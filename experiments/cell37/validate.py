"""Independent checks: original native implementation, exact cyclotomic certificates,
and high-precision amplitude arithmetic. Requires NumPy + mpmath. No large growth run.
"""
from pathlib import Path
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations,product
import ast,contextlib,gzip,hashlib,io,json,math,time
import mpmath as mp
import numpy as np
ROOT=Path(__file__).resolve().parent
start=time.perf_counter();counts=Counter();errors=Counter()
data=json.load(gzip.open(ROOT/'reference/fixture37.json.gz','rt'))
r=json.load(gzip.open(ROOT/'DCU_Mass_Cell_37_RESULTS.json.gz','rt'))
rows=[list(map(int,line.split())) for line in (ROOT/'reference/H0_s20350920.nodes.tsv').read_text().splitlines()]
assert hashlib.sha256((ROOT/'reference/H0_s20350920.nodes.tsv').read_bytes()).hexdigest()==data['full_node_table_sha256']
assert hashlib.sha256((ROOT/'reference/fixture36.json').read_bytes()).hexdigest()==data['fixture36_sha256']
ns={}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile((ROOT/'reference/native_constructor.py').read_text(),'source_native','exec'),ns)
tree=ast.parse((ROOT/'reference/DCU_Mass_Cell_25.py').read_text())
step=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='step')
ns.update(Counter=Counter,combinations=combinations,gamma=3)
exec(compile(ast.Module(body=[step],type_ignores=[]),'source_step','exec'),ns)

class PrescribedDraw:
    def __init__(self,draw):self.draw=draw;self.calls=0
    def sample(self,population,k):
        assert len(self.draw)==k and len(set(self.draw))==k and all(x in population for x in self.draw)
        self.calls+=1;return self.draw

def make_original(cp):
    st=ns['DCUStructure']()
    # Reconstruct structure; cumulative insertion charges while replaying individual
    # births are NOT the historical workload. The complete saved F/P ledgers are restored.
    for p in cp['parents'][2:]:st.add_batch([tuple(p)])
    return st

def seq_cost(st,pairs):
    seen=set(st.recorded);cost=0
    for a,b in pairs:
        for z in st.ancestors[a]|st.ancestors[b]:
            cost+=(2 if z in seen else 11)*st.path_weight(z);seen.add(z)
    return cost

for name,cp in data['worlds'].items():
    assert sum(row[5]<=cp['step'] for row in rows)==cp['n']
    assert all((None if row[1]<0 else row[1:3])==pp for row,pp in zip(rows,cp['parents']))
    st=make_original(cp)
    recorded={i for i,row in enumerate(rows[:cp['n']]) if 0<row[7]<=cp['step']}
    assert st.recorded==recorded
    assert st.chains==[row[3] for row in rows[:cp['n']]]
    counts['source_native_vertices']+=len(st)
    F,P,tau,step_index=cp['F'],cp['P'],cp['tau'],cp['step']
    for target in r['native_source_routes'][name]['rounds']:
        n=len(st);anchors=tuple(target['old_anchors'])
        source=next(z for z in range(2,n) if z not in st.recorded and z not in anchors)
        assert source==target['source']
        e=target['event'];pairs=[tuple(x) for x in e['new_pairs']]
        assert seq_cost(st,pairs)==seq_cost(st,list(reversed(pairs)))==e['cost']
        rng=PrescribedDraw(tuple(e['served']))
        F2,P2,ev=ns['step'](st,F,P,rng,0)
        assert rng.calls==1
        assert ev['cost']==e['cost'] and F2==e['next_F'] and P2==e['next_P']
        assert len(st)==e['next_n'] and [list(x) for x in ev['new_pairs']]==e['new_pairs']
        assert tuple(st.pair_to_id[tuple(sorted((source,a)))] for a in anchors)==tuple(target['new_anchors'])
        assert math.comb(F+n,3)==Fraction(e['probability']).denominator
        assert Fraction(target['future_fixed_anchor_ready_probability_upper_bound'])==Fraction(2*(n+2),(n-2)*(n-1))
        tau+=3;step_index+=1
        assert (tau,step_index)==(e['next_tau'],e['next_step'])
        F,P=F2,P2;counts['original_complete_resource_bursts']+=1
        counts['order_independent_charge_checks']+=2

# Every early service subset against actual original step, not just proposed protocol.
cp=data['worlds']['early'];base=make_original(cp)
for e in r['exact_native_rows']:
    W=tuple(e['served']);draw=W+tuple(range(len(base),len(base)+3-len(W)))
    st=deepcopy(base);rng=PrescribedDraw(draw)
    F,P,ev=ns['step'](st,cp['F'],cp['P'],rng,0)
    assert ev['cost']==e['cost'] and F==e['next_F'] and P==e['next_P'] and len(st)==e['next_n']
    assert [list(p) for p in ev['new_pairs']]==e['new_pairs']
    assert Fraction(math.comb(cp['F'],3-len(W)),math.comb(cp['F']+len(base),3))==Fraction(e['probability'])
    counts['original_complete_next_states']+=1

# Source paths checked against real parent entries, not inferred support equality.
paths=data['original_paths']
for root,tokens in paths:
    current=tokens[-1][0]
    for z,e in reversed(tokens):
        assert current==z and 0<rows[z][7]<=data['worlds']['early']['step']
        current=rows[z][1+e]
    assert current==root;counts['original_backward_paths']+=1
left=Counter(tuple(t) for _,tokens in paths[:2] for t in tokens)
right=Counter(tuple(t) for _,tokens in paths[2:] for t in tokens)
assert left==right

# Exact Q(omega) arithmetic: omega^2+omega+1=0. No floating rank or root tolerances.
def plus(x,y):return (x[0]+y[0],x[1]+y[1])
def times(x,y):
    a,b=x;c,d=y
    return (a*c-b*d,a*d+b*c-b*d)
def scale(x,c):return (x[0]*c,x[1]*c)
w=((Fraction(1),Fraction(0)),(Fraction(0),Fraction(1)),(Fraction(-1),Fraction(-1)))
zero=(Fraction(0),Fraction(0))
for d in range(3):
    assert sum(1 for j in range(3) if (j+d)%3==j)==(3 if d==0 else 0)
    s=zero
    for a in range(3):s=plus(s,w[(a*d)%3])
    assert s==((Fraction(3),Fraction(0)) if d==0 else zero)
    counts['exact_character_sums']+=1
for k,i,j,u,v in product(range(3),repeat=5):
    target=scale(w[(k*(u-v-i+j))%3],Fraction(1,9))
    renewal=zero
    for a in range(3):
        b=(a+k)%3
        renewal=plus(renewal,scale(w[(k*(u-v)+(a-b)*(i-j))%3],Fraction(1,27)))
    assert target==renewal;counts['exact_code_channel_entries']+=1

# High precision constructs matrices by component amplitudes, independently of the
# main Kronecker/projector implementation, for every full 9D effect entry.
mp.mp.dps=75
om=mp.exp(2j*mp.pi/3)
basis=list(product(range(3),repeat=2))
qmat=[];mm=[]
for k in range(3):
    M=np.zeros((9,9),complex)
    for i,(a,b) in enumerate(basis):
        for j,(c,d) in enumerate(basis):
            v=om**(k*((a-c)%3))/3 if (a-c)%3==(b-d)%3 else mp.mpc(0)
            M[i,j]=complex(v);counts['high_precision_target_entries']+=1
    qmat.append(M)
for a,b in basis:
    k=(b-a)%3;M=np.zeros((9,9),complex)
    for i,(u,v) in enumerate(basis):
        for j,(x,y) in enumerate(basis):
            value=om**(k*u+a*x-b*y)/(3*mp.sqrt(3)) if u==v else mp.mpc(0)
            M[i,j]=complex(value);counts['high_precision_renewal_entries']+=1
    mm.append(M)

def chk(x,y,label):
    z=float(np.max(abs(np.asarray(x)-np.asarray(y))));errors[label]=max(errors[label],z);assert z<3e-13
for k in range(3):
    eff=sum((M.conj().T@M for M,(a,b) in zip(mm,basis) if (b-a)%3==k),np.zeros((9,9),complex))
    chk(eff,qmat[k],'high_precision_same_effects')
for i,j in product(range(3),repeat=2):
    E=np.zeros((9,9),complex);E[4*i,4*j]=1
    for k in range(3):
        expected=qmat[k]@E@qmat[k]
        actual=sum((M@E@M.conj().T for M,(a,b) in zip(mm,basis) if (b-a)%3==k),np.zeros((9,9),complex))
        chk(expected,actual,'high_precision_code_maps');counts['high_precision_code_branch_maps']+=1
# Outside-code control exactly disjoint supports at the output.
E=np.zeros((9,9),complex);E[1,1]=1
a=sum(Q@E@Q for Q in qmat);b=sum(M@E@M.conj().T for M in mm)
chk(np.abs(np.linalg.eigvalsh(a-b)).sum()/2,r['alternative']['out_of_code_control_trace_distance'],'outside_control')
assert r['primary_compiler_success'] is False
assert r['full_operation_cost'] is None
assert all(b['complete_measurement_work'] is None for route in r['native_source_routes'].values() for b in route['rounds'])
counts['unknown_cost_not_zero']=5
summary=dict(counts=dict(counts),max_numerical_differences=dict(errors),elapsed_seconds=time.perf_counter()-start,platform='Linux',new_unconditioned_native_trajectories=0)
(ROOT/'DCU_Mass_Cell_37_TESTS.json').write_text(json.dumps(summary,indent=2)+'\n')
text='PASS: Cell37 independent finite validation.\n'+json.dumps(summary,indent=2)+'\nNumerical agreement is implementation validation, not a physical accuracy or complete native compiler claim.\n'
(ROOT/'DCU_Mass_Cell_37_TESTS.txt').write_text(text)
print(text)
