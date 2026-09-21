"""Independent finite checks; does not run old large native trajectories.
Requires NumPy, mpmath, and only the bundled sources/data. Run python validate.py.
"""
from pathlib import Path
import ast, contextlib, gzip, hashlib, io, json, math, time
from collections import Counter
from copy import deepcopy
from itertools import combinations, product
from fractions import Fraction
import mpmath as mp
import numpy as np

ROOT=Path(__file__).resolve().parent
start=time.perf_counter(); counts=Counter(); errors=Counter()
data=json.loads((ROOT/'reference/fixture36.json').read_text())
r=json.load(gzip.open(ROOT/'DCU_Mass_Cell_36_RESULTS.json.gz','rt'))
assert hashlib.sha256((ROOT/'reference/fixture36.json').read_bytes()).hexdigest()==r['fixture_sha256']

def close(x,y,label,tol=2e-12):
    error=float(np.max(abs(np.asarray(x)-np.asarray(y))))
    errors[label]=max(errors[label],error)
    assert error<tol,(label,error)

# Independent original Python class and original step (AST extraction skips old runs).
ns={'__name__':'validation_native'}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile((ROOT/'reference/native_constructor.py').read_text(),'original_native','exec'),ns)
source=(ROOT/'reference/DCU_Mass_Cell_25.py').read_text()
tree=ast.parse(source)
step=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='step')
ns.update(Counter=Counter,combinations=combinations,gamma=3)
exec(compile(ast.Module(body=[step],type_ignores=[]),'original_Cell25_step','exec'),ns)
base=ns['DCUStructure']()
state=data['primary_state']
for row in state['nodes'][2:]: base.add_batch([tuple(row[:2])])
assert list(base.chains)==[row[2] for row in state['nodes']]
class PrescribedDraw:
    def __init__(self,draw): self.draw=draw;self.calls=0
    def sample(self,population,k):
        self.calls+=1
        assert k==len(self.draw) and all(v in population for v in self.draw)
        return self.draw
for saved in r['native_transition_rows']:
    W=tuple(saved['served']);draw=W+tuple(range(len(base),len(base)+3-len(W)))
    st=deepcopy(base);rng=PrescribedDraw(draw)
    F,P,event=ns['step'](st,state['F'],state['P'],rng,state['H'])
    assert rng.calls==1
    assert (len(st),F,P)==(saved['next_n'],saved['next_F'],saved['next_P'])
    assert event['cost']==saved['cost']
    assert [list(p) for p in event['new_pairs']]==saved['new_pairs']
    assert event['relief_removed']==saved['relief_removed']
    assert set(st.recorded)=={v for p in st.parents if p is not None for v in p}
    counts['full_native_next_transitions']+=1

# Actual complete node tables: validate all projected cards, paths, epochs, and late contact.
nodes={}
for world,w in data['worlds'].items():
    p=ROOT/'reference/node_tables'/f'{world}.nodes.tsv'
    assert hashlib.sha256(p.read_bytes()).hexdigest()==w['node_table_sha256']
    rows=[list(map(int,line.split())) for line in p.read_text().splitlines()]
    assert all(row[0]==i for i,row in enumerate(rows));nodes[world]=rows
    for z,c in data['cards'][world].items():
        row=rows[int(z)]
        assert row[1:3]==c['parents'] and row[3]==c['chains']
        assert row[4:8]==[c['birth_tau'],c['birth_step'],c['write_tau'],c['write_step']]
        counts['full_source_cards']+=1
for ex in data['examples']:
    rows=nodes[ex['world']];cp=ex['checkpoint']
    assert sum(row[5]<=cp['step'] for row in rows)==cp['n']
    for root,tokens in ex['paths']:
        current=tokens[-1][0]
        # Traverse the supplied path backwards, independently of forward validation.
        for z,e in reversed(tokens):
            assert current==z and 0<rows[z][7]<=cp['step']
            current=rows[z][1+e]
        assert current==root
        counts['backward_paths']+=1
late=data['late_contact_projection']['checkpoint'];rows=nodes['H0_s20350920'][:late['n']]
assert not any(row[1:3]==[10,11] for row in rows)
used=set();stack=[10,11]
while stack:
    z=stack.pop()
    if z in used:continue
    used.add(z)
    if z>=2:stack.extend(rows[z][1:3])
cost=sum((2 if 0<rows[z][7]<=late['step'] else 11)*(2*rows[z][3]-2) for z in used)
assert cost==r['late_recording_event']['cost']==128
counts['late_native_contact_costs']+=1

# High-precision pathwise quantum oracle; explicit basis amplitudes, no 27D matrices.
mp.mp.dps=70
omega=mp.exp(2j*mp.pi/3)
roots3=[mp.mpc(1),omega,omega**2]
roots={Q:[mp.exp(2j*mp.pi*j/Q) for j in range(Q)] for Q in (27,81)}
basis=list(product(range(3),repeat=2))

def oracle(h,Q):
    # amplitude(l,meter | input k), summing path labels r and probe setting m.
    ans=np.zeros((3,3,3),float)
    for l,n,k in product(range(3),repeat=3):
        amplitude=mp.mpc(0)
        for rr,m in product(range(3),repeat=2):
            a=(rr+m)%3
            amplitude+=roots3[(k*rr+n*m-l*a)%3]*roots[Q][int(h[3*a+rr])%Q]
        ans[l,n,k]=float(abs(amplitude/9)**2)
        counts['high_precision_probabilities']+=1
    return ans

paths=data['primary']['paths'];regs=sorted({z for root,t in paths for z,e in t})
for Q in (27,81):
    for entry,mode in ((False,'primary'),(True,'entry_sensitive')):
        hby={}
        for z in regs:
            def component(side,label):
                return sum(((1,0,2)[label] if entry and e else label)
                           for _,tokens in side for zz,e in tokens if zz==z)
            hby[z]=[component(paths[:2],a)-component(paths[2:],b) for a,b in basis]
        saved=r['split_interaction_sensitivity'][str(Q)][mode]
        for z in regs:
            expected=oracle(hby[z],Q)
            close(expected,saved['single_register'][str(z)]['joint_label_probabilities'],'mp_single_hit')
        for row in saved['common_phase_residue_scan']:
            j=row['aggregate_count'];expected=oracle([j*(a-b) for a,b in basis],Q)
            close(expected,row['joint_label_probabilities'],'mp_all_residues')
        for context in saved['contexts']:
            R=context['n']+context['F'];expected=np.zeros((3,3,3));norm=Fraction(0)
            for j in range(4):
                for selected in combinations(regs,j):
                    h=np.sum([hby[z] for z in selected],axis=0) if selected else [0]*9
                    p=Fraction(math.comb(R-len(regs),3-j),math.comb(R,3));norm+=p
                    expected+=float(p)*oracle(h,Q)
            assert norm==1
            close(expected,context['one_update_joint_probabilities'],'mp_native_weighted')
            one=np.mean([oracle(hby[z],Q) for z in regs],axis=0)
            close(one,context['single_support_hit_joint_probabilities'],'mp_activity_conditioned')
            counts['weighted_quantum_contexts']+=1

# All 81 complex matrix units: coarse effects identical, but instruments differ.
def complex_array(a):
    if isinstance(a,dict) and set(a)=={'real','imag'}: return complex(a['real'],a['imag'])
    if isinstance(a,list):return [complex_array(x) for x in a]
    return a
Pi=np.asarray(complex_array(r['coarse_projectors']))
for i,j in product(range(9),repeat=2):
    # p_k(rho=|i><j|)=Pi_k[j,i]; use the spectral character construction independently.
    i1,i2=divmod(i,3);j1,j2=divmod(j,3)
    exp=[]
    for k in range(3):
        value=sum(roots3[(k*m)%3]/3 for m in range(3)
                  if (i1+m)%3==j1 and (i2+m)%3==j2)
        exp.append(complex(value))
    close(exp,Pi[:,j,i],'coarse_effect_all_matrix_units')
    counts['matrix_unit_measurement_tests']+=1

answer=dict(status='PASS',checks=dict(counts),max_errors=dict(errors),
    wall_seconds=time.perf_counter()-start,
    new_native_histories=0,new_physical_interaction_derived=False,
    independent_methods=['Original DCUStructure and original Cell25 step with explicit service draws',
       'Complete source-node tables and backwards path/ancestry traversal',
       '70-digit direct basis-amplitude summation independent of the 27D matrix circuit',
       'Exact finite-group spectral projectors on every memory matrix unit'])
(ROOT/'DCU_Mass_Cell_36_TESTS.json').write_text(json.dumps(answer,indent=2)+'\n')
print(json.dumps(answer,indent=2))
