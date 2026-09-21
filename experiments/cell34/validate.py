"""Independent finite checks for Cell34. Optional dependencies: scipy and mpmath."""
from pathlib import Path
import ast,gzip,json,math,time
from collections import Counter,defaultdict
from copy import deepcopy
from fractions import Fraction
from itertools import combinations,product
from random import Random
import numpy as np
import mpmath as mp
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parent

def load():
    with gzip.open(ROOT/'reference/DCU_Mass_Cell_26_RESULTS.json.gz','rt') as f:src=json.load(f)
    p=ROOT/'DCU_Mass_Cell_34_RESULTS.json.gz'
    if not p.exists():p=ROOT/'REPLAY_RESULTS.json.gz'
    with gzip.open(p,'rt') as f:res=json.load(f)
    return src,res

def backwards(parentage,length):
    recorded={z for ps in parentage if ps is not None for z in ps}
    out=[]
    def visit(z,reverse):
        if z<2:
            if len(reverse)==length:out.append((z,tuple(reversed(reverse))))
            return
        if len(reverse)>=length:return
        for e,par in enumerate(parentage[z]):visit(par,reverse+((z,e),))
    for end in recorded:
        if end>=2:visit(end,())
    return tuple(sorted(out))

def native_source():
    env={'__name__':'cell34_native_verifier','gamma':3}
    tree=ast.parse((ROOT/'reference/native_constructor.py').read_text())
    keep=[v for v in tree.body if isinstance(v,(ast.Import,ast.ImportFrom,ast.ClassDef,ast.FunctionDef))]
    exec(compile(ast.Module(body=keep,type_ignores=[]),'native_constructor.py','exec'),env)
    tree=ast.parse((ROOT/'reference/DCU_Mass_Cell_25.py').read_text())
    keep=[v for v in tree.body if isinstance(v,(ast.Import,ast.ImportFrom))]
    step=next(v for v in ast.walk(tree) if isinstance(v,ast.FunctionDef) and v.name=='step')
    exec(compile(ast.Module(body=keep+[step],type_ignores=[]),'original_cell25_step','exec'),env)
    return env['DCUStructure'],env['step']

class Bits:
    """Independent bitmask ancestors and SEQUENTIAL first/repeat charging."""
    def __init__(self,parents):
        self.parents=[None,None];self.bits=[1,2];self.ch=[1,1];self.used=set();self.pairs={}
        for ps in parents[2:]:self.add([tuple(ps)])
    def add(self,pairs):
        before=len(self.parents);cost=0;born=[]
        for u,v in sorted(pairs):
            assert u<before and v<before and u!=v and (u,v) not in self.pairs
            ancestors=self.bits[u]|self.bits[v]
            for z in range(before):
                if ancestors>>z&1:
                    cost+=(2 if z in self.used else 11)*(2*self.ch[z]-2)
                    self.used.add(z)
            w=len(self.parents);born.append(w);self.parents.append((u,v));self.pairs[u,v]=w
            self.ch.append(self.ch[u]+self.ch[v]);self.bits.append(ancestors|(1<<w))
        return tuple(born),cost
    def step(self,F,P,rng,H):
        n=len(self.parents);R=n+F;q=min(3,R);draw=rng.sample(range(R),q)
        W=tuple(sorted(x for x in draw if x<n));f=q-len(W);fm,pm=F-f,P+2*f
        floor=max(1,pm//6);quota=floor if floor%2==0 else floor+1
        relief=min(quota,H,fm,pm) if fm>=3 and pm>=6 else 0
        pairs=tuple(ps for ps in combinations(W,2) if ps not in self.pairs)
        born,cost=self.add(pairs)
        return fm-relief+cost,pm-relief,dict(n=n,F=F,P=P,total=R,q=q,served=W,
            forced_served=f,relief_removed=relief,new_pairs=pairs,born=born,cost=cost)

def main():
    t0=time.perf_counter();src,res=load();checks=Counter();maxerr=defaultdict(float)
    cls,step=native_source()
    for envname,env in src['environments'].items():
        parentage=[None if ps is None else tuple(ps) for ps in env['parentage']]
        n=len(parentage);F=env['checkpoint']['F'];P=env['checkpoint']['P'];H=env['H']
        base=cls()
        for ps in parentage[2:]:base.add_batch([ps])
        # Every possible maintenance subset (not a sample) is legal with this F.
        Ws=tuple(W for s in range(4) for W in combinations(range(n),s))
        incidence=np.zeros((len(Ws),n),dtype=np.int16)
        weights=[]
        for i,W in enumerate(Ws):
            incidence[i,list(W)]=1;weights.append(math.comb(F,3-len(W)))
        assert sum(weights)==math.comb(n+F,3) and min(weights)>0
        checks['complete_legal_maintenance_subsets']+=len(Ws)
        # Singleton service is replayed through the ACTUAL native function.
        class Draw:
            def __init__(self,ids):self.ids=ids
            def sample(self,population,k):assert k==3;return list(self.ids)
        for z in range(n):
            native=deepcopy(base);fn,pn,event=step(native,F,P,Draw((z,n,n+1)),H)
            assert event['served']==(z,) and event['new_pairs']==() and fn>=0 and pn>=0
            assert native.parents==base.parents
            checks['actual_native_singleton_events']+=1
        for length,panel in env['panels'].items():
            k=int(length);paths=backwards(parentage,k)
            actual=tuple((r,tuple(tuple(t) for t in ts)) for r,ts in panel['paths'])
            assert paths==actual;checks['backwards_paths']+=len(paths)
            ar0=np.zeros((len(paths),n),dtype=np.int16);ar1=ar0.copy()
            for i,(root,tokens) in enumerate(paths):
                for z,e in tokens:(ar1 if e else ar0)[i,z]=1
            s0=incidence@ar0.T;s1=incidence@ar1.T
            primary=s0+s1
            rpanel=res['environments'][envname]['panels'][length]
            for row in rpanel['pair_results']:
                i,j=row['pair']
                for secondary in (False,True):
                    # INDEPENDENT class signatures over ALL legal draws, not per-register rows.
                    table=[]
                    for a,b in product(range(3),repeat=2):
                        if secondary:
                            vals=a*s0[:,i]+(1-a if a<2 else 2)*s1[:,i]-b*s0[:,j]-(1-b if b<2 else 2)*s1[:,j]
                        else:vals=a*primary[:,i]-b*primary[:,j]
                        table.append((vals%27).astype(np.int8).tobytes())
                    got=sorted(Counter(table).values(),reverse=True)
                    expected=[1]*9 if secondary else row['protected_block_dimensions']
                    assert got==expected
                    checks['all_draw_character_classifications']+=1
            # Dense exponential generators on first shared and first separated pair.
            positive=next(row for row in rpanel['pair_results'] if row['same_support'])
            negative=next(row for row in rpanel['pair_results'] if not row['same_support'])
            I=np.eye(3);G=np.diag([0.,1.,2.]);S=np.array([[0,1,0],[1,0,0],[0,0,1.]])
            logical=np.zeros((9,3));logical[[0,4,8],range(3)]=1
            for row in (positive,negative):
                i,j=row['pair'];ea=dict(paths[i][1]);eb=dict(paths[j][1])
                for Q in (27,81):
                    for sec in (False,True):
                        phases=[]
                        for z in sorted(set(ea)|set(eb)):
                            GA=(S@G@S if sec and ea.get(z) else G) if z in ea else np.zeros((3,3))
                            GB=(S@G@S if sec and eb.get(z) else G) if z in eb else np.zeros((3,3))
                            generator=np.kron(GA,I)-np.kron(I,GB)
                            U=expm(2j*np.pi*generator/Q)
                            assert np.max(abs(U.conj().T@U-np.eye(9)))<2e-14
                            if row['same_support'] and not sec:
                                err=float(np.max(abs(U@logical-logical)))
                                maxerr['dense_protected_subspace']=max(maxerr['dense_protected_subspace'],err)
                                assert err<2e-14
                            phases.append(np.diag(U));checks['matrix_exponential_generators']+=1
                        for a,b in combinations(range(9),2):
                            coincident=all(abs(x[a]-x[b])<1e-12 for x in phases)
                            if sec or not row['same_support']:assert not coincident
        # Eight designated full continuations; compare against both constructors and saved counts.
        for index in (0,127,128,255):
            h=env['histories'][index];s=deepcopy(base);bs=Bits(parentage)
            f,p,bf,bp=F,P,F,P;ra=Random(h['seed']);rb=Random(h['seed'])
            ticks=np.zeros(n,dtype=np.int64);tau=0
            for t in range(1,2049):
                f,p,e=step(s,f,p,ra,H);bf,bp,be=bs.step(bf,bp,rb,H)
                assert e==be and (f,p)==(bf,bp) and s.parents==bs.parents
                assert s.recorded==bs.used and s.chains==bs.ch
                for z in e['served']:
                    if z<n:ticks[z]+=1
                tau+=len(e['served']);checks['independent_native_transitions']+=1
                if t%128==0:
                    assert ticks.tolist()==h['counts'][t//128]
                    checks['archived_count_checkpoints']+=1
            assert (len(s),f,p,tau)==(h['final_n'],h['final_F'],h['final_P'],h['maintenance'])
    # High-precision direct amplitudes on every overlap representative / two seed blocks.
    mp.mp.dps=70
    for envname,env in src['environments'].items():
        for length,panel in env['panels'].items():
            prs=res['environments'][envname]['panels'][length]['pair_results']
            representatives={}
            for row in prs:representatives.setdefault(row['symmetric_difference'],row['pair'])
            for pair in representatives.values():
                i,j=pair;A=panel['supports'][i];B=panel['supports'][j]
                for index in (0,128,255):
                    counts=env['histories'][index]['counts'][-1]
                    a=sum(counts[z] for z in A);b=sum(counts[z] for z in B)
                    for Q in (27,81):
                        for word in range(3):
                            probs=[]
                            for out in range(3):
                                angle=mp.mpf(a-b)/Q+mp.mpf(word-out)/3
                                amplitude=sum(mp.exp(2j*mp.pi*r*angle) for r in range(3))/3
                                truth=float(abs(amplitude)**2)
                                fast=(1+2*np.cos(2*np.pi*((a-b)/Q+(word-out)/3)))**2/9
                                maxerr['high_precision_probability']=max(maxerr['high_precision_probability'],abs(fast-truth))
                                assert abs(fast-truth)<2e-13
                                checks['70_digit_decoded_probabilities']+=1
    # Readability genuinely involves coherence: reduced states and diagonals are identical.
    ww=np.exp(2j*np.pi/3)
    states=[]
    for k in range(3):
        v=np.zeros(9,complex);v[[0,4,8]]=[ww**(k*r)/math.sqrt(3) for r in range(3)]
        rho=np.outer(v,v.conj());states.append(rho)
        reduced=np.trace(rho.reshape(3,3,3,3),axis1=1,axis2=3)
        assert np.max(abs(reduced-np.eye(3)/3))<1e-14
    for i,j in combinations(range(3),2):
        assert np.max(abs(np.diag(states[i])-np.diag(states[j])))<1e-14
        distance=sum(abs(np.linalg.eigvalsh(states[i]-states[j])))/2
        assert abs(distance-1)<1e-14;checks['orthogonal_coherence_word_pairs']+=1
    print('PASS: Cell34 independent validation')
    for k,v in checks.items():print(f'  {k}: {v}')
    for k,v in maxerr.items():print(f'  {k}: {v:.3e}')
    print('Exact integer classes checked over every legal maintenance subset, not numerical rank cutoffs.')
    print('Both source and independent constructors reproduce all designated full native transitions.')
    print('This verifies the specified attached instrument; it is not evidence for a physical particle or costed memory.')
    print(f'Validation elapsed: {time.perf_counter()-t0:.3f} s; Linux only.')
if __name__=='__main__':main()
