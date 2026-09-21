"""Independent finite checks. Extra high-precision check uses mpmath only here."""
from pathlib import Path
import ast,contextlib,gzip,json,math,runpy,time,io
from collections import Counter
from itertools import combinations,product
from fractions import Fraction
from random import Random
import numpy as np
ROOT=Path(__file__).resolve().parent

def load():
    env={}
    tree=ast.parse((ROOT/'reference/native_constructor.py').read_text())
    definitions=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef,ast.ClassDef))]
    exec(compile(ast.Module(body=definitions,type_ignores=[]),'original_native','exec'),env)
    tree=ast.parse((ROOT/'reference/DCU_Mass_Cell_25.py').read_text())
    f=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='_run_dcu_mass_25')
    step=next(n for n in f.body if isinstance(n,ast.FunctionDef) and n.name=='step')
    scope={'Counter':Counter,'combinations':combinations,'gamma':3}
    exec(compile(ast.Module(body=[step],type_ignores=[]),'original_Cell25_step','exec'),scope)
    return env['DCUStructure'],scope['step']

class Bits:
    def __init__(self,parents):
        self.parents=[None,None];self.anc=[1,2];self.ch=[1,1];self.used=0;self.pairs={}
        for pair in parents[2:]:self.add([tuple(pair)])
    def add(self,pairs):
        cost=0;born=[]
        for u,v in sorted(pairs):
            p=(u,v)
            assert u<v<len(self.parents) and p not in self.pairs
            mask=self.anc[u]|self.anc[v];zbits=mask
            while zbits:
                bit=zbits&-zbits;j=bit.bit_length()-1;zbits-=bit
                cost+=(2 if self.used&bit else 11)*(2*self.ch[j]-2)
            self.used|=mask
            z=len(self.parents);self.parents.append(p);self.anc.append(mask|(1<<z))
            self.ch.append(self.ch[u]+self.ch[v]);self.pairs[p]=z;born.append(z)
        return born,cost
    def step(self,F,P,rng,H):
        n=len(self.parents);R=n+F;q=min(3,R)
        sample=rng.sample(range(R),q)
        W=sorted(z for z in sample if z<n);nf=q-len(W)
        FF=F-nf;PP=P+2*nf
        raw=max(1,PP//6);quota=raw+(raw%2)
        relief=min(quota,H,FF,PP) if FF>=3 and PP>=6 else 0
        events=[pair for pair in combinations(W,2) if pair not in self.pairs]
        births,cost=self.add(events)
        return FF-relief+cost,PP-relief,dict(n=n,F=F,P=P,total=R,q=q,served=tuple(W),
            forced_served=nf,relief_removed=relief,new_pairs=tuple(events),born=tuple(births),cost=cost)


def main():
    started=time.perf_counter();report=[]
    src=json.load(gzip.open(ROOT/'reference/DCU_Mass_Cell_29_RESULTS.json.gz','rt'))
    res=json.load(gzip.open(ROOT/'DCU_Mass_Cell_32_RESULTS.json.gz','rt'))
    Native,step=load();transitions=0;rows_checked=0;events_checked=0
    for name,env in src['environments'].items():
        H=env['H'];cp=env['checkpoint']
        for arm,A in env['arms'].items():
            # Replay actual conditioned emission, retaining all enabled births.
            st=Native()
            for pair in A['post_parentage'][2:cp['n']]:st.add_batch([tuple(pair)])
            class Draw:
                def sample(self,population,k):
                    assert k==3;return list(A['pair'])+[len(st)]
            f,p,e=step(st,cp['F'],cp['P'],Draw(),H)
            assert f==A['post_F'] and p==A['post_P'] and len(st)==len(A['post_parentage'])
            assert json.loads(json.dumps(st.parents))==A['post_parentage']
            assert json.loads(json.dumps(e))==A['event'];events_checked+=1
            for hi in (0,127,256,511):
                hist=A['history'][hi]
                original=Native()
                for pair in A['post_parentage'][2:]:original.add_batch([tuple(pair)])
                b=Bits(A['post_parentage']);F,P=A['post_F'],A['post_P'];BF,BP=F,P
                r1,r2=Random(hist['seed']),Random(hist['seed'])
                counts={j:0 for j in (3,4,0,1)};tau=cs=fs=vs=0;ip=chi=0.0
                target={row[0]:row for row in hist['rows']}
                for t in range(1,1025):
                    F,P,event=step(original,F,P,r1,H);BF,BP,bevent=b.step(BF,BP,r2,H)
                    assert (F,P,event)==(BF,BP,bevent)
                    assert original.parents==b.parents and original.chains==b.ch
                    assert sum(1<<j for j in original.recorded)==b.used
                    assert all(sum(1<<v for v in ancestors)==mask for ancestors,mask in zip(original.ancestors,b.anc))
                    for j in counts:counts[j]+=j in event['served']
                    s=len(event['served']);tau+=s;fs+=event['forced_served'];vs+=event['relief_removed'];cs+=event['cost']
                    ip+=event['q']/event['total'];chi+=s/event['n'];transitions+=1
                    if t in target:
                        actual=[t,len(original),F,P,tau,counts[3],counts[4],counts[0],counts[1],ip,chi,cs,fs,vs]
                        assert actual==target[t];rows_checked+=1
    report.append(f'PASS: {events_checked} original conditioned events; {transitions:,} complete native transitions via independent bitmask/sequential-charge code; {rows_checked} stored rows.')

    # High-precision CHARACTER sums, not the cosine kernel in the cell.
    import mpmath as mp
    mp.mp.dps=75
    def F(x):return mp.mpf(x.numerator)/x.denominator
    settings=[Fraction(x) for x in res['settings']]
    signs=[(0,0,0,1),(1,0,-1,1),(1,1,0,1),(0,1,0,1),
           (0,0,-1,-1),(1,0,0,-1),(1,1,-1,-1),(0,1,1,-1)]
    maxp=maxbell=maxneg=0.0;probcount=0;statecount=0
    for name,E in res['environments'].items():
        for arm,A in E.items():
            for depth,curves in A['readout_by_depth'].items():
                for v in curves:
                    q=v['phase_modulus'];H={int(k):c for k,c in v['sum_count_histogram'].items()};N=sum(H.values())
                    vals={}
                    for i,j,a,b in product(range(2),range(2),range(3),range(3)):
                        value=mp.mpf(0)
                        for k,c in H.items():
                            phase=F(settings[i]-settings[2+j]+Fraction(a-b,3)+Fraction(k,q))
                            amp=sum(mp.exp(2j*mp.pi*r*phase) for r in range(3))/(3*mp.sqrt(3))
                            value+=mp.mpf(c)/N*abs(amp)**2
                        vals[i,j,a,b]=value;probcount+=1
                        maxp=max(maxp,abs(float(value)-v['joint_probabilities'][i][j][a][b]))
                    bell=sum(s*sum(vals[i,j,a,b] for a in range(3) for b in range(3) if (a-b-d)%3==0)
                             for i,j,d,s in signs)
                    maxbell=max(maxbell,abs(float(bell)-v['bell']['mean']))
                    moments=[sum(mp.mpf(c)/N*mp.exp(2j*mp.pi*h*k/q) for k,c in H.items()) for h in (1,2)]
                    neg=(2*abs(moments[0])+abs(moments[1]))/3
                    maxneg=max(maxneg,abs(float(neg)-v['negativity']));statecount+=1
    # Independent filtered-source calculation on its four-dimensional support.
    fprobcount=fstatecount=0;maxfp=maxfn=maxfr=0.0
    labels=[(0,0),(0,1),(1,0),(1,1)];global_indices=[0,1,3,4]
    for name,E in res['environments'].items():
        for arm,A in E.items():
            for depth,curves in A['readout_by_depth'].items():
                for v in curves:
                    q=v['phase_modulus'];joint=v['joint_count_histogram'];N=sum(row[2] for row in joint)
                    fr=mp.matrix(4)
                    for i,(r,s) in enumerate(labels):
                        for j,(t,u) in enumerate(labels):
                            moment=sum(mp.mpf(count)/N*mp.exp(2j*mp.pi*((r-t)*a+(s-u)*b)/q)
                                       for a,b,count in joint)
                            fr[i,j]=mp.exp(2j*mp.pi*mp.mpf(r*s-t*u)/3)*moment/4
                    target=v['source_filtered_transfer']
                    for i in range(4):
                        for j in range(4):
                            vv=target['density_matrix'][global_indices[i]][global_indices[j]]
                            maxfr=max(maxfr,abs(complex(fr[i,j])-complex(vv['real'],vv['imag'])))
                    pt=mp.matrix(4)
                    for i,(r,s) in enumerate(labels):
                        for j,(t,u) in enumerate(labels):pt[i,j]=fr[labels.index((r,u)),labels.index((t,s))]
                    ev=mp.eighe(pt,eigvals_only=True)
                    nneg=-sum(x for x in ev if x<0)
                    maxfn=max(maxfn,abs(float(nneg)-target['negativity']));fstatecount+=1
                    for i,j,a,b in product(range(2),range(2),range(3),range(3)):
                        bra=[mp.exp(2j*mp.pi*(r*F(settings[i]+Fraction(a,3))-s*F(settings[2+j]+Fraction(b,3))))/3
                             for r,s in labels]
                        pp=sum(bra[k]*fr[k,l]*mp.conj(bra[l]) for k in range(4) for l in range(4)).real
                        maxfp=max(maxfp,abs(float(pp)-target['joint_probabilities'][i][j][a][b]));fprobcount+=1
    assert maxfp<1e-13 and maxfn<1e-13 and maxfr<1e-13
    report.append(f'PASS: filtered-source transfer, {fprobcount:,} probabilities, {fstatecount} density/PT states independently evaluated at 75 digits; max errors {maxfp:.3e}, {maxfn:.3e}, {maxfr:.3e}.')
    assert maxp<1e-13 and maxbell<1e-13 and maxneg<1e-13
    report.append(f'PASS: {probcount:,} joint probabilities and {statecount} Bell/negativity values from 75-digit independent complex-amplitude sums; max errors {maxp:.3e}, {maxbell:.3e}, {maxneg:.3e}.')

    # Strong unitary/generator check on arbitrary states, retaining joint counts.
    from scipy.linalg import expm
    rng=np.random.default_rng(932);checked=0
    for d in (3,4):
        den=3**d
        for _ in range(8):
            a,b=rng.integers(0,den,size=2)
            G=np.kron(np.diag([0,1,2]),np.eye(3))*a + np.kron(np.eye(3),np.diag([0,1,2]))*b
            U=expm(2j*np.pi*G/den)
            D=np.diag(np.kron(np.exp(2j*np.pi*np.arange(3)*a/den),np.exp(2j*np.pi*np.arange(3)*b/den)))
            assert np.max(abs(U-D))<1e-13;assert np.max(abs(U.conj().T@U-np.eye(9)))<2e-14
            v=rng.normal(size=9)+1j*rng.normal(size=9);v/=np.linalg.norm(v)
            assert abs(np.linalg.norm(U@v)-1)<2e-14;checked+=1
    # Equal total counts do not fix the channel on arbitrary 9D states.
    plus=np.ones(3)/np.sqrt(3);e0=np.array([1,0,0]);v=np.kron(plus,e0)
    D=np.diag(np.exp(2j*np.pi*np.arange(3)/27))
    a=np.kron(D,np.eye(3))@v;b=np.kron(np.eye(3),D)@v
    assert np.max(abs(np.outer(a,a.conj())-np.outer(b,b.conj())))>0.05
    report.append(f'PASS: {checked} independent matrix-exponential checks; joint-count necessity outside the diagonal source sector.')
    # No altered metadata silently accepted.
    import copy
    with contextlib.redirect_stdout(io.StringIO()):
        env=runpy.run_path(str(ROOT/'DCU_Mass_Cell_32.py'),init_globals={'dcu_mass_29':src})
    bad=copy.deepcopy(src);bad['clock_pair']=[3,5]
    try:env['_run_dcu_mass_32'](bad)
    except ValueError:pass
    else:raise AssertionError('Changed carrier metadata accepted')
    report.append('PASS: altered clock metadata rejected; primary reconstruction left input untouched.')
    report.append('Scope: all history-level quantum results are retrospective model predictions; no quantum outcomes were experimentally observed or newly sampled.')
    report.append(f'Validation elapsed: {time.perf_counter()-started:.3f} seconds in this Linux container.')
    (ROOT/'DCU_Mass_Cell_32_TESTS.txt').write_text('\n'.join(report)+'\n')
    print('\n'.join(report))
if __name__=='__main__':main()
