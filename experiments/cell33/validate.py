"""Independent finite source replay, echo matrices and high-precision checks."""
from pathlib import Path
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations,product
from math import comb,sqrt
import ast,gzip,json,time,hashlib
import numpy as np
import mpmath as mp
from scipy.linalg import expm
from reproduce import definitions, normalize
ROOT=Path(__file__).resolve().parent

def run():
    t0=time.perf_counter();mp.mp.dps=65
    saved=json.load(gzip.open(ROOT/'reference/DCU_Mass_Cell_29_RESULTS.json.gz','rt'))
    result=json.load(gzip.open(ROOT/'DCU_Mass_Cell_33_RESULTS.json.gz','rt'))
    Native=definitions(ROOT/'reference/native_constructor.py')['DCUStructure']
    tree=ast.parse((ROOT/'reference/DCU_Mass_Cell_25.py').read_text())
    f=next(v for v in ast.walk(tree) if isinstance(v,ast.FunctionDef) and v.name=='step')
    env={'Counter':Counter,'combinations':combinations,'gamma':3}
    exec(compile(ast.Module(body=[f],type_ignores=[]),'<original Cell25 step>','exec'),env)
    step=env['step']
    full=0;two_step=0;maxerr=0.;matrices=0;entries=0;higherr=0.;echo_static=0
    class Draw:
        def __init__(self,served,n):self.draw=list(served)+list(range(n,n+3-len(served)))
        def sample(self,pop,k):
            assert k==3 and len(set(self.draw))==3 and all(x in pop for x in self.draw)
            return self.draw
    for en in ('no_relief','relief_H6'):
        H=0 if en=='no_relief' else 6
        for arm in ('existing_pair','recorded_pair','first_use_pair'):
            inp=saved['environments'][en]['arms'][arm]
            out=result['environments'][en][arm]
            base=Native()
            for z,pp in enumerate(inp['post_parentage'][2:],2):
                born,_=base.add_batch([tuple(pp)]);assert born==[z]
            F,P=inp['post_F'],inp['post_P'];n=len(base);R=n+F
            ex=out['exact_two_step'];groups=Counter()
            for ss,weight,cost,nn,ff,pr in ex['validation_records']:
                st=deepcopy(base);f2,p2,e=step(st,F,P,Draw(ss,n),H)
                assert (len(st),f2,p2,e['cost'])==(nn,ff,pr,cost)
                assert st.recorded==set(x for pp in st.parents[2:] for x in pp)
                assert tuple(st.parents[:n])==tuple(base.parents)
                assert weight==comb(F,3-len(ss))
                j=len(set(ss)&{3,4});groups[j,len(st)+f2]+=weight;full+=1
            assert tuple((j,rn,v) for (j,rn),v in sorted(groups.items()))==tuple(map(tuple,ex['grouped_next_pool']))
            den=comb(R,3);e1=Fraction(0);e2=Fraction(0);e12=Fraction(0)
            for (j,rn),wt in groups.items():
                # Explicit finite hypergeometric second-step tag law, separate from 6/R.
                pr2={k:Fraction(comb(2,k)*comb(rn-2,3-k),comb(rn,3)) for k in range(3)}
                assert sum(pr2.values())==1
                ej=sum(k*p for k,p in pr2.items())
                assert ej==Fraction(6,rn)
                e1+=Fraction(wt,den)*j;e2+=Fraction(wt,den)*ej;e12+=Fraction(wt,den)*j*ej
            assert e12-e1*e2==Fraction(ex['covariance']) and Fraction(ex['covariance'])<0
            two_step+=1
            # Explicit phases at high precision from the integer two-window hist.
            for depth in (3,4):
                q=3**depth;o=out['readout_by_depth'][str(depth)]
                hist=o['joint_two_window_histogram'];assert sum(x[-1] for x in hist)==512
                om=mp.matrix(9,1);om[0]=om[4]=om[8]=1/mp.sqrt(3)
                filt=mp.matrix(9,1)
                for a,b in product(range(2),repeat=2):filt[3*a+b]=mp.exp(2j*mp.pi*a*b/3)/2
                for mode in ('free','echo'):
                    for name,v in (('direct',om),('filtered',filt)):
                        rho=mp.matrix(9)
                        for a,b,c,d,count in hist:
                            aeff=a+c if mode=='free' else a-c
                            beff=b+d if mode=='free' else b-d
                            u=[mp.exp(2j*mp.pi*((i//3)*aeff+(i%3)*beff)/q) for i in range(9)]
                            for i in range(9):
                                if v[i]==0:continue
                                for j in range(9):
                                    if v[j]==0:continue
                                    rho[i,j]+=mp.mpf(count)/512*u[i]*v[i]*mp.conj(u[j]*v[j])
                        expect=o[mode][name]['density_matrix']
                        for i,j in product(range(9),repeat=2):
                            got=complex(expect[i][j]['real'],expect[i][j]['imag'])
                            higherr=max(higherr,float(abs(rho[i,j]-got)));entries+=1
                        fid=sum(mp.conj(v[i])*rho[i,j]*v[j] for i,j in product(range(9),repeat=2))
                        higherr=max(higherr,abs(float(mp.re(fid))-o[mode][name]['fidelity_to_initial']))
                        pt=mp.matrix(9)
                        for a,b,c,d in product(range(3),repeat=4):pt[3*a+b,3*c+d]=rho[3*a+d,3*c+b]
                        eigen=mp.eigsy(pt,eigvals_only=True) if all(mp.im(x)==0 for x in pt) else mp.eighe(pt,eigvals_only=True)
                        nn=-sum(min(mp.re(x),0) for x in eigen)
                        higherr=max(higherr,abs(float(nn)-o[mode][name]['negativity']));matrices+=1
                # Product-marginal channel versus all early/late pairings, independent compact expansion.
                h1=Counter();h2=Counter()
                for a,b,c,d,w in hist:h1[a,b]+=w;h2[c,d]+=w
                joint=Counter()
                for (a,b),v1 in h1.items():
                    for (c,d),v2 in h2.items():joint[(a-c)%q,(b-d)%q]+=v1*v2
                assert sum(joint.values())==512**2
                for name,v in (('direct',np.array([complex(x) for x in om]).flatten()),
                               ('filtered',np.array([complex(x) for x in filt]).flatten())):
                    rho=np.zeros((9,9),complex)
                    for (a,b),w in joint.items():
                        u=np.exp(2j*np.pi*(np.repeat(np.arange(3),3)*a+np.tile(np.arange(3),3)*b)/q)
                        x=v*u;rho+=w/512**2*np.outer(x,x.conj())
                    exp=o['independent_windows'][name]['density_matrix']
                    exp=np.array([[complex(x['real'],x['imag']) for x in row] for row in exp])
                    maxerr=max(maxerr,float(np.max(abs(rho-exp))))
    # Independent dense exponential construction of pulses and phase gates.
    J=np.eye(3)[[2,1,0]];G=np.diag([0,1,2])
    for q in (27,81):
        for a,b in product(range(q),repeat=2):
            U=J@np.diag(np.exp(2j*np.pi*np.arange(3)*b/q))@J@np.diag(np.exp(2j*np.pi*np.arange(3)*a/q))
            V=np.exp(4j*np.pi*b/q)*np.diag(np.exp(2j*np.pi*np.arange(3)*(a-b)/q))
            maxerr=max(maxerr,float(np.max(abs(U-V))))
            if a==b:
                assert np.max(abs(U-np.exp(4j*np.pi*b/q)*np.eye(3)))<2e-14;echo_static+=1
    # J=exp(-i*pi/2*(I-J)) is an independent Hermitian-generator identity.
    assert np.max(abs(expm(-1j*np.pi*(np.eye(3)-J)/2)-J))<2e-14
    for a,b in ((0,0),(1,2),(7,9),(30,21)):
        U=expm(-1j*np.pi*(np.eye(3)-J)/2)@expm(2j*np.pi*b*G/81)@expm(-1j*np.pi*(np.eye(3)-J)/2)@expm(2j*np.pi*a*G/81)
        V=np.exp(4j*np.pi*b/81)*expm(2j*np.pi*(a-b)*G/81)
        maxerr=max(maxerr,float(np.max(abs(U-V))))
    assert higherr<3e-13 and maxerr<3e-13
    print('CELL 33 — INDEPENDENT VALIDATION')
    print('Original-source complete conditional transitions checked:',full)
    print('Exact two-step covariance checks against explicit second-step tag law:',two_step)
    print('High-precision density/negativity/fidelity states:',matrices,'at',mp.mp.dps,'digits')
    print('High-precision density entries compared:',entries)
    print('Maximum high-precision discrepancy:',higherr)
    print('Exact static-noise echo instances checked:',echo_static)
    print('All 729+6561 phase pairs and four independent matrix exponentials checked.')
    print('All 24 independent-window states checked by full empirical-product histograms.')
    print('Maximum independent floating-point discrepancy:',maxerr)
    print('No Monte Carlo outcome agreement used as a pass criterion.')
    print('Runtime seconds:',time.perf_counter()-t0)
    print('PASS')
if __name__=='__main__':run()
