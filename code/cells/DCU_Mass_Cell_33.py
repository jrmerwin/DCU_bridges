# CELL 33 — fixed midpoint echo, timing correlations, and preparation transfer.
# Paste after Cell32. Reads Cell29 histories; Cell32 is an optional regression check.
# NumPy + standard library. No new Monte Carlo histories, native controls, or fits.
# ADDED ideal quantum control: J|r>=|2-r> on each qutrit at t=512 and t=1024.
# Timing is an external iteration schedule; no physical gate cost/controller is derived.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations, product
from math import comb, sqrt


def _run_dcu_mass_33(saved=None, previous=None):
    import numpy as np
    if saved is None:
        if 'dcu_mass_29' not in globals():
            raise RuntimeError('Run Cell29 first, or use the standalone Cell33 package.')
        saved = dcu_mass_29
    if previous is None:
        previous = globals().get('dcu_mass_32')
    before = deepcopy(saved)
    Q = Fraction
    times = tuple(saved['checkpoints'])
    seeds = tuple(saved['seeds'])
    arms = ('existing_pair', 'recorded_pair', 'first_use_pair')
    envs = ('no_relief', 'relief_H6')
    if (saved['horizon'], tuple(saved['clock_pair']), tuple(saved['phase_depths'])) != (1024,(3,4),(3,4)):
        raise ValueError('Frozen observation window, clock, or phase convention changed.')
    if times != (0,1,8,32,128,512,1024) or len(seeds) != 512 or tuple(saved['arm_order']) != arms:
        raise ValueError('Frozen checkpoints, cohort, or arms changed; no replacement selected.')
    mid, end = times.index(512), times.index(1024)
    col = {name: saved['history_columns'].index(name) for name in ('step','N3','N4')}
    r = np.arange(3)
    settings = (Q(0),Q(4,27),-Q(2,27),Q(2,27))
    xs, ys = np.array(list(map(float, settings[:2]))), np.array(list(map(float, settings[2:])))
    theta = xs[:,None,None,None]-ys[None,:,None,None]+(r[None,None,:,None]-r[None,None,None,:])/3
    W = np.zeros((2,2,3,3),int)
    for i,j,d,sgn in ((0,0,0,1),(1,0,-1,1),(1,1,0,1),(0,1,0,1),
                       (0,0,-1,-1),(1,0,0,-1),(1,1,-1,-1),(0,1,1,-1)):
        for a,b in product(range(3),repeat=2):
            if (a-b-d)%3==0: W[i,j,a,b]+=sgn
    projectors=np.empty((2,2,3,3,9,9),complex)
    for i,j,a,b in product(range(2),range(2),range(3),range(3)):
        av=np.exp(-2j*np.pi*r*(xs[i]+a/3))/sqrt(3)
        bv=np.exp(+2j*np.pi*r*(ys[j]+b/3))/sqrt(3)
        v=np.kron(av,bv); projectors[i,j,a,b]=np.outer(v,v.conj())
    Bop=np.einsum('xyab,xyabij->ij',W,projectors)
    F3=np.exp(2j*np.pi*r[:,None]*r[None,:]/3)/sqrt(3)
    omega=np.zeros(9,complex);omega[[0,4,8]]=1/sqrt(3)
    P01=np.diag([1,1,0])
    filtered=np.kron(P01@F3@P01,np.eye(3))@omega
    herald=float(np.vdot(filtered,filtered).real)
    assert abs(herald-4/9)<2e-14
    filtered/=sqrt(herald)
    source_vectors={'direct':omega,'filtered':filtered}
    source_states={name:np.outer(v,v.conj()) for name,v in source_vectors.items()}
    rows=np.repeat(r,3); cols=np.tile(r,3)
    X=np.roll(np.eye(3),1,axis=0)
    S=np.eye(3)[[1,0,2]]
    J=np.eye(3)[[2,1,0]]
    assert np.array_equal(X@X@S@X,J) and np.array_equal(J@J,np.eye(3))
    # Native entry reflection is S=(0 1), NOT the echo permutation J=(0 2).
    # Their active implementation as pulses remains an added instrument control.
    maxerr=Counter(); checks=Counter()

    def stats(a):
        a=np.asarray(a,float)
        return dict(n=len(a),mean=float(a.mean()),se=float(a.std(ddof=1)/sqrt(len(a))) if len(a)>1 else None)

    def neg(rho):
        ev=np.linalg.eigvalsh(rho.reshape(3,3,3,3).transpose(0,3,2,1).reshape(9,9))
        return float(-ev[ev<0].sum())

    def probabilities(rho):
        return np.real(np.einsum('ij,xyabji->xyab',rho,projectors))

    def phases(counts,q):
        return np.exp(2j*np.pi*((counts[:,0,None]*rows+counts[:,1,None]*cols)%q)/q)

    def density_report(rho,vector):
        rho=(rho+rho.conj().T)/2
        prob=probabilities(rho)
        assert np.linalg.eigvalsh(rho).min()>-2e-13
        assert abs(np.trace(rho)-1)<2e-13 and prob.min()>-2e-13
        assert np.max(abs(prob.sum(axis=(2,3))-1))<2e-13
        assert np.max(abs(prob.sum(axis=3)[:,0]-prob.sum(axis=3)[:,1]))<2e-13
        assert np.max(abs(prob.sum(axis=2)[0]-prob.sum(axis=2)[1]))<2e-13
        checks['density_matrices']+=1;checks['outcome_probabilities']+=36
        return dict(density_matrix=rho,joint_probabilities=prob,bell=float(np.sum(W*prob)),
            negativity=neg(rho),fidelity_to_initial=float(np.vdot(vector,rho@vector).real),
            purity=float(np.trace(rho@rho).real))

    def ensemble_report(counts,q):
        u=phases(counts,q);channel=u.T@u.conj()/len(u)
        assert np.max(abs(channel.diagonal()-1))<2e-14
        assert np.linalg.eigvalsh((channel+channel.conj().T)/2).min()>-2e-12
        output={}
        for name,v in source_vectors.items():
            rho=source_states[name]*channel
            report=density_report(rho,v)
            states=u*v[None,:]
            f=abs(states@v.conj())**2
            bell=np.real(np.einsum('ni,ij,nj->n',states.conj(),Bop,states))
            assert abs(f.mean()-report['fidelity_to_initial'])<4e-14
            assert abs(bell.mean()-report['bell'])<4e-14
            report['fidelity_stats']=stats(f);report['bell_stats']=stats(bell)
            # Leave-one-history-out matrices: uncertainty is history sampling, not shots.
            pure=states[:,:,None]*states[:,None,:].conj()
            loo=(len(states)*rho-pure)/(len(states)-1)
            pt=loo.reshape(-1,3,3,3,3).transpose(0,1,4,3,2).reshape(-1,9,9)
            eig=np.linalg.eigvalsh(pt);ln=-np.minimum(eig,0).sum(axis=1)
            report['negativity_jackknife_SE']=float(np.sqrt((len(states)-1)/len(states)*sum((ln-ln.mean())**2)))
            report['seed_blocks']=[]
            for sl in (slice(0,256),slice(256,512)):
                ub=u[sl];rhob=source_states[name]*(ub.T@ub.conj()/len(ub))
                report['seed_blocks'].append(dict(bell=float(np.trace(Bop@rhob).real),
                    negativity=neg(rhob),fidelity_to_initial=float(np.vdot(v,rhob@v).real)))
            output[name]=report
        z=np.exp(2j*np.pi*(counts.sum(axis=1)%q)/q)
        C1=complex(z.mean());C2=complex((z*z).mean())
        assert abs(output['direct']['negativity']-(2*abs(C1)+abs(C2))/3)<8e-14
        output['phase_moments']={'C1':C1,'C2':C2}
        output['residue_histogram']=tuple((a,b,v) for (a,b),v in sorted(Counter(
            (int(a%q),int(b%q)) for a,b in counts).items()))
        return output,channel,u

    def choose(n,k):return comb(n,k) if 0<=k<=n else 0

    def exact_two_step(inp,H):
        # Reimplement Eq.(14) only for finite conditional enumeration; validator
        # checks every resulting state against the ORIGINAL native step.
        pp=inp['post_parentage'];n=len(pp);F=int(inp['post_F']);P=int(inp['post_P']);R=F+n
        ancestors=[{0},{1}];L=[0,0];weights=[1,1];oldpairs=set()
        for z,pair in enumerate(pp[2:],2):
            a,b=map(int,pair);assert a<z and b<z and a!=b
            ancestors.append(ancestors[a]|ancestors[b]|{z})
            weights.append(weights[a]+weights[b]);L.append(2*weights[-1]-2)
            oldpairs.add(tuple(sorted((a,b))))
        recorded=set(p for pair in pp[2:] for p in pair)
        groups=Counter();records=[];subsets=0
        for s in range(4):
            wt=choose(F,3-s)
            if wt==0:continue
            for served in combinations(range(n),s):
                subsets+=1;f=3-s
                new=tuple(pair for pair in combinations(served,2) if pair not in oldpairs)
                uses=Counter(z for a,b in new for z in ancestors[a]|ancestors[b])
                cost=2*sum(L[z]*v for z,v in uses.items())+9*sum(L[z] for z in uses if z not in recorded)
                pm=P+2*f;bm=F-f
                quota=2*((max(1,pm//6)+1)//2)
                relief=min(quota,H,bm,pm) if bm>=3 and pm>=6 else 0
                Rnext=n+len(new)+F-f-relief+cost
                j=int(3 in served)+int(4 in served)
                groups[j,Rnext]+=wt
                records.append((served,wt,cost,n+len(new),F-f-relief+cost,pm-relief))
        den=comb(R,3)
        assert sum(groups.values())==den
        mass=[sum(v for (j,rn),v in groups.items() if j==k) for k in range(3)]
        byj=[sum((Q(v)*Q(6,rn) for (j,rn),v in groups.items() if j==k),Q(0))/mass[k]
             for k in range(3)]
        EJ1=sum(Q(k*mass[k],den) for k in range(3))
        EJ2=sum(Q(mass[k],den)*byj[k] for k in range(3))
        cov=sum(Q(k*mass[k],den)*byj[k] for k in range(3))-EJ1*EJ2
        assert EJ1==Q(6,R)
        checks['exact_native_service_subsets']+=subsets
        # Enumeration is new exact algebra; not extra stochastic histories.
        return dict(pool=R,subsets=subsets,first_increment_probabilities=tuple(Q(v,den) for v in mass),
            expected_second_increment_by_first=tuple(byj),expected_first_increment=EJ1,
            expected_second_increment=EJ2,covariance=cov,covariance_nonzero_exactly=(cov!=0),
            grouped_next_pool=tuple((j,rn,v) for (j,rn),v in sorted(groups.items())),
            validation_records=tuple(records))

    output={}
    for en in envs:
        output[en]={};H=0 if en=='no_relief' else 6
        for arm in arms:
            inp=saved['environments'][en]['arms'][arm]
            histories=inp['history']
            assert tuple(h['seed'] for h in histories)==seeds
            assert all(tuple(row[col['step']] for row in h['rows'])==times for h in histories)
            counts=np.array([[[int(row[col['N3']]),int(row[col['N4']])] for row in h['rows']]
                            for h in histories],dtype=np.int64)
            assert np.all(counts>=0) and np.all(np.diff(counts,axis=1)>=0)
            early=counts[:,mid];late=counts[:,end]-early
            A=early.sum(axis=1);B=late.sum(axis=1);N=len(A)
            sa,sb,saa,sbb,sab=(int(z) for z in (A.sum(),B.sum(),(A*A).sum(),(B*B).sum(),(A*B).sum()))
            va=Q(saa,N)-Q(sa,N)**2;vb=Q(sbb,N)-Q(sb,N)**2
            cv=Q(sab,N)-Q(sa*sb,N*N)
            vs=Q(int(((A+B)**2).sum()),N)-Q(sa+sb,N)**2
            vd=Q(int(((A-B)**2).sum()),N)-Q(sa-sb,N)**2
            assert vs-vd==4*cv
            co_loo=((sab-A*B)-(sa-A)*(sb-B)/(N-1))/(N-2)
            co_se=float(np.sqrt((N-1)/N*np.sum((co_loo-co_loo.mean())**2)))
            timing=dict(first_half=stats(A),second_half=stats(B),difference=stats(A-B),
                sample_covariance=float(cv*N/(N-1)),covariance_jackknife_SE=co_se,
                empirical_population_variances={'early':va,'late':vb,'free_sum':vs,'echo_difference':vd},
                empirical_population_covariance=cv,variance_identity_exact=True)
            depths={}
            for depth in (3,4):
                q=3**depth
                free,cf,uf=ensemble_report(early+late,q)
                echo,ce,ue=ensemble_report(early-late,q)
                ud,ul=phases(early,q),phases(late,q)
                # Product marginals preserve within-window two-factor correlation,
                # but REMOVE the relationship between early and late windows.
                cearly=ud.T@ud.conj()/N;clate=ul.T@ul.conj()/N
                cind=cearly*clate.conj()
                independent={}
                comparisons={}
                for name,v in source_vectors.items():
                    independent[name]=density_report(source_states[name]*cind,v)
                    delta=echo[name]['density_matrix']-independent[name]['density_matrix']
                    tdist=float(np.abs(np.linalg.eigvalsh(delta)).sum()/2)
                    sf=uf*v[None,:];se=ue*v[None,:]
                    bf=np.real(np.einsum('ni,ij,nj->n',sf.conj(),Bop,sf))
                    be=np.real(np.einsum('ni,ij,nj->n',se.conj(),Bop,se))
                    ff=abs(sf@v.conj())**2;fe=abs(se@v.conj())**2
                    comparisons[name]=dict(delta_fidelity=stats(fe-ff),delta_bell=stats(be-bf),
                        delta_negativity=echo[name]['negativity']-free[name]['negativity'],
                        echo_vs_independent_window_trace_distance=tdist)
                    # Final diagonal populations are unchanged: no energy exchange
                    # for any diagonal energy assignment in this phase basis.
                    assert np.max(abs(np.diag(echo[name]['density_matrix'])-np.diag(source_states[name])))<3e-14
                    assert np.max(abs(np.diag(free[name]['density_matrix'])-np.diag(source_states[name])))<3e-14
                # Fixed original source Bell regression, even when previous is absent.
                z=np.exp(2j*np.pi*((A+B)%q)/q)
                c1=z.mean();c2=(z*z).mean()
                p=1/9+4*np.real(c1*np.exp(2j*np.pi*theta))/27+2*np.real(c2*np.exp(4j*np.pi*theta))/27
                assert np.max(abs(free['direct']['joint_probabilities']-p))<4e-14
                if previous is not None:
                    old=previous['environments'][en][arm]['readout_by_depth']
                    old=old[depth] if depth in old else old[str(depth)]
                    endold=old[-1]
                    assert abs(free['direct']['bell']-endold['bell']['mean'])<1e-13
                    assert abs(free['direct']['negativity']-endold['negativity'])<1e-13
                    assert abs(free['filtered']['negativity']-endold['source_filtered_transfer']['negativity'])<1e-13
                    checks['Cell32_terminal_regressions']+=3
                # Full explicit J D_late J D_early tensor products, all histories.
                for k in range(N):
                    local=[]
                    for side in (0,1):
                        d1=np.diag(np.exp(2j*np.pi*r*int(early[k,side]%q)/q))
                        d2=np.diag(np.exp(2j*np.pi*r*int(late[k,side]%q)/q))
                        full=J@d2@J@d1
                        target=np.exp(4j*np.pi*int(late[k,side]%q)/q)*np.diag(
                            np.exp(2j*np.pi*r*int((early[k,side]-late[k,side])%q)/q))
                        err=float(np.max(abs(full-target)));maxerr['echo_gate_identity']=max(maxerr['echo_gate_identity'],err)
                        assert err<3e-14;local.append(full)
                    u=np.kron(*local)
                    for name,v in source_vectors.items():
                        v1=u@v;v2=ue[k]*v
                        err=float(np.max(abs(np.outer(v1,v1.conj())-np.outer(v2,v2.conj()))))
                        maxerr['tensor_echo_density']=max(maxerr['tensor_echo_density'],err)
                        assert err<3e-14
                        checks['explicit_two_qutrit_echo_states']+=1
                aphase=np.exp(2j*np.pi*A/q);bphase=np.exp(2j*np.pi*B/q)
                memory=[]
                for h in (1,2):
                    xh=aphase**h;yh=bphase.conj()**h
                    joint=(xh*yh).mean();prod=xh.mean()*yh.mean();diff=joint-prod
                    lj=(N*joint-xh*yh)/(N-1)
                    lp=(N*xh.mean()-xh)/(N-1)*(N*yh.mean()-yh)/(N-1)
                    ld=lj-lp
                    memory.append(dict(harmonic=h,joint=complex(joint),product=complex(prod),
                        difference=complex(diff),real_jackknife_SE=float(np.sqrt((N-1)/N*sum((ld.real-ld.real.mean())**2))),
                        imag_jackknife_SE=float(np.sqrt((N-1)/N*sum((ld.imag-ld.imag.mean())**2)))))
                depths[depth]=dict(free=free,echo=echo,independent_windows=independent,
                    comparisons=comparisons,temporal_character_covariance=memory,
                    joint_two_window_histogram=tuple((*key,v) for key,v in sorted(Counter(
                        tuple(map(int,(*early[k],*late[k]))) for k in range(N)).items())),
                    ideal_full_count_inverse={name:dict(negativity=neg(rho),fidelity=1.,
                        bell=float(np.trace(Bop@rho).real)) for name,rho in source_states.items()})
            output[en][arm]=dict(timing=timing,readout_by_depth=depths,
                exact_two_step=exact_two_step(inp,H),new_native_histories=0)

    # Exact finite example: at midpoint and endpoint, SAME free channels; echo
    # differs. These algebraic noise laws are not claimed as DCU trajectories.
    q=27;hist_corr_mid=Counter();hist_corr_end=Counter();hist_corr_echo=Counter()
    hist_ind_mid=Counter();hist_ind_end=Counter();hist_ind_echo=Counter()
    for k in range(q):
        hist_corr_mid[k]+=q;hist_corr_end[2*k%q]+=q;hist_corr_echo[0]+=q
    for a,b in product(range(q),repeat=2):
        hist_ind_mid[a]+=1;hist_ind_end[(a+b)%q]+=1;hist_ind_echo[(a-b)%q]+=1
    assert hist_corr_mid==hist_ind_mid and hist_corr_end==hist_ind_end
    assert hist_corr_echo==Counter({0:q*q}) and set(hist_ind_echo.values())=={q}
    counterexample=dict(modulus=q,description='one factor: early uniform k; late=k versus independent uniform late',
        midpoint_free_channels_equal_exactly=True,endpoint_free_channels_equal_exactly=True,
        correlated_echo_residues=dict(hist_corr_echo),independent_echo_residues=dict(hist_ind_echo),
        direct_echo_negativity_correlated=1.,direct_echo_negativity_independent=0.,
        source_state_trace_distance=Q(2,3),native_occurrence_claimed=False)
    assert saved==before
    ideal=float(np.trace(Bop@source_states['direct']).real)
    print('CELL 33 — FIXED MIDPOINT ECHO AND TWO-WINDOW NATIVE TIMING')
    print('512 existing histories x 6 arms; ALL reused. No new stochastic histories or native update rule.')
    print('J=(0 2) on each factor at 512 and 1024; ideal added controls, not costed native gates.')
    print('Local effective counts: early-late. No individual-history correction or fitted pulse time.')
    print('Same source, filtered preparation, phase increments 1/27 and 1/81, and Bell settings as Cell32.')
    print('\nAT 1024, 1/27, DIRECT SOURCE:')
    print(' environment / arm           fidelity free/echo   negativity free/echo      Bell free/echo')
    for en in envs:
        for arm in arms:
            z=output[en][arm]['readout_by_depth'][3];a=z['free']['direct'];b=z['echo']['direct']
            print(f' {en:12s}/{arm:15s} {a["fidelity_to_initial"]:.6f}/{b["fidelity_to_initial"]:.6f}'
                f'    {a["negativity"]:.6f}/{b["negativity"]:.6f}   {a["bell"]:+.6f}/{b["bell"]:+.6f}')
    print('\nTWO EQUAL ITERATION WINDOWS: counts, covariance, exact-next-step memory')
    for en in envs:
        for arm in arms:
            z=output[en][arm];t=z['timing'];e=z['exact_two_step']
            print(f' {en:12s}/{arm:15s} means={t["first_half"]["mean"]:.6f}, {t["second_half"]["mean"]:.6f}; '
                f'Cov={t["sample_covariance"]:+.6f} +/- {t["covariance_jackknife_SE"]:.6f}; '
                f'exact Cov(J1,J2)={float(e["covariance"]):+.9g}')
    print('\nFILTERED PREPARATION, 1/27: fraction of initial negativity, free / echo')
    for en in envs:
        for arm in arms:
            z=output[en][arm]['readout_by_depth'][3];baseline=sqrt(3)/4
            print(f' {en:12s}/{arm:15s} {z["free"]["filtered"]["negativity"]/baseline:.6f} / '
                f'{z["echo"]["filtered"]["negativity"]/baseline:.6f}; '
                f'echo-vs-independent trace distance={z["comparisons"]["filtered"]["echo_vs_independent_window_trace_distance"]:.6f}')
    print(f'\nExact finite memory counterexample: identical two free channels, echo negativity 1 vs 0; trace distance 2/3.')
    print('All final populations unchanged. This diagonal phase channel supplies no diagonal energy shift.')
    print('Echo improves state alignment without necessarily restoring entanglement; no spatial time dilation inferred.')
    print('PASS:',dict(checks));print('Max numerical errors:',dict(maxerr))
    return dict(protocol='Cell33 v1 fixed midpoint echo on all saved Cell29 histories',
        pulse_times=(512,1024),phase_depths=(3,4),primary_phase_depth=3,settings=settings,
        seeds=seeds,checkpoints=times,reused_native_histories=6*len(seeds),new_stochastic_native_updates=0,
        control='J D(late) J D(early); J=X^2 S X; both factors; ideal external iteration schedule',
        ideal_bell=ideal,herald_success_probability=herald,environments=output,
        finite_memory_counterexample=counterexample,validation_counts=dict(checks),max_numerical_errors=dict(maxerr),
        native_law_changed=False,gate_cost_derived=False,autonomous_controller_implemented=False,
        pulse_time_optimized=False,history_postselection=False,quantum_non_markovianity_claimed=False,
        spatial_or_gravity_claim=False,energy_scale_calibrated=False,neutron_search_modified=False)


dcu_mass_33 = _run_dcu_mass_33()
