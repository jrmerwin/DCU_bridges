# CELL 32 — native service histories -> a finite quantum noise channel and Bell readout.
# Paste after Cell29 (Cells30/31 and neutron-search cells are NOT prerequisites).
# NumPy + stdlib only. Uses every SAVED history; no new service law or fitted noise.
# Standalone: extract the Cell32 bundle, then run `python reproduce.py`.
# The attached complex quantum source and service-to-phase rule remain hypotheses.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from math import sqrt


def _run_dcu_mass_32(saved=None):
    import numpy as np
    if saved is None:
        if 'dcu_mass_29' not in globals():
            raise RuntimeError('Run Cell29 first, or use the standalone Cell32 reproduce.py.')
        saved = dcu_mass_29
    before = deepcopy(saved)
    if (tuple(saved['clock_pair']), tuple(saved['phase_depths']), saved['horizon']) != ((3,4),(3,4),1024):
        raise ValueError('The inherited clock, phase depths, or observation window changed.')
    cols = tuple(saved['history_columns'])
    ix = {name: cols.index(name) for name in ('step','N3','N4')}
    times = tuple(saved['checkpoints'])
    seeds = tuple(saved['seeds'])
    arms = tuple(saved['arm_order'])
    env_names = ('no_relief','relief_H6')
    depths = (3,4)
    Q = Fraction
    # Fixed for BOTH phase-depth controls. No numerical optimization.
    settings = (Q(0),Q(4,27),-Q(2,27),Q(2,27))
    x = np.array([float(v) for v in settings[:2]])
    y = np.array([float(v) for v in settings[2:]])
    rr = np.arange(3)
    theta = x[:,None,None,None]-y[None,:,None,None] + \
            (rr[None,None,:,None]-rr[None,None,None,:])/3

    def check(condition, message):
        if not condition:
            raise AssertionError(message)

    def stats(a):
        a = np.asarray(a,dtype=float)
        return dict(n=len(a),mean=float(a.mean()),
                    se=float(a.std(ddof=1)/sqrt(len(a))) if len(a)>1 else None)

    # Exact coefficients of the paper's SIGNED CGLMP-3 functional.
    W = np.zeros((2,2,3,3),dtype=int)
    tests = ((0,0,0,1),(1,0,-1,1),(1,1,0,1),(0,1,0,1),
             (0,0,-1,-1),(1,0,0,-1),(1,1,-1,-1),(0,1,1,-1))
    for i,j,d,s in tests:
        for a in range(3):
            for b in range(3):
                if (a-b-d)%3 == 0: W[i,j,a,b] += s
    def score(p): return float(np.sum(W*p))
    deterministic = Counter()
    from itertools import product
    for a0,a1,b0,b1 in product(range(3),repeat=4):
        aa,bb=(a0,a1),(b0,b1)
        deterministic[sum(int(W[i,j,aa[i],bb[j]]) for i in range(2) for j in range(2))]+=1
    check(sum(deterministic.values())==81 and max(deterministic)==2,'Local-bound enumeration')

    def probabilities(c1,c2,phase_table=theta):
        # p=1/9 + (4/27) Re[C1 exp(i2pi theta)] + (2/27) Re[C2 exp(i4pi theta)].
        return np.ones_like(phase_table)/9 + 4*np.real(c1*np.exp(2j*np.pi*phase_table))/27 \
                + 2*np.real(c2*np.exp(4j*np.pi*phase_table))/27

    ideal = probabilities(1,1)
    ideal_score = score(ideal)
    continuous = np.array([0,1/6])[:,None,None,None] - np.array([-1/12,1/12])[None,:,None,None] \
                 + (rr[None,None,:,None]-rr[None,None,None,:])/3
    source_max = 4/3 + 8*sqrt(3)/9
    check(abs(score(probabilities(1,1,continuous))-source_max)<1e-13,'Source continuous Bell regression')
    check(ideal_score>2,'Fixed legal triadic setting must violate the local bound at zero phase')
    check(abs(score(probabilities(0,0)))<1e-14,'Erased-source Bell score')

    # Local projectors independently realize the paper's x-y+(a-b)/3 convention.
    projectors = np.empty((2,2,3,3,9,9),complex)
    for i,j,a,b in product(range(2),range(2),range(3),range(3)):
        av=np.exp(-2j*np.pi*rr*(x[i]+a/3))/sqrt(3)
        bv=np.exp(+2j*np.pi*rr*(y[j]+b/3))/sqrt(3)
        v=np.kron(av,bv)
        projectors[i,j,a,b]=np.outer(v,v.conj())
    for i,j in product(range(2),repeat=2):
        check(np.max(np.abs(projectors[i,j].sum(axis=(0,1))-np.eye(9)))<2e-14,'POVM completeness')

    def density(c1,c2):
        moments={0:1,1:c1,-1:c1.conjugate(),2:c2,-2:c2.conjugate()}
        rho=np.zeros((9,9),complex)
        for r,s in product(range(3),repeat=2):rho[4*r,4*s]=moments[r-s]/3
        return rho
    rho0=density(1+0j,1+0j)
    # Source-backed transfer preparation, BEFORE clock evolution: P01 F3 P01
    # on Alice's factor. This is the manuscript's heralded branch, not a fitted state.
    F3=np.exp(2j*np.pi*rr[:,None]*rr[None,:]/3)/sqrt(3)
    P01=np.diag([1,1,0]);filter_op=P01@F3@P01
    omega=np.zeros(9,complex);omega[[0,4,8]]=1/sqrt(3)
    filtered=np.kron(filter_op,np.eye(3))@omega
    filter_success=float(np.vdot(filtered,filtered).real)
    check(abs(filter_success-4/9)<2e-14,'Source filter success probability')
    filtered/=sqrt(filter_success);rho_filter=np.outer(filtered,filtered.conj())
    reduced=np.trace(rho_filter.reshape(3,3,3,3),axis1=1,axis2=3)
    check(np.max(abs(np.linalg.eigvalsh(reduced)-[0,1/4,3/4]))<2e-14,'Source filter Schmidt spectrum')
    filter_baseline_negativity=sqrt(3)/4
    filter_baseline_bell=score(np.real(np.einsum('ij,xyabji->xyab',rho_filter,projectors)))
    rho_erased=density(0j,0j)
    check(np.max(np.abs(rho0@rho0-rho0))<1e-14,'Pure source')
    b1=4*np.sum(W*np.exp(2j*np.pi*theta))/27
    b2=2*np.sum(W*np.exp(4j*np.pi*theta))/27
    max_errors=Counter()
    counters=Counter(local_assignments=81)
    outcomes={}

    for name in env_names:
        outcomes[name]={}
        for arm in arms:
            inp=saved['environments'][name]['arms'][arm]
            histories=inp['history']
            check(tuple(h['seed'] for h in histories)==seeds,'Do not filter/reorder histories')
            # Read integer counts as integers; the ledger also contains real diagnostics.
            na=np.array([[row[ix['N3']] for row in h['rows']] for h in histories],dtype=np.int64)
            nb=np.array([[row[ix['N4']] for row in h['rows']] for h in histories],dtype=np.int64)
            check(all(tuple(row[ix['step']] for row in h['rows'])==times for h in histories),'Checkpoints changed')
            check(np.all(na>=0) and np.all(nb>=0) and np.all(np.diff(na,axis=1)>=0)
                  and np.all(np.diff(nb,axis=1)>=0),'Invalid cumulative service counts')
            records={}
            for depth in depths:
                qphase=3**depth
                one_depth=[]
                for k,t in enumerate(times):
                    aa=na[:,k]%qphase;bb=nb[:,k]%qphase
                    total=(aa+bb)%qphase
                    histogram=Counter(map(int,total))
                    joint=Counter((int(a),int(b)) for a,b in zip(aa,bb))
                    z=np.exp(2j*np.pi*total/qphase)
                    c1=complex(z.mean());c2=complex((z*z).mean())
                    c1_hist=sum(count*np.exp(2j*np.pi*j/qphase) for j,count in histogram.items())/len(z)
                    c2_hist=sum(count*np.exp(4j*np.pi*j/qphase) for j,count in histogram.items())/len(z)
                    check(abs(c1-c1_hist)<1e-14 and abs(c2-c2_hist)<1e-14,'Histogram moments')
                    rho=density(c1,c2)
                    eig=np.linalg.eigvalsh(rho)
                    check(eig.min()>-5e-14 and abs(np.trace(rho)-1)<1e-14,'Density validity')
                    pt=rho.reshape(3,3,3,3).transpose(0,3,2,1).reshape(9,9)
                    evpt=np.linalg.eigvalsh(pt)
                    neg=(2*abs(c1)+abs(c2))/3
                    check(abs(neg+sum(v for v in evpt if v<0))<2e-14,'Partial-transpose spectrum')
                    # Its spectrum is 1/3 x3, +/-|C1|/3 x2, +/-|C2|/3 x1.
                    expected=np.sort([1/3]*3+[abs(c1)/3,-abs(c1)/3]*2+[abs(c2)/3,-abs(c2)/3])
                    check(np.max(abs(evpt-expected))<2e-14,'Explicit partial-transpose eigenvalues')
                    purity=(3+4*abs(c1)**2+2*abs(c2)**2)/9
                    check(abs(np.trace(rho@rho).real-purity)<2e-14,'Purity formula')
                    marginal=np.trace(rho.reshape(3,3,3,3),axis1=1,axis2=3)
                    check(np.max(np.abs(marginal-np.eye(3)/3))<1e-14,'Local state changed')
                    p=probabilities(c1,c2)
                    dense=np.real(np.einsum('ij,xyabji->xyab',rho,projectors))
                    err=float(np.max(np.abs(dense-p)))
                    max_errors['dense_probabilities']=max(max_errors['dense_probabilities'],err)
                    check(err<4e-14,'Moment/tensor measurement mismatch')
                    check(p.min()>-4e-14 and np.max(abs(p.sum(axis=(2,3))-1))<2e-14,'Outcome normalization')
                    check(np.max(abs(p.sum(axis=3)-1/3))<2e-14 and
                          np.max(abs(p.sum(axis=2)-1/3))<2e-14,'No-signalling marginals')
                    per_history=np.real(b1*z+b2*z*z)
                    I=score(p)
                    check(abs(I-per_history.mean())<2e-14,'Linear Bell average')
                    direct=np.mean((1+2*np.cos(2*np.pi*(theta[None,...]+total[:,None,None,None,None]/qphase)))**2/27,axis=0)
                    max_errors['direct_kernel']=max(max_errors['direct_kernel'],float(np.max(abs(p-direct))))
                    check(np.max(abs(p-direct))<4e-14,'Direct history/kernel mismatch')

                    # A full 9D channel needs JOINT counts, not only their sum.
                    # Every term is a local unitary: trace preservation without a fit.
                    mixture=np.zeros((9,9),complex);completeness=np.zeros(9,float)
                    corrected=np.zeros((9,9),complex)
                    filtered_mixture=np.zeros((9,9),complex)
                    filtered_corrected=np.zeros((9,9),complex)
                    invalid_sum_only=np.zeros((9,9),complex)
                    for (a,b),count in joint.items():
                        wt=count/len(z)
                        u=np.kron(np.exp(2j*np.pi*rr*a/qphase),np.exp(2j*np.pi*rr*b/qphase))
                        r=u[:,None]*rho0*u.conj()[None,:]
                        mixture+=wt*r;completeness+=wt*abs(u)**2
                        corrected+=wt*(u.conj()[:,None]*r*u[None,:])
                        fstate=u[:,None]*rho_filter*u.conj()[None,:]
                        filtered_mixture+=wt*fstate
                        filtered_corrected+=wt*u.conj()[:,None]*fstate*u[None,:]
                        bad_u=np.kron(np.exp(2j*np.pi*rr*(a+b)/qphase),np.ones(3))
                        invalid_sum_only+=wt*bad_u[:,None]*rho_filter*bad_u.conj()[None,:]
                    check(np.max(abs(completeness-1))<2e-14,'Random-unitary trace preservation')
                    check(np.max(abs(mixture-rho))<3e-14,'Full local phase channel mismatch')
                    check(np.max(abs(corrected-rho0))<3e-14,'All-history local inverse phase control')
                    max_errors['channel_reconstruction']=max(max_errors['channel_reconstruction'],float(np.max(abs(mixture-rho))))
                    max_errors['inverse_phase']=max(max_errors['inverse_phase'],float(np.max(abs(corrected-rho0))))
                    # Transfer needs the SAME local channel and its JOINT service counts.
                    feig=np.linalg.eigvalsh(filtered_mixture)
                    check(feig.min()>-5e-14 and abs(np.trace(filtered_mixture)-1)<2e-14,'Filtered-state validity')
                    fpt=filtered_mixture.reshape(3,3,3,3).transpose(0,3,2,1).reshape(9,9)
                    fev=np.linalg.eigvalsh(fpt);fneg=float(-fev[fev<0].sum())
                    fp=np.real(np.einsum('ij,xyabji->xyab',filtered_mixture,projectors))
                    check(np.max(abs(fp.sum(axis=(2,3))-1))<3e-14 and fp.min()>-3e-14,'Filtered probabilities')
                    check(np.max(abs(fp.sum(axis=3)[:,0,:]-fp.sum(axis=3)[:,1,:]))<2e-14 and
                          np.max(abs(fp.sum(axis=2)[0,:,:]-fp.sum(axis=2)[1,:,:]))<2e-14,'Filtered no-signalling')
                    check(np.max(abs(filtered_corrected-rho_filter))<3e-14,'Filtered all-history inverse phases')
                    wrong=float(np.sum(abs(np.linalg.eigvalsh(filtered_mixture-invalid_sum_only)))/2)
                    transfer=dict(operator='(P01 F3 P01) on Alice BEFORE clock evolution',
                        herald_success_probability=filter_success,baseline_negativity=filter_baseline_negativity,
                        baseline_bell=filter_baseline_bell,density_matrix=filtered_mixture,
                        joint_probabilities=fp,bell=score(fp),negativity=fneg,
                        fraction_of_initial_negativity=fneg/filter_baseline_negativity,
                        invalid_sum_only_trace_distance=wrong,
                        local_inverse_control_max_error=float(np.max(abs(filtered_corrected-rho_filter))))
                    counters['filtered_transfer_density_matrices']+=1
                    counters['filtered_joint_probabilities']+=36
                    # Mean-phase surrogate is PURE, and NOT the adopted ensemble law.
                    meanphase=float((na[:,k]+nb[:,k]).mean()/qphase)
                    pure_mean=score(probabilities(np.exp(2j*np.pi*meanphase),np.exp(4j*np.pi*meanphase)))
                    # Monte Carlo jackknife spread; finite-ensemble negativity is nonlinear/biased.
                    l1=(len(z)*c1-z)/(len(z)-1);l2=(len(z)*c2-z*z)/(len(z)-1)
                    loo=(2*abs(l1)+abs(l2))/3
                    negse=float(np.sqrt((len(z)-1)/len(z)*np.sum((loo-loo.mean())**2)))
                    one_depth.append(dict(step=t,phase_modulus=qphase,history_count=len(z),
                        sum_count_histogram={str(j):c for j,c in sorted(histogram.items())},
                        joint_count_histogram=tuple((a,b,c) for (a,b),c in sorted(joint.items())),
                        C1=c1,C2=c2,density_matrix=rho,joint_probabilities=p,
                        bell=stats(per_history),source_upper_bound_exceeded=I>2,
                        negativity=float(neg),negativity_jackknife_SE=negse,
                        density_purity=float(purity),partial_transpose_eigenvalues=evpt,
                        mean_phase_pure_surrogate_bell=pure_mean,
                        phase_corrected_bell=ideal_score,phase_corrected_negativity=1.0,
                        erased_bell=0.0,erased_negativity=0.0,
                        individual_conditioned_state_negativity=1.0,
                        source_filtered_transfer=transfer,
                        seed_block_bell=tuple(float(per_history[i:i+256].mean()) for i in (0,256))))
                    counters['density_matrices']+=1
                    counters['joint_setting_probabilities']+=36
                    counters['saved_history_snapshots_analyzed']+=len(z)
                    counters['joint_histogram_channel_terms']+=len(joint)
                records[depth]=tuple(one_depth)
            outcomes[name][arm]=dict(readout_by_depth=records,original_recording_cost=inp['event']['cost'])
    check(saved==before,'An input changed')
    # Reuse the one-step hypergeometric law only for a ONE-step exact check.
    # The changing-history channel is NOT approximated by multiplying unconditional laws.
    one_steps=[]
    for name in env_names:
        for arm in arms:
            e=saved['environments'][name]['arms'][arm]['exact_next_step']
            pool=int(e['pool']);p=Q(3,pool);p2=Q(6,pool*(pool-1))
            weights=(1-2*p+p2,2*(p-p2),p2)
            check(sum(weights)==1 and all(v>=0 for v in weights),'One-step sampling normalization')
            moment=[]
            for h in (1,2):
                zz=np.exp(2j*np.pi*h/27)
                polynomial=1+2*float(p)*(zz-1)+float(p2)*(zz-1)**2
                direct=sum(float(v)*zz**j for j,v in enumerate(weights))
                check(abs(direct-polynomial)<2e-14,'One-step characteristic polynomial')
                moment.append(complex(direct))
            one_steps.append(dict(environment=name,arm=arm,pool=pool,increment_law=weights,
                exact_next_state_moments=tuple(moment),exact_next_bell=score(probabilities(*moment))))
    print('CELL 32 — A SMALLER VERTICAL BRIDGE: WORKLOAD -> TIMING -> QUANTUM CORRELATIONS')
    print('Saved Cell29 native histories; NO new primary native simulation or observer-identity rule.')
    print('Inherited source Omega_2 and local phase advance; classical service is NOT the source of entanglement.')
    print('Fixed CGLMP3 settings: (0, 4/27, -2/27, 2/27), unchanged across all arms/depths.')
    print(f'Pure-source score = {ideal_score:.12f}; local upper bound = 2; continuous source benchmark = {source_max:.12f}.')
    print('\n AT 128 UPDATES, PHASE INCREMENT 1/27:')
    print(' environment / arm                 CGLMP I3 +/- history MC SE    negativity +/- jackknife SE')
    for name in env_names:
        for arm in arms:
            z=outcomes[name][arm]['readout_by_depth'][3][times.index(128)]
            print(f" {name:12s}/{arm:15s} {z['bell']['mean']: .9f} +/- {z['bell']['se']:.9f}"
                  f"      {z['negativity']:.9f} +/- {z['negativity_jackknife_SE']:.9f}")
    print('\n AT 1024 UPDATES, PHASE INCREMENT 1/27:')
    for name in env_names:
        for arm in arms:
            z=outcomes[name][arm]['readout_by_depth'][3][-1]
            print(f" {name:12s}/{arm:15s} I3={z['bell']['mean']:+.9f}; negativity={z['negativity']:.9f}.")
    print('\nSOURCE-BACKED FILTER TRANSFER: P01 F3 P01 before timing; success 4/9; Schmidt squares (3/4,1/4,0).')
    for name in env_names:
        a=outcomes[name]['existing_pair']['readout_by_depth'][3][-1]['source_filtered_transfer']
        print(f" {name}: filtered negativity at 1024 = {a['negativity']:.9f}, "
              f"fraction retained = {a['fraction_of_initial_negativity']:.9f}; "
              f"wrong total-count-only channel trace distance = {a['invalid_sum_only_trace_distance']:.9f}.")
    print('Filtered-source transfer imports the declared heralded instrument operation, not a new native workload.')
    print('A failed FIXED Bell witness is not a separability or all-settings-locality verdict.')
    print('Complete erasure: I3=0, negativity=0. Log-conditioned local inverse phases: original source restored.')
    print('Inverse-phase control retains ALL histories and all measurements; it is not postselection.')
    print('No physical gate/feedback cost, energy, seconds, propagation, GR redshift, or neutron operation added.')
    print('The 1/81 control, every checkpoint and seed block, and full joint probability tables remain saved.')
    print('PASS: local-bound enumeration, source regression, density/partial-transpose spectra,')
    print('      local Kraus completeness, dense tensor outcomes, no-signalling, exact one-step laws, unchanged input.')
    return dict(protocol='Cell32 v1 fixed-source Bell response to all saved Cell29 service histories',
        settings=settings,primary_phase_depth=3,phase_depths=depths,checkpoints=times,seeds=seeds,
        ideal_bell=ideal_score,continuous_source_bell=source_max,
        filtered_preparation=dict(success_probability=filter_success,Schmidt_squares=(3/4,1/4,0),
                                  negativity=filter_baseline_negativity,applied_before_timing=True),local_deterministic_score_counts=dict(deterministic),
        bell_moment_coefficients=(complex(b1),complex(b2)),environments=outcomes,
        exact_initial_channels=one_steps,validation_counts=dict(counters),max_numerical_errors=dict(max_errors),
        new_primary_native_updates=0,reused_native_histories=6*len(seeds),
        source_instrument_presupposed=True,phase_rule_changed=False,observer_identity_needed=False,
        likelihood_fit_performed=False,measurement_settings_optimized=False,history_postselection=False,
        gate_feedback_cost_derived=False,local_time_dilation_derived=False,irreversible_decoherence_derived=False,
        neutron_search_modified=False,physical_calibration=None)


dcu_mass_32 = _run_dcu_mass_32()
