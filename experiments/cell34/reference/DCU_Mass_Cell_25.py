# CELL 25 — first vertical pilot: native maintenance -> finite phase -> joint readout.
# Paste after Cell 24. NumPy + standard library. No network, physical constants or fits.
# ADDED coupling: servicing either of two tagged siblings applies D(1/3^d) on its
# attached factor (opposite orientation in the differential control). Other native
# events leave that attached phase state unchanged. This is a passive readout/drive
# hypothesis, NOT a native energy law or a derivation of the quantum instrument.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations
from math import comb, cos, pi, sqrt
from random import Random


def _run_dcu_mass_25():
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError('Cell 25 requires NumPy in this notebook kernel.') from exc
    missing = [k for k in ('_Native', 'dcu_mass_24') if k not in globals()]
    if missing:
        raise RuntimeError('Run Cell 24 first. Missing: ' + ', '.join(missing))
    prior = dcu_mass_24
    snapshot = deepcopy(prior)
    Q = Fraction
    if (prior['local_phase_labels'] != (0, 1, 2) or
        prior['local_probability_harmonics'] != dict(constant=Q(1,3), cos1=Q(4,9), cos2=Q(2,9))):
        raise ValueError('The inherited phase/readout convention changed.')
    readout = prior['helpers']['sibling_readout']
    kernel = prior['helpers']['kernel']

    # Fix the domain before generating histories. H changes the EXISTING relief
    # rule, not a fitted quantum damping rate. Both have Gamma=3, m=0.
    steps = 4096
    checkpoints = tuple(range(0, steps+1, 64))
    environments = (('no_relief', 0), ('relief_H6', 6))
    seeds = tuple(20260920+i for i in range(128)) + tuple(20270920+i for i in range(128))
    depths = (3, 4)  # Apparatus phase resolutions, NOT spatial dimension or DAG grade.
    primary_depth = 3
    carriers = (3, 4)
    gamma = 3

    def prepare():
        # Exactly the first TWO deterministic updates of Gamma=3,m=0 genesis.
        # First q=2 services a,b; then q=3 services a,b,{a,b}.
        st = _Native()
        a, cost0 = st.add_batch([(0, 1)])
        assert a == [2] and cost0 == 0
        b, cost1 = st.add_batch([(0, 2), (1, 2)])
        assert b == [3, 4] and cost1 == 26
        fx = prior['siblings'][2]['fixture']
        assert tuple(st.parents) == fx['parentage'] and tuple(b) == fx['siblings']
        # Work is NOT reset/forgiven: carry 26 outstanding requests into observation.
        return st, cost1, 0

    def step(st, forced, relief, rng, H):
        n, total = len(st), len(st)+forced
        q = min(gamma, total)
        drawn = rng.sample(range(total), q)  # EXACT uniform sample of individual tokens.
        served = tuple(sorted(j for j in drawn if j < n))
        f = q-len(served)
        Fminus, Pminus = forced-f, relief+2*f
        quota = 2*((max(1, Pminus//6)+1)//2)
        v = min(quota, H, Fminus, Pminus) if Fminus >= gamma and Pminus >= 6 else 0
        events = [e for e in combinations(served, 2) if e not in st.pair_to_id]
        before_recorded = frozenset(st.recorded)
        hits = Counter(z for a,b in events for z in (st.ancestors[a] | st.ancestors[b]))
        independent_cost = sum((2*h + (9 if z not in before_recorded else 0))*st.path_weight(z)
                               for z,h in hits.items())
        born, cost = st.add_batch(events)
        assert cost == independent_cost and len(born) == len(events)
        assert all(e in st.pair_to_id for e in combinations(served, 2))
        afterF, afterP = Fminus-v+cost, Pminus-v
        assert afterF >= 0 and afterP >= 0
        return afterF, afterP, dict(n=n, F=forced, P=relief, total=total, q=q,
            served=served, forced_served=f, relief_removed=v, new_pairs=tuple(events),
            born=tuple(born), cost=cost)

    def mean_se(values):
        a = np.asarray(values, dtype=float)
        return dict(mean=float(a.mean()), se=float(a.std(ddof=1)/sqrt(len(a))))

    def run(seed, H, trace=False):
        st, forced, relief = prepare()
        rng = Random(seed)
        count_a = count_b = tau = sumf = sumv = sumcost = 0
        expected_count = variance_sum = variance_diff = 0.0
        per_object_ticks = Q(0)
        variance_given_service_count = 0.0
        saved, full = [], []
        # row order is public below; 'tau' is elapsed maintenance AFTER preparation.
        def save(t):
            saved.append((t, len(st), forced, relief, tau, count_a, count_b,
                          expected_count, variance_sum, variance_diff, sumcost, sumf, sumv,
                          per_object_ticks, variance_given_service_count))
        save(0)
        for t in range(1, steps+1):
            forced, relief, e = step(st, forced, relief, rng, H)
            p = e['q']/e['total']
            p2 = e['q']*(e['q']-1)/(e['total']*(e['total']-1))
            expected_count += 2*p
            variance_sum += 2*p+2*p2-4*p*p
            variance_diff += 2*(p-p2)
            count_a += int(carriers[0] in e['served'])
            count_b += int(carriers[1] in e['served'])
            selected = len(e['served'])
            tau += selected
            # Derived clock comparison, NOT a replacement for the primary drive.
            # Conditional on n,s: E[K_increment]=2s/n (hypergeometric identities).
            if selected:
                per_object_ticks += Q(selected, e['n'])
            variance_given_service_count += (selected*(2/e['n'])*(1-2/e['n'])
                                              *(e['n']-selected)/(e['n']-1))
            sumf += e['forced_served']
            sumv += e['relief_removed']
            sumcost += e['cost']
            assert forced == 26+sumcost-sumf-sumv
            assert relief == 2*sumf-sumv
            assert tau+sumf == gamma*t and count_a+count_b <= tau
            if trace:
                full.append(dict(e, observation_step=t, F_after=forced, P_after=relief,
                    n_after=len(st), tau=tau, count_a=count_a, count_b=count_b))
            if t % 64 == 0:
                save(t)
        assert st.recorded == {v for pair in st.parents[2:] for v in pair}
        return dict(seed=seed, checkpoints=tuple(saved), trace=tuple(full),
                    final_parents=tuple(st.parents) if trace else None)

    # EXACT one-service law: neither token identity nor graph distance enters it.
    st0, F0, P0 = prepare()
    n0, total0 = len(st0), len(st0)+F0
    p, p2 = Q(gamma,total0), Q(gamma*(gamma-1),total0*(total0-1))
    plus_law = {0:1-2*p+p2, 1:2*(p-p2), 2:p2}
    minus_law = {-1:p-p2, 0:1-2*(p-p2), 1:p-p2}
    draws = comb(total0, gamma)
    empirical_exact = {pair: (Counter(),Counter()) for pair in combinations(range(n0),2)}
    for draw in combinations(range(total0),gamma):
        chosen=set(draw)
        for (a,b),(plus,minus) in empirical_exact.items():
            ia,ib=int(a in chosen),int(b in chosen)
            plus[ia+ib]+=1
            minus[ia-ib]+=1
    assert all({k:Q(v,draws) for k,v in a.items()} == plus_law and
               {k:Q(v,draws) for k,v in b.items()} == minus_law
               for a,b in empirical_exact.values())
    graph=[set() for _ in st0.parents]
    for child,pair in enumerate(st0.parents[2:],2):
        for par in pair:
            graph[child].add(par); graph[par].add(child)
    distances={}
    for a in range(n0):
        dist={a:0}; queue=[a]
        for v in queue:
            for w in sorted(graph[v]):
                if w not in dist: dist[w]=dist[v]+1; queue.append(w)
        for b in range(a+1,n0): distances[a,b]=dist[b]
    assert set(distances.values()) == {1,2}
    single_step = dict(n=n0,F=F0,P=P0,q=gamma,p=p,pair_probability=p2,
        sum_increment_law=plus_law,difference_increment_law=minus_law,
        exact_draws=draws,pairs_checked=len(empirical_exact),parent_child_distances=distances,
        same_law_for_all_pairs=True,
        tagged_total_mean=2*p,tagged_total_variance=2*p+2*p2-4*p*p,
        tagged_difference_variance=2*(p-p2),two_token_covariance=p2-p*p)

    ensembles={}
    for name,H in environments:
        histories=tuple(run(seed,H,trace=(index==0)) for index,seed in enumerate(seeds))
        ensembles[name]=dict(Gamma=gamma,m=0,H=H,histories=histories)

    # Predictions given by histories, not a fitted damping curve. Terminal
    # probabilities at different checkpoints are COUNTERFACTUAL separate readouts;
    # no repeated projective measurement is inserted during the evolution.
    curves={}; summaries={}; max_tensor_error=0.0; max_rephasing_error=0.0
    for name,H in environments:
        hist=ensembles[name]['histories']
        arrays=np.asarray([r['checkpoints'] for r in hist],dtype=float)
        a=arrays[:,:,5].astype(np.int64); b=arrays[:,:,6].astype(np.int64)
        curves[name]={}
        for d in depths:
            denominator=3**d
            for mode, signed in (('sum',a+b),('difference',a-b)):
                residue=signed % denominator
                z=np.exp(2j*pi*residue/denominator)
                C1=z.mean(axis=0); C2=(z*z).mean(axis=0)
                prob=np.array([1/3+(4/9)*np.real(C1*np.exp(-2j*pi*s/3))+
                               (2/9)*np.real(C2*np.exp(-4j*pi*s/3)) for s in range(3)]).T
                direct=(1+2*np.cos(2*pi*(residue[:,:,None]/denominator-np.arange(3)/3)))**2/9
                assert np.max(np.abs(direct.mean(axis=0)-prob)) < 4e-14
                assert np.max(np.abs(prob.sum(axis=1)-1)) < 4e-14 and prob.min() > -2e-14
                # A pure trajectory has |z|=1. Loss of ENSEMBLE coherence is history averaging.
                assert np.max(np.abs(np.abs(z)-1)) < 3e-15
                purity=1/3+(4/9)*np.abs(C1)**2+(2/9)*np.abs(C2)**2
                assert purity.min() >= 1/3-2e-14 and purity.max() <= 1+2e-14
                # Deliberately WRONG mean-clock substitution, retained as a diagnostic.
                mean_phase=signed.mean(axis=0)/denominator
                mean_clock=(1+2*np.cos(2*pi*(mean_phase % 1)))**2/9
                key=f'd{d}_{mode}'
                curves[name][key]=dict(C1=tuple(map(complex,C1)),C2=tuple(map(complex,C2)),
                    probabilities=tuple(tuple(map(float,row)) for row in prob),
                    purity=tuple(map(float,purity)),P0_at_mean_phase=tuple(map(float,mean_clock)),
                    final_residue_counts=dict(Counter(map(int,residue[:,-1]))),
                    final_P0=mean_se(direct[:,-1,0]),
                    seed_block_P0=(mean_se(direct[:128,-1,0]),mean_se(direct[128:,-1,0])))
                # Compare against Cell24's actual tensor circuit at both resolutions.
                for k in (0,64,128,192):
                    aa,bb=int(a[k,-1]),int(b[k,-1])
                    xs=(Q(aa,denominator),Q(bb if mode=='sum' else -bb,denominator))
                    observed=readout(xs)['probabilities']
                    max_tensor_error=max(max_tensor_error,float(np.max(np.abs(np.array(observed)-direct[k,-1]))))
                    # History-aware inverse settings undo phases in the attached model.
                    phases0=prior['helpers']['phases']((0,1,2),xs[0])
                    phases1=prior['helpers']['phases']((0,1,2),xs[1])
                    initial=np.eye(3,dtype=complex)/sqrt(3)  # coefficients of |rr>, not a density matrix.
                    evolved=phases0[:,None]*initial*phases1[None,:]
                    inverse0=prior['helpers']['phases']((0,1,2),-xs[0])
                    inverse1=prior['helpers']['phases']((0,1,2),-xs[1])
                    restored=inverse0[:,None]*evolved*inverse1[None,:]
                    max_rephasing_error=max(max_rephasing_error,float(np.max(np.abs(restored-initial))))
                    assert Q(aa,denominator) == Q(3*aa,3*denominator)
        K=a[:,-1]+b[:,-1]; D=a[:,-1]-b[:,-1]
        compensator=arrays[:,-1,7]
        residual=K-compensator
        stat=mean_se(residual)
        summaries[name]=dict(
            final_n=mean_se(arrays[:,-1,1]),final_F=mean_se(arrays[:,-1,2]),
            maintenance_ticks=mean_se(arrays[:,-1,4]),tagged_services=mean_se(K),
            expected_tagged_services=mean_se(compensator),
            martingale_total_residual=stat,
            martingale_difference_residual=mean_se(D),
            per_object_ticks=mean_se(arrays[:,-1,13]),
            residual_from_actual_per_object_clock=mean_se(K-2*arrays[:,-1,13]),
            empirical_per_object_residual_second_moment=float(np.mean((K-2*arrays[:,-1,13])**2)),
            expected_per_object_quadratic_variation=float(arrays[:,-1,14].mean()),
            total_residual_mean_z=stat['mean']/stat['se'] if stat['se'] else None,
            empirical_total_residual_second_moment=float(np.mean(residual**2)),
            expected_total_quadratic_variation=float(arrays[:,-1,8].mean()),
            empirical_difference_second_moment=float(np.mean(D**2)),
            expected_difference_quadratic_variation=float(arrays[:,-1,9].mean()),
            state_ranges=dict(n=(int(arrays[:,-1,1].min()),int(arrays[:,-1,1].max())),
                              F=(int(arrays[:,-1,2].min()),int(arrays[:,-1,2].max()))))
    assert max_tensor_error < 6e-14 and max_rephasing_error < 6e-14

    # Same-draw passive-readout countercontrol: independently repeat one complete
    # native history per environment; it has NO phase state inside the native step.
    for name,H in environments:
        replay=run(seeds[0],H,trace=True)
        assert replay == ensembles[name]['histories'][0]
    assert prior == snapshot
    total_steps=len(environments)*len(seeds)*steps
    print('CELL 25 — VERTICAL SERVICE -> PHASE -> INTERFERENCE PILOT')
    print('Native: Gamma=3,m=0; H=0 or H=6; uniform individual-token service; all enabled births.')
    print('Prepared from TWO actual deterministic genesis updates: n=5,F=26,P=0,tau=5; siblings3,4.')
    print('Observation starts after preparation; its backlog is retained. All quoted ticks are subsequent.')
    print('ADDED: a served tagged sibling applies D(1/3^d); other events preserve the attached phase state.')
    print('d=3 primary; d=4 finer NEW SETTING in the model, not relabeling the same coarse tick.')
    print('No phase setting/terminal readout cost is supplied; counting history is not modified.')
    print(f'{len(seeds)} histories/environment, {steps} steps each; {total_steps:,} native updates retained.')
    print('\nExact initial tagged-service law (K=number of selected siblings):',plus_law)
    print('Exact initial difference law:',minus_law)
    print(f'{draws} token subsets x {len(empirical_exact)} pairs: identical laws at PC distances 1 and 2.')
    print('\n environment        mean n    mean F    mean ticks    mean tagged K    predicted mean K')
    for name,_ in environments:
        r=summaries[name]
        print(f" {name:17s} {r['final_n']['mean']:8.3f} {r['final_F']['mean']:9.3f}"
              f" {r['maintenance_ticks']['mean']:13.3f} {r['tagged_services']['mean']:16.3f}"
              f" {r['expected_tagged_services']['mean']:19.3f}")
    print('\n environment / clock    |C1|        |C2|       purity       mean P0 +/- service-MC SE')
    for name,_ in environments:
        for key,row in curves[name].items():
            pp=row['final_P0']
            print(f" {name:12s}/{key:13s} {abs(row['C1'][-1]):10.6f} {abs(row['C2'][-1]):10.6f}"
                  f" {row['purity'][-1]:11.6f} {pp['mean']:12.6f} +/- {pp['se']:.6f}")
    print('\nMean-phase substitution is NOT used; for primary sum clock:')
    for name,_ in environments:
        row=curves[name]['d3_sum']
        print(f"  {name}: E[P0]={row['final_P0']['mean']:.6f}; P0(E[phase])={row['P0_at_mean_phase'][-1]:.6f}.")
    print('Derived per-object activity clock chi=sum(s/n), NOT a new phase rule or SI calibration:')
    for name,_ in environments:
        row=summaries[name]
        print(f"  {name}: mean chi={row['per_object_ticks']['mean']:.6f}; "
              f"mean(K-2chi)={row['residual_from_actual_per_object_clock']['mean']:+.6f} "
              f"+/- {row['residual_from_actual_per_object_clock']['se']:.6f} MC SE.")
    print('Every retained history is coherent in this attached model; logged inverse phases restore P0=1.')
    print('Mixed ensemble states reflect unobserved service histories, NOT proved irreversible decoherence.')
    print('One-step equal-resolution clock noise contains no PC-distance/sector distinction; no 3D claim.')
    print('Finite-run Monte Carlo discrepancies are retained; no samples, coefficients, or phases fitted.')
    print('PASS: exact one-step enumeration, full native ledgers/bursts, tensor readout, refinement,')
    print('      pure-history checks, two exact replays, and preservation of Cell24.')
    print('No SI calibration, energy/mass prediction, physical Hamiltonian selection, or native locality added.')
    return dict(protocol='maintenance-driven passive finite-phase interferometer',
        Gamma=gamma,m=0,steps=steps,seeds=seeds,checkpoints=checkpoints,carriers=carriers,
        depths=depths,primary_depth=primary_depth,
        checkpoint_columns=('t','n','F','P','tau','count_a','count_b','expected_K',
                            'predictable_variance_K','predictable_variance_D','cost','forced_served','relief_removed',
                            'per_object_ticks','predictable_variance_K_given_s'),
        preparation=dict(parentage=tuple(st0.parents),n=5,F=26,P=0,previous_tau=5,
                         actual_genesis_updates=2,work_discarded=0),
        one_step=single_step,ensembles=ensembles,summaries=summaries,curves=curves,
        total_native_updates=total_steps,max_tensor_error=max_tensor_error,
        max_rephasing_error=max_rephasing_error,
        helpers=dict(prepare=prepare,native_step=step,run_history=run),
        stochastic_replicates_not_empirical_tests=True,
        phase_drive_is_additional_hypothesis=True,phase_retention_is_additional_hypothesis=True,
        phase_gate_recording_cost=None,feedback_into_counting=False,
        dimension_interpreted_as_spatial=False,physical_calibration=None,
        new_mass_formula=False,physical_energy_generator_selected=False)


dcu_mass_25 = _run_dcu_mass_25()
