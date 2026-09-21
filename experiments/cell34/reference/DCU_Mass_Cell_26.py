# CELL 26 — record-supported clocks: an overlap-sensitive response, not physical space.
# Paste after Cell 25. NumPy + standard library; no downloads or fitted parameters.
# ADDED readout: each served register in a fixed placement path adds the SAME 1/27
# phase unit. Two such counters drive the ORIGINAL two-sibling apparatus with +/-
# orientations. Shared-register increments cancel in its differential phase.
# This fan-out/control wiring is a proposed passive readout, NOT a new native gate,
# work charge, propagation law, quark classifier, or identification of path length with space.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations, product
from math import comb, sqrt
from random import Random


def _run_dcu_mass_26():
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError('Cell 26 needs NumPy in this notebook kernel.') from exc
    missing = [k for k in ('_Native', 'dcu_mass_24', 'dcu_mass_25') if k not in globals()]
    if missing:
        raise RuntimeError('Run Cell 25 first. Missing: ' + ', '.join(missing))
    old, instrument = dcu_mass_25, dcu_mass_24
    before = deepcopy((old, instrument))
    Q = Fraction
    if (old['Gamma'], old['m'], old['carriers'], old['primary_depth']) != (3, 0, (3, 4), 3):
        raise ValueError('Prior native/phase convention changed; no replacement chosen.')
    if instrument['local_phase_labels'] != (0, 1, 2):
        raise ValueError('Inherited instrument changed.')
    native_step = old['helpers']['native_step']
    readout = instrument['helpers']['sibling_readout']

    # Fixed scope: two previously saved first-seed histories, PAST-ONLY state at 512.
    # Enumerate ALL length-3 recorded paths (primary) and length-4 paths (transfer).
    # This is a finite retrospective test, not a preregistered particle/epoch search.
    freeze_step, horizon = 512, 2048
    record_lengths, phase_depths = (3, 4), (3, 4)
    checkpoints = tuple(range(0, horizon + 1, 128))
    seeds = tuple(20280920 + i for i in range(256))
    environments = (('no_relief', 0), ('relief_H6', 6))

    def exact_rank(matrix):
        a = [list(map(Q, row)) for row in matrix]
        if not a:
            return 0
        r = 0
        for c in range(len(a[0])):
            pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
            if pivot is None:
                continue
            a[r], a[pivot] = a[pivot], a[r]
            lead = a[r][c]
            a[r] = [v/lead for v in a[r]]
            for i in range(r+1, len(a)):
                lead = a[i][c]
                if lead:
                    a[i] = [v-lead*w for v, w in zip(a[i], a[r])]
            r += 1
            if r == len(a):
                break
        return r

    def restored_checkpoint(name, H):
        saved = old['ensembles'][name]['histories'][0]
        st, F, P = old['helpers']['prepare']()
        rng = Random(saved['seed'])
        tau = 0
        for t in range(freeze_step):
            F, P, event = native_step(st, F, P, rng, H)
            expected = saved['trace'][t]
            assert all(event[k] == expected[k] for k in event)
            assert (F, P, len(st)) == (expected['F_after'], expected['P_after'], expected['n_after'])
            tau += len(event['served'])
        next_event = saved['trace'][freeze_step]
        assert (len(st), F, P) == (next_event['n'], next_event['F'], next_event['P'])
        assert tuple(st.parents) == saved['final_parents'][:len(st)]
        return st, F, P, dict(seed=saved['seed'], after_observation_step=freeze_step,
                             n=len(st), F=F, P=P, elapsed_maintenance=tau,
                             genesis_maintenance=tau + old['preparation']['previous_tau'])

    def path_panels(st):
        paths = {0: [(0, ())], 1: [(1, ())]}
        for z, parents in enumerate(st.parents[2:], 2):
            paths[z] = [(root, tokens+((z, entry),))
                        for entry, parent in enumerate(parents)
                        for root, tokens in paths[parent]
                        if len(tokens) < max(record_lengths)]
        out = {}
        for length in record_lengths:
            selected = sorted((root, tokens) for z in sorted(st.recorded) if z >= 2
                              for root, tokens in paths[z] if len(tokens) == length)
            if len(selected) < 3:
                raise ValueError('Declared checkpoint has too few paths; no later time substituted.')
            supports = [tuple(z for z, _ in tokens) for _, tokens in selected]
            incidence = np.zeros((len(selected), len(st)), dtype=np.int64)
            for i, (root, tokens) in enumerate(selected):
                entered = root
                for z, entry in tokens:
                    assert st.parents[z][entry] == entered and z in st.recorded
                    entered = z
                assert len(set(supports[i])) == length
                incidence[i, list(supports[i])] = 1
            overlap = incidence @ incidence.T
            squared_distance = 2*length - 2*overlap
            assert np.all(squared_distance >= 0) and np.all(np.diag(squared_distance) == 0)
            pairs = tuple(combinations(range(len(selected)), 2))
            grouped = {}
            for ij in pairs:
                grouped.setdefault(int(squared_distance[ij]), []).append(ij)
            affine_rank = exact_rank((incidence[1:]-incidence[0]).tolist())
            # Two explicit information-loss/geometry diagnostics, selected lexicographically.
            alias_pair = next((ij for ij in pairs if squared_distance[ij] == 0), None)
            alias = None
            if alias_pair is not None:
                i, j = alias_pair
                a, b = selected[i][1], selected[j][1]
                registers = sorted(set(supports[i]))
                agree = 0
                for values in product(range(3), repeat=len(registers)):
                    g = dict(zip(registers, values))
                    def address(tokens):
                        return tuple((1-g[z] if g[z] in (0, 1) else 2) if entry else g[z]
                                     for z, entry in tokens)
                    agree += address(a) == address(b)
                alias = dict(pair=alias_pair, paths=(selected[i], selected[j]),
                             registers=tuple(registers), same_clock_pathwise=True,
                             address_agreement_probability=Q(agree, 3**len(registers)))
            ultrametric_counterexample = None
            for i, j, k in combinations(range(len(selected)), 3):
                dij, dik, djk = (int(squared_distance[i,j]), int(squared_distance[i,k]),
                                  int(squared_distance[j,k]))
                distances = sorted((dij, dik, djk))
                if distances[2] > distances[1]:
                    ultrametric_counterexample = dict(indices=(i,j,k),
                        squared_distances=(dij,dik,djk), paths=tuple(selected[z] for z in (i,j,k)))
                    break
            out[length] = dict(paths=tuple(selected),supports=tuple(supports),incidence=incidence,
                overlap=overlap,squared_distance=squared_distance,groups=grouped,
                unique_supports=len(set(supports)),affine_dimension=affine_rank,
                alias=alias,ultrametric_counterexample=ultrametric_counterexample)
        return out

    def choose(n, k):
        return comb(n, k) if 0 <= k <= n else 0

    def difference_law(total, q, unique_each):
        # Equal-length path difference: +1 on A\B, -1 on B\A, 0 everywhere else.
        m = unique_each
        if total < 2*m:
            raise ValueError('Support exceeds pool.')
        counts = Counter()
        for a in range(min(q, m)+1):
            for b in range(min(q-a, m)+1):
                counts[a-b] += choose(m,a)*choose(m,b)*choose(total-2*m,q-a-b)
        law = {v:Q(c,choose(total,q)) for v,c in counts.items() if c}
        assert sum(law.values()) == 1 and sum(Q(v)*p for v,p in law.items()) == 0
        return law

    def exact_checkpoint_audit(st, F, panels):
        n, total, q = len(st), len(st)+F, min(3,len(st)+F)
        p, p2 = Q(q,total), Q(q*(q-1),total*(total-1))
        a, b = p-p2, p2-p*p
        # Independently enumerate maintenance SUBSETS, with each forced-token
        # combination counted exactly. No enormous explicit forced-pool enumeration.
        representatives=[]
        for length,panel in panels.items():
            for delta,pairs in sorted(panel['groups'].items()):
                i,j=pairs[0]
                representatives.append((length,delta,i,j))
        counters=[Counter() for _ in representatives]
        raw_moments=[[0,0,0,0,0] for _ in representatives] # x,y,x^2,y^2,xy
        subsets=0
        for s in range(min(q,n)+1):
            multiplicity=choose(F,q-s)
            if not multiplicity:
                continue
            for served in combinations(range(n),s):
                subsets += 1
                W=set(served)
                for k,(length,delta,i,j) in enumerate(representatives):
                    panel=panels[length]
                    x,y=(len(W & set(panel['supports'][z])) for z in (i,j))
                    counters[k][x-y] += multiplicity
                    for z,value in enumerate((x,y,x*x,y*y,x*y)):
                        raw_moments[k][z] += multiplicity*value
        results=[]
        for k,(length,delta,i,j) in enumerate(representatives):
            distribution={v:Q(c,choose(total,q)) for v,c in counters[k].items()}
            compact=difference_law(total,q,delta//2)
            assert distribution == compact
            ex,ey,ex2,ey2,exy=[Q(z,choose(total,q)) for z in raw_moments[k]]
            shared=length-delta//2
            assert ex == ey == length*p
            assert exy-ex*ey == a*shared+b*length**2
            assert ex2+ey2-2*exy == a*delta
            results.append(dict(record_length=length,pair=(i,j),symmetric_difference=delta,
                shared_registers=shared,difference_distribution=distribution,
                conditional_covariance=exy-ex*ey,conditional_difference_variance=a*delta))
        return dict(n=n,F=F,total=total,q=q,p=p,p2=p2,activity_factor=a,
                    common_sampling_covariance=b,maintenance_subsets=subsets,tests=results)

    def mean_se(x):
        x=np.asarray(x,dtype=float)
        return dict(mean=float(x.mean()),se=float(x.std(ddof=1)/sqrt(len(x))))

    results={}
    tensor_error=0.0
    for name,H in environments:
        base,F0,P0,checkpoint=restored_checkpoint(name,H)
        panels=path_panels(base)
        exact=exact_checkpoint_audit(base,F0,panels)
        n0=len(base)
        histories=[]
        for seed in seeds:
            st,F,P=deepcopy(base),F0,P0
            rng=Random(seed)
            counts=np.zeros(n0,dtype=np.int64)
            cumulative_p=cumulative_a=0.0
            cost_sum=f_sum=v_sum=tau=0
            saved_counts=[counts.copy()]
            saved_a=[0.0]; saved_p=[0.0]
            for t in range(1,horizon+1):
                F,P,e=native_step(st,F,P,rng,H)
                total,q=e['total'],e['q']
                p=q/total; p2=q*(q-1)/(total*(total-1))
                cumulative_p += p; cumulative_a += p-p2
                for z in e['served']:
                    if z<n0:
                        counts[z]+=1
                cost_sum+=e['cost']; f_sum+=e['forced_served']; v_sum+=e['relief_removed']
                tau+=len(e['served'])
                assert F==F0+cost_sum-f_sum-v_sum and P==P0+2*f_sum-v_sum
                if t%128==0:
                    saved_counts.append(counts.copy()); saved_a.append(cumulative_a); saved_p.append(cumulative_p)
            assert tuple(st.parents[:n0])==tuple(base.parents) and tau+f_sum==3*horizon
            histories.append(dict(seed=seed,counts=tuple(tuple(map(int,v)) for v in saved_counts),
                integrated_activity=tuple(saved_a),integrated_p=tuple(saved_p),
                final_n=len(st),final_F=F,final_P=P,maintenance=tau))
        count_array=np.asarray([h['counts'] for h in histories],dtype=np.int64)
        activity=np.array([h['integrated_activity'] for h in histories])
        integrated_p=np.array([h['integrated_p'] for h in histories])
        # A derived clock-only reference eliminates the global pool factor in the
        # ENSEMBLE identity: distance^2 = 2 E[D_path^2]/E[(N_3-N_4)^2].
        # The noisy finite-sample denominator is retained, never adjusted to its mean.
        reference_square=(count_array[:,-1,3]-count_array[:,-1,4]).astype(float)**2
        panel_results={}
        for length,panel in panels.items():
            path_counts=count_array @ panel['incidence'].T
            # Equal path length makes pair differences martingales: no fitted drift.
            drift=path_counts[:,-1,:]-length*integrated_p[:,-1,None]
            groups={}
            for delta,pairs in sorted(panel['groups'].items()):
                ii,jj=np.array(pairs,dtype=int).T
                differences=path_counts[:,:,ii]-path_counts[:,:,jj]
                per_history_m2=np.mean(differences[:,-1,:].astype(float)**2,axis=1)
                a=activity[:,-1]
                estimated_distance=float(per_history_m2.mean()/a.mean())
                # Delta-method Monte Carlo SE of ratio of means; pair dependence is
                # respected by averaging within EACH independent history first.
                residual=per_history_m2-estimated_distance*a
                ratio_se=float(residual.std(ddof=1)/(sqrt(len(a))*a.mean()))
                predicted_m2=delta*a
                if reference_square.mean()>0:
                    cr=float(2*per_history_m2.mean()/reference_square.mean())
                    cr_se=float((2*per_history_m2-cr*reference_square).std(ddof=1)
                                /(sqrt(len(a))*reference_square.mean()))
                    clock_only=dict(estimate=cr,monte_carlo_SE=cr_se,reference_pair=(3,4))
                else:
                    clock_only=dict(estimate=None,monte_carlo_SE=None,reference_pair=(3,4))
                phase_results={}
                for d in phase_depths:
                    modulus=3**d
                    phases=(differences%modulus)/modulus
                    z=np.exp(2j*np.pi*phases)
                    c1=z.mean(axis=2).mean(axis=0); c2=(z*z).mean(axis=2).mean(axis=0)
                    p0=(1+2*np.cos(2*np.pi*phases))**2/9
                    curve=p0.mean(axis=2).mean(axis=0)
                    assert np.max(np.abs(curve-(1/3+4*c1.real/9+2*c2.real/9)))<6e-14
                    phase_results[d]=dict(C1=tuple(map(complex,c1)),C2=tuple(map(complex,c2)),
                        mean_P0=tuple(map(float,curve)),final_P0=mean_se(p0[:,-1,:].mean(axis=1)))
                if delta==0:
                    assert np.max(np.abs(differences))==0
                # All pairs and seed blocks remain visible; these are NOT independent pairs.
                groups[delta]=dict(number_of_path_pairs=len(pairs),representative_pair=pairs[0],
                    predicted_squared_distance=delta,estimated_squared_distance=estimated_distance,
                    clock_only_distance=clock_only,
                    monte_carlo_SE=ratio_se,observed_difference_second_moment=mean_se(per_history_m2),
                    predicted_second_moment=mean_se(predicted_m2),
                    martingale_square_residual=mean_se(per_history_m2-predicted_m2),
                    seed_block_distance_estimates=tuple(float(per_history_m2[k:k+128].mean()/a[k:k+128].mean())
                                                        for k in (0,128)),phase_results=phase_results)
                # Independent terminal application through Cell24's tensor instrument.
                i,j=pairs[0]
                for k in (0,128):
                    u,v=map(int,path_counts[k,-1,[i,j]])
                    for d in phase_depths:
                        actual=readout((Q(u,3**d),-Q(v,3**d)))['probabilities']
                        target=np.array([(1+2*np.cos(2*np.pi*float((Q(u-v,3**d)-Q(s,3))%1)))**2/9
                                         for s in range(3)])
                        tensor_error=max(tensor_error,float(np.max(np.abs(np.array(actual)-target))))
            panel_results[length]=dict(paths=panel['paths'],supports=panel['supports'],
                distinct_supports=panel['unique_supports'],affine_dimension=panel['affine_dimension'],
                squared_distance=tuple(tuple(map(int,row)) for row in panel['squared_distance']),
                alias=panel['alias'],ultrametric_counterexample=panel['ultrametric_counterexample'],
                groups=groups,mean_path_martingale_residual=mean_se(drift.mean(axis=1)))
        # Passive observational wiring never participates in the native step.
        replay=deepcopy(base); rf,rp=F0,P0; rr=Random(seeds[0])
        for _ in range(horizon):
            rf,rp,_=native_step(replay,rf,rp,rr,H)
        assert (len(replay),rf,rp)==tuple(histories[0][k] for k in ('final_n','final_F','final_P'))
        results[name]=dict(H=H,checkpoint=checkpoint,parentage=tuple(base.parents),
            exact_service_audit=exact,panels=panel_results,histories=tuple(histories),
            mean_integrated_activity=mean_se(activity[:,-1]),
            reference_clock_pair=(3,4),reference_clock_second_moment=mean_se(reference_square),
            final_population=mean_se([h['final_n'] for h in histories]))
    assert tensor_error<8e-14 and (old,instrument)==before
    print('CELL 26 — FIXED RECORD PATHS DRIVE DIFFERENTIAL CLOCKS')
    print('Same 1/27 phase unit (1/81 control); same original two-sibling readout.')
    print('A served register drives every selected path containing it; shared contributions cancel.')
    print('This readout fan-out is ADDED; it creates no native work or local communication law.')
    print(f'Two past-only checkpoints at old step {freeze_step}; {len(seeds)} new continuations each x {horizon} updates.')
    print('Primary record length 3, length 4 transfer. Phase resolution is a DIFFERENT variable.')
    print('Exact: Cov(Y_A,Y_B)=(p-p2)|A intersect B|+(p2-p^2)|A||B|.')
    print('Equal-size supports: Var(Y_A-Y_B)=(p-p2)|A symmetric-difference B|.')
    print('Accumulation: E[(N_A-N_B)^2] = |A symmetric-difference B| E[sum(p-p2)].')
    print('Clock-only identity: distance^2 = 2 E[(N_A-N_B)^2]/E[(N_3-N_4)^2]; noisy estimates saved.')
    for name,_ in environments:
        r=results[name]; cp=r['checkpoint']
        print(f'\n{name}: past checkpoint n={cp["n"]}, F={cp["F"]}, P={cp["P"]}.')
        for length,panel in r['panels'].items():
            print(f'  record length {length}: {len(panel["paths"])} paths, {panel["distinct_supports"]} support classes; '
                  f'exact affine feature dimension {panel["affine_dimension"]}.')
            print('    delta    pairs    estimated delta +/- MC SE       mean P0, 1/27')
            for delta,g in panel['groups'].items():
                print(f'    {delta:3d} {g["number_of_path_pairs"]:8d} {g["estimated_squared_distance"]:12.6f}'
                      f' +/- {g["monte_carlo_SE"]:.6f}          {g["phase_results"][3]["final_P0"]["mean"]:.6f}')
            if panel['alias'] is not None:
                aa=panel['alias']
                print(f'  Clock alias pair {aa["pair"]}: {aa["paths"]}; address agreement probability='
                      f'{aa["address_agreement_probability"]}, despite identical clocks on EVERY history.')
            if panel['ultrametric_counterexample'] is not None:
                print('  Not an ultrametric: triple squared distances',panel['ultrametric_counterexample']['squared_distances'])
    print(f'\nPASS: complete legal recorded paths; exact weighted token-subset enumeration; {2*len(seeds)*horizon:,} native updates;')
    print(f'      all prior states preserved; tensor readout error {tensor_error:.3e}; same-history passive controls.')
    print('Empirical estimates are noisy and retained without selecting successes; pairs share histories.')
    print('Operational proximity = overlap of prescribed register supports, not physical distance or 3D.')
    print('Entry orientation, prefix order, and spatial propagation are NOT supplied by this support-only clock.')
    print('No proton/neutron/u-d identification, maturity epoch, SI scale, or formation-time prediction is made.')
    return dict(protocol='record-supported passive clocks; exact overlap kernel and native continuations',
        freeze_step=freeze_step,horizon=horizon,seeds=seeds,checkpoints=checkpoints,
        record_lengths=record_lengths,phase_depths=phase_depths,primary_record_length=3,primary_phase_depth=3,
        environments=results,total_new_native_updates=2*len(seeds)*horizon,
        max_tensor_error=tensor_error,readout_fanout_added=True,native_dynamics_changed=False,
        entry_reflection_used_as_phase_sign=False,physical_geometry_claimed=False,
        spatial_dimension_identified=None,particle_identifications={},physical_calibration=None,
        helpers=dict(difference_law=difference_law,exact_rank=exact_rank))


dcu_mass_26 = _run_dcu_mass_26()
