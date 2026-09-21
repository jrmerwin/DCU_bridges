# CELL 29 — a legal comparison -> native work -> clock/readout back-action.
# Paste after Cell 28; execution needs only Cells 25/26 (and their prerequisites).
# NumPy + standard library. No downloads, external constants, particle labels or fits.
# THREE CONDITIONED LEGAL SERVICE EVENTS, not a new autonomous selection rule:
# existing pair; first canonical absent recorded/recorded pair; first canonical
# absent recorded/unrecorded pair. Continue with the UNCHANGED native sampler.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations
from math import comb, sqrt
from random import Random


def _run_dcu_mass_29():
    import numpy as np
    required = ('_Native', 'dcu_mass_24', 'dcu_mass_25', 'dcu_mass_26')
    missing = [key for key in required if key not in globals()]
    if missing:
        raise RuntimeError('Run Cells 25/26 first. Missing: ' + ', '.join(missing))
    old25, old26, instrument = dcu_mass_25, dcu_mass_26, dcu_mass_24
    before = deepcopy((old25, old26, instrument))
    if (old25['Gamma'], old25['m'], old25['carriers'], old25['primary_depth']) != (3,0,(3,4),3):
        raise ValueError('The inherited service/clock convention changed.')
    if old26['freeze_step'] != 512 or instrument['local_phase_labels'] != (0,1,2):
        raise ValueError('The source checkpoint or phase generator changed.')
    step = old25['helpers']['native_step']
    readout = instrument['helpers']['sibling_readout']
    Q = Fraction
    horizon = 1024
    checkpoints = (0,1,8,32,128,512,1024)
    seeds = tuple(20290920+i for i in range(512))
    depths = (3,4)
    clock = (3,4)
    reference = (0,1)
    environments = (('no_relief',0), ('relief_H6',6))
    arm_order = ('existing_pair', 'recorded_pair', 'first_use_pair')

    def stats(values):
        a = np.asarray(values, dtype=float)
        return dict(mean=float(a.mean()), se=float(a.std(ddof=1)/sqrt(len(a))))

    def choose(n,k):
        return comb(n,k) if 0 <= k <= n else 0

    class ConditionedDraw:
        def __init__(self, pair, n):
            self.draw = tuple(pair)+(n,)  # one actual outstanding-work token
            self.calls = 0
        def sample(self, population, k):
            assert k == 3 and self.calls == 0
            assert len(set(self.draw)) == 3 and all(v in population for v in self.draw)
            self.calls += 1
            return list(self.draw)

    def construct(parentage):
        st = _Native()
        for index, pair in enumerate(parentage[2:],2):
            born, _ = st.add_batch([tuple(pair)])
            assert born == [index]
        assert tuple(st.parents) == tuple(parentage)
        return st

    def select(st):
        # Cache labels are NOT physical types. Order by complete parentage strings.
        terms = ['a','b']
        for a,b in st.parents[2:]:
            terms.append('('+'|'.join(sorted((terms[a],terms[b])))+')')
        pool = [z for z in range(2,len(st)) if z not in clock]
        pairs = sorted((e for e in combinations(pool,2) if e not in st.pair_to_id),
                       key=lambda e:tuple(sorted(terms[z] for z in e)))
        both = [e for e in pairs if all(z in st.recorded for z in e)]
        first = [e for e in pairs if sum(z in st.recorded for z in e)==1]
        if not both or not first:
            raise ValueError('A declared intervention class is empty; no later checkpoint substituted.')
        return dict(existing_pair=(0,1), recorded_pair=both[0], first_use_pair=first[0]), \
               dict(recorded_pair_pool=len(both), first_use_pair_pool=len(first),
                    selection='first complete-parentage lexicographic pair in each status class',
                    selection_uses_mass_or_clock_effect=False)

    def one_step(st,F):
        n,R = len(st),len(st)+F
        p,p2=Q(3,R),Q(6,R*(R-1))
        law={0:1-2*p+p2,1:2*(p-p2),2:p2}
        counted=Counter(); subsets=0
        # Sum over ALL maintenance subsets, weighting their forced-token completions.
        for s in range(4):
            multiplicity=choose(F,3-s)
            if not multiplicity: continue
            for W in combinations(range(n),s):
                subsets += 1
                counted[sum(z in W for z in clock)] += multiplicity
        assert {k:Q(v,comb(R,3)) for k,v in counted.items()} == law
        assert sum(law.values()) == 1 and sum(k*v for k,v in law.items())==2*p
        C1=sum(float(v)*np.exp(2j*np.pi*k/27) for k,v in law.items())
        C2=sum(float(v)*np.exp(4j*np.pi*k/27) for k,v in law.items())
        P0=sum(float(v)*(1+2*np.cos(2*np.pi*k/27))**2/9 for k,v in law.items())
        assert abs(P0-(1/3+4*C1.real/9+2*C2.real/9))<2e-14
        return dict(n=n,F=F,pool=R,p=p,p2=p2,tagged_increment_law=law,
                    tagged_mean=2*p,mean_global_ticks=n*p,
                    expected_P0_after_one_step=float(P0),
                    weighted_subsets_checked=subsets)

    def trajectory(state,F0,P0,H,seed):
        st,F,P=deepcopy(state),F0,P0
        rng=Random(seed)
        counts=[0,0,0,0]  # 3,4,0,1: clock pair, then reference pair
        tau=sf=sv=sc=0
        integrated_p=chi=0.0
        saved=[]
        def save(t):
            saved.append((t,len(st),F,P,tau,*counts,integrated_p,chi,sc,sf,sv))
        save(0)
        for t in range(1,horizon+1):
            F,P,e=step(st,F,P,rng,H)
            for k,z in enumerate(clock+reference):counts[k]+=int(z in e['served'])
            s=len(e['served']);tau+=s
            integrated_p+=e['q']/e['total'];chi+=s/e['n']
            sf+=e['forced_served'];sv+=e['relief_removed'];sc+=e['cost']
            assert F==F0+sc-sf-sv and P==P0+2*sf-sv and tau+sf==3*t
            if t in checkpoints:save(t)
        assert tuple(st.parents[:len(state)])==tuple(state.parents)
        return dict(seed=seed,rows=tuple(saved))

    output={}; tensor_error=0.0
    for name,H in environments:
        previous=old26['environments'][name]
        base=construct(previous['parentage'])
        cp=previous['checkpoint'];F0,P0=cp['F'],cp['P'];n0=len(base)
        assert (n0,F0,P0)==(cp['n'],cp['F'],cp['P']) and F0>0
        pairs,selection=select(base)
        arms={}
        for arm in arm_order:
            pair=pairs[arm];st=deepcopy(base)
            old_anc=tuple(frozenset(a) for a in st.ancestors)
            old_recorded=frozenset(st.recorded)
            conditional=ConditionedDraw(pair,n0)
            F,P,event=step(st,F0,P0,conditional,H)
            assert conditional.calls==1 and event['served']==tuple(sorted(pair))
            assert event['forced_served']==1
            assert tuple(frozenset(st.ancestors[z]) for z in range(n0))==old_anc
            assert not (set(pair)&set(clock))
            new=arm!='existing_pair'
            assert len(event['born'])==int(new)
            if new:
                assert event['new_pairs']==(pair,)
                used=base.ancestors[pair[0]] | base.ancestors[pair[1]]
                repeat=2*sum(base.path_weight(z) for z in used)
                premium=9*sum(base.path_weight(z) for z in used if z not in old_recorded)
                assert event['cost']==repeat+premium
                assert (premium==0)==(arm=='recorded_pair')
            else:
                repeat=premium=0
                assert event['cost']==0
            exact=one_step(st,F)
            histories=tuple(trajectory(st,F,P,H,seed) for seed in seeds)
            a=np.asarray([h['rows'] for h in histories],dtype=float)
            K=a[:,:,5]+a[:,:,6];B=a[:,:,7]+a[:,:,8]
            curve={}
            for d in depths:
                x=(K.astype(np.int64)%(3**d))/(3**d)
                z=np.exp(2j*np.pi*x)
                p0=(1+2*np.cos(2*np.pi*x))**2/9
                c1=z.mean(axis=0);c2=(z*z).mean(axis=0)
                assert np.max(np.abs(p0.mean(axis=0)-(1/3+4*c1.real/9+2*c2.real/9)))<5e-14
                for i in (0,len(seeds)-1):
                    for ti in (1,len(checkpoints)-1):
                        aa,bb=int(a[i,ti,5]),int(a[i,ti,6])
                        pred=readout((Q(aa,3**d),Q(bb,3**d)))['probabilities'][0]
                        tensor_error=max(tensor_error,abs(pred-p0[i,ti]))
                curve[d]=dict(C1=tuple(map(complex,c1)),C2=tuple(map(complex,c2)),
                    P0=tuple(stats(p0[:,i]) for i in range(len(checkpoints))))
            summary={t:dict(tagged_services=stats(K[:,i]),reference_services=stats(B[:,i]),
                tagged_minus_reference=stats(K[:,i]-B[:,i]),
                predicted_tagged=stats(2*a[:,i,9]),
                tagged_minus_compensator=stats(K[:,i]-2*a[:,i,9]),
                tagged_minus_activity_clock=stats(K[:,i]-2*a[:,i,10]),
                population=stats(a[:,i,1]),forced_work=stats(a[:,i,2]),maintenance=stats(a[:,i,4]),
                seed_block_tagged=(stats(K[:256,i]),stats(K[256:,i])))
                for i,t in enumerate(checkpoints)}
            arms[arm]=dict(pair=pair,event=event,post_parentage=tuple(st.parents),post_F=F,post_P=P,
                cost_repeat=repeat,cost_first_use_premium=premium,exact_next_step=exact,
                conditional_service_event_probability=Q(F0,comb(n0+F0,3)),
                history=histories,summary=summary,curves=curve)
        R0=arms['existing_pair']['exact_next_step']['pool']
        comparisons={}
        baseline=np.asarray([h['rows'] for h in arms['existing_pair']['history']],dtype=float)
        for arm in arm_order[1:]:
            r=arms[arm];C=r['event']['cost'];R=r['exact_next_step']['pool']
            assert R==R0+C+1
            ratio=r['exact_next_step']['p']/arms['existing_pair']['exact_next_step']['p']
            assert ratio==Q(R0,R0+C+1)
            # Same exact change for EVERY old token: an immediate common-mode effect.
            shifts=tuple(r['exact_next_step']['p']-arms['existing_pair']['exact_next_step']['p']
                         for _ in range(n0))
            assert len(set(shifts))==1
            a=np.asarray([h['rows'] for h in r['history']],dtype=float)
            changes={}
            for i,t in enumerate(checkpoints):
                K=a[:,i,5]+a[:,i,6];K0=baseline[:,i,5]+baseline[:,i,6]
                p0=(1+2*np.cos(2*np.pi*(K%27)/27))**2/9
                p00=(1+2*np.cos(2*np.pi*(K0%27)/27))**2/9
                changes[t]=dict(delta_tagged_services=stats(K-K0),delta_P0=stats(p0-p00))
            comparisons[arm]=dict(initial_rate_ratio=ratio,initial_fractional_rate_change=ratio-1,
                old_token_rate_shifts=shifts,paired_changes=changes)
        output[name]=dict(H=H,checkpoint=cp,selection=selection,arms=arms,comparisons=comparisons)
    assert tensor_error<7e-14 and (old25,old26,instrument)==before
    print('CELL 29 — CONDITIONAL RECORDING BACK-ACTION, NOT A GRAVITATIONAL FIELD')
    print('Same two saved checkpoints; no search process, particle label, mass dictionary or phase rule changed.')
    print('A legal 3-slot service selects two named old parents and one forced token; all enabled pairs form.')
    print('The event is conditioned, NOT freely selected by an autonomous internal observer.')
    print('New observation starts AFTER that event; prior graph, records, work and maintenance are retained.')
    print('Clock 3,4 and reference 0,1; phase increments 1/27 (primary), 1/81 (control).')
    print(f'{len(seeds)} continuations/arm x {horizon} updates; {2*3*len(seeds)*horizon:,} native updates.')
    for name,_ in environments:
        env=output[name]
        print(f'\n{name}: pre-event n,F,P = {env["checkpoint"]["n"]}, {env["checkpoint"]["F"]}, {env["checkpoint"]["P"]}')
        print(' arm              pair          cost(repeat+premium)   post n/F     next-token rate / control')
        for arm in arm_order:
            r=env['arms'][arm];ratio=Q(1) if arm=='existing_pair' else env['comparisons'][arm]['initial_rate_ratio']
            print(f' {arm:17s} {str(r["pair"]):13s} {r["cost_repeat"]:5d}+{r["cost_first_use_premium"]:<5d}'
                  f' {len(r["post_parentage"]):4d}/{r["post_F"]:<7d} {str(ratio):>12s} = {float(ratio):.6f}')
        for t in (128,1024):
            print(f' At {t} continuation updates: mean tagged services; paired change vs existing-pair control; P0')
            for arm in arm_order:
                r=env['arms'][arm];s=r['summary'][t]['tagged_services'];p0=r['curves'][3]['P0'][checkpoints.index(t)]
                change=0 if arm=='existing_pair' else env['comparisons'][arm]['paired_changes'][t]['delta_tagged_services']['mean']
                print(f'  {arm:17s}: {s["mean"]:.6f} +/- {s["se"]:.6f}; change={change:+.6f}; '
                      f'P0={p0["mean"]:.6f} +/- {p0["se"]:.6f}')
        print(' Terminal tagged-minus-reference and tagged-minus-2chi (Monte Carlo SE):')
        for arm in arm_order:
            s=env['arms'][arm]['summary'][horizon]
            print(f'  {arm:17s}: '+ '; '.join(f'{s[k]["mean"]:+.6f} +/- {s[k]["se"]:.6f}'
                  for k in ('tagged_minus_reference','tagged_minus_activity_clock')))
    print('\nExact immediate law: p_new/p_control = R0/(R0+C+1); same for all old tokens.')
    print('Longer-run differences include native growth/queue feedback; no monotonicity claim.')
    print('Outcome probabilities are nonlinear: slower phase accumulation need not increase P0.')
    print('Both local clock pairs have equal conditional drift; no relative redshift or locality is obtained.')
    print('PASS: actual original native_step, all enabled births, exact one-step enumeration, full budget checks,')
    print(f'      tensor readout (max error {tensor_error:.3e}), native ancestry preservation, prior dictionaries unchanged.')
    return dict(protocol='v1: conditioned legal recording event drives later unchanged native phase clocks',
        horizon=horizon,seeds=seeds,checkpoints=checkpoints,phase_depths=depths,clock_pair=clock,
        reference_pair=reference,arm_order=arm_order,environments=output,
        history_columns=('step','n','F','P','tau','N3','N4','N0','N1','integrated_p','chi',
                         'construction_cost','forced_served','relief_removed'),
        total_native_updates=2*3*len(seeds)*horizon,max_tensor_error=tensor_error,
        conditioned_initial_service=True,phase_rule_unchanged=True,
        native_transition_changed=False,new_gate_cost_assumed=False,
        physical_time_calibration=None,physical_energy_calibration=None,
        gravitational_redshift_identified=False,local_signal_propagation_identified=False,
        neutron_search_modified=False,particle_identification=None)


dcu_mass_29 = _run_dcu_mass_29()
