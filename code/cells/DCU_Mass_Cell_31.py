# CELL 31 — descendant-carrier readers versus fixed anchors. Paste after Cell30.
# NumPy + stdlib; no network. Each policy shadows the SAME unchanged native history.
# Following, private sibling selection, and memory retention are ADDED observer rules.
# A shared carrier ends a DISTINCT-observer protocol, never the native construction.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations, product
from math import comb, sqrt
from random import Random
import hashlib
import json


def _run_dcu_mass_31():
    import numpy as np
    missing = [k for k in ('_Native', 'dcu_mass_25', 'dcu_mass_30') if k not in globals()]
    if missing:
        raise RuntimeError('Run Cell30 first. Missing: ' + ', '.join(missing))
    old25, old30 = dcu_mass_25, dcu_mass_30
    if (old25['Gamma'], old25['m']) != (3, 0) or old30['horizon'] != 4096:
        raise ValueError('Inherited native protocol changed.')
    if tuple(old30['receivers']) != (3, 4):
        raise ValueError('Inherited receiver convention changed.')
    native_step = old25['helpers']['native_step']
    horizon = 4096
    checkpoints = (0, 32, 128, 512, 1024, 2048, 4096)
    policies = ('fixed', 'follow_any', 'follow_unique')
    cohorts = (('replayed', tuple(old30['seeds'])),
               ('fresh', tuple(20320920 + i for i in range(128))))
    environments = (('no_relief', 0), ('relief_H6', 6))
    arms, receivers = ('recorded_pair', 'first_use_pair'), (3, 4)
    Q = Fraction

    def get(mapping, key):
        return mapping[key] if key in mapping else mapping[str(key)]

    def stats(values):
        a = np.asarray(values, dtype=float)
        return dict(n=len(a), mean=float(a.mean()) if len(a) else None,
                    se=float(a.std(ddof=1)/sqrt(len(a))) if len(a) > 1 else None)

    def choose(n, k):
        return comb(n, k) if 0 <= k <= n else 0

    def own_child_law(n, F, q, degree):
        """Incident newborn count, conditional on the current carrier's service."""
        if not (0 <= degree < n and 1 <= q <= F+n):
            raise ValueError('Invalid finite service domain.')
        absent, neutral = n-1-degree, F+degree
        denom = choose(F+n-1, q-1)
        out = {j: Q(choose(absent,j)*choose(neutral,q-1-j),denom) for j in range(q)}
        assert sum(out.values()) == 1
        return out

    def construct(parentage):
        st = _Native()
        for expected, pair in enumerate(parentage[2:], 2):
            born, _ = st.add_batch([tuple(pair)])
            assert born == [expected]
        return st

    class Reader:
        def __init__(self, st, origin, initial_known, seed, role):
            self.origin = self.active = origin
            self.ticks = self.hops = self.ambiguous = 0
            self.expected_hops_at_ticks = 0.0
            self.known = set(initial_known)
            self.alias = {z: f'r{i}' for i,z in enumerate(sorted(self.known))}
            self.initial_cards = self.cards(st, self.known)
            self.initial_active = self.alias[origin]
            self.log, self.track = [], [origin]
            # Independent observer RNG; never used by the native transition.
            # Source streams agree between B3/B4 until a pair protocol stops.
            self.rng = Random(seed*1000003 + 1009*role + origin)
            self.shell_after_first_hop_checks = 0

        def cards(self, st, vertices):
            return tuple((self.alias[z], None if st.parents[z] is None else
                          tuple(self.alias[p] for p in st.parents[z])) for z in sorted(vertices))

        def service(self, st, e, own, degree, policy):
            before = self.active
            assert before in e['served'] and e['q']==3
            self.ticks += 1
            R = e['total']; neutral = e['F']+degree
            absent = e['n']-1-degree; den = (R-1)*(R-2)
            hp = (0.0 if policy == 'fixed' else
                  1.0-neutral*(neutral-1)/den if policy == 'follow_any' else
                  2.0*absent*neutral/den)
            self.expected_hops_at_ticks += hp
            self.ambiguous += len(own) > 1
            if policy == 'follow_any' and self.hops:
                assert degree == 0
                self.shell_after_first_hop_checks += 1
            gained = set()
            for z in own:
                gained.update(st.ancestors[z]-self.known)
            self.known.update(gained)
            for z in sorted(gained):
                self.alias[z] = f'r{len(self.alias)}'
            after = before
            if policy == 'follow_any' and own:
                after = own[0] if len(own) == 1 else own[self.rng.randrange(len(own))]
            elif policy == 'follow_unique' and len(own) == 1:
                after = own[0]
            if after != before:
                assert after in e['born'] and before in st.parents[after]
                assert after not in st.recorded
                self.hops += 1; self.track.append(after)
            self.active = after
            self.log.append(dict(local_tick=self.ticks,
                carrier_before=self.alias[before], carrier_after=self.alias[after],
                new_children=tuple(self.alias[z] for z in own),
                partners=tuple(self.alias[st.parents[z][1] if st.parents[z][0] == before
                                           else st.parents[z][0]] for z in own),
                newly_known_parent_cards=self.cards(st,gained)))
            assert st.ancestors[self.active] <= self.known

        def export(self):
            return dict(origin=self.origin, active=self.active, ticks=self.ticks,
                hops=self.hops, ambiguous_services=self.ambiguous,
                expected_hops_at_ticks=self.expected_hops_at_ticks,
                initial_active_handle=self.initial_active, initial_local_cards=self.initial_cards,
                local_log=tuple(self.log), auditor_handle_map=dict(self.alias),
                auditor_track=tuple(self.track), final_known=tuple(sorted(self.known)),
                shell_after_first_hop_checks=self.shell_after_first_hop_checks)

    class Protocol:
        def __init__(self, st, source, receiver, marker, policy, seed):
            self.policy,self.source,self.receiver,self.marker = policy,source,receiver,marker
            self.readers = (Reader(st,source,st.ancestors[source]|st.ancestors[marker],seed,1),
                            Reader(st,receiver,st.ancestors[receiver],seed,2))
            self.tags = [int(marker in aa) for aa in st.ancestors]
            self.roots, self.remote_roots = {}, set()
            self.receipt = self.first_remote = self.returned = self.merge = None
            self.stopped, self.stop_step = False, None
            self.clock_p = self.clock_a = 0.0
            self.saved = []; self.shared_receipt_births = 0
            self.save(0)

        def save(self,t):
            self.saved.append(dict(step=t,stopped=self.stopped,
                ticks=tuple(r.ticks for r in self.readers),hops=tuple(r.hops for r in self.readers),
                active=tuple(r.active for r in self.readers),receipt=self.receipt is not None,
                returned=self.returned is not None,merge=self.merge is not None))

        def observe(self,st,e,t,degree,integrated_p,integrated_a):
            if self.stopped:
                return
            assert len(self.tags) == e['n']
            a,b = (r.active for r in self.readers)
            assert a != b and 0 <= a < e['n'] and 0 <= b < e['n']
            W = set(e['served'])
            newborn = tuple(zip(e['born'],e['new_pairs']))
            assert all(u < e['n'] and v < e['n'] for _,(u,v) in newborn)
            inherited = {z:self.tags[u]|self.tags[v] for z,(u,v) in newborn}
            own_a = tuple(z for z,pair in newborn if a in pair)
            own_b = tuple(z for z,pair in newborn if b in pair)
            rec = tuple(z for z in own_b if inherited[z]&1)
            # Bit2 is ancestry of a PREVIOUS, nonjointly-witnessed receiver receipt.
            ret = tuple(z for z in own_a if inherited[z]&2)
            if a in W:
                self.readers[0].service(st,e,own_a,degree[a],self.policy)
            if b in W:
                self.readers[1].service(st,e,own_b,degree[b],self.policy)
            merged = self.readers[0].active == self.readers[1].active
            self.clock_p,self.clock_a = integrated_p,integrated_a

            def card(z):
                return dict(auditor_step=t,child=z,parents=tuple(st.parents[z]),
                    active_before=(a,b),active_after=tuple(r.active for r in self.readers),
                    source_tick=self.readers[0].ticks,receiver_tick=self.readers[1].ticks,
                    source_hops=self.readers[0].hops,receiver_hops=self.readers[1].hops,
                    complete_burst_pairs=tuple(e['new_pairs']),full_burst_work=e['cost'],
                    at_coalescence=merged)

            if self.returned is None and ret:
                z=ret[0]; candidates=sorted(self.remote_roots&st.ancestors[z])
                assert candidates
                r=candidates[0]; rc=self.roots[r]
                assert rc['auditor_step']<t and self.marker in st.ancestors[r]
                assert r in self.readers[0].known and a in st.parents[z]
                self.returned=dict(card(z),receipt_child=r,
                    receipt_active_carrier=rc['active_before'][1],receipt_birth_step=rc['auditor_step'],
                    identity_verified_from_auditor_track=True,
                    tracked_identity_authenticated_from_parentage_alone=(self.policy=='fixed'))
            for z in rec:
                is_shared=z in own_a
                rc=dict(card(z),jointly_witnessed_with_source=is_shared)
                self.roots[z]=rc
                self.shared_receipt_births+=int(is_shared)
                if self.receipt is None:self.receipt=rc
                if not is_shared:
                    self.remote_roots.add(z)
                    if self.first_remote is None:self.first_remote=rc
            for z,_ in newborn:
                self.tags.append(inherited[z]|(2 if z in self.remote_roots else 0))
            assert len(self.tags)==len(st)
            assert (self.marker in self.readers[1].known)==(self.receipt is not None)
            if merged:
                common=self.readers[0].active
                assert common in e['born'] and set(st.parents[common])=={a,b}
                self.merge=dict(card(common),shared_carrier=common)
                self.stopped,self.stop_step=True,t
            if t in checkpoints or merged:self.save(t)

        def finish(self,st):
            if self.stop_step is None:self.stop_step=horizon
            for z in range(len(self.tags)):
                expected=int(self.marker in st.ancestors[z])
                expected|=2*int(bool(self.remote_roots&st.ancestors[z]))
                assert self.tags[z]==expected
            for rd in self.readers:
                assert rd.hops<=rd.ticks
                assert all(u in st.parents[v] for u,v in zip(rd.track,rd.track[1:]))
            return dict(policy=self.policy,receiver=self.receiver,source=self.source,
                receipt=self.receipt,first_remote_receipt=self.first_remote,
                returned_remote_receipt=self.returned,coalescence=self.merge,stop_step=self.stop_step,
                integrated_p_until_stop=self.clock_p,clock_difference_variance_compensator=2*self.clock_a,
                readers=tuple(rd.export() for rd in self.readers),
                all_receipt_cards=tuple(self.roots[z] for z in sorted(self.roots)),
                jointly_witnessed_receipt_births=self.shared_receipt_births,checkpoints=tuple(self.saved))

    def run_history(cp,seed,expected=None,keep_trace=False):
        st=construct(cp['initial_parentage']); start_parentage=tuple(st.parents)
        F0,P0,H=int(cp['initial_F']),int(cp['initial_P']),int(cp['H'])
        source,marker=int(cp['source']),int(cp['marker'])
        F,P=F0,P0; rng=Random(seed); degree=[0]*len(st)
        for pair in st.parents[2:]:
            for p in pair:degree[p]+=1
        protocols={(policy,b):Protocol(st,source,b,marker,policy,seed)
                   for policy in policies for b in receivers}
        integrated_p=integrated_a=0.0; tau=sf=sv=sc=0; trace=[]
        for t in range(1,horizon+1):
            F,P,e=native_step(st,F,P,rng,H)
            p=e['q']/e['total']; p2=e['q']*(e['q']-1)/(e['total']*(e['total']-1))
            integrated_p+=p;integrated_a+=p-p2
            for protocol in protocols.values():
                protocol.observe(st,e,t,degree,integrated_p,integrated_a)
            for pair in e['new_pairs']:
                for parent in pair:degree[parent]+=1
            degree.extend([0]*len(e['born']))
            sf+=e['forced_served'];sv+=e['relief_removed'];sc+=e['cost'];tau+=len(e['served'])
            assert F==F0+sc-sf-sv and P==P0+2*sf-sv and tau+sf==3*t
            if keep_trace:trace.append(dict(e,step=t,F_after=F,P_after=P))
        assert tuple(st.parents[:len(start_parentage)])==start_parentage
        out={policy:{b:protocols[policy,b].finish(st) for b in receivers} for policy in policies}
        if expected is not None:
            assert (len(st),F,P,tau,sc)==tuple(expected[k] for k in
                ('final_n','final_F','final_P','maintenance','construction_work'))
            for b in receivers:
                row=out['fixed'][b]
                for newkey,oldkey in (('receipt','receipt'),('returned_remote_receipt','returned')):
                    actual,archived=row[newkey],get(expected[oldkey],b)
                    assert (actual is None)==(archived is None)
                    if actual:
                        assert actual['child']==archived['child'] and actual['auditor_step']==archived['auditor_step']
                        assert actual['source_tick']==archived['source_ticks']
                assert row['readers'][0]['ticks']==get(expected['clocks'],source)
                assert row['readers'][1]['ticks']==get(expected['clocks'],b)
                assert not row['jointly_witnessed_receipt_births'] and row['coalescence'] is None
        return dict(seed=seed,policies=out,final_n=len(st),final_F=F,final_P=P,maintenance=tau,
            construction_work=sc,final_parents=tuple(st.parents),
            native_parentage_sha256=hashlib.sha256(json.dumps(st.parents).encode()).hexdigest(),
            global_auditor_trace=tuple(trace),old_fixed_protocol_replayed=expected is not None)

    def summarize(histories):
        result={}
        for policy in policies:
            result[policy]={}
            for b in receivers:
                rows=[h['policies'][policy][b] for h in histories]
                receipts=[r for r in rows if r['receipt'] is not None]
                returned=[r for r in rows if r['returned_remote_receipt'] is not None]
                before=[r for r in returned if not r['returned_remote_receipt']['at_coalescence']]
                merged=[r for r in rows if r['coalescence'] is not None]
                D=[r['readers'][0]['ticks']-r['readers'][1]['ticks'] for r in rows]
                sum_ticks=sum(rd['ticks'] for r in rows for rd in r['readers'])
                sum_hops=sum(rd['hops'] for r in rows for rd in r['readers'])
                result[policy][b]=dict(histories=len(rows),
                    receipts_before_coalescence=sum(not r['receipt']['at_coalescence'] for r in receipts),
                    receipts_at_coalescence=sum(r['receipt']['at_coalescence'] for r in receipts),
                    all_receipts=len(receipts),
                    first_receipts_jointly_witnessed=sum(r['receipt']['jointly_witnessed_with_source'] for r in receipts),
                    remote_returns_before_coalescence=len(before),
                    remote_returns_at_coalescence=len(returned)-len(before),all_remote_returns=len(returned),
                    coalescences=len(merged),merges_before_any_receipt=sum(r['receipt'] is None for r in merged),
                    pending_without_receipt_at_horizon=sum(r['receipt'] is None and r['coalescence'] is None for r in rows),
                    pending_after_receipt_at_horizon=sum(r['receipt'] is not None and r['returned_remote_receipt'] is None
                                                       and r['coalescence'] is None for r in rows),
                    stopped_clock_difference=stats(D),
                    clock_square_minus_compensator=stats([d*d-r['clock_difference_variance_compensator'] for d,r in zip(D,rows)]),
                    source_ticks=stats([r['readers'][0]['ticks'] for r in rows]),
                    receiver_ticks=stats([r['readers'][1]['ticks'] for r in rows]),
                    hops_per_serviced_tick=sum_hops/sum_ticks if sum_ticks else None,total_carrier_hops=sum_hops,
                    ambiguous_services=sum(rd['ambiguous_services'] for r in rows for rd in r['readers']),
                    hop_prediction_residual=stats([sum(rd['hops']-rd['expected_hops_at_ticks'] for rd in r['readers']) for r in rows]),
                    mean_terminal_track_lengths=tuple(float(np.mean([len(r['readers'][i]['auditor_track']) for r in rows])) for i in (0,1)),
                    completed_only_source_ticks=stats([r['returned_remote_receipt']['source_tick'] for r in before]),
                    common_carrier_stops_are_not_independent_roundtrips=True)
        paired={}
        for policy in policies[1:]:
            paired[policy]={}
            for b in receivers:
                vals={}
                for name,key in (('receipt','receipt'),('remote_return','returned_remote_receipt')):
                    differences=[]
                    for h in histories:
                        r,fixed=h['policies'][policy][b],h['policies']['fixed'][b]
                        e,f=r[key],fixed[key]
                        differences.append(int(e is not None and not e['at_coalescence'])-int(f is not None))
                    vals[name+'_before_coalescence_fraction_change']=stats(differences)
                paired[policy][b]=vals
        return dict(by_policy=result,paired_against_fixed=paired)

    def continuation_identity_counterexample():
        # SAME native events and SAME source-local cards, two positive-probability
        # private sibling selections. Parentage cannot authenticate that choice.
        cp=old30['environments']['no_relief']['recorded_pair']
        st=construct(cp['initial_parentage']); F,P=cp['initial_F'],cp['initial_P']
        source,marker,b=cp['source'],cp['marker'],3
        candidates=[z for z in range(len(st)) if z not in (source,marker,b)
                    and tuple(sorted((b,z))) not in st.pair_to_id]
        x,y=candidates[:2]
        protocols=[Protocol(st,source,b,marker,'follow_any',0) for _ in range(2)]
        class Choice:
            def __init__(self,index):self.index=index;self.calls=0
            def randrange(self,n):
                assert n==2 and self.calls==0
                self.calls+=1
                return self.index
        protocols[0].readers[1].rng=Choice(0)
        protocols[1].readers[1].rng=Choice(1)
        class Draw:
            def __init__(self,selected,n):
                self.draw=list(selected)+([n] if len(selected)==2 else [])
            def sample(self,population,k):
                assert k==3 and len(set(self.draw))==3 and all(z in population for z in self.draw)
                return self.draw
        events=[]; integrated_p=integrated_a=0.0
        for t in (1,2,3):
            if t==1:selected=(b,x,y)
            elif t==2:selected=(branch_v,marker)
            else:selected=(source,receipt_z)
            degree=[0]*len(st)
            for pair in st.parents[2:]:
                for z in pair:degree[z]+=1
            F,P,e=native_step(st,F,P,Draw(selected,len(st)),0)
            probability=e['q']/e['total']
            pair_probability=e['q']*(e['q']-1)/(e['total']*(e['total']-1))
            integrated_p+=probability;integrated_a+=probability-pair_probability
            for protocol in protocols:protocol.observe(st,e,t,degree,integrated_p,integrated_a)
            events.append(dict(e,step=t,F_after=F,P_after=P))
            if t==1:
                incident=[z for z,pair in zip(e['born'],e['new_pairs']) if b in pair]
                assert len(incident)==2
                branch_u,branch_v=incident
                assert protocols[0].readers[1].active==branch_u
                assert protocols[1].readers[1].active==branch_v
            elif t==2:
                receipt_z=e['born'][0]
                assert protocols[0].receipt is None and protocols[1].receipt['child']==receipt_z
        assert protocols[0].readers[0].log==protocols[1].readers[0].log
        assert protocols[0].returned is None and protocols[1].returned is not None
        assert protocols[0].merge is None and protocols[1].merge is not None
        return dict(initial_parentage=cp['initial_parentage'],source=source,receiver=b,marker=marker,
            branching_parents=(b,x,y),alternative_children=(branch_u,branch_v),
            receipt_child=receipt_z,return_child=protocols[1].returned['child'],
            source_local_log_identical=True,native_events_identical=True,
            tracked_receipt_exists=(False,True),tracked_return_exists=(False,True),
            final_carriers=tuple(tuple(r.active for r in p.readers) for p in protocols),
            original_rule_probability_each_choice=Q(1,2),
            source_local_log=tuple(protocols[0].readers[0].log),
            initial_source_cards=protocols[0].readers[0].initial_cards,
            source_handles=protocols[0].readers[0].alias,
            conditioned_native_events=tuple(events),final_parentage=tuple(st.parents),
            interpretation='A returned descendant-origin certificate does not authenticate a private continuation track.')

    def exact_audits():
        tests=[]
        for n in range(3,8):
            for F in (0,1,4):
                q=min(3,n+F)
                for degree in range(n):
                    counter=Counter()
                    for other in combinations(range(1,n+F),q-1):
                        counter[sum(degree<z<n for z in other)]+=1
                    law={j:Q(counter[j],choose(n+F-1,q-1)) for j in range(q)}
                    assert law==own_child_law(n,F,q,degree)
                    tests.append((n,F,q,degree,law))
        merge_cases=sum(a==b for a,b in product(('AB','AC'),('AB','BC')))
        assert Q(merge_cases,4)==Q(1,4)
        merger_tests=[]
        for n in range(3,9):
            for F in (0,1,4):
                R=n+F; denom=choose(R,3); weight_any=Q(0); weight_unique=Q(0)
                for W in combinations(range(R),3):
                    if 0 not in W or 1 not in W:continue
                    third=next(z for z in W if z not in (0,1))
                    weight_any+=Q(1) if third>=n else Q(1,4)
                    weight_unique+=int(third>=n)
                p2=Q(6,R*(R-1))
                expected_any=p2*(F+Q(n-2,4))/(R-2)
                expected_unique=p2*Q(F,R-2)
                assert weight_any/denom==expected_any and weight_unique/denom==expected_unique
                merger_tests.append((n,F,expected_any,expected_unique))
        return dict(conditional_child_laws=tests,conditional_child_law_cases=len(tests),
                    two_childless_carrier_merger_checks=merger_tests,
                    merger_two_parent_burst=Q(1),merger_three_parent_complete_burst=Q(1,4),
                    identity_counterexample=continuation_identity_counterexample())

    result={};total_updates=replayed_histories=0;exact=exact_audits()
    print('CELL 31 — FOLLOWING CARRIERS; NOT SPATIAL MOTION OR PARTICLE IDENTITY',flush=True)
    print('Primary follow_any: uniform incident-newborn choice with independent observer RNG.')
    print('Control follow_unique: advance only at a unique incident child. Fixed anchor unchanged.')
    print('One pre-burst carrier/tick; all native births and charges retained.')
    print('Stop distinct-reader protocols at first shared carrier; never stop native growth.')
    print('Moving identity uses auditor tracks, not authenticated native parentage alone.')
    for cohort,seeds in cohorts:
        result[cohort]={}
        for name,H in environments:
            result[cohort][name]={}
            for arm in arms:
                cp=old30['environments'][name][arm]
                snapshot=deepcopy((cp['initial_parentage'],cp['initial_F'],cp['initial_P']))
                assert H==cp['H']
                previous_by_seed={h['seed']:h for h in cp['histories']} if cohort=='replayed' else {}
                histories=tuple(run_history(cp,seed,previous_by_seed.get(seed),keep_trace=(i==0))
                                for i,seed in enumerate(seeds))
                total_updates+=len(seeds)*horizon
                replayed_histories+=sum(h['old_fixed_protocol_replayed'] for h in histories)
                assert (cp['initial_parentage'],cp['initial_F'],cp['initial_P'])==snapshot
                summary=summarize(histories)
                result[cohort][name][arm]=dict(source=cp['source'],marker=cp['marker'],H=H,
                    initial_parentage=cp['initial_parentage'],initial_F=cp['initial_F'],initial_P=cp['initial_P'],
                    histories=histories,summary=summary)
                print(f'\n{cohort}: {name}/{arm}, {len(seeds)} histories x {horizon} updates',flush=True)
                print(' policy          B   receipts(pre/merge)   returns(pre/merge)   merged  hop/tick',flush=True)
                for policy in policies:
                    for b in receivers:
                        s=summary['by_policy'][policy][b]
                        print(f' {policy:14s} {b}      {s["receipts_before_coalescence"]:3d}/{s["receipts_at_coalescence"]:<3d}'
                              f'              {s["remote_returns_before_coalescence"]:3d}/{s["remote_returns_at_coalescence"]:<3d}'
                              f'        {s["coalescences"]:3d}     {s["hops_per_serviced_tick"]:.5f}',flush=True)
    print(f'\nPASS: {replayed_histories} archived histories, both fixed receivers replayed exactly.')
    print(f'Native updates evaluated: {total_updates:,}; observer policies share each history.')
    print(f'Exact conditional-child laws: {exact["conditional_child_law_cases"]}.')
    print('Adaptive clocks: E[dN_A-dN_B|past]=0 while distinct; variance=2(p-p2).')
    print('PASS: same source-local log and native events can hide different private continuation tracks.')
    print('Continuation choices are observer rules, NOT autonomous personal identity.')
    print('No source mass formula, neutron worker, service pool, spatial metric or SI scale changed.')
    return dict(protocol='v1 descendant-carrier versus fixed reader on shared native histories',
        horizon=horizon,checkpoints=checkpoints,policies=policies,primary_policy='follow_any',
        receivers=receivers,primary_receiver=3,cohort_seeds={name:seeds for name,seeds in cohorts},
        cohorts=result,exact=exact,total_native_updates_evaluated=total_updates,
        archived_native_updates_replayed=len(cohorts[0][1])*4*horizon,
        fresh_native_updates=len(cohorts[1][1])*4*horizon,archived_histories_replayed=replayed_histories,
        observer_choice_rule_added=True,memory_retention_added=True,native_transition_changed=False,
        merger_is_competing_endpoint=True,moving_identity_authenticated_from_native_parentage=False,
        physical_worldline_or_speed_claimed=False,proper_time_claimed=False,mass_dictionary_changed=False,
        neutron_search_modified=False,helpers=dict(construct=construct,run_history=run_history,
            own_child_law=own_child_law,Reader=Reader,Protocol=Protocol,summarize=summarize))


dcu_mass_31 = _run_dcu_mass_31()
