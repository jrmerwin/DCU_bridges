# CELL 30 — event-local readers: native pulse ancestry, receipt, and return.
# Paste after Cell29. Original _Native and Cell25 native_step; NumPy + stdlib.
# ADDED access convention: a reader sees the structural ancestry of children
# made using its own anchor. No global-state access, new pools, routing law,
# chosen-bit encoding, physical position, or costed memory implementation.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations
from math import comb, sqrt
from random import Random


def _run_dcu_mass_30():
    import numpy as np
    missing = [k for k in ('_Native', 'dcu_mass_25', 'dcu_mass_29') if k not in globals()]
    if missing:
        raise RuntimeError('Run Cells25 and 29 first. Missing: ' + ', '.join(missing))
    old25, old29 = dcu_mass_25, dcu_mass_29
    if (old25['Gamma'], old25['m'], tuple(old25['carriers']), old25['primary_depth']) != (3,0,(3,4),3):
        raise ValueError('Inherited counting/phase convention changed.')
    if tuple(old29['clock_pair']) != (3,4):
        raise ValueError('Cell29 clock-pair convention changed.')
    native_step = old25['helpers']['native_step']
    previous = deepcopy(old29)
    Q = Fraction
    horizon = 4096
    seeds = tuple(20300920+i for i in range(256))  # Seed labels, not calendar dates.
    checkpoints = (0,32,128,512,1024,2048,4096)
    receivers = (3,4)  # Primary and fixed transfer reader, on the SAME histories.
    environments = (('no_relief',0),('relief_H6',6))
    arms = ('recorded_pair','first_use_pair')

    def choose(n,k):
        return comb(n,k) if 0 <= k <= n else 0

    def stats(values):
        a = np.asarray(values,dtype=float)
        if len(a) == 0:
            return dict(n=0,mean=None,se=None)
        return dict(n=len(a),mean=float(a.mean()),
                    se=float(a.std(ddof=1)/sqrt(len(a))) if len(a)>1 else None)

    def conditional_contact(total,q,carriers):
        # Conditional on the receiver token being among q selected requests.
        if not 0 <= carriers < total or not 1 <= q <= total:
            raise ValueError('Bad finite contact domain.')
        return Q(1)-Q(choose(total-1-carriers,q-1),choose(total-1,q-1))

    def construct(parentage):
        st = _Native()
        for expected,pair in enumerate(parentage[2:],2):
            born,_ = st.add_batch([tuple(pair)])
            assert born == [expected]
        return st

    def graph_distances(st,start):
        adjacency=[set() for _ in st.parents]
        for child,pair in enumerate(st.parents[2:],2):
            for parent in pair:
                adjacency[child].add(parent);adjacency[parent].add(child)
        distance={start:0};queue=[start]
        for v in queue:
            for w in sorted(adjacency[v]):
                if w not in distance:
                    distance[w]=distance[v]+1;queue.append(w)
        return distance

    def exact_initial(st,F,marker,source):
        R,q = len(st)+F,3
        p,p2 = Q(q,R),Q(q*(q-1),R*(R-1))
        h=conditional_contact(R,q,1)
        assert p*h == p2 and h == Q(2,R-1)
        counters=Counter();sample_weight=0;event_numer=cost_numer=0;subsets=0
        target=receivers[0]
        for s in range(min(q,len(st))+1):
            weight=choose(F,q-s)
            if not weight:continue
            for W in combinations(range(len(st)),s):
                subsets+=1;sample_weight+=weight
                served=set(W)
                counters[int(target in served),int(marker in served)] += weight
                if target in served and marker in served:
                    events=[pair for pair in combinations(W,2) if pair not in st.pair_to_id]
                    assert tuple(sorted((target,marker))) in events
                    event_numer+=weight
                    cost_numer+=weight*st.recording_cost(events)
        denominator=choose(R,q)
        assert sample_weight==denominator and Q(event_numer,denominator)==p2
        assert Q(sum(v for (b,m),v in counters.items() if b),denominator)==p
        distances=graph_distances(st,marker)
        strata=Counter(distances[b] for b in range(marker) if b not in st.parents[marker])
        return dict(n=len(st),F=F,total=R,source=source,marker=marker,p=p,p2=p2,
            arrival_per_receiver_tick=h,weighted_subsets=subsets,
            joint_service_law={str(k):Q(v,denominator) for k,v in counters.items()},
            mean_full_burst_cost_given_first_step_receipt=Q(cost_numer,event_numer),
            available_receiver_parent_child_distance_counts=dict(strata),
            primary_distance=distances[target],same_hazard_at_all_those_distances=True)

    def run(base,F0,P0,H,source,marker,seed,keep_trace=False):
        st,F,P=deepcopy(base),F0,P0
        rng=Random(seed)
        anchors=(source,)+receivers
        clocks={o:0 for o in anchors}
        knowledge={o:set(base.ancestors[o]) for o in anchors}
        knowledge[source].update(base.ancestors[marker])  # Source participated in emission.
        initial_frames={o:tuple(sorted(v)) for o,v in knowledge.items()}
        # Handles are assigned inside EACH frame: raw global IDs are auditor-only.
        alias={o:{z:f'r{i}' for i,z in enumerate(sorted(knowledge[o]))} for o in anchors}
        initial_cards={o:tuple((alias[o][z],None if base.parents[z] is None else
            tuple(alias[o][p] for p in base.parents[z])) for z in sorted(knowledge[o])) for o in anchors}
        assert all(marker not in knowledge[b] for b in receivers)
        # AUDITOR tag bits encode pulse ancestry and existence of ANY native receipt;
        # they never drive native_step and are not mutable particle states.
        tags=[0]*len(st);tags[marker]=1
        carrier_counts=[1,0,0]
        receipt_roots={b:set() for b in receivers};receipt_birth={b:{} for b in receivers}
        receipt={b:None for b in receivers};returned={b:None for b in receivers}
        a_receipt={b:0.0 for b in receivers};a_return={b:0.0 for b in receivers}
        v_receipt={b:0.0 for b in receivers};v_return={b:0.0 for b in receivers}
        local_logs={o:[] for o in anchors};global_trace=[];saved=[]
        sf=sv=sc=tau=0
        def save(t):
            saved.append(dict(step=t,n=len(st),F=F,P=P,tau=tau,
                source_ticks=clocks[source],receiver_ticks={b:clocks[b] for b in receivers},
                receipt={b:receipt[b] is not None for b in receivers},
                returned={b:returned[b] is not None for b in receivers},
                source_frame_size=len(knowledge[source]),
                receiver_frame_size={b:len(knowledge[b]) for b in receivers},
                pulse_carriers=carrier_counts[0],receipt_compensator=dict(a_receipt),
                return_compensator=dict(a_return)))
        save(0)
        for t in range(1,horizon+1):
            n,R=len(st),len(st)+F
            q=min(3,R);p=q/R
            for k,b in enumerate(receivers,1):
                if receipt[b] is None:
                    # Before acquisition, no receiver/carrier pair has formed.
                    rate=p*float(conditional_contact(R,q,carrier_counts[0]))
                    a_receipt[b]+=rate;v_receipt[b]+=rate*(1-rate)
                elif returned[b] is None:
                    rate=p*float(conditional_contact(R,q,carrier_counts[k]))
                    a_return[b]+=rate;v_return[b]+=rate*(1-rate)
            F,P,e=native_step(st,F,P,rng,H)
            W=set(e['served'])
            for o in anchors:clocks[o]+=int(o in W)
            born=list(zip(e['born'],e['new_pairs']))
            # Pre-burst parents only: no same-burst sibling transmits a new receipt.
            assert all(u<n and v<n for _,(u,v) in born)
            inherited={z:tags[u]|tags[v] for z,(u,v) in born}
            new_receipts={};new_returns={};receipt_bits={}
            for k,b in enumerate(receivers,1):
                # EVERY actual receipt can supply a return certificate. A source
                # cannot know which receipt was first at the distant reader.
                options=[(z,pair) for z,pair in born if b in pair and inherited[z]&1]
                for z,pair in options:
                    receipt_roots[b].add(z);receipt_birth[b][z]=t
                    receipt_bits[z]=receipt_bits.get(z,0)|(1<<k)
                if receipt[b] is None and options:
                    z,pair=options[0]
                    new_receipts[b]=(z,pair,k)
                elif receipt[b] is not None and returned[b] is None:
                    options=[(z,pair) for z,pair in born if source in pair and inherited[z]&(1<<k)]
                    if options:new_returns[b]=options[0]
            for z,_ in born:
                value=inherited[z]
                value|=receipt_bits.get(z,0)
                tags.append(value)
                for k in range(3):carrier_counts[k]+=bool(value&(1<<k))
            assert len(tags)==len(st)
            # The ONLY new access consists of this reader's own incident newborns.
            for o in anchors:
                if o not in W:continue
                own=[(z,pair) for z,pair in born if o in pair]
                gained=set()
                for z,_ in own:gained.update(st.ancestors[z]-knowledge[o])
                knowledge[o].update(gained)
                for z in sorted(gained):alias[o][z]=f'r{len(alias[o])}'
                local_logs[o].append(dict(local_tick=clocks[o],
                    new_children=tuple(alias[o][z] for z,_ in own),
                    partners=tuple(alias[o][v if u==o else u] for _,(u,v) in own),
                    newly_known_parent_cards=tuple((alias[o][z],None if st.parents[z] is None else
                        tuple(alias[o][p] for p in st.parents[z])) for z in sorted(gained))))
            def event_card(z,pair):
                return dict(auditor_step=t,child=z,parents=pair,source_ticks=clocks[source],
                    receiver_ticks={b:clocks[b] for b in receivers},
                    complete_burst_pairs=tuple(e['new_pairs']),complete_burst_cost=e['cost'],
                    queued_work_after=F,ancestry=tuple(sorted(st.ancestors[z])))
            for b,(z,pair,k) in new_receipts.items():
                assert marker in st.ancestors[z] and b in pair and marker in knowledge[b]
                assert marker in st.recorded
                receipt[b]=event_card(z,pair)
                receipt[b]['direct_marker_contact']=marker in pair
                receipt[b]['source_record_in_received_ancestry']=source in st.ancestors[z]
                peer=pair[1] if pair[0]==b else pair[0]
                assert peer in st.recorded and source in st.ancestors[peer]
            for b,(z,pair) in new_returns.items():
                proofs=sorted(receipt_roots[b]&st.ancestors[z])
                assert proofs
                r=proofs[0]
                assert receipt_birth[b][r]<t and b in st.parents[r] and marker in st.ancestors[r]
                assert source in pair and r in knowledge[source] and r in st.recorded
                returned[b]=event_card(z,pair)
                returned[b]['receipt_child']=r
                returned[b]['receipt_birth_step']=receipt_birth[b][r]
                returned[b]['direct_receipt_contact']=bool(set(pair)&receipt_roots[b])
                returned[b]['proof_is_first_receipt']=(r==receipt[b]['child'])
            for b in receivers:
                assert (marker in knowledge[b])==(receipt[b] is not None)
                if receipt[b] is not None:
                    assert bool(receipt_roots[b]&knowledge[source])==(returned[b] is not None)
            tau+=len(W);sf+=e['forced_served'];sv+=e['relief_removed'];sc+=e['cost']
            assert F==F0+sc-sf-sv and P==P0+2*sf-sv and tau+sf==3*t
            if keep_trace:global_trace.append(dict(e,step=t,F_after=F,P_after=P))
            if t in checkpoints:save(t)
        for z in range(len(st)):
            expected=int(marker in st.ancestors[z])
            for k,b in enumerate(receivers,1):
                if receipt_roots[b]&st.ancestors[z]:expected|=1<<k
            assert tags[z]==expected
        assert tuple(st.parents[:len(base)])==tuple(base.parents)
        allowed={'local_tick','new_children','partners','newly_known_parent_cards'}
        assert all(set(card)==allowed for log in local_logs.values() for card in log)
        return dict(seed=seed,clocks=clocks,receipt=receipt,returned=returned,
            receipt_compensator=a_receipt,return_compensator=a_return,
            receipt_predictable_variance=v_receipt,return_predictable_variance=v_return,
            receipt_roots={b:tuple(sorted(v)) for b,v in receipt_roots.items()},
            initial_local_parent_cards=initial_cards,auditor_handle_crosswalk=alias,
            initial_frames=initial_frames,final_frames={o:tuple(sorted(v)) for o,v in knowledge.items()},
            local_logs=local_logs,checkpoints=tuple(saved),
            final_n=len(st),final_F=F,final_P=P,maintenance=tau,construction_work=sc,
            final_parents=tuple(st.parents) if keep_trace else None,
            global_auditor_trace=tuple(global_trace))

    # Existing reflection s=(0 1): compare one persistent trit in two entry frames.
    # No arbitrary message is encoded and no register value is resampled.
    def reflect(g):return 1-g if g in (0,1) else 2
    frame_checks=[]
    for g in range(3):
        for source_entry in (0,1):
            for receiver_entry in (0,1):
                at_source=reflect(g) if source_entry else g
                at_receiver=reflect(g) if receiver_entry else g
                transported=reflect(at_receiver) if source_entry!=receiver_entry else at_receiver
                assert transported==at_source
                frame_checks.append((g,source_entry,receiver_entry,transported))
    assert len(frame_checks)==12

    results={};total_updates=0
    for name,H in environments:
        results[name]={}
        for arm in arms:
            old=old29['environments'][name]['arms'][arm]
            base=construct(old['post_parentage'])
            F0,P0=int(old['post_F']),int(old['post_P'])
            marker=len(base)-1;source=int(old['pair'][0])
            assert tuple(base.parents[marker])==tuple(old['pair']) and source not in receivers
            assert not (set(base.parents[marker])&set(receivers))
            assert marker not in base.recorded and source in base.recorded
            exact=exact_initial(base,F0,marker,source)
            histories=tuple(run(base,F0,P0,H,source,marker,s,keep_trace=(i==0))
                            for i,s in enumerate(seeds))
            total_updates+=len(seeds)*horizon
            summaries={}
            for b in receivers:
                got=[h for h in histories if h['receipt'][b] is not None]
                back=[h for h in histories if h['returned'][b] is not None]
                ticks=[h['returned'][b]['source_ticks'] for h in back]
                rec_hazard=stats([int(h['receipt'][b] is not None)-h['receipt_compensator'][b] for h in histories])
                ret_hazard=stats([int(h['returned'][b] is not None)-h['return_compensator'][b] for h in histories])
                summaries[b]=dict(receipts=len(got),returns=len(back),total_histories=len(histories),
                    not_received_by_horizon=len(histories)-len(got),not_returned_by_horizon=len(histories)-len(back),
                    receipt_fraction=len(got)/len(histories),return_fraction=len(back)/len(histories),
                    direct_receipts=sum(h['receipt'][b]['direct_marker_contact'] for h in got),
                    direct_returns=sum(h['returned'][b]['direct_receipt_contact'] for h in back),
                    returns_proving_first_receipt=sum(h['returned'][b]['proof_is_first_receipt'] for h in back),
                    source_tick_roundtrip_on_completed_only=dict(n=len(ticks),
                        minimum=min(ticks) if ticks else None,
                        median=float(np.median(ticks)) if ticks else None,
                        maximum=max(ticks) if ticks else None),
                    receiver_ticks_at_receipt_completed_only=stats([h['receipt'][b]['receiver_ticks'][b] for h in got]),
                    reception_full_burst_work_completed_only=stats([h['receipt'][b]['complete_burst_cost'] for h in got]),
                    return_full_burst_work_completed_only=stats([h['returned'][b]['complete_burst_cost'] for h in back]),
                    final_source_minus_receiver_ticks=stats([h['clocks'][source]-h['clocks'][b] for h in histories]),
                    receipt_compensator_residual=rec_hazard,return_compensator_residual=ret_hazard,
                    receipt_predictable_RMS_for_sample_mean=sqrt(sum(h['receipt_predictable_variance'][b]
                        for h in histories))/len(histories),
                    return_predictable_RMS_for_sample_mean=sqrt(sum(h['return_predictable_variance'][b]
                        for h in histories))/len(histories),
                    final_receiver_frame_size=stats([len(h['final_frames'][b]) for h in histories]))
            # Independent rerun without any observer code: same native history.
            bare=deepcopy(base);bf,bp=F0,P0;rr=Random(seeds[0])
            for _ in range(horizon):bf,bp,_=native_step(bare,bf,bp,rr,H)
            first=histories[0]
            assert (len(bare),bf,bp)==(first['final_n'],first['final_F'],first['final_P'])
            assert tuple(bare.parents)==first['final_parents']
            results[name][arm]=dict(source=source,marker=marker,receivers=receivers,H=H,
                emission_parents=tuple(base.parents[marker]),emission_cost=old['event']['cost'],
                initial_parentage=tuple(base.parents),initial_F=F0,initial_P=P0,
                exact_initial=exact,histories=histories,summary=summaries)
    # A fixed additional seed block was specified AFTER a conspicuous primary
    # receipt-compensator discrepancy. It is a diagnostic replication, not a
    # replacement panel, a changed rule, or a success-dependent stopping scheme.
    cp=results['no_relief']['recorded_pair']
    replication_seeds=tuple(20310920+i for i in range(256))
    base=construct(cp['initial_parentage'])
    replication_histories=tuple(run(base,cp['initial_F'],cp['initial_P'],0,
        cp['source'],cp['marker'],seed) for seed in replication_seeds)
    replication_summary={}
    for b in receivers:
        residual=[int(h['receipt'][b] is not None)-h['receipt_compensator'][b]
                  for h in replication_histories]
        replication_summary[b]=dict(
            receipts=sum(h['receipt'][b] is not None for h in replication_histories),
            returns=sum(h['returned'][b] is not None for h in replication_histories),
            mean_integrated_receipt_hazard=float(np.mean([h['receipt_compensator'][b]
                                                          for h in replication_histories])),
            receipt_compensator_residual=stats(residual),
            predictable_RMS_for_sample_mean=sqrt(sum(h['receipt_predictable_variance'][b]
                for h in replication_histories))/len(replication_histories))
    replication=dict(environment='no_relief',arm='recorded_pair',seeds=replication_seeds,
        reason='fixed fresh diagnostic block after large primary receipt-compensator discrepancy',
        histories=replication_histories,summary=replication_summary,
        replaces_primary=False,new_native_updates=len(replication_seeds)*horizon)
    assert old29==previous
    print('CELL 30 — EVENT-LOCAL READERS: PULSE ANCESTRY -> RECEIPT -> RETURN')
    print('Native global pool, all enabled births, original workload and relief; no new force or routing law.')
    print('Readers: fixed anchors, own maintenance clocks, incident-event ancestry cards with private handles only.')
    print('Added access convention: importing participating ancestry is a readout, not a gate-cost derivation.')
    print('Pulse = existing Cell29 emission child; its actual preparation backlog is retained.')
    print('Recipient 3 primary, recipient 4 transfer. Tags audit ancestry, not chosen-bit encoding.')
    print(f'{len(seeds)} histories x 4 settings x {horizon} updates = {total_updates:,} new native updates.')
    print('\n setting/arm                    pool     first contact per local tick    B3 receipts/returns    B4 receipts/returns')
    for name,H in environments:
        for arm in arms:
            r=results[name][arm];ex=r['exact_initial'];s3,s4=(r['summary'][b] for b in receivers)
            print(f' {name:12s}/{arm:14s} {ex["total"]:6d} {str(ex["arrival_per_receiver_tick"]):>17s}'
                  f' = {float(ex["arrival_per_receiver_tick"]):.7f}   '
                  f'{s3["receipts"]:3d}/{s3["returns"]:<3d}              {s4["receipts"]:3d}/{s4["returns"]:<3d}')
            print('   B3 source-clock round trip, COMPLETED cases only:',s3['source_tick_roundtrip_on_completed_only'])
            z=s3['final_source_minus_receiver_ticks']
            print(f'   all histories: mean source ticks - B3 ticks = {z["mean"]:+.5f} +/- {z["se"]:.5f} MC SE.')
            print('   marker-to-old-reader graph-distance strata:',ex['available_receiver_parent_child_distance_counts'])
    print('\nExact: P(first contact | receiver served, C carriers) = 1-C(R-1-C,q-1)/C(R-1,q-1).')
    print('One initial carrier: 2/(R-1); single-clock service probability is 3/R.')
    print('Every return: m in Anc(receipt), receipt in Anc(return), with strictly ordered births.')
    print('Fixed old ancestry never acquires a new pulse; event-grown reader frames sometimes do.')
    print('No spatial travel-time bound: equal current carrier counts give equal hazard at every distance.')
    print('Censored runs retained; completed-only times are NOT unconditional mean lifetimes.')
    print('All 12 persistent-trit/entry-frame cases agree after frame transport.')
    print('Fixed additional no-relief/recorded-parent seed block (primary panel NOT replaced):')
    for b,row in replication_summary.items():
        rr=row['receipt_compensator_residual']
        print(f'  B{b}: {row["receipts"]}/256 receipts, {row["returns"]}/256 returns; '
              f'mean(receipt - hazard)={rr["mean"]:+.6f} +/- {rr["se"]:.6f} MC SE.')
    print(f'Additional diagnostic updates: {replication["new_native_updates"]:,}.')
    print('PASS: exact weighted sampling; full ledgers; ancestry provenance; observer-free replays.')
    print('Local clock means have equal conditional drift. No redshift, SI clock, or photon speed assigned.')
    return dict(protocol='v1 event-local ancestry acquisition with native pulse/receipt/return witnesses',
        horizon=horizon,seeds=seeds,checkpoints=checkpoints,receivers=receivers,primary_receiver=3,
        environments=results,total_primary_native_updates=total_updates,
        diagnostic_replication=replication,total_new_native_updates=total_updates+replication['new_native_updates'],
        frame_transport_checks=frame_checks,
        native_transition_changed=False,observer_access_is_added=True,observer_memory_cost_derived=False,
        packet_is_ancestry_provenance_not_chosen_bit=True,work_completion_latency_defined=False,
        timestamp_exchange_compiled=False,spatial_frame_or_propagation_identified=False,
        calibration=None,gravitational_redshift_identified=False,neutron_search_modified=False,
        helpers=dict(conditional_contact=conditional_contact,run=run,construct=construct,exact_initial=exact_initial))


dcu_mass_30 = _run_dcu_mass_30()
