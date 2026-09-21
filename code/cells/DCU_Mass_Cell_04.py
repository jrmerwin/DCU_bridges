# CELL 4 — first contact versus actual repeated use of a fixed candidate.
# Run after Cells 1–3. Standard library only. No service clock or mass factors.
from copy import deepcopy
from fractions import Fraction
from itertools import product


def _run_dcu_mass_4():
    required = ('dcu_mass_1', 'dcu_mass_2', 'dcu_mass_3',
                '_panel', '_parents', '_term', '_Native')
    missing = [k for k in required if k not in globals()]
    if missing:
        raise RuntimeError('Run Cells 1–3 first. Missing: ' + ', '.join(missing))
    unit = dcu_mass_1['unit']
    assert unit == dcu_mass_2['unit'] == dcu_mass_3['unit'] == Fraction(90)
    motifs = _panel['motifs']
    words = tuple(''.join(w) for w in product('ab', repeat=5))
    blocks = sorted((r for r in dcu_mass_2['matched_contrasts'] if r['probe'] == 'A'),
                    key=lambda r: r['block'])
    ids = sorted({i for b in blocks for i in b['G_ids'] + b['control_ids']})
    assert len(ids) == 22 and len(words) == 32

    def mean(values):
        values = tuple(values)
        if not values:
            raise ValueError('Empty average is undefined.')
        return sum(values, Fraction(0)) / len(values)

    def measure(i, word, passes=4):
        """Fixed ports; commit a complete contact pass, then use fresh probe roots."""
        m = motifs[i]
        if len(word) != 5 or set(word) - {'a', 'b'} or passes < 2:
            raise ValueError('Use a five-letter a/b word and at least two passes.')
        if m['birth_round'] > 5:
            raise ValueError('This probe protocol is restricted to grade <= 5 supports.')
        s, mapping = _Native(), {0: 0, 1: 1}
        setup_work = 0
        for old in sorted(m['support']):
            if old >= 2:
                born, charge = s.add_batch([tuple(mapping[p] for p in _parents[old])])
                mapping[old] = born[0]
                setup_work += charge
        candidate_setup_work = setup_work
        p = s.pair_to_id[(0, 1)]
        for letter in word:
            pair = tuple(sorted((p, 'ab'.index(letter))))
            if pair not in s.pair_to_id:
                _, charge = s.add_batch([pair])
                setup_work += charge
            p = s.pair_to_id[pair]
        ports = sorted((v for v in m['support'] if v >= 2), key=_term)
        support = frozenset(mapping[v] for v in m['support'])
        assert len(support) == m['support_size'] == len(ports) + 2
        assert support | s.ancestors[p] == set(range(len(s)))
        original = tuple((s.parents[v], s.ancestors[v], s.chains[v]) for v in sorted(support))
        overlap = support & s.ancestors[p]
        exclusive = {v: s.ancestors[mapping[v]] - s.ancestors[p] for v in ports}
        repeat_prediction = {v: 2 * sum(s.path_weight(z) for z in exclusive[v]) for v in ports}
        g_nodes = {mapping[r['root']] for r in m['registry_locations'] if r['sector'] == 'G'}
        # The no-use control receives identical probes but no candidate contacts.
        cold, cp = deepcopy(s), p
        interaction_children, log = set(), []
        for k in range(passes):
            extension_work = cold_extension_work = 0
            if k:
                # Probe successors inherit ONLY the previous probe, never its
                # candidate-contact children. Extension respects primitive exchange.
                endpoint = 'ab'.index(word[0])
                born, extension_work = s.add_batch([(p, endpoint)])
                p = born[0]
                born, cold_extension_work = cold.add_batch([(cp, endpoint)])
                cp = born[0]
            assert (len(s.ancestors[p]), s.path_weight(p)) == (8 + k, 12 + 2*k)
            assert (len(cold.ancestors[cp]), cold.path_weight(cp)) == (8 + k, 12 + 2*k)
            assert support & s.ancestors[p] == support & cold.ancestors[cp] == overlap
            assert not (interaction_children & s.ancestors[p])
            contacts = []
            for old in ports:
                u = mapping[old]
                assert s.ancestors[u] - s.ancestors[p] == exclusive[old]
                if k:
                    assert s.ancestors[u] <= s.recorded
                state_before = (len(s), frozenset(s.recorded))
                baseline = s.recording_cost([(0, p)])  # Alternative, NOT committed.
                raw = s.recording_cost([(u, p)])
                net = raw - baseline
                assert state_before == (len(s), frozenset(s.recorded))
                charges = {z: (2 if z in s.recorded else 11) * s.path_weight(z)
                           for z in exclusive[old]}
                assert net == sum(charges.values())
                cold_baseline = cold.recording_cost([(0, cp)])
                cold_raw = cold.recording_cost([(u, cp)])
                newly_recorded = tuple(v for v in ports
                                       if mapping[v] in s.ancestors[u] - s.recorded)
                born, committed = s.add_batch([(u, p)])  # ACTUAL first/repeat use.
                assert committed == raw and len(born) == 1
                interaction_children.add(born[0])
                assert s.ancestors[u] <= s.recorded
                contacts.append(dict(port=old, raw=raw, baseline=baseline, net=net,
                                     cold_raw=cold_raw, cold_baseline=cold_baseline,
                                     cold_net=cold_raw-cold_baseline,
                                     G_part=sum(charges.get(z, 0) for z in g_nodes),
                                     newly_recorded_support=newly_recorded))
            assert support <= s.recorded
            assert s.recorded == {v for pair in s.parents if pair is not None for v in pair}
            assert original == tuple((s.parents[v], s.ancestors[v], s.chains[v])
                                     for v in sorted(support))
            log.append(dict(pass_number=k+1, X=mean(r['net'] for r in contacts),
                            cold_X=mean(r['cold_net'] for r in contacts),
                            G_part=mean(r['G_part'] for r in contacts),
                            extension_work=extension_work,
                            cold_extension_work=cold_extension_work, contacts=contacts))
        first_by_port = {r['port']: r['net'] for r in log[0]['contacts']}
        for k, stage in enumerate(log):
            for r in stage['contacts']:
                assert r['cold_net'] == first_by_port[r['port']]
                if k:
                    assert r['net'] == repeat_prediction[r['port']]
                    assert not r['newly_recorded_support']
        Q = log[1]['X']
        assert all(stage['X'] == Q for stage in log[1:])
        return dict(archive_motif_id=i, family=m['family'], word=word,
                    support_size=len(support), ports=tuple(ports),
                    first_X=log[0]['X'], Q=Q, startup_premium=log[0]['X']-Q,
                    Q_over_original_unit=Q/unit, candidate_setup_work=candidate_setup_work,
                    initial_probe_setup_work=setup_work-candidate_setup_work,
                    committed_work=setup_work + sum(stage['extension_work'] +
                        sum(r['raw'] for r in stage['contacts']) for stage in log),
                    passes=log)

    cases = {}
    for i in ids:
        for word in words:
            result = measure(i, word)
            old = dcu_mass_3['native'][i, word]
            assert result['first_X'] == old['X']
            assert {r['port']: r['net'] for r in result['passes'][0]['contacts']} == old['increments']
            cases[i, word] = result

    contrasts = []
    for b in blocks:
        for word in words:
            groups = [[cases[i, word] for i in b[key]] for key in ('G_ids', 'control_ids')]
            first = mean(r['first_X'] for r in groups[0]) - mean(r['first_X'] for r in groups[1])
            repeat = mean(r['Q'] for r in groups[0]) - mean(r['Q'] for r in groups[1])
            contrasts.append(dict(block=b['block'], word=word, first_delta=first,
                                  repeat_delta=repeat, startup_delta=first-repeat,
                                  Q_G=mean(r['Q'] for r in groups[0]),
                                  Q_control=mean(r['Q'] for r in groups[1])))

    n_contacts = sum(len(stage['contacts']) for r in cases.values() for stage in r['passes'])
    print('FIRST CONTACT -> REPEAT USE; SAME FIXED CANDIDATES AND PER-PORT OBSERVABLE')
    print(f'{len(cases)} histories; 4 passes each; {n_contacts:,} committed candidate contacts.')
    print('First pass exactly replays all 704 previous responses, including every port.')
    print('Passes 2, 3 and 4 agree exactly; all no-use controls retain their first response.')
    print(f'Original reporting unit remains {unit}; no physical conversion or service clock.')
    print('\nG-bearing minus matched G-free:')
    print(' block     first A    repeat A      first B    repeat B')
    for b in blocks:
        a, z = [next(r for r in contrasts if r['block'] == b['block'] and r['word'] == w)
                for w in ('aaaaa', 'bbbbb')]
        print(f" {b['block']:3d} {str(a['first_delta']):>11s} {str(a['repeat_delta']):>11s}"
              f" {str(z['first_delta']):>12s} {str(z['repeat_delta']):>11s}")
    print('\nRepeat contrast across all 32 original probe contexts:')
    print(' block       min       max    positive/zero/negative')
    for b in blocks:
        values = [r['repeat_delta'] for r in contrasts if r['block'] == b['block']]
        signs = (sum(x > 0 for x in values), sum(x == 0 for x in values), sum(x < 0 for x in values))
        print(f" {b['block']:3d} {str(min(values)):>10s} {str(max(values)):>9s}    {signs}")
    print('\nEqual-block mean contrast: first = repeat + first-use premium difference.')
    for word in ('aaaaa', 'bbbbb'):
        subset = [r for r in contrasts if r['word'] == word]
        first, repeat, startup = [mean(r[key] for r in subset)
                                  for key in ('first_delta', 'repeat_delta', 'startup_delta')]
        retained = 'undefined' if not first else f'{100*float(repeat/first):.3f}%'
        print(f'  {word}: {first} = {repeat} + {startup}; fraction retained={retained}')
        print(f"    Matched mean Q: G-bearing={mean(r['Q_G'] for r in subset)}, "
              f"G-free={mean(r['Q_control'] for r in subset)}")
    qvalues = [r['Q'] for r in cases.values()]
    print(f'\nRepeat Q range in this panel: {min(qvalues)} to {max(qvalues)}.')
    print('PASS: legal new pairs, real record-mask updates, fixed candidate ancestry,')
    print('      unchanged probe overlap, exact accounting, and repeat-use identity.')
    print('New interaction children are retained, but are not candidate ports.')
    print('Q is externally elicited work per port, NOT autonomous turnover or mass.')
    return dict(protocol='fixed-support, sequential contact passes; primitive-extended probes',
                unit=unit, passes_per_history=4, cases=cases, contrasts=contrasts,
                measure=measure)


dcu_mass_4 = _run_dcu_mass_4()
