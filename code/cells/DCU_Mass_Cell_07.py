# CELL 7 — nested-closure demand, without a threshold or resonance multiplier.
# Run after Cell 5; Cell 6 is a separate branch and is NOT required.
# Added experiment: a closure-verification schedule. Native costs are unchanged.
# Compare one shared certificate per witness vs a separate request per clique.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb


def _run_dcu_mass_7():
    needed = ('dcu_mass_5', '_panel', '_parents', '_term', '_Native')
    missing = [key for key in needed if key not in globals()]
    if missing:
        raise RuntimeError('Run Cells 1-5 first. Missing: ' + ', '.join(missing))
    catalogue = dcu_mass_5['catalogue']
    motifs = _panel['motifs']

    @lru_cache(None)
    def ancestry(v):
        if v < 2:
            return frozenset([v])
        return frozenset([v]) | ancestry(_parents[v][0]) | ancestry(_parents[v][1])

    def binary_rank(vectors):
        """Exact triangle-edge incidence rank over GF(2), not physical DOF."""
        basis = {}
        for value in vectors:
            while value:
                pivot = value.bit_length() - 1
                if pivot in basis:
                    value ^= basis[pivot]
                else:
                    basis[pivot] = value
                    break
        return len(basis)

    def inventory(i):
        m = motifs[i]
        support = set(m['support'])
        pairs = {tuple(sorted(_parents[v])): v for v in support if v >= 2}
        vertices = sorted(support - {0, 1})
        cliques = []
        # Enumerate the COMPLETE internally witnessed clique inventory in support,
        # not appearances in larger sampled hosts. No primitive core vertices.
        for k in range(3, len(vertices) + 1):
            if comb(k, 2) > len(pairs):
                break
            for core in combinations(vertices, k):
                edges = tuple(combinations(core, 2))
                if all(edge in pairs for edge in edges):
                    witnesses = tuple(pairs[edge] for edge in edges)
                    local_support = set().union(*(ancestry(v) for v in core))
                    local_support.update(witnesses)
                    assert local_support <= support
                    cliques.append(dict(core=core, witnesses=witnesses))
        assert any(set(c['core']) == set(m['core']) for c in cliques)
        counts = Counter(w for c in cliques for w in c['witnesses'])
        edges = sorted(tuple(sorted(e)) for e in m['required_edges'])
        edge_index = {edge: j for j, edge in enumerate(edges)}
        vectors = [sum(1 << edge_index[e] for e in combinations(t, 2))
                   for t in combinations(sorted(m['core']), 3)]
        rank = binary_rank(vectors)
        assert rank == len(edges) - len(m['core']) + 1
        return dict(cliques=cliques, counts=dict(counts),
                    by_order=dict(Counter(len(c['core']) for c in cliques)),
                    core_cycle_rank=rank)

    def experiment(i, info, branch):
        m = motifs[i]
        if m['birth_round'] > 5:
            raise ValueError('This fixed-probe protocol requires grade <= 5.')
        s, mapping = _Native(), {0: 0, 1: 1}
        preparation = 0
        for old in sorted(m['support']):
            if old >= 2:
                born, charge = s.add_batch([tuple(mapping[v] for v in _parents[old])])
                mapping[old] = born[0]
                preparation += charge
        p = s.pair_to_id[(0, 1)]
        for _ in range(5):
            pair = tuple(sorted((p, branch)))
            if pair not in s.pair_to_id:
                _, charge = s.add_batch([pair])
                preparation += charge
            p = s.pair_to_id[pair]
        witnesses = tuple(sorted(info['counts'], key=_term))
        support = frozenset(mapping.values())
        assert all(tuple(sorted((mapping[u], mapping[v]))) in s.pair_to_id
                   for u, v in m['required_edges'])  # Closures already exist.
        _, warm_work = s.add_batch([(mapping[w], p) for w in witnesses])
        assert support <= s.recorded
        overlap = support & s.ancestors[p]
        original_parents = tuple(s.parents[v] for v in sorted(support))
        costs = {w: 2 * sum(s.path_weight(z) for z in
                           s.ancestors[mapping[w]] - s.ancestors[p]) for w in witnesses}
        policies = {'shared': (witnesses,),
                    'per_clique': tuple(tuple(sorted(c['witnesses'], key=_term))
                                        for c in info['cliques'])}
        results = {}
        for policy, groups in policies.items():
            state, probe = deepcopy(s), p
            passes = []
            for cycle in range(2):
                contacts, extension_work = [], 0
                for group in groups:
                    # Each group gets a genuinely NEW partner; repeated witnesses
                    # across cliques never recreate an already-existing parent pair.
                    born, charge = state.add_batch([(probe, branch)])
                    probe = born[0]
                    extension_work += charge
                    assert support & state.ancestors[probe] == overlap
                    batch_quote = state.recording_cost([(mapping[w], probe) for w in group])
                    group_raw = 0
                    for w in group:
                        baseline = state.recording_cost([(0, probe)])  # Alternative arm.
                        raw = state.recording_cost([(mapping[w], probe)])
                        born, committed = state.add_batch([(mapping[w], probe)])
                        assert committed == raw and raw - baseline == costs[w]
                        group_raw += raw
                        contacts.append(dict(witness=w, raw=raw, baseline=baseline,
                                             net=raw-baseline, child=born[0]))
                    assert group_raw == batch_quote  # Native burst/sequential identity.
                total = sum(c['net'] for c in contacts)
                expected = sum(costs.values()) if policy == 'shared' else sum(
                    costs[w] * info['counts'][w] for w in witnesses)
                assert total == expected
                passes.append(dict(total=total, contacts=contacts,
                                   probe_extension_work_excluded=extension_work))
            assert passes[0]['total'] == passes[1]['total']
            assert original_parents == tuple(state.parents[v] for v in sorted(support))
            results[policy] = dict(total=passes[0]['total'],
                                   contacts_per_cycle=len(passes[0]['contacts']), passes=passes)
        return dict(preparation_work_excluded=preparation, warm_work_excluded=warm_work,
                    witness_costs=costs, policies=results)

    inventories = {r['archive_motif_id']: inventory(r['archive_motif_id']) for r in catalogue}
    cases = {(r['archive_motif_id'], label): experiment(r['archive_motif_id'],
             inventories[r['archive_motif_id']], branch)
             for r in catalogue for label, branch in (('A', 0), ('B', 1))}
    # Separate response references for TWO declared schedules, never a MeV refit.
    anchor_id = next(r['archive_motif_id'] for r in catalogue if r['name'] == 'K3-1')
    anchors = {p: cases[anchor_id, 'A']['policies'][p]['total'] for p in ('shared', 'per_clique')}
    if any(v <= 0 for v in anchors.values()):
        raise ValueError('Declared cycle reference is nonpositive; ratios are undefined.')
    summaries = []
    print('CELL 7: INTERNAL CLIQUE INVENTORY AND TWO LEGAL VERIFICATION SCHEDULES')
    print('Candidate selection unchanged. All cliques in each complete support are enumerated.')
    print('Counts are distinct structures, NOT their appearances in sampled hosts.')
    print(' case    K3/K4/K5 counts   unique witnesses  clique-witness visits  core cycle rank')
    for r in catalogue:
        i = r['archive_motif_id']
        inv = inventories[i]
        print(f" {r['name']:6s} {str(tuple(inv['by_order'].get(k,0) for k in (3,4,5))):>16s}"
              f" {len(inv['counts']):17d} {sum(inv['counts'].values()):22d} {inv['core_cycle_rank']:16d}")
    print('\nTotal repeat-work per cycle (NOT per-port Q):')
    print(' case     shared A/B       per-clique A/B     ratio shared A/B    ratio per-clique A/B')
    for r in catalogue:
        i = r['archive_motif_id']
        totals = {label: {p: cases[i,label]['policies'][p]['total'] for p in anchors}
                  for label in ('A','B')}
        ratios = {label: {p: Fraction(totals[label][p], anchors[p]) for p in anchors}
                  for label in ('A','B')}
        print(f" {r['name']:6s} {totals['A']['shared']:5d}/{totals['B']['shared']:<5d}"
              f" {totals['A']['per_clique']:8d}/{totals['B']['per_clique']:<8d}"
              f" {float(ratios['A']['shared']):9.4f}/{float(ratios['B']['shared']):.4f}"
              f" {float(ratios['A']['per_clique']):12.4f}/{float(ratios['B']['per_clique']):.4f}")
        summaries.append(dict(name=r['name'], archive_motif_id=i, totals=totals, ratios=ratios))
    contact_count = sum(len(stage['contacts']) for c in cases.values()
                        for p in c['policies'].values() for stage in p['passes'])
    print(f'\nPASS: {contact_count:,} committed verification contacts; two equal repeat cycles per schedule.')
    print('Fresh legal partners, unchanged candidate parentage and probe overlap, exact native costs.')
    print('Both schedules are PROPOSED verification protocols, not mandatory autonomous maintenance.')
    print('Repeated clique requests are actual extra operations, not independent copies of shared objects.')
    print('No registry multiplier, trigger threshold, service clock, or mass conversion was added.')
    return dict(protocol='all internal cliques; shared vs per-clique legal witness verification',
                inventories=inventories, cases=cases, summaries=summaries,
                response_anchors=anchors, committed_contacts=contact_count,
                physical_mass_calibration=None)


dcu_mass_7 = _run_dcu_mass_7()
