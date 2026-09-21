# CELL 3 — explain Block 9; audit the same X on all 32 matched chain probes.
# Run after Cells 1 and 2. Standard library only; no coefficients or units change.
from fractions import Fraction
from functools import lru_cache
from itertools import product


def _run_dcu_mass_3():
    required = ('dcu_mass_1', 'dcu_mass_2', '_panel', '_parents', '_term', '_Native')
    missing = [name for name in required if name not in globals()]
    if missing:
        raise RuntimeError('Run Cells 1 and 2 first. Missing: ' + ', '.join(missing))
    unit = dcu_mass_1['unit']
    assert unit == dcu_mass_2['unit'] == Fraction(90)
    motifs = _panel['motifs']

    def mean(values):
        values = list(values)
        if not values:
            raise ValueError('Cannot average an empty comparison group.')
        return sum(values, Fraction(0)) / len(values)

    @lru_cache(None)
    def ancestry(v):
        return (frozenset([v]) if v < 2 else
                frozenset([v]) | ancestry(_parents[v][0]) | ancestry(_parents[v][1]))

    @lru_cache(None)
    def chains(v):
        return 1 if v < 2 else sum(chains(p) for p in _parents[v])

    @lru_cache(None)
    def profile(i):
        """Candidate-only contributions, BEFORE overlap subtraction; not a mass."""
        support = set(motifs[i]['support'])
        ports = support - {0, 1}
        recorded = {p for z in support if z >= 2 for p in _parents[z]}
        assert recorded <= support
        assert all(ancestry(u) <= support for u in ports)
        rows = {}
        for z in sorted(ports):
            L = 2 * chains(z) - 2
            d = sum(z in ancestry(u) for u in ports)  # Counts ports, not paths.
            q = 2 if z in recorded else 11
            rows[z] = dict(term=_term(z), L=L, d=d, q=q,
                           k=Fraction(q * L * d, len(ports)))
        return rows

    @lru_cache(None)
    def probe_terms(word):
        """Inclusive ancestry of p0={a,b}; pj={p(j-1), word[j-1]}."""
        p = '(a|b)'
        result = {'a', 'b', p}
        for letter in word:
            p = '(' + '|'.join(sorted((p, letter))) + ')'
            result.add(p)
        return frozenset(result)

    def predict(i, word):
        P = probe_terms(word)
        return sum((r['k'] for r in profile(i).values() if r['term'] not in P),
                   Fraction(0))

    # Audit EVERY saved Cell 1/2 response and every individual port increment.
    saved = dcu_mass_1['rows'] + dcu_mass_2['rows']
    for old in saved:
        i = old['archive_motif_id']
        word = 'aaaaa' if old['probe'] == 'A' else 'bbbbb'
        assert predict(i, word) == old['X']
        P = probe_terms(word)
        for arm in old['ports']:
            predicted_port = sum(r['q'] * r['L'] for z, r in profile(i).items()
                                 if z in ancestry(arm['cache_port']) and r['term'] not in P)
            assert predicted_port == arm['increment']

    # Original exact-C blocks, with the same within-block averaging as Cell 2.
    blocks = sorted((r for r in dcu_mass_2['matched_contrasts'] if r['probe'] == 'A'),
                    key=lambda r: r['block'])
    ids = sorted({i for b in blocks for i in b['G_ids'] + b['control_ids']})
    words = [''.join(w) for w in product('ab', repeat=5)]

    # Freeze algebraic predictions for all contexts BEFORE native comparison calls.
    planned = {(i, w): predict(i, w) for i in ids for w in words}

    def native_response(i, word):
        """Same preparation and counterfactual arms as Cell 1; generalized probe word."""
        s, mapping = _Native(), {0: 0, 1: 1}
        for old in sorted(motifs[i]['support']):
            if old >= 2:
                born, _ = s.add_batch([tuple(mapping[p] for p in _parents[old])])
                mapping[old] = born[0]
        p = s.pair_to_id[(0, 1)]
        for letter in word:
            pair = tuple(sorted((p, 'ab'.index(letter))))
            if pair not in s.pair_to_id:
                s.add_batch([pair])
            p = s.pair_to_id[pair]
        height = [0, 0]
        terms = ['a', 'b']
        for pair in s.parents[2:]:
            height.append(1 + max(height[v] for v in pair))
            terms.append('(' + '|'.join(sorted(terms[v] for v in pair)) + ')')
        assert (height[p], len(s.ancestors[p]), s.path_weight(p)) == (6, 8, 12)
        assert {terms[v] for v in s.ancestors[p]} == probe_terms(word)
        assert set(mapping.values()) | s.ancestors[p] == set(range(len(s)))
        assert s.recorded == {v for pair in s.parents if pair is not None for v in pair}
        before = (len(s), frozenset(s.recorded))
        baseline = s.recording_cost([(0, p)])
        increments = {}
        for old in sorted(motifs[i]['support']):
            if old >= 2:
                u = mapping[old]
                assert u != p and tuple(sorted((u, p))) not in s.pair_to_id
                increments[old] = s.recording_cost([(u, p)]) - baseline
        assert before == (len(s), frozenset(s.recorded))
        return dict(X=mean(increments.values()), increments=increments, baseline=baseline)

    measured = {}
    for key, prediction in planned.items():
        measured[key] = native_response(*key)
        assert measured[key]['X'] == prediction, (key, prediction, measured[key])
        i, word = key
        for u, value in measured[key]['increments'].items():
            expected = sum(r['q'] * r['L'] for z, r in profile(i).items()
                           if z in ancestry(u) and r['term'] not in probe_terms(word))
            assert value == expected

    print(f'PASS: exact decomposition of {len(saved)} saved responses and all their ports.')
    print(f'PASS: {len(ids)} motifs x {len(words)} chain probes = {len(measured)} native checks;')
    print(f'      {len(ids) * 30} checks use the 30 previously untested mixed chain words.')
    print(f'All probes: grade 6, support 8, L=12. Original reporting unit remains {unit}.')

    focus = next(b for b in blocks if b['block'] == 9)
    assert focus['G_ids'] == (4855,) and focus['control_ids'] == (12628,)
    g, c = focus['G_ids'][0], focus['control_ids'][0]
    kg, kc = profile(g), profile(c)
    topology = {}
    print('\nBLOCK 9: actual cores and edge witnesses (IDs are cache labels).')
    for i in (g, c):
        m = motifs[i]
        topology[i] = dict(core=m['core'], edge_marks=m['edge_marks'])
        print(f"  {i}: core={m['core']}, S/I/G={m['full_SIG']}")
        for mark in m['edge_marks']:
            z = mark['child']
            r = profile(i)[z]
            label = mark['sector'] or 'nonregistry'
            print(f"    {tuple(mark['pair'])} -> {z}: {label}, "
                  f"L={r['L']}, q={r['q']}, d={r['d']}, k={r['k']}")

    print('\nShared ancestors with different k (G-bearing minus G-free):')
    differences = []
    for z in sorted(set(kg) & set(kc)):
        delta = kg[z]['k'] - kc[z]['k']
        if delta:
            in_A = _term(z) in probe_terms('aaaaa')
            in_B = _term(z) in probe_terms('bbbbb')
            differences.append(dict(cache_id=z, delta_k=delta, in_A=in_A, in_B=in_B))
            print(f"  {z}: {_term(z)}, L={kg[z]['L']}, "
                  f"d={kg[z]['d']} vs {kc[z]['d']}, delta k={delta}; "
                  f'in A/B={in_A}/{in_B}')
    total_g = sum((r['k'] for r in kg.values()), Fraction(0))
    total_c = sum((r['k'] for r in kc.values()), Fraction(0))
    first_g = sum((r['k'] for r in kg.values() if r['q'] == 11), Fraction(0))
    first_c = sum((r['k'] for r in kc.values() if r['q'] == 11), Fraction(0))
    print(f'  First-use contributions: {first_g} vs {first_c}; difference={first_g-first_c}.')
    print(f'  Whole support sums before overlap removal: {total_g} vs {total_c}.')
    print('  These sums are bookkeeping, NOT a new calibrated observable.')

    summaries = []
    print('\nALL NINE BLOCKS: same X, all 32 chain probes, no averaging away signs.')
    print(' block   min delta    max delta    positive/zero/negative')
    for b in blocks:
        deltas = {}
        for w in words:
            delta = (mean(measured[i, w]['X'] for i in b['G_ids']) -
                     mean(measured[i, w]['X'] for i in b['control_ids']))
            deltas[w] = delta
            if w in ('aaaaa', 'bbbbb'):
                probe = 'A' if w == 'aaaaa' else 'B'
                old = next(r for r in dcu_mass_2['matched_contrasts']
                           if (r['block'], r['probe']) == (b['block'], probe))
                assert delta == old['delta_X']
        signs = (sum(x > 0 for x in deltas.values()), sum(x == 0 for x in deltas.values()),
                 sum(x < 0 for x in deltas.values()))
        summaries.append(dict(block=b['block'], deltas=deltas, signs=signs))
        print(f" {b['block']:3d} {str(min(deltas.values())):>12s} "
              f'{str(max(deltas.values())):>12s}    {signs}')
    block9 = next(r for r in summaries if r['block'] == 9)
    assert all(delta == (Fraction(8, 9) if w[0] == 'a' else Fraction(-8, 9))
               for w, delta in block9['deltas'].items())
    print('\nBlock 9: delta=+8/9 for all 16 a-first probes; -8/9 for all 16 b-first probes.')
    print('This is an ancestry-overlap effect, not a reversal of G recording charge.')
    print('The 32 probes exhaust this chain family, NOT all legal probe shapes.')
    return dict(unit=unit, saved_response_checks=len(saved), predictions=planned,
                native=measured, profiles={i: profile(i) for i in ids},
                block_summaries=summaries, block9_topology=topology,
                block9_common_differences=differences)


dcu_mass_3 = _run_dcu_mass_3()
