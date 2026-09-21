# CELL 5 — one electron anchor, unchanged repeat-use Q, and a K5 transfer test.
# Run after Cell 4. Standard library only; no new files or network calls required.
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product


def _run_dcu_mass_5():
    required = ('dcu_mass_1', 'dcu_mass_4', '_panel', '_parents', '_term', '_reference')
    missing = [k for k in required if k not in globals()]
    if missing:
        raise RuntimeError('Run Cells 1–4 first. Missing: ' + ', '.join(missing))
    unit = dcu_mass_1['unit']
    assert unit == dcu_mass_4['unit'] == Fraction(90)
    measure = dcu_mass_4['measure']  # Unchanged native, committed-use protocol.
    motifs = _panel['motifs']
    words = tuple(''.join(w) for w in product('ab', repeat=5))

    # FIX SELECTION BEFORE MEASUREMENT. Reuse the original four candidates.
    selected = [(r['name'], r['archive_motif_id']) for r in dcu_mass_1['rows']
                if r['probe'] == 'A']
    assert [name for name, _ in selected] == ['K3-1', 'K3-2', 'K4-1', 'K4-2']
    report = next(r for r in _panel['round_reports']
                  if (r['round'], r['family']) == (12, 'K5'))
    smallest = min(motifs[i]['support_size'] for i in report['unique_motif_ids'])
    k5_ids = sorted((i for i in report['unique_motif_ids']
                    if motifs[i]['support_size'] == smallest),
                   key=lambda i: (motifs[i]['birth_round'],
                                  tuple(sorted(_term(v) for v in motifs[i]['core']))))
    assert smallest == 15 and len(k5_ids) == 3  # Original archive, not a fitted cut.
    selected += [(f'K5-{j+1}', i) for j, i in enumerate(k5_ids)]
    anchor_key = (next(i for name, i in selected if name == 'K3-1'), 'aaaaa')

    # Published central values from the 2022 CODATA listing, checked 2026-09-19.
    # https://physics.nist.gov/cuu/Constants/Table/allascii.txt
    # Fractions preserve the printed decimals; physical measurements are NOT exact.
    electron_MeV = Fraction('0.51099895069')  # m_e c^2, MeV; uncertainty 0.00000000016.
    muon_ratio = Fraction('206.7682827')      # m_mu / m_e; uncertainty 0.0000046.

    @lru_cache(None)
    def ancestry(v):
        if v < 2:
            return frozenset([v])
        return frozenset([v]) | ancestry(_parents[v][0]) | ancestry(_parents[v][1])

    @lru_cache(None)
    def chains(v):
        return 1 if v < 2 else sum(chains(p) for p in _parents[v])

    def probe_terms(word):
        p, terms = '(a|b)', {'a', 'b', '(a|b)'}
        for letter in word:
            p = '(' + '|'.join(sorted((p, letter))) + ')'
            terms.add(p)
        return terms

    # Architectural annotation, not additional workload weights.
    # Verify the RMR core operator's spectrum exactly; do NOT call it native motion.
    spectral = {}
    for n in (3, 4, 5):
        plus = Fraction(1, n-1)
        def transfer(v):
            return [-sum(v[j] for j in range(n) if j != i) / Fraction(n-1)
                    for i in range(n)]
        assert transfer([1]*n) == [-1]*n
        for j in range(n-1):
            v = [0]*n
            v[j], v[-1] = 1, -1
            assert transfer(v) == [plus*x for x in v]
        spectral[n] = dict(lambda_common=Fraction(-1), lambda_difference=plus,
                           difference_multiplicity=n-1, f_RMR=1-plus)
    registry_size = len(_reference['rows'])
    assert registry_size == 137
    rmr_muon_reference = registry_size * spectral[5]['f_RMR'] / spectral[3]['f_RMR']

    catalogue = []
    for name, i in selected:
        m = motifs[i]
        n = len(m['core'])
        assert m['family'] == f'K{n}' and m['birth_round'] <= 5
        assert {tuple(sorted(e)) for e in m['required_edges']} == set(combinations(sorted(m['core']), 2))
        assert len(m['edge_children']) == n*(n-1)//2
        for edge, child in zip(m['required_edges'], m['edge_children']):
            assert set(_parents[child]) == set(edge)
        core_support = set().union(*(ancestry(v) for v in m['core']))
        assert core_support | set(m['edge_children']) == set(m['support'])
        assert sum(m['full_SIG']) == len(m['registry_locations'])
        catalogue.append(dict(name=name, archive_motif_id=i, core=m['core'],
                              support_size=m['support_size'], grade=m['birth_round'],
                              SIG=tuple(m['full_SIG']), f_RMR=spectral[n]['f_RMR']))

    # Measure ONLY the predetermined reference, then freeze the conversion.
    cases = {anchor_key: measure(*anchor_key)}
    Q_anchor = cases[anchor_key]['Q']
    if Q_anchor <= 0:
        raise ValueError('The declared anchor is nonpositive; calibration is undefined. No replacement selected.')
    conversion = electron_MeV / Q_anchor  # Conditional MeV of rest-energy proxy per Q.
    for name, i in selected:
        for word in words:
            key = (i, word)
            if key not in cases:
                cases[key] = measure(i, word)
            result = cases[key]
            # Independent ancestry-sum audit of actual repeat-use measurements.
            ports = [v for v in motifs[i]['support'] if v >= 2]
            P = probe_terms(word)
            expected = Fraction(2 * sum(2*chains(z)-2 for u in ports
                                       for z in ancestry(u) if _term(z) not in P), len(ports))
            assert result['Q'] == expected
            assert all(s['X'] == expected for s in result['passes'][1:])
            result['Q_over_anchor'] = result['Q'] / Q_anchor
            result['rest_energy_proxy_MeV'] = conversion * result['Q']

    summaries = []
    for row in catalogue:
        i = row['archive_motif_id']
        values = [cases[i, w]['Q_over_anchor'] for w in words]
        summaries.append(dict(row, Q_A=cases[i, 'aaaaa']['Q'], Q_B=cases[i, 'bbbbb']['Q'],
                              ratio_A=cases[i, 'aaaaa']['Q_over_anchor'],
                              ratio_B=cases[i, 'bbbbb']['Q_over_anchor'],
                              ratio_min=min(values), ratio_max=max(values)))
    k5_ratios = [cases[i, w]['Q_over_anchor'] for i in k5_ids for w in words]
    k5_range = (min(k5_ratios), max(k5_ratios))
    contacts = sum(len(s['contacts']) for r in cases.values() for s in r['passes'])

    print('CELL 5: CONDITIONAL Q-ONLY MASS RATIO; ARCHITECTURAL BENCHMARK SEPARATE')
    print('Selection uses the old K3/K4 candidates and ALL smallest-support sampled K5s.')
    for r in catalogue:
        print(f"  {r['name']}: archive {r['archive_motif_id']}, core={r['core']}, "
              f"D={r['support_size']}, h={r['grade']}, S/I/G={r['SIG']}")
    print(f'\nFrozen electron anchor: K3-1/aaaaa; Q_e={Q_anchor}.')
    print(f'Electron rest energy used: {float(electron_MeV):.11f} MeV.')
    print('The anchor is not reset for probe B, mixed chains, or other instances.')
    print(f'Old reporting unit {unit} is unchanged and cancels from Q/Q_e.')
    print('\n case      Q(A)       Q(B)       ratio A    ratio B    full 32-probe ratio range')
    for r in summaries:
        print(f" {r['name']:6s} {str(r['Q_A']):>9s} {str(r['Q_B']):>10s} "
              f"{float(r['ratio_A']):11.6f} {float(r['ratio_B']):10.6f} "
              f"[{float(r['ratio_min']):.6f}, {float(r['ratio_max']):.6f}]")
    print('\nSeparate references (NOT multiplied into Q):')
    for n in (3, 4, 5):
        s = spectral[n]
        print(f"  K{n}: RMR T eigenvalues -1 (x1), {s['lambda_difference']} "
              f"(x{n-1}); declared f={s['f_RMR']}.")
    print(f'  Paper Eq. (9): 137*f5/f3 = {rmr_muon_reference} = {float(rmr_muon_reference):.6f}.')
    print(f'  CODATA muon/electron comparison target = {float(muon_ratio):.7f}.')
    print(f'  Measured K5 Q/Q_e envelope = [{float(k5_range[0]):.6f}, {float(k5_range[1]):.6f}].')
    print(f'  Conditional K5 rest-energy envelope = '
          f'[{float(electron_MeV*k5_range[0]):.6f}, {float(electron_MeV*k5_range[1]):.6f}] MeV.')
    print(f'  Fractional errors vs muon target = '
          f'[{100*float(k5_range[0]/muon_ratio-1):.6f}%, {100*float(k5_range[1]/muon_ratio-1):.6f}%].')
    print(f'\nPASS: {len(cases)} native histories; {contacts:,} committed contacts; exact repeat-Q audit.')
    print('K4 is a structural control, NOT assigned a physical vacuum mass.')
    print('This tests the chosen Q-only correspondence, NOT registry-coupled resonance or proton confinement.')
    print('Literal registry counts are not imposed active-channel counts; no 137/136, cubic, or spectral multiplier was applied.')
    return dict(protocol='Cell 4 repeat Q; original K3/K4 plus smallest-support U12 K5',
                unit=unit, anchor=dict(key=anchor_key, Q=Q_anchor,
                                      electron_rest_energy_MeV=electron_MeV,
                                      conversion_MeV_per_Q=conversion),
                source_constants=dict(adjustment='CODATA 2022', checked='2026-09-19',
                                      muon_electron_ratio=muon_ratio),
                catalogue=catalogue, summaries=summaries, cases=cases,
                K5_ratio_range=k5_range, native_contacts=contacts,
                RMR_reference=dict(spectral=spectral, muon_ratio_eq9=rmr_muon_reference,
                                   status='paper prescription; NOT a native Q factor'))


dcu_mass_5 = _run_dcu_mass_5()
