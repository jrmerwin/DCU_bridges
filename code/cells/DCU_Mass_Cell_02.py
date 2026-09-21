# CELL 2 — grade-five G-bearing/G-free comparison; Cell 1's X is unchanged.
# Run after Cell 1 in the SAME notebook. No new files or packages required.
from collections import defaultdict
from fractions import Fraction


def _run_dcu_mass_2():
    needed = ('dcu_mass_1', '_panel', '_measure', '_parents', '_term')
    absent = [name for name in needed if name not in globals()]
    if absent:
        raise RuntimeError('Run Cell 1 first. Missing: ' + ', '.join(absent))
    unit = dcu_mass_1['unit']
    if unit != Fraction(90) or dcu_mass_1['unit_source'] != 'K4-1/A':
        raise ValueError('Cell 1 reporting unit changed; expected K4-1/A = 90.')

    motifs = _panel['motifs']
    # These definitions use structure only, before any new response is measured.
    def canonical_key(i):
        return tuple(sorted(_term(v) for v in motifs[i]['core']))

    ids = {}
    for family in ('K3', 'K4'):
        report = next(r for r in _panel['round_reports']
                      if (r['round'], r['family']) == (12, family))
        ids[family] = sorted(
            (i for i in report['unique_motif_ids']
             if (motifs[i]['support_size'], motifs[i]['birth_round']) == (11, 5)),
            key=canonical_key
        )
    assert len(ids['K3']) == 281 and len(ids['K4']) == 1

    # Stronger control: exactly the SAME constituent-ancestry set C_M.
    # Core vertices/pairings can differ; this is NOT insertion/deletion of one G.
    blocks = defaultdict(lambda: {False: [], True: []})
    for i in ids['K3']:
        m = motifs[i]
        blocks[tuple(m['core_support'])][bool(m['full_SIG'][2])].append(i)
    matched = [(core, group) for core, group in blocks.items()
               if group[False] and group[True]]
    matched.sort(key=lambda item: tuple(sorted(_term(v) for v in item[0])))
    matched_g = {i for _, group in matched for i in group[True]}
    unmatched_g = [i for i in ids['K3']
                   if motifs[i]['full_SIG'][2] and i not in matched_g]

    # Check Cell 1 still replays exactly without changing its saved results.
    for old in dcu_mass_1['rows']:
        replay = _measure(motifs[old['archive_motif_id']], 'AB'.index(old['probe']))
        assert replay['X'] == old['X'] and replay['parts'] == old['parts']
        assert replay['ports'] == old['ports']

    rows = []
    for family in ('K3', 'K4'):
        for i in ids[family]:
            m = motifs[i]
            # Validate the archived constituent-ancestry field independently.
            seen, todo = set(), list(m['core'])
            while todo:
                v = todo.pop()
                if v in seen:
                    continue
                seen.add(v)
                if _parents[v] is not None:
                    todo.extend(_parents[v])
            assert seen == set(m['core_support'])
            g_locations = [x for x in m['registry_locations'] if x['sector'] == 'G']
            assert len(g_locations) == m['full_SIG'][2]
            assert m['core_SIG'][2] == 0
            assert all(x['root'] in m['edge_only'] for x in g_locations)
            near = sum(x['closure'] == 'near' for x in g_locations)
            skip = sum(x['closure'] == 'skip' for x in g_locations)
            assert (near, skip) == tuple(m['G_near_skip'])
            assert near + skip == len(g_locations)

            # Prediction BEFORE probing: G witnesses are childless here,
            # absent from both chain probes, and hit only at their own port.
            # The SAME native first-use rule applies to every sector.
            expected_G = Fraction(11 * (16 * near + 14 * skip), 9)
            for branch, probe in enumerate('AB'):
                measured = _measure(m, branch)  # Unmodified Cell 1 function.
                port_map = {p['cache_port']: p for p in measured['ports']}
                assert all(not port_map[x['root']]['target_recorded']
                           for x in g_locations)
                assert measured['parts']['G'] == expected_G
                assert sum(measured['parts'].values()) == measured['X']
                rows.append(dict(
                    measured, archive_motif_id=i, family=family, probe=probe,
                    G_count=len(g_locations), near=near, skip=skip,
                    expected_G=expected_G, X_over_original_unit=measured['X'] / unit
                ))

    lookup = {(r['archive_motif_id'], r['probe']): r for r in rows}
    def mean(values):
        values = list(values)
        return sum(values, Fraction(0)) / len(values) if values else None

    # Each matched C_M block gets equal weight; controls within a block are averaged.
    # Shared motifs/ancestry mean these are NOT independent particle statistics.
    contrasts = []
    for b, (core, group) in enumerate(matched, 1):
        for probe in 'AB':
            positive = [lookup[i, probe] for i in group[True]]
            negative = [lookup[i, probe] for i in group[False]]
            x_g, x_0 = (mean(r['X'] for r in group_rows)
                        for group_rows in (positive, negative))
            g_part = mean(r['parts']['G'] for r in positive)
            delta = x_g - x_0
            contrasts.append(dict(
                block=b, probe=probe, core_support=core,
                G_ids=tuple(group[True]), control_ids=tuple(group[False]),
                X_G_bearing=x_g, X_G_free=x_0, delta_X=delta,
                direct_G_part=g_part, other_difference=delta - g_part
            ))

    print('GRADE 5 / SUPPORT 11 / SAME NATIVE X AND SAME A/B PROBES')
    print(f'Original reporting unit remains {unit}; no mass calibration or refit.')
    print(f"K3: {sum(bool(motifs[i]['full_SIG'][2]) for i in ids['K3'])} G-bearing, "
          f"{sum(not motifs[i]['full_SIG'][2] for i in ids['K3'])} G-free.")
    print(f"K4: {len(ids['K4'])} control; "
          f"{sum(bool(motifs[i]['full_SIG'][2]) for i in ids['K4'])} G-bearing.")
    print(f'Exact-C matched blocks: {len(matched)}; '
          f'{len(matched_g)} G-bearing motifs; '
          f'{sum(len(g[False]) for _, g in matched)} G-free controls.')
    print(f'G-bearing motifs without an exact-C control: {unmatched_g}')

    print('\nWhole fixed stratum (descriptive, not the matched comparison):')
    print(' family/group   probe    n       mean X      X/90       mean G part')
    for family, has_G in (('K3', False), ('K3', True), ('K4', False)):
        for probe in 'AB':
            subset = [r for r in rows if r['family'] == family
                      and bool(r['G_count']) == has_G and r['probe'] == probe]
            if not subset:
                continue
            x = mean(r['X'] for r in subset)
            g = mean(r['parts']['G'] for r in subset)
            print(f" {family}/{'G+' if has_G else 'G-free':7s}  {probe:>3s} "
                  f'{len(subset):4d} {float(x):12.6f} {float(x/unit):10.6f} {float(g):15.6f}')

    print('\nPrimary comparison: G-bearing minus G-free, SAME constituent ancestry.')
    print(' block   G+ archive IDs    G-free archive IDs      delta A     delta B')
    for b, (_, group) in enumerate(matched, 1):
        pair = [r for r in contrasts if r['block'] == b]
        print(f" {b:3d} {str(group[True]):>16s} {str(group[False]):>22s} "
              f"{float(pair[0]['delta_X']):12.6f} {float(pair[1]['delta_X']):11.6f}")

    print('\nEqual-block means: delta X = direct G charge + difference elsewhere.')
    for probe in 'AB':
        subset = [r for r in contrasts if r['probe'] == probe]
        if not subset:
            print(f'  {probe}: no matched comparison; undefined, not zero.')
            continue
        delta = mean(r['delta_X'] for r in subset)
        g_part = mean(r['direct_G_part'] for r in subset)
        other = mean(r['other_difference'] for r in subset)
        signs = tuple(sum((r['delta_X'] > 0, r['delta_X'] == 0,
                           r['delta_X'] < 0)[j] for r in subset) for j in range(3))
        print(f'  {probe}: {delta} = {g_part} + ({other}); '
              f'delta/90={float(delta/unit):+.6f}; positive/zero/negative={signs}')

    print('\nPASS: Cell 1 replay; unchanged native cost; exact G attribution.')
    print('X minus its G-attributed part is bookkeeping, NOT a G-deletion experiment.')
    print('No G-specific coefficient, force law, service clock, or mass factor was added.')
    return dict(protocol='Cell 1 X; D=11,h=5; exact constituent-ancestry matched K3',
                unit=unit, unit_source=dcu_mass_1['unit_source'], rows=rows,
                selected_ids=ids, matched_contrasts=contrasts,
                unmatched_G_ids=unmatched_g)


dcu_mass_2 = _run_dcu_mass_2()
