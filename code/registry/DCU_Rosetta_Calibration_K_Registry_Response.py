# Calibration K — one declared RMR-inspired response on the FULL registry graph.
# Observer-side linear probe only. No native clock, graph, workload, or RNG is changed.
# State vectors are signed signals, NOT probabilities, masses, or quantum amplitudes.
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
import hashlib
import pickle


def rosetta_k_transfer(adjacency, signal):
    """T=-A D^{-1}: each source sends minus its signal / its degree to each neighbor.

    This is a proposed extension of RMR's complete-graph rule, not a derived DCU
    update. Column/source normalization is explicit. Never renormalize the output.
    """
    if len(adjacency) != len(signal):
        raise ValueError('Graph and signal sizes disagree.')
    result = [Fraction(0) for _ in adjacency]
    for source, value in enumerate(signal):
        if not value:
            continue
        if not adjacency[source]:
            raise ValueError('This rule is undefined at an isolated active source.')
        transmitted = -Fraction(value) / len(adjacency[source])
        for target in adjacency[source]:
            result[target] += transmitted
    return tuple(result)


def rosetta_k_pair_probe(adjacency, left, right):
    """Two probes with equal initial L1 norm 2: +1/-1 and +1/+1 on a named pair."""
    n = len(adjacency)
    if not 0 <= left < right < n:
        raise ValueError('Require two different valid vertices in increasing order.')
    difference = [Fraction(0) for _ in range(n)]
    common = [Fraction(0) for _ in range(n)]
    difference[left], difference[right] = Fraction(1), Fraction(-1)
    common[left] = common[right] = Fraction(1)
    y = rosetta_k_transfer(adjacency, difference)
    common_y = rosetta_k_transfer(adjacency, common)
    outside = [i for i in range(n) if i not in (left, right)]
    signed_outside_l1 = sum(abs(y[i]) for i in outside)
    eigenvalue = y[left] if all(
        y[i] == y[left] * difference[i] for i in range(n)) else None
    return {
        'difference_input': tuple(difference), 'difference_output': y,
        'common_input': tuple(common), 'common_output': common_y,
        'difference_eigenvalue': eigenvalue,
        'difference_outside_l1_over_input_l1': signed_outside_l1 / 2,
        'common_outside_l1_over_input_l1': sum(abs(common_y[i]) for i in outside) / 2,
        'difference_norm_ratio': sum(abs(a) for a in y) / 2,
    }


def rosetta_registry_response(exact, placement, response_steps=(0, 1, 2, 4, 8)):
    """Same eight near/skip pairs. Full 137-vertex response plus all matched-degree pairs.

    Neighborhood equivalence classes are discovered WITHOUT S/I/G labels. Sector
    labels are applied afterward. Exact rational arithmetic; no fitted coefficient.
    Graph restriction/renormalization to G, K8, or K3 is not performed.
    """
    if placement.get('scope') != 'read-only G structure and placement on the existing exact-closure observations':
        raise ValueError('Use the unchanged Calibration-J result.')
    steps = tuple(response_steps)
    if not steps or steps[0] != 0 or tuple(sorted(set(steps))) != steps or any(
            type(k) is not int or not 0 <= k <= 16 for k in steps):
        raise ValueError('Use increasing, unique response indices from 0 through 16, including 0.')
    ref = exact['reference']
    anatomy = placement['G_anatomy']
    before = hashlib.sha256(pickle.dumps((ref, anatomy), protocol=5)).digest()
    rows, ancestors = ref['rows'], ref['ancestors']
    n = len(rows)
    if n != 137 or [r['index'] for r in rows] != list(range(n)):
        raise ValueError('Use the complete 137-entry reference with its existing indices.')
    supports = tuple(frozenset(ancestors[r['object_id']]) for r in rows)
    if any(len(s) != 7 for s in supports):
        raise ValueError('Every reference support must have inclusive size seven.')
    # Reconstruct independently of the stored edge array and sector labels.
    adjacency = tuple(frozenset(j for j in range(n) if j != i and
        len(supports[i] & supports[j]) > 3) for i in range(n))
    if adjacency != tuple(ref['adjacency']):
        raise ValueError('Stored overlap graph disagrees with the actual reference ancestries.')
    degrees = tuple(map(len, adjacency))

    # Equal CLOSED neighborhoods identify adjacent twins without using sector names.
    neighborhood_classes = defaultdict(list)
    for i in range(n):
        neighborhood_classes[adjacency[i] | {i}].append(i)
    blocks = tuple(tuple(v) for v in neighborhood_classes.values())
    block_of = {i: k for k, block in enumerate(blocks) for i in block}
    block_report = []
    for k, block in enumerate(blocks):
        if any(j not in adjacency[i] for i, j in combinations(block, 2)):
            raise AssertionError('Closed-neighborhood class is not internally complete.')
        block_report.append({
            'block': k, 'indices': block,
            'labels': tuple(rows[i]['label'] for i in block),
            'size': len(block), 'degree': degrees[block[0]],
            'counts_SIG': tuple(sum(rows[i]['sector'] == s for i in block) for s in ('S','I','G')),
        })
    block_connections = []
    for k, a in enumerate(blocks):
        line = []
        for h, b in enumerate(blocks):
            if k == h:
                line.append('clique')
            else:
                counts = [len(adjacency[i] & set(b)) for i in a]
                if set(counts) not in ({0}, {len(b)}):
                    raise AssertionError('Between-block adjacency is not all-or-none.')
                line.append('all' if counts[0] else 'none')
        block_connections.append(tuple(line))

    by_choices = defaultdict(dict)
    for entry in anatomy:
        i = entry['index']
        key = tuple(entry['primitive_choices'])
        if rows[i]['sector'] != 'G' or entry['closure'] in by_choices[key]:
            raise ValueError('Malformed Calibration-J anatomy.')
        by_choices[key][entry['closure']] = i
    g_indices = {i for i,r in enumerate(rows) if r['sector'] == 'G'}
    if len(by_choices) != 8 or {i for v in by_choices.values() for i in v.values()} != g_indices:
        raise ValueError('Keep all sixteen G entries in their eight matched closure pairs.')
    pair_rows = []
    for choices, variants in sorted(by_choices.items()):
        if set(variants) != {'near','skip'}:
            raise ValueError('Each chain must have one near and one skip closure.')
        near, skip = variants['near'], variants['skip']
        proper_near = supports[near] - {rows[near]['object_id']}
        proper_skip = supports[skip] - {rows[skip]['object_id']}
        if len(proper_near) != 6 or proper_near != proper_skip:
            raise AssertionError('The matched roots do not have the same six proper ancestors.')
        lo, hi = sorted((near, skip))
        probe = rosetta_k_pair_probe(adjacency, lo, hi)
        common_external = tuple(sum(abs(probe['common_output'][j]) for j in range(n)
            if j not in (lo,hi) and rows[j]['sector'] == s) / 2 for s in ('S','I','G'))
        history, signal = [], probe['difference_input']
        for k in range(steps[-1]+1):
            if k in steps:
                outside = sum(abs(signal[j]) for j in range(n) if j not in (lo,hi))
                history.append({'response_step': k, 'left': signal[lo], 'right': signal[hi],
                    'outside_l1': outside, 'l1_over_initial': sum(abs(x) for x in signal)/2})
            if k != steps[-1]:
                signal = rosetta_k_transfer(adjacency, signal)
        ev = probe['difference_eigenvalue']
        # Independently compare the iterated response with its measured eigenpair.
        if ev is not None and any(h['left'] != ev**h['response_step'] or
                                 h['right'] != -ev**h['response_step'] or h['outside_l1']
                                 for h in history):
            raise AssertionError('Explicit repeated response disagrees with its eigenpair.')
        pair_rows.append({'choices': choices, 'near': near, 'skip': skip,
            'near_label': rows[near]['label'], 'skip_label': rows[skip]['label'],
            'degree': degrees[near], 'same_block': block_of[near] == block_of[skip],
            'common_external_SIG_over_input_l1': common_external,
            **probe, 'history': tuple(history)})

    # Control: ALL pairs of vertices with the SAME degree as the anatomical G pairs.
    g_degrees = {degrees[i] for i in g_indices}
    if len(g_degrees) != 1:
        raise AssertionError('Matched-degree control requires one G degree in this reference.')
    degree = next(iter(g_degrees))
    pool = [i for i,d in enumerate(degrees) if d == degree]
    controls, control_types = [], Counter()
    for a,b in combinations(pool,2):
        probe = rosetta_k_pair_probe(adjacency,a,b)
        outside = probe['difference_outside_l1_over_input_l1']
        controls.append({'pair': (a,b), 'sector_pair': tuple(sorted((rows[a]['sector'],rows[b]['sector']))),
            'eigenvalue': probe['difference_eigenvalue'], 'outside_l1_over_input_l1': outside})
        if outside == 0:
            control_types[''.join(sorted((rows[a]['sector'],rows[b]['sector'])))] += 1
    zero_leak = sum(c['outside_l1_over_input_l1'] == 0 for c in controls)

    print('REGISTRY RESPONSE — one DECLARED RMR-inspired probe, not native DCU time evolution')
    print('T_ij = -A_ij / degree(j); use the FULL 137-entry unweighted overlap graph.')
    print('Columns are source-normalized; no post-step renormalization, added edge, or fitted factor.')
    print('Signed signals are not probabilities or energies. Response step k is NOT a maintenance tick.')
    print('\nLABEL-BLIND CLOSED-NEIGHBORHOOD CLASSES (sectors counted only after grouping)')
    for b in block_report:
        print(f"  block={b['block']}: size={b['size']}; degree={b['degree']}; S/I/G={b['counts_SIG']}")
    print('Between-block adjacency:', tuple(block_connections))
    print('\nALL EIGHT ANATOMICAL NEAR/SKIP PAIRS')
    print(f"{'choices':>8} {'near / skip':>13} {'degree':>6} {'lambda-':>10} {'diff outside':>13} {'common outside':>15}")
    for p in pair_rows:
        bits=''.join('ab'[v] for v in p['choices'])
        print(f"{bits:>8} {p['near_label']+'/'+p['skip_label']:>13} {p['degree']:6d} "
              f"{str(p['difference_eigenvalue']):>10} "
              f"{str(p['difference_outside_l1_over_input_l1']):>13} "
              f"{str(p['common_outside_l1_over_input_l1']):>15}")
    print('External common-signal fractions by S/I/G:',
          sorted({tuple(map(str,p['common_external_SIG_over_input_l1'])) for p in pair_rows}))
    print('\nREPEATED DIFFERENCE RESPONSE — every pair has the displayed result when identical')
    identical = len({tuple((h['response_step'],h['left'],h['right'],h['outside_l1'],h['l1_over_initial'])
                         for h in p['history']) for p in pair_rows}) == 1
    print('  histories identical across all anatomical pairs:', identical)
    for p in pair_rows[:1] if identical else pair_rows:
        if not identical: print(' pair=',p['near_label'],p['skip_label'])
        for h in p['history']:
            print(f"  k={h['response_step']:2d}: gain={h['l1_over_initial']} "
                  f"(~{float(h['l1_over_initial']):.9g}); outside signal={h['outside_l1']}")
    print('\nMATCHED-DEGREE CONTROL — every pair, not a random or selected subset')
    print(f'  degree={degree}; eligible vertices={len(pool)}; pair differences={len(controls)}')
    print(f'  zero external difference response={zero_leak}/{len(controls)}')
    print('  zero-external sector-pair counts:', dict(sorted(control_types.items())))
    print('Same response in other sectors is retained; this is not an independence or significance test.')
    print('Local cancellation does not imply non-decaying amplitude, a physical bound state, or a mass.')
    print('The 1/2 factor, 17, 27, and 1836 are not inputs. The response rule itself is a stated hypothesis.')
    print('All native graphs, samples, lineages, work counters and random states remain unchanged.')
    if hashlib.sha256(pickle.dumps((ref, anatomy), protocol=5)).digest() != before:
        raise AssertionError('A source dictionary or anatomy record changed.')
    return {'scope': 'declared signed degree-normalized response on fixed full registry',
        'operator': 'T_ij=-A_ij/degree(j)', 'adjacency': adjacency, 'degrees': degrees,
        'neighborhood_classes': tuple(block_report), 'block_connections': tuple(block_connections),
        'near_skip_pairs': tuple(pair_rows), 'response_steps': steps,
        'matched_degree_controls': tuple(controls), 'zero_external_controls': zero_leak,
        'zero_external_control_sector_counts': dict(control_types)}


if 'rosetta_constructor_reference' not in globals() or 'rosetta_G_structure_placement' not in globals():
    raise RuntimeError('Keep the existing Calibration-H reference and Calibration-J result in this kernel.')
rosetta_registry_response_results = rosetta_registry_response(
    rosetta_constructor_reference, rosetta_G_structure_placement)
