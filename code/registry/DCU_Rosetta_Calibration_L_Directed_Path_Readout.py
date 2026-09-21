# Calibration L — a directed-path readout of the existing reference registry.
# No replacement transfer matrix, signal normalization, growth, clock, or mass law.
# x is a FORMAL marker for edge count. Polynomial coefficients are path counts.
from collections import Counter, defaultdict
import hashlib
import pickle


def rosetta_l_poly_add(left, right):
    result = [0] * max(len(left), len(right))
    for i, value in enumerate(left):
        result[i] += value
    for i, value in enumerate(right):
        result[i] += value
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def rosetta_l_poly_string(coefficients):
    terms = []
    for power, coefficient in enumerate(coefficients):
        if coefficient == 0:
            continue
        variable = '' if power == 0 else ('x' if power == 1 else f'x^{power}')
        magnitude = abs(coefficient)
        body = str(magnitude) if power == 0 else (
            variable if magnitude == 1 else f'{magnitude}*{variable}')
        if not terms:
            terms.append(('-' if coefficient < 0 else '') + body)
        else:
            terms.append((' - ' if coefficient < 0 else ' + ') + body)
    return ''.join(terms) or '0'


def rosetta_registry_directed_paths(exact, placement):
    """Path-length-resolved Ch on ALL 137 entries, with matched near/skip probes.

    P_a^a=1, P_a^b=0; P_b^a=0, P_b^b=1.
    P_{u,v}^s(x)=x*(P_u^s(x)+P_v^s(x)). No path is deduplicated merely
    because it meets another path. Summing coefficients recovers native Ch.
    G labels are consulted only AFTER all profiles and equality classes exist.
    This is a fixed-DAG readout, NOT a physical impulse-response experiment.
    """
    ref = exact['reference']
    anatomy = placement['G_anatomy']
    before = hashlib.sha256(pickle.dumps((ref, anatomy), protocol=5)).digest()
    parents = tuple(None if p is None else tuple(sorted(p)) for p in ref['parents'])
    rows = ref['rows']
    if parents[:2] != (None, None) or len(rows) != 137:
        raise ValueError('Use the existing complete two-primitive registry reference.')
    if any(r['index'] != i for i, r in enumerate(rows)):
        raise ValueError('Reference indices must retain their current ordering.')

    pa, pb, chains = [(1,), (0,)], [(0,), (1,)], [1, 1]
    ancestors = [frozenset({0}), frozenset({1})]
    seen_pairs = set()
    for z, pair in enumerate(parents[2:], 2):
        if pair is None or len(pair) != 2 or not 0 <= pair[0] < pair[1] < z or pair in seen_pairs:
            raise ValueError('Malformed, non-topological, or duplicate parent pair.')
        seen_pairs.add(pair)
        u, v = pair
        pa.append((0,) + rosetta_l_poly_add(pa[u], pa[v]))
        pb.append((0,) + rosetta_l_poly_add(pb[u], pb[v]))
        chains.append(chains[u] + chains[v])
        ancestors.append(ancestors[u] | ancestors[v] | {z})
    total = tuple(rosetta_l_poly_add(a, b) for a, b in zip(pa, pb))
    if tuple(ancestors) != tuple(ref['ancestors']):
        raise AssertionError('Parent reconstruction disagrees with stored reference ancestry.')
    if any(sum(poly) != chains[z] for z, poly in enumerate(total)):
        raise AssertionError('Path-length count does not recover native Ch.')

    # Independent explicit route enumeration: intentionally do not use a visited set.
    routes_checked = 0
    for z in range(len(parents)):
        expected = (Counter(), Counter())
        pending = [(z, 0)]
        while pending:
            node, length = pending.pop()
            if parents[node] is None:
                expected[node][length] += 1
                routes_checked += 1
            else:
                pending.extend((p, length + 1) for p in parents[node])
        actual = tuple({i: c for i, c in enumerate(p[z]) if c} for p in (pa, pb))
        if actual != tuple(dict(c) for c in expected):
            raise AssertionError('Explicit routes disagree with polynomial recurrence.')

    total_classes, named_classes = defaultdict(list), defaultdict(list)
    records = []
    for r in rows:
        z, i = r['object_id'], r['index']
        if len(ancestors[z]) != 7:
            raise ValueError('The census requires all 137 size-seven entries.')
        total_classes[total[z]].append(i)
        named_classes[(pa[z], pb[z])].append(i)
        records.append({'index': i, 'label': r['label'], 'sector': r['sector'],
                        'object_id': z, 'from_a': pa[z], 'from_b': pb[z],
                        'total': total[z], 'Ch': chains[z], 'L': 2 * chains[z] - 2})

    by_choices = defaultdict(dict)
    g_ids = {r['index'] for r in rows if r['sector'] == 'G'}
    for a in anatomy:
        key, closure, i = tuple(a['primitive_choices']), a['closure'], a['index']
        if i not in g_ids or closure not in ('near', 'skip') or closure in by_choices[key]:
            raise ValueError('Malformed existing near/skip anatomy.')
        by_choices[key][closure] = i
    if len(by_choices) != 8 or {i for v in by_choices.values() for i in v.values()} != g_ids:
        raise ValueError('Keep all eight pairs and all sixteen G entries.')

    pairs = []
    for choices, ids in sorted(by_choices.items()):
        if set(ids) != {'near', 'skip'}:
            raise ValueError('Both closures are required for each primitive-choice triple.')
        near, skip = records[ids['near']], records[ids['skip']]
        neg_skip = tuple(-c for c in skip['total'])
        difference = rosetta_l_poly_add(near['total'], neg_skip)
        pairs.append({'choices': choices, 'near': near['index'], 'skip': skip['index'],
                      'near_label': near['label'], 'skip_label': skip['label'],
                      'near_total': near['total'], 'skip_total': skip['total'],
                      'difference': difference, 'near_Ch': near['Ch'], 'skip_Ch': skip['Ch'],
                      'near_L': near['L'], 'skip_L': skip['L']})

    census = []
    for poly, ids in total_classes.items():
        if not set(ids) & g_ids:
            continue
        census.append({'total': poly, 'indices': tuple(ids),
                       'counts_SIG': tuple(sum(rows[i]['sector'] == s for i in ids)
                                           for s in ('S', 'I', 'G'))})
    resolved = []
    for i in sorted(g_ids):
        r = records[i]
        matching = named_classes[(r['from_a'], r['from_b'])]
        others = tuple(j for j in matching if j not in g_ids)
        resolved.append({'index': i, 'label': r['label'], 'from_a': r['from_a'],
                         'from_b': r['from_b'], 'nonG_matches': others,
                         'nonG_labels': tuple(rows[j]['label'] for j in others)})

    print('DIRECTED PATH READOUT — same registry and G pairs; no new evolution rule')
    print('P^s_{u,v}(x)=x*(P^s_u(x)+P^s_v(x)); x marks edge count, NOT time or frequency.')
    print('All directed primitive-rooted paths count, including different routes through shared ancestors.')
    print('\nALL EIGHT ANATOMICAL PAIRS')
    print(f"{'choices':>8} {'near/skip':>12} {'Ch near/skip':>14} {'L near/skip':>13} {'Pnear - Pskip':>22}")
    for p in pairs:
        bits = ''.join('ab'[v] for v in p['choices'])
        print(f"{bits:>8} {p['near_label']+'/'+p['skip_label']:>12} "
              f"{str((p['near_Ch'],p['skip_Ch'])):>14} {str((p['near_L'],p['skip_L'])):>13} "
              f"{rosetta_l_poly_string(p['difference']):>22}")
    print('\nTOTAL PROFILE CLASSES — all 137 entries grouped BEFORE consulting sector labels')
    for c in census:
        print(f"  P(x)={rosetta_l_poly_string(c['total'])}; Ch={sum(c['total'])}; "
              f"matching S/I/G={c['counts_SIG']}")
    print('\nPRIMITIVE-RESOLVED PROFILES — exact equality of BOTH polynomials, fixed a/b labels')
    for r in resolved:
        print(f"  {r['label']}: a=[{rosetta_l_poly_string(r['from_a'])}], "
              f"b=[{rosetta_l_poly_string(r['from_b'])}]; non-G equals={r['nonG_labels']}")
    mixed = sum(bool(r['nonG_matches']) for r in resolved)
    non_g = {j for r in resolved for j in r['nonG_matches']}
    print(f'  G entries with a non-G identical resolved profile: {mixed}/16; '
          f'distinct non-G matches={len(non_g)}.')
    print(f'Independent check: {len(parents)} node profiles, {routes_checked} explicit primitive-rooted routes.')
    print('Ch and L are native path counts, NOT newly calibrated masses or energies.')
    print('Total counts lose route identity; even the primitive-resolved profile need not recover the full DAG.')
    print('No sectors were reassigned, no favourable subset promoted, no transfer weights or targets fitted.')
    print('No new samples, lineage links, recording obligations, service, or RNG draws. Sources unchanged.')
    after = hashlib.sha256(pickle.dumps((ref, anatomy), protocol=5)).digest()
    if after != before:
        raise AssertionError('An input reference or anatomy record was changed.')
    return {'scope': 'exact primitive-rooted path-length readout, not physical dynamics',
            'polynomial_variable': 'formal edge-count marker', 'records': records,
            'near_skip_pairs': pairs, 'total_profile_census': census,
            'primitive_resolved_census': resolved, 'mixed_G_entries': mixed,
            'nonG_resolved_matches': tuple(sorted(non_g)),
            'node_profiles_verified': len(parents), 'explicit_routes_verified': routes_checked}


if not all(k in globals() for k in ('rosetta_constructor_reference', 'rosetta_G_structure_placement')):
    raise RuntimeError('Keep the original Calibration-H reference and Calibration-J anatomy in this kernel.')
rosetta_registry_path_readout = rosetta_registry_directed_paths(
    rosetta_constructor_reference, rosetta_G_structure_placement)
