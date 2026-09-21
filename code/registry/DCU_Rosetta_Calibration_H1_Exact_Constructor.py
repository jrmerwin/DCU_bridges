# Calibration H1 — compare the construction grammar, then enumerate its exhaustive closure.
# NEW reference experiment: all legal old pairs are constructed in synchronous rounds.
# This is NOT the Gamma=2 service process and has NO maintenance-clock assignment.
# All existing runs and RNGs are untouched. Standard library only.
from collections import Counter
from itertools import combinations
from fractions import Fraction
import time


def rosetta_h_choose2(n):
    return n * (n - 1) // 2


def rosetta_h_ids(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def rosetta_h_new_structure():
    return {'parents': [None, None], 'ancestors': [1, 2], 'chains': [1, 1],
            'heights': [0, 0], 'mapped7': [0, 1], 'registry_masks': [0, 0],
            'pairs': {}}


def rosetta_h_append(graph, a, b, reference):
    a, b = sorted((a, b))
    z = len(graph['parents'])
    if not 0 <= a < b < z or (a, b) in graph['pairs']:
        raise ValueError('Require distinct old parents and an unused pair.')
    ancestry = graph['ancestors'][a] | graph['ancestors'][b] | (1 << z)
    size = ancestry.bit_count()
    small = None
    if size <= 7:
        ma, mb = graph['mapped7'][a], graph['mapped7'][b]
        if ma is None or mb is None:
            raise AssertionError('A small child cannot have a large parent.')
        small = reference['pairs'][tuple(sorted((ma, mb)))]
    own = reference['index_by_object'].get(small)
    mask = graph['registry_masks'][a] | graph['registry_masks'][b]
    if own is not None:
        mask |= 1 << own
    graph['parents'].append((a, b))
    graph['ancestors'].append(ancestry)
    graph['chains'].append(graph['chains'][a] + graph['chains'][b])
    graph['heights'].append(1 + max(graph['heights'][a], graph['heights'][b]))
    graph['mapped7'].append(small)
    graph['registry_masks'].append(mask)
    graph['pairs'][a, b] = z
    return z


def rosetta_h_summary(histogram, reference, round_number=None, registry_count=137):
    n = sum(histogram.values())
    sectors = reference['sector_masks']
    return {
        'round': round_number, 'N': n, 'registry_count': registry_count,
        'registry_carriers': sum(c for m, c in histogram.items() if m),
        'registry_incidences': sum(c * m.bit_count() for m, c in histogram.items()),
        'pair_coembeddings': sum(c * rosetta_h_choose2(m.bit_count())
                                for m, c in histogram.items()),
        'sector_carriers': {s: sum(c for m, c in histogram.items() if m & bit)
                            for s, bit in sectors.items()},
        'sector_incidences': {s: sum(c * (m & bit).bit_count()
                                   for m, c in histogram.items())
                              for s, bit in sectors.items()}}


def rosetta_h_rounds5(reference, native_structure_class):
    """Independent construction vs existing DCUStructure on IDENTICAL prescribed batches.

    This tests the constructor and accounting code, not equivalence of schedulers.
    """
    graph = rosetta_h_new_structure()
    native = native_structure_class()
    rows = []
    for r in range(6):
        histogram = Counter(graph['registry_masks'])
        own = sum(a.bit_count() == 7 for a in graph['ancestors'])
        rows.append(rosetta_h_summary(histogram, reference, r, own))
        if r == 5:
            break
        n = len(graph['parents'])
        batch = [p for p in combinations(range(n), 2) if p not in graph['pairs']]
        # Independent first/repeat batch charge, including shared-factor reuse.
        hits = Counter()
        recorded = {x for ps in graph['parents'][2:] for x in ps}
        for a, b in batch:
            hits.update(rosetta_h_ids(graph['ancestors'][a] | graph['ancestors'][b]))
        expected_cost = sum(2 * (2 * graph['chains'][x] - 2) * uses
                            for x, uses in hits.items())
        expected_cost += sum(9 * (2 * graph['chains'][x] - 2)
                             for x in hits if x not in recorded)
        born, actual_cost = native.add_batch(batch)
        expected_born = [rosetta_h_append(graph, a, b, reference) for a, b in batch]
        if tuple(born) != tuple(expected_born) or actual_cost != expected_cost:
            raise AssertionError('Existing constructor disagrees on a prescribed batch.')
        if list(native.parents) != graph['parents'] or list(native.chains) != graph['chains']:
            raise AssertionError('Parentage or chain weights differ.')
        for z, actual in enumerate(native.ancestors):
            if sum(1 << x for x in actual) != graph['ancestors'][z]:
                raise AssertionError('Inclusive ancestry differs.')
        if set(native.recorded) != {x for ps in graph['parents'][2:] for x in ps}:
            raise AssertionError('Recorded masks differ.')
        rows[-1]['next_batch_work'] = actual_cost
    if [row['N'] for row in rows] != [2, 3, 5, 12, 68, 2280]:
        raise AssertionError('Unexpected exhaustive population sequence.')
    return graph, rows


def rosetta_h_stream_round6(graph, reference):
    """Enumerate every U5 pair once; obtain U6 statistics WITHOUT storing 2.6M nodes.

    U6 = two primitives plus the child of EACH unordered pair in U5, including
    children that were already born in earlier rounds. Therefore no pair is skipped.
    """
    hist = Counter({0: 2})
    sizes = Counter({1: 2})
    layer_hist = {1: Counter({0: 2})}
    for a, b in combinations(range(len(graph['parents'])), 2):
        size = (graph['ancestors'][a] | graph['ancestors'][b]).bit_count() + 1
        mask = graph['registry_masks'][a] | graph['registry_masks'][b]
        if size == 7:
            small = reference['pairs'][tuple(sorted((graph['mapped7'][a], graph['mapped7'][b])))]
            mask |= 1 << reference['index_by_object'][small]
        hist[mask] += 1
        sizes[size] += 1
        if size not in layer_hist:
            layer_hist[size] = Counter()
        layer_hist[size][mask] += 1
    return {'histogram': hist, 'size_counts': sizes, 'layer_histograms': layer_hist,
            'summary': rosetta_h_summary(hist, reference, 6)}


def rosetta_h_complete_size_cap(reference, cap=9):
    """Complete UNIVERSAL ancestry-size-bounded catalogue, not a time-truncated run.

    A child's inclusive ancestry is strictly larger than either parent's. Parents
    of size >= cap can never contribute another object of size <= cap. Pruning
    those parents therefore loses NO object in the requested finite catalogue.
    """
    if type(cap) is not int or not 7 <= cap <= 9:
        raise ValueError('This bounded experiment supports caps 7, 8, or 9.')
    graph = rosetta_h_new_structure()
    eligible, pos, pair_checks = [0, 1], 1, 0
    while pos < len(eligible):
        b = eligible[pos]
        for a in eligible[:pos]:
            pair_checks += 1
            size = (graph['ancestors'][a] | graph['ancestors'][b]).bit_count() + 1
            if size > cap:
                continue
            z = rosetta_h_append(graph, a, b, reference)
            if size < cap:
                eligible.append(z)
        pos += 1
    graph['size_cap'] = cap
    graph['pair_checks'] = pair_checks
    graph['size_counts'] = Counter(a.bit_count() for a in graph['ancestors'])
    return graph


def rosetta_h_pair_counts(histogram, sector_mask):
    result = Counter()
    for mask, count in histogram.items():
        for pair in combinations(rosetta_h_ids(mask & sector_mask), 2):
            result[pair] += count
    return result


def rosetta_h_spatial_backbone(round6, reference):
    """Reproduce the PAPER'S finite-snapshot rule, not a new dimensionality test.

    Intersect pair-coembedding supports at ancestry sizes 13,14,15 in U6; average
    their weights; retain >= the linear 95th percentile of positive weights.
    """
    counts = {size: rosetta_h_pair_counts(h, reference['sector_masks']['S'])
              for size, h in round6['layer_histograms'].items()}
    support = set(counts[13]) & set(counts[14]) & set(counts[15])
    totals = {p: sum(counts[s][p] for s in (13, 14, 15)) for p in support}
    values = sorted(totals.values())
    if not values:
        return {'support_edges': 0, 'top_edges': 0, 'top_vertices': (),
                'same_as_nine_reference_S_universals': False, 'pair_counts_by_size': counts}
    position = Fraction(95, 100) * (len(values) - 1)
    lower = position.numerator // position.denominator
    remainder = position - lower
    cut = values[lower] * (1 - remainder) + values[min(lower + 1, len(values) - 1)] * remainder
    top = {p for p, weight in totals.items() if weight >= cut}
    vertices = {v for pair in top for v in pair}
    reference_nine = {r['index'] for r in reference['rows']
                      if r['sector'] == 'S' and len(reference['adjacency'][r['index']]) == 136}
    return {'support_edges': len(support), 'top_edges': len(top), 'top_vertices': tuple(sorted(vertices)),
            'top_is_clique': len(top) == rosetta_h_choose2(len(vertices)),
            'same_as_nine_reference_S_universals': vertices == reference_nine,
            'mean_weight_threshold': cut / 3, 'pair_counts_by_size': counts,
            'persistent_weight_sums': totals, 'top_edge_set': frozenset(top)}


def rosetta_h_exact_marginal_rounds(full5, round6, reference, final_round=10):
    """Exact one- and two-registry-entry incidence counts; no later full graphs built.

    For a fixed entry r after round 5:
      d'_r = 1 + C(N,2) - C(N-d_r,2).
    The +1 is r itself; r's own parents contain no DAG-7 entry.
    For different r,s, inclusion-exclusion gives their pair-coembedding count.
    """
    if type(final_round) is not int or not 6 <= final_round <= 12:
        raise ValueError('Use 6 <= final_round <= 12 for this exact summary calculation.')
    choose = rosetta_h_choose2
    histogram = Counter(full5['registry_masks'])
    n = len(full5['parents'])
    singles = [sum(c for m, c in histogram.items() if m & (1 << i)) for i in range(137)]
    pairs = rosetta_h_pair_counts(histogram, (1 << 137) - 1)
    tags = dict(reference['sector_masks']); tags['any'] = (1 << 137) - 1
    groups = {s: sum(c for m, c in histogram.items() if m & tag) for s, tag in tags.items()}
    out = []
    for r in range(5, final_round + 1):
        row = {'round': r, 'N': n, 'registry_count': 137,
               'registry_carriers': groups['any'], 'registry_incidences': sum(singles),
               'sector_carriers': {s: groups[s] for s in ('S', 'I', 'G')},
               'pair_coembeddings': sum(pairs.values()),
               'per_entry_descendants': tuple(singles), 'pair_coembedding_counts': dict(pairs)}
        out.append(row)
        if r == 6:
            streamed = round6['summary']
            for key in ('N', 'registry_carriers', 'registry_incidences', 'sector_carriers', 'pair_coembeddings'):
                if row[key] != streamed[key]:
                    raise AssertionError('Exact recurrence disagrees with streamed round six.')
            for i, count in enumerate(singles):
                actual = sum(c for m, c in round6['histogram'].items() if m & (1 << i))
                if actual != count:
                    raise AssertionError('Single-entry recurrence disagrees with enumeration.')
            if pairs != rosetta_h_pair_counts(round6['histogram'], (1 << 137) - 1):
                raise AssertionError('Pair-entry recurrence disagrees with enumeration.')
        if r == final_round:
            break
        allpairs = choose(n)
        newpairs = {(i, j): allpairs - choose(n - singles[i]) - choose(n - singles[j])
                    + choose(n - singles[i] - singles[j] + pairs.get((i, j), 0))
                    for i, j in combinations(range(137), 2)}
        singles = [1 + allpairs - choose(n - x) for x in singles]
        groups = {s: tags[s].bit_count() + allpairs - choose(n - value)
                  for s, value in groups.items()}
        pairs, n = newpairs, 2 + allpairs
    return out


def rosetta_h_reference_experiment(reference, native_structure_class):
    started = time.perf_counter()
    full5, rounds = rosetta_h_rounds5(reference, native_structure_class)
    round6 = rosetta_h_stream_round6(full5, reference)
    bounded9 = rosetta_h_complete_size_cap(reference, 9)
    backbone = rosetta_h_spatial_backbone(round6, reference)
    later = rosetta_h_exact_marginal_rounds(full5, round6, reference, final_round=10)
    print('CONSTRUCTOR COMPARISON — independent exhaustive batches vs existing DCUStructure')
    print('PASS: five rounds; 2,280 parent lists, ancestries, chain weights and masks agree; all batch charges agree.')
    print('This is a fresh structural fixture, NOT a change to native Gamma=2 service or your saved runs.')
    print('\nEXHAUSTIVE ROUNDS — six is streamed, not stored as millions of full objects')
    print(f"{'round':>5} {'objects':>14} {'R7 roots':>9} {'any R%':>9} {'any G%':>9} {'mean R/node':>12} {'pair coemb':>13}")
    for row in rounds + [round6['summary']]:
        n = row['N']
        print(f"{row['round']:5d} {n:14d} {row['registry_count']:9d} "
              f"{100*row['registry_carriers']/n:9.5f} {100*row['sector_carriers']['G']/n:9.5f} "
              f"{row['registry_incidences']/n:12.6f} {row['pair_coembeddings']:13d}")
    print('\nCOMPLETE ancestry-size <=9 catalogue:', dict(sorted(bounded9['size_counts'].items())))
    print('Total objects:', len(bounded9['parents']), '; exact size 9:', bounded9['size_counts'][9],
          '; size-9 objects inside round six:', round6['size_counts'][9])
    print('A size bound is not a time/epoch bound. It is exact here because ancestor size strictly decreases.')
    print('\nHISTORICAL finite-round S-sector comparison:')
    print('  size-13/14/15 common coembedding edges:', backbone['support_edges'])
    print('  >=95th-percentile edges/active vertices:', backbone['top_edges'], len(backbone['top_vertices']))
    print('  same nine S registry universal nodes:', backbone['same_as_nine_reference_S_universals'])
    print('  This is concentration within selected ancestry-size slices of round six, not a force or 3D result.')
    print('\nLATER EXHAUSTIVE ROUNDS — exact marginal recurrences; later graphs NOT enumerated')
    print(f"{'round':>5} {'objects':>14} {'R roots/N':>13} {'any R%':>11} {'any G%':>11} {'mean R/node':>12}")
    for row in later:
        n = row['N']
        ns = str(n) if n < 10**14 else f'{n:.5e}'
        print(f"{row['round']:5d} {ns:>14} {137/n:13.5e} {100*row['registry_carriers']/n:11.6f} "
              f"{100*row['sector_carriers']['G']/n:11.6f} {row['registry_incidences']/n:12.6f}")
    print('Round seven would contain', later[2]['N'], 'objects: the unrestricted rule does NOT halt at six.')
    elapsed = time.perf_counter() - started
    print(f'Reference experiment elapsed: {elapsed:.2f}s. No native clocks or random draws assigned to these rounds.')
    return {'reference': reference, 'full5': full5, 'explicit_rounds': rounds,
            'round6': round6, 'bounded9': bounded9, 'spatial_backbone': backbone,
            'later_exact_marginals': later, 'elapsed_seconds': elapsed,
            'constructor_equality_checked_through_round': 5,
            'scope': 'exhaustive closure benchmark; not service-driven dynamics'}


if 'DCUStructure' not in globals() or 'rosetta_registry_anatomy_results' not in globals():
    raise RuntimeError('Use the existing kernel with DCUStructure and Calibration F results.')
rosetta_constructor_reference = rosetta_h_reference_experiment(
    rosetta_registry_anatomy_results['reference'], DCUStructure)
