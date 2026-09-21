# CELL 9 — distances from native data; distinguish dimension from metric failure.
# Run after Cells 5 and 8. Cells 6/7 are not required. Standard library only.
# Nine readouts are retained, including controls and unsuccessful conventions.
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
from itertools import combinations


def _run_dcu_mass_9():
    needed = ('dcu_mass_5', 'dcu_mass_8', '_panel', '_parents', '_reference', '_Native')
    missing = [k for k in needed if k not in globals()]
    if missing:
        raise RuntimeError('Run Cells 5 and 8 first. Missing: ' + ', '.join(missing))
    catalogue, motifs, ref = dcu_mass_5['catalogue'], _panel['motifs'], _reference
    audit = dcu_mass_8['audit_distances']  # Exact rational PSD/minor/rank checks.
    rules = (
        ('CP_hops2', 'co-parent shortest-path distance, squared (clique control)'),
        ('PC_hops2', 'parent-child shortest-path distance, squared'),
        ('Jaccard2', 'inclusive-ancestry Jaccard distance, squared'),
        ('Anc_delta', 'ancestry symmetric difference AS squared feature distance'),
        ('Anc_Ldelta', 'L-weighted symmetric difference AS squared feature distance'),
        ('Anc_Ldelta2', 'L-weighted symmetric difference AS length, then squared'),
        ('R_overlap', 'mean squared difference of all-registry overlap counts'),
        ('S_overlap', 'mean squared difference of S-only overlap counts'),
        ('R_hops', 'mean squared difference of distances to registry landmarks'),
    )

    # One common legal host: complete registry + old candidates + their mirrors.
    # This is a finite measurement frame, not the decoded cache as a universe.
    host, rmap = _Native(), {0: 0, 1: 1}
    for v, pair in enumerate(ref['parents'][2:], 2):
        born, _ = host.add_batch([tuple(rmap[p] for p in pair)])
        rmap[v] = born[0]

    @lru_cache(None)
    def embed(v, mirror=False):
        if v < 2:
            return 1-v if mirror else v
        pair = tuple(sorted(embed(p, mirror) for p in _parents[v]))
        if pair not in host.pair_to_id:
            host.add_batch([pair])
        return host.pair_to_id[pair]

    for row in catalogue:
        for v in motifs[row['archive_motif_id']]['support']:
            embed(v)
            embed(v, True)
    n = len(host)
    anc = host.ancestors
    registry = tuple(rmap[r['object_id']] for r in ref['rows'])
    spatial = tuple(rmap[r['object_id']] for r in ref['rows'] if r['sector'] == 'S')
    pc, cp = [set() for _ in range(n)], [set() for _ in range(n)]
    for child, pair in enumerate(host.parents):
        if pair is not None:
            a, b = pair
            pc[child].update(pair)
            pc[a].add(child)
            pc[b].add(child)
            cp[a].add(b)
            cp[b].add(a)
    assert n == 180 and (len(registry), len(spatial)) == (137, 81)

    def bfs(source, graph, allowed=None):
        distances, queue = {source: 0}, deque([source])
        while queue:
            v = queue.popleft()
            for w in graph[v]:
                if w not in distances and (allowed is None or w in allowed):
                    distances[w] = distances[v] + 1
                    queue.append(w)
        return distances

    @lru_cache(None)
    def hops(v, kind):
        return bfs(v, pc if kind == 'pc' else cp)

    @lru_cache(None)
    def features(v, kind):
        if kind == 'R_hops':
            return tuple(hops(v, 'pc')[r] for r in registry)
        roots = spatial if kind == 'S_overlap' else registry
        return tuple(len(anc[v] & anc[r]) for r in roots)

    @lru_cache(None)
    def pair_values(u, v):
        delta = anc[u] ^ anc[v]
        weighted = sum(host.path_weight(z) for z in delta)
        if v not in hops(u, 'cp'):
            raise ValueError('A co-parent distance is disconnected; no finite value substituted.')
        values = dict(
            CP_hops2=hops(u, 'cp')[v]**2,
            PC_hops2=hops(u, 'pc')[v]**2,
            Jaccard2=(1-Fraction(len(anc[u] & anc[v]), len(anc[u] | anc[v])))**2,
            Anc_delta=len(delta), Anc_Ldelta=weighted, Anc_Ldelta2=weighted**2)
        for kind in ('R_overlap', 'S_overlap', 'R_hops'):
            a, b = features(u, kind), features(v, kind)
            values[kind] = Fraction(sum((x-y)**2 for x, y in zip(a, b)), len(a))
        return values

    def status(result):
        return ('E' if result['positive_semidefinite'] else 'N') + str(result['rank'])

    def evaluate(vertices):
        answer = {}
        for kind, _ in rules:
            squared = [[Fraction(pair_values(u, v)[kind]) for v in vertices] for u in vertices]
            result = audit(squared)
            assert all(squared[i][j] > 0 for i, j in combinations(range(len(vertices)), 2))
            result['D2'] = squared
            result['status'] = status(result)
            answer[kind] = result
        return answer

    cases = {}
    for row in catalogue:
        m = motifs[row['archive_motif_id']]
        for edge, child in zip(m['required_edges'], m['edge_children']):
            assert set(_parents[child]) == set(edge)
        vertices = tuple(embed(v) for v in m['core'])
        support = {embed(v) for v in m['support']}
        # Shortest core paths give identical results in the minimal and common frames.
        for v in vertices:
            for kind, graph in (('pc', pc), ('cp', cp)):
                local = bfs(v, graph, support)
                assert all(local[w] == hops(v, kind)[w] for w in vertices)
        result = evaluate(vertices)
        mirrored = evaluate(tuple(embed(v, True) for v in m['core']))
        for kind, _ in rules:
            assert mirrored[kind]['D2'] == result[kind]['D2']
            reversed_result = audit([list(reversed(r)) for r in reversed(result[kind]['D2'])])
            assert status(reversed_result) == result[kind]['status']
        cases[row['name']] = dict(archive_motif_id=row['archive_motif_id'],
                                  cache_core=tuple(m['core']), readouts=result)

    # All non-clique five-point groups in the SAME already-selected eight-vertex pool.
    # These are specificity controls, not independent observations or new particles.
    pool = tuple(sorted({v for row in catalogue for v in row['core']}))
    nulls = []
    for core in combinations(pool, 5):
        missing_edges = tuple((u, v) for u, v in combinations(core, 2)
                              if embed(v) not in cp[embed(u)])
        if missing_edges:
            nulls.append(dict(cache_core=core, missing_edges=missing_edges,
                              readouts=evaluate(tuple(embed(v) for v in core))))
    assert len(pool) == 8 and len(nulls) == 30
    null_counts = {kind: dict(Counter(r['readouts'][kind]['status'] for r in nulls))
                   for kind, _ in rules}
    # A literal inherited-registry-root mask collapses these young cores to one point.
    assert all(not (anc[embed(v)] & set(registry)) for v in pool)

    # For K5, also test all five four-point subsets under each unchanged metric.
    facets = {}
    for name, case in cases.items():
        if len(case['cache_core']) == 5:
            facets[name] = {}
            for kind, _ in rules:
                matrix = case['readouts'][kind]['D2']
                facets[name][kind] = tuple(
                    dict(indices=subset, **audit([[matrix[i][j] for j in subset] for i in subset]))
                    for subset in combinations(range(5), 4))

    # Exact counterexample: the first K4 fails for parent-child hop lengths.
    square = cases['K4-1']['readouts']['PC_hops2']
    vector = (1, -1, -1, 1)
    assert all(sum(square['gram'][i][j]*vector[j] for j in range(4)) == -vector[i]
               for i in range(4))

    print('CELL 9 — NATIVE-DATA METRICS; NO EDGE-LENGTH OR MASS FIT')
    print(f'Common host: {n} objects; {len(registry)} registry roots; {len(spatial)} S roots.')
    print('E# = Euclidean, minimum dimension #; N# = non-PSD, Gram rank #.')
    print('An N result is NOT a dimension-4 Euclidean realization.')
    print('\n readout           ' + ' '.join(f'{r["name"]:>7s}' for r in catalogue))
    for kind, _ in rules:
        print(f' {kind:18s}' + ' '.join(f'{cases[r["name"]]["readouts"][kind]["status"]:>7s}'
                                     for r in catalogue))
    print('\nReadout definitions (a squared distance is NOT squared again):')
    for kind, definition in rules:
        print(f'  {kind}: {definition}')
    print('\nAll 30 non-K5 five-point controls from the old eight-vertex pool:')
    for kind, _ in rules:
        print(f'  {kind:18s}: {null_counts[kind]}')
    print('\nK5-1 squared distances from the 81 S-overlap counts:')
    matrix = cases['K5-1']['readouts']['S_overlap']['D2']
    for row in matrix:
        print('  ' + ' '.join(f'{str(x):>7s}' for x in row))
    print('K4-1 parent-child hop metric has exact negative Gram eigenpair:')
    print('  B @ (1,-1,-1,1) = -(1,-1,-1,1).')
    print('All literal inherited-registry masks of the selected core vertices are empty.')
    print('PASS: 63 candidate and 270 nonclique metric audits; all exact rational arithmetic.')
    print('PASS: primitive exchange, vertex-order invariance, and local/common core path agreement.')
    print('Rank n-1 for full weighted ancestry features is generic, NOT a K5-specific defect.')
    print('Geodesic-distance failure is not itself physical strain; no strain law was added.')
    return dict(protocol='nine native-data distance conventions on unchanged seven cores',
                rules=dict(rules), common_host_parents=tuple(host.parents),
                catalogue=cases, controls=nulls, control_counts=null_counts,
                K5_four_point_subsets=facets, literal_registry_masks_collapsed=True,
                physical_mass_calibration=None)


dcu_mass_9 = _run_dcu_mass_9()
