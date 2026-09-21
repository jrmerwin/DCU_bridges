# Calibration J — G-sector structure, triangles, and placement in the SAME exact samples.
# Read-only: no new samples, partner choices, exhaustive rounds, service steps, or clocks.
from collections import Counter, deque
from itertools import combinations, product
import hashlib
import pickle


def rosetta_j_digest(*objects):
    class Sink:
        def __init__(self):
            self.h = hashlib.sha256()
        def write(self, data):
            self.h.update(data)
            return len(data)
    sink = Sink()
    pickle.Pickler(sink, protocol=5).dump(objects)
    return sink.h.digest()


def rosetta_g_structure_placement(exact, observations):
    """Audit G in three explicitly distinct graphs using already saved parent lists.

    1. Reference registry overlap: |Anc(r) intersect Anc(s)| > 3.
    2. Undirected parent-child graph inside each seven-vertex G object.
    3. Co-parent triangles reconstructed ONLY from children inside a queried ancestry.

    A decoded-cache union is never treated as a complete universe or local environment.
    The time axis below is exhaustive round / selected lineage link, NOT native tau.
    """
    if observations.get('scope') != 'exact uniform exhaustive-population observations, not native service histories':
        raise ValueError('Use the unchanged Calibration-I observation result.')
    before = rosetta_j_digest(exact, observations)
    reference, base = exact['reference'], exact['full5']
    rr = reference['rows']
    parents = tuple(None if p is None else tuple(p) for p in observations['parents'])
    masks = tuple(observations['registry_masks'])
    heights = tuple(observations['heights'])
    if len(rr) != 137 or len(parents) != len(masks) or len(parents) != len(heights):
        raise ValueError('Incomplete reference or observation arrays.')
    if parents[:len(base['parents'])] != tuple(base['parents']):
        raise ValueError('The Calibration-I cache does not retain the Calibration-H base.')
    if parents[:2] != (None, None) or heights[:2] != (0, 0) or masks[:2] != (0, 0):
        raise ValueError('Incorrect primitive convention.')
    root_by_i = {reference['index_by_object'][m]: z
                 for z, m in enumerate(base['mapped7'])
                 if m in reference['index_by_object']}
    if set(root_by_i) != set(range(137)):
        raise ValueError('All 137 literal registry roots must exist in the exact base.')
    index_by_root = {z: i for i, z in root_by_i.items()}
    seen_pairs = set()
    for z in range(2, len(parents)):
        p = parents[z]
        if p is None or len(p) != 2:
            raise ValueError('Incomplete composite parent list.')
        a, b = p
        if not 0 <= a < b < z or p in seen_pairs:
            raise ValueError('Parent pairs must be distinct, ordered, and topological.')
        seen_pairs.add(p)
        own = (1 << index_by_root[z]) if z in index_by_root else 0
        if masks[z] != (masks[a] | masks[b] | own):
            raise ValueError('Registry mask disagrees with actual ancestry inheritance.')
        if heights[z] != 1 + max(heights[a], heights[b]):
            raise ValueError('Exhaustive birth-round grading disagrees.')

    G = tuple(i for i, row in enumerate(rr) if row['sector'] == 'G')
    Gmask = sum(1 << i for i in G)
    sector_masks = reference['sector_masks']
    adjacency = tuple(frozenset(j for j, y in enumerate(rr) if j != i and
        len(reference['ancestors'][row['object_id']] & reference['ancestors'][y['object_id']]) > 3)
        for i, row in enumerate(rr))
    if adjacency != tuple(reference['adjacency']):
        raise ValueError('Reference overlap edges disagree with reference ancestry.')
    g_adj = {i: adjacency[i].intersection(G) for i in G}
    pending, components = set(G), []
    while pending:
        todo, component = [min(pending)], set()
        while todo:
            i = todo.pop()
            if i in component:
                continue
            component.add(i)
            todo.extend(g_adj[i] - component)
        pending -= component
        components.append(tuple(sorted(component)))
    g_triangles = tuple((a, b, c) for a, b, c in combinations(G, 3)
                        if b in g_adj[a] and c in g_adj[a] and c in g_adj[b])
    component_by_i = {i: k for k, cc in enumerate(components) for i in cc}

    # Independently expand a concise candidate description and compare its exact terms.
    # This is an anatomy check, not a new generator for the observation cache.
    rules = {}
    def ref_pair(a, b):
        return reference['pairs'][tuple(sorted((a, b)))]
    x0 = ref_pair(0, 1)
    for alpha, beta, gamma in product((0, 1), repeat=3):
        x1 = ref_pair(x0, alpha)
        x2 = ref_pair(x1, beta)
        x3 = ref_pair(x2, gamma)
        for name, target, terminal_cycle in (('near', x2, 3), ('skip', x1, 4)):
            obj = ref_pair(x3, target)
            i = reference['index_by_object'].get(obj)
            if i is None or i in rules:
                raise AssertionError('The explicit G pattern descriptions are not unique.')
            rules[i] = {'primitive_choices': (alpha, beta, gamma), 'closure': name,
                        'terminal_undirected_cycle': terminal_cycle,
                        'chain_reference_objects': (x0, x1, x2, x3),
                        'target_reference_object': target}
    if set(rules) != set(G):
        raise AssertionError('Explicit chain-and-shortcut patterns do not equal the G sector.')

    profile_cache = {}
    def node_profile(root):
        if root in profile_cache:
            return profile_cache[root]
        if not 0 <= root < len(parents):
            raise ValueError('Queried root is absent from the saved cache.')
        distance, queue = {root: 0}, deque([root])
        while queue:
            z = queue.popleft()
            for p in parents[z] or ():
                if p not in distance:
                    distance[p] = distance[z] + 1
                    queue.append(p)
        observed_mask = sum(1 << i for i, g in root_by_i.items() if g in distance)
        if observed_mask != masks[root]:
            raise AssertionError('Independent ancestral traversal disagrees with registry mask.')
        present_g = tuple(i for i in G if masks[root] & (1 << i))
        overlap_tris = tuple(t for t in g_triangles if all(masks[root] & (1 << i) for i in t))
        # ONLY actual children contained in this root's own ancestry provide edges.
        cp_adj = {z: set() for z in distance}
        edge_witness = {}
        pc_adj = {z: set() for z in distance}
        for z in distance:
            if parents[z] is not None:
                a, b = parents[z]
                cp_adj[a].add(b); cp_adj[b].add(a)
                edge_witness[(a, b)] = z
                pc_adj[z].update((a, b)); pc_adj[a].add(z); pc_adj[b].add(z)
        cp_tris = tuple((a, b, c) for a in sorted(cp_adj) if a >= 2
                        for b in sorted(v for v in cp_adj[a] if v > a)
                        for c in sorted(v for v in cp_adj[a] & cp_adj[b] if v > b))
        cp_g, cp_common = [], []
        for t in cp_tris:
            a, b, c = t
            union = (masks[a] | masks[b] | masks[c]) & Gmask
            common = masks[a] & masks[b] & masks[c] & Gmask
            if union:
                record = {'core': t, 'edge_children': tuple(edge_witness[p] for p in combinations(t, 2)),
                          'G_indices': tuple(i for i in G if union & (1 << i)),
                          'common_G_indices': tuple(i for i in G if common & (1 << i))}
                cp_g.append(record)
                if common:
                    cp_common.append(record)
        terminal_triangle = False
        terminal_path = None
        if parents[root] is not None:
            a, b = parents[root]
            # Shortest undirected old-ancestry path between the two parents,
            # excluding the apex. Its two new edges form the terminal cycle.
            pd, pq = {a: 0}, deque([a])
            while pq and b not in pd:
                z = pq.popleft()
                for v in pc_adj[z]:
                    if v != root and v not in pd:
                        pd[v] = pd[z]+1; pq.append(v)
            terminal_path = pd.get(b)
            terminal_triangle = terminal_path == 1
        out = {'root': root, 'birth_round': heights[root], 'ancestry_size': len(distance),
               'registry_mask': masks[root], 'G_indices': present_g,
               'G_counts_by_component': tuple(sum(i in cc for i in present_g) for cc in components),
               'G_depths': {i: distance[root_by_i[i]] for i in present_g},
               'G_overlap_triangles': overlap_tris,
               'nonprimitive_coparent_triangles': cp_tris,
               'coparent_triangles_inheriting_G': tuple(cp_g),
               'coparent_triangles_sharing_G': tuple(cp_common),
               'apex_parent_child_triangle': terminal_triangle,
               'apex_shortest_undirected_cycle': None if terminal_path is None else terminal_path+2,
               'distance': distance}
        profile_cache[root] = out
        return out

    anatomy = []
    for i in G:
        root = root_by_i[i]; row = rr[i]; info = node_profile(root)
        if info['ancestry_size'] != 7 or info['apex_shortest_undirected_cycle'] != rules[i]['terminal_undirected_cycle']:
            raise AssertionError('Native G support or terminal cycle disagrees with pattern expansion.')
        anatomy.append({'index': i, 'label': row['label'], 'root': root,
            'component': component_by_i[i], 'support': tuple(sorted(info['distance'])),
            'parents_on_support': tuple((z, parents[z]) for z in sorted(info['distance'])),
            'shape': row['shape'], 'parent_overlap': row['parent_overlap'],
            'cross_sector_degree': tuple(sum(rr[j]['sector'] == s for j in adjacency[i]) for s in ('S','I','G')),
            'G_triangle_memberships': sum(i in t for t in g_triangles), **rules[i]})

    print('G STRUCTURE AND PLACEMENT — same exact-closure observations, no new sampling')
    print('Axes are exhaustive round and lineage link, NOT maintenance ticks.')
    print('REFERENCE G-OVERLAP GRAPH: edge iff original registry ancestries overlap in >3 vertices')
    print(f'  vertices={len(G)}; edges={sum(map(len, g_adj.values()))//2}; '
          f'components={tuple(len(c) for c in components)}; pure-G K3={len(g_triangles)}')
    for k, cc in enumerate(components):
        print(f'  component {k}: labels={tuple(rr[i]["label"] for i in cc)}; '
              f'complete={all(len(g_adj[i] & set(cc)) == len(cc)-1 for i in cc)}')
    print('  per-entry pure-G triangle counts:', dict(Counter(a['G_triangle_memberships'] for a in anatomy)))
    print('ACTUAL G-SEVEN ANATOMY: x0={a,b}; x1={x0,alpha}; x2={x1,beta}; x3={x2,gamma}')
    print('  alpha,beta,gamma are each a or b; near={x3,x2}; skip={x3,x1}.')
    print(f'{"label":>6} {"root":>6} {"a/b choices":>12} {"closure":>8} {"component":>10} {"terminal cycle":>15}')
    for a in anatomy:
        bits = ''.join('ab'[v] for v in a['primitive_choices'])
        print(f'{a["label"]:>6} {a["root"]:6d} {bits:>12} {a["closure"]:>8} '
              f'{a["component"]:10d} {a["terminal_undirected_cycle"]:15d}')

    sample_rows, sample_summary = [], []
    print('\nSAME UNIFORM SAMPLES — counts of draws, NOT disjoint particles')
    print(f'{"round":>6} {"draws":>6} {"any G":>7} {"G-overlap K3":>13} {"coparent K3":>12} {"K3 core has G":>14} {"K3 common G":>12}')
    for s in observations['sample_rounds']:
        selected = [r for r in observations['samples'] if r['round'] == s]
        if not selected:
            raise ValueError('A declared sample round has no saved draws.')
        measured = []
        for draw, r in enumerate(selected):
            p = node_profile(r['root'])
            if r['registry_mask'] != p['registry_mask'] or r['ancestry_size'] != p['ancestry_size']:
                raise AssertionError('Saved Calibration-I sample disagrees with its parents.')
            measured.append(p)
            sample_rows.append({'round': s, 'draw': draw, 'root': r['root'], 'rank': r['rank']})
        row = {'round': s, 'draws': len(selected),
               'with_G': sum(bool(p['G_indices']) for p in measured),
               'with_G_overlap_K3': sum(bool(p['G_overlap_triangles']) for p in measured),
               'with_nonprimitive_coparent_K3': sum(bool(p['nonprimitive_coparent_triangles']) for p in measured),
               'with_coparent_K3_inheriting_G': sum(bool(p['coparent_triangles_inheriting_G']) for p in measured),
               'with_coparent_K3_sharing_G': sum(bool(p['coparent_triangles_sharing_G']) for p in measured)}
        sample_summary.append(row)
        print(f'{s:6d} {row["draws"]:6d} {row["with_G"]:7d} {row["with_G_overlap_K3"]:13d} '
              f'{row["with_nonprimitive_coparent_K3"]:12d} {row["with_coparent_K3_inheriting_G"]:14d} '
              f'{row["with_coparent_K3_sharing_G"]:12d}')
    print('ALL positive G-overlap-K3 draws (root is a cache ID, not an old simulation ID):')
    for r in sample_rows:
        p = profile_cache[r['root']]
        if p['G_overlap_triangles']:
            labels = tuple(tuple(rr[i]['label'] for i in t) for t in p['G_overlap_triangles'])
            depths = tuple((rr[i]['label'], d) for i, d in p['G_depths'].items())
            print(f'  U{r["round"]} draw={r["draw"]} root={r["root"]}; ancestry={p["ancestry_size"]}; '
                  f'triangles={labels}; depths={depths}')

    lineage_rows = []
    print('\nSAME 16 G-SEEDED LINEAGES — initial G is held by actual identity')
    print(f'{"label":>6} {"initial-root depths":>21} {"G counts":>21} {"final G component counts":>25} {"shortcuts":>10}')
    transition_count = 0
    for path in observations['lineages']:
        i = path['reference_index']
        nodes = tuple(path['roots'])
        if nodes[0] != root_by_i[i] or len(path['coparents']) != len(nodes)-1:
            raise ValueError('A saved lineage has inconsistent starting identity or partners.')
        pp = [node_profile(z) for z in nodes]
        depths = tuple(p['distance'][root_by_i[i]] for p in pp)
        gains = []
        for k, q in enumerate(path['coparents']):
            p, z = nodes[k:k+2]
            if tuple(sorted((p, q['root']))) != parents[z] or heights[z] != q['round']+1:
                raise AssertionError('A saved lineage edge or exhaustive birth round changed.')
            if masks[p] & ~masks[z] or not set(pp[k]['G_overlap_triangles']) <= set(pp[k+1]['G_overlap_triangles']):
                raise AssertionError('Old inherited registry entries or fixed overlap edges were lost.')
            gain = depths[k]+1-depths[k+1]
            if gain < 0:
                raise AssertionError('Old parent route should remain available.')
            gains.append(gain)
            transition_count += 1
        row = {'reference_index': i, 'label': path['label'], 'sector': path['sector'],
               'roots': nodes, 'observation_rounds': (5,) + tuple(q['round']+1 for q in path['coparents']),
               'initial_G_or_other_depths': depths, 'depth_gains': tuple(gains),
               'G_counts': tuple(len(p['G_indices']) for p in pp),
               'G_component_counts': tuple(p['G_counts_by_component'] for p in pp),
               'G_overlap_K3_counts': tuple(len(p['G_overlap_triangles']) for p in pp),
               'G_coparent_K3_counts': tuple(len(p['coparent_triangles_inheriting_G']) for p in pp)}
        lineage_rows.append(row)
        if path['sector'] == 'G':
            print(f'{path["label"]:>6} {str(depths):>21} {str(row["G_counts"]):>21} '
                  f'{str(row["G_component_counts"][-1]):>25} {sum(g>0 for g in gains):10d}')
    g_paths = [p for p in lineage_rows if p['sector'] == 'G']
    if len(g_paths) != 16 or len(lineage_rows) != 137:
        raise ValueError('Keep all 137 original Calibration-I lineages, including all 16 G starts.')
    print(f'  All {transition_count} saved lineage links verified. No new link was chosen.')
    print('  G-seeded lineages with a later pure-G overlap K3:', sum(any(p['G_overlap_K3_counts']) for p in g_paths))
    print('  G-seeded lineages with a later internally witnessed co-parent K3 inheriting G:',
          sum(any(p['G_coparent_K3_counts']) for p in g_paths))
    print('\nGraph distinctions: overlap K3 uses registry ancestry intersections; co-parent K3 requires three witnessed comparisons.')
    print('The co-parent table reconstructs each root ancestry independently; no other cache objects supply edges.')
    print('In ambient complete U_s, all U_(s-1) objects are mutually co-parent connected; that generic completeness is NOT a particle detector.')
    print('The near/skip terminal cycles are in the UNDIRECTED parent-child graph. The directed ancestry remains acyclic.')
    print('Unchanged old G supports/edges are inherited identities, not evidence of active stabilization or a protected local position.')
    print('Distances are shortest parent steps from the changing representative; they are not spatial distances or local clock readings.')
    print('All signs, transfer dynamics, 136/27/half factors, energy and physical mass remain unassigned.')
    if rosetta_j_digest(exact, observations) != before:
        raise AssertionError('An input reference or Calibration-I observation was changed.')
    return {'scope': 'read-only G structure and placement on the existing exact-closure observations',
            'G_anatomy': anatomy, 'G_components': tuple(components),
            'G_overlap_triangles': g_triangles, 'sample_summary': sample_summary,
            'sample_rows': sample_rows, 'lineages': lineage_rows,
            'profiles': profile_cache, 'transition_count': transition_count}


if 'rosetta_constructor_reference' not in globals() or 'rosetta_exact_closure_samples' not in globals():
    raise RuntimeError('Run Calibration H1 and I first, keeping their saved results in this kernel.')
rosetta_G_structure_placement = rosetta_g_structure_placement(
    rosetta_constructor_reference, rosetta_exact_closure_samples)
