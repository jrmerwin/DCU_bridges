# Calibration F — exact historical DAG-7 registry content of the EXISTING K3 candidates.
# The 137 catalogue objects are a REFERENCE DICTIONARY, never inserted into a sampled run.
# No mass formula, force, registry coupling, new geometry, workload change, or RNG draws.
from collections import Counter
from itertools import combinations
import hashlib
import pickle


def rosetta_dag7_reference():
    """Enumerate possible two-primitive unordered objects with inclusive DAG size <=7.

    This is a finite size-bounded dictionary, NOT a claim that the full engine halts.
    Parentage keeps primitive identities; uncoloured unfolded tree shapes are used
    ONLY for the published catalogue-multiplicity partition, not identity matching.
    """
    parents = [None, None]
    ancestors = [frozenset({0}), frozenset({1})]
    shapes, terms, pairs = ['P', 'P'], ['a', 'b'], {}
    while True:
        batch = []
        for a, b in combinations(range(len(parents)), 2):
            if (a, b) in pairs:
                continue
            union = ancestors[a] | ancestors[b]
            if len(union) + 1 <= 7:
                batch.append((a, b, union))
        if not batch:
            break
        for a, b, union in batch:
            z = len(parents)
            parents.append((a, b))
            ancestors.append(union | {z})
            pairs[a, b] = z
            shapes.append('(' + '|'.join(sorted((shapes[a], shapes[b]))) + ')')
            terms.append('(' + '|'.join(sorted((terms[a], terms[b]))) + ')')

    ids = sorted((z for z in range(len(parents)) if len(ancestors[z]) == 7),
                 key=lambda z: terms[z])
    multiplicities = Counter(shapes[z] for z in ids)  # COMPLETE catalogue, not sampled frequency.
    rows = []
    for i, z in enumerate(ids):
        a, b = parents[z]
        overlap = len(ancestors[a] & ancestors[b])
        multiplicity = multiplicities[shapes[z]]
        sector = 'S' if overlap <= 3 else ('G' if multiplicity == 8 else 'I')
        rows.append({'index': i, 'label': f'R{i:03d}', 'object_id': z,
                     'term': terms[z], 'shape': shapes[z], 'parent_overlap': overlap,
                     'shape_multiplicity': multiplicity, 'sector': sector})
    sectors = {s: sum(1 << r['index'] for r in rows if r['sector'] == s)
               for s in ('S', 'I', 'G')}
    adjacency = tuple(frozenset(j for j, y in enumerate(ids)
        if y != z and len(ancestors[z] & ancestors[y]) > 3) for z in ids)
    # Source values are reference checks, never targets used to choose a sampled motif.
    if (len(parents) != 173 or len(rows) != 137 or
        Counter(r['sector'] for r in rows) != Counter({'S': 81, 'I': 40, 'G': 16}) or
        sum(map(len, adjacency)) // 2 != 5347 or
        Counter(map(len, adjacency)) != Counter({73: 126, 136: 11})):
        raise ValueError('The size-bounded dictionary does not reproduce the historical registry.')
    for r in rows:
        if r['sector'] == 'G':
            cross = tuple(sum(rows[j]['sector'] == s for j in adjacency[r['index']])
                          for s in ('S', 'I', 'G'))
            if cross != (45, 21, 7):
                raise ValueError('Historical G-sector reference connectivity disagrees.')
    single = tuple(r['index'] for r in rows if r['shape_multiplicity'] == 1)
    if len(single) != 1:
        raise ValueError('Expected one multiplicity-one catalogue entry.')
    return {'parents': tuple(parents), 'ancestors': tuple(ancestors),
            'pairs': pairs, 'terms': tuple(terms), 'rows': tuple(rows),
            'index_by_object': {z: i for i, z in enumerate(ids)},
            'sector_masks': sectors, 'adjacency': adjacency,
            'multiplicity_one_index': single[0]}


def rosetta_registry_projection(structure, end, reference):
    """Literal inherited registry membership; no local re-rooting or forced completion.

    Match actual small objects by their COMPLETE recursive parentage, then propagate
    their registry bits through the existing DAG. Each factor counts once per support.
    Native IDs 0/1 correspond to reference atoms a/b; a global atom swap just relabels
    symmetric catalogue entries. Large native objects are not coerced into registry rows.
    """
    if type(end) is not int or not 2 <= end <= len(structure):
        raise ValueError('Use a stored population prefix.')
    mapped, masks, occurrences, pairs = [0, 1], [0, 0], {}, {}
    if structure.parents[:2] != [None, None] and tuple(structure.parents[:2]) != (None, None):
        raise ValueError('Expected two original primitive roots.')
    for z in range(2, end):
        a, b = sorted(structure.parents[z])
        if not 0 <= a < b < z or (a, b) in pairs:
            raise ValueError('Invalid or repeated native parent pair.')
        union = structure.ancestors[a] | structure.ancestors[b]
        if structure.ancestors[z] != union | {z}:
            raise ValueError('The saved inclusive ancestry disagrees with parentage.')
        pairs[a, b] = z
        native_size = len(structure.ancestors[z])
        small = None
        if native_size <= 7:
            if mapped[a] is None or mapped[b] is None:
                raise ValueError('A small native object has an incompatible parent.')
            small = reference['pairs'].get(tuple(sorted((mapped[a], mapped[b]))))
            if small is None or len(reference['ancestors'][small]) != native_size:
                raise ValueError('A native small object has no exact dictionary match.')
        mapped.append(small)
        mask = masks[a] | masks[b]
        if native_size == 7:
            i = reference['index_by_object'][small]
            if i in occurrences:
                raise ValueError('Duplicate actual realization of the same native term.')
            occurrences[i] = z
            mask |= 1 << i
        masks.append(mask)
    # Independent set-membership check of propagated registry incidence.
    for z in range(end):
        literal = sum(1 << i for i, native in occurrences.items()
                      if native in structure.ancestors[z])
        if literal != masks[z]:
            raise ValueError('Registry incidence propagation failed its ancestry check.')
    return {'masks': tuple(masks), 'occurrences': occurrences,
            'mapped_small_objects': tuple(mapped), 'pairs': pairs}


def rosetta_registry_anatomy(growth, candidates, audit_population=3072):
    """Inspect all original Calibration-A triads and their exact inherited registry.

    No numerical workload selects the candidates. M is the three core objects.
    Report their separate/shared registry references and any extra DAG-7 entries
    introduced by the three edge-witness children; also retain the same frozen
    external co-parent-neighbour reference used in Calibration E.
    """
    discovery = candidates['discovery_population']
    if type(audit_population) is not int or not 2 < discovery <= audit_population:
        raise ValueError('Use the original discovery and an available later prefix.')
    if candidates.get('candidate_rule') != 'nonprimitive co-parent K3':
        raise ValueError('Use the unchanged original Calibration-A candidate catalogue.')
    branches = {b['seed']: b for b in growth['branches']}
    if (len(branches) != len(growth['branches']) or
        len({b['seed'] for b in candidates['branches']}) != len(candidates['branches']) or
        set(branches) != {b['seed'] for b in candidates['branches']}):
        raise ValueError('The distinct seed panels must match.')

    class Sink:
        def __init__(self): self.h = hashlib.sha256()
        def write(self, data): self.h.update(data); return len(data)
    def fingerprint():
        sink = Sink()
        pickle.Pickler(sink, protocol=5).dump((growth, candidates))
        return sink.h.digest()
    before = fingerprint()
    reference = rosetta_dag7_reference()
    ref_rows, sector_masks = reference['rows'], reference['sector_masks']
    singleton_bit = 1 << reference['multiplicity_one_index']
    def profile(mask):
        return tuple((mask & sector_masks[s]).bit_count() for s in ('S', 'I', 'G'))
    def indices(mask):
        return tuple(i for i in range(len(ref_rows)) if mask & (1 << i))
    def union_masks(values):
        result = 0
        for value in values: result |= value
        return result

    print('REGISTRY ANATOMY — reference dictionary vs actual sampled membership')
    print('Dictionary: 137 DAG-7 objects; S/I/G=81/40/16; overlap edges=5347.')
    print('Each reference G entry has cross-sector degree (45,21,7). No registry objects enter the simulation.')
    print(f'Original candidates/backgrounds at N={discovery}; sampled coverage ends at N={audit_population}.')
    output = []
    for saved in candidates['branches']:
        branch = branches[saved['seed']]
        run, seed = branch['run'], branch['seed']
        g, n = run.structure, len(run.structure)
        cp = branch['checkpoints'][-1]
        if (n < audit_population or branch['status'] not in ('complete','time_budget','update_budget') or
            (cp['N'], cp['updates'], cp['ticks'], cp['F']) != (n,run.updates,run.ticks,run.work)):
            raise ValueError('Use a verified snapshot beyond the common audit endpoint.')
        if (run.parameters.gamma,run.parameters.persistent_work,run.parameters.relief_cap) != (2,0,0):
            raise ValueError('Keep the unchanged Gamma=2, m=H=0 histories.')
        projection = rosetta_registry_projection(g, audit_population, reference)
        masks, occurrences = projection['masks'], projection['occurrences']
        adjacent = [set() for _ in range(discovery)]
        edge_children = {}
        for z in range(2, discovery):
            a,b = sorted(g.parents[z]); adjacent[a].add(b); adjacent[b].add(a)
            edge_children[a,b] = z
        original_cores = {(a,b,c) for a in range(2,discovery)
            for b in adjacent[a] if b > a
            for c in adjacent[a] & adjacent[b] if c > b}
        if original_cores != {tuple(r['core']) for r in saved['candidates']}:
            raise ValueError('Original candidate membership changed.')
        coverage = []
        print(f'\nseed={seed}: current N={n}; reference membership at two prefixes')
        for N in (discovery, audit_population):
            present = sum(1 << i for i,z in occurrences.items() if z < N)
            coverage.append({'N':N,'registry_indices':indices(present),'counts_SIG':profile(present)})
            print(f'  N={N}: actual registry entries={present.bit_count()}/137; '
                  f'S/I/G={profile(present)}; m=1 entry present={bool(present & singleton_bit)}')
        print('ACTUAL DAG-7 OBJECTS: native ID -> reference label / sector')
        print('  '+(', '.join(f"{z}->{ref_rows[i]['label']}/{ref_rows[i]['sector']}"
                            for i,z in sorted(occurrences.items(),key=lambda x:x[1])) or 'none'))
        print(f"{'core':>21} {'R per core':>12} {'core S/I/G':>12} {'common G':>9} "
              f"{'edge+ S/I/G':>12} {'ambient G':>10} {'m=1':>5}")
        rows=[]
        for old in saved['candidates']:
            core=tuple(old['core']); cm=tuple(masks[z] for z in core)
            ec=tuple(edge_children[a,b] for a,b in combinations(core,2))
            boundary=tuple(sorted((a,b) for a in core for b in adjacent[a] if b not in core))
            if ec != tuple(old['edge_children']) or boundary != tuple(old['external_coparent_edges']):
                raise ValueError('Candidate edge witnesses or frozen boundary changed.')
            pooled=union_masks(cm); common=cm[0]&cm[1]&cm[2]
            support=pooled|union_masks(masks[z] for z in ec)
            extra=support & ~pooled
            ambient=tuple(sorted({b for _,b in boundary}))
            ambient_mask=union_masks(masks[z] for z in ambient)
            incidence={i:tuple(core[j] for j in range(3) if cm[j] & (1<<i)) for i in indices(pooled)}
            pattern_counts={s:dict(Counter(sum(1<<j for j in range(3) if cm[j]&(1<<i))
                for i in indices(pooled & sector_masks[s]))) for s in ('S','I','G')}
            row={'core':core,'edge_children':ec,'per_core_registry':tuple(indices(x) for x in cm),
                'per_core_counts_SIG':tuple(profile(x) for x in cm),
                'union_registry':indices(pooled),'union_counts_SIG':profile(pooled),
                'all_three_registry':indices(common),'all_three_counts_SIG':profile(common),
                'edge_added_registry':indices(extra),'edge_added_counts_SIG':profile(extra),
                'ambient_nodes':ambient,'ambient_registry':indices(ambient_mask),
                'ambient_counts_SIG':profile(ambient_mask),'incidence':incidence,
                'incidence_pattern_counts_SIG':pattern_counts,
                'multiplicity_one_in_core_support':bool(pooled & singleton_bit)}
            rows.append(row)
            print(f"{str(core):>21} {str(tuple(x.bit_count() for x in cm)):>12} "
                  f"{str(profile(pooled)):>12} {profile(common)[2]:9d} "
                  f"{str(profile(extra)):>12} {profile(ambient_mask)[2]:10d} "
                  f"{str(bool(pooled & singleton_bit)):>5}")
        output.append({'seed':seed,'audit_N':audit_population,'coverage':coverage,
                       'registry_occurrences':dict(occurrences),'candidates':rows,
                       'object_registry_masks':masks})
    if fingerprint()!=before:
        raise AssertionError('A source history, candidate catalogue or RNG changed.')
    print('\nR per core counts literal inherited DAG-7 entries; core S/I/G counts their UNION once.')
    print('common G counts identical G entries inherited by ALL THREE cores, not active coupling.')
    print('edge+ lists new entries contributed by the actual triangle-closing children, outside core ancestry.')
    print('ambient G counts the UNION in the same distinct frozen external neighbours; not a vacuum or enrichment test.')
    print('m=1 is the unique catalogue multiplicity-one entry (in S), NOT the full 40-entry I sector.')
    print('Its correspondence to the RMR single interface bit remains a hypothesis, not a detector rule.')
    print('Catalogue multiplicity uses all 137 reference objects, never occurrence counts in the sampled universe.')
    print('Fixed core ancestry is immutable: fixed registry content is NOT independent evidence of mass stability.')
    print('This tests literal inherited registry content. No re-rooted copy, added registry, 27-factor, or 1/2 dynamics is assumed.')
    print('All candidates and bit-incidence patterns retained. No mass computed; all original states/RNGs unchanged.')
    return {'reference':reference,'branches':output,'discovery_population':discovery,
            'audit_population':audit_population,'observable':'literal DAG-7 registry incidence'}


if not all(k in globals() for k in ('rosetta_calibration_growth_4096','rosetta_calibration_candidates')):
    raise RuntimeError('Keep the current Calibration-C continuation and ORIGINAL Calibration-A catalogue; no new growth needed.')
rosetta_registry_anatomy_results = rosetta_registry_anatomy(
    rosetta_calibration_growth_4096, rosetta_calibration_candidates)
