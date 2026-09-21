# Calibration G — actual local two-source registry templates, not an inserted registry.
# NEW readout hypothesis: two existing objects may serve as opaque boundary atoms.
# The native graph, ancestry, charges, service and RNG are unchanged.
from collections import Counter
import hashlib
import pickle


def rosetta_local_registry_embeddings(parents, reference):
    """Find all realized size-seven constructions relative to two boundary objects.

    Every nonempty two-source construction contains the first join of the sources.
    Thus all possible source pairs are among the ACTUALLY realized parent pairs.
    For each such pair, follow only actual constructions, using the same <=7
    dictionary as Calibration F. Never complete missing constructions.
    Size seven counts two boundary objects plus five distinct internal composites;
    ancestry inside either boundary is NOT part of this local template.
    """
    parents = tuple(None if p is None else tuple(sorted(p)) for p in parents)
    if len(parents) < 3 or parents[:2] != (None, None):
        raise ValueError('Expected a two-primitive population prefix with a first birth.')
    actual_pairs = {}
    ancestry_bits = [1, 2]
    for z, p in enumerate(parents[2:], 2):
        if p is None or len(p) != 2 or not 0 <= p[0] < p[1] < z or p in actual_pairs:
            raise ValueError('Invalid, non-topological or repeated parent pair.')
        actual_pairs[p] = z
        ancestry_bits.append(ancestry_bits[p[0]] | ancestry_bits[p[1]] | (1 << z))

    occurrences, anchor_rows = [], []
    ref_pairs = reference['pairs']
    index_by_object = reference['index_by_object']
    for anchors, first_join in sorted(actual_pairs.items(), key=lambda kv: kv[1]):
        # These maps are observer bookkeeping only. They never modify the graph.
        native_to_ref = {anchors[0]: 0, anchors[1]: 1}
        ref_to_native = {0: anchors[0], 1: anchors[1]}
        known = list(anchors)
        pos = 1
        while pos < len(known):
            x = known[pos]
            # Each pair of discovered objects is queried once.
            for y in known[:pos]:
                child_ref = ref_pairs.get(tuple(sorted((native_to_ref[x], native_to_ref[y]))))
                if child_ref is None:
                    continue  # The resulting relative support would exceed seven.
                child = actual_pairs.get(tuple(sorted((x, y))))
                if child is None:
                    continue  # An unexecuted legal construction is NOT inserted.
                if child in native_to_ref:
                    if native_to_ref[child] != child_ref:
                        raise AssertionError('Local substitution produced conflicting identities.')
                    continue
                if child_ref in ref_to_native:
                    raise AssertionError('Two actual objects realize one local expression.')
                native_to_ref[child] = child_ref
                ref_to_native[child_ref] = child
                known.append(child)
            pos += 1

        nested = bool(ancestry_bits[anchors[1]] & (1 << anchors[0]))
        hits = []
        for child, child_ref in native_to_ref.items():
            if child_ref not in index_by_object:
                continue
            i = index_by_object[child_ref]
            local_support = tuple(sorted(ref_to_native[r]
                                         for r in reference['ancestors'][child_ref]))
            if len(local_support) != 7 or not set(anchors) <= set(local_support):
                raise AssertionError('A relative registry support must have exactly seven vertices.')
            # Verify every internal edge. Boundaries are deliberately opaque.
            for r in reference['ancestors'][child_ref]:
                if r < 2:
                    continue
                expected = tuple(sorted(ref_to_native[t] for t in reference['parents'][r]))
                if parents[ref_to_native[r]] != expected:
                    raise AssertionError('A local template edge is not an actual parent edge.')
            sector = reference['rows'][i]['sector']
            hits.append(i)
            occurrences.append({
                'anchors': anchors, 'root': child, 'reference_index': i,
                'reference_label': reference['rows'][i]['label'], 'sector': sector,
                'local_support': local_support, 'nested_anchors': nested,
                'primitive_anchors': sum(a < 2 for a in anchors),
                'boundary_shared_ancestors': (ancestry_bits[anchors[0]] &
                                              ancestry_bits[anchors[1]]).bit_count(),
                'root_full_ancestry_size': ancestry_bits[child].bit_count(),
            })
        anchor_rows.append({'anchors': anchors, 'first_join': first_join,
                            'bounded_realized_objects': len(native_to_ref),
                            'registry_indices': tuple(sorted(hits)),
                            'nested_anchors': nested})
    occurrences.sort(key=lambda r: (r['root'], r['anchors'], r['reference_index']))
    return {'occurrences': occurrences, 'anchor_rows': anchor_rows,
            'ancestry_bits': tuple(ancestry_bits), 'parents': parents}


def rosetta_local_registry_scan(growth, literal_results, display_examples=8):
    """Same frozen candidates; new, explicit boundary-relative registry interpretation."""
    if literal_results.get('observable') != 'literal DAG-7 registry incidence':
        raise ValueError('Use the original Calibration-F result.')
    end = literal_results['audit_population']
    discovery = literal_results['discovery_population']
    reference = literal_results['reference']
    branches = {b['seed']: b for b in growth['branches']}
    seeds = [b['seed'] for b in literal_results['branches']]
    if len(branches) != len(growth['branches']) or len(set(seeds)) != len(seeds) or set(seeds) != set(branches):
        raise ValueError('The same distinct seed panel is required.')
    if len(reference['rows']) != 137 or len(reference['parents']) != 173:
        raise ValueError('Use the unchanged Calibration-F reference catalogue.')
    if type(display_examples) is not int or display_examples < 0:
        raise ValueError('display_examples must be a nonnegative integer.')

    class Sink:
        def __init__(self): self.hash = hashlib.sha256()
        def write(self, data): self.hash.update(data); return len(data)
    def fingerprint():
        sink = Sink()
        pickle.Pickler(sink, protocol=5).dump((growth, literal_results))
        return sink.hash.digest()
    before = fingerprint()
    outputs = []
    print('LOCAL REGISTRY TEST — changed readout, unchanged native construction')
    print('Two ACTUAL objects are boundary atoms; only ACTUAL later pair constructions are followed.')
    print('Seven = two boundary objects + five internal composites. Full ancestry is never erased.')
    print(f'Original triad selection N={discovery}; common audit prefix N={end}.')

    for saved in literal_results['branches']:
        seed = saved['seed']; branch = branches[seed]; run = branch['run']; g = run.structure
        if len(g) < end or (run.parameters.gamma, run.parameters.persistent_work, run.parameters.relief_cap) != (2, 0, 0):
            raise ValueError('Keep the original Gamma=2, m=H=0 histories past the audit prefix.')
        scan = rosetta_local_registry_embeddings(g.parents[:end], reference)
        occurrences, ancestry = scan['occurrences'], scan['ancestry_bits']
        # Verify the new generalization reduces exactly to F at the original atoms.
        original = {r['reference_index']: r['root'] for r in occurrences if r['anchors'] == (0, 1)}
        if original != saved['registry_occurrences']:
            raise AssertionError('Original-primitive slice disagrees with Calibration F.')
        # Independently compare all reconstructed ancestry to the saved native sets.
        for z in range(end):
            actual = sum(1 << a for a in g.ancestors[z])
            if actual != ancestry[z]:
                raise AssertionError('Native ancestry and parent-list reconstruction disagree.')

        print(f'\nseed={seed}: local two-source catalogue occurrences')
        print(f"{'N':>6} {'source pairs':>12} {'with R7':>8} "
              f"{'R occurrences':>14} {'G occurrences':>14} {'maxR/pair':>10}")
        coverage = []
        for N in (discovery, end):
            occ = [r for r in occurrences if r['root'] < N]
            per_pair = Counter(r['anchors'] for r in occ)
            type_ids = {r['reference_index'] for r in occ}
            g_types = {r['reference_index'] for r in occ if r['sector'] == 'G'}
            row = {'N': N, 'source_pairs': N-2, 'pairs_with_registry': len(per_pair),
                   'template_types': len(type_ids), 'G_template_types': len(g_types),
                   'occurrences': len(occ), 'G_occurrences': sum(r['sector'] == 'G' for r in occ),
                   'max_registry_per_pair': max(per_pair.values(), default=0),
                   'complete_source_pairs': sum(v == 137 for v in per_pair.values()),
                   'G_by_primitive_anchor_count': dict(Counter(r['primitive_anchors'] for r in occ if r['sector'] == 'G')),
                   'G_nested_anchors': sum(r['sector'] == 'G' and r['nested_anchors'] for r in occ)}
            coverage.append(row)
            print(f"{N:6d} {N-2:12d} {len(per_pair):8d} "
                  f"{len(occ):14d} {row['G_occurrences']:14d} {row['max_registry_per_pair']:10d}")
        print('  All 137 types at ONE source pair:', coverage[-1]['complete_source_pairs'])
        print('  G occurrences by number of ORIGINAL-primitive boundary objects (0/1/2):',
              {i: coverage[-1]['G_by_primitive_anchor_count'].get(i,0) for i in (0,1,2)})
        print('  G occurrences with ancestor-related boundary objects:', coverage[-1]['G_nested_anchors'])

        print('SAME FROZEN TRIADS — counts are anchored occurrences, not independent particles')
        print(f"{'core':>21} {'R occ':>7} {'G occ':>7} {'common G':>9} "
              f"{'edge+ G':>8} {'ambient G':>10}")
        candidate_rows = []
        for old in saved['candidates']:
            core, edges = tuple(old['core']), tuple(old['edge_children'])
            cm = [ancestry[z] for z in core]
            core_union = cm[0] | cm[1] | cm[2]
            core_common = cm[0] & cm[1] & cm[2]
            edge_union = 0
            for z in edges: edge_union |= ancestry[z]
            ambient_union = 0
            for z in old['ambient_nodes']: ambient_union |= ancestry[z]
            present = tuple(j for j,r in enumerate(occurrences) if core_union & (1<<r['root']))
            g_present = tuple(j for j in present if occurrences[j]['sector']=='G')
            g_common = tuple(j for j in g_present if core_common & (1<<occurrences[j]['root']))
            g_edge = tuple(j for j,r in enumerate(occurrences) if r['sector']=='G' and
                           edge_union & (1<<r['root']) and not core_union & (1<<r['root']))
            g_ambient = tuple(j for j,r in enumerate(occurrences) if r['sector']=='G' and
                              ambient_union & (1<<r['root']))
            inc = {j: tuple(core[k] for k in range(3) if cm[k] & (1<<occurrences[j]['root'])) for j in present}
            row = {'core':core, 'edge_children':edges, 'occurrence_indices':present,
                   'G_occurrence_indices':g_present, 'common_G_occurrence_indices':g_common,
                   'G_template_indices':tuple(sorted({occurrences[j]['reference_index'] for j in g_present})),
                   'edge_added_G_occurrence_indices':g_edge, 'ambient_G_occurrence_indices':g_ambient,
                   'incidence':inc}
            candidate_rows.append(row)
            print(f"{str(core):>21} {len(present):7d} {len(g_present):7d} "
                  f"{len(g_common):9d} {len(g_edge):8d} {len(g_ambient):10d}")

        gs = [r for r in occurrences if r['sector']=='G']
        print(f'FIRST {min(len(gs),display_examples)} G-TEMPLATE OCCURRENCES (root birth order; all saved):')
        for r in gs[:display_examples]:
            print(f"  root={r['root']}, boundaries={r['anchors']}, {r['reference_label']}; "
                  f"local support={r['local_support']}; full ancestry={r['root_full_ancestry_size']}; "
                  f"nested boundaries={r['nested_anchors']}")
        if not gs: print('  none')
        outputs.append({'seed':seed, 'audit_N':end, 'coverage':coverage,
                        'occurrences':occurrences, 'anchor_rows':scan['anchor_rows'],
                        'candidates':candidate_rows})
    if fingerprint() != before:
        raise AssertionError('An input history, earlier readout or RNG changed.')
    print('\nLocal templates are relative to TWO specified boundary objects; their older histories remain native.')
    print('All eligible realized source pairs are searched; no absent construction is supplied.')
    print('Same pattern at different boundaries = different anchored occurrence, not a new particle identity.')
    print('Reference labels orient the earlier boundary as a; exchanging boundaries relabels templates, not sectors.')
    print('Coverage across several boundaries is NOT a complete 137-entry registry at any one boundary.')
    print('G labels are inherited TEMPLATE classes, not a demonstrated gravity response or full (45,21,7) connectivity.')
    print('This generalizes F\'s readout; it does not overturn F\'s sparse original-primitive result.')
    print('No mass, 136/27/half multiplier, signed transfer dynamics, new growth or changed service rule.')
    return {'observable':'actual boundary-relative two-source registry templates',
            'reference':reference, 'discovery_population':discovery,
            'audit_population':end, 'branches':outputs}


if not all(k in globals() for k in ('rosetta_calibration_growth_4096', 'rosetta_registry_anatomy_results')):
    raise RuntimeError('Keep the current growth snapshots and original Calibration-F results in this kernel.')
rosetta_local_registry_results = rosetta_local_registry_scan(
    rosetta_calibration_growth_4096, rosetta_registry_anatomy_results)
