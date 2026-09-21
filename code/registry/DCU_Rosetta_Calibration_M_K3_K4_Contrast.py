# Calibration M — K3 / K4 registry location in the SAME exact-closure observations.
# RMR motivates K3 and K4 as trial matter/vacuum templates; neither is identified
# physically by this scan. This is NOT the registry-overlap graph or an ambient clique scan.
# All edges must be witnessed by actual children in each sampled root's OWN ancestry.
from collections import Counter, defaultdict
from itertools import combinations
from statistics import median
import hashlib
import pickle


def rosetta_m_digest(*objects):
    class Sink:
        def __init__(self):
            self.h = hashlib.sha256()
        def write(self, data):
            self.h.update(data)
            return len(data)
    sink = Sink()
    pickle.Pickler(sink, protocol=5).dump(objects)
    return sink.h.digest()


def rosetta_m_cliques(adjacency, order):
    """Enumerate every order-3/4 clique once; primitive vertices are not removed here."""
    if order not in (3, 4):
        raise ValueError('This readout tests triangles and tetrahedral complete graphs only.')
    result = []
    def extend(prefix, eligible):
        need = order - len(prefix)
        if not need:
            result.append(prefix)
            return
        if len(eligible) < need:
            return
        for offset, vertex in enumerate(eligible):
            tail = tuple(w for w in eligible[offset+1:] if w in adjacency[vertex])
            extend(prefix + (vertex,), tail)
    extend((), tuple(sorted(adjacency)))
    return tuple(result)


def rosetta_registry_k3_k4_contrast(exact, observations, placement=None, max_examples=2):
    """Selection uses only internally witnessed co-parent K3/K4, NEVER sector labels.

    The primary cores contain no ORIGINAL primitive, preserving Calibration-J's
    nonprimitive-core convention. All-clique primitive-count censuses are retained.
    A motif support is union(Anc(core)) plus the actual pair-witness children.
    Registry entries are split into inherited CORE ancestry and EDGE-WITNESS-only
    additions. This is NOT a force, bound-state, physical vacuum, or mass definition.
    All K3 faces of K4 are retained and flagged per host, never deleted after results.
    Distinct motifs and (sample draw, motif) appearances have separate denominators.
    The same object may appear in many sampled ancestries. Counts are not independent
    particle counts. The decoded union is only a cache, never a background universe.
    """
    if observations.get('scope') != 'exact uniform exhaustive-population observations, not native service histories':
        raise ValueError('Use the unchanged Calibration-I observations.')
    if type(max_examples) is not int or max_examples < 0:
        raise ValueError('max_examples must be a nonnegative integer.')
    before = rosetta_m_digest(exact, observations, placement)
    ref, base = exact['reference'], exact['full5']
    rows = ref['rows']
    parents = tuple(None if q is None else tuple(sorted(q)) for q in observations['parents'])
    masks, heights = tuple(observations['registry_masks']), tuple(observations['heights'])
    if len(rows) != 137 or any(row['index'] != i for i, row in enumerate(rows)):
        raise ValueError('Use the same complete reference catalogue and row order.')
    if parents[:2] != (None, None) or len(parents) != len(masks) or len(parents) != len(heights):
        raise ValueError('Incomplete two-primitive observation cache.')
    if parents[:len(base['parents'])] != tuple(base['parents']):
        raise ValueError('Exact-closure base prefix differs from Calibration H.')
    roots = {ref['index_by_object'][m]: z for z, m in enumerate(base['mapped7'])
             if m in ref['index_by_object']}
    if set(roots) != set(range(137)):
        raise ValueError('The exhaustive base must contain all reference roots.')
    index_by_root = {z: i for i, z in roots.items()}
    sector_masks = ref['sector_masks']
    pairs = {}
    for z, pp in enumerate(parents):
        if z < 2:
            if pp is not None or masks[z] or heights[z]:
                raise ValueError('Primitive metadata disagree.')
            continue
        if pp is None or len(pp) != 2 or not 0 <= pp[0] < pp[1] < z or pp in pairs:
            raise ValueError('Invalid, duplicate, or non-topological parentage.')
        pairs[pp] = z
        own = 1 << index_by_root[z] if z in index_by_root else 0
        if masks[z] != (masks[pp[0]] | masks[pp[1]] | own):
            raise ValueError('Registry mask disagrees with parentage.')
        if heights[z] != 1 + max(heights[p] for p in pp):
            raise ValueError('Birth-round grading disagrees with parentage.')

    def sig(mask):
        return tuple((mask & sector_masks[s]).bit_count() for s in ('S', 'I', 'G'))
    def labels(mask):
        return tuple(r['label'] for r in rows if mask & (1 << r['index']))
    ancestry_cache = {}
    def ancestry(root):
        if root not in ancestry_cache:
            seen, pending = set(), [root]
            while pending:
                v = pending.pop()
                if v in seen:
                    continue
                seen.add(v)
                pending.extend(parents[v] or ())
            ancestry_cache[root] = frozenset(seen)
        return ancestry_cache[root]
    def literal(support):
        return sum(1 << i for i, root in roots.items() if root in support)

    motifs, motif_index = [], {}
    def record_motif(core, witness):
        key = (len(core), core)
        if key in motif_index:
            mi = motif_index[key]
            if motifs[mi]['edge_children'] != witness:
                raise AssertionError('One core acquired inconsistent edge witnesses.')
            return mi
        core_support = frozenset().union(*(ancestry(v) for v in core))
        edge_only = frozenset(witness) - core_support
        support = core_support | frozenset(witness)
        core_mask, full_mask = literal(core_support), literal(support)
        expected_core = 0
        for v in core:
            expected_core |= masks[v]
        if expected_core != core_mask:
            raise AssertionError('Core registry readout failed independent ancestry check.')
        for child in witness:
            if not ancestry(child) <= support:
                raise AssertionError('Motif support is not ancestry closed.')
        edge_mask = literal(edge_only)
        if core_mask & edge_mask or (core_mask | edge_mask) != full_mask:
            raise AssertionError('Core/edge registry partition failed.')
        common = (1 << 137) - 1
        for v in core:
            common &= masks[v]
        mi = len(motifs)
        motif_index[key] = mi
        motifs.append({'order': len(core), 'core': core, 'edge_children': witness,
            'core_support': tuple(sorted(core_support)), 'edge_only': tuple(sorted(edge_only)),
            'support': tuple(sorted(support)), 'support_size': len(support),
            'birth_round': max(heights[v] for v in support),
            'nested_core_pairs': sum(a in ancestry(b) or b in ancestry(a) for a,b in combinations(core,2)),
            'core_registry_mask': core_mask, 'edge_registry_mask': edge_mask,
            'registry_mask': full_mask, 'common_registry_mask': common,
            'core_SIG': sig(core_mask), 'edge_SIG': sig(edge_mask), 'full_SIG': sig(full_mask),
            'registry_locations': tuple({'index': i, 'label': row['label'], 'sector': row['sector'],
                'root': roots[i], 'role': 'core ancestry' if core_mask & (1 << i) else 'edge witness',
                'inheriting_core_members': tuple(v for v in core if masks[v] & (1 << i)),
                'witnessed_core_pair': parents[roots[i]] if roots[i] in edge_only else None}
                for i,row in enumerate(rows) if full_mask & (1 << i))})
        return mi

    hosts, by_round = [], defaultdict(list)
    for sample_index, sample in enumerate(observations['samples']):
        rd, root = sample['round'], sample['root']
        support = ancestry(root)
        if len(support) != sample['ancestry_size'] or literal(support) != sample['registry_mask']:
            raise ValueError('An observation disagrees with its actual ancestral frame.')
        adjacency = {v: set() for v in support}
        witness_by_pair = {}
        for child in support:
            if parents[child] is not None:
                a, b = parents[child]
                adjacency[a].add(b); adjacency[b].add(a)
                witness_by_pair[a,b] = child
        all_cliques = {k: rosetta_m_cliques(adjacency, k) for k in (3,4)}
        core_lists = {k: tuple(c for c in all_cliques[k] if all(v >= 2 for v in c)) for k in (3,4)}
        tetra_faces = {face for tetra in core_lists[4] for face in combinations(tetra,3)}
        records = {3: [], 4: []}
        for k in (3,4):
            for core in core_lists[k]:
                edge_children = tuple(witness_by_pair[p] for p in combinations(core,2))
                mi = record_motif(core, edge_children)
                # This role flag is local to THIS host, not a global deletion of triangles.
                records[k].append({'motif': mi, 'face_of_K4_in_host': k == 3 and core in tetra_faces,
                                   'external_registry_mask': masks[root] & ~motifs[mi]['registry_mask']})
        host = {'sample_index': sample_index, 'round': rd, 'root': root, 'rank': sample['rank'],
            'ancestry_size': len(support), 'registry_mask': masks[root], 'counts_SIG': sig(masks[root]),
            'all_core_primitive_counts': {k: dict(Counter(sum(v < 2 for v in c) for c in all_cliques[k])) for k in (3,4)},
            'motifs': records}
        hosts.append(host); by_round[rd].append(host)

    def summary(ids):
        rr = [motifs[i] for i in sorted(ids)]
        n = len(rr)
        sums = tuple(sum(r['full_SIG'][j] for r in rr) for j in range(3))
        denom = sum(sums)
        return {'motifs': n, 'median_support': median(r['support_size'] for r in rr) if n else None,
            'median_birth_round': median(r['birth_round'] for r in rr) if n else None,
            'core_has_G': sum(r['core_SIG'][2] > 0 for r in rr),
            'edge_adds_G': sum(r['edge_SIG'][2] > 0 for r in rr),
            'support_has_G': sum(r['full_SIG'][2] > 0 for r in rr),
            'sum_SIG': sums, 'mean_SIG': tuple(x/n for x in sums) if n else None,
            'pooled_SIG_fraction': tuple(x/denom for x in sums) if denom else None}
    def fmt_mean(values):
        return '--' if values is None else '/'.join(f'{v:.3f}' for v in values)
    def pct(a,b):
        return '--' if not b else f'{100*a/b:.2f}'
    reports, matched = [], []
    print('K3/K4 REGISTRY CONTRAST — trial RMR structural templates, NOT identified matter/vacuum')
    print('Same Calibration-I samples; each ancestral co-parent graph reconstructed internally.')
    print('Main cores exclude original primitives; all primitive-count censuses retained separately.')
    print('Selection ignores sector labels. All K3 are retained, including faces of K4.')
    print('\nBACKGROUND: original uniform object draws, not a proven vacuum ensemble')
    print(' round draws anyG%   mean S/I/G (distinct registry entries per sampled ancestry)')
    for rd in observations['sample_rounds']:
        hh = by_round[rd]
        vals = tuple(sum(h['counts_SIG'][j] for h in hh)/len(hh) for j in range(3))
        print(f'{rd:6d} {len(hh):5d} {pct(sum(h["counts_SIG"][2]>0 for h in hh),len(hh)):>6}   {fmt_mean(vals)}')
    print('\nDISTINCT MOTIFS WITHIN EACH ROUND SAMPLE: full support = core ancestry + edge witnesses')
    print(' round motif hosts  unique appearances medianDAG coreG edgeG       mean S/I/G')
    for rd in observations['sample_rounds']:
        hh = by_round[rd]
        unique = {k: {q['motif'] for h in hh for q in h['motifs'][k]} for k in (3,4)}
        for k in (3,4):
            row = summary(unique[k])
            row.update({'round':rd,'order':k,'unique_motif_ids':tuple(sorted(unique[k])),
                        'hosts':sum(bool(h['motifs'][k]) for h in hh),
                        'appearances':sum(len(h['motifs'][k]) for h in hh),
                        'face_appearances':sum(q['face_of_K4_in_host'] for h in hh for q in h['motifs'][k])})
            reports.append(row)
            md = '--' if row['median_support'] is None else f'{row["median_support"]:g}'
            print(f'{rd:6d} K{k:1d} {row["hosts"]:7d} {row["motifs"]:7d} {row["appearances"]:11d} '
                  f'{md:>9} {row["core_has_G"]:5d} {row["edge_adds_G"]:5d} {fmt_mean(row["mean_SIG"]):>20}')
        # Compare only equal minimal-support SIZE and same latest required birth ROUND.
        strata = defaultdict(lambda: {3: set(),4: set()})
        for k in (3,4):
            for mi in unique[k]:
                m=motifs[mi]; strata[m['support_size'],m['birth_round']][k].add(mi)
        shared = []
        for key, groups in sorted(strata.items()):
            if not groups[3] or not groups[4]:
                continue
            a,b = summary(groups[3]), summary(groups[4])
            shared.append({'support_size':key[0], 'motif_birth_round':key[1], 'K3':a,'K4':b,
                           'K3_ids':tuple(sorted(groups[3])), 'K4_ids':tuple(sorted(groups[4]))})
        matched.append({'round':rd, 'shared_strata':shared,
            'unmatched':{k:sum(len(g[k]) for g in strata.values() if not g[3] or not g[4]) for k in (3,4)}})
        if placement is not None:
            j = next(r for r in placement['sample_summary'] if r['round'] == rd)
            if (sum(bool(h['motifs'][3]) for h in hh) != j['with_nonprimitive_coparent_K3'] or
                sum(any(motifs[q['motif']]['core_SIG'][2] for q in h['motifs'][3]) for h in hh) != j['with_coparent_K3_inheriting_G']):
                raise AssertionError('Existing Calibration-J K3 result was not reproduced.')
    print('\nLATEST ROUND: equal-support-size / equal-motif-birth-round comparison; no fitted matching radius')
    print(' DAG birthRound     K3/G+       K4/G+     G% K3/K4')
    last=matched[-1]
    for s in last['shared_strata']:
        a,b=s['K3'],s['K4']
        print(f'{s["support_size"]:4d} {s["motif_birth_round"]:10d} '
              f'{a["motifs"]:6d}/{a["support_has_G"]:<3d} {b["motifs"]:6d}/{b["support_has_G"]:<3d} '
              f'{pct(a["support_has_G"],a["motifs"]):>8}/{pct(b["support_has_G"],b["motifs"])}')
    print(' Unmatched distinct motifs retained:',last['unmatched'])
    print('\nFIRST G-BEARING SUPPORTS by motif birth round then core IDs; not ranked by a physical target')
    latest_ids={k:set(r['unique_motif_ids']) for k in (3,4) for r in reports if r['round']==last['round'] and r['order']==k}
    for k in (3,4):
        positive=sorted((motifs[i] for i in latest_ids[k] if motifs[i]['full_SIG'][2]),
                        key=lambda m:(m['birth_round'],m['core']))
        print(f' K{k}: {len(positive)} positive supports; first {min(max_examples,len(positive))}:')
        for m in positive[:max_examples]:
            print(f'  core={m["core"]}; children={m["edge_children"]}; support={m["support_size"]}; birth round={m["birth_round"]}')
            for r in m['registry_locations']:
                if r['sector']=='G':
                    print(f'    {r["label"]} at root={r["root"]}: {r["role"]}; comparison={r["witnessed_core_pair"]}; core carriers={r["inheriting_core_members"]}')
    print('\ncoreG/edgeG count motifs, not registry entries. mean S/I/G counts support entries once per motif.')
    print('Core and edge-only registry sets are disjoint. A G edge witness is a registry ENTRY, not an entire G sector.')
    print('The same motif can occur in several draws and rounds; distinct/appearance counts must not be pooled as independent particles.')
    print('K4 contains four K3 faces. Face flags are retained locally; they do not select a post-hoc particle subset.')
    print('Matched strata control support size and last required construction round only; no significance or population-enrichment claim.')
    print('Motif birth round is structural grading, not time of discovery in a later host sample.')
    print('RMR K3/K4 names motivate a trial mapping to THIS co-parent readout; no physical vacuum, particle identity, energy, or mass is assigned.')
    print('No new samples, native updates, lineages, registry entries, or RNG draws. Inputs unchanged.')
    if rosetta_m_digest(exact,observations,placement)!=before:
        raise AssertionError('A source object changed.')
    return {'scope':'same-observation internally witnessed K3/K4 registry contrast; provisional templates',
            'candidate_rule':'all nonprimitive K3; retain K4 faces', 'vacuum_proxy_rule':'all nonprimitive K4',
            'edge_rule':'actual comparison child lies in each sampled root ancestry',
            'background_rule':'original uniform Calibration-I object draws',
            'motifs':motifs,'hosts':hosts,'round_reports':reports,'matched_strata':matched,
            'source_observer_seed':observations['observer_seed'],'sample_rounds':observations['sample_rounds']}


if not all(name in globals() for name in ('rosetta_constructor_reference', 'rosetta_exact_closure_samples', 'rosetta_G_structure_placement')):
    raise RuntimeError('Keep Calibration H, I and J results in this kernel; no growth is required.')
rosetta_registry_motif_contrast = rosetta_registry_k3_k4_contrast(
    rosetta_constructor_reference, rosetta_exact_closure_samples, rosetta_G_structure_placement)
