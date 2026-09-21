# Calibration N — fixed five-family registry-placement comparison.
# SAME exact-closure samples: K3, K4, K5, induced C4, and induced diamond (K4 minus one edge).
# K3/K4/K5 are RMR-motivated trial labels in this co-parent readout.
# C4 and diamond are NEW exploratory structural candidates, not identified species.
# No new samples, global pair completion, mass fitting, growth, clocks, or dynamics.
from collections import Counter, defaultdict
from itertools import combinations
from statistics import median
import hashlib
import pickle


ROSETTA_N_FAMILIES = ('K3', 'K4', 'K5', 'C4', 'diamond')


def rosetta_n_digest(*sources):
    class Sink:
        def __init__(self):
            self.h = hashlib.sha256()
        def write(self, data):
            self.h.update(data)
            return len(data)
    sink = Sink()
    pickle.Pickler(sink, protocol=5).dump(sources)
    return sink.h.digest()


def rosetta_n_cliques(adjacency, order):
    """Every clique of order 3..5; no maximality or sector-based selection."""
    if order not in (3, 4, 5):
        raise ValueError('Declared clique orders are 3, 4, and 5.')
    out = []
    def extend(prefix, eligible):
        need = order - len(prefix)
        if not need:
            out.append(prefix)
            return
        for offset, vertex in enumerate(eligible):
            if len(eligible) - offset < need:
                break
            tail = tuple(v for v in eligible[offset+1:] if v in adjacency[vertex])
            extend(prefix + (vertex,), tail)
    extend((), tuple(sorted(adjacency)))
    return tuple(out)


def rosetta_n_four_patterns(adjacency):
    """Induced C4 and diamond cores, once each, using actual common neighbors.

    Diamond: one unique shared edge, with two NONadjacent common neighbors.
    C4: two NONadjacent opposite pairs, counted only in one pair ordering.
    Missing diagonals mean absent in THIS reconstructed host, never in universal closure.
    """
    common = defaultdict(list)
    for middle in sorted(adjacency):
        for left, right in combinations(sorted(adjacency[middle]), 2):
            common[left, right].append(middle)
    output = {'C4': [], 'diamond': []}
    for (left, right), centers in sorted(common.items()):
        connected = right in adjacency[left]
        for a, b in combinations(centers, 2):
            if b in adjacency[a]:
                continue
            core = tuple(sorted((left, right, a, b)))
            if connected:
                output['diamond'].append(core)
            elif (left, right) < (a, b):
                output['C4'].append(core)
    for key, values in output.items():
        if len(values) != len(set(values)):
            raise AssertionError('A four-node core was counted more than once in one host.')
        output[key] = tuple(sorted(values))
    return output


def rosetta_registry_five_family_comparison(exact, observations, baseline, placement,
                                           max_examples=2):
    """Registry location in a fixed panel on the UNCHANGED Calibration-I observations.

    1. Reconstruct each root ancestry's co-parent graph from its own children only.
    2. Detect motifs WITHOUT registry/sector labels; main cores exclude primitives.
    3. Support = union of core ancestries plus the actual required edge children.
    4. Keep core-inherited vs edge-only registry entries, identity incidences, blocks,
       G near/skip anatomy, nesting, raw counts, and per-edge descriptive counts.
    5. Compare each family to K4 within equal support-size and last-required-round strata.

    A unique occurrence key contains (family, core, required_edges), because the SAME
    four roots can have different visible missing diagonals in different hosts.
    Each appearance keeps its own host and external-registry context. No absent child
    is inserted. No G, near/skip, workload, or target value selects a candidate.
    """
    if observations.get('scope') != 'exact uniform exhaustive-population observations, not native service histories':
        raise ValueError('Use the saved Calibration-I exact observations.')
    if baseline.get('scope') != 'same-observation internally witnessed K3/K4 registry contrast; provisional templates':
        raise ValueError('Use the saved Calibration-M comparison for regression checks.')
    if placement.get('scope') != 'read-only G structure and placement on the existing exact-closure observations':
        raise ValueError('Use the saved Calibration-J anatomy.')
    if type(max_examples) is not int or max_examples < 0:
        raise ValueError('max_examples must be a nonnegative integer.')
    before = rosetta_n_digest(exact, observations, baseline, placement)
    ref, base = exact['reference'], exact['full5']
    rows = ref['rows']; sector_masks = ref['sector_masks']
    parents = tuple(None if p is None else tuple(sorted(p)) for p in observations['parents'])
    masks, heights = tuple(observations['registry_masks']), tuple(observations['heights'])
    if len(rows) != 137 or tuple(r['index'] for r in rows) != tuple(range(137)):
        raise ValueError('Keep the complete reference registry and its index ordering.')
    if parents[:2] != (None, None) or len(parents) != len(masks) or len(parents) != len(heights):
        raise ValueError('Incomplete exact-closure arrays.')
    if parents[:len(base['parents'])] != tuple(None if p is None else tuple(sorted(p)) for p in base['parents']):
        raise ValueError('The observation cache must retain its H reference prefix.')
    roots = {ref['index_by_object'][m]: z for z,m in enumerate(base['mapped7'])
             if m in ref['index_by_object']}
    if set(roots) != set(range(137)):
        raise ValueError('All reference roots must be present.')
    index_by_root = {v:k for k,v in roots.items()}
    seen_pairs = set()
    for z,p in enumerate(parents):
        if z < 2:
            if p is not None or masks[z] or heights[z]:
                raise ValueError('Primitive convention changed.')
            continue
        if p is None or len(p)!=2 or not 0 <= p[0] < p[1] < z or p in seen_pairs:
            raise ValueError('Invalid or repeated parent construction.')
        seen_pairs.add(p)
        own = 1 << index_by_root[z] if z in index_by_root else 0
        if masks[z] != (masks[p[0]] | masks[p[1]] | own):
            raise ValueError('A saved registry mask disagrees with parentage.')
        if heights[z] != 1 + max(heights[v] for v in p):
            raise ValueError('A saved exhaustive-round grading disagrees with parentage.')

    # Block labels come from branch inheritance, never cache-ID ordering or inferred space.
    x0 = ref['pairs'][0,1]
    pa = ref['pairs'][tuple(sorted((0,x0)))]; pb = ref['pairs'][tuple(sorted((1,x0)))]
    block_by_i = {}
    for i,row in enumerate(rows):
        ss = ref['ancestors'][row['object_id']]
        aa,bb = pa in ss, pb in ss
        if not aa and not bb:
            raise AssertionError('Reference entry is outside the verified three-block partition.')
        block_by_i[i] = 'U' if aa and bb else ('A' if aa else 'B')
    block_masks = {b:sum(1<<i for i in range(137) if block_by_i[i]==b) for b in ('A','B','U')}
    if tuple(block_masks[b].bit_count() for b in ('A','B','U'))!=(63,63,11):
        raise AssertionError('Reference block counts changed.')
    closure_by_i = {a['index']:a['closure'] for a in placement['G_anatomy']}
    g_indices = {r['index'] for r in rows if r['sector']=='G'}
    if set(closure_by_i)!=g_indices or set(closure_by_i.values())!={'near','skip'}:
        raise ValueError('Keep all sixteen existing G anatomical labels.')

    def sig(mask):
        return tuple((mask & sector_masks[s]).bit_count() for s in ('S','I','G'))
    def block_sig(mask):
        return tuple((mask & block_masks[b]).bit_count() for b in ('A','B','U'))
    def indices(mask):
        while mask:
            bit=mask & -mask
            yield bit.bit_length()-1
            mask ^= bit
    ancestry_cache = {}
    def ancestry(root):
        if root not in ancestry_cache:
            seen,pending=set(),[root]
            while pending:
                v=pending.pop()
                if v in seen: continue
                seen.add(v);pending.extend(parents[v] or ())
            ancestry_cache[root]=frozenset(seen)
        return ancestry_cache[root]
    def literal(support):
        return sum(1<<i for i,r in roots.items() if r in support)

    motifs=[]; motif_index={}
    expected_edges={'K3':3,'K4':6,'K5':10,'C4':4,'diamond':5}
    def record(family,core,adj,witness):
        required=tuple(p for p in combinations(core,2) if p[1] in adj[p[0]])
        if len(required)!=expected_edges[family]:
            raise AssertionError('A candidate has the wrong induced edge count.')
        missing=tuple(p for p in combinations(core,2) if p not in required)
        key=(family,core,required)
        children=tuple(witness[e] for e in required)
        if key in motif_index:
            mi=motif_index[key]
            if motifs[mi]['edge_children']!=children: raise AssertionError('Witness identity changed.')
            return mi
        core_support=frozenset().union(*(ancestry(v) for v in core))
        edge_only=frozenset(children)-core_support
        support=core_support|frozenset(children)
        cm,em,fm=literal(core_support),literal(edge_only),literal(support)
        inherited=0; common=(1<<137)-1
        for v in core: inherited |= masks[v]; common &= masks[v]
        if inherited!=cm or cm & em or cm | em != fm:
            raise AssertionError('Registry core/edge partition failed.')
        if any(not ancestry(c)<=support for c in children):
            raise AssertionError('A required edge child is missing ancestry from the motif support.')
        # A missing diagonal cannot sneak in through a core's own ancestry.
        if any(parents[c] in missing for c in support if parents[c] is not None):
            raise AssertionError('An allegedly missing diagonal has an internal witness.')
        loc=[]
        for i in indices(fm):
            r=roots[i]
            member_indices=tuple(k for k,v in enumerate(core) if masks[v] & (1<<i))
            loc.append({'index':i,'label':rows[i]['label'],'sector':rows[i]['sector'],
                'block':block_by_i[i],'closure':closure_by_i.get(i), 'root':r,
                'role':'core ancestry' if cm & (1<<i) else 'edge witness',
                'inheriting_core_members':tuple(core[k] for k in member_indices),
                'core_incidence_positions':member_indices,
                'witnessed_core_pair':parents[r] if r in edge_only else None})
        edge_marks=[]
        for (u,v),c in zip(required,children):
            i=index_by_root.get(c)
            edge_marks.append({'pair':(u,v),'child':c,
                'label':None if i is None else rows[i]['label'],
                'sector':None if i is None else rows[i]['sector'],
                'block':None if i is None else block_by_i[i],
                'closure':closure_by_i.get(i), 'already_in_core':c in core_support})
        # Signature ignores registry labels and core ordering; it is a COARSE summary,
        # not an isomorphism certificate. A/B are exchange-normalized only in this summary.
        def role_signature(swap=False):
            return tuple(sorted((a['sector'], 'B' if swap and a['block']=='A' else
                'A' if swap and a['block']=='B' else a['block'], a['role'],
                a['closure'] or '-',len(a['inheriting_core_members'])) for a in loc))
        mi=len(motifs);motif_index[key]=mi
        motifs.append({'family':family,'order':len(core),'core':core,'required_edges':required,
            'missing_pairs_in_host':missing,'edge_children':children,'edge_count':len(children),
            'core_support':tuple(sorted(core_support)),'edge_only':tuple(sorted(edge_only)),
            'support':tuple(sorted(support)),'support_size':len(support),
            'birth_round':max(heights[v] for v in support),
            'max_core_birth_round':max(heights[v] for v in core),
            'nested_core_pairs':sum(a in ancestry(b) or b in ancestry(a) for a,b in combinations(core,2)),
            'core_registry_mask':cm,'edge_registry_mask':em,'registry_mask':fm,
            'common_registry_mask':common,'core_SIG':sig(cm),'edge_SIG':sig(em),'full_SIG':sig(fm),
            'core_blocks_ABU':block_sig(cm),'edge_blocks_ABU':block_sig(em),'full_blocks_ABU':block_sig(fm),
            'G_near_skip':tuple(sum(a['sector']=='G' and a['closure']==c for a in loc) for c in ('near','skip')),
            'direct_G_edge_count':sum(e['sector']=='G' for e in edge_marks),
            'core_incidence_SIG':tuple(sig(masks[v]) for v in core),
            'registry_locations':tuple(loc),'edge_marks':tuple(edge_marks),
            'coarse_role_signature':min(role_signature(),role_signature(True))})
        return mi

    frames={}; hosts=[]; by_round=defaultdict(list)
    old_hosts={h['sample_index']:h for h in baseline['hosts']}
    if set(old_hosts)!=set(range(len(observations['samples']))):
        raise ValueError('Calibration M and I must contain the same sample rows.')
    verified_old=0
    for si,sample in enumerate(observations['samples']):
        root=sample['root'];rd=sample['round']
        ss=ancestry(root)
        if len(ss)!=sample['ancestry_size'] or literal(ss)!=sample['registry_mask']:
            raise ValueError('Sample frame disagrees with its saved ancestry/mask.')
        if root not in frames:
            adj={v:set() for v in ss};witness={}
            for child in sorted(ss):
                if parents[child] is not None:
                    a,b=parents[child];adj[a].add(b);adj[b].add(a);witness[a,b]=child
            all_patterns={f'K{k}':rosetta_n_cliques(adj,k) for k in (3,4,5)}
            all_patterns.update(rosetta_n_four_patterns(adj))
            primary={f:tuple(c for c in cc if all(v>=2 for v in c)) for f,cc in all_patterns.items()}
            four_faces={c for k in primary['K4'] for c in combinations(k,3)}
            five_tri_faces={c for k in primary['K5'] for c in combinations(k,3)}
            five_tet_faces={c for k in primary['K5'] for c in combinations(k,4)}
            recs={f:[] for f in ROSETTA_N_FAMILIES}
            for family in ROSETTA_N_FAMILIES:
                for core in primary[family]:
                    mi=record(family,core,adj,witness)
                    if not set(motifs[mi]['support'])<=ss:
                        raise AssertionError('An outside child supplied a candidate edge.')
                    recs[family].append({'motif':mi,
                        'face_of_K4_in_host':family=='K3' and core in four_faces,
                        'contained_in_K5_in_host':(family=='K3' and core in five_tri_faces) or
                                                  (family=='K4' and core in five_tet_faces),
                        'external_registry_mask':masks[root] & ~motifs[mi]['registry_mask']})
            frames[root]={'motifs':{f:tuple(v) for f,v in recs.items()},
                'all_core_primitive_counts':{f:dict(Counter(sum(v<2 for v in c) for c in cc))
                                             for f,cc in all_patterns.items()}}
        frame=frames[root]
        host={'sample_index':si,'round':rd,'root':root,'rank':sample['rank'],
              'ancestry_size':len(ss),'registry_mask':masks[root], 'counts_SIG':sig(masks[root]),
              'motifs':frame['motifs'],'all_core_primitive_counts':frame['all_core_primitive_counts']}
        old=old_hosts[si]
        if (old['round'],old['root'],old['rank'])!=(rd,root,sample['rank']):
            raise ValueError('Calibration M refers to a different sample ordering.')
        for k in (3,4):
            newer={motifs[a['motif']]['core']:a for a in host['motifs'][f'K{k}']}
            older={baseline['motifs'][a['motif']]['core']:a for a in old['motifs'][k]}
            if set(newer)!=set(older) or frame['all_core_primitive_counts'][f'K{k}']!=old['all_core_primitive_counts'][k]:
                raise AssertionError('Existing M clique census changed.')
            for core,a in newer.items():
                b=older[core];nm=motifs[a['motif']];om=baseline['motifs'][b['motif']]
                for key in ('edge_children','support','core_support','edge_only','birth_round',
                            'core_registry_mask','edge_registry_mask','registry_mask','common_registry_mask',
                            'core_SIG','edge_SIG','full_SIG'):
                    if nm[key]!=om[key]: raise AssertionError('Existing M motif measurement changed: '+key)
                if a['face_of_K4_in_host']!=b['face_of_K4_in_host']:
                    raise AssertionError('Existing triangle/tetrahedron nesting changed.')
                verified_old+=1
        hosts.append(host);by_round[rd].append(host)

    def summary(ids):
        rr=[motifs[i] for i in sorted(ids)];n=len(rr)
        def total(key,length=3): return tuple(sum(r[key][j] for r in rr) for j in range(length))
        def avg(vals):return tuple(x/n for x in vals) if n else None
        full=total('full_SIG');cs=total('core_SIG');es=total('edge_SIG');bs=total('full_blocks_ABU')
        ng=sum(r['full_SIG'][2]>0 for r in rr)
        edge_den=sum(r['edge_count'] for r in rr);gedges=sum(r['direct_G_edge_count'] for r in rr)
        regden=sum(full)
        return {'motifs':n,'median_support':median(r['support_size'] for r in rr) if n else None,
            'median_birth_round':median(r['birth_round'] for r in rr) if n else None,
            'support_has_G':ng,'support_G_fraction':ng/n if n else None,
            'core_has_G':sum(r['core_SIG'][2]>0 for r in rr),
            'edge_adds_G':sum(r['edge_SIG'][2]>0 for r in rr),
            'distinct_G_labels':tuple(sorted({a['label'] for r in rr for a in r['registry_locations'] if a['sector']=='G'})),
            'sum_SIG':full,'mean_SIG':avg(full),'mean_core_SIG':avg(cs),'mean_edge_SIG':avg(es),
            'mean_blocks_ABU':avg(bs),'G_near_skip_occurrences':total('G_near_skip',2),
            'direct_G_edges':gedges,'edge_denominator':edge_den,
            'direct_G_edge_fraction':gedges/edge_den if edge_den else None,
            'pooled_SIG_fraction':tuple(x/regden for x in full) if regden else None,
            'registry_empty':sum(not r['registry_mask'] for r in rr),
            'birth_round_counts':dict(sorted(Counter(r['birth_round'] for r in rr).items())),
            'max_core_birth_round_counts':dict(sorted(Counter(r['max_core_birth_round'] for r in rr).items())),
            'core_old_enough_for_literal_G':sum(r['max_core_birth_round']>=min(heights[roots[i]] for i in g_indices) for r in rr),
            'distinct_coarse_role_signatures':len({r['coarse_role_signature'] for r in rr})}
    reports=[];matched=[];background=[]
    for rd in observations['sample_rounds']:
        hh=by_round[rd]
        if not hh: raise ValueError('A declared observation round has no samples.')
        background.append({'round':rd,'draws':len(hh),
            'any_G':sum(bool(h['registry_mask'] & sector_masks['G']) for h in hh),
            'mean_SIG':tuple(sum(h['counts_SIG'][j] for h in hh)/len(hh) for j in range(3))})
        unique={f:{a['motif'] for h in hh for a in h['motifs'][f]} for f in ROSETTA_N_FAMILIES}
        for f in ROSETTA_N_FAMILIES:
            row=summary(unique[f]);row.update({'round':rd,'family':f,
                'unique_motif_ids':tuple(sorted(unique[f])),
                'hosts':sum(bool(h['motifs'][f]) for h in hh),
                'appearances':sum(len(h['motifs'][f]) for h in hh),
                'face_appearances':sum(a['face_of_K4_in_host'] for h in hh for a in h['motifs'][f]),
                'inside_K5_appearances':sum(a['contained_in_K5_in_host'] for h in hh for a in h['motifs'][f])})
            reports.append(row)
        for family in ('K3','K5','C4','diamond'):
            groups=defaultdict(lambda:{family:set(),'K4':set()})
            for f in (family,'K4'):
                for mi in unique[f]:
                    mm=motifs[mi];groups[mm['support_size'],mm['birth_round']][f].add(mi)
            shared=[]
            for key,gg in sorted(groups.items()):
                if gg[family] and gg['K4']:
                    shared.append({'support_size':key[0],'motif_birth_round':key[1],
                        'candidate':summary(gg[family]),'K4':summary(gg['K4']),
                        'candidate_ids':tuple(sorted(gg[family])),'K4_ids':tuple(sorted(gg['K4']))})
            matched.append({'round':rd,'family':family,'shared_strata':tuple(shared),
                'unmatched':{f:sum(len(g[f]) for g in groups.values() if not g[family] or not g['K4'])
                             for f in (family,'K4')}})
    latest=observations['sample_rounds'][-1]
    final={r['family']:r for r in reports if r['round']==latest}
    # Cross-family identical COARSE signatures are retained, never treated as distinct masses.
    sig_classes=defaultdict(lambda:defaultdict(set))
    for f in ROSETTA_N_FAMILIES:
        for mi in final[f]['unique_motif_ids']:
            sig_classes[motifs[mi]['coarse_role_signature']][f].add(mi)
    cross={f:{'nonempty_signatures':len({motifs[i]['coarse_role_signature'] for i in final[f]['unique_motif_ids']
                                      if motifs[i]['registry_mask']}),
              'motifs_sharing_nonempty_signature_with_K4':sum(
                  bool(motifs[i]['registry_mask']) and bool(sig_classes[motifs[i]['coarse_role_signature']].get('K4'))
                  for i in final[f]['unique_motif_ids'])}
           for f in ('K3','K5','C4','diamond')}
    def fmt(x):return '--' if x is None else '/'.join(f'{v:.3f}' for v in x)
    def pct(x,n):return '--' if not n else f'{100*x/n:.2f}'
    print('FIVE-FAMILY REGISTRY CONTRAST — K3/K4/K5 plus two exploratory non-clique candidates')
    print('Same exact-closure samples. K4 is a trial vacuum template, not identified physical vacuum.')
    print('C4 = induced four-node ring; diamond = induced K4 with one edge absent in THIS host.')
    print('Original primitives excluded from main cores; all primitive-count censuses retained.')
    print(f'Calibration-M regression: {verified_old} existing sample-motif appearances agree exactly.')
    print('\nALL SAMPLED ROUNDS — no new draws; unique motifs and appearances have different denominators')
    print(' round   family hosts unique appearances medDAG coreG edgeG     G%      mean S/I/G')
    for r in reports:
        md='--' if r['median_support'] is None else f'{r["median_support"]:g}'
        print(f'{r["round"]:6d} {r["family"]:>8} {r["hosts"]:5d} {r["motifs"]:6d} {r["appearances"]:11d}'
              f' {md:>6} {r["core_has_G"]:5d} {r["edge_adds_G"]:5d}'
              f' {pct(r["support_has_G"],r["motifs"]):>6} {fmt(r["mean_SIG"]):>17}')
    print(f'\nLATEST U{latest} — location and anatomy; NOT a mass conversion or enrichment significance test')
    print(' family       core S/I/G       edge S/I/G   mean A/B/U    G near/skip  G edges/all edges')
    for f,r in final.items():
        print(f'{f:>8} {fmt(r["mean_core_SIG"]):>16} {fmt(r["mean_edge_SIG"]):>16}'
              f' {fmt(r["mean_blocks_ABU"]):>14} {str(r["G_near_skip_occurrences"]):>14}'
              f' {r["direct_G_edges"]}/{r["edge_denominator"]} ({pct(r["direct_G_edges"],r["edge_denominator"])}%)')
    print('G near/skip and edge counts count reused entry appearances across DISTINCT motifs, not independent entries.')
    print('Distinct G identities across each entire family (NOT all inside one motif):',
          {f:len(r['distinct_G_labels']) for f,r in final.items()})
    print('\nCONSTRUCTION-ROUND CONTROL: a later host does not imply a newly formed motif')
    for f,r in final.items():
        print(f' {f}: required birth rounds={r["birth_round_counts"]}; cores old enough to contain literal G='
              f'{r["core_old_enough_for_literal_G"]}/{r["motifs"]}')
    g_births={heights[roots[i]] for i in g_indices}
    print(' All literal G birth rounds:',sorted(g_births))
    print(' Every core vertex participates in an edge: motif completed by round 5 => all core vertices born by 4.')
    print(' Thus a G-free core in such an early-completed motif is forced by birth order, not a functional exclusion.')
    print('\nLATEST: equal-support-size / equal-last-required-round comparisons against K4')
    for panel in (p for p in matched if p['round']==latest):
        print(f' {panel["family"]} versus K4:')
        print('   DAG round    candidate/G+          K4/G+    G% candidate/K4')
        for s in panel['shared_strata']:
            a,b=s['candidate'],s['K4']
            print(f'  {s["support_size"]:4d} {s["motif_birth_round"]:5d}'
                  f' {a["motifs"]:9d}/{a["support_has_G"]:<5d} {b["motifs"]:9d}/{b["support_has_G"]:<5d}'
                  f' {pct(a["support_has_G"],a["motifs"]):>8}/{pct(b["support_has_G"],b["motifs"])}')
        if not panel['shared_strata']:print('   No shared stratum: absence of comparability, not equality.')
        print('   unmatched retained:',panel['unmatched'])
    print('\nCOARSE PLACEMENT SIGNATURE OVERLAP WITH K4, latest round')
    print('Signature retains sector, A/B/U block, core-or-edge role, G closure, number of inheriting cores;')
    print('A/B exchange-normalized. It discards labels and detailed core/edge arrangements; equality is NOT isomorphism.')
    for f,v in cross.items():
        print(f' {f}: nonempty signatures={v["nonempty_signatures"]}; motifs matching a nonempty K4 signature='
              f'{v["motifs_sharing_nonempty_signature_with_K4"]}/{final[f]["motifs"]}')
    print('\nNESTING, latest sample appearances (local flags, not rejected candidates)')
    for f in ('K3','K4'):
        r=final[f]
        print(f' {f}: appearances={r["appearances"]}; face of K4={r["face_appearances"]};'
              f' contained in K5={r["inside_K5_appearances"]}')
    print('\nFIRST G-BEARING SUPPORTS by required birth round, then core and edge IDs')
    for family in ROSETTA_N_FAMILIES:
        positive=sorted((motifs[i] for i in final[family]['unique_motif_ids'] if motifs[i]['full_SIG'][2]),
                        key=lambda m:(m['birth_round'],m['core'],m['required_edges']))
        print(f' {family}: {len(positive)} G-bearing supports; first {min(max_examples,len(positive))}:')
        for mm in positive[:max_examples]:
            print(f'  core={mm["core"]}; children={mm["edge_children"]}; support={mm["support_size"]}; round={mm["birth_round"]}')
            if mm['missing_pairs_in_host']: print('   unwitnessed pairs in this host:',mm['missing_pairs_in_host'])
            for a in mm['registry_locations']:
                if a['sector']=='G':
                    print(f'   {a["label"]}: {a["role"]}; block={a["block"]}; {a["closure"]};'
                          f' pair={a["witnessed_core_pair"]}; carriers={a["inheriting_core_members"]}')
    print('\nAll families fixed before this scan; C4/diamond are exploratory comparisons, not named species.')
    print('No missing diagonal is erased from the universe. An induced non-clique can change class in a larger host.')
    print('Required birth round is NOT the time a non-clique first appears or disappears in a host.')
    print('The cache is not the ambient universe. Every edge is witnessed inside the individual observed ancestry.')
    print('Registry roots and candidate motifs may overlap. Unique counts and host appearances are not independent particles.')
    print('Matched strata control support size and last required round only; no p-value, classifier, or physical separation claimed.')
    print('These finite exact-closure rounds are not established modern epochs. No mass, energy, spatial clock, or dynamics assigned.')
    if rosetta_n_digest(exact, observations, baseline, placement)!=before:
        raise AssertionError('A source result changed.')
    return {'scope':'fixed five-family internally witnessed registry comparison; structural candidates only',
        'families':ROSETTA_N_FAMILIES,'edge_rule':'individual ancestral co-parent reconstruction',
        'nonclique_rule':'induced C4 and diamond; no deletions or global absence assertions',
        'source_observer_seed':observations['observer_seed'],'sample_rounds':observations['sample_rounds'],
        'motifs':motifs,'hosts':hosts,'round_reports':reports,'background_reports':background,
        'matched_strata':matched,'coarse_signature_overlap_with_K4':cross,
        'reference_blocks':block_by_i,'verified_M_appearances':verified_old}


if not all(k in globals() for k in ('rosetta_constructor_reference','rosetta_exact_closure_samples',
                                  'rosetta_registry_motif_contrast','rosetta_G_structure_placement')):
    raise RuntimeError('Keep the existing H, I, J and M results. No growth or new sample is needed.')
rosetta_registry_candidate_panel = rosetta_registry_five_family_comparison(
    rosetta_constructor_reference, rosetta_exact_closure_samples,
    rosetta_registry_motif_contrast, rosetta_G_structure_placement)
