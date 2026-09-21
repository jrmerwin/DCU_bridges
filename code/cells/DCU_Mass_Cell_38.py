# CELL 38 — native record information and the registry's conditional-dependence graph.
# Self-contained: Python standard library + NumPy. No old notebook variables or downloads.
# Uses fixed archived parent cards and the native persistent-trit law, not new histories.
# Information is in bits; it is not heat, energy, spatial distance, or particle mass.
import json as _j38, zlib as _z38, base64 as _b38
from pathlib import Path as _Path38

_FIXTURE38 = _j38.loads(_z38.decompress(_b38.b64decode(
'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'
)))

def _run_dcu_mass_38(fixture):
    import numpy as np
    from collections import Counter
    from itertools import combinations
    from math import log2
    from functools import lru_cache
    lg3 = log2(3)
    rng = np.random.Generator(np.random.PCG64(20380302))
    checks = dict(path_checks=0, pair_distributions=0, triple_distributions=0,
                  registry_pair_distributions=0, enumerated_assignments=0)
    max_error = 0.0

    @lru_cache(None)
    def all_trits(n):
        a = np.arange(3**n, dtype=np.int64)[:, None]
        return ((a // (3**np.arange(n, dtype=np.int64))) % 3).astype(np.int8)

    def read_tokens(g, keys, tokens):
        pos = {z: i for i,z in enumerate(keys)}
        value = np.empty((len(g), len(tokens)), dtype=np.int8)
        for j,(z,e) in enumerate(tokens):
            x = g[:, pos[z]]
            value[:,j] = np.where((e == 1) & (x < 2), 1-x, x)
        return value

    def entropy(*arrays, exact=False):
        a = np.column_stack(arrays)
        code = a.astype(np.int64) @ (3**np.arange(a.shape[1], dtype=np.int64))
        _, counts = np.unique(code, return_counts=True)
        if exact:
            assert np.all(counts == counts[0]), 'Full fixed-token distribution must be uniform on its image.'
        p = counts/counts.sum()
        return float(-np.dot(p,np.log2(p)))

    def observe(tokens_list, sample_count=None):
        keys = tuple(sorted({z for ts in tokens_list for z,e in ts}))
        if sample_count is None:
            g = all_trits(len(keys)); checks['enumerated_assignments'] += len(g)
        else:
            g = rng.integers(0,3,size=(sample_count,len(keys)),dtype=np.int8)
        arrays = [read_tokens(g,keys,ts) for ts in tokens_list]
        return arrays

    def mutual(tokens_a,tokens_b, sample_count=None):
        a,b=observe([tokens_a,tokens_b],sample_count)
        exact=sample_count is None
        return entropy(a,exact=exact)+entropy(b,exact=exact)-entropy(a,b,exact=exact)

    def conditional(tokens_a,tokens_b,tokens_c):
        a,b,c=observe([tokens_a,tokens_b,tokens_c])
        return (entropy(a,c,exact=True)+entropy(b,c,exact=True)
                -entropy(c,exact=True)-entropy(a,b,c,exact=True))

    def native_paths(parents,length):
        children={z for pp in parents[2:] for z in pp}
        cache={0:[(0,())],1:[(1,())]}
        for z,pp in enumerate(parents[2:],2):
            cache[z]=[(root,tokens+((z,e),)) for e,p in enumerate(pp)
                      for root,tokens in cache[p] if len(tokens)<length]
        return sorted((root,tokens) for z in sorted(children) if z>=2
                      for root,tokens in cache[z] if len(tokens)==length)

    out={'protocol':'Cell38 record information v1','information_unit':'bits',
         'register_draws':'once per register in each independent record-value realization',
         'new_native_history':False,'physical_space_or_energy_identified':False,
         'panels':{},'registry':{},'formula_log2_3':lg3}
    print('CELL 38: RECORD INFORMATION, CONDITIONAL SHARING, AND THE ACTUAL REGISTRY')
    print('Complete addresses; fixed entry reflection; no energy or spatial valuation.')
    for env,data in fixture['record_panels'].items():
        out['panels'][env]={'checkpoint':data['checkpoint'],'lengths':{}}
        for length_string,old_panel in data['panels'].items():
            length=int(length_string)
            paths=native_paths(data['parentage'],length)
            expected=[(r,tuple(map(tuple,t))) for r,t in old_panel['paths']]
            assert paths==expected
            checks['path_checks']+=len(paths)
            tokens=[p[1] for p in paths]
            supports=[set(z for z,e in ts) for ts in tokens]
            rows=[]; representatives={}; aliases=[]
            for i,j in combinations(range(len(paths)),2):
                shared=len(supports[i]&supports[j]); target=shared*lg3
                obs=mutual(tokens[i],tokens[j]); checks['pair_distributions']+=1
                max_error=max(max_error,abs(obs-target)); assert abs(obs-target)<2e-11
                prefix=mutual(tokens[i][:1],tokens[j][:1])
                assert prefix<=obs+2e-11
                row={'pair':[i,j],'shared_registers':shared,'MI_bits':obs,
                     'MI_trits_exact':shared,'first_token_MI_bits':prefix}
                rows.append(row); representatives.setdefault(shared,(i,j))
                if supports[i]==supports[j]:
                    a,b=observe([tokens[i],tokens[j]])
                    aliases.append(dict(row,address_equality_probability=float(np.all(a==b,axis=1).mean())))
            # Deterministic structural triple selection, not favorable information values.
            selected={}; all_triple_classes=Counter()
            for i,j in combinations(range(len(paths)),2):
                common=supports[i]&supports[j]
                for k in range(len(paths)):
                    if k in (i,j): continue
                    key=(len(common),len(common-supports[k]))
                    all_triple_classes[key]+=1
                    if len(supports[i]|supports[j]|supports[k])<=8:
                        selected.setdefault(key,(i,j,k))
            triples=[]
            for key,(i,j,k) in sorted(selected.items()):
                got=conditional(tokens[i],tokens[j],tokens[k]); predicted=key[1]*lg3
                checks['triple_distributions']+=1
                assert abs(got-predicted)<3e-11; max_error=max(max_error,abs(got-predicted))
                triples.append({'indices':[i,j,k],'pair_shared':key[0],
                                'unobserved_common':key[1],'CMI_bits':got})
            sampled=[]
            for shared,(i,j) in sorted(representatives.items()):
                for shots in (4096,65536):
                    empirical=mutual(tokens[i],tokens[j],shots)
                    sampled.append({'pair':[i,j],'samples':shots,'shared':shared,
                                    'truth_bits':shared*lg3,'plugin_MI_bits':empirical,
                                    'error_bits':empirical-shared*lg3})
            out['panels'][env]['lengths'][length_string]={'paths':paths,'pairs':rows,
                'pair_overlap_counts':dict(sorted(Counter(r['shared_registers'] for r in rows).items())),
                'alias_controls':aliases,'selected_triples':triples,
                'all_ordered_conditioning_class_counts':{str(k):v for k,v in sorted(all_triple_classes.items())},
                'finite_sampling':sampled}
            print(f'  {env:10s} L={length}: {len(paths):3d} paths; {len(rows):4d} pair distributions; '
                  f'{len(triples):2d} selected conditional distributions; {len(aliases):2d} same-support aliases.')

    # All proper ancestors are recorded in the complete size-bounded source ideal.
    r=fixture['registry']; par=r['parents']; anc=[{0},{1}]
    for z,pp in enumerate(par[2:],2): anc.append(anc[pp[0]]|anc[pp[1]]|{z})
    recorded={x for pp in par[2:] for x in pp}
    c=next(z for z,pp in enumerate(par) if pp==[0,1])
    supports=[]; raw=[]
    for row in r['rows']:
        z=row['object_id']; s=anc[z]-{0,1,z}
        assert len(anc[z])==7 and len(s)==4 and s<=recorded and c in s
        supports.append(s); raw.append(tuple((v,0) for v in sorted(s)))
    n=len(supports); counts=Counter(); sector_counts={}; info=[]
    graph=[set() for _ in range(n)]
    for i,j in combinations(range(n),2):
        observed=conditional(raw[i],raw[j],((c,0),)); target=len((supports[i]&supports[j])-{c})
        checks['registry_pair_distributions']+=1
        err=abs(observed-target*lg3); max_error=max(max_error,err); assert err<3e-11
        counts[target]+=1
        sec=''.join(sorted([r['rows'][i]['sector'],r['rows'][j]['sector']]))
        sector_counts.setdefault(sec,Counter())[target]+=1
        if target:
            graph[i].add(j);graph[j].add(i)
        info.append({'i':i,'j':j,'shared_trits':target+1,'conditional_shared_trits':target})
    assert [sorted(x) for x in graph]==r['adjacency']
    universal=[i for i,g in enumerate(graph) if len(g)==n-1]
    rem=set(range(n))-set(universal); components=[]
    while rem:
        front=[min(rem)]; comp=set(front);rem.difference_update(comp)
        while front:
            v=front.pop();new=(graph[v]&rem);rem.difference_update(new);comp.update(new);front.extend(sorted(new))
        components.append(sorted(comp))
    groups=[universal]+sorted(components,key=lambda x:x[0])
    # Full lobe bundles are conditionally independent too, not only individual pairs.
    lobes=[set().union(*(supports[i] for i in comp)) for comp in groups[1:]]
    assert (lobes[0]&lobes[1])=={c}
    registry_alias=Counter(tuple(sorted(s)) for s in supports)
    out['registry']={'source':'original complete 173-object closure, 137 roots at inclusive ancestry size 7',
        'readout':'four proper nonprimitive ancestral register values, not the unresolved root',
        'per_view_entropy_bits':4*lg3,'pairs':info,'conditional_MI_trit_histogram':dict(sorted(counts.items())),
        'sector_pair_histograms':{s:dict(sorted(cn.items())) for s,cn in sorted(sector_counts.items())},
        'conditional_graph_edges':sum(map(len,graph))//2,'exact_original_adjacency_match':True,
        'block_sizes':[len(x) for x in groups],
        'block_sector_compositions':[dict(Counter(r['rows'][i]['sector'] for i in g)) for g in groups],
        'lobe_register_union_sizes':[len(s) for s in lobes],
        'conditional_information_between_full_lobe_bundles_bits':0,
        'distinct_register_bundles':len(registry_alias),
        'bundle_multiplicity_histogram':dict(Counter(registry_alias.values())),
        'universal_indices':universal,'branches':components}
    print('Registry: each view reads four ALREADY RECORDED proper ancestral registers.')
    print(f'  {n} roots; {n*(n-1)//2} pair distributions; I(Yi;Yj | gc)/log2(3) histogram={dict(sorted(counts.items()))}')
    print(f'  Positive conditional-information graph: {sum(map(len,graph))//2} edges; blocks {[len(g) for g in groups]}.')
    print('  Exact match to the old overlap graph; same-input information interpretation, not an independent derivation of its size.')
    print(f'  {len(registry_alias)} distinct register bundles for 137 distinct roots: some final-parent distinctions remain invisible.')
    out['validation']=dict(checks,max_float_entropy_identity_error=max_error)
    print('PASS:',checks,'; max entropy arithmetic residual',max_error)
    print('Finite-sample plug-in bias retained. No registry-specific trit weight, physical space, mass, heat or causal field added.')
    return out

dcu_mass_38 = _run_dcu_mass_38(_FIXTURE38)
_out38 = _Path38.cwd() / 'DCU_Cell_38_outputs'
_out38.mkdir(exist_ok=True)
(_out38 / 'DCU_Mass_Cell_38_RESULTS.json').write_text(_j38.dumps(dcu_mass_38,indent=2)+'\n',encoding='utf-8')
