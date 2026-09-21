# CELL 27 — constrained nucleon-partner search: exact nulls, then one-port substitutions.
# Paste after Cell 26 (only Cell 1 definitions/data are execution prerequisites).
# Standard library only. No network. No u/d, charge, mass, energy or epoch is assigned.
# T remains the declared registry-overlap readout, NOT the native service transition.
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import gcd


def _run_dcu_mass_27():
    need = ('_reference', '_panel', '_parents', '_Native')
    absent = [key for key in need if key not in globals()]
    if absent:
        raise RuntimeError('Run Cell 1 first. Missing: ' + ', '.join(absent))
    ref, panel, parents = _reference, _panel, _parents
    before = deepcopy((ref, panel, parents))
    Q = Fraction
    rows = ref['rows']
    N = len(rows)
    adjacency = tuple(frozenset(x) for x in ref['adjacency'])
    degree = tuple(map(len, adjacency))
    if N != 137 or Counter(r['sector'] for r in rows) != Counter(S=81, I=40, G=16):
        raise ValueError('The original registry is required; no replacement is selected.')
    c = ref['pairs'][(0, 1)]
    pa, pb = (ref['pairs'][tuple(sorted((p, c)))] for p in (0, 1))
    blocks = {'U': [], 'A': [], 'B': []}
    for i, row in enumerate(rows):
        aa = ref['ancestors'][row['object_id']]
        a, b = pa in aa, pb in aa
        assert a or b
        blocks['U' if a and b else 'A' if a else 'B'].append(i)
    assert tuple(len(blocks[x]) for x in ('U', 'A', 'B')) == (11, 63, 63)
    block_of = {i: name for name, group in blocks.items() for i in group}
    for i in range(N):
        expected = {j for j in range(N) if j != i and
                    (block_of[i] == block_of[j] or 'U' in (block_of[i], block_of[j]))}
        assert adjacency[i] == expected
    block_sectors = {name: dict(Counter(rows[i]['sector'] for i in group))
                     for name, group in blocks.items()}

    # A. Powers of the ACTUAL unweighted operator, using only exact integers.
    # M = L*T has integer entries. No underflow at k=8 is mistaken for a zero.
    L = 1
    for d in degree:
        L = L*d//gcd(L, d)
    factor = tuple(L//d for d in degree)
    assert L == 9928
    def apply_scaled(vector):
        return tuple(-sum(factor[j]*vector[j] for j in adjacency[i]) for i in range(N))
    columns = tuple(tuple(int(i == j) for i in range(N)) for j in range(N))
    powers = {0: columns}
    times = (1, 2, 4, 8)
    for k in range(1, max(times)+1):
        columns = tuple(apply_scaled(col) for col in columns)
        if k in times:
            powers[k] = columns
    sector_pairs = Counter()
    categories = Counter()
    power_counts = {k: Counter() for k in times}
    example_cross = None
    component_checks = 0
    for i, j in combinations(range(N), 2):
        bi, bj = block_of[i], block_of[j]
        if bi == bj:
            kind = 'universal_twins' if bi == 'U' else 'same_lobe'
        else:
            kind = 'opposite_lobes' if 'U' not in (bi, bj) else 'junction_lobe'
        categories[kind] += 1
        if kind == 'same_lobe':
            sector_pairs[''.join(sorted((rows[i]['sector'], rows[j]['sector'])))] += 1
        for k in times:
            x = tuple(a-b for a,b in zip(powers[k][i], powers[k][j]))
            outside = sum(abs(v) for vtx,v in enumerate(x) if vtx not in (i,j))
            if outside == 0:
                power_counts[k]['zero_external_all_pairs'] += 1
                power_counts[k][kind] += 1
            if bi == bj:
                expected = tuple((int(v == i)-int(v == j))*factor[i]**k for v in range(N))
                assert x == expected
                component_checks += N
            elif kind == 'opposite_lobes':
                # Exact decomposition: within-lobe difference + antisymmetric block mode.
                local = factor[i]**k
                block = (-(len(blocks[bi])-1)*factor[i])**k
                numerator = block-local
                assert numerator % len(blocks[bi]) == 0
                macro = numerator//len(blocks[bi])
                expected = tuple(local*(int(v == i)-int(v == j)) +
                                 macro*(int(block_of[v] == bi)-int(block_of[v] == bj))
                                 for v in range(N))
                assert x == expected and all(x[u] == 0 for u in blocks['U'])
                component_checks += N
                if example_cross is None:
                    example_cross = (i,j)
            if bi != 'U' and bj != 'U':
                assert all(x[u] == 0 for u in blocks['U'])
    assert categories == Counter(same_lobe=3906, universal_twins=55,
                                 opposite_lobes=3969, junction_lobe=1386)
    assert sector_pairs == Counter(GG=56, GI=304, GS=576, II=342, IS=1368, SS=1260)
    for k in times:
        assert power_counts[k] == Counter(zero_external_all_pairs=3961,
                                         same_lobe=3906, universal_twins=55)

    # B. Native primitive exchange, kept distinct from flavor or entry reflection.
    # This is an involution of complete parentage, NOT an edit of an existing object.
    @lru_cache(None)
    def term(v, reflected=False):
        if v < 2:
            return 'ab'[1-v if reflected else v]
        return '(' + '|'.join(sorted(term(p, reflected) for p in parents[v])) + ')'
    @lru_cache(None)
    def ancestors(v):
        if v < 2:
            return frozenset((v,))
        return frozenset((v,)) | ancestors(parents[v][0]) | ancestors(parents[v][1])
    @lru_cache(None)
    def chains(v):
        return 1 if v < 2 else sum(chains(p) for p in parents[v])
    @lru_cache(None)
    def height(v):
        return 0 if v < 2 else 1+max(height(p) for p in parents[v])
    cache_pairs = {tuple(sorted(p)):v for v,p in enumerate(parents) if v >= 2}
    cc = cache_pairs[(0,1)]
    aa,bb = (cache_pairs[tuple(sorted((p,cc)))] for p in (0,1))
    def branch(v):
        a,b = aa in ancestors(v), bb in ancestors(v)
        return 'U' if a and b else 'A' if a else 'B' if b else 'base'
    lookup = {term(v):v for v in range(len(parents))}
    motifs = panel['motifs']
    report = next(r for r in panel['round_reports'] if (r['round'],r['family']) == (12,'K3'))
    ids = tuple(report['unique_motif_ids'])
    core_lookup = {tuple(sorted(motifs[i]['core'])):i for i in ids}
    assert len(core_lookup) == len(ids)
    membership = defaultdict(set)
    hosts = {}
    for h in panel['hosts']:
        if h['round'] == 12:
            hosts[h['sample_index']] = h
            for item in h['motifs']['K3']:
                membership[item['motif']].add(h['sample_index'])
    signatures = Counter(tuple(sorted(branch(v) for v in motifs[i]['core'])) for i in ids)
    pool, misses = [], []
    for i in ids:
        m = motifs[i]
        if sorted(branch(v) for v in m['core']) != ['A','A','B']:
            continue
        for u in m['core']:
            if branch(u) != 'A':
                continue
            v = lookup.get(term(u,True))
            if v is None:
                misses.append(dict(motif=i,port=u,reason='mirror port not in decoded cache'))
                continue
            if v in m['core']:
                misses.append(dict(motif=i,port=u,reason='mirror port already a constituent'))
                continue
            core2 = tuple(sorted((set(m['core'])-{u}) | {v}))
            j = core_lookup.get(core2)
            if j is None:
                misses.append(dict(motif=i,port=u,reason='partner not in saved round-12 K3 catalogue'))
                continue
            other = motifs[j]
            matched = (m['support_size'],m['birth_round']) == (other['support_size'],other['birth_round'])
            shared_hosts = tuple(sorted(membership[i] & membership[j]))
            global_mirror = tuple(sorted(term(z,True) for z in m['core'])) == \
                            tuple(sorted(term(z) for z in core2))
            assert chains(u) == chains(v) and height(u) == height(v)
            assert sorted(branch(z) for z in core2) == ['A','B','B']
            pool.append(dict(left=i,right=j,replaced_port=u,replacement_port=v,
                fixed_ports=tuple(sorted(set(m['core'])-{u})),
                left_core=tuple(m['core']),right_core=core2,
                left_support_size=m['support_size'],right_support_size=other['support_size'],
                left_grade=m['birth_round'],right_grade=other['birth_round'],
                matched_size_grade=matched,common_host_indices=shared_hosts,
                full_primitive_mirror=global_mirror))
    matched = [r for r in pool if r['matched_size_grade']]
    primary = [r for r in matched if r['common_host_indices']]
    def order(r):
        return (r['left_grade'],r['left_support_size'],
                tuple(sorted(term(v) for v in r['left_core'])),term(r['replaced_port']))
    pool.sort(key=order); matched.sort(key=order); primary.sort(key=order)
    # Recover each complete witness support and certify actual common-host presence.
    for r in matched:
        for i in (r['left'],r['right']):
            m = motifs[i]
            assert set(map(tuple,m['required_edges'])) == set(combinations(sorted(m['core']),2))
            assert all(tuple(sorted(parents[w])) == tuple(edge)
                       for edge,w in zip(m['required_edges'],m['edge_children']))
            needed = set().union(*(ancestors(v) for v in m['core'])) | set(m['edge_children'])
            assert needed == set(m['support'])
        for h in r['common_host_indices']:
            host_support = ancestors(hosts[h]['root'])
            assert set(motifs[r['left']]['support']) | set(motifs[r['right']]['support']) <= host_support

    # C. SAME repeat-work observable Q and the old complete 32-chain probe family.
    # Every included nonprimitive support port is contacted once/pass; raw work
    # is compared with an alternative primitive/probe arm, not charged twice.
    words = tuple(''.join(w) for w in product('ab',repeat=5))
    histories = contacts = 0
    def repeat_Q(i, word):
        nonlocal histories, contacts
        m = motifs[i]
        if m['birth_round'] > 5:
            raise ValueError('Fixed grade-6 probes are not valid for this candidate.')
        st, mapping = _Native(), {0:0,1:1}
        preparation_work = 0
        for old in sorted(m['support']):
            if old >= 2:
                born,work = st.add_batch([tuple(mapping[p] for p in parents[old])])
                mapping[old] = born[0]; preparation_work += work
        p = st.pair_to_id[(0,1)]
        for letter in word:
            pair = tuple(sorted((p,'ab'.index(letter))))
            if pair not in st.pair_to_id:
                _,work=st.add_batch([pair]); preparation_work += work
            p = st.pair_to_id[pair]
        assert len(st.ancestors[p]) == 8 and st.path_weight(p) == 12
        port_ids = tuple(sorted((v for v in m['support'] if v >= 2), key=term))
        ports = tuple(mapping[v] for v in port_ids)
        support = frozenset(mapping.values())
        overlap = support & st.ancestors[p]
        old_parentage = tuple(st.parents[v] for v in sorted(support))
        _,warm_work = st.add_batch([(v,p) for v in ports])
        contacts += len(ports)
        assert support <= st.recorded
        results=[]; port_outputs=[]; excluded_extension=[]
        for _ in range(2):
            # Same first-letter successor convention as Cell 4. These grade>=7
            # successors cannot add a grade<=5 candidate to probe ancestry.
            born,work=st.add_batch([('ab'.index(word[0]),p)]); p=born[0]; excluded_extension.append(work)
            assert overlap == support & st.ancestors[p]
            values=[]
            for v in ports:
                baseline=st.recording_cost([(0,p)])
                raw=st.recording_cost([(v,p)])
                born,committed=st.add_batch([(v,p)])
                delta=raw-baseline
                assert raw == committed and delta == 2*sum(st.path_weight(z)
                        for z in st.ancestors[v]-st.ancestors[p])
                values.append(delta); contacts += 1
            results.append(Q(sum(values),len(values)))
            port_outputs.append(tuple(values))
        assert results[0] == results[1] and port_outputs[0] == port_outputs[1]
        assert old_parentage == tuple(st.parents[v] for v in sorted(support))
        histories += 1
        return dict(Q=results[0],port_ids=port_ids,ports=port_outputs[0],warm_work=warm_work,
                    preparation_work_excluded=preparation_work,
                    extension_work_excluded=tuple(excluded_extension),two_repeat_passes_equal=True)
    selected = sorted({i for r in matched for i in (r['left'],r['right'])})
    measured = {i:{word:repeat_Q(i,word) for word in words} for i in selected}
    # Original electron-anchor apparatus values are regression tests ONLY.
    assert repeat_Q(271,'aaaaa')['Q'] == Q(248,9)
    assert repeat_Q(271,'bbbbb')['Q'] == Q(308,9)
    def evaluate_pair(r):
        a,b=measured[r['left']],measured[r['right']]
        deltas={word:b[word]['Q']-a[word]['Q'] for word in words}
        even=sum(deltas.values(),Q(0))/len(words)  # Context average, NOT a new mass observable.
        if r['full_primitive_mirror']:
            for word in words:
                reflected = word.translate(str.maketrans('ab','ba'))
                assert a[word]['Q'] == b[reflected]['Q']
        return dict(deltas=deltas,delta_A=deltas['aaaaa'],delta_B=deltas['bbbbb'],
                    minimum=min(deltas.values()),maximum=max(deltas.values()),
                    signs=tuple(sum((v>0,v==0,v<0)[j] for v in deltas.values()) for j in range(3)),
                    context_mean=even,context_mean_is_mass=False)
    for r in matched:
        r['response']=evaluate_pair(r)
    primary_pattern=Counter(r['response']['signs'] for r in primary)
    primary_means=Counter(r['response']['context_mean'] for r in primary)
    other_means=Counter(r['response']['context_mean'] for r in matched if not r['common_host_indices'])
    # Checkpoint expectations validate this exact archive, not a mass target.
    assert (len(ids),len(pool),len(matched),len(primary)) == (1753,109,82,41)
    assert primary_pattern == Counter({(16,0,16):41})
    assert primary_means == Counter({Q(0):41})
    assert len(selected)==134
    assert (ref,panel,parents)==before

    print('CELL 27 — CONSTRAINED PARTNERS; NO PROTON OR NEUTRON IDENTIFICATION')
    print('Unchanged registry T=-A D^-1; response steps are NOT maintenance ticks.')
    print('Exact block sizes:',{k:len(v) for k,v in blocks.items()})
    print('Sector pair    same-lobe zero-external pairs')
    for key in ('GG','GI','GS','II','IS','SS'):
        print(f'  {key:6s} {sector_pairs[key]:8d}')
    print('Total: 3906/7875 peripheral pairs. Universal-block twins add 55 more.')
    print('\n power k    peripheral zero-external    total zero-external    junction response from peripheral dipoles')
    for k in times:
        print(f' {k:7d} {power_counts[k]["same_lobe"]:27d} '
              f'{power_counts[k]["zero_external_all_pairs"]:22d}                  exact zero')
    print('Same-lobe T^k(e_i-e_j)=73^-k(e_i-e_j), no exceptions in any sector.')
    print('Opposite-lobe differences use eigenvalues 1/73 and -62/73; no small new eigenvalue appears.')
    print('No graph, degree, sector label, response weight or coupling changed.')
    print('\nONE-PORT PRIMITIVE-EXCHANGE SEARCH (A/B are ancestry branches, NOT quark flavors)')
    print(f'{len(ids)} saved round-12 K3s; {len(pool)} AAB -> ABB one-port correspondences;')
    print(f'  {len(matched)} match complete support size and required grade; '
          f'{len(primary)} coexist in at least one actual saved host.')
    print(f'  {sum(r["full_primitive_mirror"] for r in primary)} of the co-host pairs are full primitive mirrors.')
    print(f'All {len(primary)} co-host pairs: 16 positive / 0 zero / 16 negative deltas over the same 32 probes.')
    print('Their context-averaged difference is EXACTLY zero in all cases; average is diagnostic only.')
    print('Matched pairs without a common saved host, context means:',dict(sorted(other_means.items())))
    first=primary[0]
    print('\nFirst co-host pair by grade, support size and canonical parentage (not by response):')
    print(f'  archive {first["left"]} core {first["left_core"]} -> '
          f'{first["right"]} core {first["right_core"]}')
    print(f'  only port {first["replaced_port"]} -> {first["replacement_port"]}; '
          f'D={first["left_support_size"]}, grade={first["left_grade"]}')
    for side in ('left','right'):
        i=first[side]
        print(f'  {side}: Q_A={measured[i]["aaaaa"]["Q"]}, Q_B={measured[i]["bbbbb"]["Q"]}, '
              f'core terms={tuple(term(v) for v in motifs[i]["core"])}')
    print(f'  deltas A/B={first["response"]["delta_A"]}/{first["response"]["delta_B"]}.')
    print(f'PASS: {len(list(combinations(range(N),2)))*len(times):,} exact pair-power checks; '
          f'{histories:,} native probe fixtures; {contacts:,} committed support contacts.')
    print('PASS: complete witnesses, real common-host containment, legal fresh pairs, repeat-work identity, mirror controls.')
    print('Primitive exchange is a structural correspondence, NOT a realized beta transition.')
    print('Native Q remains conditional recording work; no physical charge, mass, binding, or mature epoch supplied.')
    return dict(protocol='fixed overlap symmetry + target-free one-port branch-partner catalogue',
        blocks={k:tuple(v) for k,v in blocks.items()},block_sectors=block_sectors,
        recurrence=dict(denominator=L,steps=times,categories=dict(categories),
            zero_external_sector_pairs=dict(sector_pairs),power_counts={k:dict(v) for k,v in power_counts.items()},
            exact_component_checks=component_checks,same_lobe_eigenvalue=Q(1,73),
            cross_lobe_antisymmetric_eigenvalue=-Q(62,73),peripheral_junction_response_zero=True),
        selection=dict(round=12,K3_count=len(ids),signature_counts=dict(signatures),
            raw_one_port_pairs=len(pool),matched_pairs=len(matched),cohost_pairs=len(primary),
            failures=misses,unmatched=[r for r in pool if not r['matched_size_grade']]),
        matched_pairs=matched,primary_pairs=primary,first_example=first,
        words=words,measurements=measured,probe_fixtures=histories,committed_support_contacts=contacts,
        signed_change_patterns=dict(primary_pattern),context_mean_counts=dict(primary_means),
        noncohost_context_mean_counts=dict(other_means),
        physical_particle_identified=False,u_d_mapping=None,charge_observable=None,
        mass_calibration_changed=False,native_service_history_generated=False,
        candidate_selection_used_mass_target=False,constructor_modified=False,
        response_steps_are_native_ticks=False,source_mass_dictionary_unchanged=True)


dcu_mass_27 = _run_dcu_mass_27()
