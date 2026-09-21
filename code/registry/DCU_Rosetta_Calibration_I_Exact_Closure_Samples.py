# Calibration I — exact uniform observations and explicit lineages in exhaustive closure.
# This is NOT a faster Gamma=2 simulation. It supplies no service, backlog or clock.
# Reference U5 is reused; higher objects are decoded only when queried.
import hashlib
import math
import pickle
import random
from collections import Counter, defaultdict
from statistics import mean, median


class RosettaExactClosureReader:
    """Exact rank/unrank for U_(s+1) = {a,b} plus unordered distinct pairs in U_s.

    Ranks are local to a specified exhaustive round. Semantic object identities
    are interned parent pairs; a repeated object gets the SAME local identity,
    even when its ranks at different rounds differ. No observed subset feeds back
    into the sampling population. Every draw uses the full exact population size.
    """
    def __init__(self, exact, max_round=12):
        if type(max_round) is not int or not 5 <= max_round <= 14:
            raise ValueError('Use 5 <= max_round <= 14 for this bounded reader.')
        base = exact['full5']
        if len(base['parents']) != 2280 or len(exact['reference']['rows']) != 137:
            raise ValueError('Use the verified Calibration-H complete U5 and reference.')
        self.parents = list(base['parents'])
        self.heights = list(base['heights'])
        self.masks = list(base['registry_masks'])
        self.pairs = dict(base['pairs'])
        self.populations = [2]
        for _ in range(max_round):
            n = self.populations[-1]
            self.populations.append(2 + n * (n-1) // 2)
        self.rank_cache, self.unrank_cache = {}, {}
        for s in range(6):
            members = [z for z, h in enumerate(self.heights) if h <= s]
            if len(members) != self.populations[s]:
                raise ValueError('The supplied U5 has an inconsistent round grading.')
            ranks = [self.rank(s, z) for z in members]
            if set(ranks) != set(range(self.populations[s])):
                raise AssertionError('The base rank map is not a bijection.')

    def _round(self, s):
        if type(s) is not int or not 0 <= s < len(self.populations):
            raise ValueError('Round is outside the declared finite reader range.')

    def rank(self, s, z):
        self._round(s)
        if type(z) is not int or not 0 <= z < len(self.parents) or self.heights[z] > s:
            raise ValueError('Object does not belong to the requested exhaustive round.')
        key = (s, z)
        if key in self.rank_cache:
            return self.rank_cache[key]
        if z < 2:
            k = z
        else:
            a, b = self.parents[z]
            i, j = sorted((self.rank(s-1, a), self.rank(s-1, b)))
            k = 2 + j*(j-1)//2 + i  # colex order: all i<j, then increase j
        if not 0 <= k < self.populations[s]:
            raise AssertionError('Invalid exhaustive rank.')
        old = self.unrank_cache.get((s, k))
        if old is not None and old != z:
            raise AssertionError('Two semantic objects have the same rank.')
        self.rank_cache[key] = k
        self.unrank_cache[s, k] = z
        return k

    def _intern(self, a, b):
        a, b = sorted((a, b))
        if a == b:
            raise ValueError('Self-pairs are excluded.')
        old = self.pairs.get((a, b))
        if old is not None:
            return old
        height = 1 + max(self.heights[a], self.heights[b])
        if height <= 5:
            raise AssertionError('An object is missing from the complete U5 base.')
        z = len(self.parents)
        self.parents.append((a, b))
        self.heights.append(height)
        # All literal DAG-seven entries already exist in U5. Later roots add no own bit.
        self.masks.append(self.masks[a] | self.masks[b])
        self.pairs[a, b] = z
        return z

    def unrank(self, s, k):
        self._round(s)
        if type(k) is not int or not 0 <= k < self.populations[s]:
            raise ValueError('Rank must be an integer inside the full population.')
        old = self.unrank_cache.get((s, k))
        if old is not None:
            return old
        if k < 2:
            z = k
        else:
            t = k - 2
            j = (1 + math.isqrt(1 + 8*t)) // 2
            i = t - j*(j-1)//2
            if not 0 <= i < j < self.populations[s-1]:
                raise AssertionError('Invalid unordered-pair unranking.')
            z = self._intern(self.unrank(s-1, i), self.unrank(s-1, j))
        self.unrank_cache[s, k] = z
        prior = self.rank_cache.get((s, z))
        if prior is not None and prior != k:
            raise AssertionError('One semantic object has two ranks in a round.')
        self.rank_cache[s, z] = k
        return z

    def ancestry(self, root):
        seen, pending = set(), [root]
        while pending:
            z = pending.pop()
            if z in seen:
                continue
            seen.add(z)
            if self.parents[z] is not None:
                pending.extend(self.parents[z])
        return seen

    def uniform_other(self, s, root, rng, require_new_child=False):
        """Uniform partner in U_s minus root, with an optional new-child condition.

        If root is older than round s, a child first born at s+1 requires a partner
        from the newest shell. Rejection implements exactly that explicit condition.
        """
        excluded = self.rank(s, root)
        while True:
            k = rng.randrange(self.populations[s]-1)
            if k >= excluded:
                k += 1
            q = self.unrank(s, k)
            if q == root:
                raise AssertionError('Excluded parent was sampled.')
            if not require_new_child or max(self.heights[root], self.heights[q]) == s:
                return q, k


def rosetta_exact_closure_observations(exact, samples_per_round=128,
                                     rounds=(6, 8, 10, 12), lineage_links=5,
                                     observer_seed=20260920):
    """Observe exact exhaustive populations; do not evolve a service-driven history.

    Uniform samples are WITH replacement. Lineages choose uniformly from eligible
    partners in the FULL exhaustive population, not from previously decoded nodes.
    The 137 fixed registry starts all participate; no cost or physical target selects
    a start or a partner. Retaining the initial registry bit is inherited by definition.
    """
    if (type(samples_per_round) is not int or not 1 <= samples_per_round <= 512 or
        type(lineage_links) is not int or not 1 <= lineage_links <= 7 or
        type(observer_seed) is not int or not rounds or
        any(type(s) is not int or not 5 <= s <= 12 for s in rounds) or
        len(set(rounds)) != len(rounds)):
        raise ValueError('Use distinct rounds 5..12, 1..512 draws, and 1..7 lineage links.')
    if not callable(globals().get('rosetta_h_exact_marginal_rounds')):
        raise RuntimeError('Keep the Calibration-H1 functions in this kernel.')

    def digest():
        class Sink:
            def __init__(self): self.h = hashlib.sha256()
            def write(self, x): self.h.update(x); return len(x)
        sink = Sink()
        pickle.Pickler(sink, protocol=5).dump(exact)
        return sink.h.digest()
    before = digest()
    global_rng_before = random.getstate()
    reference, base = exact['reference'], exact['full5']
    end_round = max(max(rounds), 5+lineage_links, 6)
    reader = RosettaExactClosureReader(exact, end_round)
    expected = {r['round']: r for r in rosetta_h_exact_marginal_rounds(
        base, exact['round6'], reference, final_round=end_round)}
    sample_rng = random.Random(observer_seed)
    path_rng = random.Random(observer_seed+1)
    sector_masks = reference['sector_masks']
    roots = {reference['index_by_object'][m]: z for z, m in enumerate(base['mapped7'])
             if m in reference['index_by_object']}
    if len(roots) != 137:
        raise AssertionError('All original registry roots must be present in U5.')

    # Extra check of the existing exact counts: sector vs construction-round effect.
    cohorts = defaultdict(list)
    for i, z in roots.items():
        cohorts[reader.heights[z], reference['rows'][i]['sector']].append(i)
    control = []
    for (height, sector), indices in sorted(cohorts.items()):
        counts = {expected[10]['per_entry_descendants'][i] for i in indices} if 10 in expected else set()
        control.append({'birth_round':height, 'sector':sector, 'entries':len(indices),
                        'round10_distinct_single_counts':len(counts),
                        'round10_single_count':next(iter(counts)) if len(counts)==1 else None})
    for height in {h for h, _ in cohorts}:
        indices = [i for (h, _), ii in cohorts.items() if h == height for i in ii]
        for row in expected.values():
            if len({row['per_entry_descendants'][i] for i in indices}) != 1:
                raise AssertionError('Same-birth registry marginal counts differ.')

    print('EXACT CLOSURE OBSERVATIONS — NOT the Gamma=2 occurrence process')
    print(f'Private observer RNG seed={observer_seed}; existing native RNGs are not used.')
    print('Uniform draws WITH replacement from complete U_s; only queried ancestries are decoded.')
    print('   round draws  mean R sample/exact      any R% sample/exact     any G% sample/exact  median ancestry')
    samples, summaries = [], []
    for s in rounds:
        rows = []
        for _ in range(samples_per_round):
            rank = sample_rng.randrange(reader.populations[s])
            z = reader.unrank(s, rank)
            support = reader.ancestry(z)
            mask = reader.masks[z]
            direct = sum(1 << i for i, r in roots.items() if r in support)
            if mask != direct:
                raise AssertionError('Inherited registry mask differs from literal ancestry.')
            rows.append({'round':s, 'rank':rank, 'root':z, 'birth_round':reader.heights[z],
                         'ancestry_size':len(support), 'registry_mask':mask,
                         'registry_count':mask.bit_count(),
                         'counts_SIG':tuple((mask & sector_masks[t]).bit_count() for t in ('S','I','G'))})
        n = expected[s]['N']
        summary = {'round':s, 'draws':len(rows), 'mean_registry':mean(r['registry_count'] for r in rows),
                   'any_registry_fraction':sum(bool(r['registry_mask']) for r in rows)/len(rows),
                   'any_G_fraction':sum(r['counts_SIG'][2]>0 for r in rows)/len(rows),
                   'median_ancestry':median(r['ancestry_size'] for r in rows),
                   'exact_mean_registry':expected[s]['registry_incidences']/n,
                   'exact_any_registry_fraction':expected[s]['registry_carriers']/n,
                   'exact_any_G_fraction':expected[s]['sector_carriers']['G']/n}
        samples.extend(rows); summaries.append(summary)
        print(f"{s:8d} {len(rows):5d}  {summary['mean_registry']:8.3f}/{summary['exact_mean_registry']:<8.3f}"
              f"  {100*summary['any_registry_fraction']:8.3f}/{100*summary['exact_any_registry_fraction']:<8.3f}"
              f"  {100*summary['any_G_fraction']:8.3f}/{100*summary['exact_any_G_fraction']:<8.3f}"
              f"  {summary['median_ancestry']:15.1f}")

    paths = []
    for i in range(137):
        current = roots[i]
        nodes, coparents = [current], []
        for s in range(5, 5+lineage_links):
            q, qrank = reader.uniform_other(s, current, path_rng, require_new_child=True)
            child = reader._intern(current, q)
            if reader.heights[child] != s+1 or not reader.masks[child] & (1 << i):
                raise AssertionError('Lineage violated birth-round or inheritance identity.')
            if reader.unrank(s+1, reader.rank(s+1, child)) != child:
                raise AssertionError('Lineage rank round-trip failed.')
            coparents.append({'round':s, 'root':q, 'rank':qrank})
            nodes.append(child); current = child
        paths.append({'reference_index':i, 'label':reference['rows'][i]['label'],
                      'sector':reference['rows'][i]['sector'], 'roots':tuple(nodes),
                      'coparents':tuple(coparents),
                      'registry_counts':tuple(reader.masks[z].bit_count() for z in nodes),
                      'counts_SIG':tuple(tuple((reader.masks[z] & sector_masks[t]).bit_count()
                                              for t in ('S','I','G')) for z in nodes),
                      'ancestry_sizes':tuple(len(reader.ancestry(z)) for z in nodes)})
    print(f'\nLINEAGES: all 137 registry starts at U5; {lineage_links} genuinely new-child links each.')
    print('Partner choices are uniform among qualifying FULL-U_s partners, not among the decoded sample.')
    print('Initial entry retained in every path: exact inheritance, NOT evidence of preferential stabilization.')
    for sector in ('S','I','G'):
        pp = [p for p in paths if p['sector']==sector]
        finals = [p['registry_counts'][-1] for p in pp]
        print(f"  {sector}: starts={len(pp)}; final registry count min/median/max="
              f"{min(finals)}/{median(finals):g}/{max(finals)}")
    print('\nSAME-BIRTH MARGINAL CONTROL (existing exact recurrence):')
    for row in control:
        print(f"  first round {row['birth_round']}: sector={row['sector']}, entries={row['entries']}")
    print('Within each birth-round group, ALL entries have exactly equal descendant counts at every computed round.')
    print('This is a marginal-count statement, not equality of the full weighted/size-conditioned structures.')
    print(f'\nDecoded union contains {len(reader.parents)} objects including the 2280-object reference base.')
    print('This union is an observer cache, NOT the population or a chronological growth history.')
    print('New object IDs are cache labels; they are NOT IDs in the old sampled simulations.')
    print('Any R means at least ONE registry entry; 100.000% may be rounded and never means all 137.')
    print('No work, maintenance ticks, seconds, mass, forces, or unique particle paths are assigned.')
    if digest()!=before or random.getstate()!=global_rng_before:
        raise AssertionError('An input reference or the global random state changed.')
    return {'scope':'exact uniform exhaustive-population observations, not native service histories',
            'observer_seed':observer_seed, 'sample_rounds':tuple(rounds),
            'samples':samples, 'summaries':summaries, 'lineages':paths,
            'birth_round_control':control, 'parents':tuple(reader.parents),
            'heights':tuple(reader.heights), 'registry_masks':tuple(reader.masks),
            'populations':tuple(reader.populations)}


if 'rosetta_constructor_reference' not in globals():
    raise RuntimeError('Run Calibration H1 in this kernel first. No additional native growth is required.')
rosetta_exact_closure_samples = rosetta_exact_closure_observations(rosetta_constructor_reference)
