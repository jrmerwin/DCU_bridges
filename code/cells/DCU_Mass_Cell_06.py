# CELL 6 — registry recruitment under ONE DECLARED COUPLING EXTENSION.
# Run after Cell 5. Standard library only; earlier Q results remain unchanged.
# Native: objects, ancestry, legal recording costs. ADDED: signed transport on
# the undirected parent-child graph, and |signal| as a contact probability.
from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import fsum


def _run_dcu_mass_6():
    needed = ('dcu_mass_5', '_reference', '_Native', '_parents', '_panel', '_term')
    missing = [k for k in needed if k not in globals()]
    if missing:
        raise RuntimeError('Run Cells 1–5 first. Missing: ' + ', '.join(missing))
    catalogue = dcu_mass_5['catalogue']
    motifs, ref = _panel['motifs'], _reference
    horizons, exact_steps = (32, 64, 128), 8  # Fixed numerical windows, not clocks.

    # One COMMON legal host: full reference + all seven saved supports + mirrors.
    # Primitive-exchange closure avoids a preferred a/b preparation.
    host, rmap = _Native(), {0: 0, 1: 1}
    for old, pair in enumerate(ref['parents'][2:], 2):
        born, _ = host.add_batch([tuple(rmap[p] for p in pair)])
        rmap[old] = born[0]

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
    registry = tuple(rmap[r['object_id']] for r in ref['rows'])
    assert len(set(registry)) == len(ref['rows'])
    assert all(len(host.ancestors[v]) == 7 for v in registry)
    graph = [set() for _ in range(n)]
    for child, pair in enumerate(host.parents):
        if pair is not None:
            for parent in pair:
                graph[child].add(parent)
                graph[parent].add(child)
    graph = tuple(tuple(sorted(a)) for a in graph)
    assert sum(map(len, graph)) == 4*(n-2)
    swap = [1, 0]
    for pair in host.parents[2:]:
        swap.append(host.pair_to_id[tuple(sorted(swap[p] for p in pair))])
    assert all({swap[j] for j in graph[i]} == set(graph[swap[i]]) for i in range(n))

    def readout_kernel(branch):
        """Actual legal warming; native costs for alternative single contacts."""
        s = deepcopy(host)
        p = s.pair_to_id[(0, 1)]
        for _ in range(5):
            pair = tuple(sorted((p, branch)))
            if pair not in s.pair_to_id:
                s.add_batch([pair])
            p = s.pair_to_id[pair]
        overlap = set(range(n)) & s.ancestors[p]
        # Warm every original nonprimitive host node with the same apparatus.
        _, warm_work = s.add_batch([(v, p) for v in range(2, n)])
        assert set(range(n)) <= s.recorded
        kernel = None
        for test in range(2):
            p = s.add_batch([(p, branch)])[0][0]  # Fresh legal probe, not a duplicate.
            assert set(range(n)) & s.ancestors[p] == overlap
            baseline = s.recording_cost([(0, p)])
            costs = tuple(s.recording_cost([(v, p)]) - baseline for v in range(n))
            assert costs == tuple(2*sum(s.path_weight(z)
                                       for z in s.ancestors[v] - s.ancestors[p])
                                  for v in range(n))
            if kernel is None:
                kernel = costs
            else:
                assert costs == kernel
            # A real contact does not change the warmed host's later cost kernel.
            target = 2 if test == 0 else n-1
            _, raw = s.add_batch([(target, p)])
            assert raw - baseline == costs[target]
        return kernel, warm_work

    kernels, warming = {}, {}
    for label, branch in (('A', 0), ('B', 1)):
        kernels[label], warming[label] = readout_kernel(branch)
    assert all(kernels['A'][i] == kernels['B'][swap[i]] for i in range(n))
    # Bridge to the OLD observable: the identical kernel averaged over old ports
    # reproduces all 14 original A/B repeat-Q measurements, despite the larger host.
    for row in catalogue:
        i = row['archive_motif_id']
        ports = [embed(v) for v in motifs[i]['support'] if v >= 2]
        for label, word in (('A', 'aaaaa'), ('B', 'bbbbb')):
            assert Fraction(sum(kernels[label][v] for v in ports), len(ports)) == \
                dcu_mass_5['cases'][i, word]['Q']

    def evolve(core, adjacency):
        """Separate unit-L1 dipoles, averaged equally; never a stronger K5 input."""
        degree = tuple(map(len, adjacency))
        pairs = tuple(combinations(core, 2))
        integrated = {h: [0.0]*n for h in horizons}
        exact_reached, trials = set(), []
        max_final_signal = 0.0
        for a, b in pairs:
            x = [0.0]*n
            x[a], x[b] = 0.5, -0.5
            exact = [Fraction(0)]*n
            exact[a], exact[b] = Fraction(1, 2), Fraction(-1, 2)
            integral, reached = [0.0]*n, set()
            previous_norm = 1.0
            for t in range(horizons[-1]+1):
                norm = fsum(abs(v) for v in x)
                assert norm <= previous_norm + 2e-13
                assert abs(fsum(x)) < 2e-13  # No net/common-mode injection.
                previous_norm = norm
                if t <= exact_steps:
                    assert max(abs(float(e)-v) for e, v in zip(exact, x)) < 2e-13
                    reached.update(v for v in registry if exact[v] != 0)
                for j, value in enumerate(x):
                    integral[j] += abs(value)
                if t in horizons:
                    for j, value in enumerate(integral):
                        integrated[t][j] += value / len(pairs)
                if t == horizons[-1]:
                    max_final_signal = max(max_final_signal, norm)
                    break
                # Source normalization: outgoing negative amplitude / source degree.
                x = [-fsum(x[j]/degree[j] for j in adjacency[i]) for i in range(n)]
                if t < exact_steps:
                    exact = [-sum((exact[j]/degree[j] for j in adjacency[i]), Fraction(0))
                             for i in range(n)]
            exact_reached.update(reached)
            trials.append(dict(source_pair=(a, b), reached_by_step_8=tuple(sorted(reached)),
                               W={k: fsum(v*c for v, c in zip(integral, cost))
                                  for k, cost in kernels.items()}))
        snapshots = {}
        for h, activity in integrated.items():
            registry_activity = [activity[v] for v in registry]
            total = fsum(registry_activity)
            squares = fsum(v*v for v in registry_activity)
            snapshots[h] = dict(
                W={k: fsum(v*c for v, c in zip(activity, cost)) for k, cost in kernels.items()},
                registry_activity=total,
                effective_registry_participation=total*total/squares if squares else None,
                integrated_activity=tuple(activity))
        return dict(snapshots=snapshots, trials=trials,
                    reached_by_step_8=tuple(sorted(exact_reached)),
                    max_final_signal=max_final_signal)

    responses = {}
    for row in catalogue:
        i = row['archive_motif_id']
        core = tuple(embed(v) for v in sorted(motifs[i]['core'], key=_term))
        support = {embed(v) for v in motifs[i]['support']}
        # Isolation disables ONLY signal exchange across the support boundary.
        # Native readout costs stay identical; no physical history is deleted.
        local = tuple(tuple(j for j in neighbors if j in support) if v in support else ()
                      for v, neighbors in enumerate(graph))
        responses[row['name']] = dict(archive_motif_id=i, full=evolve(core, graph),
                                      isolated=evolve(core, local))
    anchor = responses['K3-1']['full']['snapshots'][128]['W']['A']
    if anchor <= 0:
        raise ValueError('Declared response reference is zero; no replacement chosen.')
    summaries = []
    for row in catalogue:
        response = responses[row['name']]
        full, local = (response[k]['snapshots'][128] for k in ('full', 'isolated'))
        summaries.append(dict(
            name=row['name'], W_A=full['W']['A'], W_B=full['W']['B'],
            ratio_A=full['W']['A']/anchor, ratio_B=full['W']['B']/anchor,
            gain_A=full['W']['A']/local['W']['A'], gain_B=full['W']['B']/local['W']['B'],
            N_eff=full['effective_registry_participation'],
            reached=len(response['full']['reached_by_step_8']),
            relative_addition_64_to_128={k: (full['W'][k] - response['full']['snapshots'][64]['W'][k])
                                           /full['W'][k] for k in kernels}))

    print('CELL 6 — PROPOSED COUPLING EXTENSION; NATIVE RECORDING KERNEL UNCHANGED')
    print(f'Common host: {n} objects, {len(registry)} registry roots, {sum(map(len,graph))//2} parent-child links.')
    print('Added law: T=-A D^-1 on those links; |signal| is a contact probability.')
    print('One unit-L1 balanced pulse per trial; all core pairs averaged, not summed.')
    print('Each step permits at most one contact; leftover probability means no contact.')
    print('Expected incremental work is summed through step 128; apparatus warming excluded.')
    print('No service time, field energy, sector weights, or mass multipliers were supplied.')
    print('\n case      W(A)      W(B)     ratio A   ratio B   gain A   N_eff   reached<=8')
    for r in summaries:
        print(f" {r['name']:6s} {r['W_A']:9.4f} {r['W_B']:9.4f} {r['ratio_A']:10.6f}"
              f" {r['ratio_B']:9.6f} {r['gain_A']:8.4f} {r['N_eff']:7.2f} {r['reached']:10d}")
    print('ratio = full response / K3-1 full response with readout A (fixed once).')
    print('gain = full-host / isolated-support response with the SAME readout kernel.')
    print('N_eff uses time-integrated activity, NOT independent physical degrees of freedom.')
    print('Reached counts are exact UNION counts across trials, not counts in every single pulse.')
    maximum_change = max(v for r in summaries for v in r['relative_addition_64_to_128'].values())
    print(f'Largest full-response fractional addition from steps 65–128: {maximum_change:.3e}.')
    print('This is a finite-window check, not an infinite-time tail bound.')
    print('PASS: legal host, primitive-exchange symmetry, real warming, all 14 old A/B Q values,')
    print('      exact first-eight-step transport checks, zero signed sum, and nonincreasing input budget.')
    print('Transport is floating point; first 8 steps are also checked with exact fractions.')
    print('W is expected native work under a PROPOSED scheduler, not a measured rest mass.')
    print('Cell 5 and its physical calibration remain unchanged; no MeV conversion is made here.')
    return dict(protocol='declared parent-child signed transport + probabilistic native readout',
                horizons=horizons, exact_steps=exact_steps, host_parents=tuple(host.parents),
                registry_labels={rmap[r['object_id']]: r['label'] for r in ref['rows']},
                kernels=kernels, warming_work_excluded=warming, responses=responses,
                summaries=summaries, response_anchor=dict(name='K3-1', readout='A', W=anchor),
                physical_mass_calibration=None)


dcu_mass_6 = _run_dcu_mass_6()
