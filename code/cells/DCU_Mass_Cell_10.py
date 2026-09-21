# CELL 10 — witnessed-edge 3D completion with exact algebraic certificates.
# Run after Cells 8 and 9. Standard library only. Missing distances are FREE,
# not zero; witness presence is NOT the witness object's own recorded-mask bit.
from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import sqrt


def _run_dcu_mass_10():
    needed = ('dcu_mass_8', 'dcu_mass_9', '_reference', '_parents', '_Native')
    absent = [k for k in needed if k not in globals()]
    if absent:
        raise RuntimeError('Run Cells 8 and 9 first. Missing: ' + ', '.join(absent))
    previous = dcu_mass_9
    audit = dcu_mass_8['audit_distances']
    before = deepcopy(previous)
    host = _Native()
    for expected, pair in enumerate(previous['common_host_parents'][2:], 2):
        born, _ = host.add_batch([tuple(pair)])
        assert born == [expected]
    assert len(host) == 180
    rmap = {0: 0, 1: 1}
    for v, pair in enumerate(_reference['parents'][2:], 2):
        rmap[v] = host.pair_to_id[tuple(sorted(rmap[p] for p in pair))]
    S = tuple(rmap[r['object_id']] for r in _reference['rows'] if r['sector'] == 'S')

    @lru_cache(None)
    def embed(v):
        if v < 2:
            return v
        return host.pair_to_id[tuple(sorted(embed(p) for p in _parents[v]))]

    def inputs(case):
        core = tuple(case['cache_core'])
        d = [list(map(Fraction, row)) for row in case['readouts']['S_overlap']['D2']]
        features = [tuple(len(host.ancestors[embed(v)] & host.ancestors[r]) for r in S)
                    for v in core]
        expected = [[Fraction(sum((x-y)**2 for x, y in zip(a, b)), len(S))
                     for b in features] for a in features]
        assert d == expected  # Same 81-entry readout, independently reconstructed.
        witnesses = {}
        for i, j in combinations(range(len(core)), 2):
            pair = tuple(sorted((embed(core[i]), embed(core[j]))))
            if pair in host.pair_to_id:
                witnesses[i, j] = host.pair_to_id[pair]
        mask = [[Fraction(0) if i == j else
                 d[i][j] if tuple(sorted((i, j))) in witnesses else None
                 for j in range(len(core))] for i in range(len(core))]
        result = audit(d)
        if not result['positive_semidefinite']:
            raise ValueError('The fixed S-feature matrix is not Euclidean.')
        return core, d, witnesses, mask, result

    def completion(d, missing):
        """Two tetrahedra on a shared triangle: exact certificate for K5-e."""
        u, v = missing
        a, b, c = [j for j in range(5) if j not in missing]
        A, C = d[a][b], d[a][c]
        B = (A+C-d[b][c])/2
        det = A*C-B*B
        if A <= 0 or det <= 0:
            raise ValueError('Degenerate base triangle: this certificate does not apply.')
        H = ((A, B), (B, C))
        def product(x, y):
            return sum(x[i]*H[i][j]*y[j] for i in range(2) for j in range(2))
        coefficients = {a: (Fraction(0), Fraction(0)),
                        b: (Fraction(1), Fraction(0)),
                        c: (Fraction(0), Fraction(1))}
        heights = {a: Fraction(0), b: Fraction(0), c: Fraction(0)}
        for z in (u, v):
            g0, g1 = (d[a][z]+A-d[b][z])/2, (d[a][z]+C-d[c][z])/2
            p = ((C*g0-B*g1)/det, (A*g1-B*g0)/det)
            coefficients[z] = p
            heights[z] = d[a][z] - product(p, p)
            if heights[z] <= 0:
                raise ValueError('Degenerate apex: this certificate requires positive height.')
        # These exact identities certify ALL nine retained distances in R3.
        # Only the apex-to-apex distance depends on the relative height sign.
        for i, j in combinations(range(5), 2):
            delta = tuple(x-y for x, y in zip(coefficients[i], coefficients[j]))
            if (i, j) != missing:
                assert product(delta, delta) + heights[i] + heights[j] == d[i][j]
        delta = tuple(x-y for x, y in zip(coefficients[u], coefficients[v]))
        center = product(delta, delta) + heights[u] + heights[v]
        radicand = heights[u]*heights[v]
        polynomial = (d[u][v]-center)**2 - 4*radicand
        # Strictly negative means the original length lies between, not at,
        # the two allowed squared distances center +/- 2*sqrt(radicand).
        if polynomial >= 0:
            raise ValueError('Expected the unchanged full five-point metric to have rank four.')

        options = []
        for sign in (-1, 1):
            # Numeric coordinates are illustrations/checks; feasibility above is exact.
            coordinates = []
            for z in range(5):
                p, q = coefficients[z]
                x = float(p)*sqrt(float(A)) + float(q*B)/sqrt(float(A))
                y = float(q)*sqrt(float(det/A))
                h = sqrt(float(heights[z])) * (sign if z == v else 1)
                coordinates.append((x, y, h))
            actual = [[sum((x-y)**2 for x, y in zip(p, q)) for q in coordinates]
                      for p in coordinates]
            residual = max(abs(actual[i][j]-float(d[i][j]))
                           for i, j in combinations(range(5), 2) if (i, j) != missing)
            value = float(center)-2*sign*sqrt(float(radicand))
            assert residual < 1e-10 * max(1.0, max(map(float, sum(d, []))))
            assert abs(actual[u][v]-value) < 1e-10
            options.append(dict(relative_height_sign=sign, missing_squared_distance=value,
                                coordinates_R3=coordinates, max_nine_edge_error=residual))
        return dict(missing_indices=missing, base_indices=(a, b, c), plane_Gram=H,
                    plane_coefficients=coefficients, squared_heights=heights,
                    center=center, radicand=radicand,
                    original_squared_distance=d[u][v],
                    polynomial_at_original=polynomial, options=options,
                    exact_nine_edge_certificate=True)

    candidates, releases = {}, {}
    for name, case in previous['catalogue'].items():
        core, d, witnesses, mask, result = inputs(case)
        assert len(witnesses) == len(core)*(len(core)-1)//2
        candidates[name] = dict(core=core, witnesses=witnesses, partial_D2=mask,
                                source_rank=result['rank'], fits_R3=result['fits_R3'],
                                full_metric_certificate=result)
        if len(core) == 5:
            assert result['rank'] == 4
            # Mathematical constraint release ONLY: no native object is deleted.
            releases[name] = {edge: completion(d, edge) for edge in witnesses}

    controls, native_closures = [], []
    for case in previous['controls']:
        core, d, witnesses, mask, result = inputs(case)
        missing = tuple(e for e in combinations(range(5), 2) if e not in witnesses)
        assert tuple((core[i], core[j]) for i, j in missing) == tuple(case['missing_edges'])
        assert result['rank'] == 4 and missing
        # Every absent-edge choice works. Other absent distances may be assigned
        # their source values as ONE feasible completion, not as constraints.
        certificates = [completion(d, edge) for edge in missing]
        controls.append(dict(core=core, witnesses=witnesses, partial_D2=mask,
                             fits_R3=True, certificates=certificates))
        if len(missing) == 1:
            # This is also a LEGAL native before/after experiment, not a deletion.
            state = deepcopy(host)
            i, j = missing[0]
            pair = tuple(sorted((embed(core[i]), embed(core[j]))))
            old_ancestries = [frozenset(state.ancestors[embed(v)]) for v in core]
            born, work = state.add_batch([pair])
            assert pair not in host.pair_to_id and state.parents[born[0]] == pair
            assert old_ancestries == [frozenset(state.ancestors[embed(v)]) for v in core]
            assert born[0] not in state.recorded  # PRESENT witness, still a fresh root.
            assert all(tuple(sorted((embed(core[p]), embed(core[q])))) in state.pair_to_id
                       for p, q in combinations(range(5), 2))
            native_closures.append(dict(core=core, added_cache_pair=(core[i], core[j]),
                                       new_host_child=born[0], native_work=work,
                                       before_fits_R3=True, after_fits_R3=False))

    assert previous == before  # Earlier results and graph were not changed.
    all_certificates = [c for r in releases.values() for c in r.values()] + [
        c for r in controls for c in r['certificates']]
    print('CELL 10 — WITNESSED-EDGE COMPLETION; S DISTANCES AND HOST UNCHANGED')
    print('None means unconstrained, not distance zero. Witness = actual child in host.')
    print('An existing witness need NOT itself belong to the native recorded mask.')
    print('\n case    constraints   original Gram rank   witnessed distances fit R3')
    for name, r in candidates.items():
        print(f" {name:6s} {len(r['witnesses']):11d} {r['source_rank']:20d} "
              f"{str(r['fits_R3']):>29s}")
    print(f'\nActual non-K5 controls: {len(controls)}/{len(controls)} have an exact 3D completion.')
    print(f"K5 single-constraint releases: {sum(map(len, releases.values()))}/30 are feasible.")
    print(f'Legal final-edge constructions: {len(native_closures)} local controls go feasible -> infeasible.')
    distinct_pairs = sorted({r['added_cache_pair'] for r in native_closures})
    print(f'  These are overlapping control cores, using only {len(distinct_pairs)} distinct added pairs: {distinct_pairs}.')
    print('  Native construction still accepts them; geometry is not a formation gate.')

    first_control = controls[0]  # First saved control; no gap-based selection.
    first = first_control['certificates'][0]
    u, v = first['missing_indices']
    print(f"\nFirst actual sparse control: {first_control['core']}; "
          f"absent pair {(first_control['core'][u], first_control['core'][v])}.")
    print(f"  Native squared distance when closure is required = {first['original_squared_distance']}.")
    print(f"  Allowed 3D values = {first['center']} +/- 2*sqrt({first['radicand']}).")
    print('  Numeric values:', sorted(o['missing_squared_distance'] for o in first['options']))
    print(f"  Exact closure polynomial at native value = {first['polynomial_at_original']} < 0.")
    max_error = max(o['max_nine_edge_error'] for c in all_certificates for o in c['options'])
    print(f'\nPASS: {len(all_certificates)} exact tetrahedral certificates; '
          f'{2*len(all_certificates)} explicit 3D realizations checked.')
    print(f'Largest numeric nine-edge squared-distance residual: {max_error:.3e}.')
    print('Exact arithmetic establishes feasibility/infeasibility; no optimizer was used.')
    print('Local realizations may differ between cases; no common global embedding is claimed.')
    print('Witness-to-preferred-length and dimension 3 remain stated premises.')
    print('Incompatibility is NOT yet a stress, energy, workload surcharge, or mass.')
    return dict(protocol='fixed S-overlap metric, actual witnessed-edge masks, local R3 completion',
                candidates=candidates, controls=controls, K5_constraint_releases=releases,
                legal_closures=native_closures, completion_certificate=completion,
                max_coordinate_residual=max_error, physical_mass_calibration=None)


dcu_mass_10 = _run_dcu_mass_10()
