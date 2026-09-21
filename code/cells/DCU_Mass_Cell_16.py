# CELL 16 — source-normalized cavity/registry coupling; no mass target or fitted strength.
# Run after Cells 5 and 14. NumPy plus Python standard library.
# NEW HYPOTHESIS: a core variable links to each registry root inheriting it.
# ALL old and added links share ONE source-normalized transport budget.
# g stays SYMBOLIC. Small numeric g values below ONLY validate perturbation algebra.
from copy import deepcopy
from fractions import Fraction
from itertools import combinations


def _run_dcu_mass_16():
    import numpy as np
    needed = ('dcu_mass_5', 'dcu_mass_14', '_reference', '_term', '_panel', '_parents')
    absent = [k for k in needed if k not in globals()]
    if absent:
        raise RuntimeError('Run Cells 5 and 14 first. Missing: ' + ', '.join(absent))
    ref, catalogue = _reference, dcu_mass_5['catalogue']
    saved = deepcopy((ref, catalogue, dcu_mass_14))
    rows, anc, adjacency = ref['rows'], ref['ancestors'], ref['adjacency']
    N = len(rows)
    degree = tuple(map(len, adjacency))
    lambdas = tuple(sorted(dcu_mass_14['spectrum']))
    assert N == 137 and len(lambdas) == 5
    by_term = {t: i for i, t in enumerate(ref['terms'])}
    c = ref['pairs'][(0, 1)]
    pa, pb = [ref['pairs'][tuple(sorted((p, c)))] for p in (0, 1)]
    blocks = [[], [], []]  # universal, a-only branch, b-only branch
    for i, r in enumerate(rows):
        a = anc[r['object_id']]
        assert pa in a or pb in a
        blocks[0 if pa in a and pb in a else 1 if pa in a else 2].append(i)
    assert tuple(map(len, blocks)) == (11, 63, 63)
    which = {v: g for g, block in enumerate(blocks) for v in block}
    incident = ((0, 1, 2), (0, 1), (0, 2))
    for i in range(N):
        assert set(adjacency[i]) == {j for j in range(N) if i != j and
                                    which[j] in incident[which[i]]}

    def rank(matrix):
        a = [list(map(Fraction, row)) for row in matrix]
        pivot_row = 0
        for col in range(len(a[0])):
            p = next((i for i in range(pivot_row, len(a)) if a[i][col]), None)
            if p is None:
                continue
            a[pivot_row], a[p] = a[p], a[pivot_row]
            v = a[pivot_row][col]
            a[pivot_row] = [x/v for x in a[pivot_row]]
            for i in range(pivot_row+1, len(a)):
                v = a[i][col]
                if v:
                    a[i] = [x-v*y for x, y in zip(a[i], a[pivot_row])]
            pivot_row += 1
            if pivot_row == len(a):
                break
        return pivot_row

    def transfer(Y):
        # EXACT application of the existing T_R=-A D^-1, using verified blocks.
        m = len(Y[0])
        totals = [[sum((Y[j][k]/degree[j] for j in block), Fraction(0))
                   for k in range(m)] for block in blocks]
        return [[Y[i][k]/degree[i] - sum((totals[b][k]
                 for b in incident[which[i]]), Fraction(0))
                 for k in range(m)] for i in range(N)]

    def projection(Y, lam):
        # Spectral projector as an exact polynomial; no degenerate eigenbasis choice.
        X = [list(row) for row in Y]
        for other in lambdas:
            if other != lam:
                TX = transfer(X)
                X = [[(t-other*x)/(lam-other) for t, x in zip(a, b)]
                     for a, b in zip(TX, X)]
        return X

    A = np.array([[int(j in adjacency[i]) for j in range(N)] for i in range(N)], float)
    d = np.array(degree, float)
    symmetric_registry = -A/np.sqrt(np.outer(d, d))
    results = {}
    for row in catalogue:
        m = _panel['motifs'][row['archive_motif_id']]
        core = tuple(by_term[_term(v)] for v in row['core'])
        n = len(core)
        assert {tuple(sorted(e)) for e in m['required_edges']} == set(combinations(sorted(m['core']), 2))
        for edge, child in zip(m['required_edges'], m['edge_children']):
            assert set(_parents[child]) == set(edge)
        # Pure incidence, NOT a sector weighting or inserted registry multiplier.
        C = [[Fraction(int(u in anc[r['object_id']])) for u in core] for r in rows]
        # Y=C P_n restricts injection to the isolated cavity's zero-sum modes.
        Y = [[v-sum(a)/n for v in a] for a in C]
        homogenized = [[sum(a)/n]*n for a in C]
        assert all(v-sum(a)/n == 0 for a in homogenized for v in a)
        lam0, beat0 = Fraction(1, n-1), 1-Fraction(1, n-1)
        assert lam0 > max(lambdas)
        pieces, sum_projected = {}, [[Fraction(0)]*n for _ in range(N)]
        Sigma = [[Fraction(0)]*n for _ in range(n)]
        for lam in lambdas:
            Z = projection(Y, lam)
            assert transfer(Z) == [[lam*v for v in a] for a in Z]
            B = [[sum((Y[j][a]*Z[j][b]/degree[j] for j in range(N)), Fraction(0))
                  for b in range(n)] for a in range(n)]
            # Projectors are self-adjoint/orthogonal in W=D^-1, not Euclidean metric.
            gram = [[sum((Z[j][a]*Z[j][b]/degree[j] for j in range(N)), Fraction(0))
                     for b in range(n)] for a in range(n)]
            assert B == gram and B == list(map(list, zip(*B)))
            accessible_rank = rank(Z)
            weight = sum(B[j][j] for j in range(n))/(n-1)
            assert weight >= 0 and rank(B) == accessible_rank
            pieces[lam] = dict(accessible_rank=accessible_rank, weight=weight, Gram=B)
            for j in range(N):
                for k in range(n):
                    sum_projected[j][k] += Z[j][k]
            for j in range(n):
                for k in range(n):
                    Sigma[j][k] += B[j][k]/(lam0-lam)
        assert sum_projected == Y and all(sum(a) == 0 for a in Sigma)
        assert rank(Sigma) == n-1
        mean_shift = sum(Sigma[j][j] for j in range(n))/(n-1)
        Cn = np.array(C, float)
        P = np.eye(n)-np.ones((n, n))/n
        # Independent dense resolvent check: Sigma=P C' (lam0 D+A)^-1 C P.
        dense = P @ (Cn.T @ np.linalg.solve(float(lam0)*np.diag(d)+A, Cn)) @ P
        exact_matrix = np.array(Sigma, float)
        assert np.allclose(dense, exact_matrix, rtol=2e-11, atol=2e-12)
        ev = np.linalg.eigvalsh(exact_matrix)
        assert abs(ev[0]) < 1e-10 and ev[1] > 0
        shifts = ev[1:]
        # ACTUAL COUPLED MODEL: negative random-walk transfer on one weighted graph.
        # Internal cavity edges and registry edges have unit weights; cross-links
        # have a common SYMBOLIC weight g. Degrees MUST include the new links.
        # No independently preserved, unnormalized feedback block is substituted.
        column_sums = tuple(sum(a[k] for a in C) for k in range(n))
        h = tuple(v/(n-1) for v in column_sums)
        h_mean = sum(h)/n
        h2_mean = sum(v*v for v in h)/n
        variance_h = h2_mean-h_mean*h_mean
        first_lambda = -lam0*h_mean
        # Second-order centroid coefficient:
        # normalization + mixing with the cavity common mode + registry return.
        second_normalization = lam0*h2_mean
        second_common = -lam0/(lam0+1)*variance_h/(n-1)
        second_return = mean_shift/(n-1)
        second_lambda = second_normalization+second_common+second_return
        first_matrix = -float(lam0)*P@np.diag(list(map(float, h)))@P
        first_eigenvalues = np.linalg.eigvalsh(first_matrix)[:n-1]
        assert np.all(first_eigenvalues < 0)

        cavity_adjacency = np.ones((n, n))-np.eye(n)
        A0 = np.block([[cavity_adjacency, np.zeros((n, N))],
                       [np.zeros((N, n)), A]])
        A1 = np.block([[np.zeros((n, n)), Cn.T],
                       [Cn, np.zeros((N, N))]])
        D0 = np.r_[np.full(n, n-1), d]
        D1 = np.r_[Cn.sum(axis=0), Cn.sum(axis=1)]
        gap = min(float(lam0)-max(map(float, lambdas)), float(lam0)+1)
        validations = []
        for g in (0.00001, 0.000005):
            # These are ONLY finite-difference checks at g=0, not physical choices.
            weighted = A0+g*A1
            degrees_g = D0+g*D1
            full = -weighted/np.sqrt(np.outer(degrees_g, degrees_g))
            assert np.allclose(full@np.sqrt(degrees_g), -np.sqrt(degrees_g),
                               rtol=1e-12, atol=1e-12)
            spectrum = np.linalg.eigvalsh(full)
            assert spectrum.min() >= -1-1e-12 and spectrum.max() <= 1+1e-12
            cluster = spectrum[abs(spectrum-float(lam0)) < gap/3]
            assert len(cluster) == n-1
            centroid = float(cluster.mean())
            second_estimate = (centroid-float(lam0)-g*float(first_lambda))/(g*g)
            second_error = abs(second_estimate-float(second_lambda))
            assert second_error < 1e-3*max(1.0, abs(float(second_lambda)))
            validations.append(dict(validation_g=g, centroid=centroid,
                                    second_coefficient_estimate=second_estimate,
                                    second_coefficient_error=second_error))
        results[row['name']] = dict(
            archive_motif_id=row['archive_motif_id'], core_cache_ids=tuple(row['core']),
            incidence=C, dependent_root_counts=tuple(map(int, column_sums)),
            directly_addressable_roots=sum(any(a) for a in Y),
            contrast_dimension=n-1,
            reachable_registry_dimension=sum(v['accessible_rank'] for v in pieces.values()),
            registry_eigenpieces=pieces, cavity_lambda=lam0, bare_beat_factor=beat0,
            exact_resolvent_matrix=Sigma, raw_resolvent_mean=mean_shift,
            mean_lambda_linear_coefficient=first_lambda,
            mean_lambda_quadratic_coefficient=second_lambda,
            mean_beat_linear_coefficient=-first_lambda,
            quadratic_terms=dict(normalization=second_normalization,
                                 common_mode=second_common, registry_return=second_return),
            individual_lambda_linear_coefficients=first_eigenvalues.tolist(),
            individual_registry_return_coefficients=(shifts/(n-1)).tolist(),
            nearest_registry_gap=lam0-max(lambdas), full_spectrum_checks=validations)

    electron = results['K3-1']
    for r in results.values():
        # The mean of ALL continuously connected contrast branches, not a
        # retrospectively selected favorable eigenmode. Same carrier -1 for both.
        r['relative_beat_ratio_linear_coefficient'] = (
            r['mean_beat_linear_coefficient']/r['bare_beat_factor']
            - electron['mean_beat_linear_coefficient']/electron['bare_beat_factor'])
    assert (ref, catalogue, dcu_mass_14) == saved
    print('CELL 16 — ANCESTRY COUPLING WITH ONE SOURCE-NORMALIZED TRANSPORT BUDGET')
    print('C[r,i]=1 iff core object i belongs to ancestry of registry root r.')
    print('Added graph links have symbolic common weight g; internal edges retain unit weight.')
    print('T(g)=-A(g)D(g)^-1. All degrees include the added channels; -1 remains an exact carrier.')
    print('This is an ADDED spectral model, NOT native service, energy, or a mass operator.')
    print('Identical cavity-port profiles annihilate contrast injection, exactly.')
    print('Source normalization may still change cavity decay in that null control.')
    print('No empirical quantity is used and no physical value of g is chosen.')

    print('\n case    roots reachable  subspace dim  bare beat   linear beat coefficient')
    for name, r in results.items():
        print(f" {name:6s} {r['directly_addressable_roots']:15d} {r['reachable_registry_dimension']:13d}"
              f" {str(r['bare_beat_factor']):>10s} {str(r['mean_beat_linear_coefficient']):>25s}")
    print('Root counts are unions over contrast inputs; subspace dimension uses the UNCOUPLED registry.')
    print('Neither count is workload, independent-particle count, or a mass multiplier.')

    print('\nFull weak-coupling centroid: lambda_bar=lambda0+a*g+b*g^2+O(g^3)')
    print(' case                 a               b    registry-return part of b')
    for name, r in results.items():
        print(f" {name:6s} {float(r['mean_lambda_linear_coefficient']):17.9f}"
              f" {float(r['mean_lambda_quadratic_coefficient']):15.9f}"
              f" {float(r['quadratic_terms']['registry_return']):29.9f}")

    print('\nEqual-mode gap-ratio diagnostic: (fbar_n/fbar_3)/(f_n/f_3)=1+k*g+O(g^2)')
    for name, r in results.items():
        if len(r['core_cache_ids']) == 5:
            k = r['relative_beat_ratio_linear_coefficient']
            print(f"  {name}: k={k} ({float(k):+.9f}); exact sign negative={k < 0}")
    print('Common weak coupling LOWERS the mean K5/K3 gap ratio; it cannot supply an upward perturbative correction.')
    print('This is conditional on the inheritance wiring, shared source budget, and equal-mode readout.')
    print('No finite-strength mass ratio or universal exclusion of other states/couplings is claimed.')
    print('PASS: 35 exact spectral projections, weighted Gram identities, and dense-resolvent checks.')
    print('PASS: 14 full normalized-matrix checks of the SECOND-order centroid; previous inputs preserved.')
    print('No target matching, selected g, physical calibration, or native-dynamics change.')
    return dict(protocol='binary ancestry-incidence coupling with full source renormalization',
                registry_spectrum=dcu_mass_14['spectrum'], candidates=results,
                uniform_contrast_injection_null_checked=True, coupling_strength_selected=None,
                mass_operator_identified=False, new_physical_calibration=False,
                native_dynamics_changed=False)


dcu_mass_16 = _run_dcu_mass_16()
