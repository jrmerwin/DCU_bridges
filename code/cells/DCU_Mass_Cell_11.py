# CELL 11 — quantify unavoidable witnessed-distance mismatch, not energy or mass.
# Run after Cell 10. Requires NumPy and SciPy in this notebook's environment.
# All witnessed constraints have equal weight. Missing entries are never zero-filled.
# Optimizers supply feasible UPPER bounds; exact rational algebra supplies LOWER bounds.
from copy import deepcopy
from fractions import Fraction
from itertools import combinations
from math import sqrt


def _run_dcu_mass_11():
    try:
        import numpy as np
        from scipy.optimize import least_squares
    except ImportError as exc:
        raise ImportError('Cell 11 needs NumPy and SciPy in the notebook kernel.') from exc
    needed = ('dcu_mass_8', 'dcu_mass_9', 'dcu_mass_10')
    missing = [k for k in needed if k not in globals()]
    if missing:
        raise RuntimeError('Run Cells 8–10 first. Missing: ' + ', '.join(missing))
    old9, old10 = dcu_mass_9, dcu_mass_10
    before = (deepcopy(old9), deepcopy(old10))
    audit = dcu_mass_8['audit_distances']
    complete = old10['completion_certificate']

    def positive_definite(matrix):
        # Exact LDL/Schur-complement pivots: no floating eigenvalue tolerance.
        a = [list(map(Fraction, row)) for row in matrix]
        for k in range(len(a)):
            pivot = a[k][k]
            if pivot <= 0:
                return False
            for i in range(k + 1, len(a)):
                for j in range(i, len(a)):
                    value = a[i][j] - a[i][k]*a[j][k]/pivot
                    a[i][j] = a[j][i] = value
        return True

    def lower_certificate(d, report):
        # For n=5, B has four positive eigenvalues on the centered subspace.
        # F=[e0-e4,...,e3-e4]: F'BF=A and F'F=I+11'.
        n = len(d)
        A = [[(d[i][-1] + d[j][-1] - d[i][j])/2
              for j in range(n-1)] for i in range(n-1)]
        def shifted(t):
            return [[A[i][j] - t*(1 + int(i == j))
                     for j in range(n-1)] for i in range(n-1)]
        lo = Fraction(0)
        hi = sum(report['gram'][i][i] for i in range(n))/(n-1)
        assert positive_definite(shifted(lo)) and not positive_definite(shifted(hi))
        for _ in range(42):
            mid = (lo+hi)/2
            if positive_definite(shifted(mid)):
                lo = mid
            else:
                hi = mid
        assert lo > 0 and positive_definite(shifted(lo))
        # Every centered rank<=3 realization has a centered unit null vector v.
        # lo <= v'Bv = sum_{i<j} e_ij v_i v_j and
        # sum_{i<j} v_i^2 v_j^2 <= (n-1)/(2n). Hence SSE >= 2n/(n-1)*lo^2.
        return dict(lambda_lower=lo, lambda_upper=hi,
                    exact_positive_definite_matrix=shifted(lo),
                    SSE_lower=Fraction(2*n, n-1)*lo*lo)

    def rational_realization(x, d, edges):
        # Rounding gives explicit rational 3D coordinates, then an EXACT upper bound.
        x = np.asarray(x, dtype=float)
        x = x - x.mean(axis=0)
        q = tuple(tuple(Fraction(f'{v:.12f}') for v in row) for row in x)
        residuals = {}
        for i, j in edges:
            residuals[i, j] = sum((a-b)**2 for a, b in zip(q[i], q[j])) - d[i][j]
        return dict(coordinates=q, residuals=residuals,
                    SSE=sum((r*r for r in residuals.values()), Fraction(0)))

    def fit_complete(squared):
        # Fixed objective: sum squared errors of S-derived SQUARED distances.
        d = [list(map(Fraction, row)) for row in squared]
        report = audit(d)
        if not report['positive_semidefinite']:
            raise ValueError('Expected the unchanged Euclidean S-overlap readout.')
        n = len(d)
        edges = tuple(combinations(range(n), 2))
        norm2 = sum(d[i][j]**2 for i, j in edges)
        if norm2 <= 0:
            raise ValueError('Zero normalization; no substitute scale selected.')
        if report['rank'] <= 3:
            return dict(eta2_lower=Fraction(0), eta2_upper=Fraction(0),
                        status='EXACT_ZERO', source_rank=report['rank'],
                        source_D2=d, source_certificate=report)
        if n != 5 or report['rank'] != 4:
            raise ValueError('This bound/initialization protocol is for full-rank five-point data.')
        bound = lower_certificate(d, report)
        B = np.array(report['gram'], dtype=float)
        ev, vec = np.linalg.eigh(B)
        axes = tuple(range(n-4, n))
        if any(ev[a] <= 0 for a in axes):
            raise ArithmeticError('Numerical eigendecomposition disagrees with exact rank.')
        # Four spectral starts, plus both tetrahedral completions for EACH edge.
        # Completions are ONLY initial guesses: all TEN residuals remain in the fit.
        starts = [vec[:, a]*np.sqrt(ev[list(a)]) for a in combinations(axes, 3)]
        certs = {e: complete(d, e) for e in edges}
        starts.extend(np.array(o['coordinates_R3'])
                      for c in certs.values() for o in c['options'])
        assert len(starts) == 24
        ii, jj = np.array(edges).T
        target = np.array([float(d[i][j]) for i, j in edges])
        scale = sqrt(float(norm2))
        def residual(y):
            x = y.reshape(n, 3)
            dx = x[ii] - x[jj]
            return ((dx*dx).sum(axis=1)-target)/scale
        def jacobian(y):
            x = y.reshape(n, 3)
            a = 2*(x[ii]-x[jj])/scale
            out = np.zeros((len(edges), n, 3))
            out[np.arange(len(edges)), ii] = a
            out[np.arange(len(edges)), jj] = -a
            return out.reshape(len(edges), n*3)
        trials, best = [], None
        for start_id, x0 in enumerate(starts):
            # Retain each starting realization too: failure cannot erase a valid bound.
            start = rational_realization(x0, d, edges)
            fit = least_squares(residual, x0.ravel(), jac=jacobian, method='trf',
                                loss='linear', x_scale=1.0, ftol=1e-13, xtol=1e-13,
                                gtol=1e-13, max_nfev=5000)
            if not np.all(np.isfinite(fit.x)):
                raise ArithmeticError('Nonfinite solver result; no result substituted.')
            result = rational_realization(fit.x.reshape(n, 3), d, edges)
            choice = min((start, result), key=lambda r: r['SSE'])
            if best is None or choice['SSE'] < best['SSE']:
                best = dict(choice, start_id=start_id)
            trials.append(dict(start_id=start_id, success=bool(fit.success),
                               status=int(fit.status), evaluations=int(fit.nfev),
                               optimality=float(fit.optimality),
                               eta2=float(result['SSE']/norm2), message=str(fit.message)))
        assert bound['SSE_lower'] <= best['SSE']
        # Dimensional countercontrol: all source distances fit in FOUR dimensions.
        x4 = vec[:, axes]*np.sqrt(ev[list(axes)])
        r4 = ((x4[ii]-x4[jj])**2).sum(axis=1)-target
        assert np.linalg.norm(r4)/scale < 1e-10
        return dict(status='POSITIVE_BRACKET_NOT_EXACT_MINIMUM', source_D2=d,
                    source_rank=4, normalization=norm2, lower_certificate=bound,
                    eta2_lower=bound['SSE_lower']/norm2, eta2_upper=best['SSE']/norm2,
                    best=best, trials=trials, all_edges_fitted=edges,
                    four_dimensional_relative_residual=float(np.linalg.norm(r4)/scale))

    def verify_zero(control):
        # Recheck Cell 10's exact 3D construction on EVERY currently witnessed edge.
        # At most one apex can have nonzero height in each such pair.
        c = control['certificates'][0]
        H, p, h = c['plane_Gram'], c['plane_coefficients'], c['squared_heights']
        u, v = c['missing_indices']
        assert h[u] > 0 and h[v] > 0
        checked = 0
        for i, j in control['witnesses']:
            assert (i, j) != (u, v) and h[i]*h[j] == 0
            z = tuple(a-b for a, b in zip(p[i], p[j]))
            actual = sum(z[a]*H[a][b]*z[b] for a in range(2) for b in range(2)) + h[i]+h[j]
            assert actual == control['partial_D2'][i][j]
            checked += 1
        return dict(core=control['core'], eta2_lower=Fraction(0), eta2_upper=Fraction(0),
                    status='EXACT_ZERO', checked_edges=checked)

    candidates = {name: fit_complete(row['readouts']['S_overlap']['D2'])
                  for name, row in old9['catalogue'].items()}
    controls = [verify_zero(c) for c in old10['controls']]
    source_controls = {tuple(c['cache_core']): c for c in old9['controls']}
    closures = []
    for event in old10['legal_closures']:
        source = source_controls[tuple(event['core'])]
        result = fit_complete(source['readouts']['S_overlap']['D2'])
        assert result['eta2_lower'] > 0
        closures.append(dict(core=event['core'], added_pair=event['added_cache_pair'],
                             native_formation_work=event['native_work'],
                             before_eta2=Fraction(0), after=result))
    assert (old9, old10) == before
    positive = [r for r in candidates.values() if r['source_rank'] == 4] + [r['after'] for r in closures]
    failed = sum(not t['success'] for r in positive for t in r['trials'])
    print('CELL 11 — MINIMUM RELATIVE RMS SQUARED-DISTANCE MISMATCH IN R3')
    print('eta^2 = min sum_witnesses (realized squared distance - S target)^2 / sum_witnesses target^2')
    print('Uniform edge weights. All coordinates free. No deleted or zero-filled constraints.')
    print('Exact rational lower/upper certificates; upper = a feasible fit, NOT a proved optimum.')
    print('\n case       lower eta (%)     upper eta (%)     status')
    for name, r in candidates.items():
        print(f" {name:7s} {100*sqrt(float(r['eta2_lower'])):16.8f}"
              f" {100*sqrt(float(r['eta2_upper'])):17.8f}     {r['status']}")
    print(f'\nSparse controls: {len(controls)}/{len(controls)} have EXACT minimum zero.')
    print(f"Previous legal final-edge cases: {len(closures)}/{len(closures)} now have a positive lower bound.")
    print('These are local, overlapping tests; no common global embedding is asserted.')
    print(f'{len(positive)*24} deterministic optimization starts; {failed} non-success terminations retained.')
    first = closures[0]
    print(f"\nFirst earlier closure: core {first['core']}, new pair {first['added_pair']}.")
    print(f"  Before eta=0; after eta in [{100*sqrt(float(first['after']['eta2_lower'])):.8f}%, "
          f"{100*sqrt(float(first['after']['eta2_upper'])):.8f}%].")
    print('  Native formation work is retained separately, NOT increased by this geometric diagnostic.')
    print('PASS: exact positive spectral bounds, rational 3D upper bounds, four-dimensional zero controls,')
    print('      all sparse witnessed-edge certificates, and preservation of Cells 9/10.')
    print('Percentages are relative RMS errors of SQUARED distances, not physical strain percentages.')
    print('No force, turnover, clock, energy, processing surcharge, or mass conversion was introduced.')
    return dict(protocol='uniform witnessed-edge relative squared-distance least squares in R3',
                candidates=candidates, controls=controls, legal_closure_magnitudes=closures,
                fit_complete=fit_complete, non_success_terminations=failed,
                exact_global_minimum_claimed=False, physical_mass_calibration=None)


dcu_mass_11 = _run_dcu_mass_11()
