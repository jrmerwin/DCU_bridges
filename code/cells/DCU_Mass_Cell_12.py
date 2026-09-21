# CELL 12 — quadratic-response feasibility and local spectral diagnostic.
# Run after Cells 5 and 11. NumPy only; NO new optimization, native dynamics,
# fitted coefficient, workload surcharge, or mass assignment is introduced.
from copy import deepcopy
from fractions import Fraction
from itertools import combinations


def _run_dcu_mass_12():
    import numpy as np
    needed = ('dcu_mass_5', 'dcu_mass_11')
    absent = [k for k in needed if k not in globals()]
    if absent:
        raise RuntimeError('Run Cells 5 and 11 first. Missing: ' + ', '.join(absent))
    baseline, geometry = dcu_mass_5, dcu_mass_11
    saved = (deepcopy(baseline['summaries']), deepcopy(geometry['candidates']))
    target = baseline['source_constants']['muon_electron_ratio']
    qrows = {r['name']: r for r in baseline['summaries']}
    assert baseline['anchor']['Q'] == Fraction(248, 9)
    assert geometry['candidates']['K3-1']['eta2_upper'] == 0

    # q_eff = Q/Q_e + g*eta^2. Invert for REQUIRED g only as a diagnostic.
    # These intervals are NOT used as model parameters or applied to any Q.
    sizing = []
    for name, geom in geometry['candidates'].items():
        lo, hi = geom['eta2_lower'], geom['eta2_upper']
        assert 0 <= lo <= hi
        if hi == 0:
            continue
        if lo <= 0:
            raise ValueError('Positive lower bound required for finite coefficient intervals.')
        for label in ('A', 'B'):
            q = qrows[name]['ratio_' + label]
            gap = target - q
            if gap <= 0:
                raise ValueError('This diagnostic assumes the target exceeds the baseline.')
            sizing.append(dict(name=name, probe=label, q=q,
                               eta2_bounds=(lo, hi), unit_g_ratio_bounds=(q+lo, q+hi),
                               required_g_bounds=(gap/hi, gap/lo)))

    def numeric_rank(singular, rows, columns):
        if not len(singular) or singular[0] == 0:
            return 0
        tolerance = max(rows, columns)*np.finfo(float).eps*singular[0]*100
        return int(np.sum(singular > tolerance))

    spectral = {}
    for name, geom in geometry['candidates'].items():
        d = np.asarray(geom['source_D2'], dtype=float)
        n = len(d)
        edges = tuple(combinations(range(n), 2))
        norm2 = sum(d[i, j]**2 for i, j in edges)
        if 'best' in geom:
            x = np.asarray(geom['best']['coordinates'], dtype=float)
        else:
            B = np.asarray(geom['source_certificate']['gram'], dtype=float)
            values, vectors = np.linalg.eigh(B)
            assert values.min() > -1e-10
            x = vectors[:, -3:]*np.sqrt(np.maximum(values[-3:], 0))
        x = x - x.mean(axis=0)

        def response_at(points):
            # Phi = (1/2) sum_e residual_e^2 / sum_e target_e^2.
            # This is the SAME normalized geometric objective, not native work.
            residual = np.empty(len(edges))
            R = np.zeros((len(edges), 3*n))
            omega = np.zeros((n, n))
            for k, (i, j) in enumerate(edges):
                delta = points[i]-points[j]
                residual[k] = delta@delta-d[i, j]
                R[k, 3*i:3*i+3] = 2*delta
                R[k, 3*j:3*j+3] = -2*delta
                r = residual[k]
                omega[i, i] += r
                omega[j, j] += r
                omega[i, j] -= r
                omega[j, i] -= r
            gradient = R.T@residual/norm2
            material = R.T@R/norm2
            prestress = 2*np.kron(omega, np.eye(3))/norm2
            return residual, R, gradient, material, prestress

        residual, R, grad, H0, Hp = response_at(x)
        H = H0+Hp
        left, singular, _ = np.linalg.svd(R, full_matrices=True)
        rrank = numeric_rank(singular, *R.shape)
        stress_basis = left[:, rrank:]
        stress_projection = stress_basis@(stress_basis.T@residual)
        nonzero_mismatch = geom['eta2_lower'] > 0
        relative_unbalanced = (np.linalg.norm(residual-stress_projection)/np.linalg.norm(residual)
                               if nonzero_mismatch else None)
        # Remove rigid translations and rotations, not an arbitrarily selected mode.
        rigid = np.column_stack(
            [np.tile(a, (n, 1)).ravel() for a in np.eye(3)] +
            [np.cross(np.tile(a, (n, 1)), x).ravel() for a in np.eye(3)])
        U, srigid, _ = np.linalg.svd(rigid, full_matrices=True)
        rigid_rank = numeric_rank(srigid, *rigid.shape)
        assert rigid_rank == 6
        internal = U[:, rigid_rank:]
        eigenvalues = np.linalg.eigvalsh(internal.T@H@internal)
        # Independent directional finite-difference check of the analytic Hessian.
        direction = np.arange(1, 3*n+1, dtype=float).reshape(n, 3)
        direction -= direction.mean(axis=0)
        direction /= np.linalg.norm(direction)
        step = 1e-6*max(1.0, np.linalg.norm(x))
        gp = response_at(x+step*direction)[2]
        gm = response_at(x-step*direction)[2]
        numerical = (gp-gm)/(2*step)
        analytic = H@direction.ravel()
        error = np.linalg.norm(numerical-analytic)/max(1.0, np.linalg.norm(analytic))
        assert error < 1e-7
        measured_eta2 = residual@residual/norm2
        assert abs(measured_eta2-float(geom['eta2_upper'])) < 1e-11
        spectral[name] = dict(rigidity_rank=rrank, stress_dimension=len(edges)-rrank,
                              relative_residual_outside_stress_space=relative_unbalanced,
                              gradient_norm=float(np.linalg.norm(grad)),
                              internal_curvature_eigenvalues=eigenvalues.tolist(),
                              prestress_relative_norm=float(np.linalg.norm(Hp, 2)/np.linalg.norm(H0, 2)),
                              hessian_difference_error=float(error), target_norm2=float(norm2))

    print('CELL 12 — RESPONSE-LAW SIZING; NO COEFFICIENT FIT OR NEW MASS PREDICTION')
    print('Trial form: q_eff = Q/Q_e + g*eta^2. Electron anchor has eta=0: g is unidentified.')
    print('Required g is target inversion ONLY; intervals use certified eta^2 bounds.')
    print(' case   probe       baseline q            required g interval')
    for row in sizing:
        a, b = row['required_g_bounds']
        print(f" {row['name']:6s} {row['probe']:>3s} {float(row['q']):16.9f} "
              f'[{float(a):.6e}, {float(b):.6e}]')
    for probe in ('A', 'B'):
        group = [r for r in sizing if r['probe'] == probe]
        low = max(r['required_g_bounds'][0] for r in group)
        high = min(r['required_g_bounds'][1] for r in group)
        print(f'  {probe}: common g interval for all K5s at ONE muon target is empty: {low > high}.')
    print('\nLocal curvature of Phi=(1/2)*normalized squared-distance loss:')
    print(' case    rank R   self-stress dim   min internal curvature   residual outside stress space')
    for name, row in spectral.items():
        imbalance = row['relative_residual_outside_stress_space']
        label = 'zero mismatch' if imbalance is None else f'{imbalance:.3e}'
        print(f" {name:6s} {row['rigidity_rank']:7d} {row['stress_dimension']:17d} "
              f"{row['internal_curvature_eigenvalues'][0]:24.9f} {label:>29s}")
    assert baseline['summaries'] == saved[0] and geometry['candidates'] == saved[1]
    print('\nPASS: analytic Hessians agree with directional finite differences; prior results preserved.')
    print('These are FLOATING-POINT local curvature diagnostics at saved realizations, not exact equilibria.')
    print('Phi uses the existing per-candidate normalization; it is not a shared absolute spring stiffness.')
    print('Curvature is NOT frequency without a kinetic/clock law, and NOT an inertial mass.')
    print('144 is not used; the RMR topology-only transfer matrix and its f_n are unchanged.')
    return dict(protocol='additive normalized surcharge feasibility + geometric objective Hessian',
                required_coefficient_diagnostics=sizing, local_spectral=spectral,
                selected_coefficient=None, fitted_to_mass=False,
                physical_mass_calibration_changed=False, native_dynamics_changed=False)


dcu_mass_12 = _run_dcu_mass_12()
