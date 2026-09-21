# CELL 13 — decompose the saved curvature and locate conditional service coupling.
# Run after Cells 11 and 12. NumPy only. No optimization or native history is run.
# A prospective per-witness correction is ANALYZED, not added to the DCU dynamics.
from copy import deepcopy
from itertools import combinations


def _run_dcu_mass_13():
    import numpy as np
    missing = [k for k in ('dcu_mass_11', 'dcu_mass_12') if k not in globals()]
    if missing:
        raise RuntimeError('Run Cells 11 and 12 first. Missing: ' + ', '.join(missing))
    geometry, saved_spectra = dcu_mass_11['candidates'], dcu_mass_12['local_spectral']
    before = (deepcopy(geometry), deepcopy(saved_spectra))
    results = {}

    def rank_of(singular, shape):
        if not len(singular) or singular[0] == 0:
            return 0
        tol = 100 * max(shape) * np.finfo(float).eps * singular[0]
        return int(np.sum(singular > tol))

    for name, case in geometry.items():
        target = np.asarray(case['source_D2'], dtype=float)
        n = len(target)
        edges = tuple(combinations(range(n), 2))
        C = sum(target[i, j]**2 for i, j in edges)
        if C <= 0:
            raise ValueError('The declared target normalization must be positive.')
        if 'best' in case:
            x = np.asarray(case['best']['coordinates'], dtype=float)
        else:
            B = np.asarray(case['source_certificate']['gram'], dtype=float)
            values, vectors = np.linalg.eigh(B)
            x = vectors[:, -3:] * np.sqrt(np.maximum(values[-3:], 0))
        x = x - x.mean(axis=0)
        R = np.zeros((len(edges), 3*n))
        residual = np.empty(len(edges))
        H_stress = np.zeros((3*n, 3*n))
        for k, (i, j) in enumerate(edges):
            delta = x[i]-x[j]
            residual[k] = delta @ delta-target[i, j]
            R[k, 3*i:3*i+3], R[k, 3*j:3*j+3] = 2*delta, -2*delta
            # This is r_e * Hessian(r_e) / C; it is not another fitted interaction.
            block = 2 * residual[k] * np.eye(3) / C
            a, b = slice(3*i, 3*i+3), slice(3*j, 3*j+3)
            H_stress[a, a] += block
            H_stress[b, b] += block
            H_stress[a, b] -= block
            H_stress[b, a] -= block
        H_material = R.T @ R / C
        H = H_material + H_stress
        rigid = np.column_stack(
            [np.tile(a, (n, 1)).ravel() for a in np.eye(3)] +
            [np.cross(np.tile(a, (n, 1)), x).ravel() for a in np.eye(3)])
        U, singular, _ = np.linalg.svd(rigid, full_matrices=True)
        rigid_rank = rank_of(singular, rigid.shape)
        if rigid_rank != 6:
            raise ValueError(f'{name}: expected exactly six rigid motions.')
        internal = U[:, 6:]
        H_internal = internal.T @ H @ internal
        H0_internal = internal.T @ H_material @ internal
        values, vectors = np.linalg.eigh(H_internal)
        material_values = np.linalg.eigvalsh(H0_internal)
        assert np.allclose(values, saved_spectra[name]['internal_curvature_eigenvalues'],
                           atol=2e-12, rtol=2e-10)
        assert np.isclose(C, saved_spectra[name]['target_norm2'], rtol=2e-14)
        assert values[0] > 0 and material_values[0] > 0
        rigidity_rank = rank_of(np.linalg.svd(R, compute_uv=False), R.shape)
        assert rigidity_rank == saved_spectra[name]['rigidity_rank']

        # g_e = gradient of the SINGLE edge contribution to Phi, at the saved fit.
        # The common factor/mobility of a prospective coordinate update is omitted.
        gradients = residual[:, None] * R / C
        gradient_sum = gradients.sum(axis=0)
        tensor = gradients.T @ gradients
        positive_mismatch = case['eta2_lower'] > 0
        raw_trace = float(np.trace(tensor))
        if positive_mismatch:
            trace = raw_trace
            ambient_modes = internal @ vectors
            modal_power = np.einsum('ij,ij->j', ambient_modes, tensor @ ambient_modes)
            shares = modal_power / trace
            assert np.all(shares >= -1e-14) and abs(shares.sum()-1) < 1e-10
            soft_share = float(shares[0])
            modal_shares = shares.tolist()
        else:
            # Source certificates establish exact compatibility. Residual gradients
            # of order 1e-15 in floating coordinates are roundoff, not a noise source.
            assert raw_trace < 1e-20
            trace, soft_share, modal_shares = 0.0, None, None

        # Finite enumeration checks ONLY a conditional sampling identity.
        # Test fixture: 12 requests, choose 2; unused rows have zero correction.
        # These are NOT a scientific choice of capacity or backlog for the DCU.
        population, draws = 12, 2
        marks = np.zeros((population, 3*n))
        marks[:len(edges)] = gradients
        samples = np.array([marks[list(ids)].sum(axis=0)
                            for ids in combinations(range(population), draws)])
        expected_mean = draws/population * gradient_sum
        centered = samples - samples.mean(axis=0)
        covariance = centered.T @ centered / len(samples)
        expected_covariance = (draws*(population-draws)/(population*(population-1))) * (
            tensor-np.outer(gradient_sum, gradient_sum)/population)
        assert np.allclose(samples.mean(axis=0), expected_mean, atol=1e-14, rtol=1e-10)
        assert np.allclose(covariance, expected_covariance, atol=1e-14, rtol=1e-10)

        results[name] = dict(
            n_vertices=n, edges=edges, coordinates=x.tolist(), target_norm2=float(C),
            residuals=residual.tolist(), internal_basis=internal.tolist(),
            rigidity_rank=rigidity_rank, self_stress_dimension=len(edges)-rigidity_rank,
            internal_eigenvalues=values.tolist(), material_eigenvalues=material_values.tolist(),
            unnormalized_eigenvalues=(C*values).tolist(),
            softest_curvature_change= float(values[0]/material_values[0]-1),
            H=H.tolist(), H_material=H_material.tolist(), H_prestress=H_stress.tolist(),
            edge_gradients=gradients.tolist(), gradient_sum_norm=float(np.linalg.norm(gradient_sum)),
            edge_gradient_tensor=tensor.tolist(), raw_tensor_trace=raw_trace,
            tensor_trace_using_exact_zero_certificate=trace,
            softest_mode_gradient_share=soft_share, all_mode_gradient_shares=modal_shares,
            source_positive_mismatch=bool(positive_mismatch))

    assert geometry == before[0] and saved_spectra == before[1]
    print('CELL 13 — SAVED SPECTRA, PRESTRESS, AND A CONDITIONAL SERVICE MARK')
    print('Phi, source distances, normalization, and saved coordinates are unchanged.')
    print('No new native history, coordinate trajectory, kinetic law, or mass fit.')
    print('\n case   modes    min H      max H    min H_material    softening (%)')
    for name, r in results.items():
        v, v0 = r['internal_eigenvalues'], r['material_eigenvalues']
        print(f" {name:6s} {len(v):5d} {v[0]:10.7f} {v[-1]:10.7f} {v0[0]:16.7f}"
              f" {100*r['softest_curvature_change']:16.7f}")
    print('\nCommon unnormalized loss check: min/max of C*H (NOT a new adopted energy).')
    for name, r in results.items():
        v = r['unnormalized_eigenvalues']
        print(f'  {name}: {v[0]:.9f}, {v[-1]:.9f}')
    print('\nPer-edge gradient tensor S=sum_e g_e g_e^T, g_e=r_e*grad(r_e)/C:')
    print(' case      trace S         share in softest internal mode (%)')
    for name, r in results.items():
        share = r['softest_mode_gradient_share']
        label = 'undefined (exact zero)' if share is None else f'{100*share:.9f}'
        print(f" {name:6s} {r['tensor_trace_using_exact_zero_certificate']:13.6e} {label:>42s}")
    print('\nProspective rule ONLY: a served witness applies -epsilon*g_e to attached coordinates.')
    print('The tensor is NOT a measured native noise or extra recording-work rate.')
    print('Native uniform service fixes selection probabilities; it does not supply epsilon or this rule.')
    print('PASS: all 7 saved spectra; six rigid modes; source-certified zero controls;')
    print('      conditional service mean/covariance checked by finite enumeration; old results preserved.')
    return dict(protocol='Cell 12 decomposition and per-witness gradient diagnostic', cases=results,
                new_native_history=False, new_coordinate_dynamics=False,
                mobility_selected=None, physical_mass_calibration_changed=False,
                sampling_identity='q(R-q)/(R(R-1)) * (sum g_e g_e^T - g_sum g_sum^T/R)')


dcu_mass_13 = _run_dcu_mass_13()
