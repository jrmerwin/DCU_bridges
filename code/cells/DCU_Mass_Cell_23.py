# CELL 23 — production thresholds and a finite-quantum operator audit.
# Run after Cells 18 and 19 (Cell 22 remains an unchanged frequency checkpoint).
# Requires NumPy. No internet or extra data files. ONE copy-paste notebook cell.
# Native registry data / source mass hypotheses / borrowed physics stay separate.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations
from math import cos, isclose, pi, sqrt


def _run_dcu_mass_23():
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError('Cell 23 needs NumPy in this notebook kernel.') from exc
    needed = ('_reference', 'dcu_mass_18', 'dcu_mass_19')
    missing = [k for k in needed if k not in globals()]
    if missing:
        raise RuntimeError('Run Cells 18 and 19 first. Missing: ' + ', '.join(missing))
    prior = (_reference, dcu_mass_18, dcu_mass_19)
    snapshot = deepcopy(prior)
    Q = Fraction
    N = len(_reference['rows'])
    counts = Counter(r['sector'] for r in _reference['rows'])
    if (N, counts['S'], counts['I'], counts['G']) != (137, 81, 40, 16):
        raise ValueError('The original registry changed.')
    E0 = dcu_mass_18['calibrations']['electron_rest_energy_MeV']
    if E0 != Q('0.51099895069'):
        raise ValueError('Earlier electron-energy calibration changed.')
    expected = dict(mu=Q(827, 4), tau=Q(13909, 4), proton=Q(251552, 137),
                    pion=Q(37421, 137))
    if dcu_mass_18['recipes']['screened_primary'] != expected:
        raise ValueError('Frozen screened masses changed; no substitute selected.')
    delta = dcu_mass_19['delta0']
    rL = dcu_mass_19['baryon_preflight']['r_lambda']
    if delta != Q(5, 137) or rL != 2183 + 9*delta:
        raise ValueError('Frozen screening/Lambda prescription changed.')

    # ONE NEW borrowed species mass for associated-strangeness channels ONLY.
    # PDG 2025 K+- listing p.6: OUR FIT 493.677 +/- .015 MeV.
    # https://pdg.lbl.gov/2025/listings/rpp2025-list-K-plus-minus.pdf
    # The source's unexplained kaon base 966 is NOT promoted to a prediction here.
    mK, mK_sigma = Q('493.677'), Q('0.015')
    mass_sets = {}
    for label, recipe in dcu_mass_18['recipes'].items():
        masses = {k: Q(recipe[k])*E0 for k in ('mu', 'tau', 'proton', 'pion')}
        masses.update(electron=E0, photon=Q(0), kaon=mK,
                      Lambda=(rL if label == 'screened_primary' else Q(2183))*E0)
        mass_sets[label] = masses

    # External quantum numbers are channel bookkeeping, NOT registry outputs.
    # Store electric charge, baryon number, strangeness; all selected processes
    # conserve these. These necessary checks do NOT assert a nonzero amplitude.
    species = {
        'e-': ('electron', -1, 0, 0), 'e+': ('electron', 1, 0, 0),
        'mu-': ('mu', -1, 0, 0), 'mu+': ('mu', 1, 0, 0),
        'tau-': ('tau', -1, 0, 0), 'tau+': ('tau', 1, 0, 0),
        'p': ('proton', 1, 1, 0), 'pbar': ('proton', -1, -1, 0),
        'Lambda': ('Lambda', 0, 1, -1), 'Lambdabar': ('Lambda', 0, -1, 1),
        'K+': ('kaon', 1, 0, 1), 'gamma': ('photon', 0, 0, 0),
        'pi+': ('pion', 1, 0, 0), 'pi-': ('pion', -1, 0, 0)}
    def charges(particles):
        return tuple(sum(species[p][j] for p in particles) for j in (1, 2, 3))
    def require_conservation(initial, final):
        if charges(initial) != charges(final):
            raise ValueError('Specified channel violates charge/baryon/strangeness bookkeeping.')

    # Endothermic free-particle kinematic boundaries, narrow-width approximation.
    # NOT cross sections, rates, bound-state thresholds, or beam-response fits.
    channels = (
        ('ee_to_mumu', ('e+', 'e-'), ('mu+', 'mu-')),
        ('ee_to_tautau', ('e+', 'e-'), ('tau+', 'tau-')),
        ('ee_to_ppbar', ('e+', 'e-'), ('p', 'pbar')),
        ('ee_to_Lambda_pair', ('e+', 'e-'), ('Lambda', 'Lambdabar')),
        ('gamma_p_to_p_mumu', ('gamma', 'p'), ('p', 'mu+', 'mu-')),
        ('gamma_p_to_K_Lambda', ('gamma', 'p'), ('K+', 'Lambda')),
        ('pp_to_p_K_Lambda', ('p', 'p'), ('p', 'K+', 'Lambda')))
    def threshold(masses, initial, final):
        require_conservation(initial, final)
        a, b = (masses[species[p][0]] for p in initial)
        outgoing = tuple(masses[species[p][0]] for p in final)
        M = sum(outgoing, Q(0))
        if b <= 0 or M <= a+b:
            raise ValueError('This helper is for endothermic production on a massive rest target.')
        beam_E = (M*M-a*a-b*b)/(2*b)
        beam_T = beam_E-a
        p2 = beam_E*beam_E-a*a
        assert beam_T > 0 and p2 > 0 and a*a+b*b+2*b*beam_E == M*M
        assert beam_T == (M*M-(a+b)**2)/(2*b)
        # At threshold, all final particles share the same four-velocity.
        total_E = beam_E+b
        final_E = tuple(m*total_E/M for m in outgoing)
        final_p_fractions = tuple(m/M for m in outgoing)
        assert sum(final_E) == total_E and sum(final_p_fractions) == 1
        assert all(E*E-f*f*p2 == m*m
                   for E, f, m in zip(final_E, final_p_fractions, outgoing))
        assert a*a+b*b+2*b*(beam_E-Q(1, 1000)) < M*M
        return dict(COM_energy_MeV=M, s_threshold_MeV2=M*M,
                    beam_total_energy_MeV=beam_E, beam_kinetic_energy_MeV=beam_T,
                    beam_momentum_squared_MeV2=p2,
                    kaon_mass_borrowed=any(species[p][0] == 'kaon' for p in final))

    # Predict first; benchmarks will be introduced separately below.
    thresholds = {name: {key: threshold(m, incoming, outgoing)
                        for key, incoming, outgoing in channels}
                  for name, m in mass_sets.items()}
    frozen_thresholds = deepcopy(thresholds)
    rejected = []
    for pion in ('pi+', 'pi-'):
        try:
            require_conservation((pion, 'p'), ('K+', 'Lambda'))
        except ValueError:
            rejected.append(pion)
        else:
            raise AssertionError('A charge-violating pion channel was accepted.')

    # Measured-mass BENCHMARKS, not independent measured production thresholds.
    # CODATA2022 mu/e and p/e; PDG2025 tau and Lambda listings, page1.
    # https://physics.nist.gov/cuu/Constants/Table/allascii.txt
    # https://pdg.lbl.gov/2025/listings/rpp2025-list-tau.pdf
    # https://pdg.lbl.gov/2025/listings/rpp2025-list-lambda.pdf
    measured = dict(mass_sets['screened_primary'], mu=Q('206.7682827')*E0,
                    proton=Q('1836.152673426')*E0,
                    tau=Q('1776.93'), Lambda=Q('1115.683'))
    benchmarks = {key: threshold(measured, incoming, outgoing)
                  for key, incoming, outgoing in channels}
    summaries = []
    for key, incoming, outgoing in channels:
        primary, control, benchmark = (thresholds['screened_primary'][key],
                                       thresholds['structural_control'][key], benchmarks[key])
        summaries.append(dict(channel=key, incoming=incoming, outgoing=outgoing,
            primary=primary, structural_control=control, measured_mass_benchmark=benchmark,
            COM_fractional_error=primary['COM_energy_MeV']/benchmark['COM_energy_MeV']-1,
            lab_fractional_error=primary['beam_kinetic_energy_MeV']/benchmark['beam_kinetic_energy_MeV']-1))

    def two_body(s, a, b):
        s, a, b = Q(s), Q(a), Q(b)
        if s <= 0 or min(a, b) < 0:
            raise ValueError('Need positive invariant s and nonnegative masses.')
        if s < (a+b)**2:
            return dict(allowed=False, momentum_MeV=None, beta_first=None)
        p2 = (s-(a+b)**2)*(s-(a-b)**2)/(4*s)
        EA2 = (s+a*a-b*b)**2/(4*s)
        EB2 = (s+b*b-a*a)**2/(4*s)
        assert EA2-p2 == a*a and EB2-p2 == b*b
        return dict(allowed=True, momentum_squared_MeV2=p2,
                    momentum_MeV=sqrt(float(p2)), beta_first=sqrt(float(p2/EA2)))
    # Fixed illustrative energies, not observed beam points or extra successes.
    scenarios = []
    for W in (3554, 3560, 3600):
        row = {label: two_body(Q(W)**2, m['tau'], m['tau'])
               for label, m in dict(mass_sets, measured_mass_benchmark=measured).items()}
        scenarios.append(dict(process='tau_pair', COM_energy_MeV=W, outputs=row))
    for Eg in (900, 912, 920, 1000):
        row = {label: two_body(m['proton']**2+2*m['proton']*Eg, mK, m['Lambda'])
               for label, m in dict(mass_sets, measured_mass_benchmark=measured).items()}
        scenarios.append(dict(process='gamma_p_to_K_Lambda', photon_energy_MeV=Eg, outputs=row))
    for m in mass_sets.values():
        assert two_body(4*m['tau']**2, m['tau'], m['tau'])['momentum_squared_MeV2'] == 0
        # Positive Kallen polynomial below the pseudothreshold is NOT physical.
        assert not two_body((m['Lambda']-mK)**2/2, mK, m['Lambda'])['allowed']
    for key, incoming, outgoing in channels:
        doubled = threshold({k: 2*v for k,v in mass_sets['screened_primary'].items()}, incoming, outgoing)
        assert doubled['beam_kinetic_energy_MeV'] == 2*thresholds['screened_primary'][key]['beam_kinetic_energy_MeV']

    # QUANTUM PREFLIGHT: actual overlap matrix, no native time-to-unitary assertion.
    A = np.zeros((N, N))
    for i, neighbors in enumerate(_reference['adjacency']):
        A[i, list(neighbors)] = 1.0
    degree = A.sum(axis=0)
    assert np.array_equal(A, A.T) and Counter(map(int, degree)) == Counter({73:126, 136:11})
    T = -A/degree[np.newaxis, :]
    B = -A/np.sqrt(degree[:, None]*degree[None, :])
    assert np.allclose(B, (T*np.sqrt(degree)[None, :])/np.sqrt(degree)[:, None], atol=2e-16)
    assert np.allclose(B, B.T, atol=2e-16)
    column_norm2 = (T*T).sum(axis=0)
    assert np.allclose(column_norm2, 1/degree, atol=2e-16)
    ev, vec = np.linalg.eigh(B)
    expected_ev = np.sort([-1, -62/73, 383/4964] + [1/73]*124 + [1/136]*10)
    assert np.max(np.abs(ev-expected_ev)) < 2e-14
    # Two EXPLICIT energy-map choices share eigenvectors but not cavity gaps.
    # Dimensionless tau here is NOT the native maintenance clock or seconds.
    unitary_checks = {}
    for label, energies in (('signed', ev), ('magnitude', np.abs(ev))):
        U = (vec*np.exp(-1j*energies)) @ vec.T
        half = (vec*np.exp(-0.5j*energies)) @ vec.T
        error = float(np.max(np.abs(U.conj().T@U-np.eye(N))))
        assert error < 3e-14 and np.max(np.abs(half@half-U)) < 3e-14
        unitary_checks[label] = error
    gaps = {n: dict(signed=1+Q(1,n-1), magnitude=1-Q(1,n-1)) for n in (3,4,5)}
    assert gaps[5]['signed']/gaps[3]['signed'] == Q(5,6)
    assert gaps[5]['magnitude']/gaps[3]['magnitude'] == Q(3,2)

    # Finite record-character algebra, using the DECLARED complex instrument.
    # This checks consequences of that representation, not emergence of Born's rule.
    finite = {}
    for q in (3, 9):
        j = np.arange(q)
        omega = np.exp(2j*pi/q)
        X = np.zeros((q,q), dtype=complex)
        X[(j+1)%q, j] = 1
        Z = np.diag(omega**j)
        F = omega**np.outer(j,j)/sqrt(q)
        assert np.max(np.abs(Z@X-omega*X@Z)) < 2e-14
        assert np.max(np.abs(F.conj().T@F-np.eye(q))) < 2e-14
        assert np.max(np.abs(np.abs(F)**2-1/q)) < 2e-14
        operators = np.array([np.linalg.matrix_power(X,a)@np.linalg.matrix_power(Z,b)
                              for a in range(q) for b in range(q)]).reshape(q*q,q*q)
        orth_error = float(np.max(np.abs(operators.conj()@operators.T/q-np.eye(q*q))))
        assert orth_error < 3e-14
        finite[q] = dict(dimension=q, operator_basis_size=q*q,
                         Weyl_relation='ZX=exp(2*pi*i/q)XZ', orthogonality_error=orth_error)
    q = 3
    j = np.arange(q)
    F = np.exp(2j*pi*np.outer(j,j)/q)/sqrt(q)
    e0 = np.array([1,0,0], dtype=complex)
    interference = {}
    for phase in (Q(0), Q(1,9)):
        phase_gate = np.diag(np.exp(2j*pi*j*float(phase)))
        amplitude = F.conj().T@phase_gate@F@e0
        probabilities = np.abs(amplitude)**2
        analytic = np.array([(1+2*cos(2*pi*(float(phase)-a/3)))**2/9 for a in range(3)])
        assert np.max(np.abs(probabilities-analytic)) < 1e-14
        assert abs(probabilities.sum()-1) < 1e-14
        dephased = F.conj().T@(np.eye(3)/3)@F
        assert np.max(np.abs(np.diag(dephased)-1/3)) < 1e-14
        interference[str(phase)] = tuple(map(float, probabilities))
    quantum = dict(T_is_Hermitian_in_Euclidean_basis=bool(np.allclose(T,T.T)),
        T_is_unitary=False, basis_column_norm2={str(d): Q(1,d) for d in (73,136)},
        symmetric_generator_spectrum=tuple(map(float,ev)), unitary_checks=unitary_checks,
        Kn_gap_choices=gaps, finite_character_checks=finite, interference=interference,
        phase_erased=(1/3,1/3,1/3), physical_energy_time_scale=None,
        note='Complex Hilbert representation, energy map, unitary evolution and Born readout are declared choices.')

    assert thresholds == frozen_thresholds and prior == snapshot
    print('CELL 23 — PRODUCTION THRESHOLDS + FINITE-QUANTUM PREFLIGHT')
    print('Neutrino cells untouched. Screened masses and electron unit unchanged.')
    print(f'Associated-strangeness channels borrow m_K+={mK} +/- {mK_sigma} MeV; no claim to derive base 966.')
    print('Incoming first particle is the beam; second is the rest target in lab outputs.')
    print('\n channel                         W_threshold(MeV)   beam kinetic threshold(MeV)  lab error vs mass benchmark (%)')
    for r in summaries:
        p = r['primary']
        print(f" {r['channel']:31s} {float(p['COM_energy_MeV']):16.6f}"
              f" {float(p['beam_kinetic_energy_MeV']):30.6f} {100*float(r['lab_fractional_error']):+16.8f}")
    print('Benchmarks are computed from measured masses, NOT independent measured thresholds.')
    print('Pi+ p and pi- p -> K+ Lambda rejected by charge conservation; pi0 would have the required charge.')
    print('\nIllustrative tau-pair speeds at fixed COM energies (not velocity measurements):')
    for r in scenarios[:3]:
        a,b = r['outputs']['screened_primary'],r['outputs']['measured_mass_benchmark']
        print(f"  W={r['COM_energy_MeV']} MeV: beta={a['beta_first']:.12f}, "
              f"mass-control beta={b['beta_first']:.12f}; relative difference={100*(a['beta_first']/b['beta_first']-1):+.6f}%")
    print('\nOPERATOR AUDIT: T=-A D^-1 is not a unitary step or an ordinary Hermitian matrix.')
    print('  ||T e_j||^2 = 1/d_j: 1/73 or 1/136; no Born normalization by squaring T entries.')
    print('  B=D^-1/2 T D^1/2 is symmetric; both exp(-i B tau) and exp(-i |B| tau) are unitary.')
    print('  Their K5/K3 cavity-gap ratios are 5/6 and 3/2 respectively; the energy map is extra information.')
    print('  Finite character spaces checked: q=3,9; q^2 orthogonal Weyl operators, mutually unbiased Fourier bases.')
    for phase,p in interference.items():
        print(f'  Qutrit phase {phase}: P={p}; complete phase erasure gives (1/3,1/3,1/3).')
    print('No physical phase-space area, new Bell claim, clock-to-seconds or Lorentz-violation coefficient inferred.')
    print('PASS: rational threshold conservation/scaling; legal charge bookkeeping; threshold/pseudothreshold distinction;')
    print('      actual registry similarity/spectrum; unitary alternatives; finite-character algebra; unchanged prior cells.')
    return dict(protocol='frozen masses + imported relativistic thresholds; separate finite-quantum algebra audit',
        mass_sets_MeV=mass_sets, thresholds=thresholds, measured_mass_benchmarks=benchmarks,
        measured_mass_inputs_MeV=measured, benchmark_uncertainties=dict(tau_MeV=Q('.09'),Lambda_MeV=Q('.006')),
        summaries=summaries, kinematic_scenarios=scenarios,
        kaon_input=dict(value_MeV=mK,sigma_MeV=mK_sigma,source='PDG 2025 K+- listing p.6 OUR FIT'),
        electron_energy_input_MeV=E0, rejected_charge_channels=rejected,
        quantum=quantum, dispersion_prediction=None,
        helpers=dict(threshold=threshold,two_body=two_body), independent_threshold_data_tested=False,
        native_dynamics_changed=False, quantum_postulates_derived=False, truly_blinded=False)


dcu_mass_23 = _run_dcu_mass_23()
