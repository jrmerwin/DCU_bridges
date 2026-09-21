# CELL 18 — three fixed transfers: pion branching, recoil, and spectroscopy.
# Run after Cell 17. Standard library only; no network calls or extra files.
# Imported laws: LO pseudoscalar weak decay, SR two-body kinematics, and
# nonrelativistic Coulomb reduced-mass scaling. These are NOT native DCU laws.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from math import isclose, sqrt


def _run_dcu_mass_18():
    needed = ('dcu_mass_17', '_reference')
    absent = [k for k in needed if k not in globals()]
    if absent:
        raise RuntimeError('Run Cell 17 first. Missing: ' + ', '.join(absent))
    old = dcu_mass_17
    before = deepcopy((old, _reference))
    counts = Counter(r['sector'] for r in _reference['rows'])
    N = len(_reference['rows'])
    if (N, counts['S'], counts['I'], counts['G']) != (137, 81, 40, 16):
        raise ValueError('The frozen registry changed; no replacement selected.')
    expected = {
        'structural_control': {'mu': Fraction(411, 2), 'tau': Fraction(6987, 2)},
        'screened_primary': {'mu': Fraction(827, 4), 'tau': Fraction(13909, 4)}
    }
    if old['recipes'] != expected or old['primary_recipe'] != 'screened_primary':
        raise ValueError('Cell 17 mass prescriptions changed; no silent refit.')

    # Extend the species dictionary with the EXPLICIT source formulas:
    # screen(3).pdf Eqs. (14),(15),(30),(31): p0=1836, pi0=273,
    # and the SAME +20/137 correction. No kaon or Xi fitted base is imported.
    delta0 = Fraction(counts['I'], 8*N)
    proton0 = Fraction((N-1)*3**3, 2)
    pion0 = Fraction(2*N-1)
    recipes = deepcopy(old['recipes'])
    for name, r in recipes.items():
        correction = 4*delta0 if name == 'screened_primary' else Fraction(0)
        r.update(proton=proton0+correction, pion=pion0+correction)
    assert proton0 == 1836 and pion0 == 273 and 4*delta0 == Fraction(20, 137)

    # TWO separate borrowed dimensional inputs, for DIFFERENT bridges.
    # NIST CODATA 2022: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
    electron_MeV = Fraction('0.51099895069')
    electron_sigma_MeV = Fraction('0.00000000016')
    # Parthey et al. 2011, hydrogen hyperfine-centroid 1S-2S frequency:
    # https://arxiv.org/abs/1107.3101
    hydrogen_Hz, hydrogen_sigma_Hz = Fraction(2466061413187035), Fraction(10)
    # These printed inputs are empirical measurements, NOT exact physical constants.

    def predict(r, energy_unit, hydrogen_frequency):
        """Uses source mass ratios plus the TWO declared calibrations, no targets."""
        mu, pion, proton = (Fraction(r[k]) for k in ('mu', 'pion', 'proton'))
        if not (1 < mu < pion and proton > 1 and energy_unit > 0 and hydrogen_frequency > 0):
            raise ValueError('Need ordered positive masses and positive calibrations.')
        # pi -> mu nu, massless neutrino, nonradiative two-body rest-frame kinematics.
        momentum = (pion*pion-mu*mu)/(2*pion)  # p/(m_e c)
        energy = (pion*pion+mu*mu)/(2*pion)    # E/(m_e c^2)
        kinetic = energy-mu
        beta = momentum/energy
        assert energy+momentum == pion and energy**2-momentum**2 == mu**2
        assert kinetic == (pion-mu)**2/(2*pion) and 0 < beta < 1
        # LO helicity-suppressed pseudoscalar decay: f_pi, V_ud and G_F cancel.
        pion_ratio = ((pion*pion-1)/(pion*pion-mu*mu))**2/(mu*mu)
        electron_momentum = (pion*pion-1)/(2*pion)
        assert pion_ratio == (electron_momentum/momentum)**2/(mu*mu)
        # Equal unit electric charges; gross-structure Coulomb law only.
        # nu_1S2S = A * reduced_mass/m_e. Calibrate A using H, not muonium.
        reduced_H, reduced_Mu = proton/(1+proton), mu/(1+mu)
        optical_coefficient = hydrogen_frequency/reduced_H
        muonium = optical_coefficient*reduced_Mu
        shift = hydrogen_frequency-muonium
        assert shift/hydrogen_frequency == (proton-mu)/(proton*(mu+1))
        return dict(
            pion_e_mu_ratio_LO=pion_ratio,
            recoil=dict(momentum_MeV_c=momentum*energy_unit,
                        muon_total_energy_MeV=energy*energy_unit,
                        muon_kinetic_energy_MeV=kinetic*energy_unit,
                        neutrino_energy_MeV=momentum*energy_unit,
                        beta=beta),
            spectroscopy=dict(muonium_1S2S_Hz=muonium,
                              hydrogen_minus_muonium_Hz=shift,
                              effective_Coulomb_coefficient_Hz=optical_coefficient),
            sensitivities=dict(
                dlog_p_dlog_mu=-2*mu*mu/(pion*pion-mu*mu),
                dlog_p_dlog_pion=(pion*pion+mu*mu)/(pion*pion-mu*mu),
                dlog_nu_Mu_dlog_mu=1/(mu+1),
                dlog_shift_dlog_mu=-mu/(proton-mu)-mu/(mu+1)))

    predictions = {name: predict(r, electron_MeV, hydrogen_Hz) for name, r in recipes.items()}
    frozen = deepcopy(predictions)
    for name, r in recipes.items():
        doubled_E = predict(r, 2*electron_MeV, hydrogen_Hz)
        doubled_H = predict(r, electron_MeV, 2*hydrogen_Hz)
        p = predictions[name]
        assert doubled_E['recoil']['momentum_MeV_c'] == 2*p['recoil']['momentum_MeV_c']
        assert doubled_E['recoil']['beta'] == p['recoil']['beta']
        assert doubled_E['pion_e_mu_ratio_LO'] == p['pion_e_mu_ratio_LO']
        assert doubled_H['spectroscopy']['muonium_1S2S_Hz'] == 2*p['spectroscopy']['muonium_1S2S_Hz']
        assert doubled_H['spectroscopy']['hydrogen_minus_muonium_Hz'] == 2*p['spectroscopy']['hydrogen_minus_muonium_Hz']

    # COMPARISONS ONLY — values and historical versions fixed, checked 2026-09-19.
    # PDG 2025 pion listing, pp. 1-2 and 5 (rendered pages inspected):
    # https://pdg.lbl.gov/2025/listings/rpp2025-list-pi-plus-minus.pdf
    # R is inclusive of radiative decays, whereas this pilot calculation is LO.
    observed_ratio, sigma_ratio = Fraction('0.00012327'), Fraction('0.00000023')
    observed_p, sigma_p = Fraction('29.79200'), Fraction('0.00011')
    # Meyer et al. 1999, a FROZEN historical measurement, not a claim of latest data:
    # https://arxiv.org/abs/hep-ex/9907013 ; 2 455 528 941.0(9.8) MHz.
    observed_Mu_Hz, sigma_Mu_Hz = Fraction('2455528941.0')*10**6, Fraction('9.8')*10**6
    observed_shift = hydrogen_Hz-observed_Mu_Hz
    sigma_shift = sqrt(float(hydrogen_sigma_Hz**2+sigma_Mu_Hz**2))

    # Matched measured-mass control diagnoses the limitations of the IMPORTED laws.
    # These masses enter ONLY this labeled control, never the primary predictions.
    reference_masses = dict(mu=Fraction('206.7682827'),
                            pion=Fraction('139.57039')/electron_MeV,
                            proton=Fraction('1836.152673426'))
    reference = predict(reference_masses, electron_MeV, hydrogen_Hz)
    targets = {
        'pion_e_mu_ratio_LO': (observed_ratio, float(sigma_ratio)),
        'muon_recoil_MeV_c': (observed_p, float(sigma_p)),
        'H_minus_Mu_1S2S_Hz': (observed_shift, sigma_shift),
        'muonium_1S2S_Hz': (observed_Mu_Hz, float(sigma_Mu_Hz))}
    def quantities(p):
        return dict(pion_e_mu_ratio_LO=p['pion_e_mu_ratio_LO'],
                    muon_recoil_MeV_c=p['recoil']['momentum_MeV_c'],
                    H_minus_Mu_1S2S_Hz=p['spectroscopy']['hydrogen_minus_muonium_Hz'],
                    muonium_1S2S_Hz=p['spectroscopy']['muonium_1S2S_Hz'])
    comparisons = []
    all_results = dict(predictions, measured_mass_bridge_control=reference)
    for name, p in all_results.items():
        for quantity, value in quantities(p).items():
            observed, sigma = targets[quantity]
            comparisons.append(dict(recipe=name, quantity=quantity, prediction=value,
                                    observation=observed, observation_sigma=sigma,
                                    fractional_error=value/observed-1))
    assert predictions == frozen and (old, _reference) == before

    print('CELL 18 — FIXED PRESCRIPTIONS INTO RECOIL, SPECTROSCOPY, AND PION BRANCHING')
    for name, r in recipes.items():
        print(f"  {name}: mu={r['mu']}, tau={r['tau']}, pion={r['pion']}, proton={r['proton']}")
    print('Calibrations: m_e c^2=0.51099895069 MeV (kinematics); H 1S-2S (spectroscopy).')
    print('Pion branching ratio needs NO dimensional calibration. Cell 17 clock is unchanged, unused here.')
    print('\n recipe                         p_mu (MeV/c)    H-Mu shift (THz)      R_pi e/mu')
    for name, p in all_results.items():
        print(f" {name:30s} {float(p['recoil']['momentum_MeV_c']):13.9f}"
              f" {float(p['spectroscopy']['hydrogen_minus_muonium_Hz'])/1e12:19.9f}"
              f" {float(p['pion_e_mu_ratio_LO']):17.11e}")
    print(f' Observation                    {float(observed_p):13.9f}'
          f' {float(observed_shift)/1e12:19.9f} {float(observed_ratio):17.11e}')
    print('\n recipe                         recoil error %     shift error %     R_pi error %')
    for name in all_results:
        errors = {c['quantity']: c['fractional_error'] for c in comparisons if c['recipe'] == name}
        print(f" {name:30s} {100*float(errors['muon_recoil_MeV_c']):+14.7f}"
              f" {100*float(errors['H_minus_Mu_1S2S_Hz']):+17.7f}"
              f" {100*float(errors['pion_e_mu_ratio_LO']):+17.7f}")
    primary = predictions['screened_primary']
    recoil, spec = primary['recoil'], primary['spectroscopy']
    print('\nPrimary dependent kinematic outputs (NOT separate experimental tests):')
    print(f"  muon kinetic energy={float(recoil['muon_kinetic_energy_MeV']):.9f} MeV;")
    print(f"  neutrino energy={float(recoil['neutrino_energy_MeV']):.9f} MeV; v_mu/c={float(recoil['beta']):.9f}.")
    print(f"Primary muonium line={float(spec['muonium_1S2S_Hz'])/1e12:.9f} THz;")
    print(f"  observed={float(observed_Mu_Hz)/1e12:.9f} THz; error="
          f"{1e6*float(spec['muonium_1S2S_Hz']/observed_Mu_Hz-1):+.6f} ppm.")
    print('The line and H-minus-Mu shift are ONE comparison, not independent successes.')
    print('Pion recoil contributes to standard pion-mass determinations; not independent of mass-fit data.')
    print('Muonium spectroscopy can determine m_mu/m_e; this is not a historically blind target.')
    print('LO pion ratio omits EM radiation; measured-mass LO control exposes a similar discrepancy.')
    print('Spectroscopy omits differential relativistic, recoil, QED and nuclear-size corrections.')
    print('No 1% pass count, joint significance, fitted correction, or theoretical noise allowance.')
    print('PASS: exact rational predictions, on-shell kinematics, equivalent ratios, calibration scaling, unchanged inputs.')
    return dict(protocol='source-fixed masses + borrowed kinematics/Coulomb/pseudoscalar-weak laws',
                primary_recipe='screened_primary', recipes=recipes, predictions=predictions,
                calibrations=dict(electron_rest_energy_MeV=electron_MeV,
                                  electron_rest_energy_sigma_MeV=electron_sigma_MeV,
                                  hydrogen_1S2S_Hz=hydrogen_Hz, hydrogen_sigma_Hz=hydrogen_sigma_Hz),
                comparisons=comparisons, comparison_targets=targets,
                measured_mass_control=dict(inputs=reference_masses, outputs=reference),
                prediction_function=predict, truly_blinded=False, native_dynamics_changed=False,
                radiative_corrections_applied=False, joint_significance_computed=False,
                original_Cell17_changed=False, source_checked_date='2026-09-19')


dcu_mass_18 = _run_dcu_mass_18()
