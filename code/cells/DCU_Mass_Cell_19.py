# CELL 19 — frozen angle ladder; pion-specific pointlike QED; baryon readiness.
# Run after Cell 18. Standard library only, with no network calls or new data files.
# The screening formulas are SOURCE HYPOTHESES. The QED kernel is IMPORTED physics.
# No coefficient, octant, or source formula is changed to improve a comparison.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from math import fsum, isclose, isfinite, log, log1p, pi, sqrt


def _run_dcu_mass_19():
    needed = ('dcu_mass_18', '_reference')
    absent = [key for key in needed if key not in globals()]
    if absent:
        raise RuntimeError('Run Cell 18 first. Missing: ' + ', '.join(absent))
    old = dcu_mass_18
    snapshot = deepcopy((old, _reference))
    N = len(_reference['rows'])
    sectors = Counter(r['sector'] for r in _reference['rows'])
    if (N, sectors['S'], sectors['I'], sectors['G']) != (137, 81, 40, 16):
        raise ValueError('The frozen registry changed; no replacement selected.')
    expected = {
        'structural_control': dict(mu=Fraction(411, 2), tau=Fraction(6987, 2),
                                   pion=Fraction(273), proton=Fraction(1836)),
        'screened_primary': dict(mu=Fraction(827, 4), tau=Fraction(13909, 4),
                                 pion=Fraction(37421, 137), proton=Fraction(251552, 137))}
    if old['recipes'] != expected or old['primary_recipe'] != 'screened_primary':
        raise ValueError('Cell 18 mass prescription changed; no silent refit.')

    # screen(3).pdf Eqs. (36)-(39); use equations, not rounded table entries.
    delta = Fraction(sectors['I'], 8*N)
    reactor_sin = 4*delta
    angles = dict(weak=Fraction(1, 4)-delta/2,
                  solar=Fraction(1, 4)+3*delta/2,
                  atmospheric=Fraction(1, 2)+3*delta/2,
                  reactor=reactor_sin**2)
    assert delta == Fraction(5, 137) and angles['reactor'] == Fraction(400, 18769)
    assert all(0 < value < 1 for value in angles.values())
    assert angles['atmospheric']-angles['solar'] == Fraction(1, 4)
    assert angles['solar']+3*angles['weak'] == 1
    assert reactor_sin+8*angles['weak'] == 2

    # NEW external input, not the RMR screening unit: low-energy QED alpha.
    # CODATA 2022, https://physics.nist.gov/cuu/Constants/Table/allascii.txt
    alpha_inverse, alpha_inverse_sigma = Fraction('137.035999177'), Fraction('0.000000021')
    alpha = float(1/alpha_inverse)

    def li2(x):
        """Real dilogarithm on [0,1], reflection + bounded positive series."""
        if not isfinite(x) or not 0 <= x <= 1:
            raise ValueError('Dilogarithm argument must be in [0,1].')
        if x == 0:
            return 0.0
        if x == 1:
            return pi*pi/6
        if x > 0.5:
            return pi*pi/6-log(x)*log1p(-x)-li2(1-x)
        terms, power = [], x
        for k in range(1, 10000):
            terms.append(power/(k*k))
            power *= x
            if power/((k+1)**2*(1-x)) < 1e-18:
                return fsum(terms)
        raise ArithmeticError('Dilogarithm series did not converge within its budget.')

    def qed_F(z):
        # Cirigliano-Rosell arXiv:0707.4464 Eqs. (96)-(97): pointlike inclusive kernel.
        if not isfinite(z) or not 0 < z < 1:
            raise ValueError('Radiative kernel needs 0 < squared mass ratio < 1.')
        L, M = log(z), log1p(-z)
        return (1.5*L+(13-19*z)/(8*(1-z))
                -(8-5*z)*z*L/(4*(1-z)**2)
                -(2+(1+z)*L/(1-z))*M
                -2*(1+z)*li2(1-z)/(1-z))

    def pion_prediction(masses, coupling):
        mu, pion_mass = Fraction(masses['mu']), Fraction(masses['pion'])
        if not 1 < mu < pion_mass or not isfinite(coupling) or coupling < 0:
            raise ValueError('Need ordered masses and a nonnegative finite coupling.')
        tree = ((pion_mass**2-1)/(pion_mass**2-mu**2))**2/mu**2
        correction = coupling/pi*(qed_F(float(1/pion_mass**2))
                                  -qed_F(float((mu/pion_mass)**2)))
        return dict(tree=tree, relative_pointlike_QED=correction,
                    corrected=float(tree)*(1+correction))

    mass_sets = dict(old['recipes'], measured_mass_bridge_control=old['measured_mass_control']['inputs'])
    pion = {name: pion_prediction(masses, alpha) for name, masses in mass_sets.items()}
    for name, masses in mass_sets.items():
        prior = (old['measured_mass_control']['outputs'] if name == 'measured_mass_bridge_control'
                 else old['predictions'][name])
        assert pion[name]['tree'] == prior['pion_e_mu_ratio_LO']
        assert pion_prediction(masses, 0.0)['corrected'] == float(pion[name]['tree'])
        assert isclose(pion_prediction(masses, 2*alpha)['relative_pointlike_QED'],
                       2*pion[name]['relative_pointlike_QED'], rel_tol=2e-15)
    assert isclose(li2(0.5), pi*pi/12-log(2)**2/2, rel_tol=2e-15)
    for x in (0.001, 0.1, 0.3, 0.8):
        assert isclose(li2(x)+li2(1-x), pi*pi/6-log(x)*log1p(-x), rel_tol=2e-15)

    # Lambda preflight: Eq. (20) fixes a mass; it supplies NO hadronic form factors.
    electron_MeV = old['calibrations']['electron_rest_energy_MeV']
    r_lambda, r_proton = N*sectors['G']-3**2+3**2*delta, old['recipes']['screened_primary']['proton']
    lambda_mass = r_lambda*electron_MeV
    endpoint = (r_lambda**2+1-r_proton**2)/(2*r_lambda)*electron_MeV
    preflight = dict(r_lambda=r_lambda, predicted_mass_MeV=lambda_mass,
                     electron_total_endpoint_MeV=endpoint,
                     electron_kinetic_endpoint_MeV=endpoint-electron_MeV,
                     rate=None, reason='Hadronic f_i(q^2), g_i(q^2), and V_us not supplied',
                     endpoint_compared_to_data=False,
                     Xi_bases_status='Source Eq. (24)-(25) bases were obtained from experimental masses')
    frozen = deepcopy((angles, pion, preflight))

    # COMPARISON DATA: versioned benchmarks, NOT claims of latest fit or blind targets.
    # NuFIT 6.0, Sept 2024 data, arXiv:2410.05380v2 Table 1.
    # Tuples = best fit, lower/upper 1-sigma errors, lower/upper marginal 3-sigma limits.
    fits = {
        'IC24_with_SK_NO': dict(solar=(.308,.011,.012,.275,.345),
            atmospheric=(.470,.013,.017,.435,.585), reactor=(.02215,.00058,.00056,.02030,.02388)),
        'IC19_without_SK_NO': dict(solar=(.307,.011,.012,.275,.345),
            atmospheric=(.561,.015,.012,.430,.596), reactor=(.02195,.00058,.00054,.02023,.02376)),
        'IC24_with_SK_IO': dict(solar=(.308,.011,.012,.275,.345),
            atmospheric=(.550,.015,.012,.440,.584), reactor=(.02231,.00056,.00056,.02060,.02409)),
        'IC19_without_SK_IO': dict(solar=(.308,.011,.012,.275,.345),
            atmospheric=(.562,.015,.012,.437,.597), reactor=(.02224,.00057,.00056,.02053,.02397))}
    comparisons = []
    for variant, entries in fits.items():
        for name, (center, lower, upper, low3, high3) in entries.items():
            p = float(angles[name])
            comparisons.append(dict(variant=variant, quantity=name, prediction=p,
                central=center, error_minus=lower, error_plus=upper,
                marginal_3sigma=(low3,high3), in_marginal_3sigma=low3 <= p <= high3,
                fractional_error=p/center-1, gaussian_pull=None))
    # PDG2025 Electroweak review Eq. (10.65): collider effective leptonic angle.
    # JUNO first 59.1-day result, arXiv:2511.14593v1: normal ordering, solar angle.
    simple_data = {'weak_PDG2025': ('weak', .23148, .00013),
                   'solar_JUNO_59days': ('solar', .3092, .0087)}
    simple = []
    for label, (name, center, sigma) in simple_data.items():
        p = float(angles[name])
        simple.append(dict(dataset=label, quantity=name, prediction=p, central=center,
                           sigma=sigma, fractional_error=p/center-1,
                           quoted_error_units=(p-center)/sigma))
    obs, obs_sigma = old['comparison_targets']['pion_e_mu_ratio_LO']
    pion_comparisons = {name: dict(tree_error=float(r['tree']/obs-1),
        corrected_error=r['corrected']/float(obs)-1) for name, r in pion.items()}
    lambda_observed, lambda_sigma = Fraction('1115.683'), Fraction('0.006')  # PDG2025 Lambda listing.
    lambda_check = dict(observed_MeV=lambda_observed, sigma_MeV=lambda_sigma,
                        prediction_minus_observed_keV=(lambda_mass-lambda_observed)*1000)
    assert (angles, pion, preflight) == frozen and (old, _reference) == snapshot

    print('CELL 19 — FROZEN ANGLES; PROCESS-SPECIFIC PION QED; HYPERON PREFLIGHT')
    print(f'delta0={delta}. Source Eqs. (36)-(39), no alternate angle formulas or octant flip.')
    for name, value in angles.items():
        print(f'  sin^2 {name:11s} = {str(value):>12s} = {float(value):.12f}')
    print('\nNuFIT 6.0 (Sept 2024 inputs): NO primary; IO is an unchanged-prediction sensitivity check.')
    print(' variant                 quantity      best fit     prediction   error %   in marginal 3sigma')
    for row in comparisons:
        print(f" {row['variant']:23s} {row['quantity']:11s} {row['central']:.6f} "
              f" {row['prediction']:.9f} {100*row['fractional_error']:+9.4f} {str(row['in_marginal_3sigma']):>17s}")
    print('No Gaussian pulls from the multimodal atmospheric fit; no joint-confidence claim.')
    for row in simple:
        print(f"  {row['dataset']}: {row['central']} +/- {row['sigma']}; "
              f"error={100*row['fractional_error']:+.5f}%; quoted-error units={row['quoted_error_units']:+.4f}")
    print('Exact identities: s23-s12=1/4; s12+3*s_eff=1; sin(theta13)+8*s_eff=2.')
    print('These identities are consequences, not additional independent successes. CP phase is unspecified.')
    print('\nPION ONLY: imported inclusive pointlike O(alpha) correction, recomputed with EACH mass recipe.')
    print(f'New measured input alpha^-1={alpha_inverse}; not delta0 or an RMR screening multiplier.')
    print(' recipe                         correction %        corrected R      old error %    new error %')
    for name, r in pion.items():
        c = pion_comparisons[name]
        print(f" {name:30s} {100*r['relative_pointlike_QED']:+11.7f} {r['corrected']:18.11e} "
              f"{100*c['tree_error']:+12.6f} {100*c['corrected_error']:+12.6f}")
    print(f'  Unchanged inclusive comparison: {float(obs):.8e} +/- {obs_sigma:.2e}.')
    print('Higher chiral/structure and resummed-log terms omitted, NOT covered by an assigned theory error.')
    print('This pion correction is NOT applied to the effective weak angle.')
    print('\nLAMBDA: source formula and borrowed electron unit, not rounded source prediction.')
    print(f'  Predicted mass={float(lambda_mass):.9f} MeV; residual={float(lambda_check["prediction_minus_observed_keV"]):+.6f} keV.')
    print(f'  Pure kinematic electron kinetic endpoint={float(endpoint-electron_MeV):.9f} MeV (no endpoint data test).')
    print('  Absolute hyperon rates NOT calculated: form factors/CKM missing; Xi bases are experimentally inferred in source.')
    print('PASS: exact ladder identities, prior pion replay, dilogarithm identities, alpha scaling, and unchanged inputs.')
    print('No mass/correction fit, native evolution, new weak-angle radiative prescription, or historical blindness claimed.')
    return dict(protocol='source-fixed angle ladder + imported pion-specific pointlike QED',
        primary_NuFIT_variant='IC24_with_SK_NO', angle_predictions=angles, delta0=delta,
        angle_comparisons=comparisons, supplementary_angle_comparisons=simple,
        pion_predictions=pion, pion_comparisons=pion_comparisons,
        qed_input=dict(alpha_inverse=alpha_inverse, inverse_sigma=alpha_inverse_sigma),
        baryon_preflight=preflight, lambda_mass_check=lambda_check,
        helpers=dict(dilogarithm=li2, radiative_kernel=qed_F, pion_predictor=pion_prediction),
        full_PMNS_predicted=False, hyperon_rate_predicted=False,
        higher_pion_corrections_applied=False, weak_angle_QED_factor_applied=False,
        theoretical_error_assigned=False, joint_significance_computed=False,
        native_dynamics_changed=False, prior_prescriptions_changed=False, truly_blinded=False)


dcu_mass_19 = _run_dcu_mass_19()
