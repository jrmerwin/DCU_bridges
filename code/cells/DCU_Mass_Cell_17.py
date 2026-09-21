# CELL 17 — locked mass recipes -> tau leptonic rates, with ONE clock input.
# Run after Cell 14. Standard library only. No network calls or data files needed.
# IMPORTED PHYSICS: universal leading-order V-A three-body decay law.
# This is NOT native DCU time evolution. Tau data are comparison-only.
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from math import isclose, isfinite, isqrt, log, sqrt


def _run_dcu_mass_17():
    needed = ('_reference', 'dcu_mass_14')
    missing = [k for k in needed if k not in globals()]
    if missing:
        raise RuntimeError('Run Cells 1 and 14 first. Missing: ' + ', '.join(missing))
    snapshot = deepcopy((_reference, dcu_mass_14))
    counts = Counter(r['sector'] for r in _reference['rows'])
    N = len(_reference['rows'])
    if (N, counts['S'], counts['I'], counts['G']) != (137, 81, 40, 16):
        raise ValueError('Registry changed: this frozen prescription is not applicable.')

    # SOURCE PRESCRIPTIONS, not a new derivation from registry counts.
    # leptons.pdf Eqs. 9/10: r_mu^0=3N/2; r_tau^0=17*r_mu^0.
    # screen(3).pdf Eqs. 7/9: +5/4 and -13*(5/4).
    # Do NOT substitute that paper's unexplained 205.518 / 3493.480 bases.
    f3, f5 = Fraction(1, 2), Fraction(3, 4)
    generation = Fraction(N-1, 8)
    rmu0 = N*f5/f3
    if dcu_mass_14['prescriptions']['muon_electron'] != rmu0:
        raise ValueError('The Cell 14 muon prescription changed; no silent refit.')
    rtau0 = generation*rmu0
    delta0 = Fraction(counts['I'], 8*N)
    rootG = isqrt(counts['G'])
    assert rootG**2 == counts['G']
    delta_mu = delta0*N/rootG
    B_tau = (5-1)**2 + 1 - rootG
    recipes = {
        'structural_control': {'mu': rmu0, 'tau': rtau0},
        'screened_primary': {'mu': rmu0+delta_mu, 'tau': rtau0-B_tau*delta_mu},
    }
    assert recipes['screened_primary'] == {'mu': Fraction(827, 4), 'tau': Fraction(13909, 4)}

    # SOLE measured input to the predictions: PDG 2025 muon lifetime, page 2.
    # https://pdg.lbl.gov/2025/listings/rpp2025-list-muon.pdf
    tau_mu, sigma_tau_mu = 2.1969811e-6, 2.2e-12  # seconds

    def phase(x):
        """Three-body V-A phase-space factor; x=(daughter/parent mass)^2."""
        if not isfinite(x) or not 0 <= x <= 1:
            raise ValueError('Phase argument must lie in [0,1].')
        if x == 0:
            return 1.0
        if x == 1:
            return 0.0
        return 1 - 8*x + 8*x**3 - x**4 - 12*x*x*log(x)

    def predict(recipe, clock):
        """No tau measurements or measured lepton masses enter this function."""
        mu, tau = float(recipe['mu']), float(recipe['tau'])
        if not (0 < clock and 1 < mu < tau):
            raise ValueError('Need positive clock and m_e < m_mu < m_tau.')
        normalizer = clock*mu**5*phase(1/mu**2)
        rates = {name: tau**5*phase((daughter/tau)**2)/normalizer
                 for name, daughter in (('electron', 1.0), ('muon', mu))}
        assert all(isfinite(v) and v > 0 for v in rates.values())
        return dict(mass_ratio_tau_mu=tau/mu, rate_s_inverse=rates,
                    partial_lifetime_ps={k: 1e12/v for k, v in rates.items()},
                    muon_over_electron_rate=rates['muon']/rates['electron'],
                    calibrated_coefficient_s_inverse=1/normalizer)

    # Freeze BOTH predictions before introducing comparison values.
    predictions = {name: predict(recipe, tau_mu) for name, recipe in recipes.items()}
    frozen_predictions = deepcopy(predictions)
    for name, recipe in recipes.items():
        result = predictions[name]
        mu, tau = float(recipe['mu']), float(recipe['tau'])
        for channel, daughter in (('electron', 1.0), ('muon', mu)):
            direct = (tau/mu)**5*phase((daughter/tau)**2)/phase(1/mu**2)/tau_mu
            assert isclose(direct, result['rate_s_inverse'][channel], rel_tol=2e-14)
        doubled = predict(recipe, 2*tau_mu)
        assert all(isclose(doubled['rate_s_inverse'][c], value/2, rel_tol=2e-14)
                   for c, value in result['rate_s_inverse'].items())
        assert isclose(doubled['muon_over_electron_rate'],
                       result['muon_over_electron_rate'], rel_tol=2e-14)
    assert phase(0) == 1 and phase(1) == 0

    # COMPARISON ONLY: ordinary HFLAV end-2023 fit (web report 15 Jan 2025).
    # Use the averages WITHOUT the extra unitarity constraint, not Be_univ.
    # https://hflav-eos.web.cern.ch/hflav-eos/tau/end-2023/br-fit.html
    # https://hflav-eos.web.cern.ch/hflav-eos/tau/end-2023/tau-lifetime-avg.html
    tau_tau, sigma_tau_tau = 290.29e-15, 0.53e-15
    branching = {'electron': (0.1784, 0.0004), 'muon': (0.17360, 0.00037)}
    ratio_observed, ratio_sigma = 0.9730, 0.0022  # HFLAV separately rounded fit output
    # Borrowed decay equation, f(x), and inclusive-radiation convention:
    # https://hflav-eos.web.cern.ch/hflav-eos/tau/end-2023/lepton-univ.html
    observed = {}
    for channel, (b, sb) in branching.items():
        rate = b/tau_tau
        # Marginal error only: assume lifetime and branching estimate uncorrelated.
        # Cross-channel correlations are NOT discarded into a spurious joint test.
        observed[channel] = dict(rate=rate,
            sigma=rate*sqrt((sb/b)**2+(sigma_tau_tau/tau_tau)**2))
    comparisons = []
    for name, prediction in predictions.items():
        for channel, rate in prediction['rate_s_inverse'].items():
            obs = observed[channel]
            error = rate/obs['rate'] - 1
            comparisons.append(dict(recipe=name, channel=channel, predicted_rate=rate,
                observed_rate=obs['rate'], marginal_observed_sigma=obs['sigma'],
                fractional_error=error, within_1_percent=abs(error) < 0.01,
                calibration_sigma_only=rate*sigma_tau_mu/tau_mu))
        prediction['ratio_comparison'] = dict(observed=ratio_observed, sigma=ratio_sigma,
            fractional_error=prediction['muon_over_electron_rate']/ratio_observed-1)
    # The comparison adds annotations, not changes to frozen predictions.
    assert all({k: v for k, v in p.items() if k != 'ratio_comparison'} == frozen_predictions[n]
               for n, p in predictions.items())
    assert (_reference, dcu_mass_14) == snapshot

    print('CELL 17 — FROZEN MASS -> DECAY-RATE TRANSFER; ONE MUON CLOCK CALIBRATION')
    for name, recipe in recipes.items():
        print(f"  {name}: r_mu={recipe['mu']}, r_tau={recipe['tau']}; "
              f"r_tau/r_mu={predictions[name]['mass_ratio_tau_mu']:.12f}")
    print(f'Muon lifetime input: {tau_mu:.10e} +/- {sigma_tau_mu:.1e} seconds.')
    print('Imported: universal leading-order weak decay law, including final-lepton phase space.')
    print('No measured lepton masses, tau lifetime, or tau branching values enter prediction.')
    print('\n recipe              channel   predicted(1e12/s)  observed(1e12/s)   error(%)   <1%?')
    for row in comparisons:
        print(f" {row['recipe']:19s} {row['channel']:8s} {row['predicted_rate']/1e12:17.9f}"
              f" {row['observed_rate']/1e12:19.9f} {100*row['fractional_error']:+10.5f}"
              f" {str(row['within_1_percent']):>7s}")
    print('\nPredicted inverse partial rates (ps), NOT total tau lifetime:')
    for name, p in predictions.items():
        print(f"  {name}: e={p['partial_lifetime_ps']['electron']:.9f}, "
              f"mu={p['partial_lifetime_ps']['muon']:.9f}")
    print(f'\nDependent channel-ratio check; HFLAV comparison {ratio_observed} +/- {ratio_sigma}:')
    for name, p in predictions.items():
        print(f"  {name}: Gamma_mu/Gamma_e={p['muon_over_electron_rate']:.9f}; "
              f"error={100*p['ratio_comparison']['fractional_error']:+.5f}%")
    print('The two rate checks share inputs. The ratio is NOT a third independent success.')
    print('Radiative/finite-W corrections are omitted, not fitted; no precision-theory error assigned.')
    print('Retrospective new-observable test, not a blinded discovery or native-clock derivation.')
    print('PASS: exact source arithmetic, clock scaling, two equivalent rate formulas, unchanged inputs.')
    return dict(protocol='frozen structural/screened masses + universal LO weak law',
                primary_recipe='screened_primary', recipes=recipes, predictions=predictions,
                clock_input=dict(value_s=tau_mu, sigma_s=sigma_tau_mu),
                comparisons=comparisons, observational_rates=observed,
                comparison_inputs=dict(tau_tau_s=tau_tau, sigma_tau_tau_s=sigma_tau_tau,
                                       branching_fractions=branching),
                prediction_function=predict, phase_function=phase,
                truly_blinded=False, borrowed_weak_law=True, native_dynamics_changed=False,
                theory_uncertainty_assigned=False, aggregate_significance_computed=False)


dcu_mass_17 = _run_dcu_mass_17()
