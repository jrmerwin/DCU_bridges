# CELL 20 — frozen neutrino inputs -> three-flavor VACUUM survival probabilities.
# Run after Cells 14, 18 and 19 (already present when Cell 19 is complete).
# Standard library only. No network, new fit, native evolution, or detector model.
# Source prescriptions are distinct from the IMPORTED coherent oscillation law.
from copy import deepcopy
from fractions import Fraction
from math import cos, exp, fsum, isclose, isfinite, pi, sin, sqrt
from cmath import exp as cexp


def _run_dcu_mass_20():
    needed = ('_reference', 'dcu_mass_14', 'dcu_mass_18', 'dcu_mass_19')
    missing = [k for k in needed if k not in globals()]
    if missing:
        raise RuntimeError('Run the prerequisite cells first. Missing: ' + ', '.join(missing))
    prior = (_reference, dcu_mass_14, dcu_mass_18, dcu_mass_19)
    snapshot = deepcopy(prior)
    N = len(_reference['rows'])
    prescriptions = dcu_mass_14['prescriptions']
    angles = dcu_mass_19['angle_predictions']
    electron_MeV = dcu_mass_18['calibrations']['electron_rest_energy_MeV']
    if N != 137 or electron_MeV != Fraction('0.51099895069'):
        raise ValueError('The frozen registry or electron energy input changed.')
    if (angles['solar'], angles['reactor']) != (Fraction(167, 548), Fraction(400, 18769)):
        raise ValueError('The frozen screened-angle prescription changed.')
    dimensionless_solar = Fraction(1, 4*N*(N-1)**6)
    if (prescriptions['solar_splitting_over_electron_mass_squared'] != dimensionless_solar
            or prescriptions['neutrino_splitting_ratio'] != 33):
        raise ValueError('The frozen splitting prescription changed. No silent substitution.')
    dm21 = (electron_MeV*10**6)**2 * dimensionless_solar  # eV^2, c=1 mass convention.
    dm32 = 33*dm21
    dm31 = dm21+dm32  # Never confuse the frozen R=dm32/dm21 with dm31/dm21.

    # SI definitions: h, c, elementary charge are exact, NOT fit parameters.
    # NIST: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
    h = Fraction('6.62607015e-34')       # J s
    c = Fraction(299792458)             # m / s
    elementary_charge = Fraction('1.602176634e-19')  # J / eV numerically.
    hbar_c_eVm = float(h*c/elementary_charge)/(2*pi)
    K = 1e-3/(4*hbar_c_eVm)  # Phase: K*dm^2[eV^2]*L[km]/E[MeV].
    assert isclose(K/1000, 1.2669326790391, rel_tol=1e-14)

    def model(s12sq, s13sq, solar, atmospheric31):
        s, t, a, b = map(Fraction, (s12sq, s13sq, solar, atmospheric31))
        if not (0 <= s <= 1 and 0 <= t <= 1 and 0 < a < b):
            raise ValueError('Need valid mixing fractions and normal-ordering positive gaps.')
        w = ((1-s)*(1-t), s*(1-t), t)
        assert sum(w) == 1
        coeff = (4*w[0]*w[1], 4*w[0]*w[2], 4*w[1]*w[2])
        return dict(s12sq=s, s13sq=t, dm21_eV2=a, dm31_eV2=b, dm32_eV2=b-a,
                    weights=w, coefficients=coeff,
                    gaps=(float(a), float(b), float(b-a)),
                    amplitude_float=tuple(map(float, coeff)),
                    atmospheric_amplitude=4*t*(1-t), solar_amplitude=coeff[0],
                    dm_ee_eV2=(1-s)*b+s*(b-a),
                    completely_phase_averaged=sum(p*p for p in w))

    source = model(angles['solar'], angles['reactor'], dm21, dm31)

    def probability(m, y, sigma_y=0.0):
        """P_ee at y=L[km]/E[MeV]. Optional Gaussian y spread is a sensitivity test."""
        y, sigma_y = float(y), float(sigma_y)
        if not (isfinite(y) and isfinite(sigma_y) and y >= 0 and sigma_y >= 0):
            raise ValueError('L/E and its spread must be nonnegative finite numbers.')
        if sigma_y == 0:
            deficit = fsum(A*sin(K*d*y)**2 for A, d in zip(m['amplitude_float'], m['gaps']))
        else:
            # Exact average over Gaussian y with mean y and standard deviation sigma_y.
            # NOT a detector energy kernel, wave-packet decoherence, or thermal noise.
            deficit = fsum(A*(1-cos(2*K*d*y)*exp(-2*(K*d*sigma_y)**2))/2
                           for A, d in zip(m['amplitude_float'], m['gaps']))
        result = 1-deficit
        if not -2e-13 <= result <= 1+2e-13:
            raise ArithmeticError('Probability left [0,1]; no clipping applied.')
        return result

    def at_baseline(m, L_km, E_MeV):
        if not (isfinite(L_km) and isfinite(E_MeV) and L_km >= 0 and E_MeV > 0):
            raise ValueError('Need nonnegative baseline and positive neutrino energy.')
        return probability(m, L_km/E_MeV)

    def amplitude_probability(m, y, common_mass_squared=0.0):
        # The zero in the list is a removable common phase, not a physical m1=0 test.
        masses2 = (0.0, float(m['dm21_eV2']), float(m['dm31_eV2']))
        amplitude = sum(float(w)*cexp(-2j*K*(q+common_mass_squared)*y)
                        for w, q in zip(m['weights'], masses2))
        return abs(amplitude)**2

    # All primary locations/curves fixed independently of the comparison values.
    domains = dict(km_scale=(0.0, 1.0, 2000), multiple_cycles=(0.0, 40.0, 16000))
    locations = ((1.6, 4.0), (52.5, 2.0), (52.5, 3.0),
                 (52.5, 4.0), (52.5, 6.0), (52.5, 8.0))
    period = pi/(K*float(dm21))  # R=33 makes all three phase differences commensurate.
    revival_primary = probability(source, period)
    assert abs(revival_primary-1) < 1e-14
    source_points = [at_baseline(source, L, E) for L, E in locations]
    source_frozen = deepcopy(source)

    # COMPARISON-ONLY curves from the SAME versioned NuFIT6.0 NO variants as Cell19.
    # https://arxiv.org/html/2410.05380v2 Table1; Sept2024 input data.
    # Its normal-ordering atmospheric entry is dm31, NOT dm32.
    fit_dm31 = {'IC24_with_SK_NO': Fraction('0.002513'),
                'IC19_without_SK_NO': Fraction('0.002534')}
    fit_dm21 = Fraction('0.0000749')
    empirical = {}
    for variant, atmospheric31 in fit_dm31.items():
        entries = {r['quantity']: Fraction(str(r['central']))
                   for r in dcu_mass_19['angle_comparisons'] if r['variant'] == variant}
        empirical[variant] = model(entries['solar'], entries['reactor'], fit_dm21, atmospheric31)
    primary_fit_name = 'IC24_with_SK_NO'
    if dcu_mass_19['primary_NuFIT_variant'] != primary_fit_name:
        raise ValueError('Prior primary fit variant changed; no best-match selection.')
    fit = empirical[primary_fit_name]
    # Diagnostic hybrids attribute phase/amplitude differences; never replace primary.
    hybrids = {
        'source_angles_fit_splittings': model(source['s12sq'], source['s13sq'],
                                              fit['dm21_eV2'], fit['dm31_eV2']),
        'fit_angles_source_splittings': model(fit['s12sq'], fit['s13sq'], dm21, dm31)}
    models = dict(source=source, **empirical, **hybrids)

    curves, metrics, grid_checks = {}, {}, {}
    for domain, (a, b, intervals) in domains.items():
        y = tuple(a+(b-a)*j/intervals for j in range(intervals+1))
        values = {name: tuple(probability(m, z) for z in y) for name, m in models.items()}
        curves[domain] = dict(y_km_per_MeV=y, probabilities=values)
        metric = {}
        for name in models:
            if name == primary_fit_name:
                continue
            diff = [p-q for p, q in zip(values[name], values[primary_fit_name])]
            rms = sqrt((fsum(z*z for z in diff)-(diff[0]**2+diff[-1]**2)/2)/intervals)
            j = max(range(len(y)), key=lambda k: abs(diff[k]))
            metric[name] = dict(rms_probability_difference=rms,
                grid_max_absolute_difference=abs(diff[j]), y_at_grid_max=y[j])
        metrics[domain] = metric
        # Double-grid diagnostic: no change in physical parameters or graph inputs.
        yf = tuple(a+(b-a)*j/(2*intervals) for j in range(2*intervals+1))
        df = [probability(source,z)-probability(fit,z) for z in yf]
        finer_rms = sqrt((fsum(z*z for z in df)-(df[0]**2+df[-1]**2)/2)/(2*intervals))
        finer_max = max(map(abs,df))
        grid_checks[domain] = dict(rms_change=abs(finer_rms-metric['source']['rms_probability_difference']),
                                  max_change=abs(finer_max-metric['source']['grid_max_absolute_difference']))
        assert grid_checks[domain]['rms_change'] < 2e-7 and grid_checks[domain]['max_change'] < 3e-6
    points = [dict(L_km=L, E_MeV=E, source=p,
                   empirical_control=at_baseline(fit,L,E),
                   difference=p-at_baseline(fit,L,E))
              for (L,E),p in zip(locations,source_points)]
    revival = dict(y_km_per_MeV=period, source=revival_primary,
        empirical_controls={k:probability(m,period) for k,m in empirical.items()},
        source_Gaussian_y_spread={str(f):probability(source,period,f*period) for f in (0,.01,.03)},
        spread_is_detector_response=False)

    # Additional data-summary check; NOT survival-probability data or independent of NuFIT.
    # Daya Bay final nGd sample, arXiv:2211.14988, normal ordering.
    # Note the paper reports dm32 = 2.466e-3, not dm_ee.
    daya_inputs = dict(sin2_2theta13=(0.0851,0.0024), dm32_eV2=(0.002466,0.000060))
    daya = {}
    for key,pred in (('sin2_2theta13',float(source['atmospheric_amplitude'])),
                     ('dm32_eV2',float(dm32))):
        center,sigma = daya_inputs[key]
        daya[key] = dict(prediction=pred, central=center, sigma=sigma,
                         fractional_difference=pred/center-1,
                         quoted_error_units=(pred-center)/sigma)

    # Independent representations, limits, dimensions, and exact-ratio recurrence.
    checked = 0
    for m in models.values():
        for j in range(1001):
            y = j/25  # 0..40, fixed validation grid.
            assert abs(probability(m,y)-amplitude_probability(m,y)) < 5e-13
            assert abs(amplitude_probability(m,y,.01)-probability(m,y)) < 3e-12
            checked += 1
        assert probability(m,0) == 1
    for y in (0,.125,.5,1.0,7.0,12.0,25.0,40.0):
        assert abs(probability(source,y+period)-probability(source,y)) < 2e-13
        assert at_baseline(source,y*4,4) == at_baseline(source,y*8,8)
    no13 = model(source['s12sq'],0,dm21,dm31)
    for y in (.1,1,10,20):
        expected = 1-4*float(source['s12sq']*(1-source['s12sq']))*sin(K*float(dm21)*y)**2
        assert abs(probability(no13,y)-expected) < 1e-14
    assert probability(model(0,0,dm21,dm31),12.0) == 1
    # A full common period averages every nonzero pair term to 1/2 exactly here.
    average = fsum(probability(source,period*j/4096) for j in range(4096))/4096
    assert abs(average-float(source['completely_phase_averaged'])) < 1e-13
    assert source == source_frozen and prior == snapshot

    print('CELL 20 — FROZEN INPUTS -> COHERENT VACUUM ELECTRON-SURVIVAL CURVES')
    print('Imported: unitary three-flavor relativistic phase propagation; no matter or detector law.')
    print(f"s12^2={source['s12sq']}; s13^2={source['s13sq']}; R=dm32/dm21=33.")
    for key in ('dm21_eV2','dm32_eV2','dm31_eV2','dm_ee_eV2'):
        print(f'  {key:12s} = {float(source[key]):.12e}')
    print(f'Phase constant (km/MeV): {K:.12f}; earlier electron unit reused, no new mass/length fit.')
    print('Electron-row weights:', ', '.join(f'{float(w):.12f}' for w in source['weights']))
    print(f"Solar oscillatory amplitude={float(source['solar_amplitude']):.12f}; "
          f"sin^2(2theta13)={float(source['atmospheric_amplitude']):.12f}.")
    print(f"Complete phase average={float(source['completely_phase_averaged']):.12f} (NOT a solar matter prediction).")
    print('\nExamples: ALL probabilities are calculated; control is NuFIT6.0 IC24+SK NO, NOT counts.')
    print(' L(km) E_nu(MeV)      source P_ee    fit-control P_ee   difference (percentage points)')
    for row in points:
        print(f" {row['L_km']:5.1f} {row['E_MeV']:9.1f} {row['source']:16.9f} "
              f"{row['empirical_control']:19.9f} {100*row['difference']:+24.6f}")
    print('\nFixed-grid differences from the primary FIT CURVE, no probability/energy refit:')
    for domain in domains:
        a,b,n = domains[domain]
        print(f'  L/E in [{a:g},{b:g}] km/MeV; {n+1} points:')
        for name,r in metrics[domain].items():
            print(f"    {name:30s} RMS={100*r['rms_probability_difference']:.6f} pp; "
                  f"max={100*r['grid_max_absolute_difference']:.6f} pp")
    print(f'\nExact-ratio VACUUM REVIVAL at L/E={period:.12f} km/MeV: P_ee=1.')
    for label,p in revival['empirical_controls'].items():
        print(f'  {label}: calculated P_ee at same L/E = {p:.12f}.')
    for f,p in revival['source_Gaussian_y_spread'].items():
        print(f'  Source, Gaussian L/E spread {100*float(f):g}%: mean P_ee={p:.12f}.')
    print('The Gaussian spreads are sensitivity cases, not a fitted detector or decoherence law.')
    print('\nDaya Bay final-result parameter comparison (overlaps NuFIT inputs):')
    for key,r in daya.items():
        print(f"  {key}: source={r['prediction']:.12g}, data-summary={r['central']} +/- {r['sigma']}; "
              f"difference={100*r['fractional_difference']:+.6f}%, quoted-error units={r['quoted_error_units']:+.3f}")
    print(f'PASS: {checked} amplitude/probability checks, common-mass shift, limits, recurrence,')
    print('      unit scaling, exact weight normalization, doubled grids, and unchanged inputs.')
    print('No CP phase/theta23 required for P_ee; no m1=0 measurement, mass-sum, or native-clock test.')
    print('No reactor flux, cross section, multiple-baseline weighting, efficiency, matter, or response folding.')
    print('Fit-curve agreement is NOT an independent event-spectrum validation or joint likelihood.')
    return dict(protocol='frozen screened angles and R=33 splittings + borrowed vacuum oscillations',
        source=source, comparison_models=empirical, diagnostic_hybrids=hybrids,
        phase_constant_km_per_MeV=K, hbar_c_eVm=hbar_c_eVm,
        electron_rest_energy_MeV=electron_MeV, curves=curves, metrics=metrics,
        points=points, revival=revival, daya_bay_summary_comparison=daya,
        numerical_grid_checks=grid_checks, amplitude_checks=checked,
        helpers=dict(probability=probability, at_baseline=at_baseline, model=model,
                     amplitude_probability=amplitude_probability),
        matter_included=False, detector_spectrum_compared=False, CP_phase_assumed=False,
        absolute_m1_tested=False, native_dynamics_changed=False, source_parameters_fitted=False,
        empirical_fit_curves_not_data=True, truly_blinded=False)


dcu_mass_20 = _run_dcu_mass_20()
