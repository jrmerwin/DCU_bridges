# CELL 15 — fixed global prescriptions versus dated external measurements.
# Run after Cell 14. Standard library only; no downloads, fitting, or new dynamics.
# Printed decimal measurements are NOT exact physical constants. Fractions merely
# preserve their quoted central values. Cosmological errors are posterior summaries.
from copy import deepcopy
from decimal import Decimal, localcontext
from fractions import Fraction
from math import sqrt


def _run_dcu_mass_15():
    if 'dcu_mass_14' not in globals():
        raise RuntimeError('Run Cell 14 first. No replacement formulas were selected.')
    old = dcu_mass_14
    before = deepcopy(old)
    p = old['prescriptions']
    N = sum(old['block_sizes'].values())
    I = p['matter_fraction'] * (N-1)
    F = Fraction
    assert N == 137 and I == 40
    assert p['muon_electron'] == F(3*N, 2)
    assert p['solar_splitting_over_electron_mass_squared'] == F(1, 4*N*(N-1)**6)
    assert p['neutrino_splitting_ratio'] == 33

    # Frozen observations, checked 2026-09-19. No source-table predictions reused.
    sources = {
        'CODATA2022': 'https://physics.nist.gov/cuu/Constants/Table/allascii.txt',
        'JUNO59d': 'https://arxiv.org/abs/2511.14593v1',
        'DESI2025': 'https://arxiv.org/html/2503.14738v3#S6',
        'DESI2026': 'https://arxiv.org/html/2607.27410v3#S6.SS3',
    }
    # Natural-unit notation c=1: electron mass expressed in eV.
    me, u_me = F('510998.95069'), F('0.00016')
    observations = {
        'muon_electron': dict(value=F('206.7682827'), error=F('0.0000046'),
                             source='CODATA2022', convention='rest-mass ratio'),
        'matter_fraction': dict(value=F('0.3012'), error=F('0.0079'), source='DESI2026',
                                convention='DR2 galaxy BAO + Lyman-alpha full shape; flat LCDM'),
        'solar_eV2': dict(value=F('7.50e-5'), error=F('0.12e-5'), source='JUNO59d',
                          convention='first 59.1 days; normal ordering; published three-flavor fit'),
    }
    predictions = dict(muon_electron=p['muon_electron'],
                       matter_fraction=p['matter_fraction'],
                       solar_eV2=p['solar_splitting_over_electron_mass_squared'] * me**2)
    solar_u_input = 2*predictions['solar_eV2']*u_me/me

    def comparison(prediction, observation, input_error=F(0)):
        difference = prediction-observation['value']
        # This is a standardized residual, NOT a joint likelihood or tail p-value.
        scale = sqrt(float(observation['error']**2 + input_error**2))
        return dict(prediction=prediction, **observation, difference=difference,
                    relative_difference=difference/observation['value'],
                    input_scale_error=input_error, standardized_residual=float(difference)/scale)

    comparisons = {key: comparison(value, observations[key],
                                    solar_u_input if key == 'solar_eV2' else F(0))
                   for key, value in predictions.items()}
    # Sensitivity rows overlap the primary DESI data: NEVER add them as independent tests.
    sensitivity = []
    for label, value, error, source in (
        ('2025 DR2 BAO only', '0.2975', '0.0086', 'DESI2025'),
        ('2026 Ly-alpha FS only', '0.325', '0.018', 'DESI2026'),
        ('2026 DR2 BAO+Ly-alpha FS+CMB', '0.3042', '0.0033', 'DESI2026'),
    ):
        obs = dict(value=F(value), error=F(error), source=source, convention='flat LCDM')
        sensitivity.append(dict(label=label, **comparison(predictions['matter_fraction'], obs)))

    # Verify algebraic closure on predictions. These are identities, not extra data.
    r, omega = predictions['muon_electron'], predictions['matter_fraction']
    d = predictions['solar_eV2']/me**2
    assert omega*(2*r-3) == 3*I
    assert d*8*r*(2*r-3)**6 == 3**7
    assert d == omega**7/(4*I**6*(omega+I))

    # Eliminate N using the OBSERVED muon ratio. This is a consistency diagnostic,
    # NOT a new fitted registry size or a replacement for the fixed N=137 model.
    rm = observations['muon_electron']['value']
    ur = observations['muon_electron']['error']
    Nm = 2*rm/3
    omega_from_muon = I/(Nm-1)
    solar_from_muon = me**2 * F(3**7, 8)/(rm*(2*rm-3)**6)
    # First-order input propagation; CODATA m_e/r covariance omitted. Its possible
    # effect is negligible beside the external DESI/JUNO errors used here.
    u_omega_mu = abs(2*omega_from_muon/(2*rm-3))*ur
    u_solar_mu = sqrt(float((2*solar_from_muon*u_me/me)**2 +
                           (solar_from_muon*(1/rm+12/(2*rm-3))*ur)**2))
    elimination = {
        'matter_from_muon': comparison(omega_from_muon, observations['matter_fraction'], u_omega_mu),
        'solar_from_muon': comparison(solar_from_muon, observations['solar_eV2'], F(str(u_solar_mu))),
    }

    def infer_solar_N(value):
        # Monotone inverse of the printed formula; numerical bracketing only.
        with localcontext() as ctx:
            ctx.prec = 55
            def dec(f):
                return Decimal(f.numerator)/Decimal(f.denominator)
            target, e = dec(value), dec(me)
            lo, hi = Decimal(2), Decimal(1000)
            def f(x):
                return e*e/(4*x*(x-1)**6)
            assert f(lo) > target > f(hi)
            for _ in range(160):
                mid = (lo+hi)/2
                if f(mid) > target:
                    lo = mid
                else:
                    hi = mid
            return str((lo+hi)/2)

    om, uom = (observations['matter_fraction'][key] for key in ('value', 'error'))
    so, uso = (observations['solar_eV2'][key] for key in ('value', 'error'))
    inferred = {
        'muon': dict(center=Nm, lower=2*(rm-ur)/3, upper=2*(rm+ur)/3),
        'matter': dict(center=1+I/om, lower=1+I/(om+uom), upper=1+I/(om-uom)),
        'solar': dict(center=infer_solar_N(so), lower=infer_solar_N(so+uso), upper=infer_solar_N(so-uso)),
    }
    # These intervals transform the quoted +/-1-error endpoints. They are not
    # posterior fits for N, and the numerical inverse is not another test datum.
    s = predictions['solar_eV2']
    s32, s31 = p['neutrino_splitting_ratio']*s, (p['neutrino_splitting_ratio']+1)*s
    assert s31-s == s32 and s32/s == 33
    masses = (0.0, sqrt(float(s)), sqrt(float(s31)))
    neutrinos = dict(assumptions=('normal ordering', 'm1=0', 'Delta32/Delta21=33'),
                     masses_eV=masses, sum_eV=sum(masses), delta32_eV2=s32, delta31_eV2=s31,
                     absolute_masses_measured_here=False, ratio_33_tested_here=False)
    # Supplementary limits, not independent confirmations or added likelihoods.
    # DESI2026 Table 3 uses a free nonnegative sum and three degenerate species,
    # not a likelihood evaluated at our specific hierarchical spectrum. This is
    # also different from the fixed-0.06-eV matter-fraction rows.
    neutrinos['cosmological_limit_sensitivity'] = [
        dict(model=model, source='DESI2026', upper95_eV=limit,
             hierarchy_assumption='three degenerate species; nonnegative mass-sum prior',
             difference_from_upper_eV=sum(masses)-limit)
        for model, limit in (('LCDM + free sum; DESI+CMB', 0.0592),
                             ('w0waCDM + free sum; DESI+CMB', 0.166))
    ]
    assert old == before
    print('CELL 15 — FIXED GLOBAL FORMULA BRANCH; DATED EXTERNAL COMPARISONS')
    print('N=137; no replacement formula, screening correction, or fitted theory error.')
    print('Primary: CODATA 2022, DESI DR2 BAO+Ly-alpha FS (2026), JUNO first 59.1 days.')
    print('\n quantity            prediction         observed       quoted error    relative gap     standardized')
    for key, row in comparisons.items():
        print(f" {key:18s} {float(row['prediction']):15.10g} {float(row['value']):15.10g} "
              f"{float(row['error']):14.6g} {100*float(row['relative_difference']):+12.6f}% "
              f"{row['standardized_residual']:+15.6f}")
    print('Standardized residuals are summary diagnostics, NOT extreme-tail significances.')
    print('No theoretical truncation error was provided; bare-to-measured corrections remain unspecified.')
    print('\nMatter-fraction sensitivities (overlapping data; flat LCDM, not combined):')
    for row in sensitivity:
        print(f"  {row['label']}: {float(row['value']):.4f} +/- {float(row['error']):.4f}; "
              f"standardized residual={row['standardized_residual']:+.3f}")
    print('\nN-eliminated check using observed muon ratio (not a refit):')
    for key, row in elimination.items():
        print(f"  {key}: {float(row['prediction']):.12g}; "
              f"external standardized residual={row['standardized_residual']:+.6f}")
    print('\nEffective N diagnostics; quoted +/-1-error endpoint transforms only:')
    for key, row in inferred.items():
        print(f"  {key:7s}: {float(row['center']):.9f}; "
              f"[{float(row['lower']):.9f}, {float(row['upper']):.9f}]")
    print('These are not allowed changes to the enumerated registry size.')
    print('\nConditional neutrinos; the splitting match is NOT an absolute-mass measurement:')
    print('  masses (meV):', tuple(round(1000*x, 8) for x in masses))
    print(f"  sum={1000*sum(masses):.8f} meV; Delta32={float(s32):.10g}; Delta31={float(s31):.10g} eV^2")
    for row in neutrinos['cosmological_limit_sensitivity']:
        print(f"  {row['model']}: quoted 95% upper={row['upper95_eV']:.4f} eV; "
              f"sum minus limit={row['difference_from_upper_eV']:+.8f} eV")
    print('Credible upper limits are model/prior dependent, not sharp physical cutoffs.')
    print('\nPASS: old inputs unchanged; exact formula identities; units and splitting conventions checked.')
    print('FINDING: 205.5 is not the measured muon ratio. The uncorrected exact joint branch fails.')
    print('The solar prescription and DESI-only matter prescription are compatible with these summaries.')
    print('No joint p-value, independent-confirmation count, or automatic correction was generated.')
    return dict(protocol='fixed Cell 14 source prescriptions; empirical summary consistency check',
                checked_date='2026-09-19', sources=sources, electron_input_eV=me,
                electron_uncertainty_eV=u_me, predictions=predictions, comparisons=comparisons,
                matter_sensitivity=sensitivity, eliminated_N_checks=elimination,
                effective_N_diagnostics=inferred, conditional_neutrinos=neutrinos,
                exact_joint_branch_matches_observations=False, theory_correction_applied=False,
                native_dynamics_changed=False, previous_calibrations_changed=False,
                joint_likelihood=None)


dcu_mass_15 = _run_dcu_mass_15()
