# CELL 22 — frozen atmospheric prediction on Daya Bay's RELEASED Delta-chi2 surface.
# Paste after Cell 20; Cell 21 is not an execution prerequisite. Standard library only.
# The public surface already profiles the collaboration's nuisance parameters.
# It is NOT a new detector simulation or an exact four-source-parameter likelihood.
# The two complete amplitude rows needed by our source are embedded below.
from bisect import bisect_right
from copy import deepcopy
from fractions import Fraction
from math import isclose, isfinite, sin, sqrt


def _run_dcu_mass_22():
    if 'dcu_mass_20' not in globals():
        raise RuntimeError('Run Cell 20 first. No replacement neutrino inputs selected.')
    old = dcu_mass_20
    before = deepcopy(old)
    s0 = old['source']
    s, t, d21 = (s0[k] for k in ('s12sq', 's13sq', 'dm21_eV2'))
    expected_d21 = (Fraction('510998.95069')**2 / (4*137*136**6))
    if (s, t, d21, s0['dm32_eV2']/d21) != (
            Fraction(167,548), Fraction(400,18769), expected_d21, 33):
        raise ValueError('Frozen source changed. No refitting or alternate branch allowed.')
    A = 4*t*(1-t)
    d32 = 33*d21
    dee = d32+(1-s)*d21
    assert A == s0['atmospheric_amplitude'] and dee == s0['dm_ee_eV2']
    probability = old['helpers']['probability']
    make_model = old['helpers']['model']
    K = old['phase_constant_km_per_MeV']
    # Primary prediction is fixed BEFORE consulting the likelihood or comparisons.
    frozen_source = dict(s12sq=s, s13sq=t, dm21_eV2=d21, dm32_eV2=d32,
                         amplitude=A, dm_ee_eV2=dee, R=Fraction(33))

    # Collaboration ancillary file, arXiv:2211.14988v1, header: NO, Delta m32^2.
    grid_url = ('https://arxiv.org/src/2211.14988v1/anc/'
                'DayaBay_DeltaChiSq_NO_3158days.txt')
    # Transcribed numerical excerpts from web text, NOT a byte-for-byte download:
    # 100 consecutive rows per A, delta32 = .002100 + j*.000008, j=0,...,99.
    # These are A=.083000 and A=.083500; original numeric rows 4601..4800.
    # Original web rendered lines 4608..4807 (0-based with 8 header lines).
    values_083 = '''
    40.682 39.024 37.400 35.811 34.255 32.733 31.245 29.792 28.373 26.989
    25.639 24.323 23.043 21.797 20.586 19.410 18.269 17.163 16.093 15.058
    14.058 13.094 12.165 11.272 10.415 9.593 8.808 8.058 7.345 6.667
    6.026 5.422 4.853 4.321 3.826 3.367 2.945 2.560 2.211 1.900
    1.625 1.387 1.187 1.023 0.897 0.808 0.757 0.742 0.765 0.826
    0.924 1.060 1.234 1.445 1.694 1.981 2.305 2.668 3.068 3.506
    3.983 4.498 5.050 5.641 6.270 6.937 7.643 8.387 9.169 9.990
    10.849 11.747 12.683 13.658 14.671 15.723 16.814 17.943 19.111 20.318
    21.564 22.849 24.172 25.535 26.937 28.377 29.857 31.376 32.934 34.531
    36.168 37.844 39.559 41.314 43.108 44.942 46.815 48.728 50.680 52.673
    '''
    values_0835 = '''
    39.718 38.078 36.471 34.899 33.360 31.856 30.385 28.949 27.546 26.178
    24.845 23.546 22.281 21.051 19.856 18.696 17.570 16.479 15.424 14.404
    13.419 12.469 11.555 10.676 9.833 9.025 8.253 7.517 6.817 6.153
    5.524 4.932 4.377 3.857 3.374 2.927 2.517 2.143 1.806 1.505
    1.242 1.015 0.825 0.672 0.556 0.477 0.435 0.430 0.463 0.533
    0.640 0.784 0.966 1.186 1.443 1.738 2.070 2.440 2.848 3.294
    3.777 4.299 4.858 5.455 6.090 6.764 7.475 8.225 9.012 9.838
    10.702 11.605 12.546 13.525 14.543 15.599 16.693 17.826 18.998 20.209
    21.458 22.745 24.072 25.437 26.841 28.284 29.766 31.287 32.847 34.446
    36.084 37.761 39.478 41.234 43.029 44.863 46.737 48.650 50.603 52.596
    '''
    axis_A = tuple(map(Fraction, ('0.083000','0.083500')))
    axis_d = tuple(Fraction('.002100')+j*Fraction('.000008') for j in range(100))
    grid = tuple(tuple(map(Fraction, text.split())) for text in (values_083,values_0835))
    assert all(len(row)==100 for row in grid) and axis_A[0] < A < axis_A[1]
    assert grid[0][38:40] == (Fraction('2.211'),Fraction('1.900'))
    assert grid[0][48:50] == (Fraction('.765'),Fraction('.826'))
    assert grid[1][48:50] == (Fraction('.463'),Fraction('.533'))
    assert grid[1][57:59] == (Fraction('2.440'),Fraction('2.848'))
    assert min(grid[0]) == Fraction('.742') and min(grid[1]) == Fraction('.430')

    def interpolate(amplitude, gap):
        """Exact-rational bilinear interpolation; never extrapolate or reset the minimum."""
        a, z = Fraction(amplitude), Fraction(gap)
        if not axis_A[0] <= a <= axis_A[-1] or not axis_d[0] <= z <= axis_d[-1]:
            raise ValueError('Query outside the explicitly transcribed grid excerpt.')
        j = min(bisect_right(axis_d,z)-1,len(axis_d)-2)
        fa, fz = (a-axis_A[0])/(axis_A[1]-axis_A[0]), (z-axis_d[j])/(axis_d[j+1]-axis_d[j])
        value = ((1-fa)*((1-fz)*grid[0][j]+fz*grid[0][j+1])
                 + fa*((1-fz)*grid[1][j]+fz*grid[1][j+1]))
        corners = tuple((axis_A[i],axis_d[k],grid[i][k]) for i in (0,1) for k in (j,j+1))
        assert min(c[2] for c in corners) <= value <= max(c[2] for c in corners)
        return dict(delta_chi2=value, corners=corners, fractions=(fa,fz))

    for i in range(2):
        for j in range(100):
            assert interpolate(axis_A[i],axis_d[j])['delta_chi2'] == grid[i][j]

    # Additional neighbouring amplitudes for numerical interpolation sensitivity.
    # Four consecutive mass nodes around EACH predetermined R=32,33,34 point.
    extra_nodes = (37,38,39,40,47,48,49,50,56,57,58,59)
    extra = {
        Fraction('.0825'): tuple(map(Fraction, '''
            3.064 2.704 2.382 2.096 1.143 1.157 1.209 1.298 2.630 2.985 3.378 3.809
            '''.split())),
        Fraction('.0840'): tuple(map(Fraction, '''
            1.813 1.488 1.199 .946 .207 .249 .328 .445 1.925 2.303 2.718 3.171
            '''.split()))}
    def lagrange(nodes, values, x):
        total=Fraction(0)
        for i,(xi,yi) in enumerate(zip(nodes,values)):
            weight=Fraction(1)
            for j,xj in enumerate(nodes):
                if j!=i:
                    weight *= (x-xj)/(xi-xj)
            total += weight*yi
        return total

    def cubic_check(R):
        ns = extra_nodes[4*(R-32):4*(R-32)+4]
        nodes = tuple(axis_d[j] for j in ns)
        amplitudes = tuple(map(Fraction,('.0825','.0830','.0835','.0840')))
        results=[]
        for a in amplitudes:
            if a in axis_A:
                values=tuple(grid[axis_A.index(a)][j] for j in ns)
            else:
                values=tuple(extra[a][extra_nodes.index(j)] for j in ns)
            results.append(lagrange(nodes,values,R*d21))
        return lagrange(amplitudes,results,A)

    # R=32 and 34 are FIXED sensitivity alternatives, never substituted for R=33.
    tests={}
    for R in (32,33,34):
        q=interpolate(A,R*d21)
        q.update(R=R, dm32_eV2=R*d21, dm_ee_eV2=(R+1-s)*d21,
                 cubic_interpolation=cubic_check(R))
        q['interpolation_difference']=q['delta_chi2']-q['cubic_interpolation']
        tests[R]=q
    # A physical-region statement from the published TWO-parameter surface,
    # not a posterior probability that the model is correct or a global goodness-of-fit.
    thresholds={'two_parameter_68_3':Fraction('2.30'),
                'two_parameter_95_5':Fraction('6.18')}
    for q in tests.values():
        q['in_conventional_2D_regions']={k:q['delta_chi2'] <= v for k,v in thresholds.items()}
    # Keep the collaboration's origin of Delta chi^2, not the slice minimum.
    curve_R=tuple(Fraction(30)+Fraction(j,100) for j in range(601))
    curve_q=tuple(interpolate(A,R*d21)['delta_chi2'] for R in curve_R)
    jmin=min(range(len(curve_q)),key=curve_q.__getitem__)
    slice_minimum=dict(R_grid=curve_R[jmin],delta_chi2=curve_q[jmin],
                       replaces_primary=False,profiled_over_amplitude=False)

    # The public grid uses its own solar-input constraints, not a fixed scan at our
    # exact source solar parameters. Audit this small matching difference explicitly.
    release_s,release_d=Fraction('.307'),Fraction('.0000753')
    matched_gap=dee-(1-release_s)*release_d
    models={
        'source':s0,
        'release_solar_same_dm32':make_model(release_s,t,release_d,d32+release_d),
        'release_solar_same_dmee':make_model(release_s,t,release_d,matched_gap+release_d)}
    ygrid=tuple(j/2000 for j in range(2401))  # 0..1.2 km/MeV, fixed audit domain.
    differences={name:max(abs(probability(m,y)-probability(s0,y)) for y in ygrid)
                 for name,m in models.items() if name!='source'}
    solar_matching=dict(release_s12_central=release_s,release_dm21_central=release_d,
        same_fast_frequency_dm32=matched_gap,
        same_fast_frequency_surface=interpolate(A,matched_gap),
        max_probability_difference=differences,y_max=1.2,
        four_source_parameters_conditioned_in_public_grid=False)

    # Independent match between three-frequency source and short-baseline d_ee.
    def one_fast_frequency(m,y):
        return (1-float(m['solar_amplitude'])*sin(K*float(m['dm21_eV2'])*y)**2
                  -float(m['atmospheric_amplitude'])*sin(K*float(m['dm_ee_eV2'])*y)**2)
    approximation_error=0.0
    for y in ygrid:
        error=abs(probability(s0,y)-one_fast_frequency(s0,y))
        # Exact Taylor bound because first-order offsets cancel and |(sin^2)''|<=2.
        bound=float(A*s*(1-s))*(K*float(d21)*y)**2
        assert error <= bound+2e-14
        approximation_error=max(approximation_error,error)
    approximation=dict(max_probability_error=approximation_error,
        analytic_uniform_bound=float(A*s*(1-s))*(K*float(d21)*1.2)**2,y_max=1.2)

    # Published spectral-fit summaries. These are NOT newly measured probabilities.
    # Daya Bay nGd two-frequency fit: 2211.14988 page5 Eq2 result; same sample as grid.
    # Daya Bay nH: 2406.01007v2 TableIII, A nominal; B same-data sensitivity.
    # RENO: 2412.18711v1 Eq/main result (stat/systematic errors combined in quadrature).
    summaries={
        'DayaBay_nGd_same_sample':dict(A=(.0852,.0024,.0024),dee=(.002519,.000060,.000060)),
        'DayaBay_nH_A':dict(A=(.0759,.0049,.0050),dee=(.00277,.00015,.00014)),
        'DayaBay_nH_B_same_data':dict(A=(.0776,.0053,.0053),dee=(.00280,.00014,.00014)),
        'RENO_3800days':dict(A=(.0920,sqrt(.0042**2+.0041**2),sqrt(.0044**2+.0041**2)),
            dee=(.00257,sqrt(.00011**2+.00005**2),sqrt(.00010**2+.00005**2)))}
    comparisons={}
    for name,entry in summaries.items():
        out={}
        for key,pred in (('A',float(A)),('dee',float(dee))):
            center,minus,plus=entry[key]
            out[key]=dict(prediction=pred,center=center,error_minus=minus,error_plus=plus,
                fractional_difference=pred/center-1,
                marginal_error_units=(pred-center)/(plus if pred>=center else minus))
        center,minus,plus=entry['dee']
        out['conditional_R']=dict(center=center/float(d21)-1+float(s),
            error_minus=minus/float(d21),error_plus=plus/float(d21),
            solar_inputs_fixed_to_source=True)
        comparisons[name]=out

    # Identifiability: at short baseline, changing R can be compensated by d21.
    # The following models are explicit diagnostics, NOT new source predictions.
    degeneracy={}
    for R in (32,34):
        compensated_d21=dee/(R+1-s)
        compensated=make_model(s,t,compensated_d21,(R+1)*compensated_d21)
        uncompensated=make_model(s,t,d21,(R+1)*d21)
        assert compensated['dm_ee_eV2']==dee
        degeneracy[R]=dict(compensated_dm21=compensated_d21,
            fractional_solar_change=compensated_d21/d21-1,
            max_probability_change_same_dmee=max(abs(probability(compensated,y)-probability(s0,y)) for y in ygrid),
            max_probability_change_fixed_solar=max(abs(probability(uncompensated,y)-probability(s0,y)) for y in ygrid))
    precision=dict(one_integer_frequency_step=d21,
        step_in_DayaBay_nGd_quoted_errors=float(d21)/.000060,
        conditional_three_error_separation_sigma_dee=d21/3,
        corresponding_relative_frequency_precision=(d21/3)/dee,
        independent_solar_precision_and_correlations_needed=True)

    assert old == before
    assert all(isfinite(float(q)) and q>=0 for q in curve_q)
    print('CELL 22 — RELEASED DAYA BAY SPECTRAL DELTA-CHI2, FROZEN SOURCE')
    print(f'A13={float(A):.12f}; dm32={float(d32):.12e}; dmee={float(dee):.12e} eV^2.')
    print('NO surface uses dm32, NOT dmee. 200 grid entries + 24 check entries are embedded.')
    print('Primary interpolation: bilinear; exact rational arithmetic on the published rounded numbers.')
    print('\n R     predicted dm32       released Delta-chi2    cubic check    difference')
    for R,q in tests.items():
        print(f" {R:2d} {float(q['dm32_eV2']):20.12e} {float(q['delta_chi2']):20.9f}"
              f" {float(q['cubic_interpolation']):14.9f} {float(q['interpolation_difference']):+12.9f}")
    print(f'Conventional two-parameter reference levels: {thresholds}. Not model probabilities.')
    print(f"Fixed-amplitude curve minimum on the diagnostic grid: R={float(slice_minimum['R_grid']):.2f};"
          ' R=33 is NOT replaced.')
    print('\nPublic grid retains collaboration solar priors; not the exact 4-parameter source likelihood:')
    print(f"  Same-dmee coordinate matching gives Delta-chi2={float(solar_matching['same_fast_frequency_surface']['delta_chi2']):.9f}.")
    print('  Largest probability differences from changing to the central solar inputs:',differences)
    print('  This probability check is NOT a rigorous bound on the experimental likelihood change.')
    print('\nOther spectral summaries; nGd is SAME dataset, nH B is SAME data as nH A:')
    print(' dataset                       A error units   dmee error units   conditional R +/- errors')
    for name,row in comparisons.items():
        r=row['conditional_R']
        print(f" {name:30s} {row['A']['marginal_error_units']:+10.4f} {row['dee']['marginal_error_units']:+18.4f}"
              f" {r['center']:.5f} -{r['error_minus']:.5f}/+{r['error_plus']:.5f}")
    print('Error units are one-dimensional diagnostics; no invented covariance, summed chi2, or combined significance.')
    print('\nShort-baseline mapping check:',approximation)
    print('One integer step =',float(d21),'eV^2 =',precision['step_in_DayaBay_nGd_quoted_errors'],'nGd frequency errors.')
    for R,row in degeneracy.items():
        print(f"  R={R}: holding dmee fixed needs solar change {100*float(row['fractional_solar_change']):+.4f}%; "
              f"max P change={row['max_probability_change_same_dmee']:.6g} vs {row['max_probability_change_fixed_solar']:.6g} with solar fixed.")
    print('PASS: 200 grid-node replays, no extrapolation, independent interpolation check, exact dmee identity,')
    print('      bounded two-frequency approximation, same-dmee diagnostics, and unchanged Cell20.')
    print('No new detector fit, four-parameter likelihood reconstruction, native-clock derivation, or blind claim.')
    return dict(protocol='fixed source on public DayaBay NO two-parameter profile surface',
        source=frozen_source,grid_source=grid_url,grid_A=axis_A,grid_dm32=axis_d,grid_values=grid,
        interpolation_check_nodes=extra_nodes,interpolation_check_values=extra,
        grid_transcribed_not_downloaded=True,full_grid_retrieved=False,likelihood_tests=tests,
        conventional_2D_thresholds=thresholds,curve=dict(R=curve_R,delta_chi2=curve_q),
        slice_minimum=slice_minimum,solar_matching=solar_matching,
        spectral_summaries=summaries,summary_comparisons=comparisons,
        short_baseline_approximation=approximation,degeneracy_checks=degeneracy,
        sensitivity_precision=precision,helpers=dict(surface=interpolate,one_fast_frequency=one_fast_frequency),
        joint_source_likelihood=False,collaboration_solar_constraints_retained=True,
        solar_reconditioned_to_source=False,
        new_nuisance_fit=False,combined_significance_computed=False,truly_blinded=False,
        native_dynamics_changed=False,source_parameters_changed=False)


dcu_mass_22 = _run_dcu_mass_22()
