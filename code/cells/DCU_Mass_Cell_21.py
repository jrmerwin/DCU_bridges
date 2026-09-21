# CELL 21 — frozen Cell 20 waveform -> KamLAND 2005 reconstructed detector spectrum.
# Run after Cell 20. Requires NumPy and SciPy; no downloads or additional files.
# Historical data, approximate reconstruction, NOT the collaboration likelihood.
# Source neutrino parameters are unchanged. All detector/IBD laws are external.
from copy import deepcopy
from fractions import Fraction
from math import ceil, isclose, sqrt


def _run_dcu_mass_21():
    try:
        import numpy as np
        from scipy.special import ndtr, xlogy
    except ImportError as exc:
        raise ImportError('Cell 21 requires NumPy and SciPy in the notebook kernel.') from exc
    if 'dcu_mass_20' not in globals():
        raise RuntimeError('Run Cell 20 first. No substitute neutrino prescription selected.')
    old = dcu_mass_20
    before = deepcopy(old)
    source = old['source']
    expected_dm21 = Fraction('510998.95069')**2 / (4*137*136**6)
    assert (source['s12sq'], source['s13sq']) == (Fraction(167,548), Fraction(400,18769))
    if (source['dm21_eV2'] != expected_dm21 or source['dm31_eV2'] != 34*expected_dm21
            or source['dm32_eV2'] != 33*expected_dm21):
        raise ValueError('Cell 20 source changed; this frozen comparison is not applicable.')
    K = float(old['phase_constant_km_per_MeV'])

    # Numerical transcriptions of OFFICIAL release text, not raw-byte downloads.
    # Release: https://www.awa.tohoku.ac.jp/KamLAND/datarelease/2ndresult.html
    # Event list: sort_energy.dat (all 258, published to 0.01 MeV).
    # Fissions: fission_flux_distance.dat (all 16 nonzero rows; 24 zero rows omitted).
    # Background: BG-Spectrum.dat, EVERY FIFTH row, 2.60..8.50 MeV (0.05-MeV knots).
    # Backgrounds are already on PROMPT energy scale; do not smear them again.
    events_text = """2.61 2.62 2.62 2.63 2.64 2.64 2.65 2.65 2.66 2.66 2.67 2.67 2.69 2.71 2.72 2.74 2.74 2.74 2.76 2.78 2.80 2.81 2.82 2.84 2.85 2.85 2.88 2.88 2.93 3.00 3.00 3.03 3.05 3.05 3.05 3.05 3.06 3.07 3.08 3.10 3.10 3.11 3.17 3.19 3.20 3.23 3.23 3.26 3.26 3.26 3.28 3.31 3.31 3.32 3.32 3.32 3.33 3.34 3.35 3.35 3.36 3.40 3.40 3.40 3.42 3.42 3.43 3.44 3.45 3.46 3.46 3.46 3.47 3.49 3.49 3.51 3.53 3.53 3.54 3.55 3.55 3.56 3.61 3.67 3.68 3.68 3.69 3.69 3.69 3.70 3.70 3.72 3.74 3.74 3.76 3.76 3.77 3.77 3.78 3.79 3.79 3.81 3.81 3.82 3.84 3.87 3.87 3.88 3.91 3.91 3.95 3.95 3.96 3.96 3.99 4.00 4.02 4.02 4.02 4.03 4.05 4.06 4.07 4.07 4.09 4.12 4.14 4.14 4.17 4.18 4.19 4.20 4.24 4.29 4.29 4.29 4.30 4.33 4.34 4.35 4.35 4.36 4.40 4.40 4.40 4.40 4.47 4.49 4.49 4.49 4.49 4.51 4.51 4.53 4.53 4.56 4.56 4.57 4.58 4.59 4.60 4.60 4.61 4.62 4.63 4.65 4.65 4.65 4.67 4.67 4.68 4.69 4.70 4.71 4.75 4.75 4.75 4.76 4.77 4.77 4.80 4.83 4.83 4.85 4.85 4.85 4.85 4.87 4.89 4.92 4.92 4.93 4.93 4.94 4.95 4.96 4.97 4.97 5.00 5.00 5.03 5.03 5.04 5.08 5.08 5.09 5.09 5.09 5.10 5.11 5.20 5.20 5.22 5.25 5.29 5.30 5.34 5.38 5.40 5.42 5.42 5.46 5.47 5.51 5.51 5.56 5.59 5.60 5.63 5.64 5.65 5.71 5.72 5.75 5.77 5.77 5.78 5.91 5.93 5.97 5.99 6.01 6.09 6.09 6.10 6.13 6.13 6.13 6.15 6.19 6.28 6.40 6.44 6.46 6.86 6.90 7.85 7.95"""

    fission_text = """75 100 7.607167e11 1.010018e11 3.986357e11 7.947312e10
125 150 2.515151e12 3.531436e11 1.309065e12 2.349392e11
150 175 2.636478e12 3.578637e11 1.449826e12 3.017085e11
175 200 3.805156e12 5.524455e11 2.008506e12 3.544056e11
200 225 5.450376e11 7.499991e10 3.084790e11 6.230310e10
275 300 1.824217e11 2.566876e10 1.080335e11 2.298736e10
325 350 5.389535e11 7.212189e10 2.877228e11 5.902264e10
400 425 1.151672e11 1.531188e10 6.054655e10 1.224769e10
425 450 1.548995e11 1.990010e10 7.580203e10 1.303234e10
550 575 9.577854e10 1.359798e10 4.926248e10 8.884929e9
700 725 2.222616e11 3.152956e10 1.132416e11 2.010469e10
725 750 1.013992e11 1.438409e10 5.166192e10 9.172000e9
750 775 9.155953e10 1.359493e10 5.062475e10 9.220732e9
775 800 2.774757e10 3.913515e9 1.415300e10 2.553614e9
825 850 4.134464e10 5.840020e9 2.087991e10 3.712186e9
975 1000 9.269532e10 1.314935e10 4.722713e10 8.384654e9"""

    background_text = """2.60 .472341 8.85411 4.9755
2.65 .482091 10.0101 4.437
2.70 .491893 10.7919 3.96731
2.75 .501663 9.27591 3.54834
2.80 .511368 7.2393 3.16553
2.85 .520991 4.86032 2.80678
2.90 .530523 2.44604 2.46282
2.95 .539959 1.22126 2.12841
3.00 .549295 .642395 1.80297
3.05 .558525 .349435 1.48961
3.10 .567648 .314138 1.19365
3.15 .576658 .285901 .921343
3.20 .585551 .236486 .679402
3.25 .594322 .183542 .474174
3.30 .602969 .134126 .310215
3.35 .611489 .109419 .188489
3.40 .619875 .116478 .105486
3.45 .628127 .0635336 .0539831
3.50 .636239 .0670632 .0251103
3.55 .64421 .060004 .0105654
3.60 .652034 .060004 .00401228
3.65 .65971 .0388261 .00139031
3.70 .667233 .0476502 .00048161
3.75 .674602 .0211779 .000245844
3.80 .681814 .030002 .000261231
3.85 .688865 .0388261 .000376233
3.90 .695752 .0211779 .000541198
3.95 .702474 .0141186 .000764689
4.00 .709027 .00705929 .00117707
4.05 .715408 .0105889 .00224942
4.10 .721616 .0123538 .0052975
4.15 .727648 .00352964 .0134677
4.20 .733502 .0105889 .0333577
4.25 .739176 1.62265e-14 .0771117
4.30 .744667 .0123538 .164138
4.35 .749973 .0105889 .320743
4.40 .755092 .00176482 .575616
4.45 .760023 .00705929 .950128
4.50 .764764 .00529447 1.44521
4.55 .769312 9.24184e-15 2.03001
4.60 .773664 .00705929 2.63914
4.65 .777823 .00352964 3.18331
4.70 .781783 .00176482 3.57169
4.75 .785544 0 3.73828
4.80 .789106 .00352964 3.66103
4.85 .792465 .00352964 3.36597
4.90 .795622 .00529447 2.91568
4.95 .798574 0 2.38844
5.00 .80132 .00529447 1.85731
5.05 .803862 .00705929 1.37614
5.10 .806195 .00529447 .974887
5.15 .80832 2.2998e-14 .662411
5.20 .810237 0 .432983
5.25 .811944 .00705929 .273368
5.30 .813441 .00352964 .168608
5.35 .814726 .0105889 .106323
5.40 .815801 .00176482 .071977
5.45 .816664 .00352964 .0852984
5.50 .817314 .00529447 .154246
5.55 .817753 .00352964 .32004
5.60 .81798 .00176482 .65225
5.65 .817993 1.338e-14 1.24857
5.70 .817795 .00176482 2.22053
5.75 .817383 .00352964 3.65962
5.80 .81676 .00176482 5.58659
5.85 .815925 .0105889 7.90063
5.90 .814877 .00705929 10.3553
5.95 .813619 .00705929 12.5869
6.00 .81215 .00705929 14.1997
6.05 .81047 .00352964 14.8828
6.10 .808581 .00529447 14.5106
6.15 .806483 .00352964 13.1808
6.20 .804178 .00352964 11.1747
6.25 .801665 .00705929 8.86039
6.30 .798945 .00352964 6.58508
6.35 .796021 .00352964 4.5977
6.40 .792893 .00529447 3.02216
6.45 .789562 3.27791e-14 1.87351
6.50 .786029 0 1.09665
6.55 .782296 .00352964 .606379
6.60 .778364 .00705929 .31661
6.65 .774236 .00352964 .155932
6.70 .769911 .00352964 .0723197
6.75 .765394 .00705929 .0315225
6.80 .760683 .00529447 .0128858
6.85 .755783 .00352964 .00492979
6.90 .750695 .00529447 .00176181
6.95 .74542 1.82705e-14 .000587209
7.00 .739962 .00176482 .000182277
7.05 .734321 1.86467e-14 .0000526377
7.10 .728501 .00352964 .0000141285
7.15 .722504 .00352964 .0000035004
7.20 .716332 .00529447 .000000804223
7.25 .709989 0 .000000173338
7.30 .703475 0 .0000000346252
7.35 .696794 .0105889 .00000000207755
7.40 .689949 .00352964 0
7.45 .682944 2.01515e-14 0
7.50 .67578 0 0
7.55 .668461 .0141186 0
7.60 .660991 .00176482 0
7.65 .653372 .00352964 0
7.70 .645606 .00882411 0
7.75 .637699 2.12801e-14 0
7.80 .629653 .00352964 0
7.85 .621473 .0105889 0
7.90 .613161 .00529447 0
7.95 .604721 .0105889 0
8.00 .596158 .00176482 0
8.05 .587474 .00352964 0
8.10 .578675 .00176482 0
8.15 .569765 2.27848e-14 0
8.20 .560746 .00176482 0
8.25 .551623 .00352964 0
8.30 .542403 .00352964 0
8.35 .533088 2.35372e-14 0
8.40 .523682 .00529447 0
8.45 .514192 .0105889 0
8.50 .504621 .00529447 0"""

    events = np.fromstring(events_text, sep=' ')
    flux = np.fromstring(fission_text, sep=' ').reshape(-1,6)
    bg_table = np.fromstring(background_text, sep=' ').reshape(-1,4)
    assert events.shape == (258,) and np.all(np.diff(events)>=0)
    assert (events[0], events[-1]) == (2.61,7.95)
    assert [(i,events[i]) for i in (29,50,99,129,159,199,208,230,250,257)] == [
        (29,3.00),(50,3.28),(99,3.79),(129,4.18),(159,4.59),
        (199,5.00),(208,5.10),(230,5.65),(250,6.28),(257,7.95)]
    assert flux.shape==(16,6) and np.all(flux[:,1]-flux[:,0]==25)
    assert bg_table.shape==(119,4) and np.allclose(np.diff(bg_table[:,0]),.05,atol=1e-14)
    # The release's rounded prose totals differ ~1% from the actual table sums.
    # Preserve the table. Never silently scale its individual isotope columns.
    prose_flux = np.array([1.2055e13,1.6864e12,6.4194e12,1.2143e12])
    flux_discrepancy = flux[:,2:].sum(axis=0)/prose_flux-1

    # Analysis convention fixed for every model: 14 equal prompt-energy bins.
    # This is NOT the publication's likelihood binning; no empty bins are dropped.
    edges = np.linspace(2.6,8.5,15)
    bg_totals = np.array([4.8,2.69,10.3])  # He8/Li9, accidentals, C13(alpha,n).
    protons, efficiency = 4.61e31, .898
    published_noosc, published_norm_sigma = 365.2, 23.7
    # IBD detector constants, measured laboratory values, NOT source predictions.
    # CODATA2022 https://physics.nist.gov/cuu/Constants/Table/allascii.txt
    me, mp, delta_np = .51099895069, 938.27208943, 1.29333251
    # Huber-Schwetz hep-ph/0407026 Tables 2/3, coefficient order E^0,... .
    # Pre-2011 ILL parameterization for this historical reconstruction.
    # Units: antineutrinos / fission / MeV. U238 is the original quadratic.
    isotopes = ('U235','U238','Pu239','Pu241')
    coefficients = ((3.519,-3.517,1.595,-.4171,.05004,-.002303),
                    (.976,-.162,-.0790),
                    (2.560,-2.654,1.256,-.3617,.04547,-.002143),
                    (1.487,-1.038,.4130,-.1423,.01866,-.0009229))

    def background_bins(coarse=False):
        # Piecewise-linear template integrals normalized IN the analysis window.
        # 0.10-MeV thinning is retained as an interpolation sensitivity check.
        table = bg_table[::2] if coarse else bg_table
        x = np.unique(np.r_[table[:,0],edges])
        y = np.array([np.interp(x,table[:,0],table[:,j]) for j in (1,2,3)])
        areas = (y[:,1:]+y[:,:-1])*np.diff(x)/2
        primitive = np.c_[np.zeros(3),np.cumsum(areas,axis=1)]
        bins = np.diff(primitive[:,np.searchsorted(x,edges)],axis=1)
        bins = bins/bins.sum(axis=1)[:,None]*bg_totals[:,None]
        assert np.all(bins>=0) and np.allclose(bins.sum(axis=1),bg_totals,atol=1e-12)
        return bins

    def terms(model,E,rho=0.0):
        s,t = float(model['s12sq']),float(model['s13sq'])
        a,b = float(model['dm21_eV2']),float(model['dm31_eV2'])
        w = np.array([(1-s)*(1-t),s*(1-t),t])
        if rho == 0:
            W = np.tile(w,(len(E),1))
            gaps = np.tile(np.array([a,b,b-a])[:,None],(1,len(E)))
        else:
            # Sensitivity ONLY: standard constant-density three-flavor matter law.
            # KamLAND arXiv1009.4771 uses rho=2.7g/cm^3; Ye=.5 is declared here.
            # V=-sqrt(2) G_F n_e for antineutrinos. rho/amu in cm^-3 -> m^-3.
            GF, amu_g, Ye = 1.1663787e-5, 1.66053906892e-24, .5
            V = -sqrt(2)*GF*1e-18*(rho*Ye/amu_g*1e6)*old['hbar_c_eVm']**3
            H = np.tile(np.diag([0,a,b]),(len(E),1,1))
            H += 2*E[:,None,None]*1e6*V*np.outer(np.sqrt(w),np.sqrt(w))
            ev, U = np.linalg.eigh(H)
            W = np.einsum('eij,i->ej',U,np.sqrt(w))**2
            assert np.allclose(W.sum(axis=1),1,atol=2e-14)
            gaps = np.array([ev[:,1]-ev[:,0],ev[:,2]-ev[:,0],ev[:,2]-ev[:,1]])
        A = np.array([4*W[:,0]*W[:,1],4*W[:,0]*W[:,2],4*W[:,1]*W[:,2]])
        return A, gaps

    def mean_probability(model,E,positions='uniform',rho=0.0,fast_average=False):
        if model is None:
            return np.ones((len(flux),len(E)))
        A,gaps = terms(model,E,rho)
        k = K*gaps/E
        lo,hi = flux[:,0,None,None],flux[:,1,None,None]
        L = {'low':lo,'high':hi}.get(positions,(lo+hi)/2)
        if positions not in ('uniform','midpoint','low','high'):
            raise ValueError('Unknown unresolved-baseline convention.')
        cosine = np.cos(2*k[None,:,:]*L)
        if positions == 'uniform':
            # Exact top-hat integral of cos(2kL); ASSUMPTION within 25-km bins.
            # The release does not give individual reactor locations in these bins.
            cosine *= np.sinc(k[None,:,:]*(hi-lo)/np.pi)
        if fast_average:
            cosine[:,1:,:] = 0  # Diagnostic, NOT a replacement fitted to data.
        p = 1-np.sum(A[None,:,:]*(1-cosine)/2,axis=1)
        if np.min(p)<-2e-13 or np.max(p)>1+2e-13:
            raise ArithmeticError('Probability outside [0,1]; no clipping applied.')
        return p

    def forward(model, step=.002, positions='uniform',resolution=.07,
                energy_scale=1.0,rho=0.0,fast_average=False,recoil=False,
                top_energy=12.0,return_kernel=False):
        if step<=0 or resolution<=0 or energy_scale<=0 or rho<0 or top_energy<10:
            raise ValueError('Invalid forward-integration setting.')
        threshold = delta_np+me+(delta_np+me)**2/(2*mp)
        n = ceil((top_energy-threshold)/step)
        dE = (top_energy-threshold)/n
        E = threshold+(np.arange(n)+.5)*dE
        spectra = np.array([np.exp(np.polynomial.polynomial.polyval(E,c)) for c in coefficients])
        Ee = E-delta_np
        # Strumia-Vissani astro-ph/0302055 Eq25. sigma in cm^2, energies in MeV.
        logE = np.log(E)
        sigma = 1e-43*Ee*np.sqrt(Ee*Ee-me*me)*E**(-.07056+.02018*logE-.001953*logE**3)
        prompt = Ee+me  # LO mean prompt kinematics; NOT E_prompt=E_nu.
        if recoil:
            # Independently declared first-order isotropic mean recoil estimate.
            # Not a full differential IBD generator; never replaces the primary.
            prompt -= (2*E*Ee+delta_np**2-me**2)/(2*mp)
        response = np.diff(ndtr((edges[:,None]-energy_scale*prompt)/
                                (resolution*np.sqrt(prompt))),axis=0)
        assert np.all(response>=0) and np.max(response.sum(axis=0))<=1+1e-14
        p = mean_probability(model,E,positions,rho,fast_average)
        # Flux ALREADY includes geometric dilution and total livetime.
        # Do NOT apply another 1/L^2 or multiply again by 515.1 days.
        flux_per_bin_energy = flux[:,2:]@spectra
        spectrum = np.sum(flux_per_bin_energy*p,axis=0)*sigma*protons*efficiency
        bins = response@(spectrum*dE)
        assert np.all(bins>=0) and np.all(np.isfinite(bins))
        if return_kernel:
            return bins,dict(E=E,dE=dE,response=response,unfolded_density=spectrum,
                             baseline_survival=p,flux=flux_per_bin_energy,cross_section=sigma)
        return bins

    models = dict(source=source,**old['comparison_models'],no_oscillation=None)
    fixed_signals = {name:forward(m) for name,m in models.items()}
    backgrounds = background_bins()
    bg = backgrounds.sum(axis=0)
    predictions = {name:dict(signal=s,background=bg.copy(),total=s+bg)
                   for name,s in fixed_signals.items()}
    frozen_predictions = deepcopy(predictions)
    primary_name = 'source'
    # Independent detector conventions, retained for ALL physics models alike.
    # These are sensitivity cases, never optimized/selected to repair residuals.
    variants = {
        'bin_midpoints':dict(positions='midpoint'),
        'bin_lower_endpoints':dict(positions='low'),
        'bin_upper_endpoints':dict(positions='high'),
        'energy_scale_minus2percent':dict(energy_scale=.98),
        'energy_scale_plus2percent':dict(energy_scale=1.02),
        'resolution_6p2':dict(resolution=.062),
        'resolution_7p3':dict(resolution=.073),
        'constant_rock_2p7':dict(rho=2.7),
        'first_order_mean_recoil':dict(recoil=True)}
    sensitivities = {tag:{name:forward(m,**kw)+bg for name,m in models.items()}
                     for tag,kw in variants.items()}
    # Published detailed no-oscillation expectation is an EXTERNAL normalization
    # control, not the 258 observed events. Primary retains raw published exposure.
    noosc_raw = fixed_signals['no_oscillation'].sum()
    norm_to_release = published_noosc/noosc_raw
    release_normalized = {name:s*norm_to_release+bg for name,s in fixed_signals.items()}
    averaged_fast = forward(source,fast_average=True)
    # Sensitivity, not alternate fits: does this release resolve the integer R at all?
    # Keep solar scale and angles fixed; evaluate adjacent ratios 32 and 34.
    ratio_resolution = {}
    for ratio in (32,34):
        alternative = dict(source,dm31_eV2=(ratio+1)*source['dm21_eV2'],
                           dm32_eV2=ratio*source['dm21_eV2'])
        ratio_resolution[ratio] = {}
        for position in ('uniform','midpoint'):
            delta = forward(alternative,positions=position)-forward(source,positions=position)
            ratio_resolution[ratio][position] = dict(
                max_absolute_bin_events=float(np.max(abs(delta))),
                total_absolute_bin_events=float(np.sum(abs(delta))))
    fast_position_checks = {}
    for position in ('uniform','midpoint','low','high'):
        difference = forward(source,positions=position,fast_average=True)-forward(source,positions=position)
        fast_position_checks[position] = dict(max_absolute_bin_events=float(np.max(abs(difference))),
                                            total_absolute_bin_events=float(np.sum(abs(difference))))

    # Only NOW form observed bin counts and comparison diagnostics.
    observed = np.histogram(events,edges)[0]
    assert observed.sum()==258 and len(observed)==14
    def score(expected):
        if np.any(expected<=0):
            raise ValueError('Deviance requires positive expected counts in every bin.')
        full = 2*np.sum(expected-observed+xlogy(observed,observed/expected))
        obs_sum, exp_sum = observed.sum(),expected.sum()
        rate = 2*(exp_sum-obs_sum+obs_sum*np.log(obs_sum/exp_sum))
        # Conditional multinomial shape diagnostic; NOT a fitted reactor flux.
        shape = 2*np.sum(xlogy(observed,observed/(obs_sum*expected/exp_sum)))
        assert abs(full-rate-shape)<2e-12
        return dict(total_expected=float(exp_sum),fractional_count_difference=float(exp_sum/obs_sum-1),
                    poisson_deviance=float(full),rate_component=float(rate),shape_component=float(shape))
    summaries = {name:score(p['total']) for name,p in predictions.items()}
    sensitivity_scores = {tag:{name:score(x) for name,x in p.items()} for tag,p in sensitivities.items()}
    normalized_scores = {name:score(x) for name,x in release_normalized.items()}
    # Rough released total-error scale is reported, NOT turned into a joint pull.
    # Norm/systematics/background covariance and reactor-bin uncertainty are unprofiled.

    # Validation: independent point formula, exact baseline integral, quadrature.
    Es = np.array([3.0,3.5,5.,8.,10.])
    for name,m in models.items():
        if m is None:
            continue
        centers = flux[:,:2].mean(axis=1)
        p = mean_probability(m,Es,'midpoint')
        refp = np.array([[old['helpers']['probability'](m,L/E) for E in Es] for L in centers])
        assert np.max(abs(p-refp))<2e-13
        nodes,weights = np.polynomial.legendre.leggauss(80)
        avg = np.array([[sum(w*old['helpers']['probability'](m,(a+(z+1)*(b-a)/2)/E)
                           for z,w in zip(nodes,weights))/2
                         for E in Es] for a,b in flux[:,:2]])
        assert np.max(abs(avg-mean_probability(m,Es)))<2e-12
    convergence = {}
    for name,m in models.items():
        finer = forward(m,step=.001)
        err = float(np.max(abs(finer-fixed_signals[name])))
        assert err<1e-6
        convergence[name] = err
    coarse_bg_error = float(np.max(abs(background_bins(True)-backgrounds)))
    tail_check = float(np.max(abs(forward(source,top_energy=14)-fixed_signals['source'])))
    assert tail_check<1e-6
    for name,p in predictions.items():
        assert all(np.array_equal(p[k],frozen_predictions[name][k]) for k in p)
    assert old == before
    fast_difference = dict(
        max_absolute_bin_events=float(np.max(abs(averaged_fast-fixed_signals['source']))),
        total_absolute_bin_events=float(np.sum(abs(averaged_fast-fixed_signals['source']))),
        total_rate_difference=float(averaged_fast.sum()-fixed_signals['source'].sum()))

    print('CELL 21 — KAMLAND 2005 DETECTOR-FOLDED TRANSFER, NO SOURCE FIT')
    print('Data: 258 actual prompt energies; 16 nonzero 25-km fission-flux bins; 3 background templates.')
    print('14 fixed bins, 2.6--8.5 MeV. Efficiency .898; 4.61e31 target protons.')
    print('Primary: uniform unresolved positions within each published baseline bin, vacuum, 7%/sqrt(E).')
    print('External isotope spectra and IBD law; backgrounds normalized to 4.8, 2.69, 10.3 events.')
    print(f'No-oscillation raw signal={noosc_raw:.6f}; release detailed expectation=365.2 +/- 23.7.')
    print('Fission-table/prose-total differences (%), retained:',100*flux_discrepancy)
    print('\n model                       signal   signal+BG  observed     D_Poisson   D_shape')
    for name,r in summaries.items():
        print(f" {name:27s} {fixed_signals[name].sum():8.3f} {r['total_expected']:11.3f}"
              f" {observed.sum():9d} {r['poisson_deviance']:13.5f} {r['shape_component']:10.5f}")
    print('Deviances are fixed-nuisance diagnostics, NOT official chi2/p-values or independent tests.')
    print('\nObserved bin counts:',observed.tolist())
    print('Source expected bins:',np.round(predictions['source']['total'],4).tolist())
    print('\nSource sensitivity cases (all also evaluated for controls; none selected as a fit):')
    for tag,r in sensitivity_scores.items():
        print(f"  {tag:30s}: N={r['source']['total_expected']:.6f}; D={r['source']['poisson_deviance']:.6f}")
    print(f"  release_noosc_normalization: N={normalized_scores['source']['total_expected']:.6f}; "
          f"D={normalized_scores['source']['poisson_deviance']:.6f}; scale={norm_to_release:.9f}")
    print('\nReplacing ONLY fast 31/32 interference by its phase average:')
    print('  Source max bin change:',fast_difference['max_absolute_bin_events'],'events')
    print('  Largest corresponding change in midpoint/lower/upper baseline checks:',
          max(x['max_absolute_bin_events'] for x in fast_position_checks.values()),'events')
    print('  Adjacent-R diagnostics, largest bin change (uniform / midpoint):')
    for ratio, by_position in ratio_resolution.items():
        print('    R=',ratio,by_position['uniform']['max_absolute_bin_events'],
              by_position['midpoint']['max_absolute_bin_events'])
    print('  This reconstructed spectrum does not test the exact R=33 revival.')
    print('Validation: max half-step bin change',max(convergence.values()),'events;')
    print('  background 0.05-to-0.10-MeV thinning max component/bin change',coarse_bg_error,'events;')
    print('  E_nu upper-limit 12-to-14-MeV change',tail_check,'events.')
    print('PASS: prior waveform checks, analytic/numerical baseline integrals, counts/flux orientation,')
    print('      integration refinement, background totals, deviance decomposition, and unchanged Cell 20.')
    print('Limitations: coarse baseline bins, mean prompt mapping, one Gaussian resolution, fixed background')
    print('and flux nuisances; no official detector response matrix or collaboration likelihood reproduced.')
    print('Historical dataset overlaps global fits; no blind validation or native propagation derivation claimed.')
    return dict(protocol='KamLAND 2nd result (2005): fixed vacuum source, approximate forward detector fold',
        primary_model=primary_name,source_parameters={k:source[k] for k in ('s12sq','s13sq','dm21_eV2','dm31_eV2','dm32_eV2')},
        bin_edges_prompt_MeV=edges.tolist(),observed_counts=observed.tolist(),event_prompt_MeV=events.tolist(),
        predictions={name:{k:v.tolist() for k,v in p.items()} for name,p in predictions.items()},
        summaries=summaries,background_components=backgrounds.tolist(),background_totals=bg_totals.tolist(),
        nonzero_fission_flux_bins=flux.tolist(),fission_prose_discrepancy=flux_discrepancy.tolist(),
        isotope_order=isotopes,isotope_coefficients=coefficients,background_knots=bg_table.tolist(),
        detector=dict(target_protons=protons,efficiency=efficiency,resolution_coefficient=.07,
                      primary_within_bin='uniform',primary_matter=False,primary_recoil=False,
                      noosc_release=published_noosc,noosc_release_sigma=published_norm_sigma,
                      noosc_raw=float(noosc_raw),normalization_to_release_control=float(norm_to_release)),
        sensitivity_parameters=variants,sensitivity_scores=sensitivity_scores,
        sensitivities={tag:{name:x.tolist() for name,x in p.items()} for tag,p in sensitivities.items()},
        release_normalized_control={name:x.tolist() for name,x in release_normalized.items()},
        release_normalized_scores=normalized_scores,fast_phase_average_diagnostic=fast_difference,
        fast_phase_averaged_source_signal=averaged_fast.tolist(),
        fast_baseline_position_checks=fast_position_checks,adjacent_ratio_sensitivity=ratio_resolution,
        numerical_checks=dict(energy_half_step=convergence,background_thinning=coarse_bg_error,
                              upper_energy_limit=tail_check),
        empirical_p_value_computed=False,physics_parameters_fitted=False,nuisance_parameters_fitted=False,
        approximate_detector_reconstruction=True,collaboration_likelihood_reproduced=False,
        native_dynamics_changed=False,truly_blinded=False,
        source_urls=dict(release='https://www.awa.tohoku.ac.jp/KamLAND/datarelease/2ndresult.html',
                         paper='https://arxiv.org/abs/hep-ex/0406035',
                         flux_spectra='https://arxiv.org/abs/hep-ph/0407026',
                         cross_section='https://arxiv.org/abs/astro-ph/0302055'),
        helper_forward=forward)


dcu_mass_21 = _run_dcu_mass_21()
