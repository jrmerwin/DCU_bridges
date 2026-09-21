# Cell 21 — Frozen neutrino parameters through a KamLAND detector reconstruction

## Result and status

Cell 21 evaluates the unchanged Cell 20 source through a documented, approximate reconstruction of the **second KamLAND reactor result (PRL 94, 081801, 2005)**. It is the first event-spectrum comparison in this notebook branch. It is not a reproduction of the collaboration's likelihood, a new fit, a blind prediction, or a native DCU propagation derivation.

The raw-exposure source predicts **215.405 reactor events + 17.790 background events = 233.195 candidates**, compared with 258 published candidates. The same reconstruction using the two frozen NuFIT 6.0 normal-ordering controls predicts 232.492 and 232.817 candidates. Without oscillation it predicts 381.024 candidates. No neutrino parameter, normalization, or background was optimized against those 258 events.

A fixed-nuisance Poisson-deviance diagnostic on 14 fixed prompt bins gives 19.17383 (source), 19.51219 (IC24+SK control), 19.46148 (IC19-noSK control), and 75.81426 (no oscillation). These are not confidence levels or the published KamLAND chi-square. The differences of order 0.3 between the source and global-fit controls are not evidence favoring the source. All parameters, detector conventions, and prior results remain unchanged.

The most important limitation is sensitivity: replacing the source's two fast terms by their phase average changes any primary predicted bin by at most **0.0000451 events**. Even moving all flux in each public distance bin to its center or endpoints keeps the corresponding effect below **0.00159 events per bin**. Adjacent-ratio diagnostics R=32 and R=34 are also virtually indistinguishable. This exposure/resolution/reconstruction therefore does **not** test the proposed exact R=33 revival. It primarily tests the slow oscillation scale/amplitude and the broad disappearance spectrum.

## 1. Frozen theory input

Exactly the source from Cell 20:

- sin²(theta12) = 167/548;
- sin²(theta13) = 400/18769;
- Delta m21² = (510998.95069 eV)² / [4 * 137 * 136^6];
- Delta m32² = 33 Delta m21²;
- Delta m31² = 34 Delta m21²;
- phase coefficient K = 1266.932679039099 for L in km, E in MeV, gaps in eV².

The source assignments come from the screening-paper angle equations and the selected taxonomy splitting branch. Cell 21 does not derive these assignments, change them, choose a CP phase, or test the absolute lightest mass. The primary propagation law remains the borrowed coherent three-flavor vacuum law.

The two control parameter sets are retained from Cell 20's named NuFIT 6.0 (September 2024) variants; no 'best agreement' control selection is made. KamLAND contributes to the global oscillation information, so these are not independent confirmations.

## 2. Why this historical release

The official second-result release provides an event list, distance- and isotope-resolved integrated fission fluxes, visible-energy background templates, target protons, detection efficiency, threshold, and an approximate energy resolution. It is sufficiently open to build an inspectable first forward calculation. Selection was on input availability, not on agreement with the RMR spectrum. It is not claimed to be the latest or highest-statistics reactor result.

Official release:
https://www.awa.tohoku.ac.jp/KamLAND/datarelease/2ndresult.html

Publication:
https://arxiv.org/abs/hep-ex/0406035

The publication's rendered pages 2–4 were inspected, including its selection window, systematic-error table, backgrounds, and spectrum figure. Its analysis is not equivalent to the independent fixed-nuisance diagnostic below: it uses a likelihood fit and background treatments including a floating approximately 6-MeV alpha-n component. Those procedures are not silently replaced by the present deviance.

## 3. Data provenance and transcription

The computation environment could not fetch network files directly. Numerical text was read using web retrieval and transcribed locally; these local data are not advertised as raw byte-for-byte downloaded files or checked against a remote cryptographic checksum.

### Events

https://www.awa.tohoku.ac.jp/KamLAND/datarelease/sort_energy.dat

All **258** sorted prompt energies, as published to 0.01 MeV, are included in the single cell. Prompt energy includes positron annihilation, not just kinetic energy. Count, ordering, endpoints, and ten index-specific values were checked against the displayed source. The observed histogram is:

    [31, 37, 37, 28, 40, 37, 15, 15, 11, 3, 2, 0, 2, 0]

The bins are 14 equal intervals between 2.6 and 8.5 MeV. This is an explicitly chosen analysis binning, not the collaboration's likelihood binning. Both zero bins remain in every score. Released energy rounding is retained; this is not an unbinned detector-level metrology analysis.

### Integrated fission flux

https://www.awa.tohoku.ac.jp/KamLAND/datarelease/fission_flux_distance.dat

The file has forty 25-km intervals from 0 to 1000 km. Sixteen rows are nonzero. All sixteen are transcribed; omitted zero rows contribute nothing. Isotope order is **U235, U238, Pu239, Pu241**. Each entry is integrated fissions/cm² at the detector, already including geometric dilution and livetime. Multiplying again by 1/L² or 515.1 days would be wrong.

The four sums of the table, in fissions/cm², are:

    1.19267676e13, 1.66846658e12, 6.35366737e12, 1.202152355e12

These differ by approximately -1.064%, -1.063%, -1.024%, and -1.000% from the totals in the release's prose. The cell keeps the tabulated values and exposes the discrepancy; it does not adjust isotope weights to force the prose totals.

**Unresolved positions:** the release gives distance intervals, not every reactor's precise baseline in this file. The primary adopts a uniform distribution of *received integrated fission flux* within each interval. It does not assert that reactors are physically uniformly distributed. Center and both endpoint placements are retained as alternative descriptions; these are sensitivity cases, not rigorous confidence limits or exhaustive bounds over all within-bin arrangements.

### Backgrounds

https://www.awa.tohoku.ac.jp/KamLAND/datarelease/BG-Spectrum.dat

The release tabulates arbitrary-unit shapes at 0.01-MeV spacing in prompt energy. To keep the one-cell payload readable, Cell 21 includes every fifth row from 2.60 through 8.50 MeV: **119 knots at 0.05-MeV spacing**. Linear interpolation is integrated analytically between knots and analysis-bin edges. It is normalized inside the selected prompt window to:

- He8/Li9: 4.8 events;
- accidentals: 2.69 events;
- C13(alpha,n): 10.3 events.

These are the publication/release central estimates, not values inferred from fitting the candidate spectrum. The small fast-neutron upper limit contributes no central template here. The publication quotes a total background uncertainty of approximately 7.3 events; no covariance or nuisance fit is reconstructed from that number. Backgrounds are already in visible prompt energy and are not smeared a second time.

A 0.10-MeV thinning check changes an individual background component/bin by at most 0.0371 events. This measures interpolation sensitivity, not a certified error relative to the full 0.01-MeV file or a reduction of the physical background uncertainty.

## 4. Borrowed reactor and detector model

The expected signal in prompt bin b is

    S_b = N_p * efficiency * sum_{r,a} F_{ra}
          * integral dE phi_a(E) sigma_IBD(E) R_b(E) <P_ee(E,L)>_r.

Here F_{ra} is the published integrated fission flux, including source exposure and distance dilution. N_p=4.61e31 and efficiency=0.898 are published detector inputs. The physical baseline and detector energy calibration are borrowed measured quantities, not derived native distances or clocks.

### Fission spectra

Use the pre-2011 ILL-based Huber–Schwetz parameterization from hep-ph/0407026, Appendix A Tables 2/3. The tables were inspected as rendered PDF pages 17/18. For U235, Pu239, Pu241:

    phi_a(E) = exp(sum_{k=0}^5 c_{ak} E^k),

in antineutrinos per fission per MeV. U238 uses the quadratic reproduced in Table 2. All coefficients are in the code. The fits approximate the cited measured/theoretical isotope spectra; no coefficient is fitted to KamLAND here. The original coefficient covariance and fuel/flux uncertainties are not propagated in the current diagnostic. Long-lived fuel corrections and contributions beyond the released 1000-km flux table are not individually modeled. Integration beyond the original fitted high-energy range is an extrapolation of the stated function; the extreme tail has negligible weight in this window, but this does not certify the physical spectral model there.

https://arxiv.org/abs/hep-ph/0407026

### Inverse-beta cross section

Strumia–Vissani, astro-ph/0302055 Eq. (25):

    sigma_IBD(E) = 1e-43 cm² * p_e E_e
                  * E^[-0.07056 + 0.02018 ln(E) - 0.001953 ln(E)^3],
    E_e = E - Delta_np.

All energies are in MeV. This is an external approximate total cross section, not an RMR prediction and not a freely fitted normalization. It includes the effects embodied in that published fit, but no full differential event generator is used.

https://arxiv.org/abs/astro-ph/0302055

CODATA 2022 detector-kinematic inputs:

    m_e c² = 0.51099895069 MeV,
    m_p c² = 938.27208943 MeV,
    (m_n-m_p)c² = 1.29333251 MeV.

https://physics.nist.gov/cuu/Constants/Table/allascii.txt

These laboratory masses describe the borrowed interaction/response model. They do not change the frozen neutrino masses or validate the notebook's own hadron-mass prescriptions. The threshold is [(m_n+m_e)²-m_p²]/(2m_p), expressed through these numbers.

### Energy response

Primary mean prompt energy:

    E_prompt_mean = E - Delta_np + m_e.

A Gaussian with sigma_prompt=0.07 sqrt(E_prompt_mean) MeV is integrated over each prompt bin via its cumulative distribution. This is the official release's **approximate** combined resolution; it is not a full detector response matrix or nonlinearity model. The publication's earlier/later resolution coefficients 0.073 and 0.062 are separately retained, not mixed with assumed time/flux correlations. The energy-scale +/-2% cases are sensitivity tests motivated by the quoted threshold-scale uncertainty, not a claim of a complete energy-dependent covariance.

A separate first-order isotropic mean recoil estimate replaces the mean by

    E_prompt_mean - [2E(E-Delta_np) + Delta_np² - m_e²]/(2m_p).

This is a diagnostic of the leading kinematic simplification, not a full angular/recoil treatment. It never replaces the primary after comparison.

### Baseline averaging

For k=K Delta m²/E and an interval [a,b],

    <sin²(kL)> = [1 - cos(k(a+b)) sinc(k(b-a))]/2,

where sinc(z)=sin(z)/z. This analytic integration avoids an underresolved numerical baseline grid at the fast oscillation frequencies. All three original frequencies are included, unless explicitly running the separate fast-average diagnostic. No post-hoc decoherence or thermal jitter parameter is introduced.

## 5. Primary outcomes

| Model | Reactor signal | Signal + background | Fixed-nuisance Poisson deviance | Conditional shape component |
|---|---:|---:|---:|---:|
| Frozen source | 215.405 | 233.195 | 19.17383 | 16.62421 |
| NuFIT IC24+SK NO control | 214.702 | 232.492 | 19.51219 | 16.81064 |
| NuFIT IC19-noSK NO control | 215.027 | 232.817 | 19.46148 | 16.83065 |
| No oscillation | 363.234 | 381.024 | 75.81426 | 30.95600 |

Observed: 258 candidates. The source count prediction is approximately 9.61% below the observed count. Neither its discrepancy nor its closeness to the global-fit control is called a sub-percent detector success. The source/control predictions differ by at most 0.2831 events in any of these bins. Their conditional normalized-shape total variation is 0.0031141. Those are dependent diagnostics of the same comparison, not independent confirmations.

The reconstructed **no-oscillation reactor signal** 363.234 is 0.5383% below the collaboration's detailed 365.2 +/- 23.7 estimate. This is a normalization cross-check, not validation of the oscillatory shape or of the neutrino theory. A separate control normalizes the no-oscillation calculation to 365.2 (factor 1.00541211); the source then predicts 234.361 candidates. That factor is not estimated from the 258 observations and is not used by the primary.

## 6. Statistical interpretation

For expected counts mu_b and observed counts n_b, retain

    D = 2 sum_b [mu_b - n_b + n_b ln(n_b/mu_b)],

with n ln(n/mu)=0 for n=0. The exact decomposition is

    D = D_rate + D_shape,
    D_rate = 2 [M-N + N ln(N/M)],
    D_shape = 2 sum_b n_b ln(n_b/(N mu_b/M)),

where M=sum mu_b and N=sum n_b. The shape component conditions on the observed total; it does not silently refit the primary reactor normalization. Rate and shape are not independent additional tests to count separately.

No p-value, confidence interval, Delta-chi-square significance, nuisance optimization, or combined significance with previous mass/fit comparisons is reported. Expected bins near zero make casual chi-square/dof interpretation particularly inappropriate. The collaboration's analysis profiles or floats quantities that are fixed here. Covariances, time-dependent fuel/systematics, and exact reactor baselines remain absent.

## 7. Detector and matter sensitivities

All variants were evaluated for every physics model with the same settings. None was selected because it improved the observed spectrum.

For the source, total candidates are:

- uniform within distance bins, primary: 233.195;
- midpoint locations: 233.410;
- lower/upper interval endpoints: 237.248 / 228.264;
- energy-scale -2% / +2%: 229.358 / 236.925;
- resolution coefficients .062 / .073: 233.208 / 233.189;
- constant-density matter control: 234.310;
- leading mean-recoil control: 232.593;
- published no-oscillation normalization control: 234.361.

The matter sensitivity uses the standard three-flavor Hamiltonian in the mass basis:

    H_m²(E) = diag(0,d21,d31) + 2 E V v v^T,
    v_i=sqrt(w_i), V=-sqrt(2) G_F n_e.

The sign is negative for antineutrinos. The constant rock density 2.7 g/cm³ follows the KamLAND later matter-treatment convention; Ye=.5 is explicitly assumed. G_F=1.1663787e-5 GeV^-2 and atomic mass constant 1.66053906892e-24 g are external CODATA values. Diagonalization and baseline averaging are performed without changing source vacuum parameters. This is not a fitted density profile and is not adopted after seeing residuals.

https://arxiv.org/abs/1009.4771 (discussion near Eqs. 2–6)

## 8. What happened to the revival

Cell 20's exact R=33 statement concerned coherent P_ee at one precise L/E. A prompt-energy bin in KamLAND includes multiple reactors, true energies, and detector reconstruction outcomes. It is not a single L/E state.

Under the primary fold, averaging the two atmospheric-frequency sin² terms to 1/2 changes any bin by only 4.50871e-5 events. The sum of absolute bin changes is 1.76033e-4 events. At midpoint baselines (no artificial continuous width per interval), the maximum is 0.000480857; across the midpoint and both endpoint placements it never exceeds 0.00158792 events.

Keeping angles and the solar scale fixed, changing R from 33 to 32 changes the primary bins by at most 0.00012714 events, and R=34 by at most 0.00009463. At midpoint baselines the respective maxima are 0.00107052 and 0.00093870 events. These are explicit identifiability controls, not alternate predictions selected for agreement.

Therefore this dataset/reconstruction is useful for the broad disappearance and slow solar-scale behavior, but cannot substantiate the exact integer relation or the common-phase revival. More accurate arithmetic does not recover information lost through source/detector averaging.

## 9. Execution and validation

- Replayed actual mounted Cells 01, 14, 17, 18, 19, 20 with the original source archive.
- Ran final Cell 21 twice, identical non-callable results.
- Recovered Cell 20 probabilities at 240 baseline/energy/model points independently of the vectorized folding implementation.
- Checked analytic top-hat averages with 80-point Gauss–Legendre integration of the original scalar probability function.
- Halving energy step changes bins by less than 5.7e-14 events; extending 12-MeV upper integration limit to 14 changes bins by less than 2.9e-14 events.
- Twelve independent adaptive scalar integrals agree to at most 1.21e-11 events.
- Nine constant-matter averages independently checked using matrix exponentials: max probability error 1.28e-14.
- Checked background totals, isotope orientation, threshold/response behavior, and exact deviance decomposition.
- Rejected an altered frozen splitting rather than silently accepting a new source.
- Data-development corrections were mechanical: an extra duplicate 5.09 in the first event transcription was caught by the 258-event check; a generator newline-quoting syntax error was fixed before scientific output. Neither changed model choices.

These are numerical implementation checks, not evidence that the approximate reactor flux, response, and nuisance model is experimentally exact.

## 10. Reproduction

Paste **all of DCU_Mass_Cell_21.py as one cell after Cell 20**. NumPy and SciPy must be installed in the notebook kernel. All numerical release excerpts are embedded, so no extra downloads, access keys, ROOT installation, or input-file placement is required. The cell prints the results and retains `dcu_mass_21`; its helper can produce the same documented sensitivity calculations.

The bundle includes the one-cell code, numerical output, JSON results, notes, tests, a detector-spectrum figure, and a separately inspectable transcription of the input tables. No previous model or calibration is overwritten.
