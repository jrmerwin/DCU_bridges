# Cell 15 — fixed global prescriptions: empirical consistency and scope

## What was executed

Cell 1 was replayed from the originally attached reproducibility ZIP, then Cell 14
was replayed and its entire scientific output matched the archived Cell 14 output.
Cell 15 was executed twice. Both outputs and result dictionaries agreed. A deliberate
mutation of the muon formula was rejected by the source-prescription assertions.
No local-geometry optimization, native service history, new coupling simulation,
parameter fitting, or physical recalibration was performed.

## Frozen formula branch

This is the explicitly selected cross-manuscript subset introduced in Cell 14:

- March lepton manuscript Eq. (9): r_mu = 3 N / 2.
- Unpublished taxonomy: Omega_m = I/(N-1), I=40.
- Unpublished taxonomy: Delta21/m_e^2 = 1/[4 N (N-1)^6].
- Conditional extension: R = Delta32/Delta21 = 33; normal ordering; m1=0.

N=137 is the enumerated count. The one-interface exclusion, the physical
assignments, and the depth exponent six are retained hypotheses, not newly
derived by this cell. The alternative pi-based muon expression and screening
formulas were NOT substituted. The source's printed rounded decimal predictions
were NOT used as inputs. No theoretical error bar was introduced to make the
predictions agree with data.

The taxonomy explicitly leaves radiative corrections and formula uniqueness open.
Consequently, a failure as an exact physical equality is not a theorem against an
unspecified, corrected theory. Conversely, calling a discrepancy a radiative
correction does not calculate that correction or validate a bare-to-measured map.

## External data frozen on 2026-09-19

1. NIST CODATA 2022 table, muon-electron mass ratio and electron mass energy:
   https://physics.nist.gov/cuu/Constants/Table/allascii.txt
   - r_mu = 206.7682827 +/- 0.0000046.
   - m_e c^2 = 0.51099895069 +/- 0.00000000016 MeV.
   - Natural-unit conversion used: m_e = 510998.95069 +/- 0.00016 eV.

2. JUNO Collaboration, arXiv:2511.14593v1, abstract and Results, first 59.1 days:
   https://arxiv.org/abs/2511.14593v1
   https://arxiv.org/html/2511.14593v1
   - Delta21 = (7.50 +/- 0.12) x 10^-5 eV^2, normal ordering.
   - This is JUNO's published three-flavor result, not a reanalysis. The analysis
     uses Daya Bay spectral information and external oscillation constraints.
   - No PDG/NuFIT/global-fit estimate containing the same oscillation data was added
     as an independent solar observation.

3. DESI DR2 Results IV, arXiv:2607.27410v3, Eq. (28), Table 3, sections VI.3/VI.4.2:
   https://arxiv.org/html/2607.27410v3
   - PRIMARY cosmology comparison: galaxy BAO + Lyman-alpha forest full shape,
     Omega_m = 0.3012 +/- 0.0079, flat LCDM, fixed mass sum 0.06 eV.
   - Sensitivity: Lyman-alpha forest full shape alone, 0.325 +/- 0.018.
   - Sensitivity: DESI + CMB, 0.3042 +/- 0.0033 (fixed mass sum).
   - These are posterior parameter summaries conditional on the stated cosmology,
     not model-independent direct measurements of a count fraction.

4. DESI DR2 Results II, arXiv:2503.14738v3, Eq. (17), Table 5:
   https://arxiv.org/html/2503.14738v3
   - Older BAO-only sensitivity comparison: 0.2975 +/- 0.0086 in flat LCDM.
   - These data overlap the 2026 DESI combination. They are not separate repeated
     confirmations and are not combined into a likelihood or average.

The primary cosmology row uses the 2026 DESI-only combination rather than keeping
the older, numerically closer BAO-only value as the sole comparison.

## Exact algebra and computed predictions

Using Fractions preserves the entered decimal central values; it does not make
measured constants exact. The predictions at N=137 are:

- r_mu = 205.5.
- Omega_m = 5/17 = 0.2941176470588235...
- Delta21 = 7.530538062369455 x 10^-5 eV^2.

The input electron uncertainty propagates as u(Delta21)=2 Delta21 u(m_e)/m_e,
which is negligible beside the JUNO uncertainty. Differences between the source's
rounded 7.531 x 10^-5 and the recomputed number were not silently treated as new
predictions or a formula correction.

Exact eliminated identities:

    Omega_m (2 r_mu - 3) = 3 I
    (Delta21/m_e^2) 8 r_mu (2 r_mu - 3)^6 = 3^7
    Delta21/m_e^2 = Omega_m^7 / [4 I^6 (Omega_m+I)]

These are consequences of the same prescriptions. They must not be counted as
new independent observations in addition to the original three quantities.

## Error interpretation and no joint likelihood

A standardized residual is computed as (prediction-observation) divided by the
quoted observational scale in quadrature with propagated input-scale uncertainty.
For DESI it is a Gaussian-style diagnostic using the quoted 68% posterior summary;
no full posterior or RMR-specific cosmological likelihood was evaluated. For the
muon, the very large standardized residual simply establishes that a 0.61% gap
is not explained by its experimental uncertainty. It is NOT an evaluated Gaussian
tail significance or a calibrated likelihood of the whole research program.

No joint chi-square, p-value, model evidence, or count of independent confirmations
is supplied. The primary observables originate in different experimental domains,
but they still have measurement/model conventions. Sensitivity rows explicitly
share data and are NOT added to the primary measurements.

## Cross-sector consistency without changing N

Inserting the measured muon ratio into the N-eliminated identities gives:

- effective N_mu = 137.8455218 +/- 0.0000030667;
- Omega_m implied by that identity = 0.2923003944437442;
- Delta21 implied by that identity = 7.211138288831791 x 10^-5 eV^2.

The last is approximately 2.407 quoted JUNO errors below its central value.
This is not an adopted continuous registry size or a fitted replacement model.
It demonstrates that using a common effective N to absorb the muon discrepancy
changes the neutrino prescription. The N=137 forecast remains unchanged.

The effective-N inverse from solar data uses a deterministic Decimal bisection.
The quoted intervals are monotone images of each observation's +/-one-error
endpoints, with m_e held at its listed central value. They are not posterior
confidence regions or a joint parameter fit. The m_e and mass-ratio input
covariance is not supplied; its possible effect on the cross-sector checks is
negligible beside the JUNO and DESI error scales. The root is checked by
independent floating-point evaluation as well.

## Absolute neutrino numbers: conditional, not measured here

Normal ordering, m1=0 and Delta32/Delta21=33 imply:

- m2 = 8.6778672854 meV;
- m3 = 50.6002266912 meV;
- sum = 59.2780939767 meV;
- Delta32 = 0.00248507756058 eV^2;
- Delta31 = 0.00256038294121 eV^2.

Delta31 and Delta32 are not interchanged. Solar agreement is not confirmation of
R=33, m1=0, or the absolute mass spectrum; those have not been measured in this
cell. These formulas yield an absolute scale in eV only after inputting m_e.

Supplementary comparison (not added as independent data): the 2026 DESI+CMB
analysis quotes sum < 0.0592 eV (95%) in LCDM with a free nonnegative mass sum
and three degenerate species. The conditional prediction is 0.000078094 eV above
that quoted upper bound. This near-boundary comparison is not a likelihood at
our distinct, hierarchical masses, and the paper itself warns that much of that
posterior is inconsistent with either oscillation-ordering floor. In w0waCDM,
its corresponding upper bound relaxes to 0.166 eV. Thus neither robust exclusion
nor independent confirmation follows from comparing only those marginal limits.
Do not change the cosmological model to whichever bound is favorable.

## Result

The selected uncorrected branch is NOT an exact jointly successful description
of the measured quantities, already because its muon ratio is wrong by 1.2682827.
The solar prediction and DESI-only matter fraction remain compatible with the
cited summary uncertainties. The data do not fix the cavity-registry coupling
mechanism; an independently derived correction would constitute additional work,
not a consequence of this statistical comparison.

All old dictionaries were preserved. Double replay and a deliberately altered
formula sentinel passed. The JSON encodes rational values with numerator,
denominator and a display decimal. No new geometry or native dynamics was invoked.
