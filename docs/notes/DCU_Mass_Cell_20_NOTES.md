# Cell 20 — frozen neutrino inputs to coherent vacuum survival

## Scope

This cell combines the **existing** screened-angle hypotheses from Cell19 with
the **existing** neutrino splitting prescription from Cell14. It borrows the
usual coherent, unitary, ultrarelativistic three-flavor propagation law. It does
not derive that law from native DCU service or ternary records. The physical
baseline is normal ordering. No source input is fitted, switched, or corrected.

The numerical source masses are in the convention c=1 (so the eV^2 values are
rest-energy-squared differences). The only empirical dimensional scale in the
primary neutrino parameters is the already-used electron rest energy,
0.51099895069 MeV. The exact SI h, c, and elementary-charge definitions supply
unit conversion, not an adjusted native length or clock calibration.

The result is a complete forward probability function, not detector counts.
Empirical-fit control curves are calculated from published fitted parameters;
they are not observations of probability at the example baselines/energies.
No goodness-of-fit claim or joint significance follows from curve distances.

## Unchanged source prescriptions

Screening paper, Eq.37 and Eq.39:

    s12 = sin^2(theta12) = 167/548
    s13 = sin^2(theta13) = 400/18769

Note: sin(theta13)=20/137 is squared ONCE. All digits are calculated from
fractions, not copied from rounded manuscript numbers.

Structural Taxonomy, solar mass-scale equation, plus its specified R=33 branch:

    d21 = (m_e c^2)^2/[4*137*136^6]
    d32 = 33*d21
    d31 = d32+d21 = 34*d21

The alternative angles 4/13,3/137 and alternative ratio137/4 are not substituted.
The source numerical predictions are:

    d21 = 7.530538062369455e-5 eV^2
    d32 = 2.4850775605819203e-3 eV^2
    d31 = 2.5603829412056147e-3 eV^2

A formula for splittings does not require inserting the absolute m1 into the
oscillation calculation. Setting the first propagation mass-squared to zero
is removal of a common phase; the code tests invariance under a common positive
mass-squared shift. This calculation therefore does not test m1=0, an absolute
mass sum, or the nature of Dirac versus Majorana masses.

## Borrowed probability law

The electron-row moduli squared are

    w1=(1-s12)*(1-s13), w2=s12*(1-s13), w3=s13, sum(w)=1.

Let y=L[km]/E_nu[MeV] and

    kappa = 1e-3/(4*hbar*c) = 1266.932679039099

when hbar*c is in eV m. Then

    P_ee(y)=1-4*w1*w2*sin^2(kappa*d21*y)
               -4*w1*w3*sin^2(kappa*d31*y)
               -4*w2*w3*sin^2(kappa*d32*y).

Equivalently, P_ee=|sum_i w_i exp[-2i*kappa*m_i^2*y]|^2.
It is the same for electron neutrinos and antineutrinos in the vacuum model.
Theta23 and the CP phase drop out of this disappearance channel. They have NOT
been set to convenient values to fit it, and a full predicted PMNS matrix has
not been obtained. The full-PMNS validation uses arbitrary settings only as
implementation checks. Unitary three-active-neutrino propagation and coherence
are borrowed assumptions, not new source results.

Related outputs:

    solar-term amplitude = 4*w1*w2 = 0.8117626794279422
    total fast amplitude = 4*s13*(1-s13) = 0.08343018914683618
    dm_ee^2 = (1-s12)*d31+s12*d32 = 0.002537434038716277 eV^2
    full phase average = sum_i w_i^2 = 0.5524035657126107

The full three-frequency equation is used throughout. The effective dm_ee is
reported as a derived near-reactor diagnostic, not substituted for both fast
terms at large L/E.

## A conditional exact signature: full vacuum revivals

Because d32/d21=33 and d31/d21=34, all phase differences in the squared-sine
formula are integer multiples of the solar phase. Thus

    T_y=pi/(kappa*d21)=32.92837707115688 km/MeV
    P_ee(y+T_y)=P_ee(y), P_ee(T_y)=1.

All source weights are nonzero, so this is the smallest positive full coherent
return period: the 21 phase must itself be an integer multiple of pi, and the
other frequencies are integer multiples. A common mass-squared offset changes
only the amplitude's overall phase. The source's solar/fast amplitudes determine
the pattern inside this period, but not the existence of a simultaneous revival.

This is an algebraic consequence of taking the frozen **effective** splitting
ratio as exactly33. It is not an independent physical derivation, a uniquely DCU
signature, or a measured full return. Any model with the same effective inputs
has the same waveform. If33 is instead only a structural-limit value with an
unspecified neutrino correction, the exact revival is a prediction of that limit,
not automatically of the corrected physical spectrum.

At this y the versioned fit controls give P_ee=0.9263160057610459 (IC24+SK NO)
and0.931882459557144 (IC19 noSK NO). Those are fitted-model central curves,
not measured probabilities or evidence of exclusion at any confidence level.

## Resolution sensitivity (not detector fitting)

The cell reports one optional analytical average over an ideal Gaussian y with
mean y0 and width sigma_y:

    <sin^2(k*y)>=[1-cos(2*k*y0)*exp(-2*k^2*sigma_y^2)]/2.

At the source's coherent return, sigma_y/y0=1% and3% give

    <P_ee>=0.9619269316513634, 0.9511379782266743.

The widths are declared illustrative sensitivity settings, not fitted detector
resolutions, native noise levels, physical decoherence parameters, or corrections
to the primary prediction. The Gaussian mathematical integral is over the real
line; for these positive means and widths the negative-y tail is negligible.
The formula has been independently integrated at80 digits.

The same nuisance description must be used on BOTH source and empirical curves
in any later data comparison. Resolution averaging is not a license to smooth
only the source until it appears closer to data.

## Empirical benchmarks and convention tracking

Use the same NuFIT6.0 release and NO variants as Cell19, based on data through
September2024, arXiv2410.05380v2 Table1. These are explicitly versioned benchmarks,
not claims of the newest global analysis as of execution.

- Primary IC24+SK NO: s12=.308, s13=.02215, d21=7.49e-5, d31=2.513e-3.
- Retained IC19 w/oSK NO: s12=.307, s13=.02195, d21=7.49e-5, d31=2.534e-3.

The NO atmospheric entry in that table is d31, not d32. The latter is obtained
by subtracting d21. Parameter correlations are not reconstructed from marginal
errors; no confidence bands, independent-error scans, or chi-square are used.

Two diagnostic hybrids keep the primary fit fixed and replace either its
angles or its splittings by the source values. They attribute the curve
difference, but never replace the source prediction or become extra fits.

Fixed L/E domains: [0,1] and[0,40] km/MeV, using 2001 and16001 points.
The root-mean-square discrepancy uses a trapezoidal uniform-L/E measure;
it is not flux-weighted, event-weighted, or a statistical residual. Maximum
values are grid maxima; doubled-grid checks are retained, not a proof of the
continuous supremum.

Source versus primary fit:

    domain [0,1]: RMS0.310060 percentage points, grid max0.544350 pp.
    domain [0,40]: RMS3.707516 pp, grid max8.923203 pp.

On[0,40], source angles with fit splittings differ by RMS0.344004 pp,
whereas source splittings with fit angles differ by RMS3.750270 pp.
This conditional ablation shows that accumulated phase offsets dominate that
chosen unsmeared comparison; it is not a detector-likelihood result.

DayaBay final nGd sample (arXiv2211.14988) is retained as an additional
PARAMETER-SUMMARY comparison:

    sin^2(2theta13)=0.0851 +/-0.0024
    d32=(2.466 +/-0.060)e-3 eV^2 (NO).

Here the reported gap is d32, NOT dm_ee. Source discrepancies are -1.962175%
and+0.773624%, respectively; marginal quoted-error units -0.696 and+0.318.
DayaBay also contributes to NuFIT; these checks are not independent successes.
Neither this final sample nor the source masses were historically blind.

## What is not calculated

No solar production/matter profile: a 'solar mixing angle' does not imply these
vacuum curves predict Sun-to-Earth survival. No Earth matter potential, baseline
distribution, energy-response matrix, reactor spectral flux, inverse-beta cross
section, backgrounds, acceptance, or efficiency is supplied. The illustrative
52.5km baseline is a single baseline, not a reconstruction of JUNO.
Neutrino energy is used, NOT detector prompt visible energy.

A subsequent experimental prediction must fold the same P_ee with specified
source and detector inputs, schematically:

    counts_b = normalization * sum_r integral dE
       flux_r(E)/L_r^2 * sigma(E) * efficiency(E)
       * response_b(E) * P_ee(L_r,E) + backgrounds_b.

Matter propagation must also be supplied where relevant. This cell stops
before that step; raw counts and fitted parameter curves are not interchanged.

## Sources

Source theory supplied in the conversation:
- screen(3).pdf, Eqs.(37),(39), with the source's T2 qualifications retained.
- Structural Taxonomy draft, solar splitting equation and R=33 normal-ordering branch.
- Existing executable Cells14,18,19 supply the exact frozen parameters.

External propagation/measurement sources checked for this iteration:
- PDG2025 Neutrino Masses, Mixing, and Oscillations, Eqs.(14.35)-(14.42),
  and(14.80)-(14.81). Rendered pages8,9,35 were inspected.
  https://pdg.lbl.gov/2025/reviews/rpp2025-rev-neutrino-mixing.pdf
- NuFIT6.0 arXiv2410.05380v2, Table1 and definition of dm3l.
  https://arxiv.org/html/2410.05380v2
- DayaBay Collaboration, Precision measurement of reactor antineutrino oscillation
  at kilometer-scale baselines, arXiv2211.14988, abstract/final3158-day sample.
  https://arxiv.org/abs/2211.14988
- CODATA2022 constants, electron scale and SI defining constants.
  https://physics.nist.gov/cuu/Constants/Table/allascii.txt

## Verification

The original Cells01,14,17,18,19 and source archive were replayed. Two executions
of the final cell agree in every retained non-callable result and printed output.
5005 internal probability/amplitude checks and90 independent full-PMNS checks
agree. Independent80-digit amplitudes (210 points) give maximum absolute error
8.27e-16. Both Gaussian sensitivity integrals agree at5.9e-17 or better.
Input changes are rejected; source dictionaries are preserved. Detailed output
and a development grid-refinement note are in DCU_Mass_Cell_20_TESTS.txt.
