# Cell 17 — locked mass-to-decay-rate transfer

## Outcome

Using the frozen structural mass formulas plus the previously stated screening
terms, the two calculated tau leptonic partial decay rates differ from the
selected HFLAV comparison summaries by -0.31178% (electron channel) and
-0.36592% (muon channel). The structural-only control differs by +5.18104% and
+5.18529%. No term was changed using either rate comparison.

This is a **conditional, retrospective new-observable transfer check**. It is
not a blinded prospective discovery, a native DCU decay calculation, or a
derivation of the universal weak interaction. It supports the practical use
of the selected approximate mass prescription with the explicitly borrowed
weak-decay law. It does not determine the microscopic origin of the masses or
the cause of their residual discrepancies.

## Frozen source prescriptions

1. `leptons.pdf`, Eq. (9): r_mu^0 = 3*137/2 = 205.5.
2. `leptons.pdf`, Eq. (10): r_tau^0 = 17*r_mu^0 = 3493.5.
3. `screen(3).pdf`, Eq. (7): Delta r_mu = +5/4.
4. `screen(3).pdf`, Eq. (9): Delta r_tau = -13*(5/4).

All r values are masses in electron-mass units. The primary screened recipe is
r_mu=827/4=206.75 and r_tau=13909/4=3477.25. The structural recipe is retained as
a control, not selected after a new fit.

This combines the March 14 structural formulas with the March 16 correction
terms. It does NOT adopt the screening paper's unexplained alternative bases
205.518 and 3493.480. The paper's printed final mass values are not substituted
for the resulting arithmetic. The empirical integer mapping is a source
hypothesis, not a result of the native constructor.

Both recipes are anchored to the **same one measured muon lifetime**, using
exactly the same calibration procedure. Because they assign different muon
masses, their inferred effective weak-rate coefficients differ; this is
reported, not concealed. Neither coefficient is calibrated to tau data.

## Borrowed bridge and prediction

Adopt universal leading-order V-A charged-current three-body decay, neglecting
neutrino masses:

    rate(L -> l nu nu) = K * r_L^5 * F((r_l/r_L)^2)
    F(x) = 1 - 8*x + 8*x^3 - x^4 - 12*x^2*log(x)

The units of K here are inverse seconds because the r variables are
**dimensionless** masses normalized to m_e. All physical normalization,
including the necessary hbar and weak coupling in other conventions, is
absorbed into K by the clock calibration. No physical electron-mass value is
needed in the resulting ratios.

Set

    K = 1 / [tau_mu * r_mu^5 * F(1/r_mu^2)].

Then

    rate(tau -> l nu nu) = (r_tau/r_mu)^5
                          * F((r_l/r_tau)^2) / F(1/r_mu^2) / tau_mu.

The two daughters are r_l=1 and r_l=r_mu. This is a *partial* rate. Its inverse
is a partial lifetime, not the total tau lifetime. Tau hadronic widths are not
calculated; neither the total tau lifetime nor the absolute branching
fractions are independently predicted by this cell.

The experimental comparison rate is B_l/tau_tau. Tau lifetime and branching
fractions only enter this comparison. No measured electron, muon, or tau mass
is passed to the prediction routine. The Fermi constant is not introduced as
an additional calibration; its usual extraction already uses the muon
lifetime.

The electron/muon channel ratio is also retained. It is algebraically a ratio
of the two predicted partial rates, not a third independent test. Its small
mass sensitivity makes it much less discriminating between the two mass
recipes than either partial rate. The structural control actually lies closer
to the chosen HFLAV central channel-ratio value; both are well within its
quoted uncertainty. That result is retained rather than omitted.

## External sources and measurements

Retrieved September 19, 2026, using the versions listed below. These are frozen
benchmark summaries, not live-updating inputs. Current HFLAV's tau page points
to its end-2023 report (web release 15 January 2025).

### Clock calibration

PDG 2025 muon listing, page 2, "OUR AVERAGE":

    tau_mu = (2.1969811 +/- 0.0000022) microseconds
           = 2.1969811e-6 +/- 2.2e-12 seconds.

https://pdg.lbl.gov/2025/listings/rpp2025-list-muon.pdf

The page was inspected as a rendered PDF. This input is an empirical time
calibration, not a derived native tick-to-second map.

### Tau comparison data

HFLAV end-2023 branching-fraction report, Table 3, ordinary average rows
**without the additional unitarity constraint**:

    B_e  = 0.1784  +/- 0.0004
    B_mu = 0.17360 +/- 0.00037
    B_mu/B_e = 0.9730 +/- 0.0022.

https://hflav-eos.web.cern.ch/hflav-eos/tau/end-2023/br-fit.html

The report's separate universality-improved B_e is deliberately NOT used.
Such a value incorporates the weak-universality relation being tested here.
The term "unconstrained" means no extra total-unitarity constraint; it does
not mean there are no other fit constraints or experimental correlations.
Section 4.3 has one apparent wording slip ("universality-constrained");
Sections 4.1/4.2 and Table 3 identify the actual additional constraint as
unitarity. The report calculates tests of lepton universality separately in
Section 5.

The branching table PDF pages 5 and 6 were inspected. Values are used to their
published rounding precision. The separately rounded fitted channel ratio is
not forced to equal the quotient of rounded branching fractions exactly.

HFLAV end-2023 lifetime average, Table 2:

    tau_tau = (290.29 +/- 0.53) femtoseconds.

https://hflav-eos.web.cern.ch/hflav-eos/tau/end-2023/tau-lifetime-avg.html

Current report navigation:
https://hflav.web.cern.ch/content/tau

### Physics convention

The explicit decay law, phase factor, and radiation conventions are in:

https://hflav-eos.web.cern.ch/hflav-eos/tau/end-2023/lepton-univ.html

For an additional derivation/discussion:
A. Pich, *Precision Tau Physics*, arXiv:1310.7922, Section 2:
https://arxiv.org/html/1310.7922v2

The measurements are inclusive of the radiation convention specified by the
experimental analyses. This pilot omits relative radiative and finite-W
corrections, setting their inter-channel/parent ratios to one. Calibration to
the measured muon lifetime absorbs its inclusive effective normalization, but
does not automatically include the tau-to-muon difference. No theoretical
error bar or radiative parameter has been fitted to cover a discrepancy.
Consequently this is not an exact electroweak precision prediction.

## Error interpretation and dependence

For a constant phase factor, a fractional change in r_tau/r_mu changes the
predicted rate by about five times as much: d log(rate)=5 d log(r_tau/r_mu).
Phase-space derivatives provide the explicitly calculable additional term.
Thus a sub-percent mass model need not be equally accurate, percentage for
percentage, for every nonlinear transformed observable. This transfer tests
the propagated prediction, rather than assuming an invariant error level.

Marginal errors on B_l/tau_tau use first-order propagation with zero
lifetime/branching covariance. Both channels share tau_tau and their branching
fraction estimates are correlated. No joint chi-square, combined p-value, or
independence-based significance is reported. The 1% flag is a descriptive
approximation criterion, not a statistical confidence interval, a theoretical
uncertainty, or a historical prospective preregistration.

Source papers contain the charged-lepton masses and a rough lifetime table
(screening paper Section XII). These inputs were not historically blind to
the authors or this conversation. The present *partial-rate calculation* is
new to the notebook and uses no partial-rate target for calibration. Its
success adds an observable consistency check under imported dynamics; it
cannot distinguish mechanisms that already yield the same masses and obey
the same imported weak law.

## Executed validation

- Replayed the actual original Cells 1 and 14 from their mounted files and
  source reproducibility ZIP.
- Executed Cell 17 twice; outputs and all saved non-callable fields agree.
- Recomputed all four rate predictions independently with 70-digit Decimal
  logarithms and rational-to-Decimal mass inputs. Maximum relative difference
  from the floating computation: 2.05e-16.
- Checked equivalence of direct power-ratio and calibrated-coefficient forms.
- Doubling the supplied clock halves all rates, doubles partial lifetimes,
  and leaves the predicted channel ratio unchanged.
- Confirmed phase-space endpoints F(0)=1 and F(1)=0.
- Confirmed predictions are frozen before comparison annotations and contain
  no tau-lifetime or branching-fraction input.
- Changing the old Cell 14 muon prescription causes explicit rejection,
  rather than silent acceptance of another baseline.
- Earlier registry and Cell 14 dictionaries remain unchanged.

No new native history, fitted correction, dynamical update, or geometry
reinterpretation was introduced.
