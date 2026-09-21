# Cell 18 — frozen mass prescriptions into recoil, spectroscopy, and pion branching

## Question and result

Expand Cell 17 into other observables without changing its structural or screened lepton masses. This cell adds the explicit charged-pion and proton prescriptions already in the supplied screening paper. It does not derive new species assignments or adjust masses using the comparison data.

Three branches are evaluated, all retained:

| Quantity | Screened prediction | Comparison | Screened fractional error | Structural-control error |
|---|---:|---:|---:|---:|
| Rest-pion muon recoil momentum | 29.804657018 MeV/c | 29.79200(11) MeV/c | +0.0424846% | +1.4643369% |
| H-minus-muonium 1S–2S interval difference | 10.533732607 THz | 10.532472187 THz | +0.0119670% | +0.6935465% |
| Pion e/mu inclusive rate ratio, compared with a leading-order prediction | 1.28262214387e-4 | 1.2327(23)e-4 | +4.0498210% | +2.2789958% |

These are approximation errors, not agreement within experimental uncertainty. There is no aggregate pass count or significance. The third row deliberately exposes a known limitation of omitting radiation in the imported law. All predictions and comparisons are retrospective; none is claimed to be historically blind.

## What remains fixed, and what is newly borrowed

Cell 17 source mass ratios remain exactly:

- Structural: r_mu=411/2 and r_tau=6987/2.
- Screened: r_mu=827/4 and r_tau=13909/4.

Here r=m/m_e. Add the supplied screening-paper formulas, not its rounded result columns:

- r_pi,0=273; r_pi=273+20/137=37421/137 (screen(3).pdf Eqs. 30–31).
- r_p,0=1836; r_p=1836+20/137=251552/137 (Eqs. 14–15).

The additional pion/proton assignments are pre-existing source hypotheses. Merely counting the registry does not force those physical assignments. No kaon integer 966 or Xi bases, whose source provenance is less independent, are used.

Imported physics:

1. Standard relativistic two-body kinematics; neutrino mass neglected.
2. Universal leading-order pseudoscalar leptonic decay, with the helicity-suppression factor. This is not Cell 17's three-body phase function applied incorrectly to a pion.
3. Nonrelativistic Coulomb gross-structure energies, equal unit charges, and reduced-mass scaling. This is not a new native DCU quantum/clock law.

There are two SEPARATE dimensional inputs for separate branches:

- Electron rest energy 0.51099895069(16) MeV for MeV/c and MeV recoil outputs. Its uncertainty is 0.00000000016 MeV.
- Hydrogen 1S–2S hyperfine-centroid frequency 2466061413187035(10) Hz for the optical frequency scale (Parthey et al., 2011).

The pion e/mu ratio and recoil beta=v/c require neither dimensional input. No measured tau/muon/pion/proton mass is supplied to the primary predictor. Cell 17's muon lifetime remains unchanged but is not used by this cell. No measured Rydberg constant, alpha, G_F, f_pi, or V_ud is added to the primary predictions. The hydrogen calibration absorbs an effective optical normalization, NOT the separate unknown differential bound-state corrections.

## Equations

For pi -> mu nu at rest with m_nu=0, in electron units:

    p/(m_e c) = (r_pi^2-r_mu^2)/(2 r_pi)
    E_mu/(m_e c^2) = (r_pi^2+r_mu^2)/(2 r_pi)
    K_mu/(m_e c^2) = (r_pi-r_mu)^2/(2 r_pi)
    beta_mu = (r_pi^2-r_mu^2)/(r_pi^2+r_mu^2)

The returned quantities satisfy E_mu^2-p^2=m_mu^2 and E_mu+p=m_pi exactly in these units. The primary predicts K_mu=4.123621852 MeV, beta_mu=0.271512582, and E_nu=29.804657018 MeV. These are dependent kinematic outputs, not separate experimental confirmations or new propagation-law derivations.

For the two pion leptonic modes, the leading-order ratio is

    R_pi = (1/r_mu^2) * [(r_pi^2-1)/(r_pi^2-r_mu^2)]^2.

The same pion decay constant, CKM element and weak normalization cancel in this ratio. No pion lifetime is needed. The measured ratio includes radiation; this pilot does not yet include the corresponding radiative corrections.

For 1S–2S in hydrogenic unit-charge atoms, set

    reduced_H/m_e = r_p/(1+r_p)
    reduced_Mu/m_e = r_mu/(1+r_mu)
    nu_Mu = nu_H * reduced_Mu/reduced_H.

The more mass-sensitive quantity is

    (nu_H-nu_Mu)/nu_H = (r_p-r_mu)/[r_p(1+r_mu)].

Both expressions are calculated exactly with Fractions. The predicted muonium transition is 2455.5276805803455 THz versus the frozen experimental 2455.528941000 THz, or -0.513299 ppm. The absolute optical-frequency relative error is not the same as the more discriminating H–Mu difference error, +0.0119670%. They are the same measurement in two presentations, NOT two independent successes.

The leading-order logarithmic frequency sensitivity to r_mu is 1/(1+r_mu), about 0.0048: an impressive-looking small fractional line-frequency error substantially compresses the input mass error. The logarithmic sensitivity of the H–Mu difference is -r_mu/(r_p-r_mu)-r_mu/(r_mu+1), around -1.12. The cell retains these sensitivities and the corresponding recoil sensitivities.

## Measured-mass bridge control

After freezing both recipe predictions, run precisely the same imported formulas using:

    r_mu = 206.7682827 (CODATA 2022)
    m_pi c^2 = 139.57039 MeV (PDG 2025)
    r_p = 1836.152673426 (CODATA 2022).

This is labeled `measured_mass_bridge_control`; it never supplies a value to the primary predictions. Its results are:

- R_pi=1.28334576312e-4, +4.1085230% relative to the inclusive measurement.
- p_mu=29.792140913 MeV/c, +0.0004730% relative to the specified recoil measurement.
- H–Mu difference=10.532692369 THz, +0.0020905% relative to the specified frequency difference.

This distinguishes failure of the imported leading-order approximation from disagreement introduced by the mass inputs. It does not establish that ALL of the primary error is an imported-theory effect. For example, the muonium primary is about 1260.42 MHz below the measured line, while the measured-mass control is about 220.18 MHz below it. The primary's additional mass-input error remains visible.

## Dependence, omitted physics, and epistemic limits

- PDG explicitly identifies the 29.79200 MeV/c recoil measurement as an input to precise charged-pion mass determinations. The test uses the direct recoil observation, not a synthetic momentum manufactured from measured masses, BUT it is not statistically independent of the mass data that informed the historical source paper. It is not an additional independent test of special relativity.
- The 1999 muonium optical measurement can itself be inverted to estimate the muon/electron mass ratio. It is a real optical comparison but not a historically blind or wholly independent sample relative to all muon-mass information. No joint test with mass measurements is claimed.
- The H and Mu frequencies are individually measured line intervals. The difference is computed from them, with marginal uncertainty sqrt(sigma_H^2+sigma_Mu^2); possible common metrology correlations are ignored and no precision likelihood is constructed.
- Pion radiation is not optional at the per-mille level. Cirigliano–Rosell (2007), Table 2, reports a leading relative electromagnetic contribution -3.929%, plus smaller terms. Those are known external field-theory corrections, not evidence that this residual is RMR thermal noise. They are NOT applied or fitted in this cell. Their numerical coefficients depend on the masses and conventions, so blindly importing a complete physical-mass correction factor as a universal constant would need a separate audit.
- The optical transfer omits differential relativistic, higher-recoil, QED and proton-size contributions. Measuring hydrogen fixes the normalization but cannot make all its corrections identical to muonium's. The simple formula must not be called precision spectroscopy or within measurement uncertainty.
- All three bridges import substantial established physics. They provide useful consistency/portability tests of approximate mass prescriptions, not a derivation of those external laws from the registry. Another mechanism producing the same mass inputs would give the same outputs.
- No native clock-to-second, coordinate-to-length, energy-to-work, or mode-to-mass identification is derived. No geometry, service law, initial recipe or prior result was changed.

## Sources and version convention

The sources below were checked 2026-09-19. Values are deliberately versioned historical benchmarks, not claims of newest available experiments.

### Source prescriptions

- Supplied `leptons.pdf`, Eqs. (9)–(10), retained via Cell 17.
- Supplied `screen(3).pdf`, Eqs. (7), (9), (14)–(15), (30)–(31). The paper's 205.518 and 3493.480 printed bases are not substituted.

### External calibration and comparison sources

- NIST, CODATA 2022 constants table:
  https://physics.nist.gov/cuu/Constants/Table/allascii.txt
- Parthey et al., Improved Measurement of the Hydrogen 1S–2S Transition Frequency (2011), abstract and paper:
  https://arxiv.org/abs/1107.3101
- Meyer et al., Measurement of the 1s–2s energy interval in muonium (1999), abstract: 2 455 528 941.0(9.8) MHz.
  https://arxiv.org/abs/hep-ex/9907013
- PDG 2025 charged-pion listing, rendered pp. 1, 2, 5 inspected:
  https://pdg.lbl.gov/2025/listings/rpp2025-list-pi-plus-minus.pdf
  pp. 1–2 give direct recoil 29.79200 +/- 0.00011 MeV/c and discuss its connection to mass measurements.
  p. 5 gives inclusive R_e/mu=(1.2327 +/- 0.0023)e-4. The p. 4 value 1.230e-4 is a DIFFERENT quantity, B(pi -> e nu); it is not substituted for the channel ratio.

### Imported laws and limitation diagnostics

- PDG 2025 Kinematics review, Eqs. (49.16)–(49.17); rendered page 3 inspected:
  https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf
- Cirigliano and Rosell, The Standard Model prediction for R_e/mu^(pi,K), Eq. (2):
  https://arxiv.org/abs/0707.3439
- Cirigliano and Rosell, pi/K -> e nu branching ratios to O(e^2 p^4) in Chiral Perturbation Theory, Table 2 and discussion:
  https://arxiv.org/html/0707.4464
- Cortinovis et al., Update of Muonium 1S–2S transition frequency (2023), Eqs. (1),(3),(4), and QED discussion:
  https://arxiv.org/html/2301.12883

## Executed validation

1. Replayed original mounted Cells 1, 14 and 17 using the original reproducibility archive.
2. Executed Cell 18 twice with exactly identical printed output and all non-callable result fields.
3. Exact Fraction checks verify source arithmetic, on-shell/energy identities, the two equivalent pion-ratio formulas, and the isotope-difference formula.
4. Doubling the electron energy calibration doubles dimensionful recoil and leaves beta/branching unchanged; doubling the hydrogen frequency doubles optical outputs.
5. An independent 80-digit Decimal calculation checks 12 outputs. Maximum relative discrepancy is 5.6e-79, consistent with independent decimal rounding of exact rational outputs. This tests arithmetic, not physical precision.
6. Six logarithmic derivatives checked independently by finite differences.
7. Altering Cell 17's screened muon input by 0.001 triggers explicit rejection.
8. The original registry and Cell 17 dictionary are unchanged after the experiment.

See `DCU_Mass_Cell_18_TESTS.txt` for the executed validation output.
