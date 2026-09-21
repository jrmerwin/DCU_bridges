# Cell 19 — angle-ladder transfer, pion-specific radiation, and hyperon readiness

## Scope and result

One notebook cell runs after Cell 18 using the original registry and its frozen mass dictionaries. It adds no native DCU dynamics, mass fit, adjustable screening coefficient, or thermal-noise allowance. The source's four angle formulas are distinct physical hypotheses, not consequences obtained solely by enumerating the registry. All comparisons are retrospective and use named historical releases. This is not a claim to have found the latest available global oscillation analysis.

The leading pointlike pion QED calculation brings the screened ratio from +4.049821% to -0.037993% relative to the unchanged PDG inclusive benchmark. The measured-mass control changes from +4.108523% to +0.018330%. The unscreened control moves from +2.278996% to -1.734556%. Known omitted higher-order terms exceed the remaining primary residual, so the latter is not an established precision prediction.

The angle ladder has mixed results. The effective weak-angle prediction is 0.11743% high, about 2.091 quoted observational errors. Solar agrees approximately but is 1.057% below the selected NuFIT central value and 1.441% below the additional JUNO central value. The reactor sine-squared value is 3.7845% below the selected NuFIT central value. Atmospheric mixing is explicitly in the upper octant, whereas the selected normal-ordering IC24+SK central value is in the lower octant. Its 18.03% central-value discrepancy must not be converted to a Gaussian significance using only the local one-sigma error: the quoted marginal three-sigma interval includes the upper-octant prediction. All four NuFIT variants/orderings are retained without choosing the best match.

## Source hypotheses

`screen(3).pdf`, Eqs. (36)-(39):

- delta0 = 5/137, equivalently (40/137)/8.
- weak = sin^2(theta_eff^l) = 1/4 - delta0/2 = 127/548.
- solar = sin^2(theta12) = 1/4 + 3 delta0/2 = 167/548.
- atmospheric = sin^2(theta23) = 1/2 + 3 delta0/2 = 76/137.
- sin(theta13) = 4 delta0; therefore reactor = 400/18769.

The last expression is squared once, not treated as a sine-squared input before squaring again. The older taxonomy's 4/13, 7/13 and 3/137 formulas are not substituted. The source labels the angle prescriptions partially derived (T2), with base-state and radiative mechanisms still open. No CP phase is assigned, so this does not specify an entire complex PMNS matrix.

Exact consequences checked by Fractions:

    atmospheric - solar = 1/4
    solar + 3 weak = 1
    sin(theta13) + 8 weak = 2

They are algebraic consistency identities, not additional independent observations. No statistical test of these joint identities is performed without the appropriate likelihood/covariance.

## Angle comparisons

Primary neutrino benchmark: NuFIT 6.0, data through September 2024, arXiv:2410.05380v2 Table 1, IC24 with SK atmospheric data, normal ordering. Both published atmospheric-data variants and both ordering-conditional fits are saved. Inverted-ordering comparisons are sensitivity checks; they do not switch the source hypothesis to inverted ordering.

https://arxiv.org/html/2410.05380v2

The variants differ in their IceCube and Super-Kamiokande treatment, not just an isolated SK switch. Their datasets overlap. Table intervals are one-dimensional, profiled ranges. Satisfying each marginal three-sigma range is not proof of joint compatibility, and local one-sigma errors do not encode a two-octant likelihood. No summed chi-square or pass probability is calculated.

Additional solar check: JUNO first 59.1-day result, normal ordering, sin^2(theta12)=0.3092 +/- 0.0087. This is not substituted into the NuFIT analysis, and its dependencies on external oscillation/spectral information are not ignored by calling it fully independent.

https://arxiv.org/html/2511.14593v1

Effective weak angle: PDG 2025 Electroweak Model review, Eq. (10.65), collider average 0.23148 +/- 0.00013. This is not an on-shell angle or an MSbar coupling angle, and it is not merely the LEP+SLC average. The same review's Eq. (10.63) gives a different LEP+SLC combination.

https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf

The source's expression is compared as printed to the effective observable. Its description elsewhere as a tree-level value leaves an unresolved scheme/scale matching question. Electroweak effective-angle corrections use coupling form factors and a specified renormalization convention (review Eq. (10.51)). The pion decay correction cannot be reused for this quantity. No such reuse occurs.

## Pion-specific pointlike QED calculation

Use the identical tree-level mass function as Cell 18:

    R0 = (1/r_mu^2) [(r_pi^2-1)/(r_pi^2-r_mu^2)]^2.

Borrow the fully inclusive pointlike O(e^2 p^2) electromagnetic contribution from Cirigliano and Rosell, arXiv:0707.4464, Eqs. (95)-(97):

    Delta_pt = alpha/pi [F(1/r_pi^2) - F((r_mu/r_pi)^2)]
    R_pt = R0 (1 + Delta_pt)

    F(z) = 3/2 log(z) + (13-19z)/(8(1-z))
           - (8-5z) z log(z)/(4(1-z)^2)
           - [2+(1+z)log(z)/(1-z)]log(1-z)
           - 2(1+z)Li2(1-z)/(1-z).

https://arxiv.org/html/0707.4464

This includes the pointlike virtual and real-photon effects at the stated order. It is a Standard Model/chiral expansion ingredient, not a derived RMR screening law. It depends on the process and mass ratios. The published -3.929% is not pasted in as a universal constant: the function is evaluated separately for each frozen source recipe and for the measured-mass control.

One new measured input is supplied: the low-energy alpha, from CODATA 2022 inverse alpha=137.035999177 with standard uncertainty 0.000000021. This uncertainty is recorded, not used to create a complete theoretical uncertainty. The forward alpha is evaluated as the reciprocal of the printed inverse-alpha central value. No new dimensionful calibration is needed.

https://physics.nist.gov/cuu/Constants/Table/allascii.txt

The original pion ratio comparison, (1.2327 +/- 0.0023)*10^-4, is reused from Cell 18. It is the inclusive e/mu ratio on page 5 of the PDG 2025 charged-pion listing, not the different electronic branching fraction near it.

https://pdg.lbl.gov/2025/listings/rpp2025-list-pi-plus-minus.pdf

The source's Table 2 also gives higher contributions around +0.053%, +0.073%, and +0.055%, with structure/matching qualifications. They are omitted consistently, not selectively fitted. They are each comparable to or larger than the remaining -0.038% residual. Nearer agreement of this truncated result is not proof that those contributions vanish, or that all residuals are RMR noise. The primary and measured-mass calculations remain distinguishable: their mass-input difference is not erased by correcting the shared leading approximation.

## Hyperon readiness and arithmetic check

The source's Lambda expression (Eq. (20)) is r_Lambda=2183+45/137. With the unchanged electron energy input from Cell 18 it gives

    m_Lambda c^2 = 1115.6785557269346 MeV.

The PDG 2025 comparison is 1115.683 +/- 0.006 MeV. The signed difference is -4.4442730654 keV, within the quoted six-keV error but not the source's displayed sub-keV residual. This is a re-evaluation of the explicit formula and the fixed unit, not a changed prediction.

https://pdg.lbl.gov/2025/listings/rpp2025-list-lambda.pdf

Relativistic three-body endpoint kinematics with a massless neutrino gives

    E_e,max/(m_e c^2) = (r_Lambda^2 + 1 - r_p^2)/(2 r_Lambda)
    T_e,max = 162.793565064 MeV.

This value is stored as a forward kinematic output only: there is no endpoint measurement comparison in Cell 19. It is not a decay rate or a new experimental success.

Absolute Lambda semileptonic rates require a baryon weak-current matrix element with vector and axial form factors and a flavor-changing coupling. These do not follow from supplying the masses or from importing the muon lifetime alone. Cabibbo, Swallow, Winston, arXiv:hep-ph/0307298 Sec. 2.1, Eqs. (5)-(7), describes the additional inputs. No lepton-only m^5 law is misapplied to the baryon.

https://arxiv.org/html/hep-ph/0307298

The screening paper explicitly says that its Xi integer bases 2573 and 2586 were obtained by solving from experimental masses (page 9, Eqs. (24)-(25)). Those could be used in an empirically anchored transfer with clear labeling, but not called unseen parameter-free full mass predictions. Their charge-dependent screening coefficients are separate source hypotheses.

## Validation and reproducibility

The actual original mounted Cells 01, 14, 17, and 18 were replayed with the original reproducibility archive. Cell 19 was executed twice with identical output and all non-callable fields identical. Prior input dictionaries were preserved. The three radiative predictions were recomputed independently with mpmath at 80 digits; the maximum relative disagreement with the user-cell outputs was 5.071e-17. The user cell does not require mpmath.

The real dilogarithm implementation uses its convergent power series on [0,1/2] with an explicit truncation-tail bound and reflection elsewhere; it is tested against exact special-value/reflection identities. QED alpha=0 returns the previous tree output, and the correction is linear in alpha at this order. A deliberate +0.001 change to the frozen screened muon input triggers an explicit failure instead of refitting. These checks establish implementation accuracy, not physical validity.
