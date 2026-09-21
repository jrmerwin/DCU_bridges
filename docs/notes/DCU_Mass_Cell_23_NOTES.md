# Cell 23 — production thresholds and finite-quantum bridge audit

## Decision and scope

Leave Cell 22 as the atmospheric-frequency checkpoint. It is a projected two-parameter Daya Bay compatibility test under the collaboration's solar/nuisance treatment, not a completed joint solar–atmospheric likelihood or evidence for literal integer exactness. This cell changes no neutrino parameter and performs no new oscillation fit.

The present cell opens two separate branches:

1. Frozen source mass prescriptions plus borrowed special relativity give production boundaries and near-threshold speeds.
2. Native registry adjacency plus the explicitly adopted complex record representation support finite operator checks. These distinguish nonunitary graph transfer from possible quantum generators and check finite interference.

No modified dispersion relation is derived. No new native service or construction history is simulated. Nothing in this cell identifies a recording unit with energy or a maintenance tick with a second.

## 1. Mass and channel inputs

Reused exactly from Cells 18/19, in electron-mass units:

- muon: 827/4;
- tau: 13909/4;
- proton: 1836 + 20/137;
- charged pion: 273 + 20/137 (only charge-bookkeeping negative controls here);
- Lambda: 2183 + 45/137.

Structural-only controls remain 205.5, 3493.5, 1836, 273, and 2183. The electron rest-energy unit remains 0.51099895069 MeV, with its previously declared experimental provenance. Particle–antiparticle mass equality is borrowed. Source mass dictionaries do not derive spin, charge, baryon number, strangeness, or conservation laws. The channel quantum numbers in the code are external conventional identifications.

The kaon mass is NOT part of this branch's claimed source predictions. The screening paper gives 966+45/137 but explicitly leaves the 966 base unexplained (Eq. 32 and ensuing discussion). For reactions requiring K+, the same measured kaon mass is supplied to both source and structural-control calculations:

    m_K+ c^2 = 493.677 +/- 0.015 MeV.

This is PDG 2025's `OUR FIT`, page 6 of its charged-kaon listing. The historical review text/`OUR AVERAGE` on the same listing has uncertainty .013 MeV; these are different summary conventions. The code uses .015 consistently and no mass uncertainty is fitted. The kaon is an additional species input, not a second freely fitted overall mass scale.

For reference calculations only, use NIST CODATA2022 mu/e=206.7682827 and p/e=1836.152673426, PDG2025 tau=1776.93 +/- .09 MeV, and Lambda=1115.683 +/- .006 MeV. The tau benchmark is explicitly PDG2025, not CODATA2022's older tau value; it does not change any earlier cell's comparison data.

## 2. Relativistic threshold equations

Masses below are rest energies in MeV (c=1). For an endothermic free-particle final state with total rest mass M and incident beam/target masses a,b:

    W_min = sum final masses = M,
    s = a^2+b^2+2 b E_beam (target at rest),
    E_beam,min = (M^2-a^2-b^2)/(2b),
    T_beam,min = [M^2-(a+b)^2]/(2b).

At threshold all final particles have the same laboratory velocity; they need not each be at rest in the laboratory. The code verifies their individual mass-shell equations and total four-momentum using exact rational arithmetic. These formulas describe free-particle kinematic boundaries in the narrow-width approximation, not resonances, Coulomb bound-state onsets, experimentally smeared thresholds, or a nonzero cross section.

For two final particles:

    p_*^2 = [s-(m_1+m_2)^2][s-(m_1-m_2)^2]/(4s).

The physical condition s >= (m_1+m_2)^2 is checked separately. Positivity of the Kallen polynomial below the pseudothreshold is not sufficient for physical production.

### Charge correction to the proposed pion example

K+ Lambda has total electric charge +1. With a proton target, the two-body pion channel therefore requires pi0. Neither pi+ p nor pi- p has the requisite charge. A charged-beam channel uses pi- p -> K0 Lambda instead. This cell evaluates gamma p -> K+ Lambda and pp -> p K+ Lambda, avoiding a new neutral-pion or neutral-kaon mass prescription. This is a choice of a charge-allowed alternative, not a silent claim that a charged pion has the same channel.

## 3. Executed threshold values

| Channel | Source COM boundary (MeV) | Source fixed-target beam kinetic boundary (MeV) |
|---|---:|---:|
| e+e- -> mu+mu- | 211.298066110315 | 43684.85317040625 |
| e+e- -> tau+tau- | 3553.742202573605 | 12357249.051901 |
| e+e- -> p pbar | 1876.537343707604 | 3445595.488107 |
| e+e- -> Lambda Lambdabar | 2231.357111453869 | 4871784.479824 |
| gamma p -> p mu+mu- | 1149.566738 | 235.090224 |
| gamma p -> K+ Lambda | 1609.355555726935 | 911.080831884180 |
| pp -> p K+ Lambda | 2547.624227580737 | 1582.167715757312 |

More digits and exact fractions are in RESULTS.json. COM pair-production boundaries and fixed-target values are different descriptions of the kinematics, not independent empirical tests. Comparing these values to thresholds computed from measured masses is a consistency/control comparison, NOT a new measured production-onset dataset.

For photoproduction the mass-control threshold is 911.081718760915 MeV, so the primary is approximately .887 keV lower. The inherited kaon-input uncertainty ALONE gives .025728594 MeV in the primary threshold through derivative (m_K+m_Lambda)/m_p. Thus this proximity is not a new sub-keV metrology claim. Other input errors and theoretical uncertainty have not been combined into a likelihood.

For pp associated production the control beam kinetic threshold is 1582.169629333591 MeV, approximately 1.914 keV above the primary. Close p and Lambda mass discrepancies partly cancel in these combinations. Such cancellation is not independent evidence for a mechanism.

## 4. Near-threshold amplification

For equal-mass pair creation at COM energy W:

    beta = sqrt(1-4m^2/W^2).

At W=3554 MeV, the source tau mass gives beta=.012044470133, whereas the PDG2025 central-mass calculation gives .008875974113: +35.6974% relatively. At W=3560 MeV this drops to +.953869%, and at 3600 MeV to +.126746%.

This is a SENSITIVITY example at chosen energies, not a measured-speed discrepancy. The source pair threshold is only .117797426 MeV below the central-mass boundary, or -.00331463%. The pair threshold uncertainty from the cited tau measurement alone is .18 MeV, larger than the central .14-MeV excess of 3554 MeV over 3553.86 MeV. An actual threshold scan needs beam-energy distribution, widths/radiative corrections, final-state interactions and a production amplitude. No cross section, yield, or confidence statement is inferred.

The general lesson is that relative accuracy is not invariant under transfer. Near a threshold, p scales as sqrt(excess energy), and a small absolute mass error can dominate the small excess.

## 5. Registry transfer versus a quantum generator

Construct from the ACTUAL reference adjacency A and degree matrix D:

    T = -A D^-1,
    S = D^-1/2 T D^1/2 = -D^-1/2 A D^-1/2.

T is not Hermitian in the ordinary Euclidean inner product on node coordinates. It is self-adjoint in the weighted inner product with metric D^-1, and S is its symmetric representation. This does not make T a unitary one-step propagator:

    ||T e_j||^2 = 1/d_j = 1/73 or 1/136.

Simply squaring its entries does not give a normalized transition probability per column. Squaring and renormalizing would be a different rule, not a proof of unitary evolution.

S has the same exact eigenvalue inventory checked in Cell14:

    -1 (x1), -62/73 (x1), 383/4964 (x1), 1/73 (x124), 1/136 (x10).

A Hermitian energy operator and the unitary law may be ADOPTED, e.g. H=epsilon*S. An equally well-defined but DIFFERENT energy assignment is H=epsilon*|S|, where the bars mean the spectral absolute value, not elementwise absolute values. Additive identity shifts do not affect transition probabilities. The graph spectrum alone does not determine which assignment represents physical energy. The code checks unitarity of the two mathematical propagators for an arbitrary dimensionless time, with no second/MeV map.

The distinction already matters for an isolated K_n cavity:

    signed-generator gap = 1 + 1/(n-1),
    magnitude-generator gap = 1 - 1/(n-1).

The resulting K5/K3 gap ratios are 5/6 and 3/2. The latter recovers the source's prescribed f5/f3 under a magnitude-based energy convention. The recovery is conditional on that energy convention, not a new derivation of it or of the additional factor137. A direct discrete update T^k and unitary exp(-i H t/hbar) are not interchangeable dynamics.

## 6. Finite characters and interference

Use the manuscript's finite cyclic address/character groups Z/(3^d), together with its DECLARED complex representation. In a q-dimensional computational basis:

    X|j> = |j+1 mod q>,
    Z|j> = exp(2*pi*i*j/q)|j>,
    ZX = exp(2*pi*i/q) XZ.

For q=3 and9 the cell verifies:

- unitary finite Fourier transforms;
- computational/Fourier overlaps squared exactly 1/q mathematically;
- q^2 trace-orthogonal operators X^a Z^b;
- a three-level interference circuit and its phase-erased control.

For a qutrit, apply Fourier preparation, diagonal phases (1,exp(2*pi*i*x),exp(4*pi*i*x)), then inverse Fourier readout. A finite triadic setting x=1/9 gives

    (.712386014201, .201689718788, .085924267010).

Complete dephasing between preparation and readout gives (1/3,1/3,1/3). At x=0 the coherent distribution is (1,0,0), up to numerical roundoff. These reproduce the form L(x)=(1+2cos(2*pi*x))^2/9 used in the attached instrument.

This is not a spatial wavefunction on R3, but complex state amplitudes, preparation, gates, and Born readout are still explicit representation assumptions. The q=9 Fourier transform is an ambient character-space audit, not a claim that a native F9 gate has been compiled into the construction/service grammar. No new Bell maximum or experimental quantum validation is claimed. The manuscript already gives stronger conditional Bell results for its specified tensor instrument.

Finite quantum mechanics does not have exact finite-dimensional canonical [Q,P]=i*hbar*I with hbar nonzero: taking the trace proves the obstruction. The periodic Weyl relation is an appropriate finite algebra. Its q^2 operator labels are not a measured physical phase-space area. Obtaining meters, momentum, and an action unit requires an additional valuation.

## 7. Dispersion status

A finite 137-root registry is not automatically a lattice of137 physical sites. The native process can grow indefinitely;137 is the selected ancestry-size registry. The current service law permits pairing across arbitrary existing graph distance. Neither a physical wavelength nor momentum nor local propagation scale has been supplied by that counting law.

A candidate correction such as

    E^2 = p^2 c^2 + m^2 c^4 + xi (pc)^(n+2)/E_*^n

would need a defined energy/momentum map, propagation operator, clock, scale E_*, exponent n, sign, species dependence, and any frame/direction convention. Writing xi=1/137 or5/137 does not determine these. No coefficient, length or time scale is introduced in Cell23. Discreteness does not by itself require Lorentz breaking; the cited causal-set theorem gives a distinct model where an intrinsic sprinkling cannot pick a preferred frame. This is a logical caution, not a proof of Lorentz invariance for DCU.

## 8. Rough bridge map

| Area | Current useful content | Remaining boundary |
|---|---|---|
| Mass -> kinematics | Recoil, endpoints, allowed invariant masses, production thresholds | Conservation/dispersion borrowed; yields and threshold shapes need interactions |
| Mass -> decay rates | Tau leptonic rates, pion e/mu with specified QED correction | External weak law; composite baryons need form factors and flavor couplings |
| Mass -> spectroscopy | Reduced-mass optical intervals | Differential bound-state QED, recoil, nuclear structure and source mass errors |
| Frozen neutrino parameters -> oscillations | Vacuum waveform, approximate detector fold, published atmospheric-surface point | Full joint solar/atmospheric likelihood with common priors; exact R unresolved |
| Finite records -> quantum algebra | Finite characters, interference, declared instrument's conditional Bell results | Physically selected states/generators/readout; no derivation from counts alone |
| Native geometry/time -> spacetime propagation | Ancestry/record metrics and operational service clocks | Metric/time/energy valuation and wave operator before limiting speed or LIV |

## 9. Sources and reproducibility

Author-provided sources:

- screen(3).pdf: muon/tau corrections Eqs.7,9; proton15; Lambda20; kaon32 with explicit open-base discussion.
- leptons.pdf: Eq.4 transfer operator, Eq.5 prescribed beat, Eqs.9–10 structural lepton bases.
- DCU revised manuscript 0.5: Sec.4.3 microscopic locality; Sec.5 finite address/clock characters and explicit representation choices; Sec.7.1 quantum instrument and conditional Bell scope.
- Actual Cells01,14,17,18,19 are replayed to generate inputs. Cell22 notes define the preserved frequency checkpoint; Cell22 is not an execution dependency.

External sources, read September19,2026; explicitly versioned inputs, not latest-experiment claims:

- PDG2025 Kinematics, Eqs.49.2–49.7 and49.16–49.17, plus cross-section discussion; rendered pp.1,3 inspected:
  https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf
- NIST CODATA2022:
  https://physics.nist.gov/cuu/Constants/Table/allascii.txt
- PDG2025 tau p.1, Lambda p.1, K+- p.6 (fit value), rendered pages inspected:
  https://pdg.lbl.gov/2025/listings/rpp2025-list-tau.pdf
  https://pdg.lbl.gov/2025/listings/rpp2025-list-lambda.pdf
  https://pdg.lbl.gov/2025/listings/rpp2025-list-K-plus-minus.pdf
- Farhi and Gutmann, Quantum Computation and Decision Trees, continuous-time quantum evolution on a graph is an explicitly specified Hamiltonian construction:
  https://arxiv.org/abs/quant-ph/9706062
- Gross, Hudson's Theorem for finite-dimensional quantum systems, finite phase-space/Weyl framework:
  https://arxiv.org/abs/quant-ph/0602001
- Bombelli, Henson and Sorkin, Discreteness without symmetry breaking: a theorem:
  https://arxiv.org/abs/gr-qc/0605006

Cell23 runs after Cells18/19 with NumPy. No network or extra files are required. Its output is dcu_mass_23. Exact thresholds, controls, scenarios, operator results and limits are retained. Run twice identically, with independent high-precision and matrix-exponential checks recorded in TESTS.txt. None of the numerical verification establishes physical precision or independently validates a mass formula.
