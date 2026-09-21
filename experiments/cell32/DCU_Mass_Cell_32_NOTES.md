# Cell 32 — a smaller vertical bridge: native timing to quantum correlations

## Decision and scope

Pause the persistent-observer identity/continuation branch at Cell31's explicit
counterexample. It is not a prerequisite for every vertical test. Return to the
already declared finite quantum instrument and ask how its state and measurement
statistics depend on actual native timing histories.

This is a retrospective analysis of ALL saved Cell29 histories. No new principal
native history, physical datum, unit calibration, mass fit, private service pool,
spatial metric, observer continuation, or neutron-search edit is introduced.

The primary result is a quantum channel **conditional on the existing attached
complex representation and Cell25/29 service-to-phase/phase-retention rule**.
Classical sampling does not manufacture entanglement; the entangled initial state
is supplied by the declared instrument. The other native constructions do not
otherwise disturb this spectator state under the inherited assumption.

## 1. What the papers actually support

The directly supplied September 2026 revised manuscript 0.5, *Self-Limited Growth,
Intrinsic Time, and Record Geometry in the Distinction Combinatorial Universe*,
provides:

- Sections 2.3–2.5: the precise native service/construction law, ancestry-factor
  recording work, and operational lapse. Equations (15)–(16) use **Phi squared**
  in the serviced tick rate; neither Phi nor that counter is automatically proper
  time.
- Section 6.3: passive regional labels have identical expected service per token.
  Independent local pools would change the dynamics.
- Section 7.1: the conditional shared-source state, phase family, and two global
  Bell values within the fixed instrument catalogue. These are mathematical
  results under stated representation premises, not an independent empirical
  derivation of all quantum mechanics.
- Appendix B.3: the signed CGLMP-3 functional and its continuous maximizing
  settings (0,1/6,-1/12,1/12). The pure-source maximum is
  4/3+8sqrt(3)/9 = 2.8729340511723355... .
- Section 7.1: the allowed heralded operation P01 F3 P01 has success 4/9 and
  Schmidt squares (3/4,1/4,0).
- Sections 7.2–7.3 and Appendix C: historical physical comparisons are partial;
  regional charge arithmetic is not a gravitational field or horizon. The current
  manuscript does not supply an observed local relativistic redshift or a
  Friedmann-dynamics derivation.

The related August 2026 *Combinatorial Descent of an Internally Unified
Quantum-Gravitational Equation System* is explicit about the earlier language:
its quantum-side and gravitational-regime sectors are model-internal. The latter
refers to persistent finite-capacity load. Its Sections 1.1 and 2.1 distinguish
this from empirical equivalence to physical QM/GR. It also separates lineage/
record dynamics from rendered service, coupled in one composite transition.
These are different source scopes, not interchangeable models.

The RMR screening/mass papers are not used as dynamical inputs here. Their
physical interpretations are neither erased nor imported into the native sampler.

## 2. Inputs and fixed protocol

Use `reference/DCU_Mass_Cell_29_RESULTS.json.gz` without alteration:

- Six arms: existing-pair, recorded-parent new pair, first-use new pair, in each
  of the two existing relief environments H=0 and H=6.
- 512 histories per arm, 3072 total; all seven checkpoints 0,1,8,32,128,512,1024.
- Each history retains cumulative services N3,N4 to the original sibling factors.
- The primary phase increment is 1/27; 1/81 is the existing finer-setting control.
- Construction, queues and preparation costs are those of Cell29. Changing arms
  involves the actual distinct pairs described there, not an artificial toggle
  of one charge bit on otherwise identical events.

Use the fixed legal settings

    (x0,x1,y0,y1) = (0, 4/27, -2/27, 2/27).

These approximate the paper's continuous optimum on the denominator-27 grid,
with the halfway choice for x1 rounded toward zero. The same quartet is used for
both phase increments, every history and every arm. No grid search, fitted
phase offset or outcome-dependent selection is used.

This is not preregistered: a short diagnostic calculation preceded the written
protocol. The filtered-state transfer was added after the direct-source result,
using the source's fixed operator, not an outcome-selected state. All direct-source
scientific values were preserved when the transfer was added.

## 3. The local random-unitary channel

Let D(u)=diag(1,exp(2pi i u),exp(4pi i u)). For phase modulus Q=3^d, the service
history supplies

    U_(a,b) = D(a/Q) tensor D(b/Q),
    a=N3, b=N4 modulo Q.

The joint histogram p_(a,b) determines the full channel

    E(rho) = sum_(a,b) p_(a,b) U_(a,b) rho U_(a,b)^dagger.

This is a completely positive trace-preserving map: its Kraus operators are
sqrt(p_(a,b))*U_(a,b), and their adjoint products sum to the identity. Correlations
between the two service streams remain present; independent Bernoulli noise is
not substituted. Every term is a product of local phase operations, so the
service noise does not itself create entanglement from a separable preparation.

This is the ensemble map at a specified checkpoint, not a claim of a time-
homogeneous Markov master equation for the reduced apparatus. The growing native
queues correlate service increments. Products of unconditional one-step moments
are not used to reconstruct the multi-step result.

### Two moments suffice for the diagonal shared-source state

For Omega=(|00>+|11>+|22>)/sqrt(3), K=N3+N4,

    |psi_K> = sum_r exp(2pi i r K/Q)|rr>/sqrt(3),
    C_h = E exp(2pi i h K/Q), h=1,2,
    rho = (1/3) sum_(r,s) C_(r-s) |rr><ss|,

where C_0=1 and C_(-h)=conjugate(C_h). Thus the source's complete density matrix,
not only one fringe height, follows from these two moments.

IMPORTANT: K alone does not define the channel on arbitrary two-qutrit states.
The full implementation retains the JOINT (N3,N4) histogram. The sum-only
compression is justified specifically by the |rr> support of Omega.

### Exact one-step law, not an independence assumption

At one native state R=F+n, q=3, p=3/R and p2=6/[R(R-1)], the sum increment J is
0,1,2 with probabilities

    (1-2p+p2, 2(p-p2), p2).

For z=exp(2pi i h/Q),

    E[z^J | pre-service state] = 1+2p(z-1)+p2(z-1)^2.

This verifies the immediate channel from the native sampler. The later channel
is instead computed from the full saved trajectories.

## 4. Fixed Bell readout and entanglement

For theta=x-y+(a-b)/3, the output table is

    p(a,b|x,y) = 1/9
              + (4/27) Re[C1 exp(2pi i theta)]
              + (2/27) Re[C2 exp(4pi i theta)].

It has normalized probabilities and uniform 1/3 local marginals. All four
measurement settings use the same native-history distribution. Settings are not
fed into the counting process, and no history is removed based on settings or
outcomes.

The signed functional is exactly Appendix B.3's expression. Enumerating all 81
local deterministic strategies gives scores -4 (3 strategies), -1 (48), and
2 (30). The relevant local upper bound is 2, NOT a symmetric absolute-value
bound of 2. The fixed finite setting gives pure-source I3=2.8399975540209894.
The continuous source optimum is used only as a formula regression.

### Partial-transpose diagnostic

The partial transpose of this diagonal-source density has eigenvalues

    1/3 (three times),
    +/-|C1|/3 (each twice),
    +/-|C2|/3 (each once).

Its negativity (sum of magnitudes of negative eigenvalues) is

    N(rho) = (2|C1|+|C2|)/3.

A separable state's partial transpose is positive: transposing the second
positive factor of each positive product term preserves positivity. Therefore
any strictly negative eigenvalue here certifies entanglement. In this particular
state family, C1=C2=0 also gives the explicit separable mixture
(1/3)sum_r |rr><rr|. These statements do not rely on an optimizer.

Purity is (3+4|C1|^2+2|C2|^2)/9. Each conditioned pure history has negativity 1
and purity 1; discarding the phase history generally reduces both.

The finite empirical mixture's negativity is an estimator for the population
mixture. Magnitudes are nonlinear and biased near zero; the jackknife error in
the results is a Monte Carlo diagnostic, not an experimental confidence bound.
Bell-score errors are standard errors of per-history EXPECTED scores, not
measurement shot-noise errors. No quantum detection experiment was conducted.

## 5. Results

At t=128 and increment 1/27:

| Environment | Arm | I3 +/- MC SE | Negativity |
|---|---|---|---|
| H=0 | existing pair | 2.319894878 +/- .029906993 | .923080857 |
| H=0 | recorded new pair | 2.580946071 +/- .019897365 | .947610103 |
| H=0 | first-use new pair | 2.715763219 +/- .012057190 | .967763653 |
| H=6 | existing pair | 2.030202748 +/- .042187147 | .880292393 |
| H=6 | recorded new pair | 2.211296752 +/- .036138862 | .901334829 |
| H=6 | first-use new pair | 2.421790288 +/- .026735971 | .927178518 |

The H=6 existing-pair estimate is close to the bound relative to Monte Carlo
error; it is not a high-significance population-level Bell violation. At this
fixed early checkpoint the new-record arms remain nearer the initial fringe and
retain more entanglement. This does NOT establish a general theorem that more
recording preserves entanglement or that gravity protects it.

At t=1024 all six fixed I3 values lie between -.8503 and -.7439. Yet their
negativities lie between .4394 and .5901. Failure of one fixed witness does not
show separability, nor does it prove that no other measurement violates a Bell
inequality. Measurement settings were deliberately not reoptimized.

All checkpoints, both increment depths, both 256-seed blocks, full probabilities,
and densities remain in RESULTS. No favorable time window replaces the rest.

## 6. Recoverable phase uncertainty versus intrinsic destruction

For each history, local inverse operations D(-N3/Q) and D(-N4/Q) restore the
initial source exactly. The counts are sufficient modulo Q; full unwrapped
history is not needed for this algebraic inversion. Every history is retained.
This is not a postselected Bell sample.

The operations belong to the declared triadic phase family. A physical observer
able to retain, read and act on the necessary local counts at a specified time
has not been compiled into native work; no free implementation is claimed.

Complete computational phase erasure gives I3=0 and negativity=0. Merely failing
to keep the service history yields the calculated mixed state, but does not by
itself demonstrate irreversible collapse or a fundamental decoherence process.

The pure mean-phase surrogate has negativity 1 and is not the same state as the
native-history average. It is retained only as a negative approximation control.
Clock slowing changes phase location, while history dispersion reduces ensemble
coherence. They must be distinguished before making a time-dilation claim.

## 7. Transfer to the paper's allowed filtered source

Apply P01 F3 P01 on Alice to Omega BEFORE clock evolution, with the inherited
positive-exponent F3 convention. Its unnormalized success probability is 4/9,
independent of future service and measurement settings. Conditional on success,
its nonzero coefficients are

    f_(r,s) = exp(2pi i rs/3)/2, r,s in {0,1}.

Its Schmidt squares are (3/4,1/4,0), and its initial negativity is sqrt(3)/4.
Use precisely the same joint-count channel and fixed measurement quartet. Do
not rotate into a convenient Schmidt basis without also transporting the phase
operations and measurements.

For H=0 existing-pair at t=1024, increment 1/27, its negativity is .260737815,
or .602148191 of its initial value. For H=6 it is .303557973, or .701037110.
These are calculations under a second source-backed preparation, not a second
experimental validation or a new fitted decoherence rule.

Applying D(K/Q) only to Alice was equivalent on Omega but is NOT equivalent on
this filtered state. At the two examples above, the trace distances from that
incorrect reduction are .734274489 and .602082864. The source state determines
which clock statistics may be compressed; a scalar total count is not a universal
quantum noise description.

Heralding is a declared instrument operation with a recorded prior probability.
Its physical implementation/workload is not newly derived in this cell. No
future-history/post-measurement filtering is introduced.

## 8. Consequences for the research program

A finite quantum instrument can be studied vertically without first constructing
persistent moving observers. The native timing process now supplies a full
channel whose outcome predictions, state positivity, Bell score and entanglement
are mutually consistent and independently checkable.

This is not physical GR. Conditional on a tagged-service count K the primary
phase response is still exactly K/Q, independent of the workload arm. Native
workload affects when counts accumulate and their ensemble distribution; no
additional local energy/phase-rate shift per serviced tick has been derived.
The global type-blind sampler still does not generate different conditional mean
service rates for two selected single-token clocks.

We have not established an autonomous Hamiltonian, a localized lapse, a
continuum, an SI map, permanent loss of coherence, or a physical neutron. The
observer-identity issue is parked, not declared impossible or silently repaired.

## 9. Reproduction and validation

Notebook: paste the complete Cell32 after Cell29 (or after later cells). It reads
only dcu_mass_29 and leaves it unchanged. No earlier source import is re-executed.

Standalone: run `python reproduce.py` from the extracted bundle. Python3.10+
and NumPy are sufficient. It uses the frozen full Cell29 histories in reference.
The validation script additionally uses mpmath and SciPy, not the primary cell.

The main calculation evaluates 84 direct-source density matrices and 84 filtered
ones, with 6048 joint outcome probabilities in total. It analyzes 3072 original
histories at seven checkpoints and two phase depths. Counting those views twice
does not create extra native histories.

Independent validation replays six original conditioned events and 24576 complete
native transitions using a separate bitmask/sequential-charge implementation,
checks 144 archived rows, evaluates 6048 probabilities and all 168 states at
75-digit precision, and checks 16 phase generators with independent matrix
exponentials. See TESTS.txt for actual numerical differences and runtime.

A replay from a different directory containing spaces is compared field-for-field
with the delivered results. This is a Linux reproducibility check, not a Mac
benchmark. No numerical implementation tolerance is claimed as physical precision.
