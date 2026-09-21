# Cell 25 — native maintenance to finite phase to interference

## Status

This is the first executed vertical service-to-instrument pilot in the current notebook. Native stochastic histories now supply the input to the inherited quantum phase/readout family. The event-to-phase coupling and phase-retention convention are ADDED hypotheses. The native constructor, workload and service law are unchanged. There is no new mass formula, experimental fit, spatial-dimension result, quantum-postulate derivation, time-to-second map, or inferred energy scale.

Source basis: *Self-Limited Growth, Intrinsic Time, and Record Geometry in the Distinction Combinatorial Universe*, revised manuscript 0.5 (September 2026), Sections 2.1–2.5, 3.1, 5.4, 6.3 and 7.1. The native class is the original `DCUStructure` in the structural reproduction archive. The instrument implementation and phase convention are those of Cell 24. The RMR screened mass and neutrino dictionaries are not used to fit or normalize this experiment and remain unchanged.

## 1. One explicit vertical hypothesis

Use the declared state

    |Omega_2> = (|00>+|11>+|22>)/sqrt(3)

on two siblings. The diagonal-sector, trivial-character, complex-amplitude and Born-readout assumptions remain inherited, not derived anew.

When the maintenance token of a tagged sibling is served, apply

    D(1/3^d)=diag(1, exp(2 pi i/3^d), exp(4 pi i/3^d))

to that attached factor. Untagged service, forced-work service, relief and subsequent births do not change this attached phase state. This phase-retention/spectator-state assumption is part of the hypothesis, not a theorem about every use of the record instrument. It does NOT alter the persistent ternary register values g_z. We are driving the existing phase-setting interface, not overwriting native registers.

The primary uses d=3 and the same positive orientation on both factors. The additional d=4 readout applies genuinely smaller NEW increments. Re-expressing the original increment at depth four would instead use 3/81=1/27 and leave its response unchanged; that identity is checked. Neither d=3 nor the three-outcome state is interpreted as spatial dimension.

A separately retained differential control uses the negative orientation on the second factor. No interaction between the factors beyond the declared preparation is inserted. The joint readout remains the sum of outcomes modulo three; a negative drive on factor two changes the accumulated phase to a difference of service counts.

No gate implementation or terminal readout cost is supplied by the native counting law. The phase readout is passive relative to counting: it does not insert workload, alter capacity, or draw random numbers from the native sampler. Its eventual physical implementation and resource cost remain open.

## 2. Preparation is actually reachable from genesis

For Gamma=3, m=0 and either H=0 or H=6:

1. Genesis has n=2,F=0,P=0. Both primitive tokens are selected deterministically; the child c={a,b} forms with cost zero. The cumulative maintenance count becomes 2.
2. At n=3,F=0 the three maintenance tokens are selected deterministically. All absent pairs form: {a,c} and {b,c}. Their joint workload is 26, not two separate first-use charges. The two siblings receive labels 3 and 4. Cumulative maintenance becomes 5.

Observation starts at

    n=5, F=26, P=0, prior tau=5.

The preparation's backlog is carried forward. No work is forgiven and no atypical preparation is selected. This common preparation is the same j=2 fixture as Cell 24, but here its service/genesis provenance is explicit. Two simultaneous siblings are not claimed for Gamma=2.

The attached phase is initialized on the newborn sibling factors, and observations thereafter use their maintenance tokens whether or not they are later used as parents. They remain fixed tagged objects while the universe grows.

## 3. The actual native process

Fix Gamma=3, m=0. Compare H=0 (no relief) and H=6 (the source relief cap at 2 Gamma). These are two existing constitutive parameter choices, not an inferred damping constant. Changing H affects the full feedback process: workload, births, population and later sampling. The comparison is not an intervention holding graph size and backlog fixed.

At each step, sample q=min(Gamma,F+n) distinct integer labels uniformly from range(F+n). The first n labels are maintenance tokens and the remaining F labels are individual outstanding work requests. If W is the selected maintenance set, then f=q-|W| is the selected forced count. Apply the source's rounded relief formula and construct EVERY absent pair within W. Newborns are not eligible during the same burst.

The original constructor computes the work. A separate explicit event/ancestry hit count checks it on each productive update. The ledgers satisfy, pathwise,

    F = 26 + cumulative construction work - cumulative f - cumulative relief,
    P = 2 cumulative f - cumulative relief,
    elapsed tau + cumulative f = 3 * observation steps.

There is no counterfactual suppression of unwanted births or deletion of an old object. The full source process runs regardless of the phase readout.

## 4. Fixed finite sample

- Two native parameter settings: H=0,6.
- 256 seeds each: 20260920..20261047 and 20270920..20271047, two fixed blocks of 128.
- Same seed labels reused across parameter settings as paired random-number starts. Their subsequent state-dependent draws diverge. Do not treat the two environments as independent datasets for an ad hoc combined significance.
- 4096 observation updates/history, no stopping or sample selection based on outcome.
- 2,097,152 total primary native updates, plus replay/validation histories.
- All 512 histories retained at 65 checkpoints, including zero-service periods. The first history per environment retains its full 4096-step event ledger and final parentage.
- d=3 primary and d=4 finer phase settings, same histories. Sum and difference orientations are separate declared diagnostics.

These are stochastic simulation replicates, not new blind experimental observations. No decay coefficient, phase offset, visibility function, time scale, or mass multiplier is fitted.

## 5. Exact one-service prediction

Write N=F+n, p=q/N and p2=q(q-1)/(N(N-1)). For fixed distinct tagged tokens a,b let I_a,I_b be their service indicators. Then

    E I_a=E I_b=p, E(I_a I_b)=p2.

For K=I_a+I_b,

    P(K=0)=1-2p+p2,
    P(K=1)=2(p-p2),
    P(K=2)=p2.

For D=I_a-I_b,

    P(D=+1)=P(D=-1)=p-p2,
    P(D=0)=1-2(p-p2).

At the initial observed state, N=31,q=3:

    sum law: {0:126/155, 1:28/155, 2:1/155};
    difference law: {-1:14/155, 0:127/155, 1:14/155}.

All 4495 individual-token subsets and all ten maintenance-token pairs are enumerated exactly. The pairs have parent–child distances one and two but identical service-count laws.

For an increment theta=2 pi h/3^d, its conditional harmonic multiplier is

    c_h^+ = 1-2p+p2 + 2(p-p2) exp(i theta) + p2 exp(2i theta),
    c_h^- = 1-2(p-p2)(1-cos(theta)).

These are exact CONDITIONAL one-step relations. F,n and future workload respond to sampling. We do not multiply a realized series of these factors and falsely identify it with the unconditional multi-step coherence. The experiment measures full-history phase distributions directly.

## 6. The ensemble quantum readout

Let N_a(t),N_b(t) be the cumulative tagged counts. The unwrapped phases are

    xi_+=(N_a+N_b)/3^d, xi_-=(N_a-N_b)/3^d.

Only their residues are used for trigonometric evaluation; unwrapped counts remain available. Each individual retained history remains a pure state in the adopted phase-retaining model. Terminal probabilities are

    P_s=L(xi-s/3), L(x)=(1+2 cos(2 pi x))^2/9.

Let C_h=E exp(2 pi i h xi). Then

    E P_s=1/3+(4/9) Re(C_1 exp(-2 pi i s/3))
              +(2/9) Re(C_2 exp(-4 pi i s/3)).

The ensemble density matrix in the equal-label sector has purity

    Tr(rho^2)=1/3+(4/9)|C_1|^2+(2/9)|C_2|^2.

The displayed C_h and purity are empirical estimates from the finite simulated ensemble, not exact population limits. |C_1| is a first-harmonic coherence measure, not by itself the full fringe visibility. Finite samples have a noise floor for such nonnegative magnitudes.

A low P_0 need not imply low coherence: a coherent phase can be near a fringe minimum. This is why both harmonic coherence and probabilities are retained. No destructive measurement is repeatedly performed during one history. Values at intermediate checkpoints are counterfactual terminal readouts on separately prepared realizations, not a sequence of projective collapses.

History-conditioned inverse setting operations restore the prepared phase state. This is checked with forward and inverse matrices, not only by evaluating L(0). It is a mathematical control using the simulator log and existing phase-gate family, not proof that an internal observer can recover or erase every relevant history, and not a costed native gate compiler.

## 7. Executed outcomes at 4096 updates

Means over 256 histories per environment:

| Setting | Final n | Final F | Elapsed global maintenance tau | Services to the original sibling pair |
|---|---:|---:|---:|---:|
| H=0 | 43.3867 | 1482.1758 | 594.9570 | 54.8320 |
| H=6 | 69.3047 | 2067.1250 | 702.3555 | 43.6016 |

The mean predictable tagged-service counts obtained by summing 2q/(F+n) are 55.3735 and 43.7253. Mean residuals are -0.5415 +/- 0.4668 and -0.1238 +/- 0.4220 service-Monte-Carlo SE. No run was removed for a large residual.

More relief permits more growth and aggregate maintenance in this finite panel, but the original pair receives LESS service. More aggregate activity is not the same observable as a faster clock attached to particular old objects. This is not a universal monotonicity theorem in H.

| Setting/readout | C1 magnitude | C2 magnitude | Ensemble purity | Mean P0 +/- MC SE |
|---|---:|---:|---:|---:|
| H=0, d3 sum | .141391 | .041633 | .342604 | .392601 +/- .022856 |
| H=6, d3 sum | .254287 | .092076 | .363956 | .246567 +/- .019741 |
| H=0, d4 sum | .811538 | .430716 | .667268 | .115603 +/- .010514 |
| H=6, d4 sum | .857293 | .540880 | .724990 | .072515 +/- .002596 |
| H=0, d3 difference | .262990 | .072512 | .365241 | .462174 +/- .023564 |
| H=6, d3 difference | .333787 | .024995 | .382989 | .482265 +/- .022959 |
| H=0, d4 difference | .855527 | .540516 | .723558 | .832706 +/- .012515 |
| H=6, d4 difference | .886113 | .615504 | .766497 | .863925 +/- .010220 |

The SEs quantify sampling over native histories. They are not experimental uncertainties and do not include model uncertainty or extra Born-counting noise. Seed-block summaries remain in the results; for H0,d3,sum their two P0 means are .35160 and .43360, illustrating finite-ensemble fluctuation rather than identical replication of a noisy curve.

At d3, sum orientation, replacing the distribution by its mean is markedly wrong:

    H0: E P0=.392601 versus L(E xi)=.975240;
    H6: E P0=.246567 versus L(E xi)=.027920.

No free phase-noise parameter produced this difference. It comes from the sampled native history under the declared event-to-phase coupling. It is not a demonstrated irreversible decoherence process or a derivation of thermodynamic screening residuals.

## 8. Derived global-versus-local clock diagnostic

After the first run exposed the aggregate/local difference, a diagnostic was added using an exact consequence of uniform service, without changing the primary drive or any histories:

    chi(t)=sum_{steps<=t} s_t/n_t,

where n_t is the PRE-service population. Conditional on n_t and s_t, the two tagged counts have a hypergeometric distribution with mean 2s_t/n_t. Thus

    E[N_a(t)+N_b(t)-2chi(t)]=0

at every fixed finite observation time. This is an expectation identity, NOT a pathwise equality or proof that the global activity average is relativistic proper time. It specifies the normalization rather than fitting one to the two regimes.

Results:

    H0: mean chi=27.608974, mean(K-2chi)=-.385917 +/- .422735;
    H6: mean chi=21.868351, mean(K-2chi)=-.135139 +/- .388869.

The exact rational chi per history is retained. The unnormalized global tau counts all objects' maintenance, so it does not have a universal fixed conversion to the phase count of two tagged carriers as n changes.

## 9. Spatial implication: a useful boundary, not a dimensional result

At a fixed state the one-step count kernel depends on F+n, q and the number of tagged objects, not their graph distance, age, sector, ancestry overlap or parentage. In particular,

    Cov(I_a,I_b)=p2-p^2=-q(N-q)/(N^2(N-1))

is the same for every distinct token pair. This is global without-replacement competition, not a spatial correlation function. A distance inferred solely from these equal-resolution one-step count-difference variances would assign the same separation to all distinct objects; it would not recover their parent–child metric.

This DOES NOT prove that all multi-step native responses are geometry blind: later graph/workload evolution can carry architectural information. Nor does it preclude other operational spatial observables. It tells us that instantaneous bare-token timing alone cannot identify the spatial dimension.

A future hypothesis can couple readers to actual shared record paths/supports, for which overlap enters the source covariance law. That needs a fixed support/readout definition and fresh tests. It is not silently included here.

## 10. Validation and files

The final cell was run twice with identical serialized outputs. Independent tests include:

- 16,384 complete transitions from a separate bitmask-ancestry/event-cost implementation, four histories across both seed blocks and both H settings.
- 2048 full dense nine-state circuits and eight ensemble density matrices, matching probabilities to 1.998e-15 and purity to 1.221e-15.
- 6144 independent 70-digit amplitude-derived probabilities, matching ensemble means to 6.661e-16.
- Exact initial token-subset enumeration, native ledger identities, all-enabled-burst retention, exact refinement identities, and explicit rejection of a changed inherited generator.

These validate implementation, not empirical physics. The scalar CPU run took about 24 seconds here; other machines can differ. No incomplete history was discarded.

Files: the single Jupyter cell, primary output, compressed results (all checkpoint histories and two full traces), independent test output, and these notes. No prior dictionary is overwritten.
