# Cell 24 — finite record phases, collective readout, and energy-map identifiability

## Result in one paragraph

The existing, explicitly adopted complex record instrument specifies a finite-setting interference law. Its refinement structure can distinguish generators that would be indistinguishable at a single coarse resolution. Applying that same local phase/readout to the manuscript's declared shared-source state gives an exact collective phase-addition law. Neither result selects the signed graph spectrum or its spectral absolute value as physical energy: both admit refinement-compatible finite-clock representations, with different cavity gap ratios. This is a conditional algebra and implementation study, **not an empirical quantum experiment, a derived mass hierarchy, or a new native service law**.

## Source basis and what was held fixed

1. *Self-Limited Growth, Intrinsic Time, and Record Geometry in the Distinction Combinatorial Universe*, revised manuscript 0.5 (September 2026):
   - Section 5.1: independent first-write ternary values, entry reflection, persistent shared registers; cyclic/complex frame choices are explicit.
   - Section 5.4, Eqs. (45)-(46): `[d,k]=[d+1,3k]`, pairing with record addresses, and the Prüfer 3-group phase clock.
   - Section 7.1, Eq. (55): the diagonal-sector, trivial-character state `Omega_j = sum_r |r>^tensor j / sqrt(3)`.
   - Section 7.1, Eq. (56): the specified phase-family probability `L(x)=(1+2 cos(2 pi x))^2/9`, finite triadic settings, continuous extension and scope.
   - Section 3.1, Eq. (20): at least j simultaneous siblings from one unresolved source need Gamma >= j+1. This is not an arbitrary-exact-multiplicity statement.
   - Section 2.3, Eq. (9): all absent pairs among the co-served old objects form, not only selected desired pairs.
2. Actual Cell 23, including its phase-gate realization `D(x)=diag(1,exp(2 pi i x),exp(4 pi i x))`, inverse-Fourier readout, and signed/magnitude graph-energy countermodels.
3. Actual `_Native` and `_reference` loaded by Cell 1 from the original structural reproducibility archive.

The inherited mass, threshold, screening, and neutrino dictionaries are unchanged. No experimental target, new measured constant, or fitted coefficient is used. Existing quantum-state, complex-amplitude and Born-readout assumptions are retained; they are not derived from counting in this experiment.

## 1. Finite phase law

Let x=k/3^d modulo one. In the retained character frame,

    D(x)=exp(2 pi i x G),  G=diag(0,1,2).

This sign convention follows Cell 23's phase settings. If a later physical convention used `U(t)=exp(-i H t/hbar)` with x=t/T, then H would be `-(2 pi hbar/T)G` up to a scalar shift. We do not choose such a T or identify G with physical energy here. Common scalar shifts of G do not affect the measured probabilities.

Integer labels give exact identities:

    D(x+y)=D(x)D(y),  D(x+1)=D(x),
    D(k/3^d)=D(3k/3^(d+1)).

The last equality is a change of resolution for the SAME phase element. It is not a statement that an old tick becomes three newly occurred events or that physical time changes by a factor three.

Fourier preparation from |0>, D(x), and inverse-Fourier measurement give

    P_a(x)=L(x-a/3),
    L(x)=1/3+(4/9) cos(2 pi x)+(2/9) cos(4 pi x).

The two harmonics have a 1:2 ratio in the supplied phase coordinate. This is an exact consequence of the chosen three-character family, not an observed energy ratio or a new particle assignment.

The cell checks 363 (not necessarily unique) settings across depths 1..5 and 7,380 exact pair-composition congruences across depths 1..4. Numerical matrix probabilities are separately checked against the analytic character sum.

### What refinement resolves

On every setting x=k/9, both label vectors

    G0=(0,1,2),  G1=(0,1,11)

give the same gate because 11-2=9. At x=1/27 they give different gates and probabilities:

    G0: (0.964382751440680, 0.020201271814419, 0.015415976744902)
    G1: (0.211403426955145, 0.472910401990724, 0.315686171054131)

Their total variation at that predetermined first fine setting is 0.752979324486.

This is a **coarse-sampling alias**, not two models both satisfying the entire prescribed fine family. The fine source family rejects G1. The comparison does not assert that G1 is another physical state selected by the source.

More generally, two *fixed finite real diagonal label vectors* whose gates agree at every k/3^d have identical labels; for projective equality their differences can only be a common scalar. If a relative difference a persists, equality at x=1/3^d requires a/3^d to be an integer for every d, which forces a=0. This removes the branch ambiguity within the fixed finite diagonal representation. It does not uniquely reconstruct every operator realization from probabilities alone, nor provide an isometry into graph-mode coordinates.

This phenomenon has a standard analogue in stroboscopic quantum theory: one-period eigenphases determine quasienergies only modulo the corresponding phase winding. External context: Eckardt and Anisimovas, arXiv:1502.06477v4, Section 2.3, especially Eqs. (14)-(16). That paper is context for the distinction, not a premise of our arithmetic:

https://arxiv.org/html/1502.06477v4

## 2. Shared-source collective interference

Retain the manuscript's conditional source

    |Omega_j> = (1/sqrt(3)) sum_{r=0}^2 |r>^tensor j.

This state follows only with the source's diagonal-sector and trivial-character premises. Existence of j children alone does not imply that state.

Apply the same local D(x_i) and local inverse Fourier readout on each factor. The resulting tensor amplitude is

    A(a_1,...,a_j) = 3^(-(j+1)/2)
        sum_r exp(2 pi i r [sum_i x_i - sum_i a_i/3]).

Let s=sum_i a_i mod3 and X=sum_i x_i mod1. Then

    p(a_1,...,a_j) = 3^(-(j-1)) L(X-s/3),
    P(sum outcomes=s mod3) = L(X-s/3).

For j=2, this matches Eq. (56) after the explicitly tested orientation and outcome relabeling `x_1=x, x_2=-y, a_1=-a, a_2=b`. This relabeling is not a new phase convention fitted to a result.

For x=1/27:

| j | Minimum simultaneous-source capacity | P(sum=0), one local unit total | P(sum=0), one unit at every factor |
|---|---:|---:|---:|
| 2 | 3 | 0.964382751441 | 0.863205304967 |
| 3 | 4 | 0.964382751441 | 0.712386014201 |
| 4 | 5 | 0.964382751441 | 0.535003100154 |
| 5 | 6 | 0.964382751441 | 0.356870643158 |

The all-factor protocol gives L(jx), so scanning a common local setting gives a factor j phase accumulation. But that protocol applies j local setting increments per common step. Giving all j increments to one factor gives the same collective distribution. At a fixed TOTAL setting-increment count, distributing the increments differently does not change it. The cell checks all 494 allocations of seven finite 1/27 units across j=2..5 factors.

We do NOT use x/j as a convenient equal-budget setting: for j not a power of three, this can leave the allowed finite triadic catalog. Integer allocations avoid that error.

Each local outcome marginal is 1/3, regardless of the other settings. The coherent phase signal is a joint correlation; the calculation does not create signaling through local marginals. Erasing the coherence between the three equal-label components gives a classically correlated diagonal mixture. Independent inverse-Fourier readout then produces a uniform full joint distribution and hence the uniform three-outcome sum distribution. The classical equal-label correlation before Fourier readout has not been erased by pretending the density matrix is a product state.

The comparison is an operational tensor consequence, not a new Bell maximum, a quantum advantage benchmark, or an entanglement-to-mass identification. A count of j siblings is not a K_j clique and their Hilbert dimension 3^j is not j.

### Native preparation fixture scope

Four small ideal fixtures are constructed with the unchanged `_Native` pair constructor. First prepare a legal spine to j+1 objects; take its childless apex as the common source. Then condition on all those old objects being co-served. The actual burst includes **every absent pair among the served old objects**. All extra children not using the common source are kept. This supplies j siblings from the common source without suppressing other enabled births.

Setup work and entire burst work are recorded separately from the phase experiment. These are prepared-ideal, conditional-burst checks, not sampled autonomous runs from the Gamma-specific genesis. Their use does not establish a new preparation rate. In particular, no j>=2 simultaneous-source fixture is claimed to arise under Gamma=2. The experiment assigns no native recording cost to a phase setting merely because it is a unitary in the attached instrument.

## 3. Does the clock choose the cavity energy assignment?

Cell 23 compared the symmetric graph representative

    S=-D^(-1/2) A D^(-1/2)

with its spectral absolute value. On an isolated K_n, define P0=J/n and Pc=I-P0. Then

    S_n=-P0+Pc/(n-1),
    |S_n|=P0+Pc/(n-1).

Both are symmetric, share the structural projectors, and respect vertex permutations. The signed gap is 1+1/(n-1); the magnitude gap is 1-1/(n-1).

To test whether finite-clock compatibility excludes either model, take the same algebraic normalization beta=lcm(2,3,4)=12. Then both beta*S_n and beta*|S_n| have integer eigenvalues for n=3,4,5. Therefore

    U_n(x)=exp(2 pi i beta S_n x)

and its magnitude alternative both represent Z[1/3]/Z, are unitary, obey the refinement relation, and remain permutation symmetric. Beta is a convenient derived common-period normalization, not a fitted elapsed time, energy scale, or claimed service cost. It cancels from the gap ratios.

| n | Signed gap before beta | Magnitude gap before beta | Signed return probability at x=1/27 | Magnitude return probability |
|---|---:|---:|---:|---:|
| 3 | 3/2 | 1/2 | 0.333333333333 | 0.632732523408 |
| 4 | 4/3 | 2/3 | 0.311692070720 | 0.517448787733 |
| 5 | 5/4 | 3/4 | 0.379298361349 | 0.520000000000 |

The return calculation starts at a graph-coordinate basis vector. It is a mathematical discriminator for the candidate operators, **not a compiled permitted measurement on a DCU particle**.

The K5/K3 gap ratios remain 5/6 and 3/2. The magnitude choice agrees with the separate RMR beat prescription. These clock tests do not derive why that choice, rather than the signed one, describes physical energy. The local source phase family on a qutrit does not supply the missing identification of its record states with cavity common/contrast eigenvectors.

The same compatibility construction works for the actual registry's rational eigenvalue inventory. Multiplying by 9928 gives signed labels

    (-9928,-8432,766,136,73)

and magnitude labels

    (9928,8432,766,136,73).

Both obey exact finite-clock refinement. This reuses Cell 23's actual-reference spectral audit; it does not replace the registry by K_137 or claim these are physical energy levels.

**Scope of non-selection:** We exhibit two countermodels satisfying the stated spectral, unitary, finite-clock and permutation requirements. We do not claim both are compiled native gate words, or that they satisfy some unprovided full dynamical coupling condition. Finding such a condition could select or reject either model. The current test simply shows that clock compatibility alone does not make the choice.

## 4. Clock phase versus elapsed time

Every phase-clock element x has finite order: for some d, 3^d x=0 in the group. If a map t from that group to additive real elapsed time preserved addition, then

    0=t(3^d x)=3^d t(x),

hence t(x)=0. Thus no nonzero additive real-time map exists from the wrapped phase group alone.

This is not a prohibition on physical time emerging, and not a claim the DCU has no clock. The manuscript separately supplies the cumulative service counter tau. One must retain an unwrapped occurrence count or specify an event-to-phase and event-to-duration rule. The phase class by itself has discarded winding count. The model explicitly distinguishes scalar service count from the depth-and-orientation character clock.

## 5. What this adds and what remains missing

Established here, conditional on the existing representation:

- Exact control-parameter dependence and the local 1:2 harmonic structure.
- A concrete fine-resolution discrimination of a coarse generator alias.
- Joint multi-sibling phase addition with resource-count and phase-erasure controls.
- Two inequivalent graph-energy choices survive the same finite-clock constraints.

Not established:

- A quantum instrument derived from classical object counts alone.
- A unique physical graph Hamiltonian or a compiled K5 gate.
- A map from the scalar maintenance clock to the applied phase settings.
- A frequency in Hz, an action quantum, a rest mass, or a dispersion relation.
- Agreement with any experimental qutrit or multiparty interference data.

A next specific modeling task would be to define how a serviced, depth-tagged record changes an attached phase setting, with unchanged native budgets, and then predict fringe phase versus the retained service count. That would add an event-to-phase hypothesis explicitly rather than confuse a setting scan with autonomous evolution. Any physical frequency calibration would follow that step, not select the graph-energy map retroactively.

## 6. Executed validation and artifacts

The original Cell 1 archive was materialized into the working container, and actual Cells 1,14,17,18,19,23 were replayed. Cell 24 then ran twice with identical printed output and all serializable fields. Runtime in the single-BLAS-thread test was approximately 0.4 seconds for Cell24; this is not a runtime guarantee on other computers.

Independent checks:

- 702 probabilities from direct 80-digit complex sums; maximum absolute error 9.526e-16.
- Eight full dense Kronecker-product circuits versus the tensor-contraction implementation; maximum difference 4.441e-16.
- Six SciPy matrix-exponential checks of cavity return probabilities; maximum difference 3.886e-16.
- Four independently reconstructed complete enabled bursts, including extra children and sequential native work accounting.
- 27 common-generator-offset invariance checks.
- Rejection of nontriadic exact settings, rounded float settings, and a changed inherited gap prescription.
- Preservation of the prior Cell23 and registry dictionaries.

The working development version initially constructed only the desired source pairs in the fixture. Before the final scientific run, this was corrected to form ALL pairs enabled by the proposed service realization, with all extra children and their work retained. Phase predictions did not use fixture work and were not tuned by that correction.

These are implementation tests and mathematical consequences, not experimental uncertainties or global mechanistic validation.

Files:

- DCU_Mass_Cell_24.py: one complete cell, run after Cell23; NumPy only.
- DCU_Mass_Cell_24_output.txt: executed primary output.
- DCU_Mass_Cell_24_RESULTS.json: non-callable results, Fractions encoded as strings.
- DCU_Mass_Cell_24_TESTS.txt: independent test summary.
- DCU_Mass_Cell_24_NOTES.md: this document.
