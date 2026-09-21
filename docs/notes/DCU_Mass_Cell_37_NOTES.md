# Cell 37 — joint-label compilation, and a separately scoped renewal alternative

## Result

A complete native implementation of the same-carrier nondemolition joint-label measurement has **not** been obtained. The existing factorwise quantum operations plus classical records, with no additional shared coherent resource, cannot implement that target. The obstruction is an operation-level statement, not the mistaken inference that outcome-blind workload forbids all state-sensitive measurement probabilities.

A different implementation is available at the adopted quantum-instrument level: local Fourier readout of the old memory, an independently prepared fresh equality-source pair, and outcome-conditioned phase encoding into that new pair. It realizes exactly the same labelled measurement instrument **on the protected code**, under a fixed identification of the old and new logical bases. It does not preserve the original physical carriers, is not equivalent outside the code, and does not transmit an arbitrary unknown logical superposition coherently.

Actual native bursts can create the required fresh common-source sibling pairs. Their complete construction charges are calculated here at two inherited epochs. Those charges are incurred recording obligations, NOT a derivation of all measurement, phase-encoding, controller, or factor-reattachment work. The latter quantities remain unspecified, not zero. No ideal controlled gate is declared executed merely because a birth has a calculable charge.

## 1. Sources and status of assumptions

The principal source remains *Self-Limited Growth, Intrinsic Time, and Record Geometry in the Distinction Combinatorial Universe*, revised manuscript 0.5, September 2026, supplied as `unified_deu_eq (1)(1).pdf`.

- Sections 2.1–2.4: parentage, recorded status, type-blind service, complete enabled bursts, and exact first/repeat work.
- Section 3.1, Eq.(20): at least two children of a common unresolved source require Gamma>=3. Every burst here uses Gamma=3.
- Section 5.1: persistent uniform first-write ternary values, with entry reflection; these are not arbitrary quantum output labels that may be overwritten.
- Section 7.1, Eq.(55): shared-source quantum preparation conditional on the diagonal-sector and trivial-common-character representation assumptions. Child formation alone does not prove those assumptions.
- Appendix A.2, Eq.(62): an upper bound on the entire future co-service intensity of two fixed objects.
- Section 9: the separate local factorwise edit compiler is not silently the counting-process workload.

Cell36 fixes the code, phase words, rank-one local Fourier implementation, and the coherent target projectors. The exact fixture and full source node table are retained in `reference`. Cell36's quantum apparatus is an attached representation, not an already compiled native device.

An additional source audit, `DCU_NATIVE_READER_DYNAMICS.md` (Library), Section 6, distinguishes source-backed factor-local query Kraus operators from the missing rule tying query occurrence, factor persistence and workload to global service. Its counting specialization is capacity TWO. We do not import that capacity or its phase-label convention into the capacity-three Cell36 experiment. It is corroborating evidence about the compilation boundary, not a substitute dynamical model.

Three levels remain separate:

1. **Unchanged native law:** all event probabilities, births, first writes, ledgers and maintenance counts in the resource audit.
2. **Inherited ideal instrument:** the shared-source state, the finite phase family and Cell36's chosen local measurement state update.
3. **New secondary composition:** use a successful resource burst as a readout opportunity; independently reset the new quantum pair to Omega_0; calculate the joint outcome k; apply D(k/3) to the new pair; designate it the output and attach the old storage driver. None of those trigger/reset/feed-forward/rebinding rules is inferred just from the native birth.

No geometry, physical location, energy, mass, SI unit, neutron identity, or autonomous observer is assigned.

## 2. Fixed mathematical target

Use omega=exp(2*pi*i/3) and

    Omega_k = sum_(r=0)^2 omega^(k*r) |r,r> / sqrt(3), k=0,1,2.
    C = span{|00>,|11>,|22>}.

The original four paths have supports

    A={4,7,10}, B={5,7,11}, C_path={4,7,11}, D_path={5,7,10}.

Both register and entry-token multiplicities balance. The two apparatus phases are driven by N_A+N_B and -(N_C_path+N_D_path). Thus every state in C is preserved by the fixed storage rule. Q=27 is primary and Q=81 is the existing control. The same storage connection to a NEW output pair is explicitly a reattachment premise in the renewal alternative.

The measurement target is

    I_k(rho) = Pi_k rho Pi_k,
    Pi_k = sum_(b-a=k mod3) |f_ab><f_ab|,
    |f_ab> = [sum_r omega^(-a*r)|r>/sqrt(3)] tensor
             [sum_s omega^( b*s)|s>/sqrt(3)].

Equivalently, with X|r>=|r+1 mod3>,

    Pi_k = (1/3) sum_(m=0)^2 omega^(k*m) X^m tensor X^m.

The three X powers are trace-orthogonal, Tr[(X^a)^dagger X^b]=3 delta_ab. Hence Pi_k has operator-Schmidt rank three. This is an exact rank certificate, not floating-eigenvalue thresholding.

## 3. Why the restricted same-carrier compiler fails

For product input |00>, which is inside C,

    Pi_k |00> = Omega_k/sqrt(3).

Every target outcome has probability1/3 and conditional negativity1. Each would turn a separable input into a maximally entangled two-qutrit state.

A completely refined branch of any finite protocol containing only factor-local operations, local ancillary systems initially in a product state, shared classical randomness and classical outcome feed-forward has product Kraus form A_h tensor B_h. Applied to a product input it yields a product vector; averaging unresolved branches gives a separable state. Such a protocol cannot implement the target.

This does NOT exclude an additional shared entangled ancillary source, a coherent cross-factor operation, quantum transfer, or a richer native instrument. In particular the source does provide a conditional shared-source state; the theorem does not erase that resource.

A stronger scoped bound follows for perfect local identification of the three Omega_k words. For any refined branch with reported word k, exact discrimination implies

    (A_h tensor B_h) Omega_j = 0 for j!=k.

Since |00>=(Omega_0+Omega_1+Omega_2)/sqrt(3),

    (A_h tensor B_h)|00> = (A_h tensor B_h)Omega_k/sqrt(3).

The left side is product. Thus every surviving output branch on the original two factors is product. Its squared overlap with maximally entangled Omega_k is at most1/3, by the Schmidt coefficients and Cauchy-Schwarz. Mixtures retain that bound. Cell36's local Fourier readout attains1/3; a different factor-local perfect discriminator cannot restore fidelity1 without an additional quantum resource.

These proofs are more useful than a bounded search over local gate words: they apply to every finite composition in the specified grammar. They are not universal prohibitions on native DCU quantum physics.

## 4. Logical renewal: measurement plus a fresh shared-source resource

The secondary route is intentionally a different target implementation:

1. Use the old local Fourier readout and retain outcomes a,b.
2. Compute the public label k=b-a mod3.
3. Supply a new independent Omega_0 pair from the adopted shared-source preparation.
4. Apply D(k/3) tensor I to the new pair, obtaining Omega_k.
5. Treat that pair as the new output memory. The old measured factors and their native history are not deleted.

The added quantum common-controlled-shift gate from Cell36 is not used. This does not eliminate the need for a physical readout, classical outcome communication, conditional encoding, independent source-reset rule and factor mapping.

As a channel from OLD nine-dimensional memory to NEW nine-dimensional memory, the fine Kraus maps are

    R_ab = |Omega_(b-a)><f_ab|.

Their adjoint products sum to identity. Summing over outcome pairs with fixed k gives the same POVM effect Pi_k on the whole input space.

On C,

    <f_ab|psi> = <Omega_k|psi>/sqrt(3), k=b-a,
    R_ab P_C = (1/sqrt(3)) |Omega_k><Omega_k|.

There are three pairs (a,b) per k. Therefore, for every rho supported in C,

    sum_(b-a=k) R_ab rho R_ab^dagger
       = |Omega_k><Omega_k| rho |Omega_k><Omega_k|
       = Pi_k rho Pi_k,

after the declared old/new logical-basis identification. The code checks all nine complex logical matrix units and also an input entangled with an untouched reference system. Thus this is not merely a match on three input probabilities.

An initial superposition of the phase words loses its inter-label coherence upon the measured readout, just as the target projective instrument does. For |00>, the unconditioned output is P_C/3, not |00>. Calling this coherent teleportation or unrestricted unknown-state copying would be incorrect.

The output difference outside C is retained. For |01>, the target output lies in span{|01>,|12>,|20>}, while the renewal output lies in C. The output density matrices have trace distance exactly1. Thus this is a code-restricted replacement implementation, not a realization of the complete nine-dimensional same-carrier measurement.

For each Omega_k input:

    label accuracy=1;
    old-memory fidelity=1/3, code weight=1/3, negativity=0;
    new-memory fidelity=1, code weight=1, negativity=1.

The fresh shared-source entanglement is a resource, not entanglement manufactured from classical records. Two successive ideal renewals return the same label on two successive replacement pairs.

## 5. Native resource construction and exact costs

Use the original sibling anchors(3,4), whose parents are(0,2) and(1,2). These were the apparatus anchors of Cells25–33. They are different from the endpoints10,11 of the four record paths supplying storage control.

At each of two inherited snapshots of H0_s20350920, select the smallest birth-order ID of a childless composite other than the two current anchors. This is a deterministic finite resource-selection convention, not a cost optimizer, spatial rule, physical identity, or preferred structural coupling. Select a legal service containing that source and the two anchors. All three slots service maintenance. Every absent pair among them forms.

The two children sharing the selected unresolved source are the proposed new output pair. Any additional enabled child is retained. Apply the source's conditional sibling-state representation to those two distinct factors only under the independent-reset premise above. Newborns never parent their siblings within the same burst.

The starting snapshots are:

| Snapshot | tau | Service iteration | n | F | P |
|---|---:|---:|---:|---:|---:|
| Early |128|814|19|547|4626|
| Late |4194304|649858187|10790|5503826|3890760512|

They are actual previously saved checkpoints, not freshly generated ages or a certified mature-particle epoch. The following four bursts are conditioned legal continuations, not unconditioned trajectories.

| Start snapshot | Round | Old anchors | Fresh common source | New sibling anchors | Other newborn | Repeat part | First-use premium | Total incurred work |
|---|---:|---|---:|---|---|---:|---:|---:|
| Early |1|(3,4)|8|(19,20)|none|56|54|110|
| Early |2|(19,20)|13|(21,22)|23|616|540|1156|
| Late |1|(3,4)|306|(10790,10791)|none|9432|4302|13734|
| Late |2|(10790,10791)|672|(10792,10793)|10794|21828|11394|33222|

The first round has only two births because pair(3,4) already exists as object9. In the second round, all three pairs are absent, so suppressing the anchor–anchor child would violate the constructor.

Every total is evaluated by Eq.(14): repeat charge2*L_z per event using the factor, plus a single extra9*L_z per first-written factor per burst. Ancestors are not double charged within an event, and a shared first-write premium is not charged twice across its simultaneous uses.

All three logical labels incur the same native construction charge for a given burst. The early two-round total is1266; the late total is46956. These are not universal memory-renewal costs or a theorem that cost grows monotonically with age. The selected fresh sources have different ancestries.

Since all three sampled requests are maintenance, forced service f=0. H=m=0 implies v=0. Thus each event adds its entire charge to F; none of those fresh obligations is discharged in that event. Each event adds three maintenance ticks and one service iteration; neither is a duration for the ideal quantum controls. Do not report the charge as work already completed, energy in joules, or a complete gate price.

The total work of the measurement protocol is UNKNOWN. In particular no zero is assigned for local queries, retaining/combining outcomes, conditioned phase encoding, source-state reset, controller or storage reattachment. The source distinguishes the edit compiler from the counting inventory, so adding an integer per matrix gate would be a new cost law rather than a derivation.

## 6. Readiness is neither forced nor guaranteed

For the first specified ready triple, the next-draw probability is

    1 / binomial(F+n,3).

It is1/30060260 in the early state and1/27950806457074173960 in the late state. All1160 early maintenance subsets are enumerated, including non-ready draws. The code does not force the chosen triple inside a claimed autonomous simulation.

These probabilities cannot be inverted to obtain a waiting time: native graph, workload and source eligibility change while waiting.

There is also a stronger exact availability limitation. The primary source paper's Appendix A.2 Eq.(62) bounds the entire future pair co-service intensity by

    B_Gamma(D) = 2*(D+A-1)/[(D-2)*(D-1)], A=choose(Gamma,2).

For Gamma3, this becomes2*(D+2)/[(D-2)*(D-1)]. Any future resource burst for this particular fixed two-anchor method requires those anchors to be co-served. The probability of ever receiving such a resource is no larger than their expected total co-services, hence at most this bound. The argument applies even if their pair already exists: Eq.(62) bounds intensity before the absent-pair indicator is inserted in deriving the source's later formation bound.

For the initial two anchors, the bounds are:

    early n19:   7/51 = 0.137254901960784...
    late n10790: 5396/29097933 = 0.000185442725433...

They are upper bounds, not measured success rates, and the actual resource event can be more restrictive. They remain bounds if the third source is chosen adaptively while the two anchors stay fixed. Moving/replacing anchors before readiness, allowing separate contacts, changing capacity or using an additional quantum interaction is outside this particular availability claim.

This rules out treating the displayed conditioned burst as a guaranteed repeatable memory-refresh service. It does not rule out other observer or interaction protocols.

## 7. A native trit is not a free quantum-output register

At first write the native raw trit g_z has uniform conditional law independent of the logical phase word under the passive counting attachment. With a uniform prior over k,

    p_native(k,g)=1/9.

A proposed substitution g=k gives

    p_substituted(k,g)=delta_(k,g)/3.

Both have identical uniform marginals, but their joint laws differ in total variation2/3. The first offers label-guess accuracy1/3; the second would offer1. Reproducing raw trit counts or a workload ledger does not establish that an implementation has transmitted the quantum label.

In the secondary route k is a separately declared measurement record used for feed-forward. It is not silently stored in or substituted for a native g_z. A native account of that outcome record and its conditional action is still needed.

## 8. Verdict and implications

- Same-carrier QND target: not compiled. Exact failure for independent-local operations with classical records and no new shared coherent resource.
- Code-restricted renewal instrument: exact conditional success, including correlated inputs with a reference and preservation of all three phase words on the new pair.
- Complete nine-dimensional instrument: renewal is NOT equivalent; retained trace-distance-one counterexample.
- Native creation of the fresh shared-source pair: allowed, with full construction costs, extra children and ledgers checked at early and late snapshots.
- Complete native cost, query scheduling, independent reset, feed-forward and factor reassignment: not derived.
- Autonomous repeatability: not established; the fixed-anchor resource supply has a restrictive exact all-future bound.
- Particle/mass/energy/space: not inferred. The result does not identify a neutron, and the neutron scanner is not accessed.

The principal lesson is not that records cannot influence physics. It is that a probability kernel, a quantum state-update rule, and a native resource ledger are distinct objects. They must be connected explicitly. The alternative shows that logical measurement need not require survival of the same carriers, while its resource and readiness limits prevent treating replacement as a free solution.

## 9. Execution and validation

The main cell is self-contained, uses NumPy and the standard library, and takes well below a second for this fixed calculation in the tested Linux container. It performs no new large stochastic growth run. It writes no source file, changes no earlier notebook state other than creating dcu_mass_37, and makes no network calls.

Validation uses the original DCUStructure class and the original Cell25 step, not only the new bitset implementation. It checks all1160 early successor states, all four complete resource bursts, both burst orderings for independent sequential first/repeat charge, and the complete source-node prefixes at both epochs. It verifies the four actual path entries against native cards.

An independent exact Q(omega) calculation checks243 code-channel entries. A75-digit amplitude construction checks243 target-operator entries and729 renewal-Kraus entries, plus effects and27 logical branch maps. A separate fresh-directory replay checks every discrete output exactly and float outputs to the specified numerical tolerance.

These are implementation and mathematical checks, not laboratory validation, proof that an ideal preparation is naturally realized, or completion of a native measurement compiler. The exact operator arguments carry their stated premises; numerical roundoff tolerances are not physical precision claims.
