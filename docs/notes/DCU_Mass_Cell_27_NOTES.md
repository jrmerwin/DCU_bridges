# Cell 27 — a constrained partner search, without naming a neutron prematurely

## Result and scope

The proposed late-recurrence exception does not occur in the frozen registry-overlap operator. All 3,906 same-lobe peripheral pairs remain exact eigen-differences at eigenvalue 1/73, including every GG and GS pair in that class. All equal-and-opposite peripheral dipoles have zero response on the eleven universal vertices, even when their endpoints occupy opposite lobes. No physical time is assigned to these powers.

A separate, newly specified structural partner search does find actual K3 pairs related by replacing ONE constituent with its primitive-exchanged counterpart. Of 109 correspondences in the 1,753-motif round-twelve K3 catalogue, 82 have matching complete support size and required construction grade; 41 are co-present in at least one archived host. These are structural correspondences, not quark states. Under the unchanged native repeat-work observable, every one of the 41 co-host pairs has 16 positive and 16 negative contrasts across the old 32 probes, and exactly zero context-mean contrast. No neutron or proton is identified by this experiment.

The source mass prescriptions remain untouched. In particular, the screening paper's delta(n,p)=5/2 electron units is not inserted into an operator, used as a numerical selection tolerance, or used to choose a candidate. The alternative historical 137/54 prescription is not substituted.

## 1. Source basis

- *Structure, Symmetry, and Motif Embedding of the 137-Object Distinction Registry*, manuscript 0.1, sections on the universal-junction decomposition and declared signed response. This is the source for the exact two-lobe graph, the physical scope of S/I/G labels, the 3,906 cancellation count, and the distinction between an overlap edge and an actual comparison witness.
- The original structural reproducibility archive: `DCUStructure`, the complete 137-entry reference, `reproduced_exact_observations.json.gz`, and `reproduced_five_family_full.json.gz`, loaded by original Cell1.
- Cell4's original native repeat-Q protocol, the five-letter chain family, and Cell5's fixed anchor apparatus values. No physical calibration is performed here.
- Cell26's conclusion that support-only clocks do not establish a particle identity and are blind to entry orientation.
- *Thermodynamic Screening Corrections*, section VI, as the separate paired mass hypothesis, not as a native state predicate.

External context, not an input to the computation: Borsanyi et al., *Ab initio calculation of the neutron-proton mass difference*, arXiv:1406.4088, describes the physical splitting in QCD+QED as the competition of electromagnetic and mass-isospin breaking. No QCD interaction, isospin matrix, quark charge, or empirical splitting from that work is inserted in Cell27.

## 2. Exact recurrence result

Let U,A,B have sizes 11,63,63. Their degrees are 136,73,73, and

    T = -Adjacency * inverse(Degree).

For i,j in the same peripheral lobe,

    T(e_i-e_j) = (e_i-e_j)/73,
    T^k(e_i-e_j) = (e_i-e_j)/73^k.

This follows because every vertex other than i,j sees identical outgoing columns from the two endpoints. At the endpoints themselves, the missing diagonal supplies the surviving opposite signs. Zero external response is NOT a zero vector or an element of ker(T).

The count 3,906 = 2*C(63,2) contains the following sector pairs:

    GG 56; GI 304; GS 576; II 342; IS 1368; SS 1260.

It is not a distinguished proton sector. The original structural paper explicitly presented it as a non-specificity control.

For i in A and j in B, define u_A,u_B as the unnormalized lobe indicator vectors. Then

    T^k(e_i-e_j)
      = (1/73)^k (e_i-e_j)
        + [(-62/73)^k - (1/73)^k] * (u_A-u_B)/63.

The vector is zero on U for EVERY k. Its external branch response is nonzero for k>=1. The antisymmetric mode has eigenvalue -62/73; it is not a small shift of 1/73. The response can distinguish same- versus opposite-lobe preparations, but that is a known mode decomposition, not a new neutron signature.

Within U, differences have eigenvalue 1/136, adding 55 more zero-external pairs. The total unordered registry pair count is

    3,906 same-lobe + 3,969 opposite-lobe + 55 universal-universal
      + 1,386 universal-peripheral = 9,316.

All 9,316 pairs are checked at k=1,2,4,8. There are 3,961 zero-external pairs at each of those powers. An integer matrix M=9928*T is used, retaining numerator vectors over denominator 9928^k. Nothing is called zero because a decaying floating value becomes small.

Changing which sectors are named in a preparation does not change eigenvalues of a fixed T. Changing state-dependent couplings, the observation domain, or the graph could produce a different operator, but its rule would need to be specified separately. A later host does not alter old inclusive ancestries or the fixed reference's overlap table.

## 3. Why the first structural partner map is deliberately not u/d

The reference decomposition comes from the two primitive branches:

    c={a,b}, a1={a,c}, b1={b,c}.

An object belongs to A when its inclusive ancestry contains a1 but not b1, to B for the reverse, to U when it contains both, and to 'base' otherwise. The global primitive exchange pi(a)=b, pi(b)=a extends recursively to every object. It preserves construction grade, ancestry size and path count. It preserves S/I/G classification on the complete registry.

A global relabeling is not a native weak-interaction event, nor is it an established quark-flavor operation. A/B here are literal structural branch labels. They are NOT electric charges, color charges, flavors, or physical locations.

The trial tests whether this existing doublet-like combinatorial distinction is a useful first partner relation:

1. Use every K3 in the saved round-twelve catalogue, not a selected small-response tail.
2. Select a core with branch multiset AAB.
3. Replace exactly one A constituent u by pi(u), retaining the other two constituents.
4. Require a distinct saved K3 with resulting multiset ABB and its actual three comparison witnesses.
5. Match the complete minimal support cardinality and required construction grade.
6. For the primary cohort, require both complete supports to occur in at least one SAME actual sampled host.

The test does not confuse overlap triangles with native co-parent triangles. An AAB triple is not a clique in the raw registry-overlap graph (it lacks cross-lobe edges). Native co-parent triangles CAN have such constituents because their comparison children are separate objects.

The exact-closure round index is an observation grade, not a maintenance epoch or a mature-proton age. The decoded cache is used only as a term lookup, never as an ambient universe. Co-presence is verified against the ancestor set of each actual reported host root.

## 4. Selection results and retained exclusions

- Saved round-twelve K3s: 1,753.
- One-port AAB-to-ABB correspondences found in the same catalogue: 109.
- Matched support size and required grade: 82.
- Of these, co-present in at least one actual host: 41.
- The other 41 matched pairs have no common host in this finite sample; this is not a proof of native co-occurrence impossibility.
- 27 correspondences fail the size/grade match and are retained as outside the primary matching stratum.
- Mirror-port lookup failures, duplicate-constituent failures, and absent saved partner failures are retained with reason codes.
- There are 19 full global-mirror correspondences among all 82 matched pairs; 12 among the 41 co-host primary pairs.

This is an exploratory, target-free finite-catalogue search, not a prospective particle prediction or a historically preregistered experiment. Counts embedded as assertions check replay of this exact archive. They are not physical fitting conditions. The co-host criterion is a common-environment availability criterion, not a selection by response sign.

## 5. Repeat-work protocol and transfer check

For each of the 134 distinct matched supports, use all 32 old five-letter a/b chain probes. Each probe has grade six, eight ancestry objects, and L=12. Construct each minimal support legally using the original constructor. Warm every nonprimitive support port, then perform two full repeat passes with fresh probe successors. Ports and successor orientation use Cell4's original conventions.

For support port w, the baseline-subtracted repeat cost is exactly

    q_w = 2 * sum_{z in Anc(w) minus Anc(probe)} L_z,
    Q(M,probe) = mean over all nonprimitive complete-support ports of q_w.

Each measured contact is an actual legal new pair. Baseline work is a counterfactual quote from the same pre-contact state, NOT a second committed charge. Warmup and probe-extension work are stored separately. No service history, queue waiting time, or physical action is inferred from this controlled experiment.

All 4,288 candidate/probe Q values agree with the original Cell4 measure body; all 39,584 port charges agree with an independent direct ancestry calculation. Two additional anchor regression fixtures give Q(271,A)=248/9 and Q(271,B)=308/9, without any mass conversion.

Because repeat charges depend on old ancestry and recorded status, adding passive background ancestry does not alter this Q once the candidate is fully recorded and the probe overlap is held fixed. A separate validation actually prepares both supports inside their common saved hosts, extending probes beyond host grade to avoid duplicates, and repeats both A/B arms. All 164 results equal the minimal-support Q values. These are full-background response checks, not autonomous histories.

Total primary committed nonprimitive support contacts (warm and repeated, including regression fixtures): 118,806. Total primary probe fixtures: 4,290. Extra constructor events for preparation and probe extension are not included in that contact count.

## 6. Responses and their interpretation

For every one of the 41 co-host pairs:

    Delta Q = Q(ABB) - Q(AAB)

is positive for 16 probes and negative for 16 probes, with no zero probe cases. Its average over the fixed 32-word catalogue is exactly zero. The average is a probe-symmetrization diagnostic, NOT a replacement calibrated mass observable, and the individual signs are retained.

Among the 41 size/grade-matched pairs without a shared saved host, the same mean is zero in 39 cases and is -7/5 and +7/5 in the other two. These are retained; they are not removed for failing to look like the other cases.

The first example is chosen by required grade, support size, and canonical parentage, not by closeness to a target:

    left archive 91:  core (3,4,9), branch AAB,
    right archive 17: core (3,4,10), branch ABB,
    replaced constituent 9 -> 10,
    complete support size 9, required grade 4.

In structural terms, with c={a,b}, p={a,c}, q={b,c},

    left core = {p,q,{p,c}},
    right core = {p,q,{q,c}}.

Each includes all three actual comparison-witness children. One of those witnesses happens also to be an object relevant elsewhere in the opposite support; roles are not identified merely by equal IDs. The pair is a full primitive mirror up to interchange of the two fixed constituents.

    left:  Q_A=148/7, Q_B=164/7,
    right: Q_A=164/7, Q_B=148/7,
    difference: +16/7 under A; -16/7 under B.

Thus this is a probe-orientation response, not a uniformly positive intrinsic splitting. It is not a neutron candidate with measured mass. Any positive scalar observable invariant under full primitive exchange is identical on exact mirror preparations; a fixed asymmetric external probe can still distinguish them.

## 7. Connection to the separate source mass hypothesis

The locked screening expressions imply

    r_p = 1836 + 20/137,
    r_n = r_p + 5/2,
    (r_n-r_p)/r_p = 685/503104 = 0.001361547513...

That 0.1361547513% relative increment is an algebraic consequence of the source prescription, not a tolerance chosen for this catalogue. It does not imply that a decay eigenvalue, a native work difference, or a recurrence amplitude must shift by the same relative amount. Such a relation would require an identified observable map.

The manuscript's separate K2 neutron correction does not define a native u/d state. The core/block and S/I/G labels must not be silently assigned that role. In particular, color-neutral confinement and electric neutrality should not be conflated.

## 8. Consequence for the next step

This cell narrows the search rather than identifying a neutron. The late-power overlap exception is ruled out for this fixed operator by exact algebra. A concrete one-port structural correspondence is available, but this first branch-as-type trial does not produce an orientation-independent native repeat-work splitting in the primary cohort.

A stronger follow-on needs a local state/type or interaction readout that differentiates two three-constituent states, is not just global primitive relabeling, and uses the same rules for both. Charge response and a native mass observable must be distinguished. Entry orientation and persistent records supply actual information beyond the overlap graph, but assigning them charge or energy remains a hypothesis that must be made explicit. Neither native quark flavors nor a beta-decay rule has been compiled by this test.

There is no mature-proton checkpoint and no neutron occurrence time in this cell. First availability, stochastic appearance and dynamic persistence remain different. A later deuteron claim requires binding, not mere co-presence.

## 9. Reproduction

Paste the complete `DCU_Mass_Cell_27.py` as one Jupyter cell. It can follow Cell26 but requires only Cell1's original constructor, registry and catalogue. It uses the standard library and no network. It leaves previous dictionaries unchanged and stores all matched pairs and measurements in `dcu_mass_27`.

The bundle contains the cell, execution output, exact serialized results, a partner-table CSV, this note, and the test report. It is not a replacement copy of all earlier notebook cells. Fractions in JSON and CSV are strings. The source archive remains an existing prerequisite.
