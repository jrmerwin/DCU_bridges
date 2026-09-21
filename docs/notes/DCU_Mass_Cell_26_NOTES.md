# Cell 26 — record-supported clocks and an exact overlap response

## Result and scope

A fixed placement path now controls an attached clock through the maintenance of the recorded composites on that path. Equal-length reader clocks share their common-register increments exactly. Their conditional difference variance is proportional to the symmetric difference of their register supports. This is an exact consequence of the native global sampler **and a newly stated passive readout wiring**, not an inserted spatial metric in the native dynamics.

Independent continuations from two saved native checkpoints reproduce the predicted second-moment scale across record lengths three and four. The clock readout is nevertheless blind to entry orientation: two different token paths can have identical clocks while their ternary addresses agree only with probability 1/3. The resulting feature metric is generally not an ultrametric and has affine dimension greater than three in the tested panels. These are limits of this particular readout, not a dimension measurement of physical space.

No proton, neutron, quark identity, mature-particle epoch, nuclear binding law, SI scale, or physical propagation law is inferred.

## 1. Source basis and the new hypothesis

Sources are the supplied revised DCU manuscript 0.5, especially §§2.3–2.5, 5.1–5.4, 6.3, and 7.1; and actual Cells 24–25. Cell 25 supplies its unchanged `native_step` and `prepare` functions, the original `_Native` class, and two complete first-seed traces. Cell 24 supplies the unchanged two-sibling tensor readout.

The manuscript's placement record is a primitive-rooted directed path to a **recorded endpoint**. A token `(z,e)` identifies a recorded composite z and which lower/higher parent was entered in the fixed construction order. Its persistent ternary value is read through the prescribed reflection. This representation is not an assignment of spatial coordinates.

For a selected path R, define its register support A(R) as the distinct composite IDs along that path. The reader counter N_A counts all services to these particular maintenance tokens after the observation checkpoint. There is no extra weight for depth position, sector, endpoint, path multiplicity, or physical particle label. Each counted service supplies the SAME phase increment 1/27. The 1/81 control uses a genuinely finer new increment, not a rewrite of 1/27.

For each comparison of two reader paths, their counters drive the **original two-sibling apparatus** in opposite phase orientations:

    x_1 = N_A / 3^d,  x_2 = -N_B / 3^d.

The inherited state remains Omega_2=(|00>+|11>+|22>)/sqrt(3). The quantum result is

    P_s = L((N_A-N_B)/3^d - s/3),
    L(x) = (1+2cos(2pi x))^2 / 9.

We do NOT assert that arbitrary placement paths are themselves entangled siblings. They are two control inputs to the prepared apparatus. Register maintenance drives every included reader; a shared register therefore drives both phases with opposite signs. That control fan-out is an ADDED passive readout rule. It supplies no new workload, capacity, communication time, or compiled gate implementation. Persistent g_z register values are not overwritten; their entry reflection is **not** silently interpreted as a phase sign.

All comparisons are separate counterfactual terminal readouts using the same classical service ensemble. They are not simultaneous physical apparatuses with free resources or repeated projective measurements on one quantum system. Each individual history remains coherent under the inherited retention assumption; ensemble averaging is not a proof of irreversible decoherence.

The observation counters and attached phase begin at zero at the new checkpoint. That is a declared controlled observation/preparation convention, not deletion of the native object's past or clearing its work.

## 2. Exact one-step measurement law

Let T=F+n be the total token pool, q the number sampled, and

    p=q/T,  p2=q(q-1)/(T(T-1)).

For maintenance indicators I_z, E[I_z]=p and E[I_z I_w]=p2 for z != w. Set Y_A=sum_{z in A} I_z. Counting the equal-index and distinct-index terms gives

    Cov(Y_A,Y_B | current state)
      = (p-p2) |A intersect B| + (p2-p^2) |A||B|.

The second term is the common negative without-replacement sampling contribution. It is not spatial repulsion. For supports of the same size k, the expected count difference is zero, and the common term cancels:

    Var(Y_A-Y_B | current state) = a_t |A symmetric-difference B|,
    a_t=p-p2=q(T-q)/(T(T-1)).

For phase increments, multiply this variance by 3^(-2d). Where a_t>0, a rate-normalized squared clock distance is therefore

    d_clock(A,B)^2 = 3^(2d) Var(Delta x_A - Delta x_B | state) / a_t
                  = |A symmetric-difference B|.

This is derived from the response; it is not an extra target-dependent distance formula supplied to the simulation. Nevertheless, the **choice of support-only coupling** is what makes it a set-overlap readout. This law would hold for any equal-size fixed subsets under the same exchangeable sampling; it is not special to particles, K3, or the 137-entry registry.

For two size-k supports with shared size c, let m=k-c. There are m positive, m negative, and T-2m neutral token labels. The exact one-step law is

    P(Y_A-Y_B=z) = sum_{u-v=z}
        C(m,u) C(m,v) C(T-2m,q-u-v) / C(T,q).

The code independently checks this expression against complete maintenance-subset enumeration, weighted by the number of forced-token selections. There is no sampling approximation in those finite checks.

## 3. Accumulation through an evolving native environment

The supports are fixed, all their objects already exist, and their sizes agree. Thus D_t=N_A(t)-N_B(t) is a martingale. Its increments have conditional variance a_t delta, where delta=|A symmetric-difference B| is fixed. Orthogonality of martingale increments gives, at each fixed finite horizon,

    E[D_t^2] = delta E[sum_{k<=t} a_k].

The expectation is over fresh native continuations from the fixed checkpoint. The graph, queue, and future a_k depend on past random choices. We DO NOT assume independent identical increments, condition on an entire future queue path as though it were external, or replace this identity by a pathwise equality.

The primary estimator is a ratio of ensemble means:

    delta_hat = mean(D_t^2) / mean(sum a_k).

It is not the average of per-history ratios. Within an overlap class, pairs are averaged inside each history before estimating Monte Carlo uncertainty. This respects their substantial shared-history dependence. The reported ratio standard error is the delta-method estimate using the history-wise residual U_i-delta_hat A_i. It describes simulation uncertainty, not experimental error or model uncertainty. Two blocks of 128 seeds remain separately inspectable.

### A clock-only reference, without inspecting the global queue

The original two sibling maintenance counts N_3,N_4 supply a reference difference B=N_3-N_4. Its symmetric-difference size is exactly two, so

    E[B_t^2] = 2 E[sum a_k],
    delta = 2 E[D_t^2]/E[B_t^2].

This additional derived diagnostic removes explicit F and n from the ensemble ratio. It does not require an SI constant or fit a coefficient to the geometry. It does require access to unwrapped counts and a reference pair; a single wrapped fringe does not expose this information.

Finite-reference noise is retained. The observed reference second moments are 24.94140625 versus predictable 22.38769668 (H0 checkpoint), and 14.80078125 versus 15.90575944 (H6 checkpoint). Corresponding clock-only distance estimates and their correlated ratio errors are saved, not replaced with the cleaner global-ledger normalization. This is not yet a fully compiled internal observer able to collect all those counts.

## 4. Frozen finite domain and observed results

Restore the first stored native history of each Cell25 regime, seed 20260920, after exactly 512 post-preparation iterations. The state is replayed through the actual native transition and compared with its original ledger. Only then are paths selected; no future-born object or future recorded status is used.

The two checkpoints are:

| Existing regime | Objects | Outstanding work F | Relief supply P |
|---|---:|---:|---:|
| H0 | 18 | 691 | 2842 |
| H6 | 32 | 339 | 25 |

Their different states mean this is NOT an isolated causal intervention changing H while holding the environment fixed. It tests the same response theorem in two different native environments.

Enumerate all primitive-rooted paths with exactly three tokens ending at recorded objects (primary), and all such four-token paths (transfer):

| Regime | Record length | Paths | Distinct register supports | Exact affine feature rank |
|---|---:|---:|---:|---:|
| H0 | 3 | 17 | 13 | 11 |
| H0 | 4 | 17 | 12 | 9 |
| H6 | 3 | 42 | 29 | 18 |
| H6 | 4 | 40 | 26 | 17 |

Counts retain distinct token paths, including different primitive entries with equal register supports. These panels are not random samples of all possible records or physical species.

Run 256 fresh continuations, seeds 20280920..20281175, for each checkpoint, 2048 iterations per continuation. All 1,048,576 primary native updates are retained statistically; raw per-object counts are saved at 17 checkpoints. Original recording charges, complete enabled bursts, forced work, relief, and maintenance remain unchanged. No history is dropped based on output. The same integer seeds across regimes are paired RNG starts, not evidence of independent environments.

The primary three-register results are:

| Shared registers | Predicted delta | H0 estimate +/- MC SE | H6 estimate +/- MC SE |
|---|---:|---:|---:|
| 3 | 0 | 0 | 0 |
| 2 | 2 | 2.008653 +/- .073120 | 1.946924 +/- .054392 |
| 1 | 4 | 3.929349 +/- .137194 | 4.023937 +/- .136989 |
| 0 | 6 | 5.706082 +/- .209035 | 6.087670 +/- .215280 |

For four-token paths, the respective estimates for delta=2,4,6,8 are:

    H0: 2.067146, 3.927696, 5.915289, 7.718710.
    H6: 1.898019, 4.134614, 6.224181, 8.480956.

The full errors, seed blocks, and moment residuals are in the output and RESULTS. No estimated value is required to match the exact number to numerical precision; finite stochastic fluctuations are expected and retained.

At phase increment 1/27, the three-token paths' mean P0 values, ordered by shared count 3,2,1,0, are:

    H0: 1, .602864, .467630, .409071.
    H6: 1, .671112, .530541, .456326.

Perfectly shared support gives identical counters and P0=1 on EVERY history. The other values are empirical means, not universal functions of overlap alone at arbitrary times. Finite periodic readouts can revive; these finite-window monotonic means are not a general monotonicity theorem. Larger/deeper phase denominator 81 is retained as the previously specified smaller-increment control.

The observational checkpoint was chosen after inspecting candidate path availability at several existing snapshots, before constructing the quantitative clock comparison. No preregistration or blinded source selection is claimed. The clock-only normalization was added as a derived consequence of the exact moment relation; it does not alter the coupling, paths, histories, or primary estimates.

## 5. What the readout cannot identify

### Entry blindness

In the H0 length-three panel, the two paths

    root0 -> (2,entry0) -> (3,entry1) -> (6,entry1)
    root1 -> (2,entry1) -> (3,entry1) -> (6,entry1)

have the same register support {2,3,6}. Their counters are identical on every possible service history, not just similar on average. Their first token reads g_2 versus its entry reflection. Under the source's uniform persistent ternary law, their full addresses agree only when g_2=2, with probability 1/3. The source's entry information is therefore genuine information that this support-only clock discards.

Each of the four panels contains a corresponding exact alias. The code enumerates the finite register assignments; this register-value average is distinct from the service-history ensemble. No new g values are sampled during repeated service.

### Feature dimension is not spatial dimension

The squared clock distances are squared Euclidean distances of binary support-incidence vectors. Their exact affine ranks are 11,9,18,17 in the stated order. They are NOT candidate physical spatial dimensions. They describe how many Euclidean feature coordinates are needed to preserve this particular complete finite table. No eigenvalue cutoff, low-rank fit, or scaling dimension is selected.

The resulting distance is a pseudometric on token paths (different paths can have distance zero), and a metric after quotienting by equal register support. It is generally not the manuscript's common-prefix ultrametric. For example, one H6 triple has squared distances (2,2,4), violating the ultrametric requirement that the largest distance be attained at least twice. All ordinary Euclidean triangle inequalities hold.

Clock proximity here means shared driving events. It is NOT spatial adjacency, finite-speed interaction, a force, a signal travel time, or a unique reconstruction of the underlying DAG.

## 6. Stopping point and the neutron follow-on

This completes the first record-supported clock milestone: a specified passive measurement law connects actual evolving native service to a nontrivial, exactly understood structural response, with informative blind spots. It does not complete emergent space or a microscopic energy law.

The next user-prioritized investigation is paired proton/neutron candidate identification. Existing RMR mass prescriptions are not already native uud/udd classifiers, and the current native record notebooks do not provide a verified mature-proton epoch. A reproducible structural/state predicate must be fixed before measuring formation times. The screening paper's own light-baryon prescription gives m_n-m_p=(5/2)m_e; the older taxonomy has the different expression (137/54)m_e. These are distinct hypotheses, not values to switch between after inspecting a native candidate. Neither specifies a waiting time.

Future time measurements should distinguish structural first availability, first occurrence in the stochastic ideal, and a stated dynamical-maturity criterion. Once a proton-candidate checkpoint and neutron predicate exist, the response can be measured as a distribution of additional maintenance counts, keeping incomplete/never-realized cases rather than inferring a duration from the small mass gap. A later deuteron interpretation also needs a binding observable; co-presence or immutable shared ancestry alone is insufficient to define one.

## 7. Validation and reproducibility

The original Cells 01,14,17,18,19,23,24,25 were executed from the original structural archive before Cell26. Two fresh-kernel runs of the final Cell26 produced identical outputs and all saved serializable values. Primary runtime was approximately 8.5 seconds in this container, excluding prerequisite replays; user-machine runtime can differ.

Independent checks are recorded in TESTS.txt:

- All 116 selected paths independently enumerated backwards from recorded endpoints to primitives.
- 12,954 weighted maintenance-subset visits across both record lengths; every one of the 3,942 full covariance entries checked exactly.
- Exact centered-Gram ranks independently checked with rational-domain linear algebra; metric triangle checks and entry-alias agreement checks.
- 16,384 complete native transitions checked against an independent sequential-charge bit-mask implementation, plus 128 checkpoint count/compensator checks.
- 648 outcome probabilities from 70-digit two-qutrit complex amplitudes, with maximum absolute numerical discrepancy 4.44e-16.

During development, a flat NumPy reduction over many repeated pairs exceeded an overly tight roundoff check. Replacing it with per-history means followed by the ensemble mean reduced accumulated summation roundoff. No physical parameter, path, history, or probability formula changed. Validation runtimes were shortened by using exact rational-domain rather than generic symbolic rank elimination. These were numerical implementation changes, not scientific refits.

Artifacts: one Cell26.py (requires Cell25 in the notebook), output, compressed results, test report, and these notes. The bundle is not a standalone reconstruction of all prior cells. No old artifact, mass dictionary, library file, or model law is overwritten.
