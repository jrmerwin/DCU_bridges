# New vertical results for the bridges paper

Executed completion pass 01. This report describes model calculations and synthetic measurements, not laboratory observations. The two self-contained cells require Python and NumPy. They do not read or modify the Mac neutron campaign.

## 1. Record sharing becomes information measured in bits (Cell 38)

For a fixed native record path A, let Y_A contain its complete ordered entry readouts. The source's first-write law assigns an independent uniform trit g_z to each newly recorded register. A fixed entry map permutes that trit; the value is not redrawn on subsequent service. Conditioning on the structure and known entry conventions, the whole address therefore identifies its entire register list.

Using Shannon entropy in bits:

    H(Y_A) = |A| log2(3)
    H(Y_A,Y_B) = |A union B| log2(3)
    I(Y_A;Y_B) = |A intersection B| log2(3)
    I(Y_A;Y_B | Y_C) = |(A intersection B) minus C| log2(3).

Proof: each distinct observed register contributes one independent uniform trit. Entry maps and ordered encodings are invertible on the fixed register coordinates. Substitute the appropriate union cardinalities into the entropy and conditional-information definitions. The last formula is inclusion-exclusion, not a causal screening law.

The complete distributions of all 1,913 within-panel pairs of the 116 saved Cell26 paths were enumerated. A fixed structural representative selection yielded 46 three-reader distribution checks, with at most eight involved registers. No physical spatial dimension is involved. The conditional information may vanish even though two views share records: a third full record can already reveal those shared variables.

A fixed prefix-only readout loses information relative to the full addresses. All 36 equal-support/different-entry pair controls retain full mutual information, even though their literal addresses agree only one third of the time. Thus agreement of symbol strings is different from shared information with a known frame map.

Finite-value sampling is over independently prepared register-value realizations at the fixed graphs, not repeated writes in a single immutable history. At 65,536 samples, the four length-three H0 representatives give:

| Shared registers | Exact MI (bits) | Empirical plug-in MI |
|---:|---:|---:|
| 0 | 0 | 0.00825850 |
| 1 | 1.58496250 | 1.58727949 |
| 2 | 3.16992500 | 3.17020057 |
| 3 | 4.75488750 | 4.75467199 |

The small positive disjoint-support estimate is finite-sample bias, not a new interaction. All 4,096- and 65,536-sample results, including poorer sparse distributions for longer records, remain in the results.

## 2. The registry overlap graph is exactly a conditional-information graph

This is a declared extension from one path to a bundle of already recorded ancestral cards. Reconstruct the original complete 173-object size-bounded ideal. Its 137 roots with seven inclusive ancestors have four proper nonprimitive ancestral registers each. Exclude the root itself, whose record need not yet exist, and exclude the two deterministic primitives.

Let B_i = Anc(r_i) minus {a,b,r_i}, and Y_i be the ordered values of these four registers. Each proper ancestor has a child on a path to r_i, so its record exists in the source ideal. All roots contain c={a,b}. For distinct roots, no root is an ancestor of another because both have the same inclusive ancestry size. Therefore

    I(Y_i;Y_j | g_c)
      = (|Anc(r_i) intersection Anc(r_j)| - 3) log2(3).

The previously defined overlap edge, |intersection|>3, is consequently equivalent to strictly positive conditional information. The old threshold was not changed after seeing information values.

Exact distribution checks for every one of the 9,316 root pairs give:

| Conditional common trits after g_c is known | Root pairs |
|---:|---:|
| 0 | 3,969 |
| 1 | 3,774 |
| 2 | 1,323 |
| 3 | 250 |

The positive graph has exactly 5,347 edges and is K11 join (K63 disjoint-union K63). Removing the 11 universal vertices leaves two 63-vertex components. Under the retained structural sector classification their compositions are (S,I,G)=(9,2,0), (36,19,8), and (36,19,8). Sector labels are retained from full parentage; they are not derived solely from this information graph.

The two complete lobe bundles each involve 17 distinct ancestral registers and intersect only at c. Their complete values are conditionally independent given g_c in this independent-first-write model. This is dependence among copied records, not a spatial separation, communication rate, quantum entanglement, or causal barrier.

An important limit: only 37 distinct proper-register bundles occur among the 137 roots (25 bundles with five roots each and 12 singleton bundles). Different final parent pairs can yield exactly the same observed ancestral variable set. This readout recovers the registry's coarse dependence organization but does not distinguish every root or derive all sector labels. A single four-trit marginal is uniform for every root; the graph statement concerns joint correlated readouts, not single-shot root identification.

This is an operational information interpretation of the known overlap structure. It is not an independent proof of the uniqueness, physical necessity, or empirical validity of 137. The general information identity would apply to any similarly shared independent register collection; its instantiation here has the actual registry's graph and alias structure.

## 3. Native activity can be estimated from quantum readout outcomes (Cell 39)

The source apparatus and passive service-to-phase rule are unchanged. At one fixed pre-service state, capacity three, R=F+n, p=3/R, and p2=6/[R(R-1)], the two tagged tokens have total increment K with probabilities

    P(K=0)=1-2p+p2,
    P(K=1)=2(p-p2),
    P(K=2)=p2.

The fixed three-outcome readout at offset x and phase modulus Q is

    P(s|K,x)=[1+2 cos(2 pi (x+K/Q-s/3))]^2/9.

Use all offsets x=j/9, j=0,...,8, with Q=27 primary and Q=81 control. The resulting 27-by-3 mixture matrix has rank three. Analytically, the offset Fourier components reveal C0,C1,C2; the K=0,1,2 columns form a nonsingular Vandermonde matrix because their roots of unity are distinct. The full exact mixture is identifiable.

The finite estimator uses the narrower, already specified native one-parameter relationship p2=2p^2/(3-p). It maximizes the multinomial likelihood over p in [0,3/5]; the upper limit follows from the known five-object source preparation. Estimating p is fitting an unknown observable, not changing a fundamental coupling or calibrating a desired answer. The estimator receives only outcome counts and the frozen apparatus matrix, never the true pool, its split, actual service increments, state label, or future history.

Two fresh native trajectories use H=0/6 and Random(20380301). Their states are selected at post-preparation iterations 16,64,256, after the source's five genesis maintenance ticks. The preparation retains 26 recording-work obligations. All six states are retained; the preparation itself is a seventh development control. The six states are new relative to the archived tests, but dependent checkpoints of two trajectories, not six independent worlds.

Each synthetic dataset contains 32,768 ideal shots per setting (294,912 total), with 64 datasets per state/depth. Counts are sampled directly from the exact marginal probability law using multinomial draws; hundreds of millions of native histories were NOT simulated. Every shot presumes an independent replica of the same native/prepared quantum state. This is an ideal inverse-measurement experiment, not an implemented repeatable autonomous observer.

| New state | True R | True p | Median estimated p, Q27 | Relative RMSE in p, Q27 | Relative RMSE, Q81 |
|---|---:|---:|---:|---:|---:|
| H0, iteration 16 | 39 | 0.0769231 | 0.077150 | 2.88% | 9.85% |
| H0, iteration 64 | 51 | 0.0588235 | 0.058973 | 3.66% | 11.29% |
| H0, iteration 256 | 261 | 0.0114943 | 0.011749 | 10.46% | 30.79% |
| H6, iteration 16 | 9 | 0.3333333 | 0.333313 | 0.76% | 2.01% |
| H6, iteration 64 | 150 | 0.0200000 | 0.019712 | 7.62% | 22.82% |
| H6, iteration 256 | 385 | 0.0077922 | 0.007968 | 11.90% | 39.60% |

These RMSEs compare 64 synthetic estimates with hidden native truth, separately at each state. Empirical quantiles, all outcomes and all estimates remain available. They are not laboratory error bars, uncertainty in physical DCU constants, or a claim of exact recovery on every dataset.

No zero-rate estimate occurred, but the code retains the zero/infinite-pool possibility. Exact identifiability does not ensure useful finite-shot precision: low activity and the finer phase increment make the readout less informative. The full mixture matrix condition numbers are approximately 65.5 and 595.9. These matrix diagnostics and the nonlinear one-parameter errors are distinct quantities.

Only R=F+n, or equivalently p, is identified under this experiment. F and n separately have identical output laws whenever their sum agrees and the same instrument is available. The claim does not extend to reconstructing the whole ancestry, a local gravitational potential, a time interval in seconds, or an independently verified service law.

## 4. A third vertical connection: native recording work versus information

This secondary explanatory audit was added after the first Cell38/39 results, without changing their primary inputs or estimates. It uses every productive event among the same 512 new native preparation updates; no new search or trajectory was introduced.

For a completed event burst let W_new be the set of previously unread nonprimitive registers first used in that burst. Conditional on old written values and the native event history, the source uniform first-write law gives

    H(new register values | old values, native event) = |W_new| log2(3).

This measures the entropy in newly written register values, potentially acquired by a reader given access to those cards. It does not count all information in construction choices, full path factors, or event chronology. It is not thermodynamic entropy production.

The separate native charge is

    C = 2 sum_z L_z h_z + 9 sum_(z in W_new) L_z.

There were 21 productive events in the new continuation windows:

| New trits in event | Events | Native charge range |
|---:|---:|---:|
| 0 | 3 | 12–68 |
| 1 | 17 | 48–408 |
| 2 | 1 | 144 |

Thus the same amount of newly written trit information can incur very different construction work. Positive work also occurs without a new trit. The record-value count and the path-multiplicity workload are different inventories, as in the source definitions. This explicitly excludes a universal C = constant times new Shannon bits within this measurement; it does not rule out thermodynamic interpretations with additional energy, temperature and erasure definitions.

## 5. Implementation checks and limits

The main Cell38 enumeration checked 1,913 path-pair distributions, 46 chosen triple distributions and 9,316 registry-pair distributions, processing 14,635,002 finite assignment rows. Each assignment is a register-value possibility, not a native trajectory. The largest entropy arithmetic residual was 1.78e-14 bits; exact identities are proved symbolically, not inferred from floating-point equality.

An independent validator uses finite-field ranks of affine trit maps (the reflection is 1-g modulo three), not histogram entropy. It checked all 11,229 pair identities and all 46 selected conditional identities. A separate bitmask construction reproduces the source's full 173-object parent list.

The inverse experiment checked all 2,315 possible maintenance subsets across seven states against the exact marginal increment law, weighted by the relevant forced-token multiplicities. A separate replay uses the ORIGINAL Cell25 step and native constructor for all 512 new preparation transitions. A second density-matrix construction and SciPy optimizer reproduce all 896 rate estimates: maximum p discrepancy 1.21e-7, maximum optimized-loglikelihood discrepancy 3.50e-10. These values are numerical implementation differences, not physical precision.

Known limits remain: ideal readout access/repreparation, no native apparatus resource accounting, no SI clock or energy valuation, no independent experimental data, and no completion of late-catalogue queue/RNG audits. Primary source references are the revised DCU manuscript Sections 2 and 5–7, the registry structural manuscript and original code, Cells 25/26/32, and Shannon's 1948 definition of information entropy. The per-register independence model is an inherited assumption, not derived anew here.
