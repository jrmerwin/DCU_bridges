# Cell 34: exact age and observation-window audit

This is a new read-only summary of archived inputs. No new native simulation was run.

Clock convention: tau is the cumulative number of ALL serviced maintenance tokens since the two-primitive genesis. A service iteration is not a maintenance tick. The original apparatus preparation contributes five ticks; these are included below.

Cell 34 reused the two Cell 26 starting states and their 256 continuations each. The selected paths were frozen after 512 post-preparation iterations. Each continuation lasts 2,048 additional iterations. Paths were NOT reselected from the later populations. Cell 34 itself ran zero new primary native updates.

| Regime | Starting tau | Starting objects | Terminal tau mean | Terminal tau range | Terminal objects mean | Terminal objects range |
|---|---:|---:|---:|---:|---:|---:|
| no_relief | 120 | 18 | 407.07031250 | 286–590 | 34.83203125 | 27–47 |
| relief_H6 | 177 | 32 | 554.41796875 | 410–679 | 64.21875000 | 49–76 |

Additional maintenance ticks during the 2,048-iteration continuation:
- no_relief: mean 287.07031250, range 166–470.
- relief_H6: mean 377.41796875, range 233–502.

All ages are native computational counters, not established cosmological epochs, proper times or SI seconds. The sum of updates across independent histories is computational effort, not the age of one universe.

Exact versus finite-time scope: under the unchanged two-path coupling, equal nonempty supports have a protected subspace; distinct supports have only one-dimensional common character spaces when the singleton-register actions are allowed. This classification has no age parameter. Its necessity argument is available at any finite later checkpoint with Gamma=3 and F>=3. Age can change event probabilities, finite-window coherence, the structures available for selection and the behavior of larger/different encodings. Those questions were not answered by the small fixed panels.

Source: DCU_Mass_Cell_34_BUNDLE.zip, reference/DCU_Mass_Cell_26_RESULTS.json.gz and DCU_Mass_Cell_34_RESULTS.json.gz. Starting checkpoint records agree exactly in both files. Continuation totals also satisfy the saved service/relief ledger checks.

Cell26 input SHA256: `2897eb2418affb7381fc60781c7cb246d340d5af03e84862c6bccf4f0d50f545`
Cell34 input SHA256: `1cf984f225c69cf34daf7651fb442277b1643ab0117612b7ba1e5cd8e5d625e0`
