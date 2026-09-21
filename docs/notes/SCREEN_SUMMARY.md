# Fixed K3 trial-reference screen

**CONDITIONAL_Q_RATIO_SCREEN_WITH_STRUCTURAL_TRIAL_REFERENCE**  
Protocol: `K3_TRIAL_REFERENCE_V1`. Observation cut: **140,000,000 maintenance ticks**.

The declared `P_trial = h0:K3:4,5,6` is verified and has **X(P_trial)=17/1**.
All **467 K3 occurrences / 460 ancestry-plus-actual-edge types** were screened
against that single denominator. The nearest nonself occurrence is
`h2:K3:5,6,7`, with X=16 and rho=16/17, below the reference. The smallest
positive excess is `h3:K3:0,3,7`, with delta=29/119=0.24369747899...,
versus the frozen empirical target 0.00137841946. Its excess is about
176.79486 times the target, with signed empirical gap error +175.79486.
There is no close neutron-sized proxy increment in this finite trial.
No acceptance tolerance was added, and no neutron-absence conclusion follows.

## Declaration and measurement

This is the **new structural trial declaration**, not recovery of an identified
proton. The supplied `CODEX_TRIAL_REFERENCE_AND_RUN.md` is copied unchanged beside
these files. The earlier `REFERENCE_SELECTION_REQUIRED` report and structural CSV
remain unchanged. The reference was chosen from structural data before opening
current candidate Q values: 328 rows satisfy nonprimitive, pairwise-incomparable
cores; P_trial alone has minimum complete-support size10. The prescribed order
is (support size, maximum support grade, ancestry-edge identity, host, sorted
core IDs). Readability, Q, neutron-target proximity and strict-sign behavior did
not select P_trial. Earlier aggregate response findings were already known;
this corpus is discovery data, not a newly blinded independent confirmation.

Reference certificate, in host0, seed41000920, Gamma=3, m=0, H=0:

- Cores `[4,5,6]`; actual edges `[[4,5],[4,6],[5,6]]`; witnesses `[19,252,13]`.
- Complete support `[0,1,2,3,4,5,6,13,19,252]`: ten nodes, eight nonprimitive
  ports, maximum grade4.
- Local parents `[null,null,[0,1],[0,2],[1,2],[0,3],[1,3],[5,6],[4,5],[4,6]]`;
  local cores `[4,5,6]`, local witnesses `[8,9,7]`.
- Structural completion tau9657 / iteration179958; all witnesses readable
  tau32772 / iteration924112. Outstanding first-write premium: zero.
- Ancestry-edge type:
  `871f66ee5b05e0419f65f73f284af021d8c8b408452942b0107d6b88cfb85c23`.
- Core-term-set identity:
  `2d45f823ba26ea98829efb0f54a404933f336a4e1aff15f38b91b4199ebed63e`.

The native parentage, births, first-use records and values agree with the saved
raw checkpoint. Independent set/term evaluation gives repeat contributions
`{2:0,3:0,4:8,5:0,6:12,13:40,19:32,252:44}`: total136 / eight ports =17.
The existing `repeat_vector` and a legal local `committed_repeat_check` also
return17. The saved apparatus has `root_already_outside=true` for `chain_aaaaa`;
its original extension certificate agrees with recomputation. This local
measurement verification is not new stochastic native history.

The frozen measurement is:

```text
X(M) = saved_occurrence['response']['Q169'][0]
     = saved_occurrence['response']['certified_Q169'][0][0]  (mask true)
Q(S,P) = (1/N) sum_{v in S,v>=2} 2 sum_{z in Anc(v) minus Anc(P)} L_z
L_z = 2 Ch_z - 2
```

S is the set union of complete core ancestry and all actual witnesses. N counts
**every nonprimitive complete-support port**, not three cores. Q has native
repeat-recording-work units per port; it is not an SI mass. All ports are virtually
warmed/fully recorded; fresh repeat contacts subtract the alternative primitive-arm
quote. Setup, warmup, guard construction, first-write premium, division by90 and
species/spectral multipliers do not enter X. The exported `unresolved_surcharge_exact`
is the untouched-state outstanding premium, kept separate from X.

Use the original primitive frame **0=a,1=b**, original saved support orientation,
and zero-based index0 `chain_aaaaa`. Its blueprint is
`[null,null,[0,1],[0,2],[0,3],[0,4],[0,5],[0,6]]`, root7. Order is 32 lexicographic
five-letter chains, followed by137 registry probes. All169 names, full blueprints,
certification masks and the exact mirror permutation are in `TOP_DIAGNOSTICS.json`.
Mirror-only responses obey `Q_mirror[i]=Q_original[permutation[i]]`; aaaaa0 maps to
bbbbb31. Transporting both support and apparatus restores the original response.
The mirrored arrays are mathematical diagnostic constructions. A realized
structural mirror, where noted, does not establish a transported physical g-state.

Source implementation: `dcu_search/readouts.py` (`work_coefficients`,
`repeat_vector`, `probe_blueprints`, `extension_certificate`,
`committed_repeat_check`) and `scripts/closed_core_response_v1.py`.
Source paths in files are relative to the project root `DCU_NEUTRON_SEARCH_v1/`.

## Arithmetic, ranks and complete-population counts

```text
delta_empirical = 68920973/50000000000 = 0.00137841946
rho_empirical = 50068920973/50000000000
delta_screening = 685/503104
rho_screening = 503789/503104
rho = X(M)/17
delta = rho-1
empirical_gap_error = (delta-delta_empirical)/delta_empirical
screening_gap_error = (delta-delta_screening)/delta_screening
full_ratio_error = rho/rho_empirical-1
```

Rank uses **abs(empirical_gap_error)** only. Exact fractions drive all calculations;
28-significant-digit decimals are displays. Competition rank equals one plus the
number of rows with strictly smaller exact error. Equal errors share rank, with
occurrence ID as the display tie-breaker. The screening benchmark is a separate
column, never an alternative selected per object. All467 rows remain in the CSV.

| Count | Result |
|---|---:|
| K3 occurrences / ancestry-edge types | 467 / 460 |
| Usable certified primary values | 467 |
| Undefined / negative / zero X | 0 / 0 / 0 |
| Reference-self rows | 1 |
| Primary excess positive / zero / negative | 449 / 1 / 17 |
| All edge witnesses readable / unreadable | 390 / 77 |
| Complete-support ready / unready | 390 / 77 |
| Ordered two-orientation K3 Q signatures | 459 |
| Normalized K3 coefficient signatures | 460 |
| Repeated Q-signature groups / member occurrences | 5 / 13 |
| Q groups containing different coefficient signatures | 1 |

Readiness was checked against every complete-support record, independently of the
edge-witness flag, with zero disagreements. All77 unreadable occurrences stay in
the warmed-Q screen. Primitive-containing38, composite-only429,
ancestry-comparable139, pairwise-incomparable328 and witness/core-overlap12 rows
are retained. Reference eligibility restricts only P_trial, not this population.
These K3-only signature counts are distinct from the prior four-family totals.
Repeated types or signatures are not independent confirmations. Q-signature and
coefficient-signature membership counts and stable identities remain in the CSV.

The structural `ancestry_edge_type_id` and the response export's
`response_edge_identity_id` use different hash payload schemas. Both are retained
with separate names; occurrence IDs, actual edges and full certificates bind them.
A different hash encoding was not treated as a different physical structure.

## Five leading distinct types

The literal first five distinct types in the complete ranking include P_trial
itself at rank1. That row is explicitly **REFERENCE_SELF**, not a candidate lead.
The shortlist consequently contains four nonself types. All five are singleton
ancestry types and singleton exact Q/coefficient groups; no ties occur among them
or at the fifth-type boundary. The JSON includes every selected type member and
all of their exact Q aliases. Top-five inclusion is only a review convention.

| Rank | Occurrence | X exact | rho exact | delta exact | Signed empirical gap error | Signed screening gap error |
|---:|---|---|---|---|---:|---:|
| 1 | `h0:K3:4,5,6` (SELF) | 17/1 | 1/1 | 0/1 | -1.000000000 | -1.000000000 |
| 2 | `h2:K3:5,6,7` | 16/1 | 16/17 | -1/17 | -43.674621999 | -44.203434951 |
| 3 | `h0:K3:1,3,6` | 76/5 | 76/85 | -9/85 | -77.814319598 | -78.766182911 |
| 4 | `h3:K3:0,3,7` | 148/7 | 148/119 | 29/119 | 175.794862568 | 177.985659081 |
| 5 | `h1:K3:2,3,4` | 38/3 | 38/51 | -13/51 | -185.923361996 | -188.214884786 |

Exact error fractions and full-ratio errors are in the complete CSV. A self or
same-X comparison has delta0 and signed gap error-1; it is not a neutron-sized
splitting. Here the self row is the only exact equality. The closest nonself
row's negative increment is -1/17, and every positive increment is at least29/119.
Thus the absence of a near target follows from the fixed primary values, not
from rejecting context sign reversals or applying the old LP1 selector.

## Context diagnostics at the fixed calibration

For every diagnostic probe j, keep denominator17:

```text
rho_fixed_j(M) = Q(M,j)/17
ref_response_j = Q(P_trial,j)/17
paired_contrast_j = [Q(M,j)-Q(P_trial,j)]/17
```

The reference response itself ranges from9/17 to49/34 across169 probes. It is
not recalibrated to1 per probe. The following original-orientation summaries use
all169 certified entries; JSON additionally separates the chain32 and registry137
subsets and contains the complete original and mirrored raw vectors.

| Occurrence | Fixed candidate ratio range | Paired contrast range | Contrast signs (-/0/+) |
|---|---|---|---|
| `h0:K3:4,5,6` | 9/17 to 49/34 | 0/1 to 0/1 | 0/169/0 |
| `h2:K3:5,6,7` | 12/17 to 73/34 | -5/17 to 14/17 | 9/0/160 |
| `h0:K3:1,3,6` | 16/85 to 108/85 | -69/85 to 21/85 | 168/0/1 |
| `h3:K3:0,3,7` | 36/119 to 24/17 | -11/14 to 5/7 | 88/0/81 |
| `h1:K3:2,3,4` | 16/51 to 38/51 | -29/34 to 11/51 | 168/0/1 |

All four nonself shortlisted cases reverse their contrast relative to P across
contexts. These are diagnostics; none was removed or reranked for that behavior.
For mirrored candidate supports compared with the **original reference response**
at each probe, the JSON records the same fixed denominator and the following
paired-contrast ranges and counts:

| Mirrored occurrence | Paired contrast range | Contrast signs (-/0/+) |
|---|---|---|
| `h0:K3:4,5,6` | -31/34 to 31/34 | 81/7/81 |
| `h2:K3:5,6,7` | -25/34 to 55/34 | 30/0/139 |
| `h0:K3:1,3,6` | -213/170 to 63/85 | 88/0/81 |
| `h3:K3:0,3,7` | -11/14 to 15/17 | 88/0/81 |
| `h1:K3:2,3,4` | -29/34 to 11/51 | 168/0/1 |

The saved CSV orientation classes describe each object's contrast with its **own
mirror**, not its contrast with P_trial. The historical `NO_ROBUST_WORK_LEAD`
remains separately labeled and is not an entry gate in this screen.

## Structural relations and existing controls

Native IDs are host-local. Full semantic support intersections are given below;
only host0 intersections denote literally shared native objects. The JSON includes
all selected core/witness IDs, complete compact parent cards, original support
records and semantic added/removed/shared terms, without truncation.

| Occurrence | Relation to P_trial | Shared support/core/witness terms | Witness IDs | Completion / readability tau |
|---|---|---|---|---|
| `h0:K3:4,5,6` | Reference self | 10/3/3 | `[19, 252, 13]` | 9657 / 32772 |
| `h2:K3:5,6,7` | Cross host | 7/2/1 | `[10, 24, 9]` | 271 / 652 |
| `h0:K3:1,3,6` | Same host | 5/1/0 | `[6, 8, 9]` | 29 / 415 |
| `h3:K3:0,3,7` | Cross host | 6/0/0 | `[57, 12, 11]` | 596 / 1312 |
| `h1:K3:2,3,4` | Cross host | 5/1/0 | `[5, 6, 7]` | 19 / 35 |

- `h2:K3:5,6,7` shares two constituent terms with P, but replaces P's
  `{1,2}` constituent with `{0,5}` in its own host. The changed constituents
  differ in ancestry/grade and have no common immediate parent. This does not
  satisfy the unchanged LP1 substitution. Its selected core5 is an ancestor of
  core7; it remains in the candidate pool despite failing reference eligibility.
- `h0:K3:1,3,6` shares actual support IDs `[0,1,2,3,6]` with P. It retains
  core6 and changes the other two selected cores; its witnesses `[6,8,9]`
  replace `[19,252,13]`. It is a same-host structural comparison, not a
  demonstrated local conversion or prior LP1 partner.
- `h3:K3:0,3,7` is a cross-host structure sharing six support terms and no
  core terms with P. It supplies the smallest positive primary excess.
- `h1:K3:2,3,4` is a cross-host comparable-core structure. It is structurally
  self-mirrored and has zero self-mirror Q contrast. Its stored exchanged-core
  g values differ (core3 has2; core4 has0), so that structural symmetry does not
  establish symmetry of the actual stored state.

None of the four nonself cases is P_trial's exact primitive mirror or an existing
LP1 partner. No LP1 pair IDs are attached to these five occurrences. The exact
structural mirror of P and of the other first four listed types is not realized
in any of the four final hosts; the last type is its own realized structural
mirror in host1. All selected cores lie outside the defined137-root registry
domain and have no inherited registry roots in their core ancestry. That domain
status is explicit, not a missing value converted into a physical label.

All five are existing `OLD_K3_CONTROL` formation records and were readable before
10M. The screen concerns their presence at the140M cut; it does not assert that
these nearest rows formed in a late epoch. Their genesis event IDs and first-read
counters remain in the CSV/JSON. P's completion event is the single birth
252={4,6}. Its precise served set/pre-event queues remain unresolved. Persistent
g is sampled at first write and is unchanged by reader-frame choices; Q does not
depend on g. The records supply no observed candidate-to-P transition, u/d labels,
electric-charge assignment or native decay law.

Existing controls retained: complete support/witness certificates, primitive
strata, ancestry and overlap annotations, actually unreadable cases, genesis
controls, exact mirror diagnostics and archived LP1/strict-sign results. A physical
proton classifier, independent neutron confirmation set and native mass map remain
undefined; they were not added as prerequisites for this conditional screen.

## Validation, preservation and finite scope

- All467 response-file SHA256 values match the existing response index. Exact
  original primary Q was also recomputed from each saved local parent certificate.
- All157,846 K3 orientation/probe entries agree with their certification masks and
  exact saved vectors; mirror permutation is verified. No candidate array was rebuilt.
- CSV was parsed from disk:467 unique IDs, union identical to the frozen structural
  index; all ratios, both gap errors, full-ratio errors and competition ranks were
  recomputed from the exported exact fractions. All JSON vectors have169 entries;
  all type/alias/reference IDs resolve. Ties elsewhere in the complete CSV remain ties.
- The independently checked reference agrees with its raw checkpoint, independent
  ancestry formula and actual legal warmed repeat measurement. Baseline verification:
  **36 tests passed** (`python3 neutron_search.py verify`, 1.618s test runtime).
- The native model, physical prescriptions, old source CSV, paused checkpoints and
  response arrays are unchanged. No growth, LP1 scan, or full-service audit resumed.
  Read source hashes were unchanged across export.

Partial service coverage remains **14/16 inheritance certificates** and **0/14
planned late-milestone certificates**, with four planned200M sources unavailable.
Formula/certificate checks do not complete missing queue/RNG audits or establish
an SI/cosmological clock. No measured particle mass or neutron identification is
claimed. Target proximity selected the ranking and cannot also count as an
independent prediction of that target. The negative result limits this declared
reference and warmed-Q correspondence; it does not establish neutron absence or
exclude other separately specified baryon hypotheses.

This bounded pass is complete and paused for review. The full source catalogue,
response arrays and prior controls remain available locally. The exact declaration
is preserved beside the three returned files. Reproduction from those saved local
inputs uses `python3 scripts/k3_trial_ratio_v1.py` (writes the CSV and JSON);
its code hash and input hashes are recorded in the JSON. The summary supplies
interpretation and independent audit findings, not a new observable.

<!-- FILE_MANIFEST -->

Files are plain UTF-8; no ZIP was created. CSV and JSON hashes cover whole files.
The summary body checksum covers bytes before the manifest marker, avoiding a
self-referential full-file checksum. Sizes below include the entire delivered file.

| File | Exact bytes | SHA256 scope/value |
|---|---:|---|
| `K3_RATIO_SCREEN.csv` | 691465 | Full file: `9c0617d71942ea84dbafabc33b1cea32b4ab1186cdbdfb35a22c51dcc901fe75` |
| `TOP_DIAGNOSTICS.json` | 386651 | Full file: `ee50ebab02fcaa8d61d429c45b4d2c2a8787526e0089fbf44a9f9930522dc289` |
| `CODEX_TRIAL_REFERENCE_AND_RUN.md` | 12051 | Full file: `64ba18e204c321490e38f366a87c90c196469634dc257866a23f71ab17ffcc18` |
| `SCREEN_SUMMARY.md` | 17239 | Body (16190 bytes): `4ebe7fd87226a7057afdf7e9803ad1522ef1f46ec0afbb33127e453a284560c4` |

CSV has467 data rows; JSON has five complete diagnostic objects, each with two
169-entry vectors, plus169 probe blueprints. The three requested return files are
each below5MiB and together below10MiB. Full source/version hashes are in JSON.
