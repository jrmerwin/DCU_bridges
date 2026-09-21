# K3 successor investigation: results for review

**Completed all three saved-data searches. Two near-gap structures merit
structural review, but no object is established to meet the combined small-excess,
free-state instability and decay-product criteria.** These are structural proxy
leads, not identified neutrons. The fixed P_trial comparison remains negative
for a near neutron-sized excess; the closer results use explicitly separate
comparisons with each assembly's own contained K3.

## The useful structural leads

The most relevant later-forming example is in host3:

- K3 cores **[1823,30374,42163]**, plus adjoining core **14213**.
- Four-core motif: a triangle with a tail (paw). The extra core connects only
  to core1823, through actual witness **57446={1823,14213}**. The triangle's
  three actual witnesses are46339,46340,46341.
- Complete support grows **603 -> 619 objects**, adding16 nonprimitive ports.
  This is a real support enlargement, not merely a new choice of core labels.
- Frozen warmed-Q: **22378700/601 -> 23005256/617**.
  Local fractional excess is **4625239/3451914475 = 0.001339905444789...**,
  **2.79407% below the empirical excess** and1.58952% below the separate
  screening-paper excess. Percent errors here refer to the SMALL EXCESS,
  not to the full candidate/reference ratio.
- The K3 completes at **22,309,822 ticks**; the adjoining assembly completes
  later at **31,257,756 ticks**. Thus it has the requested structural chronology.
  Core14213 already existed; the later event creates its joining witness.
- The original K3 becomes fully witness-readable at40,841,222 ticks, after
  attachment. The new witness57446 is still unreadable at 140M. Its separate
  first-write premium is84,780 native work units, excluded from Q.
- The four-core structure remains structurally present for 108,742,244 further ticks to the
  final cut. No detachment, short free lifetime, or free/bound-state comparison
  is demonstrated. Unreadability is not evidence of instability or decay.

The closest numerical local excess is host1 K3 **[12282,41181,51044] + core11959**.
Its tail connects only to12282 through witness17425. Support grows900 -> 941;
Q changes27652988/449 -> 19304328/313. Its excess
**3064507/2163846311 = 0.001416231358217...** is **2.74313% above** the empirical
excess (4.01630% above the separate screening excess). However, its tail edge
already exists at5,115,561 ticks and the K3 and whole assembly become complete
together at57,022,396. It is not a demonstrated later successor of a pre-existing
K3. It is fully readable by82,910,364 and persists structurally to140M.

All five closest local-excess examples are retained below and in
[STRUCTURAL_LEADS.csv](STRUCTURAL_LEADS.csv). The CSV also contains five
fixed-P_trial comparison controls. Selection for these review examples uses only
the already frozen exact error order; no acceptance tolerance was introduced.
Complete local parent/witness certificates for all ten examples are in
[LEAD_CERTIFICATES.json](LEAD_CERTIFICATES.json), without truncated support lists.

| Contained K3 + extra core | Added support | Local fractional excess | Empirical small-gap error | K3 / assembly completion ticks | Chronology |
|---|---:|---:|---:|---:|---|
| `h1:K3:12282,41181,51044+core:11959` | 41 | 0.001416231358 | +2.74313% | 57,022,396 / 57,022,396 | Coformed |
| `h3:K3:1823,30374,42163+core:14213` | 16 | 0.001339905445 | -2.79407% | 22,309,822 / 31,257,756 | Later assembly |
| `h1:K3:7787,8256,12981+core:825` | 6 | 0.001301770958 | -5.56061% | 3,981,504 / 40,503,508 | Later assembly |
| `h3:K3:8,12,191+core:317` | 9 | 0.001460970846 | +5.98884% | 9,931 / 290,154 | Later assembly |
| `h0:K3:13046,96868,98216+core:12961` | 22 | 0.001272237727 | -7.70315% | 129,755,357 / 129,755,357 | Coformed |

All five have unequal attachment-incidence patterns, either[1,0,0] or[0,1,0]
on the sorted triangle cores. That is a specific structural asymmetry; it does
not imply chirality, a quark/flavor assignment, or absence of graph automorphisms.
The first, second and fifth examples also illustrate why catalogue presence at
140M, exact formation time, and witness-readiness time must remain separate.

These local leads do NOT repair the fixed trial-reference screen. For example,
the first two assemblies have Q/17 about3627.951 and2193.274, respectively.
The small local ratios compare a structure to its own contained high-Q K3; no common
mass calibration makes them near-unit ratios relative to the declared P_trial.
Treating all K3s as one physical proton species would require an additional
justified mapping that this calculation neither supplies nor fits.

## 1. Complete adjoining-component test

Frozen domain: each final140M K3 plus one distinct same-host native core adjacent
to at least one of its three core vertices. This covers every induced paw,
diamond and K4 with a contained K3; larger attachments are outside this finite
pass. Full support includes every core ancestor and every actual edge witness.

| Coverage or result | Count |
|---|---:|
| K3 anchors surveyed | 467 |
| Anchors with at least one adjacent fourth core | 466 |
| Explicit zero-neighbor control | 1 |
| Anchor/attachment relations | 14,304 |
| Distinct four-core assemblies | 14,203 |
| Paw / diamond / K4 relations | 14,108 / 184 / 12 |
| Enlarged / unchanged complete supports | 14,037 / 267 |
| Positive / negative / zero local Q changes | 12,938 / 1,097 / 269 |
| Assembly completes strictly after K3 / coforms | 9,835 / 4,469 |
| All witnesses readable / unreadable at 140M | 12,352 / 1,952 |
| Certified primary measurements / resource censors | 14,304 / 0 |
| Exact empirical or screening target matches, local or fixed P | 0 |

The surveyed zero-neighbor triangle is
`h3:K3:20269,62765,100361`: each of its three cores has only the other two as
co-parent neighbors at 140M. It contributes no invented attachment row and is
retained in [ANCHOR_COVERAGE.csv](adjoining/ANCHOR_COVERAGE.csv). Its old Q is
5389768/87; support size 872; structural completion95,755,892; witnesses unread.
The branch's original report overstated responding-anchor coverage by one;
[COVERAGE_ERRATA.json](adjoining/COVERAGE_ERRATA.json) preserves the correction,
and the [clarified report](adjoining/REPORT_WITH_COVERAGE_CLARIFICATION.md)
is the report to use. No numeric result or original certificate changed.

The 92 diamonds each give two contained-K3 relations and the three K4s each give
four. These overlapping relations are not independent objects or confirmations.
All267 unchanged-support cases necessarily have unchanged Q. Two additional
support enlargements happen to have equal normalized Q. Added unnormalized work
is nonnegative, but Q is a mean: added ports can lower it. Neither that accounting
sum nor the first-write premium was substituted for the frozen observable.

The closest fixed-P_trial result is `h1:K3:2,3,4+core:5`:
Q=120/7, rho=120/119, excess1/119=0.00840336134..., about6.09637 times the
empirical target. Its empirical small-gap error is+5.09637. All59 assemblies
that specifically contain P_trial itself are retained in
[P_TRIAL_EXPANSIONS.csv](adjoining/P_TRIAL_EXPANSIONS.csv):58 paws and one diamond,
all enlarged supports,57 positive and two negative Q changes.

### Unchanged measurement and distinct denominators

```text
Q(S,P) = (1/N) sum_{v in S,v>=2} 2 sum_{z in Anc(v) minus Anc(P)} L_z
L_z = 2 Ch_z - 2
X(M) = Q(M,chain_aaaaa), original frame0=a,1=b
N = ALL nonprimitive complete-support ports
local_delta = [X(assembly)-X(contained K3)] / X(contained K3)
fixed_delta = X(assembly)/17 - 1
empirical_delta = 68920973/50000000000 = 0.00137841946
screening_delta = 685/503104
signed_gap_error = (delta-target_delta)/target_delta
```

P_trial stays `h0:K3:4,5,6`, X=17. The local denominator is fixed by actual
structural containment, not an all-reference/all-pairs fit. It is a separately
specified structural diagnostic, not a replacement proton reference. The two
benchmarks remain separate. Exact fractions, all signs, primitive/ancestry/role
overlap strata and unreadable supports are retained. No strict-sign, LP1, charge,
spin or observed-decay requirement was added as an entry gate for this screen.

The original virtually warmed, baseline-subtracted mean repeat work remains
unchanged, with no division by90, new mass factor, weighted probe average or
first-use contribution. Existing closed-family responses are reused and checked;
new paw values apply the SAME already certified operator to saved parent cards.
Only aaaaa and its bbbbb primitive-mirror transport control were evaluated for
new paw supports. No new169-context sweep, probe search or native history was run.
The mirror value is diagnostic and did not choose orientation or rank.

## 2. P4/C4/V transitions and proposed decay residues

The fixed10M cohort contains **2,989,672 connected core sets** of size2–4.
All seven saved cuts through140M are used. The existing exhaustive counters
retain the denominator; an independently checked reverse join recovers every
changed core identity without rebuilding the three-million-set atlas.

| Initial -> final induced topology | Changed core sets |
|---|---:|
| Open V -> K3 | 22 |
| P4 -> triangle with tail | 234 |
| P4 -> C4 | 128 |
| Three-leaf star -> triangle with tail | 138 |
| Triangle with tail -> diamond | 1 |
| Total | 523 |

All **759 initial C4s remain C4s**. The523 changes reduce to only **147 actual
single-child birth events**. The234 P4->paw cases reduce to **22 new K3s/events**,
the same22 initially-open-V closures;17 of those K3s are fully witness-readable
by140M and five remain unresolved. Overlapping subgraphs are not234 separate
candidate decays or independent confirmations.

The proposed structural precursor channel does occur: a P4 loses its induced
P4 label when a missing edge is added, leaving a K3 plus the fourth core. In all
**374 K3-plus-fourth-core relations**, that fourth core remains directly attached
by at least one actual edge. In 32 relations it was already in the K3's ancestry
support. There are **zero recorded detached residues, node/edge removals or
splits of an initial component** in these records.

New comparison-witness objects are initially isolated co-parent vertices because
they are born after the served old cores were selected. The full co-born lists
are preserved. That is native object creation, not demonstrated detachment or
emission of an electron/antineutrino. In 435 changed cases an all-witness-ready
flag becomes false because a new unread witness is added; no older recorded bit
was erased. This distinction prevents a bookkeeping change from being mistaken
for a decay signal.

Append-only substrate history alone does not forbid a physical interpretation
in which an active pattern decays or is replaced. However, these catalogues have
no independently defined active-particle predicate, emission rule, physical
binding classifier or free/bound-state lifetime comparison. Those stronger
criteria remain **unavailable**, not measured false for all DCU mechanisms.
Graph components can grow and merge; their connectivity is not automatically a
nuclear binding law. The saved cycle-space quotient is a GF(2) vector-space
quotient, not a spatial metric with defined locations or velocities. It cannot
supply “elsewhere” or “fast” residue motion by itself.

The temporal result concerns core sets already connected at 10M. Newly connected
or later-born sets outside that cohort are not silently treated as tested decay
controls. The separate final140M adjoining scan includes such later structures.

## 3. Two coefficient-distinct Q-alias groups

All **five members** and four pair comparisons are retained in
[ALIAS_MEMBERS.csv](aliases/ALIAS_MEMBERS.csv) and
[alias DIAGNOSTICS.json](aliases/DIAGNOSTICS.json), including both full169-entry
vectors, coefficient maps/differences, actual witnesses and complete small cards.

1. **C4(0,44,51,57) and C4(1,44,51,57), host0.** Their32-object supports
   overlap in30 objects and exchange witnesses249/315 for188/48. The four
   differing normalized coefficient terms are grade 11/12, absent from every
   frozen probe and have zero total difference. The two ordered orientation
   vectors match exactly. There is no K3 among either square's four selected
   cores, although each support contains six complete ancestral K3 certificates.
2. **K3(0,3,44), K3(1,3,44), and diamond(0,2,3,44), host0.** The two triangles
   exchange witness315 for48; the coefficient difference is ±87/5 on these
   grade12 terms, again invisible to the169 probes. The diamond instead shares
   its ENTIRE22-object support and coefficient map with K3(0,3,44). Its extra
   core2, two edges and witnesses3/10 were already present in that support.
   The extra graphlet roles add no Q excess, and closure/readability times agree.

There is **no measured orientation mismatch between alias members**. The member's
own original-versus-mirror contrast can be nonzero, but it is the SAME contrast
for every member of a group. Passive g differences do not enter Q. All five are
fully readable early structures that remain recorded at 140M, with no observed
free-state instability or isotope conversion.

For the differing coefficient maps, a separately prescribed legal probe reaching
the unresolved grade 11/12 terms could test the finite-apparatus degeneracy. It
was not constructed or used to improve a mass match here. For the identical-support
K3/diamond pair, no new probe under this same support-only Q can distinguish them;
a role-sensitive observable would be a new hypothesis requiring its own definition.
The aliases are useful structural/readout controls, not an additional measured
mass or trapped particle degree of freedom.

## Criteria assessment and validation

| Requested signature | Evidence from this pass | Status |
|---|---|---|
| K3 plus a small asymmetric adjoining motif | Complete paw/diamond/K4 census; actual extra support quantified | Structural examples found |
| Neutron-sized positive excess | Near local gaps with 2.74–2.79% small-gap errors; fixed-P result remains far | Conditional local leads only |
| A later successor to an existing K3 | Present for 9,835 relations, including second-closest local lead | Structural chronology found |
| Instability outside a bound state | No active-state/free-bound comparison supplied by these records | Not established |
| Decay leaves K3 plus detached, fast residue | Real P4->paw transitions; fourth cores remain attached; no spatial clock | Not established |
| Extra orientation response in the specified aliases | Exact ordered vectors agree | Not observed by these probes |
| An identified neutron | Combined physical evidence absent | Not established |

Target proximity was used for selection and cannot also count as an independent
prediction. These are dependent saved histories and previously inspected discovery
data. No new tolerance, energy mapping or particle assignment was introduced.
A negative fixed-proxy or alias result does not establish neutron absence.

Independent checks cover all 14,304 exported exact ratios/errors, unique IDs and
support containment, all 59 P_trial attachments from raw parentage, all 523 temporal
cases and147 birth groups, all 28 component cuts, and both alias-vector groups.
All 36 original verification tests passed. Six prepared adjoining controls passed
12 legal repeat measurements; all 196 diamond/K4 anchor relations agree with saved
closed-family primary and mirrored responses. No resource censor occurred.

Source checkpoints, response arrays, old controls and both STOP markers are
unchanged. Full-service integrity remains 14/16 inheritance and 0/14 planned late
milestone certificates, with four planned 200M snapshots unavailable. These local
checks do not complete service/queue/RNG audits. No native worker remains running.
This finite pass is complete and preserved for review and resumption.

## Files

- [STRUCTURAL_LEADS.csv](STRUCTURAL_LEADS.csv): ten compact review examples,
  all exact local/fixed ratios, chronology and explicit missing physical criteria.
- [LEAD_CERTIFICATES.json](LEAD_CERTIFICATES.json): complete parent/witness/support
  certificates for those examples; no truncated parent cards.
- [Complete numerical shortlist](adjoining/SHORTLIST.csv):278 rows, including
  every tie in the prospectively specified local/fixed and support-change strata.
- Full attachment tables: [host0](adjoining/ALL_RELATIONS_HOST_0.csv),
  [host1](adjoining/ALL_RELATIONS_HOST_1.csv), [host2](adjoining/ALL_RELATIONS_HOST_2.csv),
  [host3](adjoining/ALL_RELATIONS_HOST_3.csv). Together they preserve all 14,304 rows.
- [All523 changed core identities](temporal/CHANGED_CORE_CASES.csv),
  [full transition matrix](temporal/INITIAL_TO_FINAL_TRANSITIONS.csv),
  [birth-event certificates](temporal/EVENT_CERTIFICATES.json),
  [temporal report](temporal/REPORT.md).
- [Alias report](aliases/REPORT.md), [five alias members](aliases/ALIAS_MEMBERS.csv),
  [exact alias diagnostics](aliases/DIAGNOSTICS.json).
- [Protocol](PROTOCOL.md), [source preservation manifest](SOURCE_MANIFEST.json),
  [independent validation](INDEPENDENT_VALIDATION.json). Each branch retains its
  own frozen method, checksums, chunk progress and completed receipts.

All files are plain text/JSON/CSV except retained local compressed certificates;
no new ZIP was created. Each plain file linked above is below 5 MiB. In adjoining
full tables, certificate_relation_index is host-global: the record offset is
certificate_relation_index minus chunk.start. The compact review table supplies
both offsets explicitly. Resolve source paths relative to DCU_NEUTRON_SEARCH_v1.
<!-- REVIEW_FILE_MANIFEST -->

Delivered review-file checksums (whole files); the review body checksum excludes
this manifest to avoid a self-containing full-file digest.

| File | Exact bytes | SHA256 |
|---|---:|---|
| `STRUCTURAL_LEADS.csv` | 7361 | `c2ccb5800d6bd9e09311b50b4f4eb7604d3eccba4dd0277198a3ed44d7f42bea` |
| `LEAD_CERTIFICATES.json` | 1225453 | `69c2157638e46ea4af84069b1cf7e64fb24a68636b62019700ec14fc2ec2c155` |
| `INDEPENDENT_VALIDATION.json` | 1941 | `1b5de74d8fb7baf1054d01502991b7f758486480d4746005e77955dd7ec73479` |
| `SUCCESSOR_REVIEW.md` | 18384 | Body (17721 bytes): `33c39fe3870d019d687bbb2863eca0074605a52c6665825685905a7cf02d7cde` |
