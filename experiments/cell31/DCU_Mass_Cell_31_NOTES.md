# Cell 31 — continuing carriers, shared endpoints, and observer identity

## Result and status

The same native histories now support three specified observer readouts: fixed
anchors, following a uniformly chosen incident newborn (`follow_any`, primary),
and following only a unique incident newborn (`follow_unique`, control). Actual
parent-to-child carrier chains exist, but this finite panel does not demonstrate
more reliable communication. Following introduces shared-carrier events, and a
native ancestry certificate alone does not authenticate the private choice of
continuing branch. A constructive example has the same native events and the same
source-local record in two cases with different tracked-receiver acquisitions.

The source constructor, global sampler, all enabled births, first/repeat costs,
relief, and initial recording backlog are unchanged. The observer rule is passive;
its memory and sibling-choice machinery are NOT compiled into extra native work.
No spatial trajectory, proper time, SI scale, particle, quantum-state transport,
physical decay, or relativistic redshift is identified. The neutron search is not
read, patched, or controlled by this experiment.

## Source basis

- *Self-Limited Growth, Intrinsic Time, and Record Geometry in the Distinction
  Combinatorial Universe*, revised manuscript 0.5 (September 2026), sections
  2.1–2.5 (native objects, service, complete bursts, costs and clocks), 4.3–4.5
  (reader frames, locality and measurement back-action), and 6.3 (regional labels).
- Original `DCUStructure` and original Cell25 `step`, included unchanged in
  `reference`. The loader extracts definitions without executing their old notebooks.
- Cell29's conditioned emission checkpoints; Cell30's actual input convention and
  full fixed-anchor results. They are reference data, not a newly generated genesis.
- No RMR mass expression, empirical physical input or unrelated formula is evaluated.

The reader law and mathematical deductions below are this experiment's proposals
and consequences, not claims that the source already defines these observers.

## 1. Frozen observer experiment

All policies start with the SAME old source anchor and receiver anchor at the SAME
post-emission checkpoint. Source knows the existing pulse child through its actual
emission; it does not receive a free initial carrier jump onto that pulse. Its old
recording obligations remain. Receivers are 3 (primary) and 4 (transfer).

At a native service:
1. Count service only to the carrier that existed BEFORE the draw.
2. Read all incident newborn children and their parent cards, as in Cell30.
3. After the whole batch, optionally advance to one incident newborn.
4. Retain the accumulated clock and the union of previously obtained records.

`fixed` never advances. `follow_any` advances whenever at least one child is
available and selects uniformly among two choices. Its private generator is
independent of the native RNG. It reads neither a future event, pulse tag, target,
workload nor sector to make that choice. `follow_unique` advances only if exactly
one child is incident, remaining at its old carrier at ambiguous bursts. This is
a DIFFERENT rule/control, not an alternative selected after seeing its results.

There is one active carrier per reader at a time, not a simultaneous descendant
cloud. Unchosen children remain real native objects. All readers shadow one native
history; they do not introduce multiple real apparatuses with unaccounted capacity.

A pair protocol stops at its first shared active carrier or at the fixed horizon.
The merge event is observed and reported separately. The NATIVE HISTORY continues
to the horizon regardless. No synthetic separation or extra ownership label is
inserted after merger. This is a competing endpoint of the distinct-reader protocol,
not proof that physical observers necessarily cease to exist when they interact.

Every reader has a private log with local ticks, active-carrier handles, incident
newborns, encountered partners and newly acquired parent cards. No global iteration,
queue, total population, seed or unrelated births appear inside that log. Full
results also contain auditor fields; passing the entire dictionary to a putative
internal observer would violate this access convention.

## 2. Receipt and return conventions

A receipt is an actual incident birth r={b_t,c} of the current receiver b_t whose
parental ancestry contains the pulse m. Every valid receipt root is retained, not
only a remotely inaccessible first receipt.

A returned certificate must contain an OLD receipt in a later source-incident child.
A newborn receipt cannot be returned through a sibling in the same batch. Receipts
jointly witnessed by the source at their birth are separately marked: repeating
an already shared receipt later is not counted as a new remote acknowledgement.
The primary return statistic therefore uses previously nonjoint (remote) receipts.
This distinction was implemented before the primary run; it does not affect the
Cell30 fixed-anchor regressions, whose receipts are all nonjoint.

Receipts and returns occurring on the merging event are reported as such, not
silently discarded and not advertised as continued communication between separated
endpoints. Return timing is in the source's retained tick count. A successful-only
summary is not an unconditional mean travel time. Unfinished outcomes stay saved.

For moving readers, validating that the receipt was issued by the actual active
receiver requires its continuation track. The auditor has this log. Parentage alone
does not necessarily supply it to the source; section 6 gives an exact counterexample.

## 3. Exact adaptive-clock result

Let a_t,b_t be the current carriers selected by past information before service,
and let R=F+n, p=q/R, p2=q(q-1)/(R(R-1)). As long as the carriers are distinct,

    E[dN_a | past] = E[dN_b | past] = p,
    E[d(N_a-N_b) | past] = 0,
    E[d(N_a-N_b)^2 | past] = 2(p-p2).

These follow from uniform token sampling even though the labels vary adaptively.
If S=min(horizon, first merge), then bounded stopping and orthogonal martingale
increments give

    E[N_a(S)-N_b(S)] = 0,
    E[(N_a(S)-N_b(S))^2] = 2 E[sum_{t<=S}(p_t-p2_t)].

The merging event still has distinct pre-service carriers and is included. These
are ensemble identities, not equal tick counts on every realization. Conditioning
only on successful returns or surviving distinct pairs would introduce a different
selection and does not preserve this statement automatically.

Changing the observer's carrier therefore does not by itself create a systematic
relative mean maintenance-clock rate. A nontrivial time dilation law would need
more than a predictable choice of one globally exchangeable token.

## 4. Exact continuation and merger laws

Let d be the current carrier's co-parent degree: the number of old vertices with
which its pair already exists. Conditional on this carrier's tick, among the other
q-1 selections there are n-1-d absent-pair maintenance partners and F+d neutral
requests. The number J of incident newborns obeys

    P(J=j | own token served, state)
      = C(n-1-d,j) C(F+d,q-1-j) / C(F+n-1,q-1).

Thus follow_any has hop probability 1-P(J=0), while follow_unique has P(J=1).
The code enumerates 75 finite domains to check this law independently.

After the first follow_any hop, the active carrier is always on the childless
(unrecorded) shell before its next service. It was newborn when adopted and must
be served to get any child; as soon as that happens the observer follows a new child.
This is a consequence of the readout rule, not a matter classification. Old records
and unchosen siblings remain in the universe and in the observer's declared memory.

For two distinct childless carriers at Gamma=3, conditional on their co-service,
there is one remaining sampled request. A forced request gives one shared child and
certain merging when both follow. A third maintenance token gives two choices per
reader, with one shared choice among four combinations. Hence

    P(merge at next update | state)
      = p2 * [F + (n-2)/4] / (R-2)       [follow_any],
      = p2 * F/(R-2)                     [follow_unique at childless carriers].

These special-case laws are checked by enumeration in 18 finite domains. They
explain why merging is a structural possibility, not floating point duplication.
After merger we make no assumption about whether a physical observer should merge,
bifurcate or retain distinct internal states; the simple one-carrier representation
has reached its stated boundary.

While a reader is uninformed, the old first-pulse contact law remains

    1 - C(R-1-M,q-1)/C(R-1,q-1)

per reader tick, with M current pulse carriers. Before first acquisition, its
current carrier cannot already have a formed pair with a pulse carrier without
having observed that pair. The law contains carrier counts, not spatial distance.
Equal conditional hazards do not imply identical unconditional waiting-time laws
under different adaptive selection/stopping histories.

## 5. Execution domain and retained results

Replayed cohort: all 256 original Cell30 seeds (20300920..20301175) at four emission
states, 4096 updates/history. Fresh cohort: 128 seeds (20320920..20321047) at each
of the SAME states, with the same horizon and policies. Both blocks were fixed
before the primary run, and there was no outcome-dependent extension. The historical
replay is a new observer analysis of previously seen native histories, not independent
new physical data. The fresh cohort is new stochastic continuation, not laboratory data.

There are 1024 replayed native histories and 512 fresh histories. The total evaluated
native updates are 6,291,456, of which 4,194,304 replay Cell30 and 2,097,152 are fresh.
Each history supports six pair protocols (three rules times two receivers), sharing
one native graph and event history. Repeated labels/returns across those protocols
are not independent observations or distinct native objects.

Aggregated, DESCRIPTIVE protocol counts:

| Cohort | Rule | Pair protocols | Receipts before / at merge | Remote returns before / at merge | Merges |
|---|---|---:|---:|---:|---:|
| Replayed | fixed | 2048 | 201 / 0 | 12 / 0 | 0 |
| Replayed | follow_any | 2048 | 200 / 3 | 2 / 5 | 95 |
| Replayed | follow_unique | 2048 | 202 / 2 | 2 / 4 | 92 |
| Fresh | fixed | 1024 | 103 / 0 | 0 / 0 | 0 |
| Fresh | follow_any | 1024 | 98 / 2 | 0 / 1 | 44 |
| Fresh | follow_unique | 1024 | 96 / 2 | 1 / 0 | 45 |

87 of the primary follow_any mergers, and 41 of its fresh mergers, occur before
any receiver acquisition. Merger is not simply being relabeled a successful message.
The complete per-state/per-receiver comparisons and paired Monte Carlo errors are
in the output/results. There is no joint significance obtained by treating the two
receivers, shared source paths, policy shadows or repeated states as independent.

The maximum observed carrier sequence has seven hops, but the mean numbers of hops
are much lower (about 1–2 per reader in most panels). A hop is one native ancestry
edge, not a metre. The policy enforces at most one hop per own tick; that inequality
is not a derived physical speed limit.

No clear communication improvement is established. The primary single-carrier
rule has almost the same receipt count in the replay cohort and fewer completed
returns. The fresh block preserves the lack of a clear improvement and rarity of
returns. This does not establish that all possible mobile/continuing observer
models fail; it evaluates these specified choices and endpoints.

Example of a genuine pre-merger tracked return: replay cohort H6/first-use,
seed20301043. Source track 25->40->64->66; receiver track 3->43->57. Receipt43
predates return66={43,64}; at return the current receiver is57 and source moves
64->66, so their carriers remain distinct. The source counts11 own ticks; the
receiver has7. Full burst work at return is2332. These are not seconds, energy,
relative motion or gravitational clock comparison.

## 6. Constructive identity ambiguity

Use the actual H0/recorded emission state: source13, pulse18, receiver3. A legal
conditioned maintenance burst on 3,0,2 creates

    u19={0,3}, v20={2,3}.

The other pair {0,2} already exists; it is not recreated. Follow_any permits either
child with probability1/2. Next condition actual services that create

    r21={18,20}, then w22={13,21}.

Keep these native events IDENTICAL in two copies of the reader description.

- If the receiver chose u19, it never participates in r21. The apparent receipt
  was made by an unchosen descendant of its old anchor.
- If the receiver chose v20, r21 is its actual acquisition; it advances to21.
  The return at22 merges its carrier with the source under the simple rule.

The source's complete private local log, parent cards, own tick and final carrier
are identical in both descriptions. The native final graph, queue and all recording
charges are also identical. The only difference is the private continuation choice.

Consequently the source can recognize a descendant-origin ancestry certificate,
but it cannot infer from those cards alone which positive-probability private
continuation track was actually followed. This is an information-loss counterexample
for the present passive observer representation. It does not assert a universal
impossibility of observer identity in the DCU. Recording the choice would be extra
operational content that should be accounted for, not silently supplied by the auditor.

The counterexample and merger enumeration were appended as explanatory checks after
the first full simulation; they did not alter the frozen cohorts, policies, counts,
clock estimates or native histories. All scientific cohort data reproduced unchanged
in the subsequent complete execution of the final cell.

## 7. Reproduction and validation

The final cell ran with the original source functions, then all private logs and
certificates were checked independently. The delivered validator covers:

- 253,186 private event cards and18,432 reader frames across9,216 pair protocols;
- 1,295 actual active-receiver receipt cards,27 protocol-level returned certificates,
  and276 coalescence records (including repeated policies/receiver tests);
- 65,536 full native transitions using an independent bit-mask ancestry and sequential
  first/repeat charge implementation, across16 observer-free histories;
- all48 per-cohort/state/policy/receiver summary rows;
- 75 exact incident-child distributions,18 exact merger distributions;
- three actual conditioned events in the identical-source-view counterexample.

The validation is about implementation and stated finite mathematics. Simulated
fluctuations are not pass/fail filters; some clock residuals exceed two empirical
Monte Carlo standard errors and remain visible. They do not overturn the exact
conditional law or justify selecting additional seed blocks.

No prior dictionary, neutron search, source paper or old notebook result is modified.
The standalone folder retains source files and exact inputs. `reproduce.py` writes
REPLAY_* files rather than overwriting the supplied outcomes. Native integer data
and event decisions are compared exactly; negligible platform-dependent floating
summation differences are reported separately where applicable.

## Next well-posed step

Test a native, explicitly charged record of observer continuation, rather than
assuming a selected pointer is already a physical identity. A continuation needs
an actual record linking an observer's retained state to its next carrier. Then ask
whether a peer can verify that linkage using only received cards, and measure the
additional native opportunities and work. Do not grant free exclusive ownership,
new observer service capacity, a routing metric, or a spacetime interpretation to
make the test pass.

### Completed standalone replay
The final delivered cell was executed from a separate directory containing spaces.
Every serialized result reproduced, with zero floating differences (137.5 seconds
in this Linux runtime). This verifies independent-folder reconstruction with the
bundled sources and checkpoints. No Mac execution has been performed here.
