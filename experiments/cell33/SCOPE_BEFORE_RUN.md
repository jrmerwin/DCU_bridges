# Cell 33: fixed midpoint echo and two-window memory

Retrospective analysis of all six Cell29 ensembles, 512 histories each, through
1024 service iterations. No new native evolution. Primary pulse at 512, final
basis restoration at 1024. Original phase depths 3 and 4. No timing scan, fit,
postselection, or new load/energy rule. No script accesses the Mac search.

Proposed coherent control: J|r>=|2-r> on BOTH qutrits at midpoint and endpoint.
J = X^2 S X for the cyclic permutation X and entry reflection S=(0 1).
This matrix can be expressed in the declared finite representation; its active
implementation, schedule, memory, and physical gate cost are NOT derived.
The 512/1024 schedule is an external service-iteration schedule, not a compiled
internal clock. Sequence is J D(n_late/q) J D(n_early/q) on each factor.

Compare uncorrected free evolution with echo and retained full-count inverse.
Use the same Omega2 and heralded P01 F3 P01 state, same source Bell settings.
Measure fidelity to each initial state, negativity, Bell score, and phase moments.
Original free evolution must replay Cell32. Means and covariance of equal-duration
window counts describe the sample; counts across windows are not presumed iid.

Memory control: exact product of empirical early and late JOINT two-factor
count distributions; no resetting of real native histories. Retain cross-factor
correlations within each window. Distinguish this constructed channel from a
realizable native intervention or proof of quantum non-Markovianity.
Two native seed halves remain visible; no outcome-triggered seed expansion.

Analytic controls: equal counts in early/late windows refocus exactly; zero
counts do nothing. Two finite noise laws with identical midpoint and endpoint
free channels but different temporal pairings have distinct echo responses.
No claim that every such algebraic noise law is realized by the native process.
