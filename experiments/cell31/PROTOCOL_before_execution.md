# Cell 31 protocol fixed before the primary run

Question: Does following newly constructed descendants improve ancestry contact and
return, without adding a service budget, and when do two observer tracks cease to
have distinct carriers?

Source inputs: Cell30's four emission states, Gamma=3,m=0,H=0 or6, complete native
bursts and Cell25's unchanged native step. Primary replay: all 256 Cell30 seeds at
each state, 4096 updates. Confirmation: 128 new seeds 20320920..20321047 at each
state, same horizon. Both receivers 3/4 on the same histories. No optional sample
extension after seeing a statistic.

Observer policies (passive annotations of the SAME native histories):
- fixed: never replace the anchor (Cell30 regression).
- follow_any (primary): on a productive service of the current carrier, uniformly
  choose one of its incident newborn children; use an independent observer RNG.
  At a singleton no random draw is needed. No scoring by pulse, target, queue,
  costs, sector, or future behavior. One active carrier only.
- follow_unique (control): advance exactly when there is one incident newborn;
  retain current carrier at ambiguous two-child bursts. This avoids a tie selector
  but has a different opportunity law; it is NOT secretly pooled with primary.

All incident children are readable regardless of which is followed. Memory is the
union of their ancestries; local clock counts the pre-service current carrier.
Transition occurs once after the whole burst, no reset of ticks, no newborn service
or multihop within a burst. Source starts at its old emission-parent anchor in ALL
policies; it does not get a free initial jump onto the pulse.

Each source/receiver pair is a counterfactual observer protocol, not an army of
cost-free physical observers. Stop its protocol at first common active carrier
(coalescence) or the fixed horizon; the NATIVE PROCESS never stops for this reason.
The merge event itself is recorded. Receipts/returns on that event are separated
from receipts/returns strictly before merger. Do not classify merger as faster
spatial travel. All censored/merged cases remain; no success-conditioned mean is
called an unconditional travel time.

A receipt is any incident child of the actual current receiver whose OLD parental
ancestry contains the pulse. A return is an incident child of the actual current
source containing a PREVIOUS receipt in its ancestry. Newborn siblings cannot
return a newly born receipt. All receipt roots can certify return, not only first.
Auditor validation of moving identity uses the actual retained continuation track.
Native ancestry alone need NOT authenticate which branch was followed. Demonstrate
this distinction with a legal fork counterexample. Do not advertise an
independently readable receiver-identity certificate when only an auditor supplied
that identity. Fixed-anchor certificates retain their original scope.

Analytical targets: adaptive single-carrier tick probability q/(F+n); stopped
clock-difference variance compensator 2*sum(p-p2) while carriers distinct; exact
conditional number-of-incident-newborns law using current co-parent degree.
Movement/receipt choices never feed back into native RNG or workload. Plain logs
exclude global steps/pool/seed and use local opaque handles. Auditor data separate.

No quantum-state transport, masses, SI units, gravitational redshift, speed limit,
spatial dimension, or modified dynamics is to be inferred from carrier following.
