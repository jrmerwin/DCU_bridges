# Cell 33 — a fixed echo, native timing correlations, and limits of recovery

## Status

A fixed two-pulse permutation sequence substantially increases overlap with the
initial source. At the primary 1/27 phase increment it does NOT recover a violation
of the fixed CGLMP-3 witness at the terminal checkpoint. At the already retained
1/81 control increment it does recover a violation in all six saved ensembles.
The entanglement measure does not improve uniformly: phase realignment and
recovery of entanglement are distinct.

The native service process, recording work, six initial states and 3,072 histories
are unchanged. There are no new stochastic evolutions, native objects, mass
formulas, SI valuations, regional service pools, or neutron-search operations.
This is a retrospective, conditional quantum-instrument calculation, not an
experimental demonstration or a uniquely derived autonomous quantum controller.

## 1. Source basis

- *Self-Limited Growth, Intrinsic Time, and Record Geometry in the Distinction
  Combinatorial Universe*, revised manuscript 0.5, September 2026: Sections
  2.3–2.5 (uniform sampling, complete construction bursts, charges, clock);
  5.1/5.4 (finite cyclic/reflection representation and phase clock); 7.1 and
  Appendix B.3 (shared source, filtered preparation, signed Bell functional).
- Cell29: actual saved service counts at 0,1,8,32,128,512,1024 in all six arms,
  including the conditioned preparation backlog. Every one of 512 histories/arm
  is included. Its source code and original Cell25 step are supplied for validation.
- Cell32: full joint-count quantum channel, source Omega2, filtered source,
  fixed settings (0,4/27,-2/27,2/27), and phase increments 1/27 and 1/81.
- External motivation ONLY: E. L. Hahn, *Spin Echoes*, Physical Review 80,
  580–594 (1950), doi:10.1103/PhysRev.80.580. Pulse refocusing is a standard
  control idea. No magnetic field, nuclear spin, experimental relaxation constant,
  or phenomenological damping function is imported into this model calculation.

The RMR mass/screening prescriptions are not used to select a pulse or a result.
In particular, the screening paper's signed mass displacements are distinct from
this phase channel, which preserves diagonal populations at the final readout.

## 2. Fixed protocol and added assumptions

Let D(x)=diag(1,exp(2 pi i x),exp(4 pi i x)). The primary noise law remains
D(N3/27) tensor D(N4/27), with the existing denominator-81 control.
Define J|r>=|2-r>, which exchanges basis labels 0 and 2 and fixes 1.
This is NOT the source entry reflection S=(0 1).
For X|r>=|r+1 mod 3>, J=X^2 S X, so the matrix is expressible through the
inherited finite permutation representation. Representation compatibility is not
proof of physical gate availability. ACTIVE application of J is a new explicit
ideal control hypothesis.

Apply J to both factors immediately after iteration 512 and again immediately
after iteration 1024, before final readout. The final pulse restores the original
basis, important for the filtered preparation. No intermediate projective
measurement, history-based adjustment, timing scan or selective acceptance is used.
The same schedule applies to every arm, depth, history and input state.

Scheduling uses the external native-iteration index, NOT seconds or a demonstrated
internal reference clock. The pulses have no assigned native work or duration.
As in Cell32, the attached instrument does not feed back into the counting law.
Therefore the saved native histories remain valid counterfactual inputs under this
specified ideal control. A costed controller could change them and requires a new test.

## 3. Exact echo identity

For local first-window and second-window counts a and b, respectively:

    J D(b/Q) J D(a/Q) = exp(4 pi i b/Q) D((a-b)/Q).

The scalar factor has no effect on the density matrix. For two factors the final
channel is the mixture of product phase gates with local counts

    e3 = N3(512) - [N3(1024)-N3(512)],
    e4 = N4(512) - [N4(1024)-N4(512)].

The direct source sees their sum. The filtered preparation requires their joint
histogram; it cannot be propagated correctly with the summed count alone.

Exact recovery occurs when each factor accumulates equal counts in the two
windows, modulo Q. Equal WINDOW LENGTH does not imply equal SERVICE COUNTS.
Changing mean service rates and timing fluctuations can both prevent recovery.

For Omega2, put A=N3(512)+N4(512) and B=the corresponding second-window total.
The free and echo moments are E exp(i*omega*(A+B)) and
E exp(i*omega*(A-B)), respectively. No Gaussian approximation is used.
The identity

    Var(A+B)-Var(A-B) = 4 Cov(A,B)

is checked exactly with rational empirical moments. Positive correlation makes
the echo phase difference narrower in variance; negative correlation has the
opposite effect. Variance alone need not determine periodic coherence outside a
small-phase approximation. There is no assumed universal sign for long windows.

## 4. Readouts and uncertainty

Both preparations, phase depths and fixed Bell settings are unchanged from Cell32.
Fidelity is <psi_initial|rho|psi_initial>, with the pure-state convention (not its
square root). Negativity is the sum of negative partial-transpose eigenvalue
magnitudes. For Omega2:

    negativity = (2 |C1| + |C2|)/3,
    fidelity = (3 + 4 Re(C1) + 2 Re(C2))/9.

Fidelity and the fixed Bell witness depend on phase location. Negativity is
invariant under a known local phase rotation. An echo can align the source to its
measurement basis while leaving less entanglement than free evolution.

All means concern EXPECTED quantum outcomes under simulated histories, not
measurement-shot data. Fidelity/Bell errors are history-level standard errors.
Negativity errors are jackknife diagnostics for a nonlinear finite-mixture
estimate. The two original 256-seed blocks are separately retained. States and
policies share histories; their values are not independent physical confirmations.
No significance result is obtained by combining these dependent tests.

## 5. Temporal-information control

The independent-window surrogate is the product of the empirical early and late
JOINT two-factor count laws. It keeps within-window N3/N4 dependence but discards
cross-window dependence. Its channel is computed over ALL 512^2 early/late
pairings by factorization, with an independent full histogram check.
It is not a rerun with an actual native environment reset, and it need not be an
admissible history ensemble of the constructor. Its role is information loss.

We retain the differences between actual and product characteristic moments and
their jackknife errors. Small differences in finite empirical mixtures are not
alone promoted into proof of quantum non-Markovianity.

## 6. Exact native two-step dependence

At each of the six starting states, enumerate every maintenance subset W of size
at most three, weighted by binomial(F,3-|W|). Compute the complete enabled burst,
first/repeat workload, relief, and next pool R'. No partial burst is substituted.
With J1=number of tagged tokens 3,4 in W, uniform service at the next state gives

    E[J2 | next native state] = 6/R'.

Weighted exact rational summation gives E[J2|J1=j] and Cov(J1,J2). All six
covariances are strictly negative, checked with exact rational arithmetic. This
is a local finite-state result, not a universal sign theorem. It establishes that
these successive unconditioned tagged increments are not independent; their
dependence is mediated by the changing native state.

The full native process remains a specified classical Markov process. This
calculation makes NO claim about CP-divisibility, quantum environment memory,
information backflow, a time-homogeneous reduced master equation, or a fitted
Lindblad law. Those are distinct questions.

## 7. Exact same-channel/different-echo counterexample

On one factor let Q=27. In model A the early count k is uniform on 0..26, and the
late count equals k. In model B both counts are independent uniform variables.
At the midpoint both give the same uniform phase channel. At the endpoint their
free phases also have the same distribution: 2k is uniform modulo 27 in A, and
a+b is uniform in B. Thus both free channels agree for every input state, not
merely for one observed probability.

Echo leaves phase zero in A and a uniform phase in B. Applied to Omega2, one
output is Omega2 with negativity 1; the other is (1/3)sum |rr><rr|, with negativity
0. Their trace distance is 2/3. Enumeration uses 729 equally weighted atoms per
model. These algebraic models are NOT claimed to be native DCU history laws.

This shows why two separately known free-evolution channels are insufficient to
predict a controlled two-time experiment. One needs their temporal relationship.
The saved paired histories provide that relationship in our DCU calculation.

## 8. Energy and screening boundary

The free channel and the complete echo cycle are diagonal phase channels. Both
preserve every diagonal density-matrix entry for any input. Therefore for ANY
energy assignment diagonal in this basis, its mean energy is unchanged between
initial and final readout. No particular numerical energy scale is supplied.

The intermediate J pulse may change populations for a generic input and might
require energy/work. Its cost is not modeled. The final identity is not a proof
of thermodynamically free control.

This excludes interpreting THIS phase channel as a calculated mass/self-energy
shift or as the screening paper's thermodynamic correction. It does not exclude
an additional interaction or energy operator producing such a correction.

## 9. Reproduction

Notebook: paste the single Cell33 file after Cell32. Only the saved Cell29 state
is mandatory; an available Cell32 result supplies additional terminal regressions.
The input dictionaries are not modified. There are no network calls.

Standalone: Python 3.10+ and NumPy. Run `python reproduce.py` in this folder.
It replays Cell32 from the frozen full Cell29 histories and then Cell33. It checks
exact fields exactly and floats with a small numerical tolerance to accommodate
platform/BLAS roundoff. It writes REPLAY files, leaving reference data unchanged.
The validation script additionally requires SciPy and mpmath. Linux was tested;
no Mac runtime or physical precision is claimed.

Primary scope was recorded before first execution of the echo comparisons.
The exact two-step enumeration and same-channel counterexample were added as
explanatory audits; no pulse time, mass parameter, source, seed, or successful
subcohort was selected in response to their results. The whole exercise uses
historical data and is not described as a preregistered experimental prediction.

## 10. Complete terminal summary

### Phase increment 1/27

|Environment / arm|Fidelity free|Fidelity echo|Negativity free|Negativity echo|Bell free|Bell echo|
|---|---:|---:|---:|---:|---:|---:|
|no_relief / existing_pair|0.100638|0.550260|0.439443|0.394826|-0.743937|0.784958|
|no_relief / recorded_pair|0.083722|0.628918|0.493246|0.455146|-0.820842|1.113608|
|no_relief / first_use_pair|0.096407|0.645962|0.477603|0.495090|-0.826192|1.129475|
|relief_H6 / existing_pair|0.097689|0.672972|0.564899|0.609250|-0.850216|1.341728|
|relief_H6 / recorded_pair|0.106595|0.672854|0.583520|0.583980|-0.815555|1.326220|
|relief_H6 / first_use_pair|0.106350|0.706183|0.590002|0.619059|-0.842392|1.475126|

### Phase increment 1/81

|Environment / arm|Fidelity free|Fidelity echo|Negativity free|Negativity echo|Bell free|Bell echo|
|---|---:|---:|---:|---:|---:|---:|
|no_relief / existing_pair|0.423143|0.905648|0.893182|0.878643|0.223300|2.412713|
|no_relief / recorded_pair|0.532108|0.928112|0.910168|0.894752|0.731463|2.495852|
|no_relief / first_use_pair|0.569554|0.936354|0.901191|0.911087|0.909396|2.492837|
|relief_H6 / existing_pair|0.638631|0.943381|0.930172|0.937579|1.235113|2.606025|
|relief_H6 / recorded_pair|0.652061|0.945093|0.934095|0.934782|1.298462|2.607818|
|relief_H6 / first_use_pair|0.665517|0.952652|0.937142|0.939987|1.362223|2.637287|

### Filtered source: fraction of initial negativity, 1/27

|Environment / arm|Free|Echo|
|---|---:|---:|
|no_relief / existing_pair|0.602148|0.574275|
|no_relief / recorded_pair|0.643734|0.628047|
|no_relief / first_use_pair|0.639273|0.663298|
|relief_H6 / existing_pair|0.701037|0.723141|
|relief_H6 / recorded_pair|0.728478|0.725003|
|relief_H6 / first_use_pair|0.719658|0.734942|
