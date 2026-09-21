# Cell 11 — witnessed squared-distance mismatch: definition, proof and tests

## Scope

Input: the unchanged seven candidate cores, the thirty sparse controls, the
twenty before/after closure cases, and their S-overlap squared distances from
Cells 9–10. The physical premises remain a Euclidean three-dimensional local
realization and a preferred distance for each actually witnessed core edge.
No common global realization is asserted. The 20 closure cases overlap and
involve only two native added pairs.

The new diagnostic is
\[
 \eta_3^2 = \inf_{x_1,\ldots,x_n\in\mathbb R^3}
 \frac{\sum_{(i,j)\in E}(\|x_i-x_j\|^2-\delta_{ij})^2}
 {\sum_{(i,j)\in E}\delta_{ij}^2}.
\]
Here delta already denotes squared distance. Every witnessed edge has unit
weight in the numerator. Missing edges have no residual and no assigned zero.
The normalization is fixed by the target distances. Uniform scaling of all
squared distances leaves eta unchanged. Eta is not an average of individual
edge percentage errors and is not physical strain.

This is a selected mathematical diagnostic, not a native energy or recording
cost. Its numerator has fourth-power length scaling if a length valuation is
later imposed. Nothing in this test equates the optimizer's iterations to
native updates or maintenance ticks.

## Positive lower bound for complete five-point data

Let D be the source matrix of squared distances, J=I-11'/n, and
B=-JDJ/2. In the tested full five-point cases B is PSD and rank four.
Let ell be a positive rational lower bound on B restricted to 1-perp.

For any configuration X in R^3, center it and set G=XX'. There is a unit
vector v in 1-perp with Gv=0. Put
e_ij=||x_i-x_j||^2-delta_ij and SSE=sum_{i<j}e_ij^2.
Then
\[
 \ell\leq v^\top Bv
 = \sum_{i<j}e_{ij}v_i v_j
 \leq \sqrt{\mathrm{SSE}}
       \sqrt{\sum_{i<j}v_i^2v_j^2}.
\]
Cauchy–Schwarz also gives sum_i v_i^4 >= 1/n, so
\[
 \sum_{i<j}v_i^2v_j^2
 =\tfrac12(1-\sum_i v_i^4)\leq\frac{n-1}{2n}.
\]
Consequently
\[
 \mathrm{SSE}\geq\frac{2n}{n-1}\ell^2
 =\frac52\ell^2\quad(n=5).
\]
Dividing by the fixed target squared norm gives the certified eta-squared
lower bound used in the cell. It applies to every possible 3D realization,
not merely to the numerically found ones.

### Exact certification of ell

Let F have columns e_i-e_4, i=0,...,3. Then
F'F=I+11' and
\[
 A_{ij}=(\delta_{i4}+\delta_{j4}-\delta_{ij})/2= (F'BF)_{ij}.
\]
Thus A-ell*(I+11') positive definite is an exact sufficient certificate that
B dominates ell*I on the centered subspace. All its Schur-complement/LDL
pivots are checked with rational arithmetic. A 42-step rational bisection
brackets the smallest positive eigenvalue; its precision is numerical
bookkeeping, not a physical parameter.

## Feasible upper bounds, not proved optima

Each five-point case uses exactly 24 deterministic starting configurations:
four three-axis projections of the four-dimensional centered Gram
realization, plus the two Cell 10 tetrahedral completions for each of ten
released-edge choices. Releasing an edge generates a STARTING CONFIGURATION
only. All ten residuals are active throughout every least-squares run.

SciPy least_squares, method trf, uses the explicit analytic Jacobian and
ordinary quadratic loss. All positions are free; rigid translations and
rotations do not affect the objective. Initial and final configurations are
both retained when choosing the best feasible upper bound.

The best coordinates are rounded to 12 decimal places and interpreted as
exact rationals. The objective is then recomputed entirely with Fractions.
These rational 3D coordinates and their exact objective certify an UPPER
bound. Solver convergence or agreement of starts is never a proof of global
optimality. For K5-3, two distinct terminal objective values were encountered.

The squared-bound intervals in the JSON are exact rationals. The printed
square roots and percentages are approximate renderings.

## Zero controls

K3 and K4 zeros follow from the previously audited PSD Gram ranks <=3.
Each sparse control is rechecked against Cell 10's exact tetrahedral
certificate on every witnessed edge. The numerator is exactly zero, so its
minimum is exactly zero.

For all complete five-point data a four-dimensional coordinate realization
recovers every target distance. This is a dimensional countercontrol, not a
proposal to change the physical dimension to accommodate the candidates.

## Executed checks

- Two full executions agree in all saved result fields (excluding the callable) and printed output.
- All three K5 fits pass reverse-label and squared-distance x9 scale checks; exact lower bounds invariant, feasible upper errors <1e-13.
- Independently recomputed every rational 3D upper certificate (23 full five-point cases).
- Distinct rounded terminal objective values by selected K5: {'K5-1': 1, 'K5-2': 1, 'K5-3': 2}. No global-minimum claim from solver convergence.
- All 23 positive lower-bound matrices independently pass exact SymPy leading-principal-minor tests.

Primary execution: 4.174 seconds with NumPy 2.3.5 and
SciPy 1.17.0. The user cell requires NumPy and SciPy, not SymPy.
SymPy was used only in the additional independent exact certificate check.

23 positive full five-point problems x 24 starts = 552 solves.
No non-success termination was hidden or discarded (none occurred).
All old Cell 9/10 fields were unchanged.

## Interpretation

All tested full K5 metrics have a provably nonzero minimum under the chosen
diagnostic, while all 30 actual sparse controls have zero. K5-1 and K5-2 share
the same target geometry up to relabeling, and their bounds agree. K5-3 has a
larger bounded mismatch. This is embedding-dependent, not a universal
constant of graph order.

Nothing in these results identifies eta, eta squared, or their numerator
with a physical energy, force, native workload surcharge, or mass. In
particular, the K3 reference has eta=0, so it cannot be used to calibrate a
nonzero electron mass from this diagnostic alone. The pre-existing physical
calibration is untouched.
