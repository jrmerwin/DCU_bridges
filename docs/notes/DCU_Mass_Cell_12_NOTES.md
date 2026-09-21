# Cell 12: response-law scale and the local self-stress connection

## Inputs and scope
The saved Cell 5 response summaries and Cell 11 rational geometric certificates
were read directly. The inputs were not regenerated or changed. This is a new
calculation from those saved results, not a new native-DCU interaction history.
The provided cell runs after those two cells in the user's notebook.

The coefficient diagnostic asks whether q_eff(M) = Q(M)/Q_e + g eta_3(M)^2
could give each chosen K5 the same muon/electron target. Q_e remains 248/9.
Here g = gamma/Q_e if the additive work law were Delta Q = gamma eta_3^2.
No g was adopted or applied as a model parameter. The inverse intervals are
diagnostics of what a fit would require, not fitted predictions.

The target is the 2022 CODATA listed central muon/electron ratio 206.7682827
(standard uncertainty 0.0000046), from
https://physics.nist.gov/cuu/Constants/Table/allascii.txt ,
checked 2026-09-19. The code reuses the same saved target. This uncertainty is
not propagated into the printed g intervals; those intervals account for the
certified geometric bounds only. It is negligible for the nonintersection
between K5-1/2 and K5-3 but they must not be called experimental confidence
intervals or an exact physical measurement.

## Coefficient identifiability
The electron candidate has eta_3=0. Its anchor therefore fixes no g at all.
For L<=eta_3^2<=U, q<target and g>=0, a necessary interval for matching the
target central value is [(target-q)/U, (target-q)/L]. Disjoint intervals
exclude a common coefficient for the proposed same-species assignments under
this ansatz even without knowing the exact global minimum. Different native
state/species assignments would be a different test, not prohibited by this
diagnostic.

## A common geometric functional and its derivatives
Use exactly the existing normalized loss, with a conventional factor 1/2:
Phi(x) = r(x)'r(x)/(2 C),
r_e=||x_i-x_j||^2-delta_e, C=sum_e delta_e^2.

R is the Jacobian of r. Direct differentiation gives
gradient Phi = R' r/C,
Hessian Phi = (R'R + 2 Omega(r) tensor I_3)/C,
where Omega(r) is the graph Laplacian with signed edge weights r_e.

At an equilibrium, R' r=0; the edge residual vector then lies in the
self-stress space ker(R'). This is an algebraic feature of the hypothetical
quadratic interaction, not a newly discovered native force. In particular,
a mere error norm eta does not determine the Hessian or the direction of its
self-stress.

For complete nondegenerate core frameworks the observed Jacobian ranks are
3, 6, 9 for K3, K4, K5; self-stress-space dimensions are 0, 0, 1.
Saved K5 residuals project into that one-dimensional space with fractional
remainders below 4e-8. These are numerical stationarity checks, not exact
equilibrium proofs. The exact mismatch brackets remain those of Cell 11.

The local Hessian is projected orthogonally away from all six rigid motions.
Every tested projected spectrum is positive to numerical precision.
This is a local curvature/stability observation only: it does not prove global
minimality. K4-1 has a softer nonrigid eigenvalue than the three K5s.

C differs among candidates: this normalized functional is NOT the energy
of identical, absolute-stiffness springs across every candidate. To hypothesize
a common local edge energy one must specify its units and use the corresponding
un-normalized energy explicitly. Do not mix its coefficient with the coefficient
of the normalized diagnostic.

No kinetic matrix, clock, damping, external forcing, action/energy conversion,
native repair schedule, or geometry-dependent recording surcharge is supplied.
Hessian eigenvalues are not masses or frequencies. With a supplied kinetic
matrix M, mode frequencies would require a generalized problem H v=omega^2 M v.

## Relationship to the RMR matrix
The paper's unweighted T_n=-A(K_n)/(n-1) does not depend on x, delta or eta.
Its eigenvalues and f_n consequently do not change under a metric deformation
with the same clique graph. A geometry-dependent replacement would be an
additional operator. The number 144 is not used here.

## Executed checks
- Exact rational coefficient intervals recomputed from the saved fractions.
- All recorded geometric residual magnitudes agree with the saved feasible
  upper certificates within 1e-11 after floating-point conversion.
- Analytic Hessian-vector products agree with independently differenced gradients
  at all seven saved configurations (relative check threshold 1e-7).
- Full SVD removes exactly six rigid motions for each candidate.
- Saved summaries and geometric candidate records compare unchanged after use.
- Two executions give identical output and saved numerical results in this runtime.

Execution uses NumPy only; no nonlinear optimizer is called.
