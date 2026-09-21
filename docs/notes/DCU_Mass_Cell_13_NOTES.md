# Cell 13: saved Hessians, prestress, and a conditional service-coupling diagnostic

## What was and was not run

Inputs are the existing Cell 11 coordinates/distance certificates and Cell 12
full spectral results. No candidate selection, metric, normalization, or
coordinate fit was changed. No optimizer or native service/construction history
was run. The new computations are a decomposition of the saved Hessian and a
per-edge-gradient diagnostic.

The quadratic response functional remains a proposed interpretation:
Phi(x) = sum_e r_e(x)^2 / (2 C), where
r_e = ||x_i-x_j||^2-delta_e and C=sum_e delta_e^2.
Delta is a squared-distance target. C is fixed for each candidate but differs
between candidates. These normalized potentials do not represent an assumed
common absolute edge stiffness.

## Spectral decomposition

H_material=R'R/C; H_prestress=sum_e r_e Hessian(r_e)/C.
H=H_material+H_prestress.
The same six rigid-motion directions as Cell 12 are removed before taking
spectra. All seven saved spectra are reproduced.

Subtracting H_prestress is an algebraic comparison at the SAME coordinates.
It is not a native deletion or a new relaxation to a stress-free state. The
minimum-eigenvalue change compares the two minima and permits their eigenvectors
to differ.

For K5-1 and K5-2, minimum material curvature is about 0.0399420099 and
minimum full curvature is about 0.0360429150 (9.7618896% reduction).
For K5-3, the corresponding values are about 0.0399696454 and 0.0200435310
(49.8531178% reduction). K4-1 still has a smaller full curvature,
0.0167165508. Positive spectra are floating-point local curvature observations
at nearly stationary saved fits, not exact equilibrium or global minimum proofs.

The additional C*H output checks the effect of removing the per-candidate
normalization; it does NOT adopt a new physical potential. K4-1 remains softer
than the K5s under this check. Neither spectrum is a frequency spectrum without
a kinetic or relaxation law.

## Conditional update proposal and service moments

Define g_e=gradient(r_e^2/(2C))=r_e R_e'/C.
At an exact stationary configuration sum_e g_e=0. Individual g_e may be
nonzero for a frustrated configuration.

A POSSIBLE ADDED RULE would let service of the existing witness for e change
attached coordinates by -epsilon*g_e, with one common mobility epsilon.
Native service does not already perform this coordinate update. The analysis
neither chooses epsilon nor identifies coordinate changes with native
recording events.

Let N=F+n be the number of pooled requests, q=min(Gamma,N), and I_e indicate
selection of witness e's maintenance token in the native uniform sample of q
requests without replacement. All other requests have zero vector mark in
this conditional calculation. Define g=sum_e g_e and S=sum_e g_e g_e'.

The exact one-step conditional moments of the proposed coordinate increment are

E[Delta x | state] = -epsilon*(q/N)*g,

Cov(Delta x | state) =
epsilon^2 * q*(N-q)/(N*(N-1)) * (S-g g'/N), for N>1.

At exact force balance g=0, this covariance is proportional to S.
These follow by using E I_e=q/N and E(I_e I_f)=q(q-1)/(N(N-1)) for e!=f.
The code checks the formulas by enumerating all 66 two-request samples from
a 12-request test fixture for each candidate. The fixture sizes only validate
the identity; they are not scientific service parameters or a native history.

The native clock has E[Delta tau | state] = q*n/N. If a small displacement
has internal-mode amplitude a_j, the proposed rule linearizes to
E[Delta a_j | state] approximately -epsilon*(q/N)*lambda_j*a_j.
This is a CONDITIONAL first-order drift, not a ratio of random increments,
a pathwise time-change identity, a physical frequency, or a mass law.

A common service pool provides equal inclusion probability per maintenance
token; passive candidate labels do not impose priority or extra work.
A new feedback to F would require specified additional requests, including an
explicit legal-recording interpretation of coordinate changes. The current
calculation introduces no such requests.

## Readout of the prospective forcing directions

The source certificates give zero geometric mismatch for K3/K4. At their exact
compatible coordinates every g_e is zero, so S is exactly zero. Floating
coordinate reconstructions leave raw trace values of order 1e-30; these are
retained separately as roundoff, not interpreted as physical forcing.
Soft-mode shares for zero S are undefined, not reported as zero fractions.

For the saved K5 fits:
- K5-1/2: trace S approximately 1.005848535e-5, and the softest internal mode
  receives approximately 0.344907048% of trace S.
- K5-3: trace S approximately 2.486752393e-4, and the softest internal mode
  receives approximately 0.727032914% of trace S.

The share is v_min' S v_min / trace S, where v_min is unit-normalized in the
same Euclidean coordinate convention. It describes one-step update-covariance
shape under the stated proposal, not the fraction of steady-state energy or
a native measured noise. Relaxation, mode mixing, and subsequent state changes
would also matter in an actual trajectory.

These values are not workload, energy, or mass. No 137/136/factorial or fitted
mass coefficient is applied.

## Executed validation

1. All original spectra, normalization constants, and rigidity ranks replay.
2. Independent Laplacian and edge-block Hessian formulas agree.
3. Per-edge gradients are used to check conditional service mean/covariance
   by direct finite enumeration.
4. K3/K4 exact-zero certificates are distinguished from floating roundoff.
5. Previous input dictionaries remain unchanged.
6. Two executions give identical output and saved numerical results in this
   runtime. Cross-platform last-digit differences are possible.

Files:
- DCU_Mass_Cell_13.py: one copy-paste cell, after Cells 11 and 12.
- DCU_Mass_Cell_13_output.txt: actual output.
- DCU_Mass_Cell_13_RESULTS.json: full decompositions and per-edge directions.
