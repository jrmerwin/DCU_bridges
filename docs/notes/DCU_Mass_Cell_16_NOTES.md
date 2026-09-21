# Cell 16 — inherited cavity ports with one source-normalized transport budget

## Scope and source basis

This is a new finite spectral experiment. Inputs are the complete registry from Cell 1, the seven frozen cavity candidates from Cell 5, and the exact registry-overlap spectrum from Cell 14. Neither the geometric coordinates from Cells 8–13 nor any empirical mass, density, or neutrino target enters the calculation. The code reads only the candidate catalogue from Cell 5 and the exact spectrum from Cell 14, not their physical predictions.

Source definitions retained:

- `leptons.pdf`, Primitive III, Eqs. (4)–(5): the isolated cavity matrix T_n=-(J-I)/(n-1), its contrast eigenvalue lambda_n=1/(n-1), and the prescribed factor f_n=1-lambda_n.
- `DCU_Registry_Structure_Manuscript(1).tex`: native parentage, exact seven-vertex registry, overlap graph K_11 join (K_63 disjoint-union K_63), and the separately declared overlap readout T_R=-A_R D_R^-1.
- Cell 14: eigenvalues {-1,-62/73,1/136,1/73,383/4964} with multiplicities {1,1,10,124,1}.

NEW hypothesis: a cavity port is coupled to a registry root when that port's native object is in the root's inclusive ancestry. This uses actual ancestry, but adopting that dependency as a bidirectional signal channel is an additional interpretation. It is not an existing native recording or service rule. Co-parent cavity links, registry-overlap links, and cross-ancestry links are three different relations; this model explicitly combines them, rather than claiming they were one graph already.

## The explicit coupling

For a cavity with n core objects u_i and registry roots r, define

    C[r,i] = 1 if u_i belongs to Anc(r), else 0.

Within-cavity and within-registry link weights are one. Cross-link weight is one common SYMBOLIC positive g:

    A(g) = [[A_n, g C^T], [g C, A_R]].

With t_i=sum_r C[r,i] and s_r=sum_i C[r,i], the total degrees are

    D_c(g)[i] = (n-1) + g t_i,
    D_R(g)[r] = d_r + g s_r.

The proposed coupled operator is

    T(g) = -A(g) D(g)^-1.

The old blocks are recovered exactly at g=0. At nonzero g their degrees MUST include the new links. Retaining old denominators while adding feedback is a different operator and does not conserve the old source-normalized budget.

For positive g the graph is connected. T(g) has the exact eigenpair (-1,D(g)1); its symmetric similar matrix is -D(g)^-1/2 A(g) D(g)^-1/2 and its spectrum lies in [-1,1]. This does not make it a unitary wave evolution, a Hamiltonian, a native clock transition, or an energy operator. The eigenvalues remain dimensionless signed-transfer quantities.

No value of g is selected as physics. Equal weight g=1 is not silently imposed. The source documents do not fix the relative strength of these cross-relation links. The calculation instead derives coefficients at g=0, and establishes the leading sign without fitting g.

## Uniform-contact selection rule

Let P_n=I-11^T/n. The isolated cavity contrast family is ran(P_n). If C=b1^T, so all ports have identical registry contact profiles, C P_n=0. For that uniform-column control the cross-injection vanishes for every cavity contrast state. With equal column sums, the degree-normalized cross block also annihilates that subspace at every g.

This is a rule about identical COLUMNS (identical profiles for cavity ports), not merely a rule about the number of registry rows. It does not mean that every coupling to all registry roots vanishes. Source normalization can still change a cavity's internal decay even when the net registry signal cancels.

The actual incidence matrix is not uniform. K3-1 and K3-2 column counts are (23,23,5); K5-1 and K5-2 have (137,74,74,23,5), and K5-3 has (137,74,74,23,23). These counts follow from the frozen core objects, not from assigned particle names. In particular, the common ancestral composite {a,b} in the selected K5/K4 cores is inherited by all registry entries; that explains its 137-entry column without inserting a multiplier.

## Exact registry response and accessible subspaces

The registry operator is self-adjoint in the positive metric W=D_R^-1. Its five spectral projectors can be evaluated as exact polynomials:

    Pi_lambda = product_{mu != lambda} (T_R-mu I)/(lambda-mu).

For Y=C P_n, compute

    B_lambda = Y^T W Pi_lambda Y
             = (Pi_lambda Y)^T W (Pi_lambda Y).

This identity establishes symmetry and nonnegative quadratic forms without relying on floating eigenvectors in a degenerate eigenspace. Each accessible rank is rank(Pi_lambda Y). Their sum is the dimension of span{T_R^k C P_n v} under the UNCOUPLED registry operator, not a count of all modes reachable at finite g in the full coupled system.

The source-normalized model's return term starts at g^2. Exact block elimination gives, where the inverse exists,

    T_eff(z,g) = T_cc(g) + T_cR(g) [z I-T_RR(g)]^-1 T_Rc(g).

This follows directly from solving the registry block of the eigenvector equation. The general finite-dimensional method is the Schur-complement/Feshbach reduction; see Dusson, Sigal and Stamm, arXiv:2105.02058, Theorem 1.2 and Eq. (1.16), https://arxiv.org/html/2105.02058v1 . Its use here supplies algebra, not a physical justification of the coupling.

At g=0, the coefficient of the return term on the contrast space is

    (1/(n-1)) Sigma_n,
    Sigma_n = P_n C^T (lambda_n D_R + A_R)^-1 C P_n
            = sum_lambda B_lambda/(lambda_n-lambda).

All denominators are positive for the tested n=3,4,5, since lambda_n exceeds 383/4964. The inverse has no pole near the unperturbed positive cavity eigenvalue. Matrices and their traces are computed exactly with rational arithmetic; their individual eigenvalues are evaluated numerically only for display.

## Spectral centroid: all branches, not a favorable selected mode

The n-1 contrast eigenvalues are degenerate at g=0 and generally split when the actual port contacts differ. Let lambda_bar_n(g) be their arithmetic mean, following the complete isolated cluster near g=0. This basis-independent diagnostic gives equal weight to every branch. It is NOT a physical state selector furnished by the RMR paper.

Write h_i=t_i/(n-1), hbar=mean(h_i), and overline(h^2)=mean(h_i^2). Then

    lambda_bar_n(g) = lambda_n + a_n g + b_n g^2 + O(g^3),
    a_n = -lambda_n hbar = -sum_i t_i/[n(n-1)^2].

The first-order shift is the source-normalization effect. A registry return requires two cross-links and first enters at second order:

    b_n = lambda_n overline(h^2)
          - [lambda_n/(lambda_n+1)] Var(h)/(n-1)
          + Tr(Sigma_n)/(n-1)^2.

These are, respectively, source-normalization curvature, mixing with the cavity common mode, and registry return. One derivation expands the non-symmetric but similar T(g) directly: T1=U-T0 H and T2=T0 H^2-U H, where U=-A1 D0^-1 and H=D0^-1 D1. The isolated cluster projector is P_n on the cavity block; the common cavity eigenvalue is -1. Applying the usual second-order trace formula gives the displayed expression. A separate calculation using derivatives of the symmetric similar matrix was used for validation.

Define the diagnostic fbar_n(g)=1-lambda_bar_n(g), using the exact unit-modulus global carrier of this source-normalized graph. At zero coupling this agrees with the source f_n. Its leading coefficients are:

- K3-1/2: fbar = 1/2 + (17/4)g + O(g^2).
- K4-1: fbar = 2/3 + (145/18)g + O(g^2).
- K4-2: fbar = 2/3 + (257/36)g + O(g^2).
- K5-1/2: fbar = 3/4 + (313/80)g + O(g^2).
- K5-3: fbar = 3/4 + (331/80)g + O(g^2).

The ratio must correct the reference AND the candidate:

    [fbar_5/fbar_3]/[f_5/f_3] = 1 + k g + O(g^2),
    k = fbar_5'(0)/f_5 - fbar_3'(0)/f_3.

Exactly,

    k(K5-1)=k(K5-2)=-197/60,
    k(K5-3)=-179/60.

Thus this specific normalized coupling with this equal-mode convention lowers the ratio at sufficiently weak positive coupling. It cannot supply an upward PERTURBATIVE correction to a retained leading ratio. Multiplying by a separately assumed fixed 137 does not change this relative sign, but such multiplication is not performed in the code and is not derived by this experiment.

No claim is made about a strong-coupling regime, a differently prepared individual eigenstate, a different contact map, or a general impossibility of mass hierarchy. At finite g the paper does not supply a unique rule identifying one split branch with a particle. Choosing a favorable branch after looking at a target would not establish that rule.

## Executed results

Directly addressable registry roots / accessible contrast-driven registry dimension:

    K3-1 50 / 7; K3-2 50 / 7;
    K4-1 132 / 4; K4-2 136 / 7;
    K5-1 137 / 6; K5-2 137 / 6; K5-3 137 / 7.

The root count is a union over contrast inputs. It is not a count of 137 equally excited channels for one input. The small accessible dimension also does not in itself forbid a large coherent amplitude; the actual response matrix is the relevant computed quantity.

The averaged g^2 registry-return coefficients are approximately 0.202485 (K3), 0.182600 (K4-1), 0.225296 (K4-2), 0.173016 (K5-1/2), and 0.196882 (K5-3). These are dimensionless operator coefficients, not work or mass.

## Validation performed

1. Cell 14 replayed from the original reference generator and archived source package.
2. All 35 spectral projections checked with exact Fractions, including eigenvalue identities, their sum, and weighted orthogonal-Gram identities.
3. Seven independent exact resolvent reconstructions used 3x3 block equations, then checked each resulting vector against the original 137-row adjacency equations. Their projected matrices agree entry-for-entry with the polynomial-projector method.
4. Direct numerical dense resolvents agree with all exact small response matrices.
5. Fourteen full normalized coupled matrices, at g=10^-5 and 5x10^-6, agree with the analytic second-order centroid expansions. These g values are numerical validation points only, not selected model parameters or fits.
6. Full matrices retain the -1 carrier and spectrum within [-1,1].
7. An independent symmetric-matrix perturbation calculation reproduces all second-order coefficients; maximum absolute discrepancy about 1.71e-13.
8. Exact primitive-exchange covariance of the incidence map checked. Cavity-label permutation gives conjugate response matrices and preserves ranks/traces by construction.
9. Two executions produced identical saved numerical/rational fields and printed output in this runtime. Last displayed floating digits can differ on another platform.
10. Previous catalogue, reference, spectrum, calibrations, and native state remain unchanged.

Execution of the supplied cell took approximately 1.5 seconds in this runtime. NumPy is required; SymPy was used only in the separate validation and is not required by the user cell.

## Interpretation and next decision

The calculation supplies a genuinely explicit contact map, its exact accessible registry modes, a budget-consistent coupled operator, and a directionally restrictive perturbation. It does not supply a native strength, mass observable, unitary phase law, or particle-state selection. It does not verify the factor 137 in the source mass prescription or alter the empirically tested global formula subset.

The appropriate next missing ingredient is an independently defined cavity state and contact/phase operation. The isolated graph order alone is no longer enough once its degenerate modes split. Neither a favorable mode nor a coupling strength should be selected from the empirical discrepancy.
