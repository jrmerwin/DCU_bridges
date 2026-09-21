# Cell 14: exact global spectrum and conditional cross-sector closure

## Scope

This is a consolidation calculation, not a new native simulation, spectral mass
model, empirical fit, or dynamical coupling. It uses the complete registry from
Cell 1. Cells 6–13 are not prerequisites. Earlier geometry, clock definitions,
source files and physical calibrations are unchanged.

The code separates three objects:

1. The exactly enumerated ancestry registry and overlap graph.
2. The previously declared signed overlap operator T = -A D^(-1).
3. A fixed subset of RMR/taxonomy physical formula prescriptions.

The code does not derive (3) from (2) or identify (2) as native service.
It tests the internal algebra of one conditional joint formula subset, not the
mutual consistency of every formula/version in the supplied papers.

## Native counting does not gain a repair trigger from passive readouts

In the revised DCU counting law, Z=(X,F,P), together with fixed Gamma, m, H,
determines the distribution of service and the subsequent construction/work
update. First/repeat recorded status is determined by X. Attached coordinates,
metric mismatch, and realized ternary register values do not enter this transition.
The record instrument is expressly outcome-blind with respect to this counting law.

Consequently, starting from the same native counting state and using the same
random service choices produces the same future native counting states, regardless
of the values of passive attached readouts. This follows inductively through each
service/relief/construction update. It excludes an additional causal repair
surcharge from those readouts WITHOUT changing the law. It does not exclude
structural dependence of ordinary work (X, path weights, ancestry hits already
matter), every possible mass proxy, or an explicitly extended physical model.

Ancestry-closed ideals form a distributive lattice under set union/intersection:
I intersect (J union K) = (I intersect J) union (I intersect K).
Different path lengths into a shared vertex are not inconsistent assignments:
Ch(z) adds the parent path counts. Reusing a factor in distinct construction
events is already counted, whereas reconstructing an existing parent pair is
forbidden. Incompatible requested address values for one persistent register have
zero joint probability in the stated record instrument; the rule does not create
a second copy of that register to repair them.

## Exact spectrum of the declared overlap operator

Write O = K_u join (K_b disjoint-union K_b), with observed u=11, b=63.
The universal degree is dU=u+2b-1=136 and each branch degree is
dB=u+b-1=73. The actual adjacency is checked against these blocks entry by entry.

Within one block, every zero-sum vector is an eigenvector of T with eigenvalue
1/d for that block. This yields 10 universal modes at 1/136 and 124 branch modes
at 1/73.

On vectors constant per block (universal, branch A, branch B), the operator is

    [ -(u-1)/dU  -b/dB       -b/dB     ]
    [ -u/dU      -(b-1)/dB   0        ]
    [ -u/dU       0         -(b-1)/dB ]

with eigenvectors and eigenvalues:

    (dU,dB,dB)       -> -1
    (0,1,-1)         -> -(b-1)/dB = -62/73
    (-2b,u,u)        -> 1-(u-1)/dU-(b-1)/dB = 383/4964.

The within-block-zero-sum and block-constant spaces intersect only at zero.
Their dimensions are 134 and 3; the latter three eigenvalues are distinct.
The constructed vectors therefore form a complete eigenbasis. The program
checks every eigenvector equation with rational arithmetic against the actual
137-entry adjacency, and checks both trace(T)=0 and trace(T^2).

This spectrum is not the spectrum of K137 or of a separate K3/K5 cavity.
No entry by itself determines a particle mass or an excitation energy.

## Conditional global formula subset and N elimination

Source prescriptions held fixed:

    r_mu = m_mu/m_e = N*(3/4)/(1/2) = 3N/2
    omega = Omega_m = I/(N-1), I=40
    d_nu = Delta m21^2/m_e^2 = 1/[4 N (N-1)^6]
    R_nu = Delta m32^2/Delta m21^2 = 3^2+5^2-1=33.

The first is Eq. (9) of the supplied lepton paper. The latter cosmological and
absolute-splitting prescriptions are from the supplied unpublished taxonomy.
This subset is explicitly conditional: the taxonomy also uses a different
phase-space muon formula elsewhere. No numerical or conceptual reconciliation
of those alternate prescriptions is performed here.

Eliminating N yields:

    omega*(2 r_mu-3) = 3 I = 120
    d_nu*8 r_mu*(2 r_mu-3)^6 = 3^7 = 2187
    d_nu = omega^7/[4 I^6*(omega+I)].

These are algebraic consequences, NOT three extra confirmations of the model.
They restrict how one could change a common N while preserving these formulas;
they do not authorize fitting N away from the enumerated registry count.

The neutrino formula fixes a scale only relative to m_e. Specifying the full
spectrum additionally uses normal ordering and m1=0. Conversion to eV needs a
physical mass anchor. No experimental target or scale is read in Cell 14.

The following remain unproved physical/model identifications, not count outputs:
selection of an interface element, the physical use of N-1, the muon assignment
and selective registry coupling, sampling-to-cosmological-density mapping,
the depth exponent six and its operation, the splitting rule, and the lightest
neutrino state. Sharing integers does not eliminate these assumptions.

## Executed tests

- Reference generated from the original archived F source definitions.
- Two complete executions: identical exact result dictionaries.
- All 137 eigenvector equations checked exactly.
- Invariant-subspace dimension/spanning check, trace and squared-trace checks.
- All three formula-elimination identities checked exactly.
- Reference input dictionary compared unchanged after execution.
- Independent symbolic characteristic polynomial of the 3x3 block quotient
  and numeric full-adjacency eigenvalues checked in supplementary validation.

Only Python's standard library is required by the delivered cell. SymPy/NumPy
are used solely for the supplementary validation, not for its main certificates.
