# Cell 34 — architecture, protected subspaces, and readable distinctions

## Result in one paragraph

Under the unchanged Cell26 support-count phase attachment, a pair of nonempty
record supports has a protected subspace of dimension at least two if and only if
the supports are identical. The four saved panels contain 1,913 unordered path
pairs; 36 have identical supports and maximal protected dimension three. Every
one of those 36 pairs follows the same composite chain, differing only in its
primitive entry to genesis register 2. After quotienting by register support,
there are zero protected pairs. The three-dimensional zero-difference subspace
supports three orthogonal, coherently encoded states which the existing local
Fourier readouts distinguish perfectly through their joint outcome difference.
Thus this is a conditional readable shared-input memory, not additional particle-
specific protection. A separately labeled entry-conjugated phase coupling removes
even those degeneracies. No source mass recipe, native law, or neutron search is
changed.

## 1. Source and assumption boundaries

Source: *Self-Limited Growth, Intrinsic Time, and Record Geometry in the
Distinction Combinatorial Universe*, revised manuscript 0.5 (September 2026):
Sections 2.1–2.4 define objects, recording and global service; 5.1 defines register
paths and entry reflection; 7.1 supplies the conditional shared-source quantum
state and phase measurements. This is the same source used in Cells25–33.

Cell26 specifies the exact fan-out wiring used here: for two classical record
paths A and B, their support service counts drive a separately prepared TWO-qutrit
apparatus by D(N_A/Q) tensor D(-N_B/Q). It does NOT establish that the paths themselves
are quantum tensor factors or that this wiring/readout is physically costless.
The attached state is assumed unaffected by native events except the declared
phase drives. These are substantive scope restrictions.

Q=27 remains primary, Q=81 the existing control. Both use G=diag(0,1,2) and
D(x)=exp(2*pi*i*x*G). No rescaling, depth-dependent coefficient, state-dependent
service priority, fitted noise, echo or storage feedback is introduced.

The general decoherence-free-subspace comparison is external quantum theory:
D. A. Lidar, I. L. Chuang and K. B. Whaley, *Decoherence-Free Subspaces for Quantum
Computation*, Phys. Rev. Lett. 81,2594 (1998), DOI 10.1103/PhysRevLett.81.2594;
arXiv:quant-ph/9807004. The finite proofs below are derived directly for the
supplied DCU coupling. No Markov-semigroup assumption is imported for the reduced
quantum apparatus.

## 2. Exact noise algebra

Let I_A(z) be the indicator that register z belongs to A. Its single-service gate is

    U_z = exp[(2*pi*i/Q) H_z],
    H_z = I_A(z) (G tensor I) - I_B(z) (I tensor G).

In computational basis |r,s>, r,s in {0,1,2}, the integer character is

    h_z(r,s) = r I_A(z) - s I_B(z)  (mod Q).

A complete service burst applies the product of its U_z. A whole service history
applies a product of such commuting gates. Consequently, common character
classes give exact invariant scalar-action subspaces at every time, for every
history, with no assumptions about independence or stationarity of increments.

These classes are also NECESSARY for this all-history protection test at the
actual observation checkpoints. F is 691 or 339, n is 18 or 32, and q=3. An
all-forced draw gives identity with positive probability

    choose(F,3) / choose(F+n,3).

For every old register z, selecting only its maintenance request plus two forced
requests is a legal draw with positive probability

    choose(F,2) / choose(F+n,3).

The native constructor is idle for this singleton maintenance selection: there
is no omitted child or extra pair. Thus each individual U_z, as well as identity,
is actually in the one-step noise support. A subspace whose every pure-state
density matrix is preserved under every supported draw must see each U_z as a
common scalar; averaging cannot restore all such pure states if distinct branches
separate them. This is stronger than finding a nearly degenerate finite-sample
channel. We test exact common-character subspaces, not a recovery-assisted code.

Because Q is 27 or 81 and local eigenvalue differences are small integers, no
alias due to numerical rounding or coarse modulo-three resolution enters.

### Identical nonempty supports

Every nonzero generator is G tensor I - I tensor G. Its eigenvalues are
-2,-1,0,1,2, with respective multiplicities 1,2,3,2,1. The maximal protected space is

    C = span{|00>,|11>,|22>}.

Every U_z is the identity on C. Every density matrix supported there is preserved,
including its off-diagonal coherences, not merely its populations. The two
additional dimension-two eigenspaces acquire common global phases.

### Distinct nonempty supports

If both sets have private registers, their generators independently resolve r
and s. If one set were a proper subset of the other, a private register plus a
shared register would do the same. Thus the common eigenspaces are one-dimensional.
In our equal-length panels, distinct supports always have private registers on
both sides. Preservation of isolated computational basis states is not a protected
qubit/qutrit. This conclusion does NOT rule out partial finite-time coherence,
active error correction, alternative couplings, or larger encoded architectures.

## 3. Native architecture census

All paths are the original primitive-rooted directed placement records ending
at already-recorded objects. No later record status is imported when selecting
this panel. They are independently regenerated from the saved parentage and
also independently enumerated backwards by the validator.

| Environment | Path length | Paths | Distinct supports | Unordered path pairs | Maximal-dim-3 pairs |
|---|---:|---:|---:|---:|---:|
| H0 | 3 | 17 | 13 | 136 | 4 |
| H0 | 4 | 17 | 12 | 136 | 5 |
| H6 | 3 | 42 | 29 | 861 | 13 |
| H6 | 4 | 40 | 26 | 780 | 14 |

Total: 116 path appearances, 1,913 comparisons, 36 equal-support pairs. All 36
pairs differ only in the primitive entry through composite 2; the subsequent
composite chain and entry tokens coincide. Different panels/paths are correlated
structural cases, not independent physical observations.

After collapsing equal supports, the panels have 78,66,406,325 distinct-support
pairs (875 total). None supports a common-character subspace of dimension two.
There is no statistical long-running search needed for an exception to this
pair-level classification under the frozen law.

## 4. Readability: retained coherent information, not just an inert state

Use three codewords already expressible through the adopted source and phase
family:

    |Omega_k> = sum_(r=0)^2 exp(2*pi*i*k*r/3)|rr>/sqrt(3), k=0,1,2.

Encoding is D(k/3) on Alice BEFORE storage. It is an ideal declared preparation
operation; its physical/native work is not supplied. The three states are mutually
orthogonal, have identical computational populations, and identical single-factor
marginals I/3. Their distinction resides in joint coherence.

At the common zero-setting local Fourier measurements, the paper's convention
gives

    p(a,b|k) = 1/3 when b-a = k (mod 3), and zero otherwise.

Thus joint outcome comparison reads k exactly on an equal-support protected pair,
without the service history, without an echo and without feedback during storage.
The two local outcomes separately are uniform; no superluminal/local signaling
claim is involved. Combining the results and physically implementing preparation,
memory, and terminal measurement remain additional operational tasks.

For distinct supports and accumulated count difference D=N_A-N_B,

    p(decoded=j | encoded=k, history)
      = [1 + 2 cos(2*pi*(D/Q+(k-j)/3))]^2 / 9.

At j=k this is exactly the Cell26 kernel. ALL original overlap-class means are
reproduced, not refitted. In particular, terminal primary-length-3 correct-word
means (shared registers 3,2,1,0) are:

    H0: 1, 0.602864287, 0.467629511, 0.409071019.
    H6: 1, 0.671111700, 0.530541465, 0.456325579.

These are new interpretations/checks of the existing simulation readout, not new
independent numerical successes. All 512 histories, 17 checkpoints and both phase
moduli are retained. Under complete computational dephasing the codewords become
the same separable mixture, and this decoder succeeds with probability 1/3. That
negative control confirms that the perfect decoder uses coherent information.

The code's entire three-dimensional linear span is mathematically protected, but
we do not claim that arbitrary logical input preparation or universal logical
operations have been compiled from native constructions. Three accessible ideal
phase codewords and exact linear-subspace preservation are the demonstrated tests.

## 5. Entry information: an explicitly different coupling as sensitivity control

The native record convention includes S=(0 1), fixing label 2. Cell26's primary
phase driver ignores which parent entry the path uses. We retain that primary
rule. For a SECONDARY, separately named sensitivity test, define

    G_e = S^e G S^e,
    H_z^entry = I_A(z) G_(e_A(z)) tensor I
                - I_B(z) I tensor G_(e_B(z)).

This uses the actual permutation (eigenvalue list (1,0,2) when e=1), NOT an assumed
sign (-1)^e. Indeed S G S is not -G plus a scalar: the three sums of the original
and reflected eigenvalues are (1,1,4), not a constant.

This is a different physical coupling hypothesis while the source/preparation/
terminal measurements are held fixed. It is NOT a mere coordinate relabeling;
a global change of basis would also transform states and measurements and leave
predictions unchanged. Neither the DCU paper nor this audit selects this secondary
wiring as the correct physical one.

For a shared-support pair, the paths enter register 2 differently but subsequent
registers identically. On |rr>, the extra register-2 phase coefficients are
(-1,+1,0) (up to exchanging the two paths). They split the three code basis states.
Across all 1,913 distinct path pairs, the secondary common-character classes are
one-dimensional. No storage-protected qubit/qutrit survives this sensitivity model.

Illustratively, for the first length-three alias, the same saved histories give
terminal primary-Q=27 correct-word probabilities 0.090157660 (H0) and 0.151684098
(H6), versus exactly one under the unchanged primary coupling. These low numbers
are finite-window phase/readout results, not monotone decay laws, absolute memory
capacities, or a fitted failure threshold. The exact character splitting—not a
small numerical probability—is what removes the subspace protection.

Thus shared native ancestry alone does not guarantee the tested protection under
every state-sensitive interaction. Equality of the COUPLING seen by the encoded
states is essential. This is the key scope limit, not an excuse to rewrite the
baseline after seeing results.

## 6. What has and has not been learned

Positive: an explicitly readable, three-dimensional protected state space exists
under the existing shared-input attachment. Its protection is exact for every
native history under that fixed rule. The native inputs and their ideal quantum
readout are mutually consistent.

Negative for the proposed stronger interpretation: all pair-level protection is
explained by identical driving supports. There is no extra protection linked to
distinct record architectures, a registry sector, K3/K5 order, particle species,
energy or space. The paper's persistent trits carry distinct entry information,
but they are not established protected quantum amplitudes by this test.

As specified before execution, this is a stopping result for the two-path
support-only branch. It does not justify a longer search for rare exceptions,
and it does not identify matter or a neutron. Richer interactions, larger codes,
state-dependent energy, and physically costed preparation/readout are separate
questions. No new source assumption is imported from the older RMR mass/thermal
papers or from the unrelated typed sector simulation.

## 7. Reproduction and validation

Notebook: paste DCU_Mass_Cell_34.py after Cell26 or any later cell while dcu_mass_26
remains available. No Cell32/33/31 or neutron scanner execution is needed.
Standalone: Python3.10+ and NumPy; run `python reproduce.py` in the extracted
bundle. It reads the frozen Cell26 data, writes REPLAY files, and compares every
output field to the saved reference. No internet, old filesystem, pickle, LLM,
neural model or empirical constant is needed.

Main calculation: 7,652 exact pair/model/resolution classifications, 156,672
pathwise equal-support count snapshots, 612 old readout regressions, 1,296 dense
encoded-state tests, and 864 logical matrix-unit checks. No new primary native
simulation was run.

Independent validator additionally uses SciPy/mpmath. It checks 6,477 complete
legal maintenance subsets, all 3,826 pair/model classifications against full
burst-character tables, 50 real native singleton events, all 116 paths backwards,
132 independently exponentiated generators, and 972 decoded probabilities at
70-digit precision. It replays 16,384 complete native transitions using both the
original constructor and a separate bitmask/sequential-charge implementation,
recovering 128 archived count checkpoints. These replays are validation, not a
new statistical cohort. See TESTS.txt for actual errors and runtime.

Numerical tolerances apply only to evaluating exponentials and matrix operations.
All protection dimensions and equivalence classes use exact integer arithmetic.
No precision threshold is presented as a physical accuracy claim. Linux has been
tested; this is not a Mac benchmark.
