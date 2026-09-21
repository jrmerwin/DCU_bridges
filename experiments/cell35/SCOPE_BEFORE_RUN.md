# Cell 35: later-epoch persistence, before primary simulation

Question: Does later native age improve the frozen record-support phase memory
beyond changing how often it is serviced? Does a fixed small collective extension
have architectures unavailable to the earlier two-input test?

Baseline: September 2026 revised counting law, Gamma=3, m=0, c_first=11,
c_repeat=2, H=0 and H=6. Two-primitive genesis, complete simultaneous enabled
bursts, original relief. Do not change state or schedule using readout results.
Use exact integer uniform sampling; no mean-queue, productive-only clock or
asymptotic age substitution. Compiled implementation must be checked against
original DCUStructure and original Cell25 step given identical explicit draws.

Eight independently seeded trajectories: four per H, seeds 20350920..20350923.
The same seed integers across H are paired initial RNG states, not independent
regime replicas. Each trajectory starts at true genesis and goes to maintenance
age 4,456,448, or a resource cap of 4,000,000,000 iterations / 50,000 objects.
Censoring is retained. Epoch targets: 128, 1024, 8192, 65536, 524288, 1048576, 4194304.
Observe immediately AFTER the complete service/birth burst first crossing each
threshold; report the actual overshoot, never divide summed work by a mean rate.
No target is designated a physical matter/proton/neutron epoch.

At every epoch:
1. Enumerate all current recorded-endpoint primitive-rooted paths of exactly
   3 or 4 composite registers, retaining entry tokens. Select at most 16 unordered
   path pairs per (length,symmetric-difference size) by a separate deterministic
   selector, independent of future service. Retain full counts of eligible pairs.
2. Coverage control: sample at most 64 endpoints first recorded since the previous
   (one-eighth-age) boundary; choose two UNIFORM primitive-rooted paths per endpoint
   using integer chain counts. Keep all outcomes, including duplicates; deduplicate
   exact paths for pair sampling, keep their counts. Match path lengths; at most
   4 pairs per (length,symmetric-difference size). This is a sampled deeper/newer
   path cohort, not all future architectures.
3. Every selected pair retains Cell26/34's two-qutrit +/- coupling, phase moduli
   27 (primary) and 81 (control), and its existing codeword decoder. No echo,
   weight fit, measurement optimization, particle label or energy is introduced.
4. Read at global maintenance increments 128,512,2048,8192 and when union-register
   service count first reaches 8,32,128 (whole bursts, actual overshoots retained).
   Local-activity reads are stopped-process diagnostics, not simulated repeated
   quantum measurements. Save missing milestones, exposure and co-service counts.
   Aggregate within each native seed; do not treat all overlapping pairs as
   independent Monte Carlo samples. Report all overlap classes.
5. Availability-only census of FOUR DISTINCT support sets A,B,C,D with
   1_A+1_B=1_C+1_D, within each path-length panel. Also test entry-token multiplicity
   balance. This is a separately labeled four-input extension, not an exception
   to the two-path theorem. If found, check the inherited two-qutrit code under
   U=D((N_A+N_B)/Q) tensor D(-(N_C+N_D)/Q), no new qutrit or Gamma=5 preparation.
   Added fan-out is hypothetical; native gate/readout work is not compiled.
   Do not interpret a linear cancellation identity as particle stability.

Results support or limit only these mechanisms and covered epochs/architectures.
No exact protected quantum code has been established as a necessary condition
for matter. The typed three-sector condensation simulator in the older registry
paper is a DIFFERENT state law; it must not be imported silently.

Implementation-only timing benchmark at 65,536 ticks completed before the main
cohort. On its runtime evidence (no protection outputs inspected), the terminal
epoch was extended to 4,194,304 with a 262,144-tick follow-up. All eight initial
seeds and analysis rules remain as above. No stopped/failed history is replaced.
