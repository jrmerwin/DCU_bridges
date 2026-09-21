# Cell36: protected collective memory and an explicit state-sensitive interface

## Run

Paste the complete `DCU_Mass_Cell_36.py` as one Jupyter cell. It requires NumPy
and the Python standard library. It is self-contained, so no earlier kernel
variables, external file structure, compiler, network, model weights, or neutron
search folder are needed. Results are returned in `dcu_mass_36`.

Alternatively, from this extracted folder:

    python3 reproduce.py

This recomputes the finite analysis and compares every scientific result with the
bundled reference, then writes separate REPLAY files. The reference is unchanged.
Main runtime was about0.2s in the development Linux environment with one BLAS
thread; this is not a Mac benchmark. No multi-million-tick histories are rerun.

For independent validation (also needs mpmath):

    python3 validate.py

## What the test establishes

* The existing four-path attachment protects its code under native storage.
* The native birth/work/queue transition cannot distinguish the three phase words
  at fixed structure. This is not a claim that all quantum readouts are blind or
  that all different candidate architectures have the same workload.
* Separate local projective Fourier measurements and a coherent coarse joint
  measurement have identical reported-label probabilities but different effects
  on the stored memory.
* An explicitly ADDED probe-controlled common cyclic shift realizes a
  nondestructive measurement of the three orthogonal phase words. A separate
  ideal phase command changes the label without leaving the protected space.
* Splitting that probe interaction into two controls exposes the state to native
  phase activity between them. The completed operation can remain in the code
  space and nevertheless corrupt the label. All gap controls are retained.

This is NOT an autonomous native quantum meter, native interaction-cost
calculation, mass generation, neutron identification, binding law, spatial
manifold or SI calibration. Preparation, controlled interaction and physical
readout remain added assumptions. The existing probability rule alone does not
select a post-measurement state update.

## Contents

`DCU_Mass_Cell_36.py`: complete self-contained numerical cell.
`DCU_Mass_Cell_36_NOTES.md`: derivations, sources, status and caveats.
`DCU_Mass_Cell_36_RESULTS.json.gz`: all scientific output.
`DCU_Mass_Cell_36_output.txt`: human-readable execution output.
`validate.py`: original-model and independent high-precision checks.
`reproduce.py`: fast complete replay and field-by-field comparison.
`reference/fixture36.json`: small projection of the actual Cell35 saved inputs.
`reference/node_tables/`: complete eight original archived node tables.
`reference/native_constructor.py`, `reference/DCU_Mass_Cell_25.py`: unmodified
original model functions, read by the validator without launching old simulations.

The archive supplies no runtime environment or neural-network weights; none is
used by this calculation. All code and data are relative to the extracted folder.
