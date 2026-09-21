# DCU Cell33 — fixed quantum echo and timing memory

One notebook cell plus a standalone reproducer. No neutron files are opened.

## Existing notebook
Paste `DCU_Mass_Cell_33.py` after Cell32. It uses `dcu_mass_29`;
`dcu_mass_32`, when present, is also used for regression checks.

## Standalone
Use Python 3.10 or later with NumPy. From this folder:

```bash
python3 reproduce.py
```

Reproducing takes a few seconds in the tested Linux environment. This is not a
Mac benchmark. `REPLAY_*` files are written; frozen outputs are not overwritten.
No Internet access or prior notebook filesystem is needed.

For independent checks, with SciPy and mpmath available:

```bash
python3 validate.py
```

Read `DCU_Mass_Cell_33_NOTES.md` for the protocol, derivations, full tables,
controls, sources, and limitations. `SCOPE_BEFORE_RUN.md` records the starting
plan; subsequent exact two-step audits are identified in the notes.

This is a conditional instrument calculation using 3,072 already-saved native
histories. Pulses are new ideal controls scheduled by external service iterations,
not costed native operations. There is no physical time, energy, or mass calibration.
