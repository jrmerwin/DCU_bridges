# DCU Cell 34 — exact protected/readable-state audit

**Finding:** With the fixed Cell26 two-path phase attachment, exact nontrivial
protection occurs only for identical register supports. The 36 such path pairs
support a readable three-state phase code, but no distinct-support pair has that
protection. A separately labeled entry-sensitive coupling removes the degeneracy.
This is a conditional quantum-instrument result, not a native particle discovery.

## Run independently

Python 3.10+ and NumPy:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py
```

All data and mathematical-model source dependencies are included. No network is
used by the calculation. The install command may need internet if NumPy is not
already installed. The `REPLAY_` outputs are compared to the delivered reference;
original source/checkpoints are not overwritten. No path from the original
notebook is needed.

Optional independent checks:

```bash
python -m pip install -r requirements-validation.txt
python validate.py
```

For the existing Jupyter notebook, paste `DCU_Mass_Cell_34.py` as one cell after
Cell26 or later. It reads `dcu_mass_26` without changing it and returns
`dcu_mass_34`. It does not launch the native or neutron search processes.

## Files

- `DCU_Mass_Cell_34.py`: complete one-cell analysis.
- `DCU_Mass_Cell_34_RESULTS.json.gz`: full pair classifications, both phase moduli,
  all aggregate checkpoint curves and all pair terminal results.
- `DCU_Mass_Cell_34_output.txt`: primary output.
- `DCU_Mass_Cell_34_NOTES.md`: hypotheses, proofs, source scope and interpretation.
- `DCU_Mass_Cell_34_TESTS.txt`: independent validation actually executed.
- `SCOPE_BEFORE_RUN.md`: finite scope specified before the first calculation.
- `reference/`: unchanged native/phase sources and the original Cell26 data/notes.
- `MANIFEST.json`: hashes of the packaged release files.

Storage is passive in the counting model. Encoding and terminal readout are ideal
instrument operations with no newly derived native cost. No exact spatial or
thermodynamic conclusion is implied by the quantum-memory terminology.
