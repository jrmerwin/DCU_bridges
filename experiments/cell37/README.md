# DCU Cell 37: joint-label compilation and logical renewal

The notebook cell is completely self-contained. Paste all of `DCU_Mass_Cell_37.py` into ONE Jupyter cell and run it. It requires Python3.10+ and NumPy. No previous notebook variables, old folders, network, compiler, LLM, or neutron-search files are needed. The output is `dcu_mass_37`.

For standalone replay, from this extracted directory:

```bash
python3 reproduce.py
```

This writes `REPLAY_RESULTS.json.gz` and `REPLAY_output.txt`, leaving the delivered reference output unchanged. It compares all fields, exactly for discrete values and with a small numerical tolerance for floats.

For independent validation (NumPy and mpmath required):

```bash
python3 validate.py
```

## What is and is not implemented

The same-carrier nondemolition measurement is NOT compiled from native counting operations. The file proves/checks a restriction on independent-local quantum operations and constructs a separately labelled read-and-renew alternative on the protected logical code. Exact native costs apply to creating the fresh sibling pair, not to the complete quantum/controller operation. Missing costs are represented by `None`, not zero.

The early and late source snapshots come from the same previously simulated world. Two conditioned successful resource bursts are evaluated at each epoch; there is no new unconditioned trajectory or measured latency. Every extra birth and outstanding-work increment is retained. Instantaneous readiness and the all-future fixed-anchor bound are reported, rather than assuming a guaranteed resource event.

## Contents

- `DCU_Mass_Cell_37.py`: the single complete cell, including its compressed input fixture.
- `DCU_Mass_Cell_37_NOTES.md`: derivations, native/quantum boundaries, results and scope.
- `DCU_Mass_Cell_37_RESULTS.json.gz`: all outputs, including1160 exact first-step events.
- `DCU_Mass_Cell_37_output.txt`: the printed output.
- `DCU_Mass_Cell_37_TESTS.*`: actual independent-check results.
- `reference/`: original native constructor, original Cell25 service function, Cell36 input fixture, full source node table, and compact Cell37 inputs.
- `SCOPE_BEFORE_RUN.md`: the initial protocol; the all-future readiness bound was added later as an explanatory audit and does not change selection or simulation.
- `reproduce.py`, `validate.py`: standalone tools.

Linux was tested. No Mac runtime or physical validation is claimed.
