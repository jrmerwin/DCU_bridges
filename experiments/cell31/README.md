# Cell 31 — descendant-carrier observers

One Jupyter cell plus an offline, independent-folder reproduction. The native law
is unchanged. Three passive observer policies share each simulated history.
**No physical worldline, particle identity, spatial velocity, or redshift is claimed.**

## In the existing notebook
Paste the contents of `DCU_Mass_Cell_31.py` into ONE cell after Cell30. It reads
`_Native`, `dcu_mass_25`, and `dcu_mass_30`, and writes `dcu_mass_31`. It does not
modify the neutron worker or prior results. Requires NumPy and the standard library.

## In a separate folder, including the Mac
With Python 3.10+ and NumPy available, run:

```bash
python3 validate.py
python3 reproduce.py
```

No internet connection, external accounts, old notebook, or original absolute paths
are needed. The full main experiment took about two minutes in the development
Linux runtime; this is not a Mac performance measurement. The reproduction writes
`REPLAY_*` files and checks every serialized value against the delivered result.

`validate.py` reconstructs private logs, checks receipt and merge certificates,
replays selected native histories using independent bit-mask/sequential accounting,
and verifies the finite sampling laws and branch-identity counterexample.

## Read first
- `DCU_Mass_Cell_31_NOTES.md`: experiment, proofs, results, and physical limitations.
- `DCU_Mass_Cell_31_output.txt`: complete scientific summary, both receiver panels.
- `DCU_Mass_Cell_31_RESULTS.json.gz`: all histories, observer tracks, private event
  cards, auditor metadata, checks and summaries. Private logs and auditor fields
  are expressly separated; do not give an internal observer the whole result dict.
- `PROTOCOL_before_execution.md`: experimental design written before primary execution.

The reference folder contains the original native constructor, original Cell25
source, saved Cell29 emission states, and Cell30 checkpoint/results. The loader
extracts the original native function, verifies the conditioned emissions, and
loads saved prehistories; it does not rerun genesis or claim new mature epochs.

The random-follow policy keeps ONE active carrier. If two reader tracks adopt the
same newborn, their distinct-reader protocol stops. The entire native simulation
continues. Merger observations are retained, not counted as independent spatial
round trips. An unchosen sibling continues to exist and can propagate ancestry.

No model weights, physical constants, fitted coefficients, service routing, or
extra native phase-gate charges are used. The memory/continuation interpretation
is an explicit added readout layer, not derived physical machinery.
