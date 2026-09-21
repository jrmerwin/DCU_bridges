# Cell 35 — maintenance-age and activity-matched persistence

## What was actually tested
Eight new native worlds (Gamma=3,m=0,H=0 or6, four seeds each), from true
2-primitive genesis through 4,456,448 maintenance ticks. Record paths were
selected at 128,1024,8192,65536,524288,1048576,4194304 ticks. No age is designated a
physical matter epoch. The counting law has not been changed. The phase-memory
attachment and ideal preparation/readout remain assumptions.

The frozen two-path decoder was compared both at equal global-tick windows and
at equal union-register service counts. A separate four-path input extension
checks exact reconverging-path balance. A longitudinal confirmation uses the SAME
early support pair in its early/late environments: 32 fresh continuations at each
of16 checkpoints, 512 in total.

## Fast analysis replay (no compiler)

    python replay_readout.py

Requires NumPy. Checks the saved complete record paths, all stored counter/readout
identities, the four-path certificates, and the 512 fresh confirmation readings.
This replays analysis, NOT native generation. The reference results remain intact.

## Full reproduction / notebook
Paste the entire DCU_Mass_Cell_35.py as ONE Jupyter cell, or run:

    python DCU_Mass_Cell_35.py

Requires Python3.10+, NumPy, and a C++17 compiler (c++, clang++, or g++).
The C++ source is embedded in the single cell. On a Mac, Apple's command-line
development tools supply clang++ when installed. No LLM, model weights, network,
prior notebook definitions, or neutron-search filesystem is used.

Full execution creates `DCU_Cell35_run/` in the current working directory. It uses
2 workers by default and ~2.3GB of scratch event logs plus results. Allow several
minutes or longer depending on the machine and other work. Each native world is
bounded by4billion service iterations and50,000 objects; incomplete ages/windows
remain censored. Do not reinterpret a censored window as a stable state.

The top-level function accepts `workers=1..8` and `include_replication=False` for
just the main age ladder. Neither option changes an existing native world. The
native core and retained source/seed are checked before reusing cached files.
An interrupted complete-world generation may need regeneration; this is a bounded
experiment, not the resumable neutron scanner. No interaction with that scanner
is performed.

The package's separate `run_replication.py --run DCU_Cell35_run` repeats only the
fixed-support confirmation using the generated raw data. It is NOT required when
the single cell is run with its default include_replication=True.

## Files
- DCU_Mass_Cell_35_NOTES.md: methods, exact laws, caveats and source locations.
- DCU_Mass_Cell_35_output.txt: complete printed scientific summary.
- results/: eight full per-world analyses and parentage, selected readouts, and
  combined compressed results. No enormous raw event logs are included; full
  reproduction regenerates them deterministically with a portable integer RNG.
- replication/: all16 source-state inputs and512 fresh continuations' readings.
- reference/: unchanged original native constructor, Cell25 source, Cell34 source
  and prior exact-age audit.
- native35.cpp: compiled implementation used for this calculation.
- SCOPE_BEFORE_RUN.md: original prospective scope and implementation-only extension.
- DCU_Mass_Cell_35_TESTS.txt/json: actual validation outcomes.

The new RNG is a specified SplitMix64 stream with rejection-based bounded integer
sampling. It is not Python's old random stream. Given the same explicit sampled
tokens, complete transitions are checked against the original Python code.
No Poisson clock, float hazard, mean age inversion, or approximate queue skip
replaces native updates. The compiled code is checked for integer overflows.

## Interpretation
The two-path exact obstruction remains. Greater longevity in global ticks mostly
reflects fewer services reaching a fixed finite support. It is NOT a general
absence-of-matter theorem. Four distinct paths can preserve the code through
balanced aggregate drives, including balanced entry tokens, but that is still a
hypothetical readout wiring, not an autonomous particle, binding law or energy.

All calculations were performed in Linux. This package has not been benchmarked
on Jason's Mac. Do not equate these Gamma3 trajectories with a different search
law or with a physically calibrated cosmological epoch.
