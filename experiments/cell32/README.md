# DCU Cell 32: native timing to quantum correlations

This small branch replaces neither the core DCU nor its frozen physical dictionary.
It temporarily parks Cell31's observer-identity problem and uses the paper's
existing quantum instrument directly.

## Notebook

Paste `DCU_Mass_Cell_32.py` as one cell after Cell29 or later. It requires only
`dcu_mass_29` in the kernel and NumPy. It saves `dcu_mass_32` without modifying old
results or launching another native history campaign.

## Independent directory

With Python 3.10+ and NumPy available:

```bash
python3 reproduce.py
```

The complete saved Cell29 histories are included. No external paths, network,
notebook kernel, LLM, GPU, neutron-search process or empirical constants are needed.
Outputs use the `REPLAY_` prefix and are checked against the delivered JSON values.

For the independent checks (also requires SciPy and mpmath):

```bash
python3 validate.py
```

## Read first

- `DCU_Mass_Cell_32_NOTES.md`: source review, protocol, derivations and scope.
- `DCU_Mass_Cell_32_output.txt`: numerical summary.
- `DCU_Mass_Cell_32_RESULTS.json.gz`: all histories' histogram-derived states,
  joint probabilities, Bell scores, entanglement measures and controls.
- `DCU_Mass_Cell_32_TESTS.txt`: executed independent validation.

The original source functions are included under `reference` for verification.
The standalone primary calculation does not need to rerun the old simulations.

## Scientific boundary

The entangled source, allowed preparation/filter operations, and service-to-phase
rule are declared instrument assumptions. Native service supplies the timing
histories. This is a conditional quantum-channel calculation, not an experimental
Bell test, a derivation of physical gravity, or an SI clock calibration.
The notebook and bundle do not modify the neutron search.
