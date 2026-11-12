# Conditional thermal particle populations

The complete calibrated 17-species menu gives temperature-dependent
number, energy and entropy fractions under a declared ideal-gas
equilibrium assumption. Its massless limit has 90 fermionic and 28
bosonic degrees, with energy count 106.75. Separate controls evaluate
the actual 1,653-state capped boson Gibbs law and expose why it cannot
be treated as the uncapped continuum gas.

The [proof](PROOF.md) states the thermodynamic identities and physical
boundaries. The [receipt](results.json) contains the curves and the
comparison with the published interacting equation of state. These
conditional calculations derive no native cosmic temperature or relic
abundances.

```sh
python3 -m pytest -q model/particle_thermal_history
python3 -m model.particle_thermal_history --write
python3 -m model.particle_thermal_history --check
```
