# Finite quantum thermal relaxation

Whole-colour fermion collisions and complete capped-boson jumps relax
the actual calibrated source toward its supplied finite Gibbs state.
The [proof](PROOF.md) gives reference-safe terminal mixing bounds while
retaining every quantum bath coordinate. Executed cooling/heating
responses include full finite occupation laws and arithmetic bounds.

```sh
python3 -m pytest -q model/quantum_thermal_relaxation
python3 -m model.quantum_thermal_relaxation --write
python3 -m model.quantum_thermal_relaxation --check
```
