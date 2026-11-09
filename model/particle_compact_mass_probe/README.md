# Compact readout of the original source response

A direct algebraic dilation reads the same finite filtered probe with two
neutral qubits and four Ramsey outcomes. Its original success responses
and excitation weights agree with the selector construction; its complete
failure instrument is different. The [proof](PROOF.md) separates this
quantum-memory saving from the large unresolved control and execution
costs. The [receipt](results.json) binds the selected source and an exact
small control. It supplies no particle mass or efficient simulator.

```sh
python -m model.particle_compact_mass_probe --check
python -W error -m pytest -q model/particle_compact_mass_probe
```
