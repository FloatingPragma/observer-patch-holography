# Charged states in a bounded radial source

The same bounded radial Hamiltonian supports exactly two localized
one-body zero modes and a controlled interacting charged core. Its
original centered contact gives an actual ground doublet, positive
addition/removal poles and a large contrast in the central occupation
parity. A faithful finite Gibbs reference approximates this ground with
an explicit error and modular clock. The supplied texture and finite-box
interaction bound are explicit; a volume-uniform relativistic particle
theory requires a separate limit.

The [proof](PROOF.md) derives the radial source symmetries, filled-sea
contact, complete-Fock error bounds and actual local observable.
The [receipt](results.json) binds the source and independent controls.

```sh
python -m model.radial_interacting_core --check
python -m pytest -q model/radial_interacting_core
```
