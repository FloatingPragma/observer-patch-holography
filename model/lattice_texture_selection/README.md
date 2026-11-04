# Which spin textures cost finite energy?

On the source's selected ordered quantum ground, nonzero winding on each
admissible cubic shell costs at least `pi J kappa` per unit degree.
The bound uses actual nearest-neighbor bonds and an exact integer degree
calculation. Nonacute shells also have a positive cost, so finite energy
forces zero winding on every sufficiently distant shell. Every finite assignment has an explicit
controlled unwinding path; the result is an energy constraint, not an
exact particle label. The [proof](PROOF.md) states the lattice admissibility
condition and keeps repair-only records separate from Hamiltonian evolution.

Connecting these spin textures to the prescribed fermion mass background
requires a formation law and a common field map. No particle mass or defect
lifetime is inferred from the shell bound.

```sh
python -W error -m model.lattice_texture_selection --write
python -W error -m model.lattice_texture_selection --check
python -m pytest -q -W error model/lattice_texture_selection
```

The [receipt](results.json) retains the exact shell examples, source bond
identities, finite deformation controls and complete native repair checks.
