# The interacting ground and its first collective modes

An exact spin-coordinate transformation reduces the controlled
second-order interaction to a cube antiferromagnet. Its complete
256-state spectrum and the full-source error bound prove a unique
interacting physical ground at explicitly bounded nonzero hopping.
Its first three collective excitations occupy three distinct source
momenta. Their finite energy gap is not a particle rest mass.

The [proof](PROOF.md) retains all charged intermediate sectors, the
scalar energy shift and the distinction between the admitted hopping
coefficient and unit hopping.

```sh
python -m model.coherent_defect_wave.chiral_joint.atomic_spin_geometry --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/atomic_spin_geometry
```
