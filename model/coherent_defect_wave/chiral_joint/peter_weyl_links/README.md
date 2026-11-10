# Finite gauge links and neutral propagation

This construction retains 119 matrix-coefficient states of the native compact
gauge group. Its link matrices include every allowed transition between those
states. They obey the twelve native gauge-current identities, reverse
orientation without conjugating matter, and transport with operator norm one.

The complete three-dimensional finite source has a nonzero neutral-packet
transfer amplitude at fifth order in time, with no additional composite hopping
term. A full-Hamiltonian remainder gives strictly positive arrival probability
and a retained complementary outcome. The input is a prepared physical state;
no ground state, relativistic limit or measured particle mass is identified.

The [proof](PROOF.md) specifies the source, preparation and complete readout.
The [receipt](results.json) contains the exact intertwiners and probabilities.

From the repository root:

```sh
python -m model.coherent_defect_wave.chiral_joint.peter_weyl_links --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/peter_weyl_links
```
