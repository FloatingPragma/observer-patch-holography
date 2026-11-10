# Paired charged source spectra

The actual finite native Hamiltonian has an exact antiunitary symmetry
that exchanges source addition and removal, reflects space and reverses
the time parameter. Their source spectral measures match on the unique
interacting vacuum. The selected odd doublet carries the finite
binary-cube representation G1. The retained static reference energy is
outside this symmetry. The [proof](PROOF.md) keeps this reference boundary
explicit; continuous physical spin and physical CP or CPT are separate
claims.

Run from the repository root:

```
python -m model.coherent_defect_wave.chiral_joint.interacting_charge_conjugation --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/interacting_charge_conjugation
```
