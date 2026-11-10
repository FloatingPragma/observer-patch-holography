# Compensated transport of a native family

The five native gauge representations move under one common interacting
three-dimensional finite Hamiltonian. Each preparation includes a static
dual-charge reference, so the complete state obeys every local Gauss
constraint. The transported matter carries a quantum flux string to that
reference.

At one common source time, every representation has a strictly positive
local detector contrast between preparation on and off. The calculation
retains all quark, lepton, Higgs and link terms, with complete complementary
records and an exact bound on the omitted time-series terms. This is
compensated finite transport; no free-particle dispersion or observed mass
is identified.

The [proof](PROOF.md) specifies the reference and detector support. The
[receipt](results.json) contains the native representation dictionary and
the exact probability intervals.

```sh
python -m model.coherent_defect_wave.chiral_joint.family_transport --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/family_transport
```
