# Separating scalar motion from elastic records

The original radial-scalar detector has nonzero positive energy transfer
above a fixed source-clock frequency. An elementary frequency-cut estimate
uses the actual Wilson commutator, and complete finite energy clocks retain
an explicit original probability. The probability floor is extremely small
and its spectral window broad; it supplies no Higgs mass or surviving
particle branch. The same calculation shows that the capacity-uniform
radial detection floor also has a purely elastic counterpart.
The [proof](PROOF.md) gives the bounds and apparatus resources.

```sh
python -m model.particle_candidate_scalar_spectrum --write --check
python -m pytest -q model/particle_candidate_scalar_spectrum
```
