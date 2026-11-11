# Pion masses

The archived quark-antiquark correlator supplies an isospin-degenerate pion
triplet. With the measured rho mass and the imported pion/rho tuning ratio,
its conditional mass is **0.13799628 GeV** for π⁺, π⁰ and π⁻. The alternative
string-tension scale is also retained. Electromagnetism, unequal light-quark
masses and a controlled physical limit are required for the observed split.
See the [proof](PROOF.md) and [receipt](results.json).

```sh
python3 -m model.pion_mass_readout --write --replay-archive
python3 -m model.pion_mass_readout --check
python3 -m unittest model.pion_mass_readout.test_readout
```
