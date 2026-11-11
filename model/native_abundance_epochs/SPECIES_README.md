# From class counts to particle counts

A centre-class census fixes only sums of populations sharing a label.
The [exact calculation](SPECIES_PROOF.md) identifies those sums and gives
sharp bounds for unresolved species. Equal weights give a 75/25 quark/lepton
split in the 32-component menu and 80/20 after omitting the neutral singlet
pair. Neither weighting is a dynamical abundance. A separately declared
ideal thermal comparison includes physical helicities and distinct number
and energy weights. Applying it to native flux faces requires an implemented
species decoder and thermal preparation.

```sh
python3 -m model.native_abundance_epochs.species_identifiability --check
python3 -m pytest -q model/native_abundance_epochs/test_species_identifiability.py
```
