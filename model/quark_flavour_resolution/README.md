# Unequal quark responses from native hypercharge

A neutral two-state polarization memory splits the native up and down
removal spectra with one shared unit coupling. The charge values come from
the actual fifteen-label representation. Both polarization branches and
every original failure are retained. The complete first energy moment is
unchanged; the second moments acquire a calculable difference of one sixth.

The [proof](PROOF.md) separates this intrinsic source response from a
compensating-reference binding control. The [receipt](results.json) gives
every energy and weight. Neither unit mechanism matches both frozen quark
mass parameters through one scale, and physical quark masses require a
propagation and renormalization map.

```sh
python -m model.quark_flavour_resolution --write
python -m model.quark_flavour_resolution --check
python -m unittest model.quark_flavour_resolution.test_resolution
```
