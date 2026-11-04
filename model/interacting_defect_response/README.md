# Interacting localized fermions

A quartic occupation interaction acts on the same localized fermion modes as
the native currents and pixel-dependent mass matrix. It gives a connected
four-point function, shifts spectral poles, and admits a faithful Gibbs state
with a checked noncommuting thermodynamic response. The complete memory-record
limit is uniform in the cube radius for this localized algebra. The interaction
has exponentially decaying spatial tails; a strictly onsite density interaction
fails the exact localized-state inclusion. This finite construction supplies no
relativistic continuum or gravitational entropy law.

The [proof](PROOF.md) states the chosen interaction and its complete scope.
The [receipt](results.json) is regenerated with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.interacting_defect_response --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.interacting_defect_response --check
python -m pytest -q model/interacting_defect_response
```
