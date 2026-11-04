# Interacting readouts of a compact defect pair

The compact finite-energy texture, native currents and memory-derived pairing
operator share one 180-mode interacting core. Its nonzero spatial splitting
changes the ground sector and measured excitation energies of exact charged
four-mode slices. Complete native records approach the same nonlinear mean
dynamics with an explicit bound uniform in exterior volume. These finite
readouts require a common physical vacuum, propagation and unit map before
they can be identified with particle rest masses.

The [proof](PROOF.md) gives the exact spectrum, rational certificates and
operator inclusion. The [receipt](results.json) retains the complete
sixteen-state spectra, pole weights and source dependencies.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.compact_interacting_response --check
python -m pytest -q model/compact_interacting_response
```
