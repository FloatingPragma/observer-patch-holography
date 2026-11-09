# Spatial Fourier measurements and particle momentum

A coherent array of the selected physical probes has a complete finite
Fourier readout, including failures and the unchanged quantum reference.
Its original probabilities measure spatial phase dependence on the same
source. The [proof](PROOF.md) and [executable](momentum.py) give an explicit
instrument and its finite-aperture error without replacing the prepared
state or closing the source into a torus.

Conserved particle momentum requires a spatial action implemented on the
actual state representation and an attachment to those detector records.
Translation invariance supplies a canonical vacuum-fixing implementation;
other implementations need not fix the prepared vector. The finite
controls distinguish these cases and show why energy cooling alone does
not concentrate momentum.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m pytest -q -W error model/particle_primary_momentum
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m model.particle_primary_momentum --check
```
