# Joint energy and spatial-phase spectroscopy

The selected finite source admits a complete experiment that records
initial energy, spatial Fourier outcomes and final energy. Every failure,
energy bin and the original quantum reference is retained. The
[proof](PROOF.md) gives alias, finite-clock and source-ground errors and
a positive broad joint band for the actual plaquette aperture.

A resolved dispersion branch requires measured original probabilities
and a physical momentum attachment. The broad proved band does not
identify a particle or a rest mass. These finite comparisons need
controlled errors on the specified member, without requiring an infinite
continuum theorem.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m pytest -q -W error model/particle_primary_joint_spectroscopy
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m model.particle_primary_joint_spectroscopy --check
```
