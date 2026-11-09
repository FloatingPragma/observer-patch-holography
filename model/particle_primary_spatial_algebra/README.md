# Continuous spatial action and the momentum premise

Spatially smeared gauge-invariant charge and Wilson-loop probes give a
regular observable algebra with continuous R3 translations commuting
with the selected source-time action. Its stationary state and positive
time-energy GNS representation use the same finite Wilson/Gauss sequence.
The [proof](PROOF.md) compares bulk records with a physical infinite-lattice
algebra by an explicit boundary locality bound and identifies their
limiting cyclic word moments.

Canonical particle momentum requires translation invariance of the same
selected limiting state, which this construction does not prove. The
regular probe menu preserves the source preparation; bounded spatial
averaging need not produce a surviving field or particle response.

`spatial.py` provides exact spatial weights, rounded-translation bounds,
open-boundary interaction counts and authenticated bulk transfer plans.
The source matrices are not allocated. Independent controls show why
state invariance and the distinction from cyclic finite-box labels matter.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m pytest -q -W error model/particle_primary_spatial_algebra
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m model.particle_primary_spatial_algebra --check
```
