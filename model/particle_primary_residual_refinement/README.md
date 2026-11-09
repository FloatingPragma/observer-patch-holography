# Local moments and the interacting trial's error

Exact particle-hole and Haar contractions reduce the smallest primary
member's full-Hamiltonian trial residual from about102.09 to at most5.232
source energy units. The eleven-vector trial and its ground upper bound
are unchanged. The [proof](PROOF.md) retains the entire omitted space;
at link degree two or higher its squared residual formula is exact.

This controls the auxiliary trial's spectral proximity. A particle mass
requires an observable branch and its original selected-state weight;
neither is assigned by this residual window.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m pytest -q -W error model/particle_primary_residual_refinement
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /tmp/oph-completion-venv/bin/python -m model.particle_primary_residual_refinement --check
```
