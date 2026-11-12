# Finite atomic spectrum controls

The source-generated Coulomb potential is tested with independent changes
of box size and grid spacing. Exact reflection sectors reduce every
computation, while a maximum-principle argument proves the direction of
the fixed-spacing volume change. The [proof](PROOF.md) separates finite
spectral residuals from the continuum and physical-identification errors.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_atomic_refinement --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_atomic_refinement --check
python3 -m pytest -q model/openwave_atomic_refinement
```
