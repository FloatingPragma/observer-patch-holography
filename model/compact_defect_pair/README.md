# Finite-energy defect pairs

A compact texture on the two-spin lattice carries opposite winding numbers on
two small surrounding shells and becomes constant outside a finite cube. Its
spin-texture energy and fermion vacuum-energy difference have bounds independent
of the exterior volume. Four localized modes form an isolated band with a
controlled infinite-volume limit. An exact certificate shows that none has zero
energy, despite the two winding numbers. The same native currents and pairing
operator act on that band. Autonomous formation, stability and a particle rest-mass
interpretation require additional dynamics.

The [proof](PROOF.md) states the design, spectral and energy bounds. The
[receipt](results.json) retains every finite eigenvalue and its complete-certificate
control.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.compact_defect_pair --check
python -m pytest -q model/compact_defect_pair
```
