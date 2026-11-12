# Attraction from a native geometry response

The reciprocal scalar Hamiltonian supplies an exact matter-to-geometry
charge. A declared positive spatial elastic energy turns equal-sign
charges into an attractive force with an approximate inverse-square
profile. Exact bounds on a finite dyadic geometry carrier certify the
attraction even with the source's full nonlinear matter dependence.

This scalar extension also has a precise obstruction to universal gravity:
two preparations above its actual ground have positive excess energies
and opposite geometry-charge signs. The [proof](PROOF.md) identifies the
native coupling and the added spatial law, keeps its free stiffness and
masslessness explicit, and states the finite force bounds. It supplies no
physical Newton constant or native execution of geometry equilibration.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_native_force --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_native_force --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/openwave_native_force
```
