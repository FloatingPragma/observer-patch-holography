# Moving particles that scatter and bind

The source's critically tuned contact interaction survives its three-dimensional
two-particle continuum limit at every bounded total momentum. An explicit
principal-value proof gives uniform scattering errors, and the same coupling
and clock retain the moving bound branch. At a certified finite spacing the
scattering differs from identity by more than one throughout the selected
momentum window. This is a controlled nonrelativistic two-particle result;
a relativistic many-body field theory requires further construction.

The [proof](PROOF.md) includes all eight edge valleys, the interior coarea
estimate and the common-clock join. Its additive extension is self-adjoint
at every particle number and has a lower bound linear in particle number
at each finite spacing, but separated bound pairs make the empty Fock
vacuum fail to be a ground state. No fixed affine energy and momentum map
turns this bound branch into an exact relativistic mass shell. The
[receipt](results.json) binds the source and rational error budgets.
No observed masses or species are fitted.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m model.moving_contact_scattering --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/moving_contact_scattering
```
