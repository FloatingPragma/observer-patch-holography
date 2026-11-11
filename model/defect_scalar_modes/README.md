# Scalar modes around a compact defect

The same compact texture that traps four fermion modes also selects two scalar
bound modes. The lowest scalar mode has nearly equal coupling strengths, with a
ratio of 1.000122. Applying the full unit Yukawa interaction to the filled native
fermion sea instead places about 96.45% of its scalar response in the odd bound
mode. A variational calculation proves that this sea is not the interacting
ground. Every scalar mode and every native matter label is retained.

The [proof](PROOF.md) gives the scalar normalization, complete measurement
probabilities and strict energy bound. The [receipt](results.json) and
[mode archive](scalar_modes.npz) retain the source calculation. An interacting
particle spectrum and an observed generation count require the full joint
state and propagation calculation.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.defect_scalar_modes --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/defect_scalar_modes
```
