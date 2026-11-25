# Scalar vacuum readings on a spatial ball

The native routing word transports a marked excitation with the dispersion of a massless scalar. Its protected single-particle states have at most two Schmidt coefficients across a station region. Canonical quantization of the positive frequencies supplies a separate Gaussian model. On that model smooth momentum excitations give the conformal ball weight, while a smooth field displacement gives an additional contribution proportional to the surface area. Finite regional instruments distinguish the effects. [Construction and proofs](PROOF.md) · [Gaussian readings](measure.py) · [Native routing comparison](compiled.py)

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.vacuum_kernel.compiled
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.vacuum_kernel.measure
OPENBLAS_NUM_THREADS=1 python -m pytest -q code/matched_einstein_limit/source_implication/vacuum_kernel
```
