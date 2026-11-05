# Spatial structure and energy flow of charged composite responses

The centered cube's 138 composite spectral lines each have a certified
four-dimensional spatial response. Exact spin and charge ladders give
2,208 visible states per charge sector. Together with single-fermion
probes, these observations cover 4,032 of 11,440 states per sector.
Every spatial residue and energy-rate moment is bounded from the same
interacting source. The contact contribution is nonzero; interpreting
these finite responses as moving physical particles requires a
continuum and physical-clock construction.

The [proof](PROOF.md) gives the complete spatial matrices, exact current
equations and source assumptions. The [receipt](results.json) contains
every line and matrix entry, with no spectral cutoff. This is certified
source computation, not an executed native detector experiment.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/composite_spatial_response
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m model.composite_spatial_response --check
```

`--write` reconstructs the receipt in a fresh interpreter. The verifier
uses the retained integer witnesses and reruns their source, spectral,
rank and energy-moment checks; it requires no numerical eigensolver.
