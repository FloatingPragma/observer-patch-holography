# Scalar overlap with localized defects

A local scalar interaction couples differently to modes with different spatial
profiles. The complete four-mode band of a compact defect pair gives two doubled
response strengths with ratio about 16,385 for a specified local scalar patch.
A uniform scalar gives equal strengths. Every site insertion and its leakage
into the full exterior are retained, without family coefficients or measured
masses. The [proof](PROOF.md) and [calculation](results.json) give the finite
source result. A physical Yukawa interpretation requires a dynamical scalar
mode, its kinetic residue and a propagating field map.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m model.defect_mode_overlap --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q model/defect_mode_overlap
```
