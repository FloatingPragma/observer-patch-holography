# ROW29 independent full review

The verdict is in `REVIEW.json`; the scientific argument and scope are in
`REVIEW.md`. `MANIFEST.json` freezes the exact producer candidates, all review
files, and source/import pins. This directory does not change an official row.

From the workspace root, reproduce using one numerical thread:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 openwave-tuning/next/CHAIN4/maxwell_full_review/review.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 openwave-tuning/next/CHAIN4/maxwell_full_review/size3_check.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 openwave-tuning/next/CHAIN4/maxwell_full_review/seal.py
```

Requirements are Python, NumPy, SciPy and SymPy; exact versions are recorded in
`IMPORTS.json`. The exact fixed-frame calculation uses bounded-integer GEMM
and arbitrary-precision recombination. Floating native residuals supplement
the existing outward interval/control certificate and do not replace it.

`SIZE3_ATTEMPT1.log` preserves an initial reviewer-harness unpacking error after
the frame certificate passed. The corrected complete run is `SIZE3_RUN.log`.
No producer data or chosen frame changed in the repair.

This source has c=1 under the declared noncompact-field h-time clock. The scalar
KG source's c=1/sqrt(48) and its clock conversion are separate; a pass here is
not a proof that those sources have been unified.
