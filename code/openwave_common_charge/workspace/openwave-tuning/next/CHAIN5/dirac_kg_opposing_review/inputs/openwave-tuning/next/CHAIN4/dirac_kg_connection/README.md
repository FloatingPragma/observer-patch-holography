# Routed Dirac Klein–Gordon connection

The [proof](PROOF.md) connects the routed Dirac mass and current source to its
explicit free parameter intervention, the same h-time clock, and the free
Maxwell field's speed. Three register-derived masses, all spin columns and
flavor padding are retained. Source-specific evidence does not identify the
scalar walk or the full90 Wilson source with this construction.

Reproduce from the workspace root with Python, NumPy, SciPy and SymPy:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 openwave-tuning/next/CHAIN4/dirac_kg_connection/geometry_checks.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 openwave-tuning/next/CHAIN4/dirac_kg_connection/native_tomography.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 openwave-tuning/next/CHAIN4/dirac_kg_connection/spatial_check.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 openwave-tuning/next/CHAIN4/dirac_kg_connection/controls.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 openwave-tuning/next/CHAIN4/dirac_kg_connection/seal.py
```

The ordinary finite native runs are executed. The L=2^18 resource contract and
full tensor-field holding iterator are finite, unexecuted constructions.
`GEOMETRY_ATTEMPT1.log` records a corrected review-harness permutation-transpose
error. Producer files and selected input data are unchanged.
