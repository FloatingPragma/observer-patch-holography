# A controlled native field readout

Actual signed native currents drive a local fixed-point Poisson solver.
Its error bound includes finite relaxation, moving charges and rounding.
The complete cap output is retained. A static Fourier theorem gives a
smoothed Coulomb limit on a unit periodic three-dimensional domain, with
an O(h) bound for localized internal ports and O(h²) for an explicitly
averaged comparison source. The [proof](PROOF.md) specifies the source,
clock and refinement maps; [RESULT.json](RESULT.json) records execution.

From the workspace root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
reverse-engineering-reality/.venv/bin/python openwave-tuning/next/CHAIN/field_limit/controlled_field.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
reverse-engineering-reality/.venv/bin/python openwave-tuning/next/CHAIN/field_limit/verify.py
```

The local arithmetic is an added observation service. The physical step
is a source-executed field with a specified distributed update clock and
a controlled dynamic refinement, including its reciprocal response.
