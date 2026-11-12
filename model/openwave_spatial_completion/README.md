# Full spatial domain and trap release

The source, exact argument, complete carrier and physical boundaries are stated
in `PROOF.md`. `results.json` pins all computational inputs and records every
finite witness used by the certificate.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_spatial_completion --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_spatial_completion --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q model/openwave_spatial_completion
```

Declared source choices and empirical energy anchors retain their ancestry.
A successful finite or infinite-source theorem does not by itself identify
an observed electron or a physical electromagnetic calibration.
