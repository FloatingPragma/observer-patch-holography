# Light nuclear spectra

A common native quark action gives complete low-energy bands for three-
and four-baryon boxes. Exact local representation theory reduces spaces
of 572572 and 176729280 colour singlets to 64 and 256 low-band coordinates,
with rigorous bounds on every omitted coordinate. One range-three spin
coefficient calibrates a triton-like rest energy. The held-out helion and
alpha masses fail, so the finite boxes do not establish the physical
nuclear binding curve or halos. [Proof](PROOF.md), [results](results.json).

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m model.openwave_light_nuclear_curve --check
python -m pytest -q model/openwave_light_nuclear_curve
```
