# Checking the charged ground

An independent finite Taylor calculation and exact rational inertia counts
verify all thirteen compression intervals used to bound the actual charged
ground. The [proof](PROOF.md) checks the complete photon complement,
strict field dressing and retained-failure preparation law on the same
finite source. Its boundary reference remains fixed; the energy bounds
do not establish a relativistic mass shell.

Run from the repository root:

```
python -m model.affine_gauss_spectrum_validation --check
python -m pytest -q model/affine_gauss_spectrum_validation
```
