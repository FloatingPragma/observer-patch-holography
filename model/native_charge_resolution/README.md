# Sharper native charge energy

A nested energy window on the unchanged native source retains a certified
unnormalized charge weight above 0.419706 while narrowing the energy
interval by more than 364 times. A narrower setting retains weight above
0.31477 in a source-energy interval only 8.382e-8 wide. All complements,
reference energy and preparation failures remain explicit.

The [proof](PROOF.md) supports an explicitly refined spectral selection.
It does not relocate the frozen maximum-residue eigenvalue or change any
existing charged coupling. A full spectral-rank or separation theorem is
required to identify that original eigenvalue with the refined selection.

From the repository root:

```
python3 -m model.native_charge_resolution --check
python3 -m pytest -q -W error model/native_charge_resolution
```

The [receipt](results.json) rebuilds the exact local roots and full-source
window bounds. Its [provenance](PROVENANCE.md) includes the immutable
research comparison and the unchanged native source.
