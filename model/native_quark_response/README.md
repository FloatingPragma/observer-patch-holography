# Native quark energy records

Literal quark removal from the native interacting local reference has a
compact exact spectral calculation. The source's 140-state quark module
maps into a verified 224-state image. Its full scalar and lepton couplings
give the complete energy measures for `u^c`, `d^c` and `Q`, including exact
zero residues. The unchanged eight-cell source admits corresponding energy
bins with original physical preparation errors and all failure records
retained. The computed `u^c`/`d^c` degeneracy belongs to the equal-coupling
source; it is not a physical mass prediction.

The [proof](PROOF.md) states the actual source, static color reference and
source clock. The [receipt](results.json) records every local root and
residue, the complete interacting bin law and retained reference cost four.
This calculation supplies no quark propagation or asymptotic colored state.

From the repository root:

```sh
python -m model.native_quark_response --check
python -W error -m unittest model.native_quark_response.test_response
```
