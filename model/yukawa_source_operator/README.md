# Yukawa couplings from source records

A Yukawa coupling measures how a fermion propagation law changes when the
Higgs field changes. The native current construction fixes three allowed
interaction tensors. This package derives exact formulas that recover every
family coefficient from the corresponding Higgs variation. It also evaluates
complete, normalized three-point occupation records of the actual finite
source without measured masses or couplings. Applied to the existing declared
Higgs action, the extraction recovers its three explicit source-only trial
matrices and their singular values in its declared Higgs coordinate.

The [proof](PROOF.md) separates these two calculations. The source occupation
responses have exact values, while their identification with a Higgs field
and charged fermion propagation requires an additional source construction.
The [receipt](results.json) retains both calculations and their boundaries.

```sh
python3 -m model.yukawa_source_operator.certificate --check
python3 -m pytest -q model/yukawa_source_operator
```
