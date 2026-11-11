# Scalar pulses and native family responses

The same interacting one-cell source gives complete scalar excitation records and
different charged-scalar pulse responses for its quark and lepton sectors,
including both outcomes of an exact finite pulse.
All primitive interaction coefficients are one. The scalar normalization
comes from the actual thirty-mode ground and every accessible energy,
while the contact tensors retain all three quark colors. No measured
particle inputs enter the calculation.

The [proof](PROOF.md) distinguishes these finite responses from physical
Yukawa couplings and describes the neutral composite and a decoupled sterile
extension. The [receipt](results.json) contains the full scalar energy
records, all nonzero contact components and the complete neutral records.

```sh
python -m model.defect_family_vertices --write
python -m model.defect_family_vertices --check
python -m unittest model.defect_family_vertices.test_vertices
```
