# Finite quantum readout calibration

A complete finite quantum instrument is fixed by basis-state survival and
one coherent survival test for each nontrivial outcome block. These finite
joint tests determine its effects and Lüders state update; approximate tests
give an explicit bound valid with arbitrary untouched references and adaptive
continuations. The construction replaces an all-state invariance requirement
by a finite calibration contract. Quantum probability semantics and trusted
preparation/readback remain assumptions.

[Proof and error bound](PROOF.md) · [Exact controls](results.json)

```sh
python3 -m model.finite_instrument_calibration.certificate --check
python3 -m pytest -q model/finite_instrument_calibration
```
