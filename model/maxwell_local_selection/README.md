# Local Maxwell selection

Gauge invariance, strict face-sharing locality and the committed rotation
symmetry force every quadratic seam action, up to an additive constant,
to be a multiple of `||C A||^2`. Exact elimination classifies all ninety
allowed Hessian coefficients; no action menu or numerical fit is used.

Among real linear autonomous laws, Faraday evolution, a supplied positive
electric/magnetic energy and the specified time reversal force Ampere
evolution. Electric and magnetic normalizations, physical time and native
implementation remain inputs.
The proof distinguishes this finite result from a physical Maxwell or
photon derivation and supplies counterexamples when assumptions are dropped.

[Proof and assumptions](PROOF.md) · [Exact receipt](results.json)

```sh
python3 -m model.maxwell_local_selection.certificate --check
python3 -m pytest -q model/maxwell_local_selection
```
