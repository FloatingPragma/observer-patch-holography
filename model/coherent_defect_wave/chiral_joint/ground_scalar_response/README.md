# Scalar response of the quantum ground

The native scalar excitation number and the total Yukawa interaction
have one positive inelastic transition on the actual local ground,
near 23.4612 in the declared clock units. Exact identities reduce their
response to two energy levels. The [proof](PROOF.md) and
[receipt](results.json) retain the spin-degenerate elastic line,
every zero-weight energy outcome and both probe failures.

These are neutral spin-zero probes on the selected local ground
reference. The transition frequency has no relativistic rest-mass
identification, and spatial propagation requires a separate calculation.

Run from the repository root with the pinned `model/requirements.txt`:

```sh
python -m model.coherent_defect_wave.chiral_joint.ground_scalar_response --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/ground_scalar_response
```

`--write` regenerates the receipt using exact operator identities and
rational root isolation.
