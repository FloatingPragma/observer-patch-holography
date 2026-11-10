# Current records of a native Higgs excitation

Run from the repository root:

```sh
python -m model.native_higgs_current_response --check
pytest -q model/native_higgs_current_response
```

The [proof](PROOF.md) joins the actual quantum-Higgs preparation to
the original dressed weak-current detector on the unchanged full
native Hamiltonian. Its elementary and composite preparations both
have a strictly negative current mean in the same admitted energy
bin, with every preparation and detector failure retained.

Charge conjugation changes the detector W to Wbar. The implementation
calculates that transformation and the full charged blocks; it does
not copy the opposite-charge current value. The receipt contains
source-clock bounds, original-attempt probabilities and the explicit
remaining field/particle interpretation requirements.
