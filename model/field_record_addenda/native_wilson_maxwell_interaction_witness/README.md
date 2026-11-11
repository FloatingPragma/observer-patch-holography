# A local record that detects the original coupling

The same nonstationary native-code preparation receives an electric
kick, evolves under the actual Wilson–Maxwell source, and undergoes a
complete bounded current readout. A source-only analytic hopping bound
proves a strictly positive setting-conditioned contrast. The contrast
vanishes in the zero-coupling comparison. All native admissions,
failures, field cutoffs, setting bits and detector outcomes are retained.

[PROOF.md](PROOF.md) gives the uniform moment/remainder estimates,
literal finite-cap transfer and fixed-spacing local-record limit.
This is research evidence, with no interacting-ground, scattering,
Lorentz/ultraviolet or particle-identity promotion.

The sufficient box and photon capacity are analytic bounds; their matrices
are not numerically enumerated.

Public replay from the repository root:

```sh
.venv/bin/python -B -m model.field_record_addenda.native_wilson_maxwell_interaction_witness --check
.venv/bin/python -m pytest -q model/field_record_addenda/native_wilson_maxwell_interaction_witness -W error
```
