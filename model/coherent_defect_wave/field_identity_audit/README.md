# Charges carried by a flux object

A literal native occupation cloud carries the fifteen familiar fermion
weights, and a specified transport can carry that cloud with a protected
flux object. The [source calculation](PROOF.md) checks their actual current
intertwiner, rather than assigning charges to five face labels. The native
idle update fixes this cloud, so its weak interaction contributes no mass
correction. A chiral interacting field and measured masses require dynamics
beyond that invariant cloud. The [receipt](results.json) records exact source
controls and the complete finite transport.

```sh
python -m model.coherent_defect_wave.field_identity_audit --check
python -m pytest -q model/coherent_defect_wave/field_identity_audit
```
