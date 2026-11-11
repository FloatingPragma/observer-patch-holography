# A local field comparison on the interacting cube

The same prepared native Wilson–Maxwell source supports an explicit
bound on how a magnetic intervention changes a later electric detector.
Its physical electric field includes the charge reference and Coulomb
part. The complete finite experiment retains every preparation failure
and bounds the influence by `1.810490e-8` in the displayed setting.

The [proof](PROOF.md) derives an interacting commutator majorant on the
fixed cube, then transfers it through the actual ground preparation and
photon cutoff. Locality belongs to the reconstructed gauge-invariant
field/current algebra; a microscopic native or relativistic cone is not
proved.

Public replay from the repository root:

```sh
.venv/bin/python -B -m model.field_record_addenda.native_wilson_maxwell_locality --check
.venv/bin/python -m pytest -q model/field_record_addenda/native_wilson_maxwell_locality -W error
```
