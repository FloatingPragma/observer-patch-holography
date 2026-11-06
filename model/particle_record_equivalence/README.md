# Repair records and particle readout

Exact recorded recovery preserves the logical evolution of the native
matching code. The specified family memory also admits a four-address
encoding with exactly the same complete quantum record histories as its
twelve-address version, for compatible preparations. This supplies a concrete
finite capacity comparison. Changing pulse duty changes the action; the
finite observer readout does not provide a fitted mass correction.

[The proof](PROOF.md) states the identities and their scope.
[The receipt](results.json) binds the native source tables and numerical
readout checks. Physical spectrum identity requires a selected preparation,
action and laboratory observation map.

From the repository root:

```sh
python3 -m model.particle_record_equivalence.certificate --check
python3 -m pytest -q model/particle_record_equivalence
```
