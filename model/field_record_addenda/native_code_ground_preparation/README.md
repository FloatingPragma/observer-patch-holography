# Recorded code admission before energy agreement

This research setting keeps the completed native Wilson-Maxwell source fixed.
It records whether four trial copies entered the existing code, then compares
their energies with a Cayley transform centered on the known free baseline.
All admission and filter failures remain in the law.

The admitted code gap, rather than the much smaller gap to native complement
states, controls this instrument. The pointer needs 4096 states instead of the
parent setting's approximately 3.15e17. This is an apparatus improvement,
not a changed Hamiltonian, an unknown-ground projector or an efficient passive
repair implementation. [PROOF.md](PROOF.md) states its clock and resources.

A smaller four-slot experiment uses two copies and one 64-state pointer.
Its ground approximation is less accurate, but it certifies a
positive original-attempt stationary T111 record with all failures retained.

Public replay from the repository root:

```sh
.venv/bin/python -B -m model.field_record_addenda.native_code_ground_preparation --check
.venv/bin/python -m pytest -q model/field_record_addenda/native_code_ground_preparation -W error
```
