# Quark-type charge transport on the native interacting cube

This research calculation extends the literal compensated first-hop
experiment to all three quark-current families. It reuses the full native
Hamiltonian, its physical preparation and fixed reference. The exact
quark-spin-singlet compression gives the route effect `P0/128` for every
localized ground-doublet input, including its retained environment.

Public replay from the repository root:

```sh
.venv/bin/python -B -m model.field_record_addenda.native_quark_transport --check
.venv/bin/python -m pytest -q model/field_record_addenda/native_quark_transport -W error
```
