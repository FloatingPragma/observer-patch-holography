# Energy carried by a native quark hop

The selected first-hop response has the original quark spectrum shifted
by one unit of flux energy. Its incoming and outgoing atomic energies
are nonresonant. On the same full native Hamiltonian, a complete experiment
records a positive arrival and then a resolved energy after one source
time unit. The arrival measurement changes the state before the energy
readout; this is not a free-particle dispersion or a confinement theorem.

See [the derivation](PROOF.md) and [the exact certificate](results.json).
From the workspace root, use the RER Python environment to run
`spatial.py` for a checked replay or `--write` to regenerate it. Run
`pytest -q -Werror` on this folder for the source and normalization
controls.

Public replay from the repository root:

```sh
.venv/bin/python -B -m model.field_record_addenda.native_quark_spatial_spectrum --check
.venv/bin/python -m pytest -q model/field_record_addenda/native_quark_spatial_spectrum -W error
```
