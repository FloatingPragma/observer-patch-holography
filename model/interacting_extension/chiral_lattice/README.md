# Propagation and quantum flux

Two finite controls expose measurable propagation and gauge interaction.
The complete free spatial overlap branch approaches a relativistic cone.
A separate abelian ring evolves quantum electric flux, hopping and a loop
interaction, with every cutoff boundary and recorded outcome retained.
Its finite-time record error has an analytic cutoff bound. These controls
do not establish a common source-selected interacting chiral continuum.

The [proof](PROOF.md) gives the assumptions and recorded laws. Reproduce:

```sh
python -m model.interacting_extension.chiral_lattice.certificate --check
python -m pytest -q model/interacting_extension/chiral_lattice
```

The `light_branch`, `rotor_response` and `independent_audit` module CLIs
accept `--out DIRECTORY` for complete numerical diagnostics and arrays.
