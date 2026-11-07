# Electric motion and a chiral projector

When a Euclidean overlap projector varies with a gauge coordinate, electric
motion needs a compatible connection. Compressing the declared parent
kinetic energy also supplies a positive geometric term. This construction
checks both terms on the source projector and retains the complementary
space in a finite penalty Hamiltonian. It does not select the microscopic
kinetic law or prove a global chiral gauge theory.

The [proof](PROOF.md) states the local identities and finite-parent bound.

```sh
python -m model.interacting_extension.chiral_electric.certificate --check
python -m pytest -q model/interacting_extension/chiral_electric
```

The `covariant_kinetic` module accepts `--out DIRECTORY` for complete parent
and compressed matrices plus numerical diagnostics.
