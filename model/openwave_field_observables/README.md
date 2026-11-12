# Finite electromagnetic waves and forces

The declared cochain field produces two transverse wave polarizations and a
Coulomb-shaped force from one incidence stencil. Exact finite-occupation
Heisenberg identities reduce propagation to field-coordinate equations,
avoiding a full Fock-space simulation. Independent real-space evolution and
Green-function calculations test the same operators. A separate complete
native mass-channel calculation verifies moving massive wave records and
quantifies the finite-spacing error under boosts. A Dirichlet Green-field
consumer produces negative bound levels with exact reflection-sector
reduction and explicit spatial refinement.

The source geometry, field normalizations and mass coefficients are declared
inputs; matching these observables to physical particles requires the
corresponding source and laboratory identifications.

[Equations and scope](PROOF.md) · [Executed results](results.json)

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_field_observables --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_field_observables --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q model/openwave_field_observables
```
