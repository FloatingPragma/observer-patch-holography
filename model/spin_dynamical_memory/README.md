# Dynamical response from the spin source

The selected spin source determines an exact response equation with a memory
term that retains every complementary excitation. A complete 256-state
certificate measures that term and a noncentral time history using the source's
own ground and clock. Both conserved spin-charge channels have nonzero memory;
the full response retains memory after finite field removal. Bounded resolvent
filters transfer the response equation
to the same canonical continuum histories without an unbounded-moment shortcut.
A [coherent charge pointer](CHARGE_READOUT.md) preserves each charge class
and bounds aliases in the same continuum field.
The [proof](PROOF.md) distinguishes these source results from Lorentz covariance
and surviving continuum scattering, whose dynamical limits are unproved.

A [bounded reducing defect](RESOLVENT_CLOSURE.md) tests whether a specified
response space is preserved by the actual positive source clock. Exact dimer
controls retain all preparation and filter outcomes and distinguish a
nonreducing prepared vector from a reducing space with two energies. Two
canonical pairs and an internal U1 do not determine the source wave equation:
a local two-speed field control satisfies these inputs with a hidden force.
The criterion does not certify a vanishing continuum defect or implement an
arbitrary Gram projection as a detector.

A [source spectral hierarchy](SPECTRAL_BRANCH.md) approximates bounded
energy scores using actual resolvent moments, with finite-resolution and
boundary-weight errors. Time-translated local creators realize those moments
with explicit cutoff and control costs; Gram-certified coherent combinations
retain their original preparation probability. Exact dimer energy and
site-symmetry labels illustrate the construction. Finite Ritz nodes do not
certify a continuum pole, and no uniform growing Gram floor is inferred.

```sh
python -W error -m pytest -q model/spin_dynamical_memory
python -W error -m model.spin_dynamical_memory --write
python -W error -m model.spin_dynamical_memory --check
```
