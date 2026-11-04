# Flavor from current-dependent memory

The declared native current gates the complete repair record differently in
the up, down and charged-lepton sectors. Its stationary memory produces
different charged spectra and nontrivial real quark mixing. The calculation
retains all outcomes, all three family directions and all 45 internal rows.

The [proof](PROOF.md) classifies every family basis change preserving the same
currents and kinetic form. Such changes preserve the spectra and mixing.
All four source duty policies have exactly zero charged-current CP invariant,
including policies with an exact zero charged direction. The readable angle,
reset and feedback laws determine these dimensionless results; identifying
one such implementation with physical particle observations requires a
same-source matching theorem.

Every finite feedback history has an exact auxiliary antiunitary symmetry.
Its principal Hermitian generator has a real symmetric mass block and a
separate imaginary antisymmetric diagonal block. The complete unread
finite-memory channel has strictly positive Choi purity loss on the stated
positive-gate domain for all four duties, including arbitrary independent
initial memory at two or more steps. Exact lower bounds distinguish it from
every unitary channel. The full twelve-memory channel has no matrix
logarithm. These statements do not identify the auxiliary symmetry with
physical CP, discard the diagonal correction, or promote the stationary
coefficient catalogue to the entire correlated process.

Run the [receipt](results.json) calculation from the repository root:

```sh
python -m model.native_flavor_orbit --write
python -m model.native_flavor_orbit --check
python -m pytest -q model/native_flavor_orbit
```
