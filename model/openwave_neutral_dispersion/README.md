# Why the neutral band does not move

The periodic two-cell source has three exact symmetries combining matter
parity with inversion of its quantum links. They make every odd-sector
energy identical at all eight lattice momenta. A dressed local neutral
excitation therefore cannot move to another atomic cell at any time;
the virtual-link/profile complement is retained. This is an exact
finite-source obstruction, not a small numerical dispersion.

The [proof](PROOF.md) checks the full native link and current identities.
Run `python3 -m model.openwave_neutral_dispersion --check` and
`python3 -m pytest -q model/openwave_neutral_dispersion` from the repository
root. A larger or explicitly changed source needs its own propagation proof.
