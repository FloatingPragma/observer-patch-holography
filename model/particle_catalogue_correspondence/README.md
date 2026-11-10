# Centre classes and particle currents

A six-class defect register records colour triality and weak centre parity.
Its complete exterior-weight comparison is reproducible, including the
fifteen-state current representation used by the interacting quantum source.
The same register loses information required for a particle identification.

The exact witnesses distinguish a colour octet from a singlet with the same
centre, a bosonic diquark from a fermionic antiquark with the same full gauge
representation, and opposite hypercharges in the self-inverse class three.
These restrictions apply to a proposed transfer between constructions; they
do not prohibit declaring the larger physical representation.

Run from the repository root:

```sh
python -m model.particle_catalogue_correspondence --check
python -m pytest model/particle_catalogue_correspondence
```

[PROOF.md](PROOF.md) states the calculation and the required transfer.
[results.json](results.json) contains every weight, the exact witnesses and
the hashes of their imported operators. It supplies no measured particle
mass or physical species identification.
