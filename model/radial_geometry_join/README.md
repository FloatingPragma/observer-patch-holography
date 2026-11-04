# Charged defects change the response geometry

The central interaction energy of the bounded radial fermion defect drives
the native geometry instrument. Its actual interacting ground and charged
states give different local parity signals, producing a certified geometry
contrast greater than 0.059. The forty-pulse quantum readout retains every
failure branch, and its measurement energy injection is derived explicitly.
The [proof](PROOF.md) and [receipt](results.json) bind the same finite field
states to both the particle and geometry calculations.

The controlled observation requires prepared inputs and blank memory;
autonomous Einstein dynamics requires additional regional and continuum
premises.

```sh
python -m model.radial_geometry_join --check
python -m pytest -q -W error model/radial_geometry_join
```
