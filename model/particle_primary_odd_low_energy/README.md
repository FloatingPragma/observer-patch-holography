# Energy filters for odd source excitations

A complete coherent filter selects an odd excitation's source-energy band
with more than 3/5 original probability on the smallest member, including
clock errors. Every other selector outcome is retained. Summing those
outcomes preserves the original energy law. The same source proves a
necessary growth rate for phase dressing to remove its escaping
high-energy component. A positive fixed-energy residue, physical dispersion
and a rest mass require their own source spectral bounds.

See the [proof](PROOF.md) and [receipt](results.json).

```sh
python -m model.particle_primary_odd_low_energy --check
python -m pytest -q -W error model/particle_primary_odd_low_energy
```
