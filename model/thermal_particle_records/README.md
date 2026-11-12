# Finite thermal particle-count measurements

The same calibrated source supports complete number measurements for
its joint energy, charge, colour, helicity and momentum groups. The
instrument preserves colour coherence and includes every fermionic
occupation and all 1,653 capped boson states. Its finite thermal count
fractions keep empty events and many-particle states.

The [proof](PROOF.md) constructs the projectors and the actual boson
polarization transformation. The [receipt](results.json) records four
supplied temperatures and quantifies the initial many-excitation
complement of a one-particle detector. These are finite calibrated
measurements, without a cosmological abundance identification.

```sh
python3 -m pytest -q model/thermal_particle_records
python3 -m model.thermal_particle_records --write
python3 -m model.thermal_particle_records --check
```
