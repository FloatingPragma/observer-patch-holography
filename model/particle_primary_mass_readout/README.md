# From source transitions to mass readouts

Cooling controls the difference between an original energy transfer and
the final excitation energy of the same selected Wilson source. The
[proof](PROOF.md) gives a quantitative error without dividing by a
measurement success probability. The [receipt](results.json) binds actual
source levels and checks a mixed-state control where a transition energy
would give a false rest-energy assignment. Rest masses require the same
source's momentum and shell attachment. Declared SI comparison units are
supported explicitly; a unit choice does not provide a missing spectrum.

```sh
python -m model.particle_primary_mass_readout --check
python -W error -m pytest -q model/particle_primary_mass_readout
```
