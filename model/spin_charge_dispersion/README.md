# Bounded charged excitations with measurable energy and momentum

The ordered spin source supplies two bounded preparation branches with exact
opposite internal charges and a uniform positive probability in a joint
energy-momentum band. Their energies scale with the momentum radius on the
same field clock. The construction reuses the source's canonical fields and
spectral moment bounds; all preparation and detector complements are retained.
The [proof](PROOF.md) and exact finite source certificate identify a measurable
dispersive band. A sharp relativistic particle shell is unproved.

```sh
python -W error -m pytest -q model/spin_charge_dispersion
python -W error -m model.spin_charge_dispersion --write
python -W error -m model.spin_charge_dispersion --check
```
