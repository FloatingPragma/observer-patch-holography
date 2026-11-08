# Preparing a field at controlled positive selection

The source has an explicit family of decreasing selection fields with finite
ground-preparation and local measurement error bounds. Its order and quadratic
bond coefficients concentrate while transverse field histories keep bounded
scaled energy. [The construction](PROOF.md) supplies two canonical pairs and
their full internal circle symmetry on that source clock, with complete local
record protocols at each finite scale. It is distinct from the field defined
by a different order of limits; no sharp dispersion or physical particle
assignment follows from this result.

```sh
python -W error -m pytest -q model/spin_joint_field_limit
python -W error -m model.spin_joint_field_limit --write
python -W error -m model.spin_joint_field_limit --check
```
