# Preparing a field at controlled positive selection

The source has an explicit family of decreasing selection fields with finite
ground-preparation and local measurement error bounds. Its order and quadratic
bond coefficients concentrate while transverse field histories keep bounded
scaled energy. [The construction](PROOF.md) supplies two canonical pairs and
their full internal circle symmetry on that source clock, with complete local
record protocols at each finite scale. It is distinct from the field defined
by a different order of limits; no sharp dispersion or physical particle
assignment follows from this result.

[Source Ward identities](WARD_DYNAMICS.md) reduce both force errors to
static source scores and two inverse-energy susceptibilities. Their
filtered values give an explicit original-probability energy-window bound.
The complete finite source exhibits both exact force closure and a
polarization mismatch beside a sharp energy line; no vanishing continuum
force error follows from those controls.
For the volume-wide [uniform Ward mode](UNIFORM_MODE.md), the source first
moment and a computable finite-box order floor make its inverse-weighted
energy centre and energy-square variance converge to zero. Its original
preparation-window probability has an explicit volume-dependent lower bound.
Complete finite cube spectra retain all charge, reflection and energy
records; they cannot replace the growing-volume field or establish a
nonzero-momentum particle branch.

```sh
python -W error -m pytest -q model/spin_joint_field_limit
python -W error -m model.spin_joint_field_limit --write
python -W error -m model.spin_joint_field_limit --check
```
