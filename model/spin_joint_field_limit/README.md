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
The [two-site ground marginal](QUANTUM_FORCE_SLACK.md) and the source
variational energy bound imply `liminf(kappa_B-m_B^2)>=61/14400` on the
joint family. Thus replacing each longitudinal spin by its staggered mean
cannot close the unfiltered microscopic force equations. A corrected
infrared dispersion requires separate spectral control.
For the volume-wide [uniform Ward mode](UNIFORM_MODE.md), the source first
moment and a computable finite-box order floor make its inverse-weighted
energy centre and energy-square variance converge to zero. Its original
preparation-window probability has an explicit volume-dependent lower bound.
Complete finite cube spectra retain all charge, reflection and energy
records; they cannot replace the growing-volume field or establish a
nonzero-momentum particle branch.

The [bounded charged preparation](ARCTAN_CHARGED_BANDS.md) has two conjugate
records, each with original probability at least `1/387200` jointly in
momentum `[r/2,r]` and energy `[3r/387200,387200r]`, uniformly at small r.
It uses this positive-field ground and clock, including a finite-source
commutator error and every complementary record. The complete dimer
control retains nine preparation/energy labels. These broad response bands
supply no sharp dispersion or physical mass.

The [directed preparation](DIRECTED_CHARGED_BANDS.md) resolves each charge
into arbitrarily narrow opposite momentum caps, retaining original probability
at least `1/774400` in each charge and direction. Its spatial translation
phase has a controlled error; its energy interval remains broad. A literal
six-site periodic source verifies the charge/direction readout with all
13 outcomes and certified positive energy variance. This finite control is
not a continuum dispersion law or an identification with the odd CAR probes.

```sh
python -W error -m pytest -q model/spin_joint_field_limit
python -W error -m model.spin_joint_field_limit --write
python -W error -m model.spin_joint_field_limit --check
```
