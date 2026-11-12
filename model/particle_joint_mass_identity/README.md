# Joint particle and mass records

This package measures commuting charge, colour Casimir, helicity,
momentum and energy on the actual finite quadratic source of
[conditional_particle_spectrum](../conditional_particle_spectrum/README.md).
The probe has no species or family label. Its joint outputs distinguish
the source's conditional mass types while preserving colour degeneracy,
antiparticle modes and the quantum reference. The complete instrument
also retains vacuum and the entire many-excitation complement.

A separate exact sixteen-state calculation refines the compiled free
lepton mass instrument by spin. Forgetting spin reproduces every original
mass/charge/reservoir record. Its GeV intervals retain the one electron
anchor; the calibrated full-action coefficients retain their different
input conventions.

```
VECLIB_MAXIMUM_THREADS=1 python3 -m pytest -q model/particle_joint_mass_identity
VECLIB_MAXIMUM_THREADS=1 python3 -m model.particle_joint_mass_identity --write
VECLIB_MAXIMUM_THREADS=1 python3 -m model.particle_joint_mass_identity --check
```

[PROOF.md](PROOF.md) gives the full-domain instrument identity and the
finite-resolution boundary. Exact arithmetic resolves the optional
neutrino coefficients at momentum 1 GeV; ordinary binary64 energy
diagonalization and the wider detector bins do not. These are finite
source identities, not a derived physical particle dictionary or
independent predictions of the calibration inputs.
