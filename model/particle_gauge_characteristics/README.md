# Gauge and matter wave speeds

The native magnetic loop action and electric trace metric give the same
quadratic wave speed to all twelve gauge directions. On the declared
clock, this speed matches the Dirac principal symbol exactly when
`electric_scale * magnetic_scale = 15`. This uses no particle data.
The [proof](PROOF.md) retains gauge-null modes and gives the actual
lattice dispersion and its error. Identifying the tangent waves with
physical particles requires the interacting prepared-state response.

Run `python -m model.particle_gauge_characteristics` for the certificate,
or `python -m pytest -q -W error model/particle_gauge_characteristics`.
