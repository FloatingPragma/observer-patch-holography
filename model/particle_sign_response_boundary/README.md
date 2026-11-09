# Electric response of the native Wilson sign

An exact native gauge-coordinate family has complete squared Wilson
spectrum `(delta²:16, 1:28, 3+2delta:16)` and a unit-electric sign
derivative of norm `(1-delta)/(4delta)`. The code retains all sixty
Dirac/M15 columns and all twelve coordinate-symbol readout labels.

Small original weight in this block can suppress its low-mode count while
leaving an arbitrarily large inverse-gap derivative cost. This isolates
why the [local spectral theorem](../particle_local_wilson_density/PROOF.md)
does not by itself prove the interacting local Duhamel estimate. It is
not a failure theorem for the configured mixed state or a particle-mass
calculation. The [proof](PROOF.md) distinguishes coordinate controls,
normalizable packets, finite approximations and sufficient derivative
hypotheses.

```sh
python -m model.particle_sign_response_boundary --write
python -m model.particle_sign_response_boundary --check
python -W error -m pytest -q model/particle_sign_response_boundary
```
