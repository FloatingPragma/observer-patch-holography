# Local Wilson spectral density

The same translation-invariant gauge-coordinate state and Hamiltonian give
an all-time, volume-independent bound on the expected local weight of Wilson
modes below absolute eigenvalue1/2. The proof turns the HJL positive
decomposition into a trace estimate and instantiates its density premise
with the exact [local gauge preparation](../particle_local_gauge_state/).

For its requested roughness allowance eta, a region X satisfies
`<Tr(P_X Pi_low P_X)> <= 2 sqrt(5 eta)|X|`. The package also gives a
local Hilbert–Schmidt sign-polynomial bound and a complete finite binary
probe. Every native internal and Dirac mode is retained.

This is a local spectral-weight result, not a global Wilson gap, spatial
localization theorem, finite-duration interacting comparison or particle
mass calculation. See [the proof and remaining interface](PROOF.md).

```sh
python -m model.particle_local_wilson_density --write
python -m model.particle_local_wilson_density --check
python -W error -m pytest -q model/particle_local_wilson_density
```
