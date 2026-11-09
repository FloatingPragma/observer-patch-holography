# Energy flow on one physical source

The selected Wilson/Gauss program has an exact local energy partition and
interaction currents. Every regional energy change equals its complete
boundary flux on the same Hamiltonian clock. The construction retains
all source modes, boundary charges and the idle quantum reference, and
provides finite current-readout error bounds. A separate exact quantum-link
control exhibits nonzero energy transfer.

The same state's supported regional marginal also satisfies the finite
entropy first law under a complete controlled intervention. This marginal
entropy differs from the entropy of the retained coin and reference.
Neither result identifies proper geometry or a Lorentz stress tensor;
the [proof](PROOF.md) states the missing Einstein inputs.

```sh
python3 -m pytest -q model/particle_primary_stress -W error
python3 -m model.particle_primary_stress --write
python3 -m model.particle_primary_stress --check
```
