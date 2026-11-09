# Prepared scalar-odd energy records

A source-derived finite wavepacket can be prepared as an explicit
intervention, with the displaced source and its reference correlations
retained in idle memory. Full-space residual bounds certify its excited
energy window and a complete finite clock reads that window with original
probability greater than one quarter. The source Hamiltonian is unchanged;
this experiment differs from the original mixed-state pair detector.

The [proof](PROOF.md) states the preparation, errors and spectral boundary.
No isolated particle or GeV rest mass is identified. Run
`python -m model.particle_primary_spectral_resolution --write` and
`python -m pytest -W error model/particle_primary_spectral_resolution`.
