# Particle production during a declared expansion

The calibrated fermion source evolves under a supplied scale factor from
one to four. The executable keeps all 180 modes, computes all twelve species,
and retains the complete Gaussian states and correlated number readouts.
Exact pair phases and interval arithmetic connect the calculation to the
literal source Hamiltonian. The clock, expansion law and calibrated masses
are inputs; physical cosmological attachment is separate.

Run `python3 -m unittest model.expanding_quadratic_particles.test_production`,
then `python3 -m model.expanding_quadratic_particles --write` and
`python3 -m model.expanding_quadratic_particles --check`.
The [proof](PROOF.md) states the continuous-ramp error bound and source scope;
[results.json](results.json) contains every protocol, pair state and number law.
