# Source ground and thermal spin references

The full fermionic source selects controlled spin references at fixed
lattice spacing. Its actual ground states obey the Heisenberg ground-state
inequality with an error independent of volume. Its faithful thermal states
have a controlled spin free-energy defect and an exact modular clock fixed
by the source Hamiltonian and temperature. All charge states are retained.
The [proof](PROOF.md) distinguishes finite-volume thermal convergence from
thermodynamic pressure and ground-state limits. A unique thermodynamic
phase, physical temperature calibration and relativistic continuum require
additional results.

Run `python -m model.spin_reference_limit --check` and
`python -m pytest -q -W error model/spin_reference_limit`.
