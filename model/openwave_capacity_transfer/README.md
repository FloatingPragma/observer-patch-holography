# Constrained capacity, matter and radiation transfer

This executable attaches the retained canonical capacity scalar to **linear
Einstein gravity, conserved ideal dust and ideal radiation**. It evolves the
metric and both fluids, verifies the Hamiltonian constraint and spatial trace
as exact symbolic identities, and computes normalized linear transfer functions.

Run `python -m model.openwave_capacity_transfer --write`, then run the same
module with `--check` for the original cold comparison. Run its tests with
`python -m pytest model/openwave_capacity_transfer -q`.

For the inherited background inputs, the present comoving matter transfers at
`k/H0 = 0.1, 1, 10, 100` are approximately `0.97585, 0.97001, 0.91494, 0.90822`
times those of the same-density cosmological-constant control. These are
**conditional numerical model comparisons**, not measured matches, fitted
primordial amplitudes or predictions of a native capacity-to-scalar mechanism.

The exact constraint propagates homogeneously without resetting or projecting
the numerical solution. Tolerance refinement, earlier initialization, a
deliberately nonzero constraint and the independent Einstein–de Sitter analytic
solution test distinct failure modes. Floating ODE errors are not rigorous
interval bounds. Free-streaming neutrinos, recombination and the full Boltzmann
hierarchy remain absent. See [PROOF.md](PROOF.md) for equations and scope.
