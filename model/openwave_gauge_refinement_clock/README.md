# Actual spatial Haar blocking

Subdividing a cubic gauge source adds links and plaquettes while retaining
all Haar representations. Exact integration produces the coarse Wilson
interaction at fourth order, with coefficient `kappa^4/64` and a uniform
sixth-order remainder. Two declared steps have certified positive
coefficients. The [proof](PROOF.md) keeps the source clock and the missing
full-Hamiltonian flow explicit.

Run `python -m model.openwave_gauge_refinement_clock --write`, then
`python -m model.openwave_gauge_refinement_clock --check`.
