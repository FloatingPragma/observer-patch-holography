# Propagation under spatial refinement

The declared gauge-coordinate field has a causal influence bound that
is uniform as its assigned lattice spacing decreases. It retains the
quantum electric links, every fermion occupation and complete local
apparatus records. The [proof](PROOF.md) gives exact spacing certificates
and a sharper finite-range Wilson alternative, without identifying the
bound with a measured light speed.

State and observable convergence, Lorentz symmetry and physical particle
identification require separate arguments. Run
`python -m model.particle_physical_causality` for the certificate and
`python -W error -m pytest -q model/particle_physical_causality` for tests.
