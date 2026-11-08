# Source energy transport and a measured boost candidate

The selected eight-spin source supplies three weighted energy operators
B_i=sum_x r_i(x)e_x and their exact fluxes P_i=i[H_B,B_i]. The coordinates
r_i=+-1/2 are the source cube coordinates. These operators preserve the
same internal charge used in the finite energy readout. They supply an
operational candidate for a boost; physical momentum and Lorentz covariance
require additional source estimates.

The construction reuses the finite spin energy Ward identity and the Lean
coarse-energy and forced-balance theorems. It supplies three local flux
densities, antisymmetric stress columns and explicit source forces. Every
cell equation and its coordinate block transport is checked on the full
256-state source algebra. The created charge channels have nonzero energy
flux acceleration under this same clock.

Coarse energy totals alone lose information needed by the boost candidate.
Retaining each block's energy dipole gives exact transport of both the
weighted energy and its flux. This identifies a concrete extra observable
for a refinement-compatible spatial energy readout.

A complete native experiment prepares the same ground from blank and
controls its weighted-energy flow. Its minus-pointer probability is
certified above 5.12e-9. The source ground is therefore not invariant under
this finite candidate. Exact source algebra also gives a nonzero boost
commutator defect and a nonconserved total energy flux. Positive changes
of the candidate speed cannot remove the ground defect.

The [proof](PROOF.md) converts those defects into quantitative complete
history bounds and states the source estimates a continuum boost would
need. No physical metric, particle dispersion or symmetric Lorentz stress
is assigned by definition.

Run `python -m model.spin_boost_stress --check
model/spin_boost_stress/results.json` and
`python -W error -m pytest -q model/spin_boost_stress`.
