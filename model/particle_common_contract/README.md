# One source, state, reference and clock

The [construction contract](PROOF.md) fixes an open Wilson gauge-field
program, its Gauss-constrained polynomial preparation, a retained conjugate
reference and its Hamiltonian clock. Its finite sequence has no fitted
mass inputs. Executable checks reject mixed preparations, reference clocks,
regulators and unclassified comparison results. This resolves the choice
of construction. Its [physical preparation and boundary detector](../particle_physical_boundary_state/PROOF.md)
retain every charge outcome, with explicit nonzero Wilson-loop sectors.
The [locality consumer](../particle_primary_locality/PROOF.md) bounds actual
complete boundary/reference records under a remote plaquette intervention
on this same state and clock.
Joint physical refinement and particle identification are separate output
requirements.

Run `python -m model.particle_common_contract --write`, then
`python -m model.particle_common_contract --check`. The exporter constructs
finite specifications and verified bounds, without executing the enormous
full spatial matrix or preparation tape.
