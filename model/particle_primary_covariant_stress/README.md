# Energy with the source's rotation symmetry

Split each interior link's electric energy equally between its ends and
each plaquette's energy equally between its four corners. The resulting
local energy rotates with the physical cube, conserves the unchanged
Hamiltonian and admits complete energy/current measurements. Its difference
from the owner allocation is an exact boundary flow. Conjugate source
states have the same energy and current records. Proper geometry and
refinement-uniform stress estimates require separate proofs.

See [the derivation](PROOF.md) and [source-bound receipt](results.json).
Run `python -m model.particle_primary_covariant_stress --check` and
`python -m pytest -q -W error model/particle_primary_covariant_stress`.
