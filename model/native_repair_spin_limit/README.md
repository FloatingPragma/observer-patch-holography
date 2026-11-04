# Interacting spins under continuing native repair

The source repair acts at a fixed positive rate on every edge while the
two-flavor source approaches an interacting spin system. A local error
bound is independent of the total volume, retains complete local records,
and applies to explicit faithful source inputs with translation symmetry.
The estimate also controls a time interval that grows with the source
coupling. [The proof](PROOF.md) and [executable certificate](results.json)
specify the state, clock and record conventions.

The lattice spacing is fixed; a common relativistic limit with finite
massive and massless sectors and Einstein dynamics requires separate
derivations.

Run `python -m model.native_repair_spin_limit --check` and
`python -m pytest -q model/native_repair_spin_limit` from the repository
root. The small matrix tests retain every Fock state. The rational large
coupling examples certify analytic bounds and are not source executions.
