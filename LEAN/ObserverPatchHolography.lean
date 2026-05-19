import ObserverPatchHolography.AbstractRewriting
import ObserverPatchHolography.Primitives

/-!
# Observer-Patch Holography — Lean 4 library root

Re-exports the modules that make up the OPH formalisation effort.

**Current state is a preliminary skeleton, not a theorem-grade
formalisation of Proposition 4.2** from *Paradise as Fixed-Point
Consensus*. The `Primitives` module declares sorry-bearing signatures for
the OPH primitives (Records, Repair, Patch, Obs, Φ, gauge equivalence,
OPH-Confluence, OPH-Completeness) — these structurally depend on the
companion paper *Reality as a Consensus Protocol*.

See `README.md` and `PROOF_INDEX.md` for scope and completion tracking.
-/
