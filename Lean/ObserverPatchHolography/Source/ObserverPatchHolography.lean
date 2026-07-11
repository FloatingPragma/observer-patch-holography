import ObserverPatchHolography.AbstractRewriting
import ObserverPatchHolography.Primitives
import ObserverPatchHolography.Rule90
import ObservableNormalForms
import ObserverPatchHolography.Bridges.ObservableNormalForms

/-!
# Observer-Patch Holography — Lean 4 umbrella root

Re-exports Jonathan Hill's concrete carrier/dynamics modules, the neutral
observation-determined normal-forms proof package, and the explicit bridge
between them.

**Status: preliminary skeleton rather than theorem-grade formalisation of
Proposition 4.2** from *Paradise as Fixed-Point Consensus*. The `Primitives`
module declares sorry-bearing signatures for
the OPH primitives (Records, Repair, Patch, Obs, Φ, gauge equivalence,
OPH-Confluence, OPH-Completeness) — these structurally depend on the
companion paper *Reality as a Consensus Protocol*.

See `../README.md` and `../PROOF_INDEX.md` for scope and completion tracking.
-/
