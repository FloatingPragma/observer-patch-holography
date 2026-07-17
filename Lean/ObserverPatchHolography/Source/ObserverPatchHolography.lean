import ObserverPatchHolography.AbstractRewriting
import ObserverPatchHolography.Primitives
import ObserverPatchHolography.Rule90
import ObserverPatchHolography.BoundaryFiber
import ObserverPatchHolography.BridgeEquivalence
import ObserverPatchHolography.CapacityFixedPoint
import ObserverPatchHolography.SeedPi
import ObservableNormalForms
import ObserverPatchHolography.Bridges.ObservableNormalForms
import EventAlgebra

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

The `BridgeEquivalence`, `CapacityFixedPoint`, and `SeedPi` modules carry
the Part-A coupling-algebra layer: the bridge count/tick equivalence, the
capacity fixed-point uniqueness schema, and the CAP-P seed statement. They
formalise the algebraic layer only; the physical identities I1/I2 are
outside the formalised set.

The `EventAlgebra` library (re-exported here for convenience) is an
independent, self-contained, sorry-free development of finite-dimensional
quantum event algebras: events, states, Born weights, Lüders
conditioning, the conditional expectation onto a commutative center, the
expectation functional, and the Tsirelson bound. It deliberately imports
only Mathlib and carries no vocabulary from the rest of this repository.

See `../README.md` and `../PROOF_INDEX.md` for scope and completion tracking.
-/
