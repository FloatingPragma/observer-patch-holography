import ObserverPatchHolography.AbstractRewriting
import ObserverPatchHolography.Primitives
import ObserverPatchHolography.Rule90
import ObserverPatchHolography.BoundaryFiber
import ObserverPatchHolography.BridgeEquivalence
import ObserverPatchHolography.CapacityFixedPoint
import ObserverPatchHolography.SeedPi
import ObserverPatchHolography.CollarClause
import ObserverPatchHolography.CollarLayer
import ObserverPatchHolography.CollarStates
import ObserverPatchHolography.CollarStatesT1
import ObserverPatchHolography.CollarModularT2
import ObserverPatchHolography.CollarStatesBridge
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
module formalises the OPH primitives (Records, Repair, Patch, Obs, Φ, gauge
equivalence, OPH-Confluence, OPH-Completeness) admission-free — the three
former admissions (Lyapunov descent, termination, single-site solvability)
are discharged; these structurally depend on the companion paper *Reality
as a Consensus Protocol*.

The `CollarClause`/`CollarLayer` modules carry the issue #544 layer
separation (the overlap-consistency layer factors through the realized
constraint family; the collar clause is a declared input, not a theorem),
and the `CollarStates`/`CollarStatesT1`/`CollarModularT2` modules carry the
state-side no-gos: the current state-side axioms do not force the clause
(T0), the flux conditional expectation exists and deselects — but does not
exclude — the cross-cut coupling (T1), and the naive modular recast is
vacuous while the corrected recast buys only the diagonal clause (T2).

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
