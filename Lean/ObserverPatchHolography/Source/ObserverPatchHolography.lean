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
import ObserverPatchHolography.IcosahedralAntibridge
import ObserverPatchHolography.IcosahedralOrbitStabilizer
import ObserverPatchHolography.BridgeBoundaries
import ObserverPatchHolography.EinsteinBranch
import EventAlgebra

/-!
# Observer-Patch Holography : Lean 4 umbrella root

Re-exports Jonathan Hill's concrete carrier/dynamics modules, the neutral
observation-determined normal-forms proof package, and the explicit bridge
between them.

**Status: preliminary skeleton rather than theorem-grade formalisation of
Proposition 4.2** from *Paradise as Fixed-Point Consensus*. The `Primitives`
module formalises the OPH primitives (Records, Repair, Patch, Obs, Φ, gauge
equivalence, OPH-Confluence, OPH-Completeness) admission-free : the three
former admissions (Lyapunov descent, termination, single-site solvability)
are discharged; these structurally depend on the companion paper *Reality
as a Consensus Protocol*.

The `CollarClause`/`CollarLayer` modules carry the issue #544 layer
separation (the overlap-consistency layer factors through the realized
constraint family; the collar clause is a declared input, not a theorem),
and the `CollarStates`/`CollarStatesT1`/`CollarModularT2` modules carry the
state-side no-gos: the stated state-side axioms do not force the clause
(T0), the flux conditional expectation exists and deselects the cross-cut
coupling without excluding it (T1), and the naive modular recast is
vacuous while the corrected recast buys only the diagonal clause (T2).

The `BridgeEquivalence`, `CapacityFixedPoint`, and `SeedPi` modules carry
the Part-A coupling-algebra layer: the bridge count/tick equivalence, the
capacity fixed-point uniqueness schema, and the CAP-P seed statement. They
formalise the algebraic layer only; the physical identities I1/I2 are
outside the formalised set.

The `BridgeBoundaries` module carries the finite anti-bridges used by the
consensus correction: disjoint writes can fail to compose when a nonlinear
protected observable is omitted from the dependency graph, and equal one- and
two-bit marginals do not determine a tripartite payload. Positive local-diamond
and coherent union-collar theorems therefore keep their receipts explicit.

The `EinsteinBranch` modules carry the issue #578 algebraic/compositional
kernel: bare-tower non-entailment, common-domain typed arrows and
boundary-fibre composition, explicit null tomography, finite entropy/MaxEnt
identities, exact small-ball coefficient arithmetic, timelike and null tensor
algebra, Ward/Bianchi constancy, and strict manifest deletion logic.  The
continuum, asymptotic, and physical premises remain explicit theorem inputs;
no inhabited Einstein-admissible tower is claimed.

The `EventAlgebra` library (re-exported here for convenience) is an
independent, self-contained, sorry-free development of finite-dimensional
quantum event algebras: events, states, Born weights, Lüders
conditioning, the conditional expectation onto a commutative center, the
expectation functional, and the Tsirelson bound. It deliberately imports
only Mathlib and carries no vocabulary from the rest of this repository.

See `../README.md` and `../PROOF_INDEX.md` for scope and proof coverage.
-/
