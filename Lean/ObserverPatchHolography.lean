import ObserverPatchHolography.AbstractRewriting
import ObserverPatchHolography.CoreAxioms
import ObserverPatchHolography.Primitives
import ObserverPatchHolography.Rule90
import ObserverPatchHolography.BoundaryFiber
import ObserverPatchHolography.BridgeEquivalence
import ObserverPatchHolography.CapacityFixedPoint
import ObserverPatchHolography.CapacityClosurePrinciple
import ObserverPatchHolography.ClosurePreflight
import ObserverPatchHolography.CapacityNonidentifiability
import ObserverPatchHolography.GlobalCapacityAttachment
import ObserverPatchHolography.SeedPi
import ObserverPatchHolography.ClebschRatio
import ObserverPatchHolography.KoideCirculant
import ObserverPatchHolography.W5Stabilizer
import ObserverPatchHolography.IcosahedralAxisNoGo
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
import ObserverPatchHolography.DeSitterCapacityShock
import EventAlgebra
import ObserverPatchHolography.YangMillsLemma72
import ObserverPatchHolography.YangMillsProp81
import ObserverPatchHolography.YangMillsGap
import ObserverPatchHolography.YangMillsGapWitness
import ObserverPatchHolography.RepairGapChain
import ObserverPatchHolography.QuotientLumpability
import ObserverPatchHolography.QuotientLumpabilityRuntimeHarness
import ObserverPatchHolography.ScalarSeamRepair
import ObserverPatchHolography.EqualSeamSelection
import ObserverPatchHolography.RepairWordSchedule
import Time.TimeOrderLedger
import Time.ObserverHistory
import Time.ClockReadout
import Time.WorldlineRealization
import Time.ProperTimeCalibration
import Time.ClockComparison
import Tower
import ObserverPatchHolography.A2EndpointCommutator
import ObserverPatchHolography.DirectedSeamRepair
import ObserverPatchHolography.DirectedSeamRepairProgress
import ObserverPatchHolography.RateNonidentifiability
import ObserverPatchHolography.RateBridgeObstruction
import ObserverPatchHolography.Locality.DependencyCone
import ObserverPatchHolography.Locality.NoSignalling
import ObserverPatchHolography.Execution.AdaptiveRunCounterexamples
import ObserverPatchHolography.Execution.CumulativeAttemptCapacity
import ObserverPatchHolography.Execution.CumulativeAttemptCapacityExamples
import ObserverPatchHolography.Provenance.SemanticEventProvenance
import ObserverPatchHolography.Provenance.FiniteCausetCompiler
import ObserverPatchHolography.Provenance.FiniteCausetCoverCompiler
import ObserverPatchHolography.Provenance.MismatchProvenance
import ObserverPatchHolography.Provenance.SeamDeltaAggregation
import ObserverPatchHolography.Provenance.CausalInterval
import ObserverPatchHolography.Provenance.HistoryCausalInvariance
import ObserverPatchHolography.Provenance.QuotientInvariance
import ObserverPatchHolography.Provenance.RefinementNaturality
import Computation.FixedFederationProgress
import Computation.FixedFederationCounterexamples
import Computation.FixedFederationExecution
import Computation.FixedFederationExecutionExamples
import Computation.FixedFederationComplexity
import Computation.FixedFederationFanoutControl
import ObserverPatchHolography.Execution.RankedAttemptCapacity

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
coupling without excluding it (T1), and on the fixed 4×4 model the naive
modular recast does not exclude the non-central witness while the
corrected recast yields the diagonal clause, which coincides with the
paper's MSA form and is circular as a derivation (T2).

The `BridgeEquivalence`, `CapacityFixedPoint`, `CapacityClosurePrinciple`,
`ClosurePreflight`, and `SeedPi` modules carry
the Part-A coupling-algebra layer: the bridge count/tick equivalence, the
capacity fixed-point uniqueness schema, the typed closure principle, the
arbitrary-seed normalization control, and the CAP-P seed statement. They
formalise the algebraic layer only; the source selector, same-quantity
identity, and physical identities I1/I2 are outside the formalised set.

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

The `Time.TimeOrderLedger` module carries the completion-plan A1 ontology
boundary.  Closure witnesses, repair schedules, observer record orders,
modular parameters, worldline realizations, clock readouts, proper-time
intervals, and optional global-time functions are distinct types. The
committed source environment rejects the audited implicit coercions, and
cross-layer interpretations require explicit named realization maps. This is
type infrastructure only: the module constructs no physical clock, worldline,
global time function, or modular-time identity.

The `Time.ObserverHistory` through `Time.ClockComparison` modules carry the
bounded D1 continuation. Bare monotone clocks admit arbitrary strictly
increasing, including non-affine, regrading. A supplied affine unit-timelike
history yields chart-covariant timelike displacements and additive
dimensionless proper intervals; two ordered shared events uniquely determine
a positive affine comparison. A distinct third shared event supplies a
consistency test with no additional fit parameter. The source history,
physical instrument, refinement and SI scale remain unconstructed.

The `Tower.ConsensusTower` module carries the completion-plan A3 root
interface. It packages finite observer and record fibers, private matrix
algebras, commutative public subalgebras, selected states and generators, and
functorial refinement compatibility in one directed tower. Its constant
projective-partition adaptor proves only that an existing finite fiber can be
packaged. No nonconstant source tower, repair endpoint, causal net, geometry,
clock, continuum limit, or physical evolution is constructed.

See `../README.md` and `../PROOF_INDEX.md` for scope and proof coverage.
-/
