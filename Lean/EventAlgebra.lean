import EventAlgebra.Basic
import EventAlgebra.Lueders
import EventAlgebra.PartitionPinching
import EventAlgebra.PartitionPinchingCP
import EventAlgebra.PartitionAverage
import EventAlgebra.PartitionAverageCP
import EventAlgebra.TwoScalePublicRepair
import EventAlgebra.PublicRecordAlgebra
import EventAlgebra.NoBroadcastingAdapter
import EventAlgebra.StateExpectation
import EventAlgebra.Robertson
import EventAlgebra.Superselection
import EventAlgebra.Tsirelson
import EventAlgebra.ExpectationBound
import EventAlgebra.FiniteBornFrame
import EventAlgebra.FiniteEffectClosureBoundary
import EventAlgebra.FiniteBuschGleason
import EventAlgebra.InterlockingContexts
import EventAlgebra.FiniteWebBornNoGo
import EventAlgebra.FrequencyConcentration
import EventAlgebra.RecordMajorization
import EventAlgebra.SpectralEntropyBoundary
import EventAlgebra.SlotLocalitySurface
import EventAlgebra.QuantumAdequacySurface
import EventAlgebra.SchroedingerFrameFlow
import EventAlgebra.TsirelsonSaturation
import EventAlgebra.OperationalPhaseInstrument

/-!
# Finite projection-event calculus — umbrella root

A machine-checked Lean 4 development of finite-dimensional matrix algebras
with projection events over `Matrix (Fin n) (Fin n) ℂ`:

* `EventAlgebra.Basic` — events (Hermitian idempotents), states
  (positive trace-one matrices), Born weights with reality,
  nonnegativity, normalisation, additivity, complement bound, and
  monotonicity;
* `EventAlgebra.Lueders` — Lüders conditioning: state preservation,
  repeatability, idempotence, compatibility for commuting events, the
  commuting reduction, a typed state-update boundary, and the unguarded
  fixed-point characterisation of
  conditioning;
* `EventAlgebra.PartitionPinching` — a bundled commutant star-subalgebra and
  linear pinching map for arbitrary projective partitions: exact range and
  fixed points, positivity, unitality, trace preservation, the bimodule law,
  Hilbert--Schmidt geometry, uniqueness, and Lüders compatibility;
* `EventAlgebra.PartitionAverage` — the commutative partition span with
  closure, commutativity, and centrality theorems, and the bundled averaging
  expectation onto it: exact range, trace duality, uniqueness, tower laws
  with the pinching, Born-statistics preservation, and the
  classical-conditioning collapse;
* `EventAlgebra.PartitionAverageCP`: an explicit normalized Parseval-frame
  Kraus family whose conjugation sum is the partition average.  The companion
  `Dynamics.ChoiCPTP` module binds this Kraus receipt to the repository's
  finite CP/CPTP predicate;
* `EventAlgebra.TwoScalePublicRepair`: the decomposition, invariant public
  component, residual scaling, and exponential semigroup laws for relaxation
  toward an idempotent linear publicization map;
* `EventAlgebra.PublicRecordAlgebra`: the partition span bundled as a
  commutative star subalgebra, with an exact star-algebra equivalence to
  complex-valued functions on the nonzero projector labels;
* `EventAlgebra.NoBroadcastingAdapter`: the sharp-state no-cloning
  obstruction for a common isometric copier, plus the explicit interface for
  attaching a separate finite mixed-state no-broadcasting theorem;
* `EventAlgebra.StateExpectation` — the bundled expectation functional
  `M ↦ Tr(ρ M)`: positivity, normalisation, and its
  restriction to events (the Born weight);
* `EventAlgebra.Robertson` — the supplied-state finite Robertson inequality,
  its ordinary-commutator form, and exact noncommuting saturation and
  zero-variance controls;
* `EventAlgebra.Superselection` — the exact operational quotient induced by
  a supplied partition pinching, including invisibility of every
  cross-sector corner to the complete partition commutant;
* `EventAlgebra.Tsirelson` — the Tsirelson bound `‖S‖ ≤ 2√2` for CHSH
  tuples, proved abstractly in unital C*-rings and instantiated for the
  finite matrix algebras;
* `EventAlgebra.ExpectationBound` — the state-expectation bound
  `‖Tr(ρ M)‖ ≤ ‖M‖` for the L2 operator norm and the state-level CHSH
  corollary for projection events;
* `EventAlgebra.FiniteBornFrame`: the exact rank gap for the declared
  twelve-port qubit adapter: six context-additive parameters, a
  three-parameter coordinate slice, conditional uniqueness, and an explicit
  nonrepresentation control.  The separate exact producer and independent
  verifier certify the actual projector/Bloch matrix interpretation and the
  represented nonpositive-matrix control.
* `EventAlgebra.FiniteEffectClosureBoundary`: an exact continuous nonlinear
  binary-frame countermodel on the C1 celestial sphere, showing that range,
  antipodal normalization, and continuity do not force affine/Born form;
  after affinity is supplied, dense positivity tests force the coefficient
  into the closed unit ball.

Every lemma carries a doc-comment tag, **algebra-only** or
**trace-dependent**, separating the pure `*`-algebra layer from the
results that consume the trace pairing.

The modules distinguish algebraic assumptions from trace-dependent results
while retaining the quantum-measurement interpretation in their API names.
`EventAlgebra.FiniteBuschGleason` proves the finite effect-valuation
representation: additive [0,1] assignments on effects are exactly the
Born weights of a unique density matrix, with no continuity axiom.
`EventAlgebra.InterlockingContexts` locates the decision point: the
binary sharp web admits the cube countermodel while explicit unsharp
trine and calibration contexts exclude that response.  The companion
`EventAlgebra.FiniteWebBornNoGo` proves that the whole current finite battery
still admits a transverse cubic non-Born valuation, so it does not supply the
full-effect premise of `FiniteBuschGleason`.
`EventAlgebra.FrequencyConcentration` proves the exact frequency-operator
moment and concentration theorems with the uniqueness of the stable
frequency point and its circularity boundary stated as a theorem.
`EventAlgebra.RecordMajorization` proves that arbitrary partition pinching is
the uniform random-unitary average over all independently signed block
reflections.  `EventAlgebra.SpectralEntropyBoundary` proves that Lean's
totalized matrix logarithm cannot by itself define support-aware relative
entropy: it assigns zero to an exact orthogonal pure-state pair whose support
inclusion fails.  The support-aware entropy and majorization continuation
remain separate work.
`EventAlgebra.SlotLocalitySurface` composes the Tsirelson bound and the
no-signalling receipts on one supplied bipartite slot split: the
cross-party commutation is derived from the split, the lifted local
Kraus channel preserves the remote partial trace, and the trace pairing
of remote observables is channel-invariant.
`EventAlgebra.QuantumAdequacySurface` packages the represented public
frame and re-exports the committed quantum receipts through named
premise bundles.  `EventAlgebra.SchroedingerFrameFlow` proves the
composed Born-Lueders statement with conditioning closure and the
Schroedinger frame duality: the unique Busch-Gleason state of the
flow-shifted frame follows the propagator conjugation of the unique
self-adjoint generator of the supplied flow.
-/
