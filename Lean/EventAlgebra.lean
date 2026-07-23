import EventAlgebra.Basic
import EventAlgebra.Lueders
import EventAlgebra.PartitionPinching
import EventAlgebra.PartitionAverage
import EventAlgebra.StateExpectation
import EventAlgebra.Tsirelson
import EventAlgebra.ExpectationBound

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
* `EventAlgebra.StateExpectation` — the bundled expectation functional
  `M ↦ Tr(ρ M)`: positivity, normalisation, and its
  restriction to events (the Born weight);
* `EventAlgebra.Tsirelson` — the Tsirelson bound `‖S‖ ≤ 2√2` for CHSH
  tuples, proved abstractly in unital C*-rings and instantiated for the
  finite matrix algebras;
* `EventAlgebra.ExpectationBound` — the state-expectation bound
  `‖Tr(ρ M)‖ ≤ ‖M‖` for the L2 operator norm and the state-level CHSH
  corollary for projection events.

Every lemma carries a doc-comment tag, **algebra-only** or
**trace-dependent**, separating the pure `*`-algebra layer from the
results that consume the trace pairing.

The modules distinguish algebraic assumptions from trace-dependent results
while retaining the quantum-measurement interpretation in their API names.
-/
