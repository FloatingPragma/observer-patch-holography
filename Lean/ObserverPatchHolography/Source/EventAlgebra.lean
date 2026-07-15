import EventAlgebra.Basic
import EventAlgebra.Lueders
import EventAlgebra.CenterExpectation
import EventAlgebra.StateFromTrace
import EventAlgebra.Tsirelson

/-!
# Finite event algebras — umbrella root

A sorry-free Lean 4 development of finite-dimensional
quantum event algebras over `Matrix (Fin n) (Fin n) ℂ`:

* `EventAlgebra.Basic` — events (Hermitian idempotents), states
  (positive trace-one matrices), Born weights with reality,
  nonnegativity, normalisation, additivity, complement bound, and
  monotonicity;
* `EventAlgebra.Lueders` — Lüders conditioning: state preservation,
  repeatability, idempotence, compatibility for commuting events, the
  classical restriction, and the fixed-point characterisation of
  conditioning;
* `EventAlgebra.CenterExpectation` — the conditional expectation
  (pinching) onto the block-diagonal commutant of a partition of
  unity: projector properties, state preservation, trace
  selfadjointness, Pythagoras, contractivity, uniqueness, and
  compatibility with Lüders conditioning on commutant events;
* `EventAlgebra.StateFromTrace` — the expectation functional
  `M ↦ Tr(ρ M)`: linearity, positivity, normalisation, and its
  restriction to events (the Born weight);
* `EventAlgebra.Tsirelson` — the Tsirelson bound `‖S‖ ≤ 2√2` for CHSH
  tuples, proved abstractly in unital C*-rings and instantiated for the
  finite matrix algebras.

Every lemma carries a doc-comment tag, **algebra-only** or
**trace-dependent**, separating the pure `*`-algebra layer from the
results that consume the trace pairing.

These modules are deliberately free of any interpretational vocabulary:
plain finite-dimensional operator algebra, importing only Mathlib.
-/
