# B6 finite holonomy and character-phase packet

Issue `#682` is closed only at the finite algebraic level.

`Screen/HolonomyInterference.lean` defines endpoint-typed finite paths and
ordered group-valued transport. Its exact results are:

- reversal-compatible edge labels make reversed transport the inverse;
- the transport ratio of two paths with common endpoints equals the holonomy
  of the closed ratio loop;
- a one-dimensional unitary character turns that identity into an exact
  relative phase;
- a vertex rechart conjugates based holonomy, so its conjugacy class,
  nontriviality, and every class-function value are invariant; the raw
  holonomy is invariant for abelian groups, while every one-dimensional
  character phase is invariant even for a noncommutative group;
- an explicit four-vertex control has trivial holonomy on its two declared
  local triangular faces and nontrivial holonomy on a separate undeclared,
  unfilled loop;
- the same control gives an exact two-arm character-phase identity, an
  algebraic finite Aharonov--Bohm analogy only;
- characters of `Multiplicative (ZMod n)` take values among the `n`th roots of
  unity, and the result transports through a supplied cyclic-sector
  homomorphism.

The local-flat example is a four-vertex path/face control with exactly two
declared triangular faces and a separate undeclared, unfilled loop. It does
not claim that every possible triangle is a face or flat, and proves no
puncture or noncontractibility statement. The Aharonov--Bohm language refers
only to the algebraic two-path/character-phase analogy, not to a topological
or physical identification.
The phase theorem receives the edge labels, paths, character, cyclic sector,
and factorization as data. It constructs no observer-source connection,
physical gauge field, spacetime loop, flux, charge, laboratory arm, detector,
or interference fringe. Those attachments remain with the source/current and
physical-sector continuation.

The module contains no `sorry` or project axiom. Its printed dependencies are
the standard Lean/Mathlib principles only.
