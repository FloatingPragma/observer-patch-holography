# B6 finite holonomy and character-phase packet

Issue `#682` is closed only at the finite algebraic level.

`Screen/HolonomyInterference.lean` defines endpoint-typed finite paths and
ordered group-valued transport. Its exact results are:

- reversal-compatible edge labels make reversed transport the inverse;
- the transport ratio of two paths with common endpoints equals the holonomy
  of the closed ratio loop;
- a one-dimensional unitary character turns that identity into an exact
  relative phase;
- a vertex rechart conjugates based holonomy; the raw holonomy is invariant
  for abelian groups, while every one-dimensional character phase is
  invariant even for a noncommutative group;
- an explicit four-vertex control has trivial holonomy on every declared
  local triangular face and nontrivial holonomy on an unfilled global loop;
- the same control gives an exact two-arm finite Aharonov--Bohm phase;
- characters of `Multiplicative (ZMod n)` take values among the `n`th roots of
  unity, and the result transports through a supplied cyclic-sector
  homomorphism.

The local-flat control is a punctured finite complex: only its two declared
triangles are faces. It does not claim that every possible triangle is flat.
The phase theorem receives the edge labels, paths, character, cyclic sector,
and factorization as data. It constructs no observer-source connection,
physical gauge field, spacetime loop, flux, charge, laboratory arm, detector,
or interference fringe. Those attachments remain with the source/current and
physical-sector continuation.

The module contains no `sorry` or project axiom. Its printed dependencies are
the standard Lean/Mathlib principles only.
