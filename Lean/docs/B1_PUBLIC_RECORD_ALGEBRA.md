# B1 public record algebra boundary

Issue `#677` asks for the finite classical algebra carried by a projective
record partition and for the exact boundary between sharp no-cloning and the
mixed-state no-broadcasting theorem.

## Proved in Lean

`EventAlgebra/PublicRecordAlgebra.lean` bundles the span of a finite
projective partition as `ProjectivePartition.publicSubalgebra`. The bundled
algebra is commutative and lies in the partition commutant. Its active labels
are exactly the nonzero projectors. Coordinate extraction and synthesis are
mutual inverses, and `publicRecordFunctionEquiv` gives a star-algebra
equivalence

\[
  \mathcal D_{\mathrm{rec}} \simeq
  \{f : I_{\mathrm{active}} \to \mathbb C\}.
\]

Zero projectors are excluded from the label set because their coefficients
cannot be recovered from a matrix. This avoids adding spurious classical
outcomes.

`EventAlgebra/NoBroadcastingAdapter.lean` proves the sharp-state obstruction.
If one linear isometry copies two states using the same normalized blank, then
their inner product is idempotent. Over the complex numbers it is therefore
zero or one. The cloning equations also force every copied vector to have norm
zero or one, and overlap one forces equality. Therefore two distinct
alternatives copied by the same device are orthogonal.

## Boundary

The sharp theorem is not the general mixed-state no-broadcasting theorem.
`NoBroadcastingAdapter` exposes the exact plug-in point: an implementation
supplies a state type, a compatibility relation, a broadcastability predicate,
and the theorem that members of a broadcastable family are compatible. The
adapter derives pairwise compatibility for a supplied objective family; it
does not construct a broadcaster or prove the mixed-state implication.

The projective partition, matrix representation, and public-family
broadcastability receipt are inputs. No source theorem selects a physical
record partition, no laboratory record is identified, and no prediction is
emitted.
