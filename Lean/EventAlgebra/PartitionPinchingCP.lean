import EventAlgebra.PartitionPinching

namespace EventAlgebra

open Matrix

/-!
# Explicit Kraus data for partition pinching

The pinching map `X ↦ ∑ i, Pᵢ X Pᵢ` of a projective partition is
proven linear, unital, positive, trace preserving, and trace
self-adjoint in the parent module. This module adds the Kraus data:
the partition projectors themselves form a Kraus family whose
completeness relation is the partition completeness, and the map acts
by Kraus conjugation.  The module formalizes these two Kraus identities
but does not define or prove a finite complete-positivity predicate and does
not bundle the map as a CPTP channel.  The result is therefore reported as
explicit Kraus syntax plus the independently proved trace-preservation
identity, not as a formal CP or CPTP theorem.
-/

variable {n k : ℕ}

/-- The partition projectors are a complete Kraus family:
`∑ i, Pᵢᴴ Pᵢ = 1`. Together with the manifest Kraus form
`partitionPinching part X = ∑ i, Pᵢ X Pᵢᴴ`, this is the exact normalization
identity used by the standard Kraus argument.  No CP predicate is asserted
by this theorem. -/
theorem ProjectivePartition.kraus_complete
    (part : ProjectivePartition n k) :
    ∑ i, (part.proj i)ᴴ * part.proj i = 1 := by
  have h : ∀ i ∈ Finset.univ,
      (part.proj i)ᴴ * part.proj i = part.proj i := by
    intro i _
    rw [(part.isEvent i).1.eq]
    exact (part.isEvent i).2
  rw [Finset.sum_congr rfl h]
  exact part.complete

/-- The pinching map in explicit Kraus form: each summand is a
conjugation `Pᵢ X Pᵢᴴ` by a Kraus operator. -/
theorem partitionPinching_kraus_form
    (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionPinching part X
      = ∑ i, part.proj i * X * (part.proj i)ᴴ := by
  unfold partitionPinching
  apply Finset.sum_congr rfl
  intro i _
  rw [(part.isEvent i).1.eq]

end EventAlgebra

#print axioms EventAlgebra.ProjectivePartition.kraus_complete
#print axioms EventAlgebra.partitionPinching_kraus_form
