import EventAlgebra.PartitionPinching

namespace EventAlgebra

open Matrix

/-!
# Complete positivity of partition pinching

The pinching map `X ↦ ∑ i, Pᵢ X Pᵢ` of a projective partition is
proven linear, unital, positive, trace preserving, and trace
self-adjoint in the parent module. This module adds the Kraus normalization: the partition
projectors themselves form a Kraus family whose completeness relation is
the partition completeness, so the pinching map is completely positive
and trace preserving with the manifest Kraus representation
`{Pᵢ}`. Public quantum records therefore commit through a physical
channel, closing the channel-form gap of the thermodynamic completion
program.
-/

variable {n k : ℕ}

/-- The partition projectors are a complete Kraus family:
`∑ i, Pᵢᴴ Pᵢ = 1`. Together with the manifest Kraus form
`partitionPinching part X = ∑ i, Pᵢ X Pᵢᴴᴴ`, this exhibits pinching as
a completely positive trace-preserving channel. -/
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
