import EventAlgebra.PartitionAverage

/-!
# The public record algebra of a projective partition

The public records selected by a projective partition are the linear span of
its mutually orthogonal projectors.  `PartitionAverage.lean` proves the loose
closure and commutativity statements.  This module packages that span as a
star subalgebra and records the exact finite classical representation on the
nonzero (active) projectors.

Zero projectors are deliberately removed from the classical outcome type:
they carry no record and have no coefficient that can be recovered from the
matrix.  No claim about a source-selected physical partition is made here.
-/

namespace EventAlgebra

open Matrix

variable {n k : ℕ}

namespace ProjectivePartition

variable (part : ProjectivePartition n k)

/-- The partition span bundled as the finite public record star algebra. -/
noncomputable def publicSubalgebra :
    StarSubalgebra ℂ (Matrix (Fin n) (Fin n) ℂ) where
  carrier := part.span
  zero_mem' := Submodule.zero_mem _
  add_mem' := Submodule.add_mem _
  one_mem' := part.one_mem_span
  mul_mem' := part.mul_mem_span
  algebraMap_mem' c := by
    rw [Algebra.algebraMap_eq_smul_one]
    exact Submodule.smul_mem _ c part.one_mem_span
  star_mem' := part.conjTranspose_mem_span

@[simp]
theorem mem_publicSubalgebra_iff
    (X : Matrix (Fin n) (Fin n) ℂ) :
    X ∈ part.publicSubalgebra ↔ X ∈ part.span :=
  Iff.rfl

/-- The public record algebra is commutative. -/
theorem publicSubalgebra_mul_comm
    (X Y : part.publicSubalgebra) : X * Y = Y * X := by
  apply Subtype.ext
  exact part.span_mul_comm X.property Y.property

/-- Bundled inclusion of public records into the block-diagonal commutant. -/
theorem publicSubalgebra_le_commutant :
    part.publicSubalgebra ≤ part.commutant :=
  fun _ hX => part.span_le_commutant hX

/-- Active record labels: precisely the nonzero projectors. -/
def ActiveIndex := {i : Fin k // part.proj i ≠ 0}

noncomputable instance : Fintype part.ActiveIndex := by
  classical
  exact Fintype.ofFinset
    (Finset.univ.filter fun i => part.proj i ≠ 0) (by
      intro i
      change (i ∈ Finset.univ.filter fun j => part.proj j ≠ 0) ↔
        part.proj i ≠ 0
      simp)
noncomputable instance : DecidableEq part.ActiveIndex := by
  classical
  exact Classical.decEq _

/-- A nonzero projective event has nonzero trace. -/
theorem trace_ne_zero_of_active (i : part.ActiveIndex) :
    (part.proj i.1).trace ≠ 0 := by
  intro htr
  exact i.2 ((part.isEvent i.1).posSemidef.trace_eq_zero_iff.mp htr)

/-- The coefficient of a public matrix on an active record projector. -/
noncomputable def recordCoordinate
    (X : Matrix (Fin n) (Fin n) ℂ) (i : part.ActiveIndex) : ℂ :=
  ((part.proj i.1).trace)⁻¹ * bornWeight X (part.proj i.1)

@[simp]
theorem recordCoordinate_proj_self (i : part.ActiveIndex) :
    part.recordCoordinate (part.proj i.1) i = 1 := by
  rw [recordCoordinate, bornWeight, (part.isEvent i.1).2]
  exact inv_mul_cancel₀ (part.trace_ne_zero_of_active i)

theorem recordCoordinate_proj_ne
    (i j : part.ActiveIndex) (hij : i ≠ j) :
    part.recordCoordinate (part.proj i.1) j = 0 := by
  have hval : i.1 ≠ j.1 := by
    intro h
    apply hij
    exact Subtype.ext h
  rw [recordCoordinate, bornWeight, part.orthogonal i.1 j.1 hval,
    trace_zero, mul_zero]

@[simp]
theorem recordCoordinate_one (i : part.ActiveIndex) :
    part.recordCoordinate 1 i = 1 := by
  rw [recordCoordinate, bornWeight, one_mul]
  exact inv_mul_cancel₀ (part.trace_ne_zero_of_active i)

/-- Synthesis of a classical function as a public record matrix. -/
noncomputable def recordSynthesis
    (f : part.ActiveIndex → ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  ∑ i, f i • part.proj i.1

theorem recordSynthesis_mem (f : part.ActiveIndex → ℂ) :
    part.recordSynthesis f ∈ part.publicSubalgebra := by
  exact Submodule.sum_mem _ fun i _ =>
    Submodule.smul_mem _ _ (part.proj_mem_span i.1)

@[simp]
theorem recordSynthesis_zero :
    part.recordSynthesis (0 : part.ActiveIndex → ℂ) = 0 := by
  simp [recordSynthesis]

@[simp]
theorem recordSynthesis_add
    (f g : part.ActiveIndex → ℂ) :
    part.recordSynthesis (f + g) =
      part.recordSynthesis f + part.recordSynthesis g := by
  simp [recordSynthesis, Finset.sum_add_distrib, add_smul]

@[simp]
theorem recordSynthesis_smul
    (c : ℂ) (f : part.ActiveIndex → ℂ) :
    part.recordSynthesis (c • f) = c • part.recordSynthesis f := by
  simp [recordSynthesis, Finset.smul_sum, smul_smul]

@[simp]
theorem recordSynthesis_mul
    (f g : part.ActiveIndex → ℂ) :
    part.recordSynthesis (f * g) =
      part.recordSynthesis f * part.recordSynthesis g := by
  classical
  rw [recordSynthesis, recordSynthesis, recordSynthesis, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro i _
  rw [Finset.mul_sum, Finset.sum_eq_single i]
  · simp [(part.isEvent i.1).2, smul_smul, mul_comm]
  · intro j _ hji
    have hval : i.1 ≠ j.1 := by
      intro h
      apply hji
      exact Subtype.ext h.symm
    rw [smul_mul_assoc, mul_smul_comm,
      part.orthogonal i.1 j.1 hval]
    simp
  · intro hi
    exact absurd (Finset.mem_univ i) hi

@[simp]
theorem recordSynthesis_star
    (f : part.ActiveIndex → ℂ) :
    (part.recordSynthesis f)ᴴ = part.recordSynthesis (star f) := by
  classical
  rw [recordSynthesis, recordSynthesis, conjTranspose_sum]
  apply Finset.sum_congr rfl
  intro i _
  rw [conjTranspose_smul, (part.isEvent i.1).1.eq]
  rfl

@[simp]
theorem recordCoordinate_recordSynthesis
    (f : part.ActiveIndex → ℂ) (j : part.ActiveIndex) :
    part.recordCoordinate (part.recordSynthesis f) j = f j := by
  rw [recordSynthesis, recordCoordinate, bornWeight, Finset.sum_mul,
    trace_sum, Finset.sum_eq_single j]
  · rw [smul_mul_assoc, (part.isEvent j.1).2, trace_smul, smul_eq_mul]
    field_simp [part.trace_ne_zero_of_active j]
  · intro i _ hij
    have hval : i.1 ≠ j.1 := by
      intro h
      apply hij
      exact Subtype.ext h
    simp [part.orthogonal i.1 j.1 hval]
  · intro hj
    exact absurd (Finset.mem_univ j) hj

/-- Active projectors are linearly independent. -/
theorem active_projectors_linearIndependent :
    LinearIndependent ℂ (fun i : part.ActiveIndex => part.proj i.1) := by
  rw [Fintype.linearIndependent_iff]
  intro f hf i
  have hcoord := congr_arg (fun X => part.recordCoordinate X i) hf
  change part.recordCoordinate (part.recordSynthesis f) i =
    part.recordCoordinate 0 i at hcoord
  rw [part.recordCoordinate_recordSynthesis] at hcoord
  simpa [recordCoordinate, bornWeight] using hcoord

/-- Every public matrix is reconstructed by its active record coordinates. -/
theorem recordSynthesis_recordCoordinate
    (X : Matrix (Fin n) (Fin n) ℂ) (hX : X ∈ part.publicSubalgebra) :
    part.recordSynthesis (part.recordCoordinate X) = X := by
  have havg : partitionAverage part X = X :=
    partitionAverage_fixes part hX
  classical
  change (∑ i : part.ActiveIndex,
    (((part.proj i.1).trace)⁻¹ * bornWeight X (part.proj i.1)) •
      part.proj i.1) = X
  calc
    _ = ∑ i : Fin k,
        (((part.proj i).trace)⁻¹ * bornWeight X (part.proj i)) •
          part.proj i := by
      let term : Fin k → Matrix (Fin n) (Fin n) ℂ := fun i =>
        (((part.proj i).trace)⁻¹ * bornWeight X (part.proj i)) • part.proj i
      have hsub :
          (∑ i : part.ActiveIndex, term i.1) =
            (∑ i ∈ Finset.univ.filter (fun i => part.proj i ≠ 0), term i) := by
        simpa [ActiveIndex] using
          (Finset.sum_subtype_eq_sum_filter
            (s := Finset.univ) (p := fun i : Fin k => part.proj i ≠ 0) term)
      have hfilter :
          (∑ i ∈ Finset.univ.filter (fun i => part.proj i ≠ 0), term i) =
            (∑ i : Fin k, term i) := by
        apply Finset.sum_subset (Finset.filter_subset _ _)
        intro i _ hi
        have hzero : part.proj i = 0 := by
          by_contra hne
          exact hi (by simp [hne])
        simp [term, hzero, bornWeight]
      exact hsub.trans hfilter
    _ = partitionAverage part X := rfl
    _ = X := havg

@[simp]
theorem recordSynthesis_one :
    part.recordSynthesis (1 : part.ActiveIndex → ℂ) = 1 := by
  have h := part.recordSynthesis_recordCoordinate 1 part.one_mem_span
  have hc : part.recordCoordinate 1 = (1 : part.ActiveIndex → ℂ) := by
    funext i
    exact part.recordCoordinate_one i
  rwa [hc] at h

/-- Synthesis as a star-algebra homomorphism from the finite classical
function algebra into the public matrix algebra. -/
noncomputable def recordSynthesisStarAlgHom :
    (part.ActiveIndex → ℂ) →⋆ₐ[ℂ] part.publicSubalgebra where
  toFun f := ⟨part.recordSynthesis f, part.recordSynthesis_mem f⟩
  map_zero' := by
    apply Subtype.ext
    exact part.recordSynthesis_zero
  map_one' := by
    apply Subtype.ext
    exact part.recordSynthesis_one
  map_add' f g := by
    apply Subtype.ext
    exact part.recordSynthesis_add f g
  map_mul' f g := by
    apply Subtype.ext
    exact part.recordSynthesis_mul f g
  commutes' c := by
    apply Subtype.ext
    change part.recordSynthesis ((algebraMap ℂ (part.ActiveIndex → ℂ)) c) =
      (algebraMap ℂ (Matrix (Fin n) (Fin n) ℂ)) c
    rw [Algebra.algebraMap_eq_smul_one, Algebra.algebraMap_eq_smul_one]
    rw [← part.recordSynthesis_one, ← part.recordSynthesis_smul]
  map_star' f := by
    apply Subtype.ext
    exact (part.recordSynthesis_star f).symm

theorem recordSynthesisStarAlgHom_bijective :
    Function.Bijective part.recordSynthesisStarAlgHom := by
  constructor
  · intro f g h
    funext i
    have hval : part.recordSynthesis f = part.recordSynthesis g :=
      congr_arg Subtype.val h
    have hcoord := congr_arg (fun X => part.recordCoordinate X i) hval
    simpa using hcoord
  · intro X
    refine ⟨part.recordCoordinate X.1, ?_⟩
    apply Subtype.ext
    exact part.recordSynthesis_recordCoordinate X.1 X.2

/-- Exact finite classical representation of the public record algebra:
public matrices are star-algebra equivalent to complex-valued functions on
the active projector labels. -/
noncomputable def publicRecordFunctionEquiv :
    part.publicSubalgebra ≃⋆ₐ[ℂ] (part.ActiveIndex → ℂ) :=
  (StarAlgEquiv.ofBijective part.recordSynthesisStarAlgHom
    part.recordSynthesisStarAlgHom_bijective).symm

@[simp]
theorem publicRecordFunctionEquiv_apply
    (X : part.publicSubalgebra) (i : part.ActiveIndex) :
    part.publicRecordFunctionEquiv X i = part.recordCoordinate X.1 i := by
  let e : (part.ActiveIndex → ℂ) ≃⋆ₐ[ℂ] part.publicSubalgebra :=
    StarAlgEquiv.ofBijective part.recordSynthesisStarAlgHom
      part.recordSynthesisStarAlgHom_bijective
  have he : e (part.recordCoordinate X.1) = X := by
    apply Subtype.ext
    exact part.recordSynthesis_recordCoordinate X.1 X.2
  have hs : e.symm X = part.recordCoordinate X.1 := by
    exact e.symm_apply_eq.mpr he.symm
  exact congr_fun hs i

-- Axiom audit: all statements must remain within the standard Mathlib basis.
#print axioms ProjectivePartition.publicSubalgebra_mul_comm
#print axioms ProjectivePartition.publicSubalgebra_le_commutant
#print axioms ProjectivePartition.trace_ne_zero_of_active
#print axioms ProjectivePartition.recordCoordinate_recordSynthesis
#print axioms ProjectivePartition.active_projectors_linearIndependent
#print axioms ProjectivePartition.recordSynthesis_recordCoordinate
#print axioms ProjectivePartition.recordSynthesisStarAlgHom_bijective
#print axioms ProjectivePartition.publicRecordFunctionEquiv_apply

end ProjectivePartition

end EventAlgebra
