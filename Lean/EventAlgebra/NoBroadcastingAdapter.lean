import EventAlgebra.PublicRecordAlgebra
import Mathlib.Analysis.InnerProductSpace.TensorProduct

/-!
# Sharp no-cloning and the no-broadcasting adapter boundary

The finite public record algebra is classical, but that fact must not be
silently upgraded to the full mixed-state no-broadcasting theorem.  This file
does two precise things:

1. it proves the standard sharp-state no-cloning obstruction directly for a
   linear isometry on a Hilbert-space tensor product;
2. it declares the small adapter through which a separately formalized finite
   no-broadcasting theorem can certify pairwise compatibility of a broadcast
   family.

The adapter's implication is data supplied by an implementation.  Merely
declaring a family broadcastable does not prove the mixed-state theorem.
-/

namespace EventAlgebra

open scoped TensorProduct

variable {H : Type*}
variable [NormedAddCommGroup H] [InnerProductSpace ℂ H]

/-- A single linear-isometric device claimed to clone two sharp states using
the same normalized blank state. -/
structure SharpCloneWitness (ψ φ blank : H) where
  copier : (H ⊗[ℂ] H) →ₗᵢ[ℂ] (H ⊗[ℂ] H)
  blank_norm : ‖blank‖ = 1
  clone_psi : copier (ψ ⊗ₜ[ℂ] blank) = ψ ⊗ₜ[ℂ] ψ
  clone_phi : copier (φ ⊗ₜ[ℂ] blank) = φ ⊗ₜ[ℂ] φ

/-- Isometric cloning forces the overlap of the two sharp states to be an
idempotent scalar. -/
theorem SharpCloneWitness.overlap_idempotent
    {ψ φ blank : H} (w : SharpCloneWitness ψ φ blank) :
    inner ℂ ψ φ * inner ℂ ψ φ = inner ℂ ψ φ := by
  have h := w.copier.inner_map_map
    (ψ ⊗ₜ[ℂ] blank) (φ ⊗ₜ[ℂ] blank)
  rw [w.clone_psi, w.clone_phi, TensorProduct.inner_tmul,
    TensorProduct.inner_tmul,
    inner_self_eq_one_of_norm_eq_one w.blank_norm, mul_one] at h
  exact h

/-- Exact sharp-state no-cloning dichotomy: two states cloned by the same
isometry have overlap zero or overlap one. -/
theorem SharpCloneWitness.overlap_zero_or_one
    {ψ φ blank : H} (w : SharpCloneWitness ψ φ blank) :
    inner ℂ ψ φ = 0 ∨ inner ℂ ψ φ = 1 := by
  apply eq_zero_or_one_of_sq_eq_self
  simpa [pow_two] using w.overlap_idempotent

/-- Any vector copied from a normalized blank by the witness has norm zero or
one.  This helper makes the normalization hidden in the cloning equations
explicit; the zero vector remains a harmless degenerate case. -/
theorem SharpCloneWitness.norm_zero_or_one
    {ψ φ blank state : H} (w : SharpCloneWitness ψ φ blank)
    (hclone : w.copier (state ⊗ₜ[ℂ] blank) = state ⊗ₜ[ℂ] state) :
    ‖state‖ = 0 ∨ ‖state‖ = 1 := by
  have hnorm := congrArg norm hclone
  rw [w.copier.norm_map, TensorProduct.norm_tmul,
    TensorProduct.norm_tmul, w.blank_norm, mul_one] at hnorm
  apply eq_zero_or_one_of_sq_eq_self
  simpa [pow_two] using hnorm.symm

/-- The overlap-one branch of exact sharp cloning consists only of identical
vectors.  No separate normalization premise is needed: the cloning equations
and normalized blank force every nonzero copied vector to have norm one. -/
theorem SharpCloneWitness.eq_of_overlap_one
    {ψ φ blank : H} (w : SharpCloneWitness ψ φ blank)
    (hoverlap : inner ℂ ψ φ = 1) : ψ = φ := by
  have hψ := w.norm_zero_or_one w.clone_psi
  have hφ := w.norm_zero_or_one w.clone_phi
  have hψ_one : ‖ψ‖ = 1 := by
    rcases hψ with hzero | hone
    · have : ψ = 0 := norm_eq_zero.mp hzero
      subst ψ
      simp at hoverlap
    · exact hone
  have hφ_one : ‖φ‖ = 1 := by
    rcases hφ with hzero | hone
    · have : φ = 0 := norm_eq_zero.mp hzero
      subst φ
      simp at hoverlap
    · exact hone
  have hscaled := (inner_eq_norm_mul_iff).mp (show
    inner ℂ ψ φ = (‖ψ‖ : ℂ) * (‖φ‖ : ℂ) by
      simpa [hψ_one, hφ_one] using hoverlap)
  simpa [hψ_one, hφ_one] using hscaled

/-- In particular, alternatives whose overlap is not one can be copied as
sharp records only if they are orthogonal. -/
theorem SharpCloneWitness.orthogonal_of_overlap_ne_one
    {ψ φ blank : H} (w : SharpCloneWitness ψ φ blank)
    (hne : inner ℂ ψ φ ≠ 1) : inner ℂ ψ φ = 0 := by
  rcases w.overlap_zero_or_one with hzero | hone
  · exact hzero
  · exact (hne hone).elim

/-- Literal sharp-record no-cloning boundary: two distinct alternatives
copied by the same isometry from the same normalized blank are orthogonal. -/
theorem SharpCloneWitness.orthogonal_of_ne
    {ψ φ blank : H} (w : SharpCloneWitness ψ φ blank)
    (hne : ψ ≠ φ) : inner ℂ ψ φ = 0 := by
  rcases w.overlap_zero_or_one with hzero | hone
  · exact hzero
  · exact (hne (w.eq_of_overlap_one hone)).elim

/-- Minimal interface expected from a finite mixed-state
no-broadcasting formalization.  `Broadcastable` is deliberately abstract;
the implementation, not this record, must prove that broadcastability implies
the chosen pairwise compatibility relation. -/
structure NoBroadcastingAdapter (State : Type*) where
  Compatible : State → State → Prop
  Broadcastable : Set State → Prop
  compatible_of_broadcastable :
    ∀ {family : Set State}, Broadcastable family →
      ∀ ⦃ρ σ : State⦄, ρ ∈ family → σ ∈ family → Compatible ρ σ

/-- A public family is objective relative to an adapter exactly when the
adapter certifies that the family is broadcastable. -/
def NoBroadcastingAdapter.IsObjectiveFamily
    {State : Type*} (adapter : NoBroadcastingAdapter State)
    (family : Set State) : Prop :=
  adapter.Broadcastable family

/-- The only conclusion exposed by the adapter boundary: members of an
objective public family satisfy the implementation's compatibility relation. -/
theorem NoBroadcastingAdapter.objective_pair_compatible
    {State : Type*} (adapter : NoBroadcastingAdapter State)
    {family : Set State} (hobj : adapter.IsObjectiveFamily family)
    {ρ σ : State} (hρ : ρ ∈ family) (hσ : σ ∈ family) :
    adapter.Compatible ρ σ :=
  adapter.compatible_of_broadcastable hobj hρ hσ

#print axioms SharpCloneWitness.overlap_idempotent
#print axioms SharpCloneWitness.overlap_zero_or_one
#print axioms SharpCloneWitness.norm_zero_or_one
#print axioms SharpCloneWitness.eq_of_overlap_one
#print axioms SharpCloneWitness.orthogonal_of_overlap_ne_one
#print axioms SharpCloneWitness.orthogonal_of_ne
#print axioms NoBroadcastingAdapter.objective_pair_compatible

end EventAlgebra
