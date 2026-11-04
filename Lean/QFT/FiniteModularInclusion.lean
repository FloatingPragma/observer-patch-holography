import Mathlib.LinearAlgebra.FiniteDimensional.Basic
import Mathlib.Data.Real.Basic

/-!
# Finite-dimensional modular inclusions

An invertible linear map carrying a finite-dimensional subspace into itself
carries it onto itself. Its inverse therefore preserves the same subspace.
Applied to the vector space underlying a finite matrix subalgebra, this
excludes strict half-sided compression by modular automorphisms.

The source-specific CAR algebra and its nonzero modular leakage are checked
in `model/finite_modular_inclusion`. This module supplies the dimension
argument, without identifying a finite modular flow with null translation.
-/

namespace OPH.QFT.FiniteModularInclusion

variable {K V : Type*} [DivisionRing K] [AddCommGroup V] [Module K V]

/-- Finite dimension is needed only for the preserved subspace. -/
theorem inverse_preserves (S : Submodule K V) [FiniteDimensional K S]
    (e : V ≃ₗ[K] V) (h : ∀ x ∈ S, e x ∈ S) :
    ∀ x ∈ S, e.symm x ∈ S := by
  let f : S →ₗ[K] S :=
    { toFun := fun x => ⟨e x, h x x.property⟩
      map_add' := by intro x y; apply Subtype.ext; exact e.map_add x y
      map_smul' := by intro c x; apply Subtype.ext; exact e.map_smul c x }
  have hinj : Function.Injective f := by
    intro x y hxy
    apply Subtype.ext
    exact e.injective (congrArg Subtype.val hxy)
  have hsurj : Function.Surjective f := LinearMap.injective_iff_surjective.mp hinj
  intro x hx
  obtain ⟨y, hy⟩ := hsurj ⟨x, hx⟩
  have he : e (y : V) = x := congrArg Subtype.val hy
  rw [← he, e.symm_apply_apply]
  exact y.property

/-- A one-sided finite image inclusion is equality. -/
theorem image_eq (S : Submodule K V) [FiniteDimensional K S]
    (e : V ≃ₗ[K] V) (h : ∀ x ∈ S, e x ∈ S) :
    S.map e.toLinearMap = S := by
  apply le_antisymm
  · rintro x ⟨y, hy, rfl⟩
    exact h y hy
  · intro x hx
    exact ⟨e.symm x, inverse_preserves S e h x hx, e.apply_symm_apply x⟩

/-- Positive-time invariance entails invariance at every real parameter.
No continuity assumption is needed for this finite-dimensional statement. -/
theorem half_sided_is_two_sided (S : Submodule K V) [FiniteDimensional K S]
    (flow : ℝ → V ≃ₗ[K] V)
    (hinv : ∀ t, flow (-t) = (flow t).symm)
    (hpos : ∀ t, 0 ≤ t → ∀ x ∈ S, flow t x ∈ S) :
    ∀ t, ∀ x ∈ S, flow t x ∈ S := by
  intro t x hx
  by_cases ht : 0 ≤ t
  · exact hpos t ht x hx
  · have hneg : 0 ≤ -t := neg_nonneg.mpr (le_of_not_ge ht)
    have hi := inverse_preserves S (flow (-t)) (hpos (-t) hneg) x hx
    simpa only [← hinv (-t), neg_neg] using hi

end OPH.QFT.FiniteModularInclusion
