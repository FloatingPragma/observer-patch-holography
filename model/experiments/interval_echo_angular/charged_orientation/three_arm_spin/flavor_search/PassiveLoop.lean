import Mathlib

namespace SourceLaw.PassiveLoop

/-- If the second completed lap has only the blank component, every
code restored by an invertible first lap belongs to that component.
A nonzero scalar times a unitary is one such invertible code action. -/
theorem restored_code {V D : Type*} [AddCommGroup V] [Module ℂ V]
    [AddCommGroup D] [Module ℂ D]
    (W P : V →ₗ[ℂ] V) (C : D →ₗ[ℂ] V) (T : D ≃ₗ[ℂ] D) (a : ℂ)
    (project : ∀ v, P (P v) = P v)
    (twice : ∀ v, W (W v) = a ^ 2 • P v)
    (blank : ∀ v, W (P v) = a • P v)
    (restore : ∀ x, W (C x) = C (T x)) :
    (∀ x, P (C x) = C x) ∧ (∀ x, C (T x) = a • C x) := by
  have support (v : V) : P (W (W v)) = W (W v) := by
    rw [twice, map_smul, project]
  have code (x : D) : P (C x) = C x := by
    obtain ⟨y, rfl⟩ := T.surjective x
    obtain ⟨z, rfl⟩ := T.surjective y
    simpa only [restore] using support (C z)
  refine ⟨code, fun x => ?_⟩
  simpa only [code, restore] using blank (C x)

end SourceLaw.PassiveLoop
