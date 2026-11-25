import Mathlib

namespace SourceLaw.ChargedCharacter

/-- If a source vector and a terminal bra transform by different characters
under a common rotation, the joint detector amplitude is zero. -/
theorem selection {V : Type*} [AddCommGroup V] [Module ℂ V]
    (T : V →ₗ[ℂ] V) (read : V →ₗ[ℂ] ℂ) (x : V) (a b : ℂ)
    (src : T x = a • x) (dst : read (T x) = b * read x) (distinct : a ≠ b) :
    read x = 0 := by
  rw [src, map_smul, smul_eq_mul] at dst
  have h : (a-b) * read x = 0 := by rw [sub_mul, dst, sub_self]
  exact (mul_eq_zero.mp h).resolve_left (sub_ne_zero.mpr distinct)

/-- A reflection-symmetric three-arm moment matrix with diagonal d and
off-diagonal c has the same moment d-c on both nontrivial characters. -/
theorem two_winding_moments (d c w : ℂ) (root : w^2+w+1=0) :
    d+c*w+c*w^2=d-c ∧ d+c*w^2+c*w=d-c := by
  constructor <;> linear_combination c*root

/-- The unbiased-basis overlap inequality and the Schmidt-two fidelity limit
are incompatible with an observed match sum above five thirds. -/
theorem Schmidt_bound_excluded (a b F : ℝ) (overlap : a+b≤1+F)
    (observed : 5/3<a+b) : ¬ F≤2/3 := by linarith

end SourceLaw.ChargedCharacter
