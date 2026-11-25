import Mathlib

namespace SourceLaw.ThreeArmFlavor

/-- A linear operator commuting with the internal matrix units that take
the reference coordinate to each basis vector is scalar on the entire basis. -/
theorem scalar_from_units {V ι : Type*} [AddCommGroup V] [Module ℂ V]
    (b : Module.Basis ι ℂ V) (r : ι) (T : V →ₗ[ℂ] V)
    (M : ι → V →ₗ[ℂ] V)
    (units : ∀ i v, M i v = (b.repr v r) • b i)
    (commute : ∀ i v, T (M i v) = M i (T v)) :
    ∀ i, T (b i) = (b.repr (T (b r)) r) • b i := by
  classical
  intro i
  have h := commute i (b r)
  simpa [units] using h

/-- The interval is in the same formation units, after subtracting the
same prepared unit-filled host. An external common multiplicative unit cancels. -/
theorem mass_ratio {light heavy : ℝ} (lo : 8 ≤ light) (hi : heavy ≤ 27) :
    heavy / light ≤ 27/8 := by
  apply (div_le_iff₀ (by linarith : 0 < light)).2
  linarith

/-- Scale the three square-root energies by the smallest. Any three in
[1,2] obey the normalized bound, including repeated lines. -/
theorem square_root_triangle {x y z : ℝ} (x₁ : 1 ≤ x) (x₂ : x ≤ 2)
    (y₁ : 1 ≤ y) (y₂ : y ≤ 2) (z₁ : 1 ≤ z) (z₂ : z ≤ 2) :
    8 * (x^2+y^2+z^2) ≤ 3 * (x+y+z)^2 := by
  nlinarith [mul_nonneg (sub_nonneg.mpr x₁) (sub_nonneg.mpr x₂),
    mul_nonneg (sub_nonneg.mpr y₁) (sub_nonneg.mpr y₂),
    mul_nonneg (sub_nonneg.mpr z₁) (sub_nonneg.mpr z₂),
    sq_nonneg (x+y+z-4)]

end SourceLaw.ThreeArmFlavor
