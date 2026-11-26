import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

namespace OPH.SharedLoadReference

/-- Adjacent-charge overlap of a normalized one-load reference. -/
theorem one_load_bound (a b : ℝ) (norm : a^2 + b^2 = 1) : 2 * (a*b) ≤ 1 := by
  nlinarith [sq_nonneg (a - b)]

/-- Deficit in the two-load visibility bound, for arbitrary multiplicity-sector norms. -/
theorem two_load_deficit (a b c : ℝ) :
    (a^2 + b^2 + c^2)^2 - 2 * (a*b + b*c)^2 =
      (a^2 + c^2 - b^2)^2 + 2*b^2*(a-c)^2 := by ring

theorem two_load_bound (a b c : ℝ) (norm : a^2 + b^2 + c^2 = 1) :
    2 * (a*b + b*c)^2 ≤ 1 := by
  have h := two_load_deficit a b c
  rw [norm] at h
  nlinarith [sq_nonneg (a^2+c^2-b^2), mul_nonneg (sq_nonneg b) (sq_nonneg (a-c))]

/-- The binomial split of two independent equal references attains equality. -/
theorem balanced_two_attains (middle : ℝ) (hm : middle^2 = 1/2) :
    2 * ((1/2 : ℝ)*middle + middle*(1/2))^2 = 1 := by nlinarith

end OPH.SharedLoadReference
