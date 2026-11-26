import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

namespace CommonNucleonScale

/-- a and b are the predictions divided by their respective targets after the
same positive unit conversion. Any common relative tolerance must span both. -/
theorem necessary_tolerance (a b e : ℝ) (ha : 0 ≤ a) (hb : 0 ≤ b)
    (fitA : |a - 1| ≤ e) (fitB : |b - 1| ≤ e) :
    |b - a| ≤ e * (a + b) := by
  obtain ⟨al, au⟩ := abs_le.mp fitA
  obtain ⟨bl, bu⟩ := abs_le.mp fitB
  have aLo : 1-e ≤ a := by linarith
  have bLo : 1-e ≤ b := by linarith
  have aHi : a ≤ 1+e := by linarith
  have bHi : b ≤ 1+e := by linarith
  apply abs_le.mpr
  constructor
  · nlinarith only [mul_le_mul_of_nonneg_left bLo ha, mul_le_mul_of_nonneg_right aHi hb]
  · nlinarith only [mul_le_mul_of_nonneg_left aLo hb, mul_le_mul_of_nonneg_right bHi ha]

/-- Scaling cancels in the normalized difference and sum. -/
theorem scale_cancels (s a b : ℝ) (hs : 0 < s) :
    |s*b-s*a|/(s*a+s*b) = |b-a|/(a+b) := by
  rw [← mul_sub, ← mul_add, abs_mul, abs_of_pos hs]
  exact mul_div_mul_left _ _ hs.ne'

end CommonNucleonScale
