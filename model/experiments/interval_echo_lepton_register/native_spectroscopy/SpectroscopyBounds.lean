import Mathlib

namespace OPH.NativeSpectroscopy

/-- Positive formation lines in a span of at most 64 in energy have Koide
quotient strictly below two thirds.  The variables are positive square roots. -/
theorem koide_below_balance {a b c : ℝ} (ha : 0 < a) (hab : a ≤ b) (hbc : b ≤ c)
    (hc : c ≤ 8 * a) :
    3 * (a ^ 2 + b ^ 2 + c ^ 2) < 2 * (a + b + c) ^ 2 := by
  have middle := mul_nonpos_of_nonneg_of_nonpos (sub_nonneg.mpr hab)
    (show b + a - 4 * (a + c) ≤ 0 by linarith)
  have outer := mul_nonpos_of_nonneg_of_nonpos (show 0 ≤ c by linarith)
    (sub_nonpos.mpr hc)
  nlinarith only [middle, outer, sq_pos_of_pos ha]

/-- The local quadratic certificate for a positively weighted row bound.
The sum of these residual squares is nonnegative on every free form. -/
theorem weighted_edge_residual (a b x y : ℝ) :
    0 ≤ b ^ 2 * x ^ 2 + a ^ 2 * y ^ 2 - 2 * a * b * x * y := by
  nlinarith only [sq_nonneg (b * x - a * y)]

end OPH.NativeSpectroscopy
