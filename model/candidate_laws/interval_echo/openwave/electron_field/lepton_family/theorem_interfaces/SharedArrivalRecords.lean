import Mathlib.Data.Real.Basic
import Mathlib.Algebra.Order.Field.Basic
import Mathlib.Tactic.Linarith

/-! How accurately two accepted counts on all original starts must agree
to transfer a conditional binary detector prediction. -/
namespace OPH.SharedArrivalRecords

private theorem scaled_difference (a b x y ε : ℝ)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 1)
    (click : a * x - b * y ≤ ε)
    (quiet : b * (1 - y) - a * (1 - x) ≤ ε) : a * (x - y) ≤ ε := by
  rcases le_total b a with h | h
  · nlinarith [mul_nonneg (sub_nonneg.mpr h) hy₀]
  · nlinarith [mul_nonneg (sub_nonneg.mpr h) (sub_nonneg.mpr hy₁)]

/-- The reference admits fraction a, and the candidate admits fraction b.
Both have nonzero observed admission so their conditional read rates x and
y are defined. If the two original-start bins (arrival with a click,
arrival without a click) agree to ε apiece, the conditional error is at
most ε/a. The smaller candidate admission needs no separate divisor. -/
theorem conditional_read (a b x y ε : ℝ) (ha : 0 < a)
    (hy₀ : 0 ≤ y) (hy₁ : y ≤ 1)
    (click : |a * x - b * y| ≤ ε)
    (quiet : |a * (1 - x) - b * (1 - y)| ≤ ε) :
    |x - y| ≤ ε / a := by
  rcases abs_le.mp click with ⟨hcl, hcu⟩
  rcases abs_le.mp quiet with ⟨hql, hqu⟩
  have upper := scaled_difference a b x y ε hy₀ hy₁ hcu (by linarith)
  have lower := scaled_difference a b (1 - x) (1 - y) ε (by linarith) (by linarith) hqu
    (by simpa only [sub_sub_cancel] using (show b * y - a * x ≤ ε by linarith))
  apply abs_le.mpr
  constructor
  · have h : -(x - y) ≤ ε / a := (le_div_iff₀ ha).mpr (by nlinarith)
    linarith
  · exact (le_div_iff₀ ha).mpr (by nlinarith)

/-- With ±1 settings, an error ε in a binary agreement probability is an
error 2ε in the signed correlation. Four local comparisons must account
for the full gap between the reference Bell value and the local bound. -/
theorem four_setting_gap (r c : Fin 4 → ℝ) (ε : ℝ)
    (localBound : c 0 + c 1 + c 2 - c 3 ≤ 2)
    (records : ∀ j, |r j - c j| ≤ 2 * ε) :
    (r 0 + r 1 + r 2 - r 3 - 2) / 8 ≤ ε := by
  linarith [(abs_le.mp (records 0)).2, (abs_le.mp (records 1)).2,
    (abs_le.mp (records 2)).2, (abs_le.mp (records 3)).1]

end OPH.SharedArrivalRecords
