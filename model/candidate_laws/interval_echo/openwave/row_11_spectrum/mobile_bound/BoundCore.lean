import Mathlib.Tactic

namespace OPH.BoundCore

/-- The energy envelope for two or more separated runs has this maximum at
two runs when the contrast coefficient exceeds one half. -/
theorem separated_ceiling (d runs κ : ℝ) (hc : 2 ≤ runs) (hk : 1/2 ≤ κ) :
    (4*d - 2*runs)*κ + d + runs ≤ (4*d - 4)*κ + d + 2 := by
  nlinarith [mul_nonneg (sub_nonneg.mpr hc) (show 0 ≤ 2*κ - 1 by linarith)]

/-- If a sharp core keeps its mean energy above the complementary block,
the normalized probability of an intact core has a uniform lower bound. -/
theorem all_time_core_weight (gap coupling p : ℝ) (hd : 0 < gap)
    (hp' : p ≤ 1)
    (energy : gap^2 * (1-p)^2 ≤ 4*coupling^2*p*(1-p)) :
    gap^2 / (gap^2 + 4*coupling^2) ≤ p := by
  have den : 0 < gap^2 + 4*coupling^2 := by positivity
  apply (div_le_iff₀ den).2
  by_cases h : p = 1
  · rw [h]
    nlinarith [sq_nonneg coupling]
  · have rem : 0 < 1-p := sub_pos.mpr (lt_of_le_of_ne hp' h)
    have diff : (gap^2*(1-p) - 4*coupling^2*p)*(1-p) ≤ 0 := by nlinarith
    have h' : gap^2*(1-p) - 4*coupling^2*p ≤ 0 := by
      by_contra! pos
      exact (not_lt_of_ge diff) (mul_pos pos rem)
    nlinarith

/-- The measured charge-core bound uses integer local contact counts. -/
theorem native_floors :
    (13:ℝ)^2/(13^2+4) = 169/173 ∧ (12:ℝ)^2/(12^2+6) = 24/25 := by
  norm_num

end OPH.BoundCore
