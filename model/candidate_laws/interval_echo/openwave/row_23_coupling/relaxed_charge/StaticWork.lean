import Mathlib.Tactic

namespace IntervalEcho.RelaxedCharge

/-- One shared vacancy chooses the cheaper well. Its signed mixed energy is
the length of the overlap of the two wells' charge-shifted energy intervals. -/
theorem sharedVacancyWork (a b c d : ℝ) (hab : a ≤ b) (hcd : c ≤ d) :
    min b d - min b c - min a d + min a c = max 0 (min b d - max a c) := by
  simp only [min_def, max_def]
  split_ifs <;> linarith

/-- Connected signed work cancels for additive, separately stocked regions. -/
theorem independentStockWork (a b c d : ℝ) :
    (b + d) - (b + c) - (a + d) + (a + c) = 0 := by ring

/-- The sign contrast of a closed thermal walk factorizes at the two pins. -/
theorem walkSignFactor (a b c d : ℝ) :
    b * d - b * c - a * d + a * c = (b-a) * (d-c) := by ring

/-- Visits to both pins have positive signed weight for ordered contact gains. -/
theorem visitedWalkWeight (a b c d : ℝ) (hab : b < a) (hcd : d < c) :
    0 < b * d - b * c - a * d + a * c := by
  rw [walkSignFactor]
  exact mul_pos_of_neg_of_neg (sub_neg.mpr hab) (sub_neg.mpr hcd)

end IntervalEcho.RelaxedCharge
