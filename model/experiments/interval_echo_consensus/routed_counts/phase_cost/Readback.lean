import Mathlib

namespace OPH.TwelvePortPhaseCost

/-- Detector A clicks on a single, and B clicks on a double or with rate a on a single.
The two pooled counts (here in units of one port) determine the square sum. -/
theorem decode (singles doubles a : ℝ) :
    singles + 4 * ((a * singles + doubles) - a * singles) = singles + 4 * doubles := by
  ring

/-- In the total-twelve sector, the full quadratic reading is twelve plus twice
the double count. Occupancy and total load then also determine the vacancies. -/
theorem floor_gap (zero one two : ℕ) (ports : zero + one + two = 12)
    (load : one + 2 * two = 12) :
    zero = two ∧ one + 4 * two = 12 + 2 * two ∧
    (one + 4 * two = 12 ↔ one = 12) ∧
    (one ≠ 12 → 14 ≤ one + 4 * two) := by omega

theorem conservative_write (one two : ℕ) :
    (one + 2) + 4 * two + 2 = one + 4 * (two + 1) := by omega

end OPH.TwelvePortPhaseCost
