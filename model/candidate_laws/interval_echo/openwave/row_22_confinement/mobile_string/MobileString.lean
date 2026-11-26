import Mathlib.Tactic.Linarith
import Mathlib.Data.Real.Basic

namespace IntervalEcho.MobileString

/-- For any normalized interior trial state, changing both fixed spins changes
only its two boundary expectations. This applies to arbitrary entangled trials
and also to the one-site interior, whose two bounding contacts share a spin.
Energies are relative to the two common outer lower pins. -/
theorem work_range {S : Type*} (energy left right : S → ℝ)
    (minus plus : S) (work : ℝ)
    (leftBound : ∀ s, -1 ≤ left s ∧ left s ≤ 1)
    (rightBound : ∀ s, -1 ≤ right s ∧ right s ≤ 1)
    (minLower : ∀ s, energy minus ≤ energy s)
    (minUpper : ∀ s, energy plus - 16 * (left plus + right plus) ≤
      energy s - 16 * (left s + right s))
    (read : work = 48 + energy plus - 16 * (left plus + right plus) - energy minus) :
    16 ≤ work ∧ work ≤ 80 := by
  constructor
  · linarith [(leftBound plus).2, (rightBound plus).2, minLower plus]
  · linarith [(leftBound minus).1, (rightBound minus).1, minUpper minus]

/-- At zero transverse exchange a sharp intact word costs 32 + 8*l.
The least screened word costs 80 when there is room for a third site. -/
theorem sharp_crossover (length : ℝ) :
    (length ≤ 6 → min (32 + 8*length) 80 = 32 + 8*length) ∧
    (6 ≤ length → min (32 + 8*length) 80 = 80) := by
  constructor
  · intro h
    apply min_eq_left
    linarith
  · intro h
    apply min_eq_right
    linarith

end IntervalEcho.MobileString
