import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Tactic

namespace IntervalEcho.Publication

/-- The multiplier of any single matrix unit for a freshly flagged signed read. -/
def tick (noFlag : ℝ) (same : Bool) : ℝ :=
  noFlag + (1 - noFlag) * (if same then 1 else 0)

theorem common_rate (rate dt : ℝ) (n : ℕ) (same : Bool) :
    (tick (Real.exp (-(rate * dt))) same) ^ n =
      if same then 1 else Real.exp ((n : ℝ) * (-(rate * dt))) := by
  cases same
  · simp only [tick, Bool.false_eq_true, if_false, mul_zero, add_zero]
    exact (Real.exp_nat_mul _ _).symm
  · simp [tick]

theorem two_collars (a b t : ℝ) :
    Real.exp (-(a * t)) * Real.exp (-(b * t)) = Real.exp (-((a + b) * t)) := by
  rw [← Real.exp_add]
  congr 1
  ring

/-- Every positive-rate pinching of a balanced pair loses linearly near zero.
A certainly bright finite Hamiltonian preparation with quadratic loss cannot
give this reading at a sufficiently short nonzero wait. -/
theorem short_time_exclusion (g t B : ℝ) (hg : 0 < g) (ht : 0 < t)
    (small : g * t ≤ 1) (speed : B^2 * t < g / 4) :
    B^2 * t^2 < (1 - Real.exp (-(g*t))) / 2 := by
  have hx : 0 < g * t := mul_pos hg ht
  have hx1 : 0 < 1 + g * t := by linarith
  have he : 1 + g * t ≤ Real.exp (g*t) := by linarith [Real.add_one_le_exp (g*t)]
  have hi := one_div_le_one_div_of_le hx1 he
  rw [Real.exp_neg]
  have ha : (1 / (1 + g*t)) * (1 + g*t) = 1 := by field_simp
  have hb : 0 ≤ g*t*(1-g*t) := mul_nonneg hx.le (by linarith)
  have hc : 1 / (1 + g*t) ≤ 1 - g*t / 2 := by
    apply (div_le_iff₀ hx1).mpr
    nlinarith
  have shrink := mul_lt_mul_of_pos_right speed ht
  rw [← one_div]
  nlinarith

end IntervalEcho.Publication
