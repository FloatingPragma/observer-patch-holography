import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Linarith

/-! # Request counts for uniform pair repair

The executable instrument supplies the success law.  These arithmetic
lemmas bound its sum of reciprocal success probabilities and evaluate the
finite product multiplying its slowest waiting-time tail.
-/
namespace OPH.RoutedAttemptBudget

/-- Mean, in units of the number of possible unordered port pairs. -/
noncomputable def waitSum (gap : ℕ) : ℕ → ℝ
  | 0 => 0
  | n+1 => waitSum gap n + 1 / ((n+1 : ℝ) * (n+1+(gap : ℝ)))

theorem wait_le_balanced (gap n : ℕ) : waitSum gap n ≤ waitSum 0 n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [waitSum, Nat.cast_zero, add_zero]
    apply add_le_add ih
    apply one_div_le_one_div_of_le (by positivity)
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    linarith [show 0 ≤ (gap : ℝ) by positivity]

theorem balanced_bound (n : ℕ) : waitSum 0 (n+1) ≤ 2 - 1/(n+1 : ℝ) := by
  induction n with
  | zero => norm_num [waitSum]
  | succ n ih =>
    rw [waitSum]
    have x : 0 < (n : ℝ)+1 := by positivity
    have y : 0 < (n : ℝ)+1+1 := by positivity
    push_cast
    simp only [add_zero]
    apply (add_le_add_left ih (1 / (((n : ℝ)+1+1)*((n : ℝ)+1+1)))).trans
    field_simp
    nlinarith

theorem mean_lt_two (gap n : ℕ) : waitSum gap n < 2 := by
  cases n with
  | zero => norm_num [waitSum]
  | succ n =>
    have h := (wait_le_balanced gap (n+1)).trans (balanced_bound n)
    have pos : 0 < 1 / ((n : ℝ)+1) := by positivity
    linarith

/-- Slowest-rate geometric moment factor for stages 2 through n+1. -/
noncomputable def momentProduct (gap : ℕ) : ℕ → ℝ
  | 0 => 1
  | n+1 => momentProduct gap n *
    (((n+2 : ℝ) * (n+2+(gap : ℝ))) / ((n+1 : ℝ) * (n+3+(gap : ℝ))))

theorem moment_telescopes (gap n : ℕ) :
    momentProduct gap n = ((n+1 : ℝ) * ((gap : ℝ)+2)) / (n+2+(gap : ℝ)) := by
  induction n with
  | zero => simp only [momentProduct, Nat.cast_zero, zero_add, one_mul]
            field_simp
            ring
  | succ n ih =>
    rw [momentProduct, ih]
    push_cast
    field_simp
    ring

end OPH.RoutedAttemptBudget
