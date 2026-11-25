import Mathlib.Data.Complex.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.LinearCombination

/-! One return amplitude after applying a common phase to the marked output
of a unitary dilation, undoing the dilation and heralding its input space.
The executable pulse columns supply the effect eigenvalue on each character. -/
namespace OPH.LeptonRegisterFeedback

noncomputable def tick (p : ℝ) (z : ℂ) : ℂ := 1 + (z - 1) * (p : ℂ)

@[simp] theorem tick_re (p : ℝ) (z : ℂ) : (tick p z).re = 1 + (z.re - 1) * p := by
  simp [tick]

@[simp] theorem tick_im (p : ℝ) (z : ℂ) : (tick p z).im = z.im * p := by
  simp [tick]

theorem single_dwell_identifies (p q : ℝ) (z : ℂ) (cal : z.im ≠ 0)
    (same : (tick p z).im = (tick q z).im) : p = q := by
  simpa [tick_im, cal] using same

theorem retention (p : ℝ) (z : ℂ) (unit : z.re^2 + z.im^2 = 1) :
    (tick p z).re^2 + (tick p z).im^2 = 1 - 2*p*(1-p)*(1-z.re) := by
  rw [tick_re, tick_im]
  linear_combination p^2 * unit

theorem half_cycle_erases_phase (p : ℝ) (small : p < 1/2) :
    (tick p (-1)).im = 0 ∧ 0 < (tick p (-1)).re := by
  constructor
  · simp
  · simp only [tick_re, Complex.neg_re, Complex.one_re]
    linarith

theorem commuting_character_ticks (p : ℝ) (x y : ℂ) :
    tick p x * tick p y = tick p y * tick p x := mul_comm _ _

end OPH.LeptonRegisterFeedback
