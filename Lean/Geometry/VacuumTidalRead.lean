import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring

/-!
# What two vacuum polarizations look like to curvature readers

Work in a supplied linearized Minkowski chart of signature (-,+,+,+).
The covector is (s,0,0,-s). The amplitudes include the common second
derivative of the wave profile at the read event. These are classical
curvature observations in the chart.
-/

namespace OPH.VacuumTidalRead

inductive Axis | time | x | y | z

open Axis

noncomputable section

def wave (s : ℝ) : Axis → ℝ
  | time => s
  | z => -s
  | _ => 0

def amplitude (plus cross : ℝ) : Axis → Axis → ℝ
  | x, x => plus
  | y, y => -plus
  | x, y | y, x => cross
  | _, _ => 0

/-- Twice the linearized Riemann coefficient of a plane perturbation. -/
def curvature (s plus cross : ℝ) (a b c d : Axis) : ℝ :=
  let k := wave s
  let h := amplitude plus cross
  h a d * k b * k c + h b c * k a * k d -
    h a c * k b * k d - h b d * k a * k c

/-- Contract the first and third indices in the stated signature. -/
def twiceRicci (s plus cross : ℝ) (b d : Axis) : ℝ :=
  let r := curvature s plus cross
  (-(r time b time d)) + r x b x d + r y b y d + r z b z d

theorem ricci_vanishes (s plus cross : ℝ) (b d : Axis) :
    twiceRicci s plus cross b d = 0 := by
  cases b <;> cases d <;> simp [twiceRicci, curvature, amplitude, wave]

/-- Normalized arm accelerations for the x arm and the 45-degree x+y arm. -/
def armRead (s plus cross : ℝ) : ℝ × ℝ :=
  let r := curvature s plus cross
  (-r time x time x,
   -(r time x time x + r time x time y + r time y time x + r time y time y) / 2)

theorem arm_response (s plus cross : ℝ) :
    armRead s plus cross = (s ^ 2 * plus, s ^ 2 * cross) := by
  apply Prod.ext <;> simp [armRead, curvature, amplitude, wave] <;> ring

theorem arm_separates {s : ℝ} (hs : s ≠ 0) {p c p' c' : ℝ} :
    armRead s p c = armRead s p' c' ↔ p = p' ∧ c = c' := by
  rw [arm_response, arm_response, Prod.mk.injEq]
  have h₂ : s ^ 2 ≠ 0 := pow_ne_zero 2 hs
  simp [h₂]

end

end OPH.VacuumTidalRead
