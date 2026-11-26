import Mathlib.Tactic

/-! Randomized dark, bright and signal clicks in a four-arm receiver. -/

namespace IntervalEcho.DriftingReceiver

/-- Two response endpoints give the click rate at ideal bright-exit probability p. -/
def rate (dark bright p : ℝ) : ℝ := (1 - p) * dark + p * bright

/-- The conditional score for the true exit has zero mean for each possible
previous history, even if both response endpoints change with that history. -/
theorem concealedScore (dark bright p : ℝ) :
    rate dark bright p / 4 - (1 - p) * dark / 4 - p * bright / 4 = 0 := by
  unfold rate
  ring

/-- Inversion uses the observed dark/bright span. The counts may be arbitrary. -/
theorem invertScore {dark bright signal p allowance : ℝ} (h : 0 < bright - dark) :
    |signal - dark - p * (bright - dark)| ≤ allowance ↔
    (signal - dark - allowance) / (bright - dark) ≤ p ∧
      p ≤ (signal - dark + allowance) / (bright - dark) := by
  rw [abs_le, div_le_iff₀ h, le_div_iff₀ h]
  constructor <;> rintro ⟨hl, hu⟩ <;> constructor <;> linarith

/-- A source-label offset can counterfeit a distinct true source exit with
identical dark, bright and half controls. Both click rates are admissible. -/
theorem exposedLabel :
    rate (1/8) (5/8) (1/4) = rate (1/8) (5/8) (1/8) + 1/16 := by
  norm_num [rate]

end IntervalEcho.DriftingReceiver

