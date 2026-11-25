import Mathlib

namespace SourceLaw.ResonantTransfer

/-- The receiver ground includes a parked compensator. The opposite local
transitions have exactly the same sum for any common splitting. -/
theorem energy_shell (ω : ℝ) :
    (41 + ω) / 2 + (89 - ω) / 2 = 65 ∧
    (41 - ω) / 2 + (89 + ω) / 2 = 65 := by constructor <;> ring

/-- The local number projector has transition magnitude 1/(2ω) on each
arm. The actual stored contact has nonzero matrix element on the shell. -/
theorem native_on_shell (ω : ℝ) (gap : ω^2 = 1601) :
    -32 * (1 / (2 * ω))^2 = -(8 / 1601) ∧
    (-32 * (1 / (2 * ω))^2) * ω = -8 / ω := by
  have nonzero : ω ≠ 0 := by intro e; rw [e] at gap; norm_num at gap
  constructor <;> field_simp [nonzero] <;> nlinarith

/-- With the matching cut the calibrated positive/negative local energy
levels are intact, while the off-diagonal product coupling vanishes. -/
theorem matched_cut (ω : ℝ) : (0 : ℝ) * (1 / (2 * ω))^2 = 0 := by simp

end SourceLaw.ResonantTransfer
