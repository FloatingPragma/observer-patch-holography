import Mathlib.Tactic

namespace OPH.PooledFeedback

/-- Both published intervals cover the same surviving ensemble mean. -/
theorem common_interval {m a b r : ℝ} (hm : m ∈ Set.Icc (-1) 1)
    (ha : |a-m| ≤ r) (hb : |b-m| ≤ r) :
    m ∈ Set.Icc (max (-1) (max a b - r)) (min 1 (min a b + r)) := by
  rcases hm with ⟨hlo, hhi⟩
  rcases abs_le.mp ha with ⟨ha₁, ha₂⟩
  rcases abs_le.mp hb with ⟨hb₁, hb₂⟩
  constructor
  · apply max_le hlo
    rw [sub_le_iff_le_add]
    exact max_le (by linarith) (by linarith)
  · apply le_min hhi
    rw [← sub_le_iff_le_add]
    exact le_min (by linarith) (by linarith)

/-- An exterior sign cannot be reversed on the covered event. -/
theorem sign_sound {m a r : ℝ} (ha : |a-m| ≤ r) :
    (r < a → 0 < m) ∧ (a < -r → m < 0) := by
  rcases abs_le.mp ha with ⟨h₁, h₂⟩
  constructor <;> intro h <;> linarith

theorem separated {m a b r : ℝ} (ha : |a-m| ≤ r) (hb : |b-m| ≤ r) :
    (2*r < m → r < a ∧ r < b) ∧ (m < -2*r → a < -r ∧ b < -r) := by
  rcases abs_le.mp ha with ⟨ha₁, ha₂⟩
  rcases abs_le.mp hb with ⟨hb₁, hb₂⟩
  constructor <;> intro h <;> constructor <;> linarith

end OPH.PooledFeedback
