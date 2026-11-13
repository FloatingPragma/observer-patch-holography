import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# The actual four-corner record transport bound

The retained tensor decoder acts on co-located records. These inequalities
bound the additional error when a staggered component is transported by its
four-corner mean. They assume measured finite neighbour variations, not a
continuum metric, physical stress, or a curvature-acquisition instrument.
-/

namespace OPH.TensorRecordTransport

theorem four_corner_error {a b c d dx dy : ℝ}
    (hab : |a - b| ≤ dx) (hac : |a - c| ≤ dy)
    (had : |a - d| ≤ dx + dy) :
    |a - (a + b + c + d) / 4| ≤ (dx + dy) / 2 := by
  rcases abs_le.mp hab with ⟨hab₁, hab₂⟩
  rcases abs_le.mp hac with ⟨hac₁, hac₂⟩
  rcases abs_le.mp had with ⟨had₁, had₂⟩
  apply abs_le.mpr
  constructor <;> linarith

theorem transported_remainder {b c target ε τ : ℝ}
    (htransport : |b - c| ≤ τ) (hdecoder : |c - target| ≤ ε) :
    |b - target| ≤ ε + τ := by
  rcases abs_le.mp htransport with ⟨ht₁, ht₂⟩
  rcases abs_le.mp hdecoder with ⟨hd₁, hd₂⟩
  apply abs_le.mpr
  constructor <;> linarith

#print axioms four_corner_error
#print axioms transported_remainder

end OPH.TensorRecordTransport
