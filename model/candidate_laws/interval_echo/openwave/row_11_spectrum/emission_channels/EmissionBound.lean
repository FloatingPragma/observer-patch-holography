import Mathlib

namespace SourceLaw.FreeEmission

/-- The probability amplitude in a resolved daughter chamber is bounded at every
free exposure by the parent energy spread and the leakage of that chamber.
The gap hypothesis is on the compressed daughter Hamiltonian. -/
theorem finite_width
    {E F : Type*} [NormedAddCommGroup E] [NormedSpace ℂ E]
    [NormedAddCommGroup F] [NormedSpace ℂ F]
    (H U : E →ₗ[ℂ] E) (P : E →ₗ[ℂ] F) (K : F →ₗ[ℂ] F)
    (v : E) (energy : ℂ) (gap spread leakage : ℝ)
    (positive : 0 < gap) (commute : U (H v) = H (U v))
    (isometry : ∀ w, ‖U w‖ = ‖w‖) (contract : ∀ w, ‖P w‖ ≤ ‖w‖)
    (parent : ‖energy • v - H v‖ ≤ spread)
    (escape : ‖P (H (U v)) - K (P (U v))‖ ≤ leakage)
    (separation : ∀ y, gap * ‖y‖ ≤ ‖energy • y - K y‖) :
    ‖P (U v)‖ ≤ (spread + leakage) / gap := by
  have decompose : energy • P (U v) - K (P (U v)) =
      P (U (energy • v - H v)) + (P (H (U v)) - K (P (U v))) := by
    simp only [map_sub, map_smul, commute]
    abel
  apply (le_div_iff₀ positive).2
  rw [mul_comm]
  calc
    gap * ‖P (U v)‖ ≤ ‖energy • P (U v) - K (P (U v))‖ := separation _
    _ ≤ ‖P (U (energy • v - H v))‖ + ‖P (H (U v)) - K (P (U v))‖ := by
      rw [decompose]
      exact norm_add_le _ _
    _ ≤ ‖U (energy • v - H v)‖ + leakage := add_le_add (contract _) escape
    _ ≤ spread + leakage := by rw [isometry]; exact add_le_add parent le_rfl

end SourceLaw.FreeEmission
