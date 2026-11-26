import Mathlib

/-! Discrete convexity of any translated positive Klein–Gordon branch, including
an arbitrary linear drift. Set z = m + i c k and d = i c Δk. -/

namespace PolarBeat

theorem adj_secants (z d : ℂ) (offset drift : ℝ) :
    (offset + drift + ‖z + d‖) - (offset + ‖z‖) ≤
    (offset + 2 * drift + ‖z + 2 * d‖) - (offset + drift + ‖z + d‖) := by
  have h := norm_add_le z (z + 2 * d)
  have midpoint : z + (z + 2 * d) = 2 * (z + d) := by ring
  rw [midpoint, norm_mul] at h
  norm_num at h
  linarith

end PolarBeat
