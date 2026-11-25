import Mathlib.Tactic

namespace IntervalEcho.FiniteDensity

/-- An attractive inverse-distance energy has negative second difference even
    when its zero of energy and distance origin are free calibrations. -/
theorem attractive_inverse_curvature (u C r : ℝ) (hC : 0 < C) (hr : 0 < r) :
    (u - C / r) - 2 * (u - C / (r+1)) + (u - C / (r+2)) < 0 := by
  have h₁ : 0 < r+1 := by linarith
  have h₂ : 0 < r+2 := by linarith
  have he : (u - C/r) - 2*(u - C/(r+1)) + (u - C/(r+2)) =
      -(2*C)/(r*(r+1)*(r+2)) := by
    field_simp [ne_of_gt hr, ne_of_gt h₁, ne_of_gt h₂]
    ring
  rw [he]
  exact div_neg_of_neg_of_pos (by linarith) (by positivity)

/-- The connected partition odds of independent halves equal one. -/
theorem cut_odds (a₀ a₁ b₀ b₁ : ℝ) (ha₀ : 0 < a₀) (ha₁ : 0 < a₁)
    (hb₀ : 0 < b₀) (hb₁ : 0 < b₁) :
    ((a₁*b₁)*(a₀*b₀))/((a₁*b₀)*(a₀*b₁)) = 1 := by
  field_simp

/-- Fixing a nontrivial population correlates two distinct sites already
    in a uniform zero-interaction medium. -/
theorem fixed_population_covariance (L N : ℝ) (hs : 1 < L) (hn : 0 < N)
    (hv : N < L) :
    N*(N-1)/(L*(L-1)) - (N/L)^2 < 0 := by
  have hL : 0 < L := by linarith
  have hd : 0 < L-1 := by linarith
  have he : N*(N-1)/(L*(L-1)) - (N/L)^2 = -N*(L-N)/(L^2*(L-1)) := by
    field_simp [ne_of_gt hL, ne_of_gt hd]
    ring
  rw [he]
  exact div_neg_of_neg_of_pos (mul_neg_of_neg_of_pos (by linarith) (by linarith)) (by positivity)

end IntervalEcho.FiniteDensity
