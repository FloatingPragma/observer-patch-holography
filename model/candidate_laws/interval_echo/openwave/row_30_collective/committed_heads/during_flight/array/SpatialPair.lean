import Mathlib.Tactic.Ring

namespace OPH.WallArray

/-- A wave and its reflected branch share one second-order recurrence,
    for arbitrary complex amplitudes when the two phases are reciprocal. -/
theorem common_spatial_recurrence {R : Type*} [CommRing R]
    (r s a b : R) (h : r * s = 1) (j : ℕ) :
    (a * r^(j+2) + b * s^(j+2)) + (a * r^j + b * s^j) =
      (r+s) * (a * r^(j+1) + b * s^(j+1)) := by
  have h₁ : r*r + 1 = (r+s)*r := by rw [← h]; ring
  have h₂ : s*s + 1 = (r+s)*s := by rw [← h]; ring
  calc
    _ = (a*r^j)*(r*r+1) + (b*s^j)*(s*s+1) := by simp only [pow_succ]; ring
    _ = _ := by rw [h₁, h₂]; simp only [pow_succ]; ring

end OPH.WallArray
