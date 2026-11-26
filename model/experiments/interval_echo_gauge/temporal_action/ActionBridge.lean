import ScaledMaxwellStability

namespace OPH.TemporalActionReadout

open OPH.DiscreteCoulombGreen OPH.PositionSpaceMaxwellAction OPH.LocalFaceMaxwellAction
open OPH.ScaledMaxwellStability

/-- The four addressed contributions have precisely the deposited signs and forward times. -/
theorem half_step_polynomial (A J : ℕ → Fin 30 → ℝ) (φ ρ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    4 * stepLagrangian (1/2) A φ ρ J n =
      realSeamEnergy (electricFieldScaled (1/2) A φ n)
      - faceEnergy (faceCurvature (A (n+1)))
      + 2 * realSeamInner (J n) (A (n+1)) + 2 * realPortInner (ρ n) (φ n) := by
  unfold stepLagrangian localSourcedAction
  ring

/-- Symmetric one-word variation of a weighted local square. -/
theorem two_path_difference (weight val delta : ℤ) :
    weight * (val + delta)^2 - weight * (val - delta)^2 = 4 * weight * val * delta := by
  ring

#print axioms half_step_polynomial
#print axioms two_path_difference

end OPH.TemporalActionReadout
