import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring

namespace OPH.NativeOccupationAction

/-- The four-corner action of a symmetric binary contact, with its declared
velocity-only continuation. Here s is the measured log stay/flip ratio. -/
noncomputable def lift (c s a x y : ℝ) : ℝ :=
  c + s * (x + y - 2*x*y) + (a/2)*y*(y-1)

/-- Difference from the translation-invariant two-point action is a boundary
term at the declared mode-extremal curvature. -/
theorem boundary_term (c s x y : ℝ) :
    lift c s (4*s) x y =
      c + s*(y-x)^2 + (s*y*(y-1) - s*x*(x-1)) := by
  unfold lift
  ring

/-- Every real replacement at an interior point of a constant path has the
displayed gap, including both embedded constant modes. -/
theorem constant_path_variation (c s t z : ℝ) :
    lift c s (4*s) t z + lift c s (4*s) z t - 2*lift c s (4*s) t t =
      2*s*(z-t)^2 := by
  unfold lift
  ring

/-- At the initial-zero edge every single quadratic curvature predicts zero
in the held quarter-minus-three-quarter-halves comparison. -/
theorem quarter_half_quadratic (c s a : ℝ) :
    (lift c s 0 0 (1/4) - lift c s a 0 (1/4)) -
    (3/4)*(lift c s 0 0 (1/2) - lift c s a 0 (1/2)) = 0 := by
  unfold lift
  ring

end OPH.NativeOccupationAction
