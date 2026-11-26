import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring

namespace EnergySource

/-- The antisymmetric click on one link, given its two populations and
the real part of their coherence in the native position basis. -/
noncomputable def odd (a b z : ℝ) : ℝ := (a+b)/2-z

/-- The four-link sector has diagonal one and hopping minus one half above
its ground energy. Each undirected coherence occurs twice in the trace. -/
theorem four_link_aperture (a b c d x y z w : ℝ) :
    (odd a b x + odd b c y + odd c d z + odd d a w)/4 =
      ((a+b+c+d)-(x+y+z+w))/4 := by unfold odd; ring

/-- Composition of a source click with an independent terminated guide and
its common, calibrated inward response; silent source trials are counted. -/
theorem response_per_preparation (energy current p : ℝ) (h : p=energy/4) :
    p*(-current)+(1-p)*0 = -energy*current/4 := by rw [h]; ring

end EnergySource
