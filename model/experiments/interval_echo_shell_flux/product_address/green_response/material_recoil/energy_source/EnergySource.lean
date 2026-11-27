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

namespace EnergySource

/-- The sharp parked energy excludes the shared half unit read by the odd meter.
Together the effects give both native diagonals and both off-diagonal entries. -/
theorem two_material_energy_read (p q eNear eFar coherence : ℝ) :
    ((eNear-1/2)*p + (eFar-1/2)*q + odd p q coherence)/512 =
      (eNear*p + eFar*q - coherence)/512 := by unfold odd; ring

/-- Each side of a close pair loses the contact it supplied against the vacuum.
The comparison keeps both materials and changes only their relative placement. -/
theorem binding_contrast (k p q mismatch : ℝ) :
    (k*p^2+1/2)+(k*q^2+1/2)-(k*(p-q)^2+mismatch/2) =
    2*k*p*q+1-mismatch/2 := by ring

/-- Use the same negative signed displacement per emitted pulse on both
preparations. Binding reduces the magnitude of their inward response. -/
theorem recoil_from_binding (eNear eFar : ℝ) :
    -(eNear / 512) / 15 - (-(eFar / 512) / 15) = (eFar-eNear)/7680 := by ring

end EnergySource
