import Mathlib.Tactic

set_option autoImplicit false

namespace OPH.QFT.ElectricCutoffFloor

/-- Completing the square controls the omitted electric sector. The source
must supply the two diagonal lower bounds and the off-diagonal norm bound. -/
theorem two_block_rayleigh_lower (L b c u v : ℝ) (hgap : L < b) :
    (L - c^2 / (b-L)) * (u^2+v^2) ≤ L*u^2+b*v^2-2*c*u*v := by
  have hg : 0 < b-L := sub_pos.mpr hgap
  have hs : 0 ≤ ((b-L)*v-c*u)^2+c^2*v^2 :=
    add_nonneg (sq_nonneg _) (mul_nonneg (sq_nonneg _) (sq_nonneg _))
  have hid : (L*u^2+b*v^2-2*c*u*v) -
      (L-c^2/(b-L))*(u^2+v^2) =
      (((b-L)*v-c*u)^2+c^2*v^2)/(b-L) := by
    field_simp
    ring
  have := div_nonneg hs hg.le
  rw [← hid] at this
  linarith

theorem normalized_rayleigh_lower (L b c u v : ℝ) (hgap : L < b)
    (hnorm : u^2+v^2=1) :
    L-c^2/(b-L) ≤ L*u^2+b*v^2-2*c*u*v := by
  simpa [hnorm] using two_block_rayleigh_lower L b c u v hgap

end OPH.QFT.ElectricCutoffFloor
