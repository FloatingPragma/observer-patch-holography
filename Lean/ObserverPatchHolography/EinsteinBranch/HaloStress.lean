import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-! # Stress and circular speed on a static spherical exterior

Use geometric units and orthonormal stresses. The mass equation fixes the
density. Its radial Einstein equation and the radial conservation equation
are hypotheses about a metric and effective source on an annulus.
-/

namespace OPH.EinsteinBranch.HaloStress

/-- At vanishing radial pressure, the radial Einstein equation fixes the
square of the local circular speed. The hydrostatic equation then fixes the
tangential pressure. This statement concerns a static source. -/
theorem tangentialClosure (r m betaSq rho pTan : ℝ) (hStatic : r - 2*m ≠ 0)
    (hRadial : (r - 2*m)*betaSq = m)
    (hConserved : rho*betaSq = 2*pTan) :
    betaSq = m/(r-2*m) ∧ pTan = rho*betaSq/2 := by
  constructor
  · apply (eq_div_iff hStatic).2
    simpa [mul_comm] using hRadial
  · linarith

/-- With a fixed linear mass profile m = h r and radial pressure w rho
in the scale-free annulus, the square of the circular speed depends on w.
The density profile alone cannot pick a unique orbital response. -/
theorem pressureChangesOrbit (h w₁ w₂ x₁ x₂ : ℝ) (hh : h ≠ 0)
    (hdiff : w₁ ≠ w₂)
    (h₁ : (1-2*h)*x₁ = h*(1+w₁))
    (h₂ : (1-2*h)*x₂ = h*(1+w₂)) : x₁ ≠ x₂ := by
  intro he
  have ht : h * (1+w₁) = h * (1+w₂) := by rw [← h₁, ← h₂, he]
  exact hdiff (by linarith [mul_left_cancel₀ hh ht])

/-- The signed shear of a spherical thin shell outside the aperture, in
units of its mass divided by the aperture area. Here c is the absolute
cosine of the line-of-sight angle of its intersection with the rim. -/
theorem exteriorShellContribution (c : ℝ) (hc : 0 < c) :
    1-c-(1-c^2)/(2*c) = -(1-c)^2/(2*c) ∧ 1-c-(1-c^2)/(2*c) ≤ 0 := by
  have he : 1-c-(1-c^2)/(2*c) = -(1-c)^2/(2*c) := by
    field_simp
    ring
  refine ⟨he, ?_⟩
  rw [he]
  exact div_nonpos_of_nonpos_of_nonneg (neg_nonpos.mpr (sq_nonneg _)) (by linarith)

/-- One upper measurement of normalized halo mass bounds every point of a
common monotone profile whose growth per unit scaled radius cannot increase. -/
theorem massTransfer (f : ℝ → ℝ) (x d cap : ℝ) (hx : 0 < x) (hd : 0 < d)
    (hCap : f d ≤ cap) (hm : x ≤ d → f x ≤ f d)
    (hs : d ≤ x → f x / x ≤ f d / d) :
    f x ≤ if x ≤ d then cap else cap * (x/d) := by
  split_ifs with h
  · exact (hm h).trans hCap
  · have h₁ := (div_le_div_iff₀ hx hd).mp (hs (le_of_lt (lt_of_not_ge h)))
    have h₂ := mul_le_mul_of_nonneg_right hCap (le_of_lt hx)
    have h₃ : f x * d ≤ cap * x := h₁.trans h₂
    rw [← mul_div_assoc]
    exact (le_div_iff₀ hd).mpr h₃

end OPH.EinsteinBranch.HaloStress
