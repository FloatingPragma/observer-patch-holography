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

/-- One orbital interval at a positive scaled radius. Equal radii and repeated
measurements are allowed. Upper endpoints may be inconsistent with the lowers. -/
structure OrbitInterval where
  radius : ℝ
  lower : ℝ
  upper : ℝ
  radius_pos : 0 < radius
  lower_nonneg : 0 ≤ lower

/-- A nonnegative enclosed mass on positive radii, nondecreasing with radius
and with nonincreasing mass per unit radius. No differentiability is assumed. -/
structure OrbitalShape (m : ℝ → ℝ) : Prop where
  nonneg : ∀ x, 0 < x → 0 ≤ m x
  increasing : ∀ x y, 0 < x → x ≤ y → m x ≤ m y
  sublinear : ∀ x y, 0 < x → x ≤ y → m y / y ≤ m x / x

/-- Least possible mass forced at radius x by one orbital lower bound. -/
noncomputable def OrbitInterval.floorAt (b : OrbitInterval) (x : ℝ) : ℝ :=
  b.lower * min 1 (x / b.radius)

private theorem floor_shape (b : OrbitInterval) : OrbitalShape b.floorAt := by
  have perRadius (x : ℝ) (hx : 0 < x) :
      b.floorAt x / x = min (b.lower / x) (b.lower / b.radius) := by
    rw [OrbitInterval.floorAt, mul_min_of_nonneg _ _ b.lower_nonneg, mul_one,
      ← min_div_div_right hx.le]
    congr 1
    field_simp [hx.ne', b.radius_pos.ne']
  constructor
  · intro x hx
    exact mul_nonneg b.lower_nonneg (le_min zero_le_one (div_nonneg hx.le b.radius_pos.le))
  · intro x y _ hxy
    exact mul_le_mul_of_nonneg_left
      (min_le_min le_rfl (div_le_div_of_nonneg_right hxy b.radius_pos.le)) b.lower_nonneg
  · intro x y hx hxy
    rw [perRadius x hx, perRadius y (hx.trans_le hxy)]
    exact min_le_min (div_le_div_of_nonneg_left b.lower_nonneg hx hxy) le_rfl

private theorem floor_forced (b : OrbitInterval) {m : ℝ → ℝ} (h : OrbitalShape m)
    (hb : b.lower ≤ m b.radius) (x : ℝ) (hx : 0 < x) : b.floorAt x ≤ m x := by
  unfold OrbitInterval.floorAt
  rcases le_total b.radius x with hbx | hxb
  · rw [min_eq_left ((le_div_iff₀ b.radius_pos).2 (by simpa using hbx)), mul_one]
    exact hb.trans (h.increasing _ _ b.radius_pos hbx)
  · rw [min_eq_right ((div_le_iff₀ b.radius_pos).2 (by simpa using hxb)), ← mul_div_assoc]
    apply (div_le_iff₀ b.radius_pos).2
    exact (mul_le_mul_of_nonneg_right hb hx.le).trans
      ((div_le_div_iff₀ b.radius_pos hx).1 (h.sublinear _ _ hx hxb))

/-- The finite pointwise maximum of all truncated lower profiles, with the
zero profile for an empty list. -/
noncomputable def orbitalEnvelope : List OrbitInterval → ℝ → ℝ
  | [], _ => 0
  | b :: bs, x => max (b.floorAt x) (orbitalEnvelope bs x)

private theorem envelope_shape (bs : List OrbitInterval) : OrbitalShape (orbitalEnvelope bs) := by
  induction bs with
  | nil =>
    constructor <;> simp [orbitalEnvelope]
  | cons b bs ih =>
    have h := floor_shape b
    constructor
    · intro x hx
      exact (h.nonneg _ hx).trans (le_max_left _ _)
    · intro x y hx hxy
      exact max_le_max (h.increasing _ _ hx hxy) (ih.increasing _ _ hx hxy)
    · intro x y hx hxy
      simp only [orbitalEnvelope, ← max_div_div_right hx.le,
        ← max_div_div_right (hx.trans_le hxy).le]
      exact max_le_max (h.sublinear _ _ hx hxy) (ih.sublinear _ _ hx hxy)

private theorem envelope_le (bs : List OrbitInterval) (x u : ℝ) :
    orbitalEnvelope bs x ≤ u ↔ 0 ≤ u ∧ ∀ b ∈ bs, b.floorAt x ≤ u := by
  induction bs with
  | nil => simp [orbitalEnvelope]
  | cons b bs ih => simp [orbitalEnvelope, ih, and_left_comm]

private theorem member_le_envelope {bs : List OrbitInterval} {b : OrbitInterval}
    (hb : b ∈ bs) (x : ℝ) : b.floorAt x ≤ orbitalEnvelope bs x :=
  ((envelope_le bs x (orbitalEnvelope bs x)).1 le_rfl).2 b hb

/-- Finite orbital feasibility. Propagate every lower endpoint to every
measured radius. Exactly when all propagated lowers fit below the respective
upper endpoints, one scale-free mass profile fits every interval. The witness
is the finite maximum of the truncated linear profiles. -/
theorem finiteOrbitalFeasibility (bs : List OrbitInterval) :
    (∃ m : ℝ → ℝ, OrbitalShape m ∧
      ∀ b ∈ bs, b.lower ≤ m b.radius ∧ m b.radius ≤ b.upper) ↔
    ∀ i ∈ bs, ∀ j ∈ bs, i.lower * min 1 (j.radius / i.radius) ≤ j.upper := by
  constructor
  · rintro ⟨m, hm, fits⟩ i hi j hj
    exact (floor_forced i hm (fits i hi).1 j.radius j.radius_pos).trans (fits j hj).2
  · intro h
    refine ⟨orbitalEnvelope bs, envelope_shape bs, ?_⟩
    intro b hb
    have self : b.floorAt b.radius = b.lower := by
      simp [OrbitInterval.floorAt, div_self b.radius_pos.ne']
    have hbb : b.floorAt b.radius ≤ b.upper := h b hb b hb
    refine ⟨self ▸ member_le_envelope hb b.radius, (envelope_le bs _ _).2 ⟨?_, ?_⟩⟩
    · rw [self] at hbb
      exact b.lower_nonneg.trans hbb
    · intro i hi
      exact h i hi b hb

end OPH.EinsteinBranch.HaloStress
