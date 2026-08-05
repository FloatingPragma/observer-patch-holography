import Mathlib

/-!
# The C3-circulant Koide identity

This module formalizes the exact algebra behind the conditional charged-family
circulant receipt.  If three cosine coordinates have sum zero and square sum
`3/2`, the three signed eigenvalues

`lambda_k = a + 2 * b * c_k`

have sum `3a` and square sum `3a^2 + 6b^2`.  Their signed Koide quotient is
therefore

`1/3 + (2/3) * (b/a)^2`.

For a nonnegative modulus ratio, this quotient equals `2/3` exactly when the
ratio is `1/sqrt(2)`.

The physical square-root-mass reading additionally requires every eigenvalue
to be nonnegative.  Outside that chamber physical square roots use absolute
values, so this module deliberately makes no physical-mass or phase claim.
It also does not formalize the open source-to-charged-family attachment.
-/

namespace OPH.KoideCirculant

/-- The three signed circulant roots sum to the singlet contribution. -/
theorem root_sum
    {a b c₀ c₁ c₂ : ℝ}
    (hcos : c₀ + c₁ + c₂ = 0) :
    (a + 2 * b * c₀) + (a + 2 * b * c₁) + (a + 2 * b * c₂) =
      3 * a := by
  calc
    (a + 2 * b * c₀) + (a + 2 * b * c₁) + (a + 2 * b * c₂) =
        3 * a + 2 * b * (c₀ + c₁ + c₂) := by ring
    _ = 3 * a := by rw [hcos]; ring

/-- Their square sum separates into singlet and centered-plane powers. -/
theorem root_square_sum
    {a b c₀ c₁ c₂ : ℝ}
    (hcos : c₀ + c₁ + c₂ = 0)
    (hcosSq : c₀ ^ 2 + c₁ ^ 2 + c₂ ^ 2 = 3 / 2) :
    (a + 2 * b * c₀) ^ 2
        + (a + 2 * b * c₁) ^ 2
        + (a + 2 * b * c₂) ^ 2 =
      3 * a ^ 2 + 6 * b ^ 2 := by
  calc
    (a + 2 * b * c₀) ^ 2
          + (a + 2 * b * c₁) ^ 2
          + (a + 2 * b * c₂) ^ 2 =
        3 * a ^ 2
          + 4 * a * b * (c₀ + c₁ + c₂)
          + 4 * b ^ 2 * (c₀ ^ 2 + c₁ ^ 2 + c₂ ^ 2) := by ring
    _ = 3 * a ^ 2 + 6 * b ^ 2 := by rw [hcos, hcosSq]; ring

/-- The exact signed Koide quotient in terms of the modulus ratio. -/
theorem koide_formula
    {a b : ℝ}
    (ha : a ≠ 0) :
    (3 * a ^ 2 + 6 * b ^ 2) / (3 * a) ^ 2 =
      1 / 3 + (2 / 3) * (b / a) ^ 2 := by
  field_simp
  ring

/-- The conventional balance coordinate has the reciprocal-square-root form. -/
theorem balance_eq_inv_sqrt_two :
    Real.sqrt 2 / 2 = 1 / Real.sqrt 2 := by
  have hsqrt : Real.sqrt 2 ≠ 0 := by positivity
  field_simp
  nlinarith [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)]

/-- On the nonnegative modulus branch, Koide balance is exactly one condition
on the singlet-to-centered-plane norm ratio. -/
theorem koide_eq_two_thirds_iff
    {ratio : ℝ}
    (hratio : 0 ≤ ratio) :
    1 / 3 + (2 / 3) * ratio ^ 2 = 2 / 3 ↔
      ratio = 1 / Real.sqrt 2 := by
  have hsqrt_sq : (Real.sqrt 2) ^ 2 = 2 :=
    Real.sq_sqrt (by norm_num)
  have hsqrt_pos : 0 < Real.sqrt 2 := Real.sqrt_pos.2 (by norm_num)
  have hbalance_nonneg : 0 ≤ Real.sqrt 2 / 2 := by positivity
  have hbalance_sq : (Real.sqrt 2 / 2) ^ 2 = (1 : ℝ) / 2 := by
    nlinarith
  rw [← balance_eq_inv_sqrt_two]
  constructor
  · intro h
    have hratio_sq : ratio ^ 2 = (1 : ℝ) / 2 := by
      nlinarith
    nlinarith [sq_nonneg (ratio - Real.sqrt 2 / 2)]
  · intro h
    rw [h]
    nlinarith

/-! ## Phase cancellation: the cosine hypotheses hold for every δ

`root_sum` and `root_square_sum` take the two cosine identities as HYPOTHESES on abstract
coordinates `c₀ c₁ c₂`. Theorem 9.1 of the paper derives them for the actual cosines at angles
separated by `2π/3`, which is what makes "the phase δ cancels from Q" a theorem rather than an
assumption. The three results below discharge those hypotheses for EVERY `δ`, so the Koide
quotient is phase-independent unconditionally.

Scope unchanged: this adds no physical-mass claim, does not select the phase, and does not
supply the source-to-charged-family attachment. It closes exactly the trigonometric step. -/

/-- The three cosine coordinates at `2π/3` spacing sum to zero, for every phase `δ`. -/
theorem cos_coords_sum_zero (δ : ℝ) :
    Real.cos δ + Real.cos (δ + 2 * Real.pi / 3) + Real.cos (δ + 4 * Real.pi / 3) = 0 := by
  have h2 : Real.cos (δ + 2 * Real.pi / 3)
      = Real.cos δ * Real.cos (2 * Real.pi / 3) - Real.sin δ * Real.sin (2 * Real.pi / 3) :=
    Real.cos_add _ _
  have h4 : Real.cos (δ + 4 * Real.pi / 3)
      = Real.cos δ * Real.cos (4 * Real.pi / 3) - Real.sin δ * Real.sin (4 * Real.pi / 3) :=
    Real.cos_add _ _
  have c2 : Real.cos (2 * Real.pi / 3) = -(1 / 2) := by
    have h : (2 : ℝ) * Real.pi / 3 = Real.pi - Real.pi / 3 := by ring
    rw [h, Real.cos_pi_sub, Real.cos_pi_div_three]
  have c4 : Real.cos (4 * Real.pi / 3) = -(1 / 2) := by
    have h : (4 : ℝ) * Real.pi / 3 = Real.pi / 3 + Real.pi := by ring
    rw [h, Real.cos_add_pi, Real.cos_pi_div_three]
  have s2 : Real.sin (2 * Real.pi / 3) = Real.sqrt 3 / 2 := by
    have h : (2 : ℝ) * Real.pi / 3 = Real.pi - Real.pi / 3 := by ring
    rw [h, Real.sin_pi_sub, Real.sin_pi_div_three]
  have s4 : Real.sin (4 * Real.pi / 3) = -(Real.sqrt 3 / 2) := by
    have h : (4 : ℝ) * Real.pi / 3 = Real.pi / 3 + Real.pi := by ring
    rw [h, Real.sin_add_pi, Real.sin_pi_div_three]
  rw [h2, h4, c2, c4, s2, s4]; ring

/-- Their squares sum to `3/2`, for every phase `δ`. -/
theorem cos_coords_sq_sum (δ : ℝ) :
    Real.cos δ ^ 2 + Real.cos (δ + 2 * Real.pi / 3) ^ 2
      + Real.cos (δ + 4 * Real.pi / 3) ^ 2 = 3 / 2 := by
  have h2 : Real.cos (δ + 2 * Real.pi / 3)
      = Real.cos δ * Real.cos (2 * Real.pi / 3) - Real.sin δ * Real.sin (2 * Real.pi / 3) :=
    Real.cos_add _ _
  have h4 : Real.cos (δ + 4 * Real.pi / 3)
      = Real.cos δ * Real.cos (4 * Real.pi / 3) - Real.sin δ * Real.sin (4 * Real.pi / 3) :=
    Real.cos_add _ _
  have c2 : Real.cos (2 * Real.pi / 3) = -(1 / 2) := by
    have h : (2 : ℝ) * Real.pi / 3 = Real.pi - Real.pi / 3 := by ring
    rw [h, Real.cos_pi_sub, Real.cos_pi_div_three]
  have c4 : Real.cos (4 * Real.pi / 3) = -(1 / 2) := by
    have h : (4 : ℝ) * Real.pi / 3 = Real.pi / 3 + Real.pi := by ring
    rw [h, Real.cos_add_pi, Real.cos_pi_div_three]
  have s2 : Real.sin (2 * Real.pi / 3) = Real.sqrt 3 / 2 := by
    have h : (2 : ℝ) * Real.pi / 3 = Real.pi - Real.pi / 3 := by ring
    rw [h, Real.sin_pi_sub, Real.sin_pi_div_three]
  have s4 : Real.sin (4 * Real.pi / 3) = -(Real.sqrt 3 / 2) := by
    have h : (4 : ℝ) * Real.pi / 3 = Real.pi / 3 + Real.pi := by ring
    rw [h, Real.sin_add_pi, Real.sin_pi_div_three]
  have hsq3 : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  have pyth : Real.sin δ ^ 2 + Real.cos δ ^ 2 = 1 := Real.sin_sq_add_cos_sq δ
  rw [h2, h4, c2, c4, s2, s4]
  nlinarith [hsq3, pyth]

/-- **The phase cancels.** For every `δ`, the signed Koide quotient of the circulant roots built
from the actual cosine coordinates is `1/3 + (2/3)(b/a)^2` -- an expression in which `δ` does not
appear. This is the unconditional form of `koide_formula`: no cosine hypothesis is assumed. -/
theorem koide_quotient_phase_independent {a b : ℝ} (ha : a ≠ 0) (δ : ℝ) :
    ((a + 2 * b * Real.cos δ) ^ 2
        + (a + 2 * b * Real.cos (δ + 2 * Real.pi / 3)) ^ 2
        + (a + 2 * b * Real.cos (δ + 4 * Real.pi / 3)) ^ 2)
      / ((a + 2 * b * Real.cos δ)
          + (a + 2 * b * Real.cos (δ + 2 * Real.pi / 3))
          + (a + 2 * b * Real.cos (δ + 4 * Real.pi / 3))) ^ 2
      = 1 / 3 + (2 / 3) * (b / a) ^ 2 := by
  rw [root_sum (cos_coords_sum_zero δ),
      root_square_sum (cos_coords_sum_zero δ) (cos_coords_sq_sum δ)]
  exact koide_formula ha

end OPH.KoideCirculant
