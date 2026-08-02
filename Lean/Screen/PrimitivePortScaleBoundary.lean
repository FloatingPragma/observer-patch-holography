import Mathlib
import A5PrimitivePortPrediction

namespace OPH.PrimitivePortScaleBoundary

/-!
# Primitive-port scale boundary

This module isolates the algebraic scale freedom in the frozen primitive-port
coefficient packet.  The normalized quadratic coefficient is a separate
constant, while the quartic and sixth-order coefficients carry powers of a
declared positive-length parameter `a`.  Exact ratios remove that parameter.

These theorems are a nonidentifiability and counterfamily boundary.  They do
not construct source models at every scale, select a physical value of `a`, or
derive the metric coefficient `kappa`.  In particular, the relation
`a^2 = kappa * P * ell^2` does not select `kappa = 1`.
-/

/-- The normalized coefficient of the leading `k^2` term.  It is kept
separate from the scale-dependent corrections. -/
def normalizedLeadingCoefficient : ℝ := 1

/-- Scale-dependent coefficient of the isotropic `k^4` correction. -/
noncomputable def scaleC4 (a : ℝ) : ℝ := -(a ^ 2) / 20

/-- Scale-dependent coefficient of the isotropic `k^6` correction. -/
noncomputable def scaleB0 (a : ℝ) : ℝ := a ^ 4 / 840

/-- Scale-dependent coefficient of the anisotropic `k^6` correction. -/
noncomputable def scaleB6 (a : ℝ) : ℝ := 2 * a ^ 4 / 7875

/-- The scale-dependent coefficients are the frozen dimensionless constants
multiplied by the appropriate power of `a`. -/
theorem coefficients_match_frozen_source (a : ℝ) :
    scaleC4 a = a ^ 2 * OPH.A5PrimitivePortPrediction.C4 ∧
      scaleB0 a = a ^ 4 * OPH.A5PrimitivePortPrediction.B0 ∧
      scaleB6 a = a ^ 4 * OPH.A5PrimitivePortPrediction.B6 := by
  constructor
  · norm_num [scaleC4, OPH.A5PrimitivePortPrediction.C4]
    ring
  constructor
  · norm_num [scaleB0, OPH.A5PrimitivePortPrediction.B0]
    ring
  · norm_num [scaleB6, OPH.A5PrimitivePortPrediction.B6]
    ring

/-- Rescaling `a` by `s` rescales the quartic coefficient by `s^2`. -/
theorem scaleC4_rescale (s a : ℝ) :
    scaleC4 (s * a) = s ^ 2 * scaleC4 a := by
  unfold scaleC4
  ring

/-- Rescaling `a` by `s` rescales the isotropic sixth-order coefficient by
`s^4`. -/
theorem scaleB0_rescale (s a : ℝ) :
    scaleB0 (s * a) = s ^ 4 * scaleB0 a := by
  unfold scaleB0
  ring

/-- Rescaling `a` by `s` rescales the anisotropic sixth-order coefficient by
`s^4`. -/
theorem scaleB6_rescale (s a : ℝ) :
    scaleB6 (s * a) = s ^ 4 * scaleB6 a := by
  unfold scaleB6
  ring

/-- The normalized leading coefficient is independent of every rescaling of
`a`; no scale variable enters its definition. -/
theorem normalizedLeadingCoefficient_rescale (_s _a : ℝ) :
    normalizedLeadingCoefficient = 1 := by
  rfl

/-- The isotropic sixth-order to squared-quartic ratio is scale free. -/
theorem scaleB0_over_scaleC4_sq (a : ℝ) (ha : a ≠ 0) :
    scaleB0 a / scaleC4 a ^ 2 = 10 / 21 := by
  unfold scaleB0 scaleC4
  field_simp [ha]
  ring

/-- The anisotropic sixth-order to squared-quartic ratio is scale free. -/
theorem scaleB6_over_scaleC4_sq (a : ℝ) (ha : a ≠ 0) :
    scaleB6 a / scaleC4 a ^ 2 = 32 / 315 := by
  unfold scaleB6 scaleC4
  field_simp [ha]
  ring

/-- The two sixth-order coefficients have a scale-free ratio. -/
theorem scaleB6_over_scaleB0 (a : ℝ) (ha : a ≠ 0) :
    scaleB6 a / scaleB0 a = 16 / 75 := by
  unfold scaleB6 scaleB0
  field_simp [ha]
  ring

/-- A conditional metric scale relation determines only the dimensionless
ratio `(a / ell)^2` when the reference length is nonzero. -/
theorem metric_ratio_of_scale_relation
    (a kappa P ell : ℝ) (hell : ell ≠ 0)
    (hscale : a ^ 2 = kappa * P * ell ^ 2) :
    (a / ell) ^ 2 = kappa * P := by
  calc
    (a / ell) ^ 2 = a ^ 2 / ell ^ 2 := by ring
    _ = (kappa * P * ell ^ 2) / ell ^ 2 := by rw [hscale]
    _ = kappa * P := by field_simp [hell]

/-- Under the additional sign conditions, the conditional metric relation
has the positive-square-root form. -/
theorem metric_ratio_eq_sqrt
    (a kappa P ell : ℝ) (ha : 0 ≤ a) (hell : 0 < ell)
    (hkp : 0 ≤ kappa * P)
    (hscale : a ^ 2 = kappa * P * ell ^ 2) :
    a / ell = Real.sqrt (kappa * P) := by
  have hell0 : ell ≠ 0 := ne_of_gt hell
  have hratio := metric_ratio_of_scale_relation a kappa P ell hell0 hscale
  have hleft : 0 ≤ a / ell := div_nonneg ha (le_of_lt hell)
  have hright : 0 ≤ Real.sqrt (kappa * P) := Real.sqrt_nonneg _
  have hsqrt : (Real.sqrt (kappa * P)) ^ 2 = kappa * P :=
    Real.sq_sqrt hkp
  nlinarith

/-- Pure algebra permits a continuum of metric coefficients at fixed nonzero
`P`: choosing `a = s * ell` is compatible with `kappa = s^2 / P`.  This is
not an existence theorem for the source dynamics. -/
theorem metric_coefficient_counterfamily
    (s P ell : ℝ) (hP : P ≠ 0) :
    (s * ell) ^ 2 = (s ^ 2 / P) * P * ell ^ 2 := by
  field_simp [hP]

#print axioms coefficients_match_frozen_source
#print axioms scaleC4_rescale
#print axioms scaleB0_rescale
#print axioms scaleB6_rescale
#print axioms normalizedLeadingCoefficient_rescale
#print axioms scaleB0_over_scaleC4_sq
#print axioms scaleB6_over_scaleC4_sq
#print axioms scaleB6_over_scaleB0
#print axioms metric_ratio_of_scale_relation
#print axioms metric_ratio_eq_sqrt
#print axioms metric_coefficient_counterfamily

end OPH.PrimitivePortScaleBoundary
