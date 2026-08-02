import Mathlib

namespace OPH.VolumeReadoutBridge

/-!
# Conditional screen-volume readout algebra

The determinant kernel of the screen-to-curvature identification is exact.
In three dimensions a conformal rescaling changes the positive square-root
determinant by `exp (3 * zeta)`.  If a separate physical attachment proves
that the quotient-visible collar readout and its frozen reference are those
positive spatial volume densities, its logarithmic scalar is exactly `zeta`.

This file isolates that implication from the physical premises needed to use
it.  It does not construct a uniform-density cut, identify
`lambda * sqrt (det sigma)` with the induced spatial volume density, derive a
conserved stress tensor, prove adiabaticity or freeze-out, or attach the field
to cosmological data.  The final counterexample records why positivity of a
collar readout alone cannot replace the missing volume-ratio premise.
-/

/-- The unprojected scalar carried by a positive readout and positive frozen
reference.  Positivity is a premise of the bridge theorems rather than part of
this algebraic definition. -/
noncomputable def screenScalar (J Jbar : ℝ) : ℝ :=
  (1 / 3 : ℝ) * Real.log (J / Jbar)

/-- Three-volume density of the declared conformal scalar branch.  The common
positive factors `a^3` and `rhoRef` encode the background scale and reference
density. -/
noncomputable def conformalVolumeDensity (a rhoRef zeta : ℝ) : ℝ :=
  a ^ 3 * Real.exp (3 * zeta) * rhoRef

/-- Frozen-background three-volume density with the same scale and reference
density. -/
noncomputable def backgroundVolumeDensity (a rhoRef : ℝ) : ℝ :=
  a ^ 3 * rhoRef

/-- A collar-adapted three-metric in spatial ADM form.  The first coordinate
is normal to the two-surface, `lambda` is its positive normal lapse,
`beta0,beta1` are tangential shifts, and `s00,s01,s11` are the components of
the induced two-metric. -/
def collarMetric
    (lambda s00 s01 s11 beta0 beta1 : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  !![
    lambda ^ 2 + s00 * beta0 ^ 2 + 2 * s01 * beta0 * beta1 + s11 * beta1 ^ 2,
      s00 * beta0 + s01 * beta1,
      s01 * beta0 + s11 * beta1;
    s00 * beta0 + s01 * beta1, s00, s01;
    s01 * beta0 + s11 * beta1, s01, s11]

/-- Exact collar determinant identity.  Tangential shift cancels, leaving the
square normal lapse times the determinant of the induced two-metric. -/
theorem collarMetric_determinant
    (lambda s00 s01 s11 beta0 beta1 : ℝ) :
    Matrix.det (collarMetric lambda s00 s01 s11 beta0 beta1) =
      lambda ^ 2 * (s00 * s11 - s01 ^ 2) := by
  simp [collarMetric, Matrix.det_fin_three]
  ring

/-- Under the standard two-dimensional positive-definiteness conditions, the
collar product `lambda * sqrt(det sigma)` equals the positive square root of
the collar determinant.  Calling this a physical spatial volume density also
requires the separate metric-attachment premise described above. -/
theorem collar_density_eq_sqrt_det
    (lambda s00 s01 s11 beta0 beta1 : ℝ)
    (hlambda : 0 < lambda) (_hS00 : 0 < s00)
    (_hSigma : 0 < s00 * s11 - s01 ^ 2) :
    lambda * Real.sqrt (s00 * s11 - s01 ^ 2) =
      Real.sqrt (Matrix.det (collarMetric lambda s00 s01 s11 beta0 beta1)) := by
  rw [collarMetric_determinant]
  rw [Real.sqrt_mul (sq_nonneg lambda)]
  rw [Real.sqrt_sq_eq_abs, abs_of_pos hlambda]

/-- In three dimensions, multiplying a three-by-three matrix by
`exp (2*zeta)` multiplies its determinant by `exp (6*zeta)`.  Applying the
identity to a spatial metric requires symmetry and positive definiteness. -/
theorem conformal_metric_determinant
    (hBar : Matrix (Fin 3) (Fin 3) ℝ) (zeta : ℝ) :
    Matrix.det (Real.exp (2 * zeta) • hBar) =
      Real.exp (6 * zeta) * Matrix.det hBar := by
  rw [Matrix.det_smul]
  simp only [Fintype.card_fin]
  congr 1
  rw [show Real.exp (2 * zeta) ^ 3 =
      Real.exp (2 * zeta) * Real.exp (2 * zeta) * Real.exp (2 * zeta) by ring]
  rw [← Real.exp_add, ← Real.exp_add]
  congr 1
  ring

/-- For a positive reference determinant, the positive square-root
determinant ratio is `exp (3*zeta)`.  Its spatial-volume interpretation uses
the separate metric premise. -/
theorem conformal_metric_sqrt_det_ratio
    (hBar : Matrix (Fin 3) (Fin 3) ℝ) (zeta : ℝ)
    (hDet : 0 < Matrix.det hBar) :
    Real.sqrt (Matrix.det (Real.exp (2 * zeta) • hBar)) /
        Real.sqrt (Matrix.det hBar) = Real.exp (3 * zeta) := by
  rw [conformal_metric_determinant]
  rw [show Real.exp (6 * zeta) = Real.exp (3 * zeta) ^ 2 by
    rw [pow_two, ← Real.exp_add]
    congr 1
    ring]
  rw [Real.sqrt_mul (sq_nonneg (Real.exp (3 * zeta)))]
  rw [Real.sqrt_sq_eq_abs, abs_of_pos (Real.exp_pos _)]
  field_simp [Real.sqrt_ne_zero'.mpr hDet]

/-- With one common positive background scale and reference density, the
three-volume-density ratio is exactly the conformal factor `exp (3*zeta)`. -/
theorem conformal_volume_ratio
    (a rhoRef zeta : ℝ) (ha : 0 < a) (hrho : 0 < rhoRef) :
    conformalVolumeDensity a rhoRef zeta /
        backgroundVolumeDensity a rhoRef = Real.exp (3 * zeta) := by
  apply (div_eq_iff ?_).2
  · simp only [conformalVolumeDensity, backgroundVolumeDensity]
    rw [show 3 * zeta = zeta * 3 by ring]
    ring
  · exact mul_ne_zero (pow_ne_zero 3 ha.ne') hrho.ne'

/-- The exact conditional kernel requested by the native volume-readout
bridge: once the collar ratio equals the spatial volume ratio, the logarithmic
screen scalar is `zeta`. -/
theorem screenScalar_eq_zeta_of_conformal_volume_readout
    (J Jbar a rhoRef zeta : ℝ)
    (ha : 0 < a) (hrho : 0 < rhoRef)
    (hJ : J = conformalVolumeDensity a rhoRef zeta)
    (hJbar : Jbar = backgroundVolumeDensity a rhoRef) :
    screenScalar J Jbar = zeta := by
  rw [hJ, hJbar, screenScalar,
    conformal_volume_ratio a rhoRef zeta ha hrho, Real.log_exp]
  ring

/-- Equivalent interface form.  A producer may certify the ratio directly,
without exposing a coordinate representation of the metric. -/
theorem screenScalar_eq_zeta_of_ratio
    (J Jbar zeta : ℝ)
    (_hJ : 0 < J) (_hJbar : 0 < Jbar)
    (hRatio : J / Jbar = Real.exp (3 * zeta)) :
    screenScalar J Jbar = zeta := by
  rw [screenScalar, hRatio, Real.log_exp]
  ring

/-- Direct determinant form of the conditional identity for a
positive-determinant reference matrix.  A physical metric interpretation
requires symmetry and positive definiteness. -/
theorem screenScalar_of_conformal_metric
    (hBar : Matrix (Fin 3) (Fin 3) ℝ) (zeta : ℝ)
    (hDet : 0 < Matrix.det hBar) :
    screenScalar
        (Real.sqrt (Matrix.det (Real.exp (2 * zeta) • hBar)))
        (Real.sqrt (Matrix.det hBar)) = zeta := by
  apply screenScalar_eq_zeta_of_ratio
  · exact Real.sqrt_pos.2 (by
      rw [conformal_metric_determinant]
      exact mul_pos (Real.exp_pos _) hDet)
  · exact Real.sqrt_pos.2 hDet
  · exact conformal_metric_sqrt_det_ratio hBar zeta hDet

/-- A positive calibration ratio contributes only its logarithmic scalar.
This makes the residual reference-density ambiguity explicit. -/
theorem screenScalar_of_calibrated_ratio
    (J Jbar calibration zeta : ℝ)
    (_hJ : 0 < J) (_hJbar : 0 < Jbar)
    (hcalibration : 0 < calibration)
    (hRatio : J / Jbar = calibration * Real.exp (3 * zeta)) :
    screenScalar J Jbar = zeta + (1 / 3 : ℝ) * Real.log calibration := by
  rw [screenScalar, hRatio,
    Real.log_mul hcalibration.ne' (Real.exp_ne_zero (3 * zeta)), Real.log_exp]
  ring

/-- Pointwise screen scalar on a finite or continuum index type. -/
noncomputable def screenScalarField {I : Type*} (J Jbar : I → ℝ) : I → ℝ :=
  fun x => screenScalar (J x) (Jbar x)

/-- The exact ratio premise implies equality of the complete scalar fields. -/
theorem screenScalarField_eq_zeta_of_ratio {I : Type*}
    (J Jbar zeta : I → ℝ)
    (hJ : ∀ x, 0 < J x) (hJbar : ∀ x, 0 < Jbar x)
    (hRatio : ∀ x, J x / Jbar x = Real.exp (3 * zeta x)) :
    screenScalarField J Jbar = zeta := by
  funext x
  exact screenScalar_eq_zeta_of_ratio
    (J x) (Jbar x) (zeta x) (hJ x) (hJbar x) (hRatio x)

/-- A common positive factor cancels between two positive-density readouts.
A native producer must still prove that a physical rechart transforms both
readouts by this same density factor. -/
theorem screenScalar_common_density_factor
    (J Jbar jacobian : ℝ) (_hJ : 0 < J) (hJbar : 0 < Jbar)
    (hJacobian : 0 < jacobian) :
    screenScalar (jacobian * J) (jacobian * Jbar) = screenScalar J Jbar := by
  unfold screenScalar
  congr 2
  field_simp [hJbar.ne', hJacobian.ne']

/-- Relabeling an index type by a bijection commutes definitionally with the
scalar readout.  This is an index identity, not a physical covariance theorem
or a proof that `J` factors through the physical quotient. -/
theorem screenScalarField_reindex {I I' : Type*} (e : I ≃ I')
    (J Jbar : I → ℝ) :
    screenScalarField (fun y => J (e.symm y)) (fun y => Jbar (e.symm y)) =
      fun y => screenScalarField J Jbar (e.symm y) := by
  rfl

/-- A low-mode calibration defect disappears after any linear projector that
annihilates it.  This covers the monopole/dipole removal step without assuming
that a particular numerical projector is the physical one. -/
theorem projected_scalar_eq_projected_zeta {I : Type*}
    (projectHigh : (I → ℝ) →ₗ[ℝ] (I → ℝ))
    (q zeta lowMode : I → ℝ)
    (hq : q = zeta + lowMode)
    (hLow : projectHigh lowMode = 0) :
    projectHigh q = projectHigh zeta := by
  rw [hq, map_add, hLow, add_zero]

/-- Combined calibrated-ratio and low-mode theorem.  A reference mismatch is
harmless precisely when its logarithm lies in the declared removed subspace. -/
theorem projected_calibrated_ratio_eq_projected_zeta {I : Type*}
    (projectHigh : (I → ℝ) →ₗ[ℝ] (I → ℝ))
    (J Jbar calibration zeta : I → ℝ)
    (hJ : ∀ x, 0 < J x) (hJbar : ∀ x, 0 < Jbar x)
    (hCalibration : ∀ x, 0 < calibration x)
    (hRatio : ∀ x, J x / Jbar x = calibration x * Real.exp (3 * zeta x))
    (hLow : projectHigh (fun x => (1 / 3 : ℝ) * Real.log (calibration x)) = 0) :
    projectHigh (screenScalarField J Jbar) = projectHigh zeta := by
  apply projected_scalar_eq_projected_zeta projectHigh
      (screenScalarField J Jbar) zeta
      (fun x => (1 / 3 : ℝ) * Real.log (calibration x))
  · funext x
    exact screenScalar_of_calibrated_ratio
      (J x) (Jbar x) (calibration x) (zeta x)
      (hJ x) (hJbar x) (hCalibration x) (hRatio x)
  · exact hLow

/-- Positivity and the bare collar-product form do not determine curvature.
The unit collar and unit background are positive, but their scalar is zero and
therefore differs from the nonzero candidate `zeta = 1`. -/
theorem positive_collar_readout_does_not_force_zeta :
    ∃ (lambda sigmaDet lambdaBar sigmaDetBar zeta : ℝ),
      0 < lambda ∧ 0 < sigmaDet ∧ 0 < lambdaBar ∧ 0 < sigmaDetBar ∧
      screenScalar (lambda * Real.sqrt sigmaDet)
          (lambdaBar * Real.sqrt sigmaDetBar) ≠ zeta := by
  refine ⟨1, 1, 1, 1, 1, by norm_num, by norm_num, by norm_num, by norm_num, ?_⟩
  norm_num [screenScalar]

#print axioms conformal_volume_ratio
#print axioms collarMetric_determinant
#print axioms collar_density_eq_sqrt_det
#print axioms conformal_metric_determinant
#print axioms conformal_metric_sqrt_det_ratio
#print axioms screenScalar_eq_zeta_of_conformal_volume_readout
#print axioms screenScalar_eq_zeta_of_ratio
#print axioms screenScalar_of_conformal_metric
#print axioms screenScalar_of_calibrated_ratio
#print axioms screenScalarField_eq_zeta_of_ratio
#print axioms screenScalar_common_density_factor
#print axioms screenScalarField_reindex
#print axioms projected_scalar_eq_projected_zeta
#print axioms projected_calibrated_ratio_eq_projected_zeta
#print axioms positive_collar_readout_does_not_force_zeta

end OPH.VolumeReadoutBridge
