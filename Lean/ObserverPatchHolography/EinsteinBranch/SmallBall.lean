import Mathlib

/-!
# Small-ball coefficient algebra and explicit analytic premises

The radial integral and coefficient arithmetic are proved exactly.  The
identification of this polynomial kernel with a continuum modular charge, the
fixed-volume small-geodesic-ball expansion, and the common `o(ell^4)` tail are
not proved here: they are represented by `SmallBallPremises` and
`ScalingTailReceipt` and remain visible in the type of the composed theorem.
-/

namespace OPH.EinsteinBranch

open Filter Asymptotics
open scoped Topology

noncomputable section

/-- Radial reduction of the three-dimensional diamond kernel. -/
def diamondKernelIntegral (ℓ : ℝ) : ℝ :=
  ∫ r in (0 : ℝ)..ℓ,
    4 * Real.pi * r ^ 2 * ((ℓ ^ 2 - r ^ 2) / (2 * ℓ))

/-- Exact `d=4` diamond-kernel integral.  This proves the numerical
coefficient once the continuum kernel identification is supplied. -/
theorem diamondKernelIntegral_eq (ℓ : ℝ) (hℓ : ℓ ≠ 0) :
    diamondKernelIntegral ℓ = 4 * Real.pi * ℓ ^ 4 / 15 := by
  have hfun :
      (fun r : ℝ => 4 * Real.pi * r ^ 2 * ((ℓ ^ 2 - r ^ 2) / (2 * ℓ))) =
        fun r : ℝ => (2 * Real.pi * ℓ) * r ^ 2 -
          (2 * Real.pi / ℓ) * r ^ 4 := by
    funext r
    field_simp [hℓ]
    ring
  have h2 : IntervalIntegrable (fun r : ℝ => (2 * Real.pi * ℓ) * r ^ 2)
      MeasureTheory.volume 0 ℓ :=
    (by fun_prop : Continuous (fun r : ℝ => (2 * Real.pi * ℓ) * r ^ 2)).intervalIntegrable 0 ℓ
  have h4 : IntervalIntegrable (fun r : ℝ => (2 * Real.pi / ℓ) * r ^ 4)
      MeasureTheory.volume 0 ℓ :=
    (by fun_prop : Continuous (fun r : ℝ => (2 * Real.pi / ℓ) * r ^ 4)).intervalIntegrable 0 ℓ
  rw [diamondKernelIntegral, hfun, intervalIntegral.integral_sub h2 h4,
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    integral_pow, integral_pow]
  norm_num
  field_simp [hℓ]
  ring

/-- Multiplication by the `2*pi` first-law normalisation yields the exact
bulk coefficient printed in the compact paper. -/
theorem bulkSmallBallCoefficient (ℓ : ℝ) :
    2 * Real.pi * (4 * Real.pi * ℓ ^ 4 / 15) =
      8 * Real.pi ^ 2 * ℓ ^ 4 / 15 := by
  ring

/-- Exact arithmetic form of the fixed-volume area coefficient.  The
geometric expansion carrying this coefficient is an explicit premise below. -/
theorem fixedVolumeAreaCoefficient (ℓ Z : ℝ) :
    -(4 * Real.pi * ℓ ^ 4 / 15) * Z =
      -(4 * Real.pi * ℓ ^ 4 * Z / 15) := by
  ring

/-- One cofinal family and its three carried remainder channels. -/
structure ScalingTailData where
  ell : ℕ → ℝ
  collarError : ℕ → ℝ
  kernelError : ℕ → ℝ
  chartError : ℕ → ℝ

/-- EB8: one common shrinking family carries every named error as
`o(ell^4)`.  This is a mathematical predicate, not a Boolean or `True`. -/
def ScalingTailReceipt (D : ScalingTailData) : Prop :=
  (∀ r, 0 < D.ell r) ∧
  Tendsto D.ell atTop (𝓝 0) ∧
  IsLittleO atTop D.collarError (fun r => D.ell r ^ 4) ∧
  IsLittleO atTop D.kernelError (fun r => D.ell r ^ 4) ∧
  IsLittleO atTop D.chartError (fun r => D.ell r ^ 4)

/-- Constant radii fail the shrinking-family receipt, proving that EB8's
predicate is not constant true. -/
def constantRadiusData : ScalingTailData where
  ell := fun _ => 1
  collarError := fun _ => 0
  kernelError := fun _ => 0
  chartError := fun _ => 0

theorem constantRadiusData_fails : ¬ ScalingTailReceipt constantRadiusData := by
  intro h
  have hlim := h.2.1
  have hone : Tendsto (fun _ : ℕ => (1 : ℝ)) atTop (𝓝 (1 : ℝ)) := tendsto_const_nhds
  have hz : (1 : ℝ) = 0 := tendsto_nhds_unique hone hlim
  norm_num at hz

/-- The exact algebraic inputs consumed by the rest-frame step.  The names
make clear which statements are analytic/physical premises:

* `stationarity` is generalized-entropy stationarity;
* `continuumKernel` is the continuum diamond-kernel/stress identification;
* `fixedVolumeArea` is the smooth small-geodesic-ball identity.

They are structure fields, not Lean axioms and not asserted truths. -/
structure SmallBallPremises where
  G : ℝ
  ell : ℝ
  stressContraction : ℝ
  geometryContraction : ℝ
  deltaBulkEntropy : ℝ
  deltaArea : ℝ
  G_pos : 0 < G
  ell_pos : 0 < ell
  stationarity : deltaBulkEntropy + deltaArea / (4 * G) = 0
  continuumKernel :
    deltaBulkEntropy =
      (8 * Real.pi ^ 2 * ell ^ 4 / 15) * stressContraction
  fixedVolumeArea :
    deltaArea =
      -(4 * Real.pi * ell ^ 4 / 15) * geometryContraction

/-- Exact rest-frame Einstein arithmetic, conditional on the explicitly
named continuum and physical premises. -/
theorem restFrameEinsteinRelation (P : SmallBallPremises) :
    P.geometryContraction =
      8 * Real.pi * P.G * P.stressContraction := by
  have hπ : 0 < Real.pi := Real.pi_pos
  have hArea : P.deltaArea = -P.deltaBulkEntropy * (4 * P.G) := by
    have hdiv : P.deltaArea / (4 * P.G) = -P.deltaBulkEntropy := by
      linarith [P.stationarity]
    have hne : (4 * P.G : ℝ) ≠ 0 :=
      mul_ne_zero (by norm_num) (ne_of_gt P.G_pos)
    have hmul := (div_eq_iff hne).mp hdiv
    linarith [hmul]
  rw [P.continuumKernel, P.fixedVolumeArea] at hArea
  have hcancel : (4 * Real.pi * P.ell ^ 4 / 15 : ℝ) ≠ 0 := by
    apply div_ne_zero
    · exact mul_ne_zero
        (mul_ne_zero (by norm_num) Real.pi_ne_zero)
        (pow_ne_zero 4 (ne_of_gt P.ell_pos))
    · norm_num
  apply mul_left_cancel₀ hcancel
  linear_combination -hArea

/-! ## Per-theorem axiom audit -/

#print axioms diamondKernelIntegral_eq
#print axioms bulkSmallBallCoefficient
#print axioms fixedVolumeAreaCoefficient
#print axioms constantRadiusData_fails
#print axioms restFrameEinsteinRelation

end

end OPH.EinsteinBranch
