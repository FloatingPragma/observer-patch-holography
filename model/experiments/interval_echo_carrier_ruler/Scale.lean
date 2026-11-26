import DispersionArmingInterface
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Inverse

namespace OPH.CarrierRuler

open OPH.PrimitivePortScaleBoundary OPH.DispersionArming

/-- Conditional antisymmetric probability at either receiver. -/
noncomputable def minus (x : ℝ) : ℝ := (1 - Real.cos x) / 2
noncomputable def recoveredPhase (p : ℝ) : ℝ := Real.arccos (1 - 2 * p)

theorem undo_phase {x : ℝ} (lo : 0 ≤ x) (hi : x ≤ Real.pi) :
    recoveredPhase (minus x) = x := by
  unfold recoveredPhase minus
  rw [show 1 - 2 * ((1 - Real.cos x) / 2) = Real.cos x by ring,
    Real.arccos_cos lo hi]

/-- Taking the squared modulus of the normalized antisymmetric output on two
equal-amplitude inputs gives the observed conditional minus probability. -/
theorem signed_exit (x y : ℝ) :
    ((Real.cos x - Real.cos y) ^ 2 + (Real.sin x - Real.sin y) ^ 2) / 4 =
      minus (x - y) := by
  rw [minus, Real.cos_sub]
  nlinarith [Real.sin_sq_add_cos_sq x, Real.sin_sq_add_cos_sq y]

/-- A reference reading cancels an unknown common positive phase gain on the enrolled branch. -/
theorem ratio_from_receipts {g a L : ℝ} (hg : 0 < g) (ha : 0 < a) (hL : 0 < L)
    (sample : g * a ≤ Real.pi) (reference : g * L ≤ Real.pi) :
    recoveredPhase (minus (g * a)) / recoveredPhase (minus (g * L)) = a / L := by
  rw [undo_phase (mul_pos hg ha).le sample, undo_phase (mul_pos hg hL).le reference]
  field_simp

/-- A second phase can produce the identical first intensity outside that branch. -/
theorem folded_phase (x : ℝ) : minus (2 * Real.pi - x) = minus x := by
  simp [minus, Real.cos_sub]

/-- Simultaneously scaling both lengths and inversely scaling the gradient leaves both readings fixed. -/
theorem unit_freedom (g a L c : ℝ) (hc : 0 < c) :
    minus ((g / c) * (c * a)) = minus (g * a) ∧
    minus ((g / c) * (c * L)) = minus (g * L) := by
  constructor <;> congr 1 <;> field_simp

/-- The declared standard fills the actual thirty-seam coefficient slot. -/
theorem declared_row (L : ℝ) (hL : 0 < L) :
    0 < 3 * L / 2 ∧ scaleC4 (3 * L / 2) = -(9 * L ^ 2) / 80 ∧
    scaleB0 (3 * L / 2) = 27 * L ^ 4 / 4480 ∧
    edgeScaleB6 (3 * L / 2) = -(9 * L ^ 4) / 22400 := by
  refine ⟨by positivity, ?_, ?_, ?_⟩ <;> simp only [scaleC4, scaleB0, edgeScaleB6] <;> ring

/-- The dilation changes the physical quartic slot while preserving the dimensionless ray. -/
theorem dilation (L : ℝ) (hL : 0 < L) :
    scaleC4 (3 * L / 2) ≠ scaleC4 (2 * L) ∧
    scaleB0 (2 * L) / scaleC4 (2 * L) ^ 2 = 10 / 21 ∧
    edgeScaleB6 (2 * L) / scaleC4 (2 * L) ^ 2 = -2 / 63 ∧
    edgeScaleB6 (2 * L) / scaleB0 (2 * L) = -1 / 15 := by
  refine ⟨?_, edgeScale_ray _ (by positivity)⟩
  intro heq
  have same := scale_parameter_faithful _ _ (by positivity) (by positivity) heq
  linarith

end OPH.CarrierRuler
