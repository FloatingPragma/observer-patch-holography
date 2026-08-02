import SeamCurrentEdge30Moment

open scoped BigOperators goldenRatio

namespace OPH.SeamCurrentEdge30Remainder

open OPH.SeamCurrentCarrierQuotient
open OPH.SeamCurrentEdge30Moment

/-!
# Eighth moment for the seam-current full-symbol remainder

The frozen edge-ray source uses the second, fourth, and sixth moments.  This
separate extension computes the eighth moment required to bound the first
omitted Taylor term without changing any frozen source artifact.

The result concerns the exact finite spatial symbol.  It supplies no physical
frequency, clock, field attachment, propagation law, or comparison result.
-/

abbrev Vec3 := OPH.PrimitivePortFrameQuotient.Vec3

private theorem sqrt5_pow_two : Real.sqrt 5 ^ 2 = (5 : ℝ) := by
  norm_num

private theorem sqrt5_pow_three : Real.sqrt 5 ^ 3 = 5 * Real.sqrt 5 := by
  calc
    Real.sqrt 5 ^ 3 = Real.sqrt 5 * Real.sqrt 5 ^ 2 := by ring
    _ = 5 * Real.sqrt 5 := by rw [sqrt5_pow_two]; ring

private theorem sqrt5_pow_four : Real.sqrt 5 ^ 4 = (25 : ℝ) := by
  calc
    Real.sqrt 5 ^ 4 = (Real.sqrt 5 ^ 2) ^ 2 := by ring
    _ = 25 := by rw [sqrt5_pow_two]; norm_num

private theorem sqrt5_pow_five : Real.sqrt 5 ^ 5 = 25 * Real.sqrt 5 := by
  calc
    Real.sqrt 5 ^ 5 = Real.sqrt 5 * Real.sqrt 5 ^ 4 := by ring
    _ = 25 * Real.sqrt 5 := by rw [sqrt5_pow_four]; ring

private theorem sqrt5_pow_six : Real.sqrt 5 ^ 6 = (125 : ℝ) := by
  calc
    Real.sqrt 5 ^ 6 = (Real.sqrt 5 ^ 2) ^ 3 := by ring
    _ = 125 := by rw [sqrt5_pow_two]; norm_num

private theorem sqrt5_pow_seven : Real.sqrt 5 ^ 7 = 125 * Real.sqrt 5 := by
  calc
    Real.sqrt 5 ^ 7 = Real.sqrt 5 * Real.sqrt 5 ^ 6 := by ring
    _ = 125 * Real.sqrt 5 := by rw [sqrt5_pow_six]; ring

private theorem sqrt5_pow_eight : Real.sqrt 5 ^ 8 = (625 : ℝ) := by
  calc
    Real.sqrt 5 ^ 8 = (Real.sqrt 5 ^ 4) ^ 2 := by ring
    _ = 625 := by rw [sqrt5_pow_four]; norm_num

/-- Normalized eighth moment of the complete intrinsic seam support. -/
noncomputable def seamMoment8 (k : Vec3) : ℝ :=
  ∑ e : Fin 30,
    OPH.PrimitivePortTranslationBridge.dot (carrierSeamDifference e) k ^ 8 / 256

set_option maxHeartbeats 12000000 in
/-- Exact eighth moment, split into its isotropic radial term and the radial
multiple of the normalized spin-six invariant. -/
theorem seamMoment8_eq (k : Vec3) :
    seamMoment8 k = (10 / 3 : ℝ) * radiusSquared k ^ 4
      - (8 / 15 : ℝ) * radiusSquared k * I6 k := by
  unfold seamMoment8
  simp_rw [carrierSeamDifference_eq_endpointDifference]
  simp [radiusSquared, I6, endpointDifference, portVector, seamLeft, seamRight,
    OPH.PrimitivePortTranslationBridge.dot, Fin.sum_univ_succ]
  ring_nf
  simp only [sqrt5_pow_eight, sqrt5_pow_seven, sqrt5_pow_six,
    sqrt5_pow_five, sqrt5_pow_four, sqrt5_pow_three, sqrt5_pow_two]
  ring

end OPH.SeamCurrentEdge30Remainder

/- Axiom audit: an exact finite table and a real polynomial identity only. -/

#print axioms OPH.SeamCurrentEdge30Remainder.seamMoment8_eq
