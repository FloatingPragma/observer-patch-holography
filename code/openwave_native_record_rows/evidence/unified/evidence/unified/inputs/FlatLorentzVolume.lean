import Geometry.FlatDiamondVolume
import Mathlib.LinearAlgebra.Matrix.ToLinearEquiv
import Mathlib.MeasureTheory.Measure.Lebesgue.EqHaar
import Mathlib.MeasureTheory.Group.Prod

set_option autoImplicit false

namespace OPH.FlatLorentzVolume

open MeasureTheory Set Real
open OPH.SourceCausalBoundary OPH.FlatDiamondVolume
open scoped Matrix BigOperators

noncomputable section

def dot (x y : Space) : ℝ := x 0 * y 0 + x 1 * y 1 + x 2 * y 2
def square (z : Spacetime) : ℝ := z.1 ^ 2 - dot z.2 z.2

def boost (T τ : ℝ) (x : Space) (z : Spacetime) : Spacetime :=
  ((T * z.1 - dot x z.2) / τ,
    fun i => z.2 i + ((dot x z.2 / (T + τ) - z.1) / τ) * x i)

def boostMatrix (T τ : ℝ) (x : Space) : Matrix (Fin 4) (Fin 4) ℝ :=
  !![T / τ, -x 0 / τ, -x 1 / τ, -x 2 / τ;
     -x 0 / τ, 1 + x 0 * x 0 / (τ * (T + τ)), x 0 * x 1 / (τ * (T + τ)), x 0 * x 2 / (τ * (T + τ));
     -x 1 / τ, x 1 * x 0 / (τ * (T + τ)), 1 + x 1 * x 1 / (τ * (T + τ)), x 1 * x 2 / (τ * (T + τ));
     -x 2 / τ, x 2 * x 0 / (τ * (T + τ)), x 2 * x 1 / (τ * (T + τ)), 1 + x 2 * x 2 / (τ * (T + τ))]

theorem dot_self (x : Space) : dot x x = spatialNorm x ^ 2 := by
  unfold spatialNorm
  rw [EuclideanSpace.real_norm_sq_eq]
  simp [dot, Fin.sum_univ_succ]
  ring

theorem boost_square_identity (T τ : ℝ) (x : Space) (z : Spacetime)
    (hτ : τ ≠ 0) (hTτ : T + τ ≠ 0) :
    square (boost T τ x z) - square z =
      (T ^ 2 - τ ^ 2 - dot x x) * (z.1 - dot x z.2 / (T + τ)) ^ 2 / τ ^ 2 := by
  unfold square boost dot
  dsimp only
  field_simp
  ring

theorem boost_preserves_square {T τ : ℝ} (x : Space) (hτ : τ ≠ 0)
    (hTτ : T + τ ≠ 0) (hrel : T ^ 2 - dot x x = τ ^ 2) (z : Spacetime) :
    square (boost T τ x z) = square z := by
  have h := boost_square_identity T τ x z hτ hTτ
  have hz : T ^ 2 - τ ^ 2 - dot x x = 0 := by linarith
  rw [hz, zero_mul, zero_div] at h
  exact sub_eq_zero.mp h

set_option maxHeartbeats 800000 in
theorem boostMatrix_det (T τ : ℝ) (x : Space) (hτ : τ ≠ 0) (hTτ : T + τ ≠ 0) :
    (boostMatrix T τ x).det = (T * (T + τ) - dot x x) / (τ * (T + τ)) := by
  have h₁ : (1 : Fin 4).succAbove (2 : Fin 3) = 3 := by decide
  have h₂ : (2 : Fin 4).succAbove (2 : Fin 3) = 3 := by decide
  have h₃ : (3 : Fin 4).succAbove (2 : Fin 3) = 2 := by decide
  simp [boostMatrix, Matrix.det_succ_row_zero, Fin.sum_univ_succ, Matrix.submatrix, dot, h₁, h₂, h₃]
  field_simp
  ring

theorem boostMatrix_det_one {T τ : ℝ} (x : Space) (hτ : τ ≠ 0)
    (hTτ : T + τ ≠ 0) (hrel : T ^ 2 - dot x x = τ ^ 2) :
    (boostMatrix T τ x).det = 1 := by
  rw [boostMatrix_det T τ x hτ hTτ]
  apply (div_eq_one_iff_eq (mul_ne_zero hτ hTτ)).mpr
  nlinarith

theorem dot_le (x y : Space) : dot x y ≤ spatialNorm x * spatialNorm y := by
  have h := real_inner_le_norm (WithLp.toLp 2 x : EuclideanSpace ℝ (Fin 3))
    (WithLp.toLp 2 y)
  have hscalar (a b : ℝ) : inner ℝ a b = a * b := by
    change b * a = a * b
    ring
  simpa [PiLp.inner_apply, hscalar, Fin.sum_univ_succ, dot, spatialNorm,
    mul_comm, add_assoc] using h

theorem spatialNorm_nonneg (x : Space) : 0 ≤ spatialNorm x := norm_nonneg _

theorem boost_time_pos {T τ : ℝ} (x : Space) (hτ : 0 < τ) (hT : spatialNorm x < T)
    (z : Spacetime) (hz : spatialNorm z.2 ≤ z.1) (ht : 0 < z.1) :
    0 < (boost T τ x z).1 := by
  apply div_pos _ hτ
  have hxy := dot_le x z.2
  have hmul := mul_le_mul_of_nonneg_left hz (spatialNorm_nonneg x)
  nlinarith [mul_pos (sub_pos.mpr hT) ht]

theorem boost_zero (T τ : ℝ) (x : Space) : boost T τ x 0 = 0 := by
  simp [boost, dot]
  rfl

theorem boost_neg (T τ : ℝ) (x : Space) (z : Spacetime) :
    boost T τ x (-z) = -boost T τ x z := by
  apply Prod.ext
  · simp [boost, dot]
    ring
  · ext i
    simp [boost, dot]
    ring

/-- The explicit boost preserves future orientation as well as the quadratic form. -/
theorem boost_future_iff {T τ : ℝ} (x : Space) (hτ : 0 < τ) (hT : spatialNorm x < T)
    (hrel : T ^ 2 - dot x x = τ ^ 2) (z : Spacetime) :
    spatialNorm (boost T τ x z).2 ≤ (boost T τ x z).1 ↔ spatialNorm z.2 ≤ z.1 := by
  have hTp : 0 < T := (spatialNorm_nonneg x).trans_lt hT
  have hsq := boost_preserves_square x hτ.ne' (by linarith) hrel z
  simp only [square, dot_self] at hsq
  constructor
  · intro hb
    have hb0 : 0 ≤ (boost T τ x z).1 := (spatialNorm_nonneg _).trans hb
    have hq := sq_le_sq₀ (spatialNorm_nonneg _) hb0 |>.mpr hb
    have ht : 0 ≤ z.1 := by
      by_contra hnot
      have htneg : z.1 < 0 := lt_of_not_ge hnot
      have hzneg : spatialNorm (-z).2 ≤ (-z).1 := by
        have hn : spatialNorm (-z).2 = spatialNorm z.2 := by simp [spatialNorm]
        rw [hn]
        change spatialNorm z.2 ≤ -z.1
        nlinarith [spatialNorm_nonneg z.2]
      have hh := boost_time_pos x hτ hT (-z) hzneg (by simpa using neg_pos.mpr htneg)
      rw [boost_neg] at hh
      change 0 < -(boost T τ x z).1 at hh
      linarith
    nlinarith [spatialNorm_nonneg z.2]
  · intro hz
    have hz0 : 0 ≤ z.1 := (spatialNorm_nonneg _).trans hz
    have ht : 0 ≤ (boost T τ x z).1 := by
      change 0 ≤ (T * z.1 - dot x z.2) / τ
      apply div_nonneg _ hτ.le
      have hxy := dot_le x z.2
      have hmul := mul_le_mul_of_nonneg_left hz (spatialNorm_nonneg x)
      nlinarith [mul_nonneg (sub_nonneg.mpr hT.le) hz0]
    have hq := sq_le_sq₀ (spatialNorm_nonneg _) hz0 |>.mpr hz
    nlinarith [spatialNorm_nonneg (boost T τ x z).2]

theorem boost_tip {T τ : ℝ} (x : Space) (hτ : τ ≠ 0) (hTτ : T + τ ≠ 0)
    (hrel : T ^ 2 - dot x x = τ ^ 2) : boost T τ x (T, x) = (τ, 0) := by
  apply Prod.ext
  · change (T * T - dot x x) / τ = τ
    apply (div_eq_iff hτ).mpr
    nlinarith
  · ext i
    change x i + ((dot x x / (T + τ) - T) / τ) * x i = 0
    have hh : (dot x x / (T + τ) - T) / τ = -1 := by
      field_simp
      nlinarith
    rw [hh]
    ring

def split : (Fin 4 → ℝ) ≃ᵐ Spacetime :=
  MeasurableEquiv.piFinSuccAbove (fun _ : Fin 4 => ℝ) 0

theorem split_formula (v : Fin 4 → ℝ) : split v = (v 0, fun i => v i.succ) := by
  apply Prod.ext
  · rfl
  · ext i
    change v ((0 : Fin 4).succAbove i) = v i.succ
    simp

theorem matrix_split (T τ : ℝ) (x : Space) (v : Fin 4 → ℝ) :
    split ((boostMatrix T τ x).mulVec v) = boost T τ x (split v) := by
  simp only [split_formula]
  apply Prod.ext
  · simp [Matrix.mulVec, boostMatrix, boost, dot, Matrix.vecHead, Matrix.vecTail]
    ring
  · ext i
    fin_cases i <;> simp [Matrix.mulVec, boostMatrix, boost, dot, Matrix.vecHead, Matrix.vecTail,
      div_eq_mul_inv, mul_inv_rev] <;> ring

theorem boost_measurePreserving {T τ : ℝ} (x : Space) (hτ : τ ≠ 0)
    (hTτ : T + τ ≠ 0) (hrel : T ^ 2 - dot x x = τ ^ 2) :
    MeasurePreserving (boost T τ x) coordinateVolume coordinateVolume := by
  let A := Matrix.toLin' (boostMatrix T τ x)
  have hdet : LinearMap.det A = 1 := by
    rw [LinearMap.det_toLin']
    exact boostMatrix_det_one x hτ hTτ hrel
  have hA : MeasurePreserving A := by
    refine ⟨A.continuous_of_finiteDimensional.measurable, ?_⟩
    have hh := Real.map_linearMap_volume_pi_eq_smul_volume_pi (show LinearMap.det A ≠ 0 by simp [hdet])
    simpa only [hdet, inv_one, abs_one, ENNReal.ofReal_one, one_smul] using hh
  have hs : MeasurePreserving split := volume_preserving_piFinSuccAbove (fun _ : Fin 4 => ℝ) 0
  have heq : boost T τ x = split ∘ A ∘ split.symm := by
    funext z
    obtain ⟨v, rfl⟩ := split.surjective z
    simp only [Function.comp_apply, MeasurableEquiv.symm_apply_apply]
    exact (matrix_split T τ x v).symm
  rw [heq]
  exact hs.comp (hA.comp hs.symm)

theorem boost_sub (T τ : ℝ) (x : Space) (z w : Spacetime) :
    boost T τ x (z - w) = boost T τ x z - boost T τ x w := by
  apply Prod.ext
  · simp [boost, dot]
    ring
  · ext i
    simp [boost, dot]
    ring

/-- The preimage identity is about the actual causal diamond as a set. -/
theorem diamond_boost_preimage {T τ : ℝ} (x : Space) (hτ : 0 < τ)
    (hT : spatialNorm x < T) (hrel : T ^ 2 - dot x x = τ ^ 2) :
    diamond 1 0 (T, x) = boost T τ x ⁻¹' diamond 1 0 (τ, 0) := by
  have hTτ : T + τ ≠ 0 := by have := spatialNorm_nonneg x; linarith
  ext z
  have h₁ := boost_future_iff x hτ hT hrel z
  have h₂ := boost_future_iff x hτ hT hrel ((T, x) - z)
  rw [boost_sub, boost_tip x hτ.ne' hTτ hrel] at h₂
  simpa only [diamond, future, past, mem_preimage, mem_inter_iff, mem_setOf_eq, margin,
    one_mul, Prod.fst_zero, Prod.snd_zero, sub_zero, Prod.fst_sub, Prod.snd_sub,
    sub_nonneg] using (h₁.and h₂).symm

/-- Actual coordinate volume of a timelike diamond with arbitrary spatial tip. -/
theorem tilted_diamond_volume {T : ℝ} (x : Space) (hT : spatialNorm x < T) :
    coordinateVolume.real (diamond 1 0 (T, x)) =
      Real.pi * (T ^ 2 - spatialNorm x ^ 2) ^ 2 / 24 := by
  have hTp : 0 < T := (spatialNorm_nonneg x).trans_lt hT
  have hq : 0 < T ^ 2 - spatialNorm x ^ 2 := by
    nlinarith [spatialNorm_nonneg x]
  let τ := Real.sqrt (T ^ 2 - spatialNorm x ^ 2)
  have hτ : 0 < τ := Real.sqrt_pos.mpr hq
  have hrel : T ^ 2 - dot x x = τ ^ 2 := by rw [dot_self]; exact (Real.sq_sqrt hq.le).symm
  have hTτ : T + τ ≠ 0 := by linarith
  rw [measureReal_def, diamond_boost_preimage x hτ hT hrel,
    (boost_measurePreserving x hτ.ne' hTτ hrel).measure_preimage
      (diamond_measurable 1 0 (τ, 0)).nullMeasurableSet]
  change coordinateVolume.real (diamond 1 (0, 0) (τ, 0)) = _
  rw [full_vertical_diamond_volume hτ.le (by norm_num : (0 : ℝ) < 1)]
  have hs : τ ^ 2 = T ^ 2 - spatialNorm x ^ 2 := Real.sq_sqrt hq.le
  rw [show τ ^ 4 = (τ ^ 2) ^ 2 by ring, hs]
  ring

theorem tilted_diamond_volume_of_causal {T : ℝ} (x : Space) (hT : spatialNorm x ≤ T) :
    coordinateVolume.real (diamond 1 0 (T, x)) =
      Real.pi * (T ^ 2 - spatialNorm x ^ 2) ^ 2 / 24 := by
  rcases lt_or_eq_of_le hT with hlt | heq
  · exact tilted_diamond_volume x hlt
  · have hsub : diamond 1 0 (T, x) ⊆ {z | margin 1 0 z = 0} := by
      intro z hz
      have h₁ : 0 ≤ z.1 - spatialNorm z.2 := by
        simpa only [future, mem_setOf_eq, margin, one_mul, Prod.fst_zero, Prod.snd_zero, sub_zero] using hz.1
      have h₂ : 0 ≤ T - z.1 - spatialNorm (x - z.2) := by
        have hh : 0 ≤ 1 * (T - z.1) - spatialNorm (x - z.2) := hz.2
        simpa only [one_mul] using hh
      have htri : spatialNorm x ≤ spatialNorm (x - z.2) + spatialNorm z.2 := by
        simpa only [sub_zero] using spatialNorm_sub_triangle 0 z.2 x
      change 1 * (z.1 - 0) - spatialNorm (z.2 - 0) = 0
      simp only [one_mul, sub_zero]
      linarith
    have hnull := measure_mono_null hsub (coordinate_future_null (by norm_num : (1 : ℝ) ≠ 0) 0)
    simp [measureReal_def, hnull, heq]

theorem margin_translate (c : ℝ) (p q r : Spacetime) :
    margin c (p - r) (q - r) = margin c p q := by
  simp [margin, sub_sub_sub_cancel_right]

theorem diamond_translate_preimage (c : ℝ) (p q : Spacetime) :
    diamond c p q = (fun z => z - p) ⁻¹' diamond c 0 (q - p) := by
  ext z
  have h₁ := margin_translate c p z p
  have h₂ := margin_translate c z q p
  simp only [sub_self] at h₁
  simp only [diamond, future, past, mem_inter_iff, mem_preimage, mem_setOf_eq, h₁, h₂]

theorem coordinate_sub_preserving (p : Spacetime) :
    MeasurePreserving (fun z => z - p) coordinateVolume coordinateVolume := by
  simpa only [coordinateVolume, Prod.map_def, sub_eq_add_neg, Prod.fst_add,
    Prod.snd_add, Prod.fst_neg, Prod.snd_neg] using
    (measurePreserving_add_right (volume : Measure ℝ) (-p.1)).prod
      (measurePreserving_add_right (volume : Measure Space) (-p.2))

/-- Translated, possibly null, actual causal diamonds in coordinate volume. -/
theorem causal_diamond_volume (p q : Spacetime)
    (hpq : spatialNorm (q.2 - p.2) ≤ q.1 - p.1) :
    coordinateVolume.real (diamond 1 p q) =
      Real.pi * ((q.1 - p.1) ^ 2 - spatialNorm (q.2 - p.2) ^ 2) ^ 2 / 24 := by
  rw [measureReal_def, diamond_translate_preimage,
    (coordinate_sub_preserving p).measure_preimage (diamond_measurable 1 0 (q - p)).nullMeasurableSet]
  exact tilted_diamond_volume_of_causal (q.2 - p.2) hpq

#print axioms causal_diamond_volume
#print axioms tilted_diamond_volume
#print axioms boost_measurePreserving
#print axioms boost_future_iff
#print axioms boost_preserves_square
#print axioms boostMatrix_det_one

end
end OPH.FlatLorentzVolume
