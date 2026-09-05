import CarrierFrequencySpeed
import Mathlib.Analysis.Calculus.Deriv.MeanValue

open scoped BigOperators

namespace OPH.SeamMaxwellContinuum

open OPH.PrimitivePortTranslationBridge
open OPH.SeamCurrentEdge30Moment
open OPH.SeamCurrentHomogeneousAction
open OPH.SeamCurrentFreePhotonLift
open OPH.SeamCurrentPhotonLeptonThreshold

/-!
# Quantitative continuum control of the complete seam symbol

The exact thirty-seam symbol converges quadratically to the Euclidean
Laplacian symbol. These estimates control the cosine and sine entries of
the declared transverse Maxwell propagator. They hold at every momentum,
not just for a truncated polynomial or a sampled set of directions.

The accompanying paper constructs the real L2 field group and proves its
strong continuum limit by Plancherel and dominated convergence. Those
functional-analytic steps are paper proofs, not declarations in this file.
The common Euclidean domain, scale and reversible field evolution are
supplied; this does not construct spacetime or select physical time, fields,
current, units, or a detector from an authenticated source history.
-/

private theorem cubic_sine_lower {x : ℝ} (hx : 0 ≤ x) :
    x - x ^ 3 / 6 ≤ Real.sin x := by
  have hm : Monotone (fun y : ℝ ↦ Real.sin y - y + y ^ 3 / 6) := by
    apply monotone_of_hasDerivAt_nonneg
      (f' := fun y ↦ Real.cos y - 1 + 3 * y ^ 2 / 6)
    · intro y
      convert ((Real.hasDerivAt_sin y).sub (hasDerivAt_id y)).add
        (((hasDerivAt_id y).pow 3).div_const 6) using 1
      simp only [id_eq]
      ring
    · intro y
      change 0 ≤ Real.cos y - 1 + 3 * y ^ 2 / 6
      linarith [Real.one_sub_sq_div_two_le_cos (x := y)]
  have h : Real.sin 0 - 0 + 0 ^ 3 / 6 ≤ Real.sin x - x + x ^ 3 / 6 := hm hx
  norm_num at h
  linarith

/-- Global fourth-order one-sided cosine bound, proved by two derivative
comparisons rather than a small-argument expansion. -/
theorem cosine_quartic_upper (x : ℝ) :
    Real.cos x ≤ 1 - x ^ 2 / 2 + x ^ 4 / 24 := by
  suffices h : ∀ y : ℝ, 0 ≤ y →
      Real.cos y ≤ 1 - y ^ 2 / 2 + y ^ 4 / 24 by
    have h4 : |x| ^ 4 = x ^ 4 := by
      calc
        |x| ^ 4 = (|x| ^ 2) ^ 2 := by ring
        _ = (x ^ 2) ^ 2 := by rw [sq_abs]
        _ = x ^ 4 := by ring
    simpa only [Real.cos_abs, sq_abs, h4] using h |x| (abs_nonneg x)
  intro y hy
  have hm : MonotoneOn (fun z : ℝ ↦
      1 - z ^ 2 / 2 + z ^ 4 / 24 - Real.cos z) (Set.Ici 0) := by
    apply monotoneOn_of_hasDerivWithinAt_nonneg (convex_Ici 0)
      (f' := fun z ↦ -(2 * z / 2) + 4 * z ^ 3 / 24 + Real.sin z)
    · fun_prop
    · intro z _
      apply HasDerivAt.hasDerivWithinAt
      convert (((hasDerivAt_const z (1 : ℝ)).sub
        (((hasDerivAt_id z).pow 2).div_const 2)).add
        (((hasDerivAt_id z).pow 4).div_const 24)).sub
        (Real.hasDerivAt_cos z) using 1
      simp only [id_eq]
      ring
    · intro z hz
      have hz' : 0 ≤ z := (interior_subset hz : z ∈ Set.Ici 0)
      linarith [cubic_sine_lower hz']
  have h : 1 - (0 : ℝ) ^ 2 / 2 + 0 ^ 4 / 24 - Real.cos 0 ≤
      1 - y ^ 2 / 2 + y ^ 4 / 24 - Real.cos y := hm (by simp) hy hy
  norm_num at h
  linarith

/-- The source fourth moment in the unit-seam normalization. -/
theorem unit_seam_fourth_moment_eq (k : ThresholdVec3) :
    (∑ e : Fin 30, dot k (unitCarrierSeamDirection e) ^ 4) =
      6 * radiusSquared k ^ 2 := by
  rw [← seamMoment4_eq k]
  unfold seamMoment4
  apply Finset.sum_congr rfl
  intro e _
  unfold unitCarrierSeamDirection
  rw [dot_smul_right, OPH.SeamCurrentFreePhotonLift.dot_comm k]
  ring

/-- Global consistency bound for the full positive cosine operator. -/
theorem exact_symbol_quadratic_error {a : ℝ} (ha : a ≠ 0)
    (k : ThresholdVec3) :
    0 ≤ radiusSquared k -
        OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k ∧
    radiusSquared k -
        OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol a k ≤
      a ^ 2 * radiusSquared k ^ 2 / 20 := by
  refine ⟨sub_nonneg.mpr (exact_fz12_symbol_le_radiusSquared ha k), ?_⟩
  rw [OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol_eq_edgeCurrentCharacterSymbol
    a ha]
  have hsum :
      ∑ e : Fin 30, ((a * dot k (unitCarrierSeamDirection e)) ^ 2 / 2 -
        (a * dot k (unitCarrierSeamDirection e)) ^ 4 / 24) ≤
      ∑ e : Fin 30, (1 - Real.cos (a * dot k (unitCarrierSeamDirection e))) := by
    apply Finset.sum_le_sum
    intro e _
    linarith [cosine_quartic_upper (a * dot k (unitCarrierSeamDirection e))]
  have heq :
      (∑ e : Fin 30, ((a * dot k (unitCarrierSeamDirection e)) ^ 2 / 2 -
        (a * dot k (unitCarrierSeamDirection e)) ^ 4 / 24)) =
      a ^ 2 / 2 * (10 * radiusSquared k) -
        a ^ 4 / 24 * (6 * radiusSquared k ^ 2) := by
    rw [← unit_seam_second_moment_eq k, ← unit_seam_fourth_moment_eq k]
    rw [Finset.sum_sub_distrib, Finset.mul_sum, Finset.mul_sum]
    congr 1 <;> apply Finset.sum_congr rfl <;> intro e _ <;> ring
  rw [heq] at hsum
  have h := mul_le_mul_of_nonneg_left hsum
    (show 0 ≤ 1 / (5 * a ^ 2) from by positivity)
  have halg : (1 / (5 * a ^ 2)) *
      (a ^ 2 / 2 * (10 * radiusSquared k) -
        a ^ 4 / 24 * (6 * radiusSquared k ^ 2)) =
      radiusSquared k - a ^ 2 * radiusSquared k ^ 2 / 20 := by
    field_simp [ha]
    ring
  rw [halg] at h
  unfold edgeCurrentCharacterSymbol
  linarith

/-- Frequency error with no infrared division and no ultraviolet cutoff.
The bound is deliberately uniform rather than asymptotically optimal. -/
theorem exact_frequency_error {a : ℝ} (ha : a ≠ 0)
    (k : ThresholdVec3) :
    0 ≤ euclideanMomentumMagnitude k - photonModeFrequency a k ∧
    euclideanMomentumMagnitude k - photonModeFrequency a k ≤
      a ^ 2 * euclideanMomentumMagnitude k ^ 3 / 20 := by
  have hr := euclideanMomentumMagnitude_nonnegative k
  have hw := photonModeFrequency_nonnegative a k
  have hle := exact_fz12_frequency_le_euclideanMomentumMagnitude ha k
  refine ⟨sub_nonneg.mpr hle, ?_⟩
  have he := (exact_symbol_quadratic_error ha k).2
  rw [radiusSquared_eq_euclideanMomentumMagnitude_sq,
    ← photonModeFrequency_sq a k] at he
  by_cases hz : euclideanMomentumMagnitude k = 0
  · rw [hz] at hle ⊢
    have hwz := le_antisymm hle hw
    simp [hwz]
  · have hp : 0 < euclideanMomentumMagnitude k := lt_of_le_of_ne hr (Ne.symm hz)
    apply (mul_le_mul_iff_left₀ hp).mp
    nlinarith [mul_nonneg hw (sub_nonneg.mpr hle)]

/-- Entrywise cosine propagator control on the complete symbol. -/
theorem cosine_propagator_error {a : ℝ} (ha : a ≠ 0)
    (k : ThresholdVec3) (t : ℝ) :
    |Real.cos (t * photonModeFrequency a k) -
      Real.cos (t * euclideanMomentumMagnitude k)| ≤
      |t| * (a ^ 2 * euclideanMomentumMagnitude k ^ 3 / 20) := by
  calc
    _ ≤ |t * photonModeFrequency a k - t * euclideanMomentumMagnitude k| :=
      Real.abs_cos_sub_cos_le _ _
    _ = |t| * (euclideanMomentumMagnitude k - photonModeFrequency a k) := by
      rw [← mul_sub, abs_mul, abs_of_nonpos (sub_nonpos.mpr
        (exact_fz12_frequency_le_euclideanMomentumMagnitude ha k))]
      ring
    _ ≤ _ := mul_le_mul_of_nonneg_left (exact_frequency_error ha k).2 (abs_nonneg t)

/-- Entrywise sine propagator control, including the zero mode. -/
theorem sine_propagator_error {a : ℝ} (ha : a ≠ 0)
    (k : ThresholdVec3) (t : ℝ) :
    |Real.sin (t * photonModeFrequency a k) -
      Real.sin (t * euclideanMomentumMagnitude k)| ≤
      |t| * (a ^ 2 * euclideanMomentumMagnitude k ^ 3 / 20) := by
  calc
    _ ≤ |t * photonModeFrequency a k - t * euclideanMomentumMagnitude k| :=
      Real.abs_sin_sub_sin_le _ _
    _ = |t| * (euclideanMomentumMagnitude k - photonModeFrequency a k) := by
      rw [← mul_sub, abs_mul, abs_of_nonpos (sub_nonpos.mpr
        (exact_fz12_frequency_le_euclideanMomentumMagnitude ha k))]
      ring
    _ ≤ _ := mul_le_mul_of_nonneg_left (exact_frequency_error ha k).2 (abs_nonneg t)

end OPH.SeamMaxwellContinuum

#print axioms OPH.SeamMaxwellContinuum.cosine_quartic_upper
#print axioms OPH.SeamMaxwellContinuum.exact_symbol_quadratic_error
#print axioms OPH.SeamMaxwellContinuum.exact_frequency_error
#print axioms OPH.SeamMaxwellContinuum.cosine_propagator_error
#print axioms OPH.SeamMaxwellContinuum.sine_propagator_error
