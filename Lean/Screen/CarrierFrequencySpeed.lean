import SeamCurrentPhotonLeptonThreshold

open scoped BigOperators goldenRatio

namespace OPH.CarrierFrequencySpeed

open OPH.PrimitivePortFrameQuotient
open OPH.PrimitivePortTranslationBridge
open OPH.SeamCurrentEdge30Moment
open OPH.SeamCurrentHomogeneousAction
open OPH.SeamCurrentPhotonLeptonThreshold

/-!
# Exact speed control for normalized positive cosine carriers

This file isolates one target-free consequence of a positive finite cosine
symbol.  A weighted support whose second moment is a Euclidean tight frame has
an exact sine-feature realization.  After quadratic normalization, that
feature map is a contraction.  Its norm, the nonnegative square root of the
cosine symbol, is therefore globally one-Lipschitz.

The theorem is mathematical.  Calling the feature norm a physical frequency,
its parameter a physical time, or its Lipschitz bound a signal-speed bound
requires separate position, field, clock, wave-packet, frame, scale, and
readout attachments.  No such attachment is assumed or proved here.
-/

abbrev Vec3 := OPH.PrimitivePortTranslationBridge.Vec3

/-- The Euclidean `L2` realization of the three-coordinate carrier chart. -/
noncomputable def euclideanVec (x : Vec3) : EuclideanSpace ℝ (Fin 3) :=
  WithLp.toLp 2 x

/-- Euclidean magnitude in the selected carrier chart. -/
noncomputable def euclideanMagnitude (x : Vec3) : ℝ := ‖euclideanVec x‖

theorem euclideanMagnitude_sq (x : Vec3) :
    euclideanMagnitude x ^ 2 = dot x x := by
  unfold euclideanMagnitude euclideanVec dot
  rw [EuclideanSpace.real_norm_sq_eq]
  apply Finset.sum_congr rfl
  intro d _
  simp [sq]

/-- Positive weighted finite support with an exact Euclidean tight-frame
second moment.  `tight` is the frame constant before quadratic
normalization. -/
structure PositiveTightFrame (ι : Type*) [Fintype ι] where
  weight : ι → ℝ
  direction : ι → Vec3
  tight : ℝ
  weight_pos : ∀ i, 0 < weight i
  tight_pos : 0 < tight
  secondMoment : ∀ x : Vec3,
    (∑ i : ι, weight i * dot x (direction i) ^ 2) = tight * dot x x

variable {ι : Type*} [Fintype ι]

/-- Quadratically normalized positive cosine symbol. -/
noncomputable def cosineSymbol
    (frame : PositiveTightFrame ι) (a : ℝ) (k : Vec3) : ℝ :=
  (2 / (frame.tight * a ^ 2)) * ∑ i : ι,
    frame.weight i * (1 - Real.cos (a * dot k (frame.direction i)))

/-- Per-support sine amplitude used by the exact Hilbert-space embedding. -/
noncomputable def featureAmplitude
    (frame : PositiveTightFrame ι) (a : ℝ) (i : ι) : ℝ :=
  Real.sqrt (4 * frame.weight i / (a ^ 2 * frame.tight))

/-- Exact sine-feature realization of the normalized cosine symbol. -/
noncomputable def feature
    (frame : PositiveTightFrame ι) (a : ℝ) (k : Vec3) :
    EuclideanSpace ℝ ι :=
  WithLp.toLp 2 (fun i ↦
    featureAmplitude frame a i *
      Real.sin ((a * dot k (frame.direction i)) / 2))

/-- The feature amplitude squares to its declared positive weight factor. -/
theorem featureAmplitude_sq
    (frame : PositiveTightFrame ι) (a : ℝ) (i : ι) :
    featureAmplitude frame a i ^ 2 =
      4 * frame.weight i / (a ^ 2 * frame.tight) := by
  unfold featureAmplitude
  rw [Real.sq_sqrt]
  exact div_nonneg
    (mul_nonneg (by norm_num) (frame.weight_pos i).le)
    (mul_nonneg (sq_nonneg a) frame.tight_pos.le)

/-- The feature norm squared is exactly the normalized full cosine symbol. -/
theorem feature_norm_sq_eq_cosineSymbol
    (frame : PositiveTightFrame ι) {a : ℝ} (ha : a ≠ 0) (k : Vec3) :
    ‖feature frame a k‖ ^ 2 = cosineSymbol frame a k := by
  rw [EuclideanSpace.real_norm_sq_eq]
  unfold feature cosineSymbol
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  rw [mul_pow, featureAmplitude_sq frame a i]
  have htrig := Real.cos_two_mul_eq_one_sub
    ((a * dot k (frame.direction i)) / 2)
  have htwo : 2 * ((a * dot k (frame.direction i)) / 2) =
      a * dot k (frame.direction i) := by ring
  rw [htwo] at htrig
  field_simp [ha, ne_of_gt frame.tight_pos]
  rw [htrig]
  ring

/-- The nonnegative auxiliary frequency supplied by the exact feature norm. -/
noncomputable def frequency
    (frame : PositiveTightFrame ι) (a : ℝ) (k : Vec3) : ℝ :=
  ‖feature frame a k‖

theorem frequency_nonnegative
    (frame : PositiveTightFrame ι) (a : ℝ) (k : Vec3) :
    0 ≤ frequency frame a k := norm_nonneg _

/-- The auxiliary frequency squares to the full normalized cosine symbol. -/
theorem frequency_sq_eq_cosineSymbol
    (frame : PositiveTightFrame ι) {a : ℝ} (ha : a ≠ 0) (k : Vec3) :
    frequency frame a k ^ 2 = cosineSymbol frame a k :=
  feature_norm_sq_eq_cosineSymbol frame ha k

/-- Coordinatewise sine contraction after the exact positive amplitude is
inserted. -/
theorem feature_coordinate_difference_sq_le
    (frame : PositiveTightFrame ι) {a : ℝ} (ha : a ≠ 0)
    (k p : Vec3) (i : ι) :
    (feature frame a k i - feature frame a p i) ^ 2 ≤
      (frame.weight i / frame.tight) *
        dot (k - p) (frame.direction i) ^ 2 := by
  have hsin := Real.abs_sin_sub_sin_le
    ((a * dot k (frame.direction i)) / 2)
    ((a * dot p (frame.direction i)) / 2)
  have hsin_sq :
      (Real.sin ((a * dot k (frame.direction i)) / 2) -
          Real.sin ((a * dot p (frame.direction i)) / 2)) ^ 2 ≤
        ((a * dot k (frame.direction i)) / 2 -
          (a * dot p (frame.direction i)) / 2) ^ 2 := by
    have h := (sq_le_sq₀ (abs_nonneg _) (abs_nonneg _)).2 hsin
    simpa only [sq_abs] using h
  have hdot : dot (k - p) (frame.direction i) =
      dot k (frame.direction i) - dot p (frame.direction i) := by
    unfold dot
    simp only [Pi.sub_apply, sub_mul]
    rw [Finset.sum_sub_distrib]
  change
    (featureAmplitude frame a i *
        Real.sin ((a * dot k (frame.direction i)) / 2) -
      featureAmplitude frame a i *
        Real.sin ((a * dot p (frame.direction i)) / 2)) ^ 2 ≤ _
  calc
    _ = featureAmplitude frame a i ^ 2 *
        (Real.sin ((a * dot k (frame.direction i)) / 2) -
          Real.sin ((a * dot p (frame.direction i)) / 2)) ^ 2 := by ring
    _ ≤ featureAmplitude frame a i ^ 2 *
        ((a * dot k (frame.direction i)) / 2 -
          (a * dot p (frame.direction i)) / 2) ^ 2 :=
      mul_le_mul_of_nonneg_left hsin_sq (sq_nonneg _)
    _ = (frame.weight i / frame.tight) *
        dot (k - p) (frame.direction i) ^ 2 := by
      rw [featureAmplitude_sq frame a i, hdot]
      field_simp [ha, ne_of_gt frame.tight_pos]
      ring

/-- The complete sine-feature map is a Euclidean contraction. -/
theorem feature_dist_le
    (frame : PositiveTightFrame ι) {a : ℝ} (ha : a ≠ 0)
    (k p : Vec3) :
    ‖feature frame a k - feature frame a p‖ ≤ euclideanMagnitude (k - p) := by
  have hsum :
      (∑ i : ι, (feature frame a k i - feature frame a p i) ^ 2) ≤
        dot (k - p) (k - p) := by
    calc
      (∑ i : ι, (feature frame a k i - feature frame a p i) ^ 2) ≤
          ∑ i : ι, (frame.weight i / frame.tight) *
            dot (k - p) (frame.direction i) ^ 2 := by
        apply Finset.sum_le_sum
        intro i _
        exact feature_coordinate_difference_sq_le frame ha k p i
      _ = (1 / frame.tight) *
          ∑ i : ι, frame.weight i *
            dot (k - p) (frame.direction i) ^ 2 := by
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro i _
        ring
      _ = dot (k - p) (k - p) := by
        rw [frame.secondMoment]
        field_simp [ne_of_gt frame.tight_pos]
  apply (sq_le_sq₀ (norm_nonneg _)
    (show 0 ≤ euclideanMagnitude (k - p) from norm_nonneg _)).1
  rw [EuclideanSpace.real_norm_sq_eq, euclideanMagnitude_sq]
  simpa only [Pi.sub_apply, Real.norm_eq_abs, sq_abs] using hsum

/-- Global unit Lipschitz bound for the complete nonnegative frequency.
No Taylor expansion or angular orientation is used. -/
theorem frequency_global_one_lipschitz
    (frame : PositiveTightFrame ι) {a : ℝ} (ha : a ≠ 0)
    (k p : Vec3) :
    |frequency frame a k - frequency frame a p| ≤
      euclideanMagnitude (k - p) := by
  exact (abs_norm_sub_norm_le (feature frame a k) (feature frame a p)).trans
    (feature_dist_le frame ha k p)

/-! ## Exact FZ-12 and FZ-11 support instantiations -/

/-- The complete normalized thirty-seam support as a positive tight frame. -/
noncomputable def edge30Frame : PositiveTightFrame (Fin 30) where
  weight := fun _ ↦ 1
  direction := unitCarrierSeamDirection
  tight := 10
  weight_pos := by intro; norm_num
  tight_pos := by norm_num
  secondMoment := by
    intro x
    simpa using unit_seam_second_moment_eq x

theorem edge30_cosineSymbol_eq_fz12
    {a : ℝ} (ha : a ≠ 0) (k : Vec3) :
    cosineSymbol edge30Frame a k =
      OPH.SeamCurrentHomogeneousAction.edgeCurrentCharacterSymbol a k := by
  unfold cosineSymbol edge30Frame
    OPH.SeamCurrentHomogeneousAction.edgeCurrentCharacterSymbol
  field_simp [ha]
  ring

theorem fz12_frequency_global_one_lipschitz
    {a : ℝ} (ha : a ≠ 0) (k p : Vec3) :
    |frequency edge30Frame a k - frequency edge30Frame a p| ≤
      euclideanMagnitude (k - p) :=
  frequency_global_one_lipschitz edge30Frame ha k p

/-- Unit primitive-port directions in the same explicit frame as the
registered twelve-port carrier. -/
noncomputable def unitPortDirection (p : Fin 12) : Vec3 :=
  rawRadius⁻¹ • portVector p

/-- The complete normalized primitive-port support has tight-frame constant
four. -/
theorem unit_port_second_moment_eq (x : Vec3) :
    (∑ p : Fin 12, dot x (unitPortDirection p) ^ 2) = 4 * dot x x := by
  have hsqrt : Real.sqrt 5 ^ 2 = 5 := by norm_num
  have hrad : rawRadius ^ 2 = (5 + Real.sqrt 5) / 2 := by
    calc
      rawRadius ^ 2 = Real.goldenRatio + 2 := by
        unfold rawRadius
        exact Real.sq_sqrt (by linarith [Real.goldenRatio_pos])
      _ = (5 + Real.sqrt 5) / 2 := by
        unfold Real.goldenRatio
        ring
  unfold unitPortDirection portVector
  simp_rw [dot_smul_right]
  unfold dot
  simp [Fin.sum_univ_succ]
  field_simp [rawRadius_ne_zero]
  nlinarith

/-- The frozen primitive twelve-port support as a positive tight frame. -/
noncomputable def vertex12Frame : PositiveTightFrame (Fin 12) where
  weight := fun _ ↦ 1
  direction := unitPortDirection
  tight := 4
  weight_pos := by intro; norm_num
  tight_pos := by norm_num
  secondMoment := by
    intro x
    simpa using unit_port_second_moment_eq x

theorem vertex12_cosineSymbol_eq_fz11
    {a : ℝ} (ha : a ≠ 0) (k : Vec3) :
    cosineSymbol vertex12Frame a k =
      OPH.PrimitivePortTranslationBridge.cosineSymbol
        a k unitPortDirection := by
  unfold cosineSymbol vertex12Frame
    OPH.PrimitivePortTranslationBridge.cosineSymbol
    OPH.PrimitivePortTranslationBridge.portPhase
  field_simp [ha]
  ring

theorem fz11_frequency_global_one_lipschitz
    {a : ℝ} (ha : a ≠ 0) (k p : Vec3) :
    |frequency vertex12Frame a k - frequency vertex12Frame a p| ≤
      euclideanMagnitude (k - p) :=
  frequency_global_one_lipschitz vertex12Frame ha k p

end OPH.CarrierFrequencySpeed

/- Axiom audit: exact real inequalities and imported finite moment theorems. -/

#print axioms OPH.CarrierFrequencySpeed.frequency_global_one_lipschitz
#print axioms OPH.CarrierFrequencySpeed.fz12_frequency_global_one_lipschitz
#print axioms OPH.CarrierFrequencySpeed.fz11_frequency_global_one_lipschitz
