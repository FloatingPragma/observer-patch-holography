import SeamCurrentHomogeneousAction

open scoped BigOperators

namespace OPH.SeamCurrentDirichletGenerator

open OPH.PrimitivePortFrameQuotient
open OPH.SeamCurrentCarrierQuotient
open OPH.SeamCurrentEdge30Moment
open OPH.SeamCurrentHomogeneousAction

/-!
# The source-selected Markov and Dirichlet operator on the carrier completion

The exact seam-current quotient has a dense image in its response-Gram
rank-three completion.  The directed source seams give sixty translations of
that completion.  This file proves that the A2-natural, A3-unique move law
selects their equal average and that this average is a finite-range Markov
operator.  Its dimensionless generator `I - P` has a nonnegative local
Dirichlet square and satisfies the exact carré-du-champ identity.

This is a theorem about the response-selected carrier and its exact record
translations.  It supplies neither a physical field nor a clock, a
dimensionful length, a continuum action, or a laboratory readout.  Calling
the operator a physical spatial kinetic term requires those additional
identifications.
-/

/-- Product-coordinate presentation of the rank-three completion. -/
abbrev Vec3 := OPH.SeamCurrentEdge30Moment.Vec3

/-- Image of an exact even-sum record in the product-coordinate presentation
of the response-Gram completion. -/
noncomputable def completionPoint (p : D6Point) : Vec3 :=
  euclideanEquivVec3 (d6Position p)

/-- Translation vector in the completion supplied by one directed source
seam.  The common `rawRadius⁻¹` factor belongs to the proved response-Gram
metric; it is not a laboratory length. -/
noncomputable def completionSeamStep (e : DirectedSeam) : Vec3 :=
  rawRadius⁻¹ • directedSeamChartStep e

/-- The completion step is exactly the image of the corresponding cumulative
record displacement. -/
theorem completionSeamStep_eq_record_image (e : DirectedSeam) :
    completionSeamStep e = completionPoint (directedSeamStep e) := by
  change rawRadius⁻¹ • integerFrame (directedSeamAxisCurrent e) =
    euclideanEquivVec3
      (rawRadius⁻¹ • euclideanEquivVec3.symm
        (integerFrame (directedSeamAxisCurrent e)))
  rfl

/-- Addition of exact records becomes ordinary translation in the completion
chart. -/
theorem completionPoint_translate (e : DirectedSeam) (p : D6Point) :
    completionPoint (d6Translate (directedSeamStep e) p) =
      completionPoint p + completionSeamStep e := by
  rw [completionSeamStep_eq_record_image]
  unfold completionPoint
  rw [d6Translate, d6Position_add]
  ext d
  change
    (d6Position (directedSeamStep e)) d + (d6Position p) d =
      (d6Position p) d + (d6Position (directedSeamStep e)) d
  ring

/-- Reversal of a directed source seam reverses its completion translation. -/
theorem completionSeamStep_reverse (e : DirectedSeam) :
    completionSeamStep (reverseDirectedSeam e) = -completionSeamStep e := by
  unfold completionSeamStep directedSeamChartStep
  rw [directedSeamAxisCurrent_reverse, integerFrame_neg]
  simp

/-- Every raw directed seam chart displacement has squared norm four. -/
theorem directedSeamChartStep_norm_sq (e : DirectedSeam) :
    OPH.PrimitivePortTranslationBridge.dot
      (directedSeamChartStep e) (directedSeamChartStep e) = 4 := by
  classical
  let i : Fin 30 × Bool := directedSeamIndexEquiv.symm e
  have hi : indexedDirectedSeam i = e :=
    directedSeamIndexEquiv.apply_symm_apply e
  rw [← hi, directedSeamChartStep_indexed]
  cases i.2
  · simp only [Bool.false_eq_true, ↓reduceIte]
    exact carrierSeamDifference_norm_sq i.1
  · simp only [↓reduceIte]
    unfold OPH.PrimitivePortTranslationBridge.dot
    simp only [Pi.neg_apply, neg_mul_neg]
    exact carrierSeamDifference_norm_sq i.1

/-- The response-Gram completion step has one exact common squared norm.  The
raw integer-chart value four is therefore not the norm seen by the completion
metric. -/
theorem completionSeamStep_norm_sq (e : DirectedSeam) :
    OPH.PrimitivePortTranslationBridge.dot
      (completionSeamStep e) (completionSeamStep e) =
        2 - (2 / 5) * Real.sqrt 5 := by
  have hscale :
      OPH.PrimitivePortTranslationBridge.dot
          (completionSeamStep e) (completionSeamStep e) =
        rawRadius⁻¹ ^ 2 *
          OPH.PrimitivePortTranslationBridge.dot
            (directedSeamChartStep e) (directedSeamChartStep e) := by
    unfold completionSeamStep OPH.PrimitivePortTranslationBridge.dot
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro d _
    simp only [Pi.smul_apply, smul_eq_mul]
    ring
  have hsqrt : Real.sqrt 5 ^ 2 = 5 := by norm_num
  have hrad : rawRadius ^ 2 = (5 + Real.sqrt 5) / 2 := by
    calc
      rawRadius ^ 2 = Real.goldenRatio + 2 := by
        unfold rawRadius
        exact Real.sq_sqrt (by linarith [Real.goldenRatio_pos])
      _ = (5 + Real.sqrt 5) / 2 := by ring
  rw [hscale, directedSeamChartStep_norm_sq]
  field_simp [rawRadius_ne_zero]
  nlinarith

/-- Equivalent inverse-square-root presentation of the completion seam norm. -/
theorem completionSeamStep_norm_sq_eq_two_sub_two_div_sqrt_five
    (e : DirectedSeam) :
    OPH.PrimitivePortTranslationBridge.dot
      (completionSeamStep e) (completionSeamStep e) =
        2 - 2 / Real.sqrt 5 := by
  rw [completionSeamStep_norm_sq]
  have hsqrt : Real.sqrt 5 ^ 2 = 5 := by norm_num
  have hsqrt_ne : Real.sqrt 5 ≠ 0 := by positivity
  field_simp [hsqrt_ne]
  nlinarith

/-- The exact response-normalized seam squared norm lies strictly between one
and six fifths. -/
theorem completionSeamStep_norm_sq_bounds (e : DirectedSeam) :
    1 < OPH.PrimitivePortTranslationBridge.dot
        (completionSeamStep e) (completionSeamStep e) ∧
      OPH.PrimitivePortTranslationBridge.dot
        (completionSeamStep e) (completionSeamStep e) < 6 / 5 := by
  rw [completionSeamStep_norm_sq]
  have hsqrt_pos : 0 < Real.sqrt 5 := by positivity
  have hsqrt_sq : Real.sqrt 5 ^ 2 = 5 := by norm_num
  constructor <;> nlinarith

/-- Orientation reversal is an involutive equivalence of the exact sixty
source labels. -/
def reverseDirectedSeamEquiv : DirectedSeam ≃ DirectedSeam where
  toFun := reverseDirectedSeam
  invFun := reverseDirectedSeam
  left_inv := by
    intro e
    apply Subtype.ext
    exact Prod.swap_swap e.1
  right_inv := by
    intro e
    apply Subtype.ext
    exact Prod.swap_swap e.1

/-! ## Selected finite-range average -/

/-- Completion average associated with an arbitrary real directed-seam
weight law. -/
noncomputable def weightedCompletionAverage
    (weight : DirectedSeam → ℝ) (f : Vec3 → ℝ) : Vec3 → ℝ :=
  fun x ↦ ∑ e : DirectedSeam,
    weight e * f (x + completionSeamStep e)

/-- Equal source-counting average on the rank-three carrier completion. -/
noncomputable def completionMarkovAverage (f : Vec3 → ℝ) : Vec3 → ℝ :=
  weightedCompletionAverage (fun _ ↦ 1 / 60) f

/-- A2 naturality, A3 unique minimality, and source normalization remove every
directional coefficient from the completion average. -/
theorem a2a3_selected_completion_average_eq_markov
    (selection : A2A3DirectedSeamProjection) :
    weightedCompletionAverage selection.selected = completionMarkovAverage := by
  funext f x
  unfold weightedCompletionAverage completionMarkovAverage
  apply Finset.sum_congr rfl
  intro e _
  rw [a2a3_directed_seam_weight_eq_one_sixtieth selection e]

/-- The completion average restricts on the dense exact carrier to the same
source-counting translation operator. -/
theorem completionMarkovAverage_restricts_to_d6
    (f : Vec3 → ℝ) (p : D6Point) :
    completionMarkovAverage f (completionPoint p) =
      ∑ e : DirectedSeam, (1 / 60 : ℝ) *
        f (completionPoint (d6Translate (directedSeamStep e) p)) := by
  unfold completionMarkovAverage weightedCompletionAverage
  apply Finset.sum_congr rfl
  intro e _
  rw [completionPoint_translate]

/-- Complexification of the completion average agrees exactly with the
previously selected source-counting operator on every dense exact record. -/
theorem completionMarkovAverage_complex_restriction
    (f : Vec3 → ℝ) (p : D6Point) :
    (completionMarkovAverage f (completionPoint p) : ℂ) =
      sourceCountingEdgeCurrentOperator
        (fun q ↦ (f (completionPoint q) : ℂ)) p := by
  unfold completionMarkovAverage weightedCompletionAverage
    sourceCountingEdgeCurrentOperator weightedEdgeCurrentOperator
  push_cast
  apply Finset.sum_congr rfl
  intro e _
  rw [completionPoint_translate]

/-- The jump law is centrally symmetric.  This is the finite detailed-balance
ingredient; no measure on a physical space is assumed. -/
theorem completionMarkovAverage_eq_reverse_average
    (f : Vec3 → ℝ) (x : Vec3) :
    completionMarkovAverage f x =
      ∑ e : DirectedSeam, (1 / 60 : ℝ) *
        f (x - completionSeamStep e) := by
  unfold completionMarkovAverage weightedCompletionAverage
  calc
    (∑ e : DirectedSeam, (1 / 60 : ℝ) *
        f (x + completionSeamStep e)) =
        ∑ e : DirectedSeam, (1 / 60 : ℝ) *
          f (x + completionSeamStep (reverseDirectedSeam e)) := by
      exact (reverseDirectedSeamEquiv.bijective.sum_comp
        (fun e : DirectedSeam ↦
          (1 / 60 : ℝ) * f (x + completionSeamStep e))).symm
    _ = ∑ e : DirectedSeam, (1 / 60 : ℝ) *
        f (x - completionSeamStep e) := by
      apply Finset.sum_congr rfl
      intro e _
      rw [completionSeamStep_reverse]
      simp [sub_eq_add_neg]

/-! ## Exact detailed balance on the dense record carrier -/

/-- Following one directed seam and following its reverse undo each other on
the exact cumulative record carrier. -/
theorem d6Translate_directed_reverse_iff
    (e : DirectedSeam) (p q : D6Point) :
    d6Translate (directedSeamStep e) p = q ↔
      d6Translate (directedSeamStep (reverseDirectedSeam e)) q = p := by
  constructor
  · intro h
    ext i
    have hi := congrFun (congrArg D6Point.control h) i
    have hrev :
        (directedSeamStep (reverseDirectedSeam e)).control i =
          -(directedSeamStep e).control i := by
      exact congrFun (directedSeamAxisCurrent_reverse e) i
    simp only [d6Translate, d6Add] at hi ⊢
    rw [hrev]
    omega
  · intro h
    ext i
    have hi := congrFun (congrArg D6Point.control h) i
    have hrev :
        (directedSeamStep (reverseDirectedSeam e)).control i =
          -(directedSeamStep e).control i := by
      exact congrFun (directedSeamAxisCurrent_reverse e) i
    simp only [d6Translate, d6Add] at hi ⊢
    rw [hrev] at hi
    omega

/-- Number of directed source labels carrying exact record `p` to exact
record `q`. -/
noncomputable def d6TransitionMultiplicity (p q : D6Point) : ℕ := by
  classical
  exact ∑ e : DirectedSeam,
    if d6Translate (directedSeamStep e) p = q then 1 else 0

set_option maxHeartbeats 1000000 in
/-- Exact detailed balance for counting measure: the multiplicity of every
record transition equals that of its reverse. -/
theorem d6TransitionMultiplicity_symmetric (p q : D6Point) :
    d6TransitionMultiplicity p q = d6TransitionMultiplicity q p := by
  classical
  unfold d6TransitionMultiplicity
  calc
    (∑ e : DirectedSeam,
      if d6Translate (directedSeamStep e) p = q then 1 else 0) =
        ∑ e : DirectedSeam,
          if d6Translate
              (directedSeamStep (reverseDirectedSeam e)) p = q
            then 1 else 0 := by
      exact (reverseDirectedSeamEquiv.bijective.sum_comp
        (fun e : DirectedSeam ↦
          if d6Translate (directedSeamStep e) p = q then 1 else 0)).symm
    _ = ∑ e : DirectedSeam,
        if d6Translate (directedSeamStep e) q = p then 1 else 0 := by
      apply Finset.sum_congr rfl
      intro e _
      by_cases h : d6Translate (directedSeamStep e) q = p
      · have hr := (d6Translate_directed_reverse_iff e q p).1 h
        simp [h, hr]
      · have hr :
            d6Translate (directedSeamStep (reverseDirectedSeam e)) p ≠ q := by
          intro hr
          exact h ((d6Translate_directed_reverse_iff e q p).2 hr)
        simp [h, hr]

/-- Normalized exact transition weight on the dense record carrier. -/
noncomputable def d6TransitionWeight (p q : D6Point) : ℝ :=
  (1 / 60) * d6TransitionMultiplicity p q

/-- The normalized record transition kernel is reversible for counting
measure. -/
theorem d6TransitionWeight_symmetric (p q : D6Point) :
    d6TransitionWeight p q = d6TransitionWeight q p := by
  unfold d6TransitionWeight
  rw [d6TransitionMultiplicity_symmetric]

/-! ## Markov properties -/

/-- The equal average preserves constants exactly. -/
theorem completionMarkovAverage_const (c : ℝ) :
    completionMarkovAverage (fun _ ↦ c) = fun _ ↦ c := by
  funext x
  simp [completionMarkovAverage, weightedCompletionAverage,
    directedSeam_card]

/-- The equal average preserves pointwise nonnegativity. -/
theorem completionMarkovAverage_nonnegative
    {f : Vec3 → ℝ} (hf : ∀ x, 0 ≤ f x) :
    ∀ x, 0 ≤ completionMarkovAverage f x := by
  intro x
  unfold completionMarkovAverage weightedCompletionAverage
  exact Finset.sum_nonneg fun e _ ↦ mul_nonneg (by norm_num) (hf _)

/-- The equal average is order preserving. -/
theorem completionMarkovAverage_mono
    {f g : Vec3 → ℝ} (hfg : ∀ x, f x ≤ g x) :
    ∀ x, completionMarkovAverage f x ≤ completionMarkovAverage g x := by
  intro x
  unfold completionMarkovAverage weightedCompletionAverage
  apply Finset.sum_le_sum
  intro e _
  exact mul_le_mul_of_nonneg_left (hfg _) (by norm_num)

/-- The average commutes with all translations of the carrier completion. -/
theorem completionMarkovAverage_translation_covariant
    (f : Vec3 → ℝ) (a x : Vec3) :
    completionMarkovAverage (fun y ↦ f (a + y)) x =
      completionMarkovAverage f (a + x) := by
  unfold completionMarkovAverage weightedCompletionAverage
  apply Finset.sum_congr rfl
  intro e _
  congr 2
  apply congrArg f
  abel

/-- The finite source average preserves continuity on the carrier
completion. -/
theorem completionMarkovAverage_continuous
    {f : Vec3 → ℝ} (hf : Continuous f) :
    Continuous (completionMarkovAverage f) := by
  unfold completionMarkovAverage weightedCompletionAverage
  fun_prop

/-- The average at a point depends only on its sixty declared source
translates. -/
theorem completionMarkovAverage_local
    (f g : Vec3 → ℝ) (x : Vec3)
    (hlocal : ∀ e : DirectedSeam,
      f (x + completionSeamStep e) = g (x + completionSeamStep e)) :
    completionMarkovAverage f x = completionMarkovAverage g x := by
  unfold completionMarkovAverage weightedCompletionAverage
  apply Finset.sum_congr rfl
  intro e _
  rw [hlocal e]

/-! ## Exact Fourier character on the mathematical completion -/

/-- The exact conversion from unit seam directions to the response-Gram
completion steps.  It is a mathematical metric factor, not a physical
length. -/
noncomputable def completionStepScale : ℝ := 2 * rawRadius⁻¹

theorem completionStepScale_pos : 0 < completionStepScale := by
  unfold completionStepScale
  exact mul_pos (by norm_num) (inv_pos.mpr rawRadius_pos)

/-- Complex completion average associated with an arbitrary real directed-seam
weight law. -/
noncomputable def weightedComplexCompletionAverage
    (weight : DirectedSeam → ℝ) (f : Vec3 → ℂ) : Vec3 → ℂ :=
  fun x ↦ ∑ e : DirectedSeam,
    (weight e : ℂ) * f (x + completionSeamStep e)

/-- Complex version of the same equal source-counting completion average. -/
noncomputable def completionComplexMarkovAverage
    (f : Vec3 → ℂ) : Vec3 → ℂ :=
  weightedComplexCompletionAverage (fun _ ↦ 1 / 60) f

/-- The A2/A3-selected complex completion average is the equal source
average. -/
theorem a2a3_selected_complex_completion_average_eq_markov
    (selection : A2A3DirectedSeamProjection) :
    weightedComplexCompletionAverage selection.selected =
      completionComplexMarkovAverage := by
  funext f x
  unfold weightedComplexCompletionAverage completionComplexMarkovAverage
  apply Finset.sum_congr rfl
  intro e _
  rw [a2a3_directed_seam_weight_eq_one_sixtieth selection e]

/-- The exact completion average is the certified sixty-label chart
average at the response-Gram completion factor. -/
theorem completionComplexMarkovAverage_eq_sourceCountingChartAverage
    (f : Vec3 → ℂ) :
    completionComplexMarkovAverage f =
      sourceCountingChartAverage completionStepScale f := by
  funext x
  unfold completionComplexMarkovAverage weightedComplexCompletionAverage
    sourceCountingChartAverage
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro e _
  push_cast
  congr 1
  apply congrArg f
  congr 1
  unfold completionSeamStep completionStepScale
    unitDirectedSeamChartStep
  rw [smul_smul]
  congr 1
  ring

/-- Dimensionless complex generator `I - P` of the same selected completion
average. -/
noncomputable def completionComplexDirichletGenerator
    (f : Vec3 → ℂ) : Vec3 → ℂ :=
  fun x ↦ f x - completionComplexMarkovAverage f x

/-- Exact Fourier eigenvalue of the completion generator. -/
noncomputable def completionFourierSymbol (k : Vec3) : ℝ :=
  1 - (1 / 30) * ∑ e : Fin 30,
    Real.cos (completionStepScale *
      OPH.PrimitivePortTranslationBridge.dot k
        (unitCarrierSeamDirection e))

/-- The Fourier character of the canonical completion generator is
nonnegative. -/
theorem completionFourierSymbol_nonnegative (k : Vec3) :
    0 ≤ completionFourierSymbol k := by
  have hsum :
      (∑ e : Fin 30,
        Real.cos (completionStepScale *
          OPH.PrimitivePortTranslationBridge.dot k
            (unitCarrierSeamDirection e))) ≤ 30 := by
    calc
      (∑ e : Fin 30,
        Real.cos (completionStepScale *
          OPH.PrimitivePortTranslationBridge.dot k
            (unitCarrierSeamDirection e))) ≤
          ∑ _e : Fin 30, (1 : ℝ) := by
        apply Finset.sum_le_sum
        intro e _
        exact Real.cos_le_one _
      _ = 30 := by simp
  unfold completionFourierSymbol
  nlinarith

/-- Plane waves diagonalize the canonical completion generator itself.  No
auxiliary field operator is introduced at this mathematical level. -/
theorem completionComplexDirichletGenerator_planeWave
    (k x : Vec3) :
    completionComplexDirichletGenerator
        (OPH.PrimitivePortTranslationBridge.planeWave k) x =
      (completionFourierSymbol k : ℂ) *
        OPH.PrimitivePortTranslationBridge.planeWave k x := by
  unfold completionComplexDirichletGenerator
  rw [show completionComplexMarkovAverage
      (OPH.PrimitivePortTranslationBridge.planeWave k) x =
        sourceCountingChartAverage completionStepScale
          (OPH.PrimitivePortTranslationBridge.planeWave k) x by
    rw [completionComplexMarkovAverage_eq_sourceCountingChartAverage]]
  rw [sourceCountingChartAverage_planeWave]
  unfold completionFourierSymbol
  push_cast
  ring

/-- Under the declared quadratic normalization, the eigenvalue of the
canonical completion generator is exactly the FZ-12 edge-current character.
This is an equality of mathematical carrier symbols; physical field, time,
and scale remain outside it. -/
theorem completionFourierSymbol_normalizes_to_edgeCurrentCharacterSymbol
    (k : Vec3) :
    (6 / completionStepScale ^ 2) * completionFourierSymbol k =
      edgeCurrentCharacterSymbol completionStepScale k := by
  have hscale : completionStepScale ≠ 0 :=
    ne_of_gt completionStepScale_pos
  have hsum :
      (∑ e : Fin 30,
        (1 - Real.cos (completionStepScale *
          OPH.PrimitivePortTranslationBridge.dot k
            (unitCarrierSeamDirection e)))) =
        30 - ∑ e : Fin 30,
          Real.cos (completionStepScale *
            OPH.PrimitivePortTranslationBridge.dot k
              (unitCarrierSeamDirection e)) := by
    rw [Finset.sum_sub_distrib]
    simp
  unfold completionFourierSymbol edgeCurrentCharacterSymbol
  rw [hsum]
  field_simp [hscale]
  ring

/-- Pure coordinate dilation of the internal completion character.  The
parameter `a` is a supplied coordinate scale.  The equality theorem below
requires it to be nonzero; interpreting it as a length additionally requires
positivity.  This definition does not identify it with a physical ruler or the
result with a field frequency. -/
noncomputable def dilatedCompletionFourierSymbol
    (a : ℝ) (k : Vec3) : ℝ :=
  (6 / a ^ 2) *
    completionFourierSymbol ((a / completionStepScale) • k)

/-- The coordinate-dilated completion symbol is algebraically the edge-current
character at the supplied scale.  Physical position and
`edgeCurrentCharacterSymbol a k = omega_physical k ^ 2` remain separate
premises. -/
theorem dilatedCompletionFourierSymbol_eq_edgeCurrentCharacterSymbol
    (a : ℝ) (ha : a ≠ 0) (k : Vec3) :
    dilatedCompletionFourierSymbol a k =
      edgeCurrentCharacterSymbol a k := by
  have hstep : completionStepScale ≠ 0 :=
    ne_of_gt completionStepScale_pos
  have hangle (e : Fin 30) :
      completionStepScale *
          OPH.PrimitivePortTranslationBridge.dot
            ((a / completionStepScale) • k)
            (unitCarrierSeamDirection e) =
        a * OPH.PrimitivePortTranslationBridge.dot k
          (unitCarrierSeamDirection e) := by
    have hdot :
        OPH.PrimitivePortTranslationBridge.dot
            ((a / completionStepScale) • k)
            (unitCarrierSeamDirection e) =
          (a / completionStepScale) *
            OPH.PrimitivePortTranslationBridge.dot k
              (unitCarrierSeamDirection e) := by
      unfold OPH.PrimitivePortTranslationBridge.dot
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro d _
      simp
      ring
    rw [hdot]
    field_simp [hstep]
  have hsum :
      (∑ e : Fin 30,
        (1 - Real.cos (a *
          OPH.PrimitivePortTranslationBridge.dot k
            (unitCarrierSeamDirection e)))) =
        30 - ∑ e : Fin 30,
          Real.cos (a *
            OPH.PrimitivePortTranslationBridge.dot k
              (unitCarrierSeamDirection e)) := by
    rw [Finset.sum_sub_distrib]
    simp
  unfold dilatedCompletionFourierSymbol completionFourierSymbol
    edgeCurrentCharacterSymbol
  simp_rw [hangle]
  rw [hsum]
  field_simp [ha]
  ring

/-- Composed mathematical record: the normalized eigenvalue of the canonical
completion generator is the edge-current character, and the independently
proved finite-moment ratios give the exact edge-30 ray.  This conjunction does
not itself formalize the Taylor-derivative bridge between the two statements. -/
theorem completionDirichletGenerator_and_fz12_ray_certificate :
    (∀ k : Vec3,
      (6 / completionStepScale ^ 2) * completionFourierSymbol k =
        edgeCurrentCharacterSymbol completionStepScale k) ∧
    seamLambdaB0 / seamLambdaC4 ^ 2 = 10 / 21 ∧
    seamLambdaB6 / seamLambdaC4 ^ 2 = -2 / 63 ∧
    seamLambdaB6 / seamLambdaB0 = -1 / 15 := by
  refine ⟨completionFourierSymbol_normalizes_to_edgeCurrentCharacterSymbol, ?_⟩
  exact seam_edge30_ray

/-! ## Canonical dimensionless Dirichlet generator -/

/-- Dimensionless carrier generator selected by the exact source average. -/
noncomputable def completionDirichletGenerator
    (f : Vec3 → ℝ) : Vec3 → ℝ :=
  fun x ↦ f x - completionMarkovAverage f x

/-- Generator built directly from a supplied A2/A3 projection. -/
noncomputable def selectedCompletionGenerator
    (selection : A2A3DirectedSeamProjection)
    (f : Vec3 → ℝ) : Vec3 → ℝ :=
  fun x ↦ f x - weightedCompletionAverage selection.selected f x

/-- Every admitted A2/A3 projection gives the same dimensionless carrier
generator. -/
theorem a2a3_selected_completion_generator_eq_dirichlet
    (selection : A2A3DirectedSeamProjection) :
    selectedCompletionGenerator selection = completionDirichletGenerator := by
  funext f x
  unfold selectedCompletionGenerator completionDirichletGenerator
  rw [a2a3_selected_completion_average_eq_markov selection]

/-- Constants lie in the kernel of the selected carrier generator. -/
theorem completionDirichletGenerator_const (c : ℝ) :
    completionDirichletGenerator (fun _ ↦ c) = 0 := by
  funext x
  simp [completionDirichletGenerator, completionMarkovAverage_const]

/-- The selected dimensionless generator preserves continuity. -/
theorem completionDirichletGenerator_continuous
    {f : Vec3 → ℝ} (hf : Continuous f) :
    Continuous (completionDirichletGenerator f) := by
  unfold completionDirichletGenerator
  exact hf.sub (completionMarkovAverage_continuous hf)

/-- Positive maximum principle for the selected generator `I - P`. -/
theorem completionDirichletGenerator_nonnegative_at_global_maximum
    (f : Vec3 → ℝ) (x : Vec3)
    (hmax : ∀ y, f y ≤ f x) :
    0 ≤ completionDirichletGenerator f x := by
  unfold completionDirichletGenerator completionMarkovAverage
    weightedCompletionAverage
  have hsum :
      (∑ e : DirectedSeam, (1 / 60 : ℝ) *
        f (x + completionSeamStep e)) ≤
        ∑ _e : DirectedSeam, (1 / 60 : ℝ) * f x := by
    apply Finset.sum_le_sum
    intro e _
    exact mul_le_mul_of_nonneg_left (hmax _) (by norm_num)
  simpa [directedSeam_card] using sub_nonneg.mpr hsum

/-- Pointwise local Dirichlet square of the symmetric finite jump law. -/
noncomputable def completionCarreDuChamp
    (f : Vec3 → ℝ) (x : Vec3) : ℝ :=
  (1 / 120) * ∑ e : DirectedSeam,
    (f (x + completionSeamStep e) - f x) ^ 2

/-- The local Dirichlet square is nonnegative. -/
theorem completionCarreDuChamp_nonnegative
    (f : Vec3 → ℝ) (x : Vec3) :
    0 ≤ completionCarreDuChamp f x := by
  unfold completionCarreDuChamp
  exact mul_nonneg (by norm_num) (Finset.sum_nonneg fun e _ ↦ sq_nonneg _)

/-- Exact carré-du-champ identity for the selected dimensionless generator.
It is the local algebraic Dirichlet identity, independent of a continuum
measure or clock. -/
theorem completionDirichletGenerator_carre_du_champ
    (f : Vec3 → ℝ) (x : Vec3) :
    2 * f x * completionDirichletGenerator f x -
        completionDirichletGenerator (fun y ↦ f y ^ 2) x =
      2 * completionCarreDuChamp f x := by
  let S : ℝ := ∑ e : DirectedSeam, f (x + completionSeamStep e)
  let Q : ℝ := ∑ e : DirectedSeam, f (x + completionSeamStep e) ^ 2
  have haverage : completionMarkovAverage f x = (1 / 60) * S := by
    unfold completionMarkovAverage weightedCompletionAverage
    rw [Finset.mul_sum]
  have haverageSq :
      completionMarkovAverage (fun y ↦ f y ^ 2) x = (1 / 60) * Q := by
    unfold completionMarkovAverage weightedCompletionAverage
    rw [Finset.mul_sum]
  have hdiff :
      (∑ e : DirectedSeam,
        (f (x + completionSeamStep e) - f x) ^ 2) =
        Q - 2 * f x * S + 60 * f x ^ 2 := by
    simp_rw [sub_sq]
    simp only [Finset.sum_sub_distrib, Finset.sum_add_distrib]
    simp only [Finset.sum_const, Finset.card_univ, directedSeam_card]
    rw [← Finset.sum_mul]
    rw [← Finset.mul_sum]
    dsimp only [S, Q]
    ring
  unfold completionDirichletGenerator completionCarreDuChamp
  rw [haverage, haverageSq, hdiff]
  ring

/-! ## Axiom audit

The theorems use the finite source orbit and ordinary Mathlib algebra.  The
A2/A3 conclusion consumes the named projection premises from
`SeamCurrentHomogeneousAction`; it does not promote them to axioms or infer a
physical interpretation.
-/

#print axioms completionSeamStep_eq_record_image
#print axioms completionPoint_translate
#print axioms directedSeamChartStep_norm_sq
#print axioms completionSeamStep_norm_sq
#print axioms completionSeamStep_norm_sq_eq_two_sub_two_div_sqrt_five
#print axioms completionSeamStep_norm_sq_bounds
#print axioms completionMarkovAverage_restricts_to_d6
#print axioms completionMarkovAverage_complex_restriction
#print axioms completionMarkovAverage_eq_reverse_average
#print axioms d6Translate_directed_reverse_iff
#print axioms d6TransitionMultiplicity_symmetric
#print axioms d6TransitionWeight_symmetric
#print axioms a2a3_selected_completion_average_eq_markov
#print axioms completionMarkovAverage_const
#print axioms completionMarkovAverage_nonnegative
#print axioms completionMarkovAverage_translation_covariant
#print axioms completionMarkovAverage_continuous
#print axioms completionMarkovAverage_local
#print axioms a2a3_selected_complex_completion_average_eq_markov
#print axioms completionComplexMarkovAverage_eq_sourceCountingChartAverage
#print axioms completionFourierSymbol_nonnegative
#print axioms completionComplexDirichletGenerator_planeWave
#print axioms completionFourierSymbol_normalizes_to_edgeCurrentCharacterSymbol
#print axioms dilatedCompletionFourierSymbol_eq_edgeCurrentCharacterSymbol
#print axioms completionDirichletGenerator_and_fz12_ray_certificate
#print axioms a2a3_selected_completion_generator_eq_dirichlet
#print axioms completionDirichletGenerator_continuous
#print axioms completionDirichletGenerator_nonnegative_at_global_maximum
#print axioms completionCarreDuChamp_nonnegative
#print axioms completionDirichletGenerator_carre_du_champ

end OPH.SeamCurrentDirichletGenerator
