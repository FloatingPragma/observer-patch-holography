import SeamCurrentEdge30Moment
import PortGramA5Isometry
import ObserverPatchHolography.EqualSeamSelection

open scoped BigOperators

namespace OPH.SeamCurrentHomogeneousAction

open OPH.PortFrameGram
open OPH.PortGramA5Isometry
open OPH.PrimitivePortFrameQuotient
open OPH.RepairWordCarrierReadout
open OPH.SeamCurrentCarrierQuotient
open OPH.SeamCurrentEdge30Moment
open ObserverPatchHolography.EqualSeamSelection

/-- Unambiguous local name for the common three-coordinate response chart. -/
abbrev Vec3 := OPH.SeamCurrentEdge30Moment.Vec3

/-!
# Homogeneous action of the exact seam-current carrier

The source incidence packet has sixty directed seams.  The registered proper
carrier rotations act transitively on this finite alphabet.  A directed seam
adds its signed boundary readout to the exact `D6` record carrier.  Record
addition gives a simply transitive action on `D6Point`, acts by isometries for
the induced pullback of the response-selected rank-three Gram metric, and
extends as ordinary translation on that metric's Euclidean completion.  This
is not the usual discrete, word, or six-dimensional lattice metric on `D6`.

The finite conditional A2/A3 theorem in the second half removes the remaining
directional weight choice.  If the feasible move laws and objective are
natural under the registered rotations and A3 supplies a unique normalized
minimizer, the sixty directed seams have weight `1/60`.  Their induced
operator is therefore the source-counting homogeneous convolution on the
exact record carrier.  Feasibility, objective naturality, and unique
minimality are ordinary theorem premises, not consequences asserted merely
from the names A1, A2, and A3.

This closes a mathematical action bridge.  It does not identify the record
completion with physical position, a seam event with physical field
propagation, or the source-counting move simplex with the feasible physical
history simplex.  Those semantic and temporal-completeness statements remain
ordinary premises.  The canonical `CarrierShadow` permits inequivalent repair
maps, as proved in `EqualSeamSelection`, so the present three typed shadows do
not discharge them by themselves.
-/

/-! ## The sixty directed source seams and their proper-carrier action -/

/-- An ordered incidence seam.  Both orientations are retained as separate
source events. -/
abbrev DirectedSeam :=
  {p : Fin 12 × Fin 12 // adj p.1 p.2 = true}

instance : Nonempty DirectedSeam :=
  ⟨⟨(0, 1), by decide⟩⟩

/-- The exact directed seam alphabet has sixty elements. -/
theorem directedSeam_card : Fintype.card DirectedSeam = 60 := by
  native_decide

set_option maxHeartbeats 4000000 in
/-- Every registered proper rotation preserves incidence. -/
theorem portMap_preserves_adj :
    ∀ r : ProperRotation, ∀ i j : Fin 12,
      adj (portMap r i) (portMap r j) = adj i j := by
  native_decide

/-- Relabel a directed seam by one registered proper carrier rotation. -/
def rotateDirectedSeam (r : ProperRotation) (e : DirectedSeam) :
    DirectedSeam :=
  ⟨(portMap r e.1.1, portMap r e.1.2), by
    rw [portMap_preserves_adj]
    exact e.2⟩

set_option maxHeartbeats 8000000 in
/-- The listed proper rotations act transitively on directed seams. -/
theorem proper_rotations_transitive_on_directed_seams :
    ∀ e f : DirectedSeam,
      ∃ r : ProperRotation, rotateDirectedSeam r e = f := by
  native_decide

set_option maxHeartbeats 8000000 in
/-- At each directed seam, distinct registered rotations have distinct
images.  This finite check avoids enumerating a redundant target variable. -/
theorem proper_rotation_source_injective :
    ∀ e : DirectedSeam,
      Function.Injective (fun r : ProperRotation ↦
        rotateDirectedSeam r e) := by
  native_decide

/-- A registered row carrying one directed seam to a declared target is
unique. -/
theorem proper_rotation_target_unique
    (e f : DirectedSeam) (r s : ProperRotation)
    (hr : rotateDirectedSeam r e = f)
    (hs : rotateDirectedSeam s e = f) : r = s := by
  apply proper_rotation_source_injective e
  exact hr.trans hs.symm

/-- The listed proper rotations act simply transitively on directed seams:
for each ordered source and target seam there is exactly one registered row
carrying the first to the second. -/
theorem proper_rotations_simply_transitive_on_directed_seams
    (e f : DirectedSeam) :
    ∃! r : ProperRotation, rotateDirectedSeam r e = f := by
  obtain ⟨r, hr⟩ := proper_rotations_transitive_on_directed_seams e f
  exact ⟨r, hr, fun s hs ↦ proper_rotation_target_unique e f s r hs hr⟩

/-! ## Additive `D6` record carrier -/

/-- Zero cumulative even-sum record. -/
def d6Zero : D6Point :=
  ⟨fun _ ↦ 0, by
    refine ⟨0, ?_⟩
    simp⟩

/-- Addition of cumulative even-sum records. -/
def d6Add (p q : D6Point) : D6Point :=
  ⟨fun i ↦ p.control i + q.control i, by
    rcases p.even_total with ⟨a, ha⟩
    rcases q.even_total with ⟨b, hb⟩
    refine ⟨a + b, ?_⟩
    simp only [Finset.sum_add_distrib] at ha hb ⊢
    omega⟩

/-- Additive inverse of a cumulative even-sum record. -/
def d6Neg (p : D6Point) : D6Point :=
  ⟨fun i ↦ -p.control i, by
    rcases p.even_total with ⟨a, ha⟩
    refine ⟨-a, ?_⟩
    simp only [Finset.sum_neg_distrib] at ha ⊢
    omega⟩

/-- Difference of cumulative even-sum records. -/
def d6Sub (p q : D6Point) : D6Point := d6Add p (d6Neg q)

theorem d6Add_assoc (p q r : D6Point) :
    d6Add (d6Add p q) r = d6Add p (d6Add q r) := by
  ext i
  simp [d6Add, add_assoc]

theorem d6Add_comm (p q : D6Point) : d6Add p q = d6Add q p := by
  ext i
  simp [d6Add, add_comm]

theorem d6Add_zero (p : D6Point) : d6Add p d6Zero = p := by
  ext i
  simp [d6Add, d6Zero]

theorem d6Add_neg (p : D6Point) : d6Add p (d6Neg p) = d6Zero := by
  ext i
  simp [d6Add, d6Neg, d6Zero]

/-- Regular translation by one cumulative record. -/
def d6Translate (a p : D6Point) : D6Point := d6Add a p

/-- The additive record action is transitive. -/
theorem d6Translate_transitive (p q : D6Point) :
    ∃ a : D6Point, d6Translate a p = q := by
  refine ⟨d6Sub q p, ?_⟩
  ext i
  simp [d6Translate, d6Sub, d6Add, d6Neg]

/-- The additive record action is free. -/
theorem d6Translate_free (a b p : D6Point)
    (h : d6Translate a p = d6Translate b p) : a = b := by
  ext i
  have hi := congrFun (congrArg D6Point.control h) i
  simp only [d6Translate, d6Add] at hi
  omega

/-- The action is simply transitive: there is exactly one cumulative record
translation carrying one record point to another. -/
theorem existsUnique_d6Translate (p q : D6Point) :
    ∃! a : D6Point, d6Translate a p = q := by
  obtain ⟨a, ha⟩ := d6Translate_transitive p q
  exact ⟨a, ha, fun b hb ↦ d6Translate_free b a p (hb.trans ha.symm)⟩

/-! ## The induced pullback metric and its completion carry the action -/

/-- The exact Gram-frame position map is additive. -/
theorem d6Position_add (p q : D6Point) :
    d6Position (d6Add p q) = d6Position p + d6Position q := by
  ext d
  simp [d6Position, d6IntegerPoint, pointEuclideanFrame, pointFrame,
    integerFrame, frameMap, castIntegerControl, d6Add]
  rw [← mul_add, ← Finset.sum_add_distrib]
  apply congrArg (fun z : ℝ ↦ rawRadius⁻¹ * z)
  apply Finset.sum_congr rfl
  intro x _
  ring_nf

/-- Every exact record translation is an isometry for the response-Gram
metric pulled back along `d6Position`.  No claim is made for the standard
discrete, word, or six-dimensional lattice metric on the same records. -/
theorem d6Translate_isometry (a : D6Point) : Isometry (d6Translate a) := by
  apply Isometry.of_dist_eq
  intro p q
  calc
    dist (d6Translate a p) (d6Translate a q) =
        dist (d6Position (d6Translate a p))
          (d6Position (d6Translate a q)) :=
      (d6Position_isometry.dist_eq _ _).symm
    _ = dist (d6Position a + d6Position p)
        (d6Position a + d6Position q) := by
      simp only [d6Translate, d6Position_add]
    _ = dist (d6Position p) (d6Position q) :=
      dist_add_left (d6Position a) (d6Position p) (d6Position q)
    _ = dist p q := d6Position_isometry.dist_eq p q

/-- Translation on the Euclidean completion of the induced pullback metric
by an exact `D6` record. -/
noncomputable def completedTranslation
    (a : D6Point) (x : EuclideanVec3) : EuclideanVec3 :=
  d6Position a + x

/-- Dense-record translation and completion translation commute exactly. -/
theorem completedTranslation_intertwines (a p : D6Point) :
    completedTranslation a (d6Position p) =
      d6Position (d6Translate a p) := by
  rw [d6Translate, d6Position_add]
  rfl

/-- Exact-record translations are isometries of the Euclidean completion. -/
theorem completedTranslation_isometry (a : D6Point) :
    Isometry (completedTranslation a) := by
  apply Isometry.of_dist_eq
  intro x y
  exact dist_add_left (d6Position a) x y

/-! ## Directed seam steps in the homogeneous carrier -/

/-- Signed six-axis record created by one directed seam event. -/
def directedSeamAxisCurrent (e : DirectedSeam) : Fin 6 → ℤ :=
  signedLoadControl (directedBoundary e.1.1 e.1.2 1)

set_option maxHeartbeats 1000000 in
/-- Every directed seam event lies in the even-sum `D6` carrier. -/
theorem directedSeamAxisCurrent_even :
    ∀ e : DirectedSeam,
      EvenAxisTotal (directedSeamAxisCurrent e) := by
  rintro ⟨⟨u, v⟩, huv⟩
  have hne : u ≠ v := by
    intro huvEq
    subst v
    rw [adj_irrefl] at huv
    contradiction
  rcases lt_or_gt_of_ne hne with hlt | hgt
  · obtain ⟨e, hleft, hright⟩ :=
      (seam_table_complete u v).1 ⟨hlt, huv⟩
    change EvenAxisTotal
      (signedLoadControl (directedBoundary u v 1))
    rw [← hleft, ← hright, ← seamBoundary_atom]
    exact seamAxisCurrent_even (seamAtom e 1)
  · have hvu : adj v u = true := by
      simpa [adj_symm] using huv
    obtain ⟨e, hleft, hright⟩ :=
      (seam_table_complete v u).1 ⟨hgt, hvu⟩
    change EvenAxisTotal
      (signedLoadControl (directedBoundary u v 1))
    rw [← directedBoundary_reverse_neg u v 1,
      ← hleft, ← hright, ← seamBoundary_atom]
    exact seamAxisCurrent_even (seamAtom e (-1))

/-- Exact record displacement attached to one directed seam event. -/
def directedSeamStep (e : DirectedSeam) : D6Point :=
  ⟨directedSeamAxisCurrent e, directedSeamAxisCurrent_even e⟩

/-- Canonical smaller-to-larger orientation of a seam-table row. -/
def canonicalDirectedSeam (e : Fin 30) : DirectedSeam :=
  ⟨(seamLeft e, seamRight e), (seam_table_sound e).2⟩

/-- Reverse the orientation of a directed seam. -/
def reverseDirectedSeam (e : DirectedSeam) : DirectedSeam :=
  ⟨(e.1.2, e.1.1), by simpa [adj_symm] using e.2⟩

/-- The exact sixty-label presentation: thirty canonical seam-table rows and
their reverses. -/
def indexedDirectedSeam (e : Fin 30 × Bool) : DirectedSeam :=
  if e.2 then reverseDirectedSeam (canonicalDirectedSeam e.1)
  else canonicalDirectedSeam e.1

/-- The indexed presentation has no duplicate directed seam. -/
theorem indexedDirectedSeam_injective :
    Function.Injective indexedDirectedSeam := by
  rintro ⟨e, b⟩ ⟨f, c⟩ h
  cases b <;> cases c
  · have hp : (seamLeft e, seamRight e) =
        (seamLeft f, seamRight f) := congrArg Subtype.val h
    have hef : e = f := seam_table_injective hp
    subst f
    rfl
  · have hp : (seamLeft e, seamRight e) =
        (seamRight f, seamLeft f) := congrArg Subtype.val h
    have heLt := (seam_table_sound e).1
    have hfLt := (seam_table_sound f).1
    have hfst : seamLeft e = seamRight f := by
      simpa using congrArg Prod.fst hp
    have hsnd : seamRight e = seamLeft f := by
      simpa using congrArg Prod.snd hp
    omega
  · have hp : (seamRight e, seamLeft e) =
        (seamLeft f, seamRight f) := congrArg Subtype.val h
    have heLt := (seam_table_sound e).1
    have hfLt := (seam_table_sound f).1
    have hfst : seamRight e = seamLeft f := by
      simpa using congrArg Prod.fst hp
    have hsnd : seamLeft e = seamRight f := by
      simpa using congrArg Prod.snd hp
    omega
  · have hp : (seamRight e, seamLeft e) =
        (seamRight f, seamLeft f) := congrArg Subtype.val h
    have horiginal : (seamLeft e, seamRight e) =
        (seamLeft f, seamRight f) := by
      apply Prod.ext
      · exact congrArg Prod.snd hp
      · exact congrArg Prod.fst hp
    have hef : e = f := seam_table_injective horiginal
    subst f
    rfl

/-- The indexed presentation contains every directed seam exactly once. -/
theorem indexedDirectedSeam_bijective :
    Function.Bijective indexedDirectedSeam := by
  apply (Fintype.bijective_iff_injective_and_card indexedDirectedSeam).2
  exact ⟨indexedDirectedSeam_injective, by
    simp [directedSeam_card]⟩

/-- Exact equivalence between sixty directed labels and the incidence
subtype. -/
noncomputable def directedSeamIndexEquiv :
    (Fin 30 × Bool) ≃ DirectedSeam :=
  Equiv.ofBijective indexedDirectedSeam indexedDirectedSeam_bijective

/-- The directed-seam step agrees exactly with the independently proved
thirty-current boundary/readout composite. -/
theorem directedSeamStep_canonical_control (e : Fin 30) :
    (directedSeamStep (canonicalDirectedSeam e)).control =
      seamAxisCurrent (seamAtom e 1) := by
  change signedLoadControl
      (directedBoundary (seamLeft e) (seamRight e) 1) =
    seamAxisCurrent (seamAtom e 1)
  rw [seamAxisCurrent, seamBoundary_atom]

/-- Reversing a seam reverses its signed six-axis record. -/
theorem directedSeamAxisCurrent_reverse (e : DirectedSeam) :
    directedSeamAxisCurrent (reverseDirectedSeam e) =
      fun i ↦ -directedSeamAxisCurrent e i := by
  funext i
  simp [directedSeamAxisCurrent, reverseDirectedSeam, signedLoadControl,
    directedBoundary]
  ring_nf

/-- The integer carrier chart is additive under reversal. -/
theorem integerFrame_neg (z : Fin 6 → ℤ) :
    integerFrame (fun i ↦ -z i) = -integerFrame z := by
  funext d
  simp [integerFrame, frameMap, castIntegerControl]

/-- Raw response chart of one directed seam displacement. -/
noncomputable def directedSeamChartStep (e : DirectedSeam) : Vec3 :=
  integerFrame (directedSeamAxisCurrent e)

/-- The canonical orientation gives the source seam difference exactly. -/
theorem directedSeamChartStep_canonical (e : Fin 30) :
    directedSeamChartStep (canonicalDirectedSeam e) =
      carrierSeamDifference e := by
  unfold directedSeamChartStep carrierSeamDifference
  apply congrArg integerFrame
  exact directedSeamStep_canonical_control e

/-- The reversed orientation gives the negative source seam difference. -/
theorem directedSeamChartStep_reverse_canonical (e : Fin 30) :
    directedSeamChartStep
        (reverseDirectedSeam (canonicalDirectedSeam e)) =
      -carrierSeamDifference e := by
  unfold directedSeamChartStep
  rw [directedSeamAxisCurrent_reverse, integerFrame_neg]
  simpa [directedSeamChartStep] using
    congrArg Neg.neg (directedSeamChartStep_canonical e)

/-- Every one of the sixty directed labels is one of the thirty canonical
source differences or its negative. -/
theorem directedSeamChartStep_indexed (e : Fin 30 × Bool) :
    directedSeamChartStep (indexedDirectedSeam e) =
      if e.2 then -carrierSeamDifference e.1
      else carrierSeamDifference e.1 := by
  cases e with
  | mk i reversed =>
      cases reversed <;>
        simp [indexedDirectedSeam, directedSeamChartStep_canonical,
          directedSeamChartStep_reverse_canonical]

/-- Unit direction used by the frozen edge-current moments.  The raw
integer-record seam displacement has squared norm four, so physical scale
normalization must divide it by two before the Taylor coefficients are read. -/
noncomputable def unitCarrierSeamDirection (e : Fin 30) : Vec3 :=
  (1 / 2 : ℝ) • carrierSeamDifference e

/-- Every normalized seam direction has squared norm one. -/
theorem unitCarrierSeamDirection_norm_sq (e : Fin 30) :
    OPH.PrimitivePortTranslationBridge.dot
      (unitCarrierSeamDirection e) (unitCarrierSeamDirection e) = 1 := by
  calc
    OPH.PrimitivePortTranslationBridge.dot
        (unitCarrierSeamDirection e) (unitCarrierSeamDirection e) =
      (1 / 4 : ℝ) * OPH.PrimitivePortTranslationBridge.dot
        (carrierSeamDifference e) (carrierSeamDifference e) := by
      unfold unitCarrierSeamDirection OPH.PrimitivePortTranslationBridge.dot
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro d _
      simp
      ring
    _ = 1 := by
      rw [carrierSeamDifference_norm_sq]
      norm_num

/-- Unit-normalized chart step for one directed seam label. -/
noncomputable def unitDirectedSeamChartStep (e : DirectedSeam) : Vec3 :=
  (1 / 2 : ℝ) • directedSeamChartStep e

/-- The sixty normalized directed labels are the thirty unit seam directions
and their negatives. -/
theorem unitDirectedSeamChartStep_indexed (e : Fin 30 × Bool) :
    unitDirectedSeamChartStep (indexedDirectedSeam e) =
      if e.2 then -unitCarrierSeamDirection e.1
      else unitCarrierSeamDirection e.1 := by
  rw [unitDirectedSeamChartStep, directedSeamChartStep_indexed]
  cases e.2 <;> simp [unitCarrierSeamDirection]

/-- In the finite source chart, the canonical directed seam step is the same
carrier difference whose projective support gives the edge-30 moments. -/
theorem directedSeamStep_canonical_chart (e : Fin 30) :
    euclideanEquivVec3
        (d6Position (directedSeamStep (canonicalDirectedSeam e))) =
      rawRadius⁻¹ • carrierSeamDifference e := by
  change rawRadius⁻¹ •
      integerFrame (directedSeamStep (canonicalDirectedSeam e)).control =
    rawRadius⁻¹ • carrierSeamDifference e
  rw [directedSeamStep_canonical_control]
  rfl

/-! ## A2/A3 source-counting selection and the unique homogeneous operator -/

/-- A2 naturality and A3 uniqueness on the exact directed-seam alphabet.
Geometry fixes `rotateDirectedSeam`; callers supply only the feasible move
simplex, its natural objective, and its unique selected minimizer. -/
structure A2A3DirectedSeamProjection where
  feasible : Set (DirectedSeam → ℝ)
  objective : (DirectedSeam → ℝ) → ℝ
  feasible_natural :
    ∀ r, ∀ p ∈ feasible, (p ∘ rotateDirectedSeam r) ∈ feasible
  objective_natural :
    ∀ r, ∀ p ∈ feasible,
      objective (p ∘ rotateDirectedSeam r) = objective p
  selected : DirectedSeam → ℝ
  selected_feasible : selected ∈ feasible
  selected_minimal : ∀ p ∈ feasible, objective selected ≤ objective p
  selected_unique :
    ∀ p ∈ feasible,
      (∀ q ∈ feasible, objective p ≤ objective q) → p = selected
  selected_positive : ∀ e, 0 < selected e
  selected_normalized : (∑ e : DirectedSeam, selected e) = 1

/-- Package the source geometry and supplied A2/A3 data into the generic
natural unique projection theorem. -/
def A2A3DirectedSeamProjection.toNaturalUniqueMoveProjection
    (selection : A2A3DirectedSeamProjection) :
    NaturalUniqueMoveProjection
      (Move := DirectedSeam) (Presentation := ProperRotation) where
  feasible := selection.feasible
  objective := selection.objective
  relabel := rotateDirectedSeam
  feasible_natural := selection.feasible_natural
  objective_natural := selection.objective_natural
  selected := selection.selected
  selected_feasible := selection.selected_feasible
  selected_minimal := selection.selected_minimal
  selected_unique := selection.selected_unique
  presentation_transitive := proper_rotations_transitive_on_directed_seams
  selected_positive := selection.selected_positive
  selected_normalized := selection.selected_normalized

/-- A2 covariance, A3 uniqueness, and the exact regular source orbit force
one sixtieth on every directed seam. -/
theorem a2a3_directed_seam_weight_eq_one_sixtieth
    (selection : A2A3DirectedSeamProjection) :
    ∀ e : DirectedSeam, selection.selected e = 1 / 60 := by
  intro e
  rw [show selection.selected e =
      selection.toNaturalUniqueMoveProjection.selected e by rfl]
  rw [selected_eq_inverse_card]
  rw [directedSeam_card]
  norm_num

/-- Translation average on functions over the exact record carrier. -/
noncomputable def weightedEdgeCurrentOperator
    (weight : DirectedSeam → ℝ) (f : D6Point → ℂ) :
    D6Point → ℂ :=
  fun p ↦ ∑ e : DirectedSeam,
    (weight e : ℂ) * f (d6Translate (directedSeamStep e) p)

/-- Source-counting homogeneous edge-current operator. -/
noncomputable def sourceCountingEdgeCurrentOperator
    (f : D6Point → ℂ) : D6Point → ℂ :=
  weightedEdgeCurrentOperator (fun _ ↦ 1 / 60) f

/-- Translation covariance is exact for every weight law. -/
theorem weightedEdgeCurrentOperator_homogeneous
    (weight : DirectedSeam → ℝ) (f : D6Point → ℂ)
    (a p : D6Point) :
    weightedEdgeCurrentOperator weight
        (fun q ↦ f (d6Translate a q)) p =
      weightedEdgeCurrentOperator weight f (d6Translate a p) := by
  unfold weightedEdgeCurrentOperator
  apply Finset.sum_congr rfl
  intro e _
  congr 1
  apply congrArg f
  ext i
  simp [d6Translate, d6Add]
  omega

/-- The A2/A3-selected operator equals the parameter-free source-counting
operator on the exact `D6` record carrier. -/
theorem a2a3_selected_edge_current_operator_eq_source_counting
    (selection : A2A3DirectedSeamProjection) :
    weightedEdgeCurrentOperator selection.selected =
      sourceCountingEdgeCurrentOperator := by
  funext f p
  unfold weightedEdgeCurrentOperator sourceCountingEdgeCurrentOperator
  apply Finset.sum_congr rfl
  intro e _
  rw [a2a3_directed_seam_weight_eq_one_sixtieth selection e]

/-! ## Exact sixty-label Fourier and FZ-12 normalization bridge

This section remains on the source chart.  It proves the algebra behind the
frozen edge-current candidate without identifying that chart or operator with
a physical field.  Its directions are the raw record displacements divided by
two, because the raw displacements have squared norm four while the frozen
moments use unit directions.  The induced `d6Position` chart also contributes
the common `rawRadius⁻¹` factor proved above.  A physical scale convention must
bind both normalizations before laboratory comparison.
-/

/-- Equal source-counting average of all sixty directed seam translations in
the raw response chart. -/
noncomputable def sourceCountingChartAverage
    (a : ℝ) (f : Vec3 → ℂ) : Vec3 → ℂ :=
  fun x ↦ (1 / 60 : ℂ) * ∑ e : DirectedSeam,
    f (x + a • unitDirectedSeamChartStep e)

/-- Reindexing the sixty exact labels gives one positive and one negative
translation for each of the thirty seam-table rows. -/
theorem sourceCountingChartAverage_eq_thirty_pairs
    (a : ℝ) (f : Vec3 → ℂ) (x : Vec3) :
    sourceCountingChartAverage a f x =
      (1 / 60 : ℂ) * ∑ e : Fin 30,
        (f (x + a • unitCarrierSeamDirection e) +
          f (x - a • unitCarrierSeamDirection e)) := by
  unfold sourceCountingChartAverage
  congr 1
  calc
    (∑ e : DirectedSeam, f (x + a • unitDirectedSeamChartStep e)) =
        ∑ e : Fin 30 × Bool,
          f (x + a • unitDirectedSeamChartStep (indexedDirectedSeam e)) := by
      exact (indexedDirectedSeam_bijective.sum_comp
        (fun e : DirectedSeam ↦
          f (x + a • unitDirectedSeamChartStep e))).symm
    _ = ∑ e : Fin 30,
        (f (x + a • unitCarrierSeamDirection e) +
          f (x - a • unitCarrierSeamDirection e)) := by
      rw [Fintype.sum_prod_type]
      apply Finset.sum_congr rfl
      intro e _
      simp [unitDirectedSeamChartStep_indexed, sub_eq_add_neg, add_comm]

/-- The FZ-12 raw-chart generator, with the quadratic normalization fixed by
the sixty-label source average. -/
noncomputable def edgeCurrentGenerator
    (a : ℝ) (f : Vec3 → ℂ) : Vec3 → ℂ :=
  fun x ↦ (1 / (5 * a ^ 2) : ℂ) * ∑ e : Fin 30,
    OPH.PrimitivePortTranslationBridge.pairedDifference
      (a • unitCarrierSeamDirection e) f x

/-- Exact generator identity:
`(6/a²) (I-P) = (1/(5a²)) Σ_e pairedDifference_e`.
No nonzero-scale premise is needed for this algebraic equality. -/
theorem normalized_average_generator_eq_edgeCurrentGenerator
    (a : ℝ) (f : Vec3 → ℂ) (x : Vec3) :
    (6 / a ^ 2 : ℂ) *
        (f x - sourceCountingChartAverage a f x) =
      edgeCurrentGenerator a f x := by
  rw [sourceCountingChartAverage_eq_thirty_pairs]
  unfold edgeCurrentGenerator
  simp only [OPH.PrimitivePortTranslationBridge.pairedDifference]
  have hsumConst : (∑ _e : Fin 30, f x) = 30 * f x := by simp
  have hsumHalf :
      (∑ e : Fin 30,
        (f (x + a • unitCarrierSeamDirection e) +
          f (x - a • unitCarrierSeamDirection e)) / 2) =
        (1 / 2 : ℂ) * ∑ e : Fin 30,
          (f (x + a • unitCarrierSeamDirection e) +
            f (x - a • unitCarrierSeamDirection e)) := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro e _
    ring
  rw [Finset.sum_sub_distrib]
  rw [hsumConst, hsumHalf]
  ring

/-- Exact scalar character sampled by the normalized source generator. -/
noncomputable def edgeCurrentCharacterSymbol (a : ℝ) (k : Vec3) : ℝ :=
  (1 / (5 * a ^ 2)) * ∑ e : Fin 30,
    (1 - Real.cos (a *
      OPH.PrimitivePortTranslationBridge.dot k
        (unitCarrierSeamDirection e)))

/-- The source edge-current generator is diagonal on raw-chart plane waves,
with exactly the FZ-12 character symbol. -/
theorem edgeCurrentGenerator_planeWave
    (a : ℝ) (k x : Vec3) :
    edgeCurrentGenerator a
        (OPH.PrimitivePortTranslationBridge.planeWave k) x =
      (edgeCurrentCharacterSymbol a k : ℂ) *
        OPH.PrimitivePortTranslationBridge.planeWave k x := by
  unfold edgeCurrentGenerator edgeCurrentCharacterSymbol
  simp_rw [OPH.PrimitivePortTranslationBridge.pairedDifference_planeWave,
    OPH.PrimitivePortTranslationBridge.dot_smul_right]
  push_cast
  rw [← Finset.sum_mul]
  ring

/-- A positive/negative translation pair has the cosine character expected
from the exact plane-wave calculation. -/
theorem planeWave_translation_pair
    (a : ℝ) (k x v : Vec3) :
    OPH.PrimitivePortTranslationBridge.planeWave k (x + a • v) +
        OPH.PrimitivePortTranslationBridge.planeWave k (x - a • v) =
      ((2 * Real.cos (a *
        OPH.PrimitivePortTranslationBridge.dot k v) : ℝ) : ℂ) *
        OPH.PrimitivePortTranslationBridge.planeWave k x := by
  have h := OPH.PrimitivePortTranslationBridge.pairedDifference_planeWave
    k x (a • v)
  simp only [OPH.PrimitivePortTranslationBridge.pairedDifference,
    OPH.PrimitivePortTranslationBridge.dot_smul_right] at h
  push_cast at h ⊢
  linear_combination -2 * h

/-- The sixty-label averaging character is exactly the average of the thirty
cosines. -/
theorem sourceCountingChartAverage_planeWave
    (a : ℝ) (k x : Vec3) :
    sourceCountingChartAverage a
        (OPH.PrimitivePortTranslationBridge.planeWave k) x =
      (((1 / 30 : ℝ) * ∑ e : Fin 30,
        Real.cos (a * OPH.PrimitivePortTranslationBridge.dot k
          (unitCarrierSeamDirection e))) : ℂ) *
        OPH.PrimitivePortTranslationBridge.planeWave k x := by
  rw [sourceCountingChartAverage_eq_thirty_pairs]
  simp_rw [planeWave_translation_pair]
  push_cast
  rw [← Finset.sum_mul, ← Finset.mul_sum]
  ring

/-! ## Axiom audit

The finite orbit checks use native evaluation and therefore report the
standard native decision axiom.  The structural action and selection
theorems use ordinary theorem arguments and Mathlib only.  No project axiom or
`sorryAx` occurs.
-/

#print axioms proper_rotations_transitive_on_directed_seams
#print axioms proper_rotations_simply_transitive_on_directed_seams
#print axioms d6Translate_transitive
#print axioms d6Translate_isometry
#print axioms completedTranslation_intertwines
#print axioms directedSeamStep_canonical_chart
#print axioms a2a3_directed_seam_weight_eq_one_sixtieth
#print axioms weightedEdgeCurrentOperator_homogeneous
#print axioms a2a3_selected_edge_current_operator_eq_source_counting
#print axioms sourceCountingChartAverage_eq_thirty_pairs
#print axioms normalized_average_generator_eq_edgeCurrentGenerator
#print axioms edgeCurrentGenerator_planeWave
#print axioms sourceCountingChartAverage_planeWave

end OPH.SeamCurrentHomogeneousAction
