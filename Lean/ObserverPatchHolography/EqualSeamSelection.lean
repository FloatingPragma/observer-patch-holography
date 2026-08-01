import Mathlib
import ObserverPatchHolography.CoreAxioms
import ObserverPatchHolography.ScalarSeamRepair

namespace ObserverPatchHolography.EqualSeamSelection

set_option linter.unusedSectionVars false

open OPH.CoreAxioms
open ObserverPatchHolography.ScalarSeamRepair

/-!
# Equal-seam selection: exact conditional theorem and axiom boundary

The canonical A1--A3 structures type a repair interface, interpretation
naturality, and an information projection.  They do not constrain the repair
map itself.  The first section gives an exact shadow-level counterexample:
the same A1 boundary and the same trivial A2/A3 data admit identity and Boolean
flip repairs.

The positive theorem then isolates the smallest finite scalar packet that
does select the equal-seam generator.

* A1-R supplies the complete finite atomic seam alphabet and its endpoints.
* A2-R makes every completed scalar reconciliation local, endpoint-agreeing,
  and endpoint-sum-preserving.  `ScalarSeamRepair` then proves that each move
  is pair averaging.
* A2 naturality and A3 uniqueness make the selected move distribution
  invariant under presentation equivalence.  If that action is transitive,
  normalization forces one over the seam count.

The selected expected repair is therefore

`I - graphLaplacian / (2 * seamCount)`.

This is a theorem with ordinary hypotheses, not a new Lean axiom.  The
canonical basis does not currently provide the A1-R/A2-R fields.  On a
refined mesh, carrier A5 symmetry alone need not be transitive on all refined
edges, so an all-level use additionally needs either a proved single move
orbit or a source-derived unit-counting reference on the complete refined
move grammar.  That reference is an additional branch or source premise until
it is derived from the canonical structures; placing it in an A1/A2/A3 clause
refinement does not itself constitute such a derivation.  No physical clock or
continuum operator is selected here.
-/

noncomputable section

/-! ## Canonical-shadow counterexample -/

/-- The concrete oriented twelve-port boundary as an A1 boundary packet. -/
def concreteBoundary : BoundaryPacket (Fin 12) where
  adjacent := OPH.PortFrameGram.adj
  oriented := orientedFaces
  relabel := Equiv.refl (Fin 12)
  adjacent_iff := by
    intro i j
    rfl
  oriented_matches := by
    exact List.Perm.refl _

/-- An A1 carrier shadow whose Boolean working state is left unchanged. -/
def identityRepairCarrier : CarrierShadow (Fin 12) Bool Unit Unit where
  boundary := concreteBoundary
  readback := fun _ _ => ()
  record := fun _ => ()
  repair := id
  checkpoint := fun _ => ()

/-- The same A1 carrier shadow with Boolean flip as its repair map. -/
def flipRepairCarrier : CarrierShadow (Fin 12) Bool Unit Unit where
  boundary := concreteBoundary
  readback := fun _ _ => ()
  record := fun _ => ()
  repair := Bool.not
  checkpoint := fun _ => ()

/-- A shared trivial A2 meaning-naturality packet. -/
def trivialMeaningNaturality : MeaningNaturality Unit Unit Unit where
  interpret := fun _ _ => ()
  dataTransport := fun _ _ _ => ()
  meaningTransport := fun _ _ _ => ()
  natural := by
    intro _ _ _
    rfl

/-- A shared one-point A3 information-projection packet. -/
def trivialInfoProjection : InfoProjectionSpec Unit Unit Unit where
  feasible := fun _ => True
  restrict := fun _ _ => ()
  reference := fun _ => ()
  cover := {()}
  cover_state_determining := by
    intro s t _ _ _
    exact Subsingleton.elim s t
  weight := fun _ => 1
  weight_positive := by
    intro _ _
    norm_num
  localDivergence := fun _ _ => 0
  divergence := fun _ => 0
  divergence_is_weighted := by
    intro _
    simp
  realized := ()
  realized_feasible := trivial
  realized_minimal := by
    intro _ _
    norm_num
  realized_unique := by
    intro s _ _
    exact Subsingleton.elim s ()

/-- The typed canonical carrier shadow does not identify a repair law.  Both
carriers have the same boundary, readback, records, checkpoint interface, and
shared A2/A3 packets, while their repairs differ. -/
theorem canonical_shadow_repair_not_identified :
    identityRepairCarrier.repair ≠ flipRepairCarrier.repair := by
  intro h
  have hfalse := congrFun h false
  simp [identityRepairCarrier, flipRepairCarrier] at hfalse

/-! ## Conditional complete-grammar theorem -/

variable {ι Move Presentation : Type*}
  [DecidableEq ι] [Fintype Move] [DecidableEq Move]

/-- The finite scalar content of the proposed A1-R/A2-R clauses.  `Move` is
the complete primitive attempt alphabet; each completed move is a linear
endomorphism on the common scalar workspace.  Local support, agreement, and
sum preservation are the scalar specialization of completed reconciliation.
-/
structure CompleteScalarReconciliationGrammar where
  left : Move → ι
  right : Move → ι
  completed : Move → ((ι → ℝ) →ₗ[ℝ] (ι → ℝ))
  supported :
    ∀ e x w, w ≠ left e → w ≠ right e → completed e x w = x w
  agrees : ∀ e x, completed e x (left e) = completed e x (right e)
  preservesEndpointSum :
    ∀ e x,
      completed e x (left e) + completed e x (right e) =
        x (left e) + x (right e)

/-- A presentation-natural A3 move-selection problem.  The feasible family
and objective are invariant under every admitted A2 presentation map; the
A3 minimizer is unique; and the admitted maps act transitively on the move
alphabet.  This is the exact hypothesis bundle that turns presentation
naturality into equal weights on one move orbit. -/
structure NaturalUniqueMoveProjection where
  feasible : Set (Move → ℝ)
  objective : (Move → ℝ) → ℝ
  relabel : Presentation → Move → Move
  feasible_natural :
    ∀ g, ∀ p ∈ feasible, (p ∘ relabel g) ∈ feasible
  objective_natural :
    ∀ g, ∀ p ∈ feasible, objective (p ∘ relabel g) = objective p
  selected : Move → ℝ
  selected_feasible : selected ∈ feasible
  selected_minimal : ∀ p ∈ feasible, objective selected ≤ objective p
  selected_unique :
    ∀ p ∈ feasible,
      (∀ q ∈ feasible, objective p ≤ objective q) → p = selected
  presentation_transitive :
    ∀ e f : Move, ∃ g : Presentation, relabel g e = f
  selected_positive : ∀ e, 0 < selected e
  selected_normalized : (Finset.univ.sum selected) = 1

/-- A source-emitted unit-counting A3 projection on the *complete* primitive
move alphabet.  This packet is the constructive alternative when geometry
has several move orbits.  The reference is fixed before optimization to one
unit per registered primitive event.  The objective is nonnegative on the
feasible set and vanishes there only at that counting reference.

The canonical A1--A3 basis does not currently construct this packet on every
refinement level.  In particular, A5 invariance within each edge orbit does
not supply its cross-orbit reference field. -/
structure ExactCountingMoveProjection where
  feasible : Set (Move → ℝ)
  objective : (Move → ℝ) → ℝ
  selected : Move → ℝ
  selected_feasible : selected ∈ feasible
  selected_minimal : ∀ p ∈ feasible, objective selected ≤ objective p
  counting_feasible :
    (fun _ : Move ↦ 1 / (Fintype.card Move : ℝ)) ∈ feasible
  counting_zero :
    objective (fun _ : Move ↦ 1 / (Fintype.card Move : ℝ)) = 0
  objective_nonnegative : ∀ p ∈ feasible, 0 ≤ objective p
  zero_identifies_counting :
    ∀ p ∈ feasible, objective p = 0 →
      p = (fun _ : Move ↦ 1 / (Fintype.card Move : ℝ))

/-- A unique natural minimizer is invariant under every admitted
presentation map. -/
theorem selected_invariant
    (selection : NaturalUniqueMoveProjection
      (Move := Move) (Presentation := Presentation))
    (g : Presentation) :
    selection.selected ∘ selection.relabel g = selection.selected := by
  have hfeasible :
      (selection.selected ∘ selection.relabel g) ∈ selection.feasible :=
    selection.feasible_natural g selection.selected selection.selected_feasible
  have hminimum :
      ∀ q ∈ selection.feasible,
        selection.objective
            (selection.selected ∘ selection.relabel g) ≤
          selection.objective q := by
    intro q hq
    rw [selection.objective_natural g selection.selected
      selection.selected_feasible]
    exact selection.selected_minimal q hq
  exact selection.selected_unique
    (selection.selected ∘ selection.relabel g) hfeasible hminimum

/-- Transitivity of the admitted presentation maps makes the selected move
weight constant. -/
theorem selected_constant
    (selection : NaturalUniqueMoveProjection
      (Move := Move) (Presentation := Presentation)) :
    ∀ e f : Move, selection.selected e = selection.selected f := by
  intro e f
  obtain ⟨g, hgef⟩ := selection.presentation_transitive e f
  have hinvariant := congrFun (selected_invariant selection g) e
  simpa [Function.comp_apply, hgef] using hinvariant.symm

/-- A constant normalized weight on a nonempty finite move alphabet is the
uniform counting weight. -/
theorem constant_normalized_eq_inverse_card
    [Nonempty Move]
    (weight : Move → ℝ)
    (hconstant : ∀ e f : Move, weight e = weight f)
    (hnormalized : (Finset.univ.sum weight) = 1) :
    ∀ e : Move, weight e = 1 / (Fintype.card Move : ℝ) := by
  intro e
  have hall :
      ∀ f ∈ (Finset.univ : Finset Move), weight f = weight e := by
    intro f _
    exact hconstant f e
  have hsum :
      (Finset.univ.sum weight) = (Fintype.card Move : ℝ) * weight e := by
    rw [Finset.sum_congr rfl hall]
    simp [Finset.sum_const, Finset.card_univ]
  have hcard : (Fintype.card Move : ℝ) ≠ 0 := by
    exact_mod_cast Fintype.card_ne_zero
  rw [hnormalized] at hsum
  field_simp
  nlinarith

/-- The A2/A3 selection packet fixes the uniform counting measure on one
transitive primitive move orbit. -/
theorem selected_eq_inverse_card
    [Nonempty Move]
    (selection : NaturalUniqueMoveProjection
      (Move := Move) (Presentation := Presentation)) :
    ∀ e : Move,
      selection.selected e = 1 / (Fintype.card Move : ℝ) :=
  constant_normalized_eq_inverse_card selection.selected
    (selected_constant selection) selection.selected_normalized

/-- A nonnegative information projection onto a source-emitted unit-counting
reference selects that reference exactly.  Transitivity is not required, so
the theorem applies to a complete refined move alphabet containing several
A5 edge orbits. -/
theorem exactCounting_selected_eq_inverse_card
    (selection : ExactCountingMoveProjection (Move := Move)) :
    selection.selected =
      (fun _ : Move ↦ 1 / (Fintype.card Move : ℝ)) := by
  have hle : selection.objective selection.selected ≤ 0 := by
    have hminimal :=
      selection.selected_minimal
        (fun _ : Move ↦ 1 / (Fintype.card Move : ℝ))
        selection.counting_feasible
    rw [selection.counting_zero] at hminimal
    exact hminimal
  have hge : 0 ≤ selection.objective selection.selected :=
    selection.objective_nonnegative selection.selected
      selection.selected_feasible
  have hzero : selection.objective selection.selected = 0 :=
    le_antisymm hle hge
  exact selection.zero_identifies_counting selection.selected
    selection.selected_feasible hzero

/-- Expected completed reconciliation under the selected move law. -/
def selectedRepair
    (grammar : CompleteScalarReconciliationGrammar (ι := ι) (Move := Move))
    (selection : NaturalUniqueMoveProjection
      (Move := Move) (Presentation := Presentation)) :
    (ι → ℝ) →ₗ[ℝ] (ι → ℝ) :=
  ∑ e : Move, selection.selected e • grammar.completed e

/-- Expected completed reconciliation for a source-emitted exact counting
projection. -/
def exactCountingRepair
    (grammar : CompleteScalarReconciliationGrammar (ι := ι) (Move := Move))
    (selection : ExactCountingMoveProjection (Move := Move)) :
    (ι → ℝ) →ₗ[ℝ] (ι → ℝ) :=
  ∑ e : Move, selection.selected e • grammar.completed e

/-- A2-R's scalar support, equalizer, and conservation clauses force every
completed primitive move to be pair averaging. -/
theorem completed_eq_pairAverage
    (grammar : CompleteScalarReconciliationGrammar (ι := ι) (Move := Move))
    (e : Move) :
    grammar.completed e = pairAverage (grammar.left e) (grammar.right e) :=
  eq_pairAverage_of_supported_agreeing_sum_preserving
    (grammar.left e) (grammar.right e) (grammar.completed e)
    (grammar.supported e) (grammar.agrees e)
    (grammar.preservesEndpointSum e)

/-- Complete atomic grammar, completed local reconciliation, and one-orbit
natural A3 selection give the uniform seam-repair operator. -/
theorem selectedRepair_eq_uniformSeamRepair
    [Nonempty Move]
    (grammar : CompleteScalarReconciliationGrammar (ι := ι) (Move := Move))
    (selection : NaturalUniqueMoveProjection
      (Move := Move) (Presentation := Presentation)) :
    selectedRepair grammar selection =
      uniformSeamRepair grammar.left grammar.right := by
  rw [selectedRepair, uniformSeamRepair]
  simp_rw [completed_eq_pairAverage grammar]
  simp_rw [selected_eq_inverse_card selection]
  rw [Finset.smul_sum]

/-- The composed finite theorem.  No relaxation coefficient remains: one
completed tick is the uniform average of all primitive completed moves, so
its generator is the graph Laplacian divided by twice the move count. -/
theorem selectedRepair_eq_id_sub_graphLaplacian
    [Nonempty Move]
    (grammar : CompleteScalarReconciliationGrammar (ι := ι) (Move := Move))
    (selection : NaturalUniqueMoveProjection
      (Move := Move) (Presentation := Presentation)) :
    selectedRepair grammar selection =
      LinearMap.id -
        (1 / (2 * (Fintype.card Move : ℝ))) •
          graphLaplacian grammar.left grammar.right := by
  rw [selectedRepair_eq_uniformSeamRepair grammar selection]
  exact uniformSeamRepair_eq_id_sub_graphLaplacian
    grammar.left grammar.right

/-- Complete atomic grammar plus a source-emitted unit-counting A3 reference
gives the same exact finite Laplacian generator without assuming that A5 acts
transitively on the refined move alphabet. -/
theorem exactCountingRepair_eq_id_sub_graphLaplacian
    [Nonempty Move]
    (grammar : CompleteScalarReconciliationGrammar (ι := ι) (Move := Move))
    (selection : ExactCountingMoveProjection (Move := Move)) :
    exactCountingRepair grammar selection =
      LinearMap.id -
        (1 / (2 * (Fintype.card Move : ℝ))) •
          graphLaplacian grammar.left grammar.right := by
  have hrepair :
      exactCountingRepair grammar selection =
        uniformSeamRepair grammar.left grammar.right := by
    rw [exactCountingRepair, uniformSeamRepair]
    simp_rw [completed_eq_pairAverage grammar]
    rw [exactCounting_selected_eq_inverse_card selection]
    rw [Finset.smul_sum]
  rw [hrepair]
  exact uniformSeamRepair_eq_id_sub_graphLaplacian
    grammar.left grammar.right

/-- The base-carrier specialization: a complete transitive thirty-seam
grammar gives `I - L/60`. -/
theorem selected_thirty_seam_repair
    (grammar : CompleteScalarReconciliationGrammar
      (ι := ι) (Move := Fin 30))
    (selection : NaturalUniqueMoveProjection
      (Move := Fin 30) (Presentation := Presentation)) :
    selectedRepair grammar selection =
      LinearMap.id - (1 / 60 : ℝ) •
        graphLaplacian grammar.left grammar.right := by
  rw [selectedRepair_eq_id_sub_graphLaplacian grammar selection]
  norm_num

end

end ObserverPatchHolography.EqualSeamSelection

/-! Axiom audit: all results use ordinary theorem arguments and Mathlib. -/

#print axioms ObserverPatchHolography.EqualSeamSelection.canonical_shadow_repair_not_identified
#print axioms ObserverPatchHolography.EqualSeamSelection.selectedRepair_eq_id_sub_graphLaplacian
#print axioms ObserverPatchHolography.EqualSeamSelection.exactCountingRepair_eq_id_sub_graphLaplacian
#print axioms ObserverPatchHolography.EqualSeamSelection.selected_thirty_seam_repair
