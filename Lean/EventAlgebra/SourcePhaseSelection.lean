import EventAlgebra.SourcePhaseSelectionRelation
import EventAlgebra.SourcePhaseSelectionSemantics
import EventAlgebra.LuedersPhaseInstrument

set_option autoImplicit false

namespace EventAlgebra.SourcePhaseSelection

open Matrix
open OPH
open EventAlgebra.SourcePhaseSelectionInput
open EventAlgebra.SourcePhaseSelectionSemantics
open EventAlgebra.SourcePhaseSelectionRelation

theorem generated_effect_values_exact :
    generatedEffectValues = [.positive, .negative] := by
  decide

theorem generated_positive_effect_unique :
    ∃! value : EffectValue, upperImaginaryNegative value := by
  refine ⟨.positive, ?_, ?_⟩
  · norm_num [upperImaginaryNegative, effectMatrix]
  · intro value hvalue
    cases value with
    | positive => rfl
    | negative => norm_num [upperImaginaryNegative, effectMatrix] at hvalue

theorem first_noncommuting_effect_eq_positive :
    sourceSelectedGeneratedEffect = effectMatrix .positive := by
  rfl

theorem source_selected_generated_effect_eq_sourcePhaseLift :
    sourceSelectedGeneratedEffect = OPH.QFT.sourcePhaseLift := by
  rw [first_noncommuting_effect_eq_positive,
    OPH.QFT.sourcePhaseLift_eq_rhoYPlus]
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [effectMatrix, OPH.QFT.rhoYPlus] <;> ring

/-- The source-selected effect is exactly the phase-slot effect consumed by
the current declared Lünders instrument.  This proves effect selection only;
it does not construct or source-produce an instrument. -/
theorem source_selected_generated_effect_eq_declared_lueders_effect :
    sourceSelectedGeneratedEffect =
      committedEffectPair InstrumentContext.phase 0 := by
  rw [source_selected_generated_effect_eq_sourcePhaseLift]
  rfl

theorem pair_00_03_effect_eq_sourcePhaseLift :
    generatedEffect .pair0003 = OPH.QFT.sourcePhaseLift := by
  rw [OPH.QFT.sourcePhaseLift_eq_rhoYPlus]
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [generatedEffect, effectValueOf, effectMatrix, OPH.QFT.rhoYPlus] <;>
    ring

theorem source_phase_effect_universe_complete :
    (∀ event : PhaseEvent, event ∈ sourcePairOrder) ∧
      sourcePairOrder.Nodup ∧ sourcePairOrder.length = 12 := by
  decide

theorem phase_transition_exists_iff_enabled :
    ∀ state : PhaseState, ∀ event : PhaseEvent,
      HasSourceStepAt state event ↔ cellStatus state event = .enabled := by
  intro state event
  constructor
  · rintro ⟨step, hstate, hevent⟩
    subst state
    subst event
    simpa [SourceStep.source, SourceStep.event] using step.cell.property
  · intro henabled
    exact ⟨⟨⟨(state, event), henabled⟩⟩, rfl, rfl⟩

theorem disabled_phase_cell_has_no_source_transition :
    DisabledNoTransitionSignature := by
  intro cell hstep
  have henabled :=
    (phase_transition_exists_iff_enabled cell.val.1 cell.val.2).mp hstep
  rw [cell.property] at henabled
  cases henabled

set_option maxHeartbeats 2000000 in
theorem phase_lueders_enabled_update_eq_generated_effect :
    ∀ step : SourceStep,
      EventAlgebra.luedersUpdate
          (stateMatrix step.source) (generatedEffect step.event) =
        generatedEffect step.event := by
  rintro ⟨⟨⟨state, event⟩, enabled⟩⟩
  cases state <;> cases event
  all_goals simp [cellStatus] at enabled
  all_goals
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [SourceStep.source, SourceStep.event,
        EventAlgebra.luedersUpdate, EventAlgebra.bornWeight,
        stateMatrix, generatedEffect, effectValueOf, effectMatrix,
        Matrix.trace, Matrix.diag, Matrix.mul_apply, Fin.sum_univ_two] <;>
      ring_nf <;> simp [Complex.I_sq] <;> norm_num <;> ring

theorem localRepair_eq_generated_result
    (event : PhaseEvent) (state : PhaseState)
    (hneq : state ≠ eventState event) :
    localRepair (Carrier event) (eventSite event) (abstractState event state) =
      abstractState event (eventState event) := by
  have htrigger : LocalTrigger (eventSite event) (abstractState event state) := by
    refine ⟨(), ?_, ?_⟩
    · exact Or.inl rfl
    · simpa [edgeConsistentAt, Carrier, eventSite, abstractState] using hneq
  have hsolvable : LocallySolvable (eventSite event)
      (abstractState event state) := by
    refine ⟨eventState event, ?_⟩
    intro edge _
    cases edge
    simp [edgeConsistentAt, Carrier, eventSite, Function.update]
  have hfire := And.intro htrigger hsolvable
  unfold localRepair
  split
  · rename_i hcond
    have hchosen : Classical.choose hcond.2 = eventState event := by
      have hspec := Classical.choose_spec hcond.2
      have hone := hspec () (Or.inl rfl)
      simpa [edgeConsistentAt, Carrier, eventSite, abstractState,
        Function.update] using hone
    funext patch
    cases patch
    change Classical.choose hcond.2 = eventState event
    exact hchosen
  · rename_i hcond
    exact (hcond hfire).elim

theorem enabled_phase_transition_refines_accepted_or_equality_stutter :
    EnabledTransitionRefinementSignature := by
  intro step
  by_cases h : step.source = step.result
  · exact Or.inr (by
      funext patch
      simpa [abstractState] using h.symm)
  · apply Or.inl
    refine ⟨eventSite step.event, ?_, ?_⟩
    · exact (localRepair_eq_generated_result step.event step.source h).symm
    · intro heq
      apply h
      have hrepair :=
        localRepair_eq_generated_result step.event step.source h
      have hrecords := heq.symm.trans hrepair
      have happ := congrFun hrecords ()
      simpa [abstractState] using happ

structure ExcludedPayloadFields where
  nextStateFingerprint : String
  completenessFingerprint : String
  recordedConclusionFingerprint : String

noncomputable def semanticProjection (_ : ExcludedPayloadFields) :=
  (stateMatrix, generatedEffect, cellStatus)

theorem excluded_payload_fields_do_not_affect_semantic_projection
    (left right : ExcludedPayloadFields) :
    semanticProjection left = semanticProjection right := by
  rfl

theorem anchored_phase_effect_source_selection_candidate :
    generatedEffectValues = [.positive, .negative] ∧
    (∃! value : EffectValue, upperImaginaryNegative value) ∧
    sourceSelectedGeneratedEffect = OPH.QFT.sourcePhaseLift ∧
    sourceSelectedGeneratedEffect =
      committedEffectPair InstrumentContext.phase 0 ∧
    ((∀ event : PhaseEvent, event ∈ sourcePairOrder) ∧
      sourcePairOrder.Nodup ∧ sourcePairOrder.length = 12) ∧
    (∀ state : PhaseState, ∀ event : PhaseEvent,
      HasSourceStepAt state event ↔ cellStatus state event = .enabled) ∧
    DisabledNoTransitionSignature ∧
    EnabledTransitionRefinementSignature ∧
    EventAlgebra.SourcePhaseSelectionSemantics.sourceSelectionPacketSha256 =
      "3412b8fa528635bf670e8d3ba7a1a68558d63524a9795e4603ca5148b6617970" := by
  exact ⟨generated_effect_values_exact, generated_positive_effect_unique,
    source_selected_generated_effect_eq_sourcePhaseLift,
    source_selected_generated_effect_eq_declared_lueders_effect,
    source_phase_effect_universe_complete,
    phase_transition_exists_iff_enabled,
    disabled_phase_cell_has_no_source_transition,
    enabled_phase_transition_refines_accepted_or_equality_stutter, rfl⟩

#print axioms generated_effect_values_exact
#print axioms generated_positive_effect_unique
#print axioms first_noncommuting_effect_eq_positive
#print axioms source_selected_generated_effect_eq_sourcePhaseLift
#print axioms source_selected_generated_effect_eq_declared_lueders_effect
#print axioms pair_00_03_effect_eq_sourcePhaseLift
#print axioms source_phase_effect_universe_complete
#print axioms phase_transition_exists_iff_enabled
#print axioms disabled_phase_cell_has_no_source_transition
#print axioms phase_lueders_enabled_update_eq_generated_effect
#print axioms enabled_phase_transition_refines_accepted_or_equality_stutter
#print axioms excluded_payload_fields_do_not_affect_semantic_projection
#print axioms anchored_phase_effect_source_selection_candidate

end EventAlgebra.SourcePhaseSelection
