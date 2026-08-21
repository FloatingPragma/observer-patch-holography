import EventAlgebra.SourcePhaseSelection
import EventAlgebra.LuedersPhaseInstrument

set_option autoImplicit false

/-!
# Source-phase selection to the declared Lüders outcome map

This module joins the generated source-phase semantics to the already declared
Lüders instrument, at the matrix level only.  A generated positive phase
effect is assigned the declared phase-outcome label `0`, and a generated
negative phase effect is assigned label `1`.  For every enabled `SourceStep`,
the corresponding normalized declared Lüders outcome map is then proved to
equal the generated semantic result matrix.

The label below is an index into `committedEffectPair`; it is not a public or
recorded outcome.  Nothing here constructs a run, common preparation,
readback, provenance, custody, or physical implementation.  In particular,
the bridge does not discharge PR-03, PR-64's public-readback clause, PR-65, or
the source-selection status of the declared Lüders instrument.
-/

namespace EventAlgebra.SourcePhaseInstrumentOutcomeBridge

open EventAlgebra
open EventAlgebra.SourcePhaseSelectionInput
open EventAlgebra.SourcePhaseSelectionRelation
open EventAlgebra.SourcePhaseSelectionSemantics

noncomputable section

/-- The generated effect's index in the declared binary phase-effect pair.
This is a semantic table index, not a recorded public outcome. -/
def generatedOutcome : PhaseEvent → Fin 2
  | .pair0002 => 0
  | .pair0003 => 0
  | .pair0004 => 1
  | .pair0005 => 1
  | .pair0102 => 0
  | .pair0103 => 0
  | .pair0104 => 1
  | .pair0105 => 1
  | .pair0204 => 0
  | .pair0205 => 0
  | .pair0304 => 0
  | .pair0305 => 0

/-- The phase state selected by the generated event semantics. -/
def generatedResultState (event : PhaseEvent) : PhaseState :=
  match effectValueOf event with
  | .positive => .phasePositive
  | .negative => .phaseNegative

/-- The generated result state is the unique state whose committed bytes
match the event operation bytes. -/
theorem eventState_eq_generatedResultState (event : PhaseEvent) :
    eventState event = generatedResultState event := by
  have hmem : generatedResultState event ∈ matchingStates event := by
    cases event <;> decide
  exact ((Classical.choose_spec (matchingStates_exists_unique event)).2
    (generatedResultState event) hmem).symm

/-- The matrix of the generated result state is exactly the generated phase
effect. -/
theorem stateMatrix_generatedResultState (event : PhaseEvent) :
    stateMatrix (generatedResultState event) = generatedEffect event := by
  cases event <;> rfl

/-- Every matrix in the generated phase-state table carries an existing
state certificate.  The three real states reuse the committed web-projector
event certificates; the two phase states reuse the Pauli-Y certificates. -/
theorem stateMatrix_isState (state : PhaseState) :
    IsState (stateMatrix state) := by
  cases state with
  | orbit00 =>
      constructor
      · have heq : stateMatrix .orbit00 =
            OPH.QFT.complexifyRealMatrix
              (webContextProjector OPH.QFT.WebContext.diagonal) := by
          ext i j
          fin_cases i <;> fin_cases j <;>
            norm_num [stateMatrix, webContextProjector,
              OPH.QFT.complexifyRealMatrix, OPH.QFT.recordProjector]
        rw [heq]
        exact (webContextProjector_isEvent
          OPH.QFT.WebContext.diagonal).posSemidef
      · norm_num [stateMatrix, Matrix.trace, Matrix.diag,
          Fin.sum_univ_two]
  | orbit02 =>
      constructor
      · have heq : stateMatrix .orbit02 =
            OPH.QFT.complexifyRealMatrix (webContextProjector
              (OPH.QFT.WebContext.conjugated (2 : Fin 6))) := by
          rw [webContextProjector_conjugated, OPH.QFT.conjProjector_two]
          ext i j
          fin_cases i <;> fin_cases j <;>
            norm_num [stateMatrix, OPH.QFT.complexifyRealMatrix] <;> ring
        rw [heq]
        exact (webContextProjector_isEvent
          (OPH.QFT.WebContext.conjugated (2 : Fin 6))).posSemidef
      · norm_num [stateMatrix, Matrix.trace, Matrix.diag,
          Fin.sum_univ_two]
  | orbit04 =>
      constructor
      · have heq : stateMatrix .orbit04 =
            OPH.QFT.complexifyRealMatrix (webContextProjector
              (OPH.QFT.WebContext.conjugated (4 : Fin 6))) := by
          rw [webContextProjector_conjugated, OPH.QFT.conjProjector_four]
          ext i j
          fin_cases i <;> fin_cases j <;>
            norm_num [stateMatrix, OPH.QFT.complexifyRealMatrix] <;> ring
        rw [heq]
        exact (webContextProjector_isEvent
          (OPH.QFT.WebContext.conjugated (4 : Fin 6))).posSemidef
      · norm_num [stateMatrix, Matrix.trace, Matrix.diag,
          Fin.sum_univ_two]
  | phaseNegative =>
      have heq : stateMatrix .phaseNegative = OPH.QFT.rhoYMinus := by
        ext i j
        fin_cases i <;> fin_cases j <;>
          norm_num [stateMatrix, OPH.QFT.rhoYMinus] <;> ring
      rw [heq]
      exact OPH.QFT.rhoYMinus_isState
  | phasePositive =>
      have heq : stateMatrix .phasePositive = OPH.QFT.rhoYPlus := by
        ext i j
        fin_cases i <;> fin_cases j <;>
          norm_num [stateMatrix, OPH.QFT.rhoYPlus] <;> ring
      rw [heq]
      exact OPH.QFT.rhoYPlus_isState

/-- The generated effect is exactly the declared phase-effect entry selected
by `generatedOutcome`.  This identifies effects and outcome-map indices only;
it does not assert that an outcome was observed. -/
theorem generatedEffect_eq_committedPhaseOutcome (event : PhaseEvent) :
    generatedEffect event =
      committedEffectPair InstrumentContext.phase (generatedOutcome event) := by
  cases event <;>
    ext i j <;>
    fin_cases i <;>
    fin_cases j <;>
    norm_num [generatedEffect, generatedOutcome, effectValueOf, effectMatrix,
      committedEffectPair, committedContextEffect,
      OPH.QFT.sourcePhaseLift_eq_rhoYPlus, OPH.QFT.rhoYPlus,
      Matrix.one_apply] <;>
    ring

/-- Every generated phase effect is nonzero. -/
theorem generatedEffect_ne_zero (event : PhaseEvent) :
    generatedEffect event ≠ 0 := by
  intro hzero
  have h00 := congrFun (congrFun hzero (0 : Fin 2)) (0 : Fin 2)
  cases event <;>
    norm_num [generatedEffect, effectValueOf, effectMatrix] at h00

/-- Enabled source steps have nonzero weight for their generated effect.  This
rules out the totalized zero-weight branch of `luedersUpdate` in the bridge. -/
theorem sourceStep_generated_bornWeight_ne_zero (step : SourceStep) :
    bornWeight (stateMatrix step.source) (generatedEffect step.event) ≠ 0 := by
  intro hzero
  have hupdate :=
    EventAlgebra.SourcePhaseSelection.phase_lueders_enabled_update_eq_generated_effect
      step
  have hz : (0 : Matrix (Fin 2) (Fin 2) ℂ) = generatedEffect step.event := by
    simpa [EventAlgebra.luedersUpdate, hzero] using hupdate
  exact generatedEffect_ne_zero step.event hz.symm

/-- For every enabled source step, the normalized outcome map of the declared
Lüders phase instrument at the generated effect's declared table index equals
the generated semantic result matrix.  This is a conditional semantic bridge,
not evidence that the declared instrument was source-produced or run. -/
theorem sourceStep_normalized_lueders_outcome_eq_result (step : SourceStep) :
    (bornWeight (stateMatrix step.source)
        (committedEffectPair InstrumentContext.phase
          (generatedOutcome step.event)))⁻¹ •
      luedersPhaseInstrument.outcomeMap InstrumentContext.phase
        (generatedOutcome step.event) (stateMatrix step.source) =
      stateMatrix step.result := by
  have hevent : IsEvent (generatedEffect step.event) := by
    rw [generatedEffect_eq_committedPhaseOutcome]
    exact committedEffectPair_isEvent InstrumentContext.phase
      (generatedOutcome step.event)
  rw [luedersPhaseInstrument_outcomeMap,
    ← generatedEffect_eq_committedPhaseOutcome step.event,
    luedersOutcomeMap_normalized hevent,
    EventAlgebra.SourcePhaseSelection.phase_lueders_enabled_update_eq_generated_effect]
  rw [SourceStep.result, eventState_eq_generatedResultState,
    stateMatrix_generatedResultState]

/-- The composed stated-domain receipt: all generated matrices are states,
the generated effects occupy the declared phase table entries, the enabled
source domain has 48 cells, and every enabled step has nonzero weight and
equals the normalized output of the corresponding declared Lüders map.  The
receipt carries semantic table indices only, not public outcomes or run
custody. -/
theorem sourcePhaseInstrumentOutcomeBridge_receipt :
    (∀ state : PhaseState, IsState (stateMatrix state)) ∧
    (∀ event : PhaseEvent,
      generatedEffect event =
        committedEffectPair InstrumentContext.phase (generatedOutcome event)) ∧
    enabledCells.length = 48 ∧
    (∀ step : SourceStep,
      bornWeight (stateMatrix step.source) (generatedEffect step.event) ≠ 0 ∧
      (bornWeight (stateMatrix step.source)
          (committedEffectPair InstrumentContext.phase
            (generatedOutcome step.event)))⁻¹ •
        luedersPhaseInstrument.outcomeMap InstrumentContext.phase
          (generatedOutcome step.event) (stateMatrix step.source) =
        stateMatrix step.result) := by
  exact ⟨stateMatrix_isState, generatedEffect_eq_committedPhaseOutcome,
    enabledCellCensus, fun step =>
      ⟨sourceStep_generated_bornWeight_ne_zero step,
        sourceStep_normalized_lueders_outcome_eq_result step⟩⟩

#print axioms eventState_eq_generatedResultState
#print axioms stateMatrix_isState
#print axioms generatedEffect_eq_committedPhaseOutcome
#print axioms sourceStep_generated_bornWeight_ne_zero
#print axioms sourceStep_normalized_lueders_outcome_eq_result
#print axioms sourcePhaseInstrumentOutcomeBridge_receipt

end

end EventAlgebra.SourcePhaseInstrumentOutcomeBridge
