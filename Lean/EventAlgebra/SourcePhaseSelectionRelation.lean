import EventAlgebra.Lueders
import ObserverPatchHolography.Primitives
import EventAlgebra.SourcePhaseSelectionInput

set_option autoImplicit false

/-!
# Source-bound EventAlgebra relation for phase selection

The relation uses a
decidable search over the generated state and event arrays.
-/

namespace EventAlgebra.SourcePhaseSelectionRelation

open OPH
open EventAlgebra.SourcePhaseSelectionInput

def matchingStates (event : PhaseEvent) : Finset PhaseState :=
  Finset.univ.filter fun state => stateBytes state = operationBytes event

theorem matchingStates_card_one :
    ∀ event : PhaseEvent, (matchingStates event).card = 1 := by
  intro event
  cases event <;> decide

theorem matchingStates_eq_singleton (event : PhaseEvent) :
    ∃ state : PhaseState, matchingStates event = {state} :=
  Finset.card_eq_one.mp (matchingStates_card_one event)

theorem matchingStates_exists_unique (event : PhaseEvent) :
    ∃! state : PhaseState, state ∈ matchingStates event := by
  rcases matchingStates_eq_singleton event with ⟨state, hstate⟩
  refine ⟨state, by simp [hstate], ?_⟩
  intro other hother
  simpa [hstate] using hother

noncomputable def eventState (event : PhaseEvent) : PhaseState :=
  Classical.choose (matchingStates_exists_unique event)

theorem eventState_mem_matchingStates (event : PhaseEvent) :
    eventState event ∈ matchingStates event :=
  (Classical.choose_spec (matchingStates_exists_unique event)).1

theorem eventStateOperationBytes (event : PhaseEvent) :
    stateBytes (eventState event) = operationBytes event := by
  simpa [matchingStates] using eventState_mem_matchingStates event

noncomputable def Carrier (event : PhaseEvent) : OPHCarrier where
  Patch := Unit
  State := fun _ => PhaseState
  Edge := Unit
  src := fun _ => ()
  tgt := fun _ => ()
  Iface := fun _ => PhaseState
  projSrc := fun _ state => state
  projTgt := fun _ _ => eventState event
  weight := fun _ => 1
  dist := fun _ left right => if left = right then 0 else 1
  weight_pos := by intro _; norm_num
  dist_eq_zero := by intro _ left right; simp

noncomputable def abstractState (event : PhaseEvent) (state : PhaseState) :
    Records (Carrier event) :=
  fun _ => state

noncomputable def eventSite (event : PhaseEvent) : Site (Carrier event) := ()

structure SourceStep where
  cell : EnabledCell
  deriving DecidableEq, Fintype

def SourceStep.source (step : SourceStep) : PhaseState := step.cell.val.1
def SourceStep.event (step : SourceStep) : PhaseEvent := step.cell.val.2
noncomputable def SourceStep.result (step : SourceStep) : PhaseState :=
  eventState step.event

def HasSourceStepAt (state : PhaseState) (event : PhaseEvent) : Prop :=
  ∃ step : SourceStep, step.source = state ∧ step.event = event

def DisabledNoTransitionSignature : Prop :=
  ∀ cell : DisabledCell, ¬ HasSourceStepAt cell.val.1 cell.val.2

def EnabledTransitionRefinementSignature : Prop :=
  ∀ step : SourceStep,
    acceptedStep (Carrier step.event)
      (abstractState step.event step.source)
      (abstractState step.event step.result) ∨
    abstractState step.event step.result =
      abstractState step.event step.source

def changingEnabledCells : List EnabledCell :=
  enabledCells.filter fun cell =>
    decide (stateBytes cell.val.1 ≠ operationBytes cell.val.2)

def alreadyOperationEnabledCells : List EnabledCell :=
  enabledCells.filter fun cell =>
    decide (stateBytes cell.val.1 = operationBytes cell.val.2)

theorem changingEnabledCellCensus : changingEnabledCells.length = 36 := by decide
theorem alreadyOperationEnabledCellCensus :
    alreadyOperationEnabledCells.length = 12 := by decide

def sourceBoundStep : SourceStep where
  cell := sourceBoundEnabledCell

def sourceBoundDisabledWitness : DisabledCell := sourceBoundDisabledCell

noncomputable def sourceBoundCarrier : OPHCarrier := Carrier sourceBoundStep.event

noncomputable def sourceBoundRecord : Records sourceBoundCarrier :=
  abstractState sourceBoundStep.event sourceBoundStep.source

#print axioms matchingStates_exists_unique
#print axioms eventStateOperationBytes
#print axioms changingEnabledCellCensus
#print axioms alreadyOperationEnabledCellCensus

end EventAlgebra.SourcePhaseSelectionRelation
