import EventAlgebra.SourcePhaseInstrumentOutcomeBridge

set_option autoImplicit false

/-!
# Common-preparation hull of the generated source-phase relation

This module locates the exact common-support surface of the generated
source-phase relation.  Exactly three generated states admit every one of the
twelve generated phase events.  The declared two-dimensional run state lies
in their convex hull, with exact weights `265/537`, `136/537`, and `136/537`.
Every state in that common support assigns Born weight `1/2` to either
generated Pauli-Y orientation, and its enabled transitions inherit the
normalized declared Lueders-outcome equality proved by
`SourcePhaseInstrumentOutcomeBridge`.

The declared run state is not itself one of the generated states, and no
generated `SourceStep` starts from it.  Moreover, the generated result of
every phase event is outside the common support.  Common support therefore
means pointwise availability of separate one-step cells, not a convexly
closed state relation, reset or re-preparation rule, or sequential protocol.

This is a finite compatibility theorem, not a source-produced common
preparation.  The convex coefficients are an exact decomposition of the
already declared run state; no source rule selects or samples them.  The
twelve generated phase events are not the eight public instrument contexts,
and no instrument implementation, public outcome, readback, run binding,
provenance, custody, cross-context additivity, or physical attachment is
constructed here.  In particular, PR-03, PR-64's public-readback clause,
PR-65, and PR-52 remain open.
-/

namespace EventAlgebra.SourcePhaseCommonPreparationHull

open Matrix
open EventAlgebra
open EventAlgebra.SourcePhaseSelectionInput
open EventAlgebra.SourcePhaseSelectionRelation
open EventAlgebra.SourcePhaseSelectionSemantics
open EventAlgebra.SourcePhaseSelection
open EventAlgebra.SourcePhaseInstrumentOutcomeBridge

noncomputable section

/-- A generated state has common phase support when every generated phase
event is enabled from that same state.  This is a semantic-domain predicate,
not a source-preparation protocol. -/
def HasAllPhaseEvents (state : PhaseState) : Prop :=
  ∀ event : PhaseEvent, HasSourceStepAt state event

/-- The exact common-support set: the three real orbit states, and no
generated Pauli-Y result state, admit all twelve generated phase events. -/
theorem hasAllPhaseEvents_iff (state : PhaseState) :
    HasAllPhaseEvents state ↔
      state = .orbit00 ∨ state = .orbit02 ∨ state = .orbit04 := by
  cases state with
  | orbit00 =>
      constructor
      · intro _
        exact Or.inl rfl
      · intro _ event
        rw [phase_transition_exists_iff_enabled]
        cases event <;> rfl
  | orbit02 =>
      constructor
      · intro _
        exact Or.inr (Or.inl rfl)
      · intro _ event
        rw [phase_transition_exists_iff_enabled]
        cases event <;> rfl
  | orbit04 =>
      constructor
      · intro _
        exact Or.inr (Or.inr rfl)
      · intro _ event
        rw [phase_transition_exists_iff_enabled]
        cases event <;> rfl
  | phaseNegative =>
      constructor
      · intro hall
        have h := hall .pair0002
        rw [phase_transition_exists_iff_enabled] at h
        cases h
      · rintro (h | h | h) <;> cases h
  | phasePositive =>
      constructor
      · intro hall
        have h := hall .pair0004
        rw [phase_transition_exists_iff_enabled] at h
        cases h
      · rintro (h | h | h) <;> cases h

/-- The exact convex weights of the common-support decomposition. -/
def orbit00Weight : ℝ := 265 / 537
def orbit02Weight : ℝ := 136 / 537
def orbit04Weight : ℝ := 136 / 537

theorem commonSupportWeights_nonnegative :
    0 ≤ orbit00Weight ∧ 0 ≤ orbit02Weight ∧ 0 ≤ orbit04Weight := by
  norm_num [orbit00Weight, orbit02Weight, orbit04Weight]

theorem commonSupportWeights_sum_one :
    orbit00Weight + orbit02Weight + orbit04Weight = 1 := by
  norm_num [orbit00Weight, orbit02Weight, orbit04Weight]

/-- The declared run matrix reconstructed as a convex combination of the
three states that admit every generated phase event.  This definition is a
retrospective exact decomposition, not a source mixing operation. -/
def commonSupportMixture : Matrix (Fin 2) (Fin 2) ℂ :=
  orbit00Weight • stateMatrix .orbit00 +
    orbit02Weight • stateMatrix .orbit02 +
    orbit04Weight • stateMatrix .orbit04

/-- The common-support mixture is exactly the declared run state. -/
theorem commonSupportMixture_eq_committedRunState :
    commonSupportMixture = committedRunState := by
  rw [committedRunState_eq_literal]
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [commonSupportMixture, orbit00Weight, orbit02Weight,
      orbit04Weight, stateMatrix, Matrix.smul_apply, Complex.real_smul]

/-- The three coefficients are uniquely fixed once a normalized real
combination of the exact common-support matrices is required to equal the
declared run state.  In fact matrix equality alone forces the coefficients;
normalization is not an extra hypothesis.  This is uniqueness of a
retrospective decomposition, not source selection of the coefficients. -/
theorem commonSupportMixture_coefficients_unique
    (a b c : ℝ) (hmix :
      a • stateMatrix .orbit00 +
          b • stateMatrix .orbit02 +
          c • stateMatrix .orbit04 =
        committedRunState) :
    a = orbit00Weight ∧ b = orbit02Weight ∧ c = orbit04Weight := by
  have h00 := congrFun (congrFun hmix (0 : Fin 2)) (0 : Fin 2)
  have h00re := congrArg Complex.re h00
  have h01 := congrFun (congrFun hmix (0 : Fin 2)) (1 : Fin 2)
  have h01re := congrArg Complex.re h01
  have h11 := congrFun (congrFun hmix (1 : Fin 2)) (1 : Fin 2)
  have h11re := congrArg Complex.re h11
  norm_num [stateMatrix, committedRunState, Matrix.smul_apply,
    Complex.real_smul, Matrix.diagonal_apply] at h00re h01re h11re
  have hsqrt : 0 < OPH.QFT.sqrt3 := OPH.QFT.sqrt3_pos
  constructor
  · norm_num [orbit00Weight]
    nlinarith
  · constructor
    · norm_num [orbit02Weight]
      nlinarith
    · norm_num [orbit04Weight]
      nlinarith

/-- The declared run matrix is not itself any generated `PhaseState` matrix.
The convex-hull identity therefore does not place the hull point in the
generated relation's state type. -/
theorem no_generated_state_eq_committedRunState (state : PhaseState) :
    stateMatrix state ≠ committedRunState := by
  intro heq
  have h00 := congrFun (congrFun heq (0 : Fin 2)) (0 : Fin 2)
  cases state <;>
    norm_num [stateMatrix, committedRunState, Matrix.diagonal_apply] at h00

/-- Consequently no generated source step starts from the declared run
matrix.  A convex-lift or source mixing semantics would be additional data. -/
theorem no_sourceStep_starts_at_committedRunState :
    ¬ ∃ step : SourceStep, stateMatrix step.source = committedRunState := by
  rintro ⟨step, hstep⟩
  exact no_generated_state_eq_committedRunState step.source hstep

/-- After any generated event, the generated result state no longer admits
all twelve phase events.  Common support is therefore a pointwise one-step
property, not a reusable or sequential protocol without an added reset or
re-preparation rule. -/
theorem generatedResultState_not_hasAllPhaseEvents (event : PhaseEvent) :
    ¬ HasAllPhaseEvents (generatedResultState event) := by
  rw [hasAllPhaseEvents_iff]
  cases event <;> simp [generatedResultState, effectValueOf]

/-- Either generated Pauli-Y orientation has Born weight `1/2` on each state
in the exact common support. -/
theorem commonSupport_effectValue_weight_half
    {state : PhaseState} (hstate : HasAllPhaseEvents state)
    (value : EffectValue) :
    bornWeight (stateMatrix state) (effectMatrix value) = 1 / 2 := by
  rw [hasAllPhaseEvents_iff] at hstate
  rcases hstate with rfl | rfl | rfl <;>
    cases value <;>
    norm_num [bornWeight, stateMatrix, effectMatrix, Matrix.trace,
      Matrix.diag, Matrix.mul_apply, Fin.sum_univ_two]

/-- Consequently every generated phase event has weight `1/2` on every
common-support state. -/
theorem commonSupport_generatedEffect_weight_half
    {state : PhaseState} (hstate : HasAllPhaseEvents state)
    (event : PhaseEvent) :
    bornWeight (stateMatrix state) (generatedEffect event) = 1 / 2 := by
  exact commonSupport_effectValue_weight_half hstate (effectValueOf event)

/-- The same exact half weight through the generated event's declared phase
outcome-table entry. -/
theorem commonSupport_declaredOutcome_weight_half
    {state : PhaseState} (hstate : HasAllPhaseEvents state)
    (event : PhaseEvent) :
    bornWeight (stateMatrix state)
        (committedEffectPair InstrumentContext.phase (generatedOutcome event)) =
      1 / 2 := by
  rw [← generatedEffect_eq_committedPhaseOutcome event]
  exact commonSupport_generatedEffect_weight_half hstate event

/-- Every common-support state/event pair inherits the normalized outcome of
the declared Lueders phase instrument.  The conclusion is an equality of
semantic matrices; its table index is not a recorded public outcome. -/
theorem commonSupport_normalized_lueders_outcome_eq_result
    {state : PhaseState} (hstate : HasAllPhaseEvents state)
    (event : PhaseEvent) :
    (bornWeight (stateMatrix state)
        (committedEffectPair InstrumentContext.phase
          (generatedOutcome event)))⁻¹ •
      luedersPhaseInstrument.outcomeMap InstrumentContext.phase
        (generatedOutcome event) (stateMatrix state) =
      stateMatrix (generatedResultState event) := by
  rcases hstate event with ⟨step, hsource, hevent⟩
  have houtcome := sourceStep_normalized_lueders_outcome_eq_result step
  simpa [hsource, hevent, SourceStep.result,
    eventState_eq_generatedResultState] using houtcome

/-- The composed finite compatibility receipt.  It states the exact common
support, convex weights and run-state identity, half weights, and inherited
normalized outcomes, while carrying no operational provenance. -/
theorem sourcePhaseCommonPreparationHull_receipt :
    (∀ state : PhaseState,
      HasAllPhaseEvents state ↔
        state = .orbit00 ∨ state = .orbit02 ∨ state = .orbit04) ∧
    (0 ≤ orbit00Weight ∧ 0 ≤ orbit02Weight ∧ 0 ≤ orbit04Weight) ∧
    orbit00Weight + orbit02Weight + orbit04Weight = 1 ∧
    commonSupportMixture = committedRunState ∧
    (∀ state : PhaseState, stateMatrix state ≠ committedRunState) ∧
    (¬ ∃ step : SourceStep,
      stateMatrix step.source = committedRunState) ∧
    (∀ event : PhaseEvent,
      ¬ HasAllPhaseEvents (generatedResultState event)) ∧
    (∀ a b c : ℝ,
      a • stateMatrix .orbit00 +
            b • stateMatrix .orbit02 +
            c • stateMatrix .orbit04 =
          committedRunState →
        a = orbit00Weight ∧ b = orbit02Weight ∧ c = orbit04Weight) ∧
    (∀ (state : PhaseState), HasAllPhaseEvents state →
      ∀ event : PhaseEvent,
        bornWeight (stateMatrix state) (generatedEffect event) = 1 / 2 ∧
        (bornWeight (stateMatrix state)
            (committedEffectPair InstrumentContext.phase
              (generatedOutcome event)))⁻¹ •
          luedersPhaseInstrument.outcomeMap InstrumentContext.phase
            (generatedOutcome event) (stateMatrix state) =
          stateMatrix (generatedResultState event)) := by
  exact ⟨hasAllPhaseEvents_iff, commonSupportWeights_nonnegative,
    commonSupportWeights_sum_one, commonSupportMixture_eq_committedRunState,
    no_generated_state_eq_committedRunState,
    no_sourceStep_starts_at_committedRunState,
    generatedResultState_not_hasAllPhaseEvents,
    commonSupportMixture_coefficients_unique,
    fun state hstate event =>
      ⟨commonSupport_generatedEffect_weight_half hstate event,
        commonSupport_normalized_lueders_outcome_eq_result hstate event⟩⟩

#print axioms hasAllPhaseEvents_iff
#print axioms commonSupportMixture_eq_committedRunState
#print axioms commonSupportMixture_coefficients_unique
#print axioms no_generated_state_eq_committedRunState
#print axioms no_sourceStep_starts_at_committedRunState
#print axioms generatedResultState_not_hasAllPhaseEvents
#print axioms commonSupport_effectValue_weight_half
#print axioms commonSupport_normalized_lueders_outcome_eq_result
#print axioms sourcePhaseCommonPreparationHull_receipt

end

end EventAlgebra.SourcePhaseCommonPreparationHull
