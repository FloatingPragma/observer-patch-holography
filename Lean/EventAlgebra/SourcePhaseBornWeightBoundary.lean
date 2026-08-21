import EventAlgebra.SourceReachabilityDelimitation
import EventAlgebra.SourcePhaseSelection
import EventAlgebra.PhaseInstrumentDetermination

set_option autoImplicit false

/-!
# Source-phase Born-weight boundary

This module joins the current source-reachable class to the source-selected
two-dimensional phase effect without inventing a map between their different
carriers.  A prospective bridge is typed explicitly.  If it sends every
record-diagonal source matrix to a record-diagonal two-dimensional
preparation, then every reachable preparation assigns Born weight `1/2` to
both source-generated Pauli-Y orientations.  The same weight is assigned to
the diagonal comparator `(1/2) I`.

More exactly, the selected-effect weight of any two-dimensional state is
`1/2 - Im(rho 0 1)`.  Thus a non-half selected-effect weight requires a
nonzero imaginary off-diagonal entry.  A prospective source-to-phase bridge
whose output differs from `1/2` therefore cannot preserve record diagonality
on that input.  This is a conditional, single-shot Born-weight statement,
not a construction of PR-65 and not an operational readback no-go.  The
module types no source operation, public outcome, sequential instrument
behavior, run, provenance, custody, physical region, or validation claim.
-/

namespace EventAlgebra.SourcePhaseBornWeightBoundary

open Matrix
open OPH.QFT
open EventAlgebra.SourceReachability
open EventAlgebra.SourcePhaseSelection
open EventAlgebra.SourcePhaseSelectionSemantics
open scoped ComplexOrder

noncomputable section

abbrev PhaseMatrix := Matrix (Fin 2) (Fin 2) ℂ

/-! ## Exact selected-effect weight on two-dimensional preparations -/

/-- The lower off-diagonal entry of a two-dimensional state is the complex
conjugate of the upper entry. -/
theorem state_offdiag_conj {rho : PhaseMatrix} (hrho : IsState rho) :
    rho 1 0 = (starRingEnd ℂ) (rho 0 1) := by
  have h := congrFun (congrFun hrho.1.isHermitian.eq 1) 0
  rw [Matrix.conjTranspose_apply] at h
  simpa using h.symm

/-- The source-selected Pauli-`+Y` effect in commutator coordinates. -/
theorem state_sourceSelected_weight_commutator_form
    {rho : PhaseMatrix} (hrho : IsState rho) :
    bornWeight rho sourceSelectedGeneratedEffect =
      1 / 2 + Complex.I / 2 * (rho 0 1 - rho 1 0) := by
  have htrace := hrho.2
  rw [source_selected_generated_effect_eq_sourcePhaseLift,
    sourcePhaseLift_eq_rhoYPlus, bornWeight_fin_two]
  rw [Matrix.trace_fin_two] at htrace
  simp [rhoYPlus]
  linear_combination (1 / 2 : ℂ) * htrace

/-- The source-selected Pauli-`+Y` effect reads exactly the imaginary part of
the upper off-diagonal coordinate of any two-dimensional state. -/
theorem state_sourceSelected_weight_eq_half_sub_im
    {rho : PhaseMatrix} (hrho : IsState rho) :
    bornWeight rho sourceSelectedGeneratedEffect =
      1 / 2 - ((rho 0 1).im : ℂ) := by
  rw [state_sourceSelected_weight_commutator_form hrho,
    state_offdiag_conj hrho, Complex.sub_conj]
  push_cast
  linear_combination ((rho 0 1).im : ℂ) * Complex.I_sq

/-- The diagonal comparator with the same diagonal entries as either
source-generated Pauli-Y projector.  It is an effect, not a projection and
not an instrument. -/
def diagonalComparator : PhaseMatrix := (1 / 2 : ℝ) • (1 : PhaseMatrix)

theorem diagonalComparator_isEffect : IsEffect diagonalComparator := by
  exact isEffect_one.smul (by norm_num) (by norm_num)

/-- Every normalized record-diagonal two-dimensional preparation gives the
source-selected phase effect weight `1/2`. -/
theorem recordDiagonal_state_sourceSelected_weight_half
    {rho : PhaseMatrix} (hrho : IsState rho) (hdiag : IsRecordDiagonal rho) :
    bornWeight rho sourceSelectedGeneratedEffect = 1 / 2 := by
  have h01 : rho 0 1 = 0 := hdiag 0 1 (by decide)
  have h10 : rho 1 0 = 0 := hdiag 1 0 (by decide)
  have htrace := hrho.2
  simp [Matrix.trace, Matrix.diag, Fin.sum_univ_two] at htrace
  rw [source_selected_generated_effect_eq_sourcePhaseLift,
    sourcePhaseLift_eq_rhoYPlus]
  simp [bornWeight, rhoYPlus, Matrix.trace, Matrix.diag, Matrix.mul_apply,
    Fin.sum_univ_two, h01, h10]
  linear_combination (1 / 2 : ℂ) * htrace

/-- The opposite source-generated orientation has the same weight `1/2` on
every normalized record-diagonal preparation. -/
theorem recordDiagonal_state_negative_weight_half
    {rho : PhaseMatrix} (hrho : IsState rho) (hdiag : IsRecordDiagonal rho) :
    bornWeight rho (effectMatrix .negative) = 1 / 2 := by
  have h01 : rho 0 1 = 0 := hdiag 0 1 (by decide)
  have h10 : rho 1 0 = 0 := hdiag 1 0 (by decide)
  have htrace := hrho.2
  simp [Matrix.trace, Matrix.diag, Fin.sum_univ_two] at htrace
  simp [bornWeight, effectMatrix, Matrix.trace, Matrix.diag, Matrix.mul_apply,
    Fin.sum_univ_two, h01, h10]
  linear_combination (1 / 2 : ℂ) * htrace

/-- The diagonal comparator also has weight `1/2` on every state. -/
theorem state_diagonalComparator_weight_half
    {rho : PhaseMatrix} (hrho : IsState rho) :
    bornWeight rho diagonalComparator = 1 / 2 := by
  rw [show diagonalComparator = (1 / 2 : ℂ) • (1 : PhaseMatrix) by
    ext i j
    simp [diagonalComparator, Matrix.smul_apply, Complex.real_smul]]
  simp [bornWeight, Matrix.trace_smul, hrho.2]

/-- On a record-diagonal preparation the selected orientation, the opposite
orientation, and the diagonal comparator have equal single-shot Born
weights. -/
theorem recordDiagonal_state_phase_weights_equal
    {rho : PhaseMatrix} (hrho : IsState rho) (hdiag : IsRecordDiagonal rho) :
    bornWeight rho sourceSelectedGeneratedEffect =
        bornWeight rho (effectMatrix .negative) ∧
      bornWeight rho sourceSelectedGeneratedEffect =
        bornWeight rho diagonalComparator := by
  rw [recordDiagonal_state_sourceSelected_weight_half hrho hdiag,
    recordDiagonal_state_negative_weight_half hrho hdiag,
    state_diagonalComparator_weight_half hrho]
  exact ⟨rfl, rfl⟩

/-! ## Typed bridges from the committed source class -/

/-- A prospective preparation bridge from one committed source carrier to
the two-dimensional phase system.  The state field says only that reachable
source matrices are sent to normalized positive preparations.  It does not
assume a particular transport, source operation, public outcome, or receipt. -/
structure SourceToPhasePreparationBridge (c : Carrier) where
  prepare : Matrix c.Index c.Index ℂ → PhaseMatrix
  state_of_reachable : ∀ {M : Matrix c.Index c.Index ℂ},
    Reachable c M → IsState (prepare M)

/-- The bridge preserves record diagonality when every record-diagonal input
is mapped to a record-diagonal phase preparation. -/
def PreservesRecordDiagonality {c : Carrier}
    (B : SourceToPhasePreparationBridge c) : Prop :=
  ∀ M : Matrix c.Index c.Index ℂ,
    IsRecordDiagonal M → IsRecordDiagonal (B.prepare M)

/-- Every record-preserving bridge from a reachable committed source state
gives the source-selected phase effect weight `1/2`. -/
theorem reachable_recordPreserving_phase_weight_half
    {c : Carrier} (B : SourceToPhasePreparationBridge c)
    (hpres : PreservesRecordDiagonality B)
    {M : Matrix c.Index c.Index ℂ} (hM : Reachable c M) :
    bornWeight (B.prepare M) sourceSelectedGeneratedEffect = 1 / 2 := by
  exact recordDiagonal_state_sourceSelected_weight_half
    (B.state_of_reachable hM) (hpres M (reachable_recordDiagonal hM))

/-- A non-half selected-effect weight on a reachable source state forces the
bridge output to have a nonzero imaginary off-diagonal coordinate. -/
theorem nonhalf_weight_requires_nonzero_im
    {c : Carrier} (B : SourceToPhasePreparationBridge c)
    {M : Matrix c.Index c.Index ℂ} (hM : Reachable c M)
    (hweight : bornWeight (B.prepare M) sourceSelectedGeneratedEffect ≠ 1 / 2) :
    (B.prepare M 0 1).im ≠ 0 := by
  intro him
  apply hweight
  rw [state_sourceSelected_weight_eq_half_sub_im (B.state_of_reachable hM),
    him]
  simp

/-- Consequently, a bridge whose output has a non-half selected-effect weight
on a reachable input does not preserve record diagonality. -/
theorem nonhalf_weight_refutes_recordPreservation
    {c : Carrier} (B : SourceToPhasePreparationBridge c)
    {M : Matrix c.Index c.Index ℂ} (hM : Reachable c M)
    (hweight : bornWeight (B.prepare M) sourceSelectedGeneratedEffect ≠ 1 / 2) :
    ¬ PreservesRecordDiagonality B := by
  intro hpres
  exact hweight (reachable_recordPreserving_phase_weight_half B hpres hM)

/-! ## Load-bearing controls -/

/-- A constant diagonal control bridge on the committed 86/247 carrier.  It
is included only to show that the bridge type and record-preservation
hypothesis are jointly inhabitable; it is not source production. -/
def constantDiagonalControlBridge :
    SourceToPhasePreparationBridge Carrier.pair247 where
  prepare _ := committedRunState
  state_of_reachable _ := committedRunState_isState

theorem constantDiagonalControlBridge_preserves :
    PreservesRecordDiagonality constantDiagonalControlBridge := by
  intro M hM i j hij
  fin_cases i <;> fin_cases j <;>
    first
      | exact (hij rfl).elim
      | simp [constantDiagonalControlBridge, committedRunState]

theorem constantDiagonalControlBridge_weight_half
    {M : Matrix Carrier.pair247.Index Carrier.pair247.Index ℂ}
    (hM : Reachable Carrier.pair247 M) :
    bornWeight (constantDiagonalControlBridge.prepare M)
      sourceSelectedGeneratedEffect = 1 / 2 :=
  reachable_recordPreserving_phase_weight_half constantDiagonalControlBridge
    constantDiagonalControlBridge_preserves hM

/-- A deliberately off-diagonal control bridge.  It ignores its input and
returns the selected Pauli-Y state, so it is not a source construction; its
purpose is to prove that dropping record preservation permits a non-half
selected-effect weight. -/
def offDiagonalControlBridge :
    SourceToPhasePreparationBridge Carrier.pair247 where
  prepare _ := sourceSelectedGeneratedEffect
  state_of_reachable _ := by
    rw [source_selected_generated_effect_eq_sourcePhaseLift,
      sourcePhaseLift_eq_rhoYPlus]
    exact rhoYPlus_isState

theorem sourceSelectedGeneratedEffect_isEvent :
    IsEvent sourceSelectedGeneratedEffect := by
  rw [source_selected_generated_effect_eq_sourcePhaseLift]
  exact sourcePhaseLift_isEvent

theorem sourceSelectedGeneratedEffect_isState :
    IsState sourceSelectedGeneratedEffect := by
  rw [source_selected_generated_effect_eq_sourcePhaseLift,
    sourcePhaseLift_eq_rhoYPlus]
  exact rhoYPlus_isState

theorem offDiagonalControlBridge_weight_one
    (M : Matrix Carrier.pair247.Index Carrier.pair247.Index ℂ) :
    bornWeight (offDiagonalControlBridge.prepare M)
      sourceSelectedGeneratedEffect = 1 := by
  rw [show offDiagonalControlBridge.prepare M =
      sourceSelectedGeneratedEffect by rfl,
    bornWeight, sourceSelectedGeneratedEffect_isEvent.2,
    sourceSelectedGeneratedEffect_isState.2]

theorem offDiagonalControlBridge_not_preserving :
    ¬ PreservesRecordDiagonality offDiagonalControlBridge := by
  apply nonhalf_weight_refutes_recordPreservation
    offDiagonalControlBridge correlationState_reachable
  rw [offDiagonalControlBridge_weight_one]
  norm_num

/-- The composed stated-domain receipt: all current reachable matrices obey
the record law; every record-preserving phase bridge assigns equal half
weight to both source orientations and the diagonal comparator; and a
non-half selected-effect weight requires a nonzero imaginary off-diagonal
entry. -/
theorem sourcePhaseBornWeightBoundary_receipt :
    IsEffect diagonalComparator ∧
    (∀ {c : Carrier} {M : Matrix c.Index c.Index ℂ}, Reachable c M →
      IsRecordLaw M) ∧
    (∀ {rho : PhaseMatrix}, IsState rho → IsRecordDiagonal rho →
      bornWeight rho sourceSelectedGeneratedEffect =
          bornWeight rho (effectMatrix .negative) ∧
        bornWeight rho sourceSelectedGeneratedEffect =
        bornWeight rho diagonalComparator) ∧
    (∀ {c : Carrier} (B : SourceToPhasePreparationBridge c),
      PreservesRecordDiagonality B →
      ∀ {M : Matrix c.Index c.Index ℂ}, Reachable c M →
        bornWeight (B.prepare M) sourceSelectedGeneratedEffect = 1 / 2) ∧
    (∀ {c : Carrier} (B : SourceToPhasePreparationBridge c)
      {M : Matrix c.Index c.Index ℂ}, Reachable c M →
      bornWeight (B.prepare M) sourceSelectedGeneratedEffect ≠ 1 / 2 →
        (B.prepare M 0 1).im ≠ 0) ∧
    PreservesRecordDiagonality constantDiagonalControlBridge ∧
    ¬ PreservesRecordDiagonality offDiagonalControlBridge := by
  exact ⟨diagonalComparator_isEffect,
    reachable_isRecordLaw,
    fun hrho hdiag => recordDiagonal_state_phase_weights_equal hrho hdiag,
    fun B hpres _ hM => reachable_recordPreserving_phase_weight_half B hpres hM,
    fun B _ hM hweight => nonhalf_weight_requires_nonzero_im B hM hweight,
    constantDiagonalControlBridge_preserves,
    offDiagonalControlBridge_not_preserving⟩

end

#print axioms recordDiagonal_state_sourceSelected_weight_half
#print axioms state_offdiag_conj
#print axioms state_sourceSelected_weight_commutator_form
#print axioms state_sourceSelected_weight_eq_half_sub_im
#print axioms diagonalComparator_isEffect
#print axioms recordDiagonal_state_negative_weight_half
#print axioms state_diagonalComparator_weight_half
#print axioms recordDiagonal_state_phase_weights_equal
#print axioms reachable_recordPreserving_phase_weight_half
#print axioms nonhalf_weight_requires_nonzero_im
#print axioms nonhalf_weight_refutes_recordPreservation
#print axioms constantDiagonalControlBridge_preserves
#print axioms constantDiagonalControlBridge_weight_half
#print axioms sourceSelectedGeneratedEffect_isEvent
#print axioms sourceSelectedGeneratedEffect_isState
#print axioms offDiagonalControlBridge_weight_one
#print axioms offDiagonalControlBridge_not_preserving
#print axioms sourcePhaseBornWeightBoundary_receipt

end EventAlgebra.SourcePhaseBornWeightBoundary
