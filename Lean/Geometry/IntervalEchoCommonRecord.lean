import Geometry.CommonWorldInstrumentJoin

/-! A uniform eight-bit ticket writes one classical source rail and one repair
record. Tickets 179 through 255 flag refusal. The next three fair bits index the
instrument after the write, so every context has the same admitted fibres. -/

namespace OPH.IntervalEchoCommonRecord

open OPH.CommonWorldInstrumentJoin OPH.QFT EventAlgebra
open scoped BigOperators

def writeRail (ticket : Fin 256) : Option (Fin 2) :=
  if ticket.val < 111 then some 0 else if ticket.val < 179 then some 1 else none

def fibre (outcome : Option (Fin 2)) : Finset (Fin 256) :=
  Finset.univ.filter (fun t => writeRail t = outcome)

theorem ticket_counts :
    (fibre (some 0)).card = 111 ∧ (fibre (some 1)).card = 68 ∧
      (fibre none).card = 77 := by
  have h₀ : fibre (some 0) = Finset.Iio (111 : Fin 256) := by
    ext t
    simp only [fibre, Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_Iio]
    change (if t.val < 111 then some (0 : Fin 2) else
      if t.val < 179 then some 1 else none) = some 0 ↔ t.val < 111
    split_ifs <;> simp_all
  have h₁ : fibre (some 1) = Finset.Ico (111 : Fin 256) 179 := by
    ext t
    simp only [fibre, Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_Ico]
    change (if t.val < 111 then some (0 : Fin 2) else
      if t.val < 179 then some 1 else none) = some 1 ↔ 111 ≤ t.val ∧ t.val < 179
    split_ifs <;> simp_all
  have h₂ : fibre none = Finset.Ici (179 : Fin 256) := by
    ext t
    simp only [fibre, Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_Ici]
    change (if t.val < 111 then some (0 : Fin 2) else
      if t.val < 179 then some 1 else none) = none ↔ 179 ≤ t.val
    split_ifs <;> simp_all
    omega
  rw [h₀, h₁, h₂, Fin.card_Iio, Fin.card_Ico, Fin.card_Ici]
  decide

/-- The same recorded ticket drives both writes before the context choice. -/
def publicRecord (ticket : Fin 256) (selector : Fin 8) :=
  (writeRail ticket).map (fun rail => (rail, rail, selector))

theorem common_fibre (x : Fin 2) (c : Fin 8) :
    Finset.univ.filter (fun t => publicRecord t c = some (x, x, c)) =
      fibre (some x) := by
  ext t
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, fibre, publicRecord]
  cases writeRail t <;> simp

noncomputable def admittedMass (x : Fin 2) : ℝ :=
  ((fibre (some x)).card : ℝ) / ((fibre (some 0)).card + (fibre (some 1)).card)

theorem mass_eq_ref (x : Fin 2) : admittedMass x = runRepairLaw.ref x := by
  fin_cases x <;>
    norm_num [admittedMass, ticket_counts.1, ticket_counts.2.1, runRepairLaw, runRef]

noncomputable def prepared : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.diagonal (fun x => (admittedMass x : ℂ))

theorem prepared_eq_run : prepared = committedRunState := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [prepared, admittedMass, ticket_counts.1, ticket_counts.2.1,
      committedRunState, Matrix.diagonal]

theorem common_preparation :
    partitionPinching (diagonalPartition 2) prepared = prepared ∧
    ∀ x, (runRepairLaw.ref x : ℂ) = prepared x x := by
  rw [prepared_eq_run]
  exact ⟨diagonal_pinching_runState, runRepairLaw_ref_eq_runState⟩

/-- Put the prepared state in the common-world record itself. The instrument,
repair law, horizon and other architecture fields come from the typed witness. -/
noncomputable def sourceWorld : InstrumentedCommonWorldArchitecture :=
  { instrumentedCommittedWitness with
    state := ⟨prepared, by rw [prepared_eq_run]; exact committedRunState_isState⟩
    state_diagonal_fixed := common_preparation.1 }

theorem one_world_record :
    sourceWorld.state.1 = sourceWorld.decision.prep ∧
    (∀ x, admittedMass x = sourceWorld.repairLaw.ref x) ∧
    partitionPinching (diagonalPartition 2) sourceWorld.state.1 = sourceWorld.state.1 :=
  ⟨prepared_eq_run, mass_eq_ref, common_preparation.1⟩

#print axioms one_world_record
#print axioms common_fibre

end OPH.IntervalEchoCommonRecord
