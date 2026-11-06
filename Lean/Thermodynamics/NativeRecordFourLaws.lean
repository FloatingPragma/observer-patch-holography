import CofinalSpectralTailFamily

set_option autoImplicit false

open scoped BigOperators
open Filter Topology

namespace OPH.Thermodynamics.NativeRecord

/-!
# Four laws for a stopped native record

The record alphabet is the prefix code `0, 10, 11`. Its energy is the
number of native repair attempts in the record, multiplied by `log 2`.
Uniform two-bit tapes produce the record weights by finite counting.
The resulting full-resampling kernel has the same reference on its state
and transition sides and inhabits `FiniteFourLawCore` at multiplier one.

`model/native_record_thermodynamics` checks the binding of each observed bit
to the native endpoint source on every source state. This module proves the
finite prefix calculation and thermodynamic consequences. Energy and
inverse temperature are in record-cost units; no laboratory identification
or Gibbs law for the entire underlying source is asserted.
-/

/-- Complete stopped records, with `false` denoting zero. -/
def prefixRecord : Fin 3 → List Bool :=
  ![[false], [true, false], [true, true]]

/-- The number of native attempts recorded by the stopping rule. -/
def attemptCount (x : Fin 3) : ℕ := (prefixRecord x).length

/-- Energy is specified from attempt cost before any probability definition. -/
noncomputable def recordEnergy (x : Fin 3) : ℝ :=
  (attemptCount x : ℝ) * Real.log 2

/-- The first zero stops the first attempt; a first one consumes the second. -/
def sampleRecord (tape : Bool × Bool) : Fin 3 :=
  if tape.1 then if tape.2 then 2 else 1 else 0

/-- Count full two-bit tapes with the same stopped record. -/
def recordMultiplicity (x : Fin 3) : ℕ :=
  ((Finset.univ : Finset (Bool × Bool)).filter
    (fun tape => sampleRecord tape = x)).card

/-- Pushforward of the uniform law on four two-bit tapes. -/
noncomputable def recordWeight (x : Fin 3) : ℝ :=
  (recordMultiplicity x : ℝ) / 4

/-- The assigned record is the consumed prefix of its complete coin tape. -/
theorem sampleRecord_prefix (tape : Bool × Bool) :
    prefixRecord (sampleRecord tape) =
      if tape.1 then [tape.1, tape.2] else [tape.1] := by
  rcases tape with ⟨a, b⟩
  cases a <;> cases b <;> rfl

theorem attemptCount_eq : attemptCount = ![1, 2, 2] := by
  decide

theorem attemptCount_bounds (x : Fin 3) :
    1 ≤ attemptCount x ∧ attemptCount x ≤ 2 := by
  rw [attemptCount_eq]
  fin_cases x <;> norm_num

theorem recordMultiplicity_eq : recordMultiplicity = ![2, 1, 1] := by
  decide

theorem recordWeight_eq : recordWeight = ![1 / 2, 1 / 4, 1 / 4] := by
  funext x
  simp only [recordWeight, recordMultiplicity_eq]
  fin_cases x <;> norm_num

theorem recordWeight_pos (x : Fin 3) : 0 < recordWeight x := by
  rw [recordWeight_eq]
  fin_cases x <;> norm_num

theorem recordWeight_sum : ∑ x, recordWeight x = 1 := by
  rw [recordWeight_eq]
  norm_num [Fin.sum_univ_succ]

/-- The finite tape count equals the prefix probability from fair coins. -/
theorem recordWeight_eq_prefix_probability (x : Fin 3) :
    recordWeight x = (2 : ℝ)⁻¹ ^ attemptCount x := by
  rw [recordWeight_eq, attemptCount_eq]
  fin_cases x <;> norm_num

theorem log_two_pos : 0 < Real.log 2 := Real.log_pos (by norm_num)

theorem recordEnergy_eq :
    recordEnergy = ![Real.log 2, 2 * Real.log 2, 2 * Real.log 2] := by
  funext x
  simp only [recordEnergy, attemptCount_eq]
  fin_cases x <;> norm_num

/-- The thermal identification follows from the counted prefix probabilities. -/
theorem record_thermal (x : Fin 3) :
    -Real.log (recordWeight x) = recordEnergy x := by
  rw [recordWeight_eq_prefix_probability, Real.log_pow, Real.log_inv]
  unfold recordEnergy
  ring

theorem record_boltzmann_weight (x : Fin 3) :
    gibbsWeight recordEnergy 1 x = recordWeight x := by
  unfold gibbsWeight
  rw [← record_thermal]
  simpa using Real.exp_log (recordWeight_pos x)

theorem record_partition_one : partitionZ recordEnergy 1 = 1 := by
  unfold partitionZ
  simp_rw [record_boltzmann_weight]
  exact recordWeight_sum

/-- The complete resampling step forgets the previous record label. -/
noncomputable def recordRepairLaw : RepairLawData (Fin 3) Unit where
  ref := recordWeight
  ref_pos := recordWeight_pos
  ref_law := recordWeight_sum
  visible := fun _ => ()

/-- Exact record-cost calibration, with no remaining calibration hypothesis. -/
noncomputable def recordCalibration : CalibrationData recordRepairLaw where
  energy := recordEnergy
  beta := 1
  logZ := 0
  beta_pos := by norm_num
  thermal := by
    intro x
    simpa only [recordRepairLaw, one_mul, add_zero] using record_thermal x

/-- Every row is the independently counted prefix law. -/
theorem record_kernel_eq (x y : Fin 3) :
    FourLawSurface.repairKernel recordRepairLaw x y = recordWeight y := by
  simp [FourLawSurface.repairKernel, recordRepairLaw, heatBath, fiberMass,
    recordWeight_sum]

/-- A normalized input is sent exactly to the counted record reference. -/
theorem record_push_eq (p : Fin 3 → ℝ) (hp : ∑ x, p x = 1) :
    push p (FourLawSurface.repairKernel recordRepairLaw) = recordWeight := by
  funext y
  unfold push
  simp_rw [record_kernel_eq]
  rw [← Finset.sum_mul, hp, one_mul]

/-- All fields of the finite four-law core have concrete input data. -/
theorem record_four_laws : FiniteFourLawCore recordRepairLaw recordCalibration :=
  finiteFourLawCore recordRepairLaw recordCalibration

/-- The spectrum has a strictly positive, explicitly computed gap. -/
theorem record_energy_gap :
    recordEnergy 1 - recordEnergy 0 = Real.log 2 ∧
      0 < recordEnergy 1 - recordEnergy 0 := by
  rw [recordEnergy_eq]
  norm_num only [Matrix.cons_val_zero, Matrix.cons_val_one]
  constructor
  · ring
  · nlinarith [log_two_pos]

theorem record_energy_nondegenerate : recordEnergy 0 ≠ recordEnergy 1 := by
  have h := record_energy_gap.2
  linarith

theorem record_ground_gap (x : Fin 3) :
    recordEnergy x = Real.log 2 ∨
      Real.log 2 + Real.log 2 ≤ recordEnergy x := by
  rw [recordEnergy_eq]
  fin_cases x
  · exact Or.inl rfl
  · exact Or.inr (by ring_nf; rfl)
  · exact Or.inr (by ring_nf; rfl)

theorem record_ground_iff (x : Fin 3) :
    recordEnergy x = Real.log 2 ↔ x = 0 := by
  rw [recordEnergy_eq]
  fin_cases x <;> norm_num

theorem record_ground_nonempty :
    (Finset.univ.filter (fun x => recordEnergy x = Real.log 2)).Nonempty := by
  exact ⟨0, by simp [record_ground_iff]⟩

theorem record_ground_degeneracy :
    groundDegeneracy recordEnergy (Real.log 2) = 1 := by
  simp only [groundDegeneracy, record_ground_iff]
  decide

/-- The fixed record-cost spectrum has zero residual Gibbs entropy. -/
theorem record_third_entropy_limit :
    Tendsto (gibbsEntropy recordEnergy) atTop (𝓝 0) := by
  have h := record_four_laws.third_entropy_limit (Real.log 2) (Real.log 2)
    log_two_pos record_ground_gap record_ground_nonempty
  exact h.2.mpr record_ground_degeneracy

/-- The nonconstant cost spectrum makes its Gibbs multiplier identifiable. -/
theorem record_thermometer (beta1 beta2 : ℝ)
    (h : gibbs recordEnergy beta1 = gibbs recordEnergy beta2) : beta1 = beta2 :=
  record_four_laws.zeroth_thermometer beta1 beta2 0 1
    record_energy_nondegenerate h

/-- Any calibration with this same attempt-cost energy has multiplier one. -/
theorem record_calibration_multiplier_unique
    (C : CalibrationData recordRepairLaw) (h : C.energy = recordEnergy) :
    C.beta = 1 := by
  have hbeta := record_four_laws.zeroth_multiplier_unique C h 0 1
    record_energy_nondegenerate
  exact hbeta.symm

#print axioms recordWeight_eq_prefix_probability
#print axioms record_thermal
#print axioms record_partition_one
#print axioms record_kernel_eq
#print axioms record_push_eq
#print axioms record_four_laws
#print axioms record_energy_gap
#print axioms record_ground_degeneracy
#print axioms record_third_entropy_limit
#print axioms record_thermometer
#print axioms record_calibration_multiplier_unique

end OPH.Thermodynamics.NativeRecord
