import CofinalSpectralTailFamily

set_option autoImplicit false

open scoped BigOperators
open Filter Topology

namespace OPH.Thermodynamics.ProtectedGaugeCorner

/-!
# Finite thermodynamics of the declared protected Gauss-singlet corner

The complete four-outcome instrument in `model/protected_reference/gauge_join.py`
restricts, at each fixed protected memory character, to the replacement
channel with probabilities `(25/33, 8/33)`. Its two logical basis vectors are
the endpoint Gauss singlets. The dual reference, this reducing corner and
its admission are declared resources. No physical particle identification
is part of this statement.

The native endpoint coin supplies the integer weights 25 and 8. This module
proves the counting law, derives its diagonal kernel from all four retained
outcomes, and supplies every input of `FiniteFourLawCore`. Energy is the
same supported modular energy, in dimensionless units, with multiplier one.
Neither laboratory calibration nor a refinement-uniform low-temperature
family is inferred. The executable checks the quantum channel on every
matrix unit, including its extension by any untouched reference.
-/

/-- A uniform ticket realization of the native 25:8 endpoint coin. -/
def sampleCorner (ticket : Fin 33) : Fin 2 :=
  if ticket.val < 25 then 0 else 1

def cornerMultiplicity (x : Fin 2) : ℕ :=
  ((Finset.univ : Finset (Fin 33)).filter
    (fun ticket => sampleCorner ticket = x)).card

theorem cornerMultiplicity_eq : cornerMultiplicity = ![25, 8] := by
  decide

/-- The same source-derived native coin law, conditioned on the active pair. -/
noncomputable def cornerWeight (x : Fin 2) : ℝ :=
  (cornerMultiplicity x : ℝ) / 33

theorem cornerWeight_eq : cornerWeight = ![25 / 33, 8 / 33] := by
  funext x
  simp only [cornerWeight, cornerMultiplicity_eq]
  fin_cases x <;> norm_num

theorem cornerWeight_pos (x : Fin 2) : 0 < cornerWeight x := by
  rw [cornerWeight_eq]
  fin_cases x <;> norm_num

theorem cornerWeight_sum : ∑ x, cornerWeight x = 1 := by
  rw [cornerWeight_eq]
  norm_num [Fin.sum_univ_succ]

/-- Forward move, forward hold, reverse move and reverse hold, in that order. -/
def branchInput : Fin 4 → Fin 2 := ![1, 0, 0, 1]

def branchOutput : Fin 4 → Fin 2 := ![0, 0, 1, 1]

noncomputable def branchWeight : Fin 4 → ℝ :=
  ![25 / 33, 25 / 33, 8 / 33, 8 / 33]

/-- Sum every subnormalized retained outcome; no success conditioning. -/
noncomputable def cornerKernel (x y : Fin 2) : ℝ :=
  ∑ h : Fin 4,
    if branchInput h = x ∧ branchOutput h = y then branchWeight h else 0

theorem cornerKernel_eq (x y : Fin 2) :
    cornerKernel x y = cornerWeight y := by
  rw [cornerWeight_eq]
  fin_cases x <;> fin_cases y <;>
    norm_num [cornerKernel, branchInput, branchOutput, branchWeight,
      Fin.sum_univ_succ]

/-- Every outcome remains in the logical state space; its unread visible fibre
is the complete two-state corner. The outcome archive is traced only for
this marginal channel. -/
noncomputable def cornerRepairLaw : RepairLawData (Fin 2) Unit where
  ref := cornerWeight
  ref_pos := cornerWeight_pos
  ref_law := cornerWeight_sum
  visible := fun _ => ()

theorem corner_heatBath_eq (x y : Fin 2) :
    FourLawSurface.repairKernel cornerRepairLaw x y = cornerKernel x y := by
  rw [cornerKernel_eq]
  norm_num [FourLawSurface.repairKernel, cornerRepairLaw, heatBath, fiberMass,
    cornerWeight_eq]

theorem corner_push_eq (p : Fin 2 → ℝ) (hp : ∑ x, p x = 1) :
    push p cornerKernel = cornerWeight := by
  funext y
  unfold push
  simp_rw [cornerKernel_eq]
  rw [← Finset.sum_mul, hp, one_mul]

/-- Supported modular energy from the native integer counts. -/
noncomputable def cornerEnergy (x : Fin 2) : ℝ :=
  Real.log 33 - Real.log (cornerMultiplicity x)

theorem cornerEnergy_eq :
    cornerEnergy = ![Real.log 33 - Real.log 25, Real.log 33 - Real.log 8] := by
  funext x
  simp only [cornerEnergy, cornerMultiplicity_eq]
  fin_cases x <;> norm_num

theorem corner_thermal (x : Fin 2) :
    -Real.log (cornerWeight x) = cornerEnergy x := by
  rw [cornerWeight_eq, cornerEnergy_eq]
  fin_cases x <;> norm_num [Real.log_div]

/-- A modular-unit instance, with no claim of physical unit calibration. -/
noncomputable def cornerCalibration : CalibrationData cornerRepairLaw where
  energy := cornerEnergy
  beta := 1
  logZ := 0
  beta_pos := by norm_num
  thermal := by
    intro x
    simpa only [cornerRepairLaw, one_mul, add_zero] using corner_thermal x

/-- The finite core consumes the actual same-corner kernel and modular energy. -/
theorem corner_four_laws : FiniteFourLawCore cornerRepairLaw cornerCalibration :=
  finiteFourLawCore cornerRepairLaw cornerCalibration

/-- The finite conclusion and its actual complete-outcome kernel are joined. -/
theorem corner_same_kernel_four_laws :
    FiniteFourLawCore cornerRepairLaw cornerCalibration ∧
      FourLawSurface.repairKernel cornerRepairLaw = cornerKernel := by
  exact ⟨corner_four_laws, funext fun x => funext fun y => corner_heatBath_eq x y⟩

theorem corner_energy_gap :
    cornerEnergy 1 - cornerEnergy 0 = Real.log (25 / 8 : ℝ) := by
  rw [cornerEnergy_eq, Real.log_div (by norm_num) (by norm_num)]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one]
  ring

theorem corner_gap_pos : 0 < Real.log (25 / 8 : ℝ) :=
  Real.log_pos (by norm_num)

theorem corner_energy_nondegenerate : cornerEnergy 0 ≠ cornerEnergy 1 := by
  have h := corner_gap_pos
  rw [← corner_energy_gap] at h
  linarith

theorem corner_ground_gap (x : Fin 2) :
    cornerEnergy x = cornerEnergy 0 ∨
      cornerEnergy 0 + Real.log (25 / 8 : ℝ) ≤ cornerEnergy x := by
  fin_cases x
  · exact Or.inl rfl
  · right
    change cornerEnergy 0 + Real.log (25 / 8 : ℝ) ≤ cornerEnergy 1
    rw [← corner_energy_gap]
    linarith

theorem corner_ground_iff (x : Fin 2) :
    cornerEnergy x = cornerEnergy 0 ↔ x = 0 := by
  fin_cases x
  · simp
  · simp [Ne.symm corner_energy_nondegenerate]

theorem corner_ground_nonempty :
    (Finset.univ.filter (fun x => cornerEnergy x = cornerEnergy 0)).Nonempty := by
  exact ⟨0, by simp⟩

theorem corner_ground_degeneracy :
    groundDegeneracy cornerEnergy (cornerEnergy 0) = 1 := by
  simp only [groundDegeneracy, corner_ground_iff]
  decide

/-- This mathematical Gibbs-family limit concerns the fixed two-state
spectrum. No native variation of its multiplier or growing capacity is
implemented by the one fixed source coin. -/
theorem corner_third_entropy_limit :
    Tendsto (gibbsEntropy cornerEnergy) atTop (𝓝 0) := by
  have h := corner_four_laws.third_entropy_limit (cornerEnergy 0)
    (Real.log (25 / 8 : ℝ)) corner_gap_pos corner_ground_gap corner_ground_nonempty
  exact h.2.mpr corner_ground_degeneracy

theorem corner_thermometer (beta1 beta2 : ℝ)
    (h : gibbs cornerEnergy beta1 = gibbs cornerEnergy beta2) : beta1 = beta2 :=
  corner_four_laws.zeroth_thermometer beta1 beta2 0 1
    corner_energy_nondegenerate h

/-- Fixing the same counted modular energy fixes the Gibbs multiplier. -/
theorem corner_calibration_multiplier_unique
    (C : CalibrationData cornerRepairLaw) (h : C.energy = cornerEnergy) :
    C.beta = 1 := by
  have hbeta := corner_four_laws.zeroth_multiplier_unique C h 0 1
    corner_energy_nondegenerate
  exact hbeta.symm

#print axioms cornerMultiplicity_eq
#print axioms cornerKernel_eq
#print axioms corner_heatBath_eq
#print axioms corner_push_eq
#print axioms corner_thermal
#print axioms corner_four_laws
#print axioms corner_same_kernel_four_laws
#print axioms corner_energy_gap
#print axioms corner_third_entropy_limit
#print axioms corner_thermometer
#print axioms corner_calibration_multiplier_unique

end OPH.Thermodynamics.ProtectedGaugeCorner
