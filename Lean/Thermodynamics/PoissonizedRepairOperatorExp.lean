import PoissonizedRepair
import Mathlib.Analysis.SpecialFunctions.Exponential

/-!
# Exponential of an idempotent repair generator

This file identifies the literal Banach-algebra exponential of a bounded
idempotent endomorphism.  It is kept separate from the algebraic
`poissonizedRepair` development so that the additional completeness and
continuity hypotheses are visible at the interface.
-/

namespace OPH.Thermodynamics

section GenericIdempotent

variable {A : Type*} [NormedRing A] [NormedAlgebra ℝ A] [CompleteSpace A]

/-- The Banach-algebra exponential of a scalar multiple of an idempotent has
the exact two-point spectral form. -/
theorem normedSpace_exp_smul_idempotent
    (a : ℝ) (P : A) (hP : IsIdempotentElem P) :
    NormedSpace.exp (a • P) = 1 + (Real.exp a - 1) • P := by
  rw [NormedSpace.exp_eq_tsum ℝ]
  change (∑' n : ℕ, ((n.factorial⁻¹ : ℝ) • (a • P) ^ n)) =
    1 + (Real.exp a - 1) • P
  rw [(NormedSpace.expSeries_summable' (𝕂 := ℝ) (a • P)).tsum_eq_zero_add]
  simp only [Nat.factorial_zero, Nat.cast_one, inv_one, pow_zero, one_smul]
  congr 1
  have hfull : Summable (fun n : ℕ => ((n.factorial⁻¹ : ℝ) * a ^ n)) := by
    simpa [smul_eq_mul] using
      (NormedSpace.expSeries_summable' (𝕂 := ℝ) a)
  have hscalar :
      (∑' n : ℕ, (((n + 1).factorial⁻¹ : ℝ) * a ^ (n + 1))) =
        Real.exp a - 1 := by
    rw [Real.exp_eq_exp_ℝ, NormedSpace.exp_eq_tsum ℝ]
    change (∑' n : ℕ, (((n + 1).factorial⁻¹ : ℝ) * a ^ (n + 1))) =
      (∑' n : ℕ, ((n.factorial⁻¹ : ℝ) * a ^ n)) - 1
    rw [hfull.tsum_eq_zero_add]
    simp
  have hsummable :
      Summable (fun n : ℕ => (((n + 1).factorial⁻¹ : ℝ) * a ^ (n + 1))) := by
    simpa [mul_comm] using hfull.comp_injective
      (fun _ _ h => Nat.succ.inj h)
  rw [← hscalar, ← hsummable.tsum_smul_const P]
  apply tsum_congr
  intro n
  rw [smul_pow, hP.pow_succ_eq]
  simp [smul_smul, mul_comm]

end GenericIdempotent

section ContinuousEndomorphism

variable {M : Type*} [NormedAddCommGroup M] [NormedSpace ℝ M] [CompleteSpace M]

/-- The bounded version of the declared repair generator. -/
noncomputable def continuousRepairGenerator
    (gamma : ℝ) (E : M →L[ℝ] M) : M →L[ℝ] M :=
  gamma • (E - 1)

/-- The bounded closed-form Poissonized repair map. -/
noncomputable def continuousPoissonizedRepair
    (gamma t : ℝ) (E : M →L[ℝ] M) : M →L[ℝ] M :=
  E + Real.exp (-gamma * t) • (1 - E)

omit [CompleteSpace M] in
@[simp]
theorem continuousRepairGenerator_apply
    (gamma : ℝ) (E : M →L[ℝ] M) (x : M) :
    continuousRepairGenerator gamma E x = gamma • (E x - x) := by
  simp [continuousRepairGenerator]

omit [CompleteSpace M] in
/-- Forgetting continuity recovers the existing algebraic repair generator. -/
theorem continuousRepairGenerator_toLinearMap
    (gamma : ℝ) (E : M →L[ℝ] M) :
    (continuousRepairGenerator gamma E).toLinearMap =
      repairGenerator gamma E.toLinearMap := by
  ext x
  simp [continuousRepairGenerator_apply, repairGenerator_apply]

omit [CompleteSpace M] in
@[simp]
theorem continuousPoissonizedRepair_apply
    (gamma t : ℝ) (E : M →L[ℝ] M) (x : M) :
    continuousPoissonizedRepair gamma t E x =
      E x + Real.exp (-gamma * t) • (x - E x) := by
  simp [continuousPoissonizedRepair]

omit [CompleteSpace M] in
/-- Forgetting continuity recovers the algebraic Poissonized repair map
definitionally up to the two identity-map interfaces. -/
theorem continuousPoissonizedRepair_toLinearMap
    (gamma t : ℝ) (E : M →L[ℝ] M) :
    (continuousPoissonizedRepair gamma t E).toLinearMap =
      poissonizedRepair gamma t E.toLinearMap := by
  ext x
  simp [continuousPoissonizedRepair_apply, poissonizedRepair_apply]

/-- Literal Mathlib operator-exponential identity for a bounded idempotent
repair.  This closes the gap between the algebraic/ODE flow and
`NormedSpace.exp` in the Banach algebra of continuous endomorphisms. -/
theorem normedSpace_exp_continuousRepairGenerator
    (gamma t : ℝ) (E : M →L[ℝ] M) (hE : E.comp E = E) :
    NormedSpace.exp (t • continuousRepairGenerator gamma E) =
      E + Real.exp (-gamma * t) • (1 - E) := by
  have hE' : IsIdempotentElem E := by
    change E.comp E = E
    exact hE
  have hQ : IsIdempotentElem (1 - E) := hE'.one_sub
  have hgenerator :
      t • continuousRepairGenerator gamma E = (-gamma * t) • (1 - E) := by
    simp only [continuousRepairGenerator]
    module
  rw [hgenerator, normedSpace_exp_smul_idempotent (-gamma * t) (1 - E) hQ]
  module

/-- The explicit operator-exponential identity expressed through the bounded
closed-form repair definition. -/
theorem normedSpace_exp_continuousRepairGenerator_eq_continuousPoissonizedRepair
    (gamma t : ℝ) (E : M →L[ℝ] M) (hE : E.comp E = E) :
    NormedSpace.exp (t • continuousRepairGenerator gamma E) =
      continuousPoissonizedRepair gamma t E := by
  rw [normedSpace_exp_continuousRepairGenerator gamma t E hE]
  rfl

/-- Pointwise bridge from the literal operator exponential to the existing
algebraic Poissonized repair. -/
theorem normedSpace_exp_continuousRepairGenerator_apply
    (gamma t : ℝ) (E : M →L[ℝ] M) (hE : E.comp E = E) (x : M) :
    NormedSpace.exp (t • continuousRepairGenerator gamma E) x =
      poissonizedRepair gamma t E.toLinearMap x := by
  rw [normedSpace_exp_continuousRepairGenerator_eq_continuousPoissonizedRepair
    gamma t E hE]
  exact LinearMap.congr_fun
    (continuousPoissonizedRepair_toLinearMap gamma t E) x

end ContinuousEndomorphism

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.normedSpace_exp_smul_idempotent
#print axioms OPH.Thermodynamics.normedSpace_exp_continuousRepairGenerator
#print axioms OPH.Thermodynamics.normedSpace_exp_continuousRepairGenerator_apply
