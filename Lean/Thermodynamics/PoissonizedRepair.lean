import EventAlgebra.TwoScalePublicRepair

/-!
# Poissonized idempotent repair

For an idempotent real-linear repair `E`, when `gamma >= 0` and `t >= 0` a
Poisson clock of rate `gamma` leaves the state untouched when no repair occurs
and sends it to `E x` after one or more repairs.  Its exact finite-time map is
therefore

`T(t) = E + exp (-gamma * t) (I - E)`.

This file proves that closed form, the semigroup law, invariance of the public
component, and the kernel of the declared generator `L = gamma (E - I)` for
all real parameters as an algebraic continuation.  The literal Poisson
probability interpretation is restricted to `gamma >= 0` and `t >= 0`.
No physical clock calibration, positivity, or channel structure is inferred
from an arbitrary linear `E`.
-/

namespace OPH.Thermodynamics

variable {M : Type*} [AddCommGroup M] [Module ℝ M]

/-- The conditional-expectation generator `gamma (E - I)`. -/
def repairGenerator (gamma : ℝ) (E : M →ₗ[ℝ] M) : M →ₗ[ℝ] M :=
  gamma • (E - LinearMap.id)

@[simp]
theorem repairGenerator_apply (gamma : ℝ) (E : M →ₗ[ℝ] M) (x : M) :
    repairGenerator gamma E x = gamma • (E x - x) := by
  simp [repairGenerator]

/-- Exact repair map.  It has the literal Poisson-mixture interpretation only
for nonnegative `gamma` and nonnegative `t`; the definition is algebraically
well formed for all real parameters. -/
noncomputable def poissonizedRepair
    (gamma t : ℝ) (E : M →ₗ[ℝ] M) : M →ₗ[ℝ] M :=
  E + Real.exp (-gamma * t) • (LinearMap.id - E)

/-- Pointwise closed form, matching two-scale public relaxation. -/
theorem poissonizedRepair_apply (gamma t : ℝ) (E : M →ₗ[ℝ] M) (x : M) :
    poissonizedRepair gamma t E x =
      E x + Real.exp (-gamma * t) • (x - E x) := by
  simp [poissonizedRepair]

/-- Operator form `T(t) = E + exp(-gamma*t)(I-E)`. -/
theorem poissonizedRepair_eq_projection_plus_exp_residual
    (gamma t : ℝ) (E : M →ₗ[ℝ] M) :
    poissonizedRepair gamma t E =
      E + Real.exp (-gamma * t) • (LinearMap.id - E) :=
  rfl

@[simp]
theorem poissonizedRepair_zero (gamma : ℝ) (E : M →ₗ[ℝ] M) :
    poissonizedRepair gamma 0 E = LinearMap.id := by
  ext x
  simp [poissonizedRepair_apply]

/-- An idempotent repair is invariant after the Poissonized map. -/
theorem poissonizedRepair_average_invariant
    (gamma t : ℝ) (E : M →ₗ[ℝ] M) (hE : E.comp E = E) :
    E.comp (poissonizedRepair gamma t E) = E := by
  ext x
  rw [LinearMap.comp_apply, poissonizedRepair_apply]
  have hEE : E (E x) = E x := by
    have h := LinearMap.congr_fun hE x
    simpa [LinearMap.comp_apply] using h
  simp [hEE]

/-- Once publicized, the state stays fixed under Poissonized repair. -/
theorem poissonizedRepair_fixes_range
    (gamma t : ℝ) (E : M →ₗ[ℝ] M) (hE : E.comp E = E) :
    (poissonizedRepair gamma t E).comp E = E := by
  ext x
  rw [LinearMap.comp_apply, poissonizedRepair_apply]
  have hEE : E (E x) = E x := by
    have h := LinearMap.congr_fun hE x
    simpa [LinearMap.comp_apply] using h
  simp [hEE]

/-- Exact semigroup law for the Poissonized family. -/
theorem poissonizedRepair_add
    (gamma s t : ℝ) (E : M →ₗ[ℝ] M) (hE : E.comp E = E) :
    poissonizedRepair gamma (s + t) E =
      (poissonizedRepair gamma s E).comp (poissonizedRepair gamma t E) := by
  ext x
  rw [LinearMap.comp_apply]
  change EventAlgebra.publicRelaxTime gamma (s + t) E x =
    EventAlgebra.publicRelaxTime gamma s E
      (EventAlgebra.publicRelaxTime gamma t E x)
  exact EventAlgebra.publicRelaxTime_add gamma s t E hE x

/-- With nonzero rate, the generator vanishes exactly on the fixed space of
the repair. -/
theorem repairGenerator_eq_zero_iff
    {gamma : ℝ} (hgamma : gamma ≠ 0) (E : M →ₗ[ℝ] M) (x : M) :
    repairGenerator gamma E x = 0 ↔ E x = x := by
  rw [repairGenerator_apply]
  constructor
  · intro h
    have : E x - x = 0 := by
      exact (smul_eq_zero.mp h).resolve_left hgamma
    exact sub_eq_zero.mp this
  · intro h
    simp [h]

section Flow

variable {N : Type*} [NormedAddCommGroup N] [NormedSpace ℝ N]

/-- The closed form is differentiable and its derivative has the expected
exponentially decaying residual. -/
theorem hasDerivAt_poissonizedRepair_apply
    (gamma t : ℝ) (E : N →ₗ[ℝ] N) (x : N) :
    HasDerivAt (fun s : ℝ => poissonizedRepair gamma s E x)
      ((-gamma * Real.exp (-gamma * t)) • (x - E x)) t := by
  have hinner : HasDerivAt (fun s : ℝ => -gamma * s) (-gamma) t := by
    simpa using (hasDerivAt_id t).const_mul (-gamma)
  have hexp : HasDerivAt (fun s : ℝ => Real.exp (-gamma * s))
      (Real.exp (-gamma * t) * (-gamma)) t :=
    (Real.hasDerivAt_exp (-gamma * t)).comp t hinner
  have hscaled := hexp.smul_const (x - E x)
  have hconst : HasDerivAt (fun _s : ℝ => E x) 0 t :=
    hasDerivAt_const t (E x)
  convert hconst.add hscaled using 1
  simp [mul_comm]

/-- The Poissonized closed form solves `dT/dt = L T` for the declared
generator `L = gamma (E-I)` whenever `E` is idempotent.  Together with the
initial-value and semigroup theorems above, this is the exact exponentiated
flow statement available in the linear interface; no separate Banach-algebra
endomorphism exponential is bundled. -/
theorem hasDerivAt_poissonizedRepair_eq_generator
    (gamma t : ℝ) (E : N →ₗ[ℝ] N) (hE : E.comp E = E) (x : N) :
    HasDerivAt (fun s : ℝ => poissonizedRepair gamma s E x)
      (repairGenerator gamma E (poissonizedRepair gamma t E x)) t := by
  convert hasDerivAt_poissonizedRepair_apply gamma t E x using 1
  rw [repairGenerator_apply, poissonizedRepair_apply]
  have hEE : E (E x) = E x := by
    have h := LinearMap.congr_fun hE x
    simpa [LinearMap.comp_apply] using h
  simp [hEE]
  module

end Flow

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.poissonizedRepair_add
#print axioms OPH.Thermodynamics.repairGenerator_eq_zero_iff
#print axioms OPH.Thermodynamics.hasDerivAt_poissonizedRepair_eq_generator
