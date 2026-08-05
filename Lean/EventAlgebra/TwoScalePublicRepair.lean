import Mathlib

/-!
# Two-scale public repair

This module isolates the exact algebra of relaxation toward an idempotent
publicization map.  If `E` is idempotent, then

`publicRelax a E x = E x + a • (x - E x)`

keeps the public component `E x` fixed and multiplies the unresolved
component by `a`.  Consequently two relaxation steps compose by multiplying
their unresolved amplitudes.  The exponential specialization therefore is a
semigroup under addition of its time parameters.

The results are purely linear.  They do not identify `E` with a physical
channel, do not establish complete positivity, and do not attach the real
parameters to laboratory time or a source-derived rate.
-/

namespace EventAlgebra

section Algebraic

variable {R M : Type*} [CommRing R] [AddCommGroup M] [Module R M]

/-- Relaxation toward the range of `E`, with residual amplitude `a`. -/
def publicRelax (a : R) (E : M →ₗ[R] M) (x : M) : M :=
  E x + a • (x - E x)

@[simp]
theorem publicRelax_zero (E : M →ₗ[R] M) (x : M) :
    publicRelax 0 E x = E x := by
  simp [publicRelax]

@[simp]
theorem publicRelax_one (E : M →ₗ[R] M) (x : M) :
    publicRelax 1 E x = x := by
  simp [publicRelax]

/-- The exact public/residual decomposition. -/
theorem publicRelax_decomposition (a : R) (E : M →ₗ[R] M) (x : M) :
    publicRelax a E x = a • x + (1 - a) • E x := by
  simp only [publicRelax, smul_sub, one_smul, sub_smul]
  module

/-- Applying an idempotent publicization map after relaxation erases the
residual and returns the original public component. -/
theorem publicRelax_average_invariant (a : R) (E : M →ₗ[R] M)
    (hE : E.comp E = E) (x : M) :
    E (publicRelax a E x) = E x := by
  have hEE : E (E x) = E x := by
    have h := LinearMap.congr_fun hE x
    simpa [LinearMap.comp_apply] using h
  simp [publicRelax, hEE]

/-- The public component itself is fixed by every relaxation amplitude. -/
theorem publicRelax_fixes_range (a : R) (E : M →ₗ[R] M)
    (hE : E.comp E = E) (x : M) :
    publicRelax a E (E x) = E x := by
  have hEE : E (E x) = E x := by
    have h := LinearMap.congr_fun hE x
    simpa [LinearMap.comp_apply] using h
  simp [publicRelax, hEE]

/-- Relaxation scales the unresolved residual exactly. -/
theorem publicRelax_residual (a : R) (E : M →ₗ[R] M)
    (hE : E.comp E = E) (x : M) :
    publicRelax a E x - E (publicRelax a E x) = a • (x - E x) := by
  rw [publicRelax_average_invariant a E hE x]
  simp [publicRelax]

/-- Exact composition law: unresolved amplitudes multiply. -/
theorem publicRelax_compose (a b : R) (E : M →ₗ[R] M)
    (hE : E.comp E = E) (x : M) :
    publicRelax a E (publicRelax b E x) = publicRelax (a * b) E x := by
  have hEE : E (E x) = E x := by
    have h := LinearMap.congr_fun hE x
    simpa [LinearMap.comp_apply] using h
  simp [publicRelax, hEE, mul_smul]

/-- Any linear statistic invariant under `E` is invariant under relaxation. -/
theorem publicRelax_preserves_linear_statistic
    (a : R) (E : M →ₗ[R] M) (f : M →ₗ[R] R)
    (hf : f.comp E = f) (x : M) :
    f (publicRelax a E x) = f x := by
  have hfE : f (E x) = f x := by
    have h := LinearMap.congr_fun hf x
    simpa [LinearMap.comp_apply] using h
  simp [publicRelax, hfE]

end Algebraic

section Exponential

variable {M : Type*} [AddCommGroup M] [Module ℝ M]

/-- Exponential relaxation with declared rate `gamma` and parameter `t`. -/
noncomputable def publicRelaxTime
    (gamma t : ℝ) (E : M →ₗ[ℝ] M) (x : M) : M :=
  publicRelax (Real.exp (-gamma * t)) E x

@[simp]
theorem publicRelaxTime_zero (gamma : ℝ) (E : M →ₗ[ℝ] M) (x : M) :
    publicRelaxTime gamma 0 E x = x := by
  simp [publicRelaxTime]

/-- The exponential family is a semigroup for an idempotent `E`. -/
theorem publicRelaxTime_add (gamma s t : ℝ) (E : M →ₗ[ℝ] M)
    (hE : E.comp E = E) (x : M) :
    publicRelaxTime gamma (s + t) E x =
      publicRelaxTime gamma s E (publicRelaxTime gamma t E x) := by
  rw [publicRelaxTime, publicRelaxTime, publicRelaxTime,
    publicRelax_compose]
  · congr 2
    rw [← Real.exp_add]
    congr 1
    ring
  · exact hE

/-- The public component is time-independent along exponential relaxation. -/
theorem publicRelaxTime_average_invariant (gamma t : ℝ)
    (E : M →ₗ[ℝ] M) (hE : E.comp E = E) (x : M) :
    E (publicRelaxTime gamma t E x) = E x :=
  publicRelax_average_invariant _ E hE x

/-- Closed-form exponential residual law. -/
theorem publicRelaxTime_residual (gamma t : ℝ)
    (E : M →ₗ[ℝ] M) (hE : E.comp E = E) (x : M) :
    publicRelaxTime gamma t E x - E (publicRelaxTime gamma t E x) =
      Real.exp (-gamma * t) • (x - E x) :=
  publicRelax_residual _ E hE x

end Exponential

end EventAlgebra

#print axioms EventAlgebra.publicRelax_compose
#print axioms EventAlgebra.publicRelaxTime_add
