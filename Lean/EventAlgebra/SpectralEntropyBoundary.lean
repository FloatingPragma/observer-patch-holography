import EventAlgebra.RecordMajorization
import Mathlib.Analysis.CStarAlgebra.ContinuousFunctionalCalculus.Projection
import Mathlib.Analysis.SpecialFunctions.ContinuousFunctionalCalculus.ExpLog.Basic

/-!
# B9: support boundary for totalized matrix logarithms

Lean's real logarithm is totalized by `Real.log 0 = 0`.  Consequently the
continuous-functional-calculus logarithm of every finite projection is the
zero matrix.  The superficially plausible raw formula

`Re Tr(ρ (log ρ - log σ))`

therefore returns zero for two orthogonal pure states.  A support-aware
Umegaki divergence must instead return `∞` when `supp ρ` is not contained in
`supp σ`.

This module proves that incompatibility exactly.  It is a fail-closed B9
boundary: the specialized full-support Gibbs-family logarithms already in
the repository cannot simply be totalized to discharge B9's support-aware
finite-matrix contract.
-/

namespace EventAlgebra

open Matrix

variable {n : ℕ}

section CFC

attribute [local instance] Matrix.linftyOpNormedAddCommGroup
  Matrix.linftyOpNormedRing Matrix.linftyOpNormedAlgebra

/-- The tempting totalized matrix logarithm. -/
noncomputable def totalizedMatrixLog
    (A : Matrix (Fin n) (Fin n) ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  CFC.log A

/-- The tempting raw real-valued Umegaki formula using the totalized log. -/
noncomputable def totalizedRelativeEntropy
    (A B : Matrix (Fin n) (Fin n) ℂ) : ℝ :=
  ((A * (totalizedMatrixLog A - totalizedMatrixLog B)).trace).re

/-- Totalization collapses the logarithm of every projection to zero because
both spectral values, zero and one, have logarithm zero in Lean's total real
logarithm. -/
theorem totalizedMatrixLog_event_eq_zero
    {P : Matrix (Fin n) (Fin n) ℂ} (hP : IsEvent P) :
    totalizedMatrixLog P = 0 := by
  rw [totalizedMatrixLog, CFC.log, ← cfc_zero ℝ P]
  apply cfc_congr
  intro x hx
  have hIdem : IsIdempotentElem P := hP.2
  have hx01 : x ∈ ({0, 1} : Set ℝ) :=
    hIdem.spectrum_subset ℝ hx
  simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hx01
  rcases hx01 with rfl | rfl <;> simp

/-- In particular, the raw totalized formula gives zero between any two
projection events, regardless of their supports. -/
theorem totalizedRelativeEntropy_events_eq_zero
    {P Q : Matrix (Fin n) (Fin n) ℂ} (hP : IsEvent P) (hQ : IsEvent Q) :
    totalizedRelativeEntropy P Q = 0 := by
  rw [totalizedRelativeEntropy, totalizedMatrixLog_event_eq_zero hP,
    totalizedMatrixLog_event_eq_zero hQ, sub_self, mul_zero,
    Matrix.trace_zero]
  rfl

end CFC

/-- Kernel inclusion formulation of `supp A ≤ supp B`: every vector killed
by `B` is also killed by `A`.  This is the only support fact needed by the
counterexample. -/
def supportContained (A B : Matrix (Fin n) (Fin n) ℂ) : Prop :=
  ∀ v : Fin n → ℂ, B *ᵥ v = 0 → A *ᵥ v = 0

/-- The first coordinate basis vector. -/
def binaryBasisZero : Fin 2 → ℂ := fun i => if i = 0 then 1 else 0

/-- The second coordinate projector kills the first basis vector. -/
theorem binary_proj_one_mulVec_basisZero :
    binaryCoordinatePartition.proj 1 *ᵥ binaryBasisZero = 0 := by
  funext i
  fin_cases i <;>
    simp [binaryCoordinatePartition, binaryBasisZero, Matrix.mulVec]

/-- The first coordinate projector fixes the first basis vector. -/
theorem binary_proj_zero_mulVec_basisZero :
    binaryCoordinatePartition.proj 0 *ᵥ binaryBasisZero = binaryBasisZero := by
  funext i
  fin_cases i <;>
    simp [binaryCoordinatePartition, binaryBasisZero, Matrix.mulVec]

theorem binaryBasisZero_ne_zero : binaryBasisZero ≠ 0 := by
  intro h
  have := congr_fun h (0 : Fin 2)
  norm_num [binaryBasisZero] at this

/-- Each coordinate projection has unit trace and hence is a density
matrix, not merely an abstract idempotent. -/
theorem binaryCoordinatePartition_proj_trace (i : Fin 2) :
    (binaryCoordinatePartition.proj i).trace = 1 := by
  rw [binaryCoordinatePartition, Matrix.trace_diagonal]
  fin_cases i <;> simp [Pi.single_apply]

theorem binaryCoordinatePartition_proj_isState (i : Fin 2) :
    IsState (binaryCoordinatePartition.proj i) :=
  ⟨(binaryCoordinatePartition.isEvent i).posSemidef,
    binaryCoordinatePartition_proj_trace i⟩

/-- Exact support failure for two orthogonal pure coordinate states. -/
theorem not_supportContained_binary_orthogonal :
    ¬ supportContained (binaryCoordinatePartition.proj 0)
      (binaryCoordinatePartition.proj 1) := by
  intro h
  have hz := h binaryBasisZero binary_proj_one_mulVec_basisZero
  rw [binary_proj_zero_mulVec_basisZero] at hz
  exact binaryBasisZero_ne_zero hz

/-- The witness really is a pair of orthogonal density matrices with failed
support inclusion. -/
theorem binary_orthogonal_density_receipt :
    IsState (binaryCoordinatePartition.proj 0) ∧
      IsState (binaryCoordinatePartition.proj 1) ∧
      binaryCoordinatePartition.proj 0 * binaryCoordinatePartition.proj 1 = 0 ∧
      ¬ supportContained (binaryCoordinatePartition.proj 0)
        (binaryCoordinatePartition.proj 1) := by
  refine ⟨binaryCoordinatePartition_proj_isState 0,
    binaryCoordinatePartition_proj_isState 1, ?_,
    not_supportContained_binary_orthogonal⟩
  exact binaryCoordinatePartition.orthogonal 0 1 (by norm_num)

/-- The raw formula assigns zero to the orthogonal pure-state pair. -/
theorem totalizedRelativeEntropy_binary_orthogonal_eq_zero :
    totalizedRelativeEntropy (binaryCoordinatePartition.proj 0)
      (binaryCoordinatePartition.proj 1) = 0 :=
  totalizedRelativeEntropy_events_eq_zero
    (binaryCoordinatePartition.isEvent 0)
    (binaryCoordinatePartition.isEvent 1)

/-- **Support-awareness no-go.** Any extended-real divergence that returns
`∞` on support failure cannot equal the finite raw totalized-log formula on
the orthogonal pure-state witness. -/
theorem supportAware_not_totalizedRelativeEntropy
    (D : Matrix (Fin 2) (Fin 2) ℂ → Matrix (Fin 2) (Fin 2) ℂ → WithTop ℝ)
    (hSupport : ∀ A B, ¬ supportContained A B → D A B = ⊤) :
    D (binaryCoordinatePartition.proj 0) (binaryCoordinatePartition.proj 1)
      ≠ (totalizedRelativeEntropy (binaryCoordinatePartition.proj 0)
          (binaryCoordinatePartition.proj 1) : WithTop ℝ) := by
  rw [hSupport _ _ not_supportContained_binary_orthogonal,
    totalizedRelativeEntropy_binary_orthogonal_eq_zero]
  exact WithTop.top_ne_coe

#print axioms EventAlgebra.totalizedMatrixLog_event_eq_zero
#print axioms EventAlgebra.binary_orthogonal_density_receipt
#print axioms EventAlgebra.totalizedRelativeEntropy_binary_orthogonal_eq_zero
#print axioms EventAlgebra.not_supportContained_binary_orthogonal
#print axioms EventAlgebra.supportAware_not_totalizedRelativeEntropy

end EventAlgebra
