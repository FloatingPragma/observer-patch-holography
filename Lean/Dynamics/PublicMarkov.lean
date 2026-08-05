import EventAlgebra.PublicRecordAlgebra
import Mathlib.Data.Complex.BigOperators
import Mathlib.Topology.Connected.TotallyDisconnected

/-!
# Finite public dynamics: stochasticity and the continuity obstruction

The public record algebra of a projective partition is exactly a finite
function algebra on its active labels.  This file isolates the dynamics of
that finite classical algebra.

For complex-valued functions, `PointwiseNonnegative` means that every value
is real and has nonnegative real part.  A complex-linear map is called
positive when it preserves this cone, and unital when it fixes the constant
function `1`.  Under precisely those conventions, positive unital maps are
represented exactly by real row-stochastic kernels.

A continuous one-parameter permutation flow on a discrete finite record set
is trivial.  This is a finite public-algebra obstruction, not a claim about
continuous private matrix-algebra dynamics.  No source selects a particular
kernel, rate, clock, or private Hamiltonian here.
-/

namespace OPH.Dynamics

open scoped BigOperators

universe u

variable {ι : Type u} [Fintype ι] [DecidableEq ι]

/-! ## Positive unital maps are stochastic kernels -/

/-- The record basis function supported on the label `j`. -/
def recordBasis (j : ι) : ι → ℂ := fun i => if i = j then 1 else 0

/-- Positivity convention for the complex finite function algebra: every
coordinate is real and its real part is nonnegative. -/
def PointwiseNonnegative (f : ι → ℂ) : Prop :=
  ∀ i, 0 ≤ (f i).re ∧ (f i).im = 0

/-- A complex-linear map is positive when it preserves the declared
pointwise cone. -/
def IsPointwisePositive (T : (ι → ℂ) →ₗ[ℂ] (ι → ℂ)) : Prop :=
  ∀ f, PointwiseNonnegative f → PointwiseNonnegative (T f)

/-- A linear map on finite record functions is unital when it fixes the
constant function `1`. -/
def IsUnitalRecordMap (T : (ι → ℂ) →ₗ[ℂ] (ι → ℂ)) : Prop :=
  T 1 = 1

omit [Fintype ι] in
theorem recordBasis_pointwiseNonnegative (j : ι) :
    PointwiseNonnegative (recordBasis j) := by
  intro i
  by_cases h : i = j <;> simp [recordBasis, h]

/-- Every finite record function is the sum of its basis coordinates. -/
theorem sum_smul_recordBasis (f : ι → ℂ) :
    (∑ j, f j • recordBasis j) = f := by
  funext i
  simp [recordBasis]

@[simp]
theorem sum_recordBasis :
    (∑ j : ι, recordBasis j) = (1 : ι → ℂ) := by
  simpa using sum_smul_recordBasis (1 : ι → ℂ)

/-- The real kernel extracted from a linear map by applying it to record
basis functions.  Positivity will prove that no imaginary information was
discarded. -/
noncomputable def recordKernel (T : (ι → ℂ) →ₗ[ℂ] (ι → ℂ)) : ι → ι → ℝ :=
  fun x y => (T (recordBasis y) x).re

omit [Fintype ι] in
/-- Positivity makes every extracted kernel coefficient nonnegative. -/
theorem recordKernel_nonnegative
    (T : (ι → ℂ) →ₗ[ℂ] (ι → ℂ)) (hT : IsPointwisePositive T)
    (x y : ι) :
    0 ≤ recordKernel T x y :=
  (hT (recordBasis y) (recordBasis_pointwiseNonnegative y) x).1

omit [Fintype ι] in
/-- On a positive map, the image of a basis function is exactly the real
coefficient recorded by `recordKernel`. -/
theorem apply_recordBasis_eq_recordKernel
    (T : (ι → ℂ) →ₗ[ℂ] (ι → ℂ)) (hT : IsPointwisePositive T)
    (x y : ι) :
    T (recordBasis y) x = (recordKernel T x y : ℂ) := by
  apply Complex.ext
  · rfl
  · simpa [recordKernel] using
      (hT (recordBasis y) (recordBasis_pointwiseNonnegative y) x).2

/-- Unitality makes the extracted kernel row stochastic. -/
theorem recordKernel_row_sum
    (T : (ι → ℂ) →ₗ[ℂ] (ι → ℂ)) (hT : IsPointwisePositive T)
    (h1 : IsUnitalRecordMap T) (x : ι) :
    ∑ y, recordKernel T x y = 1 := by
  have him : ∀ y, T (recordBasis y) x = (recordKernel T x y : ℂ) :=
    fun y => apply_recordBasis_eq_recordKernel T hT x y
  have hsum : (∑ y, (recordKernel T x y : ℂ)) = 1 := by
    calc
      (∑ y, (recordKernel T x y : ℂ)) = ∑ y, T (recordBasis y) x := by
        apply Finset.sum_congr rfl
        intro y _
        exact (him y).symm
      _ = (∑ y, T (recordBasis y)) x := by rw [Finset.sum_apply]
      _ = T (∑ y, recordBasis y) x := by rw [map_sum]
      _ = T 1 x := by rw [sum_recordBasis]
      _ = 1 := by rw [h1]; rfl
  exact_mod_cast hsum

/-- Exact stochastic representation of a positive finite record map. -/
theorem apply_eq_recordKernel_sum
    (T : (ι → ℂ) →ₗ[ℂ] (ι → ℂ)) (hT : IsPointwisePositive T)
    (f : ι → ℂ) (x : ι) :
    T f x = ∑ y, (recordKernel T x y : ℂ) * f y := by
  calc
    T f x = T (∑ y, f y • recordBasis y) x := by
      rw [sum_smul_recordBasis]
    _ = (∑ y, T (f y • recordBasis y)) x := by rw [map_sum]
    _ = ∑ y, T (f y • recordBasis y) x := by rw [Finset.sum_apply]
    _ = ∑ y, f y * T (recordBasis y) x := by simp
    _ = ∑ y, f y * (recordKernel T x y : ℂ) := by
      apply Finset.sum_congr rfl
      intro y _
      rw [apply_recordBasis_eq_recordKernel T hT]
    _ = ∑ y, (recordKernel T x y : ℂ) * f y := by
      apply Finset.sum_congr rfl
      intro y _
      exact mul_comm _ _

/-- A real row-stochastic kernel in the observable convention
`(Tf)(x) = sum_y K(x,y) f(y)`. -/
def IsRowStochastic (K : ι → ι → ℝ) : Prop :=
  (∀ x y, 0 ≤ K x y) ∧ (∀ x, ∑ y, K x y = 1)

/-- The complex-linear record map defined by a real kernel. -/
noncomputable def recordMapOfKernel (K : ι → ι → ℝ) :
    (ι → ℂ) →ₗ[ℂ] (ι → ℂ) where
  toFun f x := ∑ y, (K x y : ℂ) * f y
  map_add' f g := by
    funext x
    simp [mul_add, Finset.sum_add_distrib]
  map_smul' c f := by
    funext x
    simp [Finset.mul_sum, mul_left_comm]

omit [DecidableEq ι] in
@[simp]
theorem recordMapOfKernel_apply (K : ι → ι → ℝ) (f : ι → ℂ) (x : ι) :
    recordMapOfKernel K f x = ∑ y, (K x y : ℂ) * f y :=
  rfl

omit [DecidableEq ι] in
/-- A nonnegative real kernel defines a positive record map. -/
theorem recordMapOfKernel_positive (K : ι → ι → ℝ)
    (hK : ∀ x y, 0 ≤ K x y) :
    IsPointwisePositive (recordMapOfKernel K) := by
  intro f hf x
  constructor
  · simp only [recordMapOfKernel_apply]
    rw [Complex.re_sum]
    simp only [Complex.mul_re, Complex.ofReal_re, Complex.ofReal_im,
      zero_mul, sub_zero]
    exact Finset.sum_nonneg fun y _ => mul_nonneg (hK x y) (hf y).1
  · simp only [recordMapOfKernel_apply]
    rw [Complex.im_sum]
    simp only [Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im,
      zero_mul]
    exact Finset.sum_eq_zero fun y _ => by simp [(hf y).2]

omit [DecidableEq ι] in
/-- A row-normalized real kernel defines a unital record map. -/
theorem recordMapOfKernel_unital (K : ι → ι → ℝ)
    (hK : ∀ x, ∑ y, K x y = 1) :
    IsUnitalRecordMap (recordMapOfKernel K) := by
  funext x
  simp only [recordMapOfKernel_apply, Pi.one_apply, mul_one]
  exact_mod_cast hK x

/-- Extracting the kernel of a nonnegative kernel map returns the original
real kernel. -/
theorem recordKernel_recordMapOfKernel (K : ι → ι → ℝ) :
    recordKernel (recordMapOfKernel K) = K := by
  funext x y
  simp [recordKernel, recordMapOfKernel_apply, recordBasis]

/-- Distinct real kernels define distinct record maps. -/
theorem recordMapOfKernel_injective :
    Function.Injective (recordMapOfKernel :
      (ι → ι → ℝ) → (ι → ℂ) →ₗ[ℂ] (ι → ℂ)) := by
  intro K L h
  rw [← recordKernel_recordMapOfKernel K,
    ← recordKernel_recordMapOfKernel L, h]

/-- Reconstructing a positive record map from its extracted kernel is exact. -/
theorem recordMapOfKernel_recordKernel
    (T : (ι → ℂ) →ₗ[ℂ] (ι → ℂ)) (hT : IsPointwisePositive T) :
    recordMapOfKernel (recordKernel T) = T := by
  apply LinearMap.ext
  intro f
  funext x
  exact (apply_eq_recordKernel_sum T hT f x).symm

/-- Classification theorem: positive unital complex-linear maps of a finite
record function algebra are exactly real row-stochastic matrices. -/
theorem positive_unital_iff_stochastic
    (T : (ι → ℂ) →ₗ[ℂ] (ι → ℂ)) :
    IsPointwisePositive T ∧ IsUnitalRecordMap T ↔
      ∃ K : ι → ι → ℝ, IsRowStochastic K ∧ T = recordMapOfKernel K := by
  constructor
  · rintro ⟨hT, h1⟩
    refine ⟨recordKernel T, ?_, (recordMapOfKernel_recordKernel T hT).symm⟩
    exact ⟨recordKernel_nonnegative T hT, recordKernel_row_sum T hT h1⟩
  · rintro ⟨K, ⟨hK0, hK1⟩, rfl⟩
    exact ⟨recordMapOfKernel_positive K hK0,
      recordMapOfKernel_unital K hK1⟩

/-! ## Application to the exact public record algebra -/

/-- The stochastic classification applies directly to the active-label
function algebra in `PublicRecordAlgebra.lean`.  Together with
`publicRecordFunctionEquiv`, this is the classical normal form for positive
unital public-record dynamics. -/
theorem activeRecord_positive_unital_iff_stochastic
    {n k : ℕ} (part : EventAlgebra.ProjectivePartition n k)
    (T : (part.ActiveIndex → ℂ) →ₗ[ℂ] (part.ActiveIndex → ℂ)) :
    IsPointwisePositive T ∧ IsUnitalRecordMap T ↔
      ∃ K : part.ActiveIndex → part.ActiveIndex → ℝ,
        IsRowStochastic K ∧ T = recordMapOfKernel K :=
  positive_unital_iff_stochastic T

/-! ## Continuous reversible public dynamics -/

/-- A continuous one-parameter flow of permutations of public record labels.
The group law is included in the interface although continuity and the
identity value suffice to force triviality. -/
structure ContinuousPermutationFlow (ι : Type u) [TopologicalSpace ι] where
  toPerm : ℝ → Equiv.Perm ι
  map_zero : toPerm 0 = Equiv.refl ι
  map_add : ∀ s t, toPerm (s + t) = toPerm s * toPerm t
  continuous_apply : ∀ i, Continuous fun t => toPerm t i

namespace ContinuousPermutationFlow

variable [TopologicalSpace ι] [DiscreteTopology ι]

omit [Fintype ι] [DecidableEq ι] in
/-- A continuous real one-parameter permutation flow on a discrete record
set fixes every label. -/
theorem apply_eq (F : ContinuousPermutationFlow ι) (t : ℝ) (i : ι) :
    F.toPerm t i = i := by
  have hconst : F.toPerm t i = F.toPerm 0 i :=
    PreconnectedSpace.constant (α := ℝ) (Y := ι) inferInstance
      (F.continuous_apply i)
  rw [F.map_zero] at hconst
  exact hconst

omit [Fintype ι] [DecidableEq ι] in
/-- Therefore the entire continuous permutation flow is the identity. -/
theorem toPerm_eq_refl (F : ContinuousPermutationFlow ι) (t : ℝ) :
    F.toPerm t = Equiv.refl ι := by
  ext i
  exact F.apply_eq t i

omit [Fintype ι] [DecidableEq ι] in
/-- The induced action on the finite public function algebra is pointwise
trivial. -/
theorem function_action_eq (F : ContinuousPermutationFlow ι)
    (t : ℝ) (f : ι → ℂ) :
    (fun i => f ((F.toPerm t).symm i)) = f := by
  rw [F.toPerm_eq_refl t]
  rfl

end ContinuousPermutationFlow

-- Axiom audit: declarations may use only the standard Mathlib basis.
#print axioms OPH.Dynamics.positive_unital_iff_stochastic
#print axioms OPH.Dynamics.activeRecord_positive_unital_iff_stochastic
#print axioms OPH.Dynamics.ContinuousPermutationFlow.toPerm_eq_refl
#print axioms OPH.Dynamics.ContinuousPermutationFlow.function_action_eq

end OPH.Dynamics
