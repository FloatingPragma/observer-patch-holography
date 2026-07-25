import Mathlib

/-!
# Finite-capacity curvature behind the de Sitter shock-sign diagnostic

This module formalizes the exact algebraic core of the finite-screen
capacity calculation.  Sector dimensions are positive integers in the
finite model.  The gradient and Hessian statements below concern their
positive-real analytic relaxation.

For

`A(d₁, ..., dₙ) = ∑ᵢ (dᵢ / M) log dᵢ`, with `M = ∑ᵢ dᵢ`,

the Hessian at the symmetric point `dᵢ = d` acts by

`x ↦ (x - (2/n)(∑ᵢ xᵢ) 1) / (n d²)`.

The homogeneous vector has eigenvalue `-1/(n d²)`, while every fixed-`M`
tangent vector has eigenvalue `+1/(n d²)`.  The gradient is nonzero, so this
does not state an unconstrained interior maximum.  On the positive-real
transfer path, a maximum occurs at the one-sided budget boundary when a
fraction of a fixed total capacity is transferred from the horizon sectors
to an observer: the symmetric area drops by `log(1-f)`.  Integer-valued
finite transfers are the admissible finite-model points on this analytic
path.

The final section also checks the radius cancellation in pure de Sitter:
`μ² = ((D-2)/2)(2/L)L = D-2`, equal to the `ℓ = 1` spherical Laplacian
eigenvalue.

No theorem in this module identifies the analytic area functional with a
physical horizon, the screen kinetic operator with a graph Laplacian, or the
icosahedral rotation triplet with a gauge mode.  Those attachments remain
separate premises of the paper-level shock-spectrum diagnostic.
-/

namespace OPH.DeSitterCapacityShock

open scoped BigOperators

noncomputable section

/-- Generalized entropy of a finite sector distribution with relaxed positive
sector dimensions. -/
def generalizedSectorEntropy
    {α : Type*}
    [Fintype α]
    (probability dimension : α → ℝ) : ℝ :=
  -∑ i, probability i * Real.log (probability i) +
    ∑ i, probability i * Real.log (dimension i)

/-- Dimension-proportional sector weights give the exact finite entropy
`log M`, where `M` is the total dimension. -/
theorem generalizedSectorEntropy_dimension_weights
    {α : Type*}
    [Fintype α]
    (dimension : α → ℝ)
    {M : ℝ}
    (hM : 0 < M)
    (hdimension : ∀ i, 0 < dimension i)
    (hsum : ∑ i, dimension i = M) :
    generalizedSectorEntropy (fun i => dimension i / M) dimension =
      Real.log M := by
  have hlog (i : α) :
      Real.log (dimension i / M) =
        Real.log (dimension i) - Real.log M :=
    Real.log_div (ne_of_gt (hdimension i)) (ne_of_gt hM)
  simp only [generalizedSectorEntropy]
  simp_rw [hlog]
  rw [← Finset.sum_neg_distrib, ← Finset.sum_add_distrib]
  calc
    ∑ i, (-(dimension i / M *
          (Real.log (dimension i) - Real.log M)) +
        dimension i / M * Real.log (dimension i)) =
        ∑ i, dimension i / M * Real.log M := by
          apply Finset.sum_congr rfl
          intro i hi
          ring
    _ = (∑ i, dimension i / M) * Real.log M := by
          rw [Finset.sum_mul]
    _ = Real.log M := by
          rw [← Finset.sum_div, hsum]
          field_simp [ne_of_gt hM]

/-- Action of the symmetric-point Hessian in the positive-real relaxation of
the sector dimensions. -/
def analyticHessianAction
    (n : ℕ)
    (d : ℝ)
    (x : Fin n → ℝ)
    (i : Fin n) : ℝ :=
  (x i - (2 / (n : ℝ)) * ∑ j, x j) / ((n : ℝ) * d ^ 2)

/-- A fixed-total-horizon-capacity tangent vector is an eigenvector with the
positive transverse eigenvalue. -/
theorem analyticHessianAction_of_sum_zero
    {n : ℕ}
    {d : ℝ}
    (x : Fin n → ℝ)
    (hsum : ∑ j, x j = 0)
    (i : Fin n) :
    analyticHessianAction n d x i = x i / ((n : ℝ) * d ^ 2) := by
  simp [analyticHessianAction, hsum]

/-- The homogeneous direction is an eigenvector with the opposite-sign
eigenvalue. -/
theorem analyticHessianAction_homogeneous
    {n : ℕ}
    (hn : 0 < n)
    {d c : ℝ}
    (i : Fin n) :
    analyticHessianAction n d (fun _ => c) i =
      (-1 / ((n : ℝ) * d ^ 2)) * c := by
  have hn0 : (n : ℝ) ≠ 0 := by
    exact_mod_cast (Nat.ne_of_gt hn)
  simp only [analyticHessianAction, Finset.sum_const, Finset.card_univ,
    Fintype.card_fin, nsmul_eq_mul]
  field_simp
  ring

/-- The transverse Hessian eigenvalue is positive whenever the symmetric
dimension is nonzero. -/
theorem transverseEigenvalue_pos
    {n : ℕ}
    (hn : 0 < n)
    {d : ℝ}
    (hd : d ≠ 0) :
    0 < 1 / ((n : ℝ) * d ^ 2) := by
  have hnR : (0 : ℝ) < n := by
    exact_mod_cast hn
  have hd2 : 0 < d ^ 2 := sq_pos_of_ne_zero hd
  exact one_div_pos.mpr (mul_pos hnR hd2)

/-- The homogeneous Hessian eigenvalue is negative whenever the symmetric
dimension is nonzero. -/
theorem homogeneousEigenvalue_neg
    {n : ℕ}
    (hn : 0 < n)
    {d : ℝ}
    (hd : d ≠ 0) :
    -1 / ((n : ℝ) * d ^ 2) < 0 := by
  have hnR : (0 : ℝ) < n := by
    exact_mod_cast hn
  have hd2 : 0 < d ^ 2 := sq_pos_of_ne_zero hd
  exact div_neg_of_neg_of_pos (by norm_num) (mul_pos hnR hd2)

/-- Symmetric-point gradient component of the relaxed area functional. -/
def analyticGradientComponent (n : ℕ) (d : ℝ) : ℝ :=
  1 / ((n : ℝ) * d)

/-- The symmetric point is not an unconstrained stationary point. -/
theorem analyticGradientComponent_ne_zero
    {n : ℕ}
    (hn : 0 < n)
    {d : ℝ}
    (hd : d ≠ 0) :
    analyticGradientComponent n d ≠ 0 := by
  have hn0 : (n : ℝ) ≠ 0 := by
    exact_mod_cast (Nat.ne_of_gt hn)
  exact one_div_ne_zero (mul_ne_zero hn0 hd)

/-- Area coordinate on the symmetric ray. -/
def symmetricArea (d : ℝ) : ℝ :=
  Real.log d

/-- Positive-real interpolation of the symmetric-ray area after transferring
the fraction `f` of the capacity to an external observer ledger. -/
def capacityTransferArea (d f : ℝ) : ℝ :=
  Real.log ((1 - f) * d)

/-- Exact logarithmic area loss along the positive-real uniform-transfer
interpolation.  Admissible integer depletions are special cases. -/
theorem capacityTransferArea_sub
    {d f : ℝ}
    (hd : 0 < d)
    (hf : f < 1) :
    capacityTransferArea d f - symmetricArea d = Real.log (1 - f) := by
  have hremaining : 0 < 1 - f := sub_pos.mpr hf
  rw [capacityTransferArea, symmetricArea,
    Real.log_mul (ne_of_gt hremaining) (ne_of_gt hd)]
  ring

/-- Positive transfer fractions strictly decrease the symmetric-ray area.
The baseline `f = 0` is therefore a one-sided boundary maximum on
`0 ≤ f < 1`. -/
theorem capacityTransferArea_lt_baseline
    {d f : ℝ}
    (hd : 0 < d)
    (hf0 : 0 < f)
    (hf1 : f < 1) :
    capacityTransferArea d f < symmetricArea d := by
  have hremaining : 0 < 1 - f := sub_pos.mpr hf1
  have hremaining_lt_one : 1 - f < 1 := by linarith
  have hlog : Real.log (1 - f) < 0 :=
    Real.log_neg hremaining hremaining_lt_one
  have hdiff :=
    capacityTransferArea_sub (d := d) (f := f) hd hf1
  linarith

/-- Pure-de-Sitter shock mass parameter at spacetime dimension `D` and
radius `L`, written in the form used by the Killing-horizon equation. -/
def pureDeSitterMuSquared (D : ℕ) (L : ℝ) : ℝ :=
  (((D : ℝ) - 2) / 2) * (2 / L) * L

/-- The `ℓ = 1` scalar-Laplacian eigenvalue on `S^(D-2)`. -/
def lOneSphereEigenvalue (D : ℕ) : ℝ :=
  (D : ℝ) - 2

/-- The de Sitter radius cancels exactly, leaving the `ℓ = 1` eigenvalue. -/
theorem pureDeSitterMuSquared_eq_lOne
    (D : ℕ)
    {L : ℝ}
    (hL : L ≠ 0) :
    pureDeSitterMuSquared D L = lOneSphereEigenvalue D := by
  unfold pureDeSitterMuSquared lOneSphereEigenvalue
  rw [mul_assoc, div_mul_cancel₀ 2 hL]
  ring

end

end OPH.DeSitterCapacityShock
