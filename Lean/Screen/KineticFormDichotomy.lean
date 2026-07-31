import Mathlib

namespace OPH.KineticFormDichotomy

/-!
# Kinetic-form dichotomy on the twelve-dimensional current algebra

The port-current pairing `B(f, h) = -Re tr(K(f) K(h))` on
`u(1) + su(2) + su(3)` is ad-invariant, so on each simple ideal it is a
single multiple of the Killing form.  A candidate kinetic form is therefore
fixed on the simple ideals by one Killing-relative rational coefficient per
ideal, and the su(2):su(3) coefficient ratio separates the two declared
candidates: the port-response trace form and the rank-fifteen matter trace
form.

This file checks rational arithmetic only: the Killing-relative coefficient
map `c = T / h_dual`, the two branch ratios `6` and `3/2` together with
their inequality, the frozen cofactor expansion of the matter-branch
concurrency determinant, and the row-scaling covariance of its zero locus.

Imported premises, not established here:

* the four-dimensional one-loop gauge-beta law and the declared
  `(nG, nH) = (3, 1)` coefficient column `(41/6, -19/6, -7)` from
  `RGRepresentationFrontier`;
* the port-branch indices `(T2, T3) = (2, 1/2)`, pinned from the exact
  Python certificate
  `code/angular_sprint/runtime/kinetic_form_selection_receipt.json`
  (producer `kinetic_form_selection_certificate.py`) on the declared
  charged-double-triplet current fixture (conditional input, not
  source-selected);
* the selection between the port branch and the matter branch, which stays
  open: no statement in this file selects either branch as physical.

No measured coupling value appears in this file.
-/

/-- Killing-relative su(2) coefficient of a trace form with total Dynkin
    index `T2`.  The convention is `c = T / h_dual` with su(2) dual Coxeter
    number `h_dual = 2`; equivalently `c = B / (-kappa)` with `kappa` the
    Killing form, and `T` stated in the standard `T(fund) = 1/2`
    convention. -/
def killingRelativeSU2 (T2 : ℚ) : ℚ :=
  T2 / 2

/-- Killing-relative su(3) coefficient of a trace form with total Dynkin
    index `T3`.  The convention is `c = T / h_dual` with su(3) dual Coxeter
    number `h_dual = 3`; equivalently `c = B / (-kappa)` with `kappa` the
    Killing form, and `T` stated in the standard `T(fund) = 1/2`
    convention. -/
def killingRelativeSU3 (T3 : ℚ) : ℚ :=
  T3 / 3

/-- su(2):su(3) Killing-relative coefficient ratio of the port-response
    trace form.  The indices `(T2, T3) = (2, 1/2)` are pinned from the
    exact Python certificate
    `code/angular_sprint/runtime/kinetic_form_selection_receipt.json`
    (producer `kinetic_form_selection_certificate.py`) on the declared
    charged-double-triplet current fixture: the kernel block realizes the
    complex adjoint (index `2`) and the even block the complex defining
    triplet (index `1/2`), giving `c2 = 1` and `c3 = 1/6` by exact
    Killing-Gram proportionality.  The indices enter as a conditional
    input, not as a source-selected result. -/
def portBranchRatio : ℚ :=
  killingRelativeSU2 2 / killingRelativeSU3 (1 / 2)

/-- su(2):su(3) Killing-relative coefficient ratio of the rank-fifteen
    matter trace form.  The indices `(T2, T3) = (2, 2)` are the
    `RGRepresentationFrontier` per-copy Weyl index sums
    `su2WeylIndexPerCopy` and `su3WeylIndexPerCopy`. -/
def matterBranchRatio : ℚ :=
  killingRelativeSU2 2 / killingRelativeSU3 2

/-- Exact evaluation of both branch ratios. -/
theorem branch_ratios : portBranchRatio = 6 ∧ matterBranchRatio = 3 / 2 := by
  norm_num [portBranchRatio, matterBranchRatio,
    killingRelativeSU2, killingRelativeSU3]

/-- The two candidate kinetic forms disagree on the su(2):su(3)
    Killing-relative coefficient ratio. -/
theorem branch_dichotomy : portBranchRatio ≠ matterBranchRatio := by
  norm_num [portBranchRatio, matterBranchRatio,
    killingRelativeSU2, killingRelativeSU3]

/-- Cofactor expansion of the matter-branch concurrency determinant.  The
    first column `x` is the sealed inverse-coupling vector; the second
    column is the matter kinetic-index column `k = (10/3, 2, 2)`; the third
    column is the declared `(nG, nH) = (3, 1)` one-loop coefficient column
    `b = (41/6, -19/6, -7)` from `RGRepresentationFrontier`.  This lemma
    freezes the statistic's exact cofactor coefficients before any
    comparison; no measured value appears. -/
theorem matter_branch_cofactor_form (x1 x2 x3 : ℚ) :
    Matrix.det !![x1, 10 / 3, 41 / 6; x2, 2, -19 / 6; x3, 2, -7] =
      (-23 / 3) * x1 + 37 * x2 + (-218 / 9) * x3 := by
  simp [Matrix.det_fin_three]
  ring

/-- The frozen cofactor statistic vanishes exactly when its integer
    rescaling by `-9` vanishes.  The `x` column is the sealed
    inverse-coupling vector; this restates the zero condition of
    `matter_branch_cofactor_form` with integer coefficients before any
    comparison; no measured value appears. -/
theorem matter_branch_concurrency_integer_form (x1 x2 x3 : ℚ) :
    ((-23 / 3) * x1 + 37 * x2 + (-218 / 9) * x3 = 0) ↔
      (69 * x1 - 333 * x2 + 218 * x3 = 0) := by
  constructor
  · intro h
    linear_combination (-9 : ℚ) * h
  · intro h
    linear_combination (-(1 : ℚ) / 9) * h

/-- Per-factor coupling renormalization rescales one determinant row per
    gauge factor; the determinant picks up the product of the three row
    factors. -/
theorem det_row_scaling (l1 l2 l3 : ℚ) (M : Matrix (Fin 3) (Fin 3) ℚ) :
    Matrix.det (Matrix.of fun i j => (![l1, l2, l3] i) * M i j) =
      l1 * l2 * l3 * Matrix.det M := by
  simp [Matrix.det_fin_three]
  ring

/-- Zero-set invariance of the frozen statistic under nonzero per-row
    rescaling: per-factor coupling renormalization rescales one determinant
    row per gauge factor, so the zero locus of the frozen statistic is
    normalization-covariant. -/
theorem det_row_scaling_zero_iff (l1 l2 l3 : ℚ)
    (M : Matrix (Fin 3) (Fin 3) ℚ)
    (h1 : l1 ≠ 0) (h2 : l2 ≠ 0) (h3 : l3 ≠ 0) :
    (Matrix.det (Matrix.of fun i j => (![l1, l2, l3] i) * M i j) = 0 ↔
      Matrix.det M = 0) := by
  rw [det_row_scaling]
  simp [mul_eq_zero, h1, h2, h3]

#print axioms branch_ratios
#print axioms branch_dichotomy
#print axioms matter_branch_cofactor_form
#print axioms matter_branch_concurrency_integer_form
#print axioms det_row_scaling
#print axioms det_row_scaling_zero_iff

end OPH.KineticFormDichotomy
