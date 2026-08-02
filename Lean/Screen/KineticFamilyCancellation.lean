import KineticFormDichotomy
import RGRepresentationFrontier

namespace OPH.KineticFamilyCancellation

open OPH.RGRepresentationFrontier

/-!
# Complete-family cancellation in the conditional kinetic-ray statistic

The rank-fifteen matter census has the per-copy representation-index column

`k = (10/3, 2, 2)`.

After importing the standard four-dimensional one-loop beta functional, one
complete matter family contributes `(20/9, 4/3, 4/3) = (2/3) k`.  Its
multiplicity therefore drops out of the alternating determinant
`det(x, k, b(nG,nH))` exactly.  This file proves that cancellation and exposes
the remaining scalar-count dependence as an exact plane.  At `nH = 1` the
plane specializes to the frozen #646 integer equation.

The results are conditional algebra for the declared matter-trace kinetic ray.
They do not select that ray as the physical kinetic action, derive a physical
family or scalar census, supply threshold or higher-loop corrections, or turn
the sealed inverse-coupling column into an OPH postdiction.
-/

/-- The complete-family matter-trace kinetic-index column, in the order
    `(U(1)_Y, SU(2), SU(3))`. -/
def matterKineticColumn : Fin 3 → ℚ :=
  ![10 / 3, 2, 2]

/-- The imported one-loop beta column at family parameter `nG` and scalar
    doublet parameter `nH`, in the same gauge-factor order. -/
def oneLoopBetaColumn (nG nH : ℚ) : Fin 3 → ℚ :=
  ![bY nG nH, b2 nG nH, b3 nG]

/-- A complete matter family shifts the imported one-loop beta column in the
    direction of the declared matter kinetic ray.  This is the structural
    reason its multiplicity cancels from the determinant below. -/
theorem beta_column_complete_family_decomposition (nG nH : ℚ) :
    oneLoopBetaColumn nG nH =
      oneLoopBetaColumn 0 nH + (2 * nG / 3) • matterKineticColumn := by
  funext i
  fin_cases i <;>
    simp [oneLoopBetaColumn, matterKineticColumn, bY, b2, b3] <;>
    ring

/-- The conditional scale-free determinant with a sealed first column `x`,
    the matter-trace kinetic ray as its second column, and the imported
    one-loop beta coefficients as its third column. -/
def matterBranchDet (x1 x2 x3 nG nH : ℚ) : ℚ :=
  Matrix.det !![
    x1, 10 / 3, bY nG nH;
    x2, 2, b2 nG nH;
    x3, 2, b3 nG]

/-- Exact cancellation of complete-family multiplicity.  The theorem is
    stated over rational `nG`, so in particular it covers every natural family
    count admitted by the imported one-loop formula. -/
theorem complete_family_multiplicity_cancels
    (x1 x2 x3 nG nH : ℚ) :
    matterBranchDet x1 x2 x3 nG nH =
      matterBranchDet x1 x2 x3 0 nH := by
  simp [matterBranchDet, Matrix.det_fin_three, bY, b2, b3]
  ring

/-- Exact cofactor form after complete-family cancellation.  Only the scalar
    doublet parameter remains. -/
theorem general_scalar_count_cofactor_form
    (x1 x2 x3 nG nH : ℚ) :
    matterBranchDet x1 x2 x3 nG nH =
      (-(nH + 22) / 3) * x1 +
        ((nH + 110) / 3) * x2 +
        (2 * (nH - 110) / 9) * x3 := by
  simp [matterBranchDet, Matrix.det_fin_three, bY, b2, b3]
  ring

/-- The determinant vanishes exactly on this scalar-count-dependent plane.
    No division by `nH` is used, so the statement holds for every rational
    scalar parameter. -/
theorem general_scalar_count_plane
    (x1 x2 x3 nG nH : ℚ) :
    (matterBranchDet x1 x2 x3 nG nH = 0) ↔
      (3 * (nH + 22) * x1 -
          3 * (nH + 110) * x2 -
          2 * (nH - 110) * x3 = 0) := by
  rw [general_scalar_count_cofactor_form]
  constructor
  · intro h
    linear_combination (-9 : ℚ) * h
  · intro h
    linear_combination (-(1 : ℚ) / 9) * h

/-- With one scalar doublet, every complete-family multiplicity gives the
    frozen rational cofactor form used by the #646 statistic. -/
theorem one_scalar_doublet_cofactor_form
    (x1 x2 x3 nG : ℚ) :
    matterBranchDet x1 x2 x3 nG 1 =
      (-23 / 3) * x1 + 37 * x2 + (-218 / 9) * x3 := by
  rw [general_scalar_count_cofactor_form]
  ring

/-- The one-scalar determinant for every complete-family multiplicity is
    exactly the already frozen `(nG,nH) = (3,1)` determinant.  This identifies
    the general cancellation theorem with the #646 statistic definition, not
    with a physical coupling comparison. -/
theorem one_scalar_doublet_matches_frozen_three_one
    (x1 x2 x3 nG : ℚ) :
    matterBranchDet x1 x2 x3 nG 1 =
      Matrix.det !![
        x1, 10 / 3, 41 / 6;
        x2, 2, -19 / 6;
        x3, 2, -7] := by
  calc
    matterBranchDet x1 x2 x3 nG 1 =
        (-23 / 3) * x1 + 37 * x2 + (-218 / 9) * x3 :=
      one_scalar_doublet_cofactor_form x1 x2 x3 nG
    _ = Matrix.det !![
          x1, 10 / 3, 41 / 6;
          x2, 2, -19 / 6;
          x3, 2, -7] :=
      (OPH.KineticFormDichotomy.matter_branch_cofactor_form x1 x2 x3).symm

/-- The `nH = 1` specialization of the general plane is the frozen integer
    equation `69*x1 - 333*x2 + 218*x3 = 0`, independently of `nG`. -/
theorem one_scalar_doublet_plane
    (x1 x2 x3 nG : ℚ) :
    (matterBranchDet x1 x2 x3 nG 1 = 0) ↔
      (69 * x1 - 333 * x2 + 218 * x3 = 0) := by
  rw [one_scalar_doublet_cofactor_form]
  constructor
  · intro h
    linear_combination (-9 : ℚ) * h
  · intro h
    linear_combination (-(1 : ℚ) / 9) * h

#print axioms beta_column_complete_family_decomposition
#print axioms complete_family_multiplicity_cancels
#print axioms general_scalar_count_cofactor_form
#print axioms general_scalar_count_plane
#print axioms one_scalar_doublet_cofactor_form
#print axioms one_scalar_doublet_matches_frozen_three_one
#print axioms one_scalar_doublet_plane

end OPH.KineticFamilyCancellation
