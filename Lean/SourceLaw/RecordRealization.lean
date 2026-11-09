import Mathlib.Data.Real.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Complete record realizations from reflected moments

Let the columns of `X` be source conditional expectations of accessible
record cylinders and let `W` be the source probability weight matrix. For
invertible `X` and `W`, the observed Gram `Xᵀ W X` and insertion moments
`Xᵀ W L X` recover each selective source map as `X⁻¹ L X`. The same
similarity gives every finite chronological word and its reflected moment.

The hypotheses are explicit matrix inverse identities. The executable
`model/native_record_response/short_records.py` certificate supplies their
66-dimensional native population instance; its numerical source matrices
are not kernel-evaluated in this module. No quantum tomography, source
selection, statistical error bound or physical energy interpretation follows
from these algebraic statements alone.
-/

namespace OPH.SourceLaw.RecordRealization

variable {n : Type*} [Fintype n] [DecidableEq n]
variable {X Y W Z : Matrix n n ℝ}

/-- The two-sided inverse of a record Gram matrix from explicit source inverses. -/
theorem gram_inverse
    (hXY : X * Y = 1) (hYX : Y * X = 1)
    (hWZ : W * Z = 1) (hZW : Z * W = 1) :
    (Y * Z * Y.transpose) * (X.transpose * W * X) = 1 ∧
    (X.transpose * W * X) * (Y * Z * Y.transpose) = 1 := by
  constructor
  · calc
      _ = Y * Z * (X * Y).transpose * W * X := by
        simp only [Matrix.transpose_mul, Matrix.mul_assoc]
      _ = 1 := by
        rw [hXY, Matrix.transpose_one, Matrix.mul_one]
        simp only [Matrix.mul_assoc, hZW, Matrix.one_mul, hYX]
  · calc
      _ = X.transpose * W * (X * Y) * Z * Y.transpose := by
        simp only [Matrix.mul_assoc]
      _ = X.transpose * Y.transpose := by
        rw [hXY, Matrix.mul_one]
        simp only [Matrix.mul_assoc, hWZ, Matrix.one_mul]
      _ = 1 := by rw [← Matrix.transpose_mul, hYX, Matrix.transpose_one]

/-- An inserted raw-record matrix is recovered by the reflected Gram inverse. -/
theorem recover_insertion (hXY : X * Y = 1) (hZW : Z * W = 1)
    (L : Matrix n n ℝ) :
    (Y * Z * Y.transpose) * (X.transpose * W * L * X) = Y * L * X := by
  calc
    _ = Y * Z * (X * Y).transpose * W * L * X := by
      simp only [Matrix.transpose_mul, Matrix.mul_assoc]
    _ = Y * L * X := by
      rw [hXY, Matrix.transpose_one, Matrix.mul_one]
      simp only [Matrix.mul_assoc, hZW, Matrix.one_mul]

/-- Recovery uses the ordinary nonsingular inverse of the observed Gram matrix. -/
theorem recover_from_gram
    (hXY : X * Y = 1) (hYX : Y * X = 1)
    (hWZ : W * Z = 1) (hZW : Z * W = 1) (L : Matrix n n ℝ) :
    (X.transpose * W * X)⁻¹ * (X.transpose * W * L * X) = Y * L * X := by
  rw [Matrix.inv_eq_left_inv (gram_inverse hXY hYX hWZ hZW).1]
  exact recover_insertion hXY hZW L

/-- Exact recovery determines the insertion uniquely on the represented source space. -/
theorem insertion_unique (hXY : X * Y = 1) (hZW : Z * W = 1)
    {L K : Matrix n n ℝ}
    (h : X.transpose * W * L * X = X.transpose * W * K * X) : L = K := by
  have e : Y * L * X = Y * K * X := by
    rw [← recover_insertion hXY hZW L, ← recover_insertion hXY hZW K, h]
  have e' := congrArg (fun A => X * A * Y) e
  simpa only [Matrix.mul_assoc, ← Matrix.mul_assoc X Y, hXY,
    Matrix.one_mul, Matrix.mul_one] using e'

variable {α : Type*}

/-- Chronological backward insertions, including empty histories. -/
def history (L : α → Matrix n n ℝ) (word : List α) : Matrix n n ℝ :=
  (word.map L).prod

/-- All held-out words obey the same reconstruction, at arbitrary finite length. -/
theorem recover_history (hXY : X * Y = 1) (hYX : Y * X = 1)
    (L : α → Matrix n n ℝ) (word : List α) :
    history (fun a => Y * L a * X) word = Y * history L word * X := by
  induction word with
  | nil => simp [history, hYX]
  | cons a word ih =>
      simp only [history, List.map_cons, List.prod_cons] at *
      rw [ih]
      calc
        _ = Y * L a * (X * Y) * (List.map L word).prod * X := by
          simp only [Matrix.mul_assoc]
        _ = Y * (L a * (List.map L word).prod) * X := by
          rw [hXY, Matrix.mul_one]
          simp only [Matrix.mul_assoc]

/-- Every reflected history moment follows from the recovered selective maps. -/
theorem gram_history (hXY : X * Y = 1) (hYX : Y * X = 1)
    (L : α → Matrix n n ℝ) (word : List α) :
    (X.transpose * W * X) * history (fun a => Y * L a * X) word =
      X.transpose * W * history L word * X := by
  rw [recover_history hXY hYX]
  calc
    _ = X.transpose * W * (X * Y) * history L word * X := by
      simp only [Matrix.mul_assoc]
    _ = _ := by rw [hXY, Matrix.mul_one]

#print axioms gram_inverse
#print axioms recover_from_gram
#print axioms gram_history
#print axioms recover_insertion
#print axioms insertion_unique
#print axioms recover_history

end OPH.SourceLaw.RecordRealization
