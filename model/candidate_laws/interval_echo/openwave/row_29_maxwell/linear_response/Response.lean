import Mathlib.Data.Complex.Basic
import Mathlib.LinearAlgebra.Matrix.ConjTranspose
import Mathlib.Tactic.Ring

namespace IntervalEcho.MaxwellRead

/-- The reduced receiver Gram matrix survives a source isometry, including
entangled incoming records and a source environment retained in its index. -/
theorem receiver_gram {s t r : Type*} [Fintype s] [Fintype t] [DecidableEq s]
    (U : Matrix t s ℂ) (A : Matrix s r ℂ) (hU : U.conjTranspose * U = 1) :
    (U * A).conjTranspose * (U * A) = A.conjTranspose * A := by
  rw [Matrix.conjTranspose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc U.conjTranspose,
      hU, Matrix.one_mul]

/-- Reversing the perturbation isolates a potential difference on each edge
when the reference amplitude is a real constant. -/
theorem edge_difference (g ar ai br bi : ℝ) :
    (((g + ar) * bi - ai * (g + br)) -
     ((g - ar) * (-bi) - (-ai) * (g - br))) / 2 = g * (bi - ai) := by ring

theorem face_difference (g a b c : ℝ) :
    g * (b - a) + g * (c - b) + g * (a - c) = 0 := by ring

end IntervalEcho.MaxwellRead
