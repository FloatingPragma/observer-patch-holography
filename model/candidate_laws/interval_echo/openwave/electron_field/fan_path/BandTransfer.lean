import Mathlib.Data.Complex.Basic
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

/-! # Two recorded amplitudes through a two-line spectral selection

The source and receiving matrices include the chosen clock and the spectator
vacuum read. Each of the two distinct admitted energy lines acquires a
nonzero phase during free evolution. The certified nonzero determinants of
the source and receiver therefore imply two independent accepted amplitudes
at every wait.
-/

namespace OPH.WrittenFan

theorem two_directions_survive (source receiver : Matrix (Fin 2) (Fin 2) ℂ)
    (phase : Fin 2 → ℂ) (hs : source.det ≠ 0) (hr : receiver.det ≠ 0)
    (hp : ∀ i, phase i ≠ 0) :
    (receiver * Matrix.diagonal phase * source).det ≠ 0 := by
  rw [Matrix.det_mul, Matrix.det_mul, Matrix.det_diagonal]
  exact mul_ne_zero (mul_ne_zero hr (Finset.prod_ne_zero_iff.mpr fun i _ => hp i)) hs

end OPH.WrittenFan
