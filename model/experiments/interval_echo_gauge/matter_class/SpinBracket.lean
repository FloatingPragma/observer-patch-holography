import Mathlib.Data.Complex.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Data.Fin.VecNotation
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Ring

namespace OPH.NativeMatterClass

abbrev Axis := Fin 3 → ℝ
abbrev Doublet := Fin 2 → ℂ

def cross (a b : Axis) : Axis :=
  ![a 1*b 2-a 2*b 1, a 2*b 0-a 0*b 2, a 0*b 1-a 1*b 0]

noncomputable def act (a : Axis) (v : Doublet) : Doublet :=
  ![-Complex.I/2 * ((a 2 : ℂ)*v 0 + (a 0 - Complex.I*a 1)*v 1),
    -Complex.I/2 * ((a 0 + Complex.I*a 1)*v 0 - (a 2 : ℂ)*v 1)]

/-- The addressed doublet represents the rank-three skew algebra. -/
theorem commutator (a b : Axis) (v : Doublet) :
    act a (act b v) - act b (act a v) = act (cross a b) v := by
  funext j
  fin_cases j <;>
    simp [act, cross] <;>
    ring_nf <;> simp [Complex.I_sq] <;> ring

/-- A central axis vanishes on the two active addresses. The third, odd
    address is a one-dimensional trivial summand at every axis. -/
theorem zero_axis (v : Doublet) : act 0 v = 0 := by
  funext j
  fin_cases j <;> simp [act]

end OPH.NativeMatterClass
