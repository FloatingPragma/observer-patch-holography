import Mathlib.Data.Complex.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Data.Fin.VecNotation
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Ring

namespace OPH.MarkedPortRamsey

abbrev Axis := Fin 3 → ℝ
abbrev Doublet := Fin 2 → ℂ

/-- For a unit axis this is the spinor half-turn, -i n·σ. -/
noncomputable def halfTurn (n : Axis) (v : Doublet) : Doublet :=
  ![-Complex.I * ((n 2 : ℂ)*v 0 + (n 0 - Complex.I*n 1)*v 1),
    -Complex.I * ((n 0 + Complex.I*n 1)*v 0 - (n 2 : ℂ)*v 1)]

def basis (j : Fin 2) : Doublet := fun k => if k = j then 1 else 0

/-- A balanced launch and material analyzer report the real amplitude on
    all starts, even when other doublet or orthogonal exits are populated. -/
theorem bright_sub_dark (z : ℂ) :
    (Complex.normSq (1 + z) - Complex.normSq (1 - z)) / 4 = z.re := by
  simp [Complex.normSq_apply]
  ring

theorem signed_half_turn_pair (a b : Axis) (j : Fin 2) :
    ((halfTurn b (halfTurn a (basis j))) j).re =
      -(a 0 * b 0 + a 1 * b 1 + a 2 * b 2) := by
  fin_cases j <;>
    simp [halfTurn, basis, Complex.mul_re, Complex.mul_im] <;> ring

/-- Either resolved doublet address gives the same all-start fringe. -/
theorem routed_fringe (a b : Axis) (j : Fin 2) :
    let z := halfTurn b (halfTurn a (basis j)) j
    (Complex.normSq (1 + z) - Complex.normSq (1 - z)) / 4 =
      -(a 0 * b 0 + a 1 * b 1 + a 2 * b 2) :=
  (bright_sub_dark _).trans (signed_half_turn_pair a b j)

end OPH.MarkedPortRamsey
