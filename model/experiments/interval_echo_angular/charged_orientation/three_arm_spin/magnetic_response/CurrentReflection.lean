import Mathlib.Data.Fintype.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Linarith

/-! # Current in the two reflection eigenspaces

In a real basis of a stationary winding doublet, divide the Hermitian
bond-current matrix by the imaginary unit. This gives a real skew matrix.
In the basis of the reflection fixing the measured bond, its two basis
vectors acquire opposite signs. If the current commutes with this
reflection, every current matrix element vanishes.
-/

namespace OPH.ThreeArmSpin.CurrentReflection

def parity (i : Fin 2) : ℝ := if i = 0 then 1 else -1

theorem skew_current_zero (a : Fin 2 → Fin 2 → ℝ)
    (hskew : ∀ i j, a j i = -a i j)
    (hparity : ∀ i j, a i j * parity j = parity i * a i j) :
    ∀ i j, a i j = 0 := by
  have d₀ := hskew 0 0
  have d₁ := hskew 1 1
  have u := hparity 0 1
  have l := hparity 1 0
  simp [parity] at u l
  intro i j
  fin_cases i <;> fin_cases j
  · change a 0 0 = 0
    linarith
  · change a 0 1 = 0
    linarith
  · change a 1 0 = 0
    linarith
  · change a 1 1 = 0
    linarith

end OPH.ThreeArmSpin.CurrentReflection
