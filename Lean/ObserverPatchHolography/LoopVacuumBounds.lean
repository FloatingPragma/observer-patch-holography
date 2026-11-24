import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Real.Basic
import Mathlib.Logic.Equiv.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum

/-! Finite bounds used in the two-loop vacuum comparison.
The probability and Dirichlet-form comparisons are hypotheses of the scalar
lemmas. The source calculation identifies the observables and proves them. -/

namespace ObserverPatchHolography.LoopVacuumBounds

/-- Two triangular vacancy loops acting on their three possible neighboring marks. -/
theorem shared_spoke_orders :
    (Equiv.swap (0 : Fin 3) 1) ((Equiv.swap (1 : Fin 3) 2) 1) ≠
    (Equiv.swap (1 : Fin 3) 2) ((Equiv.swap (0 : Fin 3) 1) 1) := by decide

/-- With one vacancy, a uniformly chosen root has at most inverse-size
mean binary variance for a loop click supported on that root. -/
theorem single_vacancy_variance {n : ℕ} (hn : 0 < n)
    (click vacancy : Fin n → ℝ) (h : ∀ i, click i ≤ vacancy i)
    (htotal : ∑ i, vacancy i = 1) :
    (∑ i, click i * (1 - click i)) / n ≤ 1 / (n : ℝ) := by
  apply (div_le_div_iff_of_pos_right (by exact_mod_cast hn : 0 < (n : ℝ))).mpr
  calc
    ∑ i, click i * (1 - click i) ≤ ∑ i, vacancy i :=
      Finset.sum_le_sum (fun i _ => by nlinarith [h i, sq_nonneg (click i)])
    _ = 1 := htotal

/-- For any root-supported binary read, its mean variance is bounded by the
mean vacancy density, also in a preparation with fluctuating vacancy count. -/
theorem vacancy_density_variance {n : ℕ} (hn : 0 < n)
    (click vacancy : Fin n → ℝ) (h : ∀ i, click i ≤ vacancy i) :
    (∑ i, click i * (1 - click i)) / n ≤ (∑ i, vacancy i) / (n : ℝ) := by
  apply (div_le_div_iff_of_pos_right (by exact_mod_cast hn : 0 < (n : ℝ))).mpr
  exact Finset.sum_le_sum (fun i _ => by nlinarith [h i, sq_nonneg (click i)])

/-- Transferring a Poincaré floor through two bounds on the ground density. -/
theorem density_gap_floor {c low high varH varG energyH energyG : ℝ}
    (hc : 0 ≤ c) (hl : 0 ≤ low) (hh : 0 < high)
    (hvar : varG ≤ high * varH) (hp : c * varH ≤ energyH)
    (henergy : low * energyH ≤ energyG) :
    c * low / high * varG ≤ energyG := by
  calc
    c * low / high * varG ≤ c * low / high * (high * varH) :=
      mul_le_mul_of_nonneg_left hvar (div_nonneg (mul_nonneg hc hl) hh.le)
    _ = low * (c * varH) := by field_simp
    _ ≤ low * energyH := mul_le_mul_of_nonneg_left hp hl
    _ ≤ energyG := henergy

#print axioms shared_spoke_orders
#print axioms single_vacancy_variance
#print axioms vacancy_density_variance
#print axioms density_gap_floor

end ObserverPatchHolography.LoopVacuumBounds
