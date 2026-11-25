import EquationOfState
import Mathlib.Tactic

/-! A held dimer accumulates a relative phase per work cycle. The ruler reads
three return times proportional to cycle duration. All energies here are
relative to the ground mode in the same clock units. -/
namespace OPH.MeteredEnergyCone

theorem fixed_action_gap (action t ref : ℝ) (ha : action ≠ 0) (ht : t ≠ 0)
    (hr : ref ≠ 0) : (action/t)/(action/ref) = ref/t := by field_simp

theorem fixed_duty_gap (action t ref : ℝ) (ha : action ≠ 0) (ht : t ≠ 0)
    (hr : ref ≠ 0) : ((action*t/ref)/t)/(action/ref) = 1 := by field_simp

/-- Eulerian pressure on the supplied deformation, with occupation held fixed.
The numerator differentiates the actually inferred gap times its occupation,
and the denominator differentiates the radar volume. -/
theorem work_pressure (p l : ℝ) (hl : 0 < l) :
    -(deriv (fun x : ℝ => p/x) l) / deriv (fun x : ℝ => x^3) l
    = (p/l)/(3*l^3) := by
  have h := (hasDerivAt_const l p).div (hasDerivAt_id l) (ne_of_gt hl)
  have v := (hasDerivAt_id l).pow 3
  change HasDerivAt (fun x : ℝ => p/x) _ l at h
  change HasDerivAt (fun x : ℝ => x^3) _ l at v
  rw [h.deriv, v.deriv]
  simp only [id_eq]
  field_simp
  ring

theorem constant_duty_pressure (p l : ℝ) :
    -(deriv (fun _ : ℝ => p) l) / deriv (fun x : ℝ => x^3) l = 0 := by simp

theorem compression_work (p : ℝ) : p/1 - p/2 = p/2 := by ring

section Reach
variable {X : Type*} [PseudoMetricSpace X]

/-- A path of n energy exchanges of length at most one has reach at most n
in the ruler's own metric. -/
theorem bounded_reach (walk : ℕ → X) (h : ∀ i, dist (walk i) (walk (i+1)) ≤ 1)
    (n : ℕ) : dist (walk 0) (walk n) ≤ (n : ℝ) := by
  induction n with
  | zero => simp
  | succ n ih =>
    calc
      dist (walk 0) (walk (n+1)) ≤ dist (walk 0) (walk n) +
        dist (walk n) (walk (n+1)) := dist_triangle _ _ _
      _ ≤ (n : ℝ) + 1 := add_le_add ih (h n)
      _ = _ := by norm_num

theorem disjoint_inputs (a b : X) (r s : ℝ) (h : r+s < dist a b) :
    ¬ ∃ z : X, dist a z ≤ r ∧ dist b z ≤ s := by
  rintro ⟨z, ha, hb⟩
  have t := dist_triangle a z b
  rw [dist_comm z b] at t
  linarith
end Reach

#print axioms fixed_action_gap
#print axioms fixed_duty_gap
#print axioms work_pressure
#print axioms constant_duty_pressure
#print axioms compression_work
#print axioms bounded_reach
#print axioms disjoint_inputs
end OPH.MeteredEnergyCone
