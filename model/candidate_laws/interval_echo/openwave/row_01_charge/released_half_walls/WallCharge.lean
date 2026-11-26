import Mathlib.Tactic

namespace IntervalEcho.WallCharge

def bit (a : Bool) : ℤ := if a then 1 else 0

/-- The two signed equal adjacent readings detect a reversal of alternating
order. Opposite readings give its unperturbed continuation. -/
def bond : Bool → Bool → ℤ
  | true, true => 1
  | false, false => -1
  | _, _ => 0

theorem two_site (a b : Bool) : bond a b = bit a + bit b - 1 := by
  cases a <;> cases b <;> decide

/-- Twice the total polarization charge in occupied-site units. The successor
permutation accounts for every site once at each end of its two bonds. -/
theorem ring_charge {V : Type*} [Fintype V] (succ : V ≃ V) (s : V → Bool) :
    ∑ i, bond (s i) (s (succ i)) = 2 * (∑ i, bit (s i)) - Fintype.card V := by
  simp only [two_site, Finset.sum_sub_distrib, Finset.sum_add_distrib]
  have moved : ∑ i, bit (s (succ i)) = ∑ i, bit (s i) :=
    Equiv.sum_comp succ (fun i => bit (s i))
  rw [moved]
  simp
  ring

/-- A real field linear in physical occupations can be implemented as the
difference of two native quadratic contrasts against zero and two. -/
theorem linear_voltage (load : ℝ) : load^2 - (load-2)^2 = 4*(load-1) := by ring

/-- For a sharp addition, let m count minority equal bonds. Energy conservation
gives 2*c*E[m] = E[off] and the sum of off-diagonal norms is at most n/2.
The probability of any minority bond is at most E[m]. -/
theorem two_wall_retention (c n m off bad : ℝ) (hc : 0 < c)
    (energy : 2*c*m = off) (exchange : off ≤ n/2) (markov : bad ≤ m) :
    1 - n/(4*c) ≤ 1-bad := by
  have h : m ≤ n/(4*c) := by
    apply (le_div_iff₀ (by positivity : (0:ℝ) < 4*c)).2
    linarith
  linarith

end IntervalEcho.WallCharge
