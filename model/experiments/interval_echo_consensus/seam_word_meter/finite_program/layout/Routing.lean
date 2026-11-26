import Mathlib.Algebra.Group.End
import Mathlib.Tactic

namespace OPH.ShelfRouting

variable {V : Type*} [DecidableEq V]

/-- Swaps along a path in order, and back through every edge except the last. -/
def routeWord (a z : V) : List V → Equiv.Perm V
  | [] => Equiv.swap a z
  | b :: bs => Equiv.swap a b * routeWord b z bs * Equiv.swap a b

/-- The native exchange word permutes just the required endpoint labels.
The destination is absent from the interior of a simple physical route. -/
theorem routeWord_endpoints (via : List V) (a z : V) (ha : z ≠ a) (hn : z ∉ via) :
    routeWord a z via = Equiv.swap a z := by
  induction via generalizing a with
  | nil => rfl
  | cons b bs ih =>
    have hzb : z ≠ b := fun eq => hn (List.mem_cons.mpr (Or.inl eq))
    have hn' : z ∉ bs := fun h => hn (List.mem_cons_of_mem b h)
    simp only [routeWord]
    rw [ih b hzb hn']
    have h := Equiv.swap_mul_swap_mul_swap (x := z) (y := b) (z := a) hzb ha
    rw [Equiv.swap_comm b a, Equiv.swap_comm z b] at h
    exact h

/-- One writer, two publication dimers and two three-port stations reserve
twelve separate ports. Each data digit consumes two more. -/
theorem storage_iff (owners bits : ℕ) (h : 2 ≤ owners) :
    (2*bits + 12 ≤ 12*owners) ↔ bits ≤ 6*owners - 6 := by omega

/-- The two mounted station cells each leave four dimers. Every ordinary
cell leaves six, and two of these nonstation dimers publish flags. -/
theorem construction_saturates (owners : ℕ) (h : 2 ≤ owners) :
    4+4+6*(owners-2)-2 = 6*owners-6 := by omega

#print axioms routeWord_endpoints
#print axioms storage_iff

end OPH.ShelfRouting
