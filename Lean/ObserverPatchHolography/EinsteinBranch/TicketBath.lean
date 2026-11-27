import Mathlib.Data.Complex.Basic
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Data.Fintype.Fin
import Mathlib.Data.Fintype.Prod
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring
import Mathlib.Util.AssertNoSorry

/-! A complete local bath sweep conditioned on its three outgoing bits. -/

namespace OPH.EinsteinBranch.TicketBath

abbrev Ticket := Bool × Bool × Bool
abbrev Block (ι : Type*) := Ticket → ι → Ticket → ι → ℂ

/-- Read the displaced bit and insert an independent fair bit in its place.
The uninterpreted indices may be any outside quantum system. -/
noncomputable def first {ι : Type*} (z : Bool) (ρ : Block ι) : Block ι :=
  fun x i y j => if x.1 = y.1 then ρ (z, x.2) i (z, y.2) j / 2 else 0

noncomputable def second {ι : Type*} (z : Bool) (ρ : Block ι) : Block ι :=
  fun x i y j => if x.2.1 = y.2.1 then
    ρ (x.1, z, x.2.2) i (y.1, z, y.2.2) j / 2 else 0

noncomputable def third {ι : Type*} (z : Bool) (ρ : Block ι) : Block ι :=
  fun x i y j => if x.2.2 = y.2.2 then
    ρ (x.1, x.2.1, z) i (y.1, y.2.1, z) j / 2 else 0

noncomputable def sweep {ι : Type*} (z : Ticket) (ρ : Block ι) : Block ι :=
  third z.2.2 (second z.2.1 (first z.1 ρ))

/-- Every conditional output is the uniform ticket times the full surviving
outside matrix block, including off-diagonal spectator entries. -/
theorem conditional_sweep {ι : Type*} (ρ : Block ι) (z x y : Ticket) (i j : ι) :
    sweep z ρ x i y j = if x = y then ρ z i z j / 8 else 0 := by
  rcases x with ⟨a, b, c⟩
  rcases y with ⟨a', b', c'⟩
  rcases z with ⟨u, v, w⟩
  by_cases ha : a = a' <;> by_cases hb : b = b' <;> by_cases hc : c = c' <;>
    simp [sweep, first, second, third, ha, hb, hc, div_div]
  ring

/-- The three less populous words each have a unique ticket. -/
def word : Ticket → Fin 4
  | (true, false, true) => 1
  | (true, true, false) => 2
  | (true, true, true) => 3
  | _ => 0

def ticketCount (w : Fin 4) : ℕ := if w = 0 then 5 else 1

theorem word_card (w : Fin 4) :
    (Finset.univ.filter fun x : Ticket => word x = w).card = ticketCount w := by
  fin_cases w <;> decide

/-- Discard the ticket after writing its native word. Conditioning on any bath
log leaves that word's prescribed probability times the same exterior block. -/
theorem selected_word_block {ι : Type*} (ρ : Block ι) (z : Ticket) (w : Fin 4) (i j : ι) :
    (∑ x ∈ Finset.univ.filter (fun x : Ticket => word x = w), sweep z ρ x i x j)
      = (ticketCount w : ℂ) * (ρ z i z j / 8) := by
  simp [conditional_sweep, word_card, nsmul_eq_mul]

assert_no_sorry conditional_sweep
assert_no_sorry selected_word_block
#print axioms conditional_sweep
#print axioms selected_word_block

end OPH.EinsteinBranch.TicketBath
