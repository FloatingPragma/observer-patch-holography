import Mathlib.Data.Fintype.Basic
import Mathlib.Data.Fin.Basic

namespace OPH.RepairOutputs

/-- The two possible loads at the first endpoint; the other endpoint uses the opposite coin. -/
def level (a b : Fin 3) (coin : Bool) : Nat :=
  (a.val + b.val + (if coin then 1 else 0)) / 2

def contrast (a b : Fin 3) : Int := (a.val : Int) - b.val

/-- When tracing the opposite endpoint and signed head, different input matrix
units can never interfere in a single reader's marginal. The coin is part
of the traced head and can be fixed on both sides of the matrix element. -/
theorem source_custody : ∀ a b u v : Fin 3, ∀ c : Bool,
    contrast a b = contrast u v → level a b c = level u v c → a = u ∧ b = v := by
  decide

/-- For one fixed value of the unobserved coin, the reduced matrix
coefficient is diagonal in both the incoming pair and the outgoing load.
Adding the two terms and dividing by two gives the entire endpoint channel,
including all off-diagonal input units and spectator matrix elements. -/
theorem local_coefficient (a b u v : Fin 3) (c : Bool) (i j : Nat) :
    (if contrast a b = contrast u v ∧ level a b (!c) = level u v (!c) ∧
           level a b c = i ∧ level u v c = j then 1 else 0 : Nat) =
      (if a = u ∧ b = v ∧ i = j ∧ level a b c = i then 1 else 0) := by
  apply if_congr _ rfl rfl
  constructor
  · rintro ⟨hd, ho, hi, hj⟩
    obtain ⟨rfl, rfl⟩ := source_custody a b u v (!c) hd ho
    exact ⟨rfl, rfl, hi.symm.trans hj, hi⟩
  · rintro ⟨rfl, rfl, rfl, rfl⟩
    exact ⟨rfl, rfl, rfl, rfl⟩

theorem equal_iff_even : ∀ a b : Fin 3, ∀ c : Bool,
    level a b c = level a b (!c) ↔ (a.val + b.val) % 2 = 0 := by
  decide

theorem recover_sum : ∀ a b : Fin 3, ∀ c : Bool,
    level a b c + level a b (!c) = a.val + b.val := by
  decide

/-- An all-event agreement test on both odd sum levels forces every local
coarsening of the three levels to take a single value. -/
theorem full_level_coarsening {Y : Type*} (f : Nat → Y)
    (agree : ∀ a b : Fin 3, f (level a b false) = f (level a b true)) :
    f 0 = f 1 ∧ f 1 = f 2 :=
  ⟨agree 0 1, agree 1 2⟩

/-- The sharp question on an enrolled set is independently readable exactly
when any two source labels sharing a possible output have the same answer.
This also characterizes functions constant on the generated components. -/
theorem exact_decoder_iff {X O Y : Type*} (defaultAnswer : Y)
    (reachable : X → O → Prop) (f : X → Y) :
    (∃ d : O → Y, ∀ x o, reachable x o → d o = f x) ↔
    ∀ x y o, reachable x o → reachable y o → f x = f y := by
  classical
  constructor
  · rintro ⟨d, h⟩ x y o hx hy
    exact (h x o hx).symm.trans (h y o hy)
  · intro h
    refine ⟨fun o => if ex : ∃ x, reachable x o then f (Classical.choose ex)
                              else defaultAnswer, ?_⟩
    intro x o hx
    have ex : ∃ y, reachable y o := ⟨x, hx⟩
    simp only [dif_pos ex]
    exact h _ x o (Classical.choose_spec ex) hx

end OPH.RepairOutputs
