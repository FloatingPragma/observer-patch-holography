import Std

/-!
A six-port word is represented by the bit mask of its double loads; all
other ports are vacancies. Each wing has degrees 6,6,5,5,5,5 from the ring
outward. The diagonal potential counts unlike neighbors and doubles at
private (degree five) sites. Energies below refer to the detached wing.
-/

namespace IntervalEcho.ReleasedWing

def occupied (word site : Nat) : Bool := word.testBit site

def inventory : List Nat :=
  (List.range 64).filter fun w => ((List.range 6).filter (occupied w)).length == 3

def walls (w : Nat) : Nat :=
  ((List.range 5).filter fun s => occupied w s != occupied w (s+1)).length

def privateDoubles (w : Nat) : Nat :=
  ([2,3,4,5].filter (occupied w)).length

def potential (w : Nat) : Nat := 32 * (walls w + privateDoubles w)

def twiceTrialEnergy (w : Nat) : Nat := 2 * potential w + walls w

theorem inventory_count : inventory.length = 20 := by decide +kernel

theorem contrast_histogram :
    ([64,96,128,160,192,224].map fun energy =>
      (inventory.filter fun w => potential w == energy).length) =
    [1,1,5,5,6,2] := by decide +kernel

theorem only_low_word : (inventory.filter fun w => potential w < 96) = [7] := by
  decide +kernel

theorem trial_energy :
    7 ∈ inventory ∧ potential 7 = 64 ∧ twiceTrialEnergy 7 = 129 := by decide +kernel

end IntervalEcho.ReleasedWing
