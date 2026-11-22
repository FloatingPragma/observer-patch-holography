import ObserverPatchHolography.ObservablePriority.Network
import ObserverPatchHolography.ObservablePriority.InstalledRows

namespace InstalledNetwork
open ObserverPatchSharp PriorityCellData CellBits PriorityNetwork

inductive Federation where
  | twenty | eighty | threeTwenty
  deriving DecidableEq, Repr

def size : Federation → Nat
  | .twenty => 20 | .eighty => 80 | .threeTwenty => 320

def rows (s : Federation) := rowAt (size s) InstalledRows.inventory
def board (s : Federation) := PriorityNetwork.inventory (size s) (rows s)
def chart (s : Federation) (k : Nat) := make (12*size s) (rows s k)

set_option maxRecDepth 100000
set_option maxHeartbeats 0

theorem rows_length (s : Federation) :
    (rowsFrom (size s) InstalledRows.inventory).length = size s := by cases s <;> decide

theorem admitted_rows (s : Federation) : Admitted (size s) (rows s) :=
  from_blocks (fun _ _ h _ hr => InstalledRows.every_row h hr) (rows_length s)

theorem row_numbers (s : Federation) : ∀ k < size s, (rows s k).cell = k := by
  have ok : (List.range (size s)).all (fun k => decide ((rows s k).cell = k)) = true := by
    cases s <;> decide
  exact fun k hk => of_decide_eq_true (List.all_eq_true.mp ok k (List.mem_range.mpr hk))

theorem covered_rows (s : Federation) : Covered (size s) (rows s) := by
  have ok : (List.range (size s)).all (ownerOK (size s) (rows s)) = true := by
    cases s <;> decide
  exact fun k hk => List.all_eq_true.mp ok k (List.mem_range.mpr hk)

theorem required_rows (s : Federation) : ∀ k < size s, boundaryRequired (chart s k) = true := by
  have ok : (List.range (size s)).all (fun k => boundaryRequired (chart s k)) = true := by
    cases s <;> decide
  exact fun k hk => List.all_eq_true.mp ok k (List.mem_range.mpr hk)

/-- Explicit legal complete schedule obtained from all addresses of this federation. -/
def candidate (s : Federation) : ObservablePriority.Candidate (board s) :=
  installed (covered_rows s) (locallyOrdered (admitted_rows s))

theorem rank_eq (s : Federation) : (candidate s).rank = cellRank := rfl

theorem full_window (s : Federation) (k : Nat) (hk : k < size s) :
    Ordered cellRank (candidate s).schedule ∧
    Ordered cellRank (menu (chart s k)) ∧
    PredecessorClosed cellRank (candidate s).schedule (keep k) ∧
    ((candidate s).schedule.filter (keep k)).Perm (menu (chart s k)) :=
  ⟨(candidate s).legal,locallyOrdered (admitted_rows s) k hk,closed _ _,
   restriction (covered_rows s) (locallyOrdered (admitted_rows s)) _ hk (candidate s).complete⟩

theorem minimal_window (s : Federation) (k : Nat) (hk : k < size s) (other : Edge → Bool)
    (cl : PredecessorClosed cellRank (candidate s).schedule other)
    (demand : ∀ a ∈ (chart s k).inner, other a = true) :
    ∀ b ∈ menu (chart s k), other b = true :=
  least_window k hk (required_rows s k hk) (candidate s).complete other cl demand

theorem every_full_cost (s : Federation) (k : Nat) (hk : k < size s) (full : List Edge)
    (p : full.Perm (board s)) (legal : Ordered cellRank full) :
    let c := chart s k
    Communicates c (fun L => localWord c (sharpRun cellRank full L)) 3 ∧
    ∀ b, Communicates c (fun L => localWord c (sharpRun cellRank full L)) b → 3 ≤ b :=
  physical_cost (covered_rows s) (admitted_rows s) k hk full p legal

theorem every_full_tree (s : Federation) (k : Nat) (hk : k < size s) (full : List Edge)
    (p : full.Perm (board s)) (legal : Ordered cellRank full) :
    let c := chart s k
    ((threeQuestions c cellRank).depth = 3 ∧ ∀ L,
      (threeQuestions c cellRank).run c L = localWord c (sharpRun cellRank full L)) ∧
    ∀ x : DecisionTree, (∀ L, x.run c L = localWord c (sharpRun cellRank full L)) → 3 ≤ x.depth :=
  physical_tree (covered_rows s) (admitted_rows s) k hk full p legal

theorem installed_cost (s : Federation) (k : Nat) (hk : k < size s) :
    let c := chart s k
    let full := (candidate s).schedule
    Communicates c (fun L => localWord c (sharpRun cellRank full L)) 3 ∧
    ∀ b, Communicates c (fun L => localWord c (sharpRun cellRank full L)) b → 3 ≤ b :=
  every_full_cost s k hk _ (candidate s).complete (candidate s).legal

#print axioms full_window
#print axioms minimal_window
#print axioms installed_cost
#print axioms every_full_cost
#print axioms every_full_tree
end InstalledNetwork
