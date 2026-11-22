import ObserverPatchHolography.ObservablePriority.Screen20Block0
import ObserverPatchHolography.ObservablePriority.Screen80Block0
import ObserverPatchHolography.ObservablePriority.Screen80Block1
import ObserverPatchHolography.ObservablePriority.Screen80Block2
import ObserverPatchHolography.ObservablePriority.Screen80Block3
import ObserverPatchHolography.ObservablePriority.Screen320Block0
import ObserverPatchHolography.ObservablePriority.Screen320Block1
import ObserverPatchHolography.ObservablePriority.Screen320Block2
import ObserverPatchHolography.ObservablePriority.Screen320Block3
import ObserverPatchHolography.ObservablePriority.Screen320Block4
import ObserverPatchHolography.ObservablePriority.Screen320Block5
import ObserverPatchHolography.ObservablePriority.Screen320Block6
import ObserverPatchHolography.ObservablePriority.Screen320Block7
import ObserverPatchHolography.ObservablePriority.Screen320Block8
import ObserverPatchHolography.ObservablePriority.Screen320Block9
import ObserverPatchHolography.ObservablePriority.Screen320Block10
import ObserverPatchHolography.ObservablePriority.Screen320Block11
import ObserverPatchHolography.ObservablePriority.Screen320Block12
import ObserverPatchHolography.ObservablePriority.Screen320Block13
import ObserverPatchHolography.ObservablePriority.Screen320Block14
import ObserverPatchHolography.ObservablePriority.Screen320Block15
namespace InstalledRows
open PriorityCellData CellBits ObserverPatchSharp
def inventory : List (Nat × List Row) := [(240,Screen20Block0.rows),
  (960,Screen80Block0.rows),
  (960,Screen80Block1.rows),
  (960,Screen80Block2.rows),
  (960,Screen80Block3.rows),
  (3840,Screen320Block0.rows),
  (3840,Screen320Block1.rows),
  (3840,Screen320Block2.rows),
  (3840,Screen320Block3.rows),
  (3840,Screen320Block4.rows),
  (3840,Screen320Block5.rows),
  (3840,Screen320Block6.rows),
  (3840,Screen320Block7.rows),
  (3840,Screen320Block8.rows),
  (3840,Screen320Block9.rows),
  (3840,Screen320Block10.rows),
  (3840,Screen320Block11.rows),
  (3840,Screen320Block12.rows),
  (3840,Screen320Block13.rows),
  (3840,Screen320Block14.rows),
  (3840,Screen320Block15.rows)]
def valid (item : Nat × List Row) := item.2.all (admitted item.1)
theorem verified : inventory.all valid = true := by
  simp only [inventory,valid,List.all_cons,List.all_nil,Screen20Block0.verified,Screen80Block0.verified,Screen80Block1.verified,Screen80Block2.verified,Screen80Block3.verified,Screen320Block0.verified,Screen320Block1.verified,Screen320Block2.verified,Screen320Block3.verified,Screen320Block4.verified,Screen320Block5.verified,Screen320Block6.verified,Screen320Block7.verified,Screen320Block8.verified,Screen320Block9.verified,Screen320Block10.verified,Screen320Block11.verified,Screen320Block12.verified,Screen320Block13.verified,Screen320Block14.verified,Screen320Block15.verified,Bool.and_self]
theorem every_row {n : Nat} {rows : List Row} (h : (n,rows) ∈ inventory)
    {r : Row} (hr : r ∈ rows) : admitted n r = true :=
  List.all_eq_true.mp (List.all_eq_true.mp verified _ h) _ hr
theorem costs {n : Nat} {rows : List Row} (h : (n,rows) ∈ inventory)
    {r : Row} (hr : r ∈ rows) :
    let c := make n r
    Communicates c (fun L => localWord c (sharpRun cellRank (menu c) L)) 3 ∧
    ∀ b, Communicates c (fun L => localWord c (sharpRun cellRank (menu c) L)) b → 3 ≤ b :=
  row_exact n r (every_row h hr)
#print axioms verified
#print axioms costs
end InstalledRows
