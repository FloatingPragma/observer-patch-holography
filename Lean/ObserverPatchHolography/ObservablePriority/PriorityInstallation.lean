import ObserverPatchHolography.ObservablePriority.ObservablePriority

namespace ObservablePriority
open ObserverPatchSharp

def noIncidentTies (screen : List Edge) (rank : Edge → Nat) : Bool :=
  screen.all fun a => screen.all fun b =>
    decide (a = b ∨ meets a b = false ∨ rank a ≠ rank b)

/-- Sorting constructs a complete schedule for every loop-free simple screen
    and numeric priority with no tie at an incident pair. -/
def ofRank (screen : List Edge) (rank : Edge → Nat) (unique : screen.Nodup)
    (separate : noIncidentTies screen rank = true) : Candidate screen := by
  let cmp := fun a b : Edge => decide (rank b ≤ rank a)
  let xs := screen.mergeSort cmp
  have perm : xs.Perm screen := List.mergeSort_perm screen cmp
  have monotone : xs.Pairwise (fun a b => cmp a b = true) :=
    List.pairwise_mergeSort (le := cmp) (fun a b c hab hbc =>
      decide_eq_true (Nat.le_trans (of_decide_eq_true hbc) (of_decide_eq_true hab)))
      (fun a b => by simp [cmp, Nat.le_total]) screen
  exact ⟨rank, xs, perm, (perm.nodup_iff.mpr unique).and monotone |>.imp_of_mem
    (fun {a b} ha hb ⟨hne, le⟩ => ⟨hne, fun hm => by
      have tested := List.all_eq_true.mp (List.all_eq_true.mp separate a (perm.mem_iff.mp ha)) b
        (perm.mem_iff.mp hb)
      have diff : rank a ≠ rank b := by simpa [hne, hm] using tested
      exact Nat.lt_of_le_of_ne (of_decide_eq_true le) (Ne.symm diff)⟩)⟩

/-- Install only after checking the two finite structural conditions. -/
def install? (screen : List Edge) (rank : Edge → Nat) : Option (Candidate screen) :=
  if unique : screen.Nodup then
    if separate : noIncidentTies screen rank then some (ofRank screen rank unique separate)
    else none
  else none

theorem install_iff (screen : List Edge) (rank : Edge → Nat) :
    (install? screen rank).isSome = true ↔ screen.Nodup ∧ noIncidentTies screen rank = true := by
  by_cases hu : screen.Nodup <;> cases ht : noIncidentTies screen rank <;> simp [install?, hu, ht]

#print axioms ofRank
#print axioms install_iff
end ObservablePriority
