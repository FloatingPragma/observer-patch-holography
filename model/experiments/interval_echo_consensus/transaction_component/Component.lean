import ObserverPatchHolography.Provenance.HistoryCausalInvariance
import Mathlib.Tactic.FinCases

namespace OPH.TransactionComponent

open OPH.Provenance

/-- Two collars, with boundary, three data bits, and seal on each. -/
abbrev Reg := Bool × Fin 5
abbrev Id := Fin 4
abbrev File := VersionedState Reg Bool Id
abbrev Spec := CommitSpec Reg Bool Id

def island (e : Id) : Bool := e.val % 2 == 1

def reads (e : Id) : Finset Reg := Finset.univ.filter (fun r => r.1 == island e)

def writes (e : Id) : Finset Reg :=
  (reads e).filter (fun r => if e.val < 2 then r.2.val > 0 ∧ r.2.val < 4 else r.2.val = 4)

def spec (boundary : Bool → Bool) (e : Id) : Spec where
  eventId := e
  readSet := reads e
  causalSupp := reads e
  writeSet := writes e
  supp_subset_read := Finset.Subset.rfl
  writeValue := fun _ => if e.val < 2 then boundary (island e) else true

def start (value : Reg → Bool) : File := ⟨value, fun _ => none⟩

def order (d s : Bool) : List Id :=
  (if d then [1,0] else [0,1]) ++ (if s then [3,2] else [2,3])

def program (b : Bool → Bool) (d s : Bool) : List Spec := (order d s).map (spec b)

theorem fresh (v : Reg → Bool) (b : Bool → Bool) (d s : Bool) :
    FreshStart (start v) (program b d s) := by
  intro _ _ _
  simp [start]

theorem labels (b : Bool → Bool) (d s : Bool) :
    (program b d s).map CommitSpec.eventId = order d s := by
  cases d <;> cases s <;> rfl

theorem unique_labels (b : Bool → Bool) (d s : Bool) :
    ((program b d s).map CommitSpec.eventId).Nodup := by
  rw [labels]
  cases d <;> cases s <;> decide

theorem covers (b : Bool → Bool) (d s : Bool) :
    ExecutedCarrierComplete (program b d s) := by
  intro e
  rw [labels]
  cases d <;> cases s <;> fin_cases e <;> decide

/-- The four physically decoded labels therefore instantiate the actual provenance log. -/
noncomputable def eventLog (v : Reg → Bool) (b : Bool → Bool) (d s : Bool) :
    SemanticEventLog Reg Bool Id :=
  semanticEventLogOfExecution (start v) (program b d s) (fresh v b d s)
    (unique_labels b d s) (covers b d s)

theorem support_visible (b : Bool → Bool) (e : Id) (r : Reg) :
    r ∈ (spec b e).causalSupp ↔ r.1 = island e := by simp [spec, reads]

theorem disjoint_writes (b : Bool → Bool) (e f : Id) (h : island e ≠ island f) :
    IndependentPair (spec b e) (spec b f) := by
  have hr {i : Id} {r : Reg} (hw : r ∈ writes i) := (Finset.mem_filter.mp hw).1
  have sep {i j : Id} (ne : island i ≠ island j) {r : Reg}
      (a : r ∈ reads i) (c : r ∈ reads j) : False := by
    simp [reads] at a c
    exact ne (a.symm.trans c)
  change (∀ r ∈ writes e, r ∉ writes f) ∧
    (∀ r ∈ writes e, r ∉ reads f) ∧ (∀ r ∈ writes f, r ∉ reads e)
  exact ⟨fun _ a c => sep h (hr a) (hr c), fun _ a c => sep h (hr a) c,
    fun _ a c => sep (Ne.symm h) (hr a) c⟩

theorem disjoint_commit_diamond (q : File) (b : Bool → Bool) (e f : Id)
    (h : island e ≠ island f) :
    (spec b e).apply ((spec b f).apply q) = (spec b f).apply ((spec b e).apply q) := by
  exact (CommitSpec.apply_comm (disjoint_writes b e f h).1 q).symm

/-- Each seal reads precisely a material data register and version written by its publication. -/
theorem seal_cites_data (q : File) (b : Bool → Bool) (i : Bool) :
    let d : Id := if i then 1 else 0
    let s : Id := if i then 3 else 2
    AuthenticatedDirectSemanticParent ((spec b d).commitAt q)
      ((spec b s).commitAt ((spec b d).apply q)) := by
  have witness (d s : Id) (a : Reg) (hr : a ∈ reads s) (hw : a ∈ writes d) :
      AuthenticatedDirectSemanticParent ((spec b d).commitAt q)
        ((spec b s).commitAt ((spec b d).apply q)) := by
    refine ⟨a, hr, hw, ?_, rfl⟩
    exact CommitSpec.apply_writer_mem (s := spec b d) (q := q) hw
  cases i
  · exact witness 0 2 (false, 2) (by decide) (by decide)
  · exact witness 1 3 (true, 2) (by decide) (by decide)

/-- The local decoder evaluates all three bits in a single selected bank. -/
def visible (old shadow : Fin 3 → Bool) (published : Bool) := if published then shadow else old

theorem staging_is_private (old left right : Fin 3 → Bool) :
    visible old left false = visible old right false := rfl

theorem one_selector_publishes (old : Fin 3 → Bool) (boundary : Bool) :
    visible old (fun _ => boundary) true = fun _ => boundary := rfl

#print axioms disjoint_commit_diamond
#print axioms seal_cites_data
#print axioms eventLog
end OPH.TransactionComponent
