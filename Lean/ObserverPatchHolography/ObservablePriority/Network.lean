import ObserverPatchHolography.ObservablePriority.PriorityInstallation
import ObserverPatchHolography.ObservablePriority.PriorityCellData
import Mathlib.Data.List.Dedup

namespace PriorityNetwork
open ObserverPatchSharp ObservablePriority PriorityCellData CellBits

abbrev Table := Nat → Row
def localMenu (N : Nat) (rows : Table) (k : Nat) := menu (make (12*N) (rows k))
def keep (k : Nat) (e : Edge) := decide (e.lo/12 = k ∨ e.hi/12 = k)
def inventory (N : Nat) (rows : Table) :=
  ((List.range N).flatMap (localMenu N rows)).dedup

theorem mem_inventory {N : Nat} {rows : Table} {e : Edge} :
    e ∈ inventory N rows ↔ ∃ k < N, e ∈ localMenu N rows k := by
  simp [inventory]

/-- Ports not in the row's own cell must name the identical addressed edge
    on the opposite cell's menu. Both endpoints must belong to this screen. -/
def ownerOK (N : Nat) (rows : Table) (k : Nat) : Bool :=
  (localMenu N rows k).all fun e =>
    keep k e && decide (e.lo/12 < N) && decide (e.hi/12 < N) &&
    (decide (e.lo/12 = k) || (localMenu N rows (e.lo/12)).contains e) &&
    (decide (e.hi/12 = k) || (localMenu N rows (e.hi/12)).contains e)

def Covered (N : Nat) (rows : Table) := ∀ k < N, ownerOK N rows k = true
def LocallyOrdered (N : Nat) (rows : Table) :=
  ∀ k < N, Ordered cellRank (localMenu N rows k)

private theorem owners {N : Nat} {r : Table} (h : Covered N r) {e : Edge}
    (he : e ∈ inventory N r) :
    (e.lo/12 < N ∧ e ∈ localMenu N r (e.lo/12)) ∧
    (e.hi/12 < N ∧ e ∈ localMenu N r (e.hi/12)) := by
  obtain ⟨k,hk,hl⟩ := mem_inventory.mp he
  have ok := List.all_eq_true.mp (h k hk) e hl
  simp only [Bool.and_eq_true,Bool.or_eq_true,decide_eq_true_eq,List.contains_iff_mem] at ok
  constructor
  · exact ⟨ok.1.1.1.2,ok.1.2.elim (fun eq => eq ▸ hl) id⟩
  · exact ⟨ok.1.1.2,ok.2.elim (fun eq => eq ▸ hl) id⟩

private theorem no_local_tie {xs : List Edge} {a b : Edge} (h : Ordered cellRank xs)
    (ha : a ∈ xs) (hb : b ∈ xs) (ab : a ≠ b) (meet : meets a b = true) :
    cellRank a ≠ cellRank b := by
  induction xs with
  | nil => simp at ha
  | cons hd tl ih =>
    obtain ⟨hhd,htl⟩ := List.pairwise_cons.mp h
    rcases List.mem_cons.mp ha with rfl | hat
    · exact Nat.ne_of_gt ((hhd b ((List.mem_cons.mp hb).resolve_left (Ne.symm ab))).2 meet)
    rcases List.mem_cons.mp hb with rfl | hbt
    · exact Nat.ne_of_lt ((hhd a hat).2 ((meets_symm ..).trans meet))
    · exact ih htl hat hbt

theorem ties {N : Nat} {r : Table} (h : Covered N r) (legal : LocallyOrdered N r) :
    noIncidentTies (inventory N r) cellRank = true := by
  apply List.all_eq_true.mpr
  intro a ha
  apply List.all_eq_true.mpr
  intro b hb
  apply decide_eq_true
  by_cases hab : a = b
  · exact .inl hab
  by_cases meet : meets a b = true
  · right; right
    obtain ⟨al,ar⟩ := owners h ha
    obtain ⟨bl,br⟩ := owners h hb
    have common := meet
    simp only [meets,incident,Bool.or_eq_true,decide_eq_true_eq] at common
    rcases common with (eq|eq)|(eq|eq)
    · exact no_local_tie (legal _ al.1) al.2 (eq ▸ bl.2) hab meet
    · exact no_local_tie (legal _ al.1) al.2 (eq ▸ br.2) hab meet
    · exact no_local_tie (legal _ ar.1) ar.2 (eq ▸ bl.2) hab meet
    · exact no_local_tie (legal _ ar.1) ar.2 (eq ▸ br.2) hab meet
  · exact .inr (.inl (Bool.eq_false_iff.mpr meet))

def installed {N : Nat} {r : Table} (h : Covered N r) (legal : LocallyOrdered N r) :
    Candidate (inventory N r) :=
  ofRank _ cellRank (List.nodup_dedup _) (ties h legal)

theorem installed_rank {N : Nat} {r : Table} (h : Covered N r) (l : LocallyOrdered N r) :
    (installed h l).rank = cellRank := rfl

theorem closed (k : Nat) (full : List Edge) : PredecessorClosed cellRank full (keep k) := by
  intro a _ hak b _ hm hr
  by_cases same : a.lo/12 = a.hi/12
  · simp [keep] at hak ⊢
    simp [meets,incident] at hm
    rcases hm with (eq|eq)|(eq|eq) <;> rcases hak with hk|hk <;> grind
  · simp [cellRank,same] at hr

theorem restriction {N : Nat} {r : Table} (cov : Covered N r) (legal : LocallyOrdered N r)
    (k : Nat) (hk : k < N) {full : List Edge} (p : full.Perm (inventory N r)) :
    (full.filter (keep k)).Perm (localMenu N r k) := by
  have d1 : (full.filter (keep k)).Nodup := List.Pairwise.filter _ (p.nodup_iff.mpr (List.nodup_dedup _))
  have d2 : (localMenu N r k).Nodup := (legal k hk).imp (fun h => h.1)
  have same : ∀ e, e ∈ full.filter (keep k) ↔ e ∈ localMenu N r k := by
    intro e
    simp only [List.mem_filter]
    constructor
    · rintro ⟨hm,hh⟩
      obtain ⟨hl,hr⟩ := owners cov (p.mem_iff.mp hm)
      rcases of_decide_eq_true hh with hl'|hr'
      · exact hl' ▸ hl.2
      · exact hr' ▸ hr.2
    · intro hm
      refine ⟨p.mem_iff.mpr (mem_inventory.mpr ⟨k,hk,hm⟩),?_⟩
      have ok := List.all_eq_true.mp (cov k hk) e hm
      exact (Bool.and_eq_true_iff.mp (Bool.and_eq_true_iff.mp
        (Bool.and_eq_true_iff.mp (Bool.and_eq_true_iff.mp ok).1).1).1).1
  apply List.perm_iff_count.mpr
  intro e
  simp only [d1.count,d2.count,same]

#print axioms ties
#print axioms installed
#print axioms closed
#print axioms restriction

abbrev Blocks := List (Nat × List Row)
def rowsFrom (N : Nat) (blocks : Blocks) : List Row :=
  (blocks.filter (fun entry => entry.1 == 12*N)).flatMap Prod.snd

def emptyRow : Row := ⟨0,⟨0,0,0⟩,⟨0,0,0⟩,⟨0,0,0⟩⟩
def rowAt (N : Nat) (blocks : Blocks) (k : Nat) : Row :=
  (rowsFrom N blocks)[k]?.getD emptyRow

theorem from_blocks {blocks : Blocks}
    (all : ∀ n rs, (n,rs) ∈ blocks → ∀ r ∈ rs, admitted n r = true)
    {N : Nat} (len : (rowsFrom N blocks).length = N) (k : Nat) (hk : k < N) :
    admitted (12*N) (rowAt N blocks k) = true := by
  have hlt : k < (rowsFrom N blocks).length := by rw [len]; exact hk
  have mem : rowAt N blocks k ∈ rowsFrom N blocks := by
    simp [rowAt,List.getElem?_eq_getElem hlt,List.getElem_mem]
  obtain ⟨⟨n,rs⟩,hm,hr⟩ := List.mem_flatMap.mp mem
  obtain ⟨hb,hn⟩ := List.mem_filter.mp hm
  have eq : n = 12*N := by simpa using hn
  exact eq ▸ all n rs hb _ hr

def Admitted (N : Nat) (t : Table) := ∀ k < N, admitted (12*N) (t k) = true

theorem locallyOrdered {N : Nat} {t : Table} (ad : Admitted N t) : LocallyOrdered N t :=
  fun k hk => (legalB_spec _ _).mp (Bool.and_eq_true_iff.mp (ad k hk)).2

theorem physical_cost {N : Nat} {t : Table} (cov : Covered N t) (ad : Admitted N t)
    (k : Nat) (hk : k < N) (full : List Edge) (p : full.Perm (inventory N t))
    (legal : Ordered cellRank full) :
    let c := make (12*N) (t k)
    Communicates c (fun L => localWord c (sharpRun cellRank full L)) 3 ∧
    ∀ b, Communicates c (fun L => localWord c (sharpRun cellRank full L)) b → 3 ≤ b :=
  exactly_three_full _ _ full (keep k) legal (locallyOrdered ad k hk) (closed k full)
    (restriction cov (locallyOrdered ad) k hk p) (Bool.and_eq_true_iff.mp (ad k hk)).1

theorem physical_tree {N : Nat} {t : Table} (cov : Covered N t) (ad : Admitted N t)
    (k : Nat) (hk : k < N) (full : List Edge) (p : full.Perm (inventory N t))
    (legal : Ordered cellRank full) :
    let c := make (12*N) (t k)
    ((threeQuestions c cellRank).depth = 3 ∧ ∀ L,
      (threeQuestions c cellRank).run c L = localWord c (sharpRun cellRank full L)) ∧
    ∀ x : DecisionTree, (∀ L, x.run c L = localWord c (sharpRun cellRank full L)) → 3 ≤ x.depth :=
  exactly_three_tree_full _ _ full (keep k) legal (locallyOrdered ad k hk) (closed k full)
    (restriction cov (locallyOrdered ad) k hk p) (Bool.and_eq_true_iff.mp (ad k hk)).1

def boundaryRequired (c : CellChart) : Bool :=
  c.outer.toList.all fun b => c.inner.any fun a => meets a b && decide (cellRank b < cellRank a)

/-- Every closed patch containing the thirty internal demands contains all
    three named boundary addresses, even if another subset was proposed. -/
theorem least_window {N : Nat} {t : Table} (k : Nat) (hk : k < N)
    (required : boundaryRequired (make (12*N) (t k)) = true)
    {full : List Edge} (pf : full.Perm (inventory N t)) (other : Edge → Bool)
    (cl : PredecessorClosed cellRank full other)
    (demand : ∀ a ∈ (make (12*N) (t k)).inner, other a = true) :
    ∀ b ∈ localMenu N t k, other b = true := by
  intro b hb
  have inside : ∀ e ∈ localMenu N t k, e ∈ full :=
    fun e he => pf.mem_iff.mpr (mem_inventory.mpr ⟨k,hk,he⟩)
  rcases List.mem_append.mp hb with hi|ho
  · exact demand b hi
  · obtain ⟨a,ha,test⟩ := List.any_eq_true.mp (List.all_eq_true.mp required b ho)
    obtain ⟨meets,lt⟩ := Bool.and_eq_true_iff.mp test
    exact cl a (inside a (List.mem_append_left _ ha)) (demand a ha) b (inside b hb) meets
      (of_decide_eq_true lt)

#print axioms physical_cost
#print axioms physical_tree
#print axioms least_window
end PriorityNetwork
