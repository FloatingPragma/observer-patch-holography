import ObserverPatchHolography.ObservablePriority.BoundaryInput

/-! Exact fixed-length, zero-error classical communication for a sharp cell.
    The closed menu runs from tail to head. Each witness is an explicit physical
    load, and the checker computes its internal occupancy word by strict repair. -/
namespace CellBits

open ObserverPatchSharp PatchBoundaryInput

structure Three (α : Type) where
  a : α
  b : α
  c : α
  deriving DecidableEq, Repr

def Three.toList {α : Type} (t : Three α) : List α := [t.a,t.b,t.c]
abbrev Flags := Three Bool

def flags (i : Fin 8) : Flags := ⟨i.val.testBit 0, i.val.testBit 1, i.val.testBit 2⟩

def pack (f : Flags) : Fin 8 :=
  ⟨(if f.a then 1 else 0) + (if f.b then 2 else 0) + (if f.c then 4 else 0), by
    rcases f with ⟨a,b,c⟩; cases a <;> cases b <;> cases c <;> decide⟩

theorem flags_pack (f : Flags) : flags (pack f) = f := by
  rcases f with ⟨a,b,c⟩; cases a <;> cases b <;> cases c <;> decide

theorem flags_inj : ∀ i j : Fin 8, flags i = flags j → i = j := by decide

structure CellChart where
  nPorts : Nat
  own : List Nat
  inner : List Edge
  outer : Three Edge
  boundary : Three Nat
  remote : Three Nat
  spare : Three Nat
  deriving Repr

def menu (c : CellChart) := c.inner ++ c.outer.toList
def support (c : CellChart) := c.remote.toList ++ c.spare.toList

def localView (c : CellChart) (L : Source) : Source :=
  fun p => if p ∈ c.own then L p else .one

def censusOn (ports : List Nat) (L : Source) : Three Nat :=
  ⟨ports.countP (fun p => decide (L p = .zero)), ports.countP (fun p => decide (L p = .one)),
    ports.countP (fun p => decide (L p = .two))⟩

def census (c : CellChart) (L : Source) := censusOn (List.range c.nPorts) L

def base (c : CellChart) : Source := fun p =>
  if p ∈ c.boundary.toList then .two else if p ∈ c.own then .zero else .one

def witness (c : CellChart) (i : Fin 8) : Source := fun p =>
  if p ∈ support c then
    let f := flags i
    let vacant := (p == c.remote.a && f.a) || (p == c.remote.b && f.b) ||
      (p == c.remote.c && f.c) || (p == c.spare.a && !f.a) ||
      (p == c.spare.b && !f.b) || (p == c.spare.c && !f.c)
    if vacant then .zero else .one
  else base c p

def incidenceOK (e : Edge) (x y : Nat) : Bool :=
  decide (x ≠ y) && incident e x && incident e y

def shapeOK (c : CellChart) : Bool :=
  c.inner.length == 30 && c.own.length == 12 && decide c.own.Nodup &&
  c.own.all (· < c.nPorts) && decide c.inner.Nodup && decide c.boundary.toList.Nodup &&
  c.boundary.toList.all c.own.contains && incidenceOK c.outer.a c.boundary.a c.remote.a &&
  incidenceOK c.outer.b c.boundary.b c.remote.b && incidenceOK c.outer.c c.boundary.c c.remote.c

def layoutOK (c : CellChart) : Bool :=
  shapeOK c && decide (support c).Nodup && (support c).all (· < c.nPorts) &&
  (support c).all (fun p => !c.own.contains p) &&
  c.inner.all (fun e => c.own.contains e.lo && c.own.contains e.hi) &&
  c.inner.all (fun e => !c.outer.toList.contains e)

def boundaryFlags (c : CellChart) (wins : List Edge) : Flags :=
  ⟨wins.contains c.outer.a, wins.contains c.outer.b, wins.contains c.outer.c⟩

def fromFlags (c : CellChart) (f : Flags) : List Edge :=
  (if f.a then [c.outer.a] else []) ++ (if f.b then [c.outer.b] else []) ++
    (if f.c then [c.outer.c] else [])

def localWord (c : CellChart) (s : Material) : List Bool :=
  c.inner.map (fun e => (s.collar e).isSome)

def wordFromWins (c : CellChart) (wins : List Edge) : List Bool :=
  c.inner.map wins.contains

def readFlags (c : CellChart) (word : List Bool) : Flags :=
  let used := (c.inner.zip word).filterMap (fun (e, hit) => if hit then some e else none)
  ⟨!touched used c.boundary.a, !touched used c.boundary.b, !touched used c.boundary.c⟩

/-- Eight real sharp preparations; the only varying loads are on six remote sites.
    The finite histogram check costs six ports even for the largest federation. -/
def witnessOK (c : CellChart) : Bool := decide (∀ i : Fin 8,
  let L := witness c i
  let wins := accepted (initialReady L) meets (menu c)
  readFlags c (wordFromWins c wins) = flags i ∧
  boundaryFlags c (accepted (initialReady L) meets c.outer.toList) = flags i ∧
  censusOn (support c) L = (⟨3,3,0⟩ : Three Nat))

def cellCheck (c : CellChart) : Bool := layoutOK c && witnessOK c

private theorem layout_facts {c : CellChart} (h : layoutOK c = true) :
    (support c).Nodup ∧ (∀ p ∈ support c, p < c.nPorts) ∧
    (∀ p ∈ support c, p ∉ c.own) ∧
    (∀ e ∈ c.inner, e.lo ∈ c.own ∧ e.hi ∈ c.own) ∧
    (∀ e ∈ c.inner, e ∉ c.outer.toList) := by
  simp [layoutOK, List.all_eq_true] at h
  exact ⟨h.1.1.1.1.2,h.1.1.1.2,h.1.1.2,h.1.2,h.2⟩

theorem witness_local {c : CellChart} (h : layoutOK c = true) (i : Fin 8) :
    localView c (witness c i) = localView c (base c) := by
  funext p
  by_cases hp : p ∈ c.own
  · have hs : p ∉ support c := fun hs => (layout_facts h).2.2.1 p hs hp
    simp [localView,witness,hp,hs]
  · simp [localView,hp]

private theorem counts_extend {n : Nat} {sites : List Nat} {L R : Source}
    (hd : sites.Nodup) (bound : ∀ p ∈ sites, p < n) (outside : ∀ p ∉ sites, L p = R p)
    (g : Load) (small : sites.countP (fun p => decide (L p = g)) =
      sites.countP (fun p => decide (R p = g))) :
    (List.range n).countP (fun p => decide (L p = g)) =
      (List.range n).countP (fun p => decide (R p = g)) := by
  let allPorts := List.range n
  have perm : (allPorts.filter sites.contains).Perm sites := by
    apply List.perm_iff_count.mpr
    intro p
    have dn : (allPorts.filter sites.contains).Nodup :=
      List.Pairwise.filter _ (List.nodup_range (n := n))
    rw [dn.count, hd.count]
    simp only [List.mem_filter,List.contains_iff_mem,allPorts,List.mem_range]
    by_cases hp : p ∈ sites <;> simp_all
  let left := fun p => decide (L p = g)
  let right := fun p => decide (R p = g)
  have on := (perm.countP_eq left).trans (small.trans (perm.countP_eq right).symm)
  have off : (allPorts.filter (fun p => !sites.contains p)).countP left =
      (allPorts.filter (fun p => !sites.contains p)).countP right := by
    apply List.countP_congr
    intro p hp
    have hn : p ∉ sites := by simpa using (List.mem_filter.mp hp).2
    simp [left,right,outside p hn]
  change allPorts.countP left = allPorts.countP right
  rw [List.countP_eq_countP_filter_add allPorts left sites.contains,
      List.countP_eq_countP_filter_add allPorts right sites.contains,on,off]

private theorem witness_facts {c : CellChart} (h : witnessOK c = true) (i : Fin 8) :
    readFlags c (wordFromWins c (accepted (initialReady (witness c i)) meets (menu c))) = flags i ∧
    boundaryFlags c (accepted (initialReady (witness c i)) meets c.outer.toList) = flags i ∧
    censusOn (support c) (witness c i) = ⟨3,3,0⟩ := (of_decide_eq_true h) i

theorem witness_census {c : CellChart} (h : cellCheck c = true) (i j : Fin 8) :
    census c (witness c i) = census c (witness c j) := by
  obtain ⟨hl,hw⟩ := Bool.and_eq_true_iff.mp h
  obtain ⟨hd,bound,_⟩ := layout_facts hl
  have small := (witness_facts hw i).2.2.trans (witness_facts hw j).2.2.symm
  have hout : ∀ p ∉ support c, witness c i p = witness c j p := by
    intro p hp
    simp [witness,hp]
  unfold census censusOn
  congr 1
  · exact counts_extend hd bound hout .zero (congrArg Three.a small)
  · exact counts_extend hd bound hout .one (congrArg Three.b small)
  · exact counts_extend hd bound hout .two (congrArg Three.c small)

private theorem localWord_material (c : CellChart) (L : Source) (wins : List Edge) :
    localWord c (material L wins) = wordFromWins c wins := by
  unfold localWord wordFromWins
  apply List.map_congr_left
  intro e _
  by_cases he : e ∈ wins <;> simp [material,he]

theorem recover_witness (c : CellChart) (rank : Edge → Nat) (ho : Ordered rank (menu c))
    (h : witnessOK c = true) (i : Fin 8) :
    readFlags c (localWord c (sharpRun rank (menu c) (witness c i))) = flags i := by
  rw [sharpRun_eq_material _ _ _ ho,localWord_material]
  exact (witness_facts h i).1

theorem eight_distinct_records (c : CellChart) (rank : Edge → Nat) (ho : Ordered rank (menu c))
    (h : witnessOK c = true) {i j : Fin 8}
    (eq : localWord c (sharpRun rank (menu c) (witness c i)) =
      localWord c (sharpRun rank (menu c) (witness c j))) : i = j :=
  flags_inj i j ((recover_witness c rank ho h i).symm.trans
    ((congrArg (readFlags c) eq).trans (recover_witness c rank ho h j)))

private theorem decode_touched (c : CellChart) (wins : List Edge) (inside : wins ⊆ c.outer.toList)
    (p : Nat) : touched (fromFlags c (boundaryFlags c wins)) p = touched wins p := by
  apply Bool.eq_iff_iff.mpr
  simp only [touched,List.any_eq_true]
  have hm : ∀ e, e ∈ fromFlags c (boundaryFlags c wins) ↔ e ∈ wins := by
    intro e
    simp only [fromFlags,boundaryFlags,List.mem_append]
    by_cases a : c.outer.a ∈ wins <;> by_cases b : c.outer.b ∈ wins <;>
      by_cases d : c.outer.c ∈ wins <;> simp [a,b,d] <;>
      have hh := @inside e <;> simp [Three.toList] at hh <;> grind
  simp only [hm]

def sendFlags (c : CellChart) (L : Source) : Fin 8 :=
  pack (boundaryFlags c (accepted (initialReady L) meets c.outer.toList))

def replay (c : CellChart) (rank : Edge → Nat) (self : Source) (msg : Fin 8) :=
  sharpRun rank c.inner (informed self (fromFlags c (flags msg)))

private theorem same_word (c : CellChart) {s t : Material} (h : Same c.inner s t) :
    localWord c s = localWord c t :=
  List.map_congr_left (fun e he => congrArg Option.isSome (h.1 e he))

theorem three_flags_suffice (c : CellChart) (rank : Edge → Nat) (ho : Ordered rank (menu c))
    (hl : layoutOK c = true) (L : Source) :
    Same c.inner (sharpRun rank (menu c) L) (replay c rank (localView c L) (sendFlags c L)) := by
  obtain ⟨_,_,_,ins,sep⟩ := layout_facts hl
  let wins := accepted (initialReady L) meets c.outer.toList
  have hb := boundary_flags_suffice rank c.inner c.outer.toList ho sep L
  have hi : Ordered rank c.inner := (List.pairwise_append.mp ho).1
  have touch := decode_touched c wins (fun _ h => accepted_mem h)
  have initialSame : Same c.inner (initial (informed L wins))
      (initial (informed (localView c L) (fromFlags c (boundaryFlags c wins)))) := by
    constructor
    · intro _ _; rfl
    · intro e he
      simp [initial,informed,touch,localView,(ins e he).1,(ins e he).2]
  have localSame := window_local (menu := c.inner) (run := c.inner) (fun _ h => h) initialSame
  rw [serviceFold_eq_material,serviceFold_eq_material,
      ← sharpRun_eq_material rank _ _ hi,← sharpRun_eq_material rank _ _ hi] at localSame
  unfold replay sendFlags
  rw [flags_pack]
  exact ⟨fun e he => (hb.1 e he).trans (localSame.1 e he),
    fun e he => ⟨((hb.2 e he).1).trans (localSame.2 e he).1,
      ((hb.2 e he).2).trans (localSame.2 e he).2⟩⟩

/-- The receiver may know the complete global 0/1/2 census for free. Its remaining
    input is its own preparation and one of at most 2^b classical messages. -/
def Communicates (c : CellChart) (target : Source → List Bool) (b : Nat) : Prop :=
  ∃ send : Source → Fin (2^b), ∃ receive : Source → Three Nat → Fin (2^b) → List Bool,
    ∀ L, receive (localView c L) (census c L) (send L) = target L

private theorem length_of_distinct_subset {α : Type _} [BEq α] [LawfulBEq α] {xs pool : List α}
    (h : xs.Nodup) (inside : xs ⊆ pool) : xs.length ≤ pool.length := by
  induction xs generalizing pool with
  | nil => simp
  | cons x tail ih =>
    obtain ⟨hn,hd⟩ := List.nodup_cons.mp h
    have hx := inside List.mem_cons_self
    have bound := ih (pool := pool.erase x) hd (by
      intro a ha
      have ne : a ≠ x := fun eq => hn (eq ▸ ha)
      exact (List.mem_erase_of_ne ne).mpr (inside (List.mem_cons_of_mem x ha)))
    have hp := List.length_erase_of_mem hx
    have pos := List.length_pos_of_mem hx
    simp only [List.length_cons]
    omega

private theorem eight_map_nodup {α : Type _} {f : Fin 8 → α} (inj : ∀ {i j}, f i = f j → i = j) :
    ((List.finRange 8).map f).Nodup := by
  apply List.pairwise_map.mpr
  exact (show (List.finRange 8).Nodup by decide).imp (fun hn eq => hn (inj eq))

theorem lower_bound (c : CellChart) (rank : Edge → Nat) (ho : Ordered rank (menu c))
    (h : cellCheck c = true) (b : Nat)
    (protocol : Communicates c (fun L => localWord c (sharpRun rank (menu c) L)) b) : 3 ≤ b := by
  obtain ⟨send,receive,correct⟩ := protocol
  obtain ⟨hl,hw⟩ := Bool.and_eq_true_iff.mp h
  let f := fun i => send (witness c i)
  have inj : ∀ {i j}, f i = f j → i = j := by
    intro i j eq
    apply eight_distinct_records c rank ho hw
    apply (correct (witness c i)).symm.trans
    refine Eq.trans ?_ (correct (witness c j))
    rw [witness_local hl i,witness_local hl j,witness_census h i j]
    exact congrArg _ eq
  have bound := length_of_distinct_subset (pool := List.finRange (2^b))
    (eight_map_nodup inj) (fun x _ => List.mem_finRange x)
  simp only [List.length_map,List.length_finRange] at bound
  exact (Nat.pow_le_pow_iff_right (a := 2) (by decide)).mp (show 2^3 ≤ 2^b from bound)

theorem exactly_three (c : CellChart) (rank : Edge → Nat) (ho : Ordered rank (menu c))
    (h : cellCheck c = true) :
    Communicates c (fun L => localWord c (sharpRun rank (menu c) L)) 3 ∧
    ∀ b, Communicates c (fun L => localWord c (sharpRun rank (menu c) L)) b → 3 ≤ b := by
  constructor
  · refine ⟨sendFlags c,fun self _ msg => localWord c (replay c rank self msg),fun L => ?_⟩
    exact (same_word c (three_flags_suffice c rank ho (Bool.and_eq_true_iff.mp h).1 L)).symm
  · exact lower_bound c rank ho h

/-- A lifting interface for the global predecessor closure. In a complete schedule
    the closure may have a different order from the sorted local chart. -/
theorem full_word_of_closed (c : CellChart) (rank : Edge → Nat) (full : List Edge)
    (keep : Edge → Bool) (hf : Ordered rank full) (hm : Ordered rank (menu c))
    (closed : PredecessorClosed rank full keep) (perm : (full.filter keep).Perm (menu c))
    (L : Source) : localWord c (sharpRun rank full L) = localWord c (sharpRun rank (menu c) L) := by
  apply List.map_congr_left
  intro e he
  have em := perm.mem_iff.mpr (List.mem_append_left c.outer.toList he)
  obtain ⟨hef,hek⟩ := List.mem_filter.mp em
  have pf := sharpRun_eq_of_legal_perm rank L perm (List.Pairwise.filter keep hf) hm
  exact congrArg Option.isSome ((closed_agrees rank full keep hf closed L e hef hek).trans
    (congrArg (fun s : Material => s.collar e) pf))

theorem exactly_three_full (c : CellChart) (rank : Edge → Nat) (full : List Edge)
    (keep : Edge → Bool) (hf : Ordered rank full) (hm : Ordered rank (menu c))
    (closed : PredecessorClosed rank full keep) (perm : (full.filter keep).Perm (menu c))
    (h : cellCheck c = true) :
    Communicates c (fun L => localWord c (sharpRun rank full L)) 3 ∧
    ∀ b, Communicates c (fun L => localWord c (sharpRun rank full L)) b → 3 ≤ b := by
  have eq := funext (full_word_of_closed c rank full keep hf hm closed perm)
  rw [eq]
  exact exactly_three c rank hm h

/-- At each internal node the environment may answer an arbitrary Boolean question
    about the source. A leaf reconstructs from the local input and global census. -/
inductive DecisionTree where
  | leaf (receive : Source → Three Nat → List Bool)
  | branch (question : Source → Bool) (no yes : DecisionTree)

def DecisionTree.depth : DecisionTree → Nat
  | .leaf _ => 0
  | .branch _ no yes => max no.depth yes.depth + 1

def DecisionTree.run (c : CellChart) (L : Source) : DecisionTree → List Bool
  | .leaf receive => receive (localView c L) (census c L)
  | .branch question no yes => if question L then yes.run c L else no.run c L

private def DecisionTree.outputs (self : Source) (counts : Three Nat) :
    DecisionTree → List (List Bool)
  | .leaf receive => [receive self counts]
  | .branch _ no yes => no.outputs self counts ++ yes.outputs self counts

private theorem DecisionTree.mem_outputs (t : DecisionTree) (c : CellChart) (L : Source) :
    t.run c L ∈ t.outputs (localView c L) (census c L) := by
  induction t with
  | leaf => simp [run,outputs]
  | branch q no yes hn hy =>
    cases h : q L <;> simp [run,outputs,h,hn,hy]

private theorem DecisionTree.leaf_bound (t : DecisionTree) (self : Source) (counts : Three Nat) :
    (t.outputs self counts).length ≤ 2^t.depth := by
  induction t with
  | leaf => simp [outputs,depth]
  | branch q no yes hn hy =>
    have bn := Nat.pow_le_pow_right (n := 2) (by decide) (Nat.le_max_left no.depth yes.depth)
    have byes := Nat.pow_le_pow_right (n := 2) (by decide) (Nat.le_max_right no.depth yes.depth)
    simp only [outputs,depth,List.length_append,Nat.pow_succ]
    omega

theorem decision_tree_lower (c : CellChart) (rank : Edge → Nat) (ho : Ordered rank (menu c))
    (h : cellCheck c = true) (t : DecisionTree)
    (correct : ∀ L, t.run c L = localWord c (sharpRun rank (menu c) L)) : 3 ≤ t.depth := by
  obtain ⟨hl,hw⟩ := Bool.and_eq_true_iff.mp h
  let target := fun i => localWord c (sharpRun rank (menu c) (witness c i))
  have hd := eight_map_nodup (f := target) (eight_distinct_records c rank ho hw)
  let pool := t.outputs (localView c (base c)) (census c (witness c 0))
  have inside : (List.finRange 8).map target ⊆ pool := by
    intro word hm
    obtain ⟨i,_,rfl⟩ := List.mem_map.mp hm
    have hit := t.mem_outputs c (witness c i)
    rw [correct, witness_local hl i,witness_census h i 0] at hit
    exact hit
  have cap := Nat.le_trans (length_of_distinct_subset hd inside)
    (t.leaf_bound (localView c (base c)) (census c (witness c 0)))
  simp only [List.length_map,List.length_finRange] at cap
  exact (Nat.pow_le_pow_iff_right (a := 2) (by decide)).mp (show 2^3 ≤ 2^t.depth from cap)

def threeQuestions (c : CellChart) (rank : Edge → Nat) : DecisionTree :=
  let bit := fun (L : Source) => flags (sendFlags c L)
  let last := fun a b d => DecisionTree.leaf
    (fun self _ => localWord c (replay c rank self (pack ⟨a,b,d⟩)))
  let third := fun a b => DecisionTree.branch (fun L => (bit L).c) (last a b false) (last a b true)
  let second := fun a => DecisionTree.branch (fun L => (bit L).b) (third a false) (third a true)
  DecisionTree.branch (fun L => (bit L).a) (second false) (second true)

theorem tree_upper (c : CellChart) (rank : Edge → Nat) (ho : Ordered rank (menu c))
    (hl : layoutOK c = true) : (threeQuestions c rank).depth = 3 ∧
    ∀ L, (threeQuestions c rank).run c L = localWord c (sharpRun rank (menu c) L) := by
  refine ⟨rfl,fun L => ?_⟩
  have suff := (same_word c (three_flags_suffice c rank ho hl L)).symm
  have pack_flags : ∀ m : Fin 8, pack (flags m) = m := by decide
  rw [← suff, ← pack_flags (sendFlags c L)]
  simp only [threeQuestions,DecisionTree.run]
  generalize sendFlags c L = msg
  cases h : flags msg with
  | mk a b d =>
    cases a <;> cases b <;> cases d <;> rfl

theorem exactly_three_tree_full (c : CellChart) (rank : Edge → Nat) (full : List Edge)
    (keep : Edge → Bool) (hf : Ordered rank full) (hm : Ordered rank (menu c))
    (closed : PredecessorClosed rank full keep) (perm : (full.filter keep).Perm (menu c))
    (h : cellCheck c = true) :
    ((threeQuestions c rank).depth = 3 ∧
      ∀ L, (threeQuestions c rank).run c L = localWord c (sharpRun rank full L)) ∧
    ∀ t : DecisionTree, (∀ L, t.run c L = localWord c (sharpRun rank full L)) → 3 ≤ t.depth := by
  have eq := full_word_of_closed c rank full keep hf hm closed perm
  obtain ⟨hu,hr⟩ := tree_upper c rank hm (Bool.and_eq_true_iff.mp h).1
  exact ⟨⟨hu,fun L => (hr L).trans (eq L).symm⟩,
    fun t ht => decision_tree_lower c rank hm h t (fun L => (ht L).trans (eq L))⟩

#print axioms witness_census
#print axioms eight_distinct_records
#print axioms exactly_three_full
#print axioms exactly_three_tree_full
end CellBits
