import ObserverPatchHolography.ObservablePriority.CellBits

/-! Convert the encoder's three decisions into actual reads of the completed
    signed collars on the full physical screen. -/
namespace CellFlagReadout
open ObserverPatchSharp CellBits

def readBoundary (c : CellChart) (s : Material) : Flags :=
  ⟨(s.collar c.outer.a).isSome, (s.collar c.outer.b).isSome, (s.collar c.outer.c).isSome⟩

private theorem unrelated_write (jobs : List Edge) (e : Edge) (absent : e ∉ jobs)
    (start : Material) : (jobs.foldr service start).collar e = start.collar e := by
  induction jobs with
  | nil => rfl
  | cons a tail ih =>
    have ne : e ≠ a := fun h => absent (List.mem_cons.mpr (Or.inl h))
    have abs : e ∉ tail := fun h => absent (List.mem_cons_of_mem a h)
    rw [List.foldr_cons,service]
    split <;> simp [write,ne,ih abs]

private theorem physical_outer (c : CellChart) (rank : Edge → Nat) (ord : Ordered rank (menu c))
    (hl : layoutOK c = true) (L : Source) (e : Edge) (he : e ∈ c.outer.toList) :
    ((sharpRun rank (menu c) L).collar e).isSome =
      (accepted (initialReady L) meets c.outer.toList).contains e := by
  have sep : c.inner.all (fun f => !c.outer.toList.contains f) = true :=
    (Bool.and_eq_true_iff.mp hl).2
  have absent : e ∉ c.inner := by
    intro hin
    have hn := List.all_eq_true.mp sep e hin
    simp [he] at hn
  rw [sharpRun_eq_material _ _ _ ord,← serviceFold_eq_material,menu,List.foldr_append,
      unrelated_write _ e absent _,serviceFold_eq_material]
  by_cases hit : e ∈ accepted (initialReady L) meets c.outer.toList <;> simp [material,hit]

theorem local_message_is_readout (c : CellChart) (rank : Edge → Nat) (ord : Ordered rank (menu c))
    (hl : layoutOK c = true) (L : Source) :
    pack (readBoundary c (sharpRun rank (menu c) L)) = sendFlags c L := by
  have h := physical_outer c rank ord hl L
  simp only [sendFlags,readBoundary,boundaryFlags,h c.outer.a (by simp [Three.toList]),
    h c.outer.b (by simp [Three.toList]),h c.outer.c (by simp [Three.toList])]

private theorem collar_in_closed_menu (c : CellChart) (rank : Edge → Nat) (full : List Edge)
    (keep : Edge → Bool) (hf : Ordered rank full) (hm : Ordered rank (menu c))
    (closed : PredecessorClosed rank full keep) (perm : (full.filter keep).Perm (menu c))
    (L : Source) (e : Edge) (he : e ∈ menu c) :
    (sharpRun rank full L).collar e = (sharpRun rank (menu c) L).collar e := by
  obtain ⟨ef,ek⟩ := List.mem_filter.mp (perm.mem_iff.mpr he)
  have order := sharpRun_eq_of_legal_perm rank L perm (List.Pairwise.filter keep hf) hm
  exact (closed_agrees rank full keep hf closed L e ef ek).trans (congrArg (·.collar e) order)

theorem full_message_is_readout (c : CellChart) (rank : Edge → Nat) (full : List Edge)
    (keep : Edge → Bool) (hf : Ordered rank full) (hm : Ordered rank (menu c))
    (closed : PredecessorClosed rank full keep) (perm : (full.filter keep).Perm (menu c))
    (hl : layoutOK c = true) (L : Source) :
    pack (readBoundary c (sharpRun rank full L)) = sendFlags c L := by
  have h := collar_in_closed_menu c rank full keep hf hm closed perm L
  rw [← local_message_is_readout c rank hm hl]
  simp only [readBoundary,h c.outer.a (by simp [menu,Three.toList]),
    h c.outer.b (by simp [menu,Three.toList]),h c.outer.c (by simp [menu,Three.toList])]

/-- The concrete receiver accepts the three physical boundary hit bits. Its output
    agrees at every internal signed address for every initial sharp preparation. -/
theorem replay_from_full_hits (c : CellChart) (rank : Edge → Nat) (full : List Edge)
    (keep : Edge → Bool) (hf : Ordered rank full) (hm : Ordered rank (menu c))
    (closed : PredecessorClosed rank full keep) (perm : (full.filter keep).Perm (menu c))
    (hl : layoutOK c = true) (L : Source) (e : Edge) (he : e ∈ c.inner) :
    (sharpRun rank full L).collar e = (replay c rank (localView c L)
      (pack (readBoundary c (sharpRun rank full L)))).collar e := by
  rw [full_message_is_readout c rank full keep hf hm closed perm hl L]
  exact (collar_in_closed_menu c rank full keep hf hm closed perm L e
    (List.mem_append_left _ he)).trans ((three_flags_suffice c rank hm hl L).1 e he)

#print axioms full_message_is_readout
#print axioms replay_from_full_hits
end CellFlagReadout
