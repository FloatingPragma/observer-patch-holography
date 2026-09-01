import ObserverPatchHolography.Provenance.MismatchProvenance

/-!
# History invariance of the generated causal order

An implementation history executes commit specifications against a threaded
versioned state; a specification declares its identity, read, write, and
certified support sets, and its written payload.  Executor stutters are
steps that commit nothing: the executed event list and the final state read
only the committed specifications, so stutters generate no event and no
edge.

The schedule-invariance theorems bound what an executor can do to the
generated order.  Two adjacent commits whose write sets are disjoint from
each other and from each other's certified supports commute: swapping them
preserves the final state and the full direct-parent relation on event
identifiers, and therefore the generated causal order.  Histories related
by any chain of such independent swaps carry equal labeled provenance.  A
dependent pair is the boundary: swapping a writer past a reader that
certifies the written register changes the edge set, so the independence
clauses are load-bearing.

On executed histories the append-only ancestry rank is derived rather than
declared: a direct parent always sits at a strictly earlier position of the
executed event list, so the list position witnesses the rank premise on
this branch.

A separate persistence theorem supplies the record-side reading of a
long-lived cut: a freshly injected overlap mismatch survives every
continuation whose commits avoid the overlap support, at full strength, at
every intermediate stage.  Nothing here selects a schedule, derives
fairness, assigns physical time, or identifies the informational order with
physical signal propagation.
-/

namespace OPH.Provenance

universe u v w

variable {Register : Type u} {Value : Type v} {EventId : Type w}
variable [DecidableEq Register]

/-- A commit specification: identity, declared register sets, and the
written payload.  The payload is declared per specification, so execution
threads only provenance. -/
structure CommitSpec (Register : Type u) (Value : Type v)
    (EventId : Type w) [DecidableEq Register] where
  /-- The specification's event identifier. -/
  eventId : EventId
  /-- The declared read set. -/
  readSet : Finset Register
  /-- The declared write set. -/
  writeSet : Finset Register
  /-- The certified causal read support. -/
  causalSupp : Finset Register
  /-- The causal support only cites read registers. -/
  supp_subset_read : causalSupp ⊆ readSet
  /-- The written payload on the write set. -/
  writeValue : Register → Value

namespace CommitSpec

/-- Applying a specification writes its payload and stamps its identifier
on the write set, and fixes every other register. -/
def apply (s : CommitSpec Register Value EventId)
    (q : VersionedState Register Value EventId) :
    VersionedState Register Value EventId where
  value a := if a ∈ s.writeSet then s.writeValue a else q.value a
  writer a := if a ∈ s.writeSet then some s.eventId else q.writer a

@[simp] theorem apply_value_mem {s : CommitSpec Register Value EventId}
    {q : VersionedState Register Value EventId} {a : Register}
    (ha : a ∈ s.writeSet) : (s.apply q).value a = s.writeValue a :=
  if_pos ha

@[simp] theorem apply_value_notMem {s : CommitSpec Register Value EventId}
    {q : VersionedState Register Value EventId} {a : Register}
    (ha : a ∉ s.writeSet) : (s.apply q).value a = q.value a :=
  if_neg ha

@[simp] theorem apply_writer_mem {s : CommitSpec Register Value EventId}
    {q : VersionedState Register Value EventId} {a : Register}
    (ha : a ∈ s.writeSet) : (s.apply q).writer a = some s.eventId :=
  if_pos ha

@[simp] theorem apply_writer_notMem {s : CommitSpec Register Value EventId}
    {q : VersionedState Register Value EventId} {a : Register}
    (ha : a ∉ s.writeSet) : (s.apply q).writer a = q.writer a :=
  if_neg ha

/-- The semantic commit executed by a specification at a state. -/
def commitAt (s : CommitSpec Register Value EventId)
    (q : VersionedState Register Value EventId) :
    SemanticCommit Register Value EventId where
  eventId := s.eventId
  before := q
  after := s.apply q
  readSet := s.readSet
  writeSet := s.writeSet
  causalSupp := s.causalSupp
  supp_subset_read := s.supp_subset_read
  frame_value := fun _ ha => if_neg ha
  frame_writer := fun _ ha => if_neg ha
  stamp := fun _ ha => if_pos ha

/-- Direct parenthood of an executed commit reads only the pre-commit
writers on the certified support. -/
theorem parent_commitAt_iff (c : SemanticCommit Register Value EventId)
    (s : CommitSpec Register Value EventId)
    (q : VersionedState Register Value EventId) :
    DirectSemanticParent c (s.commitAt q) ↔
      ∃ a ∈ s.causalSupp, q.writer a = some c.eventId :=
  Iff.rfl

/-- Two states agreeing on the writers of the certified support execute
the same parent verdicts. -/
theorem parent_commitAt_congr {s : CommitSpec Register Value EventId}
    {q q' : VersionedState Register Value EventId}
    (h : ∀ a ∈ s.causalSupp, q.writer a = q'.writer a)
    (c : SemanticCommit Register Value EventId) :
    DirectSemanticParent c (s.commitAt q) ↔
      DirectSemanticParent c (s.commitAt q') := by
  constructor
  · rintro ⟨a, ha, hw⟩
    exact ⟨a, ha, (h a ha).symm.trans hw⟩
  · rintro ⟨a, ha, hw⟩
    exact ⟨a, ha, (h a ha).trans hw⟩

/-- Two write-disjoint specifications apply in either order to the same
state. -/
theorem apply_comm {s t : CommitSpec Register Value EventId}
    (h : ∀ a : Register, a ∈ s.writeSet → a ∉ t.writeSet)
    (q : VersionedState Register Value EventId) :
    t.apply (s.apply q) = s.apply (t.apply q) := by
  have value_eq : ∀ a : Register,
      (t.apply (s.apply q)).value a = (s.apply (t.apply q)).value a := by
    intro a
    by_cases hs : a ∈ s.writeSet
    · have ht : a ∉ t.writeSet := h a hs
      simp [apply_value_mem, apply_value_notMem, hs, ht]
    · by_cases ht : a ∈ t.writeSet
      · simp [apply_value_mem, apply_value_notMem, hs, ht]
      · simp [apply_value_notMem, hs, ht]
  have writer_eq : ∀ a : Register,
      (t.apply (s.apply q)).writer a = (s.apply (t.apply q)).writer a := by
    intro a
    by_cases hs : a ∈ s.writeSet
    · have ht : a ∉ t.writeSet := h a hs
      simp [apply_writer_mem, apply_writer_notMem, hs, ht]
    · by_cases ht : a ∈ t.writeSet
      · simp [apply_writer_mem, apply_writer_notMem, hs, ht]
      · simp [apply_writer_notMem, hs, ht]
  cases hq : t.apply (s.apply q) with
  | mk v w =>
      cases hq' : s.apply (t.apply q) with
      | mk v' w' =>
          have hv : v = v' := by
            funext a
            have := value_eq a
            rw [hq, hq'] at this
            exact this
          have hw : w = w' := by
            funext a
            have := writer_eq a
            rw [hq, hq'] at this
            exact this
          rw [hv, hw]

end CommitSpec

/-! ## Execution of histories with stutters -/

/-- An implementation step: one commit specification or one executor
stutter. -/
inductive RunStep (Register : Type u) (Value : Type v) (EventId : Type w)
    [DecidableEq Register] where
  /-- A semantic commit step. -/
  | commit (s : CommitSpec Register Value EventId)
  /-- An executor stutter: no record changes and no event. -/
  | stutter

/-- The committed specifications of a history: stutters contribute
nothing. -/
def committedSpecs :
    List (RunStep Register Value EventId) →
      List (CommitSpec Register Value EventId)
  | [] => []
  | RunStep.commit s :: rest => s :: committedSpecs rest
  | RunStep.stutter :: rest => committedSpecs rest

/-- The state after executing a specification list. -/
def execState (q : VersionedState Register Value EventId) :
    List (CommitSpec Register Value EventId) →
      VersionedState Register Value EventId
  | [] => q
  | s :: rest => execState (s.apply q) rest

/-- The executed event list of a specification list. -/
def execLog (q : VersionedState Register Value EventId) :
    List (CommitSpec Register Value EventId) →
      List (SemanticCommit Register Value EventId)
  | [] => []
  | s :: rest => s.commitAt q :: execLog (s.apply q) rest

/-- Inserting a stutter anywhere in a history changes neither the
committed specifications, hence neither the executed event list nor the
final state. -/
theorem committedSpecs_insert_stutter
    (pre post : List (RunStep Register Value EventId)) :
    committedSpecs (pre ++ RunStep.stutter :: post) =
      committedSpecs (pre ++ post) := by
  induction pre with
  | nil => rfl
  | cons head tail ih =>
      cases head with
      | commit s => simp [committedSpecs, ih]
      | stutter => simp [committedSpecs, ih]

@[simp] theorem execState_append (q : VersionedState Register Value EventId)
    (xs ys : List (CommitSpec Register Value EventId)) :
    execState q (xs ++ ys) = execState (execState q xs) ys := by
  induction xs generalizing q with
  | nil => rfl
  | cons s rest ih => simp [execState, ih]

@[simp] theorem execLog_append (q : VersionedState Register Value EventId)
    (xs ys : List (CommitSpec Register Value EventId)) :
    execLog q (xs ++ ys) = execLog q xs ++ execLog (execState q xs) ys := by
  induction xs generalizing q with
  | nil => rfl
  | cons s rest ih => simp [execLog, execState, ih]

/-! ## The executed direct-parent relation -/

/-- The direct-parent relation of an executed history, read on event
identifiers. -/
def execParents (q : VersionedState Register Value EventId)
    (specs : List (CommitSpec Register Value EventId))
    (e f : EventId) : Prop :=
  ∃ c ∈ execLog q specs, ∃ d ∈ execLog q specs,
    c.eventId = e ∧ d.eventId = f ∧ DirectSemanticParent c d

/-- Writers in a threaded state are initial or stamped by an executed
specification. -/
theorem execState_writer_cases
    (q : VersionedState Register Value EventId)
    (specs : List (CommitSpec Register Value EventId)) (a : Register) :
    (execState q specs).writer a = q.writer a ∨
      ∃ s ∈ specs, a ∈ s.writeSet ∧
        (execState q specs).writer a = some s.eventId := by
  induction specs generalizing q with
  | nil => exact Or.inl rfl
  | cons s rest ih =>
      rcases ih (s.apply q) with h | ⟨t, ht, hat, hw⟩
      · by_cases ha : a ∈ s.writeSet
        · exact Or.inr ⟨s, List.mem_cons_self, ha, by
            simpa [execState, CommitSpec.apply_writer_mem ha] using h⟩
        · exact Or.inl (by
            simpa [execState, CommitSpec.apply_writer_notMem ha] using h)
      · exact Or.inr ⟨t, List.mem_cons_of_mem s ht, hat, hw⟩

/-- Executed events carry exactly the specification identifiers, in
order. -/
theorem execLog_ids (q : VersionedState Register Value EventId)
    (specs : List (CommitSpec Register Value EventId)) :
    (execLog q specs).map SemanticCommit.eventId =
      specs.map CommitSpec.eventId := by
  induction specs generalizing q with
  | nil => rfl
  | cons s rest ih => simp [execLog, CommitSpec.commitAt, ih]

/-! ## The append-only rank is derived on executed histories -/

/-- A history is provenance-fresh at a start state when no specification
identifier occurs as a writer of the start state. -/
def FreshStart (q : VersionedState Register Value EventId)
    (specs : List (CommitSpec Register Value EventId)) : Prop :=
  ∀ s ∈ specs, ∀ a : Register, q.writer a ≠ some s.eventId

/-- The parent of an executed commit is always stamped by a strict-prefix
specification: executed provenance is append-only. -/
theorem exec_parent_mem_prefix
    {q : VersionedState Register Value EventId}
    {xs ys : List (CommitSpec Register Value EventId)}
    {s : CommitSpec Register Value EventId}
    (hfresh : FreshStart q (xs ++ s :: ys))
    {c : SemanticCommit Register Value EventId}
    (hclog : c.eventId ∈ (xs ++ s :: ys).map CommitSpec.eventId)
    (hparent : DirectSemanticParent c (s.commitAt (execState q xs))) :
    c.eventId ∈ xs.map CommitSpec.eventId := by
  obtain ⟨a, _, hw⟩ := hparent
  have hw' : (execState q xs).writer a = some c.eventId := hw
  rcases execState_writer_cases q xs a with hinit | ⟨t, ht, _, hwt⟩
  · rw [hinit] at hw'
    obtain ⟨e, he, hid⟩ := List.mem_map.mp hclog
    exact absurd (hw'.trans (congrArg some hid.symm)) (hfresh e he a)
  · rw [hwt] at hw'
    have hid : t.eventId = c.eventId := by
      simpa using hw'
    exact hid ▸ List.mem_map_of_mem ht

/-- Identifier rank in a list: the position of the first occurrence. -/
def idRank [DecidableEq EventId] :
    List EventId → EventId → ℕ
  | [], _ => 0
  | x :: rest, e => if x = e then 0 else idRank rest e + 1

theorem idRank_append_self [DecidableEq EventId]
    {l₁ l₂ : List EventId} {e : EventId} (h : e ∉ l₁) :
    idRank (l₁ ++ e :: l₂) e = l₁.length := by
  induction l₁ with
  | nil => simp [idRank]
  | cons x rest ih =>
      have hxe : x ≠ e := fun hx => h (hx ▸ List.mem_cons_self)
      have hrest : e ∉ rest := fun hr => h (List.mem_cons_of_mem x hr)
      simp [idRank, hxe, ih hrest]

theorem idRank_lt_of_mem_prefix [DecidableEq EventId]
    {l₁ l₂ : List EventId} {e : EventId} (h : e ∈ l₁) :
    idRank (l₁ ++ l₂) e < l₁.length := by
  induction l₁ with
  | nil => exact absurd h (List.not_mem_nil)
  | cons x rest ih =>
      by_cases hxe : x = e
      · simp [idRank, hxe]
      · have hrest : e ∈ rest := by
          rcases List.mem_cons.mp h with hx | hr
          · exact absurd hx.symm hxe
          · exact hr
        have := ih hrest
        simp only [List.cons_append, idRank, hxe, if_false,
          List.length_cons]
        omega

/-- Every executed event is the commit of one specification at the state
after a strict prefix. -/
theorem execLog_mem_split
    {q : VersionedState Register Value EventId}
    {specs : List (CommitSpec Register Value EventId)}
    {d : SemanticCommit Register Value EventId}
    (hd : d ∈ execLog q specs) :
    ∃ xs s ys, specs = xs ++ s :: ys ∧
      d = CommitSpec.commitAt s (execState q xs) := by
  induction specs generalizing q with
  | nil => simp [execLog] at hd
  | cons s rest ih =>
      rcases List.mem_cons.mp hd with hds | hdrest
      · exact ⟨[], s, rest, by simp, by simpa [execState] using hds⟩
      · obtain ⟨xs, t, ys, hsplit, hcommit⟩ := ih hdrest
        exact ⟨s :: xs, t, ys, by rw [hsplit]; rfl, by
          simpa [execState] using hcommit⟩

/-- On a fresh, duplicate-free executed history the specification
position is an ancestry rank: every executed direct-parent pair is
position-ordered.  The append-only rank premise is derived on this
branch. -/
theorem execParents_rank_lt [DecidableEq EventId]
    {q : VersionedState Register Value EventId}
    {specs : List (CommitSpec Register Value EventId)}
    (hfresh : FreshStart q specs)
    (hnodup : (specs.map CommitSpec.eventId).Nodup)
    {e f : EventId} (h : execParents q specs e f) :
    idRank (specs.map CommitSpec.eventId) e <
      idRank (specs.map CommitSpec.eventId) f := by
  classical
  obtain ⟨c, hc, d, hd, hce, hdf, hparent⟩ := h
  obtain ⟨xs, s, ys, hsplit, hcommit⟩ := execLog_mem_split hd
  subst hsplit
  have hce' : c.eventId ∈
      (xs ++ s :: ys).map CommitSpec.eventId := by
    have := List.mem_map_of_mem (f := SemanticCommit.eventId) hc
    rwa [execLog_ids] at this
  have hprefix : c.eventId ∈ xs.map CommitSpec.eventId := by
    refine exec_parent_mem_prefix hfresh hce' ?_
    rw [← hcommit]
    exact hparent
  have hf : f = s.eventId := by
    rw [← hdf, hcommit]
    rfl
  have hids :
      (xs ++ s :: ys).map CommitSpec.eventId =
        xs.map CommitSpec.eventId ++
          s.eventId :: ys.map CommitSpec.eventId := by
    simp
  have hnot : s.eventId ∉ xs.map CommitSpec.eventId := by
    rw [hids] at hnodup
    have hdisj :=
      (List.nodup_append.mp hnodup).2.2
    intro hmem
    exact hdisj s.eventId hmem s.eventId List.mem_cons_self rfl
  have hlt :
      idRank ((xs ++ s :: ys).map CommitSpec.eventId) c.eventId <
        (xs.map CommitSpec.eventId).length := by
    rw [hids]
    exact idRank_lt_of_mem_prefix hprefix
  have heq :
      idRank ((xs ++ s :: ys).map CommitSpec.eventId) s.eventId =
        (xs.map CommitSpec.eventId).length := by
    rw [hids]
    exact idRank_append_self hnot
  rw [← hce, hf, heq]
  exact hlt

/-! ## Independent commutation -/

/-- A specification cites an event at a state: it carries the child
identifier and certifies a register whose writer is the cited event. -/
def SpecCites (s : CommitSpec Register Value EventId)
    (q : VersionedState Register Value EventId) (e f : EventId) : Prop :=
  s.eventId = f ∧ ∃ a ∈ s.causalSupp, q.writer a = some e

/-- The executed child-citation relation. -/
def execChildCites (q : VersionedState Register Value EventId)
    (specs : List (CommitSpec Register Value EventId))
    (e f : EventId) : Prop :=
  ∃ d ∈ execLog q specs, d.eventId = f ∧
    ∃ a ∈ d.causalSupp, d.before.writer a = some e

theorem execChildCites_nil (q : VersionedState Register Value EventId)
    (e f : EventId) : ¬ execChildCites q [] e f := by
  rintro ⟨d, hd, _⟩
  simp [execLog] at hd

theorem execChildCites_cons (q : VersionedState Register Value EventId)
    (s : CommitSpec Register Value EventId)
    (rest : List (CommitSpec Register Value EventId)) (e f : EventId) :
    execChildCites q (s :: rest) e f ↔
      SpecCites s q e f ∨ execChildCites (s.apply q) rest e f := by
  constructor
  · rintro ⟨d, hd, hdf, a, ha, hw⟩
    rcases List.mem_cons.mp hd with hds | hdrest
    · subst hds
      exact Or.inl ⟨hdf, a, ha, hw⟩
    · exact Or.inr ⟨d, hdrest, hdf, a, ha, hw⟩
  · rintro (⟨hdf, a, ha, hw⟩ | ⟨d, hd, hdf, a, ha, hw⟩)
    · exact ⟨s.commitAt q, List.mem_cons_self, hdf, a, ha, hw⟩
    · exact ⟨d, List.mem_cons_of_mem _ hd, hdf, a, ha, hw⟩

theorem execChildCites_append (q : VersionedState Register Value EventId)
    (xs ys : List (CommitSpec Register Value EventId)) (e f : EventId) :
    execChildCites q (xs ++ ys) e f ↔
      execChildCites q xs e f ∨
        execChildCites (execState q xs) ys e f := by
  constructor
  · rintro ⟨d, hd, hrest⟩
    rw [execLog_append] at hd
    rcases List.mem_append.mp hd with hx | hy
    · exact Or.inl ⟨d, hx, hrest⟩
    · exact Or.inr ⟨d, hy, hrest⟩
  · rintro (⟨d, hd, hrest⟩ | ⟨d, hd, hrest⟩)
    · exact ⟨d, by rw [execLog_append]; exact List.mem_append_left _ hd,
        hrest⟩
    · exact ⟨d, by rw [execLog_append]; exact List.mem_append_right _ hd,
        hrest⟩

/-- The executed parent relation splits into identifier presence and a
child citation. -/
theorem execParents_iff (q : VersionedState Register Value EventId)
    (specs : List (CommitSpec Register Value EventId)) (e f : EventId) :
    execParents q specs e f ↔
      e ∈ specs.map CommitSpec.eventId ∧ execChildCites q specs e f := by
  constructor
  · rintro ⟨c, hc, d, hd, hce, hdf, a, ha, hw⟩
    refine ⟨?_, d, hd, hdf, a, ha, by rw [hw, hce]⟩
    have := List.mem_map_of_mem (f := SemanticCommit.eventId) hc
    rw [execLog_ids] at this
    exact hce ▸ this
  · rintro ⟨hmem, d, hd, hdf, a, ha, hw⟩
    rw [← execLog_ids q specs] at hmem
    obtain ⟨c, hc, hce⟩ := List.mem_map.mp hmem
    exact ⟨c, hc, d, hd, hce, hdf, a, ha, by rw [hw, hce]⟩

/-- Independence of two specifications: disjoint writes, and neither
writes into the other's certified support. -/
def IndependentPair (c d : CommitSpec Register Value EventId) : Prop :=
  (∀ a ∈ c.writeSet, a ∉ d.writeSet) ∧
    (∀ a ∈ c.writeSet, a ∉ d.causalSupp) ∧
    (∀ a ∈ d.writeSet, a ∉ c.causalSupp)

theorem IndependentPair.symm {c d : CommitSpec Register Value EventId}
    (h : IndependentPair c d) : IndependentPair d c :=
  ⟨fun a ha hb => h.1 a hb ha, h.2.2, h.2.1⟩

/-- Swapping an independent adjacent pair preserves the final state. -/
theorem swap_execState {c d : CommitSpec Register Value EventId}
    (h : IndependentPair c d)
    (q : VersionedState Register Value EventId)
    (pre post : List (CommitSpec Register Value EventId)) :
    execState q (pre ++ c :: d :: post) =
      execState q (pre ++ d :: c :: post) := by
  rw [execState_append, execState_append]
  show execState (d.apply (c.apply (execState q pre))) post =
    execState (c.apply (d.apply (execState q pre))) post
  rw [CommitSpec.apply_comm h.1]

theorem specCites_congr {s : CommitSpec Register Value EventId}
    {q q' : VersionedState Register Value EventId}
    (h : ∀ a ∈ s.causalSupp, q.writer a = q'.writer a) (e f : EventId) :
    SpecCites s q e f ↔ SpecCites s q' e f := by
  constructor
  · rintro ⟨hf, a, ha, hw⟩
    exact ⟨hf, a, ha, (h a ha).symm.trans hw⟩
  · rintro ⟨hf, a, ha, hw⟩
    exact ⟨hf, a, ha, (h a ha).trans hw⟩

/-- Swapping an independent adjacent pair preserves the executed
direct-parent relation on event identifiers. -/
theorem swap_execParents {c d : CommitSpec Register Value EventId}
    (h : IndependentPair c d)
    (q : VersionedState Register Value EventId)
    (pre post : List (CommitSpec Register Value EventId)) (e f : EventId) :
    execParents q (pre ++ c :: d :: post) e f ↔
      execParents q (pre ++ d :: c :: post) e f := by
  rw [execParents_iff, execParents_iff]
  have hmem : e ∈ (pre ++ c :: d :: post).map CommitSpec.eventId ↔
      e ∈ (pre ++ d :: c :: post).map CommitSpec.eventId := by
    simp only [List.map_append, List.map_cons, List.mem_append,
      List.mem_cons]
    tauto
  have hqp := execState q pre
  have hcqp : ∀ a ∈ c.causalSupp,
      (execState q pre).writer a =
        (d.apply (execState q pre)).writer a := by
    intro a ha
    have : a ∉ d.writeSet := fun hw => h.2.2 a hw ha
    rw [CommitSpec.apply_writer_notMem this]
  have hdqp : ∀ a ∈ d.causalSupp,
      (c.apply (execState q pre)).writer a =
        (execState q pre).writer a := by
    intro a ha
    have : a ∉ c.writeSet := fun hw => h.2.1 a hw ha
    rw [CommitSpec.apply_writer_notMem this]
  have hstate : d.apply (c.apply (execState q pre)) =
      c.apply (d.apply (execState q pre)) :=
    CommitSpec.apply_comm h.1 (execState q pre)
  have hcites :
      execChildCites q (pre ++ c :: d :: post) e f ↔
        execChildCites q (pre ++ d :: c :: post) e f := by
    rw [execChildCites_append, execChildCites_append,
      execChildCites_cons, execChildCites_cons,
      execChildCites_cons, execChildCites_cons]
    show _ ∨ SpecCites c (execState q pre) e f ∨
        SpecCites d (c.apply (execState q pre)) e f ∨
        execChildCites (execState (c.apply (execState q pre)) [d]) post e f
      ↔ _ ∨ SpecCites d (execState q pre) e f ∨
        SpecCites c (d.apply (execState q pre)) e f ∨
        execChildCites (execState (d.apply (execState q pre)) [c]) post e f
    have h1 : SpecCites c (execState q pre) e f ↔
        SpecCites c (d.apply (execState q pre)) e f :=
      specCites_congr hcqp e f
    have h2 : SpecCites d (c.apply (execState q pre)) e f ↔
        SpecCites d (execState q pre) e f :=
      specCites_congr (fun a ha => (hdqp a ha)) e f
    have h3 : execState (c.apply (execState q pre)) [d] =
        execState (d.apply (execState q pre)) [c] := by
      show d.apply (c.apply (execState q pre)) =
        c.apply (d.apply (execState q pre))
      exact hstate
    rw [h3]
    tauto
  rw [hmem, hcites]

/-- One independent adjacent swap between two histories. -/
def IndependentSwap
    (l l' : List (CommitSpec Register Value EventId)) : Prop :=
  ∃ pre c d post, IndependentPair c d ∧
    l = pre ++ c :: d :: post ∧ l' = pre ++ d :: c :: post

/-- Histories related by any chain of independent adjacent swaps execute
to the same final state and the same direct-parent relation, hence the
same generated causal order. -/
theorem eqvGen_independentSwap_invariant
    {l l' : List (CommitSpec Register Value EventId)}
    (h : Relation.EqvGen IndependentSwap l l')
    (q : VersionedState Register Value EventId) :
    execState q l = execState q l' ∧
      ∀ e f : EventId, execParents q l e f ↔ execParents q l' e f := by
  induction h with
  | rel x y hxy =>
      obtain ⟨pre, c, d, post, hpair, hx, hy⟩ := hxy
      subst hx hy
      exact ⟨swap_execState hpair q pre post,
        fun e f => swap_execParents hpair q pre post e f⟩
  | refl x => exact ⟨rfl, fun _ _ => Iff.rfl⟩
  | symm x y _ ih =>
      exact ⟨ih.1.symm, fun e f => (ih.2 e f).symm⟩
  | trans x y z _ _ ih₁ ih₂ =>
      exact ⟨ih₁.1.trans ih₂.1,
        fun e f => (ih₁.2 e f).trans (ih₂.2 e f)⟩

/-- Stutters are history-transparent: inserting an executor stutter
anywhere changes neither the executed state nor the executed parent
relation. -/
theorem stutter_transparent
    (q : VersionedState Register Value EventId)
    (pre post : List (RunStep Register Value EventId)) :
    execState q (committedSpecs (pre ++ RunStep.stutter :: post)) =
        execState q (committedSpecs (pre ++ post)) ∧
      ∀ e f : EventId,
        execParents q (committedSpecs (pre ++ RunStep.stutter :: post)) e f ↔
          execParents q (committedSpecs (pre ++ post)) e f := by
  rw [committedSpecs_insert_stutter]
  exact ⟨rfl, fun _ _ => Iff.rfl⟩

/-! ## Persistence of unconsumed mismatch -/

/-- A continuation that avoids an overlap support leaves that overlap's
score unchanged. -/
theorem execState_score_eq {Overlap : Type v}
    (M : MismatchSystem Register Value EventId Overlap) {e : Overlap}
    (q : VersionedState Register Value EventId)
    (specs : List (CommitSpec Register Value EventId))
    (h : ∀ s ∈ specs, ∀ a ∈ M.support e, a ∉ s.writeSet) :
    M.score e (execState q specs) = M.score e q := by
  induction specs generalizing q with
  | nil => rfl
  | cons s rest ih =>
      have hstep : M.score e (s.apply q) = M.score e q :=
        M.score_local e (s.apply q) q
          (fun a ha =>
            CommitSpec.apply_value_notMem (h s List.mem_cons_self a ha))
      calc M.score e (execState (s.apply q) rest)
          = M.score e (s.apply q) :=
            ih (s.apply q) (fun t ht => h t (List.mem_cons_of_mem s ht))
        _ = M.score e q := hstep

/-- A standing mismatch that no continuation commit consumes persists at
full strength through every stage: the record-side content of a
long-lived cut. -/
theorem unconsumed_mismatch_persists {Overlap : Type v}
    (M : MismatchSystem Register Value EventId Overlap) {e : Overlap}
    (q : VersionedState Register Value EventId)
    (specs : List (CommitSpec Register Value EventId))
    (hpos : 0 < M.score e q)
    (h : ∀ s ∈ specs, ∀ a ∈ M.support e, a ∉ s.writeSet)
    {xs ys : List (CommitSpec Register Value EventId)}
    (hsplit : specs = xs ++ ys) :
    0 < M.score e (execState q xs) := by
  have hxs : ∀ s ∈ xs, ∀ a ∈ M.support e, a ∉ s.writeSet := by
    intro s hs
    exact h s (hsplit ▸ List.mem_append_left ys hs)
  rw [execState_score_eq M q xs hxs]
  exact hpos

/-! ## The dependent boundary

Swapping a writer past a reader that certifies the written register
changes the executed edge set: the independence clauses are
load-bearing. -/

namespace DependentSwapControl

/-- The writer specification: event `0` writes the single register. -/
def writerSpec : CommitSpec Unit Bool (Fin 2) where
  eventId := 0
  readSet := ∅
  writeSet := {()}
  causalSupp := ∅
  supp_subset_read := by decide
  writeValue _ := true

/-- The reader specification: event `1` certifies the single register and
writes nothing. -/
def readerSpec : CommitSpec Unit Bool (Fin 2) where
  eventId := 1
  readSet := {()}
  writeSet := ∅
  causalSupp := {()}
  supp_subset_read := by decide
  writeValue _ := false

/-- The blank start state. -/
def start : VersionedState Unit Bool (Fin 2) :=
  ⟨fun _ => false, fun _ => none⟩

/-- Writer before reader executes the provenance edge; reader before
writer does not.  The pair fails independence exactly because the writer
writes into the reader's certified support. -/
theorem dependent_swap_changes_edges :
    execParents start [writerSpec, readerSpec] 0 1 ∧
      ¬ execParents start [readerSpec, writerSpec] 0 1 ∧
      ¬ IndependentPair writerSpec readerSpec := by
  refine ⟨?_, ?_, ?_⟩
  · refine (execParents_iff start _ 0 1).mpr ⟨by simp [writerSpec], ?_⟩
    rw [execChildCites_cons]
    refine Or.inr ?_
    rw [execChildCites_cons]
    refine Or.inl ⟨rfl, (), by simp [readerSpec], ?_⟩
    exact CommitSpec.apply_writer_mem (by simp [writerSpec])
  · intro h
    obtain ⟨_, hcites⟩ := (execParents_iff start _ 0 1).mp h
    rw [execChildCites_cons] at hcites
    rcases hcites with ⟨_, a, _, hw⟩ | hcites
    · exact absurd hw (by simp [start])
    · rw [execChildCites_cons] at hcites
      rcases hcites with ⟨hid, _⟩ | hcites
      · exact absurd hid (by decide)
      · exact execChildCites_nil _ _ _ hcites
  · rintro ⟨_, hsupp, _⟩
    exact hsupp () (by simp [writerSpec]) (by simp [readerSpec])

end DependentSwapControl

/-! ## The diamond under two schedules

The two response specifications of the Boolean diamond are independent,
so the serialization order of the responses leaves the final state and
the full parent relation unchanged: the committed instance of history
invariance. -/

namespace DiamondScheduleWitness

/-- The injection: event `0` writes the seam register. -/
def injectionSpec : CommitSpec (Fin 4) Bool (Fin 4) where
  eventId := 0
  readSet := ∅
  writeSet := {0}
  causalSupp := ∅
  supp_subset_read := by decide
  writeValue _ := true

/-- The first response: certifies the seam, writes its own output. -/
def leftSpec : CommitSpec (Fin 4) Bool (Fin 4) where
  eventId := 1
  readSet := {0}
  writeSet := {1}
  causalSupp := {0}
  supp_subset_read := by decide
  writeValue _ := true

/-- The second response: certifies the seam, writes a disjoint output. -/
def rightSpec : CommitSpec (Fin 4) Bool (Fin 4) where
  eventId := 2
  readSet := {0}
  writeSet := {2}
  causalSupp := {0}
  supp_subset_read := by decide
  writeValue _ := true

/-- The answer: certifies both outputs. -/
def answerSpec : CommitSpec (Fin 4) Bool (Fin 4) where
  eventId := 3
  readSet := {1, 2}
  writeSet := {3}
  causalSupp := {1, 2}
  supp_subset_read := by decide
  writeValue _ := true

/-- The blank start state. -/
def start : VersionedState (Fin 4) Bool (Fin 4) :=
  ⟨fun _ => false, fun _ => none⟩

/-- The two responses are independent. -/
theorem responses_independent : IndependentPair leftSpec rightSpec := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

/-- The two serializations of the diamond execute to the same final state
and the same direct-parent relation. -/
theorem diamond_schedule_invariant :
    execState start [injectionSpec, leftSpec, rightSpec, answerSpec] =
        execState start [injectionSpec, rightSpec, leftSpec, answerSpec] ∧
      ∀ e f : Fin 4,
        execParents start
            [injectionSpec, leftSpec, rightSpec, answerSpec] e f ↔
          execParents start
            [injectionSpec, rightSpec, leftSpec, answerSpec] e f :=
  ⟨swap_execState responses_independent start [injectionSpec] [answerSpec],
    fun e f =>
      swap_execParents responses_independent start
        [injectionSpec] [answerSpec] e f⟩

/-- The invariant relation is nonempty: the injection parents the first
response in both serializations. -/
theorem injection_parents_left :
    execParents start [injectionSpec, leftSpec, rightSpec, answerSpec]
      0 1 := by
  refine (execParents_iff start _ 0 1).mpr ⟨by simp [injectionSpec], ?_⟩
  rw [execChildCites_cons]
  refine Or.inr ?_
  rw [execChildCites_cons]
  refine Or.inl ⟨rfl, 0, by simp [leftSpec], ?_⟩
  exact CommitSpec.apply_writer_mem (by simp [injectionSpec])

end DiamondScheduleWitness

#print axioms execParents_rank_lt
#print axioms swap_execState
#print axioms swap_execParents
#print axioms eqvGen_independentSwap_invariant
#print axioms stutter_transparent
#print axioms execState_score_eq
#print axioms unconsumed_mismatch_persists
#print axioms DependentSwapControl.dependent_swap_changes_edges
#print axioms DiamondScheduleWitness.diamond_schedule_invariant
#print axioms DiamondScheduleWitness.injection_parents_left

end OPH.Provenance
