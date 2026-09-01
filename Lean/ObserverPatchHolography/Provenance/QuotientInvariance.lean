import ObserverPatchHolography.Provenance.HistoryCausalInvariance

/-!
# Quotient invariance of mismatch provenance

The physical quotient exposes a declared visible register carrier.
Overlap supports and certified causal supports are quotient-visible by
discipline, so every mismatch predicate and every parent verdict reads
only the visible restriction of a versioned state: two states agreeing on
visible values and visible writers carry identical scores, identical
change, injection, and repair verdicts, and identical direct-parent
verdicts.

The hidden boundary is exact in the other direction.  A commit whose
write set misses the visible carrier fixes the visible restriction, so it
changes no visible mismatch, and along any executed history whose
identifiers are fresh it is never cited by a commit with visible
certified support: hidden or presentation-only writes create no mismatch
change and no generated edge.  Both directions are machine-checked, with
the visible-writer trace lemma carrying the run-level control.

No continuum quotient, gauge group, or physical presentation is
constructed here; the visible carrier is declared data.
-/

namespace OPH.Provenance

universe u v w x

variable {Register : Type u} {Value : Type v} {EventId : Type w}
variable {Overlap : Type x}
variable [DecidableEq Register]

/-- Agreement of two versioned states on a visible register carrier:
equal values and equal writers on every visible register. -/
def VisiblyEqual (visible : Finset Register)
    (q q' : VersionedState Register Value EventId) : Prop :=
  ∀ a ∈ visible, q.value a = q'.value a ∧ q.writer a = q'.writer a

namespace MismatchSystem

variable (M : MismatchSystem Register Value EventId Overlap)

/-- A mismatch system is quotient-visible over a carrier when every
overlap support lies inside it. -/
def VisibleOn (visible : Finset Register) : Prop :=
  ∀ e : Overlap, M.support e ⊆ visible

/-- Visibly equal states carry identical scores on every overlap of a
quotient-visible system. -/
theorem score_visible_congr {visible : Finset Register}
    (hvis : M.VisibleOn visible)
    {q q' : VersionedState Register Value EventId}
    (h : VisiblyEqual visible q q') (e : Overlap) :
    M.score e q = M.score e q' :=
  M.score_local e q q' (fun a ha => (h a (hvis e ha)).1)

/-- Change, fresh-injection, and repair verdicts of a commit read only
the visible restrictions of its snapshots. -/
theorem verdicts_visible_congr {visible : Finset Register}
    (hvis : M.VisibleOn visible)
    {c c' : SemanticCommit Register Value EventId}
    (hbefore : VisiblyEqual visible c.before c'.before)
    (hafter : VisiblyEqual visible c.after c'.after) (e : Overlap) :
    (M.Changes c e ↔ M.Changes c' e) ∧
      (M.FreshInjects c e ↔ M.FreshInjects c' e) ∧
      (M.Repairs c e ↔ M.Repairs c' e) := by
  have hb := M.score_visible_congr hvis hbefore e
  have ha := M.score_visible_congr hvis hafter e
  refine ⟨?_, ?_, ?_⟩ <;>
    simp [Changes, FreshInjects, Repairs, hb, ha]

end MismatchSystem

/-- A parent verdict reads only the visible writers when the child's
certified support is visible. -/
theorem directSemanticParent_visible_congr {visible : Finset Register}
    {c : SemanticCommit Register Value EventId}
    {d d' : SemanticCommit Register Value EventId}
    (hsupp : d.causalSupp = d'.causalSupp)
    (hvis : d.causalSupp ⊆ visible)
    (hbefore : VisiblyEqual visible d.before d'.before) :
    DirectSemanticParent c d ↔ DirectSemanticParent c d' := by
  constructor
  · rintro ⟨a, ha, hw⟩
    exact ⟨a, hsupp ▸ ha,
      ((hbefore a (hvis ha)).2).symm.trans hw⟩
  · rintro ⟨a, ha, hw⟩
    have ha' : a ∈ d.causalSupp := hsupp ▸ ha
    exact ⟨a, ha', ((hbefore a (hvis ha')).2).trans hw⟩

/-! ## The hidden boundary -/

/-- A hidden commit specification writes outside the visible carrier. -/
def HiddenSpec (visible : Finset Register)
    (s : CommitSpec Register Value EventId) : Prop :=
  ∀ a ∈ s.writeSet, a ∉ visible

/-- Applying a hidden specification fixes the visible restriction. -/
theorem hiddenSpec_visiblyEqual {visible : Finset Register}
    {s : CommitSpec Register Value EventId}
    (hs : HiddenSpec visible s)
    (q : VersionedState Register Value EventId) :
    VisiblyEqual visible (s.apply q) q := by
  intro a ha
  have hnot : a ∉ s.writeSet := fun hw => hs a hw ha
  exact ⟨CommitSpec.apply_value_notMem hnot,
    CommitSpec.apply_writer_notMem hnot⟩

/-- A hidden commit changes no visible mismatch. -/
theorem hiddenSpec_not_changes {visible : Finset Register}
    (M : MismatchSystem Register Value EventId Overlap)
    (hvis : M.VisibleOn visible)
    {s : CommitSpec Register Value EventId}
    (hs : HiddenSpec visible s)
    (q : VersionedState Register Value EventId) (e : Overlap) :
    ¬ M.Changes (s.commitAt q) e := by
  intro hchange
  exact hchange
    (M.score_visible_congr hvis (hiddenSpec_visiblyEqual hs q) e)

/-- Visible writers along an executed history are initial or stamped by
a specification writing a visible register. -/
theorem execState_visible_writer_cases {visible : Finset Register}
    (q : VersionedState Register Value EventId)
    (specs : List (CommitSpec Register Value EventId))
    {a : Register} (ha : a ∈ visible) :
    (execState q specs).writer a = q.writer a ∨
      ∃ s ∈ specs, a ∈ s.writeSet ∧ ¬ HiddenSpec visible s ∧
        (execState q specs).writer a = some s.eventId := by
  rcases execState_writer_cases q specs a with h | ⟨s, hs, hmem, hw⟩
  · exact Or.inl h
  · exact Or.inr ⟨s, hs, hmem, fun hhidden => hhidden a hmem ha, hw⟩

/-- Duplicate-free identifiers identify specifications. -/
theorem spec_eq_of_eventId_eq
    {specs : List (CommitSpec Register Value EventId)}
    (hnodup : (specs.map CommitSpec.eventId).Nodup)
    {t h : CommitSpec Register Value EventId}
    (ht : t ∈ specs) (hh : h ∈ specs)
    (hid : t.eventId = h.eventId) : t = h := by
  induction specs with
  | nil => exact absurd ht (List.not_mem_nil)
  | cons s rest ih =>
      rw [List.map_cons, List.nodup_cons] at hnodup
      rcases List.mem_cons.mp ht with rfl | htrest
      · rcases List.mem_cons.mp hh with rfl | hhrest
        · rfl
        · exact absurd (hid ▸ List.mem_map_of_mem hhrest) hnodup.1
      · rcases List.mem_cons.mp hh with rfl | hhrest
        · exact absurd (hid ▸ List.mem_map_of_mem htrest) hnodup.1
        · exact ih hnodup.2 htrest hhrest

/-- Along a fresh executed history, a hidden commit is never a generated
parent of any commit whose certified support is visible: hidden writes
are invisible to the generated order. -/
theorem hiddenSpec_no_visible_edge {visible : Finset Register}
    {q : VersionedState Register Value EventId}
    {specs : List (CommitSpec Register Value EventId)}
    (hfresh : FreshStart q specs)
    (hnodup : (specs.map CommitSpec.eventId).Nodup)
    {h : CommitSpec Register Value EventId}
    (hh : h ∈ specs) (hhidden : HiddenSpec visible h)
    (hsupp : ∀ s ∈ specs, s.causalSupp ⊆ visible)
    (f : EventId) : ¬ execParents q specs h.eventId f := by
  intro hedge
  obtain ⟨-, d, hd, -, a, ha, hw⟩ :=
    (execParents_iff q specs h.eventId f).mp hedge
  obtain ⟨xs, s, ys, hsplit, hcommit⟩ := execLog_mem_split hd
  have hsmem : s ∈ specs := by
    rw [hsplit]
    exact List.mem_append_right xs List.mem_cons_self
  have havis : a ∈ visible := by
    apply hsupp s hsmem
    have hds : d.causalSupp = s.causalSupp := by rw [hcommit]; rfl
    rwa [hds] at ha
  have hw' : (execState q xs).writer a = some h.eventId := by
    have : d.before = execState q xs := by rw [hcommit]; rfl
    rwa [this] at hw
  rcases execState_visible_writer_cases q xs havis with hinit | ⟨t, ht, hta, htvis, hwt⟩
  · rw [hinit] at hw'
    exact hfresh h hh a hw'
  · rw [hwt] at hw'
    have hid : t.eventId = h.eventId := by simpa using hw'
    have htspecs : t ∈ specs := by
      rw [hsplit]
      exact List.mem_append_left _ ht
    have hth : t = h := spec_eq_of_eventId_eq hnodup htspecs hh hid
    exact htvis (hth ▸ hhidden)

#print axioms MismatchSystem.score_visible_congr
#print axioms MismatchSystem.verdicts_visible_congr
#print axioms directSemanticParent_visible_congr
#print axioms hiddenSpec_not_changes
#print axioms hiddenSpec_no_visible_edge

end OPH.Provenance
