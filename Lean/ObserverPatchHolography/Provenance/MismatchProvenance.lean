import ObserverPatchHolography.Provenance.SemanticEventProvenance

/-!
# Mismatch provenance over versioned semantic commits

An overlap between observer patches carries a quotient-visible mismatch
score that depends only on the committed values on its support.  The score
is direction-blind: it reports disagreement and nothing else.  Direction
enters through provenance.  A commit changes an overlap's mismatch exactly
when the score differs across its snapshots, injects a fresh mismatch when
the score rises from zero, and repairs when the score drops.  A response
reads from an injecting commit when its certified causal support meets the
overlap support at a register whose pre-commit writer is that commit.

Two machine-checked boundaries anchor the design.  A commit whose write set
misses an overlap's support cannot change that overlap's score, so hidden
or presentation-only writes create no mismatch edge.  Two versioned states
with identical values carry identical scores on every overlap while their
writer provenance may differ, and a concrete two-register countermodel
realizes both parent attributions on value-identical snapshots.  A static
residual therefore determines no causal arrow; the writer field is
load-bearing.  A second countermodel shows the read-from clause is
load-bearing for the mismatch-response relation: a repair that does not
read the injected seam version acquires no injection-to-response edge.

The mismatch-response relation is contained in direct semantic parenthood,
so every mismatch edge is also a provenance edge of the generated causal
order.  No physical propagation, spacetime, or continuum object appears.
-/

namespace OPH.Provenance

universe u v w x

variable {Register : Type u} {Value : Type v} {EventId : Type w}
variable {Overlap : Type x}
variable [DecidableEq Register]

/-- A quotient-visible mismatch system: each overlap carries a register
support and a score into the natural numbers, zero meaning agreement, and
the score depends only on the committed values on the support.  Writer
provenance is invisible to the score by construction. -/
structure MismatchSystem (Register : Type u) (Value : Type v)
    (EventId : Type w) (Overlap : Type x) [DecidableEq Register] where
  /-- The register support of each overlap. -/
  support : Overlap → Finset Register
  /-- The mismatch score of each overlap at a versioned state. -/
  score : Overlap → VersionedState Register Value EventId → ℕ
  /-- The score reads only the values on the support. -/
  score_local : ∀ (e : Overlap)
    (q q' : VersionedState Register Value EventId),
    (∀ a ∈ support e, q.value a = q'.value a) → score e q = score e q'

namespace MismatchSystem

variable (M : MismatchSystem Register Value EventId Overlap)

/-- A commit changes an overlap's mismatch when the score differs across
its snapshots. -/
def Changes (c : SemanticCommit Register Value EventId) (e : Overlap) :
    Prop :=
  M.score e c.after ≠ M.score e c.before

/-- A commit injects a fresh mismatch when the score rises from zero. -/
def FreshInjects (c : SemanticCommit Register Value EventId)
    (e : Overlap) : Prop :=
  M.score e c.before = 0 ∧ 0 < M.score e c.after

/-- A commit amplifies an existing mismatch when a positive score rises. -/
def Amplifies (c : SemanticCommit Register Value EventId) (e : Overlap) :
    Prop :=
  0 < M.score e c.before ∧ M.score e c.before < M.score e c.after

/-- A commit repairs an overlap when the score drops. -/
def Repairs (d : SemanticCommit Register Value EventId) (e : Overlap) :
    Prop :=
  M.score e d.after < M.score e d.before

  /-- A response reads from a commit on an overlap when its certified causal
support meets the overlap support at a register carrying the exact version
written by that commit.  A writer identifier alone is insufficient: write
membership and value continuity are part of the witness. -/
def ReadsFrom (d c : SemanticCommit Register Value EventId) (e : Overlap) :
    Prop :=
  ∃ a ∈ d.causalSupp ∩ M.support e, a ∈ c.writeSet ∧
    d.before.writer a = some c.eventId ∧
    d.before.value a = c.after.value a

/-- The distinguished mismatch-response parent relation: a fresh injection
whose changed seam version the response reads and repairs. -/
def MismatchResponseParent (c d : SemanticCommit Register Value EventId) :
    Prop :=
  ∃ e : Overlap, M.FreshInjects c e ∧ M.ReadsFrom d c e ∧ M.Repairs d e

/-- A fresh injection changes the overlap. -/
theorem changes_of_freshInjects
    {c : SemanticCommit Register Value EventId} {e : Overlap}
    (h : M.FreshInjects c e) : M.Changes c e := by
  obtain ⟨hzero, hpos⟩ := h
  intro heq
  rw [heq, hzero] at hpos
  exact lt_irrefl 0 hpos

/-- Every mismatch-response edge is an authenticated direct semantic
parent edge. -/
theorem authenticatedDirectSemanticParent_of_mismatchResponse
    {c d : SemanticCommit Register Value EventId}
    (h : M.MismatchResponseParent c d) :
    AuthenticatedDirectSemanticParent c d := by
  obtain ⟨e, _, ⟨a, hmem, hwrite, hwriter, hvalue⟩, _⟩ := h
  exact ⟨a, (Finset.mem_inter.mp hmem).1, hwrite, hwriter, hvalue⟩

/-- Backward-compatible raw-citation consequence of an authenticated
mismatch-response edge. -/
theorem directSemanticParent_of_mismatchResponse
    {c d : SemanticCommit Register Value EventId}
    (h : M.MismatchResponseParent c d) : DirectSemanticParent c d :=
  directSemanticParent_of_authenticated
    (M.authenticatedDirectSemanticParent_of_mismatchResponse h)

/-- A commit whose write set misses the overlap support leaves that
overlap's score unchanged.  Hidden or presentation-only writes create no
mismatch change. -/
theorem not_changes_of_writes_disjoint
    {c : SemanticCommit Register Value EventId} {e : Overlap}
    (h : ∀ a ∈ M.support e, a ∉ c.writeSet) : ¬ M.Changes c e := by
  have hvalues : ∀ a ∈ M.support e, c.after.value a = c.before.value a :=
    fun a ha => c.frame_value a (h a ha)
  exact fun hchanges => hchanges (M.score_local e c.after c.before hvalues)

/-- Two versioned states with identical values carry identical scores on
every overlap.  The score is writer-blind. -/
theorem score_writer_blind
    (q q' : VersionedState Register Value EventId)
    (hvalues : ∀ a : Register, q.value a = q'.value a) (e : Overlap) :
    M.score e q = M.score e q' :=
  M.score_local e q q' (fun a _ => hvalues a)

end MismatchSystem

/-! ## Static mismatch carries no orientation

The two-register countermodel: one overlap supported on both registers,
scored by disagreement of the two values.  Two pre-commit snapshots hold
the same disagreeing values, so every overlap score agrees, while their
writer fields attribute the standing versions to the two different patch
commits.  The direct-parent verdict for the same response shape differs
between the snapshots.  A static residual determines no causal arrow. -/

namespace StaticOrientationControl

/-- The two-register disagreement system: registers, values, and event
identifiers are booleans, one overlap is supported on both registers, and
the score is one exactly when the two values differ. -/
def system : MismatchSystem Bool Bool Bool Unit where
  support _ := {false, true}
  score _ q := if q.value false = q.value true then 0 else 1
  score_local := by
    intro e q q' hvalues
    have hfalse := hvalues false (by simp)
    have htrue := hvalues true (by simp)
    simp [hfalse, htrue]

/-- The disagreeing snapshot whose standing versions are attributed to the
patch-`false` commit. -/
def leftAttributed : VersionedState Bool Bool Bool where
  value r := r
  writer _ := some false

/-- The same disagreeing values attributed to the patch-`true` commit. -/
def rightAttributed : VersionedState Bool Bool Bool where
  value r := r
  writer _ := some true

/-- A response commit shape over a declared pre-commit snapshot: it reads
and cites both registers, writes both values to `true`, and stamps itself
as writer.  Only the snapshot varies between the two instances. -/
def response (q : VersionedState Bool Bool Bool) :
    SemanticCommit Bool Bool Bool where
  eventId := true
  before := q
  after := { value := fun _ => true, writer := fun _ => some true }
  readSet := {false, true}
  writeSet := {false, true}
  causalSupp := {false, true}
  supp_subset_read := by intro a ha; exact ha
  frame_value := by intro a ha; cases a <;> simp at ha
  frame_writer := by intro a ha; cases a <;> simp at ha
  stamp := by intro a _; rfl

/-- A candidate parent commit carrying the identifier `false`.  It genuinely
writes both displayed versions, so the orientation control does not rely on
a forged writer label. -/
def candidateParent : SemanticCommit Bool Bool Bool where
  eventId := false
  before := { value := fun _ => false, writer := fun _ => none }
  after := leftAttributed
  readSet := ∅
  writeSet := {false, true}
  causalSupp := ∅
  supp_subset_read := by intro a ha; exact absurd ha (Finset.notMem_empty a)
  frame_value := by intro a ha; cases a <;> simp at ha
  frame_writer := by intro a ha; cases a <;> simp at ha
  stamp := by intro a _; rfl

/-- The two snapshots agree on every overlap score, and the score is
positive: the standing mismatch is statically identical. -/
theorem scores_agree_and_positive :
    (∀ e : Unit, system.score e leftAttributed =
      system.score e rightAttributed) ∧
    0 < system.score () leftAttributed := by
  constructor
  · intro e
    exact system.score_writer_blind leftAttributed rightAttributed
      (fun _ => rfl) e
  · simp [system, leftAttributed]

/-- The direct-parent verdict differs between the two value-identical
snapshots: the left attribution makes the candidate a direct semantic
parent of the response, the right attribution does not. -/
theorem orientation_not_determined :
    AuthenticatedDirectSemanticParent candidateParent
        (response leftAttributed) ∧
    ¬ AuthenticatedDirectSemanticParent candidateParent
        (response rightAttributed) := by
  constructor
  · exact ⟨false, by simp [response], by simp [candidateParent], rfl, rfl⟩
  · rintro ⟨a, _, _, hwriter, _⟩
    have : (some true : Option Bool) = some false := hwriter
    simp at this

end StaticOrientationControl

/-! ## The read-from clause is load-bearing

A repair that does not read the injected seam version acquires no
injection-to-response edge, whatever the scores do. -/

namespace UnconsumedInjectionControl

open StaticOrientationControl

/-- A repairing commit over the disagreeing left-attributed snapshot whose
certified causal support is empty: it repairs the overlap without citing
the injected version. -/
def blindRepair : SemanticCommit Bool Bool Bool where
  eventId := true
  before := leftAttributed
  after := { value := fun _ => true, writer := fun _ => some true }
  readSet := ∅
  writeSet := {false, true}
  causalSupp := ∅
  supp_subset_read := by intro a ha; exact absurd ha (Finset.notMem_empty a)
  frame_value := by intro a ha; cases a <;> simp at ha
  frame_writer := by intro a ha; cases a <;> simp at ha
  stamp := by intro a _; rfl

/-- The blind repair repairs the overlap yet reads from no commit on it,
so it is no mismatch response to any commit. -/
theorem repairs_without_reading :
    system.Repairs blindRepair () ∧
    (∀ c : SemanticCommit Bool Bool Bool,
      ¬ system.MismatchResponseParent c blindRepair) := by
  constructor
  · simp [MismatchSystem.Repairs, system, blindRepair, leftAttributed]
  · rintro c ⟨e, _, ⟨a, hmem, _, _, _⟩, _⟩
    have ha := (Finset.mem_inter.mp hmem).1
    exact absurd ha (Finset.notMem_empty a)

end UnconsumedInjectionControl

#print axioms MismatchSystem.directSemanticParent_of_mismatchResponse
#print axioms MismatchSystem.authenticatedDirectSemanticParent_of_mismatchResponse
#print axioms MismatchSystem.not_changes_of_writes_disjoint
#print axioms MismatchSystem.score_writer_blind
#print axioms StaticOrientationControl.scores_agree_and_positive
#print axioms StaticOrientationControl.orientation_not_determined
#print axioms UnconsumedInjectionControl.repairs_without_reading

end OPH.Provenance
