import Mathlib.Logic.Relation
import Mathlib.Data.Finset.Basic
import Time.TimeOrderLedger

/-!
# Source-derived semantic event provenance

A quotient-visible semantic register carries the identifier of its last
semantic writer.  An accepted semantic commit reads a snapshot, writes a
declared register set, stamps its own identifier on every written register,
and leaves every other register untouched.  A raw writer citation is not
yet a causal edge: an arbitrary snapshot can name an event which never
wrote the cited register.  The authenticated direct-parent relation
therefore requires a read-after-write witness consisting of all four facts:
the child certifies the register, the parent wrote it, the child's pre-state
names the parent as last writer, and the child's value is the value in the
parent's post-state.

The generated causal order is the transitive closure of authenticated
direct-parent edges.  `SemanticEventLog` remains the abstract certificate
interface and therefore carries a well-founded rank witness.  The executed
history bridge in `HistoryCausalInvariance.lean` constructs this interface
from a fresh duplicate-free threaded log and derives the rank from commit
position; producers should use that bridge rather than supply an unrelated
rank.  The closure is the least strict transitive relation containing every
authenticated edge.

The ancestry rank is a well-foundedness witness for the append-only semantic
ledger.  It is not a repair schedule, modular parameter, clock reading,
proper time, or global time function, and no such identification is made
here.  The module derives an informational order among semantic commits; it
does not identify that order with physical signal propagation, construct a
spacetime, or select a continuum limit.
-/

namespace OPH.Provenance

open OPH.TimeOrderLedger

universe u v w

/-- A quotient-visible register file together with writer provenance: each
register holds a value and the identifier of the semantic commit whose
written value is present, with `none` marking initial data. -/
structure VersionedState (Register : Type u) (Value : Type v)
    (EventId : Type w) where
  /-- The committed value visible at each register. -/
  value : Register → Value
  /-- The last semantic writer of each register, `none` for initial data. -/
  writer : Register → Option EventId

/-- An accepted semantic commit: pre- and post-commit versioned snapshots,
declared read, write, and certified causal-support register sets, and the
commit's event identifier.  The frame law fixes every register outside the
write set, and the stamp law records the commit as writer of every written
register.  Aborted proposals, retries, and executor stutters are not
commits. -/
structure SemanticCommit (Register : Type u) (Value : Type v)
    (EventId : Type w) [DecidableEq Register] where
  /-- The commit's event identifier. -/
  eventId : EventId
  /-- The pre-commit versioned snapshot. -/
  before : VersionedState Register Value EventId
  /-- The post-commit versioned snapshot. -/
  after : VersionedState Register Value EventId
  /-- The declared read set. -/
  readSet : Finset Register
  /-- The declared write set. -/
  writeSet : Finset Register
  /-- The certified causal read support: the registers on which the
  commit's enablement, accepted payload, visible label, or continuation
  genuinely depends.  Support soundness and completeness are receipt
  obligations of the producing pipeline. -/
  causalSupp : Finset Register
  /-- The causal support only cites read registers. -/
  supp_subset_read : causalSupp ⊆ readSet
  /-- Frame law: registers outside the write set keep value and writer. -/
  frame_value : ∀ a : Register, a ∉ writeSet → after.value a = before.value a
  /-- Frame law for provenance: unwritten registers keep their writer. -/
  frame_writer : ∀ a : Register, a ∉ writeSet → after.writer a = before.writer a
  /-- Stamp law: every written register records this commit as writer. -/
  stamp : ∀ a : Register, a ∈ writeSet → after.writer a = some eventId

variable {Register : Type u} {Value : Type v} {EventId : Type w}
variable [DecidableEq Register]

/-- A raw writer citation between two commits.  This predicate deliberately
does not assert that `c` wrote the cited register or that the cited value
matches `c.after`; it is useful as the negative-control boundary for
untrusted snapshots.  Generated precedence uses
`AuthenticatedDirectSemanticParent`, not this predicate. -/
def DirectSemanticParent (c d : SemanticCommit Register Value EventId) :
    Prop :=
  ∃ a ∈ d.causalSupp, d.before.writer a = some c.eventId

instance [DecidableEq EventId]
    (c d : SemanticCommit Register Value EventId) :
    Decidable (DirectSemanticParent c d) :=
  decidable_of_iff
    (∃ a ∈ d.causalSupp, d.before.writer a = some c.eventId) Iff.rfl

/-- An authenticated read-after-write witness at one register.  The event
identifier is not accepted as provenance by itself: the alleged parent must
actually write the register and the child must read the exact committed
value which appears in that parent's post-state. -/
def ReadAfterWriteAt (c d : SemanticCommit Register Value EventId)
    (a : Register) : Prop :=
  a ∈ d.causalSupp ∧ a ∈ c.writeSet ∧
    d.before.writer a = some c.eventId ∧
    d.before.value a = c.after.value a

/-- Authenticated direct semantic parenthood: at least one certified child
support register carries the exact value-version written by the parent. -/
def AuthenticatedDirectSemanticParent
    (c d : SemanticCommit Register Value EventId) : Prop :=
  ∃ a : Register, ReadAfterWriteAt c d a

instance [DecidableEq EventId] [DecidableEq Value]
    (c d : SemanticCommit Register Value EventId) :
    Decidable (AuthenticatedDirectSemanticParent c d) :=
  decidable_of_iff
    (∃ a ∈ d.causalSupp, a ∈ c.writeSet ∧
      d.before.writer a = some c.eventId ∧
      d.before.value a = c.after.value a) (by
        simp only [AuthenticatedDirectSemanticParent, ReadAfterWriteAt])

/-- Authentication implies the underlying writer citation. -/
theorem directSemanticParent_of_authenticated
    {c d : SemanticCommit Register Value EventId}
    (h : AuthenticatedDirectSemanticParent c d) :
    DirectSemanticParent c d := by
  obtain ⟨a, ha, _, hw, _⟩ := h
  exact ⟨a, ha, hw⟩

/-- An abstract semantic event-log certificate: one commit per event
identifier, coherent identifiers, and a well-founded rank increasing across
authenticated read-after-write edges.  Finiteness is supplied by a
`Fintype EventId` at finite-world use sites rather than hidden in this
structure.  The executed-history constructor derives the rank from a
threaded append-only log; arbitrary inhabitants of this abstract interface
must still provide the certificate.  The rank carries no clock, schedule,
or physical-time reading. -/
structure SemanticEventLog (Register : Type u) (Value : Type v)
    (EventId : Type w) [DecidableEq Register] where
  /-- The commit recorded under each event identifier. -/
  commitOf : EventId → SemanticCommit Register Value EventId
  /-- Identifier coherence. -/
  commitOf_eventId : ∀ e : EventId, (commitOf e).eventId = e
  /-- The declared ancestry rank. -/
  rank : EventId → ℕ
  /-- The rank strictly increases across every direct parent edge. -/
  rank_lt_of_parent : ∀ {c d : EventId},
    AuthenticatedDirectSemanticParent (commitOf c) (commitOf d) →
      rank c < rank d

namespace SemanticEventLog

variable (L : SemanticEventLog Register Value EventId)

/-- The authenticated direct parent relation of the log, read on event
identifiers.  Raw writer labels which fail write-membership or value-version
continuity do not generate edges. -/
def ParentEdge (c d : EventId) : Prop :=
  AuthenticatedDirectSemanticParent (L.commitOf c) (L.commitOf d)

instance [DecidableEq EventId] [DecidableEq Value] (c d : EventId) :
    Decidable (L.ParentEdge c d) :=
  decidable_of_iff
    (AuthenticatedDirectSemanticParent (L.commitOf c) (L.commitOf d)) Iff.rfl

/-- The generated causal precedence: the transitive closure of the direct
parent edges. -/
def GeneratedBefore (c d : EventId) : Prop :=
  Relation.TransGen L.ParentEdge c d

/-- The reflexive generated precedence. -/
def GeneratedBeforeEq (c d : EventId) : Prop :=
  c = d ∨ L.GeneratedBefore c d

/-- The ancestry rank strictly increases along the generated precedence. -/
theorem rank_lt_of_generatedBefore {c d : EventId}
    (h : L.GeneratedBefore c d) : L.rank c < L.rank d := by
  induction h with
  | single hedge => exact L.rank_lt_of_parent hedge
  | tail _ hedge ih => exact lt_trans ih (L.rank_lt_of_parent hedge)

/-- The generated precedence is irreflexive. -/
theorem generatedBefore_irrefl (c : EventId) : ¬ L.GeneratedBefore c c :=
  fun h => lt_irrefl (L.rank c) (L.rank_lt_of_generatedBefore h)

/-- The generated precedence is transitive. -/
theorem generatedBefore_trans {c d e : EventId}
    (hcd : L.GeneratedBefore c d) (hde : L.GeneratedBefore d e) :
    L.GeneratedBefore c e :=
  Relation.TransGen.trans hcd hde

/-- The generated precedence is asymmetric. -/
theorem generatedBefore_asymm {c d : EventId}
    (h : L.GeneratedBefore c d) : ¬ L.GeneratedBefore d c := by
  intro h'
  exact lt_asymm (L.rank_lt_of_generatedBefore h)
    (L.rank_lt_of_generatedBefore h')

/-- The generated precedence packaged as the A1 observer-record order type:
strict order data on event identifiers, produced from provenance rather
than declared. -/
def generatedRecordOrder : ObserverRecordOrder EventId where
  precedes := L.GeneratedBefore
  irrefl := L.generatedBefore_irrefl
  trans := fun h h' => L.generatedBefore_trans h h'

@[simp] theorem generatedRecordOrder_precedes (c d : EventId) :
    (L.generatedRecordOrder).precedes c d ↔ L.GeneratedBefore c d :=
  Iff.rfl

/-- Every direct parent edge lies in the generated precedence. -/
theorem parentEdge_generatedBefore {c d : EventId}
    (h : L.ParentEdge c d) : L.GeneratedBefore c d :=
  Relation.TransGen.single h

/-- Minimality: any transitive relation containing every direct parent
edge contains the generated precedence. -/
theorem generatedBefore_le_of_transitive
    (S : EventId → EventId → Prop)
    (hedges : ∀ c d : EventId, L.ParentEdge c d → S c d)
    (htrans : ∀ {c d e : EventId}, S c d → S d e → S c e)
    {c d : EventId} (h : L.GeneratedBefore c d) : S c d := by
  induction h with
  | single hedge => exact hedges _ _ hedge
  | tail _ hedge ih => exact htrans ih (hedges _ _ hedge)

/-- Exactness: a candidate precedence that contains every direct parent
edge and adds no comparability beyond the generated precedence equals the
generated precedence.  A precedence adapter satisfying these two clauses
verifies the source-derived order; it cannot choose a different one. -/
theorem eq_generatedBefore_of_exact
    (S : EventId → EventId → Prop)
    (hedges : ∀ c d : EventId, L.ParentEdge c d → S c d)
    (htrans : ∀ {c d e : EventId}, S c d → S d e → S c e)
    (hsupported : ∀ c d : EventId, S c d → L.GeneratedBefore c d)
    (c d : EventId) : S c d ↔ L.GeneratedBefore c d :=
  ⟨hsupported c d,
    L.generatedBefore_le_of_transitive S hedges (fun h h' => htrans h h')⟩

end SemanticEventLog

/-! ## Forged-writer negative control -/

namespace ForgedWriterControl

/-- An alleged parent whose identifier is `false` but which writes no
register. -/
def allegedParent : SemanticCommit Unit Bool Bool where
  eventId := false
  before := ⟨fun _ => false, fun _ => none⟩
  after := ⟨fun _ => false, fun _ => none⟩
  readSet := ∅
  writeSet := ∅
  causalSupp := ∅
  supp_subset_read := by simp
  frame_value := fun _ _ => rfl
  frame_writer := fun _ _ => rfl
  stamp := by simp

/-- A child snapshot forged to name the alleged parent as writer. -/
def forgedChild : SemanticCommit Unit Bool Bool where
  eventId := true
  before := ⟨fun _ => true, fun _ => some false⟩
  after := ⟨fun _ => true, fun _ => some false⟩
  readSet := {()}
  writeSet := ∅
  causalSupp := {()}
  supp_subset_read := by simp
  frame_value := fun _ _ => rfl
  frame_writer := fun _ _ => rfl
  stamp := by simp

/-- A writer identifier by itself still creates a raw citation, but it is
rejected by the authenticated relation because the alleged parent did not
write the register. -/
theorem raw_writer_label_is_not_authenticated :
    DirectSemanticParent allegedParent forgedChild ∧
      ¬ AuthenticatedDirectSemanticParent allegedParent forgedChild := by
  constructor
  · exact ⟨(), by simp [forgedChild], rfl⟩
  · rintro ⟨a, _, hwrite, _, _⟩
    exact absurd hwrite (by simp [allegedParent])

end ForgedWriterControl

#print axioms SemanticEventLog.rank_lt_of_generatedBefore
#print axioms SemanticEventLog.generatedBefore_irrefl
#print axioms SemanticEventLog.generatedBefore_asymm
#print axioms SemanticEventLog.generatedRecordOrder
#print axioms SemanticEventLog.generatedBefore_le_of_transitive
#print axioms SemanticEventLog.eq_generatedBefore_of_exact
#print axioms ForgedWriterControl.raw_writer_label_is_not_authenticated

end OPH.Provenance
