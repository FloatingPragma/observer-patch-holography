import Mathlib.Logic.Relation
import Mathlib.Data.Finset.Basic
import Time.TimeOrderLedger

/-!
# Source-derived semantic event provenance

A quotient-visible semantic register carries the identifier of its last
semantic writer.  An accepted semantic commit reads a snapshot, writes a
declared register set, stamps its own identifier on every written register,
and leaves every other register untouched.  The direct semantic parent
relation reads that provenance: `c` is a direct parent of `d` exactly when
some register in the certified causal read support of `d` still carries the
committed value written by `c` in the pre-commit snapshot of `d`.

The generated causal order is the transitive closure of the direct parent
relation.  Under one declared ancestry rank that strictly increases across
direct parent edges, the closure is irreflexive, transitive, and asymmetric,
so it inhabits the A1 `ObserverRecordOrder` type ledger with no further
choice.  The closure is also the least strict transitive relation containing
every direct parent edge: any candidate precedence that contains the edges
and adds no unsupported comparability equals the generated order exactly.

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

/-- Direct semantic parenthood between two commits: some register in the
child's certified causal support carries the parent's committed version in
the child's pre-commit snapshot.  This is version provenance, not timestamp
order. -/
def DirectSemanticParent (c d : SemanticCommit Register Value EventId) :
    Prop :=
  ∃ a ∈ d.causalSupp, d.before.writer a = some c.eventId

instance [DecidableEq EventId]
    (c d : SemanticCommit Register Value EventId) :
    Decidable (DirectSemanticParent c d) :=
  decidable_of_iff
    (∃ a ∈ d.causalSupp, d.before.writer a = some c.eventId) Iff.rfl

/-- A finite semantic event log: one commit per event identifier, coherent
identifiers, and one declared ancestry rank that strictly increases across
direct parent edges.  The rank witnesses that the ledger is append-only; it
carries no clock, schedule, or physical-time reading. -/
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
    DirectSemanticParent (commitOf c) (commitOf d) → rank c < rank d

namespace SemanticEventLog

variable (L : SemanticEventLog Register Value EventId)

/-- The direct parent relation of the log, read on event identifiers. -/
def ParentEdge (c d : EventId) : Prop :=
  DirectSemanticParent (L.commitOf c) (L.commitOf d)

instance [DecidableEq EventId] (c d : EventId) :
    Decidable (L.ParentEdge c d) :=
  decidable_of_iff
    (DirectSemanticParent (L.commitOf c) (L.commitOf d)) Iff.rfl

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

#print axioms SemanticEventLog.rank_lt_of_generatedBefore
#print axioms SemanticEventLog.generatedBefore_irrefl
#print axioms SemanticEventLog.generatedBefore_asymm
#print axioms SemanticEventLog.generatedRecordOrder
#print axioms SemanticEventLog.generatedBefore_le_of_transitive
#print axioms SemanticEventLog.eq_generatedBefore_of_exact

end OPH.Provenance
