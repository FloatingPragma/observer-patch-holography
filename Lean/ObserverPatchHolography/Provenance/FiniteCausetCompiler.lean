import ObserverPatchHolography.Provenance.SemanticEventProvenance

/-!
# Finite ranked-relation compiler into authenticated provenance

This module proves an expressivity theorem for the source-derived causal
grammar.  Given a finite event type, a decidable relation `R`, and a supplied
natural-number rank which strictly increases on `R`, it constructs semantic
commits with one register per event such that authenticated direct parenthood
is exactly `R`.  Generated precedence is therefore the transitive closure of
`R`; when `R` is already transitive, generated precedence is exactly `R`.

The construction uses event `c`'s own register as its source token.  Commit
`c` writes that register, and commit `d` certifies/read-cites it precisely when
`R c d`.  Unit values make value continuity exact while the writer stamp
authenticates the cited event identity.

When `R` is supplied as the full strict order of a poset, this generic
compiler therefore authenticates every strict predecessor as a `ParentEdge`;
it does not first reduce `R` to its Hasse/cover relation.  The separate
geometry-seeded simulation compiler uses authenticated geometric covers and
recovers the full order by transitive closure.  The companion
`FiniteCausetCoverCompiler` provides the distinct abstract strengthening that
authenticates exactly the Hasse covers and reconstructs `R` by transitive
closure.

This is a finite abstract-log compiler, not a source-selection theorem.  It
does not derive `R`, a rank, physical event identity, a signal relation,
manifoldlikeness, density, or a continuum.  Nor does it prove that the
compiled per-commit snapshots arise from one particular sequential executor;
that stronger history-custody condition remains a separate interface.
-/

namespace OPH.Provenance.FiniteCausetCompiler

universe u

variable {Event : Type u} [Fintype Event] [DecidableEq Event]
variable (R : Event → Event → Prop) [DecidableRel R]

/-- Registers cited by `d`: exactly its supplied `R`-predecessors. -/
def predecessorRegisters (d : Event) : Finset Event :=
  Finset.univ.filter fun c ↦ R c d

/-- Pre-state for `d`: every cited predecessor register names its own event
as last writer.  Uncited registers carry no writer. -/
def beforeState (d : Event) : VersionedState Event Unit Event where
  value := fun _ ↦ ()
  writer := fun c ↦ if R c d then some c else none

/-- Post-state for `d`: `d` stamps its own singleton register and leaves all
other registers unchanged. -/
def afterState (d : Event) : VersionedState Event Unit Event where
  value := fun _ ↦ ()
  writer := fun c ↦
    if c = d then some d else (beforeState R d).writer c

/-- The semantic commit compiling one event of the supplied relation. -/
def commit (d : Event) : SemanticCommit Event Unit Event where
  eventId := d
  before := beforeState R d
  after := afterState R d
  readSet := predecessorRegisters R d
  writeSet := {d}
  causalSupp := predecessorRegisters R d
  supp_subset_read := fun _ h ↦ h
  frame_value := by
    intro _ _
    rfl
  frame_writer := by
    intro c hc
    have hne : c ≠ d := by simpa using hc
    simp [afterState, hne]
  stamp := by
    intro c hc
    have hcd : c = d := by simpa using hc
    subst c
    simp [afterState]

/-- The compiler is exact already at the authenticated direct-edge level.
For a transitive input `R`, these edges are all supplied strict-order pairs;
the compiler does not reduce them to Hasse covers. -/
theorem authenticatedParent_commit_iff (c d : Event) :
    AuthenticatedDirectSemanticParent (commit R c) (commit R d) ↔ R c d := by
  constructor
  · rintro ⟨a, hasupp, hawrite, _, _⟩
    have hac : a = c := by
      simpa [commit] using hawrite
    subst a
    simpa [commit, predecessorRegisters] using hasupp
  · intro h
    refine ⟨c, ?_, ?_, ?_, ?_⟩
    · simp [commit, predecessorRegisters, h]
    · simp [commit]
    · simp [commit, beforeState, h]
    · rfl

/-- Compile a supplied ranked relation into an abstract semantic event log. -/
def log (rank : Event → ℕ)
    (rank_lt : ∀ {c d : Event}, R c d → rank c < rank d) :
    SemanticEventLog Event Unit Event where
  commitOf := commit R
  commitOf_eventId := fun _ ↦ rfl
  rank := rank
  rank_lt_of_parent := by
    intro c d h
    exact rank_lt ((authenticatedParent_commit_iff R c d).mp h)

variable (rank : Event → ℕ)
variable (rank_lt : ∀ {c d : Event}, R c d → rank c < rank d)

/-- Exact direct-edge receipt for the compiled semantic log. -/
theorem parentEdge_iff (c d : Event) :
    (log R rank rank_lt).ParentEdge c d ↔ R c d :=
  authenticatedParent_commit_iff R c d

/-- Without assuming transitivity, generated precedence is exactly the
transitive closure of the supplied relation. -/
theorem generatedBefore_iff_transGen (c d : Event) :
    (log R rank rank_lt).GeneratedBefore c d ↔
      Relation.TransGen R c d := by
  constructor
  · intro h
    induction h with
    | single hedge =>
        exact Relation.TransGen.single ((parentEdge_iff R rank rank_lt _ _).mp hedge)
    | tail _ hedge ih =>
        exact Relation.TransGen.tail ih
          ((parentEdge_iff R rank rank_lt _ _).mp hedge)
  · intro h
    induction h with
    | single hedge =>
        exact Relation.TransGen.single
          ((parentEdge_iff R rank rank_lt _ _).mpr hedge)
    | tail _ hedge ih =>
        exact Relation.TransGen.tail ih
          ((parentEdge_iff R rank rank_lt _ _).mpr hedge)

/-- If the supplied ranked relation is transitive, the complete generated
strict order is exactly that relation. -/
theorem generatedBefore_iff
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c)
    (c d : Event) :
    (log R rank rank_lt).GeneratedBefore c d ↔ R c d := by
  constructor
  · exact (log R rank rank_lt).generatedBefore_le_of_transitive R
      (fun a b h ↦ (parentEdge_iff R rank rank_lt a b).mp h) htrans
  · intro h
    exact (log R rank rank_lt).parentEdge_generatedBefore
      ((parentEdge_iff R rank rank_lt c d).mpr h)

/-- Reflexive-order version: a supplied finite ranked strict transitive
relation is realized exactly as authenticated equality-or-ancestry. -/
theorem generatedBeforeEq_iff
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c)
    (c d : Event) :
    (log R rank rank_lt).GeneratedBeforeEq c d ↔ c = d ∨ R c d := by
  simp only [SemanticEventLog.GeneratedBeforeEq]
  exact or_congr Iff.rfl (generatedBefore_iff R rank rank_lt htrans c d)

/-! ## Removing the supplied rank for finite strict transitive relations -/

/-- Intrinsic finite rank: the number of strict predecessors. -/
def predecessorRank (d : Event) : ℕ :=
  (predecessorRegisters R d).card

omit [DecidableEq Event] in
/-- In a strict transitive relation, predecessor sets grow strictly along
the relation, so predecessor count is an automatically derived rank. -/
theorem predecessorRank_lt_of_rel
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c)
    {c d : Event} (hcd : R c d) :
    predecessorRank R c < predecessorRank R d := by
  have hsub : predecessorRegisters R c ⊆ predecessorRegisters R d := by
    intro a ha
    have hac : R a c := by
      simpa [predecessorRegisters] using ha
    simpa [predecessorRegisters] using htrans hac hcd
  have hmem : c ∈ predecessorRegisters R d := by
    simp [predecessorRegisters, hcd]
  have hnot : c ∉ predecessorRegisters R c := by
    simp [predecessorRegisters, hirrefl c]
  exact Finset.card_lt_card
    ((Finset.ssubset_iff_of_subset hsub).2 ⟨c, hmem, hnot⟩)

/-- Every finite decidable strict transitive relation therefore compiles
without a supplied rank.  The rank field of the resulting abstract semantic
log is derived from predecessor cardinality. -/
def strictTransitiveLog
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c) :
    SemanticEventLog Event Unit Event :=
  log R (predecessorRank R)
    (fun h ↦ predecessorRank_lt_of_rel R hirrefl htrans h)

/-- Exact direct authenticated edges in the rank-free compiler. -/
theorem strictTransitiveLog_parentEdge_iff
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c)
    (c d : Event) :
    (strictTransitiveLog R hirrefl htrans).ParentEdge c d ↔ R c d :=
  parentEdge_iff R (predecessorRank R)
    (fun h ↦ predecessorRank_lt_of_rel R hirrefl htrans h) c d

/-- Exact generated strict order in the rank-free compiler. -/
theorem strictTransitiveLog_generatedBefore_iff
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c)
    (c d : Event) :
    (strictTransitiveLog R hirrefl htrans).GeneratedBefore c d ↔ R c d :=
  generatedBefore_iff R (predecessorRank R)
    (fun h ↦ predecessorRank_lt_of_rel R hirrefl htrans h) htrans c d

/-- Bundled exact-realization receipt.  Every finite decidable strict
transitive relation is the source-generated order of an authenticated
semantic log.  This proves grammar expressivity, not physical or dynamical
selection of the relation. -/
theorem finiteStrictTransitiveRelation_realized
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c) :
    ∃ L : SemanticEventLog Event Unit Event,
      (∀ c d, L.ParentEdge c d ↔ R c d) ∧
      (∀ c d, L.GeneratedBefore c d ↔ R c d) := by
  exact ⟨strictTransitiveLog R hirrefl htrans,
    strictTransitiveLog_parentEdge_iff R hirrefl htrans,
    strictTransitiveLog_generatedBefore_iff R hirrefl htrans⟩

/-! ## Axiom audit -/

#print axioms authenticatedParent_commit_iff
#print axioms log
#print axioms parentEdge_iff
#print axioms generatedBefore_iff_transGen
#print axioms generatedBefore_iff
#print axioms generatedBeforeEq_iff
#print axioms predecessorRank_lt_of_rel
#print axioms strictTransitiveLog
#print axioms strictTransitiveLog_parentEdge_iff
#print axioms strictTransitiveLog_generatedBefore_iff
#print axioms finiteStrictTransitiveRelation_realized

end OPH.Provenance.FiniteCausetCompiler
