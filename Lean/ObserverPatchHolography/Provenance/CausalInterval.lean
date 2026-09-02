import ObserverPatchHolography.Provenance.SemanticEventProvenance
import Mathlib.Data.Fin.VecNotation

/-!
# Source-derived causal intervals and the Boolean response diamond

For two events related by the generated precedence, the causal interval is
the set of events between them.  The interval is defined from provenance
alone: no coordinate, chart, clock, or manifold participates.  Calling its
boundary null would require a separate Lorentzian realization theorem, and
none is claimed here.

The committed four-event log realizes the minimal nontrivial interval: one
injection writes a seam register, two responses each read exactly that
register and write disjoint outputs, and one answer reads both outputs.
The generated order is the Boolean diamond: the responses are incomparable,
and the interval from injection to answer is the full four-element carrier.
The log also witnesses schedule blindness at the definitional level: the
pre-commit snapshot of the second response contains the first response's
output, yet no edge between the responses exists because the certified
causal support does not cite it.  A parent-child pair alone spans a
two-element interval, so a nontrivial diamond requires branching,
independent response, and semantic joining.
-/

namespace OPH.Provenance

universe u v w

variable {Register : Type u} {Value : Type v} {EventId : Type w}
variable [DecidableEq Register]

namespace SemanticEventLog

variable (L : SemanticEventLog Register Value EventId)

/-- The ancestry rank is monotone along the reflexive precedence. -/
theorem rank_le_of_generatedBeforeEq {a b : EventId}
    (h : L.GeneratedBeforeEq a b) : L.rank a ≤ L.rank b := by
  rcases h with rfl | h
  · exact le_refl _
  · exact le_of_lt (L.rank_lt_of_generatedBefore h)

/-- The source-derived causal interval between two events: every event
above the first and below the second in the reflexive generated
precedence. -/
def interval (a b : EventId) : Set EventId :=
  {x : EventId | L.GeneratedBeforeEq a x ∧ L.GeneratedBeforeEq x b}

/-- Every interval is finite at a finite event cutoff.  This is the local
finiteness axiom of an abstract causal set at that cutoff; it does not identify
the informational order with physical causality or establish a refinement
limit. -/
theorem interval_finite [Fintype EventId] (a b : EventId) :
    (L.interval a b).Finite :=
  Set.toFinite _

/-- At every finite cutoff, authenticated provenance therefore supplies the
order-theoretic causal-set axioms: a partial order with finite intervals.
The statement is deliberately informational; causal-set faithful embedding,
including calibrated count density, and manifoldlikeness remain separate
certificates. -/
theorem finiteCausalSetAxioms [Fintype EventId] :
    (∀ e, L.GeneratedBeforeEq e e) ∧
    (∀ a b c, L.GeneratedBeforeEq a b → L.GeneratedBeforeEq b c →
      L.GeneratedBeforeEq a c) ∧
    (∀ a b, L.GeneratedBeforeEq a b → L.GeneratedBeforeEq b a → a = b) ∧
    (∀ a b, (L.interval a b).Finite) :=
  ⟨L.generatedBeforeEq_refl,
    fun _ _ _ hab hbc ↦ L.generatedBeforeEq_trans hab hbc,
    fun _ _ hab hba ↦ L.generatedBeforeEq_antisymm hab hba,
    L.interval_finite⟩

theorem left_mem_interval {a b : EventId}
    (h : L.GeneratedBeforeEq a b) : a ∈ L.interval a b :=
  ⟨Or.inl rfl, h⟩

theorem right_mem_interval {a b : EventId}
    (h : L.GeneratedBeforeEq a b) : b ∈ L.interval a b :=
  ⟨h, Or.inl rfl⟩

/-! ## The interval-linearity obstruction

The following two predicates isolate the exact order-theoretic obstruction
seen in a chain-like event log.  They concern the source-generated
informational order only.  In particular, they make no statement about a
Hasse diagram being a physical causal set, about manifoldlikeness, or about a
continuum realization. -/

/-- Every two events in every source-derived interval are comparable.  This
is a global chain condition on generated intervals, not a dimension or
manifold criterion. -/
def IntervalLinear : Prop :=
  ∀ ⦃a b x y : EventId⦄,
    x ∈ L.interval a b → y ∈ L.interval a b →
      L.GeneratedBeforeEq x y ∨ L.GeneratedBeforeEq y x

/-- There is no incomparable fork-and-join witness: whenever `x` and `y`
both lie above `a` and below `b`, they are comparable.  The name records the
forbidden order pattern; it does not assert any graph-theoretic forest
presentation. -/
def NoForkJoin : Prop :=
  ∀ ⦃a b x y : EventId⦄,
    L.GeneratedBeforeEq a x → L.GeneratedBeforeEq a y →
    L.GeneratedBeforeEq x b → L.GeneratedBeforeEq y b →
      L.GeneratedBeforeEq x y ∨ L.GeneratedBeforeEq y x

/-- Interval linearity is exactly the absence of an incomparable
fork-and-join witness. -/
theorem intervalLinear_iff_noForkJoin :
    L.IntervalLinear ↔ L.NoForkJoin := by
  constructor
  · intro h a b x y hax hay hxb hyb
    exact h ⟨hax, hxb⟩ ⟨hay, hyb⟩
  · intro h a b x y hx hy
    exact h hx.1 hy.1 hx.2 hy.2

end SemanticEventLog

/-! ## The committed Boolean response diamond -/

namespace BooleanDiamond

/-- The initial snapshot: all four registers false, no writers. -/
def state0 : VersionedState (Fin 4) Bool (Fin 4) :=
  ⟨![false, false, false, false], ![none, none, none, none]⟩

/-- After the seam injection at register 0. -/
def state1 : VersionedState (Fin 4) Bool (Fin 4) :=
  ⟨![true, false, false, false], ![some 0, none, none, none]⟩

/-- After the first response's output at register 1. -/
def state2 : VersionedState (Fin 4) Bool (Fin 4) :=
  ⟨![true, true, false, false], ![some 0, some 1, none, none]⟩

/-- After the second response's output at register 2. -/
def state3 : VersionedState (Fin 4) Bool (Fin 4) :=
  ⟨![true, true, true, false], ![some 0, some 1, some 2, none]⟩

/-- After the answer's output at register 3. -/
def state4 : VersionedState (Fin 4) Bool (Fin 4) :=
  ⟨![true, true, true, true], ![some 0, some 1, some 2, some 3]⟩

/-- The seam injection: event 0 writes register 0 from the initial
snapshot and cites nothing. -/
def injection : SemanticCommit (Fin 4) Bool (Fin 4) where
  eventId := 0
  before := state0
  after := state1
  readSet := ∅
  writeSet := {0}
  causalSupp := ∅
  supp_subset_read := by decide
  frame_value := by decide
  frame_writer := by decide
  stamp := by decide

/-- The first response: event 1 cites exactly the injected register 0 and
writes its own output register 1. -/
def leftResponse : SemanticCommit (Fin 4) Bool (Fin 4) where
  eventId := 1
  before := state1
  after := state2
  readSet := {0}
  writeSet := {1}
  causalSupp := {0}
  supp_subset_read := by decide
  frame_value := by decide
  frame_writer := by decide
  stamp := by decide

/-- The second response: event 2 also cites exactly register 0 and writes
register 2.  Its pre-commit snapshot contains the first response's output,
but its certified support does not cite it, so the serialization leaves no
provenance edge between the responses. -/
def rightResponse : SemanticCommit (Fin 4) Bool (Fin 4) where
  eventId := 2
  before := state2
  after := state3
  readSet := {0}
  writeSet := {2}
  causalSupp := {0}
  supp_subset_read := by decide
  frame_value := by decide
  frame_writer := by decide
  stamp := by decide

/-- The answer: event 3 cites both response outputs and writes register
3. -/
def answer : SemanticCommit (Fin 4) Bool (Fin 4) where
  eventId := 3
  before := state3
  after := state4
  readSet := {1, 2}
  writeSet := {3}
  causalSupp := {1, 2}
  supp_subset_read := by decide
  frame_value := by decide
  frame_writer := by decide
  stamp := by decide

/-- The four-event diamond log with ancestry ranks 0, 1, 1, 2. -/
def log : SemanticEventLog (Fin 4) Bool (Fin 4) where
  commitOf := ![injection, leftResponse, rightResponse, answer]
  commitOf_eventId := by decide
  rank := ![0, 1, 1, 2]
  rank_lt_of_parent := by
    have h : ∀ c d : Fin 4,
        AuthenticatedDirectSemanticParent
            (![injection, leftResponse, rightResponse, answer] c)
            (![injection, leftResponse, rightResponse, answer] d) →
          (![0, 1, 1, 2] : Fin 4 → ℕ) c < ![0, 1, 1, 2] d := by decide
    intro c d
    exact h c d

/-- The direct parent edges are exactly injection-to-responses and
responses-to-answer. -/
theorem parentEdge_iff (c d : Fin 4) :
    log.ParentEdge c d ↔
      (c = 0 ∧ d = 1) ∨ (c = 0 ∧ d = 2) ∨
        (c = 1 ∧ d = 3) ∨ (c = 2 ∧ d = 3) := by
  revert c d
  decide

theorem before_left : log.GeneratedBefore 0 1 :=
  Relation.TransGen.single (by decide)

theorem before_right : log.GeneratedBefore 0 2 :=
  Relation.TransGen.single (by decide)

theorem left_before_answer : log.GeneratedBefore 1 3 :=
  Relation.TransGen.single (by decide)

theorem right_before_answer : log.GeneratedBefore 2 3 :=
  Relation.TransGen.single (by decide)

theorem injection_before_answer : log.GeneratedBefore 0 3 :=
  log.generatedBefore_trans before_left left_before_answer

/-- The two responses are incomparable: the diamond genuinely branches. -/
theorem responses_incomparable :
    ¬ log.GeneratedBefore 1 2 ∧ ¬ log.GeneratedBefore 2 1 := by
  constructor
  · intro h
    exact absurd (log.rank_lt_of_generatedBefore h) (by decide)
  · intro h
    exact absurd (log.rank_lt_of_generatedBefore h) (by decide)

/-- The interval from injection to answer is the full four-element
carrier: the Boolean response diamond. -/
theorem interval_injection_answer_eq_univ :
    log.interval 0 3 = Set.univ := by
  apply Set.eq_univ_of_forall
  intro x
  fin_cases x
  · exact ⟨Or.inl rfl, Or.inr injection_before_answer⟩
  · exact ⟨Or.inr before_left, Or.inr left_before_answer⟩
  · exact ⟨Or.inr before_right, Or.inr right_before_answer⟩
  · exact ⟨Or.inr injection_before_answer, Or.inl rfl⟩

/-- A parent-child pair alone spans a two-element interval: the interval
from the injection to the first response contains exactly the two
endpoints. -/
theorem interval_injection_left_eq_pair :
    log.interval 0 1 = {0, 1} := by
  ext x
  constructor
  · rintro ⟨_, hupper⟩
    rcases hupper with rfl | hlt
    · exact Or.inr rfl
    · left
      have hrank := log.rank_lt_of_generatedBefore hlt
      fin_cases x
      · rfl
      · exact absurd hrank (by decide)
      · exact absurd hrank (by decide)
      · exact absurd hrank (by decide)
  · rintro (rfl | rfl)
    · exact log.left_mem_interval (Or.inr before_left)
    · exact log.right_mem_interval (Or.inr before_left)

/-- The Boolean response diamond is a concrete counterexample to interval
linearity: its two response events lie in one interval and are incomparable.
-/
theorem not_intervalLinear : ¬ log.IntervalLinear := by
  intro h
  have hcomp : log.GeneratedBeforeEq 1 2 ∨
      log.GeneratedBeforeEq 2 1 :=
    h (a := 0) (b := 3) (x := 1) (y := 2)
      ⟨Or.inr before_left, Or.inr left_before_answer⟩
      ⟨Or.inr before_right, Or.inr right_before_answer⟩
  rcases hcomp with h12 | h21
  · exact responses_incomparable.1 (by
      rcases h12 with hEq | hBefore
      · omega
      · exact hBefore)
  · exact responses_incomparable.2 (by
      rcases h21 with hEq | hBefore
      · omega
      · exact hBefore)

/-- Equivalently, the Boolean response diamond contains an incomparable
fork-and-join witness. -/
theorem not_noForkJoin : ¬ log.NoForkJoin := by
  intro h
  exact not_intervalLinear
    ((log.intervalLinear_iff_noForkJoin).mpr h)

end BooleanDiamond

#print axioms SemanticEventLog.interval
#print axioms SemanticEventLog.finiteCausalSetAxioms
#print axioms SemanticEventLog.intervalLinear_iff_noForkJoin
#print axioms BooleanDiamond.log
#print axioms BooleanDiamond.parentEdge_iff
#print axioms BooleanDiamond.responses_incomparable
#print axioms BooleanDiamond.interval_injection_answer_eq_univ
#print axioms BooleanDiamond.interval_injection_left_eq_pair
#print axioms BooleanDiamond.not_intervalLinear
#print axioms BooleanDiamond.not_noForkJoin

end OPH.Provenance
