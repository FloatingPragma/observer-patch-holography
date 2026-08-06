import Tower.FixedPointEndpoint

/-!
# A full-premise inhabitant for the A4 endpoint packet

This module discharges the existence half of inherited obligation 2 of
completion-plan issue `#693`.  The closed A4 packet (issue `#699`) proves a
conditional endpoint theorem: termination, confluence, semantic fixed-point
completeness, repair-output congruence, and enabledness congruence together
force one completed public endpoint across schedules and public-equivalent
representatives.  The A4 closure records that no exhibited source branch
combines every premise.  This module exhibits one explicit finite branch
satisfying every premise of `public_endpoint_exists_unique_on_public_class`
simultaneously and consumes that theorem to state the unconditional
corollary, a statement with no remaining hypothesis.

The branch is nontrivial in the three audited directions.  Every public
class has two distinct raw representatives because the readback erases one
hidden bit (`two_raws_per_class`, `hidden_collapse`).  The rewrite system
has a genuine choice point: at the doubly broken start both labelled moves
are enabled and produce distinct successors (`choice_point_exists`), so
confluence is exercised on an actual peak rather than holding vacuously
(`confluence_nonvacuous`); the peak is closed by an explicit commutation
diamond and Newman's lemma.  The corollary
`branch_endpoint_unconditional` then states, premise-free, that every
start completes to one public endpoint independent of the schedule and of
the raw representative of its public class.

## Claim boundary

This is the existence half of obligation 2 only: a branch combining every
A4 premise exists as an explicit finite mathematical object.  The
source-selection half remains open: no theorem here states that the OPH
source produces this branch, and the branch is synthetic rather than
source-derived.  No physical world, clock, continuum limit, or
cross-regulator naturality is inferred.
-/

namespace OPH.Tower

namespace FullPremiseInhabitant

open PublicWorldPresentation
open PublicWorldPresentation.FiniteRepairSystem
open Relation

/-- Raw configurations: a public two-flag pair and one hidden bit. -/
abbrev BranchConfig := (Bool × Bool) × Bool

/-- The branch presentation.  The readback publishes the two flags and
erases the hidden bit, so every public class has exactly two raw
representatives. -/
@[reducible] def branchPresentation : PublicWorldPresentation (Bool × Bool) where
  Config := BranchConfig
  configFintype := inferInstance
  configNonempty := ⟨((false, false), false)⟩
  readback := Prod.fst

/-- Two labelled repair moves, one per public flag. -/
inductive BranchMove
  | fixLeft | fixRight
  deriving DecidableEq, Fintype

/-- Each move sets its own flag to `true`, keeps the other flag, and
transports the hidden bit unchanged. -/
def branchRepair : BranchMove → BranchConfig → BranchConfig
  | .fixLeft, ((_, b), h) => ((true, b), h)
  | .fixRight, ((a, _), h) => ((a, true), h)

/-- The labelled repair system of the branch. -/
@[reducible] def branchSystem : branchPresentation.FiniteRepairSystem where
  Move := BranchMove
  moveFintype := inferInstance
  repair := branchRepair

/-- Declared semantic consistency: both public flags read `true`. -/
def branchConsistent (x : BranchConfig) : Prop := x.1 = (true, true)

/-- `fixLeft` is enabled exactly when the first public flag is `false`; the
hidden bit never obstructs quiescence. -/
theorem fixLeft_enabled_iff (x : BranchConfig) :
    branchRepair BranchMove.fixLeft x ≠ x ↔ x.1.1 = false := by
  rcases x with ⟨⟨a, b⟩, h⟩
  cases a <;> simp [branchRepair]

/-- `fixRight` is enabled exactly when the second public flag is `false`. -/
theorem fixRight_enabled_iff (x : BranchConfig) :
    branchRepair BranchMove.fixRight x ≠ x ↔ x.1.2 = false := by
  rcases x with ⟨⟨a, b⟩, h⟩
  cases b <;> simp [branchRepair]

/-- The public effect of one move as a function of the readback alone. -/
def publicStep : BranchMove → Bool × Bool → Bool × Bool
  | .fixLeft, p => (true, p.2)
  | .fixRight, p => (p.1, true)

/-- The raw repair covers its public effect through the readback. -/
theorem readback_branchRepair (m : BranchMove) (x : BranchConfig) :
    (branchRepair m x).1 = publicStep m x.1 := by
  rcases x with ⟨⟨a, b⟩, h⟩
  cases m <;> rfl

/-- Premises four and five, the exact A4 congruence packet: the public
effect and the enabledness of every labelled move depend only on the
readback. -/
theorem branchQuotientCompatible :
    QuotientCompatible branchPresentation branchSystem :=
  { repair_respects_public := by
      intro m x y hxy
      show (branchRepair m x).1 = (branchRepair m y).1
      rw [readback_branchRepair, readback_branchRepair]
      exact congrArg (publicStep m) hxy
    enabled_respects_public := by
      intro m x y hxy
      have hpub : x.1 = y.1 := hxy
      cases m
      · show branchRepair BranchMove.fixLeft x ≠ x ↔
          branchRepair BranchMove.fixLeft y ≠ y
        rw [fixLeft_enabled_iff, fixLeft_enabled_iff, hpub]
      · show branchRepair BranchMove.fixRight x ≠ x ↔
          branchRepair BranchMove.fixRight y ≠ y
        rw [fixRight_enabled_iff, fixRight_enabled_iff, hpub] }

/-- Mismatch measure: the number of public flags reading `false`. -/
def branchMeasure (x : BranchConfig) : ℕ :=
  (if x.1.1 then 0 else 1) + (if x.1.2 then 0 else 1)

/-- Every accepted step strictly decreases the mismatch measure. -/
theorem branch_measure_desc {x y : BranchConfig}
    (h : step branchPresentation branchSystem x y) :
    branchMeasure y < branchMeasure x := by
  obtain ⟨m, rfl, hne⟩ := h
  rcases x with ⟨⟨a, b⟩, hb⟩
  cases m <;> cases a <;> cases b <;>
    simp_all [branchRepair, branchMeasure]

/-- Premise one: termination. -/
theorem branch_terminating :
    OPH.AbstractRewriting.Terminating (step branchPresentation branchSystem) :=
  terminating_of_measure branchPresentation branchSystem
    branchMeasure branch_measure_desc

/-- One labelled move, taken when enabled and skipped when quiescent, gives
a reflexive-transitive accepted reduction. -/
theorem reflTransGen_repair (m : BranchMove) (x : BranchConfig) :
    ReflTransGen (step branchPresentation branchSystem) x (branchRepair m x) := by
  by_cases hm : branchRepair m x = x
  · rw [hm]
  · exact ReflTransGen.single ⟨m, rfl, hm⟩

/-- The two moves commute on raw configurations. -/
theorem branchRepair_comm (x : BranchConfig) :
    branchRepair BranchMove.fixRight (branchRepair BranchMove.fixLeft x) =
      branchRepair BranchMove.fixLeft (branchRepair BranchMove.fixRight x) := by
  rcases x with ⟨⟨a, b⟩, h⟩
  rfl

/-- Direct diamond argument on the explicit rewrite system: equal moves join
immediately and distinct moves join through the commuted double repair.  The
mixed cases are the genuine peaks; they are the reason this proof is a
diamond argument rather than a determinism remark. -/
theorem branch_locallyConfluent :
    OPH.AbstractRewriting.LocallyConfluent
      (step branchPresentation branchSystem) := by
  rintro x y z ⟨m, rfl, _⟩ ⟨n, rfl, _⟩
  cases m <;> cases n
  · exact ⟨branchRepair BranchMove.fixLeft x,
      ReflTransGen.refl, ReflTransGen.refl⟩
  · refine ⟨branchRepair BranchMove.fixRight
        (branchRepair BranchMove.fixLeft x),
      reflTransGen_repair BranchMove.fixRight
        (branchRepair BranchMove.fixLeft x), ?_⟩
    rw [branchRepair_comm]
    exact reflTransGen_repair BranchMove.fixLeft
      (branchRepair BranchMove.fixRight x)
  · refine ⟨branchRepair BranchMove.fixRight
        (branchRepair BranchMove.fixLeft x), ?_,
      reflTransGen_repair BranchMove.fixRight
        (branchRepair BranchMove.fixLeft x)⟩
    rw [branchRepair_comm]
    exact reflTransGen_repair BranchMove.fixLeft
      (branchRepair BranchMove.fixRight x)
  · exact ⟨branchRepair BranchMove.fixRight x,
      ReflTransGen.refl, ReflTransGen.refl⟩

/-- Premise two: confluence, by Newman's lemma from termination and the
diamond argument. -/
theorem branch_confluent :
    OPH.AbstractRewriting.Confluent (step branchPresentation branchSystem) :=
  OPH.AbstractRewriting.newman_lemma _
    branch_terminating branch_locallyConfluent

/-- Premise three: semantic completeness.  Simultaneous repair fixed points
are exactly the declared consistent configurations. -/
theorem branch_completeFor :
    CompleteFor branchPresentation branchSystem branchConsistent := by
  intro x
  rcases x with ⟨⟨a, b⟩, h⟩
  constructor
  · intro hfix
    have hl : branchRepair BranchMove.fixLeft ((a, b), h) = ((a, b), h) :=
      hfix BranchMove.fixLeft
    have hr : branchRepair BranchMove.fixRight ((a, b), h) = ((a, b), h) :=
      hfix BranchMove.fixRight
    have ha : true = a := congrArg (fun c => c.1.1) hl
    have hb : true = b := congrArg (fun c => c.1.2) hr
    show ((a, b) : Bool × Bool) = (true, true)
    rw [← ha, ← hb]
  · intro hcons m
    have hab : a = true ∧ b = true := by
      have hpair : ((a, b) : Bool × Bool) = (true, true) := hcons
      simpa [Prod.ext_iff] using hpair
    obtain ⟨ha, hb⟩ := hab
    subst ha
    subst hb
    cases m <;> rfl

/-! ## Receipts: the premises are exercised on nontrivial content -/

/-- The doubly broken start with cleared hidden bit. -/
def choiceStart : BranchConfig := ((false, false), false)

/-- Receipt one: a genuine choice point.  At `choiceStart` two distinct
moves are enabled and produce distinct successors. -/
theorem choice_point_exists :
    branchRepair BranchMove.fixLeft choiceStart ≠ choiceStart ∧
      branchRepair BranchMove.fixRight choiceStart ≠ choiceStart ∧
      branchRepair BranchMove.fixLeft choiceStart ≠
        branchRepair BranchMove.fixRight choiceStart := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

/-- Receipt two: confluence is exercised on an actual peak.  The two
enabled moves at `choiceStart` give two distinct one-step successors, and
both reduce to the common repaired configuration. -/
theorem confluence_nonvacuous :
    step branchPresentation branchSystem choiceStart ((true, false), false) ∧
      step branchPresentation branchSystem choiceStart
        ((false, true), false) ∧
      (((true, false), false) : BranchConfig) ≠ ((false, true), false) ∧
      ReflTransGen (step branchPresentation branchSystem)
        ((true, false), false) ((true, true), false) ∧
      ReflTransGen (step branchPresentation branchSystem)
        ((false, true), false) ((true, true), false) := by
  refine ⟨⟨BranchMove.fixLeft, rfl, by decide⟩,
    ⟨BranchMove.fixRight, rfl, by decide⟩, by decide, ?_, ?_⟩
  · exact ReflTransGen.single ⟨BranchMove.fixRight, rfl, by decide⟩
  · exact ReflTransGen.single ⟨BranchMove.fixLeft, rfl, by decide⟩

/-- Fixed completing schedule used by the collapse receipt. -/
def collapseSchedule : List BranchMove :=
  [BranchMove.fixLeft, BranchMove.fixRight]

/-- Receipt three: hidden data collapses.  Two raw configurations differing
only in the hidden bit are distinct, share one public class, complete under
the same schedule, and land on one public endpoint. -/
theorem hidden_collapse :
    (((false, false), false) : BranchConfig) ≠ ((false, false), true) ∧
      branchPresentation.toPublicWorld ((false, false), false) =
        branchPresentation.toPublicWorld ((false, false), true) ∧
      CompletedSchedule branchPresentation branchSystem collapseSchedule
        ((false, false), false) ∧
      CompletedSchedule branchPresentation branchSystem collapseSchedule
        ((false, false), true) ∧
      branchPresentation.toPublicWorld
          (run branchPresentation branchSystem collapseSchedule
            ((false, false), false)) =
        branchPresentation.toPublicWorld
          (run branchPresentation branchSystem collapseSchedule
            ((false, false), true)) := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · intro hcontra
    exact Bool.false_ne_true (congrArg Prod.snd hcontra)
  · exact (branchPresentation.toPublicWorld_eq_iff _ _).2 rfl
  · intro m; cases m <;> rfl
  · intro m; cases m <;> rfl
  · exact (branchPresentation.toPublicWorld_eq_iff _ _).2 rfl

/-- Every public class of the branch has two distinct raw representatives:
the readback erases the hidden bit on every class, so the quotient is
nontrivial everywhere rather than at one exhibited point. -/
theorem two_raws_per_class (w : branchPresentation.PublicWorld) :
    ∃ x y : BranchConfig, x ≠ y ∧
      branchPresentation.toPublicWorld x = w ∧
      branchPresentation.toPublicWorld y = w := by
  refine Quotient.inductionOn w ?_
  intro z
  refine ⟨(z.1, false), (z.1, true), ?_, ?_, ?_⟩
  · intro hcontra
    exact Bool.false_ne_true (congrArg Prod.snd hcontra)
  · exact (branchPresentation.toPublicWorld_eq_iff _ _).2 rfl
  · exact (branchPresentation.toPublicWorld_eq_iff _ _).2 rfl

/-! ## The unconditional endpoint corollary -/

/-- The A4 endpoint theorem consumed with every premise discharged.  The
statement carries no hypothesis: every start has a completed schedule, the
completed endpoint is consistent, and every completed schedule from every
raw representative of the same public class lands on one public endpoint.
Schedule independence and representative independence hold jointly.  This
theorem is the existence half of E2 obligation 2; whether the OPH source
produces this branch is a separate open question outside this module. -/
theorem branch_endpoint_unconditional (start : BranchConfig) :
    ∃ (world : branchPresentation.PublicWorld) (schedule : List BranchMove),
      CompletedSchedule branchPresentation branchSystem schedule start ∧
      branchConsistent (run branchPresentation branchSystem schedule start) ∧
      world = branchPresentation.toPublicWorld
        (run branchPresentation branchSystem schedule start) ∧
      ∀ otherStart : BranchConfig,
        branchPresentation.toPublicWorld otherStart =
          branchPresentation.toPublicWorld start →
        ∀ other : List BranchMove,
          CompletedSchedule branchPresentation branchSystem other otherStart →
          branchPresentation.toPublicWorld
              (run branchPresentation branchSystem other otherStart) =
            world :=
  public_endpoint_exists_unique_on_public_class
    branchPresentation branchSystem
    branch_terminating branch_confluent branchQuotientCompatible
    branchConsistent branch_completeFor start

end FullPremiseInhabitant

-- Axiom audit: only Lean/Mathlib foundations are permitted here.
#print axioms OPH.Tower.FullPremiseInhabitant.branch_terminating
#print axioms OPH.Tower.FullPremiseInhabitant.branch_confluent
#print axioms OPH.Tower.FullPremiseInhabitant.branchQuotientCompatible
#print axioms OPH.Tower.FullPremiseInhabitant.branch_completeFor
#print axioms OPH.Tower.FullPremiseInhabitant.choice_point_exists
#print axioms OPH.Tower.FullPremiseInhabitant.confluence_nonvacuous
#print axioms OPH.Tower.FullPremiseInhabitant.hidden_collapse
#print axioms OPH.Tower.FullPremiseInhabitant.two_raws_per_class
#print axioms OPH.Tower.FullPremiseInhabitant.branch_endpoint_unconditional

end OPH.Tower
