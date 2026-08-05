import Tower.PublicWorldQuotient
import ObserverPatchHolography.Primitives

/-!
# Finite fixed-point endpoints for labelled repair

This module implements the endpoint half of completion-plan issue `#699`.
A `FiniteRepairSystem` consists of finitely many labelled repair maps on the
raw configurations of a `PublicWorldPresentation`.  Its accepted step relation
records only moves that actually change the configuration.

The hypotheses and their jobs are kept separate:

* termination gives a finite run reaching some fixed configuration;
* confluence makes fixed endpoints from the same start schedule-independent;
* `CompletedSchedule` says that a finite schedule has not stopped while any
  labelled repair remains enabled.

This is a terminal completion postcondition, not an infinite scheduler-
fairness claim.
Repair labels and schedules are proof witnesses.  They are not a physical
time, a worldline, or part of public-world identity.

Three finite countermodels locate necessity precisely.  The fork system is
terminating and has completed schedules but, without confluence, reaches two
different public endpoints.  The settling system is terminating and
confluent, but an incomplete truncated schedule stops before its fixed endpoint.
The representative system is terminating and confluent but cannot descend to
the quotient because a hidden representative bit controls its public endpoint.
-/

namespace OPH.Tower

open Relation

universe u v

namespace PublicWorldPresentation

variable {Signature : Type v} (P : PublicWorldPresentation Signature)

/-- A finite family of labelled repair maps on one raw presentation. -/
structure FiniteRepairSystem where
  Move : Type u
  moveFintype : Fintype Move
  repair : Move → P.Config → P.Config

attribute [instance] FiniteRepairSystem.moveFintype

namespace FiniteRepairSystem

variable (R : P.FiniteRepairSystem)

/-- An accepted step is a labelled repair that changes the configuration. -/
def step (x y : P.Config) : Prop :=
  ∃ m : R.Move, y = R.repair m x ∧ R.repair m x ≠ x

/-- A configuration is fixed when every labelled repair is quiescent. -/
def IsFixed (x : P.Config) : Prop := ∀ m : R.Move, R.repair m x = x

/-- Execute a finite repair schedule from left to right.  Stuttering labels
are allowed in a schedule but do not count as accepted steps. -/
def run : List R.Move → P.Config → P.Config
  | [], x => x
  | m :: ms, x => run ms (R.repair m x)

/-- Terminal completion: after the schedule ends, no labelled repair remains
enabled.  This is deliberately not called scheduler fairness. -/
def CompletedSchedule (schedule : List R.Move) (start : P.Config) : Prop :=
  IsFixed P R (run P R schedule start)

/-- An explicit semantic completeness interface.  The repair system is
complete for `consistent` exactly when simultaneous repair fixed points are
the declared consistent configurations.  This is never inferred from
termination or confluence. -/
def CompleteFor (consistent : P.Config → Prop) : Prop :=
  ∀ x : P.Config, IsFixed P R x ↔ consistent x

/-- Exact hypotheses needed for labelled repair to act coherently on public
quotient classes.  Output congruence makes every labelled move descend to the
quotient.  Enabledness congruence is separate and stronger: it prevents a
move from being quiescent for one representative but enabled for another. -/
structure QuotientCompatible : Prop where
  repair_respects_public : ∀ (m : R.Move) {x y : P.Config},
    P.readback x = P.readback y →
      P.readback (R.repair m x) = P.readback (R.repair m y)
  enabled_respects_public : ∀ (m : R.Move) {x y : P.Config},
    P.readback x = P.readback y →
      (R.repair m x ≠ x ↔ R.repair m y ≠ y)

/-- A quotient-compatible labelled move descends to a literal endomorphism
of the public-world quotient. -/
noncomputable def publicMove (hcompat : QuotientCompatible P R) (m : R.Move) :
    P.PublicWorld → P.PublicWorld :=
  Quotient.lift
    (fun x => P.toPublicWorld (R.repair m x))
    (by
      intro x y hxy
      exact (P.toPublicWorld_eq_iff _ _).2
        (hcompat.repair_respects_public m hxy))

@[simp]
theorem publicMove_mk (hcompat : QuotientCompatible P R) (m : R.Move)
    (x : P.Config) :
    publicMove P R hcompat m (P.toPublicWorld x) =
      P.toPublicWorld (R.repair m x) :=
  rfl

/-- Applying the same finite schedule to public-equivalent representatives
keeps the endpoints public-equivalent. -/
theorem run_respects_public (hcompat : QuotientCompatible P R)
    (schedule : List R.Move) {x y : P.Config}
    (hxy : P.readback x = P.readback y) :
    P.readback (run P R schedule x) = P.readback (run P R schedule y) := by
  induction schedule generalizing x y with
  | nil => exact hxy
  | cons m ms ih =>
      exact ih (hcompat.repair_respects_public m hxy)

/-- Completion is representative-independent when enabledness, not merely
the move output, respects public equivalence. -/
theorem completedSchedule_congr (hcompat : QuotientCompatible P R)
    (schedule : List R.Move) {x y : P.Config}
    (hxy : P.readback x = P.readback y) :
    CompletedSchedule P R schedule x ↔ CompletedSchedule P R schedule y := by
  have hend := run_respects_public P R hcompat schedule hxy
  constructor
  · intro hx m
    by_contra hy
    exact (hcompat.enabled_respects_public m hend).2 hy (hx m)
  · intro hy m
    by_contra hx
    exact (hcompat.enabled_respects_public m hend).1 hx (hy m)

/-- Under an explicit completeness receipt, schedule completion is exactly
consistency of the resulting configuration. -/
theorem completedSchedule_iff_consistent
    (consistent : P.Config → Prop) (hcomplete : CompleteFor P R consistent)
    (schedule : List R.Move) (start : P.Config) :
    CompletedSchedule P R schedule start ↔
      consistent (run P R schedule start) :=
  hcomplete _

/-- The accepted-step normal forms are exactly the simultaneous fixed points
of all labelled repairs. -/
theorem normalForm_iff_fixed (x : P.Config) :
    OPH.AbstractRewriting.IsNormalForm (step P R) x ↔ IsFixed P R x := by
  constructor
  · intro hnf m
    by_contra hne
    exact hnf (R.repair m x) ⟨m, rfl, hne⟩
  · intro hfix y hstep
    obtain ⟨m, _, hne⟩ := hstep
    exact hne (hfix m)

/-- Schedule completion is equivalently normality of its endpoint. -/
theorem completedSchedule_iff_normalForm (schedule : List R.Move)
    (start : P.Config) :
    CompletedSchedule P R schedule start ↔
      OPH.AbstractRewriting.IsNormalForm (step P R) (run P R schedule start) :=
  (normalForm_iff_fixed P R _).symm

/-- Running concatenated schedules is sequential composition. -/
theorem run_append (first second : List R.Move) (start : P.Config) :
    run P R (first ++ second) start =
      run P R second (run P R first start) := by
  induction first generalizing start with
  | nil => rfl
  | cons m ms ih =>
      simp only [List.cons_append, run]
      exact ih (R.repair m start)

/-- Every finite labelled schedule gives a reflexive-transitive accepted-step
reduction after its stuttering labels are erased. -/
theorem run_reachable (schedule : List R.Move) (start : P.Config) :
    ReflTransGen (step P R) start (run P R schedule start) := by
  induction schedule generalizing start with
  | nil => exact ReflTransGen.refl
  | cons m ms ih =>
      by_cases hmove : R.repair m start = start
      · simpa only [run, hmove] using ih start
      · exact ReflTransGen.head ⟨m, rfl, hmove⟩ (ih (R.repair m start))

/-- Every finite accepted reduction has a labelled schedule witness. -/
theorem reachable_has_schedule {start finish : P.Config}
    (h : ReflTransGen (step P R) start finish) :
    ∃ schedule : List R.Move, run P R schedule start = finish := by
  induction h with
  | refl => exact ⟨[], rfl⟩
  | tail hreach hstep ih =>
      obtain ⟨schedule, hschedule⟩ := ih
      obtain ⟨m, rfl, _⟩ := hstep
      refine ⟨schedule ++ [m], ?_⟩
      rw [run_append P R, hschedule]
      rfl

/-- A strictly decreasing natural-valued measure proves termination. -/
theorem terminating_of_measure (measure : P.Config → ℕ)
    (hdesc : ∀ {x y : P.Config}, step P R x y → measure y < measure x) :
    OPH.AbstractRewriting.Terminating (step P R) := by
  exact Subrelation.wf
    (fun {_ _} hxy => hdesc hxy)
    (InvImage.wf _ wellFounded_lt)

/-- Termination gives a reachable normal form from every raw configuration. -/
theorem exists_reachable_normal
    (hterm : OPH.AbstractRewriting.Terminating (step P R)) (start : P.Config) :
    ∃ finish : P.Config,
      ReflTransGen (step P R) start finish ∧
        OPH.AbstractRewriting.IsNormalForm (step P R) finish := by
  refine hterm.induction
    (C := fun x => ∃ finish : P.Config,
      ReflTransGen (step P R) x finish ∧
        OPH.AbstractRewriting.IsNormalForm (step P R) finish)
    start ?_
  intro x ih
  by_cases hnext : ∃ y : P.Config, step P R x y
  · obtain ⟨y, hxy⟩ := hnext
    obtain ⟨finish, hyf, hnormal⟩ := ih y hxy
    exact ⟨finish, ReflTransGen.head hxy hyf, hnormal⟩
  · exact ⟨x, ReflTransGen.refl, fun y hxy => hnext ⟨y, hxy⟩⟩

/-- Termination supplies an actually completed finite schedule from every
start; completion is not postulated as an unreachable oracle. -/
theorem exists_completedSchedule
    (hterm : OPH.AbstractRewriting.Terminating (step P R)) (start : P.Config) :
    ∃ schedule : List R.Move, CompletedSchedule P R schedule start := by
  obtain ⟨finish, hreach, hnormal⟩ :=
    exists_reachable_normal P R hterm start
  obtain ⟨schedule, hschedule⟩ := reachable_has_schedule P R hreach
  refine ⟨schedule, ?_⟩
  rw [CompletedSchedule, hschedule]
  exact (normalForm_iff_fixed P R finish).1 hnormal

/-- Raw configurations fixed by every labelled repair. -/
def FixedConfiguration := {x : P.Config // IsFixed P R x}

/-- Public fixed-point object: public classes that possess a fixed raw
representative.  It can contain more than one element when distinct basins
settle to distinct public signatures. -/
def FixedPointObject :=
  {w : P.PublicWorld //
    ∃ x : FixedConfiguration P R, P.toPublicWorld x.1 = w}

/-- The fixed-point object is inhabited whenever the repair relation
terminates and the raw presentation is inhabited. -/
theorem fixedPointObject_nonempty
    (hterm : OPH.AbstractRewriting.Terminating (step P R)) :
    Nonempty (FixedPointObject P R) := by
  obtain ⟨start⟩ := P.configNonempty
  obtain ⟨schedule, hcompleted⟩ :=
    exists_completedSchedule P R hterm start
  let finish := run P R schedule start
  exact ⟨⟨P.toPublicWorld finish, ⟨⟨finish, hcompleted⟩, rfl⟩⟩⟩

/-- Confluence plus terminal completion makes two finite schedules from the
same start agree even as raw configurations.  The public quotient conclusion
is therefore not obtained by hiding a raw schedule dependence. -/
theorem completed_runs_equal
    (hconf : OPH.AbstractRewriting.Confluent (step P R))
    {start : P.Config} {first second : List R.Move}
    (hfirst : CompletedSchedule P R first start)
    (hsecond : CompletedSchedule P R second start) :
    run P R first start = run P R second start := by
  apply OPH.AbstractRewriting.unique_normal_form (step P R) hconf
  · exact ⟨run_reachable P R first start,
      (completedSchedule_iff_normalForm P R first start).1 hfirst⟩
  · exact ⟨run_reachable P R second start,
      (completedSchedule_iff_normalForm P R second start).1 hsecond⟩

/-- The corresponding public-world classes are schedule-independent. -/
theorem completed_public_endpoints_equal
    (hconf : OPH.AbstractRewriting.Confluent (step P R))
    {start : P.Config} {first second : List R.Move}
    (hfirst : CompletedSchedule P R first start)
    (hsecond : CompletedSchedule P R second start) :
    P.toPublicWorld (run P R first start) =
      P.toPublicWorld (run P R second start) := by
  exact congrArg P.toPublicWorld
    (completed_runs_equal P R hconf hfirst hsecond)

/-- Strong public schedule independence: completed runs may start from
different raw representatives of the same public class.  Quotient
compatibility transports the first schedule to the second representative;
confluence then compares the two completed runs from that common start. -/
theorem completed_public_endpoints_equal_of_equiv
    (hconf : OPH.AbstractRewriting.Confluent (step P R))
    (hcompat : QuotientCompatible P R)
    {firstStart secondStart : P.Config}
    (hstart : P.readback firstStart = P.readback secondStart)
    {first second : List R.Move}
    (hfirst : CompletedSchedule P R first firstStart)
    (hsecond : CompletedSchedule P R second secondStart) :
    P.toPublicWorld (run P R first firstStart) =
      P.toPublicWorld (run P R second secondStart) := by
  have hfirstAtSecond : CompletedSchedule P R first secondStart :=
    (completedSchedule_congr P R hcompat first hstart).1 hfirst
  have hrun := run_respects_public P R hcompat first hstart
  calc
    P.toPublicWorld (run P R first firstStart) =
        P.toPublicWorld (run P R first secondStart) :=
      (P.toPublicWorld_eq_iff _ _).2 hrun
    _ = P.toPublicWorld (run P R second secondStart) :=
      completed_public_endpoints_equal P R hconf hfirstAtSecond hsecond

/-- Combined endpoint theorem. Termination produces a completed finite witness;
confluence and completion prove that every other completed schedule from the
same start has the same raw and hence public endpoint. -/
theorem public_endpoint_exists_unique
    (hterm : OPH.AbstractRewriting.Terminating (step P R))
    (hconf : OPH.AbstractRewriting.Confluent (step P R))
    (start : P.Config) :
    ∃ (world : P.PublicWorld) (schedule : List R.Move),
      CompletedSchedule P R schedule start ∧
      world = P.toPublicWorld (run P R schedule start) ∧
      ∀ other : List R.Move, CompletedSchedule P R other start →
        run P R other start = run P R schedule start ∧
        P.toPublicWorld (run P R other start) = world := by
  obtain ⟨schedule, hcompleted⟩ :=
    exists_completedSchedule P R hterm start
  refine ⟨P.toPublicWorld (run P R schedule start), schedule, hcompleted, rfl, ?_⟩
  intro other hother
  have heq := completed_runs_equal P R hconf hother hcompleted
  exact ⟨heq, congrArg P.toPublicWorld heq⟩

/-- Full bounded endpoint packet on one public input class.  Termination
constructs a completed schedule; completeness identifies its fixed endpoint
with the caller's semantic consistency predicate; confluence and quotient
compatibility make every completed schedule from every representative of the
same public class land in the same public world. -/
theorem public_endpoint_exists_unique_on_public_class
    (hterm : OPH.AbstractRewriting.Terminating (step P R))
    (hconf : OPH.AbstractRewriting.Confluent (step P R))
    (hcompat : QuotientCompatible P R)
    (consistent : P.Config → Prop)
    (hcomplete : CompleteFor P R consistent)
    (start : P.Config) :
    ∃ (world : P.PublicWorld) (schedule : List R.Move),
      CompletedSchedule P R schedule start ∧
      consistent (run P R schedule start) ∧
      world = P.toPublicWorld (run P R schedule start) ∧
      ∀ (otherStart : P.Config),
        P.toPublicWorld otherStart = P.toPublicWorld start →
        ∀ other : List R.Move, CompletedSchedule P R other otherStart →
          P.toPublicWorld (run P R other otherStart) = world := by
  obtain ⟨schedule, hcompleted⟩ :=
    exists_completedSchedule P R hterm start
  refine ⟨P.toPublicWorld (run P R schedule start), schedule,
    hcompleted, (hcomplete _).1 hcompleted, rfl, ?_⟩
  intro otherStart hstart other hother
  have hreadback : P.readback start = P.readback otherStart :=
    ((P.toPublicWorld_eq_iff _ _).1 hstart).symm
  exact (completed_public_endpoints_equal_of_equiv P R hconf hcompat
    hreadback hcompleted hother).symm

end FiniteRepairSystem

/-! ## Explicit handoff to the existing OPH primitive repair stack -/

namespace OPHPrimitiveEndpoint

open FiniteRepairSystem

variable (C : OPH.OPHCarrier)

/-- Explicit finiteness receipt for every dependent patch-state type. -/
class FinitePatchStates where
  stateFintype : ∀ i : C.Patch, Fintype (C.State i)

variable [FinitePatchStates C]

/-- Expose the dependent finiteness receipt to Lean's instance search.  The
`reducible` annotation is needed because `Fintype.piFintype` must see through
the dependent state family when constructing a `Fintype (OPH.Records C)`. -/
@[reducible] noncomputable instance finitePatchStateFintype
    (i : C.Patch) : Fintype (C.State i) :=
  FinitePatchStates.stateFintype i

/-- The dependent product of finite patch-state types is finite. -/
@[reducible] noncomputable instance recordsFintype :
    Fintype (OPH.Records C) :=
  @Pi.instFintype C.Patch C.State C.patchDecEq C.patchFintype
    (fun i => FinitePatchStates.stateFintype i)

/-- The actual OPH raw-record presentation: configurations are
`OPH.Records C`, and public readback is exactly `OPH.obsMap C`.  An explicit
seed certifies that the dependent record type is inhabited. -/
@[reducible] noncomputable def presentation (seed : OPH.Records C) :
    PublicWorldPresentation (OPH.Obs C) where
  Config := OPH.Records C
  configFintype := recordsFintype C
  configNonempty := ⟨seed⟩
  readback := OPH.obsMap C

/-- The labelled repair system is the repository's constructed
`OPH.localRepair`; its accepted relation is therefore `OPH.acceptedStep`. -/
@[reducible] noncomputable def localSystem (seed : OPH.Records C) :
    (presentation C seed).FiniteRepairSystem where
  Move := OPH.Site C
  moveFintype := C.patchFintype
  repair := OPH.localRepair C

/-- The A4 labelled step relation is definitionally the existing OPH
accepted asynchronous step relation. -/
theorem local_step_eq_acceptedStep (seed : OPH.Records C) :
    step (presentation C seed) (localSystem C seed) = OPH.acceptedStep C :=
  rfl

/-- The existing `termination_holds` theorem discharges termination for the
actual constructed local-repair system. -/
theorem local_terminating (seed : OPH.Records C) :
    OPH.AbstractRewriting.Terminating
      (step (presentation C seed) (localSystem C seed)) := by
  change OPH.Termination C
  exact OPH.termination_holds C

/-- Both the output and the enabledness of every constructed OPH local move
depend only on `obsMap`.  Thus the labelled local system satisfies the exact
quotient-compatibility packet, not merely congruence of the composite
choice-selected `OPH.Repair`. -/
noncomputable def localQuotientCompatible (seed : OPH.Records C) :
    QuotientCompatible (presentation C seed) (localSystem C seed) where
  repair_respects_public := by
    intro i x y hxy
    exact OPH.obsMap_localRepair_congr C hxy i
  enabled_respects_public := by
    intro i x y hxy
    exact OPH.localRepair_fire_congr C hxy i

/-- The repository's canonical recursive `OPH.Repair` reaches a simultaneous
fixed point of every actual `OPH.localRepair` move. -/
theorem rawRepair_isFixed (seed x : OPH.Records C) :
    IsFixed (presentation C seed) (localSystem C seed) (OPH.Repair C x) := by
  intro i
  by_contra hfire
  exact OPH.Repair_normalForm C x
    (OPH.localRepair C i (OPH.Repair C x)) ⟨i, rfl, hfire⟩

omit [FinitePatchStates C] in
/-- The canonical recursive repair is idempotent on raw records because the
first application is a normal form. -/
theorem rawRepair_idempotent :
    ∀ x : OPH.Records C, OPH.Repair C (OPH.Repair C x) = OPH.Repair C x := by
  intro x
  apply OPH.Repair_of_normal
  intro hfire
  obtain ⟨i, hi⟩ := hfire
  exact OPH.Repair_normalForm C x
    (OPH.localRepair C i (OPH.Repair C x)) ⟨i, rfl, hi⟩

/-- `OPH.Repair` descends to the literal public-world quotient because the
existing `repair_respects_gauge` theorem makes it a congruence for
`Setoid.ker (obsMap C)`. -/
noncomputable def publicRepair (seed : OPH.Records C) :
    (presentation C seed).PublicWorld → (presentation C seed).PublicWorld :=
  Quotient.lift
    (fun x => (presentation C seed).toPublicWorld (OPH.Repair C x))
    (by
      intro x y hxy
      apply ((presentation C seed).toPublicWorld_eq_iff _ _).2
      exact OPH.repair_respects_gauge C x y hxy)

@[simp]
theorem publicRepair_mk (seed x : OPH.Records C) :
    publicRepair C seed ((presentation C seed).toPublicWorld x) =
      (presentation C seed).toPublicWorld (OPH.Repair C x) :=
  rfl

/-- The quotient repair operator is idempotent: its image is a public
fixed-point endpoint. -/
theorem publicRepair_idempotent (seed : OPH.Records C) :
    ∀ q, publicRepair C seed (publicRepair C seed q) = publicRepair C seed q := by
  intro q
  refine Quotient.inductionOn q ?_
  intro x
  change
    publicRepair C seed
        ((presentation C seed).toPublicWorld (OPH.Repair C x)) =
      (presentation C seed).toPublicWorld (OPH.Repair C x)
  rw [publicRepair_mk, rawRepair_idempotent C x]

/-- Gauge-equivalent raw representatives have exactly the same canonical
public endpoint.  This is the representative-independence theorem absent
from same-start confluence alone. -/
theorem canonicalEndpoint_rep_independent (seed : OPH.Records C)
    {x y : OPH.Records C} (hxy : OPH.gaugeEquiv C x y) :
    (presentation C seed).toPublicWorld (OPH.Repair C x) =
      (presentation C seed).toPublicWorld (OPH.Repair C y) := by
  apply ((presentation C seed).toPublicWorld_eq_iff _ _).2
  exact OPH.repair_respects_gauge C x y hxy

/-- Every canonical OPH endpoint belongs to the A4 public fixed-point object
for the actual labelled local-repair system. -/
theorem canonicalEndpoint_mem_fixedPointObject (seed x : OPH.Records C) :
    ∃ w : FixedPointObject (presentation C seed) (localSystem C seed),
      w.1 = (presentation C seed).toPublicWorld (OPH.Repair C x) := by
  exact ⟨⟨(presentation C seed).toPublicWorld (OPH.Repair C x),
    ⟨⟨OPH.Repair C x, rawRepair_isFixed C seed x⟩, rfl⟩⟩, rfl⟩

/-- For the actual constructed local moves, explicit confluence is enough to
upgrade same-start schedule independence to independence from both schedule
and gauge-equivalent raw representative.  Confluence is deliberately a
premise: it is false for some OPH carriers. -/
theorem local_completed_endpoints_equal_of_gauge
    (seed : OPH.Records C)
    (hconf : OPH.AbstractRewriting.Confluent (OPH.acceptedStep C))
    {x y : OPH.Records C} (hxy : OPH.gaugeEquiv C x y)
    {first second : List (OPH.Site C)}
    (hfirst : CompletedSchedule (presentation C seed) (localSystem C seed)
      first x)
    (hsecond : CompletedSchedule (presentation C seed) (localSystem C seed)
      second y) :
    (presentation C seed).toPublicWorld
        (run (presentation C seed) (localSystem C seed) first x) =
      (presentation C seed).toPublicWorld
        (run (presentation C seed) (localSystem C seed) second y) := by
  exact completed_public_endpoints_equal_of_equiv
    (presentation C seed) (localSystem C seed) hconf
    (localQuotientCompatible C seed) hxy hfirst hsecond

/-! ### The frustration-free `acceptedStepLR` completeness packet -/

/-- Adapt any explicit OPH local move `lr` to the same public presentation.
Its A4 step relation is definitionally `OPH.acceptedStepLR lr`. -/
@[reducible] noncomputable def lrSystem (seed : OPH.Records C)
    (lr : C.Patch → OPH.Records C → OPH.Records C) :
    (presentation C seed).FiniteRepairSystem where
  Move := C.Patch
  moveFintype := C.patchFintype
  repair := lr

theorem lr_step_eq_acceptedStepLR (seed : OPH.Records C)
    (lr : C.Patch → OPH.Records C → OPH.Records C) :
    step (presentation C seed) (lrSystem C seed lr) = OPH.acceptedStepLR lr :=
  rfl

variable (lr : C.Patch → OPH.Records C → OPH.Records C)
variable
  (H1 : ∀ (i : C.Patch) (x : OPH.Records C) (j : C.Patch),
    j ≠ i → (lr i x) j = x j)
  (H2 : ∀ (i : C.Patch) (x : OPH.Records C),
    lr i x ≠ x ↔ ∃ e : C.Edge,
      (C.src e = i ∨ C.tgt e = i) ∧ ¬ OPH.edgeConsistentAt e x)
  (H3 : ∀ (i : C.Patch) (x : OPH.Records C),
    lr i x ≠ x → ∀ e : C.Edge,
      (C.src e = i ∨ C.tgt e = i) → OPH.edgeConsistentAt e (lr i x))

include H1 H2 H3

/-- The existing H1--H3 descent theorem supplies termination for the adapted
`acceptedStepLR` system. -/
theorem lr_terminating (seed : OPH.Records C) :
    OPH.AbstractRewriting.Terminating
      (step (presentation C seed) (lrSystem C seed lr)) := by
  change WellFounded (fun y x : OPH.Records C => OPH.acceptedStepLR lr x y)
  exact OPH.termination lr H1 H2 H3

/-- The existing H1--H3 completeness theorem supplies the explicit semantic
bridge: simultaneous repair fixed points are exactly `OPH.Consistent`
records. -/
theorem lr_completeFor (seed : OPH.Records C) :
    CompleteFor (presentation C seed) (lrSystem C seed lr)
      (OPH.Consistent C) := by
  intro x
  change (∀ i : C.Patch, lr i x = x) ↔ OPH.Consistent C x
  exact (OPH.normalForm_iff_all_quiescent lr H1 H2 H3 x).symm.trans
    (OPH.completeness lr H1 H2 H3 x)

/-- An arbitrary H1--H3 local repair needs an additional explicit gauge
packet before it can act on quotient classes.  H1--H3 alone prove descent of
the mismatch measure and semantic completeness, not quotient congruence. -/
noncomputable def lrQuotientCompatible (seed : OPH.Records C)
    (Hgauge : ∀ (i : C.Patch) (x y : OPH.Records C),
      OPH.gaugeEquiv C x y →
        OPH.gaugeEquiv C (lr i x) (lr i y))
    (Henabled : ∀ (i : C.Patch) (x y : OPH.Records C),
      OPH.gaugeEquiv C x y →
        (lr i x ≠ x ↔ lr i y ≠ y)) :
    QuotientCompatible (presentation C seed) (lrSystem C seed lr) where
  repair_respects_public := by
    intro i x y hxy
    exact Hgauge i x y hxy
  enabled_respects_public := by
    intro i x y hxy
    exact Henabled i x y hxy

/-- The complete OPH-specific conditional endpoint theorem.  H1--H3 supply
termination and `normal form ↔ OPH.Consistent`; `hconf` supplies schedule
independence; `Hgauge` and `Henabled` supply representative independence;
and `CompletedSchedule` is the explicitly named finite completion postcondition.
No physical time or infinite scheduler fairness is asserted. -/
theorem lr_public_endpoint_exists_unique_on_gauge_class
    (Hgauge : ∀ (i : C.Patch) (x y : OPH.Records C),
      OPH.gaugeEquiv C x y →
        OPH.gaugeEquiv C (lr i x) (lr i y))
    (Henabled : ∀ (i : C.Patch) (x y : OPH.Records C),
      OPH.gaugeEquiv C x y →
        (lr i x ≠ x ↔ lr i y ≠ y))
    (hconf : OPH.AbstractRewriting.Confluent (OPH.acceptedStepLR lr))
    (seed start : OPH.Records C) :
    ∃ (world : (presentation C seed).PublicWorld)
        (schedule : List C.Patch),
      CompletedSchedule (presentation C seed) (lrSystem C seed lr)
        schedule start ∧
      OPH.Consistent C
        (run (presentation C seed) (lrSystem C seed lr) schedule start) ∧
      world = (presentation C seed).toPublicWorld
        (run (presentation C seed) (lrSystem C seed lr) schedule start) ∧
      ∀ (otherStart : OPH.Records C), OPH.gaugeEquiv C otherStart start →
        ∀ other : List C.Patch,
          CompletedSchedule (presentation C seed) (lrSystem C seed lr)
            other otherStart →
          (presentation C seed).toPublicWorld
              (run (presentation C seed) (lrSystem C seed lr)
                other otherStart) = world := by
  obtain ⟨world, schedule, hdone, hconsistent, hworld, hunique⟩ :=
    public_endpoint_exists_unique_on_public_class
      (presentation C seed) (lrSystem C seed lr)
      (lr_terminating C lr H1 H2 H3 seed) hconf
      (lrQuotientCompatible C lr seed Hgauge Henabled)
      (OPH.Consistent C) (lr_completeFor C lr H1 H2 H3 seed) start
  refine ⟨world, schedule, hdone, hconsistent, hworld, ?_⟩
  intro otherStart hstart other hother
  apply hunique otherStart
  · exact ((presentation C seed).toPublicWorld_eq_iff _ _).2 hstart
  · exact hother

noncomputable instance demoFinitePatchStates :
    FinitePatchStates OPH.demoCarrier where
  stateFintype := fun _ => by
    change Fintype Bool
    infer_instance

omit [FinitePatchStates C] H1 H2 H3 in
/-- The existing two-patch copy repair is a nonempty concrete instance of
the H1--H3 termination-plus-completeness packet. -/
theorem demoLR_completeFor :
    CompleteFor
      (presentation OPH.demoCarrier (fun b => b))
      (lrSystem OPH.demoCarrier (fun b => b) OPH.demoLR)
      (OPH.Consistent OPH.demoCarrier) :=
  lr_completeFor OPH.demoCarrier OPH.demoLR OPH.demoLR_H1
    OPH.demoLR_H2 OPH.demoLR_H3 (fun b => b)

omit [FinitePatchStates C] H1 H2 H3 in
theorem demoLR_terminating :
    OPH.AbstractRewriting.Terminating
      (step
        (presentation OPH.demoCarrier (fun b => b))
        (lrSystem OPH.demoCarrier (fun b => b) OPH.demoLR)) :=
  lr_terminating OPH.demoCarrier OPH.demoLR OPH.demoLR_H1
    OPH.demoLR_H2 OPH.demoLR_H3 (fun b => b)

end OPHPrimitiveEndpoint

/-! ## Countermodels -/

namespace EndpointCountermodels

open FiniteRepairSystem

/-! ### Removing confluence: a terminating fork with two completed endpoints -/

inductive ForkConfig
  | start | left | right
  deriving DecidableEq, Fintype

inductive ForkMove
  | chooseLeft | chooseRight
  deriving DecidableEq, Fintype

/-- The two terminal fork branches have different public signatures. -/
def forkReadback : ForkConfig → Bool
  | .start | .left => false
  | .right => true

def forkPresentation : PublicWorldPresentation Bool where
  Config := ForkConfig
  configFintype := inferInstance
  configNonempty := ⟨.start⟩
  readback := forkReadback

def forkRepair : ForkMove → ForkConfig → ForkConfig
  | .chooseLeft, .start => .left
  | .chooseRight, .start => .right
  | _, x => x

def forkSystem : forkPresentation.FiniteRepairSystem where
  Move := ForkMove
  moveFintype := inferInstance
  repair := forkRepair

def forkMeasure : ForkConfig → ℕ
  | .start => 1
  | .left | .right => 0

theorem fork_measure_desc {x y : ForkConfig}
    (h : step forkPresentation forkSystem x y) :
    forkMeasure y < forkMeasure x := by
  obtain ⟨m, rfl, hne⟩ := h
  cases m <;> cases x <;> simp_all [forkSystem, forkRepair, forkMeasure]

/-- The fork is terminating. -/
theorem fork_terminating :
    OPH.AbstractRewriting.Terminating (step forkPresentation forkSystem) :=
  terminating_of_measure forkPresentation forkSystem forkMeasure fork_measure_desc

/-- Choosing the left branch completes the repair. -/
theorem fork_left_completed :
    CompletedSchedule forkPresentation forkSystem
      [ForkMove.chooseLeft] ForkConfig.start := by
  intro m
  cases m <;> rfl

/-- Choosing the right branch completes the repair. -/
theorem fork_right_completed :
    CompletedSchedule forkPresentation forkSystem
      [ForkMove.chooseRight] ForkConfig.start := by
  intro m
  cases m <;> rfl

/-- Without confluence, two completed schedules from the same start can determine
different public-world classes. -/
theorem fork_public_endpoints_differ :
    forkPresentation.toPublicWorld
        (run forkPresentation forkSystem
          [ForkMove.chooseLeft] ForkConfig.start) ≠
      forkPresentation.toPublicWorld
        (run forkPresentation forkSystem
          [ForkMove.chooseRight] ForkConfig.start) := by
  intro h
  have hs := (forkPresentation.toPublicWorld_eq_iff _ _).1 h
  exact Bool.false_ne_true hs

/-- The terminating fork is genuinely nonconfluent. -/
theorem fork_not_confluent :
    ¬ OPH.AbstractRewriting.Confluent (step forkPresentation forkSystem) := by
  intro hconf
  have heq := completed_public_endpoints_equal forkPresentation forkSystem
    hconf fork_left_completed fork_right_completed
  exact fork_public_endpoints_differ heq

/-! ### Removing completion: a terminating confluent run may stop early -/

def settlePresentation : PublicWorldPresentation Bool where
  Config := Bool
  configFintype := inferInstance
  configNonempty := ⟨false⟩
  readback := id

def settleRepair (_ : Unit) (_ : Bool) : Bool := true

def settleSystem : settlePresentation.FiniteRepairSystem where
  Move := Unit
  moveFintype := inferInstance
  repair := settleRepair

def settleMeasure (x : Bool) : ℕ := if x then 0 else 1

theorem settle_measure_desc {x y : Bool}
    (h : step settlePresentation settleSystem x y) :
    settleMeasure y < settleMeasure x := by
  obtain ⟨m, rfl, hne⟩ := h
  cases m
  cases x <;> simp_all [settleSystem, settleRepair, settleMeasure]

theorem settle_terminating :
    OPH.AbstractRewriting.Terminating (step settlePresentation settleSystem) :=
  terminating_of_measure settlePresentation settleSystem
    settleMeasure settle_measure_desc

theorem settle_locallyConfluent :
    OPH.AbstractRewriting.LocallyConfluent
      (step settlePresentation settleSystem) := by
  rintro x y z ⟨m, rfl, _⟩ ⟨n, rfl, _⟩
  cases m
  cases n
  exact ⟨settleRepair () x, ReflTransGen.refl, ReflTransGen.refl⟩

/-- The settling system is confluent as well as terminating. -/
theorem settle_confluent :
    OPH.AbstractRewriting.Confluent (step settlePresentation settleSystem) :=
  OPH.AbstractRewriting.newman_lemma (step settlePresentation settleSystem)
    settle_terminating settle_locallyConfluent

/-- The empty schedule is not completed at the unsettled state. -/
theorem settle_empty_not_completed :
    ¬ CompletedSchedule settlePresentation settleSystem [] false := by
  intro hcompleted
  have h := hcompleted ()
  exact (by decide : true ≠ false) h

/-- One settling move is completed. -/
theorem settle_one_move_completed :
    CompletedSchedule settlePresentation settleSystem [()] false := by
  intro m
  cases m
  rfl

/-- Even with termination and confluence, dropping completion lets a truncated
schedule report a nonfixed public class different from the actual endpoint. -/
theorem settle_incomplete_public_value_differs :
    settlePresentation.toPublicWorld
        (run settlePresentation settleSystem [] false) ≠
      settlePresentation.toPublicWorld
        (run settlePresentation settleSystem [()] false) := by
  intro h
  have hs := (settlePresentation.toPublicWorld_eq_iff _ _).1 h
  exact Bool.false_ne_true hs

/-! ### Removing quotient congruence: representative dependence survives
termination, confluence, and completion -/

/-- The first bit is public and the second bit is hidden. -/
def representativePresentation : PublicWorldPresentation Bool where
  Config := Bool × Bool
  configFintype := inferInstance
  configNonempty := ⟨(false, false)⟩
  readback := Prod.fst

/-- Copy the hidden bit into both positions.  This deterministic operation is
idempotent and terminating, but it lets hidden representative data determine
the eventual public bit. -/
def representativeRepair (_ : Unit) (x : Bool × Bool) : Bool × Bool :=
  (x.2, x.2)

def representativeSystem : representativePresentation.FiniteRepairSystem where
  Move := Unit
  moveFintype := inferInstance
  repair := representativeRepair

def representativeMeasure (x : Bool × Bool) : ℕ :=
  if x.1 = x.2 then 0 else 1

theorem representative_measure_desc {x y : Bool × Bool}
    (h : step representativePresentation representativeSystem x y) :
    representativeMeasure y < representativeMeasure x := by
  obtain ⟨m, rfl, hne⟩ := h
  cases m
  rcases x with ⟨a, b⟩
  cases a <;> cases b <;>
    simp_all [representativeSystem, representativeRepair, representativeMeasure]

theorem representative_terminating :
    OPH.AbstractRewriting.Terminating
      (step representativePresentation representativeSystem) :=
  terminating_of_measure representativePresentation representativeSystem
    representativeMeasure representative_measure_desc

theorem representative_locallyConfluent :
    OPH.AbstractRewriting.LocallyConfluent
      (step representativePresentation representativeSystem) := by
  rintro x y z ⟨m, rfl, _⟩ ⟨n, rfl, _⟩
  cases m
  cases n
  exact ⟨representativeRepair () x, ReflTransGen.refl, ReflTransGen.refl⟩

/-- Even the raw rewrite relation is confluent. -/
theorem representative_confluent :
    OPH.AbstractRewriting.Confluent
      (step representativePresentation representativeSystem) :=
  OPH.AbstractRewriting.newman_lemma
    (step representativePresentation representativeSystem)
    representative_terminating representative_locallyConfluent

theorem representative_left_completed :
    CompletedSchedule representativePresentation representativeSystem
      [()] (false, false) := by
  intro m
  cases m
  rfl

theorem representative_right_completed :
    CompletedSchedule representativePresentation representativeSystem
      [()] (false, true) := by
  intro m
  cases m
  rfl

/-- The two starts are the same public quotient class. -/
theorem representative_initial_public_eq :
    representativePresentation.toPublicWorld (false, false) =
      representativePresentation.toPublicWorld (false, true) :=
  (representativePresentation.toPublicWorld_eq_iff _ _).2 rfl

/-- Their completed endpoints expose opposite public bits. -/
theorem representative_completed_endpoints_differ :
    representativePresentation.toPublicWorld
        (run representativePresentation representativeSystem [()] (false, false)) ≠
      representativePresentation.toPublicWorld
        (run representativePresentation representativeSystem [()] (false, true)) := by
  intro h
  have hs := (representativePresentation.toPublicWorld_eq_iff _ _).1 h
  exact Bool.false_ne_true hs

/-- Consequently the repair cannot descend to the quotient.  This is the
precise missing hypothesis: same-start confluence does not imply that a repair
respects public/gauge-equivalent representatives. -/
theorem representative_no_descended_repair :
    ¬ ∃ f : representativePresentation.PublicWorld →
        representativePresentation.PublicWorld,
      ∀ x : Bool × Bool,
        f (representativePresentation.toPublicWorld x) =
          representativePresentation.toPublicWorld
            (representativeRepair () x) := by
  rintro ⟨f, hf⟩
  have h := congrArg f representative_initial_public_eq
  rw [hf, hf] at h
  exact representative_completed_endpoints_differ h

end EndpointCountermodels

end PublicWorldPresentation

/-! ## Typed OPH primitive adaptor at an A3 consensus-tower regulator -/

namespace ConsensusTower

open PublicWorldPresentation
open PublicWorldPresentation.FiniteRepairSystem
open PublicWorldPresentation.OPHPrimitiveEndpoint

variable {ι : Type u} [Preorder ι] (T : ConsensusTower ι) (r : ι)
variable (C : OPH.OPHCarrier) [FinitePatchStates C]

/-- The extra realization data needed to attach an OPH primitive carrier to
one regulator of an A3 `ConsensusTower`.  A3 supplies only observer and record
types.  The adaptor therefore names a raw seed and an injective encoding of
the complete OPH overlap observation into the tower's observer-indexed public
signature.  Injectivity is exactly what makes the induced public quotient
equal to OPH gauge equivalence rather than a coarser quotient. -/
structure OPHPrimitiveReadbackAdapter where
  seed : OPH.Records C
  encode : OPH.Obs C → T.PublicSignature r
  encode_injective : Function.Injective encode

variable (A : OPHPrimitiveReadbackAdapter T r C)

/-- The literal `ConsensusTower.atRegulator` presentation selected by the
adaptor. -/
@[reducible] noncomputable def primitivePresentation :
    PublicWorldPresentation (T.PublicSignature r) := by
  letI : Fintype (OPH.Records C) := recordsFintype C
  letI : Nonempty (OPH.Records C) := ⟨A.seed⟩
  exact T.atRegulator r (OPH.Records C)
    (fun x => A.encode (OPH.obsMap C x))

/-- Equality in the tower-regulator public quotient is exactly the existing
OPH `gaugeEquiv`, because the adaptor's public encoding is injective. -/
theorem primitive_toPublicWorld_eq_iff_gauge (x y : OPH.Records C) :
    (primitivePresentation T r C A).toPublicWorld x =
        (primitivePresentation T r C A).toPublicWorld y ↔
      OPH.gaugeEquiv C x y := by
  rw [(primitivePresentation T r C A).toPublicWorld_eq_iff]
  change A.encode (OPH.obsMap C x) = A.encode (OPH.obsMap C y) ↔
    OPH.gaugeEquiv C x y
  constructor
  · intro h
    change OPH.obsMap C x = OPH.obsMap C y
    exact A.encode_injective h
  · intro h
    exact congrArg A.encode h

/-- The actual constructed OPH local repair family, attached to the
selected A3 regulator presentation. -/
@[reducible] noncomputable def primitiveLocalSystem :
    (primitivePresentation T r C A).FiniteRepairSystem where
  Move := OPH.Site C
  moveFintype := C.patchFintype
  repair := OPH.localRepair C

/-- The attached labelled relation is definitionally the repository's
`OPH.acceptedStep`. -/
theorem primitive_step_eq_acceptedStep :
    step (primitivePresentation T r C A) (primitiveLocalSystem T r C A) =
      OPH.acceptedStep C :=
  rfl

/-- Existing OPH single-move congruence and firing congruence make the
tower-attached local system quotient-compatible. -/
noncomputable def primitiveLocalQuotientCompatible :
    QuotientCompatible (primitivePresentation T r C A)
      (primitiveLocalSystem T r C A) where
  repair_respects_public := by
    intro i x y hxy
    exact congrArg A.encode
      (OPH.obsMap_localRepair_congr C (A.encode_injective hxy) i)
  enabled_respects_public := by
    intro i x y hxy
    exact OPH.localRepair_fire_congr C (A.encode_injective hxy) i

/-- The canonical existing `OPH.Repair` descends to an endomorphism of the
A3 regulator's public quotient. -/
noncomputable def primitivePublicRepair :
    (primitivePresentation T r C A).PublicWorld →
      (primitivePresentation T r C A).PublicWorld :=
  Quotient.lift
    (fun x => (primitivePresentation T r C A).toPublicWorld (OPH.Repair C x))
    (by
      intro x y hxy
      apply ((primitivePresentation T r C A).toPublicWorld_eq_iff _ _).2
      exact congrArg A.encode
        (OPH.repair_respects_gauge C x y (A.encode_injective hxy)))

@[simp]
theorem primitivePublicRepair_mk (x : OPH.Records C) :
    primitivePublicRepair T r C A
        ((primitivePresentation T r C A).toPublicWorld x) =
      (primitivePresentation T r C A).toPublicWorld (OPH.Repair C x) :=
  rfl

/-- The descended tower-regulator endpoint operator is idempotent. -/
theorem primitivePublicRepair_idempotent :
    ∀ q, primitivePublicRepair T r C A
        (primitivePublicRepair T r C A q) =
      primitivePublicRepair T r C A q := by
  intro q
  refine Quotient.inductionOn q ?_
  intro x
  change primitivePublicRepair T r C A
      ((primitivePresentation T r C A).toPublicWorld (OPH.Repair C x)) =
    (primitivePresentation T r C A).toPublicWorld (OPH.Repair C x)
  rw [primitivePublicRepair_mk,
    PublicWorldPresentation.OPHPrimitiveEndpoint.rawRepair_idempotent C x]

/-- The tower-attached public fixed-point object is inhabited by the
canonical OPH repair endpoint. -/
theorem primitiveCanonicalEndpoint_mem_fixedPointObject (x : OPH.Records C) :
    ∃ w : FixedPointObject (primitivePresentation T r C A)
        (primitiveLocalSystem T r C A),
      w.1 = (primitivePresentation T r C A).toPublicWorld
        (OPH.Repair C x) := by
  have hfixed : IsFixed (primitivePresentation T r C A)
      (primitiveLocalSystem T r C A) (OPH.Repair C x) := by
    intro i
    by_contra hfire
    exact OPH.Repair_normalForm C x
      (OPH.localRepair C i (OPH.Repair C x)) ⟨i, rfl, hfire⟩
  exact ⟨⟨(primitivePresentation T r C A).toPublicWorld (OPH.Repair C x),
    ⟨⟨OPH.Repair C x, hfixed⟩, rfl⟩⟩, rfl⟩

/-- Termination of the existing constructed local repair supplies
inhabitedness independently of any confluence claim. -/
theorem primitiveFixedPointObject_nonempty :
    Nonempty (FixedPointObject (primitivePresentation T r C A)
      (primitiveLocalSystem T r C A)) := by
  apply fixedPointObject_nonempty
  change OPH.Termination C
  exact OPH.termination_holds C

/-- Attach an arbitrary explicit H1--H3 local repair family to the same A3
regulator presentation. -/
@[reducible] noncomputable def primitiveLRSystem
    (lr : C.Patch → OPH.Records C → OPH.Records C) :
    (primitivePresentation T r C A).FiniteRepairSystem where
  Move := C.Patch
  moveFintype := C.patchFintype
  repair := lr

/-- Explicit quotient packet for a caller-supplied local repair family at the
tower regulator. -/
noncomputable def primitiveLRQuotientCompatible
    (lr : C.Patch → OPH.Records C → OPH.Records C)
    (Hgauge : ∀ (i : C.Patch) (x y : OPH.Records C),
      OPH.gaugeEquiv C x y → OPH.gaugeEquiv C (lr i x) (lr i y))
    (Henabled : ∀ (i : C.Patch) (x y : OPH.Records C),
      OPH.gaugeEquiv C x y → (lr i x ≠ x ↔ lr i y ≠ y)) :
    QuotientCompatible (primitivePresentation T r C A)
      (primitiveLRSystem T r C A lr) where
  repair_respects_public := by
    intro i x y hxy
    exact congrArg A.encode (Hgauge i x y (A.encode_injective hxy))
  enabled_respects_public := by
    intro i x y hxy
    exact Henabled i x y (A.encode_injective hxy)

/-- Complete OPH/A3 endpoint packet.  Every hypothesis not supplied by H1--H3
is visible: confluence, gauge congruence, enabledness congruence, and finite
terminal completion.  The conclusion is independent of both schedule and raw
representative throughout one gauge/public class. -/
theorem primitiveLR_endpoint_exists_unique_on_gauge_class
    (lr : C.Patch → OPH.Records C → OPH.Records C)
    (H1 : ∀ (i : C.Patch) (x : OPH.Records C) (j : C.Patch),
      j ≠ i → (lr i x) j = x j)
    (H2 : ∀ (i : C.Patch) (x : OPH.Records C),
      lr i x ≠ x ↔ ∃ e : C.Edge,
        (C.src e = i ∨ C.tgt e = i) ∧ ¬ OPH.edgeConsistentAt e x)
    (H3 : ∀ (i : C.Patch) (x : OPH.Records C),
      lr i x ≠ x → ∀ e : C.Edge,
        (C.src e = i ∨ C.tgt e = i) →
          OPH.edgeConsistentAt e (lr i x))
    (Hgauge : ∀ (i : C.Patch) (x y : OPH.Records C),
      OPH.gaugeEquiv C x y → OPH.gaugeEquiv C (lr i x) (lr i y))
    (Henabled : ∀ (i : C.Patch) (x y : OPH.Records C),
      OPH.gaugeEquiv C x y → (lr i x ≠ x ↔ lr i y ≠ y))
    (hconf : OPH.AbstractRewriting.Confluent (OPH.acceptedStepLR lr))
    (start : OPH.Records C) :
    ∃ (world : (primitivePresentation T r C A).PublicWorld)
        (schedule : List C.Patch),
      CompletedSchedule (primitivePresentation T r C A)
          (primitiveLRSystem T r C A lr) schedule start ∧
      OPH.Consistent C
        (run (primitivePresentation T r C A)
          (primitiveLRSystem T r C A lr) schedule start) ∧
      world = (primitivePresentation T r C A).toPublicWorld
        (run (primitivePresentation T r C A)
          (primitiveLRSystem T r C A lr) schedule start) ∧
      ∀ (otherStart : OPH.Records C), OPH.gaugeEquiv C otherStart start →
        ∀ other : List C.Patch,
          CompletedSchedule (primitivePresentation T r C A)
              (primitiveLRSystem T r C A lr) other otherStart →
          (primitivePresentation T r C A).toPublicWorld
              (run (primitivePresentation T r C A)
                (primitiveLRSystem T r C A lr) other otherStart) = world := by
  have hterm : OPH.AbstractRewriting.Terminating
      (step (primitivePresentation T r C A)
        (primitiveLRSystem T r C A lr)) := by
    change WellFounded (fun y x : OPH.Records C => OPH.acceptedStepLR lr x y)
    exact OPH.termination lr H1 H2 H3
  have hcomplete : CompleteFor (primitivePresentation T r C A)
      (primitiveLRSystem T r C A lr) (OPH.Consistent C) := by
    intro x
    change (∀ i : C.Patch, lr i x = x) ↔ OPH.Consistent C x
    exact (OPH.normalForm_iff_all_quiescent lr H1 H2 H3 x).symm.trans
      (OPH.completeness lr H1 H2 H3 x)
  obtain ⟨world, schedule, hdone, hconsistent, hworld, hunique⟩ :=
    public_endpoint_exists_unique_on_public_class
      (primitivePresentation T r C A) (primitiveLRSystem T r C A lr)
      hterm hconf
      (primitiveLRQuotientCompatible T r C A lr Hgauge Henabled)
      (OPH.Consistent C) hcomplete start
  refine ⟨world, schedule, hdone, hconsistent, hworld, ?_⟩
  intro otherStart hstart other hother
  apply hunique otherStart
  · exact (primitive_toPublicWorld_eq_iff_gauge T r C A _ _).2 hstart
  · exact hother

end ConsensusTower

-- Axiom audit: no admission or project-specific axiom is permitted here.
#print axioms OPH.Tower.PublicWorldPresentation.FiniteRepairSystem.exists_completedSchedule
#print axioms OPH.Tower.PublicWorldPresentation.FiniteRepairSystem.fixedPointObject_nonempty
#print axioms OPH.Tower.PublicWorldPresentation.FiniteRepairSystem.completed_runs_equal
#print axioms OPH.Tower.PublicWorldPresentation.FiniteRepairSystem.public_endpoint_exists_unique
#print axioms OPH.Tower.PublicWorldPresentation.FiniteRepairSystem.public_endpoint_exists_unique_on_public_class
#print axioms OPH.Tower.PublicWorldPresentation.OPHPrimitiveEndpoint.publicRepair_idempotent
#print axioms OPH.Tower.PublicWorldPresentation.OPHPrimitiveEndpoint.lr_public_endpoint_exists_unique_on_gauge_class
#print axioms OPH.Tower.PublicWorldPresentation.EndpointCountermodels.fork_not_confluent
#print axioms OPH.Tower.PublicWorldPresentation.EndpointCountermodels.settle_confluent
#print axioms OPH.Tower.PublicWorldPresentation.EndpointCountermodels.settle_empty_not_completed
#print axioms OPH.Tower.PublicWorldPresentation.EndpointCountermodels.representative_confluent
#print axioms OPH.Tower.PublicWorldPresentation.EndpointCountermodels.representative_no_descended_repair
#print axioms OPH.Tower.ConsensusTower.primitivePublicRepair_idempotent
#print axioms OPH.Tower.ConsensusTower.primitiveLR_endpoint_exists_unique_on_gauge_class

end OPH.Tower
