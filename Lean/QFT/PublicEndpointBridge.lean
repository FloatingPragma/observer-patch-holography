import QFT.JointInstance
import Tower.FixedPointEndpoint

/-!
# The A4 public-endpoint bridge at the joint instance

This module implements the A4 realization/naturality bridge of
completion-plan issue `#693`.  It attacks bounded fragments of inherited
obligations 1 and 3: a public readback adaptor tying the A4 kernel
quotient to the bare tower's physical quotient, and raw-configuration
refinement naturality across a bare regulator pair.

`PublicReadbackAdaptor` declares, at one regulator of a bare consensus
tower, an A4 public readback on the tower's raw finite presentation
states together with both kernel inclusions: readback equality implies
physical-quotient equality and conversely.  The A4 kernel quotient of
the induced presentation is therefore the bare tower's physical quotient
through the derived equivalence `publicWorldEquivQuot`, whose two
composites are the identity (`toQuot_ofQuot`, `ofQuot_toQuot`) and whose
compatibility triangles `toQuot_mk` and `ofQuot_quotient` commute with
the two quotient maps.  The kernel identity is recorded literally as
`publicSetoid_eq_quotient_kernel`.

`PublicEndpointBridge` attaches an A4 labelled repair system to the
adaptor: the quotient-compatibility packet, and a raw normalizer that is
fixed under every labelled move, reachable by accepted steps, and
projects to the bare tower's declared normal form.  The
endpoint-compatibility theorem `endpoint_matches` shows that under
declared confluence every completed A4 schedule from every raw start
lands, through the adaptor, on the bare normal form of the start's
quotient class; `endpoint_matches_descended` restates the agreement
inside the A4 public world against the transported bare normal form
`descendedNormalForm`, and `endpoints_agree_of_quotient_eq` upgrades it
to independence from schedule and raw representative at once.
Confluence stays a premise, exactly as in the A4 packet, because it
fails for some carriers.

`RawRefinementSystem` states the obligation-3 raw-configuration
refinement maps covering the bare coarse maps.
`ReadbackRefinementNaturality` states readback naturality across one
regulator pair as the structure's law; its descended `publicCoarse`
commutes with the bare coarse maps (`publicCoarse_toQuot`) and with the
descended signatures (`publicCoarse_signature`).
`RepairRefinementNaturality` states repair and normalizer naturality,
with the schedule-level consequence `run_natural`.  For the identity
case of the committed demo stack all three laws are proved outright and
the squares are evaluated at the regulator pair `0 ≤ 1` of the joint
witness instance's bare tower.

Two inhabitant families keep the bridge away from vacuity.  The identity
adaptor and copy-repair bridge on the committed `demoTower` glue the
interfaces to `jointWitnessInstance`:
`jointWitnessInstance_endpoint_reads_selected_slice` sends a completed
A4 endpoint through the adaptor into the joint instance's source-link
readout and reproduces the selected world slice.  The hidden-bit tower
supplies the pair where the equivalence is genuine: the readback is
nonconstant, distinct raw states with equal public part share one
quotient class, and the copy repair completes to the bare normal form
through the hidden bit.

## Claim boundary

The bridge is bounded gluing of two committed packets: the A4 kernel
quotient/endpoint interfaces and the E2 joint finite instance.  No
source produces the adaptor, the repair system, or the refinement maps;
no physical world, cross-regulator continuum endpoint, or
schedule-independent consistency claim beyond the stated conditional
theorems is made; and the C2, CP/CPTP, and full obligation-2 packets
remain open E2 scope.
-/

namespace OPH.QFT

open OPH.Tower
open OPH.Tower.PublicWorldPresentation
open OPH.Tower.PublicWorldPresentation.FiniteRepairSystem
open OPH.EinsteinBranch (BareConsensusTower demoTower demoMismatch)
open Relation

universe u u' v

/-! ## Obligation-1 fragment: the public readback adaptor -/

/-- An A4 public readback adaptor at one regulator of a bare consensus
tower.  The readback is an A4 public signature on the tower's raw finite
presentation states; the two kernel laws say that its kernel is exactly
the kernel of the tower's physical quotient map.  The induced A4 kernel
quotient is then the bare physical quotient through a derived
equivalence rather than by declaration. -/
structure PublicReadbackAdaptor {ι : Type u} [Preorder ι]
    (B : BareConsensusTower ι) (r : ι) where
  /-- The public signature carrier of the A4 readback. -/
  Signature : Type v
  /-- The raw presentation carrier is inhabited. -/
  stateNonempty : Nonempty (B.State r)
  /-- The A4 public readback on raw presentation states. -/
  readback : B.State r → Signature
  /-- Kernel inclusion, forward: readback equality refines the physical
  quotient. -/
  quotient_eq_of_readback : ∀ {x y : B.State r},
    readback x = readback y → B.quotient r x = B.quotient r y
  /-- Kernel inclusion, backward: physical-quotient equality refines
  readback equality. -/
  readback_eq_of_quotient : ∀ {x y : B.State r},
    B.quotient r x = B.quotient r y → readback x = readback y

namespace PublicReadbackAdaptor

variable {ι : Type u} [Preorder ι] {B : BareConsensusTower ι} {r : ι}
variable (A : PublicReadbackAdaptor B r)

/-- The literal A4 presentation selected by the adaptor: raw
configurations are the tower's finite presentation states at the
regulator and the public readback is the declared one. -/
@[reducible] def presentation : PublicWorldPresentation A.Signature where
  Config := B.State r
  configFintype := B.stateFinite r
  configNonempty := A.stateNonempty
  readback := A.readback

/-- The declared A4 public setoid is literally the kernel of the bare
tower's physical quotient map. -/
theorem publicSetoid_eq_quotient_kernel :
    A.presentation.publicSetoid = Setoid.ker (B.quotient r) :=
  Setoid.ext fun _ _ =>
    ⟨fun h => A.quotient_eq_of_readback h,
      fun h => A.readback_eq_of_quotient h⟩

/-- Two raw states share an A4 public class exactly when they share a
bare physical quotient class. -/
theorem toPublicWorld_eq_iff_quotient_eq (x y : B.State r) :
    A.presentation.toPublicWorld x = A.presentation.toPublicWorld y ↔
      B.quotient r x = B.quotient r y := by
  rw [A.presentation.toPublicWorld_eq_iff]
  exact ⟨fun h => A.quotient_eq_of_readback h,
    fun h => A.readback_eq_of_quotient h⟩

/-- The descended map from the A4 kernel quotient onto the bare physical
quotient. -/
def toQuot : A.presentation.PublicWorld → B.Quot r :=
  Quotient.lift (B.quotient r) (fun _ _ h => A.quotient_eq_of_readback h)

/-- Compatibility triangle, raw side: the descended map commutes with
the two quotient maps. -/
@[simp] theorem toQuot_mk (x : B.State r) :
    A.toQuot (A.presentation.toPublicWorld x) = B.quotient r x :=
  rfl

/-- The descended map is injective by the backward kernel inclusion. -/
theorem toQuot_injective : Function.Injective A.toQuot := by
  intro w₁ w₂ h
  refine Quotient.inductionOn₂ w₁ w₂ ?_ h
  intro x y hxy
  exact (A.toPublicWorld_eq_iff_quotient_eq x y).2 hxy

/-- The descended map is surjective because the bare physical quotient
map is surjective. -/
theorem toQuot_surjective : Function.Surjective A.toQuot := by
  intro q
  obtain ⟨x, hx⟩ := B.quotient_surjective r q
  exact ⟨A.presentation.toPublicWorld x, hx⟩

/-- The descended map is a bijection. -/
theorem toQuot_bijective : Function.Bijective A.toQuot :=
  ⟨A.toQuot_injective, A.toQuot_surjective⟩

/-- The A4 kernel quotient IS the bare tower's physical quotient: the
derived equivalence whose forward map is the descended readback
projection. -/
noncomputable def publicWorldEquivQuot :
    A.presentation.PublicWorld ≃ B.Quot r :=
  Equiv.ofBijective A.toQuot A.toQuot_bijective

/-- The inverse leg of the equivalence. -/
noncomputable def ofQuot : B.Quot r → A.presentation.PublicWorld :=
  A.publicWorldEquivQuot.symm

/-- The two declared legs are mutually inverse, first composite. -/
@[simp] theorem toQuot_ofQuot (q : B.Quot r) :
    A.toQuot (A.ofQuot q) = q :=
  A.publicWorldEquivQuot.apply_symm_apply q

/-- The two declared legs are mutually inverse, second composite. -/
@[simp] theorem ofQuot_toQuot (w : A.presentation.PublicWorld) :
    A.ofQuot (A.toQuot w) = w :=
  A.publicWorldEquivQuot.symm_apply_apply w

/-- Compatibility triangle, quotient side: the inverse leg carries a
bare quotient class to the A4 class of any of its raw representatives. -/
theorem ofQuot_quotient (x : B.State r) :
    A.ofQuot (B.quotient r x) = A.presentation.toPublicWorld x := by
  apply A.toQuot_injective
  rw [A.toQuot_ofQuot, A.toQuot_mk]

/-- The bare tower's normal form transported through the equivalence to
an endomorphism of the A4 public world. -/
noncomputable def descendedNormalForm :
    A.presentation.PublicWorld → A.presentation.PublicWorld :=
  fun w => A.ofQuot (B.normalForm r (A.toQuot w))

/-- The transported normal form projects back to the bare normal form. -/
@[simp] theorem toQuot_descendedNormalForm (w : A.presentation.PublicWorld) :
    A.toQuot (A.descendedNormalForm w) = B.normalForm r (A.toQuot w) :=
  A.toQuot_ofQuot _

/-- The transported normal form inherits idempotence from the bare
tower law. -/
theorem descendedNormalForm_idempotent (w : A.presentation.PublicWorld) :
    A.descendedNormalForm (A.descendedNormalForm w) =
      A.descendedNormalForm w := by
  unfold descendedNormalForm
  rw [A.toQuot_ofQuot, B.normalForm_idempotent]

end PublicReadbackAdaptor

/-- Obligation-1 surface at a joint instance: a public readback adaptor
on the instance's bare source-link tower at one bare regulator. -/
abbrev JointFiniteInstance.ReadbackAdaptor {κ : Type u'} [Preorder κ]
    {ι : Type u} [Preorder ι] (J : JointFiniteInstance κ ι) (r : ι) :=
  PublicReadbackAdaptor J.bareTower r

/-! ## Endpoint compatibility: the A4 repair bridge -/

/-- An A4 endpoint bridge over a public readback adaptor: a labelled A4
repair system on the adaptor's presentation, the exact
quotient-compatibility packet, and a raw normalizer that is fixed under
every labelled move, reachable by accepted steps, and projects to the
bare tower's declared normal form.  The normalizer replaces a
termination premise: it supplies the completed witness schedule that the
endpoint theorem compares against. -/
structure PublicEndpointBridge {ι : Type u} [Preorder ι]
    (B : BareConsensusTower ι) (r : ι) (A : PublicReadbackAdaptor B r) where
  /-- The labelled A4 repair system on the adaptor's presentation. -/
  system : A.presentation.FiniteRepairSystem
  /-- The A4 quotient-compatibility packet for the repair system. -/
  compat : QuotientCompatible A.presentation system
  /-- The raw normalizer. -/
  normalize : B.State r → B.State r
  /-- The normalizer output is fixed under every labelled move. -/
  normalize_fixed : ∀ x : B.State r,
    IsFixed A.presentation system (normalize x)
  /-- The normalizer output is reachable by accepted steps. -/
  normalize_reachable : ∀ x : B.State r,
    ReflTransGen (step A.presentation system) x (normalize x)
  /-- The normalizer projects to the bare tower's declared normal
  form. -/
  normalize_matches : ∀ x : B.State r,
    B.quotient r (normalize x) = B.normalForm r (B.quotient r x)

namespace PublicEndpointBridge

variable {ι : Type u} [Preorder ι] {B : BareConsensusTower ι} {r : ι}
variable {A : PublicReadbackAdaptor B r}

/-- The bridge supplies a completed finite schedule from every raw
start, so completion is witnessed rather than postulated. -/
theorem completedSchedule_exists (Br : PublicEndpointBridge B r A)
    (start : B.State r) :
    ∃ schedule : List Br.system.Move,
      CompletedSchedule A.presentation Br.system schedule start := by
  obtain ⟨schedule, hrun⟩ :=
    reachable_has_schedule A.presentation Br.system
      (Br.normalize_reachable start)
  refine ⟨schedule, ?_⟩
  show IsFixed A.presentation Br.system
    (run A.presentation Br.system schedule start)
  rw [hrun]
  exact Br.normalize_fixed start

/-- Endpoint compatibility: under declared confluence, every completed
A4 schedule from every raw start lands, through the adaptor, on the
bare tower's normal form of the start's physical quotient class. -/
theorem endpoint_matches (Br : PublicEndpointBridge B r A)
    (hconf : OPH.AbstractRewriting.Confluent (step A.presentation Br.system))
    {start : B.State r} {schedule : List Br.system.Move}
    (hsched : CompletedSchedule A.presentation Br.system schedule start) :
    A.toQuot (A.presentation.toPublicWorld
        (run A.presentation Br.system schedule start)) =
      B.normalForm r (B.quotient r start) := by
  obtain ⟨witness, hwitness⟩ :=
    reachable_has_schedule A.presentation Br.system
      (Br.normalize_reachable start)
  have hdone : CompletedSchedule A.presentation Br.system witness start := by
    show IsFixed A.presentation Br.system
      (run A.presentation Br.system witness start)
    rw [hwitness]
    exact Br.normalize_fixed start
  have hruns :=
    completed_runs_equal A.presentation Br.system hconf hsched hdone
  rw [hruns, hwitness, A.toQuot_mk]
  exact Br.normalize_matches start

/-- Endpoint compatibility inside the A4 public world: the completed A4
public endpoint equals the bare normal form transported through the
equivalence. -/
theorem endpoint_matches_descended (Br : PublicEndpointBridge B r A)
    (hconf : OPH.AbstractRewriting.Confluent (step A.presentation Br.system))
    {start : B.State r} {schedule : List Br.system.Move}
    (hsched : CompletedSchedule A.presentation Br.system schedule start) :
    A.presentation.toPublicWorld
        (run A.presentation Br.system schedule start) =
      A.descendedNormalForm (A.presentation.toPublicWorld start) := by
  apply A.toQuot_injective
  rw [A.toQuot_descendedNormalForm]
  exact Br.endpoint_matches hconf hsched

/-- Completed endpoints are independent of both the schedule and the
raw representative of the start's physical quotient class. -/
theorem endpoints_agree_of_quotient_eq (Br : PublicEndpointBridge B r A)
    (hconf : OPH.AbstractRewriting.Confluent (step A.presentation Br.system))
    {x y : B.State r} (hq : B.quotient r x = B.quotient r y)
    {sx sy : List Br.system.Move}
    (hx : CompletedSchedule A.presentation Br.system sx x)
    (hy : CompletedSchedule A.presentation Br.system sy y) :
    A.toQuot (A.presentation.toPublicWorld
        (run A.presentation Br.system sx x)) =
      A.toQuot (A.presentation.toPublicWorld
        (run A.presentation Br.system sy y)) := by
  rw [Br.endpoint_matches hconf hx, Br.endpoint_matches hconf hy, hq]

end PublicEndpointBridge

/-! ## Obligation-3 fragment: raw refinement maps and naturality laws -/

/-- Raw-configuration refinement maps covering the bare coarse maps:
one raw map per regulator pair, functorial, and commuting with the two
physical quotient maps. -/
structure RawRefinementSystem {ι : Type u} [Preorder ι]
    (B : BareConsensusTower ι) where
  /-- The raw refinement map for one regulator pair, fine to coarse. -/
  rawCoarse : ∀ {r s : ι}, r ≤ s → B.State s → B.State r
  /-- The raw map covers the bare physical coarse map. -/
  rawCoarse_covers : ∀ {r s : ι} (hrs : r ≤ s) (x : B.State s),
    B.quotient r (rawCoarse hrs x) = B.coarse hrs (B.quotient s x)
  /-- Raw refinement is the identity on the reflexive pair. -/
  rawCoarse_refl : ∀ (r : ι) (x : B.State r),
    rawCoarse (show r ≤ r from le_rfl) x = x
  /-- Raw refinement composes along regulator triples. -/
  rawCoarse_trans : ∀ {r s t : ι} (hrs : r ≤ s) (hst : s ≤ t)
    (x : B.State t),
    rawCoarse hrs (rawCoarse hst x) = rawCoarse (hrs.trans hst) x

/-- Readback refinement naturality across one regulator pair, stated as
the structure's law: a declared signature refinement map intertwines the
two adaptor readbacks through the raw refinement map. -/
structure ReadbackRefinementNaturality {ι : Type u} [Preorder ι]
    {B : BareConsensusTower ι} (RS : RawRefinementSystem B)
    {r s : ι} (hrs : r ≤ s)
    (Ar : PublicReadbackAdaptor B r) (As : PublicReadbackAdaptor B s) where
  /-- The declared signature refinement map. -/
  signatureCoarse : As.Signature → Ar.Signature
  /-- Readback naturality: the general form of the obligation-3 law. -/
  readback_natural : ∀ x : B.State s,
    Ar.readback (RS.rawCoarse hrs x) = signatureCoarse (As.readback x)

namespace ReadbackRefinementNaturality

variable {ι : Type u} [Preorder ι] {B : BareConsensusTower ι}
variable {RS : RawRefinementSystem B} {r s : ι} {hrs : r ≤ s}
variable {Ar : PublicReadbackAdaptor B r} {As : PublicReadbackAdaptor B s}

/-- The raw refinement map descends to the A4 public worlds because
readback naturality makes it a congruence for the fine public setoid. -/
def publicCoarse (N : ReadbackRefinementNaturality RS hrs Ar As) :
    As.presentation.PublicWorld → Ar.presentation.PublicWorld :=
  Quotient.lift
    (fun x => Ar.presentation.toPublicWorld (RS.rawCoarse hrs x))
    (fun x y hxy =>
      (Ar.presentation.toPublicWorld_eq_iff _ _).2 (by
        show Ar.readback (RS.rawCoarse hrs x) =
          Ar.readback (RS.rawCoarse hrs y)
        rw [N.readback_natural, N.readback_natural]
        exact congrArg N.signatureCoarse hxy))

@[simp] theorem publicCoarse_mk
    (N : ReadbackRefinementNaturality RS hrs Ar As) (x : B.State s) :
    N.publicCoarse (As.presentation.toPublicWorld x) =
      Ar.presentation.toPublicWorld (RS.rawCoarse hrs x) :=
  rfl

/-- Naturality square, quotient side: the descended public refinement
map commutes with the bare coarse map through the two adaptors. -/
theorem publicCoarse_toQuot
    (N : ReadbackRefinementNaturality RS hrs Ar As)
    (w : As.presentation.PublicWorld) :
    Ar.toQuot (N.publicCoarse w) = B.coarse hrs (As.toQuot w) := by
  refine Quotient.inductionOn w ?_
  intro x
  exact RS.rawCoarse_covers hrs x

/-- Naturality square, signature side: the descended public refinement
map commutes with the descended public signatures. -/
theorem publicCoarse_signature
    (N : ReadbackRefinementNaturality RS hrs Ar As)
    (w : As.presentation.PublicWorld) :
    Ar.presentation.publicSignature (N.publicCoarse w) =
      N.signatureCoarse (As.presentation.publicSignature w) := by
  refine Quotient.inductionOn w ?_
  intro x
  exact N.readback_natural x

end ReadbackRefinementNaturality

/-- Repair refinement naturality across one regulator pair, stated as
the structure's law: a declared move translation intertwines the two
labelled repair systems and the two raw normalizers through the raw
refinement map. -/
structure RepairRefinementNaturality {ι : Type u} [Preorder ι]
    {B : BareConsensusTower ι} (RS : RawRefinementSystem B)
    {r s : ι} (hrs : r ≤ s)
    {Ar : PublicReadbackAdaptor B r} {As : PublicReadbackAdaptor B s}
    (Brr : PublicEndpointBridge B r Ar)
    (Brs : PublicEndpointBridge B s As) where
  /-- The declared move translation, fine to coarse. -/
  moveCoarse : Brs.system.Move → Brr.system.Move
  /-- Repair naturality: the general form of the obligation-3 law. -/
  repair_natural : ∀ (m : Brs.system.Move) (x : B.State s),
    RS.rawCoarse hrs (Brs.system.repair m x) =
      Brr.system.repair (moveCoarse m) (RS.rawCoarse hrs x)
  /-- Normalizer naturality across the regulator pair. -/
  normalize_natural : ∀ x : B.State s,
    RS.rawCoarse hrs (Brs.normalize x) =
      Brr.normalize (RS.rawCoarse hrs x)

namespace RepairRefinementNaturality

variable {ι : Type u} [Preorder ι] {B : BareConsensusTower ι}
variable {RS : RawRefinementSystem B} {r s : ι} {hrs : r ≤ s}
variable {Ar : PublicReadbackAdaptor B r} {As : PublicReadbackAdaptor B s}
variable {Brr : PublicEndpointBridge B r Ar}
variable {Brs : PublicEndpointBridge B s As}

/-- Schedule-level consequence of repair naturality: refining a whole
finite run commutes with running the translated schedule from the
refined start. -/
theorem run_natural (RN : RepairRefinementNaturality RS hrs Brr Brs)
    (schedule : List Brs.system.Move) (x : B.State s) :
    RS.rawCoarse hrs (run As.presentation Brs.system schedule x) =
      run Ar.presentation Brr.system (schedule.map RN.moveCoarse)
        (RS.rawCoarse hrs x) := by
  induction schedule generalizing x with
  | nil => rfl
  | cons m ms ih =>
      show RS.rawCoarse hrs
          (run As.presentation Brs.system ms (Brs.system.repair m x)) = _
      rw [ih (Brs.system.repair m x), RN.repair_natural m x]
      rfl

/-- The refined normalizer projects onto the bare normal form of the
coarse quotient class: normalizer naturality composed with the covering
law and the coarse endpoint law. -/
theorem normalize_quotient_square
    (RN : RepairRefinementNaturality RS hrs Brr Brs) (x : B.State s) :
    B.quotient r (RS.rawCoarse hrs (Brs.normalize x)) =
      B.normalForm r (B.coarse hrs (B.quotient s x)) := by
  rw [RN.normalize_natural, Brr.normalize_matches, RS.rawCoarse_covers]

end RepairRefinementNaturality

/-! ## The identity demo stack glued to the joint witness instance

The committed `demoTower` has identity physical quotient maps, so its
adaptor readback is the identity and its A4 kernel quotient is a
bijective copy of the physical quotient.  The copy repair on raw pairs
realizes the tower's declared normal form as an A4 endpoint, and the
completed endpoint feeds the joint witness instance's source-link
readout. -/

section DemoBridge

/-- The identity public readback adaptor on the demo tower: the demo
physical quotient map is the identity, so both kernel inclusions hold
definitionally. -/
def demoJointAdaptor (r : ℕ) : PublicReadbackAdaptor demoTower r where
  Signature := Bool × Bool
  stateNonempty := ⟨(true, true)⟩
  readback := id
  quotient_eq_of_readback := fun h => h
  readback_eq_of_quotient := fun h => h

/-- The demo readback separates raw states. -/
theorem demoJointAdaptor_readback_nonconstant :
    (demoJointAdaptor 0).readback (true, true) ≠
      (demoJointAdaptor 0).readback (false, false) := by
  intro h
  exact Bool.noConfusion (congrArg Prod.fst h)

/-- The copy repair on the demo tower's raw pairs. -/
def demoCopyRepair (x : Bool × Bool) : Bool × Bool := (x.1, x.1)

/-- The copy repair fixes a raw pair exactly when the pair agrees. -/
theorem demoCopyRepair_fix_iff (x : Bool × Bool) :
    demoCopyRepair x = x ↔ x.1 = x.2 := by
  rcases x with ⟨a, b⟩
  simp [demoCopyRepair, Prod.ext_iff]

/-- The A4 endpoint bridge on the demo adaptor: one labelled copy move,
the quotient-compatibility packet, and the copy normalizer projecting to
the demo tower's declared normal form. -/
def demoJointBridge (r : ℕ) :
    PublicEndpointBridge demoTower r (demoJointAdaptor r) where
  system :=
    { Move := Unit
      moveFintype := inferInstance
      repair := fun _ => demoCopyRepair }
  compat :=
    { repair_respects_public := by
        intro m x y hxy
        have h : x = y := hxy
        subst h
        rfl
      enabled_respects_public := by
        intro m x y hxy
        have h : x = y := hxy
        subst h
        exact Iff.rfl }
  normalize := demoCopyRepair
  normalize_fixed := fun _ _ => rfl
  normalize_reachable := fun x => by
    by_cases h : demoCopyRepair x = x
    · rw [h]
    · exact ReflTransGen.single ⟨(), rfl, h⟩
  normalize_matches := fun _ => rfl

/-- The demo bridge's accepted relation is terminating: the demo
mismatch strictly decreases along accepted copy moves. -/
theorem demoJointBridge_terminating (r : ℕ) :
    OPH.AbstractRewriting.Terminating
      (step (demoJointAdaptor r).presentation (demoJointBridge r).system) :=
  terminating_of_measure _ _ demoMismatch (by
    intro x y h
    obtain ⟨m, hy, hne⟩ := h
    subst hy
    rcases x with ⟨a, b⟩
    by_cases hab : a = b
    · subst hab
      exact absurd rfl hne
    · show demoMismatch (a, a) < demoMismatch (a, b)
      simp [demoMismatch, hab])

/-- The demo bridge's accepted relation is confluent: the single copy
move is deterministic, so Newman's lemma applies. -/
theorem demoJointBridge_confluent (r : ℕ) :
    OPH.AbstractRewriting.Confluent
      (step (demoJointAdaptor r).presentation (demoJointBridge r).system) :=
  OPH.AbstractRewriting.newman_lemma _
    (demoJointBridge_terminating r) (by
      rintro x y z ⟨m, rfl, _⟩ ⟨n, rfl, _⟩
      exact ⟨demoCopyRepair x, ReflTransGen.refl, ReflTransGen.refl⟩)

/-- The obligation-1 adaptor surface inhabited at the joint witness
instance: the demo adaptor typed against the instance's own bare
tower. -/
def jointWitnessReadbackAdaptor (r : ℕ) :
    jointWitnessInstance.ReadbackAdaptor r :=
  demoJointAdaptor r

/-- The endpoint bridge typed against the joint witness instance's bare
tower. -/
def jointWitnessEndpointBridge (r : ℕ) :
    PublicEndpointBridge jointWitnessInstance.bareTower r
      (jointWitnessReadbackAdaptor r) :=
  demoJointBridge r

/-- The A4 public endpoint feeds the joint instance: a completed A4
schedule whose start settles to the selected branch point reads out,
through the adaptor and the instance's source-link readout, the selected
world slice. -/
theorem jointWitnessInstance_endpoint_reads_selected_slice
    {start : Bool × Bool} {schedule : List Unit}
    (hsched : CompletedSchedule (demoJointAdaptor 0).presentation
      (demoJointBridge 0).system schedule start)
    (hsettle : demoTower.normalForm 0 start =
      jointWitnessInstance.square.branch.point 0) :
    jointWitnessInstance.square.readout.read 0
        ((demoJointAdaptor 0).toQuot
          ((demoJointAdaptor 0).presentation.toPublicWorld
            (run (demoJointAdaptor 0).presentation
              (demoJointBridge 0).system schedule start))) =
      jointWitnessInstance.world.slice := by
  rw [PublicEndpointBridge.endpoint_matches (demoJointBridge 0)
    (demoJointBridge_confluent 0) hsched]
  exact jointWitnessInstance.square.readout.read_of_settles_to 0
    ((demoTower.normalForm_idempotent 0 start).trans hsettle)

/-- The endpoint-to-slice law evaluated on explicit data: the broken
start `(true, false)` completes in one copy move, its endpoint projects
to the branch point, and the readout returns the selected slice. -/
theorem jointWitnessInstance_endpoint_reads_selected_slice_eval :
    jointWitnessInstance.square.readout.read 0
        ((demoJointAdaptor 0).toQuot
          ((demoJointAdaptor 0).presentation.toPublicWorld
            (run (demoJointAdaptor 0).presentation
              (demoJointBridge 0).system [()] (true, false)))) =
      jointWitnessInstance.world.slice :=
  jointWitnessInstance_endpoint_reads_selected_slice
    (fun _ => rfl) rfl

end DemoBridge

/-! ## The hidden-bit pair: a genuine kernel quotient

On the demo tower the physical quotient map is the identity, so the
adaptor equivalence is a bijective copy.  The hidden-bit tower is the
smallest self-contained pair where the equivalence does real work: raw
states carry one private bit that the physical quotient and the
nonconstant readback both erase, and the copy repair completes to the
bare normal form while transporting the hidden bit. -/

section HiddenBitBridge

/-- A bare consensus tower whose raw presentation states carry a public
two-patch pair and one hidden bit.  The physical quotient forgets the
hidden bit; mismatch, repair data, normal form, and boundary all live on
the public pair; the coarse maps are the identity. -/
def hiddenBitTower : BareConsensusTower ℕ where
  State := fun _ => (Bool × Bool) × Bool
  stateFinite := fun _ => inferInstance
  Quot := fun _ => Bool × Bool
  quotFinite := fun _ => inferInstance
  Boundary := fun _ => Bool
  quotient := fun _ => Prod.fst
  quotient_surjective := fun _ q => ⟨(q, false), rfl⟩
  mismatch := fun _ => demoMismatch
  step := fun _ x y => demoMismatch x = 1 ∧ demoMismatch y = 0
  normalForm := fun _ x => (x.1, x.1)
  consistent := fun _ x => x.1 = x.2
  boundary := fun _ x => x.1
  consistent_iff := by
    intro r x
    simp only [demoMismatch]
    by_cases h : x.1 = x.2 <;> simp [h]
  step_descends := by
    intro r x y h
    omega
  normalForm_consistent := by
    intro r x
    rfl
  normalForm_idempotent := by
    intro r x
    rfl
  boundary_normalForm := by
    intro r x
    rfl
  coarse := fun _ => id
  coarse_refl := by
    intro r x
    rfl
  coarse_trans := by
    intro r s t hrs hst x
    rfl
  normalForm_natural := by
    intro r s hrs x
    rfl

/-- The hidden-bit adaptor: the readback is the public projection, so
both kernel inclusions hold definitionally while the readback is
nonconstant and noninjective. -/
def hiddenBitAdaptor (r : ℕ) : PublicReadbackAdaptor hiddenBitTower r where
  Signature := Bool × Bool
  stateNonempty := ⟨((true, true), false)⟩
  readback := Prod.fst
  quotient_eq_of_readback := fun h => h
  readback_eq_of_quotient := fun h => h

/-- The hidden-bit readback separates public data. -/
theorem hiddenBitAdaptor_readback_nonconstant :
    (hiddenBitAdaptor 0).readback ((true, true), false) ≠
      (hiddenBitAdaptor 0).readback ((false, false), false) := by
  intro h
  exact Bool.noConfusion (congrArg Prod.fst h)

/-- The kernel quotient does real work: two raw states differing only in
the hidden bit are distinct yet share one A4 public class. -/
theorem hiddenBitAdaptor_collapses_hidden_bit :
    (((true, true), false) : (Bool × Bool) × Bool) ≠ ((true, true), true) ∧
      (hiddenBitAdaptor 0).presentation.toPublicWorld ((true, true), false) =
        (hiddenBitAdaptor 0).presentation.toPublicWorld
          ((true, true), true) :=
  ⟨fun h => Bool.noConfusion (congrArg Prod.snd h),
    ((hiddenBitAdaptor 0).presentation.toPublicWorld_eq_iff _ _).2 rfl⟩

/-- The equivalence is genuine: the raw-to-public quotient map is
noninjective while the descended map onto the bare physical quotient is
a bijection. -/
theorem hiddenBitAdaptor_equivalence_genuine :
    ¬ Function.Injective
        (hiddenBitAdaptor 0).presentation.toPublicWorld ∧
      Function.Bijective (hiddenBitAdaptor 0).toQuot :=
  ⟨fun hinj =>
      hiddenBitAdaptor_collapses_hidden_bit.1
        (hinj hiddenBitAdaptor_collapses_hidden_bit.2),
    (hiddenBitAdaptor 0).toQuot_bijective⟩

/-- The copy repair on hidden-bit raw states: repair the public pair,
transport the hidden bit. -/
def hiddenBitCopyRepair (x : (Bool × Bool) × Bool) : (Bool × Bool) × Bool :=
  ((x.1.1, x.1.1), x.2)

/-- The hidden-bit copy repair fixes a raw state exactly when its public
pair agrees; the hidden bit never obstructs quiescence. -/
theorem hiddenBitCopyRepair_fix_iff (x : (Bool × Bool) × Bool) :
    hiddenBitCopyRepair x = x ↔ x.1.1 = x.1.2 := by
  rcases x with ⟨⟨a, b⟩, c⟩
  simp [hiddenBitCopyRepair, Prod.ext_iff]

/-- The A4 endpoint bridge on the hidden-bit adaptor. -/
def hiddenBitBridge (r : ℕ) :
    PublicEndpointBridge hiddenBitTower r (hiddenBitAdaptor r) where
  system :=
    { Move := Unit
      moveFintype := inferInstance
      repair := fun _ => hiddenBitCopyRepair }
  compat :=
    { repair_respects_public := by
        intro m x y hxy
        have h : x.1 = y.1 := hxy
        show ((x.1.1, x.1.1) : Bool × Bool) = (y.1.1, y.1.1)
        rw [h]
      enabled_respects_public := by
        intro m x y hxy
        have h : x.1 = y.1 := hxy
        show hiddenBitCopyRepair x ≠ x ↔ hiddenBitCopyRepair y ≠ y
        rw [ne_eq, ne_eq, hiddenBitCopyRepair_fix_iff,
          hiddenBitCopyRepair_fix_iff, h] }
  normalize := hiddenBitCopyRepair
  normalize_fixed := fun _ _ => rfl
  normalize_reachable := fun x => by
    by_cases h : hiddenBitCopyRepair x = x
    · rw [h]
    · exact ReflTransGen.single ⟨(), rfl, h⟩
  normalize_matches := fun _ => rfl

/-- The hidden-bit accepted relation is terminating: the public mismatch
strictly decreases along accepted copy moves. -/
theorem hiddenBitBridge_terminating (r : ℕ) :
    OPH.AbstractRewriting.Terminating
      (step (hiddenBitAdaptor r).presentation (hiddenBitBridge r).system) :=
  terminating_of_measure _ _ (fun x => demoMismatch x.1) (by
    intro x y h
    obtain ⟨m, hy, hne⟩ := h
    subst hy
    rcases x with ⟨⟨a, b⟩, c⟩
    by_cases hab : a = b
    · subst hab
      exact absurd rfl hne
    · show demoMismatch (a, a) < demoMismatch (a, b)
      simp [demoMismatch, hab])

/-- The hidden-bit accepted relation is confluent by determinism and
Newman's lemma. -/
theorem hiddenBitBridge_confluent (r : ℕ) :
    OPH.AbstractRewriting.Confluent
      (step (hiddenBitAdaptor r).presentation (hiddenBitBridge r).system) :=
  OPH.AbstractRewriting.newman_lemma _
    (hiddenBitBridge_terminating r) (by
      rintro x y z ⟨m, rfl, _⟩ ⟨n, rfl, _⟩
      exact ⟨hiddenBitCopyRepair x, ReflTransGen.refl, ReflTransGen.refl⟩)

/-- The one-move schedule is completed at the broken start with hidden
bit set. -/
theorem hiddenBitBridge_completed_singleton :
    CompletedSchedule (hiddenBitAdaptor 0).presentation
      (hiddenBitBridge 0).system [()] ((true, false), true) :=
  fun _ => rfl

/-- Endpoint compatibility exercised on explicit hidden-bit data: the
completed endpoint of the broken start projects to the bare normal form
of its quotient class. -/
theorem hiddenBitBridge_endpoint_matches_eval :
    (hiddenBitAdaptor 0).toQuot
        ((hiddenBitAdaptor 0).presentation.toPublicWorld
          (run (hiddenBitAdaptor 0).presentation (hiddenBitBridge 0).system
            [()] ((true, false), true))) =
      hiddenBitTower.normalForm 0
        (hiddenBitTower.quotient 0 ((true, false), true)) :=
  PublicEndpointBridge.endpoint_matches (hiddenBitBridge 0)
    (hiddenBitBridge_confluent 0) hiddenBitBridge_completed_singleton

/-- The evaluated endpoint is the repaired public pair. -/
theorem hiddenBitBridge_endpoint_value :
    hiddenBitTower.normalForm 0
        (hiddenBitTower.quotient 0 ((true, false), true)) =
      (true, true) :=
  rfl

/-- Representative independence across the hidden bit: two distinct raw
starts differing only in the hidden bit complete to one bare endpoint. -/
theorem hiddenBitBridge_endpoints_representative_independent :
    (((true, false), false) : (Bool × Bool) × Bool) ≠
        ((true, false), true) ∧
      (hiddenBitAdaptor 0).toQuot
          ((hiddenBitAdaptor 0).presentation.toPublicWorld
            (run (hiddenBitAdaptor 0).presentation
              (hiddenBitBridge 0).system [()] ((true, false), false))) =
        (hiddenBitAdaptor 0).toQuot
          ((hiddenBitAdaptor 0).presentation.toPublicWorld
            (run (hiddenBitAdaptor 0).presentation
              (hiddenBitBridge 0).system [()] ((true, false), true))) :=
  ⟨fun h => Bool.noConfusion (congrArg Prod.snd h), rfl⟩

end HiddenBitBridge

/-! ## Refinement naturality proved outright on the identity cases -/

section NaturalityDemos

/-- The identity raw refinement system on the demo tower: the demo
coarse maps are the identity, so the covering and functoriality laws
hold definitionally. -/
def demoRawRefinement : RawRefinementSystem demoTower where
  rawCoarse := fun _ => id
  rawCoarse_covers := by
    intro r s hrs x
    rfl
  rawCoarse_refl := fun _ _ => rfl
  rawCoarse_trans := by
    intro r s t hrs hst x
    rfl

/-- Readback refinement naturality on the demo stack, proved outright
for every regulator pair with the identity signature refinement. -/
def demoJointNaturality (r s : ℕ) (hrs : r ≤ s) :
    ReadbackRefinementNaturality demoRawRefinement hrs
      (demoJointAdaptor r) (demoJointAdaptor s) where
  signatureCoarse := id
  readback_natural := fun _ => rfl

/-- The descended public refinement map at the joint witness pair
`0 ≤ 1` acts as the identity on public classes. -/
theorem demoJointNaturality_publicCoarse_mk (x : Bool × Bool) :
    (demoJointNaturality 0 1 (Nat.zero_le 1)).publicCoarse
        ((demoJointAdaptor 1).presentation.toPublicWorld x) =
      (demoJointAdaptor 0).presentation.toPublicWorld x :=
  rfl

/-- The quotient-side naturality square evaluated at the joint witness
instance's regulator pair `0 ≤ 1`. -/
theorem demoJointNaturality_square
    (w : (demoJointAdaptor 1).presentation.PublicWorld) :
    (demoJointAdaptor 0).toQuot
        ((demoJointNaturality 0 1 (Nat.zero_le 1)).publicCoarse w) =
      demoTower.coarse (Nat.zero_le 1) ((demoJointAdaptor 1).toQuot w) :=
  ReadbackRefinementNaturality.publicCoarse_toQuot _ w

/-- Repair refinement naturality on the demo stack, proved outright for
every regulator pair with the identity move translation. -/
def demoJointRepairNaturality (r s : ℕ) (hrs : r ≤ s) :
    RepairRefinementNaturality demoRawRefinement hrs
      (demoJointBridge r) (demoJointBridge s) where
  moveCoarse := id
  repair_natural := fun _ _ => rfl
  normalize_natural := fun _ => rfl

/-- Schedule-level run naturality evaluated at the joint witness
instance's regulator pair `0 ≤ 1`. -/
theorem demoJointRepairNaturality_run (schedule : List Unit)
    (x : Bool × Bool) :
    demoRawRefinement.rawCoarse (Nat.zero_le 1)
        (run (demoJointAdaptor 1).presentation
          (demoJointBridge 1).system schedule x) =
      run (demoJointAdaptor 0).presentation (demoJointBridge 0).system
        (schedule.map
          (demoJointRepairNaturality 0 1 (Nat.zero_le 1)).moveCoarse) x :=
  RepairRefinementNaturality.run_natural
    (demoJointRepairNaturality 0 1 (Nat.zero_le 1)) schedule x

/-- The identity raw refinement system on the hidden-bit tower. -/
def hiddenBitRawRefinement : RawRefinementSystem hiddenBitTower where
  rawCoarse := fun _ => id
  rawCoarse_covers := by
    intro r s hrs x
    rfl
  rawCoarse_refl := fun _ _ => rfl
  rawCoarse_trans := by
    intro r s t hrs hst x
    rfl

/-- Readback refinement naturality on the hidden-bit stack, proved
outright for every regulator pair; the readback is noninjective, so the
descended square carries a genuine quotient on both sides. -/
def hiddenBitNaturality (r s : ℕ) (hrs : r ≤ s) :
    ReadbackRefinementNaturality hiddenBitRawRefinement hrs
      (hiddenBitAdaptor r) (hiddenBitAdaptor s) where
  signatureCoarse := id
  readback_natural := fun _ => rfl

/-- The quotient-side naturality square on the hidden-bit stack at the
regulator pair `0 ≤ 1`. -/
theorem hiddenBitNaturality_square
    (w : (hiddenBitAdaptor 1).presentation.PublicWorld) :
    (hiddenBitAdaptor 0).toQuot
        ((hiddenBitNaturality 0 1 (Nat.zero_le 1)).publicCoarse w) =
      hiddenBitTower.coarse (Nat.zero_le 1)
        ((hiddenBitAdaptor 1).toQuot w) :=
  ReadbackRefinementNaturality.publicCoarse_toQuot _ w

/-- The signature-side naturality square on the hidden-bit stack at the
regulator pair `0 ≤ 1`. -/
theorem hiddenBitNaturality_signature
    (w : (hiddenBitAdaptor 1).presentation.PublicWorld) :
    (hiddenBitAdaptor 0).presentation.publicSignature
        ((hiddenBitNaturality 0 1 (Nat.zero_le 1)).publicCoarse w) =
      (hiddenBitNaturality 0 1 (Nat.zero_le 1)).signatureCoarse
        ((hiddenBitAdaptor 1).presentation.publicSignature w) :=
  ReadbackRefinementNaturality.publicCoarse_signature _ w

end NaturalityDemos

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.PublicReadbackAdaptor.publicSetoid_eq_quotient_kernel
#print axioms OPH.QFT.PublicReadbackAdaptor.toQuot_bijective
#print axioms OPH.QFT.PublicReadbackAdaptor.toQuot_ofQuot
#print axioms OPH.QFT.PublicReadbackAdaptor.ofQuot_toQuot
#print axioms OPH.QFT.PublicReadbackAdaptor.ofQuot_quotient
#print axioms OPH.QFT.PublicReadbackAdaptor.descendedNormalForm_idempotent
#print axioms OPH.QFT.PublicEndpointBridge.completedSchedule_exists
#print axioms OPH.QFT.PublicEndpointBridge.endpoint_matches
#print axioms OPH.QFT.PublicEndpointBridge.endpoint_matches_descended
#print axioms OPH.QFT.PublicEndpointBridge.endpoints_agree_of_quotient_eq
#print axioms OPH.QFT.ReadbackRefinementNaturality.publicCoarse_toQuot
#print axioms OPH.QFT.ReadbackRefinementNaturality.publicCoarse_signature
#print axioms OPH.QFT.RepairRefinementNaturality.run_natural
#print axioms OPH.QFT.RepairRefinementNaturality.normalize_quotient_square
#print axioms OPH.QFT.jointWitnessInstance_endpoint_reads_selected_slice
#print axioms OPH.QFT.jointWitnessInstance_endpoint_reads_selected_slice_eval
#print axioms OPH.QFT.hiddenBitAdaptor_collapses_hidden_bit
#print axioms OPH.QFT.hiddenBitAdaptor_equivalence_genuine
#print axioms OPH.QFT.hiddenBitBridge_endpoint_matches_eval
#print axioms OPH.QFT.hiddenBitBridge_endpoints_representative_independent
#print axioms OPH.QFT.demoJointNaturality_square
#print axioms OPH.QFT.demoJointRepairNaturality_run
#print axioms OPH.QFT.hiddenBitNaturality_square
