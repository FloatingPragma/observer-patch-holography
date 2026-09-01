import QFT.PublicEndpointBridge

/-!
# The F1 inhabited-tower premise row: one enriched common-domain tower

This module addresses row 1 of the Einstein-branch premise matrix of
completion-plan issue `#694`: the inhabited enriched common-domain tower,
owned jointly by the A3, A4, and E2 lanes.  The row closes structurally
when one constructed value carries the finite record and algebra data,
the geometry data, and the readout data over one shared bare consensus
tower.

`EnrichedCommonDomainTower` bundles exactly that: an E2 joint finite
instance (record and algebra data), an E7 event and geometry readout
fragment (geometry data), and an A4 public endpoint bridge with its
readback adaptor (readout data), all typed over the one bare tower that
the joint instance carries.  The field `fragment_carried` pins the
bundled fragment to the fragment of the joint instance's source-link
square, so the geometry channel of the bundle is the channel the square
consumes rather than an unrelated readout.

The carriage theorems read every channel off the one value: the E6
common-origin and readback laws through `joint`, branch matching and
cross-regulator coherence of the fragment through `fragment_carried`,
and the endpoint-to-slice law `endpoint_reads_selected_slice`, which
sends a completed A4 endpoint through the adaptor into the square's
observer readout under the declared confluence premise.

`demoEnrichedTower` inhabits the structure from the committed witnesses:
`jointWitnessInstance`, `demoFragment`, and the demo adaptor and bridge,
all over `demoTower`; the index types compose without re-instantiation.
The evaluated theorems exercise each channel of the one value on
explicit data, and the nontriviality receipt keeps every channel away
from constant adaptors.

## Gap rows

The named Props `ModularReadoutRow`, `StressReadoutRow`,
`EntropyReadoutRow`, `ScaleReadoutRow`, and
`EinsteinReadoutCompletionRow` state what the row lacks for F1
consumption.  Each of the four readout rows is the typed surface of one
missing readout channel; the bare language cannot type-separate the four
targets, which `readout_rows_share_one_surface` records as a definitional
equivalence.  The degeneration receipt
`structuralReadoutRow_demoTower_by_relabeling` proves that the surface is
closable on the demo tower by relabeling the committed event readout, so
closing these Props carries no physical identification and cannot
discharge the matrix rows by itself.  None of the five Props is proved
for the inhabitant here, and no `EinsteinAdmissibleTower` is assembled:
the completion row remains an F1 obligation.

## Negative receipt

`enriched_bare_not_einstein_decisive` proves at every enriched tower over
an `ℕ`-indexed bare tower that no predicate of the bare consensus
language, evaluated on the shared tower, agrees with the Einstein
equation across the geometric extensions of that same tower: every bare
tower admits an Einstein and a non-Einstein extension sharing it as
reduct.  The demo instantiation composes with the committed
counterextensions of `bare_consensus_not_einstein_complete`, which share
the inhabitant's own bare tower as reduct.  Einstein truth is therefore
not inherited from the bare consensus tower alone; it can only enter
through the readout completion and the analytic premises that remain
open.

## Claim boundary

The inhabitant certifies the structural row only: one explicit value
carries the listed interfaces over one shared bare tower at the committed
demo tier.  The modular, stress, entropy, and scale readouts, the
Einstein readout completion, and every analytic and physical premise row
of the `#694` matrix (guarded Ward limit, first law, MaxEnt envelope,
small-ball asymptotics, vacuum reference, universal coupling, absolute
scale, source-causal continuum receipts, clock network) remain open with F1.  No
source produces the bundle, no physical attachment, uniqueness, or
continuum claim is made, and closure of the stated gap Props is necessary
interface surface rather than sufficient row closure.
-/

namespace OPH.QFT

open OPH.Tower
open OPH.Tower.PublicWorldPresentation
open OPH.Tower.PublicWorldPresentation.FiniteRepairSystem
open OPH.EinsteinBranch

universe u u' v w

/-! ## The enriched common-domain tower -/

/-- The F1 row-1 bundle: one constructed value carrying the finite record
and algebra data (the E2 joint finite instance), the geometry data (the
E7 event and geometry readout fragment), and the readout data (the A4
public endpoint bridge with its readback adaptor), all over the one bare
consensus tower that the joint instance carries.

`fragment_carried` pins the bundled fragment to the fragment of the joint
instance's source-link square: the geometry channel of the bundle is the
channel the square consumes.  The bridge lives at one declared bare
regulator; its adaptor identifies the A4 public setoid with the kernel of
the shared tower's physical quotient by construction. -/
structure EnrichedCommonDomainTower (κ : Type u) [Preorder κ]
    (ι : Type u') [Preorder ι] where
  /-- The E2 joint finite instance: record and algebra data over the
  consensus tower, glued to the shared bare tower by its source-link
  square. -/
  joint : JointFiniteInstance.{u, u', v} κ ι
  /-- The E7 event and geometry readout fragment over the shared bare
  tower: the geometry data of the row. -/
  fragment : EventGeometryReadoutFragment.{u'} joint.bareTower
  /-- The bundled fragment is the fragment of the joint instance's
  source-link square. -/
  fragment_carried : fragment = joint.square.fragment
  /-- The bare regulator at which the endpoint bridge lives. -/
  bridgeRegulator : ι
  /-- The A4 public readback adaptor on the shared bare tower at the
  bridge regulator. -/
  adaptor : PublicReadbackAdaptor.{u', w} joint.bareTower bridgeRegulator
  /-- The A4 endpoint bridge over the adaptor: the readout data of the
  row. -/
  bridge : PublicEndpointBridge.{u', w} joint.bareTower bridgeRegulator adaptor

namespace EnrichedCommonDomainTower

variable {κ : Type u} [Preorder κ] {ι : Type u'} [Preorder ι]

/-- The shared bare tower of the bundle. -/
abbrev bareTower (E : EnrichedCommonDomainTower κ ι) :
    BareConsensusTower ι :=
  E.joint.bareTower

/-- Record and algebra carriage: the E6 common-origin conclusion read off
the one value.  For an active visible pair whose chart owner is the
observer's label, the event support lies below the declared observer
region, the packet record lies in the observer's accessible interface,
and the soldering coordinate is a readout of the same packet. -/
theorem common_origin (E : EnrichedCommonDomainTower κ ι)
    {i : E.joint.signature.Chart} {e : E.joint.signature.Event}
    (hi : E.joint.world.activeChart i) (he : E.joint.world.activeEvent e)
    (hv : E.joint.soldering.visible i e)
    (ho : E.joint.world.chartOwner i =
      E.joint.observer.label E.joint.regulator) :
    E.joint.net.regionLE E.joint.regulator
        (E.joint.world.supportRead (E.joint.world.packetOf e))
        (E.joint.cut.observerRegion (E.joint.world.chartOwner i)) ∧
      E.joint.world.recordRead (E.joint.world.packetOf e) ∈
        E.joint.observer.accessible E.joint.regulator ∧
      E.joint.soldering.coordinate i e =
        E.joint.world.geometryRead i (E.joint.world.packetOf e) :=
  E.joint.joint_common_origin hi he hv ho

/-- Record carriage: the observer's self-readback fixes the record of
every active visible pair owned by the observer's label. -/
theorem record_readback (E : EnrichedCommonDomainTower κ ι)
    {i : E.joint.signature.Chart} {e : E.joint.signature.Event}
    (hi : E.joint.world.activeChart i) (he : E.joint.world.activeEvent e)
    (hv : E.joint.soldering.visible i e)
    (ho : E.joint.world.chartOwner i =
      E.joint.observer.label E.joint.regulator) :
    E.joint.observer.readback E.joint.regulator
        (E.joint.world.recordRead (E.joint.world.packetOf e)) =
      E.joint.world.recordRead (E.joint.world.packetOf e) :=
  E.joint.joint_record_readback hi he hv ho

/-- Geometry carriage, event side: through `fragment_carried` the bundled
fragment's event packet at the square's branch is the encoded packet
assignment of the selected world.  The pin identifies the two packet
codomains only propositionally, so the generic statement is
heterogeneous; the inhabitant below discharges the homogeneous form,
where the codomains agree definitionally. -/
theorem fragment_event_matches (E : EnrichedCommonDomainTower κ ι)
    (r : ι) :
    HEq (E.fragment.eventRead r (E.joint.square.branch.point r))
      (E.joint.square.eventEncode r E.joint.world.packetOf) := by
  rw [E.fragment_carried]
  exact heq_of_eq (E.joint.square.event_matches r)

/-- Geometry carriage, coordinate side: through `fragment_carried` the
bundled fragment's geometry packet at the square's branch is the encoded
chart geometry readout of the selected world. -/
theorem fragment_geometry_matches (E : EnrichedCommonDomainTower κ ι)
    (r : ι) :
    HEq (E.fragment.geometryRead r (E.joint.square.branch.point r))
      (E.joint.square.geometryEncode r E.joint.world.geometryRead) := by
  rw [E.fragment_carried]
  exact heq_of_eq (E.joint.square.geometry_matches r)

/-- Geometry carriage, cross-regulator coherence: on the square's branch
the bundled fragment's event packets form a coherent family over the
packet coarse maps. -/
theorem fragment_event_coherent (E : EnrichedCommonDomainTower κ ι)
    {r s : ι} (hrs : r ≤ s) :
    E.fragment.eventCoarse hrs
        (E.fragment.eventRead s (E.joint.square.branch.point s)) =
      E.fragment.eventRead r (E.joint.square.branch.point r) :=
  E.fragment.eventRead_point_coherent E.joint.square.branch hrs

/-- Geometry carriage, cross-regulator coherence of the geometry
packets. -/
theorem fragment_geometry_coherent (E : EnrichedCommonDomainTower κ ι)
    {r s : ι} (hrs : r ≤ s) :
    E.fragment.geometryCoarse hrs
        (E.fragment.geometryRead s (E.joint.square.branch.point s)) =
      E.fragment.geometryRead r (E.joint.square.branch.point r) :=
  E.fragment.geometryRead_point_coherent E.joint.square.branch hrs

/-- Readout carriage: under declared confluence, a completed A4 schedule
whose start settles to the square's branch point reads out, through the
bundled adaptor and the square's observer readout, the selected world
slice.  This is the endpoint-to-slice law of the one value: the bridge's
endpoint agrees with the bare normal form, and the settling basin of the
branch reproduces the slice. -/
theorem endpoint_reads_selected_slice (E : EnrichedCommonDomainTower κ ι)
    (hconf : OPH.AbstractRewriting.Confluent
      (step E.adaptor.presentation E.bridge.system))
    {start : E.joint.bareTower.State E.bridgeRegulator}
    {schedule : List E.bridge.system.Move}
    (hsched : CompletedSchedule E.adaptor.presentation E.bridge.system
      schedule start)
    (hsettle : E.joint.bareTower.normalForm E.bridgeRegulator
        (E.joint.bareTower.quotient E.bridgeRegulator start) =
      E.joint.square.branch.point E.bridgeRegulator) :
    E.joint.square.readout.read E.bridgeRegulator
        (E.adaptor.toQuot (E.adaptor.presentation.toPublicWorld
          (run E.adaptor.presentation E.bridge.system schedule start))) =
      E.joint.world.slice := by
  rw [E.bridge.endpoint_matches hconf hsched]
  exact E.joint.square.readout.read_of_settles_to E.bridgeRegulator
    ((E.joint.bareTower.normalForm_idempotent E.bridgeRegulator
      (E.joint.bareTower.quotient E.bridgeRegulator start)).trans hsettle)

end EnrichedCommonDomainTower

/-! ## The inhabitant from the committed witnesses -/

/-- The row-1 inhabitant: the committed E2 joint witness instance, the
committed E7 demo fragment, and the committed A4 demo adaptor and bridge,
all over the one shared bare tower `demoTower` at bridge regulator `0`.
The fragment pinning holds definitionally, since the joint witness
square carries `demoFragment`.  The index types compose directly, so no
re-instantiation is needed. -/
noncomputable def demoEnrichedTower : EnrichedCommonDomainTower Unit ℕ where
  joint := jointWitnessInstance
  fragment := demoFragment
  fragment_carried := rfl
  bridgeRegulator := 0
  adaptor := jointWitnessReadbackAdaptor 0
  bridge := jointWitnessEndpointBridge 0

/-- The inhabitant's shared bare tower is the committed demo tower. -/
theorem demoEnrichedTower_bareTower_eq :
    demoEnrichedTower.joint.bareTower = demoTower :=
  rfl

/-- Record and algebra carriage evaluated on every witness pair of the
one value: support containment, accessible-record membership, and the
packet coordinate law hold across the whole witness world. -/
theorem demoEnrichedTower_common_origin (i e : Fin 2) :
    demoEnrichedTower.joint.net.regionLE demoEnrichedTower.joint.regulator
        (demoEnrichedTower.joint.world.supportRead
          (demoEnrichedTower.joint.world.packetOf e))
        (demoEnrichedTower.joint.cut.observerRegion
          (demoEnrichedTower.joint.world.chartOwner i)) ∧
      demoEnrichedTower.joint.world.recordRead
          (demoEnrichedTower.joint.world.packetOf e) ∈
        demoEnrichedTower.joint.observer.accessible
          demoEnrichedTower.joint.regulator ∧
      demoEnrichedTower.joint.soldering.coordinate i e =
        demoEnrichedTower.joint.world.geometryRead i
          (demoEnrichedTower.joint.world.packetOf e) :=
  jointWitnessInstance_common_origin i e

/-- Record carriage evaluated on every witness pair of the one value. -/
theorem demoEnrichedTower_record_readback (i e : Fin 2) :
    demoEnrichedTower.joint.observer.readback
        demoEnrichedTower.joint.regulator
        (demoEnrichedTower.joint.world.recordRead
          (demoEnrichedTower.joint.world.packetOf e)) =
      demoEnrichedTower.joint.world.recordRead
        (demoEnrichedTower.joint.world.packetOf e) :=
  jointWitnessInstance_record_readback i e

/-- Geometry carriage evaluated at the inhabitant, in homogeneous form:
the bundled fragment's packet codomains agree definitionally with the
square's, and the fragment matches the encoded selected world data at the
square's branch at every bare regulator. -/
theorem demoEnrichedTower_fragment_matches (r : ℕ) :
    demoEnrichedTower.fragment.eventRead r
        (demoEnrichedTower.joint.square.branch.point r) =
      demoEnrichedTower.joint.square.eventEncode r
        demoEnrichedTower.joint.world.packetOf ∧
      demoEnrichedTower.fragment.geometryRead r
          (demoEnrichedTower.joint.square.branch.point r) =
        demoEnrichedTower.joint.square.geometryEncode r
          demoEnrichedTower.joint.world.geometryRead :=
  ⟨demoEnrichedTower.joint.square.event_matches r,
    demoEnrichedTower.joint.square.geometry_matches r⟩

/-- Readout carriage evaluated at the inhabitant: the broken start
`(true, false)` completes in one copy move, its A4 endpoint projects to
the branch point, and the square's observer readout returns the selected
world slice. -/
theorem demoEnrichedTower_endpoint_reads_selected_slice :
    demoEnrichedTower.joint.square.readout.read 0
        (demoEnrichedTower.adaptor.toQuot
          (demoEnrichedTower.adaptor.presentation.toPublicWorld
            (run demoEnrichedTower.adaptor.presentation
              demoEnrichedTower.bridge.system [()] (true, false)))) =
      demoEnrichedTower.joint.world.slice :=
  demoEnrichedTower.endpoint_reads_selected_slice
    (demoJointBridge_confluent 0) (fun _ => rfl) rfl

/-- Nontriviality receipt for the inhabitant, one clause per channel: the
square's event encoding separates the selected world packets from a
constant packet map, the bundled fragment separates two consistent
quotient states, the bundled adaptor readback separates two raw states,
and the world record readout separates the two packets. -/
theorem demoEnrichedTower_nontrivial :
    demoEnrichedTower.joint.square.eventEncode 0
        demoEnrichedTower.joint.world.packetOf ≠
      demoEnrichedTower.joint.square.eventEncode 0 (fun _ => (0 : Fin 2)) ∧
      demoEnrichedTower.fragment.eventRead 0 (true, true) ≠
        demoEnrichedTower.fragment.eventRead 0 (false, false) ∧
      demoEnrichedTower.adaptor.readback (true, true) ≠
        demoEnrichedTower.adaptor.readback (false, false) ∧
      demoEnrichedTower.joint.world.recordRead (0 : Fin 2) ≠
        demoEnrichedTower.joint.world.recordRead (1 : Fin 2) :=
  ⟨jointWitnessInstance_eventEncode_separates,
    OPH.Tower.demoFragment_eventRead_separates,
    demoJointAdaptor_readback_nonconstant,
    jointWitnessInstance_record_separates⟩

/-! ## Gap rows: what the row lacks for F1 consumption

The Props below are stated and left open.  Proving them supplies typed
surface only; the F1 matrix rows demand source-derived readouts with
physical identification, which the bare language does not express.  The
degeneration receipt makes that limit explicit. -/

/-- A structural readout row over a bare consensus tower: a per-regulator
carrier family, a normal-form-invariant readout on every physical
quotient state, refinement-natural coarse maps, and a separation receipt
on the consistent sector.  This is the typed surface of one missing
readout channel of `CommonReadoutTower`; it carries no physical
identification. -/
def StructuralReadoutRow {ι : Type u'} [Preorder ι]
    (B : BareConsensusTower ι) : Prop :=
  ∃ (Carrier : ι → Type) (read : ∀ r, B.Quot r → Carrier r)
    (coarsen : ∀ r s : ι, r ≤ s → Carrier s → Carrier r),
    (∀ r q, read r (B.normalForm r q) = read r q) ∧
      (∀ (r s : ι) (hrs : r ≤ s) (q : B.Quot s),
        coarsen r s hrs (read s q) = read r (B.coarse hrs q)) ∧
      ∃ (r : ι) (q q' : B.Quot r),
        B.consistent r q ∧ B.consistent r q' ∧ read r q ≠ read r q'

/-- The modular readout the row lacks: the typed surface of the modular
channel of the Einstein readout completion on the shared bare tower.  The
name records the intended physical target only; the bare language cannot
type-separate it from the other missing channels. -/
def ModularReadoutRow {κ : Type u} [Preorder κ] {ι : Type u'} [Preorder ι]
    (E : EnrichedCommonDomainTower κ ι) : Prop :=
  StructuralReadoutRow E.joint.bareTower

/-- The stress readout the row lacks: the typed surface of the stress
channel on the shared bare tower. -/
def StressReadoutRow {κ : Type u} [Preorder κ] {ι : Type u'} [Preorder ι]
    (E : EnrichedCommonDomainTower κ ι) : Prop :=
  StructuralReadoutRow E.joint.bareTower

/-- The entropy readout the row lacks: the typed surface of the entropy
channel on the shared bare tower. -/
def EntropyReadoutRow {κ : Type u} [Preorder κ] {ι : Type u'} [Preorder ι]
    (E : EnrichedCommonDomainTower κ ι) : Prop :=
  StructuralReadoutRow E.joint.bareTower

/-- The scale readout the row lacks: the typed surface of the scale
channel on the shared bare tower. -/
def ScaleReadoutRow {κ : Type u} [Preorder κ] {ι : Type u'} [Preorder ι]
    (E : EnrichedCommonDomainTower κ ι) : Prop :=
  StructuralReadoutRow E.joint.bareTower

/-- The Einstein readout completion the row lacks: a six-readout
common-domain tower over the shared bare tower carrying the bundled
fragment's event and geometry readouts through injective encodings at
every regulator, with typed-arrow commutation, refinement naturality, and
the boundary fibre.  A witness would assemble an
`EinsteinAdmissibleTower` on the shared tower; every analytic Einstein
premise of `composedEinsteinBranch` remains a separate argument even
then. -/
def EinsteinReadoutCompletionRow {κ : Type u} [Preorder κ]
    {ι : Type u'} [Preorder ι]
    (E : EnrichedCommonDomainTower κ ι) : Prop :=
  ∃ R : CommonReadoutTower E.joint.bareTower,
    TypedArrowCommutation R ∧
      ReadoutNaturality R ∧
      BoundaryFiber E.joint.bareTower ∧
      (∃ evEnc : ∀ r : ι, E.fragment.EventPacket r → R.Event r,
        (∀ r, Function.Injective (evEnc r)) ∧
          ∀ r q, R.eventRead r q = evEnc r (E.fragment.eventRead r q)) ∧
      (∃ geoEnc : ∀ r : ι, E.fragment.GeometryPacket r → R.Geometry r,
        (∀ r, Function.Injective (geoEnc r)) ∧
          ∀ r q, R.geometryRead r q = geoEnc r (E.fragment.geometryRead r q))

/-- The four readout gap rows are one structural surface: the bare
language does not distinguish the modular, stress, entropy, and scale
targets.  Physical identification is what separates them, and it lives
outside this language with F1. -/
theorem readout_rows_share_one_surface {κ : Type u} [Preorder κ]
    {ι : Type u'} [Preorder ι] (E : EnrichedCommonDomainTower κ ι) :
    (ModularReadoutRow E ↔ StressReadoutRow E) ∧
      (StressReadoutRow E ↔ EntropyReadoutRow E) ∧
      (EntropyReadoutRow E ↔ ScaleReadoutRow E) :=
  ⟨Iff.rfl, Iff.rfl, Iff.rfl⟩

/-- Degeneration receipt for the readout surfaces: on the demo tower the
structural readout row is closable by relabeling the committed protected
bit as the readout.  A proof of the four gap rows by this route carries
no physical content, so the matrix rows for the modular, stress, entropy,
and scale readouts cannot be discharged by the structural Props alone. -/
theorem structuralReadoutRow_demoTower_by_relabeling :
    StructuralReadoutRow demoTower :=
  ⟨fun _ => Bool, fun _ q => q.1, fun _ _ _ => id,
    fun _ _ => rfl, fun _ _ _ _ => rfl,
    0, (true, true), (false, false), rfl, rfl,
    fun h => Bool.noConfusion h⟩

/-! ## Negative receipt: bare consensus does not inherit Einstein truth -/

/-- An Einstein geometric extension sharing an arbitrary bare tower as
reduct: one point, unit metric, zero curvature and stress, zero
cosmological term.  The construction uses no property of the tower. -/
def einsteinExtensionOf (B : BareConsensusTower ℕ) : GeometricExtension where
  reduct := B
  Point := Unit
  point_nonempty := ⟨()⟩
  metric := fun _ => 1
  curvature := fun _ => 0
  stress := fun _ => 0
  cosmological := 0
  coupling := 1

/-- A non-Einstein geometric extension sharing the same arbitrary bare
tower as reduct: identical except for unit curvature. -/
def nonEinsteinExtensionOf (B : BareConsensusTower ℕ) :
    GeometricExtension where
  reduct := B
  Point := Unit
  point_nonempty := ⟨()⟩
  metric := fun _ => 1
  curvature := fun _ => 1
  stress := fun _ => 0
  cosmological := 0
  coupling := 1

/-- The Einstein extension of any bare tower satisfies the Einstein
equation. -/
theorem einsteinEq_einsteinExtensionOf (B : BareConsensusTower ℕ) :
    EinsteinEq (einsteinExtensionOf B) := by
  intro p
  simp [einsteinExtensionOf]

/-- The non-Einstein extension of any bare tower violates the Einstein
equation. -/
theorem not_einsteinEq_nonEinsteinExtensionOf (B : BareConsensusTower ℕ) :
    ¬ EinsteinEq (nonEinsteinExtensionOf B) := by
  intro h
  have hp := h ()
  change (1 : ℤ) = 0 at hp
  omega

/-- Negative receipt at every enriched tower: the shared bare tower alone
does not decide that Einstein truth is inherited.  Even restricted to
geometric extensions whose reduct is the bundle's own bare tower, no
predicate of the bare consensus language agrees with the Einstein
equation, because that one tower carries both an Einstein and a
non-Einstein extension. -/
theorem enriched_bare_not_einstein_decisive {κ : Type u} [Preorder κ]
    (E : EnrichedCommonDomainTower κ ℕ) :
    ¬ ∃ decide : BareConsensusTower ℕ → Prop,
        ∀ X : GeometricExtension, X.reduct = E.joint.bareTower →
          (decide X.reduct ↔ EinsteinEq X) := by
  rintro ⟨decide, h⟩
  have h1 := h (einsteinExtensionOf E.joint.bareTower) rfl
  have h2 := h (nonEinsteinExtensionOf E.joint.bareTower) rfl
  exact not_einsteinEq_nonEinsteinExtensionOf E.joint.bareTower
    (h2.1 (h1.2 (einsteinEq_einsteinExtensionOf E.joint.bareTower)))

/-- The negative receipt instantiated at the inhabitant. -/
theorem demoEnrichedTower_not_einstein_decisive :
    ¬ ∃ decide : BareConsensusTower ℕ → Prop,
        ∀ X : GeometricExtension,
          X.reduct = demoEnrichedTower.joint.bareTower →
          (decide X.reduct ↔ EinsteinEq X) :=
  enriched_bare_not_einstein_decisive demoEnrichedTower

/-- The committed counterextensions of the non-entailment theorem share
the inhabitant's own bare tower as reduct. -/
theorem demoEnrichedTower_counterextensions_share_bareTower :
    demoEinsteinExtension.reduct = demoEnrichedTower.joint.bareTower ∧
      demoNonEinsteinExtension.reduct = demoEnrichedTower.joint.bareTower :=
  ⟨rfl, rfl⟩

/-- The committed non-entailment theorem composed at the inhabitant: no
bare-language predicate evaluated on the inhabitant's shared tower agrees
with the Einstein equation on both committed counterextensions of that
very tower.  This restates `bare_consensus_not_einstein_complete` at this
instance through `einsteinEq_demoEinsteinExtension` and
`not_einsteinEq_demoNonEinsteinExtension`. -/
theorem demoEnrichedTower_bare_not_einstein_complete :
    ¬ ∃ decide : BareConsensusTower ℕ → Prop,
        (decide demoEnrichedTower.joint.bareTower ↔
          EinsteinEq demoEinsteinExtension) ∧
        (decide demoEnrichedTower.joint.bareTower ↔
          EinsteinEq demoNonEinsteinExtension) := by
  rintro ⟨decide, h1, h2⟩
  exact not_einsteinEq_demoNonEinsteinExtension
    (h2.1 (h1.2 einsteinEq_demoEinsteinExtension))

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.EnrichedCommonDomainTower.common_origin
#print axioms OPH.QFT.EnrichedCommonDomainTower.record_readback
#print axioms OPH.QFT.EnrichedCommonDomainTower.fragment_event_matches
#print axioms OPH.QFT.EnrichedCommonDomainTower.fragment_geometry_matches
#print axioms OPH.QFT.EnrichedCommonDomainTower.fragment_event_coherent
#print axioms OPH.QFT.EnrichedCommonDomainTower.fragment_geometry_coherent
#print axioms OPH.QFT.EnrichedCommonDomainTower.endpoint_reads_selected_slice
#print axioms OPH.QFT.demoEnrichedTower_bareTower_eq
#print axioms OPH.QFT.demoEnrichedTower_common_origin
#print axioms OPH.QFT.demoEnrichedTower_record_readback
#print axioms OPH.QFT.demoEnrichedTower_fragment_matches
#print axioms OPH.QFT.demoEnrichedTower_endpoint_reads_selected_slice
#print axioms OPH.QFT.demoEnrichedTower_nontrivial
#print axioms OPH.QFT.readout_rows_share_one_surface
#print axioms OPH.QFT.structuralReadoutRow_demoTower_by_relabeling
#print axioms OPH.QFT.einsteinEq_einsteinExtensionOf
#print axioms OPH.QFT.not_einsteinEq_nonEinsteinExtensionOf
#print axioms OPH.QFT.enriched_bare_not_einstein_decisive
#print axioms OPH.QFT.demoEnrichedTower_not_einstein_decisive
#print axioms OPH.QFT.demoEnrichedTower_counterextensions_share_bareTower
#print axioms OPH.QFT.demoEnrichedTower_bare_not_einstein_complete
