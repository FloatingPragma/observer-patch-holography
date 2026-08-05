import QFT.ObserverEventWorld
import Time.ObserverHistory

/-!
# Observer event towers: precedence, selected refinement, raw-event descent

This module implements deliverables 9, 10, and 11 of completion-plan issue
`#712` over the event worlds of `QFT/ObserverEventWorld.lean`.

`EventPrecedenceAdapter` decorates one event world with an A1
`ObserverRecordOrder` on event slots and a partial record-label adapter.
The consistency law ties the two orders together: on active visible labeled
pairs, event precedence implies the chart owner's record-order precedence
of the corresponding record labels.

`SelectedEventTowerRefinement` packages one selected coarse-to-fine
refinement between two event worlds over a regulator pair `r ≤ s` of one
tower/net family: slot maps for events, charts, and packets, naturality
fields for the active masks, chart ownership, observer regions, chart
regions, event supports, packets, record readouts, and visibility, and a
declared recharting map for the geometry.  The composite theorem
`refine_common_origin` transports the common-origin conclusion: at the
refined image of coarse data, support containment, record membership, and
the recharted coordinate law hold in the fine world.

`RawEventCandidates` supplies the pre-event layer: raw candidates with a
declared semantic-coincidence setoid, the event quotient, and the descent
theorems showing that the induced packet map and every packet-derived
readout are well defined on the quotient.  The C2 three-germ control
`OPH.C2Soldering.overlapControl` connects negatively: a reflexive symmetric
relation without transitivity admits no coincidence setoid, so the quotient
constructor is closed to it.

Claim boundary: every structure is a finite conditional interface over the
committed access-cut and event-world packets.  Event precedence is order
data under the A1 type ledger; it carries no repair schedule, modular
parameter, clock reading, proper time, or global time function.  Refinement
data are declared maps; no continuum limit, convergence, spacetime, or
physical coverage is claimed, and no theorem grants any observer access to
the ambient private algebra.
-/

namespace OPH.QFT

open OPH.Tower
open OPH.TimeOrderLedger

universe u v

variable {ι : Type u} [Preorder ι]
variable {T : ConsensusTower ι} {N : FiniteCausalObserverNet T}

/-! ## Deliverable 9: event precedence over one event world -/

section Precedence

variable {r : ι} {A : ObserverAccessCut T N r} {sig : EventWorldSignature}
variable {S : OPH.C2Soldering.EventFrameSoldering sig.Event sig.Chart}

/-- An event precedence adapter over one observer quantum-event world.

The precedence carrier is the A1 `ObserverRecordOrder` type of the
time-and-order type ledger, so precedence is strict order data on event
slots and nothing else.  Under the A1 ledger discipline no identification
with any other ledger layer is included: the precedence is neither a repair
schedule (`RepairOrder`), a modular parameter (`ModularParameter`), a clock
reading (`ClockReadout`), a proper time (`ProperTime`), nor a global time
function (`GlobalTimeFunction`), and none of those types is produced or
consumed here.

`labeled` and `recordLabel` form the partial record-label adapter: a
declared mask of events carrying a record label of the tower.  The
consistency law scopes to active visible labeled pairs under one chart and
routes through the chart owner: event precedence implies the owner's
record-order precedence of the corresponding record labels. -/
structure EventPrecedenceAdapter (W : FiniteEventGeometryWorld A sig S) where
  /-- The declared strict precedence on event slots, carried by the A1
  order type. -/
  eventOrder : ObserverRecordOrder sig.Event
  /-- The declared record-label mask: which events carry a record label. -/
  labeled : sig.Event → Prop
  /-- The declared record label of a labeled event. -/
  recordLabel : ∀ e, labeled e → T.Record r
  /-- Consistency with the owner's record order: on an active visible
  labeled pair, event precedence implies record-order precedence of the
  labels under the chart owner's declared order. -/
  precedes_consistent : ∀ {i : sig.Chart} {e f : sig.Event},
    W.activeChart i → W.activeEvent e → W.activeEvent f →
    S.visible i e → S.visible i f →
    ∀ (hle : labeled e) (hlf : labeled f),
      eventOrder.precedes e f →
        (T.recordOrder r (W.chartOwner i)).precedes
          (recordLabel e hle) (recordLabel f hlf)

namespace EventPrecedenceAdapter

variable {W : FiniteEventGeometryWorld A sig S}

/-- The declared event precedence is asymmetric: irreflexivity and
transitivity of the A1 order exclude a two-way pair. -/
theorem precedes_asymm (P : EventPrecedenceAdapter W) {e f : sig.Event}
    (h : P.eventOrder.precedes e f) : ¬ P.eventOrder.precedes f e :=
  fun h' => P.eventOrder.irrefl e (P.eventOrder.trans h h')

end EventPrecedenceAdapter

end Precedence

/-! ### A witness adapter on the dimension-two event world -/

section PrecedenceWitness

open OPH.D1Time

/-- The single record label of the one-block witness tower: the sure event
of the one-block partition, certified nonzero. -/
noncomputable def witnessRecordLabel : witnessTower.Record () :=
  ⟨0, by
    intro h
    have h00 := congrFun (congrFun h 0) 0
    simp [oneBlockPartition] at h00⟩

/-- The witness precedence adapter on the dimension-two event world: the
canonical strict chain on the two event slots, with the record label
declared on the first event only.  The witness tower carries the discrete
record order, whose precedence is empty, so a labeled strictly ordered pair
would contradict the consistency law; the single-label mask keeps the
adapter exact while the event precedence stays genuinely nonempty. -/
noncomputable def witnessPrecedenceAdapter :
    EventPrecedenceAdapter witnessEventWorld where
  eventOrder := finRecordOrder 2
  labeled := fun e => e = 0
  recordLabel := fun _ _ => witnessRecordLabel
  precedes_consistent := by
    intro i e f _ _ _ _ _ hle hlf hp
    rw [show e = 0 from hle, show f = 0 from hlf] at hp
    exact absurd hp ((finRecordOrder 2).irrefl 0)

/-- The witness precedence is genuinely nonempty and one-way: the first
event strictly precedes the second, the reverse pair fails, and the first
event carries the declared record label. -/
theorem witnessPrecedenceAdapter_strict :
    witnessPrecedenceAdapter.eventOrder.precedes (0 : Fin 2) (1 : Fin 2) ∧
      ¬ witnessPrecedenceAdapter.eventOrder.precedes (1 : Fin 2) (0 : Fin 2) ∧
      witnessPrecedenceAdapter.labeled (0 : Fin 2) := by
  refine ⟨?_, ?_, rfl⟩
  · show (0 : Fin 2) < 1
    decide
  · show ¬ (1 : Fin 2) < 0
    decide

end PrecedenceWitness

/-! ## Deliverable 10: selected coarse-to-fine refinement -/

section Refinement

variable {r s : ι}

/-- One selected coarse-to-fine refinement between two event worlds over
the regulator pair `r ≤ s` of a single tower/net family.

The slot maps carry events, charts, and packets from the coarse signature
to the fine one.  The naturality fields make every declared readout commute
with the tower refinement maps: the active masks and visibility are
preserved, chart ownership refines through `observerRefine`, observer
regions, chart regions, and event supports refine through `regionRefine`,
`packetOf` commutes with the packet map, and the fine record of a mapped
packet is the `algebraRefine` image of the coarse record.

Geometry refinement carries a declared recharting map in place of literal
coordinate equality.  The identity
`Sf.coordinate (chartMap i) (eventMap e) = Sc.coordinate i e` is excluded
by design: the two solderings carry independent affine and Lorentz chart
gauges, so a refinement is free to re-express every coordinate.  The
`rechart` field names that re-expression, and `coordinate_rechart` is the
only geometric compatibility imposed. -/
structure SelectedEventTowerRefinement (hrs : r ≤ s)
    {Ac : ObserverAccessCut T N r} {Af : ObserverAccessCut T N s}
    {sigc sigf : EventWorldSignature}
    {Sc : OPH.C2Soldering.EventFrameSoldering sigc.Event sigc.Chart}
    {Sf : OPH.C2Soldering.EventFrameSoldering sigf.Event sigf.Chart}
    (Wc : FiniteEventGeometryWorld Ac sigc Sc)
    (Wf : FiniteEventGeometryWorld Af sigf Sf) where
  /-- The event slot map. -/
  eventMap : sigc.Event → sigf.Event
  /-- The chart slot map. -/
  chartMap : sigc.Chart → sigf.Chart
  /-- The packet slot map. -/
  packetMap : sigc.Packet → sigf.Packet
  /-- Active events map to active events. -/
  activeEvent_natural : ∀ e, Wc.activeEvent e → Wf.activeEvent (eventMap e)
  /-- Active charts map to active charts. -/
  activeChart_natural : ∀ i, Wc.activeChart i → Wf.activeChart (chartMap i)
  /-- Visibility is preserved on mapped pairs. -/
  visible_natural : ∀ {i e}, Sc.visible i e →
    Sf.visible (chartMap i) (eventMap e)
  /-- Chart ownership refines through the declared observer refinement. -/
  chartOwner_natural : ∀ i,
    Wf.chartOwner (chartMap i) = T.observerRefine hrs (Wc.chartOwner i)
  /-- Observer support regions refine through the declared region
  refinement. -/
  observerRegion_natural : ∀ o : T.Observer r,
    Af.observerRegion (T.observerRefine hrs o) =
      N.regionRefine hrs (Ac.observerRegion o)
  /-- Chart support regions refine through the declared region
  refinement. -/
  chartRegion_natural : ∀ i,
    Wf.chartRegion (chartMap i) = N.regionRefine hrs (Wc.chartRegion i)
  /-- The packet of a mapped event is the mapped packet. -/
  packetOf_natural : ∀ e,
    Wf.packetOf (eventMap e) = packetMap (Wc.packetOf e)
  /-- Event support readouts commute with the region refinement through the
  packet map. -/
  supportRead_natural : ∀ p,
    Wf.supportRead (packetMap p) = N.regionRefine hrs (Wc.supportRead p)
  /-- The `algebraRefine` image of the coarse record is the fine record of
  the mapped packet. -/
  recordRead_natural : ∀ p,
    T.algebraRefine hrs (Wc.recordRead p) = Wf.recordRead (packetMap p)
  /-- The declared recharting of coarse chart coordinates into the fine
  soldering. -/
  rechart : sigc.Chart → OPH.C1Lorentz.Herm2 → OPH.C1Lorentz.Herm2
  /-- On active visible coarse pairs the fine coordinate of the mapped pair
  is the recharted coarse coordinate. -/
  coordinate_rechart : ∀ {i e}, Wc.activeChart i → Wc.activeEvent e →
    Sc.visible i e →
      Sf.coordinate (chartMap i) (eventMap e) = rechart i (Sc.coordinate i e)

namespace SelectedEventTowerRefinement

variable {hrs : r ≤ s}
variable {Ac : ObserverAccessCut T N r} {Af : ObserverAccessCut T N s}
variable {sigc sigf : EventWorldSignature}
variable {Sc : OPH.C2Soldering.EventFrameSoldering sigc.Event sigc.Chart}
variable {Sf : OPH.C2Soldering.EventFrameSoldering sigf.Event sigf.Chart}
variable {Wc : FiniteEventGeometryWorld Ac sigc Sc}
variable {Wf : FiniteEventGeometryWorld Af sigf Sf}

/-- A coarse active visible pair maps to a fine active visible pair, so the
fine common-origin theorem applies at the mapped data. -/
theorem mapped_active_visible (R : SelectedEventTowerRefinement hrs Wc Wf)
    {i : sigc.Chart} {e : sigc.Event} (hi : Wc.activeChart i)
    (he : Wc.activeEvent e) (hv : Sc.visible i e) :
    Wf.activeChart (R.chartMap i) ∧ Wf.activeEvent (R.eventMap e) ∧
      Sf.visible (R.chartMap i) (R.eventMap e) :=
  ⟨R.activeChart_natural i hi, R.activeEvent_natural e he,
    R.visible_natural hv⟩

/-- Composite transport theorem: the common-origin conclusion survives the
selected refinement.  At the refined image of coarse data, the refined
event support lies below the refined owner region, the `algebraRefine`
image of the coarse record lies in the fine observer-local algebra at the
refined support meet the fine owner public algebra, and the recharted
coarse coordinate is the fine geometry readout of the mapped packet.  Each
conjunct is the corresponding fine common-origin conjunct at the mapped
pair rewritten along the naturality fields; the support conjunct also
follows from the coarse common origin through `region_refine_mono`
alone. -/
theorem refine_common_origin (R : SelectedEventTowerRefinement hrs Wc Wf)
    {i : sigc.Chart} {e : sigc.Event} (hi : Wc.activeChart i)
    (he : Wc.activeEvent e) (hv : Sc.visible i e) :
    N.regionLE s (N.regionRefine hrs (Wc.supportRead (Wc.packetOf e)))
        (N.regionRefine hrs (Ac.observerRegion (Wc.chartOwner i))) ∧
      T.algebraRefine hrs (Wc.recordRead (Wc.packetOf e)) ∈
        Af.observerLocalAlgebra (T.observerRefine hrs (Wc.chartOwner i))
            (N.regionRefine hrs (Wc.supportRead (Wc.packetOf e))) ⊓
          T.publicAlgebra s (T.observerRefine hrs (Wc.chartOwner i)) ∧
      R.rechart i (Sc.coordinate i e) =
        Wf.geometryRead (R.chartMap i) (R.packetMap (Wc.packetOf e)) := by
  obtain ⟨hfi, hfe, hfv⟩ := R.mapped_active_visible hi he hv
  obtain ⟨-, hmem, hcoord⟩ := Wf.common_origin hfi hfe hfv
  refine ⟨N.region_refine_mono hrs (Wc.common_origin hi he hv).1, ?_, ?_⟩
  · rw [R.packetOf_natural, R.supportRead_natural, R.chartOwner_natural,
      ← R.recordRead_natural] at hmem
    exact hmem
  · rw [R.packetOf_natural] at hcoord
    rw [← R.coordinate_rechart hi he hv]
    exact hcoord

end SelectedEventTowerRefinement

end Refinement

/-! ### The identity refinement witness -/

section RefinementWitness

/-- The identity refinement of the dimension-two witness world over the
single regulator: every slot map is the identity and the recharting is the
identity re-expression. -/
noncomputable def witnessIdentityRefinement :
    SelectedEventTowerRefinement (le_refl ()) witnessEventWorld
      witnessEventWorld where
  eventMap := id
  chartMap := id
  packetMap := id
  activeEvent_natural := fun _ h => h
  activeChart_natural := fun _ h => h
  visible_natural := fun h => h
  chartOwner_natural := fun _ => rfl
  observerRegion_natural := fun _ => rfl
  chartRegion_natural := fun _ => rfl
  packetOf_natural := fun _ => rfl
  supportRead_natural := fun _ => rfl
  recordRead_natural := fun _ => rfl
  rechart := fun _ x => x
  coordinate_rechart := fun _ _ _ => rfl

/-- The refinement interface is inhabited over the witness world. -/
theorem selectedEventTowerRefinement_inhabited :
    Nonempty (SelectedEventTowerRefinement (le_refl ()) witnessEventWorld
      witnessEventWorld) :=
  ⟨witnessIdentityRefinement⟩

end RefinementWitness

/-! ## Deliverable 11: raw event candidates and quotient descent -/

section RawEvents

/-- Raw event candidates: a carrier of raw event slots with a declared
semantic-coincidence setoid.  The setoid is supplied data; reflexivity,
symmetry, and transitivity are its `iseqv` obligations, and the negative
control below shows that the transitivity obligation carries content. -/
structure RawEventCandidates where
  /-- The raw candidate carrier. -/
  Raw : Type v
  /-- The declared semantic-coincidence setoid. -/
  coincidence : Setoid Raw

namespace RawEventCandidates

variable (C : RawEventCandidates)

/-- The event type of a raw candidate family: raw slots up to declared
coincidence. -/
def EventQuotient : Type v :=
  Quotient C.coincidence

/-- Descend a coincidence-invariant readout to the event quotient. -/
def descend {α : Sort*} (f : C.Raw → α)
    (hf : ∀ a b, C.coincidence.r a b → f a = f b) :
    C.EventQuotient → α :=
  Quotient.lift f hf

@[simp] theorem descend_mk {α : Sort*} (f : C.Raw → α)
    (hf : ∀ a b, C.coincidence.r a b → f a = f b) (a : C.Raw) :
    C.descend f hf (Quotient.mk C.coincidence a) = f a :=
  rfl

/-- Well-definedness of descended readouts: coincident raw candidates
receive equal values. -/
theorem descend_well_defined {α : Sort*} (f : C.Raw → α)
    (hf : ∀ a b, C.coincidence.r a b → f a = f b) {a b : C.Raw}
    (hab : C.coincidence.r a b) :
    C.descend f hf (Quotient.mk C.coincidence a) =
      C.descend f hf (Quotient.mk C.coincidence b) :=
  congrArg (C.descend f hf) (Quotient.sound hab)

end RawEventCandidates

/-- A coincidence-invariant raw packet assignment into one event-world
signature.  Invariance is the exact hypothesis under which the packet map
and every packet-derived readout descend to the event quotient. -/
structure RawPacketAssignment (C : RawEventCandidates)
    (sig : EventWorldSignature) where
  /-- The raw packet assignment. -/
  rawPacket : C.Raw → sig.Packet
  /-- Declared coincidence invariance of the raw packet assignment. -/
  packet_invariant : ∀ a b, C.coincidence.r a b → rawPacket a = rawPacket b

namespace RawPacketAssignment

variable {C : RawEventCandidates} {sig : EventWorldSignature}

/-- The descended packet map on the event quotient. -/
def inducedPacket (P : RawPacketAssignment C sig) :
    C.EventQuotient → sig.Packet :=
  C.descend P.rawPacket P.packet_invariant

@[simp] theorem inducedPacket_mk (P : RawPacketAssignment C sig)
    (a : C.Raw) :
    P.inducedPacket (Quotient.mk C.coincidence a) = P.rawPacket a :=
  rfl

/-- The induced packet map is well defined on coincidence classes. -/
theorem inducedPacket_well_defined (P : RawPacketAssignment C sig)
    {a b : C.Raw} (hab : C.coincidence.r a b) :
    P.inducedPacket (Quotient.mk C.coincidence a) =
      P.inducedPacket (Quotient.mk C.coincidence b) :=
  congrArg P.inducedPacket (Quotient.sound hab)

/-- Every packet field descends: any readout of the packet is
coincidence-invariant through the assignment, and its descended value at a
class is the readout of the induced packet. -/
theorem descend_packetField (P : RawPacketAssignment C sig) {α : Sort*}
    (g : sig.Packet → α) (q : C.EventQuotient) :
    C.descend (fun a => g (P.rawPacket a))
        (fun a b hab => congrArg g (P.packet_invariant a b hab)) q =
      g (P.inducedPacket q) := by
  refine Quotient.inductionOn q ?_
  intro a
  rfl

variable {r : ι} {A : ObserverAccessCut T N r}
variable {S : OPH.C2Soldering.EventFrameSoldering sig.Event sig.Chart}

/-- The event-support field descends through the induced packet map. -/
theorem descend_supportRead (P : RawPacketAssignment C sig)
    (W : FiniteEventGeometryWorld A sig S) (q : C.EventQuotient) :
    C.descend (fun a => W.supportRead (P.rawPacket a))
        (fun a b hab =>
          congrArg W.supportRead (P.packet_invariant a b hab)) q =
      W.supportRead (P.inducedPacket q) :=
  P.descend_packetField W.supportRead q

/-- The public-record field descends through the induced packet map. -/
theorem descend_recordRead (P : RawPacketAssignment C sig)
    (W : FiniteEventGeometryWorld A sig S) (q : C.EventQuotient) :
    C.descend (fun a => W.recordRead (P.rawPacket a))
        (fun a b hab =>
          congrArg W.recordRead (P.packet_invariant a b hab)) q =
      W.recordRead (P.inducedPacket q) :=
  P.descend_packetField W.recordRead q

/-- Each chart's geometry field descends through the induced packet map. -/
theorem descend_geometryRead (P : RawPacketAssignment C sig)
    (W : FiniteEventGeometryWorld A sig S) (i : sig.Chart)
    (q : C.EventQuotient) :
    C.descend (fun a => W.geometryRead i (P.rawPacket a))
        (fun a b hab =>
          congrArg (W.geometryRead i) (P.packet_invariant a b hab)) q =
      W.geometryRead i (P.inducedPacket q) :=
  P.descend_packetField (W.geometryRead i) q

end RawPacketAssignment

/-- Negative control, shared with the C2 germ stack: the three-germ overlap
relation `OPH.C2Soldering.overlapControl` is reflexive and symmetric and
fails transitivity, so no setoid has it as its coincidence relation and the
quotient constructor above is closed to it. -/
theorem overlapControl_no_coincidence_setoid :
    ¬ ∃ st : Setoid (Fin 3),
        ∀ a b, st.r a b ↔ OPH.C2Soldering.overlapControl a b := by
  rintro ⟨st, hst⟩
  obtain ⟨h01, h12, h02⟩ := OPH.C2Soldering.overlapControl_not_transitive
  exact h02 ((hst 0 2).mp
    (st.iseqv.trans ((hst 0 1).mpr h01) ((hst 1 2).mpr h12)))

/-- Witness raw candidates over the two-event witness world: two raw slots
per event, coincident exactly when their event components agree. -/
def witnessRawEventCandidates : RawEventCandidates where
  Raw := Fin 2 × Fin 2
  coincidence :=
    { r := fun a b => a.1 = b.1
      iseqv := ⟨fun _ => rfl, Eq.symm, Eq.trans⟩ }

/-- The witness raw packet assignment: the packet of a raw slot is the
witness packet of its event component. -/
def witnessRawPacketAssignment :
    RawPacketAssignment witnessRawEventCandidates witnessEventSignature where
  rawPacket := fun a => a.1
  packet_invariant := fun _ _ hab => hab

/-- The witness coincidence genuinely identifies distinct raw slots: the
two raw candidates of the first event are distinct, coincident, and receive
the same induced packet. -/
theorem witnessRaw_nontrivial_descent :
    ((0, 0) : Fin 2 × Fin 2) ≠ (0, 1) ∧
      witnessRawEventCandidates.coincidence.r (0, 0) (0, 1) ∧
      witnessRawPacketAssignment.inducedPacket
          (Quotient.mk witnessRawEventCandidates.coincidence (0, 0)) =
        witnessRawPacketAssignment.inducedPacket
          (Quotient.mk witnessRawEventCandidates.coincidence (0, 1)) := by
  refine ⟨by decide, rfl, ?_⟩
  exact witnessRawPacketAssignment.inducedPacket_well_defined rfl

end RawEvents

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.EventPrecedenceAdapter.precedes_asymm
#print axioms OPH.QFT.witnessPrecedenceAdapter_strict
#print axioms OPH.QFT.SelectedEventTowerRefinement.mapped_active_visible
#print axioms OPH.QFT.SelectedEventTowerRefinement.refine_common_origin
#print axioms OPH.QFT.selectedEventTowerRefinement_inhabited
#print axioms OPH.QFT.RawEventCandidates.descend_well_defined
#print axioms OPH.QFT.RawPacketAssignment.inducedPacket_well_defined
#print axioms OPH.QFT.RawPacketAssignment.descend_packetField
#print axioms OPH.QFT.RawPacketAssignment.descend_supportRead
#print axioms OPH.QFT.RawPacketAssignment.descend_recordRead
#print axioms OPH.QFT.RawPacketAssignment.descend_geometryRead
#print axioms OPH.QFT.overlapControl_no_coincidence_setoid
#print axioms OPH.QFT.witnessRaw_nontrivial_descent
