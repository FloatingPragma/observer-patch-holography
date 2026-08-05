import QFT.ObserverAccessCut
import Geometry.EventFrameSoldering

/-!
# Observer quantum-event worlds over the C2 soldering interface

This module implements the event-world layer of completion-plan issue
`#712` over the access cuts of `QFT/ObserverAccessCut.lean` and the C2
`EventFrameSoldering` interface.  An `EventWorldSignature` fixes finite
carrier types for events, charts, and packets, so every readout codomain
is independent of any quotient state.  A `FiniteEventGeometryWorld` binds
one packet per event to four declared readouts: the event's support
region, the chart support region with observer ownership, the public
record observable, and the observer-relative chart coordinate.

The common-origin theorem is the payoff: for an active visible pair the
event support lies below the owner's observer region, the packet-derived
record lies in the meet of the observer-local algebra at the event
support with the owner's public algebra, and the soldering coordinate is
a readout of the same packet.  One packet is the common source of the
support, the record, the accessibility, and the coordinate.

The witness section instantiates the world on the C2 two-event control
soldering over the dimension-two witness cut, with a packet-sensitive
scalar record readout, and evaluates the common-origin theorem on every
pair.  A second witness section instantiates the world over the
dimension-four proper-meet cut of `QFT/ObserverAccessCut.lean`: a
three-element packet carrier distinct from the event carrier, projector
records outside the scalar span, a diagonal owner region that the block
region fails to lie below, and a coordinate law factoring through a
packet-separating readout.

Claim boundary: every structure is a finite conditional interface over
the committed E1, C2, and access-cut packets.  Charts and observers stay
distinct types; ownership is declared data.  No physical region,
coverage, instrument, spacetime, causal cone, event actualization, or
continuum object is claimed, and no theorem grants any observer access
to the ambient private algebra.  E1 retains the noncommutative source
witness; source realization belongs to E2.
-/

namespace OPH.QFT

open OPH.Tower

universe u v

/-- Fixed finite carrier types for events, charts, and packets.  The
carriers are structural data, so readout codomains built from them do not
depend on any quotient state. -/
structure EventWorldSignature where
  Event : Type v
  Chart : Type v
  Packet : Type v
  eventFintype : Fintype Event
  chartFintype : Fintype Chart
  packetFintype : Fintype Packet

attribute [instance] EventWorldSignature.eventFintype
  EventWorldSignature.chartFintype EventWorldSignature.packetFintype

variable {ι : Type u} [Preorder ι]
variable {T : ConsensusTower ι} {N : FiniteCausalObserverNet T} {r : ι}

/-- A finite observer quantum-event world: one declared packet per event
determines the event support region, the public record observable, and
the chart coordinate, while declared chart ownership ties chart supports
to observer regions of one access cut.

`activeEvent` and `activeChart` are declared masks; every law is scoped
to active visible pairs.  Charts and observers remain distinct types
related only through `chartOwner`. -/
structure FiniteEventGeometryWorld (A : ObserverAccessCut T N r)
    (sig : EventWorldSignature)
    (S : OPH.C2Soldering.EventFrameSoldering sig.Event sig.Chart) where
  /-- Declared active-event mask. -/
  activeEvent : sig.Event → Prop
  /-- Declared active-chart mask. -/
  activeChart : sig.Chart → Prop
  /-- The declared packet of each event. -/
  packetOf : sig.Event → sig.Packet
  /-- The declared support region of each chart. -/
  chartRegion : sig.Chart → N.Region r
  /-- The declared owner of each chart.  Ownership is data; no claim of
  uniqueness, exclusivity, or physical attachment is included. -/
  chartOwner : sig.Chart → T.Observer r
  /-- Owner-support containment: an active chart's support region lies
  below its owner's observer region. -/
  chartRegion_le_owner : ∀ i, activeChart i →
    N.regionLE r (chartRegion i) (A.observerRegion (chartOwner i))
  /-- The packet-derived support region of an event. -/
  supportRead : sig.Packet → N.Region r
  /-- The packet-derived chart readout of the geometric coordinate. -/
  geometryRead : sig.Chart → sig.Packet → OPH.C1Lorentz.Herm2
  /-- Visibility support inclusion: an active visible event has its
  packet-derived support inside the seeing chart's support region. -/
  visible_support : ∀ {i e}, activeChart i → activeEvent e →
    S.visible i e →
      N.regionLE r (supportRead (packetOf e)) (chartRegion i)
  /-- The packet-derived coordinate law: on active visible pairs the C2
  soldering coordinate is a readout of the event's packet. -/
  coordinate_from_packet : ∀ {i e}, activeChart i → activeEvent e →
    S.visible i e →
      S.coordinate i e = geometryRead i (packetOf e)
  /-- The packet-derived public record observable. -/
  recordRead : sig.Packet → ConsensusTower.PrivateAlgebra T r
  /-- Record localization: the record of a packet lies in the local
  algebra of the packet's support region. -/
  record_local : ∀ p, recordRead p ∈ N.localAlgebra r (supportRead p)
  /-- Owner-public membership: on active visible pairs the record lies in
  the chart owner's public algebra. -/
  record_public : ∀ {i e}, activeChart i → activeEvent e →
    S.visible i e →
      recordRead (packetOf e) ∈ T.publicAlgebra r (chartOwner i)

namespace FiniteEventGeometryWorld

variable {A : ObserverAccessCut T N r} {sig : EventWorldSignature}
variable {S : OPH.C2Soldering.EventFrameSoldering sig.Event sig.Chart}

/-- The record of an active visible event is accessible to the chart
owner: owner-public membership composed with the cut's public
containment. -/
theorem record_mem_accessible (W : FiniteEventGeometryWorld A sig S)
    {i : sig.Chart} {e : sig.Event} (hi : W.activeChart i)
    (he : W.activeEvent e) (hv : S.visible i e) :
    W.recordRead (W.packetOf e) ∈ A.accessibleAlgebra (W.chartOwner i) :=
  A.public_le_accessible (W.chartOwner i) (W.record_public hi he hv)

/-- Common-origin theorem: on an active visible pair, one packet
determines everything.  The event support lies below the owner's observer
region, the record lies in the meet of the observer-local algebra at the
event support with the owner's public algebra, and the chart coordinate
is a readout of the same packet.  The proof composes the declared readout
laws `visible_support`, `chartRegion_le_owner`, `record_local`,
`record_public`, and `coordinate_from_packet`; no content beyond those
declared laws is derived. -/
theorem common_origin (W : FiniteEventGeometryWorld A sig S)
    {i : sig.Chart} {e : sig.Event} (hi : W.activeChart i)
    (he : W.activeEvent e) (hv : S.visible i e) :
    N.regionLE r (W.supportRead (W.packetOf e))
        (A.observerRegion (W.chartOwner i)) ∧
      W.recordRead (W.packetOf e) ∈
        A.observerLocalAlgebra (W.chartOwner i)
            (W.supportRead (W.packetOf e)) ⊓
          T.publicAlgebra r (W.chartOwner i) ∧
      S.coordinate i e = W.geometryRead i (W.packetOf e) := by
  have hsup := W.visible_support hi he hv
  have hle := N.regionLE_trans r hsup (W.chartRegion_le_owner i hi)
  have hpub := W.record_public hi he hv
  refine ⟨hle, StarSubalgebra.mem_inf.mpr ⟨?_, hpub⟩,
    W.coordinate_from_packet hi he hv⟩
  exact A.mem_observerLocalAlgebra.mpr
    ⟨W.record_local _, A.public_le_accessible (W.chartOwner i) hpub⟩

end FiniteEventGeometryWorld

/-! ## A witness world on the C2 control soldering

The witness reuses the dimension-two access cut and the C2 two-event
control soldering.  Both events carry distinct packets, the record
readout distinguishes them, and every pair is active and visible, so the
common-origin theorem applies everywhere. -/

section EventWorldWitness

open EventAlgebra

/-- The witness signature: two events, two charts, two packets.  The
declaration is reducible so that numeral and coercion instances on the
carriers resolve. -/
abbrev witnessEventSignature : EventWorldSignature :=
  { Event := Fin 2
    Chart := Fin 2
    Packet := Fin 2
    eventFintype := inferInstance
    chartFintype := inferInstance
    packetFintype := inferInstance }

/-- The witness event world: identity packets on the C2 control
soldering, one owned full-region chart support, a proper singleton event
support, and a packet-sensitive scalar record readout. -/
noncomputable def witnessEventWorld :
    FiniteEventGeometryWorld witnessAccessCut witnessEventSignature
      OPH.C2Soldering.controlSoldering where
  activeEvent := fun _ => True
  activeChart := fun _ => True
  packetOf := id
  chartRegion := fun _ => ({0, 1} : Finset (Fin 2))
  chartOwner := fun _ => ()
  chartRegion_le_owner := fun _ _ => Finset.Subset.rfl
  supportRead := fun _ => ({0} : Finset (Fin 2))
  geometryRead := fun i p => OPH.C2Soldering.controlSoldering.coordinate i p
  visible_support := by
    intro i e _ _ _
    exact Finset.singleton_subset_iff.mpr (by simp)
  coordinate_from_packet := fun _ _ _ => rfl
  recordRead := fun p => ((p : ℕ) + 1 : ℂ) • 1
  record_local := fun p =>
    ((diagonalPartition 2).publicSubalgebra).smul_mem (one_mem _) _
  record_public := fun _ _ _ => (oneBlockPartition_mem_iff _).mpr ⟨_, rfl⟩

/-- The witness record readout distinguishes the two packets, so the
record is packet data rather than a constant adaptor. -/
theorem witnessEventWorld_record_ne :
    witnessEventWorld.recordRead 0 ≠ witnessEventWorld.recordRead 1 := by
  intro h
  have h2 : ((((0 : Fin 2) : ℕ) + 1 : ℂ) • (1 : Matrix (Fin 2) (Fin 2) ℂ)) =
      ((((1 : Fin 2) : ℕ) + 1 : ℂ) • (1 : Matrix (Fin 2) (Fin 2) ℂ)) := h
  have h00 := congrFun (congrFun h2 0) 0
  norm_num [Matrix.smul_apply, Matrix.one_apply_eq] at h00

/-- The common-origin theorem evaluated on every witness pair: each pair
is active and visible, so support containment, record membership, and the
packet coordinate law hold across the whole witness world. -/
theorem witnessEventWorld_common_origin (i e : Fin 2) :
    witnessNet.regionLE ()
        (witnessEventWorld.supportRead (witnessEventWorld.packetOf e))
        (witnessAccessCut.observerRegion (witnessEventWorld.chartOwner i)) ∧
      witnessEventWorld.recordRead (witnessEventWorld.packetOf e) ∈
        witnessAccessCut.observerLocalAlgebra
            (witnessEventWorld.chartOwner i)
            (witnessEventWorld.supportRead (witnessEventWorld.packetOf e)) ⊓
          witnessTower.publicAlgebra () (witnessEventWorld.chartOwner i) ∧
      OPH.C2Soldering.controlSoldering.coordinate i e =
        witnessEventWorld.geometryRead i (witnessEventWorld.packetOf e) :=
  witnessEventWorld.common_origin trivial trivial trivial

end EventWorldWitness

/-! ## A dimension-four world with nondegenerate origin receipts

The witness above discharges every common-origin conjunct through the
top owner region, scalar records, and an identity packet assignment.
The world in this section sits over the dimension-four proper-meet stack
of `QFT/ObserverAccessCut.lean` and separates each conjunct from those
degenerate discharges: the tower public algebra is the full diagonal
span, so the record readout takes nonscalar public values; the owner
region is the diagonal region, which the block region fails to lie
below, so the support containment is falsifiable; the packet carrier is
a three-element type distinct from the two-element event carrier and the
packet assignment misses one packet label; the coordinate law factors
through a packet readout that separates packet labels. -/

section ProperMeetEventWorld

open EventAlgebra

/-- The dimension-four witness signature: two events, two charts, three
packets, so the packet carrier is a type distinct from the event
carrier. -/
abbrev properMeetEventSignature : EventWorldSignature :=
  { Event := Fin 2
    Chart := Fin 2
    Packet := Fin 3
    eventFintype := inferInstance
    chartFintype := inferInstance
    packetFintype := inferInstance }

/-- The event label read back from a packet label. -/
def properMeetPacketEvent (p : Fin 3) : Fin 2 :=
  if p = 0 then 0 else 1

/-- The dimension-four witness world: an embedding packet assignment
into the three-element packet carrier, the diagonal owner region,
projector records in the diagonal public algebra, and a coordinate
readout through `properMeetPacketEvent`. -/
noncomputable def properMeetEventWorld :
    FiniteEventGeometryWorld properMeetAccessCut properMeetEventSignature
      OPH.C2Soldering.controlSoldering where
  activeEvent := fun _ => True
  activeChart := fun _ => True
  packetOf := Fin.castSucc
  chartRegion := fun _ => ProperMeetRegion.diag
  chartOwner := fun _ => ()
  chartRegion_le_owner := fun _ _ => Or.inl rfl
  supportRead := fun _ => ProperMeetRegion.diag
  geometryRead := fun i p =>
    OPH.C2Soldering.controlSoldering.coordinate i (properMeetPacketEvent p)
  visible_support := fun _ _ _ => Or.inl rfl
  coordinate_from_packet := by
    intro i e _ _ _
    fin_cases e <;> rfl
  recordRead := fun p => (diagonalPartition 4).proj p.castSucc
  record_local := fun p => (diagonalPartition 4).proj_mem_span p.castSucc
  record_public := fun {i e} _ _ _ =>
    (diagonalPartition 4).proj_mem_span (Fin.castSucc e).castSucc

/-- The record readout of the first event, unfolded to the first
diagonal projector. -/
theorem properMeetEventWorld_record_eq :
    properMeetEventWorld.recordRead (properMeetEventWorld.packetOf 0) =
      (diagonalPartition 4).proj 0 := rfl

/-- Nondegeneracy receipt for the record conjunct: the record of the
first event lies outside the scalar span, so the record readout takes
nonscalar public values. -/
theorem properMeetEventWorld_record_nonscalar :
    properMeetEventWorld.recordRead (properMeetEventWorld.packetOf 0) ∉
      (oneBlockPartition 4).publicSubalgebra := by
  rw [properMeetEventWorld_record_eq]
  intro hmem
  obtain ⟨c, hc⟩ := (oneBlockPartition_mem_iff _).mp hmem
  have h00 := congrFun (congrFun hc 0) 0
  have h11 := congrFun (congrFun hc 1) 1
  simp [diagonalPartition, Matrix.smul_apply] at h00 h11
  rw [h00] at h11
  exact one_ne_zero h11

/-- The record readout separates the two events. -/
theorem properMeetEventWorld_record_ne :
    properMeetEventWorld.recordRead (properMeetEventWorld.packetOf 0) ≠
      properMeetEventWorld.recordRead (properMeetEventWorld.packetOf 1) := by
  intro h
  have h' : (diagonalPartition 4).proj 0 = (diagonalPartition 4).proj 1 := h
  have h00 := congrFun (congrFun h' 0) 0
  simp [diagonalPartition] at h00

/-- Nondegeneracy receipt for the support conjunct: the block region
fails the containment below the owner region, so the containment law is
falsifiable in this world. -/
theorem properMeetEventWorld_support_falsifiable :
    ¬ properMeetNet.regionLE () ProperMeetRegion.block
      (properMeetAccessCut.observerRegion
        (properMeetEventWorld.chartOwner 0)) :=
  fun h => ProperMeetRegion.block_not_le_diag h

/-- Nondegeneracy receipt for the coordinate conjunct: the packet
readout separates packet labels, so the coordinate law composes through
a readout with content. -/
theorem properMeetEventWorld_geometryRead_ne :
    properMeetEventWorld.geometryRead 0 0 ≠
      properMeetEventWorld.geometryRead 0 1 := by
  intro h
  have ht := congrArg Prod.fst h
  norm_num [properMeetEventWorld, properMeetPacketEvent,
    OPH.C2Soldering.controlSoldering, OPH.C2Soldering.controlEventPoint,
    OPH.C2Soldering.controlChartOffset] at ht

/-- The packet assignment misses the third packet label, so the packet
carrier is strictly larger than the event image and the assignment is no
identity map. -/
theorem properMeetEventWorld_packet_missed (e : Fin 2) :
    properMeetEventWorld.packetOf e ≠ 2 := by
  fin_cases e <;> decide

/-- The common-origin theorem evaluated on every dimension-four witness
pair. -/
theorem properMeetEventWorld_common_origin (i e : Fin 2) :
    properMeetNet.regionLE ()
        (properMeetEventWorld.supportRead (properMeetEventWorld.packetOf e))
        (properMeetAccessCut.observerRegion
          (properMeetEventWorld.chartOwner i)) ∧
      properMeetEventWorld.recordRead (properMeetEventWorld.packetOf e) ∈
        properMeetAccessCut.observerLocalAlgebra
            (properMeetEventWorld.chartOwner i)
            (properMeetEventWorld.supportRead
              (properMeetEventWorld.packetOf e)) ⊓
          properMeetTower.publicAlgebra ()
            (properMeetEventWorld.chartOwner i) ∧
      OPH.C2Soldering.controlSoldering.coordinate i e =
        properMeetEventWorld.geometryRead i
          (properMeetEventWorld.packetOf e) :=
  properMeetEventWorld.common_origin trivial trivial trivial

end ProperMeetEventWorld

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.FiniteEventGeometryWorld.record_mem_accessible
#print axioms OPH.QFT.FiniteEventGeometryWorld.common_origin
#print axioms OPH.QFT.witnessEventWorld_record_ne
#print axioms OPH.QFT.witnessEventWorld_common_origin
#print axioms OPH.QFT.properMeetEventWorld_record_nonscalar
#print axioms OPH.QFT.properMeetEventWorld_record_ne
#print axioms OPH.QFT.properMeetEventWorld_support_falsifiable
#print axioms OPH.QFT.properMeetEventWorld_geometryRead_ne
#print axioms OPH.QFT.properMeetEventWorld_packet_missed
#print axioms OPH.QFT.properMeetEventWorld_common_origin
