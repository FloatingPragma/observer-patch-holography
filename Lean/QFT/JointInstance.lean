import QFT.SourceLinkSquare
import Tower.OperationalObserver

/-!
# The first finite joint instance of the E2 integration gate

This module implements the first milestone of completion-plan issue
`#693`.  A `JointFiniteInstance` packages, over one consensus tower, a
finite causal observer net, an observer access cut, a finite observer
event world over the C2 soldering interface, an operational observer
whose accessible interface agrees with the access cut at the observer's
own label, and a source-link square against a bare consensus tower.  The
compatibility field `accessible_compatible` is the sharpest well-typed
link between the B17 observer receipt and the E6 access cut: the two
accessible interfaces are one subalgebra of the private matrix algebra
at the joint regulator.

`joint_common_origin` restates the E6 common-origin theorem through the
operational observer's own interface: for an active visible pair whose
chart owner is the observer's label, the event support lies below the
declared observer region, the packet record lands in the observer's
accessible algebra through the compatibility field, and the soldering
coordinate is a readout of the same packet.  `joint_record_readback`
adds the B17 clause on the same data: the observer's self-readback fixes
the record of every such pair, since the record is public to the owner.

The witness instance glues the committed witnesses.  The consensus side
reuses the E6 dimension-two stack of `QFT/ObserverAccessCut.lean` and
`QFT/ObserverEventWorld.lean`; the bare side reuses the demo tower, the
demo branch, and the demo source-link square of
`Tower/EventGeometryReadout.lean` and `QFT/SourceLinkSquare.lean`; the
two sides keep unrelated regulator index types, `Unit` and `ℕ`.  The
operational observer `jointWitnessObserver` is a fresh dimension-two
inhabitant over the E6 witness tower: its accessible interface is the
diagonal algebra of the access cut, its readback is the trace pairing
against the selected witness state followed by the scalar embedding,
its control is the identity intervention, and its records ride the
one-block partition of the tower.  The B17 witness observer lives over
a different tower constant and is reused nowhere here.

Separation facts keep the instance away from constant adaptors: the
observer's accessible interface sits strictly between the tower public
scalars and the private matrix algebra, the readback annihilates the
raising matrix unit while fixing every public scalar, the square's event
encoding distinguishes the selected world packets from a constant packet
map, the selected observer readout separates two bare quotient states,
and the world record readout separates the two packets.

## Claim boundary

The instance certifies finite joint inhabitation only: one explicit
tower carries all listed interfaces at once.  Source production of the
tower data, CP or CPTP structure for the declared maps, adaptive
measurement scheduling, C2 realization, and every continuum statement
remain open E2 scope.  No uniqueness, canonicity, or physical
attachment is claimed for any component.
-/

namespace OPH.QFT

open OPH.Tower (ConsensusTower OperationalObserver)
open OPH.EinsteinBranch (BareConsensusTower demoTower)

universe u u' v

/-- A finite joint instance of the E2 integration gate: one consensus
tower carrying a finite causal observer net, an observer access cut, a
finite observer event world, and an operational observer, together with
a source-link square against a bare consensus tower.

`accessible_compatible` is the joint glue: the operational observer's
accessible interface at the joint regulator is the access cut's
accessible algebra at the observer's own label.  Every component keeps
its own claim boundary; the instance adds no field beyond the
compatibility equation. -/
structure JointFiniteInstance (κ : Type u) [Preorder κ]
    (ι : Type u') [Preorder ι] where
  /-- The consensus tower carrying the joint structure. -/
  tower : ConsensusTower κ
  /-- The finite causal observer net over the tower. -/
  net : FiniteCausalObserverNet tower
  /-- The joint regulator at which the cut and the world live. -/
  regulator : κ
  /-- The observer access cut on the net at the joint regulator. -/
  cut : ObserverAccessCut tower net regulator
  /-- The event-world signature: fixed finite carriers. -/
  signature : EventWorldSignature.{v}
  /-- The C2 soldering interface on the signature carriers. -/
  soldering :
    OPH.C2Soldering.EventFrameSoldering signature.Event signature.Chart
  /-- The finite observer event world over the cut and the soldering. -/
  world : FiniteEventGeometryWorld cut signature soldering
  /-- The operational observer over the same tower. -/
  observer : OperationalObserver tower
  /-- Compatibility of the observer with the cut: the accessible
  interface at the joint regulator is the cut's accessible algebra at
  the observer's own label. -/
  accessible_compatible :
    observer.accessible regulator =
      cut.accessibleAlgebra (observer.label regulator)
  /-- The bare consensus tower of the source-link square.  Its regulator
  index type is unrelated to the consensus index type. -/
  bareTower : BareConsensusTower ι
  /-- The source-link square gluing the bare readouts to the selected
  event world. -/
  square : SourceLinkSquare bareTower world

namespace JointFiniteInstance

variable {κ : Type u} [Preorder κ] {ι : Type u'} [Preorder ι]

/-- The joint common-origin theorem: the E6 common-origin conclusion
restated through the operational observer's own interface.  For an
active visible pair whose chart owner is the observer's label, the event
support lies below the declared observer region, the packet record lies
in the observer's accessible interface by the compatibility field, and
the soldering coordinate is a readout of the same packet. -/
theorem joint_common_origin (J : JointFiniteInstance κ ι)
    {i : J.signature.Chart} {e : J.signature.Event}
    (hi : J.world.activeChart i) (he : J.world.activeEvent e)
    (hv : J.soldering.visible i e)
    (ho : J.world.chartOwner i = J.observer.label J.regulator) :
    J.net.regionLE J.regulator
        (J.world.supportRead (J.world.packetOf e))
        (J.cut.observerRegion (J.world.chartOwner i)) ∧
      J.world.recordRead (J.world.packetOf e) ∈
        J.observer.accessible J.regulator ∧
      J.soldering.coordinate i e =
        J.world.geometryRead i (J.world.packetOf e) := by
  obtain ⟨hle, hmem, hcoord⟩ := J.world.common_origin hi he hv
  refine ⟨hle, ?_, hcoord⟩
  have hacc : J.world.recordRead (J.world.packetOf e) ∈
      J.cut.accessibleAlgebra (J.world.chartOwner i) :=
    J.cut.observerLocal_le_accessible (J.world.chartOwner i)
      (J.world.supportRead (J.world.packetOf e))
      (StarSubalgebra.mem_inf.mp hmem).1
  rw [ho] at hacc
  rw [J.accessible_compatible]
  exact hacc

/-- The joint readback law: the operational observer's self-readback
fixes the record of every active visible pair whose chart owner is the
observer's label.  The record is public to the owner by the world's
owner-public law, and the B17 readback clause fixes public elements. -/
theorem joint_record_readback (J : JointFiniteInstance κ ι)
    {i : J.signature.Chart} {e : J.signature.Event}
    (hi : J.world.activeChart i) (he : J.world.activeEvent e)
    (hv : J.soldering.visible i e)
    (ho : J.world.chartOwner i = J.observer.label J.regulator) :
    J.observer.readback J.regulator
        (J.world.recordRead (J.world.packetOf e)) =
      J.world.recordRead (J.world.packetOf e) := by
  have hpub := J.world.record_public hi he hv
  rw [ho] at hpub
  exact J.observer.readback_fixes_public J.regulator _ hpub

end JointFiniteInstance

/-! ## A fresh operational observer over the E6 witness tower

The E6 witness tower is the constant one-regulator tower of the
one-block partition at dimension two, so its public algebra is the
scalar span of the identity and its records ride the one-block active
index.  The observer declares the diagonal algebra of the witness access
cut as its accessible interface, reads back through the trace pairing
against the selected witness state, controls through the identity
intervention, and predicts through the identity on the discrete record
order. -/

section JointWitness

/-- The fresh operational observer over the E6 witness tower.  The
accessible interface is the diagonal algebra of `witnessAccessCut`, the
readback is the state pairing `X ↦ Tr(ρ X) • 1` against the selected
witness state, and every remaining clause is discharged on the constant
tower data. -/
noncomputable def jointWitnessObserver : OperationalObserver witnessTower where
  label := fun _ => ()
  label_refine := by intros; rfl
  accessible := fun _ => (diagonalPartition 2).publicSubalgebra
  public_le_accessible := fun _ =>
    oneBlockPartition_publicSubalgebra_le _
  readback := fun r X => (witnessTower.state r () * X).trace • 1
  readback_mem_public := fun _ X _ =>
    ((oneBlockPartition 2).publicSubalgebra).smul_mem (one_mem _) _
  readback_fixes_public := fun r X hX => by
    obtain ⟨c, rfl⟩ := (oneBlockPartition_mem_iff X).mp hX
    rw [mul_smul_comm, mul_one, Matrix.trace_smul,
      (witnessTower.state_isState r ()).2]
    simp
  record_durable := fun _ _ => rfl
  control := fun _ _ => LinearMap.id
  control_preserves_accessible := fun _ _ _ hX => hX
  predict := fun _ x => x
  predict_no_regress := fun _ _ h => h
  predict_refine := by intros; rfl
  readout := fun r x =>
    (witnessTower.state r () * (witnessTower.recordElement r () x).1).trace
  predict_readout_congruent := fun _ _ _ h => h
  readout_eq_pairing := fun _ _ => rfl
  record_fixed_by_control := fun _ _ _ => rfl

/-- The witness readback annihilates the raising matrix unit: the
selected witness state pairs to zero against every off-diagonal matrix
unit. -/
theorem jointWitnessObserver_readback_raise :
    jointWitnessObserver.readback () FiniteCausalObserverNet.m2Raise = 0 := by
  show ((Matrix.diagonal (Pi.single 0 1) : Matrix (Fin 2) (Fin 2) ℂ) *
      FiniteCausalObserverNet.m2Raise).trace •
      (1 : Matrix (Fin 2) (Fin 2) ℂ) = 0
  have htr : ((Matrix.diagonal (Pi.single 0 1) : Matrix (Fin 2) (Fin 2) ℂ) *
      FiniteCausalObserverNet.m2Raise).trace = 0 := by
    rw [Matrix.trace_fin_two]
    simp [FiniteCausalObserverNet.m2Raise, Matrix.mul_apply]
  rw [htr, zero_smul]

/-- The witness readback differs from the identity map: it annihilates
the raising matrix unit, which is nonzero. -/
theorem jointWitnessObserver_readback_ne_id :
    jointWitnessObserver.readback () ≠ id := by
  intro h
  have hr := congrFun h FiniteCausalObserverNet.m2Raise
  rw [jointWitnessObserver_readback_raise] at hr
  have h01 : (0 : ℂ) = 1 := congrFun (congrFun hr (0 : Fin 2)) (1 : Fin 2)
  exact zero_ne_one h01

/-- Strictness through the observer, lower interposition: the tower
public scalars sit strictly inside the observer's accessible
interface. -/
theorem jointWitnessObserver_public_lt_accessible :
    witnessTower.publicAlgebra () (jointWitnessObserver.label ()) <
      jointWitnessObserver.accessible () :=
  witnessAccessCut_public_lt_accessible ()

/-- Strictness through the observer, upper interposition: the observer's
accessible interface is a proper subalgebra of the private matrix
algebra. -/
theorem jointWitnessObserver_accessible_lt_top :
    jointWitnessObserver.accessible () <
      (⊤ : StarSubalgebra ℂ
        (ConsensusTower.PrivateAlgebra witnessTower ())) :=
  witnessAccessCut_accessible_lt_top ()

/-- The witness joint instance: the E6 dimension-two stack, the fresh
operational observer, and the demo source-link square in one record.
The compatibility field holds definitionally, since the observer's
accessible interface is declared as the cut's diagonal algebra. -/
noncomputable def jointWitnessInstance : JointFiniteInstance Unit ℕ where
  tower := witnessTower
  net := witnessNet
  regulator := ()
  cut := witnessAccessCut
  signature := witnessEventSignature
  soldering := OPH.C2Soldering.controlSoldering
  world := witnessEventWorld
  observer := jointWitnessObserver
  accessible_compatible := rfl
  bareTower := demoTower
  square := demoSourceLinkSquare

/-- The joint common-origin theorem evaluated on every witness pair:
each pair is active and visible and the sole observer owns every chart,
so support containment, accessible-record membership through the
observer's interface, and the packet coordinate law hold across the
whole witness world. -/
theorem jointWitnessInstance_common_origin (i e : Fin 2) :
    jointWitnessInstance.net.regionLE jointWitnessInstance.regulator
        (jointWitnessInstance.world.supportRead
          (jointWitnessInstance.world.packetOf e))
        (jointWitnessInstance.cut.observerRegion
          (jointWitnessInstance.world.chartOwner i)) ∧
      jointWitnessInstance.world.recordRead
          (jointWitnessInstance.world.packetOf e) ∈
        jointWitnessInstance.observer.accessible
          jointWitnessInstance.regulator ∧
      jointWitnessInstance.soldering.coordinate i e =
        jointWitnessInstance.world.geometryRead i
          (jointWitnessInstance.world.packetOf e) :=
  jointWitnessInstance.joint_common_origin (i := i) (e := e)
    trivial trivial trivial rfl

/-- The joint readback law evaluated on every witness pair: the fresh
observer's readback fixes the record of each packet. -/
theorem jointWitnessInstance_record_readback (i e : Fin 2) :
    jointWitnessInstance.observer.readback jointWitnessInstance.regulator
        (jointWitnessInstance.world.recordRead
          (jointWitnessInstance.world.packetOf e)) =
      jointWitnessInstance.world.recordRead
        (jointWitnessInstance.world.packetOf e) :=
  jointWitnessInstance.joint_record_readback (i := i) (e := e)
    trivial trivial trivial rfl

/-- Separation through the instance, square side: the event encoding of
the packaged square distinguishes the selected world packets from the
constant packet map. -/
theorem jointWitnessInstance_eventEncode_separates :
    jointWitnessInstance.square.eventEncode 0
        jointWitnessInstance.world.packetOf ≠
      jointWitnessInstance.square.eventEncode 0 (fun _ => (0 : Fin 2)) :=
  demoSourceLinkSquare_eventEncode_separates

/-- Separation through the instance, readout side: the packaged selected
observer readout separates two bare quotient states. -/
theorem jointWitnessInstance_readout_separates :
    jointWitnessInstance.square.readout.read 0 (true, true) ≠
      jointWitnessInstance.square.readout.read 0 (false, false) :=
  demoObserverReadout_separates

/-- Separation through the instance, world side: the record readout of
the packaged world separates the two packets. -/
theorem jointWitnessInstance_record_separates :
    jointWitnessInstance.world.recordRead (0 : Fin 2) ≠
      jointWitnessInstance.world.recordRead (1 : Fin 2) :=
  witnessEventWorld_record_ne

/-- Nontriviality receipt for the witness joint instance: strict
interposition of the observer's accessible interface on both sides, a
readback distinct from the identity, and the three packaged separation
facts. -/
theorem jointWitnessInstance_nontrivial :
    (witnessTower.publicAlgebra () (jointWitnessObserver.label ()) <
        jointWitnessObserver.accessible ()) ∧
      (jointWitnessObserver.accessible () <
        (⊤ : StarSubalgebra ℂ
          (ConsensusTower.PrivateAlgebra witnessTower ()))) ∧
      jointWitnessObserver.readback () ≠ id ∧
      jointWitnessInstance.square.eventEncode 0
          jointWitnessInstance.world.packetOf ≠
        jointWitnessInstance.square.eventEncode 0 (fun _ => (0 : Fin 2)) ∧
      jointWitnessInstance.world.recordRead (0 : Fin 2) ≠
        jointWitnessInstance.world.recordRead (1 : Fin 2) :=
  ⟨jointWitnessObserver_public_lt_accessible,
    jointWitnessObserver_accessible_lt_top,
    jointWitnessObserver_readback_ne_id,
    jointWitnessInstance_eventEncode_separates,
    jointWitnessInstance_record_separates⟩

end JointWitness

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.JointFiniteInstance.joint_common_origin
#print axioms OPH.QFT.JointFiniteInstance.joint_record_readback
#print axioms OPH.QFT.jointWitnessObserver_readback_raise
#print axioms OPH.QFT.jointWitnessObserver_readback_ne_id
#print axioms OPH.QFT.jointWitnessObserver_public_lt_accessible
#print axioms OPH.QFT.jointWitnessObserver_accessible_lt_top
#print axioms OPH.QFT.jointWitnessInstance_common_origin
#print axioms OPH.QFT.jointWitnessInstance_record_readback
#print axioms OPH.QFT.jointWitnessInstance_eventEncode_separates
#print axioms OPH.QFT.jointWitnessInstance_readout_separates
#print axioms OPH.QFT.jointWitnessInstance_record_separates
#print axioms OPH.QFT.jointWitnessInstance_nontrivial
