import ObserverPatchHolography.EinsteinBranch.Consensus

/-!
# All-quotient event and geometry readout over the bare consensus tower

This module implements deliverables 1, 3, 5, and 6 of completion-plan issue
`#713` over the Einstein-branch `BareConsensusTower`.

`SelectedSettledBranch` packages one physical quotient point per regulator,
fixed by the normal form and coherent under the physical fine-to-coarse
maps.  Existence of such a branch for a given tower is an explicit premise,
carried by inhabitation of the structure.  `demoSelectedBranch` inhabits it
on the concrete `demoTower`, and `demoUnsettledCandidate` is a point family
on the same tower that fails the settledness field on explicit data.

`EventGeometryReadoutFragment` carries event and geometry readouts defined
on every physical quotient state, with fixed packet codomains per
regulator, normal-form invariance, functorial fine-to-coarse packet maps,
and coarse naturality.  `demoFragment` inhabits the fragment on `demoTower`
with readouts that separate two consistent quotient states.

`ReadoutFragmentPacket` is the Einstein handoff: it bundles a fragment and
nothing else; the obligations that remain with the Einstein readout
completion are recorded in its docstring and are not fields.

`boundaryFiber_of_complete_signature` derives the boundary-fibre conclusion
from an event signature that factors through the protected boundary and is
complete on consistent states.  Instantiating the signature as the
protected boundary itself is inadmissible, and
`boundaryFiber_of_complete_signature_trivial_instance` records that
degeneration as a definitional equivalence.

## Claim boundary

Every statement is a finite conditional structure over declared tower data.
The fragment supplies the event and geometry portion of the Einstein
readout interface only.  Modular, stress, entropy, and scale readouts and
every Einstein premise (first law, MaxEnt envelope, small-ball, Ward,
Bianchi, vacuum reference, universal coupling, scale) remain with issue
`#694`.  The matching square between these readouts and the E6 selected
observer event world of issue `#712` (deliverables 2 and 4 of issue `#713`)
is a separate composition and is open.  No physical spacetime, no
probability law, and no continuum claim follow.
-/

namespace OPH.Tower

open OPH.EinsteinBranch

universe u v

/-! ## Selected settled branch -/

/-- A selected settled branch of a bare consensus tower: one physical
quotient point per regulator, fixed by the normal form and coherent under
the physical fine-to-coarse maps.

Existence of a coherent settled branch is an explicit premise of the E7
source-link square; nothing in the bare tower interface forces one to
exist.  The protected boundary carries no coarse-naturality law in
`BareConsensusTower`, so this structure states no cross-regulator boundary
coherence; the provable boundary facts are
`SelectedSettledBranch.boundary_point_eq_of_settles_to` per regulator and
`SelectedSettledBranch.boundary_coarse_point` for coarsened branch
points. -/
structure SelectedSettledBranch {ι : Type u} [Preorder ι]
    (B : BareConsensusTower ι) where
  /-- The branch point at each regulator. -/
  point : ∀ r, B.Quot r
  /-- Each branch point is a fixed point of the normal form. -/
  settled : ∀ r, B.normalForm r (point r) = point r
  /-- The branch is coherent under the physical fine-to-coarse maps. -/
  coherent : ∀ {r s} (hrs : r ≤ s), B.coarse hrs (point s) = point r

namespace SelectedSettledBranch

variable {ι : Type u} [Preorder ι] {B : BareConsensusTower ι}

/-- Settledness forces mismatch zero: the branch point equals its own
normal form, and normal forms are consistent. -/
theorem mismatch_point_eq_zero (br : SelectedSettledBranch B) (r : ι) :
    B.mismatch r (br.point r) = 0 := by
  have h : B.consistent r (B.normalForm r (br.point r)) :=
    B.normalForm_consistent r (br.point r)
  rw [br.settled r] at h
  exact (B.consistent_iff r (br.point r)).1 h

/-- The branch point is consistent at every regulator. -/
theorem consistent_point (br : SelectedSettledBranch B) (r : ι) :
    B.consistent r (br.point r) :=
  (B.consistent_iff r (br.point r)).2 (br.mismatch_point_eq_zero r)

/-- Per-regulator boundary fact: the branch point carries the same
protected boundary value as every quotient state that settles to it. -/
theorem boundary_point_eq_of_settles_to (br : SelectedSettledBranch B)
    (r : ι) {q : B.Quot r} (h : B.normalForm r q = br.point r) :
    B.boundary r (br.point r) = B.boundary r q := by
  rw [← h]
  exact B.boundary_normalForm r q

/-- The provable cross-regulator boundary fact: coherence transports the
branch point, so the regulator-`r` boundary of the coarsened fine point is
the boundary of the coarse point.  A comparison between `B.Boundary s` and
`B.Boundary r` values is not well typed, and the bare tower declares no
boundary coarse-naturality law. -/
theorem boundary_coarse_point (br : SelectedSettledBranch B) {r s : ι}
    (hrs : r ≤ s) :
    B.boundary r (B.coarse hrs (br.point s)) = B.boundary r (br.point r) :=
  congrArg (B.boundary r) (br.coherent hrs)

/-- Settledness descends along every physical coarse map by normal-form
naturality alone; the coherence field is not used. -/
theorem settled_coarse_point (br : SelectedSettledBranch B) {r s : ι}
    (hrs : r ≤ s) :
    B.normalForm r (B.coarse hrs (br.point s)) =
      B.coarse hrs (br.point s) := by
  rw [← B.normalForm_natural hrs, br.settled s]

end SelectedSettledBranch

/-- The constant consistent family is a selected settled branch of the demo
tower: its points are normal-form fixed and the demo coarse maps are the
identity. -/
def demoSelectedBranch : SelectedSettledBranch demoTower where
  point := fun _ => (true, true)
  settled := fun _ => rfl
  coherent := fun _ => rfl

/-- Negative-control candidate: the constant mixed family on the demo
tower. -/
def demoUnsettledCandidate (r : ℕ) : demoTower.Quot r := (true, false)

/-- The mixed candidate fails settledness on explicit data: at regulator
`0` its normal form is `(true, true)` and differs in the second
component. -/
theorem demoUnsettledCandidate_not_settled :
    ¬ ∀ r, demoTower.normalForm r (demoUnsettledCandidate r) =
      demoUnsettledCandidate r := by
  intro h
  have h0 : ((true, true) : Bool × Bool) = (true, false) := h 0
  exact Bool.noConfusion (congrArg Prod.snd h0)

/-- No selected settled branch has the mixed candidate as its point
family. -/
theorem demoUnsettledCandidate_no_branch :
    ¬ ∃ br : SelectedSettledBranch demoTower,
      br.point = demoUnsettledCandidate := by
  rintro ⟨br, hpoint⟩
  refine demoUnsettledCandidate_not_settled fun r => ?_
  rw [← hpoint]
  exact br.settled r

/-! ## All-quotient event and geometry readout fragment -/

/-- Event and geometry readouts defined on every physical quotient state,
with fixed packet codomains per regulator, normal-form invariance,
functorial fine-to-coarse packet maps, and coarse naturality.

Variance: the consensus tower of `Tower.ConsensusTower` refines coarse to
fine, sending regulator-`r` data to regulator-`s` data when `r ≤ s`; the
bare consensus quotient and this fragment coarse-grain fine to coarse,
sending regulator-`s` data to regulator-`r` data through
`BareConsensusTower.coarse` and the packet maps. -/
structure EventGeometryReadoutFragment {ι : Type u} [Preorder ι]
    (B : BareConsensusTower ι) where
  /-- Fixed event-packet codomain at each regulator. -/
  EventPacket : ι → Type
  /-- Fixed geometry-packet codomain at each regulator. -/
  GeometryPacket : ι → Type
  /-- Event readout on every physical quotient state. -/
  eventRead : ∀ r, B.Quot r → EventPacket r
  /-- Geometry readout on every physical quotient state. -/
  geometryRead : ∀ r, B.Quot r → GeometryPacket r
  eventRead_normalForm : ∀ r q,
    eventRead r (B.normalForm r q) = eventRead r q
  geometryRead_normalForm : ∀ r q,
    geometryRead r (B.normalForm r q) = geometryRead r q
  /-- Fine-to-coarse event-packet map, in the direction of
  `BareConsensusTower.coarse`. -/
  eventCoarse : ∀ {r s}, r ≤ s → EventPacket s → EventPacket r
  /-- Fine-to-coarse geometry-packet map. -/
  geometryCoarse : ∀ {r s}, r ≤ s → GeometryPacket s → GeometryPacket r
  eventCoarse_refl : ∀ r p, eventCoarse (show r ≤ r from le_rfl) p = p
  eventCoarse_trans : ∀ {r s t} (hrs : r ≤ s) (hst : s ≤ t) p,
    eventCoarse hrs (eventCoarse hst p) = eventCoarse (hrs.trans hst) p
  geometryCoarse_refl : ∀ r p,
    geometryCoarse (show r ≤ r from le_rfl) p = p
  geometryCoarse_trans : ∀ {r s t} (hrs : r ≤ s) (hst : s ≤ t) p,
    geometryCoarse hrs (geometryCoarse hst p) =
      geometryCoarse (hrs.trans hst) p
  eventRead_natural : ∀ {r s} (hrs : r ≤ s) q,
    eventCoarse hrs (eventRead s q) = eventRead r (B.coarse hrs q)
  geometryRead_natural : ∀ {r s} (hrs : r ≤ s) q,
    geometryCoarse hrs (geometryRead s q) = geometryRead r (B.coarse hrs q)

namespace EventGeometryReadoutFragment

variable {ι : Type u} [Preorder ι] {B : BareConsensusTower ι}

/-- The event readout is a function of the normal form: quotient states
with a common normal form receive the same event packet. -/
theorem eventRead_eq_of_normalForm_eq (F : EventGeometryReadoutFragment B)
    (r : ι) {q q' : B.Quot r}
    (h : B.normalForm r q = B.normalForm r q') :
    F.eventRead r q = F.eventRead r q' := by
  rw [← F.eventRead_normalForm r q, h, F.eventRead_normalForm r q']

/-- The geometry readout is a function of the normal form. -/
theorem geometryRead_eq_of_normalForm_eq
    (F : EventGeometryReadoutFragment B) (r : ι) {q q' : B.Quot r}
    (h : B.normalForm r q = B.normalForm r q') :
    F.geometryRead r q = F.geometryRead r q' := by
  rw [← F.geometryRead_normalForm r q, h, F.geometryRead_normalForm r q']

/-- On a selected settled branch the event packets form a coherent family
over the packet coarse maps.  This equation is the all-quotient side of the
E7 source-link square; the matching side against the E6 selected event
world is a separate composition. -/
theorem eventRead_point_coherent (F : EventGeometryReadoutFragment B)
    (br : SelectedSettledBranch B) {r s : ι} (hrs : r ≤ s) :
    F.eventCoarse hrs (F.eventRead s (br.point s)) =
      F.eventRead r (br.point r) := by
  rw [F.eventRead_natural hrs, br.coherent hrs]

/-- On a selected settled branch the geometry packets form a coherent
family over the packet coarse maps. -/
theorem geometryRead_point_coherent (F : EventGeometryReadoutFragment B)
    (br : SelectedSettledBranch B) {r s : ι} (hrs : r ≤ s) :
    F.geometryCoarse hrs (F.geometryRead s (br.point s)) =
      F.geometryRead r (br.point r) := by
  rw [F.geometryRead_natural hrs, br.coherent hrs]

end EventGeometryReadoutFragment

/-! ## Einstein handoff packet -/

/-- Einstein handoff packet for the event and geometry readout fragment.

The packet certifies exactly the fragment data: fixed packet codomains, the
two readouts, normal-form invariance, functorial fine-to-coarse packet
maps, and coarse naturality over the physical quotient.

Obligations that remain with the Einstein readout completion (issue
`#694`): modular, stress, entropy, and scale readouts on the same repaired
quotient domain; the typed arrows with their commutation receipts
(`OPH.EinsteinBranch.TypedArrowCommutation`); refinement naturality for all
six readouts (`OPH.EinsteinBranch.ReadoutNaturality`); the boundary-fibre
receipt; assembly into `OPH.EinsteinBranch.EinsteinAdmissibleTower`; and
every analytic Einstein premise.  These obligations are docstring-level
boundary statements; none is a field and none is proved in this module. -/
structure ReadoutFragmentPacket {ι : Type u} [Preorder ι]
    (B : BareConsensusTower ι) where
  /-- The event and geometry fragment being handed off. -/
  fragment : EventGeometryReadoutFragment B

/-- Package an event and geometry fragment as the Einstein handoff. -/
def readoutFragmentPacket {ι : Type u} [Preorder ι]
    {B : BareConsensusTower ι} (F : EventGeometryReadoutFragment B) :
    ReadoutFragmentPacket B :=
  ⟨F⟩

/-! ## Boundary fibre from a complete event signature -/

/-- The signature factors through the protected boundary on consistent
states: equal boundary values force equal signature values. -/
def SignatureFactorsThroughBoundary {ι : Type u} [Preorder ι]
    (B : BareConsensusTower ι) {Signature : ι → Type v}
    (Sigma : ∀ r, B.Quot r → Signature r) : Prop :=
  ∀ r (x y : B.Quot r), B.consistent r x → B.consistent r y →
    B.boundary r x = B.boundary r y → Sigma r x = Sigma r y

/-- The signature is complete on consistent states: equal signature values
force equal quotient states. -/
def SignatureComplete {ι : Type u} [Preorder ι]
    (B : BareConsensusTower ι) {Signature : ι → Type v}
    (Sigma : ∀ r, B.Quot r → Signature r) : Prop :=
  ∀ r (x y : B.Quot r), B.consistent r x → B.consistent r y →
    Sigma r x = Sigma r y → x = y

/-- Boundary-fibre theorem from a complete event signature: a signature
family that factors through the protected boundary and is complete on
consistent states forces the boundary fibre of the tower.

Inadmissibility note: taking `Sigma := B.boundary` makes the completeness
hypothesis the conclusion itself
(`boundaryFiber_of_complete_signature_trivial_instance`), so the theorem is
substantive only for a signature constructed independently of the protected
boundary and proved complete on its own grounds. -/
theorem boundaryFiber_of_complete_signature {ι : Type u} [Preorder ι]
    (B : BareConsensusTower ι) {Signature : ι → Type v}
    (Sigma : ∀ r, B.Quot r → Signature r)
    (hfactor : SignatureFactorsThroughBoundary B Sigma)
    (hcomplete : SignatureComplete B Sigma) :
    BoundaryFiber B := by
  intro r x y hx hy hB
  exact hcomplete r x y hx hy (hfactor r x y hx hy hB)

/-- With the protected boundary as its own signature, the factoring
hypothesis holds trivially. -/
theorem signatureFactorsThroughBoundary_boundary {ι : Type u} [Preorder ι]
    (B : BareConsensusTower ι) :
    SignatureFactorsThroughBoundary B B.boundary :=
  fun _ _ _ _ _ hB => hB

/-- Circularity control: with the protected boundary as its own signature,
the completeness hypothesis is definitionally the boundary-fibre statement,
so `boundaryFiber_of_complete_signature` adds nothing in that instance.
The proof is `Iff.rfl`. -/
theorem boundaryFiber_of_complete_signature_trivial_instance {ι : Type u}
    [Preorder ι] (B : BareConsensusTower ι) :
    SignatureComplete B B.boundary ↔ BoundaryFiber B :=
  Iff.rfl

/-! ## Readout witness on the demo tower -/

/-- Readout witness on the demo tower: the event packet is the protected
bit and the geometry packet is a numeric weight of the same bit.  Both
readouts are normal-form invariant by construction and nonconstant across
quotient states. -/
def demoFragment : EventGeometryReadoutFragment demoTower where
  EventPacket := fun _ => Bool
  GeometryPacket := fun _ => ℕ
  eventRead := fun _ (q : Bool × Bool) => q.1
  geometryRead := fun _ (q : Bool × Bool) => if q.1 then 1 else 0
  eventRead_normalForm := fun _ _ => rfl
  geometryRead_normalForm := fun _ _ => rfl
  eventCoarse := fun _ => id
  geometryCoarse := fun _ => id
  eventCoarse_refl := fun _ _ => rfl
  eventCoarse_trans := fun _ _ _ => rfl
  geometryCoarse_refl := fun _ _ => rfl
  geometryCoarse_trans := fun _ _ _ => rfl
  eventRead_natural := fun _ _ => rfl
  geometryRead_natural := fun _ _ => rfl

/-- The witness event readout separates two explicit quotient states of the
demo tower, so the fragment interface admits nonconstant readouts. -/
theorem demoFragment_eventRead_separates :
    demoFragment.eventRead 0 (true, true) ≠
      demoFragment.eventRead 0 (false, false) := by
  intro h
  exact Bool.noConfusion h

/-- The witness geometry readout separates the same two quotient states. -/
theorem demoFragment_geometryRead_separates :
    demoFragment.geometryRead 0 (true, true) ≠
      demoFragment.geometryRead 0 (false, false) := by
  intro h
  have h10 : (1 : ℕ) = 0 := h
  exact Nat.one_ne_zero h10

/-- The separated states are both consistent, so the separation happens
inside the consistent sector of the quotient. -/
theorem demoFragment_separates_consistent_states :
    demoTower.consistent 0 (true, true) ∧
      demoTower.consistent 0 (false, false) ∧
      demoFragment.eventRead 0 (true, true) ≠
        demoFragment.eventRead 0 (false, false) :=
  ⟨rfl, rfl, demoFragment_eventRead_separates⟩

/-- Evaluating the witness fragment on the demo settled branch returns the
protected bit of the branch point at every regulator. -/
theorem demoFragment_eventRead_branch (r : ℕ) :
    demoFragment.eventRead r (demoSelectedBranch.point r) = true := rfl

/-- The Einstein handoff packet is inhabited on the demo tower. -/
def demoReadoutFragmentPacket : ReadoutFragmentPacket demoTower :=
  readoutFragmentPacket demoFragment

/-! ## Per-theorem axiom audit -/

#print axioms OPH.Tower.SelectedSettledBranch.mismatch_point_eq_zero
#print axioms OPH.Tower.SelectedSettledBranch.consistent_point
#print axioms OPH.Tower.SelectedSettledBranch.boundary_point_eq_of_settles_to
#print axioms OPH.Tower.SelectedSettledBranch.boundary_coarse_point
#print axioms OPH.Tower.SelectedSettledBranch.settled_coarse_point
#print axioms OPH.Tower.demoUnsettledCandidate_not_settled
#print axioms OPH.Tower.demoUnsettledCandidate_no_branch
#print axioms OPH.Tower.EventGeometryReadoutFragment.eventRead_eq_of_normalForm_eq
#print axioms OPH.Tower.EventGeometryReadoutFragment.geometryRead_eq_of_normalForm_eq
#print axioms OPH.Tower.EventGeometryReadoutFragment.eventRead_point_coherent
#print axioms OPH.Tower.EventGeometryReadoutFragment.geometryRead_point_coherent
#print axioms OPH.Tower.boundaryFiber_of_complete_signature
#print axioms OPH.Tower.signatureFactorsThroughBoundary_boundary
#print axioms OPH.Tower.boundaryFiber_of_complete_signature_trivial_instance
#print axioms OPH.Tower.demoFragment_eventRead_separates
#print axioms OPH.Tower.demoFragment_geometryRead_separates
#print axioms OPH.Tower.demoFragment_separates_consistent_states
#print axioms OPH.Tower.demoFragment_eventRead_branch

end OPH.Tower
