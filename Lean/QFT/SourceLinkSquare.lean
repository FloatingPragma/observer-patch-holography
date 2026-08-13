import Tower.EventGeometryReadout
import QFT.ObserverEventWorld

/-!
# The source-link square between all-quotient readouts and one selected world

This module constructs the finite source-link square between all-quotient
readouts and one selected observer world.

`WorldSlice` is the fixed readout codomain: the components of a finite
observer event world that can vary with a quotient state, namely the
active-event mask, the active-chart mask, and the packet assignment over
one fixed event-world signature.  The signature carriers are structural
data, so one slice type serves every regulator.
`FiniteEventGeometryWorld.slice` extracts the slice of a world.

`SelectedObserverReadout` is the deliverable-2 interface: a world slice
read from every physical quotient state of a bare consensus tower, with
the fixed codomain `WorldSlice sig` at every regulator, normal-form
invariance, and the reproduction law returning the slice of one declared
selected event world on a selected settled branch.  The access cut, the
tower and net data, the soldering, and the world laws are
quotient-independent structural parameters, so the slice carries the
whole quotient-varying content; a readout valued in proof-carrying world
bundles would add only constant components and proof-irrelevance
plumbing.

`SourceLinkSquare` is the deliverable-4 composition: it bundles a
selected settled branch, an all-quotient event and geometry readout
fragment, a selected observer readout, and the componentwise matching
equations stating that the fragment packets at the branch are the
encoded packet data of the selected world.  `eventRead_matches` and
`geometryRead_matches` derive the commutation of the square: along any
regulator pair the packet coarse image of the fine encoded selected
packet is the coarse encoded selected packet.  The normal-form stability
theorems reproduce the matched packets and the matched slice at the
normal form of the branch point and on its whole settling basin.

The square carries a two-tier hierarchy.  The plain square is declared
bookkeeping: nothing constrains the encodings, and
`SourceLinkSquare.ofConstantEncodings` inhabits the square for every
branch, fragment, and readout triple.  `InjectiveSourceLinkSquare` is
the determination tier: encoding injectivity is well typed and strictly
stronger than the bare matching equations, and under it
`InjectiveSourceLinkSquare.packetOf_determined` and
`InjectiveSourceLinkSquare.geometryRead_determined` recover the world
packet assignment and coordinate readout from the fragment values on
the branch.

Typing note.  A literal equation between a fragment packet and the world
packet data is not well typed: `fragment.EventPacket r` is an arbitrary
per-regulator codomain while the world data lives on the fixed signature
carriers.  The declared encodings `eventEncode` and `geometryEncode`
state equality after a declared encoding of the world data into each
regulator codomain; whether that equality carries content depends on the
encodings, which is what the injective tier certifies.  The two towers
also carry opposite index variance: the consensus tower of the event
world refines coarse to fine while the bare tower coarse-grains fine to
coarse, and the two index types are unrelated.  The selected world
therefore enters at one fixed consensus-tower regulator through its
slice and packet data alone; no equation relates the consensus
refinement maps to the bare coarse maps.

## Claim boundary

The square is finite bookkeeping tying the all-quotient readouts to one
selected observer world.  No source produces the branch, no mechanism
selects the world, no probability law is stated, and no continuum claim
follows. Einstein readouts beyond the event and geometry fragment are scoped
to gravitation composition (#729), with physical readout instrumentation
scoped separately to #737.
-/

namespace OPH.QFT

open OPH.Tower
open OPH.EinsteinBranch

universe u u' v

/-! ## World slices -/

/-- The quotient-varying slice of a finite observer event world over one
fixed signature: the active-event mask, the active-chart mask, and the
packet assignment.  The signature carriers are structural data, so the
slice type is one fixed codomain for every regulator of a bare consensus
tower. -/
structure WorldSlice (sig : EventWorldSignature) where
  /-- The active-event mask of the slice. -/
  eventMask : sig.Event → Prop
  /-- The active-chart mask of the slice. -/
  chartMask : sig.Chart → Prop
  /-- The packet assignment of the slice. -/
  packetAssign : sig.Event → sig.Packet

variable {ι : Type u} [Preorder ι]
variable {κ : Type u'} [Preorder κ]
variable {T : ConsensusTower κ} {N : FiniteCausalObserverNet T} {r0 : κ}
variable {A : ObserverAccessCut T N r0} {sig : EventWorldSignature}
variable {S : OPH.C2Soldering.EventFrameSoldering sig.Event sig.Chart}

/-- The slice of a finite observer event world: its active-event mask,
its active-chart mask, and its packet assignment.  Every other world
component (regions, ownership, algebras, soldering, laws) is
quotient-independent structural data and stays outside the slice. -/
def FiniteEventGeometryWorld.slice (W : FiniteEventGeometryWorld A sig S) :
    WorldSlice sig where
  eventMask := W.activeEvent
  chartMask := W.activeChart
  packetAssign := W.packetOf

/-! ## Deliverable 2: the selected observer readout -/

/-- A selected observer readout over one bare consensus tower: a world
slice read from every physical quotient state, with the fixed codomain
`WorldSlice sig` at every regulator, normal-form invariance, and the
reproduction law on the selected settled branch, where the readout
returns the slice of the declared selected event world.

The selected world sits at one fixed regulator of its own consensus
tower; the two regulator index types are unrelated, and the world enters
only through its slice. -/
structure SelectedObserverReadout (B : BareConsensusTower ι)
    (br : SelectedSettledBranch B)
    (W : FiniteEventGeometryWorld A sig S) where
  /-- The slice read at each physical quotient state. -/
  read : ∀ r, B.Quot r → WorldSlice sig
  /-- Normal-form invariance of the readout. -/
  read_normalForm : ∀ r q, read r (B.normalForm r q) = read r q
  /-- On the selected settled branch the readout returns the slice of
  the selected world at every regulator. -/
  read_point : ∀ r, read r (br.point r) = W.slice

namespace SelectedObserverReadout

variable {B : BareConsensusTower ι} {br : SelectedSettledBranch B}
variable {W : FiniteEventGeometryWorld A sig S}

/-- The readout is a function of the normal form: quotient states with a
common normal form receive the same slice. -/
theorem read_eq_of_normalForm_eq (R : SelectedObserverReadout B br W)
    (r : ι) {q q' : B.Quot r}
    (h : B.normalForm r q = B.normalForm r q') :
    R.read r q = R.read r q' := by
  rw [← R.read_normalForm r q, h, R.read_normalForm r q']

/-- Every quotient state that settles to the branch point reads out the
selected world slice. -/
theorem read_of_settles_to (R : SelectedObserverReadout B br W) (r : ι)
    {q : B.Quot r} (h : B.normalForm r q = br.point r) :
    R.read r q = W.slice := by
  rw [← R.read_normalForm r q, h, R.read_point r]

/-- On the branch the readout is regulator-independent: every branch
point reads out the same selected slice. -/
theorem read_point_regulator_independent
    (R : SelectedObserverReadout B br W) (r s : ι) :
    R.read r (br.point r) = R.read s (br.point s) := by
  rw [R.read_point r, R.read_point s]

end SelectedObserverReadout

/-! ## Deliverable 4: the source-link square -/

/-- The source-link square: a selected settled branch, an all-quotient
event and geometry readout fragment, a selected observer readout, and
the componentwise matching equations tying the fragment packets at the
branch to the packet data of the selected world.

`eventEncode` and `geometryEncode` are declared encodings of the world
packet data into the fragment codomains; a literal equation between
`fragment.EventPacket r` values and functions on the signature carriers
is not well typed, so the matching equations hold through the encodings.
The plain square is the declared-bookkeeping tier: nothing constrains
the encodings, `ofConstantEncodings` inhabits the square for every
branch, fragment, and readout triple, and the matching equations alone
determine nothing about the world.  The injective strengthening
`InjectiveSourceLinkSquare` is well typed and strictly stronger; under
it the matching determines the world packet data from the fragment
values on the branch (`InjectiveSourceLinkSquare.packetOf_determined`). -/
structure SourceLinkSquare (B : BareConsensusTower ι)
    (W : FiniteEventGeometryWorld A sig S) where
  /-- The selected settled branch of the bare tower. -/
  branch : SelectedSettledBranch B
  /-- The all-quotient event and geometry readout fragment. -/
  fragment : EventGeometryReadoutFragment B
  /-- The selected observer readout against the branch and the world. -/
  readout : SelectedObserverReadout B branch W
  /-- Declared encoding of the world event-packet data into the fragment
  event codomain at each regulator. -/
  eventEncode : ∀ r : ι, (sig.Event → sig.Packet) → fragment.EventPacket r
  /-- Declared encoding of the world chart-geometry data into the
  fragment geometry codomain at each regulator. -/
  geometryEncode : ∀ r : ι,
    (sig.Chart → sig.Packet → OPH.C1Lorentz.Herm2) →
      fragment.GeometryPacket r
  /-- Componentwise event matching at the branch: the fragment event
  packet at the branch point is the encoded packet assignment of the
  selected world. -/
  event_matches : ∀ r,
    fragment.eventRead r (branch.point r) = eventEncode r W.packetOf
  /-- Componentwise geometry matching at the branch: the fragment
  geometry packet at the branch point is the encoded chart geometry
  readout of the selected world. -/
  geometry_matches : ∀ r,
    fragment.geometryRead r (branch.point r) =
      geometryEncode r W.geometryRead

namespace SourceLinkSquare

variable {B : BareConsensusTower ι} {W : FiniteEventGeometryWorld A sig S}

/-- Commutation of the square, event side: along any regulator pair the
packet coarse image of the fine encoded selected packet is the coarse
encoded selected packet.  The proof composes fragment naturality with
branch coherence through the matching equations, so this is
`EventGeometryReadoutFragment.eventRead_point_coherent` transported to
the matched packets. -/
theorem eventRead_matches (Q : SourceLinkSquare B W) {r s : ι}
    (hrs : r ≤ s) :
    Q.fragment.eventCoarse hrs (Q.eventEncode s W.packetOf) =
      Q.eventEncode r W.packetOf := by
  rw [← Q.event_matches s, ← Q.event_matches r]
  exact Q.fragment.eventRead_point_coherent Q.branch hrs

/-- Commutation of the square, geometry side. -/
theorem geometryRead_matches (Q : SourceLinkSquare B W) {r s : ι}
    (hrs : r ≤ s) :
    Q.fragment.geometryCoarse hrs (Q.geometryEncode s W.geometryRead) =
      Q.geometryEncode r W.geometryRead := by
  rw [← Q.geometry_matches s, ← Q.geometry_matches r]
  exact Q.fragment.geometryRead_point_coherent Q.branch hrs

/-- Normal-form stability of the matched event packet: the fragment
event packet at the normal form of the branch point is the encoded
selected packet. -/
theorem eventRead_matches_normalForm (Q : SourceLinkSquare B W) (r : ι) :
    Q.fragment.eventRead r (B.normalForm r (Q.branch.point r)) =
      Q.eventEncode r W.packetOf := by
  rw [Q.branch.settled r]
  exact Q.event_matches r

/-- Normal-form stability of the matched geometry packet. -/
theorem geometryRead_matches_normalForm (Q : SourceLinkSquare B W)
    (r : ι) :
    Q.fragment.geometryRead r (B.normalForm r (Q.branch.point r)) =
      Q.geometryEncode r W.geometryRead := by
  rw [Q.branch.settled r]
  exact Q.geometry_matches r

/-- Normal-form stability of the matched slice: the readout at the
normal form of the branch point is the selected world slice. -/
theorem read_matches_normalForm (Q : SourceLinkSquare B W) (r : ι) :
    Q.readout.read r (B.normalForm r (Q.branch.point r)) = W.slice := by
  rw [Q.branch.settled r]
  exact Q.readout.read_point r

/-- Settling-basin form, event side: every quotient state that settles
to the branch point carries the encoded selected event packet. -/
theorem eventRead_matches_of_settles_to (Q : SourceLinkSquare B W)
    (r : ι) {q : B.Quot r} (h : B.normalForm r q = Q.branch.point r) :
    Q.fragment.eventRead r q = Q.eventEncode r W.packetOf := by
  rw [← Q.fragment.eventRead_normalForm r q, h]
  exact Q.event_matches r

/-- Settling-basin form, geometry side. -/
theorem geometryRead_matches_of_settles_to (Q : SourceLinkSquare B W)
    (r : ι) {q : B.Quot r} (h : B.normalForm r q = Q.branch.point r) :
    Q.fragment.geometryRead r q = Q.geometryEncode r W.geometryRead := by
  rw [← Q.fragment.geometryRead_normalForm r q, h]
  exact Q.geometry_matches r

/-- The square on the settling basin: every quotient state that settles
to the branch point reads out the selected slice, the encoded selected
event packet, and the encoded selected geometry packet. -/
theorem matches_of_settles_to (Q : SourceLinkSquare B W) (r : ι)
    {q : B.Quot r} (h : B.normalForm r q = Q.branch.point r) :
    Q.readout.read r q = W.slice ∧
      Q.fragment.eventRead r q = Q.eventEncode r W.packetOf ∧
      Q.fragment.geometryRead r q = Q.geometryEncode r W.geometryRead :=
  ⟨Q.readout.read_of_settles_to r h,
    Q.eventRead_matches_of_settles_to r h,
    Q.geometryRead_matches_of_settles_to r h⟩

/-- Tier-1 saturation: constant encodings at the branch values close the
matching equations for every branch, fragment, and readout triple.  The
plain square therefore carries no constraint beyond its three components;
determination lives in the injective tier below. -/
noncomputable def ofConstantEncodings (br : SelectedSettledBranch B)
    (frag : EventGeometryReadoutFragment B)
    (Rd : SelectedObserverReadout B br W) : SourceLinkSquare B W where
  branch := br
  fragment := frag
  readout := Rd
  eventEncode := fun r _ => frag.eventRead r (br.point r)
  geometryEncode := fun r _ => frag.geometryRead r (br.point r)
  event_matches := fun _ => rfl
  geometry_matches := fun _ => rfl

/-- Every branch, fragment, and readout triple inhabits the plain
square. -/
theorem nonempty_of_readout (br : SelectedSettledBranch B)
    (frag : EventGeometryReadoutFragment B)
    (Rd : SelectedObserverReadout B br W) :
    Nonempty (SourceLinkSquare B W) :=
  ⟨ofConstantEncodings br frag Rd⟩

end SourceLinkSquare

/-! ## The injective strengthening: the determination tier -/

/-- A source-link square whose declared encodings are injective at every
regulator.  Injectivity is the well-typed strengthening that the plain
square omits: it turns the matching equations from declared bookkeeping
into determination, because an injective encoding admits at most one
preimage for the fragment value at the branch. -/
structure InjectiveSourceLinkSquare (B : BareConsensusTower ι)
    (W : FiniteEventGeometryWorld A sig S) extends SourceLinkSquare B W where
  /-- The event encoding is injective at every regulator. -/
  eventEncode_injective : ∀ r, Function.Injective (eventEncode r)
  /-- The geometry encoding is injective at every regulator. -/
  geometryEncode_injective : ∀ r, Function.Injective (geometryEncode r)

namespace InjectiveSourceLinkSquare

variable {B : BareConsensusTower ι} {W W' : FiniteEventGeometryWorld A sig S}

/-- Determination, event side: through an injective encoding the fragment
event packet at the branch pins the world packet assignment.  Every packet
assignment whose encoding matches the branch value is `W.packetOf`. -/
theorem packetOf_determined (Q : InjectiveSourceLinkSquare B W) (r : ι)
    (p : sig.Event → sig.Packet)
    (h : Q.fragment.eventRead r (Q.branch.point r) = Q.eventEncode r p) :
    p = W.packetOf :=
  Q.eventEncode_injective r (h.symm.trans (Q.event_matches r))

/-- Determination, geometry side: through an injective encoding the
fragment geometry packet at the branch pins the world coordinate
readout. -/
theorem geometryRead_determined (Q : InjectiveSourceLinkSquare B W) (r : ι)
    (g : sig.Chart → sig.Packet → OPH.C1Lorentz.Herm2)
    (h : Q.fragment.geometryRead r (Q.branch.point r) =
      Q.geometryEncode r g) :
    g = W.geometryRead :=
  Q.geometryEncode_injective r (h.symm.trans (Q.geometry_matches r))

/-- Two worlds matched through one injective square carry equal world
data: if a second world's packet assignment and coordinate readout encode
to the same branch values, they equal the selected world's packet
assignment and coordinate readout on every event, chart, and packet. -/
theorem world_data_determined (Q : InjectiveSourceLinkSquare B W) (r : ι)
    (hev : Q.fragment.eventRead r (Q.branch.point r) =
      Q.eventEncode r W'.packetOf)
    (hgeo : Q.fragment.geometryRead r (Q.branch.point r) =
      Q.geometryEncode r W'.geometryRead) :
    W'.packetOf = W.packetOf ∧ W'.geometryRead = W.geometryRead :=
  ⟨Q.packetOf_determined r W'.packetOf hev,
    Q.geometryRead_determined r W'.geometryRead hgeo⟩

end InjectiveSourceLinkSquare

/-! ## The demo squares

The witnesses glue the two committed demo stacks: the bare demo tower
with its selected branch and readout fragments, and the dimension-two
witness event world over the witness access cut.  The two stacks keep
unrelated regulator index types, `ℕ` for the bare tower and `Unit` for
the witness consensus tower; the squares compose them without any shared
index, as the interface requires.  Three squares mark the tiers: the
constant-geometry square is the labeled degenerate control on its
geometry side, `demoNonconstantGeometrySquare` matches the geometry
through a nonconstant value encoding with a separation receipt, and
`demoInjectiveSourceLinkSquare` inhabits the injective determination
tier and exercises the uniqueness theorems on explicit data. -/

section DemoSquare

/-- The active demo slice: every event and chart active, identity packet
assignment.  This is the slice of the witness event world. -/
def demoActiveSlice : WorldSlice witnessEventSignature where
  eventMask := fun _ => True
  chartMask := fun _ => True
  packetAssign := id

/-- The inert demo slice: no active event or chart, constant packet
assignment.  It differs from the active slice and drives the separation
fact. -/
def demoInertSlice : WorldSlice witnessEventSignature where
  eventMask := fun _ => False
  chartMask := fun _ => False
  packetAssign := fun _ => 0

/-- The demo observer readout on the bare demo tower against the witness
event world: states with protected bit `true` read the active slice,
states with protected bit `false` read the inert slice.  The readout
depends only on the protected bit, so normal-form invariance and branch
reproduction hold definitionally. -/
noncomputable def demoObserverReadout :
    SelectedObserverReadout demoTower demoSelectedBranch
      witnessEventWorld where
  read := fun _ q => cond q.1 demoActiveSlice demoInertSlice
  read_normalForm := fun _ _ => rfl
  read_point := fun _ => rfl

/-- The demo readout separates two quotient states of the demo tower:
the two slices differ in their packet assignment at the second event. -/
theorem demoObserverReadout_separates :
    demoObserverReadout.read 0 (true, true) ≠
      demoObserverReadout.read 0 (false, false) := by
  intro h
  have h1 : (1 : Fin 2) = 0 :=
    congrFun (congrArg WorldSlice.packetAssign h) 1
  exact absurd h1 (by decide)

/-- Settling-basin evaluation on explicit data: the inconsistent state
`(true, false)` settles to the branch point and reads out the selected
world slice. -/
theorem demoObserverReadout_read_settling_state :
    demoObserverReadout.read 0 (true, false) = witnessEventWorld.slice :=
  demoObserverReadout.read_of_settles_to 0 rfl

/-- The demo source-link square: demo branch, demo fragment, demo
observer readout, an event encoding probing the packet assignment at the
second event, and the constant geometry adaptor at the matched value.
The matching equations hold definitionally: the fragment reads the
protected bit `true` and weight `1` at the branch point, the encoded
selected event packet is `true`, and the encoded geometry value is `1`.
The geometry side is the labeled degenerate control of the demo family;
`demoNonconstantGeometrySquare` below matches the geometry through a
nonconstant encoding. -/
noncomputable def demoSourceLinkSquare :
    SourceLinkSquare demoTower witnessEventWorld where
  branch := demoSelectedBranch
  fragment := demoFragment
  readout := demoObserverReadout
  eventEncode := fun _ p => decide (p 1 = 1)
  geometryEncode := fun _ _ => (1 : ℕ)
  event_matches := fun _ => rfl
  geometry_matches := fun _ => rfl

/-- Separation fact for the matching data: the demo event encoding
distinguishes the packet data of the selected world from the constant
packet map, so the matched packet carries world data rather than one
adaptor value.  The geometry encoding of the demo square is the constant
adaptor at the matched value; the event side carries the separation. -/
theorem demoSourceLinkSquare_eventEncode_separates :
    demoSourceLinkSquare.eventEncode 0 witnessEventWorld.packetOf ≠
      demoSourceLinkSquare.eventEncode 0 (fun _ => 0) := by
  intro h
  have hb : true = false := h
  exact Bool.noConfusion hb

/-- Commutation of the demo square evaluated on the explicit regulator
pair `0 ≤ 1`: the coarse image of the fine encoded selected packet is
the coarse encoded selected packet. -/
theorem demoSourceLinkSquare_eventRead_matches :
    demoSourceLinkSquare.fragment.eventCoarse (Nat.zero_le 1)
        (demoSourceLinkSquare.eventEncode 1 witnessEventWorld.packetOf) =
      demoSourceLinkSquare.eventEncode 0 witnessEventWorld.packetOf :=
  demoSourceLinkSquare.eventRead_matches (Nat.zero_le 1)

/-- The demo square on the settling basin at the explicit off-branch
state `(true, false)`: the readout returns the selected slice and both
fragment packets return the encoded selected data. -/
theorem demoSourceLinkSquare_matches_of_settling_state :
    demoSourceLinkSquare.readout.read 0 (true, false) =
        witnessEventWorld.slice ∧
      demoSourceLinkSquare.fragment.eventRead 0 (true, false) =
        demoSourceLinkSquare.eventEncode 0 witnessEventWorld.packetOf ∧
      demoSourceLinkSquare.fragment.geometryRead 0 (true, false) =
        demoSourceLinkSquare.geometryEncode 0
          witnessEventWorld.geometryRead :=
  demoSourceLinkSquare.matches_of_settles_to 0 rfl

/-! ### The nonconstant geometry encoding -/

/-- The demo geometry value encoding: read the coordinate function at the
designated chart/packet pair `(0, 1)` and inject the designated demo
value pair into `ℕ`, `0` for the zero coordinate and `1` for every other
coordinate.  On the demo's two-element value set at the designated pair,
zero and the witness coordinate, this assignment is injective. -/
noncomputable def demoGeometryValueEncode
    (g : Fin 2 → Fin 2 → OPH.C1Lorentz.Herm2) : ℕ :=
  letI := Classical.dec (g 0 1 = 0)
  if g 0 1 = 0 then 0 else 1

/-- The designated geometry value of the selected world is nonzero: at
the chart/packet pair `(0, 1)` the witness coordinate has scalar part
`1`. -/
theorem witnessEventWorld_designated_geometry_ne_zero :
    witnessEventWorld.geometryRead 0 1 ≠ 0 := by
  intro h
  have h1 := congrArg Prod.fst h
  norm_num [witnessEventWorld, OPH.C2Soldering.controlSoldering,
    OPH.C2Soldering.controlEventPoint, OPH.C2Soldering.controlChartOffset,
    Prod.fst_add] at h1

/-- The value encoding sends the selected world's coordinate readout to
`1`. -/
theorem demoGeometryValueEncode_world :
    demoGeometryValueEncode witnessEventWorld.geometryRead = 1 := by
  unfold demoGeometryValueEncode
  exact if_neg witnessEventWorld_designated_geometry_ne_zero

/-- The value encoding sends the constant zero readout to `0`. -/
theorem demoGeometryValueEncode_zero :
    demoGeometryValueEncode (fun _ _ => 0) = 0 := by
  unfold demoGeometryValueEncode
  exact if_pos rfl

/-- The demo square with a nonconstant geometry encoding: demo branch,
demo fragment, demo observer readout, the packet-probing event encoding,
and the geometry value encoding at the designated pair.  The geometry
matching carries coordinate content: the fragment weight `1` at the
branch equals the encoded coordinate value of the selected world because
the designated witness coordinate is nonzero. -/
noncomputable def demoNonconstantGeometrySquare :
    SourceLinkSquare demoTower witnessEventWorld where
  branch := demoSelectedBranch
  fragment := demoFragment
  readout := demoObserverReadout
  eventEncode := fun _ p => decide (p 1 = 1)
  geometryEncode := fun _ g => demoGeometryValueEncode g
  event_matches := fun _ => rfl
  geometry_matches := fun _ => demoGeometryValueEncode_world.symm

/-- Separation receipt for the geometry matching: the value encoding
distinguishes the selected world's coordinate readout from the constant
zero readout, so the matched geometry packet carries coordinate data of
the world rather than one adaptor value. -/
theorem demoNonconstantGeometrySquare_geometryEncode_separates :
    demoNonconstantGeometrySquare.geometryEncode 0
        witnessEventWorld.geometryRead ≠
      demoNonconstantGeometrySquare.geometryEncode 0 (fun _ _ => 0) := by
  intro h
  have h1 : demoGeometryValueEncode witnessEventWorld.geometryRead =
      demoGeometryValueEncode (fun _ _ => 0) := h
  rw [demoGeometryValueEncode_world, demoGeometryValueEncode_zero] at h1
  exact Nat.one_ne_zero h1

/-! ### The injective demo square -/

/-- An injective-codomain readout fragment on the demo tower: the event
packet codomain is the packet-assignment function space and the geometry
packet codomain is the coordinate-readout function space, so the identity
encodings into them are injective.  States with protected bit `true`
read the selected world data; states with protected bit `false` read the
inert constants. -/
noncomputable def demoInjectiveFragment :
    EventGeometryReadoutFragment demoTower where
  EventPacket := fun _ => Fin 2 → Fin 2
  GeometryPacket := fun _ => Fin 2 → Fin 2 → OPH.C1Lorentz.Herm2
  eventRead := fun _ q => cond q.1 id fun _ => 0
  geometryRead := fun _ q =>
    cond q.1 witnessEventWorld.geometryRead fun _ _ => 0
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

/-- The injective demo square: demo branch, injective-codomain fragment,
demo observer readout, and identity encodings on both sides.  The
matching equations hold definitionally and both encodings are injective,
so the square inhabits the determination tier. -/
noncomputable def demoInjectiveSourceLinkSquare :
    InjectiveSourceLinkSquare demoTower witnessEventWorld where
  branch := demoSelectedBranch
  fragment := demoInjectiveFragment
  readout := demoObserverReadout
  eventEncode := fun _ p => p
  geometryEncode := fun _ g => g
  event_matches := fun _ => rfl
  geometry_matches := fun _ => rfl
  eventEncode_injective := fun _ => Function.injective_id
  geometryEncode_injective := fun _ => Function.injective_id

/-- Determination exercised on explicit data: every packet assignment
whose encoding matches the injective fragment at the branch is the
selected world's packet assignment. -/
theorem demoInjectiveSourceLinkSquare_packetOf_determined
    (p : Fin 2 → Fin 2)
    (h : demoInjectiveSourceLinkSquare.fragment.eventRead 0
        (demoInjectiveSourceLinkSquare.branch.point 0) =
      demoInjectiveSourceLinkSquare.eventEncode 0 p) :
    p = witnessEventWorld.packetOf :=
  demoInjectiveSourceLinkSquare.packetOf_determined 0 p h

/-- The injective matching rejects the constant packet assignment: its
encoding differs from the fragment value at the branch, because
determination would force it to equal the identity packet map. -/
theorem demoInjectiveSourceLinkSquare_rejects_constant :
    demoInjectiveSourceLinkSquare.fragment.eventRead 0
        (demoInjectiveSourceLinkSquare.branch.point 0) ≠
      demoInjectiveSourceLinkSquare.eventEncode 0 (fun _ => 0) := by
  intro h
  have h1 : (fun _ : Fin 2 => (0 : Fin 2)) = witnessEventWorld.packetOf :=
    demoInjectiveSourceLinkSquare.packetOf_determined 0 _ h
  have h2 : (0 : Fin 2) = 1 := congrFun h1 1
  exact absurd h2 (by decide)

end DemoSquare

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.SelectedObserverReadout.read_eq_of_normalForm_eq
#print axioms OPH.QFT.SelectedObserverReadout.read_of_settles_to
#print axioms OPH.QFT.SelectedObserverReadout.read_point_regulator_independent
#print axioms OPH.QFT.SourceLinkSquare.eventRead_matches
#print axioms OPH.QFT.SourceLinkSquare.geometryRead_matches
#print axioms OPH.QFT.SourceLinkSquare.eventRead_matches_normalForm
#print axioms OPH.QFT.SourceLinkSquare.geometryRead_matches_normalForm
#print axioms OPH.QFT.SourceLinkSquare.read_matches_normalForm
#print axioms OPH.QFT.SourceLinkSquare.eventRead_matches_of_settles_to
#print axioms OPH.QFT.SourceLinkSquare.geometryRead_matches_of_settles_to
#print axioms OPH.QFT.SourceLinkSquare.matches_of_settles_to
#print axioms OPH.QFT.SourceLinkSquare.nonempty_of_readout
#print axioms OPH.QFT.InjectiveSourceLinkSquare.packetOf_determined
#print axioms OPH.QFT.InjectiveSourceLinkSquare.geometryRead_determined
#print axioms OPH.QFT.InjectiveSourceLinkSquare.world_data_determined
#print axioms OPH.QFT.demoObserverReadout_separates
#print axioms OPH.QFT.demoObserverReadout_read_settling_state
#print axioms OPH.QFT.demoSourceLinkSquare_eventEncode_separates
#print axioms OPH.QFT.demoSourceLinkSquare_eventRead_matches
#print axioms OPH.QFT.demoSourceLinkSquare_matches_of_settling_state
#print axioms OPH.QFT.witnessEventWorld_designated_geometry_ne_zero
#print axioms OPH.QFT.demoGeometryValueEncode_world
#print axioms OPH.QFT.demoGeometryValueEncode_zero
#print axioms OPH.QFT.demoNonconstantGeometrySquare_geometryEncode_separates
#print axioms OPH.QFT.demoInjectiveSourceLinkSquare_packetOf_determined
#print axioms OPH.QFT.demoInjectiveSourceLinkSquare_rejects_constant
