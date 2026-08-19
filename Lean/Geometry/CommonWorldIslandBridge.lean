import Geometry.CommonWorldKinematicsWitness
import Geometry.RestFiberShellTransport

set_option autoImplicit false

open scoped BigOperators

/-!
# CW1 cross-island bridge: one typed dictionary between the two witness islands (issue #740)

WHAT IS PROVED.  The committed common-world witness of
`Geometry/CommonWorldKinematicsWitness.lean` joins its kinematics island and
its screen island only by conjunction on one record.  This module adds
exactly one cross-island join: a typed dictionary between one committed
screen-derived object and the rest fiber of the record's declared frame on
the committed Hermitian Lorentz module.

The threaded screen object, named exactly.  The committed witness's screen
island carries a declared Maxwell evolution bundle and a declared
electroweak premise bundle; neither field carries the rank-three Gram
quotient itself.  The object threaded through the bridge is the committed
rank-three source Gram quotient `FrameQuotient` of
`Screen/PrimitivePortFrameQuotient.lean`: the quotient of the six-axis
signed-record control module by the radical of the finite source Gram form,
whose entries are read from the registered twelve-port incidence table on
the same committed twelve-port index type that the witness's screen join
evaluates.  This quotient is the exact committed screen object that admits
a committed isometry to a rest fiber: `frameQuotientEquivStandardRest` of
`Geometry/SpatialReadbackSoldering.lean` at the standard frame, and the
boost-shear isometry `restFiberIsometry` of
`Geometry/RestFiberShellTransport.lean` at every declared frame.

The bridge.  For every declared frame `u`, `quotientRestBridge u` is the
composition of the committed first-isomorphism equivalence
`quotientEquivVec3`, the canonical linear equivalence onto the Euclidean
completion carrier, and the committed boost-shear isometry
`restFiberIsometry u`.  Its metric identity `quotientRestBridge_metric`
states that the committed rest metric pulled back through the bridge is
exactly the committed screen quotient Gram; at the standard frame the
bridge equals the committed candidate bridge
`frameQuotientEquivStandardRest` (`quotientRestBridge_standardFrame`); the
declared radial readout of the rest fiber reads the screen Gram exactly
(`quotientRestBridge_radius`); and the bridge carries the class of every
integer record control to the committed screen chain image `sourceChain`
scaled by the committed normalization radius
(`quotientRestBridge_sourceChain`).  Every conjunct reuses a committed
theorem; nothing is reproved.

The extended witness.  `BridgedCommonWorldArchitecture` extends the
committed `CommonWorldArchitecture` by two fields: one declared linear
equivalence `bridge` from the screen Gram quotient onto the rest fiber of
the record's own frame field, and the metric clause `bridge_metric` that
pins the pullback of the rest metric to the screen quotient Gram.  The
composed receipt `bridgedCommonWorld_receipt` consumes one typed antecedent
bundle (`BridgedWorldAntecedents`: one extended record and one declared
radial flux packet) and proves, of the one record simultaneously, all eight
committed island packets of the base witness plus the bridge clause
`CarriesIslandBridge` (metric equality, zero-preservation both ways, the
matching dimension three on both sides, the radius readout, and the
commuting square with every committed oriented Lorentz chart transport) and
the bridged shell readout `CarriesBridgedShellReadout` (the committed
transported inverse-square law read on the screen Gram: at every nonzero
quotient class the strength is the charge over `4 π` times the class's own
Gram square).  The explicit inhabitant `bridgedCommittedWitness` extends
the committed inhabitant `committedWitness` by the canonical bridge, with
`bridgedCommittedWitness_extends` a definitional equality.

The load-bearing receipt.  `BareBridgedArchitecture` is the same extension
with the metric clause dropped.  The inhabitant `mismatchedExtension`
extends the same committed base witness by the doubled dictionary
`doubledBridge` (the canonical bridge composed with scaling by two), a
perfectly typed linear equivalence between the same two objects;
`mismatchedExtension_breaks_metric` exhibits a quotient class at which its
rest-metric pullback is four while the screen Gram is one.  So dropping the
metric clause readmits a mismatched pair, and
`bridge_metric_excludes_doubling` shows the retained clause excludes that
dictionary from every extended record.

DECLARED DATA.  Every base-record field stays a declared hypothesis
exactly as in the committed witness.  The bridge field of the extended
record is declared data constrained by the metric clause; the canonical
bridge is a function of the frame alone and no source theorem selects the
frame.  The flux field of the antecedent bundle is the committed
`RadialFluxData` premise packet: its strength field is the declared radial
readout (premise register row PR-18) and its normalized field is the
declared shell-flux normalization (premise register row PR-19), consumed
exactly as in the committed shell-transport claim; `unitFlux` is one
explicit demonstration instance whose parameter values select nothing.

FALSIFIER.  The module fails if the canonical bridge violates the metric
identity at some frame and pair of quotient classes, if it differs from the
committed candidate bridge at the standard frame, if the extended record
type is uninhabited over the committed inhabitant, if a claimed committed
receipt fails at the record's fields, if the doubled dictionary satisfies
the metric clause at the exhibited class, or if some clause of the composed
receipt fails at the antecedent bundle.

NONCLAIMS AND BOUNDARY (the joins that REMAIN missing; CW1 stays open).
The witness grows by exactly one join.  The bridge is a candidate local
readout dictionary between two committed mathematical objects, not a
physical identification: no observer-produced physical spacetime, rod,
or position attachment is claimed, and no absolute space exists on either
side (PR-52 open).  The bridge is not unique: composing it with any
rest-metric-preserving automorphism satisfies the same clause, and no
committed theorem selects one representative.  The remaining cross-island
joins, named exactly: no shared clock or calibration exists, since no map
identifies the Maxwell step index with the kinematic proper time and no
physical clock or laboratory unit attaches to either parameter (PR-15
open); no common action exists, since the assembled screen action and the
declared kinematic flow share no variational principle on one carrier
(PR-54 open); no common dynamics exists, since the bridge is a static
dictionary with no time index on either side and no equation of motion is
transported through it.  The twelve ports, thirty seams, and twenty faces
of the screen carrier themselves have no map to points, intervals,
or cones of the declared module: the bridge carries only the derived
rank-three control quotient, so the carrier dictionary row PR-53 stays
open.  The committed normalization radius `rawRadius` in the source-chain
identity is the committed unit-diagonal Gram normalization, not a physical
length calibration.  Nothing here closes issue #740, and no
observation-ledger row is promoted: OL-N1 stays owed.

Axiom audit.  Every proof composes committed receipts with exact finite
mathematics; the module adds no project axiom and uses no native decision
procedure.  The guard lines at the end of the file show at most `propext`,
`Classical.choice`, and `Quot.sound`.
-/

namespace OPH.CommonWorldIslandBridge

open OPH.C1Lorentz OPH.C2Soldering OPH.CommonWorld
open OPH.RestFiberShellTransport OPH.InverseSquareShellLaw
open OPH.PrimitivePortFrameQuotient
open Real Module

noncomputable section

variable {n k : ℕ}

/-! ## The typed dictionary from the screen Gram quotient to every rest fiber -/

/-- Canonical linear equivalence from the product-coordinate presentation of
the screen quotient onto the committed Euclidean completion carrier.  The
underlying function is the identity of the `WithLp` synonym. -/
def vec3LinearEquivEuclidean : Vec3 ≃ₗ[ℝ] EuclideanVec3 :=
  (WithLp.linearEquiv 2 ℝ (Fin 3 → ℝ)).symm

@[simp] theorem vec3LinearEquivEuclidean_apply (v : Vec3) (i : Fin 3) :
    (vec3LinearEquivEuclidean v) i = v i := rfl

/-- The cross-island bridge: the committed rank-three screen Gram quotient
carried onto the rest fiber of one declared frame, as the composition of the
committed first-isomorphism equivalence `quotientEquivVec3`, the canonical
carrier equivalence, and the committed boost-shear isometry
`restFiberIsometry`.  Nothing is reproved. -/
def quotientRestBridge (u : FrameHyperboloid) :
    FrameQuotient ≃ₗ[ℝ] RestSpace u :=
  (quotientEquivVec3.trans vec3LinearEquivEuclidean).trans (restFiberIsometry u)

theorem quotientRestBridge_apply (u : FrameHyperboloid) (s : FrameQuotient) :
    quotientRestBridge u s =
      restFiberIsometry u (vec3LinearEquivEuclidean (quotientEquivVec3 s)) :=
  rfl

/-- Exact metric identity of the bridge: the committed rest metric pulled
back through the bridge is the committed screen quotient Gram, at every
declared frame. -/
theorem quotientRestBridge_metric (u : FrameHyperboloid)
    (q r : FrameQuotient) :
    restMetric u (quotientRestBridge u q) (quotientRestBridge u r) =
      quotientGram q r := by
  rw [quotientRestBridge_apply, quotientRestBridge_apply,
    restFiberIsometry_metric]
  rfl

/-- At the standard frame the bridge is the committed candidate bridge of
`Geometry/SpatialReadbackSoldering.lean`: the commuting triangle with the
committed screen-to-rest identification is definitional up to the committed
boost-shear reduction. -/
theorem quotientRestBridge_standardFrame (q : FrameQuotient) :
    quotientRestBridge standardFrame q = frameQuotientEquivStandardRest q := by
  rw [quotientRestBridge_apply, restFiberIsometry_standardFrame]
  rfl

/-- The declared radial readout of the rest fiber reads the screen Gram
exactly through the bridge. -/
theorem quotientRestBridge_radius (u : FrameHyperboloid) (q : FrameQuotient) :
    restRadius u (quotientRestBridge u q) ^ 2 = quotientGram q q := by
  rw [restRadius_sq, quotientRestBridge_metric]

/-- The bridge carries the quotient class of every integer record control to
the committed screen chain image, scaled by the committed unit-diagonal
normalization radius: the commuting triangle with
`Geometry/RestFiberShellTransport.lean`'s `sourceChain`. -/
theorem quotientRestBridge_sourceChain (u : FrameHyperboloid)
    (p : IntegerFramePoint) :
    quotientRestBridge u
        (Submodule.Quotient.mk (castIntegerControl p.control)) =
      rawRadius • sourceChain u p := by
  rw [quotientRestBridge_apply]
  have hmk : quotientEquivVec3
      ((Submodule.Quotient.mk (castIntegerControl p.control) :
        FrameQuotient)) = pointFrame p :=
    quotientEquivVec3_mk (castIntegerControl p.control)
  rw [hmk]
  unfold sourceChain
  rw [← map_smul]
  congr 1
  unfold pointEuclideanFrame
  rw [smul_inv_smul₀ rawRadius_ne_zero]
  rfl

/-! ## The extended witness record -/

/-- The bridged common-world architecture: the committed two-island record
extended by one declared cross-island dictionary from the screen Gram
quotient onto the rest fiber of the record's own frame field, pinned by the
metric clause.  Every base field stays a declared hypothesis exactly as in
the committed witness. -/
structure BridgedCommonWorldArchitecture (n k : ℕ) extends
    CommonWorldArchitecture n k where
  /-- The declared cross-island dictionary. -/
  bridge : FrameQuotient ≃ₗ[ℝ] RestSpace toCommonWorldArchitecture.frame
  /-- The bridge clause: the rest metric pulled back through the dictionary
  is the committed screen quotient Gram. -/
  bridge_metric : ∀ q r : FrameQuotient,
    restMetric toCommonWorldArchitecture.frame (bridge q) (bridge r) =
      quotientGram q r

/-- The island-bridge clause of the composed receipt: metric equality, zero
preservation both ways, dimension three on both sides, the radius readout,
and the commuting square with every committed oriented Lorentz chart
transport. -/
def CarriesIslandBridge (W : BridgedCommonWorldArchitecture n k) : Prop :=
  (∀ q r : FrameQuotient,
    restMetric W.frame (W.bridge q) (W.bridge r) = quotientGram q r) ∧
  (∀ q : FrameQuotient, W.bridge q = 0 ↔ q = 0) ∧
  (finrank ℝ FrameQuotient = 3 ∧ finrank ℝ (RestSpace W.frame) = 3) ∧
  (∀ q : FrameQuotient,
    restRadius W.frame (W.bridge q) ^ 2 = quotientGram q q) ∧
  (∀ (L : OrientedLorentzEquiv) (q r : FrameQuotient),
    restMetric (L.mapFrame W.frame)
        (L.restEquiv W.frame (W.bridge q))
        (L.restEquiv W.frame (W.bridge r)) =
      quotientGram q r)

theorem carriesIslandBridge (W : BridgedCommonWorldArchitecture n k) :
    CarriesIslandBridge W :=
  ⟨W.bridge_metric,
    fun _ => W.bridge.map_eq_zero_iff,
    ⟨frameQuotient_finrank, finrank_restSpace W.frame⟩,
    fun q => by rw [restRadius_sq, W.bridge_metric],
    fun L q r => by
      rw [L.restEquiv_preserves_metric, W.bridge_metric]⟩

/-- The bridged shell readout: the committed transported inverse-square law
of `Geometry/RestFiberShellTransport.lean` read on the screen Gram through
the bridge.  The flux packet carries the declared premises PR-18 and
PR-19. -/
def CarriesBridgedShellReadout (W : BridgedCommonWorldArchitecture n k)
    (D : RadialFluxData) : Prop :=
  ∀ q : FrameQuotient, q ≠ 0 →
    restField W.frame D (W.bridge q) =
      D.charge / (4 * π * quotientGram q q)

theorem carriesBridgedShellReadout (W : BridgedCommonWorldArchitecture n k)
    (D : RadialFluxData) : CarriesBridgedShellReadout W D := by
  intro q hq
  have hv : W.bridge q ≠ 0 := fun h => hq (W.bridge.map_eq_zero_iff.mp h)
  rw [restField_inverse_square_metric W.frame D hv, W.bridge_metric]

/-! ## The composed receipt on one antecedent bundle -/

/-- The typed antecedent bundle of the bridged witness: one extended record
and one declared radial flux packet (PR-18 and PR-19, consumed exactly as
in the committed shell-transport claim). -/
structure BridgedWorldAntecedents (n k : ℕ) where
  /-- The extended common-world record. -/
  world : BridgedCommonWorldArchitecture n k
  /-- PR-18 (radial readout) and PR-19 (shell-flux normalization), as the
  committed `RadialFluxData` premise packet. -/
  flux : RadialFluxData

/-- **The bridged common-world receipt (issue #740, one added join; CW1
stays open).**  Every antecedent bundle simultaneously carries the eight
committed island packets of the base witness on its one record, the
island-bridge clause tying the screen Gram quotient to the rest fiber of
the record's one frame, and the bridged shell readout.  The two islands
share the rest-fiber object up to the bridge isometry; the joins that
remain missing (shared clock, common action, common dynamics, carrier
dictionary, physical attachment) are listed in the file header and stay on
the open register rows PR-15, PR-52, PR-53, and PR-54. -/
theorem bridgedCommonWorld_receipt (A : BridgedWorldAntecedents n k) :
    CarriesMassShell A.world.toCommonWorldArchitecture ∧
      CarriesInternalClock A.world.toCommonWorldArchitecture ∧
      CarriesFreeEvolution A.world.toCommonWorldArchitecture ∧
      CarriesCausalOrder A.world.toCommonWorldArchitecture ∧
      CarriesKinematicJoins A.world.toCommonWorldArchitecture ∧
      CarriesTemporalMaxwell A.world.toCommonWorldArchitecture ∧
      CarriesElectroweakBreaking A.world.toCommonWorldArchitecture ∧
      CarriesScreenPortMeasure A.world.toCommonWorldArchitecture ∧
      CarriesIslandBridge A.world ∧
      CarriesBridgedShellReadout A.world A.flux := by
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8⟩ :=
    commonWorldArchitecture_receipt A.world.toCommonWorldArchitecture
  exact ⟨h1, h2, h3, h4, h5, h6, h7, h8,
    carriesIslandBridge A.world,
    carriesBridgedShellReadout A.world A.flux⟩

/-! ## The explicit inhabitant extending the committed inhabitant -/

/-- The explicit extended inhabitant: the committed witness of
`Geometry/CommonWorldKinematicsWitness.lean` extended by the canonical
bridge at its own declared frame.  The parameter values select nothing. -/
def bridgedCommittedWitness : BridgedCommonWorldArchitecture 2 1 where
  toCommonWorldArchitecture := committedWitness
  bridge := quotientRestBridge committedWitness.frame
  bridge_metric := quotientRestBridge_metric committedWitness.frame

/-- The extended inhabitant restricts to the committed inhabitant
definitionally. -/
theorem bridgedCommittedWitness_extends :
    bridgedCommittedWitness.toCommonWorldArchitecture = committedWitness :=
  rfl

/-- The extended record is inhabited. -/
theorem bridgedCommonWorld_inhabited :
    Nonempty (BridgedCommonWorldArchitecture 2 1) :=
  ⟨bridgedCommittedWitness⟩

/-- The inhabitant's bridge agrees with the committed candidate bridge of
`Geometry/SpatialReadbackSoldering.lean` at the inhabitant's standard
frame. -/
theorem bridgedCommittedWitness_bridge_standard (q : FrameQuotient) :
    bridgedCommittedWitness.bridge q = frameQuotientEquivStandardRest q :=
  quotientRestBridge_standardFrame q

/-- The inhabitant's bridge carries every integer record control class to
the committed screen chain image, scaled by the committed normalization
radius. -/
theorem bridgedCommittedWitness_sourceChain (p : IntegerFramePoint) :
    bridgedCommittedWitness.bridge
        (Submodule.Quotient.mk (castIntegerControl p.control)) =
      rawRadius • sourceChain committedWitness.frame p :=
  quotientRestBridge_sourceChain committedWitness.frame p

/-- One explicit demonstration instance of the flux premise packet: unit
charge with the matching radial profile.  The values select nothing. -/
def unitFlux : RadialFluxData where
  charge := 1
  strength := fun r => 1 / (4 * π * r ^ 2)
  normalized := fun r hr => by
    rw [shellContent_eq hr, one_div,
      inv_mul_cancel₀ (by positivity : 4 * π * r ^ 2 ≠ 0)]

/-- The explicit antecedent bundle of the composed receipt. -/
def committedBridgedAntecedents : BridgedWorldAntecedents 2 1 :=
  ⟨bridgedCommittedWitness, unitFlux⟩

/-- The composed receipt holds of the one explicit inhabitant and the one
demonstration flux packet. -/
theorem bridgedCommittedWitness_receipt :
    CarriesMassShell bridgedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesInternalClock bridgedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesFreeEvolution bridgedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesCausalOrder bridgedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesKinematicJoins bridgedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesTemporalMaxwell
        bridgedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesElectroweakBreaking
        bridgedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesScreenPortMeasure
        bridgedCommittedWitness.toCommonWorldArchitecture ∧
      CarriesIslandBridge bridgedCommittedWitness ∧
      CarriesBridgedShellReadout bridgedCommittedWitness unitFlux :=
  bridgedCommonWorld_receipt committedBridgedAntecedents

/-! ## The load-bearing receipt: dropping the bridge clause readmits a
mismatched pair -/

/-- Scaling every rest-fiber vector by two, as a linear equivalence. -/
def restDoubling (u : FrameHyperboloid) : RestSpace u ≃ₗ[ℝ] RestSpace u where
  toFun v := (2 : ℝ) • v
  invFun v := (2 : ℝ)⁻¹ • v
  left_inv v := inv_smul_smul₀ two_ne_zero v
  right_inv v := smul_inv_smul₀ two_ne_zero v
  map_add' v w := smul_add (2 : ℝ) v w
  map_smul' a v := smul_comm (2 : ℝ) a v

/-- The mismatched dictionary: the canonical bridge composed with scaling by
two.  A perfectly typed linear equivalence between the same two objects. -/
def doubledBridge (u : FrameHyperboloid) :
    FrameQuotient ≃ₗ[ℝ] RestSpace u :=
  (quotientRestBridge u).trans (restDoubling u)

/-- One explicit quotient class of unit Gram square. -/
def unitQuotient : FrameQuotient := quotientEquivVec3.symm ![1, 0, 0]

theorem unitQuotient_gram : quotientGram unitQuotient unitQuotient = 1 := by
  unfold quotientGram unitQuotient
  rw [LinearEquiv.apply_symm_apply]
  unfold OPH.PrimitivePortTranslationBridge.dot
  simp [Fin.sum_univ_three]

/-- The mismatched dictionary violates the metric identity at the unit
class, at every declared frame: its rest-metric pullback there is four
while the screen Gram is one. -/
theorem doubledBridge_breaks_metric (u : FrameHyperboloid) :
    restMetric u (doubledBridge u unitQuotient)
        (doubledBridge u unitQuotient) ≠
      quotientGram unitQuotient unitQuotient := by
  have happ : doubledBridge u unitQuotient =
      (2 : ℝ) • quotientRestBridge u unitQuotient := rfl
  rw [happ, restMetric_smul_left, restMetric_smul_right,
    quotientRestBridge_metric, unitQuotient_gram]
  norm_num

/-- The same extension with the metric clause dropped. -/
structure BareBridgedArchitecture (n k : ℕ) extends
    CommonWorldArchitecture n k where
  /-- A declared cross-island dictionary with no metric clause. -/
  bridge : FrameQuotient ≃ₗ[ℝ] RestSpace toCommonWorldArchitecture.frame

/-- The readmitted mismatched pair: the bare extension of the same committed
base witness by the doubled dictionary. -/
def mismatchedExtension : BareBridgedArchitecture 2 1 where
  toCommonWorldArchitecture := committedWitness
  bridge := doubledBridge committedWitness.frame

/-- The mismatched extension's dictionary violates the dropped clause at an
exhibited class: the bare record type admits a pair the extended record
type excludes. -/
theorem mismatchedExtension_breaks_metric :
    ∃ q : FrameQuotient,
      restMetric mismatchedExtension.frame
          (mismatchedExtension.bridge q) (mismatchedExtension.bridge q) ≠
        quotientGram q q :=
  ⟨unitQuotient, doubledBridge_breaks_metric committedWitness.frame⟩

/-- The retained metric clause is load-bearing: no extended record carries
the doubled dictionary. -/
theorem bridge_metric_excludes_doubling
    (W : BridgedCommonWorldArchitecture n k)
    (h : ∀ q : FrameQuotient, W.bridge q = doubledBridge W.frame q) :
    False := by
  have hm := W.bridge_metric unitQuotient unitQuotient
  rw [h unitQuotient] at hm
  exact doubledBridge_breaks_metric W.frame hm

end

end OPH.CommonWorldIslandBridge

/- Axiom audit: committed receipts and exact finite mathematics only.
Expected axioms per line: at most `propext`, `Classical.choice`,
`Quot.sound`.  No native decision procedure is used. -/

#print axioms OPH.CommonWorldIslandBridge.quotientRestBridge_metric
#print axioms OPH.CommonWorldIslandBridge.quotientRestBridge_standardFrame
#print axioms OPH.CommonWorldIslandBridge.quotientRestBridge_radius
#print axioms OPH.CommonWorldIslandBridge.quotientRestBridge_sourceChain
#print axioms OPH.CommonWorldIslandBridge.carriesIslandBridge
#print axioms OPH.CommonWorldIslandBridge.carriesBridgedShellReadout
#print axioms OPH.CommonWorldIslandBridge.bridgedCommonWorld_receipt
#print axioms OPH.CommonWorldIslandBridge.bridgedCommittedWitness_extends
#print axioms OPH.CommonWorldIslandBridge.bridgedCommonWorld_inhabited
#print axioms OPH.CommonWorldIslandBridge.bridgedCommittedWitness_bridge_standard
#print axioms OPH.CommonWorldIslandBridge.bridgedCommittedWitness_sourceChain
#print axioms OPH.CommonWorldIslandBridge.bridgedCommittedWitness_receipt
#print axioms OPH.CommonWorldIslandBridge.unitQuotient_gram
#print axioms OPH.CommonWorldIslandBridge.doubledBridge_breaks_metric
#print axioms OPH.CommonWorldIslandBridge.mismatchedExtension_breaks_metric
#print axioms OPH.CommonWorldIslandBridge.bridge_metric_excludes_doubling
