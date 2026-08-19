import Geometry.FreeEvolutionPersistence
import Geometry.CausalOrderComposition
import QFT.ObserverAccessCut
import TemporalMaxwellEvolution
import ElectroweakBreakingComposition

set_option autoImplicit false

open scoped BigOperators

/-!
# CW1 partial witness: one inhabited two-island common-world architecture (issue #740)

WHAT IS PROVED.  One typed record, `CommonWorldArchitecture`, bundles the
shared objects of two clusters of committed lane receipts.  The kinematics
island consists of one declared positive mass parameter, one declared state
frame, and one declared flow frame on the committed Hermitian Lorentz module
`Herm2` of `Geometry/CanonicalLorentzModule.lean`, together with one declared
projective partition and state for the committed partition causal net.  The
screen island consists of one declared Maxwell evolution bundle on the
committed twelve-port, thirty-seam, twenty-face carrier and one declared
electroweak-breaking premise bundle whose base instantiates the committed
Standard-Model structure premises, including the port-dual weighting of the
same twelve ports.

For every inhabitant of the record, four kinematics receipts hold of the one
mass, frame, flow frame, partition, and state simultaneously: the kinematic
mass-shell packet of `Geometry/MassShellKinematics.lean`, the internal-clock
packet of `Geometry/InternalClockRestFrequency.lean`, the free-evolution
persistence packet of `Geometry/FreeEvolutionPersistence.lean`, and the
causal-order composition of `Geometry/CausalOrderComposition.lean`.  The
shared-object identifications are proved by `rfl`: the clock packet and the
persistence packet read the same mass field and the same frame field, and
the four-momentum of the clock packet is the four-momentum of the
persistence packet.  Three join theorems tie the packets together on the one
module: the record's four-momentum is future-causal from the origin and from
the committed causal net's bottom-region image, the clock worldline is
causally monotone in its proper-time parameter, and every nonnegative-time
free-evolution displacement is future-causal, all in the cone order of the
same committed module that carries the finite causal nets.

On the screen island, the temporal Maxwell receipt of
`Screen/TemporalMaxwellEvolution.lean` and the electroweak-breaking receipt
of `Screen/ElectroweakBreakingComposition.lean` hold of the one record
simultaneously, and one join theorem evaluates the pinned port-dual measure
of the record's Standard-Model structure premises against the record's
Maxwell Gauss load on the one committed twelve-port index type: every port
weight is `1/12` and the measure-weighted total load is one twelfth of the
plain total load, with no transport map between the two packages.

One explicit inhabitant `committedWitness` instantiates the record: unit
declared mass at the standard frame, the committed one-block partition and
diagonal witness state, the committed nonstatic Maxwell demonstration
bundle, and the committed electroweak premise instance.  The composed
receipt `commonWorldArchitecture_receipt` holds of every inhabitant, and
`committedWitness_receipt` states it of the one explicit inhabitant, so
every listed lane receipt holds of one record simultaneously.

DECLARED DATA.  Every record field is a declared hypothesis: the mass,
frames, partition, state, Maxwell histories and sources, and electroweak
parameters are supplied, never derived.  The two vertex assignments of the
causal composition stay declared witness constructions of the committed
module.

FALSIFIER.  The witness fails if any listed committed receipt fails at the
record's fields, if a claimed `rfl` identification requires a transport map,
if a join theorem fails (a positive-mass four-momentum outside the future
cone, a causally decreasing worldline, a non-causal nonnegative flow
displacement, or a port-dual weight differing from `1/12` on some port of
the record's base premises), or if the record type is uninhabited.

NONCLAIMS AND BOUNDARY (the missing joins; CW1 stays open).  This is a
partial CW1 inhabitant covering the kinematics cluster and the screen
cluster as two internally joined islands.  The witness joins exactly six
lane receipts: mass shell, internal clock, free-evolution persistence, and
causal order on the kinematics island, and temporal Maxwell evolution and
electroweak breaking on the screen island.  Every other lane is absent
from the record: gravitation and the Einstein-tensor bridge,
thermodynamics and the repair generator, matter dynamics beyond the
committed structure premises, quantization and Born readout, and the
electromagnetic observer readout appear in no field, no packet, and no
join.  The two islands are joined to
each other only by conjunction on one record.  No typed dictionary joins
the screen carrier to the Lorentz module: no map sends the twelve ports,
thirty seams, or twenty faces to points, intervals, or cones of the
declared module, and no committed theorem constrains such a map (PR-52 and
PR-53 open).  No shared clock exists: the Maxwell step index is a declared
natural-number evolution parameter, the kinematics proper time is a real
coordinate on the declared module, and no map identifies them or attaches
either to a physical clock (PR-15 open).  No common action exists: the
assembled screen action and the declared kinematic flow share no
variational principle on one carrier (PR-54 open).  No observer or
public-record surface produces or reads either island; the partition and
state fields type the committed partition net, whose causal clauses do not
consume the physical data.  The Maxwell face index type and the
reconstructed face triples of the port-dual measure are distinct committed
types with no committed dictionary; the twelve-port evaluation is the only
cross-package measured overlap proved on the screen island.  Nothing here
closes issue #740, and no observation-ledger row is promoted.

Axiom audit.  Every proof composes committed receipts with exact finite
mathematics; the module adds no project axiom and uses no native decision
procedure.  The guard lines at the end of the file show at most `propext`,
`Classical.choice`, and `Quot.sound`.
-/

namespace OPH.CommonWorld

open OPH.C1Lorentz OPH.C2Soldering OPH.CausalComposition OPH.QFT
open OPH.TemporalMaxwellEvolution OPH.ElectroweakBreakingComposition
open OPH.DiscreteCoulombGreen OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.PrimitivePortDualMeasure (portDualWeight)

noncomputable section

variable {n k : ℕ}

/-! ## The one typed architecture record -/

/-- The common-world architecture record: one simultaneously inhabitable
bundle of the shared objects of the kinematics island and the screen
island.  Every field is a declared hypothesis.  The kinematics fields are
typed on the one committed Hermitian Lorentz module `Herm2`; the screen
fields are typed on the one committed twelve-port carrier and its committed
coarse grammar. -/
structure CommonWorldArchitecture (n k : ℕ) where
  /-- The one declared positive mass parameter, shared by the mass-shell,
  internal-clock, and free-evolution packets. -/
  mass : ℝ
  /-- Declared positivity of the mass parameter. -/
  mass_pos : 0 < mass
  /-- The one declared state frame on the future unit-timelike hyperboloid
  of the committed module, shared by the mass-shell, internal-clock, and
  free-evolution packets. -/
  frame : FrameHyperboloid
  /-- The declared frame whose worldline directs the free-evolution flow. -/
  flowFrame : FrameHyperboloid
  /-- The declared projective partition of the committed partition causal
  net. -/
  partition : EventAlgebra.ProjectivePartition n k
  /-- The declared state matrix of the committed partition causal net. -/
  state : EventAlgebra.StateMatrix n
  /-- The declared screen data on the committed twelve-port, thirty-seam,
  twenty-face complex: seam and port histories, sources, the Ampere
  evolution law, the initial Gauss constraint, and the continuity
  equation. -/
  maxwell : MaxwellEvolutionBundle
  /-- The declared electroweak parameter fields over a declared instance of
  the committed Standard-Model structure premises, including the port-dual
  weighting of the same committed twelve ports. -/
  ew : EWBreakingPremiseData

/-! ## The derived lane packets and their `rfl` identifications -/

/-- The internal-clock antecedents read off the record: no second mass or
frame is introduced. -/
def clockPacket (W : CommonWorldArchitecture n k) : InternalClockAntecedents :=
  ⟨W.mass, W.mass_pos, W.frame⟩

/-- The free-evolution antecedents read off the record: no second mass,
frame, or flow frame is introduced. -/
def persistencePacket (W : CommonWorldArchitecture n k) :
    FreeEvolutionAntecedents :=
  ⟨W.mass, W.mass_pos, W.frame, W.flowFrame⟩

/-- The causal-composition antecedents read off the record. -/
def causalPacket (W : CommonWorldArchitecture n k) :
    CausalCompositionAntecedents n k :=
  ⟨W.partition, W.state⟩

/-- The clock packet and the persistence packet carry the record's one mass
field: both identifications are definitional. -/
theorem sharedMass (W : CommonWorldArchitecture n k) :
    (clockPacket W).mass = W.mass ∧
      (persistencePacket W).mass = (clockPacket W).mass :=
  ⟨rfl, rfl⟩

/-- The clock packet and the persistence packet carry the record's one
frame field: both identifications are definitional. -/
theorem sharedFrame (W : CommonWorldArchitecture n k) :
    (clockPacket W).frame = W.frame ∧
      (persistencePacket W).frame = (clockPacket W).frame :=
  ⟨rfl, rfl⟩

/-- The clock packet and the persistence packet parameterize the same
four-momentum vector of the one committed module: definitional. -/
theorem sharedMomentum (W : CommonWorldArchitecture n k) :
    fourMomentum (clockPacket W).mass (clockPacket W).frame =
      fourMomentum (persistencePacket W).mass (persistencePacket W).frame :=
  rfl

/-- The internal clock of the clock packet is the free-evolution flow phase
of the persistence packet along its own frame: the committed clock identity
instantiated at the record's shared fields. -/
theorem sharedClock (W : CommonWorldArchitecture n k) (τ : ℝ) :
    internalClock (clockPacket W).mass (clockPacket W).frame τ =
      flowPhase (persistencePacket W).frame
        (fourMomentum (persistencePacket W).mass
          (persistencePacket W).frame) τ :=
  (flowPhase_self_eq_internalClock W.mass W.frame τ).symm

/-! ## The kinematics island: four packets on one module instance -/

/-- The kinematic mass-shell packet at the record's one mass and frame. -/
def CarriesMassShell (W : CommonWorldArchitecture n k) : Prop :=
  lorentzQ (fourMomentum W.mass W.frame) = W.mass ^ 2
  ∧ (fourMomentum W.mass W.frame).1 ^ 2 =
      W.mass ^ 2 + spatialNormSq (fourMomentum W.mass W.frame).2
  ∧ W.mass ≤ (fourMomentum W.mass W.frame).1
  ∧ 0 < (fourMomentum W.mass W.frame).1
  ∧ ((fourMomentum W.mass W.frame).1 = W.mass ↔
      (fourMomentum W.mass W.frame).2 = 0)
  ∧ (∀ L : OrientedLorentzEquiv,
      lorentzQ (L (fourMomentum W.mass W.frame)) = W.mass ^ 2)
  ∧ (∀ L : OrientedLorentzEquiv,
      (L (fourMomentum W.mass W.frame)).1 ^ 2 =
        W.mass ^ 2 + spatialNormSq (L (fourMomentum W.mass W.frame)).2)
  ∧ (∀ v : Herm2, lorentzQ v = W.mass ^ 2 → 0 < v.1 →
      ∃ u : FrameHyperboloid, fourMomentum W.mass u = v)

theorem carriesMassShell (W : CommonWorldArchitecture n k) :
    CarriesMassShell W :=
  ⟨lorentzQ_fourMomentum W.mass W.frame,
    fourMomentum_time_sq_eq_mass_sq_add_spatial W.mass W.frame,
    mass_le_fourMomentum_time W.mass_pos.le W.frame,
    fourMomentum_time_pos W.mass_pos W.frame,
    fourMomentum_time_eq_mass_iff_spatial_zero W.mass_pos W.frame,
    fun L => lorentzQ_oriented_fourMomentum L W.mass W.frame,
    fun L => fourMomentum_time_sq_oriented L W.mass W.frame,
    fun _ hshell hpos => exists_frame_fourMomentum_eq W.mass_pos hshell hpos⟩

/-- The internal-clock packet at the record's one mass and frame: the body
of the committed composed receipt of
`Geometry/InternalClockRestFrequency.lean`. -/
def CarriesInternalClock (W : CommonWorldArchitecture n k) : Prop :=
  (∀ τ : ℝ,
    planeWavePhase (fourMomentum W.mass W.frame)
      (frameWorldline W.frame τ) = W.mass * τ) ∧
  (∀ τ : ℝ,
    HasDerivAt
      (fun σ : ℝ =>
        planeWavePhase (fourMomentum W.mass W.frame)
          (frameWorldline W.frame σ)) W.mass τ) ∧
  (∀ L : OrientedLorentzEquiv,
    internalClock W.mass (L.mapFrame W.frame) =
      internalClock W.mass W.frame) ∧
  (∀ T : ℝ,
    Function.Periodic (internalClock W.mass W.frame) T ↔
      ∃ m : ℤ, T = m * (2 * Real.pi / W.mass)) ∧
  (∀ (t : ℝ) (xs : Spatial) (i : Fin 3) (s : ℝ),
    HasDerivAt
      (fun σ : ℝ =>
        wavePhase (fourMomentum W.mass W.frame)
          (t, xs + σ • (Pi.single i 1 : Spatial)))
      ((fourMomentum W.mass W.frame).2 i) s) ∧
  (∀ (xs : Spatial) (t : ℝ),
    HasDerivAt
      (fun σ : ℝ => planeWavePhase (fourMomentum W.mass W.frame) (σ, xs))
      (fourMomentum W.mass W.frame).1 t) ∧
  (fourMomentum W.mass W.frame).1 ^ 2 =
    W.mass ^ 2 + spatialNormSq (fourMomentum W.mass W.frame).2

theorem carriesInternalClock (W : CommonWorldArchitecture n k) :
    CarriesInternalClock W :=
  internalClock_receipt (clockPacket W)

/-- The free-evolution persistence packet at the record's one mass, frame,
and flow frame: the body of the committed composed receipt of
`Geometry/FreeEvolutionPersistence.lean`. -/
def CarriesFreeEvolution (W : CommonWorldArchitecture n k) : Prop :=
  lorentzQ (fourMomentum W.mass W.frame) = W.mass ^ 2 ∧
  (∀ τ : ℝ,
    freeEvolution W.flowFrame τ (planeWave (fourMomentum W.mass W.frame)) =
      fun x => flowPhase W.flowFrame (fourMomentum W.mass W.frame) τ *
        planeWave (fourMomentum W.mass W.frame) x) ∧
  (∀ (τ t : ℝ) (xs : Spatial),
    HasDerivAt
      (fun σ : ℝ =>
        planeWavePhase (fourMomentum W.mass W.frame)
          ((σ, xs) + frameWorldline W.flowFrame τ))
      (fourMomentum W.mass W.frame).1 t) ∧
  (∀ (τ t : ℝ) (xs : Spatial) (i : Fin 3) (s : ℝ),
    HasDerivAt
      (fun σ : ℝ =>
        wavePhase (fourMomentum W.mass W.frame)
          ((t, xs + σ • (Pi.single i 1 : Spatial)) +
            frameWorldline W.flowFrame τ))
      ((fourMomentum W.mass W.frame).2 i) s) ∧
  (∀ (τ : ℝ) (x : Herm2),
    ‖freeEvolution W.flowFrame τ
      (planeWave (fourMomentum W.mass W.frame)) x‖ = 1) ∧
  (∀ m' : ℝ, 0 < m' → ∀ (f' : FrameHyperboloid) (τ : ℝ),
    freeEvolution W.flowFrame τ (planeWave (fourMomentum W.mass W.frame)) =
      freeEvolution W.flowFrame τ (planeWave (fourMomentum m' f')) →
    W.mass = m' ∧ W.frame = f') ∧
  (∀ ψ : Herm2 → ℂ, freeEvolution W.flowFrame 0 ψ = ψ) ∧
  (∀ (τ σ : ℝ) (ψ : Herm2 → ℂ),
    freeEvolution W.flowFrame σ (freeEvolution W.flowFrame τ ψ) =
      freeEvolution W.flowFrame (τ + σ) ψ) ∧
  (∀ τ : ℝ,
    flowPhase W.frame (fourMomentum W.mass W.frame) τ =
      internalClock W.mass W.frame τ)

theorem carriesFreeEvolution (W : CommonWorldArchitecture n k) :
    CarriesFreeEvolution W :=
  freeEvolution_receipt (persistencePacket W)

/-- The causal-order composition at the record's one partition and state,
on the same committed module: the body of the committed composed receipt of
`Geometry/CausalOrderComposition.lean`. -/
def CarriesCausalOrder (W : CommonWorldArchitecture n k) : Prop :=
  ((∀ v, causalLE v v) ∧
    (∀ u v w, causalLE u v → causalLE v w → causalLE u w) ∧
    (∀ u v, causalLE u v → causalLE v u → u = v)) ∧
  (∀ U V : RichRegion,
    richFibreNet.regionLE () U V ↔
      causalLE (richVertex U) (richVertex V)) ∧
  Function.Injective richVertex ∧
  ((IsTimelikeSep (richVertex .bot) (richVertex .r86) ∧
      causalLE (richVertex .bot) (richVertex .r86)) ∧
    (IsTimelikeSep (richVertex .r88) (richVertex .top) ∧
      causalLE (richVertex .r88) (richVertex .top)) ∧
    (IsSpacelikeSep (richVertex .r86) (richVertex .r88) ∧
      ¬ causalLE (richVertex .r86) (richVertex .r88) ∧
      ¬ causalLE (richVertex .r88) (richVertex .r86)) ∧
    (IsSpacelikeSep (richVertex .r247) (richVertex .r384) ∧
      ¬ causalLE (richVertex .r247) (richVertex .r384) ∧
      ¬ causalLE (richVertex .r384) (richVertex .r247))) ∧
  (∀ U V : RichRegion, richFibreNet.disjoint () U V →
    IsSpacelikeSep (richVertex U) (richVertex V)) ∧
  (∀ (L : OrientedLorentzEquiv) (U V : RichRegion),
    richFibreNet.regionLE () U V ↔
      causalLE (L (richVertex U)) (L (richVertex V))) ∧
  (∀ (Chart : Type*) (T : LorentzOverlapCocycle Chart) (i j : Chart)
      (U V : RichRegion),
    richFibreNet.regionLE () U V ↔
      causalLE (T.act i j (richVertex U)) (T.act i j (richVertex V))) ∧
  (∀ U V : Finset (Fin 2),
    (FiniteCausalObserverNet.partitionPublicCausalNet
        W.partition W.state).regionLE () U V ↔
      causalLE (richVertex (siteRegionLabel U))
        (richVertex (siteRegionLabel V))) ∧
  (∀ U V : Finset (Fin 2),
    (FiniteCausalObserverNet.partitionPublicCausalNet
        W.partition W.state).disjoint () U V →
      IsSpacelikeSep (richVertex (siteRegionLabel U))
        (richVertex (siteRegionLabel V))) ∧
  (∃ f : RichRegion → Herm2, f ≠ richVertex ∧
    ∀ U V : RichRegion,
      richFibreNet.regionLE () U V ↔ causalLE (f U) (f V))

theorem carriesCausalOrder (W : CommonWorldArchitecture n k) :
    CarriesCausalOrder W :=
  causalOrderComposition_receipt (causalPacket W)

/-! ## The kinematics joins: the packets meet in one cone order -/

/-- The record's four-momentum is future-causal from the origin of the one
committed module: the mass-shell packet lands inside the cone order that
carries the committed causal nets. -/
theorem fourMomentum_causal_from_origin (W : CommonWorldArchitecture n k) :
    causalLE (0 : Herm2) (fourMomentum W.mass W.frame) := by
  show IsFutureCausal (fourMomentum W.mass W.frame - 0)
  rw [sub_zero]
  exact ⟨by rw [lorentzQ_fourMomentum]; exact sq_nonneg W.mass,
    (fourMomentum_time_pos W.mass_pos W.frame).le⟩

/-- The committed causal net's bottom-region image causally precedes the
record's four-momentum: the finite net and the mass-shell packet meet in
the same module. -/
theorem richBot_causal_fourMomentum (W : CommonWorldArchitecture n k) :
    causalLE (richVertex RichRegion.bot) (fourMomentum W.mass W.frame) := by
  have h : richVertex RichRegion.bot = (0 : Herm2) := rfl
  rw [h]
  exact fourMomentum_causal_from_origin W

/-- The internal-clock worldline of the record's frame is causally monotone
in its proper-time parameter. -/
theorem frameWorldline_causal_mono (W : CommonWorldArchitecture n k)
    (τ σ : ℝ) (h : τ ≤ σ) :
    causalLE (frameWorldline W.frame τ) (frameWorldline W.frame σ) := by
  show IsFutureCausal (frameWorldline W.frame σ - frameWorldline W.frame τ)
  have hdiff : frameWorldline W.frame σ - frameWorldline W.frame τ =
      (σ - τ) • (W.frame : Herm2) := by
    show σ • (W.frame : Herm2) - τ • (W.frame : Herm2) =
      (σ - τ) • (W.frame : Herm2)
    rw [sub_smul]
  rw [hdiff]
  refine ⟨?_, ?_⟩
  · rw [lorentzQ_smul, W.frame.2.1, mul_one]
    exact sq_nonneg _
  · have ht : ((σ - τ) • (W.frame : Herm2)).1 =
        (σ - τ) * (W.frame : Herm2).1 := rfl
    rw [ht]
    exact mul_nonneg (by linarith) W.frame.2.2.le

/-- Every nonnegative-parameter displacement of the record's free-evolution
flow is future-causal: the persistence flow moves every event of the one
module into its own causal future. -/
theorem freeEvolution_flow_causal (W : CommonWorldArchitecture n k)
    (x : Herm2) (τ : ℝ) (h : 0 ≤ τ) :
    causalLE x (x + frameWorldline W.flowFrame τ) := by
  show IsFutureCausal (x + frameWorldline W.flowFrame τ - x)
  have hcancel : x + frameWorldline W.flowFrame τ - x =
      frameWorldline W.flowFrame τ := by abel
  rw [hcancel]
  refine ⟨?_, ?_⟩
  · show 0 ≤ lorentzQ (τ • (W.flowFrame : Herm2))
    rw [lorentzQ_smul, W.flowFrame.2.1, mul_one]
    exact sq_nonneg τ
  · show 0 ≤ (τ • (W.flowFrame : Herm2)).1
    have ht : (τ • (W.flowFrame : Herm2)).1 =
        τ * (W.flowFrame : Herm2).1 := rfl
    rw [ht]
    exact mul_nonneg h W.flowFrame.2.2.le

/-- The kinematics joins packaged as one clause of the composed receipt. -/
def CarriesKinematicJoins (W : CommonWorldArchitecture n k) : Prop :=
  causalLE (0 : Herm2) (fourMomentum W.mass W.frame) ∧
  causalLE (richVertex RichRegion.bot) (fourMomentum W.mass W.frame) ∧
  (∀ τ σ : ℝ, τ ≤ σ →
    causalLE (frameWorldline W.frame τ) (frameWorldline W.frame σ)) ∧
  (∀ (x : Herm2) (τ : ℝ), 0 ≤ τ →
    causalLE x (x + frameWorldline W.flowFrame τ))

theorem carriesKinematicJoins (W : CommonWorldArchitecture n k) :
    CarriesKinematicJoins W :=
  ⟨fourMomentum_causal_from_origin W,
    richBot_causal_fourMomentum W,
    fun τ σ h => frameWorldline_causal_mono W τ σ h,
    fun x τ h => freeEvolution_flow_causal W x τ h⟩

/-! ## The screen island: two packages on one committed carrier -/

/-- The temporal Maxwell packet at the record's one screen bundle: the body
of the committed composed receipt of
`Screen/TemporalMaxwellEvolution.lean`. -/
def CarriesTemporalMaxwell (W : CommonWorldArchitecture n k) : Prop :=
  (∀ m, realBoundary (electricField W.maxwell.A W.maxwell.phi m) =
    W.maxwell.rho m)
  ∧ (∀ (ρ' : ℕ → Fin 12 → ℝ) (m : ℕ),
      realBoundary (electricField W.maxwell.A W.maxwell.phi m) = ρ' m →
      (realBoundary (electricField W.maxwell.A W.maxwell.phi (m + 1)) =
          ρ' (m + 1) ↔
        ρ' (m + 1) - ρ' m + realBoundary (W.maxwell.J m) = 0))
  ∧ (∀ m, magneticField W.maxwell.A (m + 1) - magneticField W.maxwell.A m =
      -faceCurvature (electricField W.maxwell.A W.maxwell.phi m))
  ∧ (∀ m, fieldEnergy W.maxwell.A W.maxwell.phi (m + 1) =
      fieldEnergy W.maxwell.A W.maxwell.phi m -
        (1 / 2) * realSeamInner
          (electricField W.maxwell.A W.maxwell.phi m +
            electricField W.maxwell.A W.maxwell.phi (m + 1))
          (W.maxwell.J m))
  ∧ (∀ χ : ℕ → Fin 12 → ℝ,
      (∀ m, electricField (gaugeTransformA W.maxwell.A χ)
          (gaugeTransformPhi W.maxwell.phi χ) m =
        electricField W.maxwell.A W.maxwell.phi m)
      ∧ (∀ m, magneticField (gaugeTransformA W.maxwell.A χ) m =
          magneticField W.maxwell.A m)
      ∧ AmpereEvolution (gaugeTransformA W.maxwell.A χ)
          (gaugeTransformPhi W.maxwell.phi χ) W.maxwell.J
      ∧ (∀ m, fieldEnergy (gaugeTransformA W.maxwell.A χ)
          (gaugeTransformPhi W.maxwell.phi χ) m =
        fieldEnergy W.maxwell.A W.maxwell.phi m))
  ∧ (W.maxwell.J = (fun _ ↦ 0) →
      (∀ m, fieldEnergy W.maxwell.A W.maxwell.phi m =
        fieldEnergy W.maxwell.A W.maxwell.phi 0)
      ∧ (∀ m, W.maxwell.A (m + 2) - (2 : ℝ) • W.maxwell.A (m + 1) +
          W.maxwell.A m + localMaxwellOperator (W.maxwell.A (m + 1)) +
          realCoboundary (W.maxwell.phi (m + 1) - W.maxwell.phi m) = 0))
  ∧ ((∀ m, W.maxwell.A m = W.maxwell.A 0) →
      (∀ m, W.maxwell.phi m = W.maxwell.phi 0) →
      W.maxwell.J = (fun _ ↦ 0) →
      ((∑ p : Fin 12, W.maxwell.rho 0 p) = 0
      ∧ (∀ m, magneticField W.maxwell.A m = 0)
      ∧ (∀ m, electricField W.maxwell.A W.maxwell.phi m =
          realCoulombField (W.maxwell.rho 0))))

theorem carriesTemporalMaxwell (W : CommonWorldArchitecture n k) :
    CarriesTemporalMaxwell W :=
  temporalMaxwellEvolution_receipt W.maxwell

/-- The electroweak-breaking packet at the record's one premise bundle: the
body of the committed composed receipt of
`Screen/ElectroweakBreakingComposition.lean`. -/
def CarriesElectroweakBreaking (W : CommonWorldArchitecture n k) : Prop :=
  ((∀ h : Doublet,
      (∀ h' : Doublet,
        potential W.ew.lam W.ew.vev h ≤ potential W.ew.lam W.ew.vev h')
        ↔ normSqSum h = W.ew.vev ^ 2)
    ∧ normSqSum (vacuum W.ew.vev) = W.ew.vev ^ 2
    ∧ (∀ h : Doublet, normSqSum h = W.ew.vev ^ 2 → h ≠ 0))
  ∧ (∀ A : EWParam,
      infAction W.ew.gW W.ew.gY A (vacuum W.ew.vev) = 0
        ↔ ∃! t : ℝ, A = ⟨0, 0, W.ew.gY * t, W.ew.gW * t⟩)
  ∧ ((∀ u : Doublet, 0 ≤ secondVariation W.ew.lam W.ew.vev u)
    ∧ (∀ (t : ℝ) (u : Doublet),
        potential W.ew.lam W.ew.vev (vacuum W.ew.vev + t • u)
          = secondVariation W.ew.lam W.ew.vev u * t ^ 2
            + 4 * W.ew.lam * W.ew.vev * (u 1).re * normSqSum u * t ^ 3
            + W.ew.lam * normSqSum u ^ 2 * t ^ 4)
    ∧ (∀ u : Doublet,
        secondVariation W.ew.lam W.ew.vev u = 0
          ↔ ∃ A : EWParam, infAction W.ew.gW W.ew.gY A (vacuum W.ew.vev) = u)
    ∧ (∀ u : Doublet,
        secondVariation W.ew.lam W.ew.vev u = 0
          ↔ ∃! c : ℝ × ℝ × ℝ, u = goldstoneForm c))
  ∧ ((∀ A : EWParam,
        gaugeMassForm W.ew.gW W.ew.gY W.ew.vev A = 0
          ↔ infAction W.ew.gW W.ew.gY A (vacuum W.ew.vev) = 0)
    ∧ gaugeMassForm W.ew.gW W.ew.gY W.ew.vev wDirection
        = W.ew.gW ^ 2 * W.ew.vev ^ 2 / 4
    ∧ gaugeMassForm W.ew.gW W.ew.gY W.ew.vev (zDirection W.ew.gW W.ew.gY)
          / paramNormSq (zDirection W.ew.gW W.ew.gY)
        = (W.ew.gW ^ 2 + W.ew.gY ^ 2) * W.ew.vev ^ 2 / 4
    ∧ gaugeMassForm W.ew.gW W.ew.gY W.ew.vev
        (photonDirection W.ew.gW W.ew.gY) = 0
    ∧ (gaugeMassForm W.ew.gW W.ew.gY W.ew.vev wDirection
          / paramNormSq wDirection)
        / (gaugeMassForm W.ew.gW W.ew.gY W.ew.vev
            (zDirection W.ew.gW W.ew.gY)
            / paramNormSq (zDirection W.ew.gW W.ew.gY))
        = W.ew.gW ^ 2 / (W.ew.gW ^ 2 + W.ew.gY ^ 2))
  ∧ ((∀ (ψ : Doublet) (χ : ℂ),
        yukawaLine W.ew.yuk (vacuum W.ew.vev) ψ χ
          = ((W.ew.yuk * W.ew.vev : ℝ) : ℂ) * (ψ 0 * χ))
    ∧ ((∀ (ψ : Doublet) (χ : ℂ),
          yukawaLine W.ew.yuk (vacuum W.ew.vev) ψ χ = 0)
        ↔ W.ew.yuk = 0 ∨ W.ew.vev = 0)
    ∧ ∃ dRow sRow : Fin 10,
        OPH.ExteriorSelection.mem W.ew.base.selectionMask.val dRow = true
          ∧ OPH.ExteriorSelection.mem W.ew.base.selectionMask.val sRow
              = true
          ∧ OPH.ExteriorSelection.isDoublet dRow = true
          ∧ OPH.ExteriorSelection.isDoublet sRow = false
          ∧ OPH.ExteriorSelection.isColored dRow = false
          ∧ OPH.ExteriorSelection.isColored sRow = false
          ∧ OPH.ExteriorSelection.odd dRow
              = OPH.ExteriorSelection.odd sRow
          ∧ OPH.BaryonDimensionSix.left
                (OPH.MatterGrammarIndexBridge.fieldIndex dRow)
              = OPH.BaryonDimensionSix.left
                (OPH.MatterGrammarIndexBridge.fieldIndex sRow)
          ∧ (OPH.ExteriorSelection.charge dRow
                + OPH.ExteriorSelection.charge sRow
                + scalarChargeSix = 0
              ∨ OPH.ExteriorSelection.charge dRow
                + OPH.ExteriorSelection.charge sRow
                - scalarChargeSix = 0))

theorem carriesElectroweakBreaking (W : CommonWorldArchitecture n k) :
    CarriesElectroweakBreaking W :=
  electroweakBreakingComposition_receipt W.ew

/-! ## The screen join: one port-dual measure evaluates one Gauss load -/

/-- Typed identity of the screen port carrier: the Maxwell Gauss load at
any step and the pinned port-dual weight of the record's Standard-Model
structure premises are functions on the one committed twelve-port type, so
their pointwise product is well typed with no transport map, and each
weighted load equals one twelfth of the plain load. -/
theorem sharedPortEvaluation (W : CommonWorldArchitecture n k) (m : ℕ)
    (p : Fin 12) :
    ((portDualWeight W.ew.base.vertexShare W.ew.base.faceWeight p : ℚ) : ℝ) *
        W.maxwell.rho m p
      = (1 / 12 : ℝ) * W.maxwell.rho m p := by
  rw [W.ew.base.port_dual_weight_pinned p]
  norm_num

/-- The screen join clause: every port weight of the record's base premises
is pinned to `1/12`, and the measure-weighted total Gauss load of the
record's Maxwell bundle is one twelfth of the plain total load, at every
step, on the one committed port carrier. -/
def CarriesScreenPortMeasure (W : CommonWorldArchitecture n k) : Prop :=
  (∀ p : Fin 12,
    portDualWeight W.ew.base.vertexShare W.ew.base.faceWeight p = 1 / 12) ∧
  (∀ m : ℕ,
    ∑ p : Fin 12,
      ((portDualWeight W.ew.base.vertexShare W.ew.base.faceWeight p : ℚ) :
          ℝ) * W.maxwell.rho m p
      = (1 / 12 : ℝ) * ∑ p : Fin 12, W.maxwell.rho m p)

theorem carriesScreenPortMeasure (W : CommonWorldArchitecture n k) :
    CarriesScreenPortMeasure W := by
  refine ⟨fun p => W.ew.base.port_dual_weight_pinned p, fun m => ?_⟩
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun p _ => sharedPortEvaluation W m p

/-! ## The composed receipt -/

/-- **The common-world two-island receipt (issue #740, partial CW1).**
Every inhabitant of the one architecture record simultaneously carries the
kinematic mass-shell packet, the internal-clock packet, the free-evolution
persistence packet, and the causal-order composition on its one mass,
frame, flow frame, partition, and state over the one committed Lorentz
module, joined by the cone-order theorems; and the temporal Maxwell packet
and the electroweak-breaking packet on its one screen bundle and premise
bundle, joined by the port-dual measure evaluation on the one committed
twelve-port carrier.  The two islands are joined to each other only by
conjunction on the record; the missing cross-island dictionaries are listed
in the file header and stay on the open register rows PR-15, PR-52, PR-53,
and PR-54.  The record carries no gravitation, thermodynamics,
matter-dynamics, quantization, or observer-readout lane. -/
theorem commonWorldArchitecture_receipt (W : CommonWorldArchitecture n k) :
    CarriesMassShell W ∧ CarriesInternalClock W ∧ CarriesFreeEvolution W ∧
      CarriesCausalOrder W ∧ CarriesKinematicJoins W ∧
      CarriesTemporalMaxwell W ∧ CarriesElectroweakBreaking W ∧
      CarriesScreenPortMeasure W :=
  ⟨carriesMassShell W, carriesInternalClock W, carriesFreeEvolution W,
    carriesCausalOrder W, carriesKinematicJoins W,
    carriesTemporalMaxwell W, carriesElectroweakBreaking W,
    carriesScreenPortMeasure W⟩

/-! ## The explicit inhabitant -/

/-- The explicit inhabitant of the architecture record: unit declared mass
at the standard frame with the standard frame directing the flow, the
committed one-block partition on the two-dimensional event algebra with the
committed diagonal witness state, the committed nonstatic Maxwell
demonstration bundle, and the committed electroweak premise instance.  The
parameter values select nothing. -/
def committedWitness : CommonWorldArchitecture 2 1 where
  mass := 1
  mass_pos := one_pos
  frame := standardFrame
  flowFrame := standardFrame
  partition := oneBlockPartition 2
  state := witnessState
  maxwell := demoBundle
  ew := committedEWBreakingPremiseData

/-- The architecture record is inhabited. -/
theorem commonWorld_inhabited : Nonempty (CommonWorldArchitecture 2 1) :=
  ⟨committedWitness⟩

/-- Every listed lane receipt holds of the one explicit inhabitant
simultaneously. -/
theorem committedWitness_receipt :
    CarriesMassShell committedWitness ∧
      CarriesInternalClock committedWitness ∧
      CarriesFreeEvolution committedWitness ∧
      CarriesCausalOrder committedWitness ∧
      CarriesKinematicJoins committedWitness ∧
      CarriesTemporalMaxwell committedWitness ∧
      CarriesElectroweakBreaking committedWitness ∧
      CarriesScreenPortMeasure committedWitness :=
  commonWorldArchitecture_receipt committedWitness

/-- The inhabitant's screen history is nonstatic: the committed
demonstration bundle moves on its first step, so the screen island of the
witness is not the trivial evolution. -/
theorem committedWitness_screen_nonstatic :
    committedWitness.maxwell.A 1 ≠ committedWitness.maxwell.A 0 :=
  demoBundle_nonvacuous.1

end

end OPH.CommonWorld

/- Axiom audit: committed receipts and exact finite mathematics only.
Expected axioms per line: at most `propext`, `Classical.choice`,
`Quot.sound`.  No native decision procedure is used. -/

#print axioms OPH.CommonWorld.sharedMass
#print axioms OPH.CommonWorld.sharedFrame
#print axioms OPH.CommonWorld.sharedMomentum
#print axioms OPH.CommonWorld.sharedClock
#print axioms OPH.CommonWorld.carriesMassShell
#print axioms OPH.CommonWorld.carriesInternalClock
#print axioms OPH.CommonWorld.carriesFreeEvolution
#print axioms OPH.CommonWorld.carriesCausalOrder
#print axioms OPH.CommonWorld.fourMomentum_causal_from_origin
#print axioms OPH.CommonWorld.richBot_causal_fourMomentum
#print axioms OPH.CommonWorld.frameWorldline_causal_mono
#print axioms OPH.CommonWorld.freeEvolution_flow_causal
#print axioms OPH.CommonWorld.carriesKinematicJoins
#print axioms OPH.CommonWorld.carriesTemporalMaxwell
#print axioms OPH.CommonWorld.carriesElectroweakBreaking
#print axioms OPH.CommonWorld.sharedPortEvaluation
#print axioms OPH.CommonWorld.carriesScreenPortMeasure
#print axioms OPH.CommonWorld.commonWorldArchitecture_receipt
#print axioms OPH.CommonWorld.commonWorld_inhabited
#print axioms OPH.CommonWorld.committedWitness_receipt
#print axioms OPH.CommonWorld.committedWitness_screen_nonstatic
