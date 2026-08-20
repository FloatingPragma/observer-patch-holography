import ObserverPatchHolography.EinsteinBranch.RegisterSurface
import HorizonThermalitySurface
import QFT.InhabitedTowerRow

set_option autoImplicit false

/-!
# An explicit inhabitant of the Einstein-branch register surface

## Claim

`registerOfHorizonRecord D H` is an explicit value of
`OPH.EinsteinBranch.EinsteinRegisterData ℕ (MassZero Ω) TorusPoint` for every
declared repair law `D : RepairLawData Ω B` (register row PR-07) and every
horizon record `H : HorizonRecordData D` (register row PR-42).  Its typed
modular flow supply is the finite thermodynamic discharge itself: the field
`modularFlowData` is definitionally
`thermoFirstLawData Ω D.ref H.boostCharge H.centralCharge`, so the variation
space of the composed Einstein branch at this value is the simplex tangent
space `MassZero Ω` of the thermodynamics lane, and the `finiteFirstLaw` clause
of `einstein_register_composition` at this value is the statement
`HorizonSurface.entanglement_first_law_first_order D H δρ` of the
observation-ledger row OL-B3.  `demoEinsteinRegister` evaluates the
constructor at one explicit two-state repair law with a nonconstant boost
charge.

The value exists, so the twenty-five premise clauses of the register are
jointly satisfiable and the composed implication has one instance on record.
The nonvacuity receipts exercise the clauses that this data can fail on
nonzero data: the bare tower has consistent and inconsistent states, the
readouts separate two consistent states, the entropy differential and the
bulk pairing are nonzero on an explicit mass-zero variation, the chart stress
field is nonconstant across the sixteen-point torus, the vacuum-referenced
constant is nonzero, the small-ball stress contraction is nonzero, and the
scaling family carries nonzero error channels.  `split_clause_load_bearing`
shows which clause the join depends on: with a uniform reference and zero
central charge in place of the split, the same boost charge fails the
modular-generator clause on the same variation.

Four conclusions of the composed bundle hold at this value by definition and
select nothing: the Einstein clause (the chart geometry is defined as
`κ · stress + Λ · η`), the small-ball arithmetic (the ball geometry
contraction is defined as `8π · G · stress`), the universal-source clause
(the entropy stress field is the stress field), and the vacuum reference.
The degeneration receipts `einstein_clause_definitional` and
`small_ball_clause_definitional` record the first two as `rfl`.  The
content of the inhabitant is the joint satisfiability of the premise bundle
and the typed join of the first-law clause to the thermodynamic variation
space; the Einstein clause at this value is a restatement of the chart
definition, a chart violating it is rejected by the null-tomography,
vacuum-reference, Ward, or Bianchi fields, and the conclusion itself filters
nothing.

## Declared inputs

Every input is a synthetic witness at the demo tier and is named here.

* Tower and readouts (PR-28): the committed `demoTower`; six readouts
  extending the committed E7 fragment `OPH.Tower.demoFragment` (event bit and
  geometry weight) by identity and decision arrows, with identity coarse
  maps.
* Audit candidate (PR-28): `demoCandidate D H`.  Evaluated fields that the
  data could fail: mismatch before and after one copy repair on `demoTower`,
  the split residual of the horizon record, the Ward and Bianchi defects of
  the torus chart at the base point, and the shrinking family.  Evaluated
  fields that vanish by definition: the kernel, area, and vacuum residues
  compare a quantity with its own defining expression, and the MaxEnt
  envelope defect is `2π - 2π` on the committed tangent.  The chart index
  dimension is `Fintype.card (Fin 4)`.  Declared fields: boundary-fibre
  cardinality, arrow digests, leak count, tomography rank, uncovered
  timelike directions, metric compatibility defect, component count, source
  digests, and the rescaling stabilizer.
* Modular flow (PR-23): `D.ref`, `H.boostCharge`, `H.centralCharge`, and the
  split `H.split`.
* MaxEnt tangent (PR-26): the committed `passingMaxEntTangent`.
* Small-ball reals (PR-26, PR-27): Newton constant `1`, radius `1`, stress
  contraction `1`, and the geometry contraction, bulk entropy variation, and
  area variation solved from the three small-ball identities.
* Scaling family (PR-27): `ℓ_r = 1/(r+1)` with error channels `ℓ_r^5`.
* Chart (PR-24, PR-25, PR-28): the four-direction binary torus
  `Fin 4 → ZMod 2` with unit-vector steps, the diagonal stress field whose
  `(i, i)` entry reads coordinate `i + 1`, geometry `κ · stress + Λ · η` with
  `κ = 8π`, Newton constant `1`, and `Λ = 1`.

## Not proved

No physical metric, curvature, stress tensor, horizon, temperature, Newton
constant, spacetime, or causal structure is attached.  The repair law and the
horizon record are arguments; the demo instance is a finite faithful law with
an injective visible map and is not the source-generated pair (the B12 and
B20 no-gos stand).  The torus chart, the readouts, and the ball reals are
chosen witnesses.  The result is a consistency certificate for the
register surface and a typed join between the thermodynamic discharge and the
composed branch; it is not evidence that any physical tower satisfies the
premises, and the observer-to-spacetime attachment (PR-52) is untouched.
The common-readout clause of the composed conclusion is taken at one
consistent state paired with itself, the only pairing the boundary fibre of
`demoTower` admits.  Five of the six readout channels are copies of one
protected bit, and `QFT/InhabitedTowerRow.lean` records that closing its
readout rows by such relabeling carries no physical identification.
-/

namespace OPH.EinsteinBranchInhabitant

open OPH.EinsteinBranch OPH.Thermodynamics
open Filter Asymptotics
open scoped Topology

noncomputable section

/-! ## The torus chart -/

/-- The chart point type: the sixteen points of the four-direction binary
torus. -/
abbrev TorusPoint : Type := Fin 4 → ZMod 2

/-- Coordinate step in direction `i`: add the unit vector. -/
def torusStep (i : Fin 4) (p : TorusPoint) : TorusPoint := p + Pi.single i 1

/-- The chart base point. -/
def torusBase : TorusPoint := 0

/-- Adding a single coordinate value is one symmetric step or no step. -/
theorem torus_step_reach (p : TorusPoint) (i : Fin 4) (c : ZMod 2) :
    SymmReachable torusStep p (p + Pi.single i c) := by
  have hc : c = 0 ∨ c = 1 := by revert c; decide
  rcases hc with rfl | rfl
  · rw [Pi.single_zero, add_zero]
    exact Relation.ReflTransGen.refl
  · exact Relation.ReflTransGen.single (Or.inl ⟨i, rfl⟩)

/-- Every torus point is reachable from the base in the symmetric step
closure: the chart is connected. -/
theorem torus_connected (q : TorusPoint) : SymmReachable torusStep torusBase q := by
  have h1 := torus_step_reach torusBase (0 : Fin 4) (q 0)
  have h2 := torus_step_reach (torusBase + Pi.single (0 : Fin 4) (q 0)) (1 : Fin 4) (q 1)
  have h3 := torus_step_reach
    (torusBase + Pi.single (0 : Fin 4) (q 0) + Pi.single (1 : Fin 4) (q 1)) (2 : Fin 4) (q 2)
  have h4 := torus_step_reach
    (torusBase + Pi.single (0 : Fin 4) (q 0) + Pi.single (1 : Fin 4) (q 1)
      + Pi.single (2 : Fin 4) (q 2)) (3 : Fin 4) (q 3)
  have hq : torusBase + Pi.single (0 : Fin 4) (q 0) + Pi.single (1 : Fin 4) (q 1)
      + Pi.single (2 : Fin 4) (q 2) + Pi.single (3 : Fin 4) (q 3) = q := by
    have hsum := Finset.univ_sum_single q
    rw [Fin.sum_univ_four] at hsum
    rw [torusBase, zero_add]
    exact hsum
  rw [hq] at h4
  exact ((h1.trans h2).trans h3).trans h4

/-- The diagonal stress field: entry `(i, i)` reads coordinate `i + 1`. -/
def torusStress (p : TorusPoint) : Mat 3 :=
  fun i j => if i = j then ((p (i + 1)).val : ℝ) else 0

theorem torusStress_symm (p : TorusPoint) (i j : Fin 4) :
    torusStress p i j = torusStress p j i := by
  unfold torusStress
  by_cases h : i = j
  · subst h
    rfl
  · rw [if_neg h, if_neg (Ne.symm h)]

theorem fin4_succ_ne_self (i : Fin 4) : i + 1 ≠ i := by
  revert i
  decide

/-- A step in direction `i` does not move the `(i, j)` stress entries, since
they read coordinate `i + 1`. -/
theorem torusStress_step (p : TorusPoint) (i j : Fin 4) :
    torusStress (torusStep i p) i j = torusStress p i j := by
  unfold torusStress torusStep
  rw [Pi.add_apply, Pi.single_eq_of_ne (fin4_succ_ne_self i), add_zero]

/-- The discrete Ward identity of the stress field. -/
theorem torusStress_ward (p : TorusPoint) (j : Fin 4) :
    ddiv torusStep torusStress p j = 0 := by
  unfold ddiv
  exact Finset.sum_eq_zero fun i _ => by rw [torusStress_step, sub_self]

/-- The chart Newton constant. -/
def torusNewton : ℝ := 1

/-- The chart coupling, the scale identification `κ = 8π G`. -/
def torusCoupling : ℝ := 8 * Real.pi * torusNewton

/-- The vacuum-referenced metric residue. -/
def torusLambda : ℝ := 1

/-- The geometry field `κ · stress + Λ · η`. -/
def torusGeometry (p : TorusPoint) : Mat 3 :=
  fun i j => torusCoupling * torusStress p i j + torusLambda * eta 3 i j

theorem torusGeometry_symm (p : TorusPoint) (i j : Fin 4) :
    torusGeometry p i j = torusGeometry p j i := by
  unfold torusGeometry
  rw [torusStress_symm p i j, eta_symm (n := 3) i j]

/-- The discrete contracted Bianchi identity of the geometry field. -/
theorem torusGeometry_bianchi (p : TorusPoint) (j : Fin 4) :
    ddiv torusStep torusGeometry p j = 0 := by
  unfold ddiv
  exact Finset.sum_eq_zero fun i _ => by
    unfold torusGeometry
    rw [torusStress_step]
    ring

theorem quadOf_torusGeometry (p : TorusPoint) (k : OPH.EinsteinBranch.V 3) :
    quadOf (torusGeometry p) k =
      torusCoupling * quadOf (torusStress p) k + torusLambda * quadOf (eta 3) k := by
  unfold quadOf bilinOf torusGeometry
  rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun j _ => ?_
  ring

/-- Null-null matching of the geometry and stress charges. -/
theorem torus_null_tomography (p : TorusPoint) (k : OPH.EinsteinBranch.V 3)
    (hk : quadOf (eta 3) k = 0) :
    quadOf (torusGeometry p) k = torusCoupling * quadOf (torusStress p) k := by
  rw [quadOf_torusGeometry, hk, mul_zero, add_zero]

/-- The continuum chart premise bundle of `Composition.lean`, inhabited. -/
def torusChartPremises : ContinuumEinsteinPremises TorusPoint where
  step := torusStep
  base := torusBase
  geometry := torusGeometry
  stress := torusStress
  entropyStress := torusStress
  coupling := torusCoupling
  newton := torusNewton
  referenceLambda := torusLambda
  geometry_symmetric := torusGeometry_symm
  stress_symmetric := torusStress_symm
  nullBalance := torus_null_tomography
  ward := torusStress_ward
  bianchi := torusGeometry_bianchi
  connected := torus_connected
  universalCoupling := rfl
  vacuumReference := rfl
  physicalScale := rfl

/-! ## Readouts over the demo tower -/

/-- Six typed readouts over `demoTower`.  The event and geometry channels are
the committed E7 fragment readouts (protected bit, numeric weight of the
bit); the modular, entropy, and scale channels copy the weight; the stress
channel copies the bit; the arrows are identities and one decision map; the
coarse maps are identities. -/
def demoReadouts : CommonReadoutTower demoTower where
  Geometry := fun _ => ℕ
  Modular := fun _ => ℕ
  Event := fun _ => Bool
  Stress := fun _ => Bool
  Entropy := fun _ => ℕ
  Scale := fun _ => ℕ
  geometryRead := fun _ (q : Bool × Bool) => if q.1 then 1 else 0
  modularRead := fun _ (q : Bool × Bool) => if q.1 then 1 else 0
  eventRead := fun _ (q : Bool × Bool) => q.1
  stressRead := fun _ (q : Bool × Bool) => q.1
  entropyRead := fun _ (q : Bool × Bool) => if q.1 then 1 else 0
  scaleRead := fun _ (q : Bool × Bool) => if q.1 then 1 else 0
  geometryToModular := fun _ => id
  modularToEvent := fun _ n => decide (n = 1)
  eventToStress := fun _ => id
  modularToEntropy := fun _ => id
  entropyToGeometry := fun _ => id
  geometryToScale := fun _ => id
  geometryCoarse := fun _ => id
  modularCoarse := fun _ => id
  eventCoarse := fun _ => id
  stressCoarse := fun _ => id
  entropyCoarse := fun _ => id
  scaleCoarse := fun _ => id

theorem demoReadouts_arrows : TypedArrowCommutation demoReadouts := by
  refine ⟨fun _ _ => rfl, ?_, fun _ _ => rfl, fun _ _ => rfl, fun _ _ => rfl,
    fun _ _ => rfl⟩
  rintro _ ⟨b, _⟩
  cases b <;> rfl

theorem demoReadouts_natural : ReadoutNaturality demoReadouts :=
  ⟨fun _ _ => rfl, fun _ _ => rfl, fun _ _ => rfl, fun _ _ => rfl, fun _ _ => rfl,
    fun _ _ => rfl⟩

/-- The event and geometry channels agree with the committed E7 fragment. -/
theorem demoReadouts_extends_demoFragment (r : ℕ) (q : demoTower.Quot r) :
    demoReadouts.eventRead r q = OPH.Tower.demoFragment.eventRead r q ∧
      demoReadouts.geometryRead r q = OPH.Tower.demoFragment.geometryRead r q :=
  ⟨rfl, rfl⟩

/-- The finite common-domain bundle of `Consensus.lean`, inhabited. -/
def demoAdmissibleTower : EinsteinAdmissibleTower ℕ where
  consensus := demoTower
  readouts := demoReadouts
  boundaryFiber := demoTower_boundaryFiber
  typedArrows := demoReadouts_arrows
  refinementNaturality := demoReadouts_natural

/-- The Einstein readout completion row of `QFT/InhabitedTowerRow.lean`,
closed at the enriched demo tower by the readouts above with identity
encodings of the fragment packets.  This is the relabeling route that module
anticipates; it carries no physical identification. -/
theorem demoEnrichedTower_readout_completion :
    OPH.QFT.EinsteinReadoutCompletionRow OPH.QFT.demoEnrichedTower :=
  ⟨demoReadouts, demoReadouts_arrows, demoReadouts_natural, demoTower_boundaryFiber,
    ⟨fun _ => id, fun _ => Function.injective_id, fun _ _ => rfl⟩,
    ⟨fun _ => id, fun _ => Function.injective_id, fun _ _ => rfl⟩⟩

/-! ## Small-ball reals and the scaling family -/

def ballNewton : ℝ := 1
def ballRadius : ℝ := 1
def ballStress : ℝ := 1
def ballGeometry : ℝ := 8 * Real.pi * ballNewton * ballStress
def ballBulkEntropy : ℝ := (8 * Real.pi ^ 2 * ballRadius ^ 4 / 15) * ballStress
def ballArea : ℝ := -(4 * Real.pi * ballRadius ^ 4 / 15) * ballGeometry

theorem ballNewton_pos : 0 < ballNewton := one_pos
theorem ballRadius_pos : 0 < ballRadius := one_pos

theorem ball_stationarity : ballBulkEntropy + ballArea / (4 * ballNewton) = 0 := by
  unfold ballBulkEntropy ballArea ballGeometry ballNewton ballRadius ballStress
  ring

/-- The small-ball premise bundle of `SmallBall.lean`, inhabited. -/
def ballPremises : SmallBallPremises where
  G := ballNewton
  ell := ballRadius
  stressContraction := ballStress
  geometryContraction := ballGeometry
  deltaBulkEntropy := ballBulkEntropy
  deltaArea := ballArea
  G_pos := ballNewton_pos
  ell_pos := ballRadius_pos
  stationarity := ball_stationarity
  continuumKernel := rfl
  fixedVolumeArea := rfl

/-- The shrinking radii `1 / (r + 1)`. -/
def tailRadius (r : ℕ) : ℝ := 1 / ((r : ℝ) + 1)

/-- The error channels `ℓ_r ^ 5`. -/
def tailError (r : ℕ) : ℝ := tailRadius r ^ 5

theorem tailRadius_pos (r : ℕ) : 0 < tailRadius r := by
  unfold tailRadius
  positivity

theorem tailRadius_tendsto : Tendsto tailRadius atTop (𝓝 0) :=
  tendsto_one_div_add_atTop_nhds_zero_nat

theorem tailError_isLittleO : tailError =o[atTop] fun r => tailRadius r ^ 4 := by
  have h1 : tailRadius =o[atTop] (fun _ : ℕ => (1 : ℝ)) :=
    (isLittleO_one_iff ℝ).mpr tailRadius_tendsto
  have h2 : (fun r => tailRadius r ^ 4) =O[atTop] (fun r => tailRadius r ^ 4) :=
    isBigO_refl _ _
  refine (h1.mul_isBigO h2).congr (fun r => ?_) (fun r => ?_)
  · unfold tailError
    ring
  · ring

/-- One cofinal shrinking family with nonzero error channels. -/
def tailFamily : ScalingTailData where
  ell := tailRadius
  collarError := tailError
  kernelError := tailError
  chartError := tailError

theorem tailFamily_receipt : ScalingTailReceipt tailFamily :=
  ⟨tailRadius_pos, tailRadius_tendsto, tailError_isLittleO, tailError_isLittleO,
    tailError_isLittleO⟩

/-! ## The audit candidate -/

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω] {B : Type*} [DecidableEq B]

/-- The pointwise split residual of a horizon record sums to zero. -/
theorem split_residual_zero (D : RepairLawData Ω B) (H : HorizonRecordData D) :
    ∑ x, (-Real.log (D.ref x) - (2 * Real.pi * H.boostCharge x + H.centralCharge x)) = 0 :=
  Finset.sum_eq_zero fun x _ => by rw [H.split x]; ring

/-- The audit candidate, evaluated on the inhabitant where an evaluator
exists and declared otherwise (see the module header).  The mismatch, split,
Ward, Bianchi, and tail fields are evaluated on data that could fail them;
the kernel, area, and vacuum residues compare a quantity with its own
defining expression and vanish by definition. -/
def demoCandidate (D : RepairLawData Ω B) (H : HorizonRecordData D) : BranchCandidate where
  mismatchBefore := demoMismatch (true, false)
  mismatchAfter := demoMismatch (demoTower.normalForm 0 (true, false))
  boundaryFiberCardinality := 1
  arrowSourceDigest := 0
  arrowTargetDigest := 0
  targetLeakCount := 0
  eventDimension := Fintype.card (Fin 4)
  tomographyRank := 9
  entropyDefect :=
    ∑ x, (-Real.log (D.ref x) - (2 * Real.pi * H.boostCharge x + H.centralCharge x))
  envelopeDefect := passingMaxEntTangent.entropyRate - passingMaxEntTangent.multiplier
  kernelObserved := ballBulkEntropy
  kernelExpected := (8 * Real.pi ^ 2 * ballRadius ^ 4 / 15) * ballStress
  areaObserved := ballArea
  areaExpected := -(4 * Real.pi * ballRadius ^ 4 / 15) * ballGeometry
  ell := tailRadius
  totalError := tailError
  uncoveredTimelikeDirections := 0
  wardDefect := fun j => ddiv torusStep torusStress torusBase j
  bianchiDefect := fun j => ddiv torusStep torusGeometry torusBase j
  metricCompatibilityDefect := 0
  componentCount := 1
  stressSource := 0
  geometrySource := 0
  vacuumResidue :=
    torusGeometry torusBase 0 0 -
      (torusCoupling * torusStress torusBase 0 0 + torusLambda * eta 3 0 0)
  scaleStabilizer := {1}

theorem demoCandidate_strict (D : RepairLawData Ω B) (H : HorizonRecordData D) :
    StrictManifest (demoCandidate D H) Finset.univ := by
  refine ⟨rfl, fun d _ => ?_⟩
  cases d with
  | finiteRepair =>
      show demoMismatch (demoTower.normalForm 0 (true, false)) < demoMismatch (true, false)
      decide
  | boundaryFiber => rfl
  | typedCommonDomain => rfl
  | sourceFactorization => rfl
  | lorentzGeometry => exact Fintype.card_fin 4
  | nineDirectionTomography => rfl
  | entropyFirstLaw => exact split_residual_zero D H
  | maxEntEnvelope =>
      exact sub_eq_zero.mpr (maxEnt_envelope_identity _ passingMaxEntTangent_passes)
  | diamondKernel => rfl
  | fixedVolumeArea => rfl
  | uniformTail => exact ⟨tailRadius_tendsto, tailError_isLittleO⟩
  | timelikeCoverage => rfl
  | wardIdentity => exact fun j => torusStress_ward torusBase j
  | bianchiIdentity => exact fun j => torusGeometry_bianchi torusBase j
  | metricCompatibility => rfl
  | connectedComponent => rfl
  | universalCoupling => rfl
  | vacuumReference => exact sub_self _
  | physicalScale => rfl

/-! ## The register inhabitant over a horizon record -/

/-- **The Einstein-branch register surface, inhabited over any horizon
record.**  The typed modular flow supply is the thermodynamic discharge at
the declared reference of `D` with the boost and central charges of `H`; the
remaining supplies are the synthetic witnesses above. -/
def registerOfHorizonRecord (D : RepairLawData Ω B) (H : HorizonRecordData D) :
    EinsteinRegisterData ℕ (MassZero Ω) TorusPoint where
  tower := demoTower
  readouts := demoReadouts
  protectedBoundaryFiber := demoTower_boundaryFiber
  typedCommonDomain := demoReadouts_arrows
  refinementNaturality := demoReadouts_natural
  auditCandidate := demoCandidate D H
  auditManifest := Finset.univ
  strictAudit := demoCandidate_strict D H
  modularFlowData := thermoFirstLawData Ω D.ref H.boostCharge H.centralCharge
  entropyDifferentialIsModular :=
    (thermoFirstLawData_passes D.ref H.boostCharge H.centralCharge H.split).1
  modularGeneratorSplit :=
    (thermoFirstLawData_passes D.ref H.boostCharge H.centralCharge H.split).2.1
  centralChargeIsEdge :=
    (thermoFirstLawData_passes D.ref H.boostCharge H.centralCharge H.split).2.2
  maxEntTangent := passingMaxEntTangent
  maxEntNormalization := passingMaxEntTangent_passes.1
  maxEntConstraintSpeed := passingMaxEntTangent_passes.2.1
  maxEntGibbsDifferential := passingMaxEntTangent_passes.2.2
  smallBallNewton := ballNewton
  smallBallRadius := ballRadius
  restFrameStressContraction := ballStress
  restFrameGeometryContraction := ballGeometry
  bulkEntropyVariation := ballBulkEntropy
  areaVariation := ballArea
  smallBallNewton_pos := ballNewton_pos
  smallBallRadius_pos := ballRadius_pos
  generalizedEntropyStationarity := ball_stationarity
  diamondKernelIdentification := rfl
  fixedVolumeAreaExpansion := rfl
  scalingFamily := tailFamily
  commonScalingTail := tailFamily_receipt
  chartStep := torusStep
  chartBase := torusBase
  geometryReadout := torusGeometry
  stressReadout := torusStress
  entropyStressReadout := torusStress
  couplingConstant := torusCoupling
  newtonConstant := torusNewton
  vacuumLambda := torusLambda
  geometrySymmetric := torusGeometry_symm
  stressSymmetric := torusStress_symm
  nullStressTomography := torus_null_tomography
  nullTranslationWard := torusStress_ward
  bianchiConservation := torusGeometry_bianchi
  chartConnected := torus_connected
  universalSource := rfl
  vacuumReference := rfl
  scaleIdentification := rfl

/-- The typed join, at the data level: the modular flow field of the register
is the thermodynamic first-law datum on the simplex tangent space. -/
theorem registerOfHorizonRecord_modularFlowData (D : RepairLawData Ω B)
    (H : HorizonRecordData D) :
    (registerOfHorizonRecord D H).modularFlowData =
      thermoFirstLawData Ω D.ref H.boostCharge H.centralCharge :=
  rfl

/-- The composed Einstein branch has an instance: the full conclusion bundle
at the inhabitant, for every mass-zero variation of the thermodynamic
model.  The common-readout clause is instantiated at the consistent state
`(true, true)` paired with itself; the boundary fibre of `demoTower` admits
no pairing of two distinct consistent states. -/
theorem register_composition_on_horizon_variation (D : RepairLawData Ω B)
    (H : HorizonRecordData D) (δρ : MassZero Ω) :
    EinsteinCompositionConclusion (registerOfHorizonRecord D H).readouts 0
      (true, true) (true, true) (registerOfHorizonRecord D H).modularFlowData δρ
      (registerOfHorizonRecord D H).maxEntTangent
      (registerOfHorizonRecord D H).smallBallPremises
      (registerOfHorizonRecord D H).scalingFamily
      (registerOfHorizonRecord D H).chartPremises
      (registerOfHorizonRecord D H).auditCandidate :=
  einstein_register_composition (registerOfHorizonRecord D H) 0 (true, true) (true, true)
    rfl rfl rfl δρ

/-- The `finiteFirstLaw` clause of the composed branch at the inhabitant is
the OL-B3 statement `HorizonSurface.entanglement_first_law_first_order D H δρ`:
one proposition, on one variation space, reached through the composed
Einstein branch. -/
theorem composed_first_law_is_horizon_first_law (D : RepairLawData Ω B)
    (H : HorizonRecordData D) (δρ : MassZero Ω) :
    (thermoFirstLawData Ω D.ref H.boostCharge H.centralCharge).entropyDifferential δρ =
      2 * Real.pi *
          (thermoFirstLawData Ω D.ref H.boostCharge H.centralCharge).bulkPairing δρ +
        (thermoFirstLawData Ω D.ref H.boostCharge H.centralCharge).edgeDifferential δρ :=
  (register_composition_on_horizon_variation D H δρ).finiteFirstLaw

/-- The Einstein clause at the inhabitant, with the constant equal to the
declared vacuum residue, reached through the composed conclusion.  The
statement does not mention `D` or `H`, and it holds by `rfl`
(`einstein_clause_definitional`); this theorem records only that the
composed branch returns it. -/
theorem register_einstein_equation (D : RepairLawData Ω B) (H : HorizonRecordData D)
    (p : TorusPoint) (i j : Fin 4) :
    torusGeometry p i j =
      8 * Real.pi * torusNewton * torusStress p i j + torusLambda * eta 3 i j :=
  (register_composition_on_horizon_variation D H 0).einstein p i j

/-- Degeneration receipt for the Einstein clause.  At this inhabitant the
clause is the definition of `torusGeometry` with `torusCoupling` unfolded, so
it holds by `rfl` and consumes no premise of the register.  The inhabitant
certifies that the premises are jointly satisfiable; it does not exhibit the
Einstein clause as a constraint that the chart data could have violated.  A
chart violating the clause is rejected by the null-tomography,
vacuum-reference, Ward, or Bianchi fields of the register, which is the
content of `continuumEinstein_from_explicit_premises`. -/
theorem einstein_clause_definitional (p : TorusPoint) (i j : Fin 4) :
    torusGeometry p i j =
      8 * Real.pi * torusNewton * torusStress p i j + torusLambda * eta 3 i j :=
  rfl

/-- Degeneration receipt for the small-ball clause.  `ballGeometry` is defined
as `8π · ballNewton · ballStress`, so the `smallBallArithmetic` conclusion at
this inhabitant holds by `rfl`; the three small-ball identities were solved
for the ball reals rather than tested on them. -/
theorem small_ball_clause_definitional :
    ballGeometry = 8 * Real.pi * ballNewton * ballStress :=
  rfl

/-! ## The demo repair law and horizon record -/

/-- A two-state faithful reference law `(1/3, 2/3)`. -/
def demoRef (x : Fin 2) : ℝ := if x = 0 then 1 / 3 else 2 / 3

theorem demoRef_zero : demoRef 0 = 1 / 3 := by simp [demoRef]
theorem demoRef_one : demoRef 1 = 2 / 3 := by simp [demoRef]

/-- The demo repair law (PR-07 shape): the reference above with the identity
as the repaired visible datum. -/
def demoRepairLaw : RepairLawData (Fin 2) (Fin 2) where
  ref := demoRef
  ref_pos := fun x => by
    unfold demoRef
    split_ifs <;> norm_num
  ref_law := by
    rw [Fin.sum_univ_two, demoRef_zero, demoRef_one]
    norm_num
  visible := id

/-- A nonconstant boost charge. -/
def demoBoost (x : Fin 2) : ℝ := if x = 0 then 1 else 0

/-- The central charge forced by the split. -/
def demoCentral (x : Fin 2) : ℝ := -Real.log (demoRef x) - 2 * Real.pi * demoBoost x

/-- The demo horizon record (PR-42 shape) over the demo repair law. -/
def demoHorizonRecord : HorizonRecordData demoRepairLaw where
  boostCharge := demoBoost
  centralCharge := demoCentral
  split := fun x => by
    show -Real.log (demoRef x) = 2 * Real.pi * demoBoost x + demoCentral x
    unfold demoCentral
    ring
  central_measurable := fun x y h => by rw [show x = y from h]

/-- **The explicit inhabitant.** -/
def demoEinsteinRegister : EinsteinRegisterData ℕ (MassZero (Fin 2)) TorusPoint :=
  registerOfHorizonRecord demoRepairLaw demoHorizonRecord

/-! ## Nonvacuity receipts -/

/-- The mass-zero variation `(1, -1)`. -/
def demoVariation : MassZero (Fin 2) :=
  ⟨![1, -1], by
    rw [LinearMap.mem_ker, massFunctional_apply, Fin.sum_univ_two]
    simp⟩

/-- The entropy differential of the inhabitant is nonzero on the demo
variation: the first-law clause is not the identity `0 = 0`. -/
theorem demo_entropyDifferential_ne_zero :
    demoEinsteinRegister.modularFlowData.entropyDifferential demoVariation ≠ 0 := by
  have hval : demoEinsteinRegister.modularFlowData.entropyDifferential demoVariation =
      Real.log (2 / 3) - Real.log (1 / 3) := by
    show pairing (Fin 2) (fun x => -(1 + Real.log (demoRef x))) ![1, -1] = _
    rw [pairing_apply, Fin.sum_univ_two]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, demoRef_zero, demoRef_one]
    ring
  rw [hval]
  have hlt := Real.log_lt_log (by norm_num : (0 : ℝ) < 1 / 3)
    (by norm_num : (1 / 3 : ℝ) < 2 / 3)
  exact ne_of_gt (sub_pos.mpr hlt)

/-- The bulk pairing of the inhabitant is nonzero on the demo variation. -/
theorem demo_bulkPairing_ne_zero :
    demoEinsteinRegister.modularFlowData.bulkPairing demoVariation ≠ 0 := by
  have hval : demoEinsteinRegister.modularFlowData.bulkPairing demoVariation = 1 := by
    show pairing (Fin 2) demoBoost ![1, -1] = 1
    rw [pairing_apply, Fin.sum_univ_two]
    simp [demoBoost]
  rw [hval]
  exact one_ne_zero

/-- The stress field is nonconstant across the chart. -/
theorem torusStress_nonconstant :
    torusStress (Pi.single 1 1) 0 0 ≠ torusStress torusBase 0 0 := by
  have h1 : torusStress (Pi.single 1 1) 0 0 = 1 := by
    unfold torusStress
    rw [if_pos rfl]
    have h : ((Pi.single (1 : Fin 4) (1 : ZMod 2) : TorusPoint) (0 + 1)).val = 1 := by
      decide
    rw [h, Nat.cast_one]
  have h0 : torusStress torusBase 0 0 = 0 := by
    unfold torusStress torusBase
    rw [if_pos rfl]
    have h : ((0 : TorusPoint) (0 + 1)).val = 0 := by decide
    rw [h, Nat.cast_zero]
  rw [h1, h0]
  exact one_ne_zero

/-- The readouts separate two consistent quotient states. -/
theorem demoReadouts_separate :
    demoTower.consistent 0 (true, true) ∧ demoTower.consistent 0 (false, false) ∧
      demoReadouts.geometryRead 0 (true, true) ≠ demoReadouts.geometryRead 0 (false, false) :=
  ⟨rfl, rfl, fun h => Nat.one_ne_zero h⟩

/-- The nonvacuity bundle of the explicit inhabitant: the clauses of the
register that this data can fail are exercised on nonzero data.  The Einstein,
small-ball, universal-source, and vacuum-reference conclusions hold by
definition at this value and are recorded separately by the degeneration
receipts; the stress nonconstancy below certifies the Ward clause on
nonconstant data and says nothing about the Einstein clause, which is built
from the stress by definition. -/
theorem demoEinsteinRegister_nonvacuous :
    (∃ r, ∃ x y : demoEinsteinRegister.tower.Quot r,
        demoEinsteinRegister.tower.consistent r x ∧
          ¬ demoEinsteinRegister.tower.consistent r y) ∧
    demoEinsteinRegister.modularFlowData.entropyDifferential demoVariation ≠ 0 ∧
    demoEinsteinRegister.modularFlowData.bulkPairing demoVariation ≠ 0 ∧
    demoEinsteinRegister.stressReadout (Pi.single 1 1) 0 0 ≠
      demoEinsteinRegister.stressReadout demoEinsteinRegister.chartBase 0 0 ∧
    demoEinsteinRegister.vacuumLambda ≠ 0 ∧
    demoEinsteinRegister.restFrameStressContraction ≠ 0 ∧
    demoEinsteinRegister.scalingFamily.collarError 0 ≠ 0 ∧
    demoEinsteinRegister.readouts.geometryRead 0 (true, true) ≠
      demoEinsteinRegister.readouts.geometryRead 0 (false, false) :=
  ⟨demoTower_nondegenerate, demo_entropyDifferential_ne_zero, demo_bulkPairing_ne_zero,
    torusStress_nonconstant, one_ne_zero, one_ne_zero,
    by
      show tailError 0 ≠ 0
      unfold tailError tailRadius
      norm_num,
    demoReadouts_separate.2.2⟩

/-! ## Load-bearing clause -/

/-- A uniform reference. -/
def badTau (_ : Fin 2) : ℝ := 1 / 2

/-- A zero central charge. -/
def badCentral (_ : Fin 2) : ℝ := 0

/-- With a uniform reference and zero central charge in place of the split,
the same boost charge fails the modular-generator clause on the demo
variation: the uniform reference has zero modular pairing on `(1, -1)` while
`2π` times the bulk pairing is `2π`.  No horizon record exists for this
triple, since its split field would be false; the split of the horizon
record is the load-bearing input of the join. -/
theorem split_clause_load_bearing :
    ¬ FiniteFirstLawPremises (thermoFirstLawData (Fin 2) badTau demoBoost badCentral) := by
  rintro ⟨-, hsplit, -⟩
  have h := LinearMap.congr_fun hsplit demoVariation
  have hL : (thermoFirstLawData (Fin 2) badTau demoBoost badCentral).modularPairing
      demoVariation = 0 := by
    show pairing (Fin 2) (fun x => -Real.log (badTau x)) ![1, -1] = 0
    rw [pairing_apply, Fin.sum_univ_two]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, badTau]
    ring
  have hR : ((2 * Real.pi) •
        (thermoFirstLawData (Fin 2) badTau demoBoost badCentral).bulkPairing +
      (thermoFirstLawData (Fin 2) badTau demoBoost badCentral).centralPairing)
        demoVariation = 2 * Real.pi := by
    rw [LinearMap.add_apply, LinearMap.smul_apply]
    show (2 * Real.pi) • pairing (Fin 2) demoBoost ![1, -1] +
      pairing (Fin 2) badCentral ![1, -1] = 2 * Real.pi
    rw [pairing_apply, pairing_apply, Fin.sum_univ_two, Fin.sum_univ_two]
    simp [demoBoost, badCentral]
  rw [hL, hR] at h
  exact Real.pi_ne_zero (by linarith)

end

end OPH.EinsteinBranchInhabitant

/-! ## Per-theorem axiom audit -/

#print axioms OPH.EinsteinBranchInhabitant.torus_connected
#print axioms OPH.EinsteinBranchInhabitant.torusStress_ward
#print axioms OPH.EinsteinBranchInhabitant.torusGeometry_bianchi
#print axioms OPH.EinsteinBranchInhabitant.torus_null_tomography
#print axioms OPH.EinsteinBranchInhabitant.torusChartPremises
#print axioms OPH.EinsteinBranchInhabitant.demoAdmissibleTower
#print axioms OPH.EinsteinBranchInhabitant.demoEnrichedTower_readout_completion
#print axioms OPH.EinsteinBranchInhabitant.ballPremises
#print axioms OPH.EinsteinBranchInhabitant.tailFamily_receipt
#print axioms OPH.EinsteinBranchInhabitant.demoCandidate_strict
#print axioms OPH.EinsteinBranchInhabitant.registerOfHorizonRecord
#print axioms OPH.EinsteinBranchInhabitant.registerOfHorizonRecord_modularFlowData
#print axioms OPH.EinsteinBranchInhabitant.register_composition_on_horizon_variation
#print axioms OPH.EinsteinBranchInhabitant.composed_first_law_is_horizon_first_law
#print axioms OPH.EinsteinBranchInhabitant.register_einstein_equation
#print axioms OPH.EinsteinBranchInhabitant.einstein_clause_definitional
#print axioms OPH.EinsteinBranchInhabitant.small_ball_clause_definitional
#print axioms OPH.EinsteinBranchInhabitant.demoEinsteinRegister
#print axioms OPH.EinsteinBranchInhabitant.demoEinsteinRegister_nonvacuous
#print axioms OPH.EinsteinBranchInhabitant.split_clause_load_bearing
