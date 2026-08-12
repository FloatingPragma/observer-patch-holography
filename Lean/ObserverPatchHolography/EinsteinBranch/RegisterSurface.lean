import ObserverPatchHolography.EinsteinBranch.Composition
import EinsteinPremiseLink

/-!
# Register surface for the composed Einstein branch

One citable surface where every premise of the composed Einstein branch is
named in register vocabulary.  `EinsteinRegisterData` names each physical
premise of `composedEinsteinBranch` as an explicit field, and each field
docstring ties that field to a descriptively named declared premise of the
V3 program: the typed modular flow supply, the half-sided inclusion/null
translation supply, the null-stress tomography supply, the
generalized-entropy stationarity supply, the small-ball geometry supply,
and the physical tower/calibration supply.

Every theorem in this module is a composition or re-export of existing
theorems from `Composition.lean`, `Tensor.lean`, and
`Thermodynamics/EinsteinPremiseLink.lean`.  No premise is established here,
no value of `EinsteinRegisterData` is constructed, and every conclusion is
conditional on the named supplies.  The finite first-law discharge from the
thermodynamics lane acts on the thermodynamic model's own variation space;
its hypotheses do not literally instantiate `EinsteinRegisterData`, so its
two re-exports stand side by side with the register structure inside the
same namespace, and their proofs live in
`Thermodynamics/EinsteinPremiseLink.lean`.
-/

namespace OPH.EinsteinBranch

universe u v w

noncomputable section

/-- Register-vocabulary premise surface for the composed Einstein branch.
Each field names one physical premise consumed by `composedEinsteinBranch`,
and each field docstring ties the field to a descriptively named declared
premise of the V3 program.  Together with `einstein_register_composition`
this structure is one citable surface where every premise is named in
register vocabulary.  The structure asserts nothing on its own: constructing
a value requires supplying every named premise, and this module constructs
none. -/
structure EinsteinRegisterData (ι : Type u) [Preorder ι]
    (V : Type v) [AddCommGroup V] [Module ℝ V] (P : Type u) where
  /-- Physical tower/calibration supply: the bare finite consensus tower
  carrying presentations, physical quotients, mismatch and repair data, and
  protected boundary values. -/
  tower : BareConsensusTower ι
  /-- Physical tower/calibration supply: the six typed readouts emitted from
  the repaired quotient domain of the tower. -/
  readouts : CommonReadoutTower tower
  /-- Physical tower/calibration supply: protected boundary data identify
  one consistent quotient class per regulator. -/
  protectedBoundaryFiber : BoundaryFiber tower
  /-- Physical tower/calibration supply: the typed readout arrows commute on
  the common quotient domain. -/
  typedCommonDomain : TypedArrowCommutation readouts
  /-- Physical tower/calibration supply: every typed readout commutes with
  physical refinement. -/
  refinementNaturality : ReadoutNaturality readouts
  /-- Physical tower/calibration supply: raw evaluator data of the strict
  dependency audit. -/
  auditCandidate : BranchCandidate
  /-- Physical tower/calibration supply: the declared dependency manifest of
  the audit. -/
  auditManifest : Finset Dependency
  /-- Physical tower/calibration supply: strict, deletion-sensitive
  validation of every declared dependency. -/
  strictAudit : StrictManifest auditCandidate auditManifest
  /-- Typed modular flow supply: the finite first-law pairings (entropy
  differential, modular, bulk, central, edge) on the variation space. -/
  modularFlowData : FiniteFirstLawData V
  /-- Typed modular flow supply: the entropy differential equals the modular
  pairing. -/
  entropyDifferentialIsModular :
    modularFlowData.entropyDifferential = modularFlowData.modularPairing
  /-- Typed modular flow supply: the modular generator splits as `2 * pi`
  times the bulk pairing plus the central pairing. -/
  modularGeneratorSplit :
    modularFlowData.modularPairing =
      (2 * Real.pi) • modularFlowData.bulkPairing +
        modularFlowData.centralPairing
  /-- Typed modular flow supply: the central pairing equals the edge
  differential. -/
  centralChargeIsEdge :
    modularFlowData.centralPairing = modularFlowData.edgeDifferential
  /-- Generalized-entropy stationarity supply: first-order scalar data along
  the declared MaxEnt family. -/
  maxEntTangent : MaxEntTangent
  /-- Generalized-entropy stationarity supply: the MaxEnt family preserves
  normalisation at first order. -/
  maxEntNormalization : maxEntTangent.normalizationRate = 0
  /-- Generalized-entropy stationarity supply: unit constraint speed along
  the MaxEnt family. -/
  maxEntConstraintSpeed : maxEntTangent.constraintRate = 1
  /-- Generalized-entropy stationarity supply: the Gibbs differential
  equation along the MaxEnt family. -/
  maxEntGibbsDifferential :
    maxEntTangent.entropyRate =
      maxEntTangent.multiplier * maxEntTangent.constraintRate +
        maxEntTangent.normalizationMultiplier *
          maxEntTangent.normalizationRate
  /-- Small-ball geometry supply: the Newton constant of the small-ball
  chart. -/
  smallBallNewton : ℝ
  /-- Small-ball geometry supply: the geodesic-ball radius. -/
  smallBallRadius : ℝ
  /-- Small-ball geometry supply: the rest-frame stress contraction read at
  the ball centre. -/
  restFrameStressContraction : ℝ
  /-- Small-ball geometry supply: the rest-frame geometry contraction read
  at the ball centre. -/
  restFrameGeometryContraction : ℝ
  /-- Generalized-entropy stationarity supply: the first-order bulk entropy
  variation of the ball. -/
  bulkEntropyVariation : ℝ
  /-- Generalized-entropy stationarity supply: the first-order area
  variation of the ball at fixed volume. -/
  areaVariation : ℝ
  /-- Small-ball geometry supply: positivity of the Newton constant. -/
  smallBallNewton_pos : 0 < smallBallNewton
  /-- Small-ball geometry supply: positivity of the ball radius. -/
  smallBallRadius_pos : 0 < smallBallRadius
  /-- Generalized-entropy stationarity supply: the bulk entropy variation
  balances the area variation over `4 * G` at fixed volume. -/
  generalizedEntropyStationarity :
    bulkEntropyVariation + areaVariation / (4 * smallBallNewton) = 0
  /-- Typed modular flow supply: the continuum diamond-kernel identification
  of the bulk entropy variation with the modular charge of the stress. -/
  diamondKernelIdentification :
    bulkEntropyVariation =
      (8 * Real.pi ^ 2 * smallBallRadius ^ 4 / 15) *
        restFrameStressContraction
  /-- Small-ball geometry supply: the fixed-volume small-geodesic-ball area
  identity. -/
  fixedVolumeAreaExpansion :
    areaVariation =
      -(4 * Real.pi * smallBallRadius ^ 4 / 15) *
        restFrameGeometryContraction
  /-- Small-ball geometry supply: one cofinal shrinking family with its
  named remainder channels. -/
  scalingFamily : ScalingTailData
  /-- Small-ball geometry supply: the common family carries every named
  remainder as `o(ell^4)`. -/
  commonScalingTail : ScalingTailReceipt scalingFamily
  /-- Physical tower/calibration supply: coordinate steps of one discrete
  chart on a connected component. -/
  chartStep : Fin 4 → P → P
  /-- Physical tower/calibration supply: the base point of the chart. -/
  chartBase : P
  /-- Physical tower/calibration supply: the chart-level geometry readout. -/
  geometryReadout : P → Mat 3
  /-- Physical tower/calibration supply: the chart-level reconstructed
  stress readout. -/
  stressReadout : P → Mat 3
  /-- Null-stress tomography supply: the stress source feeding the entropy
  channel. -/
  entropyStressReadout : P → Mat 3
  /-- Physical tower/calibration supply: the coupling entering the null
  balance. -/
  couplingConstant : ℝ
  /-- Physical tower/calibration supply: the Newton constant of the
  chart-level scale identification. -/
  newtonConstant : ℝ
  /-- Physical tower/calibration supply: the vacuum-referenced metric
  residue. -/
  vacuumLambda : ℝ
  /-- Null-stress tomography supply: the geometry charge form is
  symmetric. -/
  geometrySymmetric : ∀ p i j, geometryReadout p i j = geometryReadout p j i
  /-- Null-stress tomography supply: the stress charge form is symmetric. -/
  stressSymmetric : ∀ p i j, stressReadout p i j = stressReadout p j i
  /-- Null-stress tomography supply: null-null matching of the geometry and
  stress charges on the null cone, the scaling-limit null-horizon
  balance. -/
  nullStressTomography : ∀ p k, quadOf (eta 3) k = 0 →
    quadOf (geometryReadout p) k =
      couplingConstant * quadOf (stressReadout p) k
  /-- Half-sided inclusion/null translation supply: the weak Ward identity
  for the reconstructed stress, the conservation shadow of the null
  translations generated by half-sided inclusion. -/
  nullTranslationWard : ∀ p j, ddiv chartStep stressReadout p j = 0
  /-- Small-ball geometry supply: the contracted Bianchi identity of the
  smooth geometric side. -/
  bianchiConservation : ∀ p j, ddiv chartStep geometryReadout p j = 0
  /-- Physical tower/calibration supply: connectivity of the chart component
  in the symmetric step closure. -/
  chartConnected : ∀ q, SymmReachable chartStep chartBase q
  /-- Null-stress tomography supply: universal coupling; the entropy channel
  and the null tomography consume one stress source. -/
  universalSource : entropyStressReadout = stressReadout
  /-- Physical tower/calibration supply: the vacuum-reference condition
  fixing the pointwise metric residue at the chart base. -/
  vacuumReference :
    geometryReadout chartBase 0 0 =
      couplingConstant * stressReadout chartBase 0 0 +
        vacuumLambda * eta 3 0 0
  /-- Physical tower/calibration supply: the independent scale
  identification of the coupling. -/
  scaleIdentification : couplingConstant = 8 * Real.pi * newtonConstant

/-- The finite first-law premise package of `Entropy.lean`, assembled from
the typed modular flow supply fields.  Pure repackaging; each conjunct is a
register field. -/
theorem EinsteinRegisterData.firstLawPremises
    {ι : Type u} [Preorder ι] {V : Type v} [AddCommGroup V] [Module ℝ V]
    {P : Type u} (D : EinsteinRegisterData ι V P) :
    FiniteFirstLawPremises D.modularFlowData :=
  ⟨D.entropyDifferentialIsModular, D.modularGeneratorSplit,
    D.centralChargeIsEdge⟩

/-- The MaxEnt envelope premise package of `Entropy.lean`, assembled from
the generalized-entropy stationarity supply fields.  Pure repackaging; each
conjunct is a register field. -/
theorem EinsteinRegisterData.maxEntPremises
    {ι : Type u} [Preorder ι] {V : Type v} [AddCommGroup V] [Module ℝ V]
    {P : Type u} (D : EinsteinRegisterData ι V P) :
    MaxEntEnvelopePremises D.maxEntTangent :=
  ⟨D.maxEntNormalization, D.maxEntConstraintSpeed, D.maxEntGibbsDifferential⟩

/-- The small-ball premise bundle of `SmallBall.lean`, assembled from the
generalized-entropy stationarity, typed modular flow, and small-ball
geometry supply fields.  Pure repackaging; each field of the bundle is a
register field. -/
def EinsteinRegisterData.smallBallPremises
    {ι : Type u} [Preorder ι] {V : Type v} [AddCommGroup V] [Module ℝ V]
    {P : Type u} (D : EinsteinRegisterData ι V P) : SmallBallPremises where
  G := D.smallBallNewton
  ell := D.smallBallRadius
  stressContraction := D.restFrameStressContraction
  geometryContraction := D.restFrameGeometryContraction
  deltaBulkEntropy := D.bulkEntropyVariation
  deltaArea := D.areaVariation
  G_pos := D.smallBallNewton_pos
  ell_pos := D.smallBallRadius_pos
  stationarity := D.generalizedEntropyStationarity
  continuumKernel := D.diamondKernelIdentification
  fixedVolumeArea := D.fixedVolumeAreaExpansion

/-- The continuum chart premise bundle of `Composition.lean`, assembled from
the null-stress tomography, half-sided inclusion/null translation,
small-ball geometry, and physical tower/calibration supply fields.  Pure
repackaging; each field of the bundle is a register field. -/
def EinsteinRegisterData.chartPremises
    {ι : Type u} [Preorder ι] {V : Type v} [AddCommGroup V] [Module ℝ V]
    {P : Type u} (D : EinsteinRegisterData ι V P) :
    ContinuumEinsteinPremises P where
  step := D.chartStep
  base := D.chartBase
  geometry := D.geometryReadout
  stress := D.stressReadout
  entropyStress := D.entropyStressReadout
  coupling := D.couplingConstant
  newton := D.newtonConstant
  referenceLambda := D.vacuumLambda
  geometry_symmetric := D.geometrySymmetric
  stress_symmetric := D.stressSymmetric
  nullBalance := D.nullStressTomography
  ward := D.nullTranslationWard
  bianchi := D.bianchiConservation
  connected := D.chartConnected
  universalCoupling := D.universalSource
  vacuumReference := D.vacuumReference
  physicalScale := D.scaleIdentification

/-- Register restatement of the composed Einstein branch.  Consuming the
register surface re-exports the full conclusion bundle
`EinsteinCompositionConclusion`: the strict dependency audit, the
common-readout composition, the typed arrows and refinement naturality, the
finite first law, the MaxEnt envelope, the exact small-ball arithmetic, the
scaling-tail receipt, the universal source, and the componentwise Einstein
equation with the vacuum-referenced constant.  This theorem is a composition
of `composedEinsteinBranch` with the field-to-premise assembly above; it
proves the implication on the admissible branch and does not prove that a
tower satisfying the premises exists. -/
theorem einstein_register_composition
    {ι : Type u} [Preorder ι] {V : Type v} [AddCommGroup V] [Module ℝ V]
    {P : Type u} (D : EinsteinRegisterData ι V P)
    (r : ι) (q q' : D.tower.Quot r)
    (hq : D.tower.consistent r q) (hq' : D.tower.consistent r q')
    (hBoundary : D.tower.boundary r q = D.tower.boundary r q')
    (δρ : V) :
    EinsteinCompositionConclusion D.readouts r q q' D.modularFlowData δρ
      D.maxEntTangent D.smallBallPremises D.scalingFamily D.chartPremises
      D.auditCandidate :=
  composedEinsteinBranch D.readouts D.auditCandidate D.auditManifest
    D.strictAudit D.protectedBoundaryFiber D.typedCommonDomain
    D.refinementNaturality r q q' hq hq' hBoundary D.modularFlowData
    D.firstLawPremises δρ D.maxEntTangent D.maxEntPremises
    D.smallBallPremises D.scalingFamily D.commonScalingTail D.chartPremises

/-- Ward/Bianchi Lambda-constancy on the register surface: the null-stress
tomography supply, the half-sided inclusion/null translation supply, the
contracted Bianchi identity, and chart connectivity force one constant
metric coefficient across the connected chart.  Re-export of
`einstein_equation_with_constant_symm` at the register fields. -/
theorem einstein_register_lambda_constancy
    {ι : Type u} [Preorder ι] {V : Type v} [AddCommGroup V] [Module ℝ V]
    {P : Type u} (D : EinsteinRegisterData ι V P) :
    ∃ Λ : ℝ, ∀ p i j,
      D.geometryReadout p i j =
        D.couplingConstant * D.stressReadout p i j + Λ * eta 3 i j :=
  einstein_equation_with_constant_symm D.chartStep D.chartBase
    D.chartConnected D.geometryReadout D.stressReadout D.couplingConstant
    D.geometrySymmetric D.stressSymmetric D.nullStressTomography
    D.bianchiConservation D.nullTranslationWard

/-- Finite first-law discharge, cited from the thermodynamics lane.  On the
finite thermodynamic model the typed modular flow supply is a theorem: the
premise package `FiniteFirstLawPremises` holds for the simplex-tangent data
under the declared split of the faithful reference into `2 * pi` times a
bulk charge plus a central charge.  This corollary re-exports
`OPH.Thermodynamics.thermoFirstLawData_passes` from
`Thermodynamics/EinsteinPremiseLink.lean`; the proof lives there.  Its
hypotheses act on the model's own variation space and do not literally
instantiate `EinsteinRegisterData`, so the re-export stands side by side
with the register structure. -/
theorem einstein_register_first_law_discharge
    {Ω : Type w} [Fintype Ω] (tau Bc Z : Ω → ℝ)
    (hsplit : ∀ x, -Real.log (tau x) = 2 * Real.pi * Bc x + Z x) :
    FiniteFirstLawPremises
      (OPH.Thermodynamics.thermoFirstLawData Ω tau Bc Z) :=
  OPH.Thermodynamics.thermoFirstLawData_passes tau Bc Z hsplit

/-- The discharged first law applied through the Einstein-branch identity:
on every mass-preserving variation of the finite thermodynamic model the
entropy differential is `2 * pi` times the bulk pairing plus the edge
differential.  Re-export of
`OPH.Thermodynamics.thermo_first_law_on_simplex_tangent`; the proof lives in
`Thermodynamics/EinsteinPremiseLink.lean`. -/
theorem einstein_register_first_law_on_model
    {Ω : Type w} [Fintype Ω] (tau Bc Z : Ω → ℝ)
    (hsplit : ∀ x, -Real.log (tau x) = 2 * Real.pi * Bc x + Z x)
    (δρ : OPH.Thermodynamics.MassZero Ω) :
    (OPH.Thermodynamics.thermoFirstLawData Ω tau Bc Z).entropyDifferential
        δρ =
      2 * Real.pi *
          (OPH.Thermodynamics.thermoFirstLawData Ω tau Bc Z).bulkPairing δρ
        + (OPH.Thermodynamics.thermoFirstLawData Ω tau Bc
            Z).edgeDifferential δρ :=
  OPH.Thermodynamics.thermo_first_law_on_simplex_tangent tau Bc Z hsplit δρ

/-! ## Per-theorem axiom audit -/

#print axioms EinsteinRegisterData.firstLawPremises
#print axioms EinsteinRegisterData.maxEntPremises
#print axioms einstein_register_composition
#print axioms einstein_register_lambda_constancy
#print axioms einstein_register_first_law_discharge
#print axioms einstein_register_first_law_on_model

end

end OPH.EinsteinBranch
