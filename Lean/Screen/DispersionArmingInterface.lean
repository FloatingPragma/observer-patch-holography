import SeamCurrentEdge30Moment
import PrimitivePortScaleBoundary

set_option autoImplicit false

namespace OPH.DispersionArming

open OPH.SeamCurrentEdge30Moment
open OPH.A5OrbitRaySeparation
open OPH.PrimitivePortScaleBoundary

/-!
# Schema-level arming checklist for the frozen dispersion row FZ-12

The frozen prediction register row FZ-12 fixes the source-seam coefficient
manifold and names the physical premises required before a comparison can be
armed.  This module mirrors those clauses as placeholder fields around the
committed coefficient and control-support objects.

## What is proved

* For every externally supplied nonzero scale, the three coefficient slots
  obey the exact committed scale-free ratios; the edge and vertex branches
  have opposite rank-six signs.
* `DispersionArmingInterface` is a schema with one placeholder per protocol
  clause.  `ArmedComparison` is definitionally just such a record paired with
  a natural-number contract-token placeholder.  Its equivalence theorem is
  packaging, not a physical eligibility theorem.
* `stipulatedInterface` constructs mock inhabitants by setting each `Prop`
  field to `True` and supplying arbitrary positive scale and lower-bound
  labels.  These inhabitants show that this permissive schema constrains
  neither label.  Propositional extensionality shows only that inhabited
  `Prop` fields carry no data separation in this representation.
* The coefficient ratios alone select no unique positive scale.

## What is not proved here

No arming occurred.  The module has stipulated full inhabitants, but no
source-derived, artifact-bound, custody-backed evidential inhabitant.  The
mock `Prop` proofs, natural-number digest, scale, and lower-bound labels do not
establish a physical action, sector attachment, gluing, frame transport,
analysis contract, source scale, source lower bound, provenance, or custody.
The schema-freedom receipts are not universal no-go theorems against a future
source derivation.  A future evidential type must carry concrete witnesses and
artifact-bound derivation and contract receipts.  The frozen FZ-12 bytes are
untouched, no comparison datum enters, and the row remains unarmed.
-/

/-! ## The edge-branch coefficient slot at a declared scale -/

/-- Edge-branch anisotropic sixth-order coefficient at a declared scale
`a`.  The committed vertex-branch counterpart is
`OPH.PrimitivePortScaleBoundary.scaleB6`; the committed dimensionless
edge value is `OPH.A5OrbitRaySeparation.B6 .edge30 = -1/12600`. -/
noncomputable def edgeScaleB6 (a : ℝ) : ℝ := -(a ^ 4) / 12600

/-- The committed seam-derived rational row has the frozen dimensionless
values. -/
theorem seam_row_values :
    seamLambdaC4 = -1 / 20 ∧ seamLambdaB0 = 1 / 840 ∧
      seamLambdaB6 = -1 / 12600 := by
  refine ⟨?_, ?_, ?_⟩ <;>
    norm_num [seamLambdaC4, seamLambdaB0, seamLambdaB6, edge30Weight]

/-- At every declared scale the three coefficient slots equal the
committed seam-derived rational row times the declared power of the
scale.  This is the committed carrier scale slot of the arming path. -/
theorem edgeScale_matches_seam_row (a : ℝ) :
    scaleC4 a = a ^ 2 * ((seamLambdaC4 : ℚ) : ℝ) ∧
      scaleB0 a = a ^ 4 * ((seamLambdaB0 : ℚ) : ℝ) ∧
      edgeScaleB6 a = a ^ 4 * ((seamLambdaB6 : ℚ) : ℝ) := by
  obtain ⟨h1, h2, h3⟩ := seam_row_values
  rw [h1, h2, h3]
  push_cast
  refine ⟨?_, ?_, ?_⟩
  · simp only [scaleC4]; ring
  · simp only [scaleB0]; ring
  · simp only [edgeScaleB6]; ring

/-- At every nonzero scale the slot ratios are exactly the frozen FZ-12
manifold ratios.  The ratios are scale free, which is the content of the
ratio-only scale-nonselection receipts below. -/
theorem edgeScale_ray (a : ℝ) (ha : a ≠ 0) :
    scaleB0 a / scaleC4 a ^ 2 = 10 / 21 ∧
      edgeScaleB6 a / scaleC4 a ^ 2 = -2 / 63 ∧
      edgeScaleB6 a / scaleB0 a = -1 / 15 := by
  refine ⟨scaleB0_over_scaleC4_sq a ha, ?_, ?_⟩
  · simp only [edgeScaleB6, scaleC4]
    field_simp
    ring
  · simp only [edgeScaleB6, scaleB0]
    field_simp
    norm_num

/-- Distinct positive scales give distinct quartic coefficients: the
declared scale is a faithful parameter of the coefficient slot. -/
theorem scale_parameter_faithful (a b : ℝ) (ha : 0 < a) (hb : 0 < b)
    (h : scaleC4 a = scaleC4 b) : a = b := by
  have h2 : a ^ 2 = b ^ 2 := by
    simp only [scaleC4] at h
    linarith
  have := abs_eq_abs.mp (by
    have := sq_abs a ▸ sq_abs b ▸ h2
    nlinarith [abs_nonneg a, abs_nonneg b, sq_abs a, sq_abs b] :
      |a| = |b|)
  rcases this with h' | h'
  · exact h'
  · nlinarith

/-- **The committed ray selects no scale.**  The frozen manifold ratios
hold at every positive scale, so no positive scale is uniquely selected by
them.  This schema supplies no scale choice; an independent future source
theorem may supply one. -/
theorem committed_ray_selects_no_scale :
    ¬ ∃! a : ℝ, 0 < a ∧
      (scaleB0 a / scaleC4 a ^ 2 = 10 / 21 ∧
        edgeScaleB6 a / scaleC4 a ^ 2 = -2 / 63 ∧
        edgeScaleB6 a / scaleB0 a = -1 / 15) := by
  rintro ⟨a, -, huniq⟩
  have h1 : (1 : ℝ) = a := huniq 1 ⟨one_pos, edgeScale_ray 1 one_ne_zero⟩
  have h2 : (2 : ℝ) = a := huniq 2 ⟨two_pos, edgeScale_ray 2 two_ne_zero⟩
  have h12 : (1 : ℝ) = 2 := h1.trans h2.symm
  norm_num at h12

/-- The edge and vertex branches keep opposite rank-six signs at every
nonzero scale, mirroring the frozen row's separation from FZ-11. -/
theorem edge_vertex_rank_six_sign_separation (a : ℝ) (ha : a ≠ 0) :
    edgeScaleB6 a < 0 ∧ 0 < scaleB6 a := by
  have h4 : 0 < a ^ 4 := by positivity
  constructor
  · simp only [edgeScaleB6]
    linarith
  · simp only [scaleB6]
    linarith

/-! ## The committed edge-orbit support certificate -/

/-- Statement of the committed edge-orbit support certificate: the
complete thirty-seam support, one common squared chord, the projective
class multiset identity with the edge orbit, and the identification of
the seam-derived coefficients with the registered edge-30 orbit row.
This is the committed counterpart of the protocol clause "the complete
edge orbit as its sole direct support". -/
abbrev EdgeOrbitSupportCertificate : Prop :=
  (Fintype.card (Fin 30) = 30) ∧
    (∀ e : Fin 30,
      OPH.PrimitivePortTranslationBridge.dot
        (carrierSeamDifference e) (carrierSeamDifference e) = 4) ∧
    ((Finset.univ.val.map seamAxisClass : Multiset (Fin 15)) =
      Finset.univ.val.map edgeAxisClass) ∧
    seamLambdaC4 = C4 ∧ seamLambdaB0 = B0 ∧ seamLambdaB6 = B6 .edge30

/-! ## The arming interface -/

/-- **Typed arming interface for row FZ-12.**  One field per clause of
the frozen comparison protocol.  The premise clauses are typed
obligations: `Prop` fields required to hold whose content the committed
corpus does not determine, and positive real data fields whose values the
committed corpus does not determine.  The determined clauses are `Prop`
fields whose content the committed corpus proves at every declared
positive scale.  An inhabitant certifies no arming. -/
structure DispersionArmingInterface where
  /-- Protocol clause "a source-derived homogeneous position action with
  the complete edge orbit as its sole direct support": the physical
  position action is homogeneous and directly supported on the complete
  edge orbit alone.  This schema does not supply or determine the
  proposition; an independent future source theorem may do so. -/
  homogeneousPositionAction : Prop
  /-- The position-action clause holds. -/
  homogeneousPositionAction_holds : homogeneousPositionAction
  /-- Protocol clause "equal source weights": the physical action weights
  the thirty directions equally.  This schema does not supply or determine
  the proposition; the committed control normalization
  `OPH.SeamCurrentEdge30Moment.edge30Weight` is equal-weight, and no
  committed theorem attaches it to a physical action. -/
  equalSourceWeights : Prop
  /-- The equal-weights clause holds. -/
  equalSourceWeights_holds : equalSourceWeights
  /-- Protocol clause "continuous field and same-operator sector
  attachment": the compared field is continuous and the comparison reads
  the same operator sector that the source theorem fixes.  This schema does
  not supply or determine the proposition. -/
  sectorAttachment : Prop
  /-- The sector-attachment clause holds. -/
  sectorAttachment_holds : sectorAttachment
  /-- Protocol clause "cofinal gluing": the finite patch embeds
  cofinally into the compared regime.  This schema does not supply or
  determine the proposition. -/
  cofinalGluing : Prop
  /-- The cofinal-gluing clause holds. -/
  cofinalGluing_holds : cofinalGluing
  /-- Protocol clause "finite scale": one declared positive carrier
  scale filling the committed coefficient slot.  The value is not supplied
  or determined by this schema (`committed_ray_selects_no_scale`,
  `interface_scale_free_parameter`). -/
  carrierScale : ℝ
  /-- The declared scale is positive. -/
  carrierScale_pos : 0 < carrierScale
  /-- Protocol clause "coherent frame and boost transport": one coherent
  frame for the comparison and a boost transport connecting it to the
  source frame.  This schema does not supply or determine the proposition. -/
  frameBoostTransport : Prop
  /-- The frame-transport clause holds. -/
  frameBoostTransport_holds : frameBoostTransport
  /-- Protocol clause "frozen nuisance and coverage rules": the readout
  rule of the comparison, with nuisance isolation and calibrated
  coverage, fixed before exposure.  This schema does not supply or determine
  the proposition. -/
  nuisanceCoverageRules : Prop
  /-- The readout-rule clause holds. -/
  nuisanceCoverageRules_holds : nuisanceCoverageRules
  /-- Separate exclusivity obligation, protocol clause "the source-derived
  positive lower bound owned by issue #664": one positive lower bound for
  the admissible scales in comparison units.  The value is not supplied or
  determined by this schema (`committed_corpus_selects_no_lower_bound`); an
  independent future source theorem may supply it. -/
  exclusivityLowerBound : ℝ
  /-- The lower bound is positive. -/
  exclusivityLowerBound_pos : 0 < exclusivityLowerBound
  /-- The declared scale lies in the admissible manifold the exclusivity
  obligation quantifies over: the bound does not exceed the declared
  scale. -/
  exclusivityLowerBound_le_scale : exclusivityLowerBound ≤ carrierScale
  /-- The clause that the lower bound is source derived, in the sense of
  the append-only clarification of the frozen row.  This schema does not
  supply or determine the proposition; a stipulated inhabitant supplies no
  source derivation, while an independent future source theorem remains
  admissible. -/
  lowerBoundSourceDerived : Prop
  /-- The source-derivation clause holds. -/
  lowerBoundSourceDerived_holds : lowerBoundSourceDerived
  /-- Determined clause: the frozen manifold ratios hold at the declared
  scale.  The committed corpus proves this field for every positive
  scale (`edgeScale_ray`). -/
  manifoldRay :
    scaleB0 carrierScale / scaleC4 carrierScale ^ 2 = 10 / 21 ∧
      edgeScaleB6 carrierScale / scaleC4 carrierScale ^ 2 = -2 / 63 ∧
      edgeScaleB6 carrierScale / scaleB0 carrierScale = -1 / 15
  /-- Determined clause: the committed edge-orbit support certificate
  (`source_seam_edge30_control_certificate`). -/
  edgeOrbitSupport : EdgeOrbitSupportCertificate
  /-- Determined clause: the coefficient slots at the declared scale equal
  the committed seam-derived rational row times the declared power of the
  scale (`edgeScale_matches_seam_row`). -/
  coefficientSlot :
    scaleC4 carrierScale = carrierScale ^ 2 * ((seamLambdaC4 : ℚ) : ℝ) ∧
      scaleB0 carrierScale = carrierScale ^ 4 * ((seamLambdaB0 : ℚ) : ℝ) ∧
      edgeScaleB6 carrierScale =
        carrierScale ^ 4 * ((seamLambdaB6 : ℚ) : ℝ)

/-- The interface inhabitant with stipulated external fields over the
committed determined content.  The premise `Prop` fields are stipulated
as `True`, and the scale and lower bound are stipulated parameters: the
definition certifies no arming, no source derivation, and no physical
scale. -/
noncomputable def stipulatedInterface (a m : ℝ) (ha : 0 < a) (hm : 0 < m)
    (hle : m ≤ a) : DispersionArmingInterface where
  homogeneousPositionAction := True
  homogeneousPositionAction_holds := trivial
  equalSourceWeights := True
  equalSourceWeights_holds := trivial
  sectorAttachment := True
  sectorAttachment_holds := trivial
  cofinalGluing := True
  cofinalGluing_holds := trivial
  carrierScale := a
  carrierScale_pos := ha
  frameBoostTransport := True
  frameBoostTransport_holds := trivial
  nuisanceCoverageRules := True
  nuisanceCoverageRules_holds := trivial
  exclusivityLowerBound := m
  exclusivityLowerBound_pos := hm
  exclusivityLowerBound_le_scale := hle
  lowerBoundSourceDerived := True
  lowerBoundSourceDerived_holds := trivial
  manifoldRay := edgeScale_ray a (ne_of_gt ha)
  edgeOrbitSupport := source_seam_edge30_control_certificate
  coefficientSlot := edgeScale_matches_seam_row a

@[simp] theorem stipulatedInterface_carrierScale (a m : ℝ) (ha : 0 < a)
    (hm : 0 < m) (hle : m ≤ a) :
    (stipulatedInterface a m ha hm hle).carrierScale = a := rfl

@[simp] theorem stipulatedInterface_exclusivityLowerBound (a m : ℝ)
    (ha : 0 < a) (hm : 0 < m) (hle : m ≤ a) :
    (stipulatedInterface a m ha hm hle).exclusivityLowerBound = m := rfl

/-- Explicit nontrivial inhabitant with both stipulated parameters equal
to the label `1`.  It certifies no arming. -/
noncomputable def unitStipulatedInterface : DispersionArmingInterface :=
  stipulatedInterface 1 1 one_pos one_pos le_rfl

theorem dispersionArmingInterface_nonempty :
    Nonempty DispersionArmingInterface := ⟨unitStipulatedInterface⟩

/-! ## Underdetermination inside the placeholder schema -/

/-- **The scale is a free parameter.**  Every positive real is realized
as the carrier scale of an interface inhabitant. -/
theorem interface_scale_free_parameter (a : ℝ) (ha : 0 < a) :
    ∃ I : DispersionArmingInterface, I.carrierScale = a :=
  ⟨stipulatedInterface a a ha ha le_rfl, rfl⟩

/-- **The exclusivity lower bound is a free parameter.**  Every positive
real is realized as the lower bound of an interface inhabitant. -/
theorem interface_lower_bound_free_parameter (m : ℝ) (hm : 0 < m) :
    ∃ I : DispersionArmingInterface,
      I.exclusivityLowerBound = m ∧ I.carrierScale = m :=
  ⟨stipulatedInterface m m hm hm le_rfl, rfl, rfl⟩

/-- **The placeholder schema does not determine the scale.**  Two mock
inhabitants agree at the exclusivity lower-bound label and differ at the
carrier-scale label.  This is not a no-go against future source derivations. -/
theorem committed_corpus_does_not_determine_scale :
    ∃ I₁ I₂ : DispersionArmingInterface,
      I₁.exclusivityLowerBound = I₂.exclusivityLowerBound ∧
        I₁.carrierScale ≠ I₂.carrierScale := by
  refine ⟨stipulatedInterface 1 1 one_pos one_pos le_rfl,
    stipulatedInterface 2 1 two_pos one_pos (by norm_num), rfl, ?_⟩
  simp only [stipulatedInterface_carrierScale]
  norm_num

/-- **The placeholder schema does not determine the lower bound.**  Two mock
inhabitants agree at the carrier-scale label and differ at the lower-bound
label.  This is not a universal source no-go. -/
theorem committed_corpus_does_not_determine_lower_bound :
    ∃ I₁ I₂ : DispersionArmingInterface,
      I₁.carrierScale = I₂.carrierScale ∧
        I₁.exclusivityLowerBound ≠ I₂.exclusivityLowerBound := by
  refine ⟨stipulatedInterface 2 1 two_pos one_pos (by norm_num),
    stipulatedInterface 2 2 two_pos two_pos le_rfl, rfl, ?_⟩
  simp only [stipulatedInterface_exclusivityLowerBound]
  norm_num

/-- **The permissive schema singles out no lower bound.**  Every positive
real occurs as the lower-bound label of a stipulated inhabitant.  This proves
only that the present record definition does not constrain the value. -/
theorem committed_corpus_selects_no_lower_bound :
    ¬ ∃! m : ℝ, 0 < m ∧
      ∃ I : DispersionArmingInterface, I.exclusivityLowerBound = m := by
  rintro ⟨m, -, huniq⟩
  have h1 : (1 : ℝ) = m :=
    huniq 1 ⟨one_pos, stipulatedInterface 1 1 one_pos one_pos le_rfl, rfl⟩
  have h2 : (2 : ℝ) = m :=
    huniq 2 ⟨two_pos, stipulatedInterface 2 2 two_pos two_pos le_rfl, rfl⟩
  have h12 : (1 : ℝ) = 2 := h1.trans h2.symm
  norm_num at h12

/-- **The required-to-hold premise propositions carry no data separation.**
Any two interface inhabitants have propositionally equal premise `Prop`
fields: each clause holds in both, so propositional extensionality
identifies them.  This reflects the choice of inhabited `Prop` placeholders;
it is not evidence that the physical clauses hold or that real custody can
only be represented by the numeric labels. -/
theorem prop_fields_carry_no_separation (I₁ I₂ : DispersionArmingInterface) :
    I₁.homogeneousPositionAction = I₂.homogeneousPositionAction ∧
      I₁.equalSourceWeights = I₂.equalSourceWeights ∧
      I₁.sectorAttachment = I₂.sectorAttachment ∧
      I₁.cofinalGluing = I₂.cofinalGluing ∧
      I₁.frameBoostTransport = I₂.frameBoostTransport ∧
      I₁.nuisanceCoverageRules = I₂.nuisanceCoverageRules ∧
      I₁.lowerBoundSourceDerived = I₂.lowerBoundSourceDerived :=
  ⟨propext ⟨fun _ => I₂.homogeneousPositionAction_holds,
      fun _ => I₁.homogeneousPositionAction_holds⟩,
    propext ⟨fun _ => I₂.equalSourceWeights_holds,
      fun _ => I₁.equalSourceWeights_holds⟩,
    propext ⟨fun _ => I₂.sectorAttachment_holds,
      fun _ => I₁.sectorAttachment_holds⟩,
    propext ⟨fun _ => I₂.cofinalGluing_holds,
      fun _ => I₁.cofinalGluing_holds⟩,
    propext ⟨fun _ => I₂.frameBoostTransport_holds,
      fun _ => I₁.frameBoostTransport_holds⟩,
    propext ⟨fun _ => I₂.nuisanceCoverageRules_holds,
      fun _ => I₁.nuisanceCoverageRules_holds⟩,
    propext ⟨fun _ => I₂.lowerBoundSourceDerived_holds,
      fun _ => I₁.lowerBoundSourceDerived_holds⟩⟩

/-! ## The post-custody contract token and the armed comparison -/

/-- Abstract token for the protocol clause "a post-custody
dataset-specific contract".  The digest stands for a contract artifact of
one eligible dataset; this schema does not bind or authenticate its value,
and an inhabitant certifies no contract. -/
structure PostCustodyDatasetContract where
  /-- Abstract digest of the dataset-specific contract.  This token does not
  determine or authenticate its value. -/
  contractDigest : ℕ

/-- Explicit nontrivial inhabitant with a stipulated nonzero digest.  It
certifies no contract. -/
def stipulatedContract : PostCustodyDatasetContract := ⟨1⟩

/-- **The contract digest is a free parameter.**  Every natural number is
realized as the digest of a token. -/
theorem contract_digest_free_parameter (d : ℕ) :
    ∃ C : PostCustodyDatasetContract, C.contractDigest = d := ⟨⟨d⟩, rfl⟩

/-- Two distinct contract stipulations, both consistent with the schema:
the token alone selects no contract artifact. -/
theorem contract_distinct_stipulations :
    ∃ C₁ C₂ : PostCustodyDatasetContract, C₁ ≠ C₂ := by
  refine ⟨⟨0⟩, ⟨1⟩, fun h => ?_⟩
  exact Nat.zero_ne_one
    (congrArg PostCustodyDatasetContract.contractDigest h)

/-- **An armed comparison, typed.**  Exactly one interface inhabitant and
one post-custody dataset contract token.  An inhabitant of this structure
built from stipulations certifies no arming; the frozen row remains
ineligible and unarmed. -/
structure ArmedComparison where
  /-- The typed arming interface. -/
  armingInterface : DispersionArmingInterface
  /-- The post-custody dataset-specific contract token. -/
  datasetContract : PostCustodyDatasetContract

/-- Explicit nontrivial inhabitant assembled from the stipulated
interface and the stipulated contract.  It certifies no arming. -/
noncomputable def stipulatedArmedComparison : ArmedComparison :=
  ⟨unitStipulatedInterface, stipulatedContract⟩

/-- **Packaging, exact form.**  This record is definitionally an interface
inhabitant paired with a contract-token placeholder; no physical eligibility
follows. -/
def armedComparisonEquiv :
    ArmedComparison ≃ DispersionArmingInterface × PostCustodyDatasetContract where
  toFun A := (A.armingInterface, A.datasetContract)
  invFun p := ⟨p.1, p.2⟩
  left_inv _ := rfl
  right_inv _ := rfl

/-- Nonemptiness of the packaging record is equivalent to nonemptiness of its
two placeholder components. -/
theorem armedComparison_constructible_iff :
    Nonempty ArmedComparison ↔
      Nonempty DispersionArmingInterface ∧
        Nonempty PostCustodyDatasetContract := by
  constructor
  · rintro ⟨A⟩
    exact ⟨⟨A.armingInterface⟩, ⟨A.datasetContract⟩⟩
  · rintro ⟨⟨I⟩, ⟨C⟩⟩
    exact ⟨⟨I, C⟩⟩

/-! ## The eligibility antecedent -/

/-- The conjunction of the protocol premise clauses of one interface
inhabitant: the eligibility antecedent of the frozen comparison protocol,
with the finite-scale and exclusivity clauses read on the data fields. -/
def EligibilityAntecedent (I : DispersionArmingInterface) : Prop :=
  I.homogeneousPositionAction ∧ I.equalSourceWeights ∧
    I.sectorAttachment ∧ I.cofinalGluing ∧ 0 < I.carrierScale ∧
    I.frameBoostTransport ∧ I.nuisanceCoverageRules ∧
    (0 < I.exclusivityLowerBound ∧
      I.exclusivityLowerBound ≤ I.carrierScale ∧
      I.lowerBoundSourceDerived)

/-- Every interface inhabitant projects proofs of its stored proposition
fields.  This tautological schema fact transfers no eligibility to a physical
comparison, because the mock inhabitants simply set those fields to `True`. -/
theorem eligibilityAntecedent_of_interface (I : DispersionArmingInterface) :
    EligibilityAntecedent I :=
  ⟨I.homogeneousPositionAction_holds, I.equalSourceWeights_holds,
    I.sectorAttachment_holds, I.cofinalGluing_holds, I.carrierScale_pos,
    I.frameBoostTransport_holds, I.nuisanceCoverageRules_holds,
    I.exclusivityLowerBound_pos, I.exclusivityLowerBound_le_scale,
    I.lowerBoundSourceDerived_holds⟩

/-- The antecedent read on the interface carried by an armed
comparison. -/
theorem armed_eligibilityAntecedent (A : ArmedComparison) :
    EligibilityAntecedent A.armingInterface :=
  eligibilityAntecedent_of_interface A.armingInterface

/-! ## Composed receipt -/

/-- Composed receipt: the coefficient identities and control support hold at
every supplied positive scale, the placeholder schema admits distinct scale
and lower-bound labels, the ratios alone select no scale, and the comparison
record is definitionally a pair of placeholders.  No evidential arming
conclusion follows. -/
theorem dispersionArmingInterface_receipt :
    (∀ a : ℝ, a ≠ 0 →
      (scaleB0 a / scaleC4 a ^ 2 = 10 / 21 ∧
        edgeScaleB6 a / scaleC4 a ^ 2 = -2 / 63 ∧
        edgeScaleB6 a / scaleB0 a = -1 / 15)) ∧
    EdgeOrbitSupportCertificate ∧
    (∃ I₁ I₂ : DispersionArmingInterface,
      I₁.exclusivityLowerBound = I₂.exclusivityLowerBound ∧
        I₁.carrierScale ≠ I₂.carrierScale) ∧
    (∃ I₁ I₂ : DispersionArmingInterface,
      I₁.carrierScale = I₂.carrierScale ∧
        I₁.exclusivityLowerBound ≠ I₂.exclusivityLowerBound) ∧
    (¬ ∃! a : ℝ, 0 < a ∧
      (scaleB0 a / scaleC4 a ^ 2 = 10 / 21 ∧
        edgeScaleB6 a / scaleC4 a ^ 2 = -2 / 63 ∧
        edgeScaleB6 a / scaleB0 a = -1 / 15)) ∧
    (Nonempty ArmedComparison ↔
      Nonempty DispersionArmingInterface ∧
        Nonempty PostCustodyDatasetContract) :=
  ⟨edgeScale_ray, source_seam_edge30_control_certificate,
    committed_corpus_does_not_determine_scale,
    committed_corpus_does_not_determine_lower_bound,
    committed_ray_selects_no_scale,
    armedComparison_constructible_iff⟩

end OPH.DispersionArming

/- Axiom audit: propext, Classical.choice, and Quot.sound only. -/

#print axioms OPH.DispersionArming.seam_row_values
#print axioms OPH.DispersionArming.edgeScale_matches_seam_row
#print axioms OPH.DispersionArming.edgeScale_ray
#print axioms OPH.DispersionArming.scale_parameter_faithful
#print axioms OPH.DispersionArming.committed_ray_selects_no_scale
#print axioms OPH.DispersionArming.edge_vertex_rank_six_sign_separation
#print axioms OPH.DispersionArming.stipulatedInterface_carrierScale
#print axioms OPH.DispersionArming.stipulatedInterface_exclusivityLowerBound
#print axioms OPH.DispersionArming.dispersionArmingInterface_nonempty
#print axioms OPH.DispersionArming.interface_scale_free_parameter
#print axioms OPH.DispersionArming.interface_lower_bound_free_parameter
#print axioms OPH.DispersionArming.committed_corpus_does_not_determine_scale
#print axioms OPH.DispersionArming.committed_corpus_does_not_determine_lower_bound
#print axioms OPH.DispersionArming.committed_corpus_selects_no_lower_bound
#print axioms OPH.DispersionArming.prop_fields_carry_no_separation
#print axioms OPH.DispersionArming.contract_digest_free_parameter
#print axioms OPH.DispersionArming.contract_distinct_stipulations
#print axioms OPH.DispersionArming.armedComparison_constructible_iff
#print axioms OPH.DispersionArming.eligibilityAntecedent_of_interface
#print axioms OPH.DispersionArming.armed_eligibilityAntecedent
#print axioms OPH.DispersionArming.dispersionArmingInterface_receipt
