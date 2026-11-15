import SeamCurrentCarrierQuotient
import A5OrbitRaySeparation

open scoped BigOperators goldenRatio

namespace OPH.SeamCurrentEdge30Moment

open OPH.PrimitivePortFrameQuotient
open OPH.SeamCurrentCarrierQuotient
open OPH.A5OrbitRaySeparation

/-!
# Exact edge-orbit moments of the complete seam-current support

The thirty source seams first enter the selected carrier through their exact
incidence boundary and the six-axis load quotient.  This file proves that the
resulting carrier differences have one common nonzero norm and that their
complete unoriented direction multiset is the edge-30 direction multiset,
including multiplicity.  It then computes the second, fourth, and sixth
normalized moments and derives the edge-30 coefficients and ratios used by
the control packet.

No seam is declared to be a physical hop.  Applying this control ray to a
physical dispersion symbol `Λ` requires two separately named premises: that
the seam differences act as physical translations and that observer readout
uses this completed carrier.  No clock, physical scale, polarization, or
laboratory field is supplied here.
-/

abbrev Vec3 := OPH.PrimitivePortFrameQuotient.Vec3

/-! ## Intrinsic seam differences and the edge-30 projective support -/

/-- The twelve registered carrier vertices in the same exact frame used by
the Gram completion.  This table is only a convenient expansion of
`rawAxis` and its antipodes; the next two theorems certify that binding. -/
noncomputable def portVector : Fin 12 → Vec3 :=
  ![
    ![0, -1, -φ],
    ![-1, -φ, 0],
    ![-φ, 0, -1],
    ![1, -φ, 0],
    ![0, 1, -φ],
    ![-φ, 0, 1],
    ![φ, 0, -1],
    ![0, -1, φ],
    ![-1, φ, 0],
    ![φ, 0, 1],
    ![1, φ, 0],
    ![0, 1, φ]
  ]

theorem portVector_positivePort (i : Fin 6) :
    portVector (positivePort i) = rawAxis i := by
  fin_cases i <;> rfl

theorem portVector_antipode (p : Fin 12) :
    portVector (OPH.PortFrameGram.antipode p) = -portVector p := by
  funext d
  fin_cases p <;> fin_cases d <;>
    simp [portVector, OPH.PortFrameGram.antipode]

/-- Carrier vector obtained intrinsically from a unit current on one source
seam, through the proved boundary and signed-load quotient. -/
noncomputable def carrierSeamDifference (e : Fin 30) : Vec3 :=
  integerFrame (seamAxisCurrent (seamAtom e 1))

/-- Direct endpoint difference in the exact carrier frame. -/
noncomputable def endpointDifference (e : Fin 30) : Vec3 :=
  portVector (seamRight e) - portVector (seamLeft e)

set_option maxHeartbeats 3000000 in
/-- The intrinsic current construction and the endpoint difference are the
same vector for all thirty seams. -/
theorem carrierSeamDifference_eq_endpointDifference (e : Fin 30) :
    carrierSeamDifference e = endpointDifference e := by
  rw [carrierSeamDifference, seamAxisCurrent_eq_table]
  funext d
  fin_cases e <;> fin_cases d <;>
    simp [endpointDifference, portVector, integerFrame, frameMap,
      castIntegerControl, seamAxisTable, seamAtom, rawAxis,
      seamLeft, seamRight, Fin.sum_univ_succ] <;>
    ring

set_option maxHeartbeats 2000000 in
/-- Every intrinsic seam difference has squared norm exactly four. -/
theorem carrierSeamDifference_norm_sq (e : Fin 30) :
    OPH.PrimitivePortTranslationBridge.dot
      (carrierSeamDifference e) (carrierSeamDifference e) = 4 := by
  rw [carrierSeamDifference_eq_endpointDifference]
  fin_cases e <;>
    simp [endpointDifference, portVector,
      OPH.PrimitivePortTranslationBridge.dot, seamLeft, seamRight,
      Fin.sum_univ_succ] <;>
    nlinarith [Real.goldenRatio_sq]

theorem carrierSeamDifference_ne_zero (e : Fin 30) :
    carrierSeamDifference e ≠ 0 := by
  intro h
  have hnorm := carrierSeamDifference_norm_sq e
  simp [h, OPH.PrimitivePortTranslationBridge.dot] at hnorm

/-- Midpoint direction of one incidence edge.  The thirty rows form the
usual edge-30 orbit before quotienting direction reversal. -/
noncomputable def edgeMidpoint (e : Fin 30) : Vec3 :=
  portVector (seamLeft e) + portVector (seamRight e)

/-- Fifteen representatives of the unoriented edge axes. -/
def representativeEdge : Fin 15 → Fin 30 :=
  ![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15]

/-- Projective class of each seam-difference direction. -/
def seamAxisClass : Fin 30 → Fin 15 :=
  ![14, 12, 9, 6, 5, 13, 10, 2, 1, 7, 3, 0, 8, 4, 0,
    11, 4, 1, 11, 8, 5, 3, 2, 7, 6, 10, 9, 13, 12, 14]

/-- Projective class of each edge-midpoint direction. -/
def edgeAxisClass : Fin 30 → Fin 15 :=
  ![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 11,
    14, 13, 8, 14, 12, 4, 10, 7, 9, 3, 6, 2, 5, 1, 0]

/-- Each projective seam class occurs exactly twice. -/
theorem seamAxisClass_fiber_two :
    ∀ i : Fin 15,
      (Finset.univ.filter (fun e : Fin 30 ↦ seamAxisClass e = i)).card = 2 := by
  decide

/-- Each projective edge-midpoint class occurs exactly twice. -/
theorem edgeAxisClass_fiber_two :
    ∀ i : Fin 15,
      (Finset.univ.filter (fun e : Fin 30 ↦ edgeAxisClass e = i)).card = 2 := by
  decide

/-- The two complete projective class multisets agree, including
multiplicity. -/
theorem projective_class_multisets_equal :
    (Finset.univ.val.map seamAxisClass : Multiset (Fin 15)) =
      Finset.univ.val.map edgeAxisClass := by
  decide

set_option maxHeartbeats 3000000 in
/-- Exact geometric realization of every seam projective class.  A common
factor `φ` relates the unit-current difference to the chosen midpoint
representative; the minus sign is discarded by projectivization. -/
theorem scaled_seamDifference_eq_neg_edgeRepresentative (e : Fin 30) :
    φ • carrierSeamDifference e =
      -edgeMidpoint (representativeEdge (seamAxisClass e)) := by
  have hsqrt : Real.sqrt 5 ^ 2 = 5 := by norm_num
  rw [carrierSeamDifference_eq_endpointDifference]
  funext d
  fin_cases e <;> fin_cases d <;>
    simp [endpointDifference, edgeMidpoint, portVector, representativeEdge,
      seamAxisClass, seamLeft, seamRight] <;>
    ring_nf <;>
    nlinarith [hsqrt]

set_option maxHeartbeats 3000000 in
/-- Exact realization of every row of the edge-midpoint orbit by its chosen
projective representative, with orientation deliberately forgotten. -/
theorem edgeMidpoint_eq_up_to_sign_representative (e : Fin 30) :
    edgeMidpoint e = edgeMidpoint (representativeEdge (edgeAxisClass e)) ∨
      edgeMidpoint e = -edgeMidpoint (representativeEdge (edgeAxisClass e)) := by
  fin_cases e <;>
    simp [edgeMidpoint, portVector, representativeEdge, edgeAxisClass,
      seamLeft, seamRight]

/-! ## Exact moment identities -/

private theorem sqrt5_pow_two : Real.sqrt 5 ^ 2 = (5 : ℝ) := by
  norm_num

private theorem sqrt5_pow_three : Real.sqrt 5 ^ 3 = 5 * Real.sqrt 5 := by
  calc
    Real.sqrt 5 ^ 3 = Real.sqrt 5 * Real.sqrt 5 ^ 2 := by ring
    _ = 5 * Real.sqrt 5 := by rw [sqrt5_pow_two]; ring

private theorem sqrt5_pow_four : Real.sqrt 5 ^ 4 = (25 : ℝ) := by
  calc
    Real.sqrt 5 ^ 4 = (Real.sqrt 5 ^ 2) ^ 2 := by ring
    _ = 25 := by rw [sqrt5_pow_two]; norm_num

private theorem sqrt5_pow_five : Real.sqrt 5 ^ 5 = 25 * Real.sqrt 5 := by
  calc
    Real.sqrt 5 ^ 5 = Real.sqrt 5 * Real.sqrt 5 ^ 4 := by ring
    _ = 25 * Real.sqrt 5 := by rw [sqrt5_pow_four]; ring

private theorem sqrt5_pow_six : Real.sqrt 5 ^ 6 = (125 : ℝ) := by
  calc
    Real.sqrt 5 ^ 6 = (Real.sqrt 5 ^ 2) ^ 3 := by ring
    _ = 125 := by rw [sqrt5_pow_two]; norm_num

/-- Squared radial coordinate in the selected carrier chart. -/
noncomputable def radiusSquared (k : Vec3) : ℝ :=
  OPH.PrimitivePortTranslationBridge.dot k k

/-- Normalized second moment of the complete intrinsic seam support. -/
noncomputable def seamMoment2 (k : Vec3) : ℝ :=
  ∑ e : Fin 30,
    OPH.PrimitivePortTranslationBridge.dot (carrierSeamDifference e) k ^ 2 / 4

/-- Normalized fourth moment of the complete intrinsic seam support. -/
noncomputable def seamMoment4 (k : Vec3) : ℝ :=
  ∑ e : Fin 30,
    OPH.PrimitivePortTranslationBridge.dot (carrierSeamDifference e) k ^ 4 / 16

/-- Normalized sixth moment of the complete intrinsic seam support. -/
noncomputable def seamMoment6 (k : Vec3) : ℝ :=
  ∑ e : Fin 30,
    OPH.PrimitivePortTranslationBridge.dot (carrierSeamDifference e) k ^ 6 / 64

/-- Icosahedral spin-six invariant, normalized to one on every unit vertex
direction.  This explicit polynomial is the exact expansion of the
vertex-orbit residual used by the finite moment certificate. -/
noncomputable def I6 (k : Vec3) : ℝ :=
  φ * (-105 * k 0 ^ 4 * k 1 ^ 2 / 16
      + 105 * k 0 ^ 4 * k 2 ^ 2 / 16
      + 105 * k 0 ^ 2 * k 1 ^ 4 / 16
      - 105 * k 0 ^ 2 * k 2 ^ 4 / 16
      - 105 * k 1 ^ 4 * k 2 ^ 2 / 16
      + 105 * k 1 ^ 2 * k 2 ^ 4 / 16)
    - 5 * k 0 ^ 6 / 16
    + 45 * k 0 ^ 4 * k 1 ^ 2 / 8
    - 15 * k 0 ^ 4 * k 2 ^ 2 / 16
    - 15 * k 0 ^ 2 * k 1 ^ 4 / 16
    - 225 * k 0 ^ 2 * k 1 ^ 2 * k 2 ^ 2 / 8
    + 45 * k 0 ^ 2 * k 2 ^ 4 / 8
    - 5 * k 1 ^ 6 / 16
    + 45 * k 1 ^ 4 * k 2 ^ 2 / 8
    - 15 * k 1 ^ 2 * k 2 ^ 4 / 16
    - 5 * k 2 ^ 6 / 16

set_option maxHeartbeats 3000000 in
/-- `I6` has the declared vertex normalization. -/
theorem I6_vertex_normalization (p : Fin 12) :
    I6 (portVector p) = radiusSquared (portVector p) ^ 3 := by
  fin_cases p <;>
    simp [I6, radiusSquared, portVector,
      OPH.PrimitivePortTranslationBridge.dot, Fin.sum_univ_succ] <;>
    ring_nf
  all_goals
    simp only [sqrt5_pow_six, sqrt5_pow_five, sqrt5_pow_four,
      sqrt5_pow_three, sqrt5_pow_two]
    ring

set_option maxHeartbeats 4000000 in
/-- Exact isotropic second moment. -/
theorem seamMoment2_eq (k : Vec3) :
    seamMoment2 k = 10 * radiusSquared k := by
  unfold seamMoment2
  simp_rw [carrierSeamDifference_eq_endpointDifference]
  simp [radiusSquared, endpointDifference, portVector, seamLeft, seamRight,
    OPH.PrimitivePortTranslationBridge.dot, Fin.sum_univ_succ]
  ring_nf
  simp only [sqrt5_pow_two]
  ring

set_option maxHeartbeats 5000000 in
/-- Exact isotropic fourth moment. -/
theorem seamMoment4_eq (k : Vec3) :
    seamMoment4 k = 6 * radiusSquared k ^ 2 := by
  unfold seamMoment4
  simp_rw [carrierSeamDifference_eq_endpointDifference]
  simp [radiusSquared, endpointDifference, portVector, seamLeft, seamRight,
    OPH.PrimitivePortTranslationBridge.dot, Fin.sum_univ_succ]
  ring_nf
  simp only [sqrt5_pow_four, sqrt5_pow_two]
  ring

set_option maxHeartbeats 8000000 in
/-- Exact sixth moment, split into the isotropic radial term and the unique
spin-six invariant fixed by the vertex normalization. -/
theorem seamMoment6_eq (k : Vec3) :
    seamMoment6 k = (30 / 7 : ℝ) * radiusSquared k ^ 3
      - (2 / 7 : ℝ) * I6 k := by
  unfold seamMoment6
  simp_rw [carrierSeamDifference_eq_endpointDifference]
  simp [radiusSquared, I6, endpointDifference, portVector, seamLeft, seamRight,
    OPH.PrimitivePortTranslationBridge.dot, Fin.sum_univ_succ]
  ring_nf
  simp only [sqrt5_pow_six, sqrt5_pow_five, sqrt5_pow_four,
    sqrt5_pow_three, sqrt5_pow_two]
  ring

/-! ## Derived edge-30 control ray and physical boundary -/

/-- Equal-weight normalization of the complete thirty-direction support. -/
def edge30Weight : ℚ := 1 / 5

/-- Fourth-order coefficient of the normalized control symbol `Λ`. -/
def seamLambdaC4 : ℚ := -edge30Weight * 6 / 24

/-- Isotropic sixth-order coefficient of `Λ`. -/
def seamLambdaB0 : ℚ := edge30Weight * (30 / 7) / 720

/-- Spin-six coefficient of `Λ`. -/
def seamLambdaB6 : ℚ := edge30Weight * (-2 / 7) / 720

/-- The coefficients derived from the seam moments are exactly the edge-30
row of the independently registered orbit table. -/
theorem seam_coefficients_eq_edge30 :
    seamLambdaC4 = C4 ∧
      seamLambdaB0 = B0 ∧
      seamLambdaB6 = B6 .edge30 := by
  norm_num [seamLambdaC4, seamLambdaB0, seamLambdaB6, edge30Weight,
    C4, B0, B6, OPH.A5PrimitivePortPrediction.C4,
    OPH.A5PrimitivePortPrediction.B0]

/-- Exact normalized edge-30 ray derived from the complete seam moments. -/
theorem seam_edge30_ray :
    seamLambdaB0 / seamLambdaC4 ^ 2 = 10 / 21 ∧
      seamLambdaB6 / seamLambdaC4 ^ 2 = -2 / 63 ∧
      seamLambdaB6 / seamLambdaB0 = -1 / 15 := by
  norm_num [seamLambdaC4, seamLambdaB0, seamLambdaB6, edge30Weight]

/-- Complete algebraic certificate tying the source-current support to the
edge-30 control packet. -/
theorem source_seam_edge30_control_certificate :
    (Fintype.card (Fin 30) = 30) ∧
    (∀ e : Fin 30,
      OPH.PrimitivePortTranslationBridge.dot
        (carrierSeamDifference e) (carrierSeamDifference e) = 4) ∧
    ((Finset.univ.val.map seamAxisClass : Multiset (Fin 15)) =
      Finset.univ.val.map edgeAxisClass) ∧
    seamLambdaC4 = C4 ∧ seamLambdaB0 = B0 ∧ seamLambdaB6 = B6 .edge30 := by
  exact ⟨by decide, carrierSeamDifference_norm_sq,
    projective_class_multisets_equal, seam_coefficients_eq_edge30⟩

/-- Physical use of the derived ray is conditional on two named attachment
premises.  Lean does not manufacture proofs of either premise. -/
def PhysicalEdge30Conclusion
    (seamDifferencesActAsPhysicalTranslations : Prop)
    (observerReadoutUsesCompletedCarrier : Prop) : Prop :=
  seamDifferencesActAsPhysicalTranslations ∧
    observerReadoutUsesCompletedCarrier ∧
    seamLambdaC4 = C4 ∧
    seamLambdaB0 = B0 ∧
    seamLambdaB6 = B6 .edge30

theorem physical_edge30_of_named_premises
    (seamDifferencesActAsPhysicalTranslations : Prop)
    (observerReadoutUsesCompletedCarrier : Prop)
    (htranslation : seamDifferencesActAsPhysicalTranslations)
    (hreadout : observerReadoutUsesCompletedCarrier) :
    PhysicalEdge30Conclusion seamDifferencesActAsPhysicalTranslations
      observerReadoutUsesCompletedCarrier := by
  rcases seam_coefficients_eq_edge30 with ⟨hC4, hB0, hB6⟩
  exact ⟨htranslation, hreadout, hC4, hB0, hB6⟩

end OPH.SeamCurrentEdge30Moment

/- Axiom audit: exact finite tables, real polynomial identities, and the
existing response-selected carrier construction only. -/

#print axioms OPH.SeamCurrentEdge30Moment.carrierSeamDifference_eq_endpointDifference
#print axioms OPH.SeamCurrentEdge30Moment.carrierSeamDifference_norm_sq
#print axioms OPH.SeamCurrentEdge30Moment.projective_class_multisets_equal
#print axioms OPH.SeamCurrentEdge30Moment.scaled_seamDifference_eq_neg_edgeRepresentative
#print axioms OPH.SeamCurrentEdge30Moment.seamMoment2_eq
#print axioms OPH.SeamCurrentEdge30Moment.seamMoment4_eq
#print axioms OPH.SeamCurrentEdge30Moment.seamMoment6_eq
#print axioms OPH.SeamCurrentEdge30Moment.source_seam_edge30_control_certificate
#print axioms OPH.SeamCurrentEdge30Moment.physical_edge30_of_named_premises
