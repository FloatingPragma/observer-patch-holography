import Geometry.SourceOrderFrameCompatibilityPacket
import Geometry.EinsteinTensorBridge
import ObserverPatchHolography.EinsteinBranch.Composition

/-!
# Source-order to conditional Einstein-form composition

This module gives a genuinely shared-carrier composition between the finite
source-order/frame packet and the existing Einstein tensor algebra.  The
rank-three source quotient supplies an exact two-sphere of unit directions;
through the proved celestial/null-ray and Pauli-coordinate equivalences those
directions contain the nine null tomography directions already used by the
Einstein branch.  Consequently a balance law stated only on those nine exact
source-unit directions suffices for the finite null-tomography step.  Ward
and contracted Bianchi identities then make the remaining metric coefficient
constant, and the explicit vacuum reference and scale identification give the
displayed Einstein-form equation.

The distinguished nine-element family is the fixed algebraic inverse image
of the Einstein branch's pre-existing coordinate tomography frame.  It is
neither selected by the generated poset/source dynamics nor invariant, as a
distinguished set, under Lorentz or `SO(3)` transformations.  For each fixed
unit direction in that family, `(1,n)` is its canonical future-null
representative.  These directions are not observed provenance links, sampled
signals, or points on an identified physical sky.

The theorem keeps the source order and tensor fields on the same finite
`Event` carrier.  The tensor calculation uses the independent source Gram
quotient's directions, not the generated order; no order-to-field or
provenance-link-to-direction selection theorem is proved.  The source stack
now independently constructs a four-dimensional `ℝ × FrameQuotient` carrier
with `(+---)` signature.  This module does not yet identify that carrier's
finite differences with its supplied `V 3`, `Mat 3`, `eta 3`, or `Fin 4`
tensor/step data.  It also does not prove the source-direction balance, Ward
or Bianchi identities, connectivity, vacuum reference, Newton-scale
identification, physical causal faithfulness, count--volume calibration,
manifoldlikeness, a refinement/cofinal limit, smoothness, or a physical
spacetime interpretation.  In particular, its Einstein-form conclusion is a
conditional finite-carrier algebraic theorem, not a derivation of a smooth
Einstein manifold.  The premise bundle also carries `universalCoupling` into
the conclusion for downstream composition, but that equality is not used by
the tensor-equation proof.
-/

namespace OPH

noncomputable section

open C1Lorentz EinsteinBranch Provenance

universe u v w z

/-! ## Source directions in Einstein coordinates -/

/-- The canonical future-null Einstein-coordinate representative `(1,n)` of
a source-unit direction. -/
def sourceNullEinsteinVector (q : SourceUnitDirection) : V 3 :=
  pauliEinsteinChart
    (celestialRepresentative (sourceUnitDirectionEquivCelestial q))

/-- Every source direction's canonical `(1,n)` representative is null in the
Einstein branch's `(-+++)` coordinate convention. -/
theorem sourceNullEinsteinVector_null (q : SourceUnitDirection) :
    quadOf (eta 3) (sourceNullEinsteinVector q) = 0 := by
  exact (lorentzQ_eq_zero_iff_einsteinQuad_eq_zero _).mp
    (celestialRepresentative (sourceUnitDirectionEquivCelestial q)).2.1

/-- The canonical representative has future time coordinate one. -/
@[simp] theorem sourceNullEinsteinVector_time (q : SourceUnitDirection) :
    sourceNullEinsteinVector q 0 = 1 := by
  rfl

/-- A future-null-vector witness for each of the nine exact Einstein-branch
tomography directions. -/
def tomographyFutureNullVector (r : Fin 9) : FutureNullVector := by
  refine ⟨pauliEinsteinChart.symm (tomographyDirections r), ?_⟩
  apply (isFutureNull_iff_einstein _).2
  refine ⟨tomographyDirections_null r, ?_⟩
  fin_cases r <;> norm_num [tomographyDirections, pauliEinsteinChart]

/-- The source-unit direction corresponding to one of the nine exact null
tomography directions. -/
def sourceTomographyDirections (r : Fin 9) : SourceUnitDirection :=
  sourceUnitDirectionEquivCelestial.symm
    (normalizedDirection (tomographyFutureNullVector r))

/-- The canonical `(1,n)` representative of each fixed, coordinate-chosen
direction is exactly the pre-existing Einstein tomography vector, not merely
the same projective ray.  The nine-element family itself is not selected by
the source order. -/
theorem sourceNullEinsteinVector_sourceTomographyDirections (r : Fin 9) :
    sourceNullEinsteinVector (sourceTomographyDirections r) =
      tomographyDirections r := by
  funext i
  fin_cases r <;> fin_cases i <;>
    simp [sourceNullEinsteinVector, sourceTomographyDirections,
      tomographyFutureNullVector, normalizedDirection, tomographyDirections,
      sourceUnitDirectionEquivCelestial, celestialRepresentative,
      pauliEinsteinChart]

/-! ## Same-event-type source-indexed Einstein premises

The tensor interface below is supplied separately from the constructed
source-order carrier. Sharing the finite `Event` type supplies indexing only;
it does not identify tensor coordinates or differences with carrier geometry.
-/

/-- Minimal premises for the finite Einstein *shape* determined by the nine
exact source-unit tomography directions.  This package deliberately has no
phantom dependency on a source-order packet and omits the fields that only
name or physically normalize the final constants: entropy stress,
universal-source equality, Newton's constant, a reference cosmological
constant, vacuum calibration, and `8πG` scale identification. -/
structure SourceDirectionEinsteinShapePremises (Event : Type w)
    [Fintype Event] where
  step : Fin 4 → Event → Event
  base : Event
  geometry : Event → Mat 3
  stress : Event → Mat 3
  coupling : ℝ
  geometry_symmetric : ∀ p i j, geometry p i j = geometry p j i
  stress_symmetric : ∀ p i j, stress p i j = stress p j i
  /-- Balance on the nine fixed algebraic inverse images of the null
  tomography frame.  It remains supplied and is not selected by the poset. -/
  sourceTomographyBalance : ∀ p r,
    quadOf (geometry p)
        (sourceNullEinsteinVector (sourceTomographyDirections r)) =
      coupling * quadOf (stress p)
        (sourceNullEinsteinVector (sourceTomographyDirections r))
  ward : ∀ p j, ddiv step stress p j = 0
  bianchi : ∀ p j, ddiv step geometry p j = 0
  connected : ∀ q, SymmReachable step base q

namespace SourceDirectionEinsteinShapePremises

variable {Event : Type w} [Fintype Event]
variable (A : SourceDirectionEinsteinShapePremises Event)

/-- The minimal nine source-direction balances are exactly the coordinate
tomography balances. -/
theorem tomographyBalance (p : Event) (r : Fin 9) :
    quadOf (A.geometry p) (tomographyDirections r) =
      A.coupling * quadOf (A.stress p) (tomographyDirections r) := by
  rw [← sourceNullEinsteinVector_sourceTomographyDirections r]
  exact A.sourceTomographyBalance p r

/-- The nine source-direction balances determine each symmetric tensor
difference modulo one metric coefficient. -/
theorem pointwiseMetricAmbiguity :
    ∀ p, ∃ lam : ℝ, ∀ i j,
      A.geometry p i j =
        A.coupling * A.stress p i j + lam * eta 3 i j := by
  intro p
  let B : Mat 3 := fun i j ↦ A.geometry p i j - A.coupling * A.stress p i j
  have hBsymm : ∀ i j, B i j = B j i := by
    intro i j
    dsimp [B]
    rw [A.geometry_symmetric p i j, A.stress_symmetric p i j]
  have hBzero : ∀ r, quadOf B (tomographyDirections r) = 0 := by
    intro r
    dsimp [B]
    rw [quadOf_sub_smul, A.tomographyBalance p r, sub_self]
  refine ⟨metricCoefficient B, fun i j ↦ ?_⟩
  have hij := nine_null_directions_determine_mod_metric B hBsymm hBzero i j
  dsimp [B] at hij
  linarith

/-- Nine source-direction balances recover the all-null balance without an
additional premise. -/
theorem nullBalance (p : Event) (k : V 3)
    (hk : quadOf (eta 3) k = 0) :
    quadOf (A.geometry p) k =
      A.coupling * quadOf (A.stress p) k := by
  obtain ⟨lam, hlam⟩ := A.pointwiseMetricAmbiguity p
  have hmatrix :
      (fun i j ↦ A.geometry p i j - A.coupling * A.stress p i j) =
        fun i j ↦ lam * eta 3 i j := by
    funext i j
    rw [hlam i j]
    ring
  have hzero :
      quadOf (fun i j ↦
        A.geometry p i j - A.coupling * A.stress p i j) k = 0 := by
    rw [hmatrix]
    have hquadzero :
        quadOf (fun _ _ ↦ (0 : ℝ)) k = 0 := by
      simp [quadOf, bilinOf]
    have hscale :
        quadOf (fun i j ↦ lam * eta 3 i j) k =
          lam * quadOf (eta 3) k := by
      have h := quadOf_sub_smul
        (fun _ _ ↦ (0 : ℝ)) (eta 3) (-lam) k
      rw [hquadzero] at h
      simpa only [zero_sub, neg_mul, neg_neg] using h
    rw [hscale, hk, mul_zero]
  have hsub := quadOf_sub_smul
    (A.geometry p) (A.stress p) A.coupling k
  rw [hzero] at hsub
  linarith

/-- Ward and Bianchi conservation make the pointwise metric coefficient one
constant on the connected finite carrier.  This is the Einstein-form shape;
its physical interpretation and continuum limit remain separate. -/
theorem constantMetricAmbiguity :
    ∃ Lambda : ℝ, ∀ p i j,
      A.geometry p i j =
        A.coupling * A.stress p i j + Lambda * eta 3 i j := by
  choose lam hlam using A.pointwiseMetricAmbiguity
  obtain ⟨Lambda, hLambda⟩ := lambda_constant_symm A.step A.base A.connected
    A.geometry A.stress lam A.coupling hlam A.bianchi A.ward
  exact ⟨Lambda, fun p i j ↦ by rw [hlam p i j, hLambda p]⟩

end SourceDirectionEinsteinShapePremises

/-- Explicit finite-carrier premises for the source-indexed Einstein
composition.  This is the existing continuum/physical premise package with
its all-null `nullBalance` field replaced by balance on the exact unit
tomography directions of the source Gram quotient. -/
structure SourceIndexedEinsteinPremises
    {Register : Type u} {Value : Type v} {Event : Type w} {Chart : Type z}
    [DecidableEq Register] [Fintype Event]
    (P : SourceOrderFrameCompatibilityPacket Register Value Event Chart) where
  step : Fin 4 → Event → Event
  base : Event
  geometry : Event → Mat 3
  stress : Event → Mat 3
  entropyStress : Event → Mat 3
  coupling : ℝ
  newton : ℝ
  referenceLambda : ℝ
  geometry_symmetric : ∀ p i j, geometry p i j = geometry p j i
  stress_symmetric : ∀ p i j, stress p i j = stress p j i
  /-- Balance on the nine exact source directions inverse-mapped from the
  Einstein branch's null-tomography frame.  This field quantifies over only
  those nine fixed directions rather than all null or all source directions,
  and it is not proved from provenance. -/
  sourceTomographyBalance : ∀ p r,
    quadOf (geometry p)
        (sourceNullEinsteinVector (sourceTomographyDirections r)) =
      coupling * quadOf (stress p)
        (sourceNullEinsteinVector (sourceTomographyDirections r))
  ward : ∀ p j, ddiv step stress p j = 0
  bianchi : ∀ p j, ddiv step geometry p j = 0
  connected : ∀ q, SymmReachable step base q
  universalCoupling : entropyStress = stress
  vacuumReference :
    geometry base 0 0 =
      coupling * stress base 0 0 + referenceLambda * eta 3 0 0
  physicalScale : coupling = 8 * Real.pi * newton

namespace SourceIndexedEinsteinPremises

variable {Register : Type u} {Value : Type v} {Event : Type w}
  {Chart : Type z}
variable [DecidableEq Register] [Fintype Event]
variable {P : SourceOrderFrameCompatibilityPacket Register Value Event Chart}
variable (A : SourceIndexedEinsteinPremises P)

/-- Source-unit balance contains the nine exact null-tomography equalities. -/
theorem tomographyBalance (p : Event) (r : Fin 9) :
    quadOf (A.geometry p) (tomographyDirections r) =
      A.coupling * quadOf (A.stress p) (tomographyDirections r) := by
  rw [← sourceNullEinsteinVector_sourceTomographyDirections r]
  exact A.sourceTomographyBalance p r

/-- At each event, the nine source-indexed balances determine the geometric
tensor modulo one multiple of the Minkowski metric. -/
theorem pointwiseMetricAmbiguity :
    ∀ p, ∃ lam : ℝ, ∀ i j,
      A.geometry p i j =
        A.coupling * A.stress p i j + lam * eta 3 i j := by
  intro p
  let B : Mat 3 := fun i j ↦ A.geometry p i j - A.coupling * A.stress p i j
  have hBsymm : ∀ i j, B i j = B j i := by
    intro i j
    dsimp [B]
    rw [A.geometry_symmetric p i j, A.stress_symmetric p i j]
  have hBzero : ∀ r, quadOf B (tomographyDirections r) = 0 := by
    intro r
    dsimp [B]
    rw [quadOf_sub_smul, A.tomographyBalance p r, sub_self]
  refine ⟨metricCoefficient B, fun i j ↦ ?_⟩
  have hij := nine_null_directions_determine_mod_metric B hBsymm hBzero i j
  dsimp [B] at hij
  linarith

/-- The nine source-direction balances recover the all-null balance consumed
by the existing Jacobson/Ward/Bianchi composition.  This is a theorem, not an
additional premise. -/
theorem nullBalance (p : Event) (k : V 3)
    (hk : quadOf (eta 3) k = 0) :
    quadOf (A.geometry p) k =
      A.coupling * quadOf (A.stress p) k := by
  obtain ⟨lam, hlam⟩ := A.pointwiseMetricAmbiguity p
  have hmatrix :
      (fun i j ↦ A.geometry p i j - A.coupling * A.stress p i j) =
        fun i j ↦ lam * eta 3 i j := by
    funext i j
    rw [hlam i j]
    ring
  have hzero :
      quadOf (fun i j ↦
        A.geometry p i j - A.coupling * A.stress p i j) k = 0 := by
    rw [hmatrix]
    have hquadzero :
        quadOf (fun _ _ ↦ (0 : ℝ)) k = 0 := by
      simp [quadOf, bilinOf]
    have hscale :
        quadOf (fun i j ↦ lam * eta 3 i j) k =
          lam * quadOf (eta 3) k := by
      have h := quadOf_sub_smul
        (fun _ _ ↦ (0 : ℝ)) (eta 3) (-lam) k
      rw [hquadzero] at h
      simpa only [zero_sub, neg_mul, neg_neg] using h
    rw [hscale, hk, mul_zero]
  have hsub := quadOf_sub_smul
    (A.geometry p) (A.stress p) A.coupling k
  rw [hzero] at hsub
  linarith

/-- Adapter to the existing explicit Einstein-branch premise type.  Its
all-null field is filled by `nullBalance`, derived from the nine source
directions; the remaining physical and conservation fields stay explicit. -/
def toContinuumEinsteinPremises : ContinuumEinsteinPremises Event where
  step := A.step
  base := A.base
  geometry := A.geometry
  stress := A.stress
  entropyStress := A.entropyStress
  coupling := A.coupling
  newton := A.newton
  referenceLambda := A.referenceLambda
  geometry_symmetric := A.geometry_symmetric
  stress_symmetric := A.stress_symmetric
  nullBalance := A.nullBalance
  ward := A.ward
  bianchi := A.bianchi
  connected := A.connected
  universalCoupling := A.universalCoupling
  vacuumReference := A.vacuumReference
  physicalScale := A.physicalScale

/-- Ward+Bianchi conservation promotes the pointwise metric ambiguity to one
constant on the connected same-event carrier. -/
theorem constantMetricAmbiguity :
    ∃ Lambda : ℝ, ∀ p i j,
      A.geometry p i j =
        A.coupling * A.stress p i j + Lambda * eta 3 i j := by
  choose lam hlam using A.pointwiseMetricAmbiguity
  obtain ⟨Lambda, hLambda⟩ := lambda_constant_symm A.step A.base A.connected
    A.geometry A.stress lam A.coupling hlam A.bianchi A.ward
  exact ⟨Lambda, fun p i j ↦ by rw [hlam p i j, hLambda p]⟩

/-- The source-indexed balance, conservation, reference, and scale premises
give the Einstein-form equation on the same event carrier. -/
theorem einsteinEquation :
    ∀ p i j,
      A.geometry p i j =
        8 * Real.pi * A.newton * A.stress p i j +
          A.referenceLambda * eta 3 i j :=
  continuumEinstein_from_explicit_premises A.toContinuumEinsteinPremises

end SourceIndexedEinsteinPremises

/-! ## Composed conclusion -/

/-- Exact minimal output obtained by adjoining the source-order/frame
consequences to the independent source-direction Einstein-shape theorem.
The source order and tensor fields share `Event`, but the order is not used
to prove the tensor identity; that missing order-to-field bridge remains
visible rather than hidden in a phantom premise dependency. -/
structure SourceOrderEinsteinShapeConclusion
    {Register : Type u} {Value : Type v} {Event : Type w} {Chart : Type z}
    [DecidableEq Register] [Fintype Event]
    (P : SourceOrderFrameCompatibilityPacket Register Value Event Chart)
    (A : SourceDirectionEinsteinShapePremises Event) : Prop where
  finiteSourceOrderFrame : P.FiniteConsequences
  sourceTomographyBalance : ∀ p r,
    quadOf (A.geometry p)
        (sourceNullEinsteinVector (sourceTomographyDirections r)) =
      A.coupling * quadOf (A.stress p)
        (sourceNullEinsteinVector (sourceTomographyDirections r))
  allNullBalance : ∀ p k, quadOf (eta 3) k = 0 →
    quadOf (A.geometry p) k = A.coupling * quadOf (A.stress p) k
  pointwiseMetricAmbiguity : ∀ p, ∃ lam : ℝ, ∀ i j,
    A.geometry p i j =
      A.coupling * A.stress p i j + lam * eta 3 i j
  constantMetricAmbiguity : ∃ Lambda : ℝ, ∀ p i j,
    A.geometry p i j =
      A.coupling * A.stress p i j + Lambda * eta 3 i j

/-- Minimal same-event-type composition in a separately supplied `3+1` tensor
interface. It removes six non-load-bearing normalization/unification fields
from the Einstein-shape derivation while making explicit that a source-order
packet and source-direction tensor data are composed, not that the former
constructs the latter. -/
theorem sourceOrderEinsteinShape_from_minimal_source_directions
    {Register : Type u} {Value : Type v} {Event : Type w} {Chart : Type z}
    [DecidableEq Register] [Fintype Event]
    (P : SourceOrderFrameCompatibilityPacket Register Value Event Chart)
    (A : SourceDirectionEinsteinShapePremises Event) :
    SourceOrderEinsteinShapeConclusion P A := by
  exact ⟨P.finiteConsequences, A.sourceTomographyBalance,
    A.nullBalance, A.pointwiseMetricAmbiguity, A.constantMetricAmbiguity⟩

/-- Exact output of the same-event-type source-order/Einstein composition.
The tensor interface remains separately supplied. -/
structure SourceOrderEinsteinConclusion
    {Register : Type u} {Value : Type v} {Event : Type w} {Chart : Type z}
    [DecidableEq Register] [Fintype Event]
    (P : SourceOrderFrameCompatibilityPacket Register Value Event Chart)
    (A : SourceIndexedEinsteinPremises P) : Prop where
  finiteSourceOrderFrame : P.FiniteConsequences
  sourceTomographyBalance : ∀ p r,
    quadOf (A.geometry p)
        (sourceNullEinsteinVector (sourceTomographyDirections r)) =
      A.coupling * quadOf (A.stress p)
        (sourceNullEinsteinVector (sourceTomographyDirections r))
  allNullBalance : ∀ p k, quadOf (eta 3) k = 0 →
    quadOf (A.geometry p) k = A.coupling * quadOf (A.stress p) k
  pointwiseMetricAmbiguity : ∀ p, ∃ lam : ℝ, ∀ i j,
    A.geometry p i j =
      A.coupling * A.stress p i j + lam * eta 3 i j
  constantMetricAmbiguity : ∃ Lambda : ℝ, ∀ p i j,
    A.geometry p i j =
      A.coupling * A.stress p i j + Lambda * eta 3 i j
  universalSource : A.entropyStress = A.stress
  einstein : ∀ p i j,
    A.geometry p i j =
      8 * Real.pi * A.newton * A.stress p i j +
        A.referenceLambda * eta 3 i j

/-- Main composition theorem.  Unlike a conjunction of unrelated source and
Einstein packages, the proof consumes the source quotient's unit-direction
space to obtain the nine null-tomography equations used in the tensor step,
and all fields live on the source packet's finite event carrier. -/
theorem sourceOrderEinstein_from_source_directions
    {Register : Type u} {Value : Type v} {Event : Type w} {Chart : Type z}
    [DecidableEq Register] [Fintype Event]
    (P : SourceOrderFrameCompatibilityPacket Register Value Event Chart)
    (A : SourceIndexedEinsteinPremises P) :
    SourceOrderEinsteinConclusion P A := by
  exact ⟨P.finiteConsequences, A.sourceTomographyBalance,
    A.nullBalance, A.pointwiseMetricAmbiguity, A.constantMetricAmbiguity,
    A.universalCoupling, A.einsteinEquation⟩

/-! ## Per-theorem axiom audit -/

#print axioms sourceNullEinsteinVector_null
#print axioms sourceNullEinsteinVector_time
#print axioms sourceNullEinsteinVector_sourceTomographyDirections
#print axioms SourceDirectionEinsteinShapePremises.tomographyBalance
#print axioms SourceDirectionEinsteinShapePremises.pointwiseMetricAmbiguity
#print axioms SourceDirectionEinsteinShapePremises.nullBalance
#print axioms SourceDirectionEinsteinShapePremises.constantMetricAmbiguity
#print axioms sourceOrderEinsteinShape_from_minimal_source_directions
#print axioms SourceIndexedEinsteinPremises.tomographyBalance
#print axioms SourceIndexedEinsteinPremises.pointwiseMetricAmbiguity
#print axioms SourceIndexedEinsteinPremises.nullBalance
#print axioms SourceIndexedEinsteinPremises.toContinuumEinsteinPremises
#print axioms SourceIndexedEinsteinPremises.constantMetricAmbiguity
#print axioms SourceIndexedEinsteinPremises.einsteinEquation
#print axioms sourceOrderEinstein_from_source_directions

end

end OPH
