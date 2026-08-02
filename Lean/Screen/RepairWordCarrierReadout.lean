import PrimitivePortFrameQuotient
import ObserverPatchHolography.A2EndpointCommutator

open scoped Topology

namespace OPH.RepairWordCarrierReadout

open ObserverPatchHolography.A2EndpointCommutator
open OPH.PrimitivePortFrameQuotient

/-!
# Universal repair words and the response-selected carrier

The conditional A2 endpoint theorem shows that accepted signed port words
factor through the free abelian group on the six antipodal axes when its
source steps, inverse rows, and endpoint diamonds are supplied. The current
simulator does not emit that complete word-action packet. The repair-response theorem
constructs a metric completion of the corresponding integer six-axis record
module.  This file proves that these are the same algebraic object before any
Cartesian coordinates are introduced.

The resulting universal record module has a dense isometric image in the
repair-selected three-dimensional carrier.  Its metric completion is
Euclidean three-space up to the canonical comparison equivalence for abstract
completions.  This is a local source-carrier theorem.  A physical position
interpretation still requires the source action to be faithful and A2-RC to
select this response topology for observer readback.  Global overlap and
refinement gluing, physical scale, time, and field attachment are not proved
here.
-/

/-- Exact additive equivalence between the universal six-axis record group
and the concrete integer-control coordinates used by the Gram completion. -/
noncomputable def recordControlEquiv :
    FreeAbelianGroup Axis ≃+ (Axis → ℤ) :=
  (FreeAbelianGroup.equivFinsupp Axis).trans
    (Finsupp.linearEquivFunOnFinite ℤ ℤ Axis).toAddEquiv

/-- The universal abelian record as the exact integer-control wrapper consumed
by the response-selected carrier construction. -/
noncomputable def recordPoint (r : FreeAbelianGroup Axis) : IntegerFramePoint :=
  ⟨recordControlEquiv r⟩

/-- Recover the universal record from its integer-control wrapper. -/
noncomputable def pointRecord (p : IntegerFramePoint) : FreeAbelianGroup Axis :=
  recordControlEquiv.symm p.control

@[simp]
theorem pointRecord_recordPoint (r : FreeAbelianGroup Axis) :
    pointRecord (recordPoint r) = r := by
  simp [pointRecord, recordPoint]

@[simp]
theorem recordPoint_pointRecord (p : IntegerFramePoint) :
    recordPoint (pointRecord p) = p := by
  ext i
  simp [pointRecord, recordPoint]

/-- The universal record group and the operational integer wrapper are
bijective before either is equipped with the response metric. -/
noncomputable def recordPointEquiv :
    FreeAbelianGroup Axis ≃ IntegerFramePoint where
  toFun := recordPoint
  invFun := pointRecord
  left_inv := pointRecord_recordPoint
  right_inv := recordPoint_pointRecord

theorem recordPoint_bijective : Function.Bijective recordPoint :=
  recordPointEquiv.bijective

theorem recordPoint_injective : Function.Injective recordPoint :=
  recordPoint_bijective.1

theorem recordPoint_surjective : Function.Surjective recordPoint :=
  recordPoint_bijective.2

/-! ## Source port loads and the signed quotient

The bounded repair source acts on twelve integer port loads, while its
primitive attempt alphabet consists of thirty seams.  Those seam attempts are
not relabeled as translations here.  Instead, the source load itself has a
canonical antipodal-odd readback.  The next theorems show that the response
carrier identifies exactly the antipodal-even integer loads and no others.
-/

/-- Integer working load on the twelve registered source ports. -/
abbrev PortLoad := Fin 12 → ℤ

/-- The antipodal-odd load readback in the six positive port coordinates. -/
def signedLoadControl (x : PortLoad) : Axis → ℤ :=
  fun i ↦ x (positivePort i) -
    x (OPH.PortFrameGram.antipode (positivePort i))

/-- Explicit section of the signed-load readback.  It places the six controls
on the registered positive representatives and zero on their antipodes. -/
def controlLoad (z : Axis → ℤ) : PortLoad :=
  ![0, 0, 0, 0, 0, 0, z 4, z 2, z 0, z 5, z 1, z 3]

theorem signedLoadControl_controlLoad (z : Axis → ℤ) :
    signedLoadControl (controlLoad z) = z := by
  funext i
  fin_cases i <;>
    simp [signedLoadControl, controlLoad, positivePort, positiveSourcePort,
      sourceToRERPort, OPH.PortFrameGram.antipode]

/-- The response-carrier integer point read directly from a twelve-port load. -/
noncomputable def loadPoint (x : PortLoad) : IntegerFramePoint :=
  ⟨signedLoadControl x⟩

theorem loadPoint_surjective : Function.Surjective loadPoint := by
  intro p
  refine ⟨controlLoad p.control, ?_⟩
  ext i
  exact congrFun (signedLoadControl_controlLoad p.control) i

/-- An integer port load is antipodally even when every positive/negative
pair has equal load. -/
def AntipodalEven (x : PortLoad) : Prop :=
  ∀ i : Axis,
    x (positivePort i) =
      x (OPH.PortFrameGram.antipode (positivePort i))

theorem signedLoadControl_eq_zero_iff (x : PortLoad) :
    signedLoadControl x = 0 ↔ AntipodalEven x := by
  constructor
  · intro h i
    have hi := congrFun h i
    simp only [signedLoadControl, Pi.zero_apply] at hi
    omega
  · intro h
    funext i
    simp [signedLoadControl, h i]

theorem signedLoadControl_sub (x y : PortLoad) :
    signedLoadControl (x - y) =
      signedLoadControl x - signedLoadControl y := by
  funext i
  simp [signedLoadControl]
  ring

/-- Two integer source loads define the same carrier record exactly when
their difference is antipodally even.  Thus the response quotient introduces
no additional integral load relation. -/
theorem loadPoint_eq_iff_antipodalEven_difference (x y : PortLoad) :
    loadPoint x = loadPoint y ↔ AntipodalEven (x - y) := by
  constructor
  · intro h
    have hc : signedLoadControl x = signedLoadControl y := by
      exact congrArg IntegerFramePoint.control h
    rw [← signedLoadControl_eq_zero_iff, signedLoadControl_sub, hc, sub_self]
  · intro h
    apply IntegerFramePoint.ext
    have hz : signedLoadControl (x - y) = 0 :=
      (signedLoadControl_eq_zero_iff (x - y)).2 h
    rw [signedLoadControl_sub] at hz
    exact sub_eq_zero.mp hz

/-- The twelve-port load quotient has dense image in the same
repair-selected carrier.  This statement ranges over the full integer load
module; a chosen protected-total sector is a separate affine subspace. -/
noncomputable def loadPosition (x : PortLoad) : EuclideanVec3 :=
  pointEuclideanFrame (loadPoint x)

theorem loadPosition_denseRange : DenseRange loadPosition := by
  apply pointEuclideanFrame_denseRange.mono
  rintro y ⟨p, rfl⟩
  obtain ⟨x, rfl⟩ := loadPoint_surjective p
  exact ⟨x, rfl⟩

/-- A word in the conditional six-axis endpoint control has the carrier record
supplied by its exponent vector. This control alphabet is not the canonical
thirty-seam repair alphabet. Order is retained by the word itself; only the
cumulative carrier coordinate factors through this map. -/
noncomputable def wordPoint (word : List SignedPort) : IntegerFramePoint :=
  recordPoint (exponentVector word)

@[simp]
theorem wordPoint_forward_reverse (i : Axis) :
    wordPoint [forward i, reverse i] = recordPoint 0 := by
  simp [wordPoint]

/-- Equal A2 exponent vectors give identical carrier records. -/
theorem wordPoint_eq_of_exponentVector_eq
    {left right : List SignedPort}
    (h : exponentVector left = exponentVector right) :
    wordPoint left = wordPoint right := by
  simp only [wordPoint, h]

/-- The cumulative carrier record forgets exactly the abelian word relation
and no additional load relation. -/
theorem wordPoint_eq_iff_exponentVector_eq
    (left right : List SignedPort) :
    wordPoint left = wordPoint right ↔
      exponentVector left = exponentVector right := by
  constructor
  · intro h
    apply recordPoint_injective
    exact h
  · exact wordPoint_eq_of_exponentVector_eq

/-- Wrapper that equips the universal free-abelian record itself with the
repair-response metric.  A wrapper prevents this metric from being confused
with a discrete word or record-count metric. -/
@[ext]
structure UniversalRecordPoint where
  record : FreeAbelianGroup Axis

/-- Exact bijection from the universal record wrapper to the existing
integer-control wrapper. -/
noncomputable def universalToInteger
    (p : UniversalRecordPoint) : IntegerFramePoint :=
  recordPoint p.record

theorem universalToInteger_injective :
    Function.Injective universalToInteger := by
  intro p q h
  ext
  exact recordPoint_injective h

theorem universalToInteger_surjective :
    Function.Surjective universalToInteger := by
  intro p
  obtain ⟨r, rfl⟩ := recordPoint_surjective p
  exact ⟨⟨r⟩, rfl⟩

/-- Pull back the source-Gram metric from the exact integer-control wrapper. -/
noncomputable instance universalRecordMetric : MetricSpace UniversalRecordPoint :=
  MetricSpace.induced universalToInteger universalToInteger_injective inferInstance

theorem universalToInteger_isometry : Isometry universalToInteger :=
  MetricSpace.isometry_induced universalToInteger universalToInteger_injective

/-- Intrinsic universal-record position in the repair-selected Euclidean
carrier.  The coordinates present the completion; they do not select it. -/
noncomputable def universalPosition
    (p : UniversalRecordPoint) : EuclideanVec3 :=
  pointEuclideanFrame (universalToInteger p)

theorem universalPosition_isometry : Isometry universalPosition :=
  pointEuclideanFrame_isometry.comp universalToInteger_isometry

/-- Every integer-control record is represented by a universal A2 record, so
the universal-record image inherits the proved dense range in the selected
carrier. -/
theorem universalPosition_denseRange : DenseRange universalPosition := by
  apply pointEuclideanFrame_denseRange.mono
  rintro y ⟨p, rfl⟩
  obtain ⟨r, rfl⟩ := universalToInteger_surjective p
  exact ⟨r, rfl⟩

/-- Euclidean three-space is an abstract metric completion of the universal
free-abelian repair-record module equipped with the response-selected Gram
metric. -/
noncomputable def universalRecordCompletion :
    AbstractCompletion UniversalRecordPoint where
  space := EuclideanVec3
  coe := universalPosition
  uniformStruct := inferInstance
  complete := inferInstance
  separation := inferInstance
  isUniformInducing := universalPosition_isometry.isUniformInducing
  dense := universalPosition_denseRange

/-- Canonical comparison equivalence from the standard metric completion of
the universal record module to the repair-selected Euclidean carrier. -/
noncomputable def universalCompletionEquivEuclidean3 :
    UniformSpace.Completion UniversalRecordPoint ≃ᵤ EuclideanVec3 :=
  AbstractCompletion.compareEquiv UniformSpace.Completion.cPkg
    universalRecordCompletion

/-- If the descended observer action represents the universal record group
faithfully, equality of visible word actions is exactly equality of the
response-carrier records.  Faithfulness is the remaining semantic gate; it is
not assumed elsewhere by this file. -/
theorem wordPoint_eq_iff_evalWord_eq_of_faithful
    {G : Type*} [Group G] (g : Axis → G)
    (φ : FreeAbelianGroup Axis →+ Additive (generatedSubgroup g))
    (hfaithful : Function.Injective φ)
    (hword : ∀ word,
      Additive.toMul (φ (exponentVector word)) = evalWord g word)
    (left right : List SignedPort) :
    wordPoint left = wordPoint right ↔
      evalWord g left = evalWord g right := by
  constructor
  · intro hp
    have he : exponentVector left = exponentVector right := by
      apply recordPoint_injective
      exact hp
    calc
      evalWord g left = Additive.toMul (φ (exponentVector left)) :=
        (hword left).symm
      _ = Additive.toMul (φ (exponentVector right)) := by rw [he]
      _ = evalWord g right := hword right
  · intro heval
    change recordPoint (exponentVector left) =
      recordPoint (exponentVector right)
    apply congrArg recordPoint
    apply hfaithful
    change Additive.toMul (φ (exponentVector left)) =
      Additive.toMul (φ (exponentVector right))
    calc
      Additive.toMul (φ (exponentVector left)) = evalWord g left := hword left
      _ = evalWord g right := heval
      _ = Additive.toMul (φ (exponentVector right)) := (hword right).symm

end OPH.RepairWordCarrierReadout

#print axioms OPH.RepairWordCarrierReadout.recordPoint_bijective
#print axioms OPH.RepairWordCarrierReadout.loadPoint_eq_iff_antipodalEven_difference
#print axioms OPH.RepairWordCarrierReadout.loadPosition_denseRange
#print axioms OPH.RepairWordCarrierReadout.universalPosition_denseRange
#print axioms OPH.RepairWordCarrierReadout.wordPoint_eq_iff_evalWord_eq_of_faithful
