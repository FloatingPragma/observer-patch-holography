import Geometry.SourceDerivedSpacetimeCarrier
import Geometry.SpatialReadbackSoldering

/-!
# Finite source-order/frame compatibility packet

This module composes three already separate finite ingredients:

* an authenticated semantic log whose reflexive generated order is attached,
  by a supplied source-native chart interface, to the declared future cone;
* the independently proved rank-three source Gram quotient; and
* one supplied oriented Lorentz transport per event, used to choose that
  event's internal observer frame and to transport the canonical
  `FrameQuotient`--`RestSpace standardFrame` bridge into its rest fiber.

The composition is deliberately noncircular.  The semantic log determines
the informational order; on a finite carrier its canonical longest-parent-
path height supplies the ordinal coordinate only inside an event placement.
It does not select the positive conversion scale, spatial event readback, or
frame transports. Conversely, the supplied frames do not add provenance
edges. An independent real axis and the exact rank-three source quotient
construct the algebraic carrier `ℝ × FrameQuotient`, its `(+---)` form,
and the source-unit future-null map. A spatial placement with source height,
the explicit scale, edge-speed,
and converse-support certificates constructs the one-chart cone attachment;
an arbitrary packet may still supply a larger finite atlas.  None of those
finite algebraic statements identifies a physical sky or signal relation.

**Nonclaims.**  This packet does not identify its algebraic unit sphere with
a physical celestial screen and proves no physical-causality theorem,
no volume or count-density law, no manifoldlikeness estimator, no
refinement/cofinal convergence theorem, and no topological, smooth-manifold,
or continuum-spacetime limit.  The exact finite order/cone embedding is
conditional on supplied spatial coordinates, an authenticated-edge speed
bound, and an explicit cone-support converse.  Event-frame transports and
any nontrivial multi-chart atlas remain supplied data.
-/

namespace OPH

noncomputable section

open C1Lorentz C2Soldering CausalComposition EventPopulation
  PrimitivePortFrameQuotient Provenance SourceDerivedSpacetime

universe u v w z

/-! ## Exact source-direction/celestial bridge -/

/-- Unit directions of the exact source Gram quotient.  This is the unit
sphere of the proved rank-three carrier, not yet a physical celestial
screen. -/
abbrev SourceUnitDirection :=
  SourceSpatialUnitDirection

/-- The unit sphere of the exact rank-three source quotient is algebraically
the C1 celestial two-sphere.  Both sides use the same proved Euclidean Gram
under `quotientEquivVec3`; no event population or spacetime is used. -/
noncomputable def sourceUnitDirectionEquivCelestial :
    SourceUnitDirection ≃ CelestialSphere where
  toFun q := by
    refine ⟨quotientEquivVec3 q.1, ?_⟩
    simpa [quotientGram, spatialNormSq,
      OPH.PrimitivePortTranslationBridge.dot, pow_two] using q.2
  invFun n := by
    refine ⟨quotientEquivVec3.symm n.1, ?_⟩
    simpa [quotientGram, spatialNormSq,
      OPH.PrimitivePortTranslationBridge.dot, pow_two] using n.2
  left_inv q := by
    apply Subtype.ext
    simp
  right_inv n := by
    apply Subtype.ext
    simp

/-- Equivalently, every source-unit direction labels exactly one intrinsic
future-null ray.  This uses the previously proved null-ray/celestial
equivalence and carries no assertion that a provenance edge is a light ray.
-/
noncomputable def sourceUnitDirectionEquivFutureNullRay :
    SourceUnitDirection ≃ FutureNullRay :=
  sourceUnitDirectionEquivCelestial.trans futureNullRayEquivCelestial.symm

/-- A finite source-order/chart compatibility package together with a supplied
oriented Lorentz transport for each event.  Each transport chooses the
event's internal frame as the image of `standardFrame`; it is not generated
by semantic provenance. -/
structure SourceOrderFrameCompatibilityPacket
    (Register : Type u) (Value : Type v) (Event : Type w) (Chart : Type z)
    [DecidableEq Register] [Fintype Event] where
  sourceCharts : SourceDerivedCausalChartInterface
    Register Value Event Chart
  standardToEventFrame : Event → OrientedLorentzEquiv

/-- An order-faithful canonical-height/positive-scale/spatial placement constructs
the active source chart packet directly, without the legacy rich-region
population map or a separately declared precedence relation. -/
noncomputable def SourceOrderFrameCompatibilityPacket.ofFaithfulPlacement
    {Register : Type u} {Value : Type v} {Event : Type w}
    [DecidableEq Register] [Fintype Event] [DecidableEq Event]
    [DecidableEq Value]
    {L : SemanticEventLog Register Value Event}
    (placement : FaithfulRankSpatialCausalPlacement L)
    (standardToEventFrame : Event → OrientedLorentzEquiv) :
    SourceOrderFrameCompatibilityPacket Register Value Event Unit where
  sourceCharts := placement.toSourceDerivedCausalChartInterface
  standardToEventFrame := standardToEventFrame

/-- Canonical algebraic specialization with the same standard frame at every
event.  This removes the frame-transport field when only existence of the
finite source-order/frame packet is required.  It is a trivial gauge choice,
not a derivation of event-dependent physical observer frames. -/
noncomputable def SourceOrderFrameCompatibilityPacket.ofFaithfulPlacementStandardFrame
    {Register : Type u} {Value : Type v} {Event : Type w}
    [DecidableEq Register] [Fintype Event] [DecidableEq Event]
    [DecidableEq Value]
    {L : SemanticEventLog Register Value Event}
    (placement : FaithfulRankSpatialCausalPlacement L) :
    SourceOrderFrameCompatibilityPacket Register Value Event Unit :=
  SourceOrderFrameCompatibilityPacket.ofFaithfulPlacement placement
    (fun _ ↦ OrientedLorentzEquiv.refl)

namespace SourceOrderFrameCompatibilityPacket

variable {Register : Type u} {Value : Type v} {Event : Type w}
  {Chart : Type z}
variable [DecidableEq Register] [Fintype Event]
variable (P : SourceOrderFrameCompatibilityPacket Register Value Event Chart)

/-- The event-local internal frame selected by the packet's supplied Lorentz
transport.  This is supplied frame data, not a result of the event order. -/
def eventFrame (e : Event) : FrameHyperboloid :=
  (P.standardToEventFrame e).mapFrame standardFrame

/-- The supplied base-chart event frames extended to all charts by the
already committed overlap cocycle.  Its underlying event-germ atlas is
definitionally the atlas carried by `sourceCharts`. -/
def frameSoldering : EventFrameSoldering Event Chart :=
  EventFrameSoldering.ofBaseFrame P.sourceCharts.atlas
    P.sourceCharts.base P.eventFrame

/-- At the declared base chart, the soldering recovers the packet's supplied
per-event frame. -/
theorem frameSoldering_frame_base (e : Event) :
    P.frameSoldering.frame P.sourceCharts.base e = P.eventFrame e := by
  apply Subtype.ext
  exact P.sourceCharts.atlas.overlap.lorentz_self
    P.sourceCharts.base (P.eventFrame e : Herm2)

/-- The source-generated reflexive informational order agrees with the
declared future-cone displacement relation in every supplied chart.  The
atlas and cone equivalence are fields of `sourceCharts`, not consequences of
the order. -/
theorem generatedBeforeEq_iff_futureCausal
    (i : Chart) (e f : Event) :
    P.sourceCharts.semanticLog.GeneratedBeforeEq e f ↔
      IsFutureCausal
        (P.sourceCharts.atlas.displacement i e f) :=
  P.sourceCharts.generatedBeforeEq_iff_displacement_futureCausal i e f

/-- The strict source-generated order agrees exactly with the non-diagonal
part of the supplied future-cone displacement relation in every chart. -/
theorem generatedBefore_iff_futureCausal_and_ne
    (i : Chart) (e f : Event) :
    P.sourceCharts.semanticLog.GeneratedBefore e f ↔
      IsFutureCausal (P.sourceCharts.atlas.displacement i e f) ∧ e ≠ f :=
  P.sourceCharts.generatedBefore_iff_displacement_futureCausal_and_ne i e f

/-- The source Gram quotient has real rank three. -/
theorem sourceFrame_finrank :
    Module.finrank ℝ FrameQuotient = 3 :=
  frameQuotient_finrank

/-- Every supplied event frame has a rank-three Lorentz-orthogonal rest
fiber. -/
theorem eventRest_finrank (e : Event) :
    Module.finrank ℝ (RestSpace (P.eventFrame e)) = 3 :=
  finrank_restSpace (P.eventFrame e)

/-- The canonical quotient-to-standard-rest bridge transported into the
supplied rest fiber at one event.  Its existence uses the explicitly
supplied `standardToEventFrame`; provenance does not select this map. -/
def frameQuotientEquivEventRest (e : Event) :
    FrameQuotient ≃ₗ[ℝ] RestSpace (P.eventFrame e) :=
  frameQuotientEquivStandardRest.trans
    ((P.standardToEventFrame e).restEquiv standardFrame)

/-- The event-local bridge preserves exactly the source quotient Gram and
the positive rest metric.  This is an algebraic metric bridge, not a physical
rod, scale, volume, or spacetime-metric identification. -/
theorem frameQuotientEquivEventRest_preserves_metric
    (e : Event) (q r : FrameQuotient) :
    restMetric (P.eventFrame e) (P.frameQuotientEquivEventRest e q)
        (P.frameQuotientEquivEventRest e r) = quotientGram q r := by
  change restMetric
      ((P.standardToEventFrame e).mapFrame standardFrame)
      ((P.standardToEventFrame e).restEquiv standardFrame
        (frameQuotientEquivStandardRest q))
      ((P.standardToEventFrame e).restEquiv standardFrame
        (frameQuotientEquivStandardRest r)) = quotientGram q r
  rw [(P.standardToEventFrame e).restEquiv_preserves_metric]
  exact frameQuotientEquivStandardRest_preserves_metric q r

/-- The exact finite consequences exposed by this composition: source order
versus the supplied future cone, rank three on both spatial carriers, the
algebraic source-unit-direction/celestial/null-ray equivalences, and the
event-local metric bridge.  It also exposes the existing algebraic soldering
package, the independently constructed source-carrier rank four and
`(+---)` axes, rank-three rest fibers, Lorentz overlap cocycles, and
conditional celestial covariance of future-null displacements.  The
nonclaims in the module header are not hidden in this conjunction. -/
def FiniteConsequences : Prop :=
  (∀ (i : Chart) (e f : Event),
      P.sourceCharts.semanticLog.GeneratedBeforeEq e f ↔
        IsFutureCausal
          (P.sourceCharts.atlas.displacement i e f)) ∧
  (∀ (i : Chart) (e f : Event),
      P.sourceCharts.semanticLog.GeneratedBefore e f ↔
        IsFutureCausal (P.sourceCharts.atlas.displacement i e f) ∧ e ≠ f) ∧
  AlgebraicSolderingConsequences P.frameSoldering ∧
  Module.finrank ℝ SourceSpacetimeCarrier = 4 ∧
  sourceLorentzQ sourceTimeAxis = 1 ∧
  (∀ i : Fin 3, sourceLorentzQ (sourceSpatialAxis i) = -1) ∧
  (∀ i : Fin 3,
    sourceLorentzB sourceTimeAxis (sourceSpatialAxis i) = 0) ∧
  (∀ i j : Fin 3,
    sourceLorentzB (sourceSpatialAxis i) (sourceSpatialAxis j) =
      if i = j then -1 else 0) ∧
  (∀ q : SourceUnitDirection,
    IsSourceFutureCausal (sourceUnitNullVector q)) ∧
  Module.finrank ℝ FrameQuotient = 3 ∧
  Nonempty (SourceUnitDirection ≃ CelestialSphere) ∧
  Nonempty (SourceUnitDirection ≃ FutureNullRay) ∧
  (∀ e : Event, Module.finrank ℝ (RestSpace (P.eventFrame e)) = 3) ∧
  (∀ (e : Event) (q r : FrameQuotient),
      restMetric (P.eventFrame e) (P.frameQuotientEquivEventRest e q)
          (P.frameQuotientEquivEventRest e r) = quotientGram q r)

/-- Every packet supplies the exact finite consequences listed in
`FiniteConsequences`.  The conjunction is an exposed interface, not an
exhaustiveness theorem about every proposition derivable from the packet. -/
theorem finiteConsequences : P.FiniteConsequences := by
  refine ⟨P.generatedBeforeEq_iff_futureCausal,
    P.generatedBefore_iff_futureCausal_and_ne,
    P.frameSoldering.algebraicConsequences,
    sourceSpacetimeCarrier_finrank, sourceLorentzQ_timeAxis,
    sourceLorentzQ_spatialAxis, sourceLorentzB_time_spatialAxis,
    sourceLorentzB_spatialAxis,
    sourceUnitNullVector_futureCausal,
    sourceFrame_finrank,
    ⟨sourceUnitDirectionEquivCelestial⟩,
    ⟨sourceUnitDirectionEquivFutureNullRay⟩,
    P.eventRest_finrank, P.frameQuotientEquivEventRest_preserves_metric⟩

#print axioms sourceUnitDirectionEquivCelestial
#print axioms sourceUnitDirectionEquivFutureNullRay
#print axioms SourceOrderFrameCompatibilityPacket.ofFaithfulPlacement
#print axioms SourceOrderFrameCompatibilityPacket.ofFaithfulPlacementStandardFrame
#print axioms SourceOrderFrameCompatibilityPacket.frameSoldering_frame_base
#print axioms SourceOrderFrameCompatibilityPacket.generatedBeforeEq_iff_futureCausal
#print axioms SourceOrderFrameCompatibilityPacket.generatedBefore_iff_futureCausal_and_ne
#print axioms SourceOrderFrameCompatibilityPacket.sourceFrame_finrank
#print axioms SourceOrderFrameCompatibilityPacket.eventRest_finrank
#print axioms SourceOrderFrameCompatibilityPacket.frameQuotientEquivEventRest_preserves_metric
#print axioms SourceOrderFrameCompatibilityPacket.finiteConsequences

end SourceOrderFrameCompatibilityPacket

/-! ## Fully constructed non-chain control packet -/

namespace SourceDerivedSpacetime.BooleanDiamondPlacement

/-- The explicit authenticated Boolean diamond, its order-faithful source-carrier
placement, one global chart, and the standard-frame gauge form a complete
finite source-order/frame packet with no additional existence field.  This
is a non-chain control, not a manifold or physical-frame construction. -/
noncomputable def sourceOrderFramePacket :
    SourceOrderFrameCompatibilityPacket (Fin 4) Bool (Fin 4) Unit :=
  SourceOrderFrameCompatibilityPacket.ofFaithfulPlacementStandardFrame
    faithfulPlacement

/-- All advertised finite carrier, cone, direction, and soldering
consequences hold on the explicit non-chain control. -/
theorem sourceOrderFramePacket_consequences :
    sourceOrderFramePacket.FiniteConsequences :=
  sourceOrderFramePacket.finiteConsequences

#print axioms sourceOrderFramePacket
#print axioms sourceOrderFramePacket_consequences

end SourceDerivedSpacetime.BooleanDiamondPlacement

end

end OPH
