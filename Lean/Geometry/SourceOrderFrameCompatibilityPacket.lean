import Geometry.EventPopulationChartInterface
import Geometry.SpatialReadbackSoldering

/-!
# Finite source-order/frame compatibility packet

This module composes three already separate finite ingredients:

* an authenticated semantic log whose reflexive generated order is attached,
  by a supplied event-population/chart interface, to the declared future cone;
* the independently proved rank-three source Gram quotient; and
* one supplied oriented Lorentz transport per event, used to choose that
  event's internal observer frame and to transport the canonical
  `FrameQuotient`--`RestSpace standardFrame` bridge into its rest fiber.

The composition is deliberately noncircular.  The semantic log determines
the informational order, but it does not construct the atlas, cone
identification, population, or frame transports.  Conversely, the supplied
frames do not add provenance edges.  The rank-three source quotient is a
spatial readback carrier and is not identified with the four-dimensional
event carrier.  Its unit Gram directions are, however, exactly equivalent
to the already constructed celestial two-sphere and hence to future-null
rays.  This is an algebraic unit-direction theorem, not a physical sky or
signal identification.

**Nonclaims.**  This packet does not identify its algebraic unit sphere with
a physical celestial screen and proves no source derivation of a
Lorentz/conformal group, no physical-causality or faithful embedding theorem,
no volume or count-density law, no manifoldlikeness or dimension estimator,
no refinement/cofinal convergence theorem, and no topological,
smooth-manifold, or continuum-spacetime limit.  Every geometric attachment
used below remains visible as supplied data.
-/

namespace OPH

noncomputable section

open C1Lorentz C2Soldering CausalComposition EventPopulation
  PrimitivePortFrameQuotient Provenance

universe u v w z

/-! ## Exact source-direction/celestial bridge -/

/-- Unit directions of the exact source Gram quotient.  This is the unit
sphere of the proved rank-three carrier, not yet a physical celestial
screen. -/
abbrev SourceUnitDirection :=
  {q : FrameQuotient // quotientGram q q = 1}

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
  sourceCharts : SourceDerivedOrderEventPopulationChartInterface
    Register Value Event Chart
  standardToEventFrame : Event → OrientedLorentzEquiv

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
  EventFrameSoldering.ofBaseFrame P.sourceCharts.geometry.atlas
    P.sourceCharts.geometry.base P.eventFrame

/-- At the declared base chart, the soldering recovers the packet's supplied
per-event frame. -/
theorem frameSoldering_frame_base (e : Event) :
    P.frameSoldering.frame P.sourceCharts.geometry.base e = P.eventFrame e := by
  apply Subtype.ext
  exact P.sourceCharts.geometry.atlas.overlap.lorentz_self
    P.sourceCharts.geometry.base (P.eventFrame e : Herm2)

/-- The source-generated reflexive informational order agrees with the
declared future-cone displacement relation in every supplied chart.  The
atlas and cone equivalence are fields of `sourceCharts`, not consequences of
the order. -/
theorem generatedBeforeEq_iff_futureCausal
    (i : Chart) (e f : Event) :
    P.sourceCharts.semanticLog.GeneratedBeforeEq e f ↔
      IsFutureCausal
        (P.sourceCharts.geometry.atlas.displacement i e f) :=
  P.sourceCharts.generatedBeforeEq_iff_displacement_futureCausal i e f

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
package: ambient rank four, rank-three rest fibers, Lorentz overlap cocycles,
and conditional celestial covariance of future-null displacements.  The
nonclaims in the module header are not hidden in this conjunction. -/
def FiniteConsequences : Prop :=
  (∀ (i : Chart) (e f : Event),
      P.sourceCharts.semanticLog.GeneratedBeforeEq e f ↔
        IsFutureCausal
          (P.sourceCharts.geometry.atlas.displacement i e f)) ∧
  AlgebraicSolderingConsequences P.frameSoldering ∧
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
    P.frameSoldering.algebraicConsequences, sourceFrame_finrank,
    ⟨sourceUnitDirectionEquivCelestial⟩,
    ⟨sourceUnitDirectionEquivFutureNullRay⟩,
    P.eventRest_finrank, P.frameQuotientEquivEventRest_preserves_metric⟩

#print axioms sourceUnitDirectionEquivCelestial
#print axioms sourceUnitDirectionEquivFutureNullRay
#print axioms SourceOrderFrameCompatibilityPacket.frameSoldering_frame_base
#print axioms SourceOrderFrameCompatibilityPacket.generatedBeforeEq_iff_futureCausal
#print axioms SourceOrderFrameCompatibilityPacket.sourceFrame_finrank
#print axioms SourceOrderFrameCompatibilityPacket.eventRest_finrank
#print axioms SourceOrderFrameCompatibilityPacket.frameQuotientEquivEventRest_preserves_metric
#print axioms SourceOrderFrameCompatibilityPacket.finiteConsequences

end SourceOrderFrameCompatibilityPacket

end

end OPH
