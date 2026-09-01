import Geometry.CelestialSoldering

/-!
# C2: observer-frame soldering and exact premise reduction

This module adds an observer-frame readback to an event-germ atlas and proves
its overlap covariance.  It also states the C2 premise reduction precisely.
The algebraic soldering packet automatically supplies affine/Lorentz cocycle,
displacement, interval, frame, celestial, and ambient module dimension/signature
facts.  It does **not** discharge the independent source-atlas, population,
separation, open-chart, physical-cone, refinement, causal, or clock receipts.

The result is intentionally an interface theorem, not a claim that a physical
event manifold has been constructed.

`Geometry.SourceDerivedSpacetimeCarrier` constructs an algebraic `1+3`
carrier and, from a certified source-height/spatial placement plus its explicit
faithfulness conditions, a finite one-chart instance. That result does
not turn the open-atlas, physical-causality, refinement, volume, clock, or
smooth-manifold receipts below into theorems.
-/

namespace OPH.C2Soldering

noncomputable section

open OPH.C1Lorentz

/-- An event-germ atlas together with a future-unit observer-frame readback
that transforms through the same overlap Lorentz map. -/
structure EventFrameSoldering (Event Chart : Type*)
    extends EventGermAtlas Event Chart where
  frame : Chart → Event → FrameHyperboloid
  frame_overlap : ∀ {i j e}, visible i e → visible j e →
    frame j e = (overlap.lorentz i j).mapFrame (frame i e)

/-- Construct all chart-frame readbacks from one base-chart frame selection.
The overlap law follows from the single Lorentz cocycle; it is not a second
independent covariance assumption. -/
def EventFrameSoldering.ofBaseFrame {Event Chart : Type*}
    (A : EventGermAtlas Event Chart) (base : Chart)
    (seedFrame : Event → FrameHyperboloid) :
    EventFrameSoldering Event Chart where
  toEventGermAtlas := A
  frame := fun i e => (A.overlap.lorentz base i).mapFrame (seedFrame e)
  frame_overlap := by
    intro i j e _ _
    apply Subtype.ext
    exact A.overlap.lorentz_cocycle base i j (seedFrame e).1

/-- The exact C2 algebraic consequences.  Celestial covariance is included
as a theorem schema only when the displacement is future null. -/
def AlgebraicSolderingConsequences {Event Chart : Type*}
    (S : EventFrameSoldering Event Chart) : Prop :=
  Module.finrank ℝ Herm2 = 4 ∧
  (∀ u : FrameHyperboloid, Module.finrank ℝ (RestSpace u) = 3) ∧
  (∀ i j k x, S.overlap.act i k x =
    S.overlap.act j k (S.overlap.act i j x)) ∧
  (∀ i j x, S.overlap.act j i (S.overlap.act i j x) = x) ∧
  (∀ (i j : Chart) (p q : Event),
    S.visible i p → S.visible i q → S.visible j p → S.visible j q →
      S.displacement j p q =
        S.overlap.lorentz i j (S.displacement i p q)) ∧
  (∀ (i j : Chart) (e : Event),
    S.visible i e → S.visible j e →
      S.frame j e =
        (S.overlap.lorentz i j).mapFrame (S.frame i e)) ∧
  (∀ (i j : Chart) (p q : Event),
    S.visible i p → S.visible i q → S.visible j p → S.visible j q →
      lorentzQ (S.displacement j p q) = lorentzQ (S.displacement i p q)) ∧
  (∀ (i j : Chart) (p q : Event)
      (_hip : S.visible i p) (_hiq : S.visible i q)
      (_hjp : S.visible j p) (_hjq : S.visible j q)
      (hi : IsFutureNull (S.displacement i p q))
      (hj : IsFutureNull (S.displacement j p q)),
    S.celestialSolder j p q hj =
      (S.overlap.lorentz i j).celestialAction
        (S.celestialSolder i p q hi))

theorem EventFrameSoldering.algebraicConsequences
    {Event Chart : Type*} (S : EventFrameSoldering Event Chart) :
    AlgebraicSolderingConsequences S := by
  refine ⟨finrank_Herm2, finrank_restSpace,
    S.overlap.act_cocycle, S.overlap.act_reverse_left, ?_, ?_, ?_, ?_⟩
  · intro i j p q hip hiq hjp hjq
    exact S.toEventGermAtlas.displacement_overlap hip hiq hjp hjq
  · intro i j e hie hje
    exact S.frame_overlap hie hje
  · intro i j p q hip hiq hjp hjq
    exact S.toEventGermAtlas.interval_overlap hip hiq hjp hjq
  · intro i j p q hip hiq hjp hjq hi hj
    exact S.toEventGermAtlas.celestialSolder_overlap hip hiq hjp hjq hi hj

/-! ## A finite nondegenerate algebraic control

This two-event/two-chart model verifies that the contract is consistent and
that its affine translation sector is nontrivial.  It is a control model, not
a source realization: its coordinates and chart offsets are declared below.
-/

def controlEventPoint (e : Fin 2) : Herm2 :=
  if e = 0 then 0 else (1, ![1, 0, 0])

def controlChartOffset (i : Fin 2) : Herm2 :=
  if i = 0 then 0 else (2, ![0, 1, 0])

def controlOverlap : LorentzOverlapCocycle (Fin 2) where
  lorentz := fun _ _ => OrientedLorentzEquiv.refl
  translation := fun i j => controlChartOffset j - controlChartOffset i
  lorentz_self := by simp [OrientedLorentzEquiv.refl]
  lorentz_cocycle := by simp [OrientedLorentzEquiv.refl]
  translation_self := by simp
  translation_cocycle := by
    intro i j k
    simp [OrientedLorentzEquiv.refl]

/-- Two declared event points in two affinely shifted charts with one standard
internal frame. -/
def controlSoldering : EventFrameSoldering (Fin 2) (Fin 2) where
  visible := fun _ _ => True
  coordinate := fun i e => controlEventPoint e + controlChartOffset i
  overlap := controlOverlap
  coordinate_overlap := by
    intro i j e _ _
    simp [controlOverlap, LorentzOverlapCocycle.act,
      OrientedLorentzEquiv.refl]
  frame := fun _ _ => standardFrame
  frame_overlap := by
    intro i j e _ _
    rfl

theorem control_chart_translation_nonzero :
    controlOverlap.translation 0 1 ≠ 0 := by
  intro h
  have ht := congrArg Prod.fst h
  norm_num [controlOverlap, controlChartOffset] at ht

theorem control_displacement_nonzero :
    controlSoldering.displacement 0 0 1 ≠ 0 := by
  intro h
  have ht := congrArg Prod.fst h
  norm_num [EventGermAtlas.displacement, controlSoldering,
    controlEventPoint, controlChartOffset] at ht

theorem control_displacement_futureNull :
    IsFutureNull (controlSoldering.displacement 0 0 1) := by
  norm_num [EventGermAtlas.displacement, controlSoldering,
    controlEventPoint, controlChartOffset, IsFutureNull, lorentzQ,
    spatialNormSq, Fin.sum_univ_succ]

#print axioms EventFrameSoldering.algebraicConsequences
#print axioms control_chart_translation_nonzero
#print axioms control_displacement_futureNull

end

end OPH.C2Soldering
