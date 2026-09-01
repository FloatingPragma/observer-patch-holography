import Geometry.CausalOrderComposition
import Geometry.EventGermDisplacement
import ObserverPatchHolography.Provenance.SemanticEventProvenance

/-!
# Source-derived finite causal-chart interface

Authenticated semantic commits are the finite event carrier and
`SemanticEventLog.GeneratedBeforeEq` is their reflexive source-generated
order. This module attaches that order to a supplied finite Lorentz-chart
package without introducing a second population map or precedence field.

The exact interface carries:

* the authenticated semantic log;
* a supplied finite `EventGermAtlas`;
* total visibility at the finite cutoff;
* base-chart separation; and
* exact agreement between generated order and the base-chart future cone.

The overlap cocycle propagates separation, order/cone agreement, interval
invariance, and parent-edge future causality to every supplied chart. The
atlas, total visibility, separation, and cone identification remain explicit
hypotheses. Nothing here constructs open chart images, topology, physical
signal causality, an operational clock, count--volume calibration,
manifoldlikeness, refinement, curvature, or a continuum limit.
-/

namespace OPH.EventPopulation

noncomputable section

open OPH.C1Lorentz OPH.C2Soldering OPH.CausalComposition OPH.Provenance

/-- The identity affine-Lorentz overlap cocycle on any chart-label type.
This small active helper is used by one-chart source-derived placements; it
does not introduce a chart, event population, or causal order. -/
def trivialOverlap (Chart : Type*) : LorentzOverlapCocycle Chart where
  lorentz := fun _ _ => OrientedLorentzEquiv.refl
  translation := fun _ _ => 0
  lorentz_self := by simp [OrientedLorentzEquiv.refl]
  lorentz_cocycle := by simp [OrientedLorentzEquiv.refl]
  translation_self := fun _ => rfl
  translation_cocycle := by
    intro i j k
    simp [OrientedLorentzEquiv.refl]

/-- A source-native finite causal-chart interface. The authenticated semantic
log indexes the complete finite event carrier directly, so there is no
separate population map and no freely chosen precedence relation. -/
structure SourceDerivedCausalChartInterface
    (Register Value Event Chart : Type*) [DecidableEq Register]
    [Fintype Event] where
  /-- The authenticated semantic event log on the event carrier. -/
  semanticLog : SemanticEventLog Register Value Event
  /-- The supplied finite Lorentz chart package. -/
  atlas : EventGermAtlas Event Chart
  /-- This finite precursor uses total chart visibility. Local/open atlas
  domains remain a continuum obligation. -/
  total : ∀ (i : Chart) (e : Event), atlas.visible i e
  /-- A base chart used to state the one load-bearing cone attachment. -/
  base : Chart
  /-- Distinct events have distinct coordinates in the base chart. -/
  separation : Function.Injective (atlas.coordinate base)
  /-- The source-generated reflexive order, not an auxiliary precedence
  field, agrees exactly with the supplied base-chart future-cone order. -/
  generated_cone_base : ∀ e f : Event,
    semanticLog.GeneratedBeforeEq e f ↔
      causalLE (atlas.coordinate base e) (atlas.coordinate base f)

namespace SourceDerivedCausalChartInterface

variable {Register Value Event Chart : Type*} [DecidableEq Register]
variable [Fintype Event]
variable (S : SourceDerivedCausalChartInterface Register Value Event Chart)

/-- Every chart coordinate is the overlap transport of the base-chart
coordinate. -/
theorem coordinate_from_base (i : Chart) (e : Event) :
    S.atlas.coordinate i e =
      S.atlas.overlap.act S.base i (S.atlas.coordinate S.base e) :=
  S.atlas.coordinate_overlap (S.total S.base e) (S.total i e)

/-- Source-event separation propagates to every supplied chart. -/
theorem separation_all (i : Chart) :
    Function.Injective (S.atlas.coordinate i) := by
  intro e f h
  apply S.separation
  have h2 :
      S.atlas.overlap.act S.base i (S.atlas.coordinate S.base e) =
        S.atlas.overlap.act S.base i (S.atlas.coordinate S.base f) := by
    rw [← S.coordinate_from_base i e, ← S.coordinate_from_base i f, h]
  have h3 := congrArg (S.atlas.overlap.act i S.base) h2
  rwa [S.atlas.overlap.act_reverse_left, S.atlas.overlap.act_reverse_left]
    at h3

/-- The exact source-generated order/cone correspondence propagates from
the base chart to every supplied chart. -/
theorem generatedBeforeEq_iff_causalLE (i : Chart) (e f : Event) :
    S.semanticLog.GeneratedBeforeEq e f ↔
      causalLE (S.atlas.coordinate i e) (S.atlas.coordinate i f) := by
  rw [S.coordinate_from_base i e, S.coordinate_from_base i f,
    causalLE_act_iff]
  exact S.generated_cone_base e f

/-- Equivalent displacement form of the source-order/cone
correspondence. -/
theorem generatedBeforeEq_iff_displacement_futureCausal
    (i : Chart) (e f : Event) :
    S.semanticLog.GeneratedBeforeEq e f ↔
      IsFutureCausal (S.atlas.displacement i e f) :=
  S.generatedBeforeEq_iff_causalLE i e f

/-- The strict authenticated precedence is exactly the irreflexive part of
the supplied cone order in every chart. This is a finite faithful-cone
statement conditional on the supplied chart and cone attachment, not a
construction of that attachment. -/
theorem generatedBefore_iff_causalLE_and_ne (i : Chart) (e f : Event) :
    S.semanticLog.GeneratedBefore e f ↔
      causalLE (S.atlas.coordinate i e) (S.atlas.coordinate i f) ∧ e ≠ f := by
  constructor
  · intro h
    exact ⟨(S.generatedBeforeEq_iff_causalLE i e f).mp (Or.inr h),
      fun hef ↦ S.semanticLog.generatedBefore_irrefl e (hef ▸ h)⟩
  · rintro ⟨hcone, hne⟩
    rcases (S.generatedBeforeEq_iff_causalLE i e f).mpr hcone with heq | hlt
    · exact absurd heq hne
    · exact hlt

/-- The strict authenticated precedence is likewise the non-diagonal part
of the future-causal displacement relation. -/
theorem generatedBefore_iff_displacement_futureCausal_and_ne
    (i : Chart) (e f : Event) :
    S.semanticLog.GeneratedBefore e f ↔
      IsFutureCausal (S.atlas.displacement i e f) ∧ e ≠ f :=
  S.generatedBefore_iff_causalLE_and_ne i e f

/-- Lorentz intervals of source-event displacements are chart invariant. -/
theorem interval_chart_invariant (i j : Chart) (e f : Event) :
    lorentzQ (S.atlas.displacement j e f) =
      lorentzQ (S.atlas.displacement i e f) :=
  S.atlas.interval_overlap (S.total i e) (S.total i f)
    (S.total j e) (S.total j f)

/-- Every authenticated direct parent is future-causal in every supplied
chart. -/
theorem parentEdge_displacement_futureCausal
    (i : Chart) {e f : Event} (h : S.semanticLog.ParentEdge e f) :
    IsFutureCausal (S.atlas.displacement i e f) :=
  (S.generatedBeforeEq_iff_displacement_futureCausal i e f).mp
    (Or.inr (S.semanticLog.parentEdge_generatedBefore h))

end SourceDerivedCausalChartInterface

-- Axiom audit: no project axiom or admission is used.
#print axioms OPH.EventPopulation.SourceDerivedCausalChartInterface.separation_all
#print axioms OPH.EventPopulation.SourceDerivedCausalChartInterface.generatedBeforeEq_iff_causalLE
#print axioms OPH.EventPopulation.SourceDerivedCausalChartInterface.generatedBefore_iff_causalLE_and_ne
#print axioms OPH.EventPopulation.SourceDerivedCausalChartInterface.parentEdge_displacement_futureCausal

end

end OPH.EventPopulation
