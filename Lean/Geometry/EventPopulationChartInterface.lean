import Geometry.CausalOrderComposition
import Geometry.EventGermDisplacement

/-!
# OL-A4: the typed event-population chart interface at the finite level

The event-manifold row of the observation ledger (OL-A4) owes, beyond the
committed soldering and causal-order receipts, the attachments that would
turn charted records into an event manifold: an event population, a chart
assignment, separation of distinct events, compatibility of the causal order
with the declared future cone, and the overlap law between charts.  This
file types those five obligations as one antecedent bundle,
`EventPopulationChartInterface`, and proves the composed consequence: every
inhabitant carries the event-manifold-shaped package at the finite level.

The bundle names, as declared fields: a finite type of events populated by
an injection from the committed rich-fibre region carrier
(`QFT/RichFibreRegionalNet.lean`); a chart assignment into the declared
Hermitian Lorentz module, packaged as a committed `EventGermAtlas`
(`Geometry/EventGermDisplacement.lean`) with total visibility, which carries
the chart coordinates, the `LorentzOverlapCocycle` transition law of
`Geometry/LorentzOverlapCocycle.lean`, and the overlap compatibility clause;
a separation clause (distinct events have distinct base-chart coordinates);
a cone clause (declared causal precedence holds exactly when the base-chart
displacement lies in the declared future cone of
`Geometry/CausalOrderComposition.lean`); and a population order clause (the
committed region order lands exactly in the declared precedence through the
population map).

The composed theorem `eventManifoldPackage` derives, from the bundle alone
and the committed receipts: the chart assignment is injective in every
chart, not only the base chart; Lorentz intervals of displacements are
chart-invariant through the committed cocycle receipt; the declared
precedence is realized in every chart as the cone order of chart
displacements; the precedence is reflexive, transitive, and antisymmetric,
inherited from the committed cone partial order; the committed region order
embeds into the cone order of populated chart images in every chart; the
exact triple-overlap cocycle holds; and the coordinate module has real
dimension four.  The chart propagation reuses the committed theorems
`causalLE_act_iff` and `richRegion_le_iff_causalLE` rather than reproving
them.

Nonvacuity: `witnessInterface` inhabits the bundle on the committed
six-element rich-fibre carrier with two affinely shifted charts, and the
cone clause is exercised on a nontrivial causal pair (bottom below window
r86, timelike separated) and on an incomparable spacelike pair (windows r86
and r88).

Load-bearing receipts: dropping the separation clause admits an inhabitant
whose chart identifies two distinct events
(`separation_clause_load_bearing`), and dropping the cone clause admits an
order-violating chart assignment, a time-reversed chart under which a
declared-precedent pair leaves the future cone
(`cone_clause_load_bearing`).

**Declared data.**  Every field of the bundle is a declared hypothesis.
The event type, the population map, the precedence relation, the charts,
and the overlap transitions are supplied, not derived.  The witness
inhabitant is a declared finite construction on committed data.

**Falsifier.**  The composition fails if some inhabitant of the bundle has
a non-injective chart in some chart label, a chart-dependent Lorentz
interval, a chart where the precedence and the cone order disagree, a
precedence violating a partial-order law, or a populated pair breaking the
committed order equivalence; it also fails if the witness violates a
clause, or if either weakened bundle refuses the exhibited counterexample.

**Nonclaims.**  No source atlas: nothing here produces charts, events, or
transitions from observer repair data.  No physical event population: the
events are a declared finite type, not records of a physical world.  No
clock, no continuum, no open-chart topology, and no four-dimensional
manifold claim beyond the finite chart package landing in the
four-dimensional declared module.  The stable-causality and open-image
conditions (PR-16) and the observer-to-physical-spacetime and causal
attachment premise (PR-52) stay open.  For the observation-ledger row
OL-A4 this interface supports at most a formal precursor partial
classification; the source-atlas, physical-event-population, and clock
attachments that row owes are not represented by any field here and stay
owed.
-/

namespace OPH.EventPopulation

noncomputable section

open OPH.C1Lorentz OPH.C2Soldering OPH.CausalComposition OPH.QFT

/-! ## The typed antecedent bundle -/

/-- The typed OL-A4 antecedent bundle: a finite event population injected
from the committed rich-fibre region carrier, a total chart assignment into
the declared Hermitian Lorentz module through a committed event-germ atlas,
a separation clause, a base-chart cone clause against the declared future
cone, and the population order clause tying the committed region order to
the declared precedence.  Every field is a declared hypothesis. -/
structure EventPopulationChartInterface (Event Chart : Type*)
    [Fintype Event] where
  /-- The event population: each committed rich-fibre region names one
  event. -/
  populate : RichRegion → Event
  /-- The population map is injective: distinct committed regions name
  distinct events. -/
  populate_injective : Function.Injective populate
  /-- The declared causal precedence on the event population. -/
  prec : Event → Event → Prop
  /-- The chart assignment into the declared Lorentz module, carried by a
  committed event-germ atlas: chart coordinates, the affine Lorentz
  overlap cocycle, and the overlap compatibility clause. -/
  atlas : EventGermAtlas Event Chart
  /-- Total visibility: the finite chart package has no visibility gaps,
  every chart sees every event. -/
  total : ∀ (i : Chart) (e : Event), atlas.visible i e
  /-- A declared base chart. -/
  base : Chart
  /-- The separation clause: distinct events have distinct base-chart
  coordinates. -/
  separation : Function.Injective (atlas.coordinate base)
  /-- The cone clause: declared precedence holds exactly when the
  base-chart displacement lies in the declared future cone. -/
  cone_base : ∀ e f : Event,
    prec e f ↔ causalLE (atlas.coordinate base e) (atlas.coordinate base f)
  /-- The population order clause: the committed region order lands
  exactly in the declared precedence through the population map. -/
  populate_le_iff : ∀ U V : RichRegion,
    RichRegion.le U V ↔ prec (populate U) (populate V)

namespace EventPopulationChartInterface

variable {Event Chart : Type*} [Fintype Event]
variable (S : EventPopulationChartInterface Event Chart)

/-- Every chart coordinate is the overlap transition of the base-chart
coordinate: the committed overlap compatibility clause at total
visibility. -/
theorem coordinate_from_base (i : Chart) (e : Event) :
    S.atlas.coordinate i e =
      S.atlas.overlap.act S.base i (S.atlas.coordinate S.base e) :=
  S.atlas.coordinate_overlap (S.total S.base e) (S.total i e)

/-- Separation propagates from the base chart to every chart: overlap
transitions are invertible by the committed reverse-transition receipt. -/
theorem separation_all (i : Chart) :
    Function.Injective (S.atlas.coordinate i) := by
  intro e f h
  apply S.separation
  have h2 : S.atlas.overlap.act S.base i (S.atlas.coordinate S.base e) =
      S.atlas.overlap.act S.base i (S.atlas.coordinate S.base f) := by
    rw [← S.coordinate_from_base i e, ← S.coordinate_from_base i f, h]
  have h3 := congrArg (S.atlas.overlap.act i S.base) h2
  rwa [S.atlas.overlap.act_reverse_left, S.atlas.overlap.act_reverse_left]
    at h3

/-- The cone clause propagates from the base chart to every chart through
the committed invariance of the cone order under affine overlap-cocycle
transitions (`causalLE_act_iff`), reused rather than reproved. -/
theorem cone_all (i : Chart) (e f : Event) :
    S.prec e f ↔
      causalLE (S.atlas.coordinate i e) (S.atlas.coordinate i f) := by
  rw [S.coordinate_from_base i e, S.coordinate_from_base i f,
    causalLE_act_iff]
  exact S.cone_base e f

/-- The declared precedence is realized in every chart as the cone order
of chart displacements. -/
theorem prec_iff_displacement_futureCausal (i : Chart) (e f : Event) :
    S.prec e f ↔ IsFutureCausal (S.atlas.displacement i e f) :=
  S.cone_all i e f

/-- Lorentz intervals of displacements are chart-invariant: the committed
interval-overlap receipt at total visibility. -/
theorem interval_chart_invariant (i j : Chart) (e f : Event) :
    lorentzQ (S.atlas.displacement j e f) =
      lorentzQ (S.atlas.displacement i e f) :=
  S.atlas.interval_overlap (S.total i e) (S.total i f)
    (S.total j e) (S.total j f)

/-- Reflexivity of the declared precedence, inherited from the committed
cone order. -/
theorem prec_refl (e : Event) : S.prec e e :=
  (S.cone_base e e).mpr (causalLE_refl _)

/-- Transitivity of the declared precedence, inherited from the committed
cone order. -/
theorem prec_trans {e f g : Event}
    (h1 : S.prec e f) (h2 : S.prec f g) : S.prec e g :=
  (S.cone_base e g).mpr
    (causalLE_trans ((S.cone_base e f).mp h1) ((S.cone_base f g).mp h2))

/-- Antisymmetry of the declared precedence, inherited from the committed
cone antisymmetry through the separation clause. -/
theorem prec_antisymm {e f : Event}
    (h1 : S.prec e f) (h2 : S.prec f e) : e = f :=
  S.separation
    (causalLE_antisymm ((S.cone_base e f).mp h1) ((S.cone_base f e).mp h2))

/-- The committed region order embeds into the cone order of populated
chart images, in every chart, through the population order clause and the
committed chart invariance. -/
theorem populate_le_cone (i : Chart) (U V : RichRegion) :
    RichRegion.le U V ↔
      causalLE (S.atlas.coordinate i (S.populate U))
        (S.atlas.coordinate i (S.populate V)) :=
  (S.populate_le_iff U V).trans (S.cone_all i _ _)

end EventPopulationChartInterface

/-! ## The composed theorem -/

/-- **The composed OL-A4 interface theorem.**  Every inhabitant of the
typed bundle carries the event-manifold-shaped package at the finite
level: the coordinate module has real dimension four; the chart atlas is
injective in every chart; the population map is injective; Lorentz
intervals of displacements are chart-invariant; the declared precedence is
realized in every chart as the cone order of chart displacements; the
precedence is reflexive, transitive, and antisymmetric, inherited from the
committed cone partial order; the committed region order embeds into the
cone order of populated chart images in every chart; and the exact
triple-overlap cocycle holds.  Every clause beyond the declared fields is
derived by reusing committed receipts. -/
theorem eventManifoldPackage {Event Chart : Type*} [Fintype Event]
    (S : EventPopulationChartInterface Event Chart) :
    Module.finrank ℝ Herm2 = 4 ∧
    (∀ i : Chart, Function.Injective (S.atlas.coordinate i)) ∧
    Function.Injective S.populate ∧
    (∀ (i j : Chart) (e f : Event),
      lorentzQ (S.atlas.displacement j e f) =
        lorentzQ (S.atlas.displacement i e f)) ∧
    (∀ (i : Chart) (e f : Event),
      S.prec e f ↔ IsFutureCausal (S.atlas.displacement i e f)) ∧
    (∀ e : Event, S.prec e e) ∧
    (∀ e f g : Event, S.prec e f → S.prec f g → S.prec e g) ∧
    (∀ e f : Event, S.prec e f → S.prec f e → e = f) ∧
    (∀ (i : Chart) (U V : RichRegion),
      RichRegion.le U V ↔
        causalLE (S.atlas.coordinate i (S.populate U))
          (S.atlas.coordinate i (S.populate V))) ∧
    (∀ (i j k : Chart) (x : Herm2),
      S.atlas.overlap.act i k x =
        S.atlas.overlap.act j k (S.atlas.overlap.act i j x)) :=
  ⟨finrank_Herm2, S.separation_all, S.populate_injective,
    S.interval_chart_invariant, S.prec_iff_displacement_futureCausal,
    S.prec_refl, fun _ _ _ h1 h2 => S.prec_trans h1 h2,
    fun _ _ h1 h2 => S.prec_antisymm h1 h2, S.populate_le_cone,
    S.atlas.overlap.act_cocycle⟩

/-! ## The finite inhabitant -/

/-- Chart offsets of the witness: the base chart is unshifted, the second
chart is shifted by one unit along the scalar axis. -/
def witnessOffset (i : Fin 2) : Herm2 :=
  if i = 0 then 0 else (1, 0)

/-- The witness overlap cocycle: identity Lorentz part, translation by the
offset difference. -/
def witnessOverlap : LorentzOverlapCocycle (Fin 2) where
  lorentz := fun _ _ => OrientedLorentzEquiv.refl
  translation := fun i j => witnessOffset j - witnessOffset i
  lorentz_self := by simp [OrientedLorentzEquiv.refl]
  lorentz_cocycle := by simp [OrientedLorentzEquiv.refl]
  translation_self := by simp
  translation_cocycle := by
    intro i j k
    simp [OrientedLorentzEquiv.refl]

/-- The witness atlas: the committed declared vertex assignment
`richVertex` of the causal-order composition, shifted per chart. -/
def witnessAtlas : EventGermAtlas RichRegion (Fin 2) where
  visible := fun _ _ => True
  coordinate := fun i U => richVertex U + witnessOffset i
  overlap := witnessOverlap
  coordinate_overlap := by
    intro i j e _ _
    simp [witnessOverlap, LorentzOverlapCocycle.act,
      OrientedLorentzEquiv.refl]

/-- **Nonvacuity.**  The typed bundle is inhabited on the committed
six-element rich-fibre carrier with two affinely shifted charts: identity
population, the committed region order as precedence, and the committed
declared vertex assignment as base chart.  The cone clause is the
committed order equivalence `richRegion_le_iff_causalLE`, reused rather
than reproved. -/
def witnessInterface : EventPopulationChartInterface RichRegion (Fin 2) where
  populate := id
  populate_injective := fun _ _ h => h
  prec := RichRegion.le
  atlas := witnessAtlas
  total := fun _ _ => trivial
  base := 0
  separation := by
    intro e f h
    apply richVertex_injective
    simpa [witnessAtlas, witnessOffset] using h
  cone_base := by
    intro e f
    simpa [witnessAtlas, witnessOffset] using richRegion_le_iff_causalLE e f
  populate_le_iff := fun _ _ => Iff.rfl

/-- The witness base chart is the committed declared vertex assignment. -/
theorem witness_coordinate_base (U : RichRegion) :
    witnessInterface.atlas.coordinate 0 U = richVertex U := by
  simp [witnessInterface, witnessAtlas, witnessOffset]

/-- The two witness charts differ: the overlap clause of the witness is
exercised on a nontrivial affine transition. -/
theorem witness_charts_differ :
    witnessInterface.atlas.coordinate 1 RichRegion.bot ≠
      witnessInterface.atlas.coordinate 0 RichRegion.bot := by
  intro h
  have ht := congrArg Prod.fst h
  norm_num [witnessInterface, witnessAtlas, witnessOffset, richVertex] at ht

/-- **The cone clause is genuinely exercised.**  The witness carries a
nontrivial causal pair: bottom precedes window r86, the base-chart
displacement is future-causal and timelike; and an incomparable pair: the
windows r86 and r88 are unrelated in both directions and spacelike
separated. -/
theorem witness_nontrivial_pairs :
    (witnessInterface.prec RichRegion.bot RichRegion.r86 ∧
      causalLE (witnessInterface.atlas.coordinate 0 RichRegion.bot)
        (witnessInterface.atlas.coordinate 0 RichRegion.r86) ∧
      IsTimelikeSep (witnessInterface.atlas.coordinate 0 RichRegion.bot)
        (witnessInterface.atlas.coordinate 0 RichRegion.r86)) ∧
    (¬ witnessInterface.prec RichRegion.r86 RichRegion.r88 ∧
      ¬ witnessInterface.prec RichRegion.r88 RichRegion.r86 ∧
      IsSpacelikeSep (witnessInterface.atlas.coordinate 0 RichRegion.r86)
        (witnessInterface.atlas.coordinate 0 RichRegion.r88)) := by
  have hb := witness_coordinate_base RichRegion.bot
  have h86 := witness_coordinate_base RichRegion.r86
  have h88 := witness_coordinate_base RichRegion.r88
  refine ⟨⟨(by decide : RichRegion.le RichRegion.bot RichRegion.r86), ?_, ?_⟩,
    (by decide : ¬ RichRegion.le RichRegion.r86 RichRegion.r88),
    (by decide : ¬ RichRegion.le RichRegion.r88 RichRegion.r86), ?_⟩
  · rw [hb, h86]
    exact richVertex_nondegenerate.1.2
  · rw [hb, h86]
    exact richVertex_nondegenerate.1.1
  · rw [h86, h88]
    exact richVertex_nondegenerate.2.2.1.1

/-! ## Load-bearing receipts -/

/-- The bundle with the separation clause removed; every other clause is
kept verbatim.  Used only to exhibit the collapse counterexample. -/
structure BundleWithoutSeparation (Event Chart : Type*) [Fintype Event] where
  populate : RichRegion → Event
  populate_injective : Function.Injective populate
  prec : Event → Event → Prop
  atlas : EventGermAtlas Event Chart
  total : ∀ (i : Chart) (e : Event), atlas.visible i e
  base : Chart
  cone_base : ∀ e f : Event,
    prec e f ↔ causalLE (atlas.coordinate base e) (atlas.coordinate base f)
  populate_le_iff : ∀ U V : RichRegion,
    RichRegion.le U V ↔ prec (populate U) (populate V)

/-- The trivial one-label overlap cocycle. -/
def trivialOverlap (Chart : Type*) : LorentzOverlapCocycle Chart where
  lorentz := fun _ _ => OrientedLorentzEquiv.refl
  translation := fun _ _ => 0
  lorentz_self := by simp [OrientedLorentzEquiv.refl]
  lorentz_cocycle := by simp [OrientedLorentzEquiv.refl]
  translation_self := fun _ => rfl
  translation_cocycle := by
    intro i j k
    simp [OrientedLorentzEquiv.refl]

/-- The collapsed atlas: a seventh event shares the chart point of the
committed bottom region. -/
def collapsedAtlas : EventGermAtlas (Option RichRegion) Unit where
  visible := fun _ _ => True
  coordinate := fun _ e => richVertex (e.getD RichRegion.bot)
  overlap := trivialOverlap Unit
  coordinate_overlap := by
    intro i j e _ _
    simp [trivialOverlap, LorentzOverlapCocycle.act,
      OrientedLorentzEquiv.refl]

/-- An inhabitant of the separation-free bundle whose chart identifies two
distinct events. -/
def separationDropWitness : BundleWithoutSeparation (Option RichRegion) Unit where
  populate := some
  populate_injective := fun _ _ h => Option.some.inj h
  prec := fun e f =>
    causalLE (richVertex (e.getD RichRegion.bot))
      (richVertex (f.getD RichRegion.bot))
  atlas := collapsedAtlas
  total := fun _ _ => trivial
  base := ()
  cone_base := fun _ _ => Iff.rfl
  populate_le_iff := fun U V => by
    simpa using richRegion_le_iff_causalLE U V

/-- **The separation clause is load-bearing.**  The separation-free bundle
admits an inhabitant in which two distinct events share one chart point,
so the population is not charted faithfully.  The full bundle excludes
this by the separation clause (`EventPopulationChartInterface.separation_all`). -/
theorem separation_clause_load_bearing :
    ∃ e f : Option RichRegion, e ≠ f ∧
      separationDropWitness.atlas.coordinate separationDropWitness.base e =
        separationDropWitness.atlas.coordinate separationDropWitness.base f :=
  ⟨none, some RichRegion.bot, by simp, rfl⟩

/-- The bundle with the cone clause removed; every other clause is kept
verbatim.  Used only to exhibit the order-violating counterexample. -/
structure BundleWithoutCone (Event Chart : Type*) [Fintype Event] where
  populate : RichRegion → Event
  populate_injective : Function.Injective populate
  prec : Event → Event → Prop
  atlas : EventGermAtlas Event Chart
  total : ∀ (i : Chart) (e : Event), atlas.visible i e
  base : Chart
  separation : Function.Injective (atlas.coordinate base)
  populate_le_iff : ∀ U V : RichRegion,
    RichRegion.le U V ↔ prec (populate U) (populate V)

/-- The time-reversed atlas: the negated committed vertex assignment. -/
def timeReversedAtlas : EventGermAtlas RichRegion Unit where
  visible := fun _ _ => True
  coordinate := fun _ U => -richVertex U
  overlap := trivialOverlap Unit
  coordinate_overlap := by
    intro i j e _ _
    simp [trivialOverlap, LorentzOverlapCocycle.act,
      OrientedLorentzEquiv.refl]

/-- An inhabitant of the cone-free bundle with an order-violating chart:
population, separation, and the population order clause all hold, and the
chart is the time-reversed committed assignment. -/
def coneDropWitness : BundleWithoutCone RichRegion Unit where
  populate := id
  populate_injective := fun _ _ h => h
  prec := RichRegion.le
  atlas := timeReversedAtlas
  total := fun _ _ => trivial
  base := ()
  separation := by
    intro e f h
    exact richVertex_injective (neg_injective h)
  populate_le_iff := fun _ _ => Iff.rfl

/-- **The cone clause is load-bearing.**  The cone-free bundle admits an
inhabitant whose chart violates the declared order: bottom precedes window
r86 in the declared precedence, and the time-reversed chart displacement
leaves the declared future cone.  The full bundle excludes this in every
chart (`EventPopulationChartInterface.cone_all`). -/
theorem cone_clause_load_bearing :
    ∃ e f : RichRegion,
      coneDropWitness.prec e f ∧
        ¬ causalLE (coneDropWitness.atlas.coordinate coneDropWitness.base e)
          (coneDropWitness.atlas.coordinate coneDropWitness.base f) := by
  refine ⟨RichRegion.bot, RichRegion.r86,
    (by decide : RichRegion.le RichRegion.bot RichRegion.r86), ?_⟩
  intro h
  have ht := h.2
  simp only [coneDropWitness, timeReversedAtlas, richVertex,
    Prod.fst_sub] at ht
  norm_num at ht

/-! ## The composed receipt -/

/-- The composed OL-A4 interface receipt: every inhabitant of the typed
bundle carries the finite event-manifold-shaped package; the bundle is
inhabited on the committed rich-fibre carrier with the cone clause
exercised on a timelike pair and a spacelike pair over two genuinely
distinct charts; dropping the separation clause admits a
two-events-one-chart collapse; and dropping the cone clause admits an
order-violating chart.  No clause attaches a source atlas, a physical
event population, a clock, a topology, or a physical spacetime. -/
theorem eventPopulationChartInterface_receipt :
    (∀ (Event Chart : Type) [Fintype Event]
        (S : EventPopulationChartInterface Event Chart),
      Module.finrank ℝ Herm2 = 4 ∧
      (∀ i : Chart, Function.Injective (S.atlas.coordinate i)) ∧
      Function.Injective S.populate ∧
      (∀ (i j : Chart) (e f : Event),
        lorentzQ (S.atlas.displacement j e f) =
          lorentzQ (S.atlas.displacement i e f)) ∧
      (∀ (i : Chart) (e f : Event),
        S.prec e f ↔ IsFutureCausal (S.atlas.displacement i e f)) ∧
      (∀ e : Event, S.prec e e) ∧
      (∀ e f g : Event, S.prec e f → S.prec f g → S.prec e g) ∧
      (∀ e f : Event, S.prec e f → S.prec f e → e = f) ∧
      (∀ (i : Chart) (U V : RichRegion),
        RichRegion.le U V ↔
          causalLE (S.atlas.coordinate i (S.populate U))
            (S.atlas.coordinate i (S.populate V))) ∧
      (∀ (i j k : Chart) (x : Herm2),
        S.atlas.overlap.act i k x =
          S.atlas.overlap.act j k (S.atlas.overlap.act i j x))) ∧
    ((witnessInterface.prec RichRegion.bot RichRegion.r86 ∧
        causalLE (witnessInterface.atlas.coordinate 0 RichRegion.bot)
          (witnessInterface.atlas.coordinate 0 RichRegion.r86) ∧
        IsTimelikeSep (witnessInterface.atlas.coordinate 0 RichRegion.bot)
          (witnessInterface.atlas.coordinate 0 RichRegion.r86)) ∧
      (¬ witnessInterface.prec RichRegion.r86 RichRegion.r88 ∧
        ¬ witnessInterface.prec RichRegion.r88 RichRegion.r86 ∧
        IsSpacelikeSep (witnessInterface.atlas.coordinate 0 RichRegion.r86)
          (witnessInterface.atlas.coordinate 0 RichRegion.r88)) ∧
      witnessInterface.atlas.coordinate 1 RichRegion.bot ≠
        witnessInterface.atlas.coordinate 0 RichRegion.bot) ∧
    (∃ e f : Option RichRegion, e ≠ f ∧
      separationDropWitness.atlas.coordinate separationDropWitness.base e =
        separationDropWitness.atlas.coordinate separationDropWitness.base f) ∧
    (∃ e f : RichRegion,
      coneDropWitness.prec e f ∧
        ¬ causalLE (coneDropWitness.atlas.coordinate coneDropWitness.base e)
          (coneDropWitness.atlas.coordinate coneDropWitness.base f)) :=
  ⟨fun _ _ _ S => eventManifoldPackage S,
    ⟨witness_nontrivial_pairs.1, witness_nontrivial_pairs.2,
      witness_charts_differ⟩,
    separation_clause_load_bearing, cone_clause_load_bearing⟩

-- Axiom audit: no project axiom or admission is used.
#print axioms OPH.EventPopulation.EventPopulationChartInterface.separation_all
#print axioms OPH.EventPopulation.EventPopulationChartInterface.cone_all
#print axioms OPH.EventPopulation.EventPopulationChartInterface.prec_iff_displacement_futureCausal
#print axioms OPH.EventPopulation.EventPopulationChartInterface.interval_chart_invariant
#print axioms OPH.EventPopulation.EventPopulationChartInterface.prec_antisymm
#print axioms OPH.EventPopulation.EventPopulationChartInterface.populate_le_cone
#print axioms OPH.EventPopulation.eventManifoldPackage
#print axioms OPH.EventPopulation.witness_coordinate_base
#print axioms OPH.EventPopulation.witness_charts_differ
#print axioms OPH.EventPopulation.witness_nontrivial_pairs
#print axioms OPH.EventPopulation.separation_clause_load_bearing
#print axioms OPH.EventPopulation.cone_clause_load_bearing
#print axioms OPH.EventPopulation.eventPopulationChartInterface_receipt

end

end OPH.EventPopulation
