import Geometry.SourceDerivedSpacetimeCarrier
import Mathlib.Tactic

namespace OPH.RecordedConeScale

open OPH.Provenance OPH.PrimitivePortFrameQuotient OPH.SourceDerivedSpacetime

variable {R V E : Type} [DecidableEq R] [DecidableEq V] [Fintype E] [DecidableEq E]

/-- Turn the tested closed lower bounds and strict upper bounds on a squared
clock scale into an actual faithful placement of the authenticated log. -/
noncomputable def placementOfBounds (L : SemanticEventLog R V E) (x : E → FrameQuotient)
    (k : ℝ) (hk : 0 < k)
    (separate : ∀ {e f}, L.sourceHeight e = L.sourceHeight f → x e = x f → e = f)
    (lower : ∀ {e f}, L.ParentEdge e f →
      quotientGram (x f-x e) (x f-x e) / ((L.sourceHeight f : ℝ)-L.sourceHeight e)^2 ≤ k)
    (upper : ∀ {e f}, L.sourceHeight e < L.sourceHeight f → ¬ L.GeneratedBefore e f →
      k < quotientGram (x f-x e) (x f-x e) / ((L.sourceHeight f : ℝ)-L.sourceHeight e)^2) :
    FaithfulRankSpatialCausalPlacement L := by
  let q : RankSpatialCausalPlacement L := {
    timeScale := Real.sqrt k
    timeScale_pos := Real.sqrt_pos.2 hk
    spatial := x
    edge_speed_sq := by
      intro e f h
      have hpos : 0 < (L.sourceHeight f : ℝ) - L.sourceHeight e := by
        apply sub_pos.mpr; exact_mod_cast L.sourceHeight_lt_of_parent h
      simpa [mul_pow, Real.sq_sqrt hk.le] using
        (div_le_iff₀ (sq_pos_of_pos hpos)).mp (lower h) }
  apply FaithfulRankSpatialCausalPlacement.ofIncomparableSpacelike q separate
  intro e f h hn
  have hpos : 0 < (L.sourceHeight f : ℝ) - L.sourceHeight e := by
    apply sub_pos.mpr; exact_mod_cast h
  simpa [q, mul_pow, Real.sq_sqrt hk.le] using
    (lt_div_iff₀ (sq_pos_of_pos hpos)).mp (upper h hn)

/-- Strict separation of any increasing-height nonancestor is necessary
for the same registered clock, since even a null pair is in the cone. -/
theorem necessaryUpper {L : SemanticEventLog R V E} (q : FaithfulRankSpatialCausalPlacement L)
    {e f : E} (hh : L.sourceHeight e < L.sourceHeight f) (hn : ¬ L.GeneratedBefore e f) :
    q.timeScale^2 < quotientGram (q.spatial f - q.spatial e) (q.spatial f - q.spatial e) /
      ((L.sourceHeight f : ℝ) - L.sourceHeight e)^2 := by
  have hpos : 0 < (L.sourceHeight f : ℝ) - L.sourceHeight e := by
    apply sub_pos.mpr; exact_mod_cast hh
  apply (lt_div_iff₀ (sq_pos_of_pos hpos)).2
  by_contra! h
  have hc : sourceCausalLE (q.eventPoint e) (q.eventPoint f) := by
    constructor
    · change 0 ≤ (q.timeScale * (L.sourceHeight f : ℝ) - q.timeScale * L.sourceHeight e) *
          (q.timeScale * (L.sourceHeight f : ℝ) - q.timeScale * L.sourceHeight e) -
          quotientGram (q.spatial f - q.spatial e) (q.spatial f - q.spatial e)
      nlinarith
    · change 0 ≤ q.timeScale * (L.sourceHeight f : ℝ) - q.timeScale * L.sourceHeight e
      rw [← mul_sub]; exact mul_nonneg q.timeScale_pos.le hpos.le
  rcases q.cone_supported e f hc with hEq | before
  · subst f; exact (lt_irrefl _ hh)
  · exact hn before

theorem necessaryLower {L : SemanticEventLog R V E} (q : RankSpatialCausalPlacement L)
    {e f : E} (hp : L.ParentEdge e f) :
    quotientGram (q.spatial f - q.spatial e) (q.spatial f - q.spatial e) /
      ((L.sourceHeight f : ℝ)-L.sourceHeight e)^2 ≤ q.timeScale^2 := by
  have hpos : 0 < (L.sourceHeight f : ℝ) - L.sourceHeight e := by
    apply sub_pos.mpr; exact_mod_cast L.sourceHeight_lt_of_parent hp
  apply (div_le_iff₀ (sq_pos_of_pos hpos)).2
  simpa only [mul_pow] using q.edge_speed_sq hp

abbrev Slot := Fin 4
abbrev Port := Fin 12
abbrev Word := Slot → Port × Port

def linked : Word := ![(0,1),(3,4),(0,10),(2,10)]
def separate : Word := ![(0,1),(3,6),(6,10),(2,10)]
def touches (w : Word) (i : Slot) (p : Port) : Bool := p == (w i).1 || p == (w i).2

/-- A later event reads a parent's port before any intervening event stamps it. -/
def parent (w : Word) (i j : Slot) : Bool :=
  decide (i < j) && (List.finRange 12).any fun p =>
    touches w i p && touches w j p && (List.finRange 4).all fun t =>
      !(decide (i < t ∧ t < j) && touches w t p)

def height (w : Word) : ℕ → Slot → ℕ
  | 0, _ => 0
  | fuel+1, j => Finset.univ.sup fun i => if parent w i j then height w fuel i + 1 else 0

def tag (w : Word) (i : Slot) := (w i, i, height w 4 i)

theorem equal_local_tags : tag linked 0 = tag separate 0 ∧ tag linked 3 = tag separate 3 ∧
    height linked 4 0 = 0 ∧ height linked 4 3 = 2 := by decide

theorem transmitted : Relation.TransGen (fun i j => parent linked i j = true) 0 3 :=
  .tail (.single (b := (2 : Slot)) (by decide)) (by decide)

theorem isolated : ¬ Relation.TransGen (fun i j => parent separate i j = true) 0 3 := by
  have noEdge : ∀ j, parent separate 0 j = false := by decide
  have none : ∀ j, ¬ Relation.TransGen (fun i k => parent separate i k = true) 0 j := by
    intro j h
    induction h with
    | single hp => simp [noEdge] at hp
    | tail _ _ ih => exact ih
  exact none 3

/-- No dimension, nonlinear placement, or choice of causal classifier can
make the two different recorded orders a common function of these local tags. -/
theorem matched_comparison {Point : Type} (place : (Port × Port) × Slot × ℕ → Point)
    (rel : Point → Point → Prop) :
    ¬ ((rel (place (tag linked 0)) (place (tag linked 3)) ↔
      Relation.TransGen (fun i j => parent linked i j = true) 0 3) ∧
      (rel (place (tag separate 0)) (place (tag separate 3)) ↔
      Relation.TransGen (fun i j => parent separate i j = true) 0 3)) := by
  rintro ⟨h, k⟩
  have ht := h.mpr transmitted
  rw [equal_local_tags.1, equal_local_tags.2.1] at ht
  exact isolated (k.mp ht)

end OPH.RecordedConeScale
