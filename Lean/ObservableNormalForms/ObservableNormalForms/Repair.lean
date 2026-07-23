import Mathlib

/-!
# Strong local repair and projection

This module formalizes the collar-section criterion.  The audit exposed a missing
edge-case hypothesis in the prose statement: the write space must be
nonempty.  Without it the repair domain `X_W × X_D` can be empty, so a repair
exists vacuously even when the collar projection is not surjective.  We prove
both the corrected theorem (in a form stronger than the manuscript because
no finiteness is needed classically) and the empty-write-space counterexample.
-/

namespace ObservableNormalForms

universe u v

section StrongRepair

variable {W : Type u} {D : Type v}

/-- A total exact repair into `R` that preserves the protected collar and
fixes states already in `R`. -/
structure StrongRepair (R : Set (W × D)) where
  toFun : W × D → {z : W × D // z ∈ R}
  preservesCollar : ∀ z : W × D, (toFun z).1.2 = z.2
  fixes : ∀ (z : W × D) (_hz : z ∈ R), (toFun z).1 = z

/-- Surjectivity of the collar projection restricted to the consistency
relation. -/
def CollarProjectionSurjective (R : Set (W × D)) : Prop :=
  Function.Surjective (fun z : {z : W × D // z ∈ R} => z.1.2)

theorem collar_projection_surjective_of_strongRepair
    [Nonempty W] {R : Set (W × D)} (repair : StrongRepair R) :
    CollarProjectionSurjective R := by
  intro d
  let w : W := Classical.choice (inferInstance : Nonempty W)
  refine ⟨repair.toFun (w, d), ?_⟩
  exact repair.preservesCollar (w, d)

noncomputable def strongRepairOfCollarProjection
    {R : Set (W × D)} (hsurj : CollarProjectionSurjective R) :
    StrongRepair R := by
  classical
  exact {
    toFun := fun z =>
      if hz : z ∈ R then
        ⟨z, hz⟩
      else
        Classical.choose (hsurj z.2)
    preservesCollar := by
      intro z
      by_cases hz : z ∈ R
      · simp [hz]
      · simp only [hz, ↓reduceDIte]
        exact Classical.choose_spec (hsurj z.2)
    fixes := by
      intro z hz
      simp [hz]
  }

/-- Collar-section criterion.  Finiteness is unnecessary for this
existence equivalence in classical Lean; nonemptiness of `W` is necessary. -/
theorem strongRepair_exists_iff_collarProjection_surjective
    [Nonempty W] (R : Set (W × D)) :
    Nonempty (StrongRepair R) ↔ CollarProjectionSurjective R := by
  constructor
  · rintro ⟨repair⟩
    exact collar_projection_surjective_of_strongRepair repair
  · intro hsurj
    exact ⟨strongRepairOfCollarProjection hsurj⟩

/-- No-repair corollary, with the necessary nonempty-write hypothesis:
one missing collar value is a certificate that no total collar-preserving
exact repair exists. -/
theorem no_strongRepair_of_missing_collar
    [Nonempty W] {R : Set (W × D)} {d : D}
    (hmissing : ¬ ∃ z : {z : W × D // z ∈ R}, z.1.2 = d) :
    IsEmpty (StrongRepair R) := by
  constructor
  intro repair
  apply hmissing
  exact (collar_projection_surjective_of_strongRepair repair) d

/-- Machine-checked edge-case counterexample to the manuscript's original
statement without `Nonempty X_W`: the domain is empty, hence a strong repair
exists vacuously, but the projection from the empty relation to `Unit` is not
surjective. -/
theorem empty_write_space_counterexample :
    Nonempty (StrongRepair (∅ : Set (Empty × Unit))) ∧
      ¬ CollarProjectionSurjective (∅ : Set (Empty × Unit)) := by
  constructor
  · let repair : StrongRepair (∅ : Set (Empty × Unit)) := {
      toFun := fun z => nomatch z.1
      preservesCollar := by
        intro z
        exact Empty.elim z.1
      fixes := by
        intro z
        exact Empty.elim z.1
    }
    exact ⟨repair⟩
  · intro hsurj
    rcases hsurj () with ⟨z, _⟩
    exact Empty.elim z.1.1

end StrongRepair

section RobustMargin

variable {X : Type u} {D : Type v}
variable [MetricSpace D]

/-- Distance from a protected collar value to the projected consistency
relation.  For a real-valued `infDist`, nonemptiness of the relation is a
load-bearing hypothesis in the results below. -/
noncomputable def repairMargin (R : Set X) (collar : X → D) (d : D) : ℝ :=
  Metric.infDist d (collar '' R)

/-- Robust no-repair margin.  If state-space distance
dominates collar distance and `R` is nonempty, distance to consistency is at
least the projected repair margin. -/
theorem robust_no_repair_margin
    [PseudoMetricSpace X] {R : Set X} {collar : X → D}
    (hR : R.Nonempty)
    (hdominates : ∀ x y : X,
      dist (collar x) (collar y) ≤ dist x y)
    (x : X) :
    repairMargin R collar (collar x) ≤ Metric.infDist x R := by
  rw [repairMargin, Metric.le_infDist hR]
  intro y hy
  have hyImage : collar y ∈ collar '' R := ⟨y, hy, rfl⟩
  exact (Metric.infDist_le_dist_of_mem hyImage).trans
    (hdominates x y)

/-- Compact nonempty projected relations give a strictly positive margin at
every collar value outside the projection. -/
theorem repairMargin_pos_of_compact
    [TopologicalSpace X] {R : Set X} {collar : X → D} {d : D}
    (hcompact : IsCompact R)
    (hR : R.Nonempty)
    (hcollar : Continuous collar)
    (hmissing : d ∉ collar '' R) :
    0 < repairMargin R collar d := by
  exact ((hcompact.image hcollar).isClosed.notMem_iff_infDist_pos
    (hR.image collar)).mp hmissing

/-- Why `R.Nonempty` is necessary with ordinary real-valued metric distance:
Mathlib, like many treatments, defines distance to the empty set as zero. -/
theorem empty_relation_repairMargin_zero
    (collar : X → D) (d : D) :
    repairMargin (∅ : Set X) collar d = 0 := by
  simp [repairMargin, Metric.infDist_empty]

end RobustMargin

end ObservableNormalForms
