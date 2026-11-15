import Geometry.GoldenSourceCountLimit

set_option autoImplicit false

/-!
# Null boundaries for the supplied flat causal relation

The metric is Euclidean in the three spatial coordinates. The measure is
the declared product coordinate volume, restricted to the same finite box
as the actual golden event counts. Fubini reduces each null cone to a
singleton in time; null-boundary hypotheses are conclusions here.
-/

namespace OPH.SourceCausalBoundary

open MeasureTheory Set Filter OPH.GoldenSourceCountLimit
open scoped Topology

noncomputable section

abbrev Space := Fin 3 → ℝ
abbrev Spacetime := ℝ × Space

def spatialNorm (x : Space) : ℝ := ‖(WithLp.toLp 2 x : EuclideanSpace ℝ (Fin 3))‖

def margin (c : ℝ) (p q : Spacetime) : ℝ :=
  c * (q.1 - p.1) - spatialNorm (q.2 - p.2)

def future (c : ℝ) (p : Spacetime) : Set Spacetime := {q | 0 ≤ margin c p q}
def past (c : ℝ) (p : Spacetime) : Set Spacetime := {q | 0 ≤ margin c q p}
def diamond (c : ℝ) (p q : Spacetime) : Set Spacetime := future c p ∩ past c q

theorem continuous_spatialNorm : Continuous spatialNorm :=
  (PiLp.continuous_toLp 2 (fun _ : Fin 3 => ℝ)).norm

theorem continuous_margin (c : ℝ) : Continuous (fun z : Spacetime × Spacetime =>
    margin c z.1 z.2) := by
  unfold margin
  exact (continuous_const.mul (continuous_snd.fst.sub continuous_fst.fst)).sub
    (continuous_spatialNorm.comp (continuous_snd.snd.sub continuous_fst.snd))

theorem continuous_margin_right (c : ℝ) (p : Spacetime) : Continuous (margin c p) :=
  (continuous_margin c).comp (continuous_const.prodMk continuous_id)

theorem continuous_margin_left (c : ℝ) (p : Spacetime) :
    Continuous (fun q => margin c q p) :=
  (continuous_margin c).comp (continuous_id.prodMk continuous_const)

theorem diamond_measurable (c : ℝ) (p q : Spacetime) : MeasurableSet (diamond c p q) :=
  (isClosed_le continuous_const (continuous_margin_right c p)).measurableSet.inter
    (isClosed_le continuous_const (continuous_margin_left c q)).measurableSet

/-- A measurable time graph has zero product mass even inside the clipped window. -/
theorem time_graph_null (T L : ℝ) (g : Space → ℝ) (hg : Measurable g) :
    spaceTimeVolume T L {z | z.1 = g z.2} = 0 := by
  unfold spaceTimeVolume
  have hm : MeasurableSet {z : Spacetime | z.1 = g z.2} :=
    measurableSet_eq_fun measurable_fst (hg.comp measurable_snd)
  rw [Measure.prod_apply_symm hm]
  have hs (x : Space) :
      (fun t : ℝ => (t, x)) ⁻¹' {z : Spacetime | z.1 = g z.2} = {g x} := by
    ext t
    simp
  simp only [hs, measure_singleton, lintegral_zero]

theorem future_null (T L : ℝ) {c : ℝ} (hc : c ≠ 0) (p : Spacetime) :
    spaceTimeVolume T L {q | margin c p q = 0} = 0 := by
  have hg : Continuous (fun x : Space => p.1 + spatialNorm (x - p.2) / c) :=
    continuous_const.add ((continuous_spatialNorm.comp
      (continuous_id.sub continuous_const)).div_const c)
  have heq : {q | margin c p q = 0} =
      {q : Spacetime | q.1 = p.1 + spatialNorm (q.2 - p.2) / c} := by
    ext q
    simp only [mem_setOf_eq, margin]
    constructor
    · intro h
      have hh : q.1 - p.1 = spatialNorm (q.2 - p.2) / c :=
        (eq_div_iff hc).mpr (by nlinarith)
      linarith
    · intro h
      have hh := (div_mul_cancel₀ (spatialNorm (q.2 - p.2)) hc)
      rw [h]
      nlinarith
  rw [heq]
  exact time_graph_null T L _ hg.measurable

theorem past_null (T L : ℝ) {c : ℝ} (hc : c ≠ 0) (p : Spacetime) :
    spaceTimeVolume T L {q | margin c q p = 0} = 0 := by
  have hg : Continuous (fun x : Space => p.1 - spatialNorm (p.2 - x) / c) :=
    continuous_const.sub ((continuous_spatialNorm.comp
      (continuous_const.sub continuous_id)).div_const c)
  have heq : {q | margin c q p = 0} =
      {q : Spacetime | q.1 = p.1 - spatialNorm (p.2 - q.2) / c} := by
    ext q
    simp only [mem_setOf_eq, margin]
    constructor
    · intro h
      have hh : p.1 - q.1 = spatialNorm (p.2 - q.2) / c :=
        (eq_div_iff hc).mpr (by nlinarith)
      linarith
    · intro h
      have hh := (div_mul_cancel₀ (spatialNorm (p.2 - q.2)) hc)
      rw [h]
      nlinarith
  rw [heq]
  exact time_graph_null T L _ hg.measurable

theorem future_frontier_null (T L : ℝ) {c : ℝ} (hc : c ≠ 0) (p : Spacetime) :
    spaceTimeVolume T L (frontier (future c p)) = 0 := by
  apply measure_mono_null _ (future_null T L hc p)
  intro q hq
  exact (frontier_le_subset_eq continuous_const (continuous_margin_right c p) hq).symm

theorem past_frontier_null (T L : ℝ) {c : ℝ} (hc : c ≠ 0) (p : Spacetime) :
    spaceTimeVolume T L (frontier (past c p)) = 0 := by
  apply measure_mono_null _ (past_null T L hc p)
  intro q hq
  exact (frontier_le_subset_eq continuous_const (continuous_margin_left c p) hq).symm

theorem diamond_frontier_null (T L : ℝ) {c : ℝ} (hc : c ≠ 0) (p q : Spacetime) :
    spaceTimeVolume T L (frontier (diamond c p q)) = 0 := by
  apply measure_mono_null (frontier_inter_subset _ _)
  exact measure_union_null
    (measure_mono_null inter_subset_left (future_frontier_null T L hc p))
    (measure_mono_null inter_subset_right (past_frontier_null T L hc q))

/-- Pairwise null separation is null for the actual product reference measure. -/
theorem pair_null (T L : ℝ) {c : ℝ} (hc : c ≠ 0) :
    ((spaceTimeVolume T L).prod (spaceTimeVolume T L))
      {z : Spacetime × Spacetime | margin c z.1 z.2 = 0} = 0 := by
  apply Measure.measure_prod_null_of_ae_null
    (isClosed_eq (continuous_margin c) continuous_const).measurableSet
  exact Eventually.of_forall (fun p => future_null T L hc p)

/-- The coincidence diagonal is contained in the null pair relation. -/
theorem diagonal_null (T L : ℝ) {c : ℝ} (hc : c ≠ 0) :
    ((spaceTimeVolume T L).prod (spaceTimeVolume T L))
      {z : Spacetime × Spacetime | z.1 = z.2} = 0 := by
  apply measure_mono_null _ (pair_null T L hc)
  intro z hz
  simp only [mem_setOf_eq] at hz ⊢
  simp [hz, margin, spatialNorm]

/-- Fixed causal diamonds satisfy the fixed-set count theorem with no supplied
null-frontier hypothesis. This statement does not identify generated intervals. -/
theorem source_diamond_count_tendsto {T L c : ℝ} (hT : 0 ≤ T) (hL : 0 < L)
    (hc : 0 < c) (p q : Spacetime) :
    Tendsto (fun n =>
      (sourceDelta (n + 2) L c * (L ^ 3 / (Nat.fib (n + 2) : ℝ) ^ 3)) *
        (eventCount (n + 2) T L (sourceDelta (n + 2) L c) (diamond c p q) : ℝ))
      atTop (𝓝 ((spaceTimeVolume T L).real (diamond c p q))) :=
  source_event_count_tendsto hT hL hc _ (diamond_measurable c p q)
    (diamond_frontier_null T L hc.ne' p q)

#print axioms diamond_frontier_null
#print axioms pair_null
#print axioms diagonal_null
#print axioms source_diamond_count_tendsto

end
end OPH.SourceCausalBoundary
