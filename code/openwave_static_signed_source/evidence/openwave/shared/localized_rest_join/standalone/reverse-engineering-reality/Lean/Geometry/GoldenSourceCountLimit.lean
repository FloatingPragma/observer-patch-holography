import Geometry.GoldenSourceAssignment
import Mathlib.MeasureTheory.Integral.DominatedConvergence
import Mathlib.MeasureTheory.Measure.Prod

set_option autoImplicit false

/-!
# Actual golden continuity-set counts

The finite counts use the exact golden coordinates and their actual equal
Lebesgue masses. Indicator convergence is derived from the proved assignment
bounds and dominated convergence; no weak-measure or count-convergence input
is supplied. Representatives need not lie in their assigned cells.

These results concern a declared mathematical population and coordinate
measure. They do not select a source population, causal read law or clock,
and do not identify generated-order interval or strict-pair counts.
-/

open MeasureTheory Set Filter
open scoped Topology BigOperators

namespace OPH.GoldenSourceCountLimit

noncomputable section

variable {X : Type*} [PseudoMetricSpace X] [MeasurableSpace X]
  {ι : Type*} [Fintype ι]

noncomputable def sampled (C : ι → Set X) (s : ι → X) (A : Set X) (x : X) : ℝ :=
  ∑ i, (C i).indicator (fun _ => A.indicator (fun _ => (1 : ℝ)) (s i)) x

omit [PseudoMetricSpace X] [MeasurableSpace X] in
theorem sampled_eq (C : ι → Set X) (s : ι → X) (A : Set X)
    (hd : Pairwise (fun i j => Disjoint (C i) (C j)))
    {x : X} {i : ι} (hi : x ∈ C i) :
    sampled C s A x = A.indicator (fun _ => (1 : ℝ)) (s i) := by
  classical
  unfold sampled
  rw [Finset.sum_eq_single i]
  · exact indicator_of_mem hi _
  · intro j hj hji
    apply indicator_of_notMem
    intro hx
    exact Set.disjoint_left.mp (hd hji) hx hi
  · simp

/-- Discontinuous counts converge when cell displacement vanishes and the
target boundary is null. The statement assumes no integral or count limit. -/
theorem partition_count_tendsto {ι : ℕ → Type*} [∀ n, Fintype (ι n)]
    (μ : Measure X) [IsFiniteMeasure μ]
    (C : (n : ℕ) → ι n → Set X) (s : (n : ℕ) → ι n → X)
    (hC : ∀ n i, MeasurableSet (C n i))
    (hd : ∀ n, Pairwise (fun i j => Disjoint (C n i) (C n j)))
    (hc : ∀ n, ∀ᵐ x ∂μ, ∃ i, x ∈ C n i)
    (H : ℕ → ℝ) (hH : Tendsto H atTop (𝓝 0))
    (he : ∀ n i x, x ∈ C n i → dist (s n i) x ≤ H n)
    (A : Set X) (hA : MeasurableSet A) (hfront : μ (frontier A) = 0) :
    Tendsto (fun n => ∑ i, μ.real (C n i) * A.indicator (fun _ => (1 : ℝ)) (s n i))
      atTop (𝓝 (μ.real A)) := by
  classical
  have hmeas (n : ℕ) : AEStronglyMeasurable (sampled (C n) (s n) A) μ := by
    unfold sampled
    convert (Finset.aestronglyMeasurable_sum Finset.univ fun i _ =>
      (aestronglyMeasurable_const : AEStronglyMeasurable
        (fun _ : X => A.indicator (fun _ => (1 : ℝ)) (s n i)) μ).indicator (hC n i)) using 1
    ext x
    simp
  have hnorm (n : ℕ) : ∀ᵐ x ∂μ, ‖sampled (C n) (s n) A x‖ ≤ (1 : ℝ) := by
    filter_upwards [hc n] with x hx
    obtain ⟨i,hi⟩ := hx
    rw [sampled_eq _ _ _ (hd n) hi]
    by_cases ha : s n i ∈ A <;> simp [ha]
  have hlim : ∀ᵐ x ∂μ, Tendsto (fun n => sampled (C n) (s n) A x) atTop
      (𝓝 (A.indicator (fun _ => (1 : ℝ)) x)) := by
    filter_upwards [ae_all_iff.mpr hc, compl_mem_ae_iff.mpr hfront] with x hx hxf
    let b : (n : ℕ) → ι n := fun n => Classical.choose (hx n)
    have hb (n : ℕ) : x ∈ C n (b n) := Classical.choose_spec (hx n)
    have ht : Tendsto (fun n => s n (b n)) atTop (𝓝 x) := by
      rw [tendsto_iff_dist_tendsto_zero]
      exact squeeze_zero (fun _ => dist_nonneg) (fun n => he n (b n) x (hb n)) hH
    have hi : ContinuousAt (A.indicator (fun _ => (1 : ℝ))) x :=
      continuous_const.continuousOn.continuousAt_indicator hxf
    convert hi.tendsto.comp ht using 1
    ext n
    exact sampled_eq _ _ _ (hd n) (hb n)
  have ht := tendsto_integral_of_dominated_convergence (fun _ => (1 : ℝ))
    hmeas (integrable_const _) hnorm hlim
  have heq (n : ℕ) : ∫ x, sampled (C n) (s n) A x ∂μ =
      ∑ i, μ.real (C n i) * A.indicator (fun _ => (1 : ℝ)) (s n i) := by
    unfold sampled
    rw [integral_finset_sum]
    · apply Finset.sum_congr rfl
      intro i hi
      rw [integral_indicator (hC n i), setIntegral_const]
      rfl
    · intro i hi
      exact (integrable_const _).indicator (hC n i)
  simpa only [heq, integral_indicator hA, setIntegral_const, smul_eq_mul, mul_one] using ht


open OPH.GoldenSourceAssignment

/-- The actual golden coordinates, in ordinary Lebesgue coordinates. -/
noncomputable def goldenPosition {n : ℕ} (L : ℝ) (b : Fin 3 → Fin (Nat.fib n)) :
    Fin 3 → ℝ := WithLp.ofLp (site L b)

theorem cube_measure_finite (L : ℝ) : IsFiniteMeasure (volume.restrict (cube L)) := by
  apply isFiniteMeasure_restrict.mpr
  unfold cube
  rw [Real.volume_pi_Ico]
  exact ENNReal.prod_ne_top (fun _ _ => ENNReal.ofReal_ne_top)

/-- Count the actual golden population labels whose injective readback is in A. -/
noncomputable def goldenCount (n : ℕ) (L : ℝ) (A : Set (Fin 3 → ℝ)) : ℕ := by
  classical
  exact (Finset.univ.filter (fun b : Fin 3 → Fin (Nat.fib n) => goldenPosition L b ∈ A)).card

/-- Actual normalized source counts converge on every measurable continuity
set. Coordinate Lebesgue measure is restricted to the fixed cube. -/
theorem golden_count_tendsto {L : ℝ} (hL : 0 < L)
    (A : Set (Fin 3 → ℝ)) (hA : MeasurableSet A)
    (hfront : (volume.restrict (cube L)) (frontier A) = 0) :
    Tendsto (fun n => (L ^ 3 / (Nat.fib (n + 2) : ℝ) ^ 3) *
      (goldenCount (n + 2) L A : ℝ))
      atTop (𝓝 ((volume.restrict (cube L)).real A)) := by
  classical
  let μ : Measure (Fin 3 → ℝ) := volume.restrict (cube L)
  letI : IsFiniteMeasure μ := cube_measure_finite L
  have hcube : MeasurableSet (cube L) := MeasurableSet.univ_pi fun _ => measurableSet_Ico
  have hsub (n : ℕ) (b : Fin 3 → Fin (Nat.fib (n + 2))) :
      cell (n + 2) L b ⊆ cube L := by
    rw [← cell_cover (by omega : 0 < n + 2) hL]
    exact subset_iUnion (cell (n + 2) L) b
  have hmass (n : ℕ) (b : Fin 3 → Fin (Nat.fib (n + 2))) :
      μ.real (cell (n + 2) L b) = L ^ 3 / (Nat.fib (n + 2) : ℝ) ^ 3 := by
    change (volume.restrict (cube L)).real _ = _
    rw [measureReal_restrict_apply (cell_measurable _ _ _), inter_eq_left.mpr (hsub n b)]
    exact cell_mass (by omega) hL.le b
  have hcover (n : ℕ) : ∀ᵐ x ∂μ, ∃ b, x ∈ cell (n + 2) L b := by
    filter_upwards [ae_restrict_mem hcube] with x hx
    rw [← cell_cover (by omega : 0 < n + 2) hL] at hx
    exact mem_iUnion.mp hx
  have hfib : Tendsto (fun n => (Nat.fib (n + 2) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp Nat.fib_add_two_strictMono.tendsto_atTop
  have ht := partition_count_tendsto μ (fun n => cell (n + 2) L)
    (fun _ => goldenPosition L) (fun n => cell_measurable (n + 2) L)
    (fun n => cell_disjoint (by omega : 0 < n + 2) hL) hcover
    (fun n => 2 * L / Nat.fib (n + 2)) (hfib.const_div_atTop (2 * L))
    (fun n b x hx => ?_) A hA hfront
  · convert ht using 1
    ext n
    simp only [goldenCount, hmass, Set.indicator_apply, mul_ite, mul_one, mul_zero]
    rw [← Finset.sum_filter]
    simp [nsmul_eq_mul, mul_comm]
  · apply (dist_pi_le_iff (by positivity : 0 ≤ 2 * L / (Nat.fib (n + 2) : ℝ))).mpr
    intro i
    have hxi := hx i (mem_univ i)
    simpa only [goldenPosition, site, WithLp.ofLp_toLp, Real.dist_eq] using
      coordinate_assignment (by omega : 0 < n + 2) (b i).isLt hL.le hxi.1 hxi.2.le

abbrev timeIndex (T δ : ℝ) := Fin (⌊T / δ⌋₊ + 1)

def timeCell (T δ : ℝ) (j : timeIndex T δ) : Set ℝ :=
  Ico ((j.val : ℝ) * δ) (min (((j.val : ℝ) + 1) * δ) T)

theorem timeCell_subset {T δ : ℝ} (hδ : 0 < δ) (j : timeIndex T δ) :
    timeCell T δ j ⊆ Ico 0 T := by
  intro t ht
  exact ⟨(mul_nonneg (Nat.cast_nonneg _) hδ.le).trans ht.1,
    ht.2.trans_le (min_le_right _ _)⟩

theorem timeCell_cover {T δ : ℝ} (hδ : 0 < δ) {t : ℝ} (ht : t ∈ Ico 0 T) :
    ∃ j : timeIndex T δ, t ∈ timeCell T δ j := by
  let k : ℕ := ⌊t / δ⌋₊
  have hkle : k ≤ ⌊T / δ⌋₊ := Nat.floor_mono (div_le_div_of_nonneg_right ht.2.le hδ.le)
  let j : timeIndex T δ := ⟨k, Nat.lt_succ_of_le hkle⟩
  refine ⟨j, ?_, lt_min ?_ ht.2⟩
  · have h := Nat.floor_le (div_nonneg ht.1 hδ.le)
    exact (le_div_iff₀ hδ).mp h
  · have h := Nat.lt_floor_add_one (t / δ)
    exact (div_lt_iff₀ hδ).mp h

theorem timeCell_disjoint {T δ : ℝ} (hδ : 0 < δ) :
    Pairwise (fun i j : timeIndex T δ => Disjoint (timeCell T δ i) (timeCell T δ j)) := by
  intro i j hij
  apply disjoint_left.mpr
  intro t hi hj
  rcases lt_or_gt_of_ne (Fin.val_ne_of_ne hij) with hlt | hgt
  · have hh : ((i.val : ℝ) + 1) * δ ≤ (j.val : ℝ) * δ := by
      apply mul_le_mul_of_nonneg_right _ hδ.le
      exact_mod_cast hlt
    exact (not_lt_of_ge hj.1) (hi.2.trans_le ((min_le_left _ _).trans hh))
  · have hh : ((j.val : ℝ) + 1) * δ ≤ (i.val : ℝ) * δ := by
      apply mul_le_mul_of_nonneg_right _ hδ.le
      exact_mod_cast hgt
    exact (not_lt_of_ge hi.1) (hj.2.trans_le ((min_le_left _ _).trans hh))

theorem timeCell_assignment {T δ : ℝ} (j : timeIndex T δ) {t : ℝ}
    (ht : t ∈ timeCell T δ j) : dist ((j.val : ℝ) * δ) t ≤ δ := by
  rw [Real.dist_eq, abs_of_nonpos (sub_nonpos.mpr ht.1)]
  have hh := ht.2.trans_le (min_le_left _ _)
  nlinarith



/-- Actual event cells, with the last temporal cell clipped at the fixed window. -/
noncomputable def eventCell (n : ℕ) (T L δ : ℝ)
    (i : timeIndex T δ × (Fin 3 → Fin (Nat.fib n))) : Set (ℝ × (Fin 3 → ℝ)) :=
  timeCell T δ i.1 ×ˢ cell n L i.2

noncomputable def eventPosition {n : ℕ} {T : ℝ} (L δ : ℝ)
    (i : timeIndex T δ × (Fin 3 → Fin (Nat.fib n))) : ℝ × (Fin 3 → ℝ) :=
  ((i.1.val : ℝ) * δ, goldenPosition L i.2)

noncomputable def spaceTimeVolume (T L : ℝ) : Measure (ℝ × (Fin 3 → ℝ)) :=
  (volume.restrict (Ico 0 T)).prod (volume.restrict (cube L))

theorem eventCell_measurable (n : ℕ) (T L δ : ℝ) (i) :
    MeasurableSet (eventCell n T L δ i) :=
  measurableSet_Ico.prod (cell_measurable n L i.2)

theorem eventCell_disjoint {n : ℕ} (hn : 0 < n) {T L δ : ℝ}
    (hL : 0 < L) (hδ : 0 < δ) :
    Pairwise (fun i j => Disjoint (eventCell n T L δ i) (eventCell n T L δ j)) := by
  intro i j hij
  apply disjoint_left.mpr
  intro x hi hj
  by_cases h : i.1 = j.1
  · have hh : i.2 ≠ j.2 := fun hh => hij (Prod.ext h hh)
    exact disjoint_left.mp (cell_disjoint hn hL hh) hi.2 hj.2
  · exact disjoint_left.mp (timeCell_disjoint hδ h) hi.1 hj.1

theorem eventCell_cover {n : ℕ} (hn : 0 < n) {T L δ : ℝ}
    (hL : 0 < L) (hδ : 0 < δ) :
    ⋃ i, eventCell n T L δ i = Ico 0 T ×ˢ cube L := by
  ext x
  constructor
  · intro hx
    obtain ⟨i, hi⟩ := mem_iUnion.mp hx
    refine ⟨timeCell_subset hδ i.1 hi.1, ?_⟩
    rw [← cell_cover hn hL]
    exact mem_iUnion.mpr ⟨i.2, hi.2⟩
  · intro hx
    obtain ⟨j,hj⟩ := timeCell_cover hδ hx.1
    have hxs := hx.2
    rw [← cell_cover hn hL] at hxs
    obtain ⟨b,hb⟩ := mem_iUnion.mp hxs
    exact mem_iUnion.mpr ⟨(j,b),hj,hb⟩

theorem eventCell_assignment {n : ℕ} (hn : 0 < n) {T L δ : ℝ}
    (hL : 0 < L) (hδ : 0 < δ) (i) {x : ℝ × (Fin 3 → ℝ)}
    (hx : x ∈ eventCell n T L δ i) :
    dist (eventPosition L δ i) x ≤ δ + 2 * L / Nat.fib n := by
  have hs : dist (goldenPosition L i.2) x.2 ≤ 2 * L / Nat.fib n := by
    apply (dist_pi_le_iff (by positivity)).mpr
    intro j
    have hj := hx.2 j (mem_univ j)
    simpa only [goldenPosition, site, WithLp.ofLp_toLp, Real.dist_eq] using
      coordinate_assignment hn (i.2 j).isLt hL.le hj.1 hj.2.le
  change max (dist ((i.1.val : ℝ) * δ) x.1) (dist (goldenPosition L i.2) x.2) ≤ _
  exact max_le ((timeCell_assignment i.1 hx.1).trans (le_add_of_nonneg_right (by positivity)))
    (hs.trans (le_add_of_nonneg_left hδ.le))

theorem spaceTimeVolume_finite (T L : ℝ) : IsFiniteMeasure (spaceTimeVolume T L) := by
  letI : IsFiniteMeasure (volume.restrict (cube L)) := cube_measure_finite L
  letI : IsFiniteMeasure (volume.restrict (Ico 0 T)) :=
    isFiniteMeasure_restrict.mpr (by simp [Real.volume_Ico])
  unfold spaceTimeVolume
  infer_instance

theorem event_weighted_count_tendsto {T L : ℝ} (hL : 0 < L)
    (δ : ℕ → ℝ) (hδ : ∀ n, 0 < δ n) (hδ0 : Tendsto δ atTop (𝓝 0))
    (A : Set (ℝ × (Fin 3 → ℝ))) (hA : MeasurableSet A)
    (hfront : spaceTimeVolume T L (frontier A) = 0) :
    Tendsto (fun n => ∑ i, (spaceTimeVolume T L).real (eventCell (n + 2) T L (δ n) i) *
      A.indicator (fun _ => (1 : ℝ)) (eventPosition L (δ n) i))
      atTop (𝓝 ((spaceTimeVolume T L).real A)) := by
  letI : IsFiniteMeasure (spaceTimeVolume T L) := spaceTimeVolume_finite T L
  have hcov (n : ℕ) : ∀ᵐ x ∂spaceTimeVolume T L, ∃ i, x ∈ eventCell (n + 2) T L (δ n) i := by
    have hh : ∀ᵐ x ∂spaceTimeVolume T L, x ∈ Ico 0 T ×ˢ cube L := by
      unfold spaceTimeVolume
      rw [Measure.prod_restrict]
      exact ae_restrict_mem (measurableSet_Ico.prod (MeasurableSet.univ_pi fun _ => measurableSet_Ico))
    filter_upwards [hh] with x hx
    rw [← eventCell_cover (by omega : 0 < n + 2) hL (hδ n)] at hx
    exact mem_iUnion.mp hx
  have hfib : Tendsto (fun n => (Nat.fib (n + 2) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp Nat.fib_add_two_strictMono.tendsto_atTop
  exact partition_count_tendsto (spaceTimeVolume T L)
    (fun n => eventCell (n + 2) T L (δ n)) (fun n => eventPosition L (δ n))
    (fun n => eventCell_measurable (n + 2) T L (δ n))
    (fun n => eventCell_disjoint (by omega : 0 < n + 2) hL (hδ n)) hcov
    (fun n => δ n + 2 * L / Nat.fib (n + 2))
    (by simpa using hδ0.add (hfib.const_div_atTop (2 * L)))
    (fun n i x hx => eventCell_assignment (by omega : 0 < n + 2) hL (hδ n) i hx)
    A hA hfront


theorem cube_volume_real {L : ℝ} (hL : 0 ≤ L) : volume.real (cube L) = L ^ 3 := by
  unfold cube Measure.real
  rw [Real.volume_pi_Ico_toReal (fun _ => hL)]
  simp

theorem timeCell_mass_le {T δ : ℝ} (hδ : 0 < δ) (j : timeIndex T δ) :
    volume.real (timeCell T δ j) ≤ δ := by
  unfold timeCell
  rw [Real.volume_real_Ico]
  apply max_le _ hδ.le
  have h := min_le_left (((j.val : ℝ) + 1) * δ) T
  nlinarith

theorem eventCell_mass_le {n : ℕ} (hn : 0 < n) {T L δ : ℝ}
    (hL : 0 < L) (hδ : 0 < δ) (i) :
    (spaceTimeVolume T L).real (eventCell n T L δ i) ≤
      δ * (L ^ 3 / (Nat.fib n : ℝ) ^ 3) := by
  have hb : cell n L i.2 ⊆ cube L := by
    rw [← cell_cover hn hL]
    exact subset_iUnion (cell n L) i.2
  unfold spaceTimeVolume eventCell
  rw [measureReal_prod_prod,
    measureReal_restrict_apply (show MeasurableSet (timeCell T δ i.1) from measurableSet_Ico),
    inter_eq_left.mpr (timeCell_subset hδ i.1),
    measureReal_restrict_apply (cell_measurable n L i.2), inter_eq_left.mpr hb,
    cell_mass hn hL.le i.2]
  exact mul_le_mul_of_nonneg_right (timeCell_mass_le hδ i.1) (by positivity)

theorem sum_eventCell_mass {n : ℕ} (hn : 0 < n) {T L δ : ℝ}
    (hT : 0 ≤ T) (hL : 0 < L) (hδ : 0 < δ) :
    (∑ i, (spaceTimeVolume T L).real (eventCell n T L δ i)) = T * L ^ 3 := by
  letI : IsFiniteMeasure (spaceTimeVolume T L) := spaceTimeVolume_finite T L
  rw [← measureReal_iUnion_fintype (eventCell_disjoint hn hL hδ) (eventCell_measurable n T L δ),
    eventCell_cover hn hL hδ]
  unfold spaceTimeVolume
  rw [measureReal_prod_prod, measureReal_restrict_apply measurableSet_Ico, inter_self,
    measureReal_restrict_apply (show MeasurableSet (cube L) from
      MeasurableSet.univ_pi fun _ => measurableSet_Ico), inter_self,
    Real.volume_real_Ico, sub_zero, max_eq_left hT, cube_volume_real hL.le]

theorem sum_uniform_event_weight {n : ℕ} (hn : 0 < n) (T L δ : ℝ) :
    (∑ _ : timeIndex T δ × (Fin 3 → Fin (Nat.fib n)),
      δ * (L ^ 3 / (Nat.fib n : ℝ) ^ 3)) =
      ((⌊T / δ⌋₊ : ℝ) + 1) * δ * L ^ 3 := by
  have hq : (Nat.fib n : ℝ) ≠ 0 := by exact_mod_cast (Nat.fib_pos.mpr hn).ne'
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_prod,
    Fintype.card_fin, Fintype.card_fun, nsmul_eq_mul, Nat.cast_mul, Nat.cast_pow,
    Nat.cast_add, Nat.cast_one]
  field_simp

/-- Equal event weights overcount the clipped last cell by at most one layer. -/
theorem event_uniform_weight_error {n : ℕ} (hn : 0 < n) {T L δ : ℝ}
    (hT : 0 ≤ T) (hL : 0 < L) (hδ : 0 < δ) (A : Set (ℝ × (Fin 3 → ℝ))) :
    |(∑ i : timeIndex T δ × (Fin 3 → Fin (Nat.fib n)), (δ * (L ^ 3 / (Nat.fib n : ℝ) ^ 3)) *
        A.indicator (fun _ => (1 : ℝ)) (eventPosition (T := T) L δ i)) -
      ∑ i, (spaceTimeVolume T L).real (eventCell n T L δ i) *
        A.indicator (fun _ => (1 : ℝ)) (eventPosition L δ i)| ≤ L ^ 3 * δ := by
  classical
  rw [← Finset.sum_sub_distrib]
  simp_rw [← sub_mul]
  calc
    _ ≤ ∑ i, |(δ * (L ^ 3 / (Nat.fib n : ℝ) ^ 3) -
        (spaceTimeVolume T L).real (eventCell n T L δ i)) *
        A.indicator (fun _ => (1 : ℝ)) (eventPosition L δ i)| :=
      Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, (δ * (L ^ 3 / (Nat.fib n : ℝ) ^ 3) -
        (spaceTimeVolume T L).real (eventCell n T L δ i)) := by
      apply Finset.sum_le_sum
      intro i hi
      have hnon : 0 ≤ δ * (L ^ 3 / (Nat.fib n : ℝ) ^ 3) -
          (spaceTimeVolume T L).real (eventCell n T L δ i) :=
        sub_nonneg.mpr (eventCell_mass_le hn hL hδ i)
      by_cases ha : eventPosition (T := T) L δ i ∈ A <;> simp [ha, abs_of_nonneg hnon, hnon]
    _ = (((⌊T / δ⌋₊ : ℝ) + 1) * δ - T) * L ^ 3 := by
      rw [Finset.sum_sub_distrib, sum_uniform_event_weight hn,
        sum_eventCell_mass hn hT hL hδ]
      ring
    _ ≤ L ^ 3 * δ := by
      have hf := (le_div_iff₀ hδ).mp (Nat.floor_le (div_nonneg hT hδ.le))
      have hh : ((⌊T / δ⌋₊ : ℝ) + 1) * δ - T ≤ δ := by nlinarith
      nlinarith [mul_le_mul_of_nonneg_right hh (pow_nonneg hL.le 3)]


/-- Inclusive layers 0 through floor(T/δ), with no discarded endpoint events. -/
noncomputable def eventCount (n : ℕ) (T L δ : ℝ)
    (A : Set (ℝ × (Fin 3 → ℝ))) : ℕ := by
  classical
  exact (Finset.univ.filter (fun i : timeIndex T δ × (Fin 3 → Fin (Nat.fib n)) =>
    eventPosition L δ i ∈ A)).card

/-- Actual raw event counts, with equal full-layer weights, converge on a
fixed continuity set. The last temporal cell discrepancy is derived. -/
theorem event_count_tendsto {T L : ℝ} (hT : 0 ≤ T) (hL : 0 < L)
    (δ : ℕ → ℝ) (hδ : ∀ n, 0 < δ n) (hδ0 : Tendsto δ atTop (𝓝 0))
    (A : Set (ℝ × (Fin 3 → ℝ))) (hA : MeasurableSet A)
    (hfront : spaceTimeVolume T L (frontier A) = 0) :
    Tendsto (fun n => (δ n * (L ^ 3 / (Nat.fib (n + 2) : ℝ) ^ 3)) *
      (eventCount (n + 2) T L (δ n) A : ℝ))
      atTop (𝓝 ((spaceTimeVolume T L).real A)) := by
  classical
  let U (n : ℕ) : ℝ := ∑ i : timeIndex T (δ n) × (Fin 3 → Fin (Nat.fib (n + 2))),
    (δ n * (L ^ 3 / (Nat.fib (n + 2) : ℝ) ^ 3)) *
      A.indicator (fun _ => (1 : ℝ)) (eventPosition L (δ n) i)
  let W (n : ℕ) : ℝ := ∑ i, (spaceTimeVolume T L).real (eventCell (n + 2) T L (δ n) i) *
      A.indicator (fun _ => (1 : ℝ)) (eventPosition L (δ n) i)
  have hw : Tendsto W atTop (𝓝 ((spaceTimeVolume T L).real A)) :=
    event_weighted_count_tendsto hL δ hδ hδ0 A hA hfront
  have hd : Tendsto (fun n => U n - W n) atTop (𝓝 0) := by
    apply tendsto_iff_norm_sub_tendsto_zero.mpr
    simp only [sub_zero, Real.norm_eq_abs]
    apply squeeze_zero (fun _ => abs_nonneg _)
      (fun n => event_uniform_weight_error (by omega : 0 < n + 2) hT hL (hδ n) A)
    simpa using (tendsto_const_nhds.mul hδ0 :
      Tendsto (fun n => L ^ 3 * δ n) atTop (𝓝 (L ^ 3 * 0)))
  have hu : Tendsto U atTop (𝓝 ((spaceTimeVolume T L).real A)) := by
    simpa only [sub_add_cancel, zero_add] using hd.add hw
  convert hu using 1
  ext n
  unfold U eventCount
  simp only [Set.indicator_apply, mul_ite, mul_one, mul_zero]
  rw [← Finset.sum_filter]
  simp [nsmul_eq_mul, mul_comm]

/-- The declared source-net layer duration a_q/c with a_q=L/sqrt(q). -/
noncomputable def sourceDelta (n : ℕ) (L c : ℝ) : ℝ :=
  (L / Real.sqrt (Nat.fib n)) / c

theorem sourceDelta_pos {n : ℕ} (hn : 0 < n) {L c : ℝ} (hL : 0 < L) (hc : 0 < c) :
    0 < sourceDelta n L c := by
  unfold sourceDelta
  exact div_pos (div_pos hL (Real.sqrt_pos.mpr (by exact_mod_cast Nat.fib_pos.mpr hn))) hc

theorem sourceDelta_tendsto (L c : ℝ) :
    Tendsto (fun n => sourceDelta (n + 2) L c) atTop (𝓝 0) := by
  have hfib : Tendsto (fun n => (Nat.fib (n + 2) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp Nat.fib_add_two_strictMono.tendsto_atTop
  have h := ((Real.tendsto_sqrt_atTop.comp hfib).const_div_atTop L).div_const c
  simpa only [sourceDelta, zero_div] using h

/-- Normalized counts for the actual golden source-net time grid and site
positions. The causal read law is not needed for this sampling theorem. -/
theorem source_event_count_tendsto {T L c : ℝ} (hT : 0 ≤ T) (hL : 0 < L) (hc : 0 < c)
    (A : Set (ℝ × (Fin 3 → ℝ))) (hA : MeasurableSet A)
    (hfront : spaceTimeVolume T L (frontier A) = 0) :
    Tendsto (fun n =>
      (sourceDelta (n + 2) L c * (L ^ 3 / (Nat.fib (n + 2) : ℝ) ^ 3)) *
        (eventCount (n + 2) T L (sourceDelta (n + 2) L c) A : ℝ))
      atTop (𝓝 ((spaceTimeVolume T L).real A)) :=
  event_count_tendsto hT hL _ (fun n => sourceDelta_pos (by omega : 0 < n + 2) hL hc)
    (sourceDelta_tendsto L c) A hA hfront


/-- Label counts are genuine point counts: the golden coordinate map is injective. -/
theorem goldenPosition_injective {n : ℕ} {L : ℝ} (hL : 0 < L) :
    Function.Injective (goldenPosition (n := n) L) := by
  intro b d h
  apply site_injective hL
  simpa only [goldenPosition, WithLp.toLp_ofLp] using congrArg (WithLp.toLp 2) h

theorem eventPosition_injective {n : ℕ} {T L δ : ℝ} (hL : 0 < L) (hδ : 0 < δ) :
    Function.Injective (eventPosition (n := n) (T := T) L δ) := by
  intro i j h
  apply Prod.ext
  · apply Fin.ext
    have he := congrArg Prod.fst h
    have hh : (i.1.val : ℝ) = (j.1.val : ℝ) := mul_right_cancel₀ hδ.ne' he
    exact_mod_cast hh
  · exact goldenPosition_injective hL (congrArg Prod.snd h)

theorem goldenCount_univ (n : ℕ) (L : ℝ) : goldenCount n L univ = Nat.fib n ^ 3 := by
  simp [goldenCount]

theorem eventCount_univ (n : ℕ) (T L δ : ℝ) :
    eventCount n T L δ univ = (⌊T / δ⌋₊ + 1) * Nat.fib n ^ 3 := by
  simp [eventCount]

/-- The actual inclusive endpoint is counted even when its clipped cell is empty. -/
theorem aligned_endpoint_layer_control :
    eventCount 2 1 1 1 univ = 2 ∧
      timeCell 1 1 ⟨1, by norm_num [timeIndex]⟩ = ∅ := by
  constructor
  · norm_num [eventCount_univ]
  · norm_num [timeCell]

end
end OPH.GoldenSourceCountLimit

#print axioms OPH.GoldenSourceCountLimit.partition_count_tendsto
#print axioms OPH.GoldenSourceCountLimit.golden_count_tendsto

#print axioms OPH.GoldenSourceCountLimit.event_uniform_weight_error
#print axioms OPH.GoldenSourceCountLimit.source_event_count_tendsto

#print axioms OPH.GoldenSourceCountLimit.eventPosition_injective
#print axioms OPH.GoldenSourceCountLimit.aligned_endpoint_layer_control
