import Geometry.SourceNetOrderLimit
import Geometry.SourceCountTransport
import Geometry.SourceCausalBoundary

set_option autoImplicit false

/-!
# Generated-order limits for the actual golden source population

The sites are the fractional golden orbit, not a replacement regular grid.
The grid is used only for the equal-volume cell assignment of
`Geometry.GoldenSourceAssignment`. The read relation is the complete-neighbour path relation at radius L/sqrt(q).
All event layers, including the clipped final layer, remain in the counts.
-/

namespace OPH.GoldenSourceCausalLimit

open MeasureTheory Set Filter
open OPH.GoldenSourceAssignment OPH.GoldenSourceCountLimit
open OPH.SourceNetCausalCone OPH.SourceNetOrderLimit OPH.SourceCausalBoundary
open OPH.SourceCountTransport
open scoped Topology BigOperators Classical

noncomputable section

abbrev Carrier := EuclideanSpace ℝ (Fin 3)

def domain (L : ℝ) : Set Carrier := {x | ∀ i, 0 ≤ x i ∧ x i < L}
def population (n : ℕ) (L : ℝ) : Set Carrier := Set.range (site (n := n) L)
def radius (n : ℕ) (L : ℝ) : ℝ := L / Real.sqrt (Nat.fib n)
def assignment (n : ℕ) (L : ℝ) : ℝ := 2 * Real.sqrt 3 * L / Nat.fib n

theorem domain_convex (L : ℝ) : Convex ℝ (domain L) := by
  intro x hx y hy a b ha hb hab i
  simpa only [PiLp.add_apply, PiLp.smul_apply, smul_eq_mul] using
    (convex_Ico (0 : ℝ) L (hx i) (hy i) ha hb hab)

theorem population_subset {n : ℕ} {L : ℝ} (hL : 0 < L) :
    population n L ⊆ domain L := by
  rintro x ⟨b, rfl⟩ i
  change 0 ≤ L * orbit (b i).val ∧ L * orbit (b i).val < L
  constructor
  · exact mul_nonneg hL.le (Int.fract_nonneg _)
  · simpa using mul_lt_mul_of_pos_left (Int.fract_lt_one _) hL

/-- The proved permutation partition gives a covering by actual source sites. -/
theorem population_covers {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 < L) :
    Covers (domain L) (population n L) (assignment n L) := by
  intro x hx
  have hxc : WithLp.ofLp x ∈ cube L := fun i _ => hx i
  rw [← cell_cover hn hL] at hxc
  obtain ⟨b, hb⟩ := mem_iUnion.mp hxc
  refine ⟨site L b, ⟨b, rfl⟩, ?_⟩
  rw [← dist_eq_norm]
  exact tensor_assignment hn hL.le b x (fun i =>
    ⟨(hb i (mem_univ i)).1, (hb i (mem_univ i)).2.le⟩)

theorem radius_pos {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 < L) :
    0 < radius n L := by
  exact div_pos hL (Real.sqrt_pos.mpr (by exact_mod_cast Nat.fib_pos.mpr hn))

theorem assignment_nonneg (n : ℕ) {L : ℝ} (hL : 0 ≤ L) :
    0 ≤ assignment n L := by unfold assignment; positivity

theorem assignment_radius_ratio {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 < L) :
    assignment n L / radius n L = 2 * Real.sqrt 3 / Real.sqrt (Nat.fib n) := by
  have hq : (0 : ℝ) < Nat.fib n := by exact_mod_cast Nat.fib_pos.mpr hn
  have hs : Real.sqrt (Nat.fib n) ≠ 0 := (Real.sqrt_pos.mpr hq).ne'
  have hsq := Real.sq_sqrt hq.le
  unfold assignment radius
  field_simp
  nlinarith

theorem assignment_ratio_tendsto {L : ℝ} (hL : 0 < L) :
    Tendsto (fun n => assignment (n + 2) L / radius (n + 2) L) atTop (𝓝 0) := by
  have hfib : Tendsto (fun n => (Nat.fib (n + 2) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp Nat.fib_add_two_strictMono.tendsto_atTop
  have heq (n : ℕ) := assignment_radius_ratio (by omega : 0 < n + 2) hL
  simp_rw [heq]
  exact (Real.tendsto_sqrt_atTop.comp hfib).const_div_atTop (2 * Real.sqrt 3)

abbrev Label (n : ℕ) (T L c : ℝ) :=
  timeIndex T (sourceDelta n L c) × (Fin 3 → Fin (Nat.fib n))

def placed {n : ℕ} {T L c : ℝ} (i : Label n T L c) : ℕ × population n L :=
  (i.1.val, ⟨site L i.2, ⟨i.2, rfl⟩⟩)

def generated {n : ℕ} {T L c : ℝ} (i j : Label n T L c) : Prop :=
  Precedes (radius n L) (placed i) (placed j)

def strictlyGenerated {n : ℕ} {T L c : ℝ} (i j : Label n T L c) : Prop :=
  i.1.val < j.1.val ∧ generated i j

def position {n : ℕ} {T L c : ℝ} (i : Label n T L c) : Spacetime :=
  eventPosition L (sourceDelta n L c) i

theorem placed_time {n : ℕ} {T L c : ℝ} (i : Label n T L c) :
    layerTime (radius n L) c (placed i).1 = (position i).1 := rfl

theorem placed_space {n : ℕ} {T L c : ℝ} (i : Label n T L c) :
    ((placed i).2 : Carrier) = WithLp.toLp 2 (position i).2 := by
  simp [placed, position, eventPosition, goldenPosition]

/-- Eventual membership in the generated order, with moving sampled events.
The null-separated case is explicitly excluded, rather than assumed reflected. -/
theorem generated_eventually_iff {T L c : ℝ} (hL : 0 < L) (hc : 0 < c)
    (e f : (n : ℕ) → Label (n + 2) T L c) {p q : Spacetime}
    (he : Tendsto (fun n => position (e n)) atTop (𝓝 p))
    (hf : Tendsto (fun n => position (f n)) atTop (𝓝 q))
    (hn : margin c p q ≠ 0) :
    ∀ᶠ n in atTop, (generated (e n) (f n) ↔ 0 ≤ margin c p q) ∧
      (strictlyGenerated (e n) (f n) ↔ 0 ≤ margin c p q) := by
  have het : Tendsto (fun n => layerTime (radius (n + 2) L) c (placed (e n)).1)
      atTop (𝓝 p.1) := by simpa only [placed_time] using (continuous_fst.tendsto p).comp he
  have hft : Tendsto (fun n => layerTime (radius (n + 2) L) c (placed (f n)).1)
      atTop (𝓝 q.1) := by simpa only [placed_time] using (continuous_fst.tendsto q).comp hf
  have hex : Tendsto (fun n => ((placed (e n)).2 : Carrier)) atTop
      (𝓝 (WithLp.toLp 2 p.2)) := by
    simpa only [placed_space] using ((PiLp.continuous_toLp 2 _).tendsto p.2).comp
      ((continuous_snd.tendsto p).comp he)
  have hfx : Tendsto (fun n => ((placed (f n)).2 : Carrier)) atTop
      (𝓝 (WithLp.toLp 2 q.2)) := by
    simpa only [placed_space] using ((PiLp.continuous_toLp 2 _).tendsto q.2).comp
      ((continuous_snd.tendsto q).comp hf)
  have hnorm : ‖(WithLp.toLp 2 q.2 : Carrier) - WithLp.toLp 2 p.2‖ =
      spatialNorm (q.2 - p.2) := rfl
  rcases lt_or_gt_of_ne hn with hn | hn
  · have hsp : c * (q.1 - p.1) <
        ‖(WithLp.toLp 2 q.2 : Carrier) - WithLp.toLp 2 p.2‖ := by
      rw [hnorm]
      exact sub_neg.mp hn
    filter_upwards [eventually_not_precedes_of_spacelike
      (fun n => population (n + 2) L) (fun n => radius (n + 2) L) hc
      (fun n => placed (e n)) (fun n => placed (f n)) het hft hex hfx hsp] with n hno
    simp [generated, strictlyGenerated, hno, hn.not_ge]
  · have hti : ‖(WithLp.toLp 2 q.2 : Carrier) - WithLp.toLp 2 p.2‖ <
        c * (q.1 - p.1) := by
      rw [hnorm]
      exact sub_pos.mp hn
    filter_upwards [eventually_precedes_of_timelike (domain_convex L)
      (fun n => population (n + 2) L) (fun n => radius (n + 2) L)
      (fun n => assignment (n + 2) L) hc
      (fun _ => radius_pos (by omega) hL) (fun _ => assignment_nonneg _ hL.le)
      (fun _ => population_subset hL) (fun _ => population_covers (by omega) hL)
      (assignment_ratio_tendsto hL) (fun n => placed (e n)) (fun n => placed (f n))
      het hft hex hfx hti] with n hy
    exact ⟨iff_of_true hy.2 hn.le, iff_of_true hy hn.le⟩

def intervalMember {n : ℕ} {T L c : ℝ} (e f i : Label n T L c) : Prop :=
  generated e i ∧ generated i f

def intervalCount {n : ℕ} {T L c : ℝ} (e f : Label n T L c) : ℕ :=
  count (intervalMember e f)

/-- Every choice of the actual assigned representatives converges to its point. -/
theorem assigned_position_tendsto {T L c : ℝ} (hL : 0 < L) (hc : 0 < c)
    {x : Spacetime} (b : (n : ℕ) → Label (n + 2) T L c)
    (hb : ∀ n, x ∈ eventCell (n + 2) T L (sourceDelta (n + 2) L c) (b n)) :
    Tendsto (fun n => position (b n)) atTop (𝓝 x) := by
  have hfib : Tendsto (fun n => (Nat.fib (n + 2) : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp Nat.fib_add_two_strictMono.tendsto_atTop
  rw [tendsto_iff_dist_tendsto_zero]
  apply squeeze_zero (fun _ => dist_nonneg)
    (fun n => eventCell_assignment (by omega : 0 < n + 2) hL
      (sourceDelta_pos (by omega) hL hc) (b n) (hb n))
  simpa using (sourceDelta_tendsto L c).add (hfib.const_div_atTop (2 * L))

theorem interval_classification_ae {T L c : ℝ} (hL : 0 < L) (hc : 0 < c)
    (e f : (n : ℕ) → Label (n + 2) T L c) {p q : Spacetime}
    (he : Tendsto (fun n => position (e n)) atTop (𝓝 p))
    (hf : Tendsto (fun n => position (f n)) atTop (𝓝 q)) :
    ∀ᵐ x ∂spaceTimeVolume T L, ∀ᶠ n in atTop, ∀ i,
      x ∈ eventCell (n + 2) T L (sourceDelta (n + 2) L c) i →
        (intervalMember (e n) (f n) i ↔ x ∈ diamond c p q) := by
  have hcov (n : ℕ) : ∀ᵐ x ∂spaceTimeVolume T L, ∃ i,
      x ∈ eventCell (n + 2) T L (sourceDelta (n + 2) L c) i :=
    eventCell_ae_cover (by omega : 0 < n + 2) hL (sourceDelta_pos (by omega) hL hc)
  filter_upwards [ae_all_iff.mpr hcov, compl_mem_ae_iff.mpr (future_null T L hc.ne' p),
    compl_mem_ae_iff.mpr (past_null T L hc.ne' q)] with x hx hpx hxq
  let b : (n : ℕ) → Label (n + 2) T L c := fun n => Classical.choose (hx n)
  have hb (n : ℕ) : x ∈ eventCell (n + 2) T L (sourceDelta (n + 2) L c) (b n) :=
    Classical.choose_spec (hx n)
  have hbt := assigned_position_tendsto hL hc b hb
  filter_upwards [generated_eventually_iff hL hc e b he hbt hpx,
    generated_eventually_iff hL hc b f hbt hf hxq] with n hleft hright
  intro i hi
  have hieq : i = b n := by
    by_contra hn
    exact disjoint_left.mp (eventCell_disjoint (by omega : 0 < n + 2) hL
      (sourceDelta_pos (by omega) hL hc) hn) hi (hb n)
  subst i
  exact hleft.1.and hright.1

/-- Moving generated intervals converge in the actual counting measure.
The reference volume is explicitly clipped to the supplied time/space box. -/
theorem generated_interval_count_tendsto {T L c : ℝ} (hT : 0 ≤ T) (hL : 0 < L)
    (hc : 0 < c) (e f : (n : ℕ) → Label (n + 2) T L c) {p q : Spacetime}
    (he : Tendsto (fun n => position (e n)) atTop (𝓝 p))
    (hf : Tendsto (fun n => position (f n)) atTop (𝓝 q)) :
    Tendsto (fun n =>
      (sourceDelta (n + 2) L c * (L ^ 3 / (Nat.fib (n + 2) : ℝ) ^ 3)) *
        (intervalCount (e n) (f n) : ℝ))
      atTop (𝓝 ((spaceTimeVolume T L).real (diamond c p q))) :=
  event_selected_count_tendsto hT hL _ (fun n => sourceDelta_pos (by omega) hL hc)
    (sourceDelta_tendsto L c) (fun n => intervalMember (e n) (f n)) _
    (diamond_measurable c p q) (interval_classification_ae hL hc e f he hf)

#print axioms population_covers
#print axioms assignment_ratio_tendsto
#print axioms generated_eventually_iff
#print axioms generated_interval_count_tendsto

end
end OPH.GoldenSourceCausalLimit
