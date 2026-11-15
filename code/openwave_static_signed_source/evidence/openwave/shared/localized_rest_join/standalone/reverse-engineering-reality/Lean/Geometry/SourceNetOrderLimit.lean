import Geometry.SourceNetLayeredOrder
import Mathlib.Analysis.SpecificLimits.Basic

set_option autoImplicit false

/-!
# Stability of the generated source-net order away from the null cone

The relation below is the actual all-neighbour path relation, proved equal in
`Geometry.SourceNetLayeredOrder` to the reflexive transitive closure of the
layered reads. A supplied
covering, read radius and model clock are retained as hypotheses. Event
positions and times may move. The conclusion is eventual exact agreement
away from the null cone, not a count limit or a measure assumption.
-/

namespace OPH.SourceNetOrderLimit

open Filter Set OPH.SourceNetCausalCone
open scoped Topology

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- Model time of an actual layer, with the supplied radius/speed clock. -/
noncomputable def layerTime (a c : ℝ) (j : ℕ) : ℝ := (j : ℝ) * (a / c)

theorem layer_gap {a c : ℝ} (hc : c ≠ 0) {j k : ℕ} (hjk : j ≤ k) :
    ((k - j : ℕ) : ℝ) * a = c * (layerTime a c k - layerTime a c j) := by
  rw [Nat.cast_sub hjk]
  unfold layerTime
  field_simp

theorem inner_layer_gap {a c h : ℝ} (ha : a ≠ 0) (hc : c ≠ 0)
    {j k : ℕ} (hjk : j ≤ k) :
    ((k - j : ℕ) : ℝ) * (a - 2 * h) =
      c * (layerTime a c k - layerTime a c j) * (1 - 2 * h / a) := by
  rw [Nat.cast_sub hjk]
  unfold layerTime
  field_simp

/-- A strictly timelike limiting pair eventually has an actual read path.
No path-existence or order-convergence premise is supplied. -/
theorem eventually_precedes_of_timelike
    {Ω : Set E} (hΩ : Convex ℝ Ω) (S : ℕ → Set E)
    (a h : ℕ → ℝ) {c : ℝ} (hc : 0 < c)
    (ha : ∀ n, 0 < a n) (hh : ∀ n, 0 ≤ h n)
    (hS : ∀ n, S n ⊆ Ω) (hcover : ∀ n, Covers Ω (S n) (h n))
    (hratio : Tendsto (fun n => h n / a n) atTop (𝓝 0))
    (e f : (n : ℕ) → ℕ × S n) {t u : ℝ} {x y : E}
    (het : Tendsto (fun n => layerTime (a n) c (e n).1) atTop (𝓝 t))
    (hft : Tendsto (fun n => layerTime (a n) c (f n).1) atTop (𝓝 u))
    (hex : Tendsto (fun n => ((e n).2 : E)) atTop (𝓝 x))
    (hfx : Tendsto (fun n => ((f n).2 : E)) atTop (𝓝 y))
    (htime : ‖y - x‖ < c * (u - t)) :
    ∀ᶠ n in atTop, (e n).1 < (f n).1 ∧ Precedes (a n) (e n) (f n) := by
  have htu : t < u := by nlinarith [norm_nonneg (y - x)]
  have htlt := het.eventually_lt hft htu
  have hdist := (hfx.sub hex).norm
  have hinner : Tendsto
      (fun n => c * (layerTime (a n) c (f n).1 - layerTime (a n) c (e n).1) *
        (1 - 2 * h n / a n)) atTop (𝓝 (c * (u - t))) := by
    have hs : Tendsto (fun n => 1 - 2 * (h n / a n)) atTop (𝓝 (1 - 2 * 0)) :=
      tendsto_const_nhds.sub (tendsto_const_nhds.mul hratio)
    have hp := ((hft.sub het).const_mul c).mul hs
    simpa only [mul_zero, sub_zero, mul_one, mul_div_assoc] using hp
  filter_upwards [htlt, hdist.eventually_lt hinner htime] with n hn hnear
  have hjk : (e n).1 < (f n).1 := by
    have hd : 0 < a n / c := div_pos (ha n) hc
    have hn' : ((e n).1 : ℝ) < ((f n).1 : ℝ) :=
      (mul_lt_mul_iff_left₀ hd).mp hn
    exact_mod_cast hn'
  refine ⟨hjk, hjk.le, covering_constructs_path hΩ (hS n) (hcover n) (hh n)
    (e n).2.property (f n).2.property (Nat.sub_pos_of_lt hjk) ?_⟩
  rw [inner_layer_gap (ha n).ne' hc.ne' hjk.le]
  exact hnear.le

omit [NormedSpace ℝ E] in
/-- A strictly spacelike or backwards limiting pair eventually has no path. -/
theorem eventually_not_precedes_of_spacelike
    (S : ℕ → Set E) (a : ℕ → ℝ) {c : ℝ} (hc : 0 < c)
    (e f : (n : ℕ) → ℕ × S n) {t u : ℝ} {x y : E}
    (het : Tendsto (fun n => layerTime (a n) c (e n).1) atTop (𝓝 t))
    (hft : Tendsto (fun n => layerTime (a n) c (f n).1) atTop (𝓝 u))
    (hex : Tendsto (fun n => ((e n).2 : E)) atTop (𝓝 x))
    (hfx : Tendsto (fun n => ((f n).2 : E)) atTop (𝓝 y))
    (hspace : c * (u - t) < ‖y - x‖) :
    ∀ᶠ n in atTop, ¬ Precedes (a n) (e n) (f n) := by
  have hd := (hfx.sub hex).norm
  have ht := (hft.sub het).const_mul c
  filter_upwards [ht.eventually_lt hd hspace] with n hn
  intro hp
  have ho := reachable_outer hp.2
  rw [layer_gap hc.ne' hp.1] at ho
  exact (not_lt_of_ge ho) hn

#print axioms eventually_precedes_of_timelike
#print axioms eventually_not_precedes_of_spacelike

end OPH.SourceNetOrderLimit
