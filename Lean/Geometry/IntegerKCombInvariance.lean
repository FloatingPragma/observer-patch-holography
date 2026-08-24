import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecificLimits.Basic

/-!
# Exact invariance receipts for the integer-k ringdown comb template

Build-stage companion to the numeric instrument in
`code/gravitation/ringdown_comb/`, for the frozen target
`falsification/frozen_targets/fz01_2026-07-17/frozen_target_integer_k_comb_2026-07-17.md`.
The abstract template has teeth `f_k = f0 + s * ln k` for natural
`k >= 2`, scale `s > 0`, and offset `f0`. In the frozen target's
declared physical reading, `s` is `G M / (c^3 g(chi))` of a Kerr
remnant scaled by the KMS bookkeeping and carries the entire mass,
spin, and redshift dependence, while `f0` is the rotation line
`m Omega_H / (2 pi)`.

Proved here, exactly and for every template: the offset-subtracted
ratio of any two teeth equals `ln a / ln b` and is therefore identical
across all templates (mass, spin, and redshift invariance of the
frozen observable); teeth are strictly monotone in `k` and sit above
the offset for `k >= 2`; the reference ladder against the `k = 2`
tooth satisfies exact rational brackets, `158/100 < ln 3 / ln 2 <
159/100`, `ln 4 / ln 2 = 2` exactly, and `232/100 < ln 5 / ln 2 <
233/100`, each bracket reduced to a natural-number power inequality;
and the KMS weight hierarchy `(k-1)/k` is strictly increasing, bounded
by one, and converges to one.

What is not proved here. The template's physical reading, that the
scale belongs to a Kerr remnant through `G M / (c^3 g(chi))` and the
offset to its rotation line, is the frozen target's declared
identification, not derived in this module. No event data enters; no
strain likelihood, prior normalization, or decision rule is
formalized; and nothing here is a registered, frozen, or scored
prediction. The join from this abstract template to a derived
strain-level likelihood is open.
-/

namespace OPH.IntegerKCombInvariance

open Filter Real

/-- Abstract integer-k comb template: positive scale `s`, offset `f0`.
In the declared physical reading the scale carries all mass, spin, and
redshift dependence and the offset is the rotation line. -/
structure CombTemplate where
  scale : ℝ
  offset : ℝ
  scale_pos : 0 < scale

/-- Explicit nontrivial inhabitant with concrete scale and offset. -/
noncomputable def referenceTemplate : CombTemplate :=
  { scale := 2, offset := 3, scale_pos := by norm_num }

noncomputable instance : Inhabited CombTemplate := ⟨referenceTemplate⟩

/-- Tooth `k` of a template: `f_k = f0 + s * ln k`. -/
noncomputable def tooth (T : CombTemplate) (k : ℕ) : ℝ :=
  T.offset + T.scale * Real.log (k : ℝ)

/-- The offset-subtracted tooth is the pure logarithmic part. -/
theorem tooth_sub_offset (T : CombTemplate) (k : ℕ) :
    tooth T k - T.offset = T.scale * Real.log (k : ℝ) := by
  simp [tooth]

/-- The frozen ratio law at template level: the offset-subtracted
ratio of any two teeth is `ln a / ln b`, with the scale cancelled. -/
theorem offsetSubtracted_ratio (T : CombTemplate) (a b : ℕ) :
    (tooth T a - T.offset) / (tooth T b - T.offset)
      = Real.log (a : ℝ) / Real.log (b : ℝ) := by
  have hs : T.scale ≠ 0 := ne_of_gt T.scale_pos
  rw [tooth_sub_offset, tooth_sub_offset, mul_div_mul_left _ _ hs]

/-- Mass, spin, and redshift invariance of the frozen observable: the
offset-subtracted ratio is identical for every template, hence for
every remnant mass, spin, and redshift in the declared reading, since
the scale carries all of that dependence and cancels. -/
theorem ratio_template_independent (T₁ T₂ : CombTemplate) (a b : ℕ) :
    (tooth T₁ a - T₁.offset) / (tooth T₁ b - T₁.offset)
      = (tooth T₂ a - T₂.offset) / (tooth T₂ b - T₂.offset) := by
  rw [offsetSubtracted_ratio, offsetSubtracted_ratio]

/-- Strict monotonicity of teeth in `k`. -/
theorem tooth_lt_tooth (T : CombTemplate) {a b : ℕ}
    (ha : 1 ≤ a) (hab : a < b) : tooth T a < tooth T b := by
  have ha0 : (0 : ℝ) < (a : ℝ) := by
    exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one ha
  have hcast : (a : ℝ) < (b : ℝ) := by exact_mod_cast hab
  have hlog : Real.log (a : ℝ) < Real.log (b : ℝ) :=
    Real.log_lt_log ha0 hcast
  have hmul := mul_lt_mul_of_pos_left hlog T.scale_pos
  unfold tooth
  linarith

/-- Every tooth with `k >= 2` sits strictly above the offset (the
rotation line of the declared reading). -/
theorem offset_lt_tooth (T : CombTemplate) {k : ℕ} (hk : 2 ≤ k) :
    T.offset < tooth T k := by
  have h1 : (1 : ℝ) < (k : ℝ) := by
    have : 1 < k := by omega
    exact_mod_cast this
  have hlog : 0 < Real.log (k : ℝ) := Real.log_pos h1
  have hmul := mul_pos T.scale_pos hlog
  unfold tooth
  linarith

theorem log_two_pos : (0 : ℝ) < Real.log 2 :=
  Real.log_pos (by norm_num)

/-- Natural-number reduction of the lower `ln 3 / ln 2` bracket. -/
theorem pow_bound_32_lower : (2 : ℕ) ^ 158 < 3 ^ 100 := by norm_num

/-- Natural-number reduction of the upper `ln 3 / ln 2` bracket. -/
theorem pow_bound_32_upper : (3 : ℕ) ^ 100 < 2 ^ 159 := by norm_num

/-- Natural-number reduction of the lower `ln 5 / ln 2` bracket. -/
theorem pow_bound_52_lower : (2 : ℕ) ^ 232 < 5 ^ 100 := by norm_num

/-- Natural-number reduction of the upper `ln 5 / ln 2` bracket. -/
theorem pow_bound_52_upper : (5 : ℕ) ^ 100 < 2 ^ 233 := by norm_num

/-- `158 ln 2 < 100 ln 3`, from `2^158 < 3^100`. -/
theorem log_ladder_32_lower :
    158 * Real.log 2 < 100 * Real.log 3 := by
  have hpow : (2 : ℝ) ^ (158 : ℕ) < (3 : ℝ) ^ (100 : ℕ) := by
    exact_mod_cast pow_bound_32_lower
  have h := Real.log_lt_log (by positivity) hpow
  rw [Real.log_pow, Real.log_pow] at h
  exact_mod_cast h

/-- `100 ln 3 < 159 ln 2`, from `3^100 < 2^159`. -/
theorem log_ladder_32_upper :
    100 * Real.log 3 < 159 * Real.log 2 := by
  have hpow : (3 : ℝ) ^ (100 : ℕ) < (2 : ℝ) ^ (159 : ℕ) := by
    exact_mod_cast pow_bound_32_upper
  have h := Real.log_lt_log (by positivity) hpow
  rw [Real.log_pow, Real.log_pow] at h
  exact_mod_cast h

/-- `232 ln 2 < 100 ln 5`, from `2^232 < 5^100`. -/
theorem log_ladder_52_lower :
    232 * Real.log 2 < 100 * Real.log 5 := by
  have hpow : (2 : ℝ) ^ (232 : ℕ) < (5 : ℝ) ^ (100 : ℕ) := by
    exact_mod_cast pow_bound_52_lower
  have h := Real.log_lt_log (by positivity) hpow
  rw [Real.log_pow, Real.log_pow] at h
  exact_mod_cast h

/-- `100 ln 5 < 233 ln 2`, from `5^100 < 2^233`. -/
theorem log_ladder_52_upper :
    100 * Real.log 5 < 233 * Real.log 2 := by
  have hpow : (5 : ℝ) ^ (100 : ℕ) < (2 : ℝ) ^ (233 : ℕ) := by
    exact_mod_cast pow_bound_52_upper
  have h := Real.log_lt_log (by positivity) hpow
  rw [Real.log_pow, Real.log_pow] at h
  exact_mod_cast h

/-- Exact rational bracket for the `k = 3` ladder entry:
`158/100 < ln 3 / ln 2 < 159/100`. -/
theorem log3_div_log2_bounds :
    (158 : ℝ) / 100 < Real.log 3 / Real.log 2 ∧
      Real.log 3 / Real.log 2 < (159 : ℝ) / 100 := by
  constructor
  · rw [div_lt_div_iff₀ (by norm_num) log_two_pos]
    have h := log_ladder_32_lower
    linarith
  · rw [div_lt_div_iff₀ log_two_pos (by norm_num)]
    have h := log_ladder_32_upper
    linarith

/-- The `k = 4` ladder entry is exactly two: `ln 4 / ln 2 = 2`. -/
theorem log4_div_log2 : Real.log 4 / Real.log 2 = 2 := by
  have h4 : (4 : ℝ) = 2 ^ (2 : ℕ) := by norm_num
  rw [h4, Real.log_pow]
  push_cast
  exact mul_div_cancel_right₀ 2 (ne_of_gt log_two_pos)

/-- Exact rational bracket for the `k = 5` ladder entry:
`232/100 < ln 5 / ln 2 < 233/100`. -/
theorem log5_div_log2_bounds :
    (232 : ℝ) / 100 < Real.log 5 / Real.log 2 ∧
      Real.log 5 / Real.log 2 < (233 : ℝ) / 100 := by
  constructor
  · rw [div_lt_div_iff₀ (by norm_num) log_two_pos]
    have h := log_ladder_52_lower
    linarith
  · rw [div_lt_div_iff₀ log_two_pos (by norm_num)]
    have h := log_ladder_52_upper
    linarith

/-- Ladder bracket at tooth level: the offset-subtracted `3 : 2` tooth
ratio of every template lies in `(158/100, 159/100)`. -/
theorem tooth_ratio_32_bounds (T : CombTemplate) :
    (158 : ℝ) / 100 < (tooth T 3 - T.offset) / (tooth T 2 - T.offset) ∧
      (tooth T 3 - T.offset) / (tooth T 2 - T.offset) < (159 : ℝ) / 100 := by
  have h := offsetSubtracted_ratio T 3 2
  push_cast at h
  rw [h]
  exact log3_div_log2_bounds

/-- Ladder identity at tooth level: the offset-subtracted `4 : 2`
tooth ratio of every template is exactly two. -/
theorem tooth_ratio_42 (T : CombTemplate) :
    (tooth T 4 - T.offset) / (tooth T 2 - T.offset) = 2 := by
  have h := offsetSubtracted_ratio T 4 2
  push_cast at h
  rw [h, log4_div_log2]

/-- Ladder bracket at tooth level: the offset-subtracted `5 : 2` tooth
ratio of every template lies in `(232/100, 233/100)`. -/
theorem tooth_ratio_52_bounds (T : CombTemplate) :
    (232 : ℝ) / 100 < (tooth T 5 - T.offset) / (tooth T 2 - T.offset) ∧
      (tooth T 5 - T.offset) / (tooth T 2 - T.offset) < (233 : ℝ) / 100 := by
  have h := offsetSubtracted_ratio T 5 2
  push_cast at h
  rw [h]
  exact log5_div_log2_bounds

/-- KMS detailed-balance weight `(k-1)/k` of the frozen secondary
structure. -/
noncomputable def kmsWeight (k : ℕ) : ℝ := ((k : ℝ) - 1) / (k : ℝ)

/-- The `k = 2` weight is one half. -/
theorem kmsWeight_two : kmsWeight 2 = 1 / 2 := by
  norm_num [kmsWeight]

/-- The KMS weight hierarchy is strictly increasing. -/
theorem kmsWeight_strictMono {a b : ℕ} (ha : 1 ≤ a) (hab : a < b) :
    kmsWeight a < kmsWeight b := by
  have ha0 : (0 : ℝ) < (a : ℝ) := by
    exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one ha
  have hcast : (a : ℝ) < (b : ℝ) := by exact_mod_cast hab
  have hb0 : (0 : ℝ) < (b : ℝ) := lt_trans ha0 hcast
  unfold kmsWeight
  rw [div_lt_div_iff₀ ha0 hb0]
  have expandA : ((a : ℝ) - 1) * (b : ℝ) = (a : ℝ) * (b : ℝ) - (b : ℝ) := by
    ring
  have expandB : ((b : ℝ) - 1) * (a : ℝ) = (a : ℝ) * (b : ℝ) - (a : ℝ) := by
    ring
  rw [expandA, expandB]
  linarith

/-- Every KMS weight is strictly below one. -/
theorem kmsWeight_lt_one {k : ℕ} (hk : 1 ≤ k) : kmsWeight k < 1 := by
  have hk0 : (0 : ℝ) < (k : ℝ) := by
    exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hk
  unfold kmsWeight
  rw [div_lt_one hk0]
  linarith

/-- The KMS weight hierarchy converges to one. -/
theorem kmsWeight_tendsto_one :
    Tendsto kmsWeight atTop (nhds 1) := by
  have hinv : Tendsto (fun k : ℕ => ((k : ℝ))⁻¹) atTop (nhds 0) :=
    tendsto_inv_atTop_zero.comp tendsto_natCast_atTop_atTop
  have hmain : Tendsto (fun k : ℕ => 1 - ((k : ℝ))⁻¹) atTop (nhds (1 - 0)) :=
    tendsto_const_nhds.sub hinv
  rw [sub_zero] at hmain
  apply hmain.congr'
  filter_upwards [eventually_ge_atTop 1] with k hk
  have hk0 : ((k : ℝ)) ≠ 0 := by
    have : (0 : ℝ) < (k : ℝ) := by
      exact_mod_cast Nat.lt_of_lt_of_le Nat.zero_lt_one hk
    exact ne_of_gt this
  simp only [kmsWeight]
  rw [sub_div, div_self hk0, one_div]

/-- Evaluation of the inhabitant: concrete tooth of the concrete
template. -/
theorem referenceTemplate_tooth_two :
    tooth referenceTemplate 2 = 3 + 2 * Real.log 2 := by
  simp [tooth, referenceTemplate]

end OPH.IntegerKCombInvariance

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.IntegerKCombInvariance.tooth_sub_offset
#print axioms OPH.IntegerKCombInvariance.offsetSubtracted_ratio
#print axioms OPH.IntegerKCombInvariance.ratio_template_independent
#print axioms OPH.IntegerKCombInvariance.tooth_lt_tooth
#print axioms OPH.IntegerKCombInvariance.offset_lt_tooth
#print axioms OPH.IntegerKCombInvariance.log_two_pos
#print axioms OPH.IntegerKCombInvariance.pow_bound_32_lower
#print axioms OPH.IntegerKCombInvariance.pow_bound_32_upper
#print axioms OPH.IntegerKCombInvariance.pow_bound_52_lower
#print axioms OPH.IntegerKCombInvariance.pow_bound_52_upper
#print axioms OPH.IntegerKCombInvariance.log_ladder_32_lower
#print axioms OPH.IntegerKCombInvariance.log_ladder_32_upper
#print axioms OPH.IntegerKCombInvariance.log_ladder_52_lower
#print axioms OPH.IntegerKCombInvariance.log_ladder_52_upper
#print axioms OPH.IntegerKCombInvariance.log3_div_log2_bounds
#print axioms OPH.IntegerKCombInvariance.log4_div_log2
#print axioms OPH.IntegerKCombInvariance.log5_div_log2_bounds
#print axioms OPH.IntegerKCombInvariance.tooth_ratio_32_bounds
#print axioms OPH.IntegerKCombInvariance.tooth_ratio_42
#print axioms OPH.IntegerKCombInvariance.tooth_ratio_52_bounds
#print axioms OPH.IntegerKCombInvariance.kmsWeight_two
#print axioms OPH.IntegerKCombInvariance.kmsWeight_strictMono
#print axioms OPH.IntegerKCombInvariance.kmsWeight_lt_one
#print axioms OPH.IntegerKCombInvariance.kmsWeight_tendsto_one
#print axioms OPH.IntegerKCombInvariance.referenceTemplate_tooth_two
