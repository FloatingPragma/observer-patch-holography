import Mathlib.Data.Finset.Prod
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-! # Correlated collar words on a finite spatial admission register

Two uniformly sampled distinct native cells carry the ordered words. Conditional
on seeing both words in a permitted correlation region, their measured charge
contrast is constant. These results account for finite admission counts and
for the cross term when measuring the union of two disjoint regions. The
routing length and any physical stress interpretation are separate inputs.
-/

namespace OPH.EinsteinBranch.CollarAdmission

open Filter
open scoped Topology BigOperators

variable {α : Type*}

/-- The weight includes the conditional modular contrast and the inverse number
of ordered distinct placements in the full host. -/
noncomputable def pairCharge (weight : ℝ) (region : Finset α) : ℝ :=
  weight * (region.offDiag.card : ℝ)

/-- Count all ordered distinct coincident admissions in an arbitrary region. -/
theorem pairCharge_count (w : ℝ) (s : Finset α) :
    pairCharge w s = w * s.card * ((s.card : ℝ) - 1) := by
  rw [pairCharge, Finset.offDiag_card, Nat.cast_sub (Nat.le_mul_self _), Nat.cast_mul]
  ring

variable [DecidableEq α]

/-- Mixed admissions contribute twice the product of the two disjoint populations. -/
theorem pairCharge_disjoint (w : ℝ) (s t : Finset α) (hd : Disjoint s t) :
    pairCharge w (s ∪ t) - pairCharge w s - pairCharge w t =
      2 * w * (s.card : ℝ) * (t.card : ℝ) := by
  simp only [pairCharge_count, Finset.card_union_of_disjoint hd, Nat.cast_add]
  ring

/-- Restrict tag correlation to members of the same ancestor. The union charge
is additive exactly up to the shared-ancestor cross population. -/
theorem pairCharge_ancestor_cross {ι : Type*} (groups : Finset ι) (w : ℝ)
    (s t : ι → Finset α) (hd : ∀ i ∈ groups, Disjoint (s i) (t i)) :
    (∑ i ∈ groups, pairCharge w (s i ∪ t i)) -
    (∑ i ∈ groups, pairCharge w (s i)) - (∑ i ∈ groups, pairCharge w (t i)) =
    ∑ i ∈ groups, 2 * w * ((s i).card : ℝ) * ((t i).card : ℝ) := by
  rw [← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun i hi => pairCharge_disjoint w (s i) (t i) (hd i hi)

/-- Without-replacement normalization expressed through the admitted density
per squared radius, microscopic site density, and total host count. -/
theorem quartic_admission_formula (n N r : ℝ) (hN : N ≠ 0) (hN1 : N ≠ 1) (hr : r ≠ 0) :
    n * (n - 1) / (N * (N - 1)) / r^4 =
    ((n / N / r^2)^2 - (n / N / r^2) * (1 / (N * r^2))) / (1 - 1/N) := by
  have hsub : N - 1 ≠ 0 := sub_ne_zero.mpr hN1
  have hdiff : 1 - 1/N ≠ 0 := by
    intro h
    have h' := sub_eq_zero.mp h
    field_simp [hN] at h'
    exact hN1 h'
  field_simp [hN, hsub, hr, hdiff]

/-- A quadratic one-slot admission fraction gives a fourth-order distinct-pair
fraction when the normalized inverse microscopic population vanishes. -/
theorem quartic_admission_limit (n N r : ℕ → ℝ) (a : ℝ)
    (hN : ∀ i, N i ≠ 0) (hN1 : ∀ i, N i ≠ 1) (hr : ∀ i, r i ≠ 0)
    (ha : Tendsto (fun i => n i / N i / r i ^ 2) atTop (𝓝 a))
    (hb : Tendsto (fun i => 1 / (N i * r i ^ 2)) atTop (𝓝 0))
    (hc : Tendsto (fun i => 1 / N i) atTop (𝓝 0)) :
    Tendsto (fun i => n i * (n i - 1) / (N i * (N i - 1)) / r i ^ 4) atTop (𝓝 (a^2)) := by
  simp_rw [quartic_admission_formula _ _ _ (hN _) (hN1 _) (hr _)]
  have hnum := (ha.pow 2).sub (ha.mul hb)
  have hden := (tendsto_const_nhds (x := (1 : ℝ))).sub hc
  simpa using hnum.div hden (by norm_num : (1 : ℝ) - 0 ≠ 0)

end OPH.EinsteinBranch.CollarAdmission
