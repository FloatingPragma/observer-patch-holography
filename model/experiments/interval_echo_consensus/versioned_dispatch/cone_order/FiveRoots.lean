import Geometry.SourceDerivedSpacetimeCarrier
import Mathlib.Tactic

namespace OPH.RecordedConeFiveRoots

abbrev Point := Fin 3 → ℝ
abbrev Root := Fin 5

def sqDistance (x c : Point) : ℝ := ∑ j, (x j - c j)^2

/-- A nonzero signed weighting of the roots annihilates every affine score.
The fourth coordinate in the lift allows the score to have a constant term. -/
private theorem affineDependence (x : Root → Point) :
    ∃ w : Root → ℝ, (∃ i, w i ≠ 0) ∧
      ∀ (a : ℝ) (b : Point), ∑ i, w i * (a + ∑ j, x i j * b j) = 0 := by
  classical
  let lift (i : Root) : Fin 4 → ℝ := Fin.cons 1 (x i)
  have dep : ¬ LinearIndependent ℝ lift := by
    intro h
    have dim := h.fintype_card_le_finrank
    simp [Root] at dim
  obtain ⟨w, annihilate, nonzero⟩ := Fintype.not_linearIndependent_iff.mp dep
  have column (k : Fin 4) : ∑ i, w i * lift i k = 0 := by
    simpa [Finset.sum_apply, Pi.smul_apply] using congrFun annihilate k
  have mass : ∑ i, w i = 0 := by simpa [lift] using column 0
  have moment (j : Fin 3) : ∑ i, w i * x i j = 0 := by
    simpa [lift] using column j.succ
  refine ⟨w, nonzero, ?_⟩
  intro a b
  simp_rw [mul_add, Finset.mul_sum]
  rw [Finset.sum_add_distrib, ← Finset.sum_mul, mass, zero_mul, zero_add, Finset.sum_comm]
  simp_rw [← mul_assoc, ← Finset.sum_mul, moment, zero_mul, Finset.sum_const_zero]

/-- Even with arbitrary real squared-radius thresholds, the closed balls of
three-space cannot realize every subset of five fixed, possibly coincident roots. -/
theorem cannotShatter (x : Root → Point) :
    ¬ ∀ S : Set Root, ∃ (c : Point) (r : ℝ),
        ∀ i, i ∈ S ↔ sqDistance (x i) c ≤ r := by
  classical
  intro shatter
  obtain ⟨w, ⟨i₀, hi₀⟩, cancels⟩ := affineDependence x
  obtain ⟨c, r, hc⟩ := shatter {i | 0 < w i}
  obtain ⟨d, s, hd⟩ := shatter {i | w i < 0}
  let gap (i : Root) := (r - sqDistance (x i) c) - (s - sqDistance (x i) d)
  have positive (i : Root) (hi : w i ≠ 0) : 0 < w i * gap i := by
    change ∀ i, 0 < w i ↔ _ at hc
    change ∀ i, w i < 0 ↔ _ at hd
    rcases lt_or_gt_of_ne hi with hn | hp
    · have inD := (hd i).mp hn
      have outC := lt_of_not_ge (fun h => hn.not_gt ((hc i).mpr h))
      exact mul_pos_of_neg_of_neg hn (by dsimp [gap]; linarith)
    · have inC := (hc i).mp hp
      have outD := lt_of_not_ge (fun h => hp.not_gt ((hd i).mpr h))
      exact mul_pos hp (by dsimp [gap]; linarith)
  have strict : 0 < ∑ i, w i * gap i := by
    apply Finset.sum_pos' _ ⟨i₀, Finset.mem_univ i₀, positive i₀ hi₀⟩
    intro i _
    by_cases hi : w i = 0
    · simp [hi]
    · exact (positive i hi).le
  have affine (i : Root) :
      gap i = (r - s + ∑ j, ((d j)^2 - (c j)^2)) +
        ∑ j, x i j * (2 * (c j - d j)) := by
    simp only [gap, sqDistance, Fin.sum_univ_succ]
    ring
  simp_rw [affine] at strict
  exact (ne_of_gt strict) (cancels _ _)

/-- For simultaneous source events the existing Lorentz cone cuts out a
closed ball of squared radius t² at every future query time t. -/
theorem sourceConeBall (x c : PrimitivePortFrameQuotient.FrameQuotient)
    (t : ℝ) (ht : 0 ≤ t) :
    SourceDerivedSpacetime.sourceCausalLE (0, x) (t, c) ↔
      sqDistance (PrimitivePortFrameQuotient.quotientEquivVec3 x)
        (PrimitivePortFrameQuotient.quotientEquivVec3 c) ≤ t^2 := by
  open PrimitivePortFrameQuotient SourceDerivedSpacetime in
    change (0 ≤ (t - 0) * (t - 0) - quotientGram (c - x) (c - x) ∧ 0 ≤ t - 0) ↔ _
  have eqGram : PrimitivePortFrameQuotient.quotientGram (c - x) (c - x) =
      sqDistance (PrimitivePortFrameQuotient.quotientEquivVec3 x)
        (PrimitivePortFrameQuotient.quotientEquivVec3 c) := by
    unfold PrimitivePortFrameQuotient.quotientGram PrimitivePortTranslationBridge.dot sqDistance
    rw [map_sub]
    apply Finset.sum_congr rfl
    intro j _
    simp only [Pi.sub_apply]
    ring
  simp [eqGram, ht, sq, sub_nonneg]

/-- A different future query event is permitted for every subset. Only
the five root positions and their simultaneous time are shared across trials. -/
theorem noCommonRootConeChart (x : Root → PrimitivePortFrameQuotient.FrameQuotient) :
    ¬ ∀ S : Set Root, ∃ (c : PrimitivePortFrameQuotient.FrameQuotient) (t : ℝ),
        0 ≤ t ∧ ∀ i, i ∈ S ↔ SourceDerivedSpacetime.sourceCausalLE (0, x i) (t, c) := by
  intro realize
  apply cannotShatter (fun i => PrimitivePortFrameQuotient.quotientEquivVec3 (x i))
  intro S
  obtain ⟨c, t, ht, hc⟩ := realize S
  exact ⟨PrimitivePortFrameQuotient.quotientEquivVec3 c, t^2,
    fun i => (hc i).trans (sourceConeBall _ _ _ ht)⟩

end OPH.RecordedConeFiveRoots
