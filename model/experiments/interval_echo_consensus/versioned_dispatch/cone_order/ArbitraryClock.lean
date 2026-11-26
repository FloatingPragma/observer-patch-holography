import Geometry.SourceDerivedSpacetimeCarrier
import Mathlib.Tactic

namespace OPH.RecordedConeArbitraryClock

abbrev Root := Fin 7
abbrev EventPoint := Fin 4 → ℝ

def score (x c : EventPoint) : ℝ :=
  (c 0 - x 0)^2 - (c 1 - x 1)^2 - (c 2 - x 2)^2 - (c 3 - x 3)^2

def future (x c : EventPoint) : Prop := x 0 ≤ c 0 ∧ 0 ≤ score x c

/-- The extra coordinate forces the weight of the latest root to vanish.
The six remaining weights annihilate all affine functions of spacetime. -/
private theorem dependence (x : Root → EventPoint) (anchor : Root) :
    ∃ w : Root → ℝ, w anchor = 0 ∧ (∃ i, w i ≠ 0) ∧
      ∀ (a : ℝ) (b : EventPoint), ∑ i, w i * (a + ∑ j, x i j * b j) = 0 := by
  classical
  let lift (i : Root) : Fin 6 → ℝ := Fin.cons (if i = anchor then 1 else 0) (Fin.cons 1 (x i))
  have dep : ¬ LinearIndependent ℝ lift := by
    intro h
    have dim := h.fintype_card_le_finrank
    simp [Root] at dim
  obtain ⟨w, total, nz⟩ := Fintype.not_linearIndependent_iff.mp dep
  have col (k : Fin 6) : ∑ i, w i * lift i k = 0 := by
    simpa [Finset.sum_apply, Pi.smul_apply] using congrFun total k
  have wa : w anchor = 0 := by simpa [lift, mul_ite] using col 0
  have mass : ∑ i, w i = 0 := by simpa [lift] using col (0 : Fin 5).succ
  have moment (j : Fin 4) : ∑ i, w i * x i j = 0 := by simpa [lift] using col j.succ.succ
  refine ⟨w, wa, nz, ?_⟩
  intro a b
  simp_rw [mul_add, Finset.mul_sum]
  rw [Finset.sum_add_distrib, ← Finset.sum_mul, mass, zero_mul, zero_add, Finset.sum_comm]
  simp_rw [← mul_assoc, ← Finset.sum_mul, moment, zero_mul, Finset.sum_const_zero]

/-- Past cones of freely timed points of 1+3 Minkowski space cannot shatter
seven fixed spacetime points. Neither root times nor sink times are prescribed. -/
theorem cannotShatter (x : Root → EventPoint) :
    ¬ ∀ S : Set Root, ∃ c : EventPoint, ∀ i, i ∈ S ↔ future (x i) c := by
  classical
  intro shatter
  obtain ⟨a, _, latest⟩ := Finset.exists_max_image Finset.univ (fun i => x i 0) Finset.univ_nonempty
  obtain ⟨w, wa, ⟨k, wk⟩, cancel⟩ := dependence x a
  obtain ⟨c, hc⟩ := shatter {i | i = a ∨ 0 < w i}
  obtain ⟨d, hd⟩ := shatter {i | i = a ∨ w i < 0}
  change ∀ i, (i = a ∨ 0 < w i) ↔ future (x i) c at hc
  change ∀ i, (i = a ∨ w i < 0) ↔ future (x i) d at hd
  have tc (i : Root) : x i 0 ≤ c 0 := (latest i (Finset.mem_univ i)).trans ((hc a).mp (.inl rfl)).1
  have td (i : Root) : x i 0 ≤ d 0 := (latest i (Finset.mem_univ i)).trans ((hd a).mp (.inl rfl)).1
  let gap (i : Root) := score (x i) c - score (x i) d
  have positive (i : Root) (hi : w i ≠ 0) : 0 < w i * gap i := by
    have ne : i ≠ a := by rintro rfl; exact hi wa
    rcases lt_or_gt_of_ne hi with hn | hp
    · have inside := ((hd i).mp (.inr hn)).2
      have outside : score (x i) c < 0 := lt_of_not_ge fun h =>
        ((hc i).mpr ⟨tc i, h⟩).elim ne hn.not_gt
      exact mul_pos_of_neg_of_neg hn (sub_neg.mpr (outside.trans_le inside))
    · have inside := ((hc i).mp (.inr hp)).2
      have outside : score (x i) d < 0 := lt_of_not_ge fun h =>
        ((hd i).mpr ⟨td i, h⟩).elim ne hp.not_gt
      exact mul_pos hp (sub_pos.mpr (outside.trans_le inside))
  have strict : 0 < ∑ i, w i * gap i := by
    apply Finset.sum_pos' _ ⟨k, Finset.mem_univ k, positive k wk⟩
    intro i _
    by_cases hi : w i = 0
    · simp [hi]
    · exact (positive i hi).le
  let offset := score 0 c - score 0 d
  let slope : EventPoint := ![2*(d 0 - c 0), 2*(c 1 - d 1), 2*(c 2 - d 2), 2*(c 3 - d 3)]
  have affine (i : Root) : gap i = offset + ∑ j, x i j * slope j := by
    simp [gap, offset, slope, score, Fin.sum_univ_succ]
    ring
  simp_rw [affine] at strict
  exact (ne_of_gt strict) (cancel _ _)

open SourceDerivedSpacetime PrimitivePortFrameQuotient Provenance

noncomputable def coord (u : SourceSpacetimeCarrier) : EventPoint :=
  ![u.1, quotientEquivVec3 u.2 0, quotientEquivVec3 u.2 1, quotientEquivVec3 u.2 2]

theorem coord_future (u v : SourceSpacetimeCarrier) :
    sourceCausalLE u v ↔ future (coord u) (coord v) := by
  change (0 ≤ (v.1-u.1)*(v.1-u.1) - quotientGram (v.2-u.2) (v.2-u.2) ∧
    0 ≤ v.1-u.1) ↔ _
  have eqGram : quotientGram (v.2-u.2) (v.2-u.2) =
      ∑ j, (quotientEquivVec3 v.2 j - quotientEquivVec3 u.2 j)^2 := by
    unfold quotientGram PrimitivePortTranslationBridge.dot
    rw [map_sub]
    simp [pow_two, Pi.sub_apply]
  rw [eqGram]
  simp [future, coord, score, Fin.sum_univ_succ, sub_nonneg, pow_two]
  constructor <;> rintro ⟨h, g⟩ <;> exact ⟨by linarith, by linarith⟩

/-- A single executed semantic log realizing these cuts has no order-faithful
map into the repository's Lorentz carrier, even with freely chosen times. -/
theorem noSpacetimeMap {R V E : Type*} [DecidableEq R] [DecidableEq V]
    [Fintype E] [DecidableEq E] (L : SemanticEventLog R V E)
    (root : Root → E) (sink : Set Root → E)
    (disjoint : ∀ i S, root i ≠ sink S)
    (cuts : ∀ i S, i ∈ S ↔ L.GeneratedBefore (root i) (sink S)) :
    ¬ ∃ place : E → SourceSpacetimeCarrier, ∀ e f,
      L.GeneratedBeforeEq e f ↔ sourceCausalLE (place e) (place f) := by
  rintro ⟨place, faithful⟩
  apply cannotShatter (fun i => coord (place (root i)))
  intro S
  refine ⟨coord (place (sink S)), fun i => ?_⟩
  have eqBefore : L.GeneratedBeforeEq (root i) (sink S) ↔ L.GeneratedBefore (root i) (sink S) :=
    or_iff_right (disjoint i S)
  exact (cuts i S).trans (eqBefore.symm.trans ((faithful _ _).trans (coord_future _ _)))

end OPH.RecordedConeArbitraryClock
