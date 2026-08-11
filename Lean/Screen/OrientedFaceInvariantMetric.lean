import OrientedFaceBracketSelector

/-!
# Carrier-induced invariant-metric phase diagram for the face discriminator

Proof-kernel companion to
`code/b14_jacobi/invariant_metric_phase.certificate.json`.  The independently
replayed Python certificate derives, from the pinned oriented-face bracket and
Reynolds channel data, the exact closed-form squared distances from the face
bracket to the classified compact families `P`, `F`, `G` under **every**
invariant carrier inner product.  The pinned carrier splits multiplicity-free
into sectors `1 + 3 + 3' + 5`, so by the certificate's commutant-dimension
receipt the complete family of invariant inner products is the positive
sector-scale cone `(alpha, beta, gamma, delta)`; the induced Hilbert--Schmidt
distances never involve `alpha` and are the three-term Laurent forms below.

This file takes those closed forms as literal definitions and **proves the
phase diagram as real quantified theorems**: `P` is strictly excluded for
every invariant metric; every sector-balanced metric (`gamma = beta`) selects
`G`; every metric with `delta <= 50*beta` and `beta <= 6*delta` selects `G`
regardless of `gamma`; and the mirror family `F` genuinely takes over at the
exact witness `(8, 1, 1)`.  The reference point `(1,1,1)` reproduces the
pinned selector distances, tying the two certificate layers together.

**Boundary.**  The closed forms themselves are certificate content: their
derivation from the pinned tensors is the replayed Python computation, not a
Lean theorem.  The comparison ranges over the compact locus classified by the
pinned `b14_compact_locus` certificate, conditional on its three named
textbook lemmas.  The nearest-point rule and the restriction to
carrier-induced metrics are declared discriminator choices.  No norm,
bracket, current, or physical gauge structure is source-selected here, and no
laboratory identification is made.
-/

namespace OPH.OrientedFaceInvariantMetric

open OPH.OrientedFaceBracketSelector

noncomputable def s5 : ℝ := Real.sqrt 5

private lemma hs5_pos : (0 : ℝ) < s5 := by
  simp [s5, Real.sqrt_pos]

private lemma hs5_sq : s5 ^ 2 = 5 := by
  simp [s5, Real.sq_sqrt]

/-- Squared distance from the face bracket to the compact plane `P` under the
sector-scale metric `(beta, gamma, delta)`; `alpha` never enters. -/
noncomputable def dP2 (β γ δ : ℝ) : ℝ :=
  (15 - 3*s5)/β + (15 + 3*s5)/γ + ((15 - 3*s5)/2) * β/δ^2 + ((15 + 3*s5)/2) * γ/δ^2

/-- Squared distance from the face bracket to the compact `su(3)+so(3)` family
`F` (the mirror cell) under the sector-scale metric. -/
noncomputable def dF2 (β γ δ : ℝ) : ℝ :=
  ((15 + 3*s5)/2) * γ/δ^2 + (15 + 3*s5)/γ + (60 + 12*s5)/(11*β)

/-- Squared distance from the face bracket to the compact `su(3)+so(3)` family
`G` under the sector-scale metric. -/
noncomputable def dG2 (β γ δ : ℝ) : ℝ :=
  ((15 - 3*s5)/2) * β/δ^2 + (15 - 3*s5)/β + (60 - 12*s5)/(11*γ)

/-- Parameterized family distance in one display. -/
noncomputable def squaredDistanceAt (β γ δ : ℝ) : CompactFamily → ℝ
  | .P => dP2 β γ δ
  | .F => dF2 β γ δ
  | .G => dG2 β γ δ

/-! ## Reference-point calibration against the pinned selector distances -/

theorem reference_P : dP2 1 1 1 = distanceP := by
  simp only [dP2, distanceP]; ring

theorem reference_F : dF2 1 1 1 = distanceF := by
  simp only [dF2, distanceF, s5]; ring

theorem reference_G : dG2 1 1 1 = distanceG := by
  simp only [dG2, distanceG, s5]; ring

/-! ## Galois structure of the phase constants

The `F` and `G` coefficient tables are `sqrt 5 |-> -sqrt 5` conjugates of one
another under the sector swap `beta <-> gamma`.  At the constant level the
conjugate pairs multiply to rationals: the quadratic-slope pair to `45`, the
reciprocal pair to `(30/11)^2`. -/

theorem gap_coefficient_norm :
    ((15 + 3*s5)/2) * ((15 - 3*s5)/2) = 45 := by
  have h := hs5_sq; nlinarith [h]

theorem reciprocal_coefficient_norm :
    ((105 + 45*s5)/11) * ((105 - 45*s5)/11) = 900/121 := by
  have h := hs5_sq; nlinarith [h]

/-! ## P is excluded for every invariant carrier metric -/

theorem dG2_lt_dP2 {β γ δ : ℝ} (hβ : 0 < β) (hγ : 0 < γ) (hδ : 0 < δ) :
    dG2 β γ δ < dP2 β γ δ := by
  have hs := hs5_pos
  have hs2 := hs5_sq
  have hgap : dP2 β γ δ - dG2 β γ δ
      = ((105 + 45*s5)/11)/γ + ((15 + 3*s5)/2) * γ/δ^2 := by
    simp only [dP2, dG2]
    field_simp
    ring
  have h1 : 0 < ((105 + 45*s5)/11)/γ := by positivity
  have h2 : 0 < ((15 + 3*s5)/2) * γ/δ^2 := by positivity
  linarith

theorem dF2_lt_dP2 {β γ δ : ℝ} (hβ : 0 < β) (hγ : 0 < γ) (hδ : 0 < δ) :
    dF2 β γ δ < dP2 β γ δ := by
  have hs := hs5_pos
  have hs2 := hs5_sq
  have h3 : 0 < 105 - 45*s5 := by nlinarith
  have h4 : 0 < 15 - 3*s5 := by nlinarith
  have hgap : dP2 β γ δ - dF2 β γ δ
      = ((105 - 45*s5)/11)/β + ((15 - 3*s5)/2) * β/δ^2 := by
    simp only [dP2, dF2]
    field_simp
    ring
  have h1 : 0 < ((105 - 45*s5)/11)/β := by positivity
  have h2 : 0 < ((15 - 3*s5)/2) * β/δ^2 := by positivity
  linarith

/-! ## Every sector-balanced invariant metric selects G -/

theorem balanced_dG2_lt_dF2 {β δ : ℝ} (hβ : 0 < β) (hδ : 0 < δ) :
    dG2 β β δ < dF2 β β δ := by
  have hs := hs5_pos
  have hgap : dF2 β β δ - dG2 β β δ = (3*s5) * β/δ^2 + ((90/11)*s5)/β := by
    simp only [dF2, dG2]
    field_simp
    ring
  have h1 : 0 < (3*s5) * β/δ^2 := by positivity
  have h2 : 0 < ((90/11)*s5)/β := by positivity
  linarith

/-- On every sector-balanced invariant carrier metric, `G` is the unique
nearest classified compact family. -/
theorem balanced_unique_nearest_G {β δ : ℝ} (hβ : 0 < β) (hδ : 0 < δ)
    (family : CompactFamily) (h : family ≠ .G) :
    dG2 β β δ < squaredDistanceAt β β δ family := by
  cases family with
  | P => exact dG2_lt_dP2 hβ hβ hδ
  | F => exact balanced_dG2_lt_dF2 hβ hδ
  | G => exact False.elim (h rfl)

/-! ## The phase box: bounded sector-3⁺ chirality keeps G selected -/

private lemma k_branch_poly {β δ : ℝ} (hβ : 0 < β) (hδ : 0 < δ)
    (h1 : δ ≤ 50*β) (h2 : β ≤ 6*δ) :
    ((15 - 3*s5)/2)*β^2 + ((105 - 45*s5)/11)*δ^2 < 28*β*δ := by
  have hs := hs5_pos
  have hs2 := hs5_sq
  have hprod : 0 ≤ (50*β - δ) * (6*δ - β) := by
    apply mul_nonneg <;> linarith
  have h301 : 0 < 903*s5 - 1715 := by nlinarith
  have hC : (105 - 45*s5)/11 - 3*((15 - 3*s5)/2)/25 < 0 := by nlinarith
  nlinarith [mul_pos hβ hδ, mul_pos hδ hδ, mul_pos hβ hβ,
    mul_nonneg (le_of_lt (mul_pos hβ hδ)) (le_of_lt h301)]

private lemma h_branch_poly {γ δ : ℝ} (_hγ : 0 < γ) (_hδ : 0 < δ) :
    28*γ*δ ≤ ((15 + 3*s5)/2)*γ^2 + ((105 + 45*s5)/11)*δ^2 := by
  have hs := hs5_pos
  have hs2 := hs5_sq
  have haF : 0 < (15 + 3*s5)/2 := by positivity
  have hdisc : 0 < 4*((15 + 3*s5)/2)*((105 + 45*s5)/11) - 784 := by nlinarith
  nlinarith [sq_nonneg (2*((15 + 3*s5)/2)*γ - 28*δ), mul_pos _hδ _hδ, haF]

/-- For every invariant carrier metric whose sector-3⁺ scale satisfies
`delta/50 <= beta <= 6*delta`, the family `G` is strictly nearer than `F`
for **all** `gamma` and `delta`. -/
theorem box_dG2_lt_dF2 {β γ δ : ℝ} (hβ : 0 < β) (hγ : 0 < γ) (hδ : 0 < δ)
    (h1 : δ ≤ 50*β) (h2 : β ≤ 6*δ) :
    dG2 β γ δ < dF2 β γ δ := by
  have hs := hs5_pos
  have hδ2 : (0:ℝ) < δ^2 := by positivity
  have hgap : dF2 β γ δ - dG2 β γ δ
      = (((15 + 3*s5)/2)*γ/δ^2 + ((105 + 45*s5)/11)/γ)
        - (((15 - 3*s5)/2)*β/δ^2 + ((105 - 45*s5)/11)/β) := by
    simp only [dF2, dG2]
    field_simp
    ring
  have hk : ((15 - 3*s5)/2)*β/δ^2 + ((105 - 45*s5)/11)/β < 28/δ := by
    rw [div_add_div _ _ (by positivity : (δ^2 : ℝ) ≠ 0) (ne_of_gt hβ),
      div_lt_div_iff₀ (by positivity) hδ]
    have := k_branch_poly hβ hδ h1 h2
    nlinarith [this, mul_pos hβ hδ]
  have hh : 28/δ ≤ ((15 + 3*s5)/2)*γ/δ^2 + ((105 + 45*s5)/11)/γ := by
    rw [div_add_div _ _ (by positivity : (δ^2 : ℝ) ≠ 0) (ne_of_gt hγ),
      div_le_div_iff₀ hδ (by positivity)]
    have := h_branch_poly hγ hδ
    nlinarith [this, mul_pos hγ hδ]
  linarith

/-- Unique nearest family on the whole phase box, for all `gamma`, `delta`. -/
theorem box_unique_nearest_G {β γ δ : ℝ} (hβ : 0 < β) (hγ : 0 < γ) (hδ : 0 < δ)
    (h1 : δ ≤ 50*β) (h2 : β ≤ 6*δ)
    (family : CompactFamily) (h : family ≠ .G) :
    dG2 β γ δ < squaredDistanceAt β γ δ family := by
  cases family with
  | P => exact dG2_lt_dP2 hβ hγ hδ
  | F => exact box_dG2_lt_dF2 hβ hγ hδ h1 h2
  | G => exact False.elim (h rfl)

/-! ## The F region is nonempty: F wins at `(8, 1, 1)` -/

theorem F_wins_at_witness : dF2 8 1 1 < dG2 8 1 1 := by
  have hs := hs5_pos
  have hs2 := hs5_sq
  have hgap : dG2 8 1 1 - dF2 8 1 1 = (3885 - 1593*s5)/88 := by
    simp only [dG2, dF2]
    norm_num
    ring
  nlinarith [hgap]

#print axioms reference_G
#print axioms dG2_lt_dP2
#print axioms dF2_lt_dP2
#print axioms balanced_unique_nearest_G
#print axioms box_unique_nearest_G
#print axioms F_wins_at_witness

end OPH.OrientedFaceInvariantMetric
