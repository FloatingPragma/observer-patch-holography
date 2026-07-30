import Mathlib

namespace ObserverPatchHolography.ScalarSeamRepair

/-!
# Scalar seam-repair algebra

This file proves the finite linear algebra used by a candidate scalar
seam-repair construction.  It does not assert that the OPH axioms select this
construction as the physical repair law.

For two endpoints, support off the endpoints, endpoint agreement after the
update, and preservation of the endpoint sum already force pair averaging.
The resulting linear map is an idempotent retraction, commutes with endpoint
exchange, and is nonidentity when the endpoints are distinct.  Uniformly
averaging these maps over a finite seam family gives the identity minus the
graph Laplacian divided by twice the number of seams.  In particular, a
thirty-seam family gives `I - L / 60`.
-/

noncomputable section

variable {ι : Type*} [DecidableEq ι]

/-- Exchange the scalar values at two endpoints and leave all other
coordinates fixed. -/
def swapEndpoints (u v : ι) (x : ι → ℝ) : ι → ℝ :=
  fun w => if w = u then x v else if w = v then x u else x w

/-- The linear update that replaces both endpoint values by their arithmetic
mean and leaves every other coordinate fixed. -/
def pairAverage (u v : ι) : (ι → ℝ) →ₗ[ℝ] (ι → ℝ) where
  toFun x w :=
    if w = u ∨ w = v then (x u + x v) / 2 else x w
  map_add' x y := by
    funext w
    by_cases h : w = u ∨ w = v
    · simp [h]
      ring
    · simp [h]
  map_smul' c x := by
    funext w
    by_cases h : w = u ∨ w = v
    · simp [h]
      ring
    · simp [h]

@[simp]
theorem pairAverage_left (u v : ι) (x : ι → ℝ) :
    pairAverage u v x u = (x u + x v) / 2 := by
  simp [pairAverage]

@[simp]
theorem pairAverage_right (u v : ι) (x : ι → ℝ) :
    pairAverage u v x v = (x u + x v) / 2 := by
  simp [pairAverage]

@[simp]
theorem pairAverage_off_seam (u v w : ι) (x : ι → ℝ)
    (hwu : w ≠ u) (hwv : w ≠ v) :
    pairAverage u v x w = x w := by
  simp [pairAverage, hwu, hwv]

/-- Pair averaging maps every input into the endpoint-agreement subspace. -/
theorem pairAverage_agrees (u v : ι) (x : ι → ℝ) :
    pairAverage u v x u = pairAverage u v x v := by
  simp

/-- Pair averaging preserves the total scalar load on the two endpoints. -/
theorem pairAverage_preserves_sum (u v : ι) (x : ι → ℝ) :
    pairAverage u v x u + pairAverage u v x v = x u + x v := by
  simp

/-- An input that already agrees at the endpoints is fixed pointwise. -/
theorem pairAverage_fixes_agreement (u v : ι) (x : ι → ℝ)
    (hagree : x u = x v) :
    pairAverage u v x = x := by
  funext w
  by_cases h : w = u ∨ w = v
  · rcases h with rfl | rfl <;> simp [hagree]
  · simp [pairAverage, h]

/-- Pair averaging is idempotent. -/
theorem pairAverage_idempotent (u v : ι) (x : ι → ℝ) :
    pairAverage u v (pairAverage u v x) = pairAverage u v x :=
  pairAverage_fixes_agreement u v (pairAverage u v x)
    (pairAverage_agrees u v x)

/-- Pair averaging commutes with exchanging the two endpoint labels. -/
theorem pairAverage_swap_equivariant (u v : ι) (x : ι → ℝ)
    (hne : u ≠ v) :
    pairAverage u v (swapEndpoints u v x) =
      swapEndpoints u v (pairAverage u v x) := by
  funext w
  by_cases hwu : w = u
  · subst w
    simp [swapEndpoints, hne.symm]
    ring
  · by_cases hwv : w = v
    · subst w
      simp [swapEndpoints, hne.symm]
      ring
    · simp [swapEndpoints, pairAverage, hwu, hwv]

/-- A genuine two-endpoint averaging update is not the identity map. -/
theorem pairAverage_ne_id (u v : ι) (hne : u ≠ v) :
    pairAverage u v ≠ LinearMap.id := by
  intro heq
  let x : ι → ℝ := fun w => if w = u then 1 else 0
  have hcoord :
      pairAverage u v x u = (LinearMap.id : (ι → ℝ) →ₗ[ℝ] (ι → ℝ)) x u := by
    rw [heq]
  norm_num [pairAverage, x, hne, hne.symm] at hcoord

/-- Support off the seam, agreement of the output endpoints, and preservation
of their sum uniquely determine pair averaging.  Linearity is carried by the
type of `T`; idempotence and swap equivariance need not be assumed because
they follow for the uniquely determined map. -/
theorem eq_pairAverage_of_supported_agreeing_sum_preserving
    (u v : ι) (T : (ι → ℝ) →ₗ[ℝ] (ι → ℝ))
    (hsupport : ∀ x w, w ≠ u → w ≠ v → T x w = x w)
    (hagrees : ∀ x, T x u = T x v)
    (hsum : ∀ x, T x u + T x v = x u + x v) :
    T = pairAverage u v := by
  ext x w
  by_cases hwu : w = u
  · subst w
    have ha := hagrees x
    have hs := hsum x
    simp only [pairAverage_left]
    linarith
  · by_cases hwv : w = v
    · subst w
      have ha := hagrees x
      have hs := hsum x
      simp only [pairAverage_right]
      linarith
    · rw [hsupport x w hwu hwv, pairAverage_off_seam u v w x hwu hwv]

/-- The single-edge graph-Laplacian contribution.  The coordinate theorem
below identifies this linear-map definition with the usual signed endpoint
difference formula. -/
def edgeLaplacian (u v : ι) : (ι → ℝ) →ₗ[ℝ] (ι → ℝ) :=
  (2 : ℝ) • (LinearMap.id - pairAverage u v)

/-- Coordinate form of the single-edge graph Laplacian. -/
theorem edgeLaplacian_apply (u v : ι) (x : ι → ℝ) (w : ι) :
    edgeLaplacian u v x w =
      (if w = u then x u - x v else 0) +
        (if w = v then x v - x u else 0) := by
  by_cases hwu : w = u
  · subst w
    by_cases huv : u = v
    · subst v
      simp [edgeLaplacian, pairAverage]
    · simp [edgeLaplacian, pairAverage, huv]
      ring
  · by_cases hwv : w = v
    · subst w
      simp [edgeLaplacian, pairAverage, hwu]
      ring
    · simp [edgeLaplacian, pairAverage, hwu, hwv]

/-- Pair averaging is the identity minus one half of the single-edge graph
Laplacian. -/
theorem pairAverage_eq_id_sub_half_edgeLaplacian (u v : ι) :
    pairAverage u v =
      LinearMap.id - (1 / 2 : ℝ) • edgeLaplacian u v := by
  rw [edgeLaplacian]
  module

variable {ε : Type*} [Fintype ε]

/-- Sum of the single-edge Laplacians for a finite indexed seam family. -/
def graphLaplacian (left right : ε → ι) :
    (ι → ℝ) →ₗ[ℝ] (ι → ℝ) :=
  ∑ e : ε, edgeLaplacian (left e) (right e)

/-- Uniform average of pair-averaging updates over a nonempty finite seam
family. -/
def uniformSeamRepair (left right : ε → ι) :
    (ι → ℝ) →ₗ[ℝ] (ι → ℝ) :=
  (1 / (Fintype.card ε : ℝ)) •
    ∑ e : ε, pairAverage (left e) (right e)

/-- Uniform seam averaging gives the identity minus the graph Laplacian
divided by twice the seam count. -/
theorem uniformSeamRepair_eq_id_sub_graphLaplacian
    [Nonempty ε] (left right : ε → ι) :
    uniformSeamRepair left right =
      LinearMap.id -
        (1 / (2 * (Fintype.card ε : ℝ))) • graphLaplacian left right := by
  rw [uniformSeamRepair, graphLaplacian]
  simp_rw [pairAverage_eq_id_sub_half_edgeLaplacian]
  rw [Finset.sum_sub_distrib]
  simp only [Finset.sum_const, Finset.card_univ]
  have hcard : (Fintype.card ε : ℝ) ≠ 0 := by
    exact_mod_cast Fintype.card_ne_zero
  have hnormalize :
      (1 / (Fintype.card ε : ℝ)) * (Fintype.card ε : ℝ) = 1 := by
    field_simp
  have hcoefficient :
      (1 / (Fintype.card ε : ℝ)) * (1 / 2 : ℝ) =
        1 / (2 * (Fintype.card ε : ℝ)) := by
    field_simp
  rw [← Nat.cast_smul_eq_nsmul ℝ]
  simp only [smul_sub, smul_smul, Finset.smul_sum]
  rw [hnormalize, one_smul]
  simp_rw [hcoefficient]

/-- For a family indexed by thirty seams, uniform pair averaging is exactly
`I - L / 60`.  The endpoint maps remain explicit inputs; this theorem does
not identify them with a physical carrier. -/
theorem uniform_thirty_seam_repair
    (left right : Fin 30 → ι) :
    uniformSeamRepair left right =
      LinearMap.id - (1 / 60 : ℝ) • graphLaplacian left right := by
  rw [uniformSeamRepair_eq_id_sub_graphLaplacian]
  norm_num

end

end ObserverPatchHolography.ScalarSeamRepair

/- Axiom audit: finite real linear algebra with standard axioms only. -/

#print axioms ObserverPatchHolography.ScalarSeamRepair.eq_pairAverage_of_supported_agreeing_sum_preserving
#print axioms ObserverPatchHolography.ScalarSeamRepair.uniform_thirty_seam_repair
