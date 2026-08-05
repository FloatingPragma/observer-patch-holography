import Geometry.CelestialNullCone

/-!
# A finite-effect closure boundary on the celestial sphere

This issue-scoped module records an exact obstruction to inferring a qubit
Born functional from binary-projector normalization, range, and continuity
alone.  On the unit sphere, the nonlinear frame function

`F(n) = (1 + n_z^3) / 2`

is continuous, lies in `[0,1]`, and satisfies `F(-n) = 1 - F(n)`, but it is
not of the affine/Born form `(1 + q · n) / 2`.  The nonrepresentation proof
uses only the three coordinate axes and the unit vector `(3/5, 0, 4/5)`.

The final result is a positive closure statement in the other direction:
once affinity is assumed, checking the centered functional on any dense set
of unit directions extends the bound to the whole sphere and forces the
coefficient vector into the closed unit ball.

Boundary: none of these statements derives affinity, a Gleason/Busch theorem,
an effect algebra, or source-produced physical effects.  The sphere and dot
product are the explicitly declared finite-dimensional C1 coordinates.
-/

namespace EventAlgebra.FiniteEffectClosureBoundary

open OPH.C1Lorentz

noncomputable section

/-- The exact nonlinear binary-projector frame function, defined on all
spatial coordinates so that continuity can be established before restriction
to the celestial unit sphere. -/
def nonlinearBinaryWeight (n : Spatial) : ℝ :=
  (1 + (n 2) ^ 3) / 2

/-- The centered affine/Born ansatz associated with a coefficient vector
`q`.  This is only a coordinate definition; no density-matrix semantics are
asserted here. -/
def affineBinaryWeight (q n : Spatial) : ℝ :=
  (1 + spatialDot q n) / 2

/-- Every coordinate of a unit spatial vector has square at most one. -/
theorem coordinate_sq_le_one {n : Spatial} (hn : spatialNormSq n = 1)
    (i : Fin 3) : (n i) ^ 2 ≤ 1 := by
  calc
    (n i) ^ 2 ≤ spatialNormSq n := by
      unfold spatialNormSq
      exact Finset.single_le_sum
        (fun j _ ↦ sq_nonneg (n j)) (Finset.mem_univ i)
    _ = 1 := hn

/-- Coordinate bounds inherited from the celestial unit-sphere equation. -/
theorem coordinate_mem_unitInterval {n : Spatial} (hn : spatialNormSq n = 1)
    (i : Fin 3) : -1 ≤ n i ∧ n i ≤ 1 := by
  have hsq := coordinate_sq_le_one hn i
  constructor <;> nlinarith [sq_nonneg (n i - 1), sq_nonneg (n i + 1)]

/-- The nonlinear frame function is globally continuous, hence its
restriction to the celestial sphere is continuous. -/
theorem continuous_nonlinearBinaryWeight : Continuous nonlinearBinaryWeight := by
  unfold nonlinearBinaryWeight
  fun_prop

/-- The explicit restriction of the nonlinear frame function to the C1
celestial sphere. -/
def celestialBinaryWeight (n : CelestialSphere) : ℝ :=
  nonlinearBinaryWeight n.1

theorem continuous_celestialBinaryWeight : Continuous celestialBinaryWeight := by
  exact continuous_nonlinearBinaryWeight.comp continuous_subtype_val

/-- On every celestial direction, the nonlinear frame function is a genuine
number in the probability interval `[0,1]`. -/
theorem nonlinearBinaryWeight_mem_Icc {n : Spatial}
    (hn : spatialNormSq n = 1) : nonlinearBinaryWeight n ∈ Set.Icc (0 : ℝ) 1 := by
  have hz := coordinate_mem_unitInterval hn (2 : Fin 3)
  have hquadLower : 0 ≤ (n 2) ^ 2 - n 2 + 1 := by
    nlinarith [sq_nonneg (n 2 - (1 / 2 : ℝ))]
  have hquadUpper : 0 ≤ (n 2) ^ 2 + n 2 + 1 := by
    nlinarith [sq_nonneg (n 2 + (1 / 2 : ℝ))]
  have hcubedLower : 0 ≤ (n 2 + 1) * ((n 2) ^ 2 - n 2 + 1) :=
    mul_nonneg (by linarith [hz.1]) hquadLower
  have hcubedUpper : 0 ≤ (1 - n 2) * ((n 2) ^ 2 + n 2 + 1) :=
    mul_nonneg (by linarith [hz.2]) hquadUpper
  constructor <;>
    simp only [nonlinearBinaryWeight] <;>
    nlinarith [hcubedLower, hcubedUpper]

/-- Exact normalization of antipodal binary contexts. -/
theorem nonlinearBinaryWeight_neg (n : Spatial) :
    nonlinearBinaryWeight (-n) = 1 - nonlinearBinaryWeight n := by
  simp only [nonlinearBinaryWeight, Pi.neg_apply]
  ring

theorem nonlinearBinaryWeight_antipodal_sum (n : Spatial) :
    nonlinearBinaryWeight n + nonlinearBinaryWeight (-n) = 1 := by
  rw [nonlinearBinaryWeight_neg]
  ring

/-- The nonlinear frame function has no affine/Born coefficient vector.
The witness is exact: the three positive coordinate axes force
`q = (0,0,1)`, while `(3/5,0,4/5)` distinguishes `z^3` from `z`. -/
theorem nonlinearBinaryWeight_not_affine :
    ¬ ∃ q : Spatial, ∀ n : Spatial, spatialNormSq n = 1 →
      nonlinearBinaryWeight n = affineBinaryWeight q n := by
  rintro ⟨q, hq⟩
  have hx := hq (![(1 : ℝ), 0, 0] : Spatial) (by
    norm_num [spatialNormSq, Fin.sum_univ_succ])
  have hy := hq (![(0 : ℝ), 1, 0] : Spatial) (by
    norm_num [spatialNormSq, Fin.sum_univ_succ])
  have hz := hq (![(0 : ℝ), 0, 1] : Spatial) (by
    norm_num [spatialNormSq, Fin.sum_univ_succ])
  have hq0 : q 0 = 0 := by
    simp [nonlinearBinaryWeight, affineBinaryWeight, spatialDot,
      Fin.sum_univ_succ, Matrix.cons_val_two] at hx
    norm_num at hx ⊢
    exact hx
  have hq1 : q 1 = 0 := by
    simp [nonlinearBinaryWeight, affineBinaryWeight, spatialDot,
      Fin.sum_univ_succ, Matrix.cons_val_two] at hy
    norm_num at hy ⊢
    exact hy
  have hq2 : q 2 = 1 := by
    simp [nonlinearBinaryWeight, affineBinaryWeight, spatialDot,
      Fin.sum_univ_succ, Matrix.cons_val_two] at hz
    linarith
  have hw := hq (![(3 / 5 : ℝ), 0, (4 / 5 : ℝ)] : Spatial) (by
    norm_num [spatialNormSq, Fin.sum_univ_succ])
  simp [nonlinearBinaryWeight, affineBinaryWeight, spatialDot,
    Fin.sum_univ_succ, Matrix.cons_val_two] at hw
  simp [hq0, hq2] at hw
  norm_num at hw

/-- Continuity of a fixed centered affine functional on the celestial
sphere. -/
theorem continuous_celestialSpatialDot (q : Spatial) :
    Continuous (fun n : CelestialSphere ↦ spatialDot q n.1) := by
  have h0 : Continuous (fun n : CelestialSphere ↦ n.1 (0 : Fin 3)) :=
    (continuous_apply (0 : Fin 3)).comp continuous_subtype_val
  have h1 : Continuous (fun n : CelestialSphere ↦ n.1 (1 : Fin 3)) :=
    (continuous_apply (1 : Fin 3)).comp continuous_subtype_val
  have h2 : Continuous (fun n : CelestialSphere ↦ n.1 (2 : Fin 3)) :=
    (continuous_apply (2 : Fin 3)).comp continuous_subtype_val
  simpa [spatialDot, Fin.sum_univ_succ] using
    (continuous_const.mul h0).add
      ((continuous_const.mul h1).add (continuous_const.mul h2))

/-- A centered affine bound verified on a dense family of celestial
directions extends to every direction.  Density is used only after the
affine functional has been supplied explicitly. -/
theorem dense_unit_tests_extend_affine_bound
    (D : Set CelestialSphere) (hDense : Dense D) (q : Spatial)
    (hTest : ∀ n ∈ D, |spatialDot q n.1| ≤ 1) :
    ∀ n : CelestialSphere, |spatialDot q n.1| ≤ 1 := by
  intro n
  have hClosed :
      IsClosed {m : CelestialSphere | |spatialDot q m.1| ≤ 1} :=
    isClosed_le (continuous_celestialSpatialDot q).abs continuous_const
  have hSubset : D ⊆ {m : CelestialSphere | |spatialDot q m.1| ≤ 1} := by
    intro m hm
    exact hTest m hm
  exact closure_minimal hSubset hClosed (hDense n)

/-- Spatial dot product with a scaled right argument. -/
theorem spatialDot_smul_right (q r : Spatial) (a : ℝ) :
    spatialDot q (a • r) = a * spatialDot q r := by
  simp only [spatialDot, Pi.smul_apply, smul_eq_mul]
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  ring

@[simp] theorem spatialDot_self (q : Spatial) :
    spatialDot q q = spatialNormSq q := by
  simp only [spatialDot, spatialNormSq, pow_two]

/-- Testing a centered affine functional on every unit direction forces its
coefficient vector into the closed Euclidean unit ball. -/
theorem all_unit_tests_force_closed_unit_ball (q : Spatial)
    (hTest : ∀ n : CelestialSphere, |spatialDot q n.1| ≤ 1) :
    spatialNormSq q ≤ 1 := by
  by_contra hNot
  have hNorm : 1 < spatialNormSq q := lt_of_not_ge hNot
  have hNormPos : 0 < spatialNormSq q := lt_trans zero_lt_one hNorm
  have hSqrtPos : 0 < Real.sqrt (spatialNormSq q) :=
    Real.sqrt_pos.2 hNormPos
  let a : ℝ := 1 / Real.sqrt (spatialNormSq q)
  let n : Spatial := a • q
  have hn : spatialNormSq n = 1 := by
    dsimp [n, a]
    rw [spatialNormSq_smul, div_pow, one_pow,
      Real.sq_sqrt hNormPos.le]
    field_simp [ne_of_gt hNormPos]
  let ns : CelestialSphere := ⟨n, hn⟩
  have hBound := hTest ns
  have hDot : spatialDot q n = Real.sqrt (spatialNormSq q) := by
    rw [show n = a • q by rfl, spatialDot_smul_right, spatialDot_self]
    dsimp [a]
    field_simp [ne_of_gt hSqrtPos]
    nlinarith [Real.sq_sqrt hNormPos.le]
  change |spatialDot q n| ≤ 1 at hBound
  rw [hDot, abs_of_nonneg (Real.sqrt_nonneg _)] at hBound
  nlinarith [Real.sq_sqrt hNormPos.le]

/-- Dense centered tests force the same closed-unit-ball constraint once
affinity has already been assumed. -/
theorem dense_unit_tests_force_closed_unit_ball
    (D : Set CelestialSphere) (hDense : Dense D) (q : Spatial)
    (hTest : ∀ n ∈ D, |spatialDot q n.1| ≤ 1) :
    spatialNormSq q ≤ 1 := by
  apply all_unit_tests_force_closed_unit_ball q
  exact dense_unit_tests_extend_affine_bound D hDense q hTest

/-- Probability-interval tests for the affine/Born ansatz on a dense family
of directions force its coefficient vector into the closed unit ball.  This
is a positivity-after-affinity result, not a derivation of affinity. -/
theorem dense_affine_probability_tests_force_closed_unit_ball
    (D : Set CelestialSphere) (hDense : Dense D) (q : Spatial)
    (hTest : ∀ n ∈ D, affineBinaryWeight q n.1 ∈ Set.Icc (0 : ℝ) 1) :
    spatialNormSq q ≤ 1 := by
  apply dense_unit_tests_force_closed_unit_ball D hDense q
  intro n hn
  have hInterval := hTest n hn
  rw [abs_le]
  simp only [affineBinaryWeight] at hInterval
  constructor <;> linarith [hInterval.1, hInterval.2]

end

end EventAlgebra.FiniteEffectClosureBoundary
