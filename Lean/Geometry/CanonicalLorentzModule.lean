import Mathlib

/-!
# C1: the canonical Hermitian Lorentz module

This file formalizes only the finite-dimensional geometry in completion lane
`C1`.  The Pauli-coordinate model `Herm2` gives unique real coordinates for
the vector space of `2 × 2` complex Hermitian matrices.  Matrix determinant is
the Lorentz quadratic form with one positive and three negative squares.  The
explicit linear equivalence to the Einstein tensor coordinates is stated in
`Geometry/EinsteinTensorBridge.lean`.

Nothing here identifies this module with physical spacetime, supplies clocks
or rods, takes a continuum limit, or derives the construction from observer
repair data.
-/

namespace OPH.C1Lorentz

noncomputable section

/-- The three real Pauli coordinates. -/
abbrev Spatial := Fin 3 → ℝ

/-- Pauli coordinates for `2 × 2` complex Hermitian matrices: one scalar
coordinate and three traceless coordinates. -/
abbrev Herm2 := ℝ × Spatial

/-- Squared Euclidean length on the spatial Pauli coordinates. -/
def spatialNormSq (x : Spatial) : ℝ := ∑ i, (x i) ^ 2

/-- The polarization of the spatial squared norm. -/
def spatialDot (x y : Spatial) : ℝ := ∑ i, x i * y i

/-- The determinant quadratic form in Pauli coordinates, with convention
`(+---)`. -/
def lorentzQ (v : Herm2) : ℝ := v.1 ^ 2 - spatialNormSq v.2

/-- The symmetric bilinear form polarizing `lorentzQ`. -/
def lorentzB (u v : Herm2) : ℝ := u.1 * v.1 - spatialDot u.2 v.2

@[simp] theorem lorentzB_self (v : Herm2) : lorentzB v v = lorentzQ v := by
  simp only [lorentzB, lorentzQ, spatialDot, spatialNormSq, pow_two]

theorem lorentzB_symm (u v : Herm2) : lorentzB u v = lorentzB v u := by
  simp only [lorentzB, spatialDot]
  congr 1
  · ring
  apply Finset.sum_congr rfl
  intro i _
  ring

theorem lorentzB_add_right (u v w : Herm2) :
    lorentzB u (v + w) = lorentzB u v + lorentzB u w := by
  simp only [lorentzB, spatialDot, Prod.fst_add, Prod.snd_add, Pi.add_apply]
  have hsum : (∑ i, u.2 i * (v.2 i + w.2 i)) =
      (∑ i, u.2 i * v.2 i) + ∑ i, u.2 i * w.2 i := by
    rw [← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro i _
    ring
  rw [hsum]
  ring

theorem lorentzB_smul_right (a : ℝ) (u v : Herm2) :
    lorentzB u (a • v) = a * lorentzB u v := by
  simp only [lorentzB, spatialDot, Prod.smul_fst, Prod.smul_snd,
    Pi.smul_apply, smul_eq_mul]
  have hsum : (∑ i, u.2 i * (a * v.2 i)) = a * ∑ i, u.2 i * v.2 i := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    ring
  rw [hsum]
  ring

theorem lorentzB_add_left (u v w : Herm2) :
    lorentzB (u + v) w = lorentzB u w + lorentzB v w := by
  rw [lorentzB_symm, lorentzB_add_right, lorentzB_symm w u,
    lorentzB_symm w v]

theorem lorentzB_smul_left (a : ℝ) (u v : Herm2) :
    lorentzB (a • u) v = a * lorentzB u v := by
  rw [lorentzB_symm, lorentzB_smul_right, lorentzB_symm v u]

/-- The standard Pauli-matrix realization
`t I + x σ₁ + y σ₂ + z σ₃`. -/
def toMatrix (v : Herm2) : Matrix (Fin 2) (Fin 2) ℂ :=
  !![(v.1 + v.2 2 : ℂ), (v.2 0 : ℂ) - (v.2 1 : ℂ) * Complex.I;
     (v.2 0 : ℂ) + (v.2 1 : ℂ) * Complex.I, (v.1 - v.2 2 : ℂ)]

theorem toMatrix_isHermitian (v : Herm2) : (toMatrix v).IsHermitian := by
  apply Matrix.IsHermitian.ext
  intro i j
  fin_cases i <;> fin_cases j <;>
    apply Complex.ext <;> simp [toMatrix]

/-- Determinant of the Pauli realization is exactly the Lorentz quadratic
form, viewed as a real complex number. -/
theorem det_toMatrix (v : Herm2) : Matrix.det (toMatrix v) = (lorentzQ v : ℂ) := by
  rw [Matrix.det_fin_two]
  simp [toMatrix, lorentzQ, spatialNormSq, Fin.sum_univ_succ]
  ring_nf
  rw [Complex.I_sq]
  ring

/-- The Pauli realization is additive. -/
theorem toMatrix_add (u v : Herm2) : toMatrix (u + v) = toMatrix u + toMatrix v := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [toMatrix] <;> ring

/-- The Pauli realization respects real scalar multiplication. -/
theorem toMatrix_smul (a : ℝ) (v : Herm2) :
    toMatrix (a • v) = a • toMatrix v := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [toMatrix] <;> ring

/-- The Pauli realization is injective. -/
theorem toMatrix_injective : Function.Injective toMatrix := by
  intro u v h
  have h00 := congr_fun (congr_fun h 0) 0
  have h11 := congr_fun (congr_fun h 1) 1
  have h01 := congr_fun (congr_fun h 0) 1
  apply Prod.ext
  · have h00re := congrArg Complex.re h00
    have h11re := congrArg Complex.re h11
    simp [toMatrix] at h00re h11re
    linarith
  · funext i
    fin_cases i
    · have := congrArg Complex.re h01
      simpa [toMatrix] using this
    · have := congrArg Complex.im h01
      simpa [toMatrix] using this
    · have h00re := congrArg Complex.re h00
      have h11re := congrArg Complex.re h11
      simp [toMatrix] at h00re h11re
      change u.2 (2 : Fin 3) = v.2 (2 : Fin 3)
      linarith

/-- Explicit coordinates extracted from an arbitrary complex `2 × 2`
matrix.  On Hermitian matrices these are inverse to `toMatrix`. -/
def matrixCoords (A : Matrix (Fin 2) (Fin 2) ℂ) : Herm2 :=
  Prod.mk (((A 0 0).re + (A 1 1).re) / 2)
    ![(A 0 1).re, -(A 0 1).im, ((A 0 0).re - (A 1 1).re) / 2]

theorem toMatrix_matrixCoords_of_isHermitian
    (A : Matrix (Fin 2) (Fin 2) ℂ) (hA : A.IsHermitian) :
    toMatrix (matrixCoords A) = A := by
  have h00 := hA.apply 0 0
  have h11 := hA.apply 1 1
  have h10 := hA.apply 1 0
  have h00im := congrArg Complex.im h00
  have h11im := congrArg Complex.im h11
  simp at h00im h11im
  have h00im0 : (A 0 0).im = 0 := by linarith [h00im]
  have h11im0 : (A 1 1).im = 0 := by linarith [h11im]
  ext i j
  fin_cases i <;> fin_cases j
  · apply Complex.ext
    · simp [matrixCoords, toMatrix]
      ring
    · simpa [matrixCoords, toMatrix] using h00im0.symm
  · apply Complex.ext <;> simp [matrixCoords, toMatrix]
  · apply Complex.ext
    · simpa [matrixCoords, toMatrix] using congrArg Complex.re h10
    · simpa [matrixCoords, toMatrix] using congrArg Complex.im h10
  · apply Complex.ext
    · simp [matrixCoords, toMatrix]
      ring
    · simpa [matrixCoords, toMatrix] using h11im0.symm

theorem matrixCoords_toMatrix (v : Herm2) : matrixCoords (toMatrix v) = v := by
  exact toMatrix_injective (toMatrix_matrixCoords_of_isHermitian
    (toMatrix v) (toMatrix_isHermitian v))

/-- Exact representation theorem: a complex `2 × 2` matrix is Hermitian iff
it has a unique real Pauli-coordinate representative. -/
theorem isHermitian_iff_existsUnique_toMatrix (A : Matrix (Fin 2) (Fin 2) ℂ) :
    A.IsHermitian ↔ ∃! v : Herm2, toMatrix v = A := by
  constructor
  · intro hA
    refine ⟨matrixCoords A, toMatrix_matrixCoords_of_isHermitian A hA, ?_⟩
    intro v hv
    exact toMatrix_injective (hv.trans (toMatrix_matrixCoords_of_isHermitian A hA).symm)
  · rintro ⟨v, rfl, _⟩
    exact toMatrix_isHermitian v

/-- The real dimension of the Hermitian `2 × 2` module is four. -/
theorem finrank_Herm2 : Module.finrank ℝ Herm2 = 4 := by
  simp [Herm2, Spatial]

/-- Constructive inertia certificate: in the canonical splitting into the
scalar line and traceless three-space, the determinant is one positive square
minus three squares. -/
theorem inertia_one_three (t : ℝ) (x : Spatial) :
    lorentzQ (t, x) = t ^ 2 - ∑ i : Fin 3, (x i) ^ 2 := rfl

theorem time_axis_positive {t : ℝ} (ht : t ≠ 0) : lorentzQ (t, 0) > 0 := by
  simp [lorentzQ, spatialNormSq, sq_pos_of_ne_zero ht]

theorem spatial_axis_negative {x : Spatial} (hx : x ≠ 0) : lorentzQ (0, x) < 0 := by
  have hnonneg : 0 ≤ spatialNormSq x := by
    exact Finset.sum_nonneg fun i _ ↦ sq_nonneg (x i)
  have hne : spatialNormSq x ≠ 0 := by
    intro hz
    apply hx
    funext i
    have hle : (x i) ^ 2 ≤ spatialNormSq x := by
      exact Finset.single_le_sum (fun j (_ : j ∈ Finset.univ) ↦ sq_nonneg (x j))
        (Finset.mem_univ i)
    have hi : (x i) ^ 2 = 0 :=
      le_antisymm (by simpa [hz] using hle) (sq_nonneg (x i))
    exact sq_eq_zero_iff.mp hi
  have hpos : 0 < spatialNormSq x := lt_of_le_of_ne hnonneg (Ne.symm hne)
  simpa [lorentzQ] using (neg_lt_zero.mpr hpos)

#print axioms det_toMatrix
#print axioms isHermitian_iff_existsUnique_toMatrix
#print axioms finrank_Herm2
#print axioms inertia_one_three
#print axioms time_axis_positive
#print axioms spatial_axis_negative

end

end OPH.C1Lorentz
