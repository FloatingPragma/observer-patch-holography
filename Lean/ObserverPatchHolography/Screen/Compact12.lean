import Mathlib

namespace OPH.Compact12

noncomputable section

abbrev I3 := Fin 3
abbrev CMat3 := Matrix I3 I3 ℂ
abbrev RMat3 := Matrix I3 I3 ℝ

/-- The compact real Lie algebra `u(3)`, realized as skew-Hermitian matrices. -/
def u3 : Submodule ℝ CMat3 where
  carrier := {A | A.conjTranspose = -A}
  zero_mem' := by simp
  add_mem' := fun {A B} hA hB ↦ by
    change A.conjTranspose = -A at hA
    change B.conjTranspose = -B at hB
    change (A + B).conjTranspose = -(A + B)
    rw [Matrix.conjTranspose_add, hA, hB]
    abel
  smul_mem' := fun r {A} hA ↦ by
    change A.conjTranspose = -A at hA
    change (r • A).conjTranspose = -(r • A)
    rw [Matrix.conjTranspose_smul, hA]
    rw [star_trivial]
    exact smul_neg r A

/-- The compact real Lie algebra `so(3)`, realized as skew-symmetric matrices. -/
def so3 : Submodule ℝ RMat3 where
  carrier := {A | A.transpose = -A}
  zero_mem' := by simp
  add_mem' := fun {A B} hA hB ↦ by
    change A.transpose = -A at hA
    change B.transpose = -B at hB
    change (A + B).transpose = -(A + B)
    rw [Matrix.transpose_add, hA, hB]
    abel
  smul_mem' := fun r {A} hA ↦ by
    change A.transpose = -A at hA
    change (r • A).transpose = -(r • A)
    rw [Matrix.transpose_smul, hA]
    simp

abbrev U3 := u3
abbrev SO3 := so3

def commC (A B : CMat3) : CMat3 := A * B - B * A

def commR (A B : RMat3) : RMat3 := A * B - B * A

theorem commC_antisymm (A B : CMat3) : commC A B = -commC B A := by
  simp [commC]

theorem commR_antisymm (A B : RMat3) : commR A B = -commR B A := by
  simp [commR]

theorem commC_jacobi (A B C : CMat3) :
    commC A (commC B C) + commC B (commC C A) + commC C (commC A B) = 0 := by
  simp [commC]
  noncomm_ring

theorem commR_jacobi (A B C : RMat3) :
    commR A (commR B C) + commR B (commR C A) + commR C (commR A B) = 0 := by
  simp [commR]
  noncomm_ring

theorem commC_add_left (A B C : CMat3) :
    commC (A + B) C = commC A C + commC B C := by
  unfold commC
  noncomm_ring

theorem commR_add_left (A B C : RMat3) :
    commR (A + B) C = commR A C + commR B C := by
  unfold commR
  noncomm_ring

theorem commC_add_right (A B C : CMat3) :
    commC A (B + C) = commC A B + commC A C := by
  unfold commC
  noncomm_ring

theorem commR_add_right (A B C : RMat3) :
    commR A (B + C) = commR A B + commR A C := by
  unfold commR
  noncomm_ring

theorem commC_self (A : CMat3) : commC A A = 0 := by simp [commC]

theorem commR_self (A : RMat3) : commR A A = 0 := by simp [commR]

theorem commC_leibniz (A B C : CMat3) :
    commC A (commC B C) = commC (commC A B) C + commC B (commC A C) := by
  unfold commC
  noncomm_ring

theorem commR_leibniz (A B C : RMat3) :
    commR A (commR B C) = commR (commR A B) C + commR B (commR A C) := by
  unfold commR
  noncomm_ring

theorem commC_smul_right (r : ℝ) (A B : CMat3) :
    commC A (r • B) = r • commC A B := by
  simp only [commC]
  rw [Algebra.mul_smul_comm, Algebra.smul_mul_assoc]
  exact (smul_sub r (A * B) (B * A)).symm

theorem commR_smul_right (r : ℝ) (A B : RMat3) :
    commR A (r • B) = r • commR A B := by
  simp only [commR]
  rw [Algebra.mul_smul_comm, Algebra.smul_mul_assoc]
  exact (smul_sub r (A * B) (B * A)).symm

def uBracket (A B : U3) : U3 := ⟨commC A.1 B.1, by
  change (A.1 * B.1 - B.1 * A.1).conjTranspose = -(A.1 * B.1 - B.1 * A.1)
  rw [Matrix.conjTranspose_sub, Matrix.conjTranspose_mul,
    Matrix.conjTranspose_mul, A.2, B.2]
  noncomm_ring⟩

def soBracket (A B : SO3) : SO3 := ⟨commR A.1 B.1, by
  change (A.1 * B.1 - B.1 * A.1).transpose = -(A.1 * B.1 - B.1 * A.1)
  rw [Matrix.transpose_sub, Matrix.transpose_mul,
    Matrix.transpose_mul, A.2, B.2]
  noncomm_ring⟩

/-- The direct-sum compact algebra `u(3) ⊕ so(3)`. -/
abbrev Compact := U3 × SO3

def bracket (X Y : Compact) : Compact :=
  (uBracket X.1 Y.1, soBracket X.2 Y.2)

instance compactBracket : Bracket Compact Compact := ⟨bracket⟩

theorem bracket_antisymm (X Y : Compact) : bracket X Y = -bracket Y X := by
  apply Prod.ext
  · apply Subtype.ext
    change commC X.1.1 Y.1.1 = -commC Y.1.1 X.1.1
    exact commC_antisymm _ _
  · apply Subtype.ext
    change commR X.2.1 Y.2.1 = -commR Y.2.1 X.2.1
    exact commR_antisymm _ _

theorem bracket_jacobi (X Y Z : Compact) :
    bracket X (bracket Y Z) + bracket Y (bracket Z X) + bracket Z (bracket X Y) = 0 := by
  apply Prod.ext
  · apply Subtype.ext
    change commC X.1.1 (commC Y.1.1 Z.1.1) +
      commC Y.1.1 (commC Z.1.1 X.1.1) + commC Z.1.1 (commC X.1.1 Y.1.1) = 0
    exact commC_jacobi _ _ _
  · apply Subtype.ext
    change commR X.2.1 (commR Y.2.1 Z.2.1) +
      commR Y.2.1 (commR Z.2.1 X.2.1) + commR Z.2.1 (commR X.2.1 Y.2.1) = 0
    exact commR_jacobi _ _ _

instance compactLieRing : LieRing Compact where
  add_lie X Y Z := by
    apply Prod.ext
    · apply Subtype.ext
      change commC (X.1.1 + Y.1.1) Z.1.1 =
        commC X.1.1 Z.1.1 + commC Y.1.1 Z.1.1
      exact commC_add_left _ _ _
    · apply Subtype.ext
      change commR (X.2.1 + Y.2.1) Z.2.1 =
        commR X.2.1 Z.2.1 + commR Y.2.1 Z.2.1
      exact commR_add_left _ _ _
  lie_add X Y Z := by
    apply Prod.ext
    · apply Subtype.ext
      change commC X.1.1 (Y.1.1 + Z.1.1) =
        commC X.1.1 Y.1.1 + commC X.1.1 Z.1.1
      exact commC_add_right _ _ _
    · apply Subtype.ext
      change commR X.2.1 (Y.2.1 + Z.2.1) =
        commR X.2.1 Y.2.1 + commR X.2.1 Z.2.1
      exact commR_add_right _ _ _
  lie_self X := by
    apply Prod.ext
    · apply Subtype.ext
      change commC X.1.1 X.1.1 = 0
      exact commC_self _
    · apply Subtype.ext
      change commR X.2.1 X.2.1 = 0
      exact commR_self _
  leibniz_lie X Y Z := by
    apply Prod.ext
    · apply Subtype.ext
      change commC X.1.1 (commC Y.1.1 Z.1.1) =
        commC (commC X.1.1 Y.1.1) Z.1.1 + commC Y.1.1 (commC X.1.1 Z.1.1)
      exact commC_leibniz _ _ _
    · apply Subtype.ext
      change commR X.2.1 (commR Y.2.1 Z.2.1) =
        commR (commR X.2.1 Y.2.1) Z.2.1 + commR Y.2.1 (commR X.2.1 Z.2.1)
      exact commR_leibniz _ _ _

instance compactLieAlgebra : LieAlgebra ℝ Compact where
  lie_smul r X Y := by
    apply Prod.ext
    · apply Subtype.ext
      change commC X.1.1 (r • Y.1.1) = r • commC X.1.1 Y.1.1
      exact commC_smul_right _ _ _
    · apply Subtype.ext
      change commR X.2.1 (r • Y.2.1) = r • commR X.2.1 Y.2.1
      exact commR_smul_right _ _ _

/-! Explicit coordinates certify the real dimension `9 + 3 = 12`. -/

abbrev UCoord := (I3 → ℝ) × (I3 → ℂ)
abbrev SOCoord := I3 → ℝ

def uCoord (A : U3) : UCoord :=
  (![(A.1 0 0).im, (A.1 1 1).im, (A.1 2 2).im],
   ![A.1 0 1, A.1 0 2, A.1 1 2])

def uMatrix (p : UCoord) : CMat3 :=
  !![Complex.I * p.1 0, p.2 0, p.2 1;
     -star (p.2 0), Complex.I * p.1 1, p.2 2;
     -star (p.2 1), -star (p.2 2), Complex.I * p.1 2]

theorem uMatrix_mem (p : UCoord) : uMatrix p ∈ u3 := by
  change (uMatrix p).conjTranspose = -uMatrix p
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [uMatrix, Matrix.conjTranspose_apply]

def soCoord (A : SO3) : SOCoord := ![A.1 0 1, A.1 0 2, A.1 1 2]

def soMatrix (p : SOCoord) : RMat3 :=
  !![0, p 0, p 1; -p 0, 0, p 2; -p 1, -p 2, 0]

theorem soMatrix_mem (p : SOCoord) : soMatrix p ∈ so3 := by
  change (soMatrix p).transpose = -soMatrix p
  ext i j
  fin_cases i <;> fin_cases j <;> simp [soMatrix, Matrix.transpose_apply]

theorem eq_I_mul_im_of_star_eq_neg {z : ℂ} (h : star z = -z) :
    Complex.I * (z.im : ℂ) = z := by
  apply Complex.ext
  · have hr := congrArg Complex.re h
    simp at hr ⊢
    linarith
  · simp

theorem uMatrix_uCoord (A : U3) :
    uMatrix (uCoord A) = A.1 := by
  have h := A.2
  change A.1.conjTranspose = -A.1 at h
  have h00 : star (A.1 0 0) = -A.1 0 0 := by
    simpa [Matrix.conjTranspose_apply] using congrFun (congrFun h 0) 0
  have h11 : star (A.1 1 1) = -A.1 1 1 := by
    simpa [Matrix.conjTranspose_apply] using congrFun (congrFun h 1) 1
  have h22 : star (A.1 2 2) = -A.1 2 2 := by
    simpa [Matrix.conjTranspose_apply] using congrFun (congrFun h 2) 2
  have h10 : star (A.1 0 1) = -A.1 1 0 := by
    simpa [Matrix.conjTranspose_apply] using congrFun (congrFun h 1) 0
  have h20 : star (A.1 0 2) = -A.1 2 0 := by
    simpa [Matrix.conjTranspose_apply] using congrFun (congrFun h 2) 0
  have h21 : star (A.1 1 2) = -A.1 2 1 := by
    simpa [Matrix.conjTranspose_apply] using congrFun (congrFun h 2) 1
  ext i j
  fin_cases i <;> fin_cases j
  · simpa [uMatrix, uCoord] using eq_I_mul_im_of_star_eq_neg h00
  · simp [uMatrix, uCoord]
  · simp [uMatrix, uCoord]
  · simp [uMatrix, uCoord, h10]
  · simpa [uMatrix, uCoord] using eq_I_mul_im_of_star_eq_neg h11
  · simp [uMatrix, uCoord]
  · simp [uMatrix, uCoord, h20]
  · simp [uMatrix, uCoord, h21]
  · simpa [uMatrix, uCoord] using eq_I_mul_im_of_star_eq_neg h22

theorem uCoord_uMatrix (p : UCoord) :
    uCoord ⟨uMatrix p, uMatrix_mem p⟩ = p := by
  rcases p with ⟨d, z⟩
  apply Prod.ext
  · funext i
    fin_cases i <;> simp [uCoord, uMatrix]
  · funext i
    fin_cases i <;> simp [uCoord, uMatrix]

theorem soMatrix_soCoord (A : SO3) :
    soMatrix (soCoord A) = A.1 := by
  have h := A.2
  change A.1.transpose = -A.1 at h
  have h00 := congrFun (congrFun h 0) 0
  have h11 := congrFun (congrFun h 1) 1
  have h22 := congrFun (congrFun h 2) 2
  have h10 := congrFun (congrFun h 1) 0
  have h20 := congrFun (congrFun h 2) 0
  have h21 := congrFun (congrFun h 2) 1
  simp [Matrix.transpose_apply] at h00 h11 h22 h10 h20 h21
  ext i j
  fin_cases i <;> fin_cases j <;> simp [soMatrix, soCoord] <;>
    linarith

theorem soCoord_soMatrix (p : SOCoord) :
    soCoord ⟨soMatrix p, soMatrix_mem p⟩ = p := by
  funext i
  fin_cases i <;> simp [soCoord, soMatrix]

def uCoordLinear : U3 →ₗ[ℝ] UCoord where
  toFun := uCoord
  map_add' A B := by
    apply Prod.ext <;> funext i <;> fin_cases i <;> simp [uCoord]
  map_smul' r A := by
    apply Prod.ext
    · funext i
      fin_cases i
      · change (r • A.1 0 0).im = r * (A.1 0 0).im
        simpa only [smul_eq_mul] using (Complex.smul_im r (A.1 0 0))
      · change (r • A.1 1 1).im = r * (A.1 1 1).im
        simpa only [smul_eq_mul] using (Complex.smul_im r (A.1 1 1))
      · change (r • A.1 2 2).im = r * (A.1 2 2).im
        simpa only [smul_eq_mul] using (Complex.smul_im r (A.1 2 2))
    · funext i
      fin_cases i <;> simp [uCoord]

def soCoordLinear : SO3 →ₗ[ℝ] SOCoord where
  toFun := soCoord
  map_add' A B := by
    funext i
    fin_cases i <;> simp [soCoord]
  map_smul' r A := by
    funext i
    fin_cases i <;> simp [soCoord]

def uCoordEquiv : U3 ≃ₗ[ℝ] UCoord :=
  LinearEquiv.ofBijective uCoordLinear ⟨
    fun A B hab ↦ Subtype.ext <| by
      change uCoord A = uCoord B at hab
      rw [← uMatrix_uCoord A, ← uMatrix_uCoord B, hab],
    fun p ↦ ⟨⟨uMatrix p, uMatrix_mem p⟩, uCoord_uMatrix p⟩⟩

def soCoordEquiv : SO3 ≃ₗ[ℝ] SOCoord :=
  LinearEquiv.ofBijective soCoordLinear ⟨
    fun A B hab ↦ Subtype.ext <| by
      change soCoord A = soCoord B at hab
      rw [← soMatrix_soCoord A, ← soMatrix_soCoord B, hab],
    fun p ↦ ⟨⟨soMatrix p, soMatrix_mem p⟩, soCoord_soMatrix p⟩⟩

theorem finrank_u3 : Module.finrank ℝ U3 = 9 := by
  rw [uCoordEquiv.finrank_eq]
  rw [Module.finrank_prod, Module.finrank_pi, Module.finrank_pi_fintype]
  norm_num [Complex.finrank_real_complex]

theorem finrank_so3 : Module.finrank ℝ SO3 = 3 := by
  rw [soCoordEquiv.finrank_eq]
  simp [SOCoord]

theorem finrank_compact : Module.finrank ℝ Compact = 12 := by
  rw [Module.finrank_prod, finrank_u3, finrank_so3]

/-- A concrete symmetric, traceless quadrupole on the three coordinate axes. -/
def H : RMat3 := !![(1 : ℝ), 0, 0; 0, -1, 0; 0, 0, 0]

/-- A concrete infinitesimal rotation in the 0--1 plane. -/
def K : RMat3 := !![(0 : ℝ), 1, 0; -1, 0, 0; 0, 0, 0]

theorem H_symmetric : H.transpose = H := by
  ext i j
  fin_cases i <;> fin_cases j <;> norm_num [H, Matrix.transpose_apply]

theorem H_traceless : H.trace = 0 := by
  norm_num [Matrix.trace, H, Fin.sum_univ_succ]

theorem H_noncentral : H * K - K * H ≠ 0 := by
  intro h
  have h01 := congrFun (congrFun h (0 : I3)) (1 : I3)
  norm_num [H, K, Matrix.mul_apply, Fin.sum_univ_succ] at h01

/-- Multiplication by `i` turns the symmetric observable `H` into a `u(3)` generator. -/
def iH : U3 := ⟨fun i j ↦ Complex.I * H i j, by
  change Matrix.conjTranspose (fun i j ↦ Complex.I * H i j : CMat3) =
    -(fun i j ↦ Complex.I * H i j : CMat3)
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [H, Matrix.conjTranspose_apply]⟩

/-- The same real rotation matrix is skew-Hermitian, hence lies in `u(3)`. -/
def kU : U3 := ⟨fun i j ↦ (K i j : ℂ), by
  change Matrix.conjTranspose (fun i j ↦ (K i j : ℂ) : CMat3) =
    -(fun i j ↦ (K i j : ℂ) : CMat3)
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [K, Matrix.conjTranspose_apply]⟩

theorem iH_noncentral_in_u3 : uBracket iH kU ≠ 0 := by
  intro h
  have h01 := congrFun (congrFun (congrArg Subtype.val h) (0 : I3)) (1 : I3)
  norm_num [uBracket, commC, iH, kU, H, K, Matrix.mul_apply,
    Fin.sum_univ_succ] at h01

#print axioms commC_jacobi
#print axioms commR_jacobi
#print axioms bracket_antisymm
#print axioms bracket_jacobi
#print axioms compactLieRing
#print axioms compactLieAlgebra
#print axioms uCoordEquiv
#print axioms soCoordEquiv
#print axioms finrank_u3
#print axioms finrank_so3
#print axioms finrank_compact
#print axioms H_symmetric
#print axioms H_traceless
#print axioms H_noncentral
#print axioms iH_noncentral_in_u3

end

end OPH.Compact12
