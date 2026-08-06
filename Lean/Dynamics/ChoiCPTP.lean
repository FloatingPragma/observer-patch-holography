import Dynamics.ConditionalExpectationGenerator
import EventAlgebra.PartitionAverageCP
import EventAlgebra.TwoScalePublicRepair

/-!
# Complete positivity and the CPTP discharge for the B2 channels

The B2 packet supplies explicit Kraus data and trace preservation for the
partition-averaging and partition-pinching channels, and it records the
absence of a formal complete-positivity predicate as an explicit nonclaim.
This module defines that predicate and discharges the CPTP claims.

* `IsCompletelyPositive` is amplified positivity: for every ancilla size
  `m`, the block-entrywise amplification `amplify m Phi` (the map
  `id_m tensor Phi` written on `(Fin m x Fin n)`-indexed matrices) sends
  positive-semidefinite matrices to positive-semidefinite matrices.
* `isCompletelyPositive_of_kraus` proves the Kraus criterion: a map with a
  finite Kraus family `Phi X = sum_c K_c X K_c^H` is completely positive,
  by amplifying each Kraus operator to `1 tensor K_c` and conjugating.
* `IsCPTP` bundles complete positivity with trace preservation, with
  closure lemmas for the identity, composition, sums, nonnegative real
  scalings, and convex combinations.
* The **discharge**: `partitionAverage_isCPTP` and
  `partitionPinching_isCPTP` consume the committed Kraus identities
  (`partitionAverage_kraus_form`, `partitionPinching_kraus_form`) and the
  committed trace identities (`trace_partitionAverage`,
  `trace_partitionPinching`); nothing is redefined.
* `relaxationChannel` is the B2 relaxation-semigroup member written as the
  convex combination `exp(-gamma t) id + (1 - exp(-gamma t)) E` of the
  identity and the pinching; `relaxationChannel_apply_eq_publicRelaxTime`
  proves pointwise agreement with the committed `publicRelaxTime`, and
  `relaxationChannel_isCPTP` proves every member with `gamma, t >= 0` is
  CPTP.
* `choiMatrix` is the finite Choi matrix; complete positivity forces it to
  be positive semidefinite (`choiMatrix_posSemidef`).
* Negative control: the transpose map on `M_2` is positive and trace
  preserving and fails complete positivity
  (`transposeMap_not_completelyPositive`), refuted through an exact 4x4
  computation showing its Choi matrix has the negative expectation `-2` on
  an explicit vector.
-/

namespace OPH.Dynamics

open Matrix
open scoped ComplexOrder

variable {n k : ℕ}

/-! ## Blocks, amplification, and the complete-positivity predicate -/

/-- The `(a, b)` block of a `(Fin m x Fin n)`-indexed matrix. -/
def blockOf {m : ℕ} (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ)
    (a b : Fin m) : CMat n :=
  Matrix.of fun i j => X (a, i) (b, j)

@[simp]
theorem blockOf_apply {m : ℕ} (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ)
    (a b : Fin m) (i j : Fin n) :
    blockOf X a b i j = X (a, i) (b, j) :=
  rfl

/-- The `m`-fold matrix amplification `id_m tensor Phi`, acting blockwise on
`(Fin m x Fin n)`-indexed matrices. -/
def amplify (m : ℕ) (Φ : CMat n →ₗ[ℂ] CMat n) :
    Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ →ₗ[ℂ]
      Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ where
  toFun X := Matrix.of fun p q => Φ (blockOf X p.1 q.1) p.2 q.2
  map_add' X Y := by
    ext p q
    show Φ (blockOf (X + Y) p.1 q.1) p.2 q.2 =
      Φ (blockOf X p.1 q.1) p.2 q.2 + Φ (blockOf Y p.1 q.1) p.2 q.2
    rw [show blockOf (X + Y) p.1 q.1 =
      blockOf X p.1 q.1 + blockOf Y p.1 q.1 from rfl, map_add,
      Matrix.add_apply]
  map_smul' c X := by
    ext p q
    show Φ (blockOf (c • X) p.1 q.1) p.2 q.2 =
      c • Φ (blockOf X p.1 q.1) p.2 q.2
    rw [show blockOf (c • X) p.1 q.1 = c • blockOf X p.1 q.1 from rfl,
      map_smul, Matrix.smul_apply]

@[simp]
theorem amplify_apply {m : ℕ} (Φ : CMat n →ₗ[ℂ] CMat n)
    (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ) (p q : Fin m × Fin n) :
    amplify m Φ X p q = Φ (blockOf X p.1 q.1) p.2 q.2 :=
  rfl

/-- **Complete positivity**, defined as amplified positivity: every matrix
amplification of the map preserves positive semidefiniteness. -/
def IsCompletelyPositive (Φ : CMat n →ₗ[ℂ] CMat n) : Prop :=
  ∀ (m : ℕ) (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ),
    X.PosSemidef → (amplify m Φ X).PosSemidef

/-- Positivity of a matrix map, with no amplification. -/
def IsPositiveMap (Φ : CMat n →ₗ[ℂ] CMat n) : Prop :=
  ∀ X : CMat n, X.PosSemidef → (Φ X).PosSemidef

/-- Trace preservation of a matrix map. -/
def IsTracePreserving (Φ : CMat n →ₗ[ℂ] CMat n) : Prop :=
  ∀ X : CMat n, (Φ X).trace = X.trace

/-! ## The Kraus criterion -/

/-- The amplified Kraus operator `1_m tensor K`. -/
def oneKron (m : ℕ) (K : CMat n) :
    Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ :=
  Matrix.of fun p q => if p.1 = q.1 then K p.2 q.2 else 0

@[simp]
theorem oneKron_apply (m : ℕ) (K : CMat n) (p q : Fin m × Fin n) :
    oneKron m K p q = if p.1 = q.1 then K p.2 q.2 else 0 :=
  rfl

theorem oneKron_conjTranspose (m : ℕ) (K : CMat n) :
    (oneKron m K)ᴴ = oneKron m Kᴴ := by
  ext p q
  simp only [Matrix.conjTranspose_apply, oneKron_apply]
  by_cases h : p.1 = q.1
  · rw [if_pos h.symm, if_pos h]
  · rw [if_neg (h ∘ Eq.symm), star_zero, if_neg h]

/-- Left multiplication by `1_m tensor K` acts as `K` on each block row. -/
theorem oneKron_mul_apply (m : ℕ) (K : CMat n)
    (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ)
    (a : Fin m) (i : Fin n) (q : Fin m × Fin n) :
    (oneKron m K * X) (a, i) q = (K * blockOf X a q.1) i q.2 := by
  rw [Matrix.mul_apply, Matrix.mul_apply, Fintype.sum_prod_type]
  rw [Finset.sum_eq_single a
    (fun r1 _ hr1 => Finset.sum_eq_zero fun r2 _ => by
      simp [oneKron, Ne.symm hr1])
    (fun h => absurd (Finset.mem_univ a) h)]
  refine Finset.sum_congr rfl fun r2 _ => ?_
  simp [oneKron, blockOf]

/-- Right multiplication by `(1_m tensor K)^H` acts as `K^H` on each block
column. -/
theorem mul_oneKron_conjTranspose_apply (m : ℕ) (K : CMat n)
    (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ)
    (p : Fin m × Fin n) (b : Fin m) (j : Fin n) :
    (X * (oneKron m K)ᴴ) p (b, j) = (blockOf X p.1 b * Kᴴ) p.2 j := by
  rw [oneKron_conjTranspose, Matrix.mul_apply, Matrix.mul_apply,
    Fintype.sum_prod_type]
  rw [Finset.sum_eq_single b
    (fun s1 _ hs1 => Finset.sum_eq_zero fun s2 _ => by
      simp [oneKron, hs1])
    (fun h => absurd (Finset.mem_univ b) h)]
  refine Finset.sum_congr rfl fun s2 _ => ?_
  simp [oneKron, blockOf]

/-- Conjugation by an amplified Kraus operator is blockwise conjugation by
the original Kraus operator. -/
theorem oneKron_conj_apply (m : ℕ) (K : CMat n)
    (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ)
    (a b : Fin m) (i j : Fin n) :
    (oneKron m K * X * (oneKron m K)ᴴ) (a, i) (b, j) =
      (K * blockOf X a b * Kᴴ) i j := by
  have hblock : blockOf (X * (oneKron m K)ᴴ) a b = blockOf X a b * Kᴴ := by
    ext p q
    rw [blockOf_apply, mul_oneKron_conjTranspose_apply]
  rw [mul_assoc, oneKron_mul_apply, hblock, ← mul_assoc]

/-- A map with a finite Kraus family amplifies to a sum of conjugations by
the amplified Kraus operators. -/
theorem amplify_eq_sum_oneKron_conj {ι : Type*} [Fintype ι]
    (Φ : CMat n →ₗ[ℂ] CMat n) (K : ι → CMat n)
    (hK : ∀ X, Φ X = ∑ c, K c * X * (K c)ᴴ) (m : ℕ)
    (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ) :
    amplify m Φ X = ∑ c, oneKron m (K c) * X * (oneKron m (K c))ᴴ := by
  ext p q
  obtain ⟨a, i⟩ := p
  obtain ⟨b, j⟩ := q
  rw [amplify_apply, hK, Matrix.sum_apply, Matrix.sum_apply]
  exact Finset.sum_congr rfl fun c _ =>
    (oneKron_conj_apply m (K c) X a b i j).symm

/-- **Kraus criterion.** A map admitting a finite Kraus family is
completely positive: every amplification is a sum of conjugations, and
conjugation preserves positive semidefiniteness. -/
theorem isCompletelyPositive_of_kraus {ι : Type*} [Fintype ι]
    (Φ : CMat n →ₗ[ℂ] CMat n) (K : ι → CMat n)
    (hK : ∀ X, Φ X = ∑ c, K c * X * (K c)ᴴ) :
    IsCompletelyPositive Φ := by
  intro m X hX
  rw [amplify_eq_sum_oneKron_conj Φ K hK m X]
  exact Matrix.posSemidef_sum _ fun c _ =>
    hX.mul_mul_conjTranspose_same _

/-- Reindexing equivalence for the one-dimensional amplification. -/
def oneProdEquiv (n : ℕ) : Fin 1 × Fin n ≃ Fin n where
  toFun := Prod.snd
  invFun i := (0, i)
  left_inv p := by
    obtain ⟨a, i⟩ := p
    obtain rfl : a = 0 := Subsingleton.elim a 0
    rfl
  right_inv i := rfl

/-- A completely positive map is in particular positive, through the
trivial one-dimensional amplification. -/
theorem IsCompletelyPositive.isPositiveMap {Φ : CMat n →ₗ[ℂ] CMat n}
    (hΦ : IsCompletelyPositive Φ) : IsPositiveMap Φ := by
  intro X hX
  have hsub : (X.submatrix (Prod.snd : Fin 1 × Fin n → Fin n)
      Prod.snd).PosSemidef := hX.submatrix _
  have hamp := hΦ 1 _ hsub
  have hEq : amplify 1 Φ (X.submatrix Prod.snd Prod.snd) =
      (Φ X).submatrix (Prod.snd : Fin 1 × Fin n → Fin n) Prod.snd := by
    ext p q
    obtain ⟨a, i⟩ := p
    obtain ⟨b, j⟩ := q
    rfl
  rw [hEq] at hamp
  have hamp' : ((Φ X).submatrix ⇑(oneProdEquiv n) ⇑(oneProdEquiv n)).PosSemidef :=
    hamp
  exact (Matrix.posSemidef_submatrix_equiv (oneProdEquiv n)).mp hamp'

/-! ## CPTP and closure lemmas -/

/-- **CPTP**: completely positive and trace preserving. -/
structure IsCPTP (Φ : CMat n →ₗ[ℂ] CMat n) : Prop where
  completelyPositive : IsCompletelyPositive Φ
  tracePreserving : IsTracePreserving Φ

theorem amplify_id (m : ℕ)
    (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ) :
    amplify m (LinearMap.id : CMat n →ₗ[ℂ] CMat n) X = X := by
  ext p q
  obtain ⟨a, i⟩ := p
  obtain ⟨b, j⟩ := q
  rfl

theorem amplify_comp (m : ℕ) (Φ Ψ : CMat n →ₗ[ℂ] CMat n)
    (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ) :
    amplify m (Φ ∘ₗ Ψ) X = amplify m Φ (amplify m Ψ X) := by
  ext p q
  obtain ⟨a, i⟩ := p
  obtain ⟨b, j⟩ := q
  rfl

theorem amplify_add (m : ℕ) (Φ Ψ : CMat n →ₗ[ℂ] CMat n)
    (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ) :
    amplify m (Φ + Ψ) X = amplify m Φ X + amplify m Ψ X := by
  ext p q
  obtain ⟨a, i⟩ := p
  obtain ⟨b, j⟩ := q
  rfl

theorem amplify_smul (m : ℕ) (c : ℂ) (Φ : CMat n →ₗ[ℂ] CMat n)
    (X : Matrix (Fin m × Fin n) (Fin m × Fin n) ℂ) :
    amplify m (c • Φ) X = c • amplify m Φ X := by
  ext p q
  obtain ⟨a, i⟩ := p
  obtain ⟨b, j⟩ := q
  rfl

theorem isCompletelyPositive_id :
    IsCompletelyPositive (LinearMap.id : CMat n →ₗ[ℂ] CMat n) := by
  intro m X hX
  rwa [amplify_id]

theorem IsCompletelyPositive.comp {Φ Ψ : CMat n →ₗ[ℂ] CMat n}
    (hΦ : IsCompletelyPositive Φ) (hΨ : IsCompletelyPositive Ψ) :
    IsCompletelyPositive (Φ ∘ₗ Ψ) := by
  intro m X hX
  rw [amplify_comp]
  exact hΦ m _ (hΨ m X hX)

theorem IsCompletelyPositive.add {Φ Ψ : CMat n →ₗ[ℂ] CMat n}
    (hΦ : IsCompletelyPositive Φ) (hΨ : IsCompletelyPositive Ψ) :
    IsCompletelyPositive (Φ + Ψ) := by
  intro m X hX
  rw [amplify_add]
  exact (hΦ m X hX).add (hΨ m X hX)

theorem IsCompletelyPositive.smul_ofReal {Φ : CMat n →ₗ[ℂ] CMat n}
    (hΦ : IsCompletelyPositive Φ) {r : ℝ} (hr : 0 ≤ r) :
    IsCompletelyPositive ((r : ℂ) • Φ) := by
  intro m X hX
  rw [amplify_smul]
  exact (hΦ m X hX).smul (Complex.zero_le_real.mpr hr)

theorem isCPTP_id : IsCPTP (LinearMap.id : CMat n →ₗ[ℂ] CMat n) :=
  ⟨isCompletelyPositive_id, fun _ => rfl⟩

theorem IsCPTP.comp {Φ Ψ : CMat n →ₗ[ℂ] CMat n}
    (hΦ : IsCPTP Φ) (hΨ : IsCPTP Ψ) : IsCPTP (Φ ∘ₗ Ψ) :=
  ⟨hΦ.completelyPositive.comp hΨ.completelyPositive, fun X => by
    rw [LinearMap.comp_apply, hΦ.tracePreserving, hΨ.tracePreserving]⟩

/-- Convex combinations of CPTP maps with nonnegative real weights summing
to one are CPTP. -/
theorem IsCPTP.convexCombination {Φ Ψ : CMat n →ₗ[ℂ] CMat n}
    (hΦ : IsCPTP Φ) (hΨ : IsCPTP Ψ) {a b : ℝ}
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1) :
    IsCPTP ((a : ℂ) • Φ + (b : ℂ) • Ψ) := by
  refine ⟨(hΦ.completelyPositive.smul_ofReal ha).add
    (hΨ.completelyPositive.smul_ofReal hb), fun X => ?_⟩
  have hone : (a : ℂ) + (b : ℂ) = 1 := by
    rw [← Complex.ofReal_add, hab, Complex.ofReal_one]
  rw [LinearMap.add_apply, LinearMap.smul_apply, LinearMap.smul_apply,
    Matrix.trace_add, Matrix.trace_smul, Matrix.trace_smul,
    hΦ.tracePreserving, hΨ.tracePreserving, ← add_smul, hone, one_smul]

/-! ## The B2 discharge -/

/-- **CPTP discharge for partition averaging.** The committed Parseval-frame
Kraus family and the committed trace identity assemble into the formal CPTP
statement for the B2 publicization channel. -/
theorem partitionAverage_isCPTP (part : EventAlgebra.ProjectivePartition n k) :
    IsCPTP (EventAlgebra.partitionAverageLinearMap part) := by
  refine ⟨isCompletelyPositive_of_kraus _
    (fun c : Fin k × Fin n × Fin n =>
      EventAlgebra.partitionAverageKraus part c.1 c.2.1 c.2.2) ?_,
    fun X => EventAlgebra.trace_partitionAverage part X⟩
  intro X
  rw [EventAlgebra.partitionAverageLinearMap_apply,
    EventAlgebra.partitionAverage_kraus_form part X]
  simp only [Fintype.sum_prod_type]

/-- **CPTP discharge for partition pinching.** The partition projectors are
the committed Kraus family; trace preservation is the committed identity. -/
theorem partitionPinching_isCPTP (part : EventAlgebra.ProjectivePartition n k) :
    IsCPTP (EventAlgebra.partitionPinchingLinearMap part) := by
  refine ⟨isCompletelyPositive_of_kraus _ part.proj ?_,
    fun X => EventAlgebra.trace_partitionPinching part X⟩
  intro X
  rw [EventAlgebra.partitionPinchingLinearMap_apply]
  exact EventAlgebra.partitionPinching_kraus_form part X

/-! ## The relaxation semigroup members are CPTP -/

/-- The B2 relaxation-semigroup member at rate `gamma` and time `t`,
written as the convex combination of the identity and the pinching with
weights `exp(-gamma t)` and `1 - exp(-gamma t)`. -/
noncomputable def relaxationChannel (gamma t : ℝ)
    (part : EventAlgebra.ProjectivePartition n k) : CMat n →ₗ[ℂ] CMat n :=
  ((Real.exp (-gamma * t) : ℝ) : ℂ) • LinearMap.id +
    ((1 - Real.exp (-gamma * t) : ℝ) : ℂ) •
      EventAlgebra.partitionPinchingLinearMap part

/-- Real scalars act on complex matrices through the coercion. -/
theorem ofReal_smul_matrix (r : ℝ) (X : CMat n) : (r : ℂ) • X = r • X := by
  ext i j
  rw [Matrix.smul_apply, Matrix.smul_apply, smul_eq_mul, Complex.real_smul]

/-- The committed pinching, bundled with real scalars so that it can feed
the committed relaxation semigroup `publicRelaxTime`. The underlying
function is `EventAlgebra.partitionPinching`; both module laws are imported
from the committed complex-linear bundle. -/
noncomputable def partitionPinchingRealLinearMap
    (part : EventAlgebra.ProjectivePartition n k) : CMat n →ₗ[ℝ] CMat n where
  toFun := EventAlgebra.partitionPinching part
  map_add' X Y := by
    simpa only [EventAlgebra.partitionPinchingLinearMap_apply] using
      (EventAlgebra.partitionPinchingLinearMap part).map_add X Y
  map_smul' r X := by
    show EventAlgebra.partitionPinching part (r • X) =
      r • EventAlgebra.partitionPinching part X
    rw [← ofReal_smul_matrix r X, ← ofReal_smul_matrix r
      (EventAlgebra.partitionPinching part X)]
    simpa only [EventAlgebra.partitionPinchingLinearMap_apply] using
      (EventAlgebra.partitionPinchingLinearMap part).map_smul
        ((r : ℝ) : ℂ) X

@[simp]
theorem partitionPinchingRealLinearMap_apply
    (part : EventAlgebra.ProjectivePartition n k) (X : CMat n) :
    partitionPinchingRealLinearMap part X =
      EventAlgebra.partitionPinching part X :=
  rfl

/-- The channel agrees pointwise with the committed exponential relaxation
`publicRelaxTime` applied to the pinching: this is the literal B2 operator
identity presenting the semigroup member as a convex combination of the
identity and the conditional expectation. -/
theorem relaxationChannel_apply_eq_publicRelaxTime (gamma t : ℝ)
    (part : EventAlgebra.ProjectivePartition n k) (X : CMat n) :
    relaxationChannel gamma t part X =
      EventAlgebra.publicRelaxTime gamma t
        (partitionPinchingRealLinearMap part) X := by
  rw [EventAlgebra.publicRelaxTime, EventAlgebra.publicRelax_decomposition,
    relaxationChannel]
  simp only [LinearMap.add_apply, LinearMap.smul_apply, LinearMap.id_apply,
    partitionPinchingRealLinearMap_apply,
    EventAlgebra.partitionPinchingLinearMap_apply]
  rw [ofReal_smul_matrix, ofReal_smul_matrix]
  rfl

/-- **CPTP discharge for the relaxation semigroup.** For nonnegative rate
and time, the semigroup member is a convex combination of the identity and
the pinching, hence CPTP. -/
theorem relaxationChannel_isCPTP {gamma t : ℝ}
    (hgamma : 0 ≤ gamma) (ht : 0 ≤ t)
    (part : EventAlgebra.ProjectivePartition n k) :
    IsCPTP (relaxationChannel gamma t part) := by
  have hexp : Real.exp (-gamma * t) ≤ 1 := by
    rw [Real.exp_le_one_iff, neg_mul]
    exact neg_nonpos.mpr (mul_nonneg hgamma ht)
  exact isCPTP_id.convexCombination (partitionPinching_isCPTP part)
    (Real.exp_nonneg _) (sub_nonneg.mpr hexp) (by ring)

/-! ## Choi matrices -/

/-- The unnormalized maximally entangled matrix `sum_{a b} e_ab tensor
e_ab`, presented as an outer product. -/
noncomputable def maxEntangled (n : ℕ) : Matrix (Fin n × Fin n) (Fin n × Fin n) ℂ :=
  Matrix.vecMulVec (fun p => if p.1 = p.2 then 1 else 0)
    (star fun p : Fin n × Fin n => if p.1 = p.2 then 1 else 0)

theorem maxEntangled_posSemidef (n : ℕ) : (maxEntangled n).PosSemidef :=
  Matrix.posSemidef_vecMulVec_self_star _

theorem blockOf_maxEntangled (a b : Fin n) :
    blockOf (maxEntangled n) a b = EventAlgebra.matrixUnit a b := by
  ext i j
  by_cases hai : a = i <;> by_cases hbj : b = j <;>
    simp [maxEntangled, Matrix.vecMulVec_apply, Pi.star_apply,
      EventAlgebra.matrixUnit_apply, hai, hbj]

/-- The finite Choi matrix of a matrix map. -/
noncomputable def choiMatrix (Φ : CMat n →ₗ[ℂ] CMat n) :
    Matrix (Fin n × Fin n) (Fin n × Fin n) ℂ :=
  amplify n Φ (maxEntangled n)

theorem choiMatrix_apply (Φ : CMat n →ₗ[ℂ] CMat n) (p q : Fin n × Fin n) :
    choiMatrix Φ p q = Φ (EventAlgebra.matrixUnit p.1 q.1) p.2 q.2 := by
  rw [choiMatrix, amplify_apply, blockOf_maxEntangled]

/-- **Choi direction of Choi's theorem.** Complete positivity forces the
Choi matrix to be positive semidefinite. -/
theorem choiMatrix_posSemidef {Φ : CMat n →ₗ[ℂ] CMat n}
    (hΦ : IsCompletelyPositive Φ) : (choiMatrix Φ).PosSemidef :=
  hΦ n (maxEntangled n) (maxEntangled_posSemidef n)

/-! ## Negative control: the transpose map -/

/-- The transpose map, bundled as a complex-linear map. -/
def transposeMap (n : ℕ) : CMat n →ₗ[ℂ] CMat n where
  toFun := Matrix.transpose
  map_add' X Y := Matrix.transpose_add X Y
  map_smul' c X := Matrix.transpose_smul c X

@[simp]
theorem transposeMap_apply (X : CMat n) : transposeMap n X = Xᵀ :=
  rfl

theorem transposeMap_isPositiveMap : IsPositiveMap (transposeMap n) :=
  fun _ hX => hX.transpose

theorem transposeMap_isTracePreserving : IsTracePreserving (transposeMap n) :=
  fun X => Matrix.trace_transpose X

/-- The witness vector `e_(0,1) - e_(1,0)` refuting positivity of the
transpose Choi matrix. -/
def transposeWitness : Fin 2 × Fin 2 → ℂ := fun p =>
  if p = (0, 1) then 1 else if p = (1, 0) then -1 else 0

theorem choiMatrix_transposeMap_apply (p q : Fin 2 × Fin 2) :
    choiMatrix (transposeMap 2) p q =
      if p.1 = q.2 ∧ q.1 = p.2 then 1 else 0 := by
  rw [choiMatrix_apply, transposeMap_apply, Matrix.transpose_apply,
    EventAlgebra.matrixUnit_apply]

/-- The exact 4x4 computation: the transpose Choi matrix has expectation
`-2` on the witness vector. -/
theorem transposeWitness_expectation :
    star transposeWitness ⬝ᵥ
      (choiMatrix (transposeMap 2) *ᵥ transposeWitness) = -2 := by
  simp only [dotProduct, Matrix.mulVec, choiMatrix_transposeMap_apply,
    transposeWitness, Pi.star_apply, Fintype.sum_prod_type,
    Fin.sum_univ_two]
  norm_num

/-- **Negative control.** The transpose map on `M_2` fails complete
positivity: its Choi matrix has a negative expectation, contradicting the
positivity forced by any amplification. -/
theorem transposeMap_not_completelyPositive :
    ¬ IsCompletelyPositive (transposeMap 2) := by
  intro h
  have hval := (choiMatrix_posSemidef h).dotProduct_mulVec_nonneg
    transposeWitness
  rw [transposeWitness_expectation] at hval
  have hre := (Complex.nonneg_iff.mp hval).1
  norm_num at hre

/-- The transpose map on `M_2` is positive and trace preserving and fails
complete positivity: positivity does not imply complete positivity, so the
CPTP discharge above carries strictly more content than positivity plus
trace preservation. -/
theorem transposeMap_positive_tracePreserving_not_CP :
    IsPositiveMap (transposeMap 2) ∧ IsTracePreserving (transposeMap 2) ∧
      ¬ IsCompletelyPositive (transposeMap 2) :=
  ⟨transposeMap_isPositiveMap, transposeMap_isTracePreserving,
    transposeMap_not_completelyPositive⟩

end OPH.Dynamics

#print axioms OPH.Dynamics.isCompletelyPositive_of_kraus
#print axioms OPH.Dynamics.IsCompletelyPositive.isPositiveMap
#print axioms OPH.Dynamics.IsCPTP.convexCombination
#print axioms OPH.Dynamics.partitionAverage_isCPTP
#print axioms OPH.Dynamics.partitionPinching_isCPTP
#print axioms OPH.Dynamics.relaxationChannel_apply_eq_publicRelaxTime
#print axioms OPH.Dynamics.relaxationChannel_isCPTP
#print axioms OPH.Dynamics.choiMatrix_posSemidef
#print axioms OPH.Dynamics.transposeMap_not_completelyPositive
