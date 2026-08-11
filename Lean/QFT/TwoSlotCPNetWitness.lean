import QFT.CPRestrictionNet

/-!
# The two-slot inhabitant of the conditional-expectation net

This module inhabits the redesigned `CPRegionalNet` interface with the
designated observer pair of the committed source-operator payload.  The
ambient algebra consists of the matrices indexed by the joint slot
space `Fin 13 × Fin 14`;
the four regions form the diamond `bot < left, right < top` with

* `localAlgebra left` the left slot, equal to the image of observer
  86's source-generated algebra under the slot embedding;
* `localAlgebra right` the right slot, equal to the image of observer
  88's source-generated algebra;
* `localAlgebra bot` the scalars and `localAlgebra top` everything;
* locality on the declared disjoint pair by slot commutation;
* conditional expectations: normalized full trace at `bot`, the two
  partial-trace slot expectations, and the identity at `top`, with
  membership, fixing, positivity, trace, and tower receipts;
* the genuine joint coverage law: the two slots generate the top
  algebra.

The checkpoint local-channel receipt holds in Heisenberg form: pinching
by observer 86's lifted checkpoint projectors fixes every right-slot
observable pointwise.  The redesign rationale is instantiated: the left
regional algebra carries an embedded full 13-by-13 matrix algebra, so
no unital algebra homomorphism from it reaches the scalars, and the
star-homomorphic restriction interface is uninhabitable over this
diamond.

**Boundary.**  The slot assembly and region labels are declared
postprocessors over the post-hoc payload transcriptions; the witness
carries no tower stage, spacetime, channel semantics, instrument, clock,
or physical claim, and the extraction is ineligible as validation.
Issue #692 remains open for source selection of the joint carrier and a
coherent CP-net correlation or descent receipt.  Closed issue #712
consumes the older finite-net interface and does not supply that bridge.
-/

namespace OPH.QFT

open Matrix
open Kronecker
open scoped ComplexOrder

noncomputable section

/-! ## The mirror slot expectation -/

variable {α β : Type*} [Fintype α] [Fintype β] [DecidableEq α] [DecidableEq β]

/-- Partial trace over the second factor. -/
def ptraceSnd (M : Matrix (α × β) (α × β) ℂ) : Matrix α α ℂ :=
  fun a a' => ∑ b, M (a, b) (a', b)

omit [DecidableEq α] [DecidableEq β] in
/-- The partial trace over the right slot reproduces the trace. -/
theorem ptraceSnd_trace (M : Matrix (α × β) (α × β) ℂ) :
    (ptraceSnd M).trace = M.trace := by
  unfold ptraceSnd Matrix.trace
  rw [Fintype.sum_prod_type]
  rfl

omit [Fintype α] [DecidableEq α] in
/-- The partial trace of a left-slot element. -/
theorem ptraceSnd_slotLeft (A : Matrix α α ℂ) :
    ptraceSnd (A ⊗ₖ (1 : Matrix β β ℂ)) = (Fintype.card β : ℂ) • A := by
  ext a a'
  unfold ptraceSnd
  simp only [Matrix.kroneckerMap_apply, Matrix.smul_apply, smul_eq_mul]
  have hred : (∑ b : β, A a a' * (1 : Matrix β β ℂ) b b) =
      (Fintype.card β : ℂ) * A a a' := by
    simp [Matrix.one_apply_eq, Finset.sum_const, Finset.card_univ,
      nsmul_eq_mul, mul_comm]
  exact hred

omit [Fintype α] [DecidableEq β] in
/-- The partial trace of a right-slot element. -/
theorem ptraceSnd_slotRight (B : Matrix β β ℂ) :
    ptraceSnd ((1 : Matrix α α ℂ) ⊗ₖ B) = B.trace • (1 : Matrix α α ℂ) := by
  ext a a'
  unfold ptraceSnd
  simp only [Matrix.kroneckerMap_apply, Matrix.smul_apply, smul_eq_mul]
  have htr : B.trace = ∑ b : β, B b b := rfl
  have hred : (∑ b : β, (1 : Matrix α α ℂ) a a' * B b b) =
      B.trace * (1 : Matrix α α ℂ) a a' := by
    rw [htr, Finset.sum_mul]
    exact Finset.sum_congr rfl fun b _ => mul_comm _ _
  exact hred

/-- The conditional expectation onto the left slot. -/
def leftSlotExpectation :
    Matrix (α × β) (α × β) ℂ →ₗ[ℂ] Matrix (α × β) (α × β) ℂ where
  toFun M := (Fintype.card β : ℂ)⁻¹ • (ptraceSnd M ⊗ₖ (1 : Matrix β β ℂ))
  map_add' M N := by
    unfold ptraceSnd
    rw [← smul_add, ← Matrix.add_kronecker]
    congr 2
    ext a a'
    simp [Finset.sum_add_distrib]
  map_smul' c M := by
    unfold ptraceSnd
    rw [RingHom.id_apply, smul_comm]
    congr 1
    rw [← Matrix.smul_kronecker]
    congr 1
    ext a a'
    simp [Finset.mul_sum]

omit [DecidableEq β] in
/-- The right-slot conjugation summand in explicit entry form. -/
theorem slotRight_conj_apply (Kk : Matrix β β ℂ)
    (M : Matrix (α × β) (α × β) ℂ) (z z' : α × β) :
    (((1 : Matrix α α ℂ) ⊗ₖ Kk) * M * ((1 : Matrix α α ℂ) ⊗ₖ Kk)ᴴ) z z'
      = ∑ c, ∑ c', Kk z.2 c * M (z.1, c) (z'.1, c')
          * star (Kk z'.2 c') := by
  rw [Matrix.conjTranspose_kronecker, Matrix.conjTranspose_one]
  rw [Matrix.mul_apply]
  have hmid : ∀ w : α × β,
      (((1 : Matrix α α ℂ) ⊗ₖ Kk) * M) z w
        = ∑ c, Kk z.2 c * M (z.1, c) w := by
    intro w
    rw [Matrix.mul_apply, Fintype.sum_prod_type]
    rw [Finset.sum_eq_single z.1
      (fun a _ ha => by
        apply Finset.sum_eq_zero
        intro c _
        simp [Matrix.kroneckerMap_apply, Ne.symm ha])
      (fun h => absurd (Finset.mem_univ _) h)]
    apply Finset.sum_congr rfl
    intro c _
    simp [Matrix.kroneckerMap_apply]
  calc ∑ w : α × β, (((1 : Matrix α α ℂ) ⊗ₖ Kk) * M) z w *
        ((1 : Matrix α α ℂ) ⊗ₖ Kkᴴ) w z'
      = ∑ w : α × β, (∑ c, Kk z.2 c * M (z.1, c) w) *
          ((1 : Matrix α α ℂ) ⊗ₖ Kkᴴ) w z' :=
        Finset.sum_congr rfl fun w _ => by rw [hmid w]
    _ = ∑ c', (∑ c, Kk z.2 c * M (z.1, c) (z'.1, c')) *
          star (Kk z'.2 c') := by
        rw [Fintype.sum_prod_type]
        rw [Finset.sum_eq_single z'.1
          (fun a _ ha => by
            apply Finset.sum_eq_zero
            intro c' _
            simp [Matrix.kroneckerMap_apply, ha])
          (fun h => absurd (Finset.mem_univ _) h)]
        apply Finset.sum_congr rfl
        intro c' _
        simp [Matrix.kroneckerMap_apply, Matrix.conjTranspose_apply]
    _ = ∑ c', ∑ c, Kk z.2 c * M (z.1, c) (z'.1, c') * star (Kk z'.2 c') := by
        refine Finset.sum_congr rfl fun c' _ => ?_
        rw [Finset.sum_mul]
    _ = ∑ c, ∑ c', Kk z.2 c * M (z.1, c) (z'.1, c') * star (Kk z'.2 c') :=
        Finset.sum_comm

/-- The matrix-unit Kraus identity of the left expectation. -/
theorem sum_single_conj_eq_slotLeft_ptraceSnd
    (M : Matrix (α × β) (α × β) ℂ) :
    ∑ p : β × β,
        ((1 : Matrix α α ℂ) ⊗ₖ Matrix.single p.1 p.2 (1 : ℂ)) * M *
          ((1 : Matrix α α ℂ) ⊗ₖ Matrix.single p.1 p.2 (1 : ℂ))ᴴ =
      ptraceSnd M ⊗ₖ (1 : Matrix β β ℂ) := by
  ext z z'
  rw [Matrix.sum_apply]
  have hterm : ∀ p : β × β,
      (((1 : Matrix α α ℂ) ⊗ₖ Matrix.single p.1 p.2 (1 : ℂ)) * M *
        ((1 : Matrix α α ℂ) ⊗ₖ Matrix.single p.1 p.2 (1 : ℂ))ᴴ) z z' =
      if p.1 = z.2 ∧ p.1 = z'.2 then M (z.1, p.2) (z'.1, p.2) else 0 := by
    intro p
    rw [slotRight_conj_apply]
    rw [Finset.sum_eq_single p.2
      (fun c _ hc => by
        have hzero : Matrix.single p.1 p.2 (1 : ℂ) z.2 c = 0 := by
          simp [Matrix.single, Ne.symm hc]
        simp [hzero])
      (fun h => absurd (Finset.mem_univ _) h)]
    rw [Finset.sum_eq_single p.2
      (fun c' _ hc' => by
        have hzero : Matrix.single p.1 p.2 (1 : ℂ) z'.2 c' = 0 := by
          simp [Matrix.single, Ne.symm hc']
        simp [hzero])
      (fun h => absurd (Finset.mem_univ _) h)]
    by_cases h1 : p.1 = z.2 <;> by_cases h2 : p.1 = z'.2 <;>
      simp [Matrix.single, h1, h2]
  rw [Finset.sum_congr rfl fun p _ => hterm p]
  by_cases hz : z.2 = z'.2
  · simp [Fintype.sum_prod_type, hz, Finset.sum_ite_eq', ptraceSnd,
      Matrix.kroneckerMap_apply, Matrix.one_apply]
  · have hfalse : ∀ p1 : β, ¬ (p1 = z.2 ∧ p1 = z'.2) := by
      intro p1 hp
      exact hz (hp.1 ▸ hp.2 ▸ rfl)
    simp [hfalse, ptraceSnd, Matrix.kroneckerMap_apply, hz]

/-- The left expectation preserves positive semidefiniteness. -/
theorem leftSlotExpectation_posSemidef [Nonempty β]
    {M : Matrix (α × β) (α × β) ℂ} (hM : M.PosSemidef) :
    (leftSlotExpectation (α := α) (β := β) M).PosSemidef := by
  have h : leftSlotExpectation (α := α) (β := β) M =
      (Fintype.card β : ℂ)⁻¹ • ∑ p : β × β,
        ((1 : Matrix α α ℂ) ⊗ₖ Matrix.single p.1 p.2 (1 : ℂ)) * M *
          ((1 : Matrix α α ℂ) ⊗ₖ Matrix.single p.1 p.2 (1 : ℂ))ᴴ := by
    rw [sum_single_conj_eq_slotLeft_ptraceSnd]
    rfl
  rw [h]
  have hnn : (0 : ℂ) ≤ (Fintype.card β : ℂ)⁻¹ := by
    rw [← Complex.ofReal_natCast, ← Complex.ofReal_inv]
    exact_mod_cast inv_nonneg.mpr (Nat.cast_nonneg _)
  exact (Matrix.posSemidef_sum _ fun p _ =>
    hM.mul_mul_conjTranspose_same _).smul hnn


/-! ## Slot-expectation receipts -/

omit [Fintype α] [DecidableEq α] in
/-- The left expectation fixes the left slot pointwise. -/
theorem leftSlotExpectation_fixes_left [Nonempty β] (A : Matrix α α ℂ) :
    leftSlotExpectation (A ⊗ₖ (1 : Matrix β β ℂ)) =
      A ⊗ₖ (1 : Matrix β β ℂ) := by
  have hcard : (Fintype.card β : ℂ) ≠ 0 :=
    Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  have h : leftSlotExpectation (A ⊗ₖ (1 : Matrix β β ℂ)) =
      (Fintype.card β : ℂ)⁻¹ •
        (ptraceSnd (A ⊗ₖ (1 : Matrix β β ℂ)) ⊗ₖ (1 : Matrix β β ℂ)) := rfl
  rw [h, ptraceSnd_slotLeft, Matrix.smul_kronecker, smul_smul,
    inv_mul_cancel₀ hcard, one_smul]

omit [DecidableEq α] in
/-- The left expectation preserves the trace. -/
theorem leftSlotExpectation_trace [Nonempty β]
    (M : Matrix (α × β) (α × β) ℂ) :
    (leftSlotExpectation M).trace = M.trace := by
  have hcard : (Fintype.card β : ℂ) ≠ 0 :=
    Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  have h : leftSlotExpectation (α := α) (β := β) M =
      (Fintype.card β : ℂ)⁻¹ • (ptraceSnd M ⊗ₖ (1 : Matrix β β ℂ)) := rfl
  rw [h, Matrix.trace_smul, Matrix.trace_kronecker, Matrix.trace_one,
    ptraceSnd_trace, smul_eq_mul,
    mul_comm M.trace ((Fintype.card β : ℕ) : ℂ), ← mul_assoc,
    inv_mul_cancel₀ hcard, one_mul]

/-- The scalar conditional expectation: the normalized full trace. -/
def scalarExpectation (γ : Type*) [Fintype γ] [DecidableEq γ] :
    Matrix γ γ ℂ →ₗ[ℂ] Matrix γ γ ℂ where
  toFun M := (Fintype.card γ : ℂ)⁻¹ • M.trace • (1 : Matrix γ γ ℂ)
  map_add' M N := by
    rw [Matrix.trace_add, add_smul, smul_add]
  map_smul' c M := by
    rw [Matrix.trace_smul, RingHom.id_apply, smul_comm c, smul_assoc]

variable {γ : Type*} [Fintype γ] [DecidableEq γ]

theorem scalarExpectation_mem_bot (M : Matrix γ γ ℂ) :
    scalarExpectation γ M ∈ (⊥ : StarSubalgebra ℂ (Matrix γ γ ℂ)) := by
  rw [StarSubalgebra.mem_bot]
  refine ⟨(Fintype.card γ : ℂ)⁻¹ * M.trace, ?_⟩
  have h : scalarExpectation γ M =
      (Fintype.card γ : ℂ)⁻¹ • M.trace • (1 : Matrix γ γ ℂ) := rfl
  rw [Algebra.algebraMap_eq_smul_one, h, smul_smul]

theorem scalarExpectation_fixes [Nonempty γ] (M : Matrix γ γ ℂ)
    (hM : M ∈ (⊥ : StarSubalgebra ℂ (Matrix γ γ ℂ))) :
    scalarExpectation γ M = M := by
  rw [StarSubalgebra.mem_bot] at hM
  obtain ⟨c, rfl⟩ := hM
  have hcard : (Fintype.card γ : ℂ) ≠ 0 :=
    Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  have htr : (algebraMap ℂ (Matrix γ γ ℂ) c).trace =
      c * (Fintype.card γ : ℂ) := by
    rw [Algebra.algebraMap_eq_smul_one, Matrix.trace_smul, Matrix.trace_one,
      smul_eq_mul]
  have h : scalarExpectation γ (algebraMap ℂ (Matrix γ γ ℂ) c) =
      (Fintype.card γ : ℂ)⁻¹ •
        (algebraMap ℂ (Matrix γ γ ℂ) c).trace • (1 : Matrix γ γ ℂ) := rfl
  rw [h, htr, Algebra.algebraMap_eq_smul_one, smul_smul]
  congr 1
  rw [mul_comm c, ← mul_assoc, inv_mul_cancel₀ hcard, one_mul]

theorem scalarExpectation_posSemidef [Nonempty γ] {M : Matrix γ γ ℂ}
    (hM : M.PosSemidef) : (scalarExpectation γ M).PosSemidef := by
  have h0 : scalarExpectation γ M =
      (Fintype.card γ : ℂ)⁻¹ • M.trace • (1 : Matrix γ γ ℂ) := rfl
  rw [h0, smul_smul]
  have hnn : (0 : ℂ) ≤ (Fintype.card γ : ℂ)⁻¹ * M.trace := by
    refine mul_nonneg ?_ hM.trace_nonneg
    rw [← Complex.ofReal_natCast, ← Complex.ofReal_inv]
    exact_mod_cast inv_nonneg.mpr (Nat.cast_nonneg _)
  exact Matrix.PosSemidef.one.smul hnn

theorem scalarExpectation_trace (M : Matrix γ γ ℂ) :
    (scalarExpectation γ M).trace = M.trace := by
  have h : scalarExpectation γ M =
      (Fintype.card γ : ℂ)⁻¹ • M.trace • (1 : Matrix γ γ ℂ) := rfl
  rw [h, Matrix.trace_smul, Matrix.trace_smul, Matrix.trace_one,
    smul_eq_mul, smul_eq_mul]
  by_cases hc : Fintype.card γ = 0
  · have hempty : IsEmpty γ := Fintype.card_eq_zero_iff.mp hc
    have hM0 : M.trace = 0 := by
      rw [Matrix.trace]
      simp
    rw [hM0, hc]
    simp
  · have hcard : (Fintype.card γ : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hc
    rw [mul_comm M.trace, ← mul_assoc, inv_mul_cancel₀ hcard, one_mul]

/-- The scalar expectation depends on a matrix only through its trace. -/
theorem scalarExpectation_of_trace_eq {M N : Matrix γ γ ℂ}
    (h : M.trace = N.trace) :
    scalarExpectation γ M = scalarExpectation γ N := by
  have hM : scalarExpectation γ M =
      (Fintype.card γ : ℂ)⁻¹ • M.trace • (1 : Matrix γ γ ℂ) := rfl
  have hN : scalarExpectation γ N =
      (Fintype.card γ : ℂ)⁻¹ • N.trace • (1 : Matrix γ γ ℂ) := rfl
  rw [hM, hN, h]

/-! ## The diamond region system -/

/-- The four regions of the two-slot diamond. -/
inductive TwoSlotRegion : Type
  | bot : TwoSlotRegion
  | left : TwoSlotRegion
  | right : TwoSlotRegion
  | top : TwoSlotRegion
deriving DecidableEq, Fintype

namespace TwoSlotRegion

/-- The diamond order, as a Boolean relation. -/
def le : TwoSlotRegion → TwoSlotRegion → Bool
  | bot, _ => true
  | _, top => true
  | left, left => true
  | right, right => true
  | _, _ => false

/-- The diamond meet. -/
def meet : TwoSlotRegion → TwoSlotRegion → TwoSlotRegion
  | top, V => V
  | U, top => U
  | left, left => left
  | right, right => right
  | _, _ => bot

/-- The declared disjoint pair. -/
def disj : TwoSlotRegion → TwoSlotRegion → Bool
  | left, right => true
  | right, left => true
  | _, _ => false

end TwoSlotRegion

/-- The left expectation scalarises the right slot: a right-slot
observable is replaced by its normalized trace. -/
theorem leftSlotExpectation_scalarises_right {α β : Type*} [Fintype α]
    [Fintype β] [DecidableEq α] [DecidableEq β] [Nonempty β]
    (B : Matrix β β ℂ) :
    leftSlotExpectation ((1 : Matrix α α ℂ) ⊗ₖ B) =
      ((Fintype.card β : ℂ)⁻¹ * B.trace) • 1 := by
  have h : leftSlotExpectation ((1 : Matrix α α ℂ) ⊗ₖ B) =
      (Fintype.card β : ℂ)⁻¹ •
        (ptraceSnd ((1 : Matrix α α ℂ) ⊗ₖ B) ⊗ₖ (1 : Matrix β β ℂ)) := rfl
  rw [h, ptraceSnd_slotRight, Matrix.smul_kronecker, smul_smul,
    Matrix.one_kronecker_one]

/-- The left expectation is idempotent. -/
theorem leftSlotExpectation_idem {α β : Type*} [Fintype α] [Fintype β]
    [DecidableEq α] [DecidableEq β] [Nonempty β]
    (M : Matrix (α × β) (α × β) ℂ) :
    leftSlotExpectation (α := α) (β := β)
        (leftSlotExpectation M) = leftSlotExpectation M := by
  have h : leftSlotExpectation (α := α) (β := β) M =
      (Fintype.card β : ℂ)⁻¹ •
        (ptraceSnd M ⊗ₖ (1 : Matrix β β ℂ)) := rfl
  rw [h, map_smul, leftSlotExpectation_fixes_left, ← h]

namespace TwoSlotRegion

theorem le_refl' : ∀ U : TwoSlotRegion, le U U = true := by decide
theorem le_trans' : ∀ U V W : TwoSlotRegion,
    le U V = true → le V W = true → le U W = true := by decide
theorem le_antisymm' : ∀ U V : TwoSlotRegion,
    le U V = true → le V U = true → U = V := by decide
theorem meet_le_left' : ∀ U V : TwoSlotRegion,
    le (meet U V) U = true := by decide
theorem meet_le_right' : ∀ U V : TwoSlotRegion,
    le (meet U V) V = true := by decide
theorem le_meet' : ∀ W U V : TwoSlotRegion,
    le W U = true → le W V = true → le W (meet U V) = true := by decide
theorem disj_symm' : ∀ U V : TwoSlotRegion,
    disj U V = true → disj V U = true := by decide
theorem disj_irrefl' : ∀ U : TwoSlotRegion, ¬ disj U U = true := by decide

end TwoSlotRegion

/-! ## The two-slot witness -/

/-- The ambient index of the designated pair. -/
abbrev PairIndex := Fin 13 × Fin 14

/-- The regional algebra assignment of the diamond. -/
def twoSlotAlgebra : TwoSlotRegion → StarSubalgebra ℂ (Matrix PairIndex PairIndex ℂ)
  | TwoSlotRegion.bot => ⊥
  | TwoSlotRegion.left => (slotLeft (Fin 14) (α := Fin 13)).range
  | TwoSlotRegion.right => (slotRight (Fin 13) (β := Fin 14)).range
  | TwoSlotRegion.top => ⊤

/-- The conditional-expectation assignment of the diamond. -/
def twoSlotExpect : TwoSlotRegion →
    Matrix PairIndex PairIndex ℂ →ₗ[ℂ] Matrix PairIndex PairIndex ℂ
  | TwoSlotRegion.bot => scalarExpectation PairIndex
  | TwoSlotRegion.left => leftSlotExpectation
  | TwoSlotRegion.right => rightSlotExpectation
  | TwoSlotRegion.top => LinearMap.id

/-- Left-slot membership of the left expectation value. -/
theorem leftSlotExpectation_mem (M : Matrix PairIndex PairIndex ℂ) :
    leftSlotExpectation M ∈ (slotLeft (Fin 14) (α := Fin 13)).range := by
  have h : leftSlotExpectation (α := Fin 13) (β := Fin 14) M =
      (Fintype.card (Fin 14) : ℂ)⁻¹ •
        (ptraceSnd M ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ)) := rfl
  rw [h]
  exact SMulMemClass.smul_mem _ ⟨ptraceSnd M, rfl⟩

/-- Right-slot membership of the right expectation value. -/
theorem rightSlotExpectation_mem (M : Matrix PairIndex PairIndex ℂ) :
    rightSlotExpectation M ∈ (slotRight (Fin 13) (β := Fin 14)).range := by
  have h : rightSlotExpectation (α := Fin 13) (β := Fin 14) M =
      (Fintype.card (Fin 13) : ℂ)⁻¹ •
        ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ OPH.Locality.ptraceFst M) := rfl
  rw [h]
  exact SMulMemClass.smul_mem _ ⟨OPH.Locality.ptraceFst M, rfl⟩

/-- **The two-slot inhabitant.**  The diamond of the designated pair
inhabits the conditional-expectation net interface. -/
def twoSlotNet : CPRegionalNet PairIndex where
  Region := TwoSlotRegion
  regionFintype := inferInstance
  regionNonempty := ⟨TwoSlotRegion.bot⟩
  regionLE U V := TwoSlotRegion.le U V = true
  regionLE_refl := TwoSlotRegion.le_refl'
  regionLE_trans {U V W} := TwoSlotRegion.le_trans' U V W
  regionLE_antisymm {U V} := TwoSlotRegion.le_antisymm' U V
  overlap := TwoSlotRegion.meet
  overlap_le_left := TwoSlotRegion.meet_le_left'
  overlap_le_right := TwoSlotRegion.meet_le_right'
  le_overlap {W U V} := TwoSlotRegion.le_meet' W U V
  disjoint U V := TwoSlotRegion.disj U V = true
  disjoint_symm {U V} := TwoSlotRegion.disj_symm' U V
  disjoint_irrefl := TwoSlotRegion.disj_irrefl'
  localAlgebra := twoSlotAlgebra
  isotony {U V} hUV := by
    cases U <;> cases V <;> (try (simp [TwoSlotRegion.le] at hUV))
    · exact le_refl _
    · exact bot_le
    · exact bot_le
    · exact bot_le
    · exact le_refl _
    · exact le_top
    · exact le_refl _
    · exact le_top
    · exact le_refl _
  locality {U V} hUV := by
    cases U <;> cases V <;> (try (simp [TwoSlotRegion.disj] at hUV))
    · exact fun X hX Y hY => slot_commute hX hY
    · exact fun X hX Y hY => (slot_commute hY hX).symm
  expect := twoSlotExpect
  expect_mem U X := by
    cases U
    · exact scalarExpectation_mem_bot X
    · exact leftSlotExpectation_mem X
    · exact rightSlotExpectation_mem X
    · trivial
  expect_fixes U X hX := by
    cases U
    · exact scalarExpectation_fixes X hX
    · obtain ⟨A, rfl⟩ := hX
      exact leftSlotExpectation_fixes_left A
    · obtain ⟨B, rfl⟩ := hX
      exact rightSlotExpectation_fixes_right B
    · rfl
  expect_posSemidef U {X} hX := by
    cases U
    · exact scalarExpectation_posSemidef hX
    · exact leftSlotExpectation_posSemidef hX
    · exact rightSlotExpectation_posSemidef hX
    · exact hX
  expect_trace U X := by
    cases U
    · exact scalarExpectation_trace X
    · exact leftSlotExpectation_trace X
    · exact rightSlotExpectation_trace X
    · rfl
  expect_tower {U V} hUV X := by
    cases U <;> cases V <;> (try (simp [TwoSlotRegion.le] at hUV))
    · exact scalarExpectation_of_trace_eq (scalarExpectation_trace X)
    · exact scalarExpectation_of_trace_eq (leftSlotExpectation_trace X)
    · exact scalarExpectation_of_trace_eq (rightSlotExpectation_trace X)
    · rfl
    · exact leftSlotExpectation_idem X
    · rfl
    · exact rightSlotExpectation_idem X
    · rfl
    · rfl
  generating := {TwoSlotRegion.left, TwoSlotRegion.right}
  coverage := by
    have hset : (⋃ U ∈ ({TwoSlotRegion.left, TwoSlotRegion.right} :
        Finset TwoSlotRegion),
          (twoSlotAlgebra U : Set (Matrix PairIndex PairIndex ℂ))) =
        (((slotLeft (Fin 14) (α := Fin 13)).range :
            Set (Matrix PairIndex PairIndex ℂ)) ∪
          ((slotRight (Fin 13) (β := Fin 14)).range :
            Set (Matrix PairIndex PairIndex ℂ))) := by
      ext x
      simp only [Finset.mem_insert, Finset.mem_singleton, Set.mem_iUnion,
        Set.mem_union, exists_prop]
      constructor
      · rintro ⟨U, (rfl | rfl), hx⟩
        · exact Or.inl hx
        · exact Or.inr hx
      · rintro (hx | hx)
        · exact ⟨TwoSlotRegion.left, Or.inl rfl, hx⟩
        · exact ⟨TwoSlotRegion.right, Or.inr rfl, hx⟩
    rw [hset]
    exact slot_ranges_generate_top

/-! ## Marginal nonuniqueness control -/

/-- A nonzero off-diagonal correlation matrix unit on the two-slot
ambient space.  It is invisible to both one-slot conditional
expectations. -/
def offDiagonalCorrelation : Matrix PairIndex PairIndex ℂ :=
  Matrix.single ((0, 0) : PairIndex) ((1, 1) : PairIndex) 1

/-- The correlation matrix unit is genuinely nonzero. -/
theorem offDiagonalCorrelation_ne_zero : offDiagonalCorrelation ≠ 0 := by
  intro h
  have hentry := congrFun
    (congrFun h ((0, 0) : PairIndex)) ((1, 1) : PairIndex)
  simp [offDiagonalCorrelation] at hentry

/-- The right-factor partial trace of the correlation unit vanishes. -/
theorem offDiagonalCorrelation_ptraceSnd :
    ptraceSnd offDiagonalCorrelation = 0 := by
  ext a a'
  simp [ptraceSnd, offDiagonalCorrelation, Matrix.single]

/-- The left-factor partial trace of the correlation unit vanishes. -/
theorem offDiagonalCorrelation_ptraceFst :
    OPH.Locality.ptraceFst offDiagonalCorrelation = 0 := by
  ext b b'
  simp [OPH.Locality.ptraceFst, offDiagonalCorrelation, Matrix.single]

/-- The left-slot conditional expectation erases the correlation unit. -/
theorem leftSlotExpectation_offDiagonalCorrelation :
    leftSlotExpectation offDiagonalCorrelation = 0 := by
  change (14 : ℂ)⁻¹ •
    (ptraceSnd offDiagonalCorrelation ⊗ₖ
      (1 : Matrix (Fin 14) (Fin 14) ℂ)) = 0
  rw [offDiagonalCorrelation_ptraceSnd]
  simp

/-- The right-slot conditional expectation erases the correlation unit. -/
theorem rightSlotExpectation_offDiagonalCorrelation :
    rightSlotExpectation offDiagonalCorrelation = 0 := by
  change (13 : ℂ)⁻¹ •
    ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
      OPH.Locality.ptraceFst offDiagonalCorrelation) = 0
  rw [offDiagonalCorrelation_ptraceFst]
  simp

/-- Both regional conditional expectations erase the same nonzero
correlation unit.  Consequently the pair of regional expectations is not
jointly injective: coverage by the two factor algebras does not imply unique
gluing or reconstruction from the two marginals. -/
theorem slotExpectations_not_jointly_injective :
    ¬ Function.Injective
      (fun M : Matrix PairIndex PairIndex ℂ =>
        (leftSlotExpectation M, rightSlotExpectation M)) := by
  intro hinj
  have hpair :
      (leftSlotExpectation offDiagonalCorrelation,
          rightSlotExpectation offDiagonalCorrelation) =
        (leftSlotExpectation 0, rightSlotExpectation 0) := by
    rw [leftSlotExpectation_offDiagonalCorrelation,
      rightSlotExpectation_offDiagonalCorrelation]
    simp
  exact offDiagonalCorrelation_ne_zero (hinj hpair)

/-! ## Local-channel and rationale receipts -/

/-- **Local-channel receipt.**  Pinching by observer 86's lifted
checkpoint projectors fixes every right-slot observable pointwise: the
channel acts inside the left region and leaves the disjoint region's
observables exactly unchanged in Heisenberg form. -/
theorem checkpoint_pinch_fixes_right (B : Matrix (Fin 14) (Fin 14) ℂ) :
    ∑ v : Fin 4,
        (chk86 v ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ)) *
          ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ B) *
            (chk86 v ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ))ᴴ =
      (1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ B := by
  have hterm : ∀ v : Fin 4,
      (chk86 v ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ)) *
        ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ B) *
          (chk86 v ⊗ₖ (1 : Matrix (Fin 14) (Fin 14) ℂ))ᴴ =
      chk86 v ⊗ₖ B := by
    intro v
    rw [Matrix.conjTranspose_kronecker, Matrix.conjTranspose_one,
      chk86_conjTranspose, ← Matrix.mul_kronecker_mul,
      ← Matrix.mul_kronecker_mul, mul_one, one_mul, mul_one, chk86_idem]
  rw [Finset.sum_congr rfl fun v _ => hterm v]
  have hsum : (∑ v : Fin 4, chk86 v ⊗ₖ B) =
      (∑ v : Fin 4, chk86 v) ⊗ₖ B := by
    ext z z'
    rw [Matrix.sum_apply]
    simp [Matrix.kroneckerMap_apply, Matrix.sum_apply, Finset.sum_mul]
  rw [hsum, chk86_sum]

/-- Every checkpoint projector is proper: none vanishes and none is the
identity, so the pinch channel is a non-identity local operation. -/
theorem chk86_proper : ∀ v : Fin 4,
    (∃ s : Fin 13, obs86.tbl s 1 = (v : ℕ)) ∧
      (∃ s : Fin 13, obs86.tbl s 1 ≠ (v : ℕ)) := by decide

/-- The slot embedding into the left regional algebra, as an algebra
homomorphism. -/
def slotLeftIntoRange :
    Matrix (Fin 13) (Fin 13) ℂ →ₐ[ℂ]
      (twoSlotAlgebra TwoSlotRegion.left) where
  toFun A := ⟨slotLeft (Fin 14) A, ⟨A, rfl⟩⟩
  map_one' := Subtype.ext (map_one (slotLeft (Fin 14)))
  map_mul' A B := Subtype.ext (map_mul (slotLeft (Fin 14)) A B)
  map_zero' := Subtype.ext (map_zero (slotLeft (Fin 14)))
  map_add' A B := Subtype.ext (map_add (slotLeft (Fin 14)) A B)
  commutes' c := Subtype.ext (AlgHomClass.commutes (slotLeft (Fin 14)) c)

/-- **Redesign rationale, instantiated.**  No unital algebra
homomorphism carries the left regional algebra to the scalars, so the
star-homomorphic restriction interface admits no restriction from the
left region to the scalar bottom of this diamond. -/
theorem twoSlot_left_no_scalar_hom
    (φ : (twoSlotAlgebra TwoSlotRegion.left) →ₐ[ℂ] ℂ) : False :=
  no_scalar_restriction_of_matrix_factor (by norm_num)
    (twoSlotAlgebra TwoSlotRegion.left) slotLeftIntoRange φ

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.ptraceSnd_trace
#print axioms OPH.QFT.sum_single_conj_eq_slotLeft_ptraceSnd
#print axioms OPH.QFT.leftSlotExpectation_posSemidef
#print axioms OPH.QFT.leftSlotExpectation_fixes_left
#print axioms OPH.QFT.leftSlotExpectation_trace
#print axioms OPH.QFT.scalarExpectation_fixes
#print axioms OPH.QFT.scalarExpectation_posSemidef
#print axioms OPH.QFT.twoSlotNet
#print axioms OPH.QFT.checkpoint_pinch_fixes_right
#print axioms OPH.QFT.twoSlot_left_no_scalar_hom
#print axioms OPH.QFT.ptraceSnd_slotRight
#print axioms OPH.QFT.leftSlotExpectation_scalarises_right
#print axioms OPH.QFT.chk86_proper
#print axioms OPH.QFT.offDiagonalCorrelation_ne_zero
#print axioms OPH.QFT.leftSlotExpectation_offDiagonalCorrelation
#print axioms OPH.QFT.rightSlotExpectation_offDiagonalCorrelation
#print axioms OPH.QFT.slotExpectations_not_jointly_injective
