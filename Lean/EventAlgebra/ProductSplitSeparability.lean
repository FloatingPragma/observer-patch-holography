import EventAlgebra.SourceReachabilityDelimitation

set_option autoImplicit false

/-!
# Product-split separability and the slot-local CHSH bound

Issue #730 OL-C3, continuation of the record-diagonal delimitation of
`EventAlgebra.SourceReachabilityDelimitation`.  That module proves the
classical CHSH bound `2` for reachable record-diagonal states against four
JOINTLY record-diagonal readouts, and the obstruction that no
record-diagonality-preserving state map reaches the declared Bell state.
Its boundary clause states that a slot-local separability theorem for
arbitrary settings is not proved there.  This module supplies that
theorem on the PR-44 product split, together with the counter-witness
that fixes which hypothesis carries the bound.

## The split

The carrier is `α × β` for finite types `α`, `β`; the two slots are the
fixed embeddings `slotLeft β : A ↦ A ⊗ₖ 1` and `slotRight α : B ↦ 1 ⊗ₖ B`
of `QFT.JointSlotFactorisation`, the declared PR-44 tensor split.  The
record basis of the product carrier is the product basis.  A state is
record-diagonal when it is diagonal in this basis (`IsRecordDiagonal`),
and satisfies the record law when in addition its diagonal is a
probability vector (`IsRecordLaw`); both predicates are those of the
source-reachability module.

## Results

1. **Separability.**  Every record-diagonal state on `α × β` is the
   explicit convex combination of the product point states
   `single a a 1 ⊗ₖ single b b 1` with weights its diagonal entries
   (`recordDiagonal_eq_sum_single_kronecker`).  With the record law the
   weights are nonnegative reals summing to one, so the state is a finite
   convex combination of product states of the two slots
   (`IsSlotSeparable`, `productDiagonal_separable`).  Separability is
   stated relative to the PR-44 split and to no other factorization.
2. **Slot-local bounded settings.**  `IsUnitInterval A` is the Loewner
   condition `-1 ≤ A ≤ 1` for a Hermitian `A`; every dichotomic
   observable `1 - 2P` of a projection event satisfies it
   (`isUnitInterval_of_isEvent`), and the diagonal entries of such an
   `A` are reals in `[-1, 1]` (`IsUnitInterval.diag_re_abs_le_one`).
   The slot-local CHSH operator `slotChsh A₀ A₁ B₀ B₁` is the CHSH
   combination of `slotLeft β Aᵢ` and `slotRight α Bⱼ`; the witness
   operator `chshProd` of `EventAlgebra.TsirelsonSaturation` is an
   instance (`chshProd_eq_slotChsh`).
3. **The slot-local bound.**  For every record-law state on `α × β` and
   every four unit-interval settings, two in each slot, the trace of the
   state against the slot-local CHSH operator is real with modulus at
   most `2` (`productDiagonal_slotLocal_chsh_le_two`,
   `productDiagonal_slotLocal_chsh_norm_le_two`).  No commutation
   hypothesis inside a slot is assumed: `A₀ A₁ ≠ A₁ A₀` is allowed, and
   the committed Pauli pair is such a pair
   (`pauliZ_pauliX_not_commute`, `within_wing_noncommutation_bounded`).
   On the committed pair carriers `Fin 13 × Fin 13` and `Fin 13 × Fin 14`
   every reachable state is separable across the split and obeys the
   bound against all slot-local unit-interval settings
   (`reachable_slotLocal_chsh_le_two`); this strengthens the
   jointly-diagonal clause of the source-reachability receipt.
4. **The counter-witness.**  Pairwise cross-wing commutation is not a
   sufficient hypothesis.  On `Fin 2 × Fin 2` the product-diagonal pure
   state `|00⟩⟨00|` (`ket00`), conjugated settings
   `Uᴴ (slotLeft Aᵢ) U`, `Uᴴ (slotRight Bⱼ) U` with `U = CNOT · (H ⊗ 1)`
   (`bellUnitary`) and the committed Bell settings, give four
   unit-interval observables that commute pairwise across the wings and
   whose CHSH combination has trace exactly `2√2` on `ket00`
   (`crossWing_commutation_not_sufficient`).  The conjugated left
   setting `conjA₁` is in neither slot (`conjA₁_not_slotLeft`,
   `conjA₁_not_slotRight`).  The hypothesis that carries the bound in
   clause 3 is slot-locality relative to the record split; commutation
   alone does not carry it.
5. **Residual.**  For record-law states on the PR-44 split, the only
   unit-interval readouts whose CHSH value can exceed `2` are readouts
   that are not slot-local with respect to the split
   (`slotLocal_excess_requires_nonRecordLaw`); on the committed pair
   carriers such an excess certifies non-reachability
   (`slotLocal_excess_not_reachable`).  Producing a slot-local violation
   requires a state that is not record-diagonal, the route delimited by
   `bellState_target_requires_nondiagonal_production`.
6. **Composed receipt.**  `productSplitSeparability_receipt` conjoins
   the clauses from a typed antecedent bundle `ProductSplitAntecedent`;
   `bellSettingsAntecedent` and `countedAntecedent` are inhabitants.

## Hypotheses, falsifier, nonclaims

Hypotheses: the PR-44 split (the two slot embeddings on the product
carrier), the PR-02 trace-pairing readout of states against observables,
and the record law of the state.  The reachable-class clauses consume
the closure theorem of the source-reachability module and its
hypotheses.

Falsifier: a record-law state on some `α × β` and four unit-interval
settings, two in each slot, with CHSH trace of modulus above `2`; or a
failure of any receipt clause.  The counter-witness clause is falsified
if its conjugated settings fail cross-wing commutation, fail the
unit-interval condition, or give a trace other than `2√2`.

Nonclaims: no source entangled preparation, no two-wing instrument, no
physical region, spacelike separation, or laboratory channel is produced
(PR-52 is open; OL-C3 realization clause).  The split is the declared
PR-44 datum.  The theorem is a separability bound for one declared
tensor split; it is not a hidden-variable theorem and not an
all-architectures no-go.  No statement about PR-04 or complex phase
structure is made.

## Tagging convention

As in `EventAlgebra.Basic`: **algebra-only** statements consume no trace
pairing; **trace-dependent** statements pass through `Tr(ρ M)`.
-/

namespace EventAlgebra.ProductSplitSeparability

open Matrix
open Kronecker
open OPH.QFT
open EventAlgebra.SourceReachability
open EventAlgebra.TsirelsonSaturation
open scoped ComplexOrder

noncomputable section

variable {α β : Type*} [Fintype α] [Fintype β] [DecidableEq α] [DecidableEq β]
variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-! ## Traces against record-diagonal states -/

omit [DecidableEq ι] in
/-- **Trace-dependent.**  Against a record-diagonal matrix the trace
pairing reads only the diagonal of the observable. -/
theorem trace_mul_of_recordDiagonal {ρ M : Matrix ι ι ℂ}
    (h : IsRecordDiagonal ρ) :
    (ρ * M).trace = ∑ i, ρ i i * M i i := by
  unfold Matrix.trace Matrix.diag
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Matrix.mul_apply]
  exact Finset.sum_eq_single i
    (fun j _ hji => by rw [h i j (Ne.symm hji), zero_mul])
    (fun hi => absurd (Finset.mem_univ i) hi)

/-- **Algebra-only.**  The product of a left-slot and a right-slot
element is the Kronecker product of the factors. -/
theorem slotLeft_mul_slotRight (A : Matrix α α ℂ) (B : Matrix β β ℂ) :
    slotLeft β A * slotRight α B = A ⊗ₖ B := by
  show (A ⊗ₖ (1 : Matrix β β ℂ)) * ((1 : Matrix α α ℂ) ⊗ₖ B) = A ⊗ₖ B
  rw [← Matrix.mul_kronecker_mul, mul_one, one_mul]

omit [DecidableEq α] [DecidableEq β] in
/-- **Trace-dependent.**  Against a record-diagonal state on the split,
a product observable is read through the diagonals of its factors. -/
theorem trace_mul_kronecker_of_recordDiagonal
    {ρ : Matrix (α × β) (α × β) ℂ} (h : IsRecordDiagonal ρ)
    (A : Matrix α α ℂ) (B : Matrix β β ℂ) :
    (ρ * (A ⊗ₖ B)).trace = ∑ p : α × β, ρ p p * (A p.1 p.1 * B p.2 p.2) := by
  rw [trace_mul_of_recordDiagonal h]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [Matrix.kroneckerMap_apply]

/-! ## Separability across the split -/

/-- **Algebra-only.**  The product-basis decomposition: a record-diagonal
matrix on the split is the sum of the product point matrices
`single a a 1 ⊗ₖ single b b 1` weighted by its diagonal entries. -/
theorem recordDiagonal_eq_sum_single_kronecker
    {ρ : Matrix (α × β) (α × β) ℂ} (h : IsRecordDiagonal ρ) :
    ρ = ∑ p : α × β, ρ p p •
      (Matrix.single p.1 p.1 (1 : ℂ) ⊗ₖ Matrix.single p.2 p.2 (1 : ℂ)) := by
  ext q r
  rw [Matrix.sum_apply]
  simp only [Matrix.smul_apply, Matrix.single_kronecker_single, one_mul,
    Prod.mk.eta, smul_eq_mul]
  by_cases hqr : q = r
  · subst hqr
    rw [Finset.sum_eq_single q
      (fun p _ hpq => by
        rw [Matrix.single_apply_of_ne p p (1 : ℂ) q q (fun hc => hpq hc.1),
          mul_zero])
      (fun hq => absurd (Finset.mem_univ q) hq),
      Matrix.single_apply_same, mul_one]
  · rw [h q r hqr]
    refine (Finset.sum_eq_zero fun p _ => ?_).symm
    rw [Matrix.single_apply_of_ne p p (1 : ℂ) q r
      (fun hc => hqr (hc.1.symm.trans hc.2)), mul_zero]

/-- **Separability relative to the PR-44 split.**  A matrix on `α × β`
is slot-separable when it is a finite convex combination of product
states `σᵢ ⊗ₖ τᵢ`, each factor positive semidefinite of trace one on its
slot.  The split is the declared PR-44 datum; separability is asserted
relative to it and to no other factorization. -/
def IsSlotSeparable (ρ : Matrix (α × β) (α × β) ℂ) : Prop :=
  ∃ (k : ℕ) (w : Fin k → ℝ) (σ : Fin k → Matrix α α ℂ)
    (τ : Fin k → Matrix β β ℂ),
    (∀ i, 0 ≤ w i) ∧ ∑ i, w i = 1 ∧
    (∀ i, (σ i).PosSemidef ∧ (σ i).trace = 1) ∧
    (∀ i, (τ i).PosSemidef ∧ (τ i).trace = 1) ∧
    ρ = ∑ i, (w i : ℂ) • (σ i ⊗ₖ τ i)

/-- **Algebra-only.**  The point matrix `single a a 1` is a state of its
slot: positive semidefinite with trace one. -/
theorem single_self_isState (a : α) :
    (Matrix.single a a (1 : ℂ)).PosSemidef ∧
      (Matrix.single a a (1 : ℂ)).trace = 1 := by
  refine ⟨?_, Matrix.trace_single_eq_same a 1⟩
  rw [← Matrix.diagonal_single]
  refine Matrix.posSemidef_diagonal_iff.mpr fun j => ?_
  rw [Pi.single_apply]
  split_ifs
  · exact zero_le_one
  · exact le_rfl

/-- **Product-diagonal states are separable.**  Every record-law state
on the split is slot-separable: the product point states with the
diagonal weights form the convex decomposition. -/
theorem productDiagonal_separable {ρ : Matrix (α × β) (α × β) ℂ}
    (h : IsRecordLaw ρ) : IsSlotSeparable ρ := by
  obtain ⟨w, hw0, hw1, rfl⟩ := h
  let e : α × β ≃ Fin (Fintype.card (α × β)) := Fintype.equivFin (α × β)
  have hd : IsRecordDiagonal (Matrix.diagonal fun p : α × β => (w p : ℂ)) :=
    fun i j hij => Matrix.diagonal_apply_ne _ hij
  refine ⟨Fintype.card (α × β), fun i => w (e.symm i),
    fun i => Matrix.single (e.symm i).1 (e.symm i).1 1,
    fun i => Matrix.single (e.symm i).2 (e.symm i).2 1,
    fun i => hw0 _, ?_, fun i => single_self_isState _,
    fun i => single_self_isState _, ?_⟩
  · rw [Equiv.sum_comp e.symm w]
    exact hw1
  · rw [recordDiagonal_eq_sum_single_kronecker hd]
    refine (Finset.sum_congr rfl fun p _ => ?_).trans
      (Equiv.sum_comp e.symm (fun p : α × β => (w p : ℂ) •
        (Matrix.single p.1 p.1 (1 : ℂ) ⊗ₖ Matrix.single p.2 p.2 (1 : ℂ)))).symm
    rw [Matrix.diagonal_apply_eq]

/-! ## Unit-interval observables -/

/-- **Algebra-only.**  The Loewner unit interval: `A` is Hermitian with
`-1 ≤ A ≤ 1`, stated as positive semidefiniteness of `1 - A` and
`1 + A`. -/
def IsUnitInterval (A : Matrix ι ι ℂ) : Prop :=
  A.IsHermitian ∧ (1 - A).PosSemidef ∧ (1 + A).PosSemidef

omit [DecidableEq ι] in
/-- **Algebra-only.**  A Hermitian idempotent is positive semidefinite. -/
theorem posSemidef_of_hermitian_idempotent {P : Matrix ι ι ℂ}
    (hP : P.IsHermitian) (hPP : P * P = P) : P.PosSemidef := by
  have h : Pᴴ * P = P := by rw [hP.eq, hPP]
  exact h ▸ Matrix.posSemidef_conjTranspose_mul_self P

/-- **Algebra-only.**  The dichotomic observable `1 - 2P` of a Hermitian
idempotent `P` lies in the unit interval. -/
theorem isUnitInterval_one_sub_two_smul {P : Matrix ι ι ℂ}
    (hP : P.IsHermitian) (hPP : P * P = P) :
    IsUnitInterval (1 - 2 • P) := by
  have hQ : ((1 : Matrix ι ι ℂ) - P).IsHermitian := by
    show ((1 : Matrix ι ι ℂ) - P)ᴴ = 1 - P
    rw [Matrix.conjTranspose_sub, Matrix.conjTranspose_one, hP.eq]
  have hQQ : ((1 : Matrix ι ι ℂ) - P) * (1 - P) = 1 - P := by
    rw [sub_mul, one_mul, mul_sub, mul_one, hPP]
    abel
  refine ⟨?_, ?_, ?_⟩
  · show ((1 : Matrix ι ι ℂ) - 2 • P)ᴴ = 1 - 2 • P
    rw [Matrix.conjTranspose_sub, Matrix.conjTranspose_one,
      Matrix.conjTranspose_nsmul, hP.eq]
  · have h : (1 : Matrix ι ι ℂ) - (1 - 2 • P) = P + P := by
      rw [two_smul]; abel
    rw [h]
    exact (posSemidef_of_hermitian_idempotent hP hPP).add
      (posSemidef_of_hermitian_idempotent hP hPP)
  · have h : (1 : Matrix ι ι ℂ) + (1 - 2 • P) = (1 - P) + (1 - P) := by
      rw [two_smul]; abel
    rw [h]
    exact (posSemidef_of_hermitian_idempotent hQ hQQ).add
      (posSemidef_of_hermitian_idempotent hQ hQQ)

/-- **Algebra-only.**  The dichotomic observable of a projection event is
a unit-interval observable. -/
theorem isUnitInterval_of_isEvent {n : ℕ} {P : Matrix (Fin n) (Fin n) ℂ}
    (hP : IsEvent P) : IsUnitInterval (1 - 2 • P) :=
  isUnitInterval_one_sub_two_smul hP.1 hP.2

omit [Fintype ι] in
/-- **Algebra-only.**  The identity is a unit-interval observable. -/
theorem isUnitInterval_one : IsUnitInterval (1 : Matrix ι ι ℂ) := by
  refine ⟨Matrix.isHermitian_one, ?_, ?_⟩
  · rw [sub_self]
    exact Matrix.PosSemidef.zero
  · exact Matrix.PosSemidef.one.add Matrix.PosSemidef.one

/-- **Algebra-only.**  The image of a dichotomic observable of a Hermitian
idempotent under a star algebra homomorphism is a unit-interval
observable; this covers the slot embeddings. -/
theorem isUnitInterval_map_one_sub_two_smul {κ : Type*} [Fintype κ]
    [DecidableEq κ] (φ : Matrix ι ι ℂ →⋆ₐ[ℂ] Matrix κ κ ℂ)
    {P : Matrix ι ι ℂ} (hP : P.IsHermitian) (hPP : P * P = P) :
    IsUnitInterval (φ (1 - 2 • P)) := by
  rw [map_sub, map_one, map_nsmul]
  apply isUnitInterval_one_sub_two_smul
  · show (φ P)ᴴ = φ P
    rw [← Matrix.star_eq_conjTranspose, ← map_star,
      Matrix.star_eq_conjTranspose, hP.eq]
  · rw [← map_mul, hPP]

omit [Fintype ι] in
/-- **Algebra-only.**  The diagonal entries of a unit-interval observable
are real. -/
theorem IsUnitInterval.diag_im_zero {A : Matrix ι ι ℂ} (h : IsUnitInterval A)
    (i : ι) : (A i i).im = 0 := by
  have h1 : star (A i i) = A i i := h.1.apply i i
  have h2 := congrArg Complex.im h1
  rw [Complex.star_def, Complex.conj_im] at h2
  linarith

omit [Fintype ι] in
/-- **Algebra-only.**  A diagonal entry of a unit-interval observable
equals the coercion of its real part. -/
theorem IsUnitInterval.diag_eq_re {A : Matrix ι ι ℂ} (h : IsUnitInterval A)
    (i : ι) : A i i = ((A i i).re : ℂ) := by
  apply Complex.ext
  · rw [Complex.ofReal_re]
  · rw [Complex.ofReal_im]
    exact h.diag_im_zero i

omit [Fintype ι] in
/-- **Algebra-only.**  The diagonal entries of a unit-interval observable
lie in `[-1, 1]`. -/
theorem IsUnitInterval.diag_re_abs_le_one {A : Matrix ι ι ℂ}
    (h : IsUnitInterval A) (i : ι) : |(A i i).re| ≤ 1 := by
  have h1 : (0 : ℂ) ≤ ((1 : Matrix ι ι ℂ) - A) i i := h.2.1.diag_nonneg
  have h2 : (0 : ℂ) ≤ ((1 : Matrix ι ι ℂ) + A) i i := h.2.2.diag_nonneg
  rw [Matrix.sub_apply, Matrix.one_apply_eq] at h1
  rw [Matrix.add_apply, Matrix.one_apply_eq] at h2
  have h1' := (Complex.nonneg_iff.mp h1).1
  have h2' := (Complex.nonneg_iff.mp h2).1
  rw [Complex.sub_re, Complex.one_re] at h1'
  rw [Complex.add_re, Complex.one_re] at h2'
  rw [abs_le]
  constructor <;> linarith

/-- **Algebra-only.**  Unitary conjugation preserves the unit interval. -/
theorem IsUnitInterval.conj {A U : Matrix ι ι ℂ} (h : IsUnitInterval A)
    (hU : Uᴴ * U = 1) : IsUnitInterval (Uᴴ * A * U) := by
  refine ⟨?_, ?_, ?_⟩
  · show (Uᴴ * A * U)ᴴ = Uᴴ * A * U
    rw [Matrix.conjTranspose_mul, Matrix.conjTranspose_mul,
      Matrix.conjTranspose_conjTranspose, h.1.eq, Matrix.mul_assoc]
  · have hs : (1 : Matrix ι ι ℂ) - Uᴴ * A * U = Uᴴ * (1 - A) * U := by
      rw [Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_one, hU]
    rw [hs]
    exact h.2.1.conjTranspose_mul_mul_same U
  · have hs : (1 : Matrix ι ι ℂ) + Uᴴ * A * U = Uᴴ * (1 + A) * U := by
      rw [Matrix.mul_add, Matrix.add_mul, Matrix.mul_one, hU]
    rw [hs]
    exact h.2.2.conjTranspose_mul_mul_same U

/-! ## The slot-local CHSH operator -/

/-- **Algebra-only.**  The CHSH combination of four operators on one
carrier, in the committed sign convention. -/
def chshCombination (X₀ X₁ Y₀ Y₁ : Matrix ι ι ℂ) : Matrix ι ι ℂ :=
  X₀ * Y₀ + X₀ * Y₁ + X₁ * Y₀ - X₁ * Y₁

/-- **Algebra-only.**  The slot-local CHSH operator on the PR-44 split:
the CHSH combination of the slot-lifted settings `Aᵢ ⊗ₖ 1` and
`1 ⊗ₖ Bⱼ`. -/
def slotChsh (A₀ A₁ : Matrix α α ℂ) (B₀ B₁ : Matrix β β ℂ) :
    Matrix (α × β) (α × β) ℂ :=
  chshCombination (slotLeft β A₀) (slotLeft β A₁) (slotRight α B₀)
    (slotRight α B₁)

/-- **Algebra-only.**  The witness CHSH operator of
`EventAlgebra.TsirelsonSaturation` is the slot-local CHSH operator of the
committed Bell settings. -/
theorem chshProd_eq_slotChsh :
    chshProd = slotChsh (1 - 2 • pA₀) (1 - 2 • pA₁) (1 - 2 • pB₀) (1 - 2 • pB₁) := by
  have hL : ∀ P : Matrix (Fin 2) (Fin 2) ℂ,
      slotLeft (Fin 2) ((1 : Matrix (Fin 2) (Fin 2) ℂ) - 2 • P) =
        1 - 2 • slotLeft (Fin 2) P := fun P => by
    rw [map_sub, map_one, map_nsmul]
  have hR : ∀ P : Matrix (Fin 2) (Fin 2) ℂ,
      slotRight (Fin 2) ((1 : Matrix (Fin 2) (Fin 2) ℂ) - 2 • P) =
        1 - 2 • slotRight (Fin 2) P := fun P => by
    rw [map_sub, map_one, map_nsmul]
  unfold slotChsh chshCombination
  rw [hL, hL, hR, hR]
  rfl

/-- **Trace-dependent.**  Against a record-diagonal state the slot-local
CHSH trace is the diagonal-weighted sum of the pointwise CHSH
combinations of the setting diagonals. -/
theorem trace_slotChsh_of_recordDiagonal {ρ : Matrix (α × β) (α × β) ℂ}
    (h : IsRecordDiagonal ρ) (A₀ A₁ : Matrix α α ℂ) (B₀ B₁ : Matrix β β ℂ) :
    (ρ * slotChsh A₀ A₁ B₀ B₁).trace =
      ∑ p : α × β, ρ p p *
        (A₀ p.1 p.1 * B₀ p.2 p.2 + A₀ p.1 p.1 * B₁ p.2 p.2 +
          A₁ p.1 p.1 * B₀ p.2 p.2 - A₁ p.1 p.1 * B₁ p.2 p.2) := by
  unfold slotChsh chshCombination
  simp only [slotLeft_mul_slotRight]
  rw [mul_sub, mul_add, mul_add, Matrix.trace_sub, Matrix.trace_add,
    Matrix.trace_add, trace_mul_kronecker_of_recordDiagonal h,
    trace_mul_kronecker_of_recordDiagonal h,
    trace_mul_kronecker_of_recordDiagonal h,
    trace_mul_kronecker_of_recordDiagonal h,
    ← Finset.sum_add_distrib, ← Finset.sum_add_distrib,
    ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun p _ => ?_
  ring

/-! ## The slot-local bound -/

/-- **Trace-dependent.**  For a record-law state and four unit-interval
slot-local settings, the CHSH trace is a real number of modulus at most
`2`. -/
theorem slotChsh_trace_real {ρ : Matrix (α × β) (α × β) ℂ}
    (hρ : IsRecordLaw ρ) {A₀ A₁ : Matrix α α ℂ} {B₀ B₁ : Matrix β β ℂ}
    (hA₀ : IsUnitInterval A₀) (hA₁ : IsUnitInterval A₁)
    (hB₀ : IsUnitInterval B₀) (hB₁ : IsUnitInterval B₁) :
    ∃ s : ℝ, |s| ≤ 2 ∧ (ρ * slotChsh A₀ A₁ B₀ B₁).trace = (s : ℂ) := by
  obtain ⟨w, hw0, hw1, rfl⟩ := hρ
  have hd : IsRecordDiagonal (Matrix.diagonal fun p : α × β => (w p : ℂ)) :=
    fun i j hij => Matrix.diagonal_apply_ne _ hij
  refine ⟨∑ p : α × β, w p *
    ((A₀ p.1 p.1).re * (B₀ p.2 p.2).re + (A₀ p.1 p.1).re * (B₁ p.2 p.2).re +
      (A₁ p.1 p.1).re * (B₀ p.2 p.2).re - (A₁ p.1 p.1).re * (B₁ p.2 p.2).re),
    ?_, ?_⟩
  · exact diagonal_chsh_le_two w (fun p => (A₀ p.1 p.1).re)
      (fun p => (A₁ p.1 p.1).re) (fun p => (B₀ p.2 p.2).re)
      (fun p => (B₁ p.2 p.2).re) hw0 hw1
      (fun p => hA₀.diag_re_abs_le_one p.1) (fun p => hA₁.diag_re_abs_le_one p.1)
      (fun p => hB₀.diag_re_abs_le_one p.2) (fun p => hB₁.diag_re_abs_le_one p.2)
  · rw [trace_slotChsh_of_recordDiagonal hd]
    push_cast
    refine Finset.sum_congr rfl fun p _ => ?_
    rw [Matrix.diagonal_apply_eq, ← hA₀.diag_eq_re, ← hA₁.diag_eq_re,
      ← hB₀.diag_eq_re, ← hB₁.diag_eq_re]

/-- **Trace-dependent.**  **The slot-local CHSH bound.**  Every
record-law state on the PR-44 split, read against every four
unit-interval settings, two in each slot, has a real CHSH trace of
modulus at most `2`.  No commutation inside a slot is assumed. -/
theorem productDiagonal_slotLocal_chsh_le_two {ρ : Matrix (α × β) (α × β) ℂ}
    (hρ : IsRecordLaw ρ) {A₀ A₁ : Matrix α α ℂ} {B₀ B₁ : Matrix β β ℂ}
    (hA₀ : IsUnitInterval A₀) (hA₁ : IsUnitInterval A₁)
    (hB₀ : IsUnitInterval B₀) (hB₁ : IsUnitInterval B₁) :
    ((ρ * slotChsh A₀ A₁ B₀ B₁).trace).im = 0 ∧
      |((ρ * slotChsh A₀ A₁ B₀ B₁).trace).re| ≤ 2 := by
  obtain ⟨s, hs, heq⟩ := slotChsh_trace_real hρ hA₀ hA₁ hB₀ hB₁
  rw [heq, Complex.ofReal_im, Complex.ofReal_re]
  exact ⟨rfl, hs⟩

/-- **Trace-dependent.**  Norm form of the slot-local CHSH bound. -/
theorem productDiagonal_slotLocal_chsh_norm_le_two
    {ρ : Matrix (α × β) (α × β) ℂ} (hρ : IsRecordLaw ρ)
    {A₀ A₁ : Matrix α α ℂ} {B₀ B₁ : Matrix β β ℂ}
    (hA₀ : IsUnitInterval A₀) (hA₁ : IsUnitInterval A₁)
    (hB₀ : IsUnitInterval B₀) (hB₁ : IsUnitInterval B₁) :
    ‖(ρ * slotChsh A₀ A₁ B₀ B₁).trace‖ ≤ 2 := by
  obtain ⟨s, hs, heq⟩ := slotChsh_trace_real hρ hA₀ hA₁ hB₀ hB₁
  rw [heq, Complex.norm_real, Real.norm_eq_abs]
  exact hs

/-- **Algebra-only.**  A record-diagonal state (positive semidefinite,
trace one) satisfies the record law; this is the state-form entry point
of the bound. -/
theorem isRecordLaw_of_recordDiagonal_state {ρ : Matrix ι ι ℂ}
    (hpsd : ρ.PosSemidef) (htr : ρ.trace = 1) (hd : IsRecordDiagonal ρ) :
    IsRecordLaw ρ := by
  refine ⟨fun i => (ρ i i).re,
    fun i => (Complex.nonneg_iff.mp hpsd.diag_nonneg).1, ?_, ?_⟩
  · have h := congrArg Complex.re htr
    simp only [Matrix.trace, Matrix.diag, Complex.re_sum, Complex.one_re] at h
    exact h
  · ext i j
    by_cases hij : i = j
    · subst hij
      rw [Matrix.diagonal_apply_eq]
      apply Complex.ext
      · rw [Complex.ofReal_re]
      · rw [Complex.ofReal_im]
        exact (Complex.nonneg_iff.mp hpsd.diag_nonneg).2.symm
    · rw [Matrix.diagonal_apply_ne _ hij]
      exact hd i j hij

/-- **Trace-dependent.**  State form of the slot-local bound: positive
semidefinite, trace one, record-diagonal. -/
theorem recordDiagonal_state_slotLocal_chsh_norm_le_two
    {ρ : Matrix (α × β) (α × β) ℂ} (hpsd : ρ.PosSemidef) (htr : ρ.trace = 1)
    (hd : IsRecordDiagonal ρ)
    {A₀ A₁ : Matrix α α ℂ} {B₀ B₁ : Matrix β β ℂ}
    (hA₀ : IsUnitInterval A₀) (hA₁ : IsUnitInterval A₁)
    (hB₀ : IsUnitInterval B₀) (hB₁ : IsUnitInterval B₁) :
    ‖(ρ * slotChsh A₀ A₁ B₀ B₁).trace‖ ≤ 2 :=
  productDiagonal_slotLocal_chsh_norm_le_two
    (isRecordLaw_of_recordDiagonal_state hpsd htr hd) hA₀ hA₁ hB₀ hB₁

/-! ## Noncommutation inside a wing -/

/-- **Algebra-only.**  The committed left-wing settings `Z` and `X` do not
commute. -/
theorem pauliZ_pauliX_not_commute :
    Robertson.pauliZ * Robertson.pauliX ≠ Robertson.pauliX * Robertson.pauliZ := by
  intro h
  have h01 := congrFun (congrFun h 0) 1
  norm_num [Robertson.pauliZ, Robertson.pauliX, Matrix.mul_apply,
    Fin.sum_univ_two] at h01

/-- **Trace-dependent.**  Noncommutation inside a wing does not lift a
record-law state above the bound: the committed Bell settings, with the
noncommuting left pair `Z`, `X`, read any record-law state on
`Fin 2 × Fin 2` to a CHSH trace of modulus at most `2`. -/
theorem within_wing_noncommutation_bounded
    {ρ : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ} (hρ : IsRecordLaw ρ) :
    Robertson.pauliZ * Robertson.pauliX ≠ Robertson.pauliX * Robertson.pauliZ ∧
      ‖(ρ * chshProd).trace‖ ≤ 2 := by
  refine ⟨pauliZ_pauliX_not_commute, ?_⟩
  rw [chshProd_eq_slotChsh]
  exact productDiagonal_slotLocal_chsh_norm_le_two hρ
    (isUnitInterval_of_isEvent pA₀_isEvent) (isUnitInterval_of_isEvent pA₁_isEvent)
    (isUnitInterval_of_isEvent pB₀_isEvent) (isUnitInterval_of_isEvent pB₁_isEvent)

/-! ## The committed pair carriers -/

/-- **Trace-dependent.**  On the 86/247 pair carrier every reachable state
obeys the slot-local bound against all unit-interval settings. -/
theorem reachable247_slotLocal_chsh_le_two
    {M : Matrix PairIndex247 PairIndex247 ℂ} (h : Reachable .pair247 M)
    {A₀ A₁ B₀ B₁ : Matrix (Fin 13) (Fin 13) ℂ}
    (hA₀ : IsUnitInterval A₀) (hA₁ : IsUnitInterval A₁)
    (hB₀ : IsUnitInterval B₀) (hB₁ : IsUnitInterval B₁) :
    ‖(M * slotChsh A₀ A₁ B₀ B₁).trace‖ ≤ 2 :=
  productDiagonal_slotLocal_chsh_norm_le_two (reachable_isRecordLaw h)
    hA₀ hA₁ hB₀ hB₁

/-- **Trace-dependent.**  On the 86/88 pair carrier every reachable state
obeys the slot-local bound against all unit-interval settings. -/
theorem reachable88_slotLocal_chsh_le_two
    {M : Matrix PairIndex PairIndex ℂ} (h : Reachable .pair88 M)
    {A₀ A₁ : Matrix (Fin 13) (Fin 13) ℂ} {B₀ B₁ : Matrix (Fin 14) (Fin 14) ℂ}
    (hA₀ : IsUnitInterval A₀) (hA₁ : IsUnitInterval A₁)
    (hB₀ : IsUnitInterval B₀) (hB₁ : IsUnitInterval B₁) :
    ‖(M * slotChsh A₀ A₁ B₀ B₁).trace‖ ≤ 2 :=
  productDiagonal_slotLocal_chsh_norm_le_two (reachable_isRecordLaw h)
    hA₀ hA₁ hB₀ hB₁

/-- **Trace-dependent.**  **The committed-carrier strengthening.**  Every
reachable state on the two committed pair carriers is separable across
the PR-44 split and obeys the CHSH bound `2` against all slot-local
unit-interval settings, jointly record-diagonal or not. -/
theorem reachable_slotLocal_chsh_le_two :
    (∀ M : Matrix PairIndex247 PairIndex247 ℂ, Reachable .pair247 M →
      IsSlotSeparable M ∧
      ∀ A₀ A₁ B₀ B₁ : Matrix (Fin 13) (Fin 13) ℂ,
        IsUnitInterval A₀ → IsUnitInterval A₁ →
        IsUnitInterval B₀ → IsUnitInterval B₁ →
        ‖(M * slotChsh A₀ A₁ B₀ B₁).trace‖ ≤ 2) ∧
    (∀ M : Matrix PairIndex PairIndex ℂ, Reachable .pair88 M →
      IsSlotSeparable M ∧
      ∀ (A₀ A₁ : Matrix (Fin 13) (Fin 13) ℂ) (B₀ B₁ : Matrix (Fin 14) (Fin 14) ℂ),
        IsUnitInterval A₀ → IsUnitInterval A₁ →
        IsUnitInterval B₀ → IsUnitInterval B₁ →
        ‖(M * slotChsh A₀ A₁ B₀ B₁).trace‖ ≤ 2) :=
  ⟨fun _ h => ⟨productDiagonal_separable (reachable_isRecordLaw h),
      fun _ _ _ _ hA₀ hA₁ hB₀ hB₁ =>
        reachable247_slotLocal_chsh_le_two h hA₀ hA₁ hB₀ hB₁⟩,
    fun _ h => ⟨productDiagonal_separable (reachable_isRecordLaw h),
      fun _ _ _ _ hA₀ hA₁ hB₀ hB₁ =>
        reachable88_slotLocal_chsh_le_two h hA₀ hA₁ hB₀ hB₁⟩⟩

/-- **Trace-dependent.**  The committed counted state obeys the
slot-local bound against all unit-interval settings; this strengthens
`counted_state_diagonal_chsh_le_two` from jointly diagonal readouts to
all slot-local readouts. -/
theorem counted_state_slotLocal_chsh_le_two
    {A₀ A₁ B₀ B₁ : Matrix (Fin 13) (Fin 13) ℂ}
    (hA₀ : IsUnitInterval A₀) (hA₁ : IsUnitInterval A₁)
    (hB₀ : IsUnitInterval B₀) (hB₁ : IsUnitInterval B₁) :
    ‖(correlationState * slotChsh A₀ A₁ B₀ B₁).trace‖ ≤ 2 :=
  reachable247_slotLocal_chsh_le_two correlationState_reachable hA₀ hA₁ hB₀ hB₁

/-! ## The counter-witness: cross-wing commutation is not sufficient -/

/-- The product-diagonal pure state `|00⟩⟨00|` on the two-qubit split. -/
def ket00 : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  Matrix.single ((0 : Fin 2), (0 : Fin 2)) ((0 : Fin 2), (0 : Fin 2)) 1

/-- **Algebra-only.**  `|00⟩⟨00|` satisfies the record law. -/
theorem ket00_isRecordLaw : IsRecordLaw ket00 := by
  refine ⟨fun p => if p = ((0 : Fin 2), (0 : Fin 2)) then 1 else 0,
    fun p => ?_, ?_, ?_⟩
  · dsimp only
    split_ifs <;> norm_num
  · simp
  · have hfun : (fun p : Fin 2 × Fin 2 =>
        (((if p = ((0 : Fin 2), (0 : Fin 2)) then (1 : ℝ) else 0) : ℝ) : ℂ)) =
        Pi.single ((0 : Fin 2), (0 : Fin 2)) (1 : ℂ) := by
      funext p
      rw [Pi.single_apply]
      split_ifs <;> simp
    rw [hfun, Matrix.diagonal_single]
    rfl

/-- The integer core of the basis change `CNOT · (H ⊗ 1)`: the entry at
`((a, b), (c, d))` is `(-1)^{ac}` when `b = d ⊕ a` and `0` otherwise. -/
def bellCore : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  Matrix.of fun p q =>
    if p.1 = 0 then (if p.2 = q.2 then 1 else 0)
    else (if p.2 = q.2 then 0 else (if q.1 = 0 then 1 else -1))

/-- The basis change `U = CNOT · (H ⊗ 1)`, carrying `|00⟩` to the Bell
vector `(|00⟩ + |11⟩)/√2`. -/
def bellUnitary : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  invSqrt2 • bellCore

theorem invSqrt2_mul_self : invSqrt2 * invSqrt2 = (2 : ℂ)⁻¹ := by
  unfold invSqrt2
  rw [← Complex.ofReal_mul, ← mul_inv,
    Real.mul_self_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
  push_cast
  rfl

theorem star_invSqrt2_eq : star invSqrt2 = invSqrt2 := by
  unfold invSqrt2
  rw [RCLike.star_def, Complex.conj_ofReal]

theorem bellCore_conjTranspose_mul :
    bellCoreᴴ * bellCore = (2 : ℂ) • (1 : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) := by
  ext ⟨a, b⟩ ⟨c, d⟩
  fin_cases a <;> fin_cases b <;> fin_cases c <;> fin_cases d <;>
    simp only [Matrix.mul_apply, Fintype.sum_prod_type, Fin.sum_univ_two] <;>
    norm_num [bellCore, Matrix.conjTranspose_apply]

theorem bellCore_mul_conjTranspose :
    bellCore * bellCoreᴴ = (2 : ℂ) • (1 : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) := by
  ext ⟨a, b⟩ ⟨c, d⟩
  fin_cases a <;> fin_cases b <;> fin_cases c <;> fin_cases d <;>
    simp only [Matrix.mul_apply, Fintype.sum_prod_type, Fin.sum_univ_two] <;>
    norm_num [bellCore, Matrix.conjTranspose_apply]

theorem bellCore_ket00 :
    bellCore * ket00 * bellCoreᴴ = OPH.Dynamics.maxEntangled 2 := by
  ext ⟨a, b⟩ ⟨c, d⟩
  fin_cases a <;> fin_cases b <;> fin_cases c <;> fin_cases d <;>
    simp only [Matrix.mul_apply, Fintype.sum_prod_type, Fin.sum_univ_two] <;>
    norm_num [bellCore, ket00, OPH.Dynamics.maxEntangled,
      Matrix.conjTranspose_apply, Matrix.vecMulVec_apply, Matrix.single_apply]

/-- **Algebra-only.**  `Uᴴ U = 1`. -/
theorem bellUnitary_conjTranspose_mul : bellUnitaryᴴ * bellUnitary = 1 := by
  unfold bellUnitary
  rw [Matrix.conjTranspose_smul, star_invSqrt2_eq, Matrix.smul_mul,
    Matrix.mul_smul, smul_smul, bellCore_conjTranspose_mul, smul_smul,
    invSqrt2_mul_self, inv_mul_cancel₀ two_ne_zero, one_smul]

/-- **Algebra-only.**  `U Uᴴ = 1`. -/
theorem bellUnitary_mul_conjTranspose : bellUnitary * bellUnitaryᴴ = 1 := by
  unfold bellUnitary
  rw [Matrix.conjTranspose_smul, star_invSqrt2_eq, Matrix.smul_mul,
    Matrix.mul_smul, smul_smul, bellCore_mul_conjTranspose, smul_smul,
    invSqrt2_mul_self, inv_mul_cancel₀ two_ne_zero, one_smul]

/-- **Algebra-only.**  `U |00⟩⟨00| Uᴴ` is the declared Bell state. -/
theorem bellUnitary_ket00 : bellUnitary * ket00 * bellUnitaryᴴ = bellStateProd := by
  unfold bellUnitary
  rw [Matrix.conjTranspose_smul, star_invSqrt2_eq, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.smul_mul, smul_smul, bellCore_ket00,
    invSqrt2_mul_self]
  rfl

/-- **Algebra-only.**  Conjugates multiply through the sandwich when
`U Uᴴ = 1`. -/
theorem conj_mul_conj {U X Y : Matrix ι ι ℂ} (hU : U * Uᴴ = 1) :
    (Uᴴ * X * U) * (Uᴴ * Y * U) = Uᴴ * (X * Y) * U := by
  calc (Uᴴ * X * U) * (Uᴴ * Y * U) = Uᴴ * X * (U * Uᴴ) * Y * U := by
        simp only [Matrix.mul_assoc]
    _ = Uᴴ * (X * Y) * U := by
        rw [hU, Matrix.mul_one]
        simp only [Matrix.mul_assoc]

/-- **Algebra-only.**  The CHSH combination of conjugated operators is the
conjugated CHSH combination. -/
theorem chshCombination_conj {U X₀ X₁ Y₀ Y₁ : Matrix ι ι ℂ} (hU : U * Uᴴ = 1) :
    chshCombination (Uᴴ * X₀ * U) (Uᴴ * X₁ * U) (Uᴴ * Y₀ * U) (Uᴴ * Y₁ * U) =
      Uᴴ * chshCombination X₀ X₁ Y₀ Y₁ * U := by
  unfold chshCombination
  rw [conj_mul_conj hU, conj_mul_conj hU, conj_mul_conj hU, conj_mul_conj hU]
  simp only [Matrix.mul_add, Matrix.add_mul, Matrix.mul_sub, Matrix.sub_mul]

/-- The conjugated left setting from `Z`. -/
def conjA₀ : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  bellUnitaryᴴ * slotLeft (Fin 2) (1 - 2 • pA₀) * bellUnitary

/-- The conjugated left setting from `X`. -/
def conjA₁ : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  bellUnitaryᴴ * slotLeft (Fin 2) (1 - 2 • pA₁) * bellUnitary

/-- The conjugated right setting from `(Z + X)/√2`. -/
def conjB₀ : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  bellUnitaryᴴ * slotRight (Fin 2) (1 - 2 • pB₀) * bellUnitary

/-- The conjugated right setting from `(Z - X)/√2`. -/
def conjB₁ : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  bellUnitaryᴴ * slotRight (Fin 2) (1 - 2 • pB₁) * bellUnitary

/-- **Algebra-only.**  Each conjugated left setting commutes with each
conjugated right setting, by the slot-disjointness theorem
`OPH.QFT.slot_commute` transported through the unitary. -/
theorem conj_cross_commute (X Y : Matrix (Fin 2) (Fin 2) ℂ) :
    (bellUnitaryᴴ * slotLeft (Fin 2) X * bellUnitary) *
        (bellUnitaryᴴ * slotRight (Fin 2) Y * bellUnitary) =
      (bellUnitaryᴴ * slotRight (Fin 2) Y * bellUnitary) *
        (bellUnitaryᴴ * slotLeft (Fin 2) X * bellUnitary) := by
  have hc : slotLeft (Fin 2) X * slotRight (Fin 2) Y =
      slotRight (Fin 2) Y * slotLeft (Fin 2) X :=
    (slot_commute ⟨X, rfl⟩ ⟨Y, rfl⟩).eq
  rw [conj_mul_conj bellUnitary_mul_conjTranspose,
    conj_mul_conj bellUnitary_mul_conjTranspose, hc]

theorem conjA₀_conjB₀_comm : conjA₀ * conjB₀ = conjB₀ * conjA₀ := conj_cross_commute _ _
theorem conjA₀_conjB₁_comm : conjA₀ * conjB₁ = conjB₁ * conjA₀ := conj_cross_commute _ _
theorem conjA₁_conjB₀_comm : conjA₁ * conjB₀ = conjB₀ * conjA₁ := conj_cross_commute _ _
theorem conjA₁_conjB₁_comm : conjA₁ * conjB₁ = conjB₁ * conjA₁ := conj_cross_commute _ _

theorem conjA₀_isUnitInterval : IsUnitInterval conjA₀ :=
  (isUnitInterval_map_one_sub_two_smul _ pA₀_isEvent.1 pA₀_isEvent.2).conj
    bellUnitary_conjTranspose_mul

theorem conjA₁_isUnitInterval : IsUnitInterval conjA₁ :=
  (isUnitInterval_map_one_sub_two_smul _ pA₁_isEvent.1 pA₁_isEvent.2).conj
    bellUnitary_conjTranspose_mul

theorem conjB₀_isUnitInterval : IsUnitInterval conjB₀ :=
  (isUnitInterval_map_one_sub_two_smul _ pB₀_isEvent.1 pB₀_isEvent.2).conj
    bellUnitary_conjTranspose_mul

theorem conjB₁_isUnitInterval : IsUnitInterval conjB₁ :=
  (isUnitInterval_map_one_sub_two_smul _ pB₁_isEvent.1 pB₁_isEvent.2).conj
    bellUnitary_conjTranspose_mul

/-- **Trace-dependent.**  On the product-diagonal state `|00⟩⟨00|` the
conjugated settings reach the CHSH value `2√2` exactly, by trace
cyclicity and the committed pairing `bellStateProd_chshProd_trace`. -/
theorem ket00_conj_chsh_trace :
    (ket00 * chshCombination conjA₀ conjA₁ conjB₀ conjB₁).trace =
      ((2 * Real.sqrt 2 : ℝ) : ℂ) := by
  unfold conjA₀ conjA₁ conjB₀ conjB₁
  rw [chshCombination_conj bellUnitary_mul_conjTranspose,
    show chshCombination (slotLeft (Fin 2) (1 - 2 • pA₀))
        (slotLeft (Fin 2) (1 - 2 • pA₁)) (slotRight (Fin 2) (1 - 2 • pB₀))
        (slotRight (Fin 2) (1 - 2 • pB₁)) = chshProd from chshProd_eq_slotChsh.symm]
  simp only [← Matrix.mul_assoc]
  rw [Matrix.trace_mul_comm]
  simp only [← Matrix.mul_assoc]
  rw [bellUnitary_ket00]
  exact bellStateProd_chshProd_trace

/-- **Algebra-only.**  A left-slot element vanishes off the right
diagonal. -/
theorem slotLeft_apply_of_ne {A : Matrix α α ℂ} {a a' : α} {b b' : β}
    (h : b ≠ b') : slotLeft β A (a, b) (a', b') = 0 := by
  show (A ⊗ₖ (1 : Matrix β β ℂ)) (a, b) (a', b') = 0
  rw [Matrix.kronecker_apply, Matrix.one_apply_ne h, mul_zero]

/-- **Algebra-only.**  A right-slot element has left-diagonal blocks
independent of the left index. -/
theorem slotRight_apply_diag {B : Matrix β β ℂ} (a a' : α) (b b' : β) :
    slotRight α B (a, b) (a, b') = slotRight α B (a', b) (a', b') := by
  show ((1 : Matrix α α ℂ) ⊗ₖ B) (a, b) (a, b') =
    ((1 : Matrix α α ℂ) ⊗ₖ B) (a', b) (a', b')
  rw [Matrix.kronecker_apply, Matrix.kronecker_apply, Matrix.one_apply_eq,
    Matrix.one_apply_eq]

theorem conjA₁_entry_00_01 :
    conjA₁ ((0 : Fin 2), (0 : Fin 2)) ((0 : Fin 2), (1 : Fin 2)) = 1 := by
  unfold conjA₁ bellUnitary
  rw [pA₁_dichotomic]
  show (((invSqrt2 • bellCore)ᴴ * (Robertson.pauliX ⊗ₖ (1 : Matrix (Fin 2) (Fin 2) ℂ))) *
    (invSqrt2 • bellCore)) _ _ = 1
  rw [Matrix.conjTranspose_smul, star_invSqrt2_eq, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.smul_mul, smul_smul, invSqrt2_mul_self]
  norm_num [bellCore, Robertson.pauliX, Matrix.mul_apply, Fintype.sum_prod_type,
    Fin.sum_univ_two, Matrix.conjTranspose_apply, Matrix.kroneckerMap_apply,
    Matrix.one_apply]

theorem conjA₁_entry_10_11 :
    conjA₁ ((1 : Fin 2), (0 : Fin 2)) ((1 : Fin 2), (1 : Fin 2)) = -1 := by
  unfold conjA₁ bellUnitary
  rw [pA₁_dichotomic]
  show (((invSqrt2 • bellCore)ᴴ * (Robertson.pauliX ⊗ₖ (1 : Matrix (Fin 2) (Fin 2) ℂ))) *
    (invSqrt2 • bellCore)) _ _ = -1
  rw [Matrix.conjTranspose_smul, star_invSqrt2_eq, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.smul_mul, smul_smul, invSqrt2_mul_self]
  norm_num [bellCore, Robertson.pauliX, Matrix.mul_apply, Fintype.sum_prod_type,
    Fin.sum_univ_two, Matrix.conjTranspose_apply, Matrix.kroneckerMap_apply,
    Matrix.one_apply]

/-- **Algebra-only.**  The conjugated left setting `conjA₁` is not a
left-slot element. -/
theorem conjA₁_not_slotLeft (A : Matrix (Fin 2) (Fin 2) ℂ) :
    slotLeft (Fin 2) A ≠ conjA₁ := by
  intro hA
  have h := congrFun (congrFun hA ((0 : Fin 2), (0 : Fin 2))) ((0 : Fin 2), (1 : Fin 2))
  rw [slotLeft_apply_of_ne (by decide), conjA₁_entry_00_01] at h
  exact zero_ne_one h

/-- **Algebra-only.**  The conjugated left setting `conjA₁` is not a
right-slot element. -/
theorem conjA₁_not_slotRight (B : Matrix (Fin 2) (Fin 2) ℂ) :
    slotRight (Fin 2) B ≠ conjA₁ := by
  intro hB
  have h0 := congrFun (congrFun hB ((0 : Fin 2), (0 : Fin 2))) ((0 : Fin 2), (1 : Fin 2))
  have h1 := congrFun (congrFun hB ((1 : Fin 2), (0 : Fin 2))) ((1 : Fin 2), (1 : Fin 2))
  rw [conjA₁_entry_00_01] at h0
  rw [conjA₁_entry_10_11, slotRight_apply_diag 1 0 0 1, h0] at h1
  norm_num at h1

/-- **The counter-witness.**  There exist a record-law (product-diagonal)
state on the two-qubit split and four unit-interval observables that
commute pairwise across the wings, with CHSH trace exactly `2√2`; one of
the observables is in neither slot.  The hypothesis that carries the
slot-local bound is slot-locality relative to the record split;
cross-wing commutation alone does not carry it. -/
theorem crossWing_commutation_not_sufficient :
    ∃ (ρ X₀ X₁ Y₀ Y₁ : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ),
      IsRecordLaw ρ ∧
      IsUnitInterval X₀ ∧ IsUnitInterval X₁ ∧
      IsUnitInterval Y₀ ∧ IsUnitInterval Y₁ ∧
      X₀ * Y₀ = Y₀ * X₀ ∧ X₀ * Y₁ = Y₁ * X₀ ∧
      X₁ * Y₀ = Y₀ * X₁ ∧ X₁ * Y₁ = Y₁ * X₁ ∧
      (ρ * chshCombination X₀ X₁ Y₀ Y₁).trace = ((2 * Real.sqrt 2 : ℝ) : ℂ) ∧
      (∀ A : Matrix (Fin 2) (Fin 2) ℂ, slotLeft (Fin 2) A ≠ X₁) ∧
      (∀ B : Matrix (Fin 2) (Fin 2) ℂ, slotRight (Fin 2) B ≠ X₁) :=
  ⟨ket00, conjA₀, conjA₁, conjB₀, conjB₁, ket00_isRecordLaw,
    conjA₀_isUnitInterval, conjA₁_isUnitInterval, conjB₀_isUnitInterval,
    conjB₁_isUnitInterval, conjA₀_conjB₀_comm, conjA₀_conjB₁_comm,
    conjA₁_conjB₀_comm, conjA₁_conjB₁_comm, ket00_conj_chsh_trace,
    conjA₁_not_slotLeft, conjA₁_not_slotRight⟩

/-! ## The residual -/

/-- **Trace-dependent.**  A slot-local CHSH excess above `2` certifies
that the matrix is not a record law: on the PR-44 split, the only
unit-interval readouts that can exceed `2` on a record-law state are
readouts that are not slot-local with respect to the split. -/
theorem slotLocal_excess_requires_nonRecordLaw {ρ : Matrix (α × β) (α × β) ℂ}
    {A₀ A₁ : Matrix α α ℂ} {B₀ B₁ : Matrix β β ℂ}
    (hA₀ : IsUnitInterval A₀) (hA₁ : IsUnitInterval A₁)
    (hB₀ : IsUnitInterval B₀) (hB₁ : IsUnitInterval B₁)
    (hex : 2 < ‖(ρ * slotChsh A₀ A₁ B₀ B₁).trace‖) : ¬ IsRecordLaw ρ :=
  fun h => absurd (productDiagonal_slotLocal_chsh_norm_le_two h hA₀ hA₁ hB₀ hB₁)
    (not_le.mpr hex)

/-- **Trace-dependent.**  For a state (positive semidefinite, trace one),
a slot-local excess certifies off-diagonal record coherence. -/
theorem slotLocal_excess_requires_nonDiagonal {ρ : Matrix (α × β) (α × β) ℂ}
    (hpsd : ρ.PosSemidef) (htr : ρ.trace = 1)
    {A₀ A₁ : Matrix α α ℂ} {B₀ B₁ : Matrix β β ℂ}
    (hA₀ : IsUnitInterval A₀) (hA₁ : IsUnitInterval A₁)
    (hB₀ : IsUnitInterval B₀) (hB₁ : IsUnitInterval B₁)
    (hex : 2 < ‖(ρ * slotChsh A₀ A₁ B₀ B₁).trace‖) : ¬ IsRecordDiagonal ρ :=
  fun hd => slotLocal_excess_requires_nonRecordLaw hA₀ hA₁ hB₀ hB₁ hex
    (isRecordLaw_of_recordDiagonal_state hpsd htr hd)

/-- **Trace-dependent.**  On the committed pair carriers a slot-local
excess certifies non-reachability: no state of the committed
source-reachable class supports it. -/
theorem slotLocal_excess_not_reachable :
    (∀ (M : Matrix PairIndex247 PairIndex247 ℂ)
        (A₀ A₁ B₀ B₁ : Matrix (Fin 13) (Fin 13) ℂ),
      IsUnitInterval A₀ → IsUnitInterval A₁ →
      IsUnitInterval B₀ → IsUnitInterval B₁ →
      2 < ‖(M * slotChsh A₀ A₁ B₀ B₁).trace‖ → ¬ Reachable .pair247 M) ∧
    (∀ (M : Matrix PairIndex PairIndex ℂ)
        (A₀ A₁ : Matrix (Fin 13) (Fin 13) ℂ) (B₀ B₁ : Matrix (Fin 14) (Fin 14) ℂ),
      IsUnitInterval A₀ → IsUnitInterval A₁ →
      IsUnitInterval B₀ → IsUnitInterval B₁ →
      2 < ‖(M * slotChsh A₀ A₁ B₀ B₁).trace‖ → ¬ Reachable .pair88 M) :=
  ⟨fun _ _ _ _ _ hA₀ hA₁ hB₀ hB₁ hex h =>
      slotLocal_excess_requires_nonRecordLaw hA₀ hA₁ hB₀ hB₁ hex
        (reachable_isRecordLaw h),
    fun _ _ _ _ _ hA₀ hA₁ hB₀ hB₁ hex h =>
      slotLocal_excess_requires_nonRecordLaw hA₀ hA₁ hB₀ hB₁ hex
        (reachable_isRecordLaw h)⟩

/-! ## The composed receipt -/

/-- **The typed antecedent bundle.**  One record-law state on the PR-44
split `α × β` and four unit-interval settings, two in each slot.  No
commutation inside a slot is required. -/
structure ProductSplitAntecedent (α β : Type*) [Fintype α] [Fintype β]
    [DecidableEq α] [DecidableEq β] where
  /-- The state on the split. -/
  state : Matrix (α × β) (α × β) ℂ
  /-- The state satisfies the record law. -/
  law : IsRecordLaw state
  /-- First left-slot setting. -/
  A₀ : Matrix α α ℂ
  /-- Second left-slot setting. -/
  A₁ : Matrix α α ℂ
  /-- First right-slot setting. -/
  B₀ : Matrix β β ℂ
  /-- Second right-slot setting. -/
  B₁ : Matrix β β ℂ
  /-- `A₀` lies in the unit interval. -/
  unitA₀ : IsUnitInterval A₀
  /-- `A₁` lies in the unit interval. -/
  unitA₁ : IsUnitInterval A₁
  /-- `B₀` lies in the unit interval. -/
  unitB₀ : IsUnitInterval B₀
  /-- `B₁` lies in the unit interval. -/
  unitB₁ : IsUnitInterval B₁

/-- **Nonvacuity.**  The two-qubit bundle: `|00⟩⟨00|` with the committed
Bell settings `Z`, `X`, `(Z ± X)/√2`. -/
def bellSettingsAntecedent : ProductSplitAntecedent (Fin 2) (Fin 2) where
  state := ket00
  law := ket00_isRecordLaw
  A₀ := 1 - 2 • pA₀
  A₁ := 1 - 2 • pA₁
  B₀ := 1 - 2 • pB₀
  B₁ := 1 - 2 • pB₁
  unitA₀ := isUnitInterval_of_isEvent pA₀_isEvent
  unitA₁ := isUnitInterval_of_isEvent pA₁_isEvent
  unitB₀ := isUnitInterval_of_isEvent pB₀_isEvent
  unitB₁ := isUnitInterval_of_isEvent pB₁_isEvent

/-- **Nonvacuity.**  The committed 86/247 bundle: the counted state with
the dichotomic observables of two committed observer-86 checkpoint
projectors on the left slot and the identity on the right slot. -/
def countedAntecedent : ProductSplitAntecedent (Fin 13) (Fin 13) where
  state := correlationState
  law := isRecordLaw_correlationState
  A₀ := 1 - 2 • chk86 0
  A₁ := 1 - 2 • chk86 1
  B₀ := 1
  B₁ := 1
  unitA₀ := isUnitInterval_one_sub_two_smul (chk86_conjTranspose 0) (chk86_idem 0)
  unitA₁ := isUnitInterval_one_sub_two_smul (chk86_conjTranspose 1) (chk86_idem 1)
  unitB₀ := isUnitInterval_one
  unitB₁ := isUnitInterval_one

/-- **The product-split separability receipt.**  For every antecedent
bundle on every PR-44 split, one conjunction:

1. **Separability**: the state is a finite convex combination of product
   states of the two slots;
2. **Reality**: the slot-local CHSH trace is real;
3. **Bound**: its real part has modulus at most `2`, and so does its
   norm;
4. **Counter-witness**: a record-law state and four pairwise
   cross-commuting unit-interval observables reach `2√2`, one of them in
   neither slot;
5. `2 < 2√2`.

Boundary, stated exactly: the split is the declared PR-44 datum; the
bound covers slot-local settings with respect to that split and says
nothing about readouts outside both slots; no source entangled
preparation, instrument, or physical region is produced (PR-52 is open;
OL-C3 realization clause). -/
theorem productSplitSeparability_receipt (D : ProductSplitAntecedent α β) :
    IsSlotSeparable D.state ∧
    ((D.state * slotChsh D.A₀ D.A₁ D.B₀ D.B₁).trace).im = 0 ∧
    |((D.state * slotChsh D.A₀ D.A₁ D.B₀ D.B₁).trace).re| ≤ 2 ∧
    ‖(D.state * slotChsh D.A₀ D.A₁ D.B₀ D.B₁).trace‖ ≤ 2 ∧
    (∃ (ρ X₀ X₁ Y₀ Y₁ : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ),
      IsRecordLaw ρ ∧
      IsUnitInterval X₀ ∧ IsUnitInterval X₁ ∧
      IsUnitInterval Y₀ ∧ IsUnitInterval Y₁ ∧
      X₀ * Y₀ = Y₀ * X₀ ∧ X₀ * Y₁ = Y₁ * X₀ ∧
      X₁ * Y₀ = Y₀ * X₁ ∧ X₁ * Y₁ = Y₁ * X₁ ∧
      (ρ * chshCombination X₀ X₁ Y₀ Y₁).trace = ((2 * Real.sqrt 2 : ℝ) : ℂ) ∧
      (∀ A : Matrix (Fin 2) (Fin 2) ℂ, slotLeft (Fin 2) A ≠ X₁) ∧
      (∀ B : Matrix (Fin 2) (Fin 2) ℂ, slotRight (Fin 2) B ≠ X₁)) ∧
    (2 : ℝ) < 2 * Real.sqrt 2 :=
  ⟨productDiagonal_separable D.law,
    (productDiagonal_slotLocal_chsh_le_two D.law D.unitA₀ D.unitA₁ D.unitB₀
      D.unitB₁).1,
    (productDiagonal_slotLocal_chsh_le_two D.law D.unitA₀ D.unitA₁ D.unitB₀
      D.unitB₁).2,
    productDiagonal_slotLocal_chsh_norm_le_two D.law D.unitA₀ D.unitA₁ D.unitB₀
      D.unitB₁,
    crossWing_commutation_not_sufficient,
    classical_bound_lt_saturation⟩

/-- **Nonvacuity of the receipt.**  The receipt is instantiated at the
two-qubit bundle and at the committed 86/247 bundle. -/
theorem productSplitSeparability_receipt_nonvacuous :
    Nonempty (ProductSplitAntecedent (Fin 2) (Fin 2)) ∧
    Nonempty (ProductSplitAntecedent (Fin 13) (Fin 13)) :=
  ⟨⟨bellSettingsAntecedent⟩, ⟨countedAntecedent⟩⟩

end

end EventAlgebra.ProductSplitSeparability

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
-- No native_decide is used anywhere in this module.
#print axioms EventAlgebra.ProductSplitSeparability.trace_mul_of_recordDiagonal
#print axioms EventAlgebra.ProductSplitSeparability.recordDiagonal_eq_sum_single_kronecker
#print axioms EventAlgebra.ProductSplitSeparability.productDiagonal_separable
#print axioms EventAlgebra.ProductSplitSeparability.isUnitInterval_of_isEvent
#print axioms EventAlgebra.ProductSplitSeparability.IsUnitInterval.diag_re_abs_le_one
#print axioms EventAlgebra.ProductSplitSeparability.IsUnitInterval.conj
#print axioms EventAlgebra.ProductSplitSeparability.chshProd_eq_slotChsh
#print axioms EventAlgebra.ProductSplitSeparability.trace_slotChsh_of_recordDiagonal
#print axioms EventAlgebra.ProductSplitSeparability.productDiagonal_slotLocal_chsh_le_two
#print axioms EventAlgebra.ProductSplitSeparability.productDiagonal_slotLocal_chsh_norm_le_two
#print axioms EventAlgebra.ProductSplitSeparability.recordDiagonal_state_slotLocal_chsh_norm_le_two
#print axioms EventAlgebra.ProductSplitSeparability.within_wing_noncommutation_bounded
#print axioms EventAlgebra.ProductSplitSeparability.reachable_slotLocal_chsh_le_two
#print axioms EventAlgebra.ProductSplitSeparability.counted_state_slotLocal_chsh_le_two
#print axioms EventAlgebra.ProductSplitSeparability.ket00_isRecordLaw
#print axioms EventAlgebra.ProductSplitSeparability.bellUnitary_conjTranspose_mul
#print axioms EventAlgebra.ProductSplitSeparability.bellUnitary_ket00
#print axioms EventAlgebra.ProductSplitSeparability.ket00_conj_chsh_trace
#print axioms EventAlgebra.ProductSplitSeparability.conjA₁_not_slotLeft
#print axioms EventAlgebra.ProductSplitSeparability.conjA₁_not_slotRight
#print axioms EventAlgebra.ProductSplitSeparability.crossWing_commutation_not_sufficient
#print axioms EventAlgebra.ProductSplitSeparability.slotLocal_excess_requires_nonRecordLaw
#print axioms EventAlgebra.ProductSplitSeparability.slotLocal_excess_not_reachable
#print axioms EventAlgebra.ProductSplitSeparability.productSplitSeparability_receipt
#print axioms EventAlgebra.ProductSplitSeparability.productSplitSeparability_receipt_nonvacuous
