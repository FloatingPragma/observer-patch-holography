import EventAlgebra.ExpectationBound
import EventAlgebra.SlotLocalitySurface
import EventAlgebra.Robertson
import Dynamics.ChoiCPTP
import QFT.SourceCorrelationCapstone

set_option autoImplicit false

/-!
# Tsirelson saturation on the committed slot interface, and the diagonal CHSH delimitation

Issue #730 OL-C3.  Two finite receipts on committed carriers.

1. **Saturation witness.**  A declared four-dimensional state — the
   normalized `OPH.Dynamics.maxEntangled 2`, the one entangled object of the
   committed corpus, carried to `Fin 4` along the fixed index equivalence
   `finProdFinEquiv` — together with four declared projection events (the
   spectral projections of the two Pauli directions on the left factor and
   of the two `π/4`-rotated directions on the right factor) attains the
   CHSH expectation value `2√2` exactly (`chsh_saturation`).  The four
   factor events inhabit the committed PR-44 interface
   `EventAlgebra.SlotLocalityInterface` at `da = db = 2` with the singleton
   identity Kraus family (`bellSlotInterface`); the committed norm bound of
   `EventAlgebra.QuantumSurface.slot_locality_receipts` and the committed
   state-expectation bound `matrix_state_tsirelson_bound_of_events` are
   cited, not re-proved.  The packaged statement is
   `tsirelson_saturation_receipt`: upper bound, exact attainment, and the
   identification of the attaining observable with the slot-lifted CHSH
   operator of the interface witness.  The Bell state and all four setting
   matrices are real; no Pauli-Y or genuinely complex phase effect enters
   this witness.

2. **Diagonal delimitation.**  For any convex weights and any four
   `[-1,1]`-valued diagonal observables on any finite carrier, the CHSH
   expectation has absolute value at most `2` (`diagonal_chsh_le_two`,
   `diagonal_state_chsh_le_two`), and in particular the committed counted
   state `OPH.QFT.correlationState` obeys the classical bound `2` against
   every record-diagonal CHSH readout
   (`counted_state_diagonal_chsh_le_two`).

What IS proved here: finite matrix algebra over `ℂ` on the carriers
`Fin 2`, `Fin 2 × Fin 2`, `Fin 4`, and `Fin 13 × Fin 13` — event
certificates, a state certificate (positive semidefinite, trace one), the
exact expectation value `2√2`, and the exact classical bound `2`.  All of
it is premise-free finite mathematics given the committed modules it
imports.

What is NOT proved (open premises, by register row):

* the saturating state and the four events are DECLARED data of this
  module; no theorem produces them from committed source runs — source
  entangled preparation, physical regions, spacelike separation, and
  instruments stay open (register row PR-52, and the OL-C3 realization
  clause of issue #730);
* the bipartite slot split consumed by the witness is the supplied PR-44
  datum; nothing here produces the split;
* the delimitation clause states only that one fixed diagonal state with four
  jointly record-diagonal bounded readouts cannot exceed `2`.  It does not
  cover arbitrary context-dependent classical experiments and it does not
  identify the missing resource with PR-04 or with genuinely complex phase
  structure.  This real Bell witness needs joint coherence/entanglement and
  suitable noncommuting settings; source production of those objects remains
  open.

## Tagging convention

As in `EventAlgebra.Basic`: **algebra-only** statements consume no trace
pairing; **trace-dependent** statements pass through `Tr(ρ M)`.
-/

namespace EventAlgebra.TsirelsonSaturation

open Matrix
open Kronecker
open OPH.QFT
open scoped ComplexOrder
open scoped Matrix.Norms.L2Operator

noncomputable section

/-! ## Scalar facts about `√2` -/

/-- The complex coercion of `(√2)⁻¹`, the single irrational scalar of the
witness. -/
def invSqrt2 : ℂ := ((Real.sqrt 2)⁻¹ : ℝ)

private theorem sqrt2_mul_self : Real.sqrt 2 * Real.sqrt 2 = 2 :=
  Real.mul_self_sqrt (by norm_num)

private theorem sqrt2_ne_zero : Real.sqrt 2 ≠ 0 :=
  ne_of_gt (Real.sqrt_pos.mpr (by norm_num))

private theorem invSqrt2_sq_mul_two : invSqrt2 * invSqrt2 * 2 = 1 := by
  have h : ((Real.sqrt 2)⁻¹ * (Real.sqrt 2)⁻¹ * 2 : ℝ) = 1 := by
    field_simp
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
  calc invSqrt2 * invSqrt2 * 2
      = (((Real.sqrt 2)⁻¹ * (Real.sqrt 2)⁻¹ * 2 : ℝ) : ℂ) := by
        unfold invSqrt2; push_cast; ring
    _ = 1 := by rw [h]; norm_num

private theorem invSqrt2_mul_two : invSqrt2 * 2 = ((Real.sqrt 2 : ℝ) : ℂ) := by
  have h : ((Real.sqrt 2)⁻¹ * 2 : ℝ) = Real.sqrt 2 := by
    rw [inv_mul_eq_div, div_eq_iff sqrt2_ne_zero]
    exact sqrt2_mul_self.symm
  calc invSqrt2 * 2 = (((Real.sqrt 2)⁻¹ * 2 : ℝ) : ℂ) := by
        unfold invSqrt2; push_cast; ring
    _ = ((Real.sqrt 2 : ℝ) : ℂ) := by rw [h]

/-! ## The two factor algebras: involutions and events

The left-party observables are the committed Pauli matrices
`EventAlgebra.Robertson.pauliZ` and `pauliX`; the right-party observables
are the two `π/4`-rotated combinations `(Z ± X)/√2`.  All four are
self-adjoint involutions, so `(1 - A)/2` is a projection event in each
case.  Everything in this section is **algebra-only**. -/

private theorem pauliZ_mul_self :
    Robertson.pauliZ * Robertson.pauliZ = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [Robertson.pauliZ, Matrix.mul_apply, Fin.sum_univ_two,
      Matrix.one_apply]

private theorem pauliX_mul_self :
    Robertson.pauliX * Robertson.pauliX = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [Robertson.pauliX, Matrix.mul_apply, Fin.sum_univ_two,
      Matrix.one_apply]

private theorem pauliZ_mul_pauliX_anti :
    Robertson.pauliZ * Robertson.pauliX +
      Robertson.pauliX * Robertson.pauliZ = 0 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [Robertson.pauliZ, Robertson.pauliX, Matrix.mul_apply,
      Fin.sum_univ_two, Matrix.add_apply, Matrix.zero_apply]

private theorem pauliZX_eq_neg :
    Robertson.pauliZ * Robertson.pauliX =
      -(Robertson.pauliX * Robertson.pauliZ) :=
  eq_neg_of_add_eq_zero_left pauliZ_mul_pauliX_anti

/-- **Algebra-only.** The first right-party observable `(Z + X)/√2`. -/
def bObs₀ : Matrix (Fin 2) (Fin 2) ℂ :=
  invSqrt2 • (Robertson.pauliZ + Robertson.pauliX)

/-- **Algebra-only.** The second right-party observable `(Z - X)/√2`. -/
def bObs₁ : Matrix (Fin 2) (Fin 2) ℂ :=
  invSqrt2 • (Robertson.pauliZ - Robertson.pauliX)

private theorem star_invSqrt2 : star invSqrt2 = invSqrt2 := by
  unfold invSqrt2
  rw [RCLike.star_def, Complex.conj_ofReal]

/-- **Algebra-only.** `(Z + X)/√2` is Hermitian. -/
theorem bObs₀_isHermitian : bObs₀.IsHermitian := by
  show bObs₀ᴴ = bObs₀
  unfold bObs₀
  rw [Matrix.conjTranspose_smul, star_invSqrt2, Matrix.conjTranspose_add,
    Robertson.pauliZ_isHermitian.eq, Robertson.pauliX_isHermitian.eq]

/-- **Algebra-only.** `(Z - X)/√2` is Hermitian. -/
theorem bObs₁_isHermitian : bObs₁.IsHermitian := by
  show bObs₁ᴴ = bObs₁
  unfold bObs₁
  rw [Matrix.conjTranspose_smul, star_invSqrt2, Matrix.conjTranspose_sub,
    Robertson.pauliZ_isHermitian.eq, Robertson.pauliX_isHermitian.eq]

/-- **Algebra-only.** `(Z + X)/√2` is an involution. -/
theorem bObs₀_mul_self : bObs₀ * bObs₀ = 1 := by
  unfold bObs₀
  have h : (Robertson.pauliZ + Robertson.pauliX) *
      (Robertson.pauliZ + Robertson.pauliX) =
        (2 : ℂ) • (1 : Matrix (Fin 2) (Fin 2) ℂ) := by
    rw [add_mul, mul_add, mul_add, pauliZ_mul_self, pauliX_mul_self,
      pauliZX_eq_neg, two_smul]
    abel
  rw [smul_mul_assoc, mul_smul_comm, smul_smul, h, smul_smul,
    invSqrt2_sq_mul_two, one_smul]

/-- **Algebra-only.** `(Z - X)/√2` is an involution. -/
theorem bObs₁_mul_self : bObs₁ * bObs₁ = 1 := by
  unfold bObs₁
  have h : (Robertson.pauliZ - Robertson.pauliX) *
      (Robertson.pauliZ - Robertson.pauliX) =
        (2 : ℂ) • (1 : Matrix (Fin 2) (Fin 2) ℂ) := by
    rw [sub_mul, mul_sub, mul_sub, pauliZ_mul_self, pauliX_mul_self,
      pauliZX_eq_neg, two_smul]
    abel
  rw [smul_mul_assoc, mul_smul_comm, smul_smul, h, smul_smul,
    invSqrt2_sq_mul_two, one_smul]

/-- **Algebra-only.** The sum of the rotated observables is `√2 · Z`. -/
theorem bObs_sum :
    bObs₀ + bObs₁ = ((Real.sqrt 2 : ℝ) : ℂ) • Robertson.pauliZ := by
  unfold bObs₀ bObs₁
  rw [← smul_add]
  have h : (Robertson.pauliZ + Robertson.pauliX) +
      (Robertson.pauliZ - Robertson.pauliX) = (2 : ℂ) • Robertson.pauliZ := by
    rw [two_smul]; abel
  rw [h, smul_smul, invSqrt2_mul_two]

/-- **Algebra-only.** The difference of the rotated observables is
`√2 · X`. -/
theorem bObs_diff :
    bObs₀ - bObs₁ = ((Real.sqrt 2 : ℝ) : ℂ) • Robertson.pauliX := by
  unfold bObs₀ bObs₁
  rw [← smul_sub]
  have h : (Robertson.pauliZ + Robertson.pauliX) -
      (Robertson.pauliZ - Robertson.pauliX) = (2 : ℂ) • Robertson.pauliX := by
    rw [two_smul]; abel
  rw [h, smul_smul, invSqrt2_mul_two]

/-- **Algebra-only.** The spectral projection `(1 - A)/2` of a self-adjoint
involution `A` is an event. -/
theorem isEvent_half_one_sub {k : ℕ} {A : Matrix (Fin k) (Fin k) ℂ}
    (hsa : A.IsHermitian) (hinv : A * A = 1) :
    IsEvent ((2⁻¹ : ℂ) • (1 - A)) := by
  constructor
  · show ((2⁻¹ : ℂ) • (1 - A))ᴴ = (2⁻¹ : ℂ) • (1 - A)
    rw [Matrix.conjTranspose_smul, Matrix.conjTranspose_sub,
      Matrix.conjTranspose_one, hsa.eq]
    congr 1
    simp
  · rw [smul_mul_assoc, mul_smul_comm, smul_smul]
    have h : ((1 : Matrix (Fin k) (Fin k) ℂ) - A) * (1 - A) =
        (2 : ℂ) • ((1 : Matrix (Fin k) (Fin k) ℂ) - A) := by
      rw [sub_mul, one_mul, mul_sub, mul_one, hinv, two_smul]
      abel
    rw [h, smul_smul]
    have h2 : (2⁻¹ * 2⁻¹ * 2 : ℂ) = 2⁻¹ := by norm_num
    rw [h2]

/-- **Algebra-only.** The dichotomic observable of the spectral projection
recovers the involution: `1 - 2 · ((1 - A)/2) = A`. -/
theorem one_sub_two_smul_half {k : ℕ} (A : Matrix (Fin k) (Fin k) ℂ) :
    (1 : Matrix (Fin k) (Fin k) ℂ) -
      2 • ((2⁻¹ : ℂ) • ((1 : Matrix (Fin k) (Fin k) ℂ) - A)) = A := by
  rw [two_smul, ← add_smul]
  have h : ((2 : ℂ)⁻¹ + (2 : ℂ)⁻¹) = 1 := by norm_num
  rw [h, one_smul, sub_sub_cancel]

/-- **Algebra-only.** The left-party event with observable `Z`. -/
def pA₀ : Matrix (Fin 2) (Fin 2) ℂ := (2⁻¹ : ℂ) • (1 - Robertson.pauliZ)

/-- **Algebra-only.** The left-party event with observable `X`. -/
def pA₁ : Matrix (Fin 2) (Fin 2) ℂ := (2⁻¹ : ℂ) • (1 - Robertson.pauliX)

/-- **Algebra-only.** The right-party event with observable `(Z + X)/√2`. -/
def pB₀ : Matrix (Fin 2) (Fin 2) ℂ := (2⁻¹ : ℂ) • (1 - bObs₀)

/-- **Algebra-only.** The right-party event with observable `(Z - X)/√2`. -/
def pB₁ : Matrix (Fin 2) (Fin 2) ℂ := (2⁻¹ : ℂ) • (1 - bObs₁)

/-- **Algebra-only.** `pA₀` is a projection event. -/
theorem pA₀_isEvent : IsEvent pA₀ :=
  isEvent_half_one_sub Robertson.pauliZ_isHermitian pauliZ_mul_self

/-- **Algebra-only.** `pA₁` is a projection event. -/
theorem pA₁_isEvent : IsEvent pA₁ :=
  isEvent_half_one_sub Robertson.pauliX_isHermitian pauliX_mul_self

/-- **Algebra-only.** `pB₀` is a projection event. -/
theorem pB₀_isEvent : IsEvent pB₀ :=
  isEvent_half_one_sub bObs₀_isHermitian bObs₀_mul_self

/-- **Algebra-only.** `pB₁` is a projection event. -/
theorem pB₁_isEvent : IsEvent pB₁ :=
  isEvent_half_one_sub bObs₁_isHermitian bObs₁_mul_self

/-- **Algebra-only.** The dichotomic observable of `pA₀` is `Z`. -/
theorem pA₀_dichotomic :
    (1 : Matrix (Fin 2) (Fin 2) ℂ) - 2 • pA₀ = Robertson.pauliZ :=
  one_sub_two_smul_half Robertson.pauliZ

/-- **Algebra-only.** The dichotomic observable of `pA₁` is `X`. -/
theorem pA₁_dichotomic :
    (1 : Matrix (Fin 2) (Fin 2) ℂ) - 2 • pA₁ = Robertson.pauliX :=
  one_sub_two_smul_half Robertson.pauliX

/-- **Algebra-only.** The dichotomic observable of `pB₀` is `(Z + X)/√2`. -/
theorem pB₀_dichotomic :
    (1 : Matrix (Fin 2) (Fin 2) ℂ) - 2 • pB₀ = bObs₀ :=
  one_sub_two_smul_half bObs₀

/-- **Algebra-only.** The dichotomic observable of `pB₁` is `(Z - X)/√2`. -/
theorem pB₁_dichotomic :
    (1 : Matrix (Fin 2) (Fin 2) ℂ) - 2 • pB₁ = bObs₁ :=
  one_sub_two_smul_half bObs₁

/-! ## The PR-44 interface witness -/

/-- **The interface witness.**  The four declared factor events inhabit the
committed PR-44 interface `EventAlgebra.SlotLocalityInterface` at
`da = db = 2`, with the singleton identity Kraus family.  The events are
declared data of this module, produced by no committed source theorem. -/
def bellSlotInterface : SlotLocalityInterface 2 2 where
  A₀ := pA₀
  A₁ := pA₁
  eventA₀ := pA₀_isEvent
  eventA₁ := pA₁_isEvent
  B₀ := pB₀
  B₁ := pB₁
  eventB₀ := pB₀_isEvent
  eventB₁ := pB₁_isEvent
  krausCard := 1
  kraus := fun _ => 1
  krausComplete := by simp

/-- **Algebra-only.** The slot-lifted CHSH operator of the interface
witness, on the joint carrier `Fin 2 × Fin 2` — verbatim the operator of
clause 1 of the committed `slot_locality_receipts`. -/
def chshProd : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  (1 - 2 • slotLeft (Fin 2) bellSlotInterface.A₀) *
      (1 - 2 • slotRight (Fin 2) bellSlotInterface.B₀) +
    (1 - 2 • slotLeft (Fin 2) bellSlotInterface.A₀) *
      (1 - 2 • slotRight (Fin 2) bellSlotInterface.B₁) +
    (1 - 2 • slotLeft (Fin 2) bellSlotInterface.A₁) *
      (1 - 2 • slotRight (Fin 2) bellSlotInterface.B₀) -
    (1 - 2 • slotLeft (Fin 2) bellSlotInterface.A₁) *
      (1 - 2 • slotRight (Fin 2) bellSlotInterface.B₁)

/-- **Algebra-only.** The committed norm-form Tsirelson bound of
`EventAlgebra.QuantumSurface.slot_locality_receipts`, cited at the
interface witness: `‖chshProd‖ ≤ 2√2` in the L2 operator norm. -/
theorem chshProd_norm_le : ‖chshProd‖ ≤ 2 * Real.sqrt 2 :=
  (QuantumSurface.slot_locality_receipts bellSlotInterface).1

/-- **Algebra-only.** Structural form of the witness CHSH operator:
`chshProd = √2 · (Z ⊗ Z + X ⊗ X)` written through the slot embeddings.
The cross-party products collapse by the committed kronecker calculus. -/
theorem chshProd_eq :
    chshProd = ((Real.sqrt 2 : ℝ) : ℂ) •
      (slotLeft (Fin 2) Robertson.pauliZ *
          slotRight (Fin 2) Robertson.pauliZ +
        slotLeft (Fin 2) Robertson.pauliX *
          slotRight (Fin 2) Robertson.pauliX) := by
  have hL : ∀ P : Matrix (Fin 2) (Fin 2) ℂ,
      slotLeft (Fin 2) ((1 : Matrix (Fin 2) (Fin 2) ℂ) - 2 • P) =
        1 - 2 • slotLeft (Fin 2) P := fun P => by
    rw [map_sub, map_one, map_nsmul]
  have hR : ∀ P : Matrix (Fin 2) (Fin 2) ℂ,
      slotRight (Fin 2) ((1 : Matrix (Fin 2) (Fin 2) ℂ) - 2 • P) =
        1 - 2 • slotRight (Fin 2) P := fun P => by
    rw [map_sub, map_one, map_nsmul]
  have hA₀ : (1 : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) -
      2 • slotLeft (Fin 2) bellSlotInterface.A₀ =
        slotLeft (Fin 2) Robertson.pauliZ := by
    rw [show bellSlotInterface.A₀ = pA₀ from rfl, ← hL, pA₀_dichotomic]
  have hA₁ : (1 : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) -
      2 • slotLeft (Fin 2) bellSlotInterface.A₁ =
        slotLeft (Fin 2) Robertson.pauliX := by
    rw [show bellSlotInterface.A₁ = pA₁ from rfl, ← hL, pA₁_dichotomic]
  have hB₀ : (1 : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) -
      2 • slotRight (Fin 2) bellSlotInterface.B₀ =
        slotRight (Fin 2) bObs₀ := by
    rw [show bellSlotInterface.B₀ = pB₀ from rfl, ← hR, pB₀_dichotomic]
  have hB₁ : (1 : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) -
      2 • slotRight (Fin 2) bellSlotInterface.B₁ =
        slotRight (Fin 2) bObs₁ := by
    rw [show bellSlotInterface.B₁ = pB₁ from rfl, ← hR, pB₁_dichotomic]
  unfold chshProd
  rw [hA₀, hA₁, hB₀, hB₁,
    ← mul_add (slotLeft (Fin 2) Robertson.pauliZ), add_sub_assoc,
    ← mul_sub (slotLeft (Fin 2) Robertson.pauliX),
    ← map_add (slotRight (Fin 2)), ← map_sub (slotRight (Fin 2)),
    bObs_sum, bObs_diff, map_smul, map_smul,
    mul_smul_comm, mul_smul_comm, ← smul_add]

/-! ## The declared state on the joint carrier -/

/-- **The declared state**, product-carrier form: the normalized
`OPH.Dynamics.maxEntangled 2`.  DECLARED, not source-produced: no committed
theorem prepares it (register row PR-52; OL-C3 realization clause). -/
def bellStateProd : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ :=
  (2⁻¹ : ℂ) • OPH.Dynamics.maxEntangled 2

/-- **Algebra-only.** The declared state is positive semidefinite. -/
theorem bellStateProd_posSemidef : bellStateProd.PosSemidef := by
  have h : (0 : ℂ) ≤ (2⁻¹ : ℂ) := by
    rw [show (2⁻¹ : ℂ) = ((2⁻¹ : ℝ) : ℂ) by norm_num]
    exact_mod_cast (by norm_num : (0 : ℝ) ≤ 2⁻¹)
  exact (OPH.Dynamics.maxEntangled_posSemidef 2).smul h

private theorem maxEntangled_two_trace :
    (OPH.Dynamics.maxEntangled 2).trace = 2 := by
  norm_num [OPH.Dynamics.maxEntangled, Matrix.trace, Matrix.diag,
    Matrix.vecMulVec_apply, Fintype.sum_prod_type, Fin.sum_univ_two,
    Pi.star_apply]

/-- **Trace-dependent.** The declared state has trace one. -/
theorem bellStateProd_trace : bellStateProd.trace = 1 := by
  unfold bellStateProd
  rw [Matrix.trace_smul, maxEntangled_two_trace, smul_eq_mul]
  norm_num

private theorem slot_pair_kronecker (A B : Matrix (Fin 2) (Fin 2) ℂ) :
    slotLeft (Fin 2) A * slotRight (Fin 2) B = A ⊗ₖ B := by
  show (A ⊗ₖ (1 : Matrix (Fin 2) (Fin 2) ℂ)) *
      ((1 : Matrix (Fin 2) (Fin 2) ℂ) ⊗ₖ B) = A ⊗ₖ B
  rw [← Matrix.mul_kronecker_mul, mul_one, one_mul]

private theorem pairing_ZZ :
    (OPH.Dynamics.maxEntangled 2 *
      (Robertson.pauliZ ⊗ₖ Robertson.pauliZ)).trace = 2 := by
  norm_num [OPH.Dynamics.maxEntangled, Matrix.trace, Matrix.diag,
    Matrix.mul_apply, Matrix.vecMulVec_apply, Matrix.kroneckerMap_apply,
    Fintype.sum_prod_type, Fin.sum_univ_two, Pi.star_apply,
    Robertson.pauliZ]

private theorem pairing_XX :
    (OPH.Dynamics.maxEntangled 2 *
      (Robertson.pauliX ⊗ₖ Robertson.pauliX)).trace = 2 := by
  norm_num [OPH.Dynamics.maxEntangled, Matrix.trace, Matrix.diag,
    Matrix.mul_apply, Matrix.vecMulVec_apply, Matrix.kroneckerMap_apply,
    Fintype.sum_prod_type, Fin.sum_univ_two, Pi.star_apply,
    Robertson.pauliX]

/-- **Trace-dependent.** The exact product-carrier pairing: the declared
state pairs with the witness CHSH operator to `2√2` exactly. -/
theorem bellStateProd_chshProd_trace :
    (bellStateProd * chshProd).trace = ((2 * Real.sqrt 2 : ℝ) : ℂ) := by
  rw [chshProd_eq, slot_pair_kronecker, slot_pair_kronecker]
  unfold bellStateProd
  rw [smul_mul_assoc, mul_smul_comm, smul_smul, Matrix.trace_smul,
    mul_add, Matrix.trace_add, pairing_ZZ, pairing_XX, smul_eq_mul]
  push_cast
  ring

/-! ## Transport to `Fin 4` along the fixed index equivalence

The committed state-expectation vocabulary (`EventAlgebra.expectation`,
`IsState`, `IsEvent`, `matrix_state_tsirelson_bound_of_events`) is typed
over `Fin n`; the slot machinery is typed over `Fin 2 × Fin 2`.  The fixed
equivalence `finProdFinEquiv` transports between the two carriers; every
transfer lemma below is elementary reindexing. -/

/-- The fixed index equivalence between the joint slot carrier and
`Fin 4`. -/
def carrierEquiv : Fin 2 × Fin 2 ≃ Fin 4 := finProdFinEquiv

/-- Reindexing along the fixed carrier equivalence. -/
def toFour (M : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) :
    Matrix (Fin 4) (Fin 4) ℂ :=
  Matrix.reindex carrierEquiv carrierEquiv M

private theorem toFour_map_mul
    (M N : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) :
    toFour M * toFour N = toFour (M * N) := by
  simp only [toFour, Matrix.reindex_apply, Matrix.submatrix_mul_equiv]

private theorem toFour_one : toFour 1 = 1 := by
  simp only [toFour, Matrix.reindex_apply]
  exact Matrix.submatrix_one_equiv carrierEquiv.symm

private theorem toFour_conjTranspose
    (M : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) :
    (toFour M)ᴴ = toFour Mᴴ := by
  simp only [toFour, Matrix.reindex_apply, Matrix.conjTranspose_submatrix]

private theorem toFour_add
    (M N : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) :
    toFour (M + N) = toFour M + toFour N := by
  ext i j
  simp [toFour, Matrix.reindex_apply, Matrix.submatrix_apply]

private theorem toFour_sub
    (M N : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) :
    toFour (M - N) = toFour M - toFour N := by
  ext i j
  simp [toFour, Matrix.reindex_apply, Matrix.submatrix_apply]

private theorem toFour_nsmul (n : ℕ)
    (M : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) :
    toFour (n • M) = n • toFour M := by
  ext i j
  simp only [toFour, Matrix.reindex_apply, Matrix.submatrix_apply,
    Matrix.smul_apply]

private theorem toFour_trace
    (M : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) :
    (toFour M).trace = M.trace := by
  unfold toFour Matrix.trace Matrix.diag
  simp only [Matrix.reindex_apply, Matrix.submatrix_apply]
  exact Equiv.sum_comp carrierEquiv.symm (fun p => M p p)

private theorem toFour_posSemidef
    {M : Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ} (h : M.PosSemidef) :
    (toFour M).PosSemidef := by
  unfold toFour
  rw [Matrix.reindex_apply]
  exact (Matrix.posSemidef_submatrix_equiv carrierEquiv.symm).mpr h

/-- **The declared state on `Fin 4`**: the reindexed normalized
`maxEntangled 2`.  DECLARED, not source-produced (PR-52; OL-C3 realization
clause). -/
def bellState : Matrix (Fin 4) (Fin 4) ℂ := toFour bellStateProd

/-- **Trace-dependent.** The declared state is a state: positive
semidefinite with trace one. -/
theorem bellState_isState : IsState bellState := by
  refine ⟨toFour_posSemidef bellStateProd_posSemidef, ?_⟩
  show (toFour bellStateProd).trace = 1
  rw [toFour_trace, bellStateProd_trace]

/-- **Algebra-only.** The reindexed slot-lifted left event `pA₀ ⊗ 1`. -/
def qA₀ : Matrix (Fin 4) (Fin 4) ℂ := toFour (slotLeft (Fin 2) pA₀)

/-- **Algebra-only.** The reindexed slot-lifted left event `pA₁ ⊗ 1`. -/
def qA₁ : Matrix (Fin 4) (Fin 4) ℂ := toFour (slotLeft (Fin 2) pA₁)

/-- **Algebra-only.** The reindexed slot-lifted right event `1 ⊗ pB₀`. -/
def qB₀ : Matrix (Fin 4) (Fin 4) ℂ := toFour (slotRight (Fin 2) pB₀)

/-- **Algebra-only.** The reindexed slot-lifted right event `1 ⊗ pB₁`. -/
def qB₁ : Matrix (Fin 4) (Fin 4) ℂ := toFour (slotRight (Fin 2) pB₁)

private theorem isEvent_toFour_lift {k : ℕ}
    (φ : Matrix (Fin k) (Fin k) ℂ →⋆ₐ[ℂ]
      Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ)
    {P : Matrix (Fin k) (Fin k) ℂ} (hP : IsEvent P) :
    IsEvent (toFour (φ P)) := by
  constructor
  · show (toFour (φ P))ᴴ = toFour (φ P)
    rw [toFour_conjTranspose]
    congr 1
    rw [← Matrix.star_eq_conjTranspose, ← map_star]
    congr 1
    rw [Matrix.star_eq_conjTranspose, hP.1.eq]
  · rw [toFour_map_mul, ← map_mul, hP.2]

/-- **Algebra-only.** `qA₀` is a projection event of `Fin 4`. -/
theorem qA₀_isEvent : IsEvent qA₀ := isEvent_toFour_lift _ pA₀_isEvent

/-- **Algebra-only.** `qA₁` is a projection event of `Fin 4`. -/
theorem qA₁_isEvent : IsEvent qA₁ := isEvent_toFour_lift _ pA₁_isEvent

/-- **Algebra-only.** `qB₀` is a projection event of `Fin 4`. -/
theorem qB₀_isEvent : IsEvent qB₀ := isEvent_toFour_lift _ pB₀_isEvent

/-- **Algebra-only.** `qB₁` is a projection event of `Fin 4`. -/
theorem qB₁_isEvent : IsEvent qB₁ := isEvent_toFour_lift _ pB₁_isEvent

private theorem q_commute_left_right (P Q : Matrix (Fin 2) (Fin 2) ℂ) :
    toFour (slotLeft (Fin 2) P) * toFour (slotRight (Fin 2) Q) =
      toFour (slotRight (Fin 2) Q) * toFour (slotLeft (Fin 2) P) := by
  rw [toFour_map_mul, toFour_map_mul]
  congr 1
  exact (slot_commute ⟨P, rfl⟩ ⟨Q, rfl⟩).eq

/-- **Algebra-only.** Cross-party commutation on `Fin 4`, from the
committed slot-disjointness theorem `OPH.QFT.slot_commute` — not from
entry computation. -/
theorem qA₀_qB₀_comm : qA₀ * qB₀ = qB₀ * qA₀ := q_commute_left_right pA₀ pB₀

/-- **Algebra-only.** Cross-party commutation on `Fin 4` (see
`qA₀_qB₀_comm`). -/
theorem qA₀_qB₁_comm : qA₀ * qB₁ = qB₁ * qA₀ := q_commute_left_right pA₀ pB₁

/-- **Algebra-only.** Cross-party commutation on `Fin 4` (see
`qA₀_qB₀_comm`). -/
theorem qA₁_qB₀_comm : qA₁ * qB₀ = qB₀ * qA₁ := q_commute_left_right pA₁ pB₀

/-- **Algebra-only.** Cross-party commutation on `Fin 4` (see
`qA₀_qB₀_comm`). -/
theorem qA₁_qB₁_comm : qA₁ * qB₁ = qB₁ * qA₁ := q_commute_left_right pA₁ pB₁

/-- **Algebra-only.** The CHSH observable on `Fin 4`, in the committed
event convention `1 - 2P` — the exact operator shape consumed by
`matrix_state_tsirelson_bound_of_events`. -/
def chshObs : Matrix (Fin 4) (Fin 4) ℂ :=
  (1 - 2 • qA₀) * (1 - 2 • qB₀) + (1 - 2 • qA₀) * (1 - 2 • qB₁) +
    (1 - 2 • qA₁) * (1 - 2 • qB₀) - (1 - 2 • qA₁) * (1 - 2 • qB₁)

/-- **Algebra-only.** The `Fin 4` CHSH observable is the reindexed witness
CHSH operator of the PR-44 interface: the attaining observable and the
interface observable are the same operator up to the fixed carrier
equivalence. -/
theorem toFour_chshProd_eq : toFour chshProd = chshObs := by
  show toFour
      ((1 - 2 • slotLeft (Fin 2) pA₀) * (1 - 2 • slotRight (Fin 2) pB₀) +
        (1 - 2 • slotLeft (Fin 2) pA₀) * (1 - 2 • slotRight (Fin 2) pB₁) +
        (1 - 2 • slotLeft (Fin 2) pA₁) * (1 - 2 • slotRight (Fin 2) pB₀) -
        (1 - 2 • slotLeft (Fin 2) pA₁) * (1 - 2 • slotRight (Fin 2) pB₁)) =
    chshObs
  simp only [toFour_sub, toFour_add, ← toFour_map_mul, toFour_one,
    toFour_nsmul]
  rfl

/-! ## The saturation theorems -/

/-- **Trace-dependent.** **Exact Tsirelson saturation.**  The declared
state gives CHSH expectation exactly `2√2` on the committed expectation
functional.  The state and events are declared data (PR-52; OL-C3
realization clause); the equality itself is exact finite matrix algebra. -/
theorem chsh_saturation :
    expectation bellState chshObs = ((2 * Real.sqrt 2 : ℝ) : ℂ) := by
  have h1 : expectation bellState chshObs =
      (bellState * chshObs).trace := rfl
  rw [h1, ← toFour_chshProd_eq]
  show (toFour bellStateProd * toFour chshProd).trace = _
  rw [toFour_map_mul, toFour_trace]
  exact bellStateProd_chshProd_trace

/-- **Trace-dependent.** Norm form of the saturation: the modulus of the
CHSH expectation of the declared state equals the committed upper bound
`2√2`. -/
theorem chsh_saturation_norm :
    ‖expectation bellState chshObs‖ = 2 * Real.sqrt 2 := by
  rw [chsh_saturation, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg (by positivity)]

/-- **Trace-dependent.** The committed state-level upper bound
`matrix_state_tsirelson_bound_of_events`, cited at the four reindexed
witness events: every state on `Fin 4` obeys `‖⟨CHSH⟩‖ ≤ 2√2`. -/
theorem chsh_expectation_bound {ρ : Matrix (Fin 4) (Fin 4) ℂ}
    (hρ : IsState ρ) :
    ‖expectation ρ chshObs‖ ≤ 2 * Real.sqrt 2 := by
  unfold chshObs
  exact matrix_state_tsirelson_bound_of_events hρ qA₀ qA₁ qB₀ qB₁
    qA₀_isEvent qA₁_isEvent qB₀_isEvent qB₁_isEvent
    qA₀_qB₀_comm qA₀_qB₁_comm qA₁_qB₀_comm qA₁_qB₁_comm

/-- **The composed OL-C3 saturation receipt.**  One conjunction whose
clauses share the same interface witness, the same events, and the same
declared state:

1. **Upper bound at the witness** — the committed norm-form Tsirelson
   bound of `slot_locality_receipts`, cited at `bellSlotInterface`;
2. **Upper bound for every state** — the committed
   `matrix_state_tsirelson_bound_of_events` at the reindexed witness
   events;
3. **Exact attainment** — the declared state is a state and its CHSH
   expectation equals `2√2` exactly, in value and in modulus;
4. **Interface identification** — the attaining observable is the
   reindexed slot-lifted CHSH operator of the PR-44 interface witness.

Boundary, stated exactly: the state and events are DECLARED data of this
module; no committed theorem prepares the state from source runs.  Source
entangled preparation, physical regions, spacelike separation, and
instruments stay open (register row PR-52; OL-C3 realization clause).
The PR-44 split is a supplied datum. -/
theorem tsirelson_saturation_receipt :
    (‖chshProd‖ ≤ 2 * Real.sqrt 2) ∧
    (∀ ρ : Matrix (Fin 4) (Fin 4) ℂ, IsState ρ →
      ‖expectation ρ chshObs‖ ≤ 2 * Real.sqrt 2) ∧
    (IsState bellState ∧
      expectation bellState chshObs = ((2 * Real.sqrt 2 : ℝ) : ℂ) ∧
      ‖expectation bellState chshObs‖ = 2 * Real.sqrt 2) ∧
    chshObs = toFour chshProd :=
  ⟨chshProd_norm_le, fun _ hρ => chsh_expectation_bound hρ,
    ⟨bellState_isState, chsh_saturation, chsh_saturation_norm⟩,
    toFour_chshProd_eq.symm⟩

/-! ## The classical diagonal delimitation

A state diagonal in the record basis, read against any four diagonal
`[-1,1]`-valued observables, cannot exceed the classical CHSH value `2`.
This is a common-static-record statement: the same diagonal state and four
simultaneously defined diagonal functions appear in one CHSH expression.  It
does not classify arbitrary context-dependent operational models.  In
particular, it makes no necessity or sufficiency claim about PR-04 or complex
phase structure; the saturating Bell construction above is entirely real. -/

/-- **Premise-free real arithmetic.** The pointwise classical CHSH bound:
for reals of modulus at most one,
`|a₀b₀ + a₀b₁ + a₁b₀ - a₁b₁| ≤ 2`. -/
theorem chsh_real_le_two {a₀ a₁ b₀ b₁ : ℝ}
    (ha₀ : |a₀| ≤ 1) (ha₁ : |a₁| ≤ 1) (hb₀ : |b₀| ≤ 1) (hb₁ : |b₁| ≤ 1) :
    |a₀ * b₀ + a₀ * b₁ + a₁ * b₀ - a₁ * b₁| ≤ 2 := by
  have key : a₀ * b₀ + a₀ * b₁ + a₁ * b₀ - a₁ * b₁ =
      a₀ * (b₀ + b₁) + a₁ * (b₀ - b₁) := by ring
  rw [key]
  have h1 : |a₀ * (b₀ + b₁) + a₁ * (b₀ - b₁)| ≤
      |b₀ + b₁| + |b₀ - b₁| := by
    calc |a₀ * (b₀ + b₁) + a₁ * (b₀ - b₁)|
        ≤ |a₀ * (b₀ + b₁)| + |a₁ * (b₀ - b₁)| := abs_add_le _ _
      _ = |a₀| * |b₀ + b₁| + |a₁| * |b₀ - b₁| := by rw [abs_mul, abs_mul]
      _ ≤ 1 * |b₀ + b₁| + 1 * |b₀ - b₁| := by gcongr
      _ = |b₀ + b₁| + |b₀ - b₁| := by ring
  have hb₀' := abs_le.mp hb₀
  have hb₁' := abs_le.mp hb₁
  have h2 : |b₀ + b₁| + |b₀ - b₁| ≤ 2 := by
    rcases abs_cases (b₀ + b₁) with ⟨e1, _⟩ | ⟨e1, _⟩ <;>
      rcases abs_cases (b₀ - b₁) with ⟨e2, _⟩ | ⟨e2, _⟩ <;>
      rw [e1, e2] <;> linarith [hb₀'.1, hb₀'.2, hb₁'.1, hb₁'.2]
  linarith

/-- **Premise-free real arithmetic.** The weighted classical CHSH bound:
convex weights against pointwise-bounded readouts stay within `2`. -/
theorem diagonal_chsh_le_two {ι : Type*} [Fintype ι]
    (w a₀ a₁ b₀ b₁ : ι → ℝ)
    (hw : ∀ i, 0 ≤ w i) (hsum : ∑ i, w i = 1)
    (ha₀ : ∀ i, |a₀ i| ≤ 1) (ha₁ : ∀ i, |a₁ i| ≤ 1)
    (hb₀ : ∀ i, |b₀ i| ≤ 1) (hb₁ : ∀ i, |b₁ i| ≤ 1) :
    |∑ i, w i *
        (a₀ i * b₀ i + a₀ i * b₁ i + a₁ i * b₀ i - a₁ i * b₁ i)| ≤ 2 := by
  calc |∑ i, w i *
        (a₀ i * b₀ i + a₀ i * b₁ i + a₁ i * b₀ i - a₁ i * b₁ i)|
      ≤ ∑ i, |w i *
        (a₀ i * b₀ i + a₀ i * b₁ i + a₁ i * b₀ i - a₁ i * b₁ i)| :=
        Finset.abs_sum_le_sum_abs _ _
    _ = ∑ i, w i *
        |a₀ i * b₀ i + a₀ i * b₁ i + a₁ i * b₀ i - a₁ i * b₁ i| := by
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [abs_mul, abs_of_nonneg (hw i)]
    _ ≤ ∑ i, w i * 2 :=
        Finset.sum_le_sum fun i _ =>
          mul_le_mul_of_nonneg_left
            (chsh_real_le_two (ha₀ i) (ha₁ i) (hb₀ i) (hb₁ i)) (hw i)
    _ = 2 := by rw [← Finset.sum_mul, hsum, one_mul]

/-- **Trace-dependent.** The matrix form of the diagonal delimitation: a
diagonal state with convex weights, read against four diagonal
`[-1,1]`-valued observables in the CHSH combination, has expectation of
modulus at most `2`. -/
theorem diagonal_state_chsh_le_two {ι : Type*} [Fintype ι] [DecidableEq ι]
    (w a₀ a₁ b₀ b₁ : ι → ℝ)
    (hw : ∀ i, 0 ≤ w i) (hsum : ∑ i, w i = 1)
    (ha₀ : ∀ i, |a₀ i| ≤ 1) (ha₁ : ∀ i, |a₁ i| ≤ 1)
    (hb₀ : ∀ i, |b₀ i| ≤ 1) (hb₁ : ∀ i, |b₁ i| ≤ 1) :
    ‖(Matrix.diagonal (fun i => (w i : ℂ)) *
        (Matrix.diagonal (fun i => (a₀ i : ℂ)) *
            Matrix.diagonal (fun i => (b₀ i : ℂ)) +
          Matrix.diagonal (fun i => (a₀ i : ℂ)) *
            Matrix.diagonal (fun i => (b₁ i : ℂ)) +
          Matrix.diagonal (fun i => (a₁ i : ℂ)) *
            Matrix.diagonal (fun i => (b₀ i : ℂ)) -
          Matrix.diagonal (fun i => (a₁ i : ℂ)) *
            Matrix.diagonal (fun i => (b₁ i : ℂ)))).trace‖ ≤ 2 := by
  have hred : (Matrix.diagonal (fun i => (w i : ℂ)) *
      (Matrix.diagonal (fun i => (a₀ i : ℂ)) *
          Matrix.diagonal (fun i => (b₀ i : ℂ)) +
        Matrix.diagonal (fun i => (a₀ i : ℂ)) *
          Matrix.diagonal (fun i => (b₁ i : ℂ)) +
        Matrix.diagonal (fun i => (a₁ i : ℂ)) *
          Matrix.diagonal (fun i => (b₀ i : ℂ)) -
        Matrix.diagonal (fun i => (a₁ i : ℂ)) *
          Matrix.diagonal (fun i => (b₁ i : ℂ)))).trace =
      ((∑ i, w i *
        (a₀ i * b₀ i + a₀ i * b₁ i + a₁ i * b₀ i - a₁ i * b₁ i) : ℝ) : ℂ) := by
    simp only [Matrix.diagonal_mul_diagonal, Matrix.diagonal_add,
      Matrix.diagonal_sub, Matrix.trace_diagonal]
    push_cast
    ring_nf
  rw [hred, Complex.norm_real, Real.norm_eq_abs]
  exact diagonal_chsh_le_two w a₀ a₁ b₀ b₁ hw hsum ha₀ ha₁ hb₀ hb₁

/-- **Trace-dependent.** **The counted-state delimitation.**  The committed
counted state `OPH.QFT.correlationState` — diagonal in the record basis by
construction — obeys the classical CHSH bound `2` against every four
diagonal `[-1,1]`-valued readouts.  The committed diagonal record surface
therefore does not attain the Tsirelson value `2√2` through this fixed-state,
jointly diagonal readout model.  Only that scoped diagonality is insufficient;
no implication involving PR-04 or genuinely complex phase structure is
asserted. -/
theorem counted_state_diagonal_chsh_le_two
    (a₀ a₁ b₀ b₁ : Fin 13 × Fin 13 → ℝ)
    (ha₀ : ∀ p, |a₀ p| ≤ 1) (ha₁ : ∀ p, |a₁ p| ≤ 1)
    (hb₀ : ∀ p, |b₀ p| ≤ 1) (hb₁ : ∀ p, |b₁ p| ≤ 1) :
    ‖(OPH.QFT.correlationState *
        (Matrix.diagonal (fun p => (a₀ p : ℂ)) *
            Matrix.diagonal (fun p => (b₀ p : ℂ)) +
          Matrix.diagonal (fun p => (a₀ p : ℂ)) *
            Matrix.diagonal (fun p => (b₁ p : ℂ)) +
          Matrix.diagonal (fun p => (a₁ p : ℂ)) *
            Matrix.diagonal (fun p => (b₀ p : ℂ)) -
          Matrix.diagonal (fun p => (a₁ p : ℂ)) *
            Matrix.diagonal (fun p => (b₁ p : ℂ)))).trace‖ ≤ 2 := by
  have hcs : OPH.QFT.correlationState =
      Matrix.diagonal fun p =>
        (((OPH.QFT.jointCount p : ℝ) / 32 : ℝ) : ℂ) := by
    unfold OPH.QFT.correlationState
    congr 1
    funext p
    push_cast
    ring
  rw [hcs]
  refine diagonal_state_chsh_le_two _ a₀ a₁ b₀ b₁
    (fun p => by positivity) ?_ ha₀ ha₁ hb₀ hb₁
  rw [← Finset.sum_div,
    show (∑ p : Fin 13 × Fin 13, (OPH.QFT.jointCount p : ℝ)) =
      ((∑ p : Fin 13 × Fin 13, OPH.QFT.jointCount p : ℕ) : ℝ) by
        push_cast; rfl,
    OPH.QFT.jointCount_total]
  norm_num

end

end EventAlgebra.TsirelsonSaturation

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
-- No native_decide is used anywhere in this module.
#print axioms EventAlgebra.TsirelsonSaturation.bObs₀_isHermitian
#print axioms EventAlgebra.TsirelsonSaturation.bObs₁_isHermitian
#print axioms EventAlgebra.TsirelsonSaturation.bObs₀_mul_self
#print axioms EventAlgebra.TsirelsonSaturation.bObs₁_mul_self
#print axioms EventAlgebra.TsirelsonSaturation.bObs_sum
#print axioms EventAlgebra.TsirelsonSaturation.bObs_diff
#print axioms EventAlgebra.TsirelsonSaturation.isEvent_half_one_sub
#print axioms EventAlgebra.TsirelsonSaturation.one_sub_two_smul_half
#print axioms EventAlgebra.TsirelsonSaturation.pA₀_isEvent
#print axioms EventAlgebra.TsirelsonSaturation.pA₁_isEvent
#print axioms EventAlgebra.TsirelsonSaturation.pB₀_isEvent
#print axioms EventAlgebra.TsirelsonSaturation.pB₁_isEvent
#print axioms EventAlgebra.TsirelsonSaturation.pA₀_dichotomic
#print axioms EventAlgebra.TsirelsonSaturation.pA₁_dichotomic
#print axioms EventAlgebra.TsirelsonSaturation.pB₀_dichotomic
#print axioms EventAlgebra.TsirelsonSaturation.pB₁_dichotomic
#print axioms EventAlgebra.TsirelsonSaturation.chshProd_norm_le
#print axioms EventAlgebra.TsirelsonSaturation.chshProd_eq
#print axioms EventAlgebra.TsirelsonSaturation.bellStateProd_posSemidef
#print axioms EventAlgebra.TsirelsonSaturation.bellStateProd_trace
#print axioms EventAlgebra.TsirelsonSaturation.bellStateProd_chshProd_trace
#print axioms EventAlgebra.TsirelsonSaturation.bellState_isState
#print axioms EventAlgebra.TsirelsonSaturation.qA₀_isEvent
#print axioms EventAlgebra.TsirelsonSaturation.qA₁_isEvent
#print axioms EventAlgebra.TsirelsonSaturation.qB₀_isEvent
#print axioms EventAlgebra.TsirelsonSaturation.qB₁_isEvent
#print axioms EventAlgebra.TsirelsonSaturation.qA₀_qB₀_comm
#print axioms EventAlgebra.TsirelsonSaturation.qA₀_qB₁_comm
#print axioms EventAlgebra.TsirelsonSaturation.qA₁_qB₀_comm
#print axioms EventAlgebra.TsirelsonSaturation.qA₁_qB₁_comm
#print axioms EventAlgebra.TsirelsonSaturation.toFour_chshProd_eq
#print axioms EventAlgebra.TsirelsonSaturation.chsh_saturation
#print axioms EventAlgebra.TsirelsonSaturation.chsh_saturation_norm
#print axioms EventAlgebra.TsirelsonSaturation.chsh_expectation_bound
#print axioms EventAlgebra.TsirelsonSaturation.tsirelson_saturation_receipt
#print axioms EventAlgebra.TsirelsonSaturation.chsh_real_le_two
#print axioms EventAlgebra.TsirelsonSaturation.diagonal_chsh_le_two
#print axioms EventAlgebra.TsirelsonSaturation.diagonal_state_chsh_le_two
#print axioms EventAlgebra.TsirelsonSaturation.counted_state_diagonal_chsh_le_two
