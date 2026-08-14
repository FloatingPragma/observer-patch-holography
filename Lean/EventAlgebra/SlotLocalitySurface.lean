import EventAlgebra.Tsirelson
import QFT.JointSlotFactorisation

set_option autoImplicit false

namespace EventAlgebra

/-!
# The slot-locality surface: Tsirelson and no-signalling on one split

One typed interface for the observation-ledger rows OL-C3 (entanglement
to the Tsirelson bound) and OL-C4 (no-signalling and locality receipts),
V3 issue #730.  `SlotLocalityInterface` bundles the register row PR-44
(the supplied bipartite slot split with four factor events and one
Kraus-complete left-local family); the joint carrier, the trace pairing
of the conclusions, and the effect reading are the represented surface
of register row PR-02.  The one conclusion theorem
`QuantumSurface.slot_locality_receipts` states, on the SAME joint
carrier through the SAME slot embeddings `slotLeft`/`slotRight` and the
SAME Kraus family:

1. the Tsirelson bound for the four slot-lifted events, with the
   cross-party commutation DERIVED from the slot split
   (`OPH.QFT.slot_commute`) rather than supplied as a hypothesis;
2. exact no-signalling: the lifted Kraus channel preserves the remote
   partial trace for every joint matrix;
3. Born-pairing invariance: the trace pairing of every right-slot
   observable against any joint matrix is invariant under the left
   channel, connecting the locality receipt to the PR-02 trace-pairing
   vocabulary.

Composition boundaries, stated exactly: the Tsirelson clause is an
upper bound with no attainment claim and no entangled-state
realization; the no-signalling clauses are finite algebraic statements
about one declared tensor split, with no physical spacelike separation,
laboratory channel, or instrument attached (register row PR-52 records
that open premise); the bipartite split itself is the declared PR-44
factorization, produced by no theorem of this module.  Everything is
finite-dimensional matrix algebra over the L2 operator norm.
-/

open Matrix
open OPH.QFT
open Kronecker

/-- **Register row PR-44: the bipartite slot split with a local
operation.**  The supplied data of one bipartite factorization: two
factor events per party on the two factor algebras, and one
Kraus-complete family on the left factor.  The split itself enters the
conclusions as the fixed slot embeddings `slotLeft`/`slotRight` on the
joint carrier `Fin da × Fin db`, the declared PR-44 factorization.
Events and trace pairings are read on the represented matrix algebra of
register row PR-02.  No theorem of this module produces the split, the
events, or the Kraus family from source data. -/
structure SlotLocalityInterface (da db : ℕ) where
  /-- First left-party factor event (PR-44 supplied datum). -/
  A₀ : Matrix (Fin da) (Fin da) ℂ
  /-- Second left-party factor event (PR-44 supplied datum). -/
  A₁ : Matrix (Fin da) (Fin da) ℂ
  /-- `A₀` is a projection event of the left factor algebra. -/
  eventA₀ : IsEvent A₀
  /-- `A₁` is a projection event of the left factor algebra. -/
  eventA₁ : IsEvent A₁
  /-- First right-party factor event (PR-44 supplied datum). -/
  B₀ : Matrix (Fin db) (Fin db) ℂ
  /-- Second right-party factor event (PR-44 supplied datum). -/
  B₁ : Matrix (Fin db) (Fin db) ℂ
  /-- `B₀` is a projection event of the right factor algebra. -/
  eventB₀ : IsEvent B₀
  /-- `B₁` is a projection event of the right factor algebra. -/
  eventB₁ : IsEvent B₁
  /-- Size of the left-local Kraus family (PR-44 supplied datum). -/
  krausCard : ℕ
  /-- The left-local Kraus family (PR-44 supplied datum). -/
  kraus : Fin krausCard → Matrix (Fin da) (Fin da) ℂ
  /-- Kraus completeness: the family is trace preserving,
  `∑ Kᵢᴴ Kᵢ = 1`. -/
  krausComplete : ∑ i, (kraus i)ᴴ * kraus i = 1

namespace QuantumSurface

section SlotLocality

open scoped Matrix.Norms.L2Operator

variable {α β : Type*} [Fintype α] [Fintype β] [DecidableEq α] [DecidableEq β]

/-- **Trace-dependent.**  The trace pairing of a joint matrix against a
right-slot observable factors through the partial trace over the left
slot: `Tr(M · (1 ⊗ B)) = Tr(ptraceFst M · B)`.  This is the bridge from
the locality receipts to the PR-02 trace-pairing vocabulary. -/
theorem trace_mul_slotRight (M : Matrix (α × β) (α × β) ℂ)
    (B : Matrix β β ℂ) :
    (M * slotRight α B).trace = (OPH.Locality.ptraceFst M * B).trace := by
  have hentry : ∀ z : α × β,
      (M * slotRight α B) z z = ∑ b' : β, M z (z.1, b') * B b' z.2 := by
    intro z
    show (M * ((1 : Matrix α α ℂ) ⊗ₖ B)) z z = _
    rw [Matrix.mul_apply, Fintype.sum_prod_type]
    rw [Finset.sum_eq_single z.1
      (fun a' _ ha' => Finset.sum_eq_zero fun b' _ => by
        simp [Matrix.kroneckerMap_apply, ha'])
      (fun h => absurd (Finset.mem_univ _) h)]
    refine Finset.sum_congr rfl fun b' _ => ?_
    simp [Matrix.kroneckerMap_apply]
  unfold Matrix.trace Matrix.diag
  calc ∑ z : α × β, (M * slotRight α B) z z
      = ∑ z : α × β, ∑ b' : β, M z (z.1, b') * B b' z.2 :=
        Finset.sum_congr rfl fun z _ => hentry z
    _ = ∑ b : β, ∑ b' : β, (∑ a : α, M (a, b) (a, b')) * B b' b := by
        rw [Fintype.sum_prod_type, Finset.sum_comm]
        refine Finset.sum_congr rfl fun b _ => ?_
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun b' _ => ?_
        rw [Finset.sum_mul]
    _ = ∑ b : β, (OPH.Locality.ptraceFst M * B) b b := by
        refine Finset.sum_congr rfl fun b _ => ?_
        rw [Matrix.mul_apply]
        rfl

/-- **The slot-locality receipts, composed.**  For every PR-44
interface on the declared slot split of the joint carrier
`Fin da × Fin db`, one conjunction whose clauses share the same
embeddings and the same Kraus family:

1. **Tsirelson bound.**  The four slot-lifted events give dichotomic
   observables whose CHSH combination has L2 operator norm at most
   `2√2`; the cross-party commutation is derived from the slot split
   (`slot_commute`), consuming no supplied commutation hypothesis.  The
   bound is an inequality with no attainment claim and no
   entangled-state realization.
2. **No-signalling.**  Conjugating any joint matrix by the lifted
   left-local Kraus family leaves the partial trace over the left slot
   exactly invariant.
3. **Born-pairing invariance.**  The trace pairing of every right-slot
   observable against any joint matrix is invariant under the left
   channel: the PR-02 pairing read by the remote party sees no change.

Finite/modal boundary: all clauses are exact finite matrix algebra on
one declared tensor split; no physical spacelike separation, laboratory
channel, or instrument enters (register row PR-52 stays open), and the
split is the declared PR-44 datum. -/
theorem slot_locality_receipts {da db : ℕ} [NeZero da] [NeZero db]
    (S : SlotLocalityInterface da db) :
    (‖(1 - 2 • slotLeft (Fin db) S.A₀) * (1 - 2 • slotRight (Fin da) S.B₀) +
        (1 - 2 • slotLeft (Fin db) S.A₀) * (1 - 2 • slotRight (Fin da) S.B₁) +
        (1 - 2 • slotLeft (Fin db) S.A₁) * (1 - 2 • slotRight (Fin da) S.B₀) -
        (1 - 2 • slotLeft (Fin db) S.A₁) * (1 - 2 • slotRight (Fin da) S.B₁)‖
      ≤ 2 * Real.sqrt 2) ∧
    (∀ M : Matrix (Fin da × Fin db) (Fin da × Fin db) ℂ,
      OPH.Locality.ptraceFst
          (∑ i, slotLeft (Fin db) (S.kraus i) * M *
            (slotLeft (Fin db) (S.kraus i))ᴴ) =
        OPH.Locality.ptraceFst M) ∧
    (∀ (M : Matrix (Fin da × Fin db) (Fin da × Fin db) ℂ)
        (B : Matrix (Fin db) (Fin db) ℂ),
      ((∑ i, slotLeft (Fin db) (S.kraus i) * M *
          (slotLeft (Fin db) (S.kraus i))ᴴ) * slotRight (Fin da) B).trace =
        (M * slotRight (Fin da) B).trace) := by
  haveI : Nontrivial (Matrix (Fin da × Fin db) (Fin da × Fin db) ℂ) := by
    have h0a : 0 < da := Nat.pos_of_ne_zero (NeZero.ne da)
    have h0b : 0 < db := Nat.pos_of_ne_zero (NeZero.ne db)
    refine ⟨0, 1, fun h => ?_⟩
    have h00 := Matrix.ext_iff.mpr h (⟨0, h0a⟩, ⟨0, h0b⟩) (⟨0, h0a⟩, ⟨0, h0b⟩)
    rw [Matrix.zero_apply, Matrix.one_apply_eq] at h00
    exact zero_ne_one h00
  -- Dichotomic observables push through the slot star embeddings.
  have hlift : ∀ {k : ℕ}
      (φ : Matrix (Fin k) (Fin k) ℂ →⋆ₐ[ℂ]
        Matrix (Fin da × Fin db) (Fin da × Fin db) ℂ)
      (P : Matrix (Fin k) (Fin k) ℂ), IsEvent P →
      IsSelfAdjoint (1 - 2 • φ P) ∧ (1 - 2 • φ P) * (1 - 2 • φ P) = 1 := by
    intro k φ P hP
    obtain ⟨hsa, hinv⟩ := hP.one_sub_two_smul_involution
    have hmap : φ (1 - 2 • P) = 1 - 2 • φ P := by
      rw [map_sub, map_one, map_nsmul]
    refine ⟨?_, ?_⟩
    · rw [← hmap]
      show star (φ (1 - 2 • P)) = φ (1 - 2 • P)
      rw [← map_star, hsa.star_eq]
    · rw [← hmap, ← map_mul, hinv, map_one]
  -- The no-signalling clause, proved once and shared with clause 3.
  have hns : ∀ M : Matrix (Fin da × Fin db) (Fin da × Fin db) ℂ,
      OPH.Locality.ptraceFst
          (∑ i, slotLeft (Fin db) (S.kraus i) * M *
            (slotLeft (Fin db) (S.kraus i))ᴴ) =
        OPH.Locality.ptraceFst M := by
    intro M
    have hlift_eq : ∀ i, slotLeft (Fin db) (S.kraus i) =
        OPH.Locality.liftLeft (S.kraus i) := fun i => slotLeft_eq_liftLeft _
    simp_rw [hlift_eq]
    exact OPH.Locality.ptraceFst_local_kraus Finset.univ S.kraus
      S.krausComplete M
  refine ⟨?_, hns, ?_⟩
  · -- Clause 1: the Tsirelson bound with derived commutation.
    obtain ⟨saA₀, invA₀⟩ := hlift (slotLeft (Fin db)) S.A₀ S.eventA₀
    obtain ⟨saA₁, invA₁⟩ := hlift (slotLeft (Fin db)) S.A₁ S.eventA₁
    obtain ⟨saB₀, invB₀⟩ := hlift (slotRight (Fin da)) S.B₀ S.eventB₀
    obtain ⟨saB₁, invB₁⟩ := hlift (slotRight (Fin da)) S.B₁ S.eventB₁
    have hc : ∀ (A : Matrix (Fin da) (Fin da) ℂ)
        (B : Matrix (Fin db) (Fin db) ℂ),
        slotLeft (Fin db) A * slotRight (Fin da) B =
          slotRight (Fin da) B * slotLeft (Fin db) A := fun A B =>
      (slot_commute ⟨A, rfl⟩ ⟨B, rfl⟩).eq
    exact tsirelson_bound _ _ _ _ saA₀ saA₁ saB₀ saB₁
      invA₀ invA₁ invB₀ invB₁
      (one_sub_two_smul_commute (hc S.A₀ S.B₀))
      (one_sub_two_smul_commute (hc S.A₀ S.B₁))
      (one_sub_two_smul_commute (hc S.A₁ S.B₀))
      (one_sub_two_smul_commute (hc S.A₁ S.B₁))
  · -- Clause 3: Born-pairing invariance through the partial trace.
    intro M B
    calc ((∑ i, slotLeft (Fin db) (S.kraus i) * M *
            (slotLeft (Fin db) (S.kraus i))ᴴ) * slotRight (Fin da) B).trace
        = (OPH.Locality.ptraceFst
            (∑ i, slotLeft (Fin db) (S.kraus i) * M *
              (slotLeft (Fin db) (S.kraus i))ᴴ) * B).trace :=
          trace_mul_slotRight _ B
      _ = (OPH.Locality.ptraceFst M * B).trace := by rw [hns M]
      _ = (M * slotRight (Fin da) B).trace := (trace_mul_slotRight M B).symm

end SlotLocality

end QuantumSurface

end EventAlgebra

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms EventAlgebra.QuantumSurface.trace_mul_slotRight
#print axioms EventAlgebra.QuantumSurface.slot_locality_receipts
