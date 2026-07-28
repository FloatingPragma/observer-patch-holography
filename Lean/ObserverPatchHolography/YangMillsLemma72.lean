import Mathlib

/-!
# Yang–Mills finite repair-gap — Lemma 7.2 (keystone: uniform hidden fiber)

Formalisation of the **finite, real** §7.2 keystone of B. Müller, *Explaining
the Yang–Mills Mass Gap with Observer-Patch Repair Dynamics* (r1515):

> **Lemma 7.2.** A positive-semidefinite relaxation `D` on `L²(F)` over a
> uniform hidden fiber `F` (`|F| ≥ 2`) that commutes with the **full** symmetric
> group `S_F` and has kernel exactly the constants is a **strictly positive
> scalar** multiple of `I − E_F`, where `E_F` is the mean projection.

Proved concretely via the commutant of the permutation representation
(two-valued form, `commutant_perm_two_valued`) — no algebraically-closed-field
Schur needed.

## Honest scope (read this first)

This is an **implication**, nothing more. `hComm` (full-`S_F` symmetry) and
`hKer` (kernel = constants) are the paper's *physical modelling assumptions*
about the relaxation operator's structure — hypotheses here, discharged
**nowhere** in Lean. This proves "uniform-fiber structure ⇒ scalar relaxation";
it says **nothing** about `Δ_YM`. The continuum certificate (Müller's
**Assumption 9.2**: Schwinger convergence, reflection positivity,
Osterwalder–Schrader reconstruction, non-triviality) is the genuine open
problem and is **untouched, unassumed, unclaimed**.

## Integration

Extracted verbatim from the retired single-file artifact `RepairGap.lean`
(Part I; now consolidated into these modules — see git history) so the assembly
module `ObserverPatchHolography.YangMillsGap` can discharge its `lemma_7_2`
obligation by direct import:
```
exact ObserverPatchHolography.YangMillsLemma72.lemma_7_2 hF D hPSD hComm hKer
```
(the two `EF` definitions are token-identical, hence definitionally equal).

SCOPE, verbatim: Machine-checked: the finite representation gap Δ_rep ≥ c_* > 0
(Lemma 7.2 / Lemma 7.4 / Prop 8.1 / Thm 7.3 assembly) and the conditional
reduction "Assumption 9.2 + finite gap ⇒ Δ_YM ≥ c_*". Assumption 9.2 itself is
stated as an explicit hypothesis and is not touched.
-/

open scoped Matrix
open Matrix

namespace ObserverPatchHolography.YangMillsLemma72

variable {F : Type*} [Fintype F] [DecidableEq F]

/-- `E_F` : orthogonal projection onto constants (the mean operator), as a matrix.
    Every entry is `1/|F|`. -/
noncomputable def EF (F : Type*) [Fintype F] : Matrix F F ℝ :=
  (Fintype.card F : ℝ)⁻¹ • Matrix.of (fun _ _ => (1 : ℝ))

omit [DecidableEq F] in
@[simp] lemma EF_apply (i j : F) : EF F i j = (Fintype.card F : ℝ)⁻¹ := by
  simp [EF, Matrix.smul_apply, Matrix.of_apply, smul_eq_mul]

omit [Fintype F] in
/-- **2-transitivity of the full symmetric group** (explicit witness).
    Given two distinct source points and two distinct target points, there is a permutation
    taking the first pair to the second.  Built as a two-swap composition, so no appeal to
    the packaged `isMultiplyPretransitive` machinery is needed. -/
theorem exists_perm_maps_two {i j i' j' : F} (hij : i ≠ j) (hij' : i' ≠ j') :
    ∃ σ : Equiv.Perm F, σ i = i' ∧ σ j = j' := by
  refine ⟨Equiv.swap (Equiv.swap i i' j) j' * Equiv.swap i i', ?_, ?_⟩
  · -- σ i = i'
    rw [Equiv.Perm.mul_apply, Equiv.swap_apply_left]
    -- goal: swap (swap i i' j) j' i' = i'
    refine Equiv.swap_apply_of_ne_of_ne ?_ hij'
    -- goal: i' ≠ swap i i' j
    intro h
    apply hij
    have h2 : Equiv.swap i i' i' = Equiv.swap i i' (Equiv.swap i i' j) :=
      congrArg (Equiv.swap i i') h
    rwa [Equiv.swap_apply_right, Equiv.swap_apply_self] at h2
  · -- σ j = j'
    rw [Equiv.Perm.mul_apply, Equiv.swap_apply_left]

/-- **Conjugation identity.**  A matrix commuting with every permutation matrix is invariant
    under simultaneous permutation of its row and column indices:
    `M (σ i) (σ j) = M i j`.  Proof: push the commutation `σP·M = M·σP` through `· *ᵥ eₖ`
    and read off columns. -/
theorem perm_conj_invariant (M : Matrix F F ℝ)
    (hM : ∀ σ : Equiv.Perm F, Commute (σ.permMatrix ℝ) M) :
    ∀ (σ : Equiv.Perm F) (i j : F), M (σ i) (σ j) = M i j := by
  intro σ i j
  -- `Commute a b` is definitionally `a * b = b * a`
  have hcomm : σ.permMatrix ℝ * M = M * σ.permMatrix ℝ := hM σ
  -- apply `· *ᵥ Pi.single (σ j) 1` to both sides
  have h : (σ.permMatrix ℝ * M) *ᵥ Pi.single (σ j) (1 : ℝ)
         = (M * σ.permMatrix ℝ) *ᵥ Pi.single (σ j) (1 : ℝ) := by rw [hcomm]
  rw [← mulVec_mulVec, ← mulVec_mulVec, permMatrix_mulVec, permMatrix_mulVec] at h
  -- h : (M *ᵥ e) ∘ σ = M *ᵥ (e ∘ σ)      where  e = Pi.single (σ j) 1
  -- simplify  e ∘ σ = Pi.single j 1   (σ injective)
  have hcomp : (Pi.single (σ j) (1 : ℝ)) ∘ σ = Pi.single j (1 : ℝ) := by
    funext l
    simp only [Function.comp_apply, Pi.single_apply, σ.injective.eq_iff]
  rw [hcomp, mulVec_single_one, mulVec_single_one] at h
  -- h : M.col (σ j) ∘ σ = M.col j
  have := congrFun h i
  simpa only [Function.comp_apply, Matrix.col_apply] using this

/-- **Commutant of the permutation representation, concretely.**
    A real matrix commuting with every permutation matrix is two-valued:
    a common diagonal value `a` and a common off-diagonal value `b`.
    (Diagonal-constant from transitivity; off-diagonal-constant from 2-transitivity,
    which needs `2 ≤ |F|`.) -/
theorem commutant_perm_two_valued
    (hF : 2 ≤ Fintype.card F) (M : Matrix F F ℝ)
    (hM : ∀ σ : Equiv.Perm F, Commute (σ.permMatrix ℝ) M) :
    ∃ a b : ℝ, ∀ i j, M i j = if i = j then a else b := by
  have hconj := perm_conj_invariant M hM
  obtain ⟨p, q, hpq⟩ := Fintype.one_lt_card_iff.mp (by omega : 1 < Fintype.card F)
  refine ⟨M p p, M p q, fun i j => ?_⟩
  by_cases h : i = j
  · subst h
    rw [if_pos rfl]
    have := hconj (Equiv.swap p i) p p
    rwa [Equiv.swap_apply_left] at this
  · rw [if_neg h]
    obtain ⟨σ, hσp, hσq⟩ := exists_perm_maps_two hpq h
    have := hconj σ p q
    rwa [hσp, hσq] at this

/-- **LEMMA 7.2 (keystone).**  Let `D` be positive semidefinite (this already bundles
    Hermitian), commuting with the full symmetric group's permutation action, with kernel
    exactly the constants.  Then `D = c_F · (I − E_F)` for some `c_F > 0`.

    Roles of the hypotheses (all genuinely used, none vacuous):
    * `hComm`  ⇒ `D` is two-valued: `D i j = if i=j then a else b`   (Schur / commutant step);
    * `hKer`   ⇒ constants ∈ ker gives `a + (|F|−1)·b = 0`, and ker = constants gives `b ≠ 0`;
    * `hPSD`   ⇒ the quadratic form on `e_p − e_q` equals `2(a−b) ≥ 0`, fixing the sign.
    Setting `c_F := a − b` yields `c_F > 0` and `D = c_F · (1 − E_F)` entrywise.

    HONEST: `hComm` and `hKer` are the paper's physical modelling assumptions, hypotheses
    here — this proves "uniform-fiber structure ⇒ scalar relaxation", nothing about `Δ_YM`. -/
theorem lemma_7_2
    (hF : 2 ≤ Fintype.card F) (D : Matrix F F ℝ)
    (hPSD  : D.PosSemidef)
    (hComm : ∀ σ : Equiv.Perm F, Commute (σ.permMatrix ℝ) D)
    (hKer  : ∀ v : F → ℝ, D.mulVec v = 0 ↔ ∃ c : ℝ, v = fun _ => c) :
    ∃ cF : ℝ, 0 < cF ∧ D = cF • (1 - EF F) := by
  -- Step 0: two distinct fiber points and the two-valued form of D.
  obtain ⟨p, q, hpq⟩ := Fintype.one_lt_card_iff.mp (by omega : 1 < Fintype.card F)
  obtain ⟨a, b, hab⟩ := commutant_perm_two_valued hF D hComm
  have hN0 : (Fintype.card F : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
  -- Step 1: constants in the kernel  ⇒  a + (N-1)·b = 0.
  have hconst_ker : D.mulVec (fun _ => (1 : ℝ)) = 0 := (hKer _).mpr ⟨1, rfl⟩
  have hrow : a + ((Fintype.card F : ℝ) - 1) * b = 0 := by
    have hp0 : ∑ j, D p j = 0 := by
      have h := congrFun hconst_ker p
      simpa [Matrix.mulVec, dotProduct, mul_one] using h
    have hsum : ∑ j, D p j = a + ((Fintype.card F : ℝ) - 1) * b := by
      have hshift : ∀ j, D p j = b + (if p = j then a - b else 0) := by
        intro j; rw [hab p j]; split_ifs <;> ring
      simp_rw [hshift]
      rw [Finset.sum_add_distrib, Finset.sum_ite_eq]
      simp only [Finset.mem_univ, if_true, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
      ring
    rw [hsum] at hp0; linarith
  -- Step 2: positive-semidefiniteness on  e_p − e_q  gives  0 ≤ a − b.
  have hquad : (0 : ℝ) ≤ a - b := by
    have hpos := hPSD.dotProduct_mulVec_nonneg ((Pi.single p 1 : F → ℝ) - (Pi.single q 1 : F → ℝ))
    have hstar : star ((Pi.single p 1 : F → ℝ) - (Pi.single q 1 : F → ℝ))
               = (Pi.single p 1 : F → ℝ) - (Pi.single q 1 : F → ℝ) := by
      funext k; simp [Pi.star_apply, star_trivial]
    rw [hstar] at hpos
    have hval : ((Pi.single p 1 : F → ℝ) - (Pi.single q 1 : F → ℝ)) ⬝ᵥ
                (D *ᵥ ((Pi.single p 1 : F → ℝ) - (Pi.single q 1 : F → ℝ))) = 2 * (a - b) := by
      rw [mulVec_sub, mulVec_single_one, mulVec_single_one,
          sub_dotProduct, dotProduct_sub, dotProduct_sub,
          single_dotProduct, single_dotProduct, single_dotProduct, single_dotProduct]
      simp only [Matrix.col_apply, one_mul]
      rw [hab p p, hab p q, hab q p, hab q q,
          if_pos rfl, if_pos rfl, if_neg hpq, if_neg (Ne.symm hpq)]
      ring
    rw [hval] at hpos; linarith
  -- Step 3: ker = constants  ⇒  b ≠ 0  (else D = 0, contradicting kernel exactness).
  have hbne : b ≠ 0 := by
    intro hb0
    have ha0 : a = 0 := by rw [hb0] at hrow; simpa using hrow
    have hD0 : D.mulVec ((Pi.single p 1 : F → ℝ)) = 0 := by
      have hDeq : D = 0 := by ext i j; rw [hab i j, hb0, ha0]; simp
      rw [hDeq]; exact Matrix.zero_mulVec _
    obtain ⟨c, hc⟩ := (hKer ((Pi.single p 1 : F → ℝ))).mp hD0
    have h1 : ((Pi.single p 1 : F → ℝ)) p = c := by rw [hc]
    have h2 : ((Pi.single p 1 : F → ℝ)) q = c := by rw [hc]
    rw [Pi.single_eq_same] at h1
    rw [Pi.single_eq_of_ne (Ne.symm hpq)] at h2
    exact one_ne_zero (h1.trans h2.symm)
  -- Step 4: assemble.  cF := a − b > 0, and D = cF·(1 − E_F) entrywise.
  refine ⟨a - b, ?_, ?_⟩
  · -- 0 < a − b
    have hne : a - b ≠ 0 := by
      have hcf : a - b = (Fintype.card F : ℝ) * (-b) := by linear_combination hrow
      rw [hcf]; exact mul_ne_zero hN0 (neg_ne_zero.mpr hbne)
    exact lt_of_le_of_ne hquad (Ne.symm hne)
  · -- D = (a − b) • (1 − E_F)
    have ha_eq : a = (1 - (Fintype.card F : ℝ)) * b := by linear_combination hrow
    ext i j
    rw [hab i j]
    simp only [Matrix.smul_apply, Matrix.sub_apply, Matrix.one_apply, EF_apply, smul_eq_mul]
    by_cases h : i = j
    · rw [if_pos h, if_pos h, ha_eq]; field_simp; ring
    · rw [if_neg h, if_neg h, ha_eq]; field_simp; ring

/-! ## Axiom self-audit (build-log visible)

Expected report for every theorem below: exactly
`[propext, Classical.choice, Quot.sound]` — no `sorryAx`, no project axiom. -/

#print axioms exists_perm_maps_two
#print axioms perm_conj_invariant
#print axioms commutant_perm_two_valued
#print axioms lemma_7_2

end ObserverPatchHolography.YangMillsLemma72
