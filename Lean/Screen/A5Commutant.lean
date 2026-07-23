import Mathlib
import A5PortModule
import PortFrameGram

namespace OPH.A5Commutant

open OPH.A5PortModule (P60 pAct)
open OPH.PortFrameGram (adj antipode adj_irrefl antipode_not_adjacent
  antipode_fixed_point_free)

/-! # Equivariant endomorphisms of the port module: the commutant, exactly

Issue #568.  The listed icosahedral action partitions ordered port pairs
into four orbits: diagonal, adjacent, distance-two, and antipodal
(`pair_transport`, kernel-checked).  Consequently the algebra of
`ℚ`-linear endomorphisms commuting with every listed rotation is exactly
the four-dimensional span of the four orbital matrices — identity,
adjacency, distance-two, antipode (`commutant_decomposition`,
`orbitals_independent`, `orbitals_equivariant`).  The matrix-entry
invariance used throughout is proven equivalent to commutation with the
linear action (`equivariant_iff_commutes`), so no postulate hides in the
executable definition.

Dimension four matches the rational isotypic structure
`P12 ≅ 1 ⊕ 5 ⊕ (3 ⊕ 3')` over `ℚ`: two rational irreducibles with
endomorphism field `ℚ` and one with endomorphism field `ℚ(√5)`, giving
`1 + 1 + 2 = 4` (`Screen/A5CharacterField.lean` carries the Galois pairing
of `3` and `3'`).

BOUNDARY.  This file is coefficient algebra on the explicit sixty-element
list; physical currents, global descent, and matter realization are
outside it, per the lane boundaries recorded in the companion modules. -/

/-- Matrix-entry equivariance under every listed rotation. -/
def Equivariant (M : Fin 12 → Fin 12 → ℚ) : Prop :=
  ∀ g ∈ P60, ∀ i j : Fin 12, M (g i) (g j) = M i j

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Boolean form: every listed rotation preserves port adjacency. -/
theorem P60_preserves_adj_bool :
    (P60.all fun g => (List.finRange 12).all fun i =>
      (List.finRange 12).all fun j => adj (g i) (g j) == adj i j) = true := by
  decide

/-- Every listed rotation preserves port adjacency, at the `Fin 12`
level. -/
theorem P60_preserves_adj :
    ∀ g ∈ P60, ∀ i j : Fin 12, adj (g i) (g j) = adj i j := by
  intro g hg i j
  have h := List.all_eq_true.mp P60_preserves_adj_bool g hg
  have hi := List.all_eq_true.mp h i (List.mem_finRange i)
  have hj := List.all_eq_true.mp hi j (List.mem_finRange j)
  exact eq_of_beq hj

set_option maxHeartbeats 2000000 in
set_option maxRecDepth 16384 in
/-- Boolean form: every listed rotation commutes with the antipode. -/
theorem P60_antipode_bool :
    (P60.all fun g => (List.finRange 12).all fun i =>
      g (antipode i) == antipode (g i)) = true := by
  decide

/-- Every listed rotation commutes with the antipode, at the `Fin 12`
level. -/
theorem P60_antipode_equivariant :
    ∀ g ∈ P60, ∀ i : Fin 12, g (antipode i) = antipode (g i) := by
  intro g hg i
  have h := List.all_eq_true.mp P60_antipode_bool g hg
  exact eq_of_beq (List.all_eq_true.mp h i (List.mem_finRange i))

/-- Canonical representative pair of an ordered port pair: `(0, 0)` on the
diagonal, `(0, 1)` for adjacent pairs, `(0, 11)` for antipodal pairs,
`(0, 5)` for distance-two pairs. -/
def rep (i j : Fin 12) : Fin 12 :=
  if i = j then 0 else if adj i j then 1 else if j = antipode i then 11
  else 5

set_option maxHeartbeats 8000000 in
set_option maxRecDepth 16384 in
/-- Boolean form: every ordered pair is carried by a listed rotation to
its canonical representative. -/
theorem pair_transport_bool :
    ((List.finRange 12).all fun i => (List.finRange 12).all fun j =>
      P60.any fun g => g i == 0 && g j == rep i j) = true := by
  decide

/-- Pair transport: every ordered pair is carried by a listed rotation to
its canonical representative — the four orbit classes, kernel-checked. -/
theorem pair_transport :
    ∀ i j : Fin 12, ∃ g ∈ P60, g i = 0 ∧ g j = rep i j := by
  intro i j
  have h := List.all_eq_true.mp pair_transport_bool i (List.mem_finRange i)
  have hj := List.all_eq_true.mp h j (List.mem_finRange j)
  obtain ⟨g, hg, hgb⟩ := List.any_eq_true.mp hj
  rw [Bool.and_eq_true] at hgb
  exact ⟨g, hg, eq_of_beq hgb.1, eq_of_beq hgb.2⟩

/-- The identity orbital. -/
def dOrb (i j : Fin 12) : ℚ := if i = j then 1 else 0

/-- The adjacency orbital. -/
def aOrb (i j : Fin 12) : ℚ := if adj i j then 1 else 0

/-- The antipode orbital. -/
def pOrb (i j : Fin 12) : ℚ := if j = antipode i then 1 else 0

/-- The distance-two orbital. -/
def nOrb (i j : Fin 12) : ℚ :=
  if i ≠ j ∧ adj i j = false ∧ j ≠ antipode i then 1 else 0

theorem equivariant_dOrb : Equivariant dOrb := by
  intro g _ i j
  unfold dOrb
  by_cases h : i = j
  · rw [if_pos (by rw [h]), if_pos h]
  · rw [if_neg (fun he => h (g.injective he)), if_neg h]

theorem equivariant_aOrb : Equivariant aOrb := by
  intro g hg i j
  unfold aOrb
  rw [P60_preserves_adj g hg i j]

theorem equivariant_pOrb : Equivariant pOrb := by
  intro g hg i j
  have hiff : g j = antipode (g i) ↔ j = antipode i := by
    rw [← P60_antipode_equivariant g hg i]
    exact ⟨fun h => g.injective h, fun h => by rw [h]⟩
  unfold pOrb
  by_cases h : j = antipode i
  · rw [if_pos (hiff.mpr h), if_pos h]
  · rw [if_neg (fun hc => h (hiff.mp hc)), if_neg h]

theorem equivariant_nOrb : Equivariant nOrb := by
  intro g hg i j
  have h1 : g i = g j ↔ i = j := ⟨fun h => g.injective h, fun h => by rw [h]⟩
  have h2 := P60_preserves_adj g hg i j
  have h3 : g j = antipode (g i) ↔ j = antipode i := by
    rw [← P60_antipode_equivariant g hg i]
    exact ⟨fun h => g.injective h, fun h => by rw [h]⟩
  unfold nOrb
  by_cases hc : i ≠ j ∧ adj i j = false ∧ j ≠ antipode i
  · rw [if_pos ⟨fun he => hc.1 (h1.mp he), by rw [h2]; exact hc.2.1,
      fun he => hc.2.2 (h3.mp he)⟩, if_pos hc]
  · rw [if_neg (fun hgc => hc ⟨fun he => hgc.1 (h1.mpr he),
      by rw [← h2]; exact hgc.2.1, fun he => hgc.2.2 (h3.mpr he)⟩),
      if_neg hc]

/-- Every equivariant matrix is the four-parameter orbital combination,
with coefficients read off at the four canonical entries. -/
theorem commutant_decomposition (M : Fin 12 → Fin 12 → ℚ)
    (hM : Equivariant M) :
    ∀ i j : Fin 12,
      M i j = M 0 0 * dOrb i j + M 0 1 * aOrb i j +
        M 0 5 * nOrb i j + M 0 11 * pOrb i j := by
  intro i j
  obtain ⟨g, hg, hgi, hgj⟩ := pair_transport i j
  have hMij : M i j = M 0 (rep i j) := by
    rw [← hM g hg i j, hgi, hgj]
  by_cases h1 : i = j
  · subst h1
    have ha := adj_irrefl i
    have hp : i ≠ antipode i :=
      fun h => antipode_fixed_point_free i h.symm
    rw [hMij]
    simp [rep, dOrb, aOrb, nOrb, pOrb, ha, hp]
  · by_cases h2 : adj i j = true
    · have hp : j ≠ antipode i := by
        intro h
        rw [h] at h2
        rw [antipode_not_adjacent i] at h2
        exact Bool.false_ne_true h2
      rw [hMij]
      simp [rep, dOrb, aOrb, nOrb, pOrb, h1, h2, hp]
    · by_cases h3 : j = antipode i
      · subst h3
        have ha := antipode_not_adjacent i
        have hne : i ≠ antipode i :=
          fun h => antipode_fixed_point_free i h.symm
        rw [hMij]
        simp [rep, dOrb, aOrb, nOrb, pOrb, hne, ha]
      · rw [hMij]
        simp [rep, dOrb, aOrb, nOrb, pOrb, h1, h2, h3]

/-- The four orbitals are linearly independent: the commutant has
dimension exactly four. -/
theorem orbitals_independent (a b c d : ℚ)
    (h : ∀ i j : Fin 12,
      a * dOrb i j + b * aOrb i j + c * nOrb i j + d * pOrb i j = 0) :
    a = 0 ∧ b = 0 ∧ c = 0 ∧ d = 0 := by
  have e1 : adj (0 : Fin 12) 0 = false := by decide
  have e2 : adj (0 : Fin 12) 1 = true := by decide
  have e3 : adj (0 : Fin 12) 5 = false := by decide
  have e4 : adj (0 : Fin 12) 11 = false := by decide
  have e5 : antipode (0 : Fin 12) = 11 := by decide
  have f1 : (0 : Fin 12) ≠ 1 := by decide
  have f2 : (0 : Fin 12) ≠ 5 := by decide
  have f3 : (0 : Fin 12) ≠ 11 := by decide
  have f4 : (1 : Fin 12) ≠ 11 := by decide
  have f5 : (5 : Fin 12) ≠ 11 := by decide
  have h00 := h 0 0
  have h01 := h 0 1
  have h05 := h 0 5
  have h011 := h 0 11
  simp [dOrb, aOrb, nOrb, pOrb, e1, e2, e3, e4, e5, f1, f2, f3, f4,
    f5] at h00 h01 h05 h011
  exact ⟨h00, h01, h05, h011⟩

/-- Matrix action as a linear endomorphism of the port module. -/
def matL (M : Fin 12 → Fin 12 → ℚ) : (Fin 12 → ℚ) →ₗ[ℚ] (Fin 12 → ℚ) where
  toFun v i := ∑ j, M i j * v j
  map_add' v w := by
    funext i
    show ∑ j, M i j * (v j + w j) =
      (∑ j, M i j * v j) + ∑ j, M i j * w j
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun j _ => by ring
  map_smul' c v := by
    funext i
    show ∑ j, M i j * (c * v j) = c * ∑ j, M i j * v j
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun j _ => by ring

/-- Entry equivariance is exactly commutation with the listed linear
action: nothing is hidden in the executable definition. -/
theorem equivariant_iff_commutes (M : Fin 12 → Fin 12 → ℚ) :
    Equivariant M ↔
      ∀ g ∈ P60, ∀ v : Fin 12 → ℚ,
        matL M (pAct g v) = pAct g (matL M v) := by
  constructor
  · intro hM g hg v
    have hinv : ∀ x : Fin 12, g⁻¹ (g x) = x :=
      fun x => Equiv.symm_apply_apply g x
    have hinv' : ∀ x : Fin 12, g (g⁻¹ x) = x :=
      fun x => Equiv.apply_symm_apply g x
    funext i
    show ∑ j, M i j * v (g⁻¹ j) = ∑ j, M (g⁻¹ i) j * v j
    calc ∑ j, M i j * v (g⁻¹ j)
        = ∑ j, M i (g j) * v (g⁻¹ (g j)) :=
          (Equiv.sum_comp g fun j => M i j * v (g⁻¹ j)).symm
      _ = ∑ j, M (g⁻¹ i) j * v j := by
          refine Finset.sum_congr rfl fun j _ => ?_
          rw [hinv j]
          have := hM g hg (g⁻¹ i) j
          rw [hinv' i] at this
          rw [this]
  · intro h g hg i j
    have hinv : ∀ x : Fin 12, g⁻¹ (g x) = x :=
      fun x => Equiv.symm_apply_apply g x
    have hinv' : ∀ x : Fin 12, g (g⁻¹ x) = x :=
      fun x => Equiv.apply_symm_apply g x
    have hbasis := congrFun (h g hg fun x => if x = j then (1 : ℚ) else 0)
      (g i)
    have hL : (matL M (pAct g fun x => if x = j then (1 : ℚ) else 0)) (g i)
        = M (g i) (g j) := by
      show (∑ k, M (g i) k * if g⁻¹ k = j then (1 : ℚ) else 0)
        = M (g i) (g j)
      have hcond : ∀ k : Fin 12,
          (M (g i) k * if g⁻¹ k = j then (1 : ℚ) else 0)
            = if k = g j then M (g i) k else 0 := by
        intro k
        by_cases hk : k = g j
        · rw [if_pos (by rw [hk, hinv]), if_pos hk, mul_one]
        · rw [if_neg (fun hc => hk (by rw [← hc, hinv'])), if_neg hk,
            mul_zero]
      rw [Finset.sum_congr rfl fun k _ => hcond k, Finset.sum_ite_eq']
      simp
    have hR : (pAct g (matL M fun x => if x = j then (1 : ℚ) else 0)) (g i)
        = M i j := by
      show (∑ k, M (g⁻¹ (g i)) k * if k = j then (1 : ℚ) else 0) = M i j
      rw [hinv i]
      have hcond : ∀ k : Fin 12,
          (M i k * if k = j then (1 : ℚ) else 0)
            = if k = j then M i k else 0 := by
        intro k
        by_cases hk : k = j
        · rw [if_pos hk, if_pos hk, mul_one]
        · rw [if_neg hk, if_neg hk, mul_zero]
      rw [Finset.sum_congr rfl fun k _ => hcond k, Finset.sum_ite_eq']
      simp
    rw [hL, hR] at hbasis
    exact hbasis

end OPH.A5Commutant

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.A5Commutant.P60_preserves_adj
#print axioms OPH.A5Commutant.P60_antipode_equivariant
#print axioms OPH.A5Commutant.pair_transport
#print axioms OPH.A5Commutant.commutant_decomposition
#print axioms OPH.A5Commutant.orbitals_independent
#print axioms OPH.A5Commutant.equivariant_iff_commutes
