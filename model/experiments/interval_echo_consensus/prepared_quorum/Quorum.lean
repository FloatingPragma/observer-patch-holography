import Mathlib.Data.Finset.Card
import Mathlib.Data.Fintype.Card
import Mathlib.Tactic

namespace OPH.PreparedQuorum

/-- Any two quorums at the registered threshold share a nonfaulty signer. -/
theorem honest_intersection {n f q : ℕ} (a b bad : Finset (Fin n))
    (ha : q ≤ a.card) (hb : q ≤ b.card) (hf : bad.card ≤ f)
    (threshold : n + f + 1 ≤ 2 * q) :
    ∃ i, i ∈ a ∧ i ∈ b ∧ i ∉ bad := by
  by_contra! hn
  have hsub : a ∩ b ⊆ bad := by
    intro i hi
    exact hn i (Finset.mem_inter.mp hi).1 (Finset.mem_inter.mp hi).2
  have hi := Finset.card_le_card hsub
  have hu : (a ∪ b).card ≤ n := by simpa using (a ∪ b).card_le_univ
  have he := Finset.card_union_add_card_inter a b
  omega

/-- After excluding all faulty signers, a quorum is still available. -/
theorem honest_available {n f q : ℕ} (bad : Finset (Fin n))
    (hf : bad.card ≤ f) (h : q + f ≤ n) :
    q ≤ ((Finset.univ : Finset (Fin n)) \ bad).card := by
  rw [Finset.card_sdiff]
  simp
  omega

#print axioms honest_intersection
#print axioms honest_available
end OPH.PreparedQuorum
