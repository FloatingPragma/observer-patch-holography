import ObserverPatchHolography.IcosahedralAntibridge
import Mathlib.GroupTheory.SpecificGroups.Alternating
import Mathlib.GroupTheory.GroupAction.Quotient

/-!
# Icosahedral orbit–stabilizer receipts for the `A5` vertex action

Companion to `IcosahedralAntibridge`.  That file builds the concrete
transitive `A5`-action on the twelve-point homogeneous space
`DefectSpace = A5 ⧸ H5` and proves the anti-bridge no-go.  This file adds
the finite **orbit / stabilizer / faithfulness** computations that the
`A5` formalization lane requires
(`FloatingPragma/observer-patch-holography#568`, acceptance item
"finite orbit, stabilizer, character, rank, and kernel computations are
checked by the prover"):

* `stabilizer_coe_one_eq_H5` — the stabilizer of the base vertex `⟦1⟧`
  is exactly `H5` (so vertex stabilizers are the order-5 rotation
  subgroups), and `card_stabilizer_coe_one = 5`.
* `orbit_stabilizer_defectSpace` — the orbit–stabilizer identity on the
  twelve vertices: `|A5| = 12 * 5`, i.e. `60 = 12 * 5`.
* `defectSpace_action_faithful` — the kernel of the action is trivial
  (`H5.normalCore = ⊥`), because `A5` is simple and `H5` is a proper
  subgroup.  This is the "kernel computation" the lane asks for.

Claim discipline: these are finite facts about the abstract icosahedral
model `A5 ⧸ H5`.  They do not by themselves establish any physical
current, compact Lie-algebra classification, cover, or matter claim; the
`#568` lane keeps those as explicit theorem boundaries and open items.
-/

namespace OPH

open MulAction

/-- The stabilizer of the base vertex `⟦1⟧ : A5 ⧸ H5` is exactly `H5`.
    Physically: the vertex stabilizers are the order-5 rotation
    subgroups. -/
theorem stabilizer_coe_one_eq_H5 :
    stabilizer A5 ((1 : A5) : DefectSpace) = H5 := by
  simp [stabilizer_quotient]

/-- The base vertex stabilizer has order 5. -/
theorem card_stabilizer_coe_one :
    Nat.card (stabilizer A5 ((1 : A5) : DefectSpace)) = 5 := by
  rw [stabilizer_coe_one_eq_H5, card_H5]

/-- **Orbit–stabilizer on the twelve vertices.**  The order of `A5`
    factors as `(number of vertices) * (order of a vertex stabilizer)`:
    `60 = 12 * 5`.  Uses transitivity (orbit `= univ`, size 12) and the
    order-5 vertex stabilizer. -/
theorem orbit_stabilizer_defectSpace :
    Nat.card A5 = Nat.card DefectSpace * Nat.card (stabilizer A5 ((1 : A5) : DefectSpace)) := by
  rw [card_defectSpace, card_stabilizer_coe_one, nat_card_A5]

/-- **The vertex action is faithful.**  The kernel of the `A5`-action on
    the twelve vertices is the normal core of `H5`, which is trivial
    because `A5` is simple and `H5 ≠ ⊤`.  Hence `A5` embeds into the
    permutations of the twelve icosahedral vertices. -/
theorem defectSpace_action_faithful : H5.normalCore = ⊥ := by
  rcases (IsSimpleGroup.eq_bot_or_eq_top_of_normal (H5.normalCore)
      (Subgroup.normalCore_normal H5)) with h | h
  · exact h
  · exfalso
    have hle : H5.normalCore ≤ H5 := Subgroup.normalCore_le H5
    rw [h, top_le_iff] at hle
    have h60 : Nat.card H5 = 60 := by
      rw [hle]; simpa using (Subgroup.card_top (G := A5)).trans nat_card_A5
    rw [card_H5] at h60
    omega

/-! ## Axiom audit

Expected footprint for every entry: `[propext, Classical.choice,
Quot.sound]`.  No `sorry`, no `native_decide`, no extra axioms. -/

#print axioms stabilizer_coe_one_eq_H5
#print axioms card_stabilizer_coe_one
#print axioms orbit_stabilizer_defectSpace
#print axioms defectSpace_action_faithful

end OPH
