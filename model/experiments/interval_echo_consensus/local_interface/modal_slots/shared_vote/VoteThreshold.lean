import Mathlib.Tactic

namespace OPH.ModalQuorum

/-- For a binary-valued field on 96 packets, a tie goes to its lower value.
A packet carrying zero can switch the field in one edit precisely at 48. -/
theorem one_more_changes (high : ℕ) :
    (¬ 48 < high ∧ 48 < high + 1) ↔ high = 48 := by omega

/-- With all shared packets initially zero, each physical private packet edit
can retire at most one high vote of its owner. -/
theorem private_cost (privateVotes edits remaining : ℕ)
    (bounded : privateVotes ≤ remaining + edits) (edge : remaining = 48) :
    privateVotes - 48 ≤ edits := by omega

theorem shared_pair_cost (a b : ℕ) (ha : 71 ≤ 48 + a) (hb : 68 ≤ 48 + b) :
    43 ≤ a + b := by omega

end OPH.ModalQuorum

