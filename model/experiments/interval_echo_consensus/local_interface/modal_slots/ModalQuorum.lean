import Mathlib.Tactic

namespace OPH.ModalQuorum

/-- Field histogram when every private packet reports its carrier value.
The remaining votes can depend on the whole joint state. -/
def fieldVotes (privateCount : ℕ) (target : ℕ) (others : ℕ → ℕ) (value : ℕ) : ℕ :=
  (if value = target then privateCount else 0) + others value

/-- Every competitor loses, including when the remaining packets are
chosen adversarially and correlate across observers. -/
theorem unique_mode (privateCount shared target : ℕ) (others : ℕ → ℕ)
    (hmajor : shared < privateCount) (hcap : ∀ x, others x ≤ shared)
    (x : ℕ) (hx : x ≠ target) :
    fieldVotes privateCount target others x < fieldVotes privateCount target others target := by
  simp only [fieldVotes, hx, if_false, zero_add, if_true]
  exact lt_of_le_of_lt (hcap x) (lt_of_lt_of_le hmajor (Nat.le_add_right _ _))

/-- On each full 96-node support all fields of every supplied triple have a
unique literal mode. The arbitrary other votes need only obey their actual
overlap bound. -/
theorem native_three_supports (a b c : ℕ) (u v w : ℕ → ℕ)
    (hu : ∀ x, u x ≤ 25) (hv : ∀ x, v x ≤ 28) (hw : ∀ x, w x ≤ 3) :
    (∀ x, x ≠ a → fieldVotes 71 a u x < fieldVotes 71 a u a) ∧
    (∀ x, x ≠ b → fieldVotes 68 b v x < fieldVotes 68 b v b) ∧
    (∀ x, x ≠ c → fieldVotes 93 c w x < fieldVotes 93 c w c) :=
  ⟨unique_mode 71 25 a u (by decide) hu,
   unique_mode 68 28 b v (by decide) hv,
   unique_mode 93 3 c w (by decide) hw⟩

end OPH.ModalQuorum
