import Mathlib.Analysis.Complex.Basic
import Mathlib.Tactic.Linarith

/-!
# Directed winding of a binary written head

The executable attachment enumerates the two unit-contrast reversal sectors
of the registered seam and checks the signed endpoint current. A completed
clockwise transport on the closed guide takes the all-zero head word to
the all-one word. The operator below describes the branch heralded by a
successful signed hop at every contact.
-/

namespace OPH.Dynamics.WrittenHeadWinding

/-- The material and head differences checked on each native matrix entry
preserve the two endpoint charges, in the displayed orientation. -/
theorem endpoint_balance (a b c a' b' c' : ℤ)
    (ha : a' - a = c' - c) (hb : b - b' = c' - c) :
    a' - c' = a - c ∧ b' + c' = b + c := by
  constructor <;> linarith

/-- Coordinates of the two equal-position winding words. -/
def forward (x : ℂ × ℂ) : ℂ × ℂ := (0, x.1)

theorem twice_forward (x : ℂ × ℂ) : forward (forward x) = 0 := rfl

/-- Neither sign can be the scalar action of an oriented completed transport
on a nonzero state in this closure sector. -/
theorem no_signed_eigenvector (x : ℂ × ℂ)
    (h : forward x = x ∨ forward x = -x) : x = 0 := by
  rcases x with ⟨a, b⟩
  rcases h with h | h
  · have ha : 0 = a := congrArg Prod.fst h
    have hb : a = b := congrArg Prod.snd h
    simp [← hb, ← ha]
  · have ha : (0 : ℂ) = -a := congrArg Prod.fst h
    have hb : a = -b := congrArg Prod.snd h
    have hz : a = 0 := by simpa using congrArg Neg.neg ha.symm
    have hy : b = 0 := by simpa [hz] using congrArg Neg.neg hb.symm
    simp [hz, hy]

/-- For normalized amplitudes, the expectation of the successful oriented
winding has magnitude at most one half. Unobserved orthogonal components
only reduce the available squared norm. -/
theorem winding_overlap (a b : ℂ) (h : ‖a‖ ^ 2 + ‖b‖ ^ 2 ≤ 1) :
    ‖star b * a‖ ≤ 1 / 2 := by
  rw [norm_mul, norm_star]
  nlinarith [sq_nonneg (‖a‖ - ‖b‖)]

end OPH.Dynamics.WrittenHeadWinding
