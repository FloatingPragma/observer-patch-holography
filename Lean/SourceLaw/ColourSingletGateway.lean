import Mathlib.Data.Real.Basic
import Mathlib.Logic.Relation
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Currents at colour-singlet gateways

The cell projector and the gateway incidences are evaluated by the single-load
experiment. These lemmas give the algebraic steps from its gateway Gram matrices
to vanishing probability current. Spatial real and imaginary amplitudes obey
the same real Gram equations.
-/

namespace SourceLaw.ColourSingletGateway

/-- A three-column gateway Gram matrix with diagonal greater than two and
off-diagonal signs has trivial kernel over the reals. Apply to both quadratures. -/
theorem independent_gateway (s a b c : ℝ) (u v w : Bool) (scale : 2 < s)
    (first : s * a + (if u then b else -b) + (if v then c else -c) = 0)
    (second : (if u then a else -a) + s * b + (if w then c else -c) = 0)
    (third : (if v then a else -a) + (if w then b else -b) + s * c = 0) :
    a = 0 ∧ b = 0 ∧ c = 0 := by
  have sum : a ^ 2 + b ^ 2 + c ^ 2 = 0 := by
    by_contra nonzero
    have positive : 0 < a ^ 2 + b ^ 2 + c ^ 2 := lt_of_le_of_ne
      (add_nonneg (add_nonneg (sq_nonneg a) (sq_nonneg b)) (sq_nonneg c)) (Ne.symm nonzero)
    have excess := mul_pos (sub_pos.mpr scale) positive
    cases u <;> cases v <;> cases w <;> simp at first second third <;>
      nlinarith only [excess, congrArg (a * ·) first, congrArg (b * ·) second,
        congrArg (c * ·) third, sq_nonneg (a+b), sq_nonneg (a-b), sq_nonneg (a+c),
        sq_nonneg (a-c), sq_nonneg (b+c), sq_nonneg (b-c)]
  constructor <;> [skip; constructor] <;>
    nlinarith [sq_nonneg a, sq_nonneg b, sq_nonneg c]

/-- Real proportionality of the endpoint amplitudes makes the imaginary
hopping current vanish, including where either amplitude vanishes. -/
theorem constrained_current (a b c d α : ℝ) (real_part : α * a = c)
    (imag_part : α * b = d) : a * d - b * c = 0 := by
  rw [← real_part, ← imag_part]
  ring

/-- A forced zero propagates along any finite sequence of equalities up to
sign. This certifies connected components touching a full-rank endpoint. -/
theorem linked_zero {E : Type*} (link : E → E → Prop) (coef : E → ℝ)
    (carry : ∀ i j, link i j → coef i = coef j ∨ coef i = -coef j)
    {i j : E} (path : Relation.ReflTransGen link i j) (forced : coef i = 0) :
    coef j = 0 := by
  induction path with
  | refl => exact forced
  | tail h edge ih =>
    rcases carry _ _ edge with eq | eq <;> rw [ih] at eq <;> linarith

end SourceLaw.ColourSingletGateway
