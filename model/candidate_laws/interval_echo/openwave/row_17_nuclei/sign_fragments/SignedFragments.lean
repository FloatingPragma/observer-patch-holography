import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Analysis.Complex.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

namespace IntervalEcho.SignedFragments

abbrev T := Fin 3

private def pos (a : T) : T := max a 1
private def neg (a : T) : T := min a 1
private def idEntry (a b c d : T) : ℤ := if a = c ∧ b = d then 1 else 0
private def exchange (a b c d : T) : ℤ := idEntry a b c d - idEntry b a c d
private def contrast (a b c d : T) : ℤ :=
  idEntry a b c d * ((a.val : ℤ) - b.val)^2
private def unlike (a b : T) : Bool := (a == 0 && b == 2) || (a == 2 && b == 0)
private def excessEntry (h : T → T → T → T → ℤ) (a b c d : T) : ℤ :=
  h (pos a) (pos b) (pos c) (pos d) * idEntry (neg a) (neg b) (neg c) (neg d)
private def vacancyEntry (h : T → T → T → T → ℤ) (a b c d : T) : ℤ :=
  h (neg a) (neg b) (neg c) (neg d) * idEntry (pos a) (pos b) (pos c) (pos d)

/-- Exchange is written in twice-energy units; contrast is its unscaled diagonal. -/
theorem all_contact_entries : ∀ a b c d : T,
    exchange a b c d - excessEntry exchange a b c d - vacancyEntry exchange a b c d
       = (if unlike a b then -idEntry a b c d - idEntry b a c d else 0) ∧
    contrast a b c d - excessEntry contrast a b c d - vacancyEntry contrast a b c d
       = (if unlike a b then 2 * idEntry a b c d else 0) := by decide +kernel

theorem embedded_word_recovery : ∀ a : T, (pos a).val + (neg a).val = a.val + 1 := by
  decide +kernel

theorem referenced_site : ∀ a : T, (a.val : ℤ)^2 - 1 =
  ((pos a).val : ℤ)^2 - 1 + (((neg a).val : ℤ)^2 - 1) := by decide +kernel

/-- Apply once to real and once to imaginary parts of each unlike swap pair. -/
theorem pair_energy (k x y : ℝ) :
    (2*k-1/2)*(x^2+y^2)-x*y =
    (2*k-1)*(x^2+y^2)+(x-y)^2/2 := by ring

theorem joined_dominates {I : Type*} [Fintype I] (k : ℝ) (hk : 1/2 ≤ k)
    (x y : I → ℝ) :
    0 ≤ ∑ i, ((2*k-1/2)*((x i)^2+(y i)^2)-(x i)*(y i)) := by
  apply Finset.sum_nonneg
  intro i _
  rw [pair_energy]
  have h : 0 ≤ 2*k-1 := by linarith
  exact add_nonneg (mul_nonneg h (add_nonneg (sq_nonneg _) (sq_nonneg _)))
    (div_nonneg (sq_nonneg _) (by norm_num))

/-- A normalized collision-free competitor has overlap weight at most 1-c
    with the separated ground. Its orthogonal weight costs at least the gap. -/
theorem collision_penalty {ground joined split gap c weight : ℝ}
    (g : 0 ≤ gap) (outside : weight ≤ 1-c)
    (spectral : ground + gap*(1-weight) ≤ split) (dominates : split ≤ joined) :
    ground + gap*c ≤ joined := by
  have h : c ≤ 1-weight := by linarith
  linarith [mul_le_mul_of_nonneg_left h g]

#print axioms all_contact_entries
#print axioms embedded_word_recovery
#print axioms referenced_site
#print axioms pair_energy
#print axioms joined_dominates
#print axioms collision_penalty

end IntervalEcho.SignedFragments
