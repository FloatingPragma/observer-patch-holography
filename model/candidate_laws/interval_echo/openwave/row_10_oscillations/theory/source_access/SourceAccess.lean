import Mathlib.Data.Complex.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Algebra.BigOperators.Ring.Finset

/-! Which entries of an input can contribute to a fixed complete history of sharp records? -/
namespace OPH.NeutralSourceAccess

open scoped BigOperators

variable {I T : Type*} [Fintype I]
abbrev Kernel (I : Type*) := I → I → ℂ

def Preserves (tag : I → T) (k : Kernel I) : Prop :=
  ∀ i j, tag i ≠ tag j → k i j = 0

noncomputable def comp (k l : Kernel I) : Kernel I := fun i j => ∑ v, k i v * l v j

theorem comp_preserves (tag : I → T) (k l : Kernel I)
    (hk : Preserves tag k) (hl : Preserves tag l) : Preserves tag (comp k l) := by
  intro i j hij
  apply Finset.sum_eq_zero
  intro v _
  by_cases h : tag i = tag v
  · rw [hl v j (fun e => hij (h.trans e)), mul_zero]
  · rw [hk i v h, zero_mul]

noncomputable def identityK [DecidableEq I] : Kernel I := fun i j => if i = j then 1 else 0

omit [Fintype I] in
theorem sharp_preserves [DecidableEq I] (tag : I → T) (weight : I → ℂ) :
    Preserves tag (fun i j => if i = j then weight i else 0) := by
  intro i j hij
  exact if_neg (fun h => hij (congrArg tag h))

noncomputable def history [DecidableEq I] : List (Kernel I) → Kernel I
  | [] => identityK
  | k :: rest => comp k (history rest)

theorem history_preserves [DecidableEq I] (tag : I → T) (ks : List (Kernel I))
    (allowed : ∀ k ∈ ks, Preserves tag k) : Preserves tag (history ks) := by
  induction ks with
  | nil => exact sharp_preserves tag (fun _ => 1)
  | cons k rest ih =>
    exact comp_preserves tag k _ (allowed k List.mem_cons_self)
      (ih (fun l hl => allowed l (List.mem_cons_of_mem k hl)))

/-- Probability weights for one fine-grained branch, followed by sharp final records
with a common acceptance weight. The earlier reads and controls may vary by branch. -/
noncomputable def weight (action input : Kernel I) (finish : I → ℂ) : ℂ :=
  ∑ y, ∑ i, ∑ j, finish y * action y i * input i j * star (action y j)

theorem same_block_weight (tag : I → T) (action r s : Kernel I) (finish : I → ℂ)
    (ha : Preserves tag action) (agree : ∀ i j, tag i = tag j → r i j = s i j) :
    weight action r finish = weight action s finish := by
  unfold weight
  apply Finset.sum_congr rfl
  intro y _
  apply Finset.sum_congr rfl
  intro i _
  apply Finset.sum_congr rfl
  intro j _
  by_cases h : tag y = tag i
  · by_cases h' : tag y = tag j
    · rw [agree i j (h.symm.trans h')]
    · simp [ha y j h']
  · simp [ha y i h]

theorem branch_equality [DecidableEq I] (tag : I → T) (ks : List (Kernel I))
    (allowed : ∀ k ∈ ks, Preserves tag k) (r s : Kernel I) (finish : I → ℂ)
    (agree : ∀ i j, tag i = tag j → r i j = s i j) :
    weight (history ks) r finish = weight (history ks) s finish :=
  same_block_weight tag _ r s finish (history_preserves tag ks allowed) agree

/-- The sum of the two absolute quadrature responses at relative angle 120°
controls the length of the original coherence, with arbitrary orientation. -/
theorem two_axes (x y u : ℝ) (angle : 3 * y^2 = (x + 2*u)^2) :
    3 * (x^2 + y^2) ≤ 4 * (|x| + |u|)^2 := by
  have h := le_abs_self (x*u)
  rw [abs_mul] at h
  nlinarith [sq_abs x, sq_abs u, mul_nonneg (abs_nonneg x) (abs_nonneg u)]

end OPH.NeutralSourceAccess

#print axioms OPH.NeutralSourceAccess.comp_preserves
#print axioms OPH.NeutralSourceAccess.sharp_preserves
#print axioms OPH.NeutralSourceAccess.history_preserves
#print axioms OPH.NeutralSourceAccess.same_block_weight
#print axioms OPH.NeutralSourceAccess.branch_equality

#print axioms OPH.NeutralSourceAccess.two_axes
