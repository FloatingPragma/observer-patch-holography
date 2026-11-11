import Dynamics.ProtectedCharge

/-!
# Conjugation symmetry of finite record dynamics

An equivariant finite kernel preserves a conjugation-even state, and an odd
observable has zero expectation in that state. The native source must supply
the equivariance and preparation premises. No record label is identified with
physical baryon number by these linear identities.
-/

namespace OPH.Dynamics

open scoped BigOperators

variable {State : Type*} [Fintype State]

/-- A distribution invariant under the declared record permutation. -/
def ConjugationEven (c : Equiv.Perm State) (p : State → ℝ) : Prop :=
  ∀ x, p (c x) = p x

/-- A kernel intertwines conjugation on its input and output records. -/
def ConjugationEquivariant (c : Equiv.Perm State)
    (kernel : State → State → ℝ) : Prop :=
  ∀ x y, kernel (c x) (c y) = kernel x y

/-- Complete finite evolution preserves the even preparation premise. -/
theorem kernelPush_conjugationEven (c : Equiv.Perm State)
    (kernel : State → State → ℝ) (p : State → ℝ)
    (hK : ConjugationEquivariant c kernel) (hp : ConjugationEven c p) :
    ConjugationEven c (kernelPush kernel p) := by
  intro y
  change (∑ x, p x * kernel x (c y)) = ∑ x, p x * kernel x y
  calc
    (∑ x, p x * kernel x (c y)) =
        ∑ x, p (c x) * kernel (c x) (c y) := (Equiv.sum_comp c _).symm
    _ = ∑ x, p x * kernel x y := by
      apply Finset.sum_congr rfl
      intro x _
      rw [hp x, hK x y]

/-- Every finite number of complete steps preserves conjugation symmetry. -/
theorem iterate_kernelPush_conjugationEven (c : Equiv.Perm State)
    (kernel : State → State → ℝ) (p : State → ℝ)
    (hK : ConjugationEquivariant c kernel) (hp : ConjugationEven c p)
    (n : ℕ) :
    ConjugationEven c ((kernelPush kernel)^[n] p) := by
  induction n with
  | zero => simpa using hp
  | succ n ih =>
    rw [Function.iterate_succ_apply']
    exact kernelPush_conjugationEven c kernel _ hK ih

/-- Odd record observables have zero mean under an even preparation. -/
theorem odd_chargeExpectation_zero (c : Equiv.Perm State)
    (p observable : State → ℝ) (hp : ConjugationEven c p)
    (hodd : ∀ x, observable (c x) = -observable x) :
    chargeExpectation p observable = 0 := by
  have hsum : chargeExpectation p observable = -chargeExpectation p observable := by
    unfold chargeExpectation
    calc
      (∑ x, p x * observable x) = ∑ x, p (c x) * observable (c x) :=
        (Equiv.sum_comp c _).symm
      _ = ∑ x, -(p x * observable x) := by
        apply Finset.sum_congr rfl
        intro x _
        rw [hp x, hodd x, mul_neg]
      _ = -(∑ x, p x * observable x) := by rw [Finset.sum_neg_distrib]
  linarith

end OPH.Dynamics

#print axioms OPH.Dynamics.kernelPush_conjugationEven
#print axioms OPH.Dynamics.iterate_kernelPush_conjugationEven
#print axioms OPH.Dynamics.odd_chargeExpectation_zero
