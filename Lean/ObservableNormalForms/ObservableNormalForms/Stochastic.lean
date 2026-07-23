import Mathlib

/-!
# Finite Markov drift certificates

This module adds a finite-state stochastic certificate adjacent to the
manuscript's schedule discussion.  It proves the standard affine drift
iteration for the finite Markov operator and a one-time finite Markov tail
bound.  No conditional-expectation or martingale machinery is used.
-/

namespace ObservableNormalForms

universe u

/-- A finite Markov transition kernel represented by real probabilities. -/
structure FiniteMarkovKernel (S : Type u) [Fintype S] where
  probability : S → S → ℝ
  probability_nonneg : ∀ x y : S, 0 ≤ probability x y
  probability_sum_one : ∀ x : S, ∑ y : S, probability x y = 1

namespace FiniteMarkovKernel

variable {S : Type u} [Fintype S]

/-- The Markov operator on real-valued observables. -/
def apply (K : FiniteMarkovKernel S) (V : S → ℝ) (x : S) : ℝ :=
  ∑ y : S, K.probability x y * V y

theorem apply_mono (K : FiniteMarkovKernel S)
    {V W : S → ℝ} (hVW : ∀ x : S, V x ≤ W x) (x : S) :
    K.apply V x ≤ K.apply W x := by
  apply Finset.sum_le_sum
  intro y _
  exact mul_le_mul_of_nonneg_left (hVW y) (K.probability_nonneg x y)

theorem apply_affine (K : FiniteMarkovKernel S)
    (a b : ℝ) (V : S → ℝ) (x : S) :
    K.apply (fun y => a * V y + b) x = a * K.apply V x + b := by
  simp only [apply]
  calc
    (∑ y : S, K.probability x y * (a * V y + b))
        = ∑ y : S,
            (a * (K.probability x y * V y) + K.probability x y * b) := by
          apply Finset.sum_congr rfl
          intro y _
          ring
    _ = a * (∑ y : S, K.probability x y * V y) +
          (∑ y : S, K.probability x y) * b := by
          rw [Finset.sum_add_distrib, Finset.mul_sum, Finset.sum_mul]
    _ = a * (∑ y : S, K.probability x y * V y) + b := by
          rw [K.probability_sum_one]
          ring

/-- `iterateExpectation n V x` is the expected value of `V` after `n`
kernel steps when the initial state is `x`. -/
def iterateExpectation (K : FiniteMarkovKernel S) :
    ℕ → (S → ℝ) → S → ℝ
  | 0, V => V
  | n + 1, V => K.apply (K.iterateExpectation n V)

/-- Geometric-sum form of the finite Markov affine-drift iteration. -/
theorem drift_iteration_geomSum
    (K : FiniteMarkovKernel S)
    {V : S → ℝ} {κ ε : ℝ}
    (hκ : 0 ≤ κ)
    (hdrift : ∀ x : S, K.apply V x ≤ κ * V x + ε)
    (n : ℕ) (x : S) :
    K.iterateExpectation n V x ≤
      κ ^ n * V x + ε * ∑ i ∈ Finset.range n, κ ^ i := by
  induction n generalizing x with
  | zero => simp [iterateExpectation]
  | succ n ih =>
      have hmono :
          K.apply (K.iterateExpectation n V) x ≤
            K.apply (fun y =>
              κ ^ n * V y + ε * ∑ i ∈ Finset.range n, κ ^ i) x :=
        K.apply_mono (fun y => ih y) x
      calc
        K.iterateExpectation (n + 1) V x
            = K.apply (K.iterateExpectation n V) x := rfl
        _ ≤ K.apply (fun y =>
              κ ^ n * V y + ε * ∑ i ∈ Finset.range n, κ ^ i) x := hmono
        _ = κ ^ n * K.apply V x +
              ε * ∑ i ∈ Finset.range n, κ ^ i := by
              rw [K.apply_affine]
        _ ≤ κ ^ n * (κ * V x + ε) +
              ε * ∑ i ∈ Finset.range n, κ ^ i :=
              add_le_add
                (mul_le_mul_of_nonneg_left (hdrift x) (pow_nonneg hκ n)) le_rfl
        _ = κ ^ (n + 1) * V x +
              ε * ∑ i ∈ Finset.range (n + 1), κ ^ i := by
              rw [Finset.sum_range_succ]
              ring

/-- Requested closed form: if `0 ≤ κ < 1`, then
`E[V(X_n)] ≤ κ^n V(x) + ε(1-κ^n)/(1-κ)`. -/
theorem finite_markov_drift_iteration
    (K : FiniteMarkovKernel S)
    {V : S → ℝ} {κ ε : ℝ}
    (hκ0 : 0 ≤ κ)
    (hκ1 : κ < 1)
    (hdrift : ∀ x : S, K.apply V x ≤ κ * V x + ε)
    (n : ℕ) (x : S) :
    K.iterateExpectation n V x ≤
      κ ^ n * V x + ε * ((1 - κ ^ n) / (1 - κ)) := by
  calc
    K.iterateExpectation n V x ≤
        κ ^ n * V x + ε * ∑ i ∈ Finset.range n, κ ^ i :=
      K.drift_iteration_geomSum hκ0 hdrift n x
    _ = κ ^ n * V x + ε * ((1 - κ ^ n) / (1 - κ)) := by
      rw [geom_sum_eq hκ1.ne]
      have hfrac :
          (κ ^ n - 1) / (κ - 1) = (1 - κ ^ n) / (1 - κ) := by
        rw [show 1 - κ ^ n = -(κ ^ n - 1) by ring,
          show 1 - κ = -(κ - 1) by ring, neg_div_neg_eq]
      rw [hfrac]

end FiniteMarkovKernel

section FiniteTail

variable {S : Type u} [Fintype S]

def finiteExpectation (μ V : S → ℝ) : ℝ :=
  ∑ x : S, μ x * V x

def finiteEventMass (μ : S → ℝ) (p : S → Prop) [DecidablePred p] : ℝ :=
  ∑ x ∈ Finset.univ.filter p, μ x

/-- One-time Markov inequality for a finite distribution.  Normalization of
`μ` is not required for the algebraic inequality itself; nonnegativity is. -/
theorem finite_markov_tail_bound
    (μ V : S → ℝ)
    (hμ : ∀ x : S, 0 ≤ μ x)
    (hV : ∀ x : S, 0 ≤ V x)
    {a : ℝ} (ha : 0 < a) :
    finiteEventMass μ (fun x => a ≤ V x) ≤
      finiteExpectation μ V / a := by
  classical
  apply (le_div_iff₀ ha).2
  rw [mul_comm]
  calc
    a * finiteEventMass μ (fun x => a ≤ V x)
        = ∑ x : S, a * (if a ≤ V x then μ x else 0) := by
          rw [finiteEventMass, Finset.mul_sum, Finset.sum_filter]
          simp
    _ ≤ ∑ x : S, μ x * V x := by
      apply Finset.sum_le_sum
      intro x _
      by_cases hx : a ≤ V x
      · simp only [hx, if_true]
        nlinarith [hμ x]
      · simp only [hx, if_false, mul_zero]
        exact mul_nonneg (hμ x) (hV x)
    _ = finiteExpectation μ V := rfl

end FiniteTail

end ObservableNormalForms
