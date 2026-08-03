import Mathlib
import ObservableNormalForms.Exact

/-!
# Finite Markov drift certificates

This module adds a finite-state stochastic certificate adjacent to the
manuscript's schedule discussion.  It proves the standard affine drift
iteration for the finite Markov operator and a one-time finite Markov tail
bound.  No conditional-expectation or martingale machinery is used.

Beyond the drift iteration, the module completes `prop:markov-receipt`:

* an explicit finite path law (`pathWeight` on tuples `Fin n → S`) with the
  `n`-step endpoint distribution proved to be its endpoint marginal;
* the identity between expectation under the endpoint distribution and the
  operator iterate, which turns the drift iterate into a statement about an
  actual distribution;
* the positive-support layer: supported steps, supported paths, and
  observation preservation along every supported path;
* the packaged one-time endpoint receipt `finite_markov_endpoint_receipt`;
* the stationary and occupation wrappers of `cor:markov-occupation`;
* a persistent-noise negative control and support-checker rejection tests.

Everything is finite: sums over `Fintype`, no measure theory.
-/

namespace ObservableNormalForms

universe u v

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

/-! ## Finite event-mass toolbox

Small algebraic lemmas about `finiteEventMass` needed to phrase tail and
settling events for one distribution. -/

section EventMass

variable {S : Type u} [Fintype S]

theorem finiteEventMass_congr (μ : S → ℝ) {p q : S → Prop}
    [DecidablePred p] [DecidablePred q] (h : ∀ x, p x ↔ q x) :
    finiteEventMass μ p = finiteEventMass μ q := by
  unfold finiteEventMass
  congr 1
  exact Finset.filter_congr fun x _ => h x

theorem finiteEventMass_mono (μ : S → ℝ) (hμ : ∀ x : S, 0 ≤ μ x)
    {p q : S → Prop} [DecidablePred p] [DecidablePred q]
    (h : ∀ x, p x → q x) :
    finiteEventMass μ p ≤ finiteEventMass μ q := by
  unfold finiteEventMass
  apply Finset.sum_le_sum_of_subset_of_nonneg
  · intro x hx
    rw [Finset.mem_filter] at hx ⊢
    exact ⟨hx.1, h x hx.2⟩
  · intro x _ _
    exact hμ x

theorem finiteEventMass_add_compl (μ : S → ℝ) (p : S → Prop) [DecidablePred p] :
    finiteEventMass μ p + finiteEventMass μ (fun x => ¬ p x) = ∑ x : S, μ x := by
  unfold finiteEventMass
  exact Finset.sum_filter_add_sum_filter_not _ _ _

end EventMass

/-! ## Final state of a listed path -/

section PathTarget

variable {α : Type u}

/-- Final state of a listed path: the start point when the list is empty,
otherwise the last listed state. -/
def pathTarget (x : α) : List α → α
  | [] => x
  | y :: l => pathTarget y l

theorem pathTarget_mem_cons : ∀ (l : List α) (x : α), pathTarget x l ∈ x :: l := by
  intro l
  induction l with
  | nil => intro x; simp [pathTarget]
  | cons y t ih =>
      intro x
      have h := ih y
      simp only [pathTarget]
      exact List.mem_cons_of_mem x h

end PathTarget

/-! ## The `n`-step endpoint distribution and the finite path law -/

namespace FiniteMarkovKernel

variable {S : Type u} [Fintype S] [DecidableEq S]

/-- The law of the chain after `n` kernel steps from `x`, by first-step
decomposition.  `endpointDistribution 0 x` is the Dirac mass at `x`. -/
def endpointDistribution (K : FiniteMarkovKernel S) : ℕ → S → S → ℝ
  | 0, x, y => if y = x then 1 else 0
  | n + 1, x, y => ∑ z : S, K.probability x z * K.endpointDistribution n z y

theorem endpointDistribution_nonneg (K : FiniteMarkovKernel S) :
    ∀ (n : ℕ) (x y : S), 0 ≤ K.endpointDistribution n x y := by
  intro n
  induction n with
  | zero =>
      intro x y
      simp only [endpointDistribution]
      split <;> norm_num
  | succ n ih =>
      intro x y
      simp only [endpointDistribution]
      exact Finset.sum_nonneg fun z _ =>
        mul_nonneg (K.probability_nonneg x z) (ih z y)

theorem endpointDistribution_sum_one (K : FiniteMarkovKernel S) :
    ∀ (n : ℕ) (x : S), ∑ y : S, K.endpointDistribution n x y = 1 := by
  intro n
  induction n with
  | zero =>
      intro x
      simp [endpointDistribution]
  | succ n ih =>
      intro x
      simp only [endpointDistribution]
      rw [Finset.sum_comm]
      calc ∑ z : S, ∑ y : S, K.probability x z * K.endpointDistribution n z y
          = ∑ z : S, K.probability x z * ∑ y : S, K.endpointDistribution n z y := by
            simp [Finset.mul_sum]
        _ = ∑ z : S, K.probability x z := by simp [ih]
        _ = 1 := K.probability_sum_one x

/-- **Kernel-law wrapper.**  Expectation of `V` under the `n`-step endpoint
distribution equals the operator iterate `iterateExpectation`.  This is the
bridge that turns the algebraic drift iterate into a statement about an
actual finite distribution. -/
theorem finiteExpectation_endpointDistribution
    (K : FiniteMarkovKernel S) (V : S → ℝ) :
    ∀ (n : ℕ) (x : S),
      finiteExpectation (K.endpointDistribution n x) V =
        K.iterateExpectation n V x := by
  intro n
  induction n with
  | zero =>
      intro x
      simp [finiteExpectation, endpointDistribution, iterateExpectation, ite_mul]
  | succ n ih =>
      intro x
      calc finiteExpectation (K.endpointDistribution (n + 1) x) V
          = ∑ y : S, (∑ z : S, K.probability x z * K.endpointDistribution n z y)
              * V y := by
            simp only [finiteExpectation, endpointDistribution]
        _ = ∑ y : S, ∑ z : S,
              K.probability x z * K.endpointDistribution n z y * V y := by
            simp [Finset.sum_mul]
        _ = ∑ z : S, ∑ y : S,
              K.probability x z * K.endpointDistribution n z y * V y :=
            Finset.sum_comm
        _ = ∑ z : S, K.probability x z *
              ∑ y : S, K.endpointDistribution n z y * V y := by
            simp [Finset.mul_sum, mul_assoc]
        _ = ∑ z : S, K.probability x z * K.iterateExpectation n V z := by
            refine Finset.sum_congr rfl fun z _ => ?_
            rw [← ih z, finiteExpectation]
        _ = K.iterateExpectation (n + 1) V x := rfl

/-- Product weight of the `n`-step path `p` started at `x`: the probability
that the chain visits exactly the listed states in order. -/
def pathWeight (K : FiniteMarkovKernel S) : (n : ℕ) → S → (Fin n → S) → ℝ
  | 0, _, _ => 1
  | n + 1, x, p => K.probability x (p 0) * K.pathWeight n (p 0) (Fin.tail p)

/-- Endpoint of an `n`-step path started at `x`. -/
def pathEndpoint : (n : ℕ) → S → (Fin n → S) → S
  | 0, x, _ => x
  | n + 1, _, p => p (Fin.last n)

omit [DecidableEq S] in
theorem pathWeight_cons (K : FiniteMarkovKernel S) (n : ℕ) (x z : S)
    (q : Fin n → S) :
    K.pathWeight (n + 1) x (Fin.cons z q) =
      K.probability x z * K.pathWeight n z q := by
  simp [pathWeight, Fin.tail_cons]

omit [Fintype S] [DecidableEq S] in
theorem pathEndpoint_cons (n : ℕ) (x z : S) (q : Fin n → S) :
    pathEndpoint (n + 1) x (Fin.cons z q) = pathEndpoint n z q := by
  cases n with
  | zero =>
      simp [pathEndpoint, show (Fin.last 0 : Fin 1) = 0 from rfl]
  | succ m =>
      simp [pathEndpoint, ← Fin.succ_last, Fin.cons_succ]

/-- Splitting an `(n+1)`-tuple path into its first state and its tail. -/
private def finPiConsEquiv (S : Type u) (n : ℕ) :
    S × (Fin n → S) ≃ (Fin (n + 1) → S) where
  toFun zq := Fin.cons zq.1 zq.2
  invFun p := (p 0, Fin.tail p)
  left_inv := fun zq => by simp
  right_inv := fun p => Fin.cons_self_tail p

omit [DecidableEq S] in
private theorem sum_pi_fin_succ {n : ℕ} (f : (Fin (n + 1) → S) → ℝ) :
    ∑ p : Fin (n + 1) → S, f p =
      ∑ z : S, ∑ q : Fin n → S, f (Fin.cons z q) := by
  calc ∑ p : Fin (n + 1) → S, f p
      = ∑ zq : S × (Fin n → S), f (Fin.cons zq.1 zq.2) :=
        ((finPiConsEquiv S n).sum_comp f).symm
    _ = ∑ z : S, ∑ q : Fin n → S, f (Fin.cons z q) := by
        rw [Fintype.sum_prod_type]

omit [DecidableEq S] in
theorem pathWeight_nonneg (K : FiniteMarkovKernel S) :
    ∀ (n : ℕ) (x : S) (p : Fin n → S), 0 ≤ K.pathWeight n x p := by
  intro n
  induction n with
  | zero => intro x p; simp [pathWeight]
  | succ n ih =>
      intro x p
      simp only [pathWeight]
      exact mul_nonneg (K.probability_nonneg x (p 0)) (ih (p 0) (Fin.tail p))

omit [DecidableEq S] in
/-- The path law is a probability law on `n`-step paths. -/
theorem pathWeight_sum_one (K : FiniteMarkovKernel S) :
    ∀ (n : ℕ) (x : S), ∑ p : Fin n → S, K.pathWeight n x p = 1 := by
  intro n
  induction n with
  | zero =>
      intro x
      rw [Fintype.sum_unique]
      rfl
  | succ n ih =>
      intro x
      rw [sum_pi_fin_succ]
      calc ∑ z : S, ∑ q : Fin n → S, K.pathWeight (n + 1) x (Fin.cons z q)
          = ∑ z : S, ∑ q : Fin n → S, K.probability x z * K.pathWeight n z q := by
            refine Finset.sum_congr rfl fun z _ => Finset.sum_congr rfl fun q _ => ?_
            rw [pathWeight_cons]
        _ = ∑ z : S, K.probability x z * ∑ q : Fin n → S, K.pathWeight n z q := by
            simp [Finset.mul_sum]
        _ = ∑ z : S, K.probability x z := by simp [ih]
        _ = 1 := K.probability_sum_one x

/-- **The finite path law.**  The `n`-step endpoint distribution is the
endpoint marginal of the product path weights: nothing about `X_n` is assumed
beyond the kernel entries along explicit finite paths. -/
theorem endpointDistribution_eq_pathWeight_sum (K : FiniteMarkovKernel S) :
    ∀ (n : ℕ) (x y : S),
      K.endpointDistribution n x y =
        ∑ p : Fin n → S,
          if pathEndpoint n x p = y then K.pathWeight n x p else 0 := by
  intro n
  induction n with
  | zero =>
      intro x y
      simp only [endpointDistribution, pathEndpoint, pathWeight]
      rw [Fintype.sum_unique]
      by_cases h : y = x
      · subst h; simp
      · rw [if_neg h, if_neg fun hxy => h hxy.symm]
  | succ n ih =>
      intro x y
      calc K.endpointDistribution (n + 1) x y
          = ∑ z : S, K.probability x z * K.endpointDistribution n z y := by
            simp only [endpointDistribution]
        _ = ∑ z : S, K.probability x z *
              ∑ q : Fin n → S,
                if pathEndpoint n z q = y then K.pathWeight n z q else 0 := by
            simp only [ih]
        _ = ∑ z : S, ∑ q : Fin n → S,
              if pathEndpoint n z q = y
                then K.probability x z * K.pathWeight n z q else 0 := by
            refine Finset.sum_congr rfl fun z _ => ?_
            rw [Finset.mul_sum]
            refine Finset.sum_congr rfl fun q _ => ?_
            rw [mul_ite, mul_zero]
        _ = ∑ p : Fin (n + 1) → S,
              if pathEndpoint (n + 1) x p = y
                then K.pathWeight (n + 1) x p else 0 := by
            rw [sum_pi_fin_succ]
            refine Finset.sum_congr rfl fun z _ =>
              Finset.sum_congr rfl fun q _ => ?_
            rw [pathEndpoint_cons, pathWeight_cons]

/-! ## Positive support: supported steps, supported paths, observation -/

/-- Every step of the listed path has positive kernel probability. -/
def SupportedPath (K : FiniteMarkovKernel S) : S → List S → Prop
  | _, [] => True
  | x, y :: l => 0 < K.probability x y ∧ K.SupportedPath y l

/-- Positive endpoint mass is witnessed by an explicit supported path: the
event `X_n = y` cannot receive mass except through states the kernel can
actually visit. -/
theorem exists_supportedPath_of_endpointDistribution_pos
    (K : FiniteMarkovKernel S) :
    ∀ {n : ℕ} {x y : S}, 0 < K.endpointDistribution n x y →
      ∃ l : List S, l.length = n ∧ K.SupportedPath x l ∧ pathTarget x l = y := by
  intro n
  induction n with
  | zero =>
      intro x y h
      simp only [endpointDistribution] at h
      by_cases hyx : y = x
      · exact ⟨[], rfl, trivial, by simp [pathTarget, hyx]⟩
      · rw [if_neg hyx] at h
        exact absurd h (lt_irrefl 0)
  | succ n ih =>
      intro x y h
      simp only [endpointDistribution] at h
      have h' : ∑ z : S, (0 : ℝ) <
          ∑ z : S, K.probability x z * K.endpointDistribution n z y := by
        simpa using h
      obtain ⟨z, -, hz⟩ := Finset.exists_lt_of_sum_lt h'
      have hzpos : 0 < K.probability x z ∧ 0 < K.endpointDistribution n z y := by
        rcases mul_pos_iff.mp hz with h1 | h1
        · exact h1
        · exact absurd h1.1 (not_lt.mpr (K.probability_nonneg x z))
      obtain ⟨l, hlen, hpath, hlast⟩ := ih hzpos.2
      exact ⟨z :: l, by simp [hlen], ⟨hzpos.1, hpath⟩, hlast⟩

omit [DecidableEq S] in
/-- **Observation preservation along supported paths.**  If every positive-
probability step of the kernel is a permitted rewrite or a stutter, and the
rewrite relation preserves the observation, then every state visited by a
supported path carries the initial observation. -/
theorem observation_eq_of_mem_supportedPath
    {𝓑 : Type v} {r : S → S → Prop} {B : S → 𝓑}
    (K : FiniteMarkovKernel S)
    (hsup : ∀ ⦃a b : S⦄, 0 < K.probability a b → r a b ∨ b = a)
    (hobs : ObservationPreserving r B) :
    ∀ {l : List S} {x : S}, K.SupportedPath x l → ∀ y ∈ l, B y = B x := by
  intro l
  induction l with
  | nil =>
      intro x _ y hy
      simp at hy
  | cons z t ih =>
      intro x hpath y hy
      obtain ⟨hstep, htail⟩ := hpath
      have hBz : B z = B x := by
        rcases hsup hstep with hr | heq
        · exact (hobs hr).symm
        · rw [heq]
      rcases List.mem_cons.mp hy with h | h
      · rw [h, hBz]
      · rw [ih htail y h, hBz]

/-- Endpoint form: any state with positive `n`-step mass carries the initial
observation. -/
theorem observation_eq_of_endpointDistribution_pos
    {𝓑 : Type v} {r : S → S → Prop} {B : S → 𝓑}
    (K : FiniteMarkovKernel S)
    (hsup : ∀ ⦃a b : S⦄, 0 < K.probability a b → r a b ∨ b = a)
    (hobs : ObservationPreserving r B)
    {n : ℕ} {x y : S} (h : 0 < K.endpointDistribution n x y) :
    B y = B x := by
  obtain ⟨l, -, hpath, htarget⟩ :=
    K.exists_supportedPath_of_endpointDistribution_pos h
  have hy : y ∈ x :: l := by
    have hmem := pathTarget_mem_cons l x
    rw [htarget] at hmem
    exact hmem
  rcases List.mem_cons.mp hy with h' | h'
  · rw [h']
  · exact K.observation_eq_of_mem_supportedPath hsup hobs hpath y h'

/-! ## The one-time endpoint tail bound and the packaged receipt -/

/-- Drift-to-tail composition on the actual endpoint distribution:
`P_x[V(X_n) > δ] ≤ (κ^n V(x) + ξ(1-κ^n)/(1-κ)) / δ`. -/
theorem endpointDistribution_tail_bound
    (K : FiniteMarkovKernel S) {V : S → ℝ} {κ ξ δ : ℝ}
    (hV : ∀ x : S, 0 ≤ V x) (hκ0 : 0 ≤ κ) (hκ1 : κ < 1)
    (hdrift : ∀ x : S, K.apply V x ≤ κ * V x + ξ)
    (hδ : 0 < δ) (n : ℕ) (x : S) :
    finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y) ≤
      (κ ^ n * V x + ξ * ((1 - κ ^ n) / (1 - κ))) / δ := by
  have h1 : finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y) ≤
      finiteEventMass (K.endpointDistribution n x) (fun y => δ ≤ V y) :=
    finiteEventMass_mono _ (K.endpointDistribution_nonneg n x)
      fun y hy => le_of_lt hy
  have h2 := finite_markov_tail_bound (K.endpointDistribution n x) V
    (K.endpointDistribution_nonneg n x) hV hδ
  have h3 : finiteExpectation (K.endpointDistribution n x) V ≤
      κ ^ n * V x + ξ * ((1 - κ ^ n) / (1 - κ)) := by
    rw [K.finiteExpectation_endpointDistribution]
    exact K.finite_markov_drift_iteration hκ0 hκ1 hdrift n x
  calc finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y)
      ≤ finiteEventMass (K.endpointDistribution n x) (fun y => δ ≤ V y) := h1
    _ ≤ finiteExpectation (K.endpointDistribution n x) V / δ := h2
    _ ≤ (κ ^ n * V x + ξ * ((1 - κ ^ n) / (1 - κ))) / δ := by gcongr

/-- **Finite Markov drift receipt** (`prop:markov-receipt`), packaged.

Under the affine drift condition, positive-support inclusion in permitted
rewrites or stutters, and observation preservation of the rewrite relation,
the `n`-step endpoint distribution satisfies, for `p := (κ^n V(x) +
ξ(1-κ^n)/(1-κ))/δ`:

1. tail bound `P_x[V(X_n) > δ] ≤ p`;
2. settling mass `P_x[V(X_n) ≤ δ] ≥ 1 - p`, the `(δ,0)`-settling receipt
   consumed by `cor:probabilistic`;
3. observation preservation on the endpoint support:
   `P_x[X_n = y] > 0 → B(y) = B(x)`.

**Boundary.** This is a ONE-TIME ENDPOINT bound: it constrains the law of
`X_n` for each fixed `n` separately, and by (3) the observation along
supported paths.  It does NOT bound the probability that the whole
trajectory stays in the `δ`-tube over a time window; see the persistent-
noise control `PersistentNoiseControl` below for a kernel where every
one-time bound holds while window confinement fails badly. -/
theorem finite_markov_endpoint_receipt
    {𝓑 : Type v}
    (K : FiniteMarkovKernel S) (r : S → S → Prop) (B : S → 𝓑)
    {V : S → ℝ} {κ ξ δ : ℝ}
    (hV : ∀ x : S, 0 ≤ V x) (hκ0 : 0 ≤ κ) (hκ1 : κ < 1)
    (hdrift : ∀ x : S, K.apply V x ≤ κ * V x + ξ)
    (hsup : ∀ ⦃a b : S⦄, 0 < K.probability a b → r a b ∨ b = a)
    (hobs : ObservationPreserving r B)
    (hδ : 0 < δ) (n : ℕ) (x : S) :
    finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y) ≤
        (κ ^ n * V x + ξ * ((1 - κ ^ n) / (1 - κ))) / δ ∧
    1 - (κ ^ n * V x + ξ * ((1 - κ ^ n) / (1 - κ))) / δ ≤
        finiteEventMass (K.endpointDistribution n x) (fun y => V y ≤ δ) ∧
    ∀ y : S, 0 < K.endpointDistribution n x y → B y = B x := by
  have htail := K.endpointDistribution_tail_bound hV hκ0 hκ1 hdrift hδ n x
  refine ⟨htail, ?_, fun y hy =>
    K.observation_eq_of_endpointDistribution_pos hsup hobs hy⟩
  have hcompl : finiteEventMass (K.endpointDistribution n x) (fun y => V y ≤ δ) +
      finiteEventMass (K.endpointDistribution n x) (fun y => ¬ V y ≤ δ) = 1 := by
    have h := finiteEventMass_add_compl
      (K.endpointDistribution n x) (fun y => V y ≤ δ)
    rw [K.endpointDistribution_sum_one n x] at h
    exact h
  have hcongr : finiteEventMass (K.endpointDistribution n x)
      (fun y => ¬ V y ≤ δ) =
      finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y) :=
    finiteEventMass_congr _ fun y => not_le
  linarith [htail, hcompl, hcongr]

/-! ## Stationary and occupation wrappers (`cor:markov-occupation`) -/

/-- A stationary distribution for a finite kernel. -/
structure IsStationary (K : FiniteMarkovKernel S) (π : S → ℝ) : Prop where
  nonneg : ∀ x : S, 0 ≤ π x
  sum_one : ∑ x : S, π x = 1
  fixed : ∀ y : S, ∑ x : S, π x * K.probability x y = π y

omit [DecidableEq S] in
theorem finiteExpectation_apply_stationary
    {K : FiniteMarkovKernel S} {π : S → ℝ} (hπ : K.IsStationary π)
    (V : S → ℝ) :
    finiteExpectation π (K.apply V) = finiteExpectation π V := by
  unfold finiteExpectation FiniteMarkovKernel.apply
  calc ∑ x : S, π x * ∑ y : S, K.probability x y * V y
      = ∑ x : S, ∑ y : S, π x * K.probability x y * V y := by
        simp [Finset.mul_sum, mul_assoc]
    _ = ∑ y : S, ∑ x : S, π x * K.probability x y * V y := Finset.sum_comm
    _ = ∑ y : S, (∑ x : S, π x * K.probability x y) * V y := by
        simp [Finset.sum_mul]
    _ = ∑ y : S, π y * V y := by
        refine Finset.sum_congr rfl fun y _ => ?_
        rw [hπ.fixed y]

omit [DecidableEq S] in
/-- Stationary expectation bound `E_π V ≤ ξ/(1-κ)`. -/
theorem stationary_expectation_bound
    {K : FiniteMarkovKernel S} {π : S → ℝ} (hπ : K.IsStationary π)
    {V : S → ℝ} {κ ξ : ℝ} (hκ1 : κ < 1)
    (hdrift : ∀ x : S, K.apply V x ≤ κ * V x + ξ) :
    finiteExpectation π V ≤ ξ / (1 - κ) := by
  have h1κ : 0 < 1 - κ := by linarith
  have hstep : finiteExpectation π V ≤ κ * finiteExpectation π V + ξ := by
    calc finiteExpectation π V
        = finiteExpectation π (K.apply V) :=
          (finiteExpectation_apply_stationary hπ V).symm
      _ ≤ ∑ x : S, π x * (κ * V x + ξ) := by
          unfold finiteExpectation
          exact Finset.sum_le_sum fun x _ =>
            mul_le_mul_of_nonneg_left (hdrift x) (hπ.nonneg x)
      _ = κ * (∑ x : S, π x * V x) + ξ * ∑ x : S, π x := by
          rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
          exact Finset.sum_congr rfl fun x _ => by ring
      _ = κ * finiteExpectation π V + ξ := by
          simp only [finiteExpectation, hπ.sum_one, mul_one]
  rw [le_div_iff₀ h1κ]
  nlinarith [hstep]

omit [DecidableEq S] in
/-- Stationary tail bound `π{V > δ} ≤ ξ/((1-κ)δ)`. -/
theorem stationary_tail_bound
    {K : FiniteMarkovKernel S} {π : S → ℝ} (hπ : K.IsStationary π)
    {V : S → ℝ} {κ ξ δ : ℝ}
    (hV : ∀ x : S, 0 ≤ V x) (hκ1 : κ < 1)
    (hdrift : ∀ x : S, K.apply V x ≤ κ * V x + ξ)
    (hδ : 0 < δ) :
    finiteEventMass π (fun y => δ < V y) ≤ ξ / ((1 - κ) * δ) := by
  have h1κ : 0 < 1 - κ := by linarith
  have h1 : finiteEventMass π (fun y => δ < V y) ≤
      finiteEventMass π (fun y => δ ≤ V y) :=
    finiteEventMass_mono _ hπ.nonneg fun y hy => le_of_lt hy
  have h2 := finite_markov_tail_bound π V hπ.nonneg hV hδ
  have h3 := stationary_expectation_bound hπ hκ1 hdrift
  calc finiteEventMass π (fun y => δ < V y)
      ≤ finiteEventMass π (fun y => δ ≤ V y) := h1
    _ ≤ finiteExpectation π V / δ := h2
    _ ≤ (ξ / (1 - κ)) / δ := by gcongr
    _ = ξ / ((1 - κ) * δ) := by rw [div_div]

omit [Fintype S] [DecidableEq S] in
private theorem geom_sum_le_one_div {κ : ℝ} (hκ0 : 0 ≤ κ) (hκ1 : κ < 1)
    (N : ℕ) : ∑ n ∈ Finset.range N, κ ^ n ≤ 1 / (1 - κ) := by
  have h1κ : 0 < 1 - κ := by linarith
  rw [geom_sum_eq hκ1.ne]
  rw [show (κ ^ N - 1) / (κ - 1) = (1 - κ ^ N) / (1 - κ) by
    rw [show 1 - κ ^ N = -(κ ^ N - 1) by ring,
      show 1 - κ = -(κ - 1) by ring, neg_div_neg_eq]]
  have hpow : (0 : ℝ) ≤ κ ^ N := pow_nonneg hκ0 N
  gcongr
  linarith

/-- **Occupation receipt** (summed form): expected number of the first `N`
endpoint laws giving mass to the excursion event is bounded by
`V(x)/(δ(1-κ)) + N ξ/(δ(1-κ))`.

**Boundary.**  The left side is a SUM OF ONE-TIME endpoint masses (equally:
`N` times the excursion probability at an independent uniformly sampled
time).  It is NOT the probability of an excursion-free window; when `ξ > 0`
the bound grows linearly in `N` and never certifies infinite-horizon
pathwise confinement. -/
theorem occupation_tail_bound
    (K : FiniteMarkovKernel S) {V : S → ℝ} {κ ξ δ : ℝ}
    (hV : ∀ x : S, 0 ≤ V x) (hκ0 : 0 ≤ κ) (hκ1 : κ < 1) (hξ : 0 ≤ ξ)
    (hdrift : ∀ x : S, K.apply V x ≤ κ * V x + ξ)
    (hδ : 0 < δ) (N : ℕ) (x : S) :
    ∑ n ∈ Finset.range N,
        finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y) ≤
      V x / (δ * (1 - κ)) + N * (ξ / (δ * (1 - κ))) := by
  have h1κ : 0 < 1 - κ := by linarith
  have hterm : ∀ n : ℕ,
      finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y) ≤
        κ ^ n * (V x / δ) + ξ / ((1 - κ) * δ) := by
    intro n
    have hb := K.endpointDistribution_tail_bound hV hκ0 hκ1 hdrift hδ n x
    have hfrac : ξ * ((1 - κ ^ n) / (1 - κ)) ≤ ξ / (1 - κ) := by
      have hpow : (0 : ℝ) ≤ κ ^ n := pow_nonneg hκ0 n
      rw [← mul_div_assoc]
      gcongr
      nlinarith
    calc finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y)
        ≤ (κ ^ n * V x + ξ * ((1 - κ ^ n) / (1 - κ))) / δ := hb
      _ ≤ (κ ^ n * V x + ξ / (1 - κ)) / δ := by gcongr
      _ = κ ^ n * (V x / δ) + ξ / ((1 - κ) * δ) := by
          rw [add_div, mul_div_assoc, div_div]
  calc ∑ n ∈ Finset.range N,
        finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y)
      ≤ ∑ n ∈ Finset.range N, (κ ^ n * (V x / δ) + ξ / ((1 - κ) * δ)) :=
        Finset.sum_le_sum fun n _ => hterm n
    _ = (∑ n ∈ Finset.range N, κ ^ n) * (V x / δ) +
          N * (ξ / ((1 - κ) * δ)) := by
        rw [Finset.sum_add_distrib, ← Finset.sum_mul, Finset.sum_const,
          Finset.card_range, nsmul_eq_mul]
    _ ≤ (1 / (1 - κ)) * (V x / δ) + N * (ξ / ((1 - κ) * δ)) := by
        have hVδ : 0 ≤ V x / δ := div_nonneg (hV x) (le_of_lt hδ)
        have hgeom := geom_sum_le_one_div hκ0 hκ1 N
        gcongr
    _ = V x / (δ * (1 - κ)) + N * (ξ / (δ * (1 - κ))) := by
        rw [one_div, inv_mul_eq_div, div_div, mul_comm (1 - κ) δ]

/-- Occupation receipt, averaged form: the excursion frequency over the
first `N` times (equally, the excursion probability at an independent
uniformly sampled time in `{0, …, N-1}`) is at most
`V(x)/(Nδ(1-κ)) + ξ/(δ(1-κ))`. -/
theorem occupation_tail_bound_average
    (K : FiniteMarkovKernel S) {V : S → ℝ} {κ ξ δ : ℝ}
    (hV : ∀ x : S, 0 ≤ V x) (hκ0 : 0 ≤ κ) (hκ1 : κ < 1) (hξ : 0 ≤ ξ)
    (hdrift : ∀ x : S, K.apply V x ≤ κ * V x + ξ)
    (hδ : 0 < δ) {N : ℕ} (hN : 0 < N) (x : S) :
    (∑ n ∈ Finset.range N,
        finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y)) / N ≤
      V x / (N * δ * (1 - κ)) + ξ / (δ * (1 - κ)) := by
  have h1κ : 0 < 1 - κ := by linarith
  have hNR : (0 : ℝ) < N := by exact_mod_cast hN
  have hδ0 : δ ≠ 0 := ne_of_gt hδ
  have hκ0' : (1 : ℝ) - κ ≠ 0 := ne_of_gt h1κ
  have hN0 : (N : ℝ) ≠ 0 := ne_of_gt hNR
  have h := K.occupation_tail_bound hV hκ0 hκ1 hξ hdrift hδ N x
  calc (∑ n ∈ Finset.range N,
        finiteEventMass (K.endpointDistribution n x) (fun y => δ < V y)) / N
      ≤ (V x / (δ * (1 - κ)) + N * (ξ / (δ * (1 - κ)))) / N := by gcongr
    _ = V x / (N * δ * (1 - κ)) + ξ / (δ * (1 - κ)) := by
        field_simp

/-! ## Window confinement (the pathwise event) -/

/-- Probability that the first `n` sampled states all lie in `A`: the
pathwise (window) confinement event, as a finite recursion. -/
def confinedMass (K : FiniteMarkovKernel S) (A : S → Prop) [DecidablePred A] :
    ℕ → S → ℝ
  | 0, _ => 1
  | n + 1, x => ∑ y : S, if A y then K.probability x y * K.confinedMass A n y else 0

omit [DecidableEq S] in
/-- The window-confinement mass is the path-law mass of the all-coordinates
event: `confinedMass` really is the probability of the pathwise event, under
the same finite path law whose endpoint marginal is `endpointDistribution`. -/
theorem confinedMass_eq_pathWeight_sum
    (K : FiniteMarkovKernel S) (A : S → Prop) [DecidablePred A] :
    ∀ (n : ℕ) (x : S),
      K.confinedMass A n x =
        ∑ p : Fin n → S,
          if (∀ i, A (p i)) then K.pathWeight n x p else 0 := by
  intro n
  induction n with
  | zero =>
      intro x
      simp only [confinedMass]
      rw [Fintype.sum_unique, if_pos (fun i => i.elim0)]
      rfl
  | succ n ih =>
      intro x
      rw [sum_pi_fin_succ]
      simp only [confinedMass]
      refine Finset.sum_congr rfl fun z _ => ?_
      have hcond : ∀ q : Fin n → S,
          (∀ i, A ((Fin.cons z q : Fin (n + 1) → S) i)) ↔
            (A z ∧ ∀ i, A (q i)) := by
        intro q
        constructor
        · intro h
          exact ⟨by simpa using h 0,
            fun i => by simpa [Fin.cons_succ] using h i.succ⟩
        · rintro ⟨hz, hq⟩ i
          induction i using Fin.cases with
          | zero => simpa using hz
          | succ j => simpa [Fin.cons_succ] using hq j
      have hrw : (∑ q : Fin n → S,
            if (∀ i, A ((Fin.cons z q : Fin (n + 1) → S) i))
              then K.pathWeight (n + 1) x (Fin.cons z q) else 0)
          = ∑ q : Fin n → S,
              if A z then
                (if (∀ i, A (q i))
                  then K.probability x z * K.pathWeight n z q else 0)
              else 0 := by
        refine Finset.sum_congr rfl fun q _ => ?_
        rw [pathWeight_cons, if_congr (hcond q) rfl rfl, ite_and]
      rw [hrw]
      by_cases hz : A z
      · simp only [if_pos hz]
        rw [ih z, Finset.mul_sum]
        refine Finset.sum_congr rfl fun q _ => ?_
        rw [mul_ite, mul_zero]
      · simp [hz]

end FiniteMarkovKernel

/-! ## Persistent-noise negative control

A kernel satisfying the drift condition with `ξ > 0` for which every
one-time endpoint bound holds, with excursion probability exactly `1/2` at
each single time `n ≥ 1`, while the probability of an excursion-free
window of length `N` is `(1/2)^N → 0`.  Machine-checked witness that the
endpoint receipt does NOT imply pathwise confinement, exactly as the
manuscript's limitation clause states. -/

namespace PersistentNoiseControl

open FiniteMarkovKernel

/-- Fair-coin kernel on `Bool`: at every step, move to `true` or `false`
with probability `1/2` each, independently of the current state. -/
noncomputable def noisyKernel : FiniteMarkovKernel Bool where
  probability _ _ := 1 / 2
  probability_nonneg := fun _ _ => by norm_num
  probability_sum_one := fun _ => by
    rw [Fintype.sum_bool]
    norm_num

/-- Escape indicator: `1` on the escaped state `true`, `0` on `false`. -/
def escapeIndicator (b : Bool) : ℝ := if b then 1 else 0

theorem escapeIndicator_nonneg : ∀ b : Bool, 0 ≤ escapeIndicator b := by
  intro b
  unfold escapeIndicator
  split <;> norm_num

/-- The drift condition holds with `κ = 0` and persistent noise `ξ = 1/2`. -/
theorem noisy_drift :
    ∀ x : Bool, noisyKernel.apply escapeIndicator x ≤
      0 * escapeIndicator x + 1 / 2 := by
  intro x
  simp only [FiniteMarkovKernel.apply, noisyKernel]
  rw [Fintype.sum_bool]
  simp [escapeIndicator]

/-- For every `n ≥ 1` the endpoint law is uniform on `Bool`. -/
theorem noisy_endpoint_uniform :
    ∀ n : ℕ, n ≠ 0 → ∀ x y : Bool,
      noisyKernel.endpointDistribution n x y = 1 / 2 := by
  intro n
  induction n with
  | zero => intro h; exact absurd rfl h
  | succ n ih =>
      intro _ x y
      rcases Nat.eq_zero_or_pos n with hn | hn
      · subst hn
        simp only [endpointDistribution, noisyKernel]
        rw [Fintype.sum_bool]
        cases y <;> norm_num
      · have hne : n ≠ 0 := Nat.pos_iff_ne_zero.mp hn
        simp only [endpointDistribution]
        rw [Fintype.sum_bool, ih hne true y, ih hne false y]
        simp only [noisyKernel]
        norm_num

/-- **Nonvanishing one-time excursion probability.**  At every single time
`n ≥ 1`, `P_x[V(X_n) > 1/2] = 1/2` exactly, so the one-time endpoint bound
holds at a value that is constant in `n`. -/
theorem noisy_endpoint_tail_const (n : ℕ) (hn : n ≠ 0) (x : Bool) :
    finiteEventMass (noisyKernel.endpointDistribution n x)
        (fun y => (1 : ℝ) / 2 < escapeIndicator y) = 1 / 2 := by
  unfold finiteEventMass
  rw [Finset.sum_filter, Fintype.sum_bool]
  rw [if_pos (by norm_num [escapeIndicator] :
        (1 : ℝ) / 2 < escapeIndicator true)]
  rw [if_neg (by norm_num [escapeIndicator] :
        ¬ (1 : ℝ) / 2 < escapeIndicator false)]
  rw [noisy_endpoint_uniform n hn x true, add_zero]

/-- **Window confinement collapses.**  The probability that the first `N`
states all stay inside the `δ = 1/2` tube is exactly `(1/2)^N`. -/
theorem noisy_confinedMass :
    ∀ (n : ℕ) (x : Bool),
      noisyKernel.confinedMass (fun b => escapeIndicator b ≤ 1 / 2) n x =
        (1 / 2) ^ n := by
  intro n
  induction n with
  | zero =>
      intro x
      simp [FiniteMarkovKernel.confinedMass]
  | succ n ih =>
      intro x
      simp only [FiniteMarkovKernel.confinedMass]
      rw [Fintype.sum_bool]
      rw [if_neg (by norm_num [escapeIndicator] :
            ¬ escapeIndicator true ≤ (1 : ℝ) / 2)]
      rw [if_pos (by norm_num [escapeIndicator] :
            escapeIndicator false ≤ (1 : ℝ) / 2)]
      rw [ih false]
      simp only [noisyKernel]
      rw [pow_succ]
      ring

/-- Persistent noise permits repeated excursions: window confinement decays
to zero even though every one-time endpoint bound holds. -/
theorem noisy_confinement_vanishes (x : Bool) :
    Filter.Tendsto
      (fun n => noisyKernel.confinedMass
        (fun b => escapeIndicator b ≤ 1 / 2) n x)
      Filter.atTop (nhds 0) := by
  simp only [noisy_confinedMass]
  exact tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num) (by norm_num)

/-- Concrete form: below any positive threshold after finitely many steps. -/
theorem noisy_confinement_eventually_below {c : ℝ} (hc : 0 < c) (x : Bool) :
    ∃ N : ℕ, noisyKernel.confinedMass
      (fun b => escapeIndicator b ≤ 1 / 2) N x < c := by
  obtain ⟨N, hN⟩ := exists_pow_lt_of_lt_one hc (by norm_num : (1 : ℝ) / 2 < 1)
  exact ⟨N, by rw [noisy_confinedMass]; exact hN⟩

end PersistentNoiseControl

/-! ## Support-checker audit

The positive-support hypothesis of the receipt is a genuine checker: it must
REJECT declarations whose rewrite relation does not actually cover the
kernel's moving support.  Both rejection tests and a passing control are
theorems. -/

namespace SupportAudit

open FiniteMarkovKernel

/-- Deterministic swap kernel on `Bool`: always moves to the other state. -/
def swapKernel : FiniteMarkovKernel Bool where
  probability x y := if y = !x then 1 else 0
  probability_nonneg := fun x y => by
    split <;> norm_num
  probability_sum_one := fun x => by
    rw [Fintype.sum_bool]
    cases x <;> norm_num

/-- Pure-stutter kernel on `Bool`: never moves. -/
def stutterKernel : FiniteMarkovKernel Bool where
  probability x y := if y = x then 1 else 0
  probability_nonneg := fun x y => by
    split <;> norm_num
  probability_sum_one := fun x => by
    rw [Fintype.sum_bool]
    cases x <;> norm_num

/-- **Rejection test (empty declaration).**  The empty rewrite relation does
NOT pass the support check against a kernel that actually moves. -/
theorem empty_support_declaration_fails :
    ¬ ∀ ⦃a b : Bool⦄,
        0 < PersistentNoiseControl.noisyKernel.probability a b →
          (fun _ _ => False) a b ∨ b = a := by
  intro h
  have hmove : (0 : ℝ) <
      PersistentNoiseControl.noisyKernel.probability false true := by
    simp only [PersistentNoiseControl.noisyKernel]
    norm_num
  rcases h hmove with hF | hEq
  · exact hF
  · exact Bool.noConfusion hEq

/-- **Rejection test (disconnected declaration).**  A nonempty rewrite
relation touching only a self-loop at `true`, disconnected from the swap
kernel's actual support, does NOT pass the support check. -/
theorem disconnected_support_declaration_fails :
    ¬ ∀ ⦃a b : Bool⦄, 0 < swapKernel.probability a b →
        (fun a b => a = true ∧ b = true) a b ∨ b = a := by
  intro h
  have hmove : (0 : ℝ) < swapKernel.probability false true := by
    simp only [swapKernel]
    norm_num
  rcases h hmove with ⟨hfa, _⟩ | hEq
  · exact Bool.noConfusion hfa
  · exact Bool.noConfusion hEq

/-- **Passing control (stutter support).**  The pure-stutter kernel passes
the support check even with the empty rewrite relation: all its support is
stuttering. -/
theorem stutter_support_declaration_holds :
    ∀ ⦃a b : Bool⦄, 0 < stutterKernel.probability a b →
      (fun _ _ => False) a b ∨ b = a := by
  intro a b h
  right
  by_contra hne
  have hzero : stutterKernel.probability a b = 0 := by
    simp [stutterKernel, hne]
  rw [hzero] at h
  exact lt_irrefl 0 h

/-- End-to-end: on the stutter kernel, the receipt's observation clause
holds for the IDENTITY observation: every state with positive endpoint mass
is literally the initial state's observation class. -/
example (n : ℕ) (x y : Bool)
    (hy : 0 < stutterKernel.endpointDistribution n x y) :
    (id : Bool → Bool) y = id x :=
  stutterKernel.observation_eq_of_endpointDistribution_pos
    (r := fun _ _ => False)
    stutter_support_declaration_holds
    (fun _ _ h => h.elim)
    hy

end SupportAudit

end ObservableNormalForms
