import Mathlib

namespace OPH.DiscreteRefinement

/-!
# Exact discrete-refinement lemmas

This file formalizes the inexpensive kernel of the issue-656 theorem packet.
The analytic DR-2 step is represented by `bi_refinement_periodic_factor_constant`:
a continuous normalized factor invariant under two shifts whose ratio is
irrational is constant.  The physical existence of those two shifts and their
action on one isolated covariance ray is deliberately not asserted here.
-/

/-- The nonnegative integer lattice consequence of one exact scale relation. -/
theorem lattice_law_nat {M N : Type*} [CommMonoid M] [CommMonoid N]
    (F : M → N) (B k : M) (lambda : N)
    (hscale : ∀ x, F (B * x) = lambda * F x) :
    ∀ m : ℕ, F (B ^ m * k) = lambda ^ m * F k := by
  intro m
  induction m with
  | zero => simp
  | succ m ih =>
      calc
        F (B ^ (m + 1) * k) = F (B * (B ^ m * k)) := by
          congr 1
          rw [pow_succ]
          ac_rfl
        _ = lambda * F (B ^ m * k) := hscale _
        _ = lambda * (lambda ^ m * F k) := by rw [ih]
        _ = lambda ^ (m + 1) * F k := by
          rw [pow_succ]
          ac_rfl

/-- A periodic factor gives the exact additive-log version of DR-1. -/
theorem periodic_normal_form_scales (theta a : ℝ) (p : ℝ → ℝ)
    (hperiod : ∀ t, p (t + a) = p t) (t : ℝ) :
    Real.exp (-theta * (t + a)) * p (t + a) =
      Real.exp (-theta * a) * (Real.exp (-theta * t) * p t) := by
  rw [hperiod]
  rw [show -theta * (t + a) = (-theta * a) + (-theta * t) by ring]
  rw [Real.exp_add]
  ring

/-- The additive subgroup of exact periods of a real-valued function. -/
def periodSubgroup (p : ℝ → ℝ) : AddSubgroup ℝ where
  carrier := {a | ∀ x, p (x + a) = p x}
  zero_mem' := by simp
  add_mem' := by
    intro a b ha hb x
    calc
      p (x + (a + b)) = p ((x + a) + b) := by ring_nf
      _ = p (x + a) := hb (x + a)
      _ = p x := ha x
  neg_mem' := by
    intro a ha x
    have h := ha (x - a)
    have h' : p x = p (x - a) := by simpa using h
    simpa [sub_eq_add_neg] using h'.symm

/-- Invariance under each of two shifts extends to their generated subgroup. -/
theorem closure_pair_le_periods (p : ℝ → ℝ) {a b : ℝ}
    (ha : ∀ x, p (x + a) = p x) (hb : ∀ x, p (x + b) = p x) :
    AddSubgroup.closure ({a, b} : Set ℝ) ≤ periodSubgroup p := by
  rw [AddSubgroup.closure_le]
  intro x hx
  rcases hx with (rfl | hx)
  · exact ha
  · simpa using hx ▸ hb

/-- A continuous function invariant under a dense subgroup of shifts is constant. -/
theorem continuous_dense_periods_constant (p : ℝ → ℝ) (hp : Continuous p)
    (S : AddSubgroup ℝ) (hS : Dense (S : Set ℝ))
    (hperiods : S ≤ periodSubgroup p) :
    ∃ A : ℝ, ∀ x, p x = A := by
  refine ⟨p 0, ?_⟩
  intro x
  let f : ℝ → ℝ := fun s => p (x + s)
  let g : ℝ → ℝ := fun _ => p x
  have hf : Continuous f := hp.comp (continuous_const.add continuous_id)
  have hg : Continuous g := continuous_const
  have heq : Set.EqOn f g (S : Set ℝ) := by
    intro s hs
    exact hperiods hs x
  have hfg : f = g := Continuous.ext_on hS hf hg heq
  have h_at := congrFun hfg (-x)
  have hzero : p 0 = p x := by simpa [f, g] using h_at
  exact hzero.symm

/-- Iteration of a multiplicative shift relation. -/
theorem shift_multiplier_iterate (p : ℝ → ℝ) (b mu : ℝ)
    (hshift : ∀ x, p (x + b) = mu * p x) :
    ∀ (n : ℕ) (x : ℝ), p (x + n * b) = mu ^ n * p x := by
  intro n x
  induction n with
  | zero => simp
  | succ n ih =>
      simp only [Nat.cast_succ]
      calc
        p (x + (n + 1) * b) = p ((x + n * b) + b) := by
          congr 1
          ring
        _ = mu * p (x + n * b) := hshift _
        _ = mu * (mu ^ n * p x) := by rw [ih]
        _ = mu ^ (n + 1) * p x := by rw [pow_succ]; ring

/-- A uniformly positive and bounded profile cannot acquire a multiplier
strictly larger than one under an exact shift. -/
theorem shift_multiplier_not_gt_one_of_bounds (p : ℝ → ℝ) (b mu c C : ℝ)
    (hc : 0 < c) (hmu : 0 < mu)
    (hbounds : ∀ x, c ≤ p x ∧ p x ≤ C)
    (hshift : ∀ x, p (x + b) = mu * p x) :
    ¬1 < mu := by
  intro hmuOne
  obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt (C / c) hmuOne
  have hlarge : C < mu ^ n * c := (div_lt_iff₀ hc).mp hn
  have hlower : mu ^ n * c ≤ mu ^ n * p 0 :=
    mul_le_mul_of_nonneg_left (hbounds 0).1 (pow_nonneg hmu.le n)
  have hiter := shift_multiplier_iterate p b mu hshift n 0
  have hupper := (hbounds (n * b)).2
  simp only [zero_add] at hiter
  rw [hiter] at hupper
  exact (not_lt_of_ge hupper) (hlarge.trans_le hlower)

/-- Compact bounded positivity forces the multiplier in the second DR-2
shift to be one.  Periodicity and continuity supply these bounds in the
analytic proof; the bounds are explicit premises in this kernel lemma. -/
theorem positive_bounded_shift_multiplier_eq_one (p : ℝ → ℝ) (b mu c C : ℝ)
    (hc : 0 < c) (hmu : 0 < mu)
    (hbounds : ∀ x, c ≤ p x ∧ p x ≤ C)
    (hshift : ∀ x, p (x + b) = mu * p x) :
    mu = 1 := by
  rcases lt_trichotomy mu 1 with hlt | heq | hgt
  · have hback : ∀ x, p (x + (-b)) = mu⁻¹ * p x := by
      intro x
      apply (eq_inv_mul_iff_mul_eq₀ hmu.ne').2
      have shifted := hshift (x - b)
      simpa [sub_eq_add_neg] using shifted.symm
    have hinv : 1 < mu⁻¹ := (one_lt_inv₀ hmu).2 hlt
    exact False.elim
      (shift_multiplier_not_gt_one_of_bounds p (-b) mu⁻¹ c C hc (inv_pos.mpr hmu)
        hbounds hback hinv)
  · exact heq
  · exact False.elim
      (shift_multiplier_not_gt_one_of_bounds p b mu c C hc hmu hbounds hshift hgt)

/-- DR-2, normalized form: two incommensurate exact periods force the
periodic factor to be constant. -/
theorem bi_refinement_periodic_factor_constant (p : ℝ → ℝ) (hp : Continuous p)
    {a b : ℝ} (hirr : Irrational (a / b))
    (ha : ∀ x, p (x + a) = p x) (hb : ∀ x, p (x + b) = p x) :
    ∃ A : ℝ, ∀ x, p x = A := by
  apply continuous_dense_periods_constant p hp (AddSubgroup.closure ({a, b} : Set ℝ))
  · exact dense_addSubgroupClosure_pair_iff.mpr hirr
  · exact closure_pair_le_periods p ha hb

/-- DR-2 compatibility and rigidity on a normalized factor.  The second raw
refinement first acts with a positive multiplier `mu`.  Uniform bounds force
that multiplier to one; irrational shift ratio and continuity then force the
factor to be constant. -/
theorem bi_refinement_multiplier_and_shape (p : ℝ → ℝ) (hp : Continuous p)
    {a b mu c C : ℝ} (hirr : Irrational (a / b))
    (hc : 0 < c) (hmu : 0 < mu)
    (hbounds : ∀ x, c ≤ p x ∧ p x ≤ C)
    (ha : ∀ x, p (x + a) = p x)
    (hb : ∀ x, p (x + b) = mu * p x) :
    mu = 1 ∧ ∃ A : ℝ, ∀ x, p x = A := by
  have hmuOne := positive_bounded_shift_multiplier_eq_one p b mu c C hc hmu hbounds hb
  refine ⟨hmuOne, ?_⟩
  apply bi_refinement_periodic_factor_constant p hp hirr ha
  intro x
  simpa [hmuOne] using hb x

/-- DR-2 collapse of a supplied one-ratio normal form after the normalized
factor has two incommensurate exact periods. -/
theorem bi_refinement_normal_form_is_pure_power (F p : ℝ → ℝ) (theta : ℝ)
    (hp : Continuous p) {a b : ℝ} (hirr : Irrational (a / b))
    (ha : ∀ x, p (x + a) = p x) (hb : ∀ x, p (x + b) = p x)
    (hF : ∀ x, F x = Real.exp (-theta * x) * p x) :
    ∃ A : ℝ, ∀ x, F x = A * Real.exp (-theta * x) := by
  obtain ⟨A, hA⟩ := bi_refinement_periodic_factor_constant p hp hirr ha hb
  refine ⟨A, fun x => ?_⟩
  rw [hF, hA]
  ring

/-- Unique factorization rules out a positive binary/ternary power resonance. -/
theorem two_pow_ne_three_pow (p q : ℕ) (hp : 0 < p) (hq : 0 < q) :
    2 ^ q ≠ 3 ^ p := by
  intro equality
  have hdivThree : 3 ∣ 2 ^ q := by
    rw [equality]
    exact dvd_pow_self 3 hp.ne'
  have hdivTwo : 2 ∣ 3 ^ p := by
    rw [← equality]
    exact dvd_pow_self 2 hq.ne'
  have impossible : 3 ∣ 2 ∧ 2 ∣ 3 := by
    exact ⟨
      Nat.Prime.dvd_of_dvd_pow (by norm_num) hdivThree,
      Nat.Prime.dvd_of_dvd_pow (by norm_num) hdivTwo,
    ⟩
  norm_num at impossible

/-- A homomorphism from a finite group to strictly positive real
multiplication is trivial.  This is the finite-group-to-scale no-go. -/
theorem finite_group_to_positive_reals_trivial {G : Type*} [Group G] [Finite G]
    (phi : G →* ℝ) (hpositive : ∀ g, 0 < phi g) :
    ∀ g, phi g = 1 := by
  intro g
  have horder : orderOf g ≠ 0 := (isOfFinOrder_of_finite g).orderOf_pos.ne'
  have hpow : (phi g) ^ orderOf g = 1 := by
    rw [← map_pow, pow_orderOf_eq_one, map_one]
  exact (pow_eq_one_iff_of_nonneg (hpositive g).le horder).mp hpow

/-- Frequency-n icosahedral divisibility-tower counts. -/
def faceCount (n : ℕ) : ℕ := 20 * n ^ 2
def edgeCount (n : ℕ) : ℕ := 30 * n ^ 2
def vertexCount (n : ℕ) : ℕ := 10 * n ^ 2 + 2

theorem tower_counts_one :
    (vertexCount 1, edgeCount 1, faceCount 1) = (12, 30, 20) := by
  norm_num [vertexCount, edgeCount, faceCount]

theorem tower_counts_two :
    (vertexCount 2, edgeCount 2, faceCount 2) = (42, 120, 80) := by
  norm_num [vertexCount, edgeCount, faceCount]

theorem tower_counts_three :
    (vertexCount 3, edgeCount 3, faceCount 3) = (92, 270, 180) := by
  norm_num [vertexCount, edgeCount, faceCount]

theorem tower_counts_six :
    (vertexCount 6, edgeCount 6, faceCount 6) = (362, 1080, 720) := by
  norm_num [vertexCount, edgeCount, faceCount]

/-- Same-parent barycentric coordinates for the divisibility mesh. -/
structure Barycentric where
  i : ℕ
  j : ℕ
  k : ℕ
  deriving DecidableEq

/-- Denominator multiplication on one parent-face coordinate system. -/
def refine (m : ℕ) (x : Barycentric) : Barycentric :=
  ⟨m * x.i, m * x.j, m * x.k⟩

theorem refinement_composes (m n : ℕ) (x : Barycentric) :
    refine m (refine n x) = refine (m * n) x := by
  cases x
  simp [refine, Nat.mul_assoc]

theorem binary_ternary_square (x : Barycentric) :
    refine 2 (refine 3 x) = refine 3 (refine 2 x) := by
  rw [refinement_composes, refinement_composes]

#print axioms lattice_law_nat
#print axioms periodic_normal_form_scales
#print axioms bi_refinement_periodic_factor_constant
#print axioms positive_bounded_shift_multiplier_eq_one
#print axioms bi_refinement_multiplier_and_shape
#print axioms bi_refinement_normal_form_is_pure_power
#print axioms two_pow_ne_three_pow
#print axioms finite_group_to_positive_reals_trivial
#print axioms tower_counts_six
#print axioms binary_ternary_square

end OPH.DiscreteRefinement
