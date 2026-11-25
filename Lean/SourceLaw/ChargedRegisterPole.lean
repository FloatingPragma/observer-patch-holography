import Mathlib

/-!
# Shared spectral weight in a charged register

An energy detector has a common outcome on three equally likely preparations.
Their weights at that outcome follow from the two-tip block of the material
compact projector. The return experiment supplies the physical source and the
finite filter. These bounds concern discrimination of one preparation from a
single energy reading.
-/

namespace OPH.SourceLaw.ChargedRegisterPole

/-- A common spectral outcome with probability at least one ninth in
each of three channels prevents perfect single-trial discrimination. -/
theorem shared_outcome_error {a b c : ℝ} (ha : 1 / 9 ≤ a)
    (hb : 1 / 9 ≤ b) (hc : 1 / 9 ≤ c) :
    1 - (a + b + c - max a (max b c)) / 3 ≤ 25 / 27 := by
  rcases le_total b c with h | h
  · rw [max_eq_right h]
    rcases le_total a c with h' | h'
    · rw [max_eq_right h']; linarith
    · rw [max_eq_left h']; linarith
  · rw [max_eq_left h]
    rcases le_total a b with h' | h'
    · rw [max_eq_right h']; linarith
    · rw [max_eq_left h']; linarith

/-- The two-tip projector supplies weight (2-cos t)/9 to each character. -/
theorem two_tip_weight (t : ℝ) :
    1 / 9 ≤ (2 - Real.cos t) / 9 ∧ (2 - Real.cos t) / 9 ≤ 1 / 3 := by
  constructor
  · linarith [Real.cos_le_one t]
  · linarith [Real.neg_one_le_cos t]

/-- The three equally weighted material preparations have common-line
weights summing to two thirds. A correct guess on that line collects a numerator
of at most one third before averaging over the equal priors. -/
theorem marked_edge_ceiling {a b c success : ℝ}
    (sum : a + b + c = 2 / 3) (ha : a ≤ 1 / 3)
    (hb : b ≤ 1 / 3) (hc : c ≤ 1 / 3)
    (hguess : success ≤ 1 - (a + b + c - max a (max b c)) / 3) :
    success ≤ 8 / 9 := by
  have hm : max a (max b c) ≤ 1 / 3 := max_le ha (max_le hb hc)
  linarith

/-- A positive filter with unit weight at the reference pole and tail
bounded by epsilon encloses the true pole weight, including shot error. -/
theorem measured_weight_interval {w tail eps shot estimate : ℝ}
    (ht0 : 0 ≤ tail) (ht : tail ≤ eps)
    (low : w + tail - shot ≤ estimate) (high : estimate ≤ w + tail + shot) :
    estimate - shot - eps ≤ w ∧ w ≤ estimate + shot := by
  constructor <;> linarith

/-- On ordered singular rails the physical attenuation lies in the unit
interval and leaves an equal complex amplitude on both qubit coordinates. -/
theorem attenuated_qubit {s₁ s₂ : ℝ} (h₁ : 0 ≤ s₁) (h₂ : 0 < s₂)
    (h : s₁ ≤ s₂) (z₁ z₂ : ℂ) :
    0 ≤ s₁ / s₂ ∧ s₁ / s₂ ≤ 1 ∧
      ((s₁ : ℂ) * z₁, ((s₁ / s₂ : ℝ) : ℂ) * ((s₂ : ℂ) * z₂)) =
      (((min s₁ s₂ : ℝ) : ℂ) * z₁, ((min s₁ s₂ : ℝ) : ℂ) * z₂) := by
  refine ⟨div_nonneg h₁ h₂.le, (div_le_one h₂).2 h, ?_⟩
  rw [min_eq_left h, ← mul_assoc, ← Complex.ofReal_mul, div_mul_cancel₀ s₁ h₂.ne']

/-- A transitive port-permutation action with nonzero phases has no
nonzero invariant subspace of words vanishing at an observed gateway. -/
theorem magnetic_gateway_exclusion {I G : Type*} (gateway : I)
    (site : G → I → I) (phase : G → I → ℂ)
    (turn : G → (I → ℂ) → (I → ℂ)) (W : Set (I → ℂ))
    (reach : ∀ i, ∃ g, site g i = gateway)
    (unit : ∀ g i, phase g i ≠ 0)
    (monomial : ∀ g x i, turn g x (site g i) = phase g i * x i)
    (invariant : ∀ g x, x ∈ W → turn g x ∈ W)
    (dark : ∀ x ∈ W, x gateway = 0) {x : I → ℂ} (hx : x ∈ W) :
    x = 0 := by
  funext i
  obtain ⟨g, hg⟩ := reach i
  have h := dark (turn g x) (invariant g x hx)
  rw [← hg, monomial] at h
  exact (mul_eq_zero.mp h).resolve_left (unit g i)

/-- At a transcendental phase coordinate, a ratio of polynomial values
can equal a base-field scalar only if the two polynomials are proportional. -/
theorem phase_ratio_ne_base_scalar {K L : Type*} [Field K] [Field L] [Algebra K L]
    {z : L} (hz : Transcendental K z) {p q : Polynomial K} (a : K)
    (hp : p ≠ Polynomial.C a * q) (hq : Polynomial.aeval z q ≠ 0) :
    Polynomial.aeval z p / Polynomial.aeval z q ≠ algebraMap K L a := by
  intro h
  apply hp
  apply (transcendental_iff_injective.mp hz)
  simpa using (div_eq_iff hq).mp h

#print axioms phase_ratio_ne_base_scalar
#print axioms attenuated_qubit
#print axioms magnetic_gateway_exclusion
#print axioms shared_outcome_error
#print axioms two_tip_weight
#print axioms marked_edge_ceiling
#print axioms measured_weight_interval

end OPH.SourceLaw.ChargedRegisterPole
