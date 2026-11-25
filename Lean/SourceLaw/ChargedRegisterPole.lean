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

#print axioms shared_outcome_error
#print axioms two_tip_weight
#print axioms marked_edge_ceiling
#print axioms measured_weight_interval

end OPH.SourceLaw.ChargedRegisterPole
