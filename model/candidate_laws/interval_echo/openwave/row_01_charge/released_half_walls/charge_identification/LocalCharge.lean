import Mathlib.Tactic

namespace IntervalEcho.LocalCharge

/-- Padded walks first write a vertex word from zero, use an edge, then return
to zero. The residual is the measured total minus the independent integer stock. -/
theorem padded_edge_identity (left right : ℝ) (prefixLeft prefixRight suffix edge idle pad : ℝ)
    (out : left = prefixLeft + edge + suffix + (pad - 1) * idle)
    (back : right = prefixRight + suffix + pad * idle) :
    edge - (prefixRight - prefixLeft) = left - right + idle := by
  rw [out, back]
  ring

/-- Padding and two measured calibrations bound each residual after a gradient
is removed. All three error bounds use ordinary integrated-load units. -/
theorem calibrated_edge_error (edge gradient left right idle δ ε : ℝ)
    (identity : edge - gradient = left - right + idle)
    (hl : |left| ≤ δ) (hr : |right| ≤ δ) (hz : |idle| ≤ ε) :
    |edge - gradient| ≤ 2 * δ + ε := by
  rw [identity]
  calc
    |left - right + idle| ≤ |left - right| + |idle| := abs_add_le ..
    _ ≤ |left| + |right| + |idle| := by gcongr; exact abs_sub _ _
    _ ≤ 2 * δ + ε := by linarith

/-- A gradient in the local-word graph telescopes on any closed stock. -/
theorem whole_ring_error {V : Type*} [Fintype V] (next : V ≃ V)
    (reading standard pathSum : V → ℝ) (ε : ℝ)
    (localBound : ∀ j, |(reading j - standard j) - (pathSum (next j) - pathSum j)| ≤ ε) :
    |∑ j, reading j - ∑ j, standard j| ≤ Fintype.card V * ε := by
  have perm : ∑ j, pathSum (next j) = ∑ j, pathSum j := Equiv.sum_comp next pathSum
  let d := fun j => (reading j - standard j) - (pathSum (next j) - pathSum j)
  have h : ∑ j, d j = ∑ j, reading j - ∑ j, standard j := by
    simp only [d, Finset.sum_sub_distrib, perm, sub_self, sub_zero]
  rw [← h]
  calc
    |∑ j, d j| ≤ ∑ j, |d j| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _ : V, ε := Finset.sum_le_sum fun j _ => localBound j
    _ = _ := by simp

/-- Equal local word histograms give equal integrated readings for *any*
possibly nonlinear lookup table evaluated at those sharp words. -/
theorem equal_patch {W : Type*} [Fintype W] (left right : W → ℕ) (f : W → ℝ)
    (same : ∀ w, left w = right w) :
    ∑ w, (left w : ℝ) * f w = ∑ w, (right w : ℝ) * f w := by
  simp_rw [same]

/-- Integrated stock error, residual vacuum response, and asymmetry of the
two tested wall patches have independent contributions to a half-charge bound. -/
theorem half_response_error (left right rest inserted δvac δsum δpair : ℝ)
    (stock : |left + right + rest - inserted| ≤ δsum)
    (vacuum : |rest| ≤ δvac) (matched : |left - right| ≤ δpair) :
    |left - inserted / 2| ≤ (δsum + δvac + δpair) / 2 := by
  rw [abs_le] at *
  constructor <;> linarith

theorem half_response (left right rest inserted : ℝ)
    (stock : left + right + rest = inserted) (vacuum : rest = 0)
    (matched : left = right) : left = inserted / 2 ∧ right = inserted / 2 := by
  constructor <;> linarith

end IntervalEcho.LocalCharge

#print axioms IntervalEcho.LocalCharge.padded_edge_identity
#print axioms IntervalEcho.LocalCharge.calibrated_edge_error
#print axioms IntervalEcho.LocalCharge.whole_ring_error
#print axioms IntervalEcho.LocalCharge.equal_patch
#print axioms IntervalEcho.LocalCharge.half_response_error
#print axioms IntervalEcho.LocalCharge.half_response
