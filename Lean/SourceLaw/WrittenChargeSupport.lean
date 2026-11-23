import Mathlib

/-! A zero of a positive contact and ring Laplacian with a nonnegative
occupation cost is constant along every enabled transition and vanishes on
every state with positive cost. Real and imaginary parts give the same
conclusion for complex amplitudes. -/
namespace SourceLaw.WrittenChargeSupport

theorem zero_energy_support {V E : Type*} [Fintype V] [Fintype E]
    (left right : E → V) (weight : E → ℝ) (cost amp : V → ℝ)
    (hw : ∀ e, 0 < weight e) (hc : ∀ v, 0 ≤ cost v)
    (zero : (∑ e, weight e * (amp (left e) - amp (right e))^2) +
      (∑ v, cost v * (amp v)^2) = 0) :
    (∀ e, amp (left e) = amp (right e)) ∧
      (∀ v, 0 < cost v → amp v = 0) := by
  classical
  have he (e : E) : 0 ≤ weight e * (amp (left e) - amp (right e))^2 :=
    mul_nonneg (hw e).le (sq_nonneg _)
  have hv (v : V) : 0 ≤ cost v * (amp v)^2 := mul_nonneg (hc v) (sq_nonneg _)
  have hE : 0 ≤ ∑ e, weight e * (amp (left e) - amp (right e))^2 :=
    Finset.sum_nonneg fun e _ => he e
  have hV : 0 ≤ ∑ v, cost v * (amp v)^2 := Finset.sum_nonneg fun v _ => hv v
  constructor
  · intro e
    have le := Finset.single_le_sum (f := fun i => weight i * (amp (left i) - amp (right i))^2)
      (fun i _ => he i) (Finset.mem_univ e)
    have h0 : weight e * (amp (left e) - amp (right e))^2 = 0 := le_antisymm (by linarith) (he e)
    have hsq := (mul_eq_zero.mp h0).resolve_left (ne_of_gt (hw e))
    exact sub_eq_zero.mp (sq_eq_zero_iff.mp hsq)
  · intro v positive
    have le := Finset.single_le_sum (f := fun i => cost i * (amp i)^2)
      (fun i _ => hv i) (Finset.mem_univ v)
    have h0 : cost v * (amp v)^2 = 0 := le_antisymm (by linarith) (hv v)
    have hsq := (mul_eq_zero.mp h0).resolve_left (ne_of_gt positive)
    exact sq_eq_zero_iff.mp hsq

end SourceLaw.WrittenChargeSupport
