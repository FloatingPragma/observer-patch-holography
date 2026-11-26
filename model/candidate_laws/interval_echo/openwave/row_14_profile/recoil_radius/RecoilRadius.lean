import Mathlib

namespace SignedCentroid

/-- Two excesses and two vacancies have equal-magnitude opposite pair
centroid displacements. Their signed second moment records pair sizes. -/
theorem coordinate_identity (a b c d : ℝ) :
    let x := (a+b+c+d)/4
    2*((a-x)^2 + (b-x)^2 - (c-x)^2 - (d-x)^2) =
      (a-b)^2 - (c-d)^2 := by
  dsimp
  ring

theorem summed_identity {n : ℕ} (a b c d : Fin n → ℝ) :
    let x := fun i => (a i + b i + c i + d i) / 4
    2 * (∑ i, ((a i-x i)^2 + (b i-x i)^2 - (c i-x i)^2 - (d i-x i)^2)) =
      ∑ i, ((a i-b i)^2 - (c i-d i)^2) := by
  dsimp
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun i _ => coordinate_identity (a i) (b i) (c i) (d i)

theorem inverse_order {neutral positive : ℝ} (hn : 0 < neutral)
    (hp : neutral < positive) : 1 / positive < 1 / neutral := by
  exact one_div_lt_one_div_of_lt hn hp

/-- Two independent stocks with impulses f u and (1-f) u. Their
minimal quadratic work follows by completing this square. -/
theorem split_work (p v f : ℝ) :
    (p+v)*(f^2*p + (1-f)^2*v) = p*v + ((p+v)*f-v)^2 := by ring

end SignedCentroid
