import Mathlib.Tactic

namespace OPH.FreeEnd

/-- If a dichotomic analyzer and its evolved observable have a scalar
anticommutator, projecting onto the prepared sign fixes the second sign's
mean. For a negative sign, apply this to the negatively oriented analyzer. -/
theorem signed_sandwich {A : Type*} [Ring A] (p y t c : A)
    (hs : p*y = p) (hs' : y*p = p) (hp : p*p = p)
    (hac : y*t+t*y = c+c) (hc : p*c = c*p) :
    p*t*p+p*t*p = c*p+c*p := by
  calc
    _ = (p*y)*t*p+p*t*(y*p) := by rw [hs, hs']
    _ = p*(y*t+t*y)*p := by noncomm_ring
    _ = p*c*p+p*c*p := by rw [hac]; noncomm_ring
    _ = c*p+c*p := by rw [hc, mul_assoc, hp]

/-- A fixed-parity herald annihilates the compressed odd current over the
complex matrices, where an additive element equal to its negative is zero. -/
theorem odd_compression {A : Type*} [Ring A] (p b q : A)
    (hleft : p*q = p) (hright : q*p = p) (hodd : q*b = -(b*q)) :
    p*b*p = -(p*b*p) := by
  calc
    _ = (p*q)*b*p := by rw [hleft]
    _ = p*(q*b)*p := by rw [mul_assoc p q b]
    _ = -(p*b*(q*p)) := by rw [hodd]; noncomm_ring
    _ = _ := by rw [hright]

/-- The single-mode sector can carry only difference frequencies in
a band of width four hopping units, independent of its absolute energy. -/
theorem one_mode_width (j h a b : ℝ)
    (ha : 2*(j-h) ≤ a ∧ a ≤ 2*(j+h)) (hb : 2*(j-h) ≤ b ∧ b ≤ 2*(j+h)) :
    |a-b| ≤ 4*h := by
  rcases ha with ⟨_, _⟩
  rcases hb with ⟨_, _⟩
  exact abs_le.mpr ⟨by linarith, by linarith⟩

/-- A product state has fidelity at most one half with the two equal
Schmidt branches of routing an excitation coherent with its idle copy.
Here a,b and c,d bound absolute overlaps on the two sides of the cut. -/
theorem switched_product_fidelity (a b c d : ℝ)
    (hx : a^2+b^2 ≤ 1) (hy : c^2+d^2 ≤ 1) :
    (a*c+b*d)^2/2 ≤ 1/2 := by
  have hc : (a*c+b*d)^2 ≤ (a^2+b^2)*(c^2+d^2) := by
    nlinarith [sq_nonneg (a*d-b*c)]
  have hb : (a^2+b^2)*(c^2+d^2) ≤ 1*1 :=
    mul_le_mul hx hy (by positivity) (by norm_num)
  linarith

#print axioms switched_product_fidelity

#print axioms signed_sandwich
#print axioms odd_compression
#print axioms one_mode_width

end OPH.FreeEnd
