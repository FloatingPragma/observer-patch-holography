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

#print axioms signed_sandwich
#print axioms odd_compression
#print axioms one_mode_width

end OPH.FreeEnd
