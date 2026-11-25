import Mathlib

namespace SourceLaw.WheelPole

/-- The five hopping amplitudes into the central port cancel exactly. -/
theorem dark_gateway (z : ℂ) (period : z^5 = 1) (nontrivial : z ≠ 1) :
    1 + z + z^2 + z^3 + z^4 = 0 := by
  have prod : (z-1) * (1 + z + z^2 + z^3 + z^4) = 0 := by
    linear_combination period
  exact (mul_eq_zero.mp prod).resolve_left (sub_ne_zero.mpr nontrivial)

/-- Clockwise and counterclockwise rim neighbors enter with the same
native coefficient. Their sum acts diagonally on either circular wave. -/
theorem rim_energy (z w a : ℂ) :
    (51/2 : ℂ)*a - (z*a + w*a)/2 = ((51-z-w)/2)*a := by ring

/-- No nontrivial fifth root equals its own inverse. The opposite characters
therefore label distinct waves, even though the real rim kernel pairs their energies. -/
theorem opposite_distinct (z w : ℂ) (period : z^5 = 1) (unit : z*w = 1)
    (nontrivial : z ≠ 1) : z ≠ w := by
  intro equal
  have sq : z^2 = 1 := by simpa [← equal, pow_two] using unit
  apply nontrivial
  calc
    z = z * (z^2)^2 := by rw [sq]; ring
    _ = z^5 := by ring
    _ = 1 := period

end SourceLaw.WheelPole
