import Mathlib

namespace SourceLaw.DaughterShell

/-- Each ordinary neutral word has at most three neighbours of weight seven
and all its other neighbours have weight one. -/
theorem ordinary_row (potential degree highNeighbours : ℝ)
    (hv : potential ≤ 112) (hd : degree ≤ 12) (hh : highNeighbours ≤ 3) :
    potential + degree / 2 + (7 * highNeighbours + degree - highNeighbours) / 2 ≤ 133 := by
  linarith

/-- A high neutral word has ten different-load swaps and at most two high neighbours. -/
theorem high_row (highNeighbours : ℝ) (hh : highNeighbours ≤ 2) :
    128 + 10 / 2 + (7 * highNeighbours + 10 - highNeighbours) / 14 ≤ 942 / 7 := by
  linarith

/-- The weakest of the three finite parent-band floors suffices on all boards. -/
theorem energy_deficit (parent charged neutral : ℝ)
    (hp : 397 / 2 ≤ parent) (hc : charged ≤ 209 / 4) (hn : neutral ≤ 942 / 7) :
    327 / 28 ≤ parent - (charged + neutral) := by
  linarith

/-- A stationary transfer between parent and independent-outgoing energy bases
has zero amplitude in every certified upper-band/daughter channel. -/
theorem stationary_amplitude (parent charged neutral : ℝ) (amplitude : ℂ)
    (hp : 397 / 2 ≤ parent) (hc : charged ≤ 209 / 4) (hn : neutral ≤ 942 / 7)
    (stationary : ((charged + neutral : ℝ) : ℂ) * amplitude = amplitude * (parent : ℂ)) :
    amplitude = 0 := by
  have h : charged + neutral ≠ parent := by
    have := energy_deficit parent charged neutral hp hc hn
    linarith
  have hc' : ((charged + neutral : ℝ) : ℂ) ≠ (parent : ℂ) := by exact_mod_cast h
  have zero : (((charged + neutral : ℝ) : ℂ) - (parent : ℂ)) * amplitude = 0 := by
    rw [sub_mul, stationary, mul_comm amplitude, sub_self]
  exact (mul_eq_zero.mp zero).resolve_left (sub_ne_zero.mpr hc')

end SourceLaw.DaughterShell
