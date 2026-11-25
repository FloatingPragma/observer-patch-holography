import Mathlib.Tactic

namespace MaterialRecoil

/-- Averaging three independent fair contact switches. Signs describe the
    coherence of one test load in the left/right material basis. -/
theorem switch_average (x y z : ℝ) :
    (1 + x + y + z + x*y + x*z + y*z + x*y*z)/8 =
    ((1+x)/2)*((1+y)/2)*((1+z)/2) := by ring

theorem active_phase_filter (a b c : Bool) :
    let sign : Bool → ℝ := fun hit => if hit then 1 else -1
    ((1+sign a)/2)*((1+sign b)/2)*((1+sign c)/2) =
      if a && b && c then 1 else 0 := by
  cases a <;> cases b <;> cases c <;> norm_num

/-- The free flight time is a half-period. Equal mirror weights divide
    the native displacement response by two. -/
theorem current_to_load (near far : ℝ) :
    let d : ℝ := 65/2
    d/(2*(1+d^2))*(far-near) = -(520/4229)*((near-far)/8) := by
  dsimp
  ring

/-- The two bare diagonal energies on a board of arbitrary native degree
    differ by one common constant. The off-diagonal swap is shared. -/
theorem bare_shift (k deg : ℝ) :
    ((6-deg)*k*2^2 + deg*(k*2^2+1/2)) -
    ((6-deg)*k + deg*(k+1/2)) = 18*k := by ring

end MaterialRecoil
