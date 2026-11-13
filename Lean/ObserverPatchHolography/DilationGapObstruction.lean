import Mathlib.Data.Real.Basic
import Mathlib.Tactic

/-! Necessary energy inequalities for a common vacuum and clock.

The first theorem consumes a nonzero finite-energy excitation whose energy
can be scaled by every positive real factor while the same positive gap
continues to apply. A vacuum-preserving unitary dilation covariance supplies
that hypothesis for a flat-space conformal vacuum theory. No such hypothesis
is asserted for a compact spatial cylinder or a scale-breaking vacuum.

The second theorem records the different obstruction of a fixed dimensionless
source gap under inverse-spacing physical-clock conversion. Neither theorem
constructs a quantum field theory or supplies Yang--Mills continuum premises.
-/
namespace OPH.DilationGapObstruction

theorem scaled_finite_excitation_excludes_positive_gap
    (gap energy : ℝ) (hgap : 0 < gap) (henergy : 0 < energy)
    (hscaled : ∀ scale : ℝ, 0 < scale → gap ≤ scale * energy) : False := by
  have hscale : 0 < gap / (2 * energy) := div_pos hgap (by positivity)
  have hbound := hscaled (gap / (2 * energy)) hscale
  have hvalue : gap / (2 * energy) * energy = gap / 2 := by
    field_simp
  rw [hvalue] at hbound
  linarith

theorem inverse_spacing_gap_exceeds_finite_target
    (gap spacing target : ℝ) (hspacing : 0 < spacing) (htarget : 0 < target)
    (hsmall : spacing < gap / target) : target < gap / spacing := by
  apply (lt_div_iff₀ hspacing).2
  have hproduct := (lt_div_iff₀ htarget).1 hsmall
  nlinarith

#print axioms scaled_finite_excitation_excludes_positive_gap
#print axioms inverse_spacing_gap_exceeds_finite_target

end OPH.DilationGapObstruction
