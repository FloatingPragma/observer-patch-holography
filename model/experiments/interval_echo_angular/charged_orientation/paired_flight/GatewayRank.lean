import Mathlib

/-! Finite receiver determinants and a localized phase on an extended orbital. -/

namespace SourceLaw.PairedFlight

/-- A single scalar interface contributes a rank-one matrix between two
launch and two receiving amplitudes at a fixed frequency. -/
theorem one_contact (a b c d : ℂ) : (a*c)*(b*d) - (a*d)*(b*c) = 0 := by ring

/-- Two equal rails transmit both independent components whenever the
single-rail amplitude is nonzero. -/
theorem paired_contacts {z : ℂ} (h : z ≠ 0) : z*z - 0*0 ≠ 0 := by
  simpa using mul_ne_zero h h

/-- At any untouched occupied address, preserving the selected orbital
identifies its final logical amplitudes with the initial ones. -/
theorem unchanged_site {w a b : ℂ} (hw : w ≠ 0) (h : w*a = w*b) : a = b :=
  mul_left_cancel₀ hw h

/-- An addressed pi phase on one component of an equal two-rail input
moves the displayed probability into other standing modes. -/
theorem local_phase_leakage (w : ℝ) :
    1 - (1-w)^2 - w^2 = 2*w*(1-w) := by ring

end SourceLaw.PairedFlight
