import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith

/-! Probability accounting for a uniformly addressed three-site spin transfer.
The source allocation bounds the occupation and electric readings together.
All weights refer to one original sign bank, including its refused starts. -/
namespace OPH.EndpointArrival

/-- The two-link and three-site refusal allocation preserves this much of the
original preparation. The electric contribution has additional slack. -/
theorem admission (occupation casimir deficit accepted : ℝ)
    (hc : 0 ≤ casimir)
    (budget : occupation / 16 + 5 * casimir / 18 ≤ deficit)
    (refusal : 1 - accepted ≤ 3 * occupation / 16 + 5 * casimir / 9) :
    1 - 3 * deficit ≤ accepted ∧
      1 - 3 * deficit + 5 * casimir / 18 ≤ accepted := by
  constructor <;> linarith

/-- A centered second moment bounds the discarded tails. A finite reader
can additionally lose the specified fraction of all original starts. -/
theorem recorded_band (deficit accepted tail reported : ℝ)
    (ha : 1 - 3 * deficit ≤ accepted)
    (ht : 900 * tail ≤ 248) (hr : accepted - tail - 1 / 32 ≤ reported) :
    1 - 3 * deficit - 248 / 900 - 1 / 32 ≤ reported := by
  linarith

end OPH.EndpointArrival
