import Time.ClockComparison
import Mathlib.Tactic

/-! A clock with two material strokes per complete oscillator period is
sampled on one parity of the stroke count. Two common crests fix its
comparison map. A separate finite counting result expresses the missing
multiplicity when angular cells and a local clock supply three scale powers. -/
namespace OPH.ScreenDensityClock

open OPH.D1Time
noncomputable section

/-- Two counted marker strokes and three counted oscillator flips precede launch. -/
def comparison (phase : ℝ) : PositiveAffineClockMap :=
  ⟨1/2, 2 - phase/2, by norm_num⟩

theorem clock_at_common_crest (phase periods : ℝ) :
    (comparison phase).map (2 + 2*periods + phase) = 3 + periods := by
  simp only [comparison, PositiveAffineClockMap.map]
  ring

/-- The observations at the two anchors identify the map used to predict
any subsequent same-phase crest. -/
theorem anchors_identify {phase i j : ℝ} (h : i < j) :
    comparison phase = PositiveAffineClockMap.throughTwoPoints
      (2+2*i+phase) (2+2*j+phase) (3+i) (3+j) (by linarith) (by linarith) :=
  PositiveAffineClockMap.throughTwoPoints_unique _ _ _ _ _ _ _
    (clock_at_common_crest _ _) (clock_at_common_crest _ _)

theorem missed_period (phase periods : ℝ) :
    (3+periods-1) - (comparison phase).map (2+2*periods+phase) = -1 := by
  rw [clock_at_common_crest]
  ring

/-- For a fixed conformal crest interval, a constant proper-clock cadence
at each end gives the measured signal period ratio. -/
theorem plateau_period_ratio {k qe qr : ℝ} (hk : 0 < k) (he : 0 < qe) :
    (k*qr)/(k*qe)=qr/qe := by field_simp

/-- Angular enrollment gives q² publishers and q strokes per publisher.
At unit reference normalization a further multiplicity is determined by
the fourth-power count criterion. -/
theorem required_multiplicity {q w : ℝ} (hq : 0 < q) :
    (q^2*q)*w=q^4 ↔ w=q := by
  rw [show q^4=(q^2*q)*q by ring]
  exact mul_right_inj' (ne_of_gt (mul_pos (sq_pos_of_pos hq) hq))

theorem angular_census_shortfall {q : ℝ} (hq : 1 < q) : q^2*q < q^4 := by
  have h : 0 < q^2*q := by positivity
  calc
    q^2*q = (q^2*q)*1 := by ring
    _ < (q^2*q)*q := mul_lt_mul_of_pos_left hq h
    _ = q^4 := by ring

#print axioms clock_at_common_crest
#print axioms anchors_identify
#print axioms missed_period
#print axioms plateau_period_ratio
#print axioms required_multiplicity
#print axioms angular_census_shortfall

end
end OPH.ScreenDensityClock
