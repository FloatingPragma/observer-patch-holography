import Mathlib.Data.Real.Basic
import Mathlib.Algebra.Order.Field.Basic
import Mathlib.Tactic.Linarith

/-! A common neutral probe acts as a T₃ + b Y on the calibrated weak frame.
The three readings use the exterior weights of masks 24, 15 and 23. -/
namespace OPH.EndpointNeutralDirection

/-- A resolved singlet and an independently prepared upper doublet member
determine both coefficients and predict the lower doublet reading. -/
theorem reconstruct (a b s u charged : ℝ)
    (hs : s = b) (hu : 2 * u = a - b) (hc : 2 * charged = -a - b) :
    a = s + 2 * u ∧ b = s ∧ charged = -s - u := by
  constructor
  · linarith
  constructor <;> linarith

/-- An upper-member null fixes the direction up to the nonzero singlet scale. -/
theorem electromagnetic_line (a b s : ℝ) (hs : s = b) (hu : a - b = 0) :
    (a, b) = (s, s) := by
  apply Prod.ext <;> dsimp <;> linarith

/-- The two generators agree on the singlet and differ on the upper
doublet member for any resolved singlet response. -/
theorem singlet_ambiguity (s : ℝ) (hs : s ≠ 0) :
    (s - s) / 2 ≠ (0 - s) / 2 := by
  intro h
  apply hs
  linarith

/-- A two-sided interval on the upper member gives the same signed
coefficient bound even when zero is inside the experimental interval. -/
theorem coefficient_error (a b u ε : ℝ) (hu : 2 * u = a - b)
    (lower : -ε ≤ u) (upper : u ≤ ε) :
    -(2 * ε) ≤ a - b ∧ a - b ≤ 2 * ε := by
  constructor <;> linarith

/-- An independent nonzero singlet sensitivity turns the count interval
into a bound on the relative admixture of the two generators. -/
theorem relative_error (a b u ε floor : ℝ) (hu : 2 * u = a - b)
    (lower : -ε ≤ u) (upper : u ≤ ε) (hfloor : 0 < floor) (hs : floor ≤ |b|) :
    |a / b - 1| ≤ 2 * ε / floor := by
  have hb : 0 < |b| := lt_of_lt_of_le hfloor hs
  have hdiff : |a - b| ≤ 2 * ε := abs_le.mpr (coefficient_error a b u ε hu lower upper)
  rw [div_sub_one (abs_pos.mp hb), abs_div]
  calc
    |a - b| / |b| ≤ (2 * ε) / |b| := div_le_div_of_nonneg_right hdiff hb.le
    _ ≤ 2 * ε / floor := div_le_div_of_nonneg_left (by linarith) hfloor hs

end OPH.EndpointNeutralDirection
