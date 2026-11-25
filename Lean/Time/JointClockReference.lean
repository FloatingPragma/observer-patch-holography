import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum

/-!
# When the modular rates of two classical energy labels factor

A two-hand reference diagonal in the product energy basis has positive
weights a,b,c,d in row-major order. On an Alice transition the logarithmic
gap is log(a/c) or log(b/d), according to Bob's energy label. The statement
concerns the continuous rates; finite time samples can alias a phase.
-/
namespace OPH.Time.JointClockReference

theorem gaps_equal_iff {a b c d : ℝ} (ha : 0 < a) (hb : 0 < b)
    (hc : 0 < c) (hd : 0 < d) :
    Real.log (a / c) = Real.log (b / d) ↔ a*d = b*c := by
  have h : Real.log (a/c) = Real.log (b/d) ↔ a/c = b/d :=
    ⟨Real.log_injOn_pos (div_pos ha hc) (div_pos hb hd), congrArg Real.log⟩
  exact h.trans (div_eq_div_iff hc.ne' hd.ne')

def productOfMarginals (a b c d : ℝ) : Prop :=
  a = (a+b)*(a+c) ∧ b = (a+b)*(b+d) ∧
  c = (c+d)*(a+c) ∧ d = (c+d)*(b+d)

theorem determinant_iff_product {a b c d : ℝ} (ht : a+b+c+d=1) :
    a*d=b*c ↔ productOfMarginals a b c d := by
  constructor
  · intro h
    unfold productOfMarginals
    constructor
    · nlinarith [congrArg (a * ·) ht]
    constructor
    · nlinarith [congrArg (b * ·) ht]
    constructor
    · nlinarith [congrArg (c * ·) ht]
    · nlinarith [congrArg (d * ·) ht]
  · rintro ⟨h, _⟩
    nlinarith [congrArg (a * ·) ht]

theorem common_rate_iff_product {a b c d : ℝ} (ha : 0 < a) (hb : 0 < b)
    (hc : 0 < c) (hd : 0 < d) (ht : a+b+c+d=1) :
    Real.log (a/c) = Real.log (b/d) ↔ productOfMarginals a b c d :=
  (gaps_equal_iff ha hb hc hd).trans (determinant_iff_product ht)

theorem exchange_rates_differ :
    Real.log ((5/12 : ℝ)/(1/4)) ≠ Real.log ((1/4 : ℝ)/(1/12)) := by
  intro h
  have := (gaps_equal_iff (a:=(5/12 : ℝ)) (b:=1/4) (c:=1/4) (d:=1/12)
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)).mp h
  norm_num at this

end OPH.Time.JointClockReference
