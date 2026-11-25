import ObserverPatchHolography.EinsteinBranch.CollarAdmission
import Mathlib.Tactic.NormNum
import Mathlib.Util.AssertNoSorry

/-! Exact without-replacement normalization after the five-sector station count.
Here b = H/m is the measured greatest station distance per subdivision frequency
and t = 1/R is the inverse outer radius in station hops. The sector count and
the path bounds for b are derived in the station_clock experiment. -/

namespace OPH.EinsteinBranch.CollarAdmission

open Filter
open scoped Topology

/-- The five corner sectors, at offsets zero, one, one, two and two,
counted out to even distance 2k from a single starting triangle. -/
def evenCornerBall (k : ℕ) : ℕ :=
  (k + 1)^2 + 2 * k * (k + 1) + 2 * k^2

theorem corner_annulus (k : ℕ) :
    evenCornerBall (2 * k) = evenCornerBall k + 15*k^2 + 4*k := by
  unfold evenCornerBall
  ring

theorem station_pair_finite (R b : ℝ) (hr : R ≠ 0) (hb : b ≠ 0)
    (hn : 5 * R^4 ≠ 1) :
    let n := (15 / 16 : ℝ)*R^2 + R
    let N := 5*R^4
    let r := 2/(b*R)
    n*(n-1)/(N*(N-1))/r^4 =
      (15/16 + 1/R)*(15/16 + 1/R - (1/R)^2)*b^4 /
        (400*(1 - (1/R)^4/5)) := by
  dsimp
  have hd : (1 - (1 / R)^4 / 5 : ℝ) ≠ 0 := by
    intro h
    apply hn
    field_simp [hr] at h
    linarith
  have hn' : 5 * R^4 - 1 ≠ 0 := sub_ne_zero.mpr hn
  field_simp [hr, hb, hd, hn']
  ring

theorem station_pair_limit (t b : ℕ → ℝ) (ht : Tendsto t atTop (𝓝 0))
    (hb : Tendsto b atTop (𝓝 6)) :
    Tendsto (fun i => (15/16 + t i)*(15/16 + t i - (t i)^2)*(b i)^4 /
      (400*(1 - (t i)^4/5))) atTop (𝓝 (729/256 : ℝ)) := by
  have h := (tendsto_const_nhds (x := (15 / 16 : ℝ))).add ht
  have numerator := (h.mul (h.sub (ht.pow 2))).mul (hb.pow 4)
  have denominator := (tendsto_const_nhds (x := (400 : ℝ))).mul
    ((tendsto_const_nhds (x := (1 : ℝ))).sub ((ht.pow 4).div_const 5))
  convert numerator.div denominator (by norm_num) using 1 <;> norm_num

assert_no_sorry corner_annulus
#print axioms corner_annulus
assert_no_sorry station_pair_finite
assert_no_sorry station_pair_limit
#print axioms station_pair_finite
#print axioms station_pair_limit

end OPH.EinsteinBranch.CollarAdmission
