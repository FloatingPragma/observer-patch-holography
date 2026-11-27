import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith

namespace IntervalEcho.PolarSpectralRead

/-- From an enrolled finite positive spectral measure, the measured acceptance
of a window bounds the mass of its selected pole. All off-line terms enter. -/
theorem residue_lower {ι : Type*} [DecidableEq ι] (s : Finset ι)
    (w f : ι → ℝ) (i : ι) (hi : i ∈ s) (p noise leak : ℝ)
    (hw : ∀ j ∈ s, 0 ≤ w j) (total : ∑ j ∈ s, w j ≤ 1)
    (transmission : f i ≤ 1) (hleak : 0 ≤ leak)
    (outside : ∀ j ∈ s.erase i, f j ≤ leak)
    (meter : p ≤ (∑ j ∈ s, w j * f j) + noise) :
    p - noise - leak ≤ w i := by
  have member : ∀ j ∈ s.erase i, j ∈ s := fun _ h => Finset.mem_of_mem_erase h
  have background : (∑ j ∈ s.erase i, w j * f j) ≤ (∑ j ∈ s.erase i, w j) * leak := by
    rw [Finset.sum_mul]
    exact Finset.sum_le_sum fun j hj => mul_le_mul_of_nonneg_left (outside j hj) (hw j (member j hj))
  have mass := Finset.sum_erase_add s w hi
  have spec := Finset.sum_erase_add s (fun j => w j * f j) hi
  have rest : ∑ j ∈ s.erase i, w j ≤ 1 := by linarith [hw i hi]
  have small := mul_le_mul_of_nonneg_right rest hleak
  have target := mul_le_mul_of_nonneg_left transmission (hw i hi)
  linarith

/-- A measured lower bound on the total pole weight controls its squared
contributions to the long-time return, including coincident frequencies. -/
theorem squared_weight_floor {ι : Type*} (s : Finset ι) (w : ι → ℝ)
    (mass : ℝ) (hn : 0 ≤ mass) (hm : mass ≤ ∑ j ∈ s, w j) :
    mass ^ 2 ≤ (s.card : ℝ) * ∑ j ∈ s, w j ^ 2 := by
  have hsum : 0 ≤ (∑ j ∈ s, w j) + mass := by linarith
  calc
    mass ^ 2 ≤ (∑ j ∈ s, w j) ^ 2 := by
      nlinarith [mul_nonneg (sub_nonneg.mpr hm) hsum]
    _ ≤ _ := by simpa using Finset.sum_mul_sq_le_sq_mul_sq s (fun _ => (1 : ℝ)) w

end IntervalEcho.PolarSpectralRead
