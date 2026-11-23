import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Algebra.Order.Floor.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity

/-! # First variations from accessible tile multiplicities

A full tile of area x has entropy contrast log(1+x). Physically admissible
interpretations of the tile alphabet and of area are inputs to the comparison.
-/

namespace OPH.EinsteinBranch.TileCapacity

/-- Losing the partial tile at each of two endpoints changes their difference
by at most one tile. -/
theorem floor_contrast (a b : ℝ) :
    |(a - b) - ((⌊a⌋ : ℝ) - (⌊b⌋ : ℝ))| ≤ 1 := by
  have ha := Int.floor_le a
  have hb := Int.floor_le b
  have ha' := Int.lt_floor_add_one a
  have hb' := Int.lt_floor_add_one b
  apply abs_le.mpr
  constructor <;> linarith

/-- The area-normalized log multiplicity differs from one by at most x. -/
theorem normalized_capacity (x : ℝ) (hx : 0 < x) :
    0 ≤ 1 - Real.log (1+x)/x ∧ 1 - Real.log (1+x)/x ≤ x := by
  have hpos : 0 < 1+x := by linarith
  have upper := Real.log_le_sub_one_of_pos hpos
  have lower := Real.one_sub_inv_le_log_of_pos hpos
  have hxlow : x / (1+x) ≤ Real.log (1+x) := by
    calc
      x / (1+x) = 1 - (1+x)⁻¹ := by field_simp; ring
      _ ≤ _ := lower
  have hx2 : x*(1-x) ≤ x / (1+x) := by
    apply (le_div_iff₀ hpos).2
    nlinarith [mul_nonneg hx.le (sq_nonneg x)]
  constructor
  · have h := (div_le_iff₀ hx).2 (show Real.log (1+x) ≤ 1*x by linarith)
    linarith
  · have h := (le_div_iff₀ hx).2 (show (1-x)*x ≤ Real.log (1+x) by linarith [hx2.trans hxlow])
    linarith

/-- Uniform first-variation bound for a pair of areas, including the partial
tiles and the entropy of the unequal accessible alphabets. -/
theorem transfer_residual (x a b : ℝ) (hx : 0 < x) :
    |(a-b) - ((⌊a/x⌋ : ℝ) - (⌊b/x⌋ : ℝ)) * Real.log (1+x)| ≤
      (|a-b| + 1) * x := by
  let L := Real.log (1+x)
  let n : ℝ := (⌊a/x⌋ : ℝ) - (⌊b/x⌋ : ℝ)
  have hL : 0 ≤ L := Real.log_nonneg (by linarith)
  have hLx : L ≤ x := by
    have h := Real.log_le_sub_one_of_pos (show 0 < 1+x by linarith)
    dsimp [L]; linarith
  have hround : |(a-b)/x - n| ≤ 1 := by
    simpa only [sub_div, n] using floor_contrast (a/x) (b/x)
  obtain ⟨hlow,hhigh⟩ := normalized_capacity x hx
  calc
    |(a-b) - n * L| =
        |(a-b) * (1-L/x) + ((a-b)/x-n) * L| := by congr 1; ring
    _ ≤ |(a-b) * (1-L/x)| + |((a-b)/x-n) * L| := abs_add_le _ _
    _ = |a-b| * (1-L/x) + |(a-b)/x-n| * L := by
      rw [abs_mul, abs_mul, abs_of_nonneg hlow, abs_of_nonneg hL]
    _ ≤ |a-b| * x + 1*x :=
      add_le_add (mul_le_mul_of_nonneg_left hhigh (abs_nonneg _))
        (mul_le_mul hround hLx hL (by norm_num))
    _ = (|a-b| + 1) * x := by ring

/-- For a shared scalar bulk readout the spread in held area slopes directly
lower-bounds the sum of the stationarity errors in the two frames. -/
theorem shared_bulk_frame_spread (bulk a b : ℝ) :
    |a-b| ≤ |bulk+a| + |bulk+b| := by
  calc
    |a-b| = |(bulk+a)-(bulk+b)| := by congr 1; ring
    _ ≤ _ := abs_sub _ _

/-- At a uniform word equilibrium, the interior probability of each word is
its positive exterior completion count divided by the total count. The
coarse Shannon derivative on a unit transfer cancels the log-capacity slope. -/
theorem uniform_word_transfer (m n z : ℝ) (hm : 0 < m) (hn : 0 < n) (hz : 0 < z) :
    Real.log (m/z) - Real.log (n/z) + (Real.log n - Real.log m) = 0 := by
  rw [Real.log_div hm.ne' hz.ne', Real.log_div hn.ne' hz.ne']
  ring

/-- An aperture doubling as accessible volume quadruples admits no positive
isotropic length ratio preserving both a three-dimensional volume and area. -/
theorem no_isotropic_aperture_rescaling (h : ℝ) (hpos : 0 < h)
    (volume : 4*h^3 = 1) (area : 2*h^2 = 1) : False := by
  have hhalf : 2*h = 1 := by
    nlinarith [congrArg (fun t : ℝ => 2*h*t) area]
  nlinarith [congrArg (fun t : ℝ => h*t) hhalf]

end OPH.EinsteinBranch.TileCapacity
