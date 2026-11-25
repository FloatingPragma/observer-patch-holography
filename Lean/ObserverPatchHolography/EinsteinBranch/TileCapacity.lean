import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Algebra.Order.Floor.Ring
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import Mathlib.Util.AssertNoSorry

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

open scoped BigOperators in
/-- On any finite word algebra, normalized fixed-volume transfers annihilate
the constant and volume parts of the entropy-area gradient. Uniform errors
in the completion-count calibration and the reference log weights bound the
entire first variation, with no dependence on the number of words. -/
theorem fixed_volume_capacity_bound {ι : Type*} [Fintype ι]
    (v volume area count prob : ι → ℝ) (Z η c volCoeff ea ep : ℝ)
    (hZ : 0 < Z) (hd : ∀ i, 0 < count i) (hp : ∀ i, 0 < prob i)
    (hnorm : ∑ i, v i = 0) (hvol : ∑ i, v i * volume i = 0)
    (harea : ∀ i, |Real.log (count i) - η * area i - c - volCoeff * volume i| ≤ ea)
    (href : ∀ i, |Real.log (prob i / (count i / Z))| ≤ ep) :
    |∑ i, v i * (-Real.log (prob i) + η * area i)| ≤
      (∑ i, |v i|) * (ea + ep) := by
  let e := fun i => Real.log (count i) - η * area i - c - volCoeff * volume i
  let b := fun i => Real.log (prob i / (count i / Z))
  have decomp (i : ι) :
      v i * (-Real.log (prob i) + η * area i) =
        v i * (Real.log Z - c) - (v i * volume i) * volCoeff + v i * (-(e i + b i)) := by
    dsimp [b, e]
    rw [Real.log_div (hp i).ne' (div_pos (hd i) hZ).ne', Real.log_div (hd i).ne' hZ.ne']
    ring
  have each (i : ι) : |v i * (-(e i + b i))| ≤ |v i| * (ea + ep) := by
    rw [abs_mul, abs_neg]
    exact mul_le_mul_of_nonneg_left ((abs_add_le _ _).trans
      (add_le_add (harea i) (href i))) (abs_nonneg _)
  simp_rw [decomp]
  rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.sum_mul,
      ← Finset.sum_mul, hnorm, hvol]
  simp only [zero_mul, sub_zero, zero_add]
  calc
    _ ≤ ∑ i, |v i * (-(e i + b i))| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, |v i| * (ea + ep) := Finset.sum_le_sum (fun i _ => each i)
    _ = _ := by rw [Finset.sum_mul]

open scoped BigOperators in
/-- On a faithful finite record distribution all signed probability tangents
are feasible locally in the simplex. If all normalized fixed-mean-volume
transfers are admitted, stationarity is equivalent to a three-record test.
For affine surveyed area use score i = η * area i - log (prob i).
Each test uses the same two anchor records, of unequal measured volume. -/
theorem fixed_volume_record_test {ι : Type*} [Fintype ι]
    (volume score : ι → ℝ) (a b : ι) (hab : volume b ≠ volume a) :
    (∀ v : ι → ℝ, (∑ i, v i) = 0 → (∑ i, v i * volume i) = 0 →
      (∑ i, v i * score i) = 0) ↔
    ∀ i, (volume b - volume a) * score i =
      (volume b - volume i) * score a + (volume i - volume a) * score b := by
  classical
  constructor
  · intro h i
    let v : ι → ℝ := fun j =>
      (if j = i then volume b - volume a else 0) -
      (if j = a then volume b - volume i else 0) -
      (if j = b then volume i - volume a else 0)
    have transfer (f : ι → ℝ) :
        (∑ j, v j * f j) = (volume b - volume a) * f i -
          (volume b - volume i) * f a - (volume i - volume a) * f b := by
      simp [v, sub_mul, ite_mul, Finset.sum_sub_distrib]
    have hn := transfer (fun _ => 1)
    have hv := transfer volume
    have hs := transfer score
    simp only [mul_one] at hn
    have hvol : (∑ j, v j * volume j) = 0 := by rw [hv]; ring
    have hz := h v (by linarith only [hn]) hvol
    linarith only [hs, hz]
  · intro h v hn hv
    have he : (volume b - volume a) * (∑ i, v i * score i) =
        (∑ i, v i) * (volume b * score a - volume a * score b) +
        (∑ i, v i * volume i) * (score b - score a) := by
      rw [Finset.mul_sum, Finset.sum_mul, Finset.sum_mul, ← Finset.sum_add_distrib]
      apply Finset.sum_congr rfl
      intro i _
      calc
        _ = v i * ((volume b - volume a) * score i) := by ring
        _ = _ := by rw [h i]; ring
    rw [hn, hv] at he
    exact (mul_eq_zero.mp (by simpa using he)).resolve_left (sub_ne_zero.mpr hab)

assert_no_sorry fixed_volume_record_test
#print axioms fixed_volume_record_test

/-- An aperture doubling as accessible volume quadruples admits no positive
isotropic length ratio preserving both a three-dimensional volume and area. -/
theorem no_isotropic_aperture_rescaling (h : ℝ) (hpos : 0 < h)
    (volume : 4*h^3 = 1) (area : 2*h^2 = 1) : False := by
  have hhalf : 2*h = 1 := by
    nlinarith [congrArg (fun t : ℝ => 2*h*t) area]
  nlinarith [congrArg (fun t : ℝ => h*t) hhalf]

/-- Splitting the middle occupancy equally between its two neighbouring
classes preserves the mean cut load. These two classes have equal word
multiplicities by complementing a subset. Their loss of log multiplicity
pays a strictly positive interior entropy slope at the uniform word reference. -/
theorem balanced_class_spread (r : ℕ) (hr : 0 < r) :
    Real.log (((2*r).choose r : ℕ) : ℝ) -
        (Real.log (((2*r).choose (r-1) : ℕ) : ℝ) +
         Real.log (((2*r).choose (r+1) : ℕ) : ℝ))/2 > 0 := by
  have hs : (2*r).choose (r-1) = (2*r).choose (r+1) :=
    Nat.choose_symm_of_eq_add (by omega)
  have hrel := Nat.choose_succ_right_eq (2*r) r
  rw [show 2*r-r=r by omega] at hrel
  have hp : 0 < (2*r).choose (r+1) := Nat.choose_pos (by omega)
  have hlt : (2*r).choose (r+1) < (2*r).choose r := by
    by_contra! h
    have hle := Nat.mul_le_mul_right r h
    nlinarith
  have hlog := Real.log_lt_log (show 0 < (((2*r).choose (r+1) : ℕ) : ℝ) by exact_mod_cast hp)
    (show (((2*r).choose (r+1) : ℕ) : ℝ) < (((2*r).choose r : ℕ) : ℝ) by exact_mod_cast hlt)
  rw [hs]; linarith

/-- At one host resolution, independent area and volume calibration
constants also cancel from the ratio between two homothetic measurements. -/
theorem two_radii_shape_ratio (a v : ℝ) (hv : v ≠ 0) :
    (a/2)^3/(v/4)^2 = 2*(a^3/v^2) := by
  field_simp
  ring

end OPH.EinsteinBranch.TileCapacity
