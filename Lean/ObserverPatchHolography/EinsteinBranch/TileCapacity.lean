import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
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

namespace PartitionRepair

open scoped BigOperators

/-- A fixed graph region chooses the visible endpoints of each bond and includes
or excludes each cell's volume bit. Slots on distinct bonds or cells are independent
only in the reference, not in the varied preparation. -/
inductive Slot where
  | hiddenBond | bothEnds | leftEnd | rightEnd | hiddenCell | volumeBit
  deriving DecidableEq

abbrev tickets : Slot → ℕ
  | .hiddenCell | .volumeBit => 3
  | _ => 5

abbrev letters : Slot → ℕ
  | .hiddenBond | .hiddenCell => 1
  | .bothEnds => 5
  | .leftEnd | .rightEnd => 3
  | .volumeBit => 2

/-- The five ordered legal bond words, in the repair cycle's own order. -/
def endpoints : Fin 5 → Fin 3 × Fin 3
  | 0 => (0, 0)
  | 1 => (1, 1)
  | 2 => (1, 2)
  | 3 => (2, 2)
  | _ => (2, 1)

/-- One of three equiprobable cell tickets is occupied, giving the repair's
two-to-one stationary odds. The bond ticket is exactly the legal word. -/
def observe : (s : Slot) → Fin (tickets s) → Fin (letters s)
  | .hiddenBond, _ | .hiddenCell, _ => 0
  | .bothEnds, t => t
  | .leftEnd, t => (endpoints t).1
  | .rightEnd, t => (endpoints t).2
  | .volumeBit, t => if t = 2 then 1 else 0

def aperture : (s : Slot) → Fin (letters s) → ℕ
  | .leftEnd, x | .rightEnd, x => if x = 0 then 0 else 1
  | _, _ => 0

def occupancy : (s : Slot) → Fin (letters s) → ℕ
  | .volumeBit, x => x.val
  | _, _ => 0

def baseMultiplicity : Slot → ℕ
  | .hiddenBond => 5
  | .hiddenCell => 3
  | .volumeBit => 2
  | _ => 1

def count (s : Slot) (a : Fin (letters s)) : ℕ :=
  ∑ t : Fin (tickets s), if observe s t = a then 1 else 0

/-- The fibers of the two actual endpoint projections have sizes 1,2,2.
Multiplication by two per occupied volume bit clears its reciprocal weight. -/
theorem fiber_identity (s : Slot) (a : Fin (letters s)) :
    0 < count s a ∧
    count s a * 2 ^ occupancy s a = baseMultiplicity s * 2 ^ aperture s a := by
  cases s <;> revert a <;> decide

noncomputable def weight (s : Slot) (a : Fin (letters s)) : ℝ :=
  count s a / (tickets s : ℝ)

/-- The microscopic commands use the same two fair bits as the executable
repair. The two middle commands hold the cyclic bond word. -/
def advanceBond (i : Fin 5) (command : Fin 4) : Fin 5 :=
  if command = 0 then ⟨(i.val + 4) % 5, Nat.mod_lt _ (by decide)⟩
  else if command = 3 then ⟨(i.val + 1) % 5, Nat.mod_lt _ (by decide)⟩ else i

def advanceVolume (v : Fin 2) (command : Fin 4) : Fin 2 :=
  if v = 0 then (if command = 3 then 1 else 0)
  else if 2 ≤ command.val then 1 else 0

/-- Uniform legal bond words and the counted volume marginal are invariant
under the actual four commands. Every bond communicates with every other
in two steps and every volume bit with every other in one. -/
theorem repair_selects_weights :
    (∀ j : Fin 5, (∑ i : Fin 5, ∑ c : Fin 4, if advanceBond i c = j then 1 else 0) = 4) ∧
    (∀ b : Fin 2, (∑ v : Fin 2, count .volumeBit v *
      (∑ c : Fin 4, if advanceVolume v c = b then 1 else 0)) = 4 * count .volumeBit b) ∧
    (∀ i j : Fin 5, ∃ a b : Fin 4, advanceBond (advanceBond i a) b = j) ∧
    (∀ i j : Fin 2, ∃ c : Fin 4, advanceVolume i c = j) := by decide

noncomputable def baseLog (s : Slot) : ℝ :=
  Real.log (baseMultiplicity s) - Real.log (tickets s)

private theorem weight_pos (s : Slot) (a : Fin (letters s)) : 0 < weight s a := by
  have ht : 0 < tickets s := by cases s <;> decide
  exact div_pos (by exact_mod_cast (fiber_identity s a).1) (by exact_mod_cast ht)

private theorem log_weight (s : Slot) (a : Fin (letters s)) :
    Real.log (weight s a) = baseLog s + Real.log 2 *
      ((aperture s a : ℝ) - (occupancy s a : ℝ)) := by
  obtain ⟨hc, he⟩ := fiber_identity s a
  have ht : 0 < tickets s := by cases s <;> decide
  have hb : 0 < baseMultiplicity s := by cases s <;> decide
  have h : (count s a : ℝ) * (2 : ℝ) ^ occupancy s a =
      (baseMultiplicity s : ℝ) * (2 : ℝ) ^ aperture s a := by exact_mod_cast he
  have hl := congrArg Real.log h
  rw [Real.log_mul (by positivity) (by positivity),
      Real.log_mul (by positivity) (by positivity), Real.log_pow, Real.log_pow] at hl
  unfold weight baseLog
  rw [Real.log_div (by positivity) (by positivity)]
  linarith

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

abbrev RegionalWord (slots : ι → Slot) := (i : ι) → Fin (letters (slots i))

noncomputable def reference (slots : ι → Slot) (w : RegionalWord slots) : ℝ :=
  ∏ i, weight (slots i) (w i)

/-- The regional product counts the mass of every whole global ticket word
whose restriction is the requested regional record. In particular a different
cut uses the same microscopic word distribution. -/
theorem reference_from_tickets (slots : ι → Slot) (w : RegionalWord slots) :
    reference slots w = ∑ t : (i : ι) → Fin (tickets (slots i)),
      ∏ i, (if observe (slots i) (t i) = w i then 1 else 0 : ℝ) / tickets (slots i) := by
  rw [← Fintype.prod_sum (fun i (t : Fin (tickets (slots i))) =>
    (if observe (slots i) t = w i then 1 else 0 : ℝ) / tickets (slots i))]
  apply Finset.prod_congr rfl
  intro i _
  simp only [weight, count, div_eq_mul_inv, ← Finset.sum_mul]
  congr 1
  push_cast
  rfl

theorem reference_normalized (slots : ι → Slot) : ∑ w, reference slots w = 1 := by
  unfold reference
  rw [← Fintype.prod_sum]
  have h (s : Slot) : ∑ a, weight s a = 1 := by
    have hn : ∑ a, count s a = tickets s := by cases s <;> decide
    have ht : tickets s ≠ 0 := by cases s <;> decide
    simp [weight, div_eq_mul_inv, ← Finset.sum_mul, ← Nat.cast_sum, hn, ht]
  simp [h]

noncomputable def totalArea (slots : ι → Slot) (w : RegionalWord slots) : ℝ :=
  ∑ i, (aperture (slots i) (w i) : ℝ)

noncomputable def totalVolume (slots : ι → Slot) (w : RegionalWord slots) : ℝ :=
  ∑ i, (occupancy (slots i) (w i) : ℝ)

omit [DecidableEq ι] in
theorem regional_log_law (slots : ι → Slot) (w : RegionalWord slots) :
    Real.log (reference slots w) = (∑ i, baseLog (slots i)) +
      Real.log 2 * (totalArea slots w - totalVolume slots w) := by
  unfold reference totalArea totalVolume
  rw [Real.log_prod (fun i _ => (weight_pos _ (w i)).ne')]
  simp_rw [log_weight]
  rw [Finset.sum_add_distrib, ← Finset.mul_sum, Finset.sum_sub_distrib]

/-- Differentiate the actual regional Shannon entropy. The words along the
path can be arbitrarily correlated. Normalization of the tangent and its
point of departure are the only restrictions. This gives the zero derivative
of entropy plus the calibrated area at fixed mean volume. -/
theorem shannon_first_variation (slots : ι → Slot)
    (p : RegionalWord slots → ℝ → ℝ) (v : RegionalWord slots → ℝ)
    (at_reference : ∀ w, p w 0 = reference slots w)
    (along : ∀ w, HasDerivAt (p w) (v w) 0) (normalized : ∑ w, v w = 0) :
    HasDerivAt (fun t => ∑ w, Real.negMulLog (p w t))
      (Real.log 2 * ((∑ w, v w * totalVolume slots w) -
        (∑ w, v w * totalArea slots w))) 0 := by
  have pos (w : RegionalWord slots) : 0 < p w 0 := by
    rw [at_reference]
    exact Finset.prod_pos fun i _ => weight_pos _ (w i)
  have hs := HasDerivAt.fun_sum (u := Finset.univ) (fun w _ =>
    (Real.hasDerivAt_negMulLog (pos w).ne').comp 0 (along w))
  have term (w : RegionalWord slots) :
      (-Real.log (p w 0) - 1) * v w =
        v w * (-(∑ i, baseLog (slots i)) - 1) +
        Real.log 2 * (v w * totalVolume slots w - v w * totalArea slots w) := by
    rw [at_reference, regional_log_law]
    ring
  simpa only [Function.comp_def, term, Finset.sum_add_distrib, ← Finset.sum_mul,
    ← Finset.mul_sum, normalized, zero_mul, zero_add, Finset.sum_sub_distrib] using hs

/-- At fixed first variation of the region's independently counted volume, the
differentiable sum of the Shannon entropy and the mean calibrated aperture
has exactly zero derivative. -/
theorem fixed_volume_stationarity (slots : ι → Slot)
    (p : RegionalWord slots → ℝ → ℝ) (v : RegionalWord slots → ℝ)
    (at_reference : ∀ w, p w 0 = reference slots w)
    (along : ∀ w, HasDerivAt (p w) (v w) 0) (normalized : ∑ w, v w = 0)
    (fixed_volume : ∑ w, v w * totalVolume slots w = 0) :
    HasDerivAt (fun t => (∑ w, Real.negMulLog (p w t)) +
      Real.log 2 * ∑ w, p w t * totalArea slots w) 0 0 := by
  have hS := shannon_first_variation slots p v at_reference along normalized
  have hA := HasDerivAt.fun_sum (u := Finset.univ) (fun w _ =>
    (along w).mul_const (totalArea slots w))
  simpa only [fixed_volume, zero_sub, mul_neg, neg_add_cancel] using
    hS.add (hA.const_mul (Real.log 2))

assert_no_sorry repair_selects_weights
assert_no_sorry fixed_volume_stationarity
#print axioms fixed_volume_stationarity
assert_no_sorry reference_from_tickets
assert_no_sorry shannon_first_variation
#print axioms shannon_first_variation

end PartitionRepair

end OPH.EinsteinBranch.TileCapacity
