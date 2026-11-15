import Geometry.SourcePopulationQuadrature
import Mathlib.NumberTheory.Real.GoldenRatio
import Mathlib.Data.Nat.ModEq
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import Mathlib.MeasureTheory.Function.LocallyIntegrable
import Mathlib.Analysis.SpecificLimits.Basic

set_option autoImplicit false

/-!
# Golden-source population assignment

The points are the actual fractional parts `b φ - floor (b φ)` used by
the conservative source-record construction. Fibonacci arithmetic proves
their proximity to a permuted equal grid; the grid does not replace the
population. Cell assignment, unlike nearest-point covering, may assign a
point outside its cell. No source-production, read-law or clock is derived.
-/

namespace OPH.GoldenSourceAssignment

open scoped BigOperators NNReal ENNReal

/-- Exact one-dimensional source readback before the length scale. -/
noncomputable def orbit (b : ℕ) : ℝ := Int.fract ((b : ℝ) * Real.goldenRatio)

/-- The equal-grid cell index, with the same integer residue as the paper. -/
def residue (p q b : ℕ) : ℕ := (b * p) % q

/-- Multiplication by the next Fibonacci number permutes all cell indices. -/
theorem residue_injective {p q : ℕ} (hc : q.Coprime p) :
    Function.Injective (fun b : Fin q => residue p q b.val) := by
  intro b d h
  have hm : b.val ≡ d.val [MOD q] := Nat.ModEq.cancel_right_of_coprime hc h
  apply Fin.ext
  simpa only [Nat.ModEq, Nat.mod_eq_of_lt b.isLt, Nat.mod_eq_of_lt d.isLt] using hm

/-- Every equal-grid cell is assigned once, rather than merely being near some site. -/
theorem residue_bijective {p q : ℕ} (hq : 0 < q) (hc : q.Coprime p) :
    Function.Bijective (fun b : Fin q =>
      (⟨residue p q b.val, Nat.mod_lt _ hq⟩ : Fin q)) := by
  apply Function.Injective.bijective_of_finite
  intro b d h
  exact residue_injective hc (congrArg Fin.val h)

theorem residue_pos {p q b : ℕ} (hc : q.Coprime p)
    (hb : 0 < b) (hbq : b < q) : 0 < residue p q b := by
  by_contra h
  have hz : residue p q b = 0 := Nat.eq_zero_of_not_pos h
  have hm : b ≡ 0 [MOD q] := Nat.ModEq.cancel_right_of_coprime hc (by
    simpa [Nat.ModEq, residue] using hz)
  have : b = 0 := by simpa only [Nat.ModEq, Nat.mod_eq_of_lt hbq, Nat.zero_mod] using hm
  omega

/-- The approximation error times the denominator is below one.
This derives the estimate from Fibonacci identities, not a supplied Diophantine bound. -/
theorem fibonacci_error {n : ℕ} (hn : 0 < n) :
    |(Nat.fib (n + 1) : ℝ) - Real.goldenRatio * Nat.fib n| <
      1 / (Nat.fib n : ℝ) := by
  have hq : (0 : ℝ) < Nat.fib n := by exact_mod_cast Nat.fib_pos.mpr hn
  have hf : (Nat.fib n : ℝ) < Real.goldenRatio ^ n := by
    have hid := Real.fib_succ_sub_goldenConj_mul_fib n
    have hmono : (Nat.fib n : ℝ) ≤ Nat.fib (n + 1) := by
      exact_mod_cast (Nat.fib_le_fib_succ (n := n))
    have hneg := mul_neg_of_neg_of_pos Real.goldenConj_neg hq
    linarith
  rw [Real.fib_succ_sub_goldenRatio_mul_fib, abs_pow,
    abs_of_neg Real.goldenConj_neg, ← Real.inv_goldenRatio, inv_pow, one_div]
  exact inv_strictAnti₀ hq hf

/-- Every orbit point before the denominator moves by less than one grid spacing. -/
theorem fibonacci_scaled_error {n b : ℕ} (hn : 0 < n) (hb : b < Nat.fib n) :
    |(b : ℝ) * Real.goldenRatio - (b : ℝ) * Nat.fib (n + 1) / Nat.fib n| <
      1 / (Nat.fib n : ℝ) := by
  have hq : (0 : ℝ) < Nat.fib n := by exact_mod_cast Nat.fib_pos.mpr hn
  have hb' : (b : ℝ) < Nat.fib n := by exact_mod_cast hb
  have he := fibonacci_error hn
  have herr : |(b : ℝ) * Real.goldenRatio - (b : ℝ) * Nat.fib (n + 1) / Nat.fib n| =
      ((b : ℝ) / Nat.fib n) * |(Nat.fib (n + 1) : ℝ) - Real.goldenRatio * Nat.fib n| := by
    rw [← abs_of_nonneg (div_nonneg (Nat.cast_nonneg b) hq.le), ← abs_mul]
    apply abs_eq_abs.mpr
    right
    field_simp
    ring
  rw [herr]
  have hfactor : (b : ℝ) / Nat.fib n < 1 := (div_lt_one hq).mpr hb'
  calc
    _ ≤ 1 * |(Nat.fib (n + 1) : ℝ) - Real.goldenRatio * Nat.fib n| :=
      mul_le_mul_of_nonneg_right hfactor.le (abs_nonneg _)
    _ < _ := by simpa using he

/-- Integer quotient plus fractional residue, with real division on the left. -/
theorem quotient_residue {p q b : ℕ} (hq : 0 < q) :
    (b : ℝ) * p / q = ((b * p / q : ℕ) : ℝ) + (residue p q b : ℝ) / q := by
  have hid : ((b * p / q : ℕ) : ℝ) * q + (residue p q b : ℝ) = (b : ℝ) * p := by
    exact_mod_cast (show b * p / q * q + residue p q b = b * p by
      simpa only [residue, Nat.mul_comm q] using Nat.div_add_mod (b * p) q)
  have hq' : (q : ℝ) ≠ 0 := by exact_mod_cast hq.ne'
  field_simp
  linarith

/-- No wrapping occurs: the rational and irrational readbacks have the same floor. -/
theorem golden_floor {n b : ℕ} (hn : 0 < n) (hb : b < Nat.fib n) :
    ⌊(b : ℝ) * Real.goldenRatio⌋ = (b * Nat.fib (n + 1) / Nat.fib n : ℕ) := by
  by_cases hb0 : b = 0
  · simp [hb0]
  have hq : 0 < Nat.fib n := Nat.fib_pos.mpr hn
  have hq' : (0 : ℝ) < Nat.fib n := by exact_mod_cast hq
  have hk : 1 ≤ residue (Nat.fib (n + 1)) (Nat.fib n) b :=
    residue_pos (Nat.fib_coprime_fib_succ n) (Nat.pos_of_ne_zero hb0) hb
  have hk' : (1 : ℝ) ≤ residue (Nat.fib (n + 1)) (Nat.fib n) b := by exact_mod_cast hk
  have hkq : (residue (Nat.fib (n + 1)) (Nat.fib n) b : ℝ) + 1 ≤ Nat.fib n := by
    exact_mod_cast (Nat.mod_lt (b * Nat.fib (n + 1)) hq)
  have he := abs_lt.mp (fibonacci_scaled_error hn hb)
  rw [quotient_residue hq] at he
  apply Int.floor_eq_iff.mpr
  simp only [Int.cast_natCast]
  have hlo := (div_le_div_iff_of_pos_right hq').mpr hk'
  have hhi : (residue (Nat.fib (n + 1)) (Nat.fib n) b : ℝ) / Nat.fib n +
      1 / Nat.fib n ≤ 1 := by
    rw [← add_div, div_le_one hq']
    exact hkq
  constructor <;> linarith

/-- Exact affine offset from the permuted grid: the floor term has been eliminated. -/
theorem golden_orbit_eq {n b : ℕ} (hn : 0 < n) (hb : b < Nat.fib n) :
    orbit b = (residue (Nat.fib (n + 1)) (Nat.fib n) b : ℝ) / Nat.fib n +
      ((b : ℝ) * Real.goldenRatio - (b : ℝ) * Nat.fib (n + 1) / Nat.fib n) := by
  unfold orbit Int.fract
  rw [golden_floor hn hb]
  have hid := quotient_residue (p := Nat.fib (n + 1)) (b := b) (Nat.fib_pos.mpr hn)
  simp only [Int.cast_natCast]
  linarith

theorem golden_grid_error {n b : ℕ} (hn : 0 < n) (hb : b < Nat.fib n) :
    |orbit b - (residue (Nat.fib (n + 1)) (Nat.fib n) b : ℝ) / Nat.fib n| <
      1 / (Nat.fib n : ℝ) := by
  rw [golden_orbit_eq hn hb, add_sub_cancel_left]
  exact fibonacci_scaled_error hn hb

/-- Irrationality prevents duplicate sites, independently of their cell labels. -/
theorem orbit_injective : Function.Injective orbit := by
  intro b d h
  by_contra hbd
  have hdiff : (b : ℤ) - d ≠ 0 := by
    rw [sub_ne_zero]
    exact_mod_cast hbd
  have hirr := (Real.goldenRatio_irrational.intCast_mul hdiff).ne_int
    (⌊(b : ℝ) * Real.goldenRatio⌋ - ⌊(d : ℝ) * Real.goldenRatio⌋)
  apply hirr
  unfold orbit Int.fract at h
  push_cast
  nlinarith

/-- The source sample's displacement throughout its assigned interval.
The factor two accounts for both orbit-to-grid and within-cell displacement. -/
theorem coordinate_assignment {n b : ℕ} (hn : 0 < n) (hb : b < Nat.fib n)
    {L x : ℝ} (hL : 0 ≤ L)
    (hlo : L * (residue (Nat.fib (n + 1)) (Nat.fib n) b : ℝ) / Nat.fib n ≤ x)
    (hhi : x ≤ L * ((residue (Nat.fib (n + 1)) (Nat.fib n) b : ℝ) + 1) /
      Nat.fib n) :
    |L * orbit b - x| ≤ 2 * L / Nat.fib n := by
  have he := abs_le.mp (golden_grid_error hn hb).le
  have hlo' := mul_le_mul_of_nonneg_left he.1 hL
  have hhi' := mul_le_mul_of_nonneg_left he.2 hL
  apply abs_le.mpr
  constructor <;> (ring_nf at hlo' hhi' hlo hhi ⊢; linarith)

/-- The three-dimensional source sites in the Euclidean repair carrier. -/
noncomputable def site {n : ℕ} (L : ℝ) (b : Fin 3 → Fin (Nat.fib n)) :
    EuclideanSpace ℝ (Fin 3) := WithLp.toLp 2 (fun i => L * orbit (b i).val)

/-- The pointwise assignment estimate used by the equal-cell count law. -/
theorem tensor_assignment {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 ≤ L)
    (b : Fin 3 → Fin (Nat.fib n)) (x : EuclideanSpace ℝ (Fin 3))
    (hx : ∀ i, L * (residue (Nat.fib (n + 1)) (Nat.fib n) (b i).val : ℝ) /
      Nat.fib n ≤ x i ∧
      x i ≤ L * ((residue (Nat.fib (n + 1)) (Nat.fib n) (b i).val : ℝ) + 1) /
        Nat.fib n) :
    dist (site L b) x ≤ 2 * Real.sqrt 3 * L / Nat.fib n := by
  have hq : (0 : ℝ) < Nat.fib n := by exact_mod_cast Nat.fib_pos.mpr hn
  have hcoord (i : Fin 3) : |(site L b) i - x i| ≤ 2 * L / Nat.fib n :=
    coordinate_assignment hn (b i).isLt hL (hx i).1 (hx i).2
  rw [dist_eq_norm]
  have hsum : ‖site L b - x‖ ^ 2 ≤ 3 * (2 * L / Nat.fib n) ^ 2 := by
    rw [EuclideanSpace.real_norm_sq_eq]
    calc
      _ ≤ ∑ _ : Fin 3, (2 * L / Nat.fib n) ^ 2 := by
        apply Finset.sum_le_sum
        intro i _
        simpa only [PiLp.sub_apply, sq_abs] using
          pow_le_pow_left₀ (abs_nonneg ((site L b) i - x i)) (hcoord i) 2
      _ = _ := by simp
  have hsqrt : (Real.sqrt 3) ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  have hright : 0 ≤ 2 * Real.sqrt 3 * L / Nat.fib n := by positivity
  have hrightsq : (2 * Real.sqrt 3 * L / Nat.fib n) ^ 2 =
      3 * (2 * L / Nat.fib n) ^ 2 := by
    rw [div_pow, mul_pow, mul_pow, hsqrt]
    ring
  nlinarith [hrightsq]

/-- The population has no duplicate Euclidean positions at a positive scale. -/
theorem site_injective {n : ℕ} {L : ℝ} (hL : 0 < L) :
    Function.Injective (site (n := n) L) := by
  intro b d h
  funext i
  apply Fin.ext
  apply orbit_injective
  have hi := congrArg (fun x : EuclideanSpace ℝ (Fin 3) => x i) h
  exact (mul_left_cancel₀ hL.ne') hi

/-- Exact cardinality of the source-position image. -/
theorem site_card {n : ℕ} {L : ℝ} (hL : 0 < L) :
    Fintype.card (Set.range (site (n := n) L)) = (Nat.fib n) ^ 3 := by
  classical
  rw [← Fintype.card_congr (Equiv.ofInjective _ (site_injective hL))]
  simp

open Set MeasureTheory

/-- Half-open permutation cells give a literal partition; omitted outer faces are null. -/
noncomputable def axisCell (n : ℕ) (L : ℝ) (b : Fin (Nat.fib n)) : Set ℝ :=
  Ico (L * (residue (Nat.fib (n + 1)) (Nat.fib n) b.val : ℝ) / Nat.fib n)
    (L * ((residue (Nat.fib (n + 1)) (Nat.fib n) b.val : ℝ) + 1) / Nat.fib n)

theorem axisCell_disjoint {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 < L) :
    Pairwise (fun b d : Fin (Nat.fib n) => Disjoint (axisCell n L b) (axisCell n L d)) := by
  intro b d hbd
  have hk : residue (Nat.fib (n + 1)) (Nat.fib n) b.val ≠
      residue (Nat.fib (n + 1)) (Nat.fib n) d.val :=
    fun h => hbd (residue_injective (Nat.fib_coprime_fib_succ n) h)
  have hq : (0 : ℝ) < Nat.fib n := by exact_mod_cast Nat.fib_pos.mpr hn
  apply Set.disjoint_left.mpr
  intro x hb hd
  rcases lt_or_gt_of_ne hk with h | h
  · have h' : (residue (Nat.fib (n + 1)) (Nat.fib n) b.val : ℝ) + 1 ≤
        residue (Nat.fib (n + 1)) (Nat.fib n) d.val := by exact_mod_cast h
    have := div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left h' hL.le) hq.le
    exact (hb.2.trans_le (this.trans hd.1)).false
  · have h' : (residue (Nat.fib (n + 1)) (Nat.fib n) d.val : ℝ) + 1 ≤
        residue (Nat.fib (n + 1)) (Nat.fib n) b.val := by exact_mod_cast h
    have := div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left h' hL.le) hq.le
    exact (hd.2.trans_le (this.trans hb.1)).false

/-- The cells cover the half-open spatial axis exactly, including its origin. -/
theorem axisCell_cover {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 < L) :
    ⋃ b : Fin (Nat.fib n), axisCell n L b = Ico 0 L := by
  have hq : 0 < Nat.fib n := Nat.fib_pos.mpr hn
  have hq' : (0 : ℝ) < Nat.fib n := by exact_mod_cast hq
  ext x
  constructor
  · rintro ⟨_, ⟨b, rfl⟩, hx⟩
    have hk : (residue (Nat.fib (n + 1)) (Nat.fib n) b.val : ℝ) + 1 ≤ Nat.fib n := by
      exact_mod_cast Nat.mod_lt (b.val * Nat.fib (n + 1)) hq
    refine ⟨(by positivity : 0 ≤ L *
      (residue (Nat.fib (n + 1)) (Nat.fib n) b.val : ℝ) / Nat.fib n).trans hx.1, ?_⟩
    have ht : L * ((residue (Nat.fib (n + 1)) (Nat.fib n) b.val : ℝ) + 1) /
        Nat.fib n ≤ L := (div_le_iff₀ hq').mpr (mul_le_mul_of_nonneg_left hk hL.le)
    exact hx.2.trans_le ht
  · intro hx
    let k : ℕ := ⌊(Nat.fib n : ℝ) * x / L⌋₊
    have hz : 0 ≤ (Nat.fib n : ℝ) * x / L := div_nonneg (mul_nonneg hq'.le hx.1) hL.le
    have hk : k < Nat.fib n := (Nat.floor_lt hz).mpr (by
      apply (div_lt_iff₀ hL).mpr
      nlinarith [hx.2])
    obtain ⟨b, hb⟩ := (residue_bijective hq (Nat.fib_coprime_fib_succ n)).2 ⟨k, hk⟩
    have hbk : residue (Nat.fib (n + 1)) (Nat.fib n) b.val = k := congrArg Fin.val hb
    apply mem_iUnion.mpr
    refine ⟨b, ?_⟩
    change L * _ / _ ≤ x ∧ x < L * (_ + 1) / _
    rw [hbk]
    have hlo := Nat.floor_le hz
    have hhi := Nat.lt_floor_add_one ((Nat.fib n : ℝ) * x / L)
    change (k : ℝ) ≤ _ at hlo
    change _ < (k : ℝ) + 1 at hhi
    constructor
    · apply (div_le_iff₀ hq').mpr
      have := (le_div_iff₀ hL).mp hlo
      nlinarith
    · apply (lt_div_iff₀ hq').mpr
      have := (div_lt_iff₀ hL).mp hhi
      nlinarith

/-- Rectangular cells in ordinary Lebesgue coordinates. -/
noncomputable def cell (n : ℕ) (L : ℝ) (b : Fin 3 → Fin (Nat.fib n)) : Set (Fin 3 → ℝ) :=
  Set.pi univ (fun i => axisCell n L (b i))

noncomputable def cube (L : ℝ) : Set (Fin 3 → ℝ) := Set.pi univ (fun _ => Ico 0 L)

theorem cell_measurable (n : ℕ) (L : ℝ) (b : Fin 3 → Fin (Nat.fib n)) :
    MeasurableSet (cell n L b) := by
  exact MeasurableSet.univ_pi (fun _ => measurableSet_Ico)

theorem cell_disjoint {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 < L) :
    Pairwise (fun b d : Fin 3 → Fin (Nat.fib n) => Disjoint (cell n L b) (cell n L d)) := by
  intro b d hbd
  obtain ⟨i, hi⟩ := Function.ne_iff.mp hbd
  apply Set.disjoint_left.mpr
  intro x hb hd
  exact Set.disjoint_left.mp (axisCell_disjoint hn hL hi) (hb i (mem_univ i)) (hd i (mem_univ i))

theorem cell_cover {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 < L) :
    ⋃ b : Fin 3 → Fin (Nat.fib n), cell n L b = cube L := by
  unfold cell cube
  rw [Set.iUnion_univ_pi]
  simp only [axisCell_cover hn hL]

/-- The exact Lebesgue mass is derived from the rectangular partition. -/
theorem cell_mass {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 ≤ L)
    (b : Fin 3 → Fin (Nat.fib n)) :
    volume.real (cell n L b) = L ^ 3 / (Nat.fib n : ℝ) ^ 3 := by
  unfold cell axisCell Measure.real
  rw [Real.volume_pi_Ico_toReal]
  · simp only [show ∀ i : Fin 3,
      L * ((residue (Nat.fib (n + 1)) (Nat.fib n) (b i).val : ℝ) + 1) / Nat.fib n -
      L * (residue (Nat.fib (n + 1)) (Nat.fib n) (b i).val : ℝ) / Nat.fib n = L / Nat.fib n
      from fun _ => by ring, Finset.prod_const, Finset.card_univ, Fintype.card_fin, div_pow]
  · intro i
    have hq : (0 : ℝ) < Nat.fib n := by exact_mod_cast Nat.fib_pos.mpr hn
    apply div_le_div_of_nonneg_right _ hq.le
    nlinarith

theorem cell_measure_finite (n : ℕ) (L : ℝ) (b : Fin 3 → Fin (Nat.fib n)) :
    IsFiniteMeasure (volume.restrict (cell n L b)) := by
  apply isFiniteMeasure_restrict.mpr
  unfold cell axisCell
  rw [Real.volume_pi_Ico]
  exact ENNReal.prod_ne_top (fun _ _ => ENNReal.ofReal_ne_top)

/-- Actual golden-source Lebesgue quadrature on the cube, with no supplied
partition, cell masses, assignment bound or integral-error hypothesis.
The integral uses ordinary Cartesian Lebesgue volume; `toLp 2` identifies
those coordinates with the Euclidean repair carrier. -/
theorem golden_quadrature {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 < L)
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    (f : EuclideanSpace ℝ (Fin 3) → E) {K : ℝ≥0} (hK : LipschitzWith K f) :
    ‖(∑ b : Fin 3 → Fin (Nat.fib n),
      (L ^ 3 / (Nat.fib n : ℝ) ^ 3) • f (site L b)) -
      ∫ x in cube L, f (WithLp.toLp 2 x)‖ ≤
      L ^ 3 * (K : ℝ) * (2 * Real.sqrt 3 * L / Nat.fib n) := by
  classical
  let F : (Fin 3 → ℝ) → E := fun x => f (WithLp.toLp 2 x)
  have hcont : Continuous F := hK.continuous.comp (PiLp.continuous_toLp 2 _)
  have hfcube : IntegrableOn F (cube L) volume := by
    apply (hcont.integrableOn_Icc (a := fun _ => 0) (b := fun _ => L)).mono_set
    intro x hx
    exact ⟨fun i => (hx i (mem_univ i)).1, fun i => (hx i (mem_univ i)).2.le⟩
  have hfcell (b : Fin 3 → Fin (Nat.fib n)) :
      Integrable F (volume.restrict (cell n L b)) := by
    apply hfcube.mono_set
    rw [← cell_cover hn hL]
    exact subset_iUnion (cell n L) b
  letI hfinite (b : Fin 3 → Fin (Nat.fib n)) :
      IsFiniteMeasure (volume.restrict (cell n L b)) := cell_measure_finite n L b
  have hsum : (∑ b : Fin 3 → Fin (Nat.fib n), volume.restrict (cell n L b)) =
      volume.restrict (cube L) := by
    rw [← Measure.sum_fintype, ← Measure.restrict_iUnion (cell_disjoint hn hL)
      (cell_measurable n L), cell_cover hn hL]
  have he (b : Fin 3 → Fin (Nat.fib n)) :
      ‖(L ^ 3 / (Nat.fib n : ℝ) ^ 3) • f (site L b) -
        ∫ x, F x ∂volume.restrict (cell n L b)‖ ≤
      (L ^ 3 / (Nat.fib n : ℝ) ^ 3) * (K : ℝ) *
        (2 * Real.sqrt 3 * L / Nat.fib n) := by
    rw [← cell_mass hn hL.le b, ← measureReal_restrict_apply_univ,
      ← integral_const, ← integral_sub (integrable_const _) (hfcell b)]
    have hb : ∀ᵐ x ∂volume.restrict (cell n L b),
        ‖f (site L b) - F x‖ ≤ (K : ℝ) * (2 * Real.sqrt 3 * L / Nat.fib n) := by
      filter_upwards [ae_restrict_mem (cell_measurable n L b)] with x hx
      have hass : dist (site L b) (WithLp.toLp 2 x) ≤
          2 * Real.sqrt 3 * L / Nat.fib n := by
        apply tensor_assignment hn hL.le b (WithLp.toLp 2 x)
        intro i
        exact ⟨(hx i (mem_univ i)).1, (hx i (mem_univ i)).2.le⟩
      exact (hK.norm_sub_le (site L b) (WithLp.toLp 2 x)).trans
        (mul_le_mul_of_nonneg_left hass K.coe_nonneg)
    simpa only [mul_assoc, mul_comm, mul_left_comm] using norm_integral_le_of_norm_le_const hb
  have hq : (Nat.fib n : ℝ) ≠ 0 := by exact_mod_cast (Nat.fib_pos.mpr hn).ne'
  change ‖_ - ∫ x, F x ∂volume.restrict (cube L)‖ ≤ _
  rw [← hsum, integral_finset_sum_measure (fun b _ => hfcell b), ← Finset.sum_sub_distrib]
  calc
    _ ≤ ∑ _ : Fin 3 → Fin (Nat.fib n),
        (L ^ 3 / (Nat.fib n : ℝ) ^ 3) * (K : ℝ) *
          (2 * Real.sqrt 3 * L / Nat.fib n) :=
      (norm_sum_le _ _).trans (Finset.sum_le_sum fun b _ => he b)
    _ = _ := by
      simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fun,
        Fintype.card_fin, nsmul_eq_mul, Nat.cast_pow]
      field_simp

/-- Refinement convergence for the actual growing Fibonacci source population. -/
theorem golden_quadrature_tendsto {L : ℝ} (hL : 0 < L)
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    (f : EuclideanSpace ℝ (Fin 3) → E) {K : ℝ≥0} (hK : LipschitzWith K f) :
    Filter.Tendsto (fun n => ∑ b : Fin 3 → Fin (Nat.fib (n + 2)),
      (L ^ 3 / (Nat.fib (n + 2) : ℝ) ^ 3) • f (site L b))
      Filter.atTop (nhds (∫ x in cube L, f (WithLp.toLp 2 x))) := by
  apply tendsto_iff_norm_sub_tendsto_zero.mpr
  apply squeeze_zero (fun _ => norm_nonneg _) (fun n =>
    golden_quadrature (by omega : 0 < n + 2) hL f hK)
  have hfib : Filter.Tendsto (fun n => (Nat.fib (n + 2) : ℝ)) Filter.atTop Filter.atTop :=
    tendsto_natCast_atTop_atTop.comp Nat.fib_add_two_strictMono.tendsto_atTop
  have hH := hfib.const_div_atTop (2 * Real.sqrt 3 * L)
  simpa only [mul_zero] using
    (tendsto_const_nhds.mul hH : Filter.Tendsto
      (fun n => (L ^ 3 * (K : ℝ)) * (2 * Real.sqrt 3 * L / Nat.fib (n + 2)))
      Filter.atTop (nhds ((L ^ 3 * (K : ℝ)) * 0)))

/-- Half-open boundary ownership leaves the closed-cube integral unchanged. -/
theorem cube_ae_closed (L : ℝ) :
    cube L =ᵐ[volume] Icc (fun _ : Fin 3 => (0 : ℝ)) (fun _ => L) :=
  Measure.univ_pi_Ico_ae_eq_Icc

theorem golden_quadrature_closed {n : ℕ} (hn : 0 < n) {L : ℝ} (hL : 0 < L)
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    (f : EuclideanSpace ℝ (Fin 3) → E) {K : ℝ≥0} (hK : LipschitzWith K f) :
    ‖(∑ b : Fin 3 → Fin (Nat.fib n),
      (L ^ 3 / (Nat.fib n : ℝ) ^ 3) • f (site L b)) -
      ∫ x in Icc (fun _ : Fin 3 => (0 : ℝ)) (fun _ => L), f (WithLp.toLp 2 x)‖ ≤
      L ^ 3 * (K : ℝ) * (2 * Real.sqrt 3 * L / Nat.fib n) := by
  rw [← setIntegral_congr_set (cube_ae_closed L)]
  exact golden_quadrature hn hL f hK

/-- A sample can lie outside its assigned cell: treating these as Voronoi
cells or deleting the orbit-to-grid term is incorrect. -/
theorem sample_outside_assigned_cell :
    orbit 2 < (1 : ℝ) / 3 ∧ axisCell 4 1 ⟨2, by norm_num⟩ = Ico (1 / 3 : ℝ) (2 / 3) := by
  have hf := golden_floor (n := 4) (b := 2) (by norm_num) (by norm_num)
  norm_num at hf
  have hs : Real.sqrt 5 < 7 / 3 := (Real.sqrt_lt (by norm_num) (by norm_num)).mpr (by norm_num)
  constructor
  · unfold orbit Int.fract
    norm_num only [Nat.cast_ofNat]
    rw [hf]
    change 2 * ((1 + Real.sqrt 5) / 2) - 3 < 1 / 3
    linarith
  · norm_num [axisCell, residue]

/-- One cell width is not a valid global assignment bound for this family. -/
theorem one_cell_width_counterexample :
    (3 / 5 : ℝ) ∈ axisCell 4 1 ⟨2, by norm_num⟩ ∧
      (1 : ℝ) / 3 < |orbit 2 - 3 / 5| := by
  have hf := golden_floor (n := 4) (b := 2) (by norm_num) (by norm_num)
  norm_num at hf
  have hs : Real.sqrt 5 < 34 / 15 := (Real.sqrt_lt (by norm_num) (by norm_num)).mpr (by norm_num)
  constructor
  · norm_num [axisCell, residue]
  · unfold orbit Int.fract
    norm_num only [Nat.cast_ofNat]
    rw [hf]
    norm_num only [Int.cast_ofNat]
    have hneg : (2 : ℝ) * Real.goldenRatio - 3 - 3 / 5 < 0 := by
      change 2 * ((1 + Real.sqrt 5) / 2) - 3 - 3 / 5 < 0
      linarith
    rw [abs_of_neg hneg]
    change 1 / 3 < -(2 * ((1 + Real.sqrt 5) / 2) - 3 - 3 / 5)
    linarith

#print axioms golden_floor
#print axioms residue_bijective
#print axioms tensor_assignment
#print axioms site_card
#print axioms golden_quadrature
#print axioms golden_quadrature_tendsto
#print axioms golden_quadrature_closed
#print axioms one_cell_width_counterexample

end OPH.GoldenSourceAssignment
