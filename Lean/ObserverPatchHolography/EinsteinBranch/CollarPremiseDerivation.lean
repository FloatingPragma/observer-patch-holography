import ObserverPatchHolography.EinsteinBranch.PerCutCollarComposition

/-!
# Per-cut collar premises from independence and cut geometry

`PerCutCollarComposition.lean` derives the deep-regime enclosed-mass law from
a `PerCutCollarModel` whose two load-bearing inputs are premises: additivity
of the per-cut defect variance in the source (`V_additive`) and the linear
cut count `n * r`.  This module derives both one level deeper.

* `AtomicDefectModel` packages a probability space, a family of atomic defect
  contributions `X i : Ω → ℝ` with finite second moments, pairwise
  independence, a common per-atom variance, and a unit mass per atom.
  `compoundVariance_eq_sum` is Mathlib's
  `ProbabilityTheory.IndepFun.variance_sum` applied to the model: the
  variance of the summed defect over any finite sub-family equals the sum of
  the per-atom variances.  `compoundVariance_disjUnion` specializes this to
  disjoint-union composition: compound variances are additive, so the
  `V_additive` premise shape is a theorem of any atomic independent-defect
  model.  `compoundVariance_eq_slope_mul_mass` gives `V = slope * mass` on
  the modeled masses, with `slope = atomVariance / unitMass`.
* `toPerCutCollarModel` extends the measured compound variances linearly to a
  full `PerCutCollarModel`; `toPerCutCollarModel_V_measured` pins the
  extension to the measured value at every modeled mass,
  `toPerCutCollarModel_V_additive_measured` identifies its additivity on
  modeled masses with the measured disjoint-union additivity, and
  `toDeepRegimeLaw` composes through the landed chain (`perCutLaw`) to a
  `DeepRegimeLaw` with dictionary constant `G * n ^ 2 * slope`
  (`toDeepRegimeLaw_dictionary`).
* `SeparatedDenseCutFamily` models a settled-cut family as a sequence of
  radii with positive first cut, consecutive cuts at least `s` apart
  (`separated`), and every length-`s` interval to the right of a nonnegative
  point meeting a cut (`dense`), with `s > 0`.  The premises force the exact
  progression `cut k = cut 0 + k * s` (`cut_eq`).  `count r` is the number of
  cuts in `(0, r]` (`count_spec`, `indexSet_eq`, `ncard_indexSet`); the exact
  sandwich is `r / s - 1 ≤ count r ≤ r / s + 1` (`count_lower`,
  `count_upper`), the count is monotone (`count_mono`), and
  `count r / r → 1 / s` (`count_div_tendsto`).  With per-cut amplitude
  `a ≥ 0` the counted mass `count r * a` deviates from the linear law
  `(r / s) * a` by at most one amplitude (`count_linear_deviation`): the
  linear-cut-count premise holds exactly up to a single-cut boundary term,
  with density `n = 1 / s`.
* `composedCollarModel`, `composedDeepLaw`, `composed_deep_law`,
  `composed_a0_unique`: an atomic independent-defect model on an
  `s`-separated `s`-dense cut family induces the deep-regime law with unique
  constant `a0 = G * (1 / s) ^ 2 * slope`, exact on the idealized linear
  count and within one cut amplitude of the counted enclosed mass at every
  radius.
* Load-bearing receipts and inhabitants: on the fair-coin space
  `headsIndicator` has variance `1 / 4` exactly, and the perfectly correlated
  pair violates variance additivity, `Var (X + X) = 1` against `1 / 2`
  (`correlated_pair_breaks_variance_additivity`), so independence is
  necessary.  `gappyCut` is `1`-separated with positive first cut but not
  `1`-dense, and its cut count within radius `3` is `1`, violating the
  sandwich lower bound `2` (`density_necessary`), so density is necessary.
  `coinPairModel` (two genuinely independent coin atoms on a product space)
  and `unitCutFamily` (cuts at every positive integer, `s = 1`) are full
  explicit inhabitants, composed in `coinPair_unit_dictionary` to the deep
  law with constant `1 / 4`.

What is not proved here: independence of the atomic defect contributions,
their common per-atom variance, and the separated-dense cut geometry are
premises of the two structures; none is derived from the axioms of the
event-algebra corpus.  Equal per-atom variance is weaker than identical
distribution and is the only distributional input used.  The variance
function of the induced `PerCutCollarModel` is the linear extension of the
measured compound variances; between modeled masses it is fixed by extension,
not by measurement.  The passage from the counted cut number to the idealized
linear count is proved as a bounded-error statement, one cut amplitude at
every radius, together with the asymptotic density `count r / r → 1 / s`; no
stronger limit is claimed.  No physical value of `s`, `slope`, or `a0` is
asserted, and no comparison with data is made here.
-/

namespace OPH.EinsteinBranch

open MeasureTheory ProbabilityTheory Filter
open scoped ENNReal

noncomputable section

/-! ## Atomic independent-defect models -/

/-- An atomic independent-defect model for the collar defect of a compound
source.  `X i : Ω → ℝ` is the defect contribution of atom `i` on the
probability space `(Ω, μ)`; the contributions have finite second moments
(`memLp`), are pairwise independent (`indep`), and share the per-atom
variance `atomVariance > 0` (`atomVariance_eq`, weaker than identical
distribution).  Each atom carries the mass `unitMass > 0`, so a finite atom
family `s` models the source mass `s.card * unitMass`. -/
structure AtomicDefectModel where
  /-- Sample space of the collar noise. -/
  Ω : Type*
  /-- Measurable structure on the sample space. -/
  [mΩ : MeasurableSpace Ω]
  /-- Law of the collar noise. -/
  μ : Measure Ω
  /-- The law is a probability measure. -/
  isProb : IsProbabilityMeasure μ
  /-- Atom index type. -/
  ι : Type*
  /-- Defect contribution of each atom. -/
  X : ι → Ω → ℝ
  /-- Finite second moments. -/
  memLp : ∀ i, MemLp (X i) 2 μ
  /-- Pairwise independence of the atomic contributions. -/
  indep : Pairwise fun i j => IndepFun (X i) (X j) μ
  /-- Common per-atom variance. -/
  atomVariance : ℝ
  /-- Every atom has the common variance. -/
  atomVariance_eq : ∀ i, Var[X i; μ] = atomVariance
  /-- The atoms are not deterministic. -/
  atomVariance_pos : 0 < atomVariance
  /-- Mass carried by a single atom. -/
  unitMass : ℝ
  /-- The unit mass is positive. -/
  unitMass_pos : 0 < unitMass

namespace AtomicDefectModel

variable (M : AtomicDefectModel)

/-- Compound defect of a finite atom family: the summed contribution. -/
def compound (s : Finset M.ι) : M.Ω → ℝ := ∑ i ∈ s, M.X i

/-- Variance of the compound defect of a finite atom family. -/
def compoundVariance (s : Finset M.ι) : ℝ := Var[M.compound s; M.μ]

/-- Modeled mass of a finite atom family: atom count times unit mass. -/
def mass (s : Finset M.ι) : ℝ := s.card * M.unitMass

/-- Per-unit-mass variance slope. -/
def slope : ℝ := M.atomVariance / M.unitMass

/-- The variance slope is strictly positive. -/
theorem slope_pos : 0 < M.slope := div_pos M.atomVariance_pos M.unitMass_pos

/-- Variance additivity from independence: the variance of the compound
defect over any finite sub-family is the sum of the per-atom variances.
This is Mathlib's `ProbabilityTheory.IndepFun.variance_sum` applied to the
model; no measure theory is reproved here. -/
theorem compoundVariance_eq_sum (s : Finset M.ι) :
    M.compoundVariance s = ∑ i ∈ s, Var[M.X i; M.μ] := by
  unfold compoundVariance compound
  exact IndepFun.variance_sum (fun i _ => M.memLp i) (M.indep.set_pairwise ↑s)

/-- Disjoint-union composition of atom families makes the compound variance
additive.  This is the `V_additive` premise shape of `PerCutCollarModel`,
derived from probabilistic independence instead of assumed. -/
theorem compoundVariance_disjUnion (s t : Finset M.ι) (hst : Disjoint s t) :
    M.compoundVariance (s.disjUnion t hst) =
      M.compoundVariance s + M.compoundVariance t := by
  rw [M.compoundVariance_eq_sum, M.compoundVariance_eq_sum,
    M.compoundVariance_eq_sum, Finset.sum_disjUnion]

/-- Modeled mass is additive under disjoint union of atom families. -/
theorem mass_disjUnion (s t : Finset M.ι) (hst : Disjoint s t) :
    M.mass (s.disjUnion t hst) = M.mass s + M.mass t := by
  unfold mass
  rw [Finset.card_disjUnion]
  push_cast
  ring

/-- Equal-variance atoms give the linear variance law on modeled masses:
the compound variance is the slope times the modeled mass. -/
theorem compoundVariance_eq_slope_mul_mass (s : Finset M.ι) :
    M.compoundVariance s = M.slope * M.mass s := by
  rw [M.compoundVariance_eq_sum s,
    Finset.sum_congr rfl fun i _ => M.atomVariance_eq i,
    Finset.sum_const, nsmul_eq_mul]
  unfold slope mass
  have h : M.unitMass ≠ 0 := M.unitMass_pos.ne'
  field_simp

/-- The per-cut collar model induced by an atomic model: the variance
function is the linear extension `V m = slope * m` of the measured compound
variances, and every premise field of `PerCutCollarModel` on `V` is
discharged by a proof rather than assumed.  The collar geometry data `G` and
`n` are threaded as hypotheses. -/
def toPerCutCollarModel (G n : ℝ) (hG : 0 < G) (hn : 0 < n) :
    PerCutCollarModel where
  G := G
  G_pos := hG
  n := n
  n_pos := hn
  V := fun m => M.slope * m
  V_additive := fun _ _ _ _ => by ring
  V_mono := fun _ _ _ h12 => mul_le_mul_of_nonneg_left h12 M.slope_pos.le
  V_nonneg := fun _ hm => mul_nonneg M.slope_pos.le hm.le
  V_nondeg := ⟨1, one_pos, by simpa using M.slope_pos⟩

/-- The induced variance function agrees with the measured compound variance
at every modeled mass. -/
theorem toPerCutCollarModel_V_measured (G n : ℝ) (hG : 0 < G) (hn : 0 < n)
    (s : Finset M.ι) :
    (M.toPerCutCollarModel G n hG hn).V (M.mass s) = M.compoundVariance s :=
  (M.compoundVariance_eq_slope_mul_mass s).symm

/-- On modeled masses the additivity of the induced variance function is the
measured disjoint-union additivity of `compoundVariance_disjUnion`, so the
`V_additive` field of the induced model is a theorem of independence at every
modeled mass, not only a property of the linear extension. -/
theorem toPerCutCollarModel_V_additive_measured (G n : ℝ) (hG : 0 < G)
    (hn : 0 < n) (s t : Finset M.ι) (hst : Disjoint s t) :
    (M.toPerCutCollarModel G n hG hn).V (M.mass s + M.mass t) =
      M.compoundVariance s + M.compoundVariance t := by
  rw [← M.mass_disjUnion s t hst, M.toPerCutCollarModel_V_measured G n hG hn,
    M.compoundVariance_disjUnion s t hst]

/-- The induced model has variance slope `slope` on all positive masses. -/
theorem toPerCutCollarModel_varianceSlope (G n : ℝ) (hG : 0 < G)
    (hn : 0 < n) :
    ∀ m : ℝ, 0 < m → (M.toPerCutCollarModel G n hG hn).V m = M.slope * m :=
  fun _ _ => rfl

/-- The deep-regime law induced by an atomic model through the landed chain
of `PerCutCollarComposition.lean`. -/
def toDeepRegimeLaw (G n : ℝ) (hG : 0 < G) (hn : 0 < n) : DeepRegimeLaw :=
  perCutLaw (M.toPerCutCollarModel G n hG hn)

/-- Dictionary for the induced deep-regime law: it realizes the committed
profile with the positive constant `a0 = G * n ^ 2 * slope`. -/
theorem toDeepRegimeLaw_dictionary (G n : ℝ) (hG : 0 < G) (hn : 0 < n) :
    0 < G * n ^ 2 * M.slope ∧
      ∀ Mb r : ℝ, 0 < Mb → 0 < r →
        (M.toDeepRegimeLaw G n hG hn).M Mb r =
          r * Real.sqrt (Mb * (G * n ^ 2 * M.slope) / G) :=
  perCut_a0_dictionary (M.toPerCutCollarModel G n hG hn) M.slope M.slope_pos
    (fun _ _ => rfl)

end AtomicDefectModel

/-! ## Separated dense cut families -/

/-- A settled-cut family around a point source: a sequence of cut radii with
positive innermost cut, consecutive cuts at least `s` apart (`separated`),
and every interval `(a, a + s]` with `a ≥ 0` meeting a cut (`dense`), for a
fixed separation scale `s > 0`. -/
structure SeparatedDenseCutFamily where
  /-- Separation and density scale. -/
  s : ℝ
  /-- The scale is positive. -/
  s_pos : 0 < s
  /-- Radius of the `k`-th cut. -/
  cut : ℕ → ℝ
  /-- The innermost cut has positive radius. -/
  first_pos : 0 < cut 0
  /-- Consecutive cuts are at least `s` apart. -/
  separated : ∀ k, cut k + s ≤ cut (k + 1)
  /-- Every interval `(a, a + s]` with `a ≥ 0` contains a cut. -/
  dense : ∀ a : ℝ, 0 ≤ a → ∃ k, a < cut k ∧ cut k ≤ a + s

namespace SeparatedDenseCutFamily

variable (C : SeparatedDenseCutFamily)

/-- The cut radii are strictly increasing. -/
theorem cut_strictMono : StrictMono C.cut := by
  apply strictMono_nat_of_lt_succ
  intro k
  have h := C.separated k
  have hs := C.s_pos
  linarith

/-- Every cut has positive radius. -/
theorem cut_pos (k : ℕ) : 0 < C.cut k :=
  lt_of_lt_of_le C.first_pos (C.cut_strictMono.monotone (Nat.zero_le k))

/-- Density places the innermost cut within `s` of the source. -/
theorem first_le : C.cut 0 ≤ C.s := by
  obtain ⟨k, _, h2⟩ := C.dense 0 le_rfl
  have h3 : C.cut 0 ≤ C.cut k := C.cut_strictMono.monotone (Nat.zero_le k)
  linarith

/-- Density bounds consecutive gaps above by `s`. -/
theorem gap_le (k : ℕ) : C.cut (k + 1) ≤ C.cut k + C.s := by
  obtain ⟨j, h1, h2⟩ := C.dense (C.cut k) (C.cut_pos k).le
  have hkj : k < j := C.cut_strictMono.lt_iff_lt.mp h1
  have h3 : C.cut (k + 1) ≤ C.cut j :=
    C.cut_strictMono.monotone (by omega : k + 1 ≤ j)
  linarith

/-- Separation and density together force the exact arithmetic progression:
`cut k = cut 0 + k * s`. -/
theorem cut_eq (k : ℕ) : C.cut k = C.cut 0 + k * C.s := by
  induction k with
  | zero => simp
  | succ n ih =>
    have h1 := C.separated n
    have h2 := C.gap_le n
    have h3 : C.cut (n + 1) = C.cut n + C.s := le_antisymm h2 h1
    rw [h3, ih]
    push_cast
    ring

/-- Number of cuts with radius in `(0, r]`.  Since all cut radii are
positive and strictly increasing, this is the number of indices `k` with
`cut k ≤ r`; `count_spec` and `ncard_indexSet` prove exactly that. -/
def count (r : ℝ) : ℕ :=
  if C.cut 0 ≤ r then ⌊(r - C.cut 0) / C.s⌋₊ + 1 else 0

/-- `count` counts: cut `k` lies within radius `r` exactly when
`k < count r`. -/
theorem count_spec (r : ℝ) (k : ℕ) : C.cut k ≤ r ↔ k < C.count r := by
  unfold count
  split_ifs with h0
  · rw [Nat.lt_succ_iff,
      Nat.le_floor_iff (div_nonneg (by linarith) C.s_pos.le),
      le_div_iff₀ C.s_pos, C.cut_eq k]
    constructor <;> intro h <;> linarith
  · constructor
    · intro h
      exact absurd ((C.cut_strictMono.monotone (Nat.zero_le k)).trans h) h0
    · intro h
      exact absurd h (Nat.not_lt_zero k)

/-- The indices of cuts within radius `r` form the initial segment of length
`count r`. -/
theorem indexSet_eq (r : ℝ) :
    {k : ℕ | C.cut k ≤ r} = Set.Iio (C.count r) := by
  ext k
  simpa [Set.mem_Iio] using C.count_spec r k

/-- `count r` is the cardinality of the set of cuts within radius `r`. -/
theorem ncard_indexSet (r : ℝ) :
    {k : ℕ | C.cut k ≤ r}.ncard = C.count r := by
  rw [C.indexSet_eq r]
  exact Set.ncard_Iio_nat _

/-- Lower half of the exact sandwich: `r / s - 1 ≤ count r`, at every
radius. -/
theorem count_lower (r : ℝ) : r / C.s - 1 ≤ (C.count r : ℝ) := by
  have hs := C.s_pos
  unfold count
  split_ifs with h0
  · have h1 : (r - C.cut 0) / C.s < (⌊(r - C.cut 0) / C.s⌋₊ : ℝ) + 1 :=
      Nat.lt_floor_add_one _
    have h2 : r / C.s - 1 ≤ (r - C.cut 0) / C.s := by
      rw [sub_div]
      have h3 : C.cut 0 / C.s ≤ 1 := (div_le_one hs).mpr C.first_le
      linarith
    push_cast
    linarith
  · have h1 : r < C.s := lt_of_lt_of_le (not_le.mp h0) C.first_le
    have h2 : r / C.s < 1 := (div_lt_one hs).mpr h1
    push_cast
    linarith

/-- Upper half of the exact sandwich: `count r ≤ r / s + 1`. -/
theorem count_upper {r : ℝ} (hr : 0 < r) : (C.count r : ℝ) ≤ r / C.s + 1 := by
  have hs := C.s_pos
  unfold count
  split_ifs with h0
  · have h1 : (⌊(r - C.cut 0) / C.s⌋₊ : ℝ) ≤ (r - C.cut 0) / C.s :=
      Nat.floor_le (div_nonneg (by linarith) hs.le)
    have h2 : (r - C.cut 0) / C.s ≤ r / C.s := by
      rw [sub_div]
      have h3 : 0 ≤ C.cut 0 / C.s := div_nonneg C.first_pos.le hs.le
      linarith
    push_cast
    linarith
  · have h1 : 0 < r / C.s := div_pos hr hs
    push_cast
    linarith

/-- The cut count does not decrease with the radius. -/
theorem count_mono {r1 r2 : ℝ} (h : r1 ≤ r2) : C.count r1 ≤ C.count r2 := by
  by_contra hlt
  push Not at hlt
  have h1 : C.cut (C.count r2) ≤ r1 := (C.count_spec r1 _).mpr hlt
  have h2 : ¬ C.cut (C.count r2) ≤ r2 :=
    fun hc => lt_irrefl _ ((C.count_spec r2 _).mp hc)
  exact h2 (h1.trans h)

/-- Asymptotic cut density: `count r / r → 1 / s` as `r → ∞`. -/
theorem count_div_tendsto :
    Tendsto (fun r : ℝ => (C.count r : ℝ) / r) atTop (nhds (1 / C.s)) := by
  have hs := C.s_pos
  have hinv : Tendsto (fun r : ℝ => 1 / r) atTop (nhds 0) := by
    simpa only [one_div] using tendsto_inv_atTop_zero
  have hlow : Tendsto (fun r : ℝ => 1 / C.s - 1 / r) atTop (nhds (1 / C.s)) := by
    simpa using tendsto_const_nhds.sub hinv
  have hhigh : Tendsto (fun r : ℝ => 1 / C.s + 1 / r) atTop (nhds (1 / C.s)) := by
    simpa using tendsto_const_nhds.add hinv
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hlow hhigh ?_ ?_
  · filter_upwards [eventually_gt_atTop (0 : ℝ)] with r hr
    have h1 := C.count_lower r
    have hr0 : r ≠ 0 := hr.ne'
    have hs0 : C.s ≠ 0 := hs.ne'
    have key : (r / C.s - 1) / r = 1 / C.s - 1 / r := by
      field_simp
    have hsub : ((C.count r : ℝ) - (r / C.s - 1)) / r
        = (C.count r : ℝ) / r - (r / C.s - 1) / r := sub_div _ _ _
    have hnn : 0 ≤ ((C.count r : ℝ) - (r / C.s - 1)) / r :=
      div_nonneg (by linarith) hr.le
    rw [← key]
    linarith
  · filter_upwards [eventually_gt_atTop (0 : ℝ)] with r hr
    have h2 := C.count_upper hr
    have hr0 : r ≠ 0 := hr.ne'
    have hs0 : C.s ≠ 0 := hs.ne'
    have key : (r / C.s + 1) / r = 1 / C.s + 1 / r := by
      field_simp
    have hsub : ((r / C.s + 1) - (C.count r : ℝ)) / r
        = (r / C.s + 1) / r - (C.count r : ℝ) / r := sub_div _ _ _
    have hnn : 0 ≤ ((r / C.s + 1) - (C.count r : ℝ)) / r :=
      div_nonneg (by linarith) hr.le
    rw [← key]
    linarith

/-- Bounded boundary term: with per-cut amplitude `a ≥ 0`, the counted mass
`count r * a` deviates from the linear law `(r / s) * a` by at most one
amplitude, at every positive radius. -/
theorem count_linear_deviation {a : ℝ} (ha : 0 ≤ a) {r : ℝ} (hr : 0 < r) :
    |(C.count r : ℝ) * a - r / C.s * a| ≤ a := by
  have e1 : (r / C.s - 1) * a ≤ (C.count r : ℝ) * a :=
    mul_le_mul_of_nonneg_right (C.count_lower r) ha
  have e2 : (C.count r : ℝ) * a ≤ (r / C.s + 1) * a :=
    mul_le_mul_of_nonneg_right (C.count_upper hr) ha
  rw [abs_le]
  constructor <;> nlinarith

end SeparatedDenseCutFamily

/-! ## Composition: atomic model on a cut family -/

/-- Per-cut collar model composed from an atomic independent-defect model and
a separated dense cut family: cut density `n = 1 / s`, variance slope from
the atomic model. -/
def composedCollarModel (M : AtomicDefectModel) (C : SeparatedDenseCutFamily)
    (G : ℝ) (hG : 0 < G) : PerCutCollarModel :=
  M.toPerCutCollarModel G (1 / C.s) hG (by have := C.s_pos; positivity)

/-- Deep-regime law composed from an atomic model and a cut family through
the landed chain. -/
def composedDeepLaw (M : AtomicDefectModel) (C : SeparatedDenseCutFamily)
    (G : ℝ) (hG : 0 < G) : DeepRegimeLaw :=
  perCutLaw (composedCollarModel M C G hG)

/-- Composed premise discharge: the law induced by an atomic model on an
`s`-separated `s`-dense cut family realizes the committed deep profile with
positive constant `a0 = G * (1 / s) ^ 2 * slope`, exactly on the idealized
linear count; and at every radius its enclosed mass is within one cut
amplitude of the counted enclosed mass `count r * cutAmplitude`. -/
theorem composed_deep_law (M : AtomicDefectModel) (C : SeparatedDenseCutFamily)
    (G : ℝ) (hG : 0 < G) :
    0 < G * (1 / C.s) ^ 2 * M.slope ∧
      (∀ Mb r : ℝ, 0 < Mb → 0 < r →
        (composedDeepLaw M C G hG).M Mb r =
          r * Real.sqrt (Mb * (G * (1 / C.s) ^ 2 * M.slope) / G)) ∧
      ∀ Mb r : ℝ, 0 < Mb → 0 < r →
        |(C.count r : ℝ) * cutAmplitude (composedCollarModel M C G hG) Mb -
            (composedDeepLaw M C G hG).M Mb r| ≤
          cutAmplitude (composedCollarModel M C G hG) Mb := by
  obtain ⟨hpos, hdict⟩ := perCut_a0_dictionary (composedCollarModel M C G hG)
    M.slope M.slope_pos (fun _ _ => rfl)
  refine ⟨hpos, hdict, ?_⟩
  intro Mb r hMb hr
  have hM : (composedDeepLaw M C G hG).M Mb r =
      r / C.s * cutAmplitude (composedCollarModel M C G hG) Mb := by
    show 1 / C.s * r * cutAmplitude (composedCollarModel M C G hG) Mb =
      r / C.s * cutAmplitude (composedCollarModel M C G hG) Mb
    ring
  rw [hM]
  exact C.count_linear_deviation (Real.sqrt_nonneg _) hr

/-- Uniqueness of the composed constant: any constant realizing the
committed profile for the composed law equals `G * (1 / s) ^ 2 * slope`. -/
theorem composed_a0_unique (M : AtomicDefectModel) (C : SeparatedDenseCutFamily)
    (G : ℝ) (hG : 0 < G) (a0 : ℝ) (ha0 : 0 < a0)
    (hchar : ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      (composedDeepLaw M C G hG).M Mb r =
        r * Real.sqrt (Mb * a0 / (composedDeepLaw M C G hG).G)) :
    a0 = G * (1 / C.s) ^ 2 * M.slope :=
  perCut_a0_unique (composedCollarModel M C G hG) M.slope M.slope_pos
    (fun _ _ => rfl) a0 ha0 hchar

/-! ## The fair-coin space and the independence receipt -/

/-- Fair-coin measure on `Bool`: equal weight `1 / 2` on each outcome. -/
def coinMeasure : Measure Bool :=
  (2 : ℝ≥0∞)⁻¹ • (Measure.dirac false + Measure.dirac true)

instance : IsProbabilityMeasure coinMeasure := by
  refine ⟨?_⟩
  simp only [coinMeasure, Measure.smul_apply, Measure.add_apply, measure_univ,
    smul_eq_mul]
  rw [one_add_one_eq_two]
  exact ENNReal.inv_mul_cancel two_ne_zero ENNReal.ofNat_ne_top

/-- Integration over the fair coin is the two-point average. -/
theorem coinMeasure_integral (f : Bool → ℝ) :
    ∫ b, f b ∂coinMeasure = (f false + f true) / 2 := by
  unfold coinMeasure
  rw [integral_smul_measure,
    integral_add_measure Integrable.of_finite Integrable.of_finite,
    integral_dirac, integral_dirac]
  simp only [ENNReal.toReal_inv, ENNReal.toReal_ofNat, smul_eq_mul]
  ring

/-- Defect contribution of a single coin atom: `1` on heads, `0` on tails. -/
def headsIndicator : Bool → ℝ := fun b => if b then 1 else 0

/-- The coin atom has mean `1 / 2` exactly. -/
theorem headsIndicator_mean :
    ∫ b, headsIndicator b ∂coinMeasure = 1 / 2 := by
  rw [coinMeasure_integral]
  norm_num [headsIndicator]

/-- The coin atom has variance `1 / 4` exactly. -/
theorem headsIndicator_variance : Var[headsIndicator; coinMeasure] = 1 / 4 := by
  rw [variance_eq_integral Measurable.of_discrete.aemeasurable,
    coinMeasure_integral, headsIndicator_mean]
  norm_num [headsIndicator]

/-- Variance of a self-sum: doubling a single source quadruples the
variance, for any random variable on any measure space. -/
theorem variance_self_add_self {Ω : Type*} [MeasurableSpace Ω]
    (μ : Measure Ω) (X : Ω → ℝ) : Var[X + X; μ] = 4 * Var[X; μ] := by
  have h : X + X = fun ω => 2 * X ω := by
    funext ω
    simp [two_mul]
  rw [h, variance_const_mul]
  norm_num

/-- Load-bearing receipt for independence: the perfectly correlated pair
`(headsIndicator, headsIndicator)` on the fair coin violates variance
additivity of the sum.  The self-sum has variance `1`, amplitude-quadrature
additivity would demand `1 / 2`, and the two differ.  A two-atom family with
both contributions equal to `headsIndicator` satisfies every field of
`AtomicDefectModel` except `indep`, and the conclusion of
`compoundVariance_eq_sum` fails for it: independence is necessary. -/
theorem correlated_pair_breaks_variance_additivity :
    Var[headsIndicator + headsIndicator; coinMeasure] = 1 ∧
      Var[headsIndicator; coinMeasure] + Var[headsIndicator; coinMeasure]
        = 1 / 2 ∧
      Var[headsIndicator + headsIndicator; coinMeasure] ≠
        Var[headsIndicator; coinMeasure] + Var[headsIndicator; coinMeasure] := by
  have h1 : Var[headsIndicator + headsIndicator; coinMeasure] = 1 := by
    rw [variance_self_add_self, headsIndicator_variance]
    norm_num
  have h2 : Var[headsIndicator; coinMeasure] + Var[headsIndicator; coinMeasure]
      = 1 / 2 := by
    rw [headsIndicator_variance]
    norm_num
  refine ⟨h1, h2, ?_⟩
  rw [h1, h2]
  norm_num

/-! ## The density receipt -/

/-- Sparse cut sequence at radii `1, 4, 7, ...`: `1`-separated with positive
first cut, but not `1`-dense. -/
def gappyCut : ℕ → ℝ := fun k => 1 + 3 * (k : ℝ)

/-- The sparse sequence has a positive first cut. -/
theorem gappyCut_first_pos : 0 < gappyCut 0 := by
  norm_num [gappyCut]

/-- The sparse sequence is `1`-separated. -/
theorem gappyCut_separated (k : ℕ) : gappyCut k + 1 ≤ gappyCut (k + 1) := by
  unfold gappyCut
  push_cast
  linarith

/-- The sparse sequence is not `1`-dense: no cut lies in `(1, 2]`. -/
theorem gappyCut_not_dense :
    ¬ ∀ a : ℝ, 0 ≤ a → ∃ k, a < gappyCut k ∧ gappyCut k ≤ a + 1 := by
  intro hdense
  obtain ⟨k, h1, h2⟩ := hdense 1 zero_le_one
  unfold gappyCut at h1 h2
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · norm_num at h1
  · have hk1 : (1 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk
    linarith

/-- Exactly one sparse cut lies within radius `3`. -/
theorem gappyCut_indexSet : {k : ℕ | gappyCut k ≤ 3} = Set.Iio 1 := by
  ext k
  simp only [Set.mem_setOf_eq, Set.mem_Iio, gappyCut]
  constructor
  · intro h
    by_contra hk
    push Not at hk
    have hk1 : (1 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk
    linarith
  · intro hk
    have hk0 : k = 0 := by omega
    subst hk0
    norm_num

/-- Load-bearing receipt for density: the sparse sequence satisfies every
`SeparatedDenseCutFamily` premise at `s = 1` except density, its cut count
within radius `3` is `1`, and `1` violates the sandwich lower bound
`3 / 1 - 1 = 2`: density is necessary for `count_lower`. -/
theorem density_necessary :
    0 < gappyCut 0 ∧ (∀ k, gappyCut k + 1 ≤ gappyCut (k + 1)) ∧
      (¬ ∀ a : ℝ, 0 ≤ a → ∃ k, a < gappyCut k ∧ gappyCut k ≤ a + 1) ∧
      {k : ℕ | gappyCut k ≤ 3}.ncard = 1 ∧
      ¬ ((3 : ℝ) / 1 - 1 ≤ 1) := by
  refine ⟨gappyCut_first_pos, gappyCut_separated, gappyCut_not_dense, ?_,
    by norm_num⟩
  rw [gappyCut_indexSet]
  exact Set.ncard_Iio_nat 1

/-! ## Explicit inhabitants -/

/-- Explicit inhabitant of `SeparatedDenseCutFamily`: cuts at every positive
integer radius, `s = 1`. -/
def unitCutFamily : SeparatedDenseCutFamily where
  s := 1
  s_pos := one_pos
  cut := fun k => (k : ℝ) + 1
  first_pos := by norm_num
  separated := fun k => by push_cast; linarith
  dense := fun a ha => by
    refine ⟨⌊a⌋₊, ?_, ?_⟩
    · exact Nat.lt_floor_add_one a
    · have h := Nat.floor_le ha
      linarith

/-- Exact count receipt for the inhabitant: two cuts lie within radius `2`. -/
theorem unitCutFamily_count_two : unitCutFamily.count 2 = 2 := by
  unfold SeparatedDenseCutFamily.count
  norm_num [unitCutFamily]

/-- Explicit inhabitant of `AtomicDefectModel`: two genuinely independent
coin atoms on the product of two fair coins, each contributing
`headsIndicator` of its own coordinate, with per-atom variance `1 / 4` and
unit mass `1`. -/
def coinPairModel : AtomicDefectModel where
  Ω := Fin 2 → Bool
  mΩ := inferInstance
  μ := Measure.pi fun _ => coinMeasure
  isProb := inferInstance
  ι := Fin 2
  X := fun i ω => headsIndicator (ω i)
  memLp := fun i => by
    have h : MemLp (headsIndicator ∘ Function.eval i) 2
        (Measure.pi fun _ : Fin 2 => coinMeasure) :=
      MemLp.of_discrete.comp_measurePreserving
        (measurePreserving_eval (fun _ : Fin 2 => coinMeasure) i)
    exact h
  indep := fun i j hij => by
    have h : iIndepFun (fun (k : Fin 2) (ω : Fin 2 → Bool) => headsIndicator (ω k))
        (Measure.pi fun _ : Fin 2 => coinMeasure) :=
      iIndepFun_pi (X := fun _ : Fin 2 => headsIndicator)
        fun _ => Measurable.of_discrete.aemeasurable
    exact h.indepFun hij
  atomVariance := 1 / 4
  atomVariance_eq := fun i => by
    have h := (measurePreserving_eval (fun _ : Fin 2 => coinMeasure) i).variance_fun_comp
      (f := headsIndicator) Measurable.of_discrete.aemeasurable
    rw [headsIndicator_variance] at h
    exact h
  atomVariance_pos := by norm_num
  unitMass := 1
  unitMass_pos := one_pos

/-- Measured two-atom receipt, through the derived additivity: the compound
variance of both coin atoms is `1 / 2`, the slope `1 / 4` times the modeled
mass `2`. -/
theorem coinPairModel_compoundVariance_univ :
    coinPairModel.compoundVariance (Finset.univ : Finset (Fin 2)) = 1 / 2 := by
  rw [AtomicDefectModel.compoundVariance_eq_slope_mul_mass]
  have hslope : coinPairModel.slope = 1 / 4 := by
    show (1 / 4 : ℝ) / 1 = 1 / 4
    norm_num
  have hmass : coinPairModel.mass (Finset.univ : Finset (Fin 2)) = 2 := by
    show ((Finset.univ : Finset (Fin 2)).card : ℝ) * 1 = 2
    simp
  rw [hslope, hmass]
  norm_num

/-- Fully explicit composed inhabitant: the two-coin atomic model on the
unit cut family at `G = 1` induces the deep-regime law with constant
`a0 = 1 / 4` exactly. -/
theorem coinPair_unit_dictionary :
    ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      (composedDeepLaw coinPairModel unitCutFamily 1 one_pos).M Mb r =
        r * Real.sqrt (Mb * (1 / 4)) := by
  intro Mb r hMb hr
  have h := (composed_deep_law coinPairModel unitCutFamily 1 one_pos).2.1
    Mb r hMb hr
  have hslope : coinPairModel.slope = 1 / 4 := by
    show (1 / 4 : ℝ) / 1 = 1 / 4
    norm_num
  have hs : unitCutFamily.s = 1 := rfl
  rw [h, hslope, hs]
  norm_num

/-! ## Per-theorem axiom audit -/

#print axioms AtomicDefectModel.slope_pos
#print axioms AtomicDefectModel.compoundVariance_eq_sum
#print axioms AtomicDefectModel.compoundVariance_disjUnion
#print axioms AtomicDefectModel.mass_disjUnion
#print axioms AtomicDefectModel.compoundVariance_eq_slope_mul_mass
#print axioms AtomicDefectModel.toPerCutCollarModel_V_measured
#print axioms AtomicDefectModel.toPerCutCollarModel_V_additive_measured
#print axioms AtomicDefectModel.toPerCutCollarModel_varianceSlope
#print axioms AtomicDefectModel.toDeepRegimeLaw_dictionary
#print axioms SeparatedDenseCutFamily.cut_strictMono
#print axioms SeparatedDenseCutFamily.cut_pos
#print axioms SeparatedDenseCutFamily.first_le
#print axioms SeparatedDenseCutFamily.gap_le
#print axioms SeparatedDenseCutFamily.cut_eq
#print axioms SeparatedDenseCutFamily.count_spec
#print axioms SeparatedDenseCutFamily.indexSet_eq
#print axioms SeparatedDenseCutFamily.ncard_indexSet
#print axioms SeparatedDenseCutFamily.count_lower
#print axioms SeparatedDenseCutFamily.count_upper
#print axioms SeparatedDenseCutFamily.count_mono
#print axioms SeparatedDenseCutFamily.count_div_tendsto
#print axioms SeparatedDenseCutFamily.count_linear_deviation
#print axioms composed_deep_law
#print axioms composed_a0_unique
#print axioms coinMeasure_integral
#print axioms headsIndicator_mean
#print axioms headsIndicator_variance
#print axioms variance_self_add_self
#print axioms correlated_pair_breaks_variance_additivity
#print axioms gappyCut_first_pos
#print axioms gappyCut_separated
#print axioms gappyCut_not_dense
#print axioms gappyCut_indexSet
#print axioms density_necessary
#print axioms unitCutFamily_count_two
#print axioms coinPairModel_compoundVariance_univ
#print axioms coinPair_unit_dictionary

end

end OPH.EinsteinBranch
