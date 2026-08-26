import Geometry.SourceClockRateAlongWorldlines

set_option autoImplicit false

open scoped BigOperators

/-!
# The proper-length-clocked source chain along a seam-step worldline
(issues 736, 739)

STATUS.  Enrichment module on the source clock and duration row.  The
comparison module `Geometry/SourceClockRateAlongWorldlines.lean` places
index accrual and proper-length accrual side by side as two declared
readings of the abstract chain step index; the present module declares
the proper-length reading as a clocking rule with an integer output, the
clocked index `c(n) = ⌊L(n) / u⌋₊` (one chain step each time the
accumulated proper length `L(n)` crosses a multiple of a declared
internal unit `u > 0`), composes the abstract chain trajectory with it,
and proves the exact bounds, the limit rate, the return-count
conversion, the refinement invariance, and the scope statements.  The
unit `u`, the step unit `τ`, and the clocking rule are declared; the
module selects none of them from the source, and it discharges no row.

WHAT IS PROVED.

(1) Declared clocking.  `accumulatedLength τ w n = properLength n
(generatedPath τ w)` and `clockedIndex τ u w n = ⌊accumulatedLength τ w n
/ u⌋₊`.  The accumulated length is nonnegative and monotone
(`accumulatedLength_nonneg`, `accumulatedLength_mono`), so the clocked
index is monotone (`clockedIndex_mono`).  When the proper length of step
`n` is at most `u` the clocked index advances by `0` or `1` at that step
(`clockedIndex_succ_of_step_le`); at every unit `0 ≤ τ ≤ u` the step
length of a generated path is at most `τ` (`generated_step_le`), so the
increment lies in `{0, 1}` at every step of every seam-step worldline
(`clockedIndex_unit_increments`).  On the resting worldline with `u = τ`
the clocked index is the index, `c(n) = n` (`clockedIndex_rest`): the
index rule is the clocked rule at rest, exactly.

(2) Uniform worldline.  With `n = (k + 1) p + j`, `0 ≤ j ≤ k`, the
crossing count of `uniformWorldline e k` over `n` steps is `p` when `j =
0` and `p + 1` otherwise (`crossingCount_uniform_general`), so the
accumulated length satisfies `n d τ - (τ - √(τ² - 4)) ≤ L(n) ≤ n d τ` with
`d = dilationFactor τ k` (`accumulatedLength_uniform_le`,
`accumulatedLength_uniform_ge`), with equality on full periods
(`accumulatedLength_uniform_period`).  The clocked index satisfies
`n d τ / u - (τ - √(τ² - 4)) / u - 1 < c(n) ≤ n d τ / u`
(`clockedIndex_uniform_bounds`) and on full periods the two-sided bound
`n d τ / u - 1 < c(n) ≤ n d τ / u` (`clockedIndex_uniform_period_bounds`).
Hence `c(n) / n → d τ / u` (`clockedIndex_uniform_rate`), for `u = τ`
the rate is exactly `d` (`clockedIndex_uniform_rate_unit`), `d < 1` at
every timelike unit (`dilationFactor_lt_one`, cited) and the rate at rest
is `1` (`clockedIndex_rest_rate`).

(3) The clocked chain.  `clockedTrajectory ω c = ω ∘ c` and the clocked
return count `clockedReturns i ω c n`, the number of index steps `j < n`
at which the clock advances and the chain enters state `i`.  For a clock
with `c 0 = 0` and unit increments the clocked return count over `n`
index steps equals `returnCount i ω (c n)`, the chain's own count over
`c n` chain steps (`clockedReturns_eq_returnCount`), and the composed
count is block-additive across index windows through the chain's
`returnCount_block_additive` (`returnCount_clocked_block_additive`).  Any
quantity with a constant rate `a > 0` per chain step converts with the
bound of (2) (`perChainStep_conversion`); under the declared reading of
the stationary mass `chainPi 0 = 1 / meanReturn 0` (`kac_identity`) as
the expected returns per chain step, the expected returns over `n` index
steps on the uniform worldline lie in
`(chainPi 0 * (n d τ / u - 1), chainPi 0 * n d τ / u]` on full periods
(`expectedReturns_uniform_period_bounds`) and the expected returns per
index step tend to `(d τ / u) / meanReturn 0`
(`expectedReturns_uniform_rate`), the reciprocal of
`returnIndexSteps τ u k = meanReturn 0 * u / (d τ)`, which equals
`(61511 / 7155) * u / (d τ)` (`returnIndexSteps_value`) and agrees with
`lengthRuleReturnIndexSteps τ k` at `u = τ`
(`returnIndexSteps_unit`).

(4) Refinement and the join.  The clocked index depends on the path only
through the proper length (`clockedIndexPath_congr`), so it is invariant
under the declared midpoint refinement of the Lorentz-module image at the
corresponding index `2 n` (`clockedIndexPath_refine`, through
`properLength_refine`), and so is the clocked return count
(`clockedReturns_refine`).  The join's clock advances per index at a
declared `stepDuration` (`joined_clock_advance`, cited), that is through
`stepTime δ n = n δ`; on the uniform worldline with `k ≥ 1` no `δ`
reproduces the accumulated proper length at every index
(`no_stepDuration_matches_uniform`), while at `k = 0` the constant
`δ = √(τ² - 4)` does (`stepDuration_matches_uniform_zero`): the clocked
chain is a declared enrichment of the join in which the proper-time
argument `stepTime J.stepDuration n` of the phase would have to be
replaced by the accumulated proper length.  Two units give two clocked
indices on one worldline (`unit_not_selected`).

(5) Scope.  On the resting worldline with `u = τ` the clocked chain and
the index chain coincide at every index, and on every moving uniform
worldline at every timelike unit `c(n) < n` for every `n ≥ 1`
(`clocked_vs_index_scope`).

PRIOR WORK.  `Geometry/SourceClockRateAlongWorldlines.lean` proves the
exact difference of the two accrual rules, `properLength_generated_counts`,
`properLength_restWorldline`, `properLength_uniform`, `dilationFactor` and
its bounds, and the declared readings `indexRuleReturn` and
`lengthRuleReturnIndexSteps` (real numbers); no integer clock and no
composition with a trajectory appear there.  `QFT/SourceRecurrenceClock.lean`
proves `kac_identity`, `meanReturn_values`, and
`returnCount_block_additive` on the abstract step index; the composition
`returnCount i ω (c n)` appears here for the first time.
`Geometry/ProperTimeInternalAction.lean` proves `properLength_refine`.
`Geometry/CommonWorldMaxwellClockJoin.lean` proves `joined_clock_advance`
and `joined_stepDuration_not_forced`.  No prior module defines a
floor-clocked index from a proper length, nor a clocked return count.

ROWS TOUCHED.  The source clock and duration row (`u`, `τ`, and the
clocking rule are declared); the coupled-action row (the clocked chain is
a declared enrichment); the physical spacetime attachment row (the
Lorentz module is the declared image of the port class); the laboratory
clock and energy calibration import (no unit is attached to `u` or `τ`);
the light-signal row and the gravitation-route energy identification are
named and untouched.  None discharged (`dischargedRows_empty`).

NEGATIVES CITED.  The Legendre non-identifiability at scope
(`legendre_scope_cited`, re-cited): realized histories select no
accrual shape, so the clocking rule is a declared enrichment.  The rate
non-identifiability at scope (`rate_scope_cited`, re-cited): the declared
reading of an abstract transition system fixes no first-locking count
and forbids nothing about the clocked chain.

CONVENTIONS.  Signature `(+---)`; the generated path advances `τ` per
step in the scalar coordinate; `properLength` sums the clamped roots of
the forward Lorentz squares; `⌊·⌋₊` is the natural floor (zero below
zero); a window of `n` index steps is `0, …, n - 1`; `returnCount`
counts visits at chain steps `1, …, m`; `clockedReturns` counts index
steps `1, …, n` at which the clock advances into a visit.

FALSIFIER.  The module is wrong if the clocked index on the resting
worldline at `u = τ` differs from the index at some `n`, if the clocked
index on a uniform worldline exceeds `n d τ / u` at some `n`, if the
clocked return count differs from the chain's return count at the
clocked index for some clock with unit increments, or if the clocked
index changes under one midpoint refinement at the doubled index.

Axiom audit.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`; no `sorry`, no
`native_decide`, no project axiom.
-/

namespace OPH.ProperLengthClockedChain

open OPH.C1Lorentz (Herm2 lorentzQ)
open OPH.ProperTimeInternalAction (properLength refine properLength_refine)
open OPH.WorldlineHopTransport
open OPH.SourceClockRateAlongWorldlines
open OPH.InternalEnergyInertia (OpenRow)
open OPH.QFT.SourceRecurrenceClock (returnCount returnCount_succ returnCount_block_additive
  meanReturn meanReturn_values meanReturn_pos chainPi chainPi_pos kac_identity)

noncomputable section

/-! ## (1) The declared clocking rule -/

/-- The clocked index of a `Herm2` path at a declared internal unit `u`:
the number of complete units of `u` contained in the proper length over
the first `n` steps.  The rule is declared. -/
def clockedIndexPath (u : ℝ) (x : ℕ → Herm2) (n : ℕ) : ℕ := ⌊properLength n x / u⌋₊

/-- The accumulated proper length of a seam-step worldline at unit `τ`. -/
def accumulatedLength (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) : ℝ :=
  properLength n (generatedPath τ w)

/-- The clocked index along a seam-step worldline (declared): the internal
process advances one chain step each time the accumulated proper length
crosses a multiple of `u`. -/
def clockedIndex (τ u : ℝ) (w : SeamStepWorldline) (n : ℕ) : ℕ :=
  ⌊accumulatedLength τ w n / u⌋₊

theorem clockedIndex_eq_path (τ u : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    clockedIndex τ u w n = clockedIndexPath u (generatedPath τ w) n := rfl

theorem properLength_succ (x : ℕ → Herm2) (n : ℕ) :
    properLength (n + 1) x = properLength n x + Real.sqrt (lorentzQ (x (n + 1) - x n)) := by
  unfold properLength
  rw [Finset.sum_range_succ]

theorem properLength_nonneg (x : ℕ → Herm2) (n : ℕ) : 0 ≤ properLength n x :=
  Finset.sum_nonneg fun _ _ => Real.sqrt_nonneg _

theorem properLength_mono (x : ℕ → Herm2) (n : ℕ) :
    properLength n x ≤ properLength (n + 1) x := by
  rw [properLength_succ]
  linarith [Real.sqrt_nonneg (lorentzQ (x (n + 1) - x n))]

theorem accumulatedLength_nonneg (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    0 ≤ accumulatedLength τ w n :=
  properLength_nonneg _ n

theorem accumulatedLength_mono (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    accumulatedLength τ w n ≤ accumulatedLength τ w (n + 1) :=
  properLength_mono _ n

theorem clockedIndexPath_mono (u : ℝ) (hu : 0 < u) (x : ℕ → Herm2) (n : ℕ) :
    clockedIndexPath u x n ≤ clockedIndexPath u x (n + 1) :=
  Nat.floor_mono (div_le_div_of_nonneg_right (properLength_mono x n) hu.le)

/-- **Monotone.**  The clocked index never decreases. -/
theorem clockedIndex_mono (τ u : ℝ) (hu : 0 < u) (w : SeamStepWorldline) :
    Monotone (clockedIndex τ u w) :=
  monotone_nat_of_le_succ fun n => clockedIndexPath_mono u hu _ n

/-- **Unit increments.**  When the proper length of step `n` is at most
`u`, the clocked index advances by `0` or `1` at that step. -/
theorem clockedIndexPath_succ_of_step_le (u : ℝ) (hu : 0 < u) (x : ℕ → Herm2) (n : ℕ)
    (hstep : Real.sqrt (lorentzQ (x (n + 1) - x n)) ≤ u) :
    clockedIndexPath u x (n + 1) = clockedIndexPath u x n ∨
      clockedIndexPath u x (n + 1) = clockedIndexPath u x n + 1 := by
  have hle := clockedIndexPath_mono u hu x n
  have hup : clockedIndexPath u x (n + 1) ≤ clockedIndexPath u x n + 1 := by
    unfold clockedIndexPath
    rw [← Nat.floor_add_one (div_nonneg (properLength_nonneg x n) hu.le)]
    apply Nat.floor_mono
    rw [properLength_succ, add_div]
    have : Real.sqrt (lorentzQ (x (n + 1) - x n)) / u ≤ 1 := by
      rw [div_le_one hu]
      exact hstep
    linarith
  omega

theorem clockedIndex_succ_of_step_le (τ u : ℝ) (hu : 0 < u) (w : SeamStepWorldline) (n : ℕ)
    (hstep : Real.sqrt (lorentzQ (generatedPath τ w (n + 1) - generatedPath τ w n)) ≤ u) :
    clockedIndex τ u w (n + 1) = clockedIndex τ u w n ∨
      clockedIndex τ u w (n + 1) = clockedIndex τ u w n + 1 :=
  clockedIndexPath_succ_of_step_le u hu _ n hstep

theorem stepNormSq_nonneg (s : SeamStep) : 0 ≤ stepNormSq s := by
  cases s <;> simp [stepNormSq, seamNormSq]

/-- The proper length of one step of a generated path at unit `τ ≥ 0` is
at most `τ`. -/
theorem generated_step_le (τ : ℝ) (hτ : 0 ≤ τ) (w : SeamStepWorldline) (n : ℕ) :
    Real.sqrt (lorentzQ (generatedPath τ w (n + 1) - generatedPath τ w n)) ≤ τ := by
  rw [sqrt_lorentzQ_generated_step]
  calc Real.sqrt (τ ^ 2 - stepNormSq (w.steps n)) ≤ Real.sqrt (τ ^ 2) :=
        Real.sqrt_le_sqrt (by linarith [stepNormSq_nonneg (w.steps n)])
    _ = τ := Real.sqrt_sq hτ

/-- **Unit increments on every worldline** at every unit `0 ≤ τ ≤ u`. -/
theorem clockedIndex_unit_increments (τ u : ℝ) (hτ : 0 ≤ τ) (hu : τ ≤ u) (hu0 : 0 < u)
    (w : SeamStepWorldline) (n : ℕ) :
    clockedIndex τ u w (n + 1) = clockedIndex τ u w n ∨
      clockedIndex τ u w (n + 1) = clockedIndex τ u w n + 1 :=
  clockedIndex_succ_of_step_le τ u hu0 w n (le_trans (generated_step_le τ hτ w n) hu)

/-- **The index rule is the clocked rule at rest.**  On the resting
worldline with `u = τ > 0` the clocked index is the index, exactly. -/
theorem clockedIndex_rest (τ : ℝ) (hτ : 0 < τ) (p : Fin 12) (n : ℕ) :
    clockedIndex τ τ (restWorldline p) n = n := by
  unfold clockedIndex accumulatedLength
  rw [properLength_restWorldline τ hτ.le, mul_div_assoc, div_self hτ.ne', mul_one,
    Nat.floor_natCast]

/-! ## (2) The uniform worldline: bounds and the limit rate -/

/-- The crossing count of the uniform worldline over `(k + 1) p + j` steps,
`j ≤ k`: `p` full periods contribute `p`, a partial period contributes one
crossing when nonempty. -/
theorem crossingCount_uniform_general (e : Fin 30) (k p j : ℕ) (hj : j ≤ k) :
    crossingCount (uniformWorldline e k).steps ((k + 1) * p + j) =
      p + if j = 0 then 0 else 1 := by
  rw [crossingCount_block_additive, crossingCount_uniform]
  rcases j with _ | j
  · simp [crossingCount]
  · simp only [Nat.succ_ne_zero, if_false]
    congr 1
    exact crossingCount_uniform_block e k p j (by omega)

/-- The accumulated proper length of the uniform worldline at a general
index `(k + 1) p + j`, `j ≤ k`. -/
theorem accumulatedLength_uniform_general (τ : ℝ) (hτ : 0 ≤ τ) (e : Fin 30) (k p j : ℕ)
    (hj : j ≤ k) :
    accumulatedLength τ (uniformWorldline e k) ((k + 1) * p + j) =
      (p : ℝ) * ((k : ℝ) * τ + Real.sqrt (τ ^ 2 - 4)) +
        if j = 0 then 0 else Real.sqrt (τ ^ 2 - 4) + ((j : ℝ) - 1) * τ := by
  unfold accumulatedLength
  rw [properLength_generated_counts τ hτ]
  have hc := crossingCount_uniform_general e k p j hj
  have hr := crossingCount_add_restCount (uniformWorldline e k).steps ((k + 1) * p + j)
  rw [hc] at hr ⊢
  rcases j with _ | j
  · simp only [if_true, Nat.add_zero] at hr ⊢
    have hrest : restCount (uniformWorldline e k).steps ((k + 1) * p) = k * p := by
      have : (k + 1) * p = k * p + p := by ring
      omega
    rw [hrest]
    push_cast
    ring
  · simp only [Nat.succ_ne_zero, if_false] at hr ⊢
    have hrest : restCount (uniformWorldline e k).steps ((k + 1) * p + (j + 1)) =
        k * p + j := by
      have : (k + 1) * p = k * p + p := by ring
      omega
    rw [hrest]
    push_cast
    ring

theorem sqrt_sub_four_le (τ : ℝ) (hτ : 0 ≤ τ) : Real.sqrt (τ ^ 2 - 4) ≤ τ := by
  calc Real.sqrt (τ ^ 2 - 4) ≤ Real.sqrt (τ ^ 2) := Real.sqrt_le_sqrt (by linarith)
    _ = τ := Real.sqrt_sq hτ

theorem dilationFactor_mul_unit (τ : ℝ) (hτ : 0 < τ) (k : ℕ) :
    dilationFactor τ k * τ = ((k : ℝ) * τ + Real.sqrt (τ ^ 2 - 4)) / ((k : ℝ) + 1) := by
  unfold dilationFactor
  have hk : ((k : ℝ) + 1) ≠ 0 := by positivity
  field_simp

/-- The deficit of the accumulated length below `n d τ` at a general index:
zero on full periods, `(k + 1 - j) / (k + 1)` times `τ - √(τ² - 4)` inside
a period. -/
theorem uniform_deficit (τ : ℝ) (hτ : 0 < τ) (e : Fin 30) (k p j : ℕ) (hj : j ≤ k) :
    (((k + 1) * p + j : ℕ) : ℝ) * (dilationFactor τ k * τ) -
        accumulatedLength τ (uniformWorldline e k) ((k + 1) * p + j) =
      if j = 0 then 0 else
        (((k : ℝ) + 1 - j) / ((k : ℝ) + 1)) * (τ - Real.sqrt (τ ^ 2 - 4)) := by
  rw [accumulatedLength_uniform_general τ hτ.le e k p j hj, dilationFactor_mul_unit τ hτ]
  have hk : ((k : ℝ) + 1) ≠ 0 := by positivity
  rcases j with _ | j
  · simp only [if_true]
    push_cast
    field_simp
    ring
  · simp only [Nat.succ_ne_zero, if_false]
    push_cast
    field_simp
    ring

/-- **Upper bound.**  `L(n) ≤ n d τ` at every index. -/
theorem accumulatedLength_uniform_le (τ : ℝ) (hτ : 0 < τ) (e : Fin 30) (k n : ℕ) :
    accumulatedLength τ (uniformWorldline e k) n ≤ (n : ℝ) * (dilationFactor τ k * τ) := by
  have hn : n = (k + 1) * (n / (k + 1)) + n % (k + 1) := (Nat.div_add_mod n (k + 1)).symm
  have hj : n % (k + 1) ≤ k := Nat.lt_succ_iff.mp (Nat.mod_lt n (Nat.succ_pos k))
  have hd := uniform_deficit τ hτ e k (n / (k + 1)) (n % (k + 1)) hj
  rw [← hn] at hd
  have hS := sqrt_sub_four_le τ hτ.le
  have hnn : (0 : ℝ) ≤ ((k : ℝ) + 1 - (n % (k + 1) : ℕ)) / ((k : ℝ) + 1) := by
    apply div_nonneg _ (by positivity)
    have : ((n % (k + 1) : ℕ) : ℝ) ≤ k := by exact_mod_cast hj
    linarith
  split_ifs at hd with h0
  · linarith
  · nlinarith [mul_nonneg hnn (sub_nonneg.mpr hS)]

/-- **Lower bound.**  `n d τ - (τ - √(τ² - 4)) ≤ L(n)` at every index. -/
theorem accumulatedLength_uniform_ge (τ : ℝ) (hτ : 0 < τ) (e : Fin 30) (k n : ℕ) :
    (n : ℝ) * (dilationFactor τ k * τ) - (τ - Real.sqrt (τ ^ 2 - 4)) ≤
      accumulatedLength τ (uniformWorldline e k) n := by
  have hn : n = (k + 1) * (n / (k + 1)) + n % (k + 1) := (Nat.div_add_mod n (k + 1)).symm
  have hj : n % (k + 1) ≤ k := Nat.lt_succ_iff.mp (Nat.mod_lt n (Nat.succ_pos k))
  have hd := uniform_deficit τ hτ e k (n / (k + 1)) (n % (k + 1)) hj
  rw [← hn] at hd
  have hS := sqrt_sub_four_le τ hτ.le
  have hle1 : ((k : ℝ) + 1 - (n % (k + 1) : ℕ)) / ((k : ℝ) + 1) ≤ 1 := by
    rw [div_le_one (by positivity)]
    have : (0 : ℝ) ≤ ((n % (k + 1) : ℕ) : ℝ) := by positivity
    linarith
  split_ifs at hd with h0
  · linarith
  · nlinarith [mul_le_mul_of_nonneg_right hle1 (sub_nonneg.mpr hS)]

/-- On full periods the accumulated length is exactly `n d τ`. -/
theorem accumulatedLength_uniform_period (τ : ℝ) (hτ : 0 < τ) (e : Fin 30) (k p : ℕ) :
    accumulatedLength τ (uniformWorldline e k) ((k + 1) * p) =
      (((k + 1) * p : ℕ) : ℝ) * (dilationFactor τ k * τ) := by
  have hd := uniform_deficit τ hτ e k p 0 (Nat.zero_le k)
  simp only [if_true, Nat.add_zero] at hd
  linarith

/-- **Two-sided bound on the clocked index** at every index of the uniform
worldline: `n d τ / u - (τ - √(τ² - 4)) / u - 1 < c(n) ≤ n d τ / u`. -/
theorem clockedIndex_uniform_bounds (τ u : ℝ) (hτ : 0 < τ) (hu : 0 < u) (e : Fin 30)
    (k n : ℕ) :
    (n : ℝ) * (dilationFactor τ k * τ) / u - (τ - Real.sqrt (τ ^ 2 - 4)) / u - 1 <
        (clockedIndex τ u (uniformWorldline e k) n : ℝ) ∧
      (clockedIndex τ u (uniformWorldline e k) n : ℝ) ≤
        (n : ℝ) * (dilationFactor τ k * τ) / u := by
  unfold clockedIndex
  have hL0 := accumulatedLength_nonneg τ (uniformWorldline e k) n
  have hle := accumulatedLength_uniform_le τ hτ e k n
  have hge := accumulatedLength_uniform_ge τ hτ e k n
  constructor
  · have h1 := Nat.lt_floor_add_one (accumulatedLength τ (uniformWorldline e k) n / u)
    have h2 : ((n : ℝ) * (dilationFactor τ k * τ) - (τ - Real.sqrt (τ ^ 2 - 4))) / u ≤
        accumulatedLength τ (uniformWorldline e k) n / u :=
      div_le_div_of_nonneg_right hge hu.le
    rw [sub_div] at h2
    linarith
  · exact le_trans (Nat.floor_le (div_nonneg hL0 hu.le)) (div_le_div_of_nonneg_right hle hu.le)

/-- **Exact two-sided bound on full periods**: `n d τ / u - 1 < c(n) ≤ n d τ / u`
at `n = (k + 1) p`. -/
theorem clockedIndex_uniform_period_bounds (τ u : ℝ) (hτ : 0 < τ) (hu : 0 < u) (e : Fin 30)
    (k p : ℕ) :
    (((k + 1) * p : ℕ) : ℝ) * (dilationFactor τ k * τ) / u - 1 <
        (clockedIndex τ u (uniformWorldline e k) ((k + 1) * p) : ℝ) ∧
      (clockedIndex τ u (uniformWorldline e k) ((k + 1) * p) : ℝ) ≤
        (((k + 1) * p : ℕ) : ℝ) * (dilationFactor τ k * τ) / u := by
  unfold clockedIndex
  rw [accumulatedLength_uniform_period τ hτ e k p]
  have hnn : 0 ≤ (((k + 1) * p : ℕ) : ℝ) * (dilationFactor τ k * τ) / u := by
    rw [← accumulatedLength_uniform_period τ hτ e k p]
    exact div_nonneg (accumulatedLength_nonneg τ _ _) hu.le
  constructor
  · have h := Nat.lt_floor_add_one ((((k + 1) * p : ℕ) : ℝ) * (dilationFactor τ k * τ) / u)
    linarith
  · exact Nat.floor_le hnn

/-- **The clocked rate on the uniform worldline**: `c(n) / n → d τ / u`. -/
theorem clockedIndex_uniform_rate (τ u : ℝ) (hτ : 0 < τ) (hu : 0 < u) (e : Fin 30) (k : ℕ) :
    Filter.Tendsto (fun n : ℕ => (clockedIndex τ u (uniformWorldline e k) n : ℝ) / n)
      Filter.atTop (nhds (dilationFactor τ k * τ / u)) := by
  set r : ℝ := dilationFactor τ k * τ / u with hr
  set D : ℝ := (τ - Real.sqrt (τ ^ 2 - 4)) / u + 1 with hD
  have hlow : Filter.Tendsto (fun n : ℕ => r - D * (1 / (n : ℝ))) Filter.atTop (nhds r) := by
    have h := (tendsto_one_div_atTop_nhds_zero_nat (𝕜 := ℝ)).const_mul D
    have h2 := (tendsto_const_nhds (x := r)).sub h
    simpa using h2
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' hlow tendsto_const_nhds ?_ ?_
  · filter_upwards [Filter.eventually_ge_atTop 1] with n hn
    have hn' : (0 : ℝ) < n := by exact_mod_cast hn
    have hb := (clockedIndex_uniform_bounds τ u hτ hu e k n).1
    rw [le_div_iff₀ hn']
    have : (n : ℝ) * (dilationFactor τ k * τ) / u = n * r := by rw [hr]; ring
    rw [this] at hb
    have h2 : (r - D * (1 / (n : ℝ))) * n = r * n - D := by
      field_simp
    rw [h2, hD]
    linarith
  · filter_upwards [Filter.eventually_ge_atTop 1] with n hn
    have hn' : (0 : ℝ) < n := by exact_mod_cast hn
    have hb := (clockedIndex_uniform_bounds τ u hτ hu e k n).2
    rw [div_le_iff₀ hn']
    have : (n : ℝ) * (dilationFactor τ k * τ) / u = n * r := by rw [hr]; ring
    rw [this] at hb
    linarith

/-- **At `u = τ` the clocked rate is the dilation factor**, exactly. -/
theorem clockedIndex_uniform_rate_unit (τ : ℝ) (hτ : 0 < τ) (e : Fin 30) (k : ℕ) :
    Filter.Tendsto (fun n : ℕ => (clockedIndex τ τ (uniformWorldline e k) n : ℝ) / n)
      Filter.atTop (nhds (dilationFactor τ k)) := by
  have h := clockedIndex_uniform_rate τ τ hτ hτ e k
  rwa [mul_div_assoc, div_self hτ.ne', mul_one] at h

/-- **At rest the clocked rate is one.** -/
theorem clockedIndex_rest_rate (τ : ℝ) (hτ : 0 < τ) (p : Fin 12) :
    Filter.Tendsto (fun n : ℕ => (clockedIndex τ τ (restWorldline p) n : ℝ) / n)
      Filter.atTop (nhds 1) := by
  refine (tendsto_const_nhds (x := (1 : ℝ))).congr' ?_
  filter_upwards [Filter.eventually_ge_atTop 1] with n hn
  have hn' : (n : ℝ) ≠ 0 := by
    have : (0 : ℝ) < n := by exact_mod_cast hn
    exact this.ne'
  rw [clockedIndex_rest τ hτ p n, div_self hn']

/-! ## (3) The clocked chain -/

/-- The clocked trajectory (declared): the abstract chain state read at the
clocked index. -/
def clockedTrajectory (ω : ℕ → Fin 2) (c : ℕ → ℕ) : ℕ → Fin 2 := ω ∘ c

/-- The clocked return count: the number of index steps `j + 1 ≤ n` at
which the clock advances by one chain step and the chain enters state
`i` at that step. -/
def clockedReturns (i : Fin 2) (ω : ℕ → Fin 2) (c : ℕ → ℕ) : ℕ → ℕ
  | 0 => 0
  | n + 1 => clockedReturns i ω c n +
      if c (n + 1) = c n + 1 ∧ ω (c (n + 1)) = i then 1 else 0

/-- A clock with `c 0 = 0` and unit increments. -/
structure UnitClock (c : ℕ → ℕ) : Prop where
  zero : c 0 = 0
  step : ∀ n, c (n + 1) = c n ∨ c (n + 1) = c n + 1

theorem UnitClock.mono {c : ℕ → ℕ} (hc : UnitClock c) : Monotone c :=
  monotone_nat_of_le_succ fun n => by rcases hc.step n with h | h <;> omega

/-- **Returns are counted at chain steps.**  For a unit clock the clocked
return count over `n` index steps equals the chain's return count over
`c n` chain steps, exactly. -/
theorem clockedReturns_eq_returnCount (i : Fin 2) (ω : ℕ → Fin 2) (c : ℕ → ℕ)
    (hc : UnitClock c) (n : ℕ) :
    clockedReturns i ω c n = returnCount i ω (c n) := by
  induction n with
  | zero => simp [clockedReturns, hc.zero]
  | succ n ih =>
    simp only [clockedReturns]
    rcases hc.step n with h | h
    · rw [h, ih]
      simp
    · rw [h, returnCount_succ, ih]
      simp

/-- **Block additivity of the composed count** across index windows: the
count over `n + m` index steps is the count over the first `n` plus the
chain's count of the shifted trajectory over the `c (n + m) - c n` chain
steps completed in the second window. -/
theorem returnCount_clocked_block_additive (i : Fin 2) (ω : ℕ → Fin 2) (c : ℕ → ℕ)
    (hmono : Monotone c) (n m : ℕ) :
    returnCount i ω (c (n + m)) =
      returnCount i ω (c n) + returnCount i (fun j => ω (c n + j)) (c (n + m) - c n) := by
  have h : c (n + m) = c n + (c (n + m) - c n) :=
    (Nat.add_sub_cancel' (hmono (Nat.le_add_right n m))).symm
  conv_lhs => rw [h]
  exact returnCount_block_additive i ω (c n) (c (n + m) - c n)

theorem clockedReturns_block_additive (i : Fin 2) (ω : ℕ → Fin 2) (c : ℕ → ℕ)
    (hc : UnitClock c) (n m : ℕ) :
    clockedReturns i ω c (n + m) =
      clockedReturns i ω c n + returnCount i (fun j => ω (c n + j)) (c (n + m) - c n) := by
  rw [clockedReturns_eq_returnCount i ω c hc, clockedReturns_eq_returnCount i ω c hc,
    returnCount_clocked_block_additive i ω c hc.mono]

theorem clockedIndex_zero (τ u : ℝ) (w : SeamStepWorldline) : clockedIndex τ u w 0 = 0 := by
  simp [clockedIndex, accumulatedLength, properLength]

/-- The clocked index along a worldline is a unit clock at every unit
`0 ≤ τ ≤ u`. -/
theorem clockedIndex_unitClock (τ u : ℝ) (hτ : 0 ≤ τ) (hu : τ ≤ u) (hu0 : 0 < u)
    (w : SeamStepWorldline) : UnitClock (clockedIndex τ u w) :=
  ⟨clockedIndex_zero τ u w, clockedIndex_unit_increments τ u hτ hu hu0 w⟩

/-- The clocked return count along a worldline (declared): the chain's
return count at the clocked index. -/
def worldlineReturns (i : Fin 2) (ω : ℕ → Fin 2) (τ u : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    ℕ :=
  returnCount i ω (clockedIndex τ u w n)

theorem worldlineReturns_eq_clockedReturns (i : Fin 2) (ω : ℕ → Fin 2) (τ u : ℝ) (hτ : 0 ≤ τ)
    (hu : τ ≤ u) (hu0 : 0 < u) (w : SeamStepWorldline) (n : ℕ) :
    worldlineReturns i ω τ u w n = clockedReturns i ω (clockedIndex τ u w) n :=
  (clockedReturns_eq_returnCount i ω _ (clockedIndex_unitClock τ u hτ hu hu0 w) n).symm

/-- **Generic conversion.**  Any quantity with a constant rate `a > 0` per
chain step, evaluated at a chain-step count `c` with `A - 1 < c ≤ A`,
lies in `(a (A - 1), a A]`. -/
theorem perChainStep_conversion (a A : ℝ) (ha : 0 < a) (c : ℕ) (h1 : A - 1 < c)
    (h2 : (c : ℝ) ≤ A) :
    a * (A - 1) < a * c ∧ a * c ≤ a * A :=
  ⟨mul_lt_mul_of_pos_left h1 ha, mul_le_mul_of_nonneg_left h2 ha.le⟩

/-- The expected number of returns to the committed stable class over `n`
index steps (declared reading): the stationary mass `chainPi 0` per chain
step, `1 / meanReturn 0` by `kac_identity`, times the clocked index.  The
reading of the stationary mass as an expected return rate is declared;
no probability measure on trajectories is constructed here. -/
def expectedReturns (τ u : ℝ) (w : SeamStepWorldline) (n : ℕ) : ℝ :=
  chainPi 0 * clockedIndex τ u w n

theorem chainPi_eq_inv_meanReturn : chainPi 0 = 1 / meanReturn 0 := by
  rw [kac_identity, one_div_one_div]

/-- **Expected returns on full periods**, in `(chainPi 0 (n d τ / u - 1),
chainPi 0 · n d τ / u]`. -/
theorem expectedReturns_uniform_period_bounds (τ u : ℝ) (hτ : 0 < τ) (hu : 0 < u) (e : Fin 30)
    (k p : ℕ) :
    chainPi 0 * ((((k + 1) * p : ℕ) : ℝ) * (dilationFactor τ k * τ) / u - 1) <
        expectedReturns τ u (uniformWorldline e k) ((k + 1) * p) ∧
      expectedReturns τ u (uniformWorldline e k) ((k + 1) * p) ≤
        chainPi 0 * ((((k + 1) * p : ℕ) : ℝ) * (dilationFactor τ k * τ) / u) :=
  perChainStep_conversion (chainPi 0) _ (chainPi_pos 0) _
    (clockedIndex_uniform_period_bounds τ u hτ hu e k p).1
    (clockedIndex_uniform_period_bounds τ u hτ hu e k p).2

/-- **Expected returns per index step** tend to `(d τ / u) / meanReturn 0`. -/
theorem expectedReturns_uniform_rate (τ u : ℝ) (hτ : 0 < τ) (hu : 0 < u) (e : Fin 30)
    (k : ℕ) :
    Filter.Tendsto (fun n : ℕ => expectedReturns τ u (uniformWorldline e k) n / n)
      Filter.atTop (nhds ((dilationFactor τ k * τ / u) / meanReturn 0)) := by
  have h := (clockedIndex_uniform_rate τ u hτ hu e k).const_mul (chainPi 0)
  rw [chainPi_eq_inv_meanReturn] at h
  have hfun : (fun n : ℕ => expectedReturns τ u (uniformWorldline e k) n / n) =
      fun n : ℕ => 1 / meanReturn 0 * ((clockedIndex τ u (uniformWorldline e k) n : ℝ) / n) := by
    funext n
    unfold expectedReturns
    rw [chainPi_eq_inv_meanReturn]
    ring
  rw [hfun]
  have : (1 / meanReturn 0) * (dilationFactor τ k * τ / u) =
      (dilationFactor τ k * τ / u) / meanReturn 0 := by ring
  rwa [this] at h

/-- The mean return read in index steps on the uniform worldline (declared
reading): `meanReturn 0` chain steps are `meanReturn 0` units of proper
length `u`, and the uniform worldline accrues `d τ` per index step. -/
def returnIndexSteps (τ u : ℝ) (k : ℕ) : ℝ := meanReturn 0 * u / (dilationFactor τ k * τ)

theorem returnIndexSteps_value (τ u : ℝ) (k : ℕ) :
    returnIndexSteps τ u k = (61511 / 7155 : ℝ) * u / (dilationFactor τ k * τ) := by
  unfold returnIndexSteps
  rw [meanReturn_values.1]

/-- At `u = τ` the reading agrees with `lengthRuleReturnIndexSteps`. -/
theorem returnIndexSteps_unit (τ : ℝ) (hτ : 0 < τ) (k : ℕ) :
    returnIndexSteps τ τ k = lengthRuleReturnIndexSteps τ k := by
  unfold returnIndexSteps lengthRuleReturnIndexSteps
  rw [mul_comm (dilationFactor τ k) τ, ← div_div, mul_div_assoc, div_self hτ.ne', mul_one]

/-- The limit rate of expected returns per index step is the reciprocal of
`returnIndexSteps`. -/
theorem expectedReturns_rate_eq_inv (τ u : ℝ) (hτ : 2 < τ) (hu : 0 < u) (k : ℕ) :
    (dilationFactor τ k * τ / u) / meanReturn 0 = 1 / returnIndexSteps τ u k := by
  unfold returnIndexSteps
  have h1 := (dilationFactor_pos τ hτ k).ne'
  have h2 := (meanReturn_pos 0).ne'
  have h3 : τ ≠ 0 := by linarith
  field_simp

/-! ## (4) Refinement invariance, the join, and the unit -/

/-- The clocked index depends on the path only through the proper length. -/
theorem clockedIndexPath_congr (u : ℝ) (x y : ℕ → Herm2) (n : ℕ)
    (h : properLength n x = properLength n y) :
    clockedIndexPath u x n = clockedIndexPath u y n := by
  unfold clockedIndexPath
  rw [h]

/-- **Refinement invariance.**  The clocked index of the refined path at
the doubled index equals the clocked index of the path
(`properLength_refine`). -/
theorem clockedIndexPath_refine (u : ℝ) (x : ℕ → Herm2) (n : ℕ) :
    clockedIndexPath u (refine x) (2 * n) = clockedIndexPath u x n := by
  unfold clockedIndexPath
  rw [properLength_refine]

theorem clockedIndex_refine (τ u : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    clockedIndexPath u (refine (generatedPath τ w)) (2 * n) = clockedIndex τ u w n :=
  clockedIndexPath_refine u _ n

/-- The clocked return count is refinement invariant at the doubled index. -/
theorem clockedReturns_refine (i : Fin 2) (ω : ℕ → Fin 2) (τ u : ℝ) (w : SeamStepWorldline)
    (n : ℕ) :
    returnCount i ω (clockedIndexPath u (refine (generatedPath τ w)) (2 * n)) =
      worldlineReturns i ω τ u w n := by
  unfold worldlineReturns
  rw [clockedIndex_refine]

theorem uniformSteps_zero_ne_rest (e : Fin 30) (k : ℕ) : uniformSteps e k 0 ≠ .rest := by
  rw [uniformSteps_ne_rest_iff]
  simp

theorem uniformSteps_one_rest (e : Fin 30) (k : ℕ) (hk : 1 ≤ k) : uniformSteps e k 1 = .rest := by
  by_contra h
  have h' := (uniformSteps_ne_rest_iff e k 1).mp h
  rw [Nat.mod_eq_of_lt (by omega)] at h'
  omega

/-- **The join's per-index clock does not reproduce the accumulated proper
length** of a uniform worldline with a rest between crossings: the join
advances its phase through `stepTime δ n = n δ` (`joined_clock_advance`),
and no `δ` matches `L(1) = √(τ² - 4)` and `L(2) = √(τ² - 4) + τ` at a
timelike unit.  The field of the join that the clocked chain would
replace is the proper-time argument `stepTime J.stepDuration n` of the
phase, by `accumulatedLength τ w n`; this is a declared enrichment. -/
theorem no_stepDuration_matches_uniform (τ : ℝ) (hτ : 2 < τ) (e : Fin 30) (k : ℕ) (hk : 1 ≤ k) :
    ¬ ∃ δ : ℝ, ∀ n : ℕ,
      OPH.CommonWorldInstrumentJoin.stepTime δ n = accumulatedLength τ (uniformWorldline e k) n := by
  rintro ⟨δ, hδ⟩
  have h1 := hδ 1
  have h2 := hδ 2
  unfold OPH.CommonWorldInstrumentJoin.stepTime accumulatedLength at h1 h2
  rw [properLength_generated_counts τ (by linarith)] at h1 h2
  have hs0 : (uniformWorldline e k).steps 0 ≠ .rest := uniformSteps_zero_ne_rest e k
  have hs1 : (uniformWorldline e k).steps 1 = .rest := uniformSteps_one_rest e k hk
  simp [crossingCount, restCount, hs0, hs1] at h1 h2
  have := sqrt_sub_four_lt τ hτ
  linarith

/-- At `k = 0` (a crossing at every step) the constant `δ = √(τ² - 4)`
reproduces the accumulated proper length at every index. -/
theorem stepDuration_matches_uniform_zero (τ : ℝ) (hτ : 0 ≤ τ) (e : Fin 30) :
    ∃ δ : ℝ, ∀ n : ℕ,
      OPH.CommonWorldInstrumentJoin.stepTime δ n = accumulatedLength τ (uniformWorldline e 0) n := by
  refine ⟨Real.sqrt (τ ^ 2 - 4), fun n => ?_⟩
  have h := accumulatedLength_uniform_general τ hτ e 0 n 0 le_rfl
  simp only [Nat.one_mul, if_true, Nat.cast_zero, zero_mul, zero_add, add_zero] at h
  rw [h]
  unfold OPH.CommonWorldInstrumentJoin.stepTime
  ring

/-- **The unit is not selected.**  Two units give two clocked indices on
one worldline: at `u = τ` the resting worldline reads `1` at index one,
at `u = 2 τ` it reads `0`. -/
theorem unit_not_selected (τ : ℝ) (hτ : 0 < τ) (p : Fin 12) :
    ∃ u u' : ℝ, 0 < u ∧ 0 < u' ∧ u ≠ u' ∧
      clockedIndex τ u (restWorldline p) 1 ≠ clockedIndex τ u' (restWorldline p) 1 := by
  refine ⟨τ, 2 * τ, hτ, by positivity, by linarith, ?_⟩
  rw [clockedIndex_rest τ hτ p 1]
  unfold clockedIndex accumulatedLength
  rw [properLength_restWorldline τ hτ.le]
  have : ((1 : ℕ) : ℝ) * τ / (2 * τ) = 1 / 2 := by
    field_simp
    ring
  rw [this, Nat.floor_eq_zero.mpr (by norm_num)]
  exact one_ne_zero

/-! ## (5) Scope -/

/-- **Scope.**  On the resting worldline with `u = τ` the clocked chain and
the index chain coincide at every index; on every uniform worldline at
every timelike unit the clocked index is strictly below the index at
every `n ≥ 1`. -/
theorem clocked_vs_index_scope (τ : ℝ) (hτ : 2 < τ) (p : Fin 12) (e : Fin 30) (k : ℕ) :
    (∀ n, clockedIndex τ τ (restWorldline p) n = n) ∧
      ∀ n, 1 ≤ n → clockedIndex τ τ (uniformWorldline e k) n < n := by
  have hτ0 : 0 < τ := by linarith
  refine ⟨clockedIndex_rest τ hτ0 p, fun n hn => ?_⟩
  have hb := (clockedIndex_uniform_bounds τ τ hτ0 hτ0 e k n).2
  rw [mul_div_assoc, mul_div_assoc, div_self hτ0.ne', mul_one] at hb
  have hd := dilationFactor_lt_one τ hτ k
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have : (clockedIndex τ τ (uniformWorldline e k) n : ℝ) < n := by nlinarith
  exact_mod_cast this

/-- The clocked trajectory on the resting worldline at `u = τ` is the chain
trajectory itself. -/
theorem clockedTrajectory_rest (τ : ℝ) (hτ : 0 < τ) (p : Fin 12) (ω : ℕ → Fin 2) :
    clockedTrajectory ω (clockedIndex τ τ (restWorldline p)) = ω := by
  funext n
  simp [clockedTrajectory, clockedIndex_rest τ hτ p n]

/-! ## Negatives cited and rows touched -/

/-- The Legendre non-identifiability at scope, re-cited: realized histories
select no accrual shape, so the clocking rule is a declared enrichment. -/
theorem legendre_scope_recited :
    (¬ ∃ vel : ℝ → ℝ → ℝ, OPH.Variational.SolvesMomentum
        OPH.Variational.chainLogLagrangian vel) ∧
      OPH.Variational.chainCurvedLagrangian 1 ≠ OPH.Variational.chainCurvedLagrangian 2 :=
  legendre_scope_cited

/-- The rate non-identifiability at scope, re-cited: the declared reading of
an abstract transition system fixes no first-locking count and forbids
nothing about the clocked chain. -/
theorem rate_scope_recited {S : Type} [DecidableEq S] (T : S → S) (x : S) (n : ℕ)
    (hlock : OPH.RateNonidentifiability.FirstLock T n x) (hpos : 0 < n)
    (rate : OPH.RateNonidentifiability.DeclaredReading S → ℕ)
    (hrate : ∀ s : ℕ, OPH.RateNonidentifiability.FirstLock
      (OPH.RateNonidentifiability.stutterStep T s)
      (rate (OPH.RateNonidentifiability.projReading T s)) (x, Fin.last s)) :
    False :=
  rate_scope_cited T x n hlock hpos rate hrate

/-- The rows this module touches, as register labels: the source clock and
duration row (`u`, `τ`, and the clocking rule declared), the coupled
action (the clocked chain is a declared enrichment), and the physical
spacetime attachment row.  The laboratory clock and energy calibration
import, the light-signal row, and the gravitation-route energy
identification carry no `OpenRow` label.  A label is not a discharge. -/
def touchedRows : List OpenRow :=
  [OpenRow.sourceClock, OpenRow.coupledAction, OpenRow.spacetimeAttachment]

/-- The rows this module discharges: none. -/
def dischargedRows : List OpenRow := []

theorem dischargedRows_empty : dischargedRows = [] := rfl

end

/-! ## Axiom receipts -/

#print axioms clockedIndex_mono
#print axioms clockedIndex_unit_increments
#print axioms clockedIndex_rest
#print axioms clockedIndex_uniform_bounds
#print axioms clockedIndex_uniform_period_bounds
#print axioms clockedIndex_uniform_rate
#print axioms clockedIndex_uniform_rate_unit
#print axioms clockedIndex_rest_rate
#print axioms clockedReturns_eq_returnCount
#print axioms returnCount_clocked_block_additive
#print axioms clockedReturns_block_additive
#print axioms perChainStep_conversion
#print axioms expectedReturns_uniform_period_bounds
#print axioms expectedReturns_uniform_rate
#print axioms returnIndexSteps_value
#print axioms returnIndexSteps_unit
#print axioms expectedReturns_rate_eq_inv
#print axioms clockedIndexPath_refine
#print axioms clockedReturns_refine
#print axioms no_stepDuration_matches_uniform
#print axioms stepDuration_matches_uniform_zero
#print axioms unit_not_selected
#print axioms clocked_vs_index_scope
#print axioms clockedTrajectory_rest
#print axioms dischargedRows_empty

end OPH.ProperLengthClockedChain
