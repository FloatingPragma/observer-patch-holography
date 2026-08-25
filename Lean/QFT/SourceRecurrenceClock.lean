import Mathlib
import MixingChainRealization
import PhysicalCalibrationImport
import ObserverPatchHolography.RateBridgeObstruction

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# A source-selected recurrence clock precursor

ROW TOUCHED.  Every statement in this module touches the source
clock/duration row (issue 736) and no other of the three open rows.  It
touches neither the physical spacetime attachment row nor the coupled
action row, and it discharges none of the three: the duration of one
chain step stays a declared tick, exactly as in
`PhysicalCalibrationImport.labFrequency_not_forced`.

WHAT IS PROVED.

(a) Delimitation on the repair layer.  The accepted-repair relation
`OPH.acceptedStep` admits no cycle (`acceptedStep_no_cycle`): every
accepted step strictly lowers `OPH.mismatchCount`, the fact
`RateBridgeObstruction.acceptedStep_changes_obs` rests on, so a
transitive chain of accepted steps never returns to its start.  The
choice-canonical iteration `RateBridgeObstruction.stepOnce` admits no
periodic orbit of positive period other than a fixed point
(`periodic_iterate_is_fixed`), and its iterate sequence is eventually
constant (`iterate_eventually_constant`).  A periodic clock cannot be
hosted by accepted repair steps.  This locates where a clock cannot come
from; it says nothing against a clock hosted by the source chain.

(b) Construction on the committed stationary source chain
`OPH.Thermodynamics.mixingChain` with stationary law
`mixingChainStationary`.  The first-return law of each state is written
in closed form from the kernel (`returnMass`): return at step one with
the stay weight, return at step `n + 2` after an excursion of `n` stays
in the other state.  It is a probability law (`returnMass_hasSum_one`),
its mean is finite and exact (`meanReturn_hasSum`), and Kac's identity
holds for the committed chain with the exact rationals
`61511/7155` and `61511/54356` (`kac_identity`, `meanReturn_values`).
The two mean return times have the source-determined ratio
`pi_1 / pi_0 = 54356/7155` (`meanReturn_ratio`), and this ratio is
invariant under any common rescaling of the declared tick
(`meanReturn_ratio_tick_free`): rescaling invariance means that for every
pair of calibrations the ratio of laboratory return durations is the same
number, so the ratio carries no tick.

(c) Boundaries.  (i) The return process is recurrent with a non-degenerate
return-time law: return times `1` and `2` both have positive probability
(`returnMass_one_pos`, `returnMass_two_pos`), no single period exists
(`no_single_period`), and no return time has probability one
(`returnMass_lt_one`).  (ii) One chain step stays a declared tick: two
ticks give two laboratory durations for the same mean return
(`returnDuration_not_forced`, `return_tick_declaration_not_forced`).
(iii) A return count along a realized trajectory is an integer-valued
count of distinct visit events: it is bounded by the number of steps,
grows by at most one per step, and is additive over adjacent disjoint
time blocks (`returnCount_block_additive`).  This is the property the
overlapping-window incidence statistic of `QFT/SourceClockCandidate.lean`
lacks, where consecutive windows share an edge, so the total `197` counts edge
incidences with each shared edge counted in two windows and is no count
of distinct events.

(d) Scope.  Kac's identity is proved for the committed two-state chain by
direct geometric-series computation over the reals with exact rational
values; no general Kac theorem is claimed.  `returnMass` is defined as the
excursion product of the kernel with no process object constructed, so its
reading as a first-hitting law is definitional and Kac's identity is a
closed-form series identity on that definition.

WHAT IS NOT PROVED.  No physical duration is identified: `returnMass`
lives on the abstract step index, and every laboratory reading is
conditional on a declared `ClockCalibration`.  The two standing
negatives are cited at their exact scope.
`RateNonidentifiability` is an abstract-transition-system result whose
bridge to the repair layer fails (`RateBridgeObstruction`); it forbids no
repair-layer clock and is not used here to forbid one; part (a) uses only
the concrete descent fact.  `RealizedHistoryLegendreNoGo` concerns the
velocity curvature of a Lagrangian shape and is untouched: no Lagrangian
enrichment appears in this module.  No frozen prediction, no premise
discharge, no simulation run.  The construction is a precursor on the
time row: a source-selected recurrent process with a tick-free ratio.
The source clock/duration row stays an import.
-/

namespace OPH.QFT.SourceRecurrenceClock

open OPH.Thermodynamics
open OPH.PhysicalCalibrationImport
open OPH OPH.AbstractRewriting OPH.RateNonidentifiability OPH.RateBridgeObstruction
open Function Relation

noncomputable section

/-! ## (a) Delimitation: accepted repair hosts no periodic clock -/

section RepairLayer

variable {C : OPHCarrier}

/-- A transitive chain of accepted repair steps strictly lowers the
mismatch count, by the same descent fact
`RateBridgeObstruction.acceptedStep_changes_obs` rests on. -/
theorem mismatchCount_lt_of_transGen {x y : Records C}
    (h : TransGen (acceptedStep C) x y) :
    mismatchCount y < mismatchCount x := by
  induction h with
  | single hs =>
      obtain ⟨i, rfl, hf⟩ := hs
      exact mismatchCount_localRepair_lt C i _ hf
  | tail _ hs ih =>
      obtain ⟨i, rfl, hf⟩ := hs
      exact lt_trans (mismatchCount_localRepair_lt C i _ hf) ih

/-- **No cycle.**  The accepted-repair relation admits no transitive
chain from a record back to itself. -/
theorem acceptedStep_no_cycle (x : Records C) :
    ¬ TransGen (acceptedStep C) x x :=
  fun h => lt_irrefl _ (mismatchCount_lt_of_transGen h)

/-- One `stepOnce` move never raises the mismatch count. -/
theorem mismatchCount_stepOnce_le (y : Records C) :
    mismatchCount (stepOnce C y) ≤ mismatchCount y := by
  by_cases h : Locked (stepOnce C) y
  · rw [show stepOnce C y = y from h]
  · exact le_of_lt (mismatchCount_stepOnce_lt C h)

/-- Iterates of `stepOnce` never raise the mismatch count. -/
theorem mismatchCount_iterate_le (y : Records C) (k : ℕ) :
    mismatchCount ((stepOnce C)^[k] y) ≤ mismatchCount y := by
  induction k with
  | zero => simp
  | succ k ih =>
      rw [Function.iterate_succ_apply']
      exact le_trans (mismatchCount_stepOnce_le _) ih

/-- **No periodic orbit.**  A record returned to by a positive number of
`stepOnce` iterations is locked: a periodic point of the repair
iteration is a fixed point. -/
theorem locked_of_iterate_eq_self (x : Records C) {d : ℕ} (hd : 0 < d)
    (h : (stepOnce C)^[d] x = x) : Locked (stepOnce C) x := by
  by_contra hx
  obtain ⟨k, rfl⟩ : ∃ k, d = k + 1 := ⟨d - 1, by omega⟩
  have h1 : mismatchCount ((stepOnce C)^[k + 1] x) ≤ mismatchCount (stepOnce C x) := by
    rw [Function.iterate_succ_apply]
    exact mismatchCount_iterate_le _ _
  have h2 := mismatchCount_stepOnce_lt C hx
  rw [h] at h1
  omega

/-- A periodic point of the repair iteration is stationary: every iterate
equals the point itself.  No orbit of positive period other than a fixed
point exists, so the repair iteration hosts no periodic clock. -/
theorem periodic_iterate_is_fixed (x : Records C) {d : ℕ} (hd : 0 < d)
    (h : (stepOnce C)^[d] x = x) : ∀ n : ℕ, (stepOnce C)^[n] x = x :=
  fun n => (locked_of_iterate_eq_self x hd h).iterate n

/-- Some iterate within `mismatchCount x` steps is locked, by
`RateBridgeObstruction.mismatchCount_iterate_add_le`. -/
theorem exists_locked_iterate (x : Records C) :
    ∃ j : ℕ, j ≤ mismatchCount x ∧ Locked (stepOnce C) ((stepOnce C)^[j] x) := by
  by_contra hcon
  have hcon' : ∀ j, j ≤ mismatchCount x → ¬ Locked (stepOnce C) ((stepOnce C)^[j] x) :=
    fun j hj hl => hcon ⟨j, hj, hl⟩
  have := mismatchCount_iterate_add_le C x (mismatchCount x + 1)
    (fun j hj => hcon' j (by omega))
  omega

/-- **Eventually constant.**  The iterate sequence of the repair
iteration is constant from some index on. -/
theorem iterate_eventually_constant (x : Records C) :
    ∃ N : ℕ, ∀ n : ℕ, N ≤ n → (stepOnce C)^[n] x = (stepOnce C)^[N] x := by
  obtain ⟨j, _, hj⟩ := exists_locked_iterate x
  refine ⟨j, fun n hn => ?_⟩
  obtain ⟨t, rfl⟩ : ∃ t, n = t + j := ⟨n - j, by omega⟩
  rw [Function.iterate_add_apply]
  exact hj.iterate t

end RepairLayer

/-! ## (b) Construction: the first-return law of the committed source chain -/

section SourceChain

/-- The committed kernel `mixingChain` over the reals. -/
def chainP (i j : Fin 2) : ℝ := ((mixingChain i j : ℚ) : ℝ)

/-- The committed stationary law `mixingChainStationary` over the reals. -/
def chainPi (i : Fin 2) : ℝ := ((mixingChainStationary i : ℚ) : ℝ)

/-- The other state of the two-state chain. -/
def other (i : Fin 2) : Fin 2 := ![1, 0] i

theorem chainP_pos (i j : Fin 2) : 0 < chainP i j := by
  simp only [chainP]
  exact_mod_cast mixingChain_entries_pos i j

theorem chainPi_pos (i : Fin 2) : 0 < chainPi i := by
  simp only [chainPi]
  exact_mod_cast mixingChainStationary_pos i

/-- The stay weight of either state is strictly below one. -/
theorem chainP_stay_lt_one (i : Fin 2) : chainP i i < 1 := by
  fin_cases i <;> norm_num [chainP, mixingChain]

/-- Row stochasticity in the two-state form. -/
theorem chainP_row (i : Fin 2) : chainP i i + chainP i (other i) = 1 := by
  fin_cases i <;> norm_num [chainP, other, mixingChain]

/-- Leaving the other state has weight one minus staying there. -/
theorem one_sub_stay (i : Fin 2) :
    1 - chainP (other i) (other i) = chainP (other i) i := by
  fin_cases i <;> norm_num [chainP, other, mixingChain]

/-- **The first-return law of state `i`**, in closed form from the kernel:
no return at step zero, return at step one with the stay weight, return at
step `n + 2` after leaving, staying `n` times in the other state, and
coming back.  This is the excursion product of the two-state chain; it is
a source-selected object because every factor is a committed kernel
entry.  No path-space measure and no Markov process object is built here:
`returnMass` is defined as the excursion product of the kernel, its reading
as a first-hitting law is definitional, and Kac's identity below is a
closed-form series identity on this definition. -/
def returnMass (i : Fin 2) : ℕ → ℝ
  | 0 => 0
  | 1 => chainP i i
  | n + 2 => chainP i (other i) * chainP (other i) (other i) ^ n * chainP (other i) i

@[simp] theorem returnMass_zero (i : Fin 2) : returnMass i 0 = 0 := rfl

@[simp] theorem returnMass_one (i : Fin 2) : returnMass i 1 = chainP i i := rfl

@[simp] theorem returnMass_add_two (i : Fin 2) (n : ℕ) :
    returnMass i (n + 2) =
      chainP i (other i) * chainP (other i) (other i) ^ n * chainP (other i) i := rfl

theorem returnMass_nonneg (i : Fin 2) (n : ℕ) : 0 ≤ returnMass i n := by
  match n with
  | 0 => simp
  | 1 => exact le_of_lt (chainP_pos _ _)
  | n + 2 =>
      rw [returnMass_add_two]
      exact mul_nonneg (mul_nonneg (le_of_lt (chainP_pos _ _))
        (pow_nonneg (le_of_lt (chainP_pos _ _)) n)) (le_of_lt (chainP_pos _ _))

/-- The excursion tail sums to the leave weight: a geometric series with
ratio the stay weight of the other state. -/
theorem returnMass_tail_hasSum (i : Fin 2) :
    HasSum (fun n : ℕ => returnMass i (n + 2)) (chainP i (other i)) := by
  have hr0 : 0 ≤ chainP (other i) (other i) := le_of_lt (chainP_pos _ _)
  have hr1 : chainP (other i) (other i) < 1 := chainP_stay_lt_one _
  have key : HasSum
      (fun n : ℕ => chainP i (other i) * chainP (other i) (other i) ^ n * chainP (other i) i)
      (chainP i (other i) * (1 - chainP (other i) (other i))⁻¹ * chainP (other i) i) :=
    ((hasSum_geometric_of_lt_one hr0 hr1).mul_left _).mul_right _
  rw [one_sub_stay, mul_assoc, inv_mul_cancel₀ (chainP_pos _ _).ne', mul_one] at key
  exact key

/-- **The first-return law is a probability law**: return is certain, so
the defined return law is proper from either state (part (d): the reading
as a first-hitting law is definitional). -/
theorem returnMass_hasSum_one (i : Fin 2) : HasSum (returnMass i) 1 := by
  have h := (hasSum_nat_add_iff (f := returnMass i) 2).1 (returnMass_tail_hasSum i)
  have hsum : chainP i (other i) + ∑ k ∈ Finset.range 2, returnMass i k = 1 := by
    simp [Finset.sum_range_succ]
    linarith [chainP_row i]
  rwa [hsum] at h

/-- The excursion part of the mean return time, by the two geometric
series `∑ n r^n = r / (1 - r)^2` and `∑ r^n = (1 - r)⁻¹`. -/
theorem meanReturn_tail_hasSum (i : Fin 2) :
    HasSum (fun n : ℕ => ((n : ℝ) + 2) * returnMass i (n + 2))
      (chainP i (other i) * (chainP (other i) (other i) /
          (1 - chainP (other i) (other i)) ^ 2) * chainP (other i) i
        + 2 * (chainP i (other i) * (1 - chainP (other i) (other i))⁻¹ * chainP (other i) i)) := by
  have hr0 : 0 ≤ chainP (other i) (other i) := le_of_lt (chainP_pos _ _)
  have hr1 : chainP (other i) (other i) < 1 := chainP_stay_lt_one _
  have hnorm : ‖chainP (other i) (other i)‖ < 1 := by
    rw [Real.norm_eq_abs, abs_lt]
    constructor <;> linarith
  have h1 := ((hasSum_coe_mul_geometric_of_norm_lt_one hnorm).mul_left
    (chainP i (other i))).mul_right (chainP (other i) i)
  have h2 := (((hasSum_geometric_of_lt_one hr0 hr1).mul_left
    (chainP i (other i))).mul_right (chainP (other i) i)).mul_left 2
  convert h1.add h2 using 1
  funext n
  rw [returnMass_add_two]
  ring

end SourceChain

/-! ## (b) Kac's identity for the committed chain -/

section Kac

/-- **Mean return time, summed.**  The series `∑ n * returnMass i n`
converges to `1 / pi_i`, with `pi` the committed stationary law: Kac's
identity for the committed two-state chain, checked on the exact
rationals of the kernel after the geometric-series computation. -/
theorem meanReturn_hasSum (i : Fin 2) :
    HasSum (fun n : ℕ => (n : ℝ) * returnMass i n) (1 / chainPi i) := by
  have htail : HasSum (fun n : ℕ => ((n + 2 : ℕ) : ℝ) * returnMass i (n + 2))
      (chainP i (other i) * (chainP (other i) (other i) /
          (1 - chainP (other i) (other i)) ^ 2) * chainP (other i) i
        + 2 * (chainP i (other i) * (1 - chainP (other i) (other i))⁻¹
            * chainP (other i) i)) := by
    have h := meanReturn_tail_hasSum i
    convert h using 1
    funext n
    push_cast
    ring
  have h := (hasSum_nat_add_iff (f := fun n : ℕ => (n : ℝ) * returnMass i n) 2).1 htail
  have hval : chainP i (other i) * (chainP (other i) (other i) /
          (1 - chainP (other i) (other i)) ^ 2) * chainP (other i) i
        + 2 * (chainP i (other i) * (1 - chainP (other i) (other i))⁻¹
            * chainP (other i) i)
        + ∑ k ∈ Finset.range 2, (k : ℝ) * returnMass i k = 1 / chainPi i := by
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, returnMass_zero,
      returnMass_one]
    fin_cases i <;>
      norm_num [chainP, chainPi, other, mixingChain, mixingChainStationary]
  rwa [hval] at h

/-- The mean return time of state `i`: the first moment of its
first-return law, on the abstract step index. -/
def meanReturn (i : Fin 2) : ℝ := ∑' n : ℕ, (n : ℝ) * returnMass i n

/-- **Kac's identity on the committed chain**: mean return time equals the
reciprocal stationary mass. -/
theorem kac_identity (i : Fin 2) : meanReturn i = 1 / chainPi i :=
  (meanReturn_hasSum i).tsum_eq

/-- The exact rational values: `61511/7155` steps for the committed stable
class (index `0`, checkpoint class `3`) and `61511/54356` steps for the
committed unstable class (index `1`, checkpoint class `2`); the index-to-class
map is `protectedRecordLabel = ![3, 2]` of
`Thermodynamics/MixingChainRealization.lean`. -/
theorem meanReturn_values :
    meanReturn 0 = 61511 / 7155 ∧ meanReturn 1 = 61511 / 54356 := by
  constructor <;>
    · rw [kac_identity]
      norm_num [chainPi, mixingChainStationary]

theorem meanReturn_pos (i : Fin 2) : 0 < meanReturn i := by
  rw [kac_identity]
  exact one_div_pos.2 (chainPi_pos i)

/-- **The source-determined ratio.**  The mean return times stand in the
ratio `pi_1 / pi_0 = 54356/7155`: a dimensionless number fixed by the
committed stationary law alone.  The ratio is exact on the committed
literals of one pinned run (row counts `1431` and `508` of the committed
kernel); no sampling uncertainty is modeled. -/
theorem meanReturn_ratio :
    meanReturn 0 / meanReturn 1 = chainPi 1 / chainPi 0 ∧
      meanReturn 0 / meanReturn 1 = 54356 / 7155 := by
  obtain ⟨h0, h1⟩ := meanReturn_values
  rw [h0, h1]
  constructor <;> norm_num [chainPi, mixingChainStationary]

/-- Under any declared tick the ratio of laboratory return durations is
the same source-determined number. -/
theorem labReturn_ratio (cal : ClockCalibration) :
    cal.labSeconds (meanReturn 0) / cal.labSeconds (meanReturn 1)
      = chainPi 1 / chainPi 0 := by
  simp only [ClockCalibration.labSeconds_def]
  rw [mul_div_mul_right _ _ cal.tau_pos.ne']
  exact meanReturn_ratio.1

/-- **Rescaling invariance of the ratio.**  For every pair of declared
ticks the ratio of laboratory return durations coincides: the ratio
carries no tick.  Rescaling invariance means exactly this and nothing
more; each individual duration keeps its tick. -/
theorem meanReturn_ratio_tick_free (cal cal' : ClockCalibration) :
    cal.labSeconds (meanReturn 0) / cal.labSeconds (meanReturn 1)
      = cal'.labSeconds (meanReturn 0) / cal'.labSeconds (meanReturn 1) := by
  rw [labReturn_ratio, labReturn_ratio]

end Kac

/-! ## (c) Boundaries -/

section Boundaries

/-- (i) Return at step one has positive probability. -/
theorem returnMass_one_pos (i : Fin 2) : 0 < returnMass i 1 := by
  rw [returnMass_one]
  exact chainP_pos i i

/-- (i) Return at step two has positive probability. -/
theorem returnMass_two_pos (i : Fin 2) : 0 < returnMass i 2 := by
  have h := returnMass_add_two i 0
  rw [pow_zero, mul_one] at h
  rw [h]
  exact mul_pos (chainP_pos _ _) (chainP_pos _ _)

/-- (i) **Recurrent with a non-degenerate return-time law.**  No return
time carries the whole law: for every candidate period `d` some other return time has positive
probability, so no single period exists. -/
theorem no_single_period (i : Fin 2) :
    ¬ ∃ d : ℕ, ∀ n : ℕ, n ≠ d → returnMass i n = 0 := by
  rintro ⟨d, hd⟩
  by_cases h1 : d = 1
  · subst h1
    exact (returnMass_two_pos i).ne' (hd 2 (by norm_num))
  · exact (returnMass_one_pos i).ne' (hd 1 (Ne.symm h1))

/-- (i) No return time has probability one. -/
theorem returnMass_lt_one (i : Fin 2) (d : ℕ) : returnMass i d < 1 := by
  by_cases h1 : d = 1
  · subst h1
    exact lt_hasSum (returnMass_hasSum_one i) 1 (fun j _ => returnMass_nonneg i j) 2
      (by norm_num) (returnMass_two_pos i)
  · exact lt_hasSum (returnMass_hasSum_one i) d (fun j _ => returnMass_nonneg i j) 1
      (Ne.symm h1) (returnMass_one_pos i)

/-- (ii) **One chain step stays a declared tick.**  Two distinct declared
ticks give two distinct laboratory durations for one chain step; nothing
in the chain selects between them. -/
theorem stepDuration_not_forced (cal cal' : ClockCalibration)
    (hne : cal.tau ≠ cal'.tau) :
    cal.labSeconds 1 ≠ cal'.labSeconds 1 := by
  simp only [ClockCalibration.labSeconds_def, one_mul]
  exact hne

/-- (ii) Mirror of `labFrequency_not_forced` for the mean return: two
distinct declared ticks give two distinct laboratory return durations. -/
theorem returnDuration_not_forced (cal cal' : ClockCalibration)
    (hne : cal.tau ≠ cal'.tau) (i : Fin 2) :
    cal.labSeconds (meanReturn i) ≠ cal'.labSeconds (meanReturn i) := by
  simp only [ClockCalibration.labSeconds_def]
  intro h
  exact hne (mul_left_cancel₀ (meanReturn_pos i).ne' h)

/-- (ii) Two-instance receipt on the committed values: the unit tick and
the double tick give distinct return durations for either state. -/
theorem return_tick_declaration_not_forced :
    unitTick.tau ≠ doubleTick.tau ∧
      unitTick.labSeconds (meanReturn 0) ≠ doubleTick.labSeconds (meanReturn 0) ∧
      unitTick.labSeconds (meanReturn 1) ≠ doubleTick.labSeconds (meanReturn 1) := by
  have hne : unitTick.tau ≠ doubleTick.tau := by norm_num
  exact ⟨hne, returnDuration_not_forced _ _ hne 0, returnDuration_not_forced _ _ hne 1⟩

end Boundaries

/-! ## (c iii) A return count is a distinct-event count

Row statement.  The return count touches the source clock/duration row
only: it orders returns on the abstract step index and attaches no
duration, no spacetime point, and no coupled action.  It discharges none
of the three open rows. -/

section ReturnCount

/-- The number of visits to state `i` at steps `1, …, n` of a realized
trajectory `ω`: each step index contributes one when the trajectory sits
at `i` and zero otherwise. -/
def returnCount (i : Fin 2) (ω : ℕ → Fin 2) : ℕ → ℕ
  | 0 => 0
  | n + 1 => returnCount i ω n + if ω (n + 1) = i then 1 else 0

@[simp] theorem returnCount_zero (i : Fin 2) (ω : ℕ → Fin 2) :
    returnCount i ω 0 = 0 := rfl

theorem returnCount_succ (i : Fin 2) (ω : ℕ → Fin 2) (n : ℕ) :
    returnCount i ω (n + 1) = returnCount i ω n + if ω (n + 1) = i then 1 else 0 := rfl

/-- Each step adds at most one return. -/
theorem returnCount_succ_le (i : Fin 2) (ω : ℕ → Fin 2) (n : ℕ) :
    returnCount i ω (n + 1) ≤ returnCount i ω n + 1 := by
  rw [returnCount_succ]
  split_ifs <;> omega

/-- The return count is monotone in the step index. -/
theorem returnCount_mono (i : Fin 2) (ω : ℕ → Fin 2) (n : ℕ) :
    returnCount i ω n ≤ returnCount i ω (n + 1) := by
  rw [returnCount_succ]
  omega

/-- At most one return per step: the count is bounded by the step index. -/
theorem returnCount_le (i : Fin 2) (ω : ℕ → Fin 2) (n : ℕ) :
    returnCount i ω n ≤ n := by
  induction n with
  | zero => simp
  | succ n ih => exact le_trans (returnCount_succ_le i ω n) (by omega)

/-- **Distinct-event additivity.**  The count over `n + m` steps is the
count over the first `n` steps plus the count of the shifted trajectory
over the next `m` steps: adjacent disjoint time blocks add with no shared
index.  The overlapping-window incidence statistic of
`QFT/SourceClockCandidate.lean` has no such decomposition, since
consecutive windows share an edge. -/
theorem returnCount_block_additive (i : Fin 2) (ω : ℕ → Fin 2) (n m : ℕ) :
    returnCount i ω (n + m) =
      returnCount i ω n + returnCount i (fun k => ω (n + k)) m := by
  induction m with
  | zero => simp
  | succ m ih =>
      rw [← Nat.add_assoc, returnCount_succ, ih, returnCount_succ, Nat.add_assoc,
        Nat.add_assoc n m 1]

end ReturnCount

end

/-! ## Axiom receipts -/

#print axioms acceptedStep_no_cycle
#print axioms periodic_iterate_is_fixed
#print axioms iterate_eventually_constant
#print axioms returnMass_hasSum_one
#print axioms meanReturn_hasSum
#print axioms kac_identity
#print axioms meanReturn_values
#print axioms meanReturn_ratio
#print axioms meanReturn_ratio_tick_free
#print axioms no_single_period
#print axioms returnMass_lt_one
#print axioms returnDuration_not_forced
#print axioms return_tick_declaration_not_forced
#print axioms returnCount_le
#print axioms returnCount_block_additive

end OPH.QFT.SourceRecurrenceClock
