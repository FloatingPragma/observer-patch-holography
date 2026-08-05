import Time.TimeOrderLedger

/-!
# D1: finite ranked observer histories

This module provides the first bounded D1 construction above the A1 type
ledger.  A `RankedHistory` is an observer record order together with an
explicit natural-number rank that strictly increases on that order.  Such a
rank produces a monotone real clock coordinate.  The rank is supplied data:
no record order is claimed to admit one, and no physical duration is inferred.

The canonical `Fin n` chain is a nonvacuous finite model.  A separate discrete
control shows why the existence of a `ClockReadout` alone does not establish
that any two records are ordered or that any clock has ticked.
-/

namespace OPH.D1Time

open OPH.TimeOrderLedger

universe u

/-- A strict observer history with an explicitly supplied natural rank. -/
structure RankedHistory (Record : Type u) where
  order : ObserverRecordOrder Record
  rank : Record → ℕ
  rank_strict : ∀ {first second}, order.precedes first second →
    rank first < rank second

namespace RankedHistory

/-- The natural rank, cast to the reals, is a monotone clock coordinate. -/
def clock {Record : Type u} (history : RankedHistory Record) :
    ClockReadout history.order where
  read record := history.rank record
  strictlyMonotone := by
    intro first second hBefore
    exact_mod_cast history.rank_strict hBefore

@[simp] theorem clock_read {Record : Type u} (history : RankedHistory Record)
    (record : Record) : history.clock.read record = history.rank record := rfl

end RankedHistory

/-- The canonical strict order on a finite chain. -/
def finRecordOrder (n : ℕ) : ObserverRecordOrder (Fin n) where
  precedes first second := first < second
  irrefl := by simp
  trans := by
    intro first second third hFirstSecond hSecondThird
    exact hFirstSecond.trans hSecondThird

@[simp] theorem finRecordOrder_precedes_iff {n : ℕ} (first second : Fin n) :
    (finRecordOrder n).precedes first second ↔ first < second := Iff.rfl

/-- Every canonical finite chain has its literal index as a rank. -/
def finRankedHistory (n : ℕ) : RankedHistory (Fin n) where
  order := finRecordOrder n
  rank record := record.1
  rank_strict := by simp [finRecordOrder]

@[simp] theorem finRankedHistory_rank {n : ℕ} (record : Fin n) :
    (finRankedHistory n).rank record = record.1 := rfl

/-- The three-record control contains two genuine successive precedences. -/
theorem threeRecord_control :
    (finRankedHistory 3).order.precedes 0 1 ∧
      (finRankedHistory 3).order.precedes 1 2 := by
  change (0 : Fin 3) < 1 ∧ (1 : Fin 3) < 2
  decide

/-- On a discrete record order, a constant function satisfies the clock
interface vacuously.  This is a negative control: clock-interface existence
does not by itself imply a tick or a nontrivial history. -/
def discreteConstantClock (Record : Type u) (value : ℝ) :
    ClockReadout (ObserverRecordOrder.discrete Record) where
  read := fun _ => value
  strictlyMonotone := by simp [ObserverRecordOrder.discrete]

theorem discreteConstantClock_not_injective :
    ¬ Function.Injective (discreteConstantClock (Fin 2) 0).read := by
  intro hInjective
  have h : (0 : Fin 2) = 1 := hInjective rfl
  omega

#print axioms RankedHistory.clock
#print axioms threeRecord_control
#print axioms discreteConstantClock_not_injective

end OPH.D1Time
