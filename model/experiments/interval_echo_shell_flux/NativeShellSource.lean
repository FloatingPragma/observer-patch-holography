import Mathlib.Algebra.BigOperators.Field
import Mathlib.Tactic.Linarith

open scoped BigOperators

namespace OPH.NativeShellSource

variable {V : Type*} [Fintype V] [DecidableEq V]

/-- One tagged episode uses a material flight or a fresh emission at the origin.
The flight subkernel includes rejected contacts and empty offers. -/
def step (o : V) (flight : V → V → ℚ) (collect : V → ℚ) (i j : V) : ℚ :=
  flight i j + if j = o then collect i else 0

/-- Net material outflow at a vertex in units of transported tokens per appointment. -/
def outflow (pi : V → ℚ) (flight : V → V → ℚ) (v : V) : ℚ :=
  (∑ w, pi v * flight v w) - (∑ w, pi w * flight w v)

/-- Stationarity of the complete process and actual flight/collector accounting
give the source independently of shell areas or the angular current distribution. -/
theorem stationary_outflow (o : V) (pi collect : V → ℚ) (flight : V → V → ℚ)
    (hrow : ∀ v, (∑ w, flight v w) + collect v = 1)
    (hstat : ∀ v, (∑ w, pi w * step o flight collect w v) = pi v)
    (v : V) :
    outflow pi flight v =
      (if v = o then ∑ w, pi w * collect w else 0) - pi v * collect v := by
  have hout : (∑ w, pi v * flight v w) + pi v * collect v = pi v := by
    rw [← Finset.mul_sum, ← mul_add, hrow, mul_one]
  have hin := hstat v
  by_cases hv : v = o
  · simp only [step, hv, if_pos, mul_add, Finset.sum_add_distrib] at hin
    subst v
    simp [outflow]
    linarith
  · simp only [step, if_neg hv, add_zero] at hin
    simp only [outflow, if_neg hv]
    linarith

/-- The vertex part of the finite layered Gauss premise on every constrained
layer with no collector. The emitter is the sole internal source. -/
theorem steady_until_collector (o : V) (pi collect : V → ℚ)
    (flight : V → V → ℚ) (layer : V → ℤ) (depth : ℤ)
    (hrow : ∀ v, (∑ w, flight v w) + collect v = 1)
    (hstat : ∀ v, (∑ w, pi w * step o flight collect w v) = pi v)
    (external : ∀ v, layer v ≤ depth → collect v = 0)
    (v : V) (hv : layer v ≤ depth) :
    outflow pi flight v = if v = o then ∑ w, pi w * collect w else 0 := by
  rw [stationary_outflow o pi collect flight hrow hstat, external v hv, mul_zero, sub_zero]

end OPH.NativeShellSource
