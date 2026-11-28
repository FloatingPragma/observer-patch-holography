import Mathlib.Tactic

/-!
# Reading an oriented material current with a binary outside head

The two fixed exterior arrows make a full face reversal available exactly when
the shared arrow points to the destination. A quarter material pulse on an
enabled source/destination pair gives destination amplitude (b+i*a)/√2 up to
common phase; the three-quarter pulse gives (b-i*a)/√2. The physical reader
compares that outside head with its initial direction. These identities extend
over unresolved words or reference outcomes by finite addition. The executable
checks the face predicate on the actual 1045 preparation.
-/

namespace SourceLaw.ExteriorHeadCurrent

noncomputable def shortProbability (ar ai br bi : ℝ) : ℝ :=
  ((br - ai)^2 + (bi + ar)^2) / 2

noncomputable def longProbability (ar ai br bi : ℝ) : ℝ :=
  ((br + ai)^2 + (bi - ar)^2) / 2

/-- The imaginary part in the oriented material continuity equation. -/
def orientedCurrent (t ar ai br bi : ℝ) : ℝ := 2 * t * (ar * bi - ai * br)

/-- Counts on one enabled pair, with arbitrary complex source and destination
amplitudes. Frozen words give identical readings in both settings. -/
theorem pair_current (t ar ai br bi : ℝ) :
    t * (shortProbability ar ai br bi - longProbability ar ai br bi) =
      orientedCurrent t ar ai br bi := by
  unfold shortProbability longProbability orientedCurrent
  ring

/-- A fair setting bit scores a changed head in the short arm or a preserved
head in the long arm. All success probabilities here may be unnormalized. -/
theorem agreement_current {ι : Type*} [Fintype ι] (t mass : ℝ)
    (short long current : ι → ℝ) (perPair : ∀ i, t * (short i - long i) = current i) :
    2 * t * (((∑ i, short i) + (mass - ∑ i, long i)) / 2) - t * mass =
      ∑ i, current i := by
  calc
    _ = t * ((∑ i, short i) - ∑ i, long i) := by ring
    _ = ∑ i, t * (short i - long i) := by rw [← Finset.sum_sub_distrib, Finset.mul_sum]
    _ = _ := Finset.sum_congr rfl (fun i _ => perPair i)

#print axioms pair_current
#print axioms agreement_current

/-- An origin read after the inverse rest preparation has amplitude only
on the calibrated rest line in the route energy basis. This read removes
every other line of an arbitrary coherent filter. Apply pointwise to each
spectator and reference coordinate. -/
theorem collect_after_filter {ι : Type*} [Fintype ι] [DecidableEq ι]
    (rest : ι) (collect filter input : ι → ℂ)
    (other : ∀ i, i ≠ rest → collect i = 0) :
    ∑ i, collect i * (filter i * input i) = filter rest * (collect rest * input rest) := by
  rw [Finset.sum_eq_single rest]
  · ring
  · intro i _ hi
    simp [other i hi]
  · simp

#print axioms collect_after_filter

end SourceLaw.ExteriorHeadCurrent
