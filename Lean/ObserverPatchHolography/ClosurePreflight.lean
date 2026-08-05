import Mathlib
import ObserverPatchHolography.CapacityClosurePrinciple

/-!
# Source-only closure preflight lemmas

Two exact distinctions used by the issue-708 inventory are formalized here.

* A multiplicative readback normalized at any positive seed fixes that seed
  uniquely when the exponent differs from one.  The fixed point therefore
  does not select the normalization constant by itself.
* Two readouts may share one codomain while failing the commuting-square
  obligation.  A common result type is weaker than a same-quantity theorem.

The file supplies no physical attachment, candidate selector, or numerical
prediction.
-/

namespace OPH.ClosurePreflight

open Function

/-! ## Seeded normalization -/

/-- Multiplicative readback normalized at an arbitrary positive seed. -/
noncomputable def seededReadback (seed exponent value : ℝ) : ℝ :=
  seed * Real.exp (exponent * Real.log (value / seed))

/-- Log coordinates turn the seeded readback into multiplication by the
    declared exponent. -/
theorem seededReadback_log_coords
    {seed : ℝ} (hseed : 0 < seed) (exponent value : ℝ) :
    Real.log (seededReadback seed exponent value / seed) =
      exponent * Real.log (value / seed) := by
  unfold seededReadback
  rw [mul_div_cancel_left₀ _ hseed.ne', Real.log_exp]

/-- A positive seeded readback with exponent other than one fixes exactly the
    supplied seed. -/
theorem seededReadback_fixedPt_iff
    {seed exponent value : ℝ}
    (hseed : 0 < seed)
    (hexponent : exponent ≠ 1)
    (hvalue : 0 < value) :
    seededReadback seed exponent value = value ↔ value = seed := by
  have hratio : 0 < value / seed := div_pos hvalue hseed
  constructor
  · intro hfixed
    have hlog :
        Real.log (value / seed) =
          exponent * Real.log (value / seed) := by
      have hcoords := seededReadback_log_coords hseed exponent value
      rw [hfixed] at hcoords
      exact hcoords
    have hzero : Real.log (value / seed) = 0 := by
      have hfactor :
          (exponent - 1) * Real.log (value / seed) = 0 := by
        linear_combination -hlog
      rcases mul_eq_zero.mp hfactor with hexponentZero | hlogZero
      · exact absurd (by linarith : exponent = 1) hexponent
      · exact hlogZero
    have hratioOne : value / seed = 1 :=
      Real.eq_one_of_pos_of_log_eq_zero hratio hzero
    field_simp at hratioOne
    exact hratioOne
  · rintro rfl
    unfold seededReadback
    rw [div_self hseed.ne', Real.log_one, mul_zero, Real.exp_zero, mul_one]

/-! ## Same codomain is not a commuting square -/

/-- An outer readout into the shared codomain `ℝ`. -/
def outerUnitRead (_state : Unit) : ℝ := 0

/-- An inner readout into the same shared codomain `ℝ`. -/
def innerUnitRead (_state : Unit) : ℝ := 1

/-- Sharing one codomain and an equivalence of state types does not force the
    readback square to commute. -/
theorem sameCodomain_noncommuting_witness :
    ∃ (outerRead innerRead : Unit → ℝ) (identify : Unit ≃ Unit),
      ¬ ∀ outer, innerRead (identify outer) = outerRead outer := by
  refine ⟨outerUnitRead, innerUnitRead, Equiv.refl Unit, ?_⟩
  intro commutes
  have contradiction := commutes Unit.unit
  norm_num [outerUnitRead, innerUnitRead] at contradiction

/-- In particular, the concrete shared-codomain readouts above cannot be
    packaged as an `EquivalentReadback` with the identity state map. -/
theorem unitReadouts_fail_sameQuantity :
    ¬ ∀ outer,
      innerUnitRead ((Equiv.refl Unit) outer) = outerUnitRead outer := by
  intro commutes
  have contradiction := commutes Unit.unit
  norm_num [outerUnitRead, innerUnitRead] at contradiction

-- Axiom audit: only Mathlib's standard classical axioms may appear.
#print axioms seededReadback_log_coords
#print axioms seededReadback_fixedPt_iff
#print axioms sameCodomain_noncommuting_witness
#print axioms unitReadouts_fail_sameQuantity

end OPH.ClosurePreflight
