import Mathlib.Algebra.BigOperators.Field
import Mathlib.Tactic

set_option autoImplicit false

/-!
# Ground-energy excess controls complete downward-energy records

The weights are unconditioned detector probabilities, including their original
preparation weights. A row may lose weight to a retained rejection outcome;
it may never acquire weight by conditioning on a successful event.

For a stationary finite source, insert its energy eigenvalues, diagonal state
weights and the squared matrix elements of a complete contraction instrument.
The theorem needs a lower energy bound and a mean energy-excess bound. It
does not require a spectral gap. The attachment of these arrays to a source,
stationarity, apparatus implementation and continuum convergence are separate
premises.
-/

namespace OPH.QFT.GroundRecordTail

open scoped BigOperators

variable {ι : Type*} [Fintype ι]

/-- Total raw weight of records that lower source energy by at least `η`. -/
noncomputable def downwardWeight (E : ι → ℝ) (w : ι → ι → ℝ) (η : ℝ) : ℝ :=
  ∑ i, ∑ j, if E j - E i ≤ -η then w i j else 0

/-- Energy-weighted detector estimate before division by the threshold. -/
theorem threshold_mul_downwardWeight_le
    (E p : ι → ℝ) (w : ι → ι → ℝ) (E₀ η : ℝ)
    (hfloor : ∀ i, E₀ ≤ E i)
    (hw : ∀ i j, 0 ≤ w i j)
    (hrow : ∀ i, ∑ j, w i j ≤ p i) :
    η * downwardWeight E w η ≤ ∑ i, p i * (E i - E₀) := by
  classical
  unfold downwardWeight
  rw [Finset.mul_sum]
  apply Finset.sum_le_sum
  intro i _
  calc
    η * (∑ j, if E j - E i ≤ -η then w i j else 0)
        ≤ (∑ j, w i j) * (E i - E₀) := by
      rw [Finset.mul_sum, Finset.sum_mul]
      apply Finset.sum_le_sum
      intro j _
      split_ifs with h
      · have hη : η ≤ E i - E₀ := by linarith [hfloor j]
        nlinarith [hw i j]
      · simpa using mul_nonneg (hw i j) (sub_nonneg.mpr (hfloor i))
    _ ≤ p i * (E i - E₀) :=
      mul_le_mul_of_nonneg_right (hrow i) (sub_nonneg.mpr (hfloor i))

/-- A gap-free bound on original downward-event probabilities. -/
theorem downwardWeight_le_excess_div
    (E p : ι → ℝ) (w : ι → ι → ℝ) (E₀ η ε : ℝ)
    (hη : 0 < η) (hfloor : ∀ i, E₀ ≤ E i)
    (hw : ∀ i j, 0 ≤ w i j)
    (hrow : ∀ i, ∑ j, w i j ≤ p i)
    (hexcess : (∑ i, p i * (E i - E₀)) ≤ ε) :
    downwardWeight E w η ≤ ε / η := by
  apply (le_div_iff₀ hη).mpr
  rw [mul_comm]
  exact le_trans (threshold_mul_downwardWeight_le E p w E₀ η hfloor hw hrow)
    hexcess

/-- An exact ground state has no downward events at any positive threshold. -/
theorem downwardWeight_eq_zero_of_zero_excess
    (E p : ι → ℝ) (w : ι → ι → ℝ) (E₀ η : ℝ)
    (hη : 0 < η) (hfloor : ∀ i, E₀ ≤ E i)
    (hw : ∀ i j, 0 ≤ w i j)
    (hrow : ∀ i, ∑ j, w i j ≤ p i)
    (hexcess : (∑ i, p i * (E i - E₀)) ≤ 0) :
    downwardWeight E w η = 0 := by
  apply le_antisymm
  · simpa using downwardWeight_le_excess_div E p w E₀ η 0 hη hfloor hw hrow hexcess
  · classical
    unfold downwardWeight
    apply Finset.sum_nonneg
    intro i _
    apply Finset.sum_nonneg
    intro j _
    split_ifs
    · exact hw i j
    · exact le_refl 0

end OPH.QFT.GroundRecordTail
