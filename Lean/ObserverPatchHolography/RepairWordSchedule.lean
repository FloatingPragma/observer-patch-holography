import Mathlib.InformationTheory.KullbackLeibler.Basic
import Mathlib.Probability.Distributions.Uniform

/-!
# Finite repair-word schedule selected by information projection

This file isolates one conditional finite consequence of the proposed A1-R
and A2-R clauses in `docs/CANONICAL_REPAIR_LAW_RFC.md`.

Suppose A1-R has produced the complete, quotient-deduplicated type `Move` of
primitive observer-visible repair attempts. Suppose A2-R has proved temporal
completeness in the strong finite form that every probability law on every
word space `Fin k → Move` is feasible. If A3 minimizes Kullback-Leibler
divergence to source counting measure on that full simplex, Gibbs' equality
case forces the selected law to be uniform.

The result does not require a transitive geometric action on `Move`. Several
geometric move orbits receive probability in proportion to the number of
actual primitive event identities they contain. Raw presentation duplicates
cannot change the law because they are absent from the quotient event type.

These are ordinary conditional theorems. They do not prove that the current
three-axiom basis supplies a complete move grammar, that all finite words are
temporally admissible, or that source counting is the A3 reference. They do
not adopt A1-R or A2-R, identify a physical clock, or prove refinement
naturality.
-/

namespace ObserverPatchHolography.RepairWordSchedule

open MeasureTheory ProbabilityTheory InformationTheory
open scoped ENNReal

universe u

variable (Move : Type u) [Fintype Move] [Nonempty Move]

/-- A length-`k` word in the quotient-visible primitive attempt alphabet. -/
abbrev Word (k : ℕ) := Fin k → Move

/-- Source counting measure on the complete finite word space. -/
noncomputable def uniformWord (k : ℕ) : PMF (Word Move k) :=
  PMF.uniformOfFinset Finset.univ Finset.univ_nonempty

/-- The one-event source counting measure. -/
noncomputable def uniformMove : PMF Move :=
  PMF.uniformOfFinset Finset.univ Finset.univ_nonempty

/-- A PMF regarded as a measure on the discrete measurable space. Keeping the
measurable structure explicit prevents any hidden sigma-algebra selector. -/
noncomputable def discreteMeasure {α : Type*} (p : PMF α) : @Measure α ⊤ :=
  @PMF.toMeasure α ⊤ p

instance discreteMeasure_isProbabilityMeasure {α : Type*} (p : PMF α) :
    @IsProbabilityMeasure α ⊤ (discreteMeasure p) := by
  change @IsProbabilityMeasure α ⊤ (@PMF.toMeasure α ⊤ p)
  infer_instance

/-- A3's finite-word optimization packet after A1-R has fixed the complete
event identities and A2-R has certified the full word simplex. The universal
quantifier in `selected_minimal` is the load-bearing temporal-freedom premise:
there is no state-dependent deletion or cross-position coupling left in the
feasible family. -/
structure FullWordKLProjection (k : ℕ) where
  selected : PMF (Word Move k)
  selected_minimal :
    ∀ p : PMF (Word Move k),
      klDiv (discreteMeasure selected) (discreteMeasure (uniformWord Move k)) ≤
        klDiv (discreteMeasure p) (discreteMeasure (uniformWord Move k))

/-- The discrete-measure embedding of PMFs is injective. -/
theorem discreteMeasure_injective {α : Type*} :
    Function.Injective (@discreteMeasure α) := by
  exact @PMF.toMeasure_injective α ⊤ inferInstance

/-- Gibbs' equality case removes every finite-word schedule parameter: the
selected law is exactly counting measure on the complete word set. -/
theorem selected_eq_uniformWord (k : ℕ)
    (selection : FullWordKLProjection Move k) :
    selection.selected = uniformWord Move k := by
  have hle := selection.selected_minimal (uniformWord Move k)
  have hzero :
      klDiv (discreteMeasure selection.selected)
        (discreteMeasure (uniformWord Move k)) = 0 := by
    rw [klDiv_self] at hle
    exact bot_unique hle
  have hmeasure :
      discreteMeasure selection.selected =
        discreteMeasure (uniformWord Move k) :=
    klDiv_eq_zero_iff.mp hzero
  exact discreteMeasure_injective hmeasure

/-- Every word has the same exact counting weight. -/
theorem uniformWord_apply (k : ℕ) (w : Word Move k) :
    uniformWord Move k w =
      ((Fintype.card Move : ℝ≥0∞) ^ k)⁻¹ := by
  simp [uniformWord]

/-- The one-event counting law is `1 / |Move|`. -/
theorem uniformMove_apply (m : Move) :
    uniformMove Move m = (Fintype.card Move : ℝ≥0∞)⁻¹ := by
  simp [uniformMove]

/-- Every quotient-visible primitive identity has positive one-event weight. -/
theorem uniformMove_positive (m : Move) :
    0 < uniformMove Move m := by
  rw [uniformMove_apply]
  exact ENNReal.inv_pos.mpr (ENNReal.natCast_ne_top (Fintype.card Move))

/-- The full-word counting reference factorizes exactly into the product of
its one-event counting weights. Thus the selected finite-dimensional family
is IID once the full-simplex premise has been proved at every length. -/
theorem uniformWord_factorizes (k : ℕ) (w : Word Move k) :
    uniformWord Move k w = ∏ i : Fin k, uniformMove Move (w i) := by
  rw [uniformWord_apply]
  simp [uniformMove_apply, Finset.prod_const, Fintype.card_fin,
    ENNReal.inv_pow]

/-- A3 therefore selects the exact product law on every certified full word
simplex. -/
theorem selected_factorizes (k : ℕ)
    (selection : FullWordKLProjection Move k) (w : Word Move k) :
    selection.selected w = ∏ i : Fin k, uniformMove Move (w i) := by
  rw [selected_eq_uniformWord, uniformWord_factorizes]

/-- The selected finite-word law has full support. -/
theorem selected_full_support (k : ℕ)
    (selection : FullWordKLProjection Move k) (w : Word Move k) :
    0 < selection.selected w := by
  rw [selected_eq_uniformWord, uniformWord_apply]
  exact ENNReal.inv_pos.mpr
    (ENNReal.pow_ne_top (ENNReal.natCast_ne_top (Fintype.card Move)))

/-- Presentation relabeling cannot change a counting-selected word weight.
Only the quotient-visible event identities occur in this statement; a raw
duplicate presentation is not an additional event. -/
theorem uniformWord_presentation_invariant (k : ℕ) (e : Move ≃ Move)
    (w : Word Move k) :
    uniformWord Move k (e ∘ w) = uniformWord Move k w := by
  rw [uniformWord_apply, uniformWord_apply]

/-- The source-counting probability of any declared subset of primitive
events is its exact fraction of the complete event identity set. This fixes
relative weights across several geometric orbits without assuming that one
geometric group acts transitively on the entire alphabet. -/
theorem uniformMove_subset_mass (s : Finset Move) :
    ∑ m ∈ s, uniformMove Move m =
      (s.card : ℝ≥0∞) / (Fintype.card Move : ℝ≥0∞) := by
  simp [uniformMove_apply, div_eq_mul_inv]

/-- One-step prefix consistency: summing over every possible final event in
a length-`k+1` counting word recovers the weight of its length-`k` prefix. -/
theorem uniformWord_prefix_consistent (k : ℕ) (w : Word Move k) :
    ∑ m : Move, uniformWord Move (k + 1) (Fin.snoc w m) =
      uniformWord Move k w := by
  simp_rw [uniformWord_apply]
  simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  have hcard : (Fintype.card Move : ℝ≥0∞) ≠ 0 := by
    exact_mod_cast Fintype.card_ne_zero
  have hcardTop : (Fintype.card Move : ℝ≥0∞) ≠ ⊤ :=
    ENNReal.natCast_ne_top (Fintype.card Move)
  have hpow0 : (Fintype.card Move : ℝ≥0∞) ^ k ≠ 0 :=
    pow_ne_zero k hcard
  have hpowTop : (Fintype.card Move : ℝ≥0∞) ^ k ≠ ⊤ :=
    ENNReal.pow_ne_top hcardTop
  rw [pow_succ, ENNReal.mul_inv (Or.inl hpow0) (Or.inl hpowTop)]
  calc
    (Fintype.card Move : ℝ≥0∞) *
          (((Fintype.card Move : ℝ≥0∞) ^ k)⁻¹ *
            (Fintype.card Move : ℝ≥0∞)⁻¹) =
        ((Fintype.card Move : ℝ≥0∞) ^ k)⁻¹ *
          ((Fintype.card Move : ℝ≥0∞) *
            (Fintype.card Move : ℝ≥0∞)⁻¹) := by ring
    _ = ((Fintype.card Move : ℝ≥0∞) ^ k)⁻¹ := by
      rw [ENNReal.mul_inv_cancel hcard hcardTop, mul_one]

/-- A3-selected laws on consecutive certified full word simplexes have the
same one-step prefix consistency. No separate Markov or stationarity premise
is used after full-simplex minimality has been established at both lengths. -/
theorem selected_prefix_consistent (k : ℕ)
    (selectionK : FullWordKLProjection Move k)
    (selectionK1 : FullWordKLProjection Move (k + 1))
    (w : Word Move k) :
    ∑ m : Move, selectionK1.selected (Fin.snoc w m) =
      selectionK.selected w := by
  rw [selected_eq_uniformWord, selected_eq_uniformWord]
  exact uniformWord_prefix_consistent Move k w

/-! ## Axiom audit

Expected output: standard Mathlib axioms only. In particular there is no
project-level axiom and no `sorryAx`.
-/

#print axioms selected_eq_uniformWord
#print axioms selected_factorizes
#print axioms selected_full_support
#print axioms uniformMove_subset_mass
#print axioms uniformWord_prefix_consistent
#print axioms selected_prefix_consistent

end ObserverPatchHolography.RepairWordSchedule
