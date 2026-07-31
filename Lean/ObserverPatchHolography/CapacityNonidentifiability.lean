import Mathlib

/-!
# Bounded capacity-indexed counterfamily

This module formalises the all-rung arithmetic theorem behind the bounded
issue 551 counterfamily. The source packet fixes a twenty-four-record base.
Reversible identity, copy erasure, a two-class cap, and hidden spectator
multiplicity preserve base agreement, positivity, and the carrier bound while
giving different fixed sets. The complete three-axiom packet lift across
capacity is outside this structure.

The theorem separates four quantities:

* trial carrier dimension `D`;
* input coordinate `NIn D = log D`;
* correctable public-record count `M0`;
* output coordinate `NOut M0 = log M0`.

Typed self-identification implies equality of two readings once they have
been mapped to the same quantity. It does not construct or select either map.
No cosmological value or physical horizon identification appears here.
-/

namespace OPH.CapacityNonidentifiability

/-- The public dimension at extension rung `k`. -/
def publicDimension (k : ℕ) : ℕ := 24 * k

/-- A raw carrier with an additional spectator multiplicity `s`. -/
def rawDimension (k s : ℕ) : ℕ := 24 * k * s

/-- Input coordinate. Its argument remains explicitly typed as a dimension. -/
noncomputable def NIn (D : ℕ) : ℝ := Real.log D

/-- Output coordinate. Its argument remains explicitly typed as a record count. -/
noncomputable def NOut (M0 : ℕ) : ℝ := Real.log M0

/-- Equal integer readings give equal logarithmic coordinates. -/
theorem coordinates_eq_of_readings_eq {D M0 : ℕ} (h : D = M0) :
    NIn D = NOut M0 := by
  simp [NIn, NOut, h]

/-- Two differently typed descriptions together with their maps to one shared
    quantity. The `sameQuantity` field is the strange-loop identification
    obligation. -/
structure TypedSelfReadback (Simulator Simulated Quantity : Type*) where
  simulatorState : Simulator
  simulatedState : Simulated
  simulatorRead : Simulator → Quantity
  simulatedRead : Simulated → Quantity
  sameQuantity :
    simulatorRead simulatorState = simulatedRead simulatedState

/-- Once both sides are proved to read the same quantity, equality follows. -/
theorem strangeLoop_readings_eq
    {Simulator Simulated Quantity : Type*}
    (readback : TypedSelfReadback Simulator Simulated Quantity) :
    readback.simulatorRead readback.simulatorState =
      readback.simulatedRead readback.simulatedState :=
  readback.sameQuantity

/-- Reversible continuation retains every public record. -/
def identityM (k : ℕ) : ℕ := 24 * k

/-- Copy collapse retains one twenty-four-record class. -/
def erasureM (_k : ℕ) : ℕ := 24

/-- The two-class continuation retains at most two record classes. -/
def cappedM (k : ℕ) : ℕ := 24 * min k 2

/-- Hidden spectators do not change the public record count. -/
def spectatorM (k : ℕ) : ℕ := 24 * k

theorem identity_fixed (k : ℕ) :
    identityM k = publicDimension k := rfl

theorem erasure_fixed_iff {k : ℕ} (_hk : 0 < k) :
    erasureM k = publicDimension k ↔ k = 1 := by
  simp only [erasureM, publicDimension]
  omega

theorem capped_fixed_iff {k : ℕ} (hk : 0 < k) :
    cappedM k = publicDimension k ↔ k = 1 ∨ k = 2 := by
  simp only [cappedM, publicDimension]
  by_cases h : k ≤ 2
  · rw [Nat.min_eq_left h]
    constructor
    · intro _
      omega
    · intro _
      rfl
  · rw [Nat.min_eq_right (by omega)]
    constructor <;> intro h'
    · omega
    · omega

theorem spectator_fixed_iff {k s : ℕ} (hk : 0 < k) :
    rawDimension k s = spectatorM k ↔ s = 1 := by
  constructor
  · intro h
    have h' : (24 * k) * s = (24 * k) * 1 := by
      simpa [rawDimension, spectatorM, Nat.mul_assoc] using h
    exact Nat.mul_left_cancel (by omega) h'
  · rintro rfl
    simp [rawDimension, spectatorM]

/-- The exact common antecedent used by the countermodels: agreement at the
    twenty-four-record base, positivity, and the carrier bound. -/
structure DeclaredCapacityCompletion where
  M0 : ℕ → ℕ
  baseAgreement : M0 1 = 24
  positive : ∀ k, 0 < k → 0 < M0 k
  carrierBound : ∀ k, 0 < k → M0 k ≤ publicDimension k

def identityCompletion : DeclaredCapacityCompletion where
  M0 := identityM
  baseAgreement := rfl
  positive := by
    intro k hk
    simp [identityM]
    omega
  carrierBound := by
    intro k _hk
    simp [identityM, publicDimension]

def erasureCompletion : DeclaredCapacityCompletion where
  M0 := erasureM
  baseAgreement := rfl
  positive := by
    intro _k _hk
    simp [erasureM]
  carrierBound := by
    intro k hk
    simp [erasureM, publicDimension]
    omega

def cappedCompletion : DeclaredCapacityCompletion where
  M0 := cappedM
  baseAgreement := rfl
  positive := by
    intro k hk
    simp only [cappedM]
    have hmin : 0 < min k 2 := by omega
    omega
  carrierBound := by
    intro k _hk
    simp only [cappedM, publicDimension]
    exact Nat.mul_le_mul_left 24 (Nat.min_le_left k 2)

def FixedAt (completion : DeclaredCapacityCompletion) (k : ℕ) : Prop :=
  completion.M0 k = publicDimension k

/-- Identity and erasure satisfy the same declared base and bound but disagree
    at the second rung. -/
theorem sameAntecedent_differentFixedSets :
    FixedAt identityCompletion 2 ∧ ¬ FixedAt erasureCompletion 2 := by
  constructor
  · rfl
  · norm_num [FixedAt, erasureCompletion, erasureM, publicDimension]

/-- The capped completion supplies an exact multiple-zero control. -/
theorem capped_two_zeros_and_third_not :
    FixedAt cappedCompletion 1 ∧
      FixedAt cappedCompletion 2 ∧
      ¬ FixedAt cappedCompletion 3 := by
  norm_num [FixedAt, cappedCompletion, cappedM, publicDimension]

/-- A completion has one positive fixed rung. -/
def HasUniquePositiveFixedRung
    (completion : DeclaredCapacityCompletion) : Prop :=
  ∃ k, 0 < k ∧ FixedAt completion k ∧
    ∀ j, 0 < j → FixedAt completion j → j = k

theorem identity_has_no_unique_positive_fixed_rung :
    ¬ HasUniquePositiveFixedRung identityCompletion := by
  rintro ⟨k, hk, _hfixed, hunique⟩
  have h1 : (1 : ℕ) = k :=
    hunique 1 (by omega) (identity_fixed 1)
  have h2 : (2 : ℕ) = k :=
    hunique 2 (by omega) (identity_fixed 2)
  omega

/-- The declared antecedent does not force every admissible completion to have
    a unique slack zero. A continuation selector would be an additional
    source law. -/
theorem boundedCompletionClass_doesNotForceUniqueZero :
    ¬ ∀ completion : DeclaredCapacityCompletion,
      HasUniquePositiveFixedRung completion := by
  intro h
  exact identity_has_no_unique_positive_fixed_rung (h identityCompletion)

/-! ## Complete-lift layer

The executable complete-lift receipt transports terminal fibers, atom maps,
public sections, histories, joint kernels, meaning maps, feasible sets, and
extension and refinement controls across the generation-register family. The
statements below are its arithmetic skeleton. The oscillatory direction is
exhibited and excluded by the executable A2 extension-square control; the
saturation theorem covers the source-closed reading, and the class theorem
covers any admissibility reading that retains the reversible identity
completion. -/

/-- Parity-oscillation continuation: copy collapse at even rungs only. -/
def oscillationM (k : ℕ) : ℕ := if k % 2 = 1 then 24 * k else 24

theorem oscillation_fixed_iff_odd {k : ℕ} (hk : 0 < k) :
    oscillationM k = publicDimension k ↔ k % 2 = 1 := by
  unfold oscillationM publicDimension
  split <;> omega

/-- The oscillatory direction satisfies the same declared antecedent. -/
def oscillationCompletion : DeclaredCapacityCompletion where
  M0 := oscillationM
  baseAgreement := rfl
  positive := by
    intro k hk
    unfold oscillationM
    split <;> omega
  carrierBound := by
    intro k hk
    unfold oscillationM publicDimension
    split <;> omega

/-- Oscillation and identity share the antecedent and disagree at rung two. -/
theorem oscillation_identity_differentFixedSets :
    FixedAt identityCompletion 2 ∧ ¬ FixedAt oscillationCompletion 2 := by
  constructor
  · rfl
  · norm_num [FixedAt, oscillationCompletion, oscillationM, publicDimension]

/-- Source-closed reading: a completion that saturates every rung has a
    degenerate slack zero set, so no unique positive fixed rung exists. -/
theorem sourceClosed_no_unique_positive_fixed_rung
    (completion : DeclaredCapacityCompletion)
    (hsat : ∀ k, 0 < k → completion.M0 k = publicDimension k) :
    ¬ HasUniquePositiveFixedRung completion := by
  rintro ⟨k, hk, _hfixed, hunique⟩
  have h1 : (1 : ℕ) = k := hunique 1 (by omega) (hsat 1 (by omega))
  have h2 : (2 : ℕ) = k := hunique 2 (by omega) (hsat 2 (by omega))
  omega

/-- The reversible identity completion saturates every rung. -/
theorem identity_saturates (k : ℕ) (_hk : 0 < k) :
    identityCompletion.M0 k = publicDimension k := rfl

/-- Complete-class theorem: any admissibility reading of the complete lift
    that retains the reversible identity completion does not entail a unique
    slack zero. The source-closed and widened executable readings both retain
    it, so neither selects a capacity; a rung selector is an additional
    source law. -/
theorem completeClass_doesNotEntailUniqueZero
    (Admissible : DeclaredCapacityCompletion → Prop)
    (hid : Admissible identityCompletion) :
    ¬ ∀ completion, Admissible completion →
        HasUniquePositiveFixedRung completion := by
  intro h
  exact identity_has_no_unique_positive_fixed_rung (h identityCompletion hid)

-- Axiom audit.
#print axioms coordinates_eq_of_readings_eq
#print axioms strangeLoop_readings_eq
#print axioms erasure_fixed_iff
#print axioms capped_fixed_iff
#print axioms spectator_fixed_iff
#print axioms sameAntecedent_differentFixedSets
#print axioms capped_two_zeros_and_third_not
#print axioms boundedCompletionClass_doesNotForceUniqueZero
#print axioms oscillation_fixed_iff_odd
#print axioms oscillation_identity_differentFixedSets
#print axioms sourceClosed_no_unique_positive_fixed_rung
#print axioms completeClass_doesNotEntailUniqueZero

end OPH.CapacityNonidentifiability
