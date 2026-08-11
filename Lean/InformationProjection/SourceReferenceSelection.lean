import InformationProjection.ReferenceNormalForm
import InformationProjection.SourceHistoryPacket
import FiniteConditionalRepair

/-!
# Conditional realization and boundary of the history reference

Issue B7 (#683).  The log-transition-action representation theorem fixes
its reference freedom by the target-relabeling normal form: the
step-uniform reference is the unique row-stochastic kernel invariant
under every relabeling of transition targets
(`ReferenceNormalForm.unique_invariant_reference`).  That normal form is
a convention about kernels.  This module proves a conditional realization
of that kernel and records two non-uniqueness boundaries; it does not derive
the required reference data from the OPH source.

**Conditional realization theorem.**  If the *trivial visible datum* and
the *counting reference* are supplied, their conditional-resampling repair
kernel is exactly the uniform transition kernel.  Composed with the normal
form, this realizes the transition part of the representation theorem's
step-uniform reference.  Axiom 3 conditions relative to a supplied exact
reference, however, so this equality does not prove that the source selects
the counting reference or the trivial datum.

Two separating controls show that changing the supplied data can change the
answer.  The identity visible datum turns the repair kernel into the identity
kernel, and a biased mass on two points weights the rows.  They do not prove
that the inputs are uniquely selected.  Indeed every positive constant
rescaling of the counting reference gives the same conditional kernel, and
two distinct positive normalized initial laws give distinct path references
despite sharing the same uniform transition kernel.

**Multiplier boundary.**  The committed history packet matches its Gibbs
multiplier to the declared empirical mean action.  The sensitivity receipt
here proves the narrower two-point statement: tilts at two distinct
multipliers are normalized positive Gibbs laws on the committed two-state
chain and have different mean actions.  This is evidence that the mean
responds to the multiplier, not a nonidentifiability theorem.  The existing
intermediate-value receipt matches a multiplier to the declared empirical
mean, and `SourceHistoryPacket` proves the positive exponential parameter is
unique at that supplied target.  No theorem here proves that the source
selects the constraint observable or level.

**Boundary.**  Taking the counting mass and the trivial datum as the inputs
is a declared identification.  The theorem selects neither the initial path
law nor an absolute normalization of the reference mass.  No physical unit,
clock, current, amplitude, field, stationary-phase rule, or continuum object
is constructed, and the multiplier's constraint level stays a declared
source-derived postprocessing literal rather than a packet-internal function.
These statements therefore do not by themselves discharge the source-
reference, real-enrichment, or physical-attachment obligations of issue
#683.
-/

namespace OPH.InformationProjection

open OPH.Thermodynamics
open OPH.InformationProjection.ReferenceNormalForm

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]

noncomputable section

/-- The counting reference: unit source mass on every configuration. -/
def countingRef (Ω : Type*) : Ω → ℝ := fun _ => 1

/-- The trivial visible datum: no information is retained by repair. -/
def trivialDatum (Ω : Type*) : Ω → Unit := fun _ => ()

omit [DecidableEq Ω] in
/-- **Conditional uniform-kernel realization.**  The repair kernel of the
supplied trivial visible datum under the supplied counting reference is the
uniform transition kernel. -/
theorem heatBath_counting_trivial_eq_uniform :
    heatBath (countingRef Ω) (trivialDatum Ω) = uniformKernel Ω := by
  funext x y
  unfold heatBath fiberMass countingRef trivialDatum uniformKernel
  simp [Finset.filter_true_of_mem, Finset.card_univ]

omit [DecidableEq Ω] in
/-- The composed conditional statement: the repair-of-trivial-datum kernel
for the supplied counting reference is the unique target-relabel-invariant
row-stochastic kernel.  The statement does not produce its two inputs or an
initial path law. -/
theorem reference_realized_under_counting_trivial_inputs [Nonempty Ω] :
    heatBath (countingRef Ω) (trivialDatum Ω) = uniformKernel Ω ∧
      ∀ P : Ω → Ω → ℝ, RowStochastic P →
        (TargetRelabelInvariant P ↔ P = uniformKernel Ω) :=
  ⟨heatBath_counting_trivial_eq_uniform,
    fun P hP => relabel_invariant_iff_uniform P hP⟩

/-! ## Controls and non-uniqueness boundaries -/

/-- With the identity visible datum the repair kernel is the identity
kernel. -/
theorem heatBath_counting_id (x y : Bool) :
    heatBath (countingRef Bool) (id : Bool → Bool) x y =
      if y = x then 1 else 0 := by
  unfold heatBath fiberMass countingRef
  cases x <;> cases y <;>
    norm_num [Finset.sum_filter, Fintype.sum_bool] <;> decide

/-- **Datum control.**  A nontrivial visible datum breaks
target-relabel invariance: retaining the full configuration makes the
repair kernel the identity, which the target swap distinguishes. -/
theorem nontrivial_datum_not_invariant :
    ¬ TargetRelabelInvariant
      (heatBath (countingRef Bool) (id : Bool → Bool)) := by
  intro h
  have hswap := h (Equiv.swap true false) true true
  rw [heatBath_counting_id, heatBath_counting_id] at hswap
  simp [Equiv.swap_apply_left] at hswap

/-- The biased two-point reference. -/
def biasedRef : Bool → ℝ := fun x => if x then 2 else 1

/-- **Reference control.**  A non-counting reference breaks
target-relabel invariance: the biased mass weights the rows of the
trivial-datum repair kernel. -/
theorem noncounting_reference_not_invariant :
    ¬ TargetRelabelInvariant
      (heatBath biasedRef (trivialDatum Bool)) := by
  intro h
  have hswap := h (Equiv.swap true false) true true
  unfold heatBath fiberMass biasedRef trivialDatum at hswap
  simp [Equiv.swap_apply_left] at hswap
  norm_num at hswap

/-- A constant rescaling of the counting reference. -/
def scaledCountingRef (c : ℝ) (Ω : Type*) : Ω → ℝ := fun _ => c

omit [DecidableEq Ω] in
/-- **Reference-scale non-uniqueness.**  Every positive constant rescaling of
the counting reference gives the same uniform conditional-resampling kernel.
Consequently the kernel equality does not select an absolute reference-mass
normalization. -/
theorem heatBath_scaledCounting_trivial_eq_uniform [Nonempty Ω]
    {c : ℝ} (hc : 0 < c) :
    heatBath (scaledCountingRef c Ω) (trivialDatum Ω) = uniformKernel Ω := by
  funext x y
  have hcard : (Fintype.card Ω : ℝ) ≠ 0 := by
    exact_mod_cast Fintype.card_ne_zero
  have hc0 : c ≠ 0 := ne_of_gt hc
  unfold heatBath fiberMass scaledCountingRef trivialDatum uniformKernel
  simp [Finset.filter_true_of_mem, Finset.card_univ]
  field_simp [hcard, hc0]

/-- First positive normalized initial law for the path-reference control. -/
def initialLawLeft : Bool → ℝ := fun x => if x then 1 / 3 else 2 / 3

/-- Second positive normalized initial law for the path-reference control. -/
def initialLawRight : Bool → ℝ := fun x => if x then 2 / 3 else 1 / 3

theorem initialLawLeft_pos : ∀ x, 0 < initialLawLeft x := by
  intro x
  cases x <;> norm_num [initialLawLeft]

theorem initialLawRight_pos : ∀ x, 0 < initialLawRight x := by
  intro x
  cases x <;> norm_num [initialLawRight]

theorem initialLawLeft_sum : ∑ x, initialLawLeft x = 1 := by
  norm_num [initialLawLeft, Fintype.sum_bool]

theorem initialLawRight_sum : ∑ x, initialLawRight x = 1 := by
  norm_num [initialLawRight, Fintype.sum_bool]

/-- The all-false one-step path used to distinguish the two references. -/
def allFalsePath : PathSpace Bool 1 := fun _ => false

/-- **Initial-law non-uniqueness.**  Two distinct positive normalized initial
laws combined with the same uniform transition kernel give distinct
step-uniform path references.  Selecting the uniform transition kernel alone
therefore does not select the whole history reference. -/
theorem uniform_transition_does_not_determine_path_reference :
    stepUniformRef initialLawLeft 1 ≠ stepUniformRef initialLawRight 1 := by
  intro h
  have hp := congrFun h allFalsePath
  norm_num [stepUniformRef, initialLawLeft, initialLawRight, allFalsePath] at hp

/-! ## Two-point multiplier sensitivity -/

/-- The exact rational tilt of the committed chain law at tilt weight
`c` (the weight is `exp (-lambda)` of the real packet). -/
def tiltAt (c : ℚ) (g : Fin 8) : ℚ :=
  sourceTauChainQ g * c ^ sourceAction g /
    (∑ h, sourceTauChainQ h * c ^ sourceAction h)

/-- Positivity of the exact tilt at every positive weight. -/
theorem tiltAt_pos {c : ℚ} (hc : 0 < c) (g : Fin 8) : 0 < tiltAt c g := by
  unfold tiltAt
  refine div_pos (mul_pos (sourceTauChainQ_pos g) (pow_pos hc _)) ?_
  refine Finset.sum_pos (fun h _ => mul_pos (sourceTauChainQ_pos h)
    (pow_pos hc _)) ⟨0, Finset.mem_univ 0⟩

/-- Normalization of the exact tilt at every positive weight. -/
theorem tiltAt_sum {c : ℚ} (hc : 0 < c) : ∑ g, tiltAt c g = 1 := by
  unfold tiltAt
  rw [← Finset.sum_div]
  refine div_self (ne_of_gt ?_)
  refine Finset.sum_pos (fun h _ => mul_pos (sourceTauChainQ_pos h)
    (pow_pos hc _)) ⟨0, Finset.mem_univ 0⟩

/-- **Two-point multiplier sensitivity.**  The tilts at two distinct weights
are both strictly positive normalized laws on the committed chain, and their
mean actions differ.  This theorem is compatible with uniqueness at a
supplied target mean; it proves neither that uniqueness nor source selection
of the target mean. -/
theorem committed_tilts_have_distinct_mean_actions :
    (∀ g, 0 < tiltAt 1 g) ∧ (∑ g, tiltAt 1 g = 1) ∧
      (∀ g, 0 < tiltAt (1/2) g) ∧ (∑ g, tiltAt (1/2) g = 1) ∧
      (∑ g, tiltAt 1 g * (sourceAction g : ℚ)) ≠
        (∑ g, tiltAt (1/2) g * (sourceAction g : ℚ)) := by
  have hS0 : sourceAction 0 = 0 := rfl
  have hS1 : sourceAction 1 = 1 := rfl
  have hS2 : sourceAction 2 = 2 := rfl
  have hS3 : sourceAction 3 = 1 := rfl
  have hS4 : sourceAction 4 = 1 := rfl
  have hS5 : sourceAction 5 = 2 := rfl
  have hS6 : sourceAction 6 = 1 := rfl
  have hS7 : sourceAction 7 = 0 := rfl
  refine ⟨fun g => tiltAt_pos one_pos g, tiltAt_sum one_pos,
    fun g => tiltAt_pos (by norm_num) g, tiltAt_sum (by norm_num), ?_⟩
  simp only [tiltAt, Fin.sum_univ_eight, sourceTauChainQ_apply_0,
    sourceTauChainQ_apply_1, sourceTauChainQ_apply_2,
    sourceTauChainQ_apply_3, sourceTauChainQ_apply_4,
    sourceTauChainQ_apply_5, sourceTauChainQ_apply_6,
    sourceTauChainQ_apply_7, hS0, hS1, hS2, hS3, hS4, hS5, hS6, hS7]
  norm_num

end

end OPH.InformationProjection

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.InformationProjection.committed_tilts_have_distinct_mean_actions
#print axioms OPH.InformationProjection.heatBath_counting_trivial_eq_uniform
#print axioms OPH.InformationProjection.reference_realized_under_counting_trivial_inputs
#print axioms OPH.InformationProjection.nontrivial_datum_not_invariant
#print axioms OPH.InformationProjection.noncounting_reference_not_invariant
#print axioms OPH.InformationProjection.heatBath_scaledCounting_trivial_eq_uniform
#print axioms OPH.InformationProjection.uniform_transition_does_not_determine_path_reference
