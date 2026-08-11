import InformationProjection.ReferenceNormalForm
import InformationProjection.SourceHistoryPacket
import FiniteConditionalRepair

/-!
# Source selection of the history reference

Issue B7 (#683).  The log-transition-action representation theorem fixes
its reference freedom by the target-relabeling normal form: the
step-uniform reference is the unique row-stochastic kernel invariant
under every relabeling of transition targets
(`ReferenceNormalForm.unique_invariant_reference`).  That normal form is
a convention about kernels; this module upgrades it to a source
product.

**Selection theorem.**  The Axiom-3 conditional-resampling repair kernel
of the *trivial visible datum* under the *counting reference* is exactly
the uniform kernel: repairing with no retained visible information, with
every configuration carrying unit source mass, is the step-uniform
reference.  Composed with the normal form, the representation theorem's
reference is therefore the repair law's own output rather than a free
choice.

Two controls prove both inputs load-bearing.  A nontrivial visible datum
breaks target-relabel invariance (the identity datum turns the repair
kernel into the identity kernel), and a non-counting reference breaks it
(a biased mass on two points weights the rows).  So the trivial datum
and the counting mass are exactly the source data that produce the
normal form.

**Multiplier boundary.**  The committed history packet matches its Gibbs
multiplier to the declared empirical mean action.  The nonidentifiability
receipt here proves that the packet's internal receipts cannot fix the
multiplier: the tilts at two distinct multipliers are both normalized
positive Gibbs laws on the committed two-state chain and they differ, so
only the declared constraint level separates them.

**Boundary.**  The counting reference is the source counting measure of
the committed repair-word packet; taking it as the mass law is a declared
identification of that committed input, and the trivial datum is the
no-visible-information case of the Axiom-3 fibre structure.  No physical
unit, clock, current, amplitude, field, or continuum object is
constructed, and the multiplier's constraint level stays a source
literal rather than a packet-internal function.  Issue #683 is closed bounded
under its scoped exits: the finite gaps are discharged here and in the
committed packet modules, and the physical attachments carry committed
obstructions with their production owned by the named successors.
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
/-- **Source selection.**  The Axiom-3 repair kernel of the trivial
visible datum under the counting reference is the uniform kernel. -/
theorem heatBath_counting_trivial_eq_uniform :
    heatBath (countingRef Ω) (trivialDatum Ω) = uniformKernel Ω := by
  funext x y
  unfold heatBath fiberMass countingRef trivialDatum uniformKernel
  simp [Finset.filter_true_of_mem, Finset.card_univ]

omit [DecidableEq Ω] in
/-- The composed statement: the repair-of-trivial-datum kernel is the
unique target-relabel-invariant row-stochastic kernel, so the
representation theorem's step-uniform reference is a source product of
the Axiom-3 repair law. -/
theorem reference_source_selected [Nonempty Ω] :
    heatBath (countingRef Ω) (trivialDatum Ω) = uniformKernel Ω ∧
      ∀ P : Ω → Ω → ℝ, RowStochastic P →
        (TargetRelabelInvariant P ↔ P = uniformKernel Ω) :=
  ⟨heatBath_counting_trivial_eq_uniform,
    fun P hP => relabel_invariant_iff_uniform P hP⟩

/-! ## Controls: both inputs are load-bearing -/

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

/-! ## The multiplier is not packet-determined -/

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

/-- **Multiplier nonidentifiability.**  The tilts at two distinct
weights are both strictly positive normalized laws on the committed
chain, and their mean actions differ, so the packet's internal receipts
admit every weight and only the declared constraint level separates
them. -/
theorem multiplier_not_packet_determined :
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
#print axioms OPH.InformationProjection.multiplier_not_packet_determined
#print axioms OPH.InformationProjection.heatBath_counting_trivial_eq_uniform
#print axioms OPH.InformationProjection.reference_source_selected
#print axioms OPH.InformationProjection.nontrivial_datum_not_invariant
#print axioms OPH.InformationProjection.noncounting_reference_not_invariant
