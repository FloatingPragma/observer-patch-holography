import FiniteConditionalRepair
import ObservableNormalForms.ConditionalResampling
import EventAlgebra.InstrumentSelectionByCertainStates

/-!
# Coherence of finite conditional repair

The heat-bath probabilities equal the conditional-resampling probabilities
for the same positive weights and observation map. This identity transfers
kernel recognition, relative-entropy contraction, and movement of point laws.
The statements concern finite kernels; physical identification of the weights,
observation map, and repair dynamics is outside their scope.
-/

namespace OPH.FiniteConditionalRepairBridge

open ObservableNormalForms OPH.Thermodynamics

variable {Ω B : Type*} [Fintype Ω] [DecidableEq Ω] [DecidableEq B]

/-- Package the reference weights and visible datum as a weighted observation. -/
def toObservation (π : Ω → ℝ) (b : Ω → B) (hπ : ∀ x, 0 < π x) :
    FiniteWeightedObservation Ω B := ⟨π, hπ, b⟩

omit [DecidableEq Ω] in
/-- The two conditional-repair definitions have exactly the same probabilities. -/
theorem heatBath_eq_conditionalResamplingKernel
    (π : Ω → ℝ) (b : Ω → B) (hπ : ∀ x, 0 < π x) :
    heatBath π b = (toObservation π b hπ).conditionalResamplingKernel.probability := by
  rfl

/-- The heat-bath kernel with its row-normalization and nonnegativity proofs. -/
noncomputable def heatBathKernel (π : Ω → ℝ) (b : Ω → B)
    (hπ : ∀ x, 0 < π x) : FiniteMarkovKernel Ω where
  probability := heatBath π b
  probability_nonneg := by
    intro x y
    unfold heatBath
    split_ifs
    · exact div_nonneg (hπ y).le (fiberMass_pos hπ x).le
    · exact le_rfl
  probability_sum_one := heatBath_row_sum hπ

/-- Equality of probabilities identifies the packaged Markov kernels. -/
theorem heatBathKernel_eq_conditionalResamplingKernel
    (π : Ω → ℝ) (b : Ω → B) (hπ : ∀ x, 0 < π x) :
    heatBathKernel π b hπ = (toObservation π b hπ).conditionalResamplingKernel := by
  apply FiniteWeightedObservation.finiteMarkovKernel_ext
  exact heatBath_eq_conditionalResamplingKernel π b hπ

/-- Fiber support, fiber-constant rows, and weighted balance recognize heat bath. -/
theorem heatBathKernel_eq_iff_recognition
    (π : Ω → ℝ) (b : Ω → B) (hπ : ∀ x, 0 < π x) (K : FiniteMarkovKernel Ω) :
    K = heatBathKernel π b hπ ↔
      (toObservation π b hπ).FiberSupported K ∧
      (toObservation π b hπ).FiberRowConstant K ∧
      (toObservation π b hπ).WeightedDetailedBalance K := by
  rw [heatBathKernel_eq_conditionalResamplingKernel]
  exact (toObservation π b hπ).kernel_eq_conditionalResamplingKernel_iff_recognition K

/-- Conditional resampling decreases relative entropy for nonnegative input mass. -/
theorem conditionalResamplingKernel_secondLaw (M : FiniteWeightedObservation Ω B)
    (p : Ω → ℝ) (hp : ∀ x, 0 ≤ p x) :
    kl (push p M.conditionalResamplingKernel.probability) M.weight ≤ kl p M.weight := by
  exact heatBath_secondLaw (b := M.observe) M.weight_pos p hp

/-- A point law moves if its observation fiber contains a distinct state. -/
theorem push_single_conditionalResampling_ne_single
    (M : FiniteWeightedObservation Ω B) {x y : Ω} (hxy : x ≠ y)
    (hb : M.observe y = M.observe x) :
    push (Pi.single x (1 : ℝ)) M.conditionalResamplingKernel.probability ≠
      Pi.single x 1 := by
  exact EventAlgebra.push_single_heatBath_ne_single M.weight_pos hxy hb

section Controls

private def flat : FiniteWeightedObservation Bool Unit :=
  ⟨fun _ => 1, fun _ => by norm_num, fun _ => ()⟩

private def tilted : FiniteWeightedObservation Bool Unit :=
  ⟨fun x => if x then 2 else 1, fun x => by cases x <;> norm_num, fun _ => ()⟩

private def separated : FiniteWeightedObservation Bool Bool :=
  ⟨fun _ => 1, fun _ => by norm_num, id⟩

private noncomputable def identityKernel : FiniteMarkovKernel Bool where
  probability x y := if y = x then 1 else 0
  probability_nonneg x y := by split_ifs <;> norm_num
  probability_sum_one x := by cases x <;> simp

private noncomputable def uniformKernel : FiniteMarkovKernel Bool where
  probability _ _ := 1 / 2
  probability_nonneg _ _ := by norm_num
  probability_sum_one _ := by norm_num [Fintype.sum_bool]

/-- Without fiber support, row constancy and balance allow cross-fiber mixing. -/
theorem support_required :
    separated.FiberRowConstant uniformKernel ∧
    separated.WeightedDetailedBalance uniformKernel ∧
    ¬ separated.FiberSupported uniformKernel ∧
    uniformKernel ≠ separated.conditionalResamplingKernel := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro x y _; rfl
  · intro x y; norm_num [separated, uniformKernel]
  · intro h
    have := h (x := false) (y := true) (by decide)
    norm_num [uniformKernel] at this
  · intro h
    have := congrArg (fun K : FiniteMarkovKernel Bool => K.probability false true) h
    norm_num [uniformKernel, FiniteWeightedObservation.conditionalResamplingKernel,
      FiniteWeightedObservation.transition, separated] at this

/-- Without row constancy, support and balance allow the identity kernel. -/
theorem row_constancy_required :
    flat.FiberSupported identityKernel ∧
    flat.WeightedDetailedBalance identityKernel ∧
    ¬ flat.FiberRowConstant identityKernel ∧
    identityKernel ≠ flat.conditionalResamplingKernel := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro x y h; exact (h rfl).elim
  · intro x y; cases x <;> cases y <;> norm_num [flat, identityKernel]
  · intro h
    have := congrFun (h (x := false) (x' := true) rfl) false
    norm_num [identityKernel] at this
  · intro h
    have := congrArg (fun K : FiniteMarkovKernel Bool => K.probability false false) h
    norm_num [identityKernel, FiniteWeightedObservation.conditionalResamplingKernel,
      FiniteWeightedObservation.transition, FiniteWeightedObservation.fiberMass,
      FiniteWeightedObservation.fiber, flat, Fintype.sum_bool] at this

/-- Without weighted balance, support and row constancy allow the wrong weights. -/
theorem balance_required :
    tilted.FiberSupported uniformKernel ∧
    tilted.FiberRowConstant uniformKernel ∧
    ¬ tilted.WeightedDetailedBalance uniformKernel ∧
    uniformKernel ≠ tilted.conditionalResamplingKernel := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro x y h; exact (h rfl).elim
  · intro x y _; rfl
  · intro h
    have := h false true
    norm_num [tilted, uniformKernel] at this
  · intro h
    have := congrArg (fun K : FiniteMarkovKernel Bool => K.probability false false) h
    norm_num [uniformKernel, FiniteWeightedObservation.conditionalResamplingKernel,
      FiniteWeightedObservation.transition, FiniteWeightedObservation.fiberMass,
      FiniteWeightedObservation.fiber, tilted, Fintype.sum_bool] at this

/-- The heat-bath probabilities depend on the prescribed reference weights. -/
theorem reference_weights_matter :
    heatBath tilted.weight tilted.observe ≠ heatBath flat.weight flat.observe := by
  intro h
  have := congrFun (congrFun h false) false
  norm_num [heatBath, fiberMass, tilted, flat, Fintype.sum_bool] at this

/-- Zero weights cannot replace positivity in row normalization. -/
theorem zero_weights_not_stochastic :
    (∑ y : Bool, heatBath (fun _ : Bool => (0 : ℝ)) (fun _ => ()) false y) ≠ 1 := by
  norm_num [heatBath]

/-- Signed input mass can violate the entropy inequality. -/
theorem nonnegative_mass_required :
    ¬ (∀ x : Bool, 0 ≤ (if x then (0 : ℝ) else -1)) ∧
    kl (push (fun x : Bool => if x then 0 else -1)
      flat.conditionalResamplingKernel.probability) flat.weight >
      kl (fun x : Bool => if x then 0 else -1) flat.weight := by
  constructor
  · intro h
    have := h false
    norm_num at this
  · norm_num [kl, klTerm, push, FiniteWeightedObservation.conditionalResamplingKernel,
      FiniteWeightedObservation.transition, FiniteWeightedObservation.fiberMass,
      FiniteWeightedObservation.fiber, flat, Fintype.sum_bool, Real.log_div,
      Real.log_neg_eq_log]
    exact Real.log_pos (by norm_num)

/-- Distinct states alone do not force movement when observation fibers are singletons. -/
theorem singleton_fiber_fixes_point :
    (false : Bool) ≠ true ∧ separated.observe false ≠ separated.observe true ∧
    push (Pi.single false (1 : ℝ)) separated.conditionalResamplingKernel.probability =
      Pi.single false 1 := by
  refine ⟨by decide, by decide, ?_⟩
  funext y
  cases y <;>
    norm_num [push, FiniteWeightedObservation.conditionalResamplingKernel,
      FiniteWeightedObservation.transition, FiniteWeightedObservation.fiberMass,
      FiniteWeightedObservation.fiber, separated, Fintype.sum_bool]
  decide

/-- A singleton state space supplies no distinct same-fiber witness and its point law is fixed. -/
theorem distinct_state_required :
    push (Pi.single () (1 : ℝ)) (heatBath (fun _ : Unit => 1) (fun _ => ())) =
      Pi.single () 1 := by
  funext y
  cases y
  norm_num [push, heatBath, fiberMass]

end Controls

#print axioms heatBath_eq_conditionalResamplingKernel
#print axioms heatBathKernel_eq_conditionalResamplingKernel
#print axioms heatBathKernel_eq_iff_recognition
#print axioms conditionalResamplingKernel_secondLaw
#print axioms push_single_conditionalResampling_ne_single
#print axioms support_required
#print axioms row_constancy_required
#print axioms balance_required
#print axioms reference_weights_matter
#print axioms zero_weights_not_stochastic
#print axioms nonnegative_mass_required
#print axioms singleton_fiber_fixes_point
#print axioms distinct_state_required

end OPH.FiniteConditionalRepairBridge
