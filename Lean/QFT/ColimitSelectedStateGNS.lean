import QFT.ColimitGNSRepresentation
import EventAlgebra.ExpectationBound

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# A selected-state GNS representation of the completed observer colimit

A `ConsensusTower` already carries one normalized finite density matrix for
each observer, and `CoherentObserverFamily` selects compatible observers
through the regulator system.  This module uses exactly those two committed
inputs to descend the finite trace expectations to the filtered colimit,
extend the resulting bounded functional to the C*-completion, and feed it to
the conditional GNS adapter.

The construction does not select a coherent observer family from a bare
tower.  It does not prove faithfulness, construct a physical vacuum, attach
spacetime or dynamics, or establish a physical QFT, ground-state condition,
or time-slice theorem.  The concrete witness below is the declared finite
`witnessTower` receipt, whose regulator and observer types are both `Unit`.
-/

namespace OPH.QFT

open OPH.Tower
open UniformSpace

open scoped ComplexOrder InnerProductSpace Matrix.Norms.L2Operator
  OPH.QFT.CompletionStar

universe u

variable {ι : Type u} [Preorder ι]
variable {T : ConsensusTower ι}

/-! ## The compatible finite expectation on the algebraic colimit -/

/-- The trace expectation of the observer selected by `ξ` on one tower germ. -/
private noncomputable def selectedGermExpectation
    (ξ : CoherentObserverFamily T) (p : TowerGerm T) : ℂ :=
  EventAlgebra.expectation (T.state p.1 (ξ.observer p.1)) p.2

/-- State naturality and observer coherence make the selected expectation
constant on the germ relation. -/
private theorem selectedGermExpectation_sound
    (ξ : CoherentObserverFamily T) {p q : TowerGerm T}
    (h : germRel T p q) :
    selectedGermExpectation ξ p = selectedGermExpectation ξ q := by
  obtain ⟨t, hp, hq, e⟩ := h
  have hp_state := T.state_natural hp (ξ.observer p.1) p.2
  have hq_state := T.state_natural hq (ξ.observer q.1) q.2
  rw [ξ.coherent hp] at hp_state
  rw [ξ.coherent hq] at hq_state
  change
    (T.state p.1 (ξ.observer p.1) * p.2).trace =
      (T.state q.1 (ξ.observer q.1) * q.2).trace
  calc
    (T.state p.1 (ξ.observer p.1) * p.2).trace =
        (T.state t (ξ.observer t) * T.algebraRefine hp p.2).trace :=
      hp_state.symm
    _ = (T.state t (ξ.observer t) * T.algebraRefine hq q.2).trace := by rw [e]
    _ = (T.state q.1 (ξ.observer q.1) * q.2).trace := hq_state

/-- The compatible selected state, descended to the algebraic filtered
colimit as a complex-linear functional. -/
noncomputable def selectedColimitExpectationLinearMap
    (ξ : CoherentObserverFamily T) : TowerColimit T →ₗ[ℂ] ℂ where
  toFun := Quotient.lift (selectedGermExpectation ξ)
    (fun _ _ h => selectedGermExpectation_sound ξ h)
  map_add' x y := by
    obtain ⟨r, X, Y, rfl, rfl⟩ := exists_rep₂ x y
    rw [colimitMk_add]
    exact map_add (EventAlgebra.stateExpectationLinearMap
      (T.state r (ξ.observer r))) X Y
  map_smul' c x := by
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    rw [colimitMk_smul]
    exact map_smul (EventAlgebra.stateExpectationLinearMap
      (T.state r (ξ.observer r))) c X

/-- On a stage image the descended functional is the selected finite trace
expectation, with no limiting ambiguity. -/
@[simp]
theorem selectedColimitExpectationLinearMap_colimitMk
    (ξ : CoherentObserverFamily T) (r : ι)
    (X : ConsensusTower.PrivateAlgebra T r) :
    selectedColimitExpectationLinearMap ξ (colimitMk X) =
      EventAlgebra.expectation (T.state r (ξ.observer r)) X :=
  rfl

/-- The descended selected-state expectation has operator norm at most one
for the committed colimit norm. -/
theorem norm_selectedColimitExpectationLinearMap_le
    (ξ : CoherentObserverFamily T) (x : TowerColimit T) :
    ‖selectedColimitExpectationLinearMap ξ x‖ ≤ ‖x‖ := by
  obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
  obtain ⟨s, hrs, hs⟩ := exists_stage_norm_eq X
  letI : NeZero (T.dim s) :=
    ⟨Nat.ne_of_gt (dim_pos_of_observer T s ⟨ξ.observer s⟩)⟩
  have hstate := T.state_natural hrs (ξ.observer r) X
  rw [ξ.coherent hrs] at hstate
  calc
    ‖selectedColimitExpectationLinearMap ξ (colimitMk X)‖ =
        ‖EventAlgebra.expectation
          (T.state s (ξ.observer s)) (T.algebraRefine hrs X)‖ := by
      rw [selectedColimitExpectationLinearMap_colimitMk]
      exact congrArg norm hstate.symm
    _ ≤ ‖T.algebraRefine hrs X‖ :=
      EventAlgebra.norm_expectation_le_l2_opNorm
        (T.state_isState s (ξ.observer s)) _
    _ = colimitNorm (colimitMk X) := hs
    _ = ‖(colimitMk X : TowerColimit T)‖ := rfl

/-- The descended functional as a continuous linear map. -/
noncomputable def selectedColimitExpectation
    (ξ : CoherentObserverFamily T) : TowerColimit T →L[ℂ] ℂ :=
  (selectedColimitExpectationLinearMap ξ).mkContinuous 1 fun x => by
    simpa using norm_selectedColimitExpectationLinearMap_le ξ x

@[simp]
theorem selectedColimitExpectation_colimitMk
    (ξ : CoherentObserverFamily T) (r : ι)
    (X : ConsensusTower.PrivateAlgebra T r) :
    selectedColimitExpectation ξ (colimitMk X) =
      EventAlgebra.expectation (T.state r (ξ.observer r)) X :=
  rfl

/-! ## Extension to the completed colimit -/

/-- The uniformly continuous extension of the selected colimit expectation
to the completed colimit. -/
private noncomputable def selectedCompletionExpectationFun
    (ξ : CoherentObserverFamily T) : ColimitCompletion T → ℂ :=
  Completion.extension (selectedColimitExpectation ξ)

private theorem selectedCompletionExpectationFun_coe
    (ξ : CoherentObserverFamily T) (x : TowerColimit T) :
    selectedCompletionExpectationFun ξ (x : ColimitCompletion T) =
      selectedColimitExpectation ξ x := by
  exact Completion.extension_coe (selectedColimitExpectation ξ).uniformContinuous x

/-- The completion extension remains complex-linear. -/
noncomputable def selectedCompletionExpectationLinearMap
    (ξ : CoherentObserverFamily T) : ColimitCompletion T →ₗ[ℂ] ℂ where
  toFun := selectedCompletionExpectationFun ξ
  map_add' x y := by
    induction x, y using Completion.induction_on₂ with
    | hp =>
      exact isClosed_eq
        (Completion.continuous_extension.comp continuous_add)
        ((Completion.continuous_extension.comp continuous_fst).add
          (Completion.continuous_extension.comp continuous_snd))
    | ih a b =>
      rw [← Completion.coe_add]
      simp only [selectedCompletionExpectationFun_coe]
      exact map_add (selectedColimitExpectation ξ) a b
  map_smul' c x := by
    induction x using Completion.induction_on with
    | hp =>
      exact isClosed_eq
        (Completion.continuous_extension.comp (continuous_const_smul c))
        ((continuous_const_smul c).comp Completion.continuous_extension)
    | ih a =>
      rw [← Completion.coe_smul]
      simp only [selectedCompletionExpectationFun_coe]
      exact map_smul (selectedColimitExpectation ξ) c a

/-- The completion extension evaluates on every stage exactly as the finite
selected density state. -/
@[simp]
theorem selectedCompletionExpectationLinearMap_stage
    (ξ : CoherentObserverFamily T) (r : ι)
    (X : ConsensusTower.PrivateAlgebra T r) :
    selectedCompletionExpectationLinearMap ξ (stageToCompletion T r X) =
      EventAlgebra.expectation (T.state r (ξ.observer r)) X := by
  rw [stageToCompletion_apply]
  exact selectedCompletionExpectationFun_coe ξ (colimitMk X)

/- The completion order is deliberately local, matching the conditional GNS
adapter and preventing a spectral-order instance from leaking globally. -/
noncomputable local instance selectedCompletionPartialOrder :
    PartialOrder (ColimitCompletion T) :=
  CStarAlgebra.spectralOrder (ColimitCompletion T)

local instance selectedCompletionStarOrderedRing :
    StarOrderedRing (ColimitCompletion T) :=
  CStarAlgebra.spectralOrderedRing (ColimitCompletion T)

/-- The selected expectation is nonnegative on algebraic star-squares. -/
private theorem selectedColimitExpectation_star_mul_self_nonneg
    (ξ : CoherentObserverFamily T) (x : TowerColimit T) :
    0 ≤ selectedColimitExpectation ξ (star x * x) := by
  obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
  rw [colimitMk_star, colimitMk_mul, selectedColimitExpectation_colimitMk]
  apply EventAlgebra.expectation_nonneg (T.state_isState r (ξ.observer r)).1
  rw [Matrix.star_eq_conjTranspose]
  exact Matrix.posSemidef_conjTranspose_mul_self X

/-- Positivity of algebraic star-squares extends to every star-square in the
completion by density and continuity. -/
private theorem selectedCompletionExpectation_star_mul_self_nonneg
    (ξ : CoherentObserverFamily T) (x : ColimitCompletion T) :
    0 ≤ selectedCompletionExpectationLinearMap ξ (star x * x) := by
  induction x using Completion.induction_on with
  | hp =>
    exact isClosed_le continuous_const
      (Completion.continuous_extension.comp
        (continuous_mul.comp
          (completion_map_star_continuous.prodMk continuous_id)))
  | ih a =>
    rw [completion_star_coe, ← Completion.coe_mul]
    change 0 ≤ selectedCompletionExpectationFun ξ
      ((star a * a : TowerColimit T) : ColimitCompletion T)
    rw [selectedCompletionExpectationFun_coe]
    exact selectedColimitExpectation_star_mul_self_nonneg ξ a

/-- The normalized positive functional canonically induced by the states
along the explicitly supplied coherent observer family. -/
noncomputable def selectedColimitFunctional
    (ξ : CoherentObserverFamily T) : ColimitCompletion T →ₚ[ℂ] ℂ :=
  PositiveLinearMap.mk₀ (selectedCompletionExpectationLinearMap ξ) fun x hx => by
    obtain ⟨y, rfl⟩ := CStarAlgebra.nonneg_iff_eq_star_mul_self.mp hx
    exact selectedCompletionExpectation_star_mul_self_nonneg ξ y

/-- Stage evaluation for the selected normalized positive functional. -/
@[simp]
theorem selectedColimitFunctional_stage
    (ξ : CoherentObserverFamily T) (r : ι)
    (X : ConsensusTower.PrivateAlgebra T r) :
    selectedColimitFunctional ξ (stageToCompletion T r X) =
      EventAlgebra.expectation (T.state r (ξ.observer r)) X :=
  selectedCompletionExpectationLinearMap_stage ξ r X

/-- The selected positive functional is normalized. -/
theorem selectedColimitFunctional_normalized
    (ξ : CoherentObserverFamily T) :
    IsNormalizedColimitFunctional (selectedColimitFunctional ξ) := by
  change selectedColimitFunctional ξ (1 : ColimitCompletion T) = 1
  rw [← map_one (stageToCompletion T (towerBase T))]
  rw [selectedColimitFunctional_stage]
  exact EventAlgebra.expectation_one
    (T.state_isState (towerBase T) (ξ.observer (towerBase T)))

/-! ## The committed finite witness -/

/-- The actual GNS Hilbert space induced by the committed selected density
state of `witnessTower`.  This is a finite declared witness, not a physical
vacuum sector or a spacetime QFT. -/
abbrev WitnessColimitGNS :=
  ColimitGNS (selectedColimitFunctional witnessObserverFamily)

/-- The witness GNS unit class is normalized by the constructed functional. -/
theorem norm_witnessColimitGNSUnitClass :
    ‖colimitGNSUnitClass (selectedColimitFunctional witnessObserverFamily)‖ = 1 :=
  norm_colimitGNSUnitClass _
    (selectedColimitFunctional_normalized witnessObserverFamily)

/-- The constructed witness Hilbert carrier is nonzero: its distinguished
unit class has norm one.  This is mathematical nontriviality, not a physical
vacuum claim. -/
theorem witnessColimitGNSUnitClass_ne_zero :
    colimitGNSUnitClass
      (selectedColimitFunctional witnessObserverFamily) ≠ 0 := by
  intro h
  have hnorm := norm_witnessColimitGNSUnitClass
  rw [h, norm_zero] at hnorm
  exact zero_ne_one hnorm

/-- The committed separating projector remains a norm-one element of the
source completion underlying the witness GNS construction.  No faithfulness
of its represented image is asserted. -/
theorem witnessColimitGNS_source_nontrivial_receipt :
    ‖stageToCompletion witnessTower () e0‖ = 1 :=
  witnessTower_norm_stageToCompletion_e0

end OPH.QFT

-- Axiom audit: all declarations must remain on the standard Mathlib basis.
#print axioms OPH.QFT.selectedColimitExpectationLinearMap_colimitMk
#print axioms OPH.QFT.norm_selectedColimitExpectationLinearMap_le
#print axioms OPH.QFT.selectedCompletionExpectationLinearMap_stage
#print axioms OPH.QFT.selectedColimitFunctional
#print axioms OPH.QFT.selectedColimitFunctional_stage
#print axioms OPH.QFT.selectedColimitFunctional_normalized
#print axioms OPH.QFT.norm_witnessColimitGNSUnitClass
#print axioms OPH.QFT.witnessColimitGNSUnitClass_ne_zero
#print axioms OPH.QFT.witnessColimitGNS_source_nontrivial_receipt
