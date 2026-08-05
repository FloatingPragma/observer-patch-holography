import Mathlib.Analysis.SpecialFunctions.Complex.CircleAddChar
import EventAlgebra.PublicRecordAlgebra
import QuantumMatterIntegration

/-!
# B10 partition-center readout and finite central action

This module closes three algebraic composition gaps in the supplied B10
interface without adding a source or laboratory identification.

* `partitionCenterAdaptor` is the partition-average map into the commutative
  public-record algebra for the supplied partition.  It factors through block
  pinching but is not identified with block pinching.  An explicit one-block
  control shows that the two maps can differ.
* `mappedCentralAction` is the supplied-weight sixth-root character action on
  the component-labelled finite product of ten nonzero mapped projector
  ranges.  Its kernel is exactly the declared diagonal `ZMod 6` kernel.
* `selectedMappedMatter` attaches the supplied anomaly-free selection mask
  to that same action.  The selected support is one of the two exterior
  parity rows, is action-invariant, and still detects exactly the same
  kernel.

The finite product, partition, component map, weights, and selection mask
remain supplied finite data.  No edge object or edge-to-partition map is part
of the interface.  Nothing here selects the data as physical matter,
constructs a source-produced instrument, derives the action by conjugation on
an ambient module, or identifies the character action with a laboratory gauge
transformation.  Exterior-degree parity is not identified with physical
fermion parity.
-/

namespace OPH.QuantumMatterIntegration

open Matrix

variable {n k : ℕ}

/-! ## The commutative partition-center readout -/

/-- The typed readout into the commutative public-record algebra relative to
the supplied partition.  Its underlying matrix is partition averaging, not
block pinching. -/
noncomputable def partitionCenterAdaptor
    (bridge : FiniteQuantumMatterBridge n k) :
    Matrix (Fin n) (Fin n) ℂ →ₗ[ℂ]
      bridge.partition.publicSubalgebra where
  toFun X := ⟨EventAlgebra.partitionAverage bridge.partition X,
    EventAlgebra.partitionAverage_mem_span bridge.partition X⟩
  map_add' X Y := by
    apply Subtype.ext
    exact EventAlgebra.partitionAverage_add bridge.partition X Y
  map_smul' c X := by
    apply Subtype.ext
    exact EventAlgebra.partitionAverage_smul bridge.partition c X

@[simp]
theorem partitionCenterAdaptor_apply
    (bridge : FiniteQuantumMatterBridge n k)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    (partitionCenterAdaptor bridge X : Matrix (Fin n) (Fin n) ℂ) =
      EventAlgebra.partitionAverage bridge.partition X :=
  rfl

/-- Public readback factors through the declared block readout.  Block
pinching may retain a noncommutative matrix block, after which averaging keeps
only its public central record.  This theorem contains no edge object or
edge-to-partition identification. -/
theorem partitionCenterAdaptor_after_blockReadout
    (bridge : FiniteQuantumMatterBridge n k)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionCenterAdaptor bridge (bridge.blockReadout X) =
      partitionCenterAdaptor bridge X := by
  apply Subtype.ext
  rw [partitionCenterAdaptor_apply, partitionCenterAdaptor_apply,
    bridge.blockReadout_eq_pinching,
    EventAlgebra.partitionPinchingLinearMap_apply]
  exact EventAlgebra.partitionAverage_partitionPinching bridge.partition X

/-- The adaptor is onto the public-record algebra because partition averaging
fixes every public record. -/
theorem partitionCenterAdaptor_surjective
    (bridge : FiniteQuantumMatterBridge n k) :
    Function.Surjective (partitionCenterAdaptor bridge) := by
  intro Y
  refine ⟨Y.1, ?_⟩
  apply Subtype.ext
  exact EventAlgebra.partitionAverage_fixes bridge.partition Y.2

/-- Every adapted public record commutes with every block-preserving
observable.  This is the precise finite center property used here. -/
theorem partitionCenterAdaptor_commutes_with_block
    (bridge : FiniteQuantumMatterBridge n k)
    (X C : Matrix (Fin n) (Fin n) ℂ)
    (hC : C ∈ bridge.partition.commutant) :
    (partitionCenterAdaptor bridge X : Matrix (Fin n) (Fin n) ℂ) * C =
      C * (partitionCenterAdaptor bridge X : Matrix (Fin n) (Fin n) ℂ) := by
  exact bridge.partition.mul_comm_of_mem_span_of_mem_commutant
    (partitionCenterAdaptor bridge X).2 hC

/-! A two-dimensional single-block control makes the distinction between
pinching and public averaging executable inside Lean. -/

/-- The single rank-two block is a projective partition with one outcome. -/
noncomputable def rankTwoSingleBlockPartition :
    EventAlgebra.ProjectivePartition 2 1 where
  proj _ := 1
  isEvent _ := by
    constructor <;> simp [Matrix.IsHermitian]
  orthogonal i j hij := by
    exact False.elim (hij (Subsingleton.elim i j))
  complete := by simp

/-- A nonzero within-block off-diagonal matrix. -/
def rankTwoWithinBlockCoherence : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.single 0 1 1

theorem rankTwoWithinBlockCoherence_ne_zero :
    rankTwoWithinBlockCoherence ≠ 0 := by
  intro h
  have h01 := congrArg (fun M : Matrix (Fin 2) (Fin 2) ℂ ↦ M 0 1) h
  simp [rankTwoWithinBlockCoherence] at h01

/-- Pinching keeps the within-block coherence in the one-block control. -/
theorem rankTwo_pinching_keeps_withinBlock :
    EventAlgebra.partitionPinching rankTwoSingleBlockPartition
        rankTwoWithinBlockCoherence = rankTwoWithinBlockCoherence := by
  exact EventAlgebra.partitionPinching_fixes rankTwoSingleBlockPartition
    (fun _ ↦ by simp [rankTwoSingleBlockPartition])

/-- Public averaging erases the same within-block coherence because only the
scalar record of the rank-two block is public. -/
theorem rankTwo_average_erases_withinBlock :
    EventAlgebra.partitionAverage rankTwoSingleBlockPartition
        rankTwoWithinBlockCoherence = 0 := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    simp [EventAlgebra.partitionAverage, EventAlgebra.bornWeight,
      rankTwoSingleBlockPartition, rankTwoWithinBlockCoherence,
      Matrix.trace]

/-- Mutation control: block pinching and the commutative public readout are
not interchangeable for a nontrivial block. -/
theorem rankTwo_pinching_ne_average :
    EventAlgebra.partitionPinching rankTwoSingleBlockPartition
        rankTwoWithinBlockCoherence ≠
      EventAlgebra.partitionAverage rankTwoSingleBlockPartition
        rankTwoWithinBlockCoherence := by
  rw [rankTwo_pinching_keeps_withinBlock,
    rankTwo_average_erases_withinBlock]
  exact rankTwoWithinBlockCoherence_ne_zero

/-! ## The explicit finite central action -/

/-- The range of the mapped component projector.  Nonzero range is required
by `FiniteQuantumMatterBridge.componentSector_nonzero`. -/
noncomputable def mappedSectorRange
    (bridge : FiniteQuantumMatterBridge n k) (i : Fin 10) :
    Submodule ℂ (Fin n → ℂ) :=
  LinearMap.range
    (Matrix.mulVecLin
      (bridge.partition.proj (bridge.sectorOfComponent i)))

/-- Component-labelled finite dependent product of the ten mapped projector
ranges.  It is a finite direct-sum model, but no explicit direct-sum
equivalence is asserted here.  The component label is retained even when two
rows carry the same supplied weight. -/
noncomputable abbrev MappedMatterModule
    (bridge : FiniteQuantumMatterBridge n k) :=
  (i : Fin 10) → mappedSectorRange bridge i

/-- The standard faithful sixth-root character evaluated on a mapped
component weight. -/
noncomputable def mappedCentralPhase
    (bridge : FiniteQuantumMatterBridge n k)
    (c : OPH.TraceBalancedKernel.C) (i : Fin 10) : ℂ :=
  ZMod.stdAddChar
    (OPH.Z6Descent.phase c
      (bridge.sectorWeight (bridge.sectorOfComponent i)))

/-- The diagonal character action imposed through the supplied sector
weights.  This is not an ambient-projector conjugation action. -/
noncomputable def mappedCentralAction
    (bridge : FiniteQuantumMatterBridge n k)
    (c : OPH.TraceBalancedKernel.C) :
    MappedMatterModule bridge →ₗ[ℂ] MappedMatterModule bridge where
  toFun v i := mappedCentralPhase bridge c i • v i
  map_add' v w := by
    funext i
    exact smul_add _ _ _
  map_smul' a v := by
    funext i
    simp only [Pi.smul_apply, RingHom.id_apply, smul_smul]
    rw [mul_comm]

/-- The central phase is additive in the central parameter. -/
theorem phase_add :
    ∀ (c d : OPH.TraceBalancedKernel.C) (w : OPH.Z6Descent.W),
      OPH.Z6Descent.phase (c + d) w =
        OPH.Z6Descent.phase c w + OPH.Z6Descent.phase d w := by
  decide +kernel

theorem mappedCentralPhase_zero
    (bridge : FiniteQuantumMatterBridge n k) (i : Fin 10) :
    mappedCentralPhase bridge 0 i = 1 := by
  simp [mappedCentralPhase, OPH.Z6Descent.phase]

theorem mappedCentralPhase_add
    (bridge : FiniteQuantumMatterBridge n k)
    (c d : OPH.TraceBalancedKernel.C) (i : Fin 10) :
    mappedCentralPhase bridge (c + d) i =
      mappedCentralPhase bridge c i * mappedCentralPhase bridge d i := by
  rw [mappedCentralPhase, mappedCentralPhase, mappedCentralPhase,
    phase_add]
  exact AddChar.map_add_eq_mul _ _ _

/-- The zero central parameter acts as the identity. -/
theorem mappedCentralAction_zero
    (bridge : FiniteQuantumMatterBridge n k) :
    mappedCentralAction bridge 0 = LinearMap.id := by
  apply LinearMap.ext
  intro v
  funext i
  simp [mappedCentralAction, mappedCentralPhase_zero]

/-- Addition of central parameters composes their actions. -/
theorem mappedCentralAction_add
    (bridge : FiniteQuantumMatterBridge n k)
    (c d : OPH.TraceBalancedKernel.C) :
    mappedCentralAction bridge (c + d) =
      (mappedCentralAction bridge c).comp
        (mappedCentralAction bridge d) := by
  apply LinearMap.ext
  intro v
  funext i
  simp only [mappedCentralAction, LinearMap.coe_mk, AddHom.coe_mk,
    LinearMap.comp_apply, mappedCentralPhase_add]
  exact SemigroupAction.mul_smul _ _ _

/-- The action of the additive inverse is a two-sided linear inverse. -/
theorem mappedCentralAction_neg_comp
    (bridge : FiniteQuantumMatterBridge n k)
    (c : OPH.TraceBalancedKernel.C) :
    (mappedCentralAction bridge (-c)).comp
        (mappedCentralAction bridge c) = LinearMap.id ∧
      (mappedCentralAction bridge c).comp
        (mappedCentralAction bridge (-c)) = LinearMap.id := by
  constructor
  · rw [← mappedCentralAction_add, neg_add_cancel,
      mappedCentralAction_zero]
  · rw [← mappedCentralAction_add, add_neg_cancel,
      mappedCentralAction_zero]

/-- Every mapped projector range contains a nonzero vector. -/
theorem mappedSectorRange_nontrivial
    (bridge : FiniteQuantumMatterBridge n k) (i : Fin 10) :
    ∃ v : mappedSectorRange bridge i, v ≠ 0 := by
  let P := bridge.partition.proj (bridge.sectorOfComponent i)
  have hP : P ≠ 0 := bridge.componentSector_nonzero i
  obtain ⟨a, b, hab⟩ : ∃ a b, P a b ≠ 0 := by
    by_contra h
    push Not at h
    apply hP
    ext x y
    exact h x y
  let e : Fin n → ℂ := Pi.single b 1
  let y : Fin n → ℂ := P *ᵥ e
  have hya : y a = P a b := by
    simp [y, e]
  have hy : y ≠ 0 := by
    intro hzero
    have ha := congrFun hzero a
    rw [hya] at ha
    exact hab ha
  refine ⟨⟨y, ?_⟩, ?_⟩
  · exact ⟨e, rfl⟩
  · intro hzero
    apply hy
    exact congrArg Subtype.val hzero

/-- Operator identity is equivalent to zero phase on every mapped component.
This step uses only faithfulness of the standard sixth-root character and
nonzero mapped projector ranges; it does not use the component kernel table. -/
theorem mappedCentralAction_eq_id_iff_component_phases_zero
    (bridge : FiniteQuantumMatterBridge n k)
    (c : OPH.TraceBalancedKernel.C) :
    mappedCentralAction bridge c = LinearMap.id ↔
      ∀ i : Fin 10,
        OPH.Z6Descent.phase c
          (bridge.sectorWeight (bridge.sectorOfComponent i)) = 0 := by
  constructor
  · intro hAction i
    obtain ⟨v, hv⟩ := mappedSectorRange_nontrivial bridge i
    let x : MappedMatterModule bridge := Pi.single i v
    have hx := congrArg
      (fun F : MappedMatterModule bridge →ₗ[ℂ] MappedMatterModule bridge ↦
        (F x) i) hAction
    have hsmul : mappedCentralPhase bridge c i • v = (1 : ℂ) • v := by
      simpa [mappedCentralAction, x] using hx
    have hphaseScalar : mappedCentralPhase bridge c i = 1 :=
      smul_left_injective ℂ hv hsmul
    apply ZMod.injective_stdAddChar
    simpa [mappedCentralPhase] using hphaseScalar
  · intro hPhases
    apply LinearMap.ext
    intro v
    funext i
    have hphase := hPhases i
    simp [mappedCentralAction, mappedCentralPhase, hphase]

/-- The actual action kernel on the mapped projector-range module is exactly
the declared six-element tensor-character kernel.  The proof composes the
independent componentwise phase-detection lemma above with the kernel-checked
exterior weight arithmetic. -/
theorem mappedCentralAction_eq_id_iff
    (bridge : FiniteQuantumMatterBridge n k)
    (c : OPH.TraceBalancedKernel.C) :
    mappedCentralAction bridge c = LinearMap.id ↔
      OPH.TraceBalancedKernel.tensorCharFun c = 0 :=
  (mappedCentralAction_eq_id_iff_component_phases_zero bridge c).trans
    (kernel_on_mapped_component_weights bridge c)

/-- The kernel of the central-parameter action has six elements. -/
theorem mappedCentralAction_kernel_card
    (bridge : FiniteQuantumMatterBridge n k) :
    Set.ncard {c : OPH.TraceBalancedKernel.C |
      mappedCentralAction bridge c = LinearMap.id} = 6 := by
  classical
  have hset :
      {c : OPH.TraceBalancedKernel.C |
        mappedCentralAction bridge c = LinearMap.id} =
      {c : OPH.TraceBalancedKernel.C |
        OPH.TraceBalancedKernel.tensorCharFun c = 0} := by
    ext c
    exact mappedCentralAction_eq_id_iff bridge c
  rw [hset]
  rw [Set.ncard_eq_toFinset_card]
  simpa only [Set.Finite.toFinset_setOf] using
    OPH.TraceBalancedKernel.kernel_card

/-! ## Attaching the anomaly-free parity selection -/

/-- The selected component support inside the mapped matter module. -/
def selectedMappedMatter
    (bridge : FiniteQuantumMatterBridge n k) :
    Submodule ℂ (MappedMatterModule bridge) where
  carrier := {v | ∀ i : Fin 10,
    OPH.ExteriorSelection.mem bridge.selectionMask.val i = false →
      v i = 0}
  zero_mem' := by simp
  add_mem' := by
    intro v w hv hw i hi
    simp [hv i hi, hw i hi]
  smul_mem' := by
    intro a v hv i hi
    simp [hv i hi]

/-- The selected support is one of the two exterior-degree parity rows. -/
theorem selectedMappedMatter_support_is_parity
    (bridge : FiniteQuantumMatterBridge n k) :
    (∀ i : Fin 10,
      OPH.ExteriorSelection.mem bridge.selectionMask.val i =
        !OPH.ExteriorSelection.odd i) ∨
    (∀ i : Fin 10,
      OPH.ExteriorSelection.mem bridge.selectionMask.val i =
        OPH.ExteriorSelection.odd i) := by
  rcases bridge_selection_is_parity_sector bridge with hEven | hOdd
  · left
    intro i
    rw [hEven]
    exact OPH.ExteriorSelection.evenMask_is_even_sector i
  · right
    intro i
    rw [hOdd]
    exact OPH.ExteriorSelection.oddMask_is_odd_sector i

/-- The diagonal central action preserves the attached selected support. -/
theorem mappedCentralAction_preserves_selected
    (bridge : FiniteQuantumMatterBridge n k)
    (c : OPH.TraceBalancedKernel.C)
    {v : MappedMatterModule bridge}
    (hv : v ∈ selectedMappedMatter bridge) :
    mappedCentralAction bridge c v ∈ selectedMappedMatter bridge := by
  intro i hi
  simp [mappedCentralAction, hv i hi]

/-- A selected component with a nonzero projector supplies a nonzero vector
in the selected submodule. -/
theorem selectedMappedMatter_nontrivial_of_mem
    (bridge : FiniteQuantumMatterBridge n k)
    (i : Fin 10)
    (hi : OPH.ExteriorSelection.mem bridge.selectionMask.val i = true) :
    ∃ v : MappedMatterModule bridge,
      v ∈ selectedMappedMatter bridge ∧ v ≠ 0 := by
  obtain ⟨u, hu⟩ := mappedSectorRange_nontrivial bridge i
  let v : MappedMatterModule bridge := Pi.single i u
  refine ⟨v, ?_, ?_⟩
  · intro j hj
    by_cases hji : j = i
    · subst j
      simp [hi] at hj
    · simp [v, hji]
  · intro hv
    have hvi := congrArg (fun x : MappedMatterModule bridge ↦ x i) hv
    exact hu (by simpa [v] using hvi)

/-- The attached selected support is inhabited for every bridge satisfying
the exhaustive chiral anomaly-free selection contract. -/
theorem selectedMappedMatter_nontrivial
    (bridge : FiniteQuantumMatterBridge n k) :
    ∃ v : MappedMatterModule bridge,
      v ∈ selectedMappedMatter bridge ∧ v ≠ 0 := by
  rcases bridge_selection_is_parity_sector bridge with hEven | hOdd
  · apply selectedMappedMatter_nontrivial_of_mem bridge (2 : Fin 10)
    simpa only [hEven] using
      (show OPH.ExteriorSelection.mem OPH.ExteriorSelection.evenMask
          (2 : Fin 10) = true by decide)
  · apply selectedMappedMatter_nontrivial_of_mem bridge (0 : Fin 10)
    simpa only [hOdd] using
      (show OPH.ExteriorSelection.mem OPH.ExteriorSelection.oddMask
          (0 : Fin 10) = true by decide)

theorem kernel_on_even_selected_components :
    ∀ c : OPH.TraceBalancedKernel.C,
      (∀ i : Fin 10,
        OPH.ExteriorSelection.mem OPH.ExteriorSelection.evenMask i = true →
        OPH.Z6Descent.phase c (componentWeight i) = 0) ↔
      OPH.TraceBalancedKernel.tensorCharFun c = 0 := by
  decide +kernel

theorem kernel_on_odd_selected_components :
    ∀ c : OPH.TraceBalancedKernel.C,
      (∀ i : Fin 10,
        OPH.ExteriorSelection.mem OPH.ExteriorSelection.oddMask i = true →
        OPH.Z6Descent.phase c (componentWeight i) = 0) ↔
      OPH.TraceBalancedKernel.tensorCharFun c = 0 := by
  decide +kernel

theorem no_even_selected_universal_minus_one :
    ∀ c : OPH.TraceBalancedKernel.C,
      ¬ (∀ i : Fin 10,
        OPH.ExteriorSelection.mem OPH.ExteriorSelection.evenMask i = true →
        OPH.Z6Descent.phase c (componentWeight i) = 3) := by
  decide +kernel

theorem no_odd_selected_universal_minus_one :
    ∀ c : OPH.TraceBalancedKernel.C,
      ¬ (∀ i : Fin 10,
        OPH.ExteriorSelection.mem OPH.ExteriorSelection.oddMask i = true →
        OPH.Z6Descent.phase c (componentWeight i) = 3) := by
  decide +kernel

/-- Either admissible parity row detects exactly the same central kernel
after attachment to the mapped-sector action. -/
theorem kernel_on_selected_mapped_components
    (bridge : FiniteQuantumMatterBridge n k)
    (c : OPH.TraceBalancedKernel.C) :
    (∀ i : Fin 10,
      OPH.ExteriorSelection.mem bridge.selectionMask.val i = true →
      OPH.Z6Descent.phase c
        (bridge.sectorWeight (bridge.sectorOfComponent i)) = 0) ↔
      OPH.TraceBalancedKernel.tensorCharFun c = 0 := by
  rcases bridge_selection_is_parity_sector bridge with hEven | hOdd
  · simpa only [hEven, bridge.sectorWeight_component] using
      kernel_on_even_selected_components c
  · simpa only [hOdd, bridge.sectorWeight_component] using
      kernel_on_odd_selected_components c

/-- Attaching the parity support does not turn gauge-center arithmetic into
fermion parity: no central parameter acts by minus one on every selected
component in either admissible row. -/
theorem no_selected_central_parameter_realizes_parity_sign
    (bridge : FiniteQuantumMatterBridge n k) :
    ∀ c : OPH.TraceBalancedKernel.C,
      ¬ (∀ i : Fin 10,
        OPH.ExteriorSelection.mem bridge.selectionMask.val i = true →
        OPH.Z6Descent.phase c
          (bridge.sectorWeight (bridge.sectorOfComponent i)) = 3) := by
  intro c
  rcases bridge_selection_is_parity_sector bridge with hEven | hOdd
  · simpa only [hEven, bridge.sectorWeight_component] using
      no_even_selected_universal_minus_one c
  · simpa only [hOdd, bridge.sectorWeight_component] using
      no_odd_selected_universal_minus_one c

/-- Exact kernel statement for the action restricted by the attached parity
mask.  It quantifies over the selected submodule rather than silently
identifying the selection mask with a separate table. -/
theorem selectedMappedCentralAction_eq_id_iff
    (bridge : FiniteQuantumMatterBridge n k)
    (c : OPH.TraceBalancedKernel.C) :
    (∀ v : MappedMatterModule bridge,
      v ∈ selectedMappedMatter bridge →
        mappedCentralAction bridge c v = v) ↔
      OPH.TraceBalancedKernel.tensorCharFun c = 0 := by
  constructor
  · intro hAction
    apply (kernel_on_selected_mapped_components bridge c).mp
    intro i hi
    obtain ⟨v, hv⟩ := mappedSectorRange_nontrivial bridge i
    let x : MappedMatterModule bridge := Pi.single i v
    have hxmem : x ∈ selectedMappedMatter bridge := by
      intro j hj
      by_cases hji : j = i
      · subst j
        simp [hi] at hj
      · simp [x, hji]
    have hx := congrArg (fun y : MappedMatterModule bridge ↦ y i)
      (hAction x hxmem)
    have hsmul : mappedCentralPhase bridge c i • v = (1 : ℂ) • v := by
      simpa [mappedCentralAction, x] using hx
    have hphaseScalar : mappedCentralPhase bridge c i = 1 :=
      smul_left_injective ℂ hv hsmul
    apply ZMod.injective_stdAddChar
    simpa [mappedCentralPhase] using hphaseScalar
  · intro hKernel v hv
    funext i
    by_cases hi : OPH.ExteriorSelection.mem bridge.selectionMask.val i = true
    · have hphase :=
        (kernel_on_selected_mapped_components bridge c).mpr hKernel i hi
      simp [mappedCentralAction, mappedCentralPhase, hphase]
    · have hiFalse :
          OPH.ExteriorSelection.mem bridge.selectionMask.val i = false :=
        Bool.eq_false_of_not_eq_true hi
      simp [mappedCentralAction, hv i hiFalse]

#print axioms partitionCenterAdaptor_after_blockReadout
#print axioms partitionCenterAdaptor_surjective
#print axioms partitionCenterAdaptor_commutes_with_block
#print axioms rankTwo_pinching_ne_average
#print axioms mappedCentralAction_zero
#print axioms mappedCentralAction_add
#print axioms mappedCentralAction_neg_comp
#print axioms mappedCentralAction_eq_id_iff_component_phases_zero
#print axioms mappedCentralAction_eq_id_iff
#print axioms mappedCentralAction_kernel_card
#print axioms selectedMappedMatter_support_is_parity
#print axioms mappedCentralAction_preserves_selected
#print axioms selectedMappedMatter_nontrivial
#print axioms kernel_on_selected_mapped_components
#print axioms no_selected_central_parameter_realizes_parity_sign
#print axioms selectedMappedCentralAction_eq_id_iff

end OPH.QuantumMatterIntegration
