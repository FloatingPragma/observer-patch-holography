import EventAlgebra.Superselection
import ExteriorComponentBridge
import Z6Descent

/-!
# Typed integration of finite sectors, exterior components, and Z6 descent

The existing finite theorems concern three different supplied objects: a
projective sector partition, a five-mode exterior algebra, and a table of
central characters.  They can be composed only after explicit maps identify
their labels.  `FiniteQuantumMatterBridge` records exactly those maps and the
declared block-preserving readout.  The commutative public record algebra is
the smaller partition span and its canonical readout is partition averaging;
neither is identified with the block pinching below.  The theorems transport
the operational block quotient and the exact weight-kernel arithmetic through
the maps rather than conjoining unrelated conclusions.

The bridge is an interface, not a source theorem.  In particular, it does not
derive a physical sector partition, select the exterior carrier as matter,
select one parity representative, select a physical global gauge-group form,
or identify any laboratory charge.
-/

namespace OPH.QuantumMatterIntegration

open Matrix

/-- Central weight canonically read from an exterior component row.  Color
triality is its color exterior degree modulo three, weak duality is its weak
degree modulo two, and the first coordinate is the frozen integer charge
modulo six. -/
def componentWeight (i : Fin 10) : OPH.Z6Descent.W :=
  ((OPH.ExteriorSelection.charge i : ZMod 6),
    ((OPH.ExteriorComponentBridge.componentDegree i).1.val : ZMod 3),
    ((OPH.ExteriorComponentBridge.componentDegree i).2.val : ZMod 2))

/-- In Standard-Model field order, the even exterior rows reproduce the five
matter weights used by the Z6 descent theorem. -/
theorem even_component_weights_eq_matterWeights :
    [componentWeight 3, componentWeight 2, componentWeight 4,
      componentWeight 8, componentWeight 9] =
      OPH.Z6Descent.matterWeights := by
  decide +kernel

/-- The common central stabilizer of all ten exterior component rows is
exactly the six-element tensor-character kernel.  Including conjugate rows
does not enlarge or shrink it. -/
theorem kernel_on_exterior_component_weights :
    ∀ c : OPH.TraceBalancedKernel.C,
      (∀ i : Fin 10, OPH.Z6Descent.phase c (componentWeight i) = 0) ↔
        OPH.TraceBalancedKernel.tensorCharFun c = 0 := by
  decide +kernel

/-- Mutation control: adjoining one fractionally charged singlet to the ten
exterior rows collapses their common central kernel to the identity. -/
theorem fractional_singlet_mutation_collapses_component_kernel :
    (Finset.univ.filter fun c : OPH.TraceBalancedKernel.C ↦
      (∀ i : Fin 10,
        OPH.Z6Descent.phase c (componentWeight i) = 0) ∧
      OPH.Z6Descent.phase c (1, 0, 0) = 0).card = 1 := by
  decide +kernel

/-- The explicit data required to identify partition sectors with the ten
exterior-component rows and their central weights.  Mapped projectors must be
nonzero.  The map is not assumed injective: rows carrying the same operational
sector may coincide, while `sectorWeight_component` prevents rows with
different central weights from being silently identified. -/
structure FiniteQuantumMatterBridge (n k : ℕ) where
  partition : EventAlgebra.ProjectivePartition n k
  blockReadout :
    Matrix (Fin n) (Fin n) ℂ →ₗ[ℂ] Matrix (Fin n) (Fin n) ℂ
  blockReadout_eq_pinching :
    blockReadout = EventAlgebra.partitionPinchingLinearMap partition
  sectorOfComponent : Fin 10 → Fin k
  componentSector_nonzero :
    ∀ i : Fin 10, partition.proj (sectorOfComponent i) ≠ 0
  sectorWeight : Fin k → OPH.Z6Descent.W
  sectorWeight_component :
    ∀ i : Fin 10, sectorWeight (sectorOfComponent i) = componentWeight i
  selectionMask : Fin 1024
  selection_nonempty : selectionMask.val ≠ 0
  selection_chiral :
    OPH.ExteriorSelection.chiral selectionMask.val = true
  selection_anomalyFree :
    OPH.ExteriorSelection.anomalyFree selectionMask.val = true

/-! ## Exact synthetic nonvacuity control -/

/-- A coordinate projector used only to show that the bridge interface is
mathematically inhabited. -/
def coordinateProjector (i : Fin 10) : Matrix (Fin 10) (Fin 10) ℂ :=
  Matrix.single i i 1

/-- Every coordinate projector in the control is nonzero. -/
theorem coordinateProjector_ne_zero (i : Fin 10) :
    coordinateProjector i ≠ 0 := by
  intro h
  have hii := congrArg (fun M : Matrix (Fin 10) (Fin 10) ℂ ↦ M i i) h
  simp [coordinateProjector] at hii

/-- The ten coordinate lines form an exact projective partition. -/
noncomputable def coordinatePartition :
    EventAlgebra.ProjectivePartition 10 10 where
  proj := coordinateProjector
  isEvent i := by
    constructor
    · rw [Matrix.IsHermitian]
      simp [coordinateProjector]
    · simp [coordinateProjector]
  orthogonal i j hij := by
    simpa [coordinateProjector] using
      (Matrix.single_mul_single_of_ne (c := (1 : ℂ)) i i j hij (1 : ℂ))
  complete := by
    simpa [coordinateProjector] using
      (Matrix.sum_single_one (m := Fin 10) (α := ℂ))

/-- A synthetic ten-sector bridge built from coordinate lines.  It is a
nonvacuity control for the interface, not a source or physical realization. -/
noncomputable def coordinateBridgeControl :
    FiniteQuantumMatterBridge 10 10 where
  partition := coordinatePartition
  blockReadout := EventAlgebra.partitionPinchingLinearMap coordinatePartition
  blockReadout_eq_pinching := rfl
  sectorOfComponent := id
  componentSector_nonzero := coordinateProjector_ne_zero
  sectorWeight := componentWeight
  sectorWeight_component _ := rfl
  selectionMask := ⟨OPH.ExteriorSelection.evenMask, by decide⟩
  selection_nonempty := by decide
  selection_chiral := OPH.ExteriorSelection.parity_sectors_survive.1
  selection_anomalyFree := OPH.ExteriorSelection.parity_sectors_survive.2.1

/-- Mutation control for the readout boundary: a sector diagonal is retained
by pinching and therefore is not operationally invisible. -/
theorem coordinate_diagonal_not_partitionOffDiagonal (i : Fin 10) :
    ¬ EventAlgebra.IsPartitionOffDiagonal coordinatePartition
      (coordinateProjector i) := by
  intro h
  unfold EventAlgebra.IsPartitionOffDiagonal at h
  have hfix :
      EventAlgebra.partitionPinching coordinatePartition
          (coordinateProjector i) = coordinateProjector i := by
    exact EventAlgebra.partitionPinching_fixes coordinatePartition
      (fun j ↦ coordinatePartition.proj_commute i j)
  rw [hfix] at h
  exact coordinateProjector_ne_zero i h

/-- A nonzero cross-sector coherence for the synthetic control. -/
def coordinateCoherence : Matrix (Fin 10) (Fin 10) ℂ :=
  Matrix.single 0 1 1

/-- The operational-invisibility kernel of the control is nontrivial: it
contains this explicit nonzero cross-sector coherence. -/
theorem coordinate_nonzero_offDiagonal_control :
    coordinateCoherence ≠ 0 ∧
      EventAlgebra.IsPartitionOffDiagonal coordinatePartition
        coordinateCoherence := by
  constructor
  · intro h
    have h01 := congrArg
      (fun M : Matrix (Fin 10) (Fin 10) ℂ ↦ M 0 1) h
    simp [coordinateCoherence] at h01
  · have hcorner :
        EventAlgebra.partitionCorner coordinatePartition 0 1
            coordinateCoherence = coordinateCoherence := by
      simp [EventAlgebra.partitionCorner, coordinatePartition,
        coordinateProjector, coordinateCoherence]
    unfold EventAlgebra.IsPartitionOffDiagonal
    rw [← hcorner]
    exact EventAlgebra.partitionPinching_partitionCorner_eq_zero
      coordinatePartition (by decide) coordinateCoherence

variable {n k : ℕ}

/-- The declared block readout has exactly the pinching fibres: equality of
block readouts is equivalent to equality of every trace statistic against a
matrix in the sector-preserving commutant.  This is distinct from the
commutative partition-average readout. -/
theorem declaredBlockReadout_eq_iff_operationallyEquivalent
    (bridge : FiniteQuantumMatterBridge n k)
    (X Y : Matrix (Fin n) (Fin n) ℂ) :
    bridge.blockReadout X = bridge.blockReadout Y ↔
      EventAlgebra.PartitionOperationallyEquivalent
        bridge.partition X Y := by
  rw [bridge.blockReadout_eq_pinching]
  exact
    (EventAlgebra.partitionOperationallyEquivalent_iff_pinching_eq
      bridge.partition X Y).symm

/-- Mapped rows with different supplied sector labels have orthogonal nonzero
partition projectors. -/
theorem mapped_component_sectors_orthogonal
    (bridge : FiniteQuantumMatterBridge n k)
    {i j : Fin 10}
    (hsector : bridge.sectorOfComponent i ≠ bridge.sectorOfComponent j) :
    bridge.partition.proj (bridge.sectorOfComponent i) *
        bridge.partition.proj (bridge.sectorOfComponent j) = 0 := by
  exact bridge.partition.orthogonal _ _ hsector

/-- A corner between two differently labelled mapped sectors is invisible to
every matrix in the supplied partition commutant. -/
theorem mapped_cross_component_corner_invisible
    (bridge : FiniteQuantumMatterBridge n k)
    {i j : Fin 10}
    (hsector : bridge.sectorOfComponent i ≠ bridge.sectorOfComponent j)
    (X C : Matrix (Fin n) (Fin n) ℂ)
    (hC : C ∈ bridge.partition.commutant) :
    (EventAlgebra.partitionCorner bridge.partition
        (bridge.sectorOfComponent i) (bridge.sectorOfComponent j) X * C).trace = 0 := by
  have hoff : EventAlgebra.IsPartitionOffDiagonal bridge.partition
      (EventAlgebra.partitionCorner bridge.partition
        (bridge.sectorOfComponent i) (bridge.sectorOfComponent j) X) := by
    unfold EventAlgebra.IsPartitionOffDiagonal
    exact EventAlgebra.partitionPinching_partitionCorner_eq_zero
      bridge.partition hsector X
  exact EventAlgebra.trace_mul_eq_zero_of_partitionOffDiagonal
    bridge.partition hoff C hC

/-- After the declared component-to-sector identification, the common kernel
of the supplied mapped-sector weights is exactly the Z6 tensor-character
kernel.  No group action on the partition projectors is asserted here. -/
theorem kernel_on_mapped_component_weights
    (bridge : FiniteQuantumMatterBridge n k) :
    ∀ c : OPH.TraceBalancedKernel.C,
      (∀ i : Fin 10,
        OPH.Z6Descent.phase c
          (bridge.sectorWeight (bridge.sectorOfComponent i)) = 0) ↔
        OPH.TraceBalancedKernel.tensorCharFun c = 0 := by
  intro c
  simpa only [bridge.sectorWeight_component] using
    kernel_on_exterior_component_weights c

/-- The typed bridge's admissible nonempty chiral anomaly-free selection is
one of the two conjugate fermionic-parity rows selected by the exhaustive
scan.  The theorem does not choose between them. -/
theorem bridge_selection_is_parity_sector
    (bridge : FiniteQuantumMatterBridge n k) :
    bridge.selectionMask.val = OPH.ExteriorSelection.evenMask ∨
      bridge.selectionMask.val = OPH.ExteriorSelection.oddMask :=
  OPH.ExteriorSelection.selection_unique bridge.selectionMask
    bridge.selection_nonempty bridge.selection_chiral
    bridge.selection_anomalyFree

#print axioms componentWeight
#print axioms even_component_weights_eq_matterWeights
#print axioms kernel_on_exterior_component_weights
#print axioms fractional_singlet_mutation_collapses_component_kernel
#print axioms coordinatePartition
#print axioms coordinateBridgeControl
#print axioms coordinate_diagonal_not_partitionOffDiagonal
#print axioms coordinate_nonzero_offDiagonal_control
#print axioms declaredBlockReadout_eq_iff_operationallyEquivalent
#print axioms mapped_component_sectors_orthogonal
#print axioms mapped_cross_component_corner_invisible
#print axioms kernel_on_mapped_component_weights
#print axioms bridge_selection_is_parity_sector

end OPH.QuantumMatterIntegration
