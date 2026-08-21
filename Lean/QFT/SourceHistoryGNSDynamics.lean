import InformationProjection.SourceHistoryPacket
import QFT.ColimitSelectedStateGNS
import QFT.TowerAnchoredDiamond
import EventAlgebra.SchroedingerFrameFlow
import Mathlib.RingTheory.SimpleRing.Matrix

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Source-history state and dynamics on a selected-state GNS space

This module attaches the exact eight-history packet of
`InformationProjection.SourceHistoryPacket` to one finite quantum carrier.
The empirical path law becomes a faithful diagonal density matrix on
`M₈(ℂ)`, while the exact repair-move count becomes a positive diagonal
Hamiltonian.  Both objects inhabit the same one-regulator `ConsensusTower`.
The tower-selected state then constructs the completed-colimit GNS Hilbert
space on which the Hamiltonian, its commutator generator, and its exact
unitary flow are represented.

The diagonal entries retain their source receipts: the density is the exact
normalization of the 1754 extracted path windows, and the Hamiltonian is the
number of changes along each three-state path.  The selected energy is
exactly `197 / 1754`.  An off-diagonal matrix unit has nonzero commutator
response, so the quantum dynamics does not collapse to the diagonal source
statistics.

Claim boundary.  The packet is a representation-level extraction from one
retained hash-pinned bounded run.  Its repair-count action and empirical law
were not fixed together by a prospective source contract.  Promoting the
diagonal action to an operator on all of `M₈(ℂ)` is the finite quantum
extension made here; the source data only fixes its diagonal.  The real flow
parameter has no physical clock calibration.  No physical vacuum, Lorentzian
localization, spectrum condition, continuum or infinite-volume limit,
renormalization flow, fields, particles, scattering, detector readback, or
physical time-slice theorem is constructed.
-/

namespace OPH.QFT.SourceHistoryGNSDynamics

open Matrix
open OPH.Tower
open OPH.QFT
open OPH.Dynamics
open OPH.InformationProjection
open EventAlgebra.QuantumSurface

open scoped ComplexOrder InnerProductSpace Matrix.Norms.L2Operator
  OPH.QFT.CompletionStar

noncomputable section

/-! ## The source-counted finite carrier -/

abbrev HistoryMatrix := Matrix (Fin 8) (Fin 8) ℂ

/-- The exact empirical history law, embedded as a diagonal density matrix. -/
def sourceHistoryDensity : HistoryMatrix :=
  Matrix.diagonal fun g => (sourceTauEmpR g : ℂ)

/-- The repair-move action, embedded as a diagonal finite Hamiltonian. -/
def sourceHistoryHamiltonian : HistoryMatrix :=
  Matrix.diagonal fun g => (sourceActionR g : ℂ)

/-- An off-diagonal observable used to certify genuinely quantum response. -/
def sourceHistoryProbe : HistoryMatrix :=
  Matrix.single (0 : Fin 8) (1 : Fin 8) 1

/-- The commutator response stored in the tower generator. -/
def sourceHistoryProbeGenerator : HistoryMatrix :=
  (-Complex.I) •
    (sourceHistoryHamiltonian * sourceHistoryProbe -
      sourceHistoryProbe * sourceHistoryHamiltonian)

theorem sourceHistoryDensity_isState :
    EventAlgebra.IsState sourceHistoryDensity := by
  constructor
  · refine Matrix.PosSemidef.diagonal ?_
    intro g
    exact Complex.zero_le_real.mpr (le_of_lt (sourceTauEmpR_pos g))
  · rw [sourceHistoryDensity, Matrix.trace_diagonal,
      ← Complex.ofReal_sum, sourceTauEmpR_sum, Complex.ofReal_one]

/-- The density is faithful at the finite diagonal level: every source-history
weight is strictly positive. -/
theorem sourceHistoryDensity_diagonal_pos (g : Fin 8) :
    0 < (sourceHistoryDensity g g).re := by
  simpa [sourceHistoryDensity] using sourceTauEmpR_pos g

theorem sourceHistoryHamiltonian_posSemidef :
    sourceHistoryHamiltonian.PosSemidef := by
  refine Matrix.PosSemidef.diagonal ?_
  intro g
  exact Complex.zero_le_real.mpr (by
    simp [sourceActionR])

theorem sourceHistoryHamiltonian_isSelfAdjoint :
    IsSelfAdjoint sourceHistoryHamiltonian :=
  sourceHistoryHamiltonian_posSemidef.isHermitian.isSelfAdjoint

/-- Each Hamiltonian eigenvalue is exactly the extracted path's state-change
count. -/
theorem sourceHistoryHamiltonian_eq_changes (g : Fin 8) :
    sourceHistoryHamiltonian g g =
      (((((if sourceState0 g = sourceState1 g then 0 else 1) +
        (if sourceState1 g = sourceState2 g then 0 else 1) : ℕ) : ℝ)) : ℂ) := by
  rw [sourceHistoryHamiltonian, Matrix.diagonal_apply_eq]
  simp only [sourceActionR]
  exact congrArg (fun n : ℕ => (((n : ℝ)) : ℂ))
    (sourceAction_eq_changes g)

/-- Each density diagonal entry is the exact normalization of the retained
source window count. -/
theorem sourceHistoryDensity_eq_windowCount (g : Fin 8) :
    sourceHistoryDensity g g =
      (sourceWindowCount g : ℂ) / 1754 := by
  rw [sourceHistoryDensity, Matrix.diagonal_apply_eq]
  simp only [sourceTauEmpR]
  have hr : (sourceTauEmpQ g : ℝ) =
      (sourceWindowCount g : ℝ) / 1754 := by
    rw [sourceTauEmpQ_eq_counts]
    norm_num
  rw [hr]
  norm_num

/-- The source density and action Hamiltonian commute exactly. -/
theorem sourceHistoryDensity_commutes_hamiltonian :
    Commute sourceHistoryDensity sourceHistoryHamiltonian := by
  show sourceHistoryDensity * sourceHistoryHamiltonian =
    sourceHistoryHamiltonian * sourceHistoryDensity
  rw [sourceHistoryDensity, sourceHistoryHamiltonian,
    Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  congr 1
  funext g
  ring

/-- The selected state's Hamiltonian expectation is the exact empirical mean
repair action of the retained source packet. -/
theorem sourceHistoryDensity_meanEnergy :
    EventAlgebra.expectation sourceHistoryDensity sourceHistoryHamiltonian =
      (197 / 1754 : ℂ) := by
  change (sourceHistoryDensity * sourceHistoryHamiltonian).trace = _
  rw [sourceHistoryDensity, sourceHistoryHamiltonian,
    Matrix.diagonal_mul_diagonal, Matrix.trace_diagonal]
  simp only [← Complex.ofReal_mul]
  rw [← Complex.ofReal_sum]
  change ((meanAction sourceTauEmpR sourceActionR : ℝ) : ℂ) = _
  rw [sourceTauEmpR_meanAction]
  norm_num

/-- The Hamiltonian is not scalar: histories with zero and one repair moves
have different diagonal energies. -/
theorem sourceHistoryHamiltonian_not_scalar :
    ¬ ∃ c : ℂ, sourceHistoryHamiltonian = c • 1 := by
  rintro ⟨c, hc⟩
  have h00 := congrFun (congrFun hc (0 : Fin 8)) (0 : Fin 8)
  have h11 := congrFun (congrFun hc (1 : Fin 8)) (1 : Fin 8)
  norm_num [sourceHistoryHamiltonian, sourceActionR, sourceAction] at h00 h11
  exact one_ne_zero (h11.trans h00.symm)

theorem sourceHistoryHamiltonian_ne_zero :
    sourceHistoryHamiltonian ≠ 0 := by
  intro h
  have h11 := congrFun (congrFun h (1 : Fin 8)) (1 : Fin 8)
  norm_num [sourceHistoryHamiltonian, sourceActionR, sourceAction] at h11

/-- The diagonal action produces an exact nonzero response on an off-diagonal
matrix unit.  This is an algebraic quantum extension of the source packet,
not a source-produced coherence. -/
theorem sourceHistoryProbeGenerator_entry :
    sourceHistoryProbeGenerator (0 : Fin 8) (1 : Fin 8) = Complex.I := by
  simp [sourceHistoryProbeGenerator, sourceHistoryHamiltonian,
    sourceHistoryProbe, Matrix.mul_diagonal,
    sourceActionR, sourceAction]

theorem sourceHistoryProbeGenerator_ne_zero :
    sourceHistoryProbeGenerator ≠ 0 := by
  intro h
  have h01 := congrFun (congrFun h (0 : Fin 8)) (1 : Fin 8)
  rw [sourceHistoryProbeGenerator_entry, Matrix.zero_apply] at h01
  exact Complex.I_ne_zero h01

/-! ## The one-regulator tower and selected-state GNS representation -/

def sourceHistoryState : EventAlgebra.StateMatrix 8 :=
  ⟨sourceHistoryDensity, sourceHistoryDensity_isState⟩

def sourceHistoryGenerator : HistoryMatrix →ₗ[ℂ] HistoryMatrix :=
  (-Complex.I) •
    (LinearMap.mulLeft ℂ sourceHistoryHamiltonian -
      LinearMap.mulRight ℂ sourceHistoryHamiltonian)

def sourceHistoryTower : ConsensusTower Unit := {
  ConsensusTower.constantConsensusTower
      (oneBlockPartition 8) sourceHistoryState with
  generator := fun _ _ => sourceHistoryGenerator
  generator_natural := by intros; rfl }

def observerFamily : CoherentObserverFamily sourceHistoryTower where
  observer := fun _ => ()
  coherent := fun _ => rfl

@[simp]
theorem sourceHistoryTower_state :
    sourceHistoryTower.state () () = sourceHistoryDensity :=
  rfl

@[simp]
theorem sourceHistoryTower_generator (X : HistoryMatrix) :
    sourceHistoryTower.generator () () X =
      (-Complex.I) •
        (sourceHistoryHamiltonian * X - X * sourceHistoryHamiltonian) :=
  rfl

abbrev selectedFunctional :=
  selectedColimitFunctional observerFamily

abbrev SourceHistoryColimitGNS := ColimitGNS selectedFunctional

abbrev SourceHistoryColimitGNSOperators :=
  ColimitGNSOperators selectedFunctional

noncomputable local instance sourceHistoryCompletionPartialOrder :
    PartialOrder (ColimitCompletion sourceHistoryTower) :=
  CStarAlgebra.spectralOrder (ColimitCompletion sourceHistoryTower)

local instance sourceHistoryCompletionStarOrderedRing :
    StarOrderedRing (ColimitCompletion sourceHistoryTower) :=
  CStarAlgebra.spectralOrderedRing (ColimitCompletion sourceHistoryTower)

def cyclicUnit : SourceHistoryColimitGNS :=
  colimitGNSUnitClass selectedFunctional

def stageRepresentation :
    HistoryMatrix →⋆ₐ[ℂ] SourceHistoryColimitGNSOperators :=
  (colimitGNSRepresentation selectedFunctional).comp
    (stageToCompletion sourceHistoryTower ())

def completedHamiltonian : ColimitCompletion sourceHistoryTower :=
  stageToCompletion sourceHistoryTower () sourceHistoryHamiltonian

def representedHamiltonian : SourceHistoryColimitGNSOperators :=
  stageRepresentation sourceHistoryHamiltonian

def sourceHistoryPropagator (t : ℝ) : HistoryMatrix :=
  stonePropagator sourceHistoryHamiltonian t

def sourceHistoryHeisenberg (t : ℝ) (X : HistoryMatrix) : HistoryMatrix :=
  sourceHistoryPropagator t * X * sourceHistoryPropagator (-t)

def representedPropagator (t : ℝ) : SourceHistoryColimitGNSOperators :=
  stageRepresentation (sourceHistoryPropagator t)

def representedProbeGenerator : SourceHistoryColimitGNSOperators :=
  stageRepresentation sourceHistoryProbeGenerator

theorem selectedFunctional_normalized :
    IsNormalizedColimitFunctional selectedFunctional :=
  selectedColimitFunctional_normalized observerFamily

theorem norm_cyclicUnit : ‖cyclicUnit‖ = 1 :=
  norm_colimitGNSUnitClass selectedFunctional selectedFunctional_normalized

theorem cyclicUnit_ne_zero : cyclicUnit ≠ 0 := by
  intro h
  have hn := norm_cyclicUnit
  rw [h, norm_zero] at hn
  exact zero_ne_one hn

theorem stageRepresentation_injective :
    Function.Injective stageRepresentation := by
  letI : Nontrivial SourceHistoryColimitGNS :=
    ⟨⟨cyclicUnit, 0, cyclicUnit_ne_zero⟩⟩
  exact RingHom.injective stageRepresentation.toRingHom

theorem representedHamiltonian_ne_zero : representedHamiltonian ≠ 0 := by
  exact fun h => sourceHistoryHamiltonian_ne_zero
    (stageRepresentation_injective (by
      simpa [representedHamiltonian] using h))

theorem representedHamiltonian_isSelfAdjoint :
    IsSelfAdjoint representedHamiltonian :=
  sourceHistoryHamiltonian_isSelfAdjoint.map stageRepresentation

theorem representedProbeGenerator_ne_zero :
    representedProbeGenerator ≠ 0 := by
  exact fun h => sourceHistoryProbeGenerator_ne_zero
    (stageRepresentation_injective (by
      simpa [representedProbeGenerator] using h))

/-- The selected cyclic vector recovers the exact source mean action as its
Hamiltonian matrix coefficient. -/
theorem cyclicUnit_hamiltonian_matrixCoefficient :
    ⟪cyclicUnit, representedHamiltonian cyclicUnit⟫_ℂ =
      (197 / 1754 : ℂ) := by
  rw [cyclicUnit, representedHamiltonian, stageRepresentation]
  change ⟪colimitGNSUnitClass selectedFunctional,
    colimitGNSRepresentation selectedFunctional completedHamiltonian
      (colimitGNSUnitClass selectedFunctional)⟫_ℂ = _
  rw [colimitGNSUnitClass_expectation,
    completedHamiltonian, selectedColimitFunctional_stage]
  exact sourceHistoryDensity_meanEnergy

/-- The cyclic vector is not a zero-energy vector for this source-counted
Hamiltonian.  Accordingly it is not promoted to a physical vacuum. -/
theorem representedHamiltonian_cyclicUnit_ne_zero :
    representedHamiltonian cyclicUnit ≠ 0 := by
  intro h
  have hmc := cyclicUnit_hamiltonian_matrixCoefficient
  rw [h, inner_zero_right] at hmc
  norm_num at hmc

theorem represented_tower_generator_intertwines (X : HistoryMatrix) :
    stageRepresentation (sourceHistoryTower.generator () () X) =
      (-Complex.I) •
        (representedHamiltonian * stageRepresentation X -
          stageRepresentation X * representedHamiltonian) := by
  rw [sourceHistoryTower_generator, map_smul, map_sub, map_mul, map_mul]
  rfl

/-! ## Exact represented unitary dynamics -/

theorem sourceHistoryPropagator_zero : sourceHistoryPropagator 0 = 1 :=
  stonePropagator_zero sourceHistoryHamiltonian

theorem sourceHistoryPropagator_add (s t : ℝ) :
    sourceHistoryPropagator (s + t) =
      sourceHistoryPropagator s * sourceHistoryPropagator t :=
  stonePropagator_add sourceHistoryHamiltonian s t

theorem sourceHistoryPropagator_unitary (t : ℝ) :
    sourceHistoryPropagator t ∈ unitary HistoryMatrix :=
  stonePropagator_mem_unitary sourceHistoryHamiltonian
    sourceHistoryHamiltonian_isSelfAdjoint t

theorem sourceHistoryDensity_stationary (t : ℝ) :
    sourceHistoryHeisenberg t sourceHistoryDensity = sourceHistoryDensity := by
  have hgen : Commute
      (t • ((-Complex.I) • sourceHistoryHamiltonian))
      sourceHistoryDensity :=
    (sourceHistoryDensity_commutes_hamiltonian.symm.smul_left
      (-Complex.I)).smul_left t
  have hU : Commute (sourceHistoryPropagator t) sourceHistoryDensity := by
    unfold sourceHistoryPropagator stonePropagator
    exact hgen.exp_left
  unfold sourceHistoryHeisenberg
  calc
    sourceHistoryPropagator t * sourceHistoryDensity *
          sourceHistoryPropagator (-t) =
        sourceHistoryDensity *
          (sourceHistoryPropagator t * sourceHistoryPropagator (-t)) := by
      rw [hU.eq]
      simp only [mul_assoc]
    _ = sourceHistoryDensity := by
      rw [show sourceHistoryPropagator t * sourceHistoryPropagator (-t) = 1
        from stonePropagator_mul_neg sourceHistoryHamiltonian t, mul_one]

theorem representedPropagator_zero : representedPropagator 0 = 1 := by
  rw [representedPropagator, sourceHistoryPropagator_zero, map_one]

theorem representedPropagator_add (s t : ℝ) :
    representedPropagator (s + t) =
      representedPropagator s * representedPropagator t := by
  simp only [representedPropagator, sourceHistoryPropagator_add, map_mul]

theorem representedPropagator_unitary (t : ℝ) :
    representedPropagator t ∈ unitary SourceHistoryColimitGNSOperators :=
  Unitary.map_mem stageRepresentation (sourceHistoryPropagator_unitary t)

theorem representedHeisenberg_intertwines (t : ℝ) (X : HistoryMatrix) :
    stageRepresentation (sourceHistoryHeisenberg t X) =
      representedPropagator t * stageRepresentation X *
        representedPropagator (-t) := by
  simp only [sourceHistoryHeisenberg, representedPropagator, map_mul]

theorem representedDensity_stationary (t : ℝ) :
    stageRepresentation (sourceHistoryHeisenberg t sourceHistoryDensity) =
      stageRepresentation sourceHistoryDensity := by
  rw [sourceHistoryDensity_stationary]

/-- One bundled receipt for the source-history finite dynamical GNS rung. -/
theorem sourceHistoryFiniteGNSAttachment :
    EventAlgebra.IsState sourceHistoryDensity ∧
      (∀ g : Fin 8, 0 < (sourceHistoryDensity g g).re) ∧
      sourceHistoryHamiltonian.PosSemidef ∧
      IsSelfAdjoint representedHamiltonian ∧
      representedHamiltonian ≠ 0 ∧
      ‖cyclicUnit‖ = 1 ∧
      cyclicUnit ≠ 0 ∧
      ⟪cyclicUnit, representedHamiltonian cyclicUnit⟫_ℂ =
        (197 / 1754 : ℂ) ∧
      representedHamiltonian cyclicUnit ≠ 0 ∧
      representedProbeGenerator ≠ 0 ∧
      (∀ t : ℝ, representedPropagator t ∈
        unitary SourceHistoryColimitGNSOperators) ∧
      (∀ t : ℝ, ∀ X : HistoryMatrix,
        stageRepresentation (sourceHistoryHeisenberg t X) =
          representedPropagator t * stageRepresentation X *
            representedPropagator (-t)) ∧
      (∀ t : ℝ,
        stageRepresentation
            (sourceHistoryHeisenberg t sourceHistoryDensity) =
          stageRepresentation sourceHistoryDensity) := by
  exact ⟨sourceHistoryDensity_isState,
    sourceHistoryDensity_diagonal_pos,
    sourceHistoryHamiltonian_posSemidef,
    representedHamiltonian_isSelfAdjoint,
    representedHamiltonian_ne_zero,
    norm_cyclicUnit,
    cyclicUnit_ne_zero,
    cyclicUnit_hamiltonian_matrixCoefficient,
    representedHamiltonian_cyclicUnit_ne_zero,
    representedProbeGenerator_ne_zero,
    representedPropagator_unitary,
    representedHeisenberg_intertwines,
    representedDensity_stationary⟩

end

end OPH.QFT.SourceHistoryGNSDynamics

-- Axiom audit: all declarations must remain on the standard Mathlib basis.
#print axioms OPH.QFT.SourceHistoryGNSDynamics.sourceHistoryDensity_isState
#print axioms OPH.QFT.SourceHistoryGNSDynamics.sourceHistoryDensity_meanEnergy
#print axioms OPH.QFT.SourceHistoryGNSDynamics.sourceHistoryHamiltonian_eq_changes
#print axioms OPH.QFT.SourceHistoryGNSDynamics.sourceHistoryProbeGenerator_ne_zero
#print axioms OPH.QFT.SourceHistoryGNSDynamics.stageRepresentation_injective
#print axioms OPH.QFT.SourceHistoryGNSDynamics.cyclicUnit_hamiltonian_matrixCoefficient
#print axioms OPH.QFT.SourceHistoryGNSDynamics.representedHeisenberg_intertwines
#print axioms OPH.QFT.SourceHistoryGNSDynamics.sourceHistoryFiniteGNSAttachment
