import QFT.ColimitSelectedStateGNS
import QFT.FiniteTwoSiteIsingField
import QFT.TowerAnchoredDiamond
import Mathlib.RingTheory.SimpleRing.Matrix

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# The finite Ising dynamics on its selected-state colimit GNS space

This module closes one precise composition gap between the existing finite
two-site Ising benchmark and the selected-state GNS construction.  It
reindexes the two-qubit algebra onto the rigid `Fin 4` tower carrier, installs
the Ising ground density and commutator generator in one declared one-regulator
`ConsensusTower`, transports the genuine two-slot conditional-expectation
diamond onto that same carrier, and represents the resulting Hamiltonian and
flow on the tower-selected GNS Hilbert space.  The exact object is a finite
regional diamond plus co-located dynamics on the ambient matrix algebra, not
a dynamics-preserved regional net: the interacting first-order response of a
left observable is proved to leave the left regional algebra.

The load-bearing controls are internal to that one carrier: the cyclic unit
class is normalized and has zero energy, the represented Hamiltonian and the
represented left-spin generator are nonzero, the finite unitary group and
Heisenberg flow intertwine exactly with the GNS representation, and deleting
the interaction makes the represented left-spin generator vanish.  The
ground state is *not* preserved by the left/right conditional expectations;
that explicit incompatibility prevents this finite bridge from being read as
a state-compatible physical regional QFT.

The state, Hamiltonian, regulator, slot factorization, and real flow parameter
remain declared rather than source-selected or physically calibrated.  The
ground sector is degenerate.  No Lorentzian localization, spectrum condition,
continuum or infinite-volume limit, renormalization group, fields, particles,
scattering, detector readback, or physical time-slice theorem is claimed.
-/

namespace OPH.QFT.FiniteTwoSiteIsingColimitGNS

open Matrix
open Kronecker
open OPH.Tower
open OPH.QFT
open EventAlgebra
open EventAlgebra.ProductSplitSeparability
open EventAlgebra.Robertson
open OPH.QFT.FiniteTwoSiteIsingField

open scoped ComplexOrder InnerProductSpace Matrix.Norms.L2Operator
  OPH.QFT.CompletionStar

noncomputable section

/-! ## One common finite carrier -/

abbrev StageMatrix := Matrix (Fin 4) (Fin 4) ℂ

/-- The canonical lexicographic identification of the two spin slots with
the rigid four-dimensional tower index. -/
def carrierEquiv : IsingIndex ≃ Fin 4 :=
  finProdFinEquiv

/-- Reindexing the two-qubit matrix algebra onto the tower carrier. -/
def stageEquiv : IsingMatrix ≃⋆ₐ[ℂ] StageMatrix :=
  StarAlgEquiv.ofAlgEquiv (Matrix.reindexAlgEquiv ℂ ℂ carrierEquiv)
    (fun M => by
      simp only [Matrix.reindexAlgEquiv_apply, Matrix.star_eq_conjTranspose]
      exact reindex_conjTranspose carrierEquiv M)

@[simp]
theorem stageEquiv_apply (M : IsingMatrix) :
    stageEquiv M = Matrix.reindex carrierEquiv carrierEquiv M :=
  rfl

/-- The Ising Hamiltonian on the exact tower carrier. -/
def stageHamiltonian : StageMatrix :=
  stageEquiv isingHamiltonian

/-- The chosen Ising ground density on the exact tower carrier. -/
def stageGroundDensity : StageMatrix :=
  stageEquiv isingGroundDensity

/-- The left Pauli observable on the exact tower carrier. -/
def stageLeftPauliX : StageMatrix :=
  stageEquiv leftPauliX

/-- The right Pauli observable on the exact tower carrier. -/
def stageRightPauliX : StageMatrix :=
  stageEquiv rightPauliX

/-- A right-local spin observable used to expose the state/expectation
compatibility boundary. -/
def rawRightPauliZ : IsingMatrix :=
  slotRight Spin pauliZ

/-- The same right-local spin observable on the tower carrier. -/
def stageRightPauliZ : StageMatrix :=
  stageEquiv rawRightPauliZ

/-- The nonzero Ising commutator response on the exact tower carrier. -/
def stageLeftPauliXGenerator : StageMatrix :=
  stageEquiv isingLeftPauliXGenerator

/-- The interaction-deleted Hamiltonian on the exact tower carrier. -/
def stageUncoupledHamiltonian : StageMatrix :=
  stageEquiv uncoupledHamiltonian

/-- The reindexed chosen density is a certified finite state. -/
def stageGroundState : StateMatrix 4 :=
  ⟨stageGroundDensity, by
    constructor
    · exact posSemidef_reindex_equiv carrierEquiv
        isingGroundDensity_posSemidef
    · rw [stageGroundDensity, stageEquiv_apply,
        trace_reindex_equiv carrierEquiv, isingGroundDensity_trace]⟩

/-- The Ising commutator derivation stored in the tower generator field. -/
def stageGenerator : StageMatrix →ₗ[ℂ] StageMatrix :=
  (-Complex.I) •
    (LinearMap.mulLeft ℂ stageHamiltonian -
      LinearMap.mulRight ℂ stageHamiltonian)

/-- A one-regulator tower whose selected state and generator are the finite
Ising state and Hamiltonian commutator on the same `M₄(ℂ)` carrier. -/
def isingTower : ConsensusTower Unit := {
  ConsensusTower.constantConsensusTower
      (oneBlockPartition 4) stageGroundState with
  generator := fun _ _ => stageGenerator
  generator_natural := by intros; rfl }

/-- The unique coherent observer choice in the one-regulator tower. -/
def observerFamily : CoherentObserverFamily isingTower where
  observer := fun _ => ()
  coherent := fun _ => rfl

@[simp]
theorem isingTower_state :
    isingTower.state () () = stageGroundDensity :=
  rfl

@[simp]
theorem isingTower_generator (X : StageMatrix) :
    isingTower.generator () () X =
      (-Complex.I) • (stageHamiltonian * X - X * stageHamiltonian) := by
  rfl

/-! ## The two-slot conditional-expectation net on the same carrier -/

/-- Regional algebras for the two Ising spin slots before reindexing. -/
def rawSlotAlgebra : TwoSlotRegion →
    StarSubalgebra ℂ IsingMatrix
  | .bot => ⊥
  | .left => (slotLeft Spin (α := Spin)).range
  | .right => (slotRight Spin (β := Spin)).range
  | .top => ⊤

/-- Conditional expectations for the two Ising spin slots before
reindexing. -/
def rawSlotExpectation : TwoSlotRegion → IsingMatrix →ₗ[ℂ] IsingMatrix
  | .bot => scalarExpectation IsingIndex
  | .left => leftSlotExpectation
  | .right => rightSlotExpectation
  | .top => LinearMap.id

private theorem raw_left_expectation_mem (M : IsingMatrix) :
    leftSlotExpectation M ∈
      (slotLeft Spin (α := Spin)).range := by
  have h : leftSlotExpectation (α := Spin) (β := Spin) M =
      (Fintype.card Spin : ℂ)⁻¹ •
        (ptraceSnd M ⊗ₖ (1 : Matrix Spin Spin ℂ)) := rfl
  rw [h]
  exact SMulMemClass.smul_mem _ ⟨ptraceSnd M, rfl⟩

private theorem raw_right_expectation_mem (M : IsingMatrix) :
    rightSlotExpectation M ∈
      (slotRight Spin (β := Spin)).range := by
  have h : rightSlotExpectation (α := Spin) (β := Spin) M =
      (Fintype.card Spin : ℂ)⁻¹ •
        ((1 : Matrix Spin Spin ℂ) ⊗ₖ OPH.Locality.ptraceFst M) := rfl
  rw [h]
  exact SMulMemClass.smul_mem _ ⟨OPH.Locality.ptraceFst M, rfl⟩

/-- The genuine scalar/left/right/top conditional-expectation diamond for
the two Ising spin slots.  Its left and right algebras jointly generate the
whole two-qubit algebra. -/
def rawSlotNet : CPRegionalNet IsingIndex where
  Region := TwoSlotRegion
  regionFintype := inferInstance
  regionNonempty := ⟨TwoSlotRegion.bot⟩
  regionLE U V := TwoSlotRegion.le U V = true
  regionLE_refl := TwoSlotRegion.le_refl'
  regionLE_trans {U V W} := TwoSlotRegion.le_trans' U V W
  regionLE_antisymm {U V} := TwoSlotRegion.le_antisymm' U V
  overlap := TwoSlotRegion.meet
  overlap_le_left := TwoSlotRegion.meet_le_left'
  overlap_le_right := TwoSlotRegion.meet_le_right'
  le_overlap {W U V} := TwoSlotRegion.le_meet' W U V
  disjoint U V := TwoSlotRegion.disj U V = true
  disjoint_symm {U V} := TwoSlotRegion.disj_symm' U V
  disjoint_irrefl := TwoSlotRegion.disj_irrefl'
  localAlgebra := rawSlotAlgebra
  isotony {U V} hUV := by
    cases U <;> cases V <;> (try (simp [TwoSlotRegion.le] at hUV))
    · exact le_refl _
    · exact bot_le
    · exact bot_le
    · exact bot_le
    · exact le_refl _
    · exact le_top
    · exact le_refl _
    · exact le_top
    · exact le_refl _
  locality {U V} hUV := by
    cases U <;> cases V <;> (try (simp [TwoSlotRegion.disj] at hUV))
    · exact fun X hX Y hY => slot_commute hX hY
    · exact fun X hX Y hY => (slot_commute hY hX).symm
  expect := rawSlotExpectation
  expect_mem U X := by
    cases U
    · exact scalarExpectation_mem_bot X
    · exact raw_left_expectation_mem X
    · exact raw_right_expectation_mem X
    · trivial
  expect_fixes U X hX := by
    cases U
    · exact scalarExpectation_fixes X hX
    · obtain ⟨A, rfl⟩ := hX
      exact leftSlotExpectation_fixes_left A
    · obtain ⟨B, rfl⟩ := hX
      exact rightSlotExpectation_fixes_right B
    · rfl
  expect_posSemidef U {X} hX := by
    cases U
    · exact scalarExpectation_posSemidef hX
    · exact leftSlotExpectation_posSemidef hX
    · exact rightSlotExpectation_posSemidef hX
    · exact hX
  expect_trace U X := by
    cases U
    · exact scalarExpectation_trace X
    · exact leftSlotExpectation_trace X
    · exact rightSlotExpectation_trace X
    · rfl
  expect_tower {U V} hUV X := by
    cases U <;> cases V <;> (try (simp [TwoSlotRegion.le] at hUV))
    · exact scalarExpectation_of_trace_eq (scalarExpectation_trace X)
    · exact scalarExpectation_of_trace_eq (leftSlotExpectation_trace X)
    · exact scalarExpectation_of_trace_eq (rightSlotExpectation_trace X)
    · rfl
    · exact leftSlotExpectation_idem X
    · rfl
    · exact rightSlotExpectation_idem X
    · rfl
    · rfl
  generating := {.left, .right}
  coverage := by
    have hset : (⋃ U ∈ ({TwoSlotRegion.left, TwoSlotRegion.right} :
        Finset TwoSlotRegion), (rawSlotAlgebra U : Set IsingMatrix)) =
        (((slotLeft Spin (α := Spin)).range : Set IsingMatrix) ∪
          ((slotRight Spin (β := Spin)).range : Set IsingMatrix)) := by
      ext X
      simp only [Finset.mem_insert, Finset.mem_singleton, Set.mem_iUnion,
        Set.mem_union, exists_prop]
      constructor
      · rintro ⟨U, (rfl | rfl), hX⟩
        · exact Or.inl hX
        · exact Or.inr hX
      · rintro (hX | hX)
        · exact ⟨TwoSlotRegion.left, Or.inl rfl, hX⟩
        · exact ⟨TwoSlotRegion.right, Or.inr rfl, hX⟩
    rw [hset]
    exact slot_ranges_generate_top

/-- The two-slot diamond transported onto the tower's exact `Fin 4` carrier. -/
def slotNet : CPRegionalNet (Fin 4) :=
  rawSlotNet.transport stageEquiv
    (fun hX => posSemidef_reindex_equiv carrierEquiv hX)
    (fun hY => posSemidef_reindex_equiv carrierEquiv.symm hY)
    (fun X => trace_reindex_equiv carrierEquiv X)

theorem stageLeftPauliX_mem_left :
    stageLeftPauliX ∈ slotNet.localAlgebra .left := by
  exact ⟨leftPauliX, ⟨pauliX, rfl⟩, rfl⟩

theorem stageRightPauliX_mem_right :
    stageRightPauliX ∈ slotNet.localAlgebra .right := by
  exact ⟨rightPauliX, ⟨pauliX, rfl⟩, rfl⟩

theorem stageHamiltonian_mem_top :
    stageHamiltonian ∈ slotNet.localAlgebra .top := by
  exact ⟨isingHamiltonian, trivial, rfl⟩

/-! ## Exact regional controls before representation -/

theorem raw_left_expectation_rightPauliZ_eq_zero :
    rawSlotNet.expect .left rawRightPauliZ = 0 := by
  change leftSlotExpectation
      ((1 : Matrix Spin Spin ℂ) ⊗ₖ pauliZ) = 0
  rw [leftSlotExpectation_scalarises_right]
  norm_num [pauliZ, Matrix.trace, Fin.sum_univ_two]

theorem slotNet_left_expectation_rightPauliZ_eq_zero :
    slotNet.expect .left stageRightPauliZ = 0 := by
  change stageEquiv
      (rawSlotNet.expect .left (stageEquiv.symm stageRightPauliZ)) = 0
  rw [stageRightPauliZ, StarAlgEquiv.symm_apply_apply,
    raw_left_expectation_rightPauliZ_eq_zero, map_zero]

theorem ground_expectation_rawRightPauliZ :
    (isingGroundDensity * rawRightPauliZ).trace = 1 := by
  norm_num [isingGroundDensity, ket00, rawRightPauliZ,
    slotRight, pauliZ, Matrix.mul_apply, Matrix.kroneckerMap_apply,
    Fintype.sum_prod_type, Fin.sum_univ_two, Matrix.single_apply,
    Matrix.trace]

theorem ground_expectation_stageRightPauliZ :
    EventAlgebra.expectation stageGroundDensity stageRightPauliZ = 1 := by
  change (stageGroundDensity * stageRightPauliZ).trace = 1
  rw [stageGroundDensity, stageRightPauliZ, ← map_mul, stageEquiv_apply,
    trace_reindex_equiv carrierEquiv]
  exact ground_expectation_rawRightPauliZ

/-- Boundary receipt: the left conditional expectation is trace preserving
but not preserving for the selected pure ground state.  Hence the regional
diamond and selected GNS sector are composed on one carrier, but they do not
form a state-compatible conditional-expectation net. -/
theorem selectedGroundState_not_leftExpectationInvariant :
    EventAlgebra.expectation stageGroundDensity
        (slotNet.expect .left stageRightPauliZ) ≠
      EventAlgebra.expectation stageGroundDensity stageRightPauliZ := by
  rw [slotNet_left_expectation_rightPauliZ_eq_zero,
    ground_expectation_stageRightPauliZ]
  change EventAlgebra.stateExpectationLinearMap stageGroundDensity 0 ≠ 1
  rw [map_zero]
  exact zero_ne_one

/-- The interacting left-spin response fails to commute with a right-spin
observable; the coupling spreads the response beyond the left factor at
first order. -/
theorem isingLeftGenerator_rightPauliX_commutator_entry :
    (isingLeftPauliXGenerator * rightPauliX -
        rightPauliX * isingLeftPauliXGenerator)
      ((0 : Spin), (0 : Spin)) ((1 : Spin), (1 : Spin)) =
        2 * Complex.I := by
  norm_num [isingLeftPauliXGenerator, isingHamiltonian, isingZZ,
    leftPauliX, rightPauliX, slotLeft, slotRight, pauliZ, pauliX,
    Matrix.mul_apply, Matrix.kroneckerMap_apply,
    Fintype.sum_prod_type, Fin.sum_univ_two]
  ring

theorem isingLeftPauliXGenerator_not_mem_left :
    isingLeftPauliXGenerator ∉
      (slotLeft Spin (α := Spin)).range := by
  intro hleft
  have hcomm : Commute isingLeftPauliXGenerator rightPauliX :=
    slot_commute hleft ⟨pauliX, rfl⟩
  have hzero :
      (isingLeftPauliXGenerator * rightPauliX -
          rightPauliX * isingLeftPauliXGenerator)
        ((0 : Spin), (0 : Spin)) ((1 : Spin), (1 : Spin)) = 0 := by
    rw [hcomm.eq, sub_self]
    rfl
  rw [isingLeftGenerator_rightPauliX_commutator_entry] at hzero
  exact (mul_ne_zero (by norm_num : (2 : ℂ) ≠ 0) Complex.I_ne_zero)
    hzero

theorem stageLeftPauliXGenerator_not_mem_left :
    stageLeftPauliXGenerator ∉ slotNet.localAlgebra .left := by
  intro hmem
  obtain ⟨X, hX, hEq⟩ := hmem
  apply isingLeftPauliXGenerator_not_mem_left
  have hsource : X = isingLeftPauliXGenerator := by
    apply stageEquiv.injective
    simpa only [stageLeftPauliXGenerator] using hEq
  simpa only [rawSlotAlgebra, hsource] using hX

/-! ## The selected-state GNS representation -/

/-- The positive normalized functional selected by the Ising tower and its
unique coherent observer family. -/
abbrev selectedFunctional :=
  selectedColimitFunctional observerFamily

/-- The actual Hilbert completion produced by the selected Ising state. -/
abbrev IsingColimitGNS := ColimitGNS selectedFunctional

/-- Bounded operators on the selected Ising GNS Hilbert space. -/
abbrev IsingColimitGNSOperators := ColimitGNSOperators selectedFunctional

/- The GNS adapter uses the canonical spectral order locally. -/
noncomputable local instance isingCompletionPartialOrder :
    PartialOrder (ColimitCompletion isingTower) :=
  CStarAlgebra.spectralOrder (ColimitCompletion isingTower)

local instance isingCompletionStarOrderedRing :
    StarOrderedRing (ColimitCompletion isingTower) :=
  CStarAlgebra.spectralOrderedRing (ColimitCompletion isingTower)

/-- The cyclic unit class.  It is not called a physical vacuum. -/
def cyclicUnit : IsingColimitGNS :=
  colimitGNSUnitClass selectedFunctional

/-- The exact composite representation from the common `M₄(ℂ)` stage
through the colimit completion to bounded operators on the selected GNS
Hilbert space. -/
def stageRepresentation : StageMatrix →⋆ₐ[ℂ] IsingColimitGNSOperators :=
  (colimitGNSRepresentation selectedFunctional).comp
    (stageToCompletion isingTower ())

/-- The represented image of a region in the two-slot diamond. -/
def representedLocalAlgebra (U : slotNet.Region) :
    StarSubalgebra ℂ IsingColimitGNSOperators :=
  (slotNet.localAlgebra U).map stageRepresentation

/-- The completed-colimit Ising Hamiltonian. -/
def completedHamiltonian : ColimitCompletion isingTower :=
  stageToCompletion isingTower () stageHamiltonian

/-- The Hamiltonian acting on the selected-state GNS Hilbert space. -/
def representedHamiltonian : IsingColimitGNSOperators :=
  stageRepresentation stageHamiltonian

/-- The represented finite propagator.  Its parameter is not asserted to be
physical time. -/
def representedPropagator (t : ℝ) : IsingColimitGNSOperators :=
  stageRepresentation (stageEquiv (isingPropagator t))

/-- The represented nonzero response of the left Pauli observable. -/
def representedLeftPauliXGenerator : IsingColimitGNSOperators :=
  stageRepresentation stageLeftPauliXGenerator

/-- The represented interaction-deleted Hamiltonian. -/
def representedUncoupledHamiltonian : IsingColimitGNSOperators :=
  stageRepresentation stageUncoupledHamiltonian

theorem stageHamiltonian_isSelfAdjoint : IsSelfAdjoint stageHamiltonian :=
  isingHamiltonian_isSelfAdjoint.map stageEquiv

theorem stageHamiltonian_mul_self :
    stageHamiltonian * stageHamiltonian = stageHamiltonian := by
  rw [stageHamiltonian, ← map_mul, isingHamiltonian_mul_self]

theorem stageGroundDensity_mul_hamiltonian :
    stageGroundDensity * stageHamiltonian = 0 := by
  rw [stageGroundDensity, stageHamiltonian, ← map_mul,
    isingGroundDensity_mul_hamiltonian, map_zero]

theorem stageHamiltonian_mul_groundDensity :
    stageHamiltonian * stageGroundDensity = 0 := by
  rw [stageGroundDensity, stageHamiltonian, ← map_mul,
    isingHamiltonian_mul_groundDensity, map_zero]

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

/-- The selected finite density assigns zero energy to the attached
Hamiltonian. -/
theorem selectedFunctional_completedHamiltonian :
    selectedFunctional completedHamiltonian = 0 := by
  rw [completedHamiltonian, selectedColimitFunctional_stage]
  change (stageGroundDensity * stageHamiltonian).trace = 0
  rw [stageGroundDensity_mul_hamiltonian, Matrix.trace_zero]

/-- The represented Hamiltonian belongs to the represented top region of the
same two-slot net. -/
theorem representedHamiltonian_mem_top :
    representedHamiltonian ∈ representedLocalAlgebra .top :=
  StarSubalgebra.mem_map.mpr ⟨stageHamiltonian,
    stageHamiltonian_mem_top, rfl⟩

/-- The represented left and right Pauli observables retain their declared
regional memberships. -/
theorem represented_pauli_memberships :
    stageRepresentation stageLeftPauliX ∈ representedLocalAlgebra .left ∧
      stageRepresentation stageRightPauliX ∈ representedLocalAlgebra .right :=
  ⟨StarSubalgebra.mem_map.mpr ⟨stageLeftPauliX,
      stageLeftPauliX_mem_left, rfl⟩,
    StarSubalgebra.mem_map.mpr ⟨stageRightPauliX,
      stageRightPauliX_mem_right, rfl⟩⟩

/-- The declared disjoint left and right regional images commute in the GNS
representation. -/
theorem represented_left_right_locality
    {A B : IsingColimitGNSOperators}
    (hA : A ∈ representedLocalAlgebra .left)
    (hB : B ∈ representedLocalAlgebra .right) : Commute A B := by
  obtain ⟨A₀, hA₀, rfl⟩ := StarSubalgebra.mem_map.mp hA
  obtain ⟨B₀, hB₀, rfl⟩ := StarSubalgebra.mem_map.mp hB
  exact (slotNet.locality (U := .left) (V := .right) (by rfl)
    A₀ hA₀ B₀ hB₀).map stageRepresentation

/-- The attached Hamiltonian annihilates the cyclic unit class.  The proof
uses its star-square, not merely the vanishing one-point expectation. -/
theorem representedHamiltonian_cyclicUnit_eq_zero :
    representedHamiltonian cyclicUnit = 0 := by
  have hstarSquare :
      star completedHamiltonian * completedHamiltonian =
        completedHamiltonian := by
    rw [completedHamiltonian, ← map_star, ← map_mul]
    have hstar : star stageHamiltonian = stageHamiltonian :=
      stageHamiltonian_isSelfAdjoint
    exact congrArg (stageToCompletion isingTower ()) (by
      calc
        star stageHamiltonian * stageHamiltonian =
            stageHamiltonian * stageHamiltonian :=
          congrArg (fun Y => Y * stageHamiltonian) hstar
        _ = stageHamiltonian := stageHamiltonian_mul_self)
  have hn : ‖selectedFunctional.toPreGNS completedHamiltonian‖ = 0 := by
    rw [PositiveLinearMap.preGNS_norm_def,
      PositiveLinearMap.ofPreGNS_toPreGNS, hstarSquare,
      selectedFunctional_completedHamiltonian]
    simp
  apply norm_eq_zero.mp
  rw [representedHamiltonian, stageRepresentation, cyclicUnit]
  change ‖(colimitGNSRepresentation selectedFunctional completedHamiltonian)
      (colimitGNSUnitClass selectedFunctional)‖ = 0
  rw [colimitGNSRepresentation_apply_unitClass,
    UniformSpace.Completion.norm_coe]
  exact hn

theorem representedHamiltonian_isSelfAdjoint :
    IsSelfAdjoint representedHamiltonian :=
  stageHamiltonian_isSelfAdjoint.map stageRepresentation

theorem representedHamiltonian_mul_self :
    representedHamiltonian * representedHamiltonian =
      representedHamiltonian := by
  rw [representedHamiltonian, ← map_mul, stageHamiltonian_mul_self]

theorem stageRepresentation_injective :
    Function.Injective stageRepresentation := by
  letI : Nontrivial IsingColimitGNS :=
    ⟨⟨cyclicUnit, 0, cyclicUnit_ne_zero⟩⟩
  exact RingHom.injective stageRepresentation.toRingHom

theorem stageHamiltonian_ne_zero : stageHamiltonian ≠ 0 := by
  intro hH
  have hH₀ : isingHamiltonian = 0 := by
    have hEq : stageEquiv isingHamiltonian = stageEquiv 0 := by
      simpa only [stageHamiltonian, map_zero] using hH
    exact stageEquiv.injective hEq
  apply isingLeftPauliXGenerator_ne_zero
  simp [isingLeftPauliXGenerator, hH₀]

/-- The Hamiltonian survives the selected-state GNS quotient; its zero action
on the cyclic ground vector is therefore not a collapsed representation. -/
theorem representedHamiltonian_ne_zero : representedHamiltonian ≠ 0 := by
  exact fun h => stageHamiltonian_ne_zero
    (stageRepresentation_injective (by
      simpa [representedHamiltonian] using h))

theorem stageLeftPauliXGenerator_ne_zero :
    stageLeftPauliXGenerator ≠ 0 := by
  intro hG
  apply isingLeftPauliXGenerator_ne_zero
  have hEq : stageEquiv isingLeftPauliXGenerator = stageEquiv 0 := by
    simpa only [stageLeftPauliXGenerator, map_zero] using hG
  exact stageEquiv.injective hEq

/-- The interaction response also survives representation. -/
theorem representedLeftPauliXGenerator_ne_zero :
    representedLeftPauliXGenerator ≠ 0 := by
  exact fun h => stageLeftPauliXGenerator_ne_zero
    (stageRepresentation_injective (by
      simpa [representedLeftPauliXGenerator] using h))

/-- The represented interaction response leaves the represented left factor.
This is the finite first-order spreading receipt, not a relativistic
propagation statement. -/
theorem representedLeftPauliXGenerator_not_mem_left :
    representedLeftPauliXGenerator ∉ representedLocalAlgebra .left := by
  intro hmem
  obtain ⟨X, hX, hEq⟩ := StarSubalgebra.mem_map.mp hmem
  apply stageLeftPauliXGenerator_not_mem_left
  have hsource : X = stageLeftPauliXGenerator := by
    apply stageRepresentation_injective
    simpa only [representedLeftPauliXGenerator] using hEq
  simpa only [hsource] using hX

theorem isingTower_generator_leftPauliX :
    isingTower.generator () () stageLeftPauliX =
      stageLeftPauliXGenerator := by
  rw [isingTower_generator, stageHamiltonian, stageLeftPauliX,
    stageLeftPauliXGenerator, isingLeftPauliXGenerator,
    ← map_mul, ← map_mul, ← map_sub, ← map_smul]

/-- The generator stored in the tower is exactly the represented Hamiltonian
commutator after applying the selected GNS representation. -/
theorem represented_tower_generator_intertwines (X : StageMatrix) :
    stageRepresentation (isingTower.generator () () X) =
      (-Complex.I) •
        (representedHamiltonian * stageRepresentation X -
          stageRepresentation X * representedHamiltonian) := by
  rw [isingTower_generator, map_smul, map_sub, map_mul, map_mul]
  rfl

/-! ## The represented unitary group and Heisenberg dynamics -/

theorem representedPropagator_zero : representedPropagator 0 = 1 := by
  rw [representedPropagator, isingPropagator_zero, map_one, map_one]

theorem representedPropagator_add (s t : ℝ) :
    representedPropagator (s + t) =
      representedPropagator s * representedPropagator t := by
  simp only [representedPropagator, isingPropagator_add, map_mul]

theorem representedPropagator_unitary (t : ℝ) :
    representedPropagator t ∈ unitary IsingColimitGNSOperators := by
  exact Unitary.map_mem stageRepresentation
    (Unitary.map_mem stageEquiv (isingPropagator_unitary t))

/-- The original finite Heisenberg flow is represented by exact conjugation
on the selected-state GNS Hilbert space. -/
theorem representedHeisenberg_intertwines (t : ℝ) (X : IsingMatrix) :
    stageRepresentation (stageEquiv (isingHeisenberg t X)) =
      representedPropagator t * stageRepresentation (stageEquiv X) *
        representedPropagator (-t) := by
  simp only [isingHeisenberg, representedPropagator, map_mul]

theorem representedGroundDensity_stationary (t : ℝ) :
    stageRepresentation
        (stageEquiv (isingHeisenberg t isingGroundDensity)) =
      stageRepresentation stageGroundDensity := by
  rw [isingGroundDensity_stationary, stageGroundDensity]

/-- The interaction-deletion control has zero represented response for the
same left Pauli observable. -/
theorem represented_uncoupled_left_generator_eq_zero :
    (-Complex.I) •
      (representedUncoupledHamiltonian *
          stageRepresentation stageLeftPauliX -
        stageRepresentation stageLeftPauliX *
          representedUncoupledHamiltonian) = 0 := by
  have h := congrArg
    (fun Y : IsingMatrix => stageRepresentation (stageEquiv Y))
    uncoupledLeftPauliXGenerator_eq_zero
  simpa only [representedUncoupledHamiltonian,
    stageUncoupledHamiltonian, stageLeftPauliX, map_smul, map_sub, map_mul,
    map_zero] using h

theorem represented_uncoupled_flow_static (t : ℝ) (X : IsingMatrix) :
    stageRepresentation (stageEquiv (uncoupledHeisenberg t X)) =
      stageRepresentation (stageEquiv X) := by
  rw [uncoupledHeisenberg_eq]

/-! ## Bundled composition receipt -/

/-- One theorem collecting the exact finite QFT-composition rung: a
normalized nonzero GNS sector, a nonzero self-adjoint Hamiltonian with a
zero-energy cyclic vector, a nonzero interaction response that leaves the
left factor, exact unitary group/Heisenberg transport, a static deletion
control, and the explicit failure of selected-state compatibility for the
regional conditional expectation. -/
theorem finiteSelectedGNSHamiltonianAttachment :
    ‖cyclicUnit‖ = 1 ∧
      cyclicUnit ≠ 0 ∧
      representedHamiltonian ≠ 0 ∧
      IsSelfAdjoint representedHamiltonian ∧
      representedHamiltonian cyclicUnit = 0 ∧
      representedHamiltonian ∈ representedLocalAlgebra .top ∧
      representedLeftPauliXGenerator ≠ 0 ∧
      representedLeftPauliXGenerator ∉ representedLocalAlgebra .left ∧
      (∀ t : ℝ, representedPropagator t ∈
        unitary IsingColimitGNSOperators) ∧
      (∀ t : ℝ, ∀ X : IsingMatrix,
        stageRepresentation (stageEquiv (isingHeisenberg t X)) =
          representedPropagator t * stageRepresentation (stageEquiv X) *
            representedPropagator (-t)) ∧
      ((-Complex.I) •
        (representedUncoupledHamiltonian *
            stageRepresentation stageLeftPauliX -
          stageRepresentation stageLeftPauliX *
            representedUncoupledHamiltonian) = 0) ∧
      EventAlgebra.expectation stageGroundDensity
          (slotNet.expect .left stageRightPauliZ) ≠
        EventAlgebra.expectation stageGroundDensity stageRightPauliZ := by
  exact ⟨norm_cyclicUnit, cyclicUnit_ne_zero,
    representedHamiltonian_ne_zero, representedHamiltonian_isSelfAdjoint,
    representedHamiltonian_cyclicUnit_eq_zero,
    representedHamiltonian_mem_top,
    representedLeftPauliXGenerator_ne_zero,
    representedLeftPauliXGenerator_not_mem_left,
    representedPropagator_unitary,
    representedHeisenberg_intertwines,
    represented_uncoupled_left_generator_eq_zero,
    selectedGroundState_not_leftExpectationInvariant⟩

end

end OPH.QFT.FiniteTwoSiteIsingColimitGNS

-- Axiom audit: all declarations must remain on the standard Mathlib basis.
#print axioms OPH.QFT.FiniteTwoSiteIsingColimitGNS.stageGroundState
#print axioms OPH.QFT.FiniteTwoSiteIsingColimitGNS.isingTower_generator
#print axioms OPH.QFT.FiniteTwoSiteIsingColimitGNS.rawSlotNet
#print axioms OPH.QFT.FiniteTwoSiteIsingColimitGNS.slotNet
#print axioms OPH.QFT.FiniteTwoSiteIsingColimitGNS.stageRepresentation_injective
#print axioms OPH.QFT.FiniteTwoSiteIsingColimitGNS.representedHamiltonian_cyclicUnit_eq_zero
#print axioms OPH.QFT.FiniteTwoSiteIsingColimitGNS.representedLeftPauliXGenerator_not_mem_left
#print axioms OPH.QFT.FiniteTwoSiteIsingColimitGNS.representedHeisenberg_intertwines
#print axioms OPH.QFT.FiniteTwoSiteIsingColimitGNS.selectedGroundState_not_leftExpectationInvariant
#print axioms OPH.QFT.FiniteTwoSiteIsingColimitGNS.finiteSelectedGNSHamiltonianAttachment
