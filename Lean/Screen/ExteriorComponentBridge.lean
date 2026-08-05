import Mathlib.LinearAlgebra.ExteriorAlgebra.Basis
import ExteriorSelection

/-!
# Five-mode exterior basis and the finite matter component table

This module connects the component table in `Screen.ExteriorSelection` to an
actual Mathlib exterior algebra.  The one-particle carrier is the supplied
five-dimensional split coordinate space.  Mathlib gives its exterior algebra
a basis indexed by all subsets of the five coordinate modes.  Counting color
and weak modes partitions those 32 basis labels into the twelve bidegrees
`(c,w)` with `0 <= c <= 3` and `0 <= w <= 2`.  Removing the vacuum `(0,0)`
and top `(3,2)` leaves exactly the ten rows of the existing selection table.

The same module exposes the exact square-zero and anticommutation relations
of exterior creation.  These are algebraic exclusion statements inside the
declared exterior module.  They do not select this module as physical matter
and do not prove spin--statistics, locality, a particle spectrum, or a
laboratory charge assignment.
-/

namespace OPH.ExteriorComponentBridge

open scoped BigOperators

abbrev Mode := Fin 5
abbrev FiveModeCarrier := Mode → ℂ
abbrev FiveModeExterior := ExteriorAlgebra ℂ FiveModeCarrier
abbrev Bidegree := Fin 4 × Fin 3

/-- The coordinate basis of the supplied five-mode carrier. -/
noncomputable def fiveModeBasis : Module.Basis Mode ℂ FiveModeCarrier :=
  Pi.basisFun ℂ Mode

/-- Mathlib's exterior basis, indexed by the 32 subsets of five modes. -/
noncomputable def fiveModeExteriorBasis :
    Module.Basis (Finset Mode) ℂ FiveModeExterior :=
  fiveModeBasis.ExteriorAlgebra

/-- Color degree of an exterior basis label: modes `0,1,2`. -/
def colorDegree (s : Finset Mode) : ℕ :=
  (s.filter fun i ↦ i.val < 3).card

/-- Weak degree of an exterior basis label: modes `3,4`. -/
def weakDegree (s : Finset Mode) : ℕ :=
  (s.filter fun i ↦ 3 ≤ i.val).card

/-- Bidegree as a pair of natural numbers, used directly on basis labels. -/
def basisBidegree (s : Finset Mode) : ℕ × ℕ :=
  (colorDegree s, weakDegree s)

/-- The ten non-vacuum, non-top bidegrees in exactly the row order of
`Screen.ExteriorSelection`. -/
def componentDegree : Fin 10 → Bidegree := ![
  (1, 0), (0, 1), (2, 0), (1, 1), (0, 2),
  (1, 2), (2, 1), (3, 0), (2, 2), (3, 1)
]

def vacuumDegree : Bidegree := (0, 0)
def topDegree : Bidegree := (3, 2)

/-- Number of exterior-basis labels in a declared natural bidegree. -/
def bidegreeCount (c w : ℕ) : ℕ :=
  ((Finset.univ : Finset (Finset Mode)).filter
    fun s ↦ basisBidegree s = (c, w)).card

/-- The actual exterior basis has 32 labels. -/
theorem exterior_basis_label_count : Fintype.card (Finset Mode) = 32 := by
  norm_num [Fintype.card_finset]

/-- Every five-mode exterior basis label has one of the twelve possible
color/weak bidegrees. -/
theorem basis_bidegree_exhaustive :
    ∀ s : Finset Mode,
      basisBidegree s = (0, 0) ∨ basisBidegree s = (1, 0) ∨
      basisBidegree s = (2, 0) ∨ basisBidegree s = (3, 0) ∨
      basisBidegree s = (0, 1) ∨ basisBidegree s = (1, 1) ∨
      basisBidegree s = (2, 1) ∨ basisBidegree s = (3, 1) ∨
      basisBidegree s = (0, 2) ∨ basisBidegree s = (1, 2) ∨
      basisBidegree s = (2, 2) ∨ basisBidegree s = (3, 2) := by
  decide +kernel

/-- The twelve bidegree multiplicities are exactly the binomial dimensions
of `Lambda^c C^3 tensor Lambda^w C^2`. -/
theorem bidegree_count_table :
    bidegreeCount 0 0 = 1 ∧ bidegreeCount 1 0 = 3 ∧
    bidegreeCount 2 0 = 3 ∧ bidegreeCount 3 0 = 1 ∧
    bidegreeCount 0 1 = 2 ∧ bidegreeCount 1 1 = 6 ∧
    bidegreeCount 2 1 = 6 ∧ bidegreeCount 3 1 = 2 ∧
    bidegreeCount 0 2 = 1 ∧ bidegreeCount 1 2 = 3 ∧
    bidegreeCount 2 2 = 3 ∧ bidegreeCount 3 2 = 1 := by
  decide +kernel

/-- The ten component-degree rows are distinct and cover every bidegree
except the exterior vacuum and top line. -/
theorem componentDegree_exact_nontrivial_menu :
    Function.Injective componentDegree ∧
      ∀ d : Bidegree,
        d = vacuumDegree ∨ d = topDegree ∨
          ∃ i : Fin 10, componentDegree i = d := by
  decide +kernel

/-- The basis-label multiplicity of every component row equals its frozen
color-dimension times weak-dimension entry. -/
theorem component_dimension_binding :
    ∀ i : Fin 10,
      bidegreeCount (componentDegree i).1.val (componentDegree i).2.val =
        Int.toNat
          (OPH.ExteriorSelection.colorDim i *
            OPH.ExteriorSelection.weakDim i) := by
  decide +kernel

/-- The frozen integer charge is the additive exterior weight `-2c+3w`. -/
theorem component_charge_binding :
    ∀ i : Fin 10,
      OPH.ExteriorSelection.charge i =
        -2 * (componentDegree i).1.val +
          3 * (componentDegree i).2.val := by
  decide +kernel

/-- The frozen parity bit is the parity of the total exterior degree. -/
theorem component_parity_binding :
    ∀ i : Fin 10,
      OPH.ExteriorSelection.odd i =
        decide (((componentDegree i).1.val +
          (componentDegree i).2.val) % 2 = 1) := by
  decide +kernel

/-- Frozen charge conjugation is complement in the three color and two weak
modes. -/
theorem component_conjugation_binding :
    ∀ i : Fin 10,
      componentDegree (OPH.ExteriorSelection.conj i) =
        (⟨3 - (componentDegree i).1.val, by omega⟩,
          ⟨2 - (componentDegree i).2.val, by omega⟩) := by
  decide +kernel

/-- Exterior creation by the same vector twice vanishes exactly. -/
theorem creation_square_zero (v : FiveModeCarrier) (ω : FiveModeExterior) :
    ExteriorAlgebra.ι ℂ v * (ExteriorAlgebra.ι ℂ v * ω) = 0 := by
  rw [← mul_assoc, ExteriorAlgebra.ι_sq_zero, zero_mul]

/-- Exterior creation generators anticommute exactly. -/
theorem creation_generators_anticommute (v w : FiveModeCarrier) :
    ExteriorAlgebra.ι ℂ v * ExteriorAlgebra.ι ℂ w =
      -(ExteriorAlgebra.ι ℂ w * ExteriorAlgebra.ι ℂ v) :=
  eq_neg_of_add_eq_zero_left (ExteriorAlgebra.ι_add_mul_swap v w)

/-- The anticommutation law persists when the creation generators act on an
arbitrary exterior state. -/
theorem creation_actions_anticommute
    (v w : FiveModeCarrier) (ω : FiveModeExterior) :
    ExteriorAlgebra.ι ℂ v * (ExteriorAlgebra.ι ℂ w * ω) =
      -(ExteriorAlgebra.ι ℂ w * (ExteriorAlgebra.ι ℂ v * ω)) := by
  rw [← mul_assoc, creation_generators_anticommute, neg_mul, mul_assoc]

#print axioms exterior_basis_label_count
#print axioms basis_bidegree_exhaustive
#print axioms bidegree_count_table
#print axioms componentDegree_exact_nontrivial_menu
#print axioms component_dimension_binding
#print axioms component_charge_binding
#print axioms component_parity_binding
#print axioms component_conjugation_binding
#print axioms creation_square_zero
#print axioms creation_generators_anticommute
#print axioms creation_actions_anticommute

end OPH.ExteriorComponentBridge
