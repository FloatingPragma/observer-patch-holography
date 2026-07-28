import Mathlib
import A5PortAction

namespace OPH.CommonEWOrderUnit

open OPH.A5PortAction (app comp perms)

/-! # Finite local screen/electroweak order-unit precursor

The twelve-port invariant unit and the four-copy weak multiplicity unit each
define a one-dimensional ordered line.  This file checks the actual three
carrier refinement maps and proves uniqueness of the positive unital linear
map between the two line coordinates.

Boundary: these are abstract finite order-unit lines.  No theorem here
identifies them as one physical load, selects a scalar, consumes a capacity,
or proves the N-dependent electroweak bridge. -/

/-- The constant order unit on the twelve screen ports. -/
def screenUnit : List Int := List.replicate 12 1

/-- A port permutation fixes the constant screen order unit. -/
def fixesScreenUnit (p : List Nat) : Bool :=
  p.length == 12 &&
    (List.range 12).all fun i =>
      screenUnit.getD (app p i) 0 == screenUnit.getD i 0

/-- Every listed orientation-preserving port rotation fixes the screen unit. -/
theorem all_rotations_fix_screen_unit :
    perms.all fixesScreenUnit = true := by
  decide

/-- The pinned carrier refinement map from `r0` to `r1`. -/
def refine01 : List Nat :=
  [1, 0, 3, 2, 6, 7, 4, 5, 11, 10, 9, 8]

/-- The pinned carrier refinement map from `r1` to `r2`. -/
def refine12 : List Nat :=
  [1, 5, 6, 2, 9, 0, 3, 10, 4, 8, 11, 7]

/-- The pinned direct carrier refinement map from `r0` to `r2`. -/
def refine02 : List Nat :=
  [5, 1, 2, 6, 3, 10, 9, 0, 7, 11, 8, 4]

/-- The three serialized refinement maps are permutations of the twelve ports.

The carrier artifact and `A5PortAction` use different, isomorphic port
labelings.  The certificate checks membership in the carrier artifact's
60-element rotation group.  This theorem checks the exact serialized maps in
their native labeling, without conflating the two tables. -/
theorem refinement_maps_are_permutations :
    refine01.length = 12 ∧
      (List.range 12).all refine01.contains = true ∧
      refine12.length = 12 ∧
      (List.range 12).all refine12.contains = true ∧
      refine02.length = 12 ∧
      (List.range 12).all refine02.contains = true := by
  decide

/-- The direct refinement map is exactly the composite through `r1`. -/
theorem refinement_cocycle :
    comp refine12 refine01 = refine02 := by
  decide

/-- Every pinned refinement map fixes the screen order unit. -/
theorem refinement_maps_fix_screen_unit :
    fixesScreenUnit refine01 = true ∧
      fixesScreenUnit refine12 = true ∧
      fixesScreenUnit refine02 = true := by
  decide

/-- Three quark-colour copies and one lepton copy give four weak doublets. -/
def weakDoubletCopies : Nat := 3 + 1

theorem weak_doublet_copy_count : weakDoubletCopies = 4 := by
  decide

/-- Coordinate on the one-dimensional normalized screen order-unit line. -/
abbrev ScreenLine := ℚ

/-- Coordinate on the one-dimensional normalized weak order-unit line. -/
abbrev EWLine := ℚ

/-- The normalized local line isomorphism. -/
def localIntertwiner : ScreenLine →ₗ[ℚ] EWLine :=
  LinearMap.id

/-- The line isomorphism is unital. -/
theorem local_intertwiner_unital :
    localIntertwiner 1 = 1 := by
  rfl

/-- The line isomorphism preserves the positive ray. -/
theorem local_intertwiner_positive (x : ℚ) (hx : 0 ≤ x) :
    0 ≤ localIntertwiner x := by
  simpa [localIntertwiner] using hx

/-- A unital linear map between the two one-dimensional coordinates is unique. -/
theorem unique_unital_intertwiner
    (f : ScreenLine →ₗ[ℚ] EWLine)
    (hf : f 1 = 1) :
    f = localIntertwiner := by
  apply LinearMap.ext
  intro x
  calc
    f x = f (x • (1 : ℚ)) := by simp
    _ = x • f 1 := by rw [map_smul]
    _ = x := by simp [hf]
    _ = localIntertwiner x := by rfl

/-- Normalized traces agree on the two line coordinates. -/
theorem normalized_trace_preserved (x : ScreenLine) :
    localIntertwiner x = x := by
  rfl

#print axioms all_rotations_fix_screen_unit
#print axioms refinement_maps_are_permutations
#print axioms refinement_cocycle
#print axioms refinement_maps_fix_screen_unit
#print axioms weak_doublet_copy_count
#print axioms local_intertwiner_unital
#print axioms local_intertwiner_positive
#print axioms unique_unital_intertwiner
#print axioms normalized_trace_preserved

end OPH.CommonEWOrderUnit
