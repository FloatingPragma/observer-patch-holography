import Mathlib

/-!
# The finite twelve-port Born-frame rank gap

This issue-scoped module formalizes the exact linear core of B11.  The
declared qubit adapter has six antipodal binary contexts.  Context
normalization therefore has six free weights.  Its trace-one Hermitian/Born
slice is the three-parameter family cut out by three displayed golden-ratio
relations.

The main negative theorem constructs a normalized, nonnegative context weight
with no trace-one Hermitian representation.  Thus this finite context family
does not prove a Gleason or Busch theorem.  When a representation exists, the
three coordinates are unique.

Boundary: this file does not assert that the qubit adapter is a source-produced
public quantum instrument and does not identify its coordinates with a
physical density matrix.  `HermitianRepresentable` is membership in the
declared three-coordinate image, and `DensityRepresentable` adds the declared
Bloch quadratic inequality.  The separate exact producer and independent
verifier check the actual projector coordinates, matrix interpretation, and
nonpositive-matrix witness; that operator bridge is not formalized here.
-/

namespace EventAlgebra.FiniteBornFrame

/-- The six antipodal measurement axes. -/
abbrev Axis := Fin 6

/-- A port is one of two outcomes on one of the six axes. -/
abbrev Port := Bool × Axis

/-- A real weight on the twelve labeled projectors. -/
abbrev Weight := Port → ℝ

/-- The opposite rank-one projector in the same binary context. -/
def opposite (p : Port) : Port := (!p.1, p.2)

theorem opposite_involutive (p : Port) : opposite (opposite p) = p := by
  rcases p with ⟨b, i⟩
  cases b <;> rfl

/-- Build all twelve weights from one probability per antipodal context. -/
def fromAxisWeights (a : Axis → ℝ) : Weight
  | (false, i) => a i
  | (true, i) => 1 - a i

/-- Normalization/finite additivity on each available binary context. -/
def ContextAdditive (w : Weight) : Prop :=
  ∀ i : Axis, w (false, i) + w (true, i) = 1

theorem fromAxisWeights_contextAdditive (a : Axis → ℝ) :
    ContextAdditive (fromAxisWeights a) := by
  intro i
  simp [fromAxisWeights]

/-- Every context-additive weight is exactly its six-coordinate
parameterization. -/
theorem contextAdditive_eq_fromAxisWeights {w : Weight}
    (h : ContextAdditive w) :
    w = fromAxisWeights (fun i ↦ w (false, i)) := by
  funext p
  rcases p with ⟨b, i⟩
  cases b with
  | false => rfl
  | true =>
      dsimp [fromAxisWeights]
      linarith [h i]

theorem fromAxisWeights_injective : Function.Injective fromAxisWeights := by
  intro a b h
  funext i
  have hi := congrFun h (false, i)
  simpa [fromAxisWeights] using hi

/-- Exact six-parameter classification of normalized context weights. -/
theorem contextAdditive_unique_parameterization {w : Weight}
    (h : ContextAdditive w) :
    ∃! a : Axis → ℝ, w = fromAxisWeights a := by
  refine ⟨fun i ↦ w (false, i), contextAdditive_eq_fromAxisWeights h, ?_⟩
  intro a ha
  apply fromAxisWeights_injective
  rw [← ha]
  exact contextAdditive_eq_fromAxisWeights h

/-- Finite admissibility tests only the twelve listed projector weights. -/
def Admissible (w : Weight) : Prop :=
  ContextAdditive w ∧ ∀ p : Port, 0 ≤ w p ∧ w p ≤ 1

/-- Centered coordinates on one representative from each antipodal pair. -/
def centeredWeight (w : Weight) (i : Axis) : ℝ :=
  2 * w (false, i) - 1

/-- Convert six centered values to a normalized twelve-port weight. -/
noncomputable def fromCentered (c : Axis → ℝ) : Weight :=
  fromAxisWeights (fun i ↦ (1 + c i) / 2)

theorem centeredWeight_fromCentered (c : Axis → ℝ) :
    centeredWeight (fromCentered c) = c := by
  funext i
  simp [centeredWeight, fromCentered, fromAxisWeights]
  ring

/-- The exact golden-ratio coefficient used by the declared frame. -/
noncomputable def golden : ℝ := (1 + Real.sqrt 5) / 2

/-- Coordinates of a trace-one Hermitian qubit functional, with the common
axis normalization absorbed into `(x,y,z)`. -/
@[ext] structure Coordinates where
  x : ℝ
  y : ℝ
  z : ℝ

/-- Centered weights emitted by the six representative icosahedral axes,
ordered compatibly with `PortFrameGram` and antipode `i ↦ 11-i`. -/
noncomputable def frameCentered (q : Coordinates) : Axis → ℝ :=
  ![q.y + golden * q.z,
    -q.y + golden * q.z,
    golden * q.x + q.z,
    -golden * q.x + q.z,
    q.x + golden * q.y,
    q.x - golden * q.y]

/-- The three exact linear frame relations on centered weights. -/
noncomputable def FrameRelations (c : Axis → ℝ) : Prop :=
  c 0 + c 1 = golden * (c 2 + c 3) ∧
  c 2 - c 3 = golden * (c 4 + c 5) ∧
  c 4 - c 5 = golden * (c 0 - c 1)

theorem frameRelations_frameCentered (q : Coordinates) :
    FrameRelations (frameCentered q) := by
  dsimp [FrameRelations, frameCentered]
  constructor
  · ring
  constructor <;> ring

/-- Exact tomography from centered weights. -/
noncomputable def reconstruct (c : Axis → ℝ) : Coordinates where
  x := (c 4 + c 5) / 2
  y := (c 0 - c 1) / 2
  z := (c 2 + c 3) / 2

theorem reconstruct_frameCentered (q : Coordinates) :
    reconstruct (frameCentered q) = q := by
  apply Coordinates.ext
  · simp [reconstruct, frameCentered]
  · simp [reconstruct, frameCentered]
  · simp [reconstruct, frameCentered]
    ring

theorem frameCentered_injective : Function.Injective frameCentered := by
  intro q r h
  calc
    q = reconstruct (frameCentered q) := (reconstruct_frameCentered q).symm
    _ = reconstruct (frameCentered r) := congrArg reconstruct h
    _ = r := reconstruct_frameCentered r

/-- The three relations are sufficient: the displayed tomography formula
reconstructs every centered frame value. -/
theorem frameCentered_reconstruct_of_relations (c : Axis → ℝ)
    (h : FrameRelations c) :
    frameCentered (reconstruct c) = c := by
  rcases h with ⟨h01, h23, h45⟩
  funext i
  fin_cases i
  · dsimp [frameCentered, reconstruct]
    nlinarith [h01]
  · dsimp [frameCentered, reconstruct]
    nlinarith [h01]
  · dsimp [frameCentered, reconstruct]
    nlinarith [h23]
  · dsimp [frameCentered, reconstruct]
    nlinarith [h23]
  · dsimp [frameCentered, reconstruct]
    nlinarith [h45]
  · dsimp [frameCentered, reconstruct]
    nlinarith [h45]

/-- Complete exact characterization of the three-dimensional Hermitian/Born
slice inside the six-dimensional context-weight space. -/
theorem exists_frameCentered_iff_frameRelations (c : Axis → ℝ) :
    (∃ q : Coordinates, c = frameCentered q) ↔ FrameRelations c := by
  constructor
  · rintro ⟨q, rfl⟩
    exact frameRelations_frameCentered q
  · intro h
    exact ⟨reconstruct c, (frameCentered_reconstruct_of_relations c h).symm⟩

/-- Membership in the declared three-coordinate Hermitian/Born image.  This
file does not construct the corresponding two-by-two matrix. -/
noncomputable def HermitianRepresentable (w : Weight) : Prop :=
  ∃ q : Coordinates, centeredWeight w = frameCentered q

theorem hermitianRepresentable_iff_frameRelations (w : Weight) :
    HermitianRepresentable w ↔ FrameRelations (centeredWeight w) :=
  exists_frameCentered_iff_frameRelations (centeredWeight w)

/-- Informational completeness: if a trace-one Hermitian representation
exists, its three coordinates are unique. -/
theorem hermitianRepresentation_unique {w : Weight} {q r : Coordinates}
    (hq : centeredWeight w = frameCentered q)
    (hr : centeredWeight w = frameCentered r) : q = r :=
  frameCentered_injective (hq.symm.trans hr)

/-- The declared density-positivity quadratic form in the absorbed coordinates. -/
noncomputable def densityNormSq (q : Coordinates) : ℝ :=
  (golden + 2) * (q.x ^ 2 + q.y ^ 2 + q.z ^ 2)

/-- The coordinate-level density predicate: membership in the declared image
plus the declared qubit positive-cone inequality.  The matrix/PSD bridge is an
external exact certificate boundary. -/
noncomputable def DensityRepresentable (w : Weight) : Prop :=
  ∃ q : Coordinates,
    centeredWeight w = frameCentered q ∧ densityNormSq q ≤ 1

theorem densityRepresentable_hermitianRepresentable {w : Weight}
    (h : DensityRepresentable w) : HermitianRepresentable w := by
  rcases h with ⟨q, hq, _⟩
  exact ⟨q, hq⟩

theorem densityRepresentation_unique {w : Weight} {q r : Coordinates}
    (hq : centeredWeight w = frameCentered q ∧ densityNormSq q ≤ 1)
    (hr : centeredWeight w = frameCentered r ∧ densityNormSq r ≤ 1) :
    q = r :=
  hermitianRepresentation_unique hq.1 hr.1

/-- An exact context-admissible counterexample: one binary context is
deterministic and the other five are fair. -/
def nonBornCentered (i : Axis) : ℝ := if i = 0 then 1 else 0

noncomputable def nonBornWeight : Weight := fromCentered nonBornCentered

theorem nonBornWeight_admissible : Admissible nonBornWeight := by
  constructor
  · exact fromAxisWeights_contextAdditive _
  · rintro ⟨b, i⟩
    fin_cases i <;> cases b <;>
      norm_num [nonBornWeight, nonBornCentered, fromCentered, fromAxisWeights]

/-- The admissible counterexample violates the first golden frame relation,
so it has no trace-one Hermitian representation. -/
theorem nonBornWeight_not_hermitianRepresentable :
    ¬ HermitianRepresentable nonBornWeight := by
  rw [hermitianRepresentable_iff_frameRelations]
  change ¬ FrameRelations (centeredWeight (fromCentered nonBornCentered))
  rw [centeredWeight_fromCentered]
  intro h
  have hfirst := h.1
  have h2 : (2 : Axis) ≠ 0 := by decide
  have h3 : (3 : Axis) ≠ 0 := by decide
  norm_num [nonBornCentered, h2, h3] at hfirst

/-- Consequently finite context additivity and positivity on only the twelve
listed projectors do not force a density representation. -/
theorem exists_admissible_not_densityRepresentable :
    ∃ w : Weight, Admissible w ∧ ¬ DensityRepresentable w := by
  refine ⟨nonBornWeight, nonBornWeight_admissible, ?_⟩
  intro h
  exact nonBornWeight_not_hermitianRepresentable
    (densityRepresentable_hermitianRepresentable h)

#print axioms opposite_involutive
#print axioms fromAxisWeights_contextAdditive
#print axioms contextAdditive_eq_fromAxisWeights
#print axioms fromAxisWeights_injective
#print axioms contextAdditive_unique_parameterization
#print axioms centeredWeight_fromCentered
#print axioms frameRelations_frameCentered
#print axioms reconstruct_frameCentered
#print axioms frameCentered_injective
#print axioms frameCentered_reconstruct_of_relations
#print axioms exists_frameCentered_iff_frameRelations
#print axioms hermitianRepresentable_iff_frameRelations
#print axioms hermitianRepresentation_unique
#print axioms densityRepresentable_hermitianRepresentable
#print axioms densityRepresentation_unique
#print axioms nonBornWeight_admissible
#print axioms nonBornWeight_not_hermitianRepresentable
#print axioms exists_admissible_not_densityRepresentable

end EventAlgebra.FiniteBornFrame
