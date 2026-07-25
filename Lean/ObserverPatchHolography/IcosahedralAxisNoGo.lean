import Mathlib

/-!
# Narrow icosahedral residual-axis no-go

This module is the exact finite Lean companion to
`code/particles/runs/flavor/icosahedral_axis_angle_spectrum.json`.

The scope is deliberately narrow.  The machine receipt constructs the
antipode-identified axes of the canonical real three-dimensional regular
icosahedron: six five-fold vertex axes, ten three-fold face axes, and fifteen
two-fold edge axes.  Its pairwise acute-angle spectrum is represented here by
the exact `cos²` coordinates in `ℚ(√5)` and their multiplicities.

The comparison coordinate is the declared compare-only value
`|V_us| = 0.2250 = 9/40`.  For an acute angle, `sin² = 1 - cos²`; the theorem
`every_axis_sine_sq_exceeds_cabibbo` proves with exact algebraic inequalities
that every entry in the finite menu has a larger sine-square coordinate.

This does **not** exclude arbitrary `A5` flavor models.  In particular it says
nothing about spinorial or other representations, symmetry-breaking
corrections, or general overlap geometry.  It excludes only direct equality
between the Cabibbo angle and an acute Euclidean angle in this declared
31-axis real-three-dimensional menu.
-/

namespace OPH.IcosahedralAxisNoGo

/-- The three canonical unoriented residual-axis families. -/
inductive AxisFamily
  | fivefold
  | threefold
  | twofold
  deriving DecidableEq, Repr

/-- The exact antipode-identified family sizes constructed in the receipt. -/
def axisCount : AxisFamily → Nat
  | .fivefold => 6
  | .threefold => 10
  | .twofold => 15

/-- The receipt contains exactly `6 + 10 + 15 = 31` unoriented axes. -/
theorem axis_menu_count :
    axisCount .fivefold + axisCount .threefold + axisCount .twofold = 31 := by
  decide

/-- The six unordered/cross-family pair types in the exact spectrum. -/
inductive PairFamily
  | fivefoldFivefold
  | fivefoldThreefold
  | fivefoldTwofold
  | threefoldThreefold
  | threefoldTwofold
  | twofoldTwofold
  deriving DecidableEq, Repr

/-- The number of unordered axis pairs in each family row. -/
def expectedPairCount : PairFamily → Nat
  | .fivefoldFivefold => 15
  | .fivefoldThreefold => 60
  | .fivefoldTwofold => 90
  | .threefoldThreefold => 45
  | .threefoldTwofold => 150
  | .twofoldTwofold => 105

/-- An exact coordinate `a + b√5` in the real embedding of `ℚ(√5)`. -/
structure QsqrtFive where
  rational : ℚ
  sqrtFive : ℚ
  deriving DecidableEq, Repr

namespace QsqrtFive

/-- The real embedding used for the exact cosine-square comparison. -/
noncomputable def eval (x : QsqrtFive) : ℝ :=
  x.rational + x.sqrtFive * Real.sqrt 5

end QsqrtFive

/-- One distinct exact cosine-square value and its receipt multiplicity. -/
structure SpectrumEntry where
  family : PairFamily
  cosSq : QsqrtFive
  multiplicity : Nat
  deriving DecidableEq, Repr

/--
The complete sixteen-row exact spectrum from the receipt.

Repeated quadratic coordinates in different pair families are retained:
their multiplicities certify different finite rows of the construction.
-/
def spectrumEntry : Fin 16 → SpectrumEntry :=
  ![
    ⟨.fivefoldFivefold, ⟨1 / 5, 0⟩, 15⟩,
    ⟨.fivefoldThreefold, ⟨1 / 3, 2 / 15⟩, 30⟩,
    ⟨.fivefoldThreefold, ⟨1 / 3, -2 / 15⟩, 30⟩,
    ⟨.fivefoldTwofold, ⟨1 / 2, 1 / 10⟩, 30⟩,
    ⟨.fivefoldTwofold, ⟨1 / 2, -1 / 10⟩, 30⟩,
    ⟨.fivefoldTwofold, ⟨0, 0⟩, 30⟩,
    ⟨.threefoldThreefold, ⟨5 / 9, 0⟩, 15⟩,
    ⟨.threefoldThreefold, ⟨1 / 9, 0⟩, 30⟩,
    ⟨.threefoldTwofold, ⟨1 / 2, 1 / 6⟩, 30⟩,
    ⟨.threefoldTwofold, ⟨1 / 3, 0⟩, 60⟩,
    ⟨.threefoldTwofold, ⟨1 / 2, -1 / 6⟩, 30⟩,
    ⟨.threefoldTwofold, ⟨0, 0⟩, 30⟩,
    ⟨.twofoldTwofold, ⟨3 / 8, 1 / 8⟩, 30⟩,
    ⟨.twofoldTwofold, ⟨1 / 4, 0⟩, 30⟩,
    ⟨.twofoldTwofold, ⟨3 / 8, -1 / 8⟩, 30⟩,
    ⟨.twofoldTwofold, ⟨0, 0⟩, 15⟩
  ]

/-- Sum the exact receipt multiplicities for one pair-family row. -/
def emittedPairCount (family : PairFamily) : Nat :=
  ∑ i : Fin 16,
    if (spectrumEntry i).family = family then
      (spectrumEntry i).multiplicity
    else
      0

/--
All six spectrum rows have exactly the pair count dictated by the
`6 + 10 + 15` axis-family sizes.
-/
theorem spectrum_multiplicities :
    emittedPairCount .fivefoldFivefold =
        expectedPairCount .fivefoldFivefold ∧
      emittedPairCount .fivefoldThreefold =
        expectedPairCount .fivefoldThreefold ∧
      emittedPairCount .fivefoldTwofold =
        expectedPairCount .fivefoldTwofold ∧
      emittedPairCount .threefoldThreefold =
        expectedPairCount .threefoldThreefold ∧
      emittedPairCount .threefoldTwofold =
        expectedPairCount .threefoldTwofold ∧
      emittedPairCount .twofoldTwofold =
        expectedPairCount .twofoldTwofold := by
  native_decide

/-- Across the six rows the receipt accounts for all 465 axis pairs. -/
theorem spectrum_total_multiplicity :
    ∑ i : Fin 16, (spectrumEntry i).multiplicity = 465 := by
  native_decide

private theorem sqrt_five_sq : (Real.sqrt 5) ^ 2 = 5 :=
  Real.sq_sqrt (by norm_num)

private theorem sqrt_five_nonneg : 0 ≤ Real.sqrt 5 :=
  Real.sqrt_nonneg 5

private theorem sqrt_five_lt_nine_fourths :
    Real.sqrt 5 < (9 : ℝ) / 4 := by
  nlinarith [sqrt_five_sq]

/-- Every emitted exact coordinate is a valid cosine square. -/
theorem spectrum_coordinates_in_unit_interval (i : Fin 16) :
    0 ≤ (spectrumEntry i).cosSq.eval ∧
      (spectrumEntry i).cosSq.eval ≤ 1 := by
  fin_cases i <;>
    norm_num [spectrumEntry, QsqrtFive.eval] <;>
    constructor <;>
    nlinarith [sqrt_five_sq, sqrt_five_nonneg,
      sqrt_five_lt_nine_fourths]

/--
Exact finite comparison with `|V_us| = 9/40`: every nonzero acute angle in
the declared real-three-dimensional residual-axis menu has sine square
strictly larger than `(9/40)²`.

On `[0, π/2]` sine is increasing, so this is the algebraic coordinate form of
the narrow statement that every emitted acute axis angle exceeds the declared
Cabibbo comparison angle.
-/
theorem every_axis_sine_sq_exceeds_cabibbo (i : Fin 16) :
    ((9 : ℝ) / 40) ^ 2 <
      1 - (spectrumEntry i).cosSq.eval := by
  fin_cases i <;>
    simp [spectrumEntry, QsqrtFive.eval] <;>
    nlinarith [sqrt_five_sq, sqrt_five_nonneg,
      sqrt_five_lt_nine_fourths]

/-- The declared Cabibbo sine-square coordinate occurs nowhere in the menu. -/
theorem cabibbo_sine_sq_not_in_axis_menu (i : Fin 16) :
    1 - (spectrumEntry i).cosSq.eval ≠ ((9 : ℝ) / 40) ^ 2 :=
  ne_of_gt (every_axis_sine_sq_exceeds_cabibbo i)

end OPH.IcosahedralAxisNoGo
