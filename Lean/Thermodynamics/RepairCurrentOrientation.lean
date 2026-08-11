import Mathlib

/-!
# The repair-load irreversibility current and its orientation

Literal mirror of `docs/REPAIR_CURRENT_PAYLOAD.json` of `oph-physics-sim`
(payload schema v2, sha256
`a9214bd110a1e8808c16202f6970dd90fd0e6c63ac6aa58acdc181b2da39d2c1`,
extracted from the retained preregistered bundle
`runs/b12_prereg_16k_20260806` by `scripts/extract_repair_current.py`):
the exact ordered transition counts of the committed run's
transition-history windows, projected onto the declared
`repair_load_bucket` coordinate, with the two declared orientation
invariants.

The mathematical content is small and exact: the ordered count table is
not symmetric, the designated pair `3 -> 4` carries `1343` forward counts
and `0` backward counts, the designated 3-cycle `(3,4,5)` has forward
count product `1239691068` and backward count product `0`, and reversing
every counted window transposes the table, so both invariants change sign
under time reversal of the counted order.  The generic lemmas state the
reversal behaviour for arbitrary tables; the literal receipts instantiate
them on the committed counts.  The designation rule itself is certified
by kernel decides: `1343` and `1239691068` are the exact maxima of the
pair asymmetry and the cycle gap, each attained only on the designated
orbit, whose lexicographically least representative is the designated
pair or cycle by inspection of the listed orbit.  The transposition of
the table under window reversal is an identity of ordered recounting,
recorded by the extractor as reversal semantics rather than as a data
check.

**Boundary.**  These are ordered counts of one committed bounded run
under one declared quotient.  No physical arrow of time, continuum
current, thermodynamic-limit statement, or laboratory claim is made; the
physical clock is owned by E5.  The orientation bit extracted here is
consumed by the B13 orientation-selection packet
(`QFT/SourceOrientedCompletion.lean`) under a declared typed convention
stated there.
-/

set_option maxRecDepth 8192

namespace OPH.Thermodynamics

/-- The exact ordered repair-load transition counts of the committed run:
row = source bucket, column = target bucket. -/
def repairCounts : Fin 8 → Fin 8 → ℕ := fun a =>
  ![![23552, 0, 0, 0, 0, 0, 0, 0],
    ![0, 2, 32, 0, 1, 5, 0, 3],
    ![0, 0, 0, 167, 0, 5, 0, 0],
    ![0, 28, 0, 349, 1343, 0, 46, 0],
    ![0, 0, 135, 0, 0, 1554, 0, 144],
    ![0, 1, 0, 594, 261, 227, 1246, 3],
    ![1024, 10, 5, 0, 205, 93, 78, 71],
    ![0, 2, 0, 283, 23, 81, 116, 55]] a

/-- The counts account for every counted transition of the run. -/
theorem repairCounts_total :
    (∑ a : Fin 8, ∑ b : Fin 8, repairCounts a b) = 31744 := by decide

/-- The ordered table is not symmetric: the counted order carries
information. -/
theorem repairCounts_not_symmetric :
    ¬ (∀ a b : Fin 8, repairCounts a b = repairCounts b a) := by decide

/-! ## Generic reversal behaviour -/

/-- The signed pair asymmetry of an ordered count table. -/
def pairAsymmetry (C : Fin 8 → Fin 8 → ℕ) (a b : Fin 8) : ℤ :=
  (C a b : ℤ) - (C b a : ℤ)

/-- The forward count product around an ordered 3-cycle. -/
def cycleForward (C : Fin 8 → Fin 8 → ℕ) (a b c : Fin 8) : ℕ :=
  C a b * C b c * C c a

/-- Reversal of the counted order is transposition of the table. -/
def reversal (C : Fin 8 → Fin 8 → ℕ) : Fin 8 → Fin 8 → ℕ :=
  fun a b => C b a

theorem reversal_pairAsymmetry (C : Fin 8 → Fin 8 → ℕ) (a b : Fin 8) :
    pairAsymmetry (reversal C) a b = -pairAsymmetry C a b := by
  simp [pairAsymmetry, reversal]

theorem reversal_cycleForward (C : Fin 8 → Fin 8 → ℕ) (a b c : Fin 8) :
    cycleForward (reversal C) a b c = cycleForward C c b a := by
  simp only [cycleForward, reversal]
  ring

theorem reversal_involutive (C : Fin 8 → Fin 8 → ℕ) :
    reversal (reversal C) = C := rfl

/-! ## The designated literal receipts -/

/-- The designated pair `3 -> 4`: 1343 forward, 0 backward. -/
theorem designatedPair_counts :
    repairCounts 3 4 = 1343 ∧ repairCounts 4 3 = 0 := by decide

/-- The designated pair asymmetry is exactly `1343`. -/
theorem designatedPair_asymmetry :
    pairAsymmetry repairCounts 3 4 = 1343 := by decide

/-- The designated 3-cycle `(3,4,5)`: forward product `1239691068`,
backward product `0`. -/
theorem designatedCycle_products :
    cycleForward repairCounts 3 4 5 = 1239691068 ∧
      cycleForward repairCounts 5 4 3 = 0 := by decide

/-- `1343` is the exact maximum absolute pair asymmetry of the table. -/
theorem designatedPair_maximal :
    ∀ a b : Fin 8, (pairAsymmetry repairCounts a b).natAbs ≤ 1343 := by
  decide

/-- The maximum pair asymmetry is attained only on the designated pair
and its reverse; `(3, 4)` is the lexicographically least maximizer. -/
theorem designatedPair_attainment :
    ∀ a b : Fin 8, (pairAsymmetry repairCounts a b).natAbs = 1343 →
      (a = 3 ∧ b = 4) ∨ (a = 4 ∧ b = 3) := by decide

/-- The signed cycle gap of an ordered count table. -/
def cycleGap (C : Fin 8 → Fin 8 → ℕ) (a b c : Fin 8) : ℤ :=
  (cycleForward C a b c : ℤ) - (cycleForward C c b a : ℤ)

/-- `1239691068` is the exact maximum absolute cycle gap of the table. -/
theorem designatedCycle_maximal :
    ∀ a b c : Fin 8,
      (cycleGap repairCounts a b c).natAbs ≤ 1239691068 := by decide

/-- The maximum cycle gap is attained only on the six representations of
the designated cycle `{3, 4, 5}`; `(3, 4, 5)` is the lexicographically
least maximizer. -/
theorem designatedCycle_attainment :
    ∀ a b c : Fin 8,
      (cycleGap repairCounts a b c).natAbs = 1239691068 →
        (a = 3 ∧ b = 4 ∧ c = 5) ∨ (a = 4 ∧ b = 5 ∧ c = 3) ∨
        (a = 5 ∧ b = 3 ∧ c = 4) ∨ (a = 5 ∧ b = 4 ∧ c = 3) ∨
        (a = 4 ∧ b = 3 ∧ c = 5) ∨ (a = 3 ∧ b = 5 ∧ c = 4) := by decide

/-- **The source orientation bit**: the committed counted order runs the
designated cycle forward.  Time reversal flips it. -/
def repairOrientationBit : Bool :=
  decide (cycleForward repairCounts 5 4 3 < cycleForward repairCounts 3 4 5)

theorem repairOrientationBit_true : repairOrientationBit = true := by decide

/-- Reversing the counted order flips the orientation bit. -/
theorem reversal_flips_orientation :
    decide (cycleForward (reversal repairCounts) 5 4 3 <
        cycleForward (reversal repairCounts) 3 4 5) = false := by decide

/-! ## Negative control: a reversible table carries no orientation -/

/-- A detailed-balanced control table: symmetric counts. -/
def reversibleControl : Fin 8 → Fin 8 → ℕ := fun a b =>
  (a : ℕ) + (b : ℕ) + 1

theorem reversibleControl_symmetric :
    ∀ a b : Fin 8, reversibleControl a b = reversibleControl b a := by
  intro a b
  simp [reversibleControl]
  omega

/-- The control's cycle products agree in both orders: no orientation is
selected. -/
theorem reversibleControl_no_orientation :
    cycleForward reversibleControl 3 4 5 =
      cycleForward reversibleControl 5 4 3 := by decide

end OPH.Thermodynamics

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.Thermodynamics.repairCounts_total
#print axioms OPH.Thermodynamics.repairCounts_not_symmetric
#print axioms OPH.Thermodynamics.designatedPair_asymmetry
#print axioms OPH.Thermodynamics.designatedCycle_products
#print axioms OPH.Thermodynamics.designatedPair_maximal
#print axioms OPH.Thermodynamics.designatedPair_attainment
#print axioms OPH.Thermodynamics.designatedCycle_maximal
#print axioms OPH.Thermodynamics.designatedCycle_attainment
#print axioms OPH.Thermodynamics.repairOrientationBit_true
#print axioms OPH.Thermodynamics.reversal_flips_orientation
#print axioms OPH.Thermodynamics.reversibleControl_no_orientation
