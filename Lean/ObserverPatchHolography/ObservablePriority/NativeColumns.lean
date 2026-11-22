import ObserverPatchHolography.ObservablePriority.StrictQuantum

/-! Source-level sparse arithmetic for the 99-state real native echo.
Blank rows are 0…8. The written rows start at 9, with coin first, signed
contrast second and the output load pair last, in Python product order.
No floating-point calculation or matrix-exponential identity is assumed. -/

namespace StrictQuantum.Native
open StrictQuantum Load
open scoped BigOperators
noncomputable section

abbrev Input := Load × Load

def n : Load → Nat
  | empty => 0 | single => 1 | double => 2

/-- Exactly the Python product ordering of SRC and LABELS in law.py. -/
def row (j : Input) (c : Bool) : Nat :=
  let a := n j.1
  let b := n j.2
  let parity := (a+b)%2
  let half := (a+b)/2
  let l := half + c.toNat*parity
  let r := half + (!c).toNat*parity
  9 * (1 + 5*c.toNat + (a+2-b)) + 3*l + r

def strict (r : Nat) : Bool :=
  9 ≤ r && r % 9 == 4 && ((r/9-1)%5 == 0 || (r/9-1)%5 == 4)

def flag (j : Input) := (sign j.1 j.2).isSome

/-- The two entries of sqrt(2)W; the row number is global in the 99-space. -/
def J (r : Nat) (j : Input) : ℤ :=
  (if row j false = r then 1 else 0) + (if row j true = r then 1 else 0)

theorem written_range : ∀ j c, 9 ≤ row j c ∧ row j c < 99 := by decide
theorem rows_separate : ∀ j k c d, row j c = row k d ↔ j = k ∧ c = d := by decide
theorem strict_on_column : ∀ j c, strict (row j c) = flag j := by decide
theorem column_overlap : ∀ i j c, J (row j c) i = if i = j then 1 else 0 := by decide
theorem hit_row : ∀ j s c, sign j.1 j.2 = some s →
    row j c = 9 * (1 + 5*c.toNat + (if s then 4 else 0)) + 4 := by decide

/-- Integer numerator of a post-selection write entry, in units of α. -/
def hitEntry (r : Nat) (j : Input) : ℤ := if strict r then J r j else 0

theorem hitEntry_exact (r : Nat) (j : Input) :
    hitEntry r j = if flag j then J r j else 0 := by
  by_cases h : row j false = r
  · subst r; simp [hitEntry, strict_on_column]
  · by_cases h' : row j true = r
    · subst r; simp [hitEntry, strict_on_column]
    · simp [hitEntry, J, h, h']

/-- Integer numerators for the quiet branch echoᵀ(write − hit). The blank
output block is scaled by 1/2. The occupied block is scaled by α/2.
These are the Wᵀ and I − WWᵀ blocks of the source's self-adjoint echo. -/
def quietBlank2 (i j : Input) : ℤ :=
  ∑ c : Bool, if strict (row j c) then 0 else J (row j c) i

def quietWritten2 (r : Nat) (j : Input) : ℤ :=
  ∑ c : Bool, if strict (row j c) then 0 else
    (if row j c = r then 2 else 0) - ∑ i : Input, J r i * J (row j c) i

theorem quietBlank_exact (i j : Input) :
    quietBlank2 i j = if flag j then 0 else if i = j then 2 else 0 := by
  by_cases hf : flag j <;> by_cases hi : i = j <;>
    simp [quietBlank2, strict_on_column, column_overlap, hf, hi]

theorem quietWritten_exact (r : Nat) (j : Input) : quietWritten2 r j = 0 := by
  by_cases hf : flag j
  · simp [quietWritten2, strict_on_column, hf]
  · simp only [quietWritten2, strict_on_column, Bool.not_eq_true] at *
    simp only [hf, Bool.false_eq_true, ↓reduceIte, column_overlap,
       mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, Fintype.sum_bool]
    simp only [J]
    split <;> split <;> omega

/-- Both local branch coefficients used by the global recursive evaluator.
A matching coin basis ket is chosen for a hit. For quiet the output is the
unchanged blank basis ket; all its written-block entries vanish above. -/
def scalar (j : Input) (wanted : Bool) : ℝ :=
  if wanted then alpha * (hitEntry (row j false) j : ℤ)
  else (quietBlank2 j j : ℤ) / 2

theorem scalar_eq (j : Input) (wanted : Bool) :
    scalar j wanted =
      if (sign j.1 j.2).isSome = wanted then
        if (sign j.1 j.2).isSome then alpha else 1 else 0 := by
  change scalar j wanted = if flag j = wanted then (if flag j then alpha else 1) else 0
  cases wanted <;> cases h : flag j <;>
    simp [scalar, hitEntry_exact, quietBlank_exact, column_overlap, h]

variable {E V : Type} [DecidableEq V]

/-- Evaluate the native integer matrix branch numerators at each service.
No precomputed matching, endpoint isometry or branch probability is an input. -/
def matrixCoefficient (q : Screen E V) (M : E → Bool) : List E → Loads V → ℝ
  | [], _ => 1
  | e :: es, x =>
    scalar (x (q.left e), x (q.right e)) (M e) *
      matrixCoefficient q M es (fire q e x)

theorem matrixCoefficient_eq (q : Screen E V) (M : E → Bool) (es : List E) (x : Loads V) :
    matrixCoefficient q M es x = branchCoefficient q M es x := by
  induction es generalizing x with
  | nil => rfl
  | cons e es ih =>
    rw [matrixCoefficient, scalar_eq, branchCoefficient, ih]
    cases h : (tag q e x).isSome <;> cases m : M e <;>
      simp_all [tag]

theorem matrix_outcome_projection (q : Screen E V) (M : E → Bool) (es : List E) (x : Loads V) :
    matrixCoefficient q M es x =
      if flags (scan q es x).2 = es.map M then weight (scan q es x).2 else 0 := by
  rw [matrixCoefficient_eq, coefficient_scan]

#print axioms written_range
#print axioms hit_row
#print axioms hitEntry_exact
#print axioms quietBlank_exact
#print axioms quietWritten_exact
#print axioms matrix_outcome_projection
#print axioms StrictQuantum.physical_injective
#print axioms StrictQuantum.kraus_eq_V_project
#print axioms StrictQuantum.branch_probability

end
end StrictQuantum.Native
