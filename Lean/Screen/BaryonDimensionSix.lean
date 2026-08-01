import Mathlib

namespace OPH.BaryonDimensionSix

/-!
# Conditional dimension-six baryon-operator census

Finite arithmetic for the issue #641 packet.  The field order is

`Q, L, u_R, d_R, e_R, Q_bar, L_bar, u_R_bar, d_R_bar, e_R_bar`.

The arrays below freeze the declared one-generation Weyl grammar.  Hypercharge
is stored as `q6 = 6Y`, baryon number as `b3 = 3B`, and lepton number as an
integer.  The exhaustive theorem classifies ordered four-field multisets that
pass the exact abelian charge, color-triality, weak-parity, and Lorentz-
chirality tests.  The accompanying Python certificate expands every epsilon
pairing in an exterior algebra and proves that each surviving field multiset
has a one-dimensional nonzero contraction span.

BOUNDARY.  Baryon and lepton number are declared labels here.  This module does
not derive a physical matter attachment, Wilson coefficients, proton decay,
proton stability, a QCD matrix element, or a lifetime.
-/

/-- Six times hypercharge. -/
def q6 : Fin 10 → ℤ := ![1, -3, 4, -2, -6, -1, 3, -4, 2, 6]

/-- Three times the declared accidental baryon number. -/
def b3 : Fin 10 → ℤ := ![1, 0, 1, 1, 0, -1, 0, -1, -1, 0]

/-- Declared accidental lepton number. -/
def ell : Fin 10 → ℤ := ![0, 1, 0, 0, 1, 0, -1, 0, 0, -1]

/-- Color triality, with `-1` for an antifundamental. -/
def triality : Fin 10 → ℤ := ![1, 0, 1, 1, 0, -1, 0, -1, -1, 0]

/-- Indicator for a weak doublet. -/
def weak : Fin 10 → ℤ := ![1, 1, 0, 0, 0, 1, 1, 0, 0, 0]

/-- Indicator for an undotted left-Weyl spinor. -/
def left : Fin 10 → ℤ := ![1, 1, 0, 0, 0, 0, 0, 1, 1, 1]

def row (a b c d : Fin 10) : List (Fin 10) := [a, b, c, d]

def sumRow (f : Fin 10 → ℤ) (a b c d : Fin 10) : ℤ :=
  f a + f b + f c + f d

/-- Canonical nondecreasing field order for multiset enumeration. -/
def ordered (a b c d : Fin 10) : Prop :=
  a.val ≤ b.val ∧ b.val ≤ c.val ∧ c.val ≤ d.val

/-- The exact finite eligibility predicate used by the source certificate. -/
def admissible (a b c d : Fin 10) : Prop :=
  ordered a b c d ∧
    sumRow b3 a b c d ≠ 0 ∧
    sumRow q6 a b c d = 0 ∧
    sumRow triality a b c d % 3 = 0 ∧
    sumRow weak a b c d % 2 = 0 ∧
    sumRow left a b c d % 2 = 0

/-- The four positive-baryon rows followed by their Hermitian conjugates. -/
def admittedRows : List (List (Fin 10)) :=
  [ [0, 0, 0, 1],
    [0, 0, 2, 4],
    [0, 1, 2, 3],
    [2, 2, 3, 4],
    [5, 5, 5, 6],
    [5, 5, 7, 9],
    [5, 6, 7, 8],
    [7, 7, 8, 9] ]

/-- A finite wrapper keeps exhaustive evaluation to one `Fintype` quantifier. -/
structure Quadruple where
  a : Fin 10
  b : Fin 10
  c : Fin 10
  d : Fin 10
deriving DecidableEq, Fintype

def Quadruple.admissible (x : Quadruple) : Prop :=
  OPH.BaryonDimensionSix.admissible x.a x.b x.c x.d

def Quadruple.row (x : Quadruple) : List (Fin 10) :=
  OPH.BaryonDimensionSix.row x.a x.b x.c x.d

/-- Computable form of the eligibility predicate used for exhaustive checking. -/
def Quadruple.admissibleB (x : Quadruple) : Bool :=
  decide (x.a.val ≤ x.b.val) &&
    decide (x.b.val ≤ x.c.val) &&
    decide (x.c.val ≤ x.d.val) &&
    (sumRow b3 x.a x.b x.c x.d != 0) &&
    (sumRow q6 x.a x.b x.c x.d == 0) &&
    (sumRow triality x.a x.b x.c x.d % 3 == 0) &&
    (sumRow weak x.a x.b x.c x.d % 2 == 0) &&
    (sumRow left x.a x.b x.c x.d % 2 == 0)

/-- The executable classifier is definitionally tied to the stated predicate. -/
theorem admissibleB_iff (x : Quadruple) : x.admissibleB = true ↔ x.admissible := by
  simp [Quadruple.admissibleB, Quadruple.admissible,
    OPH.BaryonDimensionSix.admissible, ordered]
  tauto

/-- Computable implication: every admitted row has the displayed charge law. -/
def Quadruple.chargeLawB (x : Quadruple) : Bool :=
  !x.admissibleB ||
    ((sumRow q6 x.a x.b x.c x.d == 0) &&
      (Int.natAbs (sumRow b3 x.a x.b x.c x.d) == 3) &&
      (Int.natAbs (sumRow ell x.a x.b x.c x.d) == 1) &&
      (sumRow b3 x.a x.b x.c x.d - 3 * sumRow ell x.a x.b x.c x.d == 0))

/-- Joint finite check, split below into ten kernel-sized prefix blocks. -/
def Quadruple.packetB (x : Quadruple) : Bool :=
  (x.admissibleB == admittedRows.contains x.row) && x.chargeLawB

structure Triple where
  b : Fin 10
  c : Fin 10
  d : Fin 10
deriving DecidableEq, Fintype

def withPrefix (a : Fin 10) (x : Triple) : Quadruple := ⟨a, x.b, x.c, x.d⟩

set_option maxRecDepth 8192 in
theorem packet_prefix_0 : ∀ x : Triple, (withPrefix 0 x).packetB = true := by decide +kernel
set_option maxRecDepth 8192 in
theorem packet_prefix_1 : ∀ x : Triple, (withPrefix 1 x).packetB = true := by decide +kernel
set_option maxRecDepth 8192 in
theorem packet_prefix_2 : ∀ x : Triple, (withPrefix 2 x).packetB = true := by decide +kernel
set_option maxRecDepth 8192 in
theorem packet_prefix_3 : ∀ x : Triple, (withPrefix 3 x).packetB = true := by decide +kernel
set_option maxRecDepth 8192 in
theorem packet_prefix_4 : ∀ x : Triple, (withPrefix 4 x).packetB = true := by decide +kernel
set_option maxRecDepth 8192 in
theorem packet_prefix_5 : ∀ x : Triple, (withPrefix 5 x).packetB = true := by decide +kernel
set_option maxRecDepth 8192 in
theorem packet_prefix_6 : ∀ x : Triple, (withPrefix 6 x).packetB = true := by decide +kernel
set_option maxRecDepth 8192 in
theorem packet_prefix_7 : ∀ x : Triple, (withPrefix 7 x).packetB = true := by decide +kernel
set_option maxRecDepth 8192 in
theorem packet_prefix_8 : ∀ x : Triple, (withPrefix 8 x).packetB = true := by decide +kernel
set_option maxRecDepth 8192 in
theorem packet_prefix_9 : ∀ x : Triple, (withPrefix 9 x).packetB = true := by decide +kernel

theorem packet_all (x : Quadruple) : x.packetB = true := by
  rcases x with ⟨a, b, c, d⟩
  fin_cases a
  · exact packet_prefix_0 ⟨b, c, d⟩
  · exact packet_prefix_1 ⟨b, c, d⟩
  · exact packet_prefix_2 ⟨b, c, d⟩
  · exact packet_prefix_3 ⟨b, c, d⟩
  · exact packet_prefix_4 ⟨b, c, d⟩
  · exact packet_prefix_5 ⟨b, c, d⟩
  · exact packet_prefix_6 ⟨b, c, d⟩
  · exact packet_prefix_7 ⟨b, c, d⟩
  · exact packet_prefix_8 ⟨b, c, d⟩
  · exact packet_prefix_9 ⟨b, c, d⟩

/-- Exhaustive classification of all `C(13,4) = 715` ordered multisets. -/
theorem exhaustive_admissible_quadruples (x : Quadruple) :
    x.admissibleB = admittedRows.contains x.row := by
  have h := packet_all x
  simp only [Quadruple.packetB, Bool.and_eq_true, beq_iff_eq] at h
  exact h.1

theorem exhaustive_admissible_rows (a b c d : Fin 10) :
    (Quadruple.mk a b c d).admissibleB = admittedRows.contains (row a b c d) :=
  exhaustive_admissible_quadruples ⟨a, b, c, d⟩

/-- Every admitted row changes `B` and `L` by the same signed unit. -/
theorem admitted_charge_law_quadruples (x : Quadruple) : x.chargeLawB = true := by
  have h := packet_all x
  simp only [Quadruple.packetB, Bool.and_eq_true, beq_iff_eq] at h
  exact h.2

theorem admitted_charge_law (a b c d : Fin 10) :
    (Quadruple.mk a b c d).chargeLawB = true :=
  admitted_charge_law_quadruples ⟨a, b, c, d⟩

/-- The diagonal quotient generator fixes every field separately. -/
def diagonalPhase (i : Fin 10) : ℤ :=
  (2 * triality i + 3 * weak i + q6 i) % 6

theorem every_field_descends_through_diagonal_Z6 :
    ∀ i : Fin 10, diagonalPhase i = 0 := by
  decide

/-- Hermitian conjugation on the frozen field order. -/
def conjugate : Fin 10 → Fin 10 := ![5, 6, 7, 8, 9, 0, 1, 2, 3, 4]

theorem conjugation_involutive : ∀ i : Fin 10, conjugate (conjugate i) = i := by
  decide

/-- Two-index alternating symbol used by the component witnesses. -/
def eps2 (i j : Fin 2) : ℤ :=
  if i = 0 ∧ j = 1 then 1 else if i = 1 ∧ j = 0 then -1 else 0

/-- One fixed nonzero orientation of the three-index alternating symbol. -/
def eps3 (i j k : Fin 3) : ℤ :=
  if (i, j, k) = (0, 1, 2) then 1 else 0

/-- A nonzero tensor factor occurs in the displayed `QQQL` contraction.
The exterior-algebra certificate checks its fully aggregated coefficient. -/
theorem qqql_epsilon_factor_nonzero :
    eps3 0 1 2 * eps2 0 1 * eps2 0 1 * eps2 0 1 * eps2 0 1 ≠ 0 := by
  decide

/-- A nonzero tensor factor occurs in the displayed `DUUE` contraction.
The exterior-algebra certificate includes identical-`u_R` Pauli signs and
checks that the fully aggregated one-generation polynomial remains nonzero. -/
theorem duue_epsilon_factor_nonzero :
    eps3 0 1 2 * eps2 0 1 * eps2 0 1 ≠ 0 := by
  decide

#print axioms exhaustive_admissible_quadruples
#print axioms exhaustive_admissible_rows
#print axioms admissibleB_iff
#print axioms admitted_charge_law_quadruples
#print axioms admitted_charge_law
#print axioms every_field_descends_through_diagonal_Z6
#print axioms conjugation_involutive
#print axioms qqql_epsilon_factor_nonzero
#print axioms duue_epsilon_factor_nonzero

end OPH.BaryonDimensionSix
