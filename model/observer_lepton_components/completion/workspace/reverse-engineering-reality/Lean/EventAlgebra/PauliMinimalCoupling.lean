import EventAlgebra.Robertson
import Mathlib.Tactic.Module

/-!
# Pauli minimal-coupling square: the exact algebraic origin of gyromagnetic ratio two

This module proves, over an arbitrary complex algebra `R` (noncommutativity
permitted and required for nonzero content), the exact Clifford-square identity
for the Pauli contraction of a declared triple `X : Fin 3 → R`:

  `(σ · X)² = (X · X) • 1 + I • σ · (X × X)`

where `(X × X) k` is the Levi-Civita contraction of `X i * X j`.  Under the
declared commutation hypothesis `[X i, X j] = (I * q * ε i j k) • B k`
(`q : ℂ` a declared coupling scalar, `B : Fin 3 → R` a declared triple), the
cross term collapses exactly, `(X × X) k = (I * q) • B k`, and the square
becomes the Pauli identity

  `(σ · X)² = (X · X) • 1 - q • (σ · B)`

with coefficient exactly `q`.  The degenerate receipt shows the spin term
vanishes identically for commuting `X` (the `B = 0` branch of the hypothesis
class), a concrete inhabitant (`R = Matrix (Fin 2) (Fin 2) ℂ`, `X = B` the
committed Pauli matrices themselves, `q = 2`) shows the hypothesis class is
nonempty with nonzero spin term, and a discrimination lemma shows the
inhabitant pins the coefficient: no `c ≠ q` satisfies the same identity.  One
composed receipt derives every clause from a single typed antecedent bundle.

## Formal precursor clause (gyromagnetic reading; report, not theorem)

In the standard convention where the orbital sector of a minimally coupled
charge couples with weight `q/(2m)` per unit orbital angular momentum,
dividing the Pauli square by `2m` yields the spin coupling `(q/(2m)) • (σ · B)`
against the spin operator `σ/2`: the derived moment per unit spin angular
momentum is `q/m`, twice the orbital weight `q/(2m)`.  That factor two is the
statement `g = 2`, and this module locates it entirely in the Clifford square
above: the coefficient of `σ · B` is exactly `q` (never `q/2`).  No physical
field, electron, mass value, gauge potential, Dirac equation, or measured
moment enters; the reading is a labeling of the algebraic identity, recorded
here as a formal precursor for the magnetic-moment observation row.

## Declared data and nonclaims

The commutation hypothesis, the coupling scalar `q`, the triples `X` and `B`,
and the carrier algebra `R` are declared inputs of every statement; none is
derived from observer repair data.  No gauge potential is constructed (the
source gauge-field, current, and action attachment premise is open), no
relativistic or Dirac-sector derivation is performed, no anomalous-moment
(`g - 2`) content exists here, and no laboratory unit, calibration anchor, or
measured moment is attached (the clock/energy calibration,
observer-to-physical-spacetime, and source gauge-field premises are open).
-/

namespace EventAlgebra.PauliMinimalCoupling

open Finset

/-! ## The Pauli family and the Levi-Civita symbol -/

/-- The three committed Pauli matrices of `EventAlgebra.Robertson` as one
indexed family: `pauli 0 = X`, `pauli 1 = Y`, `pauli 2 = Z`. -/
def pauli : Fin 3 → Matrix (Fin 2) (Fin 2) ℂ :=
  ![Robertson.pauliX, Robertson.pauliY, Robertson.pauliZ]

/-- The Levi-Civita symbol on `Fin 3`, valued in `ℂ`. -/
def eps : Fin 3 → Fin 3 → Fin 3 → ℂ
  | 0, 1, 2 => 1
  | 1, 2, 0 => 1
  | 2, 0, 1 => 1
  | 2, 1, 0 => -1
  | 1, 0, 2 => -1
  | 0, 2, 1 => -1
  | _, _, _ => 0

variable {R : Type*} [Ring R] [Algebra ℂ R]

/-! ## The sigma contraction, the dot square, and the epsilon cross square -/

/-- The Pauli contraction `σ · X` of a declared triple `X : Fin 3 → R`: the
`2 × 2` matrix over `R` whose `(a, b)` entry is `∑ k, (pauli k) a b • X k`,
the scalar Pauli entries acting through the `ℂ`-algebra structure of `R`. -/
def sigmaDot (X : Fin 3 → R) : Matrix (Fin 2) (Fin 2) R :=
  Matrix.of fun a b => ∑ k, pauli k a b • X k

/-- The scalar square `X · X = ∑ k, X k * X k`. -/
def dotSelf (X : Fin 3 → R) : R := ∑ k, X k * X k

/-- The epsilon cross square `(X × X) k = ∑ i j, ε i j k • (X i * X j)`.
For commuting entries this vanishes; in general it collects the
commutators. -/
def crossSelf (X : Fin 3 → R) (k : Fin 3) : R :=
  ∑ i, ∑ j, eps i j k • (X i * X j)

/-- Closed form of the sigma contraction. -/
theorem sigmaDot_eq (X : Fin 3 → R) :
    sigmaDot X =
      Matrix.of
        ![![X 2, X 0 - Complex.I • X 1],
          ![X 0 + Complex.I • X 1, -X 2]] := by
  ext a b
  fin_cases a <;> fin_cases b <;>
    simp [sigmaDot, pauli, Robertson.pauliX, Robertson.pauliY,
      Robertson.pauliZ, Fin.sum_univ_three, sub_eq_add_neg, neg_smul]

/-- Closed form of the epsilon cross square: the three commutators. -/
theorem crossSelf_eq (X : Fin 3 → R) :
    crossSelf X =
      ![X 1 * X 2 - X 2 * X 1, X 2 * X 0 - X 0 * X 2,
        X 0 * X 1 - X 1 * X 0] := by
  funext k
  fin_cases k <;>
    (simp [crossSelf, eps, Fin.sum_univ_three, sub_eq_add_neg, neg_smul];
      try abel)

/-! ## Clause 1: the exact Clifford square identity -/

/-- **Clifford square.** For every declared triple `X` over every complex
algebra `R`, `(σ · X)² = (X · X) • 1 + I • σ · (X × X)` exactly.  No
commutation hypothesis enters. -/
theorem sigmaDot_mul_self (X : Fin 3 → R) :
    sigmaDot X * sigmaDot X =
      dotSelf X • (1 : Matrix (Fin 2) (Fin 2) R) +
        Complex.I • sigmaDot (crossSelf X) := by
  rw [sigmaDot_eq X, sigmaDot_eq (crossSelf X), crossSelf_eq X]
  ext a b
  fin_cases a <;> fin_cases b <;>
    simp [Matrix.mul_apply, Fin.sum_univ_two, dotSelf,
      Fin.sum_univ_three, mul_add, add_mul, mul_sub, sub_mul,
      smul_smul, Complex.I_mul_I, smul_sub, smul_add,
      sub_neg_eq_add, mul_neg, smul_eq_mul] <;>
    abel

/-! ## Clause 2: the commutation hypothesis collapses the cross square -/

/-- The declared commutation hypothesis: `X i * X j - X j * X i` equals the
epsilon-weighted combination `∑ k, (I * q * ε i j k) • B k` for every pair of
indices.  `q : ℂ` is the declared coupling scalar acting through the
`ℂ`-algebra structure; `X` and `B` are declared triples in `R`.  Nothing here
derives the hypothesis; consumers declare it. -/
def CommutationHypothesis (q : ℂ) (X B : Fin 3 → R) : Prop :=
  ∀ i j : Fin 3,
    X i * X j - X j * X i = ∑ k, (Complex.I * q * eps i j k) • B k

/-- **Cross collapse.** Under the commutation hypothesis the epsilon cross
square is exactly `(X × X) k = (I * q) • B k`. -/
theorem crossSelf_of_comm {q : ℂ} {X B : Fin 3 → R}
    (h : CommutationHypothesis q X B) (k : Fin 3) :
    crossSelf X k = (Complex.I * q) • B k := by
  have h12 := h 1 2
  have h20 := h 2 0
  have h01 := h 0 1
  simp only [eps, Fin.sum_univ_three, mul_one, mul_zero,
    zero_smul, add_zero, zero_add] at h12 h20 h01
  rw [crossSelf_eq]
  fin_cases k
  · simpa using h12
  · simpa using h20
  · simpa using h01

/-! ## Clause 3: the Pauli square identity with coefficient exactly `q` -/

/-- The sigma contraction is `ℂ`-homogeneous in the triple. -/
theorem sigmaDot_smul (c : ℂ) (X : Fin 3 → R) :
    sigmaDot (fun k => c • X k) = c • sigmaDot X := by
  ext a b
  simp only [sigmaDot, Matrix.of_apply, Matrix.smul_apply, smul_sum]
  exact Finset.sum_congr rfl fun k _ => (smul_comm c (pauli k a b) (X k)).symm

/-- **Pauli square.** Under the declared commutation hypothesis,
`(σ · X)² = (X · X) • 1 - q • (σ · B)` with coefficient exactly `q`: the
spin coupling inherits the full declared coupling scalar, with no factor
one-half. -/
theorem pauliSquare {q : ℂ} {X B : Fin 3 → R}
    (h : CommutationHypothesis q X B) :
    sigmaDot X * sigmaDot X =
      dotSelf X • (1 : Matrix (Fin 2) (Fin 2) R) - q • sigmaDot B := by
  have hc : crossSelf X = fun k => (Complex.I * q) • B k :=
    funext fun k => crossSelf_of_comm h k
  rw [sigmaDot_mul_self, hc, sigmaDot_smul, smul_smul, ← mul_assoc,
    Complex.I_mul_I, neg_one_mul, neg_smul, ← sub_eq_add_neg]

/-! ## Clause 4: the degenerate receipt -/

/-- For a commuting triple the epsilon cross square vanishes identically. -/
theorem crossSelf_eq_zero_of_commuting {X : Fin 3 → R}
    (h : ∀ i j, X i * X j = X j * X i) : crossSelf X = 0 := by
  rw [crossSelf_eq]
  funext k
  fin_cases k <;> simp [h 1 2, h 2 0, h 0 1]

/-- The sigma contraction of the zero triple is zero. -/
theorem sigmaDot_zero : sigmaDot (0 : Fin 3 → R) = 0 := by
  ext a b
  simp [sigmaDot]

/-- **Degenerate receipt.** For a commuting triple the spin term vanishes
identically and the Clifford square is the pure scalar square. -/
theorem sigmaDot_mul_self_of_commuting {X : Fin 3 → R}
    (h : ∀ i j, X i * X j = X j * X i) :
    sigmaDot X * sigmaDot X =
      dotSelf X • (1 : Matrix (Fin 2) (Fin 2) R) := by
  rw [sigmaDot_mul_self, crossSelf_eq_zero_of_commuting h, sigmaDot_zero,
    smul_zero, add_zero]

/-- A commuting triple satisfies the commutation hypothesis with `B = 0` for
every declared coupling: commuting `X` is exactly the `B = 0` branch. -/
theorem commutationHypothesis_of_commuting (q : ℂ) {X : Fin 3 → R}
    (h : ∀ i j, X i * X j = X j * X i) :
    CommutationHypothesis q X (0 : Fin 3 → R) := by
  intro i j
  simp [h i j]

/-! ## Clause 5: a concrete inhabitant with nonzero spin term -/

/-- **Inhabitant.** The committed Pauli matrices themselves, viewed as a
noncommuting triple in the complex algebra `Matrix (Fin 2) (Fin 2) ℂ`,
satisfy the commutation hypothesis with `B = X = pauli` and coupling `q = 2`:
`[σ_i, σ_j] = 2 I ε_{ijk} σ_k`.  The hypothesis class is nonempty. -/
theorem pauli_commutationHypothesis :
    CommutationHypothesis (R := Matrix (Fin 2) (Fin 2) ℂ) 2 pauli pauli := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    · ext a b
      fin_cases a <;> fin_cases b <;>
        simp [pauli, Robertson.pauliX, Robertson.pauliY, Robertson.pauliZ,
          Matrix.sub_apply, Matrix.smul_apply,
          Matrix.sum_apply, Fin.sum_univ_three, eps,
          smul_eq_mul] <;>
        try ring_nf <;>
        try norm_num [Complex.I_sq]

/-- The single typed antecedent bundle: a carrier triple, a field triple, a
coupling scalar, and the declared commutation hypothesis tying them. -/
structure MinimalCouplingData (R : Type*) [Ring R] [Algebra ℂ R] where
  /-- The declared coupling scalar. -/
  q : ℂ
  /-- The declared carrier triple (kinetic side of the contraction). -/
  X : Fin 3 → R
  /-- The declared commutator-target triple (field side). -/
  B : Fin 3 → R
  /-- The declared commutation hypothesis `[X i, X j] = (I q ε i j k) • B k`. -/
  comm : CommutationHypothesis q X B

/-- The concrete inhabitant packaged as a bundle: the Pauli triple over
`Matrix (Fin 2) (Fin 2) ℂ` with `q = 2` and `B = X`. -/
def pauliTriple : MinimalCouplingData (Matrix (Fin 2) (Fin 2) ℂ) :=
  ⟨2, pauli, pauli, pauli_commutationHypothesis⟩

/-- The inhabitant's spin term is nonzero: `q • (σ · B)` with `q = 2` and
`B = pauli` has a nonvanishing entry, so the Pauli square identity carries
nonzero spin content on this inhabitant. -/
theorem pauliTriple_spin_term_ne_zero :
    pauliTriple.q • sigmaDot pauliTriple.B ≠ 0 := by
  intro hzero
  have h := congrArg (fun M => M 0 0 0 0) hzero
  norm_num [pauliTriple, sigmaDot, pauli, Robertson.pauliX, Robertson.pauliY,
    Robertson.pauliZ, Fin.sum_univ_three, Matrix.cons_val_two,
    Matrix.tail_cons, Matrix.head_cons] at h

/-- **Coefficient discrimination.** On the inhabitant the Pauli square pins
the coupling coefficient: any scalar `c` satisfying the same identity equals
`2` exactly.  In particular the identity with coefficient `q / 2` or `2 q`
fails on the inhabitant, which is the algebraic content of the factor-two
(gyromagnetic) reading. -/
theorem pauliTriple_coefficient_unique (c : ℂ)
    (hc : sigmaDot pauliTriple.X * sigmaDot pauliTriple.X =
      dotSelf pauliTriple.X • (1 : Matrix (Fin 2) (Fin 2) (Matrix (Fin 2) (Fin 2) ℂ)) -
        c • sigmaDot pauliTriple.B) : c = 2 := by
  have h2 := pauliSquare pauliTriple.comm
  have hdiff : c • sigmaDot pauliTriple.B =
      pauliTriple.q • sigmaDot pauliTriple.B :=
    sub_right_injective (hc.symm.trans h2)
  have hentry := congrArg (fun M => M 0 0 0 0) hdiff
  simpa [pauliTriple, sigmaDot, pauli, Robertson.pauliX, Robertson.pauliY,
    Robertson.pauliZ, Fin.sum_univ_three] using hentry

/-! ## Composed receipt from the single typed antecedent bundle -/

/-- **Composed receipt.** From one typed antecedent bundle: the exact
Clifford square, the exact cross collapse, the Pauli square with coefficient
exactly `q`, and the degenerate branch for commuting carriers. -/
theorem minimalCoupling_receipt (D : MinimalCouplingData R) :
    (sigmaDot D.X * sigmaDot D.X =
        dotSelf D.X • (1 : Matrix (Fin 2) (Fin 2) R) +
          Complex.I • sigmaDot (crossSelf D.X)) ∧
    (∀ k, crossSelf D.X k = (Complex.I * D.q) • D.B k) ∧
    (sigmaDot D.X * sigmaDot D.X =
        dotSelf D.X • (1 : Matrix (Fin 2) (Fin 2) R) - D.q • sigmaDot D.B) ∧
    ((∀ i j, D.X i * D.X j = D.X j * D.X i) →
        sigmaDot D.X * sigmaDot D.X =
          dotSelf D.X • (1 : Matrix (Fin 2) (Fin 2) R)) :=
  ⟨sigmaDot_mul_self D.X, crossSelf_of_comm D.comm, pauliSquare D.comm,
    fun h => sigmaDot_mul_self_of_commuting h⟩

#print axioms sigmaDot_mul_self
#print axioms crossSelf_of_comm
#print axioms pauliSquare
#print axioms sigmaDot_mul_self_of_commuting
#print axioms commutationHypothesis_of_commuting
#print axioms pauli_commutationHypothesis
#print axioms pauliTriple_spin_term_ne_zero
#print axioms pauliTriple_coefficient_unique
#print axioms minimalCoupling_receipt

end EventAlgebra.PauliMinimalCoupling
