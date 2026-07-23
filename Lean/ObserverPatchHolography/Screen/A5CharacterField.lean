import Mathlib

namespace OPH.A5CharacterField

/-! # The A5 character field ℚ(√5): Galois stability and the centre dimension list

Issue #605.  The trichotomy proof's rationality lemma states that the centre
of the compact algebra is a ℚ-rational `A5`-submodule (the Lie algebra of a
torus whose integral cocharacter lattice is preserved by the action), and
that the Galois conjugation of ℚ(√5) then swaps the two three-dimensional
irreducibles, forcing equal multiplicity in the centre.  The
torus/cocharacter step is the declared classical hypothesis, exactly as the
issue permits; this file proves the Galois-stability half and the resulting
dimension list.

DECLARED INPUT.  The `A5` character table on the five conjugacy classes
(identity, double transpositions, three-cycles, and the two five-cycle
classes), rows ordered `1, 3, 3', 4, 5`:

    1  :  1   1   1   1        1
    3  :  3  -1   0   φ        φ'
    3' :  3  -1   0   φ'       φ
    4  :  4   0   1  -1       -1
    5  :  5   1  -1   0        0

with `φ = (1+√5)/2` and `φ' = (1-√5)/2`.  The table is carried in DOUBLED
form so every entry lies in `ℤ√5`: doubling is additive and injective, so
every multiplicity statement is unchanged, and the Galois conjugation of
ℚ(√5) restricts to the star involution `a + b√5 ↦ a - b√5` of `ℤ√5`
(`Zsqrtd.im_star`).

RESULTS.

* `star_chi2_one` / `star_chi2_two`: the conjugation swaps the doubled
  characters of `3` and `3'` and fixes the other three rows
  (kernel-`decide` on the concrete table).
* `multiplicities_equal_of_galoisStable`: a module whose combined character
  is Galois-stable carries `3` and `3'` with equal multiplicity.  The proof
  evaluates the imaginary part at a five-cycle class.
* `dim_eq_char_at_identity`: twice the dimension is the doubled combined
  character at the identity class (consistency of the dimension readout).
* `centreDim_cases` / `centreDim_mem_trichotomy_list` /
  `centreDim_ne_two`: inside the twelve-port module `1 + 3 + 3' + 5`
  (multiplicity bounds 1, 1, 1, 0, 1), a Galois-stable submodule has
  dimension in `{0, 1, 5, 6, 7, 11, 12}`; dimension two is impossible.

The statements map to the rationality lemma of the compact paper's
trichotomy subsection; `Screen/A5OPH.lean` (Part VI header) records the
lemma as the declared input this file discharges up to the torus step.
No physical content is at stake in this receipt lane. -/

/-- The doubled `A5` character table.  Row index: `0 ↦ 1`, `1 ↦ 3`,
`2 ↦ 3'`, `3 ↦ 4`, `4 ↦ 5`.  Column index: `0` identity, `1` double
transpositions, `2` three-cycles, `3` and `4` the two five-cycle classes.
Declared input of this module. -/
def chi2 : Fin 5 → Fin 5 → ℤ√5 := fun i c =>
  match i, c with
  | 0, 0 => ⟨2, 0⟩
  | 0, 1 => ⟨2, 0⟩
  | 0, 2 => ⟨2, 0⟩
  | 0, 3 => ⟨2, 0⟩
  | 0, 4 => ⟨2, 0⟩
  | 1, 0 => ⟨6, 0⟩
  | 1, 1 => ⟨-2, 0⟩
  | 1, 2 => ⟨0, 0⟩
  | 1, 3 => ⟨1, 1⟩
  | 1, 4 => ⟨1, -1⟩
  | 2, 0 => ⟨6, 0⟩
  | 2, 1 => ⟨-2, 0⟩
  | 2, 2 => ⟨0, 0⟩
  | 2, 3 => ⟨1, -1⟩
  | 2, 4 => ⟨1, 1⟩
  | 3, 0 => ⟨8, 0⟩
  | 3, 1 => ⟨0, 0⟩
  | 3, 2 => ⟨2, 0⟩
  | 3, 3 => ⟨-2, 0⟩
  | 3, 4 => ⟨-2, 0⟩
  | 4, 0 => ⟨10, 0⟩
  | 4, 1 => ⟨2, 0⟩
  | 4, 2 => ⟨-2, 0⟩
  | 4, 3 => ⟨0, 0⟩
  | 4, 4 => ⟨0, 0⟩

/-- The irreducible dimensions `1, 3, 3, 4, 5` in row order. -/
def dims : Fin 5 → ℕ := ![1, 3, 3, 4, 5]

/-- Doubled combined character of a module with multiplicity vector `m`. -/
def combined (m : Fin 5 → ℕ) (c : Fin 5) : ℤ√5 :=
  (m 0 : ℤ√5) * chi2 0 c + (m 1 : ℤ√5) * chi2 1 c + (m 2 : ℤ√5) * chi2 2 c
    + (m 3 : ℤ√5) * chi2 3 c + (m 4 : ℤ√5) * chi2 4 c

/-- Galois stability of a class function: the star involution (the
nontrivial element of `Gal(ℚ(√5)/ℚ)` restricted to `ℤ√5`) fixes every
value. -/
def GaloisStable (f : Fin 5 → ℤ√5) : Prop := ∀ c, star (f c) = f c

/-- The conjugation carries the doubled character of `3` to that of `3'`. -/
theorem star_chi2_one : ∀ c, star (chi2 1 c) = chi2 2 c := by decide

/-- The conjugation carries the doubled character of `3'` to that of `3`. -/
theorem star_chi2_two : ∀ c, star (chi2 2 c) = chi2 1 c := by decide

/-- The rows `1`, `4`, `5` are fixed by the conjugation. -/
theorem star_chi2_rational : ∀ c, star (chi2 0 c) = chi2 0 c
    ∧ star (chi2 3 c) = chi2 3 c ∧ star (chi2 4 c) = chi2 4 c := by decide

/-- Galois-stability half of the rationality lemma: a module whose combined
character is fixed by the conjugation carries the two three-dimensional
irreducibles with equal multiplicity.  Evaluate the imaginary part at the
first five-cycle class: it equals `m 1 - m 2`, and stability forces it to
vanish. -/
theorem multiplicities_equal_of_galoisStable (m : Fin 5 → ℕ)
    (h : GaloisStable (combined m)) : m 1 = m 2 := by
  have him : (combined m 3).im = (m 1 : ℤ) - (m 2 : ℤ) := by
    simp [combined, chi2]
    ring
  have h3 := congrArg Zsqrtd.im (h 3)
  rw [Zsqrtd.im_star, him] at h3
  omega

/-- Dimension readout consistency: the doubled combined character at the
identity class is twice the module dimension. -/
theorem dim_eq_char_at_identity (m : Fin 5 → ℕ) :
    (combined m 0).re
      = 2 * ((m 0) * dims 0 + (m 1) * dims 1 + (m 2) * dims 2
        + (m 3) * dims 3 + (m 4) * dims 4 : ℕ) := by
  simp [combined, chi2, dims]
  ring

/-- The dimension readout with the row dimensions written out. -/
theorem dims_readout (m : Fin 5 → ℕ) :
    m 0 * dims 0 + m 1 * dims 1 + m 2 * dims 2 + m 3 * dims 3
      + m 4 * dims 4
      = m 0 + 3 * m 1 + 3 * m 2 + 4 * m 3 + 5 * m 4 := by
  simp [dims]
  ring

/-- Inside the twelve-port module `1 + 3 + 3' + 5` the multiplicities are
bounded by `1, 1, 1, 0, 1`; with the Galois constraint `m 1 = m 2` the
dimension takes one of seven values. -/
theorem centreDim_cases (m : Fin 5 → ℕ) (h0 : m 0 ≤ 1) (h1 : m 1 ≤ 1)
    (h2 : m 2 ≤ 1) (h3 : m 3 = 0) (h4 : m 4 ≤ 1) (hg : m 1 = m 2) :
    let d := m 0 + 3 * m 1 + 3 * m 2 + 4 * m 3 + 5 * m 4
    d = 0 ∨ d = 1 ∨ d = 5 ∨ d = 6 ∨ d = 7 ∨ d = 11 ∨ d = 12 := by
  intro d
  have e0 : m 0 = 0 ∨ m 0 = 1 := by omega
  have e2 : m 2 = 0 ∨ m 2 = 1 := by omega
  have e4 : m 4 = 0 ∨ m 4 = 1 := by omega
  rcases e0 with e0 | e0 <;> rcases e2 with e2 | e2 <;>
    rcases e4 with e4 | e4 <;>
    simp [d, e0, e2, e4, h3, hg]

/-- The centre dimension of a Galois-stable submodule of the twelve-port
module lies in the subset-sum list of the trichotomy enumeration. -/
theorem centreDim_mem_trichotomy_list (m : Fin 5 → ℕ)
    (hstable : GaloisStable (combined m)) (h0 : m 0 ≤ 1) (h1 : m 1 ≤ 1)
    (h2 : m 2 ≤ 1) (h3 : m 3 = 0) (h4 : m 4 ≤ 1) :
    m 0 * dims 0 + m 1 * dims 1 + m 2 * dims 2 + m 3 * dims 3
      + m 4 * dims 4 ∈ ({0, 1, 5, 6, 7, 11, 12} : Finset ℕ) := by
  have hg := multiplicities_equal_of_galoisStable m hstable
  have hcases := centreDim_cases m h0 h1 h2 h3 h4 hg
  rw [dims_readout]
  simp only [Finset.mem_insert, Finset.mem_singleton]
  exact hcases

/-- Dimension two is impossible for a Galois-stable submodule of the
twelve-port module: the case excluded by the trichotomy's centre argument. -/
theorem centreDim_ne_two (m : Fin 5 → ℕ)
    (hstable : GaloisStable (combined m)) (h0 : m 0 ≤ 1) (h1 : m 1 ≤ 1)
    (h2 : m 2 ≤ 1) (h3 : m 3 = 0) (h4 : m 4 ≤ 1) :
    m 0 * dims 0 + m 1 * dims 1 + m 2 * dims 2 + m 3 * dims 3
      + m 4 * dims 4 ≠ 2 := by
  have hg := multiplicities_equal_of_galoisStable m hstable
  have hcases := centreDim_cases m h0 h1 h2 h3 h4 hg
  rw [dims_readout]
  rcases hcases with h | h | h | h | h | h | h <;> omega

end OPH.A5CharacterField

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.A5CharacterField.star_chi2_one
#print axioms OPH.A5CharacterField.star_chi2_two
#print axioms OPH.A5CharacterField.star_chi2_rational
#print axioms OPH.A5CharacterField.multiplicities_equal_of_galoisStable
#print axioms OPH.A5CharacterField.dim_eq_char_at_identity
#print axioms OPH.A5CharacterField.centreDim_cases
#print axioms OPH.A5CharacterField.centreDim_mem_trichotomy_list
#print axioms OPH.A5CharacterField.centreDim_ne_two
