import Mathlib

namespace OPH.A5SixAxes

/-! # The six-axis A5 action and the dimension-six branch module skeleton

Issue #604.  The compact-Lie trichotomy excludes the centre case with
semisimple part of dimension six by a module argument: the two `su(2)`
ideals would have to carry `1 + 5`, and the five-dimensional summand is
irreducible, so no decomposition into two three-dimensional invariant
summands exists.  This file supplies the concrete six-point `A5` action and
the module-theoretic half of that argument.

CONTENT.

* `t`, `s`, `L60`: the icosahedral six-axis action realized as
  `PSL(2, F5)` on the projective line `{0, 1, 2, 3, 4, ∞}` (index `5` is
  the point at infinity), with `t : z ↦ z + 1` and `s : z ↦ -1/z`.  `L60`
  lists all sixty group elements as explicit permutations of `Fin 6`.
* Kernel-`decide` facts: `L60` has sixty distinct members, contains the
  identity and both generators, is closed under inverse and product, and
  acts 2-transitively (`two_transitive`).  The quadratic checks are
  kernel-decided on the raw value tables `rowF`/`rowI` against the
  precomputed index tables `mulT`/`invT`, then transported to the listed
  permutations through the pointwise bridge `el_apply`; `Equiv.Perm`
  equality and multiplication are never kernel-decided.
* `V5`: the sum-zero hyperplane of the permutation module `Fin 6 → ℚ`,
  invariant under the action (`V5_invariant`), of dimension five
  (`finrank_V5`).
* `no_three_dim_invariant`: GIVEN that `V5` is irreducible under the
  listed action, no invariant subspace of the permutation module has
  dimension three.  The proof is the dimension count
  `finrank (W ⊔ V5) + finrank (W ⊓ V5) = 3 + 5` with `finrank (W ⊔ V5) ≤ 6`,
  so `W ⊓ V5` is a nonzero invariant subspace of `V5`, hence all of `V5`,
  hence `finrank W ≥ 5`.  This excludes the two-summand decomposition
  demanded by the dimension-six branch, since a summand of a
  three-plus-three splitting is a three-dimensional invariant subspace.
* `V5_irreducible`: the irreducibility hypothesis is discharged.  The
  engine is the kernel-checked sharp fiber count `count_fibers` of the
  listed action (stabilizer size ten, sharp 2-transitive fiber size two);
  averaging any nonzero sum-zero vector over the stabilizer cosets of a
  coordinate where it does not vanish produces every difference vector
  `e k - e l` inside an invariant subspace, and these span `V5`.
* `no_three_dim_invariant_unconditional`, `no_three_plus_three_split`:
  the dimension-six branch exclusion with the hypothesis discharged.

BOUNDARY.  The conditional form `no_three_dim_invariant` is retained so the
module argument and the irreducibility input remain separately auditable.
No physical content is at stake in this receipt lane. -/

/-- Build a permutation of `Fin 6` from its value vector and inverse
vector; the two inverse laws are discharged by `decide` at each literal. -/
def perm (f g : Fin 6 → Fin 6)
    (h₁ : Function.LeftInverse g f := by decide)
    (h₂ : Function.RightInverse g f := by decide) : Equiv.Perm (Fin 6) :=
  ⟨f, g, h₁, h₂⟩

/-- Generator `t : z ↦ z + 1` on `P¹(F₅)`, the five-cycle fixing `∞`. -/
def t : Equiv.Perm (Fin 6) := perm ![1, 2, 3, 4, 0, 5] ![4, 0, 1, 2, 3, 5]

/-- Generator `s : z ↦ -1/z` on `P¹(F₅)`, the involution `(0 ∞)(1 4)`. -/
def s : Equiv.Perm (Fin 6) := perm ![5, 4, 2, 3, 1, 0] ![5, 4, 2, 3, 1, 0]

/-- All sixty elements of the generated copy of `A5 ≅ PSL(2, F5)` acting on
the six axes. -/
def L60 : List (Equiv.Perm (Fin 6)) := [
    perm ![0, 1, 2, 3, 4, 5] ![0, 1, 2, 3, 4, 5],
    perm ![0, 1, 4, 5, 2, 3] ![0, 1, 4, 5, 2, 3],
    perm ![0, 2, 1, 3, 5, 4] ![0, 2, 1, 3, 5, 4],
    perm ![0, 2, 5, 4, 1, 3] ![0, 4, 1, 5, 3, 2],
    perm ![0, 3, 4, 2, 5, 1] ![0, 5, 3, 1, 2, 4],
    perm ![0, 3, 5, 1, 4, 2] ![0, 3, 5, 1, 4, 2],
    perm ![0, 4, 1, 5, 3, 2] ![0, 2, 5, 4, 1, 3],
    perm ![0, 4, 3, 2, 1, 5] ![0, 4, 3, 2, 1, 5],
    perm ![0, 5, 2, 4, 3, 1] ![0, 5, 2, 4, 3, 1],
    perm ![0, 5, 3, 1, 2, 4] ![0, 3, 4, 2, 5, 1],
    perm ![1, 0, 2, 5, 4, 3] ![1, 0, 2, 5, 4, 3],
    perm ![1, 0, 4, 3, 2, 5] ![1, 0, 4, 3, 2, 5],
    perm ![1, 2, 0, 5, 3, 4] ![2, 0, 1, 4, 5, 3],
    perm ![1, 2, 3, 4, 0, 5] ![4, 0, 1, 2, 3, 5],
    perm ![1, 3, 2, 4, 5, 0] ![5, 0, 2, 1, 3, 4],
    perm ![1, 3, 5, 0, 2, 4] ![3, 0, 4, 1, 5, 2],
    perm ![1, 4, 0, 3, 5, 2] ![2, 0, 5, 3, 1, 4],
    perm ![1, 4, 5, 2, 0, 3] ![4, 0, 3, 5, 1, 2],
    perm ![1, 5, 3, 0, 4, 2] ![3, 0, 5, 2, 4, 1],
    perm ![1, 5, 4, 2, 3, 0] ![5, 0, 3, 4, 2, 1],
    perm ![2, 0, 1, 4, 5, 3] ![1, 2, 0, 5, 3, 4],
    perm ![2, 0, 5, 3, 1, 4] ![1, 4, 0, 3, 5, 2],
    perm ![2, 1, 0, 4, 3, 5] ![2, 1, 0, 4, 3, 5],
    perm ![2, 1, 3, 5, 0, 4] ![4, 1, 0, 2, 5, 3],
    perm ![2, 3, 1, 5, 4, 0] ![5, 2, 0, 1, 4, 3],
    perm ![2, 3, 4, 0, 1, 5] ![3, 4, 0, 1, 2, 5],
    perm ![2, 4, 3, 0, 5, 1] ![3, 5, 0, 2, 1, 4],
    perm ![2, 4, 5, 1, 3, 0] ![5, 3, 0, 4, 1, 2],
    perm ![2, 5, 0, 3, 4, 1] ![2, 5, 0, 3, 4, 1],
    perm ![2, 5, 4, 1, 0, 3] ![4, 3, 0, 5, 2, 1],
    perm ![3, 0, 4, 1, 5, 2] ![1, 3, 5, 0, 2, 4],
    perm ![3, 0, 5, 2, 4, 1] ![1, 5, 3, 0, 4, 2],
    perm ![3, 1, 2, 0, 5, 4] ![3, 1, 2, 0, 5, 4],
    perm ![3, 1, 5, 4, 2, 0] ![5, 1, 4, 0, 3, 2],
    perm ![3, 2, 1, 0, 4, 5] ![3, 2, 1, 0, 4, 5],
    perm ![3, 2, 4, 5, 1, 0] ![5, 4, 1, 0, 2, 3],
    perm ![3, 4, 0, 1, 2, 5] ![2, 3, 4, 0, 1, 5],
    perm ![3, 4, 2, 5, 0, 1] ![4, 5, 2, 0, 1, 3],
    perm ![3, 5, 0, 2, 1, 4] ![2, 4, 3, 0, 5, 1],
    perm ![3, 5, 1, 4, 0, 2] ![4, 2, 5, 0, 3, 1],
    perm ![4, 0, 1, 2, 3, 5] ![1, 2, 3, 4, 0, 5],
    perm ![4, 0, 3, 5, 1, 2] ![1, 4, 5, 2, 0, 3],
    perm ![4, 1, 0, 2, 5, 3] ![2, 1, 3, 5, 0, 4],
    perm ![4, 1, 5, 3, 0, 2] ![4, 1, 5, 3, 0, 2],
    perm ![4, 2, 3, 1, 5, 0] ![5, 3, 1, 2, 0, 4],
    perm ![4, 2, 5, 0, 3, 1] ![3, 5, 1, 4, 0, 2],
    perm ![4, 3, 0, 5, 2, 1] ![2, 5, 4, 1, 0, 3],
    perm ![4, 3, 2, 1, 0, 5] ![4, 3, 2, 1, 0, 5],
    perm ![4, 5, 1, 3, 2, 0] ![5, 2, 4, 3, 0, 1],
    perm ![4, 5, 2, 0, 1, 3] ![3, 4, 2, 5, 0, 1],
    perm ![5, 0, 2, 1, 3, 4] ![1, 3, 2, 4, 5, 0],
    perm ![5, 0, 3, 4, 2, 1] ![1, 5, 4, 2, 3, 0],
    perm ![5, 1, 3, 2, 4, 0] ![5, 1, 3, 2, 4, 0],
    perm ![5, 1, 4, 0, 3, 2] ![3, 1, 5, 4, 2, 0],
    perm ![5, 2, 0, 1, 4, 3] ![2, 3, 1, 5, 4, 0],
    perm ![5, 2, 4, 3, 0, 1] ![4, 5, 1, 3, 2, 0],
    perm ![5, 3, 0, 4, 1, 2] ![2, 4, 5, 1, 3, 0],
    perm ![5, 3, 1, 2, 0, 4] ![4, 2, 3, 1, 5, 0],
    perm ![5, 4, 1, 0, 2, 3] ![3, 2, 4, 5, 1, 0],
    perm ![5, 4, 2, 3, 1, 0] ![5, 4, 2, 3, 1, 0],  ]

/-! ## Raw value tables

Kernel `decide` over `Equiv.Perm (Fin 6)` is intractable at quadratic
scale: deciding equality or membership forces each product `g * h`
through the group-instance tower and the `Equiv` inverse-law proofs, and
`mul_closed` alone visits `3600` products against `60` candidates.  The
combinatorial content is therefore kernel-checked on raw `Fin 6` value
vectors with precomputed index tables, and transported to the listed
permutations through the pointwise bridge `el_apply`.  The tables were
computed offline from `L60` and are verified here by `decide`:
`rowF_mul` and `rowF_inv` are exhaustive checks, not trusted inputs. -/

def rowF : Fin 60 → Fin 6 → Fin 6 :=
  ![![0, 1, 2, 3, 4, 5],
    ![0, 1, 4, 5, 2, 3],
    ![0, 2, 1, 3, 5, 4],
    ![0, 2, 5, 4, 1, 3],
    ![0, 3, 4, 2, 5, 1],
    ![0, 3, 5, 1, 4, 2],
    ![0, 4, 1, 5, 3, 2],
    ![0, 4, 3, 2, 1, 5],
    ![0, 5, 2, 4, 3, 1],
    ![0, 5, 3, 1, 2, 4],
    ![1, 0, 2, 5, 4, 3],
    ![1, 0, 4, 3, 2, 5],
    ![1, 2, 0, 5, 3, 4],
    ![1, 2, 3, 4, 0, 5],
    ![1, 3, 2, 4, 5, 0],
    ![1, 3, 5, 0, 2, 4],
    ![1, 4, 0, 3, 5, 2],
    ![1, 4, 5, 2, 0, 3],
    ![1, 5, 3, 0, 4, 2],
    ![1, 5, 4, 2, 3, 0],
    ![2, 0, 1, 4, 5, 3],
    ![2, 0, 5, 3, 1, 4],
    ![2, 1, 0, 4, 3, 5],
    ![2, 1, 3, 5, 0, 4],
    ![2, 3, 1, 5, 4, 0],
    ![2, 3, 4, 0, 1, 5],
    ![2, 4, 3, 0, 5, 1],
    ![2, 4, 5, 1, 3, 0],
    ![2, 5, 0, 3, 4, 1],
    ![2, 5, 4, 1, 0, 3],
    ![3, 0, 4, 1, 5, 2],
    ![3, 0, 5, 2, 4, 1],
    ![3, 1, 2, 0, 5, 4],
    ![3, 1, 5, 4, 2, 0],
    ![3, 2, 1, 0, 4, 5],
    ![3, 2, 4, 5, 1, 0],
    ![3, 4, 0, 1, 2, 5],
    ![3, 4, 2, 5, 0, 1],
    ![3, 5, 0, 2, 1, 4],
    ![3, 5, 1, 4, 0, 2],
    ![4, 0, 1, 2, 3, 5],
    ![4, 0, 3, 5, 1, 2],
    ![4, 1, 0, 2, 5, 3],
    ![4, 1, 5, 3, 0, 2],
    ![4, 2, 3, 1, 5, 0],
    ![4, 2, 5, 0, 3, 1],
    ![4, 3, 0, 5, 2, 1],
    ![4, 3, 2, 1, 0, 5],
    ![4, 5, 1, 3, 2, 0],
    ![4, 5, 2, 0, 1, 3],
    ![5, 0, 2, 1, 3, 4],
    ![5, 0, 3, 4, 2, 1],
    ![5, 1, 3, 2, 4, 0],
    ![5, 1, 4, 0, 3, 2],
    ![5, 2, 0, 1, 4, 3],
    ![5, 2, 4, 3, 0, 1],
    ![5, 3, 0, 4, 1, 2],
    ![5, 3, 1, 2, 0, 4],
    ![5, 4, 1, 0, 2, 3],
    ![5, 4, 2, 3, 1, 0]]

def rowI : Fin 60 → Fin 6 → Fin 6 :=
  ![![0, 1, 2, 3, 4, 5],
    ![0, 1, 4, 5, 2, 3],
    ![0, 2, 1, 3, 5, 4],
    ![0, 4, 1, 5, 3, 2],
    ![0, 5, 3, 1, 2, 4],
    ![0, 3, 5, 1, 4, 2],
    ![0, 2, 5, 4, 1, 3],
    ![0, 4, 3, 2, 1, 5],
    ![0, 5, 2, 4, 3, 1],
    ![0, 3, 4, 2, 5, 1],
    ![1, 0, 2, 5, 4, 3],
    ![1, 0, 4, 3, 2, 5],
    ![2, 0, 1, 4, 5, 3],
    ![4, 0, 1, 2, 3, 5],
    ![5, 0, 2, 1, 3, 4],
    ![3, 0, 4, 1, 5, 2],
    ![2, 0, 5, 3, 1, 4],
    ![4, 0, 3, 5, 1, 2],
    ![3, 0, 5, 2, 4, 1],
    ![5, 0, 3, 4, 2, 1],
    ![1, 2, 0, 5, 3, 4],
    ![1, 4, 0, 3, 5, 2],
    ![2, 1, 0, 4, 3, 5],
    ![4, 1, 0, 2, 5, 3],
    ![5, 2, 0, 1, 4, 3],
    ![3, 4, 0, 1, 2, 5],
    ![3, 5, 0, 2, 1, 4],
    ![5, 3, 0, 4, 1, 2],
    ![2, 5, 0, 3, 4, 1],
    ![4, 3, 0, 5, 2, 1],
    ![1, 3, 5, 0, 2, 4],
    ![1, 5, 3, 0, 4, 2],
    ![3, 1, 2, 0, 5, 4],
    ![5, 1, 4, 0, 3, 2],
    ![3, 2, 1, 0, 4, 5],
    ![5, 4, 1, 0, 2, 3],
    ![2, 3, 4, 0, 1, 5],
    ![4, 5, 2, 0, 1, 3],
    ![2, 4, 3, 0, 5, 1],
    ![4, 2, 5, 0, 3, 1],
    ![1, 2, 3, 4, 0, 5],
    ![1, 4, 5, 2, 0, 3],
    ![2, 1, 3, 5, 0, 4],
    ![4, 1, 5, 3, 0, 2],
    ![5, 3, 1, 2, 0, 4],
    ![3, 5, 1, 4, 0, 2],
    ![2, 5, 4, 1, 0, 3],
    ![4, 3, 2, 1, 0, 5],
    ![5, 2, 4, 3, 0, 1],
    ![3, 4, 2, 5, 0, 1],
    ![1, 3, 2, 4, 5, 0],
    ![1, 5, 4, 2, 3, 0],
    ![5, 1, 3, 2, 4, 0],
    ![3, 1, 5, 4, 2, 0],
    ![2, 3, 1, 5, 4, 0],
    ![4, 5, 1, 3, 2, 0],
    ![2, 4, 5, 1, 3, 0],
    ![4, 2, 3, 1, 5, 0],
    ![3, 2, 4, 5, 1, 0],
    ![5, 4, 2, 3, 1, 0]]

def invT : Fin 60 → Fin 60 :=
  ![0, 1, 2, 6, 9, 5, 3, 7, 8, 4, 10, 11, 20, 40, 50, 30, 21, 41, 31, 51, 12, 16, 22, 42, 54, 36, 38, 56, 28, 46, 15, 18, 32, 53, 34, 58, 25, 49, 26, 45, 13, 17, 23, 43, 57, 39, 29, 47, 55, 37, 14, 19, 52, 33, 24, 48, 27, 44, 35, 59]

def mulT : Fin 60 → Fin 60 → Fin 60 :=
  ![![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
    ![1, 0, 6, 7, 8, 9, 2, 3, 4, 5, 11, 10, 16, 17, 19, 18, 12, 13, 15, 14, 40, 41, 42, 43, 48, 49, 45, 44, 46, 47, 50, 51, 53, 52, 58, 59, 54, 55, 56, 57, 20, 21, 22, 23, 27, 26, 28, 29, 24, 25, 30, 31, 33, 32, 36, 37, 38, 39, 34, 35],
    ![2, 3, 0, 1, 5, 4, 8, 9, 6, 7, 20, 21, 22, 23, 24, 25, 28, 29, 26, 27, 10, 11, 12, 13, 14, 15, 18, 19, 16, 17, 31, 30, 34, 35, 32, 33, 38, 39, 36, 37, 50, 51, 54, 55, 52, 53, 56, 57, 59, 58, 40, 41, 44, 45, 42, 43, 46, 47, 49, 48],
    ![3, 2, 8, 9, 6, 7, 0, 1, 5, 4, 21, 20, 28, 29, 27, 26, 22, 23, 25, 24, 50, 51, 54, 55, 59, 58, 53, 52, 56, 57, 40, 41, 45, 44, 49, 48, 42, 43, 46, 47, 10, 11, 12, 13, 19, 18, 16, 17, 14, 15, 31, 30, 35, 34, 38, 39, 36, 37, 32, 33],
    ![4, 5, 7, 6, 3, 2, 9, 8, 1, 0, 30, 31, 36, 37, 35, 34, 38, 39, 32, 33, 41, 40, 46, 47, 44, 45, 49, 48, 42, 43, 21, 20, 25, 24, 26, 27, 28, 29, 22, 23, 51, 50, 56, 57, 59, 58, 54, 55, 52, 53, 11, 10, 14, 15, 16, 17, 12, 13, 18, 19],
    ![5, 4, 9, 8, 1, 0, 7, 6, 3, 2, 31, 30, 38, 39, 33, 32, 36, 37, 34, 35, 51, 50, 56, 57, 52, 53, 58, 59, 54, 55, 11, 10, 15, 14, 18, 19, 16, 17, 12, 13, 41, 40, 46, 47, 48, 49, 42, 43, 44, 45, 21, 20, 24, 25, 28, 29, 22, 23, 26, 27],
    ![6, 7, 1, 0, 9, 8, 4, 5, 2, 3, 40, 41, 42, 43, 48, 49, 46, 47, 45, 44, 11, 10, 16, 17, 19, 18, 15, 14, 12, 13, 51, 50, 58, 59, 53, 52, 56, 57, 54, 55, 30, 31, 36, 37, 33, 32, 38, 39, 35, 34, 20, 21, 27, 26, 22, 23, 28, 29, 25, 24],
    ![7, 6, 4, 5, 2, 3, 1, 0, 9, 8, 41, 40, 46, 47, 44, 45, 42, 43, 49, 48, 30, 31, 36, 37, 35, 34, 32, 33, 38, 39, 20, 21, 26, 27, 25, 24, 22, 23, 28, 29, 11, 10, 16, 17, 14, 15, 12, 13, 19, 18, 51, 50, 59, 58, 56, 57, 54, 55, 53, 52],
    ![8, 9, 3, 2, 7, 6, 5, 4, 0, 1, 50, 51, 54, 55, 59, 58, 56, 57, 53, 52, 21, 20, 28, 29, 27, 26, 25, 24, 22, 23, 41, 40, 49, 48, 45, 44, 46, 47, 42, 43, 31, 30, 38, 39, 35, 34, 36, 37, 33, 32, 10, 11, 19, 18, 12, 13, 16, 17, 15, 14],
    ![9, 8, 5, 4, 0, 1, 3, 2, 7, 6, 51, 50, 56, 57, 52, 53, 54, 55, 58, 59, 31, 30, 38, 39, 33, 32, 34, 35, 36, 37, 10, 11, 18, 19, 15, 14, 12, 13, 16, 17, 21, 20, 28, 29, 24, 25, 22, 23, 27, 26, 41, 40, 48, 49, 46, 47, 42, 43, 45, 44],
    ![10, 11, 12, 13, 19, 18, 16, 17, 14, 15, 0, 1, 2, 3, 8, 9, 6, 7, 5, 4, 22, 23, 20, 21, 28, 29, 27, 26, 24, 25, 53, 52, 50, 51, 54, 55, 58, 59, 57, 56, 42, 43, 40, 41, 45, 44, 48, 49, 46, 47, 32, 33, 31, 30, 34, 35, 39, 38, 36, 37],
    ![11, 10, 16, 17, 14, 15, 12, 13, 19, 18, 1, 0, 6, 7, 4, 5, 2, 3, 9, 8, 42, 43, 40, 41, 46, 47, 44, 45, 48, 49, 32, 33, 30, 31, 36, 37, 34, 35, 39, 38, 22, 23, 20, 21, 26, 27, 24, 25, 28, 29, 53, 52, 51, 50, 58, 59, 57, 56, 54, 55],
    ![12, 13, 10, 11, 18, 19, 14, 15, 16, 17, 22, 23, 20, 21, 28, 29, 24, 25, 27, 26, 0, 1, 2, 3, 8, 9, 5, 4, 6, 7, 52, 53, 54, 55, 50, 51, 57, 56, 58, 59, 32, 33, 34, 35, 31, 30, 39, 38, 37, 36, 42, 43, 45, 44, 40, 41, 48, 49, 47, 46],
    ![13, 12, 14, 15, 16, 17, 10, 11, 18, 19, 23, 22, 24, 25, 26, 27, 20, 21, 29, 28, 32, 33, 34, 35, 37, 36, 30, 31, 39, 38, 42, 43, 44, 45, 47, 46, 40, 41, 48, 49, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 52, 53, 55, 54, 57, 56, 58, 59, 50, 51],
    ![14, 15, 13, 12, 17, 16, 18, 19, 10, 11, 32, 33, 34, 35, 37, 36, 39, 38, 30, 31, 23, 22, 24, 25, 26, 27, 29, 28, 20, 21, 43, 42, 47, 46, 44, 45, 48, 49, 40, 41, 52, 53, 57, 56, 55, 54, 58, 59, 51, 50, 0, 1, 4, 5, 2, 3, 6, 7, 9, 8],
    ![15, 14, 18, 19, 10, 11, 13, 12, 17, 16, 33, 32, 39, 38, 31, 30, 34, 35, 36, 37, 52, 53, 57, 56, 51, 50, 54, 55, 58, 59, 0, 1, 5, 4, 9, 8, 2, 3, 6, 7, 23, 22, 24, 25, 28, 29, 20, 21, 26, 27, 43, 42, 46, 47, 48, 49, 40, 41, 44, 45],
    ![16, 17, 11, 10, 15, 14, 19, 18, 12, 13, 42, 43, 40, 41, 46, 47, 48, 49, 44, 45, 1, 0, 6, 7, 4, 5, 9, 8, 2, 3, 33, 32, 36, 37, 30, 31, 39, 38, 34, 35, 53, 52, 58, 59, 51, 50, 57, 56, 55, 54, 22, 23, 26, 27, 20, 21, 24, 25, 29, 28],
    ![17, 16, 19, 18, 12, 13, 11, 10, 15, 14, 43, 42, 48, 49, 45, 44, 40, 41, 47, 46, 53, 52, 58, 59, 55, 54, 50, 51, 57, 56, 22, 23, 27, 26, 29, 28, 20, 21, 24, 25, 1, 0, 6, 7, 8, 9, 2, 3, 4, 5, 33, 32, 37, 36, 39, 38, 34, 35, 30, 31],
    ![18, 19, 15, 14, 11, 10, 17, 16, 13, 12, 52, 53, 57, 56, 51, 50, 58, 59, 54, 55, 33, 32, 39, 38, 31, 30, 36, 37, 34, 35, 1, 0, 9, 8, 5, 4, 6, 7, 2, 3, 43, 42, 48, 49, 46, 47, 40, 41, 45, 44, 23, 22, 28, 29, 24, 25, 20, 21, 27, 26],
    ![19, 18, 17, 16, 13, 12, 15, 14, 11, 10, 53, 52, 58, 59, 55, 54, 57, 56, 50, 51, 43, 42, 48, 49, 45, 44, 47, 46, 40, 41, 23, 22, 29, 28, 27, 26, 24, 25, 20, 21, 33, 32, 39, 38, 37, 36, 34, 35, 31, 30, 1, 0, 8, 9, 6, 7, 2, 3, 5, 4],
    ![20, 21, 22, 23, 27, 26, 28, 29, 24, 25, 2, 3, 0, 1, 6, 7, 8, 9, 4, 5, 12, 13, 10, 11, 16, 17, 19, 18, 14, 15, 45, 44, 40, 41, 42, 43, 49, 48, 47, 46, 54, 55, 50, 51, 53, 52, 59, 58, 56, 57, 34, 35, 30, 31, 32, 33, 37, 36, 38, 39],
    ![21, 20, 28, 29, 24, 25, 22, 23, 27, 26, 3, 2, 8, 9, 5, 4, 0, 1, 7, 6, 54, 55, 50, 51, 56, 57, 52, 53, 59, 58, 34, 35, 31, 30, 38, 39, 32, 33, 37, 36, 12, 13, 10, 11, 18, 19, 14, 15, 16, 17, 45, 44, 41, 40, 49, 48, 47, 46, 42, 43],
    ![22, 23, 20, 21, 26, 27, 24, 25, 28, 29, 12, 13, 10, 11, 16, 17, 14, 15, 19, 18, 2, 3, 0, 1, 6, 7, 4, 5, 8, 9, 44, 45, 42, 43, 40, 41, 47, 46, 49, 48, 34, 35, 32, 33, 30, 31, 37, 36, 39, 38, 54, 55, 53, 52, 50, 51, 59, 58, 57, 56],
    ![23, 22, 24, 25, 28, 29, 20, 21, 26, 27, 13, 12, 14, 15, 18, 19, 10, 11, 17, 16, 34, 35, 32, 33, 39, 38, 31, 30, 37, 36, 54, 55, 52, 53, 57, 56, 50, 51, 59, 58, 2, 3, 0, 1, 5, 4, 8, 9, 6, 7, 44, 45, 43, 42, 47, 46, 49, 48, 40, 41],
    ![24, 25, 23, 22, 29, 28, 26, 27, 20, 21, 34, 35, 32, 33, 39, 38, 37, 36, 31, 30, 13, 12, 14, 15, 18, 19, 17, 16, 10, 11, 55, 54, 57, 56, 52, 53, 59, 58, 50, 51, 44, 45, 47, 46, 43, 42, 49, 48, 41, 40, 2, 3, 5, 4, 0, 1, 8, 9, 7, 6],
    ![25, 24, 26, 27, 20, 21, 23, 22, 29, 28, 35, 34, 37, 36, 30, 31, 32, 33, 38, 39, 44, 45, 47, 46, 41, 40, 42, 43, 49, 48, 2, 3, 4, 5, 7, 6, 0, 1, 8, 9, 13, 12, 14, 15, 16, 17, 10, 11, 18, 19, 55, 54, 56, 57, 59, 58, 50, 51, 52, 53],
    ![26, 27, 25, 24, 21, 20, 29, 28, 23, 22, 44, 45, 47, 46, 41, 40, 49, 48, 42, 43, 35, 34, 37, 36, 30, 31, 38, 39, 32, 33, 3, 2, 7, 6, 4, 5, 8, 9, 0, 1, 55, 54, 59, 58, 56, 57, 50, 51, 53, 52, 13, 12, 16, 17, 14, 15, 10, 11, 19, 18],
    ![27, 26, 29, 28, 23, 22, 25, 24, 21, 20, 45, 44, 49, 48, 43, 42, 47, 46, 40, 41, 55, 54, 59, 58, 53, 52, 57, 56, 50, 51, 13, 12, 17, 16, 19, 18, 14, 15, 10, 11, 35, 34, 37, 36, 39, 38, 32, 33, 30, 31, 3, 2, 6, 7, 8, 9, 0, 1, 4, 5],
    ![28, 29, 21, 20, 25, 24, 27, 26, 22, 23, 54, 55, 50, 51, 56, 57, 59, 58, 52, 53, 3, 2, 8, 9, 5, 4, 7, 6, 0, 1, 35, 34, 38, 39, 31, 30, 37, 36, 32, 33, 45, 44, 49, 48, 41, 40, 47, 46, 43, 42, 12, 13, 18, 19, 10, 11, 14, 15, 17, 16],
    ![29, 28, 27, 26, 22, 23, 21, 20, 25, 24, 55, 54, 59, 58, 53, 52, 50, 51, 57, 56, 45, 44, 49, 48, 43, 42, 40, 41, 47, 46, 12, 13, 19, 18, 17, 16, 10, 11, 14, 15, 3, 2, 8, 9, 6, 7, 0, 1, 5, 4, 35, 34, 39, 38, 37, 36, 32, 33, 31, 30],
    ![30, 31, 36, 37, 33, 32, 38, 39, 35, 34, 4, 5, 7, 6, 1, 0, 9, 8, 2, 3, 46, 47, 41, 40, 42, 43, 48, 49, 44, 45, 15, 14, 11, 10, 16, 17, 18, 19, 13, 12, 56, 57, 51, 50, 58, 59, 52, 53, 54, 55, 25, 24, 20, 21, 26, 27, 23, 22, 28, 29],
    ![31, 30, 38, 39, 35, 34, 36, 37, 33, 32, 5, 4, 9, 8, 3, 2, 7, 6, 0, 1, 56, 57, 51, 50, 54, 55, 59, 58, 52, 53, 25, 24, 21, 20, 28, 29, 26, 27, 23, 22, 46, 47, 41, 40, 49, 48, 44, 45, 42, 43, 15, 14, 10, 11, 18, 19, 13, 12, 16, 17],
    ![32, 33, 34, 35, 31, 30, 39, 38, 37, 36, 14, 15, 13, 12, 10, 11, 18, 19, 16, 17, 24, 25, 23, 22, 20, 21, 28, 29, 26, 27, 5, 4, 0, 1, 2, 3, 9, 8, 7, 6, 57, 56, 52, 53, 54, 55, 51, 50, 58, 59, 47, 46, 42, 43, 44, 45, 41, 40, 48, 49],
    ![33, 32, 39, 38, 37, 36, 34, 35, 31, 30, 15, 14, 18, 19, 17, 16, 13, 12, 11, 10, 57, 56, 52, 53, 58, 59, 55, 54, 51, 50, 47, 46, 43, 42, 48, 49, 44, 45, 41, 40, 24, 25, 23, 22, 29, 28, 26, 27, 20, 21, 5, 4, 1, 0, 9, 8, 7, 6, 2, 3],
    ![34, 35, 32, 33, 30, 31, 37, 36, 39, 38, 24, 25, 23, 22, 20, 21, 26, 27, 28, 29, 14, 15, 13, 12, 10, 11, 16, 17, 18, 19, 4, 5, 2, 3, 0, 1, 7, 6, 9, 8, 47, 46, 44, 45, 42, 43, 41, 40, 49, 48, 57, 56, 54, 55, 52, 53, 51, 50, 59, 58],
    ![35, 34, 37, 36, 39, 38, 32, 33, 30, 31, 25, 24, 26, 27, 29, 28, 23, 22, 21, 20, 47, 46, 44, 45, 49, 48, 43, 42, 41, 40, 57, 56, 55, 54, 59, 58, 52, 53, 51, 50, 14, 15, 13, 12, 17, 16, 18, 19, 10, 11, 4, 5, 3, 2, 7, 6, 9, 8, 0, 1],
    ![36, 37, 30, 31, 32, 33, 35, 34, 38, 39, 46, 47, 41, 40, 42, 43, 44, 45, 48, 49, 4, 5, 7, 6, 1, 0, 2, 3, 9, 8, 14, 15, 16, 17, 11, 10, 13, 12, 18, 19, 25, 24, 26, 27, 20, 21, 23, 22, 29, 28, 56, 57, 58, 59, 51, 50, 52, 53, 55, 54],
    ![37, 36, 35, 34, 38, 39, 30, 31, 32, 33, 47, 46, 44, 45, 49, 48, 41, 40, 43, 42, 25, 24, 26, 27, 29, 28, 21, 20, 23, 22, 56, 57, 59, 58, 55, 54, 51, 50, 52, 53, 4, 5, 7, 6, 3, 2, 9, 8, 1, 0, 14, 15, 17, 16, 13, 12, 18, 19, 11, 10],
    ![38, 39, 31, 30, 34, 35, 33, 32, 36, 37, 56, 57, 51, 50, 54, 55, 52, 53, 59, 58, 5, 4, 9, 8, 3, 2, 0, 1, 7, 6, 24, 25, 28, 29, 21, 20, 23, 22, 26, 27, 15, 14, 18, 19, 10, 11, 13, 12, 17, 16, 46, 47, 49, 48, 41, 40, 44, 45, 43, 42],
    ![39, 38, 33, 32, 36, 37, 31, 30, 34, 35, 57, 56, 52, 53, 58, 59, 51, 50, 55, 54, 15, 14, 18, 19, 17, 16, 11, 10, 13, 12, 46, 47, 48, 49, 43, 42, 41, 40, 44, 45, 5, 4, 9, 8, 1, 0, 7, 6, 3, 2, 24, 25, 29, 28, 23, 22, 26, 27, 21, 20],
    ![40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 6, 7, 1, 0, 2, 3, 4, 5, 8, 9, 16, 17, 11, 10, 12, 13, 14, 15, 19, 18, 26, 27, 20, 21, 22, 23, 25, 24, 29, 28, 36, 37, 30, 31, 32, 33, 35, 34, 38, 39, 58, 59, 50, 51, 53, 52, 55, 54, 56, 57],
    ![41, 40, 46, 47, 48, 49, 42, 43, 44, 45, 7, 6, 4, 5, 9, 8, 1, 0, 3, 2, 36, 37, 30, 31, 38, 39, 33, 32, 35, 34, 58, 59, 51, 50, 56, 57, 53, 52, 55, 54, 16, 17, 11, 10, 15, 14, 19, 18, 12, 13, 26, 27, 21, 20, 25, 24, 29, 28, 22, 23],
    ![42, 43, 40, 41, 45, 44, 48, 49, 46, 47, 16, 17, 11, 10, 12, 13, 19, 18, 14, 15, 6, 7, 1, 0, 2, 3, 8, 9, 4, 5, 27, 26, 22, 23, 20, 21, 29, 28, 25, 24, 58, 59, 53, 52, 50, 51, 55, 54, 57, 56, 36, 37, 32, 33, 30, 31, 35, 34, 39, 38],
    ![43, 42, 48, 49, 46, 47, 40, 41, 45, 44, 17, 16, 19, 18, 15, 14, 11, 10, 13, 12, 58, 59, 53, 52, 57, 56, 51, 50, 55, 54, 36, 37, 33, 32, 39, 38, 30, 31, 35, 34, 6, 7, 1, 0, 9, 8, 4, 5, 2, 3, 27, 26, 23, 22, 29, 28, 25, 24, 20, 21],
    ![44, 45, 47, 46, 43, 42, 49, 48, 41, 40, 26, 27, 25, 24, 23, 22, 29, 28, 20, 21, 37, 36, 35, 34, 32, 33, 39, 38, 30, 31, 17, 16, 13, 12, 14, 15, 19, 18, 11, 10, 59, 58, 55, 54, 57, 56, 53, 52, 50, 51, 7, 6, 2, 3, 4, 5, 1, 0, 8, 9],
    ![45, 44, 49, 48, 41, 40, 47, 46, 43, 42, 27, 26, 29, 28, 21, 20, 25, 24, 22, 23, 59, 58, 55, 54, 50, 51, 56, 57, 53, 52, 7, 6, 3, 2, 8, 9, 4, 5, 1, 0, 37, 36, 35, 34, 38, 39, 30, 31, 32, 33, 17, 16, 12, 13, 19, 18, 11, 10, 14, 15],
    ![46, 47, 41, 40, 49, 48, 44, 45, 42, 43, 36, 37, 30, 31, 38, 39, 35, 34, 33, 32, 7, 6, 4, 5, 9, 8, 3, 2, 1, 0, 59, 58, 56, 57, 51, 50, 55, 54, 53, 52, 26, 27, 25, 24, 21, 20, 29, 28, 23, 22, 16, 17, 15, 14, 11, 10, 19, 18, 13, 12],
    ![47, 46, 44, 45, 42, 43, 41, 40, 49, 48, 37, 36, 35, 34, 32, 33, 30, 31, 39, 38, 26, 27, 25, 24, 23, 22, 20, 21, 29, 28, 16, 17, 14, 15, 13, 12, 11, 10, 19, 18, 7, 6, 4, 5, 2, 3, 1, 0, 9, 8, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50],
    ![48, 49, 43, 42, 47, 46, 45, 44, 40, 41, 58, 59, 53, 52, 57, 56, 55, 54, 51, 50, 17, 16, 19, 18, 15, 14, 13, 12, 11, 10, 37, 36, 39, 38, 33, 32, 35, 34, 30, 31, 27, 26, 29, 28, 23, 22, 25, 24, 21, 20, 6, 7, 9, 8, 1, 0, 4, 5, 3, 2],
    ![49, 48, 45, 44, 40, 41, 43, 42, 47, 46, 59, 58, 55, 54, 50, 51, 53, 52, 56, 57, 27, 26, 29, 28, 21, 20, 22, 23, 25, 24, 6, 7, 8, 9, 3, 2, 1, 0, 4, 5, 17, 16, 19, 18, 12, 13, 11, 10, 15, 14, 37, 36, 38, 39, 35, 34, 30, 31, 33, 32],
    ![50, 51, 54, 55, 52, 53, 56, 57, 59, 58, 8, 9, 3, 2, 0, 1, 5, 4, 6, 7, 28, 29, 21, 20, 22, 23, 24, 25, 27, 26, 18, 19, 10, 11, 12, 13, 15, 14, 17, 16, 38, 39, 31, 30, 34, 35, 33, 32, 36, 37, 49, 48, 40, 41, 45, 44, 43, 42, 46, 47],
    ![51, 50, 56, 57, 59, 58, 54, 55, 52, 53, 9, 8, 5, 4, 7, 6, 3, 2, 1, 0, 38, 39, 31, 30, 36, 37, 35, 34, 33, 32, 49, 48, 41, 40, 46, 47, 45, 44, 43, 42, 28, 29, 21, 20, 25, 24, 27, 26, 22, 23, 18, 19, 11, 10, 15, 14, 17, 16, 12, 13],
    ![52, 53, 57, 56, 55, 54, 58, 59, 51, 50, 18, 19, 15, 14, 13, 12, 17, 16, 10, 11, 39, 38, 33, 32, 34, 35, 37, 36, 31, 30, 29, 28, 23, 22, 24, 25, 27, 26, 21, 20, 48, 49, 43, 42, 47, 46, 45, 44, 40, 41, 9, 8, 0, 1, 5, 4, 3, 2, 6, 7],
    ![53, 52, 58, 59, 51, 50, 57, 56, 55, 54, 19, 18, 17, 16, 11, 10, 15, 14, 12, 13, 48, 49, 43, 42, 40, 41, 46, 47, 45, 44, 9, 8, 1, 0, 6, 7, 5, 4, 3, 2, 39, 38, 33, 32, 36, 37, 31, 30, 34, 35, 29, 28, 22, 23, 27, 26, 21, 20, 24, 25],
    ![54, 55, 50, 51, 53, 52, 59, 58, 56, 57, 28, 29, 21, 20, 22, 23, 27, 26, 24, 25, 8, 9, 3, 2, 0, 1, 6, 7, 5, 4, 19, 18, 12, 13, 10, 11, 17, 16, 15, 14, 49, 48, 45, 44, 40, 41, 43, 42, 47, 46, 38, 39, 34, 35, 31, 30, 33, 32, 37, 36],
    ![55, 54, 59, 58, 56, 57, 50, 51, 53, 52, 29, 28, 27, 26, 25, 24, 21, 20, 23, 22, 49, 48, 45, 44, 47, 46, 41, 40, 43, 42, 38, 39, 35, 34, 37, 36, 31, 30, 33, 32, 8, 9, 3, 2, 7, 6, 5, 4, 0, 1, 19, 18, 13, 12, 17, 16, 15, 14, 10, 11],
    ![56, 57, 51, 50, 58, 59, 52, 53, 54, 55, 38, 39, 31, 30, 36, 37, 33, 32, 35, 34, 9, 8, 5, 4, 7, 6, 1, 0, 3, 2, 48, 49, 46, 47, 41, 40, 43, 42, 45, 44, 18, 19, 15, 14, 11, 10, 17, 16, 13, 12, 28, 29, 25, 24, 21, 20, 27, 26, 23, 22],
    ![57, 56, 52, 53, 54, 55, 51, 50, 58, 59, 39, 38, 33, 32, 34, 35, 31, 30, 37, 36, 18, 19, 15, 14, 13, 12, 10, 11, 17, 16, 28, 29, 24, 25, 23, 22, 21, 20, 27, 26, 9, 8, 5, 4, 0, 1, 3, 2, 7, 6, 48, 49, 47, 46, 43, 42, 45, 44, 41, 40],
    ![58, 59, 53, 52, 50, 51, 55, 54, 57, 56, 48, 49, 43, 42, 40, 41, 45, 44, 46, 47, 19, 18, 17, 16, 11, 10, 12, 13, 15, 14, 8, 9, 6, 7, 1, 0, 3, 2, 5, 4, 29, 28, 27, 26, 22, 23, 21, 20, 25, 24, 39, 38, 36, 37, 33, 32, 31, 30, 35, 34],
    ![59, 58, 55, 54, 57, 56, 53, 52, 50, 51, 49, 48, 45, 44, 47, 46, 43, 42, 41, 40, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 19, 18, 17, 16, 13, 12, 15, 14, 11, 10, 8, 9, 7, 6, 3, 2, 5, 4, 1, 0]]

theorem L60_length : L60.length = 60 := rfl

/-- The `i`-th listed element. -/
def el (i : Fin 60) : Equiv.Perm (Fin 6) :=
  L60.get (Fin.cast L60_length.symm i)

theorem el_mem (i : Fin 60) : el i ∈ L60 := L60.get_mem _

theorem mem_iff_el {g : Equiv.Perm (Fin 6)} (hg : g ∈ L60) :
    ∃ i : Fin 60, g = el i := by
  obtain ⟨n, hn⟩ := List.mem_iff_get.mp hg
  exact ⟨Fin.cast L60_length n, by rw [← hn]; rfl⟩

set_option maxHeartbeats 2000000 in
/-- Pointwise bridge: the listed permutations evaluate by the forward
table. -/
theorem el_apply : ∀ (i : Fin 60) (x : Fin 6), el i x = rowF i x := by
  decide

set_option maxHeartbeats 2000000 in
/-- Pointwise bridge for inverses. -/
theorem el_inv_apply :
    ∀ (i : Fin 60) (x : Fin 6), (el i)⁻¹ x = rowI i x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band0 :
    ∀ (i j : Fin 60) (x : Fin 6), i.val < 5 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band1 :
    ∀ (i j : Fin 60) (x : Fin 6), 5 ≤ i.val → i.val < 10 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band2 :
    ∀ (i j : Fin 60) (x : Fin 6), 10 ≤ i.val → i.val < 15 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band3 :
    ∀ (i j : Fin 60) (x : Fin 6), 15 ≤ i.val → i.val < 20 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band4 :
    ∀ (i j : Fin 60) (x : Fin 6), 20 ≤ i.val → i.val < 25 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band5 :
    ∀ (i j : Fin 60) (x : Fin 6), 25 ≤ i.val → i.val < 30 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band6 :
    ∀ (i j : Fin 60) (x : Fin 6), 30 ≤ i.val → i.val < 35 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band7 :
    ∀ (i j : Fin 60) (x : Fin 6), 35 ≤ i.val → i.val < 40 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band8 :
    ∀ (i j : Fin 60) (x : Fin 6), 40 ≤ i.val → i.val < 45 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band9 :
    ∀ (i j : Fin 60) (x : Fin 6), 45 ≤ i.val → i.val < 50 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band10 :
    ∀ (i j : Fin 60) (x : Fin 6), 50 ≤ i.val → i.val < 55 →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
private theorem rowF_mul_band11 :
    ∀ (i j : Fin 60) (x : Fin 6), 55 ≤ i.val →
      rowF i (rowF j x) = rowF (mulT i j) x := by decide

/-- The forward table composes by the precomputed index table: the
exhaustive `3600 × 6` check, recombined from the twelve five-row bands.
A single kernel `decide` over all sixty rows exceeds this box's memory
guard (measured > 12 GB resident, > 540 s); each five-row band stays under
8.4 GB and 20 s of kernel time, and the band cache is freed between
declarations. -/
theorem rowF_mul :
    ∀ (i j : Fin 60) (x : Fin 6),
      rowF i (rowF j x) = rowF (mulT i j) x := by
  intro i j x
  rcases Nat.lt_or_ge i.val 5 with h | h0
  · exact rowF_mul_band0 i j x h
  rcases Nat.lt_or_ge i.val 10 with h | h1
  · exact rowF_mul_band1 i j x h0 h
  rcases Nat.lt_or_ge i.val 15 with h | h2
  · exact rowF_mul_band2 i j x h1 h
  rcases Nat.lt_or_ge i.val 20 with h | h3
  · exact rowF_mul_band3 i j x h2 h
  rcases Nat.lt_or_ge i.val 25 with h | h4
  · exact rowF_mul_band4 i j x h3 h
  rcases Nat.lt_or_ge i.val 30 with h | h5
  · exact rowF_mul_band5 i j x h4 h
  rcases Nat.lt_or_ge i.val 35 with h | h6
  · exact rowF_mul_band6 i j x h5 h
  rcases Nat.lt_or_ge i.val 40 with h | h7
  · exact rowF_mul_band7 i j x h6 h
  rcases Nat.lt_or_ge i.val 45 with h | h8
  · exact rowF_mul_band8 i j x h7 h
  rcases Nat.lt_or_ge i.val 50 with h | h9
  · exact rowF_mul_band9 i j x h8 h
  rcases Nat.lt_or_ge i.val 55 with h | h10
  · exact rowF_mul_band10 i j x h9 h
  exact rowF_mul_band11 i j x h10

set_option maxHeartbeats 2000000 in
/-- The inverse table is the forward table at the precomputed inverse
index. -/
theorem rowF_inv :
    ∀ (i : Fin 60) (x : Fin 6), rowI i x = rowF (invT i) x := by decide

theorem length_L60 : L60.length = 60 ∧ L60.Nodup := by
  refine ⟨L60_length, ?_⟩
  exact List.Nodup.of_map (fun e => (List.finRange 6).map fun x => e x)
    (by decide)

theorem one_mem : (1 : Equiv.Perm (Fin 6)) ∈ L60 := by
  have h : el 0 = 1 := by
    apply Equiv.ext
    intro x
    rw [el_apply, Equiv.Perm.one_apply]
    revert x
    decide
  exact h ▸ el_mem 0

theorem t_mem : t ∈ L60 := by
  have h : el 13 = t := by
    apply Equiv.ext
    intro x
    rw [el_apply]
    revert x
    decide
  exact h ▸ el_mem 13

theorem s_mem : s ∈ L60 := by
  have h : el 59 = s := by
    apply Equiv.ext
    intro x
    rw [el_apply]
    revert x
    decide
  exact h ▸ el_mem 59

theorem inv_closed : ∀ g ∈ L60, g⁻¹ ∈ L60 := by
  intro g hg
  obtain ⟨i, rfl⟩ := mem_iff_el hg
  have h : (el i)⁻¹ = el (invT i) := by
    apply Equiv.ext
    intro x
    rw [el_inv_apply, el_apply]
    exact rowF_inv i x
  rw [h]
  exact el_mem _

theorem mul_closed : ∀ g ∈ L60, ∀ h ∈ L60, g * h ∈ L60 := by
  intro g hg h hh
  obtain ⟨i, rfl⟩ := mem_iff_el hg
  obtain ⟨j, rfl⟩ := mem_iff_el hh
  have hij : el i * el j = el (mulT i j) := by
    apply Equiv.ext
    intro x
    rw [Equiv.Perm.mul_apply, el_apply, el_apply, el_apply]
    exact rowF_mul i j x
  rw [hij]
  exact el_mem _

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- The listed action is 2-transitive on the six axes. -/
theorem two_transitive :
    ∀ i j k l : Fin 6, i ≠ j → k ≠ l → ∃ g ∈ L60, g i = k ∧ g j = l := by
  have key : ∀ i j k l : Fin 6, i ≠ j → k ≠ l →
      ∃ m : Fin 60, rowF m i = k ∧ rowF m j = l := by decide
  intro i j k l hij hkl
  obtain ⟨m, h1, h2⟩ := key i j k l hij hkl
  exact ⟨el m, el_mem m, by rw [el_apply]; exact h1,
    by rw [el_apply]; exact h2⟩

/-! ## The permutation module and the dimension count -/

/-- Coordinate-sum functional of the permutation module. -/
def sumF : (Fin 6 → ℚ) →ₗ[ℚ] ℚ where
  toFun v := ∑ i, v i
  map_add' v w := by simp [Finset.sum_add_distrib]
  map_smul' c v := by simp [Finset.mul_sum]

/-- The sum-zero hyperplane. -/
def V5 : Submodule ℚ (Fin 6 → ℚ) := LinearMap.ker sumF

/-- A permutation acts linearly on the module by right composition with its
inverse. -/
def pAct (g : Equiv.Perm (Fin 6)) : (Fin 6 → ℚ) →ₗ[ℚ] (Fin 6 → ℚ) where
  toFun v i := v (g⁻¹ i)
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

/-- Invariance under every listed group element. -/
def InvariantUnder (W : Submodule ℚ (Fin 6 → ℚ)) : Prop :=
  ∀ g ∈ L60, ∀ w ∈ W, pAct g w ∈ W

theorem V5_invariant : InvariantUnder V5 := by
  intro g _ w hw
  rw [V5, LinearMap.mem_ker] at hw ⊢
  show ∑ i, w (g⁻¹ i) = 0
  rw [Equiv.sum_comp g⁻¹ w]
  exact hw

theorem finrank_V5 : Module.finrank ℚ V5 = 5 := by
  have hsurj : Function.Surjective sumF := by
    intro q
    refine ⟨fun i => if i = 0 then q else 0, ?_⟩
    simp [sumF, Finset.sum_ite_eq']
  have hrank := LinearMap.finrank_range_add_finrank_ker sumF
  rw [LinearMap.range_eq_top.mpr hsurj] at hrank
  simp only [finrank_top, Module.finrank_self, Module.finrank_pi,
    Fintype.card_fin] at hrank
  show Module.finrank ℚ (LinearMap.ker sumF) = 5
  omega

/-- The dimension-six branch module skeleton: if `V5` is irreducible under
the listed action, the permutation module has no three-dimensional
invariant subspace; in particular `1 + 5` admits no decomposition into two
three-dimensional invariant summands. -/
theorem no_three_dim_invariant
    (hirr : ∀ U : Submodule ℚ (Fin 6 → ℚ),
      U ≤ V5 → InvariantUnder U → U = ⊥ ∨ U = V5)
    (W : Submodule ℚ (Fin 6 → ℚ)) (hW : InvariantUnder W)
    (h3 : Module.finrank ℚ W = 3) : False := by
  have hinf_inv : InvariantUnder (W ⊓ V5) := by
    intro g hg w hw
    exact ⟨hW g hg w hw.1, V5_invariant g hg w hw.2⟩
  have hsum := Submodule.finrank_sup_add_finrank_inf_eq W V5
  have hsup : Module.finrank ℚ ↥(W ⊔ V5) ≤ 6 := by
    have h := Submodule.finrank_le (W ⊔ V5)
    simpa only [Module.finrank_pi, Fintype.card_fin] using h
  have hV5 := finrank_V5
  have h2 : 2 ≤ Module.finrank ℚ ↥(W ⊓ V5) := by omega
  rcases hirr (W ⊓ V5) inf_le_right hinf_inv with hbot | heq
  · rw [hbot] at h2
    simp at h2
  · have hle : V5 ≤ W := by
      rw [← heq]
      exact inf_le_left
    have hmono : Module.finrank ℚ V5 ≤ Module.finrank ℚ W :=
      Submodule.finrank_mono hle
    omega

/-! ## Irreducibility of `V5` over `ℚ`

The declared hypothesis of `no_three_dim_invariant` is discharged here.
The engine is the sharp fiber count of the listed action at the base point
`0`: the sixty elements distribute over the joint conditions `g 0 = k`,
`g j = i` with multiplicities forced by sharp 2-transitivity — ten on the
compatible diagonal (`i = k`, `j = 0`), zero on the two mixed cases, two
otherwise.  A nonzero vector of an invariant subspace is first transported
by a listed element so that it does not vanish at `0` (`exists_to_zero`);
averaging it over a stabilizer coset then yields every difference vector
`e k - e l` inside the subspace, and these span `V5`.  Fixing the base
point and using only forward applications keeps the kernel obligations
small. -/

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Sharp fiber count of the listed action at base point `0`,
kernel-checked over all `6³` index choices. -/
theorem count_fibers :
    ∀ k i j : Fin 6,
      L60.countP (fun g => decide (g 0 = k ∧ g j = i)) =
        if i = k then (if j = 0 then 10 else 0)
        else (if j = 0 then 0 else 2) := by
  decide

/-- Every point is carried to the base point `0` by a listed element. -/
theorem exists_to_zero : ∀ a : Fin 6, ∃ g ∈ L60, g a = 0 := by decide

/-- Stabilizer-coset average of `u` with target `k`, written as a single
list sum with indicator terms so that membership in an invariant submodule
is termwise. -/
def avg (u : Fin 6 → ℚ) (k : Fin 6) : Fin 6 → ℚ :=
  (L60.map (fun g => if g 0 = k then pAct g u else 0)).sum

theorem list_sum_mem (U : Submodule ℚ (Fin 6 → ℚ)) :
    ∀ l : List (Fin 6 → ℚ), (∀ x ∈ l, x ∈ U) → l.sum ∈ U := by
  intro l
  induction l with
  | nil => intro _; rw [List.sum_nil]; exact U.zero_mem
  | cons x t ih =>
      intro h
      rw [List.sum_cons]
      exact U.add_mem (h x (by simp)) (ih fun y hy => h y (by simp [hy]))

theorem avg_mem (U : Submodule ℚ (Fin 6 → ℚ)) (hinv : InvariantUnder U)
    {u : Fin 6 → ℚ} (hu : u ∈ U) (k : Fin 6) : avg u k ∈ U := by
  refine list_sum_mem U _ ?_
  intro x hx
  obtain ⟨g, hg, rfl⟩ := List.mem_map.mp hx
  split_ifs
  · exact hinv g hg u hu
  · exact U.zero_mem

/-- Pointwise value of an indicator list sum, expressed through the fiber
counts of the underlying list. -/
theorem sum_indicator_apply (S : List (Equiv.Perm (Fin 6))) (u : Fin 6 → ℚ)
    (k i : Fin 6) :
    (S.map (fun g => if g 0 = k then pAct g u else 0)).sum i =
      ∑ j : Fin 6,
        (S.countP (fun g => decide (g 0 = k ∧ g j = i)) : ℚ) * u j := by
  induction S with
  | nil => simp
  | cons g S ih =>
      rw [List.map_cons, List.sum_cons, Pi.add_apply, ih]
      simp only [List.countP_cons]
      push_cast
      rw [Finset.sum_congr rfl
        (fun j _ => add_mul (S.countP (fun g => decide (g 0 = k ∧ g j = i)) : ℚ)
          _ (u j)), Finset.sum_add_distrib]
      have hterm : ((if g 0 = k then pAct g u else 0) : Fin 6 → ℚ) i =
          ∑ j : Fin 6,
            (if decide (g 0 = k ∧ g j = i) = true then (1 : ℚ) else 0) * u j := by
        simp only [decide_eq_true_eq]
        by_cases h : g 0 = k
        · have hcond : ∀ j : Fin 6,
              (if g 0 = k ∧ g j = i then (1 : ℚ) else 0) * u j =
                if g⁻¹ i = j then u j else 0 := by
            intro j
            by_cases hj : g⁻¹ i = j
            · subst hj
              rw [if_pos ⟨h, Equiv.apply_symm_apply g i⟩, if_pos rfl, one_mul]
            · have hne : g j ≠ i := fun he =>
                hj (by rw [← he]; exact Equiv.symm_apply_apply g j)
              rw [if_neg fun hc => hne hc.2, if_neg hj, zero_mul]
          rw [Finset.sum_congr rfl fun j _ => hcond j, Finset.sum_ite_eq]
          simp [h]
          rfl
        · simp [h]
      rw [hterm]
      ring

/-- Value of the stabilizer-coset average on a sum-zero vector. -/
theorem avg_apply {u : Fin 6 → ℚ} (hsum : ∑ j : Fin 6, u j = 0) (k i : Fin 6) :
    avg u k i = if i = k then 10 * u 0 else -2 * u 0 := by
  rw [avg, sum_indicator_apply]
  simp only [count_fibers k i]
  by_cases h : i = k
  · simp only [h, if_true]
    have : ∀ j : Fin 6,
        (((if j = 0 then 10 else 0 : ℕ)) : ℚ) * u j =
          if j = 0 then 10 * u j else 0 := by
      intro j; by_cases hj : j = 0 <;> simp [hj]
    rw [Finset.sum_congr rfl fun j _ => this j, Finset.sum_ite_eq']
    simp
  · simp only [h, if_false]
    have : ∀ j : Fin 6,
        (((if j = 0 then 0 else 2 : ℕ)) : ℚ) * u j =
          2 * u j - (if j = 0 then 2 * u j else 0) := by
      intro j; by_cases hj : j = 0 <;> simp [hj]
    rw [Finset.sum_congr rfl fun j _ => this j, Finset.sum_sub_distrib,
      ← Finset.mul_sum, hsum, Finset.sum_ite_eq']
    simp

/-- The explicit difference vector `e k - e l`. -/
def dvec (k l : Fin 6) : Fin 6 → ℚ :=
  fun i => (if i = k then 1 else 0) - (if i = l then 1 else 0)

theorem avg_sub_avg {u : Fin 6 → ℚ} (hsum : ∑ j : Fin 6, u j = 0)
    (k l : Fin 6) :
    avg u k - avg u l = (12 * u 0) • dvec k l := by
  funext i
  rw [Pi.sub_apply, avg_apply hsum, avg_apply hsum, Pi.smul_apply, dvec,
    smul_eq_mul]
  split_ifs <;> ring

theorem dvec_mem (U : Submodule ℚ (Fin 6 → ℚ)) (hU : U ≤ V5)
    (hinv : InvariantUnder U) {u : Fin 6 → ℚ} (hu : u ∈ U)
    (h0 : u 0 ≠ 0) (k l : Fin 6) : dvec k l ∈ U := by
  have hsum : ∑ j : Fin 6, u j = 0 := by
    have h := hU hu
    rw [V5, LinearMap.mem_ker] at h
    exact h
  have h12 : (12 : ℚ) * u 0 ≠ 0 := mul_ne_zero (by norm_num) h0
  have hmem : avg u k - avg u l ∈ U :=
    U.sub_mem (avg_mem U hinv hu k) (avg_mem U hinv hu l)
  rw [avg_sub_avg hsum] at hmem
  have hsc := U.smul_mem (12 * u 0)⁻¹ hmem
  rwa [smul_smul, inv_mul_cancel₀ h12, one_smul] at hsc

/-- `V5` is irreducible over `ℚ` under the listed action: the declared
hypothesis of `no_three_dim_invariant`, discharged. -/
theorem V5_irreducible :
    ∀ U : Submodule ℚ (Fin 6 → ℚ),
      U ≤ V5 → InvariantUnder U → U = ⊥ ∨ U = V5 := by
  intro U hU hinv
  rcases eq_or_ne U ⊥ with hbot | hbot
  · exact Or.inl hbot
  · refine Or.inr (le_antisymm hU ?_)
    obtain ⟨w, hw, hw0⟩ := (Submodule.ne_bot_iff U).mp hbot
    obtain ⟨a, hwa⟩ : ∃ a, w a ≠ 0 := Function.ne_iff.mp hw0
    obtain ⟨g, hg, hga⟩ := exists_to_zero a
    set u : Fin 6 → ℚ := pAct g w with hu_def
    have hu : u ∈ U := hinv g hg w hw
    have ha : u 0 ≠ 0 := by
      have hgi : g⁻¹ 0 = a := by
        rw [← hga]
        exact Equiv.symm_apply_apply g a
      show w (g⁻¹ 0) ≠ 0
      rw [hgi]
      exact hwa
    intro v hv
    have hvs : ∑ j : Fin 6, v j = 0 := by
      rw [V5, LinearMap.mem_ker] at hv
      exact hv
    have hrepr : v = ∑ j : Fin 6, v j • dvec j 0 := by
      funext i
      rw [Finset.sum_apply]
      simp only [Pi.smul_apply, dvec, smul_eq_mul]
      have expand : ∀ j : Fin 6,
          v j * ((if i = j then (1 : ℚ) else 0) - (if i = 0 then 1 else 0)) =
            (if i = j then v j else 0) - (if i = 0 then v j else 0) := by
        intro j; split_ifs <;> ring
      rw [Finset.sum_congr rfl fun j _ => expand j, Finset.sum_sub_distrib,
        Finset.sum_ite_eq]
      by_cases h0 : i = 0 <;> simp [h0, hvs]
    rw [hrepr]
    exact Submodule.sum_mem U fun j _ =>
      U.smul_mem (v j) (dvec_mem U hU hinv hu ha j 0)

/-- Unconditional dimension-six branch exclusion: the permutation module of
the six-axis action has no three-dimensional invariant subspace. -/
theorem no_three_dim_invariant_unconditional
    (W : Submodule ℚ (Fin 6 → ℚ)) (hW : InvariantUnder W)
    (h3 : Module.finrank ℚ W = 3) : False :=
  no_three_dim_invariant V5_irreducible W hW h3

/-- Paper-facing form of the dimension-six branch: `1 ⊕ 5` admits no
decomposition into two three-dimensional invariant summands. -/
theorem no_three_plus_three_split :
    ¬ ∃ W₁ W₂ : Submodule ℚ (Fin 6 → ℚ),
        InvariantUnder W₁ ∧ InvariantUnder W₂ ∧
        Module.finrank ℚ W₁ = 3 ∧ Module.finrank ℚ W₂ = 3 ∧
        W₁ ⊓ W₂ = ⊥ ∧ W₁ ⊔ W₂ = ⊤ := by
  rintro ⟨W₁, _W₂, h₁, _, d₁, _⟩
  exact no_three_dim_invariant_unconditional W₁ h₁ d₁

end OPH.A5SixAxes

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.A5SixAxes.length_L60
#print axioms OPH.A5SixAxes.two_transitive
#print axioms OPH.A5SixAxes.mul_closed
#print axioms OPH.A5SixAxes.inv_closed
#print axioms OPH.A5SixAxes.V5_invariant
#print axioms OPH.A5SixAxes.finrank_V5
#print axioms OPH.A5SixAxes.no_three_dim_invariant
#print axioms OPH.A5SixAxes.count_fibers
#print axioms OPH.A5SixAxes.V5_irreducible
#print axioms OPH.A5SixAxes.no_three_dim_invariant_unconditional
#print axioms OPH.A5SixAxes.no_three_plus_three_split
