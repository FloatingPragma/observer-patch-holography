import CarrierModeEquivariance
import Geometry.ScreenCarrierMapCandidate

open scoped BigOperators Matrix

namespace OPH.GoldenSectorCharacters

open OPH.A5PortAction
open OPH.LocalFaceMaxwellAction
open OPH.ScaledMaxwellStability
open OPH.CarrierModeOscillators
open OPH.CarrierModeEquivariance
open OPH.ScreenCarrierMapCandidate

/-!
# The golden sector split: two three-dimensional eigenspaces of the face
normal operator with Galois-conjugate characters

STATUS.  Kernel `decide` checks on committed finite tables (the sixty
listed port permutations of `A5PortAction.perms`, the twenty oriented
faces, the integer face normal matrix `faceNormalZ`, the golden projector
multiple `projGoldenZ`, its image `goldenImageZ`, and two new `ℤ[φ]`
matrices listed below), transported to exact real linear algebra through
`evalPhi`.  The carrier, its face orientation, and the operator
`N = C Cᵀ` are declared in the imported files.  No register row is
discharged.

WHAT IS PROVED.
1. Exact golden eigenprojectors over `ℤ[φ]`.  With `φ² = φ + 1`,
   `√5 = 2φ - 1`, the two roots of `x² - 6x + 4 = 0` are
   `λ₊ = 3 + √5 = 2 + 2φ` and `λ₋ = 3 - √5 = 4 - 2φ`.  The matrices
   `goldenPlusZ` and `goldenMinusZ` (entries in `Zphi = ℤ × ℤ`, the pair
   `(a, b)` meaning `a + bφ`) are twenty times the spectral projectors of
   `N` onto the eigenspaces `λ₊`, `λ₋`: they sum to twice `projGoldenZ`
   (`plus_add_minus`), their product vanishes (`plus_mul_minus`), each
   squares to twenty times itself (`plus_idem`, `minus_idem`), `N` acts
   on each as the eigenvalue (`normal_plus`, `normal_minus`), both are
   symmetric (`plus_symm`, `minus_symm`), and `goldenMinusZ` is the
   entrywise Galois conjugate `a + bφ ↦ (a + b) - bφ` of `goldenPlusZ`
   (`minus_eq_conj_plus`).  The integer identity
   `5 · goldenPlusZ = (2φ - 1) · goldenImageZ - (6φ - 8) · projGoldenZ`
   (`five_plus_decomp`) exhibits `goldenPlusZ` as the closed-form
   projector `projGolden · (N - λ₋) / (λ₊ - λ₋)` with
   `λ₊ - λ₋ = 2√5 = 4φ - 2`.
2. Traces.  `trace goldenPlusZ = trace goldenMinusZ = 60 = 20 · 3`
   (`plus_trace`, `minus_trace`), so each projector `goldenPlusZ / 20`,
   `goldenMinusZ / 20` has trace `3`.  Inference outside the theorems:
   a symmetric idempotent has rank equal to its trace, so each golden
   eigenspace is three-dimensional.
3. Characters.  `traceZphi Q p = ∑ f, sign_f · Q (p f) f` is the trace of
   the induced signed face permutation against `Q`.  For every listed
   row the values on `goldenPlusZ` are `60, -20, 0` (that is `20 · 3`,
   `20 · (-1)`, `20 · 0`) on the elements of order `1, 2, 3` and either
   `20 φ` or `20 (1 - φ)` on the elements of order five
   (`chi_plus_values`); the same holds on `goldenMinusZ`
   (`chi_minus_values`).  The two values are Galois conjugate on every
   element (`chi_conj`), agree exactly on the elements of order other
   than five (`chi_differ_iff`), and are exchanged by squaring an element
   of order five (`chi_square`): `χ₊(g²) = χ₋(g)`.  The twenty-four
   elements of order five split into twelve with `χ₊ = φ` and twelve with
   `χ₊ = 1 - φ` (`five_class_counts`).  The values are class functions:
   `χ₊(q p q⁻¹) = χ₊(p)` for every pair of rows (`chi_plus_class`).
4. Norms.  `∑_g χ₊(g)² = ∑_g χ₋(g)² = 60` in `ℤ[φ]`, stated with the
   scale `20²` (`character_norms_golden`).  The exact arithmetic is
   `12 φ² + 12 (1 - φ)² + 15 + 20 · 0 + 9 = 12 · 3 + 24 = 60` since
   `φ² + (1 - φ)² = 3`.  Inference outside the theorems: by the
   orthogonality relations, norm `60 = |A5|` reads as irreducibility, so
   the golden sector is `3 ⊕ 3'` with the eigenvalue `λ₊ = 2φ²` carrying
   the character `(3, -1, 0, φ, 1 - φ)` and the eigenvalue
   `λ₋ = 2/φ² = 3 - √5` carrying its Galois conjugate; the eigenvalue is
   conjugated along with the character, since `λ₋` is the conjugate of
   `λ₊` in `ℤ[φ]` (`lamMinus_eq_conj`).
5. Equivariance.  Every listed row fixes both projector tables entrywise,
   `goldenPlusZ (p i) (p j) = goldenPlusZ i j` (`plus_equivariant`,
   `minus_equivariant`), and the induced action matrix `faceActZ p`
   commutes with each (`faceAct_comm_plus`, `faceAct_comm_minus`).
6. Real transport.  `goldenPlusR = evalPhi ∘ goldenPlusZ / 20` and
   `goldenMinusR` are real matrices with `goldenPlusR + goldenMinusR =
   projGoldenR` (`plusR_add_minusR`), `faceNormalR * goldenPlusR =
   (3 + √5) • goldenPlusR` and `faceNormalR * goldenMinusR = (3 - √5) •
   goldenMinusR` (`normalR_plus`, `normalR_minus`), and traces `3`
   (`plusR_trace`, `minusR_trace`).

PRIOR WORK.  `Screen/ScaledMaxwellStability.lean` lists `projGoldenZ`
and `goldenImageZ` and proves the golden sector identities over `ℤ`
(`projGolden_idem`, `projGolden_image`, `projGolden_commute`,
`projGolden_quadratic`, `projGolden_symm`) and over `ℝ`
(`golden_idem`, `golden_commute`, `golden_quadratic_mat`), and exhibits
one seam-level eigenvector of `CᵀC` at `3 + √5` (`goldenMode_eigen`).
`Screen/CarrierModeEquivariance.lean` lists the induced face action of
every row (`faceAct`, `facePerm`, `faceSign`, `faceActZ`), proves the
composition law (`facePerm_comp`), the inverse (`inv_spec`), element
orders and class counts (`elemOrder_spec`, `order_counts`), projector
equivariance (`projectors_equivariant`), and the characters of the five
integer projectors with their norms (`projector_characters`,
`character_norms`); on the golden sector it obtains the sum
`6, -2, 0, 1` with norm `120` and states the split `3 ⊕ 3'` and its
assignment to the two eigenvalues as an observation.
`Screen/A5CharacterField.lean` treats the Galois conjugation of the
`A5` character field on abstract multiplicity vectors (`GaloisStable`).
`Geometry/ScreenCarrierMapCandidate.lean` defines `Zphi`, `zmul`, and
`evalPhi` with the evaluation lemmas used in section 6.  The content of
this file beyond the prior modules is: the two `ℤ[φ]` projector tables,
their projector, eigenvalue, symmetry, and conjugation identities, the
per-projector characters with the order-five class split, the Galois
conjugation of the characters and of the eigenvalues, the class function
property, the norms `60`, the equivariance of each projector, and the
real transport.

ROWS TOUCHED (none discharged).  Source clock and duration row: untouched
(no step appears).  Physical spacetime attachment row: the listed group is
an incidence automorphism group of the declared carrier; its
identification with a physical rotation group is open.  Light-signal row:
the identification of a golden mode with a physical oscillation is open;
the theorems here sort the two golden eigenspaces by symmetry type.
Coupled-action row: untouched.  Laboratory clock and energy calibration
import: untouched.  Gravitation-route energy identification: untouched.

NEGATIVES CITED.  Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`), at scope only:
realized histories select no velocity curvature or Legendre map, so every
shape of the declared evolution whose operator `N` appears here is a
declared enrichment; the spectral statements concern the declared
operator.

CONVENTIONS.  `Zphi = ℤ × ℤ` with `(a, b) = a + bφ`, `φ² = φ + 1`,
`√5 = 2φ - 1`; conjugation `zconj (a, b) = (a + b, -b)`; integer scaling
`zint n (a, b) = (n a, n b)`.  `λ₊ = (2, 2)`, `λ₋ = (4, -2)`.  Rows of
`perms` act on faces by `facePerm` with `faceSign = 1`
(`faceSign_eq_one`); `traceZphi Q p = ∑ f, zint (faceSign p f)
(Q (facePerm p f) f)`.  Products of `Zphi` matrices are written out as
`∑ k, zmul (A i k) (B k j)`; products with an integer matrix as
`∑ k, zint (N i k) (A k j)`.  The two order-five classes are identified
by the value of `χ₊`.

FALSIFIER.  An entry of `goldenPlusZ + goldenMinusZ` off `2 projGoldenZ`,
a nonzero entry of the product, an entry of `N goldenPlusZ` off
`λ₊ goldenPlusZ`, a row whose trace against `goldenPlusZ` is off the
listed values, an order-five row with `χ₊(g²) ≠ χ₋(g)`, a class count
other than `12, 12`, or a norm other than `24000 = 60 · 20²` would make
the corresponding kernel check fail.

Axiom audit.  The `#print axioms` lines at the end of the file show at
most `propext`, `Classical.choice`, and `Quot.sound`.
-/

/-! ## 1. The two golden projector tables over `ℤ[φ]` -/

/-- Integer scaling in `ℤ[φ]`. -/
def zint (n : ℤ) (x : Zphi) : Zphi := (n * x.1, n * x.2)

/-- Galois conjugation `φ ↦ 1 - φ` on `ℤ[φ]`: `a + bφ ↦ (a + b) - bφ`. -/
def zconj (x : Zphi) : Zphi := (x.1 + x.2, -x.2)

/-- `λ₊ = 3 + √5 = 2 + 2φ`. -/
def lamPlus : Zphi := (2, 2)

/-- `λ₋ = 3 - √5 = 4 - 2φ`. -/
def lamMinus : Zphi := (4, -2)

theorem lamMinus_eq_conj : lamMinus = zconj lamPlus := by decide

/-- The rational part of twenty times the `λ₊` projector. -/
def plusA : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![3, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -3],
      ![1, 3, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -1, -1, -1, -3, -1],
      ![1, 1, 3, 1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -3, -1, -1],
      ![1, 1, 1, 3, 1, -1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -3, 1, -1, -1, -1],
      ![1, 1, 1, 1, 3, -1, -1, -1, 1, -1, 1, -1, 1, 1, -3, -1, 1, -1, -1, -1],
      ![1, 1, 1, -1, -1, 3, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -3, -1, -1, -1],
      ![1, 1, -1, 1, -1, 1, 3, 1, -1, -1, 1, 1, -1, -3, 1, -1, -1, 1, -1, -1],
      ![1, 1, -1, -1, -1, 1, 1, 3, -1, 1, -1, 1, -3, -1, 1, 1, -1, 1, -1, -1],
      ![1, -1, 1, -1, 1, 1, -1, -1, 3, 1, -1, -3, 1, 1, -1, 1, -1, -1, 1, -1],
      ![1, -1, 1, -1, -1, 1, -1, 1, 1, 3, -3, -1, -1, 1, 1, 1, -1, -1, 1, -1],
      ![-1, 1, -1, 1, 1, -1, 1, -1, -1, -3, 3, 1, 1, -1, -1, -1, 1, 1, -1, 1],
      ![-1, 1, -1, 1, -1, -1, 1, 1, -3, -1, 1, 3, -1, -1, 1, -1, 1, 1, -1, 1],
      ![-1, -1, 1, 1, 1, -1, -1, -3, 1, -1, 1, -1, 3, 1, -1, -1, 1, -1, 1, 1],
      ![-1, -1, 1, -1, 1, -1, -3, -1, 1, 1, -1, -1, 1, 3, -1, 1, 1, -1, 1, 1],
      ![-1, -1, -1, -1, -3, 1, 1, 1, -1, 1, -1, 1, -1, -1, 3, 1, -1, 1, 1, 1],
      ![-1, -1, -1, -3, -1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 3, -1, 1, 1, 1],
      ![-1, -1, -1, 1, 1, -3, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, 3, 1, 1, 1],
      ![-1, -1, -3, -1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 3, 1, 1],
      ![-1, -3, -1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 3, 1],
      ![-3, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3]]

/-- The `φ` part of twenty times the `λ₊` projector. -/
def plusB : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![0, -2, -2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0],
      ![-2, 0, 0, -2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 2],
      ![-2, 0, 0, 0, -2, 0, 0, 0, -2, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 2],
      ![0, -2, 0, 0, -2, 0, 0, 0, 0, 2, -2, 0, 0, 0, 2, 0, 0, 0, 2, 0],
      ![0, 0, -2, -2, 0, 0, 0, 2, 0, 0, 0, 0, -2, 0, 0, 2, 0, 2, 0, 0],
      ![-2, 0, 0, 0, 0, 0, 0, -2, 0, -2, 2, 0, 2, 0, 0, 0, 0, 0, 0, 2],
      ![0, -2, 0, 0, 0, 0, 0, -2, 2, 0, 0, -2, 2, 0, 0, 0, 0, 0, 2, 0],
      ![0, 0, 0, 0, 2, -2, -2, 0, 0, 0, 0, 0, 0, 2, -2, 0, 2, 0, 0, 0],
      ![0, 0, -2, 0, 0, 0, 2, 0, 0, -2, 2, 0, 0, -2, 0, 0, 0, 2, 0, 0],
      ![0, 0, 0, 2, 0, -2, 0, 0, -2, 0, 0, 2, 0, 0, 0, -2, 2, 0, 0, 0],
      ![0, 0, 0, -2, 0, 2, 0, 0, 2, 0, 0, -2, 0, 0, 0, 2, -2, 0, 0, 0],
      ![0, 0, 2, 0, 0, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 0, 0, -2, 0, 0],
      ![0, 0, 0, 0, -2, 2, 2, 0, 0, 0, 0, 0, 0, -2, 2, 0, -2, 0, 0, 0],
      ![0, 2, 0, 0, 0, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 0, 0, 0, -2, 0],
      ![0, 0, 2, 2, 0, 0, 0, -2, 0, 0, 0, 0, 2, 0, 0, -2, 0, -2, 0, 0],
      ![0, 2, 0, 0, 2, 0, 0, 0, 0, -2, 2, 0, 0, 0, -2, 0, 0, 0, -2, 0],
      ![2, 0, 0, 0, 0, 0, 0, 2, 0, 2, -2, 0, -2, 0, 0, 0, 0, 0, 0, -2],
      ![2, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, -2, 0, 0, -2, 0, 0, 0, 0, -2],
      ![2, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, 0, 0, -2],
      ![0, 2, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, -2, 0]]

/-- The rational part of twenty times the `λ₋` projector. -/
def minusA : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![3, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, -3],
      ![-1, 3, 1, -1, 1, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -3, 1],
      ![-1, 1, 3, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, 1, 1, -1, -1, -3, -1, 1],
      ![1, -1, 1, 3, -1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -3, 1, -1, 1, -1],
      ![1, 1, -1, -1, 3, -1, -1, 1, 1, -1, 1, -1, -1, 1, -3, 1, 1, 1, -1, -1],
      ![-1, 1, 1, -1, -1, 3, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -3, -1, -1, 1],
      ![1, -1, -1, 1, -1, 1, 3, -1, 1, -1, 1, -1, 1, -3, 1, -1, -1, 1, 1, -1],
      ![1, 1, -1, -1, 1, -1, -1, 3, -1, 1, -1, 1, -3, 1, -1, 1, 1, 1, -1, -1],
      ![1, -1, -1, -1, 1, 1, 1, -1, 3, -1, 1, -3, 1, -1, -1, 1, -1, 1, 1, -1],
      ![1, -1, 1, 1, -1, -1, -1, 1, -1, 3, -3, 1, -1, 1, 1, -1, 1, -1, 1, -1],
      ![-1, 1, -1, -1, 1, 1, 1, -1, 1, -3, 3, -1, 1, -1, -1, 1, -1, 1, -1, 1],
      ![-1, 1, 1, 1, -1, -1, -1, 1, -3, 1, -1, 3, -1, 1, 1, -1, 1, -1, -1, 1],
      ![-1, -1, 1, 1, -1, 1, 1, -3, 1, -1, 1, -1, 3, -1, 1, -1, -1, -1, 1, 1],
      ![-1, 1, 1, -1, 1, -1, -3, 1, -1, 1, -1, 1, -1, 3, -1, 1, 1, -1, -1, 1],
      ![-1, -1, 1, 1, -3, 1, 1, -1, -1, 1, -1, 1, 1, -1, 3, -1, -1, -1, 1, 1],
      ![-1, 1, -1, -3, 1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 3, -1, 1, -1, 1],
      ![1, -1, -1, 1, 1, -3, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 3, 1, 1, -1],
      ![1, -1, -3, -1, 1, -1, 1, 1, 1, -1, 1, -1, -1, -1, -1, 1, 1, 3, 1, -1],
      ![1, -3, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 3, -1],
      ![-3, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, 3]]

/-- The `φ` part of twenty times the `λ₋` projector. -/
def minusB : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of
    ![![0, 2, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, -2, 0],
      ![2, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, 0, 0, -2],
      ![2, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, -2, 0, 0, -2, 0, 0, 0, 0, -2],
      ![0, 2, 0, 0, 2, 0, 0, 0, 0, -2, 2, 0, 0, 0, -2, 0, 0, 0, -2, 0],
      ![0, 0, 2, 2, 0, 0, 0, -2, 0, 0, 0, 0, 2, 0, 0, -2, 0, -2, 0, 0],
      ![2, 0, 0, 0, 0, 0, 0, 2, 0, 2, -2, 0, -2, 0, 0, 0, 0, 0, 0, -2],
      ![0, 2, 0, 0, 0, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 0, 0, 0, -2, 0],
      ![0, 0, 0, 0, -2, 2, 2, 0, 0, 0, 0, 0, 0, -2, 2, 0, -2, 0, 0, 0],
      ![0, 0, 2, 0, 0, 0, -2, 0, 0, 2, -2, 0, 0, 2, 0, 0, 0, -2, 0, 0],
      ![0, 0, 0, -2, 0, 2, 0, 0, 2, 0, 0, -2, 0, 0, 0, 2, -2, 0, 0, 0],
      ![0, 0, 0, 2, 0, -2, 0, 0, -2, 0, 0, 2, 0, 0, 0, -2, 2, 0, 0, 0],
      ![0, 0, -2, 0, 0, 0, 2, 0, 0, -2, 2, 0, 0, -2, 0, 0, 0, 2, 0, 0],
      ![0, 0, 0, 0, 2, -2, -2, 0, 0, 0, 0, 0, 0, 2, -2, 0, 2, 0, 0, 0],
      ![0, -2, 0, 0, 0, 0, 0, -2, 2, 0, 0, -2, 2, 0, 0, 0, 0, 0, 2, 0],
      ![0, 0, -2, -2, 0, 0, 0, 2, 0, 0, 0, 0, -2, 0, 0, 2, 0, 2, 0, 0],
      ![0, -2, 0, 0, -2, 0, 0, 0, 0, 2, -2, 0, 0, 0, 2, 0, 0, 0, 2, 0],
      ![-2, 0, 0, 0, 0, 0, 0, -2, 0, -2, 2, 0, 2, 0, 0, 0, 0, 0, 0, 2],
      ![-2, 0, 0, 0, -2, 0, 0, 0, -2, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 2],
      ![-2, 0, 0, -2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 2],
      ![0, -2, -2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0]]

/-- Twenty times the spectral projector of `N` onto the eigenvalue `λ₊`,
entries in `ℤ[φ]`. -/
def goldenPlusZ : Matrix (Fin 20) (Fin 20) Zphi :=
  Matrix.of fun i j ↦ (plusA i j, plusB i j)

/-- Twenty times the spectral projector of `N` onto the eigenvalue `λ₋`,
entries in `ℤ[φ]`. -/
def goldenMinusZ : Matrix (Fin 20) (Fin 20) Zphi :=
  Matrix.of fun i j ↦ (minusA i j, minusB i j)

/-! ## 2. Projector identities (kernel decide) -/

theorem plus_add_minus :
    ∀ i j : Fin 20, goldenPlusZ i j + goldenMinusZ i j = (2 * projGoldenZ i j, 0) := by
  decide

theorem minus_eq_conj_plus : ∀ i j : Fin 20, goldenMinusZ i j = zconj (goldenPlusZ i j) := by
  decide

theorem plus_symm : ∀ i j : Fin 20, goldenPlusZ i j = goldenPlusZ j i := by decide
theorem minus_symm : ∀ i j : Fin 20, goldenMinusZ i j = goldenMinusZ j i := by decide

theorem plus_trace : (∑ i : Fin 20, goldenPlusZ i i) = (60, 0) := by decide
theorem minus_trace : (∑ i : Fin 20, goldenMinusZ i i) = (60, 0) := by decide

set_option maxRecDepth 16384 in
theorem plus_mul_minus :
    ∀ i j : Fin 20, (∑ k : Fin 20, zmul (goldenPlusZ i k) (goldenMinusZ k j)) = 0 := by
  decide

set_option maxRecDepth 16384 in
theorem plus_idem :
    ∀ i j : Fin 20, (∑ k : Fin 20, zmul (goldenPlusZ i k) (goldenPlusZ k j)) =
      zint 20 (goldenPlusZ i j) := by
  decide

set_option maxRecDepth 16384 in
theorem minus_idem :
    ∀ i j : Fin 20, (∑ k : Fin 20, zmul (goldenMinusZ i k) (goldenMinusZ k j)) =
      zint 20 (goldenMinusZ i j) := by
  decide

set_option maxRecDepth 16384 in
/-- `N goldenPlusZ = λ₊ goldenPlusZ`. -/
theorem normal_plus :
    ∀ i j : Fin 20, (∑ k : Fin 20, zint (faceNormalZ i k) (goldenPlusZ k j)) =
      zmul lamPlus (goldenPlusZ i j) := by
  decide

set_option maxRecDepth 16384 in
/-- `N goldenMinusZ = λ₋ goldenMinusZ`. -/
theorem normal_minus :
    ∀ i j : Fin 20, (∑ k : Fin 20, zint (faceNormalZ i k) (goldenMinusZ k j)) =
      zmul lamMinus (goldenMinusZ i j) := by
  decide

/-- Closed form: `5 · goldenPlusZ = (2φ - 1) · goldenImageZ - (6φ - 8) · projGoldenZ`,
the projector `projGolden · (N - λ₋) / (λ₊ - λ₋)` with `1/√5 = (2φ - 1)/5`. -/
theorem five_plus_decomp :
    ∀ i j : Fin 20, zint 5 (goldenPlusZ i j) =
      zsub (zmul (-1, 2) (goldenImageZ i j, 0)) (zmul (-8, 6) (projGoldenZ i j, 0)) := by
  decide

/-- Closed form for the conjugate: `5 · goldenMinusZ = (1 - 2φ) · goldenImageZ - (-2 - 6φ) · projGoldenZ`. -/
theorem five_minus_decomp :
    ∀ i j : Fin 20, zint 5 (goldenMinusZ i j) =
      zsub (zmul (1, -2) (goldenImageZ i j, 0)) (zmul (-2, -6) (projGoldenZ i j, 0)) := by
  decide

/-! ## 3. Characters of the two golden eigenspaces -/

/-- Trace of the induced signed face permutation of row `p` against a
`ℤ[φ]` face matrix `Q`: `∑ f, sign_f · Q (p f) f`. -/
def traceZphi (Q : Matrix (Fin 20) (Fin 20) Zphi) (p : List Nat) : Zphi :=
  ∑ f : Fin 20, zint (faceSign p f) (Q (facePerm p f) f)

/-- **(3)** Character values on the `λ₊` eigenspace, scaled by `20`:
`3` on the identity, `-1` on order two, `0` on order three, and `φ` or
`1 - φ` on order five. -/
theorem chi_plus_values : ∀ p ∈ perms,
    (elemOrder p = 1 → traceZphi goldenPlusZ p = (60, 0)) ∧
    (elemOrder p = 2 → traceZphi goldenPlusZ p = (-20, 0)) ∧
    (elemOrder p = 3 → traceZphi goldenPlusZ p = 0) ∧
    (elemOrder p = 5 →
      traceZphi goldenPlusZ p = (0, 20) ∨ traceZphi goldenPlusZ p = (20, -20)) := by
  decide +kernel

/-- **(3)** Character values on the `λ₋` eigenspace, scaled by `20`. -/
theorem chi_minus_values : ∀ p ∈ perms,
    (elemOrder p = 1 → traceZphi goldenMinusZ p = (60, 0)) ∧
    (elemOrder p = 2 → traceZphi goldenMinusZ p = (-20, 0)) ∧
    (elemOrder p = 3 → traceZphi goldenMinusZ p = 0) ∧
    (elemOrder p = 5 →
      traceZphi goldenMinusZ p = (0, 20) ∨ traceZphi goldenMinusZ p = (20, -20)) := by
  decide +kernel

/-- **(3)** The character of the `λ₋` eigenspace is the Galois conjugate of
the character of the `λ₊` eigenspace on every listed row. -/
theorem chi_conj : ∀ p ∈ perms,
    traceZphi goldenMinusZ p = zconj (traceZphi goldenPlusZ p) := by
  decide +kernel

/-- **(3)** The two characters agree exactly on the rows of order other
than five. -/
theorem chi_differ_iff : ∀ p ∈ perms,
    (traceZphi goldenPlusZ p = traceZphi goldenMinusZ p ↔ elemOrder p ≠ 5) := by
  decide +kernel

/-- **(3)** Squaring an element of order five exchanges the two classes:
`χ₊(g²) = χ₋(g)`, the Galois conjugation `√5 ↦ -√5` realised inside the group. -/
theorem chi_square : ∀ p ∈ perms, elemOrder p = 5 →
    traceZphi goldenPlusZ (comp p p) = traceZphi goldenMinusZ p := by
  decide +kernel

/-- **(3)** The twenty-four rows of order five split into twelve with
`χ₊ = φ` and twelve with `χ₊ = 1 - φ` (each scaled by `20`). -/
theorem five_class_counts :
    (perms.filter fun p => traceZphi goldenPlusZ p = (0, 20)).length = 12 ∧
    (perms.filter fun p => traceZphi goldenPlusZ p = (20, -20)).length = 12 ∧
    (perms.filter fun p => elemOrder p = 5 ∧ traceZphi goldenPlusZ p = (0, 20)).length = 12 ∧
    (perms.filter fun p => elemOrder p = 5 ∧ traceZphi goldenPlusZ p = (20, -20)).length = 12 := by
  decide +kernel

/-- **(4)** Character norms `∑_g χ±(g)² = 60`, at scale `20²`: the sum of
`zmul (χ g) (χ g)` over the sixty rows is `(24000, 0) = 60 · 20² + 0 φ`. -/
theorem character_norms_golden :
    (perms.map fun p => zmul (traceZphi goldenPlusZ p) (traceZphi goldenPlusZ p)).sum =
      (24000, 0) ∧
    (perms.map fun p => zmul (traceZphi goldenMinusZ p) (traceZphi goldenMinusZ p)).sum =
      (24000, 0) := by
  decide +kernel

/-! ## 4. Equivariance of the two projectors -/

theorem golden_checks : ∀ p ∈ perms,
    (projCheck plusA p && projCheck plusB p && projCheck minusA p && projCheck minusB p) =
      true := by
  decide +kernel

/-- **(5)** Every listed row fixes the `λ₊` projector table entrywise. -/
theorem plus_equivariant (p : List Nat) (hp : p ∈ perms) (i j : Fin 20) :
    goldenPlusZ (facePerm p i) (facePerm p j) = goldenPlusZ i j := by
  have h := golden_checks p hp
  simp only [Bool.and_eq_true] at h
  obtain ⟨⟨⟨hA, hB⟩, _⟩, _⟩ := h
  simp only [goldenPlusZ, Matrix.of_apply]
  rw [projCheck_entry _ p hA i j, projCheck_entry _ p hB i j]

/-- **(5)** Every listed row fixes the `λ₋` projector table entrywise. -/
theorem minus_equivariant (p : List Nat) (hp : p ∈ perms) (i j : Fin 20) :
    goldenMinusZ (facePerm p i) (facePerm p j) = goldenMinusZ i j := by
  have h := golden_checks p hp
  simp only [Bool.and_eq_true] at h
  obtain ⟨⟨_, hA⟩, hB⟩ := h
  simp only [goldenMinusZ, Matrix.of_apply]
  rw [projCheck_entry _ p hA i j, projCheck_entry _ p hB i j]

theorem zint_one (x : Zphi) : zint 1 x = x := by
  obtain ⟨a, b⟩ := x
  simp [zint]

theorem zint_zero (x : Zphi) : zint 0 x = 0 := by
  obtain ⟨a, b⟩ := x
  simp [zint]

/-- Left multiplication by the action matrix moves the row index. -/
theorem faceAct_mul_left (p : List Nat) (hp : p ∈ perms) (Q : Matrix (Fin 20) (Fin 20) Zphi)
    (i j : Fin 20) :
    (∑ k : Fin 20, zint (faceActZ p i k) (Q k j)) = Q (facePerm p i) j := by
  simp only [faceActZ, Matrix.of_apply]
  rw [Finset.sum_eq_single (facePerm p i)]
  · rw [if_pos rfl, faceSign_eq_one p hp, zint_one]
  · intro k _ hk
    rw [if_neg (Ne.symm hk), zint_zero]
  · intro h; exact absurd (Finset.mem_univ _) h

/-- The face permutation applied to the inverse image of `j` returns `j`. -/
theorem facePerm_symm_apply (p : List Nat) (hp : p ∈ perms) (j : Fin 20) :
    facePerm p ((faceEquiv p hp).symm j) = j := by
  rw [← faceEquiv_apply p hp, Equiv.apply_symm_apply]

/-- Right multiplication by the action matrix moves the column index through
the inverse of the face permutation. -/
theorem faceAct_mul_right (p : List Nat) (hp : p ∈ perms) (Q : Matrix (Fin 20) (Fin 20) Zphi)
    (i j : Fin 20) :
    (∑ k : Fin 20, zint (faceActZ p k j) (Q i k)) = Q i ((faceEquiv p hp).symm j) := by
  simp only [faceActZ, Matrix.of_apply]
  rw [Finset.sum_eq_single ((faceEquiv p hp).symm j)]
  · rw [if_pos (facePerm_symm_apply p hp j), faceSign_eq_one p hp, zint_one]
  · intro k _ hk
    have hne : facePerm p k ≠ j := by
      intro h
      apply hk
      rw [← faceEquiv_apply p hp] at h
      rw [← h, Equiv.symm_apply_apply]
    rw [if_neg hne, zint_zero]
  · intro h; exact absurd (Finset.mem_univ _) h

/-- **(5)** The action matrix of every listed row commutes with the `λ₊`
projector: `faceActZ p · goldenPlusZ = goldenPlusZ · faceActZ p`. -/
theorem faceAct_comm_plus (p : List Nat) (hp : p ∈ perms) (i j : Fin 20) :
    (∑ k : Fin 20, zint (faceActZ p i k) (goldenPlusZ k j)) =
      ∑ k : Fin 20, zint (faceActZ p k j) (goldenPlusZ i k) := by
  rw [faceAct_mul_left p hp, faceAct_mul_right p hp]
  have h := plus_equivariant p hp i ((faceEquiv p hp).symm j)
  rw [facePerm_symm_apply p hp j] at h
  exact h

/-- **(5)** The action matrix of every listed row commutes with the `λ₋`
projector. -/
theorem faceAct_comm_minus (p : List Nat) (hp : p ∈ perms) (i j : Fin 20) :
    (∑ k : Fin 20, zint (faceActZ p i k) (goldenMinusZ k j)) =
      ∑ k : Fin 20, zint (faceActZ p k j) (goldenMinusZ i k) := by
  rw [faceAct_mul_left p hp, faceAct_mul_right p hp]
  have h := minus_equivariant p hp i ((faceEquiv p hp).symm j)
  rw [facePerm_symm_apply p hp j] at h
  exact h

/-- The trace against a listed row drops the sign, which is `1`. -/
theorem traceZphi_eq (Q : Matrix (Fin 20) (Fin 20) Zphi) (p : List Nat) (hp : p ∈ perms) :
    traceZphi Q p = ∑ f : Fin 20, Q (facePerm p f) f := by
  unfold traceZphi
  refine Finset.sum_congr rfl fun f _ ↦ ?_
  rw [faceSign_eq_one p hp, zint_one]

/-- **(3)** Class function: `χ₊(q p q⁻¹) = χ₊(p)` for every pair of listed
rows, from the composition law and the equivariance of the projector. -/
theorem chi_plus_class (p : List Nat) (hp : p ∈ perms) (q : List Nat) (hq : q ∈ perms) :
    traceZphi goldenPlusZ (comp (comp q p) (invPerm q)) = traceZphi goldenPlusZ p := by
  obtain ⟨hq', hqq', hq'q⟩ := inv_spec q hq
  have hqp := comp_mem q hq p hp
  have hall := comp_mem _ hqp _ hq'
  rw [traceZphi_eq _ _ hall, traceZphi_eq _ _ hp]
  have hcomp : ∀ f, facePerm (comp (comp q p) (invPerm q)) f =
      facePerm q (facePerm p (facePerm (invPerm q) f)) := by
    intro f
    rw [facePerm_comp _ hqp _ hq', facePerm_comp q hq p hp]
  simp only [hcomp]
  rw [← Fintype.sum_equiv (faceEquiv q hq) (fun g ↦ goldenPlusZ (facePerm p g) g)
    (fun f ↦ goldenPlusZ (facePerm q (facePerm p (facePerm (invPerm q) f))) f)]
  intro g
  rw [faceEquiv_apply]
  have hback : facePerm (invPerm q) (facePerm q g) = g := by
    rw [← facePerm_comp _ hq' q hq, hq'q, facePerm_id]
  rw [hback, plus_equivariant q hq]

/-- **(3)** Class function for the conjugate character. -/
theorem chi_minus_class (p : List Nat) (hp : p ∈ perms) (q : List Nat) (hq : q ∈ perms) :
    traceZphi goldenMinusZ (comp (comp q p) (invPerm q)) = traceZphi goldenMinusZ p := by
  have hall := comp_mem _ (comp_mem q hq p hp) _ (inv_spec q hq).1
  rw [chi_conj _ hall, chi_conj p hp, chi_plus_class p hp q hq]

/-! ## 5. Real transport through `evalPhi` -/

noncomputable section

theorem evalPhi_zero : evalPhi (0 : Zphi) = 0 := by
  simp [evalPhi]

theorem evalPhi_sum {ι : Type*} (s : Finset ι) (g : ι → Zphi) :
    evalPhi (∑ i ∈ s, g i) = ∑ i ∈ s, evalPhi (g i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [evalPhi_zero]
  | insert a s ha ih => rw [Finset.sum_insert ha, Finset.sum_insert ha, evalPhi_add, ih]

theorem evalPhi_zint (n : ℤ) (x : Zphi) : evalPhi (zint n x) = (n : ℝ) * evalPhi x := by
  obtain ⟨a, b⟩ := x
  simp only [evalPhi, zint]
  push_cast
  ring

theorem evalPhi_lamPlus : evalPhi lamPlus = 3 + Real.sqrt 5 := by
  rw [three_add_sqrt5_eq_two_goldenRatio_sq, Real.goldenRatio_sq]
  simp only [evalPhi, lamPlus]
  push_cast
  ring

theorem evalPhi_lamMinus : evalPhi lamMinus = 3 - Real.sqrt 5 := by
  have h : Real.sqrt 5 = 2 * Real.goldenRatio - 1 := by
    unfold Real.goldenRatio; ring
  rw [h]
  simp only [evalPhi, lamMinus]
  push_cast
  ring

/-- The real spectral projector of `N` for the eigenvalue `3 + √5`. -/
def goldenPlusR : Matrix (Fin 20) (Fin 20) ℝ :=
  Matrix.of fun i j ↦ evalPhi (goldenPlusZ i j) / 20

/-- The real spectral projector of `N` for the eigenvalue `3 - √5`. -/
def goldenMinusR : Matrix (Fin 20) (Fin 20) ℝ :=
  Matrix.of fun i j ↦ evalPhi (goldenMinusZ i j) / 20

/-- **(6)** The two real projectors resolve the golden projector. -/
theorem plusR_add_minusR : goldenPlusR + goldenMinusR = projGoldenR := by
  ext i j
  simp only [goldenPlusR, goldenMinusR, projGoldenR, scaledProj, Matrix.add_apply,
    Matrix.of_apply, Matrix.smul_apply, castZ_apply, smul_eq_mul]
  rw [← add_div, ← evalPhi_add, plus_add_minus i j]
  simp only [evalPhi]
  push_cast
  ring

/-- Entry form of a product with the integer face normal. -/
theorem normalR_mul_entry (Q : Matrix (Fin 20) (Fin 20) Zphi) (i j : Fin 20) :
    (faceNormalR * Matrix.of fun i j ↦ evalPhi (Q i j) / 20) i j =
      evalPhi (∑ k : Fin 20, zint (faceNormalZ i k) (Q k j)) / 20 := by
  simp only [Matrix.mul_apply, faceNormalR, castZ_apply, Matrix.of_apply]
  rw [evalPhi_sum, Finset.sum_div]
  refine Finset.sum_congr rfl fun k _ ↦ ?_
  rw [evalPhi_zint]
  ring

/-- **(6)** `N goldenPlusR = (3 + √5) goldenPlusR`. -/
theorem normalR_plus : faceNormalR * goldenPlusR = (3 + Real.sqrt 5) • goldenPlusR := by
  ext i j
  unfold goldenPlusR
  rw [normalR_mul_entry, normal_plus i j, evalPhi_zmul, evalPhi_lamPlus]
  simp only [Matrix.smul_apply, Matrix.of_apply, smul_eq_mul]
  ring

/-- **(6)** `N goldenMinusR = (3 - √5) goldenMinusR`. -/
theorem normalR_minus : faceNormalR * goldenMinusR = (3 - Real.sqrt 5) • goldenMinusR := by
  ext i j
  unfold goldenMinusR
  rw [normalR_mul_entry, normal_minus i j, evalPhi_zmul, evalPhi_lamMinus]
  simp only [Matrix.smul_apply, Matrix.of_apply, smul_eq_mul]
  ring

/-- Entry form of a product of two `ℤ[φ]` tables. -/
theorem mulR_entry (Q Q' : Matrix (Fin 20) (Fin 20) Zphi) (i j : Fin 20) :
    ((Matrix.of fun i j ↦ evalPhi (Q i j) / 20) * Matrix.of fun i j ↦ evalPhi (Q' i j) / 20) i j =
      evalPhi (∑ k : Fin 20, zmul (Q i k) (Q' k j)) / 400 := by
  simp only [Matrix.mul_apply, Matrix.of_apply]
  rw [evalPhi_sum, Finset.sum_div]
  refine Finset.sum_congr rfl fun k _ ↦ ?_
  rw [evalPhi_zmul]
  ring

/-- **(6)** Idempotence and orthogonality of the real projectors. -/
theorem plusR_idem : goldenPlusR * goldenPlusR = goldenPlusR := by
  ext i j
  unfold goldenPlusR
  rw [mulR_entry, plus_idem i j, evalPhi_zint]
  simp only [Matrix.of_apply]
  push_cast
  ring

theorem minusR_idem : goldenMinusR * goldenMinusR = goldenMinusR := by
  ext i j
  unfold goldenMinusR
  rw [mulR_entry, minus_idem i j, evalPhi_zint]
  simp only [Matrix.of_apply]
  push_cast
  ring

theorem plusR_mul_minusR : goldenPlusR * goldenMinusR = 0 := by
  ext i j
  unfold goldenPlusR goldenMinusR
  rw [mulR_entry, plus_mul_minus i j, evalPhi_zero]
  simp

/-- **(6)** Traces `3`: each golden eigenspace is three-dimensional (rank
equals trace for a symmetric idempotent, an inference outside this file). -/
theorem plusR_trace : Matrix.trace goldenPlusR = 3 := by
  simp only [Matrix.trace, Matrix.diag, goldenPlusR, Matrix.of_apply]
  rw [← Finset.sum_div, ← evalPhi_sum, plus_trace]
  simp only [evalPhi]
  push_cast
  ring

theorem minusR_trace : Matrix.trace goldenMinusR = 3 := by
  simp only [Matrix.trace, Matrix.diag, goldenMinusR, Matrix.of_apply]
  rw [← Finset.sum_div, ← evalPhi_sum, minus_trace]
  simp only [evalPhi]
  push_cast
  ring

theorem plusR_symm : goldenPlusRᵀ = goldenPlusR := by
  ext i j
  simp only [Matrix.transpose_apply, goldenPlusR, Matrix.of_apply]
  rw [plus_symm i j]

theorem minusR_symm : goldenMinusRᵀ = goldenMinusR := by
  ext i j
  simp only [Matrix.transpose_apply, goldenMinusR, Matrix.of_apply]
  rw [minus_symm i j]

end

end OPH.GoldenSectorCharacters

#print axioms OPH.GoldenSectorCharacters.plus_add_minus
#print axioms OPH.GoldenSectorCharacters.plus_mul_minus
#print axioms OPH.GoldenSectorCharacters.plus_idem
#print axioms OPH.GoldenSectorCharacters.normal_plus
#print axioms OPH.GoldenSectorCharacters.normal_minus
#print axioms OPH.GoldenSectorCharacters.minus_eq_conj_plus
#print axioms OPH.GoldenSectorCharacters.five_plus_decomp
#print axioms OPH.GoldenSectorCharacters.plus_trace
#print axioms OPH.GoldenSectorCharacters.chi_plus_values
#print axioms OPH.GoldenSectorCharacters.chi_minus_values
#print axioms OPH.GoldenSectorCharacters.chi_conj
#print axioms OPH.GoldenSectorCharacters.chi_differ_iff
#print axioms OPH.GoldenSectorCharacters.chi_square
#print axioms OPH.GoldenSectorCharacters.five_class_counts
#print axioms OPH.GoldenSectorCharacters.chi_plus_class
#print axioms OPH.GoldenSectorCharacters.character_norms_golden
#print axioms OPH.GoldenSectorCharacters.plus_equivariant
#print axioms OPH.GoldenSectorCharacters.faceAct_comm_plus
#print axioms OPH.GoldenSectorCharacters.faceAct_comm_minus
#print axioms OPH.GoldenSectorCharacters.plusR_add_minusR
#print axioms OPH.GoldenSectorCharacters.normalR_plus
#print axioms OPH.GoldenSectorCharacters.normalR_minus
#print axioms OPH.GoldenSectorCharacters.plusR_trace
#print axioms OPH.GoldenSectorCharacters.plusR_idem
#print axioms OPH.GoldenSectorCharacters.plusR_mul_minusR
