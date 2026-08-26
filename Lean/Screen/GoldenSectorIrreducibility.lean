import GoldenSectorCharacters

open scoped BigOperators Matrix

namespace OPH.GoldenSectorIrreducibility

open OPH.A5PortAction
open OPH.LocalFaceMaxwellAction
open OPH.ScaledMaxwellStability
open OPH.CarrierModeEquivariance
open OPH.ScreenCarrierMapCandidate
open OPH.GoldenSectorCharacters

/-!
# Irreducibility of the two golden pieces by the Burnside span criterion

STATUS.  Kernel `decide` checks on committed finite tables over `ℤ[φ]`
(a left inverse `B` of the first three columns `A` of each golden
projector table, nine listed group elements, the nine transported
`3 × 3` action tables `M`, and an inverse certificate `N` with
`N · M = e · 1`), transported to exact real linear algebra through
`evalPhi` and to the language of submodules through `Matrix.toLin'`.
The carrier, its face orientation, the operator `N = C Cᵀ`, and the two
golden projector tables are declared in the imported files.  No register
row is discharged.

WHAT IS PROVED.  For each of the two certificates `plusCert`
(`goldenPlusZ`) and `minusCert` (`goldenMinusZ`), with
`V = Fin 20 → ℝ`, `P = projR C` the real projector (`goldenPlusR`,
`goldenMinusR`), `W = range (toLin' P)`, and `gR p = castZ (faceActZ p)`:
1. Invariance.  `W` is invariant under every listed row
   (`W_invariant`), from the commutation `gR p * P = P * gR p`
   (`gR_comm_projR`), itself the real transport of `faceAct_comm_plus`,
   `faceAct_comm_minus`.
2. Rank three.  `Module.finrank ℝ W = 3` (`finrank_W`).  Route: the
   real transports `AR` (`20 × 3`, the first three columns of the
   projector table) and `BR` (`3 × 20`, the listed table `B` over `20`)
   satisfy `BR * AR = 1` (`BR_mul_AR`) and `AR * BR = P`
   (`AR_mul_BR`), both kernel checks over `ℤ[φ]` (`BA`, `AB`), so `W`
   equals the range of the injective map `toLin' AR` from `Fin 3 → ℝ`
   (`W_eq_range_AR`).  The trace identity of Mathlib is not used; the
   factorization is explicit.
3. Burnside span.  `S = span {gR p * P : p ∈ perms}` inside
   `Matrix (Fin 20) (Fin 20) ℝ` has `finrank S = 9` (`finrank_S`).
   Upper bound: every element of `S` equals `P * X * P`
   (`sandwich_of_mem_S`), so the linear map `Phi X = BR * X * AR` into
   `Matrix (Fin 3) (Fin 3) ℝ` is injective on `S` (`Phi_injOn_S`,
   `finrank_S_le`).  Lower bound: the nine transported tables
   `gt b = BR * gR (elems b) * AR`, with `elems` the nine listed rows,
   are the real transports of the rows of `M` (`gt_eq`), and the
   certificate `N · M = e · 1`, `e ≠ 0`, writes every `3 × 3` real
   matrix as an explicit combination of them (`span_gt_eq_top`); hence
   `Phi` maps `S` onto the whole nine-dimensional target
   (`map_Phi_S_eq_top`, `finrank_S_ge`).  Method: the nine coordinate
   functionals are the entries of the transported table `BR * X * AR`
   rather than nine raw matrix positions; the determinant is replaced by
   the two-sided inverse certificate `N`, which avoids a `9 × 9`
   determinant in the kernel.  Surjectivity onto the endomorphisms of
   `W` in the basis `AR`: `Phi_S_surjective`.
4. Irreducibility.  For every submodule `U ≤ W` with
   `gR p *ᵥ u ∈ U` for all listed `p` and all `u ∈ U`, `U = ⊥ ∨ U = W`
   (`irreducible`; instances `goldenPlus_irreducible`,
   `goldenMinus_irreducible`).  Proof: the image `U' = U.map (toLin' BR)`
   in `Fin 3 → ℝ` is stable under every transported table, hence, by
   item 3, under every `3 × 3` real matrix (`stab`,
   `span_gt_le_stab`); a nonzero vector of `U'` is then carried to every
   vector, so `U' = ⊥` or `U' = ⊤`, and `U` is recovered from `U'` by
   `AR` since `AR * BR = P` fixes `W` pointwise (`projR_mulVec_of_mem_W`).
5. Both sectors.  The certificate for `goldenMinusZ` is the entrywise
   Galois conjugate of the certificate for `goldenPlusZ` (same nine
   rows, same column choice, conjugate tables, `eMinus = zconj ePlus`,
   `MMinusRows = zconj MPlusRows` entrywise, `conj_certificates`); the
   real theorems are instantiated separately from the same abstract
   development, so both sectors carry the theorem and no conjugate
   corollary is inferred.

Inference outside the theorems (item 5 of the lane).  `S` equal to the
full endomorphism algebra of `W` over `ℝ` persists under extension of
scalars to `ℂ` (a spanning set stays spanning), so each golden piece is
absolutely irreducible; together with the characters of
`GoldenSectorCharacters` (`chi_plus_values`, `chi_minus_values`) the
two pieces read as the two Galois-conjugate three-dimensional
irreducibles of the abstract `A5` character table.  Neither statement is
a theorem of this file.

PRIOR WORK.  `Screen/GoldenSectorCharacters.lean` lists the two
projector tables and proves their projector, eigenvalue, equivariance,
commutation (`faceAct_comm_plus`, `faceAct_comm_minus`), character, and
norm identities, with irreducibility stated there as an inference from
the character norm `60`; this file replaces that inference by a theorem.
`Screen/CarrierModeEquivariance.lean` lists `faceActZ`, `facePerm`,
`faceSign_eq_one`, and the integer characters.  `Screen/A5OPH.lean`
records the irreducible dimension multiset `{1, 3, 3, 5}` of the
abstract group as a subset-sum table (`SubmoduleDims`,
`three_appears_twice`) with no matrix representation attached.
`Screen/A5Commutant.lean` computes the four-dimensional commutant of
the port permutation module over `ℚ` (`orbitals_independent`).
`QFT/GaugeIrreducibleBorn.lean` proves a `2 × 2` complex Schur form
(`commutant_scalar`) for a different pair of generators.
`Geometry/ScreenCarrierMapCandidate.lean` supplies `Zphi`, `zmul`,
`evalPhi`, `evalPhi_ne_zero`.  No prior module states the Burnside span
dimension or the irreducibility of a golden piece as a theorem about
submodules of `Fin 20 → ℝ`.  The content of this file beyond the prior
modules is: the explicit rank-three factorization, the Burnside span
dimension `9`, the surjection onto the transported endomorphisms, and the
irreducibility theorem for both sectors.

ROWS TOUCHED (none discharged).  Source clock and duration row: untouched
(no step appears).  Physical spacetime attachment row: the listed group is
an incidence automorphism group of the declared carrier; the theorems
here classify its action on the two golden eigenspaces as irreducible,
and the identification of that group with a physical rotation group is
open.  Light-signal row: the identification of a golden mode with a
physical oscillation is open; the theorems here show that each golden
eigenspace carries no invariant proper subspace, so no further
symmetry-based splitting of a golden mode triple is available.
Coupled-action row: untouched.  Laboratory clock and energy calibration
import: untouched.  Gravitation-route energy identification: untouched.

NEGATIVES CITED.  Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`), at scope only:
realized histories select no velocity curvature or Legendre map, so the
operator `N` whose eigenspaces are classified here is part of a declared
evolution; the representation-theoretic statements concern the declared
carrier action.

CONVENTIONS.  `Zphi = ℤ × ℤ` with `(a, b) = a + bφ`, `φ² = φ + 1`;
`zint n` is integer scaling, `zconj (a, b) = (a + b, -b)`.  The first
three columns of a projector table are `colsZ Q i a = Q i a` for
`a : Fin 3`.  Rows of `perms` act on `Fin 20 → ℝ` by `gR p = castZ
(faceActZ p)`, with `(gR p *ᵥ v) i = v (facePerm p i)` on the listed rows
(`faceSign_eq_one`).  Matrix positions `(c, a) : Fin 3 × Fin 3` are
listed at index `3c + a` (`pairIdx`).  The transported action table of
row `b` is `gt b = BR * gR (elems b) * AR`, and `M b (c, a)` is its
`ℤ[φ]` entry (`M_spec` at scale `20`).

FALSIFIER.  A `B` table with `B · A ≠ 20 · 1` or `A · B ≠ Q`, a
transported table entry off the listed `M` row, a product `N · M` off
`e · 1`, or `e = 0` would make the corresponding kernel check fail; a
listed row with `gR p * P ≠ P * gR p` would break `W_invariant`.

Axiom audit.  The `#print axioms` lines at the end of the file show at
most `propext`, `Classical.choice`, and `Quot.sound`.
-/

/-! ## 1. Tables -/

/-- The first three columns of a `Fin 20 × Fin 20` table. -/
def colsZ (Q : Matrix (Fin 20) (Fin 20) Zphi) : Matrix (Fin 20) (Fin 3) Zphi :=
  Matrix.of fun i a ↦ Q i ⟨a.val, by omega⟩

/-- The index `3c + a` of a matrix position. -/
def pairIdx (q : Fin 3 × Fin 3) : Fin 9 := ⟨q.1.val * 3 + q.2.val, by omega⟩

/-- The nine listed rows of `perms` used in the span certificate (rows
`0, 1, 2, 5, 6, 7, 10, 11, 12` of the list). -/
def elems9 : Fin 9 → List Nat :=
  ![[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    [0, 2, 4, 1, 6, 8, 3, 5, 10, 7, 9, 11],
    [0, 3, 1, 6, 2, 7, 4, 9, 5, 10, 8, 11],
    [1, 0, 3, 2, 7, 6, 5, 4, 9, 8, 11, 10],
    [1, 2, 0, 5, 3, 4, 7, 8, 6, 11, 9, 10],
    [1, 3, 7, 0, 5, 9, 2, 6, 11, 4, 8, 10],
    [2, 0, 1, 4, 5, 3, 8, 6, 7, 10, 11, 9],
    [2, 1, 5, 0, 8, 7, 4, 3, 11, 6, 10, 9],
    [2, 4, 0, 8, 1, 6, 5, 10, 3, 11, 7, 9]]

/-- The listed rows indexed by matrix positions. -/
def elems (b : Fin 3 × Fin 3) : List Nat := elems9 (pairIdx b)

theorem elems_mem : ∀ b : Fin 3 × Fin 3, elems b ∈ perms := by decide

/-- Left inverse table for the `λ₊` sector: `BPlus · A = 20 · 1₃` and `A · BPlus = goldenPlusZ` with `A` the first three columns of `goldenPlusZ`. -/
def BPlus : Matrix (Fin 3) (Fin 20) Zphi :=
  Matrix.of
    ![![(1, 0), (0, 0), (0, 0), (-1, 1), (-1, 1), (1, -2), (0, -1), (2, 0), (0, -1), (2, 0), (-2, 0), (0, 1), (-2, 0), (0, 1), (1, -1), (1, -1), (-1, 2), (0, 0), (0, 0), (-1, 0)],
      ![(0, 0), (1, 0), (0, 0), (1, -1), (1, 0), (-1, 0), (0, -1), (0, 1), (-1, 0), (-1, 1), (1, -1), (1, 0), (0, -1), (0, 1), (-1, 0), (-1, 1), (1, 0), (0, 0), (-1, 0), (0, 0)],
      ![(0, 0), (0, 0), (1, 0), (1, 0), (1, -1), (-1, 0), (-1, 0), (-1, 1), (0, -1), (0, 1), (0, -1), (0, 1), (1, -1), (1, 0), (-1, 1), (-1, 0), (1, 0), (-1, 0), (0, 0), (0, 0)]]

/-- The nine transported action tables of the `λ₊` sector, one row per listed element of `elems9`, columns indexed by the nine matrix positions `(c, a) ↦ 3c + a`, at scale `1/20`. -/
def MPlusRows : Matrix (Fin 9) (Fin 9) Zphi :=
  Matrix.of
    ![![(1, 0), (0, 0), (0, 0), (0, 0), (1, 0), (0, 0), (0, 0), (0, 0), (1, 0)],
      ![(0, 0), (-1, 1), (1, 0), (1, 0), (1, -1), (0, 0), (0, 0), (1, 0), (0, 0)],
      ![(0, 0), (1, 0), (-1, 1), (0, 0), (0, 0), (1, 0), (1, 0), (0, 0), (1, -1)],
      ![(0, 0), (1, 0), (0, -1), (1, 0), (0, 0), (0, -1), (0, 0), (0, 0), (-1, 0)],
      ![(1, 0), (0, 0), (1, -2), (0, 0), (0, 0), (-1, 0), (0, 0), (1, 0), (-1, 0)],
      ![(-1, 1), (0, 0), (-2, 0), (1, -1), (1, 0), (1, -1), (1, 0), (0, 0), (0, -1)],
      ![(1, 0), (1, -2), (0, 0), (0, 0), (-1, 0), (1, 0), (0, 0), (-1, 0), (0, 0)],
      ![(0, 0), (0, -1), (-1, 1), (1, 0), (0, -1), (1, -1), (0, 0), (-1, 0), (1, 0)],
      ![(0, 0), (0, -1), (1, 0), (0, 0), (-1, 0), (0, 0), (1, 0), (0, -1), (0, 0)]]

/-- Certificate of invertibility of `MPlusRows`: `NPlusRows · MPlusRows = ePlus · 1₉`. -/
def NPlusRows : Matrix (Fin 9) (Fin 9) Zphi :=
  Matrix.of
    ![![(-2, -4), (0, 0), (0, 0), (-1, -3), (-4, -5), (4, 7), (-1, -1), (4, 7), (-4, -7)],
      ![(-1, -1), (2, 3), (-2, -3), (-2, -2), (-2, -4), (4, 6), (1, 1), (2, 3), (-2, -3)],
      ![(-1, -1), (-2, -3), (2, 3), (-1, -1), (0, -1), (2, 3), (0, 0), (4, 6), (-4, -6)],
      ![(0, 0), (-2, -4), (0, 0), (-1, -1), (4, 7), (-4, -7), (-1, -3), (-4, -5), (4, 7)],
      ![(-1, -1), (-2, -3), (2, 3), (0, 0), (4, 6), (-4, -6), (-1, -1), (0, -1), (2, 3)],
      ![(1, 1), (0, 1), (-2, -3), (1, 1), (2, 3), (-2, -3), (-2, -2), (-2, -4), (4, 6)],
      ![(0, 0), (0, 0), (-2, -4), (3, 3), (0, 1), (-4, -5), (1, 3), (-4, -7), (2, 3)],
      ![(1, 1), (-2, -3), (0, 1), (2, 4), (0, 0), (-2, -4), (1, 1), (-2, -3), (2, 3)],
      ![(-1, -1), (2, 3), (-2, -3), (1, 3), (0, -1), (0, -1), (2, 2), (-4, -6), (2, 4)]]

/-- `ePlus = (-4, -6)`, the scalar of the certificate. -/
def ePlus : Zphi := (-4, -6)

/-- Left inverse table for the `λ₋` sector: `BMinus · A = 20 · 1₃` and `A · BMinus = goldenMinusZ` with `A` the first three columns of `goldenMinusZ`. -/
def BMinus : Matrix (Fin 3) (Fin 20) Zphi :=
  Matrix.of
    ![![(1, 0), (0, 0), (0, 0), (0, -1), (0, -1), (-1, 2), (-1, 1), (2, 0), (-1, 1), (2, 0), (-2, 0), (1, -1), (-2, 0), (1, -1), (0, 1), (0, 1), (1, -2), (0, 0), (0, 0), (-1, 0)],
      ![(0, 0), (1, 0), (0, 0), (0, 1), (1, 0), (-1, 0), (-1, 1), (1, -1), (-1, 0), (0, -1), (0, 1), (1, 0), (-1, 1), (1, -1), (-1, 0), (0, -1), (1, 0), (0, 0), (-1, 0), (0, 0)],
      ![(0, 0), (0, 0), (1, 0), (1, 0), (0, 1), (-1, 0), (-1, 0), (0, -1), (-1, 1), (1, -1), (-1, 1), (1, -1), (0, 1), (1, 0), (0, -1), (-1, 0), (1, 0), (-1, 0), (0, 0), (0, 0)]]

/-- The nine transported action tables of the `λ₋` sector, one row per listed element of `elems9`, columns indexed by the nine matrix positions `(c, a) ↦ 3c + a`, at scale `1/20`. -/
def MMinusRows : Matrix (Fin 9) (Fin 9) Zphi :=
  Matrix.of
    ![![(1, 0), (0, 0), (0, 0), (0, 0), (1, 0), (0, 0), (0, 0), (0, 0), (1, 0)],
      ![(0, 0), (0, -1), (1, 0), (1, 0), (0, 1), (0, 0), (0, 0), (1, 0), (0, 0)],
      ![(0, 0), (1, 0), (0, -1), (0, 0), (0, 0), (1, 0), (1, 0), (0, 0), (0, 1)],
      ![(0, 0), (1, 0), (-1, 1), (1, 0), (0, 0), (-1, 1), (0, 0), (0, 0), (-1, 0)],
      ![(1, 0), (0, 0), (-1, 2), (0, 0), (0, 0), (-1, 0), (0, 0), (1, 0), (-1, 0)],
      ![(0, -1), (0, 0), (-2, 0), (0, 1), (1, 0), (0, 1), (1, 0), (0, 0), (-1, 1)],
      ![(1, 0), (-1, 2), (0, 0), (0, 0), (-1, 0), (1, 0), (0, 0), (-1, 0), (0, 0)],
      ![(0, 0), (-1, 1), (0, -1), (1, 0), (-1, 1), (0, 1), (0, 0), (-1, 0), (1, 0)],
      ![(0, 0), (-1, 1), (1, 0), (0, 0), (-1, 0), (0, 0), (1, 0), (-1, 1), (0, 0)]]

/-- Certificate of invertibility of `MMinusRows`: `NMinusRows · MMinusRows = eMinus · 1₉`. -/
def NMinusRows : Matrix (Fin 9) (Fin 9) Zphi :=
  Matrix.of
    ![![(-6, 4), (0, 0), (0, 0), (-4, 3), (-9, 5), (11, -7), (-2, 1), (11, -7), (-11, 7)],
      ![(-2, 1), (5, -3), (-5, 3), (-4, 2), (-6, 4), (10, -6), (2, -1), (5, -3), (-5, 3)],
      ![(-2, 1), (-5, 3), (5, -3), (-2, 1), (-1, 1), (5, -3), (0, 0), (10, -6), (-10, 6)],
      ![(0, 0), (-6, 4), (0, 0), (-2, 1), (11, -7), (-11, 7), (-4, 3), (-9, 5), (11, -7)],
      ![(-2, 1), (-5, 3), (5, -3), (0, 0), (10, -6), (-10, 6), (-2, 1), (-1, 1), (5, -3)],
      ![(2, -1), (1, -1), (-5, 3), (2, -1), (5, -3), (-5, 3), (-4, 2), (-6, 4), (10, -6)],
      ![(0, 0), (0, 0), (-6, 4), (6, -3), (1, -1), (-9, 5), (4, -3), (-11, 7), (5, -3)],
      ![(2, -1), (-5, 3), (1, -1), (6, -4), (0, 0), (-6, 4), (2, -1), (-5, 3), (5, -3)],
      ![(-2, 1), (5, -3), (-5, 3), (4, -3), (-1, 1), (-1, 1), (4, -2), (-10, 6), (6, -4)]]

/-- `eMinus = (-10, 6)`, the scalar of the certificate. -/
def eMinus : Zphi := (-10, 6)

/-- The transported action table of the listed row `elems b` against
`B` and the first three columns of `Q`, at scale `20`: entry `(c, a)` is
`∑ k, B c k · Q (facePerm (elems b) k) a`. -/
def tableM (Q : Matrix (Fin 20) (Fin 20) Zphi) (B : Matrix (Fin 3) (Fin 20) Zphi)
    (b q : Fin 3 × Fin 3) : Zphi :=
  ∑ k : Fin 20, zmul (B q.1 k) (colsZ Q (facePerm (elems b) k) q.2)

/-- A golden sector certificate: the projector table, its commutation
with every listed row, a left inverse of its first three columns, the
nine transported action tables, and their inverse certificate. -/
structure GoldenCert where
  Q : Matrix (Fin 20) (Fin 20) Zphi
  comm : ∀ p ∈ perms, ∀ i j : Fin 20,
    (∑ k : Fin 20, zint (faceActZ p i k) (Q k j)) = ∑ k : Fin 20, zint (faceActZ p k j) (Q i k)
  B : Matrix (Fin 3) (Fin 20) Zphi
  BA : ∀ c a : Fin 3, (∑ k : Fin 20, zmul (B c k) (colsZ Q k a)) = if c = a then (20, 0) else 0
  AB : ∀ i j : Fin 20, (∑ a : Fin 3, zmul (colsZ Q i a) (B a j)) = Q i j
  M : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) Zphi
  M_spec : ∀ b q : Fin 3 × Fin 3, tableM Q B b q = zint 20 (M b q)
  N : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) Zphi
  e : Zphi
  e_ne : e ≠ 0
  NM : ∀ q q' : Fin 3 × Fin 3, (∑ b : Fin 3 × Fin 3, zmul (N q b) (M b q')) =
    if q = q' then e else 0

/-- The `λ₊` tables indexed by positions. -/
def MPlus : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) Zphi :=
  Matrix.of fun b q ↦ MPlusRows (pairIdx b) (pairIdx q)
def NPlus : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) Zphi :=
  Matrix.of fun b q ↦ NPlusRows (pairIdx b) (pairIdx q)
def MMinus : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) Zphi :=
  Matrix.of fun b q ↦ MMinusRows (pairIdx b) (pairIdx q)
def NMinus : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) Zphi :=
  Matrix.of fun b q ↦ NMinusRows (pairIdx b) (pairIdx q)

theorem plus_BA : ∀ c a : Fin 3,
    (∑ k : Fin 20, zmul (BPlus c k) (colsZ goldenPlusZ k a)) = if c = a then (20, 0) else 0 := by
  decide

set_option maxRecDepth 16384 in
theorem plus_AB : ∀ i j : Fin 20,
    (∑ a : Fin 3, zmul (colsZ goldenPlusZ i a) (BPlus a j)) = goldenPlusZ i j := by
  decide

theorem plus_M_spec : ∀ b q : Fin 3 × Fin 3,
    tableM goldenPlusZ BPlus b q = zint 20 (MPlus b q) := by
  decide +kernel

theorem plus_NM : ∀ q q' : Fin 3 × Fin 3,
    (∑ b : Fin 3 × Fin 3, zmul (NPlus q b) (MPlus b q')) = if q = q' then ePlus else 0 := by
  decide

theorem minus_BA : ∀ c a : Fin 3,
    (∑ k : Fin 20, zmul (BMinus c k) (colsZ goldenMinusZ k a)) = if c = a then (20, 0) else 0 := by
  decide

set_option maxRecDepth 16384 in
theorem minus_AB : ∀ i j : Fin 20,
    (∑ a : Fin 3, zmul (colsZ goldenMinusZ i a) (BMinus a j)) = goldenMinusZ i j := by
  decide

theorem minus_M_spec : ∀ b q : Fin 3 × Fin 3,
    tableM goldenMinusZ BMinus b q = zint 20 (MMinus b q) := by
  decide +kernel

theorem minus_NM : ∀ q q' : Fin 3 × Fin 3,
    (∑ b : Fin 3 × Fin 3, zmul (NMinus q b) (MMinus b q')) = if q = q' then eMinus else 0 := by
  decide

/-- **(5)** The minus certificate is the entrywise Galois conjugate of the
plus certificate. -/
theorem conj_certificates :
    (∀ c k, BMinus c k = zconj (BPlus c k)) ∧ (∀ b q, MMinusRows b q = zconj (MPlusRows b q)) ∧
    (∀ b q, NMinusRows b q = zconj (NPlusRows b q)) ∧ eMinus = zconj ePlus := by
  decide

/-- The `λ₊` certificate. -/
def plusCert : GoldenCert where
  Q := goldenPlusZ
  comm := fun p hp i j ↦ faceAct_comm_plus p hp i j
  B := BPlus
  BA := plus_BA
  AB := plus_AB
  M := MPlus
  M_spec := plus_M_spec
  N := NPlus
  e := ePlus
  e_ne := by decide
  NM := plus_NM

/-- The `λ₋` certificate. -/
def minusCert : GoldenCert where
  Q := goldenMinusZ
  comm := fun p hp i j ↦ faceAct_comm_minus p hp i j
  B := BMinus
  BA := minus_BA
  AB := minus_AB
  M := MMinus
  M_spec := minus_M_spec
  N := NMinus
  e := eMinus
  e_ne := by decide
  NM := minus_NM

/-! ## 2. Real transport -/

noncomputable section

/-- The real projector of a certificate, `evalPhi Q / 20`. -/
def projR (C : GoldenCert) : Matrix (Fin 20) (Fin 20) ℝ :=
  Matrix.of fun i j ↦ evalPhi (C.Q i j) / 20

theorem projR_plus : projR plusCert = goldenPlusR := rfl
theorem projR_minus : projR minusCert = goldenMinusR := rfl

/-- The real transport of the first three projector columns. -/
def AR (C : GoldenCert) : Matrix (Fin 20) (Fin 3) ℝ :=
  Matrix.of fun i a ↦ evalPhi (colsZ C.Q i a)

/-- The real transport of the left inverse table, at scale `1/20`. -/
def BR (C : GoldenCert) : Matrix (Fin 3) (Fin 20) ℝ :=
  Matrix.of fun c k ↦ evalPhi (C.B c k) / 20

/-- The real action matrix of a listed row. -/
def gR (p : List Nat) : Matrix (Fin 20) (Fin 20) ℝ := castZ (faceActZ p)

/-- **(2)** `BR * AR = 1`. -/
theorem BR_mul_AR (C : GoldenCert) : BR C * AR C = 1 := by
  ext c a
  simp only [Matrix.mul_apply, BR, AR, Matrix.of_apply, Matrix.one_apply]
  have h : (∑ k : Fin 20, evalPhi (C.B c k) / 20 * evalPhi (colsZ C.Q k a)) =
      evalPhi (∑ k : Fin 20, zmul (C.B c k) (colsZ C.Q k a)) / 20 := by
    rw [evalPhi_sum, Finset.sum_div]
    refine Finset.sum_congr rfl fun k _ ↦ ?_
    rw [evalPhi_zmul]; ring
  rw [h, C.BA c a]
  split_ifs
  · simp [evalPhi]
  · simp [evalPhi_zero]

/-- **(2)** `AR * BR = P`. -/
theorem AR_mul_BR (C : GoldenCert) : AR C * BR C = projR C := by
  ext i j
  simp only [Matrix.mul_apply, BR, AR, projR, Matrix.of_apply]
  rw [← C.AB i j, evalPhi_sum, Finset.sum_div]
  refine Finset.sum_congr rfl fun a _ ↦ ?_
  rw [evalPhi_zmul]; ring

theorem projR_idem (C : GoldenCert) : projR C * projR C = projR C := by
  rw [← AR_mul_BR, Matrix.mul_assoc, ← Matrix.mul_assoc (BR C), BR_mul_AR, Matrix.one_mul]

theorem projR_mul_AR (C : GoldenCert) : projR C * AR C = AR C := by
  rw [← AR_mul_BR, Matrix.mul_assoc, BR_mul_AR, Matrix.mul_one]

theorem BR_mul_projR (C : GoldenCert) : BR C * projR C = BR C := by
  rw [← AR_mul_BR, ← Matrix.mul_assoc, BR_mul_AR, Matrix.one_mul]

/-- **(1)** Real commutation of every listed row with the projector. -/
theorem gR_comm_projR (C : GoldenCert) (p : List Nat) (hp : p ∈ perms) :
    gR p * projR C = projR C * gR p := by
  ext i j
  simp only [Matrix.mul_apply, gR, projR, castZ_apply, Matrix.of_apply]
  have h1 : (∑ k : Fin 20, ((faceActZ p i k : ℤ) : ℝ) * (evalPhi (C.Q k j) / 20)) =
      evalPhi (∑ k : Fin 20, zint (faceActZ p i k) (C.Q k j)) / 20 := by
    rw [evalPhi_sum, Finset.sum_div]
    refine Finset.sum_congr rfl fun k _ ↦ ?_
    rw [evalPhi_zint]; ring
  have h2 : (∑ k : Fin 20, evalPhi (C.Q i k) / 20 * ((faceActZ p k j : ℤ) : ℝ)) =
      evalPhi (∑ k : Fin 20, zint (faceActZ p k j) (C.Q i k)) / 20 := by
    rw [evalPhi_sum, Finset.sum_div]
    refine Finset.sum_congr rfl fun k _ ↦ ?_
    rw [evalPhi_zint]; ring
  rw [h1, h2, C.comm p hp i j]

/-- The action of a listed row on the transported columns is a row
permutation. -/
theorem gR_mul_AR (C : GoldenCert) (p : List Nat) (hp : p ∈ perms) :
    gR p * AR C = Matrix.of fun i a ↦ evalPhi (colsZ C.Q (facePerm p i) a) := by
  ext i a
  simp only [Matrix.mul_apply, gR, AR, castZ_apply, faceActZ, Matrix.of_apply]
  simp only [Int.cast_ite, Int.cast_zero, ite_mul, zero_mul, Finset.sum_ite_eq, Finset.mem_univ,
    if_true, faceSign_eq_one p hp, Int.cast_one, one_mul]

/-- The transported action table of a listed position row. -/
def gt (C : GoldenCert) (b : Fin 3 × Fin 3) : Matrix (Fin 3) (Fin 3) ℝ :=
  BR C * gR (elems b) * AR C

/-- **(3)** The transported tables are the real transports of `M`. -/
theorem gt_eq (C : GoldenCert) (b : Fin 3 × Fin 3) :
    gt C b = Matrix.of fun c a ↦ evalPhi (C.M b (c, a)) := by
  ext c a
  unfold gt
  rw [Matrix.mul_assoc, gR_mul_AR C _ (elems_mem b)]
  simp only [Matrix.mul_apply, BR, Matrix.of_apply]
  have h : (∑ k : Fin 20, evalPhi (C.B c k) / 20 * evalPhi (colsZ C.Q (facePerm (elems b) k) a)) =
      evalPhi (tableM C.Q C.B b (c, a)) / 20 := by
    unfold tableM
    rw [evalPhi_sum, Finset.sum_div]
    refine Finset.sum_congr rfl fun k _ ↦ ?_
    rw [evalPhi_zmul]; ring
  rw [h, C.M_spec b (c, a), evalPhi_zint]
  push_cast
  ring

/-! ## 3. The golden piece `W` and its dimension -/

/-- The golden piece: the range of the real projector. -/
def W (C : GoldenCert) : Submodule ℝ (Fin 20 → ℝ) := LinearMap.range (Matrix.toLin' (projR C))

theorem mem_W_iff (C : GoldenCert) (w : Fin 20 → ℝ) : w ∈ W C ↔ (projR C).mulVec w = w := by
  constructor
  · rintro ⟨v, rfl⟩
    rw [Matrix.toLin'_apply, Matrix.mulVec_mulVec, projR_idem]
  · intro h
    exact ⟨w, by rw [Matrix.toLin'_apply, h]⟩

/-- The projector fixes `W` pointwise. -/
theorem projR_mulVec_of_mem_W (C : GoldenCert) {w : Fin 20 → ℝ} (hw : w ∈ W C) :
    (projR C).mulVec w = w := (mem_W_iff C w).1 hw

/-- **(2)** `W` is the range of the injective map `toLin' AR`. -/
theorem W_eq_range_AR (C : GoldenCert) : W C = LinearMap.range (Matrix.toLin' (AR C)) := by
  apply le_antisymm
  · rintro _ ⟨v, rfl⟩
    refine ⟨(BR C).mulVec v, ?_⟩
    rw [Matrix.toLin'_apply, Matrix.toLin'_apply, Matrix.mulVec_mulVec, AR_mul_BR]
  · rintro _ ⟨v, rfl⟩
    refine ⟨(AR C).mulVec v, ?_⟩
    rw [Matrix.toLin'_apply, Matrix.toLin'_apply, Matrix.mulVec_mulVec, projR_mul_AR]

theorem toLin'_AR_injective (C : GoldenCert) : Function.Injective (Matrix.toLin' (AR C)) := by
  intro v v' h
  have h' := congrArg (fun x ↦ (BR C).mulVec x) h
  simp only [Matrix.toLin'_apply, Matrix.mulVec_mulVec, BR_mul_AR, Matrix.one_mulVec] at h'
  exact h'

/-- **(2)** `Module.finrank ℝ W = 3`. -/
theorem finrank_W (C : GoldenCert) : Module.finrank ℝ (W C) = 3 := by
  rw [W_eq_range_AR, LinearMap.finrank_range_of_inj (toLin'_AR_injective C), Module.finrank_fin_fun]

/-- **(1)** `W` is invariant under every listed row. -/
theorem W_invariant (C : GoldenCert) (p : List Nat) (hp : p ∈ perms) {w : Fin 20 → ℝ}
    (hw : w ∈ W C) : (gR p).mulVec w ∈ W C := by
  rw [mem_W_iff, Matrix.mulVec_mulVec, ← gR_comm_projR C p hp, ← Matrix.mulVec_mulVec,
    projR_mulVec_of_mem_W C hw]

/-! ## 4. The Burnside span -/

/-- The generator `gR p * P` of the Burnside span. -/
def gP (C : GoldenCert) (p : List Nat) : Matrix (Fin 20) (Fin 20) ℝ := gR p * projR C

/-- The Burnside span: the real span of the sixty matrices `gR p * P`. -/
def S (C : GoldenCert) : Submodule ℝ (Matrix (Fin 20) (Fin 20) ℝ) :=
  Submodule.span ℝ {X | ∃ p ∈ perms, X = gP C p}

/-- The transport of a `20 × 20` matrix to the basis `AR` of `W`. -/
def Phi (C : GoldenCert) : Matrix (Fin 20) (Fin 20) ℝ →ₗ[ℝ] Matrix (Fin 3) (Fin 3) ℝ where
  toFun X := BR C * X * AR C
  map_add' X Y := by simp only [Matrix.mul_add, Matrix.add_mul]
  map_smul' c X := by simp only [Matrix.mul_smul, Matrix.smul_mul, RingHom.id_apply]

theorem Phi_apply (C : GoldenCert) (X : Matrix (Fin 20) (Fin 20) ℝ) :
    Phi C X = BR C * X * AR C := rfl

/-- **(3)** Every element of the Burnside span is sandwiched by `P`. -/
theorem sandwich_of_mem_S (C : GoldenCert) {X : Matrix (Fin 20) (Fin 20) ℝ} (hX : X ∈ S C) :
    projR C * X * projR C = X := by
  induction hX using Submodule.span_induction with
  | mem X hX =>
    obtain ⟨p, hp, rfl⟩ := hX
    unfold gP
    rw [← Matrix.mul_assoc, ← gR_comm_projR C p hp, Matrix.mul_assoc, projR_idem,
      Matrix.mul_assoc, projR_idem]
  | zero => simp
  | add X Y _ _ hX hY => rw [Matrix.mul_add, Matrix.add_mul, hX, hY]
  | smul c X _ hX => rw [Matrix.mul_smul, Matrix.smul_mul, hX]

/-- **(3)** `Phi` is injective on the Burnside span. -/
theorem Phi_injOn_S (C : GoldenCert) {X : Matrix (Fin 20) (Fin 20) ℝ} (hX : X ∈ S C)
    (h : Phi C X = 0) : X = 0 := by
  have h2 : AR C * Phi C X * BR C = 0 := by rw [h]; simp
  rw [Phi_apply, ← Matrix.mul_assoc, ← Matrix.mul_assoc, AR_mul_BR, Matrix.mul_assoc,
    Matrix.mul_assoc, AR_mul_BR, ← Matrix.mul_assoc, sandwich_of_mem_S C hX] at h2
  exact h2

theorem finrank_matrix_three : Module.finrank ℝ (Matrix (Fin 3) (Fin 3) ℝ) = 9 := by
  rw [Module.finrank_matrix]
  simp

/-- **(3)** `finrank S ≤ 9`. -/
theorem finrank_S_le (C : GoldenCert) : Module.finrank ℝ (S C) ≤ 9 := by
  have hinj : Function.Injective ((Phi C).domRestrict (S C)) := by
    rw [← LinearMap.ker_eq_bot, LinearMap.ker_eq_bot']
    intro x hx
    apply Subtype.ext
    exact Phi_injOn_S C x.2 hx
  rw [← finrank_matrix_three]
  exact LinearMap.finrank_le_finrank_of_injective hinj

/-- The real transport of the inverse certificate. -/
def NR (C : GoldenCert) (q b : Fin 3 × Fin 3) : ℝ := evalPhi (C.N q b)

theorem eR_ne (C : GoldenCert) : evalPhi C.e ≠ 0 := evalPhi_ne_zero C.e_ne

/-- Real form of the certificate `N · M = e · 1`. -/
theorem NR_mul_M (C : GoldenCert) (q q' : Fin 3 × Fin 3) :
    (∑ b : Fin 3 × Fin 3, NR C q b * evalPhi (C.M b q')) = if q = q' then evalPhi C.e else 0 := by
  unfold NR
  have h : (∑ b : Fin 3 × Fin 3, evalPhi (C.N q b) * evalPhi (C.M b q')) =
      evalPhi (∑ b : Fin 3 × Fin 3, zmul (C.N q b) (C.M b q')) := by
    rw [evalPhi_sum]
    exact Finset.sum_congr rfl fun b _ ↦ (evalPhi_zmul _ _).symm
  rw [h, C.NM q q']
  split_ifs <;> simp [evalPhi_zero]

/-- **(3)** Every `3 × 3` real matrix is an explicit combination of the
nine transported tables. -/
theorem eq_sum_gt (C : GoldenCert) (Y : Matrix (Fin 3) (Fin 3) ℝ) :
    Y = ∑ b : Fin 3 × Fin 3,
      ((∑ q : Fin 3 × Fin 3, Y q.1 q.2 * NR C q b) / evalPhi C.e) • gt C b := by
  ext c a
  simp only [Matrix.sum_apply, Matrix.smul_apply, gt_eq, Matrix.of_apply, smul_eq_mul]
  have h : (∑ b : Fin 3 × Fin 3, (∑ q : Fin 3 × Fin 3, Y q.1 q.2 * NR C q b) / evalPhi C.e *
      evalPhi (C.M b (c, a))) =
      (∑ q : Fin 3 × Fin 3, Y q.1 q.2 * ∑ b : Fin 3 × Fin 3, NR C q b * evalPhi (C.M b (c, a))) /
        evalPhi C.e := by
    simp only [Finset.sum_div, div_mul_eq_mul_div, Finset.sum_mul, Finset.mul_sum]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun q _ ↦ Finset.sum_congr rfl fun b _ ↦ ?_
    ring
  rw [h]
  simp only [NR_mul_M, mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true]
  rw [mul_div_assoc, div_self (eR_ne C), mul_one]

/-- **(3)** The nine transported tables span all `3 × 3` real matrices. -/
theorem span_gt_eq_top (C : GoldenCert) : Submodule.span ℝ (Set.range (gt C)) = ⊤ := by
  rw [eq_top_iff]
  intro Y _
  rw [eq_sum_gt C Y]
  exact Submodule.sum_mem _ fun b _ ↦ Submodule.smul_mem _ _ (Submodule.subset_span ⟨b, rfl⟩)

theorem Phi_gP (C : GoldenCert) (p : List Nat) (_hp : p ∈ perms) :
    Phi C (gP C p) = BR C * gR p * AR C := by
  rw [Phi_apply, gP, Matrix.mul_assoc, Matrix.mul_assoc, projR_mul_AR, ← Matrix.mul_assoc]

theorem gP_mem_S (C : GoldenCert) (p : List Nat) (hp : p ∈ perms) : gP C p ∈ S C :=
  Submodule.subset_span ⟨p, hp, rfl⟩

/-- **(3)** `Phi` maps the Burnside span onto all `3 × 3` real matrices. -/
theorem map_Phi_S_eq_top (C : GoldenCert) : (S C).map (Phi C) = ⊤ := by
  rw [eq_top_iff, ← span_gt_eq_top, Submodule.span_le]
  rintro _ ⟨b, rfl⟩
  exact ⟨gP C (elems b), gP_mem_S C _ (elems_mem b), by
    rw [Phi_gP C _ (elems_mem b)]; rfl⟩

/-- **(3)** Surjectivity onto the transported endomorphisms of `W`. -/
theorem Phi_S_surjective (C : GoldenCert) (Y : Matrix (Fin 3) (Fin 3) ℝ) :
    ∃ X ∈ S C, Phi C X = Y := by
  have h : Y ∈ (S C).map (Phi C) := by rw [map_Phi_S_eq_top]; trivial
  obtain ⟨X, hX, rfl⟩ := h
  exact ⟨X, hX, rfl⟩

/-- **(3)** `9 ≤ finrank S`. -/
theorem finrank_S_ge (C : GoldenCert) : 9 ≤ Module.finrank ℝ (S C) := by
  have h := Submodule.finrank_map_le (Phi C) (S C)
  rw [map_Phi_S_eq_top, finrank_top, finrank_matrix_three] at h
  exact h

/-- **(3)** `finrank S = 9`: the Burnside span of the golden piece is the
full endomorphism algebra of the three-dimensional piece. -/
theorem finrank_S (C : GoldenCert) : Module.finrank ℝ (S C) = 9 :=
  le_antisymm (finrank_S_le C) (finrank_S_ge C)

/-! ## 5. Irreducibility -/

/-- The matrices preserving a submodule of `Fin 3 → ℝ`. -/
def stab (U' : Submodule ℝ (Fin 3 → ℝ)) : Submodule ℝ (Matrix (Fin 3) (Fin 3) ℝ) where
  carrier := {Y | ∀ u ∈ U', Y.mulVec u ∈ U'}
  add_mem' := by
    intro Y Z hY hZ u hu
    rw [Matrix.add_mulVec]
    exact U'.add_mem (hY u hu) (hZ u hu)
  zero_mem' := by
    intro u hu
    rw [Matrix.zero_mulVec]
    exact U'.zero_mem
  smul_mem' := by
    intro c Y hY u hu
    rw [Matrix.smul_mulVec]
    exact U'.smul_mem c (hY u hu)

theorem mem_stab {U' : Submodule ℝ (Fin 3 → ℝ)} {Y : Matrix (Fin 3) (Fin 3) ℝ} :
    Y ∈ stab U' ↔ ∀ u ∈ U', Y.mulVec u ∈ U' := Iff.rfl

/-- The transported image of an invariant submodule is stable under every
transported table. -/
theorem span_gt_le_stab (C : GoldenCert) (U : Submodule ℝ (Fin 20 → ℝ)) (hU : U ≤ W C)
    (hinv : ∀ p ∈ perms, ∀ u ∈ U, (gR p).mulVec u ∈ U) :
    Submodule.span ℝ (Set.range (gt C)) ≤ stab (U.map (Matrix.toLin' (BR C))) := by
  rw [Submodule.span_le]
  rintro _ ⟨b, rfl⟩ _ ⟨u, hu, rfl⟩
  refine ⟨(gR (elems b)).mulVec u, hinv _ (elems_mem b) u hu, ?_⟩
  simp only [Matrix.toLin'_apply, gt, Matrix.mulVec_mulVec]
  rw [Matrix.mul_assoc, AR_mul_BR, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec,
    projR_mulVec_of_mem_W C (hU hu), Matrix.mulVec_mulVec]

/-- A submodule of `Fin 3 → ℝ` stable under every matrix is trivial. -/
theorem trivial_of_stab_top (U' : Submodule ℝ (Fin 3 → ℝ)) (h : stab U' = ⊤) :
    U' = ⊥ ∨ U' = ⊤ := by
  by_cases hbot : U' = ⊥
  · exact Or.inl hbot
  right
  obtain ⟨u, hu, hne⟩ := Submodule.exists_mem_ne_zero_of_ne_bot hbot
  rw [eq_top_iff]
  intro v _
  have hs : dotProduct u u ≠ 0 := by
    rw [Ne, dotProduct_self_eq_zero]; exact hne
  let Y : Matrix (Fin 3) (Fin 3) ℝ := Matrix.of fun i j ↦ v i * u j / dotProduct u u
  have hY : Y ∈ stab U' := by rw [h]; trivial
  have hYu : Y.mulVec u = v := by
    funext i
    have hs' : (∑ j, u j * u j) ≠ 0 := hs
    simp only [Matrix.mulVec, dotProduct, Y, Matrix.of_apply]
    rw [show (∑ j, v i * u j / (∑ j, u j * u j) * u j) = v i * (∑ j, u j * u j) / ∑ j, u j * u j by
      rw [Finset.mul_sum, Finset.sum_div]
      exact Finset.sum_congr rfl fun j _ ↦ by ring]
    rw [mul_div_assoc, div_self hs', mul_one]
  rw [← hYu]
  exact hY u hu

/-- **(4)** Irreducibility: a submodule of the golden piece invariant under
every listed row is `⊥` or the whole piece. -/
theorem irreducible (C : GoldenCert) (U : Submodule ℝ (Fin 20 → ℝ)) (hU : U ≤ W C)
    (hinv : ∀ p ∈ perms, ∀ u ∈ U, (gR p).mulVec u ∈ U) : U = ⊥ ∨ U = W C := by
  have hstab : stab (U.map (Matrix.toLin' (BR C))) = ⊤ := by
    rw [eq_top_iff, ← span_gt_eq_top C]
    exact span_gt_le_stab C U hU hinv
  rcases trivial_of_stab_top _ hstab with hb | ht
  · left
    rw [eq_bot_iff]
    intro u hu
    have hBu : (BR C).mulVec u = 0 := by
      have : Matrix.toLin' (BR C) u ∈ U.map (Matrix.toLin' (BR C)) := ⟨u, hu, rfl⟩
      rw [hb, Submodule.mem_bot] at this
      exact this
    have : u = 0 := by
      rw [← projR_mulVec_of_mem_W C (hU hu), ← AR_mul_BR, ← Matrix.mulVec_mulVec, hBu,
        Matrix.mulVec_zero]
    rw [this]; exact Submodule.zero_mem _
  · right
    refine le_antisymm hU ?_
    intro w hw
    have : Matrix.toLin' (BR C) w ∈ U.map (Matrix.toLin' (BR C)) := by rw [ht]; trivial
    obtain ⟨u, hu, huw⟩ := this
    have hw' : w = u := by
      rw [← projR_mulVec_of_mem_W C hw, ← projR_mulVec_of_mem_W C (hU hu), ← AR_mul_BR,
        ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec]
      simp only [Matrix.toLin'_apply] at huw
      rw [huw]
    rw [hw']; exact hu

/-! ## 6. The two golden pieces -/

/-- The `λ₊` piece: the range of `goldenPlusR`. -/
def WPlus : Submodule ℝ (Fin 20 → ℝ) := LinearMap.range (Matrix.toLin' goldenPlusR)
/-- The `λ₋` piece: the range of `goldenMinusR`. -/
def WMinus : Submodule ℝ (Fin 20 → ℝ) := LinearMap.range (Matrix.toLin' goldenMinusR)

theorem WPlus_eq : WPlus = W plusCert := rfl
theorem WMinus_eq : WMinus = W minusCert := rfl

theorem finrank_WPlus : Module.finrank ℝ WPlus = 3 := finrank_W plusCert
theorem finrank_WMinus : Module.finrank ℝ WMinus = 3 := finrank_W minusCert

theorem WPlus_invariant (p : List Nat) (hp : p ∈ perms) {w : Fin 20 → ℝ} (hw : w ∈ WPlus) :
    (castZ (faceActZ p)).mulVec w ∈ WPlus := W_invariant plusCert p hp hw
theorem WMinus_invariant (p : List Nat) (hp : p ∈ perms) {w : Fin 20 → ℝ} (hw : w ∈ WMinus) :
    (castZ (faceActZ p)).mulVec w ∈ WMinus := W_invariant minusCert p hp hw

/-- The Burnside span of the `λ₊` piece has dimension nine. -/
theorem finrank_SPlus :
    Module.finrank ℝ (Submodule.span ℝ
      {X : Matrix (Fin 20) (Fin 20) ℝ | ∃ p ∈ perms, X = castZ (faceActZ p) * goldenPlusR}) = 9 :=
  finrank_S plusCert
/-- The Burnside span of the `λ₋` piece has dimension nine. -/
theorem finrank_SMinus :
    Module.finrank ℝ (Submodule.span ℝ
      {X : Matrix (Fin 20) (Fin 20) ℝ | ∃ p ∈ perms, X = castZ (faceActZ p) * goldenMinusR}) = 9 :=
  finrank_S minusCert

/-- **(4)** The `λ₊` golden piece is irreducible under the listed group. -/
theorem goldenPlus_irreducible (U : Submodule ℝ (Fin 20 → ℝ)) (hU : U ≤ WPlus)
    (hinv : ∀ p ∈ perms, ∀ u ∈ U, (castZ (faceActZ p)).mulVec u ∈ U) : U = ⊥ ∨ U = WPlus :=
  irreducible plusCert U hU hinv

/-- **(4)** The `λ₋` golden piece is irreducible under the listed group. -/
theorem goldenMinus_irreducible (U : Submodule ℝ (Fin 20 → ℝ)) (hU : U ≤ WMinus)
    (hinv : ∀ p ∈ perms, ∀ u ∈ U, (castZ (faceActZ p)).mulVec u ∈ U) : U = ⊥ ∨ U = WMinus :=
  irreducible minusCert U hU hinv

end

end OPH.GoldenSectorIrreducibility

#print axioms OPH.GoldenSectorIrreducibility.elems_mem
#print axioms OPH.GoldenSectorIrreducibility.conj_certificates
#print axioms OPH.GoldenSectorIrreducibility.BR_mul_AR
#print axioms OPH.GoldenSectorIrreducibility.AR_mul_BR
#print axioms OPH.GoldenSectorIrreducibility.gR_comm_projR
#print axioms OPH.GoldenSectorIrreducibility.W_invariant
#print axioms OPH.GoldenSectorIrreducibility.finrank_W
#print axioms OPH.GoldenSectorIrreducibility.sandwich_of_mem_S
#print axioms OPH.GoldenSectorIrreducibility.Phi_injOn_S
#print axioms OPH.GoldenSectorIrreducibility.span_gt_eq_top
#print axioms OPH.GoldenSectorIrreducibility.map_Phi_S_eq_top
#print axioms OPH.GoldenSectorIrreducibility.Phi_S_surjective
#print axioms OPH.GoldenSectorIrreducibility.finrank_S
#print axioms OPH.GoldenSectorIrreducibility.irreducible
#print axioms OPH.GoldenSectorIrreducibility.finrank_WPlus
#print axioms OPH.GoldenSectorIrreducibility.finrank_WMinus
#print axioms OPH.GoldenSectorIrreducibility.finrank_SPlus
#print axioms OPH.GoldenSectorIrreducibility.finrank_SMinus
#print axioms OPH.GoldenSectorIrreducibility.goldenPlus_irreducible
#print axioms OPH.GoldenSectorIrreducibility.goldenMinus_irreducible
