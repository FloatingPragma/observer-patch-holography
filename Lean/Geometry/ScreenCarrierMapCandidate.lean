import PortFrameGram
import SeamCurrentCarrierQuotient
import LocalFaceMaxwellAction
import A5PortAction
import DiscreteRefinement
import Geometry.CanonicalLorentzModule

set_option autoImplicit false

namespace OPH.ScreenCarrierMapCandidate

open OPH.PortFrameGram OPH.SeamCurrentCarrierQuotient OPH.LocalFaceMaxwellAction
open OPH.DiscreteRefinement (Barycentric)

/-!
# A candidate carrier map from the screen complex to spatial rays

The committed screen combinatorics consist of the twelve-port incidence
structure (`PortFrameGram.neighbors`), the thirty-seam endpoint table
(`SeamCurrentCarrierQuotient.seamLeft` and `seamRight`), the twenty
committed faces (`LocalFaceMaxwellAction.faceVertices`, pinned row for row
to the committed `CoreAxioms.orientedFaces` packet), and the sixty listed
incidence automorphisms (`A5PortAction.perms`).  The 2026-08-24 deep audit
(finding F5, completion item 2) records that the common-world
Maxwell-clock object is a formal same-index product because no map sends
ports, seams, or faces of the screen complex to geometric objects of the
canonical Lorentz module.  This file constructs the first CANDIDATE such
map, with every combinatorial fidelity clause proved against the committed
tables rather than stipulated.  It addresses the first half of the audit's
completion item 2 for issue #740.  It does not construct the second half
(a gauge-invariant interaction with source current and independent
relative normalization), and it does not discharge premise row PR-53.

WHAT IS PROVED.

* Geometric target.  Twelve exact vectors in `ℤ[φ]³`, the icosahedral
  cyclic family `(0, ±1, ±φ)`, assigned to the twelve ports
  (`candidateRayZ`) and read as rays (nonzero vectors without
  normalization).  All 144 pairwise `ℤ[φ]` inner products are computed
  exactly (`candidate_dot_table`): `2 + φ` on the diagonal, `φ` on
  committed adjacency, `-φ` at graph distance two, `-(2 + φ)` at the
  antipode.  These reproduce the committed scaled Gram matrix
  `PortFrameGram.g5` up to the single overall squared-length
  normalization (`candidate_gram_bridge`).
* Fidelity clauses.  Injectivity on ports, on seams (unordered
  endpoint-ray pairs), and on faces (unordered ray triples); endpoint and
  boundary fidelity by construction against the committed tables; and the
  two-sided adjacency correspondence: two ports share a committed seam if
  and only if their rays have exact inner product `φ`
  (`PortCarrierCandidate.sharesSeam_iff`), together with the real-number
  version (`candidateRay_adjacent_iff`).  Both directions of each
  correspondence are proved.
* Distinct rays.  Non-antipodal distinct ports receive non-parallel real
  vectors, witnessed by exact nonzero cross products
  (`candidateRay_not_parallel`); antipodal ports receive exactly opposite
  vectors (`candidate_antipode`); all twelve rays are pairwise distinct
  as rays (`candidateRay_pairwise_distinct_rays`).
* Symmetry fidelity.  The committed `A5` data is combinatorial
  (permutation lists without matrices), so this file constructs the exact
  rotation matrices itself: `rotAZ` and `rotBZ` have entries in `ℤ[φ]`,
  equal twice a proper rotation (`MᵀM = 4·1` and `det M = 8`), and
  realize the listed permutations `genA` and `genB` through the map with
  the declared factor two (`candidate_equivariant_rotA`,
  `candidate_equivariant_rotB`), hence ray-level equivariance
  (`rotA_ray_equivariant`, `rotB_ray_equivariant`).  An explicit word
  certificate (`perms_generated_by_two_rotations`) shows the two
  permutations generate the entire committed sixty-element list, so
  equivariance is proved for a generating set.
* Refinement compatibility, at the level the committed notion supports.
  The committed refinement (`DiscreteRefinement.refine`) acts on
  same-parent barycentric face coordinates by denominator multiplication.
  The carrier extension `baryCarrierZ` maps a barycentric point to the
  exact integer combination of its face's three vertex rays; `refine m`
  scales the image vector by `m` (`baryCarrier_refine`), hence fixes
  every mesh ray (`baryCarrier_refine_sameRay`); corner points map to
  the port rays (`baryCarrier_corner_first` and companions); and points
  with a nonzero coordinate sum map to nonzero vectors
  (`baryCarrier_ne_zero`).
* Non-forcing.  The negated assignment `negatedRayZ` satisfies the same
  proved clauses (`negatedCandidate`), and its rays are opposite to the
  canonical rays (`negated_opposite_ray`).  The canonical assignment is
  therefore a declared selection among valid equivariant embeddings, not
  a forced identification.

WHAT IS NOT PROVED HERE.  No metric calibration, unit, or physical scale
is attached to any ray; rays carry direction data only.  No causal or
temporal content is constructed.  No compatibility with the committed
Maxwell dynamics, the seam-current coupling, or any interaction term is
stated.  No observer readout is derived.  No theorem identifies a port
with a physical direction; premise row PR-53 is not consumed and is not
discharged.  The committed corpus carries no frequency-`n` incidence
complex (only the counts in `DiscreteRefinement.tower_counts_*` and the
same-parent barycentric coordinates), so refinement compatibility across
committed scales beyond the barycentric level is open.  Source selection
among the valid candidates is open.  The structure
`PortCarrierCandidate` is constrained by exact equations on the
committed tables; it is not a stipulable schema, and its two exhibited
inhabitants show that the proved clauses do not force a unique map.
-/

/-! ## Exact scalars: `ℤ[φ]` as integer pairs

`(a, b)` denotes `a + b·φ` with `φ² = φ + 1`.  Addition is componentwise
and is taken from the product instance; multiplication, integer scaling,
negation, and subtraction are explicit definitions so that every proof
step reduces definitionally. -/

/-- `ℤ[φ]` as pairs `(a, b)` meaning `a + b·φ`. -/
abbrev Zphi := ℤ × ℤ

/-- Multiplication in `ℤ[φ]`.  Re-derivation of the rule: with
`φ² = φ + 1`, `(a + bφ)(c + dφ) = ac + (ad + bc)φ + bd·φ²
= (ac + bd) + (ad + bc + bd)φ`. -/
def zmul (x y : Zphi) : Zphi :=
  (x.1 * y.1 + x.2 * y.2, x.1 * y.2 + x.2 * y.1 + x.2 * y.2)

/-- Scaling of a `ℤ[φ]` element by a natural number. -/
def zsc (n : ℕ) (x : Zphi) : Zphi := ((n : ℤ) * x.1, (n : ℤ) * x.2)

/-- Negation in `ℤ[φ]`. -/
def zneg (x : Zphi) : Zphi := (-x.1, -x.2)

/-- Subtraction in `ℤ[φ]`. -/
def zsub (x y : Zphi) : Zphi := (x.1 - y.1, x.2 - y.2)

theorem zmul_zero (x : Zphi) : zmul x 0 = 0 := by
  obtain ⟨a, b⟩ := x
  apply Prod.ext
  · show a * 0 + b * 0 = 0
    ring
  · show a * 0 + b * 0 + b * 0 = 0
    ring

theorem zmul_add (x y z : Zphi) : zmul x (y + z) = zmul x y + zmul x z := by
  obtain ⟨a, b⟩ := x
  obtain ⟨c, d⟩ := y
  obtain ⟨e, f⟩ := z
  apply Prod.ext
  · show a * (c + e) + b * (d + f) = a * c + b * d + (a * e + b * f)
    ring
  · show a * (d + f) + b * (c + e) + b * (d + f) =
      a * d + b * c + b * d + (a * f + b * e + b * f)
    ring

theorem zmul_zsc (n : ℕ) (x y : Zphi) : zmul x (zsc n y) = zsc n (zmul x y) := by
  obtain ⟨a, b⟩ := x
  obtain ⟨c, d⟩ := y
  apply Prod.ext
  · show a * ((n : ℤ) * c) + b * ((n : ℤ) * d) = (n : ℤ) * (a * c + b * d)
    ring
  · show a * ((n : ℤ) * d) + b * ((n : ℤ) * c) + b * ((n : ℤ) * d) =
      (n : ℤ) * (a * d + b * c + b * d)
    ring

theorem zsc_add (n : ℕ) (x y : Zphi) : zsc n (x + y) = zsc n x + zsc n y := by
  obtain ⟨a, b⟩ := x
  obtain ⟨c, d⟩ := y
  apply Prod.ext
  · show (n : ℤ) * (a + c) = (n : ℤ) * a + (n : ℤ) * c
    ring
  · show (n : ℤ) * (b + d) = (n : ℤ) * b + (n : ℤ) * d
    ring

theorem zsc_zsc (m n : ℕ) (x : Zphi) : zsc m (zsc n x) = zsc (m * n) x := by
  obtain ⟨a, b⟩ := x
  apply Prod.ext
  · show (m : ℤ) * ((n : ℤ) * a) = ((m * n : ℕ) : ℤ) * a
    push_cast
    ring
  · show (m : ℤ) * ((n : ℤ) * b) = ((m * n : ℕ) : ℤ) * b
    push_cast
    ring

theorem zsc_one (x : Zphi) : zsc 1 x = x := by
  obtain ⟨a, b⟩ := x
  apply Prod.ext
  · show ((1 : ℕ) : ℤ) * a = a
    push_cast
    ring
  · show ((1 : ℕ) : ℤ) * b = b
    push_cast
    ring

theorem zsc_zero (x : Zphi) : zsc 0 x = 0 := by
  obtain ⟨a, b⟩ := x
  apply Prod.ext
  · show ((0 : ℕ) : ℤ) * a = 0
    push_cast
    ring
  · show ((0 : ℕ) : ℤ) * b = 0
    push_cast
    ring

theorem zsc_add_left (m n : ℕ) (x : Zphi) :
    zsc (m + n) x = zsc m x + zsc n x := by
  obtain ⟨a, b⟩ := x
  apply Prod.ext
  · show ((m + n : ℕ) : ℤ) * a = (m : ℤ) * a + (n : ℤ) * a
    push_cast
    ring
  · show ((m + n : ℕ) : ℤ) * b = (m : ℤ) * b + (n : ℤ) * b
    push_cast
    ring

theorem zsub_eq_zero {x y : Zphi} (h : zsub x y = 0) : x = y := by
  obtain ⟨a, b⟩ := x
  obtain ⟨c, d⟩ := y
  have h1 : a - c = 0 := congrArg Prod.fst h
  have h2 : b - d = 0 := congrArg Prod.snd h
  have ha : a = c := by omega
  have hb : b = d := by omega
  rw [ha, hb]

/-! ## Exact vectors, inner products, cross products, and matrices -/

/-- A spatial vector with exact `ℤ[φ]` coordinates. -/
abbrev VecZ := Fin 3 → Zphi

/-- Componentwise natural scaling of an exact vector. -/
def vsc (n : ℕ) (v : VecZ) : VecZ := fun k => zsc n (v k)

/-- Componentwise negation of an exact vector. -/
def vneg (v : VecZ) : VecZ := fun k => zneg (v k)

/-- Exact `ℤ[φ]` inner product of two vectors. -/
def dotZ (u v : VecZ) : Zphi :=
  zmul (u 0) (v 0) + zmul (u 1) (v 1) + zmul (u 2) (v 2)

/-- Exact cross product; a nonzero value witnesses non-parallelism of the
real evaluations. -/
def crossZ (u v : VecZ) : VecZ
  | 0 => zsub (zmul (u 1) (v 2)) (zmul (u 2) (v 1))
  | 1 => zsub (zmul (u 2) (v 0)) (zmul (u 0) (v 2))
  | 2 => zsub (zmul (u 0) (v 1)) (zmul (u 1) (v 0))

/-- Exact matrix action on a vector. -/
def matVecZ (M : Fin 3 → Fin 3 → Zphi) (v : VecZ) : VecZ :=
  fun r => zmul (M r 0) (v 0) + zmul (M r 1) (v 1) + zmul (M r 2) (v 2)

/-- Exact determinant of a three-by-three `ℤ[φ]` matrix. -/
def detZ (M : Fin 3 → Fin 3 → Zphi) : Zphi :=
  zsub (zmul (M 0 0) (zsub (zmul (M 1 1) (M 2 2)) (zmul (M 1 2) (M 2 1))))
      (zmul (M 0 1) (zsub (zmul (M 1 0) (M 2 2)) (zmul (M 1 2) (M 2 0)))) +
    zmul (M 0 2) (zsub (zmul (M 1 0) (M 2 1)) (zmul (M 1 1) (M 2 0)))

theorem dotZ_zero_right (u : VecZ) : dotZ u 0 = 0 := by
  show zmul (u 0) 0 + zmul (u 1) 0 + zmul (u 2) 0 = 0
  rw [zmul_zero, zmul_zero, zmul_zero, add_zero, add_zero]

theorem dotZ_add_right (u v w : VecZ) :
    dotZ u (v + w) = dotZ u v + dotZ u w := by
  show zmul (u 0) (v 0 + w 0) + zmul (u 1) (v 1 + w 1) + zmul (u 2) (v 2 + w 2) =
    (zmul (u 0) (v 0) + zmul (u 1) (v 1) + zmul (u 2) (v 2)) +
      (zmul (u 0) (w 0) + zmul (u 1) (w 1) + zmul (u 2) (w 2))
  rw [zmul_add, zmul_add, zmul_add]
  abel

theorem dotZ_vsc_right (n : ℕ) (u v : VecZ) :
    dotZ u (vsc n v) = zsc n (dotZ u v) := by
  show zmul (u 0) (zsc n (v 0)) + zmul (u 1) (zsc n (v 1)) +
      zmul (u 2) (zsc n (v 2)) =
    zsc n (zmul (u 0) (v 0) + zmul (u 1) (v 1) + zmul (u 2) (v 2))
  rw [zmul_zsc, zmul_zsc, zmul_zsc, zsc_add, zsc_add]

/-! ## The candidate port assignment

Re-derivation of the target values, recorded against transcription error:
in the cyclic icosahedral family `(0, ±1, ±φ)` every vector has squared
length `1 + φ² = 2 + φ`; two adjacent vertices, for example `(0, 1, φ)`
and `(1, φ, 0)`, have inner product `φ`; a vertex and its antipode have
inner product `-(2 + φ)`; the remaining pairs, for example `(0, 1, φ)`
and `(0, 1, -φ)`, have inner product `1 - φ² = -φ`.  The assignment below
was searched against the committed adjacency `PortFrameGram.neighbors`
under the committed antipode `i ↦ 11 - i`, and the theorems that follow
re-verify every claimed value by kernel computation. -/

/-- The candidate ray of each port: an exact icosahedral direction vector.
Port `11 - i` receives the negated vector of port `i`. -/
def candidateRayZ : Fin 12 → VecZ
  | 0 => ![(0, 0), (1, 0), (0, 1)]
  | 1 => ![(1, 0), (0, 1), (0, 0)]
  | 2 => ![(0, 1), (0, 0), (1, 0)]
  | 3 => ![(-1, 0), (0, 1), (0, 0)]
  | 4 => ![(0, 0), (-1, 0), (0, 1)]
  | 5 => ![(0, 1), (0, 0), (-1, 0)]
  | 6 => ![(0, -1), (0, 0), (1, 0)]
  | 7 => ![(0, 0), (1, 0), (0, -1)]
  | 8 => ![(1, 0), (0, -1), (0, 0)]
  | 9 => ![(0, -1), (0, 0), (-1, 0)]
  | 10 => ![(-1, 0), (0, -1), (0, 0)]
  | _ => ![(0, 0), (-1, 0), (0, -1)]

/-- The exact Gram value the candidate must reproduce at each pair:
`2 + φ` on the diagonal, `φ` on committed adjacency, `-(2 + φ)` at the
committed antipode, `-φ` otherwise. -/
def gramTargetZ (i j : Fin 12) : Zphi :=
  if i = j then (2, 1)
  else if adj i j then (0, 1)
  else if j = antipode i then (-2, -1)
  else (0, -1)

theorem gramTarget_diag (i : Fin 12) : gramTargetZ i i = (2, 1) := by
  simp [gramTargetZ]

theorem gramTarget_offdiag_ne :
    ∀ i j : Fin 12, i ≠ j → gramTargetZ i j ≠ (2, 1) := by decide

theorem gramTarget_phi_iff_adj :
    ∀ i j : Fin 12, gramTargetZ i j = (0, 1) ↔ adj i j = true := by decide

set_option maxHeartbeats 1000000 in
set_option maxRecDepth 8192 in
/-- All 144 exact inner products of the candidate rays match the target
table.  This is the load-bearing computation of the file. -/
theorem candidate_dot_table :
    ∀ i j : Fin 12, dotZ (candidateRayZ i) (candidateRayZ j) = gramTargetZ i j := by
  decide

/-- Conversion of `2·(a + bφ)` to the `ℤ(√5)` pair basis of
`PortFrameGram`: `2(a + bφ) = (2a + b) + b·√5`. -/
def twoInSqrt5 (x : Zphi) : ℤ × ℤ := (2 * x.1 + x.2, x.2)

set_option maxHeartbeats 1000000 in
set_option maxRecDepth 8192 in
/-- The candidate reproduces the committed scaled Gram matrix `g5` exactly,
up to the single overall squared-length normalization.  Re-derivation of
the constant, recorded against transcription error: `g5 = 5·G` with
`G i j = dot i j / (2 + φ)`, so the claim `5·dot = (2 + φ)·g5·(1/5)·5`,
cleared of halves by doubling, reads `5·(2·dot) = (5 + √5)·g5`; equality
of the constants reduces at the diagonal pair to `10·(2 + φ) = 5·(5 + √5)`,
equivalently `2φ - 1 = √5`, which holds for `φ = (1 + √5)/2`.
Both sides below are exact `ℤ(√5)` products in the committed `mulZ5`. -/
theorem candidate_gram_bridge :
    ∀ i j : Fin 12,
      mulZ5 (5, 0) (twoInSqrt5 (dotZ (candidateRayZ i) (candidateRayZ j))) =
        mulZ5 (5, 1) (g5 i j) := by
  decide

/-- The candidate intertwines the committed antipode with exact vector
negation. -/
theorem candidate_antipode :
    ∀ i : Fin 12, candidateRayZ (antipode i) = vneg (candidateRayZ i) := by
  decide

/-- Every candidate ray has a nonzero exact coordinate. -/
theorem candidate_coord_ne_zero :
    ∀ i : Fin 12, ∃ k : Fin 3, candidateRayZ i k ≠ 0 := by decide

set_option maxHeartbeats 1000000 in
set_option maxRecDepth 8192 in
/-- Distinct non-antipodal ports have an exact nonzero cross product.
Antipodal ports are excluded because their vectors are negatives of each
other, hence parallel: the exhibited assignment carries the committed Gram
value `-1` at the antipode as exact vector negation
(`candidate_antipode`), so its antipodal pairs are exactly opposite.  That
every table-matching assignment shares this necessity is a Cauchy-Schwarz
equality argument recorded here as prose, without a named theorem. -/
theorem candidate_cross_ne_zero :
    ∀ i j : Fin 12, i ≠ j → j ≠ antipode i →
      ∃ k : Fin 3, crossZ (candidateRayZ i) (candidateRayZ j) k ≠ 0 := by
  decide

/-! ## The two exact generator rotations

The committed `A5` data (`A5PortAction.perms`) is combinatorial: sixty
permutation rows without matrices.  The route taken here is therefore the
constructive one: two rotation matrices are built from scratch with
`ℤ[φ]` entries and proved to realize two listed permutations through the
candidate map; a word certificate then shows those two permutations
generate the whole listed group. -/

/-- Twice the rotation by `2π/3` about the axis `(1, 1, 1)`: the doubled
coordinate cycle `(x, y, z) ↦ (z, x, y)`. -/
def rotAZ : Fin 3 → Fin 3 → Zphi
  | 0, 2 => (2, 0)
  | 1, 0 => (2, 0)
  | 2, 1 => (2, 0)
  | _, _ => (0, 0)

/-- Twice a rotation by `2π/5` about a five-fold icosahedral axis, with
entries in `{±1, ±φ, ±(φ - 1)}`.  The true rotation is `rotBZ / 2`. -/
def rotBZ : Fin 3 → Fin 3 → Zphi
  | 0, 0 => (-1, 1)
  | 0, 1 => (0, 1)
  | 0, _ => (-1, 0)
  | 1, 0 => (0, -1)
  | 1, 1 => (1, 0)
  | 1, _ => (-1, 1)
  | _, 0 => (1, 0)
  | _, 1 => (-1, 1)
  | _, _ => (0, 1)

/-- The port permutation realized by `rotAZ` (row 10 of the committed
list). -/
def genA : Fin 12 → Fin 12 := ![2, 0, 1, 4, 5, 3, 8, 6, 7, 10, 11, 9]

/-- The port permutation realized by `rotBZ` (row 1 of the committed
list), of order five. -/
def genB : Fin 12 → Fin 12 := ![0, 2, 4, 1, 6, 8, 3, 5, 10, 7, 9, 11]

theorem genA_listed :
    (List.finRange 12).map (fun k => (genA k : ℕ)) ∈ OPH.A5PortAction.perms := by
  decide

theorem genB_listed :
    (List.finRange 12).map (fun k => (genB k : ℕ)) ∈ OPH.A5PortAction.perms := by
  decide

/-- `rotAZᵀ · rotAZ = 4·1`: the halved matrix is orthogonal. -/
theorem rotA_orthogonal :
    ∀ r c : Fin 3,
      zmul (rotAZ 0 r) (rotAZ 0 c) + zmul (rotAZ 1 r) (rotAZ 1 c) +
          zmul (rotAZ 2 r) (rotAZ 2 c) =
        if r = c then (4, 0) else 0 := by
  decide

/-- `rotBZᵀ · rotBZ = 4·1`: the halved matrix is orthogonal. -/
theorem rotB_orthogonal :
    ∀ r c : Fin 3,
      zmul (rotBZ 0 r) (rotBZ 0 c) + zmul (rotBZ 1 r) (rotBZ 1 c) +
          zmul (rotBZ 2 r) (rotBZ 2 c) =
        if r = c then (4, 0) else 0 := by
  decide

/-- `det rotAZ = 8`.  Re-derivation of the sign convention: for a doubled
matrix `M = 2R` in three dimensions, `det M = 2³ · det R`, so `det M = 8`
certifies `det R = +1`, a proper rotation and not a reflection. -/
theorem rotA_det : detZ rotAZ = (8, 0) := by decide

/-- `det rotBZ = 8`: the halved matrix is a proper rotation
(`det R = +1`), by the same doubling argument as `rotA_det`. -/
theorem rotB_det : detZ rotBZ = (8, 0) := by decide

/-- `rotAZ` realizes the listed permutation `genA` through the candidate
map, with the declared doubling factor. -/
theorem candidate_equivariant_rotA :
    ∀ (i : Fin 12) (k : Fin 3),
      matVecZ rotAZ (candidateRayZ i) k =
        zmul (2, 0) (candidateRayZ (genA i) k) := by
  decide

set_option maxHeartbeats 1000000 in
/-- `rotBZ` realizes the listed permutation `genB` through the candidate
map, with the declared doubling factor. -/
theorem candidate_equivariant_rotB :
    ∀ (i : Fin 12) (k : Fin 3),
      matVecZ rotBZ (candidateRayZ i) k =
        zmul (2, 0) (candidateRayZ (genB i) k) := by
  decide

/-- Apply a word in the two generators to a port, letters read left to
right (`true` applies `genA`, `false` applies `genB`). -/
def wordApply (w : List Bool) (k : Fin 12) : Fin 12 :=
  w.foldl (fun x b => if b then genA x else genB x) k

/-- One generator word per committed permutation row, in row order. -/
def genWords : List (List Bool) := [
  [],
  [false],
  [false, false, false, false],
  [false, false],
  [false, false, false],
  [false, true, true],
  [true, true],
  [false, false, true, true],
  [true, true, false, true],
  [false, false, false, true, true],
  [true],
  [false, true],
  [true, true, false],
  [false, false, true],
  [false, false, false, true],
  [true, false, false, false],
  [true, false, true, true],
  [false, true, false, false, false],
  [true, true, false, false, true, true],
  [false, false, true, false, false, false],
  [true, false],
  [false, true, false],
  [true, true, false, false],
  [false, false, true, false],
  [false, false, false, true, false],
  [false, true, false, true],
  [true, false, true],
  [false, false, true, false, true],
  [true, true, false, false, true],
  [false, false, false, true, false, true],
  [true, false, false],
  [true, true, false, false, false],
  [false, true, false, false],
  [false, false, false, true, false, false],
  [false, false, true, false, false],
  [true, false, false, true, true],
  [false, true, false, false, true, true],
  [true, false, true, false, true],
  [false, false, true, false, false, true, true],
  [true, true, false, false, true, false, true],
  [true, false, false, true],
  [true, false, true, false],
  [false, true, false, false, true],
  [true, true, false, false, true, false],
  [false, false, true, false, false, true],
  [true, false, true, false, false, false],
  [true, false, false, true, false, false],
  [true, false, true, false, false, true, true],
  [false, true, false, false, true, false, false],
  [true, false, true, false, false, true, false, true],
  [true, false, false, true, false],
  [true, false, true, false, false],
  [false, true, false, false, true, false],
  [true, true, false, false, true, false, false],
  [false, false, true, false, false, true, false],
  [true, false, false, true, false, true],
  [false, true, false, false, true, false, true],
  [true, false, true, false, false, true],
  [true, false, true, false, false, true, false, false],
  [true, false, true, false, false, true, false]]

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Word certificate: every committed permutation row is the fold of an
explicit word in `genA` and `genB`, so the two realized permutations
generate the whole committed sixty-element list. -/
theorem perms_generated_by_two_rotations :
    OPH.A5PortAction.perms =
      genWords.map (fun w => (List.finRange 12).map fun k => (wordApply w k : ℕ)) := by
  decide

/-! ## The candidate structure

Fields are exact equations on the committed tables.  The structure is not
a stipulable schema: an inhabitant must reproduce all 144 committed Gram
values, up to the single squared-length normalization of the bridge
theorem, and both generator actions.  Two distinct inhabitants are
exhibited below (`canonicalCandidate` and `negatedCandidate`), so the
clauses select a class of maps and do not force a unique one. -/

/-- A candidate carrier assignment of exact rays to the twelve ports,
constrained by the committed Gram table and by equivariance under the two
exact generator rotations.  Candidate status: inhabitants are declared
selections among valid equivariant embeddings; no field attaches units,
scale, causal content, or dynamics. -/
structure PortCarrierCandidate where
  /-- The exact ray assigned to each port. -/
  portMap : Fin 12 → VecZ
  /-- All 144 exact inner products match the committed Gram target. -/
  dot_table : ∀ i j : Fin 12, dotZ (portMap i) (portMap j) = gramTargetZ i j
  /-- `rotAZ` realizes `genA` through the map, with the doubling factor. -/
  equivariant_genA : ∀ (i : Fin 12) (k : Fin 3),
    matVecZ rotAZ (portMap i) k = zmul (2, 0) (portMap (genA i) k)
  /-- `rotBZ` realizes `genB` through the map, with the doubling factor. -/
  equivariant_genB : ∀ (i : Fin 12) (k : Fin 3),
    matVecZ rotBZ (portMap i) k = zmul (2, 0) (portMap (genB i) k)

/-- Port injectivity holds for every inhabitant: the diagonal Gram value
`2 + φ` is attained at no off-diagonal pair. -/
theorem PortCarrierCandidate.portMap_injective (c : PortCarrierCandidate) :
    Function.Injective c.portMap := by
  intro i j h
  by_contra hij
  have hdot := c.dot_table i j
  rw [← h] at hdot
  rw [c.dot_table i i, gramTarget_diag] at hdot
  exact gramTarget_offdiag_ne i j hij hdot.symm

/-- Adjacency correspondence for every inhabitant, both directions: two
ports are committed neighbors exactly when their rays have exact inner
product `φ`. -/
theorem PortCarrierCandidate.adjacent_iff (c : PortCarrierCandidate)
    (i j : Fin 12) :
    adj i j = true ↔ dotZ (c.portMap i) (c.portMap j) = (0, 1) := by
  rw [c.dot_table i j]
  exact (gramTarget_phi_iff_adj i j).symm

/-- Candidate seam image: the unordered pair of the two endpoint rays of
the committed seam. -/
def PortCarrierCandidate.seamImage (c : PortCarrierCandidate) (e : Fin 30) :
    Sym2 VecZ :=
  s(c.portMap (seamLeft e), c.portMap (seamRight e))

/-- Endpoint fidelity, recorded for inspection: the seam image is exactly
the image of the committed endpoint pair.  The equation holds by
construction; the committed content behind it is `seam_table_sound` and
`seam_table_complete`. -/
theorem PortCarrierCandidate.seamImage_endpoints (c : PortCarrierCandidate)
    (e : Fin 30) :
    c.seamImage e = s(c.portMap (seamLeft e), c.portMap (seamRight e)) := rfl

/-- Seam injectivity for every inhabitant: distinct committed seams have
distinct unordered endpoint-ray pairs. -/
theorem PortCarrierCandidate.seamImage_injective (c : PortCarrierCandidate) :
    Function.Injective c.seamImage := by
  intro e e' h
  rw [PortCarrierCandidate.seamImage, PortCarrierCandidate.seamImage,
    Sym2.eq_iff] at h
  rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩
  · have hl := c.portMap_injective h1
    have hr := c.portMap_injective h2
    apply seam_table_injective
    show (seamLeft e, seamRight e) = (seamLeft e', seamRight e')
    rw [hl, hr]
  · have hl := c.portMap_injective h1
    have hr := c.portMap_injective h2
    have s1 := (seam_table_sound e).1
    have s2 := (seam_table_sound e').1
    rw [hl, hr] at s1
    exact absurd s2 (lt_asymm s1)

/-- The unordered port triple of a committed face. -/
def facePortTriple (f : Fin 20) : Multiset (Fin 12) :=
  {(faceVertices f).1, (faceVertices f).2.1, (faceVertices f).2.2}

/-- Distinct committed faces have distinct unordered port triples. -/
theorem facePortTriple_injective :
    ∀ f g : Fin 20, facePortTriple f = facePortTriple g → f = g := by decide

/-- Candidate face image: the unordered triple of the three port rays of
the committed face. -/
def PortCarrierCandidate.faceImage (c : PortCarrierCandidate) (f : Fin 20) :
    Multiset VecZ :=
  (facePortTriple f).map c.portMap

/-- Boundary fidelity, recorded for inspection: the face image is exactly
the image of the committed face-port triple.  The committed content
behind the triple is `faceVertices_matches_committed`, which pins
`faceVertices` to the `CoreAxioms.orientedFaces` packet. -/
theorem PortCarrierCandidate.faceImage_eq_triple (c : PortCarrierCandidate)
    (f : Fin 20) :
    c.faceImage f =
      {c.portMap (faceVertices f).1, c.portMap (faceVertices f).2.1,
        c.portMap (faceVertices f).2.2} := by
  simp [PortCarrierCandidate.faceImage, facePortTriple]

/-- Face injectivity for every inhabitant. -/
theorem PortCarrierCandidate.faceImage_injective (c : PortCarrierCandidate) :
    Function.Injective c.faceImage := fun f g h =>
  facePortTriple_injective f g (Multiset.map_injective c.portMap_injective h)

/-- Seam correspondence for every inhabitant, both directions: two ports
share a committed seam exactly when their rays have exact inner product
`φ`.  The forward direction uses `seam_table_sound`; the reverse uses
`seam_table_complete`. -/
theorem PortCarrierCandidate.sharesSeam_iff (c : PortCarrierCandidate)
    (i j : Fin 12) :
    (∃ e : Fin 30,
        (seamLeft e = i ∧ seamRight e = j) ∨
          (seamLeft e = j ∧ seamRight e = i)) ↔
      dotZ (c.portMap i) (c.portMap j) = (0, 1) := by
  rw [← c.adjacent_iff i j]
  constructor
  · rintro ⟨e, ⟨hl, hr⟩ | ⟨hl, hr⟩⟩
    · rw [← hl, ← hr]
      exact (seam_table_sound e).2
    · rw [← hl, ← hr, adj_symm]
      exact (seam_table_sound e).2
  · intro hadj
    rcases lt_trichotomy i j with hlt | heq | hgt
    · obtain ⟨e, he1, he2⟩ := (seam_table_complete i j).mp ⟨hlt, hadj⟩
      exact ⟨e, Or.inl ⟨he1, he2⟩⟩
    · subst heq
      rw [adj_irrefl] at hadj
      exact absurd hadj Bool.false_ne_true
    · obtain ⟨e, he1, he2⟩ :=
        (seam_table_complete j i).mp ⟨hgt, by rw [adj_symm]; exact hadj⟩
      exact ⟨e, Or.inr ⟨he1, he2⟩⟩

/-- The canonical candidate: the explicit assignment `candidateRayZ`. -/
def canonicalCandidate : PortCarrierCandidate where
  portMap := candidateRayZ
  dot_table := candidate_dot_table
  equivariant_genA := candidate_equivariant_rotA
  equivariant_genB := candidate_equivariant_rotB

/-- Port injectivity of the canonical candidate, as a named corollary. -/
theorem candidateRayZ_injective : Function.Injective candidateRayZ :=
  canonicalCandidate.portMap_injective

/-! ## The negated candidate: the fidelity clauses do not force the map -/

/-- The negated assignment: every port ray replaced by its negative. -/
def negatedRayZ (i : Fin 12) : VecZ := vneg (candidateRayZ i)

set_option maxHeartbeats 1000000 in
set_option maxRecDepth 8192 in
/-- The negated assignment reproduces the same committed Gram table. -/
theorem negated_dot_table :
    ∀ i j : Fin 12, dotZ (negatedRayZ i) (negatedRayZ j) = gramTargetZ i j := by
  decide

theorem negated_equivariant_rotA :
    ∀ (i : Fin 12) (k : Fin 3),
      matVecZ rotAZ (negatedRayZ i) k = zmul (2, 0) (negatedRayZ (genA i) k) := by
  decide

set_option maxHeartbeats 1000000 in
theorem negated_equivariant_rotB :
    ∀ (i : Fin 12) (k : Fin 3),
      matVecZ rotBZ (negatedRayZ i) k = zmul (2, 0) (negatedRayZ (genB i) k) := by
  decide

/-- The negated candidate: a second inhabitant satisfying every clause,
so the clauses declare a selection and do not force the canonical map. -/
def negatedCandidate : PortCarrierCandidate where
  portMap := negatedRayZ
  dot_table := negated_dot_table
  equivariant_genA := negated_equivariant_rotA
  equivariant_genB := negated_equivariant_rotB

/-- The two exhibited inhabitants differ at every port. -/
theorem negated_ne_canonical : ∀ i : Fin 12, negatedRayZ i ≠ candidateRayZ i := by
  decide

/-! ## Real evaluation: rays in the spatial slice of the Lorentz module -/

noncomputable section

/-- Evaluation of `ℤ[φ]` in `ℝ`. -/
def evalPhi (x : Zphi) : ℝ := (x.1 : ℝ) + (x.2 : ℝ) * Real.goldenRatio

theorem evalPhi_zmul (x y : Zphi) :
    evalPhi (zmul x y) = evalPhi x * evalPhi y := by
  obtain ⟨a, b⟩ := x
  obtain ⟨c, d⟩ := y
  show ((a * c + b * d : ℤ) : ℝ) +
      ((a * d + b * c + b * d : ℤ) : ℝ) * Real.goldenRatio =
    ((a : ℝ) + (b : ℝ) * Real.goldenRatio) *
      ((c : ℝ) + (d : ℝ) * Real.goldenRatio)
  have h : Real.goldenRatio ^ 2 = Real.goldenRatio + 1 := Real.goldenRatio_sq
  push_cast
  linear_combination (-(b : ℝ) * (d : ℝ)) * h

theorem evalPhi_add (x y : Zphi) : evalPhi (x + y) = evalPhi x + evalPhi y := by
  obtain ⟨a, b⟩ := x
  obtain ⟨c, d⟩ := y
  show ((a + c : ℤ) : ℝ) + ((b + d : ℤ) : ℝ) * Real.goldenRatio =
    ((a : ℝ) + (b : ℝ) * Real.goldenRatio) +
      ((c : ℝ) + (d : ℝ) * Real.goldenRatio)
  push_cast
  ring

theorem evalPhi_zsub (x y : Zphi) :
    evalPhi (zsub x y) = evalPhi x - evalPhi y := by
  obtain ⟨a, b⟩ := x
  obtain ⟨c, d⟩ := y
  show ((a - c : ℤ) : ℝ) + ((b - d : ℤ) : ℝ) * Real.goldenRatio =
    ((a : ℝ) + (b : ℝ) * Real.goldenRatio) -
      ((c : ℝ) + (d : ℝ) * Real.goldenRatio)
  push_cast
  ring

theorem evalPhi_zneg (x : Zphi) : evalPhi (zneg x) = -evalPhi x := by
  obtain ⟨a, b⟩ := x
  show ((-a : ℤ) : ℝ) + ((-b : ℤ) : ℝ) * Real.goldenRatio =
    -((a : ℝ) + (b : ℝ) * Real.goldenRatio)
  push_cast
  ring

theorem evalPhi_zsc (n : ℕ) (x : Zphi) :
    evalPhi (zsc n x) = (n : ℝ) * evalPhi x := by
  obtain ⟨a, b⟩ := x
  show (((n : ℤ) * a : ℤ) : ℝ) + (((n : ℤ) * b : ℤ) : ℝ) * Real.goldenRatio =
    (n : ℝ) * ((a : ℝ) + (b : ℝ) * Real.goldenRatio)
  push_cast
  ring

theorem evalPhi_phi : evalPhi ((0 : ℤ), (1 : ℤ)) = Real.goldenRatio := by
  show ((0 : ℤ) : ℝ) + ((1 : ℤ) : ℝ) * Real.goldenRatio = Real.goldenRatio
  push_cast
  ring

theorem evalPhi_two : evalPhi ((2 : ℤ), (0 : ℤ)) = 2 := by
  show ((2 : ℤ) : ℝ) + ((0 : ℤ) : ℝ) * Real.goldenRatio = 2
  push_cast
  ring

/-- A nonzero `ℤ[φ]` element evaluates to a nonzero real, by the
irrationality of `φ`. -/
theorem evalPhi_ne_zero {x : Zphi} (hx : x ≠ 0) : evalPhi x ≠ 0 := by
  obtain ⟨a, b⟩ := x
  intro h
  have hphi : (a : ℝ) + (b : ℝ) * Real.goldenRatio = 0 := h
  by_cases hb : b = 0
  · subst hb
    push_cast at hphi
    have ha : a = 0 := by exact_mod_cast (by linarith : (a : ℝ) = 0)
    exact hx (by rw [ha]; rfl)
  · apply Real.goldenRatio_irrational
    refine ⟨(-a : ℚ) / (b : ℚ), ?_⟩
    have hbR : (b : ℝ) ≠ 0 := Int.cast_ne_zero.mpr hb
    have hg : Real.goldenRatio = (1 + Real.sqrt 5) / 2 := rfl
    rw [hg] at hphi
    push_cast
    rw [hg]
    field_simp
    linarith

theorem evalPhi_injective : Function.Injective evalPhi := by
  intro x y h
  by_contra hne
  have hz : zsub x y ≠ 0 := fun hz0 => hne (zsub_eq_zero hz0)
  exact evalPhi_ne_zero hz (by rw [evalPhi_zsub, h, sub_self])

/-- Componentwise evaluation into the spatial slice of the canonical
Lorentz module. -/
def evalVec (v : VecZ) : OPH.C1Lorentz.Spatial := fun k => evalPhi (v k)

theorem evalVec_injective : Function.Injective evalVec := by
  intro u v h
  funext k
  exact evalPhi_injective (congrFun h k)

theorem evalVec_vneg (v : VecZ) : evalVec (vneg v) = -evalVec v := by
  funext k
  show evalPhi (zneg (v k)) = -evalPhi (v k)
  exact evalPhi_zneg (v k)

theorem evalVec_vsc (n : ℕ) (v : VecZ) :
    evalVec (vsc n v) = (n : ℝ) • evalVec v := by
  funext k
  show evalPhi (zsc n (v k)) = (n : ℝ) * evalPhi (v k)
  exact evalPhi_zsc n (v k)

/-- The exact inner product evaluates to the committed spatial pairing of
the Lorentz module. -/
theorem evalPhi_dotZ (u v : VecZ) :
    evalPhi (dotZ u v) = OPH.C1Lorentz.spatialDot (evalVec u) (evalVec v) := by
  unfold dotZ OPH.C1Lorentz.spatialDot
  rw [Fin.sum_univ_three, evalPhi_add, evalPhi_add, evalPhi_zmul, evalPhi_zmul,
    evalPhi_zmul]
  rfl

/-- Real cross product on the spatial slice. -/
def realCross (x y : OPH.C1Lorentz.Spatial) : OPH.C1Lorentz.Spatial
  | 0 => x 1 * y 2 - x 2 * y 1
  | 1 => x 2 * y 0 - x 0 * y 2
  | 2 => x 0 * y 1 - x 1 * y 0

theorem realCross_smul_self (x : OPH.C1Lorentz.Spatial) (c : ℝ) (k : Fin 3) :
    realCross x (c • x) k = 0 := by
  fin_cases k
  · show x 1 * (c * x 2) - x 2 * (c * x 1) = 0
    ring
  · show x 2 * (c * x 0) - x 0 * (c * x 2) = 0
    ring
  · show x 0 * (c * x 1) - x 1 * (c * x 0) = 0
    ring

theorem evalPhi_crossZ (u v : VecZ) (k : Fin 3) :
    evalPhi (crossZ u v k) = realCross (evalVec u) (evalVec v) k := by
  fin_cases k
  · show evalPhi (zsub (zmul (u 1) (v 2)) (zmul (u 2) (v 1))) =
      evalVec u 1 * evalVec v 2 - evalVec u 2 * evalVec v 1
    rw [evalPhi_zsub, evalPhi_zmul, evalPhi_zmul]
    rfl
  · show evalPhi (zsub (zmul (u 2) (v 0)) (zmul (u 0) (v 2))) =
      evalVec u 2 * evalVec v 0 - evalVec u 0 * evalVec v 2
    rw [evalPhi_zsub, evalPhi_zmul, evalPhi_zmul]
    rfl
  · show evalPhi (zsub (zmul (u 0) (v 1)) (zmul (u 1) (v 0))) =
      evalVec u 0 * evalVec v 1 - evalVec u 1 * evalVec v 0
    rw [evalPhi_zsub, evalPhi_zmul, evalPhi_zmul]
    rfl

/-- The candidate ray of a port as a real spatial vector.  A ray carries
direction data only; no unit or scale is attached. -/
def candidateRay (i : Fin 12) : OPH.C1Lorentz.Spatial :=
  evalVec (candidateRayZ i)

theorem candidateRay_ne_zero (i : Fin 12) : candidateRay i ≠ 0 := by
  obtain ⟨k, hk⟩ := candidate_coord_ne_zero i
  intro h
  apply evalPhi_ne_zero hk
  have hcomp := congrFun h k
  simpa [candidateRay, evalVec] using hcomp

/-- Non-parallelism for distinct non-antipodal ports: no real scalar
carries one candidate ray to the other. -/
theorem candidateRay_not_parallel (i j : Fin 12) (hij : i ≠ j)
    (hanti : j ≠ antipode i) (c : ℝ) :
    c • candidateRay i ≠ candidateRay j := by
  intro hpar
  obtain ⟨k, hk⟩ := candidate_cross_ne_zero i j hij hanti
  apply evalPhi_ne_zero hk
  have hbridge : evalPhi (crossZ (candidateRayZ i) (candidateRayZ j) k) =
      realCross (candidateRay i) (candidateRay j) k :=
    evalPhi_crossZ _ _ _
  rw [hbridge, ← hpar]
  exact realCross_smul_self (candidateRay i) c k

/-- A nonzero vector and its negative lie on distinct rays. -/
theorem not_sameRay_neg_self {x : OPH.C1Lorentz.Spatial} (hx : x ≠ 0) :
    ¬SameRay ℝ x (-x) := by
  intro h
  obtain ⟨r, s, hr, hs, hrs⟩ := h.exists_pos hx (neg_ne_zero.mpr hx)
  rw [smul_neg] at hrs
  have hsum : (r + s) • x = 0 := by
    rw [add_smul, hrs]
    abel
  rcases smul_eq_zero.mp hsum with h0 | h0
  · exact absurd h0 (ne_of_gt (add_pos hr hs))
  · exact hx h0

/-- The candidate antipode relation at the real level. -/
theorem candidateRay_antipode (i : Fin 12) :
    candidateRay (antipode i) = -candidateRay i := by
  show evalVec (candidateRayZ (antipode i)) = -candidateRay i
  rw [candidate_antipode i]
  exact evalVec_vneg _

/-- All twelve candidate rays are pairwise distinct as rays: for distinct
ports the `SameRay` relation fails.  Non-antipodal pairs are non-parallel
outright; antipodal pairs are parallel with opposite orientation, exactly
as forced by the committed Gram value `-1` at the antipode. -/
theorem candidateRay_pairwise_distinct_rays (i j : Fin 12) (hij : i ≠ j) :
    ¬SameRay ℝ (candidateRay i) (candidateRay j) := by
  by_cases hanti : j = antipode i
  · subst hanti
    rw [candidateRay_antipode]
    exact not_sameRay_neg_self (candidateRay_ne_zero i)
  · intro h
    obtain ⟨r, s, hr, hs, hrs⟩ :=
      h.exists_pos (candidateRay_ne_zero i) (candidateRay_ne_zero j)
    apply candidateRay_not_parallel i j hij hanti (r / s)
    have hs' : s ≠ 0 := ne_of_gt hs
    calc (r / s) • candidateRay i
        = s⁻¹ • (r • candidateRay i) := by rw [div_eq_inv_mul, mul_smul]
      _ = s⁻¹ • (s • candidateRay j) := by rw [hrs]
      _ = candidateRay j := by rw [smul_smul, inv_mul_cancel₀ hs', one_smul]

/-- The candidate ray embedded in the canonical Hermitian Lorentz module,
with zero scalar component. -/
def candidateRayHerm (i : Fin 12) : OPH.C1Lorentz.Herm2 :=
  ((0 : ℝ), candidateRay i)

/-- Every candidate ray is a spacelike direction of the canonical Lorentz
module: the determinant quadratic form is strictly negative on it.  This
consumes the committed inertia certificate `spatial_axis_negative` and
adds no metric calibration. -/
theorem candidateRayHerm_spacelike (i : Fin 12) :
    OPH.C1Lorentz.lorentzQ (candidateRayHerm i) < 0 :=
  OPH.C1Lorentz.spatial_axis_negative (candidateRay_ne_zero i)

/-- Real adjacency correspondence, both directions: two ports are
committed neighbors exactly when their real rays have spatial pairing
`φ`.  The reverse direction uses the irrationality of `φ` through
`evalPhi_injective`. -/
theorem candidateRay_adjacent_iff (i j : Fin 12) :
    adj i j = true ↔
      OPH.C1Lorentz.spatialDot (candidateRay i) (candidateRay j) =
        Real.goldenRatio := by
  have hd := candidate_dot_table i j
  have hbridge : OPH.C1Lorentz.spatialDot (candidateRay i) (candidateRay j) =
      evalPhi (dotZ (candidateRayZ i) (candidateRayZ j)) :=
    (evalPhi_dotZ _ _).symm
  constructor
  · intro h
    rw [hbridge, hd, (gramTarget_phi_iff_adj i j).mpr h, evalPhi_phi]
  · intro h
    rw [hbridge, hd] at h
    apply (gramTarget_phi_iff_adj i j).mp
    apply evalPhi_injective
    rw [h, evalPhi_phi]

/-- Ray-level equivariance for the first generator: the rotation
`rotAZ / 2` carries the ray of port `i` to the ray of port `genA i`. -/
theorem rotA_ray_equivariant (i : Fin 12) :
    SameRay ℝ (evalVec (matVecZ rotAZ (candidateRayZ i)))
      (candidateRay (genA i)) := by
  have hfun : evalVec (matVecZ rotAZ (candidateRayZ i)) =
      (2 : ℝ) • candidateRay (genA i) := by
    funext k
    show evalPhi (matVecZ rotAZ (candidateRayZ i) k) =
      (2 : ℝ) * evalPhi (candidateRayZ (genA i) k)
    rw [candidate_equivariant_rotA i k, evalPhi_zmul, evalPhi_two]
  rw [hfun]
  exact (SameRay.sameRay_nonneg_smul_right _ (by norm_num)).symm

/-- Ray-level equivariance for the second generator. -/
theorem rotB_ray_equivariant (i : Fin 12) :
    SameRay ℝ (evalVec (matVecZ rotBZ (candidateRayZ i)))
      (candidateRay (genB i)) := by
  have hfun : evalVec (matVecZ rotBZ (candidateRayZ i)) =
      (2 : ℝ) • candidateRay (genB i) := by
    funext k
    show evalPhi (matVecZ rotBZ (candidateRayZ i) k) =
      (2 : ℝ) * evalPhi (candidateRayZ (genB i) k)
    rw [candidate_equivariant_rotB i k, evalPhi_zmul, evalPhi_two]
  rw [hfun]
  exact (SameRay.sameRay_nonneg_smul_right _ (by norm_num)).symm

/-- The negated candidate's rays are opposite to the canonical rays:
the two inhabitants induce distinct ray assignments at every port. -/
theorem negated_opposite_ray (i : Fin 12) :
    ¬SameRay ℝ (evalVec (negatedRayZ i)) (candidateRay i) := by
  have h : evalVec (negatedRayZ i) = -candidateRay i := evalVec_vneg _
  rw [h]
  intro hsr
  exact not_sameRay_neg_self (candidateRay_ne_zero i) hsr.symm

end

/-! ## Refinement compatibility at the committed barycentric level

The committed refinement notion (`DiscreteRefinement.refine`) multiplies
same-parent barycentric coordinates by a denominator.  The carrier
extension below maps a barycentric point of a committed face to the exact
integer combination of the face's three vertex rays.  The committed
corpus carries no frequency-`n` incidence complex, so this is the level
at which refinement compatibility is provable; compatibility with a
committed refined complex is open because no such complex is committed. -/

/-- Carrier extension of a port assignment to barycentric points of a
committed face: the exact integer combination of the three vertex rays. -/
def baryCarrierZ (pm : Fin 12 → VecZ) (f : Fin 20) (x : Barycentric) : VecZ :=
  vsc x.i (pm (faceVertices f).1) + vsc x.j (pm (faceVertices f).2.1) +
    vsc x.k (pm (faceVertices f).2.2)

/-- The committed refinement acts on the carrier image by exact scaling:
`refine m` multiplies the image vector by `m`.  This holds for every port
assignment. -/
theorem baryCarrier_refine (pm : Fin 12 → VecZ) (m : ℕ) (f : Fin 20)
    (x : Barycentric) :
    baryCarrierZ pm f (OPH.DiscreteRefinement.refine m x) =
      vsc m (baryCarrierZ pm f x) := by
  funext k
  show zsc (m * x.i) (pm (faceVertices f).1 k) +
      zsc (m * x.j) (pm (faceVertices f).2.1 k) +
      zsc (m * x.k) (pm (faceVertices f).2.2 k) =
    zsc m (zsc x.i (pm (faceVertices f).1 k) +
      zsc x.j (pm (faceVertices f).2.1 k) +
      zsc x.k (pm (faceVertices f).2.2 k))
  rw [zsc_add, zsc_add, zsc_zsc, zsc_zsc, zsc_zsc]

/-- Corner fidelity: the first barycentric corner maps to the first port
ray of the face. -/
theorem baryCarrier_corner_first (pm : Fin 12 → VecZ) (f : Fin 20) :
    baryCarrierZ pm f ⟨1, 0, 0⟩ = pm (faceVertices f).1 := by
  funext k
  show zsc 1 (pm (faceVertices f).1 k) + zsc 0 (pm (faceVertices f).2.1 k) +
      zsc 0 (pm (faceVertices f).2.2 k) = pm (faceVertices f).1 k
  rw [zsc_one, zsc_zero, zsc_zero, add_zero, add_zero]

/-- Corner fidelity for the second corner. -/
theorem baryCarrier_corner_second (pm : Fin 12 → VecZ) (f : Fin 20) :
    baryCarrierZ pm f ⟨0, 1, 0⟩ = pm (faceVertices f).2.1 := by
  funext k
  show zsc 0 (pm (faceVertices f).1 k) + zsc 1 (pm (faceVertices f).2.1 k) +
      zsc 0 (pm (faceVertices f).2.2 k) = pm (faceVertices f).2.1 k
  rw [zsc_one, zsc_zero, zsc_zero, add_zero, zero_add]

/-- Corner fidelity for the third corner. -/
theorem baryCarrier_corner_third (pm : Fin 12 → VecZ) (f : Fin 20) :
    baryCarrierZ pm f ⟨0, 0, 1⟩ = pm (faceVertices f).2.2 := by
  funext k
  show zsc 0 (pm (faceVertices f).1 k) + zsc 0 (pm (faceVertices f).2.1 k) +
      zsc 1 (pm (faceVertices f).2.2 k) = pm (faceVertices f).2.2 k
  rw [zsc_one, zsc_zero, zsc_zero, zero_add, zero_add]

/-- The exact face-center vector of a committed face under the canonical
candidate. -/
def faceCenterZ (f : Fin 20) : VecZ :=
  candidateRayZ (faceVertices f).1 + candidateRayZ (faceVertices f).2.1 +
    candidateRayZ (faceVertices f).2.2

set_option maxHeartbeats 1000000 in
set_option maxRecDepth 8192 in
/-- Each face vertex pairs with its face center to `2 + 3φ`.
Re-derivation: the three ports of a committed face are mutually adjacent,
so the pairing of a vertex with the vertex sum is
`(2 + φ) + φ + φ = 2 + 3φ`. -/
theorem faceCenter_dot_vertices :
    ∀ f : Fin 20,
      dotZ (faceCenterZ f) (candidateRayZ (faceVertices f).1) = (2, 3) ∧
        dotZ (faceCenterZ f) (candidateRayZ (faceVertices f).2.1) = (2, 3) ∧
          dotZ (faceCenterZ f) (candidateRayZ (faceVertices f).2.2) = (2, 3) := by
  decide

/-- The face-center pairing of a barycentric image is the coordinate sum
times `2 + 3φ`. -/
theorem faceCenter_dot_bary (f : Fin 20) (x : Barycentric) :
    dotZ (faceCenterZ f) (baryCarrierZ candidateRayZ f x) =
      zsc (x.i + x.j + x.k) (2, 3) := by
  unfold baryCarrierZ
  rw [dotZ_add_right, dotZ_add_right, dotZ_vsc_right, dotZ_vsc_right,
    dotZ_vsc_right]
  have h := faceCenter_dot_vertices f
  rw [h.1, h.2.1, h.2.2, zsc_add_left, zsc_add_left]

/-- Nondegeneracy: a barycentric point with nonzero coordinate sum maps
to a nonzero exact vector, so every mesh point of every committed face
determines a ray. -/
theorem baryCarrier_ne_zero (f : Fin 20) (x : Barycentric)
    (hx : x.i + x.j + x.k ≠ 0) :
    baryCarrierZ candidateRayZ f x ≠ 0 := by
  intro h
  have hdot := faceCenter_dot_bary f x
  rw [h, dotZ_zero_right] at hdot
  have h1 : (0 : ℤ) = ((x.i + x.j + x.k : ℕ) : ℤ) * 2 :=
    congrArg Prod.fst hdot
  have : x.i + x.j + x.k = 0 := by omega
  exact hx this

noncomputable section

/-- Refinement fixes every mesh ray: the real image of a refined
barycentric point lies on the same ray as the image of the original
point.  This holds for every port assignment and every denominator,
including zero, where the scaled image is the zero vector and `SameRay`
holds degenerately. -/
theorem baryCarrier_refine_sameRay (pm : Fin 12 → VecZ) (m : ℕ) (f : Fin 20)
    (x : Barycentric) :
    SameRay ℝ (evalVec (baryCarrierZ pm f x))
      (evalVec (baryCarrierZ pm f (OPH.DiscreteRefinement.refine m x))) := by
  rw [baryCarrier_refine, evalVec_vsc]
  exact SameRay.sameRay_nonneg_smul_right _ (Nat.cast_nonneg m)

end

end OPH.ScreenCarrierMapCandidate

/- Axiom audit: standard axioms only (`propext`, `Classical.choice`,
`Quot.sound`); no `sorry`, no `native_decide`, no project axiom. -/

#print axioms OPH.ScreenCarrierMapCandidate.zmul_zero
#print axioms OPH.ScreenCarrierMapCandidate.zmul_add
#print axioms OPH.ScreenCarrierMapCandidate.zmul_zsc
#print axioms OPH.ScreenCarrierMapCandidate.zsc_add
#print axioms OPH.ScreenCarrierMapCandidate.zsc_zsc
#print axioms OPH.ScreenCarrierMapCandidate.zsc_one
#print axioms OPH.ScreenCarrierMapCandidate.zsc_zero
#print axioms OPH.ScreenCarrierMapCandidate.zsc_add_left
#print axioms OPH.ScreenCarrierMapCandidate.zsub_eq_zero
#print axioms OPH.ScreenCarrierMapCandidate.dotZ_zero_right
#print axioms OPH.ScreenCarrierMapCandidate.dotZ_add_right
#print axioms OPH.ScreenCarrierMapCandidate.dotZ_vsc_right
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_zmul
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_add
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_zsub
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_zneg
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_zsc
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_phi
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_two
#print axioms OPH.ScreenCarrierMapCandidate.evalVec_vneg
#print axioms OPH.ScreenCarrierMapCandidate.evalVec_vsc
#print axioms OPH.ScreenCarrierMapCandidate.realCross_smul_self
#print axioms OPH.ScreenCarrierMapCandidate.candidate_dot_table
#print axioms OPH.ScreenCarrierMapCandidate.candidate_gram_bridge
#print axioms OPH.ScreenCarrierMapCandidate.candidate_antipode
#print axioms OPH.ScreenCarrierMapCandidate.candidate_coord_ne_zero
#print axioms OPH.ScreenCarrierMapCandidate.candidate_cross_ne_zero
#print axioms OPH.ScreenCarrierMapCandidate.gramTarget_diag
#print axioms OPH.ScreenCarrierMapCandidate.gramTarget_offdiag_ne
#print axioms OPH.ScreenCarrierMapCandidate.gramTarget_phi_iff_adj
#print axioms OPH.ScreenCarrierMapCandidate.genA_listed
#print axioms OPH.ScreenCarrierMapCandidate.genB_listed
#print axioms OPH.ScreenCarrierMapCandidate.rotA_orthogonal
#print axioms OPH.ScreenCarrierMapCandidate.rotB_orthogonal
#print axioms OPH.ScreenCarrierMapCandidate.rotA_det
#print axioms OPH.ScreenCarrierMapCandidate.rotB_det
#print axioms OPH.ScreenCarrierMapCandidate.candidate_equivariant_rotA
#print axioms OPH.ScreenCarrierMapCandidate.candidate_equivariant_rotB
#print axioms OPH.ScreenCarrierMapCandidate.perms_generated_by_two_rotations
#print axioms OPH.ScreenCarrierMapCandidate.PortCarrierCandidate.portMap_injective
#print axioms OPH.ScreenCarrierMapCandidate.PortCarrierCandidate.adjacent_iff
#print axioms OPH.ScreenCarrierMapCandidate.PortCarrierCandidate.seamImage_endpoints
#print axioms OPH.ScreenCarrierMapCandidate.PortCarrierCandidate.seamImage_injective
#print axioms OPH.ScreenCarrierMapCandidate.facePortTriple_injective
#print axioms OPH.ScreenCarrierMapCandidate.PortCarrierCandidate.faceImage_eq_triple
#print axioms OPH.ScreenCarrierMapCandidate.PortCarrierCandidate.faceImage_injective
#print axioms OPH.ScreenCarrierMapCandidate.PortCarrierCandidate.sharesSeam_iff
#print axioms OPH.ScreenCarrierMapCandidate.canonicalCandidate
#print axioms OPH.ScreenCarrierMapCandidate.candidateRayZ_injective
#print axioms OPH.ScreenCarrierMapCandidate.negated_dot_table
#print axioms OPH.ScreenCarrierMapCandidate.negated_equivariant_rotA
#print axioms OPH.ScreenCarrierMapCandidate.negated_equivariant_rotB
#print axioms OPH.ScreenCarrierMapCandidate.negatedCandidate
#print axioms OPH.ScreenCarrierMapCandidate.negated_ne_canonical
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_ne_zero
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_injective
#print axioms OPH.ScreenCarrierMapCandidate.evalVec_injective
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_dotZ
#print axioms OPH.ScreenCarrierMapCandidate.evalPhi_crossZ
#print axioms OPH.ScreenCarrierMapCandidate.candidateRay_ne_zero
#print axioms OPH.ScreenCarrierMapCandidate.candidateRay_not_parallel
#print axioms OPH.ScreenCarrierMapCandidate.not_sameRay_neg_self
#print axioms OPH.ScreenCarrierMapCandidate.candidateRay_antipode
#print axioms OPH.ScreenCarrierMapCandidate.candidateRay_pairwise_distinct_rays
#print axioms OPH.ScreenCarrierMapCandidate.candidateRayHerm_spacelike
#print axioms OPH.ScreenCarrierMapCandidate.candidateRay_adjacent_iff
#print axioms OPH.ScreenCarrierMapCandidate.rotA_ray_equivariant
#print axioms OPH.ScreenCarrierMapCandidate.rotB_ray_equivariant
#print axioms OPH.ScreenCarrierMapCandidate.negated_opposite_ray
#print axioms OPH.ScreenCarrierMapCandidate.baryCarrier_refine
#print axioms OPH.ScreenCarrierMapCandidate.baryCarrier_corner_first
#print axioms OPH.ScreenCarrierMapCandidate.baryCarrier_corner_second
#print axioms OPH.ScreenCarrierMapCandidate.baryCarrier_corner_third
#print axioms OPH.ScreenCarrierMapCandidate.faceCenter_dot_vertices
#print axioms OPH.ScreenCarrierMapCandidate.faceCenter_dot_bary
#print axioms OPH.ScreenCarrierMapCandidate.baryCarrier_ne_zero
#print axioms OPH.ScreenCarrierMapCandidate.baryCarrier_refine_sameRay
