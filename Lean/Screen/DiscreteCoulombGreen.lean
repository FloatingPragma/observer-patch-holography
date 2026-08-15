import SeamU1HolonomyClassification
import ObserverPatchHolography.ScalarSeamRepair

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.DiscreteCoulombGreen

open OPH.SeamCurrentCarrierQuotient
open OPH.SeamU1HolonomyClassification

/-!
# The discrete Green function of the committed seam Laplacian (issue #733)

WHAT IS PROVED.  Premise-free finite mathematics on the committed
thirty-seam table `seamLeft`/`seamRight` of twelve ports, over `ℚ` and
over its `ℝ` scalar extension.  The port coboundary
`(rationalCoboundary φ) e = φ (seamRight e) - φ (seamLeft e)` is the
exact transpose of the committed incidence boundary
`rationalSeamBoundary` for the equal-weight pairings
(`coboundary_boundary_adjoint`).  The vertex Laplacian is the
composition `rationalSeamBoundary ∘ₗ rationalCoboundary`; its matrix is
derived from the committed table alone
(`incidence_product_eq_laplacianEntries`,
`rationalLaplacian_eq_mulVec`), and every port has degree five
(`laplacianEntries_diag_five`).  The explicit rational Green matrix
`greenMatrix` satisfies the exact identities
`laplacianMatrix * greenMatrix = 1 - (1/12) • allOnesQ`
(`laplacian_mul_green`), symmetry (`greenMatrix_symm`), zero row sums
(`greenMatrix_row_sum`), and the reversed product identity
(`green_mul_laplacian`).  The kernel of the Laplacian is exactly the
constant loads, by the energy identity and constancy along the
committed spanning tree `treeSeams` (`laplacian_eq_zero_iff`).  For
every neutral rational load the seam field
`coulombField ρ = rationalCoboundary (greenMatrix.mulVec ρ)` solves the
discrete Gauss problem (`coulombField_gauss`), is orthogonal to every
source-free cycle (`coulombField_orthogonal_cycles`), and is the unique
orthogonal and unique energy-minimizing Gauss solution — the Thomson
variational principle (`thomson_energy_decomposition`,
`thomson_unique_orthogonal`, `thomson_minimum`).  This discharges the
canonical-choice question that the committed `DiscreteGauss` module
explicitly left open: among all Gauss solutions of a neutral load there
is exactly one distinguished by cycle orthogonality, equivalently by
minimal seam energy.  A worked receipt computes the exact potential of
the adjacent-port dipole across seam `0` grouped by the two-step hop
classification of the committed table, which coincides with graph
distance (`tableDistance_*`, `dipole_receipt_values`,
`dipole_value_display`, `dipole_value_strict_ordering`), and the
single-source Green values by distance are `7/36, 1/90, -7/180, -1/18`,
strictly decreasing (`green_single_source`,
`green_distance_strictly_decreasing`).  Finally, the committed uniform
seam-repair operator of `ScalarSeamRepair` acts on real port loads as
`1 - (1/60) • L` for the real scalar extension of this same Laplacian
(`uniformSeamRepair_eq_id_sub_realLaplacian`), and its fixed loads are
exactly the kernel of the Laplacian, exactly the constants, exactly the
loads orthogonal to the whole boundary range — the Gauss-obstruction
directions (`uniformSeamRepair_fixed_iff_constant`,
`constant_iff_orthogonal_to_boundary_range`).

WHAT IS NOT PROVED.  No theorem identifies a rational or real port load
with a physical charge density, a seam current with an electromagnetic
field or flux, the Green matrix with a laboratory potential, or the
seam graph with physical space; those attachments are exactly the open
premise-register rows PR-53 and PR-54 and are not consumed here.  The
word "Coulomb" names the canonical discrete Gauss solution of this
finite table and nothing else.  No continuum limit, no dynamics, no
action principle (see `PositionSpaceMaxwellAction` for the sourced
quadratic form), and no physical readout is stated.  Every theorem in
this file is premise-free exact mathematics; no premise-register row
and no field of any premise bundle is consumed.

Axiom audit.  Every proof composes committed results with finite
`decide` checks and exact rational and real linear algebra; the module
adds no project axiom and uses no native decision procedure.  The audit
lines at the end of the file show at most `propext`, `Classical.choice`,
and `Quot.sound`.
-/

/-! ## Integer incidence data derived from the committed seam table

The incidence entries are defined from `seamLeft`/`seamRight` only.  The
two explicit integer matrices below (`laplacianEntriesZ`, `greenNumZ`)
are pinned to the committed table by the kernel-checked `decide`
theorems that follow them; nothing about the seam graph is assumed. -/

/-- Signed integer incidence of seam `e` at port `p`: the committed
orientation credits `seamRight` and debits `seamLeft`, matching
`rationalDirectedBoundary`. -/
def incidenceZ (e : Fin 30) (p : Fin 12) : ℤ :=
  (if p = seamRight e then 1 else 0) - (if p = seamLeft e then 1 else 0)

/-- Explicit integer Laplacian table.  The next theorem proves it equal
to the incidence product derived from the committed seam table. -/
def laplacianEntriesZ : Matrix (Fin 12) (Fin 12) ℤ :=
  Matrix.of
    ![![5, -1, -1, -1, -1, 0, -1, 0, 0, 0, 0, 0],
      ![-1, 5, -1, -1, 0, -1, 0, -1, 0, 0, 0, 0],
      ![-1, -1, 5, 0, -1, -1, 0, 0, -1, 0, 0, 0],
      ![-1, -1, 0, 5, 0, 0, -1, -1, 0, -1, 0, 0],
      ![-1, 0, -1, 0, 5, 0, -1, 0, -1, 0, -1, 0],
      ![0, -1, -1, 0, 0, 5, 0, -1, -1, 0, 0, -1],
      ![-1, 0, 0, -1, -1, 0, 5, 0, 0, -1, -1, 0],
      ![0, -1, 0, -1, 0, -1, 0, 5, 0, -1, 0, -1],
      ![0, 0, -1, 0, -1, -1, 0, 0, 5, 0, -1, -1],
      ![0, 0, 0, -1, 0, 0, -1, -1, 0, 5, -1, -1],
      ![0, 0, 0, 0, -1, 0, -1, 0, -1, -1, 5, -1],
      ![0, 0, 0, 0, 0, -1, 0, -1, -1, -1, -1, 5]]

/-- One hundred eighty times the Green matrix.  All later rational and
real Green identities are exact casts of the integer facts about this
table. -/
def greenNumZ : Matrix (Fin 12) (Fin 12) ℤ :=
  Matrix.of
    ![![35, 2, 2, 2, 2, -7, 2, -7, -7, -7, -7, -10],
      ![2, 35, 2, 2, -7, 2, -7, 2, -7, -7, -10, -7],
      ![2, 2, 35, -7, 2, 2, -7, -7, 2, -10, -7, -7],
      ![2, 2, -7, 35, -7, -7, 2, 2, -10, 2, -7, -7],
      ![2, -7, 2, -7, 35, -7, 2, -10, 2, -7, 2, -7],
      ![-7, 2, 2, -7, -7, 35, -10, 2, 2, -7, -7, 2],
      ![2, -7, -7, 2, 2, -10, 35, -7, -7, 2, 2, -7],
      ![-7, 2, -7, 2, -10, 2, -7, 35, -7, 2, -7, 2],
      ![-7, -7, 2, -10, 2, 2, -7, -7, 35, -7, 2, 2],
      ![-7, -7, -10, 2, -7, -7, 2, 2, -7, 35, 2, 2],
      ![-7, -10, -7, -7, 2, -7, 2, -7, 2, 2, 35, 2],
      ![-10, -7, -7, -7, -7, 2, -7, 2, 2, 2, 2, 35]]

set_option maxRecDepth 16384 in
/-- The explicit Laplacian table is the incidence product of the
committed seam table: nothing about the graph is assumed. -/
theorem incidence_product_eq_laplacianEntries :
    ∀ p q : Fin 12,
      (∑ e : Fin 30, incidenceZ e p * incidenceZ e q) =
        laplacianEntriesZ p q := by
  decide

set_option maxRecDepth 16384 in
/-- Exact integer core of the Green identity:
`L * (180 G) = 180·1 - 15·(all ones)` entrywise. -/
theorem laplacian_green_product :
    ∀ p q : Fin 12,
      (∑ k : Fin 12, laplacianEntriesZ p k * greenNumZ k q) =
        if p = q then 165 else -15 := by
  decide

/-- The integer Green table is symmetric. -/
theorem greenNumZ_symm : ∀ p q : Fin 12, greenNumZ p q = greenNumZ q p := by
  decide

/-- Every row of the integer Green table sums to zero. -/
theorem greenNumZ_row_sum :
    ∀ p : Fin 12, (∑ q : Fin 12, greenNumZ p q) = 0 := by
  decide

/-- Every port of the committed seam table has degree five. -/
theorem laplacianEntries_diag_five :
    ∀ p : Fin 12, laplacianEntriesZ p p = 5 := by
  decide

/-! ## Hop classification of port pairs from the committed table -/

/-- Two ports bound a common seam of the committed table. -/
def seamAdj (u v : Fin 12) : Prop :=
  ∃ e : Fin 30,
    (seamLeft e = u ∧ seamRight e = v) ∨ (seamLeft e = v ∧ seamRight e = u)

instance : ∀ u v : Fin 12, Decidable (seamAdj u v) := fun u v ↦ by
  unfold seamAdj; infer_instance

/-- Two-step hop classification of an ordered port pair: `0` for equal
ports, `1` for a common seam, `2` for a common neighbour without a
common seam, `3` otherwise.  Classes `0`, `1`, `2` are exact graph
distances by construction; `tableDistance_three_witness` shows class
`3` pairs are joined in three hops, so the classification is the graph
distance of the committed table. -/
def tableDistance (u v : Fin 12) : Fin 4 :=
  if u = v then 0
  else if seamAdj u v then 1
  else if ∃ w : Fin 12, seamAdj u w ∧ seamAdj w v then 2
  else 3

set_option maxRecDepth 16384 in
/-- Every class-`3` pair is one seam plus a class-`2` pair, so class `3`
is exact graph distance three. -/
theorem tableDistance_three_witness :
    ∀ p q : Fin 12, tableDistance p q = 3 →
      ∃ w : Fin 12, seamAdj p w ∧ tableDistance w q = 2 := by
  decide

/-- Integer Green values by hop class. -/
def distGreenNum : Fin 4 → ℤ := ![35, 2, -7, -10]

set_option maxRecDepth 16384 in
/-- Single-source structure of the integer Green table: every entry is
determined by the hop class of its port pair. -/
theorem greenNumZ_distance :
    ∀ p q : Fin 12, greenNumZ p q = distGreenNum (tableDistance p q) := by
  decide

/-- Hop classes of every port relative to port `0`. -/
theorem distFromZero :
    ∀ p : Fin 12,
      tableDistance 0 p = ![0, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3] p := by
  decide

/-- Hop classes of every port relative to port `1`. -/
theorem distFromOne :
    ∀ p : Fin 12,
      tableDistance 1 p = ![1, 0, 1, 1, 2, 1, 2, 1, 2, 2, 3, 2] p := by
  decide

/-- Integer dipole table: `180` times the dipole potential as a function
of the hop classes relative to the two poles.  Unrealized class pairs
carry `0`, which keeps the table antisymmetric. -/
def dipoleNum : Fin 4 → Fin 4 → ℤ :=
  ![![0, 33, 0, 0],
    ![-33, 0, 9, 0],
    ![0, -9, 0, 3],
    ![0, 0, -3, 0]]

set_option maxRecDepth 16384 in
/-- The integer dipole numerators across seam `0` (poles `0` and `1`)
are classified exactly by the hop classes relative to the two poles. -/
theorem dipole_numerator_table :
    ∀ p : Fin 12,
      greenNumZ p 0 - greenNumZ p 1 =
        dipoleNum (tableDistance 0 p) (tableDistance 1 p) := by
  decide

/-- The integer dipole table is antisymmetric under exchanging the two
poles. -/
theorem dipoleNum_antisymm :
    ∀ a b : Fin 4, dipoleNum b a = -dipoleNum a b := by
  decide

/-! ## The rational coboundary and its adjointness -/

/-- Port coboundary on rational loads: the exact transpose of the
committed `rationalSeamBoundary` sign convention, crediting `seamRight`
and debiting `seamLeft`. -/
def rationalCoboundary : RationalPortLoad →ₗ[ℚ] RationalSeamCurrent where
  toFun φ := fun e ↦ φ (seamRight e) - φ (seamLeft e)
  map_add' φ ψ := by
    funext e
    simp only [Pi.add_apply]
    ring
  map_smul' a φ := by
    funext e
    simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
    ring

theorem rationalCoboundary_apply (φ : RationalPortLoad) (e : Fin 30) :
    rationalCoboundary φ e = φ (seamRight e) - φ (seamLeft e) := rfl

/-- Equal-weight pairing of two rational port loads. -/
def portInner (x y : RationalPortLoad) : ℚ := ∑ p : Fin 12, x p * y p

/-- Equal-weight pairing of two rational seam currents. -/
def seamInner (c d : RationalSeamCurrent) : ℚ := ∑ e : Fin 30, c e * d e

/-- Seam energy: the equal-weight sum of squared seam values. -/
def seamEnergy (c : RationalSeamCurrent) : ℚ := ∑ e : Fin 30, c e ^ 2

theorem seamInner_self_eq_energy (c : RationalSeamCurrent) :
    seamInner c c = seamEnergy c :=
  Finset.sum_congr rfl fun e _ ↦ (pow_two (c e)).symm

theorem seamEnergy_nonneg (c : RationalSeamCurrent) : 0 ≤ seamEnergy c :=
  Finset.sum_nonneg fun e _ ↦ sq_nonneg (c e)

theorem seamEnergy_eq_zero_iff (c : RationalSeamCurrent) :
    seamEnergy c = 0 ↔ c = 0 := by
  unfold seamEnergy
  rw [Finset.sum_eq_zero_iff_of_nonneg fun e _ ↦ sq_nonneg (c e)]
  constructor
  · intro h
    funext e
    exact pow_eq_zero_iff two_ne_zero |>.mp (h e (Finset.mem_univ e))
  · rintro rfl e _
    simp

/-- **Exact transposition.**  The coboundary and the committed incidence
boundary are adjoint for the equal-weight pairings; this pins the sign
convention of `rationalCoboundary` to `rationalSeamBoundary`. -/
theorem coboundary_boundary_adjoint (φ : RationalPortLoad)
    (c : RationalSeamCurrent) :
    seamInner (rationalCoboundary φ) c =
      portInner φ (rationalSeamBoundary c) := by
  unfold seamInner portInner
  have hb : ∀ p : Fin 12, rationalSeamBoundary c p =
      ∑ e : Fin 30, rationalDirectedBoundary (seamLeft e) (seamRight e)
        (c e) p := fun p ↦ rfl
  simp only [hb, rationalDirectedBoundary, Finset.mul_sum, mul_sub]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  rw [rationalCoboundary_apply, sub_mul, Finset.sum_sub_distrib]
  congr 1
  · simp only [mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ,
      if_true]
  · simp only [mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ,
      if_true]

/-! ## The vertex Laplacian as a composition and as a matrix -/

/-- The vertex Laplacian: committed boundary after coboundary. -/
def rationalLaplacian : RationalPortLoad →ₗ[ℚ] RationalPortLoad :=
  rationalSeamBoundary ∘ₗ rationalCoboundary

theorem rationalLaplacian_apply (φ : RationalPortLoad) :
    rationalLaplacian φ = rationalSeamBoundary (rationalCoboundary φ) := rfl

/-- Energy identity: the Laplacian pairing is the seam energy of the
coboundary. -/
theorem laplacian_energy_identity (φ : RationalPortLoad) :
    portInner φ (rationalLaplacian φ) =
      seamEnergy (rationalCoboundary φ) := by
  rw [rationalLaplacian_apply, ← coboundary_boundary_adjoint,
    seamInner_self_eq_energy]

/-- Rational incidence entries. -/
def incidenceQ (e : Fin 30) (p : Fin 12) : ℚ :=
  (if p = seamRight e then 1 else 0) - (if p = seamLeft e then 1 else 0)

theorem incidenceQ_eq_cast (e : Fin 30) (p : Fin 12) :
    incidenceQ e p = ((incidenceZ e p : ℤ) : ℚ) := by
  unfold incidenceQ incidenceZ
  split_ifs <;> norm_num

/-- The Laplacian matrix, derived from the committed incidence only. -/
def laplacianMatrix : Matrix (Fin 12) (Fin 12) ℚ :=
  Matrix.of fun p q ↦ ∑ e : Fin 30, incidenceQ e p * incidenceQ e q

theorem laplacianMatrix_apply (p q : Fin 12) :
    laplacianMatrix p q = ∑ e : Fin 30, incidenceQ e p * incidenceQ e q :=
  rfl

theorem laplacianMatrix_eq_cast (p q : Fin 12) :
    laplacianMatrix p q = ((laplacianEntriesZ p q : ℤ) : ℚ) := by
  rw [laplacianMatrix_apply, ← incidence_product_eq_laplacianEntries p q,
    Int.cast_sum]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  rw [Int.cast_mul, incidenceQ_eq_cast, incidenceQ_eq_cast]

theorem incidenceQ_sum_apply (e : Fin 30) (φ : RationalPortLoad) :
    (∑ q : Fin 12, incidenceQ e q * φ q) =
      φ (seamRight e) - φ (seamLeft e) := by
  unfold incidenceQ
  simp only [sub_mul, Finset.sum_sub_distrib, ite_mul, one_mul, zero_mul,
    Finset.sum_ite_eq', Finset.mem_univ, if_true]

/-- The composed Laplacian is multiplication by the incidence-derived
matrix. -/
theorem rationalLaplacian_eq_mulVec (φ : RationalPortLoad) :
    rationalLaplacian φ = laplacianMatrix.mulVec φ := by
  funext p
  have hlhs : rationalLaplacian φ p =
      ∑ e : Fin 30, incidenceQ e p * (φ (seamRight e) - φ (seamLeft e)) := by
    have hb : rationalLaplacian φ p =
        ∑ e : Fin 30, rationalDirectedBoundary (seamLeft e) (seamRight e)
          (rationalCoboundary φ e) p := rfl
    rw [hb]
    refine Finset.sum_congr rfl fun e _ ↦ ?_
    rw [rationalCoboundary_apply]
    unfold rationalDirectedBoundary incidenceQ
    split_ifs <;> ring
  have hmv : laplacianMatrix.mulVec φ p =
      ∑ q : Fin 12, laplacianMatrix p q * φ q := rfl
  rw [hlhs, hmv]
  symm
  calc ∑ q : Fin 12, laplacianMatrix p q * φ q
      = ∑ q : Fin 12, ∑ e : Fin 30, incidenceQ e p * incidenceQ e q * φ q := by
        refine Finset.sum_congr rfl fun q _ ↦ ?_
        rw [laplacianMatrix_apply, Finset.sum_mul]
    _ = ∑ e : Fin 30, ∑ q : Fin 12, incidenceQ e p * incidenceQ e q * φ q :=
        Finset.sum_comm
    _ = ∑ e : Fin 30, incidenceQ e p * (φ (seamRight e) - φ (seamLeft e)) := by
        refine Finset.sum_congr rfl fun e _ ↦ ?_
        rw [← incidenceQ_sum_apply e φ, Finset.mul_sum]
        refine Finset.sum_congr rfl fun q _ ↦ ?_
        ring

/-! ## The explicit rational Green matrix and its identities -/

/-- The rational Green matrix: the integer table over `180`. -/
def greenMatrix : Matrix (Fin 12) (Fin 12) ℚ :=
  Matrix.of fun p q ↦ ((greenNumZ p q : ℤ) : ℚ) / 180

theorem greenMatrix_apply (p q : Fin 12) :
    greenMatrix p q = ((greenNumZ p q : ℤ) : ℚ) / 180 := rfl

/-- The all-ones rational matrix. -/
def allOnesQ : Matrix (Fin 12) (Fin 12) ℚ := Matrix.of fun _ _ ↦ 1

/-- **(G1)** The exact Green identity:
`L * G = 1 - (1/12) • (all ones)`. -/
theorem laplacian_mul_green :
    laplacianMatrix * greenMatrix = 1 - ((1 : ℚ) / 12) • allOnesQ := by
  ext p q
  rw [Matrix.mul_apply]
  have hterm : ∀ k : Fin 12, laplacianMatrix p k * greenMatrix k q =
      ((laplacianEntriesZ p k * greenNumZ k q : ℤ) : ℚ) / 180 := by
    intro k
    rw [laplacianMatrix_eq_cast, greenMatrix_apply, Int.cast_mul]
    ring
  rw [Finset.sum_congr rfl fun k _ ↦ hterm k, ← Finset.sum_div,
    ← Int.cast_sum, laplacian_green_product p q]
  by_cases hpq : p = q <;>
    simp [hpq, Matrix.sub_apply, allOnesQ] <;>
    norm_num

/-- **(G2)** The Green matrix is symmetric. -/
theorem greenMatrix_symm : greenMatrixᵀ = greenMatrix := by
  ext p q
  rw [Matrix.transpose_apply, greenMatrix_apply, greenMatrix_apply,
    greenNumZ_symm q p]

/-- **(G2)** Every row of the Green matrix sums to zero. -/
theorem greenMatrix_row_sum (p : Fin 12) :
    (∑ q : Fin 12, greenMatrix p q) = 0 := by
  simp only [greenMatrix_apply]
  rw [← Finset.sum_div, ← Int.cast_sum, greenNumZ_row_sum p]
  norm_num

/-- The incidence-derived Laplacian matrix is symmetric. -/
theorem laplacianMatrix_symm : laplacianMatrixᵀ = laplacianMatrix := by
  ext p q
  rw [Matrix.transpose_apply, laplacianMatrix_apply, laplacianMatrix_apply]
  exact Finset.sum_congr rfl fun e _ ↦ mul_comm _ _

/-- The reversed exact Green identity, by symmetry of both factors. -/
theorem green_mul_laplacian :
    greenMatrix * laplacianMatrix = 1 - ((1 : ℚ) / 12) • allOnesQ := by
  have h : (laplacianMatrix * greenMatrix)ᵀ =
      greenMatrix * laplacianMatrix := by
    rw [Matrix.transpose_mul, greenMatrix_symm, laplacianMatrix_symm]
  rw [← h, laplacian_mul_green, Matrix.transpose_sub,
    Matrix.transpose_smul, Matrix.transpose_one]
  rfl

theorem allOnesQ_mulVec (x : RationalPortLoad) :
    allOnesQ.mulVec x = fun _ ↦ ∑ q : Fin 12, x q := by
  funext p
  have h : allOnesQ.mulVec x p = ∑ q : Fin 12, allOnesQ p q * x q := rfl
  rw [h]
  simp [allOnesQ]

/-! ## (G3) The kernel of the Laplacian is exactly the constants -/

/-- Constancy along the committed spanning tree `treeSeams`: a load
whose coboundary vanishes on the eleven tree seams is constant.  The
eleven hypotheses extracted below are exactly the tree seams of the
committed `SeamU1HolonomyClassification`, which span all twelve
ports. -/
theorem constant_of_coboundary_vanishing_on_tree (φ : RationalPortLoad)
    (h : ∀ e ∈ treeSeams, rationalCoboundary φ e = 0) :
    ∀ p : Fin 12, φ p = φ 0 := by
  have h0 : φ 1 - φ 0 = 0 := h 0 (by decide)
  have h1 : φ 2 - φ 0 = 0 := h 1 (by decide)
  have h2 : φ 3 - φ 0 = 0 := h 2 (by decide)
  have h3 : φ 4 - φ 0 = 0 := h 3 (by decide)
  have h4 : φ 6 - φ 0 = 0 := h 4 (by decide)
  have h7 : φ 5 - φ 1 = 0 := h 7 (by decide)
  have h8 : φ 7 - φ 1 = 0 := h 8 (by decide)
  have h11 : φ 8 - φ 2 = 0 := h 11 (by decide)
  have h14 : φ 9 - φ 3 = 0 := h 14 (by decide)
  have h17 : φ 10 - φ 4 = 0 := h 17 (by decide)
  have h20 : φ 11 - φ 5 = 0 := h 20 (by decide)
  have e1 : φ 1 = φ 0 := by linarith
  have e2 : φ 2 = φ 0 := by linarith
  have e3 : φ 3 = φ 0 := by linarith
  have e4 : φ 4 = φ 0 := by linarith
  have e5 : φ 5 = φ 0 := by linarith
  have e6 : φ 6 = φ 0 := by linarith
  have e7 : φ 7 = φ 0 := by linarith
  have e8 : φ 8 = φ 0 := by linarith
  have e9 : φ 9 = φ 0 := by linarith
  have e10 : φ 10 = φ 0 := by linarith
  have e11 : φ 11 = φ 0 := by linarith
  intro p
  fin_cases p
  · rfl
  · exact e1
  · exact e2
  · exact e3
  · exact e4
  · exact e5
  · exact e6
  · exact e7
  · exact e8
  · exact e9
  · exact e10
  · exact e11

/-- Vanishing coboundary characterizes the constant loads. -/
theorem rationalCoboundary_eq_zero_iff (φ : RationalPortLoad) :
    rationalCoboundary φ = 0 ↔ ∃ c : ℚ, φ = fun _ ↦ c := by
  constructor
  · intro hd
    exact ⟨φ 0, funext fun p ↦
      constant_of_coboundary_vanishing_on_tree φ
        (fun e _ ↦ congrFun hd e) p⟩
  · rintro ⟨c, rfl⟩
    funext e
    simp [rationalCoboundary_apply]

/-- **(G3)** The kernel of the Laplacian is exactly the constants: the
hard direction uses the energy identity, so a kernel load has vanishing
coboundary, hence is constant along the committed tree. -/
theorem laplacian_eq_zero_iff (φ : RationalPortLoad) :
    rationalLaplacian φ = 0 ↔ ∃ c : ℚ, φ = fun _ ↦ c := by
  constructor
  · intro hL
    have hE : seamEnergy (rationalCoboundary φ) = 0 := by
      rw [← laplacian_energy_identity, hL]
      unfold portInner
      simp
    exact (rationalCoboundary_eq_zero_iff φ).mp
      ((seamEnergy_eq_zero_iff _).mp hE)
  · rintro ⟨c, rfl⟩
    rw [rationalLaplacian_apply,
      (rationalCoboundary_eq_zero_iff _).mpr ⟨c, rfl⟩, map_zero]

/-! ## (G4) The canonical Coulomb solution of the Gauss problem -/

/-- The canonical seam field of a rational port load: coboundary of the
Green potential.  The word Coulomb names this canonical discrete Gauss
solution only. -/
def coulombField (ρ : RationalPortLoad) : RationalSeamCurrent :=
  rationalCoboundary (greenMatrix.mulVec ρ)

/-- **(G4)** For every neutral load, the Coulomb field solves the
discrete Gauss problem of the committed boundary. -/
theorem coulombField_gauss (ρ : RationalPortLoad)
    (hρ : rationalPortTotal ρ = 0) :
    rationalSeamBoundary (coulombField ρ) = ρ := by
  have h1 : rationalSeamBoundary (coulombField ρ) =
      rationalLaplacian (greenMatrix.mulVec ρ) := rfl
  rw [h1, rationalLaplacian_eq_mulVec, Matrix.mulVec_mulVec,
    laplacian_mul_green, Matrix.sub_mulVec, Matrix.one_mulVec,
    Matrix.smul_mulVec, allOnesQ_mulVec]
  have htot : (∑ q : Fin 12, ρ q) = 0 := hρ
  funext p
  simp [htot]

/-! ## (G5) The Thomson variational principle -/

/-- The Coulomb field is orthogonal to every source-free cycle, for any
load: `(d (G ρ), c) = (G ρ, ∂ c) = 0`. -/
theorem coulombField_orthogonal_cycles (ρ : RationalPortLoad)
    (c : RationalSeamCurrent) (hc : rationalSeamBoundary c = 0) :
    seamInner (coulombField ρ) c = 0 := by
  unfold coulombField
  rw [coboundary_boundary_adjoint, hc]
  unfold portInner
  simp

/-- **(G5)** Energy decomposition: every Gauss solution splits as the
Coulomb field plus a cycle, with exactly additive seam energy.  The
cycle fact is the committed `gauss_solution_iff_difference_is_cycle`. -/
theorem thomson_energy_decomposition (ρ : RationalPortLoad)
    (hρ : rationalPortTotal ρ = 0) (E : RationalSeamCurrent)
    (hE : rationalSeamBoundary E = ρ) :
    seamEnergy E =
      seamEnergy (coulombField ρ) + seamEnergy (E - coulombField ρ) := by
  have hFg : rationalSeamBoundary (coulombField ρ) = ρ :=
    coulombField_gauss ρ hρ
  have hker : E - coulombField ρ ∈ LinearMap.ker rationalSeamBoundary :=
    (OPH.DiscreteGauss.gauss_solution_iff_difference_is_cycle hFg).mp hE
  have horth : seamInner (coulombField ρ) (E - coulombField ρ) = 0 :=
    coulombField_orthogonal_cycles ρ _ (LinearMap.mem_ker.mp hker)
  unfold seamEnergy
  have hsplit : ∀ e : Fin 30, E e ^ 2 =
      coulombField ρ e ^ 2 + (E - coulombField ρ) e ^ 2 +
        2 * (coulombField ρ e * (E - coulombField ρ) e) := by
    intro e
    simp only [Pi.sub_apply]
    ring
  rw [Finset.sum_congr rfl fun e _ ↦ hsplit e, Finset.sum_add_distrib,
    Finset.sum_add_distrib, ← Finset.mul_sum]
  have hzero : (∑ e : Fin 30, coulombField ρ e * (E - coulombField ρ) e)
      = 0 := horth
  rw [hzero]
  ring

/-- **(G5)** Uniqueness: the Coulomb field is the only Gauss solution
orthogonal to every source-free cycle. -/
theorem thomson_unique_orthogonal (ρ : RationalPortLoad)
    (hρ : rationalPortTotal ρ = 0) (E : RationalSeamCurrent)
    (hE : rationalSeamBoundary E = ρ)
    (horth : ∀ c : RationalSeamCurrent,
      rationalSeamBoundary c = 0 → seamInner E c = 0) :
    E = coulombField ρ := by
  have hFg : rationalSeamBoundary (coulombField ρ) = ρ :=
    coulombField_gauss ρ hρ
  have hker : E - coulombField ρ ∈ LinearMap.ker rationalSeamBoundary :=
    (OPH.DiscreteGauss.gauss_solution_iff_difference_is_cycle hFg).mp hE
  have hcyc : rationalSeamBoundary (E - coulombField ρ) = 0 :=
    LinearMap.mem_ker.mp hker
  have h1 : seamInner E (E - coulombField ρ) = 0 := horth _ hcyc
  have h2 : seamInner (coulombField ρ) (E - coulombField ρ) = 0 :=
    coulombField_orthogonal_cycles ρ _ hcyc
  have hdiff : seamEnergy (E - coulombField ρ) = 0 := by
    rw [← seamInner_self_eq_energy]
    have hexp : seamInner (E - coulombField ρ) (E - coulombField ρ) =
        seamInner E (E - coulombField ρ) -
          seamInner (coulombField ρ) (E - coulombField ρ) := by
      unfold seamInner
      rw [← Finset.sum_sub_distrib]
      refine Finset.sum_congr rfl fun e _ ↦ ?_
      simp only [Pi.sub_apply]
      ring
    rw [hexp, h1, h2, sub_zero]
  have := (seamEnergy_eq_zero_iff _).mp hdiff
  exact sub_eq_zero.mp this

/-- **(G5)** Minimality: among all Gauss solutions of a neutral load,
the Coulomb field has the strictly least seam energy, with equality
exactly at the Coulomb field itself. -/
theorem thomson_minimum (ρ : RationalPortLoad)
    (hρ : rationalPortTotal ρ = 0) (E : RationalSeamCurrent)
    (hE : rationalSeamBoundary E = ρ) :
    seamEnergy (coulombField ρ) ≤ seamEnergy E ∧
      (seamEnergy E = seamEnergy (coulombField ρ) ↔ E = coulombField ρ) := by
  have hdec := thomson_energy_decomposition ρ hρ E hE
  constructor
  · rw [hdec]
    exact le_add_of_nonneg_right (seamEnergy_nonneg _)
  · constructor
    · intro heq
      have hzero : seamEnergy (E - coulombField ρ) = 0 := by linarith
      exact sub_eq_zero.mp ((seamEnergy_eq_zero_iff _).mp hzero)
    · rintro rfl
      rfl

/-! ## (G6) The exact worked dipole receipt -/

/-- The adjacent-port dipole across seam `0`: unit load at port `0`,
opposite unit load at port `1`. -/
def dipoleLoad : RationalPortLoad := ![1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

theorem dipoleLoad_neutral : rationalPortTotal dipoleLoad = 0 := by
  have h : rationalPortTotal dipoleLoad = ∑ p : Fin 12, dipoleLoad p := rfl
  rw [h]
  simp [dipoleLoad, Fin.sum_univ_succ]

/-- The exact Green potential of the dipole. -/
def dipolePotential : RationalPortLoad := greenMatrix.mulVec dipoleLoad

theorem dipolePotential_apply (p : Fin 12) :
    dipolePotential p = greenMatrix p 0 - greenMatrix p 1 := by
  have h : dipolePotential p = ∑ q : Fin 12, greenMatrix p q * dipoleLoad q :=
    rfl
  rw [h]
  simp [dipoleLoad, Fin.sum_univ_succ]
  ring

/-- Exact rational dipole values by hop-class pair. -/
def dipoleValue (a b : Fin 4) : ℚ := ((dipoleNum a b : ℤ) : ℚ) / 180

/-- **(G6)** The dipole potential at every port is classified exactly by
the hop classes relative to the two poles. -/
theorem dipole_receipt_values (p : Fin 12) :
    dipolePotential p = dipoleValue (tableDistance 0 p) (tableDistance 1 p) := by
  rw [dipolePotential_apply, greenMatrix_apply, greenMatrix_apply,
    div_sub_div_same, ← Int.cast_sub, dipole_numerator_table p]
  rfl

/-- **(G6)** The exact rational dipole values: `11/60` at the positive
pole, `1/20` on the positive near ring, `1/60` at the positive far
port, `0` on the two equatorial classes, and the negatives on the
mirror classes. -/
theorem dipole_value_display :
    dipoleValue 0 1 = 11 / 60 ∧ dipoleValue 1 0 = -(11 / 60) ∧
      dipoleValue 1 2 = 1 / 20 ∧ dipoleValue 2 1 = -(1 / 20) ∧
      dipoleValue 2 3 = 1 / 60 ∧ dipoleValue 3 2 = -(1 / 60) ∧
      dipoleValue 1 1 = 0 ∧ dipoleValue 2 2 = 0 := by
  have h01 : dipoleNum 0 1 = 33 := by decide
  have h10 : dipoleNum 1 0 = -33 := by decide
  have h12 : dipoleNum 1 2 = 9 := by decide
  have h21 : dipoleNum 2 1 = -9 := by decide
  have h23 : dipoleNum 2 3 = 3 := by decide
  have h32 : dipoleNum 3 2 = -3 := by decide
  have h11 : dipoleNum 1 1 = 0 := by decide
  have h22 : dipoleNum 2 2 = 0 := by decide
  unfold dipoleValue
  rw [h01, h10, h12, h21, h23, h32, h11, h22]
  norm_num

/-- **(G6)** The dipole table is antisymmetric under pole exchange. -/
theorem dipoleValue_antisymm (a b : Fin 4) :
    dipoleValue b a = -dipoleValue a b := by
  unfold dipoleValue
  rw [dipoleNum_antisymm a b]
  push_cast
  ring

/-- **(G6)** The exact strict ordering of the dipole values, decreasing
with hop distance from the positive pole. -/
theorem dipole_value_strict_ordering :
    (11 : ℚ) / 60 > 1 / 20 ∧ (1 : ℚ) / 20 > 1 / 60 ∧ (1 : ℚ) / 60 > 0 ∧
      (0 : ℚ) > -(1 / 60) ∧ -(1 / 60 : ℚ) > -(1 / 20) ∧
      -(1 / 20 : ℚ) > -(11 / 60) := by
  norm_num

/-- Exact rational single-source Green values by hop class. -/
def greenDistanceValue (d : Fin 4) : ℚ := ((distGreenNum d : ℤ) : ℚ) / 180

/-- **(G6)** Single-source receipt: every Green entry is the exact
value of its hop class. -/
theorem green_single_source (p q : Fin 12) :
    greenMatrix p q = greenDistanceValue (tableDistance p q) := by
  rw [greenMatrix_apply, greenNumZ_distance p q]
  rfl

/-- **(G6)** The four exact single-source values. -/
theorem greenDistanceValue_display :
    greenDistanceValue 0 = 7 / 36 ∧ greenDistanceValue 1 = 1 / 90 ∧
      greenDistanceValue 2 = -(7 / 180) ∧ greenDistanceValue 3 = -(1 / 18) := by
  have hd0 : distGreenNum 0 = 35 := by decide
  have hd1 : distGreenNum 1 = 2 := by decide
  have hd2 : distGreenNum 2 = -7 := by decide
  have hd3 : distGreenNum 3 = -10 := by decide
  unfold greenDistanceValue
  rw [hd0, hd1, hd2, hd3]
  norm_num

/-- **(G6)** The single-source Green values are strictly decreasing in
hop distance. -/
theorem green_distance_strictly_decreasing :
    greenDistanceValue 0 > greenDistanceValue 1 ∧
      greenDistanceValue 1 > greenDistanceValue 2 ∧
      greenDistanceValue 2 > greenDistanceValue 3 := by
  obtain ⟨h0, h1, h2, h3⟩ := greenDistanceValue_display
  rw [h0, h1, h2, h3]
  norm_num

/-! ## Real scalar extension of the same operators

The real layer repeats the exact rational identities after scalar
extension; every matrix fact is a cast of the same integer tables, so
this is the same committed Laplacian and Green matrix over `ℝ`. -/

noncomputable section RealLayer

/-- Real port coboundary: the same committed sign convention over `ℝ`. -/
def realCoboundary : (Fin 12 → ℝ) →ₗ[ℝ] (Fin 30 → ℝ) where
  toFun φ := fun e ↦ φ (seamRight e) - φ (seamLeft e)
  map_add' φ ψ := by
    funext e
    simp only [Pi.add_apply]
    ring
  map_smul' a φ := by
    funext e
    simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
    ring

theorem realCoboundary_apply (φ : Fin 12 → ℝ) (e : Fin 30) :
    realCoboundary φ e = φ (seamRight e) - φ (seamLeft e) := rfl

/-- Real incidence boundary: the same committed sign convention over
`ℝ`. -/
def realBoundary : (Fin 30 → ℝ) →ₗ[ℝ] (Fin 12 → ℝ) where
  toFun c := fun p ↦ ∑ e : Fin 30,
    ((if p = seamRight e then c e else 0) - (if p = seamLeft e then c e else 0))
  map_add' c d := by
    funext p
    simp only [Pi.add_apply]
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun e _ ↦ ?_
    split_ifs <;> ring
  map_smul' a c := by
    funext p
    simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun e _ ↦ ?_
    split_ifs <;> ring

theorem realBoundary_apply (c : Fin 30 → ℝ) (p : Fin 12) :
    realBoundary c p = ∑ e : Fin 30,
      ((if p = seamRight e then c e else 0) -
        (if p = seamLeft e then c e else 0)) := rfl

/-- The real boundary conserves total load. -/
theorem realBoundary_total (c : Fin 30 → ℝ) :
    (∑ p : Fin 12, realBoundary c p) = 0 := by
  simp only [realBoundary_apply]
  rw [Finset.sum_comm]
  refine Finset.sum_eq_zero fun e _ ↦ ?_
  rw [Finset.sum_sub_distrib]
  simp [Finset.sum_ite_eq', Finset.mem_univ]

/-- Real vertex Laplacian: real boundary after real coboundary. -/
def realLaplacian : (Fin 12 → ℝ) →ₗ[ℝ] (Fin 12 → ℝ) :=
  realBoundary ∘ₗ realCoboundary

theorem realLaplacian_apply (φ : Fin 12 → ℝ) :
    realLaplacian φ = realBoundary (realCoboundary φ) := rfl

/-- Equal-weight pairing of real port loads. -/
def realPortInner (x y : Fin 12 → ℝ) : ℝ := ∑ p : Fin 12, x p * y p

/-- Equal-weight pairing of real seam fields. -/
def realSeamInner (c d : Fin 30 → ℝ) : ℝ := ∑ e : Fin 30, c e * d e

/-- Real seam energy. -/
def realSeamEnergy (c : Fin 30 → ℝ) : ℝ := ∑ e : Fin 30, c e ^ 2

theorem realSeamInner_self_eq_energy (c : Fin 30 → ℝ) :
    realSeamInner c c = realSeamEnergy c :=
  Finset.sum_congr rfl fun e _ ↦ (pow_two (c e)).symm

theorem realSeamEnergy_nonneg (c : Fin 30 → ℝ) : 0 ≤ realSeamEnergy c :=
  Finset.sum_nonneg fun e _ ↦ sq_nonneg (c e)

theorem realSeamEnergy_eq_zero_iff (c : Fin 30 → ℝ) :
    realSeamEnergy c = 0 ↔ c = 0 := by
  unfold realSeamEnergy
  rw [Finset.sum_eq_zero_iff_of_nonneg fun e _ ↦ sq_nonneg (c e)]
  constructor
  · intro h
    funext e
    exact pow_eq_zero_iff two_ne_zero |>.mp (h e (Finset.mem_univ e))
  · rintro rfl e _
    simp

theorem realPortInner_self_eq_zero_iff (x : Fin 12 → ℝ) :
    realPortInner x x = 0 ↔ x = 0 := by
  unfold realPortInner
  rw [Finset.sum_eq_zero_iff_of_nonneg fun p _ ↦ mul_self_nonneg (x p)]
  constructor
  · intro h
    funext p
    exact mul_self_eq_zero.mp (h p (Finset.mem_univ p))
  · rintro rfl p _
    simp

/-- Real adjointness: the same exact transposition over `ℝ`. -/
theorem realCoboundary_boundary_adjoint (φ : Fin 12 → ℝ)
    (c : Fin 30 → ℝ) :
    realSeamInner (realCoboundary φ) c = realPortInner φ (realBoundary c) := by
  unfold realSeamInner realPortInner
  simp only [realBoundary_apply, Finset.mul_sum, mul_sub]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  rw [realCoboundary_apply, sub_mul, Finset.sum_sub_distrib]
  congr 1
  · simp only [mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ,
      if_true]
  · simp only [mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ,
      if_true]

theorem realLaplacian_energy_identity (φ : Fin 12 → ℝ) :
    realPortInner φ (realLaplacian φ) =
      realSeamEnergy (realCoboundary φ) := by
  rw [realLaplacian_apply, ← realCoboundary_boundary_adjoint,
    realSeamInner_self_eq_energy]

/-- Real tree constancy: the same committed tree seams. -/
theorem real_constant_of_coboundary_vanishing_on_tree (φ : Fin 12 → ℝ)
    (h : ∀ e ∈ treeSeams, realCoboundary φ e = 0) :
    ∀ p : Fin 12, φ p = φ 0 := by
  have h0 : φ 1 - φ 0 = 0 := h 0 (by decide)
  have h1 : φ 2 - φ 0 = 0 := h 1 (by decide)
  have h2 : φ 3 - φ 0 = 0 := h 2 (by decide)
  have h3 : φ 4 - φ 0 = 0 := h 3 (by decide)
  have h4 : φ 6 - φ 0 = 0 := h 4 (by decide)
  have h7 : φ 5 - φ 1 = 0 := h 7 (by decide)
  have h8 : φ 7 - φ 1 = 0 := h 8 (by decide)
  have h11 : φ 8 - φ 2 = 0 := h 11 (by decide)
  have h14 : φ 9 - φ 3 = 0 := h 14 (by decide)
  have h17 : φ 10 - φ 4 = 0 := h 17 (by decide)
  have h20 : φ 11 - φ 5 = 0 := h 20 (by decide)
  have e1 : φ 1 = φ 0 := by linarith
  have e2 : φ 2 = φ 0 := by linarith
  have e3 : φ 3 = φ 0 := by linarith
  have e4 : φ 4 = φ 0 := by linarith
  have e5 : φ 5 = φ 0 := by linarith
  have e6 : φ 6 = φ 0 := by linarith
  have e7 : φ 7 = φ 0 := by linarith
  have e8 : φ 8 = φ 0 := by linarith
  have e9 : φ 9 = φ 0 := by linarith
  have e10 : φ 10 = φ 0 := by linarith
  have e11 : φ 11 = φ 0 := by linarith
  intro p
  fin_cases p
  · rfl
  · exact e1
  · exact e2
  · exact e3
  · exact e4
  · exact e5
  · exact e6
  · exact e7
  · exact e8
  · exact e9
  · exact e10
  · exact e11

theorem realCoboundary_eq_zero_iff (φ : Fin 12 → ℝ) :
    realCoboundary φ = 0 ↔ ∃ c : ℝ, φ = fun _ ↦ c := by
  constructor
  · intro hd
    exact ⟨φ 0, funext fun p ↦
      real_constant_of_coboundary_vanishing_on_tree φ
        (fun e _ ↦ congrFun hd e) p⟩
  · rintro ⟨c, rfl⟩
    funext e
    simp [realCoboundary_apply]

/-- The kernel of the real Laplacian is exactly the constants. -/
theorem realLaplacian_eq_zero_iff (φ : Fin 12 → ℝ) :
    realLaplacian φ = 0 ↔ ∃ c : ℝ, φ = fun _ ↦ c := by
  constructor
  · intro hL
    have hE : realSeamEnergy (realCoboundary φ) = 0 := by
      rw [← realLaplacian_energy_identity, hL]
      unfold realPortInner
      simp
    exact (realCoboundary_eq_zero_iff φ).mp
      ((realSeamEnergy_eq_zero_iff _).mp hE)
  · rintro ⟨c, rfl⟩
    rw [realLaplacian_apply,
      (realCoboundary_eq_zero_iff _).mpr ⟨c, rfl⟩, map_zero]

/-! ### Real matrix layer: the same integer tables over `ℝ` -/

/-- Real incidence entries. -/
def incidenceR (e : Fin 30) (p : Fin 12) : ℝ :=
  (if p = seamRight e then 1 else 0) - (if p = seamLeft e then 1 else 0)

theorem incidenceR_eq_cast (e : Fin 30) (p : Fin 12) :
    incidenceR e p = ((incidenceZ e p : ℤ) : ℝ) := by
  unfold incidenceR incidenceZ
  split_ifs <;> norm_num

/-- The real Laplacian matrix, from the committed incidence only. -/
def laplacianMatrixR : Matrix (Fin 12) (Fin 12) ℝ :=
  Matrix.of fun p q ↦ ∑ e : Fin 30, incidenceR e p * incidenceR e q

theorem laplacianMatrixR_apply (p q : Fin 12) :
    laplacianMatrixR p q = ∑ e : Fin 30, incidenceR e p * incidenceR e q :=
  rfl

theorem laplacianMatrixR_eq_cast (p q : Fin 12) :
    laplacianMatrixR p q = ((laplacianEntriesZ p q : ℤ) : ℝ) := by
  rw [laplacianMatrixR_apply, ← incidence_product_eq_laplacianEntries p q,
    Int.cast_sum]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  rw [Int.cast_mul, incidenceR_eq_cast, incidenceR_eq_cast]

theorem incidenceR_sum_apply (e : Fin 30) (φ : Fin 12 → ℝ) :
    (∑ q : Fin 12, incidenceR e q * φ q) =
      φ (seamRight e) - φ (seamLeft e) := by
  unfold incidenceR
  simp only [sub_mul, Finset.sum_sub_distrib, ite_mul, one_mul, zero_mul,
    Finset.sum_ite_eq', Finset.mem_univ, if_true]

/-- The real Laplacian is multiplication by the same matrix over `ℝ`. -/
theorem realLaplacian_eq_mulVec (φ : Fin 12 → ℝ) :
    realLaplacian φ = laplacianMatrixR.mulVec φ := by
  funext p
  have hlhs : realLaplacian φ p =
      ∑ e : Fin 30, incidenceR e p * (φ (seamRight e) - φ (seamLeft e)) := by
    rw [realLaplacian_apply, realBoundary_apply]
    refine Finset.sum_congr rfl fun e _ ↦ ?_
    rw [realCoboundary_apply]
    unfold incidenceR
    split_ifs <;> ring
  have hmv : laplacianMatrixR.mulVec φ p =
      ∑ q : Fin 12, laplacianMatrixR p q * φ q := rfl
  rw [hlhs, hmv]
  symm
  calc ∑ q : Fin 12, laplacianMatrixR p q * φ q
      = ∑ q : Fin 12, ∑ e : Fin 30, incidenceR e p * incidenceR e q * φ q := by
        refine Finset.sum_congr rfl fun q _ ↦ ?_
        rw [laplacianMatrixR_apply, Finset.sum_mul]
    _ = ∑ e : Fin 30, ∑ q : Fin 12, incidenceR e p * incidenceR e q * φ q :=
        Finset.sum_comm
    _ = ∑ e : Fin 30, incidenceR e p * (φ (seamRight e) - φ (seamLeft e)) := by
        refine Finset.sum_congr rfl fun e _ ↦ ?_
        rw [← incidenceR_sum_apply e φ, Finset.mul_sum]
        refine Finset.sum_congr rfl fun q _ ↦ ?_
        ring

/-- The real Green matrix: the same integer table over `180`. -/
def greenMatrixR : Matrix (Fin 12) (Fin 12) ℝ :=
  Matrix.of fun p q ↦ ((greenNumZ p q : ℤ) : ℝ) / 180

theorem greenMatrixR_apply (p q : Fin 12) :
    greenMatrixR p q = ((greenNumZ p q : ℤ) : ℝ) / 180 := rfl

/-- The real Green matrix is the cast of the rational one: the same
matrix. -/
theorem greenMatrixR_eq_cast (p q : Fin 12) :
    greenMatrixR p q = ((greenMatrix p q : ℚ) : ℝ) := by
  rw [greenMatrixR_apply, greenMatrix_apply]
  push_cast
  ring

/-- The all-ones real matrix. -/
def allOnesR : Matrix (Fin 12) (Fin 12) ℝ := Matrix.of fun _ _ ↦ 1

/-- The exact real Green identity. -/
theorem laplacian_mul_green_real :
    laplacianMatrixR * greenMatrixR = 1 - ((1 : ℝ) / 12) • allOnesR := by
  ext p q
  rw [Matrix.mul_apply]
  have hterm : ∀ k : Fin 12, laplacianMatrixR p k * greenMatrixR k q =
      ((laplacianEntriesZ p k * greenNumZ k q : ℤ) : ℝ) / 180 := by
    intro k
    rw [laplacianMatrixR_eq_cast, greenMatrixR_apply, Int.cast_mul]
    ring
  rw [Finset.sum_congr rfl fun k _ ↦ hterm k, ← Finset.sum_div,
    ← Int.cast_sum, laplacian_green_product p q]
  by_cases hpq : p = q <;>
    simp [hpq, Matrix.sub_apply, allOnesR] <;>
    norm_num

theorem greenMatrixR_symm : greenMatrixRᵀ = greenMatrixR := by
  ext p q
  rw [Matrix.transpose_apply, greenMatrixR_apply, greenMatrixR_apply,
    greenNumZ_symm q p]

theorem greenMatrixR_row_sum (p : Fin 12) :
    (∑ q : Fin 12, greenMatrixR p q) = 0 := by
  simp only [greenMatrixR_apply]
  rw [← Finset.sum_div, ← Int.cast_sum, greenNumZ_row_sum p]
  norm_num

theorem laplacianMatrixR_symm : laplacianMatrixRᵀ = laplacianMatrixR := by
  ext p q
  rw [Matrix.transpose_apply, laplacianMatrixR_apply, laplacianMatrixR_apply]
  exact Finset.sum_congr rfl fun e _ ↦ mul_comm _ _

/-- The reversed exact real Green identity. -/
theorem green_mul_laplacian_real :
    greenMatrixR * laplacianMatrixR = 1 - ((1 : ℝ) / 12) • allOnesR := by
  have h : (laplacianMatrixR * greenMatrixR)ᵀ =
      greenMatrixR * laplacianMatrixR := by
    rw [Matrix.transpose_mul, greenMatrixR_symm, laplacianMatrixR_symm]
  rw [← h, laplacian_mul_green_real, Matrix.transpose_sub,
    Matrix.transpose_smul, Matrix.transpose_one]
  rfl

theorem allOnesR_mulVec (x : Fin 12 → ℝ) :
    allOnesR.mulVec x = fun _ ↦ ∑ q : Fin 12, x q := by
  funext p
  have h : allOnesR.mulVec x p = ∑ q : Fin 12, allOnesR p q * x q := rfl
  rw [h]
  simp [allOnesR]

/-- Symmetry of the real Green operator for the equal-weight port
pairing. -/
theorem realPortInner_green_mulVec (x y : Fin 12 → ℝ) :
    realPortInner (greenMatrixR.mulVec x) y =
      realPortInner x (greenMatrixR.mulVec y) := by
  unfold realPortInner
  have hl : ∀ p : Fin 12, greenMatrixR.mulVec x p * y p =
      ∑ q : Fin 12, greenMatrixR p q * x q * y p := by
    intro p
    have h : greenMatrixR.mulVec x p = ∑ q : Fin 12, greenMatrixR p q * x q :=
      rfl
    rw [h, Finset.sum_mul]
  have hr : ∀ p : Fin 12, x p * greenMatrixR.mulVec y p =
      ∑ q : Fin 12, x p * (greenMatrixR p q * y q) := by
    intro p
    have h : greenMatrixR.mulVec y p = ∑ q : Fin 12, greenMatrixR p q * y q :=
      rfl
    rw [h, Finset.mul_sum]
  rw [Finset.sum_congr rfl fun p _ ↦ hl p,
    Finset.sum_congr rfl fun p _ ↦ hr p, Finset.sum_comm]
  refine Finset.sum_congr rfl fun q _ ↦ Finset.sum_congr rfl fun p _ ↦ ?_
  have hsym : greenMatrixR q p = greenMatrixR p q := by
    rw [greenMatrixR_apply, greenMatrixR_apply, greenNumZ_symm q p]
  rw [hsym]
  ring

/-! ## (G7) The repair bridge to the committed uniform seam repair -/

open ObserverPatchHolography.ScalarSeamRepair in
/-- The committed thirty-seam graph Laplacian of `ScalarSeamRepair`,
instantiated on the committed seam table, is exactly the real vertex
Laplacian of this file: the same committed operator. -/
theorem graphLaplacian_eq_realLaplacian :
    graphLaplacian seamLeft seamRight = realLaplacian := by
  refine LinearMap.ext fun x ↦ ?_
  funext p
  have hg : graphLaplacian seamLeft seamRight x p =
      ∑ e : Fin 30, edgeLaplacian (seamLeft e) (seamRight e) x p := by
    unfold graphLaplacian
    rw [LinearMap.sum_apply]
    exact Finset.sum_apply p Finset.univ _
  rw [hg, realLaplacian_apply, realBoundary_apply]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  rw [edgeLaplacian_apply, realCoboundary_apply]
  split_ifs <;> ring

open ObserverPatchHolography.ScalarSeamRepair in
/-- **(G7)** The committed uniform seam repair acts on real port loads
as `1 - (1/60) • L` for the same committed Laplacian.  Boundary: this
is the same committed operator identity; no additional physical claim
is made. -/
theorem uniformSeamRepair_eq_id_sub_realLaplacian :
    uniformSeamRepair seamLeft seamRight =
      LinearMap.id - (1 / 60 : ℝ) • realLaplacian := by
  rw [uniform_thirty_seam_repair, graphLaplacian_eq_realLaplacian]

open ObserverPatchHolography.ScalarSeamRepair in
/-- **(G7)** Fixed loads of the committed uniform repair are exactly
the kernel of the Laplacian. -/
theorem uniformSeamRepair_fixed_iff (x : Fin 12 → ℝ) :
    uniformSeamRepair seamLeft seamRight x = x ↔ realLaplacian x = 0 := by
  rw [uniformSeamRepair_eq_id_sub_realLaplacian]
  simp only [LinearMap.sub_apply, LinearMap.id_apply, LinearMap.smul_apply]
  constructor
  · intro h
    have h60 : (1 / 60 : ℝ) • realLaplacian x = 0 := sub_eq_self.mp h
    rcases smul_eq_zero.mp h60 with h1 | h2
    · norm_num at h1
    · exact h2
  · intro h
    rw [h, smul_zero, sub_zero]

open ObserverPatchHolography.ScalarSeamRepair in
/-- **(G7)** Fixed loads of the committed uniform repair are exactly
the constant loads. -/
theorem uniformSeamRepair_fixed_iff_constant (x : Fin 12 → ℝ) :
    uniformSeamRepair seamLeft seamRight x = x ↔ ∃ c : ℝ, x = fun _ ↦ c :=
  (uniformSeamRepair_fixed_iff x).trans (realLaplacian_eq_zero_iff x)

/-- **(G7)** The kernel of the Laplacian consists exactly of the loads
orthogonal to the entire boundary range: the Gauss-obstruction
directions.  Over `ℚ` the committed range theorem states the same
obstruction as the zero-total hyperplane. -/
theorem constant_iff_orthogonal_to_boundary_range (x : Fin 12 → ℝ) :
    (∀ c : Fin 30 → ℝ, realPortInner x (realBoundary c) = 0) ↔
      realLaplacian x = 0 := by
  constructor
  · intro h
    have hd : realSeamEnergy (realCoboundary x) = 0 := by
      rw [← realSeamInner_self_eq_energy, realCoboundary_boundary_adjoint]
      exact h (realCoboundary x)
    have hzero : realCoboundary x = 0 := (realSeamEnergy_eq_zero_iff _).mp hd
    rw [realLaplacian_apply, hzero, map_zero]
  · intro h c
    obtain ⟨k, hk⟩ := (realLaplacian_eq_zero_iff x).mp h
    have hzero : realCoboundary x = 0 :=
      (realCoboundary_eq_zero_iff x).mpr ⟨k, hk⟩
    rw [← realCoboundary_boundary_adjoint, hzero]
    unfold realSeamInner
    simp

end RealLayer

end OPH.DiscreteCoulombGreen

/- Axiom audit: finite decidable table checks, exact rational and real
linear algebra, and the committed incidence, Gauss, and seam-repair
receipts only.  Expected axioms per line: at most `propext`,
`Classical.choice`, `Quot.sound`.  No native decision procedure is
used. -/

#print axioms OPH.DiscreteCoulombGreen.incidence_product_eq_laplacianEntries
#print axioms OPH.DiscreteCoulombGreen.laplacian_green_product
#print axioms OPH.DiscreteCoulombGreen.greenNumZ_symm
#print axioms OPH.DiscreteCoulombGreen.greenNumZ_row_sum
#print axioms OPH.DiscreteCoulombGreen.laplacianEntries_diag_five
#print axioms OPH.DiscreteCoulombGreen.tableDistance_three_witness
#print axioms OPH.DiscreteCoulombGreen.greenNumZ_distance
#print axioms OPH.DiscreteCoulombGreen.distFromZero
#print axioms OPH.DiscreteCoulombGreen.distFromOne
#print axioms OPH.DiscreteCoulombGreen.dipole_numerator_table
#print axioms OPH.DiscreteCoulombGreen.dipoleNum_antisymm
#print axioms OPH.DiscreteCoulombGreen.coboundary_boundary_adjoint
#print axioms OPH.DiscreteCoulombGreen.laplacian_energy_identity
#print axioms OPH.DiscreteCoulombGreen.rationalLaplacian_eq_mulVec
#print axioms OPH.DiscreteCoulombGreen.laplacian_mul_green
#print axioms OPH.DiscreteCoulombGreen.greenMatrix_symm
#print axioms OPH.DiscreteCoulombGreen.greenMatrix_row_sum
#print axioms OPH.DiscreteCoulombGreen.green_mul_laplacian
#print axioms OPH.DiscreteCoulombGreen.constant_of_coboundary_vanishing_on_tree
#print axioms OPH.DiscreteCoulombGreen.rationalCoboundary_eq_zero_iff
#print axioms OPH.DiscreteCoulombGreen.laplacian_eq_zero_iff
#print axioms OPH.DiscreteCoulombGreen.coulombField_gauss
#print axioms OPH.DiscreteCoulombGreen.coulombField_orthogonal_cycles
#print axioms OPH.DiscreteCoulombGreen.thomson_energy_decomposition
#print axioms OPH.DiscreteCoulombGreen.thomson_unique_orthogonal
#print axioms OPH.DiscreteCoulombGreen.thomson_minimum
#print axioms OPH.DiscreteCoulombGreen.dipoleLoad_neutral
#print axioms OPH.DiscreteCoulombGreen.dipolePotential_apply
#print axioms OPH.DiscreteCoulombGreen.dipole_receipt_values
#print axioms OPH.DiscreteCoulombGreen.dipole_value_display
#print axioms OPH.DiscreteCoulombGreen.dipoleValue_antisymm
#print axioms OPH.DiscreteCoulombGreen.dipole_value_strict_ordering
#print axioms OPH.DiscreteCoulombGreen.green_single_source
#print axioms OPH.DiscreteCoulombGreen.greenDistanceValue_display
#print axioms OPH.DiscreteCoulombGreen.green_distance_strictly_decreasing
#print axioms OPH.DiscreteCoulombGreen.realBoundary_total
#print axioms OPH.DiscreteCoulombGreen.realCoboundary_boundary_adjoint
#print axioms OPH.DiscreteCoulombGreen.realLaplacian_energy_identity
#print axioms OPH.DiscreteCoulombGreen.realCoboundary_eq_zero_iff
#print axioms OPH.DiscreteCoulombGreen.realLaplacian_eq_zero_iff
#print axioms OPH.DiscreteCoulombGreen.realLaplacian_eq_mulVec
#print axioms OPH.DiscreteCoulombGreen.greenMatrixR_eq_cast
#print axioms OPH.DiscreteCoulombGreen.laplacian_mul_green_real
#print axioms OPH.DiscreteCoulombGreen.green_mul_laplacian_real
#print axioms OPH.DiscreteCoulombGreen.greenMatrixR_symm
#print axioms OPH.DiscreteCoulombGreen.greenMatrixR_row_sum
#print axioms OPH.DiscreteCoulombGreen.realPortInner_green_mulVec
#print axioms OPH.DiscreteCoulombGreen.graphLaplacian_eq_realLaplacian
#print axioms OPH.DiscreteCoulombGreen.uniformSeamRepair_eq_id_sub_realLaplacian
#print axioms OPH.DiscreteCoulombGreen.uniformSeamRepair_fixed_iff
#print axioms OPH.DiscreteCoulombGreen.uniformSeamRepair_fixed_iff_constant
#print axioms OPH.DiscreteCoulombGreen.constant_iff_orthogonal_to_boundary_range
