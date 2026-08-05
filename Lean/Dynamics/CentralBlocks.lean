import Mathlib.Algebra.Star.NonUnitalSubalgebra
import Mathlib.Algebra.Star.Prod
import Mathlib.Algebra.Star.StarProjection
import Mathlib.Algebra.Algebra.Prod
import Mathlib.Algebra.GroupWithZero.Idempotent
import Mathlib.Topology.Connected.Clopen
import Dynamics.StoneConverse

/-!
# Central-block decomposition and control

`StoneConverse.lean` closes the private converse on one full finite matrix
block.  This module carries the bounded central-block step of the remaining
gap of issue 679: decomposition along central star projections and full
control of the two-block direct sum.

1. A central star projection of a complex star algebra cuts the algebra
   into two complementary corners, and the cut is a star-algebra
   equivalence onto the direct sum of the corners
   (`centralCutEquiv`).  Star automorphisms carry central star projections
   to central star projections (`IsCentralStarProjection.map`) and carry
   the corner of a projection into the corner of its image
   (`map_mem_centralCorner`), so automorphisms act on blocks by relabeling
   projections.
2. On the two-block algebra `Matrix ι ι ℂ × Matrix κ κ ℂ` the center is
   exactly the span of the two minimal central projections `(1, 0)` and
   `(0, 1)`, which are orthogonal, self-adjoint, and sum to one
   (`twoBlock_mem_center_iff`, `blockProjection_add`,
   `blockProjection_mul`).  Every central idempotent is one of `0`, `1`,
   `(1, 0)`, `(0, 1)` (`twoBlock_central_idempotent_cases`).
3. Every star automorphism of the two-block algebra either fixes both
   minimal central projections or swaps them
   (`starAutomorphism_blockProjectionFst_dichotomy`).  The swap branch is
   realized on the square two-block algebra by the block-swap equivalence
   (`prodSwapStarAlgEquiv`, `prodSwap_blockProjectionFst`), so the
   dichotomy is sharp.
4. A pointwise-continuous real-parameter family of star automorphisms
   through the identity fixes each block at every parameter, because the
   orbit of a minimal central projection is a continuous path inside the
   two-point set of minimal central projections
   (`continuousFlow_fixes_blockProjectionFst`).
5. The two-block Stone converse (`twoBlockStoneConverse_of_continuous`):
   every pointwise-continuous real-parameter star-automorphism group of
   the two-block algebra is blockwise conjugation by the coherent unitary
   families of one time-independent self-adjoint Hamiltonian per block.
   Block fixing is derived, so no block-preservation hypothesis is taken.

The results cover the two-block direct sum exactly.  For arbitrary finite
central-block sums, the spectral decomposition of the center of a general
finite-dimensional star subalgebra into minimal central projections and
the induction over block counts are separate obligations of issue 679.
No physical clock, source-selected Hamiltonian, or physical time
identification is claimed; the parameter is a real number.
-/

namespace OPH.Dynamics

universe u v

/-! ## Central star projections and the corner cut -/

section CentralCut

variable {A : Type u}

/-- A central star projection: a self-adjoint idempotent commuting with
every element of the algebra. -/
structure IsCentralStarProjection [Mul A] [Star A] (p : A) : Prop
    extends IsStarProjection p where
  central : ∀ x : A, p * x = x * p

variable [Ring A] [StarRing A] [Algebra ℂ A] [StarModule ℂ A]

omit [Algebra ℂ A] [StarModule ℂ A] in
/-- The complement of a central star projection is a central star
projection. -/
theorem IsCentralStarProjection.one_sub {p : A}
    (hp : IsCentralStarProjection p) : IsCentralStarProjection (1 - p) where
  toIsStarProjection := hp.toIsStarProjection.one_sub
  central x := by rw [sub_mul, mul_sub, one_mul, mul_one, hp.central x]

omit [StarModule ℂ A] in
/-- A star-algebra automorphism carries a central star projection to a
central star projection. -/
theorem IsCentralStarProjection.map {p : A} (hp : IsCentralStarProjection p)
    (F : A ≃⋆ₐ[ℂ] A) : IsCentralStarProjection (F p) where
  toIsStarProjection := hp.toIsStarProjection.map F
  central y := by
    have hy : F (F.symm y) = y := F.apply_symm_apply y
    calc F p * y = F p * F (F.symm y) := by rw [hy]
      _ = F (p * F.symm y) := (map_mul F _ _).symm
      _ = F (F.symm y * p) := by rw [hp.central]
      _ = F (F.symm y) * F p := map_mul F _ _
      _ = y * F p := by rw [hy]

omit [Algebra ℂ A] [StarModule ℂ A] in
/-- Corner multiplicativity of a central star projection: cutting is
compatible with products. -/
theorem IsCentralStarProjection.mul_cut {p : A}
    (hp : IsCentralStarProjection p) (x y : A) :
    (p * x) * (p * y) = p * (x * y) := by
  calc (p * x) * (p * y) = p * (x * (p * y)) := by rw [mul_assoc]
    _ = p * ((x * p) * y) := by rw [← mul_assoc x p y]
    _ = p * ((p * x) * y) := by rw [← hp.central x]
    _ = p * (p * (x * y)) := by rw [mul_assoc p x y]
    _ = (p * p) * (x * y) := by rw [← mul_assoc]
    _ = p * (x * y) := by rw [hp.isIdempotentElem.eq]

/-- The corner of the algebra cut by a central star projection, as a
non-unital star subalgebra.  Its multiplicative unit is the projection
itself, so the ambient unit is outside the proper corner. -/
def centralCorner (p : A) (hp : IsCentralStarProjection p) :
    NonUnitalStarSubalgebra ℂ A where
  carrier := {x : A | p * x = x}
  add_mem' := fun {a b} ha hb => by
    rw [Set.mem_setOf_eq] at ha hb ⊢
    rw [mul_add, ha, hb]
  zero_mem' := by
    rw [Set.mem_setOf_eq, mul_zero]
  mul_mem' := fun {a b} ha hb => by
    rw [Set.mem_setOf_eq] at ha hb ⊢
    rw [← mul_assoc, ha]
  smul_mem' := fun c {x} hx => by
    rw [Set.mem_setOf_eq] at hx ⊢
    rw [mul_smul_comm, hx]
  star_mem' := fun {a} ha => by
    rw [Set.mem_setOf_eq] at ha ⊢
    have h1 : star a * p = star a := by
      have h2 := congrArg star ha
      rwa [star_mul, hp.isSelfAdjoint.star_eq] at h2
    rw [hp.central (star a), h1]

omit [StarModule ℂ A] in
@[simp] theorem mem_centralCorner_iff {p : A} (hp : IsCentralStarProjection p)
    {x : A} : x ∈ centralCorner p hp ↔ p * x = x := Iff.rfl

omit [StarModule ℂ A] in
/-- Cutting any element by a central star projection lands in its
corner. -/
theorem IsCentralStarProjection.mul_mem_corner {p : A}
    (hp : IsCentralStarProjection p) (x : A) :
    p * x ∈ centralCorner p hp := by
  rw [mem_centralCorner_iff, ← mul_assoc, hp.isIdempotentElem.eq]

omit [StarModule ℂ A] in
/-- A star-algebra automorphism carries the corner of a central star
projection into the corner of the image projection: blocks travel with
their projections. -/
theorem map_mem_centralCorner {p : A} (hp : IsCentralStarProjection p)
    (F : A ≃⋆ₐ[ℂ] A) {x : A} (hx : x ∈ centralCorner p hp) :
    F x ∈ centralCorner (F p) (hp.map F) := by
  rw [mem_centralCorner_iff] at hx ⊢
  rw [← map_mul, hx]

/-- Central-cut factorization: a central star projection factors the
algebra as the direct sum of the two complementary corners, through a
star-algebra equivalence. -/
def centralCutEquiv (p : A) (hp : IsCentralStarProjection p) :
    A ≃⋆ₐ[ℂ] (centralCorner p hp × centralCorner (1 - p) hp.one_sub) where
  toFun x := (⟨p * x, hp.mul_mem_corner x⟩,
    ⟨(1 - p) * x, hp.one_sub.mul_mem_corner x⟩)
  invFun uv := (uv.1 : A) + (uv.2 : A)
  left_inv x := by
    show (p * x) + ((1 - p) * x) = x
    rw [← add_mul]
    have h : p + (1 - p) = (1 : A) := by abel
    rw [h, one_mul]
  right_inv uv := by
    obtain ⟨⟨u, hu⟩, ⟨v, hv⟩⟩ := uv
    rw [mem_centralCorner_iff] at hu hv
    have hpv : p * v = 0 := by
      have h0 : p * (1 - p) = 0 := by
        rw [mul_sub, mul_one, hp.isIdempotentElem.eq, sub_self]
      rw [← hv, ← mul_assoc, h0, zero_mul]
    have hqu : (1 - p) * u = 0 := by
      rw [sub_mul, one_mul, hu, sub_self]
    refine Prod.ext (Subtype.ext ?_) (Subtype.ext ?_)
    · show p * (u + v) = u
      rw [mul_add, hu, hpv, add_zero]
    · show (1 - p) * (u + v) = v
      rw [mul_add, hqu, hv, zero_add]
  map_mul' x y := by
    refine Prod.ext (Subtype.ext ?_) (Subtype.ext ?_)
    · show p * (x * y) = (p * x) * (p * y)
      exact (hp.mul_cut x y).symm
    · show (1 - p) * (x * y) = ((1 - p) * x) * ((1 - p) * y)
      exact (hp.one_sub.mul_cut x y).symm
  map_add' x y := by
    refine Prod.ext (Subtype.ext ?_) (Subtype.ext ?_)
    · show p * (x + y) = p * x + p * y
      exact mul_add p x y
    · show (1 - p) * (x + y) = (1 - p) * x + (1 - p) * y
      exact mul_add _ x y
  map_smul' c x := by
    refine Prod.ext (Subtype.ext ?_) (Subtype.ext ?_)
    · show p * (c • x) = c • (p * x)
      exact mul_smul_comm c p x
    · show (1 - p) * (c • x) = c • ((1 - p) * x)
      exact mul_smul_comm c _ x
  map_star' x := by
    refine Prod.ext (Subtype.ext ?_) (Subtype.ext ?_)
    · show p * star x = star (p * x)
      rw [star_mul, hp.isSelfAdjoint.star_eq]
      exact hp.central (star x)
    · show (1 - p) * star x = star ((1 - p) * x)
      rw [star_mul, hp.one_sub.isSelfAdjoint.star_eq]
      exact hp.one_sub.central (star x)

end CentralCut

/-! ## The two-block algebra and its minimal central projections -/

section TwoBlockDefs

variable (ι : Type u) (κ : Type v) [Fintype ι] [DecidableEq ι]
  [Fintype κ] [DecidableEq κ]

/-- The two-block algebra: the direct sum of two full finite matrix
blocks. -/
abbrev TwoBlockAlgebra := Matrix ι ι ℂ × Matrix κ κ ℂ

/-- The minimal central projection carried by the first block. -/
def blockProjectionFst : TwoBlockAlgebra ι κ := (1, 0)

/-- The minimal central projection carried by the second block. -/
def blockProjectionSnd : TwoBlockAlgebra ι κ := (0, 1)

omit [DecidableEq κ] in
/-- The first block projection is a central star projection. -/
theorem isCentralStarProjection_blockProjectionFst :
    IsCentralStarProjection (blockProjectionFst ι κ) where
  toIsStarProjection :=
    { isIdempotentElem := by
        show ((1, 0) : TwoBlockAlgebra ι κ) * (1, 0) = (1, 0)
        exact Prod.ext (one_mul 1) (mul_zero 0)
      isSelfAdjoint := by
        show star ((1, 0) : TwoBlockAlgebra ι κ) = (1, 0)
        exact Prod.ext (star_one _) (star_zero _) }
  central x := by
    refine Prod.ext ?_ ?_
    · show 1 * x.1 = x.1 * 1
      rw [one_mul, mul_one]
    · show 0 * x.2 = x.2 * 0
      rw [zero_mul, mul_zero]

omit [DecidableEq ι] in
/-- The second block projection is a central star projection. -/
theorem isCentralStarProjection_blockProjectionSnd :
    IsCentralStarProjection (blockProjectionSnd ι κ) where
  toIsStarProjection :=
    { isIdempotentElem := by
        show ((0, 1) : TwoBlockAlgebra ι κ) * (0, 1) = (0, 1)
        exact Prod.ext (mul_zero 0) (one_mul 1)
      isSelfAdjoint := by
        show star ((0, 1) : TwoBlockAlgebra ι κ) = (0, 1)
        exact Prod.ext (star_zero _) (star_one _) }
  central x := by
    refine Prod.ext ?_ ?_
    · show 0 * x.1 = x.1 * 0
      rw [zero_mul, mul_zero]
    · show 1 * x.2 = x.2 * 1
      rw [one_mul, mul_one]

omit [Fintype ι] [Fintype κ] in
/-- The two block projections sum to the unit. -/
theorem blockProjection_add :
    blockProjectionFst ι κ + blockProjectionSnd ι κ = 1 := by
  refine Prod.ext ?_ ?_
  · show (1 : Matrix ι ι ℂ) + 0 = 1
    rw [add_zero]
  · show (0 : Matrix κ κ ℂ) + 1 = 1
    rw [zero_add]

/-- The two block projections are orthogonal. -/
theorem blockProjection_mul :
    blockProjectionFst ι κ * blockProjectionSnd ι κ = 0 := by
  refine Prod.ext ?_ ?_
  · show (1 : Matrix ι ι ℂ) * 0 = 0
    rw [mul_zero]
  · show (0 : Matrix κ κ ℂ) * 1 = 0
    rw [zero_mul]

omit [Fintype ι] [Fintype κ] in
/-- The scalar span of the block projections in pair form. -/
theorem blockCombination_eq_pair (a b : ℂ) :
    a • blockProjectionFst ι κ + b • blockProjectionSnd ι κ =
      ((a • 1 : Matrix ι ι ℂ), (b • 1 : Matrix κ κ ℂ)) := by
  refine Prod.ext ?_ ?_
  · show a • (1 : Matrix ι ι ℂ) + b • 0 = a • 1
    rw [smul_zero, add_zero]
  · show a • (0 : Matrix κ κ ℂ) + b • 1 = b • 1
    rw [smul_zero, zero_add]

end TwoBlockDefs

/-! ## Center classification of the two-block algebra -/

section TwoBlockCenter

variable {ι : Type u} {κ : Type v} [Fintype ι] [DecidableEq ι]
  [Fintype κ] [DecidableEq κ]

private theorem scalar_eq_smul_one {n : Type*} [Fintype n] [DecidableEq n]
    (c : ℂ) : Matrix.scalar n c = c • (1 : Matrix n n ℂ) := by
  rw [Matrix.scalar_apply, Matrix.smul_one_eq_diagonal]

/-- The center of the two-block algebra is the scalar span of the two
minimal central projections. -/
theorem twoBlock_mem_center_iff {z : TwoBlockAlgebra ι κ} :
    (∀ x, z * x = x * z) ↔
      ∃ a b : ℂ,
        z = a • blockProjectionFst ι κ + b • blockProjectionSnd ι κ := by
  constructor
  · intro hz
    have h1 : ∀ m : Matrix ι ι ℂ, z.1 * m = m * z.1 := fun m =>
      congrArg Prod.fst (hz (m, 0))
    have h2 : ∀ m : Matrix κ κ ℂ, z.2 * m = m * z.2 := fun m =>
      congrArg Prod.snd (hz (0, m))
    obtain ⟨a, ha⟩ := Matrix.mem_range_scalar_of_commute_single
      (M := z.1) (fun i j _ => (h1 (Matrix.single i j 1)).symm)
    obtain ⟨b, hb⟩ := Matrix.mem_range_scalar_of_commute_single
      (M := z.2) (fun i j _ => (h2 (Matrix.single i j 1)).symm)
    refine ⟨a, b, ?_⟩
    rw [blockCombination_eq_pair]
    refine Prod.ext ?_ ?_
    · show z.1 = a • 1
      rw [← scalar_eq_smul_one, ha]
    · show z.2 = b • 1
      rw [← scalar_eq_smul_one, hb]
  · rintro ⟨a, b, rfl⟩ x
    rw [blockCombination_eq_pair]
    refine Prod.ext ?_ ?_
    · show (a • (1 : Matrix ι ι ℂ)) * x.1 = x.1 * (a • 1)
      rw [smul_mul_assoc, one_mul, mul_smul_comm, mul_one]
    · show (b • (1 : Matrix κ κ ℂ)) * x.2 = x.2 * (b • 1)
      rw [smul_mul_assoc, one_mul, mul_smul_comm, mul_one]

private theorem smul_one_idempotent_scalar {n : Type u} [Fintype n]
    [DecidableEq n] [Nonempty n] {a : ℂ}
    (h : (a • 1 : Matrix n n ℂ) * (a • 1) = a • 1) : a = 0 ∨ a = 1 := by
  have h2 : (a * a) • (1 : Matrix n n ℂ) = a • 1 := by
    rw [← smul_smul]
    calc a • (a • (1 : Matrix n n ℂ))
        = a • ((1 : Matrix n n ℂ) * (a • 1)) := by rw [one_mul]
      _ = (a • (1 : Matrix n n ℂ)) * (a • 1) :=
          (smul_mul_assoc a 1 (a • 1)).symm
      _ = a • 1 := h
  obtain ⟨i⟩ := ‹Nonempty n›
  have h3 := congrArg (fun m : Matrix n n ℂ => m i i) h2
  have haa : a * a = a := by
    simpa [Matrix.smul_apply, Matrix.one_apply_eq] using h3
  exact IsIdempotentElem.iff_eq_zero_or_one.mp haa

/-- Every central idempotent of the two-block algebra is zero, one, or one
of the two minimal central projections. -/
theorem twoBlock_central_idempotent_cases [Nonempty ι] [Nonempty κ]
    {z : TwoBlockAlgebra ι κ} (hcentral : ∀ x, z * x = x * z)
    (hidem : IsIdempotentElem z) :
    z = 0 ∨ z = 1 ∨ z = blockProjectionFst ι κ ∨
      z = blockProjectionSnd ι κ := by
  obtain ⟨a, b, hz⟩ := twoBlock_mem_center_iff.mp hcentral
  rw [blockCombination_eq_pair] at hz
  subst hz
  have h1 : (a • (1 : Matrix ι ι ℂ)) * (a • 1) = a • 1 :=
    congrArg Prod.fst hidem
  have h2 : (b • (1 : Matrix κ κ ℂ)) * (b • 1) = b • 1 :=
    congrArg Prod.snd hidem
  obtain ha | ha := smul_one_idempotent_scalar h1 <;>
    obtain hb | hb := smul_one_idempotent_scalar h2 <;> subst ha <;> subst hb
  · exact Or.inl (Prod.ext (zero_smul ℂ 1) (zero_smul ℂ 1))
  · exact Or.inr (Or.inr (Or.inr
      (Prod.ext (zero_smul ℂ 1) (one_smul ℂ 1))))
  · exact Or.inr (Or.inr (Or.inl
      (Prod.ext (one_smul ℂ 1) (zero_smul ℂ 1))))
  · exact Or.inr (Or.inl (Prod.ext (one_smul ℂ 1) (one_smul ℂ 1)))

omit [Fintype ι] [Fintype κ] [DecidableEq κ] in
theorem blockProjectionFst_ne_zero [Nonempty ι] :
    blockProjectionFst ι κ ≠ 0 := by
  intro h
  obtain ⟨i⟩ := ‹Nonempty ι›
  have h2 := congrArg (fun z : TwoBlockAlgebra ι κ => z.1 i i) h
  simp [blockProjectionFst, Matrix.one_apply_eq] at h2

omit [Fintype ι] [Fintype κ] in
theorem blockProjectionFst_ne_one [Nonempty κ] :
    blockProjectionFst ι κ ≠ 1 := by
  intro h
  obtain ⟨j⟩ := ‹Nonempty κ›
  have h2 := congrArg (fun z : TwoBlockAlgebra ι κ => z.2 j j) h
  simp [blockProjectionFst, Matrix.one_apply_eq] at h2

omit [Fintype ι] [Fintype κ] in
theorem blockProjectionFst_ne_blockProjectionSnd [Nonempty ι] :
    blockProjectionFst ι κ ≠ blockProjectionSnd ι κ := by
  intro h
  obtain ⟨i⟩ := ‹Nonempty ι›
  have h2 := congrArg (fun z : TwoBlockAlgebra ι κ => z.1 i i) h
  simp [blockProjectionFst, blockProjectionSnd,
    Matrix.one_apply_eq] at h2

end TwoBlockCenter

/-! ## Block permutation dichotomy and continuous block fixing -/

section TwoBlockControl

variable {ι : Type u} {κ : Type v} [Fintype ι] [DecidableEq ι]
  [Fintype κ] [DecidableEq κ]

/-- Every star automorphism of the two-block algebra fixes the first block
projection or carries it to the second: single automorphisms act on the
block labels by a permutation. -/
theorem starAutomorphism_blockProjectionFst_dichotomy [Nonempty ι]
    [Nonempty κ]
    (F : TwoBlockAlgebra ι κ ≃⋆ₐ[ℂ] TwoBlockAlgebra ι κ) :
    F (blockProjectionFst ι κ) = blockProjectionFst ι κ ∨
      F (blockProjectionFst ι κ) = blockProjectionSnd ι κ := by
  have hp := (isCentralStarProjection_blockProjectionFst ι κ).map F
  obtain h0 | h1 | he₁ | he₂ :=
    twoBlock_central_idempotent_cases hp.central hp.isIdempotentElem
  · exfalso
    have h2 : blockProjectionFst ι κ = 0 := by
      have h3 := congrArg F.symm h0
      rwa [F.symm_apply_apply, map_zero] at h3
    exact blockProjectionFst_ne_zero h2
  · exfalso
    have h2 : blockProjectionFst ι κ = 1 := by
      have h3 := congrArg F.symm h1
      rwa [F.symm_apply_apply, map_one] at h3
    exact blockProjectionFst_ne_one h2
  · exact Or.inl he₁
  · exact Or.inr he₂

/-- The block-swap star-algebra equivalence of a product of two star
algebras. -/
def prodSwapStarAlgEquiv (M : Type u) (N : Type v)
    [Add M] [Mul M] [SMul ℂ M] [Star M]
    [Add N] [Mul N] [SMul ℂ N] [Star N] :
    (M × N) ≃⋆ₐ[ℂ] (N × M) where
  toFun := Prod.swap
  invFun := Prod.swap
  left_inv := Prod.swap_swap
  right_inv := Prod.swap_swap
  map_mul' _ _ := rfl
  map_add' _ _ := rfl
  map_smul' _ _ := rfl
  map_star' _ := rfl

/-- On the square two-block algebra the swap automorphism realizes the
permuting branch of the block dichotomy: it carries the first block
projection to the second, so single automorphisms need not fix blocks. -/
theorem prodSwap_blockProjectionFst (n : Type u) [Fintype n]
    [DecidableEq n] :
    prodSwapStarAlgEquiv (Matrix n n ℂ) (Matrix n n ℂ)
      (blockProjectionFst n n) = blockProjectionSnd n n := rfl

/-- A pointwise-continuous family of star automorphisms through the
identity fixes the first block projection at every parameter: the orbit is
a continuous path inside the two-point set of minimal central projections,
and the parameter line is connected. -/
theorem continuousFlow_fixes_blockProjectionFst [Nonempty ι] [Nonempty κ]
    (α : ℝ → (TwoBlockAlgebra ι κ ≃⋆ₐ[ℂ] TwoBlockAlgebra ι κ))
    (hzero : α 0 (blockProjectionFst ι κ) = blockProjectionFst ι κ)
    (hcont : Continuous fun t : ℝ => α t (blockProjectionFst ι κ)) :
    ∀ t : ℝ, α t (blockProjectionFst ι κ) = blockProjectionFst ι κ := by
  set f : ℝ → TwoBlockAlgebra ι κ :=
    fun t => α t (blockProjectionFst ι κ) with hf
  have hdich : ∀ t, f t = blockProjectionFst ι κ ∨
      f t = blockProjectionSnd ι κ :=
    fun t => starAutomorphism_blockProjectionFst_dichotomy (α t)
  have hS : IsClosed (f ⁻¹' {blockProjectionFst ι κ}) :=
    isClosed_singleton.preimage hcont
  have hT : IsClosed (f ⁻¹' {blockProjectionSnd ι κ}) :=
    isClosed_singleton.preimage hcont
  have hcompl : (f ⁻¹' {blockProjectionFst ι κ})ᶜ =
      f ⁻¹' {blockProjectionSnd ι κ} := by
    ext t
    simp only [Set.mem_compl_iff, Set.mem_preimage, Set.mem_singleton_iff]
    constructor
    · intro hne
      exact (hdich t).resolve_left hne
    · intro he h1
      exact blockProjectionFst_ne_blockProjectionSnd (h1.symm.trans he)
  have hclopen : IsClopen (f ⁻¹' {blockProjectionFst ι κ}) := by
    refine ⟨hS, ?_⟩
    rw [← hcompl] at hT
    exact isClosed_compl_iff.mp hT
  obtain hempty | huniv := isClopen_iff.mp hclopen
  · exfalso
    have h0 : (0 : ℝ) ∈ f ⁻¹' {blockProjectionFst ι κ} := by
      rw [Set.mem_preimage, Set.mem_singleton_iff]
      exact hzero
    rw [hempty] at h0
    exact Set.notMem_empty 0 h0
  · intro t
    have ht : t ∈ f ⁻¹' {blockProjectionFst ι κ} :=
      huniv ▸ Set.mem_univ t
    rwa [Set.mem_preimage, Set.mem_singleton_iff] at ht

end TwoBlockControl

/-! ## Block restriction of block-fixing automorphisms -/

section BlockRestriction

variable {ι : Type u} {κ : Type v} [Fintype ι] [DecidableEq ι]
  [Fintype κ] [DecidableEq κ]

omit [DecidableEq κ] in
/-- A block-fixing automorphism carries the first block into itself. -/
theorem apply_fst_pair
    (F : TwoBlockAlgebra ι κ ≃⋆ₐ[ℂ] TwoBlockAlgebra ι κ)
    (hF : F (blockProjectionFst ι κ) = blockProjectionFst ι κ)
    (x : Matrix ι ι ℂ) :
    F (x, 0) = ((F (x, 0)).1, 0) := by
  have hsplit : ((x, 0) : TwoBlockAlgebra ι κ) =
      (x, 0) * blockProjectionFst ι κ := by
    refine Prod.ext ?_ ?_
    · show x = x * 1
      rw [mul_one]
    · show (0 : Matrix κ κ ℂ) = 0 * 0
      rw [mul_zero]
  calc F (x, 0) = F ((x, 0) * blockProjectionFst ι κ) := by rw [← hsplit]
    _ = F (x, 0) * blockProjectionFst ι κ := by rw [map_mul, hF]
    _ = ((F (x, 0)).1, 0) := by
        refine Prod.ext ?_ ?_
        · show (F (x, 0)).1 * 1 = (F (x, 0)).1
          rw [mul_one]
        · show (F (x, 0)).2 * 0 = 0
          rw [mul_zero]

omit [DecidableEq ι] in
/-- A second-block-fixing automorphism carries the second block into
itself. -/
theorem apply_snd_pair
    (F : TwoBlockAlgebra ι κ ≃⋆ₐ[ℂ] TwoBlockAlgebra ι κ)
    (hF : F (blockProjectionSnd ι κ) = blockProjectionSnd ι κ)
    (y : Matrix κ κ ℂ) :
    F (0, y) = (0, (F (0, y)).2) := by
  have hsplit : ((0, y) : TwoBlockAlgebra ι κ) =
      (0, y) * blockProjectionSnd ι κ := by
    refine Prod.ext ?_ ?_
    · show (0 : Matrix ι ι ℂ) = 0 * 0
      rw [mul_zero]
    · show y = y * 1
      rw [mul_one]
  calc F (0, y) = F ((0, y) * blockProjectionSnd ι κ) := by rw [← hsplit]
    _ = F (0, y) * blockProjectionSnd ι κ := by rw [map_mul, hF]
    _ = (0, (F (0, y)).2) := by
        refine Prod.ext ?_ ?_
        · show (F (0, y)).1 * 0 = 0
          rw [mul_zero]
        · show (F (0, y)).2 * 1 = (F (0, y)).2
          rw [mul_one]

omit [DecidableEq κ] in
/-- Fixing the first block projection passes to the inverse
automorphism. -/
theorem symm_fixes_blockProjectionFst
    (F : TwoBlockAlgebra ι κ ≃⋆ₐ[ℂ] TwoBlockAlgebra ι κ)
    (hF : F (blockProjectionFst ι κ) = blockProjectionFst ι κ) :
    F.symm (blockProjectionFst ι κ) = blockProjectionFst ι κ := by
  conv_lhs => rw [← hF]
  rw [F.symm_apply_apply]

omit [DecidableEq ι] in
/-- Fixing the second block projection passes to the inverse
automorphism. -/
theorem symm_fixes_blockProjectionSnd
    (F : TwoBlockAlgebra ι κ ≃⋆ₐ[ℂ] TwoBlockAlgebra ι κ)
    (hF : F (blockProjectionSnd ι κ) = blockProjectionSnd ι κ) :
    F.symm (blockProjectionSnd ι κ) = blockProjectionSnd ι κ := by
  conv_lhs => rw [← hF]
  rw [F.symm_apply_apply]

/-- The first-block restriction of a block-fixing star automorphism of the
two-block algebra. -/
def blockRestrictionFst
    (F : TwoBlockAlgebra ι κ ≃⋆ₐ[ℂ] TwoBlockAlgebra ι κ)
    (hF : F (blockProjectionFst ι κ) = blockProjectionFst ι κ) :
    Matrix ι ι ℂ ≃⋆ₐ[ℂ] Matrix ι ι ℂ where
  toFun x := (F (x, 0)).1
  invFun x := (F.symm (x, 0)).1
  left_inv x := by
    show (F.symm ((F (x, 0)).1, 0)).1 = x
    rw [← apply_fst_pair F hF x, F.symm_apply_apply]
  right_inv x := by
    show (F ((F.symm (x, 0)).1, 0)).1 = x
    rw [← apply_fst_pair F.symm (symm_fixes_blockProjectionFst F hF) x,
      F.apply_symm_apply]
  map_mul' x y := by
    have hxy : ((x * y : Matrix ι ι ℂ), (0 : Matrix κ κ ℂ)) =
        ((x, 0) : TwoBlockAlgebra ι κ) * (y, 0) := by
      rw [Prod.mk_mul_mk, mul_zero]
    show (F (x * y, 0)).1 = (F (x, 0)).1 * (F (y, 0)).1
    rw [hxy, map_mul, Prod.fst_mul]
  map_add' x y := by
    have hxy : ((x + y : Matrix ι ι ℂ), (0 : Matrix κ κ ℂ)) =
        ((x, 0) : TwoBlockAlgebra ι κ) + (y, 0) := by
      rw [Prod.mk_add_mk, add_zero]
    show (F (x + y, 0)).1 = (F (x, 0)).1 + (F (y, 0)).1
    rw [hxy, map_add, Prod.fst_add]
  map_smul' c x := by
    have hcx : ((c • x : Matrix ι ι ℂ), (0 : Matrix κ κ ℂ)) =
        c • ((x, 0) : TwoBlockAlgebra ι κ) := by
      rw [Prod.smul_mk, smul_zero]
    show (F (c • x, 0)).1 = c • (F (x, 0)).1
    rw [hcx, map_smul]
    rfl
  map_star' x := by
    have hsx : ((star x : Matrix ι ι ℂ), (0 : Matrix κ κ ℂ)) =
        star ((x, 0) : TwoBlockAlgebra ι κ) := by
      refine Prod.ext ?_ ?_
      · rfl
      · show (0 : Matrix κ κ ℂ) = star 0
        rw [star_zero]
    show (F (star x, 0)).1 = star ((F (x, 0)).1)
    rw [hsx, map_star]
    rfl

/-- The second-block restriction of a block-fixing star automorphism of
the two-block algebra. -/
def blockRestrictionSnd
    (F : TwoBlockAlgebra ι κ ≃⋆ₐ[ℂ] TwoBlockAlgebra ι κ)
    (hF : F (blockProjectionSnd ι κ) = blockProjectionSnd ι κ) :
    Matrix κ κ ℂ ≃⋆ₐ[ℂ] Matrix κ κ ℂ where
  toFun y := (F (0, y)).2
  invFun y := (F.symm (0, y)).2
  left_inv y := by
    show (F.symm (0, (F (0, y)).2)).2 = y
    rw [← apply_snd_pair F hF y, F.symm_apply_apply]
  right_inv y := by
    show (F (0, (F.symm (0, y)).2)).2 = y
    rw [← apply_snd_pair F.symm (symm_fixes_blockProjectionSnd F hF) y,
      F.apply_symm_apply]
  map_mul' x y := by
    have hxy : ((0 : Matrix ι ι ℂ), (x * y : Matrix κ κ ℂ)) =
        ((0, x) : TwoBlockAlgebra ι κ) * (0, y) := by
      rw [Prod.mk_mul_mk, mul_zero]
    show (F (0, x * y)).2 = (F (0, x)).2 * (F (0, y)).2
    rw [hxy, map_mul, Prod.snd_mul]
  map_add' x y := by
    have hxy : ((0 : Matrix ι ι ℂ), (x + y : Matrix κ κ ℂ)) =
        ((0, x) : TwoBlockAlgebra ι κ) + (0, y) := by
      rw [Prod.mk_add_mk, add_zero]
    show (F (0, x + y)).2 = (F (0, x)).2 + (F (0, y)).2
    rw [hxy, map_add, Prod.snd_add]
  map_smul' c y := by
    have hcy : ((0 : Matrix ι ι ℂ), (c • y : Matrix κ κ ℂ)) =
        c • ((0, y) : TwoBlockAlgebra ι κ) := by
      rw [Prod.smul_mk, smul_zero]
    show (F (0, c • y)).2 = c • (F (0, y)).2
    rw [hcy, map_smul]
    rfl
  map_star' y := by
    have hsy : ((0 : Matrix ι ι ℂ), (star y : Matrix κ κ ℂ)) =
        star ((0, y) : TwoBlockAlgebra ι κ) := by
      refine Prod.ext ?_ ?_
      · show (0 : Matrix ι ι ℂ) = star 0
        rw [star_zero]
      · rfl
    show (F (0, star y)).2 = star ((F (0, y)).2)
    rw [hsy, map_star]
    rfl

end BlockRestriction

/-! ## The two-block Stone converse -/

section TwoBlockStone

variable {ι : Type u} {κ : Type v} [Fintype ι] [DecidableEq ι]
  [Fintype κ] [DecidableEq κ]

/-- Two-block Stone converse from pointwise continuity alone.  Every
pointwise-continuous real-parameter star-automorphism group of the direct
sum of two full finite matrix blocks fixes each block and is blockwise
conjugation by the coherent unitary families of one time-independent
self-adjoint Hamiltonian per block.  Block fixing is derived from
continuity and connectedness of the parameter line, so no
block-preservation hypothesis is taken. -/
theorem twoBlockStoneConverse_of_continuous [Nonempty ι] [Nonempty κ]
    (α : ℝ → (TwoBlockAlgebra ι κ ≃⋆ₐ[ℂ] TwoBlockAlgebra ι κ))
    (hzero : ∀ x, α 0 x = x)
    (hgroup : ∀ s t : ℝ, ∀ x, α (s + t) x = α s (α t x))
    (hcont : ∀ x, Continuous fun t : ℝ => α t x) :
    ∃ (Ham₁ : Matrix ι ι ℂ) (Ham₂ : Matrix κ κ ℂ),
      IsSelfAdjoint Ham₁ ∧ IsSelfAdjoint Ham₂ ∧
        ∀ (t : ℝ) (x : TwoBlockAlgebra ι κ),
          α t x =
            (stonePropagator Ham₁ t * x.1 * stonePropagator Ham₁ (-t),
              stonePropagator Ham₂ t * x.2 * stonePropagator Ham₂ (-t)) := by
  have hfix : ∀ t, α t (blockProjectionFst ι κ) = blockProjectionFst ι κ :=
    continuousFlow_fixes_blockProjectionFst α (hzero _) (hcont _)
  have hfix₂ : ∀ t, α t (blockProjectionSnd ι κ) =
      blockProjectionSnd ι κ := by
    intro t
    have hsum : blockProjectionSnd ι κ = 1 - blockProjectionFst ι κ := by
      rw [← blockProjection_add ι κ]
      abel
    rw [hsum, map_sub, map_one, hfix t]
  obtain ⟨Ham₁, hsa₁, hconj₁⟩ :=
    finiteStoneConverse_of_continuous
      (fun t => blockRestrictionFst (α t) (hfix t))
      (fun x => by
        show (α 0 (x, 0)).1 = x
        rw [hzero])
      (fun s t x => by
        show (α (s + t) (x, 0)).1 = (α s ((α t (x, 0)).1, 0)).1
        conv_lhs => rw [hgroup s t (x, 0),
          apply_fst_pair (α t) (hfix t) x])
      (fun x => (hcont (x, 0)).fst)
  obtain ⟨Ham₂, hsa₂, hconj₂⟩ :=
    finiteStoneConverse_of_continuous
      (fun t => blockRestrictionSnd (α t) (hfix₂ t))
      (fun y => by
        show (α 0 (0, y)).2 = y
        rw [hzero])
      (fun s t y => by
        show (α (s + t) (0, y)).2 = (α s (0, (α t (0, y)).2)).2
        conv_lhs => rw [hgroup s t (0, y),
          apply_snd_pair (α t) (hfix₂ t) y])
      (fun y => (hcont (0, y)).snd)
  refine ⟨Ham₁, Ham₂, hsa₁, hsa₂, ?_⟩
  intro t x
  have hb1 : (α t (x.1, 0)).1 =
      stonePropagator Ham₁ t * x.1 * stonePropagator Ham₁ (-t) :=
    hconj₁ t x.1
  have hb2 : (α t (0, x.2)).2 =
      stonePropagator Ham₂ t * x.2 * stonePropagator Ham₂ (-t) :=
    hconj₂ t x.2
  have hxsplit : x = ((x.1, 0) : TwoBlockAlgebra ι κ) + (0, x.2) := by
    refine Prod.ext ?_ ?_
    · show x.1 = x.1 + 0
      rw [add_zero]
    · show x.2 = 0 + x.2
      rw [zero_add]
  conv_lhs => rw [hxsplit]
  rw [map_add, apply_fst_pair (α t) (hfix t) x.1,
    apply_snd_pair (α t) (hfix₂ t) x.2, hb1, hb2]
  refine Prod.ext ?_ ?_
  · show stonePropagator Ham₁ t * x.1 * stonePropagator Ham₁ (-t) + 0 = _
    rw [add_zero]
  · show 0 + stonePropagator Ham₂ t * x.2 * stonePropagator Ham₂ (-t) = _
    rw [zero_add]

end TwoBlockStone

#print axioms OPH.Dynamics.IsCentralStarProjection.map
#print axioms OPH.Dynamics.centralCutEquiv
#print axioms OPH.Dynamics.map_mem_centralCorner
#print axioms OPH.Dynamics.twoBlock_mem_center_iff
#print axioms OPH.Dynamics.twoBlock_central_idempotent_cases
#print axioms OPH.Dynamics.starAutomorphism_blockProjectionFst_dichotomy
#print axioms OPH.Dynamics.prodSwap_blockProjectionFst
#print axioms OPH.Dynamics.continuousFlow_fixes_blockProjectionFst
#print axioms OPH.Dynamics.blockRestrictionFst
#print axioms OPH.Dynamics.twoBlockStoneConverse_of_continuous

end OPH.Dynamics
