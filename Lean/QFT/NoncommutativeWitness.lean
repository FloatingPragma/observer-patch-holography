import QFT.ObserverAccessCut

/-!
# A noncommutative region-separating causal-net witness

This module addresses the E1 witness obligation of completion-plan issue
`#692` on the committed `FiniteCausalObserverNet` interface.  It supplies an
obstruction pair and a positive inhabitant, and the two together locate the
exact boundary that the supplied star-homomorphic restriction retractions
impose on noncommutative regional algebras.

The obstruction: for `n ≥ 2` the full matrix algebra `Mₙ(ℂ)` admits no
unital algebra homomorphism into the scalars (`no_scalar_hom_of_matrix`).
The matrix ring over `ℂ` is simple, so such a homomorphism would be
injective, and the off-diagonal matrix unit is a nonzero element whose
square is zero, while the scalars carry no nonzero element with square
zero.  The elementary special case: a unital algebra homomorphism into
the scalars maps units to nonzero scalars, and the images of an
anticommuting pair of units would be nonzero complex numbers whose
product equals its own negation
(`no_scalar_character_of_anticommuting_units`).  At the interface level,
`FiniteCausalObserverNet.no_matrix_factor_over_scalar_restriction` shows
that in every inhabitant at a regulator of nonzero dimension, a region
whose local algebra contains a unital copy of `Mₙ(ℂ)` with `n ≥ 2`
admits no region below it whose local algebra is the scalars, and
`FiniteCausalObserverNet.no_commutative_restriction_of_anticommuting_units`
records the elementary anticommuting-unit form against commutative
drops.  Every full matrix factor of dimension at least two above a
scalar overlap region is therefore impossible on this interface: the
standard tensor-factor picture with trivial overlap is excluded by the
restriction field, and the committed extra receipts are strictly
stronger than ordinary isotony.

The witness: `characterCausalNet` inhabits the interface over a
one-regulator dimension-four constant tower with scalar public algebra.  It
declares three regions, a bottom overlap region carrying the scalar
algebra and two declared-disjoint maximal regions carrying the block
algebras of matrices `diag(a, a, A)` and `diag(B, b, b)` with arbitrary
two-by-two blocks `A` and `B`.  Both regional algebras are noncommutative
(`characterCausalNet_left_noncommutative`,
`characterCausalNet_right_noncommutative`), neither is contained in the
other (`characterCausalNet_region_separating`), their elements commute
pairwise as the locality field records, and each admits a character onto
the scalars through its one-dimensional block, which supplies the
star-homomorphic restriction retractions the interface demands.  The
obstruction explains the block shape: a noncommutative regional algebra
compatible with a scalar drop must carry a one-dimensional block, and the
witness realizes exactly that minimal shape.

Claim boundary: everything here is finite linear algebra over one declared
regulator with the zero tower generator.  The regional repair maps are the
identity; the interface laws hold for it, and nontrivial channel or
instrument semantics belong to E2.  The declared overlap algebra is a
proper subalgebra of the set intersection of the two block algebras
(`overlapAlgebra_lt_regional_intersection` records the inclusion and
exhibits the gap); the
interface requires only the isotony inclusions, and the declared bottom
algebra is the scalar span characterized by `oneBlockPartition_mem_iff`.
No coverage law beyond the declared three-element region poset, no
physical causality, no CP or CPTP structure, no source realization, and no
continuum statement is claimed.  E1's genuine coverage semantics and the
regional-factor receipts for a source-produced net are open.
-/

namespace OPH.QFT

open Matrix
open OPH.Tower
open EventAlgebra

universe u

variable {ι : Type u} [Preorder ι]

/-! ## The scalar-character obstruction -/

/-- No unital complex-algebra homomorphism into the scalars accepts an
anticommuting pair of units.  The homomorphism sends units to nonzero
scalars, while the anticommutation relation forces the product of the two
image scalars to vanish. -/
theorem no_scalar_character_of_anticommuting_units {A : Type*} [Ring A]
    [Algebra ℂ A] (φ : A →ₐ[ℂ] ℂ) {u v : A} (hu : IsUnit u) (hv : IsUnit v)
    (huv : u * v = -(v * u)) : False := by
  have hmap : φ u * φ v = -(φ v * φ u) := by
    rw [← map_mul, huv, map_neg, map_mul]
  rw [mul_comm (φ v) (φ u)] at hmap
  have hzero : φ u * φ v = 0 :=
    add_self_eq_zero.mp (by linear_combination hmap)
  rcases mul_eq_zero.mp hzero with h | h
  · exact (hu.map φ).ne_zero h
  · exact (hv.map φ).ne_zero h

/-- No unital complex-algebra homomorphism carries a full matrix algebra
of dimension at least two into the scalars.  The matrix ring over `ℂ` is
simple, so the kernel of such a homomorphism is trivial and the
homomorphism is injective; the off-diagonal matrix unit is a nonzero
element whose square is zero, and its image would be a nonzero scalar
whose square is zero. -/
theorem no_scalar_hom_of_matrix (n : ℕ) (hn : 2 ≤ n)
    (φ : Matrix (Fin n) (Fin n) ℂ →ₐ[ℂ] ℂ) : False := by
  have h0 : 0 < n := by omega
  have h1 : 1 < n := by omega
  haveI : Nonempty (Fin n) := ⟨⟨0, h0⟩⟩
  have hne : (⟨1, h1⟩ : Fin n) ≠ ⟨0, h0⟩ := by
    simp [Fin.ext_iff]
  have hEE : Matrix.single (⟨0, h0⟩ : Fin n) (⟨1, h1⟩ : Fin n) (1 : ℂ) *
      Matrix.single ⟨0, h0⟩ ⟨1, h1⟩ 1 = 0 := by
    simp [hne]
  have hφE : φ (Matrix.single ⟨0, h0⟩ ⟨1, h1⟩ 1) = 0 := by
    have h := map_mul φ (Matrix.single (⟨0, h0⟩ : Fin n) ⟨1, h1⟩ (1 : ℂ))
      (Matrix.single ⟨0, h0⟩ ⟨1, h1⟩ 1)
    rw [hEE, map_zero] at h
    exact mul_self_eq_zero.mp h.symm
  have hinj : Function.Injective φ := RingHom.injective φ.toRingHom
  have hE0 : Matrix.single (⟨0, h0⟩ : Fin n) (⟨1, h1⟩ : Fin n) (1 : ℂ) = 0 :=
    hinj (by rw [hφE, map_zero])
  have h01 := congrFun (congrFun hE0 ⟨0, h0⟩) ⟨1, h1⟩
  simp [Matrix.single_apply_same] at h01

/-- Interface-level obstruction.  In any finite causal observer net, a
declared region inclusion `V ≤ U` with commutative local algebra at `V`
excludes anticommuting pairs of units from the local algebra at `U`: the
supplied restriction is a unital star homomorphism into a commutative
algebra, and the images of the two units would be units whose product both
equals and negates itself.  With `V` a scalar-algebra region this excludes
every full matrix factor at `U`, so the committed restriction receipts are
strictly stronger than ordinary isotony and rule out the tensor-factor
picture with trivial overlap. -/
theorem FiniteCausalObserverNet.no_commutative_restriction_of_anticommuting_units
    {T : ConsensusTower ι} (N : FiniteCausalObserverNet T) {r : ι}
    {U V : N.Region r} (hVU : N.regionLE r V U) (i : Fin (T.dim r))
    (hcomm : ∀ X Y : N.localAlgebra r V, X * Y = Y * X)
    {u v : N.localAlgebra r U} (hu : IsUnit u) (hv : IsUnit v)
    (huv : u * v = -(v * u)) : False := by
  set φ := N.restrict r hVU
  have hmap : φ u * φ v = -(φ v * φ u) := by
    rw [← map_mul, huv, map_neg, map_mul]
  rw [hcomm (φ v) (φ u)] at hmap
  have hval : ((φ u * φ v : N.localAlgebra r V) :
      ConsensusTower.PrivateAlgebra T r) =
      -((φ u * φ v : N.localAlgebra r V) :
        ConsensusTower.PrivateAlgebra T r) :=
    congrArg Subtype.val hmap
  have hcoe : ((φ u * φ v : N.localAlgebra r V) :
      ConsensusTower.PrivateAlgebra T r) = 0 := by
    ext a b
    have hab := congrFun (congrFun hval a) b
    rw [Matrix.neg_apply] at hab
    have h2 : ((φ u * φ v : N.localAlgebra r V) :
        ConsensusTower.PrivateAlgebra T r) a b +
        ((φ u * φ v : N.localAlgebra r V) :
          ConsensusTower.PrivateAlgebra T r) a b = 0 := by
      linear_combination hab
    simpa using add_self_eq_zero.mp h2
  have hzero : (φ u * φ v : N.localAlgebra r V) = 0 := Subtype.ext hcoe
  have hunit : IsUnit (φ u * φ v) := (hu.map φ).mul (hv.map φ)
  rw [hzero, isUnit_zero_iff] at hunit
  have h01 : ((0 : N.localAlgebra r V) :
      ConsensusTower.PrivateAlgebra T r) i i =
      ((1 : N.localAlgebra r V) :
        ConsensusTower.PrivateAlgebra T r) i i := by
    rw [hunit]
  simp at h01

/-- Interface-level matrix-factor obstruction.  In any finite causal
observer net at a regulator of nonzero dimension, a declared region
inclusion `V ≤ U` whose lower local algebra consists of scalar multiples
of the identity excludes every unital copy of a full matrix algebra of
dimension at least two from the local algebra at `U`.  The supplied
restriction retraction composes with the copy to a unital algebra
homomorphism from a simple matrix ring into a nontrivial target whose
elements are scalars: the composite is injective by simplicity, while
the off-diagonal matrix unit maps to a scalar multiple of the identity
whose square vanishes, hence to zero.  In particular the strict drop of
a full matrix factor onto a scalar overlap region is impossible, and a
unital star-embedding of `M₂` above a scalar region is excluded through
its underlying algebra homomorphism. -/
theorem FiniteCausalObserverNet.no_matrix_factor_over_scalar_restriction
    {T : ConsensusTower ι} (N : FiniteCausalObserverNet T) {r : ι}
    {U V : N.Region r} (hVU : N.regionLE r V U) (i : Fin (T.dim r))
    {n : ℕ} (hn : 2 ≤ n)
    (ψ : Matrix (Fin n) (Fin n) ℂ →ₐ[ℂ] N.localAlgebra r U)
    (hscalar : ∀ X : N.localAlgebra r V, ∃ c : ℂ,
      (X : ConsensusTower.PrivateAlgebra T r) = c • 1) : False := by
  have h0 : 0 < n := by omega
  have h1 : 1 < n := by omega
  haveI : Nonempty (Fin n) := ⟨⟨0, h0⟩⟩
  haveI : Nontrivial (N.localAlgebra r V) := by
    refine nontrivial_of_ne 0 1 fun h01 => ?_
    have hval := congrArg Subtype.val h01
    have hii := congrFun (congrFun hval i) i
    simp at hii
  set ρ : Matrix (Fin n) (Fin n) ℂ →ₐ[ℂ] N.localAlgebra r V :=
    ((N.restrict r hVU).toAlgHom.comp ψ) with hρdef
  have hne : (⟨1, h1⟩ : Fin n) ≠ ⟨0, h0⟩ := by
    simp [Fin.ext_iff]
  have hEE : Matrix.single (⟨0, h0⟩ : Fin n) (⟨1, h1⟩ : Fin n) (1 : ℂ) *
      Matrix.single ⟨0, h0⟩ ⟨1, h1⟩ 1 = 0 := by
    simp [hne]
  obtain ⟨c, hc⟩ := hscalar (ρ (Matrix.single ⟨0, h0⟩ ⟨1, h1⟩ 1))
  have hsq : ρ (Matrix.single ⟨0, h0⟩ ⟨1, h1⟩ 1) *
      ρ (Matrix.single ⟨0, h0⟩ ⟨1, h1⟩ 1) = 0 := by
    rw [← map_mul, hEE, map_zero]
  have hval : (ρ (Matrix.single ⟨0, h0⟩ ⟨1, h1⟩ 1) :
      ConsensusTower.PrivateAlgebra T r) *
      (ρ (Matrix.single ⟨0, h0⟩ ⟨1, h1⟩ 1) :
        ConsensusTower.PrivateAlgebra T r) = 0 := by
    have h := congrArg Subtype.val hsq
    simpa using h
  rw [hc, smul_mul_smul_comm, one_mul] at hval
  have hii := congrFun (congrFun hval i) i
  rw [Matrix.smul_apply, Matrix.zero_apply, Matrix.one_apply_eq,
    smul_eq_mul, mul_one] at hii
  have hc0 : c = 0 := mul_self_eq_zero.mp hii
  have hρE : ρ (Matrix.single ⟨0, h0⟩ ⟨1, h1⟩ 1) = 0 := by
    apply Subtype.ext
    rw [hc, hc0, zero_smul]
    simp
  have hinj : Function.Injective ρ := RingHom.injective ρ.toRingHom
  have hE0 : Matrix.single (⟨0, h0⟩ : Fin n) (⟨1, h1⟩ : Fin n) (1 : ℂ) = 0 :=
    hinj (by rw [hρE, map_zero])
  have h01 := congrFun (congrFun hE0 ⟨0, h0⟩) ⟨1, h1⟩
  simp [Matrix.single_apply_same] at h01

/-! ## The two block algebras inside `M₄(ℂ)`

`leftBlockAlgebra` holds the matrices `diag(a, a, A)`: the scalar `a` on
the first two diagonal entries and an arbitrary two-by-two block `A` on
the last two coordinates.  `rightBlockAlgebra` holds `diag(B, b, b)`.
Both carry a one-dimensional block, hence a character, and their elements
commute pairwise. -/

/-- The dimension-four matrix unit supported at entry `(i, j)`. -/
def blockUnit (i j : Fin 4) : Matrix (Fin 4) (Fin 4) ℂ :=
  fun k l => if k = i ∧ l = j then 1 else 0

/-- The block algebra of matrices `diag(a, a, A)`: scalar upper block,
free lower two-by-two block, vanishing mixed blocks. -/
def leftBlockAlgebra : StarSubalgebra ℂ (Matrix (Fin 4) (Fin 4) ℂ) where
  carrier := {X | X 0 1 = 0 ∧ X 1 0 = 0 ∧ X 1 1 = X 0 0 ∧
    X 0 2 = 0 ∧ X 0 3 = 0 ∧ X 1 2 = 0 ∧ X 1 3 = 0 ∧
    X 2 0 = 0 ∧ X 2 1 = 0 ∧ X 3 0 = 0 ∧ X 3 1 = 0}
  zero_mem' := by simp [Set.mem_setOf_eq]
  one_mem' := by simp [Set.mem_setOf_eq]
  add_mem' := by
    intro X Y hX hY
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ := hX
    obtain ⟨g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11⟩ := hY
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
      simp [Matrix.add_apply, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11,
        g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11]
  mul_mem' := by
    intro X Y hX hY
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ := hX
    obtain ⟨g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11⟩ := hY
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
      simp [Matrix.mul_apply, Fin.sum_univ_four, h1, h2, h3, h4, h5, h6, h7,
        h8, h9, h10, h11, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11]
  algebraMap_mem' := by
    intro c
    simp [Set.mem_setOf_eq, Algebra.algebraMap_eq_smul_one,
      Matrix.smul_apply]
  star_mem' := by
    intro X hX
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ := hX
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
      simp [h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11]

theorem mem_leftBlockAlgebra_iff {X : Matrix (Fin 4) (Fin 4) ℂ} :
    X ∈ leftBlockAlgebra ↔
      (X 0 1 = 0 ∧ X 1 0 = 0 ∧ X 1 1 = X 0 0 ∧
        X 0 2 = 0 ∧ X 0 3 = 0 ∧ X 1 2 = 0 ∧ X 1 3 = 0 ∧
        X 2 0 = 0 ∧ X 2 1 = 0 ∧ X 3 0 = 0 ∧ X 3 1 = 0) :=
  Iff.rfl

/-- The block algebra of matrices `diag(B, b, b)`: free upper two-by-two
block, scalar lower block, vanishing mixed blocks. -/
def rightBlockAlgebra : StarSubalgebra ℂ (Matrix (Fin 4) (Fin 4) ℂ) where
  carrier := {X | X 2 3 = 0 ∧ X 3 2 = 0 ∧ X 3 3 = X 2 2 ∧
    X 0 2 = 0 ∧ X 0 3 = 0 ∧ X 1 2 = 0 ∧ X 1 3 = 0 ∧
    X 2 0 = 0 ∧ X 2 1 = 0 ∧ X 3 0 = 0 ∧ X 3 1 = 0}
  zero_mem' := by simp [Set.mem_setOf_eq]
  one_mem' := by simp [Set.mem_setOf_eq]
  add_mem' := by
    intro X Y hX hY
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ := hX
    obtain ⟨g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11⟩ := hY
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
      simp [Matrix.add_apply, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11,
        g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11]
  mul_mem' := by
    intro X Y hX hY
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ := hX
    obtain ⟨g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11⟩ := hY
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
      simp [Matrix.mul_apply, Fin.sum_univ_four, h1, h2, h3, h4, h5, h6, h7,
        h8, h9, h10, h11, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11]
  algebraMap_mem' := by
    intro c
    simp [Set.mem_setOf_eq, Algebra.algebraMap_eq_smul_one,
      Matrix.smul_apply]
  star_mem' := by
    intro X hX
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ := hX
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
      simp [h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11]

theorem mem_rightBlockAlgebra_iff {X : Matrix (Fin 4) (Fin 4) ℂ} :
    X ∈ rightBlockAlgebra ↔
      (X 2 3 = 0 ∧ X 3 2 = 0 ∧ X 3 3 = X 2 2 ∧
        X 0 2 = 0 ∧ X 0 3 = 0 ∧ X 1 2 = 0 ∧ X 1 3 = 0 ∧
        X 2 0 = 0 ∧ X 2 1 = 0 ∧ X 3 0 = 0 ∧ X 3 1 = 0) :=
  Iff.rfl

/-- Elementwise commutation across the two block algebras: each algebra is
scalar on the coordinates where the other has its free block. -/
theorem leftBlock_commute_rightBlock {X Y : Matrix (Fin 4) (Fin 4) ℂ}
    (hX : X ∈ leftBlockAlgebra) (hY : Y ∈ rightBlockAlgebra) :
    Commute X Y := by
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ :=
    mem_leftBlockAlgebra_iff.mp hX
  obtain ⟨g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11⟩ :=
    mem_rightBlockAlgebra_iff.mp hY
  show X * Y = Y * X
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_four, h1, h2, h3, h4, h5, h6, h7,
      h8, h9, h10, h11, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11] <;>
    ring

/-! ## Block membership receipts for the matrix units -/

theorem blockUnit_two_three_mem_left : blockUnit 2 3 ∈ leftBlockAlgebra := by
  simp [mem_leftBlockAlgebra_iff, blockUnit]

theorem blockUnit_three_two_mem_left : blockUnit 3 2 ∈ leftBlockAlgebra := by
  simp [mem_leftBlockAlgebra_iff, blockUnit]

theorem blockUnit_zero_one_mem_right : blockUnit 0 1 ∈ rightBlockAlgebra := by
  simp [mem_rightBlockAlgebra_iff, blockUnit]

theorem blockUnit_one_zero_mem_right : blockUnit 1 0 ∈ rightBlockAlgebra := by
  simp [mem_rightBlockAlgebra_iff, blockUnit]

theorem blockUnit_two_three_not_mem_right :
    blockUnit 2 3 ∉ rightBlockAlgebra := by
  intro h
  have h23 := (mem_rightBlockAlgebra_iff.mp h).1
  simp [blockUnit] at h23

theorem blockUnit_zero_one_not_mem_left :
    blockUnit 0 1 ∉ leftBlockAlgebra := by
  intro h
  have h01 := (mem_leftBlockAlgebra_iff.mp h).1
  simp [blockUnit] at h01

/-- The left block algebra is noncommutative: the two lower-block matrix
units have distinct products. -/
theorem leftBlockAlgebra_noncommutative :
    ∃ X Y : leftBlockAlgebra,
      ¬ Commute (X : Matrix (Fin 4) (Fin 4) ℂ)
        (Y : Matrix (Fin 4) (Fin 4) ℂ) := by
  refine ⟨⟨blockUnit 2 3, blockUnit_two_three_mem_left⟩,
    ⟨blockUnit 3 2, blockUnit_three_two_mem_left⟩, fun h => ?_⟩
  have h22 := congrFun (congrFun h.eq 2) 2
  simp [blockUnit, Matrix.mul_apply] at h22

/-- The right block algebra is noncommutative: the two upper-block matrix
units have distinct products. -/
theorem rightBlockAlgebra_noncommutative :
    ∃ X Y : rightBlockAlgebra,
      ¬ Commute (X : Matrix (Fin 4) (Fin 4) ℂ)
        (Y : Matrix (Fin 4) (Fin 4) ℂ) := by
  refine ⟨⟨blockUnit 0 1, blockUnit_zero_one_mem_right⟩,
    ⟨blockUnit 1 0, blockUnit_one_zero_mem_right⟩, fun h => ?_⟩
  have h00 := congrFun (congrFun h.eq 0) 0
  simp [blockUnit, Matrix.mul_apply] at h00

/-- Region separation, one direction: the lower-block raising unit lies in
the left algebra and outside the right algebra. -/
theorem leftBlock_not_le_rightBlock :
    ¬ leftBlockAlgebra ≤ rightBlockAlgebra := by
  intro hle
  exact blockUnit_two_three_not_mem_right (hle blockUnit_two_three_mem_left)

/-- Region separation, converse direction: the upper-block raising unit
lies in the right algebra and outside the left algebra. -/
theorem rightBlock_not_le_leftBlock :
    ¬ rightBlockAlgebra ≤ leftBlockAlgebra := by
  intro hle
  exact blockUnit_zero_one_not_mem_left (hle blockUnit_zero_one_mem_right)

/-- The projector onto the two upper coordinates: a joint member of both
block algebras outside the scalar span.  The declared overlap algebra of
the witness net below is therefore a proper subalgebra of the set
intersection of the two regional algebras. -/
def upperPairProjector : Matrix (Fin 4) (Fin 4) ℂ :=
  blockUnit 0 0 + blockUnit 1 1

/-- Every scalar public element lies in both regional block algebras: the
declared overlap algebra is contained in the set intersection of the two
regional algebras. -/
theorem overlapAlgebra_le_regional_intersection :
    ∀ X ∈ (oneBlockPartition 4).publicSubalgebra,
      X ∈ leftBlockAlgebra ∧ X ∈ rightBlockAlgebra := fun _ hX =>
  ⟨oneBlockPartition_publicSubalgebra_le _ hX,
    oneBlockPartition_publicSubalgebra_le _ hX⟩

/-- The declared overlap algebra sits strictly inside the set intersection
of the two regional block algebras: the inclusion holds, and the projector
onto the two upper coordinates is a joint member outside the scalar
span. -/
theorem overlapAlgebra_lt_regional_intersection :
    (∀ X ∈ (oneBlockPartition 4).publicSubalgebra,
        X ∈ leftBlockAlgebra ∧ X ∈ rightBlockAlgebra) ∧
      ∃ X, X ∈ leftBlockAlgebra ∧ X ∈ rightBlockAlgebra ∧
        X ∉ (oneBlockPartition 4).publicSubalgebra := by
  refine ⟨overlapAlgebra_le_regional_intersection,
    upperPairProjector, ?_, ?_, ?_⟩
  · simp [mem_leftBlockAlgebra_iff, upperPairProjector, blockUnit,
      Matrix.add_apply]
  · simp [mem_rightBlockAlgebra_iff, upperPairProjector, blockUnit,
      Matrix.add_apply]
  · intro hmem
    obtain ⟨c, hc⟩ := (oneBlockPartition_mem_iff _).mp hmem
    have h00 := congrFun (congrFun hc 0) 0
    have h22 := congrFun (congrFun hc 2) 2
    simp [upperPairProjector, blockUnit, Matrix.add_apply,
      Matrix.smul_apply] at h00 h22
    rw [h00] at h22
    exact one_ne_zero h22

/-! ## The characters of the two block algebras -/

theorem smul_one_mem_oneBlock (c : ℂ) :
    c • (1 : Matrix (Fin 4) (Fin 4) ℂ) ∈
      (oneBlockPartition 4).publicSubalgebra :=
  (oneBlockPartition_mem_iff _).mpr ⟨c, rfl⟩

/-- The character of the left block algebra: evaluation of the scalar
block, returned as a scalar public element.  Multiplicativity holds
because the first row of a left-block matrix vanishes away from the
corner entry. -/
noncomputable def leftBlockCharacter :
    leftBlockAlgebra →⋆ₐ[ℂ] (oneBlockPartition 4).publicSubalgebra where
  toFun X := ⟨(X : Matrix (Fin 4) (Fin 4) ℂ) 0 0 • 1,
    smul_one_mem_oneBlock _⟩
  map_zero' := by
    apply Subtype.ext
    simp
  map_one' := by
    apply Subtype.ext
    simp
  map_add' X Y := by
    apply Subtype.ext
    simp [Matrix.add_apply, add_smul]
  map_mul' X Y := by
    apply Subtype.ext
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ :=
      mem_leftBlockAlgebra_iff.mp X.2
    simp [Matrix.mul_apply, Fin.sum_univ_four, h1, h4, h5, smul_smul,
      mul_comm]
  commutes' c := by
    apply Subtype.ext
    simp [Algebra.algebraMap_eq_smul_one, Matrix.smul_apply]
  map_star' X := by
    apply Subtype.ext
    simp [star_smul]

/-- The character of the right block algebra: evaluation of its scalar
block at the third diagonal entry. -/
noncomputable def rightBlockCharacter :
    rightBlockAlgebra →⋆ₐ[ℂ] (oneBlockPartition 4).publicSubalgebra where
  toFun X := ⟨(X : Matrix (Fin 4) (Fin 4) ℂ) 2 2 • 1,
    smul_one_mem_oneBlock _⟩
  map_zero' := by
    apply Subtype.ext
    simp
  map_one' := by
    apply Subtype.ext
    simp
  map_add' X Y := by
    apply Subtype.ext
    simp [Matrix.add_apply, add_smul]
  map_mul' X Y := by
    apply Subtype.ext
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ :=
      mem_rightBlockAlgebra_iff.mp X.2
    simp [Matrix.mul_apply, Fin.sum_univ_four, h1, h8, h9, smul_smul,
      mul_comm]
  commutes' c := by
    apply Subtype.ext
    simp [Algebra.algebraMap_eq_smul_one, Matrix.smul_apply]
  map_star' X := by
    apply Subtype.ext
    simp [star_smul]

/-! ## The three-region poset -/

/-- Region labels for the witness net: a bottom overlap region and two
maximal block regions. -/
inductive CharacterRegion : Type
  | bot
  | left
  | right
  deriving DecidableEq, Fintype

namespace CharacterRegion

/-- Region inclusion: reflexive, and the bottom region is below every
region. -/
def le (U V : CharacterRegion) : Prop := U = V ∨ U = bot

instance : DecidableRel le := fun U V =>
  decidable_of_iff (U = V ∨ U = bot) Iff.rfl

/-- Declared overlap: the common region on the diagonal and the bottom
region for distinct pairs. -/
def meet (U V : CharacterRegion) : CharacterRegion := if U = V then U else bot

/-- Declared disjointness: exactly the two maximal block regions. -/
def separated (U V : CharacterRegion) : Prop :=
  (U = left ∧ V = right) ∨ (U = right ∧ V = left)

instance : DecidableRel separated := fun U V =>
  decidable_of_iff ((U = left ∧ V = right) ∨ (U = right ∧ V = left)) Iff.rfl

theorem le_refl : ∀ U, le U U := by decide

theorem le_trans : ∀ U V W, le U V → le V W → le U W := by decide

theorem le_antisymm : ∀ U V, le U V → le V U → U = V := by decide

theorem meet_le_left : ∀ U V, le (meet U V) U := by decide

theorem meet_le_right : ∀ U V, le (meet U V) V := by decide

theorem le_meet : ∀ W U V, le W U → le W V → le W (meet U V) := by decide

theorem separated_symm : ∀ U V, separated U V → separated V U := by decide

theorem separated_irrefl : ∀ U, ¬ separated U U := by decide

end CharacterRegion

/-! ## The dimension-four witness tower -/

/-- The selected witness state: the pure state supported on the first
basis direction, certified through the diagonal partition's first
projector. -/
noncomputable def characterState : StateMatrix 4 :=
  ⟨Matrix.diagonal (Pi.single 0 1),
    ((diagonalPartition 4).isEvent 0).posSemidef, by
      rw [Matrix.trace_diagonal]
      simp [Pi.single_apply]⟩

/-- The witness tower: the constant one-regulator tower of the one-block
partition at dimension four.  The public algebra is the scalar span and
the generator is zero. -/
noncomputable def characterTower : ConsensusTower Unit :=
  ConsensusTower.constantConsensusTower (oneBlockPartition 4) characterState

/-- The regional algebra assignment: scalars at the bottom region and the
two block algebras at the maximal regions. -/
noncomputable def characterRegionAlgebra :
    CharacterRegion → StarSubalgebra ℂ (Matrix (Fin 4) (Fin 4) ℂ)
  | .bot => (oneBlockPartition 4).publicSubalgebra
  | .left => leftBlockAlgebra
  | .right => rightBlockAlgebra

/-- The supplied restriction system: identities on equal regions and the
two block characters on the drops to the bottom region.  The remaining
constructor pairs are excluded by the region order. -/
noncomputable def characterRestrict :
    ∀ U V : CharacterRegion, CharacterRegion.le V U →
      (characterRegionAlgebra U →⋆ₐ[ℂ] characterRegionAlgebra V)
  | .bot, .bot, _ => StarAlgHom.id ℂ _
  | .left, .left, _ => StarAlgHom.id ℂ _
  | .right, .right, _ => StarAlgHom.id ℂ _
  | .left, .bot, _ => leftBlockCharacter
  | .right, .bot, _ => rightBlockCharacter
  | .bot, .left, h => absurd h (by decide)
  | .bot, .right, h => absurd h (by decide)
  | .left, .right, h => absurd h (by decide)
  | .right, .left, h => absurd h (by decide)

/-- Diagonal evaluation of a scalar public element reproduces the element.
The restriction-retraction receipts of the witness net consume this
identity. -/
theorem scalar_corner_eval (X : (oneBlockPartition 4).publicSubalgebra)
    (i : Fin 4) :
    (X : Matrix (Fin 4) (Fin 4) ℂ) i i • (1 : Matrix (Fin 4) (Fin 4) ℂ) =
      (X : Matrix (Fin 4) (Fin 4) ℂ) := by
  obtain ⟨c, hc⟩ := (oneBlockPartition_mem_iff _).mp X.2
  rw [← hc]
  simp [Matrix.smul_apply]

/-! ## The witness net -/

/-- The character causal net: a finite causal observer net over the
dimension-four constant tower whose two declared-disjoint maximal regions
carry noncommutative, mutually noncontaining block algebras, with
character restrictions onto the scalar bottom region.  The repair maps
are the identity; the interface laws hold for it, and nontrivial repair
semantics belong to E2. -/
noncomputable def characterCausalNet : FiniteCausalObserverNet characterTower where
  Region := fun _ => CharacterRegion
  regionFintype := fun _ => inferInstance
  regionNonempty := fun _ => ⟨CharacterRegion.bot⟩
  regionLE := fun _ => CharacterRegion.le
  regionLE_refl := fun _ U => CharacterRegion.le_refl U
  regionLE_trans := fun _ {U V W} hUV hVW =>
    CharacterRegion.le_trans U V W hUV hVW
  regionLE_antisymm := fun _ {U V} hUV hVU =>
    CharacterRegion.le_antisymm U V hUV hVU
  overlap := fun _ => CharacterRegion.meet
  overlap_le_left := fun _ U V => CharacterRegion.meet_le_left U V
  overlap_le_right := fun _ U V => CharacterRegion.meet_le_right U V
  le_overlap := fun _ {W U V} hWU hWV =>
    CharacterRegion.le_meet W U V hWU hWV
  disjoint := fun _ => CharacterRegion.separated
  disjoint_symm := fun _ {U V} h => CharacterRegion.separated_symm U V h
  disjoint_irrefl := fun _ U => CharacterRegion.separated_irrefl U
  localAlgebra := fun _ => characterRegionAlgebra
  isotony := by
    intro r U V hUV
    rcases hUV with rfl | rfl
    · exact le_rfl
    · exact oneBlockPartition_publicSubalgebra_le _
  locality := by
    intro r U V hUV X Y
    rcases hUV with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
    · exact leftBlock_commute_rightBlock X.2 Y.2
    · exact (leftBlock_commute_rightBlock Y.2 X.2).symm
  restrict := fun _ {U V} h => characterRestrict U V h
  restrict_refl := by
    intro r U X
    cases U <;> rfl
  restrict_trans := by
    intro r U V W hVU hWV X
    cases U <;> cases V <;> cases W <;>
      first
        | rfl
        | exact absurd hVU (by decide)
        | exact absurd hWV (by decide)
  restrict_inclusion := by
    intro r U V hUV X
    cases U with
    | bot =>
        cases V with
        | bot => rfl
        | left =>
            apply Subtype.ext
            exact scalar_corner_eval X 0
        | right =>
            apply Subtype.ext
            exact scalar_corner_eval X 2
    | left =>
        cases V with
        | bot => exact absurd hUV (by decide)
        | left => rfl
        | right => exact absurd hUV (by decide)
    | right =>
        cases V with
        | bot => exact absurd hUV (by decide)
        | left => exact absurd hUV (by decide)
        | right => rfl
  regionRefine := fun _ U => U
  region_refine_refl := by intros; rfl
  region_refine_trans := by intros; rfl
  region_refine_mono := by intros; assumption
  overlap_natural := by intros; rfl
  disjoint_natural := by intros; assumption
  localAlgebra_natural := by intros; assumption
  repair := fun _ _ => LinearMap.id
  repair_idempotent := by intros; rfl
  repair_fixes_region := by intros; rfl
  repair_fixes_disjoint := by intros; rfl
  repair_natural := by intros; rfl

/-! ## Witness receipts -/

/-- The two maximal regions are declared disjoint. -/
theorem characterCausalNet_disjoint_left_right :
    characterCausalNet.disjoint () CharacterRegion.left
      CharacterRegion.right :=
  Or.inl ⟨rfl, rfl⟩

/-- Elementwise commutation across the declared-disjoint pair,
instantiated from the structure's locality field. -/
theorem characterCausalNet_locality_receipt
    (X : characterCausalNet.localAlgebra () CharacterRegion.left)
    (Y : characterCausalNet.localAlgebra () CharacterRegion.right) :
    Commute (X : ConsensusTower.PrivateAlgebra characterTower ())
      (Y : ConsensusTower.PrivateAlgebra characterTower ()) :=
  characterCausalNet.locality () characterCausalNet_disjoint_left_right X Y

/-- E1 receipt: the left regional algebra of the witness net is
noncommutative. -/
theorem characterCausalNet_left_noncommutative :
    ∃ X Y : characterCausalNet.localAlgebra () CharacterRegion.left,
      ¬ Commute (X : ConsensusTower.PrivateAlgebra characterTower ())
        (Y : ConsensusTower.PrivateAlgebra characterTower ()) :=
  leftBlockAlgebra_noncommutative

/-- E1 receipt: the right regional algebra of the witness net is
noncommutative. -/
theorem characterCausalNet_right_noncommutative :
    ∃ X Y : characterCausalNet.localAlgebra () CharacterRegion.right,
      ¬ Commute (X : ConsensusTower.PrivateAlgebra characterTower ())
        (Y : ConsensusTower.PrivateAlgebra characterTower ()) :=
  rightBlockAlgebra_noncommutative

/-- E1 receipt: the two disjoint regional algebras separate the regions;
neither contains the other. -/
theorem characterCausalNet_region_separating :
    ¬ characterCausalNet.localAlgebra () CharacterRegion.left ≤
        characterCausalNet.localAlgebra () CharacterRegion.right ∧
      ¬ characterCausalNet.localAlgebra () CharacterRegion.right ≤
        characterCausalNet.localAlgebra () CharacterRegion.left :=
  ⟨leftBlock_not_le_rightBlock, rightBlock_not_le_leftBlock⟩

/-- The declared overlap of the two maximal regions is the bottom
region. -/
theorem characterCausalNet_overlap_left_right :
    characterCausalNet.overlap () CharacterRegion.left
        CharacterRegion.right = CharacterRegion.bot := rfl

/-- E1 receipt: the algebra of the declared overlap of the two disjoint
block regions is exactly the scalar public span of the tower. -/
theorem characterCausalNet_meet_scalar :
    characterCausalNet.localAlgebra ()
        (characterCausalNet.overlap () CharacterRegion.left
          CharacterRegion.right) =
      (oneBlockPartition 4).publicSubalgebra := rfl

/-- The obstruction instantiated on the witness: the left region sits
above the scalar bottom region, so its algebra contains no anticommuting
pair of units.  Its noncommutativity receipts are carried by non-unit
elements, exactly as the character obstruction requires. -/
theorem characterCausalNet_left_no_anticommuting_units
    (u v : characterCausalNet.localAlgebra () CharacterRegion.left)
    (hu : IsUnit u) (hv : IsUnit v) : u * v ≠ -(v * u) := by
  intro huv
  exact characterCausalNet.no_commutative_restriction_of_anticommuting_units
    (Or.inr rfl) (0 : Fin 4)
    (fun X Y => (oneBlockPartition 4).publicSubalgebra_mul_comm X Y)
    hu hv huv

/-- The obstruction instantiated on the witness, right region. -/
theorem characterCausalNet_right_no_anticommuting_units
    (u v : characterCausalNet.localAlgebra () CharacterRegion.right)
    (hu : IsUnit u) (hv : IsUnit v) : u * v ≠ -(v * u) := by
  intro huv
  exact characterCausalNet.no_commutative_restriction_of_anticommuting_units
    (Or.inr rfl) (0 : Fin 4)
    (fun X Y => (oneBlockPartition 4).publicSubalgebra_mul_comm X Y)
    hu hv huv

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.no_scalar_character_of_anticommuting_units
#print axioms OPH.QFT.no_scalar_hom_of_matrix
#print axioms OPH.QFT.FiniteCausalObserverNet.no_commutative_restriction_of_anticommuting_units
#print axioms OPH.QFT.FiniteCausalObserverNet.no_matrix_factor_over_scalar_restriction
#print axioms OPH.QFT.mem_leftBlockAlgebra_iff
#print axioms OPH.QFT.mem_rightBlockAlgebra_iff
#print axioms OPH.QFT.blockUnit_two_three_mem_left
#print axioms OPH.QFT.blockUnit_three_two_mem_left
#print axioms OPH.QFT.blockUnit_zero_one_mem_right
#print axioms OPH.QFT.blockUnit_one_zero_mem_right
#print axioms OPH.QFT.smul_one_mem_oneBlock
#print axioms OPH.QFT.scalar_corner_eval
#print axioms OPH.QFT.characterCausalNet_overlap_left_right
#print axioms OPH.QFT.leftBlock_commute_rightBlock
#print axioms OPH.QFT.leftBlockAlgebra_noncommutative
#print axioms OPH.QFT.rightBlockAlgebra_noncommutative
#print axioms OPH.QFT.leftBlock_not_le_rightBlock
#print axioms OPH.QFT.rightBlock_not_le_leftBlock
#print axioms OPH.QFT.overlapAlgebra_le_regional_intersection
#print axioms OPH.QFT.overlapAlgebra_lt_regional_intersection
#print axioms OPH.QFT.characterCausalNet_disjoint_left_right
#print axioms OPH.QFT.characterCausalNet_locality_receipt
#print axioms OPH.QFT.characterCausalNet_left_noncommutative
#print axioms OPH.QFT.characterCausalNet_right_noncommutative
#print axioms OPH.QFT.characterCausalNet_region_separating
#print axioms OPH.QFT.characterCausalNet_meet_scalar
#print axioms OPH.QFT.characterCausalNet_left_no_anticommuting_units
#print axioms OPH.QFT.characterCausalNet_right_no_anticommuting_units
