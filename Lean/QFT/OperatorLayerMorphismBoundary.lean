import QFT.JoinNetMorphism

/-!
# The operator-layer morphism boundary over the committed carriers

`QFT.JoinNetMorphism` constructs the record-layer net morphisms of the
hinge-fibered carrier join and records one residual: the operator-layer
net morphism stays open, delimited by two negatives (the gluing
submodule is not multiplicatively closed, and the committed left
regional algebra of the 86/247 net is noncommutative while the record
layer is commutative).  This module makes that residual exact.  On the
committed 86/247 objects (the regional algebras `twoSlotAlgebra247`,
the net `twoSlotNet247`, the record layers `recordLayer247`, and the
record expectations `recordExpect247`), the boundary between the record
layer and the operator layer is an asymmetry theorem with a proved
positive complement.

**Part 1: the embedding direction exists.**  The record function
carrier embeds into the operator layer as a unital star algebra
homomorphism, the diagonal embedding `diagonalStarHom`.  It is
injective, it lands region by region exactly in the committed regional
algebras (`diagonalStarHom_mem_regionMap_iff`: the image of a record
function lies in the committed regional algebra of a region precisely
when the function lies in the record layer of that region), and it
intertwines the committed conditional expectations with the committed
record-layer expectations on every region
(`diagonalStarHom_expect_intertwine`).

**Part 2: the projection direction caps at the conditional-expectation
grade.**  The committed conditional expectation onto the left slot is
unital (`leftSlotExpectation_one`), positive
(`leftSlotExpectation_posSemidef`, committed), trace-preserving
(`leftSlotExpectation_trace`, committed), and restricts to the
committed record-layer expectation on diagonals
(`twoSlotExpect247_diagonal`, committed).  It is not multiplicative:
the committed right-slot checkpoint-style projector
`rightProjWitness = slotRight (single 0 0 1)` is idempotent, and the
expectation of its square differs from the square of its expectation
(`leftSlotExpectation_not_multiplicative`).  Positive unital
trace-preserving expectation is therefore the exact attainable grade of
operator-to-record maps; the homomorphism grade is not attainable.

**Part 3: the projection direction admits no homomorphism, in two
grades.**  Grade (a), absolute: no unital algebra homomorphism exists
from the committed left regional algebra into the record function
carrier at all (`twoSlot247_left_no_record_algHom`), and none exists
from the ambient joint algebra either
(`twoSlot247_ambient_no_record_algHom`).  Hypotheses: the map is a
unital `ℂ`-algebra homomorphism; nothing else.  Each coordinate
evaluation of such a map is a scalar character of an algebra containing
a full 13-by-13 matrix factor, which the committed
`no_scalar_restriction_of_matrix_factor` rules out.  Grade (b),
compatibility-relative: even a bare function, with no linearity,
additivity, or unitality assumed, cannot be multiplicative on the
committed left regional algebra and restrict to the committed
record-layer morphism on diagonals
(`no_multiplicative_record_projection247`).  Hypotheses, named exactly:
(H1) multiplicativity on pairs from the committed left regional algebra
`twoSlotAlgebra247 TwoSlotRegion.left`; (H2) the diagonal compatibility
clause `φ (Matrix.diagonal v) = v` for every `v` in the committed left
record layer `recordLayer247 TwoSlotRegion.left`.  The mechanism is the
collapse lemma `multiplicative_collapse_on_left`: the committed
noncommuting witness pair of `twoSlotAlgebra247_left_noncommutative`
has both of its products diagonal (the two committed record projectors
`Pi.single 0 1` and `Pi.single 1 1` pulled back along the first
coordinate), any map multiplicative on the left algebra into any
commutative multiplicative structure identifies the two products, and
the record-layer morphism separates them.  Clause (H2) is load-bearing
and not vacuous: dropping it admits trivial multiplicative maps (the
constant-one map), and the positive complement of Part 2 shows the
expectation grade satisfies the record restriction, so the no-go
isolates multiplicativity as the exact failing property.

**The composed receipt** (`operatorLayerMorphismBoundary_receipt`)
states the embedding existence with exact regional membership and
expectation intertwining, the positive complement (positivity, trace
preservation, unitality, record restriction, and the multiplicativity
failure at the committed right-slot witness), the
compatibility-relative no-go, and the absolute no-go, as one
conjunction typed against the committed premise instance
`supportDisjointPairPremises`, whose fields are the registered premise
rows PR-32 (Cartesian joint carrier), PR-33 (diamond region map), and
PR-34 (slot split).

**What stays open for the structural field-theory line.**  This module
closes the question "which morphisms exist beyond the record layer"
for maps between the committed operator layer and the committed record
layer on one carrier: embedding at the star-homomorphism grade,
projection at the conditional-expectation grade, and nothing stronger.
It does not construct morphisms between the two committed regional
nets themselves at the operator layer, a functorial net-to-net or
tower-to-tower morphism over the join, or any spacetime-covariant
structure; those stay open under the physical-attachment premise
PR-52 and the source-selected time-slice premise PR-58.

**What is NOT claimed.**  No spacetime covariance, no continuum limit,
no physical attachment, no causal reading, and no uniqueness of the
diagonal embedding is claimed.  The 86/247 carrier and its committed
regional algebras are the exact objects of every statement; no
synthetic stand-in is introduced.
-/

set_option autoImplicit false

namespace OPH.QFT

open Matrix
open Kronecker
open scoped ComplexOrder

noncomputable section

/-! ## Part 1: the diagonal star embedding of the record carrier -/

/-- The diagonal embedding of the record function carrier into the
operator layer, as a unital star algebra homomorphism. -/
def diagonalStarHom (γ : Type*) [Fintype γ] [DecidableEq γ] :
    (γ → ℂ) →⋆ₐ[ℂ] Matrix γ γ ℂ :=
  { Matrix.diagonalAlgHom ℂ with
    map_star' := fun v => by
      show Matrix.diagonal (star v) = star (Matrix.diagonal v)
      rw [Matrix.star_eq_conjTranspose, Matrix.diagonal_conjTranspose] }

theorem diagonalStarHom_apply (γ : Type*) [Fintype γ] [DecidableEq γ]
    (v : γ → ℂ) : diagonalStarHom γ v = Matrix.diagonal v := rfl

/-- The diagonal embedding is injective. -/
theorem diagonalStarHom_injective (γ : Type*) [Fintype γ] [DecidableEq γ] :
    Function.Injective (diagonalStarHom γ) :=
  fun _ _ hab => Matrix.diagonal_injective hab

/-- **Exact regional membership of the embedding.**  The image of a
record function lies in the committed regional algebra of a region
precisely when the function lies in the record layer of that region. -/
theorem diagonalStarHom_mem_regionMap_iff (U : TwoSlotRegion)
    (v : CarrierFn247) :
    diagonalStarHom PairIndex247 v ∈ twoSlotAlgebra247 U ↔
      v ∈ recordLayer247 U :=
  diagonal_mem_twoSlotAlgebra247_iff U v

/-- **The embedding intertwines the committed conditional expectations
with the committed record-layer expectations on every region.** -/
theorem diagonalStarHom_expect_intertwine (U : TwoSlotRegion)
    (v : CarrierFn247) :
    twoSlotExpect247 U (diagonalStarHom PairIndex247 v) =
      diagonalStarHom PairIndex247 (recordExpect247 U v) :=
  twoSlotExpect247_diagonal U v

/-! ## Part 2: the positive complement at the projection direction -/

variable {α β : Type*} [Fintype α] [Fintype β] [DecidableEq α] [DecidableEq β]

omit [Fintype α] in
/-- The committed left-slot conditional expectation is unital. -/
theorem leftSlotExpectation_one [Nonempty β] :
    leftSlotExpectation (α := α) (β := β) 1 = 1 := by
  have hpt : ptraceSnd (1 : Matrix (α × β) (α × β) ℂ) =
      (Fintype.card β : ℂ) • (1 : Matrix α α ℂ) := by
    have h := ptraceSnd_slotLeft (β := β) (1 : Matrix α α ℂ)
    rwa [Matrix.one_kronecker_one] at h
  have hdef : leftSlotExpectation (α := α) (β := β) 1 =
      (Fintype.card β : ℂ)⁻¹ •
        (ptraceSnd (1 : Matrix (α × β) (α × β) ℂ) ⊗ₖ
          (1 : Matrix β β ℂ)) := rfl
  rw [hdef, hpt, Matrix.smul_kronecker, Matrix.one_kronecker_one, smul_smul,
    inv_mul_cancel₀ (Nat.cast_ne_zero.mpr Fintype.card_ne_zero), one_smul]

/-- The committed right-slot projector witness: observer 247's first
checkpoint-style diagonal projector, lifted into the right slot of the
86/247 carrier through the committed slot split. -/
def rightProjWitness : Matrix PairIndex247 PairIndex247 ℂ :=
  slotRight (Fin 13) (Matrix.single (0 : Fin 13) 0 (1 : ℂ))

theorem rightProjWitness_eq :
    rightProjWitness =
      (1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
        Matrix.single (0 : Fin 13) 0 (1 : ℂ) := rfl

/-- The witness lies in the committed right regional algebra. -/
theorem rightProjWitness_mem :
    rightProjWitness ∈ twoSlotAlgebra247 TwoSlotRegion.right :=
  ⟨Matrix.single (0 : Fin 13) 0 (1 : ℂ), rfl⟩

/-- The witness is idempotent. -/
theorem rightProjWitness_idem :
    rightProjWitness * rightProjWitness = rightProjWitness := by
  unfold rightProjWitness
  rw [← map_mul, Matrix.single_mul_single_same, one_mul]

/-- The committed left expectation scalarises the witness. -/
theorem leftSlotExpectation_rightProjWitness :
    leftSlotExpectation (α := Fin 13) (β := Fin 13) rightProjWitness =
      ((13 : ℂ))⁻¹ • (1 : Matrix PairIndex247 PairIndex247 ℂ) := by
  have hdef : leftSlotExpectation (α := Fin 13) (β := Fin 13)
      rightProjWitness =
      (Fintype.card (Fin 13) : ℂ)⁻¹ •
        (ptraceSnd rightProjWitness ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ)) :=
    rfl
  rw [hdef, rightProjWitness_eq, ptraceSnd_slotRight,
    Matrix.trace_single_eq_same, one_smul, Matrix.one_kronecker_one,
    Fintype.card_fin]
  norm_num

/-- **The committed left expectation is not multiplicative.**  On the
committed right-slot projector witness, the expectation of the square
differs from the square of the expectation: the positive unital
trace-preserving grade of the projection direction does not extend to
the homomorphism grade. -/
theorem leftSlotExpectation_not_multiplicative :
    leftSlotExpectation (α := Fin 13) (β := Fin 13)
        (rightProjWitness * rightProjWitness) ≠
      leftSlotExpectation (α := Fin 13) (β := Fin 13) rightProjWitness *
        leftSlotExpectation (α := Fin 13) (β := Fin 13) rightProjWitness := by
  rw [rightProjWitness_idem, leftSlotExpectation_rightProjWitness]
  intro h
  have h2 : ((13 : ℂ)⁻¹ • (1 : Matrix PairIndex247 PairIndex247 ℂ)) *
      ((13 : ℂ)⁻¹ • (1 : Matrix PairIndex247 PairIndex247 ℂ)) =
      ((13 : ℂ)⁻¹ * (13 : ℂ)⁻¹) •
        (1 : Matrix PairIndex247 PairIndex247 ℂ) := by
    rw [Matrix.smul_mul, Matrix.mul_smul, one_mul, smul_smul]
  rw [h2] at h
  have h00 := congrFun (congrFun h ((0 : Fin 13), (0 : Fin 13)))
    ((0 : Fin 13), (0 : Fin 13))
  rw [Matrix.smul_apply, Matrix.smul_apply, Matrix.one_apply_eq,
    smul_eq_mul, smul_eq_mul, mul_one, mul_one] at h00
  have h13 : (13 : ℂ) ≠ 0 := by norm_num
  have h169 : (13 : ℂ)⁻¹ ≠ 0 := inv_ne_zero h13
  have hcancel := mul_left_cancel₀ h169
    (h00.symm.trans (mul_one (13 : ℂ)⁻¹).symm)
  rw [inv_eq_one_div] at hcancel
  field_simp at hcancel
  norm_num at hcancel

/-! ## Part 3, grade (b): the multiplicative collapse and the
compatibility-relative no-go -/

/-- The committed record projector at a left label: the one-hot record
function of the label, pulled back along the first coordinate.  It lies
in the committed left record layer. -/
def leftRecordProjector (a : Fin 13) : CarrierFn247 :=
  fun p => (Pi.single a (1 : ℂ) : Fin 13 → ℂ) p.1

theorem leftRecordProjector_mem (a : Fin 13) :
    leftRecordProjector a ∈ recordLayer247 TwoSlotRegion.left :=
  ⟨(Pi.single a (1 : ℂ) : Fin 13 → ℂ), rfl⟩

/-- The committed noncommuting witness pair of the left regional
algebra has both of its products diagonal: a left-slot diagonal
projector pulled back along the first coordinate. -/
theorem slotLeft_single_diag247 (a : Fin 13) :
    slotLeft (Fin 13) (α := Fin 13) (Matrix.single a a (1 : ℂ)) =
      Matrix.diagonal (leftRecordProjector a) := by
  show Matrix.single a a (1 : ℂ) ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ) = _
  rw [← Matrix.diagonal_single, diagonal_kronecker_one]
  rfl

/-- **The multiplicative collapse on the committed left algebra.**  Any
map from the operator layer into any commutative multiplicative
structure that is multiplicative on pairs from the committed left
regional algebra identifies the two committed record projectors: the
noncommuting witness pair has products
`diagonal (Pi.single 0 1 ∘ fst)` and `diagonal (Pi.single 1 1 ∘ fst)`,
and the map cannot tell them apart. -/
theorem multiplicative_collapse_on_left {R : Type*} [CommMonoid R]
    (φ : Matrix PairIndex247 PairIndex247 ℂ → R)
    (hmul : ∀ X Y : Matrix PairIndex247 PairIndex247 ℂ,
      X ∈ twoSlotAlgebra247 TwoSlotRegion.left →
      Y ∈ twoSlotAlgebra247 TwoSlotRegion.left →
      φ (X * Y) = φ X * φ Y) :
    φ (Matrix.diagonal (leftRecordProjector 0)) =
      φ (Matrix.diagonal (leftRecordProjector 1)) := by
  have hX : slotLeft (Fin 13) (α := Fin 13)
      (Matrix.single (0 : Fin 13) 1 (1 : ℂ)) ∈
        twoSlotAlgebra247 TwoSlotRegion.left :=
    ⟨Matrix.single (0 : Fin 13) 1 (1 : ℂ), rfl⟩
  have hY : slotLeft (Fin 13) (α := Fin 13)
      (Matrix.single (1 : Fin 13) 0 (1 : ℂ)) ∈
        twoSlotAlgebra247 TwoSlotRegion.left :=
    ⟨Matrix.single (1 : Fin 13) 0 (1 : ℂ), rfl⟩
  have hXY : slotLeft (Fin 13) (α := Fin 13)
      (Matrix.single (0 : Fin 13) 1 (1 : ℂ)) *
      slotLeft (Fin 13) (α := Fin 13)
        (Matrix.single (1 : Fin 13) 0 (1 : ℂ)) =
      Matrix.diagonal (leftRecordProjector 0) := by
    rw [← map_mul, Matrix.single_mul_single_same, one_mul,
      slotLeft_single_diag247]
  have hYX : slotLeft (Fin 13) (α := Fin 13)
      (Matrix.single (1 : Fin 13) 0 (1 : ℂ)) *
      slotLeft (Fin 13) (α := Fin 13)
        (Matrix.single (0 : Fin 13) 1 (1 : ℂ)) =
      Matrix.diagonal (leftRecordProjector 1) := by
    rw [← map_mul, Matrix.single_mul_single_same, one_mul,
      slotLeft_single_diag247]
  calc φ (Matrix.diagonal (leftRecordProjector 0))
      = φ (slotLeft (Fin 13) (α := Fin 13)
            (Matrix.single (0 : Fin 13) 1 (1 : ℂ)) *
          slotLeft (Fin 13) (α := Fin 13)
            (Matrix.single (1 : Fin 13) 0 (1 : ℂ))) := by rw [hXY]
    _ = φ (slotLeft (Fin 13) (α := Fin 13)
            (Matrix.single (0 : Fin 13) 1 (1 : ℂ))) *
        φ (slotLeft (Fin 13) (α := Fin 13)
            (Matrix.single (1 : Fin 13) 0 (1 : ℂ))) := hmul _ _ hX hY
    _ = φ (slotLeft (Fin 13) (α := Fin 13)
            (Matrix.single (1 : Fin 13) 0 (1 : ℂ))) *
        φ (slotLeft (Fin 13) (α := Fin 13)
            (Matrix.single (0 : Fin 13) 1 (1 : ℂ))) := mul_comm _ _
    _ = φ (slotLeft (Fin 13) (α := Fin 13)
            (Matrix.single (1 : Fin 13) 0 (1 : ℂ)) *
          slotLeft (Fin 13) (α := Fin 13)
            (Matrix.single (0 : Fin 13) 1 (1 : ℂ))) :=
        (hmul _ _ hY hX).symm
    _ = φ (Matrix.diagonal (leftRecordProjector 1)) := by
        rw [hYX]

/-- **The compatibility-relative no-go.**  No map from the committed
operator layer into the record function carrier is multiplicative on
the committed left regional algebra (H1) and restricts to the committed
record-layer morphism on the committed left record layer (H2).  No
linearity, additivity, unitality, or star compatibility is assumed. -/
theorem no_multiplicative_record_projection247 :
    ¬ ∃ φ : Matrix PairIndex247 PairIndex247 ℂ → CarrierFn247,
      (∀ X Y : Matrix PairIndex247 PairIndex247 ℂ,
        X ∈ twoSlotAlgebra247 TwoSlotRegion.left →
        Y ∈ twoSlotAlgebra247 TwoSlotRegion.left →
        φ (X * Y) = φ X * φ Y) ∧
      (∀ v ∈ recordLayer247 TwoSlotRegion.left,
        φ (Matrix.diagonal v) = v) := by
  rintro ⟨φ, hmul, hdiag⟩
  have hcol := multiplicative_collapse_on_left φ hmul
  rw [hdiag _ (leftRecordProjector_mem 0),
    hdiag _ (leftRecordProjector_mem 1)] at hcol
  have hval : (Pi.single (0 : Fin 13) (1 : ℂ) : Fin 13 → ℂ) 0 =
      (Pi.single (1 : Fin 13) (1 : ℂ) : Fin 13 → ℂ) 0 :=
    congrFun hcol ((0 : Fin 13), (0 : Fin 13))
  rw [Pi.single_eq_same,
    Pi.single_eq_of_ne (by decide : (0 : Fin 13) ≠ 1)] at hval
  exact one_ne_zero hval

/-! ## Part 3, grade (a): the absolute no-go for unital algebra
homomorphisms -/

/-- The committed slot embedding into the left regional algebra of the
86/247 net, as an algebra homomorphism. -/
def slotLeftIntoRange247 :
    Matrix (Fin 13) (Fin 13) ℂ →ₐ[ℂ]
      (twoSlotAlgebra247 TwoSlotRegion.left) where
  toFun A := ⟨slotLeft (Fin 13) A, ⟨A, rfl⟩⟩
  map_one' := Subtype.ext (map_one (slotLeft (Fin 13)))
  map_mul' A B := Subtype.ext (map_mul (slotLeft (Fin 13)) A B)
  map_zero' := Subtype.ext (map_zero (slotLeft (Fin 13)))
  map_add' A B := Subtype.ext (map_add (slotLeft (Fin 13)) A B)
  commutes' c := Subtype.ext (AlgHomClass.commutes (slotLeft (Fin 13)) c)

/-- **The absolute no-go, left regional algebra.**  No unital algebra
homomorphism carries the committed left regional algebra of the 86/247
net into the record function carrier, with or without any diagonal
compatibility clause: each coordinate evaluation would be a scalar
character of an algebra containing a full 13-by-13 matrix factor. -/
theorem twoSlot247_left_no_record_algHom
    (φ : (twoSlotAlgebra247 TwoSlotRegion.left) →ₐ[ℂ] CarrierFn247) :
    False :=
  no_scalar_restriction_of_matrix_factor (by norm_num)
    (twoSlotAlgebra247 TwoSlotRegion.left) slotLeftIntoRange247
    ((Pi.evalAlgHom ℂ (fun _ : PairIndex247 => ℂ)
      ((0 : Fin 13), (0 : Fin 13))).comp φ)

/-- **The absolute no-go, ambient joint algebra.**  No unital algebra
homomorphism carries the full committed joint algebra of the 86/247
carrier into the record function carrier. -/
theorem twoSlot247_ambient_no_record_algHom
    (φ : Matrix PairIndex247 PairIndex247 ℂ →ₐ[ℂ] CarrierFn247) :
    False :=
  no_scalar_restriction_of_matrix_factor (by norm_num)
    (Matrix PairIndex247 PairIndex247 ℂ)
    (slotLeft (Fin 13) (α := Fin 13)).toAlgHom
    ((Pi.evalAlgHom ℂ (fun _ : PairIndex247 => ℂ)
      ((0 : Fin 13), (0 : Fin 13))).comp φ)

/-! ## The composed receipt -/

/-- **The operator-layer morphism boundary receipt for the committed
86/247 carrier.**  The antecedent bundle is the committed premise
instance `supportDisjointPairPremises`, whose fields carry the
registered premise rows PR-32 (Cartesian joint carrier), PR-33 (diamond
region map), and PR-34 (slot split); every regional clause is typed
through its region map.  Clauses: the record carrier embeds into the
operator layer as an injective unital star algebra homomorphism whose
image lies in the committed regional algebra of a region exactly when
the record function lies in the record layer of that region, and which
intertwines the committed conditional expectations with the record
expectations; the committed left expectation is positive,
trace-preserving, and unital, and restricts to the record-layer
expectation on diagonals; the committed right-slot projector witness
breaks its multiplicativity; no map multiplicative on the committed
left regional algebra restricts to the record-layer morphism on
diagonals; and no unital algebra homomorphism carries the committed
left regional algebra into the record carrier at all.  The boundary is
the asymmetry: embedding at the star-homomorphism grade, projection at
the positive unital trace-preserving expectation grade, and nothing
stronger.  Net-to-net and tower morphism structure at the operator
layer, spacetime covariance, and physical attachment stay open under
PR-52 and PR-58. -/
theorem operatorLayerMorphismBoundary_receipt :
    (∃ ι : CarrierFn247 →⋆ₐ[ℂ]
        Matrix supportDisjointPairPremises.Carrier
          supportDisjointPairPremises.Carrier ℂ,
      Function.Injective ι ∧
      (∀ (U : TwoSlotRegion) (v : CarrierFn247),
        ι v ∈ supportDisjointPairPremises.regionMap U ↔
          v ∈ recordLayer247 U) ∧
      (∀ (U : TwoSlotRegion) (v : CarrierFn247),
        twoSlotNet247.expect U (ι v) = ι (recordExpect247 U v))) ∧
    (∀ M : Matrix PairIndex247 PairIndex247 ℂ, M.PosSemidef →
      (twoSlotNet247.expect TwoSlotRegion.left M).PosSemidef) ∧
    (∀ M : Matrix PairIndex247 PairIndex247 ℂ,
      (twoSlotNet247.expect TwoSlotRegion.left M).trace = M.trace) ∧
    (twoSlotNet247.expect TwoSlotRegion.left 1 = 1) ∧
    (∀ (U : TwoSlotRegion) (v : CarrierFn247),
      twoSlotNet247.expect U (Matrix.diagonal v) =
        Matrix.diagonal (recordExpect247 U v)) ∧
    (∃ M ∈ supportDisjointPairPremises.regionMap TwoSlotRegion.right,
      twoSlotNet247.expect TwoSlotRegion.left (M * M) ≠
        twoSlotNet247.expect TwoSlotRegion.left M *
          twoSlotNet247.expect TwoSlotRegion.left M) ∧
    (¬ ∃ φ : Matrix PairIndex247 PairIndex247 ℂ → CarrierFn247,
      (∀ X Y : Matrix PairIndex247 PairIndex247 ℂ,
        X ∈ supportDisjointPairPremises.regionMap TwoSlotRegion.left →
        Y ∈ supportDisjointPairPremises.regionMap TwoSlotRegion.left →
        φ (X * Y) = φ X * φ Y) ∧
      (∀ v ∈ recordLayer247 TwoSlotRegion.left,
        φ (Matrix.diagonal v) = v)) ∧
    (∀ _φ : (supportDisjointPairPremises.regionMap TwoSlotRegion.left)
        →ₐ[ℂ] CarrierFn247, False) :=
  ⟨⟨diagonalStarHom PairIndex247, diagonalStarHom_injective PairIndex247,
      diagonalStarHom_mem_regionMap_iff, diagonalStarHom_expect_intertwine⟩,
    fun _ hM => leftSlotExpectation_posSemidef hM,
    leftSlotExpectation_trace,
    leftSlotExpectation_one,
    twoSlotExpect247_diagonal,
    ⟨rightProjWitness, rightProjWitness_mem,
      leftSlotExpectation_not_multiplicative⟩,
    no_multiplicative_record_projection247,
    twoSlot247_left_no_record_algHom⟩

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
-- Expected axioms per declaration: propext, Classical.choice, Quot.sound.
#print axioms OPH.QFT.diagonalStarHom
#print axioms OPH.QFT.diagonalStarHom_injective
#print axioms OPH.QFT.diagonalStarHom_mem_regionMap_iff
#print axioms OPH.QFT.diagonalStarHom_expect_intertwine
#print axioms OPH.QFT.leftSlotExpectation_one
#print axioms OPH.QFT.rightProjWitness_mem
#print axioms OPH.QFT.rightProjWitness_idem
#print axioms OPH.QFT.leftSlotExpectation_rightProjWitness
#print axioms OPH.QFT.leftSlotExpectation_not_multiplicative
#print axioms OPH.QFT.slotLeft_single_diag247
#print axioms OPH.QFT.multiplicative_collapse_on_left
#print axioms OPH.QFT.no_multiplicative_record_projection247
#print axioms OPH.QFT.slotLeftIntoRange247
#print axioms OPH.QFT.twoSlot247_left_no_record_algHom
#print axioms OPH.QFT.twoSlot247_ambient_no_record_algHom
#print axioms OPH.QFT.operatorLayerMorphismBoundary_receipt
