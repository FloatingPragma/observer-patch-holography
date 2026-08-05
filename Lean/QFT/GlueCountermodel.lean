import QFT.ObserverAccessCut

/-!
# The accessible-glue countermodel

This module discharges the deferral recorded in `QFT/ObserverAccessCut.lean`:
whether ambient unique descent could force the accessible-glue receipt stayed
open there.  The answer is negative.  This file exhibits one finite causal
observer net, one access cut, and one declared two-member subregion family
with ambient unique descent on which a compatible family of accessible
sections glues to an inaccessible element.  `AccessibleGlueClosure` is
therefore logically independent of the ambient `ObserverNetDescent` packet:
the packet holds for the declared families while the closure receipt fails.

The obstruction that shapes the model: whenever the two cover restrictions
and the shared overlap restriction satisfy the linear inclusion-exclusion
identity `glue = s₁ + s₂ - s_overlap`, the glue is an algebraic combination
of the sections and one declared restriction of a section, so every
restriction-stable accessible algebra captures it and no countermodel
exists.  Symmetric two-block diagonal nets all satisfy that identity.  The
model below breaks it with an asymmetric pick: on `M₄(ℂ)` the region `u1`
carries the two-block diagonal algebra `diag(a, a, b, b)` restricted from
the top by the entry pair `(0, 2)`, the region `u2` carries the three-block
diagonal algebra `diag(a, b, c, c)` restricted by the entry triple
`(0, 1, 3)`, and the overlap is the scalar bottom region evaluated at entry
`0`.  The two restrictions read all four diagonal entries, so descent is
unique; they share only entry `0`, and the single overlap equation matches
the single shared entry, so descent is total.  The compatible accessible
family `s₁ = diag(0, 0, 1, 1)`, `s₂ = diag(0, 1, 0, 0)` has the unique glue
`diag(0, 1, 1, 0)`, and that glue violates the defining equation
`X 2 2 = X 3 3` of the accessible algebra `diag(a, b, c, c)`.

The access cut declares the top region as observer region and the
three-block algebra as accessible algebra.  Both cover restrictions map the
accessible algebra into itself and all bottom restrictions are scalar, so
the cut satisfies `restrict_preserves` while the glue escapes.  The receipt
`glue_not_accessible`, the refutation `not_accessibleGlueClosure`, and the
packaged statement `descent_holds_closure_fails` are the deliverables.

Claim boundary: every construction here is a finite conditional structure
over the committed E1 interface, in the sense of `QFT/ObserverAccessCut.lean`.
The countermodel certifies a logical independence between two receipt types
and claims no physical region, coverage, instrument, or continuum object.
-/

namespace OPH.QFT

open Matrix
open OPH.Tower
open EventAlgebra

/-! ## Dimension-four diagonal helpers -/

/-- A scalar multiple of the identity is the constant diagonal matrix. -/
theorem smul_one_eq_diagonal_four (c : ℂ) :
    c • (1 : Matrix (Fin 4) (Fin 4) ℂ) = Matrix.diagonal ![c, c, c, c] := by
  have hvec : ![c, c, c, c] = fun _ : Fin 4 => c := by
    funext i; fin_cases i <;> rfl
  rw [hvec]
  ext i j
  by_cases hij : i = j
  · subst hij
    simp [Matrix.smul_apply]
  · simp [Matrix.smul_apply, Matrix.one_apply_ne hij,
      Matrix.diagonal_apply_ne _ hij]

theorem diagonal_zero_vector_four :
    Matrix.diagonal ![(0 : ℂ), 0, 0, 0] = 0 := by
  have hvec : ![(0 : ℂ), 0, 0, 0] = (fun _ : Fin 4 => (0 : ℂ)) := by
    funext i; fin_cases i <;> rfl
  rw [hvec]
  exact Matrix.diagonal_zero

theorem diagonal_one_vector_four :
    Matrix.diagonal ![(1 : ℂ), 1, 1, 1] = 1 := by
  have hvec : ![(1 : ℂ), 1, 1, 1] = (fun _ : Fin 4 => (1 : ℂ)) := by
    funext i; fin_cases i <;> rfl
  rw [hvec]
  exact Matrix.diagonal_one

/-! ## The two declared middle algebras -/

/-- The two-block diagonal algebra `diag(a, a, b, b)` inside `M₄(ℂ)`, the
regional algebra of the first cover region. -/
noncomputable def gluePairAlgebra :
    StarSubalgebra ℂ (Matrix (Fin 4) (Fin 4) ℂ) where
  carrier := {X | ∃ a b : ℂ, X = Matrix.diagonal ![a, a, b, b]}
  zero_mem' := ⟨0, 0, diagonal_zero_vector_four.symm⟩
  one_mem' := ⟨1, 1, diagonal_one_vector_four.symm⟩
  add_mem' := by
    intro X Y hX hY
    obtain ⟨a, b, rfl⟩ := hX
    obtain ⟨c, d, rfl⟩ := hY
    refine ⟨a + c, b + d, ?_⟩
    rw [Matrix.diagonal_add]
    congr 1
    funext i
    fin_cases i <;> simp
  mul_mem' := by
    intro X Y hX hY
    obtain ⟨a, b, rfl⟩ := hX
    obtain ⟨c, d, rfl⟩ := hY
    refine ⟨a * c, b * d, ?_⟩
    rw [Matrix.diagonal_mul_diagonal]
    congr 1
    funext i
    fin_cases i <;> simp
  algebraMap_mem' := fun c =>
    ⟨c, c, by rw [Algebra.algebraMap_eq_smul_one, smul_one_eq_diagonal_four]⟩
  star_mem' := by
    intro X hX
    obtain ⟨a, b, rfl⟩ := hX
    refine ⟨star a, star b, ?_⟩
    rw [Matrix.star_eq_conjTranspose, Matrix.diagonal_conjTranspose]
    congr 1
    funext i
    fin_cases i <;> simp

theorem mem_gluePairAlgebra_iff {X : Matrix (Fin 4) (Fin 4) ℂ} :
    X ∈ gluePairAlgebra ↔ ∃ a b : ℂ, X = Matrix.diagonal ![a, a, b, b] :=
  Iff.rfl

/-- The three-block diagonal algebra `diag(a, b, c, c)` inside `M₄(ℂ)`, the
regional algebra of the second cover region and the declared accessible
algebra of the countermodel cut. -/
noncomputable def glueAccessibleAlgebra :
    StarSubalgebra ℂ (Matrix (Fin 4) (Fin 4) ℂ) where
  carrier := {X | ∃ a b c : ℂ, X = Matrix.diagonal ![a, b, c, c]}
  zero_mem' := ⟨0, 0, 0, diagonal_zero_vector_four.symm⟩
  one_mem' := ⟨1, 1, 1, diagonal_one_vector_four.symm⟩
  add_mem' := by
    intro X Y hX hY
    obtain ⟨a, b, c, rfl⟩ := hX
    obtain ⟨d, e, f, rfl⟩ := hY
    refine ⟨a + d, b + e, c + f, ?_⟩
    rw [Matrix.diagonal_add]
    congr 1
    funext i
    fin_cases i <;> simp
  mul_mem' := by
    intro X Y hX hY
    obtain ⟨a, b, c, rfl⟩ := hX
    obtain ⟨d, e, f, rfl⟩ := hY
    refine ⟨a * d, b * e, c * f, ?_⟩
    rw [Matrix.diagonal_mul_diagonal]
    congr 1
    funext i
    fin_cases i <;> simp
  algebraMap_mem' := fun c =>
    ⟨c, c, c, by rw [Algebra.algebraMap_eq_smul_one, smul_one_eq_diagonal_four]⟩
  star_mem' := by
    intro X hX
    obtain ⟨a, b, c, rfl⟩ := hX
    refine ⟨star a, star b, star c, ?_⟩
    rw [Matrix.star_eq_conjTranspose, Matrix.diagonal_conjTranspose]
    congr 1
    funext i
    fin_cases i <;> simp

theorem mem_glueAccessibleAlgebra_iff {X : Matrix (Fin 4) (Fin 4) ℂ} :
    X ∈ glueAccessibleAlgebra ↔
      ∃ a b c : ℂ, X = Matrix.diagonal ![a, b, c, c] :=
  Iff.rfl

theorem smul_one_mem_glueAccessibleAlgebra (c : ℂ) :
    c • (1 : Matrix (Fin 4) (Fin 4) ℂ) ∈ glueAccessibleAlgebra :=
  ⟨c, c, c, smul_one_eq_diagonal_four c⟩

/-! ## The corner character and the two cover restrictions -/

theorem gluePairAlgebra_row_zero :
    ∀ X ∈ gluePairAlgebra, ∀ j : Fin 4, j ≠ 0 → X 0 j = 0 := by
  intro X hX j hj
  obtain ⟨a, b, rfl⟩ := hX
  exact Matrix.diagonal_apply_ne _ (Ne.symm hj)

theorem glueAccessibleAlgebra_row_zero :
    ∀ X ∈ glueAccessibleAlgebra, ∀ j : Fin 4, j ≠ 0 → X 0 j = 0 := by
  intro X hX j hj
  obtain ⟨a, b, c, rfl⟩ := hX
  exact Matrix.diagonal_apply_ne _ (Ne.symm hj)

theorem glueDiagonalAlgebra_row_zero :
    ∀ X ∈ (diagonalPartition 4).publicSubalgebra,
      ∀ j : Fin 4, j ≠ 0 → X 0 j = 0 := by
  intro X hX j hj
  obtain ⟨c, rfl⟩ := diagonalPartition_mem_iff.mp hX
  exact Matrix.diagonal_apply_ne _ (Ne.symm hj)

/-- The corner character on a subalgebra with vanishing off-diagonal first
row: evaluation of the `(0, 0)` entry, returned as a scalar element of the
one-block public algebra.  All three bottom restrictions of the
countermodel net are instances of this map. -/
noncomputable def cornerCharacter
    (S : StarSubalgebra ℂ (Matrix (Fin 4) (Fin 4) ℂ))
    (hrow : ∀ X ∈ S, ∀ j : Fin 4, j ≠ 0 → X 0 j = 0) :
    S →⋆ₐ[ℂ] (oneBlockPartition 4).publicSubalgebra where
  toFun X := ⟨(X : Matrix (Fin 4) (Fin 4) ℂ) 0 0 • 1,
    properMeet_smul_one_mem _⟩
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
    have hprod : ((X : Matrix (Fin 4) (Fin 4) ℂ) * Y) 0 0 =
        (X : Matrix (Fin 4) (Fin 4) ℂ) 0 0 *
          (Y : Matrix (Fin 4) (Fin 4) ℂ) 0 0 := by
      rw [Matrix.mul_apply, Fin.sum_univ_four,
        hrow _ X.2 1 (by decide), hrow _ X.2 2 (by decide),
        hrow _ X.2 3 (by decide)]
      ring
    show ((X * Y : S) : Matrix (Fin 4) (Fin 4) ℂ) 0 0 • 1 = _
    simp only [MulMemClass.coe_mul]
    rw [hprod]
    simp [smul_smul, Algebra.mul_smul_comm, mul_comm]
  commutes' c := by
    apply Subtype.ext
    simp [Algebra.algebraMap_eq_smul_one, Matrix.smul_apply]
  map_star' X := by
    apply Subtype.ext
    show (star (X : Matrix (Fin 4) (Fin 4) ℂ)) 0 0 • 1 = _
    simp [Matrix.star_apply, star_smul]

/-- Product entries of two members of the diagonal algebra multiply
entrywise on the diagonal. -/
theorem diagonalMember_mul_entry {X Y : Matrix (Fin 4) (Fin 4) ℂ}
    (hX : X ∈ (diagonalPartition 4).publicSubalgebra)
    (hY : Y ∈ (diagonalPartition 4).publicSubalgebra) (i : Fin 4) :
    (X * Y) i i = X i i * Y i i := by
  obtain ⟨c, rfl⟩ := diagonalPartition_mem_iff.mp hX
  obtain ⟨d, rfl⟩ := diagonalPartition_mem_iff.mp hY
  rw [Matrix.diagonal_mul_diagonal]
  simp

/-- The first cover restriction: the entry pick `(0, 2)` of a diagonal
matrix, filling the two-block pattern `diag(a, a, b, b)`. -/
noncomputable def glueTopToU1 :
    (diagonalPartition 4).publicSubalgebra →⋆ₐ[ℂ] gluePairAlgebra where
  toFun X := ⟨Matrix.diagonal
      ![(X : Matrix (Fin 4) (Fin 4) ℂ) 0 0,
        (X : Matrix (Fin 4) (Fin 4) ℂ) 0 0,
        (X : Matrix (Fin 4) (Fin 4) ℂ) 2 2,
        (X : Matrix (Fin 4) (Fin 4) ℂ) 2 2],
    ⟨(X : Matrix (Fin 4) (Fin 4) ℂ) 0 0,
      (X : Matrix (Fin 4) (Fin 4) ℂ) 2 2, rfl⟩⟩
  map_zero' := by
    apply Subtype.ext
    simp [diagonal_zero_vector_four]
  map_one' := by
    apply Subtype.ext
    simp [Matrix.one_apply_eq, diagonal_one_vector_four]
  map_add' X Y := by
    apply Subtype.ext
    simp only [AddMemClass.coe_add, Matrix.add_apply]
    rw [Matrix.diagonal_add]
    congr 1
    funext i
    fin_cases i <;> simp
  map_mul' X Y := by
    apply Subtype.ext
    simp only [MulMemClass.coe_mul, diagonalMember_mul_entry X.2 Y.2]
    rw [Matrix.diagonal_mul_diagonal]
    congr 1
    funext i
    fin_cases i <;> simp
  commutes' c := by
    apply Subtype.ext
    simp [Algebra.algebraMap_eq_smul_one, Matrix.smul_apply,
      Matrix.one_apply_eq, smul_one_eq_diagonal_four]
  map_star' X := by
    apply Subtype.ext
    simp only [StarMemClass.coe_star, Matrix.star_apply]
    rw [Matrix.star_eq_conjTranspose, Matrix.diagonal_conjTranspose]
    congr 1
    funext i
    fin_cases i <;> simp

/-- The second cover restriction: the entry pick `(0, 1, 3)` of a diagonal
matrix, filling the three-block pattern `diag(a, b, c, c)`. -/
noncomputable def glueTopToU2 :
    (diagonalPartition 4).publicSubalgebra →⋆ₐ[ℂ] glueAccessibleAlgebra where
  toFun X := ⟨Matrix.diagonal
      ![(X : Matrix (Fin 4) (Fin 4) ℂ) 0 0,
        (X : Matrix (Fin 4) (Fin 4) ℂ) 1 1,
        (X : Matrix (Fin 4) (Fin 4) ℂ) 3 3,
        (X : Matrix (Fin 4) (Fin 4) ℂ) 3 3],
    ⟨(X : Matrix (Fin 4) (Fin 4) ℂ) 0 0,
      (X : Matrix (Fin 4) (Fin 4) ℂ) 1 1,
      (X : Matrix (Fin 4) (Fin 4) ℂ) 3 3, rfl⟩⟩
  map_zero' := by
    apply Subtype.ext
    simp [diagonal_zero_vector_four]
  map_one' := by
    apply Subtype.ext
    simp [Matrix.one_apply_eq, diagonal_one_vector_four]
  map_add' X Y := by
    apply Subtype.ext
    simp only [AddMemClass.coe_add, Matrix.add_apply]
    rw [Matrix.diagonal_add]
    congr 1
    funext i
    fin_cases i <;> simp
  map_mul' X Y := by
    apply Subtype.ext
    simp only [MulMemClass.coe_mul, diagonalMember_mul_entry X.2 Y.2]
    rw [Matrix.diagonal_mul_diagonal]
    congr 1
    funext i
    fin_cases i <;> simp
  commutes' c := by
    apply Subtype.ext
    simp [Algebra.algebraMap_eq_smul_one, Matrix.smul_apply,
      Matrix.one_apply_eq, smul_one_eq_diagonal_four]
  map_star' X := by
    apply Subtype.ext
    simp only [StarMemClass.coe_star, Matrix.star_apply]
    rw [Matrix.star_eq_conjTranspose, Matrix.diagonal_conjTranspose]
    congr 1
    funext i
    fin_cases i <;> simp

/-! ## The countermodel region diamond -/

/-- Region labels: a scalar bottom region, the two incomparable cover
regions, and the diagonal top region. -/
inductive GlueRegion : Type
  | bot
  | u1
  | u2
  | top
  deriving DecidableEq, Fintype

namespace GlueRegion

/-- Diamond order: reflexive, the bottom region below everything, the top
region above everything. -/
def le (U V : GlueRegion) : Prop := U = V ∨ U = bot ∨ V = top

instance : DecidableRel le := fun U V =>
  decidable_of_iff (U = V ∨ U = bot ∨ V = top) Iff.rfl

/-- Diamond meet. -/
def meet : GlueRegion → GlueRegion → GlueRegion
  | .top, V => V
  | .bot, _ => .bot
  | .u1, .top => .u1
  | .u1, .u1 => .u1
  | .u1, .u2 => .bot
  | .u1, .bot => .bot
  | .u2, .top => .u2
  | .u2, .u2 => .u2
  | .u2, .u1 => .bot
  | .u2, .bot => .bot

theorem le_refl : ∀ U, le U U := by decide

theorem le_trans : ∀ U V W, le U V → le V W → le U W := by decide

theorem le_antisymm : ∀ U V, le U V → le V U → U = V := by decide

theorem meet_le_left : ∀ U V, le (meet U V) U := by decide

theorem meet_le_right : ∀ U V, le (meet U V) V := by decide

theorem le_meet : ∀ W U V, le W U → le W V → le W (meet U V) := by decide

/-- The two cover regions are incomparable, so the diamond is not a
chain. -/
theorem u1_not_le_u2 : ¬ le u1 u2 := by decide

end GlueRegion

/-- Entry accessor pinning the index type `Fin 4`: elements of the tower's
private algebra carry the regulator-indexed dimension in their type, and
this wrapper lets literal indices elaborate against it. -/
def glueEntry (X : Matrix (Fin 4) (Fin 4) ℂ) (i j : Fin 4) : ℂ := X i j

@[simp]
theorem glueEntry_apply (X : Matrix (Fin 4) (Fin 4) ℂ) (i j : Fin 4) :
    glueEntry X i j = X i j :=
  rfl

/-- The regional algebra assignment of the countermodel net. -/
noncomputable def glueRegionAlgebra :
    GlueRegion → StarSubalgebra ℂ (Matrix (Fin 4) (Fin 4) ℂ)
  | .bot => (oneBlockPartition 4).publicSubalgebra
  | .u1 => gluePairAlgebra
  | .u2 => glueAccessibleAlgebra
  | .top => (diagonalPartition 4).publicSubalgebra

/-- The supplied restriction system: identities on equal regions, the two
entry-pick restrictions from the top, and the corner characters onto the
bottom region. -/
noncomputable def glueRestrict :
    ∀ U V : GlueRegion, GlueRegion.le V U →
      (glueRegionAlgebra U →⋆ₐ[ℂ] glueRegionAlgebra V)
  | .bot, .bot, _ => StarAlgHom.id ℂ _
  | .u1, .u1, _ => StarAlgHom.id ℂ _
  | .u2, .u2, _ => StarAlgHom.id ℂ _
  | .top, .top, _ => StarAlgHom.id ℂ _
  | .u1, .bot, _ => cornerCharacter gluePairAlgebra gluePairAlgebra_row_zero
  | .u2, .bot, _ =>
      cornerCharacter glueAccessibleAlgebra glueAccessibleAlgebra_row_zero
  | .top, .bot, _ =>
      cornerCharacter ((diagonalPartition 4).publicSubalgebra)
        glueDiagonalAlgebra_row_zero
  | .top, .u1, _ => glueTopToU1
  | .top, .u2, _ => glueTopToU2
  | .bot, .u1, h => absurd h (by decide)
  | .bot, .u2, h => absurd h (by decide)
  | .bot, .top, h => absurd h (by decide)
  | .u1, .u2, h => absurd h (by decide)
  | .u1, .top, h => absurd h (by decide)
  | .u2, .u1, h => absurd h (by decide)
  | .u2, .top, h => absurd h (by decide)

/-- The countermodel tower: the constant tower of the one-block partition
on dimension four, so the public algebra is the scalar span of the
identity and sits inside every regional algebra. -/
noncomputable def glueTower : ConsensusTower Unit :=
  ConsensusTower.constantConsensusTower (oneBlockPartition 4) properMeetState

/-- The countermodel net over the diamond: scalar bottom, the two-block
and three-block middle regions, and the diagonal top, with the asymmetric
entry-pick restrictions.  The declared disjointness relation is empty and
the repair maps are the identity. -/
noncomputable def glueNet : FiniteCausalObserverNet glueTower where
  Region := fun _ => GlueRegion
  regionFintype := fun _ => inferInstance
  regionNonempty := fun _ => ⟨GlueRegion.bot⟩
  regionLE := fun _ => GlueRegion.le
  regionLE_refl := fun _ U => GlueRegion.le_refl U
  regionLE_trans := fun _ {U V W} hUV hVW =>
    GlueRegion.le_trans U V W hUV hVW
  regionLE_antisymm := fun _ {U V} hUV hVU =>
    GlueRegion.le_antisymm U V hUV hVU
  overlap := fun _ => GlueRegion.meet
  overlap_le_left := fun _ U V => GlueRegion.meet_le_left U V
  overlap_le_right := fun _ U V => GlueRegion.meet_le_right U V
  le_overlap := fun _ {W U V} hWU hWV => GlueRegion.le_meet W U V hWU hWV
  disjoint := fun _ _ _ => False
  disjoint_symm := fun _ {_ _} h => h.elim
  disjoint_irrefl := fun _ _ h => h
  localAlgebra := fun _ => glueRegionAlgebra
  isotony := by
    intro r U V hUV
    rcases hUV with rfl | rfl | rfl
    · exact le_rfl
    · exact oneBlockPartition_publicSubalgebra_le _
    · cases U with
      | bot => exact oneBlockPartition_publicSubalgebra_le _
      | u1 =>
          intro X hX
          obtain ⟨a, b, rfl⟩ := hX
          exact diagonalPartition_mem_iff.mpr ⟨![a, a, b, b], rfl⟩
      | u2 =>
          intro X hX
          obtain ⟨a, b, c, rfl⟩ := hX
          exact diagonalPartition_mem_iff.mpr ⟨![a, b, c, c], rfl⟩
      | top => exact le_rfl
  locality := fun _ {_ _} h => h.elim
  restrict := fun _ {U V} h => glueRestrict U V h
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
        | u1 =>
            apply Subtype.ext
            exact oneBlockPartition_corner_eval X (0 : Fin 4)
        | u2 =>
            apply Subtype.ext
            exact oneBlockPartition_corner_eval X (0 : Fin 4)
        | top =>
            apply Subtype.ext
            exact oneBlockPartition_corner_eval X (0 : Fin 4)
    | u1 =>
        cases V with
        | bot => exact absurd hUV (by decide)
        | u1 => rfl
        | u2 => exact absurd hUV (by decide)
        | top =>
            apply Subtype.ext
            obtain ⟨a, b, h⟩ := X.2
            show Matrix.diagonal ![glueEntry X.1 0 0, glueEntry X.1 0 0,
                glueEntry X.1 2 2, glueEntry X.1 2 2] = X.1
            rw [h]
            simp
    | u2 =>
        cases V with
        | bot => exact absurd hUV (by decide)
        | u1 => exact absurd hUV (by decide)
        | u2 => rfl
        | top =>
            apply Subtype.ext
            obtain ⟨a, b, c, h⟩ := X.2
            show Matrix.diagonal ![glueEntry X.1 0 0, glueEntry X.1 1 1,
                glueEntry X.1 3 3, glueEntry X.1 3 3] = X.1
            rw [h]
            simp
    | top =>
        cases V with
        | bot => exact absurd hUV (by decide)
        | u1 => exact absurd hUV (by decide)
        | u2 => exact absurd hUV (by decide)
        | top => rfl
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

/-- Decidable equality transported along the definitional identity of the
net's region type with `GlueRegion`, registered so that finite-set
literals and membership tests over `glueNet.Region r` elaborate. -/
instance glueNetRegionDecidableEq (r : Unit) :
    DecidableEq (glueNet.Region r) :=
  inferInstanceAs (DecidableEq GlueRegion)

/-! ## The countermodel access cut -/

/-- The countermodel access cut: the observer region is the top region and
the accessible algebra is the three-block algebra `diag(a, b, c, c)`.  Both
entry-pick restrictions map the accessible algebra into itself and every
bottom restriction is scalar, so restriction stability holds without any
appeal to the cover structure. -/
noncomputable def glueAccessCut : ObserverAccessCut glueTower glueNet () where
  observerRegion := fun _ => GlueRegion.top
  accessibleAlgebra := fun _ => glueAccessibleAlgebra
  accessible_le_region := by
    intro o X hX
    obtain ⟨a, b, c, rfl⟩ := hX
    exact diagonalPartition_mem_iff.mpr ⟨![a, b, c, c], rfl⟩
  public_le_accessible := fun _ =>
    oneBlockPartition_publicSubalgebra_le _
  restrict_preserves := by
    intro o U V hU hVU X hX
    cases U with
    | bot =>
        cases V with
        | bot => exact hX
        | u1 =>
            exact absurd hVU
              (by decide : ¬ GlueRegion.le GlueRegion.u1 GlueRegion.bot)
        | u2 =>
            exact absurd hVU
              (by decide : ¬ GlueRegion.le GlueRegion.u2 GlueRegion.bot)
        | top =>
            exact absurd hVU
              (by decide : ¬ GlueRegion.le GlueRegion.top GlueRegion.bot)
    | u1 =>
        cases V with
        | bot => exact smul_one_mem_glueAccessibleAlgebra _
        | u1 => exact hX
        | u2 =>
            exact absurd hVU
              (by decide : ¬ GlueRegion.le GlueRegion.u2 GlueRegion.u1)
        | top =>
            exact absurd hVU
              (by decide : ¬ GlueRegion.le GlueRegion.top GlueRegion.u1)
    | u2 =>
        cases V with
        | bot => exact smul_one_mem_glueAccessibleAlgebra _
        | u1 =>
            exact absurd hVU
              (by decide : ¬ GlueRegion.le GlueRegion.u1 GlueRegion.u2)
        | u2 => exact hX
        | top =>
            exact absurd hVU
              (by decide : ¬ GlueRegion.le GlueRegion.top GlueRegion.u2)
    | top =>
        cases V with
        | bot => exact smul_one_mem_glueAccessibleAlgebra _
        | u1 =>
            exact ⟨glueEntry X.1 0 0, glueEntry X.1 0 0,
              glueEntry X.1 2 2, rfl⟩
        | u2 =>
            exact ⟨glueEntry X.1 0 0, glueEntry X.1 1 1,
              glueEntry X.1 3 3, rfl⟩
        | top => exact hX

/-! ## The declared family and its unique descent -/

/-- The declared two-member subregion family `{u1, u2}` of the top
region. -/
noncomputable def gluePairCover : glueNet.FiniteCover () GlueRegion.top where
  regions := {GlueRegion.u1, GlueRegion.u2}
  nonempty := ⟨GlueRegion.u1, Finset.mem_insert_self _ _⟩
  subregion := fun _ _ => Or.inr (Or.inr rfl)

/-- Ambient unique descent for the declared family: the two entry picks
`(0, 2)` and `(0, 1, 3)` jointly read every diagonal entry, so the glue is
unique, and the single overlap equation at entry `0` matches the single
shared pick, so every compatible family glues. -/
theorem gluePairCover_hasUniqueDescent : gluePairCover.HasUniqueDescent := by
  intro F
  have hu1 : GlueRegion.u1 ∈ gluePairCover.regions := by decide
  have hu2 : GlueRegion.u2 ∈ gluePairCover.regions := by decide
  obtain ⟨a, b, hs1⟩ :=
    mem_gluePairAlgebra_iff.mp (F.localSection GlueRegion.u1 hu1).2
  obtain ⟨c, d, e, hs2⟩ :=
    mem_glueAccessibleAlgebra_iff.mp (F.localSection GlueRegion.u2 hu2).2
  have hcompat :=
    F.overlap_compatible GlueRegion.u1 hu1 GlueRegion.u2 hu2
  have hcorner :
      glueEntry (F.localSection GlueRegion.u1 hu1).1 0 0 •
          (1 : Matrix (Fin 4) (Fin 4) ℂ) =
        glueEntry (F.localSection GlueRegion.u2 hu2).1 0 0 •
          (1 : Matrix (Fin 4) (Fin 4) ℂ) :=
    congrArg Subtype.val hcompat
  have hac : a = c := by
    rw [hs1, hs2] at hcorner
    have h := congrFun (congrFun hcorner 0) 0
    simpa using h
  refine ⟨⟨Matrix.diagonal ![a, d, b, e],
    diagonalPartition_mem_iff.mpr ⟨![a, d, b, e], rfl⟩⟩, ?_, ?_⟩
  · intro U hU
    fin_cases hU
    · apply Subtype.ext
      show Matrix.diagonal ![a, a, b, b] = _
      exact hs1.symm
    · apply Subtype.ext
      show Matrix.diagonal ![a, d, e, e] = _
      rw [hac]
      exact hs2.symm
  · intro Y hY
    obtain ⟨v, hv⟩ := diagonalPartition_mem_iff.mp Y.2
    have h1 := hY GlueRegion.u1 hu1
    have h2 := hY GlueRegion.u2 hu2
    have e00 : glueEntry Y.1 0 0 =
        glueEntry (F.localSection GlueRegion.u1 hu1).1 0 0 :=
      congrArg (fun Z => glueEntry Z.1 0 0) h1
    have e22 : glueEntry Y.1 2 2 =
        glueEntry (F.localSection GlueRegion.u1 hu1).1 2 2 :=
      congrArg (fun Z => glueEntry Z.1 2 2) h1
    have e11 : glueEntry Y.1 1 1 =
        glueEntry (F.localSection GlueRegion.u2 hu2).1 1 1 :=
      congrArg (fun Z => glueEntry Z.1 1 1) h2
    have e33 : glueEntry Y.1 3 3 =
        glueEntry (F.localSection GlueRegion.u2 hu2).1 3 3 :=
      congrArg (fun Z => glueEntry Z.1 3 3) h2
    rw [hs1] at e00 e22
    rw [hs2] at e11 e33
    rw [← hv] at e00 e11 e22 e33
    simp at e00 e11 e22 e33
    apply Subtype.ext
    rw [← hv]
    show Matrix.diagonal v = Matrix.diagonal ![a, d, b, e]
    congr 1
    funext i
    fin_cases i
    · simpa using e00
    · simpa using e11
    · simpa using e22
    · simpa using e33

/-- The declared descent packet: the pair family at the top region and the
singleton identity families elsewhere. -/
noncomputable def glueDescent : glueNet.ObserverNetDescent where
  declaredCover {r} {W} C :=
    match W, C with
    | .top, C => C = gluePairCover
    | .bot, C => C = glueNet.singletonCover r .bot
    | .u1, C => C = glueNet.singletonCover r .u1
    | .u2, C => C = glueNet.singletonCover r .u2
  cover_exists := by
    intro r W
    cases W with
    | bot => exact ⟨glueNet.singletonCover r .bot, rfl⟩
    | u1 => exact ⟨glueNet.singletonCover r .u1, rfl⟩
    | u2 => exact ⟨glueNet.singletonCover r .u2, rfl⟩
    | top => exact ⟨gluePairCover, rfl⟩
  unique_descent := by
    intro r W C hC
    cases W with
    | bot =>
        subst hC
        exact glueNet.singletonCover_hasUniqueDescent r _
    | u1 =>
        subst hC
        exact glueNet.singletonCover_hasUniqueDescent r _
    | u2 =>
        subst hC
        exact glueNet.singletonCover_hasUniqueDescent r _
    | top =>
        subst hC
        exact gluePairCover_hasUniqueDescent

/-! ## The compatible accessible family and its inaccessible glue -/

/-- The first section: `diag(0, 0, 1, 1)`, a member of the two-block
regional algebra and of the accessible algebra. -/
noncomputable def glueSectionU1 : glueNet.localAlgebra () GlueRegion.u1 :=
  ⟨Matrix.diagonal ![0, 0, 1, 1], ⟨0, 1, rfl⟩⟩

/-- The second section: `diag(0, 1, 0, 0)`, a member of the three-block
regional algebra, which is the accessible algebra itself. -/
noncomputable def glueSectionU2 : glueNet.localAlgebra () GlueRegion.u2 :=
  ⟨Matrix.diagonal ![0, 1, 0, 0], ⟨0, 1, 0, rfl⟩⟩

/-- The compatible family on the declared pair cover: both sections
evaluate to `0` at the shared corner entry. -/
noncomputable def glueFamily : gluePairCover.CompatibleFamily where
  localSection := fun U hU =>
    match U, hU with
    | .u1, _ => glueSectionU1
    | .u2, _ => glueSectionU2
    | .bot, hU => absurd hU (by decide)
    | .top, hU => absurd hU (by decide)
  overlap_compatible := by
    intro U hU V hV
    fin_cases hU <;> fin_cases hV <;> rfl

/-- The unique ambient glue of the family: `diag(0, 1, 1, 0)`, read off
entrywise from the two picks. -/
noncomputable def glueGlobalSection : glueNet.localAlgebra () GlueRegion.top :=
  ⟨Matrix.diagonal ![0, 1, 1, 0],
    diagonalPartition_mem_iff.mpr ⟨![0, 1, 1, 0], rfl⟩⟩

/-- The glue receipt: `diag(0, 1, 1, 0)` restricts to both declared
sections. -/
theorem gluePairCover_globalizes :
    gluePairCover.Globalizes glueFamily glueGlobalSection := by
  intro U hU
  fin_cases hU <;> exact Subtype.ext rfl

/-- Both sections of the family are accessible to every observer. -/
theorem glueFamily_sections_accessible (o : glueTower.Observer ()) :
    ∀ U (hU : U ∈ gluePairCover.regions),
      (glueFamily.localSection U hU :
          ConsensusTower.PrivateAlgebra glueTower ()) ∈
        glueAccessCut.accessibleAlgebra o := by
  intro U hU
  fin_cases hU
  · exact ⟨0, 0, 1, rfl⟩
  · exact ⟨0, 1, 0, rfl⟩

/-- The countermodel receipt: the unique ambient glue of the accessible
family violates the defining equation `X 2 2 = X 3 3` of the accessible
algebra. -/
theorem glue_not_accessible (o : glueTower.Observer ()) :
    (glueGlobalSection : ConsensusTower.PrivateAlgebra glueTower ()) ∉
      glueAccessCut.accessibleAlgebra o := by
  intro hmem
  obtain ⟨a, b, c, h⟩ := hmem
  have h22 : (1 : ℂ) = c := congrFun (congrFun h (2 : Fin 4)) (2 : Fin 4)
  have h33 : (0 : ℂ) = c := congrFun (congrFun h (3 : Fin 4)) (3 : Fin 4)
  exact one_ne_zero (h22.trans h33.symm)

/-- The selected glue of the descent packet on the declared family is the
exhibited global section. -/
theorem glueDescent_glue_eq :
    glueDescent.glue gluePairCover rfl glueFamily = glueGlobalSection :=
  (glueDescent.glue_unique gluePairCover rfl glueFamily glueGlobalSection
    gluePairCover_globalizes).symm

/-! ## Refutation of the closure receipt and the independence statement -/

/-- The accessible-glue receipt fails on the countermodel cut for every
observer: the declared pair family, the accessible compatible family, and
the inaccessible glue instantiate the negation. -/
theorem not_accessibleGlueClosure (o : glueTower.Observer ()) :
    ¬ glueAccessCut.AccessibleGlueClosure o glueDescent := by
  intro G
  exact glue_not_accessible o
    (G.glue_accessible (Or.inl rfl) gluePairCover rfl glueFamily
      (glueFamily_sections_accessible o) glueGlobalSection
      gluePairCover_globalizes)

/-- Independence, packaged: the declared family carries ambient unique
descent while the accessible-glue receipt fails for every observer of the
countermodel cut.  Together with `accessibleGlueClosure_singleton` in
`QFT/ObserverAccessCut.lean`, which inhabits the receipt on the singleton
packet of every cut, this places `AccessibleGlueClosure` strictly between
the ambient descent packet and its negation: neither receipt implies the
other. -/
theorem descent_holds_closure_fails :
    gluePairCover.HasUniqueDescent ∧
      ∀ o : glueTower.Observer (),
        ¬ glueAccessCut.AccessibleGlueClosure o glueDescent :=
  ⟨gluePairCover_hasUniqueDescent, not_accessibleGlueClosure⟩

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.gluePairCover_hasUniqueDescent
#print axioms OPH.QFT.gluePairCover_globalizes
#print axioms OPH.QFT.glueFamily_sections_accessible
#print axioms OPH.QFT.glue_not_accessible
#print axioms OPH.QFT.glueDescent_glue_eq
#print axioms OPH.QFT.not_accessibleGlueClosure
#print axioms OPH.QFT.descent_holds_closure_fails
