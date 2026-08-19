import QFT.CarrierJoinTransport

/-!
# Record-layer net morphisms for the hinge-fibered carrier join

`QFT.CarrierJoinTransport` embeds the diagonal (record) function layers of
the two committed carriers into one 338-dimensional pushout with intertwined
anchoring and walk-step transport, and its committed boundary states that no
regional-net morphism and no tower morphism is constructed.  This module
supplies the net morphism at the record layer and proves exactly where the
operator-layer structure fails to follow.

**Part 1: the record-layer regional net on each committed carrier.**  The
committed regional nets are the two-slot diamonds `twoSlotNet` (on the
86/88 carrier `PairIndex`) and `twoSlotNet247` (on the 86/247 carrier
`PairIndex247`), bundled as the committed premise instances
`designatedPairPremises` and `supportDisjointPairPremises` (premise rows
PR-32, PR-33, PR-34).  For each diamond region, the record layer of the
committed regional algebra is defined as a typed function subspace of the
carrier: constants at the bottom, first-coordinate functions on the left
slot, second-coordinate functions on the right slot, everything at the top.
The characterization theorems (`diagonal_mem_twoSlotAlgebra_iff`,
`diagonal_mem_twoSlotAlgebra247_iff`) prove that these subspaces are exactly
the diagonal restrictions of the committed regional algebras: a diagonal
matrix lies in the committed regional algebra of a region precisely when its
diagonal function lies in the record layer of that region.  Isotony of both
record-layer nets is proved, with dimensions 1, 13, 14, 182 on the 86/88
side and 1, 13, 13, 169 on the 86/247 side.

**Part 2: the induced net on the join.**  The induced region family
`JoinRegion` carries one region per committed region of each carrier (the
image of its record layer under the committed embedding) plus the hinge (the
13-dimensional shared observer-86 image).  The family is isotone for the
declared order `JoinRegion.le` (a decidable partial order proved by kernel
computation), embedding images of nested regions nest
(`embed88_image_nested`, `embed247_image_nested`), and the hinge lies below
the left slot image on both sides, through either embedding
(`hingeLayer_le_map_embed88_left`, `hingeLayer_eq_247`).  The embeddings are
injective, so every induced region has the dimension of its record layer.

**Part 3: the morphism theorems.**  Each committed embedding maps the
record-layer net of its carrier to the induced net on the join region by
region and functorially: the image assignment is definitional
(`embed88_maps_recordLayer`), inclusions are preserved, and dimensions are
preserved.  The committed conditional expectations restrict to the record
layer exactly: on every region, the committed expectation of the committed
net carries diagonal matrices to diagonal matrices, and the induced map on
diagonal functions is the typed record-layer expectation
(`twoSlotExpect_diagonal`, `twoSlotExpect247_diagonal`), which lands in the
record layer of its region and fixes it pointwise.  The matched pair of
left-slot record expectations descends to the join as one hinge projection
`joinHingeProject` that commutes with both embeddings, projects onto the
hinge region, fixes it, and preserves the hinge restriction.  The committed
join evolution restricts region-compatibly wherever the committed evolution
does: every induced region is carried onto itself by every join step
(`joinStep_maps_joinRegionLayer`), mirroring the committed
`stepEvolve_maps_regions` at the record layer.

**Part 4: what is NOT preserved, exactly.**  Two sharp negatives delimit
the record layer.  First, the gluing submodule is not multiplicatively
closed (`hingeGraph_not_mul_closed`): an explicit hinge generator has a
componentwise square outside the graph, so the quotient multiplication is
ill-defined (`join_quotient_mul_ill_defined`) and the join carries no
algebra structure through this pushout.  The operator-layer net morphism
therefore cannot be obtained by this construction; that residual stays
open.  Second, the record layer is commutative (`recordDiagonal_commute`), while
the committed left regional algebra of the 86/247 net is noncommutative
(`twoSlotAlgebra247_left_noncommutative`): the commutant and locality
structure of the committed operator-layer net (which separates disjoint
slot pairs from noncommuting region pairs) collapses at the record layer,
where every pair of regions commutes elementwise.

**The composed receipt** (`joinNetMorphism_receipt`) states the
characterization, isotony, morphism, expectation-restriction, hinge-descent,
evolution-compatibility, and delimitation clauses as one conjunction, typed
against the committed premise instances `designatedPairPremises` and
`supportDisjointPairPremises`, whose fields are the registered premise rows
PR-32 (Cartesian joint carrier), PR-33 (diamond region map), and PR-34
(slot split).  Physical time, causality, and spacetime attachment stay open
under PR-52; the source-selected time-slice content of PR-58 stays open.

**What is NOT claimed.**  Record layer only: the operator-level regional
algebras do not transport through the join, and the full operator-algebra
net morphism stays open, delimited by the negatives above.  No tower-level
morphism, no spacetime, no physical time, no causal reading, no continuum
limit, and no uniqueness of the induced net is claimed.  The declared order
`JoinRegion.le` is one declared choice of induced region order; the
fiber-uniform hinge sections stay a declared normalization; the rich-fibre
and source regional nets live on other carriers and are outside this
module's scope.  OL-C6 stays partial with the net-morphism residual
narrowed to the operator layer.
-/

set_option autoImplicit false
set_option maxRecDepth 65536

namespace OPH.QFT

open Matrix
open Kronecker

noncomputable section

/-! ## Typed generators of the record layers -/

/-- The constant-function inclusion of the scalars into a function
carrier. -/
def constLift (γ : Type*) : ℂ →ₗ[ℂ] (γ → ℂ) where
  toFun c := fun _ => c
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

theorem constLift_apply (γ : Type*) (c : ℂ) (p : γ) : constLift γ c p = c := rfl

/-- Pullback of functions along a map of index sets. -/
def pullbackFn {γ δ : Type*} (f : γ → δ) : (δ → ℂ) →ₗ[ℂ] (γ → ℂ) where
  toFun h := h ∘ f
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

theorem pullbackFn_apply {γ δ : Type*} (f : γ → δ) (h : δ → ℂ) (p : γ) :
    pullbackFn f h p = h (f p) := rfl

theorem constLift_injective (γ : Type*) [Nonempty γ] :
    Function.Injective (constLift γ) := by
  intro a b hab
  exact congrFun hab (Classical.choice ‹Nonempty γ›)

theorem pullbackFn_injective {γ δ : Type*} (f : γ → δ)
    (hf : Function.Surjective f) : Function.Injective (pullbackFn f) := by
  intro a b hab
  funext d
  obtain ⟨c, rfl⟩ := hf d
  exact congrFun hab c

/-! ## The record-layer regional nets of the committed carriers -/

/-- The record layer of the committed 86/88 regional net: for each region
of the committed diamond, the diagonal (record) function subspace of its
committed regional algebra. -/
def recordLayer88 : TwoSlotRegion → Submodule ℂ CarrierFn88
  | TwoSlotRegion.bot => LinearMap.range (constLift PairIndex)
  | TwoSlotRegion.left =>
      LinearMap.range (pullbackFn (Prod.fst : PairIndex → Fin 13))
  | TwoSlotRegion.right =>
      LinearMap.range (pullbackFn (Prod.snd : PairIndex → Fin 14))
  | TwoSlotRegion.top => ⊤

/-- The record layer of the committed 86/247 regional net. -/
def recordLayer247 : TwoSlotRegion → Submodule ℂ CarrierFn247
  | TwoSlotRegion.bot => LinearMap.range (constLift PairIndex247)
  | TwoSlotRegion.left =>
      LinearMap.range (pullbackFn (Prod.fst : PairIndex247 → Fin 13))
  | TwoSlotRegion.right =>
      LinearMap.range (pullbackFn (Prod.snd : PairIndex247 → Fin 13))
  | TwoSlotRegion.top => ⊤

/-- **Isotony of the 86/88 record-layer net.** -/
theorem recordLayer88_isotone {U V : TwoSlotRegion}
    (hUV : TwoSlotRegion.le U V = true) :
    recordLayer88 U ≤ recordLayer88 V := by
  cases U <;> cases V <;> (try (simp [TwoSlotRegion.le] at hUV))
  · exact le_rfl
  · rintro v ⟨c, rfl⟩
    exact ⟨fun _ => c, rfl⟩
  · rintro v ⟨c, rfl⟩
    exact ⟨fun _ => c, rfl⟩
  · exact le_top
  · exact le_rfl
  · exact le_top
  · exact le_rfl
  · exact le_top
  · exact le_rfl

/-- **Isotony of the 86/247 record-layer net.** -/
theorem recordLayer247_isotone {U V : TwoSlotRegion}
    (hUV : TwoSlotRegion.le U V = true) :
    recordLayer247 U ≤ recordLayer247 V := by
  cases U <;> cases V <;> (try (simp [TwoSlotRegion.le] at hUV))
  · exact le_rfl
  · rintro v ⟨c, rfl⟩
    exact ⟨fun _ => c, rfl⟩
  · rintro v ⟨c, rfl⟩
    exact ⟨fun _ => c, rfl⟩
  · exact le_top
  · exact le_rfl
  · exact le_top
  · exact le_rfl
  · exact le_top
  · exact le_rfl

/-! ## The record layers are exactly the diagonal restrictions of the
committed regional algebras -/

/-- A diagonal Kronecker factor against the identity is the diagonal of
the first-coordinate pullback. -/
theorem diagonal_kronecker_one {α β : Type*} [DecidableEq α] [DecidableEq β]
    (d : α → ℂ) :
    (Matrix.diagonal d) ⊗ₖ (1 : Matrix β β ℂ) =
      Matrix.diagonal (fun p : α × β => d p.1) := by
  rw [← Matrix.diagonal_one, Matrix.diagonal_kronecker_diagonal]
  congr 1
  funext p
  exact mul_one _

/-- The identity Kronecker factor against a diagonal is the diagonal of
the second-coordinate pullback. -/
theorem one_kronecker_diagonal {α β : Type*} [DecidableEq α] [DecidableEq β]
    (d : β → ℂ) :
    (1 : Matrix α α ℂ) ⊗ₖ Matrix.diagonal d =
      Matrix.diagonal (fun p : α × β => d p.2) := by
  rw [← Matrix.diagonal_one, Matrix.diagonal_kronecker_diagonal]
  congr 1
  funext p
  exact one_mul _

/-- Diagonal membership in the scalar bottom algebra is constancy. -/
theorem diagonal_mem_bot_iff {γ : Type*} [Fintype γ] [DecidableEq γ]
    (v : γ → ℂ) :
    Matrix.diagonal v ∈ (⊥ : StarSubalgebra ℂ (Matrix γ γ ℂ)) ↔
      ∃ c : ℂ, (fun _ : γ => c) = v := by
  constructor
  · intro hv
    obtain ⟨c, hc⟩ := StarSubalgebra.mem_bot.mp hv
    refine ⟨c, ?_⟩
    rw [Algebra.algebraMap_eq_smul_one, ← Matrix.diagonal_one,
      ← Matrix.diagonal_smul] at hc
    have hfun := Matrix.diagonal_injective hc
    funext p
    have hp := congrFun hfun p
    simpa using hp
  · rintro ⟨c, rfl⟩
    refine StarSubalgebra.mem_bot.mpr ⟨c, ?_⟩
    rw [Algebra.algebraMap_eq_smul_one, ← Matrix.diagonal_one,
      ← Matrix.diagonal_smul]
    congr 1
    funext p
    simp

/-- Diagonal membership in the left slot algebra is first-coordinate
dependence. -/
theorem diagonal_mem_slotLeft_iff {α β : Type*} [Fintype α] [Fintype β]
    [DecidableEq α] [DecidableEq β] (v : α × β → ℂ) :
    Matrix.diagonal v ∈ (slotLeft β (α := α)).range ↔
      ∃ h : α → ℂ, (fun p : α × β => h p.1) = v := by
  constructor
  · rintro ⟨A, hA⟩
    refine ⟨fun a => A a a, ?_⟩
    funext p
    have hp : (A ⊗ₖ (1 : Matrix β β ℂ)) p p = Matrix.diagonal v p p :=
      congrFun (congrFun hA p) p
    rw [Matrix.kroneckerMap_apply, Matrix.diagonal_apply_eq,
      Matrix.one_apply_eq, mul_one] at hp
    exact hp
  · rintro ⟨h, rfl⟩
    refine ⟨Matrix.diagonal h, ?_⟩
    show Matrix.diagonal h ⊗ₖ (1 : Matrix β β ℂ) = _
    exact diagonal_kronecker_one h

/-- Diagonal membership in the right slot algebra is second-coordinate
dependence. -/
theorem diagonal_mem_slotRight_iff {α β : Type*} [Fintype α] [Fintype β]
    [DecidableEq α] [DecidableEq β] (v : α × β → ℂ) :
    Matrix.diagonal v ∈ (slotRight α (β := β)).range ↔
      ∃ h : β → ℂ, (fun p : α × β => h p.2) = v := by
  constructor
  · rintro ⟨B, hB⟩
    refine ⟨fun b => B b b, ?_⟩
    funext p
    have hp : ((1 : Matrix α α ℂ) ⊗ₖ B) p p = Matrix.diagonal v p p :=
      congrFun (congrFun hB p) p
    rw [Matrix.kroneckerMap_apply, Matrix.diagonal_apply_eq,
      Matrix.one_apply_eq, one_mul] at hp
    exact hp
  · rintro ⟨h, rfl⟩
    refine ⟨Matrix.diagonal h, ?_⟩
    show (1 : Matrix α α ℂ) ⊗ₖ Matrix.diagonal h = _
    exact one_kronecker_diagonal h

/-- **The 86/88 record layer is exactly the diagonal restriction of the
committed regional algebra**, region by region. -/
theorem diagonal_mem_twoSlotAlgebra_iff (U : TwoSlotRegion)
    (v : CarrierFn88) :
    Matrix.diagonal v ∈ twoSlotAlgebra U ↔ v ∈ recordLayer88 U := by
  cases U
  · constructor
    · intro hv
      obtain ⟨c, hc⟩ := (diagonal_mem_bot_iff v).mp hv
      exact ⟨c, hc⟩
    · rintro ⟨c, rfl⟩
      exact (diagonal_mem_bot_iff _).mpr ⟨c, rfl⟩
  · constructor
    · intro hv
      obtain ⟨h, hh⟩ := (diagonal_mem_slotLeft_iff v).mp hv
      exact ⟨h, hh⟩
    · rintro ⟨h, rfl⟩
      exact (diagonal_mem_slotLeft_iff _).mpr ⟨h, rfl⟩
  · constructor
    · intro hv
      obtain ⟨h, hh⟩ := (diagonal_mem_slotRight_iff v).mp hv
      exact ⟨h, hh⟩
    · rintro ⟨h, rfl⟩
      exact (diagonal_mem_slotRight_iff _).mpr ⟨h, rfl⟩
  · constructor <;> intro <;> trivial

/-- **The 86/247 record layer is exactly the diagonal restriction of the
committed regional algebra**, region by region. -/
theorem diagonal_mem_twoSlotAlgebra247_iff (U : TwoSlotRegion)
    (v : CarrierFn247) :
    Matrix.diagonal v ∈ twoSlotAlgebra247 U ↔ v ∈ recordLayer247 U := by
  cases U
  · constructor
    · intro hv
      obtain ⟨c, hc⟩ := (diagonal_mem_bot_iff v).mp hv
      exact ⟨c, hc⟩
    · rintro ⟨c, rfl⟩
      exact (diagonal_mem_bot_iff _).mpr ⟨c, rfl⟩
  · constructor
    · intro hv
      obtain ⟨h, hh⟩ := (diagonal_mem_slotLeft_iff v).mp hv
      exact ⟨h, hh⟩
    · rintro ⟨h, rfl⟩
      exact (diagonal_mem_slotLeft_iff _).mpr ⟨h, rfl⟩
  · constructor
    · intro hv
      obtain ⟨h, hh⟩ := (diagonal_mem_slotRight_iff v).mp hv
      exact ⟨h, hh⟩
    · rintro ⟨h, rfl⟩
      exact (diagonal_mem_slotRight_iff _).mpr ⟨h, rfl⟩
  · constructor <;> intro <;> trivial

/-! ## Record-layer dimensions -/

/-- The record-layer dimension table of the 86/88 net. -/
def recordDim88 : TwoSlotRegion → ℕ
  | TwoSlotRegion.bot => 1
  | TwoSlotRegion.left => 13
  | TwoSlotRegion.right => 14
  | TwoSlotRegion.top => 182

/-- The record-layer dimension table of the 86/247 net. -/
def recordDim247 : TwoSlotRegion → ℕ
  | TwoSlotRegion.bot => 1
  | TwoSlotRegion.left => 13
  | TwoSlotRegion.right => 13
  | TwoSlotRegion.top => 169

theorem recordLayer88_finrank (U : TwoSlotRegion) :
    Module.finrank ℂ (recordLayer88 U) = recordDim88 U := by
  cases U
  · rw [show recordLayer88 TwoSlotRegion.bot =
        LinearMap.range (constLift PairIndex) from rfl,
      LinearMap.finrank_range_of_inj (constLift_injective PairIndex)]
    exact Module.finrank_self ℂ
  · rw [show recordLayer88 TwoSlotRegion.left =
        LinearMap.range (pullbackFn (Prod.fst : PairIndex → Fin 13)) from rfl,
      LinearMap.finrank_range_of_inj
        (pullbackFn_injective _ (fun i => ⟨(i, 0), rfl⟩))]
    simp [Module.finrank_fintype_fun_eq_card, recordDim88]
  · rw [show recordLayer88 TwoSlotRegion.right =
        LinearMap.range (pullbackFn (Prod.snd : PairIndex → Fin 14)) from rfl,
      LinearMap.finrank_range_of_inj
        (pullbackFn_injective _ (fun j => ⟨(0, j), rfl⟩))]
    simp [Module.finrank_fintype_fun_eq_card, recordDim88]
  · rw [show recordLayer88 TwoSlotRegion.top =
        (⊤ : Submodule ℂ CarrierFn88) from rfl, finrank_top]
    simp [Module.finrank_fintype_fun_eq_card, recordDim88]

theorem recordLayer247_finrank (U : TwoSlotRegion) :
    Module.finrank ℂ (recordLayer247 U) = recordDim247 U := by
  cases U
  · rw [show recordLayer247 TwoSlotRegion.bot =
        LinearMap.range (constLift PairIndex247) from rfl,
      LinearMap.finrank_range_of_inj (constLift_injective PairIndex247)]
    exact Module.finrank_self ℂ
  · rw [show recordLayer247 TwoSlotRegion.left =
        LinearMap.range (pullbackFn (Prod.fst : PairIndex247 → Fin 13))
        from rfl,
      LinearMap.finrank_range_of_inj
        (pullbackFn_injective _ (fun i => ⟨(i, 0), rfl⟩))]
    simp [Module.finrank_fintype_fun_eq_card, recordDim247]
  · rw [show recordLayer247 TwoSlotRegion.right =
        LinearMap.range (pullbackFn (Prod.snd : PairIndex247 → Fin 13))
        from rfl,
      LinearMap.finrank_range_of_inj
        (pullbackFn_injective _ (fun j => ⟨(0, j), rfl⟩))]
    simp [Module.finrank_fintype_fun_eq_card, recordDim247]
  · rw [show recordLayer247 TwoSlotRegion.top =
        (⊤ : Submodule ℂ CarrierFn247) from rfl, finrank_top]
    simp [Module.finrank_fintype_fun_eq_card, recordDim247]

/-! ## The induced regions on the join -/

/-- The induced region family of the join: one region per committed
region of each carrier, plus the observer-86 hinge. -/
inductive JoinRegion : Type
  | side88 : TwoSlotRegion → JoinRegion
  | side247 : TwoSlotRegion → JoinRegion
  | hinge : JoinRegion
deriving DecidableEq, Fintype

namespace JoinRegion

/-- The declared induced order: the committed diamond order within each
side, and the hinge below the left slot and the top of both sides. -/
def le : JoinRegion → JoinRegion → Bool
  | side88 U, side88 V => TwoSlotRegion.le U V
  | side247 U, side247 V => TwoSlotRegion.le U V
  | hinge, hinge => true
  | hinge, side88 V => TwoSlotRegion.le TwoSlotRegion.left V
  | hinge, side247 V => TwoSlotRegion.le TwoSlotRegion.left V
  | _, _ => false

theorem le_refl' : ∀ R : JoinRegion, le R R = true := by decide

theorem le_trans' : ∀ R S T : JoinRegion,
    le R S = true → le S T = true → le R T = true := by decide

theorem le_antisymm' : ∀ R S : JoinRegion,
    le R S = true → le S R = true → R = S := by decide

end JoinRegion

/-- The induced regional subspace assignment on the join: the embedding
image of each record layer, region by region, plus the hinge image. -/
def joinRegionLayer : JoinRegion → Submodule ℂ JoinCarrier
  | JoinRegion.side88 U => (recordLayer88 U).map (embed88 : CarrierFn88 →ₗ[ℂ] JoinCarrier)
  | JoinRegion.side247 U => (recordLayer247 U).map (embed247 : CarrierFn247 →ₗ[ℂ] JoinCarrier)
  | JoinRegion.hinge => LinearMap.range (embed88 ∘ₗ sect88)

/-- The fiber-uniform hinge section lands in the left-slot record layer
of the 86/88 carrier. -/
theorem sect88_mem_recordLayer_left (h : HingeFn) :
    sect88 h ∈ recordLayer88 TwoSlotRegion.left :=
  ⟨fun i => h i / 14, rfl⟩

/-- The fiber-uniform hinge section lands in the left-slot record layer
of the 86/247 carrier. -/
theorem sect247_mem_recordLayer_left (h : HingeFn) :
    sect247 h ∈ recordLayer247 TwoSlotRegion.left :=
  ⟨fun i => h i / 13, rfl⟩

/-- The hinge region lies below the 86/88 left-slot image. -/
theorem hingeLayer_le_map_embed88_left :
    LinearMap.range (embed88 ∘ₗ sect88) ≤
      (recordLayer88 TwoSlotRegion.left).map embed88 := by
  rintro z ⟨h, rfl⟩
  exact ⟨sect88 h, sect88_mem_recordLayer_left h, rfl⟩

/-- The hinge image is the same subspace through either embedding. -/
theorem hingeLayer_eq_247 :
    LinearMap.range (embed88 ∘ₗ sect88) =
      LinearMap.range (embed247 ∘ₗ sect247) := by
  apply le_antisymm
  · rintro z ⟨h, rfl⟩
    exact ⟨h, (embed_agree_on_hinge h).symm⟩
  · rintro z ⟨h, rfl⟩
    exact ⟨h, embed_agree_on_hinge h⟩

/-- The hinge region lies below the 86/247 left-slot image. -/
theorem hingeLayer_le_map_embed247_left :
    LinearMap.range (embed88 ∘ₗ sect88) ≤
      (recordLayer247 TwoSlotRegion.left).map embed247 := by
  rw [hingeLayer_eq_247]
  rintro z ⟨h, rfl⟩
  exact ⟨sect247 h, sect247_mem_recordLayer_left h, rfl⟩

/-- **Nesting of embedded regions, 86/88 side.**  Images of nested
record layers nest. -/
theorem embed88_image_nested {U V : TwoSlotRegion}
    (hUV : TwoSlotRegion.le U V = true) :
    (recordLayer88 U).map embed88 ≤ (recordLayer88 V).map embed88 :=
  Submodule.map_mono (recordLayer88_isotone hUV)

/-- **Nesting of embedded regions, 86/247 side.** -/
theorem embed247_image_nested {U V : TwoSlotRegion}
    (hUV : TwoSlotRegion.le U V = true) :
    (recordLayer247 U).map embed247 ≤ (recordLayer247 V).map embed247 :=
  Submodule.map_mono (recordLayer247_isotone hUV)

/-- **Isotony of the induced net on the join.** -/
theorem joinRegionLayer_isotone : ∀ R S : JoinRegion,
    JoinRegion.le R S = true → joinRegionLayer R ≤ joinRegionLayer S := by
  intro R S hRS
  cases R <;> cases S <;> (try (simp [JoinRegion.le] at hRS))
  · exact embed88_image_nested hRS
  · exact embed247_image_nested hRS
  · exact le_trans hingeLayer_le_map_embed88_left (embed88_image_nested hRS)
  · exact le_trans hingeLayer_le_map_embed247_left (embed247_image_nested hRS)
  · exact le_rfl

/-- **The embeddings are region-by-region net morphisms**: the image
assignment of the induced net is definitional on each side. -/
theorem embed88_maps_recordLayer (U : TwoSlotRegion) :
    (recordLayer88 U).map embed88 = joinRegionLayer (JoinRegion.side88 U) :=
  rfl

theorem embed247_maps_recordLayer (U : TwoSlotRegion) :
    (recordLayer247 U).map embed247 = joinRegionLayer (JoinRegion.side247 U) :=
  rfl

/-- Injective embeddings preserve region dimensions. -/
theorem finrank_map_embed88 (p : Submodule ℂ CarrierFn88) :
    Module.finrank ℂ (p.map embed88) = Module.finrank ℂ p :=
  ((Submodule.equivMapOfInjective embed88 embed88_injective p).finrank_eq).symm

theorem finrank_map_embed247 (p : Submodule ℂ CarrierFn247) :
    Module.finrank ℂ (p.map embed247) = Module.finrank ℂ p :=
  ((Submodule.equivMapOfInjective embed247 embed247_injective p).finrank_eq).symm

theorem joinRegionLayer_side88_finrank (U : TwoSlotRegion) :
    Module.finrank ℂ (joinRegionLayer (JoinRegion.side88 U)) =
      recordDim88 U := by
  rw [← embed88_maps_recordLayer, finrank_map_embed88, recordLayer88_finrank]

theorem joinRegionLayer_side247_finrank (U : TwoSlotRegion) :
    Module.finrank ℂ (joinRegionLayer (JoinRegion.side247 U)) =
      recordDim247 U := by
  rw [← embed247_maps_recordLayer, finrank_map_embed247,
    recordLayer247_finrank]

theorem joinRegionLayer_hinge_finrank :
    Module.finrank ℂ (joinRegionLayer JoinRegion.hinge) = 13 :=
  hingeImage_finrank

/-! ## The committed conditional expectations at the record layer -/

/-- The record-layer average: the normalized total, as a constant
function.  This is the record layer of the committed scalar
expectation. -/
def recordAvg (γ : Type*) [Fintype γ] : (γ → ℂ) →ₗ[ℂ] (γ → ℂ) where
  toFun v := fun _ => (∑ p, v p) / (Fintype.card γ : ℂ)
  map_add' u v := by
    funext q
    show (∑ p, (u p + v p)) / _ = _
    rw [Finset.sum_add_distrib, add_div]
    rfl
  map_smul' c v := by
    funext q
    show (∑ p, c * v p) / _ = _
    rw [← Finset.mul_sum, mul_div_assoc]
    rfl

/-- Right-slot marginalization of the 86/88 carrier. -/
def margRight88 : CarrierFn88 →ₗ[ℂ] (Fin 14 → ℂ) where
  toFun v := fun j => ∑ i : Fin 13, v (i, j)
  map_add' u v := by
    funext j
    simp [Finset.sum_add_distrib]
  map_smul' c v := by
    funext j
    simp [Finset.mul_sum]

/-- The fiber-uniform right-slot section of the 86/88 carrier. -/
def sectRight88 : (Fin 14 → ℂ) →ₗ[ℂ] CarrierFn88 where
  toFun h := fun p => h p.2 / 13
  map_add' u v := by
    funext p
    simp [add_div]
  map_smul' c h := by
    funext p
    simp [mul_div_assoc]

/-- Right-slot marginalization of the 86/247 carrier. -/
def margRight247 : CarrierFn247 →ₗ[ℂ] (Fin 13 → ℂ) where
  toFun v := fun j => ∑ i : Fin 13, v (i, j)
  map_add' u v := by
    funext j
    simp [Finset.sum_add_distrib]
  map_smul' c v := by
    funext j
    simp [Finset.mul_sum]

/-- The fiber-uniform right-slot section of the 86/247 carrier. -/
def sectRight247 : (Fin 13 → ℂ) →ₗ[ℂ] CarrierFn247 where
  toFun h := fun p => h p.2 / 13
  map_add' u v := by
    funext p
    simp [add_div]
  map_smul' c h := by
    funext p
    simp [mul_div_assoc]

/-- The record-layer conditional expectations of the 86/88 net. -/
def recordExpect88 : TwoSlotRegion → (CarrierFn88 →ₗ[ℂ] CarrierFn88)
  | TwoSlotRegion.bot => recordAvg PairIndex
  | TwoSlotRegion.left => sect88 ∘ₗ marg88
  | TwoSlotRegion.right => sectRight88 ∘ₗ margRight88
  | TwoSlotRegion.top => LinearMap.id

/-- The record-layer conditional expectations of the 86/247 net. -/
def recordExpect247 : TwoSlotRegion → (CarrierFn247 →ₗ[ℂ] CarrierFn247)
  | TwoSlotRegion.bot => recordAvg PairIndex247
  | TwoSlotRegion.left => sect247 ∘ₗ marg247
  | TwoSlotRegion.right => sectRight247 ∘ₗ margRight247
  | TwoSlotRegion.top => LinearMap.id

/-- The committed left-slot expectation carries diagonals to diagonals.
The diagonal images of the committed partial traces are the committed
`ptraceSnd_diagonal` and `ptraceFst_diagonal` of
`QFT.SourceCorrelationCapstone`. -/
theorem leftSlotExpectation_diagonal {α β : Type*} [Fintype α] [Fintype β]
    [DecidableEq α] [DecidableEq β] (v : α × β → ℂ) :
    leftSlotExpectation (Matrix.diagonal v) =
      Matrix.diagonal
        ((Fintype.card β : ℂ)⁻¹ • fun p : α × β => ∑ b, v (p.1, b)) := by
  have h : leftSlotExpectation (Matrix.diagonal v) =
      (Fintype.card β : ℂ)⁻¹ •
        (ptraceSnd (Matrix.diagonal v) ⊗ₖ (1 : Matrix β β ℂ)) := rfl
  rw [h, ptraceSnd_diagonal, diagonal_kronecker_one, Matrix.diagonal_smul]

/-- The committed right-slot expectation carries diagonals to diagonals. -/
theorem rightSlotExpectation_diagonal {α β : Type*} [Fintype α] [Fintype β]
    [DecidableEq α] [DecidableEq β] (v : α × β → ℂ) :
    rightSlotExpectation (Matrix.diagonal v) =
      Matrix.diagonal
        ((Fintype.card α : ℂ)⁻¹ • fun p : α × β => ∑ a, v (a, p.2)) := by
  have h : rightSlotExpectation (Matrix.diagonal v) =
      (Fintype.card α : ℂ)⁻¹ •
        ((1 : Matrix α α ℂ) ⊗ₖ OPH.Locality.ptraceFst (Matrix.diagonal v)) :=
    rfl
  rw [h, ptraceFst_diagonal, one_kronecker_diagonal, Matrix.diagonal_smul]

/-- The committed scalar expectation carries diagonals to diagonals. -/
theorem scalarExpectation_diagonal {γ : Type*} [Fintype γ] [DecidableEq γ]
    (v : γ → ℂ) :
    scalarExpectation γ (Matrix.diagonal v) =
      Matrix.diagonal
        ((Fintype.card γ : ℂ)⁻¹ • fun _ : γ => ∑ p, v p) := by
  have h : scalarExpectation γ (Matrix.diagonal v) =
      (Fintype.card γ : ℂ)⁻¹ •
        (Matrix.diagonal v).trace • (1 : Matrix γ γ ℂ) := rfl
  rw [h, Matrix.trace_diagonal, ← Matrix.diagonal_one, ← Matrix.diagonal_smul,
    ← Matrix.diagonal_smul]
  congr 1
  funext p
  simp

/-- **The committed 86/88 conditional expectations restrict to the record
layer.**  On every region, the committed expectation of the committed net
maps diagonal matrices to diagonal matrices, acting as the typed
record-layer expectation. -/
theorem twoSlotExpect_diagonal (U : TwoSlotRegion) (v : CarrierFn88) :
    twoSlotExpect U (Matrix.diagonal v) =
      Matrix.diagonal (recordExpect88 U v) := by
  cases U
  · rw [show twoSlotExpect TwoSlotRegion.bot = scalarExpectation PairIndex
      from rfl, scalarExpectation_diagonal]
    congr 1
    funext p
    show (Fintype.card PairIndex : ℂ)⁻¹ * (∑ q, v q) =
      (∑ q, v q) / (Fintype.card PairIndex : ℂ)
    rw [div_eq_mul_inv, mul_comm]
  · rw [show twoSlotExpect TwoSlotRegion.left =
      (leftSlotExpectation : Matrix PairIndex PairIndex ℂ →ₗ[ℂ]
        Matrix PairIndex PairIndex ℂ) from rfl, leftSlotExpectation_diagonal]
    congr 1
    funext p
    show (Fintype.card (Fin 14) : ℂ)⁻¹ * (∑ b, v (p.1, b)) =
      (∑ b, v (p.1, b)) / 14
    rw [div_eq_mul_inv, mul_comm, Fintype.card_fin]
    norm_num
  · rw [show twoSlotExpect TwoSlotRegion.right =
      (rightSlotExpectation : Matrix PairIndex PairIndex ℂ →ₗ[ℂ]
        Matrix PairIndex PairIndex ℂ) from rfl, rightSlotExpectation_diagonal]
    congr 1
    funext p
    show (Fintype.card (Fin 13) : ℂ)⁻¹ * (∑ a, v (a, p.2)) =
      (∑ a, v (a, p.2)) / 13
    rw [div_eq_mul_inv, mul_comm, Fintype.card_fin]
    norm_num
  · rfl

/-- **The committed 86/247 conditional expectations restrict to the
record layer.** -/
theorem twoSlotExpect247_diagonal (U : TwoSlotRegion) (v : CarrierFn247) :
    twoSlotExpect247 U (Matrix.diagonal v) =
      Matrix.diagonal (recordExpect247 U v) := by
  cases U
  · rw [show twoSlotExpect247 TwoSlotRegion.bot =
      scalarExpectation PairIndex247 from rfl, scalarExpectation_diagonal]
    congr 1
    funext p
    show (Fintype.card PairIndex247 : ℂ)⁻¹ * (∑ q, v q) =
      (∑ q, v q) / (Fintype.card PairIndex247 : ℂ)
    rw [div_eq_mul_inv, mul_comm]
  · rw [show twoSlotExpect247 TwoSlotRegion.left =
      (leftSlotExpectation : Matrix PairIndex247 PairIndex247 ℂ →ₗ[ℂ]
        Matrix PairIndex247 PairIndex247 ℂ) from rfl,
      leftSlotExpectation_diagonal]
    congr 1
    funext p
    show (Fintype.card (Fin 13) : ℂ)⁻¹ * (∑ b, v (p.1, b)) =
      (∑ b, v (p.1, b)) / 13
    rw [div_eq_mul_inv, mul_comm, Fintype.card_fin]
    norm_num
  · rw [show twoSlotExpect247 TwoSlotRegion.right =
      (rightSlotExpectation : Matrix PairIndex247 PairIndex247 ℂ →ₗ[ℂ]
        Matrix PairIndex247 PairIndex247 ℂ) from rfl,
      rightSlotExpectation_diagonal]
    congr 1
    funext p
    show (Fintype.card (Fin 13) : ℂ)⁻¹ * (∑ a, v (a, p.2)) =
      (∑ a, v (a, p.2)) / 13
    rw [div_eq_mul_inv, mul_comm, Fintype.card_fin]
    norm_num
  · rfl

/-- The record-layer expectations land in their record layers. -/
theorem recordExpect88_mem (U : TwoSlotRegion) (v : CarrierFn88) :
    recordExpect88 U v ∈ recordLayer88 U := by
  cases U
  · exact ⟨(∑ p, v p) / (Fintype.card PairIndex : ℂ), rfl⟩
  · exact ⟨fun i => marg88 v i / 14, rfl⟩
  · exact ⟨fun j => margRight88 v j / 13, rfl⟩
  · trivial

theorem recordExpect247_mem (U : TwoSlotRegion) (v : CarrierFn247) :
    recordExpect247 U v ∈ recordLayer247 U := by
  cases U
  · exact ⟨(∑ p, v p) / (Fintype.card PairIndex247 : ℂ), rfl⟩
  · exact ⟨fun i => marg247 v i / 13, rfl⟩
  · exact ⟨fun j => margRight247 v j / 13, rfl⟩
  · trivial

/-- The record-layer expectations fix their record layers pointwise. -/
theorem recordExpect88_fixes (U : TwoSlotRegion) {v : CarrierFn88}
    (hv : v ∈ recordLayer88 U) : recordExpect88 U v = v := by
  cases U
  · obtain ⟨c, rfl⟩ := hv
    funext p
    show (∑ _q : PairIndex, c) / (Fintype.card PairIndex : ℂ) = c
    rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    exact mul_div_cancel_left₀ c
      (Nat.cast_ne_zero.mpr Fintype.card_ne_zero)
  · obtain ⟨h, rfl⟩ := hv
    funext p
    show (∑ _j : Fin 14, h p.1) / 14 = h p.1
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    norm_num
  · obtain ⟨h, rfl⟩ := hv
    funext p
    show (∑ _i : Fin 13, h p.2) / 13 = h p.2
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    norm_num
  · rfl

theorem recordExpect247_fixes (U : TwoSlotRegion) {v : CarrierFn247}
    (hv : v ∈ recordLayer247 U) : recordExpect247 U v = v := by
  cases U
  · obtain ⟨c, rfl⟩ := hv
    funext p
    show (∑ _q : PairIndex247, c) / (Fintype.card PairIndex247 : ℂ) = c
    rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    exact mul_div_cancel_left₀ c
      (Nat.cast_ne_zero.mpr Fintype.card_ne_zero)
  · obtain ⟨h, rfl⟩ := hv
    funext p
    show (∑ _j : Fin 13, h p.1) / 13 = h p.1
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    norm_num
  · obtain ⟨h, rfl⟩ := hv
    funext p
    show (∑ _i : Fin 13, h p.2) / 13 = h p.2
    rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    norm_num
  · rfl

/-! ## The hinge projection descends to the join -/

theorem hingeGraph_hingeProject_invariant :
    hingeGraph ≤ hingeGraph.comap
      ((sect88 ∘ₗ marg88).prodMap (sect247 ∘ₗ marg247)) := by
  rintro x ⟨h, rfl⟩
  rw [Submodule.mem_comap, hingeGraph_gen_apply, LinearMap.prodMap_apply]
  refine ⟨h, ?_⟩
  rw [hingeGraph_gen_apply]
  have h1 : (sect88 ∘ₗ marg88) (sect88 h) = sect88 h := by
    rw [LinearMap.comp_apply, marg88_sect88]
  have h2 : (sect247 ∘ₗ marg247) (-(sect247 h)) = -(sect247 h) := by
    rw [map_neg, LinearMap.comp_apply, marg247_sect247]
  rw [h1, h2]

/-- The join hinge projection: the descent of the matched pair of
left-slot record expectations. -/
def joinHingeProject : JoinCarrier →ₗ[ℂ] JoinCarrier :=
  hingeGraph.mapQ hingeGraph
    ((sect88 ∘ₗ marg88).prodMap (sect247 ∘ₗ marg247))
    hingeGraph_hingeProject_invariant

/-- **The hinge projection restricts to the committed left expectation at
the record layer through the 86/88 embedding.** -/
theorem joinHingeProject_embed88 (v : CarrierFn88) :
    joinHingeProject (embed88 v) =
      embed88 (recordExpect88 TwoSlotRegion.left v) := by
  rw [embed88_apply]
  show hingeGraph.mapQ hingeGraph _ _ (Submodule.Quotient.mk (v, 0)) = _
  rw [Submodule.mapQ_apply, LinearMap.prodMap_apply, map_zero, embed88_apply]
  rfl

/-- **The hinge projection restricts to the committed left expectation at
the record layer through the 86/247 embedding.** -/
theorem joinHingeProject_embed247 (w : CarrierFn247) :
    joinHingeProject (embed247 w) =
      embed247 (recordExpect247 TwoSlotRegion.left w) := by
  rw [embed247_apply]
  show hingeGraph.mapQ hingeGraph _ _ (Submodule.Quotient.mk (0, w)) = _
  rw [Submodule.mapQ_apply, LinearMap.prodMap_apply, map_zero, embed247_apply]
  rfl

/-- The hinge projection lands in the hinge region. -/
theorem joinHingeProject_mem_hinge (z : JoinCarrier) :
    joinHingeProject z ∈ joinRegionLayer JoinRegion.hinge := by
  obtain ⟨⟨v, w⟩, rfl⟩ := Submodule.Quotient.mk_surjective hingeGraph z
  show hingeGraph.mapQ hingeGraph _ _ (Submodule.Quotient.mk (v, w)) ∈ _
  rw [Submodule.mapQ_apply, LinearMap.prodMap_apply]
  have hsplit : (Submodule.Quotient.mk
      ((sect88 ∘ₗ marg88) v, (sect247 ∘ₗ marg247) w) : JoinCarrier) =
      embed88 (sect88 (marg88 v)) + embed247 (sect247 (marg247 w)) := by
    rw [embed88_apply, embed247_apply, ← Submodule.Quotient.mk_add,
      Prod.mk_add_mk, add_zero, zero_add]
    rfl
  rw [hsplit, ← embed_agree_on_hinge (marg247 w)]
  refine ⟨marg88 v + marg247 w, ?_⟩
  rw [LinearMap.comp_apply, map_add, map_add]

/-- The hinge projection fixes the hinge region pointwise. -/
theorem joinHingeProject_fixes_hinge (h : HingeFn) :
    joinHingeProject (embed88 (sect88 h)) = embed88 (sect88 h) := by
  rw [joinHingeProject_embed88]
  congr 1
  show sect88 (marg88 (sect88 h)) = sect88 h
  rw [marg88_sect88]

/-- The hinge projection is idempotent. -/
theorem joinHingeProject_idem (z : JoinCarrier) :
    joinHingeProject (joinHingeProject z) = joinHingeProject z := by
  obtain ⟨h, hh⟩ := joinHingeProject_mem_hinge z
  rw [← hh]
  exact joinHingeProject_fixes_hinge h

/-- The hinge projection preserves the hinge restriction of the join. -/
theorem hingeRestrict_joinHingeProject (z : JoinCarrier) :
    hingeRestrict (joinHingeProject z) = hingeRestrict z := by
  obtain ⟨⟨v, w⟩, rfl⟩ := Submodule.Quotient.mk_surjective hingeGraph z
  show hingeRestrict (hingeGraph.mapQ hingeGraph _ _
    (Submodule.Quotient.mk (v, w))) = _
  rw [Submodule.mapQ_apply, LinearMap.prodMap_apply, hingeRestrict_mk,
    hingeRestrict_mk]
  show marg88 (sect88 (marg88 v)) + marg247 (sect247 (marg247 w)) = _
  rw [marg88_sect88, marg247_sect247]

/-! ## The join evolution restricts region-compatibly -/

/-- Precomposition with the inverse of a permutation, as a linear map on
functions. -/
def permAct {γ : Type*} (e : γ ≃ γ) : (γ → ℂ) →ₗ[ℂ] (γ → ℂ) :=
  pullbackFn (e.symm : γ → γ)

theorem permAct_apply {γ : Type*} (e : γ ≃ γ) (h : γ → ℂ) (p : γ) :
    permAct e h p = h (e.symm p) := rfl

theorem permAct_range_top {γ : Type*} (e : γ ≃ γ) :
    LinearMap.range (permAct e) = ⊤ := by
  rw [LinearMap.range_eq_top]
  intro v
  refine ⟨fun p => v (e p), ?_⟩
  funext p
  rw [permAct_apply]
  exact congrArg v (Equiv.apply_symm_apply e p)

theorem evolveFn88_comp_constLift (t : Fin 31) :
    (evolveFn88 t) ∘ₗ constLift PairIndex = constLift PairIndex := by
  refine LinearMap.ext fun c => ?_
  funext p
  rfl

theorem evolveFn88_comp_pullFst (t : Fin 31) :
    (evolveFn88 t) ∘ₗ pullbackFn (Prod.fst : PairIndex → Fin 13) =
      pullbackFn (Prod.fst : PairIndex → Fin 13) ∘ₗ
        permAct (walkSwapLeft t.castSucc) := by
  refine LinearMap.ext fun h => ?_
  funext p
  show h ((walkStepEquiv88 t).symm p).1 =
    h ((walkSwapLeft t.castSucc).symm p.1)
  rw [walkStepEquiv88_symm_apply]

theorem evolveFn88_comp_pullSnd (t : Fin 31) :
    (evolveFn88 t) ∘ₗ pullbackFn (Prod.snd : PairIndex → Fin 14) =
      pullbackFn (Prod.snd : PairIndex → Fin 14) ∘ₗ
        permAct (walkSwap88 t) := by
  refine LinearMap.ext fun h => ?_
  funext p
  show h ((walkStepEquiv88 t).symm p).2 = h ((walkSwap88 t).symm p.2)
  rw [walkStepEquiv88_symm_apply]

theorem evolveFn88_range_top (t : Fin 31) :
    LinearMap.range (evolveFn88 t) = ⊤ := by
  rw [LinearMap.range_eq_top]
  intro v
  refine ⟨fun p => v (walkStepEquiv88 t p), ?_⟩
  funext p
  rw [evolveFn88_apply]
  exact congrArg v (Equiv.apply_symm_apply _ p)

theorem evolveFn247_comp_constLift (t : Fin 31) :
    (evolveFn247 t) ∘ₗ constLift PairIndex247 = constLift PairIndex247 := by
  refine LinearMap.ext fun c => ?_
  funext p
  rfl

theorem evolveFn247_comp_pullFst (t : Fin 31) :
    (evolveFn247 t) ∘ₗ pullbackFn (Prod.fst : PairIndex247 → Fin 13) =
      pullbackFn (Prod.fst : PairIndex247 → Fin 13) ∘ₗ
        permAct (walkSwapLeft t.castSucc) := by
  refine LinearMap.ext fun h => ?_
  funext p
  show h ((walkStepEquiv t.castSucc).symm p).1 =
    h ((walkSwapLeft t.castSucc).symm p.1)
  rw [walkStepEquiv247_symm_apply]

theorem evolveFn247_comp_pullSnd (t : Fin 31) :
    (evolveFn247 t) ∘ₗ pullbackFn (Prod.snd : PairIndex247 → Fin 13) =
      pullbackFn (Prod.snd : PairIndex247 → Fin 13) ∘ₗ
        permAct (walkSwapRight t.castSucc) := by
  refine LinearMap.ext fun h => ?_
  funext p
  show h ((walkStepEquiv t.castSucc).symm p).2 =
    h ((walkSwapRight t.castSucc).symm p.2)
  rw [walkStepEquiv247_symm_apply]

theorem evolveFn247_range_top (t : Fin 31) :
    LinearMap.range (evolveFn247 t) = ⊤ := by
  rw [LinearMap.range_eq_top]
  intro v
  refine ⟨fun p => v (walkStepEquiv t.castSucc p), ?_⟩
  funext p
  rw [evolveFn247_apply]
  exact congrArg v (Equiv.apply_symm_apply _ p)

/-- **The committed 86/88 walk transport preserves every record layer.** -/
theorem evolveFn88_map_recordLayer (t : Fin 31) (U : TwoSlotRegion) :
    (recordLayer88 U).map (evolveFn88 t) = recordLayer88 U := by
  cases U
  · show Submodule.map _ (LinearMap.range (constLift PairIndex)) = _
    rw [← LinearMap.range_comp, evolveFn88_comp_constLift]
    rfl
  · show Submodule.map _
      (LinearMap.range (pullbackFn (Prod.fst : PairIndex → Fin 13))) = _
    rw [← LinearMap.range_comp, evolveFn88_comp_pullFst,
      LinearMap.range_comp, permAct_range_top, Submodule.map_top]
    rfl
  · show Submodule.map _
      (LinearMap.range (pullbackFn (Prod.snd : PairIndex → Fin 14))) = _
    rw [← LinearMap.range_comp, evolveFn88_comp_pullSnd,
      LinearMap.range_comp, permAct_range_top, Submodule.map_top]
    rfl
  · show Submodule.map _ (⊤ : Submodule ℂ CarrierFn88) = _
    rw [Submodule.map_top, evolveFn88_range_top]
    rfl

/-- **The committed 86/247 walk transport preserves every record
layer.** -/
theorem evolveFn247_map_recordLayer (t : Fin 31) (U : TwoSlotRegion) :
    (recordLayer247 U).map (evolveFn247 t) = recordLayer247 U := by
  cases U
  · show Submodule.map _ (LinearMap.range (constLift PairIndex247)) = _
    rw [← LinearMap.range_comp, evolveFn247_comp_constLift]
    rfl
  · show Submodule.map _
      (LinearMap.range (pullbackFn (Prod.fst : PairIndex247 → Fin 13))) = _
    rw [← LinearMap.range_comp, evolveFn247_comp_pullFst,
      LinearMap.range_comp, permAct_range_top, Submodule.map_top]
    rfl
  · show Submodule.map _
      (LinearMap.range (pullbackFn (Prod.snd : PairIndex247 → Fin 13))) = _
    rw [← LinearMap.range_comp, evolveFn247_comp_pullSnd,
      LinearMap.range_comp, permAct_range_top, Submodule.map_top]
    rfl
  · show Submodule.map _ (⊤ : Submodule ℂ CarrierFn247) = _
    rw [Submodule.map_top, evolveFn247_range_top]
    rfl

theorem joinStep_comp_embed88 (t : Fin 31) :
    (joinStep t) ∘ₗ embed88 = embed88 ∘ₗ evolveFn88 t :=
  LinearMap.ext (joinStep_embed88 t)

theorem joinStep_comp_embed247 (t : Fin 31) :
    (joinStep t) ∘ₗ embed247 = embed247 ∘ₗ evolveFn247 t :=
  LinearMap.ext (joinStep_embed247 t)

theorem joinStep_comp_hinge (t : Fin 31) :
    (joinStep t) ∘ₗ (embed88 ∘ₗ sect88) =
      (embed88 ∘ₗ sect88) ∘ₗ hingeStep t := by
  refine LinearMap.ext fun h => ?_
  show joinStep t (embed88 (sect88 h)) = embed88 (sect88 (hingeStep t h))
  rw [joinStep_embed88, evolveFn88_sect88]

theorem hingeStep_range_top (t : Fin 31) :
    LinearMap.range (hingeStep t) = ⊤ := by
  rw [LinearMap.range_eq_top]
  intro v
  refine ⟨fun i => v (walkSwapLeft t.castSucc i), ?_⟩
  funext i
  rw [hingeStep_apply]
  exact congrArg v (Equiv.apply_symm_apply _ i)

/-- **The committed join evolution restricts region-compatibly.**  Every
induced region of the join is carried onto itself by every join step,
mirroring the committed `stepEvolve_maps_regions` at the record layer. -/
theorem joinStep_maps_joinRegionLayer (t : Fin 31) (R : JoinRegion) :
    (joinRegionLayer R).map (joinStep t) = joinRegionLayer R := by
  cases R with
  | side88 U =>
      show ((recordLayer88 U).map embed88).map (joinStep t) = _
      rw [← Submodule.map_comp, joinStep_comp_embed88, Submodule.map_comp,
        evolveFn88_map_recordLayer]
      rfl
  | side247 U =>
      show ((recordLayer247 U).map embed247).map (joinStep t) = _
      rw [← Submodule.map_comp, joinStep_comp_embed247, Submodule.map_comp,
        evolveFn247_map_recordLayer]
      rfl
  | hinge =>
      show (LinearMap.range (embed88 ∘ₗ sect88)).map (joinStep t) = _
      rw [← LinearMap.range_comp, joinStep_comp_hinge,
        LinearMap.range_comp, hingeStep_range_top, Submodule.map_top]
      rfl

/-! ## What is NOT preserved: the sharp negatives -/

/-- **The gluing submodule is not multiplicatively closed.**  The hinge
generator of the constant hinge datum has a componentwise square outside
the graph: the two fiber normalizations 14 and 13 square incompatibly. -/
theorem hingeGraph_not_mul_closed :
    ∃ x : CarrierFn88 × CarrierFn247, x ∈ hingeGraph ∧ x * x ∉ hingeGraph := by
  refine ⟨(sect88 fun _ => 1, -(sect247 fun _ => 1)),
    ⟨fun _ => 1, rfl⟩, ?_⟩
  rintro ⟨g, hg⟩
  rw [hingeGraph_gen_apply] at hg
  have h1 : sect88 g = (sect88 fun _ => 1) * (sect88 fun _ => 1) :=
    congrArg Prod.fst hg
  have hg0 : g 0 / 14 = 1 / 14 * (1 / 14) := by
    have := congrFun h1 ((0 : Fin 13), (0 : Fin 14))
    rw [sect88_apply] at this
    exact this
  have hgval : g 0 = 1 / 14 := by
    have h14 : (14 : ℂ) ≠ 0 := by norm_num
    rw [div_eq_iff h14] at hg0
    rw [hg0]
    norm_num
  have h2 : -(sect247 g) =
      (-(sect247 fun _ => 1)) * (-(sect247 fun _ => 1)) :=
    congrArg Prod.snd hg
  have h2p := congrFun h2 ((0 : Fin 13), (0 : Fin 13))
  have hL : (-(sect247 g)) ((0 : Fin 13), (0 : Fin 13)) = -(g 0 / 13) := rfl
  have hR : ((-(sect247 fun _ => 1)) * (-(sect247 fun _ => 1)))
      ((0 : Fin 13), (0 : Fin 13)) = (-(1 / 13 : ℂ)) * (-(1 / 13 : ℂ)) := rfl
  rw [hL, hR, hgval] at h2p
  norm_num at h2p

/-- **The quotient multiplication is ill-defined on the join.**  Two
representatives of the same join class have componentwise products in
different classes, so the pushout carries no algebra structure and the
operator-layer net morphism cannot be obtained through it. -/
theorem join_quotient_mul_ill_defined :
    ∃ x y : CarrierFn88 × CarrierFn247,
      (Submodule.Quotient.mk x : JoinCarrier) = Submodule.Quotient.mk y ∧
      (Submodule.Quotient.mk (x * x) : JoinCarrier) ≠
        Submodule.Quotient.mk (y * y) := by
  obtain ⟨x, hx, hxx⟩ := hingeGraph_not_mul_closed
  refine ⟨x, 0, ?_, ?_⟩
  · rw [Submodule.Quotient.eq]
    simpa using hx
  · intro h
    rw [Submodule.Quotient.eq] at h
    apply hxx
    simpa using h

/-- **The record layer is commutative.**  All diagonal matrices commute,
so the commutant and locality distinctions of the committed operator
layer are invisible at the record layer. -/
theorem recordDiagonal_commute {γ : Type*} [Fintype γ] [DecidableEq γ]
    (v w : γ → ℂ) :
    Matrix.diagonal v * Matrix.diagonal w =
      Matrix.diagonal w * Matrix.diagonal v := by
  rw [Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  congr 1
  funext i
  exact mul_comm _ _

/-- **The committed left regional algebra is noncommutative at the
operator layer.**  Two left-slot matrix units fail to commute, while
their record-layer counterparts (all diagonal functions) commute by
`recordDiagonal_commute`: the operator-layer structure erased by the record
restriction is nonempty. -/
theorem twoSlotAlgebra247_left_noncommutative :
    ∃ X Y : Matrix PairIndex247 PairIndex247 ℂ,
      X ∈ twoSlotAlgebra247 TwoSlotRegion.left ∧
      Y ∈ twoSlotAlgebra247 TwoSlotRegion.left ∧
      X * Y ≠ Y * X := by
  refine ⟨slotLeft (Fin 13) (Matrix.single 0 1 1),
    slotLeft (Fin 13) (Matrix.single 1 0 1),
    ⟨Matrix.single 0 1 1, rfl⟩, ⟨Matrix.single 1 0 1, rfl⟩, ?_⟩
  intro hcomm
  rw [← map_mul, ← map_mul, Matrix.single_mul_single_same,
    Matrix.single_mul_single_same] at hcomm
  have h00 := congrFun (congrFun hcomm ((0 : Fin 13), (0 : Fin 13)))
    ((0 : Fin 13), (0 : Fin 13))
  have hL : (slotLeft (Fin 13) (Matrix.single (0 : Fin 13) 0 (1 * 1 : ℂ)))
      ((0 : Fin 13), (0 : Fin 13)) ((0 : Fin 13), (0 : Fin 13)) = 1 := by
    show Matrix.single (0 : Fin 13) 0 (1 * 1 : ℂ) 0 0 *
      (1 : Matrix (Fin 13) (Fin 13) ℂ) 0 0 = 1
    simp [Matrix.single]
  have hR : (slotLeft (Fin 13) (Matrix.single (1 : Fin 13) 1 (1 * 1 : ℂ)))
      ((0 : Fin 13), (0 : Fin 13)) ((0 : Fin 13), (0 : Fin 13)) = 0 := by
    show Matrix.single (1 : Fin 13) 1 (1 * 1 : ℂ) 0 0 *
      (1 : Matrix (Fin 13) (Fin 13) ℂ) 0 0 = 0
    simp [Matrix.single]
  rw [hL, hR] at h00
  exact one_ne_zero h00

/-! ## The composed receipt -/

/-- **The record-layer net morphism receipt for the carrier join.**  The
antecedent bundle is the pair of committed premise instances
`designatedPairPremises` and `supportDisjointPairPremises`, whose fields
carry the registered premise rows PR-32 (Cartesian joint carrier), PR-33
(diamond region map), and PR-34 (slot split); the receipt consumes them
through the committed region maps.  Clauses: the record layers are
exactly the diagonal restrictions of the committed regional algebras on
both carriers; both record-layer nets and the induced net on the join
are isotone; both committed embeddings are region-by-region net
morphisms preserving inclusions and dimensions, with the hinge region of
dimension 13; the committed conditional expectations of both nets
restrict to the record layer on every region, landing in and fixing the
record layers; the matched left-slot pair descends to the join hinge
projection commuting with both embeddings and preserving the hinge
restriction; the committed join evolution carries every induced region
onto itself; the gluing submodule is not multiplicatively closed, so the
join carries no quotient algebra structure; and the committed left
regional algebra is noncommutative while the record layer is
commutative.  Physical time, causality, and spacetime attachment stay
open under PR-52; the source-selected time-slice content of PR-58 stays
open; the operator-layer net morphism stays open, delimited by the two
negative clauses. -/
theorem joinNetMorphism_receipt :
    (∀ (U : TwoSlotRegion) (v : CarrierFn88),
      Matrix.diagonal v ∈ designatedPairPremises.regionMap U ↔
        v ∈ recordLayer88 U) ∧
    (∀ (U : TwoSlotRegion) (w : CarrierFn247),
      Matrix.diagonal w ∈ supportDisjointPairPremises.regionMap U ↔
        w ∈ recordLayer247 U) ∧
    (∀ U V : TwoSlotRegion, TwoSlotRegion.le U V = true →
      recordLayer88 U ≤ recordLayer88 V) ∧
    (∀ U V : TwoSlotRegion, TwoSlotRegion.le U V = true →
      recordLayer247 U ≤ recordLayer247 V) ∧
    (∀ R S : JoinRegion, JoinRegion.le R S = true →
      joinRegionLayer R ≤ joinRegionLayer S) ∧
    (∀ U : TwoSlotRegion,
      (recordLayer88 U).map embed88 = joinRegionLayer (JoinRegion.side88 U)) ∧
    (∀ U : TwoSlotRegion,
      (recordLayer247 U).map embed247 =
        joinRegionLayer (JoinRegion.side247 U)) ∧
    (∀ U : TwoSlotRegion,
      Module.finrank ℂ (joinRegionLayer (JoinRegion.side88 U)) =
        recordDim88 U) ∧
    (∀ U : TwoSlotRegion,
      Module.finrank ℂ (joinRegionLayer (JoinRegion.side247 U)) =
        recordDim247 U) ∧
    Module.finrank ℂ (joinRegionLayer JoinRegion.hinge) = 13 ∧
    (∀ (U : TwoSlotRegion) (v : CarrierFn88),
      twoSlotNet.expect U (Matrix.diagonal v) =
        Matrix.diagonal (recordExpect88 U v)) ∧
    (∀ (U : TwoSlotRegion) (w : CarrierFn247),
      twoSlotNet247.expect U (Matrix.diagonal w) =
        Matrix.diagonal (recordExpect247 U w)) ∧
    (∀ (U : TwoSlotRegion) (v : CarrierFn88),
      recordExpect88 U v ∈ recordLayer88 U) ∧
    (∀ U : TwoSlotRegion, ∀ v ∈ recordLayer88 U, recordExpect88 U v = v) ∧
    (∀ (U : TwoSlotRegion) (w : CarrierFn247),
      recordExpect247 U w ∈ recordLayer247 U) ∧
    (∀ U : TwoSlotRegion, ∀ w ∈ recordLayer247 U, recordExpect247 U w = w) ∧
    (∀ v : CarrierFn88, joinHingeProject (embed88 v) =
      embed88 (recordExpect88 TwoSlotRegion.left v)) ∧
    (∀ w : CarrierFn247, joinHingeProject (embed247 w) =
      embed247 (recordExpect247 TwoSlotRegion.left w)) ∧
    (∀ z : JoinCarrier,
      joinHingeProject z ∈ joinRegionLayer JoinRegion.hinge) ∧
    (∀ z : JoinCarrier,
      hingeRestrict (joinHingeProject z) = hingeRestrict z) ∧
    (∀ (t : Fin 31) (R : JoinRegion),
      (joinRegionLayer R).map (joinStep t) = joinRegionLayer R) ∧
    (∃ x : CarrierFn88 × CarrierFn247,
      x ∈ hingeGraph ∧ x * x ∉ hingeGraph) ∧
    (∃ X Y : Matrix PairIndex247 PairIndex247 ℂ,
      X ∈ twoSlotAlgebra247 TwoSlotRegion.left ∧
      Y ∈ twoSlotAlgebra247 TwoSlotRegion.left ∧ X * Y ≠ Y * X) ∧
    (∀ v w : CarrierFn247,
      Matrix.diagonal v * Matrix.diagonal w =
        Matrix.diagonal w * Matrix.diagonal v) :=
  ⟨diagonal_mem_twoSlotAlgebra_iff, diagonal_mem_twoSlotAlgebra247_iff,
    fun _ _ => recordLayer88_isotone, fun _ _ => recordLayer247_isotone,
    joinRegionLayer_isotone,
    embed88_maps_recordLayer, embed247_maps_recordLayer,
    joinRegionLayer_side88_finrank, joinRegionLayer_side247_finrank,
    joinRegionLayer_hinge_finrank,
    twoSlotExpect_diagonal, twoSlotExpect247_diagonal,
    recordExpect88_mem, fun U _ hv => recordExpect88_fixes U hv,
    recordExpect247_mem, fun U _ hw => recordExpect247_fixes U hw,
    joinHingeProject_embed88, joinHingeProject_embed247,
    joinHingeProject_mem_hinge, hingeRestrict_joinHingeProject,
    joinStep_maps_joinRegionLayer,
    hingeGraph_not_mul_closed, twoSlotAlgebra247_left_noncommutative,
    recordDiagonal_commute⟩

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
-- Expected axioms per declaration: propext, Classical.choice, Quot.sound.
#print axioms OPH.QFT.recordLayer88_isotone
#print axioms OPH.QFT.recordLayer247_isotone
#print axioms OPH.QFT.diagonal_mem_bot_iff
#print axioms OPH.QFT.diagonal_mem_slotLeft_iff
#print axioms OPH.QFT.diagonal_mem_slotRight_iff
#print axioms OPH.QFT.diagonal_mem_twoSlotAlgebra_iff
#print axioms OPH.QFT.diagonal_mem_twoSlotAlgebra247_iff
#print axioms OPH.QFT.recordLayer88_finrank
#print axioms OPH.QFT.recordLayer247_finrank
#print axioms OPH.QFT.JoinRegion.le_refl'
#print axioms OPH.QFT.JoinRegion.le_trans'
#print axioms OPH.QFT.JoinRegion.le_antisymm'
#print axioms OPH.QFT.joinRegionLayer_isotone
#print axioms OPH.QFT.embed88_image_nested
#print axioms OPH.QFT.embed247_image_nested
#print axioms OPH.QFT.embed88_maps_recordLayer
#print axioms OPH.QFT.embed247_maps_recordLayer
#print axioms OPH.QFT.hingeLayer_eq_247
#print axioms OPH.QFT.joinRegionLayer_side88_finrank
#print axioms OPH.QFT.joinRegionLayer_side247_finrank
#print axioms OPH.QFT.joinRegionLayer_hinge_finrank
#print axioms OPH.QFT.ptraceSnd_diagonal
#print axioms OPH.QFT.ptraceFst_diagonal
#print axioms OPH.QFT.leftSlotExpectation_diagonal
#print axioms OPH.QFT.rightSlotExpectation_diagonal
#print axioms OPH.QFT.scalarExpectation_diagonal
#print axioms OPH.QFT.twoSlotExpect_diagonal
#print axioms OPH.QFT.twoSlotExpect247_diagonal
#print axioms OPH.QFT.recordExpect88_mem
#print axioms OPH.QFT.recordExpect247_mem
#print axioms OPH.QFT.recordExpect88_fixes
#print axioms OPH.QFT.recordExpect247_fixes
#print axioms OPH.QFT.joinHingeProject
#print axioms OPH.QFT.joinHingeProject_embed88
#print axioms OPH.QFT.joinHingeProject_embed247
#print axioms OPH.QFT.joinHingeProject_mem_hinge
#print axioms OPH.QFT.joinHingeProject_fixes_hinge
#print axioms OPH.QFT.joinHingeProject_idem
#print axioms OPH.QFT.hingeRestrict_joinHingeProject
#print axioms OPH.QFT.evolveFn88_map_recordLayer
#print axioms OPH.QFT.evolveFn247_map_recordLayer
#print axioms OPH.QFT.joinStep_maps_joinRegionLayer
#print axioms OPH.QFT.hingeGraph_not_mul_closed
#print axioms OPH.QFT.join_quotient_mul_ill_defined
#print axioms OPH.QFT.recordDiagonal_commute
#print axioms OPH.QFT.twoSlotAlgebra247_left_noncommutative
#print axioms OPH.QFT.joinNetMorphism_receipt
