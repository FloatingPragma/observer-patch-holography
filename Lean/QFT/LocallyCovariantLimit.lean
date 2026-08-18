import Mathlib.CategoryTheory.Functor.Basic
import QFT.ObserverAccessCut

/-!
# The E3 admissible event-region category and its covariant observable functor

This module implements the first deliverable of completion-plan issue
`#700`: the admissible event-region and causal-embedding category over
one finite causal observer net, the covariant functor from that category
into regulator-indexed observable algebras, and one inhabited nontrivial
example.

An *admissible event region* is a declared region of the E1 net at a
declared regulator: admissibility means exactly membership in the net's
own region system, so no new region type is introduced.  A *causal
embedding* between admissible event regions is a regulator refinement
together with a declared region inclusion of the refined source region
into the target region.  Composition is the E1 refinement-transitivity
receipt, so the category structure is a repackaging of fields the net
already carries, not a new assumption.

`observableFunctor` sends an admissible event region to its declared
local observable algebra and a causal embedding to the composite of the
A3 algebra-refinement map with the isotony inclusion.  Functoriality is
exactly the statement that refining along a composite embedding agrees
with refining in two steps; both identities reduce to the tower's
`algebra_refine_refl` and `algebra_refine_trans` receipts, so the
functor laws are the finite refinement-independence receipt in
categorical form.

The witness example `supportGradedNet` is a finite causal observer net
over the E6 witness tower whose regional algebras genuinely vary with
the region: the empty region carries the scalar algebra and every
nonempty region carries the full diagonal algebra.  Several committed
net inhabitants already assign region-dependent algebras — `glueNet`,
`characterCausalNet`, `properMeetNet`, `simEarnedNet`, and
`sourceRegionalNet` all do — so region dependence is not what is new
here; what this example adds is the nonvacuity receipt phrased through
the functor.  `observableFunctor_obj_strict` and
`observableFunctor_map_not_surjective` certify that the functor image of
the embedding `bottomToTop` is a proper inclusion, and
`supportGradedNet_has_disjoint_pair` keeps a genuinely declared disjoint
pair with commuting algebras in the example.

The E3 filtered-colimit layer is implemented at the end of this module.
`TowerColimit` is the filtered-colimit observable algebra of one
consensus tower over its directed regulator preorder: germs of private
observables, identified when a declared common refinement equalizes
them.  That identification is exactly the quotient by the
refinement-null kernel (`colimitMk_eq_zero_iff`), and no bonding map is
assumed injective anywhere; the collapsing tower witness shows the
non-injective regime is genuinely inhabited, and
`algebraRefine_injective_of_observer` shows the corpus confines it to
observer-free regulators.  The colimit carries the star `ℂ`-algebra
structure, the canonical cocone `colimitMkHom`, and the concrete
universal property (`colimitDesc`, `colimitHom_ext`).
`limitLocalAlgebra` is the limit local algebra of an admissible event
region, with limit-level isotony (`limitLocalAlgebra_mono`), locality
(`limitLocalAlgebra_commute`), and refinement independence
(`limitLocalAlgebra_refine`).  `IsCauchyEmbedding` is the Cauchy
embedding class, defined through the net's own contravariant
restriction maps, and `limitLocalAlgebra_eq_of_isCauchyEmbedding` is
the time-slice property at the limit.  `publicObservableFunctor` with
`publicRecordInclusion` is the public-record subfunctor of the
observable functor, and `limitPublicAlgebra` is its commutative limit
public algebra.

## Claim boundary

This module constructs the finite-level category and functor, and the
E3 filtered-colimit layer over them: the filtered-colimit observable
algebra with its universal property, limit-level isotony, locality and
refinement independence, the Cauchy embedding class with the limit
time-slice property, and the public-record subfunctor with its
naturality receipt and commutative limit public algebra.  The colimit
is algebraic: no C*-norm, limit norm, or completion is constructed.  No
limit-level contravariant restriction map is defined: the corpus
carries no naturality receipt relating `restrict` to `algebraRefine`,
and the germ-wise candidate map is well-defined only given such a
receipt, so the limit-level restriction is left open rather than
constructed.  The public-record subfunctor takes a
`CoherentObserverFamily` as an input; every committed constant-observer
tower carries one.  No CP, CPTP, or adaptive physical locality property
is inferred from B2 Kraus normalisation or from the E1 finite
interface.  The causal embeddings are declared order-theoretic data of
the E1 net; no Lorentzian or spacetime reading is attached.
-/

namespace OPH.QFT

open OPH.Tower
open CategoryTheory

universe u v

variable {ι : Type u} [Preorder ι] {T : ConsensusTower ι}

/-- An admissible event region: a declared region of the net at a
declared regulator.  Admissibility is membership in the net's own region
system; no new region type is introduced. -/
structure EventRegion (N : FiniteCausalObserverNet T) : Type u where
  /-- The regulator at which the region is declared. -/
  regulator : ι
  /-- The declared region. -/
  region : N.Region regulator

/-- A causal embedding of admissible event regions: a regulator
refinement together with a declared inclusion of the refined source
region into the target region.  Both fields are proofs, so parallel
causal embeddings are equal and the category below is thin. -/
structure CausalEmbedding {N : FiniteCausalObserverNet T}
    (U V : EventRegion N) : Type where
  /-- The regulator refinement. -/
  le : U.regulator ≤ V.regulator
  /-- The refined source region is declared below the target region. -/
  emb : N.regionLE V.regulator (N.regionRefine le U.region) V.region

instance {N : FiniteCausalObserverNet T} (U V : EventRegion N) :
    Subsingleton (CausalEmbedding U V) :=
  ⟨fun f g => by cases f; cases g; rfl⟩

/-- The admissible event-region and causal-embedding category.  Identity
and composition are the net's region-refinement reflexivity and
transitivity receipts; the associativity and unit laws hold by proof
irrelevance since the category is thin. -/
instance eventRegionCategory {N : FiniteCausalObserverNet T} :
    Category (EventRegion N) where
  Hom U V := CausalEmbedding U V
  id U :=
    ⟨le_refl U.regulator, by
      rw [N.region_refine_refl]
      exact N.regionLE_refl U.regulator U.region⟩
  comp {U V W} f g :=
    ⟨f.le.trans g.le, by
      rw [← N.region_refine_trans f.le g.le U.region]
      exact N.regionLE_trans W.regulator
        (N.region_refine_mono g.le f.emb) g.emb⟩
  id_comp _ := Subsingleton.elim _ _
  comp_id _ := Subsingleton.elim _ _
  assoc _ _ _ := Subsingleton.elim _ _

/-- A regulator-indexed observable algebra: a star subalgebra of the
tower's private matrix algebra at one regulator.  This is the target
object type of the covariant observable functor. -/
structure ObservableAlgebra (T : ConsensusTower ι) : Type u where
  /-- The regulator carrying the algebra. -/
  regulator : ι
  /-- The observable algebra at that regulator. -/
  algebra : StarSubalgebra ℂ (ConsensusTower.PrivateAlgebra T regulator)

/-- The category of regulator-indexed observable algebras with star
algebra homomorphisms as morphisms. -/
noncomputable instance observableAlgebraCategory :
    Category (ObservableAlgebra T) where
  Hom A B := A.algebra →⋆ₐ[ℂ] B.algebra
  id A := StarAlgHom.id ℂ A.algebra
  comp f g := g.comp f
  id_comp _ := rfl
  comp_id _ := rfl
  assoc _ _ _ := rfl

/-- The locally covariant observable functor at the finite level: an
admissible event region goes to its declared local observable algebra;
a causal embedding goes to the A3 algebra-refinement map followed by the
E1 isotony inclusion.  The functor laws are the tower's refinement
reflexivity and transitivity receipts, i.e. the finite
refinement-independence statement in categorical form. -/
noncomputable def observableFunctor (N : FiniteCausalObserverNet T) :
    EventRegion N ⥤ ObservableAlgebra T where
  obj U := ⟨U.regulator, N.localAlgebra U.regulator U.region⟩
  map {U V} f :=
    (N.localInclusion (CausalEmbedding.emb f)).comp
      (N.localRefine (CausalEmbedding.le f) U.region)
  map_id U := by
    refine StarAlgHom.ext fun X => Subtype.ext ?_
    exact T.algebra_refine_refl U.regulator X.1
  map_comp {U V W} f g := by
    refine StarAlgHom.ext fun X => Subtype.ext ?_
    exact (T.algebra_refine_trans (CausalEmbedding.le f)
      (CausalEmbedding.le g) X.1).symm

/-- The functor image of a causal embedding is definitionally the A3
algebra-refinement map followed by the E1 isotony inclusion; at the
matrix level the inclusion adds nothing. -/
theorem observableFunctor_map_eq {N : FiniteCausalObserverNet T}
    {U V : EventRegion N} (f : U ⟶ V) :
    (observableFunctor N).map f =
      (N.localInclusion (CausalEmbedding.emb f)).comp
        (N.localRefine (CausalEmbedding.le f) U.region) :=
  rfl

/-! ## A support-graded inhabitant with genuinely varying local algebras

Region-dependent local algebras are already committed on other net
inhabitants; what no committed net carries is a machine-checked proper
inclusion between two of its regional algebras, which is what the
observable functor needs to be certified nonconstant.  The
support-graded net below is the
finite nonvacuity model for E3: regions are subsets of a two-site label
set, and the regional algebra is the diagonal algebra of matrices whose
diagonal is constant outside the region.  The empty region carries
exactly the scalars and any nonempty region carries the full diagonal
algebra, so the observable functor takes a proper inclusion on the
embedding of the empty region into a nonempty one. -/

section SupportGradedExample

/-- The least site outside `U`, defaulting to site `0` when `U` is the
full site set.  It is the canonical representative at which a
compression reads off the common outside value. -/
def offSite (U : Finset (Fin 2)) : Fin 2 :=
  if 0 ∈ U then 1 else 0

/-- Any site outside `U` certifies that the canonical outside site is
outside `U` as well. -/
theorem offSite_not_mem :
    ∀ (U : Finset (Fin 2)) (i : Fin 2), i ∉ U → offSite U ∉ U := by
  decide

/-- Outside a nonempty subset of the two-site label set there is at most
one site, so any outside site is the canonical one. -/
theorem eq_offSite_of_nonempty :
    ∀ (U : Finset (Fin 2)) (i : Fin 2),
      U.Nonempty → i ∉ U → i = offSite U := by
  decide

/-- If the canonical site outside `W` escapes an ambient `V ⊇ W`, the
two canonical outside sites agree. -/
theorem offSite_eq_of_not_mem :
    ∀ (V W : Finset (Fin 2)),
      W ⊆ V → offSite W ∉ V → offSite V = offSite W := by
  decide

/-- The support-graded diagonal algebra: diagonal matrices whose
diagonal entries are constant outside `U`.  The empty region carries
exactly the scalar multiples of the identity; a nonempty region of the
two-site label set carries the full diagonal algebra. -/
noncomputable def supportDiagonalAlgebra (U : Finset (Fin 2)) :
    StarSubalgebra ℂ (Matrix (Fin 2) (Fin 2) ℂ) where
  carrier := {X | ∃ f : Fin 2 → ℂ,
    X = Matrix.diagonal f ∧ ∀ i ∉ U, ∀ j ∉ U, f i = f j}
  mul_mem' := by
    rintro X Y ⟨f, rfl, hf⟩ ⟨g, rfl, hg⟩
    refine ⟨f * g, Matrix.diagonal_mul_diagonal f g, fun i hi j hj => ?_⟩
    simp only [Pi.mul_apply]
    rw [hf i hi j hj, hg i hi j hj]
  one_mem' := ⟨fun _ => 1, Matrix.diagonal_one.symm, fun _ _ _ _ => rfl⟩
  add_mem' := by
    rintro X Y ⟨f, rfl, hf⟩ ⟨g, rfl, hg⟩
    refine ⟨f + g, Matrix.diagonal_add f g, fun i hi j hj => ?_⟩
    simp only [Pi.add_apply]
    rw [hf i hi j hj, hg i hi j hj]
  zero_mem' := ⟨0, Matrix.diagonal_zero.symm, fun _ _ _ _ => rfl⟩
  algebraMap_mem' := fun c => by
    refine ⟨fun _ => c, ?_, fun _ _ _ _ => rfl⟩
    ext i j
    simp [Matrix.algebraMap_eq_diagonal, Matrix.diagonal_apply]
  star_mem' := by
    rintro X ⟨f, rfl, hf⟩
    refine ⟨star f, Matrix.diagonal_conjTranspose f, fun i hi j hj => ?_⟩
    simp only [Pi.star_apply]
    rw [hf i hi j hj]

theorem mem_supportDiagonalAlgebra {U : Finset (Fin 2)}
    {X : Matrix (Fin 2) (Fin 2) ℂ} :
    X ∈ supportDiagonalAlgebra U ↔ ∃ f : Fin 2 → ℂ,
      X = Matrix.diagonal f ∧ ∀ i ∉ U, ∀ j ∉ U, f i = f j :=
  Iff.rfl

/-- Isotony of the support-graded assignment: enlarging the region
weakens the outside-constancy constraint. -/
theorem supportDiagonalAlgebra_mono {U V : Finset (Fin 2)}
    (hUV : U ⊆ V) :
    supportDiagonalAlgebra U ≤ supportDiagonalAlgebra V := by
  rintro X ⟨f, rfl, hf⟩
  exact ⟨f, rfl, fun i hi j hj =>
    hf i (fun h => hi (hUV h)) j (fun h => hj (hUV h))⟩

/-- The support compression at `V`: keep the diagonal on `V`, replace
every outside diagonal entry by the entry at the canonical outside site,
and kill all off-diagonal entries. -/
noncomputable def supportCompress (V : Finset (Fin 2))
    (X : Matrix (Fin 2) (Fin 2) ℂ) : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.diagonal fun i =>
    if i ∈ V then X i i else X (offSite V) (offSite V)

theorem supportCompress_mem (V : Finset (Fin 2))
    (X : Matrix (Fin 2) (Fin 2) ℂ) :
    supportCompress V X ∈ supportDiagonalAlgebra V :=
  ⟨_, rfl, fun i hi j hj => by rw [if_neg hi, if_neg hj]⟩

/-- Compression at an ambient region fixes every member of a smaller
support-graded algebra. -/
theorem supportCompress_eq_self {U V : Finset (Fin 2)}
    {X : Matrix (Fin 2) (Fin 2) ℂ}
    (hX : X ∈ supportDiagonalAlgebra U) (hUV : U ⊆ V) :
    supportCompress V X = X := by
  obtain ⟨f, rfl, hf⟩ := hX
  ext i j
  by_cases hij : i = j
  · subst hij
    rw [supportCompress, Matrix.diagonal_apply_eq]
    by_cases hiV : i ∈ V
    · rw [if_pos hiV]
    · rw [if_neg hiV, Matrix.diagonal_apply_eq, Matrix.diagonal_apply_eq]
      exact hf (offSite V)
        (fun h => offSite_not_mem V i hiV (hUV h)) i (fun h => hiV (hUV h))
  · rw [supportCompress, Matrix.diagonal_apply_ne _ hij,
      Matrix.diagonal_apply_ne _ hij]

/-- On the two-site label set, compression at any nonempty region fixes
every diagonal matrix: the outside of a nonempty region has at most one
site, which is the canonical one. -/
theorem supportCompress_eq_self_of_nonempty {U : Finset (Fin 2)}
    (hU : U.Nonempty) {X : Matrix (Fin 2) (Fin 2) ℂ}
    (hX : ∃ f : Fin 2 → ℂ, X = Matrix.diagonal f) :
    supportCompress U X = X := by
  obtain ⟨f, rfl⟩ := hX
  ext i j
  by_cases hij : i = j
  · subst hij
    rw [supportCompress, Matrix.diagonal_apply_eq]
    by_cases hi : i ∈ U
    · rw [if_pos hi]
    · rw [if_neg hi, ← eq_offSite_of_nonempty U i hU hi]
  · rw [supportCompress, Matrix.diagonal_apply_ne _ hij,
      Matrix.diagonal_apply_ne _ hij]

/-- Two nested compressions collapse to the inner one. -/
theorem supportCompress_compress {V W : Finset (Fin 2)} (hWV : W ⊆ V)
    (X : Matrix (Fin 2) (Fin 2) ℂ) :
    supportCompress W (supportCompress V X) = supportCompress W X := by
  unfold supportCompress
  congr 1
  funext i
  by_cases hiW : i ∈ W
  · simp [hiW, hWV hiW, Matrix.diagonal_apply_eq]
  · by_cases how : offSite W ∈ V
    · simp [hiW, how, Matrix.diagonal_apply_eq]
    · simp [hiW, how, offSite_eq_of_not_mem V W hWV how,
        Matrix.diagonal_apply_eq]

theorem supportCompress_add (V : Finset (Fin 2))
    (X Y : Matrix (Fin 2) (Fin 2) ℂ) :
    supportCompress V (X + Y) =
      supportCompress V X + supportCompress V Y := by
  unfold supportCompress
  rw [Matrix.diagonal_add]
  congr 1
  funext i
  by_cases hi : i ∈ V <;> simp [hi]

theorem supportCompress_smul (V : Finset (Fin 2)) (c : ℂ)
    (X : Matrix (Fin 2) (Fin 2) ℂ) :
    supportCompress V (c • X) = c • supportCompress V X := by
  unfold supportCompress
  rw [← Matrix.diagonal_smul]
  congr 1
  funext i
  by_cases hi : i ∈ V <;> simp [hi, Matrix.smul_apply]

theorem supportCompress_mul (V : Finset (Fin 2))
    {X Y : Matrix (Fin 2) (Fin 2) ℂ}
    (hX : ∃ f : Fin 2 → ℂ, X = Matrix.diagonal f)
    (hY : ∃ g : Fin 2 → ℂ, Y = Matrix.diagonal g) :
    supportCompress V (X * Y) =
      supportCompress V X * supportCompress V Y := by
  obtain ⟨f, rfl⟩ := hX
  obtain ⟨g, rfl⟩ := hY
  unfold supportCompress
  rw [Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  congr 1
  funext i
  by_cases hi : i ∈ V <;> simp [hi]

theorem supportCompress_star (V : Finset (Fin 2))
    (X : Matrix (Fin 2) (Fin 2) ℂ) :
    supportCompress V (star X) = star (supportCompress V X) := by
  ext i j
  by_cases hij : i = j
  · subst hij
    by_cases hi : i ∈ V <;>
      simp [supportCompress, Matrix.star_apply, hi]
  · simp [supportCompress, Matrix.star_apply,
      Matrix.diagonal_apply_ne _ hij, Matrix.diagonal_apply_ne _ (Ne.symm hij)]

/-- The support restriction: compression bundled as a star algebra
homomorphism from the ambient support-graded algebra to a declared
subregion's algebra.  It retracts the isotony inclusion. -/
noncomputable def supportRestrict {U V : Finset (Fin 2)} (_hVU : V ⊆ U) :
    supportDiagonalAlgebra U →⋆ₐ[ℂ] supportDiagonalAlgebra V where
  toFun X := ⟨supportCompress V X.1, supportCompress_mem V X.1⟩
  map_one' := Subtype.ext
    (supportCompress_eq_self (one_mem (supportDiagonalAlgebra V))
      (subset_refl V))
  map_mul' X Y := by
    have hX : ∃ f : Fin 2 → ℂ, X.1 = Matrix.diagonal f := by
      obtain ⟨f, hf, -⟩ := X.2; exact ⟨f, hf⟩
    have hY : ∃ g : Fin 2 → ℂ, Y.1 = Matrix.diagonal g := by
      obtain ⟨g, hg, -⟩ := Y.2; exact ⟨g, hg⟩
    exact Subtype.ext (supportCompress_mul V hX hY)
  map_zero' := Subtype.ext
    (supportCompress_eq_self (zero_mem (supportDiagonalAlgebra V))
      (subset_refl V))
  map_add' X Y := Subtype.ext (supportCompress_add V X.1 Y.1)
  commutes' c := Subtype.ext
    (supportCompress_eq_self
      ((supportDiagonalAlgebra V).algebraMap_mem c) (subset_refl V))
  map_star' X := Subtype.ext (supportCompress_star V X.1)

/-- The support repair: compression bundled as a complex-linear
idempotent on the full private matrix algebra. -/
noncomputable def supportRepairLinearMap (U : Finset (Fin 2)) :
    Matrix (Fin 2) (Fin 2) ℂ →ₗ[ℂ] Matrix (Fin 2) (Fin 2) ℂ where
  toFun := supportCompress U
  map_add' X Y := supportCompress_add U X Y
  map_smul' c X := supportCompress_smul U c X

/-- Diagonal matrices commute; this is the locality receipt of the
support-graded net, proved in the plain matrix language so that no
instance mismatch with the tower's private-algebra spelling occurs. -/
theorem diagonal_commute (f g : Fin 2 → ℂ) :
    Matrix.diagonal f * Matrix.diagonal g =
      Matrix.diagonal g * Matrix.diagonal f := by
  rw [Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  congr 1
  funext i
  exact mul_comm (f i) (g i)

/-- The support-graded finite causal observer net over the E6 witness
tower: the first committed net inhabitant whose regional algebras
genuinely vary with the region. -/
noncomputable def supportGradedNet : FiniteCausalObserverNet witnessTower where
  Region := fun _ => Finset (Fin 2)
  regionFintype := fun _ => by classical exact inferInstance
  regionNonempty := fun _ => ⟨∅⟩
  regionLE := fun _ U V => U ⊆ V
  regionLE_refl := fun _ _ => Finset.Subset.rfl
  regionLE_trans := by
    intro r U V W hUV hVW
    exact fun x hx => hVW (hUV hx)
  regionLE_antisymm := by
    intro r U V hUV hVU
    exact hUV.antisymm hVU
  overlap := fun _ U V => U ∩ V
  overlap_le_left := fun _ _ _ => Finset.inter_subset_left
  overlap_le_right := fun _ _ _ => Finset.inter_subset_right
  le_overlap := fun _ _ _ _ hWU hWV => Finset.subset_inter hWU hWV
  disjoint := fun _ U V => U.Nonempty ∧ V.Nonempty ∧ Disjoint U V
  disjoint_symm := by
    intro r U V h
    exact ⟨h.2.1, h.1, h.2.2.symm⟩
  disjoint_irrefl := by
    intro r U h
    obtain ⟨x, hx⟩ := h.1
    exact (Finset.disjoint_left.mp h.2.2 hx) hx
  localAlgebra := fun _ U => supportDiagonalAlgebra U
  isotony := fun _ {_ _} hUV => supportDiagonalAlgebra_mono hUV
  locality := by
    intro r U V _ X Y
    obtain ⟨f, hf, -⟩ := X.2
    obtain ⟨g, hg, -⟩ := Y.2
    show X.1 * Y.1 = Y.1 * X.1
    rw [hf, hg]
    exact diagonal_commute f g
  restrict := fun _ {_ _} hVU => supportRestrict hVU
  restrict_refl := fun _ U X =>
    Subtype.ext (supportCompress_eq_self X.2 (subset_refl U))
  restrict_trans := fun _ {_ _ _} _ hWV X =>
    Subtype.ext (supportCompress_compress hWV X.1)
  restrict_inclusion := fun _ {U V} _ X =>
    Subtype.ext (supportCompress_eq_self X.2 (subset_refl U))
  regionRefine := fun _ U => U
  region_refine_refl := by intros; rfl
  region_refine_trans := by intros; rfl
  region_refine_mono := by intros; assumption
  overlap_natural := by intros; rfl
  disjoint_natural := by intros; assumption
  localAlgebra_natural := by intros; assumption
  repair := fun _ U => supportRepairLinearMap U
  repair_idempotent := fun _ U X =>
    supportCompress_compress (subset_refl U) X
  repair_fixes_region := fun _ U X =>
    supportCompress_eq_self X.2 (subset_refl U)
  repair_fixes_disjoint := by
    intro r U V hUV X
    obtain ⟨f, hf, -⟩ := X.2
    exact supportCompress_eq_self_of_nonempty hUV.1 ⟨f, hf⟩
  repair_natural := by intros; rfl

/-- The support-graded net keeps a genuinely declared pair of distinct,
nonempty disjoint region labels. -/
theorem supportGradedNet_has_disjoint_pair :
    supportGradedNet.disjoint () ({0} : Finset (Fin 2))
      ({1} : Finset (Fin 2)) := by
  simp [supportGradedNet]

/-- The first diagonal projector, the separating observable of the
example. -/
noncomputable def e0 : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.diagonal fun i => if i = 0 then 1 else 0

theorem e0_mem_singleton :
    e0 ∈ supportDiagonalAlgebra ({0} : Finset (Fin 2)) :=
  ⟨_, rfl, fun i hi j hj => by
    fin_cases i <;> fin_cases j <;> simp_all⟩

theorem e0_not_mem_empty :
    e0 ∉ supportDiagonalAlgebra (∅ : Finset (Fin 2)) := by
  rintro ⟨f, hf, hc⟩
  have h0 : (1 : ℂ) = f 0 := by
    have := congrArg (fun M => M 0 0) hf
    simpa [e0, Matrix.diagonal_apply_eq] using this
  have h1 : (0 : ℂ) = f 1 := by
    have := congrArg (fun M => M 1 1) hf
    simpa [e0, Matrix.diagonal_apply_eq] using this
  have : (1 : ℂ) = 0 := by
    rw [h0, h1]
    exact hc 0 (by simp) 1 (by simp)
  exact one_ne_zero this

/-- Strict interposition of the support grading: the empty region's
algebra is properly below a nonempty region's algebra. -/
theorem supportDiagonalAlgebra_bot_lt :
    supportDiagonalAlgebra (∅ : Finset (Fin 2)) <
      supportDiagonalAlgebra ({0} : Finset (Fin 2)) := by
  refine lt_of_le_of_ne
    (supportDiagonalAlgebra_mono (Finset.empty_subset _)) (fun h => ?_)
  have := e0_mem_singleton
  rw [← h] at this
  exact e0_not_mem_empty this

/-- The bottom admissible event region of the example: the empty region
at the sole regulator. -/
noncomputable def bottomRegion : EventRegion supportGradedNet :=
  ⟨(), (∅ : Finset (Fin 2))⟩

/-- The top admissible event region of the example: the first singleton
region at the sole regulator. -/
noncomputable def topRegion : EventRegion supportGradedNet :=
  ⟨(), ({0} : Finset (Fin 2))⟩

/-- The inhabiting causal embedding of the example. -/
noncomputable def bottomToTop : bottomRegion ⟶ topRegion :=
  ⟨le_refl (), Finset.empty_subset _⟩

/-- Nontriviality of the example, object side: the observable functor
sends the bottom region to a proper subalgebra of the top region's
image. -/
theorem observableFunctor_obj_strict :
    ((observableFunctor supportGradedNet).obj bottomRegion).algebra <
      ((observableFunctor supportGradedNet).obj topRegion).algebra :=
  supportDiagonalAlgebra_bot_lt

/-- The functor image of the inhabiting causal embedding, read back as
a star algebra homomorphism between the two support-graded algebras. -/
noncomputable def bottomToTopHom :
    supportDiagonalAlgebra (∅ : Finset (Fin 2)) →⋆ₐ[ℂ]
      supportDiagonalAlgebra ({0} : Finset (Fin 2)) :=
  (observableFunctor supportGradedNet).map bottomToTop

/-- Nontriviality of the example, morphism side: the functor image of
the inhabiting causal embedding is not surjective, since the first
diagonal projector has no preimage in the scalar algebra. -/
theorem observableFunctor_map_not_surjective :
    ¬ Function.Surjective bottomToTopHom := by
  intro hs
  obtain ⟨X, hX⟩ := hs ⟨e0, e0_mem_singleton⟩
  have hval : e0 = X.1 := by
    have := congrArg Subtype.val hX
    exact this.symm
  exact e0_not_mem_empty (hval ▸ X.2)

end SupportGradedExample

/- The quotient algebra below sits inside the full Mathlib star-algebra
hierarchy, so instance searches that end at the colimit's own instances
first wander the C*-algebra branches; the default synthesis budget is too
small for that detour. -/
set_option synthInstance.maxHeartbeats 400000

/-! ## The filtered-colimit observable algebra

The colimit is taken over the tower's regulator preorder, which is filtered
by `indexNonempty` and `commonRefinement`.  The construction never assumes
that a bonding map `algebraRefine` is injective: a germ is identified with
zero exactly when some declared refinement kills it, which is the quotient
by the refinement-null kernel (`colimitMk_eq_zero_iff`). -/

section TowerColimitCore

/-- Declared common refinements extend to any three regulators. -/
theorem exists_common_refinement₃ (T : ConsensusTower ι) (a b c : ι) :
    ∃ t : ι, a ≤ t ∧ b ≤ t ∧ c ≤ t := by
  obtain ⟨u, hau, hbu⟩ := T.commonRefinement a b
  obtain ⟨t, hut, hct⟩ := T.commonRefinement u c
  exact ⟨t, hau.trans hut, hbu.trans hut, hct⟩

/-- Declared common refinements extend to any four regulators. -/
theorem exists_common_refinement₄ (T : ConsensusTower ι) (a b c d : ι) :
    ∃ t : ι, a ≤ t ∧ b ≤ t ∧ c ≤ t ∧ d ≤ t := by
  obtain ⟨u, hau, hbu, hcu⟩ := exists_common_refinement₃ T a b c
  obtain ⟨t, hut, hdt⟩ := T.commonRefinement u d
  exact ⟨t, hau.trans hut, hbu.trans hut, hcu.trans hut, hdt⟩

/-- Equality after refinement persists under every further refinement; this
is `algebra_refine_trans` in the form the colimit uses. -/
theorem algebraRefine_eq_of_eq {r s a : ι} {hr : r ≤ a} {hs : s ≤ a}
    {X : ConsensusTower.PrivateAlgebra T r}
    {Y : ConsensusTower.PrivateAlgebra T s}
    (e : T.algebraRefine hr X = T.algebraRefine hs Y) {t : ι} (hat : a ≤ t) :
    T.algebraRefine (hr.trans hat) X = T.algebraRefine (hs.trans hat) Y := by
  rw [← T.algebra_refine_trans hr hat X, ← T.algebra_refine_trans hs hat Y, e]

/-- A germ representative for the filtered colimit: a regulator together
with a private observable declared at that regulator. -/
def TowerGerm (T : ConsensusTower ι) : Type u :=
  Σ r : ι, ConsensusTower.PrivateAlgebra T r

/-- Two germ representatives are identified when some declared common
refinement maps them to the same observable.  Injectivity of the bonding
maps is not assumed anywhere: a representative that a later bonding map
collapses is identified with the collapsed value. -/
def germRel (T : ConsensusTower ι) (p q : TowerGerm T) : Prop :=
  ∃ (t : ι) (hp : p.1 ≤ t) (hq : q.1 ≤ t),
    T.algebraRefine hp p.2 = T.algebraRefine hq q.2

theorem germRel_refl (p : TowerGerm T) : germRel T p p :=
  ⟨p.1, le_rfl, le_rfl, rfl⟩

theorem germRel_symm {p q : TowerGerm T} (h : germRel T p q) : germRel T q p := by
  obtain ⟨t, hp, hq, e⟩ := h
  exact ⟨t, hq, hp, e.symm⟩

theorem germRel_trans {p q r : TowerGerm T}
    (h₁ : germRel T p q) (h₂ : germRel T q r) : germRel T p r := by
  obtain ⟨a, hpa, hqa, ea⟩ := h₁
  obtain ⟨b, hqb, hrb, eb⟩ := h₂
  obtain ⟨t, hat, hbt⟩ := T.commonRefinement a b
  refine ⟨t, hpa.trans hat, hrb.trans hbt, ?_⟩
  calc T.algebraRefine (hpa.trans hat) p.2
      = T.algebraRefine (hqa.trans hat) q.2 := algebraRefine_eq_of_eq ea hat
    _ = T.algebraRefine (hqb.trans hbt) q.2 := rfl
    _ = T.algebraRefine (hrb.trans hbt) r.2 := algebraRefine_eq_of_eq eb hbt

/-- The germ setoid of one consensus tower. -/
instance germSetoid (T : ConsensusTower ι) : Setoid (TowerGerm T) :=
  ⟨germRel T, ⟨germRel_refl, germRel_symm, germRel_trans⟩⟩

/-- The carrier of the filtered-colimit observable algebra: germs of private
observables under the refinement-generated identification. -/
def TowerColimit (T : ConsensusTower ι) : Type u :=
  Quotient (germSetoid T)

/-- The class of a private observable in the filtered-colimit algebra. -/
def colimitMk {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    TowerColimit T :=
  Quotient.mk (germSetoid T) ⟨r, X⟩

/-- Refinement independence at the element level: reading an observable at
any declared refinement does not change its colimit class.  This is the
cocone receipt of the colimit and is `algebra_refine_trans` in quotient
form. -/
theorem colimitMk_refine {r s : ι} (hrs : r ≤ s)
    (X : ConsensusTower.PrivateAlgebra T r) :
    colimitMk (T.algebraRefine hrs X) = colimitMk X :=
  Quotient.sound ⟨s, le_rfl, hrs, by rw [T.algebra_refine_refl]⟩

/-- Two private observables share a colimit class exactly when some declared
common refinement equalizes them. -/
theorem colimitMk_eq_colimitMk_iff {r s : ι}
    {X : ConsensusTower.PrivateAlgebra T r}
    {Y : ConsensusTower.PrivateAlgebra T s} :
    colimitMk (T := T) X = colimitMk Y ↔
      ∃ (t : ι) (hr : r ≤ t) (hs : s ≤ t),
        T.algebraRefine hr X = T.algebraRefine hs Y :=
  ⟨fun h => Quotient.exact h,
    fun ⟨t, hr, hs, e⟩ => Quotient.sound ⟨t, hr, hs, e⟩⟩

/-- The colimit is jointly exhausted by the net's own data: every element is
the class of some finite-level observable. -/
theorem colimitMk_surjective (x : TowerColimit T) :
    ∃ (r : ι) (X : ConsensusTower.PrivateAlgebra T r), x = colimitMk X :=
  Quotient.inductionOn x fun p => ⟨p.1, p.2, rfl⟩

/-- Representatives can be pushed above any given regulator. -/
theorem exists_rep_above (x : TowerColimit T) (r₀ : ι) :
    ∃ (t : ι) (_ : r₀ ≤ t) (X : ConsensusTower.PrivateAlgebra T t),
      x = colimitMk X := by
  obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
  obtain ⟨t, hrt, hr₀t⟩ := T.commonRefinement r r₀
  exact ⟨t, hr₀t, T.algebraRefine hrt X, (colimitMk_refine hrt X).symm⟩

/-- Any two colimit elements admit representatives at one common regulator. -/
theorem exists_rep₂ (x y : TowerColimit T) :
    ∃ (r : ι) (X Y : ConsensusTower.PrivateAlgebra T r),
      x = colimitMk X ∧ y = colimitMk Y := by
  obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
  obtain ⟨t, hrt, Y, rfl⟩ := exists_rep_above y r
  exact ⟨t, T.algebraRefine hrt X, Y, (colimitMk_refine hrt X).symm, rfl⟩

/-- Any three colimit elements admit representatives at one common
regulator. -/
theorem exists_rep₃ (x y z : TowerColimit T) :
    ∃ (r : ι) (X Y Z : ConsensusTower.PrivateAlgebra T r),
      x = colimitMk X ∧ y = colimitMk Y ∧ z = colimitMk Z := by
  obtain ⟨r, X, Y, rfl, rfl⟩ := exists_rep₂ x y
  obtain ⟨t, hrt, Z, rfl⟩ := exists_rep_above z r
  exact ⟨t, T.algebraRefine hrt X, T.algebraRefine hrt Y, Z,
    (colimitMk_refine hrt X).symm, (colimitMk_refine hrt Y).symm, rfl⟩

/-- A declared common refinement of the regulators of two germ
representatives. -/
noncomputable def germUpper (p q : TowerGerm T) : ι :=
  (T.commonRefinement p.1 q.1).choose

theorem germUpper_left (p q : TowerGerm T) : p.1 ≤ germUpper p q :=
  (T.commonRefinement p.1 q.1).choose_spec.1

theorem germUpper_right (p q : TowerGerm T) : q.1 ≤ germUpper p q :=
  (T.commonRefinement p.1 q.1).choose_spec.2

/-- Combine two germ representatives at a declared common refinement with a
regulator-wise binary operation. -/
noncomputable def germLift₂
    (F : ∀ r : ι, ConsensusTower.PrivateAlgebra T r →
      ConsensusTower.PrivateAlgebra T r → ConsensusTower.PrivateAlgebra T r)
    (p q : TowerGerm T) : TowerGerm T :=
  ⟨germUpper p q,
    F (germUpper p q) (T.algebraRefine (germUpper_left p q) p.2)
      (T.algebraRefine (germUpper_right p q) q.2)⟩

/-- Any regulator-wise binary operation that commutes with the bonding maps
descends to germ classes. -/
theorem germLift₂_sound
    (F : ∀ r : ι, ConsensusTower.PrivateAlgebra T r →
      ConsensusTower.PrivateAlgebra T r → ConsensusTower.PrivateAlgebra T r)
    (hF : ∀ {r s : ι} (hrs : r ≤ s)
      (X Y : ConsensusTower.PrivateAlgebra T r),
      T.algebraRefine hrs (F r X Y) =
        F s (T.algebraRefine hrs X) (T.algebraRefine hrs Y))
    {p p' q q' : TowerGerm T} (hp : germRel T p p') (hq : germRel T q q') :
    germRel T (germLift₂ F p q) (germLift₂ F p' q') := by
  obtain ⟨a, hpa, hp'a, ea⟩ := hp
  obtain ⟨b, hqb, hq'b, eb⟩ := hq
  obtain ⟨t, h1, h2, h3, h4⟩ :=
    exists_common_refinement₄ T (germUpper p q) (germUpper p' q') a b
  refine ⟨t, h1, h2, ?_⟩
  show T.algebraRefine h1
      (F (germUpper p q) (T.algebraRefine (germUpper_left p q) p.2)
        (T.algebraRefine (germUpper_right p q) q.2)) =
    T.algebraRefine h2
      (F (germUpper p' q') (T.algebraRefine (germUpper_left p' q') p'.2)
        (T.algebraRefine (germUpper_right p' q') q'.2))
  rw [hF h1, hF h2,
    T.algebra_refine_trans (germUpper_left p q) h1,
    T.algebra_refine_trans (germUpper_right p q) h1,
    T.algebra_refine_trans (germUpper_left p' q') h2,
    T.algebra_refine_trans (germUpper_right p' q') h2,
    algebraRefine_eq_of_eq ea h3, algebraRefine_eq_of_eq eb h4]

/-- Apply a regulator-wise unary operation to a germ representative. -/
def germLift₁
    (F : ∀ r : ι, ConsensusTower.PrivateAlgebra T r →
      ConsensusTower.PrivateAlgebra T r)
    (p : TowerGerm T) : TowerGerm T :=
  ⟨p.1, F p.1 p.2⟩

/-- Any regulator-wise unary operation that commutes with the bonding maps
descends to germ classes. -/
theorem germLift₁_sound
    (F : ∀ r : ι, ConsensusTower.PrivateAlgebra T r →
      ConsensusTower.PrivateAlgebra T r)
    (hF : ∀ {r s : ι} (hrs : r ≤ s) (X : ConsensusTower.PrivateAlgebra T r),
      T.algebraRefine hrs (F r X) = F s (T.algebraRefine hrs X))
    {p q : TowerGerm T} (h : germRel T p q) :
    germRel T (germLift₁ F p) (germLift₁ F q) := by
  obtain ⟨t, hp, hq, e⟩ := h
  refine ⟨t, hp, hq, ?_⟩
  show T.algebraRefine hp (F p.1 p.2) = T.algebraRefine hq (F q.1 q.2)
  rw [hF hp, hF hq, e]

/-- The refinement-compatible operations descend to the quotient; addition
and multiplication push both representatives to a declared common
refinement first.  The star involution and the scalar action are kept as
raw functions rather than `Star`/`SMul` instances: those two classes occur
as *arguments* of `StarModule` and the star-subalgebra machinery, so their
only registered instance paths must be the canonical bundles below
(`StarRing`, `Algebra`), or instance unification fails on a
non-definitional diamond. -/
noncomputable instance : Add (TowerColimit T) :=
  ⟨Quotient.map₂ (germLift₂ fun _ X Y => X + Y)
    (fun _ _ hp _ _ hq =>
      germLift₂_sound _ (fun h X Y => map_add (T.algebraRefine h) X Y) hp hq)⟩

noncomputable instance : Mul (TowerColimit T) :=
  ⟨Quotient.map₂ (germLift₂ fun _ X Y => X * Y)
    (fun _ _ hp _ _ hq =>
      germLift₂_sound _ (fun h X Y => map_mul (T.algebraRefine h) X Y) hp hq)⟩

noncomputable instance : Neg (TowerColimit T) :=
  ⟨Quotient.map (germLift₁ fun _ X => -X)
    (fun _ _ h =>
      germLift₁_sound _ (fun hrs X => map_neg (T.algebraRefine hrs) X) h)⟩

/-- The star involution of colimit classes, regulator by regulator; wired
into the `StarRing` instance below rather than a bare `Star` instance. -/
noncomputable def colimitStar : TowerColimit T → TowerColimit T :=
  Quotient.map (germLift₁ fun _ X => star X)
    (fun _ _ h =>
      germLift₁_sound _ (fun hrs X => map_star (T.algebraRefine hrs) X) h)

/-- The complex scalar action on colimit classes, regulator by regulator;
wired into the `Algebra` instance below rather than a bare `SMul`
instance. -/
noncomputable def colimitSMul (c : ℂ) : TowerColimit T → TowerColimit T :=
  Quotient.map (germLift₁ fun _ X => c • X)
    (fun _ _ h =>
      germLift₁_sound _ (fun hrs X => map_smul (T.algebraRefine hrs) c X) h)

/-- The anchor regulator used to inject constants into the colimit. -/
noncomputable def towerBase (T : ConsensusTower ι) : ι :=
  T.indexNonempty.some

noncomputable instance : Zero (TowerColimit T) :=
  ⟨colimitMk (0 : ConsensusTower.PrivateAlgebra T (towerBase T))⟩

noncomputable instance : One (TowerColimit T) :=
  ⟨colimitMk (1 : ConsensusTower.PrivateAlgebra T (towerBase T))⟩

theorem colimitMk_add {r : ι} (X Y : ConsensusTower.PrivateAlgebra T r) :
    colimitMk X + colimitMk Y = colimitMk (X + Y) := by
  refine Quotient.sound ⟨germUpper ⟨r, X⟩ ⟨r, Y⟩, le_rfl,
    germUpper_left ⟨r, X⟩ ⟨r, Y⟩, ?_⟩
  show T.algebraRefine le_rfl
      (T.algebraRefine (germUpper_left ⟨r, X⟩ ⟨r, Y⟩) X +
        T.algebraRefine (germUpper_right ⟨r, X⟩ ⟨r, Y⟩) Y) =
    T.algebraRefine (germUpper_left ⟨r, X⟩ ⟨r, Y⟩) (X + Y)
  rw [T.algebra_refine_refl, map_add]

theorem colimitMk_mul {r : ι} (X Y : ConsensusTower.PrivateAlgebra T r) :
    colimitMk X * colimitMk Y = colimitMk (X * Y) := by
  refine Quotient.sound ⟨germUpper ⟨r, X⟩ ⟨r, Y⟩, le_rfl,
    germUpper_left ⟨r, X⟩ ⟨r, Y⟩, ?_⟩
  show T.algebraRefine le_rfl
      (T.algebraRefine (germUpper_left ⟨r, X⟩ ⟨r, Y⟩) X *
        T.algebraRefine (germUpper_right ⟨r, X⟩ ⟨r, Y⟩) Y) =
    T.algebraRefine (germUpper_left ⟨r, X⟩ ⟨r, Y⟩) (X * Y)
  rw [T.algebra_refine_refl, map_mul]

theorem colimitMk_neg {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    -colimitMk X = colimitMk (-X) :=
  rfl

theorem colimitMk_zero (r : ι) :
    colimitMk (0 : ConsensusTower.PrivateAlgebra T r) = (0 : TowerColimit T) := by
  obtain ⟨t, h1, h2⟩ := T.commonRefinement r (towerBase T)
  refine Quotient.sound ⟨t, h1, h2, ?_⟩
  show T.algebraRefine h1 (0 : ConsensusTower.PrivateAlgebra T r) =
    T.algebraRefine h2 (0 : ConsensusTower.PrivateAlgebra T (towerBase T))
  rw [map_zero, map_zero]

theorem colimitMk_one (r : ι) :
    colimitMk (1 : ConsensusTower.PrivateAlgebra T r) = (1 : TowerColimit T) := by
  obtain ⟨t, h1, h2⟩ := T.commonRefinement r (towerBase T)
  refine Quotient.sound ⟨t, h1, h2, ?_⟩
  show T.algebraRefine h1 (1 : ConsensusTower.PrivateAlgebra T r) =
    T.algebraRefine h2 (1 : ConsensusTower.PrivateAlgebra T (towerBase T))
  rw [map_one, map_one]

noncomputable instance : AddCommGroup (TowerColimit T) where
  add := (· + ·)
  zero := 0
  neg := Neg.neg
  nsmul := nsmulRec
  zsmul := zsmulRec
  add_assoc x y z := by
    obtain ⟨r, X, Y, Z, rfl, rfl, rfl⟩ := exists_rep₃ x y z
    show colimitMk X + colimitMk Y + colimitMk Z =
      colimitMk X + (colimitMk Y + colimitMk Z)
    simp only [colimitMk_add]
    rw [add_assoc]
  zero_add x := by
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    show (0 : TowerColimit T) + colimitMk X = colimitMk X
    rw [← colimitMk_zero r, colimitMk_add, zero_add]
  add_zero x := by
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    show colimitMk X + (0 : TowerColimit T) = colimitMk X
    rw [← colimitMk_zero r, colimitMk_add, add_zero]
  add_comm x y := by
    obtain ⟨r, X, Y, rfl, rfl⟩ := exists_rep₂ x y
    show colimitMk X + colimitMk Y = colimitMk Y + colimitMk X
    simp only [colimitMk_add]
    rw [add_comm]
  neg_add_cancel x := by
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    show -colimitMk X + colimitMk X = (0 : TowerColimit T)
    rw [colimitMk_neg, colimitMk_add, neg_add_cancel, colimitMk_zero]

noncomputable instance : Ring (TowerColimit T) :=
  { (inferInstance : AddCommGroup (TowerColimit T)) with
    mul := (· * ·)
    one := 1
    mul_assoc := fun x y z => by
      obtain ⟨r, X, Y, Z, rfl, rfl, rfl⟩ := exists_rep₃ x y z
      show colimitMk X * colimitMk Y * colimitMk Z =
        colimitMk X * (colimitMk Y * colimitMk Z)
      simp only [colimitMk_mul]
      rw [mul_assoc]
    one_mul := fun x => by
      obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
      show (1 : TowerColimit T) * colimitMk X = colimitMk X
      rw [← colimitMk_one r, colimitMk_mul, one_mul]
    mul_one := fun x => by
      obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
      show colimitMk X * (1 : TowerColimit T) = colimitMk X
      rw [← colimitMk_one r, colimitMk_mul, mul_one]
    left_distrib := fun x y z => by
      obtain ⟨r, X, Y, Z, rfl, rfl, rfl⟩ := exists_rep₃ x y z
      show colimitMk X * (colimitMk Y + colimitMk Z) =
        colimitMk X * colimitMk Y + colimitMk X * colimitMk Z
      simp only [colimitMk_add, colimitMk_mul]
      rw [left_distrib]
    right_distrib := fun x y z => by
      obtain ⟨r, X, Y, Z, rfl, rfl, rfl⟩ := exists_rep₃ x y z
      show (colimitMk X + colimitMk Y) * colimitMk Z =
        colimitMk X * colimitMk Z + colimitMk Y * colimitMk Z
      simp only [colimitMk_add, colimitMk_mul]
      rw [right_distrib]
    zero_mul := fun x => by
      obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
      show (0 : TowerColimit T) * colimitMk X = (0 : TowerColimit T)
      rw [← colimitMk_zero r, colimitMk_mul, zero_mul]
    mul_zero := fun x => by
      obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
      show colimitMk X * (0 : TowerColimit T) = (0 : TowerColimit T)
      rw [← colimitMk_zero r, colimitMk_mul, mul_zero] }

/-- The complex algebra structure of the colimit.  The scalar action is
the raw `colimitSMul` and `algebraMap` lands at the anchor regulator, so
the only registered instance path to `•` runs through this instance. -/
noncomputable instance instTowerColimitAlgebra : Algebra ℂ (TowerColimit T) where
  smul := colimitSMul
  algebraMap :=
    { toFun := fun c =>
        colimitMk (algebraMap ℂ (ConsensusTower.PrivateAlgebra T (towerBase T)) c)
      map_one' := by rw [map_one, colimitMk_one]
      map_mul' := fun a b => by rw [map_mul, ← colimitMk_mul]
      map_zero' := by rw [map_zero, colimitMk_zero]
      map_add' := fun a b => by rw [map_add, ← colimitMk_add] }
  commutes' c x := by
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    obtain ⟨t, h1, h2⟩ := T.commonRefinement (towerBase T) r
    show colimitMk (algebraMap ℂ (ConsensusTower.PrivateAlgebra T (towerBase T)) c) *
        colimitMk X =
      colimitMk X *
        colimitMk (algebraMap ℂ (ConsensusTower.PrivateAlgebra T (towerBase T)) c)
    rw [← colimitMk_refine h1 (algebraMap ℂ _ c), ← colimitMk_refine h2 X,
      AlgHomClass.commutes, colimitMk_mul, colimitMk_mul, Algebra.commutes]
  smul_def' c x := by
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    obtain ⟨t, h1, h2⟩ := T.commonRefinement (towerBase T) r
    show colimitSMul c (colimitMk X) =
      colimitMk (algebraMap ℂ (ConsensusTower.PrivateAlgebra T (towerBase T)) c) *
        colimitMk X
    have hs : colimitSMul c (colimitMk X) = colimitMk (c • X) := rfl
    rw [hs, ← colimitMk_refine h1 (algebraMap ℂ _ c), ← colimitMk_refine h2 X,
      AlgHomClass.commutes, colimitMk_mul, ← Algebra.smul_def,
      ← map_smul (T.algebraRefine h2), colimitMk_refine]

theorem colimitMk_smul {r : ι} (c : ℂ)
    (X : ConsensusTower.PrivateAlgebra T r) :
    c • colimitMk X = colimitMk (c • X) :=
  rfl

/-- The star ring structure of the colimit; the only registered instance
path to `star` runs through this instance. -/
noncomputable instance instTowerColimitStarRing : StarRing (TowerColimit T) where
  star := colimitStar
  star_involutive x := by
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    show colimitStar (colimitStar (colimitMk X)) = colimitMk X
    have hs : ∀ Y : ConsensusTower.PrivateAlgebra T r,
        colimitStar (colimitMk Y) = colimitMk (star Y) := fun _ => rfl
    rw [hs, hs, star_star]
  star_mul x y := by
    obtain ⟨r, X, Y, rfl, rfl⟩ := exists_rep₂ x y
    show colimitStar (colimitMk X * colimitMk Y) =
      colimitStar (colimitMk Y) * colimitStar (colimitMk X)
    have hs : ∀ Z : ConsensusTower.PrivateAlgebra T r,
        colimitStar (colimitMk Z) = colimitMk (star Z) := fun _ => rfl
    rw [colimitMk_mul, hs, hs, hs, colimitMk_mul, star_mul]
  star_add x y := by
    obtain ⟨r, X, Y, rfl, rfl⟩ := exists_rep₂ x y
    show colimitStar (colimitMk X + colimitMk Y) =
      colimitStar (colimitMk X) + colimitStar (colimitMk Y)
    have hs : ∀ Z : ConsensusTower.PrivateAlgebra T r,
        colimitStar (colimitMk Z) = colimitMk (star Z) := fun _ => rfl
    rw [colimitMk_add, hs, hs, hs, colimitMk_add, star_add]

theorem colimitMk_star {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    star (colimitMk X) = colimitMk (star X) :=
  rfl

instance : StarModule ℂ (TowerColimit T) :=
  ⟨fun c x => by
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    show colimitStar (colimitSMul c (colimitMk X)) =
      colimitSMul (star c) (colimitStar (colimitMk X))
    have hs : ∀ Z : ConsensusTower.PrivateAlgebra T r,
        colimitStar (colimitMk Z) = colimitMk (star Z) := fun _ => rfl
    have hm : ∀ (a : ℂ) (Z : ConsensusTower.PrivateAlgebra T r),
        colimitSMul a (colimitMk Z) = colimitMk (a • Z) := fun _ _ => rfl
    rw [hm, hs, hs, hm, star_smul]⟩

theorem colimitMk_algebraMap (r : ι) (c : ℂ) :
    colimitMk (algebraMap ℂ (ConsensusTower.PrivateAlgebra T r) c) =
      algebraMap ℂ (TowerColimit T) c := by
  obtain ⟨t, h1, h2⟩ := T.commonRefinement r (towerBase T)
  refine Quotient.sound ⟨t, h1, h2, ?_⟩
  show T.algebraRefine h1 (algebraMap ℂ (ConsensusTower.PrivateAlgebra T r) c) =
    T.algebraRefine h2
      (algebraMap ℂ (ConsensusTower.PrivateAlgebra T (towerBase T)) c)
  rw [AlgHomClass.commutes, AlgHomClass.commutes]

/-- The canonical cocone map: the inclusion of one regulator's private
algebra into the filtered-colimit observable algebra, as a star algebra
homomorphism. -/
noncomputable def colimitMkHom (T : ConsensusTower ι) (r : ι) :
    ConsensusTower.PrivateAlgebra T r →⋆ₐ[ℂ] TowerColimit T where
  toFun := colimitMk
  map_one' := colimitMk_one r
  map_mul' X Y := (colimitMk_mul X Y).symm
  map_zero' := colimitMk_zero r
  map_add' X Y := (colimitMk_add X Y).symm
  commutes' c := colimitMk_algebraMap r c
  map_star' _ := rfl

section UniversalProperty

variable {B : Type v} [Semiring B] [Algebra ℂ B] [Star B]

/-- The universal property of the filtered colimit, existence half: a
refinement-compatible family of star algebra homomorphisms out of the
tower's levels descends to the colimit. -/
noncomputable def colimitDesc
    (g : ∀ r : ι, ConsensusTower.PrivateAlgebra T r →⋆ₐ[ℂ] B)
    (hg : ∀ {r s : ι} (hrs : r ≤ s) (X : ConsensusTower.PrivateAlgebra T r),
      g s (T.algebraRefine hrs X) = g r X) :
    TowerColimit T →⋆ₐ[ℂ] B where
  toFun := Quotient.lift (fun p : TowerGerm T => g p.1 p.2) (by
    rintro p q ⟨t, hp, hq, e⟩
    show g p.1 p.2 = g q.1 q.2
    rw [← hg hp p.2, ← hg hq q.2, e])
  map_one' := map_one (g (towerBase T))
  map_mul' x y := by
    obtain ⟨r, X, Y, rfl, rfl⟩ := exists_rep₂ x y
    rw [colimitMk_mul]
    exact map_mul (g r) X Y
  map_zero' := map_zero (g (towerBase T))
  map_add' x y := by
    obtain ⟨r, X, Y, rfl, rfl⟩ := exists_rep₂ x y
    rw [colimitMk_add]
    exact map_add (g r) X Y
  commutes' c := by
    rw [← colimitMk_algebraMap (towerBase T) c]
    exact (g (towerBase T)).commutes c
  map_star' x := by
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    rw [colimitMk_star]
    exact map_star (g r) X

@[simp]
theorem colimitDesc_mk
    (g : ∀ r : ι, ConsensusTower.PrivateAlgebra T r →⋆ₐ[ℂ] B)
    (hg : ∀ {r s : ι} (hrs : r ≤ s) (X : ConsensusTower.PrivateAlgebra T r),
      g s (T.algebraRefine hrs X) = g r X)
    {r : ι} (X : ConsensusTower.PrivateAlgebra T r) :
    colimitDesc g hg (colimitMk X) = g r X :=
  rfl

/-- The factorization triangle of the universal property: descending along
a compatible family and precomposing with the canonical cocone map
recovers the family. -/
theorem colimitDesc_comp_colimitMkHom
    (g : ∀ r : ι, ConsensusTower.PrivateAlgebra T r →⋆ₐ[ℂ] B)
    (hg : ∀ {r s : ι} (hrs : r ≤ s) (X : ConsensusTower.PrivateAlgebra T r),
      g s (T.algebraRefine hrs X) = g r X) (r : ι) :
    (colimitDesc g hg).comp (colimitMkHom T r) = g r :=
  StarAlgHom.ext fun _ => rfl

/-- The universal property, uniqueness half: colimit homomorphisms agreeing
on every finite-level class are equal. -/
theorem colimitHom_ext {φ ψ : TowerColimit T →⋆ₐ[ℂ] B}
    (h : ∀ (r : ι) (X : ConsensusTower.PrivateAlgebra T r),
      φ (colimitMk X) = ψ (colimitMk X)) : φ = ψ :=
  StarAlgHom.ext fun x => by
    obtain ⟨r, X, rfl⟩ := colimitMk_surjective x
    exact h r X

end UniversalProperty

/-- **The refinement-null kernel receipt.**  A finite-level observable dies
in the colimit exactly when some declared refinement already kills it: the
construction quotients by precisely the refinement-null kernel, so
collapsing bonding maps are absorbed and nothing else is identified with
zero. -/
theorem colimitMk_eq_zero_iff {r : ι}
    (X : ConsensusTower.PrivateAlgebra T r) :
    colimitMk X = (0 : TowerColimit T) ↔
      ∃ (t : ι) (hrt : r ≤ t), T.algebraRefine hrt X = 0 := by
  constructor
  · intro h
    obtain ⟨t, h1, h2, e⟩ := Quotient.exact h
    refine ⟨t, h1, ?_⟩
    rw [e]
    exact map_zero _
  · rintro ⟨t, hrt, e⟩
    rw [← colimitMk_refine hrt X, e, colimitMk_zero]

/-- The outside proposal's inclusion regime, recovered as a theorem rather
than an assumption: when the bonding maps out of a regulator are injective,
nothing declared at that regulator is lost at the limit.  The collapsing
tower below shows this can genuinely fail, so it is not built into the
construction. -/
theorem colimitMk_injective_of_refine_injective (r : ι)
    (hinj : ∀ (s : ι) (hrs : r ≤ s),
      Function.Injective (T.algebraRefine hrs)) :
    Function.Injective
      (colimitMk (T := T) (r := r)) := by
  intro X Y h
  obtain ⟨t, h1, h2, e⟩ := Quotient.exact h
  exact hinj t h1 e

end TowerColimitCore

/-! ## Limit local algebras: isotony, locality, refinement independence -/

section LimitLocal

variable {N : FiniteCausalObserverNet T}

/-- Eventual regional membership: the colimit element is the class of an
observable declared, at some refinement of the region's regulator, inside
the refined region's local algebra. -/
def MemLimitLocal (N : FiniteCausalObserverNet T) (U : EventRegion N)
    (x : TowerColimit T) : Prop :=
  ∃ (t : ι) (ht : U.regulator ≤ t)
    (X : ConsensusTower.PrivateAlgebra T t),
    X ∈ N.localAlgebra t (N.regionRefine ht U.region) ∧ x = colimitMk X

/-- Two eventually-regional elements admit regional representatives at one
common refinement of the region's regulator. -/
theorem MemLimitLocal.exists_pair {U : EventRegion N} {x y : TowerColimit T}
    (hx : MemLimitLocal N U x) (hy : MemLimitLocal N U y) :
    ∃ (t : ι) (ht : U.regulator ≤ t)
      (X Y : ConsensusTower.PrivateAlgebra T t),
      X ∈ N.localAlgebra t (N.regionRefine ht U.region) ∧
        Y ∈ N.localAlgebra t (N.regionRefine ht U.region) ∧
          x = colimitMk X ∧ y = colimitMk Y := by
  obtain ⟨t₁, h₁, X, hX, rfl⟩ := hx
  obtain ⟨t₂, h₂, Y, hY, rfl⟩ := hy
  obtain ⟨t, ht₁, ht₂⟩ := T.commonRefinement t₁ t₂
  have hX' := N.localAlgebra_natural ht₁ _ X hX
  rw [N.region_refine_trans] at hX'
  have hY' := N.localAlgebra_natural ht₂ _ Y hY
  rw [N.region_refine_trans] at hY'
  exact ⟨t, h₁.trans ht₁, T.algebraRefine ht₁ X, T.algebraRefine ht₂ Y,
    hX', hY', (colimitMk_refine ht₁ X).symm, (colimitMk_refine ht₂ Y).symm⟩

/-- The limit local algebra of an admissible event region: the star
subalgebra of colimit classes eventually carried by the region's declared
refinements.  This is the regional face of deliverable 1: it is built from
the net's own `localAlgebra_natural` and `region_refine_trans` receipts,
with no injectivity assumed anywhere. -/
noncomputable def limitLocalAlgebra (N : FiniteCausalObserverNet T)
    (U : EventRegion N) : StarSubalgebra ℂ (TowerColimit T) where
  carrier := {x | MemLimitLocal N U x}
  mul_mem' := by
    intro x y hx hy
    obtain ⟨t, ht, X, Y, hX, hY, rfl, rfl⟩ := MemLimitLocal.exists_pair hx hy
    exact ⟨t, ht, X * Y, mul_mem hX hY, colimitMk_mul X Y⟩
  add_mem' := by
    intro x y hx hy
    obtain ⟨t, ht, X, Y, hX, hY, rfl, rfl⟩ := MemLimitLocal.exists_pair hx hy
    exact ⟨t, ht, X + Y, add_mem hX hY, colimitMk_add X Y⟩
  one_mem' :=
    ⟨U.regulator, le_rfl, 1, one_mem _, (colimitMk_one _).symm⟩
  zero_mem' :=
    ⟨U.regulator, le_rfl, 0, zero_mem _, (colimitMk_zero _).symm⟩
  algebraMap_mem' := fun c =>
    ⟨U.regulator, le_rfl, algebraMap ℂ _ c,
      StarSubalgebra.algebraMap_mem _ c, (colimitMk_algebraMap _ c).symm⟩
  star_mem' := by
    rintro x ⟨t, ht, X, hX, rfl⟩
    exact ⟨t, ht, star X, star_mem hX, rfl⟩

/-- Finite-level regional observables land in the limit local algebra of
their event region: the colimit inclusion restricts regionally. -/
theorem colimitMk_mem_limitLocalAlgebra (U : EventRegion N)
    (X : N.localAlgebra U.regulator U.region) :
    colimitMk (X : ConsensusTower.PrivateAlgebra T U.regulator) ∈
      limitLocalAlgebra N U := by
  refine ⟨U.regulator, le_rfl, X.1, ?_, rfl⟩
  rw [N.region_refine_refl]
  exact X.2

/-- **Limit-level isotony.**  Along any causal embedding the limit local
algebras are nested.  The proof is the finite isotony receipt pushed
through `region_refine_mono` at a common refinement. -/
theorem limitLocalAlgebra_mono {U V : EventRegion N} (f : U ⟶ V) :
    limitLocalAlgebra N U ≤ limitLocalAlgebra N V := by
  rintro x ⟨t, ht, X, hX, rfl⟩
  obtain ⟨u, htu, hVu⟩ := T.commonRefinement t V.regulator
  have hX' := N.localAlgebra_natural htu _ X hX
  rw [N.region_refine_trans] at hX'
  have hle := N.region_refine_mono hVu (CausalEmbedding.emb f)
  rw [N.region_refine_trans] at hle
  exact ⟨u, hVu, T.algebraRefine htu X, N.isotony u hle hX',
    (colimitMk_refine htu X).symm⟩

/-- **Limit-level locality.**  Colimit classes carried by declared-disjoint
regions at a common regulator commute in the filtered-colimit algebra: the
finite locality receipt survives the limit through `disjoint_natural`. -/
theorem limitLocalAlgebra_commute {r : ι} {U V : N.Region r}
    (hUV : N.disjoint r U V) {x y : TowerColimit T}
    (hx : x ∈ limitLocalAlgebra N ⟨r, U⟩)
    (hy : y ∈ limitLocalAlgebra N ⟨r, V⟩) :
    Commute x y := by
  obtain ⟨t₁, h₁, X, hX, rfl⟩ := hx
  obtain ⟨t₂, h₂, Y, hY, rfl⟩ := hy
  obtain ⟨t, ht₁, ht₂⟩ := T.commonRefinement t₁ t₂
  have hX' := N.localAlgebra_natural ht₁ _ X hX
  rw [N.region_refine_trans] at hX'
  have hY' := N.localAlgebra_natural ht₂ _ Y hY
  rw [N.region_refine_trans] at hY'
  have hd := N.disjoint_natural (h₁.trans ht₁) hUV
  have hc := N.commute_of_disjoint hd ⟨_, hX'⟩ ⟨_, hY'⟩
  rw [← colimitMk_refine ht₁ X, ← colimitMk_refine ht₂ Y]
  have : colimitMk (T.algebraRefine ht₁ X) * colimitMk (T.algebraRefine ht₂ Y) =
      colimitMk (T.algebraRefine ht₂ Y) * colimitMk (T.algebraRefine ht₁ X) := by
    rw [colimitMk_mul, colimitMk_mul, hc]
  exact this

/-- **Refinement independence at the limit.**  Refining an event region
along the regulator preorder does not change its limit local algebra; with
`colimitMk_refine` this is the limit-level refinement-independence
statement. -/
theorem limitLocalAlgebra_refine {r s : ι} (hrs : r ≤ s) (U : N.Region r) :
    limitLocalAlgebra N ⟨s, N.regionRefine hrs U⟩ =
      limitLocalAlgebra N ⟨r, U⟩ := by
  apply le_antisymm
  · rintro x ⟨t, ht, X, hX, rfl⟩
    refine ⟨t, hrs.trans ht, X, ?_, rfl⟩
    rw [N.region_refine_trans hrs ht U] at hX
    exact hX
  · rintro x ⟨t, ht, X, hX, rfl⟩
    obtain ⟨u, htu, hsu⟩ := T.commonRefinement t s
    have hX' := N.localAlgebra_natural htu _ X hX
    rw [N.region_refine_trans] at hX'
    refine ⟨u, hsu, T.algebraRefine htu X, ?_, (colimitMk_refine htu X).symm⟩
    rw [N.region_refine_trans hrs hsu U]
    exact hX'

/-- The colimit is a cocone over the finite-level observable functor: the
functor's image of a causal embedding — the A3 refinement followed by the
isotony inclusion, which is `(observableFunctor N).map f` by
`observableFunctor_map_eq` — is invisible at the limit. -/
theorem colimitMk_observableFunctor_map {U V : EventRegion N} (f : U ⟶ V)
    (X : N.localAlgebra U.regulator U.region) :
    colimitMk (((N.localInclusion (CausalEmbedding.emb f)).comp
        (N.localRefine (CausalEmbedding.le f) U.region) X :
        N.localAlgebra V.regulator V.region) :
          ConsensusTower.PrivateAlgebra T V.regulator) =
      colimitMk (X : ConsensusTower.PrivateAlgebra T U.regulator) :=
  colimitMk_refine (CausalEmbedding.le f) X.1

end LimitLocal

/-! ## The Cauchy embedding class and the time-slice property -/

section TimeSlice

variable {N : FiniteCausalObserverNet T}

/-- The inclusion receipt of a causal embedding at any declared refinement
of its target regulator: the refined source region is declared below the
refined target region. -/
theorem causalEmbedding_refine_le {U V : EventRegion N} (f : U ⟶ V)
    {t : ι} (ht : V.regulator ≤ t) :
    N.regionLE t (N.regionRefine ((CausalEmbedding.le f).trans ht) U.region)
      (N.regionRefine ht V.region) := by
  have h := N.region_refine_mono ht (CausalEmbedding.emb f)
  rwa [N.region_refine_trans] at h

/-- The net's contravariant restriction retracts the bundled isotony
inclusion; this is `restrict_inclusion` spelled with `localInclusion`. -/
theorem restrict_localInclusion {r : ι} {U V : N.Region r}
    (hUV : N.regionLE r U V) (X : N.localAlgebra r U) :
    N.restrict r hUV (N.localInclusion hUV X) = X :=
  N.restrict_inclusion r hUV X

/-- **The Cauchy embedding class.**  A causal embedding is Cauchy when, at
every declared refinement of its target regulator, the net's own
contravariant restriction map inverts the isotony inclusion from the right
as well as from the left.  The left inverse is the `restrict_inclusion`
receipt every embedding enjoys; being Cauchy is exactly the extra right
inverse, so the class is carried by the restriction maps and not by any
injectivity assumption. -/
def IsCauchyEmbedding {U V : EventRegion N} (f : U ⟶ V) : Prop :=
  ∀ (t : ι) (ht : V.regulator ≤ t)
    (Y : N.localAlgebra t (N.regionRefine ht V.region)),
    N.localInclusion (causalEmbedding_refine_le f ht)
      (N.restrict t (causalEmbedding_refine_le f ht) Y) = Y

/-- The Cauchy receipt is equivalent to eventual reverse isotony: beyond
the target regulator, the refined target algebra is contained in the
refined source algebra. -/
theorem isCauchyEmbedding_iff {U V : EventRegion N} (f : U ⟶ V) :
    IsCauchyEmbedding f ↔
      ∀ (t : ι) (ht : V.regulator ≤ t),
        N.localAlgebra t (N.regionRefine ht V.region) ≤
          N.localAlgebra t
            (N.regionRefine ((CausalEmbedding.le f).trans ht) U.region) := by
  constructor
  · intro hC t ht Z hZ
    have h := hC t ht ⟨Z, hZ⟩
    have hval : (N.restrict t (causalEmbedding_refine_le f ht) ⟨Z, hZ⟩ :
        ConsensusTower.PrivateAlgebra T t) = Z := congrArg Subtype.val h
    rw [← hval]
    exact (N.restrict t (causalEmbedding_refine_le f ht) ⟨Z, hZ⟩).2
  · intro hle t ht Y
    have hmem : (Y : ConsensusTower.PrivateAlgebra T t) ∈
        N.localAlgebra t
          (N.regionRefine ((CausalEmbedding.le f).trans ht) U.region) :=
      hle t ht Y.2
    have hY : Y = N.localInclusion (causalEmbedding_refine_le f ht)
        ⟨Y.1, hmem⟩ := Subtype.ext rfl
    rw [hY, restrict_localInclusion]

/-- **The time-slice property at the limit.**  A Cauchy causal embedding
induces an equality of limit local algebras, not merely the isotony
inclusion: at the limit the embedded region already carries every
observable of its target. -/
theorem limitLocalAlgebra_eq_of_isCauchyEmbedding {U V : EventRegion N}
    (f : U ⟶ V) (hf : IsCauchyEmbedding f) :
    limitLocalAlgebra N U = limitLocalAlgebra N V := by
  refine le_antisymm (limitLocalAlgebra_mono f) ?_
  rintro x ⟨t, ht, Y, hY, rfl⟩
  exact ⟨t, (CausalEmbedding.le f).trans ht, Y,
    (isCauchyEmbedding_iff f).mp hf t ht hY, rfl⟩

end TimeSlice

/-! ## Cauchy witnesses in the support-graded example -/

section CauchyWitness

/-- The full two-site region of the support-graded example. -/
noncomputable def fullRegion : EventRegion supportGradedNet :=
  ⟨(), ({0, 1} : Finset (Fin 2))⟩

/-- The causal embedding of the first singleton region into the full
region. -/
noncomputable def topToFull : topRegion ⟶ fullRegion :=
  ⟨le_refl (), Finset.singleton_subset_iff.mpr (Finset.mem_insert_self 0 {1})⟩

/-- On the two-site label set any nonempty region already carries the full
diagonal algebra, so enlarging the first singleton region to the full
region adds no observables. -/
theorem supportDiagonalAlgebra_full_le_singleton :
    supportDiagonalAlgebra ({0, 1} : Finset (Fin 2)) ≤
      supportDiagonalAlgebra ({0} : Finset (Fin 2)) := by
  rintro X ⟨f, rfl, -⟩
  refine ⟨f, rfl, ?_⟩
  intro i hi j hj
  fin_cases i <;> fin_cases j <;> simp_all

/-- **Nonvacuity of the Cauchy class.**  The singleton-to-full embedding of
the support-graded example is Cauchy while relating two distinct declared
regions. -/
theorem topToFull_isCauchyEmbedding : IsCauchyEmbedding topToFull := by
  rw [isCauchyEmbedding_iff]
  intro t ht
  exact supportDiagonalAlgebra_full_le_singleton

/-- **Nontriviality of the Cauchy class.**  The bottom-to-top embedding is
not Cauchy: the separating projector lies in the top algebra but not in the
bottom algebra, so no restriction can split that inclusion. -/
theorem bottomToTop_not_isCauchyEmbedding :
    ¬ IsCauchyEmbedding bottomToTop := by
  intro h
  have hle := (isCauchyEmbedding_iff bottomToTop).mp h () (le_refl ())
  exact e0_not_mem_empty (hle e0_mem_singleton)

/-- The time-slice receipt exercised on the example: the singleton and the
full region share one limit local algebra. -/
theorem limitLocalAlgebra_top_eq_full :
    limitLocalAlgebra supportGradedNet topRegion =
      limitLocalAlgebra supportGradedNet fullRegion :=
  limitLocalAlgebra_eq_of_isCauchyEmbedding topToFull topToFull_isCauchyEmbedding

end CauchyWitness

/-! ## The public-record subfunctor -/

section PublicRecord

variable {N : FiniteCausalObserverNet T}

/-- A refinement-coherent observer family: one observer label per
regulator, compatible with the tower's observer refinement.  Every
committed witness tower has constant observers, so carries one. -/
structure CoherentObserverFamily (T : ConsensusTower ι) : Type u where
  /-- The selected observer at each regulator. -/
  observer : ∀ r : ι, T.Observer r
  /-- The selection is natural for the declared observer refinements. -/
  coherent : ∀ {r s : ι} (hrs : r ≤ s),
    T.observerRefine hrs (observer r) = observer s

/-- The public-record part of a regional algebra: the meet of the region's
local algebra with the selected observer's commutative public algebra. -/
noncomputable def publicLocalAlgebra (N : FiniteCausalObserverNet T)
    (ω : CoherentObserverFamily T) (U : EventRegion N) :
    StarSubalgebra ℂ (ConsensusTower.PrivateAlgebra T U.regulator) :=
  N.localAlgebra U.regulator U.region ⊓
    T.publicAlgebra U.regulator (ω.observer U.regulator)

/-- The action of a causal embedding on public-record parts: the A3
refinement map restricted to the meet, using `localAlgebra_natural` with
isotony on the regional side and `public_natural` with the family's
coherence on the public side. -/
noncomputable def publicLocalRefine (N : FiniteCausalObserverNet T)
    (ω : CoherentObserverFamily T) {U V : EventRegion N} (f : U ⟶ V) :
    publicLocalAlgebra N ω U →⋆ₐ[ℂ] publicLocalAlgebra N ω V where
  toFun X := ⟨T.algebraRefine (CausalEmbedding.le f) X.1, by
    obtain ⟨h1, h2⟩ := StarSubalgebra.mem_inf.mp X.2
    refine StarSubalgebra.mem_inf.mpr ⟨?_, ?_⟩
    · exact N.isotony V.regulator (CausalEmbedding.emb f)
        (N.localAlgebra_natural (CausalEmbedding.le f) U.region X.1 h1)
    · rw [← ω.coherent (CausalEmbedding.le f)]
      exact T.public_natural (CausalEmbedding.le f)
        (ω.observer U.regulator) X.1 h2⟩
  map_one' := Subtype.ext (map_one (T.algebraRefine (CausalEmbedding.le f)))
  map_mul' X Y :=
    Subtype.ext (map_mul (T.algebraRefine (CausalEmbedding.le f)) X.1 Y.1)
  map_zero' := Subtype.ext (map_zero (T.algebraRefine (CausalEmbedding.le f)))
  map_add' X Y :=
    Subtype.ext (map_add (T.algebraRefine (CausalEmbedding.le f)) X.1 Y.1)
  commutes' c :=
    Subtype.ext (AlgHomClass.commutes (T.algebraRefine (CausalEmbedding.le f)) c)
  map_star' X :=
    Subtype.ext (map_star (T.algebraRefine (CausalEmbedding.le f)) X.1)

/-- **The public-record subfunctor.**  Each admissible event region goes to
the public part of its local algebra; a causal embedding acts by the
restricted refinement map.  The functor laws are again the tower's
`algebra_refine_refl` and `algebra_refine_trans` receipts. -/
noncomputable def publicObservableFunctor (N : FiniteCausalObserverNet T)
    (ω : CoherentObserverFamily T) : EventRegion N ⥤ ObservableAlgebra T where
  obj U := ⟨U.regulator, publicLocalAlgebra N ω U⟩
  map {U V} f := publicLocalRefine N ω f
  map_id U := by
    refine StarAlgHom.ext fun X => Subtype.ext ?_
    exact T.algebra_refine_refl U.regulator X.1
  map_comp {U V W} f g := by
    refine StarAlgHom.ext fun X => Subtype.ext ?_
    exact (T.algebra_refine_trans (CausalEmbedding.le f)
      (CausalEmbedding.le g) X.1).symm

/-- The subfunctor receipt: the public-record functor includes into the
observable functor by the componentwise star-subalgebra inclusion, and the
inclusion is natural in the causal embedding. -/
noncomputable def publicRecordInclusion (N : FiniteCausalObserverNet T)
    (ω : CoherentObserverFamily T) :
    publicObservableFunctor N ω ⟶ observableFunctor N where
  app _ := StarSubalgebra.inclusion inf_le_left
  naturality _ _ _ := StarAlgHom.ext fun _ => Subtype.ext rfl

/-- The subfunctor component at an event region is precisely the
star-subalgebra inclusion of the public part. -/
theorem publicRecordInclusion_app_eq (N : FiniteCausalObserverNet T)
    (ω : CoherentObserverFamily T) (U : EventRegion N) :
    (publicRecordInclusion N ω).app U =
      (StarSubalgebra.inclusion inf_le_left :
        publicLocalAlgebra N ω U →⋆ₐ[ℂ]
          N.localAlgebra U.regulator U.region) :=
  rfl

/-- With `publicRecordInclusion_app_eq`, every component of the subfunctor
inclusion is injective, so the public-record functor is a genuine
subobject of the observable functor at each admissible event region. -/
theorem publicRecordInclusion_app_injective (N : FiniteCausalObserverNet T)
    (ω : CoherentObserverFamily T) (U : EventRegion N) :
    Function.Injective
      (StarSubalgebra.inclusion
        (inf_le_left :
          publicLocalAlgebra N ω U ≤ N.localAlgebra U.regulator U.region)) :=
  StarSubalgebra.inclusion_injective inf_le_left

/-- Eventual public membership: the colimit element is the class of an
observable declared in the selected observer's public algebra at some
regulator. -/
def MemLimitPublic (ω : CoherentObserverFamily T) (x : TowerColimit T) : Prop :=
  ∃ (r : ι) (X : ConsensusTower.PrivateAlgebra T r),
    X ∈ T.publicAlgebra r (ω.observer r) ∧ x = colimitMk X

/-- Two eventually-public elements admit public representatives at one
common regulator, by `public_natural` and the family's coherence. -/
theorem MemLimitPublic.exists_pair {ω : CoherentObserverFamily T}
    {x y : TowerColimit T}
    (hx : MemLimitPublic ω x) (hy : MemLimitPublic ω y) :
    ∃ (t : ι) (X Y : ConsensusTower.PrivateAlgebra T t),
      X ∈ T.publicAlgebra t (ω.observer t) ∧
        Y ∈ T.publicAlgebra t (ω.observer t) ∧
          x = colimitMk X ∧ y = colimitMk Y := by
  obtain ⟨r₁, X, hX, rfl⟩ := hx
  obtain ⟨r₂, Y, hY, rfl⟩ := hy
  obtain ⟨t, h₁, h₂⟩ := T.commonRefinement r₁ r₂
  have hX' := T.public_natural h₁ (ω.observer r₁) X hX
  rw [ω.coherent h₁] at hX'
  have hY' := T.public_natural h₂ (ω.observer r₂) Y hY
  rw [ω.coherent h₂] at hY'
  exact ⟨t, T.algebraRefine h₁ X, T.algebraRefine h₂ Y, hX', hY',
    (colimitMk_refine h₁ X).symm, (colimitMk_refine h₂ Y).symm⟩

/-- The limit public-record algebra of a coherent observer family: colimit
classes eventually carried by the family's commutative public algebras. -/
noncomputable def limitPublicAlgebra (ω : CoherentObserverFamily T) :
    StarSubalgebra ℂ (TowerColimit T) where
  carrier := {x | MemLimitPublic ω x}
  mul_mem' := by
    intro x y hx hy
    obtain ⟨t, X, Y, hX, hY, rfl, rfl⟩ := MemLimitPublic.exists_pair hx hy
    exact ⟨t, X * Y, mul_mem hX hY, colimitMk_mul X Y⟩
  add_mem' := by
    intro x y hx hy
    obtain ⟨t, X, Y, hX, hY, rfl, rfl⟩ := MemLimitPublic.exists_pair hx hy
    exact ⟨t, X + Y, add_mem hX hY, colimitMk_add X Y⟩
  one_mem' :=
    ⟨towerBase T, 1, one_mem _, (colimitMk_one _).symm⟩
  zero_mem' :=
    ⟨towerBase T, 0, zero_mem _, (colimitMk_zero _).symm⟩
  algebraMap_mem' := fun c =>
    ⟨towerBase T, algebraMap ℂ _ c, StarSubalgebra.algebraMap_mem _ c,
      (colimitMk_algebraMap _ c).symm⟩
  star_mem' := by
    rintro x ⟨r, X, hX, rfl⟩
    exact ⟨r, star X, star_mem hX, rfl⟩

/-- **Limit-level public commutativity.**  The tower's finite
`public_mul_comm` receipt survives the colimit: the limit public-record
algebra is commutative. -/
theorem limitPublicAlgebra_commute (ω : CoherentObserverFamily T)
    {x y : TowerColimit T}
    (hx : x ∈ limitPublicAlgebra ω) (hy : y ∈ limitPublicAlgebra ω) :
    Commute x y := by
  obtain ⟨t, X, Y, hX, hY, rfl, rfl⟩ := MemLimitPublic.exists_pair hx hy
  have hc := T.public_mul_comm t (ω.observer t) ⟨X, hX⟩ ⟨Y, hY⟩
  have hval : X * Y = Y * X := congrArg Subtype.val hc
  have : colimitMk X * colimitMk Y = colimitMk Y * colimitMk X := by
    rw [colimitMk_mul, colimitMk_mul, hval]
  exact this

/-- Finite-level public-record observables land in the limit public
algebra: the subfunctor is compatible with the colimit. -/
theorem colimitMk_mem_limitPublicAlgebra (ω : CoherentObserverFamily T)
    (U : EventRegion N) (X : publicLocalAlgebra N ω U) :
    colimitMk (X : ConsensusTower.PrivateAlgebra T U.regulator) ∈
      limitPublicAlgebra ω :=
  ⟨U.regulator, X.1, (StarSubalgebra.mem_inf.mp X.2).2, rfl⟩

end PublicRecord

/-! ## Nonvacuity of the public-record subfunctor on the witness tower -/

section PublicRecordWitness

/-- The constant observer family of the witness tower. -/
noncomputable def witnessObserverFamily : CoherentObserverFamily witnessTower where
  observer := fun _ => ()
  coherent := fun _ => rfl

/-- **Nonvacuity and properness of the subfunctor.**  At the top region of
the support-graded example over the witness tower, the public-record part
is properly below the full regional algebra: the separating projector is
regional but not public. -/
theorem publicLocalAlgebra_lt_localAlgebra :
    publicLocalAlgebra supportGradedNet witnessObserverFamily topRegion <
      supportGradedNet.localAlgebra topRegion.regulator topRegion.region := by
  refine lt_of_le_of_ne inf_le_left fun h => ?_
  have hmem : e0 ∈
      publicLocalAlgebra supportGradedNet witnessObserverFamily topRegion := by
    rw [h]
    exact e0_mem_singleton
  obtain ⟨c, hc⟩ := (oneBlockPartition_mem_iff e0).mp hmem.2
  have h00 := congrFun (congrFun hc 0) 0
  have h11 := congrFun (congrFun hc 1) 1
  simp [e0, Matrix.smul_apply] at h00 h11
  rw [h00] at h11
  exact one_ne_zero h11

end PublicRecordWitness

/-! ## The collapsing tower: no inclusion model can express the corpus -/

section CollapsingWitness

instance : Subsingleton (Matrix (Fin 0) (Fin 0) ℂ) :=
  ⟨fun M N => by
    ext i j
    exact i.elim0⟩

/-- The unique star algebra homomorphism into the trivial zero-by-zero
matrix algebra, in which `1 = 0`. -/
noncomputable def collapseHom :
    Matrix (Fin 1) (Fin 1) ℂ →⋆ₐ[ℂ] Matrix (Fin 0) (Fin 0) ℂ where
  toFun _ := 0
  map_one' := Subsingleton.elim _ _
  map_mul' _ _ := Subsingleton.elim _ _
  map_zero' := Subsingleton.elim _ _
  map_add' _ _ := Subsingleton.elim _ _
  commutes' _ := Subsingleton.elim _ _
  map_star' _ := Subsingleton.elim _ _

/-- **A collapsing consensus tower.**  Dimension one at the coarse
regulator, dimension zero at the fine one, no observers; the sole
nontrivial bonding map collapses everything.  Every structure field is
discharged from the corpus's own interface, so the corpus genuinely
admits collapsing towers, and the outside proposal's inclusion model
cannot express this inhabitant. -/
noncomputable def collapsingTower : ConsensusTower Bool where
  indexNonempty := ⟨false⟩
  commonRefinement := fun r s => ⟨true, Bool.le_true r, Bool.le_true s⟩
  dim := fun b => cond b 0 1
  Observer := fun _ => PEmpty
  observerFinite := fun _ => inferInstance
  Record := fun _ => PUnit
  recordFinite := fun _ => inferInstance
  recordOrder := fun _ o => o.elim
  publicAlgebra := fun _ o => o.elim
  public_mul_comm := fun _ o => o.elim
  recordElement := fun _ o => o.elim
  state := fun _ o => o.elim
  state_isState := fun _ o => o.elim
  generator := fun _ o => o.elim
  observerRefine := fun _ o => o.elim
  recordRefine := fun _ _ => PUnit.unit
  algebraRefine := @fun r s h =>
    match r, s, h with
    | false, false, _ => StarAlgHom.id ℂ (Matrix (Fin 1) (Fin 1) ℂ)
    | false, true, _ => collapseHom
    | true, true, _ => StarAlgHom.id ℂ (Matrix (Fin 0) (Fin 0) ℂ)
    | true, false, h => absurd h (by decide)
  observer_refine_refl := fun _ o => o.elim
  observer_refine_trans := fun _ _ o => o.elim
  record_refine_refl := fun _ x => Subsingleton.elim _ _
  record_refine_trans := fun _ _ _ => Subsingleton.elim _ _
  algebra_refine_refl := by
    intro r X
    cases r <;> rfl
  algebra_refine_trans := by
    intro r s t hrs hst X
    cases r <;> cases s <;> cases t <;>
      first
        | rfl
        | exact Subsingleton.elim _ _
        | exact absurd hrs (by decide)
        | exact absurd hst (by decide)
  record_order_natural := fun _ o => o.elim
  public_natural := fun _ o => o.elim
  record_element_natural := fun _ o => o.elim
  state_natural := fun _ o => o.elim
  generator_natural := fun _ o => o.elim

/-- The collapsing bonding map is not injective: this tower is not a tower
of inclusions, which is exactly the case the outside proposal's model
assumed away. -/
theorem collapsingTower_refine_not_injective :
    ¬ Function.Injective
      (collapsingTower.algebraRefine (show (false : Bool) ≤ true by decide)) := by
  intro hinj
  have hcol :
      collapsingTower.algebraRefine (show (false : Bool) ≤ true by decide)
          (0 : ConsensusTower.PrivateAlgebra collapsingTower false) =
        collapsingTower.algebraRefine (show (false : Bool) ≤ true by decide)
          (1 : ConsensusTower.PrivateAlgebra collapsingTower false) :=
    Subsingleton.elim (α := Matrix (Fin 0) (Fin 0) ℂ) _ _
  have h01 := hinj hcol
  have := congrFun (congrFun h01 ⟨0, by decide⟩) ⟨0, by decide⟩
  simp at this

/-- **The refinement-null kernel is genuinely nontrivial in the collapsing
regime.**  The coarse unit is a nonzero finite-level observable whose
colimit class vanishes; `colimitMk_eq_zero_iff` names exactly which germs
die, so the quotient construction absorbs this tower where an inclusion
model cannot. -/
theorem collapsingTower_one_refinementNull :
    colimitMk (1 : ConsensusTower.PrivateAlgebra collapsingTower false) =
        (0 : TowerColimit collapsingTower) ∧
      (1 : ConsensusTower.PrivateAlgebra collapsingTower false) ≠ 0 := by
  constructor
  · rw [colimitMk_eq_zero_iff]
    exact ⟨true, by decide,
      Subsingleton.elim (α := Matrix (Fin 0) (Fin 0) ℂ) _ _⟩
  · intro h
    have := congrFun (congrFun h ⟨0, by decide⟩) ⟨0, by decide⟩
    simp at this

/-- Downstream of a total collapse the whole colimit algebra is trivial:
the limit object records faithfully that the tower's data dies. -/
theorem collapsingTower_colimit_trivial :
    (1 : TowerColimit collapsingTower) = 0 := by
  rw [← colimitMk_one false]
  exact collapsingTower_one_refinementNull.1

/-- In the corpus, a regulator that carries even one observer carries a
selected density state, whose unit trace forces a positive matrix
dimension: total collapse cannot happen where an observer watches. -/
theorem dim_pos_of_observer (T : ConsensusTower ι) (r : ι)
    (ho : Nonempty (T.Observer r)) : 0 < T.dim r := by
  obtain ⟨o⟩ := ho
  rcases Nat.eq_zero_or_pos (T.dim r) with h | h
  · exfalso
    have hst := (T.state_isState r o).2
    have h0 : (T.state r o).trace = 0 := by
      rw [Matrix.trace]
      apply Finset.sum_eq_zero
      intro i _
      exact absurd i.2 (by omega)
    rw [hst] at h0
    exact one_ne_zero h0
  · exact h

/-- **Where observers exist, the corpus forces injective bonding maps.**
If the target regulator carries an observer, its selected state forces a
positive dimension, matrix algebras over `ℂ` are simple, and a unital ring
homomorphism out of a simple ring into a nontrivial ring is injective.  So
a nontrivial refinement-null kernel is confined to observer-free
regulators — exactly where `collapsingTower` lives — and the colimit
construction above assumes neither regime. -/
theorem algebraRefine_injective_of_observer (T : ConsensusTower ι)
    {r s : ι} (hrs : r ≤ s) (ho : Nonempty (T.Observer s)) :
    Function.Injective (T.algebraRefine hrs) := by
  have hs := dim_pos_of_observer T s ho
  rcases Nat.eq_zero_or_pos (T.dim r) with h0 | hr
  · intro X Y _
    ext i j
    exact absurd i.2 (by omega)
  · haveI : Nonempty (Fin (T.dim r)) := ⟨⟨0, hr⟩⟩
    haveI : Nontrivial (ConsensusTower.PrivateAlgebra T s) := by
      refine ⟨0, 1, fun h => ?_⟩
      have := congrFun (congrFun h ⟨0, hs⟩) ⟨0, hs⟩
      simp at this
    have hinj := RingHom.injective
      (RingHomClass.toRingHom (T.algebraRefine hrs))
    exact fun X Y h => hinj h

end CollapsingWitness


end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.eventRegionCategory
#print axioms OPH.QFT.observableAlgebraCategory
#print axioms OPH.QFT.observableFunctor
#print axioms OPH.QFT.observableFunctor_map_eq
#print axioms OPH.QFT.supportGradedNet
#print axioms OPH.QFT.supportGradedNet_has_disjoint_pair
#print axioms OPH.QFT.supportDiagonalAlgebra_bot_lt
#print axioms OPH.QFT.observableFunctor_obj_strict
#print axioms OPH.QFT.observableFunctor_map_not_surjective
#print axioms OPH.QFT.colimitMkHom
#print axioms OPH.QFT.colimitMk_refine
#print axioms OPH.QFT.colimitMk_eq_colimitMk_iff
#print axioms OPH.QFT.colimitDesc_mk
#print axioms OPH.QFT.colimitDesc_comp_colimitMkHom
#print axioms OPH.QFT.colimitHom_ext
#print axioms OPH.QFT.colimitMk_eq_zero_iff
#print axioms OPH.QFT.colimitMk_injective_of_refine_injective
#print axioms OPH.QFT.colimitMk_mem_limitLocalAlgebra
#print axioms OPH.QFT.limitLocalAlgebra_mono
#print axioms OPH.QFT.limitLocalAlgebra_commute
#print axioms OPH.QFT.limitLocalAlgebra_refine
#print axioms OPH.QFT.colimitMk_observableFunctor_map
#print axioms OPH.QFT.isCauchyEmbedding_iff
#print axioms OPH.QFT.limitLocalAlgebra_eq_of_isCauchyEmbedding
#print axioms OPH.QFT.topToFull_isCauchyEmbedding
#print axioms OPH.QFT.bottomToTop_not_isCauchyEmbedding
#print axioms OPH.QFT.limitLocalAlgebra_top_eq_full
#print axioms OPH.QFT.publicObservableFunctor
#print axioms OPH.QFT.publicRecordInclusion
#print axioms OPH.QFT.publicRecordInclusion_app_eq
#print axioms OPH.QFT.publicRecordInclusion_app_injective
#print axioms OPH.QFT.limitPublicAlgebra_commute
#print axioms OPH.QFT.colimitMk_mem_limitPublicAlgebra
#print axioms OPH.QFT.publicLocalAlgebra_lt_localAlgebra
#print axioms OPH.QFT.collapsingTower_refine_not_injective
#print axioms OPH.QFT.collapsingTower_one_refinementNull
#print axioms OPH.QFT.collapsingTower_colimit_trivial
#print axioms OPH.QFT.dim_pos_of_observer
#print axioms OPH.QFT.algebraRefine_injective_of_observer
