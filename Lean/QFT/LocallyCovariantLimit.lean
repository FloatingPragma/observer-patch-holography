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

## Claim boundary

This module constructs the finite-level category and functor only.  No
filtered-colimit observable algebra, no limit-level isotony, locality,
or refinement-independence theorem, no time-slice property, no Cauchy
embedding class, no public-record subfunctor, and no naturality
statement for one is constructed here; those are the declared later E3
deliverables.  No CP, CPTP, or adaptive physical locality property is
inferred from B2 Kraus normalisation or from the E1 finite interface.
The causal embeddings are declared order-theoretic data of the E1 net;
no Lorentzian or spacetime reading is attached.
-/

namespace OPH.QFT

open OPH.Tower
open CategoryTheory

universe u

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
