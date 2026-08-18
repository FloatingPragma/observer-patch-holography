import QFT.TripleCarrierJoin
import QFT.StructuralNetAdequacySurface

/-!
# A typed hinge-fibered join of the two committed carriers, and evolved-data time slices

`QFT.TowerAnchoredCorrelation` proves the exact observer-86 marginal hinge
agreement between the committed 86/88 carrier (ambient dimension 182) and the
committed 86/247 carrier (ambient dimension 169); `QFT.TripleCarrierJoin`
supplies a common product carrier with marginal channels.  Neither module
constructs a typed object into which **both** committed carriers embed at
their diagonal (record) function layer.  This module constructs one, and
settles the region-family question of the corrected
one-step evolution interface of `QFT.PathTimeSliceInterface`.

**Part 1: the hinge-fibered join (finite linear algebra on committed data
only).**  The carriers are taken at their diagonal (record) layer: complex
functions on the committed 182-point and 169-point ambient index sets.  Every
committed counted state is diagonal, so both counted states and the shared
observer-86 marginal live exactly on this layer.  The observer-86 hinge is the
13-dimensional function space on the shared label alphabet; it maps into each
carrier by the fiber-uniform section of slot-marginalization.  The join is
the pushout: the product of the two carriers modulo the antidiagonal graph of
the two hinge sections.

* Both embeddings are typed linear maps and injective; they agree exactly on
  the hinge (`embed_agree_on_hinge`).
* The join has dimension `182 + 169 - 13 = 338`, strictly below the direct-sum
  dimension 351, and the intersection of the two embedded images is exactly
  the 13-dimensional hinge image (`embed_range_inter`): the join genuinely
  glues along the hinge rather than placing the carriers side by side.
* The hinge restriction map retracts the join onto the hinge; under it, the
  two transported counted states restrict to the **same** hinge datum, the
  committed `count86` occupation law (`joinStates_marginal_compatible`), whose
  diagonal matrix is the committed `marginal86State`.
* The committed ambient anchoring equivalences `ambientEquiv` and
  `ambientEquiv247` act on the diagonal layer exactly as the committed basis
  reindexings, and the join transports along them with commuting squares on
  both embeddings (`joinAnchor_embed88`, `joinAnchor_embed247`).
* The committed walk-step transports (`pair88StepEvolve` on the 86/88 side and
  `stepEvolve` on the 86/247 side, over the 31 actual adjacent source
  transitions, no cyclic edge) act on the diagonal layer as the function-level
  evolutions constructed here, share the observer-86 hinge action, and descend
  to one join evolution intertwined by both embeddings and by the hinge
  restriction (`joinStep_embed88`, `joinStep_embed247`,
  `hingeRestrict_joinStep`).

**Part 2: evolved-data time slices on the corrected interface.**  The
corrected `TimeIndexedNetEvolution` interface generates the next regional
algebra from the evolved admitted data.  Two exact results delimit and then
inhabit its region-family freedom:

* **Delimitation.**  Any evolution family that maps each regional algebra of
  the committed 86/247 diamond onto itself forces a step-invariant region
  family; in particular every net-compatible slice family over the committed
  walk evolution `stepEvolve` on cyclic time `Fin 32` is constant
  (`stepEvolve_slice_constant`).  For an arbitrary evolution family, the
  bottom and top regions are rigid at every step, so a region family can
  change only by exchanging the two slot regions
  (`slice_change_forces_slot_exchange`).
* **Existence.**  `slotAlternatingEvolution` inhabits the corrected interface
  with a genuinely time-dependent region family: the slice alternates between
  the left and right slot regions, the per-step automorphism is the committed
  walk-step evolution composed with the committed internal `slotSwap`, and the
  admitted slice data are the committed source generators of observer 86
  (left slot, even steps) and observer 247 (right slot, odd steps).  The
  transport-aware generation clause holds at every step.  The former
  forced-step-invariance reading is therefore wrong in both directions: it
  fails for slot-exchanging evolutions and holds exactly for region-preserving
  ones.

**What is NOT claimed.**  The join is a finite typed linear object on the
diagonal layer: no star-algebra amalgamation, free product, or C*-pushout is
constructed, the quotient carries no order or positivity structure, and the
transported counted states remain source-produced data of their own carriers.
The join is not spacetime covariance, not a physical Cauchy surface, and not a
continuum limit; no uniqueness of the join is claimed.  The hinge sections are
the fiber-uniform lifts, a declared normalization.  In Part 2 the slot
exchange `slotSwap` is a constructed internal algebra symmetry, not a source
transition, and the `stepEvolve` factor carries the committed declared
terminal-to-initial bookkeeping closure at index 31; cyclic time, the slot
assembly, and the label conventions stay declared postprocessors over
post-hoc transcriptions (premise rows PR-32, PR-33, PR-34 as bundled in the
adequacy surface).  Physical time, causality, and spacetime attachment stay
open under PR-52; the source-selected time-slice content of PR-58 stays open:
this module types and delimits the finite interface, it does not discharge
the premise.
-/

set_option autoImplicit false
set_option maxRecDepth 65536

namespace OPH.QFT

open Matrix
open Kronecker

noncomputable section

/-! ## The diagonal-layer carriers and the observer-86 hinge -/

/-- The observer-86 hinge space: complex functions on the shared
thirteen-label alphabet.  This is the layer on which the committed
`count86` marginal law lives. -/
abbrev HingeFn := Fin 13 → ℂ

/-- The 86/88 diagonal-layer carrier: complex functions on the committed
182-point ambient index set. -/
abbrev CarrierFn88 := PairIndex → ℂ

/-- The 86/247 diagonal-layer carrier: complex functions on the committed
169-point ambient index set. -/
abbrev CarrierFn247 := PairIndex247 → ℂ

/-- Marginalization of the 86/88 carrier onto the observer-86 hinge:
sum over the observer-88 fiber. -/
def marg88 : CarrierFn88 →ₗ[ℂ] HingeFn where
  toFun v := fun i => ∑ j : Fin 14, v (i, j)
  map_add' u v := by
    funext i
    simp [Finset.sum_add_distrib]
  map_smul' c v := by
    funext i
    simp [Finset.mul_sum]

/-- Marginalization of the 86/247 carrier onto the observer-86 hinge:
sum over the observer-247 fiber. -/
def marg247 : CarrierFn247 →ₗ[ℂ] HingeFn where
  toFun v := fun i => ∑ j : Fin 13, v (i, j)
  map_add' u v := by
    funext i
    simp [Finset.sum_add_distrib]
  map_smul' c v := by
    funext i
    simp [Finset.mul_sum]

theorem marg88_apply (v : CarrierFn88) (i : Fin 13) :
    marg88 v i = ∑ j : Fin 14, v (i, j) := rfl

theorem marg247_apply (v : CarrierFn247) (i : Fin 13) :
    marg247 v i = ∑ j : Fin 13, v (i, j) := rfl

/-- The fiber-uniform section of `marg88`: spread each hinge value evenly
over the 14-point observer-88 fiber. -/
def sect88 : HingeFn →ₗ[ℂ] CarrierFn88 where
  toFun h := fun p => h p.1 / 14
  map_add' u v := by
    funext p
    simp [add_div]
  map_smul' c h := by
    funext p
    simp [mul_div_assoc]

/-- The fiber-uniform section of `marg247`: spread each hinge value evenly
over the 13-point observer-247 fiber. -/
def sect247 : HingeFn →ₗ[ℂ] CarrierFn247 where
  toFun h := fun p => h p.1 / 13
  map_add' u v := by
    funext p
    simp [add_div]
  map_smul' c h := by
    funext p
    simp [mul_div_assoc]

theorem sect88_apply (h : HingeFn) (p : PairIndex) :
    sect88 h p = h p.1 / 14 := rfl

theorem sect247_apply (h : HingeFn) (p : PairIndex247) :
    sect247 h p = h p.1 / 13 := rfl

/-- The uniform lift is a section: marginalizing it recovers the hinge
datum exactly. -/
theorem marg88_sect88 (h : HingeFn) : marg88 (sect88 h) = h := by
  funext i
  rw [marg88_apply]
  simp only [sect88_apply]
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  push_cast
  field_simp

/-- The uniform lift is a section on the 86/247 side. -/
theorem marg247_sect247 (h : HingeFn) : marg247 (sect247 h) = h := by
  funext i
  rw [marg247_apply]
  simp only [sect247_apply]
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  push_cast
  field_simp

theorem sect88_injective : Function.Injective sect88 :=
  Function.LeftInverse.injective marg88_sect88

theorem sect247_injective : Function.Injective sect247 :=
  Function.LeftInverse.injective marg247_sect247

/-! ## The hinge-fibered join carrier -/

/-- The gluing submodule: the antidiagonal graph of the two hinge
sections inside the product of the two carriers. -/
def hingeGraph : Submodule ℂ (CarrierFn88 × CarrierFn247) :=
  LinearMap.range (sect88.prod (-sect247))

theorem hingeGraph_gen_apply (h : HingeFn) :
    (sect88.prod (-sect247)) h = (sect88 h, -(sect247 h)) := rfl

/-- **The join carrier**: the pushout of the two committed diagonal-layer
carriers over the observer-86 hinge. -/
abbrev JoinCarrier := (CarrierFn88 × CarrierFn247) ⧸ hingeGraph

/-- The typed embedding of the 86/88 carrier into the join. -/
def embed88 : CarrierFn88 →ₗ[ℂ] JoinCarrier :=
  hingeGraph.mkQ.comp (LinearMap.inl ℂ CarrierFn88 CarrierFn247)

/-- The typed embedding of the 86/247 carrier into the join. -/
def embed247 : CarrierFn247 →ₗ[ℂ] JoinCarrier :=
  hingeGraph.mkQ.comp (LinearMap.inr ℂ CarrierFn88 CarrierFn247)

theorem embed88_apply (v : CarrierFn88) :
    embed88 v = Submodule.Quotient.mk (v, 0) := rfl

theorem embed247_apply (w : CarrierFn247) :
    embed247 w = Submodule.Quotient.mk (0, w) := rfl

/-- **Hinge agreement.**  The two embeddings agree exactly on the
observer-86 hinge: the lift of any hinge datum lands on the same join
element through either carrier. -/
theorem embed_agree_on_hinge (h : HingeFn) :
    embed88 (sect88 h) = embed247 (sect247 h) := by
  rw [embed88_apply, embed247_apply, Submodule.Quotient.eq]
  refine ⟨h, ?_⟩
  rw [hingeGraph_gen_apply]
  simp [Prod.mk_sub_mk]

/-- The 86/88 embedding is injective. -/
theorem embed88_injective : Function.Injective embed88 := by
  intro u v huv
  rw [embed88_apply, embed88_apply, Submodule.Quotient.eq] at huv
  obtain ⟨h, hh⟩ := huv
  rw [hingeGraph_gen_apply] at hh
  have hsnd : -(sect247 h) = 0 := by
    have := congrArg Prod.snd hh
    simpa using this
  have hzero : h = 0 := by
    apply sect247_injective
    rw [map_zero]
    exact neg_eq_zero.mp hsnd
  have hfst : sect88 h = u - v := by
    have := congrArg Prod.fst hh
    simpa using this
  rw [hzero, map_zero] at hfst
  exact sub_eq_zero.mp hfst.symm

/-- The 86/247 embedding is injective. -/
theorem embed247_injective : Function.Injective embed247 := by
  intro u v huv
  rw [embed247_apply, embed247_apply, Submodule.Quotient.eq] at huv
  obtain ⟨h, hh⟩ := huv
  rw [hingeGraph_gen_apply] at hh
  have hfst : sect88 h = 0 := by
    have := congrArg Prod.fst hh
    simpa using this
  have hzero : h = 0 := by
    apply sect88_injective
    rw [map_zero]
    exact hfst
  have hsnd : -(sect247 h) = u - v := by
    have := congrArg Prod.snd hh
    simpa using this
  rw [hzero, map_zero, neg_zero] at hsnd
  exact sub_eq_zero.mp hsnd.symm

/-! ## The hinge restriction of the join -/

theorem hingeGraph_le_ker :
    hingeGraph ≤ LinearMap.ker (marg88.coprod marg247) := by
  rintro x ⟨h, rfl⟩
  rw [LinearMap.mem_ker, hingeGraph_gen_apply, LinearMap.coprod_apply,
    map_neg, marg88_sect88, marg247_sect247]
  exact add_neg_cancel h

/-- The hinge restriction of the join: the descent of the summed
marginalizations.  It is well defined exactly because the two hinge
sections are marginal-normalized. -/
def hingeRestrict : JoinCarrier →ₗ[ℂ] HingeFn :=
  hingeGraph.liftQ (marg88.coprod marg247) hingeGraph_le_ker

theorem hingeRestrict_mk (x : CarrierFn88 × CarrierFn247) :
    hingeRestrict (Submodule.Quotient.mk x) = marg88 x.1 + marg247 x.2 := by
  show (marg88.coprod marg247) x = _
  rw [LinearMap.coprod_apply]

/-- Restricting an embedded 86/88 element to the hinge is
marginalization. -/
theorem hingeRestrict_embed88 (v : CarrierFn88) :
    hingeRestrict (embed88 v) = marg88 v := by
  rw [embed88_apply, hingeRestrict_mk]
  simp

/-- Restricting an embedded 86/247 element to the hinge is
marginalization. -/
theorem hingeRestrict_embed247 (w : CarrierFn247) :
    hingeRestrict (embed247 w) = marg247 w := by
  rw [embed247_apply, hingeRestrict_mk]
  simp

/-- The hinge restriction retracts the join onto the hinge. -/
theorem hingeRestrict_retraction (h : HingeFn) :
    hingeRestrict (embed88 (sect88 h)) = h := by
  rw [hingeRestrict_embed88, marg88_sect88]

/-! ## The committed counted states on the join -/

/-- The committed 86/88 counted state, read at the diagonal layer. -/
def stateFn88 : CarrierFn88 := fun p => (jointCount88 p : ℂ) / 32

/-- The committed 86/247 counted state, read at the diagonal layer. -/
def stateFn247 : CarrierFn247 := fun p => (jointCount p : ℂ) / 32

/-- The committed observer-86 marginal law, read at the hinge layer. -/
def hingeFn86 : HingeFn := fun i => (count86 i : ℂ) / 32

/-- The diagonal matrix of `stateFn88` is the committed counted 86/88
state. -/
theorem diagonal_stateFn88 : Matrix.diagonal stateFn88 = correlationState88 :=
  rfl

/-- The diagonal matrix of `stateFn247` is the committed counted 86/247
state. -/
theorem diagonal_stateFn247 : Matrix.diagonal stateFn247 = correlationState :=
  rfl

/-- The diagonal matrix of `hingeFn86` is the committed observer-86
marginal state. -/
theorem diagonal_hingeFn86 : Matrix.diagonal hingeFn86 = marginal86State :=
  rfl

/-- Count-level marginalization of the committed 86/88 occupation table
onto the shared observer-86 counts. -/
theorem jointCount88_marginal (i : Fin 13) :
    ∑ j : Fin 14, jointCount88 (i, j) = count86 i := by
  unfold jointCount88 count86
  rw [Finset.card_filter, Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro t _
  by_cases h : obs86.path t = i
  · rw [if_pos h, Finset.sum_eq_single (obs88.path t)]
    · rw [if_pos]
      exact Prod.ext h rfl
    · intro j _ hj
      rw [if_neg]
      intro hEq
      exact hj (congrArg Prod.snd hEq).symm
    · intro habs
      exact absurd (Finset.mem_univ _) habs
  · rw [if_neg h]
    apply Finset.sum_eq_zero
    intro j _
    rw [if_neg]
    intro hEq
    exact h (congrArg Prod.fst hEq)

/-- The committed observer-86 counts total the 32 source steps. -/
theorem count86_total : ∑ i : Fin 13, count86 i = 32 := by decide

/-- The hinge datum is a normalized law: its entries sum to one. -/
theorem hingeFn86_sum_one : ∑ i : Fin 13, hingeFn86 i = 1 := by
  unfold hingeFn86
  rw [← Finset.sum_div]
  rw [show (∑ i : Fin 13, (count86 i : ℂ)) =
      ((∑ i : Fin 13, count86 i : ℕ) : ℂ) by push_cast; rfl]
  rw [count86_total]
  norm_num

/-- Marginalizing the 86/88 counted state recovers the committed hinge
datum. -/
theorem marg88_stateFn88 : marg88 stateFn88 = hingeFn86 := by
  funext i
  rw [marg88_apply]
  unfold stateFn88 hingeFn86
  rw [← Finset.sum_div]
  congr 1
  rw [show (∑ j : Fin 14, (jointCount88 (i, j) : ℂ)) =
      ((∑ j : Fin 14, jointCount88 (i, j) : ℕ) : ℂ) by push_cast; rfl]
  rw [jointCount88_marginal i]

/-- Marginalizing the 86/247 counted state recovers the committed hinge
datum. -/
theorem marg247_stateFn247 : marg247 stateFn247 = hingeFn86 := by
  funext i
  rw [marg247_apply]
  unfold stateFn247 hingeFn86
  rw [← Finset.sum_div]
  congr 1
  rw [show (∑ j : Fin 13, (jointCount (i, j) : ℂ)) =
      ((∑ j : Fin 13, jointCount (i, j) : ℕ) : ℂ) by push_cast; rfl]
  rw [jointCount_marginal_left i]

/-- **Marginal compatibility on the join.**  The two committed counted
states transport to join elements whose hinge restrictions coincide
exactly, both equal to the committed observer-86 marginal law. -/
theorem joinStates_marginal_compatible :
    hingeRestrict (embed88 stateFn88) = hingeFn86 ∧
    hingeRestrict (embed247 stateFn247) = hingeFn86 ∧
    hingeRestrict (embed88 stateFn88) =
      hingeRestrict (embed247 stateFn247) := by
  refine ⟨?_, ?_, ?_⟩
  · rw [hingeRestrict_embed88, marg88_stateFn88]
  · rw [hingeRestrict_embed247, marg247_stateFn247]
  · rw [hingeRestrict_embed88, marg88_stateFn88,
      hingeRestrict_embed247, marg247_stateFn247]

/-! ## Dimension and nondegeneracy of the join -/

theorem hingeGraph_generator_injective :
    Function.Injective (sect88.prod (-sect247)) := by
  intro a b hab
  apply sect88_injective
  have := congrArg Prod.fst hab
  simpa using this

/-- The gluing submodule has the hinge dimension 13. -/
theorem hingeGraph_finrank : Module.finrank ℂ hingeGraph = 13 := by
  show Module.finrank ℂ (LinearMap.range (sect88.prod (-sect247))) = 13
  rw [LinearMap.finrank_range_of_inj hingeGraph_generator_injective]
  simp [Module.finrank_fintype_fun_eq_card]

theorem carrierProd_finrank :
    Module.finrank ℂ (CarrierFn88 × CarrierFn247) = 351 := by
  rw [Module.finrank_prod]
  simp [Module.finrank_fintype_fun_eq_card]

/-- **The join dimension.**  The join carrier has dimension exactly
`182 + 169 - 13`: the two ambient dimensions minus the hinge
dimension. -/
theorem joinCarrier_finrank : Module.finrank ℂ JoinCarrier = 338 := by
  have h := Submodule.finrank_quotient_add_finrank hingeGraph
  rw [hingeGraph_finrank, carrierProd_finrank] at h
  have h2 : Module.finrank ℂ JoinCarrier + 13 = 351 := h
  omega

theorem joinCarrier_finrank_formula : (338 : ℕ) = 182 + 169 - 13 := by
  norm_num

/-- The join is strictly smaller than the hinge-ignoring direct sum. -/
theorem joinCarrier_lt_direct_sum :
    Module.finrank ℂ JoinCarrier <
      Module.finrank ℂ (CarrierFn88 × CarrierFn247) := by
  rw [joinCarrier_finrank, carrierProd_finrank]
  norm_num

/-- The two embedded images overlap nontrivially: the join is not a
direct sum with disjoint images. -/
theorem embed_ranges_overlap_ne_bot :
    LinearMap.range embed88 ⊓ LinearMap.range embed247 ≠ ⊥ := by
  rw [Submodule.ne_bot_iff]
  refine ⟨embed88 (sect88 fun _ => 1), ?_, ?_⟩
  · exact Submodule.mem_inf.mpr
      ⟨⟨sect88 fun _ => 1, rfl⟩,
        ⟨sect247 fun _ => 1, (embed_agree_on_hinge fun _ => 1).symm⟩⟩
  · intro hzero
    have h0 : embed88 (sect88 fun _ => 1) = embed88 0 := by
      rw [map_zero]; exact hzero
    have h1 : sect88 (fun _ => 1) = 0 := embed88_injective h0
    have h2 : (fun _ : Fin 13 => (1 : ℂ)) = 0 := by
      apply sect88_injective
      rw [map_zero]; exact h1
    exact one_ne_zero (congrFun h2 0)

/-- **The hinge image is exactly the intersection.**  A join element lies
in both embedded images precisely when it is an embedded hinge datum. -/
theorem embed_range_inter :
    LinearMap.range embed88 ⊓ LinearMap.range embed247 =
      LinearMap.range (embed88 ∘ₗ sect88) := by
  apply le_antisymm
  · intro z hz
    obtain ⟨⟨v, rfl⟩, ⟨w, hw⟩⟩ := Submodule.mem_inf.mp hz
    rw [embed88_apply, embed247_apply, Submodule.Quotient.eq] at hw
    obtain ⟨h, hh⟩ := hw
    rw [hingeGraph_gen_apply] at hh
    have hfst : sect88 h = -v := by
      have := congrArg Prod.fst hh
      simpa using this
    refine ⟨-h, ?_⟩
    rw [LinearMap.comp_apply, map_neg, hfst, neg_neg]
  · rintro z ⟨h, rfl⟩
    exact Submodule.mem_inf.mpr
      ⟨⟨sect88 h, rfl⟩, ⟨sect247 h, (embed_agree_on_hinge h).symm⟩⟩

/-- The hinge image inside the join has the hinge dimension 13. -/
theorem hingeImage_finrank :
    Module.finrank ℂ (LinearMap.range (embed88 ∘ₗ sect88)) = 13 := by
  have hinj : Function.Injective (embed88 ∘ₗ sect88) := fun a b hab =>
    sect88_injective (embed88_injective hab)
  rw [LinearMap.finrank_range_of_inj hinj]
  simp [Module.finrank_fintype_fun_eq_card]

/-! ## The committed anchoring maps at the diagonal layer -/

/-- The 86/88 anchoring at the function layer: precomposition with the
committed basis equivalence `pairFin`. -/
def anchorFn88 : CarrierFn88 ≃ₗ[ℂ] (Fin 182 → ℂ) where
  toFun v := v ∘ pairFin.symm
  invFun w := w ∘ pairFin
  map_add' _ _ := rfl
  map_smul' _ _ := rfl
  left_inv v := by
    funext p
    simp
  right_inv w := by
    funext i
    simp

/-- The 86/247 anchoring at the function layer: precomposition with the
committed basis equivalence `pairFin247`. -/
def anchorFn247 : CarrierFn247 ≃ₗ[ℂ] (Fin 169 → ℂ) where
  toFun v := v ∘ pairFin247.symm
  invFun w := w ∘ pairFin247
  map_add' _ _ := rfl
  map_smul' _ _ := rfl
  left_inv v := by
    funext p
    simp
  right_inv w := by
    funext i
    simp

/-- **Anchoring square, 86/88 side.**  On diagonal matrices, the committed
ambient anchoring star-equivalence acts exactly as the function-layer
anchoring. -/
theorem ambientEquiv_diagonal (v : CarrierFn88) :
    ambientEquiv (Matrix.diagonal v) = Matrix.diagonal (anchorFn88 v) := by
  rw [ambientEquiv_apply, Matrix.reindex_apply,
    Matrix.submatrix_diagonal_equiv]
  rfl

/-- **Anchoring square, 86/247 side.** -/
theorem ambientEquiv247_diagonal (v : CarrierFn247) :
    ambientEquiv247 (Matrix.diagonal v) =
      Matrix.diagonal (anchorFn247 v) := by
  rw [ambientEquiv247_apply, Matrix.reindex_apply,
    Matrix.submatrix_diagonal_equiv]
  rfl

/-- The anchored hinge sections. -/
def sectA88 : HingeFn →ₗ[ℂ] (Fin 182 → ℂ) :=
  anchorFn88.toLinearMap.comp sect88

def sectA247 : HingeFn →ₗ[ℂ] (Fin 169 → ℂ) :=
  anchorFn247.toLinearMap.comp sect247

/-- The anchored gluing submodule. -/
def hingeGraphA : Submodule ℂ ((Fin 182 → ℂ) × (Fin 169 → ℂ)) :=
  LinearMap.range (sectA88.prod (-sectA247))

/-- The anchored join carrier over the committed `Fin 182` and `Fin 169`
ambient index sets. -/
abbrev JoinCarrierA := ((Fin 182 → ℂ) × (Fin 169 → ℂ)) ⧸ hingeGraphA

def embedA88 : (Fin 182 → ℂ) →ₗ[ℂ] JoinCarrierA :=
  hingeGraphA.mkQ.comp (LinearMap.inl ℂ (Fin 182 → ℂ) (Fin 169 → ℂ))

def embedA247 : (Fin 169 → ℂ) →ₗ[ℂ] JoinCarrierA :=
  hingeGraphA.mkQ.comp (LinearMap.inr ℂ (Fin 182 → ℂ) (Fin 169 → ℂ))

theorem hingeGraph_map_anchor :
    hingeGraph.map
        ((anchorFn88.prodCongr anchorFn247 :
          (CarrierFn88 × CarrierFn247) ≃ₗ[ℂ]
            (Fin 182 → ℂ) × (Fin 169 → ℂ)) :
          (CarrierFn88 × CarrierFn247) →ₗ[ℂ]
            (Fin 182 → ℂ) × (Fin 169 → ℂ)) = hingeGraphA := by
  show Submodule.map _ (LinearMap.range (sect88.prod (-sect247))) = _
  rw [← LinearMap.range_comp]
  rfl

/-- The join transports along the committed anchoring equivalences. -/
def joinAnchor : JoinCarrier ≃ₗ[ℂ] JoinCarrierA :=
  Submodule.Quotient.equiv hingeGraph hingeGraphA
    (anchorFn88.prodCongr anchorFn247) hingeGraph_map_anchor

/-- The join anchoring on class representatives. -/
theorem joinAnchor_mk (x : CarrierFn88 × CarrierFn247) :
    joinAnchor (Submodule.Quotient.mk x) =
      Submodule.Quotient.mk (anchorFn88 x.1, anchorFn247 x.2) := by
  unfold joinAnchor
  rw [Submodule.Quotient.equiv_apply, Submodule.mapQ_apply]
  rfl

/-- **Anchoring commuting square on the 86/88 embedding.** -/
theorem joinAnchor_embed88 (v : CarrierFn88) :
    joinAnchor (embed88 v) = embedA88 (anchorFn88 v) := by
  rw [embed88_apply, joinAnchor_mk]
  show Submodule.Quotient.mk (anchorFn88 v, anchorFn247 0) =
    Submodule.Quotient.mk (anchorFn88 v, 0)
  rw [map_zero]

/-- **Anchoring commuting square on the 86/247 embedding.** -/
theorem joinAnchor_embed247 (w : CarrierFn247) :
    joinAnchor (embed247 w) = embedA247 (anchorFn247 w) := by
  rw [embed247_apply, joinAnchor_mk]
  show Submodule.Quotient.mk (anchorFn88 0, anchorFn247 w) =
    Submodule.Quotient.mk (0, anchorFn247 w)
  rw [map_zero]

/-! ## The committed walk-step transports at the diagonal layer

The transports run over the 31 actual adjacent source transitions
(`Fin 31`): no terminal-to-initial cyclic edge enters the join. -/

/-- The observer-86 hinge action of one committed walk step. -/
def hingeStep (t : Fin 31) : HingeFn →ₗ[ℂ] HingeFn where
  toFun h := h ∘ (walkSwapLeft t.castSucc).symm
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

theorem hingeStep_apply (t : Fin 31) (h : HingeFn) (i : Fin 13) :
    hingeStep t h i = h ((walkSwapLeft t.castSucc).symm i) := rfl

/-- The 86/88 diagonal-layer action of one committed walk step. -/
def evolveFn88 (t : Fin 31) : CarrierFn88 →ₗ[ℂ] CarrierFn88 where
  toFun v := v ∘ (walkStepEquiv88 t).symm
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

theorem evolveFn88_apply (t : Fin 31) (v : CarrierFn88) (p : PairIndex) :
    evolveFn88 t v p = v ((walkStepEquiv88 t).symm p) := rfl

/-- The 86/247 diagonal-layer action of one committed walk step. -/
def evolveFn247 (t : Fin 31) : CarrierFn247 →ₗ[ℂ] CarrierFn247 where
  toFun v := v ∘ (walkStepEquiv t.castSucc).symm
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

theorem evolveFn247_apply (t : Fin 31) (v : CarrierFn247)
    (p : PairIndex247) :
    evolveFn247 t v p = v ((walkStepEquiv t.castSucc).symm p) := rfl

theorem walkStepEquiv88_symm_apply (t : Fin 31) (p : PairIndex) :
    (walkStepEquiv88 t).symm p =
      ((walkSwapLeft t.castSucc).symm p.1, (walkSwap88 t).symm p.2) := by
  rcases p with ⟨a, b⟩
  rfl

theorem walkStepEquiv247_symm_apply (t : Fin 31) (p : PairIndex247) :
    (walkStepEquiv t.castSucc).symm p =
      ((walkSwapLeft t.castSucc).symm p.1,
        (walkSwapRight t.castSucc).symm p.2) := by
  rcases p with ⟨a, b⟩
  rfl

/-- **The committed 86/88 walk transport is the diagonal-layer action.**
The exact commuting square between the committed matrix-level evolution
and the function-level evolution. -/
theorem pair88StepEvolve_diagonal (t : Fin 31) (v : CarrierFn88) :
    pair88StepEvolve t (Matrix.diagonal v) =
      Matrix.diagonal (evolveFn88 t v) := by
  show Matrix.reindex (walkStepEquiv88 t) (walkStepEquiv88 t)
      (Matrix.diagonal v) = _
  rw [Matrix.reindex_apply, Matrix.submatrix_diagonal_equiv]
  rfl

/-- **The committed 86/247 walk transport is the diagonal-layer action.** -/
theorem stepEvolve_diagonal (t : Fin 31) (v : CarrierFn247) :
    stepEvolve t.castSucc (Matrix.diagonal v) =
      Matrix.diagonal (evolveFn247 t v) := by
  rw [stepEvolve_apply, Matrix.reindex_apply,
    Matrix.submatrix_diagonal_equiv]
  rfl

/-- The hinge sections intertwine the hinge action with each carrier
action: both sides evolve the shared observer-86 labels by the same
committed transposition. -/
theorem evolveFn88_sect88 (t : Fin 31) (h : HingeFn) :
    evolveFn88 t (sect88 h) = sect88 (hingeStep t h) := by
  funext p
  rw [evolveFn88_apply, walkStepEquiv88_symm_apply, sect88_apply,
    sect88_apply, hingeStep_apply]

theorem evolveFn247_sect247 (t : Fin 31) (h : HingeFn) :
    evolveFn247 t (sect247 h) = sect247 (hingeStep t h) := by
  funext p
  rw [evolveFn247_apply, walkStepEquiv247_symm_apply, sect247_apply,
    sect247_apply, hingeStep_apply]

/-- Marginalization intertwines each carrier action with the hinge
action. -/
theorem marg88_evolveFn88 (t : Fin 31) (v : CarrierFn88) :
    marg88 (evolveFn88 t v) = hingeStep t (marg88 v) := by
  funext i
  rw [marg88_apply, hingeStep_apply, marg88_apply]
  refine Fintype.sum_equiv (walkSwap88 t).symm _ _ ?_
  intro j
  rw [evolveFn88_apply, walkStepEquiv88_symm_apply]

theorem marg247_evolveFn247 (t : Fin 31) (v : CarrierFn247) :
    marg247 (evolveFn247 t v) = hingeStep t (marg247 v) := by
  funext i
  rw [marg247_apply, hingeStep_apply, marg247_apply]
  refine Fintype.sum_equiv (walkSwapRight t.castSucc).symm _ _ ?_
  intro j
  rw [evolveFn247_apply, walkStepEquiv247_symm_apply]

theorem hingeGraph_evolve_invariant (t : Fin 31) :
    hingeGraph ≤ hingeGraph.comap
      ((evolveFn88 t).prodMap (evolveFn247 t)) := by
  rintro x ⟨h, rfl⟩
  rw [Submodule.mem_comap, hingeGraph_gen_apply, LinearMap.prodMap_apply,
    map_neg, evolveFn88_sect88, evolveFn247_sect247]
  exact ⟨hingeStep t h, rfl⟩

/-- The join evolution: the descent of the two committed diagonal-layer
walk transports.  It exists exactly because both sides share the
observer-86 hinge action. -/
def joinStep (t : Fin 31) : JoinCarrier →ₗ[ℂ] JoinCarrier :=
  hingeGraph.mapQ hingeGraph ((evolveFn88 t).prodMap (evolveFn247 t))
    (hingeGraph_evolve_invariant t)

/-- **Transport commuting square on the 86/88 embedding.** -/
theorem joinStep_embed88 (t : Fin 31) (v : CarrierFn88) :
    joinStep t (embed88 v) = embed88 (evolveFn88 t v) := by
  rw [embed88_apply, embed88_apply]
  show hingeGraph.mapQ hingeGraph _ _ (Submodule.Quotient.mk (v, 0)) = _
  rw [Submodule.mapQ_apply]
  rfl

/-- **Transport commuting square on the 86/247 embedding.** -/
theorem joinStep_embed247 (t : Fin 31) (w : CarrierFn247) :
    joinStep t (embed247 w) = embed247 (evolveFn247 t w) := by
  rw [embed247_apply, embed247_apply]
  show hingeGraph.mapQ hingeGraph _ _ (Submodule.Quotient.mk (0, w)) = _
  rw [Submodule.mapQ_apply]
  rfl

/-- The hinge restriction intertwines the join evolution with the hinge
action. -/
theorem hingeRestrict_joinStep (t : Fin 31) (z : JoinCarrier) :
    hingeRestrict (joinStep t z) = hingeStep t (hingeRestrict z) := by
  obtain ⟨⟨v, w⟩, rfl⟩ := Submodule.Quotient.mk_surjective hingeGraph z
  show hingeRestrict (hingeGraph.mapQ hingeGraph _ _
      (Submodule.Quotient.mk (v, w))) = _
  rw [Submodule.mapQ_apply, hingeRestrict_mk, hingeRestrict_mk,
    LinearMap.prodMap_apply]
  show marg88 (evolveFn88 t v) + marg247 (evolveFn247 t w) = _
  rw [marg88_evolveFn88, marg247_evolveFn247, map_add]

/-! ## The join receipt -/

/-- **The typed carrier join, in one receipt.**  Both committed carriers
embed injectively into one explicit finite join; the embeddings agree
exactly on the observer-86 hinge; the join dimension is exactly the two
ambient dimensions minus the hinge dimension; the intersection of the
embedded images is exactly the hinge image; the transported committed
counted states restrict to the same committed hinge law; and both
embeddings intertwine the committed anchoring equivalences and the
committed walk-step transports through the commuting squares proved
above. -/
theorem carrierJoin_receipt :
    Function.Injective embed88 ∧
    Function.Injective embed247 ∧
    (∀ h : HingeFn, embed88 (sect88 h) = embed247 (sect247 h)) ∧
    Module.finrank ℂ JoinCarrier = 338 ∧
    (338 : ℕ) = 182 + 169 - 13 ∧
    Module.finrank ℂ JoinCarrier <
      Module.finrank ℂ (CarrierFn88 × CarrierFn247) ∧
    (LinearMap.range embed88 ⊓ LinearMap.range embed247 ≠ ⊥) ∧
    (LinearMap.range embed88 ⊓ LinearMap.range embed247 =
      LinearMap.range (embed88 ∘ₗ sect88)) ∧
    (hingeRestrict (embed88 stateFn88) = hingeFn86 ∧
      hingeRestrict (embed247 stateFn247) = hingeFn86) ∧
    (∀ v, joinAnchor (embed88 v) = embedA88 (anchorFn88 v)) ∧
    (∀ w, joinAnchor (embed247 w) = embedA247 (anchorFn247 w)) ∧
    (∀ (t : Fin 31) (v), joinStep t (embed88 v) =
      embed88 (evolveFn88 t v)) ∧
    (∀ (t : Fin 31) (w), joinStep t (embed247 w) =
      embed247 (evolveFn247 t w)) ∧
    (∀ (t : Fin 31) (z), hingeRestrict (joinStep t z) =
      hingeStep t (hingeRestrict z)) :=
  ⟨embed88_injective, embed247_injective, embed_agree_on_hinge,
    joinCarrier_finrank, joinCarrier_finrank_formula,
    joinCarrier_lt_direct_sum, embed_ranges_overlap_ne_bot,
    embed_range_inter,
    ⟨joinStates_marginal_compatible.1, joinStates_marginal_compatible.2.1⟩,
    joinAnchor_embed88, joinAnchor_embed247,
    joinStep_embed88, joinStep_embed247, hingeRestrict_joinStep⟩

/-! ## Evolved-data time slices: rigidity of the region family -/

/-- A matrix star-automorphism carries the scalars onto the scalars. -/
theorem starEquiv_map_bot {γ : Type*} [Fintype γ] [DecidableEq γ]
    (e : Matrix γ γ ℂ ≃⋆ₐ[ℂ] Matrix γ γ ℂ) :
    StarSubalgebra.map
        (e : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix γ γ ℂ)
        (⊥ : StarSubalgebra ℂ (Matrix γ γ ℂ)) = ⊥ := by
  apply le_antisymm
  · rintro x ⟨y, hy, rfl⟩
    obtain ⟨c, rfl⟩ := StarSubalgebra.mem_bot.mp hy
    exact ⟨c, (AlgHomClass.commutes
      (e : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix γ γ ℂ) c).symm⟩
  · intro x hx
    obtain ⟨c, rfl⟩ := StarSubalgebra.mem_bot.mp hx
    exact ⟨algebraMap ℂ _ c, ⟨c, rfl⟩,
      AlgHomClass.commutes (e : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix γ γ ℂ) c⟩

/-- A matrix star-automorphism maps a subalgebra onto the scalars only
if the subalgebra is the scalars. -/
theorem starEquiv_map_eq_bot_iff {γ : Type*} [Fintype γ] [DecidableEq γ]
    (e : Matrix γ γ ℂ ≃⋆ₐ[ℂ] Matrix γ γ ℂ)
    (A : StarSubalgebra ℂ (Matrix γ γ ℂ)) :
    StarSubalgebra.map (e : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix γ γ ℂ) A = ⊥ ↔
      A = ⊥ := by
  constructor
  · intro hmap
    apply le_antisymm _ bot_le
    intro a ha
    have hmem : e a ∈ StarSubalgebra.map
        (e : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix γ γ ℂ) A := ⟨a, ha, rfl⟩
    rw [hmap] at hmem
    obtain ⟨c, hc⟩ := StarSubalgebra.mem_bot.mp hmem
    have he : e (algebraMap ℂ (Matrix γ γ ℂ) c) =
        algebraMap ℂ (Matrix γ γ ℂ) c :=
      AlgHomClass.commutes (e : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix γ γ ℂ) c
    have hax : e a = e (algebraMap ℂ (Matrix γ γ ℂ) c) := by
      rw [he]
      exact hc.symm
    exact ⟨c, (e.injective hax).symm⟩
  · rintro rfl
    exact starEquiv_map_bot e

/-- A matrix star-automorphism carries the full algebra onto the full
algebra. -/
theorem starEquiv_map_top {γ : Type*} [Fintype γ] [DecidableEq γ]
    (e : Matrix γ γ ℂ ≃⋆ₐ[ℂ] Matrix γ γ ℂ) :
    StarSubalgebra.map
        (e : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix γ γ ℂ)
        (⊤ : StarSubalgebra ℂ (Matrix γ γ ℂ)) = ⊤ := by
  rw [← StarAlgHom.range_eq_map_top, StarSubalgebra.eq_top_iff]
  intro x
  exact ⟨e.symm x, e.apply_symm_apply x⟩

/-- A matrix star-automorphism maps a subalgebra onto the full algebra
only if the subalgebra is the full algebra. -/
theorem starEquiv_map_eq_top_iff {γ : Type*} [Fintype γ] [DecidableEq γ]
    (e : Matrix γ γ ℂ ≃⋆ₐ[ℂ] Matrix γ γ ℂ)
    (A : StarSubalgebra ℂ (Matrix γ γ ℂ)) :
    StarSubalgebra.map (e : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix γ γ ℂ) A = ⊤ ↔
      A = ⊤ := by
  constructor
  · intro hmap
    rw [StarSubalgebra.eq_top_iff]
    intro x
    have hmem : e x ∈ StarSubalgebra.map
        (e : Matrix γ γ ℂ →⋆ₐ[ℂ] Matrix γ γ ℂ) A := by
      rw [hmap]
      trivial
    obtain ⟨a, ha, hax⟩ := hmem
    have : a = x := e.injective hax
    rwa [← this]
  · rintro rfl
    exact starEquiv_map_top e

/-- **Region-preserving evolutions force step-invariant regions.**  If
every step automorphism of an inhabitant of the corrected interface over
the committed 86/247 diamond maps each regional algebra onto itself,
then the region family cannot depend on time along the step map. -/
theorem region_preserving_forces_slice_invariance
    (E : TimeIndexedNetEvolution PairIndex247 twoSlotNet247)
    (hpres : ∀ (t : E.Time) (U : TwoSlotRegion),
      StarSubalgebra.map
          (E.evolve t : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
            Matrix PairIndex247 PairIndex247 ℂ)
          (twoSlotAlgebra247 U) = twoSlotAlgebra247 U) :
    ∀ t : E.Time, E.slice (E.step t) = E.slice t := by
  intro t
  apply twoSlotAlgebra247_injective
  have h : StarSubalgebra.map
      (E.evolve t : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
        Matrix PairIndex247 PairIndex247 ℂ)
      (twoSlotAlgebra247 (E.slice t)) =
      twoSlotAlgebra247 (E.slice (E.step t)) := E.evolve_net_compatible t
  rw [hpres t (E.slice t)] at h
  exact h.symm

/-- **The committed walk evolution forces a constant region family.**
Any slice family over cyclic time `Fin 32` that is net-compatible with
the committed walk-step evolution is constant.  The constant-region
inhabitant of the corrected interface is therefore the only region shape
the committed walk evolution admits. -/
theorem stepEvolve_slice_constant
    (slice : Fin 32 → TwoSlotRegion)
    (hcompat : ∀ t : Fin 32,
      StarSubalgebra.map
          (stepEvolve t : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
            Matrix PairIndex247 PairIndex247 ℂ)
          (twoSlotAlgebra247 (slice t)) = twoSlotAlgebra247 (slice (t + 1))) :
    ∀ t : Fin 32, slice t = slice 0 := by
  have hstep : ∀ t : Fin 32, slice (t + 1) = slice t := by
    intro t
    apply twoSlotAlgebra247_injective
    rw [← hcompat t, stepEvolve_maps_regions t (slice t)]
  have hnat : ∀ n : ℕ, ∀ hn : n < 32, slice ⟨n, hn⟩ = slice 0 := by
    intro n
    induction n with
    | zero => intro hn; rfl
    | succ m ih =>
        intro hn
        have hm : m < 32 := Nat.lt_of_succ_lt hn
        have hadd : (⟨m, hm⟩ : Fin 32) + 1 = ⟨m + 1, hn⟩ := by
          apply Fin.ext
          simp [Fin.add_def]
          omega
        rw [← hadd, hstep, ih hm]
  intro t
  exact hnat t.val t.isLt

/-- The bottom region is rigid under any evolution family. -/
theorem timeSlice_bot_rigid
    (E : TimeIndexedNetEvolution PairIndex247 twoSlotNet247) (t : E.Time) :
    E.slice t = TwoSlotRegion.bot ↔
      E.slice (E.step t) = TwoSlotRegion.bot := by
  have h : StarSubalgebra.map
      (E.evolve t : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
        Matrix PairIndex247 PairIndex247 ℂ)
      (twoSlotAlgebra247 (E.slice t)) =
      twoSlotAlgebra247 (E.slice (E.step t)) := E.evolve_net_compatible t
  constructor
  · intro hb
    apply twoSlotAlgebra247_injective
    rw [← h, hb]
    exact starEquiv_map_bot (E.evolve t)
  · intro hb
    apply twoSlotAlgebra247_injective
    apply (starEquiv_map_eq_bot_iff (E.evolve t)
      (twoSlotAlgebra247 (E.slice t))).mp
    rw [h, hb]
    rfl

/-- The top region is rigid under any evolution family. -/
theorem timeSlice_top_rigid
    (E : TimeIndexedNetEvolution PairIndex247 twoSlotNet247) (t : E.Time) :
    E.slice t = TwoSlotRegion.top ↔
      E.slice (E.step t) = TwoSlotRegion.top := by
  have h : StarSubalgebra.map
      (E.evolve t : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
        Matrix PairIndex247 PairIndex247 ℂ)
      (twoSlotAlgebra247 (E.slice t)) =
      twoSlotAlgebra247 (E.slice (E.step t)) := E.evolve_net_compatible t
  constructor
  · intro hb
    apply twoSlotAlgebra247_injective
    rw [← h, hb]
    exact starEquiv_map_top (E.evolve t)
  · intro hb
    apply twoSlotAlgebra247_injective
    apply (starEquiv_map_eq_top_iff (E.evolve t)
      (twoSlotAlgebra247 (E.slice t))).mp
    rw [h, hb]
    rfl

/-- **A time-dependent region family can only exchange the two slots.**
If the region changes across one step of any inhabitant of the corrected
interface, the change is from one slot region to the other. -/
theorem slice_change_forces_slot_exchange
    (E : TimeIndexedNetEvolution PairIndex247 twoSlotNet247) (t : E.Time)
    (hne : E.slice (E.step t) ≠ E.slice t) :
    (E.slice t = TwoSlotRegion.left ∧
      E.slice (E.step t) = TwoSlotRegion.right) ∨
    (E.slice t = TwoSlotRegion.right ∧
      E.slice (E.step t) = TwoSlotRegion.left) := by
  have hU : ∀ x : TwoSlotRegion,
      x = TwoSlotRegion.bot ∨ x = TwoSlotRegion.left ∨
      x = TwoSlotRegion.right ∨ x = TwoSlotRegion.top := by
    intro x
    cases x with
    | bot => exact Or.inl rfl
    | left => exact Or.inr (Or.inl rfl)
    | right => exact Or.inr (Or.inr (Or.inl rfl))
    | top => exact Or.inr (Or.inr (Or.inr rfl))
  rcases hU (E.slice t) with hb | hl | hr | ht
  · exact absurd (((timeSlice_bot_rigid E t).mp hb).trans hb.symm) hne
  · rcases hU (E.slice (E.step t)) with hb' | hl' | hr' | ht'
    · have hx : TwoSlotRegion.left = TwoSlotRegion.bot :=
        hl.symm.trans ((timeSlice_bot_rigid E t).mpr hb')
      exact TwoSlotRegion.noConfusion hx
    · exact absurd (hl'.trans hl.symm) hne
    · exact Or.inl ⟨hl, hr'⟩
    · have hx : TwoSlotRegion.left = TwoSlotRegion.top :=
        hl.symm.trans ((timeSlice_top_rigid E t).mpr ht')
      exact TwoSlotRegion.noConfusion hx
  · rcases hU (E.slice (E.step t)) with hb' | hl' | hr' | ht'
    · have hx : TwoSlotRegion.right = TwoSlotRegion.bot :=
        hr.symm.trans ((timeSlice_bot_rigid E t).mpr hb')
      exact TwoSlotRegion.noConfusion hx
    · exact Or.inr ⟨hr, hl'⟩
    · exact absurd (hr'.trans hr.symm) hne
    · have hx : TwoSlotRegion.right = TwoSlotRegion.top :=
        hr.symm.trans ((timeSlice_top_rigid E t).mpr ht')
      exact TwoSlotRegion.noConfusion hx
  · exact absurd (((timeSlice_top_rigid E t).mp ht).trans ht.symm) hne

/-! ## Evolved-data time slices: the slot-alternating inhabitant -/

/-- The alternating region family: the left slot at even steps, the
right slot at odd steps. -/
def alternatingSlice (t : Fin 32) : TwoSlotRegion :=
  if t.val % 2 = 0 then TwoSlotRegion.left else TwoSlotRegion.right

/-- One step of cyclic time exchanges the two slots of the alternating
family; the even length 32 closes the cycle consistently. -/
theorem alternatingSlice_step : ∀ t : Fin 32,
    alternatingSlice (t + 1) = (alternatingSlice t).swap := by
  decide

/-- The alternating family is genuinely time-dependent. -/
theorem alternatingSlice_not_constant :
    alternatingSlice 1 ≠ alternatingSlice 0 := by
  decide

/-- The slot-alternating evolution: the committed walk-step reindexing
followed by the committed internal slot exchange.  The walk factor
carries the committed transitions (including the declared index-31
bookkeeping closure); the exchange factor is a constructed internal
algebra symmetry, not a source transition. -/
def alternatingEvolve (t : Fin 32) :
    Matrix PairIndex247 PairIndex247 ℂ ≃⋆ₐ[ℂ]
      Matrix PairIndex247 PairIndex247 ℂ :=
  (stepEvolve t).trans slotSwap

theorem alternatingEvolve_toHom (t : Fin 32) :
    (alternatingEvolve t : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
        Matrix PairIndex247 PairIndex247 ℂ) =
      (slotSwap : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
          Matrix PairIndex247 PairIndex247 ℂ).comp
        (stepEvolve t : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
          Matrix PairIndex247 PairIndex247 ℂ) :=
  DFunLike.ext _ _ fun _ => rfl

/-- The slot-alternating evolution maps each regional algebra of the
committed diamond onto the algebra of the swapped region. -/
theorem alternatingEvolve_map (t : Fin 32) (U : TwoSlotRegion) :
    StarSubalgebra.map
        (alternatingEvolve t : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
          Matrix PairIndex247 PairIndex247 ℂ)
        (twoSlotAlgebra247 U) = twoSlotAlgebra247 U.swap := by
  rw [alternatingEvolve_toHom, ← StarSubalgebra.map_map,
    stepEvolve_maps_regions t U, slotSwap_map_localAlgebra U]

/-- Observer 247's committed source generators, lifted into the right
slot of the committed carrier. -/
def liftedRightSourceGenerators : Set (Matrix PairIndex247 PairIndex247 ℂ) :=
  (slotRight (Fin 13) (β := Fin 13) : Matrix (Fin 13) (Fin 13) ℂ →⋆ₐ[ℂ]
    Matrix PairIndex247 PairIndex247 ℂ) '' obs247.generators

/-- The lifted right source generators generate exactly the right
regional algebra, by the committed source-generation theorem. -/
theorem liftedRightSourceGenerators_generate_right :
    StarAlgebra.adjoin ℂ liftedRightSourceGenerators =
      (slotRight (Fin 13) (β := Fin 13)).range := by
  unfold liftedRightSourceGenerators
  rw [← StarAlgHom.map_adjoin]
  exact obs247_lifted_eq_rightSlot247

/-- The alternating admitted slice data: the committed observer-86
source generators in the left slot at even steps, the committed
observer-247 source generators in the right slot at odd steps. -/
def alternatingAdmitted (t : Fin 32) :
    Set (Matrix PairIndex247 PairIndex247 ℂ) :=
  if t.val % 2 = 0 then liftedSourceGenerators
  else liftedRightSourceGenerators

/-- **Transport-aware generation of the alternating candidate.**  At
every step, evolving the admitted slice data first and generating then
yields exactly the next (exchanged) slot algebra. -/
theorem alternatingAdmitted_generates (t : Fin 32) :
    StarAlgebra.adjoin ℂ
        ((alternatingEvolve t : Matrix PairIndex247 PairIndex247 ℂ →
            Matrix PairIndex247 PairIndex247 ℂ) '' alternatingAdmitted t) =
      twoSlotAlgebra247 (alternatingSlice (t + 1)) := by
  change StarAlgebra.adjoin ℂ
      ((alternatingEvolve t : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
          Matrix PairIndex247 PairIndex247 ℂ) '' alternatingAdmitted t) = _
  rw [← StarAlgHom.map_adjoin, alternatingSlice_step]
  unfold alternatingAdmitted
  by_cases h : t.val % 2 = 0
  · rw [if_pos h, liftedSourceGenerators_generate_left]
    have hs : alternatingSlice t = TwoSlotRegion.left := by
      unfold alternatingSlice
      rw [if_pos h]
    rw [hs]
    exact alternatingEvolve_map t TwoSlotRegion.left
  · rw [if_neg h, liftedRightSourceGenerators_generate_right]
    have hs : alternatingSlice t = TwoSlotRegion.right := by
      unfold alternatingSlice
      rw [if_neg h]
    rw [hs]
    exact alternatingEvolve_map t TwoSlotRegion.right

/-- **The slot-alternating inhabitant.**  The corrected one-step
interface is inhabited over the committed 86/247 diamond with a
genuinely time-dependent region family: cyclic time `Fin 32`, the
alternating left/right slice, the walk-step-plus-exchange evolution, and
the committed per-observer source generators as admitted slice data.
Generation is proved only after applying each one-step automorphism.
This proves the possibility the withdrawn forced-step-invariance no-go
left open; the delimitation theorems above prove its exact boundary. -/
def slotAlternatingEvolution :
    TimeIndexedNetEvolution PairIndex247 twoSlotNet247 where
  Time := Fin 32
  timeFintype := inferInstance
  step := fun t => t + 1
  slice := fun t => (show twoSlotNet247.Region from alternatingSlice t)
  admitted := alternatingAdmitted
  evolve := alternatingEvolve
  evolve_net_compatible := fun t => by
    show StarSubalgebra.map _ (twoSlotAlgebra247 (alternatingSlice t)) =
      twoSlotAlgebra247 (alternatingSlice (t + 1))
    rw [alternatingSlice_step]
    exact alternatingEvolve_map t (alternatingSlice t)
  admitted_mem := fun t => by
    show alternatingAdmitted t ⊆
      (twoSlotAlgebra247 (alternatingSlice t) :
        Set (Matrix PairIndex247 PairIndex247 ℂ))
    unfold alternatingAdmitted alternatingSlice
    by_cases h : t.val % 2 = 0
    · rw [if_pos h, if_pos h]
      rintro X ⟨A, _, rfl⟩
      exact ⟨A, rfl⟩
    · rw [if_neg h, if_neg h]
      rintro X ⟨B, _, rfl⟩
      exact ⟨B, rfl⟩
  slice_generates := fun t => by
    show StarAlgebra.adjoin ℂ _ =
      twoSlotAlgebra247 (alternatingSlice (t + 1))
    exact alternatingAdmitted_generates t

/-- The slot-alternating inhabitant starts on the left slot region. -/
theorem slotAlternatingEvolution_slice_zero :
    slotAlternatingEvolution.slice (0 : Fin 32) = TwoSlotRegion.left := rfl

/-- The slot-alternating inhabitant is on the right slot region at time
one. -/
theorem slotAlternatingEvolution_slice_one :
    slotAlternatingEvolution.slice (1 : Fin 32) = TwoSlotRegion.right := rfl

/-- **Genuine time dependence.**  The region family of the
slot-alternating inhabitant is not step-invariant. -/
theorem slotAlternatingEvolution_not_step_invariant :
    slotAlternatingEvolution.slice
        (slotAlternatingEvolution.step (0 : Fin 32)) ≠
      slotAlternatingEvolution.slice (0 : Fin 32) := by
  intro h
  have h' : TwoSlotRegion.right = TwoSlotRegion.left := h
  exact TwoSlotRegion.noConfusion h'

/-- **The evolved-data time-slice result, in one receipt.**  The
corrected interface admits a genuinely time-dependent region family
(the slot-alternating inhabitant); every net-compatible slice family
over the committed walk evolution is constant; and any one-step region
change of any inhabitant is an exchange of the two slot regions. -/
theorem evolvedTimeSlice_receipt :
    (slotAlternatingEvolution.slice
        (slotAlternatingEvolution.step (0 : Fin 32)) ≠
      slotAlternatingEvolution.slice (0 : Fin 32)) ∧
    (∀ (slice : Fin 32 → TwoSlotRegion),
      (∀ t : Fin 32,
        StarSubalgebra.map
            (stepEvolve t : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
              Matrix PairIndex247 PairIndex247 ℂ)
            (twoSlotAlgebra247 (slice t)) =
          twoSlotAlgebra247 (slice (t + 1))) →
      ∀ t : Fin 32, slice t = slice 0) ∧
    (∀ (E : TimeIndexedNetEvolution PairIndex247 twoSlotNet247)
        (t : E.Time), E.slice (E.step t) ≠ E.slice t →
      (E.slice t = TwoSlotRegion.left ∧
        E.slice (E.step t) = TwoSlotRegion.right) ∨
      (E.slice t = TwoSlotRegion.right ∧
        E.slice (E.step t) = TwoSlotRegion.left)) :=
  ⟨slotAlternatingEvolution_not_step_invariant,
    stepEvolve_slice_constant,
    slice_change_forces_slot_exchange⟩

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
-- Expected axioms per declaration: propext, Classical.choice, Quot.sound.
#print axioms OPH.QFT.marg88
#print axioms OPH.QFT.marg247
#print axioms OPH.QFT.sect88
#print axioms OPH.QFT.sect247
#print axioms OPH.QFT.marg88_sect88
#print axioms OPH.QFT.marg247_sect247
#print axioms OPH.QFT.embed_agree_on_hinge
#print axioms OPH.QFT.embed88_injective
#print axioms OPH.QFT.embed247_injective
#print axioms OPH.QFT.hingeRestrict
#print axioms OPH.QFT.hingeRestrict_embed88
#print axioms OPH.QFT.hingeRestrict_embed247
#print axioms OPH.QFT.hingeRestrict_retraction
#print axioms OPH.QFT.jointCount88_marginal
#print axioms OPH.QFT.count86_total
#print axioms OPH.QFT.hingeFn86_sum_one
#print axioms OPH.QFT.marg88_stateFn88
#print axioms OPH.QFT.marg247_stateFn247
#print axioms OPH.QFT.joinStates_marginal_compatible
#print axioms OPH.QFT.hingeGraph_finrank
#print axioms OPH.QFT.joinCarrier_finrank
#print axioms OPH.QFT.joinCarrier_lt_direct_sum
#print axioms OPH.QFT.embed_ranges_overlap_ne_bot
#print axioms OPH.QFT.embed_range_inter
#print axioms OPH.QFT.hingeImage_finrank
#print axioms OPH.QFT.ambientEquiv_diagonal
#print axioms OPH.QFT.ambientEquiv247_diagonal
#print axioms OPH.QFT.joinAnchor
#print axioms OPH.QFT.joinAnchor_embed88
#print axioms OPH.QFT.joinAnchor_embed247
#print axioms OPH.QFT.pair88StepEvolve_diagonal
#print axioms OPH.QFT.stepEvolve_diagonal
#print axioms OPH.QFT.evolveFn88_sect88
#print axioms OPH.QFT.evolveFn247_sect247
#print axioms OPH.QFT.marg88_evolveFn88
#print axioms OPH.QFT.marg247_evolveFn247
#print axioms OPH.QFT.joinStep
#print axioms OPH.QFT.joinStep_embed88
#print axioms OPH.QFT.joinStep_embed247
#print axioms OPH.QFT.hingeRestrict_joinStep
#print axioms OPH.QFT.carrierJoin_receipt
#print axioms OPH.QFT.starEquiv_map_bot
#print axioms OPH.QFT.starEquiv_map_eq_bot_iff
#print axioms OPH.QFT.starEquiv_map_top
#print axioms OPH.QFT.starEquiv_map_eq_top_iff
#print axioms OPH.QFT.region_preserving_forces_slice_invariance
#print axioms OPH.QFT.stepEvolve_slice_constant
#print axioms OPH.QFT.timeSlice_bot_rigid
#print axioms OPH.QFT.timeSlice_top_rigid
#print axioms OPH.QFT.slice_change_forces_slot_exchange
#print axioms OPH.QFT.alternatingSlice_step
#print axioms OPH.QFT.alternatingEvolve_map
#print axioms OPH.QFT.liftedRightSourceGenerators_generate_right
#print axioms OPH.QFT.alternatingAdmitted_generates
#print axioms OPH.QFT.slotAlternatingEvolution
#print axioms OPH.QFT.slotAlternatingEvolution_slice_zero
#print axioms OPH.QFT.slotAlternatingEvolution_slice_one
#print axioms OPH.QFT.slotAlternatingEvolution_not_step_invariant
#print axioms OPH.QFT.evolvedTimeSlice_receipt
