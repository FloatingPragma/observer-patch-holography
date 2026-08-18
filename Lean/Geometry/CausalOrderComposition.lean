import Geometry.MassShellKinematics
import QFT.FiniteCausalObserverNet
import QFT.RichFibreRegionalNet

/-!
# OL-A3: causal-order composition of committed finite nets with the declared Lorentz cone

The declared Hermitian Lorentz module of `Geometry/CanonicalLorentzModule.lean`
carries the relation `causalLE`: the coordinate difference has nonnegative
Lorentz square and nonnegative scalar coordinate.  This file proves that
relation is a partial order on the module and composes the committed finite
causal orders with it.

Main composition: the region lattice of the committed rich-fibre observer net
(`QFT/RichFibreRegionalNet.lean`, the second preregistered E1 run) embeds in
the cone order through the declared vertex assignment `richVertex`.  The
theorem `richRegion_le_iff_causalLE` proves the equivalence in both
directions on the committed six-element carrier: one region is below another
in the committed net exactly when the interval between their images is
future-causal in the declared module.  The embedding is injective, the image
contains two timelike-related pairs and two spacelike-unrelated pairs
(`richVertex_nondegenerate`), declared net disjointness lands in spacelike
separation (`richDisj_spacelike`), and composing with every oriented Lorentz
map, and with every affine overlap-cocycle transition, preserves the
equivalence.  The committed partition net
(`QFT/FiniteCausalObserverNet.lean`, `partitionPublicCausalNet`) composes
through the label map `siteRegionLabel` with the same receipts; its
four-element carrier affords exactly one incomparable doubleton
(`partitionNet_unique_incomparable_pair`), so one unordered spacelike pair is
the maximum any faithful embedding of that carrier can produce.

Negative controls: projecting the vertex assignment to the scalar axis
collapses the committed order to a chain and breaks the equivalence
(`timeProjection_not_faithful`), and the Lorentz square alone, with the
orientation clause dropped, fails antisymmetry (`q_alone_not_antisymm`).

**Declared data.**  The Lorentz module is declared.  The vertex assignments
`richVertex` and `siteRegionLabel` are declared witness constructions: the
composed receipt exhibits a second distinct assignment with the same
equivalence, so the committed finite data does not force the embedding.  The
regions and the disjointness relations are declared fields of the committed
nets, backed in the rich-fibre case by the preregistered source-disjointness
gate on the window literals.

**Falsifier.**  The composition fails if either direction of the order
equivalence fails at any pair of committed regions, if the vertex assignment
identifies two distinct committed regions, if a declared-disjoint pair of
regions has causally related images, if some oriented Lorentz map or affine
cocycle transition breaks the equivalence, or if the cone relation violates a
partial-order law.

**Nonclaims.**  No physical causal order, spacetime, clock, scheduler, or
observer-frame attachment is claimed.  The observer-to-physical-spacetime
and causal attachment premise (PR-52), the stable-causality and open-image
conditions (PR-16), and the finite subsystem-split interface (PR-44) stay
open.  The composition is a typed identification of committed finite orders
with one declared cone order through declared vertex maps; it does not
derive the module, the maps, or any continuum limit from observer data, and
the observation-ledger row OL-A3 stays a structural partial with this
composed receipt attached.
-/

namespace OPH.CausalComposition

noncomputable section

open OPH.C1Lorentz OPH.C2Soldering OPH.QFT

/-! ## The causal order of the declared Lorentz module -/

/-- Future-causal vectors of the declared `(+---)` module: nonnegative
Lorentz square and nonnegative scalar coordinate. -/
def IsFutureCausal (v : Herm2) : Prop := 0 ≤ lorentzQ v ∧ 0 ≤ v.1

/-- The causal relation of the declared module: the difference is
future-causal. -/
def causalLE (u v : Herm2) : Prop := IsFutureCausal (v - u)

/-- Strict timelike separation. -/
def IsTimelikeSep (u v : Herm2) : Prop := 0 < lorentzQ (v - u)

/-- Spacelike separation. -/
def IsSpacelikeSep (u v : Herm2) : Prop := lorentzQ (v - u) < 0

theorem spatialNormSq_nonneg (x : Spatial) : 0 ≤ spatialNormSq x :=
  Finset.sum_nonneg fun i _ => sq_nonneg (x i)

theorem lorentzQ_zero : lorentzQ (0 : Herm2) = 0 := by
  simp [lorentzQ, spatialNormSq]

theorem isFutureCausal_zero : IsFutureCausal 0 :=
  ⟨lorentzQ_zero.ge, le_rfl⟩

/-- The Cauchy-Schwarz inequality on the three Pauli coordinates, through
the Lagrange identity. -/
theorem spatialDot_sq_le_normSq_mul (x y : Spatial) :
    spatialDot x y ^ 2 ≤ spatialNormSq x * spatialNormSq y := by
  simp only [spatialDot, spatialNormSq, Fin.sum_univ_three]
  nlinarith [sq_nonneg (x 0 * y 1 - x 1 * y 0), sq_nonneg (x 0 * y 2 - x 2 * y 0),
    sq_nonneg (x 1 * y 2 - x 2 * y 1)]

/-- On two future-causal vectors the spatial pairing is dominated by the
product of the scalar coordinates. -/
theorem spatialDot_le_time_mul {a b : Herm2}
    (ha : IsFutureCausal a) (hb : IsFutureCausal b) :
    spatialDot a.2 b.2 ≤ a.1 * b.1 := by
  obtain ⟨hqa, hta⟩ := ha
  obtain ⟨hqb, htb⟩ := hb
  have hsa : spatialNormSq a.2 ≤ a.1 ^ 2 := by
    simp only [lorentzQ] at hqa
    linarith
  have hsb : spatialNormSq b.2 ≤ b.1 ^ 2 := by
    simp only [lorentzQ] at hqb
    linarith
  have hcs := spatialDot_sq_le_normSq_mul a.2 b.2
  have hsq : spatialDot a.2 b.2 ^ 2 ≤ (a.1 * b.1) ^ 2 := by
    have hna := spatialNormSq_nonneg a.2
    have hnb := spatialNormSq_nonneg b.2
    nlinarith
  by_contra hlt
  have hlt' : a.1 * b.1 < spatialDot a.2 b.2 := not_le.mp hlt
  have hab : 0 ≤ a.1 * b.1 := mul_nonneg hta htb
  have hpos : 0 < spatialDot a.2 b.2 - a.1 * b.1 := by linarith
  have hpos2 : 0 < spatialDot a.2 b.2 + a.1 * b.1 := by linarith
  nlinarith [mul_pos hpos hpos2]

theorem spatialNormSq_add (x y : Spatial) :
    spatialNormSq (x + y) =
      spatialNormSq x + 2 * spatialDot x y + spatialNormSq y := by
  simp only [spatialNormSq, spatialDot, Fin.sum_univ_three, Pi.add_apply]
  ring

/-- The future cone is closed under addition: the reverse triangle
inequality of the declared module. -/
theorem IsFutureCausal.add {a b : Herm2}
    (ha : IsFutureCausal a) (hb : IsFutureCausal b) :
    IsFutureCausal (a + b) := by
  have hdot := spatialDot_le_time_mul ha hb
  obtain ⟨hqa, hta⟩ := ha
  obtain ⟨hqb, htb⟩ := hb
  refine ⟨?_, ?_⟩
  · have hq : lorentzQ (a + b) =
        lorentzQ a + lorentzQ b + 2 * (a.1 * b.1 - spatialDot a.2 b.2) := by
      simp only [lorentzQ, Prod.fst_add, Prod.snd_add, spatialNormSq_add]
      ring
    rw [hq]
    linarith
  · simp only [Prod.fst_add]
    linarith

theorem causalLE_refl (v : Herm2) : causalLE v v := by
  simpa [causalLE, sub_self] using isFutureCausal_zero

theorem causalLE_trans {u v w : Herm2}
    (h1 : causalLE u v) (h2 : causalLE v w) : causalLE u w := by
  simpa [sub_add_sub_cancel] using h2.add h1

/-- A future-causal vector with vanishing scalar coordinate is zero. -/
theorem futureCausal_eq_zero_of_time_eq_zero {v : Herm2}
    (h : IsFutureCausal v) (ht : v.1 = 0) : v = 0 := by
  obtain ⟨hq, _⟩ := h
  simp only [lorentzQ, ht] at hq
  have hs0 : spatialNormSq v.2 = 0 := by
    have := spatialNormSq_nonneg v.2
    nlinarith
  have hv2 : v.2 = 0 := by
    by_contra hne
    exact absurd hs0 (ne_of_gt (spatialNormSq_pos hne))
  exact Prod.ext_iff.mpr ⟨ht, hv2⟩

theorem causalLE_antisymm {u v : Herm2}
    (h1 : causalLE u v) (h2 : causalLE v u) : u = v := by
  have ht1 : 0 ≤ (v - u).1 := h1.2
  have ht2 : 0 ≤ (u - v).1 := h2.2
  rw [Prod.fst_sub] at ht1 ht2
  have hz : (v - u).1 = 0 := by
    rw [Prod.fst_sub]
    linarith
  have hzero : v - u = 0 := futureCausal_eq_zero_of_time_eq_zero h1 hz
  exact (sub_eq_zero.mp hzero).symm

/-- The causal relation of the declared module is a partial order. -/
theorem causalLE_isPartialOrder :
    (∀ v, causalLE v v) ∧
    (∀ u v w, causalLE u v → causalLE v w → causalLE u w) ∧
    (∀ u v, causalLE u v → causalLE v u → u = v) :=
  ⟨causalLE_refl, fun _ _ _ h1 h2 => causalLE_trans h1 h2,
    fun _ _ h1 h2 => causalLE_antisymm h1 h2⟩

theorem lorentzQ_neg (v : Herm2) : lorentzQ (-v) = lorentzQ v := by
  have h := lorentzQ_smul (-1) v
  simpa [neg_one_smul] using h

theorem IsSpacelikeSep.symm {u v : Herm2}
    (h : IsSpacelikeSep u v) : IsSpacelikeSep v u := by
  unfold IsSpacelikeSep at h ⊢
  rw [← neg_sub v u, lorentzQ_neg]
  exact h

theorem not_causalLE_of_spacelike {u v : Herm2}
    (h : IsSpacelikeSep u v) : ¬ causalLE u v :=
  fun hc => absurd hc.1 (not_le.mpr h)

/-- Spacelike-separated points are causally unrelated in both directions. -/
theorem spacelike_not_related {u v : Herm2} (h : IsSpacelikeSep u v) :
    ¬ causalLE u v ∧ ¬ causalLE v u :=
  ⟨not_causalLE_of_spacelike h, not_causalLE_of_spacelike h.symm⟩

/-! ## Oriented-chart invariance of the cone order -/

/-- The inverse of an oriented Lorentz map, oriented through the committed
equivalence clauses. -/
def orientedSymm (L : OrientedLorentzEquiv) : OrientedLorentzEquiv where
  toLinearEquiv := L.toLinearEquiv.symm
  map_lorentzB := by
    intro u v
    conv_rhs => rw [← L.toLinearEquiv.apply_symm_apply u,
      ← L.toLinearEquiv.apply_symm_apply v]
    exact (L.map_lorentzB _ _).symm
  map_futureNull_iff := by
    intro v
    conv_rhs => rw [← L.toLinearEquiv.apply_symm_apply v]
    exact (L.map_futureNull_iff _).symm
  map_futureUnit_iff := by
    intro v
    conv_rhs => rw [← L.toLinearEquiv.apply_symm_apply v]
    exact (L.map_futureUnit_iff _).symm

/-- Every oriented Lorentz map preserves the future causal cone. -/
theorem futureCausal_map (L : OrientedLorentzEquiv) {v : Herm2}
    (h : IsFutureCausal v) : IsFutureCausal (L v) := by
  obtain ⟨hq, ht⟩ := h
  rcases eq_or_ne v 0 with rfl | hv
  · have h0 : L (0 : Herm2) = 0 := L.toLinearEquiv.map_zero
    rw [h0]
    exact isFutureCausal_zero
  · have htpos : 0 < v.1 := by
      rcases lt_or_eq_of_le ht with h' | h'
      · exact h'
      · exact absurd (futureCausal_eq_zero_of_time_eq_zero ⟨hq, ht⟩ h'.symm) hv
    rcases lt_or_eq_of_le hq with hqpos | hqzero
    · -- timelike branch: normalize to the frame hyperboloid
      set c := Real.sqrt (lorentzQ v) with hc
      have hcpos : 0 < c := Real.sqrt_pos.mpr hqpos
      have hcne : c ≠ 0 := ne_of_gt hcpos
      have hcsq : c ^ 2 = lorentzQ v := Real.sq_sqrt (le_of_lt hqpos)
      have hunit : IsFutureUnitTimelike (c⁻¹ • v) := by
        constructor
        · rw [lorentzQ_smul, ← hcsq, inv_pow,
            inv_mul_cancel₀ (pow_ne_zero 2 hcne)]
        · have hfst : (c⁻¹ • v).1 = c⁻¹ * v.1 := rfl
          rw [hfst]
          exact mul_pos (inv_pos.mpr hcpos) htpos
      have hLu := (L.map_futureUnit_iff (c⁻¹ • v)).mpr hunit
      rw [OrientedLorentzEquiv.map_smul] at hLu
      have htu : 0 < c⁻¹ * (L v).1 := hLu.2
      have hcinv : 0 < c⁻¹ := inv_pos.mpr hcpos
      have htLv : 0 ≤ (L v).1 := by nlinarith
      exact ⟨by rw [lorentzQ_oriented]; exact hq, htLv⟩
    · -- null branch
      have hnull : IsFutureNull v := ⟨hqzero.symm, htpos⟩
      have hLn := (L.map_futureNull_iff v).mpr hnull
      exact ⟨by rw [lorentzQ_oriented]; exact hq, le_of_lt hLn.2⟩

/-- Oriented Lorentz maps preserve the future cone in both directions. -/
theorem map_futureCausal_iff (L : OrientedLorentzEquiv) (v : Herm2) :
    IsFutureCausal (L v) ↔ IsFutureCausal v := by
  constructor
  · intro h
    have h2 := futureCausal_map (orientedSymm L) h
    rwa [show (orientedSymm L) (L v) = v from
      L.toLinearEquiv.symm_apply_apply v] at h2
  · exact futureCausal_map L

/-- Chart invariance of the cone order under every oriented Lorentz map. -/
theorem causalLE_map_iff (L : OrientedLorentzEquiv) (u v : Herm2) :
    causalLE (L u) (L v) ↔ causalLE u v := by
  unfold causalLE
  rw [← OrientedLorentzEquiv.map_sub L v u]
  exact map_futureCausal_iff L (v - u)

/-- Translation invariance of the cone order. -/
theorem causalLE_translate (w u v : Herm2) :
    causalLE (u + w) (v + w) ↔ causalLE u v := by
  unfold causalLE
  rw [add_sub_add_right_eq_sub]

/-- Chart invariance under every affine overlap-cocycle transition. -/
theorem causalLE_act_iff {Chart : Type*} (T : LorentzOverlapCocycle Chart)
    (i j : Chart) (u v : Herm2) :
    causalLE (T.act i j u) (T.act i j v) ↔ causalLE u v := by
  unfold causalLE
  rw [T.act_sub_act]
  exact map_futureCausal_iff (T.lorentz i j) (v - u)

/-! ## Cone embeddings of a finite declared order -/

/-- A typed composition of one finite declared order-with-disjointness into
the causal order of the declared Lorentz module: the order equivalence holds
in both directions and declared disjointness lands in spacelike
separation. -/
structure ConeEmbedding {P : Type*} (ord disj : P → P → Prop) where
  /-- The declared vertex assignment: each committed node maps to one point
  of the declared module. -/
  vertex : P → Herm2
  /-- The two-directional order equivalence on the committed carrier. -/
  le_iff : ∀ U V, ord U V ↔ causalLE (vertex U) (vertex V)
  /-- Declared disjointness maps to spacelike separation. -/
  disjoint_spacelike : ∀ U V, disj U V → IsSpacelikeSep (vertex U) (vertex V)

namespace ConeEmbedding

variable {P : Type*} {ord disj : P → P → Prop}

/-- Antisymmetry of the committed order makes every cone embedding
injective: the receipt that no two committed regions collapse to one module
point. -/
theorem vertex_injective (E : ConeEmbedding ord disj)
    (hanti : ∀ U V, ord U V → ord V U → U = V) :
    Function.Injective E.vertex := by
  intro U V h
  refine hanti U V ((E.le_iff U V).mpr ?_) ((E.le_iff V U).mpr ?_)
  · rw [h]
    exact causalLE_refl _
  · rw [h]
    exact causalLE_refl _

/-- Composing a cone embedding with any oriented Lorentz map yields a cone
embedding: chart invariance of the composition. -/
def mapOriented (E : ConeEmbedding ord disj) (L : OrientedLorentzEquiv) :
    ConeEmbedding ord disj where
  vertex U := L (E.vertex U)
  le_iff U V := (E.le_iff U V).trans (causalLE_map_iff L _ _).symm
  disjoint_spacelike U V h := by
    have hE := E.disjoint_spacelike U V h
    unfold IsSpacelikeSep at hE ⊢
    rwa [← OrientedLorentzEquiv.map_sub, lorentzQ_oriented]

/-- Translating a cone embedding yields a cone embedding. -/
def translate (E : ConeEmbedding ord disj) (w : Herm2) :
    ConeEmbedding ord disj where
  vertex U := E.vertex U + w
  le_iff U V := (E.le_iff U V).trans (causalLE_translate w _ _).symm
  disjoint_spacelike U V h := by
    have hE := E.disjoint_spacelike U V h
    unfold IsSpacelikeSep at hE ⊢
    rwa [add_sub_add_right_eq_sub]

end ConeEmbedding

/-! ## The committed rich-fibre region lattice in the cone -/

/-- The declared vertex assignment for the committed rich-fibre region
lattice: bottom at the origin, the four windows at scalar coordinate two
with four pairwise distinct unit spatial directions, top at scalar
coordinate four on the axis.  This assignment is a declared witness, chosen
rather than forced by the committed data. -/
def richVertex : RichRegion → Herm2
  | .bot => (0, 0)
  | .r86 => (2, Pi.single 0 1)
  | .r88 => (2, Pi.single 0 (-1))
  | .r247 => (2, Pi.single 1 1)
  | .r384 => (2, Pi.single 1 (-1))
  | .top => (4, 0)

/-- Decidability of the committed rich-region disjointness, supplied here
because the committed file declares `disj` as a plain definition. -/
instance (U V : RichRegion) : Decidable (RichRegion.disj U V) :=
  decidable_of_iff
    (U ≠ V ∧ U ≠ .bot ∧ U ≠ .top ∧ V ≠ .bot ∧ V ≠ .top) Iff.rfl

set_option maxHeartbeats 1000000 in
/-- The composition theorem for the committed rich-fibre net: one region
precedes another in the committed finite causal order exactly when the
interval of their images is future-causal in the declared Lorentz module.
Both directions are proved on the whole committed carrier. -/
theorem richRegion_le_iff_causalLE (U V : RichRegion) :
    RichRegion.le U V ↔ causalLE (richVertex U) (richVertex V) := by
  cases U <;> cases V <;>
    simp [RichRegion.le, richVertex, causalLE, IsFutureCausal, lorentzQ,
      spatialNormSq, Fin.sum_univ_three, Prod.mk_sub_mk, Pi.sub_apply,
      Pi.single_apply] <;>
    norm_num

/-- Declared disjointness of the committed rich-fibre net (distinct
windows) lands in spacelike separation. -/
theorem richDisj_spacelike (U V : RichRegion) (h : RichRegion.disj U V) :
    IsSpacelikeSep (richVertex U) (richVertex V) := by
  cases U <;> cases V <;>
    first
      | exact absurd h (by decide)
      | simp [IsSpacelikeSep, richVertex, lorentzQ, spatialNormSq,
          Fin.sum_univ_three, Prod.mk_sub_mk, Pi.sub_apply,
          Pi.single_apply] <;> norm_num

/-- The cone embedding of the committed rich-fibre region lattice. -/
def richFibreConeEmbedding : ConeEmbedding RichRegion.le RichRegion.disj where
  vertex := richVertex
  le_iff := richRegion_le_iff_causalLE
  disjoint_spacelike := richDisj_spacelike

/-- The declared vertex assignment is injective on the committed carrier. -/
theorem richVertex_injective : Function.Injective richVertex :=
  richFibreConeEmbedding.vertex_injective RichRegion.le_antisymm

/-- Nondegeneracy of the composed image: two timelike-related pairs and two
spacelike-unrelated pairs, on four distinct doubletons of the committed
carrier.  This blocks a collapse of the image to a chain or to a point. -/
theorem richVertex_nondegenerate :
    (IsTimelikeSep (richVertex .bot) (richVertex .r86) ∧
      causalLE (richVertex .bot) (richVertex .r86)) ∧
    (IsTimelikeSep (richVertex .r88) (richVertex .top) ∧
      causalLE (richVertex .r88) (richVertex .top)) ∧
    (IsSpacelikeSep (richVertex .r86) (richVertex .r88) ∧
      ¬ causalLE (richVertex .r86) (richVertex .r88) ∧
      ¬ causalLE (richVertex .r88) (richVertex .r86)) ∧
    (IsSpacelikeSep (richVertex .r247) (richVertex .r384) ∧
      ¬ causalLE (richVertex .r247) (richVertex .r384) ∧
      ¬ causalLE (richVertex .r384) (richVertex .r247)) := by
  have h1 : IsSpacelikeSep (richVertex .r86) (richVertex .r88) := by
    simp [IsSpacelikeSep, richVertex, lorentzQ, spatialNormSq,
      Fin.sum_univ_three, Prod.mk_sub_mk, Pi.sub_apply, Pi.single_apply]
    norm_num
  have h2 : IsSpacelikeSep (richVertex .r247) (richVertex .r384) := by
    simp [IsSpacelikeSep, richVertex, lorentzQ, spatialNormSq,
      Fin.sum_univ_three, Prod.mk_sub_mk, Pi.sub_apply, Pi.single_apply]
    norm_num
  refine ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩, ⟨h1, spacelike_not_related h1⟩,
    ⟨h2, spacelike_not_related h2⟩⟩ <;>
    simp [IsTimelikeSep, causalLE, IsFutureCausal, richVertex, lorentzQ,
      spatialNormSq, Prod.mk_sub_mk, Pi.single_apply] <;>
    norm_num

/-- The composition stated against the committed net itself: the region
order of `richFibreNet` at its regulator agrees with the cone order of the
images. -/
theorem richFibreNet_causal_composition (U V : RichRegion) :
    richFibreNet.regionLE () U V ↔
      causalLE (richVertex U) (richVertex V) :=
  richRegion_le_iff_causalLE U V

/-- The committed net's declared disjointness lands in spacelike
separation. -/
theorem richFibreNet_disjoint_spacelike (U V : RichRegion)
    (h : richFibreNet.disjoint () U V) :
    IsSpacelikeSep (richVertex U) (richVertex V) :=
  richDisj_spacelike U V h

/-! ## The committed partition net in the cone -/

/-- The label map from the committed partition-net carrier (subsets of the
two-site label set) into the rich-fibre region lattice. -/
def siteRegionLabel (U : Finset (Fin 2)) : RichRegion :=
  if (0 : Fin 2) ∈ U then
    (if (1 : Fin 2) ∈ U then .top else .r86)
  else
    (if (1 : Fin 2) ∈ U then .r88 else .bot)

/-- The label map is a two-directional order embedding of the committed
partition-net order into the rich-fibre region order. -/
theorem siteRegionLabel_le_iff :
    ∀ U V : Finset (Fin 2),
      U ⊆ V ↔ RichRegion.le (siteRegionLabel U) (siteRegionLabel V) := by
  decide

theorem siteRegionLabel_injective :
    ∀ U V : Finset (Fin 2), siteRegionLabel U = siteRegionLabel V → U = V := by
  decide

/-- The committed partition-net disjointness maps to rich-region
disjointness. -/
theorem siteRegionLabel_disj :
    ∀ U V : Finset (Fin 2),
      (U.Nonempty ∧ V.Nonempty ∧ Disjoint U V) →
        RichRegion.disj (siteRegionLabel U) (siteRegionLabel V) := by
  decide

/-- Sharpness of the partition carrier: the doubleton `{0}`, `{1}` is
incomparable in both directions, and every incomparable pair is that
doubleton in some order, so the committed order has exactly one unordered
incomparable pair and one spacelike doubleton is the maximum any faithful
embedding of that carrier can produce. -/
theorem partitionNet_unique_incomparable_pair :
    (¬ ({0} : Finset (Fin 2)) ⊆ {1} ∧ ¬ ({1} : Finset (Fin 2)) ⊆ {0}) ∧
    ∀ U V : Finset (Fin 2), ¬ U ⊆ V → ¬ V ⊆ U →
      (U = {0} ∧ V = {1}) ∨ (U = {1} ∧ V = {0}) := by
  decide

variable {n k : ℕ}

/-- The composition theorem for the committed partition net: the region
order of `partitionPublicCausalNet` agrees with the cone order of the
composed images, in both directions on the whole committed carrier. -/
theorem partitionNet_causal_composition
    (part : EventAlgebra.ProjectivePartition n k)
    (rho : EventAlgebra.StateMatrix n) (U V : Finset (Fin 2)) :
    (FiniteCausalObserverNet.partitionPublicCausalNet part rho).regionLE
        () U V ↔
      causalLE (richVertex (siteRegionLabel U))
        (richVertex (siteRegionLabel V)) :=
  (siteRegionLabel_le_iff U V).trans (richRegion_le_iff_causalLE _ _)

/-- The committed partition net's declared disjointness lands in spacelike
separation. -/
theorem partitionNet_disjoint_spacelike
    (part : EventAlgebra.ProjectivePartition n k)
    (rho : EventAlgebra.StateMatrix n) (U V : Finset (Fin 2))
    (h : (FiniteCausalObserverNet.partitionPublicCausalNet part
      rho).disjoint () U V) :
    IsSpacelikeSep (richVertex (siteRegionLabel U))
      (richVertex (siteRegionLabel V)) :=
  richDisj_spacelike _ _ (siteRegionLabel_disj U V h)

/-- The cone embedding of the committed partition net. -/
def partitionNetConeEmbedding
    (part : EventAlgebra.ProjectivePartition n k)
    (rho : EventAlgebra.StateMatrix n) :
    ConeEmbedding
      (fun U V : Finset (Fin 2) =>
        (FiniteCausalObserverNet.partitionPublicCausalNet part
          rho).regionLE () U V)
      (fun U V : Finset (Fin 2) =>
        (FiniteCausalObserverNet.partitionPublicCausalNet part
          rho).disjoint () U V) where
  vertex U := richVertex (siteRegionLabel U)
  le_iff := partitionNet_causal_composition part rho
  disjoint_spacelike := partitionNet_disjoint_spacelike part rho

/-! ## Negative controls -/

/-- **Negative control.**  Projecting the declared vertex assignment to the
scalar axis collapses the four windows to one point, so the chain-collapsed
assignment breaks the order equivalence.  A green order equivalence
therefore cannot come from a chain image. -/
theorem timeProjection_not_faithful :
    ¬ ∀ U V : RichRegion,
        RichRegion.le U V ↔
          causalLE ((richVertex U).1, 0) ((richVertex V).1, 0) := by
  intro h
  have h1 : causalLE ((richVertex RichRegion.r86).1, (0 : Spatial))
      ((richVertex RichRegion.r88).1, (0 : Spatial)) := causalLE_refl _
  exact (by decide : ¬ RichRegion.le .r86 .r88) ((h _ _).mpr h1)

/-- **Negative control.**  With the orientation clause dropped, the Lorentz
square alone fails antisymmetry, so the scalar-coordinate clause of
`IsFutureCausal` is load-bearing. -/
theorem q_alone_not_antisymm :
    ∃ u v : Herm2, u ≠ v ∧
      0 ≤ lorentzQ (v - u) ∧ 0 ≤ lorentzQ (u - v) := by
  refine ⟨(0, 0), (1, 0), ?_, ?_, ?_⟩
  · intro hEq
    have h1 := congrArg Prod.fst hEq
    norm_num at h1
  · simp [lorentzQ, spatialNormSq]
  · simp [lorentzQ, spatialNormSq]

/-! ## The composed receipt -/

/-- The single typed antecedent bundle of the causal-order composition: the
declared projective partition and state generating the committed partition
net.  The rich-fibre net carries no free antecedent; its payload and region
lattice are committed.  The committed partition net's `regionLE` and
`disjoint` fields do not depend on the partition or the state, so this
bundle types the receipt against the committed construction without the
causal clauses consuming the physical data; the degeneracy is disclosed in
the committed net's own docstring.  The vertex assignments `richVertex` and
`siteRegionLabel` are declared witness constructions named by the receipt,
and the receipt's final clause exhibits a second distinct witness, so the
committed data does not force the embedding. -/
structure CausalCompositionAntecedents (n k : ℕ) where
  /-- The declared projective partition of the committed partition net. -/
  partition : EventAlgebra.ProjectivePartition n k
  /-- The declared state matrix of the committed partition net. -/
  state : EventAlgebra.StateMatrix n

/-- The composed causal-order receipt.  From the single declared antecedent
bundle: the cone relation of the declared module is a partial order; the
committed rich-fibre region order agrees with the cone order of the declared
images in both directions; the image is injective and nondegenerate, with
two timelike-related and two spacelike-unrelated pairs; declared
disjointness lands in spacelike separation; the equivalence is invariant
under every oriented Lorentz map and under every affine overlap-cocycle
transition; the committed partition net composes through the declared label
map with the same order equivalence and disjointness receipt; and a second
distinct vertex assignment satisfies the same equivalence, so the embedding
is a declared choice.  No clause attaches a physical causal order,
spacetime, clock, or observer frame. -/
theorem causalOrderComposition_receipt
    (P : CausalCompositionAntecedents n k) :
    ((∀ v, causalLE v v) ∧
      (∀ u v w, causalLE u v → causalLE v w → causalLE u w) ∧
      (∀ u v, causalLE u v → causalLE v u → u = v)) ∧
    (∀ U V : RichRegion,
      richFibreNet.regionLE () U V ↔
        causalLE (richVertex U) (richVertex V)) ∧
    Function.Injective richVertex ∧
    ((IsTimelikeSep (richVertex .bot) (richVertex .r86) ∧
        causalLE (richVertex .bot) (richVertex .r86)) ∧
      (IsTimelikeSep (richVertex .r88) (richVertex .top) ∧
        causalLE (richVertex .r88) (richVertex .top)) ∧
      (IsSpacelikeSep (richVertex .r86) (richVertex .r88) ∧
        ¬ causalLE (richVertex .r86) (richVertex .r88) ∧
        ¬ causalLE (richVertex .r88) (richVertex .r86)) ∧
      (IsSpacelikeSep (richVertex .r247) (richVertex .r384) ∧
        ¬ causalLE (richVertex .r247) (richVertex .r384) ∧
        ¬ causalLE (richVertex .r384) (richVertex .r247))) ∧
    (∀ U V : RichRegion, richFibreNet.disjoint () U V →
      IsSpacelikeSep (richVertex U) (richVertex V)) ∧
    (∀ (L : OrientedLorentzEquiv) (U V : RichRegion),
      richFibreNet.regionLE () U V ↔
        causalLE (L (richVertex U)) (L (richVertex V))) ∧
    (∀ (Chart : Type*) (T : LorentzOverlapCocycle Chart) (i j : Chart)
        (U V : RichRegion),
      richFibreNet.regionLE () U V ↔
        causalLE (T.act i j (richVertex U)) (T.act i j (richVertex V))) ∧
    (∀ U V : Finset (Fin 2),
      (FiniteCausalObserverNet.partitionPublicCausalNet
          P.partition P.state).regionLE () U V ↔
        causalLE (richVertex (siteRegionLabel U))
          (richVertex (siteRegionLabel V))) ∧
    (∀ U V : Finset (Fin 2),
      (FiniteCausalObserverNet.partitionPublicCausalNet
          P.partition P.state).disjoint () U V →
        IsSpacelikeSep (richVertex (siteRegionLabel U))
          (richVertex (siteRegionLabel V))) ∧
    (∃ f : RichRegion → Herm2, f ≠ richVertex ∧
      ∀ U V : RichRegion,
        richFibreNet.regionLE () U V ↔ causalLE (f U) (f V)) := by
  refine ⟨causalLE_isPartialOrder,
    richFibreNet_causal_composition,
    richVertex_injective,
    richVertex_nondegenerate,
    richFibreNet_disjoint_spacelike,
    fun L U V => (richRegion_le_iff_causalLE U V).trans
      (causalLE_map_iff L _ _).symm,
    fun _Chart T i j U V => (richRegion_le_iff_causalLE U V).trans
      (causalLE_act_iff T i j _ _).symm,
    partitionNet_causal_composition P.partition P.state,
    partitionNet_disjoint_spacelike P.partition P.state,
    fun U => richVertex U + (1, 0), ?_, ?_⟩
  · intro hEq
    have h0 := congrFun hEq RichRegion.bot
    have h1 := congrArg Prod.fst h0
    simp [richVertex] at h1
  · intro U V
    exact (richRegion_le_iff_causalLE U V).trans
      (causalLE_translate (1, 0) _ _).symm

-- Axiom audit: no project axiom or admission is used.
#print axioms OPH.CausalComposition.causalLE_isPartialOrder
#print axioms OPH.CausalComposition.spatialDot_le_time_mul
#print axioms OPH.CausalComposition.map_futureCausal_iff
#print axioms OPH.CausalComposition.causalLE_map_iff
#print axioms OPH.CausalComposition.causalLE_act_iff
#print axioms OPH.CausalComposition.ConeEmbedding.vertex_injective
#print axioms OPH.CausalComposition.richRegion_le_iff_causalLE
#print axioms OPH.CausalComposition.richDisj_spacelike
#print axioms OPH.CausalComposition.richVertex_injective
#print axioms OPH.CausalComposition.richVertex_nondegenerate
#print axioms OPH.CausalComposition.richFibreNet_causal_composition
#print axioms OPH.CausalComposition.siteRegionLabel_le_iff
#print axioms OPH.CausalComposition.partitionNet_unique_incomparable_pair
#print axioms OPH.CausalComposition.partitionNet_causal_composition
#print axioms OPH.CausalComposition.partitionNet_disjoint_spacelike
#print axioms OPH.CausalComposition.timeProjection_not_faithful
#print axioms OPH.CausalComposition.q_alone_not_antisymm
#print axioms OPH.CausalComposition.causalOrderComposition_receipt

end

end OPH.CausalComposition
