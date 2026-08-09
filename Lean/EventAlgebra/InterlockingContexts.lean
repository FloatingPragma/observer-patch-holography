import EventAlgebra.FiniteEffectClosureBoundary

/-!
# Interlocking context webs and the Born threshold theorems

This issue-scoped module packages the committed binary-projector countermodel
and an explicit unsharp extension as exact statements about webs of
measurement contexts in the Pauli coordinate model `Herm2`.

A context is a finite family of effects summing to the identity.  A web is a
set of contexts.  An assignment gives one real value per context slot and is
normalized within every context of the web.  Interlocking is the relation
identifying equal effects across context slots, and noncontextuality asks the
assignment to respect interlocking; this is proved equivalent to factoring
through a single function on effects.

Negative threshold: on the web of all sharp binary contexts (antipodal
rank-one projector pairs), the committed cube response
`nonlinearBinaryWeight` from `EventAlgebra.FiniteEffectClosureBoundary`
defines a probability-valued noncontextual assignment with no affine/Born
coefficient vector.  Binary sharp interlocking alone underdetermines Born
form; the web packaging here consumes the committed declarations.

Partial threshold: adjoining finitely many explicit unsharp contexts excludes
the committed `z`-cubic response but is not decisive for Born form.  The
exhibited trine context has exact rational members
`(3/4) P(u1)`, `(5/8) P(u2)`, `(5/8) P(u3)` along the unit directions
`u1 = (0,0,1)`, `u2 = (4/5,0,-3/5)`, `u3 = (-4/5,0,-3/5)`.  The members sum
to the identity, are effects, are unsharp, and pairwise noncommute under the
declared Pauli matrix realization `toMatrix`; each fact is verified exactly.
Auxiliary calibration and splitting contexts force every noncontextual
normalized assignment on the extended web to extend additively across the
exhibited decompositions, with trine values `(3/4) v(P(u1))`,
`(5/8) v(P(u2))`, `(5/8) v(P(u3))`; this is the additivity hypothesis of a
finite Busch argument on only those generated decompositions.  Trine
normalization then pins one exact affine relation on three binary values, and
the cube response violates that pin exactly: its forced trine total is
`31/25`.  A transverse cubic response survives every current context.

Boundary: no affinity derivation is proved from this finite web.  The full
effect-algebra theorem in `EventAlgebra.FiniteBuschGleason` represents an
additive valuation on *all* matrix effects, but a normalized assignment on
the contexts declared here does not supply that hypothesis or a bridge from
these `Herm2` slots to a full matrix-effect valuation.  The proposed finite
bridge is recorded as `FiniteBuschGleasonInterface`; the transverse planar
counterexample in `EventAlgebra.FiniteWebBornNoGo` proves that this interface
is false for the current web.  The committed dense-test closure in
`FiniteEffectClosureBoundary` remains a positivity-after-affinity statement.
The rank-gap module
`EventAlgebra.FiniteBornFrame` (B11) treats six fixed axes; this module
quantifies over all antipodal axes and adds unsharp contexts.
-/

namespace EventAlgebra.InterlockingContexts

open OPH.C1Lorentz
open EventAlgebra.FiniteEffectClosureBoundary

noncomputable section

/-! ## Effects, contexts, webs, assignments -/

/-- Exact effect predicate in Pauli coordinates: both eigenvalues
`t ± |s|` of `t·1 + s·σ` lie in `[0,1]`, stated through squares so that
every verification is rational arithmetic. -/
def IsEffect (E : Herm2) : Prop :=
  0 ≤ E.1 ∧ E.1 ≤ 1 ∧ spatialNormSq E.2 ≤ E.1 ^ 2 ∧
    spatialNormSq E.2 ≤ (1 - E.1) ^ 2

/-- The identity effect. -/
def idEffect : Herm2 := (1, 0)

/-- A context: a finite family of effects summing to the identity.  Slots
are indexed, so one effect may occupy several slots of one context. -/
structure Context where
  size : ℕ
  member : Fin size → Herm2
  member_isEffect : ∀ j, IsEffect (member j)
  sum_eq_id : ∑ j, member j = idEffect

/-- A web of contexts. -/
abbrev Web := Set Context

/-- An assignment on a web: one real value per context slot, normalized
within every context of the web.  Values on contexts outside the web are
unconstrained. -/
structure Assignment (W : Web) where
  value : (C : Context) → Fin C.size → ℝ
  normalized : ∀ C, C ∈ W → ∑ j, value C j = 1

/-- Interlocking: two slots of two contexts carry the same effect. -/
def Interlocks (C C' : Context) (j : Fin C.size) (j' : Fin C'.size) : Prop :=
  C.member j = C'.member j'

/-- Noncontextuality: the assignment respects interlocking across the web,
including repeated slots within one context. -/
def Noncontextual (W : Web) (v : Assignment W) : Prop :=
  ∀ C, C ∈ W → ∀ C', C' ∈ W → ∀ j j',
    Interlocks C C' j j' → v.value C j = v.value C' j'

/-- Factoring through effects: a single function on effects reproduces every
slot value on the web. -/
def FactorsThroughEffects (W : Web) (v : Assignment W) : Prop :=
  ∃ f : Herm2 → ℝ, ∀ C, C ∈ W → ∀ j, v.value C j = f (C.member j)

/-- The factoring characterization of noncontextuality. -/
theorem noncontextual_iff_factorsThroughEffects (W : Web) (v : Assignment W) :
    Noncontextual W v ↔ FactorsThroughEffects W v := by
  constructor
  · intro hnc
    classical
    refine ⟨fun E =>
      if h : ∃ p : (C : Context) × Fin C.size, p.1 ∈ W ∧ p.1.member p.2 = E
      then v.value h.choose.1 h.choose.2 else 0, ?_⟩
    intro C hC j
    have hEx : ∃ p : (C' : Context) × Fin C'.size,
        p.1 ∈ W ∧ p.1.member p.2 = C.member j := ⟨⟨C, j⟩, hC, rfl⟩
    simp only [dif_pos hEx]
    exact hnc C hC hEx.choose.1 hEx.choose_spec.1 j hEx.choose.2
      hEx.choose_spec.2.symm
  · rintro ⟨f, hf⟩ C hC C' hC' j j' hjj'
    rw [hf C hC j, hf C' hC' j']
    exact congrArg f hjj'

/-- Context normalization transported through a factoring function. -/
theorem factored_context_sum {W : Web} (v : Assignment W) {f : Herm2 → ℝ}
    (hf : ∀ C, C ∈ W → ∀ j, v.value C j = f (C.member j))
    {C : Context} (hC : C ∈ W) :
    ∑ j, f (C.member j) = 1 := by
  rw [← v.normalized C hC]
  exact Finset.sum_congr rfl fun j _ => (hf C hC j).symm

/-! ## Sharp projectors and scaled projectors -/

/-- The scaled rank-one projector `c · P(n)` in Pauli coordinates. -/
def scaledProj (c : ℝ) (n : Spatial) : Herm2 := (c / 2, (c / 2) • n)

/-- The sharp rank-one projector onto direction `n`. -/
def proj (n : Spatial) : Herm2 := scaledProj 1 n

/-- The calibration grain effect `(1/8) P(n)`. -/
def grain (n : Spatial) : Herm2 := scaledProj (1 / 8) n

theorem spatialNormSq_neg (n : Spatial) :
    spatialNormSq (-n) = spatialNormSq n := by
  unfold spatialNormSq
  apply Finset.sum_congr rfl
  intro i _
  simp [Pi.neg_apply]

/-- Every scaled projector with scale in `[0,1]` along a unit direction is
an effect. -/
theorem isEffect_scaledProj {n : Spatial} (hn : spatialNormSq n = 1)
    {c : ℝ} (h0 : 0 ≤ c) (h1 : c ≤ 1) : IsEffect (scaledProj c n) := by
  have hs : spatialNormSq ((c / 2) • n) = (c / 2) ^ 2 := by
    rw [spatialNormSq_smul, hn, mul_one]
  refine ⟨by dsimp [scaledProj]; linarith, by dsimp [scaledProj]; linarith,
    ?_, ?_⟩
  · exact hs.le
  · dsimp [scaledProj]
    rw [hs]
    nlinarith

theorem isEffect_proj {n : Spatial} (hn : spatialNormSq n = 1) :
    IsEffect (proj n) :=
  isEffect_scaledProj hn zero_le_one le_rfl

theorem proj_injective : Function.Injective proj := by
  intro n m h
  have h2 : (1 / 2 : ℝ) • n = (1 / 2 : ℝ) • m := congrArg Prod.snd h
  exact smul_right_injective Spatial (by norm_num) h2

open Classical in
/-- The direction of a sharp projector, extracted by choice; junk value `0`
away from the projector image. -/
def projDirection (E : Herm2) : Spatial :=
  if h : ∃ n, spatialNormSq n = 1 ∧ E = proj n then h.choose else 0

theorem projDirection_proj {n : Spatial} (hn : spatialNormSq n = 1) :
    projDirection (proj n) = n := by
  have hEx : ∃ m, spatialNormSq m = 1 ∧ proj n = proj m := ⟨n, hn, rfl⟩
  unfold projDirection
  rw [dif_pos hEx]
  exact (proj_injective hEx.choose_spec.2).symm

/-! ## The web of all sharp binary contexts -/

/-- The sharp binary context of an antipodal projector pair. -/
def binContext (n : Spatial) (hn : spatialNormSq n = 1) : Context where
  size := 2
  member := ![proj n, proj (-n)]
  member_isEffect := by
    intro j
    fin_cases j
    · exact isEffect_proj hn
    · exact isEffect_proj (by rw [spatialNormSq_neg, hn])
  sum_eq_id := by
    rw [Fin.sum_univ_two]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one,
      proj, scaledProj, idEffect, Prod.mk_add_mk, Prod.mk.injEq]
    refine ⟨by norm_num, ?_⟩
    funext i
    simp only [Pi.add_apply, Pi.smul_apply, Pi.neg_apply, Pi.zero_apply,
      smul_eq_mul]
    ring

/-- The web of all sharp binary contexts in dimension two. -/
def binaryWeb : Web :=
  {C | ∃ n : Spatial, ∃ hn : spatialNormSq n = 1, C = binContext n hn}

theorem binContext_mem_binaryWeb (n : Spatial) (hn : spatialNormSq n = 1) :
    binContext n hn ∈ binaryWeb := ⟨n, hn, rfl⟩

/-! ## Negative threshold: the cube response on the binary web -/

/-- The cube response read as a function on effects, through the projector
direction. -/
def cubeEffectResponse (E : Herm2) : ℝ :=
  nonlinearBinaryWeight (projDirection E)

/-- The cube assignment on the binary web: the committed nonlinear frame
function, packaged per context slot. -/
def cubeAssignment : Assignment binaryWeb where
  value := fun C j => cubeEffectResponse (C.member j)
  normalized := by
    rintro C ⟨n, hn, rfl⟩
    have hn' : spatialNormSq (-n) = 1 := by rw [spatialNormSq_neg, hn]
    show (∑ j : Fin 2, cubeEffectResponse ((binContext n hn).member j)) = 1
    rw [Fin.sum_univ_two]
    show cubeEffectResponse (proj n) + cubeEffectResponse (proj (-n)) = 1
    simp only [cubeEffectResponse, projDirection_proj hn,
      projDirection_proj hn']
    exact nonlinearBinaryWeight_antipodal_sum n

/-- The cube assignment factors through effects, hence is noncontextual on
the binary web. -/
theorem cubeAssignment_noncontextual : Noncontextual binaryWeb cubeAssignment :=
  (noncontextual_iff_factorsThroughEffects binaryWeb cubeAssignment).mpr
    ⟨cubeEffectResponse, fun _ _ _ => rfl⟩

/-- Every cube slot value on the binary web is a probability. -/
theorem cubeAssignment_value_mem_Icc (n : Spatial)
    (hn : spatialNormSq n = 1) (j : Fin 2) :
    cubeAssignment.value (binContext n hn) j ∈ Set.Icc (0 : ℝ) 1 := by
  have hn' : spatialNormSq (-n) = 1 := by rw [spatialNormSq_neg, hn]
  fin_cases j
  · show cubeEffectResponse (proj n) ∈ Set.Icc (0 : ℝ) 1
    simp only [cubeEffectResponse, projDirection_proj hn]
    exact nonlinearBinaryWeight_mem_Icc hn
  · show cubeEffectResponse (proj (-n)) ∈ Set.Icc (0 : ℝ) 1
    simp only [cubeEffectResponse, projDirection_proj hn']
    exact nonlinearBinaryWeight_mem_Icc hn'

/-- Negative threshold: the web of all sharp binary contexts does not
determine Born form.  The cube response is a noncontextual per-context
normalized assignment on that web whose binary values match no affine/Born
coefficient vector. -/
theorem binary_sharp_web_underdetermines_born :
    ∃ v : Assignment binaryWeb, Noncontextual binaryWeb v ∧
      ¬ ∃ q : Spatial, ∀ (n : Spatial) (hn : spatialNormSq n = 1),
        v.value (binContext n hn) (0 : Fin 2) = affineBinaryWeight q n := by
  refine ⟨cubeAssignment, cubeAssignment_noncontextual, ?_⟩
  rintro ⟨q, hq⟩
  apply nonlinearBinaryWeight_not_affine
  refine ⟨q, fun n hn => ?_⟩
  have hval : cubeAssignment.value (binContext n hn) (0 : Fin 2)
      = nonlinearBinaryWeight n := by
    show cubeEffectResponse (proj n) = nonlinearBinaryWeight n
    simp only [cubeEffectResponse, projDirection_proj hn]
  rw [← hval]
  exact hq n hn

/-! ## The explicit unsharp battery -/

/-- First trine axis: the pole. -/
def u1 : Spatial := ![0, 0, 1]

/-- Second trine axis: a rational tilt. -/
def u2 : Spatial := ![4 / 5, 0, -(3 / 5)]

/-- Third trine axis: the mirrored rational tilt. -/
def u3 : Spatial := ![-(4 / 5), 0, -(3 / 5)]

theorem u1_unit : spatialNormSq u1 = 1 := by
  norm_num [u1, spatialNormSq, Fin.sum_univ_succ]

theorem u2_unit : spatialNormSq u2 = 1 := by
  norm_num [u2, spatialNormSq, Fin.sum_univ_succ]

theorem u3_unit : spatialNormSq u3 = 1 := by
  norm_num [u3, spatialNormSq, Fin.sum_univ_succ]

/-- First trine member `(3/4) P(u1)`. -/
def trineE1 : Herm2 := scaledProj (3 / 4) u1

/-- Second trine member `(5/8) P(u2)`. -/
def trineE2 : Herm2 := scaledProj (5 / 8) u2

/-- Third trine member `(5/8) P(u3)`. -/
def trineE3 : Herm2 := scaledProj (5 / 8) u3

/-- The unsharp three-outcome trine context: exact rational effects with
weighted axes summing to zero, so the members sum to the identity. -/
def trineContext : Context where
  size := 3
  member := ![trineE1, trineE2, trineE3]
  member_isEffect := by
    intro j
    fin_cases j
    · exact isEffect_scaledProj u1_unit (by norm_num) (by norm_num)
    · exact isEffect_scaledProj u2_unit (by norm_num) (by norm_num)
    · exact isEffect_scaledProj u3_unit (by norm_num) (by norm_num)
  sum_eq_id := by
    rw [Fin.sum_univ_three]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
      Matrix.cons_val_two, Matrix.tail_cons, trineE1, trineE2, trineE3,
      scaledProj, idEffect, Prod.mk_add_mk, Prod.mk.injEq]
    refine ⟨by norm_num, ?_⟩
    funext i
    fin_cases i <;>
      norm_num [u1, u2, u3, Pi.add_apply, Pi.smul_apply, Pi.zero_apply,
        smul_eq_mul]

/-- Calibration context along `n`: eight grain copies and the antipodal
projector. -/
def repContext (n : Spatial) (hn : spatialNormSq n = 1) : Context where
  size := 9
  member := ![grain n, grain n, grain n, grain n, grain n, grain n, grain n,
    grain n, proj (-n)]
  member_isEffect := by
    intro j
    have hg : IsEffect (grain n) :=
      isEffect_scaledProj hn (by norm_num) (by norm_num)
    have hp : IsEffect (proj (-n)) :=
      isEffect_proj (by rw [spatialNormSq_neg, hn])
    fin_cases j <;> first | exact hg | exact hp
  sum_eq_id := by
    simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero,
      Matrix.cons_val_zero, Matrix.cons_val_succ, grain, proj, scaledProj,
      idEffect, Prod.mk_add_mk, Prod.mk.injEq]
    refine ⟨by norm_num, ?_⟩
    funext i
    simp only [Pi.add_apply, Pi.smul_apply, Pi.neg_apply, Pi.zero_apply,
      smul_eq_mul]
    ring

/-- Splitting context for the polar trine member: `(3/4) P(n)`, two grains,
and the antipodal projector. -/
def splitContextA (n : Spatial) (hn : spatialNormSq n = 1) : Context where
  size := 4
  member := ![scaledProj (3 / 4) n, grain n, grain n, proj (-n)]
  member_isEffect := by
    intro j
    have hg : IsEffect (grain n) :=
      isEffect_scaledProj hn (by norm_num) (by norm_num)
    have hp : IsEffect (proj (-n)) :=
      isEffect_proj (by rw [spatialNormSq_neg, hn])
    have hE : IsEffect (scaledProj (3 / 4) n) :=
      isEffect_scaledProj hn (by norm_num) (by norm_num)
    fin_cases j <;> first | exact hE | exact hg | exact hp
  sum_eq_id := by
    simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero,
      Matrix.cons_val_zero, Matrix.cons_val_succ, grain, proj, scaledProj,
      idEffect, Prod.mk_add_mk, Prod.mk.injEq]
    refine ⟨by norm_num, ?_⟩
    funext i
    simp only [Pi.add_apply, Pi.smul_apply, Pi.neg_apply, Pi.zero_apply,
      smul_eq_mul]
    ring

/-- Splitting context for a tilted trine member: `(5/8) P(n)`, three grains,
and the antipodal projector. -/
def splitContextB (n : Spatial) (hn : spatialNormSq n = 1) : Context where
  size := 5
  member := ![scaledProj (5 / 8) n, grain n, grain n, grain n, proj (-n)]
  member_isEffect := by
    intro j
    have hg : IsEffect (grain n) :=
      isEffect_scaledProj hn (by norm_num) (by norm_num)
    have hp : IsEffect (proj (-n)) :=
      isEffect_proj (by rw [spatialNormSq_neg, hn])
    have hE : IsEffect (scaledProj (5 / 8) n) :=
      isEffect_scaledProj hn (by norm_num) (by norm_num)
    fin_cases j <;> first | exact hE | exact hg | exact hp
  sum_eq_id := by
    simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero,
      Matrix.cons_val_zero, Matrix.cons_val_succ, grain, proj, scaledProj,
      idEffect, Prod.mk_add_mk, Prod.mk.injEq]
    refine ⟨by norm_num, ?_⟩
    funext i
    simp only [Pi.add_apply, Pi.smul_apply, Pi.neg_apply, Pi.zero_apply,
      smul_eq_mul]
    ring

/-- The extended web: all sharp binary contexts together with the trine and
its calibration and splitting contexts. -/
def extendedWeb : Web :=
  binaryWeb ∪
    {trineContext,
      repContext u1 u1_unit, repContext u2 u2_unit, repContext u3 u3_unit,
      splitContextA u1 u1_unit,
      splitContextB u2 u2_unit, splitContextB u3 u3_unit}

theorem binaryWeb_subset_extendedWeb : binaryWeb ⊆ extendedWeb :=
  Set.subset_union_left

theorem trineContext_mem_extendedWeb : trineContext ∈ extendedWeb :=
  Set.mem_union_right _ (by simp)

theorem repContext_u1_mem_extendedWeb :
    repContext u1 u1_unit ∈ extendedWeb :=
  Set.mem_union_right _ (by simp)

theorem repContext_u2_mem_extendedWeb :
    repContext u2 u2_unit ∈ extendedWeb :=
  Set.mem_union_right _ (by simp)

theorem repContext_u3_mem_extendedWeb :
    repContext u3 u3_unit ∈ extendedWeb :=
  Set.mem_union_right _ (by simp)

theorem splitContextA_u1_mem_extendedWeb :
    splitContextA u1 u1_unit ∈ extendedWeb :=
  Set.mem_union_right _ (by simp)

theorem splitContextB_u2_mem_extendedWeb :
    splitContextB u2 u2_unit ∈ extendedWeb :=
  Set.mem_union_right _ (by simp)

theorem splitContextB_u3_mem_extendedWeb :
    splitContextB u3 u3_unit ∈ extendedWeb :=
  Set.mem_union_right _ (by simp)

/-! ## Exact unsharpness and pairwise noncommutation receipts -/

/-- Each trine member is unsharp: its Pauli realization is not idempotent. -/
theorem trine_members_unsharp :
    toMatrix trineE1 * toMatrix trineE1 ≠ toMatrix trineE1 ∧
    toMatrix trineE2 * toMatrix trineE2 ≠ toMatrix trineE2 ∧
    toMatrix trineE3 * toMatrix trineE3 ≠ toMatrix trineE3 := by
  refine ⟨fun h => ?_, fun h => ?_, fun h => ?_⟩ <;>
  · have h00 := congrArg (fun M : Matrix (Fin 2) (Fin 2) ℂ => M 0 0) h
    simp only [Matrix.mul_apply, Fin.sum_univ_two] at h00
    norm_num [toMatrix, trineE1, trineE2, trineE3, scaledProj, u1, u2, u3,
      Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons,
      Pi.smul_apply, smul_eq_mul] at h00

/-- The trine members pairwise noncommute in the Pauli realization. -/
theorem trine_members_pairwise_noncommuting :
    toMatrix trineE1 * toMatrix trineE2 ≠ toMatrix trineE2 * toMatrix trineE1 ∧
    toMatrix trineE1 * toMatrix trineE3 ≠ toMatrix trineE3 * toMatrix trineE1 ∧
    toMatrix trineE2 * toMatrix trineE3 ≠ toMatrix trineE3 * toMatrix trineE2 := by
  refine ⟨fun h => ?_, fun h => ?_, fun h => ?_⟩ <;>
  · have h01 := congrArg (fun M : Matrix (Fin 2) (Fin 2) ℂ => M 0 1) h
    simp only [Matrix.mul_apply, Fin.sum_univ_two] at h01
    norm_num [toMatrix, trineE1, trineE2, trineE3, scaledProj, u1, u2, u3,
      Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons,
      Pi.smul_apply, smul_eq_mul] at h01

/-! ## Positive threshold: forced additive extension and the cube exclusion -/

/-- Positive threshold, forcing form: every noncontextual assignment on the
extended web extends additively across the exhibited decompositions.  Each
trine value equals the exact scaled binary value of its axis. -/
theorem noncontextual_forces_trine_values (v : Assignment extendedWeb)
    (hnc : Noncontextual extendedWeb v) :
    v.value trineContext (0 : Fin 3)
        = (3 / 4) * v.value (binContext u1 u1_unit) (0 : Fin 2) ∧
    v.value trineContext (1 : Fin 3)
        = (5 / 8) * v.value (binContext u2 u2_unit) (0 : Fin 2) ∧
    v.value trineContext (2 : Fin 3)
        = (5 / 8) * v.value (binContext u3 u3_unit) (0 : Fin 2) := by
  obtain ⟨f, hf⟩ :=
    (noncontextual_iff_factorsThroughEffects extendedWeb v).mp hnc
  have hbin : ∀ (n : Spatial) (hn : spatialNormSq n = 1),
      f (proj n) + f (proj (-n)) = 1 := by
    intro n hn
    have hmem := binaryWeb_subset_extendedWeb (binContext_mem_binaryWeb n hn)
    have h : (∑ j : Fin 2, f ((binContext n hn).member j)) = 1 :=
      factored_context_sum v hf hmem
    rw [Fin.sum_univ_two] at h
    exact h
  have hrep : ∀ (n : Spatial) (hn : spatialNormSq n = 1),
      repContext n hn ∈ extendedWeb →
      8 * f (grain n) + f (proj (-n)) = 1 := by
    intro n hn hmem
    have h : (∑ j : Fin 9, f ((repContext n hn).member j)) = 1 :=
      factored_context_sum v hf hmem
    simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero, repContext,
      Matrix.cons_val_zero, Matrix.cons_val_succ] at h
    linarith
  have hsplitA : ∀ (n : Spatial) (hn : spatialNormSq n = 1),
      splitContextA n hn ∈ extendedWeb →
      f (scaledProj (3 / 4) n) + 2 * f (grain n) + f (proj (-n)) = 1 := by
    intro n hn hmem
    have h : (∑ j : Fin 4, f ((splitContextA n hn).member j)) = 1 :=
      factored_context_sum v hf hmem
    simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero, splitContextA,
      Matrix.cons_val_zero, Matrix.cons_val_succ] at h
    linarith
  have hsplitB : ∀ (n : Spatial) (hn : spatialNormSq n = 1),
      splitContextB n hn ∈ extendedWeb →
      f (scaledProj (5 / 8) n) + 3 * f (grain n) + f (proj (-n)) = 1 := by
    intro n hn hmem
    have h : (∑ j : Fin 5, f ((splitContextB n hn).member j)) = 1 :=
      factored_context_sum v hf hmem
    simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, add_zero, splitContextB,
      Matrix.cons_val_zero, Matrix.cons_val_succ] at h
    linarith
  have hb1 := hbin u1 u1_unit
  have hb2 := hbin u2 u2_unit
  have hb3 := hbin u3 u3_unit
  have hr1 := hrep u1 u1_unit repContext_u1_mem_extendedWeb
  have hr2 := hrep u2 u2_unit repContext_u2_mem_extendedWeb
  have hr3 := hrep u3 u3_unit repContext_u3_mem_extendedWeb
  have hA1 := hsplitA u1 u1_unit splitContextA_u1_mem_extendedWeb
  have hB2 := hsplitB u2 u2_unit splitContextB_u2_mem_extendedWeb
  have hB3 := hsplitB u3 u3_unit splitContextB_u3_mem_extendedWeb
  have hv1 : v.value trineContext (0 : Fin 3) = f (scaledProj (3 / 4) u1) :=
    hf trineContext trineContext_mem_extendedWeb (0 : Fin 3)
  have hv2 : v.value trineContext (1 : Fin 3) = f (scaledProj (5 / 8) u2) :=
    hf trineContext trineContext_mem_extendedWeb (1 : Fin 3)
  have hv3 : v.value trineContext (2 : Fin 3) = f (scaledProj (5 / 8) u3) :=
    hf trineContext trineContext_mem_extendedWeb (2 : Fin 3)
  have hw1 : v.value (binContext u1 u1_unit) (0 : Fin 2) = f (proj u1) :=
    hf (binContext u1 u1_unit)
      (binaryWeb_subset_extendedWeb (binContext_mem_binaryWeb u1 u1_unit))
      (0 : Fin 2)
  have hw2 : v.value (binContext u2 u2_unit) (0 : Fin 2) = f (proj u2) :=
    hf (binContext u2 u2_unit)
      (binaryWeb_subset_extendedWeb (binContext_mem_binaryWeb u2 u2_unit))
      (0 : Fin 2)
  have hw3 : v.value (binContext u3 u3_unit) (0 : Fin 2) = f (proj u3) :=
    hf (binContext u3 u3_unit)
      (binaryWeb_subset_extendedWeb (binContext_mem_binaryWeb u3 u3_unit))
      (0 : Fin 2)
  refine ⟨?_, ?_, ?_⟩
  · rw [hv1, hw1]; linarith
  · rw [hv2, hw2]; linarith
  · rw [hv3, hw3]; linarith

/-- Trine normalization pins an exact affine relation on three binary
values of any noncontextual assignment on the extended web. -/
theorem trine_pins_binary_values (v : Assignment extendedWeb)
    (hnc : Noncontextual extendedWeb v) :
    (3 / 4) * v.value (binContext u1 u1_unit) (0 : Fin 2)
      + (5 / 8) * v.value (binContext u2 u2_unit) (0 : Fin 2)
      + (5 / 8) * v.value (binContext u3 u3_unit) (0 : Fin 2) = 1 := by
  obtain ⟨h1, h2, h3⟩ := noncontextual_forces_trine_values v hnc
  have hsum : v.value trineContext (0 : Fin 3)
      + v.value trineContext (1 : Fin 3)
      + v.value trineContext (2 : Fin 3) = 1 := by
    have h : (∑ j : Fin 3, v.value trineContext j) = 1 :=
      v.normalized trineContext trineContext_mem_extendedWeb
    rw [Fin.sum_univ_three] at h
    exact h
  linarith

/-- The exact cube values on the three trine axes. -/
theorem cube_values_on_trine_axes :
    nonlinearBinaryWeight u1 = 1 ∧
    nonlinearBinaryWeight u2 = 49 / 125 ∧
    nonlinearBinaryWeight u3 = 49 / 125 := by
  have h1 : u1 2 = 1 := rfl
  have h2 : u2 2 = -(3 / 5) := rfl
  have h3 : u3 2 = -(3 / 5) := rfl
  refine ⟨?_, ?_, ?_⟩ <;>
    simp only [nonlinearBinaryWeight, h1, h2, h3] <;> norm_num

/-- The forced trine total of the cube response, computed exactly. -/
theorem cube_forced_trine_total :
    (3 / 4) * nonlinearBinaryWeight u1
      + (5 / 8) * nonlinearBinaryWeight u2
      + (5 / 8) * nonlinearBinaryWeight u3 = 31 / 25 := by
  obtain ⟨h1, h2, h3⟩ := cube_values_on_trine_axes
  rw [h1, h2, h3]
  norm_num

/-- Positive threshold, exclusion form: the cube response has no
noncontextual extension to the extended web.  Its binary restrictions force
the trine values, and the forced total `31/25` violates the three-element
normalization. -/
theorem cube_response_no_noncontextual_extension :
    ¬ ∃ v : Assignment extendedWeb, Noncontextual extendedWeb v ∧
      ∀ C, C ∈ binaryWeb → ∀ j : Fin C.size,
        v.value C j = cubeAssignment.value C j := by
  rintro ⟨v, hnc, hres⟩
  have hpin := trine_pins_binary_values v hnc
  obtain ⟨hc1, hc2, hc3⟩ := cube_values_on_trine_axes
  have hval : ∀ (n : Spatial) (hn : spatialNormSq n = 1),
      v.value (binContext n hn) (0 : Fin 2) = nonlinearBinaryWeight n := by
    intro n hn
    rw [hres (binContext n hn) (binContext_mem_binaryWeb n hn) (0 : Fin 2)]
    show cubeEffectResponse (proj n) = nonlinearBinaryWeight n
    simp only [cubeEffectResponse, projDirection_proj hn]
  rw [hval u1 u1_unit, hval u2 u2_unit, hval u3 u3_unit, hc1, hc2, hc3]
    at hpin
  norm_num at hpin

/-! ## Bridge interface for the finite Busch-Gleason lane -/

/-- Proposed consumption interface for a finite Busch-Gleason lane.  This
`Prop` asserts that every noncontextual assignment on the
extended web restricts on the binary subweb to an affine/Born functional
with coefficient vector in the closed unit ball, in the committed
`affineBinaryWeight` coordinates.  Together with
`binary_sharp_web_underdetermines_born` and
`cube_response_no_noncontextual_extension`, it was intended to locate the
Born decision point at the unsharp interlocking added by `extendedWeb`.
However, `EventAlgebra.FiniteWebBornNoGo` constructs a probability-valued
transverse cubic assignment and proves this `Prop` false.  The full
representation theorem in `EventAlgebra.FiniteBuschGleason` therefore cannot
be consumed without strengthening the web hypotheses. -/
def FiniteBuschGleasonInterface : Prop :=
  ∀ v : Assignment extendedWeb, Noncontextual extendedWeb v →
    ∃ q : Spatial, spatialNormSq q ≤ 1 ∧
      ∀ (n : Spatial) (hn : spatialNormSq n = 1),
        v.value (binContext n hn) (0 : Fin 2) = affineBinaryWeight q n

#print axioms noncontextual_iff_factorsThroughEffects
#print axioms cubeAssignment_noncontextual
#print axioms cubeAssignment_value_mem_Icc
#print axioms binary_sharp_web_underdetermines_born
#print axioms trine_members_unsharp
#print axioms trine_members_pairwise_noncommuting
#print axioms noncontextual_forces_trine_values
#print axioms trine_pins_binary_values
#print axioms cube_response_no_noncontextual_extension

end

end EventAlgebra.InterlockingContexts
