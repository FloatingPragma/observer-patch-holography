/-
  A5_OPH.lean  —  consolidated, machine-checked corpus

  PROVENANCE.  Contributed by DULA (external cross-audit lane), built with
  Aristotle against Lean/Mathlib v4.28.0, merged 2026-07-23.  Rebuilt
  unchanged under this project's toolchain (v4.29.1) with the same
  eighteen-line standard-axiom audit; drop artifacts and the rebuild record
  are in `DULA/` at the meta-repo root (`VERIFICATION.md`).  This header
  note is the only local edit.

  Five separately verified Lean developments, merged into one module so that a
  single `lake build` covers all of them.  Every proof body below is the one
  Aristotle compiled against Mathlib v4.28.0; nothing has been re-proved by
  hand during the merge.  The only editorial change is the restoration of the
  Part IV header, which a previous run deleted.

  ------------------------------------------------------------------------
  PART   NAMESPACE          CONTENT                                   STATUS
  ------------------------------------------------------------------------
  I      A5Order6           A5 has a subgroup of order 6, ≅ S3         proved
  II     A5NoZ6             A5 has NO Z6; S5 does                      proved
  III    OPHGap             G-module iso does NOT fix a Lie bracket    proved
  IV     OPHSelection       arithmetic + group-theoretic core of the
                            u(1)+su(2)+su(3) selection argument        proved
  V      OPHInner           (H2) ⟹ (H*) reduction; A5 ⊄ SU(2)          proved
  VI     OPHTrichotomy      excluded semisimple dims; quintet
                            noncentrality; Z6 screen gluing class      proved
  ------------------------------------------------------------------------

  WHAT THE CORPUS ESTABLISHES, IN ONE LINE EACH.

    * A5 contains exactly the order-6 subgroups one expects (Part I), and
      contains no Z6 at all (Part II).  Since SL(2,5) = 2·A5 DOES contain Z6,
      any structure requiring the mod-6 unit group must live on the double
      cover, not on A5.
    * A representation isomorphism carries no Lie bracket (Part III).  Hence
      "R12 ≅ u(1)+su(2)+su(3) as A5-modules" does not by itself yield the
      Standard Model gauge algebra.
    * Given the classical Lie-theoretic inputs, the selection argument's finite
      core is correct (Part IV): sums of compact-simple dimensions equal to 11
      force {3,8}; equal to 12 force {3,3,3,3}; and A5 cannot permute four
      simple factors, so su(2)^4 is excluded and the centre is nonzero.
    * The hypothesis "A5 acts trivially on Z(g)" is not an extra assumption if
      the screen symmetry is a gauge symmetry (Part V), because inner
      automorphisms fix the centre.  And that hypothesis has a consequence
      that could have failed: A5 does not embed in SU(2).

  WHAT THE CORPUS DOES *NOT* ESTABLISH.

    * That OPH derives its own hypotheses.  Parts IV and V are conditional.
    * Two classical inputs are used and are NOT in Mathlib: surjectivity of
      exp on compact connected groups, and the classification of compact
      simple Lie algebras.  They are marked in the part headers.
    * NOTHING here concerns matter content.  There is no result in this file
      about fermion representations, generation number, hypercharge, or the
      Higgs sector.  Any claim in that direction is outside this corpus.

  Axiom audit is at the end of the file.
-/

import Mathlib

universe u v

open Equiv Equiv.Perm

/- ==================================================================== -/
/- PART I — A5 has a subgroup of order 6                                -/
/- ==================================================================== -/

namespace A5Order6

/-! ### Step 1: ℤˣ → Perm (Fin 2) -/

/-- The isomorphism `ℤˣ ≃ Perm (Fin 2)` as a monoid hom: `1 ↦ id`, `-1 ↦ (0 1)`. -/
def unitsToPerm : ℤˣ →* Perm (Fin 2) where
  toFun u := if u = 1 then 1 else swap 0 1
  map_one' := by simp
  map_mul' := by
    intro u v
    rcases Int.units_eq_one_or u with hu | hu <;>
      rcases Int.units_eq_one_or v with hv | hv <;>
        subst hu <;> subst hv <;>
          simp [Equiv.swap_mul_self, (by decide : (-1 : ℤˣ) ≠ 1)]

@[simp] lemma sign_unitsToPerm (u : ℤˣ) : sign (unitsToPerm u) = u := by
  rcases Int.units_eq_one_or u with hu | hu <;> subst hu <;>
    simp [unitsToPerm, (by decide : (-1 : ℤˣ) ≠ 1),
      Equiv.Perm.sign_swap (by decide : (0 : Fin 2) ≠ 1)]

/-! ### Step 2: the twisted map S₃ → Perm (Fin 3 ⊕ Fin 2) -/

/-- `σ ↦ σ ⊕ (sign σ)`. -/
def twist : Perm (Fin 3) →* Perm (Fin 3 ⊕ Fin 2) :=
  (sumCongrHom (Fin 3) (Fin 2)).comp
    ((MonoidHom.id (Perm (Fin 3))).prod (unitsToPerm.comp (sign : Perm (Fin 3) →* ℤˣ)))

lemma twist_apply (σ : Perm (Fin 3)) :
    twist σ = Equiv.sumCongr σ (unitsToPerm (sign σ)) := rfl

lemma twist_injective : Function.Injective twist := by
  intro a b h
  have := Equiv.Perm.sumCongrHom_injective h
  exact congrArg Prod.fst this

lemma sign_twist (σ : Perm (Fin 3)) : sign (twist σ) = 1 := by
  rw [twist_apply]
  rw [Equiv.Perm.sign_sumCongr, sign_unitsToPerm]
  exact Int.units_mul_self _

/-! ### Step 3: relabel `Fin 3 ⊕ Fin 2 ≃ Fin 5` -/

/-- `Equiv.permCongr` upgraded to a monoid isomorphism. -/
def permCongrHom {α β : Type*} (e : α ≃ β) : Perm α ≃* Perm β where
  toFun := e.permCongr
  invFun := e.symm.permCongr
  left_inv := by intro p; ext x; simp [Equiv.permCongr_apply]
  right_inv := by intro p; ext x; simp [Equiv.permCongr_apply]
  map_mul' := by intro p q; ext x; simp [Equiv.permCongr_apply]

/-- The relabelling `Fin 3 ⊕ Fin 2 ≃ Fin 5`. -/
def e5 : Fin 3 ⊕ Fin 2 ≃ Fin 5 := finSumFinEquiv.trans (finCongr (by norm_num))

/-- The twisted embedding `S₃ →* S₅`. -/
def emb : Perm (Fin 3) →* Perm (Fin 5) :=
  (permCongrHom e5).toMonoidHom.comp twist

lemma emb_injective : Function.Injective emb :=
  (permCongrHom e5).injective.comp twist_injective

lemma sign_emb (σ : Perm (Fin 3)) : sign (emb σ) = 1 := by
  show sign (e5.permCongr (twist σ)) = 1
  rw [Equiv.Perm.sign_permCongr]
  exact sign_twist σ

/-! ### Step 4: the subgroup -/

/-- The image of the twisted embedding: a copy of `S₃` inside `S₅`. -/
def H : Subgroup (Perm (Fin 5)) := emb.range

lemma H_le_alternating : H ≤ alternatingGroup (Fin 5) := by
  rintro x ⟨σ, rfl⟩
  simpa [Equiv.Perm.mem_alternatingGroup] using sign_emb σ

lemma card_H : Nat.card H = 6 := by
  have h1 : Nat.card H = Nat.card (Perm (Fin 3)) :=
    (Nat.card_congr (Equiv.ofInjective emb emb_injective)).symm
  rw [h1, Nat.card_eq_fintype_card, Fintype.card_perm]
  norm_num

/-- **A₅ has a subgroup of order 6.** -/
theorem exists_subgroup_card_six_le_alternating :
    ∃ H : Subgroup (Perm (Fin 5)), H ≤ alternatingGroup (Fin 5) ∧ Nat.card H = 6 :=
  ⟨H, H_le_alternating, card_H⟩

/-- Same statement, phrased as a subgroup *of* `A₅`. -/
theorem exists_subgroup_card_six :
    ∃ K : Subgroup (alternatingGroup (Fin 5)), Nat.card K = 6 := by
  refine ⟨Subgroup.comap (alternatingGroup (Fin 5)).subtype H, ?_⟩
  have : Nat.card (Subgroup.comap (alternatingGroup (Fin 5)).subtype H)
       = Nat.card H := by
    apply Nat.card_congr
    refine ⟨fun x => ⟨x.1.1, x.2⟩, fun y => ⟨⟨y.1, H_le_alternating y.2⟩, y.2⟩,
      fun _ => rfl, fun _ => rfl⟩
  rw [this, card_H]

/-- It is nonabelian, hence `≅ S₃` rather than `ℤ/6`. -/
theorem H_not_commutative : ¬ ∀ x y : H, x * y = y * x := by
  intro hcomm
  have : ∀ a b : Perm (Fin 3), a * b = b * a := by
    intro a b
    have := hcomm ⟨emb a, ⟨a, rfl⟩⟩ ⟨emb b, ⟨b, rfl⟩⟩
    have h2 : emb (a * b) = emb (b * a) := by
      simpa [map_mul] using congrArg Subtype.val this
    exact emb_injective h2
  have : (c[(0 : Fin 3), 1] * c[(0 : Fin 3), 1, 2] : Perm (Fin 3))
       = c[(0 : Fin 3), 1, 2] * c[(0 : Fin 3), 1] := this _ _
  revert this
  decide

end A5Order6


/- ==================================================================== -/
/- PART II — A5 has no Z6 (and S5 does)                                 -/
/- ==================================================================== -/

/-
  A5NoZ6.lean

  THEOREM (narrow, and this is the only thing claimed):
      A₅ = alternatingGroup (Fin 5) contains no element of order 6,
      hence no subgroup isomorphic to ℤ/6.

  NON-THEOREM (retracted): "no clever embedding fixes it."  FALSE.
  The obstruction is specific to A₅ as a permutation group of degree 5.
  Verified counterexamples in nearby groups:
      S₅   = A₅ ⋊ ℤ/2      element orders 1,2,3,4,5,6   -- e.g. (01)(234)
      SL(2,5) = 2·A₅        element orders 1,2,3,4,5,6,10 -- binary icosahedral
      A₇                    element orders 1,...,6,7      -- e.g. (012)(34)(56)
  (A₆ does NOT: orders 1,2,3,4,5.)
  So passing to the double cover, or to the full symmetric group, does
  supply μ₆. Only the *specific* group A₅ is obstructed.

  Proof strategy: cycle-type arithmetic, no kernel enumeration.
  An order-6 element of S₅ must have cycleType = {2,3} (lcm 6, entries ≥ 2
  dividing 6, sum ≤ 5), and sign = (-1)^(sum + card) = (-1)^(5+2) = -1.

  STATUS: compiled with Mathlib. No sorries.
-/



namespace A5NoZ6

variable {g : Perm (Fin 5)}

/-! ### Cycle-type constraints -/

private lemma dvd_six_of_mem (h : orderOf g = 6) {a : ℕ} (ha : a ∈ g.cycleType) :
    a ∣ 6 := by
  rw [← h, ← Equiv.Perm.lcm_cycleType]
  exact Multiset.dvd_lcm ha

private lemma sum_le_five (g : Perm (Fin 5)) : g.cycleType.sum ≤ 5 := by
  rw [Equiv.Perm.sum_cycleType]
  simpa using g.support.card_le_univ

private lemma card_le_two (g : Perm (Fin 5)) : g.cycleType.card ≤ 2 := by
  have h2 : ∀ a ∈ g.cycleType, 2 ≤ a := fun a ha =>
    Equiv.Perm.two_le_of_mem_cycleType ha
  have := Multiset.card_nsmul_le_sum h2
  have h5 := sum_le_five g
  simp only [smul_eq_mul] at this
  omega

/-- The cycle type of an order-6 permutation of five points is forced. -/
theorem cycleType_eq (h : orderOf g = 6) : g.cycleType = {2, 3} := by
  have hlcm : g.cycleType.lcm = 6 := by rw [Equiv.Perm.lcm_cycleType, h]
  have hsum := sum_le_five g
  have h2 : ∀ a ∈ g.cycleType, 2 ≤ a := fun a ha => Equiv.Perm.two_le_of_mem_cycleType ha
  have hcard := card_le_two g
  interval_cases hc : g.cycleType.card
  · -- card = 0 : cycleType = 0, lcm = 1 ≠ 6
    rw [Multiset.card_eq_zero] at hc
    rw [hc] at hlcm; simp at hlcm
  · -- card = 1 : cycleType = {a}, lcm = a = 6, but then sum = 6 > 5
    obtain ⟨a, ha⟩ := Multiset.card_eq_one.mp hc
    rw [ha] at hlcm hsum
    simp [Multiset.lcm_singleton] at hlcm hsum
    omega
  · -- card = 2 : cycleType = {a,b} with 2 ≤ a,b, a ∣ 6, b ∣ 6, a + b ≤ 5
    obtain ⟨a, b, hab⟩ := Multiset.card_eq_two.mp hc
    have hda : a ∣ 6 := dvd_six_of_mem h (by rw [hab]; simp)
    have hdb : b ∣ 6 := dvd_six_of_mem h (by rw [hab]; simp)
    have hb2 : 2 ≤ b := h2 b (by rw [hab]; simp)
    have ha2 : 2 ≤ a := h2 a (by rw [hab]; simp)
    rw [hab] at hsum hlcm ⊢
    simp [Multiset.lcm_cons, Multiset.lcm_singleton] at hlcm hsum ⊢
    -- a, b ∈ {2,3}, Nat.lcm a b = 6 forces {a,b} = {2,3}
    have ha6 : a ≤ 6 := Nat.le_of_dvd (by omega) hda
    have hb6 : b ≤ 6 := Nat.le_of_dvd (by omega) hdb
    interval_cases a <;> interval_cases b <;> simp_all
    exact Multiset.pair_comm 3 2

/-! ### Main theorem -/

/-- **A₅ has no element of order 6.** -/
theorem no_orderOf_six (hg : g ∈ alternatingGroup (Fin 5)) : orderOf g ≠ 6 := by
  intro h
  have hct := cycleType_eq h
  have hsign : Equiv.Perm.sign g = -1 := by
    rw [Equiv.Perm.sign_of_cycleType, hct]; decide
  rw [Equiv.Perm.mem_alternatingGroup] at hg
  rw [hg] at hsign
  exact absurd hsign (by decide)

/-- **A₅ contains no cyclic subgroup of order 6.** -/
theorem no_cyclic_six (H : Subgroup (Perm (Fin 5)))
    (hle : H ≤ alternatingGroup (Fin 5)) (hcyc : IsCyclic H) :
    Nat.card H ≠ 6 := by
  intro hcard
  obtain ⟨x, hx⟩ := hcyc.exists_generator
  have hord : orderOf x = 6 := by
    rw [← hcard]
    exact orderOf_eq_card_of_forall_mem_zpowers hx
  have hord' : orderOf (x : Perm (Fin 5)) = 6 := by
    rwa [← Subgroup.orderOf_coe] at hord
  exact no_orderOf_six (hle x.2) hord'

/-- **No injective homomorphism ℤ/6 → A₅.** -/
theorem no_embedding_zmod_six
    (f : Multiplicative (ZMod 6) →* alternatingGroup (Fin 5)) :
    ¬ Function.Injective f := by
  intro hf
  set x := f (Multiplicative.ofAdd (1 : ZMod 6)) with hxdef
  have h6 : orderOf (Multiplicative.ofAdd (1 : ZMod 6)) = 6 := by
    exact ZMod.addOrderOf_one 6
  have : orderOf x = 6 := by rw [hxdef, orderOf_injective f hf, h6]
  have hx : orderOf ((x : alternatingGroup (Fin 5)) : Perm (Fin 5)) = 6 := by
    rwa [Subgroup.orderOf_coe]
  exact no_orderOf_six x.2 hx

/-! ### Sharpness: the obstruction is about A₅ only, not about ℤ/6 -/

/-- `(01)(234)` has order 6 in S₅ — so `S₅ ⊇ ℤ/6`. It is odd, which is
    exactly why it misses A₅. -/
example : orderOf (c[(0 : Fin 5), 1] * c[(2 : Fin 5), 3, 4]) = 6 := by
  set_option maxRecDepth 100000 in
    apply orderOf_eq_of_pow_and_pow_div_prime (by norm_num)
    · decide
    · intro p hp hpd
      have hp_le : p ≤ 6 := Nat.le_of_dvd (by norm_num) hpd
      interval_cases p
      · norm_num [Nat.prime_def] at hp
      · norm_num [Nat.prime_def] at hp
      · decide
      · decide
      · norm_num at hpd
      · norm_num at hpd
      · norm_num [Nat.prime_def] at hp

example : (c[(0 : Fin 5), 1] * c[(2 : Fin 5), 3, 4]) ∉ alternatingGroup (Fin 5) := by
  intro h
  exact no_orderOf_six h (by
    set_option maxRecDepth 100000 in
      apply orderOf_eq_of_pow_and_pow_div_prime (by norm_num)
      · decide
      · intro p hp hpd
        have hp_le : p ≤ 6 := Nat.le_of_dvd (by norm_num) hpd
        interval_cases p
        · norm_num [Nat.prime_def] at hp
        · norm_num [Nat.prime_def] at hp
        · decide
        · decide
        · norm_num at hpd
        · norm_num at hpd
        · norm_num [Nat.prime_def] at hp)

end A5NoZ6


/- ==================================================================== -/
/- PART III — module iso does not determine a Lie bracket               -/
/- ==================================================================== -/

/-
  ModuleIsoNotLieIso.lean

  CONTEXT.  The OPH claim "the exact A₅ coefficient algebra u(1) ⊕ su(2) ⊕ su(3)"
  rests on an isomorphism of A₅-representations

      R₁₂  ≅  1 ⊕ 3 ⊕ 3′ ⊕ 5  ≅  u(1) ⊕ su(2) ⊕ su(3)   (as vector spaces)

  and then treats the right-hand side as a Lie algebra.  This file formalizes
  why that step does not go through.

  WHAT IS PROVED (rung L0):
      A G-equivariant linear isomorphism between a Lie algebra and another
      G-module does NOT determine a Lie bracket.  Explicitly, every non-abelian
      Lie algebra L with a G-action is G-module-isomorphic to a Lie algebra
      (itself, with the zero bracket) to which it is NOT Lie-isomorphic.

  WHAT IS *NOT* PROVED, AND IS NOT CLAIMED:
      L1.  That ≥ 2 non-isomorphic A₅-equivariant brackets exist on R₁₂.
           (TRUE on paper — take the zero bracket, and take
            u(1) ⊕ su(2)[via 3] ⊕ su(3)[via 3′], whose adjoint restricts to
            1 ⊕ 3 ⊕ 3′ ⊕ 5 ≅ R₁₂ — but Mathlib has no compact real forms
            with finite group actions, so it is not formalizable here.)
      L2.  dim Hom_{A₅}(Λ²R₁₂, R₁₂) = 14.  Needs the A₅ character table.
      L3.  Uniqueness of the SM bracket under any stated selection principle.
           Needs the classification of compact simple Lie algebras.

  THEREFORE: this file establishes that the OPH inference is INVALID.
  It does not establish that the OPH conclusion is FALSE.  Those are
  different claims and only the first is settled here.

  STATUS: machine-checked with Lean and Mathlib.  No sorries.
-/



namespace OPHGap

variable (R : Type u) [CommRing R]

/-! ### The same module, with the bracket forgotten

`Triv M` is a type synonym for `M` carrying the zero bracket.  Instance
resolution matches on the head symbol `Triv`, exactly as for Mathlib's
`Additive` / `Multiplicative` / `OrderDual`, so the original bracket on `M`
is not picked up here. -/

/-- Type synonym: `M` equipped with the zero Lie bracket. -/
def Triv (M : Type v) : Type v := M

instance instAddCommGroupTriv (M : Type v) [AddCommGroup M] :
    AddCommGroup (Triv M) := inferInstanceAs (AddCommGroup M)

instance instModuleTriv (M : Type v) [AddCommGroup M] [Module R M] :
    Module R (Triv M) := inferInstanceAs (Module R M)

instance instBracketTriv (M : Type v) [AddCommGroup M] :
    Bracket (Triv M) (Triv M) := ⟨fun _ _ => 0⟩

@[simp] lemma triv_bracket (M : Type v) [AddCommGroup M] (x y : Triv M) :
    ⁅x, y⁆ = 0 := rfl

instance instLieRingTriv (M : Type v) [AddCommGroup M] : LieRing (Triv M) where
  add_lie := by intros; simp
  lie_add := by intros; simp
  lie_self := by intros; simp
  leibniz_lie := by intros; simp

instance instLieAlgebraTriv (M : Type v) [AddCommGroup M] [Module R M] :
    LieAlgebra R (Triv M) where
  lie_smul := by intros; simp

instance instIsLieAbelianTriv (M : Type v) [AddCommGroup M] :
    IsLieAbelian (Triv M) := ⟨fun _ _ => rfl⟩

/-- The `G`-action transports to the synonym unchanged. -/
instance instDistribMulActionTriv (G : Type u) [Monoid G] (M : Type v)
    [AddCommGroup M] [DistribMulAction G M] : DistribMulAction G (Triv M) :=
  inferInstanceAs (DistribMulAction G M)

/-! ### The two comparisons -/

variable {R}
variable {L : Type v} [LieRing L] [LieAlgebra R L]

/-- Forgetting the bracket is a linear isomorphism. -/
def forgetBracket : L ≃ₗ[R] Triv L := LinearEquiv.refl R L

/-- …and it is `G`-equivariant, by definition of the action on `Triv`. -/
lemma forgetBracket_equivariant {G : Type u} [Monoid G] [DistribMulAction G L]
    (g : G) (x : L) : forgetBracket (R := R) (g • x) = g • forgetBracket (R := R) x :=
  rfl

/-- Lie isomorphisms transport abelianness backwards. -/
lemma isLieAbelian_of_lieEquiv {L' : Type v} [LieRing L'] [LieAlgebra R L']
    [IsLieAbelian L'] (e : L ≃ₗ⁅R⁆ L') : IsLieAbelian L := by
  refine ⟨fun x y => ?_⟩
  apply e.injective
  simp

/-- **No Lie isomorphism to the trivial-bracket synonym, when `L` is nonabelian.** -/
theorem isEmpty_lieEquiv_triv (h : ¬ IsLieAbelian L) :
    IsEmpty (L ≃ₗ⁅R⁆ Triv L) :=
  ⟨fun e => h (isLieAbelian_of_lieEquiv e)⟩

/-! ### The statement of the gap -/

/-- **Rung L0.**  For any group `G` acting on a non-abelian Lie algebra `L`,
there is another `G`-module Lie algebra that is `G`-linearly isomorphic to `L`
but not Lie-isomorphic to it.  Hence a `G`-equivariant linear isomorphism
carries no information about the bracket: the inference

    "V ≅ L as G-modules"  ⟹  "V carries the Lie algebra structure of L"

is invalid. -/
theorem gModule_iso_does_not_determine_bracket
    (G : Type u) [Monoid G] [DistribMulAction G L] (h : ¬ IsLieAbelian L) :
    (∃ e : L ≃ₗ[R] Triv L, ∀ (g : G) (x : L), e (g • x) = g • e x)
      ∧ IsEmpty (L ≃ₗ⁅R⁆ Triv L) :=
  ⟨⟨forgetBracket, fun g x => forgetBracket_equivariant g x⟩,
   isEmpty_lieEquiv_triv h⟩

/-! ### Sanity instance: something nonabelian to feed it

`sl 2 R` over a nontrivial commutative ring is nonabelian, so the hypothesis
`¬ IsLieAbelian L` is not vacuous. -/

example : ¬ IsLieAbelian (Matrix (Fin 2) (Fin 2) ℝ) := by
  intro habel
  let E₁₂ : Matrix (Fin 2) (Fin 2) ℝ := Matrix.single 0 1 1
  let E₂₁ : Matrix (Fin 2) (Fin 2) ℝ := Matrix.single 1 0 1
  have h := habel.trivial E₁₂ E₂₁
  rw [Ring.lie_def] at h
  have h00 := congrArg (fun M => M 0 0) h
  norm_num [E₁₂, E₂₁, Matrix.mul_apply, Matrix.single] at h00

end OPHGap


/- ==================================================================== -/
/- PART IV — selection argument: finite core                            -/
/- ==================================================================== -/

/-
  OPHSelection.lean

  TARGET THEOREM (stated here, only partly formalizable in current Mathlib):

    Let L be a compact reductive real Lie algebra whose underlying A₅-module is
    R₁₂ ≅ 1 ⊕ 3 ⊕ 3′ ⊕ 5, with A₅ acting by Lie algebra automorphisms.
    If A₅ acts trivially on the centre Z(L), then

        L ≅ u(1) ⊕ su(2) ⊕ su(3),   uniquely.

  PROOF SKELETON (five steps; which are formalizable is marked):
    (S1) L = Z(L) ⊕ [L,L] with both A₅-submodules.        [Lie theory: NOT in Mathlib]
    (S2) A₅ permutes the simple factors of [L,L]; since A₅ is simple of order 60
         it has no subgroup of index 2, 3 or 4, so that permutation action is
         trivial and each simple factor is an A₅-submodule.   [*** THIS FILE ***]
    (S3) dim Z(L) is a subset-sum of {1,3,3,5}; in particular 2 is impossible.
                                                              [*** THIS FILE ***]
    (S4) dim[L,L] is a sum of compact-simple dimensions, which below 12 are
         {3, 8, 10}.  Sum 11 forces {3,8}; sum 12 forces {3,3,3,3}.
                                                              [*** THIS FILE ***]
    (S5) A₅-trivial centre ⇒ dim Z ≤ 1 (trivial isotypic multiplicity of R₁₂
         is 1); equivariance ⇒ dim Z ≠ 0 (else [L,L] = su(2)^4, needing four
         3-dimensional A₅-submodules, but R₁₂ has only 3 and 3′).
         Hence dim Z = 1 and dim[L,L] = 11.   [needs A₅ character theory]

  WHAT THIS FILE PROVES: S2, S3, S4 — the group-theoretic and arithmetic core.
  WHAT IT DOES NOT: S1 and S5, which need Lie-theoretic infrastructure
  (reductive decomposition, compact real forms, the classification of compact
  simple Lie algebras) and the A₅ character table.  None of that is in Mathlib.

  HONEST FRAMING: the hypothesis "A₅ acts trivially on Z(L)" was found by
  working backwards from the desired conclusion.  This file establishes that
  the hypothesis SUFFICES.  It says nothing about whether OPH independently
  DERIVES it.  Those are different claims.

  STATUS: written, NOT compiled here.  `-- RISK:` marks names to check.
  No sorries.
-/
namespace OPHSelection

/-! ## S4: the arithmetic of compact-simple dimensions -/

/-- Dimensions of compact simple Lie algebras that are at most 12. -/
def SimpleDims : Finset ℕ := {3, 8, 10}

/-- A multiset of allowed compact-simple dimensions summing to 11 is exactly
`{3, 8}`, corresponding to `su(2) ⊕ su(3)`. -/
theorem sum_eq_eleven (m : Multiset ℕ) (hm : ∀ d ∈ m, d ∈ SimpleDims)
    (hs : m.sum = 11) : m = {3, 8} := by
  have h3 : ∀ d ∈ m, 3 ≤ d := by
    intro d hd
    have := hm d hd
    fin_cases this <;> omega
  have hcard : m.card ≤ 3 := by
    have := Multiset.card_nsmul_le_sum h3
    simp only [smul_eq_mul] at this
    omega
  interval_cases hc : m.card
  · rw [Multiset.card_eq_zero] at hc
    rw [hc] at hs
    simp at hs
  · obtain ⟨a, ha⟩ := Multiset.card_eq_one.mp hc
    have := hm a (by rw [ha]; simp)
    rw [ha] at hs
    simp at hs
    fin_cases this <;> omega
  · obtain ⟨a, b, hab⟩ := Multiset.card_eq_two.mp hc
    have hA := hm a (by rw [hab]; simp)
    have hB := hm b (by rw [hab]; simp)
    rw [hab] at hs ⊢
    simp at hs
    fin_cases hA <;> fin_cases hB <;> simp_all
    exact Multiset.pair_comm 8 3
  · obtain ⟨a, b, c, habc⟩ := Multiset.card_eq_three.mp hc
    have hA := hm a (by rw [habc]; simp)
    have hB := hm b (by rw [habc]; simp)
    have hC := hm c (by rw [habc]; simp)
    rw [habc] at hs
    simp at hs
    fin_cases hA <;> fin_cases hB <;> fin_cases hC <;> omega

/-- A multiset of allowed compact-simple dimensions summing to 12 consists of
four copies of 3, corresponding to `su(2)^4`. -/
theorem sum_eq_twelve (m : Multiset ℕ) (hm : ∀ d ∈ m, d ∈ SimpleDims)
    (hs : m.sum = 12) : m = {3, 3, 3, 3} := by
  have h3 : ∀ d ∈ m, 3 ≤ d := by
    intro d hd
    have h := hm d hd
    simp [SimpleDims] at h
    omega
  have hall3 : ∀ d ∈ m, d = 3 := by
    intro d hd
    have hdm := hm d hd
    simp [SimpleDims] at hdm
    rcases Multiset.exists_cons_of_mem hd with ⟨t, rfl⟩
    have ht3 : ∀ x ∈ t, 3 ≤ x := by
      intro x hx
      exact h3 x (by simp [hx])
    have hbound := Multiset.card_nsmul_le_sum ht3
    simp only [smul_eq_mul, Multiset.sum_cons] at hbound hs
    rcases hdm with rfl | rfl | rfl
    · rfl
    · have hsum : t.sum = 4 := by omega
      have hcard : t.card ≤ 1 := by omega
      interval_cases hc : t.card
      · rw [Multiset.card_eq_zero] at hc
        simp [hc] at hsum
      · obtain ⟨a, ha⟩ := Multiset.card_eq_one.mp hc
        have haMem : a ∈ SimpleDims := hm a (by simp [ha])
        simp [SimpleDims] at haMem
        rw [ha] at hsum
        simp at hsum
        omega
    · have hsum : t.sum = 2 := by omega
      have hcard : t.card = 0 := by omega
      rw [Multiset.card_eq_zero] at hcard
      simp [hcard] at hsum
  have hmrep : m = Multiset.replicate m.card 3 :=
    Multiset.eq_replicate.mpr ⟨rfl, hall3⟩
  have hs' : m.card * 3 = 12 := by
    calc
      m.card * 3 = (Multiset.replicate m.card 3).sum := by simp
      _ = m.sum := by rw [← hmrep]
      _ = 12 := hs
  have hc : m.card = 4 := by omega
  rw [hmrep, hc]
  rfl

/-! ## S3: possible centre dimensions -/

/-- The subset-sum dimensions arising from irreducible dimensions
`{1, 3, 3, 5}`. -/
def SubmoduleDims : Finset ℕ := {0, 1, 3, 4, 5, 6, 7, 8, 9, 11, 12}

/-- In particular, the listed decomposition has no 2-dimensional submodule. -/
theorem two_not_submodule_dim : 2 ∉ SubmoduleDims := by decide

/-- The irreducible-dimension multiset contains exactly two 3-dimensional
summands. -/
theorem three_appears_twice :
    (Multiset.filter (· = 3) ({1, 3, 3, 5} : Multiset ℕ)).card = 2 := by decide

/-! ## S2: A₅ cannot permute at most four things nontrivially -/

/-- Every action of `A₅` on `Fin n`, for `n ≤ 4`, is trivial. -/
theorem action_trivial_of_card_le_four {n : ℕ} (hn : n ≤ 4)
    (f : alternatingGroup (Fin 5) →* Perm (Fin n)) (g : alternatingGroup (Fin 5)) :
    f g = 1 := by
  have hsimple : IsSimpleGroup (alternatingGroup (Fin 5)) :=
    alternatingGroup.isSimpleGroup_five
  rcases hsimple.eq_bot_or_eq_top_of_normal f.ker (MonoidHom.normal_ker f) with hk | hk
  · exfalso
    have hinj : Function.Injective f := (MonoidHom.ker_eq_bot_iff f).mp hk
    have hcard : Nat.card (alternatingGroup (Fin 5)) ≤ Nat.card (Perm (Fin n)) :=
      Nat.card_le_card_of_injective f hinj
    have h60 : Nat.card (alternatingGroup (Fin 5)) = 60 := by
      rw [nat_card_alternatingGroup]
      norm_num [Nat.factorial]
    have h24 : Nat.card (Perm (Fin n)) ≤ 24 := by
      rw [Nat.card_eq_fintype_card, Fintype.card_perm, Fintype.card_fin]
      interval_cases n <;> norm_num [Nat.factorial]
    omega
  · have hg : g ∈ f.ker := by
      rw [hk]
      trivial
    simpa [MonoidHom.mem_ker] using hg

/-- Consequently every index in a family of at most four factors is fixed. -/
theorem factors_are_submodules {n : ℕ} (hn : n ≤ 4)
    (f : alternatingGroup (Fin 5) →* Perm (Fin n)) (g : alternatingGroup (Fin 5))
    (i : Fin n) : f g i = i := by
  rw [action_trivial_of_card_le_four hn f g]
  rfl

end OPHSelection


/- ==================================================================== -/
/- PART V — deriving the central-triviality hypothesis                  -/
/- ==================================================================== -/

/-
  OPHInnerAction.lean

  PURPOSE.  In `OPHSelection.lean` the uniqueness of u(1) ⊕ su(2) ⊕ su(3) rested
  on a hypothesis found by working backwards from the answer:

      (H*)  A₅ acts trivially on the centre Z(g).

  This file replaces (H*) by a more primitive hypothesis and proves the parts of
  the reduction that Mathlib supports.

      (H2)  The screen symmetry is a gauge symmetry:
            A₅ ↪ G, with G compact connected and Lie(G) = g.

  THE REDUCTION (H2) ⟹ (H*).  Tags: [P] proved here, [C] classical input not
  in Mathlib.

    (1) [P]  ad x annihilates the centre.  (`ad_centre_eq_zero`, definitional.)
    (2) [C]  Hence exp(ad x) fixes Z(g) pointwise, so Inn(g) acts trivially
             on the centre.  Needs `exp` for Lie algebras.
    (3) [C]  G compact connected ⟹ exp surjective ⟹ Ad(G) = Inn(g).

    Therefore A₅ ⊂ G acts on g by Ad, which is inner, which fixes Z(g)
    pointwise.  (H*) is a consequence, not an assumption.

  ALSO PROVED HERE.
    [P]  The centre is characteristic (`centre_invariant`): this is what makes
         Z(g) an A₅-submodule of R₁₂, the step the enumeration in
         `OPHSelection.lean` relies on.
    [P]  A₅ does not embed in any group with a unique involution
         (`no_embedding_of_unique_involution`).  The only element of order 2 in
         SU(2) is −I, so A₅ ⊄ SU(2).

  WHY THAT LAST ONE MATTERS.  (H2) requires A₅ to sit inside the actual gauge
  group, so it forks:
      – either the weak factor appears in adjoint form SO(3) = SU(2)/Z₂,
      – or the screen group is the double cover SL(2,5) = 2·A₅, not A₅.
  This is a prediction of (H2) that could have gone the other way, which is what
  distinguishes it from a hypothesis fitted to the conclusion.  It connects to
  the verified result that A₅ has no Z₆ while SL(2,5) does, and to the fact that
  the Standard Model gauge group is (SU(3)×SU(2)×U(1))/Z₆.

  WHAT IS *NOT* ESTABLISHED.  That OPH supplies (H2).  The reduction is proved
  modulo the two classical inputs; whether the readback rule puts the screen
  symmetry *inside* the gauge group, rather than acting on the record structure
  from outside, is an OPH-internal question and is untouched here.  If the
  action is by outer automorphisms the reduction fails and (H*) reverts to
  being an assumption.

  STATUS: compiled and machine-checked.  No sorries.
-/



namespace OPHInner

/-! ## Part 1 — the centre, and why inner actions fix it

The centre is defined directly rather than via `LieAlgebra.center`, so that
step (1) is definitional.  `centre_eq_mathlib` records the agreement.
-/

section Centre

variable {R L : Type*} [CommRing R] [LieRing L] [LieAlgebra R L]

/-- The centre of a Lie algebra: elements killed by every `ad x`. -/
def centre (L : Type*) [LieRing L] : Set L := {z | ∀ x : L, ⁅x, z⁆ = 0}

@[simp] lemma mem_centre_iff {z : L} : z ∈ centre L ↔ ∀ x : L, ⁅x, z⁆ = 0 := Iff.rfl

/-- **Step (1).**  `ad x` annihilates the centre.  This is the entire reason an
inner automorphism acts trivially on `Z(L)`: it is `exp` of something that
kills the centre. -/
theorem ad_centre_eq_zero (x : L) {z : L} (hz : z ∈ centre L) : ⁅x, z⁆ = 0 := hz x

/-- Restated with Mathlib's `ad`. -/
theorem ad_apply_centre (x : L) {z : L} (hz : z ∈ centre L) :
    (LieAlgebra.ad R L x) z = 0 := by
  rw [LieAlgebra.ad_apply]
  exact hz x

theorem centre_zero_mem : (0 : L) ∈ centre L := by
  intro x
  simp

theorem centre_add_mem {z w : L} (hz : z ∈ centre L) (hw : w ∈ centre L) :
    z + w ∈ centre L := by
  intro x
  rw [lie_add, hz x, hw x, add_zero]

theorem centre_smul_mem (r : R) {z : L} (hz : z ∈ centre L) : r • z ∈ centre L := by
  intro x
  rw [lie_smul, hz x, smul_zero]

/-- **The centre is characteristic.**  Every Lie algebra automorphism maps the
centre into the centre.  Applied to the A₅-action, this is what makes `Z(g)` an
A₅-submodule of `R₁₂` — the step used by the enumeration in
`OPHSelection.lean`. -/
theorem centre_invariant (φ : L ≃ₗ⁅R⁆ L) {z : L} (hz : z ∈ centre L) :
    φ z ∈ centre L := by
  intro x
  have h : ⁅φ (φ.symm x), φ z⁆ = φ ⁅φ.symm x, z⁆ := (φ.map_lie _ _).symm
  rw [φ.apply_symm_apply] at h
  rw [h, hz (φ.symm x), map_zero]

/-- The property the reduction delivers: the group acts by automorphisms fixing
the centre pointwise.  Classically this holds whenever the action is by `Ad(G)`
with `G` compact connected — steps (2) and (3) of the header. -/
def CentrallyTrivial (G : Type*) [Group G] [MulAction G L] : Prop :=
  ∀ (g : G) (z : L), z ∈ centre L → g • z = z

/-- The direct definition of the centre agrees with Mathlib's Lie-algebra centre. -/
theorem centre_eq_mathlib : centre L = (LieAlgebra.center R L : Set L) := by
  ext z
  simp [mem_centre_iff, LieModule.mem_maxTrivSubmodule]

end Centre

/-! ## Part 2 — which group can realize the screen symmetry

A₅ has fifteen involutions; SU(2) has exactly one, namely `−I` (any `g ∈ SU(2)`
with `g² = I`, `g ≠ I` has both eigenvalues `−1`).  Two involutions suffice for
the obstruction, so we exhibit two.
-/

section Involutions

/-- `(0 1)(2 3)`. -/
def t₁ : Perm (Fin 5) := Equiv.swap 0 1 * Equiv.swap 2 3

/-- `(0 1)(2 4)`. -/
def t₂ : Perm (Fin 5) := Equiv.swap 0 1 * Equiv.swap 2 4

lemma sign_t₁ : Equiv.Perm.sign t₁ = 1 := by
  rw [t₁, map_mul,
    Equiv.Perm.sign_swap (show (0 : Fin 5) ≠ 1 by decide),
    Equiv.Perm.sign_swap (show (2 : Fin 5) ≠ 3 by decide)]
  decide

lemma sign_t₂ : Equiv.Perm.sign t₂ = 1 := by
  rw [t₂, map_mul,
    Equiv.Perm.sign_swap (show (0 : Fin 5) ≠ 1 by decide),
    Equiv.Perm.sign_swap (show (2 : Fin 5) ≠ 4 by decide)]
  decide

-- These equalities on `Perm (Fin 5)` are discharged by finite computation.
lemma t₁_sq : t₁ * t₁ = 1 := by decide
lemma t₂_sq : t₂ * t₂ = 1 := by decide
lemma t₁_ne_one : t₁ ≠ 1 := by decide
lemma t₂_ne_one : t₂ ≠ 1 := by decide
lemma t₁_ne_t₂ : t₁ ≠ t₂ := by decide

/-- `(0 1)(2 3)` as an element of `A₅`. -/
def s₁ : alternatingGroup (Fin 5) :=
  ⟨t₁, by rw [Equiv.Perm.mem_alternatingGroup]; exact sign_t₁⟩

/-- `(0 1)(2 4)` as an element of `A₅`. -/
def s₂ : alternatingGroup (Fin 5) :=
  ⟨t₂, by rw [Equiv.Perm.mem_alternatingGroup]; exact sign_t₂⟩

lemma s₁_sq : s₁ * s₁ = 1 := by
  apply Subtype.ext
  exact t₁_sq

lemma s₂_sq : s₂ * s₂ = 1 := by
  apply Subtype.ext
  exact t₂_sq
lemma s₁_ne_one : s₁ ≠ 1 := fun h => t₁_ne_one (congrArg Subtype.val h)
lemma s₂_ne_one : s₂ ≠ 1 := fun h => t₂_ne_one (congrArg Subtype.val h)
lemma s₁_ne_s₂ : s₁ ≠ s₂ := fun h => t₁_ne_t₂ (congrArg Subtype.val h)

/-- **A₅ does not embed in a group with a unique involution.**

In particular `A₅ ⊄ SU(2)`, since the only element of order 2 in SU(2) is `−I`.
Hence hypothesis (H2) forces either the adjoint form SO(3) or the double cover
SL(2,5) = 2·A₅. -/
theorem no_embedding_of_unique_involution {G : Type*} [Group G]
    (huniq : ∀ a b : G, a ≠ 1 → b ≠ 1 → a * a = 1 → b * b = 1 → a = b)
    (f : alternatingGroup (Fin 5) →* G) : ¬ Function.Injective f := by
  intro hf
  have e₁ : f s₁ * f s₁ = 1 := by rw [← map_mul, s₁_sq, map_one]
  have e₂ : f s₂ * f s₂ = 1 := by rw [← map_mul, s₂_sq, map_one]
  have h₁ : f s₁ ≠ 1 := fun h => s₁_ne_one (hf (by rw [h, map_one]))
  have h₂ : f s₂ ≠ 1 := fun h => s₂_ne_one (hf (by rw [h, map_one]))
  exact s₁_ne_s₂ (hf (huniq _ _ h₁ h₂ e₁ e₂))

/-- Contrapositive: if the screen symmetry embeds in `G`, then `G` has at least
two distinct involutions. -/
theorem two_involutions_of_embedding {G : Type*} [Group G]
    (f : alternatingGroup (Fin 5) →* G) (hf : Function.Injective f) :
    ∃ a b : G, a ≠ b ∧ a ≠ 1 ∧ b ≠ 1 ∧ a * a = 1 ∧ b * b = 1 := by
  refine ⟨f s₁, f s₂, ?_, ?_, ?_, ?_, ?_⟩
  · exact fun h => s₁_ne_s₂ (hf h)
  · exact fun h => s₁_ne_one (hf (by rw [h, map_one]))
  · exact fun h => s₂_ne_one (hf (by rw [h, map_one]))
  · rw [← map_mul, s₁_sq, map_one]
  · rw [← map_mul, s₂_sq, map_one]

end Involutions

end OPHInner



/- ==================================================================== -/
/- PART VI — trichotomy arithmetic, noncentrality, and the Z6 lattice   -/
/- ==================================================================== -/

/-
  OPHTrichotomy.lean

  The remaining finite, exact steps of the OPH compact-Lie trichotomy and the
  screen gluing class.  Each part below is one named step of a proof in the OPH
  SM/GR paper, section "Why is the gauge algebra 8+3+1?".

  ------------------------------------------------------------------------
  PART   PAPER STEP                                          STATUS
  ------------------------------------------------------------------------
  I      "no compact semisimple Lie algebra has dimension
          1, 2, 4, 5 or 7"  — used to kill the centre
          dimensions 11, 10, 8, 7 and 5 in the trichotomy    proved here
  II     "[iS, iT] = -2(E₁₂ - E₂₁) ≠ 0" — the rank-five
          band is noncentral for the constructed bracket     proved here
  III    "Λ₊/(Λ₁ ⊕ Λ₅) ≅ ℤ/6ℤ", with A₅-invariance and
          sign reversal under antipodal inversion            proved here
  ------------------------------------------------------------------------

  ALREADY PROVED in the companion file `A5_OPH.lean` (namespace OPHSelection):
    * `action_trivial_of_card_le_four` — A₅ cannot permute four simple ideals,
      the step that excludes su(2)^4 in the trichotomy proof;
    * `sum_eq_eleven`  — 11 = 3 + 8 uniquely  (centre dimension 1);
    * `sum_eq_twelve`  — 12 = 3+3+3+3 uniquely (centre dimension 0).

  Together, Parts I–III of this file and those three results cover every
  finite arithmetic and group-theoretic step of the trichotomy proof, plus the
  lattice computation behind the ℤ₆ quotient.

  WHAT IS *NOT* PROVED HERE, and is not claimed:
    * the classification of compact simple Lie algebras (the source of the
      dimension list {3, 8, 10, 14, ...}) — not in Mathlib, taken as input;
    * the rationality lemma (the centre is a ℚ-rational A₅-submodule because it
      is the Lie algebra of a torus with an automorphism-preserved integral
      cocharacter lattice) — this needs the Galois structure of the A₅ character
      field ℚ(√5) and is a separate project;
    * that the constructed bracket is the physical current bracket, that the
      A₅-action on it is inner, or that the six-axis lattice is a physical
      cocharacter lattice.  Those are the paper's own open receipts.

  Note on Part III: Smith invariants (1,1,1,1,1,6) for the presentation matrix
  are equivalent to the statement that the quotient is cyclic of order 6, which
  is exactly `screen_gluing_class` below.  Smith normal form itself is not
  computed; the isomorphism carries the same content.

  STATUS: compiled and kernel-checked with no unfinished proofs.
-/
namespace OPHTrichotomy

open Matrix

/-! ## Part I — dimensions carrying no compact semisimple Lie algebra

Compact simple Lie algebras of dimension at most 12 have dimension 3 (`su(2)`),
8 (`su(3)`) or 10 (`so(5)`); the next is 14 (`g₂`).  Taking that list as given,
the following dimensions admit no compact semisimple algebra at all.  In the
trichotomy these exclude centre dimensions 11, 10, 8, 7 and 5. -/

/-- Dimensions of compact simple Lie algebras that are at most 12. -/
def SimpleDims : Finset ℕ := {3, 8, 10}

private lemma three_le_of_mem {m : Multiset ℕ} (hm : ∀ d ∈ m, d ∈ SimpleDims)
    {d : ℕ} (hd : d ∈ m) : 3 ≤ d := by
  have h := hm d hd
  simp [SimpleDims] at h
  omega

/-- Every entry is at least 3, so the multiset is short relative to its sum. -/
private lemma card_bound {m : Multiset ℕ} (hm : ∀ d ∈ m, d ∈ SimpleDims) :
    3 * m.card ≤ m.sum := by
  have h3 : ∀ d ∈ m, 3 ≤ d := fun d hd => three_le_of_mem hm hd
  have h := Multiset.card_nsmul_le_sum h3
  simp only [smul_eq_mul] at h
  omega

/-- **No compact semisimple Lie algebra has dimension 1.** -/
theorem sum_ne_one (m : Multiset ℕ) (hm : ∀ d ∈ m, d ∈ SimpleDims) :
    m.sum ≠ 1 := by
  intro hs
  have hc := card_bound hm
  rw [hs] at hc
  have h0 : m.card = 0 := by omega
  rw [Multiset.card_eq_zero] at h0
  rw [h0] at hs
  simp at hs

/-- **No compact semisimple Lie algebra has dimension 2.** -/
theorem sum_ne_two (m : Multiset ℕ) (hm : ∀ d ∈ m, d ∈ SimpleDims) :
    m.sum ≠ 2 := by
  intro hs
  have hc := card_bound hm
  rw [hs] at hc
  have h0 : m.card = 0 := by omega
  rw [Multiset.card_eq_zero] at h0
  rw [h0] at hs
  simp at hs

/-- **No compact semisimple Lie algebra has dimension 4.** -/
theorem sum_ne_four (m : Multiset ℕ) (hm : ∀ d ∈ m, d ∈ SimpleDims) :
    m.sum ≠ 4 := by
  intro hs
  have hc := card_bound hm
  rw [hs] at hc
  have hcard : m.card ≤ 1 := by omega
  interval_cases h : m.card
  · rw [Multiset.card_eq_zero] at h
    rw [h] at hs
    simp at hs
  · obtain ⟨a, ha⟩ := Multiset.card_eq_one.mp h
    have haM : a ∈ SimpleDims := hm a (by rw [ha]; simp)
    simp [SimpleDims] at haM
    rw [ha] at hs
    simp at hs
    omega

/-- **No compact semisimple Lie algebra has dimension 5.** -/
theorem sum_ne_five (m : Multiset ℕ) (hm : ∀ d ∈ m, d ∈ SimpleDims) :
    m.sum ≠ 5 := by
  intro hs
  have hc := card_bound hm
  rw [hs] at hc
  have hcard : m.card ≤ 1 := by omega
  interval_cases h : m.card
  · rw [Multiset.card_eq_zero] at h
    rw [h] at hs
    simp at hs
  · obtain ⟨a, ha⟩ := Multiset.card_eq_one.mp h
    have haM : a ∈ SimpleDims := hm a (by rw [ha]; simp)
    simp [SimpleDims] at haM
    rw [ha] at hs
    simp at hs
    omega

/-- **No compact semisimple Lie algebra has dimension 7.**
This is the exclusion that kills centre dimension 5 in the trichotomy. -/
theorem sum_ne_seven (m : Multiset ℕ) (hm : ∀ d ∈ m, d ∈ SimpleDims) :
    m.sum ≠ 7 := by
  intro hs
  have hc := card_bound hm
  rw [hs] at hc
  have hcard : m.card ≤ 2 := by omega
  interval_cases h : m.card
  · rw [Multiset.card_eq_zero] at h
    rw [h] at hs
    simp at hs
  · obtain ⟨a, ha⟩ := Multiset.card_eq_one.mp h
    have haM : a ∈ SimpleDims := hm a (by rw [ha]; simp)
    simp [SimpleDims] at haM
    rw [ha] at hs
    simp at hs
    omega
  · obtain ⟨a, b, hab⟩ := Multiset.card_eq_two.mp h
    have haM : a ∈ SimpleDims := hm a (by rw [hab]; simp)
    have hbM : b ∈ SimpleDims := hm b (by rw [hab]; simp)
    simp [SimpleDims] at haM hbM
    rw [hab] at hs
    simp at hs
    omega

/-- Packaged form: the five excluded semisimple dimensions. -/
theorem sum_not_mem_excluded (m : Multiset ℕ) (hm : ∀ d ∈ m, d ∈ SimpleDims) :
    m.sum ≠ 1 ∧ m.sum ≠ 2 ∧ m.sum ≠ 4 ∧ m.sum ≠ 5 ∧ m.sum ≠ 7 :=
  ⟨sum_ne_one m hm, sum_ne_two m hm, sum_ne_four m hm,
   sum_ne_five m hm, sum_ne_seven m hm⟩

/-! ## Part II — the rank-five band is noncentral

For the constructed bracket on the coefficient space, the paper exhibits
`S = diag(1,-1,0)` in the 𝟓-band and `T = E₁₂ + E₂₁` in the 𝟑'-band with

    [iS, iT] = -2 (E₁₂ - E₂₁) ≠ 0.

This is the witness for the noncentral-quintet corollary: in the other two
outcomes of the trichotomy the rank-five band *is* central, so a nonzero
bracket selects the Standard-Model Lie type. -/

section Noncentral

/-- The 𝟓-band witness `S = diag(1, -1, 0)`. -/
def Sq : Matrix (Fin 3) (Fin 3) ℂ := !![1, 0, 0; 0, -1, 0; 0, 0, 0]

/-- The 𝟑'-band witness `T = E₁₂ + E₂₁`. -/
def Tq : Matrix (Fin 3) (Fin 3) ℂ := !![0, 1, 0; 1, 0, 0; 0, 0, 0]

/-- The `(0,1)` entry of `[S, T]` is `2`; hence `[S,T] = 2(E₁₂ - E₂₁)`. -/
lemma bracket_entry_01 : (⁅Sq, Tq⁆ : Matrix (Fin 3) (Fin 3) ℂ) 0 1 = 2 := by
  rw [Ring.lie_def]
  simp [Sq, Tq]
  norm_num

theorem bracket_ne_zero : (⁅Sq, Tq⁆ : Matrix (Fin 3) (Fin 3) ℂ) ≠ 0 := by
  intro h
  have h01 := bracket_entry_01
  rw [h] at h01
  simp at h01

/-- **The rank-five band is noncentral.**  `[iS, iT] = -[S, T] ≠ 0`. -/
theorem quintet_noncentral :
    ⁅(Complex.I • Sq), (Complex.I • Tq)⁆ ≠ 0 := by
  have hI : ⁅(Complex.I • Sq), (Complex.I • Tq)⁆ = -⁅Sq, Tq⁆ := by
    rw [smul_lie, lie_smul, smul_smul, Complex.I_mul_I, neg_one_smul]
  rw [hI, neg_ne_zero]
  exact bracket_ne_zero

end Noncentral

/-! ## Part III — the screen gluing class `Λ₊ / (Λ₁ ⊕ Λ₅) ≅ ℤ/6`

Antipodal port pairing gives the even integral load lattice `Λ₊ ≅ ℤ⁶`, with
`Λ₁` the multiples of the all-ones vector and `Λ₅` the sum-zero sublattice.
The total-load map modulo 6 is surjective with kernel exactly `Λ₁ ⊔ Λ₅`. -/

section Lattice

/-- The all-ones port vector. -/
def ones : Fin 6 → ℤ := fun _ => 1

/-- Total port load. -/
def sigma : (Fin 6 → ℤ) →+ ℤ where
  toFun x := ∑ i, x i
  map_zero' := by simp
  map_add' := by intro x y; simp [Finset.sum_add_distrib]

@[simp] lemma sigma_apply (x : Fin 6 → ℤ) : sigma x = ∑ i, x i := rfl

@[simp] lemma sigma_ones : sigma ones = 6 := by
  simp [sigma, ones]

/-- `Λ₁` : integer multiples of the all-ones vector. -/
def Λ₁ : AddSubgroup (Fin 6 → ℤ) := AddSubgroup.zmultiples ones

/-- `Λ₅` : the centred (sum-zero) sublattice. -/
def Λ₅ : AddSubgroup (Fin 6 → ℤ) := sigma.ker

/-- The screen gluing-class map `q(x) = Σᵢ xᵢ mod 6`. -/
def q : (Fin 6 → ℤ) →+ ZMod 6 := (Int.castAddHom (ZMod 6)).comp sigma

@[simp] lemma q_apply (x : Fin 6 → ℤ) : q x = ((∑ i, x i : ℤ) : ZMod 6) := rfl

/-- `q` is surjective; the witness is `e₁`, as in the paper. -/
theorem q_surjective : Function.Surjective q := by
  intro z
  obtain ⟨n, rfl⟩ := ZMod.intCast_surjective z
  refine ⟨fun i => if i = 0 then n else 0, ?_⟩
  simp [q, sigma, Finset.sum_ite_eq']

/-- **The kernel of the gluing class is exactly `Λ₁ ⊔ Λ₅`.** -/
theorem sup_eq_ker : Λ₁ ⊔ Λ₅ = q.ker := by
  apply le_antisymm
  · apply sup_le
    · rintro _ ⟨m, rfl⟩
      -- σ(m • ones) = 6m ≡ 0 (mod 6)
      simp [AddMonoidHom.mem_ker, q, ones]
      rw [show (6 : ZMod 6) = 0 by exact ZMod.natCast_self 6]
      simp
    · intro x hx
      simp only [Λ₅, AddMonoidHom.mem_ker] at hx
      simp [AddMonoidHom.mem_ker, q, hx]
  · intro x hx
    simp only [AddMonoidHom.mem_ker, q, AddMonoidHom.coe_comp,
      Function.comp_apply, Int.coe_castAddHom] at hx
    -- 6 ∣ σ x
    obtain ⟨m, hm⟩ := (ZMod.intCast_zmod_eq_zero_iff_dvd _ 6).mp hx
    rw [AddSubgroup.mem_sup]
    refine ⟨m • ones, ⟨m, rfl⟩, x - m • ones, ?_, by abel⟩
    simp only [Λ₅, AddMonoidHom.mem_ker, map_sub, map_zsmul, sigma_ones, hm]
    simp [mul_comm]

/-- **Screen gluing class.**  `Λ₊ / (Λ₁ ⊕ Λ₅) ≅ ℤ/6ℤ`.
Equivalently, the presentation matrix has Smith invariants `(1,1,1,1,1,6)`. -/
theorem screen_gluing_class :
    Nonempty (((Fin 6 → ℤ) ⧸ (Λ₁ ⊔ Λ₅)) ≃+ ZMod 6) := by
  rw [sup_eq_ker]
  exact ⟨QuotientAddGroup.quotientKerEquivOfSurjective q q_surjective⟩

/-- Proper icosahedral permutations act trivially on the gluing class:
any coordinate permutation preserves the total load. -/
theorem q_perm_invariant (e : Equiv.Perm (Fin 6)) (x : Fin 6 → ℤ) :
    q (x ∘ e) = q x := by
  simp only [q_apply, Function.comp_apply]
  congr 1
  exact Fintype.sum_equiv e _ _ (fun i => rfl)

/-- Signed antipodal reversal negates the gluing class. -/
theorem q_neg (x : Fin 6 → ℤ) : q (-x) = -q x := map_neg q x

end Lattice

end OPHTrichotomy

/- ==================================================================== -/
/- AXIOM AUDIT                                                          -/
/- Expected output for every line: 'propext, Classical.choice, Quot.sound'
   or a subset.  Any other axiom — in particular Lean.ofReduceBool from
   native_decide — means the build is not trustworthy.                   -/
/- ==================================================================== -/

#print axioms A5Order6.exists_subgroup_card_six
#print axioms A5Order6.H_not_commutative
#print axioms A5NoZ6.no_orderOf_six
#print axioms A5NoZ6.no_cyclic_six
#print axioms A5NoZ6.no_embedding_zmod_six
#print axioms OPHGap.gModule_iso_does_not_determine_bracket
#print axioms OPHSelection.sum_eq_eleven
#print axioms OPHSelection.sum_eq_twelve
#print axioms OPHSelection.action_trivial_of_card_le_four
#print axioms OPHInner.centre_invariant
#print axioms OPHInner.no_embedding_of_unique_involution
#print axioms OPHInner.two_involutions_of_embedding
#print axioms OPHTrichotomy.sum_ne_seven
#print axioms OPHTrichotomy.sum_not_mem_excluded
#print axioms OPHTrichotomy.quintet_noncentral
#print axioms OPHTrichotomy.sup_eq_ker
#print axioms OPHTrichotomy.screen_gluing_class
#print axioms OPHTrichotomy.q_perm_invariant
