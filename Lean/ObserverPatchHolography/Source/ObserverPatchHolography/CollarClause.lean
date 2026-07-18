import Mathlib
import ObserverPatchHolography.Primitives

/-!
# Issue #544 (central-interface collar clause): what the overlap-consistency
# layer can and cannot force — a machine-checked layer separation

Issue #544 asks whether overlap-consistent repair on the declared
fixed-cutoff branch *forces* the central-interface collar clause (the
declared-branch input stated with Axiom 3 of the compact paper,
`par:cicclause`): non-central cross-cut couplings excluded, not merely not
selected. The sharp failure mode is pinned in
`code/collar_alignment/test_msa_characterizations.py::test_descent_invariant_but_noncentral_interface_breaks_alignment`:
a `K`-invariant but **non-central** cross-cut coupling that respects the
boundary symmetry (so it changes no constraint) yet breaks entropic
alignment. As that test's docstring puts it: *invariance alone is not the
clause; centrality is.*

## What this module proves (the layer separation)

This module makes one precise, admission-free contribution to that
question: **the overlap-consistency repair layer, as formalised in
`Primitives.lean`, factors through the realized constraint family alone.**
Fix a repair base (patches, per-patch state spaces, edges, endpoints) and
vary only the interface package (alphabets, projections, weights,
distances). If two packages induce the same edge-consistency predicate on
every record (`SameConstraintFamily`), then:

* the firing trigger, single-site solvability, the constructed local move
  and the composite `Repair` operator are **identical**
  (`SameConstraintFamily.localRepair_eq`, `SameConstraintFamily.Repair_eq`
  — on-the-nose equality of operators, not mere similarity);
* the accepted-step relation, normal forms, broken-edge count and the
  consistent set are identical
  (`acceptedStep_iff`, `normalForm_iff`, `mismatchCount_eq`,
  `consistent_iff`);
* the abstract `H2`/`H3` local-law *statements* of the
  `LocalRepairDynamics` layer transfer verbatim (`H2_iff`, `H3_iff`; `H1`
  does not mention interfaces at all).

Meanwhile the interface **energetics** are not determined by the constraint
family: `fluxKit`/`coupledKit` below share one constraint family (so their
repair operators are *equal*, `demoKits_repair_eq`) yet have different
mismatch potentials `Φ` on the same broken record (`demoKits_phi_ne`) — the
carrier-level echo of the pinned test, where the extra coupling moves the
weights (alignment functional) while respecting every constraint.

## Consequence for #544, stated honestly

Within the formalised overlap-consistency layer, the collar clause is not
forceable **because it is not expressible**: the layer's acceptance data
(edge-consistency) cannot distinguish a central interface energy from an
invariant-but-non-central one when both realize the same constraint family.
Any derivation of the clause must therefore consume structure *beyond* this
layer — the Axiom-3 density family, the `C*(K̂_Σ)` center and flux functions
`π_L(Z(C*(K̂_Σ)))`, the touched-overlap acceptance contract's algebraic
content, or refinement closure — exactly the machinery in which #543 located
the clause.

## Declared scope (what this does NOT do)

This is a **layer result, not a resolution of #544**. It does not derive
the clause (that would need the operator-algebraic layer formalised), and it
does not prove independence for the paper's admissible class (the paper's
axiom system is strictly richer than the Lean core: Axiom-3 densities,
refinement closure and the collar C*-machinery are not formalised here).
The corpus labeling therefore stands unchanged: the clause remains an
explicit axiom-level declared-branch input (`par:cicclause`), its open
status recorded in `rem:msascope` and the D2 ledger rows — this module
sharpens *where* a derivation would have to live, not whether one exists.
No `sorry`, no `native_decide`, no new axiom.
-/

namespace OPH

/-! ## Interface kits over a fixed repair base -/

/-- An interface package over a fixed repair base (patches `P`, state spaces
    `S`, edges `E`, endpoints `src`/`tgt`): the interface alphabets,
    projections, weights and distances — everything an `OPHCarrier` adds on
    top of the base. Two kits over one base model "the same system with
    different interface structure/energetics". -/
structure InterfaceKit (P : Type) (S : P → Type) (E : Type)
    (src tgt : E → P) where
  Iface : E → Type
  projSrc : (e : E) → S (src e) → Iface e
  projTgt : (e : E) → S (tgt e) → Iface e
  weight : E → NNReal
  dist : (e : E) → Iface e → Iface e → NNReal
  weight_pos : ∀ e : E, 0 < weight e
  dist_eq_zero : ∀ (e : E) (a b : Iface e), dist e a b = 0 ↔ a = b

variable {P : Type} [Fintype P] [DecidableEq P] {S : P → Type}
  {E : Type} [Fintype E] {src tgt : E → P}

/-- The carrier assembled from the base and an interface kit. Reducible so
    that `Records`/`Site`/`Edge` of two kits over one base are transparently
    the same types. -/
@[reducible] def InterfaceKit.carrier (K : InterfaceKit P S E src tgt) :
    OPHCarrier where
  Patch := P
  State := S
  Edge := E
  src := src
  tgt := tgt
  Iface := K.Iface
  projSrc := K.projSrc
  projTgt := K.projTgt
  weight := K.weight
  dist := K.dist
  weight_pos := K.weight_pos
  dist_eq_zero := K.dist_eq_zero

variable (K₁ K₂ : InterfaceKit P S E src tgt)

/-- Two interface kits realize the **same constraint family** when they
    induce the same edge-consistency predicate on every record. This is the
    formal counterpart of a coupling change that "respects the boundary
    symmetry": it moves interface structure and energetics without moving a
    single constraint. -/
def SameConstraintFamily : Prop :=
  ∀ (e : E) (x : Records K₁.carrier),
    edgeConsistentAt (C := K₁.carrier) e x ↔ edgeConsistentAt (C := K₂.carrier) e x

namespace SameConstraintFamily

variable {K₁ K₂}

/-- The firing trigger is a function of the constraint family. -/
theorem trigger_iff (h : SameConstraintFamily K₁ K₂) (i : P)
    (x : Records K₁.carrier) :
    LocalTrigger (C := K₁.carrier) i x ↔ LocalTrigger (C := K₂.carrier) i x :=
  exists_congr fun e => and_congr Iff.rfl (not_congr (h e x))

/-- Single-site transactional solvability is a function of the constraint
    family (pointwise in the candidate repair). -/
theorem solvesAt_iff (h : SameConstraintFamily K₁ K₂) (i : P)
    (x : Records K₁.carrier) (s : S i) :
    SolvesAt (C := K₁.carrier) i x s ↔ SolvesAt (C := K₂.carrier) i x s :=
  forall_congr' fun e => imp_congr Iff.rfl (h e (Function.update x i s))

/-- Local satisfiability is a function of the constraint family. -/
theorem locallySolvable_iff (h : SameConstraintFamily K₁ K₂) (i : P)
    (x : Records K₁.carrier) :
    LocallySolvable (C := K₁.carrier) i x ↔ LocallySolvable (C := K₂.carrier) i x :=
  exists_congr (h.solvesAt_iff i x)

private theorem choose_eq_of_pred_iff' {α : Sort*} {p q : α → Prop}
    (hpq : ∀ a, p a ↔ q a) (hp : ∃ a, p a) (hq : ∃ a, q a) :
    Classical.choose hp = Classical.choose hq := by
  have hh : p = q := funext fun a => propext (hpq a)
  subst hh
  rfl

/-- **The constructed local move is determined by the constraint family** —
    equal operators, not merely equivalent behaviour. The trigger, the
    solvability predicate and the classically chosen repaired state all read
    only edge-consistency data, which the two kits share. -/
theorem localRepair_eq (h : SameConstraintFamily K₁ K₂) :
    localRepair K₁.carrier = localRepair K₂.carrier := by
  funext i x
  by_cases h₁ : LocalTrigger (C := K₁.carrier) i x ∧ LocallySolvable (C := K₁.carrier) i x
  · have h₂ : LocalTrigger (C := K₂.carrier) i x ∧ LocallySolvable (C := K₂.carrier) i x :=
      ⟨(h.trigger_iff i x).mp h₁.1, (h.locallySolvable_iff i x).mp h₁.2⟩
    rw [localRepair_of_fire K₁.carrier i x h₁, localRepair_of_fire K₂.carrier i x h₂]
    exact congrArg (Function.update x i)
      (choose_eq_of_pred_iff' (h.solvesAt_iff i x) h₁.2 h₂.2)
  · have h₂ : ¬ (LocalTrigger (C := K₂.carrier) i x ∧ LocallySolvable (C := K₂.carrier) i x) :=
      fun hc => h₁ ⟨(h.trigger_iff i x).mpr hc.1, (h.locallySolvable_iff i x).mpr hc.2⟩
    rw [localRepair_of_quiescent K₁.carrier i x h₁, localRepair_of_quiescent K₂.carrier i x h₂]

/-- The broken-edge set is determined by the constraint family (even though
    the *distances* defining it differ, their zero loci agree). -/
theorem brokenSet_eq (h : SameConstraintFamily K₁ K₂) (x : Records K₁.carrier) :
    brokenSet (C := K₁.carrier) x = brokenSet (C := K₂.carrier) x := by
  ext e
  rw [mem_brokenSet_iff_not_consistent, mem_brokenSet_iff_not_consistent]
  exact not_congr (h e x)

/-- The broken-edge count is determined by the constraint family. -/
theorem mismatchCount_eq (h : SameConstraintFamily K₁ K₂) (x : Records K₁.carrier) :
    mismatchCount (C := K₁.carrier) x = mismatchCount (C := K₂.carrier) x :=
  congrArg Finset.card (h.brokenSet_eq x)

/-- The accepted-step relation is determined by the constraint family. -/
theorem acceptedStep_iff (h : SameConstraintFamily K₁ K₂)
    (x y : Records K₁.carrier) :
    acceptedStep K₁.carrier x y ↔ acceptedStep K₂.carrier x y :=
  exists_congr fun i => by rw [h.localRepair_eq]; exact Iff.rfl

/-- Normal forms are determined by the constraint family. -/
theorem normalForm_iff (h : SameConstraintFamily K₁ K₂) (x : Records K₁.carrier) :
    NormalForm K₁.carrier x ↔ NormalForm K₂.carrier x :=
  forall_congr' fun y => not_congr (h.acceptedStep_iff x y)

/-- The consistent set is determined by the constraint family. -/
theorem consistent_iff (h : SameConstraintFamily K₁ K₂) (x : Records K₁.carrier) :
    Consistent K₁.carrier x ↔ Consistent K₂.carrier x := by
  rw [consistent_iff_edgeConsistent, consistent_iff_edgeConsistent]
  exact forall_congr' fun e => h e x

private theorem Repair_eq_aux (h : SameConstraintFamily K₁ K₂) :
    ∀ (n : ℕ) (x : Records K₁.carrier), mismatchCount (C := K₁.carrier) x = n →
      Repair K₁.carrier x = Repair K₂.carrier x := by
  intro n
  induction' n using Nat.strong_induction_on with n ih
  intro x hn
  by_cases hx : ∃ i : Site K₁.carrier, localRepair K₁.carrier i x ≠ x
  · have hx' : ∃ i : Site K₂.carrier, localRepair K₂.carrier i x ≠ x := by
      obtain ⟨i, hi⟩ := hx
      exact ⟨i, by rw [← h.localRepair_eq]; exact hi⟩
    have hsite : Classical.choose hx = Classical.choose hx' :=
      choose_eq_of_pred_iff' (fun i => by rw [h.localRepair_eq]; exact Iff.rfl) hx hx'
    rw [Repair_of_fire K₁.carrier x hx, Repair_of_fire K₂.carrier x hx', ← hsite,
      ← h.localRepair_eq]
    exact ih (mismatchCount (C := K₁.carrier)
        (localRepair K₁.carrier (Classical.choose hx) x))
      (hn ▸ mismatchCount_localRepair_lt K₁.carrier (Classical.choose hx) x
        (Classical.choose_spec hx))
      _ rfl
  · have hx' : ¬ ∃ i : Site K₂.carrier, localRepair K₂.carrier i x ≠ x := fun hc => by
      obtain ⟨i, hi⟩ := hc
      exact hx ⟨i, by rw [h.localRepair_eq]; exact hi⟩
    rw [Repair_of_normal K₁.carrier x hx, Repair_of_normal K₂.carrier x hx']

/-- **The composite repair operator is determined by the constraint
    family.** Two interface packages realizing the same constraints yield
    the *same* `Repair` — the overlap-consistency dynamics cannot see any
    interface structure beyond the constraints, central or otherwise. -/
theorem Repair_eq (h : SameConstraintFamily K₁ K₂) :
    Repair K₁.carrier = Repair K₂.carrier :=
  funext fun x => Repair_eq_aux h (mismatchCount (C := K₁.carrier) x) x rfl

/-- The abstract `H2` trigger law of the `LocalRepairDynamics` layer
    transfers verbatim between kits with one constraint family. (`H1` does
    not mention interfaces at all, so it is literally the same statement.) -/
theorem H2_iff (h : SameConstraintFamily K₁ K₂)
    (lr : P → Records K₁.carrier → Records K₁.carrier) :
    (∀ (i : P) (x : Records K₁.carrier),
        lr i x ≠ x ↔
          ∃ e : E, (src e = i ∨ tgt e = i) ∧ ¬ edgeConsistentAt (C := K₁.carrier) e x) ↔
    (∀ (i : P) (x : Records K₁.carrier),
        lr i x ≠ x ↔
          ∃ e : E, (src e = i ∨ tgt e = i) ∧ ¬ edgeConsistentAt (C := K₂.carrier) e x) :=
  forall_congr' fun _ => forall_congr' fun x => iff_congr Iff.rfl
    (exists_congr fun e => and_congr Iff.rfl (not_congr (h e x)))

/-- The abstract `H3` repair law transfers verbatim between kits with one
    constraint family. -/
theorem H3_iff (h : SameConstraintFamily K₁ K₂)
    (lr : P → Records K₁.carrier → Records K₁.carrier) :
    (∀ (i : P) (x : Records K₁.carrier),
        lr i x ≠ x → ∀ e : E, (src e = i ∨ tgt e = i) →
          edgeConsistentAt (C := K₁.carrier) e (lr i x)) ↔
    (∀ (i : P) (x : Records K₁.carrier),
        lr i x ≠ x → ∀ e : E, (src e = i ∨ tgt e = i) →
          edgeConsistentAt (C := K₂.carrier) e (lr i x)) :=
  forall_congr' fun i => forall_congr' fun x => imp_congr Iff.rfl
    (forall_congr' fun e => imp_congr Iff.rfl (h e (lr i x)))

end SameConstraintFamily

/-! ## The pinned failure mode, echoed at carrier level

The test suite's sharp counterexample is a coupling that respects every
constraint (boundary symmetry) but moves the energetics off the flux
sectors. Its carrier-level echo: `fluxKit` charges a broken overlap unit
cost (interface energy a function of the "flux" mismatch alone), while
`coupledKit` adds an extra symmetric cross-cut cost on the same broken
configurations. Same constraint family — provably the same `Repair`
operator — but different mismatch potential `Φ` on the same broken record.
The repair layer is blind to precisely the structure the collar clause
speaks about. -/

/-- Flux-only interface energy on a two-patch, one-edge base: unit cost on
    the broken overlap. -/
def fluxKit : InterfaceKit Bool (fun _ => Bool) Unit (fun _ => false) (fun _ => true) where
  Iface := fun _ => Bool
  projSrc := fun _ s => s
  projTgt := fun _ s => s
  weight := fun _ => 1
  dist := fun _ a b => if a = b then 0 else 1
  weight_pos := fun _ => one_pos
  dist_eq_zero := by
    intro _ a b
    by_cases hab : a = b
    · rw [if_pos hab]; exact ⟨fun _ => hab, fun _ => rfl⟩
    · rw [if_neg hab]; exact ⟨fun h1 => absurd h1 one_ne_zero, fun h2 => absurd h2 hab⟩

/-- The same constraint family carrying an extra symmetric cross-cut cost:
    every broken configuration is charged double. Constraint-preserving
    (same zero locus), energetics-shifting — the carrier-level analogue of
    the invariant-but-non-central coupling class. -/
def coupledKit : InterfaceKit Bool (fun _ => Bool) Unit (fun _ => false) (fun _ => true) where
  Iface := fun _ => Bool
  projSrc := fun _ s => s
  projTgt := fun _ s => s
  weight := fun _ => 1
  dist := fun _ a b => if a = b then 0 else 2
  weight_pos := fun _ => one_pos
  dist_eq_zero := by
    intro _ a b
    by_cases hab : a = b
    · rw [if_pos hab]; exact ⟨fun _ => hab, fun _ => rfl⟩
    · rw [if_neg hab]
      exact ⟨fun h1 => absurd h1 (by norm_num), fun h2 => absurd h2 hab⟩

/-- The two kits realize the same constraint family: their projections are
    identical, so edge-consistency is literally the same predicate. -/
theorem demoKits_same : SameConstraintFamily fluxKit coupledKit :=
  fun _ _ => Iff.rfl

/-- **The repair operators are equal on the nose** — the extra coupling is
    invisible to overlap-consistent repair. -/
theorem demoKits_repair_eq : Repair fluxKit.carrier = Repair coupledKit.carrier :=
  demoKits_same.Repair_eq

/-- **The energetics differ where the repair cannot look**: on the broken
    identity record the two mismatch potentials disagree (`1 ≠ 2`). Together
    with `demoKits_repair_eq` this is the layer separation in one concrete
    instance: constraint-respecting coupling changes `Φ` (the alignment
    functional's analogue) while every repair-layer notion stays fixed. -/
theorem demoKits_phi_ne :
    Φ fluxKit.carrier (fun b => b) ≠ Φ coupledKit.carrier (fun b => b) := by
  have h1 : Φ fluxKit.carrier (fun b => b) = 1 := by
    unfold Φ
    rw [Fintype.sum_unique]
    show (1 : NNReal) * (if (false : Bool) = true then 0 else 1) = 1
    rw [if_neg (by decide : ¬ (false : Bool) = true), mul_one]
  have h2 : Φ coupledKit.carrier (fun b => b) = 2 := by
    unfold Φ
    rw [Fintype.sum_unique]
    show (1 : NNReal) * (if (false : Bool) = true then 0 else 2) = 2
    rw [if_neg (by decide : ¬ (false : Bool) = true), one_mul]
  rw [h1, h2]
  norm_num

/-! ### Axiom audit — the layer-separation results are admission-free. -/
#print axioms SameConstraintFamily.localRepair_eq
#print axioms SameConstraintFamily.Repair_eq
#print axioms SameConstraintFamily.acceptedStep_iff
#print axioms SameConstraintFamily.normalForm_iff
#print axioms SameConstraintFamily.consistent_iff
#print axioms SameConstraintFamily.mismatchCount_eq
#print axioms SameConstraintFamily.H2_iff
#print axioms SameConstraintFamily.H3_iff
#print axioms demoKits_same
#print axioms demoKits_repair_eq
#print axioms demoKits_phi_ne

end OPH
