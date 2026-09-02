import Mathlib.Order.Interval.Finset.Basic
import ObserverPatchHolography.Provenance.FiniteCausetCompiler

/-!
# Finite cover-relation compiler into authenticated provenance

This module strengthens `FiniteCausetCompiler` at the direct-edge boundary.
For a finite decidable irreflexive transitive relation `R`, the existing
compiler can authenticate every strict `R`-predecessor as a direct
`ParentEdge`. Here the authenticated direct edges are instead exactly the
Hasse/cover relation of `R`, while generated precedence remains exactly `R`.

The proof separates two claims:

* finite strict order is the transitive closure of its cover relation; and
* the existing provenance compiler realizes that cover relation exactly at
  the authenticated direct-parent level.

This matches the logical shape used by the geometry-seeded causal-set control:
only cover links need be encoded as authenticated reads; transitive causal
order is reconstructed from them. It is still an expressivity theorem for a
supplied relation. It does not derive `R`, select a physical event population,
construct a threaded OPH execution, supply a count--volume law, or establish
manifoldlikeness or a continuum limit.
-/

namespace OPH.Provenance.FiniteCausetCompiler

universe u

variable {Event : Type u} [Fintype Event] [DecidableEq Event]
variable (R : Event → Event → Prop) [DecidableRel R]

/-- The strict cover relation of `R`: `b` covers `a` when `R a b` and there
is no `R`-intermediate event. -/
def CoverRelation (a b : Event) : Prop :=
  R a b ∧ ∀ ⦃c : Event⦄, R a c → ¬ R c b

instance coverRelationDecidable : DecidableRel (CoverRelation R) :=
  fun _ _ => by
    unfold CoverRelation
    infer_instance

omit [Fintype Event] [DecidableEq Event] [DecidableRel R] in
/-- A cover edge is in particular an edge of the supplied strict relation. -/
theorem coverRelation_rel {a b : Event} (h : CoverRelation R a b) : R a b :=
  h.1

/-- Every finite decidable irreflexive transitive relation is exactly the
transitive closure of its cover relation.

A temporary partial-order structure with `a ≤ b` defined as `a = b ∨ R a b`
lets us use Mathlib's locally-finite theorem `lt_iff_transGen_covBy`; the
result is then translated back to the relation-only interface used by OPH. -/
theorem rel_iff_transGen_coverRelation
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c)
    (a b : Event) :
    R a b ↔ Relation.TransGen (CoverRelation R) a b := by
  letI : PartialOrder Event := {
    le a b := a = b ∨ R a b
    le_refl a := Or.inl rfl
    le_trans a b c hab hbc := by
      rcases hab with rfl | hab
      · exact hbc
      rcases hbc with rfl | hbc
      · exact Or.inr hab
      · exact Or.inr (htrans hab hbc)
    le_antisymm a b hab hba := by
      rcases hab with hab | hab
      · exact hab
      rcases hba with hba | hba
      · exact hba.symm
      · exact False.elim (hirrefl a (htrans hab hba))
  }
  letI : DecidableLE Event := fun a b => by
    change Decidable (a = b ∨ R a b)
    infer_instance
  letI : DecidableLT Event := decidableLTOfDecidableLE
  letI : LocallyFiniteOrder Event := Fintype.toLocallyFiniteOrder

  have hlt : ∀ x y : Event, x < y ↔ R x y := by
    intro x y
    rw [lt_iff_le_not_ge]
    change ((x = y ∨ R x y) ∧ ¬ (y = x ∨ R y x)) ↔ R x y
    constructor
    · rintro ⟨hxy, hyx⟩
      rcases hxy with hxy | hxy
      · subst y
        exact False.elim (hyx (Or.inl rfl))
      · exact hxy
    · intro hxy
      refine ⟨Or.inr hxy, ?_⟩
      rintro (hyx | hyx)
      · subst y
        exact hirrefl x hxy
      · exact hirrefl x (htrans hxy hyx)

  have hcov : ∀ x y : Event, x ⋖ y ↔ CoverRelation R x y := by
    intro x y
    constructor
    · intro h
      refine ⟨(hlt x y).mp h.1, ?_⟩
      intro z hxz hzy
      exact h.2 ((hlt x z).mpr hxz) ((hlt z y).mpr hzy)
    · rintro ⟨hxy, hno⟩
      refine ⟨(hlt x y).mpr hxy, ?_⟩
      intro z hxz hzy
      exact hno ((hlt x z).mp hxz) ((hlt z y).mp hzy)

  constructor
  · intro hab
    have hchain : Relation.TransGen (fun x y : Event => x ⋖ y) a b :=
      (lt_iff_transGen_covBy).mp ((hlt a b).mpr hab)
    clear hab
    induction hchain with
    | single hcover =>
        exact Relation.TransGen.single ((hcov _ _).mp hcover)
    | tail _ hcover ih =>
        exact Relation.TransGen.tail ih ((hcov _ _).mp hcover)
  · intro hchain
    induction hchain with
    | single hcover =>
        exact hcover.1
    | tail _ hcover ih =>
        exact htrans ih hcover.1

omit [DecidableEq Event] in
/-- The predecessor-cardinality rank of the full relation increases on every
cover edge, so the cover compiler needs no independently supplied rank. -/
theorem predecessorRank_lt_of_cover
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c)
    {a b : Event} (h : CoverRelation R a b) :
    predecessorRank R a < predecessorRank R b :=
  predecessorRank_lt_of_rel R hirrefl htrans h.1

/-- Compile only the cover edges of a finite strict transitive relation into
an authenticated semantic log. The well-founded rank is derived from the
full relation's predecessor cardinality; it is not additional input. -/
def coverLog
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c) :
    SemanticEventLog Event Unit Event :=
  log (CoverRelation R) (predecessorRank R)
    (fun h => predecessorRank_lt_of_cover R hirrefl htrans h)

/-- Direct authenticated parenthood in `coverLog` is exactly the Hasse/cover
relation, not the full transitive order. -/
theorem coverLog_parentEdge_iff
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c)
    (a b : Event) :
    (coverLog R hirrefl htrans).ParentEdge a b ↔ CoverRelation R a b :=
  parentEdge_iff (CoverRelation R) (predecessorRank R)
    (fun h => predecessorRank_lt_of_cover R hirrefl htrans h) a b

/-- Transitive closure of the authenticated cover edges reconstructs exactly
the supplied strict relation. -/
theorem coverLog_generatedBefore_iff
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c)
    (a b : Event) :
    (coverLog R hirrefl htrans).GeneratedBefore a b ↔ R a b :=
  (generatedBefore_iff_transGen
      (CoverRelation R) (predecessorRank R)
      (fun h => predecessorRank_lt_of_cover R hirrefl htrans h) a b).trans
    (rel_iff_transGen_coverRelation R hirrefl htrans a b).symm

/-- Bundled strengthening of finite-poset grammar expressivity: every finite
decidable strict transitive relation has an authenticated semantic log whose
direct parents are exactly its covers and whose generated precedence is
exactly the original relation. This remains an abstract-log existence receipt,
not a source-selection or execution-custody theorem. -/
theorem finiteStrictTransitiveRelation_cover_realized
    (hirrefl : ∀ a : Event, ¬ R a a)
    (htrans : ∀ {a b c : Event}, R a b → R b c → R a c) :
    ∃ L : SemanticEventLog Event Unit Event,
      (∀ a b, L.ParentEdge a b ↔ CoverRelation R a b) ∧
      (∀ a b, L.GeneratedBefore a b ↔ R a b) := by
  exact ⟨coverLog R hirrefl htrans,
    coverLog_parentEdge_iff R hirrefl htrans,
    coverLog_generatedBefore_iff R hirrefl htrans⟩

/-! ## Axiom audit -/

#print axioms rel_iff_transGen_coverRelation
#print axioms predecessorRank_lt_of_cover
#print axioms coverLog
#print axioms coverLog_parentEdge_iff
#print axioms coverLog_generatedBefore_iff
#print axioms finiteStrictTransitiveRelation_cover_realized

end OPH.Provenance.FiniteCausetCompiler
