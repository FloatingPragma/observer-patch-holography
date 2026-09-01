import ObserverPatchHolography.Provenance.CausalInterval

/-!
# Consequences of a certified refinement order map

A refinement receipt between two semantic event logs supplies a map on
event identifiers and the substantive edge-preservation clause: every fine
direct-parent edge lands in the coarse reflexive generated precedence.  The
module does not derive that clause from register provenance or prove that an
arbitrary map is a semantic refinement.  Conditional on the supplied
certificate, closure preservation and interval inclusion follow.

The committed witness coarsens the four-event Boolean diamond onto a
three-event chain by merging the two incomparable responses into one
coarse response.  The reversed assignment on the same logs violates the
edge clause, so it inhabits no refinement receipt: the clause is
load-bearing, and a refinement cannot invert generated precedence.

Interval equality, surjectivity onto coarse intervals, and any continuum
refinement limit are not supplied.
-/

namespace OPH.Provenance

universe u₁ v₁ w₁ u₂ v₂ w₂

variable {RegisterF : Type u₁} {ValueF : Type v₁} {EventIdF : Type w₁}
variable {RegisterC : Type u₂} {ValueC : Type v₂} {EventIdC : Type w₂}
variable [DecidableEq RegisterF] [DecidableEq RegisterC]

/-- A certified candidate refinement map between two semantic event logs.
The `edge_natural` field is a load-bearing supplied certificate, not a
derived theorem that the map is a physical or provenance refinement. -/
structure LogRefinement
    (Lf : SemanticEventLog RegisterF ValueF EventIdF)
    (Lc : SemanticEventLog RegisterC ValueC EventIdC) where
  /-- The coarse image of each fine event identifier. -/
  map : EventIdF → EventIdC
  /-- Every fine direct-parent edge lands in the coarse reflexive
  precedence; merged endpoints are declared coincident. -/
  edge_natural : ∀ {e f : EventIdF}, Lf.ParentEdge e f →
    Lc.GeneratedBeforeEq (map e) (map f)

namespace LogRefinement

variable {Lf : SemanticEventLog RegisterF ValueF EventIdF}
variable {Lc : SemanticEventLog RegisterC ValueC EventIdC}

/-- Conditional closure preservation: the fine generated precedence maps
into the coarse reflexive precedence under the supplied edge certificate. -/
theorem generatedBefore_maps_under_edge_certificate
    (R : LogRefinement Lf Lc)
    {e f : EventIdF} (h : Lf.GeneratedBefore e f) :
    Lc.GeneratedBeforeEq (R.map e) (R.map f) := by
  induction h with
  | single hedge => exact R.edge_natural hedge
  | tail _ hedge ih =>
      exact Lc.generatedBeforeEq_trans ih (R.edge_natural hedge)

/-- Conditional preservation of reflexive generated precedence. -/
theorem generatedBeforeEq_maps_under_edge_certificate
    (R : LogRefinement Lf Lc)
    {e f : EventIdF} (h : Lf.GeneratedBeforeEq e f) :
    Lc.GeneratedBeforeEq (R.map e) (R.map f) := by
  rcases h with rfl | h
  · exact Or.inl rfl
  · exact R.generatedBefore_maps_under_edge_certificate h

/-- Conditional interval inclusion: the image of a fine causal interval
lies in the coarse interval of the image endpoints.  Surjectivity and
interval equality are not consequences. -/
theorem interval_maps_under_edge_certificate
    (R : LogRefinement Lf Lc) {a b x : EventIdF}
    (hx : x ∈ Lf.interval a b) :
    R.map x ∈ Lc.interval (R.map a) (R.map b) :=
  ⟨R.generatedBeforeEq_maps_under_edge_certificate hx.1,
    R.generatedBeforeEq_maps_under_edge_certificate hx.2⟩

end LogRefinement

/-! ## The diamond-to-chain witness -/

namespace DiamondChainRefinement

open BooleanDiamond

/-- The three-event coarse chain: injection, merged response, answer,
each citing its predecessor's register. -/
def chainState0 : VersionedState (Fin 3) Bool (Fin 3) :=
  ⟨![false, false, false], ![none, none, none]⟩

def chainState1 : VersionedState (Fin 3) Bool (Fin 3) :=
  ⟨![true, false, false], ![some 0, none, none]⟩

def chainState2 : VersionedState (Fin 3) Bool (Fin 3) :=
  ⟨![true, true, false], ![some 0, some 1, none]⟩

def chainState3 : VersionedState (Fin 3) Bool (Fin 3) :=
  ⟨![true, true, true], ![some 0, some 1, some 2]⟩

/-- The coarse injection writes register 0. -/
def coarseInjection : SemanticCommit (Fin 3) Bool (Fin 3) where
  eventId := 0
  before := chainState0
  after := chainState1
  readSet := ∅
  writeSet := {0}
  causalSupp := ∅
  supp_subset_read := by decide
  frame_value := by decide
  frame_writer := by decide
  stamp := by decide

/-- The merged coarse response cites register 0 and writes register 1. -/
def coarseResponse : SemanticCommit (Fin 3) Bool (Fin 3) where
  eventId := 1
  before := chainState1
  after := chainState2
  readSet := {0}
  writeSet := {1}
  causalSupp := {0}
  supp_subset_read := by decide
  frame_value := by decide
  frame_writer := by decide
  stamp := by decide

/-- The coarse answer cites register 1 and writes register 2. -/
def coarseAnswer : SemanticCommit (Fin 3) Bool (Fin 3) where
  eventId := 2
  before := chainState2
  after := chainState3
  readSet := {1}
  writeSet := {2}
  causalSupp := {1}
  supp_subset_read := by decide
  frame_value := by decide
  frame_writer := by decide
  stamp := by decide

/-- The three-event chain log. -/
def chainLog3 : SemanticEventLog (Fin 3) Bool (Fin 3) where
  commitOf := ![coarseInjection, coarseResponse, coarseAnswer]
  commitOf_eventId := by decide
  rank := ![0, 1, 2]
  rank_lt_of_parent := by
    have h : ∀ c d : Fin 3,
        AuthenticatedDirectSemanticParent
            (![coarseInjection, coarseResponse, coarseAnswer] c)
            (![coarseInjection, coarseResponse, coarseAnswer] d) →
          (![0, 1, 2] : Fin 3 → ℕ) c < ![0, 1, 2] d := by decide
    intro c d
    exact h c d

theorem chain3_edge_01 : chainLog3.ParentEdge 0 1 := by decide

theorem chain3_edge_12 : chainLog3.ParentEdge 1 2 := by decide

/-- The coincidence-merging coarse map: both diamond responses land on
the coarse response. -/
def diamondToChain : Fin 4 → Fin 3 := ![0, 1, 1, 2]

/-- The diamond refines onto the chain: every diamond edge lands in the
chain's reflexive precedence, with the two responses merged. -/
def refinement : LogRefinement BooleanDiamond.log chainLog3 where
  map := diamondToChain
  edge_natural := by
    intro e f hedge
    have hcase := (BooleanDiamond.parentEdge_iff e f).mp hedge
    rcases hcase with ⟨he, hf⟩ | ⟨he, hf⟩ | ⟨he, hf⟩ | ⟨he, hf⟩ <;>
      subst he <;> subst hf
    · exact Or.inr (Relation.TransGen.single chain3_edge_01)
    · exact Or.inr (Relation.TransGen.single chain3_edge_01)
    · exact Or.inr (Relation.TransGen.single chain3_edge_12)
    · exact Or.inr (Relation.TransGen.single chain3_edge_12)

/-- The refined diamond interval lands in the full coarse chain
interval. -/
theorem diamond_interval_refines {x : Fin 4}
    (hx : x ∈ BooleanDiamond.log.interval 0 3) :
    diamondToChain x ∈ chainLog3.interval 0 2 :=
  refinement.interval_maps_under_edge_certificate hx

/-- The reversed assignment inhabits no refinement receipt: it would
carry the first diamond edge against the chain's precedence. -/
theorem reversed_map_not_natural :
    ¬ ∃ R : LogRefinement BooleanDiamond.log chainLog3,
      R.map = ![2, 1, 1, 0] := by
  rintro ⟨R, hmap⟩
  have hedge : BooleanDiamond.log.ParentEdge 0 1 := by decide
  have h := R.edge_natural hedge
  rw [hmap] at h
  rcases h with heq | hlt
  · exact absurd heq (by decide)
  · have := chainLog3.rank_lt_of_generatedBefore hlt
    exact absurd this (by decide)

end DiamondChainRefinement

#print axioms LogRefinement.generatedBefore_maps_under_edge_certificate
#print axioms LogRefinement.interval_maps_under_edge_certificate
#print axioms DiamondChainRefinement.refinement
#print axioms DiamondChainRefinement.diamond_interval_refines
#print axioms DiamondChainRefinement.reversed_map_not_natural

end OPH.Provenance
