import ObservableNormalForms.Exact

/-!
# Observation-relative endpoint uniqueness

This module separates ordinary same-source rewrite confluence from the
observer-facing question relevant to normal forms: do any two sources with
the same protected observation have equivalent normal endpoints?

The main theorem gives an exact characterization.  Under observation
preservation and completeness of the declared consistent set, endpoint
uniqueness modulo an arbitrary relation is equivalent to that observation
identifying consistent states modulo the same relation.  Weak normalization
is needed only for existence of endpoints, not for the equivalence itself.
-/

namespace ObservableNormalForms

open Relation

universe u v

section Characterization

variable {Q : Type u} {𝓑 : Type v}

/-- The protected observation identifies consistent states modulo `E`.

No algebraic laws are imposed on `E`: applications may instantiate it with
equality, an equivalence relation, or any coarser observational relation. -/
def BoundaryIdentifiesModulo
    (C : Set Q) (B : Q → 𝓑) (E : Q → Q → Prop) : Prop :=
  ∀ ⦃c d : Q⦄, c ∈ C → d ∈ C → B c = B d → E c d

/-- Any normal endpoints reached from possibly different sources with the
same protected observation are equivalent modulo `E`.

This is deliberately a cross-source property.  Ordinary confluence compares
two reductions from one source and does not, by itself, imply this property. -/
def ObserverEndpointUniqueModulo
    (r : Q → Q → Prop) (B : Q → 𝓑) (E : Q → Q → Prop) : Prop :=
  ∀ ⦃x y nx ny : Q⦄,
    B x = B y →
    ReflTransGen r x nx →
    ReflTransGen r y ny →
    IsNormalForm r nx →
    IsNormalForm r ny →
    E nx ny

/-- Observer-facing endpoint uniqueness is exactly injectivity of the
protected observation on consistent normal forms, modulo the chosen relation.

Observation preservation transports equality of the source observations to
the endpoints.  Completeness supplies both directions between consistency and
normality; its reverse direction proves the converse by taking zero-step
reductions from two consistent states. -/
theorem boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo
    {r : Q → Q → Prop} {C : Set Q} {B : Q → 𝓑} {E : Q → Q → Prop}
    (hobs : ObservationPreserving r B)
    (hcomplete : CompleteFor r C) :
    BoundaryIdentifiesModulo C B E ↔
      ObserverEndpointUniqueModulo r B E := by
  constructor
  · intro hident x y nx ny hxy hxnx hyny hnx hny
    apply hident ((hcomplete nx).mp hnx) ((hcomplete ny).mp hny)
    have hnxObs : B x = B nx :=
      observation_eq_of_reflTransGen hobs hxnx
    have hnyObs : B y = B ny :=
      observation_eq_of_reflTransGen hobs hyny
    exact hnxObs.symm.trans (hxy.trans hnyObs)
  · intro hunique c d hc hd hcd
    exact hunique hcd ReflTransGen.refl ReflTransGen.refl
      ((hcomplete c).mpr hc) ((hcomplete d).mpr hd)

/-- With weak normalization, the characterization also produces an
equivalent pair of normal endpoints for every pair of equally observed
sources.  The preceding theorem controls every other choice of endpoints. -/
theorem exists_equivalent_observer_endpoints
    {r : Q → Q → Prop} {C : Set Q} {B : Q → 𝓑} {E : Q → Q → Prop}
    (hobs : ObservationPreserving r B)
    (hcomplete : CompleteFor r C)
    (hident : BoundaryIdentifiesModulo C B E)
    (hwn : ∀ x : Q, WeaklyNormalizing r x)
    {x y : Q} (hxy : B x = B y) :
    ∃ nx ny : Q,
      ReflTransGen r x nx ∧ IsNormalForm r nx ∧
      ReflTransGen r y ny ∧ IsNormalForm r ny ∧ E nx ny := by
  rcases hwn x with ⟨nx, hxnx, hnx⟩
  rcases hwn y with ⟨ny, hyny, hny⟩
  refine ⟨nx, ny, hxnx, hnx, hyny, hny, ?_⟩
  exact (boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo
    hobs hcomplete).mp hident hxy hxnx hyny hnx hny

end Characterization

/-! ### A complete two-bit repair system -/

namespace TwoBitRepair

/-- The first bit is protected; the second bit is a repairable defect flag. -/
abbrev State := Bool × Bool

def observe (q : State) : Bool := q.1

def coarseObserve (_q : State) : Unit := ()

def consistent : Set State := {q | q.2 = false}

/-- A defective state repairs in one step without changing its protected bit. -/
def step (x y : State) : Prop :=
  x.2 = true ∧ y = (x.1, false)

theorem step_observationPreserving :
    ObservationPreserving step observe := by
  intro x y hxy
  rcases hxy with ⟨_, rfl⟩
  rfl

theorem step_completeFor_consistent :
    CompleteFor step consistent := by
  intro x
  rcases x with ⟨b, d⟩
  cases d <;> simp [IsNormalForm, step, consistent]

theorem boundary_identifies_consistent :
    BoundaryIdentifiesModulo consistent observe Eq := by
  intro c d hc hd hobs
  rcases c with ⟨bc, dc⟩
  rcases d with ⟨bd, dd⟩
  simp only [consistent, Set.mem_setOf_eq] at hc hd
  simp only [observe] at hobs
  subst dc
  subst dd
  subst bd
  rfl

theorem consistent_fiber_singleton (b : Bool) :
    fiber consistent observe b = {(b, false)} := by
  ext x
  rcases x with ⟨bx, dx⟩
  cases bx <;> cases dx <;> cases b <;>
    simp [fiber, consistent, observe]

theorem weaklyNormalizing (x : State) :
    WeaklyNormalizing step x := by
  rcases x with ⟨b, d⟩
  cases d
  · exact ⟨(b, false), ReflTransGen.refl,
      (step_completeFor_consistent (b, false)).mpr rfl⟩
  · exact ⟨(b, false), ReflTransGen.single ⟨rfl, rfl⟩,
      (step_completeFor_consistent (b, false)).mpr rfl⟩

/-- The fine observation produces genuinely unique endpoints for all equally
observed sources. -/
theorem observerEndpointUnique :
    ObserverEndpointUniqueModulo step observe Eq :=
  (boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo
    step_observationPreserving step_completeFor_consistent).mp
    boundary_identifies_consistent

theorem confluent_on_observation_fiber (b : Bool) :
    ConfluentOnFiber step observe b :=
  singleton_fiber_weakNormalization_confluent
    step_observationPreserving step_completeFor_consistent
    (consistent_fiber_singleton b) (fun x _ => weaklyNormalizing x)

/-- The repair relation is same-source confluent even when viewed through the
coarse observation below.  Combined with
`coarse_observerEndpointUnique_fails`, this separates ordinary confluence from
cross-source observer-facing endpoint uniqueness. -/
theorem coarse_confluent :
    ConfluentOnFiber step coarseObserve () := by
  intro x y z _ hxy hxz
  exact confluent_on_observation_fiber (observe x) rfl hxy hxz

/-- Weak normalization makes the endpoint guarantee nonvacuous for every pair
of sources with the same protected bit. -/
theorem exists_equal_endpoints
    {x y : State} (hxy : observe x = observe y) :
    ∃ nx ny : State,
      ReflTransGen step x nx ∧ IsNormalForm step nx ∧
      ReflTransGen step y ny ∧ IsNormalForm step ny ∧ nx = ny :=
  exists_equivalent_observer_endpoints
    step_observationPreserving step_completeFor_consistent
    boundary_identifies_consistent weaklyNormalizing hxy

/-- A boundary that discards the protected bit fails to identify even the
consistent states.  Thus completeness and normalization cannot compensate for
an information-deficient observation. -/
theorem coarse_boundary_does_not_identify :
    ¬ BoundaryIdentifiesModulo consistent coarseObserve Eq := by
  intro h
  have hbad := h (c := (false, false)) (d := (true, false))
    (by rfl) (by rfl) rfl
  exact Bool.false_ne_true (congrArg Prod.fst hbad)

theorem coarse_observerEndpointUnique_fails :
    ¬ ObserverEndpointUniqueModulo step coarseObserve Eq := by
  intro h
  have hbad := h (x := (false, false)) (y := (true, false))
    (nx := (false, false)) (ny := (true, false))
    rfl ReflTransGen.refl ReflTransGen.refl
    ((step_completeFor_consistent (false, false)).mpr rfl)
    ((step_completeFor_consistent (true, false)).mpr rfl)
  exact Bool.false_ne_true (congrArg Prod.fst hbad)

end TwoBitRepair

end ObservableNormalForms
