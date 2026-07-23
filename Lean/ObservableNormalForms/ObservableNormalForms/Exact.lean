import Mathlib

/-!
# Exact observable fibers

This module formalizes the exact set-theoretic and rewriting layer of
*Observation-Determined Normal Forms*.  Its definitions depend only on the mathematical
interfaces declared in the manuscript.

The main results correspond to the universal-property theorem, the
reachable-fiber proposition, the empty/singleton corollary, and the audited
terminal alternative.  The partial normalizer is genuinely proof carrying:
successful values inhabit the consistent subtype, while `none` is proved
equivalent to emptiness of the corresponding observable fiber.
-/

namespace ObservableNormalForms

open Relation

universe u v

section Fibers

variable {Q : Type u} {𝓑 : Type v}
variable {C : Set Q} {B : Q → 𝓑}

/-- The consistent states carrying observation `b`. -/
def fiber (C : Set Q) (B : Q → 𝓑) (b : 𝓑) : Set Q :=
  {q | q ∈ C ∧ B q = b}

/-- Observable determination means injectivity of the protected observation
map after restricting it to the consistent set. -/
def ObservableDetermined (C : Set Q) (B : Q → 𝓑) : Prop :=
  Set.InjOn B C

/-- The four clauses of the universal property of a canonical partial
normalizer.  `Option.none` is the paper's proved obstruction value `⊥`. -/
structure IsPartialNormalizer (C : Set Q) (B : Q → 𝓑)
    (N : Q → Option {q : Q // q ∈ C}) : Prop where
  sound : ∀ {x : Q} {c : {q : Q // q ∈ C}}, N x = some c → B c.1 = B x
  obstructed : ∀ x : Q, N x = none ↔ ¬ (fiber C B (B x)).Nonempty
  fixes : ∀ (c : Q) (hc : c ∈ C), N c = some ⟨c, hc⟩
  boundary_extensional : ∀ x y : Q, B x = B y → N x = N y

/-- Choice-based implementation of the canonical partial normalizer.  The
subsequent theorem proves that under observable determination this choice is
canonical and is the unique function satisfying `IsPartialNormalizer`. -/
noncomputable def canonicalPartialNormalizer (C : Set Q) (B : Q → 𝓑) :
    Q → Option {q : Q // q ∈ C} := by
  classical
  exact fun x =>
    if h : (fiber C B (B x)).Nonempty then
      some ⟨Classical.choose h, (Classical.choose_spec h).1⟩
    else
      none

theorem canonicalPartialNormalizer_spec
    (hinj : ObservableDetermined C B) :
    IsPartialNormalizer C B (canonicalPartialNormalizer C B) := by
  classical
  constructor
  · intro x c hN
    unfold canonicalPartialNormalizer at hN
    split at hN
    next h =>
      have hc : c = ⟨Classical.choose h, (Classical.choose_spec h).1⟩ :=
        Option.some.inj hN.symm
      subst c
      exact (Classical.choose_spec h).2
    next h => simp at hN
  · intro x
    simp [canonicalPartialNormalizer]
  · intro c hc
    have hfiber : (fiber C B (B c)).Nonempty := ⟨c, hc, rfl⟩
    rw [canonicalPartialNormalizer, dif_pos hfiber]
    apply congrArg some
    apply Subtype.ext
    exact hinj (Classical.choose_spec hfiber).1 hc (Classical.choose_spec hfiber).2
  · intro x y hxy
    by_cases hx : (fiber C B (B x)).Nonempty
    · have hy : (fiber C B (B y)).Nonempty := by
        rcases hx with ⟨c, hc, hcx⟩
        exact ⟨c, hc, hcx.trans hxy⟩
      rw [canonicalPartialNormalizer, dif_pos hx,
        canonicalPartialNormalizer, dif_pos hy]
      apply congrArg some
      apply Subtype.ext
      apply hinj (Classical.choose_spec hx).1 (Classical.choose_spec hy).1
      exact (Classical.choose_spec hx).2.trans
        ((Classical.choose_spec hy).2.trans hxy.symm).symm
    · have hy : ¬ (fiber C B (B y)).Nonempty := by
        intro hy
        rcases hy with ⟨c, hc, hcy⟩
        exact hx ⟨c, hc, hcy.trans hxy.symm⟩
      simp [canonicalPartialNormalizer, hx, hy]

theorem partialNormalizer_unique
    (hinj : ObservableDetermined C B)
    {N : Q → Option {q : Q // q ∈ C}}
    (hN : IsPartialNormalizer C B N) :
    N = canonicalPartialNormalizer C B := by
  classical
  funext x
  by_cases hx : (fiber C B (B x)).Nonempty
  · rcases hx with ⟨c, hc, hBc⟩
    have hNne : N x ≠ none := by
      intro hnone
      exact ((hN.obstructed x).mp hnone) ⟨c, hc, hBc⟩
    cases hNx : N x with
    | none => exact absurd hNx hNne
    | some d =>
        let hxfiber : (fiber C B (B x)).Nonempty := ⟨c, hc, hBc⟩
        rw [canonicalPartialNormalizer, dif_pos hxfiber]
        apply congrArg some
        apply Subtype.ext
        apply hinj d.2 (Classical.choose_spec hxfiber).1
        exact (hN.sound hNx).trans (Classical.choose_spec hxfiber).2.symm
  · have hnone : N x = none := (hN.obstructed x).mpr hx
    simp [canonicalPartialNormalizer, hx, hnone]

/-- Universal-property theorem, clauses (a) ↔ (b): observable determination is
equivalent to existence of a unique proof-carrying partial normalizer. -/
theorem observableDetermined_iff_unique_partialNormalizer :
    ObservableDetermined C B ↔
      ∃! N : Q → Option {q : Q // q ∈ C}, IsPartialNormalizer C B N := by
  constructor
  · intro hinj
    refine ⟨canonicalPartialNormalizer C B,
      canonicalPartialNormalizer_spec hinj, ?_⟩
    intro N hN
    exact partialNormalizer_unique hinj hN
  · rintro ⟨N, hN, _⟩
    intro c hc d hd hB
    have hEq := hN.boundary_extensional c d hB
    rw [hN.fixes c hc, hN.fixes d hd] at hEq
    exact congrArg Subtype.val (Option.some.inj hEq)

end Fibers

section Rewriting

variable {Q : Type u} {𝓑 : Type v}

/-- A state is a normal form when it has no outgoing rewrite step. -/
def IsNormalForm (r : Q → Q → Prop) (x : Q) : Prop :=
  ∀ y : Q, ¬ r x y

/-- Observation preservation for a one-step rewrite relation. -/
def ObservationPreserving (r : Q → Q → Prop) (B : Q → 𝓑) : Prop :=
  ∀ ⦃x y : Q⦄, r x y → B x = B y

/-- Completeness identifies the normal forms exactly with the declared
consistent set. -/
def CompleteFor (r : Q → Q → Prop) (C : Set Q) : Prop :=
  ∀ x : Q, IsNormalForm r x ↔ x ∈ C

/-- Weak normalization of one state. -/
def WeaklyNormalizing (r : Q → Q → Prop) (x : Q) : Prop :=
  ∃ nf : Q, ReflTransGen r x nf ∧ IsNormalForm r nf

/-- Confluence restricted to states with a fixed observation. -/
def ConfluentOnFiber (r : Q → Q → Prop) (B : Q → 𝓑) (b : 𝓑) : Prop :=
  ∀ ⦃x y z : Q⦄, B x = b → ReflTransGen r x y → ReflTransGen r x z →
    ∃ w : Q, ReflTransGen r y w ∧ ReflTransGen r z w

theorem observation_eq_of_reflTransGen
    {r : Q → Q → Prop} {B : Q → 𝓑}
    (hobs : ObservationPreserving r B)
    {x y : Q} (hxy : ReflTransGen r x y) :
    B x = B y := by
  induction hxy with
  | refl => rfl
  | tail _ hyz ih => exact ih.trans (hobs hyz)

/-- Reachable-fiber proposition: a reachable normal form belongs to the
consistent observable fiber of its source. -/
theorem reachable_normalForm_mem_fiber
    {r : Q → Q → Prop} {C : Set Q} {B : Q → 𝓑}
    (hobs : ObservationPreserving r B)
    (hcomplete : CompleteFor r C)
    {x nf : Q} (hred : ReflTransGen r x nf)
    (hnf : IsNormalForm r nf) :
    nf ∈ fiber C B (B x) := by
  exact ⟨(hcomplete nf).mp hnf,
    (observation_eq_of_reflTransGen hobs hred).symm⟩

/-- Empty/singleton corollary, clause (i): an empty consistent fiber cannot be hidden
behind a nonterminating search; no state in that observation class reaches a
normal form at all. -/
theorem empty_fiber_no_reachable_normalForm
    {r : Q → Q → Prop} {C : Set Q} {B : Q → 𝓑} {b : 𝓑}
    (hobs : ObservationPreserving r B)
    (hcomplete : CompleteFor r C)
    (hempty : ¬ (fiber C B b).Nonempty)
    {x : Q} (hx : B x = b) :
    ¬ WeaklyNormalizing r x := by
  rintro ⟨nf, hred, hnf⟩
  apply hempty
  have hmem := reachable_normalForm_mem_fiber hobs hcomplete hred hnf
  exact ⟨nf, hmem.1, hmem.2.trans hx⟩

/-- Empty/singleton corollary, clause (ii): a singleton consistent fiber forces every
reachable normal form in that observation class to be its unique point. -/
theorem singleton_fiber_forces_normalForm
    {r : Q → Q → Prop} {C : Set Q} {B : Q → 𝓑} {b : 𝓑} {c : Q}
    (hobs : ObservationPreserving r B)
    (hcomplete : CompleteFor r C)
    (hsingle : fiber C B b = {c})
    {x nf : Q} (hx : B x = b)
    (hred : ReflTransGen r x nf)
    (hnf : IsNormalForm r nf) :
    nf = c := by
  have hmem := reachable_normalForm_mem_fiber hobs hcomplete hred hnf
  have : nf ∈ fiber C B b := ⟨hmem.1, hmem.2.trans hx⟩
  rw [hsingle] at this
  simpa using this

/-- Empty/singleton corollary, clause (iii): weak normalization plus a singleton
consistent fiber implies confluence on that fiber.  This is the corrected
fiberwise version; weak normalization is not assumed on an empty fiber. -/
theorem singleton_fiber_weakNormalization_confluent
    {r : Q → Q → Prop} {C : Set Q} {B : Q → 𝓑} {b : 𝓑} {c : Q}
    (hobs : ObservationPreserving r B)
    (hcomplete : CompleteFor r C)
    (hsingle : fiber C B b = {c})
    (hwn : ∀ x : Q, B x = b → WeaklyNormalizing r x) :
    ConfluentOnFiber r B b := by
  intro x y z hx hxy hxz
  have hy : B y = b :=
    (observation_eq_of_reflTransGen hobs hxy).symm.trans hx
  have hz : B z = b :=
    (observation_eq_of_reflTransGen hobs hxz).symm.trans hx
  rcases hwn y hy with ⟨nfy, hynf, hnfy⟩
  rcases hwn z hz with ⟨nfz, hznf, hnfz⟩
  have hnfy_eq : nfy = c :=
    singleton_fiber_forces_normalForm hobs hcomplete hsingle hy hynf hnfy
  have hnfz_eq : nfz = c :=
    singleton_fiber_forces_normalForm hobs hcomplete hsingle hz hznf hnfz
  subst nfy
  subst nfz
  exact ⟨c, hynf, hznf⟩

/-- Certificate returned by the audited-terminal theorem below.  It records
not only quiescence, but also preservation of the input observation and the
precise condition under which the terminal state is semantically consistent. -/
structure AuditedTerminalCertificate
    (r : Q → Q → Prop) (C : Set Q) (B : Q → 𝓑) (x t : Q) : Prop where
  reaches : ReflTransGen r x t
  terminal : IsNormalForm r t
  observation : B t = B x
  consistent_iff_realizable :
    t ∈ C ↔ (fiber C B (B x)).Nonempty

/-- Audited-terminal alternative.  If rewriting preserves observations, all
states weakly normalize, semantic consistency implies quiescence, and the
observation map separates all terminal states, then every state reaches a
unique terminal state.  That terminal is consistent exactly when its input
observation has a consistent realization.

This theorem intentionally does not identify `C` with all terminal states:
unrealizable observation fibers may settle to an audited but semantically
inconsistent terminal state. -/
theorem existsUnique_auditedTerminal
    {r : Q → Q → Prop} {C : Set Q} {B : Q → 𝓑}
    (hobs : ObservationPreserving r B)
    (hwn : ∀ x : Q, WeaklyNormalizing r x)
    (hconsistentTerminal : C ⊆ {t : Q | IsNormalForm r t})
    (hterminalInjective : Set.InjOn B {t : Q | IsNormalForm r t})
    (x : Q) :
    ∃! t : Q, AuditedTerminalCertificate r C B x t := by
  rcases hwn x with ⟨t, hxt, htNormal⟩
  have htObs : B t = B x :=
    (observation_eq_of_reflTransGen hobs hxt).symm
  have htConsistent : t ∈ C ↔ (fiber C B (B x)).Nonempty := by
    constructor
    · intro htC
      exact ⟨t, htC, htObs⟩
    · rintro ⟨c, hcC, hcObs⟩
      have htc : t = c := by
        apply hterminalInjective htNormal (hconsistentTerminal hcC)
        exact htObs.trans hcObs.symm
      simpa [htc] using hcC
  refine ⟨t, ⟨hxt, htNormal, htObs, htConsistent⟩, ?_⟩
  intro u hu
  exact (hterminalInjective htNormal hu.terminal
    (htObs.trans hu.observation.symm)).symm

/-! ### Observation-leak counterexample -/

namespace ObservationLeakCounterexample

/-- Three states: an input, its intended same-observation terminal, and a
different-observation terminal reached by the leaking rewrite. -/
inductive State where
  | start
  | intended
  | escaped
  deriving DecidableEq, Fintype

open State

def observe : State → Bool
  | start => false
  | intended => false
  | escaped => true

def consistent : Set State := {intended, escaped}

def step : State → State → Prop :=
  fun x y => x = start ∧ y = escaped

theorem step_complete_for_consistent : CompleteFor step consistent := by
  intro x
  cases x <;> simp [IsNormalForm, step, consistent]

theorem intended_fiber_singleton :
    fiber consistent observe false = {intended} := by
  ext x
  cases x <;> simp [fiber, consistent, observe]

theorem start_reaches_escaped : ReflTransGen step start escaped :=
  ReflTransGen.single ⟨rfl, rfl⟩

theorem escaped_normal : IsNormalForm step escaped := by
  simp [IsNormalForm, step]

theorem step_not_observationPreserving :
    ¬ ObservationPreserving step observe := by
  intro h
  have hfalse := h (show step start escaped from ⟨rfl, rfl⟩)
  simp [observe] at hfalse

/-- Explicit separation requested by the audit: even with exact
normal-form completeness and a singleton consistent fiber at `false`, a
rewrite from that observation class reaches the other terminal when
observation preservation is dropped. -/
theorem singleton_fiber_insufficient_without_observation_preservation :
    CompleteFor step consistent ∧
      fiber consistent observe false = {intended} ∧
      ReflTransGen step start escaped ∧
      IsNormalForm step escaped ∧
      escaped ≠ intended ∧
      ¬ ObservationPreserving step observe := by
  exact ⟨step_complete_for_consistent, intended_fiber_singleton,
    start_reaches_escaped, escaped_normal, by decide,
    step_not_observationPreserving⟩

end ObservationLeakCounterexample

end Rewriting

end ObservableNormalForms
