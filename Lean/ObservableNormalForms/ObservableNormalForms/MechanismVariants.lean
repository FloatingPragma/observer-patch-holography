import ObservableNormalForms.Exact
import ObservableNormalForms.Repair

/-!
# Conditional mechanism-variant comparison

This module provides a domain-neutral kernel for comparing finite mechanism
variants.  It separates executable traces, target and validity predicates,
protected behavior, authority quotients, and support-indexed repair relations.

The fixed-relation results below are thin wrappers around `StrongRepair`.
They do not turn a static section into a mechanism transition or establish
strategy preservation.  Pairwise behavior comparison is relative to an
explicit observation, protected set, and admissibility predicate.  A family
minimum additionally requires an explicitly declared comparison family.
-/

namespace ObservableNormalForms.MechanismVariants

open Relation

universe u v w x y z

/-- Variant-indexed dynamics with an explicit trace carrier.  Every trace is
certified to start in an initial state and to be reachable by the declared
step relation. -/
structure Dynamics (Variant : Type u) where
  State : Variant → Type v
  Step : (i : Variant) → State i → State i → Prop
  Initial : (i : Variant) → State i → Prop
  Trace : Variant → Type w
  start : (i : Variant) → Trace i → State i
  finish : (i : Variant) → Trace i → State i
  startsInitial : ∀ (i : Variant) (t : Trace i), Initial i (start i t)
  traceReachable : ∀ (i : Variant) (t : Trace i),
    ReflTransGen (Step i) (start i t) (finish i t)

abbrev State {Variant : Type u} (M : Dynamics Variant) := M.State

abbrev Step {Variant : Type u} (M : Dynamics Variant) := M.Step

abbrev Initial {Variant : Type u} (M : Dynamics Variant) := M.Initial

abbrev Trace {Variant : Type u} (M : Dynamics Variant) := M.Trace

/-- Input-relative reachability remains a mechanism-level relation. -/
def Reachable {Variant : Type u} (M : Dynamics Variant)
    (i : Variant) (x y : M.State i) : Prop :=
  ReflTransGen (M.Step i) x y

/-- Every data-carrying trace projects to the ordinary reachability relation. -/
theorem trace_to_reflTransGen {Variant : Type u} (M : Dynamics Variant)
    (i : Variant) (t : M.Trace i) :
    Reachable M i (M.start i t) (M.finish i t) :=
  M.traceReachable i t

/-- Direct reuse of the exact rewriting lemma for a state observation that is
preserved by each variant's step relation.  This theorem does not apply to a
path-sensitive trace observation without a separate bridge. -/
theorem stateObservation_eq_of_trace
    {Variant : Type u} {O : Type x} (M : Dynamics Variant)
    (observeState : (i : Variant) → M.State i → O)
    (hpreserves : ∀ i, ObservableNormalForms.ObservationPreserving
      (M.Step i) (observeState i))
    (i : Variant) (t : M.Trace i) :
    observeState i (M.start i t) = observeState i (M.finish i t) := by
  exact ObservableNormalForms.observation_eq_of_reflTransGen
    (hpreserves i) (M.traceReachable i t)

abbrev TargetBad {Variant : Type u} (M : Dynamics Variant) :=
  (i : Variant) → M.Trace i → Prop

abbrev ValidTrace {Variant : Type u} (M : Dynamics Variant) :=
  (i : Variant) → M.Trace i → Prop

abbrev CoreOK {Variant : Type u} (M : Dynamics Variant) :=
  (i : Variant) → M.Trace i → Prop

/-- A trace is admissible exactly when it satisfies the declared validity and
core predicates and avoids the selected target.  Initiality and reachability
are carried by `Dynamics`. -/
def Admissible {Variant : Type u} (M : Dynamics Variant)
    (target : TargetBad M) (valid : ValidTrace M) (core : CoreOK M)
    (i : Variant) (t : M.Trace i) : Prop :=
  valid i t ∧ core i t ∧ ¬ target i t

abbrev Observation {Variant : Type u} (M : Dynamics Variant) (O : Type x) :=
  (i : Variant) → M.Trace i → O

/-- Behaviors realized by traces satisfying one explicit admissibility
profile. -/
def RealizedBehavior {Variant : Type u} {O : Type x}
    (M : Dynamics Variant)
    (admissible : (i : Variant) → M.Trace i → Prop)
    (observe : Observation M O) (i : Variant) : Set O :=
  {o | ∃ t : M.Trace i, admissible i t ∧ observe i t = o}

/-- Protected behaviors not realized by a variant under one fixed profile. -/
def BehaviorCut {Variant : Type u} {O : Type x}
    (M : Dynamics Variant) (protectedSet : Set O)
    (admissible : (i : Variant) → M.Trace i → Prop)
    (observe : Observation M O) (i : Variant) : Set O :=
  protectedSet \ RealizedBehavior M admissible observe i

/-- Observation preservation means that the selected variant loses none of
the ex ante protected behaviors. -/
def ObservationPreserving {Variant : Type u} {O : Type x}
    (M : Dynamics Variant) (protectedSet : Set O)
    (admissible : (i : Variant) → M.Trace i → Prop)
    (observe : Observation M O) (i : Variant) : Prop :=
  BehaviorCut M protectedSet admissible observe i = ∅

/-- Descriptive preorder induced by inclusion of behavior cuts. -/
def BehaviorLE {Variant : Type u} {O : Type x}
    (cut : Variant → Set O) (i j : Variant) : Prop :=
  cut i ⊆ cut j

/-- Strict pairwise behavior order. -/
def BehaviorLT {Variant : Type u} {O : Type x}
    (cut : Variant → Set O) (i j : Variant) : Prop :=
  BehaviorLE cut i j ∧ ¬ BehaviorLE cut j i

theorem behaviorLT_of_subset_of_witness
    {Variant : Type u} {O : Type x} {cut : Variant → Set O}
    {i j : Variant} (hsub : BehaviorLE cut i j)
    {o : O} (hi : o ∉ cut i) (hj : o ∈ cut j) :
    BehaviorLT cut i j := by
  refine ⟨hsub, ?_⟩
  intro hreverse
  exact hi (hreverse hj)

/-- A family minimum is meaningful only relative to the supplied family and
requires comparison with every declared member. -/
def IsBehaviorMinimum {Variant : Type u} {O : Type x}
    (family : Set Variant) (cut : Variant → Set O) (i : Variant) : Prop :=
  i ∈ family ∧ ∀ j ∈ family, BehaviorLE cut i j

abbrev RawCaps (Variant : Type u) (Capability : Type v) :=
  Variant → Set Capability

abbrev AuthorityView (Capability : Type v) (Authority : Type w) :=
  Capability → Authority

def ViewedAuthority {Variant : Type u} {Capability : Type v}
    {Authority : Type w} (raw : RawCaps Variant Capability)
    (view : AuthorityView Capability Authority) (i : Variant) : Set Authority :=
  view '' raw i

/-- No authority class visible through `view` is added by `candidate` relative
to `baseline`.  Identity view recovers the raw-capability reading. -/
def ZeroNewAuthority {Variant : Type u} {Capability : Type v}
    {Authority : Type w} (raw : RawCaps Variant Capability)
    (view : AuthorityView Capability Authority)
    (baseline candidate : Variant) : Prop :=
  ViewedAuthority raw view candidate ⊆ ViewedAuthority raw view baseline

abbrev StrategicOK (Variant : Type u) := Variant → Prop

abbrev Supports (Variant : Type u) (Support : Type v) :=
  Variant → Set Support

/-- A relation may depend on both the mechanism variant and the chosen write
support, with support-dependent write and collar carriers. -/
abbrev Relation {Variant : Type u} {Support : Type v}
    (W : Variant → Support → Type w)
    (D : Variant → Support → Type x) :=
  (i : Variant) → (s : Support) → Set (W i s × D i s)

/-- Evidence that a support-dependent relation is exactly the range of an
encoding of admissible mechanism traces, rather than an unrelated free set. -/
structure RelationEncoding {Variant : Type u} (M : Dynamics Variant)
    (admissible : (i : Variant) → M.Trace i → Prop)
    {Support : Type v} (W : Variant → Support → Type w)
    (D : Variant → Support → Type x) (R : Relation W D)
    (i : Variant) (s : Support) where
  encode : {t : M.Trace i // admissible i t} → W i s × D i s
  exactRange : R i s = Set.range encode

abbrev SupportCost (Variant : Type u) (Support : Type v) (Cost : Type w) :=
  Variant → Support → Cost

/-- Typed provenance for one mechanism comparison.  Target, observation,
authority quotient, selected support, support-cost source, and comparison
family travel together.  The cost field records provenance only: neither
`cut` nor `PairwiseLT` depends on it. -/
structure ComparisonPolicy {Variant : Type u} (M : Dynamics Variant)
    (O : Type x) (Capability : Type y) (Authority : Type z)
    (Support : Type v) (Cost : Type w) where
  target : TargetBad M
  valid : ValidTrace M
  core : CoreOK M
  observe : Observation M O
  protectedSet : Set O
  rawCaps : RawCaps Variant Capability
  authorityView : AuthorityView Capability Authority
  authorityBaseline : Variant
  supports : Supports Variant Support
  selectedSupport : Variant → Support
  supportCost : SupportCost Variant Support Cost
  family : Set Variant

namespace ComparisonPolicy

def admissible {Variant : Type u} {M : Dynamics Variant}
    {O : Type x} {Capability : Type y} {Authority : Type z}
    {Support : Type v} {Cost : Type w}
    (p : ComparisonPolicy M O Capability Authority Support Cost)
    (i : Variant) (t : M.Trace i) : Prop :=
  Admissible M p.target p.valid p.core i t

def cut {Variant : Type u} {M : Dynamics Variant}
    {O : Type x} {Capability : Type y} {Authority : Type z}
    {Support : Type v} {Cost : Type w}
    (p : ComparisonPolicy M O Capability Authority Support Cost)
    (i : Variant) : Set O :=
  BehaviorCut M p.protectedSet p.admissible p.observe i

def Eligible {Variant : Type u} {M : Dynamics Variant}
    {O : Type x} {Capability : Type y} {Authority : Type z}
    {Support : Type v} {Cost : Type w}
    (p : ComparisonPolicy M O Capability Authority Support Cost)
    (i : Variant) : Prop :=
  i ∈ p.family ∧
    p.selectedSupport i ∈ p.supports i ∧
    ZeroNewAuthority p.rawCaps p.authorityView p.authorityBaseline i

/-- Selected support cost is exposed for audit provenance, but does not enter
the behavior order or repair-width statements. -/
def selectedCost {Variant : Type u} {M : Dynamics Variant}
    {O : Type x} {Capability : Type y} {Authority : Type z}
    {Support : Type v} {Cost : Type w}
    (p : ComparisonPolicy M O Capability Authority Support Cost)
    (i : Variant) : Cost :=
  p.supportCost i (p.selectedSupport i)

def PairwiseLT {Variant : Type u} {M : Dynamics Variant}
    {O : Type x} {Capability : Type y} {Authority : Type z}
    {Support : Type v} {Cost : Type w}
    (p : ComparisonPolicy M O Capability Authority Support Cost)
    (i j : Variant) : Prop :=
  p.Eligible i ∧ p.Eligible j ∧ BehaviorLT p.cut i j

end ComparisonPolicy

/-- Coverage of a collar declared before inspecting a repair result. -/
def AdmissibleCollar {W : Type u} {D : Type v}
    (R : Set (W × D)) (collar : Set D) : Prop :=
  ∀ d, d ∈ collar → ∃ w, (w, d) ∈ R

/-- Coverage of the full declared collar is exactly surjectivity of the
restricted relation's collar projection.  This bridge only repackages
witnesses and introduces no choice. -/
theorem collarProjectionSurjective_iff_full_admissibleCollar
    {W : Type u} {D : Type v} (R : Set (W × D)) :
    CollarProjectionSurjective R ↔ AdmissibleCollar R Set.univ := by
  constructor
  · intro hsurj d _
    rcases hsurj d with ⟨z, hz⟩
    refine ⟨z.1.1, ?_⟩
    rw [← hz]
    simpa only [Prod.eta] using z.2
  · intro hcoverage d
    rcases hcoverage d (Set.mem_univ d) with ⟨w, hw⟩
    exact ⟨⟨(w, d), hw⟩, rfl⟩

/-- Surjectivity of the collar component of the declared trace encoding. -/
def EncodedCollarSurjective
    {Variant : Type u} {M : Dynamics Variant}
    {Support : Type v} {W : Variant → Support → Type w}
    {D : Variant → Support → Type x} {R : Relation W D}
    {admissible : (i : Variant) → M.Trace i → Prop}
    {i : Variant} {s : Support}
    (E : RelationEncoding M admissible W D R i s) : Prop :=
  Function.Surjective (fun t => (E.encode t).2)

/-- `exactRange` turns full fixed-relation collar coverage into surjectivity of
the actual admissible-trace encoding.  The proof is witness-preserving and
uses no choice. -/
theorem admissibleCollar_iff_encodedCollarSurjective
    {Variant : Type u} {M : Dynamics Variant}
    {Support : Type v} {W : Variant → Support → Type w}
    {D : Variant → Support → Type x} {R : Relation W D}
    {admissible : (i : Variant) → M.Trace i → Prop}
    {i : Variant} {s : Support}
    (E : RelationEncoding M admissible W D R i s) :
    AdmissibleCollar (R i s) Set.univ ↔ EncodedCollarSurjective E := by
  constructor
  · intro hcoverage d
    rcases hcoverage d (Set.mem_univ d) with ⟨w, hw⟩
    rw [E.exactRange] at hw
    rcases hw with ⟨t, ht⟩
    exact ⟨t, congrArg Prod.snd ht⟩
  · intro hsurj d _
    rcases hsurj d with ⟨t, ht⟩
    refine ⟨(E.encode t).1, ?_⟩
    rw [E.exactRange]
    exact ⟨t, Prod.ext rfl ht⟩

/-- Direct fixed-relation reuse of the OPH collar-section criterion.  The
result is static coverage only. -/
theorem strongRepair_exists_iff_full_admissibleCollar
    {W : Type u} {D : Type v} [Nonempty W] (R : Set (W × D)) :
    Nonempty (StrongRepair R) ↔ AdmissibleCollar R Set.univ := by
  exact (strongRepair_exists_iff_collarProjection_surjective R).trans
    (collarProjectionSurjective_iff_full_admissibleCollar R)

/-- Actual proof composition from a declared exact trace range, through
encoded collar surjectivity, to existence of a fixed-relation `StrongRepair`.
This does not claim same-initial repair, all-path completeness, or full
`Repair` width. -/
theorem strongRepair_exists_iff_encodedCollarSurjective
    {Variant : Type u} {M : Dynamics Variant}
    {Support : Type v} {W : Variant → Support → Type w}
    {D : Variant → Support → Type x} {R : Relation W D}
    {admissible : (i : Variant) → M.Trace i → Prop}
    {i : Variant} {s : Support} [Nonempty (W i s)]
    (E : RelationEncoding M admissible W D R i s) :
    Nonempty (StrongRepair (R i s)) ↔ EncodedCollarSurjective E := by
  exact (strongRepair_exists_iff_full_admissibleCollar (R i s)).trans
    (admissibleCollar_iff_encodedCollarSurjective E)

/-- A concrete hole in a fixed relation rules out total exact
collar-preserving repair for that relation. -/
theorem no_strongRepair_of_missing_admissibleCollar
    {W : Type u} {D : Type v} [Nonempty W]
    {R : Set (W × D)} {d : D} (hmissing : ¬ ∃ w, (w, d) ∈ R) :
    IsEmpty (StrongRepair R) := by
  apply no_strongRepair_of_missing_collar (d := d)
  rintro ⟨z, hz⟩
  apply hmissing
  refine ⟨z.1.1, ?_⟩
  rw [← hz]
  simpa only [Prod.eta] using z.2

end ObservableNormalForms.MechanismVariants
