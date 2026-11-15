import Geometry.SourceTemporalConnected

/-!
# A maximal preservation guard with a universal completing word

The guard admits a proposed seam exactly when retained observations and
the proposed scalar state distinguish the source record domain. It
does not query a successful route or a realized payload. Every finite
connected support has one completing word which the guard accepts from
every protected checkpoint whose record determines the current receiver
scalar. Initial sampling and the retained updates maintain this hypothesis.
Other proposals, including rejected ones,
are part of an attempt history. A proposal probability law is separate.
-/

set_option autoImplicit false

namespace OPH.SourceTemporalGuard
noncomputable section
open OPH.SourceTemporalTomography OPH.SourceTemporalConnected OPH.SourceTemporalNative
open ObserverPatchHolography.ScalarSeamRepair

variable {P X : Type*} [DecidableEq P]

structure Checkpoint (X P : Type*) where
  state : X → P → ℝ
  record : X → List ℝ

def Protected (c : Checkpoint X P) : Prop :=
  ∀ x y, c.record x = c.record y → c.state x = c.state y → x = y

def HasRoot (root : P) (c : Checkpoint X P) : Prop :=
  ∀ x y, c.record x = c.record y → c.state x root = c.state y root

def advance (root : P) (c : Checkpoint X P) (e : P × P) : Checkpoint X P where
  state x := pairAverage e.1 e.2 (c.state x)
  record x := pairAverage e.1 e.2 (c.state x) root :: c.record x

def run (root : P) : List (P × P) → Checkpoint X P → Checkpoint X P
  | [], c => c
  | e::es, c => run root es (advance root c e)

def Complete (root : P) (c : Checkpoint X P) (es : List (P × P)) : Prop :=
  Function.Injective (run root es c).record

def guard (root : P) (c : Checkpoint X P) (e : P × P) : Checkpoint X P := by
  classical
  exact if Protected (advance root c e) then advance root c e else c

def guardedRun (root : P) : List (P × P) → Checkpoint X P → Checkpoint X P
  | [], c => c
  | e::es, c => guardedRun root es (guard root c e)

theorem advance_has_root (root : P) (c : Checkpoint X P) (e : P × P) :
    HasRoot root (advance root c e) := by
  intro x y h
  exact (List.cons.inj h).1

/-- The extra receiver sample is determined by the proposed state.
Admission therefore needs only the old record and proposed native state. -/
theorem advance_protected_iff (root : P) (c : Checkpoint X P) (e : P × P) :
    Protected (advance root c e) ↔
      ∀ x y, c.record x = c.record y →
        pairAverage e.1 e.2 (c.state x) = pairAverage e.1 e.2 (c.state y) → x = y := by
  constructor
  · intro h x y hr hs
    exact h x y (congrArg₂ List.cons (congrFun hs root) hr) hs
  · intro h x y hr hs
    exact h x y (List.cons.inj hr).2 hs

def linearCheckpoint {V : Type*} [AddCommGroup V] [Module ℝ V] {m : ℕ}
    (L : Fin m → V →ₗ[ℝ] ℝ) (state : V →ₗ[ℝ] (P → ℝ)) : Checkpoint V P where
  state := state
  record x := List.ofFn (fun i => L i x)

/-- Exact matrix implementation of the abstract guard on a linear domain.
The test is computed from response maps, without realized source values. -/
theorem linear_protected_iff {V : Type*} [AddCommGroup V] [Module ℝ V] {m : ℕ}
    (L : Fin m → V →ₗ[ℝ] ℝ) (state : V →ₗ[ℝ] (P → ℝ)) :
    Protected (linearCheckpoint L state) ↔
      ∀ z, (∀ i, L i z = 0) → state z = 0 → z = 0 := by
  constructor
  · intro h z hz hs
    apply h z 0
    · apply congrArg List.ofFn
      funext i
      simp [hz i]
    · change state z = state 0
      simpa only [map_zero] using hs
  · intro h x y hr hs
    have he : (fun i => L i x) = (fun i => L i y) := by simpa [linearCheckpoint] using hr
    apply sub_eq_zero.mp
    apply h (x-y)
    · intro i
      simp [congrFun he i]
    · change state x = state y at hs
      rw [map_sub,hs,sub_self]

theorem run_record_extends (root : P) (es : List (P × P)) (c : Checkpoint X P)
    (x y : X) (h : (run root es c).record x = (run root es c).record y) :
    c.record x = c.record y := by
  induction es generalizing c with
  | nil => exact h
  | cons e es ih => exact (List.cons.inj (ih _ h)).2

theorem run_equal_data (root : P) (es : List (P × P)) (c : Checkpoint X P)
    (x y : X) (hr : c.record x = c.record y) (hs : c.state x = c.state y) :
    (run root es c).record x = (run root es c).record y := by
  induction es generalizing c with
  | nil => exact hr
  | cons e es ih =>
    apply ih
    · simp only [advance,hr,hs]
    · simp only [advance,hs]

/-- An erasing first move cannot be part of a complete observation word. -/
theorem complete_first_protected (root : P) (c : Checkpoint X P) (e : P × P)
    (es : List (P × P)) (h : Complete root c (e::es)) : Protected (advance root c e) := by
  intro x y hr hs
  exact h (run_equal_data root es (advance root c e) x y hr hs)

theorem guardedRun_complete_word (root : P) (c : Checkpoint X P) (es : List (P × P))
    (h : Complete root c es) : guardedRun root es c = run root es c := by
  classical
  induction es generalizing c with
  | nil => rfl
  | cons e es ih =>
    have hp := complete_first_protected root c e es h
    simpa only [guardedRun,guard,if_pos hp,run] using ih (advance root c e) h

def TraceAgrees (root : P) : List (P × P) → (P → ℝ) → (P → ℝ) → Prop
  | [], x, y => x root = y root
  | e::es, x, y => x root = y root ∧
      TraceAgrees root es (pairAverage e.1 e.2 x) (pairAverage e.1 e.2 y)

theorem trace_head (root : P) (es : List (P × P)) (x y : P → ℝ)
    (h : TraceAgrees root es x y) : x root = y root := by
  cases es with
  | nil => exact h
  | cons => exact h.1

theorem trace_append (root : P) (es fs : List (P × P)) (x y : P → ℝ)
    (h : TraceAgrees root (es ++ fs) x y) :
    TraceAgrees root fs (wordMap es x) (wordMap es y) := by
  induction es generalizing x y with
  | nil => exact h
  | cons e es ih => exact ih _ _ h.2

theorem trace_plan_samples {root : P} {K : Set P} (plan : Plan root K) (x y : P → ℝ)
    (h : TraceAgrees root (lower plan) x y) : (execute plan x).2 = (execute plan y).2 := by
  induction plan generalizing x y with
  | done => rfl
  | step K u vs hf hn hk he next ih =>
    have ht := trace_append root (route u vs) (lower next) x y h
    have hh := trace_head root (lower next) _ _ ht
    simp only [execute,List.cons.injEq]
    exact ⟨hh,ih _ _ ht⟩

theorem run_trace (root : P) (es : List (P × P)) (c : Checkpoint X P)
    (hc : HasRoot root c) (x y : X)
    (h : (run root es c).record x = (run root es c).record y) :
    TraceAgrees root es (c.state x) (c.state y) := by
  induction es generalizing c with
  | nil => exact hc x y h
  | cons e es ih =>
    exact ⟨hc x y (run_record_extends root (e::es) c x y h),
      ih (advance root c e) (advance_has_root root c e) h⟩

/-- The topology-derived word completes every protected checkpoint.
Sampling every mean retains the phase samples. -/
theorem plan_completes {root : P} (plan : Plan root {root}) (c : Checkpoint X P)
    (hp : Protected c) (hr : HasRoot root c) : Complete root c (lower plan) := by
  intro x y h
  have ht := run_trace root (lower plan) c hr x y h
  apply hp x y (run_record_extends root (lower plan) c x y h)
  apply completion_injective plan
  exact congrArg₂ List.cons (trace_head root (lower plan) _ _ ht)
    (trace_plan_samples plan _ _ ht)

theorem guard_preserves (root : P) (c : Checkpoint X P) (e : P × P)
    (hp : Protected c) (hr : HasRoot root c) :
    Protected (guard root c e) ∧ HasRoot root (guard root c e) := by
  classical
  by_cases h : Protected (advance root c e)
  · simp only [guard,if_pos h]
    exact ⟨h,advance_has_root root c e⟩
  · simpa only [guard,if_neg h] using And.intro hp hr

theorem guardedRun_preserves (root : P) (es : List (P × P)) (c : Checkpoint X P)
    (hp : Protected c) (hr : HasRoot root c) :
    Protected (guardedRun root es c) ∧ HasRoot root (guardedRun root es c) := by
  induction es generalizing c with
  | nil => exact ⟨hp,hr⟩
  | cons e es ih =>
    have h := guard_preserves root c e hp hr
    exact ih _ h.1 h.2

theorem guardedRun_append (root : P) (es fs : List (P × P)) (c : Checkpoint X P) :
    guardedRun root (es ++ fs) c = guardedRun root fs (guardedRun root es c) := by
  induction es generalizing c with
  | nil => rfl
  | cons e es ih => exact ih _

theorem guardedRun_record_extends (root : P) (es : List (P × P)) (c : Checkpoint X P)
    (x y : X) (h : (guardedRun root es c).record x = (guardedRun root es c).record y) :
    c.record x = c.record y := by
  classical
  induction es generalizing c with
  | nil => exact h
  | cons e es ih =>
    have he := ih (guard root c e) h
    by_cases hp : Protected (advance root c e)
    · simp only [guard,if_pos hp] at he
      exact (List.cons.inj he).2
    · simpa only [guard,if_neg hp] using he

/-- An occurrence of the completing word inside arbitrary proposal noise
forces complete receiver information. No failed prefix is discarded. -/
theorem completing_occurrence {root : P} (plan : Plan root {root})
    (before after : List (P × P)) (c : Checkpoint X P) (hp : Protected c) (hr : HasRoot root c) :
    Function.Injective (guardedRun root (before ++ lower plan ++ after) c).record := by
  have hi := guardedRun_preserves root before c hp hr
  have hc := plan_completes plan (guardedRun root before c) hi.1 hi.2
  rw [guardedRun_append,guardedRun_append,guardedRun_complete_word root _ _ hc]
  intro x y h
  exact hc (guardedRun_record_extends root after _ x y h)

/-- The universal completing word and its bound are derived from topology.
Probability or fairness assumptions can be applied to proposal occurrences
after this deterministic guarantee; they are not smuggled into the guard. -/
theorem connected_guard_completion [Fintype P] (G : SimpleGraph P)
    (hc : G.Connected) (root : P) :
    ∃ es : List (P × P), Supported G es ∧
      2 * es.length + Fintype.card P ≤ Fintype.card P ^ 2 ∧
      ∀ c : Checkpoint X P, Protected c → HasRoot root c →
        ∀ before after, Function.Injective
          (guardedRun root (before ++ es ++ after) c).record := by
  obtain ⟨plan,hs,_⟩ := connected_native_completion G hc root
  refine ⟨lower plan,hs,?_,?_⟩
  · simpa only [lower_length] using (complete_work_bound plan).2
  · intro c hp hr before after
    exact completing_occurrence plan before after c hp hr

end
end OPH.SourceTemporalGuard
