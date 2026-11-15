import Geometry.SourceTemporalNative
import Geometry.SourceSeamPathTomography

/-!
# Native calibrated paths and continuation without loss

A receiver can extend a calibrated connected region by sweeping one new
vertex through a simple path of calibrated vertices. The endpoint
sample recovers the new value and the evolved region. This is a construction
from scalar seam means, with no record payload in the schedule. Paths and
their supported incidence are explicit inputs; no metric menu is selected.
-/

set_option autoImplicit false

namespace OPH.SourceTemporalTomography
noncomputable section
open OPH.SourceTemporalNative OPH.SourceSeamPathTomography
open ObserverPatchHolography.ScalarSeamRepair

universe u
variable {P : Type u} [DecidableEq P]

def route (u : P) : List P → List (P × P)
  | [] => []
  | v::vs => (u,v)::route v vs

def destination (u : P) : List P → P
  | [] => u
  | v::vs => destination v vs

theorem route_length (u : P) (vs : List P) : (route u vs).length = vs.length := by
  induction vs generalizing u with
  | nil => rfl
  | cons v vs ih => simp [route,ih]

/-- The abstract sweep endpoint is the actual local value of the native word.
Simple paths prevent an earlier source from reappearing as an unseen relay. -/
theorem native_endpoint (u : P) (vs : List P) (h : (u::vs).Nodup) (x : P → ℝ) :
    wordMap (route u vs) x (destination u vs) = endpoint (x u) (vs.map x) := by
  induction vs generalizing u x with
  | nil => rfl
  | cons v vs ih =>
    have hu := (List.nodup_cons.mp h).1
    have hv := (List.nodup_cons.mp (List.nodup_cons.mp h).2).1
    have hm : vs.map (pairAverage u v x) = vs.map x := by
      apply List.map_congr_left
      intro p hp
      exact pairAverage_off_seam u v p x
        (fun he => hu (by simp [← he,hp])) (fun he => hv (he ▸ hp))
    simpa only [route,wordMap,LinearMap.comp_apply,destination,ih v
      (List.nodup_cons.mp h).2,hm,pairAverage_right,List.map_cons,endpoint]

/-- Only the endpoint sample and calibrated relay values enter
the inverse. The unknown remote source value is the conclusion. -/
theorem native_path_decode (u : P) (vs : List P) (h : (u::vs).Nodup) (x : P → ℝ) :
    invert (vs.map x) (wordMap (route u vs) x (destination u vs)) = x u := by
  rw [native_endpoint u vs h x,invert_endpoint]

theorem native_path_identifies (u : P) (vs : List P) (h : (u::vs).Nodup)
    (x y : P → ℝ) (hr : ∀ p ∈ vs, x p = y p)
    (hs : wordMap (route u vs) x (destination u vs) =
      wordMap (route u vs) y (destination u vs)) : x u = y u := by
  have hm : vs.map x = vs.map y := List.map_congr_left hr
  rw [← native_path_decode u vs h x,hs,hm,native_path_decode u vs h y]

/-- Sweeping never changes a vertex outside its declared path. -/
theorem route_off_path (u p : P) (vs : List P) (hp : p ∉ u::vs) (x : P → ℝ) :
    wordMap (route u vs) x p = x p := by
  induction vs generalizing u x with
  | nil => rfl
  | cons v vs ih =>
    simp only [List.mem_cons,not_or] at hp
    rw [route,wordMap,LinearMap.comp_apply,ih v (by simp [hp])]
    exact pairAverage_off_seam u v p x hp.1 hp.2.1

theorem route_agreement (u : P) (vs : List P) (x y : P → ℝ)
    (h : ∀ p ∈ u::vs, x p = y p) :
    ∀ p ∈ u::vs, wordMap (route u vs) x p = wordMap (route u vs) y p := by
  induction vs generalizing u x y with
  | nil => simpa [route,wordMap] using h
  | cons v vs ih =>
    have hm : ∀ p ∈ v::vs, pairAverage u v x p = pairAverage u v y p := by
      intro p hp
      by_cases he : p=u ∨ p=v
      · simp [pairAverage,he,h u (by simp),h v (by simp)]
      · simp [pairAverage,he,h p (by simp only [List.mem_cons] at hp ⊢; tauto)]
    intro p hp
    simp only [route,wordMap,LinearMap.comp_apply]
    by_cases hv : p ∈ v::vs
    · exact ih v _ _ hm p hv
    · have he : p=u := by simpa [hv] using hp
      subst p
      rw [route_off_path v u vs hv,route_off_path v u vs hv]
      simp [h u (by simp),h v (by simp)]

/-- A calibrated phase preserves full distinguishability when its previous
record is retained. Equality of the evolved whole state cannot hide a lost
initial distinction. -/
theorem phase_preserves (u : P) (vs : List P) (h : (u::vs).Nodup)
    (x y : P → ℝ) (hr : ∀ p ∈ vs, x p = y p)
    (hs : wordMap (route u vs) x = wordMap (route u vs) y) : x = y := by
  have hu := native_path_identifies u vs h x y hr (congrFun hs _)
  funext p
  by_cases hp : p ∈ u::vs
  · rcases List.mem_cons.mp hp with rfl | hp
    · exact hu
    · exact hr p hp
  · simpa only [route_off_path u p vs hp] using congrFun hs p

/-- After one endpoint sample the whole enlarged calibrated region has
identical current values on an observation fiber. -/
theorem phase_calibrates (K : Set P) (u : P) (vs : List P) (h : (u::vs).Nodup)
    (hk : ∀ p ∈ vs, p ∈ K) (x y : P → ℝ) (hr : ∀ p ∈ K, x p = y p)
    (hs : wordMap (route u vs) x (destination u vs) =
      wordMap (route u vs) y (destination u vs)) :
    ∀ p ∈ insert u K, wordMap (route u vs) x p = wordMap (route u vs) y p := by
  have hu := native_path_identifies u vs h x y (fun p hp => hr p (hk p hp)) hs
  have ha : ∀ p ∈ u::vs, x p = y p := by
    intro p hp
    rcases List.mem_cons.mp hp with rfl | hp
    · exact hu
    · exact hr p (hk p hp)
  intro p hp
  by_cases hm : p ∈ u::vs
  · exact route_agreement u vs x y ha p hm
  · rw [route_off_path u p vs hm,route_off_path u p vs hm]
    exact hr p (by simpa [show p ≠ u from fun he => hm (by simp [he])] using hp)

/-- A finite topology certificate. Each new vertex has a simple path to the
same receiver through calibrated vertices. Termination requires coverage of
every port. No scalar values or decoder correctness occur in this type. -/
inductive Plan (root : P) : Set P → Type u
  | done (K : Set P) (covered : ∀ p, p ∈ K) : Plan root K
  | step (K : Set P) (u : P) (vs : List P)
      (fresh : u ∉ K) (simple : (u::vs).Nodup)
      (calibrated : ∀ p ∈ vs, p ∈ K) (ends : destination u vs = root)
      (next : Plan root (insert u K)) : Plan root K

def execute {root : P} {K : Set P} : Plan root K → (P → ℝ) → (P → ℝ) × List ℝ
  | .done _ _, x => (x,[])
  | .step _ u vs _ _ _ _ next, x =>
      let y := wordMap (route u vs) x
      let rest := execute next y
      (rest.1,y root::rest.2)

def means {root : P} {K : Set P} : Plan root K → ℕ
  | .done _ _ => 0
  | .step _ _ vs _ _ _ _ next => vs.length + means next

def phases {root : P} {K : Set P} : Plan root K → ℕ
  | .done _ _ => 0
  | .step _ _ _ _ _ _ _ next => 1 + phases next

theorem execute_sample_count {root : P} {K : Set P} (plan : Plan root K) (x : P → ℝ) :
    (execute plan x).2.length = phases plan := by
  induction plan generalizing x with
  | done => rfl
  | step K u vs hf hn hk he next ih => simpa [execute,phases,Nat.add_comm] using ih _

/-- Complete native tomography on a topology-certified calibration plan.
All new observations are at the one receiver; all updates are native means. -/
theorem execute_identifies {root : P} {K : Set P} (plan : Plan root K)
    (x y : P → ℝ) (hk : ∀ p ∈ K, x p = y p)
    (hs : (execute plan x).2 = (execute plan y).2) : x = y := by
  induction plan generalizing x y with
  | done K hc => exact funext (fun p => hk p (hc p))
  | step K u vs hf hn hv he next ih =>
    have hh : wordMap (route u vs) x root = wordMap (route u vs) y root ∧
        (execute next (wordMap (route u vs) x)).2 =
        (execute next (wordMap (route u vs) y)).2 := by simpa [execute] using hs
    have hp : wordMap (route u vs) x (destination u vs) =
        wordMap (route u vs) y (destination u vs) := by simpa [he] using hh.1
    have hnxt := phase_calibrates K u vs hn hv x y hk hp
    exact phase_preserves u vs hn x y (fun p h => hk p (hv p h)) (ih _ _ hnxt hh.2)

def completionSamples {root : P} (plan : Plan root {root}) (x : P → ℝ) : List ℝ :=
  x root :: (execute plan x).2

theorem completion_injective {root : P} (plan : Plan root {root}) :
    Function.Injective (completionSamples plan) := by
  intro x y h
  have hh : x root = y root ∧ (execute plan x).2 = (execute plan y).2 := by
    simpa [completionSamples] using h
  exact execute_identifies plan x y (fun p hp => by
    have he : p = root := hp
    simpa [he] using hh.1) hh.2

/-- Exact continuation criterion on an arbitrary preparation domain: past
records together with the remaining scalar state must distinguish the
meaning. A topology-certified native completion discharges sufficiency.
The criterion does not choose which record meanings must be preserved. -/
theorem completion_iff_preserved {X O Value : Type*} {root : P}
    (plan : Plan root {root}) (record : X → O) (state : X → P → ℝ) (meaning : X → Value) :
    (∀ x y, record x = record y →
      completionSamples plan (state x) = completionSamples plan (state y) →
      meaning x = meaning y) ↔
    (∀ x y, record x = record y → state x = state y → meaning x = meaning y) := by
  constructor
  · intro h x y hr hs
    exact h x y hr (congrArg (completionSamples plan) hs)
  · intro h x y hr hs
    exact h x y hr (completion_injective plan hs)

def lower {root : P} {K : Set P} : Plan root K → List (P × P)
  | .done _ _ => []
  | .step _ u vs _ _ _ _ next => route u vs ++ lower next

theorem wordMap_append (es fs : List (P × P)) (x : P → ℝ) :
    wordMap (es ++ fs) x = wordMap fs (wordMap es x) := by
  induction es generalizing x with
  | nil => rfl
  | cons e es ih => simpa only [List.cons_append,wordMap,LinearMap.comp_apply] using ih _

theorem execute_native {root : P} {K : Set P} (plan : Plan root K) (x : P → ℝ) :
    (execute plan x).1 = wordMap (lower plan) x := by
  induction plan generalizing x with
  | done => rfl
  | step K u vs hf hn hk he next ih => simp only [execute,lower,wordMap_append,ih]

theorem lower_length {root : P} {K : Set P} (plan : Plan root K) :
    (lower plan).length = means plan := by
  induction plan with
  | done => rfl
  | step K u vs hf hn hk he next ih => simp only [lower,means,List.length_append,route_length,ih]

def knownCount (K : Set P) : ℕ := K.ncard

theorem knownCount_insert [Fintype P] (K : Set P) (u : P) (hu : u ∉ K) :
    knownCount (insert u K) = knownCount K + 1 := by
  exact Set.ncard_insert_of_notMem hu

theorem phase_count [Fintype P] {root : P} {K : Set P} (plan : Plan root K) :
    phases plan + knownCount K = Fintype.card P := by
  classical
  induction plan with
  | done K hc =>
    have he : K = Set.univ := Set.eq_univ_of_forall hc
    simp [phases,knownCount,he,Nat.card_eq_fintype_card]
  | step K u vs hf hn hk he next ih =>
    rw [knownCount_insert K u hf] at ih
    simp only [phases]
    omega

/-- Each phase uses at most the current number of calibrated vertices;
summing that bound gives the triangular native work bound. -/
theorem triangular_work [Fintype P] {root : P} {K : Set P} (plan : Plan root K) :
    2 * means plan + knownCount K ^ 2 + phases plan ≤ Fintype.card P ^ 2 := by
  classical
  induction plan with
  | done K hc =>
    have he : K = Set.univ := Set.eq_univ_of_forall hc
    simp [means,phases,knownCount,he,Nat.card_eq_fintype_card]
  | step K u vs hf hn hk he next ih =>
    have hl : vs.length ≤ knownCount K := by
      have hsub : (vs.toFinset : Set P) ⊆ K := by
        intro p hp
        exact hk p (List.mem_toFinset.mp hp)
      have hc := Set.ncard_le_ncard hsub
      simpa only [Set.ncard_coe_finset,List.toFinset_card_of_nodup
        (List.nodup_cons.mp hn).2,knownCount] using hc
    rw [knownCount_insert K u hf] at ih
    simp only [means,phases]
    nlinarith

theorem complete_work_bound [Fintype P] {root : P} (plan : Plan root {root}) :
    phases plan + 1 = Fintype.card P ∧
      2 * means plan + Fintype.card P ≤ Fintype.card P ^ 2 := by
  classical
  have hc : knownCount ({root} : Set P) = 1 := by simp [knownCount]
  have hp := phase_count plan
  have hw := triangular_work plan
  rw [hc] at hp hw
  constructor
  · exact hp
  · nlinarith

/-- Finite source data generate the maximal continuation constraint through
the joint kernel of retained samples and the current physical state. -/
def joint {V I : Type*} [AddCommGroup V] [Module ℝ V]
    (L : I → V →ₗ[ℝ] ℝ) (state : V →ₗ[ℝ] (P → ℝ)) :
    I ⊕ P → V →ₗ[ℝ] ℝ :=
  Sum.elim L (fun p => (LinearMap.proj p).comp state)

theorem completion_iff_joint_span {V I : Type*} [AddCommGroup V] [Module ℝ V]
    [Fintype I] [Fintype P] {root : P} (plan : Plan root {root})
    (L : I → V →ₗ[ℝ] ℝ) (state : V →ₗ[ℝ] (P → ℝ)) (f : V →ₗ[ℝ] ℝ) :
    (∀ x y, (∀ i, L i x = L i y) →
      completionSamples plan (state x) = completionSamples plan (state y) → f x = f y) ↔
    f ∈ Submodule.span ℝ (Set.range (joint L state)) := by
  rw [← SourceTemporalObservation.observable_iff_span]
  constructor
  · intro h x y hxy
    apply h x y (fun i => hxy (.inl i))
    apply congrArg (completionSamples plan)
    exact funext (fun p => hxy (.inr p))
  · intro h x y hxy hs
    have he := completion_injective plan hs
    apply h x y
    intro i
    cases i with
    | inl i => exact hxy i
    | inr p => exact congrFun he p

end
end OPH.SourceTemporalTomography
