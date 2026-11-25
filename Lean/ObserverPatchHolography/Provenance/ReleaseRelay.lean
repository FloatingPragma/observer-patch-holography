import ObservableNormalForms.SchedulerClassObstructions

/-! A measured two-register relay writes a network bit or clears a private
pending flag. Its entire row-support class inhabits the registered scheduler
interface. A separate recurrence describes the arrival sets after acknowledgements. -/
namespace OPH.ReleaseRelay

open ObservableNormalForms
open ObservableNormalForms.ProtectedObstructions
open scoped BigOperators
open ObservableNormalForms.ProtectedObstructions.SchedulerClass
open ObservableNormalForms.ProtectedObstructions.NonidentityExactness

def Menu (q r : State) : Prop :=
  r = q ∨ r = (!q.1, q.2) ∨ r = (q.1, false)

theorem menu_iff (q r : State) :
    Menu q r ↔ mismatch r ≤ mismatch q ∧ LocalMove q r := by
  obtain ⟨a,b⟩ := q
  obtain ⟨c,d⟩ := r
  cases a <;> cases b <;> cases c <;> cases d <;>
    norm_num [Menu, mismatch, LocalMove]

theorem cost_le {q r : State} (h : Menu q r) : stepCost q r ≤ 1 := by
  obtain ⟨a,b⟩ := q
  obtain ⟨c,d⟩ := r
  cases a <;> cases b <;> cases c <;> cases d <;>
    norm_num [Menu, stepCost, privateCost, sharedCost] at *

noncomputable def scheduler (K : FiniteMarkovKernel State)
    (h : ∀ q r, 0 < K.probability q r → Menu q r) : AdmissibleScheduler where
  kernel := K
  mismatch_nonincrease := fun q r hr => ((menu_iff q r).mp (h q r hr)).1
  locality := fun q r hr => ((menu_iff q r).mp (h q r hr)).2
  capacity_bound := by
    intro q
    calc
      _ ≤ ∑ r, K.probability q r := by
        apply Finset.sum_le_sum
        intro r _
        by_cases hp : 0 < K.probability q r
        · exact mul_le_of_le_one_right (K.probability_nonneg q r) (cost_le (h q r hp))
        · have hz := le_antisymm (le_of_not_gt hp) (K.probability_nonneg q r)
          simp [hz]
      _ = 1 := K.probability_sum_one q

theorem admissible_support (A : AdmissibleScheduler) {q r : State}
    (h : 0 < A.kernel.probability q r) : Menu q r :=
  (menu_iff q r).mpr ⟨A.mismatch_nonincrease q r h, A.locality q r h⟩

variable {V : Type*}

def sweep (adj : V → V → Prop) (lit : Set V) : Set V :=
  {v | v ∈ lit ∨ ∃ u ∈ lit, adj u v}

inductive Reached (adj : V → V → Prop) (seed : Set V) : ℕ → V → Prop where
  | start {v} (h : v ∈ seed) : Reached adj seed 0 v
  | wait {n v} (h : Reached adj seed n v) : Reached adj seed (n+1) v
  | cross {n u v} (h : Reached adj seed n u) (edge : adj u v) : Reached adj seed (n+1) v

theorem sweep_cone (adj : V → V → Prop) (seed : Set V) (n : ℕ) (v : V) :
    v ∈ (sweep adj)^[n] seed ↔ Reached adj seed n v := by
  induction n generalizing v with
  | zero => exact ⟨Reached.start, fun h => by cases h; assumption⟩
  | succ n ih =>
    rw [Function.iterate_succ_apply']
    constructor
    · rintro (h | ⟨u,h,hu⟩)
      · exact .wait ((ih v).mp h)
      · exact .cross ((ih u).mp h) hu
    · intro h
      cases h with
      | wait hp => exact Or.inl ((ih _).mpr hp)
      | cross hp he => exact Or.inr ⟨_,(ih _).mpr hp,he⟩

/-- A common two-valued preparation enters only through the differences
between its two conditional means at each probe. -/
theorem fair_cov (a b c d : ℝ) :
    (a*c + b*d)/2 - ((a+b)/2)*((c+d)/2) = (a-b)*(c-d)/4 := by ring

/-- A common deadline clips each increasing conditional mean. Its pairwise
difference product, used in the independent-copy covariance formula, can
only decrease. Both caps may differ. -/
theorem clipped_product {a b c d : ℝ} (hab : a ≤ b) (hcd : c ≤ d) (u v : ℝ) :
    0 ≤ (min b u - min a u)*(min d v - min c v) ∧
    (min b u - min a u)*(min d v - min c v) ≤ (b-a)*(d-c) := by
  have diff {x y : ℝ} (h : x ≤ y) (z : ℝ) :
      0 ≤ min y z - min x z ∧ min y z - min x z ≤ y-x := by
    rcases le_total y z with hy | hy
    · rw [min_eq_left hy, min_eq_left (h.trans hy)]
      constructor <;> linarith
    · rw [min_eq_right hy]
      rcases le_total x z with hx | hx
      · rw [min_eq_left hx]
        constructor <;> linarith
      · rw [min_eq_right hx]
        constructor <;> linarith
  obtain ⟨h₁,h₂⟩ := diff hab u
  obtain ⟨h₃,h₄⟩ := diff hcd v
  exact ⟨mul_nonneg h₁ h₃, mul_le_mul h₂ h₄ h₃ (sub_nonneg.mpr hab)⟩

#print axioms menu_iff
#print axioms cost_le
#print axioms scheduler
#print axioms admissible_support
#print axioms sweep_cone
#print axioms fair_cov
#print axioms clipped_product

end OPH.ReleaseRelay
