import ObservableNormalForms.ObserverConfluence
import Mathlib.Tactic

/-! A routed strict write consumes one vacancy and one double and publishes
one mark. The native commutator tests and routing construction supply this
quotient map. -/
namespace RoutedCounts
open ObservableNormalForms Relation

structure State where
  vacancies : Nat
  doubles : Nat
  marks : Nat
  deriving DecidableEq

inductive Step : State → State → Prop
  | write (v d m : Nat) : Step ⟨v+1,d+1,m⟩ ⟨v,d,m+1⟩

def ready (x : State) := min x.vacancies x.doubles
def invariant (x : State) := (x.vacancies+x.marks, x.doubles+x.marks)
def consistent : Set State := {x | x.vacancies = 0 ∨ x.doubles = 0}

theorem preserves : ObservationPreserving Step invariant := by
  intro x y h
  cases h
  simp [invariant, Nat.add_right_comm]
  simp [Nat.add_assoc]

theorem strict_descent {x y : State} (h : Step x y) : ready y < ready x := by
  cases h
  simp [ready]

theorem complete : CompleteFor Step consistent := by
  rintro ⟨v,d,m⟩
  change (∀ z, ¬ Step ⟨v,d,m⟩ z) ↔ v = 0 ∨ d = 0
  constructor
  · intro h
    cases v with
    | zero => simp
    | succ v =>
      cases d with
      | zero => simp
      | succ d => exact False.elim (h _ (Step.write v d m))
  · rintro h z hz
    cases hz
    omega

theorem diamond {x y z : State} (a : Step x y) (b : Step x z) :
    ∃ w, ReflTransGen Step y w ∧ ReflTransGen Step z w := by
  cases a
  cases b
  exact ⟨_,.refl,.refl⟩

theorem distinguishes : BoundaryIdentifiesModulo consistent invariant Eq := by
  rintro ⟨v,d,m⟩ ⟨w,e,n⟩ a b h
  change v = 0 ∨ d = 0 at a
  change w = 0 ∨ e = 0 at b
  have hv := congrArg Prod.fst h
  have hd := congrArg Prod.snd h
  change v+m=w+n at hv
  change d+m=e+n at hd
  have hm : m=n := by omega
  cases hm
  congr 1 <;> omega

theorem unique_normal_endpoint : ObserverEndpointUniqueModulo Step invariant Eq :=
  (boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo preserves complete).mp distinguishes

theorem normalizes (x : State) : WeaklyNormalizing Step x := by
  obtain ⟨v,d,m⟩ := x
  induction v generalizing d m with
  | zero =>
    exact ⟨_, .refl, (complete _).mpr (Or.inl rfl)⟩
  | succ v ih =>
    cases d with
    | zero =>
      exact ⟨_, .refl, (complete _).mpr (Or.inr rfl)⟩
    | succ d =>
      obtain ⟨y, hy, hn⟩ := ih d (m+1)
      exact ⟨y, (ReflTransGen.single (Step.write v d m)).trans hy, hn⟩

/-- Global quadratic optimality among all ternary words of the same length and load. -/
theorem minimum_quadratic_iff (v s d : Nat) :
    (v = 0 ∨ d = 0) ↔
    ∀ v' s' d' : Nat, v'+s'+d' = v+s+d → s'+2*d' = s+2*d →
      s+4*d ≤ s'+4*d' := by
  constructor
  · intro h v' s' d' hn hl
    omega
  · intro h
    by_cases hv : v = 0
    · exact Or.inl hv
    · by_cases hd : d = 0
      · exact Or.inr hd
      · have hw := h (v-1) (s+2) (d-1) (by omega) (by omega)
        omega

#print axioms minimum_quadratic_iff
#print axioms preserves
#print axioms complete
#print axioms strict_descent
#print axioms diamond
#print axioms normalizes
#print axioms unique_normal_endpoint

end RoutedCounts
