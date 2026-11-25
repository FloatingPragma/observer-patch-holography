import ObserverPatchHolography.A2EndpointCommutator
import Mathlib.Tactic

/-! Position descent for the finite port-word carrier and angular moments of
a transverse analyzer. The frame, path and encounter probability are supplied. -/
namespace OPH.ScreenAngularTransfer

open ObserverPatchHolography.A2EndpointCommutator

abbrev Position := Fin 3 → ZMod 5

def translation (v : Position) : Equiv.Perm Position where
  toFun x := x+v
  invFun x := x-v
  left_inv x := add_sub_cancel_right x v
  right_inv x := sub_add_cancel x v

def step {F : Type*} (v : Position) (s : Position × F) : Position × F :=
  (s.1+v,s.2)

theorem descends {F : Type*} (v : Position) :
    Descends (Prod.fst : Position × F → Position) (step v) (translation v) :=
  fun _ => rfl

theorem all_diamonds (v w : Position) {F : Type*} (s : Position × F) :
    step v (step w s) = step w (step v s) := by
  simp only [step, add_right_comm]

/-- The endpoint theorem's surjective quotient can retain any internal fibre. -/
theorem source_fifteen (v : Axis → Position) :
    FifteenCommute (fun i => translation (v i)) := by
  apply quotient_fifteen_commute (q := (Prod.fst : Position × Unit → Position))
    (T := fun i => step (v i))
  · intro x; exact ⟨(x,()),rfl⟩
  · intro i; exact descends (v i)
  · intro i j _ x; rw [all_diamonds]

section Angular
variable {S : Type*} [Fintype S]

/-- The raw single-attempt intensity responds to the incident monopole and
directional second moment. The longitudinal rejection is outside this count. -/
theorem raw_moment (c mu : S → ℝ) :
    (∑ s, c s * ((1+(mu s)^2)/2)) =
      ((∑ s, c s) + (∑ s, c s*(mu s)^2))/2 := by
  simp_rw [show ∀ s, c s*((1+(mu s)^2)/2) = (c s+c s*(mu s)^2)/2 by intro s; ring]
  rw [← Finset.sum_div, Finset.sum_add_distrib]

theorem orbit_cancellation (c mu : S → ℝ) (h₀ : ∑ s, c s = 0)
    (h₂ : ∑ s, c s*(mu s)^2 = 0) :
    ∑ s, c s*((1+(mu s)^2)/2) = 0 := by rw [raw_moment,h₀,h₂]; norm_num

/-- Each isotropic orbit has second moment one third along the selected axis. -/
theorem isotropic_survival (c mu : S → ℝ) (h₀ : ∑ s, c s = 1)
    (h₂ : ∑ s, c s*(mu s)^2 = 1/3) :
    ∑ s, c s*((1+(mu s)^2)/2) = 2/3 := by rw [raw_moment,h₀,h₂]; norm_num

theorem conditional_gain (b : ℝ) (hb : 0 ≤ b) :
    b/(b+(1-b)*(2/3)) = 3*b/(2+b) := by
  have h : 0 < 2+b := by linarith
  have h' : 0 < b+(1-b)*(2/3) := by linarith
  apply (div_eq_div_iff (ne_of_gt h') (ne_of_gt h)).2
  ring
end Angular

#print axioms descends
#print axioms all_diamonds
#print axioms source_fifteen
#print axioms raw_moment
#print axioms orbit_cancellation
#print axioms isotropic_survival
#print axioms conditional_gain
end OPH.ScreenAngularTransfer
