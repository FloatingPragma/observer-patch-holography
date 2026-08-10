import InformationProjection.LogTransitionAction

/-!
# The step-uniform reference as a gauge normal form

Issue B7 (#683).  The log-transition-action representation theorem
(`markov_path_law_eq_gibbs`) exhibits every strictly positive Markov path law
as the Gibbs tilt of the *declared* step-uniform reference.  This module
removes one layer of that declaredness by characterizing the reference
exactly: among row-stochastic kernels, the uniform kernel is the unique one
that is invariant under relabeling of transition targets and the unique one
with row-constant transition weight; among strictly positive row-stochastic
kernels it is also the unique one whose log-transition step action is
constant.  Each characterization is an iff with
the nondegeneracy premises stated, and a two-state nonuniform control shows
the invariance premise is load-bearing.

Consequently the (reference, action, multiplier) representation of the
committed law with a target-relabeling-invariant Markov reference is unique:
the reference must be `stepUniformRef` and the action is then the
log-transition action up to the already-characterized gauge orbit.

**Boundary.**  This is a representation-level normal-form theorem: it fixes
the reference gauge by an invariance property, exactly as a coordinate
convention is fixed.  It is not a source selection.  No theorem here claims
that the OPH source produces the relabeling-invariance principle, a physical
action scale, a clock, a current, or an interference rule; those remain the
open B7 attachments.
-/

namespace OPH.InformationProjection.ReferenceNormalForm

open OPH.InformationProjection

variable {Ω : Type*} [Fintype Ω]

/-- Row-stochasticity of a transition kernel. -/
def RowStochastic (P : Ω → Ω → ℝ) : Prop :=
  ∀ x, ∑ y : Ω, P x y = 1

/-- Invariance of a kernel under every relabeling of transition targets. -/
def TargetRelabelInvariant (P : Ω → Ω → ℝ) : Prop :=
  ∀ (σ : Equiv.Perm Ω) (x y : Ω), P x (σ y) = P x y

/-- The uniform kernel. -/
noncomputable def uniformKernel (Ω : Type*) [Fintype Ω] : Ω → Ω → ℝ :=
  fun _ _ => ((Fintype.card Ω : ℝ))⁻¹

omit [Fintype Ω] in
/-- Target-relabeling invariance forces each row to be constant. -/
theorem row_const_of_relabel_invariant (P : Ω → Ω → ℝ)
    (h : TargetRelabelInvariant P) (x y y' : Ω) : P x y = P x y' := by
  classical
  simpa using h (Equiv.swap y y') x y'

/-- **Reference normal form.**  A row-stochastic kernel is invariant under
every relabeling of transition targets iff it is the uniform kernel. -/
theorem relabel_invariant_iff_uniform [Nonempty Ω] (P : Ω → Ω → ℝ)
    (hrow : RowStochastic P) :
    TargetRelabelInvariant P ↔ P = uniformKernel Ω := by
  constructor
  · intro h
    funext x y
    have hconst : ∀ y', P x y' = P x y := fun y' =>
      row_const_of_relabel_invariant P h x y' y
    have hsum := hrow x
    rw [Finset.sum_congr rfl (fun y' _ => hconst y')] at hsum
    rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul] at hsum
    have hcard : (0 : ℝ) < (Fintype.card Ω : ℝ) := by
      exact_mod_cast Fintype.card_pos
    unfold uniformKernel
    field_simp
    linarith [hsum]
  · rintro rfl
    intro σ x y
    rfl

/-- Constant step weight also forces the uniform kernel. -/
theorem row_const_iff_uniform [Nonempty Ω] (P : Ω → Ω → ℝ)
    (hrow : RowStochastic P) :
    (∀ x y y', P x y = P x y') ↔ P = uniformKernel Ω := by
  constructor
  · intro hconst
    exact (relabel_invariant_iff_uniform P hrow).mp fun σ x y => hconst x (σ y) y
  · rintro rfl
    exact fun x y y' => rfl

/-- Constant log-transition step action forces the uniform kernel on
strictly positive kernels: the step actions `-log (P x y)` all agree iff the
kernel is uniform. -/
theorem constant_step_action_iff_uniform [Nonempty Ω] (P : Ω → Ω → ℝ)
    (hrow : RowStochastic P) (hpos : ∀ x y, 0 < P x y) :
    (∃ c : ℝ, ∀ x y, -Real.log (P x y) = c) ↔ P = uniformKernel Ω := by
  constructor
  · rintro ⟨c, hc⟩
    have hconst : ∀ x y y', P x y = P x y' := by
      intro x y y'
      have h1 := hc x y
      have h2 := hc x y'
      have hlog : Real.log (P x y) = Real.log (P x y') := by linarith
      exact Real.log_injOn_pos (Set.mem_Ioi.mpr (hpos x y))
        (Set.mem_Ioi.mpr (hpos x y')) hlog
    exact (row_const_iff_uniform P hrow).mp hconst
  · rintro rfl
    refine ⟨Real.log (Fintype.card Ω : ℝ), fun x y => ?_⟩
    unfold uniformKernel
    rw [Real.log_inv]
    ring

/-- The step-uniform reference of the representation theorem is exactly the
Markov path law of the unique relabeling-invariant row-stochastic kernel. -/
theorem stepUniformRef_is_normal_form [Nonempty Ω] (pi : Ω → ℝ) (n : ℕ)
    (P : Ω → Ω → ℝ) (hrow : RowStochastic P)
    (hinv : TargetRelabelInvariant P) :
    markovPathLaw pi P n = stepUniformRef pi n := by
  rw [(relabel_invariant_iff_uniform P hrow).mp hinv,
    stepUniformRef_eq_uniform_markov]
  rfl

/-- **Unique invariant-reference representation.**  Every
target-relabeling-invariant row-stochastic Markov reference has the
step-uniform path law; combined with `markov_path_law_eq_gibbs` and
`action_unique_up_to_gauge`, the (reference, action, multiplier) data of the
representation theorem carries no residual reference freedom beyond the
characterized action gauge orbit.  The invariance premise is the declared
gauge condition of the module boundary, not a source-derived law. -/
theorem unique_invariant_reference [Nonempty Ω] (pi : Ω → ℝ) (n : ℕ)
    (Q : Ω → Ω → ℝ) (hrow : RowStochastic Q)
    (hinv : TargetRelabelInvariant Q) :
    markovPathLaw pi Q n = stepUniformRef pi n :=
  stepUniformRef_is_normal_form pi n Q hrow hinv

/-! ## Negative control

A two-state kernel with distinct positive rows summing to one is
row-stochastic and strictly positive but neither relabeling-invariant nor of
constant step action: the invariance premise is load-bearing, not decorative. -/

/-- The biased two-state control kernel. -/
noncomputable def biasedKernel : Bool → Bool → ℝ :=
  fun _ y => if y then (2/3 : ℝ) else (1/3 : ℝ)

theorem biasedKernel_row_stochastic : RowStochastic biasedKernel := by
  intro x
  simp [biasedKernel]
  norm_num

theorem biasedKernel_pos : ∀ x y, 0 < biasedKernel x y := by
  intro x y
  by_cases h : y <;> simp [biasedKernel, h]

theorem biasedKernel_not_invariant : ¬ TargetRelabelInvariant biasedKernel := by
  intro h
  have := h (Equiv.swap true false) true true
  simp [Equiv.swap_apply_left, biasedKernel] at this
  norm_num at this

theorem biasedKernel_not_uniform : biasedKernel ≠ uniformKernel Bool := by
  intro h
  have := congrFun (congrFun h true) true
  simp [biasedKernel, uniformKernel] at this
  norm_num at this

#print axioms relabel_invariant_iff_uniform
#print axioms constant_step_action_iff_uniform
#print axioms unique_invariant_reference
#print axioms biasedKernel_not_invariant

end OPH.InformationProjection.ReferenceNormalForm
