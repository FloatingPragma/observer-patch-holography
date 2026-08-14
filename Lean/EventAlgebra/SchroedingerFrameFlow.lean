import EventAlgebra.QuantumAdequacySurface

set_option autoImplicit false

namespace EventAlgebra

/-!
# The Schroedinger frame duality and the composed Born-Lueders statement

The composition module for the observation-ledger rows OL-C1 (Born
probabilities and Lueders conditioning) and OL-C2 (Schroedinger
dynamics), V3 issue #730, extending the surface of
`EventAlgebra.QuantumAdequacySurface` inside its namespaces.  Two
statements make the measurement surface and the dynamics surface share
one object, the represented frame state:

* `QuantumSurface.born_lueders_composed`: for a represented public
  frame (register rows PR-02 and PR-03) there is a unique state
  reproducing every frame weight in Born form AND carrying the Lueders
  receipts (conditioned state, repeatability, idempotence, composition
  for commuting events, order exchange) on the SAME state; and
  conditioning is closed: `QuantumSurface.conditionedFrame` packages
  the conditioned valuation as a represented public frame again, and
  `QuantumSurface.frameState_conditionedFrame` identifies its unique
  Busch-Gleason state with `luedersUpdate` of the original frame state.
* `QuantumSurface.schroedinger_frame_duality`: for the same frame
  bundled with a supplied pointwise continuous star-automorphism flow
  (register row PR-43, the field `flow` of `FlowedPublicFrame`), there
  is a self-adjoint generator, unique up to one additive real scalar,
  whose propagator conjugation represents the flow, and for every
  parameter value the Heisenberg-shifted frame valuation
  `E ↦ value (flow t E)` is again an additive effect valuation whose
  unique Busch-Gleason state is the propagator conjugate
  `U(-t) ρ U(t)` of the frame state: the Schroedinger-von Neumann
  evolution OF THE REPRESENTED STATE.

Composition boundaries, stated exactly: the flow is the supplied PR-43
premise, and no theorem here derives Schroedinger dynamics from source
data; no Hamiltonian is source-selected; the real parameter is a real
number, with no physical clock attached (register row PR-15 records
that open premise).  The frame representation and additivity stay the
declared rows PR-02 and PR-03.  Everything is finite-dimensional
matrix algebra.
-/

open Matrix
open scoped ComplexOrder

/-- **Register rows PR-02, PR-03, and PR-43: the flowed public
frame.**  The represented public frame of PR-02 and PR-03 (field
`frame`), bundled with a supplied real-parameter star-automorphism
group of the represented algebra, pointwise continuous in the
parameter.  The flow is the register row PR-43: a supplied premise,
produced by no theorem of this module.  The real parameter is not a
physical time; register row PR-15 records the missing clock
attachment. -/
structure FlowedPublicFrame (d : ℕ) where
  /-- The represented public frame (register rows PR-02 and PR-03). -/
  frame : RepresentedPublicFrame d
  /-- **Register row PR-43.**  The supplied star-automorphism flow of
  the represented algebra. -/
  flow : ℝ → (Matrix (Fin d) (Fin d) ℂ ≃⋆ₐ[ℂ] Matrix (Fin d) (Fin d) ℂ)
  /-- The flow starts at the identity. -/
  flow_zero : ∀ x, flow 0 x = x
  /-- The flow is a one-parameter group. -/
  flow_group : ∀ s t : ℝ, ∀ x, flow (s + t) x = flow s (flow t x)
  /-- The flow is pointwise continuous in the parameter. -/
  flow_continuous : ∀ x, Continuous fun t : ℝ => flow t x

namespace QuantumSurface

open OPH.Dynamics

variable {d : ℕ}

/-! ## OL-C1: the composed Born-Lueders statement -/

/-- **Born and Lueders on one state.**  Under PR-02 and PR-03 exactly
one state reproduces every frame weight in Born form, and that same
state carries the Lueders conditioning receipts: for every event of
nonvanishing Born weight, the conditioned matrix is a state assigning
its event weight one, conditioning is idempotent, sequential
conditioning on commuting events composes to the product event, and
commuting events condition in either order.  The uniqueness clause
consumes only the representation clauses; the conditioning clauses hold
for the represented state by the Lueders module. -/
theorem born_lueders_composed (F : RepresentedPublicFrame d) :
    ∃! ρ : Matrix (Fin d) (Fin d) ℂ,
      IsState ρ ∧
      (∀ E, IsEffect E → F.value E = (bornWeight ρ E).re) ∧
      (∀ P : Matrix (Fin d) (Fin d) ℂ, IsEvent P → bornWeight ρ P ≠ 0 →
        IsState (luedersUpdate ρ P) ∧
        bornWeight (luedersUpdate ρ P) P = 1 ∧
        luedersUpdate (luedersUpdate ρ P) P = luedersUpdate ρ P) ∧
      (∀ P Q : Matrix (Fin d) (Fin d) ℂ, IsEvent P → P * Q = Q * P →
        bornWeight ρ P ≠ 0 →
        luedersUpdate (luedersUpdate ρ P) Q = luedersUpdate ρ (P * Q)) ∧
      (∀ P Q : Matrix (Fin d) (Fin d) ℂ, IsEvent P → IsEvent Q →
        P * Q = Q * P → bornWeight ρ P ≠ 0 → bornWeight ρ Q ≠ 0 →
        luedersUpdate (luedersUpdate ρ P) Q =
          luedersUpdate (luedersUpdate ρ Q) P) := by
  refine ⟨frameState F,
    ⟨frameState_isState F, fun E hE => born_represents F hE,
      fun P hP hw => ⟨luedersUpdate_isState (frameState_isState F) hP hw,
        bornWeight_luedersUpdate_self hP hw, luedersUpdate_idem hP hw⟩,
      fun P Q hP hc hw =>
        luedersUpdate_luedersUpdate_of_commute hP hc hw,
      fun P Q hP hQ hc hwP hwQ =>
        luedersUpdate_comm hP hQ hc hwP hwQ⟩, ?_⟩
  rintro ρ' ⟨hs', hrep', -, -, -⟩
  exact state_ext_of_forall_effect_re hs' (frameState_isState F)
    fun E hE => by rw [← hrep' E hE, born_represents F hE]

/-- **Conditioning closure of the frame.**  The conditioned valuation
`E ↦ Re Tr(luedersUpdate ρ P · E)` of a represented public frame is
again a represented public frame: the PR-02/PR-03 bundle is closed
under Lueders conditioning on events of nonvanishing weight. -/
noncomputable def conditionedFrame (F : RepresentedPublicFrame d)
    {P : Matrix (Fin d) (Fin d) ℂ} (hP : IsEvent P)
    (hw : bornWeight (frameState F) P ≠ 0) : RepresentedPublicFrame d where
  value E := (bornWeight (luedersUpdate (frameState F) P) E).re
  nonneg {_E} hE := (isEffectValuation_born (lueders_state F hP hw)).nonneg hE
  le_one {_E} hE := (isEffectValuation_born (lueders_state F hP hw)).le_one hE
  map_one := (isEffectValuation_born (lueders_state F hP hw)).map_one
  additive {_E _E'} hE hE' hEE' :=
    (isEffectValuation_born (lueders_state F hP hw)).additive hE hE' hEE'

/-- **Born after Lueders.**  The unique Busch-Gleason state of the
conditioned frame IS the Lueders update of the original frame state:
the representation theorem and the conditioning rule commute on the
represented surface. -/
theorem frameState_conditionedFrame (F : RepresentedPublicFrame d)
    {P : Matrix (Fin d) (Fin d) ℂ} (hP : IsEvent P)
    (hw : bornWeight (frameState F) P ≠ 0) :
    frameState (conditionedFrame F hP hw) =
      luedersUpdate (frameState F) P := by
  refine state_ext_of_forall_effect_re
    (frameState_isState (conditionedFrame F hP hw))
    (lueders_state F hP hw) fun E hE => ?_
  have h := born_represents (conditionedFrame F hP hw) hE
  simpa [conditionedFrame] using h.symm

/-! ## OL-C2: the Schroedinger frame duality -/

section Propagator

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The adjoint of the Stone propagator is the propagator at the
negated parameter: `U(t)* = U(-t)`.  Glue between unitarity and the
one-parameter group law. -/
theorem stonePropagator_star (Ham : Matrix ι ι ℂ)
    (hHam : IsSelfAdjoint Ham) (t : ℝ) :
    star (stonePropagator Ham t) = stonePropagator Ham (-t) := by
  have hu := stonePropagator_mem_unitary Ham hHam t
  have h1 : star (stonePropagator Ham t) * stonePropagator Ham t = 1 :=
    (Unitary.mem_iff.mp hu).1
  have h2 : stonePropagator Ham t * stonePropagator Ham (-t) = 1 := by
    rw [← stonePropagator_add, add_neg_cancel, stonePropagator_zero]
  calc star (stonePropagator Ham t)
      = star (stonePropagator Ham t) *
          (stonePropagator Ham t * stonePropagator Ham (-t)) := by
        rw [h2, mul_one]
    _ = (star (stonePropagator Ham t) * stonePropagator Ham t) *
          stonePropagator Ham (-t) := by rw [mul_assoc]
    _ = stonePropagator Ham (-t) := by rw [h1, one_mul]

/-- The propagator at a parameter and at its negation multiply to the
identity. -/
theorem stonePropagator_mul_neg (Ham : Matrix ι ι ℂ) (t : ℝ) :
    stonePropagator Ham t * stonePropagator Ham (-t) = 1 := by
  rw [← stonePropagator_add, add_neg_cancel, stonePropagator_zero]

end Propagator

/-- **The Schroedinger frame duality.**  For a flowed public frame
(register rows PR-02, PR-03, PR-43) there is a self-adjoint generator
such that: the supplied flow is the propagator conjugation of that
generator; for every parameter value the Heisenberg-shifted frame
valuation `E ↦ value (flow t E)` is again an additive effect valuation,
the propagator conjugate `U(-t) ρ U(t)` of the frame state is a state,
and the unique Busch-Gleason state of the shifted valuation IS that
conjugate: the Schroedinger-von Neumann evolution of the represented
state; and the generator is unique up to one additive real scalar
multiple of the identity among self-adjoint generators representing the
same flow.  The flow stays the supplied PR-43 premise: nothing here
derives Schroedinger dynamics, no Hamiltonian is source-selected, and
the real parameter carries no physical clock (register row PR-15). -/
theorem schroedinger_frame_duality [NeZero d] (W : FlowedPublicFrame d) :
    ∃ Ham : Matrix (Fin d) (Fin d) ℂ, IsSelfAdjoint Ham ∧
      (∀ (t : ℝ) (x : Matrix (Fin d) (Fin d) ℂ),
        W.flow t x = stonePropagator Ham t * x * stonePropagator Ham (-t)) ∧
      (∀ t : ℝ,
        IsEffectValuation (fun E => W.frame.value (W.flow t E)) ∧
        IsState (stonePropagator Ham (-t) * frameState W.frame *
          stonePropagator Ham t) ∧
        gleasonState (fun E => W.frame.value (W.flow t E)) =
          stonePropagator Ham (-t) * frameState W.frame *
            stonePropagator Ham t) ∧
      (∀ H' : Matrix (Fin d) (Fin d) ℂ, IsSelfAdjoint H' →
        (∀ (t : ℝ) (x : Matrix (Fin d) (Fin d) ℂ),
          W.flow t x = stonePropagator H' t * x * stonePropagator H' (-t)) →
        ∃ r : ℝ, Ham = H' + (r : ℂ) • 1) := by
  haveI : Nonempty (Fin d) := ⟨⟨0, Nat.pos_of_ne_zero (NeZero.ne d)⟩⟩
  obtain ⟨Ham, hHam, hconj⟩ := finiteStoneConverse_of_continuous
    W.flow W.flow_zero W.flow_group W.flow_continuous
  refine ⟨Ham, hHam, hconj, ?_, ?_⟩
  · intro t
    set U : ℝ → Matrix (Fin d) (Fin d) ℂ := stonePropagator Ham with hU
    have hstar : ∀ s : ℝ, (U s)ᴴ = U (-s) := fun s => by
      rw [← Matrix.star_eq_conjTranspose]
      exact stonePropagator_star Ham hHam s
    -- The flow preserves effects: unitary conjugation of the interval.
    have heff : ∀ {E : Matrix (Fin d) (Fin d) ℂ}, IsEffect E →
        IsEffect (W.flow t E) := by
      intro E hE
      constructor
      · rw [hconj t E]
        have h := hE.1.mul_mul_conjTranspose_same (U t)
        rwa [hstar t] at h
      · have hone : (1 : Matrix (Fin d) (Fin d) ℂ) - W.flow t E =
            W.flow t (1 - E) := by rw [map_sub, map_one]
        rw [hone, hconj t (1 - E)]
        have h := hE.2.mul_mul_conjTranspose_same (U t)
        rwa [hstar t] at h
    -- The shifted valuation is an additive effect valuation.
    have hval : IsEffectValuation (fun E => W.frame.value (W.flow t E)) := by
      constructor
      · intro E hE
        exact W.frame.nonneg (heff hE)
      · intro E hE
        exact W.frame.le_one (heff hE)
      · show W.frame.value (W.flow t 1) = 1
        rw [map_one]
        exact W.frame.map_one
      · intro E E' hE hE' hEE'
        show W.frame.value (W.flow t (E + E')) =
          W.frame.value (W.flow t E) + W.frame.value (W.flow t E')
        rw [map_add]
        exact W.frame.additive (heff hE) (heff hE')
          (by rw [← map_add]; exact heff hEE')
    -- The propagator conjugate of the frame state is a state.
    have hstate : IsState (U (-t) * frameState W.frame * U t) := by
      constructor
      · have h := (frameState_isState W.frame).1.mul_mul_conjTranspose_same
          (U (-t))
        rwa [hstar (-t), neg_neg] at h
      · rw [Matrix.trace_mul_cycle, stonePropagator_mul_neg, one_mul]
        exact (frameState_isState W.frame).2
    refine ⟨hval, hstate, ?_⟩
    -- Both sides represent the shifted valuation; uniqueness closes.
    refine state_ext_of_forall_effect_re (gleasonState_isState hval)
      hstate fun E hE => ?_
    have h1 := gleasonState_represents hval hE
    rw [← h1]
    have h2 := born_represents W.frame (heff hE)
    rw [h2]
    have hbw : bornWeight (frameState W.frame) (W.flow t E) =
        bornWeight (U (-t) * frameState W.frame * U t) E := by
      rw [hconj t E]
      show (frameState W.frame * (U t * E * U (-t))).trace =
        (U (-t) * frameState W.frame * U t * E).trace
      calc (frameState W.frame * (U t * E * U (-t))).trace
          = ((frameState W.frame * U t * E) * U (-t)).trace := by
            rw [← mul_assoc, ← mul_assoc]
        _ = (U (-t) * (frameState W.frame * U t * E)).trace :=
            Matrix.trace_mul_comm _ _
        _ = (U (-t) * frameState W.frame * U t * E).trace := by
            rw [← mul_assoc, ← mul_assoc]
    rw [hbw]
  · intro H' hH' hflow'
    exact stoneFlow_generator_unique_up_to_real_scalar hHam hH'
      fun t x => by rw [← hconj t x, hflow' t x]

end QuantumSurface

end EventAlgebra

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms EventAlgebra.QuantumSurface.born_lueders_composed
#print axioms EventAlgebra.QuantumSurface.conditionedFrame
#print axioms EventAlgebra.QuantumSurface.frameState_conditionedFrame
#print axioms EventAlgebra.QuantumSurface.stonePropagator_star
#print axioms EventAlgebra.QuantumSurface.stonePropagator_mul_neg
#print axioms EventAlgebra.QuantumSurface.schroedinger_frame_duality
