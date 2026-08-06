import Variational.DiscreteNoether
import Variational.FiniteRealTransfer
import MixingChainRealization

namespace OPH.Variational

open OPH.Thermodynamics

/-!
# The finite Legendre bridge between the Lagrangian and Hamiltonian faces

The committed variational stack carries the Lagrangian face of the
recovered mechanics: `localAction` sums a two-point Lagrangian over
consecutive record pairs, `stationary_localAction_discreteEulerLagrange`
gives the discrete Euler-Lagrange equation at interior junctions, and the
discrete Noether theorem transports a constant segment current along
minimizing chains.  The Stone converse ladder of `Dynamics.StoneConverse`
forces the Hamiltonian face on the private block.  This module proves the
finite Legendre bridge between the two faces, in four exact pieces.

1.  The finite Legendre transform.  The momentum relation is the
    second-slot derivative packet of the two-point Lagrangian, in the
    committed `HasDerivAt` style, and the transform along a declared
    velocity solver is `H (x, p) = p * vel x p - L x (vel x p)`.  On the
    quadratic class `L (x, v) = a * v ^ 2 / 2 - V x` with `0 < a` the
    velocity slot is strictly convex, the momentum relation is linear
    with exact solver `p / a`, the transform evaluates to
    `p ^ 2 / (2 * a) + V x`, the transform applied twice returns the
    quadratic Lagrangian, and the Fenchel-Young inequality holds with
    equality exactly at the solved velocity.  At `a = 0` the transform
    collapses: no velocity solves the momentum relation at any nonzero
    momentum, every velocity solves it at momentum zero, and no velocity
    solver exists.

2.  The equivalence theorem.  At a junction of the committed stack, the
    discrete Euler-Lagrange sum vanishes exactly when the Legendre
    momenta satisfy the discrete Hamilton equations: the position update
    solves the velocity from the momentum, and the momentum update
    applies the force at the junction record.  The quadratic class gets
    the explicit flow map `(x, p) ↦ (x + p / a, p - V' (x + p / a))`,
    the strictly convex case is proved with the momentum map strictly
    monotone in the velocity slot and a declared solver section, both
    forms compose with the committed stationarity theorem, and the
    committed three-record witness chain of `FiniteRealTransfer`
    discharges every premise with exact values.

3.  The two-faces statement on the realized chain.  For a strictly
    positive transition kernel on any state type, a path minimizes the
    log-transition action within a family of paths exactly when it is a
    most probable path of that family.  The committed two-state mixing
    chain of `MixingChainRealization` realizes the statement in the
    committed `localAction` form: the log-transition action of a chain
    path equals the local action of one explicit two-point Lagrangian at
    the embedded records, so the single-site discrete Euler-Lagrange
    minimizers of that action are exactly the most probable paths of the
    realized chain.  The chain is reversible with exact receipts: its
    stationary edge flow is symmetric with cross value `535 / 61511`,
    and the chain is the row normalization of that symmetric flow by the
    stationary law.

4.  The Noether correspondence at witness level.  The committed constant
    chain current of the free translation witness equals the Legendre
    segment momentum on every segment, the transformed Hamiltonian of
    the free quadratic class is conserved by the discrete Hamilton flow,
    and the flow reproduces the committed witness chain from its initial
    record and momentum with exact conserved values `2` and `1`.

Claim boundary.  Finite exact statements over the committed packets: no
continuum limit and no physical unit attachment.  The selection of the
source action stays with issue B7 and the field-theoretic composition
with issue E9.  The two-faces statement identifies action minimizers
with most probable paths over the realized chain; a derivation that
produces the realized chain from a Hamiltonian flow is the open
remainder of the two-faces unification, and no theorem here claims it.
-/

variable {N : ℕ}

/-! ## The finite Legendre transform -/

/-- The momentum relation of the finite Legendre transform: the momentum
`p` is the second-slot derivative of the two-point Lagrangian at the pair
`(x, v)`, in the committed derivative-packet style of
`stationary_localAction_discreteEulerLagrange`. -/
def MomentumRelation (L : ℝ → ℝ → ℝ) (x v p : ℝ) : Prop :=
  HasDerivAt (fun w => L x w) p v

/-- The finite Legendre transform of a two-point Lagrangian along a
declared velocity solver: `H (x, p) = p * vel x p - L x (vel x p)`. -/
noncomputable def legendreTransform (L : ℝ → ℝ → ℝ) (vel : ℝ → ℝ → ℝ) :
    ℝ → ℝ → ℝ :=
  fun x p => p * vel x p - L x (vel x p)

/-- A velocity solver for `L`: at every base point and momentum, the
solved velocity satisfies the momentum relation. -/
def SolvesMomentum (L : ℝ → ℝ → ℝ) (vel : ℝ → ℝ → ℝ) : Prop :=
  ∀ x p, MomentumRelation L x (vel x p) p

/-- The quadratic Lagrangian class, velocity reading: stiffness `a` on
the second slot, potential `V` on the first.  Strict convexity in the
velocity slot is the condition `0 < a`. -/
noncomputable def quadLagrangian (a : ℝ) (V : ℝ → ℝ) (x v : ℝ) : ℝ :=
  a * v ^ 2 / 2 - V x

/-- The exact velocity solver of the quadratic class: the momentum
relation `a * v = p` is linear with solution `v = p / a`. -/
noncomputable def quadVelocitySolver (a : ℝ) : ℝ → ℝ → ℝ :=
  fun _x p => p / a

/-- The Legendre-transformed Hamiltonian of the quadratic class. -/
noncomputable def quadHamiltonian (a : ℝ) (V : ℝ → ℝ) (x p : ℝ) : ℝ :=
  p ^ 2 / (2 * a) + V x

/-- The second-slot derivative packet of the quadratic class: the
momentum at velocity `v` is `a * v`. -/
theorem quad_momentum_relation (a : ℝ) (V : ℝ → ℝ) (x v : ℝ) :
    MomentumRelation (quadLagrangian a V) x v (a * v) := by
  have h2 : HasDerivAt (fun w : ℝ => w ^ 2) (2 * v) v := by
    simpa using hasDerivAt_pow 2 v
  have h3 : HasDerivAt (fun w : ℝ => a * w ^ 2 / 2 - V x)
      (a * (2 * v) / 2) v :=
    ((h2.const_mul a).div_const 2).sub_const (V x)
  have hval : a * (2 * v) / 2 = a * v := by ring
  rw [hval] at h3
  exact h3

/-- Momenta of the quadratic class are unique: any packet value at
velocity `v` equals `a * v`. -/
theorem quad_momentum_unique (a : ℝ) (V : ℝ → ℝ) (x v p : ℝ)
    (h : MomentumRelation (quadLagrangian a V) x v p) : p = a * v :=
  HasDerivAt.unique h (quad_momentum_relation a V x v)

/-- Strict convexity of the velocity slot on the quadratic class: for
`0 < a` the map `v ↦ a * v ^ 2 / 2 - V x` is strictly convex on the
whole line. -/
theorem quadLagrangian_strictConvexOn (a : ℝ) (V : ℝ → ℝ) (ha : 0 < a)
    (x : ℝ) :
    StrictConvexOn ℝ Set.univ fun v : ℝ => quadLagrangian a V x v := by
  refine ⟨convex_univ, fun v _ w _ hvw b c hb hc hbc => ?_⟩
  simp only [quadLagrangian, smul_eq_mul]
  have hc' : c = 1 - b := by linarith
  subst hc'
  nlinarith [mul_pos (mul_pos hb hc)
    (pow_two_pos_of_ne_zero (sub_ne_zero.mpr hvw)), ha]

/-- The momentum map of the quadratic class is strictly monotone for
`0 < a`: the derivative form of strict convexity in the velocity
slot. -/
theorem quad_momentum_map_strictMono (a : ℝ) (ha : 0 < a) :
    StrictMono fun v : ℝ => a * v :=
  fun _ _ h => mul_lt_mul_of_pos_left h ha

/-- Under strict convexity the velocity solving the momentum relation is
unique: two solutions at the same momentum coincide. -/
theorem quad_velocity_unique (a : ℝ) (V : ℝ → ℝ) (ha : a ≠ 0)
    (x v₁ v₂ p : ℝ)
    (h₁ : MomentumRelation (quadLagrangian a V) x v₁ p)
    (h₂ : MomentumRelation (quadLagrangian a V) x v₂ p) : v₁ = v₂ := by
  have e₁ := quad_momentum_unique a V x v₁ p h₁
  have e₂ := quad_momentum_unique a V x v₂ p h₂
  exact mul_left_cancel₀ ha (e₁.symm.trans e₂)

/-- The exact solver solves: `p / a` satisfies the momentum relation of
the quadratic class at every base point and momentum. -/
theorem quad_solver_solves (a : ℝ) (V : ℝ → ℝ) (ha : a ≠ 0) :
    SolvesMomentum (quadLagrangian a V) (quadVelocitySolver a) := by
  intro x p
  have h := quad_momentum_relation a V x (p / a)
  have hval : a * (p / a) = p := by field_simp
  rwa [hval] at h

/-- **The finite Legendre transform of the quadratic class.**  Along the
exact solver, the transform of `a * v ^ 2 / 2 - V x` is
`p ^ 2 / (2 * a) + V x`. -/
theorem quad_legendreTransform (a : ℝ) (V : ℝ → ℝ) (ha : a ≠ 0) :
    legendreTransform (quadLagrangian a V) (quadVelocitySolver a)
      = quadHamiltonian a V := by
  funext x p
  simp only [legendreTransform, quadVelocitySolver, quadLagrangian,
    quadHamiltonian]
  field_simp
  ring

/-- The transformed Hamiltonian lies in the quadratic class with
inverted stiffness and negated potential. -/
theorem quad_hamiltonian_in_class (a : ℝ) (V : ℝ → ℝ) (ha : a ≠ 0) :
    quadHamiltonian a V = quadLagrangian a⁻¹ fun t => -V t := by
  funext x p
  simp only [quadHamiltonian, quadLagrangian]
  field_simp
  ring

/-- The solver of the transformed Hamiltonian solves its momentum
relation, so the second application of the transform runs through an
admissible solver. -/
theorem quad_hamiltonian_solver_solves (a : ℝ) (V : ℝ → ℝ) (ha : a ≠ 0) :
    SolvesMomentum (quadHamiltonian a V) (quadVelocitySolver a⁻¹) := by
  have h := quad_solver_solves a⁻¹ (fun t => -V t) (inv_ne_zero ha)
  rwa [← quad_hamiltonian_in_class a V ha] at h

/-- **Involutivity on the quadratic class.**  Applying the finite
Legendre transform twice, each time along the exact solver of the
transformed function, returns the quadratic Lagrangian. -/
theorem quad_legendre_involutive (a : ℝ) (V : ℝ → ℝ) (ha : a ≠ 0) :
    legendreTransform (quadHamiltonian a V) (quadVelocitySolver a⁻¹)
      = quadLagrangian a V := by
  rw [quad_hamiltonian_in_class a V ha,
    quad_legendreTransform a⁻¹ (fun t => -V t) (inv_ne_zero ha)]
  funext x v
  simp only [quadHamiltonian, quadLagrangian]
  field_simp
  ring

/-- **Fenchel-Young on the quadratic class.**  The transform dominates
every Legendre value: `p * w - L x w ≤ H x p` for every velocity `w`. -/
theorem quad_fenchel_young (a : ℝ) (V : ℝ → ℝ) (ha : 0 < a)
    (x p w : ℝ) :
    p * w - quadLagrangian a V x w ≤ quadHamiltonian a V x p := by
  have key : quadHamiltonian a V x p - (p * w - quadLagrangian a V x w)
      = a / 2 * (w - p / a) ^ 2 := by
    simp only [quadHamiltonian, quadLagrangian]
    field_simp
    ring
  have h0 : 0 ≤ a / 2 * (w - p / a) ^ 2 := by positivity
  linarith [key, h0]

/-- Fenchel-Young equality holds exactly at the solved velocity: the
maximizing `w` is the unique solution `p / a` of the momentum
relation. -/
theorem quad_fenchel_young_eq_iff (a : ℝ) (V : ℝ → ℝ) (ha : 0 < a)
    (x p w : ℝ) :
    p * w - quadLagrangian a V x w = quadHamiltonian a V x p
      ↔ w = p / a := by
  have key : quadHamiltonian a V x p - (p * w - quadLagrangian a V x w)
      = a / 2 * (w - p / a) ^ 2 := by
    simp only [quadHamiltonian, quadLagrangian]
    field_simp
    ring
  constructor
  · intro h
    have h0 : a / 2 * (w - p / a) ^ 2 = 0 := by linarith [key]
    have ha2 : a / 2 ≠ 0 := by positivity
    have h1 : (w - p / a) ^ 2 = 0 := (mul_eq_zero.mp h0).resolve_left ha2
    have h2 : w - p / a = 0 := pow_eq_zero_iff two_ne_zero |>.mp h1
    linarith [h2]
  · intro h
    subst h
    have h0 : quadHamiltonian a V x p
        - (p * (p / a) - quadLagrangian a V x (p / a)) = 0 := by
      rw [key]
      simp
    linarith [h0]

/-! ## The degenerate control `a = 0` -/

/-- At `a = 0` every velocity satisfies the momentum relation at
momentum zero: the momentum map is constant and the transform loses its
solved velocity. -/
theorem quad_degenerate_every_velocity (V : ℝ → ℝ) (x v : ℝ) :
    MomentumRelation (quadLagrangian 0 V) x v 0 := by
  have h := quad_momentum_relation 0 V x v
  rwa [zero_mul] at h

/-- At `a = 0` no velocity satisfies the momentum relation at any
nonzero momentum. -/
theorem quad_degenerate_no_solution (V : ℝ → ℝ) (x p : ℝ) (hp : p ≠ 0) :
    ¬ ∃ v : ℝ, MomentumRelation (quadLagrangian 0 V) x v p := by
  rintro ⟨v, hv⟩
  have h := quad_momentum_unique 0 V x v p hv
  rw [zero_mul] at h
  exact hp h

/-- **The degenerate collapse.**  At `a = 0` no velocity solver exists,
so the finite Legendre transform has no admissible solver input on the
degenerate quadratic class. -/
theorem quad_degenerate_no_solver (V : ℝ → ℝ) :
    ¬ ∃ vel : ℝ → ℝ → ℝ, SolvesMomentum (quadLagrangian 0 V) vel := by
  rintro ⟨vel, hvel⟩
  exact quad_degenerate_no_solution V 0 1 one_ne_zero ⟨vel 0 1, hvel 0 1⟩

/-! ## The equivalence theorem on the committed stack -/

/-- The quadratic two-point Lagrangian of the committed stack: the
velocity of a segment is its record difference, so the second slot
carries the later record.  `witnessL` of `FiniteRealTransfer` is the
instance `a = 2`, `V = 0`. -/
noncomputable def quadTwoPoint (a : ℝ) (V : ℝ → ℝ) (x y : ℝ) : ℝ :=
  quadLagrangian a V x (y - x)

/-- The committed free witness Lagrangian is the quadratic two-point
class at stiffness `2` and zero potential. -/
theorem witnessL_eq_quadTwoPoint :
    witnessL = quadTwoPoint 2 fun _ => 0 := by
  funext x y
  simp only [witnessL, quadTwoPoint, quadLagrangian]
  ring

/-- The Legendre momentum of segment `n`: the second-slot derivative of
the quadratic two-point Lagrangian on the segment, `a` times the record
difference.  With the unit translation generator this is the committed
`segmentMomentum` of the free class. -/
noncomputable def quadSegmentMomentum (a : ℝ) (γ : Fin (N + 1) → ℝ)
    (n : Fin N) : ℝ :=
  a * (γ n.succ - γ n.castSucc)

/-- One step of the discrete Hamilton flow of the quadratic class: the
position advances by the velocity solved from the momentum, and the
momentum advances by the force at the advanced position. -/
noncomputable def quadHamiltonStep (a : ℝ) (V' : ℝ → ℝ) (q : ℝ × ℝ) :
    ℝ × ℝ :=
  (q.1 + q.2 / a, q.2 - V' (q.1 + q.2 / a))

/-- Second-slot derivative packet of the quadratic two-point Lagrangian:
the segment momentum. -/
theorem quadTwoPoint_hasDerivAt_second (a : ℝ) (V : ℝ → ℝ) (x y : ℝ) :
    HasDerivAt (fun w => quadTwoPoint a V x w) (a * (y - x)) y := by
  have h1 : HasDerivAt (fun w : ℝ => w - x) 1 y :=
    (hasDerivAt_id y).sub_const x
  have h2 : HasDerivAt (fun w : ℝ => (w - x) ^ 2) (2 * (y - x)) y := by
    simpa using h1.pow 2
  have h3 : HasDerivAt (fun w : ℝ => a * (w - x) ^ 2 / 2 - V x)
      (a * (2 * (y - x)) / 2) y :=
    ((h2.const_mul a).div_const 2).sub_const (V x)
  have hval : a * (2 * (y - x)) / 2 = a * (y - x) := by ring
  rw [hval] at h3
  exact h3

/-- First-slot derivative packet of the quadratic two-point Lagrangian:
the negated segment momentum plus the force `-V'`. -/
theorem quadTwoPoint_hasDerivAt_first (a : ℝ) (V V' : ℝ → ℝ) (y x : ℝ)
    (hV : HasDerivAt V (V' x) x) :
    HasDerivAt (fun w => quadTwoPoint a V w y)
      (-(a * (y - x)) - V' x) x := by
  have h1 : HasDerivAt (fun w : ℝ => y - w) (-1) x :=
    (hasDerivAt_id x).const_sub y
  have h2 : HasDerivAt (fun w : ℝ => (y - w) ^ 2)
      (2 * (y - x) * -1) x := by
    simpa using h1.pow 2
  have h3 : HasDerivAt (fun w : ℝ => a * (y - w) ^ 2 / 2 - V w)
      (a * (2 * (y - x) * -1) / 2 - V' x) x :=
    ((h2.const_mul a).div_const 2).sub hV
  have hval : a * (2 * (y - x) * -1) / 2 - V' x
      = -(a * (y - x)) - V' x := by ring
  rw [hval] at h3
  have hfun : (fun w : ℝ => a * (y - w) ^ 2 / 2 - V w)
      = fun w => quadTwoPoint a V w y := by
    funext w
    simp only [quadTwoPoint, quadLagrangian]
  rwa [hfun] at h3

/-- **The equivalence theorem, quadratic class.**  At a junction of the
committed stack, the discrete Euler-Lagrange sum of the quadratic
two-point Lagrangian vanishes exactly when one step of the discrete
Hamilton flow carries the junction state `(position, momentum)` of the
incoming segment to the junction state of the outgoing segment: the
position update solves the velocity from the momentum, and the momentum
update applies the force at the junction record. -/
theorem quad_discreteEulerLagrange_iff_hamiltonStep (a : ℝ) (ha : a ≠ 0)
    (V' : ℝ → ℝ) (γ : Fin (N + 1) → ℝ) {k m : Fin N}
    (hkm : k.succ = m.castSucc) :
    quadSegmentMomentum a γ k
        + (-(quadSegmentMomentum a γ m) - V' (γ m.castSucc)) = 0
      ↔ quadHamiltonStep a V' (γ k.castSucc, quadSegmentMomentum a γ k)
          = (γ m.castSucc, quadSegmentMomentum a γ m) := by
  have hpt : γ k.succ = γ m.castSucc := by rw [hkm]
  have hpos : γ k.castSucc + quadSegmentMomentum a γ k / a
      = γ m.castSucc := by
    rw [quadSegmentMomentum, mul_div_cancel_left₀ _ ha]
    rw [← hpt]
    ring
  have hstep : quadHamiltonStep a V'
      (γ k.castSucc, quadSegmentMomentum a γ k)
      = (γ m.castSucc,
          quadSegmentMomentum a γ k - V' (γ m.castSucc)) := by
    unfold quadHamiltonStep
    rw [hpos]
  rw [hstep, Prod.mk.injEq]
  constructor
  · intro h
    exact ⟨rfl, by linarith⟩
  · rintro ⟨-, h⟩
    linarith

/-- **Equivalence composed with the committed stack, quadratic class.**
A history that minimizes the local action of the quadratic two-point
Lagrangian among single-site variations at a junction satisfies the
discrete Hamilton equations there: the flow step carries the incoming
junction state to the outgoing junction state.  This composes
`stationary_localAction_discreteEulerLagrange` with the derivative
packets of the quadratic class and the junction equivalence. -/
theorem quad_stationary_hamiltonStep (a : ℝ) (ha : a ≠ 0)
    (V V' : ℝ → ℝ) (hV : ∀ t, HasDerivAt V (V' t) t)
    (γ : Fin (N + 1) → ℝ) {k m : Fin N} (hkm : k.succ = m.castSucc)
    (hmin : ∀ x, localAction (quadTwoPoint a V) γ
      ≤ localAction (quadTwoPoint a V) (Function.update γ k.succ x)) :
    quadHamiltonStep a V' (γ k.castSucc, quadSegmentMomentum a γ k)
      = (γ m.castSucc, quadSegmentMomentum a γ m) := by
  have h₂ := quadTwoPoint_hasDerivAt_second a V (γ k.castSucc) (γ k.succ)
  have h₁ := quadTwoPoint_hasDerivAt_first a V V' (γ m.succ)
    (γ m.castSucc) (hV (γ m.castSucc))
  have hEL := stationary_localAction_discreteEulerLagrange
    (quadTwoPoint a V) γ hkm _ _ h₂ h₁ hmin
  exact (quad_discreteEulerLagrange_iff_hamiltonStep a ha V' γ hkm).mp hEL

/-! ## The general strictly convex case -/

/-- A solver section of a strictly monotone momentum map is a two-sided
inverse: solving after evaluating returns the velocity. -/
theorem solver_section_retracts (d₂ vel : ℝ → ℝ → ℝ)
    (hmono : ∀ x, StrictMono fun v => d₂ x v)
    (hsec : ∀ x p, d₂ x (vel x p) = p) (x v : ℝ) :
    vel x (d₂ x v) = v :=
  (hmono x).injective (hsec x (d₂ x v))

/-- **The equivalence theorem, strictly convex case.**  For momentum
maps `d₁`, `d₂` with a strictly monotone second slot and a declared
solver section, the discrete Euler-Lagrange sum at a junction vanishes
exactly when the discrete Hamilton equations hold there: the junction
record is the velocity solved from the incoming momentum, and the
incoming momentum balances the first-slot force evaluated at the
velocity solved from the outgoing momentum. -/
theorem discreteEulerLagrange_iff_hamilton_strictConvex
    (d₁ d₂ vel : ℝ → ℝ → ℝ)
    (hmono : ∀ x, StrictMono fun v => d₂ x v)
    (hsec : ∀ x p, d₂ x (vel x p) = p)
    (γ : Fin (N + 1) → ℝ) {k m : Fin N} (hkm : k.succ = m.castSucc) :
    d₂ (γ k.castSucc) (γ k.succ) + d₁ (γ m.castSucc) (γ m.succ) = 0
      ↔ (γ m.castSucc
            = vel (γ k.castSucc) (d₂ (γ k.castSucc) (γ k.succ))
          ∧ d₂ (γ k.castSucc) (γ k.succ)
              + d₁ (γ m.castSucc)
                  (vel (γ m.castSucc)
                    (d₂ (γ m.castSucc) (γ m.succ))) = 0) := by
  have hpt : γ k.succ = γ m.castSucc := by rw [hkm]
  have hk := solver_section_retracts d₂ vel hmono hsec (γ k.castSucc)
    (γ k.succ)
  have hm := solver_section_retracts d₂ vel hmono hsec (γ m.castSucc)
    (γ m.succ)
  rw [hk, hm]
  constructor
  · intro h
    exact ⟨hpt.symm, h⟩
  · rintro ⟨-, h⟩
    exact h

/-- **Strictly convex case composed with the committed stack.**  A
history that minimizes the local action among single-site variations at
a junction, for a differentiable two-point Lagrangian whose second-slot
momentum map is strictly monotone with a declared solver section,
satisfies the discrete Hamilton equations at that junction. -/
theorem stationary_hamilton_strictConvex
    (L : ℝ → ℝ → ℝ) (d₁ d₂ vel : ℝ → ℝ → ℝ)
    (hL : ∀ p : ℝ × ℝ, HasFDerivAt (Function.uncurry L)
      ((d₁ p.1 p.2) • (ContinuousLinearMap.fst ℝ ℝ ℝ)
        + (d₂ p.1 p.2) • (ContinuousLinearMap.snd ℝ ℝ ℝ)) p)
    (hmono : ∀ x, StrictMono fun v => d₂ x v)
    (hsec : ∀ x p, d₂ x (vel x p) = p)
    (γ : Fin (N + 1) → ℝ) {k m : Fin N} (hkm : k.succ = m.castSucc)
    (hmin : ∀ x, localAction L γ
      ≤ localAction L (Function.update γ k.succ x)) :
    γ m.castSucc = vel (γ k.castSucc) (d₂ (γ k.castSucc) (γ k.succ))
      ∧ d₂ (γ k.castSucc) (γ k.succ)
          + d₁ (γ m.castSucc)
              (vel (γ m.castSucc) (d₂ (γ m.castSucc) (γ m.succ))) = 0 := by
  have hd₂ : HasDerivAt (fun x => L (γ k.castSucc) x)
      (d₂ (γ k.castSucc) (γ k.succ)) (γ k.succ) := by
    have hL' := hL (γ k.castSucc, γ k.succ)
    have hg : HasDerivAt (fun x => ((γ k.castSucc, x) : ℝ × ℝ))
        ((0, 1) : ℝ × ℝ) (γ k.succ) :=
      (hasDerivAt_const _ _).prodMk (hasDerivAt_id _)
    have h := hL'.comp_hasDerivAt (γ k.succ) hg
    simpa [ContinuousLinearMap.add_apply,
      ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst',
      ContinuousLinearMap.coe_snd', smul_eq_mul,
      Function.uncurry] using h
  have hd₁ : HasDerivAt (fun x => L x (γ m.succ))
      (d₁ (γ m.castSucc) (γ m.succ)) (γ m.castSucc) := by
    have hL' := hL (γ m.castSucc, γ m.succ)
    have hg : HasDerivAt (fun x => ((x, γ m.succ) : ℝ × ℝ))
        ((1, 0) : ℝ × ℝ) (γ m.castSucc) :=
      (hasDerivAt_id _).prodMk (hasDerivAt_const _ _)
    have h := hL'.comp_hasDerivAt (γ m.castSucc) hg
    simpa [ContinuousLinearMap.add_apply,
      ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst',
      ContinuousLinearMap.coe_snd', smul_eq_mul,
      Function.uncurry] using h
  have hEL := stationary_localAction_discreteEulerLagrange
    L γ hkm _ _ hd₂ hd₁ hmin
  exact (discreteEulerLagrange_iff_hamilton_strictConvex d₁ d₂ vel hmono
    hsec γ hkm).mp hEL

/-- The quadratic two-point class realizes the strictly convex
hypotheses: for `0 < a` its momentum map is strictly monotone in the
second slot, and `vel x p = x + p / a` is a solver section. -/
theorem quadTwoPoint_realizes_strictConvex (a : ℝ) (ha : 0 < a) :
    (∀ x : ℝ, StrictMono fun y : ℝ => a * (y - x))
      ∧ ∀ x p : ℝ, a * (x + p / a - x) = p := by
  constructor
  · intro x y z h
    exact mul_lt_mul_of_pos_left (by linarith) ha
  · intro x p
    field_simp
    ring

/-! ## Exact witness on the committed three-site chain -/

/-- **The equivalence witness on the committed chain.**  The affine
minimizer of `FiniteRealTransfer` satisfies the discrete Hamilton
equations at its interior junction: the premise is the committed real
local minimality `witness_real_local_min`, transported through the
identification of `witnessL` with the quadratic class, and the
conclusion is one exact step of the discrete Hamilton flow. -/
theorem witness_stationary_hamiltonStep :
    quadHamiltonStep 2 (fun _ => 0)
        (witnessEmb 0 (0 : Fin 2).castSucc,
          quadSegmentMomentum 2 (witnessEmb 0) 0)
      = (witnessEmb 0 (1 : Fin 2).castSucc,
          quadSegmentMomentum 2 (witnessEmb 0) 1) := by
  refine quad_stationary_hamiltonStep 2 two_ne_zero (fun _ => 0)
    (fun _ => 0) (fun t => hasDerivAt_const t 0) (witnessEmb 0)
    (k := 0) (m := 1) (by decide) ?_
  intro x
  have h := witness_real_local_min x
  rwa [witnessL_eq_quadTwoPoint] at h

/-- Exact values of the witness bridge: both segment momenta equal `2`,
the flow step carries `(0, 2)` to `(1, 2)`, and the transformed
Hamiltonian reads `1` on the witness state. -/
theorem witness_hamilton_values :
    quadSegmentMomentum 2 (witnessEmb 0) 0 = 2
      ∧ quadSegmentMomentum 2 (witnessEmb 0) 1 = 2
      ∧ quadHamiltonStep 2 (fun _ => 0) ((0 : ℝ), (2 : ℝ)) = (1, 2)
      ∧ quadHamiltonian 2 (fun _ => 0) 0 2 = 1 := by
  have e0 : witnessEmb 0 (0 : Fin 3) = 0 := by simp [witnessEmb]
  have e1 : witnessEmb 0 (1 : Fin 3) = 1 := by simp [witnessEmb]
  have e2 : witnessEmb 0 (2 : Fin 3) = 2 := by simp [witnessEmb]
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [quadSegmentMomentum, Fin.succ_zero_eq_one, Fin.castSucc_zero,
      e0, e1]
    norm_num
  · rw [quadSegmentMomentum, Fin.succ_one_eq_two, Fin.castSucc_one,
      e1, e2]
    norm_num
  · norm_num [quadHamiltonStep, Prod.mk.injEq]
  · norm_num [quadHamiltonian]

/-! ## The two-faces statement on the realized chain -/

/-- The transition weight of a path through a kernel: the product of
the kernel over consecutive record pairs, the conditional path
probability given the initial record. -/
noncomputable def pathTransitionWeight {Ω : Type*} (P : Ω → Ω → ℝ)
    (s : Fin (N + 1) → Ω) : ℝ :=
  ∏ n : Fin N, P (s n.castSucc) (s n.succ)

/-- The log-transition action of a path: the sum of the negated kernel
logarithms over consecutive record pairs. -/
noncomputable def logTransitionAction {Ω : Type*} (P : Ω → Ω → ℝ)
    (s : Fin (N + 1) → Ω) : ℝ :=
  ∑ n : Fin N, -Real.log (P (s n.castSucc) (s n.succ))

/-- For a strictly positive kernel the path weight is the exponential of
the negated log-transition action. -/
theorem pathTransitionWeight_eq_exp_neg {Ω : Type*} (P : Ω → Ω → ℝ)
    (hP : ∀ i j, 0 < P i j) (s : Fin (N + 1) → Ω) :
    pathTransitionWeight P s
      = Real.exp (-logTransitionAction P s) := by
  unfold pathTransitionWeight logTransitionAction
  have hflip : -(∑ n : Fin N, -Real.log (P (s n.castSucc) (s n.succ)))
      = ∑ n : Fin N, Real.log (P (s n.castSucc) (s n.succ)) := by
    simp
  rw [hflip, Real.exp_sum]
  exact Finset.prod_congr rfl fun n _ => (Real.exp_log (hP _ _)).symm

/-- **The two-faces comparison.**  For a strictly positive kernel, one
path has log-transition action at or below another exactly when its
transition weight is at or above the other's. -/
theorem logTransitionAction_le_iff {Ω : Type*} (P : Ω → Ω → ℝ)
    (hP : ∀ i j, 0 < P i j) (s t : Fin (N + 1) → Ω) :
    logTransitionAction P s ≤ logTransitionAction P t
      ↔ pathTransitionWeight P t ≤ pathTransitionWeight P s := by
  rw [pathTransitionWeight_eq_exp_neg P hP,
    pathTransitionWeight_eq_exp_neg P hP, Real.exp_le_exp,
    neg_le_neg_iff]

/-- **Minimizers are most probable paths, single-site form.**  A path
minimizes the log-transition action among all single-site variations
exactly when it maximizes the transition weight among them.  The
single-site quantifier is the discrete Euler-Lagrange minimality of the
committed stack, read on the chain's path space. -/
theorem logTransitionAction_site_min_iff {Ω : Type*} (P : Ω → Ω → ℝ)
    (hP : ∀ i j, 0 < P i j) (s : Fin (N + 1) → Ω) :
    (∀ (i : Fin (N + 1)) (x : Ω), logTransitionAction P s
        ≤ logTransitionAction P (Function.update s i x))
      ↔ ∀ (i : Fin (N + 1)) (x : Ω),
          pathTransitionWeight P (Function.update s i x)
            ≤ pathTransitionWeight P s :=
  forall₂_congr fun _ _ => logTransitionAction_le_iff P hP s _

/-- Minimizers are most probable paths, global form: minimality of the
log-transition action over every path of the same length is maximality
of the transition weight over every path of the same length. -/
theorem logTransitionAction_global_min_iff {Ω : Type*} (P : Ω → Ω → ℝ)
    (hP : ∀ i j, 0 < P i j) (s : Fin (N + 1) → Ω) :
    (∀ t : Fin (N + 1) → Ω, logTransitionAction P s
        ≤ logTransitionAction P t)
      ↔ ∀ t : Fin (N + 1) → Ω,
          pathTransitionWeight P t ≤ pathTransitionWeight P s :=
  forall_congr' fun t => logTransitionAction_le_iff P hP s t

/-! ## The committed realized chain -/

/-- The real transition weights of the committed two-state mixing
chain. -/
noncomputable def chainWeight (i j : Fin 2) : ℝ := (mixingChain i j : ℝ)

/-- Every transition weight of the committed chain is strictly
positive. -/
theorem chainWeight_pos (i j : Fin 2) : 0 < chainWeight i j := by
  have h := mixingChain_entries_pos i j
  unfold chainWeight
  exact_mod_cast h

/-- The embedding of chain paths into real records: state `i` sits at
the real value `i`. -/
def chainEmb (s : Fin (N + 1) → Fin 2) : Fin (N + 1) → ℝ :=
  fun n => ((s n : ℕ) : ℝ)

/-- The log-transition two-point Lagrangian of the committed chain: the
bilinear interpolation of the negated log-weight table through the
embedded states. -/
noncomputable def chainLogLagrangian (x y : ℝ) : ℝ :=
  (1 - x) * (1 - y) * -Real.log (chainWeight 0 0)
    + (1 - x) * y * -Real.log (chainWeight 0 1)
    + x * (1 - y) * -Real.log (chainWeight 1 0)
    + x * y * -Real.log (chainWeight 1 1)

/-- Corner evaluation: at embedded states the bilinear Lagrangian reads
the negated log-weight table exactly. -/
theorem chainLogLagrangian_corner (i j : Fin 2) :
    chainLogLagrangian ((i : ℕ) : ℝ) ((j : ℕ) : ℝ)
      = -Real.log (chainWeight i j) := by
  fin_cases i <;> fin_cases j <;>
    norm_num [chainLogLagrangian]

/-- **The log-transition action is a committed local action.**  On
embedded chain paths, the local action of the bilinear log-transition
Lagrangian equals the log-transition action of the chain path.  This
puts the realized chain's path cost in the exact form of the committed
variational stack. -/
theorem localAction_chainEmb (s : Fin (N + 1) → Fin 2) :
    localAction chainLogLagrangian (chainEmb s)
      = logTransitionAction chainWeight s := by
  unfold localAction logTransitionAction chainEmb
  exact Finset.sum_congr rfl fun n _ =>
    chainLogLagrangian_corner (s n.castSucc) (s n.succ)

/-- **The two-faces statement on the committed chain, single-site
form.**  A chain path minimizes the committed-form local action of the
log-transition Lagrangian among all single-site chain variations exactly
when it is a most probable path among them. -/
theorem two_faces_mixing_chain (s : Fin (N + 1) → Fin 2) :
    (∀ (i : Fin (N + 1)) (x : Fin 2),
        localAction chainLogLagrangian (chainEmb s)
          ≤ localAction chainLogLagrangian
              (chainEmb (Function.update s i x)))
      ↔ ∀ (i : Fin (N + 1)) (x : Fin 2),
          pathTransitionWeight chainWeight (Function.update s i x)
            ≤ pathTransitionWeight chainWeight s := by
  simp only [localAction_chainEmb]
  exact logTransitionAction_site_min_iff chainWeight chainWeight_pos s

/-- The two-faces statement on the committed chain, global form. -/
theorem two_faces_mixing_chain_global (s : Fin (N + 1) → Fin 2) :
    (∀ t : Fin (N + 1) → Fin 2,
        localAction chainLogLagrangian (chainEmb s)
          ≤ localAction chainLogLagrangian (chainEmb t))
      ↔ ∀ t : Fin (N + 1) → Fin 2,
          pathTransitionWeight chainWeight t
            ≤ pathTransitionWeight chainWeight s := by
  simp only [localAction_chainEmb]
  exact logTransitionAction_global_min_iff chainWeight chainWeight_pos s

/-! ## Reversibility receipts of the committed chain -/

/-- The stationary edge flow of the committed chain: the stationary law
times the transition matrix. -/
def mixingFlow (i j : Fin 2) : ℚ :=
  mixingChainStationary i * mixingChain i j

/-- **Detailed balance, exact.**  The stationary edge flow of the
committed chain is symmetric: the chain is reversible. -/
theorem mixingFlow_symmetric (i j : Fin 2) :
    mixingFlow i j = mixingFlow j i := by
  fin_cases i <;> fin_cases j <;>
    norm_num [mixingFlow, mixingChainStationary, mixingChain]

/-- The exact cross value of the stationary edge flow. -/
theorem mixingFlow_cross_value : mixingFlow 0 1 = 535 / 61511 := by
  norm_num [mixingFlow, mixingChainStationary, mixingChain]

/-- The committed chain is the row normalization of its symmetric
stationary edge flow: the generator reading of reversibility at the
finite exact level. -/
theorem mixingChain_eq_flow_normalized (i j : Fin 2) :
    mixingChain i j = mixingFlow i j / mixingChainStationary i := by
  fin_cases i <;> fin_cases j <;>
    norm_num [mixingFlow, mixingChainStationary, mixingChain]

/-! ## The Noether correspondence at witness level -/

/-- The committed Noether chain current of the free translation witness
equals the Legendre segment momentum of the quadratic class at stiffness
`2`, on every path and segment. -/
theorem noether_current_eq_legendre_momentum (γ : Fin (N + 1) → ℝ)
    (n : Fin N) :
    segmentMomentum freeDTwo unitTranslationGenerator γ n
      = quadSegmentMomentum 2 γ n := by
  simp [segmentMomentum, freeDTwo, unitTranslationGenerator,
    quadSegmentMomentum]

/-- The discrete Hamilton flow of the free quadratic class conserves the
momentum readout: the Legendre image of the translation Noether
current. -/
theorem free_hamiltonStep_conserves_momentum (a : ℝ) (q : ℝ × ℝ) :
    (quadHamiltonStep a (fun _ => 0) q).2 = q.2 := by
  simp [quadHamiltonStep]

/-- The discrete Hamilton flow of the free quadratic class conserves the
transformed Hamiltonian readout. -/
theorem free_hamiltonStep_conserves_hamiltonian (a : ℝ) (q : ℝ × ℝ) :
    quadHamiltonian a (fun _ => 0) (quadHamiltonStep a (fun _ => 0) q).1
        (quadHamiltonStep a (fun _ => 0) q).2
      = quadHamiltonian a (fun _ => 0) q.1 q.2 := by
  simp [quadHamiltonStep, quadHamiltonian]

/-- The Legendre transform of the free witness class: the instance
`a = 2`, `V = 0` of the quadratic transform. -/
theorem witness_legendre_hamiltonian :
    legendreTransform (quadLagrangian 2 fun _ => 0)
        (quadVelocitySolver 2)
      = quadHamiltonian 2 fun _ => 0 :=
  quad_legendreTransform 2 (fun _ => 0) two_ne_zero

/-- **The Noether correspondence on the quadratic witness.**  On the
committed affine free path, the constant chain current equals the
Legendre segment momentum on every segment, the discrete Hamilton flow
reproduces the witness chain from its initial record and momentum, and
the momentum and Hamiltonian readouts are conserved with exact values
`2` and `1`. -/
theorem witness_noether_correspondence :
    (∀ n : Fin 2,
        segmentMomentum freeDTwo unitTranslationGenerator
            affineFreePath3 n
          = quadSegmentMomentum 2 affineFreePath3 n)
      ∧ quadHamiltonStep 2 (fun _ => 0)
          (affineFreePath3 0, quadSegmentMomentum 2 affineFreePath3 0)
          = (affineFreePath3 1, quadSegmentMomentum 2 affineFreePath3 1)
      ∧ quadHamiltonStep 2 (fun _ => 0)
          (affineFreePath3 1, quadSegmentMomentum 2 affineFreePath3 1)
          = (affineFreePath3 2, quadSegmentMomentum 2 affineFreePath3 2)
      ∧ quadSegmentMomentum 2 affineFreePath3 0 = 2
      ∧ quadHamiltonian 2 (fun _ => 0) (affineFreePath3 0)
          (quadSegmentMomentum 2 affineFreePath3 0) = 1
      ∧ quadHamiltonian 2 (fun _ => 0) (affineFreePath3 1)
          (quadSegmentMomentum 2 affineFreePath3 1) = 1 := by
  refine ⟨fun n => noether_current_eq_legendre_momentum _ n,
    ?_, ?_, ?_, ?_, ?_⟩ <;>
    norm_num [quadSegmentMomentum, quadHamiltonStep, quadHamiltonian,
      affineFreePath3, Prod.mk.injEq, Fin.castSucc_zero,
      Fin.succ_zero_eq_one, Fin.castSucc_one, Fin.succ_one_eq_two]

/-! ## Composite receipt -/

/-- **The Legendre bridge receipt.**  The committed chain's stationary
edge flow is symmetric, the single-site minimizers of the committed-form
log-transition action are exactly the most probable paths of the
realized chain at every length, and the committed witness minimizer
satisfies one exact step of the discrete Hamilton flow at its interior
junction. -/
theorem legendreBridge_receipt :
    (∀ i j, mixingFlow i j = mixingFlow j i)
      ∧ (∀ (M : ℕ) (s : Fin (M + 1) → Fin 2),
          (∀ (i : Fin (M + 1)) (x : Fin 2),
              localAction chainLogLagrangian (chainEmb s)
                ≤ localAction chainLogLagrangian
                    (chainEmb (Function.update s i x)))
            ↔ ∀ (i : Fin (M + 1)) (x : Fin 2),
                pathTransitionWeight chainWeight (Function.update s i x)
                  ≤ pathTransitionWeight chainWeight s)
      ∧ quadHamiltonStep 2 (fun _ => 0)
          (witnessEmb 0 (0 : Fin 2).castSucc,
            quadSegmentMomentum 2 (witnessEmb 0) 0)
          = (witnessEmb 0 (1 : Fin 2).castSucc,
              quadSegmentMomentum 2 (witnessEmb 0) 1) :=
  ⟨mixingFlow_symmetric, fun _ s => two_faces_mixing_chain s,
    witness_stationary_hamiltonStep⟩

end OPH.Variational

#print axioms OPH.Variational.quad_legendreTransform
#print axioms OPH.Variational.quad_legendre_involutive
#print axioms OPH.Variational.quad_solver_solves
#print axioms OPH.Variational.quad_hamiltonian_solver_solves
#print axioms OPH.Variational.quad_velocity_unique
#print axioms OPH.Variational.quadLagrangian_strictConvexOn
#print axioms OPH.Variational.quad_fenchel_young
#print axioms OPH.Variational.quad_fenchel_young_eq_iff
#print axioms OPH.Variational.quad_degenerate_every_velocity
#print axioms OPH.Variational.quad_degenerate_no_solution
#print axioms OPH.Variational.quad_degenerate_no_solver
#print axioms OPH.Variational.witnessL_eq_quadTwoPoint
#print axioms OPH.Variational.quad_discreteEulerLagrange_iff_hamiltonStep
#print axioms OPH.Variational.quad_stationary_hamiltonStep
#print axioms OPH.Variational.solver_section_retracts
#print axioms OPH.Variational.discreteEulerLagrange_iff_hamilton_strictConvex
#print axioms OPH.Variational.stationary_hamilton_strictConvex
#print axioms OPH.Variational.quadTwoPoint_realizes_strictConvex
#print axioms OPH.Variational.witness_stationary_hamiltonStep
#print axioms OPH.Variational.witness_hamilton_values
#print axioms OPH.Variational.pathTransitionWeight_eq_exp_neg
#print axioms OPH.Variational.logTransitionAction_le_iff
#print axioms OPH.Variational.logTransitionAction_site_min_iff
#print axioms OPH.Variational.localAction_chainEmb
#print axioms OPH.Variational.two_faces_mixing_chain
#print axioms OPH.Variational.two_faces_mixing_chain_global
#print axioms OPH.Variational.mixingFlow_symmetric
#print axioms OPH.Variational.mixingFlow_cross_value
#print axioms OPH.Variational.mixingChain_eq_flow_normalized
#print axioms OPH.Variational.noether_current_eq_legendre_momentum
#print axioms OPH.Variational.free_hamiltonStep_conserves_momentum
#print axioms OPH.Variational.free_hamiltonStep_conserves_hamiltonian
#print axioms OPH.Variational.witness_legendre_hamiltonian
#print axioms OPH.Variational.witness_noether_correspondence
#print axioms OPH.Variational.legendreBridge_receipt
