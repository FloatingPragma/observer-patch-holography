import Geometry.InternalEnergyInertia
import Geometry.InternalClockRestFrequency

set_option autoImplicit false

open scoped BigOperators

/-!
# Proper-time internal action: refinement invariance and the proper-time
principle select slope one
(issues 736, 739)

STATUS.  Selection theorem for the open slope of
`Geometry/InternalEnergyInertia.lean`.  That module proves the inertial
coefficient `m + E` under declared shape A and `m` under declared shape B
and names the slope family `m + λ E` with a per-step additive term as the
open selection.  This module declares two principles (a midpoint
refinement of the step index and a proper-time principle for the internal
process) and proves that, given them, the slope is one.  Every action
shape below is a declared enrichment; the selection is proved relative to
the declared principles and to nothing else.

WHAT IS PROVED.

(1) Refinement kinematics.  `refine x` is the declared midpoint refinement
of a step-indexed path (`refine x (2n) = x n`,
`refine x (2n+1) = (x n + x (n+1)) / 2`).  `properLength M x` is the sum
over the `M` window steps of the real square root of the committed Lorentz
square of the forward increment.  `properLength_refine` proves
`properLength (2M) (refine x) = properLength M x` for every path (each half
increment carries one quarter of the Lorentz square, the root gives one
half, two halves give one; on a step of negative Lorentz square both sides
read the clamped root zero).  `clockAction_refine` proves
`clockAction (2M) (refine x) = clockAction M x / 2` for the committed
quadratic functional, and the step count doubles by definition.

(2) Internal-action family.  `internalAction a b E M x =
E * (a * properLength M x + b * M)`.  `refinementInvariant_iff_b_zero`:
for `E ≠ 0`, invariance under the declared refinement on every timelike
window is equivalent to `b = 0`.  The additive per-step term of shape B
fails refinement invariance; the proper-length term passes it.

(3) Proper-time principle (declared).  `properTimeInternalAction E M x =
E * properLength M x`: the internal process accrues action `E` per unit
proper length.  `restPhase_step` and `restPhase_window_properLength`
connect it to the declared internal clock of
`Geometry/InternalClockRestFrequency.lean`: the rest phase increment over a
proper-time step equals the mass parameter times the step, and over the
sampled frame worldline the accumulated rest phase equals the mass
parameter times the proper length of the sampled path.  The declared
internal clock is the `b = 0`, `a = 1` member with `E` read as the mass
parameter.

(4) Length-form composite.  `lengthAction m E M F x =
(m + E) * properLength M x + lengthForceTerm M F x`, with the declared
potential pairing `lengthForceTerm M F x = forceTerm M F x / 2`.
`lengthAction_hasDerivAt`: on a timelike window and for every variation
`η`, the directional derivative of the length action along `η` at the path
is `lengthFirstVariation m E M F x η`, the explicit sum
`(m + E) * ∑ B(unitTangent x k, η (k+1) - η k) + lengthForceTerm M F η`.
`lengthFirstVariation_interior` telescopes it to the interior nodes;
`lengthStationary_iff_eom`: vanishing of the first variation on every
endpoint-fixed variation is equivalent to
`(m + E) • (unitTangent x (j+1) - unitTangent x j) = impulse F j` at every
interior node, by the nondegeneracy of the committed pairing.  The inertial
coefficient of the refinement-invariant form is `m + E`
(`length_inertial_coefficient`).  Invariance is proved under the declared
midpoint refinement only; no general reparametrization is defined.

(5) Einbein relation.  On a window with unit-Lorentz-square increments
`properLength = M` (`properLength_unit_increments`), the length first
variation is one half of the shape-A first variation
(`lengthFirstVariation_unit_eq_half_composite`), and the two equations of
motion coincide termwise (`length_eom_unit_iff_composite_eom`).  Nothing
is claimed off the unit gauge.

(6) Conclusion (`slope_selection`).  For the slope family
`slopeAction lam b m E M F x = m * properLength + internalAction lam b E +
lengthForceTerm`: refinement invariance of the internal part on timelike
windows is equivalent to `b = 0`; given `b = 0`, agreement of the internal
part with the proper-time principle on timelike windows is equivalent to
`lam = 1`; and at `lam = 1, b = 0` the family action is the length action,
whose equation of motion carries the coefficient `m + E`.  Every other
slope is an internal action that is not `E` per unit proper length
(`slope_ne_one_not_properTime`).

ROWS TOUCHED.  Coupled-action row (a declared length-form composite on the
kinematics island; the potential is step indexed with no worldline-field
coupling).  Source clock and duration row (the window index is the
declared step index; the refinement is a declared reindexing, and the
identification of proper length with a laboratory duration is the
laboratory clock and energy calibration import).  Physical spacetime
attachment row (`Herm2` is the declared Lorentz module).  Laboratory clock
and energy calibration import (the modular energy ledger is identified with
a physical internal energy only through it).  The mass parameter `m` is
declared.  This module discharges none of these rows.

NEGATIVES CITED.  The Legendre non-identifiability of
`Variational/RealizedHistoryLegendreNoGo.lean` at its scope: realized
source histories select no velocity curvature or Legendre map, so the
length form is a declared enrichment, selected among declared enrichments
by a declared invariance principle (`legendre_scope_cited`).  The abstract
rate non-identifiability (its repair-layer bridge fails per
`ObserverPatchHolography/RateBridgeObstruction.lean`) forbids nothing
about a source-hosted process and is neither used nor contradicted.

CONVENTIONS.  Signature `(+---)`; `Herm2 = ℝ × (Fin 3 → ℝ)`;
`lorentzQ v = v.1 ^ 2 - |v.2| ^ 2`; forward differences
`x (k + 1) - x k`; `Real.sqrt` clamps negative arguments to zero, so
`properLength` is the proper length on timelike-or-null windows and a
clamped sum otherwise; `unitTangent x k = (1 / sqrt (lorentzQ Δ_k)) • Δ_k`;
`impulse F j = F j - F (j + 1)` (negative forward difference, the declared
orientation of the potential, inherited); `lengthForceTerm = forceTerm / 2`
so that the length-form pairing has polarization constant one.

FALSIFIER.  The module is wrong if `properLength` or `clockAction` fail
the displayed refinement identities, if some `b ≠ 0` passes refinement
invariance at nonzero `E`, if the displayed derivative differs from the
true directional derivative on a timelike window, if stationarity differs
from the displayed equation of motion, or if the unit-gauge relation to
shape A misses a factor.

Axiom audit.  No project axiom, no native decision procedure; the guard
lines at the end show at most `propext`, `Classical.choice`, `Quot.sound`.
-/

namespace OPH.ProperTimeInternalAction

open OPH.C1Lorentz OPH.C2Soldering OPH.CommonWorldJointAction
open OPH.InternalEnergyInertia OPH.Variational

noncomputable section

/-! ## (1) Refinement kinematics -/

/-- The declared midpoint refinement of a step-indexed path: even nodes
copy the path, odd nodes take the midpoint of the adjacent nodes. -/
def refine (x : ℕ → Herm2) (n : ℕ) : Herm2 :=
  if n % 2 = 0 then x (n / 2)
  else (1 / 2 : ℝ) • (x (n / 2) + x (n / 2 + 1))

theorem refine_even (x : ℕ → Herm2) (n : ℕ) : refine x (2 * n) = x n := by
  unfold refine
  have h1 : (2 * n) % 2 = 0 := by omega
  have h2 : (2 * n) / 2 = n := by omega
  rw [if_pos h1, h2]

theorem refine_odd (x : ℕ → Herm2) (n : ℕ) :
    refine x (2 * n + 1) = (1 / 2 : ℝ) • (x n + x (n + 1)) := by
  unfold refine
  have h1 : (2 * n + 1) % 2 ≠ 0 := by omega
  have h2 : (2 * n + 1) / 2 = n := by omega
  rw [if_neg h1, h2]

/-- Each half increment of the refined path is one half of the parent
increment: the first half. -/
theorem refine_increment_first (x : ℕ → Herm2) (n : ℕ) :
    refine x (2 * n + 1) - refine x (2 * n) =
      (1 / 2 : ℝ) • (x (n + 1) - x n) := by
  rw [refine_odd, refine_even]
  module

/-- The second half. -/
theorem refine_increment_second (x : ℕ → Herm2) (n : ℕ) :
    refine x (2 * n + 1 + 1) - refine x (2 * n + 1) =
      (1 / 2 : ℝ) • (x (n + 1) - x n) := by
  rw [show 2 * n + 1 + 1 = 2 * (n + 1) by ring, refine_even, refine_odd]
  module

/-- A sum over a doubled window splits into even and odd terms. -/
theorem sum_range_two_mul (f : ℕ → ℝ) (M : ℕ) :
    ∑ n ∈ Finset.range (2 * M), f n =
      ∑ n ∈ Finset.range M, (f (2 * n) + f (2 * n + 1)) := by
  induction M with
  | zero => simp
  | succ M ih =>
    rw [show 2 * (M + 1) = 2 * M + 1 + 1 by ring, Finset.sum_range_succ,
      Finset.sum_range_succ, ih, Finset.sum_range_succ]
    ring

/-- The declared proper-length functional on the window `0, …, M`: the sum
of the real square roots of the committed Lorentz squares of the forward
increments.  On timelike-or-null increments this is the proper length; on a
negative Lorentz square the root clamps to zero. -/
def properLength (M : ℕ) (x : ℕ → Herm2) : ℝ :=
  ∑ k ∈ Finset.range M, Real.sqrt (lorentzQ (x (k + 1) - x k))

/-- The root of one quarter of a Lorentz square is one half of the root. -/
theorem sqrt_lorentzQ_half (v : Herm2) :
    Real.sqrt (lorentzQ ((1 / 2 : ℝ) • v)) = Real.sqrt (lorentzQ v) / 2 := by
  rw [lorentzQ_smul, Real.sqrt_mul (by positivity)]
  have : Real.sqrt ((1 / 2 : ℝ) ^ 2) = 1 / 2 := by
    rw [Real.sqrt_sq (by norm_num)]
  rw [this]
  ring

/-- **Refinement invariance of proper length.**  The proper length of the
refined path on the doubled window equals the proper length of the path. -/
theorem properLength_refine (M : ℕ) (x : ℕ → Herm2) :
    properLength (2 * M) (refine x) = properLength M x := by
  unfold properLength
  rw [sum_range_two_mul]
  refine Finset.sum_congr rfl fun n _ => ?_
  rw [refine_increment_first, refine_increment_second, sqrt_lorentzQ_half]
  ring

/-- **Refinement halves the committed quadratic functional.**  The step
count doubles and the clock action halves. -/
theorem clockAction_refine (M : ℕ) (x : ℕ → Herm2) :
    clockAction (2 * M) (refine x) = clockAction M x / 2 := by
  unfold clockAction
  rw [sum_range_two_mul, Finset.sum_div]
  refine Finset.sum_congr rfl fun n _ => ?_
  rw [refine_increment_first, refine_increment_second, lorentzQ_smul]
  ring


/-! ## (2) The internal-action family and refinement invariance -/

/-- A timelike window: every forward increment on the window has a
strictly positive Lorentz square. -/
def TimelikeWindow (M : ℕ) (x : ℕ → Herm2) : Prop :=
  ∀ k, k < M → 0 < lorentzQ (x (k + 1) - x k)

/-- The declared internal-action family: the ledger energy `E` multiplies
a proper-length term with weight `a` and a per-step additive term with
weight `b`.  Shape B of `Geometry/InternalEnergyInertia.lean` is the
member `a = 0, b = 1`; the proper-time form is `a = 1, b = 0`. -/
def internalAction (a b E : ℝ) (M : ℕ) (x : ℕ → Herm2) : ℝ :=
  E * (a * properLength M x + b * M)

/-- Refinement invariance of an internal action on timelike windows: the
action on the doubled window of the refined path equals the action on the
window. -/
def RefinementInvariant (a b E : ℝ) : Prop :=
  ∀ M : ℕ, ∀ x : ℕ → Herm2, TimelikeWindow M x →
    internalAction a b E (2 * M) (refine x) = internalAction a b E M x

/-- The rest path: the scalar coordinate counts the step, the spatial
coordinates vanish.  Every increment has Lorentz square one. -/
def restPath (n : ℕ) : Herm2 := ((n : ℝ), (0 : Spatial))

theorem restPath_increment (n : ℕ) :
    restPath (n + 1) - restPath n = ((1 : ℝ), (0 : Spatial)) := by
  unfold restPath
  ext <;> simp

theorem restPath_lorentzQ (n : ℕ) :
    lorentzQ (restPath (n + 1) - restPath n) = 1 := by
  rw [restPath_increment]
  simp [lorentzQ, spatialNormSq]

theorem restPath_timelike (M : ℕ) : TimelikeWindow M restPath := by
  intro k _
  rw [restPath_lorentzQ]
  exact one_pos

/-- The proper length of the rest path on the window is the step count. -/
theorem restPath_properLength (M : ℕ) : properLength M restPath = M := by
  unfold properLength
  rw [Finset.sum_congr rfl fun k _ => by rw [restPath_lorentzQ, Real.sqrt_one]]
  simp

/-- **Refinement invariance selects the proper-length form.**  For a
nonzero ledger energy, the internal action is invariant under the declared
refinement on every timelike window exactly when the per-step additive
weight vanishes.  The forward direction uses the rest path on the window
of one step: the invariance equation `E (a + 2 b) = E (a + b)` gives
`b = 0`.  The reverse direction is `properLength_refine`. -/
theorem refinementInvariant_iff_b_zero (a b E : ℝ) (hE : E ≠ 0) :
    RefinementInvariant a b E ↔ b = 0 := by
  constructor
  · intro h
    have h1 := h 1 restPath (restPath_timelike 1)
    unfold internalAction at h1
    rw [properLength_refine, restPath_properLength] at h1
    push_cast at h1
    have : E * b = 0 := by linarith
    rcases mul_eq_zero.mp this with h0 | h0
    · exact absurd h0 hE
    · exact h0
  · intro hb M x _
    unfold internalAction
    rw [hb, properLength_refine]
    ring

/-- Shape B's additive per-step term (`a = 0, b = 1`) fails refinement
invariance at every nonzero ledger energy. -/
theorem additive_term_not_refinementInvariant (E : ℝ) (hE : E ≠ 0) :
    ¬ RefinementInvariant 0 1 E := by
  rw [refinementInvariant_iff_b_zero 0 1 E hE]
  exact one_ne_zero

/-- The proper-length term (`a = 1, b = 0`) is refinement invariant. -/
theorem properLength_term_refinementInvariant (E : ℝ) :
    RefinementInvariant 1 0 E := by
  intro M x _
  unfold internalAction
  rw [properLength_refine]
  ring

/-! ## (3) The proper-time principle (declared) and the internal clock -/

/-- The declared proper-time principle: the internal process accrues
action `E` per unit proper length, so the internal action is the ledger
energy times the proper length of the window.  This is the member
`a = 1, b = 0` of the declared family. -/
def properTimeInternalAction (E : ℝ) (M : ℕ) (x : ℕ → Herm2) : ℝ :=
  E * properLength M x

theorem properTimeInternalAction_eq (E : ℝ) (M : ℕ) (x : ℕ → Herm2) :
    properTimeInternalAction E M x = internalAction 1 0 E M x := by
  unfold properTimeInternalAction internalAction
  ring

/-- The rest phase of the declared internal clock of
`Geometry/InternalClockRestFrequency.lean`, as a function of the
proper-time coordinate along the frame worldline. -/
def restPhase (mass : ℝ) (frame : FrameHyperboloid) (τ : ℝ) : ℝ :=
  planeWavePhase (fourMomentum mass frame) (frameWorldline frame τ)

/-- The rest phase increment over a proper-time step equals the mass
parameter times the step: the declared internal clock accrues phase at
rate `mass` per unit proper time. -/
theorem restPhase_step (mass : ℝ) (frame : FrameHyperboloid) (τ δ : ℝ) :
    restPhase mass frame (τ + δ) - restPhase mass frame τ = mass * δ := by
  unfold restPhase
  rw [planeWavePhase_frameWorldline, planeWavePhase_frameWorldline]
  ring

/-- The frame worldline sampled at proper-time steps of declared size
`δ`. -/
def sampledWorldline (frame : FrameHyperboloid) (δ : ℝ) (n : ℕ) : Herm2 :=
  frameWorldline frame (n * δ)

theorem sampledWorldline_increment (frame : FrameHyperboloid) (δ : ℝ)
    (n : ℕ) :
    sampledWorldline frame δ (n + 1) - sampledWorldline frame δ n =
      δ • (frame : Herm2) := by
  unfold sampledWorldline frameWorldline
  rw [← sub_smul]
  push_cast
  ring_nf

theorem sampledWorldline_lorentzQ (frame : FrameHyperboloid) (δ : ℝ)
    (n : ℕ) :
    lorentzQ (sampledWorldline frame δ (n + 1) - sampledWorldline frame δ n)
      = δ ^ 2 := by
  rw [sampledWorldline_increment, lorentzQ_smul, frame.2.1, mul_one]

/-- The proper length of the sampled frame worldline is the step count
times the step size, for a nonnegative step size. -/
theorem sampledWorldline_properLength (frame : FrameHyperboloid) {δ : ℝ}
    (hδ : 0 ≤ δ) (M : ℕ) :
    properLength M (sampledWorldline frame δ) = M * δ := by
  unfold properLength
  rw [Finset.sum_congr rfl fun k _ => by
    rw [sampledWorldline_lorentzQ, Real.sqrt_sq hδ]]
  simp

/-- **The declared internal clock is the proper-time form.**  Over the
sampled frame worldline the accumulated rest phase equals the mass
parameter times the proper length of the sampled path: the phase is the
internal action `properTimeInternalAction mass` with the ledger energy
read as the declared mass parameter. -/
theorem restPhase_window_properLength (mass : ℝ) (frame : FrameHyperboloid)
    {δ : ℝ} (hδ : 0 ≤ δ) (M : ℕ) :
    restPhase mass frame (M * δ) - restPhase mass frame 0 =
      properTimeInternalAction mass M (sampledWorldline frame δ) := by
  unfold properTimeInternalAction
  rw [sampledWorldline_properLength frame hδ M]
  have h := restPhase_step mass frame 0 (M * δ)
  rw [zero_add] at h
  rw [h]


/-! ## (4) The length-form composite and its first variation -/

/-- Exact expansion of the Lorentz square along a line. -/
theorem lorentzQ_add_smul (u v : Herm2) (t : ℝ) :
    lorentzQ (u + t • v) =
      lorentzQ u + t * (2 * lorentzB u v) + t ^ 2 * lorentzQ v := by
  rw [lorentzQ_add, lorentzB_smul_right, lorentzQ_smul]
  ring

/-- The quadratic polynomial along the line has derivative `2 B u v` at
the origin. -/
theorem hasDerivAt_lorentzQ_line (u v : Herm2) :
    HasDerivAt (fun t : ℝ => lorentzQ (u + t • v)) (2 * lorentzB u v) 0 := by
  have hfun : (fun t : ℝ => lorentzQ (u + t • v)) =
      fun t : ℝ => lorentzQ u + t * (2 * lorentzB u v) + t ^ 2 * lorentzQ v := by
    funext t
    exact lorentzQ_add_smul u v t
  rw [hfun]
  have h1 : HasDerivAt (fun t : ℝ => t * (2 * lorentzB u v)) (2 * lorentzB u v) 0 := by
    simpa using (hasDerivAt_id (0 : ℝ)).mul_const (2 * lorentzB u v)
  have h2 : HasDerivAt (fun t : ℝ => t ^ 2 * lorentzQ v)
      ((2 : ℕ) * (0 : ℝ) ^ (2 - 1) * lorentzQ v) 0 :=
    (hasDerivAt_pow 2 (0 : ℝ)).mul_const (lorentzQ v)
  have h := (h1.const_add (lorentzQ u)).add h2
  convert h using 1
  simp

/-- **Derivative of the root of the Lorentz square along a line.**  At a
vector of strictly positive Lorentz square the derivative at the origin is
the pairing divided by the root. -/
theorem hasDerivAt_sqrt_lorentzQ_line (u v : Herm2) (hu : 0 < lorentzQ u) :
    HasDerivAt (fun t : ℝ => Real.sqrt (lorentzQ (u + t • v)))
      (lorentzB u v / Real.sqrt (lorentzQ u)) 0 := by
  have hne : lorentzQ (u + (0 : ℝ) • v) ≠ 0 := by
    rw [zero_smul, add_zero]; exact ne_of_gt hu
  have h := (hasDerivAt_lorentzQ_line u v).sqrt hne
  convert h using 1
  rw [zero_smul, add_zero]
  have hs : Real.sqrt (lorentzQ u) ≠ 0 := by
    exact ne_of_gt (Real.sqrt_pos.mpr hu)
  field_simp

/-- The declared unit tangent at step `k`: the forward increment divided by
the root of its Lorentz square.  On a step of nonpositive Lorentz square the
root clamps to zero and the division reads zero. -/
def unitTangent (x : ℕ → Herm2) (k : ℕ) : Herm2 :=
  (1 / Real.sqrt (lorentzQ (x (k + 1) - x k))) • (x (k + 1) - x k)

theorem lorentzB_unitTangent (x : ℕ → Herm2) (k : ℕ) (w : Herm2) :
    lorentzB (unitTangent x k) w =
      lorentzB (x (k + 1) - x k) w / Real.sqrt (lorentzQ (x (k + 1) - x k)) := by
  unfold unitTangent
  rw [lorentzB_smul_left]
  ring

/-- The declared potential pairing of the length form: the inherited
`forceTerm` normalized by one half, so that the pairing carries the
polarization constant one of the linearized length form. -/
def lengthForceTerm (M : ℕ) (F x : ℕ → Herm2) : ℝ := forceTerm M F x / 2

theorem forceTerm_smul (M : ℕ) (F η : ℕ → Herm2) (t : ℝ) :
    forceTerm M F (t • η) = t * forceTerm M F η := by
  unfold forceTerm
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun k _ => ?_
  have hstep : (t • η) (k + 1) - (t • η) k = t • (η (k + 1) - η k) := by
    simp only [Pi.smul_apply]
    rw [smul_sub]
  rw [hstep, lorentzB_smul_right]
  ring

theorem lengthForceTerm_line (M : ℕ) (F x η : ℕ → Herm2) (t : ℝ) :
    lengthForceTerm M F (x + t • η) =
      lengthForceTerm M F x + t * lengthForceTerm M F η := by
  unfold lengthForceTerm
  rw [forceTerm_add, forceTerm_smul]
  ring

/-- The declared length-form composite: the rest parameter plus the ledger
energy multiplies the proper length, plus the declared potential pairing. -/
def lengthAction (m E : ℝ) (M : ℕ) (F x : ℕ → Herm2) : ℝ :=
  (m + E) * properLength M x + lengthForceTerm M F x

/-- The explicit first variation of the length form: the unit tangents
paired with the variation increments, weighted by `m + E`, plus the
potential pairing of the variation. -/
def lengthFirstVariation (m E : ℝ) (M : ℕ) (F x η : ℕ → Herm2) : ℝ :=
  (m + E) * (∑ k ∈ Finset.range M, lorentzB (unitTangent x k) (η (k + 1) - η k))
    + lengthForceTerm M F η

theorem line_increment (x η : ℕ → Herm2) (t : ℝ) (k : ℕ) :
    (x + t • η) (k + 1) - (x + t • η) k =
      (x (k + 1) - x k) + t • (η (k + 1) - η k) := by
  simp only [Pi.add_apply, Pi.smul_apply]
  rw [smul_sub]
  abel

/-- **Directional derivative of the length action.**  On a timelike window
the derivative at the path along any variation is the explicit first
variation. -/
theorem lengthAction_hasDerivAt (m E : ℝ) (M : ℕ) (F x η : ℕ → Herm2)
    (hx : TimelikeWindow M x) :
    HasDerivAt (fun t : ℝ => lengthAction m E M F (x + t • η))
      (lengthFirstVariation m E M F x η) 0 := by
  have hfun : (fun t : ℝ => lengthAction m E M F (x + t • η)) =
      fun t : ℝ => (m + E) * (∑ k ∈ Finset.range M,
        Real.sqrt (lorentzQ ((x (k + 1) - x k) + t • (η (k + 1) - η k))))
        + (lengthForceTerm M F x + t * lengthForceTerm M F η) := by
    funext t
    unfold lengthAction properLength
    rw [lengthForceTerm_line]
    congr 2
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [line_increment]
  rw [hfun]
  have hsum : HasDerivAt (fun t : ℝ => ∑ k ∈ Finset.range M,
      Real.sqrt (lorentzQ ((x (k + 1) - x k) + t • (η (k + 1) - η k))))
      (∑ k ∈ Finset.range M, lorentzB (x (k + 1) - x k) (η (k + 1) - η k) /
        Real.sqrt (lorentzQ (x (k + 1) - x k))) 0 := by
    exact HasDerivAt.fun_sum (u := Finset.range M)
      (A := fun k t => Real.sqrt (lorentzQ ((x (k + 1) - x k) + t • (η (k + 1) - η k))))
      (A' := fun k => lorentzB (x (k + 1) - x k) (η (k + 1) - η k) /
        Real.sqrt (lorentzQ (x (k + 1) - x k)))
      (fun k hk => hasDerivAt_sqrt_lorentzQ_line _ _ (hx k (Finset.mem_range.mp hk)))
  have hlin : HasDerivAt (fun t : ℝ => lengthForceTerm M F x + t * lengthForceTerm M F η)
      (lengthForceTerm M F η) 0 := by
    simpa using ((hasDerivAt_id (0 : ℝ)).mul_const (lengthForceTerm M F η)).const_add
      (lengthForceTerm M F x)
  have h := (hsum.const_mul (m + E)).add hlin
  convert h using 1
  unfold lengthFirstVariation
  congr 2
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [lorentzB_unitTangent]


/-- A pairing sum against variation increments is one half of the inherited
`forceTerm` at the pairing sequence. -/
theorem sum_pairing_eq_half_forceTerm (M : ℕ) (G η : ℕ → Herm2) :
    (∑ k ∈ Finset.range M, lorentzB (G k) (η (k + 1) - η k)) =
      forceTerm M G η / 2 := by
  unfold forceTerm
  rw [Finset.sum_div]
  refine Finset.sum_congr rfl fun k _ => ?_
  ring

/-- The residual of the length-form equation of motion at the interior
node `j + 1`. -/
def lengthResidual (m E : ℝ) (F x : ℕ → Herm2) (j : ℕ) : Herm2 :=
  (m + E) • (unitTangent x j - unitTangent x (j + 1)) + impulse F j

/-- **Interior form of the length first variation.**  For an endpoint-fixed
variation on the window `0, …, N + 1`, the first variation is the sum over
the interior nodes of the residual paired with the variation. -/
theorem lengthFirstVariation_interior (m E : ℝ) (N : ℕ) (F x η : ℕ → Herm2)
    (h0 : η 0 = 0) (hN : η (N + 1) = 0) :
    lengthFirstVariation m E (N + 1) F x η =
      ∑ k ∈ Finset.range N, lorentzB (lengthResidual m E F x k) (η (k + 1)) := by
  unfold lengthFirstVariation lengthForceTerm
  rw [sum_pairing_eq_half_forceTerm, forceTerm_interior N (unitTangent x) η h0 hN,
    forceTerm_interior N F η h0 hN, Finset.sum_div, Finset.sum_div, Finset.mul_sum,
    ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  unfold lengthResidual impulse
  rw [lorentzB_add_left, lorentzB_smul_left]
  ring

/-- Stationarity of the length form on the window `0, …, M`: the first
variation vanishes on every variation fixed at the two endpoints.  By
`lengthAction_hasDerivAt` this is the vanishing of the directional
derivative on a timelike window. -/
def LengthStationary (m E : ℝ) (M : ℕ) (F x : ℕ → Herm2) : Prop :=
  ∀ η : ℕ → Herm2, η 0 = 0 → η M = 0 → lengthFirstVariation m E M F x η = 0

theorem lengthResidual_eq_zero_iff (m E : ℝ) (F x : ℕ → Herm2) (j : ℕ) :
    lengthResidual m E F x j = 0 ↔
      (m + E) • (unitTangent x (j + 1) - unitTangent x j) = impulse F j := by
  unfold lengthResidual
  rw [show unitTangent x j - unitTangent x (j + 1) =
      -(unitTangent x (j + 1) - unitTangent x j) by abel, smul_neg,
    neg_add_eq_sub, sub_eq_zero, eq_comm]

/-- Length stationarity gives the equation of motion at every interior
node: the single-node variation isolates the residual, and nondegeneracy of
the committed pairing kills it. -/
theorem lengthStationary_eom {m E : ℝ} {N : ℕ} {F x : ℕ → Herm2}
    (hs : LengthStationary m E (N + 1) F x) {j : ℕ} (hj : j < N) :
    (m + E) • (unitTangent x (j + 1) - unitTangent x j) = impulse F j := by
  have key : ∀ w : Herm2, lorentzB (lengthResidual m E F x j) w = 0 := by
    intro w
    let η : ℕ → Herm2 := fun n => if n = j + 1 then w else 0
    have h0 : η 0 = 0 := if_neg (by omega)
    have hN1 : η (N + 1) = 0 := if_neg (by omega)
    have hfv := hs η h0 hN1
    rw [lengthFirstVariation_interior m E N F x η h0 hN1] at hfv
    have hterm : ∀ k ∈ Finset.range N,
        lorentzB (lengthResidual m E F x k) (η (k + 1)) =
          if k = j then lorentzB (lengthResidual m E F x j) w else 0 := by
      intro k _
      by_cases hkj : k = j
      · subst hkj
        rw [if_pos rfl]
        have hv : η (k + 1) = w := if_pos rfl
        rw [hv]
      · rw [if_neg hkj]
        have hv : η (k + 1) = 0 := if_neg (by omega)
        rw [hv, lorentzB_zero_right]
    rw [Finset.sum_congr rfl hterm, Finset.sum_ite_eq',
      if_pos (Finset.mem_range.mpr hj)] at hfv
    exact hfv
  exact (lengthResidual_eq_zero_iff m E F x j).mp (lorentzB_nondegenerate key)

/-- The equation of motion at every interior node gives length
stationarity: the interior first variation vanishes termwise. -/
theorem lengthStationary_of_eom {m E : ℝ} {N : ℕ} {F x : ℕ → Herm2}
    (heom : ∀ j, j < N →
      (m + E) • (unitTangent x (j + 1) - unitTangent x j) = impulse F j) :
    LengthStationary m E (N + 1) F x := by
  intro η h0 hN
  rw [lengthFirstVariation_interior m E N F x η h0 hN]
  refine Finset.sum_eq_zero fun k hk => ?_
  rw [(lengthResidual_eq_zero_iff m E F x k).mpr (heom k (Finset.mem_range.mp hk)),
    lorentzB_zero_left]

/-- **Headline (given the declared length form).**  On a positive window,
length stationarity under endpoint-fixed variations is equivalent to
`(m + E) • (unitTangent x (j + 1) - unitTangent x j) = impulse F j` at every
interior node. -/
theorem lengthStationary_iff_eom (m E : ℝ) (N : ℕ) (F x : ℕ → Herm2) :
    LengthStationary m E (N + 1) F x ↔
      ∀ j, j < N →
        (m + E) • (unitTangent x (j + 1) - unitTangent x j) = impulse F j :=
  ⟨fun hs _ hj => lengthStationary_eom hs hj, lengthStationary_of_eom⟩

/-- **The inertial coefficient of the refinement-invariant form is
`m + E`.**  (Refinement means the declared midpoint refinement.)  On a
timelike window the derivative of the length action along every
endpoint-fixed variation is the first variation, and its vanishing on
all such variations is the displayed equation of motion with coefficient
`inertialCoefficient m E = m + E` on the unit-tangent difference; the
coefficient exceeds its zero-energy value by exactly `E`. -/
theorem length_inertial_coefficient (m E : ℝ) (N : ℕ) (F x : ℕ → Herm2)
    (hx : TimelikeWindow (N + 1) x) :
    (∀ η : ℕ → Herm2, HasDerivAt (fun t : ℝ => lengthAction m E (N + 1) F (x + t • η))
        (lengthFirstVariation m E (N + 1) F x η) 0) ∧
    (LengthStationary m E (N + 1) F x ↔
      ∀ j, j < N →
        inertialCoefficient m E • (unitTangent x (j + 1) - unitTangent x j) =
          impulse F j) ∧
    inertialCoefficient m E - inertialCoefficient m 0 = E := by
  refine ⟨fun η => lengthAction_hasDerivAt m E (N + 1) F x η hx,
    lengthStationary_iff_eom m E N F x, ?_⟩
  unfold inertialCoefficient
  ring


/-! ## (5) Einbein relation: the unit gauge -/

/-- On a window with unit-Lorentz-square increments the proper length is
the step count. -/
theorem properLength_unit_increments {M : ℕ} {x : ℕ → Herm2}
    (hu : ∀ k, k < M → lorentzQ (x (k + 1) - x k) = 1) :
    properLength M x = (M : ℝ) := by
  unfold properLength
  rw [Finset.sum_congr rfl fun k hk => by
    rw [hu k (Finset.mem_range.mp hk), Real.sqrt_one]]
  simp

/-- On a unit-Lorentz-square step the unit tangent is the increment. -/
theorem unitTangent_unit {x : ℕ → Herm2} {k : ℕ}
    (hk : lorentzQ (x (k + 1) - x k) = 1) :
    unitTangent x k = x (k + 1) - x k := by
  unfold unitTangent
  rw [hk, Real.sqrt_one, div_one, one_smul]

/-- **Einbein relation at the level of first variations.**  On a window
with unit-Lorentz-square increments the length first variation is one half
of the shape-A first variation of `Geometry/InternalEnergyInertia.lean`,
for every variation.  The factor one half is the ratio of the polarization
constants of the root and of the square at a unit argument. -/
theorem lengthFirstVariation_unit_eq_half_composite (m E : ℝ) (M : ℕ)
    (F x η : ℕ → Herm2) (hu : ∀ k, k < M → lorentzQ (x (k + 1) - x k) = 1) :
    lengthFirstVariation m E M F x η = compositeFirstVariation m E M F x η / 2 := by
  unfold lengthFirstVariation compositeFirstVariation lengthForceTerm
    clockFirstVariation
  rw [Finset.sum_congr rfl fun k hk => by
    rw [unitTangent_unit (hu k (Finset.mem_range.mp hk))]]
  rw [add_div, mul_div_assoc, Finset.sum_div]
  congr 2
  refine Finset.sum_congr rfl fun k _ => ?_
  ring

/-- **Einbein relation at the level of equations of motion.**  On a window
with unit-Lorentz-square increments the length-form equation of motion at
every interior node is the shape-A equation of motion
`(m + E) • secondDifference x j = impulse F j` at every interior node. -/
theorem length_eom_unit_iff_composite_eom (m E : ℝ) (N : ℕ) (F x : ℕ → Herm2)
    (hu : ∀ k, k < N + 1 → lorentzQ (x (k + 1) - x k) = 1) :
    (∀ j, j < N →
        (m + E) • (unitTangent x (j + 1) - unitTangent x j) = impulse F j) ↔
      ∀ j, j < N → (m + E) • secondDifference x j = impulse F j := by
  have hT : ∀ j, j < N →
      unitTangent x (j + 1) - unitTangent x j = secondDifference x j := by
    intro j hj
    rw [unitTangent_unit (hu (j + 1) (by omega)), unitTangent_unit (hu j (by omega))]
    rfl
  constructor
  · intro h j hj
    rw [← hT j hj]
    exact h j hj
  · intro h j hj
    rw [hT j hj]
    exact h j hj

/-- On the unit gauge the length action and shape A agree up to the
declared potential normalization: the proper length equals the clock
action and the potential pairing is one half of the inherited one. -/
theorem lengthAction_unit (m E : ℝ) (M : ℕ) (F x : ℕ → Herm2)
    (hu : ∀ k, k < M → lorentzQ (x (k + 1) - x k) = 1) :
    lengthAction m E M F x = (m + E) * clockAction M x + forceTerm M F x / 2 := by
  unfold lengthAction lengthForceTerm
  rw [properLength_unit_increments hu, clockAction_unit_increments hu]

/-! ## (6) Conclusion: the slope selection -/

/-- The declared slope family: the rest parameter multiplies the proper
length, the ledger energy enters through the declared internal-action
family with slope `lam` and per-step weight `b`, and the declared
potential pairing is added.  Shape A of
`Geometry/InternalEnergyInertia.lean` sits at `lam = 1, b = 0` on the unit
gauge and shape B at `lam = 0, b = 1`. -/
def slopeAction (lam b m E : ℝ) (M : ℕ) (F x : ℕ → Herm2) : ℝ :=
  m * properLength M x + internalAction lam b E M x + lengthForceTerm M F x

/-- Agreement of an internal action with the declared proper-time
principle on every timelike window. -/
def ProperTimeAgreement (lam b E : ℝ) : Prop :=
  ∀ M : ℕ, ∀ x : ℕ → Herm2, TimelikeWindow M x →
    internalAction lam b E M x = properTimeInternalAction E M x

/-- At `b = 0` and nonzero ledger energy, the proper-time principle holds
on every timelike window exactly when the slope is one: the rest path on
one step reads `E * lam = E`. -/
theorem properTimeAgreement_iff_slope_one (lam E : ℝ) (hE : E ≠ 0) :
    ProperTimeAgreement lam 0 E ↔ lam = 1 := by
  constructor
  · intro h
    have h1 := h 1 restPath (restPath_timelike 1)
    unfold internalAction properTimeInternalAction at h1
    rw [restPath_properLength] at h1
    push_cast at h1
    have : E * (lam - 1) = 0 := by linarith
    rcases mul_eq_zero.mp this with h0 | h0
    · exact absurd h0 hE
    · linarith
  · intro hl M x _
    unfold internalAction properTimeInternalAction
    rw [hl]
    ring

/-- Every slope other than one is an internal action that is not `E` per
unit proper length: a timelike window witnesses the disagreement. -/
theorem slope_ne_one_not_properTime (lam E : ℝ) (hE : E ≠ 0) (hl : lam ≠ 1) :
    ∃ M : ℕ, ∃ x : ℕ → Herm2, TimelikeWindow M x ∧
      internalAction lam 0 E M x ≠ properTimeInternalAction E M x := by
  by_contra hcon
  push Not at hcon
  exact hl ((properTimeAgreement_iff_slope_one lam E hE).mp
    fun M x hx => hcon M x hx)

theorem slopeAction_one_zero (m E : ℝ) (M : ℕ) (F x : ℕ → Herm2) :
    slopeAction 1 0 m E M F x = lengthAction m E M F x := by
  unfold slopeAction lengthAction internalAction
  ring

/-- At zero per-step weight the slope family is the length action with
ledger energy rescaled by the slope; so on a timelike window the member
with slope `lam` has inertial coefficient `m + lam * E` by
`length_inertial_coefficient m (lam * E)`. -/
theorem slopeAction_zero_eq_lengthAction (lam m E : ℝ) (M : ℕ) (F x : ℕ → Herm2) :
    slopeAction lam 0 m E M F x = lengthAction m (lam * E) M F x := by
  unfold slopeAction lengthAction internalAction
  ring

/-- The general-slope inertial coefficient: the member of the family with
slope `lam` and zero per-step weight is stationary (first variation
vanishing on every endpoint-fixed variation, an algebraic statement valid
for every path with the clamped unit tangent)
exactly when `(m + lam * E) • (unitTangent x (j + 1) - unitTangent x j) =
impulse F j` at every interior node. -/
theorem slopeAction_inertial_coefficient (lam m E : ℝ) (N : ℕ) (F x : ℕ → Herm2) :
    LengthStationary m (lam * E) (N + 1) F x ↔
      ∀ j, j < N →
        (m + lam * E) • (unitTangent x (j + 1) - unitTangent x j) = impulse F j :=
  lengthStationary_iff_eom m (lam * E) N F x

/-- **Conclusion: refinement invariance and the proper-time principle
select slope one.**  For a nonzero ledger energy: (i) refinement invariance
of the internal part on timelike windows holds exactly when the per-step
weight vanishes; (ii) given that, agreement with the declared proper-time
principle on timelike windows holds exactly when the slope is one, which
is the uniqueness of the declared member within the family (the principle
is stated as `E * properLength`, the `lam = 1` member, so (ii) says no
other member agrees with it on all timelike windows);
(iii) at slope one and zero per-step weight the family action is the
length action; (iv) on a timelike window the length action is stationary
exactly when `(m + E) • (unitTangent x (j + 1) - unitTangent x j) =
impulse F j` at every interior node, the inertial coefficient being
`m + E`.  The refinement, the proper-time principle, the potential
pairing, and the mass parameter are declared; the selection given them is
proved. -/
theorem slope_selection (m E : ℝ) (hE : E ≠ 0) :
    (∀ lam b : ℝ, RefinementInvariant lam b E ↔ b = 0) ∧
    (∀ lam : ℝ, ProperTimeAgreement lam 0 E ↔ lam = 1) ∧
    (∀ M : ℕ, ∀ F x : ℕ → Herm2,
      slopeAction 1 0 m E M F x = lengthAction m E M F x) ∧
    (∀ N : ℕ, ∀ F x : ℕ → Herm2, TimelikeWindow (N + 1) x →
      (LengthStationary m E (N + 1) F x ↔
        ∀ j, j < N →
          inertialCoefficient m E • (unitTangent x (j + 1) - unitTangent x j) =
            impulse F j)) :=
  ⟨fun lam b => refinementInvariant_iff_b_zero lam b E hE,
    fun lam => properTimeAgreement_iff_slope_one lam E hE,
    fun M F x => slopeAction_one_zero m E M F x,
    fun N F x hx => (length_inertial_coefficient m E N F x hx).2.1⟩

/-! ## Negatives cited and rows touched -/

/-- The Legendre non-identifiability at its scope, re-cited from
`Geometry/InternalEnergyInertia.lean`: realized source histories select no
velocity curvature or Legendre map, so the length form and every member of
the slope family are declared enrichments; the selection above is among
declared enrichments by a declared invariance principle. -/
theorem legendre_scope_cited :
    (¬ ∃ vel : ℝ → ℝ → ℝ, SolvesMomentum chainLogLagrangian vel) ∧
      chainCurvedLagrangian 1 ≠ chainCurvedLagrangian 2 :=
  legendre_nonidentifiability_cited

/-- The rows this module touches, as register labels: the coupled action
(declared length form), the source clock and duration row (declared step
index and refinement), and the physical spacetime attachment row (declared
Lorentz module).  The laboratory clock and energy calibration import is
named in the header and carries no `OpenRow` label, so it is absent from
this list.  A label is not a discharge. -/
def touchedRows : List OpenRow :=
  [OpenRow.coupledAction, OpenRow.sourceClock, OpenRow.spacetimeAttachment]

/-- The rows this module discharges: none. -/
def dischargedRows : List OpenRow := []

theorem dischargedRows_empty : dischargedRows = [] := rfl

end

end OPH.ProperTimeInternalAction

/- Axiom audit: expected at most `propext`, `Classical.choice`, `Quot.sound`
per line; no native decision procedure. -/

#print axioms OPH.ProperTimeInternalAction.properLength_refine
#print axioms OPH.ProperTimeInternalAction.clockAction_refine
#print axioms OPH.ProperTimeInternalAction.refinementInvariant_iff_b_zero
#print axioms OPH.ProperTimeInternalAction.additive_term_not_refinementInvariant
#print axioms OPH.ProperTimeInternalAction.properLength_term_refinementInvariant
#print axioms OPH.ProperTimeInternalAction.restPhase_step
#print axioms OPH.ProperTimeInternalAction.restPhase_window_properLength
#print axioms OPH.ProperTimeInternalAction.hasDerivAt_sqrt_lorentzQ_line
#print axioms OPH.ProperTimeInternalAction.lengthAction_hasDerivAt
#print axioms OPH.ProperTimeInternalAction.lengthFirstVariation_interior
#print axioms OPH.ProperTimeInternalAction.lengthStationary_iff_eom
#print axioms OPH.ProperTimeInternalAction.length_inertial_coefficient
#print axioms OPH.ProperTimeInternalAction.properLength_unit_increments
#print axioms OPH.ProperTimeInternalAction.lengthFirstVariation_unit_eq_half_composite
#print axioms OPH.ProperTimeInternalAction.length_eom_unit_iff_composite_eom
#print axioms OPH.ProperTimeInternalAction.lengthAction_unit
#print axioms OPH.ProperTimeInternalAction.properTimeAgreement_iff_slope_one
#print axioms OPH.ProperTimeInternalAction.slope_ne_one_not_properTime
#print axioms OPH.ProperTimeInternalAction.slope_selection
#print axioms OPH.ProperTimeInternalAction.slopeAction_zero_eq_lengthAction
#print axioms OPH.ProperTimeInternalAction.slopeAction_inertial_coefficient
#print axioms OPH.ProperTimeInternalAction.legendre_scope_cited
#print axioms OPH.ProperTimeInternalAction.dischargedRows_empty
