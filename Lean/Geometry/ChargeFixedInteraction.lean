import Geometry.CarrierDynamicsCompatibility
import Geometry.CommonWorldJointAction
import Geometry.MassShellKinematics
import Geometry.FreeEvolutionPersistence

set_option autoImplicit false

open scoped BigOperators

namespace OPH.ChargeFixedInteraction

open OPH.ScreenCarrierMapCandidate
open OPH.SeamCurrentCarrierQuotient
open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.CarrierDynamicsCompatibility
open OPH.CommonWorldJointAction
open OPH.C1Lorentz (Spatial Herm2 spatialDot lorentzB lorentzQ)

/-!
# Charge-fixed interaction normalization and the coupled worldline equation

`Geometry/CarrierDynamicsCompatibility.lean` declares the interaction
candidate `interactionCandidate κ h A φ x N`, pairing the scaled seam
electric field with the worldline increment through the embedded seam
directions, and proves its relative normalization `κ` non-forced.  This
file varies the coupled action, the committed scaled window action plus
the candidate, in both sectors, and shows that one declared identification
fixes `κ`.  Every statement concerns a candidate coupled action on the
finite thirty-seam complex.

WHAT IS PROVED.

* Exact seam geometry.  The incidence boundary of the embedded seam
  direction field at every port is the constant `6 - 2φ` times the
  candidate ray of that port (`seamBoundaryVecZ_eq`, kernel computation
  over `ℤ[φ]`; real form `realBoundary_seamPaired`).  The candidate ray of
  a seam endpoint pairs with the seam direction to exactly `2` at the
  larger endpoint and `-2` at the smaller one (`ray_dot_seamVector_right`,
  `ray_dot_seamVector_left`).  The seam-crossing constant `12 - 4φ` is
  twice the boundary constant and is positive (`seamCrossingConst_pos`).

* (a) Field sector.  The coupled action equals the committed window action
  at augmented sources plus an endpoint term depending on the seam
  potential only through its two window endpoint values
  (`coupledAction_reduction`).  The augmented sources are the declared seam
  current plus the induced current `J_int κ m = (κ / h²) • (j (m+1) - j m)`
  and the declared port load plus the induced load
  `ρ_int κ n = -(κ / h) • ∂ (j n)`, where `j n` is the worldline seam
  current, the worldline increment paired with each embedded seam
  direction (`worldlineSeamCurrent`).  The field-sector first variation is
  therefore the committed Euler-Lagrange form at the augmented sources
  (`coupledAction_field_variation_interior`): stationarity in the seam
  potential is the committed scaled Ampere update at the augmented current
  (`coupled_stationary_A_iff_ampere`, `ampereResidual_augmented`), and
  stationarity in the port potential is the committed Gauss constraint at
  the augmented load (`coupled_stationary_phi_iff_gauss`); both together
  are `coupled_field_equations`.  The induced sources satisfy the committed
  scaled continuity equation identically (`induced_continuity`), and the
  induced load is neutral at every step (`inducedLoad_total`).

* (b) Charge fixing, conditional.  The induced load at port `p` across
  step `n` is `-(κ / h)` times the crossing factor, the boundary constant
  times the candidate ray of `p` paired with the increment
  (`inducedLoad_eq_crossingFactor`).  The port-charge identification is
  the declared clause "induced load at `p` equals the declared charge
  `q`" (`PortChargeIdentification`).  At nonzero step and nonzero crossing
  factor the identification holds exactly when
  `κ = -(q h) / crossingFactor` (`identification_iff_kappa`), so the
  normalization is unique under the identification
  (`kappa_unique_of_identification`).  On the seam-crossing worldline the
  crossing factor at the larger endpoint is `12 - 4φ`, computed exactly
  from the `ℤ[φ]` seam table and nonzero, and the identification reads
  `κ = -(q h) / (12 - 4φ)` (`seam_crossing_charge_fixed`); the smaller
  endpoint then carries `-q` (`seam_crossing_left_charge`).  Without the
  identification the normalization stays free: the augmented Gauss
  constraint is satisfiable at every `κ` by a declared load, and the
  committed non-forcing holds unchanged
  (`normalization_free_without_identification`).

* (c) Worldline sector.  With the declared weight `m` on the committed
  worldline action, the total action moves under an endpoint-fixed
  worldline variation by the coupled residual paired against the interior
  variation values plus the weighted clock action of the variation
  (`totalAction_worldline_expansion`).  Stationarity is the vanishing of
  the residual at every interior node (`worldline_stationary_iff`), and
  the residual vanishes exactly when `2 m` times the forward second
  difference of the path equals the purely spatial vector `κ` times the
  forward difference of the embedded seam electric field
  (`coupled_worldline_equation`, components in
  `coupled_worldline_equation_components`, momentum form in
  `coupled_momentum_law`).  At zero coupling and unit weight this is the
  committed equal-increment clause
  (`worldline_stationary_zero_coupling_iff_uniform`).

SIGN AND FRAME CONVENTIONS, RE-DERIVED.  Committed conventions: the scaled
electric field is `E n = -(A (n+1) - A n) / h - d (φ n)`; the port boundary
`∂` is the transpose of the port coboundary `d`, so `∂ c` at a port is the
net signed inflow, the continuum analogue of minus the divergence; the
committed Gauss constraint is `∂ E = ρ`; the committed Ampere update is
`E (m+1) - E m = h • (Cᵀ B (m+1) - J m)`; the committed Lagrangian carries
`+ h ⟨J, A⟩ + h ⟨ρ, φ⟩`.  Independent continuum check of the induced
sources: with a coupling density `κ E · v`, the seam-potential variation
of `E` is `-δȦ`, and one integration by parts turns `-κ ∫ δȦ · v` into
`+κ ∫ δA · v̇`, an additional source current `κ v̇` entering with the same
sign as `J`; the port-potential variation of `E` is `-∇ δφ`, and one
integration by parts turns `-κ ∫ ∇δφ · v` into `+κ ∫ δφ ∇·v`, an
additional load `κ ∇·v = -κ ∂ v` in the committed inflow convention.
With `v = j / h` and `v̇ = (j (m+1) - j m) / h²` this is `J_int` and
`ρ_int` exactly.  The induced load is the bound charge of a polarization
`κ v` in the committed sign convention, which is why it is neutral.
Because `∂` after `d` is the positive graph Laplacian (degree `5` at
every port), the committed constraint `∂ E = ρ` corresponds to
`div E = -ρ` in continuum terms, so a committed positive port load carries
the sign of a continuum negative charge density; the sign of the fixed
normalization below is stated relative to that convention.
Worldline sign: the committed quadratic form is `(+---)`, so the spatial
part of the worldline density is `-m |Δx|²`; the continuum
Euler-Lagrange equation of `-m |ẋ|² + κ E · ẋ` is `-2 m ẍ + κ Ė = 0`, and
the discrete law `2 m Δ²x = κ ΔE` matches it, with the scalar coordinate
free of any force.

WHAT IS DECLARED AND WHAT IS DERIVED.  Declared: the coupling shape (the
electric field paired with the increment through the candidate
embedding), the candidate embedding itself, the weight `m` on the
worldline action, the worldline action shape, the step `h`, the window,
and the port-charge identification.  Derived from these: the augmented
sources, the field equations, the continuity identity, the neutrality of
the induced load, the unique fixed value of `κ` under the identification,
the coupled residual, and the equation of motion.

WHAT IS NOT SUPPLIED.  The force term is `κ` times the step difference of
the embedded electric field, a polarization-type term produced by the
declared coupling shape; it differs from the Lorentz force
`q (E + v × B)`, no magnetic term appears, and no theorem identifies the
fixed `κ` with a monopole charge: the induced load is neutral, and the
identification names the value at one port of a neutral load.  A
minimal-coupling candidate pairing the seam potential with the increment
is a live alternative shape and is outside this file.  No physical
spacetime attachment, no calibration, no unit, no readout, no observer
map.  The weight `m` is a declared parameter of the worldline action
shape on the coupled-action row, the same declared mass parameter as in
`Geometry/MassShellKinematics.lean` (`discreteMomentum_frame`); it is
attached to no energy readout and touches neither the physical spacetime
attachment row nor the source clock and duration row.  The fixed value
depends on the step through `κ / h` (`chargeFixedKappa_step_scaling`), so
no step is selected and the source clock and duration row is untouched.
The coupled-action row of the committed register stays open: the charge
identification is declared, and `RealizedHistoryLegendreNoGo` shows
realized source histories do not select the velocity curvature or the
Legendre map, so the worldline action shape and its weight are declared
enrichments.  `RateNonidentifiability` is an abstract-transition-system
result whose bridge to the repair layer fails (`RateBridgeObstruction`),
so nothing here cites it for or against a repair-layer clock.  The
declared charge is a free datum (`declared_charge_unconstrained`).

FALSIFIER.  The module fails if a port row of the boundary table
disagrees with the constant `6 - 2φ`, if an endpoint pairing differs from
`±2`, if the reduction identity misses a term, if the induced sources
violate continuity, if the fixed value of `κ` differs from
`-(q h) / (12 - 4φ)` on the seam crossing, or if the worldline residual
disagrees with the exact expansion.
-/

/-! ## Exact port boundary of the embedded seam directions -/

/-- The exact `ℤ[φ]` port boundary of the embedded seam direction field:
at port `p`, the signed sum over seams of the seam direction vector, with
the committed incidence sign (`+` at the larger endpoint, `-` at the
smaller endpoint). -/
def seamBoundaryVecZ (p : Fin 12) (k : Fin 3) : Zphi :=
  ∑ e : Fin 30,
    ((if p = seamRight e then seamVectorZ e k else 0) -
      (if p = seamLeft e then seamVectorZ e k else 0))

/-- The committed geometric constant `6 - 2φ` of the candidate embedding. -/
def boundaryConstZ : Zphi := ((6 : ℤ), (-2 : ℤ))

/-- Exact table identity: the port boundary of the embedded seam direction
field is the constant `6 - 2φ` times the candidate ray of the port, at
every port and coordinate. -/
theorem seamBoundaryVecZ_eq :
    ∀ (p : Fin 12) (k : Fin 3),
      seamBoundaryVecZ p k = zmul boundaryConstZ (candidateRayZ p k) := by
  decide

theorem boundaryConstZ_ne_zero : boundaryConstZ ≠ 0 := by decide

/-- The self-pairing of the candidate ray with the embedded direction of a
seam it terminates: exactly `2` at the larger endpoint, at every seam. -/
theorem ray_dot_seamVector_right :
    ∀ e : Fin 30,
      dotZ (candidateRayZ (seamRight e)) (seamVectorZ e) = ((2 : ℤ), (0 : ℤ)) := by
  decide

/-- The pairing at the smaller endpoint: exactly `-2`, at every seam. -/
theorem ray_dot_seamVector_left :
    ∀ e : Fin 30,
      dotZ (candidateRayZ (seamLeft e)) (seamVectorZ e) = ((-2 : ℤ), (0 : ℤ)) := by
  decide

/-- The seam-crossing constant `12 - 4φ`: twice the boundary constant. -/
def seamCrossingConstZ : Zphi := ((12 : ℤ), (-4 : ℤ))

theorem seamCrossingConstZ_eq : seamCrossingConstZ = zmul boundaryConstZ (2, 0) := by
  decide

theorem seamCrossingConstZ_ne_zero : seamCrossingConstZ ≠ 0 := by decide

noncomputable section

/-! ## Real layer of the boundary identity -/

theorem evalPhi_zero : evalPhi 0 = 0 := by
  simp [evalPhi]

theorem evalPhi_sub (x y : Zphi) : evalPhi (x - y) = evalPhi x - evalPhi y := by
  obtain ⟨a, b⟩ := x
  obtain ⟨c, d⟩ := y
  show ((a - c : ℤ) : ℝ) + ((b - d : ℤ) : ℝ) * Real.goldenRatio =
    ((a : ℝ) + (b : ℝ) * Real.goldenRatio) -
      ((c : ℝ) + (d : ℝ) * Real.goldenRatio)
  push_cast
  ring

/-- Evaluation of `ℤ[φ]` as an additive homomorphism. -/
def evalPhiHom : Zphi →+ ℝ where
  toFun := evalPhi
  map_zero' := evalPhi_zero
  map_add' := evalPhi_add

theorem evalPhiHom_apply (x : Zphi) : evalPhiHom x = evalPhi x := rfl

/-- The real value of the committed geometric constant: `6 - 2φ`. -/
def boundaryConst : ℝ := evalPhi boundaryConstZ

theorem boundaryConst_eq : boundaryConst = 6 - 2 * Real.goldenRatio := by
  show ((6 : ℤ) : ℝ) + ((-2 : ℤ) : ℝ) * Real.goldenRatio = 6 - 2 * Real.goldenRatio
  push_cast
  ring

theorem boundaryConst_pos : 0 < boundaryConst := by
  rw [boundaryConst_eq]
  linarith [Real.goldenRatio_lt_two]

theorem boundaryConst_ne_zero : boundaryConst ≠ 0 := boundaryConst_pos.ne'

/-- The real seam-crossing constant `12 - 4φ`. -/
def seamCrossingConst : ℝ := evalPhi seamCrossingConstZ

theorem seamCrossingConst_eq : seamCrossingConst = 2 * boundaryConst := by
  rw [seamCrossingConst, seamCrossingConstZ_eq, evalPhi_zmul, boundaryConst]
  rw [mul_comm]
  congr 1
  show ((2 : ℤ) : ℝ) + ((0 : ℤ) : ℝ) * Real.goldenRatio = 2
  push_cast
  ring

theorem seamCrossingConst_pos : 0 < seamCrossingConst := by
  rw [seamCrossingConst_eq]
  linarith [boundaryConst_pos]

theorem seamCrossingConst_ne_zero : seamCrossingConst ≠ 0 :=
  seamCrossingConst_pos.ne'

/-- Coordinate form of the boundary identity over `ℝ`: the incidence
boundary of the embedded seam direction field at port `p` is the
boundary constant times the candidate ray of `p`. -/
theorem seamVector_boundary_coord (p : Fin 12) (i : Fin 3) :
    (∑ e : Fin 30,
      ((if p = seamRight e then seamVector e i else 0) -
        (if p = seamLeft e then seamVector e i else 0))) =
      boundaryConst * candidateRay p i := by
  have hterm : ∀ e : Fin 30,
      ((if p = seamRight e then seamVector e i else 0) -
        (if p = seamLeft e then seamVector e i else 0)) =
      evalPhiHom ((if p = seamRight e then seamVectorZ e i else 0) -
        (if p = seamLeft e then seamVectorZ e i else 0)) := by
    intro e
    rw [evalPhiHom_apply, evalPhi_sub, apply_ite evalPhi, apply_ite evalPhi,
      evalPhi_zero]
    rfl
  rw [Finset.sum_congr rfl fun e _ ↦ hterm e, ← map_sum]
  change evalPhi (seamBoundaryVecZ p i) = _
  rw [seamBoundaryVecZ_eq, evalPhi_zmul]
  rfl


/-! ## The worldline seam current and its induced sources -/

/-- The worldline seam current across step `n`: the spatial increment of
the worldline paired with each embedded seam direction.  This is the seam
field the interaction candidate pairs against the electric field. -/
def worldlineSeamCurrent (x : ℕ → Herm2) (n : ℕ) : Fin 30 → ℝ :=
  fun e ↦ spatialDot (seamVector e) ((x (n + 1)).2 - (x n).2)

/-- The interaction candidate is the window sum of the seam pairing of the
scaled electric field against the worldline seam current. -/
theorem interactionCandidate_eq_inner (κ h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (M : ℕ) :
    interactionCandidate κ h A φ x M =
      κ * ∑ n ∈ Finset.range M,
        realSeamInner (electricFieldScaled h A φ n) (worldlineSeamCurrent x n) :=
  rfl

theorem worldlineSeamCurrent_add (x η : ℕ → Herm2) (n : ℕ) :
    worldlineSeamCurrent (x + η) n =
      worldlineSeamCurrent x n + worldlineSeamCurrent η n := by
  funext e
  simp only [worldlineSeamCurrent, Pi.add_apply, Prod.snd_add, spatialDot,
    Pi.sub_apply]
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  ring

/-- The port boundary of a seam field of the form "embedded seam direction
paired with a fixed spatial vector": the boundary constant times the
candidate ray of the port paired with that vector. -/
theorem realBoundary_seamPaired (v : Spatial) (p : Fin 12) :
    realBoundary (fun e ↦ spatialDot (seamVector e) v) p =
      boundaryConst * spatialDot (candidateRay p) v := by
  rw [realBoundary_apply]
  have hterm : ∀ e : Fin 30,
      ((if p = seamRight e then spatialDot (seamVector e) v else 0) -
        (if p = seamLeft e then spatialDot (seamVector e) v else 0)) =
      ∑ i : Fin 3, ((if p = seamRight e then seamVector e i else 0) -
        (if p = seamLeft e then seamVector e i else 0)) * v i := by
    intro e
    split_ifs <;> simp [spatialDot]
  rw [Finset.sum_congr rfl fun e _ ↦ hterm e, Finset.sum_comm]
  simp only [← Finset.sum_mul, seamVector_boundary_coord]
  unfold spatialDot
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  ring

/-- The port boundary of the worldline seam current. -/
theorem realBoundary_worldlineSeamCurrent (x : ℕ → Herm2) (n : ℕ) (p : Fin 12) :
    realBoundary (worldlineSeamCurrent x n) p =
      boundaryConst * spatialDot (candidateRay p) ((x (n + 1)).2 - (x n).2) :=
  realBoundary_seamPaired _ p

/-- The induced seam current of the interaction candidate at step `m`:
`κ / h²` times the forward step difference of the worldline seam
current.  It enters the scaled Ampere update added to the declared seam
current. -/
def inducedCurrent (κ h : ℝ) (x : ℕ → Herm2) (m : ℕ) : Fin 30 → ℝ :=
  (κ / h ^ 2) • (worldlineSeamCurrent x (m + 1) - worldlineSeamCurrent x m)

/-- The induced port load of the interaction candidate at step `n`: minus
`κ / h` times the port boundary of the worldline seam current.  It enters
the Gauss constraint added to the declared port load. -/
def inducedLoad (κ h : ℝ) (x : ℕ → Herm2) (n : ℕ) : Fin 12 → ℝ :=
  -((κ / h) • realBoundary (worldlineSeamCurrent x n))

/-- Port value of the induced load: minus `κ / h` times the boundary
constant times the candidate ray of the port paired with the worldline
increment. -/
theorem inducedLoad_apply (κ h : ℝ) (x : ℕ → Herm2) (n : ℕ) (p : Fin 12) :
    inducedLoad κ h x n p =
      -(κ / h) * (boundaryConst *
        spatialDot (candidateRay p) ((x (n + 1)).2 - (x n).2)) := by
  simp only [inducedLoad, Pi.neg_apply, Pi.smul_apply, smul_eq_mul,
    realBoundary_worldlineSeamCurrent]
  ring

/-- The induced load is neutral: its total over the twelve ports vanishes
at every step and every normalization. -/
theorem inducedLoad_total (κ h : ℝ) (x : ℕ → Herm2) (n : ℕ) :
    (∑ p : Fin 12, inducedLoad κ h x n p) = 0 := by
  simp only [inducedLoad, Pi.neg_apply, Pi.smul_apply, smul_eq_mul,
    Finset.sum_neg_distrib, ← Finset.mul_sum, realBoundary_total, mul_zero,
    neg_zero]

/-- The induced sources satisfy the committed scaled continuity equation
identically. -/
theorem induced_continuity (κ h : ℝ) (hh : h ≠ 0) (x : ℕ → Herm2) (n : ℕ) :
    inducedLoad κ h x (n + 1) - inducedLoad κ h x n +
      h • realBoundary (inducedCurrent κ h x n) = 0 := by
  unfold inducedLoad inducedCurrent
  rw [map_smul, map_sub]
  funext p
  simp only [Pi.add_apply, Pi.sub_apply, Pi.neg_apply, Pi.smul_apply,
    smul_eq_mul, Pi.zero_apply]
  field_simp
  ring

theorem portInner_neg_left (x y : Fin 12 → ℝ) :
    realPortInner (-x) y = -realPortInner x y := by
  unfold realPortInner
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun p _ ↦ ?_
  simp

theorem portInner_smul_left (c : ℝ) (x y : Fin 12 → ℝ) :
    realPortInner (c • x) y = c * realPortInner x y := by
  unfold realPortInner
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun p _ ↦ ?_
  simp [mul_assoc]

theorem portInner_add_left (x y z : Fin 12 → ℝ) :
    realPortInner (x + y) z = realPortInner x z + realPortInner y z := by
  rw [realPortInner_comm, portInner_add_right, realPortInner_comm z x,
    realPortInner_comm z y]


/-! ## The coupled action and its reduction to augmented sources -/

/-- The coupled field action: the committed scaled window action plus the
interaction candidate over the same `N + 1` steps. -/
def coupledAction (κ h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (x : ℕ → Herm2) : ℝ :=
  windowAction h N A φ ρ J + interactionCandidate κ h A φ x (N + 1)

/-- The endpoint term of the reduction: `κ / h` times the seam pairing of
the worldline seam current against the seam potential, at the first node
minus at the last node.  It depends on the seam potential only through
its two window endpoint values and does not depend on the port
potential. -/
def endpointTerm (κ h : ℝ) (x : ℕ → Herm2) (N : ℕ) (A : ℕ → Fin 30 → ℝ) : ℝ :=
  (κ / h) * (realSeamInner (worldlineSeamCurrent x 0) (A 0) -
    realSeamInner (worldlineSeamCurrent x (N + 1)) (A (N + 1)))

/-- Source augmentation of the committed window action: adding a seam
current and a port load adds the window sum of their pairings against the
forward seam potential and the port potential. -/
theorem windowAction_augment (h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (ρ ρ' : ℕ → Fin 12 → ℝ) (J J' : ℕ → Fin 30 → ℝ) :
    windowAction h N A φ (ρ + ρ') (J + J') = windowAction h N A φ ρ J +
      ∑ n ∈ Finset.range (N + 1),
        (h * realSeamInner (J' n) (A (n + 1)) + h * realPortInner (ρ' n) (φ n)) := by
  unfold windowAction
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun n _ ↦ ?_
  unfold stepLagrangian localSourcedAction
  simp only [Pi.add_apply]
  rw [seamInner_add_left, portInner_add_left]
  ring

/-- Per-step summation by parts: the interaction density at step `n` is
the augmented-source density plus a telescoping difference. -/
theorem interaction_step_parts (κ h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (n : ℕ) :
    κ * realSeamInner (electricFieldScaled h A φ n) (worldlineSeamCurrent x n) =
      (h * realSeamInner (inducedCurrent κ h x n) (A (n + 1)) +
        h * realPortInner (inducedLoad κ h x n) (φ n)) +
      ((κ / h) * realSeamInner (worldlineSeamCurrent x n) (A n) -
        (κ / h) * realSeamInner (worldlineSeamCurrent x (n + 1)) (A (n + 1))) := by
  rw [realSeamInner_comm]
  have hk := kinetic_pairing h hh (worldlineSeamCurrent x n) A φ n
  unfold inducedCurrent inducedLoad
  rw [seamInner_smul_left, realSeamInner_sub_left, portInner_neg_left,
    portInner_smul_left]
  have hX : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  simp only [div_eq_mul_inv, pow_two, mul_inv]
  linear_combination (κ * h⁻¹) * hk -
    (κ * realSeamInner (worldlineSeamCurrent x n) (electricFieldScaled h A φ n) +
      κ * h⁻¹ * (realSeamInner (worldlineSeamCurrent x (n + 1)) (A (n + 1)) -
        realSeamInner (worldlineSeamCurrent x n) (A (n + 1)))) * hX

/-- **Reduction identity.**  The coupled action equals the committed window
action at the augmented sources (declared port load plus induced load,
declared seam current plus induced current) plus the endpoint term. -/
theorem coupledAction_reduction (κ h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (x : ℕ → Herm2) :
    coupledAction κ h N A φ ρ J x =
      windowAction h N A φ (ρ + inducedLoad κ h x) (J + inducedCurrent κ h x) +
        endpointTerm κ h x N A := by
  unfold coupledAction
  rw [windowAction_augment, interactionCandidate_eq_inner, Finset.mul_sum,
    Finset.sum_congr rfl fun n _ ↦ interaction_step_parts κ h hh A φ x n,
    Finset.sum_add_distrib, Finset.sum_range_sub']
  unfold endpointTerm
  ring

/-- The endpoint term is invariant under seam-potential variations vanishing
at both window endpoints. -/
theorem endpointTerm_fixed (κ h : ℝ) (x : ℕ → Herm2) (N : ℕ)
    (A a : ℕ → Fin 30 → ℝ) (ha0 : a 0 = 0) (haN : a (N + 1) = 0) :
    endpointTerm κ h x N (A + a) = endpointTerm κ h x N A := by
  unfold endpointTerm
  simp only [Pi.add_apply, ha0, haN, add_zero]


/-! ## (a) Field-sector variation of the coupled action -/

/-- **Exact field-sector expansion.**  For seam-potential variations
vanishing at both window endpoints and arbitrary port-potential
variations, the coupled action moves by the committed first variation at
the augmented sources plus the committed quadratic remainder. -/
theorem coupledAction_field_expansion (κ h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (x : ℕ → Herm2) (a : ℕ → Fin 30 → ℝ)
    (f : ℕ → Fin 12 → ℝ) (ha0 : a 0 = 0) (haN : a (N + 1) = 0) :
    coupledAction κ h N (A + a) (φ + f) ρ J x =
      coupledAction κ h N A φ ρ J x +
        firstVariation h N A φ (ρ + inducedLoad κ h x) (J + inducedCurrent κ h x)
          a f +
        quadraticRemainder h N a f := by
  rw [coupledAction_reduction κ h hh, coupledAction_reduction κ h hh,
    endpointTerm_fixed κ h x N A a ha0 haN, windowAction_expansion]
  ring

/-- The field-sector first variation in Euler-Lagrange form: the scaled
Ampere residual at the augmented current paired against the interior
seam variation, minus `h` times the Gauss residual at the augmented load
paired against the port variation. -/
theorem coupledAction_field_variation_interior (κ h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (x : ℕ → Herm2) (a : ℕ → Fin 30 → ℝ)
    (f : ℕ → Fin 12 → ℝ) (ha0 : a 0 = 0) (haN : a (N + 1) = 0) :
    coupledAction κ h N (A + a) (φ + f) ρ J x =
      coupledAction κ h N A φ ρ J x +
        ((∑ m ∈ Finset.range N,
          realSeamInner (ampereResidual h A φ (J + inducedCurrent κ h x) m)
            (a (m + 1))) -
          h * ∑ n ∈ Finset.range (N + 1),
            realPortInner (realBoundary (electricFieldScaled h A φ n) -
              (ρ n + inducedLoad κ h x n)) (f n)) +
        quadraticRemainder h N a f := by
  rw [coupledAction_field_expansion κ h hh N A φ ρ J x a f ha0 haN,
    firstVariation_interior h hh N A φ _ _ a f ha0 haN]
  rfl

/-- The Ampere residual at the augmented current is the committed residual
plus `κ / h` times the forward difference of the worldline seam current:
the source term of the interaction candidate. -/
theorem ampereResidual_augmented (κ h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (x : ℕ → Herm2) (m : ℕ) :
    ampereResidual h A φ (J + inducedCurrent κ h x) m =
      ampereResidual h A φ J m +
        (κ / h) • (worldlineSeamCurrent x (m + 1) - worldlineSeamCurrent x m) := by
  unfold ampereResidual inducedCurrent
  simp only [Pi.add_apply]
  funext e
  simp only [Pi.sub_apply, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  field_simp
  ring

/-- **Stationarity in the seam potential is the scaled Ampere update at the
augmented current.**  The coupled action is stationary under endpoint-fixed
seam-potential variations exactly when the committed scaled Ampere residual
vanishes at every interior step with the declared seam current replaced by
the declared current plus the induced current. -/
theorem coupled_stationary_A_iff_ampere (κ h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (x : ℕ → Herm2) :
    (∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
      coupledAction κ h N (A + a) φ ρ J x =
        coupledAction κ h N A φ ρ J x + quadraticRemainder h N a 0) ↔
    (∀ m, m < N → ampereResidual h A φ (J + inducedCurrent κ h x) m = 0) := by
  rw [← action_stationary_A_iff_ampere h hh N A φ (ρ + inducedLoad κ h x)
    (J + inducedCurrent κ h x)]
  refine forall_congr' fun a ↦ forall_congr' fun ha0 ↦ forall_congr' fun haN ↦ ?_
  rw [coupledAction_reduction κ h hh, coupledAction_reduction κ h hh,
    endpointTerm_fixed κ h x N A a ha0 haN]
  constructor <;> intro hs <;> linarith

/-- **Stationarity in the port potential is the Gauss constraint at the
augmented load.**  The coupled action is stationary under port-potential
variations exactly when the committed Gauss constraint holds at every
window step with the declared load replaced by the declared load plus the
induced load. -/
theorem coupled_stationary_phi_iff_gauss (κ h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (x : ℕ → Herm2) :
    (∀ f : ℕ → Fin 12 → ℝ,
      coupledAction κ h N A (φ + f) ρ J x =
        coupledAction κ h N A φ ρ J x + quadraticRemainder h N 0 f) ↔
    (∀ n, n < N + 1 →
      realBoundary (electricFieldScaled h A φ n) = ρ n + inducedLoad κ h x n) := by
  have key := action_stationary_phi_iff_gauss h hh N A φ (ρ + inducedLoad κ h x)
    (J + inducedCurrent κ h x)
  simp only [Pi.add_apply] at key
  rw [← key]
  refine forall_congr' fun f ↦ ?_
  rw [coupledAction_reduction κ h hh, coupledAction_reduction κ h hh]
  constructor <;> intro hs <;> linarith

/-- The coupled field equations: the scaled Ampere update at the augmented
current together with the Gauss constraint at the augmented load, in the
committed sign convention `E (m+1) - E m = h • (Cᵀ B (m+1) - J_total m)`
and `∂ E n = ρ_total n`. -/
theorem coupled_field_equations (κ h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (x : ℕ → Herm2) :
    ((∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
        coupledAction κ h N (A + a) φ ρ J x =
          coupledAction κ h N A φ ρ J x + quadraticRemainder h N a 0) ∧
      (∀ f : ℕ → Fin 12 → ℝ,
        coupledAction κ h N A (φ + f) ρ J x =
          coupledAction κ h N A φ ρ J x + quadraticRemainder h N 0 f)) ↔
    ((∀ m, m < N →
        electricFieldScaled h A φ (m + 1) - electricFieldScaled h A φ m =
          h • (localMaxwellOperator (A (m + 1)) - (J m + inducedCurrent κ h x m))) ∧
      (∀ n, n < N + 1 →
        realBoundary (electricFieldScaled h A φ n) = ρ n + inducedLoad κ h x n)) := by
  rw [coupled_stationary_A_iff_ampere κ h hh, coupled_stationary_phi_iff_gauss κ h hh]
  refine and_congr (forall_congr' fun m ↦ imp_congr Iff.rfl ?_) Iff.rfl
  unfold ampereResidual
  rw [sub_eq_zero]
  rfl


/-! ## (b) Charge fixing of the relative normalization -/

/-- The crossing factor of a port and a spatial increment: the boundary
constant times the candidate ray of the port paired with the increment.
Exact geometric data of the candidate embedding; no unit attached. -/
def crossingFactor (p : Fin 12) (v : Spatial) : ℝ :=
  boundaryConst * spatialDot (candidateRay p) v

theorem inducedLoad_eq_crossingFactor (κ h : ℝ) (x : ℕ → Herm2) (n : ℕ)
    (p : Fin 12) :
    inducedLoad κ h x n p =
      -(κ / h) * crossingFactor p ((x (n + 1)).2 - (x n).2) :=
  inducedLoad_apply κ h x n p

/-- The declared port-charge identification: the induced load at port `p`
across step `n` equals the declared charge `q`.  This clause is a
declaration, not a theorem; the port load of the committed Gauss
constraint is identified with a charge only here.  The identification
concerns the induced load only: it is step-local and proportional to the
increment (a resting worldline induces no load, `inducedLoad_rest`), so
the identified value is a polarization-type port load carried for the
crossing step; the declared load at the port is taken as zero here or
added separately in the augmented constraint. -/
def PortChargeIdentification (κ h q : ℝ) (x : ℕ → Herm2) (n : ℕ)
    (p : Fin 12) : Prop :=
  inducedLoad κ h x n p = q

/-- The normalization value selected by the identification: minus the
declared charge times the step, divided by the crossing factor. -/
def chargeFixedKappa (q h : ℝ) (p : Fin 12) (v : Spatial) : ℝ :=
  -(q * h) / crossingFactor p v

/-- **Charge fixing.**  At nonzero step and nonzero crossing factor, the
port-charge identification holds exactly when the relative normalization
equals the charge-fixed value: the identification fixes `κ` uniquely. -/
theorem identification_iff_kappa (κ h q : ℝ) (hh : h ≠ 0) (x : ℕ → Herm2)
    (n : ℕ) (p : Fin 12)
    (hc : crossingFactor p ((x (n + 1)).2 - (x n).2) ≠ 0) :
    PortChargeIdentification κ h q x n p ↔
      κ = chargeFixedKappa q h p ((x (n + 1)).2 - (x n).2) := by
  unfold PortChargeIdentification chargeFixedKappa
  rw [inducedLoad_eq_crossingFactor]
  constructor
  · intro hq
    rw [← hq]
    field_simp
  · intro hk
    rw [hk]
    field_simp

/-- Uniqueness: two normalizations satisfying the same identification at a
nonzero crossing factor coincide. -/
theorem kappa_unique_of_identification (κ₁ κ₂ h q : ℝ) (hh : h ≠ 0)
    (x : ℕ → Herm2) (n : ℕ) (p : Fin 12)
    (hc : crossingFactor p ((x (n + 1)).2 - (x n).2) ≠ 0)
    (h₁ : PortChargeIdentification κ₁ h q x n p)
    (h₂ : PortChargeIdentification κ₂ h q x n p) : κ₁ = κ₂ := by
  rw [identification_iff_kappa κ₁ h q hh x n p hc] at h₁
  rw [identification_iff_kappa κ₂ h q hh x n p hc] at h₂
  rw [h₁, h₂]

/-- The seam-crossing worldline: the spatial part moves by the embedded
direction of seam `e` across the first step and rests afterwards. -/
def seamCrossingWorldline (e : Fin 30) : ℕ → Herm2 :=
  fun n ↦ ((0 : ℝ), if n = 0 then (0 : Spatial) else seamVector e)

theorem seamCrossing_increment (e : Fin 30) :
    (seamCrossingWorldline e (0 + 1)).2 - (seamCrossingWorldline e 0).2 =
      seamVector e := by
  show (if (1 : ℕ) = 0 then (0 : Spatial) else seamVector e) -
    (if (0 : ℕ) = 0 then (0 : Spatial) else seamVector e) = seamVector e
  norm_num

/-- Real form of the endpoint pairings: `2` at the larger endpoint. -/
theorem ray_seamVector_right (e : Fin 30) :
    spatialDot (candidateRay (seamRight e)) (seamVector e) = 2 := by
  unfold candidateRay seamVector
  rw [← evalPhi_dotZ, ray_dot_seamVector_right]
  show ((2 : ℤ) : ℝ) + ((0 : ℤ) : ℝ) * Real.goldenRatio = 2
  push_cast
  ring

/-- Real form of the endpoint pairings: `-2` at the smaller endpoint. -/
theorem ray_seamVector_left (e : Fin 30) :
    spatialDot (candidateRay (seamLeft e)) (seamVector e) = -2 := by
  unfold candidateRay seamVector
  rw [← evalPhi_dotZ, ray_dot_seamVector_left]
  show ((-2 : ℤ) : ℝ) + ((0 : ℤ) : ℝ) * Real.goldenRatio = -2
  push_cast
  ring

theorem crossingFactor_right (e : Fin 30) :
    crossingFactor (seamRight e) (seamVector e) = seamCrossingConst := by
  rw [crossingFactor, ray_seamVector_right, seamCrossingConst_eq]
  ring

theorem crossingFactor_left (e : Fin 30) :
    crossingFactor (seamLeft e) (seamVector e) = -seamCrossingConst := by
  rw [crossingFactor, ray_seamVector_left, seamCrossingConst_eq]
  ring

/-- **Seam-crossing charge fixing.**  For the worldline crossing seam `e`
across the first step, the identification of the induced load at the
larger endpoint with a declared charge `q` holds exactly when
`κ = -(q h) / (12 - 4φ)`, with the constant computed exactly from the
`ℤ[φ]` seam data and proved nonzero. -/
theorem seam_crossing_charge_fixed (κ h q : ℝ) (hh : h ≠ 0) (e : Fin 30) :
    PortChargeIdentification κ h q (seamCrossingWorldline e) 0 (seamRight e) ↔
      κ = -(q * h) / seamCrossingConst := by
  rw [identification_iff_kappa κ h q hh _ 0 _
    (by rw [seamCrossing_increment, crossingFactor_right];
        exact seamCrossingConst_ne_zero)]
  unfold chargeFixedKappa
  rw [seamCrossing_increment, crossingFactor_right]

/-- The induced load of the seam crossing is a neutral pair at the seam
endpoints: the value at the smaller endpoint is the negative of the value
at the larger endpoint. -/
theorem seam_crossing_dipole (κ h : ℝ) (e : Fin 30) :
    inducedLoad κ h (seamCrossingWorldline e) 0 (seamLeft e) =
      -inducedLoad κ h (seamCrossingWorldline e) 0 (seamRight e) := by
  rw [inducedLoad_eq_crossingFactor, inducedLoad_eq_crossingFactor,
    seamCrossing_increment, crossingFactor_left, crossingFactor_right]
  ring

/-- Under the identification the smaller endpoint carries `-q`. -/
theorem seam_crossing_left_charge (κ h q : ℝ) (e : Fin 30)
    (hid : PortChargeIdentification κ h q (seamCrossingWorldline e) 0
      (seamRight e)) :
    inducedLoad κ h (seamCrossingWorldline e) 0 (seamLeft e) = -q := by
  rw [seam_crossing_dipole]
  unfold PortChargeIdentification at hid
  rw [hid]

/-- The identification fixes the ratio `κ / h` only: scaling the step
scales the fixed normalization.  No step is selected. -/
theorem chargeFixedKappa_step_scaling (q h c : ℝ) (p : Fin 12) (v : Spatial) :
    chargeFixedKappa q (c * h) p v = c * chargeFixedKappa q h p v := by
  unfold chargeFixedKappa
  ring

/-- The fixed normalization carries exactly the declared charge datum: at
nonzero step and crossing factor, equal fixed values force equal declared
charges. -/
theorem chargeFixedKappa_injective (q₁ q₂ h : ℝ) (hh : h ≠ 0) (p : Fin 12)
    (v : Spatial) (hc : crossingFactor p v ≠ 0)
    (heq : chargeFixedKappa q₁ h p v = chargeFixedKappa q₂ h p v) : q₁ = q₂ := by
  unfold chargeFixedKappa at heq
  have h1 := (div_left_inj' hc).mp heq
  have h2 : q₁ * h = q₂ * h := by linarith
  exact mul_right_cancel₀ hh h2

/-- **Without the identification the normalization stays free.**  For every
normalization, the augmented Gauss constraint is satisfiable by a declared
port load (a free declared load absorbs any residual: a trivial absorption
recorded for completeness, carrying no independent freeness content); and
the committed
non-forcing (gauge invariance at every value, separation on the pulse
configuration) holds unchanged. -/
theorem normalization_free_without_identification (h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) :
    (∀ κ : ℝ, ∃ ρ : ℕ → Fin 12 → ℝ, ∀ n,
      realBoundary (electricFieldScaled h A φ n) = ρ n + inducedLoad κ h x n) ∧
    ((∀ (κ h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ χ : ℕ → Fin 12 → ℝ)
        (x : ℕ → Herm2) (N : ℕ),
      interactionCandidate κ h (gaugeTransformA A χ)
          (gaugeTransformPhiScaled h φ χ) x N =
        interactionCandidate κ h A φ x N) ∧
      interactionCandidate 1 1 pulseHistory (fun _ ↦ 0) stepWorldline 1 ≠
        interactionCandidate 0 1 pulseHistory (fun _ ↦ 0) stepWorldline 1) := by
  refine ⟨fun κ ↦ ⟨fun n ↦ realBoundary (electricFieldScaled h A φ n) -
    inducedLoad κ h x n, fun n ↦ ?_⟩, interaction_normalization_not_forced⟩
  exact (sub_add_cancel _ _).symm


/-! ## (c) Worldline sector: the coupled equation of motion -/

theorem spatialDot_sub_right (u v w : Spatial) :
    spatialDot u (v - w) = spatialDot u v - spatialDot u w := by
  unfold spatialDot
  rw [← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  simp only [Pi.sub_apply]
  ring

theorem spatialDot_sub_left (u v w : Spatial) :
    spatialDot (u - v) w = spatialDot u w - spatialDot v w := by
  unfold spatialDot
  rw [← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  simp only [Pi.sub_apply]
  ring

theorem spatialDot_smul_left (c : ℝ) (u w : Spatial) :
    spatialDot (c • u) w = c * spatialDot u w := by
  unfold spatialDot
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

theorem spatialDot_zero_right (u : Spatial) : spatialDot u 0 = 0 := by
  simp [spatialDot]

/-- The embedded seam electric field at step `n`: the scaled seam electric
field pushed to a spatial vector through the candidate seam directions. -/
def embeddedField (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    Spatial :=
  ∑ e : Fin 30, electricFieldScaled h A φ n e • seamVector e

/-- The seam pairing of the electric field against the worldline seam
current is the spatial pairing of the embedded field against the
worldline increment. -/
theorem seamInner_worldline_eq_dot (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (n : ℕ) :
    realSeamInner (electricFieldScaled h A φ n) (worldlineSeamCurrent x n) =
      spatialDot (embeddedField h A φ n) ((x (n + 1)).2 - (x n).2) := by
  unfold realSeamInner worldlineSeamCurrent embeddedField spatialDot
  simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul, Finset.mul_sum,
    Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun i _ ↦ Finset.sum_congr rfl fun e _ ↦ ?_
  ring

/-- The interaction candidate is additive in the worldline. -/
theorem interaction_worldline_linear (κ h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x η : ℕ → Herm2) (M : ℕ) :
    interactionCandidate κ h A φ (x + η) M =
      interactionCandidate κ h A φ x M + interactionCandidate κ h A φ η M := by
  simp only [interactionCandidate_eq_inner, worldlineSeamCurrent_add,
    realSeamInner_add_right, Finset.sum_add_distrib, mul_add]

/-- Interior form of the interaction in an endpoint-fixed worldline
variation: the sum over interior nodes of the earlier-minus-later
difference of the embedded field paired with the spatial variation. -/
theorem interaction_worldline_interior (κ h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (η : ℕ → Herm2) (h0 : η 0 = 0) (hN : η (N + 1) = 0) :
    interactionCandidate κ h A φ η (N + 1) =
      κ * ∑ k ∈ Finset.range N,
        spatialDot (embeddedField h A φ k - embeddedField h A φ (k + 1))
          (η (k + 1)).2 := by
  rw [interactionCandidate_eq_inner]
  congr 1
  simp only [seamInner_worldline_eq_dot, spatialDot_sub_right,
    Finset.sum_sub_distrib, spatialDot_sub_left]
  have hS1 : (∑ n ∈ Finset.range (N + 1),
      spatialDot (embeddedField h A φ n) (η (n + 1)).2) =
      ∑ k ∈ Finset.range N, spatialDot (embeddedField h A φ k) (η (k + 1)).2 := by
    rw [Finset.sum_range_succ, hN]
    simp only [Prod.snd_zero, spatialDot_zero_right, add_zero]
  have hS2 : (∑ n ∈ Finset.range (N + 1),
      spatialDot (embeddedField h A φ n) (η n).2) =
      ∑ k ∈ Finset.range N,
        spatialDot (embeddedField h A φ (k + 1)) (η (k + 1)).2 := by
    rw [Finset.sum_range_succ', h0]
    simp only [Prod.snd_zero, spatialDot_zero_right, add_zero]
  rw [hS1, hS2]

/-- The total action: the coupled field action plus the mass-weighted
worldline action.  The weight `m` is a declared enrichment of the
committed worldline action; the value `m = 1` is the committed joint
action plus the interaction candidate. -/
def totalAction (κ h m : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (x : ℕ → Herm2) : ℝ :=
  coupledAction κ h N A φ ρ J x + m * clockAction (N + 1) x

/-- At unit weight and the record's declared sources the total action is
the committed joint action plus the interaction candidate. -/
theorem totalAction_one (κ : ℝ)
    (Jr : OPH.CommonWorldMaxwellClockJoin.MaxwellClockJoinedArchitecture) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) :
    totalAction κ Jr.scaled.h 1 N A φ Jr.scaled.rho Jr.scaled.J x =
      jointAction Jr N A φ x + interactionCandidate κ Jr.scaled.h A φ x (N + 1) := by
  unfold totalAction coupledAction jointAction
  ring

/-- The coupled worldline residual at the interior node `k + 1`: twice the
declared weight times the difference of consecutive increments, earlier
minus later (the negative of the second difference), plus the purely
spatial vector `κ` times the forward difference of the embedded field. -/
def worldlineResidual (κ h m : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (x : ℕ → Herm2) (k : ℕ) : Herm2 :=
  (2 * m) • ((x (k + 1) - x k) - (x (k + 1 + 1) - x (k + 1))) +
    ((0 : ℝ), κ • (embeddedField h A φ (k + 1) - embeddedField h A φ k))

/-- The residual pairs against a vector as the weighted clock term plus the
interaction term. -/
theorem lorentzB_worldlineResidual (κ h m : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (k : ℕ) (w : Herm2) :
    lorentzB (worldlineResidual κ h m A φ x k) w =
      m * (2 * lorentzB ((x (k + 1) - x k) - (x (k + 1 + 1) - x (k + 1))) w) +
        κ * spatialDot (embeddedField h A φ k - embeddedField h A φ (k + 1)) w.2 := by
  unfold worldlineResidual
  rw [OPH.C1Lorentz.lorentzB_add_left, OPH.C1Lorentz.lorentzB_smul_left]
  have hsp : lorentzB (((0 : ℝ),
      κ • (embeddedField h A φ (k + 1) - embeddedField h A φ k)) : Herm2) w =
      κ * spatialDot (embeddedField h A φ k - embeddedField h A φ (k + 1)) w.2 := by
    unfold lorentzB
    rw [spatialDot_smul_left, spatialDot_sub_left, spatialDot_sub_left]
    ring
  rw [hsp]
  ring

/-- **Exact worldline-sector expansion.**  For variations vanishing at the
window endpoints, the total action moves by the residual paired against
the interior variation values plus the weighted clock action of the
variation. -/
theorem totalAction_worldline_expansion (κ h m : ℝ) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (x η : ℕ → Herm2) (h0 : η 0 = 0) (hN : η (N + 1) = 0) :
    totalAction κ h m N A φ ρ J (x + η) =
      totalAction κ h m N A φ ρ J x +
        (∑ k ∈ Finset.range N,
          lorentzB (worldlineResidual κ h m A φ x k) (η (k + 1))) +
        m * clockAction (N + 1) η := by
  unfold totalAction coupledAction
  rw [interaction_worldline_linear, clockAction_expansion,
    clockFirstVariation_interior N x η h0 hN,
    interaction_worldline_interior κ h N A φ η h0 hN]
  simp only [lorentzB_worldlineResidual, Finset.sum_add_distrib, ← Finset.mul_sum]
  ring


/-- Stationarity of the total action under worldline variations vanishing
at the nodes `0` and `N + 1`: the action moves by exactly the weighted
clock action of the variation. -/
def WorldlineStationary (κ h m : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (x : ℕ → Herm2) : Prop :=
  ∀ η : ℕ → Herm2, η 0 = 0 → η (N + 1) = 0 →
    totalAction κ h m N A φ ρ J (x + η) =
      totalAction κ h m N A φ ρ J x + m * clockAction (N + 1) η

/-- **Worldline stationarity is the vanishing of the coupled residual at
every interior node.**  The single-node variation isolates one residual,
and nondegeneracy of the committed Lorentz pairing forces it to zero. -/
theorem worldline_stationary_iff (κ h m : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (x : ℕ → Herm2) :
    WorldlineStationary κ h m N A φ ρ J x ↔
      ∀ k, k < N → worldlineResidual κ h m A φ x k = 0 := by
  constructor
  · intro hs k hk
    apply lorentzB_nondegenerate
    intro w
    let η : ℕ → Herm2 := fun n ↦ if n = k + 1 then w else 0
    have h0 : η 0 = 0 := if_neg (by omega)
    have hN1 : η (N + 1) = 0 := if_neg (by omega)
    have hexp := totalAction_worldline_expansion κ h m N A φ ρ J x η h0 hN1
    rw [hs η h0 hN1] at hexp
    have hsum : (∑ j ∈ Finset.range N,
        lorentzB (worldlineResidual κ h m A φ x j) (η (j + 1))) = 0 := by
      linarith
    have hterm : ∀ j ∈ Finset.range N,
        lorentzB (worldlineResidual κ h m A φ x j) (η (j + 1)) =
          if j = k then lorentzB (worldlineResidual κ h m A φ x k) w else 0 := by
      intro j _
      by_cases hjk : j = k
      · subst hjk
        rw [if_pos rfl]
        have hv : η (j + 1) = w := if_pos rfl
        rw [hv]
      · rw [if_neg hjk]
        have hv : η (j + 1) = 0 := if_neg (by omega)
        rw [hv, lorentzB_zero_right]
    rw [Finset.sum_congr rfl hterm, Finset.sum_ite_eq',
      if_pos (Finset.mem_range.mpr hk)] at hsum
    exact hsum
  · intro hr η h0 hN
    rw [totalAction_worldline_expansion κ h m N A φ ρ J x η h0 hN]
    have hz : (∑ k ∈ Finset.range N,
        lorentzB (worldlineResidual κ h m A φ x k) (η (k + 1))) = 0 := by
      refine Finset.sum_eq_zero fun k hk ↦ ?_
      rw [hr k (Finset.mem_range.mp hk), lorentzB_zero_left]
    rw [hz, add_zero]

/-- **The coupled equation of motion.**  The residual vanishes exactly when
twice the declared weight times the forward second difference of the path
equals the purely spatial vector `κ` times the forward difference of the
embedded seam electric field. -/
theorem coupled_worldline_equation (κ h m : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (k : ℕ) :
    worldlineResidual κ h m A φ x k = 0 ↔
      (2 * m) • ((x (k + 1 + 1) - x (k + 1)) - (x (k + 1) - x k)) =
        ((0 : ℝ), κ • (embeddedField h A φ (k + 1) - embeddedField h A φ k)) := by
  unfold worldlineResidual
  have hneg : (2 * m) • ((x (k + 1 + 1) - x (k + 1)) - (x (k + 1) - x k)) =
      -((2 * m) • ((x (k + 1) - x k) - (x (k + 1 + 1) - x (k + 1)))) := by
    rw [← smul_neg, neg_sub]
  rw [hneg, neg_eq_iff_add_eq_zero]

/-- Component form: the scalar coordinate has vanishing weighted second
difference; the spatial coordinates obey the field-difference law. -/
theorem coupled_worldline_equation_components (κ h m : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (k : ℕ) :
    worldlineResidual κ h m A φ x k = 0 ↔
      ((2 * m) * ((x (k + 1 + 1) - x (k + 1)) - (x (k + 1) - x k)).1 = 0 ∧
        (2 * m) • ((x (k + 1 + 1) - x (k + 1)) - (x (k + 1) - x k)).2 =
          κ • (embeddedField h A φ (k + 1) - embeddedField h A φ k)) := by
  rw [coupled_worldline_equation, Prod.ext_iff]
  exact Iff.rfl

/-- Zero-coupling regression: the residual is the weighted committed
second-difference clause. -/
theorem worldlineResidual_zero_coupling (h m : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (k : ℕ) :
    worldlineResidual 0 h m A φ x k =
      (2 * m) • ((x (k + 1) - x k) - (x (k + 1 + 1) - x (k + 1))) := by
  unfold worldlineResidual
  rw [zero_smul]
  exact add_zero _

/-- At zero coupling and unit weight, worldline stationarity is the
committed equal-increment clause of the joint action. -/
theorem worldline_stationary_zero_coupling_iff_uniform (h : ℝ) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (x : ℕ → Herm2) :
    WorldlineStationary 0 h 1 N A φ ρ J x ↔
      ∀ k, k < N → x (k + 1 + 1) - x (k + 1) = x (k + 1) - x k := by
  rw [worldline_stationary_iff]
  refine forall_congr' fun k ↦ forall_congr' fun _ ↦ ?_
  rw [worldlineResidual_zero_coupling, smul_eq_zero, sub_eq_zero]
  constructor
  · rintro (h2 | hd)
    · norm_num at h2
    · exact hd.symm
  · intro hd
    exact Or.inr hd.symm

/-- The discrete momentum of the weighted worldline action: twice the
weight times the increment. -/
def discreteMomentum (m : ℝ) (x : ℕ → Herm2) (k : ℕ) : Herm2 :=
  (2 * m) • (x (k + 1) - x k)

/-- On a frame increment the discrete momentum is twice the declared-mass
momentum of `Geometry/MassShellKinematics.lean`: the inertial coefficient
of the declared quadratic action is `2 m`, so the weight `m` equals half
the inertial coefficient and the parameter matching `fourMomentum` is
`2 m`. -/
theorem discreteMomentum_frame (m : ℝ) (frame : OPH.C1Lorentz.FrameHyperboloid)
    (x : ℕ → Herm2) (k : ℕ) (hx : x (k + 1) - x k = (frame : Herm2)) :
    discreteMomentum m x k = (2 : ℝ) • OPH.C1Lorentz.fourMomentum m frame := by
  unfold discreteMomentum OPH.C1Lorentz.fourMomentum
  rw [hx, smul_smul]

/-- Momentum form of the coupled equation of motion: the forward difference
of the discrete momentum equals the purely spatial vector `κ` times the
forward difference of the embedded field. -/
theorem coupled_momentum_law (κ h m : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (k : ℕ) :
    worldlineResidual κ h m A φ x k = 0 ↔
      discreteMomentum m x (k + 1) - discreteMomentum m x k =
        ((0 : ℝ), κ • (embeddedField h A φ (k + 1) - embeddedField h A φ k)) := by
  rw [coupled_worldline_equation]
  unfold discreteMomentum
  rw [smul_sub]

/-! ## (d) Boundary -/

/-- A resting worldline induces no load: zero spatial increment gives zero
induced load at every port, so the identified value is step-local and
increment-proportional. -/
theorem inducedLoad_rest (κ h : ℝ) (x : ℕ → Herm2) (n : ℕ)
    (hx : (x (n + 1)).2 = (x n).2) (p : Fin 12) : inducedLoad κ h x n p = 0 := by
  rw [inducedLoad_apply, hx, sub_self, spatialDot_zero_right]
  ring

/-- The declared charge is a free datum of the identification: every real
value is realized by some normalization on the seam crossing, at nonzero
step.  Nothing in the module selects `q`. -/
theorem declared_charge_unconstrained (h : ℝ) (hh : h ≠ 0) (e : Fin 30) :
    ∀ q : ℝ, ∃ κ : ℝ,
      PortChargeIdentification κ h q (seamCrossingWorldline e) 0 (seamRight e) := by
  intro q
  exact ⟨-(q * h) / seamCrossingConst,
    (seam_crossing_charge_fixed _ h q hh e).mpr rfl⟩

end

end OPH.ChargeFixedInteraction

/- Axiom audit: standard axioms only (`propext`, `Classical.choice`,
`Quot.sound`); no `sorry`, no `native_decide`, no project axiom. -/

#print axioms OPH.ChargeFixedInteraction.seamBoundaryVecZ_eq
#print axioms OPH.ChargeFixedInteraction.realBoundary_seamPaired
#print axioms OPH.ChargeFixedInteraction.coupledAction_reduction
#print axioms OPH.ChargeFixedInteraction.coupledAction_field_variation_interior
#print axioms OPH.ChargeFixedInteraction.ampereResidual_augmented
#print axioms OPH.ChargeFixedInteraction.coupled_stationary_A_iff_ampere
#print axioms OPH.ChargeFixedInteraction.coupled_stationary_phi_iff_gauss
#print axioms OPH.ChargeFixedInteraction.coupled_field_equations
#print axioms OPH.ChargeFixedInteraction.induced_continuity
#print axioms OPH.ChargeFixedInteraction.inducedLoad_total
#print axioms OPH.ChargeFixedInteraction.identification_iff_kappa
#print axioms OPH.ChargeFixedInteraction.kappa_unique_of_identification
#print axioms OPH.ChargeFixedInteraction.seam_crossing_charge_fixed
#print axioms OPH.ChargeFixedInteraction.seam_crossing_left_charge
#print axioms OPH.ChargeFixedInteraction.seamCrossingConst_ne_zero
#print axioms OPH.ChargeFixedInteraction.chargeFixedKappa_step_scaling
#print axioms OPH.ChargeFixedInteraction.chargeFixedKappa_injective
#print axioms OPH.ChargeFixedInteraction.normalization_free_without_identification
#print axioms OPH.ChargeFixedInteraction.totalAction_worldline_expansion
#print axioms OPH.ChargeFixedInteraction.worldline_stationary_iff
#print axioms OPH.ChargeFixedInteraction.coupled_worldline_equation
#print axioms OPH.ChargeFixedInteraction.coupled_worldline_equation_components
#print axioms OPH.ChargeFixedInteraction.worldline_stationary_zero_coupling_iff_uniform
#print axioms OPH.ChargeFixedInteraction.discreteMomentum_frame
#print axioms OPH.ChargeFixedInteraction.coupled_momentum_law
#print axioms OPH.ChargeFixedInteraction.inducedLoad_rest
#print axioms OPH.ChargeFixedInteraction.declared_charge_unconstrained
