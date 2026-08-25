import Geometry.CommonWorldJointAction
import CapFirstLaw
import PhysicalCalibrationImport
import Variational.RealizedHistoryLegendreNoGo

set_option autoImplicit false

open scoped BigOperators

/-!
# Internal energy as inertia: a conditional finite precursor (issues 736, 739)

STATUS.  Precursor of the dynamical half of the mass-energy identity,
conditional on one declared composite action shape.  The kinematic half is
committed in `Geometry/MassShellKinematics.lean` (the scalar coordinate of
the declared-mass momentum equals the mass parameter at rest).  Nothing here
derives the identity; the module locates its dynamical content in one shape
selection and proves the consequence of that selection exactly.

WHAT IS PROVED.

(1) Internal energy ledger.  `modularEnergy p tau` is the modular energy
`⟨K⟩_p = ∑ p x * (-log tau x)` of `Thermodynamics/CapFirstLaw.lean`, bound
directly to that module's objects: `cap_firstLaw_modular` restates the exact
cap first law with the ledger name, and `modularEnergy_nonneg` proves the
ledger nonnegative for a nonnegative state and a sub-unit reference.  The
binding is light (one definition and two lemmas), so the abstract-ledger
fallback is unused.

(2) Composite shape A (declared).  For a declared rest parameter `m`, an
internal energy `E`, and a declared step-indexed potential `F : ℕ → Herm2`,
`compositeAction m E M F x = (m + E) * clockAction M x + forceTerm M F x`,
where `clockAction` is the committed quadratic worldline functional of
`Geometry/CommonWorldJointAction.lean` in the `(+---)` convention and
`forceTerm` pairs the potential with the path increment, normalized by the
polarization constant `2` of the committed functional.  The expansion is
exact (`compositeAction_expansion`), and `compositeStationary_iff_eom`
proves stationarity under endpoint-fixed variations equivalent to the
discrete equation of motion at every interior node:
`(m + E) • secondDifference x j = impulse F j`, with `impulse F j = F j -
F (j + 1)` the negative forward difference of the declared potential at
`j`, a declared orientation of the potential.  The inertial
coefficient multiplying the second difference is exactly `m + E`
(`composite_inertial_coefficient`): conditional on shape A, internal energy
contributes to inertia with coefficient one, in the dimensionless module
convention where the invariant-speed conversion is one.  Rest-frame
corollary: `compositeFourMomentum m E standardFrame` has scalar coordinate
`m + E` (`compositeFourMomentum_rest`), through the committed kinematic
identities.

(3) Non-forcing: composite shape B (declared).  `additiveAction m E M F x =
m * clockAction M x + E * M + forceTerm M F x` enters the same ledger as an
additive per-step term.  Shape B is stationary-well-posed in the committed
sense (`additiveStationary_iff_eom`) with inertial coefficient `m`.  The two
shapes coincide at zero internal energy (`shapes_agree_zero_energy`) and on
every path with unit-Lorentz-square increments, in particular on the joined
record worldline at step duration one (`shapes_agree_unit_worldline`), yet
their equations of motion differ: the explicit history `parabolaPath` with
the explicit potential `parabolaPotential` is stationary for shape A at
`m = 1, E = 1` and fails stationarity for shape B at the same data
(`shape_selection_discriminated`).  The contribution of internal energy to
inertia is therefore selected by the shape; the selection is the open
physical content.  Shape B carries no `E`-dependent dynamics (its additive
term is path independent, so it is shape A at zero energy plus a
constant); the open selection is the slope family `m + λ E`, of which
shapes A and B are the endpoints `λ = 1` and `λ = 0`.

(4) Declared-enrichment status.  The committed realized-history law does
not select a velocity curvature or Legendre map
(`Variational/RealizedHistoryLegendreNoGo.lean`, cited as
`legendre_nonidentifiability_cited`), so the multiplicative shape A and the
additive shape B are both declared enrichments; no theorem here reads either
off a history.

(5) Unit dictionary.  The ledger converts to joules through the declared
tick of `Thermodynamics/PhysicalCalibrationImport.lean`
(`composite_energy_joules`), and two distinct ticks give distinct joule
readings for one nonzero composite energy (`composite_joules_not_forced`).

ROWS TOUCHED.  This module touches the coupled-action row (a declared
composite shape on the kinematics island of the direct-sum action; the
force term is a step-indexed potential, with no worldline-field coupling to
the Maxwell sector) and the internal-energy identification (the ledger is
the modular energy; its identification with a physical internal energy is
the energy-clock receipt, open).  It does NOT discharge the source clock or
duration row (the window index is the declared step index; the step
duration stays declared), the physical spacetime attachment row (`Herm2` is
the declared Lorentz module; the carrier map stays open), or any derived
mass (`m` is a declared parameter).

NEGATIVES CITED.  Only the Legendre no-go is invoked, at its scope: realized
source histories do not select the velocity curvature or Legendre map, so
every Lagrangian shape below is declared.  No clock claim is made, so the
abstract-transition-system rate result is neither invoked nor contradicted.

CONVENTIONS.  The rest parameter `m` is read in tick-energy units by
declaration (the invariant-speed conversion is one), which is what makes
the sum `m + E` with the ledger in nats meaningful.  Signature `(+---)`; `secondDifference x j = (x (j + 2) -
x (j + 1)) - (x (j + 1) - x j)`; `impulse F j = F j - F (j + 1)`, so a
constant potential is pure gauge (`forceTerm_const`); the factor `2` in
`forceTerm` is the polarization constant of `lorentzQ` and is a
normalization of the declared potential.

FALSIFIER.  The module fails if the exact expansion misses a term, if
stationarity of either shape differs from its displayed equation of motion,
if the two shapes disagree at zero energy or on the unit worldline, if the
parabola history fails shape-A stationarity or satisfies shape-B
stationarity, or if the rest-frame scalar coordinate differs from `m + E`.

Axiom audit.  No project axiom, no native decision procedure; the guard
lines at the end show at most `propext`, `Classical.choice`, `Quot.sound`.
-/

namespace OPH.InternalEnergyInertia

open OPH.C1Lorentz OPH.C2Soldering OPH.CommonWorldJointAction
open OPH.Thermodynamics OPH.PhysicalCalibrationImport OPH.Variational

noncomputable section

/-! ## (1) The internal energy ledger, bound to the cap first law -/

/-- The modular energy `⟨K⟩_p` with `K = -log tau`: the internal energy
ledger of a finite thermodynamic system carried along the worldline.  The
expression is the one appearing in `cap_firstLaw_exact`. -/
def modularEnergy {Ω : Type*} [Fintype Ω] (p tau : Ω → ℝ) : ℝ :=
  ∑ x, p x * (-Real.log (tau x))

/-- The exact cap first law in ledger names: the entropy difference to the
reference equals the ledger difference minus the relative entropy. -/
theorem cap_firstLaw_modular {Ω : Type*} [Fintype Ω] [DecidableEq Ω]
    (p tau : Ω → ℝ) (hτ : ∀ x, 0 < tau x) :
    shannon p - shannon tau =
      modularEnergy p tau - modularEnergy tau tau - kl p tau :=
  cap_firstLaw_exact p tau hτ

/-- The ledger is nonnegative for a nonnegative state and a reference with
values in `(0, 1]`. -/
theorem modularEnergy_nonneg {Ω : Type*} [Fintype Ω] (p tau : Ω → ℝ)
    (hp : ∀ x, 0 ≤ p x) (hτ0 : ∀ x, 0 < tau x) (hτ1 : ∀ x, tau x ≤ 1) :
    0 ≤ modularEnergy p tau := by
  unfold modularEnergy
  refine Finset.sum_nonneg fun x _ => mul_nonneg (hp x) ?_
  rw [neg_nonneg]
  exact Real.log_nonpos (hτ0 x).le (hτ1 x)

/-- The ledger vanishes for the reference concentrated on a point distribution and
read at itself only when that mass is one; the general zero clause: a state
supported where the reference equals one has zero ledger. -/
theorem modularEnergy_zero_of_support {Ω : Type*} [Fintype Ω]
    (p tau : Ω → ℝ) (h : ∀ x, p x ≠ 0 → tau x = 1) :
    modularEnergy p tau = 0 := by
  unfold modularEnergy
  refine Finset.sum_eq_zero fun x _ => ?_
  by_cases hx : p x = 0
  · rw [hx, zero_mul]
  · rw [h x hx, Real.log_one, neg_zero, mul_zero]

/-! ## (2) The declared potential term -/

/-- The path whose increments are the declared potential: its partial
sums.  It lets the potential term reuse the committed first-variation
algebra of the worldline functional. -/
def potentialPath (F : ℕ → Herm2) (k : ℕ) : Herm2 :=
  ∑ i ∈ Finset.range k, F i

theorem potentialPath_increment (F : ℕ → Herm2) (k : ℕ) :
    potentialPath F (k + 1) - potentialPath F k = F k := by
  unfold potentialPath
  rw [Finset.sum_range_succ]
  abel

/-- The declared potential term: the step-indexed potential paired with the
path increment at every window step, normalized by the polarization
constant `2` of the committed quadratic functional.  The normalization is a
convention on the declared potential. -/
def forceTerm (M : ℕ) (F x : ℕ → Herm2) : ℝ :=
  ∑ k ∈ Finset.range M, 2 * lorentzB (F k) (x (k + 1) - x k)

/-- The potential term is the committed first variation of the worldline
functional at the potential path. -/
theorem forceTerm_eq_clockFirstVariation (M : ℕ) (F x : ℕ → Herm2) :
    forceTerm M F x = clockFirstVariation M (potentialPath F) x := by
  unfold forceTerm clockFirstVariation
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [potentialPath_increment]

/-- The potential term is additive in the path. -/
theorem forceTerm_add (M : ℕ) (F x η : ℕ → Herm2) :
    forceTerm M F (x + η) = forceTerm M F x + forceTerm M F η := by
  unfold forceTerm
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  have hstep : (x + η) (k + 1) - (x + η) k =
      (x (k + 1) - x k) + (η (k + 1) - η k) := by
    simp only [Pi.add_apply]
    abel
  rw [hstep, lorentzB_add_right]
  ring

/-- Interior form of the potential term at an endpoint-fixed variation: the
backward difference of the potential paired with the variation at the
interior nodes. -/
theorem forceTerm_interior (N : ℕ) (F η : ℕ → Herm2) (h0 : η 0 = 0)
    (hN : η (N + 1) = 0) :
    forceTerm (N + 1) F η =
      ∑ k ∈ Finset.range N, 2 * lorentzB (F k - F (k + 1)) (η (k + 1)) := by
  rw [forceTerm_eq_clockFirstVariation,
    clockFirstVariation_interior N (potentialPath F) η h0 hN]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [potentialPath_increment, potentialPath_increment]

/-- A constant potential is pure gauge: its term telescopes to the window
endpoints and moves no endpoint-fixed variation. -/
theorem forceTerm_const (M : ℕ) (c : Herm2) (x : ℕ → Herm2) :
    forceTerm M (fun _ => c) x = 2 * lorentzB c (x M - x 0) := by
  unfold forceTerm
  have hterm : ∀ k ∈ Finset.range M,
      2 * lorentzB c (x (k + 1) - x k) =
        2 * lorentzB c (x (k + 1)) - 2 * lorentzB c (x k) := by
    intro k _
    rw [lorentzB_sub_right]
    ring
  rw [Finset.sum_congr rfl hterm,
    Finset.sum_range_sub (fun n => 2 * lorentzB c (x n)) M, lorentzB_sub_right]
  ring

/-! ## (3) Composite shape A: internal energy multiplies the worldline
functional -/

/-- The second difference of a path at the interior node `j + 1`. -/
def secondDifference (x : ℕ → Herm2) (j : ℕ) : Herm2 :=
  (x (j + 1 + 1) - x (j + 1)) - (x (j + 1) - x j)

/-- The impulse of the declared potential at the interior node `j + 1`: the
negative forward difference `F j - F (j + 1)` at `j`, a declared
orientation of the potential. -/
def impulse (F : ℕ → Herm2) (j : ℕ) : Herm2 := F j - F (j + 1)

/-- The inertial coefficient of shape A. -/
def inertialCoefficient (m E : ℝ) : ℝ := m + E

/-- **Composite shape A (declared).**  The rest parameter plus the internal
energy multiplies the committed worldline functional; the declared
potential term is added.  The shape is a declared enrichment
(`legendre_nonidentifiability_cited`). -/
def compositeAction (m E : ℝ) (M : ℕ) (F x : ℕ → Herm2) : ℝ :=
  (m + E) * clockAction M x + forceTerm M F x

/-- The first variation of shape A. -/
def compositeFirstVariation (m E : ℝ) (M : ℕ) (F x η : ℕ → Herm2) : ℝ :=
  (m + E) * clockFirstVariation M x η + forceTerm M F η

/-- Exact expansion of shape A: displaced action equals action, plus first
variation, plus the scaled clock action of the displacement. -/
theorem compositeAction_expansion (m E : ℝ) (M : ℕ) (F x η : ℕ → Herm2) :
    compositeAction m E M F (x + η) =
      compositeAction m E M F x + compositeFirstVariation m E M F x η +
        (m + E) * clockAction M η := by
  unfold compositeAction compositeFirstVariation
  rw [clockAction_expansion, forceTerm_add]
  ring

/-- Stationarity of shape A on the window `0, …, M`: every endpoint-fixed
variation moves the action by exactly the scaled quadratic remainder. -/
def CompositeStationary (m E : ℝ) (M : ℕ) (F x : ℕ → Herm2) : Prop :=
  ∀ η : ℕ → Herm2, η 0 = 0 → η M = 0 →
    compositeAction m E M F (x + η) =
      compositeAction m E M F x + (m + E) * clockAction M η

theorem compositeStationary_firstVariation_zero {m E : ℝ} {M : ℕ}
    {F x : ℕ → Herm2} (hs : CompositeStationary m E M F x) {η : ℕ → Herm2}
    (h0 : η 0 = 0) (hM : η M = 0) :
    compositeFirstVariation m E M F x η = 0 := by
  have hexp := compositeAction_expansion m E M F x η
  have hst := hs η h0 hM
  linarith

/-- The equation-of-motion residual at the interior node `k + 1`. -/
def eomResidual (m E : ℝ) (F x : ℕ → Herm2) (k : ℕ) : Herm2 :=
  impulse F k - (m + E) • secondDifference x k

theorem residual_eq (m E : ℝ) (F x : ℕ → Herm2) (k : ℕ) :
    (m + E) • ((x (k + 1) - x k) - (x (k + 1 + 1) - x (k + 1))) +
        (F k - F (k + 1)) = eomResidual m E F x k := by
  unfold eomResidual impulse secondDifference
  rw [show (x (k + 1) - x k) - (x (k + 1 + 1) - x (k + 1)) =
      -((x (k + 1 + 1) - x (k + 1)) - (x (k + 1) - x k)) by abel, smul_neg]
  abel

/-- Interior form of the first variation of shape A: the residual paired
with the variation at every interior node. -/
theorem compositeFirstVariation_interior (m E : ℝ) (N : ℕ)
    (F x η : ℕ → Herm2) (h0 : η 0 = 0) (hN : η (N + 1) = 0) :
    compositeFirstVariation m E (N + 1) F x η =
      ∑ k ∈ Finset.range N, 2 * lorentzB (eomResidual m E F x k) (η (k + 1)) := by
  unfold compositeFirstVariation
  rw [clockFirstVariation_interior N x η h0 hN, forceTerm_interior N F η h0 hN,
    Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [← residual_eq, lorentzB_add_left, lorentzB_smul_left]
  ring

/-- Shape A stationarity forces the discrete equation of motion at every
interior node: the variation supported at one node paired against the
residual vanishes, and nondegeneracy kills the residual. -/
theorem compositeStationary_eom {m E : ℝ} {N : ℕ} {F x : ℕ → Herm2}
    (hs : CompositeStationary m E (N + 1) F x) {j : ℕ} (hj : j < N) :
    (m + E) • secondDifference x j = impulse F j := by
  have key : ∀ w : Herm2, lorentzB (eomResidual m E F x j) w = 0 := by
    intro w
    let η : ℕ → Herm2 := fun n => if n = j + 1 then w else 0
    have h0 : η 0 = 0 := if_neg (by omega)
    have hN1 : η (N + 1) = 0 := if_neg (by omega)
    have hfv := compositeStationary_firstVariation_zero hs h0 hN1
    rw [compositeFirstVariation_interior m E N F x η h0 hN1] at hfv
    have hterm : ∀ k ∈ Finset.range N,
        2 * lorentzB (eomResidual m E F x k) (η (k + 1)) =
          if k = j then 2 * lorentzB (eomResidual m E F x j) w else 0 := by
      intro k _
      by_cases hkj : k = j
      · subst hkj
        rw [if_pos rfl]
        have hv : η (k + 1) = w := if_pos rfl
        rw [hv]
      · rw [if_neg hkj]
        have hv : η (k + 1) = 0 := if_neg (by omega)
        rw [hv, lorentzB_zero_right, mul_zero]
    rw [Finset.sum_congr rfl hterm, Finset.sum_ite_eq',
      if_pos (Finset.mem_range.mpr hj)] at hfv
    linarith
  have hD := lorentzB_nondegenerate key
  unfold eomResidual at hD
  exact (sub_eq_zero.mp hD).symm

/-- The discrete equation of motion at every interior node gives shape A
stationarity: the first variation vanishes termwise. -/
theorem compositeStationary_of_eom {m E : ℝ} {N : ℕ} {F x : ℕ → Herm2}
    (heom : ∀ j, j < N → (m + E) • secondDifference x j = impulse F j) :
    CompositeStationary m E (N + 1) F x := by
  intro η h0 hN
  have hfv : compositeFirstVariation m E (N + 1) F x η = 0 := by
    rw [compositeFirstVariation_interior m E N F x η h0 hN]
    refine Finset.sum_eq_zero fun k hk => ?_
    have hres : eomResidual m E F x k = 0 := by
      unfold eomResidual
      rw [heom k (Finset.mem_range.mp hk), sub_self]
    rw [hres, lorentzB_zero_left, mul_zero]
  have hexp := compositeAction_expansion m E (N + 1) F x η
  rw [hexp, hfv, add_zero]

/-- **Headline (conditional on shape A).**  On a positive window, shape A
stationarity under endpoint-fixed variations is equivalent to the discrete
equation of motion `(m + E) • secondDifference x j = impulse F j` at every
interior node.  The inertial coefficient multiplying the second difference
is exactly `m + E`. -/
theorem compositeStationary_iff_eom (m E : ℝ) (N : ℕ) (F x : ℕ → Herm2) :
    CompositeStationary m E (N + 1) F x ↔
      ∀ j, j < N → (m + E) • secondDifference x j = impulse F j :=
  ⟨fun hs _ hj => compositeStationary_eom hs hj, compositeStationary_of_eom⟩

/-- The same equivalence with the coefficient named: internal energy enters
the inertial coefficient with slope one. -/
theorem composite_inertial_coefficient (m E : ℝ) (N : ℕ) (F x : ℕ → Herm2) :
    (CompositeStationary m E (N + 1) F x ↔
      ∀ j, j < N → inertialCoefficient m E • secondDifference x j = impulse F j) ∧
    inertialCoefficient m E - inertialCoefficient m 0 = E := by
  refine ⟨compositeStationary_iff_eom m E N F x, ?_⟩
  unfold inertialCoefficient
  ring

/-! ## (4) Free shape-A worldlines and the rest-frame corollary -/

theorem impulse_zero (j : ℕ) : impulse (fun _ => (0 : Herm2)) j = 0 := by
  unfold impulse
  simp

/-- An equal-increment path has vanishing second difference at every
interior node of the window. -/
theorem secondDifference_uniform {M : ℕ} {x : ℕ → Herm2}
    (hu : ∀ k, k < M → x (k + 1) - x k = x 1 - x 0) {j : ℕ}
    (hj : j + 1 < M) : secondDifference x j = 0 := by
  unfold secondDifference
  rw [hu (j + 1) hj, hu j (by omega), sub_self]

/-- The joined record's worldline is stationary for the free shape A at
every rest parameter and every internal energy: the free equation of
motion is the equal-increment clause, which is energy-blind. -/
theorem joinedPath_composite_free_stationary
    (J : OPH.CommonWorldMaxwellClockJoin.MaxwellClockJoinedArchitecture)
    (m E : ℝ) (N : ℕ) :
    CompositeStationary m E (N + 1) (fun _ => 0) (joinedPath J) := by
  refine compositeStationary_of_eom fun j hj => ?_
  rw [secondDifference_uniform (M := N + 1) (fun k _ => joinedPath_uniform J k)
    (by omega), smul_zero, impulse_zero]

/-- The composite four-momentum: the committed declared-mass momentum at
the parameter `m + E`. -/
def compositeFourMomentum (m E : ℝ) (frame : FrameHyperboloid) : Herm2 :=
  fourMomentum (m + E) frame

/-- **Rest-frame identity at the composite parameter.**  The composite
momentum is declared at parameter `m + E`, so at the standard frame its
scalar coordinate is `m + E` by the committed kinematic identity of
`Geometry/MassShellKinematics.lean`; this is the committed identity at the
composite parameter, and `composite_momentum_of_increment` is the only
bridge from the shape-A dynamics to it. -/
theorem compositeFourMomentum_rest (m E : ℝ) :
    (compositeFourMomentum m E standardFrame).1 = m + E :=
  standardFrame_time_eq_mass (m + E)

/-- The composite momentum lies on the shell of the composite parameter. -/
theorem compositeFourMomentum_shell (m E : ℝ) (frame : FrameHyperboloid) :
    lorentzQ (compositeFourMomentum m E frame) = (m + E) ^ 2 :=
  lorentzQ_fourMomentum (m + E) frame

/-- The composite mass-shell identity in coordinates. -/
theorem compositeFourMomentum_time_sq (m E : ℝ) (frame : FrameHyperboloid) :
    (compositeFourMomentum m E frame).1 ^ 2 =
      (m + E) ^ 2 + spatialNormSq (compositeFourMomentum m E frame).2 :=
  fourMomentum_time_sq_eq_mass_sq_add_spatial (m + E) frame

/-- Rest characterization for a positive composite parameter. -/
theorem compositeFourMomentum_rest_iff {m E : ℝ} (h : 0 < m + E)
    (frame : FrameHyperboloid) :
    (compositeFourMomentum m E frame).1 = m + E ↔
      (compositeFourMomentum m E frame).2 = 0 :=
  fourMomentum_time_eq_mass_iff_spatial_zero h frame

/-- Link between shape A and the composite momentum: the inertial
coefficient times an increment along a frame direction is the increment
scale times the composite momentum. -/
theorem composite_momentum_of_increment (m E δ : ℝ) (frame : FrameHyperboloid) :
    (m + E) • (δ • (frame : Herm2)) = δ • compositeFourMomentum m E frame := by
  unfold compositeFourMomentum fourMomentum
  rw [smul_comm]

/-- On the joined record's worldline the coefficient times the increment is
the declared step duration times the composite momentum at the record's
frame. -/
theorem joinedPath_composite_momentum
    (J : OPH.CommonWorldMaxwellClockJoin.MaxwellClockJoinedArchitecture)
    (m E : ℝ) (n : ℕ) :
    (m + E) • (joinedPath J (n + 1) - joinedPath J n) =
      J.stepDuration • compositeFourMomentum m E J.frame := by
  rw [joinedPath_increment, composite_momentum_of_increment]

/-! ## (5) Composite shape B: internal energy as an additive per-step term -/

/-- **Composite shape B (declared).**  The rest parameter alone multiplies
the worldline functional; the ledger enters as the additive term `E * M`,
accrued once per window step; the potential term is unchanged.  Also a
declared enrichment. -/
def additiveAction (m E : ℝ) (M : ℕ) (F x : ℕ → Herm2) : ℝ :=
  m * clockAction M x + E * (M : ℝ) + forceTerm M F x

theorem additiveAction_eq_composite (m E : ℝ) (M : ℕ) (F x : ℕ → Herm2) :
    additiveAction m E M F x = compositeAction m 0 M F x + E * (M : ℝ) := by
  unfold additiveAction compositeAction
  ring

/-- Stationarity of shape B: every endpoint-fixed variation moves the action
by exactly the rest-parameter-scaled quadratic remainder. -/
def AdditiveStationary (m E : ℝ) (M : ℕ) (F x : ℕ → Herm2) : Prop :=
  ∀ η : ℕ → Herm2, η 0 = 0 → η M = 0 →
    additiveAction m E M F (x + η) =
      additiveAction m E M F x + m * clockAction M η

/-- Shape B stationarity is shape A stationarity at zero internal energy:
the additive term cancels from every variation. -/
theorem additiveStationary_iff_composite_zero (m E : ℝ) (M : ℕ)
    (F x : ℕ → Herm2) :
    AdditiveStationary m E M F x ↔ CompositeStationary m 0 M F x := by
  unfold AdditiveStationary CompositeStationary
  refine forall_congr' fun η => imp_congr Iff.rfl (imp_congr Iff.rfl ?_)
  rw [additiveAction_eq_composite, additiveAction_eq_composite, add_zero]
  constructor <;> intro h <;> linarith

/-- **Shape B is stationary-well-posed with inertial coefficient `m`.** -/
theorem additiveStationary_iff_eom (m E : ℝ) (N : ℕ) (F x : ℕ → Herm2) :
    AdditiveStationary m E (N + 1) F x ↔
      ∀ j, j < N → m • secondDifference x j = impulse F j := by
  rw [additiveStationary_iff_composite_zero, compositeStationary_iff_eom]
  simp only [add_zero]

/-- The two shapes coincide at zero internal energy. -/
theorem shapes_agree_zero_energy (m : ℝ) (M : ℕ) (F x : ℕ → Herm2) :
    compositeAction m 0 M F x = additiveAction m 0 M F x := by
  unfold compositeAction additiveAction
  ring

/-- A path with unit-Lorentz-square increments has clock action equal to the
window length. -/
theorem clockAction_unit_increments {M : ℕ} {x : ℕ → Herm2}
    (hu : ∀ k, k < M → lorentzQ (x (k + 1) - x k) = 1) :
    clockAction M x = (M : ℝ) := by
  unfold clockAction
  rw [Finset.sum_congr rfl fun k hk => hu k (Finset.mem_range.mp hk)]
  simp

/-- The two shapes coincide on every path with unit-Lorentz-square
increments. -/
theorem shapes_agree_unit_worldline (m E : ℝ) (M : ℕ) (F x : ℕ → Herm2)
    (hu : ∀ k, k < M → lorentzQ (x (k + 1) - x k) = 1) :
    compositeAction m E M F x = additiveAction m E M F x := by
  unfold compositeAction additiveAction
  rw [clockAction_unit_increments hu]
  ring

/-- At step duration one the joined record's worldline has unit-Lorentz-
square increments. -/
theorem joinedPath_unit_increments
    (J : OPH.CommonWorldMaxwellClockJoin.MaxwellClockJoinedArchitecture)
    (hδ : J.stepDuration = 1) (k : ℕ) :
    lorentzQ (joinedPath J (k + 1) - joinedPath J k) = 1 := by
  rw [joinedPath_increment, hδ, one_smul]
  exact J.frame.2.1

/-- The two shapes coincide on the committed joined witness's worldline for
every rest parameter, internal energy, window, and potential. -/
theorem shapes_agree_committed_worldline (m E : ℝ) (M : ℕ) (F : ℕ → Herm2) :
    compositeAction m E M F
        (joinedPath OPH.CommonWorldMaxwellClockJoin.committedJoinedWitness) =
      additiveAction m E M F
        (joinedPath OPH.CommonWorldMaxwellClockJoin.committedJoinedWitness) :=
  shapes_agree_unit_worldline m E M F _ fun k _ =>
    joinedPath_unit_increments
      OPH.CommonWorldMaxwellClockJoin.committedJoinedWitness rfl k

/-! ## (6) The discriminating receipt: shape selection is the content -/

/-- The explicit history `x j = (j ^ 2, 0)`: constant second difference
`(2, 0)`. -/
def parabolaPath (j : ℕ) : Herm2 := (((j : ℝ) ^ 2), (0 : Spatial))

/-- The explicit potential `F j = (-4 j, 0)`: constant impulse `(4, 0)`. -/
def parabolaPotential (j : ℕ) : Herm2 := ((-4 * (j : ℝ)), (0 : Spatial))

theorem parabolaPath_secondDifference (j : ℕ) :
    secondDifference parabolaPath j = ((2 : ℝ), (0 : Spatial)) := by
  unfold secondDifference parabolaPath
  refine Prod.ext ?_ ?_
  · simp only [Prod.fst_sub]
    push_cast
    ring
  · funext i
    simp

theorem parabolaPotential_impulse (j : ℕ) :
    impulse parabolaPotential j = ((4 : ℝ), (0 : Spatial)) := by
  unfold impulse parabolaPotential
  refine Prod.ext ?_ ?_
  · simp only [Prod.fst_sub]
    push_cast
    ring
  · funext i
    simp

/-- At `m = 1, E = 1` the parabola history with the parabola potential is
stationary for shape A on every window. -/
theorem parabola_shapeA_stationary (N : ℕ) :
    CompositeStationary 1 1 (N + 1) parabolaPotential parabolaPath := by
  refine compositeStationary_of_eom fun j _ => ?_
  rw [parabolaPath_secondDifference, parabolaPotential_impulse]
  refine Prod.ext ?_ ?_
  · simp only [Prod.smul_fst, smul_eq_mul]
    norm_num
  · funext i
    simp

/-- At the same data the parabola history fails shape B stationarity on the
window `0, …, 2`: shape B demands `1 • (2, 0) = (4, 0)`. -/
theorem parabola_shapeB_not_stationary :
    ¬ AdditiveStationary 1 1 2 parabolaPotential parabolaPath := by
  intro hs
  have h := (additiveStationary_iff_eom 1 1 1 parabolaPotential
    parabolaPath).mp hs 0 Nat.zero_lt_one
  rw [parabolaPath_secondDifference, parabolaPotential_impulse] at h
  have h1 := congrArg Prod.fst h
  simp only [Prod.smul_fst, smul_eq_mul] at h1
  norm_num at h1

/-- At zero internal energy the parabola history also fails shape A on the
window `0, …, 2`: the internal energy is what makes shape A stationary. -/
theorem parabola_shapeA_zero_energy_not_stationary :
    ¬ CompositeStationary 1 0 2 parabolaPotential parabolaPath := by
  intro hs
  have h := (compositeStationary_iff_eom 1 0 1 parabolaPotential
    parabolaPath).mp hs 0 Nat.zero_lt_one
  rw [parabolaPath_secondDifference, parabolaPotential_impulse] at h
  have h1 := congrArg Prod.fst h
  simp only [Prod.smul_fst, smul_eq_mul] at h1
  norm_num at h1

/-- **Non-forcing receipt.**  One explicit history and one explicit
potential: stationary for shape A at `m = 1, E = 1`; not stationary for
shape B at the same data; not stationary for shape A at `E = 0`.  Both
shapes agree at zero energy and on the committed witness worldline.  The
contribution of internal energy to inertia is selected by the shape. -/
theorem shape_selection_discriminated :
    CompositeStationary 1 1 2 parabolaPotential parabolaPath ∧
    ¬ AdditiveStationary 1 1 2 parabolaPotential parabolaPath ∧
    ¬ CompositeStationary 1 0 2 parabolaPotential parabolaPath ∧
    (∀ (m : ℝ) (M : ℕ) (F x : ℕ → Herm2),
      compositeAction m 0 M F x = additiveAction m 0 M F x) ∧
    (∀ (m E : ℝ) (M : ℕ) (F : ℕ → Herm2),
      compositeAction m E M F
          (joinedPath OPH.CommonWorldMaxwellClockJoin.committedJoinedWitness) =
        additiveAction m E M F
          (joinedPath OPH.CommonWorldMaxwellClockJoin.committedJoinedWitness)) :=
  ⟨parabola_shapeA_stationary 1, parabola_shapeB_not_stationary,
    parabola_shapeA_zero_energy_not_stationary, shapes_agree_zero_energy,
    shapes_agree_committed_worldline⟩

/-! ## (7) Declared-enrichment status, unit dictionary, row register -/

/-- The Legendre no-go at its scope, two conjuncts re-exported: the
canonical extension of the realized source action has no velocity solver,
and the regular enrichments at curvatures `1` and `2` are distinct
Lagrangians.  Their agreement on every realized history is the first
conjunct of the cited receipt and is used from there.  Shapes A and B
above are therefore declared enrichments; no theorem reads either off a
history. -/
theorem legendre_nonidentifiability_cited :
    (¬ ∃ vel : ℝ → ℝ → ℝ, SolvesMomentum chainLogLagrangian vel) ∧
      chainCurvedLagrangian 1 ≠ chainCurvedLagrangian 2 :=
  ⟨realizedHistory_legendre_nonidentifiability_receipt.2.1,
    realizedHistory_legendre_nonidentifiability_receipt.2.2.2.1⟩

/-- The composite energy converts to joules through the declared tick,
additively in the rest parameter and the ledger. -/
theorem composite_energy_joules (a : SIAnchors) (cal : ClockCalibration)
    (m E : ℝ) :
    energyJoules a cal (m + E) = energyJoules a cal m + energyJoules a cal E := by
  unfold energyJoules
  ring

/-- Two distinct declared ticks give distinct joule readings of one nonzero
composite energy: the unit conversion is an import. -/
theorem composite_joules_not_forced (a : SIAnchors)
    (cal cal' : ClockCalibration) (hne : cal.tau ≠ cal'.tau) {m E : ℝ}
    (h : m + E ≠ 0) :
    energyJoules a cal (m + E) ≠ energyJoules a cal' (m + E) := by
  intro heq
  have key : hbarOverTau a cal = hbarOverTau a cal' := mul_left_cancel₀ h heq
  unfold hbarOverTau at key
  rw [div_eq_div_iff cal.tau_pos.ne' cal'.tau_pos.ne'] at key
  have hc : (a.planckConstant : ℝ) / (2 * Real.pi) ≠ 0 :=
    (div_pos (planckConstant_real_pos a) (mul_pos two_pos Real.pi_pos)).ne'
  exact hne (mul_left_cancel₀ hc key).symm

/-- The three open rows on time, energy, and mass. -/
inductive OpenRow : Type
  /-- Source clock and duration. -/
  | sourceClock
  /-- Physical spacetime attachment. -/
  | spacetimeAttachment
  /-- Coupled action. -/
  | coupledAction

/-- The row this module touches: the coupled action, through a declared
composite shape.  The list is a register label, not a discharge. -/
def touchedRows : List OpenRow := [OpenRow.coupledAction]

/-- The rows this module discharges: none. -/
def dischargedRows : List OpenRow := []

theorem dischargedRows_empty : dischargedRows = [] := rfl

end

end OPH.InternalEnergyInertia

/- Axiom audit: expected at most `propext`, `Classical.choice`, `Quot.sound`
per line; no native decision procedure. -/

#print axioms OPH.InternalEnergyInertia.cap_firstLaw_modular
#print axioms OPH.InternalEnergyInertia.modularEnergy_nonneg
#print axioms OPH.InternalEnergyInertia.forceTerm_interior
#print axioms OPH.InternalEnergyInertia.forceTerm_const
#print axioms OPH.InternalEnergyInertia.compositeAction_expansion
#print axioms OPH.InternalEnergyInertia.compositeFirstVariation_interior
#print axioms OPH.InternalEnergyInertia.compositeStationary_eom
#print axioms OPH.InternalEnergyInertia.compositeStationary_of_eom
#print axioms OPH.InternalEnergyInertia.compositeStationary_iff_eom
#print axioms OPH.InternalEnergyInertia.composite_inertial_coefficient
#print axioms OPH.InternalEnergyInertia.joinedPath_composite_free_stationary
#print axioms OPH.InternalEnergyInertia.compositeFourMomentum_rest
#print axioms OPH.InternalEnergyInertia.compositeFourMomentum_shell
#print axioms OPH.InternalEnergyInertia.compositeFourMomentum_rest_iff
#print axioms OPH.InternalEnergyInertia.joinedPath_composite_momentum
#print axioms OPH.InternalEnergyInertia.additiveStationary_iff_eom
#print axioms OPH.InternalEnergyInertia.shapes_agree_zero_energy
#print axioms OPH.InternalEnergyInertia.shapes_agree_unit_worldline
#print axioms OPH.InternalEnergyInertia.shapes_agree_committed_worldline
#print axioms OPH.InternalEnergyInertia.parabola_shapeA_stationary
#print axioms OPH.InternalEnergyInertia.parabola_shapeB_not_stationary
#print axioms OPH.InternalEnergyInertia.parabola_shapeA_zero_energy_not_stationary
#print axioms OPH.InternalEnergyInertia.shape_selection_discriminated
#print axioms OPH.InternalEnergyInertia.legendre_nonidentifiability_cited
#print axioms OPH.InternalEnergyInertia.composite_energy_joules
#print axioms OPH.InternalEnergyInertia.composite_joules_not_forced
