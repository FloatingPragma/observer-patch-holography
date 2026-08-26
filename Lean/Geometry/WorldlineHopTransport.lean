import Geometry.PortChargeMinimalCoupling

set_option autoImplicit false

open scoped BigOperators

namespace OPH.WorldlineHopTransport

open OPH.ScreenCarrierMapCandidate
open OPH.SeamCurrentCarrierQuotient
open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.CarrierDynamicsCompatibility
open OPH.ChargeFixedInteraction
open OPH.CommonWorldJointAction
open OPH.PortChargeMinimalCoupling
open OPH.C1Lorentz (Spatial Herm2 spatialDot spatialNormSq lorentzB lorentzQ)

/-!
# Transport of a Lorentz-module worldline to a hopping port path

STATUS.  Candidate module on the coupled-action row.  The committed
E-paired interaction (`Geometry/ChargeFixedInteraction.lean`) is sourced
by a `Herm2` worldline; the minimal-coupling alternative
(`Geometry/PortChargeMinimalCoupling.lean`) is sourced by a declared port
hop.  This module declares the class of worldlines that step along the
embedded seam directions of the carrier map candidate and transports each
of them to a hopping port path, so that the two source descriptions are
compared step by step on every such worldline, and joins the transported
monopole route to the committed clock action.  Every numeric value is
exact in `ℤ[φ]` or `ℝ`; the remaining statements are identities, an
inequality, and an inequation on the finite complex.

WHAT IS PROVED.

1. Declared seam-step worldlines.  A seam-step worldline is a start port
   and a step sequence, each step a rest or a signed seam index, with
   every seam step admissible at the port reached (`SeamStepWorldline`).
   It generates a `Herm2` path whose time coordinate advances by a
   declared unit `τ` per step and whose spatial coordinate is the
   cumulative sum of the signed seam vectors (`generatedPath`), and a
   hopping path whose consecutive ports are equal or joined by the seam
   stepped (`hoppingPath`, adjacency `hop_of_admissible`).  The exact seam
   vector is the ray difference of the seam endpoints
   (`seamVectorZ_eq_ray_difference`, definitional).  The port reached at
   step `n` is read off the spatial coordinate through the ray map:
   `candidateRayZ (port n) - candidateRayZ start` is the cumulative
   spatial coordinate in `ℤ[φ]³` (`ray_port_eq_spatialZ`), and the real
   form holds for the generated path (`ray_port_eq_spatial`); the rays
   are injective (`candidateRayZ_pairwise_distinct`, the corpus theorem
   restated), so the port is the unique
   port whose ray difference matches the spatial coordinate
   (`port_determined_by_spatial`).
2. Source transport.  At every step of every seam-step worldline, at both
   endpoints of the seam stepped, the E-paired induced load equals
   `-(κ / h) (12 - 4φ)` times the forward step difference of the unit
   hopping load of the transported path (`transport_load_endpoints`); at a
   rest step both vanish (`transport_load_rest`).  At the charge-fixed
   normalization `κ = -(q h) / (12 - 4φ)` and `h ≠ 0` the endpoint loads equal the
   hopping step differences of charge `q` at every step
   (`transport_load_charge_fixed`).  The agreement is endpoint-local, as
   exhibited in `PortChargeMinimalCoupling.bridge_off_endpoint_exhibit`;
   no endpoint agreement is claimed elsewhere.
3. Current transport.  The E-paired induced current at step `m` is
   `κ / h²` times the step difference of the worldline seam current
   (`inducedCurrent_eq_step_difference`).  On a generated path the
   worldline seam current at step `n` is the exact `ℤ[φ]` pairing of each
   seam vector with the signed seam vector stepped
   (`worldlineSeamCurrent_generated`); at the stepped seam it is
   `± 4`, the seam-vector norm squared `(4, 0)` computed by kernel
   evaluation (`stepNormSqZ_forward`, `worldlineSeamCurrent_at_step`);
   it is supported off the stepped seam (`seamGram_off_diagonal_exhibit`,
   value `2` at the seam pair `(0, 1)`), unlike the hopping current, and
   the antipodal seam carries the same vector, so the value `4` appears
   there as well.  The
   hopping current is `-(q / h)` times the projection of the worldline seam
   current onto the stepped seam divided by that norm
   (`hoppingCurrent_eq_projection`).
4. Joint action.  At `h ≠ 0`, the declared transported action is the monopole coupled
   action of the transported hopping path plus the committed clock action
   of the generated path (`transportedAction`).  Its field-sector
   stationarity is the committed scaled Ampere update and Gauss constraint
   at the augmented hopping sources (`transported_field_equations`).  The
   clock action of the generated path over `M` steps is the sum of
   `τ² - 4` at seam steps and `τ²` at rest steps
   (`clockAction_generated`); with no rest step in the window it is
   `M (τ² - 4)` (`clockAction_generated_no_rest`).  A seam step is timelike
   exactly when `4 < τ²` (`seam_step_timelike_iff`); the threshold is the
   exact `ℤ[φ]` number `(4, 0)` for the unit squared, the unit `2` itself.
5. Non-forcing.  At `h ≠ 0`, two distinct declared units give two distinct generated
   paths with one hopping path and one family of induced sources
   (`two_units_two_paths`, `induced_sources_forget_unit`); the field-sector
   equations of the transported action are one and the same at every unit
   (`unit_unconstrained_by_field_sector`).

ROWS TOUCHED.  The coupled-action row (the seam-step class, the transport,
the time unit, the coupling shapes, and the joint action are declared
here); the physical spacetime attachment row (the rays are the declared
carrier map candidate; no attachment of a port to a spacetime point is
supplied); the source clock and duration row (the unit `τ` and the step
`h` are declared, none is selected); the light-signal row (no signal
propagation is attached to the timelike threshold); the laboratory clock
and energy calibration import (no unit, calibration, or readout is
attached to `τ`, `q`, or the clock action); the gravitation-route energy
identification (no identification of the clock action with an energy is
made).  The module discharges none of these rows.

NEGATIVES CITED.  The Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`): realized histories
select no velocity curvature or Legendre map, so the clock action shape,
the coupling shapes, and the transport are declared enrichments; cited at
scope only.

CONVENTIONS.  Signature `(+---)`; `Herm2 = ℝ × (Fin 3 → ℝ)`;
`lorentzQ v = v.1 ^ 2 - |v.2| ^ 2`; forward differences throughout.
Seam orientation from `seamLeft e` to `seamRight e`; a forward step across
`e` adds `seamVector e`, a backward step subtracts it; the hopping current
of a forward hop is `-(q / h)` on `e`.  `ℤ[φ]` pairs `(a, b)` denote
`a + bφ`; `seamVectorZ e = candidateRayZ (seamRight e) - candidateRayZ
(seamLeft e)` componentwise, with the positive sign.  The window of the
transported action is `N + 1` steps for both sectors.

FALSIFIER.  The module is wrong if the cumulative ray identity fails at
some admissible step, if some seam vector has `ℤ[φ]` norm squared other
than `(4, 0)`, if the endpoint transport constant differs from `12 - 4φ`,
or if the transported field equations differ from `monopole_field_equations`.

Axiom audit.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`; no `sorry`, no
`native_decide`, no project axiom.
-/

/-! ## Declared steps and seam-step worldlines -/

/-- A declared step: a rest, a forward crossing of a seam (from the smaller
to the larger endpoint), or a backward crossing. -/
inductive SeamStep
  | rest
  | forward (e : Fin 30)
  | backward (e : Fin 30)
  deriving DecidableEq

/-- The port reached from `u` by a step. -/
def stepTarget (u : Fin 12) : SeamStep → Fin 12
  | .rest => u
  | .forward e => seamRight e
  | .backward e => seamLeft e

/-- Admissibility of a step at a port: a forward crossing starts at the
smaller endpoint, a backward crossing at the larger one. -/
def StepAdmissible (u : Fin 12) : SeamStep → Prop
  | .rest => True
  | .forward e => u = seamLeft e
  | .backward e => u = seamRight e

/-- The port sequence generated by a start port and a step sequence. -/
def portSeq (start : Fin 12) (steps : ℕ → SeamStep) : ℕ → Fin 12
  | 0 => start
  | n + 1 => stepTarget (portSeq start steps n) (steps n)

/-- Declared seam-step worldline: a start port and a step sequence, every
step admissible at the port reached. -/
structure SeamStepWorldline where
  /-- The start port. -/
  start : Fin 12
  /-- The step sequence. -/
  steps : ℕ → SeamStep
  /-- Every step is admissible at the port reached. -/
  adm : ∀ n, StepAdmissible (portSeq start steps n) (steps n)

/-- The port reached at step `n`. -/
def SeamStepWorldline.port (w : SeamStepWorldline) : ℕ → Fin 12 :=
  portSeq w.start w.steps

theorem SeamStepWorldline.port_zero (w : SeamStepWorldline) : w.port 0 = w.start := rfl

theorem SeamStepWorldline.port_succ (w : SeamStepWorldline) (n : ℕ) :
    w.port (n + 1) = stepTarget (w.port n) (w.steps n) := rfl

/-- The exact `ℤ[φ]` vector of a step: zero at rest, the seam vector at a
forward crossing, its negative at a backward crossing. -/
def stepVectorZ : SeamStep → VecZ
  | .rest => 0
  | .forward e => seamVectorZ e
  | .backward e => vneg (seamVectorZ e)

/-- The exact seam vector is the ray of the larger endpoint minus the ray
of the smaller endpoint, with the positive sign; the identity is
definitional. -/
theorem seamVectorZ_eq_ray_difference (e : Fin 30) :
    seamVectorZ e = candidateRayZ (seamRight e) - candidateRayZ (seamLeft e) := by
  funext k
  rfl

theorem vneg_eq_neg (v : VecZ) : vneg v = -v := by
  funext k
  rfl

/-- The ray difference across an admissible step is the step vector. -/
theorem ray_step (u : Fin 12) (s : SeamStep) (hs : StepAdmissible u s) :
    candidateRayZ (stepTarget u s) - candidateRayZ u = stepVectorZ s := by
  cases s with
  | rest => exact sub_self _
  | forward e =>
    have h : u = seamLeft e := hs
    subst h
    exact (seamVectorZ_eq_ray_difference e).symm
  | backward e =>
    have h : u = seamRight e := hs
    subst h
    show candidateRayZ (seamLeft e) - candidateRayZ (seamRight e) = vneg (seamVectorZ e)
    rw [vneg_eq_neg, seamVectorZ_eq_ray_difference, neg_sub]

/-- An admissible step is a hop. -/
theorem hop_of_admissible (u : Fin 12) (s : SeamStep) (hs : StepAdmissible u s) :
    Hop u (stepTarget u s) := by
  cases s with
  | rest => exact Or.inl rfl
  | forward e =>
    have h : u = seamLeft e := hs
    exact Or.inr (Or.inl ⟨e, h.symm, rfl⟩)
  | backward e =>
    have h : u = seamRight e := hs
    exact Or.inr (Or.inr ⟨e, h.symm, rfl⟩)

/-- The transported hopping path of a seam-step worldline. -/
def hoppingPath (w : SeamStepWorldline) : HoppingPath where
  port := w.port
  hop := fun n ↦ hop_of_admissible _ _ (w.adm n)

theorem hoppingPath_port (w : SeamStepWorldline) (n : ℕ) :
    (hoppingPath w).port n = w.port n := rfl

/-- The exact cumulative spatial coordinate after `n` steps. -/
def spatialAtZ (w : SeamStepWorldline) (n : ℕ) : VecZ :=
  ∑ k ∈ Finset.range n, stepVectorZ (w.steps k)

theorem spatialAtZ_zero (w : SeamStepWorldline) : spatialAtZ w 0 = 0 := rfl

theorem spatialAtZ_succ (w : SeamStepWorldline) (n : ℕ) :
    spatialAtZ w (n + 1) = spatialAtZ w n + stepVectorZ (w.steps n) :=
  Finset.sum_range_succ _ _

/-- **Cumulative ray identity.**  The ray of the port reached at step `n`
minus the ray of the start port is the exact cumulative spatial coordinate,
at every step. -/
theorem ray_port_eq_spatialZ (w : SeamStepWorldline) (n : ℕ) :
    candidateRayZ (w.port n) - candidateRayZ w.start = spatialAtZ w n := by
  induction n with
  | zero => rw [w.port_zero, spatialAtZ_zero, sub_self]
  | succ n ih =>
    rw [spatialAtZ_succ, ← ih, w.port_succ, ← ray_step (w.port n) (w.steps n) (w.adm n)]
    abel

/-- The twelve candidate rays are pairwise distinct: the corpus theorem
`OPH.ScreenCarrierMapCandidate.candidateRayZ_injective`, restated in the
pointwise form used below. -/
theorem candidateRayZ_pairwise_distinct :
    ∀ i j : Fin 12, candidateRayZ i = candidateRayZ j → i = j :=
  fun _ _ h ↦ OPH.ScreenCarrierMapCandidate.candidateRayZ_injective h

/-- The port reached is the unique port whose ray difference from the start
ray is the cumulative spatial coordinate. -/
theorem port_determined_by_spatial (w : SeamStepWorldline) (n : ℕ) (p : Fin 12)
    (hp : candidateRayZ p - candidateRayZ w.start = spatialAtZ w n) : p = w.port n := by
  apply candidateRayZ_pairwise_distinct
  have h := ray_port_eq_spatialZ w n
  rw [← h] at hp
  exact sub_left_injective hp

/-! ## Real layer: the generated `Herm2` path -/

noncomputable section

/-- Evaluation of exact vectors as an additive homomorphism. -/
def evalVecHom : VecZ →+ Spatial where
  toFun := evalVec
  map_zero' := by
    funext k
    exact evalPhi_zero
  map_add' := fun u v ↦ by
    funext k
    exact evalPhi_add (u k) (v k)

theorem evalVecHom_apply (v : VecZ) : evalVecHom v = evalVec v := rfl

/-- The real step vector. -/
def stepVector (s : SeamStep) : Spatial := evalVec (stepVectorZ s)

theorem stepVector_rest : stepVector .rest = 0 := evalVecHom.map_zero

theorem stepVector_forward (e : Fin 30) : stepVector (.forward e) = seamVector e := rfl

theorem stepVector_backward (e : Fin 30) : stepVector (.backward e) = -seamVector e :=
  evalVec_vneg _

/-- The real cumulative spatial coordinate after `n` steps. -/
def spatialAt (w : SeamStepWorldline) (n : ℕ) : Spatial :=
  ∑ k ∈ Finset.range n, stepVector (w.steps k)

theorem spatialAt_eq_eval (w : SeamStepWorldline) (n : ℕ) :
    spatialAt w n = evalVec (spatialAtZ w n) := by
  unfold spatialAt spatialAtZ
  rw [← evalVecHom_apply, map_sum]
  rfl

theorem spatialAt_succ (w : SeamStepWorldline) (n : ℕ) :
    spatialAt w (n + 1) = spatialAt w n + stepVector (w.steps n) :=
  Finset.sum_range_succ _ _

/-- The generated `Herm2` path at the declared time unit `τ`: the time
coordinate advances by `τ` per step, the spatial coordinate is the
cumulative sum of the signed seam vectors.  The unit is declared. -/
def generatedPath (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) : Herm2 :=
  (τ * n, spatialAt w n)

theorem generatedPath_fst (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    (generatedPath τ w n).1 = τ * n := rfl

theorem generatedPath_snd (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    (generatedPath τ w n).2 = spatialAt w n := rfl

/-- The spatial increment of the generated path across step `n` is the
step vector, at every unit. -/
theorem generatedPath_increment (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    (generatedPath τ w (n + 1)).2 - (generatedPath τ w n).2 = stepVector (w.steps n) := by
  rw [generatedPath_snd, generatedPath_snd, spatialAt_succ, add_sub_cancel_left]

/-- **Real cumulative ray identity.**  The candidate ray of the port reached
minus the ray of the start port is the spatial coordinate of the generated
path, at every step and every unit. -/
theorem ray_port_eq_spatial (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) :
    candidateRay (w.port n) - candidateRay w.start = (generatedPath τ w n).2 := by
  rw [generatedPath_snd, spatialAt_eq_eval, ← ray_port_eq_spatialZ]
  unfold candidateRay
  rw [← evalVecHom_apply, ← evalVecHom_apply, ← evalVecHom_apply, map_sub]

theorem spatialDot_neg_right (u v : Spatial) : spatialDot u (-v) = -spatialDot u v := by
  unfold spatialDot
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  simp

/-! ## Source transport: the induced load at the seam endpoints -/

/-- The seam stepped at step `n`, when the step is a crossing. -/
def steppedSeam : SeamStep → Option (Fin 30)
  | .rest => none
  | .forward e => some e
  | .backward e => some e

/-- Forward step difference of the hopping load on the transported path at
a forward crossing: `q` at the larger endpoint, `-q` at the smaller one. -/
theorem hoppingLoad_difference_forward (q : ℝ) (w : SeamStepWorldline) (n : ℕ)
    (e : Fin 30) (hs : w.steps n = .forward e) (p : Fin 12) :
    (hoppingLoad q (hoppingPath w) (n + 1) - hoppingLoad q (hoppingPath w) n) p =
      q * (portIndicator (seamRight e) p - portIndicator (seamLeft e) p) := by
  have hl : w.port n = seamLeft e := by
    have h := w.adm n
    rw [hs] at h
    exact h
  have hr : w.port (n + 1) = seamRight e := by
    rw [w.port_succ, hs, hl]
    rfl
  simp only [Pi.sub_apply, hoppingLoad_apply, hoppingPath_port, hl, hr, portIndicator]
  split_ifs <;> ring

/-- Forward step difference of the hopping load at a backward crossing:
`-q` at the larger endpoint, `q` at the smaller one. -/
theorem hoppingLoad_difference_backward (q : ℝ) (w : SeamStepWorldline) (n : ℕ)
    (e : Fin 30) (hs : w.steps n = .backward e) (p : Fin 12) :
    (hoppingLoad q (hoppingPath w) (n + 1) - hoppingLoad q (hoppingPath w) n) p =
      -(q * (portIndicator (seamRight e) p - portIndicator (seamLeft e) p)) := by
  have hr : w.port n = seamRight e := by
    have h := w.adm n
    rw [hs] at h
    exact h
  have hl : w.port (n + 1) = seamLeft e := by
    rw [w.port_succ, hs, hr]
    rfl
  simp only [Pi.sub_apply, hoppingLoad_apply, hoppingPath_port, hl, hr, portIndicator]
  split_ifs <;> ring

/-- At a rest step the hopping load does not move. -/
theorem hoppingLoad_difference_rest (q : ℝ) (w : SeamStepWorldline) (n : ℕ)
    (hs : w.steps n = .rest) :
    hoppingLoad q (hoppingPath w) (n + 1) - hoppingLoad q (hoppingPath w) n = 0 := by
  have h : w.port (n + 1) = w.port n := by
    rw [w.port_succ, hs]
    rfl
  unfold hoppingLoad
  rw [hoppingPath_port, hoppingPath_port, h, sub_self]

/-- **Endpoint transport at every step.**  At both endpoints of the seam
stepped at step `n`, the E-paired induced load of the generated path equals
`-(κ / h) (12 - 4φ)` times the forward step difference of the unit hopping
load of the transported path.  This generalizes
`PortChargeMinimalCoupling.bridge_endpoints` from the one-step crossing to
every crossing step of every seam-step worldline; the agreement is
endpoint-local, as exhibited by `bridge_off_endpoint_exhibit`. -/
theorem transport_load_endpoints (κ h τ : ℝ) (w : SeamStepWorldline) (n : ℕ)
    (e : Fin 30) (hs : steppedSeam (w.steps n) = some e) (p : Fin 12)
    (hp : p = seamLeft e ∨ p = seamRight e) :
    inducedLoad κ h (generatedPath τ w) n p =
      -(κ / h) * seamCrossingConst *
        (hoppingLoad 1 (hoppingPath w) (n + 1) - hoppingLoad 1 (hoppingPath w) n) p := by
  have hlt := (seam_table_sound e).1
  rw [inducedLoad_eq_crossingFactor, generatedPath_increment]
  rcases hstep : w.steps n with _ | e' | e'
  · rw [hstep] at hs
    exact absurd hs (by simp [steppedSeam])
  · rw [hstep] at hs
    have he : e' = e := by simpa [steppedSeam] using hs
    subst he
    rw [stepVector_forward, hoppingLoad_difference_forward 1 w n e' hstep]
    unfold portIndicator
    rcases hp with rfl | rfl
    · rw [crossingFactor_left, if_neg hlt.ne, if_pos rfl]
      ring
    · rw [crossingFactor_right, if_pos rfl, if_neg hlt.ne']
      ring
  · rw [hstep] at hs
    have he : e' = e := by simpa [steppedSeam] using hs
    subst he
    rw [stepVector_backward, hoppingLoad_difference_backward 1 w n e' hstep]
    unfold crossingFactor portIndicator
    rw [spatialDot_neg_right]
    rcases hp with rfl | rfl
    · rw [ray_seamVector_left, if_neg hlt.ne, if_pos rfl, seamCrossingConst_eq]
      ring
    · rw [ray_seamVector_right, if_pos rfl, if_neg hlt.ne', seamCrossingConst_eq]
      ring

/-- At a rest step the induced load and the hopping step difference both
vanish at every port. -/
theorem transport_load_rest (κ h τ : ℝ) (w : SeamStepWorldline) (n : ℕ)
    (hs : w.steps n = .rest) (p : Fin 12) :
    inducedLoad κ h (generatedPath τ w) n p = 0 ∧
      (hoppingLoad 1 (hoppingPath w) (n + 1) - hoppingLoad 1 (hoppingPath w) n) p = 0 := by
  refine ⟨?_, by rw [hoppingLoad_difference_rest 1 w n hs]; rfl⟩
  apply inducedLoad_rest
  rw [generatedPath_snd, generatedPath_snd, spatialAt_succ, hs, stepVector_rest, add_zero]

/-- **Charge-fixed transport.**  At nonzero `h` and the charge-fixed normalization
`κ = -(q h) / (12 - 4φ)` the E-paired endpoint loads equal the hopping step
differences of charge `q` at every crossing step. -/
theorem transport_load_charge_fixed (q h τ : ℝ) (hh : h ≠ 0) (w : SeamStepWorldline)
    (n : ℕ) (e : Fin 30) (hs : steppedSeam (w.steps n) = some e) (p : Fin 12)
    (hp : p = seamLeft e ∨ p = seamRight e) :
    inducedLoad (-(q * h) / seamCrossingConst) h (generatedPath τ w) n p =
      (hoppingLoad q (hoppingPath w) (n + 1) - hoppingLoad q (hoppingPath w) n) p := by
  rw [transport_load_endpoints _ h τ w n e hs p hp]
  have hc := seamCrossingConst_ne_zero
  have hq : (hoppingLoad q (hoppingPath w) (n + 1) - hoppingLoad q (hoppingPath w) n) p =
      q * (hoppingLoad 1 (hoppingPath w) (n + 1) - hoppingLoad 1 (hoppingPath w) n) p := by
    simp only [Pi.sub_apply, hoppingLoad_apply]
    split_ifs <;> ring
  rw [hq]
  field_simp

/-! ## Current transport: the worldline seam current of a generated path -/

/-- The E-paired induced current at step `m` is `κ / h²` times the forward
step difference of the worldline seam current; the identity is
definitional. -/
theorem inducedCurrent_eq_step_difference (κ h : ℝ) (x : ℕ → Herm2) (m : ℕ) :
    inducedCurrent κ h x m =
      (κ / h ^ 2) • (worldlineSeamCurrent x (m + 1) - worldlineSeamCurrent x m) := rfl

end

/-- The exact `ℤ[φ]` pairing of the seam vector of `e'` with the vector of
a step. -/
def stepGramZ (s : SeamStep) (e' : Fin 30) : Zphi :=
  dotZ (seamVectorZ e') (stepVectorZ s)

/-- The `ℤ[φ]` norm squared of every seam vector is exactly `4`. -/
theorem stepNormSqZ_forward : ∀ e : Fin 30, stepGramZ (.forward e) e = ((4 : ℤ), (0 : ℤ)) := by
  decide

/-- The pairing of a seam vector with its negative is exactly `-4`. -/
theorem stepNormSqZ_backward :
    ∀ e : Fin 30, stepGramZ (.backward e) e = ((-4 : ℤ), (0 : ℤ)) := by
  decide

/-- The norm squared of the negated seam vector is exactly `4`. -/
theorem vneg_seamVectorZ_normSq :
    ∀ e : Fin 30, dotZ (vneg (seamVectorZ e)) (vneg (seamVectorZ e)) = ((4 : ℤ), (0 : ℤ)) := by
  decide

theorem stepGramZ_rest (e' : Fin 30) : stepGramZ .rest e' = 0 := dotZ_zero_right _

/-- The worldline seam current is supported off the stepped seam: the
pairing of the seam vectors of seams `1` and `0`, which share the port
`0`, is exactly `2`. -/
theorem seamGram_off_diagonal_exhibit :
    stepGramZ (.forward 0) 1 = ((2 : ℤ), (0 : ℤ)) ∧ (1 : Fin 30) ≠ 0 := by
  decide

noncomputable section

theorem evalPhi_int (a : ℤ) : evalPhi ((a, (0 : ℤ)) : Zphi) = (a : ℝ) := by
  show (a : ℝ) + ((0 : ℤ) : ℝ) * Real.goldenRatio = a
  push_cast
  ring

/-- **Worldline seam current of a generated path.**  At step `n` and seam
`e'` it is the exact pairing of the seam vector of `e'` with the step
vector, evaluated in `ℝ`; the unit `τ` does not enter. -/
theorem worldlineSeamCurrent_generated (τ : ℝ) (w : SeamStepWorldline) (n : ℕ)
    (e' : Fin 30) :
    worldlineSeamCurrent (generatedPath τ w) n e' = evalPhi (stepGramZ (w.steps n) e') := by
  unfold worldlineSeamCurrent stepGramZ
  rw [generatedPath_increment, evalPhi_dotZ]
  rfl

/-- The real seam-vector norm squared, `4`. -/
def seamNormSq : ℝ := 4

theorem seamNormSq_eq_eval : seamNormSq = evalPhi ((4 : ℤ), (0 : ℤ)) := by
  rw [evalPhi_int]
  rfl

/-- At the stepped seam the worldline seam current is `4` at a forward
crossing and `-4` at a backward crossing. -/
theorem worldlineSeamCurrent_at_step (τ : ℝ) (w : SeamStepWorldline) (n : ℕ) (e : Fin 30) :
    (w.steps n = .forward e → worldlineSeamCurrent (generatedPath τ w) n e = seamNormSq) ∧
    (w.steps n = .backward e → worldlineSeamCurrent (generatedPath τ w) n e = -seamNormSq) := by
  constructor
  · intro hs
    rw [worldlineSeamCurrent_generated, hs, stepNormSqZ_forward, evalPhi_int]
    rfl
  · intro hs
    rw [worldlineSeamCurrent_generated, hs, stepNormSqZ_backward, evalPhi_int]
    show ((-4 : ℤ) : ℝ) = -(4 : ℝ)
    push_cast
    ring

/-- At a forward crossing of seam `0` the worldline seam current at seam `1`
is `2`, while the hopping current at seam `1` vanishes: the worldline seam
current is not supported on the single seam stepped. -/
theorem worldlineSeamCurrent_off_seam_exhibit (q h τ : ℝ) (w : SeamStepWorldline) (n : ℕ)
    (hs : w.steps n = .forward 0) :
    worldlineSeamCurrent (generatedPath τ w) n 1 = 2 ∧
      hoppingCurrent q h (hoppingPath w) n 1 = 0 := by
  constructor
  · rw [worldlineSeamCurrent_generated, hs, seamGram_off_diagonal_exhibit.1, evalPhi_int]
    rfl
  · have hl : w.port n = seamLeft 0 := by
      have h := w.adm n
      rw [hs] at h
      exact h
    have hr : w.port (n + 1) = seamRight 0 := by
      rw [w.port_succ, hs, hl]
      rfl
    rw [hoppingCurrent_forward q h (hoppingPath w) n 0 hl hr]
    simp

/-- Projection of a seam field onto one seam, divided by the seam-vector
norm squared. -/
def seamProjection (J : Fin 30 → ℝ) (e : Fin 30) : Fin 30 → ℝ :=
  fun e' ↦ if e' = e then J e / seamNormSq else 0

/-- **Hopping current as a projection.**  At every crossing step the hopping
current of the transported path is `-(q / h)` times the projection of the
worldline seam current onto the stepped seam divided by the seam norm. -/
theorem hoppingCurrent_eq_projection (q h τ : ℝ) (w : SeamStepWorldline) (n : ℕ)
    (e : Fin 30) (hs : steppedSeam (w.steps n) = some e) :
    hoppingCurrent q h (hoppingPath w) n =
      -((q / h) • seamProjection (worldlineSeamCurrent (generatedPath τ w) n) e) := by
  rcases hstep : w.steps n with _ | e' | e'
  · rw [hstep] at hs
    exact absurd hs (by simp [steppedSeam])
  · rw [hstep] at hs
    have he : e' = e := by simpa [steppedSeam] using hs
    subst he
    have hl : w.port n = seamLeft e' := by
      have h := w.adm n
      rw [hstep] at h
      exact h
    have hr : w.port (n + 1) = seamRight e' := by
      rw [w.port_succ, hstep, hl]
      rfl
    rw [hoppingCurrent_forward q h (hoppingPath w) n e' hl hr]
    funext x
    unfold seamProjection
    rw [(worldlineSeamCurrent_at_step τ w n e').1 hstep]
    unfold seamNormSq
    simp only [Pi.neg_apply, Pi.smul_apply, smul_eq_mul]
    split_ifs <;> ring
  · rw [hstep] at hs
    have he : e' = e := by simpa [steppedSeam] using hs
    subst he
    have hr : w.port n = seamRight e' := by
      have h := w.adm n
      rw [hstep] at h
      exact h
    have hl : w.port (n + 1) = seamLeft e' := by
      rw [w.port_succ, hstep, hr]
      rfl
    rw [hoppingCurrent_backward q h (hoppingPath w) n e' hr hl]
    funext x
    unfold seamProjection
    rw [(worldlineSeamCurrent_at_step τ w n e').2 hstep]
    unfold seamNormSq
    simp only [Pi.neg_apply, Pi.smul_apply, smul_eq_mul]
    split_ifs <;> ring

/-! ## Joint action: transported monopole route plus the committed clock action -/

/-- The declared transported action: the monopole coupled action of the
transported hopping path plus the committed clock action of the generated
path, both over the window of `N + 1` steps. -/
def transportedAction (q h τ : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) : ℝ :=
  monopoleCoupledAction q h N A φ ρ J (hoppingPath w) +
    clockAction (N + 1) (generatedPath τ w)

/-- **Field-sector stationarity of the transported action** at `h ≠ 0` is the committed
scaled Ampere update at the augmented current and the committed Gauss
constraint at the augmented load, with the hopping sources of the
transported path; the clock term does not enter. -/
theorem transported_field_equations (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) :
    ((∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
        transportedAction q h τ N (A + a) φ ρ J w =
          transportedAction q h τ N A φ ρ J w + quadraticRemainder h N a 0) ∧
      (∀ f : ℕ → Fin 12 → ℝ,
        transportedAction q h τ N A (φ + f) ρ J w =
          transportedAction q h τ N A φ ρ J w + quadraticRemainder h N 0 f)) ↔
    ((∀ m, m < N →
        electricFieldScaled h A φ (m + 1) - electricFieldScaled h A φ m =
          h • (localMaxwellOperator (A (m + 1)) -
            (J m + hoppingCurrent q h (hoppingPath w) m))) ∧
      (∀ n, n < N + 1 →
        realBoundary (electricFieldScaled h A φ n) =
          ρ n + hoppingLoad q (hoppingPath w) n)) := by
  rw [← monopole_field_equations q h hh N A φ ρ J (hoppingPath w)]
  unfold transportedAction
  refine and_congr (forall_congr' fun a ↦ imp_congr Iff.rfl (imp_congr Iff.rfl ?_))
    (forall_congr' fun f ↦ ?_)
  · rw [add_right_comm _ (clockAction (N + 1) (generatedPath τ w)) (quadraticRemainder h N a 0)]
    exact add_left_inj _
  · rw [add_right_comm _ (clockAction (N + 1) (generatedPath τ w)) (quadraticRemainder h N 0 f)]
    exact add_left_inj _

/-! ## The clock action of a generated path -/

/-- The real norm squared of a step vector: `0` at rest, `4` at a crossing. -/
def stepNormSq : SeamStep → ℝ
  | .rest => 0
  | .forward _ => seamNormSq
  | .backward _ => seamNormSq

theorem spatialNormSq_eq_dot (v : Spatial) : spatialNormSq v = spatialDot v v := by
  unfold spatialNormSq spatialDot
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  ring

theorem spatialNormSq_stepVector (s : SeamStep) :
    spatialNormSq (stepVector s) = stepNormSq s := by
  rw [spatialNormSq_eq_dot]
  unfold stepVector
  rw [← evalPhi_dotZ]
  cases s with
  | rest =>
    show evalPhi (dotZ 0 0) = 0
    rw [dotZ_zero_right, evalPhi_zero]
  | forward e =>
    show evalPhi (stepGramZ (.forward e) e) = seamNormSq
    rw [stepNormSqZ_forward, evalPhi_int]
    rfl
  | backward e =>
    show evalPhi (dotZ (vneg (seamVectorZ e)) (vneg (seamVectorZ e))) = seamNormSq
    rw [vneg_seamVectorZ_normSq, evalPhi_int]
    rfl

/-- The Lorentz square of one step of the generated path: the unit squared
minus the step norm squared. -/
theorem lorentzQ_generated_step (τ : ℝ) (w : SeamStepWorldline) (k : ℕ) :
    lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k) =
      τ ^ 2 - stepNormSq (w.steps k) := by
  unfold lorentzQ
  rw [Prod.fst_sub, Prod.snd_sub, generatedPath_fst, generatedPath_fst,
    generatedPath_increment, spatialNormSq_stepVector]
  push_cast
  ring

/-- **Clock action of a generated path.**  Over `M` steps it is the sum of
`τ² - 4` at crossing steps and `τ²` at rest steps. -/
theorem clockAction_generated (τ : ℝ) (w : SeamStepWorldline) (M : ℕ) :
    clockAction M (generatedPath τ w) =
      ∑ k ∈ Finset.range M, (τ ^ 2 - stepNormSq (w.steps k)) := by
  unfold clockAction
  exact Finset.sum_congr rfl fun k _ ↦ lorentzQ_generated_step τ w k

theorem stepNormSq_of_ne_rest (s : SeamStep) (hs : s ≠ .rest) : stepNormSq s = 4 := by
  cases s with
  | rest => exact absurd rfl hs
  | forward _ => rfl
  | backward _ => rfl

/-- With no rest step in the window the clock action is the step count times
the Lorentz square of one seam step, `M (τ² - 4)`. -/
theorem clockAction_generated_no_rest (τ : ℝ) (w : SeamStepWorldline) (M : ℕ)
    (hM : ∀ k, k < M → w.steps k ≠ .rest) :
    clockAction M (generatedPath τ w) = (M : ℝ) * (τ ^ 2 - 4) := by
  rw [clockAction_generated]
  rw [Finset.sum_congr rfl fun k hk ↦ by
    rw [stepNormSq_of_ne_rest _ (hM k (Finset.mem_range.mp hk))]]
  rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]

/-- The exact threshold for the unit squared: the seam-vector norm squared
`(4, 0)` in `ℤ[φ]`, the unit `2` itself. -/
def unitSquaredThresholdZ : Zphi := ((4 : ℤ), (0 : ℤ))

theorem unitSquaredThresholdZ_eq_normSq (e : Fin 30) :
    unitSquaredThresholdZ = stepGramZ (.forward e) e := (stepNormSqZ_forward e).symm

theorem evalPhi_threshold : evalPhi unitSquaredThresholdZ = (2 : ℝ) ^ 2 := by
  unfold unitSquaredThresholdZ
  rw [evalPhi_int]
  norm_num

/-- **Timelike threshold.**  A crossing step of the generated path is
timelike (positive Lorentz square) exactly when `4 < τ²`; at `τ² = 4` it
is null, below it is spacelike.  Which regime holds is fixed by the declared
unit alone. -/
theorem seam_step_timelike_iff (τ : ℝ) (w : SeamStepWorldline) (k : ℕ)
    (hk : w.steps k ≠ .rest) :
    0 < lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k) ↔
      evalPhi unitSquaredThresholdZ < τ ^ 2 := by
  rw [lorentzQ_generated_step, stepNormSq_of_ne_rest _ hk, evalPhi_threshold]
  constructor <;> intro h <;> linarith

theorem seam_step_timelike_of_two_lt (τ : ℝ) (hτ : 2 < τ) (w : SeamStepWorldline) (k : ℕ)
    (hk : w.steps k ≠ .rest) :
    0 < lorentzQ (generatedPath τ w (k + 1) - generatedPath τ w k) := by
  rw [seam_step_timelike_iff τ w k hk, evalPhi_threshold]
  nlinarith

/-! ## Non-forcing of the unit -/

/-- **Two units, two paths, one transport.**  Distinct declared units give
distinct generated paths whose spatial coordinates agree at every step; the
transported hopping path does not see the unit. -/
theorem two_units_two_paths (τ₁ τ₂ : ℝ) (hτ : τ₁ ≠ τ₂) (w : SeamStepWorldline) :
    generatedPath τ₁ w ≠ generatedPath τ₂ w ∧
      ∀ n, (generatedPath τ₁ w n).2 = (generatedPath τ₂ w n).2 := by
  refine ⟨fun heq ↦ hτ ?_, fun n ↦ rfl⟩
  have h := congrArg (fun x ↦ (x 1).1) heq
  simp only [generatedPath_fst, Nat.cast_one, mul_one] at h
  exact h

/-- The E-paired induced sources of a generated path do not see the unit:
the worldline seam current, the induced load, and the induced current
coincide at any two units. -/
theorem induced_sources_forget_unit (κ h τ₁ τ₂ : ℝ) (w : SeamStepWorldline) :
    worldlineSeamCurrent (generatedPath τ₁ w) = worldlineSeamCurrent (generatedPath τ₂ w) ∧
    inducedLoad κ h (generatedPath τ₁ w) = inducedLoad κ h (generatedPath τ₂ w) ∧
    inducedCurrent κ h (generatedPath τ₁ w) = inducedCurrent κ h (generatedPath τ₂ w) := by
  have hw : worldlineSeamCurrent (generatedPath τ₁ w) =
      worldlineSeamCurrent (generatedPath τ₂ w) := by
    funext n e
    unfold worldlineSeamCurrent
    rw [generatedPath_increment, generatedPath_increment]
  refine ⟨hw, ?_, ?_⟩
  · funext n
    unfold inducedLoad
    rw [hw]
  · funext m
    unfold inducedCurrent
    rw [hw]

/-- **The unit is unconstrained by the field sector.**  At `h ≠ 0`, at every declared
unit the field-sector stationarity of the transported action is the same
pair of committed equations at the transported hopping sources, which do
not contain the unit. -/
theorem unit_unconstrained_by_field_sector (q h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) :
    ∀ τ : ℝ,
    ((∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
        transportedAction q h τ N (A + a) φ ρ J w =
          transportedAction q h τ N A φ ρ J w + quadraticRemainder h N a 0) ∧
      (∀ f : ℕ → Fin 12 → ℝ,
        transportedAction q h τ N A (φ + f) ρ J w =
          transportedAction q h τ N A φ ρ J w + quadraticRemainder h N 0 f)) ↔
    ((∀ m, m < N →
        electricFieldScaled h A φ (m + 1) - electricFieldScaled h A φ m =
          h • (localMaxwellOperator (A (m + 1)) -
            (J m + hoppingCurrent q h (hoppingPath w) m))) ∧
      (∀ n, n < N + 1 →
        realBoundary (electricFieldScaled h A φ n) =
          ρ n + hoppingLoad q (hoppingPath w) n)) :=
  fun τ ↦ transported_field_equations q h τ hh N A φ ρ J w

/-! ## The one-step crossing as a seam-step worldline -/

/-- The seam-step worldline crossing seam `e` forward at step `0` and
resting afterwards. -/
def crossingWorldline (e : Fin 30) : SeamStepWorldline where
  start := seamLeft e
  steps := fun n ↦ if n = 0 then .forward e else .rest
  adm := by
    intro n
    cases n with
    | zero => exact rfl
    | succ n => exact trivial

theorem crossingWorldline_port (e : Fin 30) (n : ℕ) :
    (crossingWorldline e).port n = (crossingPath e).port n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    rw [SeamStepWorldline.port_succ, ih]
    cases n with
    | zero => rfl
    | succ n => rfl

theorem crossingWorldline_spatial (e : Fin 30) (n : ℕ) :
    spatialAt (crossingWorldline e) n = if n = 0 then 0 else seamVector e := by
  induction n with
  | zero => rfl
  | succ n ih =>
    rw [spatialAt_succ, ih]
    cases n with
    | zero =>
      show (0 : Spatial) + stepVector (.forward e) = seamVector e
      rw [stepVector_forward, zero_add]
    | succ n =>
      show seamVector e + stepVector .rest = seamVector e
      rw [stepVector_rest, add_zero]

/-- At unit `0` the generated path of the crossing worldline is the
seam-crossing worldline of `ChargeFixedInteraction`, and its transported
hopping path has the ports of the crossing path of
`PortChargeMinimalCoupling`: `bridge_endpoints` is the step-`0` case of
`transport_load_endpoints`. -/
theorem crossingWorldline_generated (e : Fin 30) :
    generatedPath 0 (crossingWorldline e) = seamCrossingWorldline e ∧
      (hoppingPath (crossingWorldline e)).port = (crossingPath e).port := by
  refine ⟨funext fun n ↦ Prod.ext ?_ ?_, funext fun n ↦ crossingWorldline_port e n⟩
  · show (0 : ℝ) * n = 0
    ring
  · exact crossingWorldline_spatial e n

end

end OPH.WorldlineHopTransport

/- Axiom audit: standard axioms only (`propext`, `Classical.choice`,
`Quot.sound`); no `sorry`, no `native_decide`, no project axiom. -/

#print axioms OPH.WorldlineHopTransport.seamVectorZ_eq_ray_difference
#print axioms OPH.WorldlineHopTransport.hop_of_admissible
#print axioms OPH.WorldlineHopTransport.ray_port_eq_spatialZ
#print axioms OPH.WorldlineHopTransport.ray_port_eq_spatial
#print axioms OPH.WorldlineHopTransport.candidateRayZ_pairwise_distinct
#print axioms OPH.WorldlineHopTransport.port_determined_by_spatial
#print axioms OPH.WorldlineHopTransport.transport_load_endpoints
#print axioms OPH.WorldlineHopTransport.transport_load_rest
#print axioms OPH.WorldlineHopTransport.transport_load_charge_fixed
#print axioms OPH.WorldlineHopTransport.inducedCurrent_eq_step_difference
#print axioms OPH.WorldlineHopTransport.stepNormSqZ_forward
#print axioms OPH.WorldlineHopTransport.seamGram_off_diagonal_exhibit
#print axioms OPH.WorldlineHopTransport.worldlineSeamCurrent_generated
#print axioms OPH.WorldlineHopTransport.worldlineSeamCurrent_at_step
#print axioms OPH.WorldlineHopTransport.worldlineSeamCurrent_off_seam_exhibit
#print axioms OPH.WorldlineHopTransport.hoppingCurrent_eq_projection
#print axioms OPH.WorldlineHopTransport.transported_field_equations
#print axioms OPH.WorldlineHopTransport.clockAction_generated
#print axioms OPH.WorldlineHopTransport.clockAction_generated_no_rest
#print axioms OPH.WorldlineHopTransport.seam_step_timelike_iff
#print axioms OPH.WorldlineHopTransport.seam_step_timelike_of_two_lt
#print axioms OPH.WorldlineHopTransport.two_units_two_paths
#print axioms OPH.WorldlineHopTransport.induced_sources_forget_unit
#print axioms OPH.WorldlineHopTransport.unit_unconstrained_by_field_sector
#print axioms OPH.WorldlineHopTransport.crossingWorldline_generated
