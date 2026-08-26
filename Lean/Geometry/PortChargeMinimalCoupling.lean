import Geometry.ChargeFixedInteraction

set_option autoImplicit false

open scoped BigOperators

namespace OPH.PortChargeMinimalCoupling

open OPH.ScreenCarrierMapCandidate
open OPH.SeamCurrentCarrierQuotient
open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.CarrierDynamicsCompatibility
open OPH.ChargeFixedInteraction
open OPH.C1Lorentz (Spatial Herm2 spatialDot lorentzB lorentzQ)

/-!
# Minimal coupling of a declared hopping port charge on the committed carrier

STATUS.  Candidate module on the coupled-action row.  The committed
E-paired interaction candidate (`Geometry/ChargeFixedInteraction.lean`)
induces polarization-type sources and carries no monopole charge.  This
module builds a minimal-coupling alternative of the shape named there, with
the `Herm2` worldline increment replaced by a declared port hop: a declared
point charge hopping along the ports of the thirty-seam complex, minimally
coupled to the seam
potential and the port potential with the source weights of the committed
window action.  Every statement is an exact identity on the finite complex.

WHAT IS PROVED.

1. Declared point charge.  A hopping path is a port path whose consecutive
   ports are equal or joined by a seam (`HoppingPath`, proof field
   `hop`).  Its induced seam current at step `n` is `-(q / h)` times the
   oriented indicator of the seam hopped (`hoppingCurrent`), zero at a
   resting step (`hoppingCurrent_rest`); its induced port load is `q`
   times the indicator of the occupied port (`hoppingLoad`).  The sources
   satisfy the committed scaled continuity equation
   `ρ (n+1) - ρ n + h • ∂ (J n) = 0` identically (`hopping_continuity`),
   and the total load is `q` at every step (`hoppingLoad_total`).
2. Minimal coupling.  The interaction is the window sum of the seam pairing
   of the forward seam potential against the hopping current plus the port
   pairing of the port potential against the hopping load, with the weight
   `h` of the committed window action (`interactionA`).  The coupled
   action equals the committed window action at the augmented sources
   exactly, with no endpoint term (`monopoleCoupledAction_eq_augmented`),
   so its field equations are the committed scaled Ampere update and the
   committed Gauss constraint at the augmented sources
   (`monopole_stationary_A_iff_ampere`, `monopole_stationary_phi_iff_gauss`,
   `monopole_field_equations`).
3. Gauge change and conservation.  For a general seam current and port
   load, the source pairing moves under the committed scaled gauge
   transformation by the window sum of the gauge function at the forward
   step paired against the continuity residual, plus the endpoint
   difference of the gauge function paired against the load
   (`sourcePairing_gauge`).  Invariance under every gauge function
   vanishing at both window endpoints holds exactly when the continuity
   residual vanishes at every interior step
   (`sourcePairing_gauge_invariant_iff`).  For the hopping charge the
   change is `q` times the difference of the gauge function at the
   occupied ports at the two window endpoints (`interactionA_gauge`), and
   vanishes for endpoint-vanishing gauge functions
   (`interactionA_gauge_invariant`).
4. Field energy balance at a hop.  Along the scaled Ampere evolution
   driven by the hopping current the scaled staggered energy moves per
   step by
   `(q / 2)` times the pairing of `E n + E (n+1)` against the oriented
   hop indicator (`hopping_work_energy`): zero at a resting step
   (`hopping_work_energy_rest`), `(q / 2) (E n e + E (n+1) e)` for a
   forward hop across seam `e` (`hopping_work_energy_forward`), the
   negative of that for a backward hop (`hopping_work_energy_backward`).
   No energy of the charge and no equation of motion for the path are
   defined; the identity is the committed field energy balance specialized
   to the hopping current.
5. Bridge to the E-paired candidate.  On the seam-crossing worldline of
   seam `e` and the crossing path of `e`, the E-paired induced load at the
   two endpoints of `e` equals `-(κ / h) (12 - 4φ)` times the forward step
   difference of the unit hopping load (`bridge_endpoints`).  At the
   charge-fixed normalization `κ = -(q h) / (12 - 4φ)` the E-paired load
   equals the step difference of the hopping load of charge `q` at both
   endpoints (`bridge_endpoints_charge_fixed`), and endpoint agreement at
   the larger endpoint holds exactly at that normalization
   (`bridge_iff_charge_fixed`).  The agreement is endpoint-local: at seam
   `0` the E-paired load at port `10` equals `(κ / h) (12 - 4φ)`, nonzero
   at nonzero `κ`, while the hopping step difference at port `10` is zero
   (`bridge_off_endpoint_exhibit`).  So at the exhibited seam the
   polarization load is the step difference of the monopole load at the
   seam endpoints and differs from it away from them.
6. Non-forcing.  Every charge value satisfies the continuity equation and
   the endpoint-vanishing gauge invariance
   (`charge_unconstrained_by_field_sector`); two distinct charges give two
   distinct hopping loads and two distinct total loads
   (`two_charges_two_sources`).

ROWS TOUCHED.  The coupled-action row (the hopping path, the charge, the
coupling shape, and the source weights are declared here); the physical
spacetime attachment row (no attachment of a port to a spacetime point is
supplied); the source clock and duration row (the step `h` is declared,
the current scales as `q / h`, no step is selected); the laboratory clock
and energy calibration import (no unit, no calibration, no readout is
attached to `q` or to the energy identity); the gravitation-route energy
identification (the scaled staggered form is the committed field quadratic
form, and no identification with a gravitational or inertial energy is
made).  The module discharges none of these rows.

NEGATIVES CITED.  The Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`): realized source
histories select no velocity curvature or Legendre map, so the coupling
shape and the window action shape are declared enrichments.  The abstract
rate non-identifiability: its repair-layer bridge fails per
`RateBridgeObstruction`, and it forbids nothing about a source-hosted
process; nothing here cites it for or against a clock.

BOUNDARY.  The hopping path is declared and is not derived from the
`Herm2` worldline of the committed clock action; no variational equation
for the path exists in this module, since the path ranges over a discrete
set of port sequences and no stationarity clause is defined for it.  The
physical identification of a port charge lives on the coupled-action row
and the physical spacetime attachment row.

CONVENTIONS.  Signature `(+---)`; `Herm2 = ℝ × (Fin 3 → ℝ)`;
`lorentzQ v = v.1 ^ 2 - |v.2| ^ 2`; forward differences throughout.
Seam orientation: from the smaller endpoint `seamLeft e` to the larger
endpoint `seamRight e`; a forward hop across `e` carries indicator `+1`
on `e`, a backward hop `-1`.  Port coboundary `d φ e = φ (right) - φ (left)`;
port boundary `∂` is its transpose, so `∂ c` at a port is the net signed
inflow.  Scaled electric field `E n = -(A (n+1) - A n) / h - d (φ n)`.
Committed Gauss constraint `∂ E = ρ`; committed Ampere update
`E (m+1) - E m = h • (Cᵀ B (m+1) - J m)`; committed continuity
`ρ (n+1) - ρ n + h • ∂ (J n) = 0`.  Committed source weights of the window
action: `+ h ⟨J n, A (n+1)⟩ + h ⟨ρ n, φ n⟩`.  Scaled gauge transformation
`A n ↦ A n + d (χ n)`, `φ n ↦ φ n - (χ (n+1) - χ n) / h`.  Because `∂`
is the net inflow, the analogue of minus the divergence, a charge `q`
moving along the seam orientation carries current `-(q / h)` on that seam
in the committed continuity equation; the factor `1 / h` makes the
equation hold at every step, and at `h = 1` the current is `-q` times the
oriented indicator.  The sign of the energy identity in item 4 is stated
relative to this convention.

FALSIFIER.  The module is wrong if the hopping sources violate the
committed continuity equation, if the coupled action differs from the
window action at the augmented sources, if the gauge change misses a
term, if the energy identity disagrees with `energy_balance_scaled`, or
if the endpoint pairing constant differs from `12 - 4φ`.

Axiom audit.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`; no `sorry`, no
`native_decide`, no project axiom.
-/

/-! ## Declared hopping path -/

/-- Declared hop relation: two ports are equal or joined by a seam in
either direction. -/
def Hop (u v : Fin 12) : Prop :=
  u = v ∨ (∃ e : Fin 30, seamLeft e = u ∧ seamRight e = v) ∨
    (∃ e : Fin 30, seamRight e = u ∧ seamLeft e = v)

/-- Declared hopping path: a port sequence whose consecutive ports hop. -/
structure HoppingPath where
  /-- The occupied port at each step. -/
  port : ℕ → Fin 12
  /-- Consecutive ports are equal or joined by a seam. -/
  hop : ∀ n, Hop (port n) (port (n + 1))

/-- Oriented seam indicator of a hop from `u` to `v`: `+1` on the seam
oriented `u → v`, `-1` on the seam oriented `v → u`, `0` elsewhere. -/
def hopIndicator (u v : Fin 12) : Fin 30 → ℝ :=
  fun e ↦ if seamLeft e = u ∧ seamRight e = v then 1
    else if seamRight e = u ∧ seamLeft e = v then -1 else 0

/-- Port indicator of the occupied port. -/
def portIndicator (u : Fin 12) : Fin 12 → ℝ :=
  fun p ↦ if p = u then 1 else 0

theorem hopIndicator_self (u : Fin 12) : hopIndicator u u = 0 := by
  funext e
  have hlt := (seam_table_sound e).1
  unfold hopIndicator
  split_ifs with h1 h2
  · exact absurd (h1.1.trans h1.2.symm) hlt.ne
  · exact absurd (h2.1.trans h2.2.symm) hlt.ne'
  · rfl

theorem hopIndicator_forward (e₀ : Fin 30) :
    hopIndicator (seamLeft e₀) (seamRight e₀) =
      fun e ↦ if e = e₀ then (1 : ℝ) else 0 := by
  funext e
  have hlt := (seam_table_sound e).1
  have hlt₀ := (seam_table_sound e₀).1
  unfold hopIndicator
  by_cases he : e = e₀
  · subst he
    simp
  · rw [if_neg he]
    split_ifs with h1 h2
    · exact absurd (seam_table_injective (Prod.ext h1.1 h1.2)) he
    · rw [← h2.1, ← h2.2] at hlt₀
      exact absurd hlt (lt_asymm hlt₀)
    · rfl

theorem hopIndicator_backward (e₀ : Fin 30) :
    hopIndicator (seamRight e₀) (seamLeft e₀) =
      fun e ↦ if e = e₀ then (-1 : ℝ) else 0 := by
  funext e
  have hlt := (seam_table_sound e).1
  have hlt₀ := (seam_table_sound e₀).1
  unfold hopIndicator
  by_cases he : e = e₀
  · subst he
    simp [hlt.ne']
  · rw [if_neg he]
    split_ifs with h1 h2
    · rw [← h1.1, ← h1.2] at hlt₀
      exact absurd hlt (lt_asymm hlt₀)
    · exact absurd (seam_table_injective (Prod.ext h2.2 h2.1)) he
    · rfl

/-- Port boundary of a single-seam field. -/
theorem realBoundary_single (e₀ : Fin 30) (c : ℝ) (p : Fin 12) :
    realBoundary (fun e ↦ if e = e₀ then c else 0) p =
      (if p = seamRight e₀ then c else 0) - (if p = seamLeft e₀ then c else 0) := by
  rw [realBoundary_apply]
  have h1 : ∀ e : Fin 30, (if p = seamRight e then (if e = e₀ then c else 0) else 0) =
      if e = e₀ then (if p = seamRight e₀ then c else 0) else 0 := by
    intro e
    by_cases he : e = e₀
    · subst he; simp
    · simp [he]
  have h2 : ∀ e : Fin 30, (if p = seamLeft e then (if e = e₀ then c else 0) else 0) =
      if e = e₀ then (if p = seamLeft e₀ then c else 0) else 0 := by
    intro e
    by_cases he : e = e₀
    · subst he; simp
    · simp [he]
  rw [Finset.sum_sub_distrib, Finset.sum_congr rfl fun e _ ↦ h1 e,
    Finset.sum_congr rfl fun e _ ↦ h2 e, Finset.sum_ite_eq', Finset.sum_ite_eq']
  simp

/-- The port boundary of a hop indicator is the difference of the port
indicators of the target and the source port. -/
theorem realBoundary_hopIndicator (u v : Fin 12) (huv : Hop u v) (p : Fin 12) :
    realBoundary (hopIndicator u v) p = portIndicator v p - portIndicator u p := by
  unfold portIndicator
  rcases huv with rfl | ⟨e₀, rfl, rfl⟩ | ⟨e₀, rfl, rfl⟩
  · rw [hopIndicator_self, map_zero, Pi.zero_apply, sub_self]
  · rw [hopIndicator_forward, realBoundary_single]
  · rw [hopIndicator_backward, realBoundary_single]
    split_ifs <;> ring

noncomputable section

/-! ## The hopping sources -/

/-- Declared hopping port load: `q` times the indicator of the occupied
port. -/
def hoppingLoad (q : ℝ) (γ : HoppingPath) (n : ℕ) : Fin 12 → ℝ :=
  q • portIndicator (γ.port n)

/-- Declared hopping seam current: `-(q / h)` times the oriented indicator
of the seam hopped across step `n`.  The sign is the one under which the
committed scaled continuity equation `ρ (n+1) - ρ n + h • ∂ (J n) = 0`
holds, since the committed port boundary `∂` is the net signed inflow
(the analogue of minus the divergence); the factor `1 / h` makes the
equation hold at every step.  At `h = 1` it is `-q` times the
indicator. -/
def hoppingCurrent (q h : ℝ) (γ : HoppingPath) (n : ℕ) : Fin 30 → ℝ :=
  -((q / h) • hopIndicator (γ.port n) (γ.port (n + 1)))

theorem hoppingLoad_apply (q : ℝ) (γ : HoppingPath) (n : ℕ) (p : Fin 12) :
    hoppingLoad q γ n p = if p = γ.port n then q else 0 := by
  unfold hoppingLoad portIndicator
  simp only [Pi.smul_apply, smul_eq_mul]
  split_ifs <;> ring

/-- Total load `q` at every step: the hopping charge is a monopole. -/
theorem hoppingLoad_total (q : ℝ) (γ : HoppingPath) (n : ℕ) :
    (∑ p : Fin 12, hoppingLoad q γ n p) = q := by
  simp only [hoppingLoad_apply, Finset.sum_ite_eq', Finset.mem_univ, if_true]

/-- A resting step induces no current. -/
theorem hoppingCurrent_rest (q h : ℝ) (γ : HoppingPath) (n : ℕ)
    (hn : γ.port (n + 1) = γ.port n) : hoppingCurrent q h γ n = 0 := by
  unfold hoppingCurrent
  rw [hn, hopIndicator_self, smul_zero, neg_zero]

theorem hoppingCurrent_forward (q h : ℝ) (γ : HoppingPath) (n : ℕ) (e : Fin 30)
    (hl : γ.port n = seamLeft e) (hr : γ.port (n + 1) = seamRight e) :
    hoppingCurrent q h γ n = fun e' ↦ if e' = e then -(q / h) else 0 := by
  unfold hoppingCurrent
  rw [hl, hr, hopIndicator_forward]
  funext e'
  simp only [Pi.neg_apply, Pi.smul_apply, smul_eq_mul]
  split_ifs <;> ring

theorem hoppingCurrent_backward (q h : ℝ) (γ : HoppingPath) (n : ℕ) (e : Fin 30)
    (hl : γ.port n = seamRight e) (hr : γ.port (n + 1) = seamLeft e) :
    hoppingCurrent q h γ n = fun e' ↦ if e' = e then q / h else 0 := by
  unfold hoppingCurrent
  rw [hl, hr, hopIndicator_backward]
  funext e'
  simp only [Pi.neg_apply, Pi.smul_apply, smul_eq_mul]
  split_ifs <;> ring

/-- The port boundary of the hopping current: `-(q / h)` times the
difference of the port indicators of the next and the current port. -/
theorem realBoundary_hoppingCurrent (q h : ℝ) (γ : HoppingPath) (n : ℕ) :
    realBoundary (hoppingCurrent q h γ n) =
      -((q / h) • (portIndicator (γ.port (n + 1)) - portIndicator (γ.port n))) := by
  unfold hoppingCurrent
  rw [map_neg, map_smul]
  congr 2
  funext p
  rw [realBoundary_hopIndicator _ _ (γ.hop n) p]
  rfl

/-- **Continuity.**  The hopping sources satisfy the committed scaled
continuity equation `ρ (n+1) - ρ n + h • ∂ (J n) = 0` identically. -/
theorem hopping_continuity (q h : ℝ) (hh : h ≠ 0) (γ : HoppingPath) (n : ℕ) :
    hoppingLoad q γ (n + 1) - hoppingLoad q γ n +
      h • realBoundary (hoppingCurrent q h γ n) = 0 := by
  rw [realBoundary_hoppingCurrent, smul_neg, smul_smul, mul_div_cancel₀ q hh]
  unfold hoppingLoad
  rw [smul_sub]
  abel

/-! ## Minimal coupling and the augmented field equations -/

/-- The source pairing of the committed window action: the window sum of
`h ⟨J n, A (n+1)⟩ + h ⟨ρ n, φ n⟩`. -/
def sourcePairing (h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) : ℝ :=
  ∑ n ∈ Finset.range (N + 1),
    (h * realSeamInner (J n) (A (n + 1)) + h * realPortInner (ρ n) (φ n))

/-- Declared minimal coupling of the hopping charge: the source pairing at
the hopping sources, with the committed source weights. -/
def interactionA (q h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (γ : HoppingPath) : ℝ :=
  sourcePairing h N A φ (hoppingLoad q γ) (hoppingCurrent q h γ)

/-- The coupled action: the committed window action plus the minimal
coupling. -/
def monopoleCoupledAction (q h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (γ : HoppingPath) : ℝ :=
  windowAction h N A φ ρ J + interactionA q h N A φ γ

theorem windowAction_augment_sourcePairing (h : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (ρ ρ' : ℕ → Fin 12 → ℝ) (J J' : ℕ → Fin 30 → ℝ) :
    windowAction h N A φ (ρ + ρ') (J + J') =
      windowAction h N A φ ρ J + sourcePairing h N A φ ρ' J' :=
  windowAction_augment h N A φ ρ ρ' J J'

/-- **Coupled action equals the window action at the augmented sources.**
No endpoint term appears. -/
theorem monopoleCoupledAction_eq_augmented (q h : ℝ) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (γ : HoppingPath) :
    monopoleCoupledAction q h N A φ ρ J γ =
      windowAction h N A φ (ρ + hoppingLoad q γ) (J + hoppingCurrent q h γ) := by
  rw [windowAction_augment_sourcePairing]
  rfl

/-- Stationarity in the seam potential is the committed scaled Ampere
update at the augmented current. -/
theorem monopole_stationary_A_iff_ampere (q h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (γ : HoppingPath) :
    (∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
      monopoleCoupledAction q h N (A + a) φ ρ J γ =
        monopoleCoupledAction q h N A φ ρ J γ + quadraticRemainder h N a 0) ↔
    (∀ m, m < N → ampereResidual h A φ (J + hoppingCurrent q h γ) m = 0) := by
  simp only [monopoleCoupledAction_eq_augmented]
  exact action_stationary_A_iff_ampere h hh N A φ _ _

/-- Stationarity in the port potential is the committed Gauss constraint
at the augmented load. -/
theorem monopole_stationary_phi_iff_gauss (q h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (γ : HoppingPath) :
    (∀ f : ℕ → Fin 12 → ℝ,
      monopoleCoupledAction q h N A (φ + f) ρ J γ =
        monopoleCoupledAction q h N A φ ρ J γ + quadraticRemainder h N 0 f) ↔
    (∀ n, n < N + 1 →
      realBoundary (electricFieldScaled h A φ n) = ρ n + hoppingLoad q γ n) := by
  simp only [monopoleCoupledAction_eq_augmented]
  exact action_stationary_phi_iff_gauss h hh N A φ _ _

/-- The coupled field equations: the committed scaled Ampere update at the
augmented current and the committed Gauss constraint at the augmented
load. -/
theorem monopole_field_equations (q h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (γ : HoppingPath) :
    ((∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
        monopoleCoupledAction q h N (A + a) φ ρ J γ =
          monopoleCoupledAction q h N A φ ρ J γ + quadraticRemainder h N a 0) ∧
      (∀ f : ℕ → Fin 12 → ℝ,
        monopoleCoupledAction q h N A (φ + f) ρ J γ =
          monopoleCoupledAction q h N A φ ρ J γ + quadraticRemainder h N 0 f)) ↔
    ((∀ m, m < N →
        electricFieldScaled h A φ (m + 1) - electricFieldScaled h A φ m =
          h • (localMaxwellOperator (A (m + 1)) - (J m + hoppingCurrent q h γ m))) ∧
      (∀ n, n < N + 1 →
        realBoundary (electricFieldScaled h A φ n) = ρ n + hoppingLoad q γ n)) := by
  rw [monopole_stationary_A_iff_ampere q h hh, monopole_stationary_phi_iff_gauss q h hh]
  refine and_congr (forall_congr' fun m ↦ imp_congr Iff.rfl ?_) Iff.rfl
  unfold ampereResidual
  rw [sub_eq_zero]
  rfl

/-! ## Gauge change of the source pairing and conservation -/

theorem portInner_sub_right (x y z : Fin 12 → ℝ) :
    realPortInner x (y - z) = realPortInner x y - realPortInner x z := by
  rw [realPortInner_comm, portInner_sub_left, realPortInner_comm y x,
    realPortInner_comm z x]

theorem portInner_smul_right (c : ℝ) (x y : Fin 12 → ℝ) :
    realPortInner x (c • y) = c * realPortInner x y := by
  rw [realPortInner_comm, portInner_smul_left, realPortInner_comm]

/-- The committed scaled continuity residual `ρ (n+1) - ρ n + h • ∂ (J n)`. -/
def continuityResidual (h : ℝ) (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (n : ℕ) : Fin 12 → ℝ :=
  ρ (n + 1) - ρ n + h • realBoundary (J n)

/-- Per-step gauge change of the source density. -/
theorem sourceStep_gauge (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ χ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (n : ℕ) :
    h * realSeamInner (J n) (gaugeTransformA A χ (n + 1)) +
      h * realPortInner (ρ n) (gaugeTransformPhiScaled h φ χ n) =
    (h * realSeamInner (J n) (A (n + 1)) + h * realPortInner (ρ n) (φ n)) +
      realPortInner (χ (n + 1)) (continuityResidual h ρ J n) +
      (realPortInner (χ n) (ρ n) - realPortInner (χ (n + 1)) (ρ (n + 1))) := by
  unfold gaugeTransformA gaugeTransformPhiScaled continuityResidual
  rw [realSeamInner_add_right, portInner_sub_right, portInner_smul_right,
    portInner_sub_right, portInner_add_right, portInner_sub_right,
    portInner_smul_right, realSeamInner_comm (J n) (realCoboundary (χ (n + 1))),
    realCoboundary_boundary_adjoint, realPortInner_comm (ρ n) (χ (n + 1)),
    realPortInner_comm (ρ n) (χ n)]
  have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  linear_combination (realPortInner (χ n) (ρ n) - realPortInner (χ (n + 1)) (ρ n)) * hinv

/-- **Gauge change of the source pairing.**  Under the committed scaled
gauge transformation the source pairing moves by the window sum of the
gauge function at the forward step paired against the continuity
residual, plus the endpoint difference of the gauge function paired
against the load. -/
theorem sourcePairing_gauge (h : ℝ) (hh : h ≠ 0) (N : ℕ) (A : ℕ → Fin 30 → ℝ)
    (φ χ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) :
    sourcePairing h N (gaugeTransformA A χ) (gaugeTransformPhiScaled h φ χ) ρ J =
      sourcePairing h N A φ ρ J +
        (∑ n ∈ Finset.range (N + 1),
          realPortInner (χ (n + 1)) (continuityResidual h ρ J n)) +
        (realPortInner (χ 0) (ρ 0) - realPortInner (χ (N + 1)) (ρ (N + 1))) := by
  unfold sourcePairing
  rw [Finset.sum_congr rfl fun n _ ↦ sourceStep_gauge h hh A φ χ ρ J n,
    Finset.sum_add_distrib, Finset.sum_add_distrib, Finset.sum_range_sub']

/-- **Gauge invariance under endpoint-vanishing gauge functions is the
continuity equation on the interior steps.**  For a general seam current
and port load, the source pairing is invariant under every gauge function
vanishing at both window endpoints exactly when the continuity residual
vanishes at every step `n < N`. -/
theorem sourcePairing_gauge_invariant_iff (h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) :
    (∀ χ : ℕ → Fin 12 → ℝ, χ 0 = 0 → χ (N + 1) = 0 →
      sourcePairing h N (gaugeTransformA A χ) (gaugeTransformPhiScaled h φ χ) ρ J =
        sourcePairing h N A φ ρ J) ↔
    (∀ n, n < N → continuityResidual h ρ J n = 0) := by
  constructor
  · intro hinv n hn
    set r := continuityResidual h ρ J n with hr
    let χ : ℕ → Fin 12 → ℝ := fun m ↦ if m = n + 1 then r else 0
    have h0 : χ 0 = 0 := by simp [χ]
    have hN : χ (N + 1) = 0 := by
      have : N ≠ n := by omega
      simp [χ, this]
    have hg := hinv χ h0 hN
    rw [sourcePairing_gauge h hh, h0, hN, portInner_zero_left, portInner_zero_left,
      sub_zero, add_zero] at hg
    have hsum : (∑ m ∈ Finset.range (N + 1),
        realPortInner (χ (m + 1)) (continuityResidual h ρ J m)) =
        realPortInner r r := by
      have hterm : ∀ m, realPortInner (χ (m + 1)) (continuityResidual h ρ J m) =
          if m = n then realPortInner r r else 0 := by
        intro m
        by_cases hm : m = n
        · subst hm; simp [χ, hr]
        · simp [χ, hm, portInner_zero_left]
      rw [Finset.sum_congr rfl fun m _ ↦ hterm m, Finset.sum_ite_eq']
      simp [Finset.mem_range, Nat.lt_succ_of_lt hn]
    rw [hsum] at hg
    have hzero : realPortInner r r = 0 := by linarith
    exact (realPortInner_self_eq_zero_iff r).mp hzero
  · intro hres χ h0 hN
    rw [sourcePairing_gauge h hh, h0, hN, portInner_zero_left, portInner_zero_left,
      sub_zero, add_zero, Finset.sum_range_succ, hN, portInner_zero_left, add_zero]
    have hz : (∑ n ∈ Finset.range N,
        realPortInner (χ (n + 1)) (continuityResidual h ρ J n)) = 0 := by
      refine Finset.sum_eq_zero fun n hn ↦ ?_
      rw [hres n (Finset.mem_range.mp hn), portInner_zero_right]
    rw [hz, add_zero]

theorem portInner_hoppingLoad (χ : Fin 12 → ℝ) (q : ℝ) (γ : HoppingPath) (n : ℕ) :
    realPortInner χ (hoppingLoad q γ n) = q * χ (γ.port n) := by
  unfold realPortInner
  simp only [hoppingLoad_apply, mul_ite, mul_zero, Finset.sum_ite_eq',
    Finset.mem_univ, if_true]
  ring

/-- **Gauge change of the minimal coupling.**  The interaction moves by
`q` times the gauge function at the occupied port at the first window
node minus at the last window node. -/
theorem interactionA_gauge (q h : ℝ) (hh : h ≠ 0) (N : ℕ) (A : ℕ → Fin 30 → ℝ)
    (φ χ : ℕ → Fin 12 → ℝ) (γ : HoppingPath) :
    interactionA q h N (gaugeTransformA A χ) (gaugeTransformPhiScaled h φ χ) γ =
      interactionA q h N A φ γ +
        q * (χ 0 (γ.port 0) - χ (N + 1) (γ.port (N + 1))) := by
  unfold interactionA
  rw [sourcePairing_gauge h hh, portInner_hoppingLoad, portInner_hoppingLoad]
  have hz : (∑ n ∈ Finset.range (N + 1), realPortInner (χ (n + 1))
      (continuityResidual h (hoppingLoad q γ) (hoppingCurrent q h γ) n)) = 0 := by
    refine Finset.sum_eq_zero fun n _ ↦ ?_
    unfold continuityResidual
    rw [hopping_continuity q h hh γ n, portInner_zero_right]
  rw [hz]
  ring

/-- Gauge invariance of the minimal coupling under every gauge function
vanishing at both window endpoints. -/
theorem interactionA_gauge_invariant (q h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ χ : ℕ → Fin 12 → ℝ) (γ : HoppingPath)
    (h0 : χ 0 = 0) (hN : χ (N + 1) = 0) :
    interactionA q h N (gaugeTransformA A χ) (gaugeTransformPhiScaled h φ χ) γ =
      interactionA q h N A φ γ := by
  rw [interactionA_gauge q h hh, h0, hN]
  simp

/-! ## Work and energy for the hopping charge -/

theorem seamInner_single (x : Fin 30 → ℝ) (e₀ : Fin 30) (c : ℝ) :
    realSeamInner x (fun e ↦ if e = e₀ then c else 0) = c * x e₀ := by
  unfold realSeamInner
  simp only [mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true]
  ring

/-- **Field energy balance at a hop.**  Along the scaled Ampere evolution
driven by the hopping current, the scaled staggered energy moves per step by
`(q / 2)` times the pairing of `E n + E (n+1)` against the oriented hop
indicator.  No energy of the charge is defined here; the statement is the
committed field energy balance specialized to the hopping current, with
the sign following from `energy_balance_scaled` and the current sign
convention. -/
theorem hopping_work_energy (q h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (γ : HoppingPath)
    (hAmp : AmpereEvolutionScaled h A φ (hoppingCurrent q h γ)) (n : ℕ) :
    fieldEnergyScaled h A φ (n + 1) = fieldEnergyScaled h A φ n +
      (q / 2) * realSeamInner
        (electricFieldScaled h A φ n + electricFieldScaled h A φ (n + 1))
        (hopIndicator (γ.port n) (γ.port (n + 1))) := by
  rw [energy_balance_scaled h hh A φ _ hAmp n]
  unfold hoppingCurrent
  rw [seamInner_neg_right, seamInner_smul_right]
  field_simp
  ring

/-- A resting step leaves the scaled staggered energy unchanged. -/
theorem hopping_work_energy_rest (q h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (γ : HoppingPath)
    (hAmp : AmpereEvolutionScaled h A φ (hoppingCurrent q h γ)) (n : ℕ)
    (hn : γ.port (n + 1) = γ.port n) :
    fieldEnergyScaled h A φ (n + 1) = fieldEnergyScaled h A φ n := by
  rw [hopping_work_energy q h hh A φ γ hAmp n, hn, hopIndicator_self,
    realSeamInner_zero_right, mul_zero, add_zero]

/-- A forward hop across seam `e` moves the energy by
`(q / 2) (E n e + E (n+1) e)`. -/
theorem hopping_work_energy_forward (q h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (γ : HoppingPath)
    (hAmp : AmpereEvolutionScaled h A φ (hoppingCurrent q h γ)) (n : ℕ)
    (e : Fin 30) (hl : γ.port n = seamLeft e) (hr : γ.port (n + 1) = seamRight e) :
    fieldEnergyScaled h A φ (n + 1) = fieldEnergyScaled h A φ n +
      (q / 2) * (electricFieldScaled h A φ n e + electricFieldScaled h A φ (n + 1) e) := by
  rw [hopping_work_energy q h hh A φ γ hAmp n, hl, hr, hopIndicator_forward,
    seamInner_single]
  simp only [Pi.add_apply]
  ring

/-- A backward hop across seam `e` moves the energy by
`-(q / 2) (E n e + E (n+1) e)`. -/
theorem hopping_work_energy_backward (q h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (γ : HoppingPath)
    (hAmp : AmpereEvolutionScaled h A φ (hoppingCurrent q h γ)) (n : ℕ)
    (e : Fin 30) (hl : γ.port n = seamRight e) (hr : γ.port (n + 1) = seamLeft e) :
    fieldEnergyScaled h A φ (n + 1) = fieldEnergyScaled h A φ n -
      (q / 2) * (electricFieldScaled h A φ n e + electricFieldScaled h A φ (n + 1) e) := by
  rw [hopping_work_energy q h hh A φ γ hAmp n, hl, hr, hopIndicator_backward,
    seamInner_single]
  simp only [Pi.add_apply]
  ring

/-! ## Bridge to the E-paired candidate -/

/-- Declared crossing path of seam `e`: the smaller endpoint at step `0`,
the larger endpoint afterwards. -/
def crossingPath (e : Fin 30) : HoppingPath where
  port := fun n ↦ if n = 0 then seamLeft e else seamRight e
  hop := by
    intro n
    unfold Hop
    by_cases hn : n = 0
    · subst hn
      exact Or.inr (Or.inl ⟨e, by simp, by simp⟩)
    · exact Or.inl (by simp [hn])

theorem crossingPath_zero (e : Fin 30) : (crossingPath e).port 0 = seamLeft e := rfl

theorem crossingPath_one (e : Fin 30) : (crossingPath e).port 1 = seamRight e := rfl

/-- Forward step difference of the hopping load on the crossing path. -/
theorem crossingPath_load_difference (q : ℝ) (e : Fin 30) (p : Fin 12) :
    (hoppingLoad q (crossingPath e) 1 - hoppingLoad q (crossingPath e) 0) p =
      q * (portIndicator (seamRight e) p - portIndicator (seamLeft e) p) := by
  simp only [Pi.sub_apply, hoppingLoad_apply, crossingPath_zero, crossingPath_one,
    portIndicator]
  split_ifs <;> ring

/-- **Endpoint bridge.**  At both endpoints of seam `e`, the E-paired
induced load of the seam-crossing worldline equals `-(κ / h) (12 - 4φ)`
times the forward step difference of the unit hopping load of the
crossing path. -/
theorem bridge_endpoints (κ h : ℝ) (e : Fin 30) (p : Fin 12)
    (hp : p = seamLeft e ∨ p = seamRight e) :
    inducedLoad κ h (seamCrossingWorldline e) 0 p =
      -(κ / h) * seamCrossingConst *
        (hoppingLoad 1 (crossingPath e) 1 - hoppingLoad 1 (crossingPath e) 0) p := by
  have hlt := (seam_table_sound e).1
  rw [inducedLoad_eq_crossingFactor, seamCrossing_increment,
    crossingPath_load_difference]
  unfold portIndicator
  rcases hp with rfl | rfl
  · rw [crossingFactor_left, if_neg hlt.ne, if_pos rfl]
    ring
  · rw [crossingFactor_right, if_pos rfl, if_neg hlt.ne']
    ring

/-- At the charge-fixed normalization the E-paired load equals the step
difference of the hopping load of charge `q` at both endpoints. -/
theorem bridge_endpoints_charge_fixed (q h : ℝ) (hh : h ≠ 0) (e : Fin 30)
    (p : Fin 12) (hp : p = seamLeft e ∨ p = seamRight e) :
    inducedLoad (-(q * h) / seamCrossingConst) h (seamCrossingWorldline e) 0 p =
      (hoppingLoad q (crossingPath e) 1 - hoppingLoad q (crossingPath e) 0) p := by
  rw [bridge_endpoints _ h e p hp, crossingPath_load_difference,
    crossingPath_load_difference]
  have hc := seamCrossingConst_ne_zero
  field_simp

/-- Endpoint agreement at the larger endpoint holds exactly at the
charge-fixed normalization `κ = -(q h) / (12 - 4φ)`: this is the committed
charge fixing restated. -/
theorem bridge_iff_charge_fixed (κ h q : ℝ) (hh : h ≠ 0) (e : Fin 30) :
    inducedLoad κ h (seamCrossingWorldline e) 0 (seamRight e) =
        (hoppingLoad q (crossingPath e) 1 - hoppingLoad q (crossingPath e) 0)
          (seamRight e) ↔
      κ = -(q * h) / seamCrossingConst := by
  have hlt := (seam_table_sound e).1
  rw [← seam_crossing_charge_fixed κ h q hh e]
  unfold PortChargeIdentification
  rw [crossingPath_load_difference]
  unfold portIndicator
  rw [if_pos rfl, if_neg hlt.ne']
  simp

theorem ray_dot_seamVector_ten_zero :
    dotZ (candidateRayZ 10) (seamVectorZ 0) = ((-2 : ℤ), (0 : ℤ)) := by
  decide

theorem seam_zero_endpoints : seamLeft 0 = 0 ∧ seamRight 0 = 1 := by decide

/-- **Off-endpoint exhibit.**  At seam `0`, port `10` is not an endpoint;
there the E-paired load equals `(κ / h) (12 - 4φ)`, while the hopping
step difference vanishes.  The polarization load is the step difference of
the monopole load at the seam endpoints only. -/
theorem bridge_off_endpoint_exhibit (κ h q : ℝ) :
    inducedLoad κ h (seamCrossingWorldline 0) 0 10 = (κ / h) * seamCrossingConst ∧
    (hoppingLoad q (crossingPath 0) 1 - hoppingLoad q (crossingPath 0) 0) 10 = 0 := by
  constructor
  · rw [inducedLoad_eq_crossingFactor, seamCrossing_increment]
    unfold crossingFactor candidateRay seamVector
    rw [← evalPhi_dotZ, ray_dot_seamVector_ten_zero, seamCrossingConst_eq]
    show -(κ / h) * (boundaryConst * (((-2 : ℤ) : ℝ) + ((0 : ℤ) : ℝ) * Real.goldenRatio)) =
      κ / h * (2 * boundaryConst)
    push_cast
    ring
  · rw [crossingPath_load_difference, seam_zero_endpoints.1, seam_zero_endpoints.2]
    unfold portIndicator
    simp

/-- At nonzero `κ` the off-endpoint E-paired load is nonzero. -/
theorem bridge_off_endpoint_ne_zero (κ h : ℝ) (hk : κ ≠ 0) (hh : h ≠ 0) :
    inducedLoad κ h (seamCrossingWorldline 0) 0 10 ≠ 0 := by
  rw [(bridge_off_endpoint_exhibit κ h 0).1]
  exact mul_ne_zero (div_ne_zero hk hh) seamCrossingConst_ne_zero

/-! ## Non-forcing of the charge -/

/-- **The charge is unconstrained by the field sector.**  Every charge
value gives sources satisfying the committed continuity equation, an
interaction invariant under every endpoint-vanishing gauge function, and
an augmented Gauss constraint satisfiable by a declared load. -/
theorem charge_unconstrained_by_field_sector (h : ℝ) (hh : h ≠ 0) (γ : HoppingPath) :
    ∀ q : ℝ,
      (∀ n, continuityResidual h (hoppingLoad q γ) (hoppingCurrent q h γ) n = 0) ∧
      (∀ (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ χ : ℕ → Fin 12 → ℝ),
        χ 0 = 0 → χ (N + 1) = 0 →
        interactionA q h N (gaugeTransformA A χ) (gaugeTransformPhiScaled h φ χ) γ =
          interactionA q h N A φ γ) ∧
      (∀ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ), ∃ ρ : ℕ → Fin 12 → ℝ, ∀ n,
        realBoundary (electricFieldScaled h A φ n) = ρ n + hoppingLoad q γ n) := by
  intro q
  refine ⟨fun n ↦ hopping_continuity q h hh γ n,
    fun N A φ χ h0 hN ↦ interactionA_gauge_invariant q h hh N A φ χ γ h0 hN,
    fun A φ ↦ ⟨fun n ↦ realBoundary (electricFieldScaled h A φ n) - hoppingLoad q γ n,
      fun n ↦ (sub_add_cancel _ _).symm⟩⟩

/-- Two distinct charges give two distinct hopping loads and two distinct
total loads on the same path. -/
theorem two_charges_two_sources (q₁ q₂ : ℝ) (hq : q₁ ≠ q₂) (γ : HoppingPath) :
    hoppingLoad q₁ γ ≠ hoppingLoad q₂ γ ∧
    (∑ p : Fin 12, hoppingLoad q₁ γ 0 p) ≠ (∑ p : Fin 12, hoppingLoad q₂ γ 0 p) := by
  refine ⟨fun heq ↦ hq ?_, ?_⟩
  · have := congrFun (congrFun heq 0) (γ.port 0)
    simpa [hoppingLoad_apply] using this
  · rw [hoppingLoad_total, hoppingLoad_total]
    exact hq

end

end OPH.PortChargeMinimalCoupling

/- Axiom audit: standard axioms only (`propext`, `Classical.choice`,
`Quot.sound`); no `sorry`, no `native_decide`, no project axiom. -/

#print axioms OPH.PortChargeMinimalCoupling.realBoundary_hopIndicator
#print axioms OPH.PortChargeMinimalCoupling.hopping_continuity
#print axioms OPH.PortChargeMinimalCoupling.hoppingLoad_total
#print axioms OPH.PortChargeMinimalCoupling.hoppingCurrent_rest
#print axioms OPH.PortChargeMinimalCoupling.monopoleCoupledAction_eq_augmented
#print axioms OPH.PortChargeMinimalCoupling.monopole_stationary_A_iff_ampere
#print axioms OPH.PortChargeMinimalCoupling.monopole_stationary_phi_iff_gauss
#print axioms OPH.PortChargeMinimalCoupling.monopole_field_equations
#print axioms OPH.PortChargeMinimalCoupling.sourcePairing_gauge
#print axioms OPH.PortChargeMinimalCoupling.sourcePairing_gauge_invariant_iff
#print axioms OPH.PortChargeMinimalCoupling.interactionA_gauge
#print axioms OPH.PortChargeMinimalCoupling.interactionA_gauge_invariant
#print axioms OPH.PortChargeMinimalCoupling.hopping_work_energy
#print axioms OPH.PortChargeMinimalCoupling.hopping_work_energy_rest
#print axioms OPH.PortChargeMinimalCoupling.hopping_work_energy_forward
#print axioms OPH.PortChargeMinimalCoupling.hopping_work_energy_backward
#print axioms OPH.PortChargeMinimalCoupling.bridge_endpoints
#print axioms OPH.PortChargeMinimalCoupling.bridge_endpoints_charge_fixed
#print axioms OPH.PortChargeMinimalCoupling.bridge_iff_charge_fixed
#print axioms OPH.PortChargeMinimalCoupling.bridge_off_endpoint_exhibit
#print axioms OPH.PortChargeMinimalCoupling.bridge_off_endpoint_ne_zero
#print axioms OPH.PortChargeMinimalCoupling.charge_unconstrained_by_field_sector
#print axioms OPH.PortChargeMinimalCoupling.two_charges_two_sources
