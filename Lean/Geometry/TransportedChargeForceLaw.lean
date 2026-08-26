import Geometry.WorldlineHopTransport

set_option autoImplicit false

open scoped BigOperators

namespace OPH.TransportedChargeForceLaw

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
open OPH.WorldlineHopTransport
open OPH.C1Lorentz (Spatial Herm2 spatialDot spatialNormSq lorentzB lorentzQ)

/-!
# A discrete force law for the transported charge

STATUS.  Candidate module on the coupled-action row.  The transported
action of `Geometry/WorldlineHopTransport.lean` is the monopole coupled
action of the hopping path plus the committed clock action of the
generated path.  The path ranges over a discrete set of step sequences,
so the worldline sector has no real first variation; this module declares
the class of closed two-step variations of a seam-step worldline (the
step pair `(s_j, s_{j+1})` is replaced by an admissible pair with the same
endpoint, every other step is kept), computes the exact action difference
under every such variation with the potentials held fixed, and reads off
the stationarity condition as an exact balance on the scaled electric
seam field.  Every statement is an identity or an equivalence on the
finite complex; the field-sector parts of the action are unchanged by
construction.

WHAT IS PROVED.

1. Closed two-step variations.  `replaceTwo` replaces the steps `j` and
   `j + 1` of a step sequence.  A closed two-step variation of a seam-step
   worldline at `j` (`ClosedTwoStepVariation`) is an admissible pair of
   steps whose two-step target is the port `w.port (j + 2)`.  The varied
   worldline (`ClosedTwoStepVariation.worldline`) is again a seam-step
   worldline; its ports agree with the original ones at every index other
   than `j + 1` (`varied_port_eq`), and at `j + 1` its port is the declared
   intermediate port (`varied_port_mid`).
2. Exact action difference.  With `h ≠ 0`, the minimal-coupling
   interaction of a hopping path is the window sum of the per-step
   potential pairing `q (h φ n (u) - ⟨hopIndicator u v, A (n+1)⟩)` at the
   occupied port `u` and the hop `u → v` (`interactionA_eq_sum_hopPairing`).
   Under a closed two-step variation inside the window the transported
   action moves by the clock difference, `(|s_j|² + |s_{j+1}|²) -
   (|s'|² + |s''|²)` with `|rest|² = 0` and `|crossing|² = 4`, plus `q`
   times the potential pairing along the alternative two-step route minus
   along the original route (`transportedAction_closed_variation`); the
   committed window action cancels identically and the declared unit `τ`
   drops out of the difference.
3. Stationarity and local minimality.  `SingleStepStationary` (every
   closed two-step variation inside the window leaves the action fixed)
   and `SingleStepLocalMin` (none lowers it) are declared.  For the
   exchange of a forward crossing of `e` followed by a rest against a rest
   followed by the same crossing, the clock difference vanishes and the
   action difference is `q h E_{j+1}(e)`, the scaled electric seam field
   on `e` at the step of the delayed crossing (`exchange_action_difference`);
   with `q ≠ 0` and `h ≠ 0` the exchange is stationary exactly when
   `E_{j+1}(e) = 0` (`exchange_stationary_iff`).  The backward exchange
   gives `-q h E_{j+1}(e)` (`exchange_backward_action_difference`).  The
   balance is pure field: the clock cost of the two orders is equal, so
   the stationarity condition is the vanishing of the electromotive
   pairing `(A (j+2) e - A (j+1) e) + h (φ (j+1) (right) - φ (j+1) (left))`
   (`exchange_balance`).
4. Round trips.  Replacing two rests at the smaller endpoint of `e` by the
   round trip forward-then-backward across `e` moves the action by
   `-8 - q h E_{j+1}(e)` (`roundTrip_action_difference`): the clock term
   `2 (τ² - 4) - 2 τ²` and `q` times the time difference of `A` on `e`
   plus the port-potential difference, the discrete work of the electric
   field around the two-step loop.  The round trip lowers the action
   exactly when `-8 < q h E_{j+1}(e)` (`roundTrip_lowers_iff`); at zero
   field it lowers the action by `8`, so rest is not a local minimum of
   the transported action in the committed `(+---)` signature
   (`rest_not_localMin_zero_field`).
5. Static fields.  With `A` and `φ` constant in the step index the
   interaction difference of every closed two-step variation is `q` times
   the port-potential difference at the two intermediate ports minus `q`
   times the circulation of `A` around the port loop
   (`static_interaction_difference`, with `hopIndicator_symm` for the
   return legs); the two-step round trip has zero
   circulation (`roundTrip_circulation_zero`), so its interaction
   difference is `q h (φ (right) - φ (left))` and vanishes on an
   equipotential seam, while the total difference is `-8`
   (`static_roundTrip_difference`, `static_rest_interaction_stationary`).
6. Non-forcing.  Two distinct charges give two distinct exchange balances
   whenever the seam field is nonzero (`two_charges_two_balances`); every
   closed two-step difference is the same at any two declared units
   (`difference_forgets_unit`); the potentials, the charge, the unit, and
   the class of variations are declared, and no continuum limit is taken.

PRIOR WORK.  `ChargeFixedInteraction.coupled_worldline_equation` is the
worldline equation of the E-paired route, obtained from real variations of
a `Herm2` path; it keeps its own sector and is not restated here.
`PortChargeMinimalCoupling.hopping_work_energy` is the field energy
balance at a hop along the Ampere evolution; the electric pairing that
appears here is the same seam field `electricFieldScaled`.
`CommonWorldJointAction.stationarity_not_minimality` exhibits the clock
action lowered by a spatial excursion; `rest_not_localMin_zero_field` is
the port-hopping form of that observation on the transported action.
`WorldlineHopTransport.clockAction_generated` supplies the clock action of
the generated path; `interactionA_gauge` supplies the gauge behaviour of
the interaction.  The single-site variation vocabulary of
`Lean/Variational/` (real single-site replacements of a path) is the
continuum-valued analogue; the variations here range over the finite step
alphabet, so no real first variation exists and the difference is exact.

ROWS TOUCHED.  The coupled-action row (the transported action, the
variation class, and the balance are declared here); the source clock and
duration row (the unit `τ` and the step `h` are declared; `τ` drops out of
every difference); the physical spacetime attachment row (no attachment
of a port to a spacetime point is supplied); the light-signal row (no
signal propagation is attached); the laboratory clock and energy
calibration import (no unit or readout is attached to `q`, `h`, `τ`, or
the seam field); the gravitation-route energy identification (no
identification of the clock term with an energy is made).  The module
discharges none of these rows.

NEGATIVES CITED.  The Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`): realized histories
select no velocity curvature or Legendre map, so the clock action shape,
the coupling shape, and the variation class are declared enrichments;
cited at scope only.

CONVENTIONS.  Signature `(+---)`; `lorentzQ v = v.1 ^ 2 - |v.2| ^ 2`; a
crossing step of the generated path has Lorentz square `τ² - 4`, a rest
`τ²`.  Seam orientation from `seamLeft e` to `seamRight e`; the hopping
current of a forward hop is `-(q / h)` on `e`; the interaction weights
are `h ⟨J n, A (n+1)⟩ + h ⟨ρ n, φ n⟩`.  Scaled electric field
`E n = -(A (n+1) - A n) / h - d (φ n)` with `d φ e = φ (right) - φ (left)`.
Differences are written varied minus original.

FALSIFIER.  The module is wrong if some closed two-step variation moves
the committed window action, if the exchange difference carries a clock
term, or if the round-trip clock difference differs from `-8`.

Axiom audit.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`; no `sorry`, no
`native_decide`, no project axiom.
-/

/-! ## Closed two-step variations of a seam-step worldline -/

/-- Replace the steps `j` and `j + 1` of a step sequence. -/
def replaceTwo (steps : ℕ → SeamStep) (j : ℕ) (s' s'' : SeamStep) : ℕ → SeamStep :=
  fun n ↦ if n = j then s' else if n = j + 1 then s'' else steps n

theorem replaceTwo_at (steps : ℕ → SeamStep) (j : ℕ) (s' s'' : SeamStep) :
    replaceTwo steps j s' s'' j = s' := by
  simp [replaceTwo]

theorem replaceTwo_succ (steps : ℕ → SeamStep) (j : ℕ) (s' s'' : SeamStep) :
    replaceTwo steps j s' s'' (j + 1) = s'' := by
  simp [replaceTwo]

theorem replaceTwo_other (steps : ℕ → SeamStep) (j : ℕ) (s' s'' : SeamStep) (n : ℕ)
    (h1 : n ≠ j) (h2 : n ≠ j + 1) : replaceTwo steps j s' s'' n = steps n := by
  simp [replaceTwo, h1, h2]

/-- Declared closed two-step variation of a seam-step worldline at step
`j`: an admissible alternative step at the port `w.port j`, an admissible
alternative step at the port reached, with the same two-step endpoint
`w.port (j + 2)`. -/
structure ClosedTwoStepVariation (w : SeamStepWorldline) (j : ℕ) where
  /-- The alternative first step. -/
  s' : SeamStep
  /-- The alternative second step. -/
  s'' : SeamStep
  /-- The first alternative step is admissible at the port reached at `j`. -/
  adm1 : StepAdmissible (w.port j) s'
  /-- The second alternative step is admissible at the intermediate port. -/
  adm2 : StepAdmissible (stepTarget (w.port j) s') s''
  /-- The alternative route ends at the original port `w.port (j + 2)`. -/
  closed : stepTarget (stepTarget (w.port j) s') s'' = w.port (j + 2)

namespace ClosedTwoStepVariation

variable {w : SeamStepWorldline} {j : ℕ} (v : ClosedTwoStepVariation w j)

/-- The intermediate port of the alternative route. -/
def mid : Fin 12 := stepTarget (w.port j) v.s'

/-- The varied step sequence. -/
def steps : ℕ → SeamStep := replaceTwo w.steps j v.s' v.s''

/-- The varied port sequence, before the admissibility proof. -/
def portRaw : ℕ → Fin 12 := portSeq w.start v.steps

theorem portRaw_succ (n : ℕ) : v.portRaw (n + 1) = stepTarget (v.portRaw n) (v.steps n) := rfl

/-- Below the varied pair the ports agree with the original ones. -/
theorem portRaw_le (n : ℕ) (hn : n ≤ j) : v.portRaw n = w.port n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    rw [portRaw_succ, ih (by omega), w.port_succ]
    unfold steps
    rw [replaceTwo_other _ _ _ _ _ (by omega) (by omega)]

/-- At `j + 1` the varied port is the intermediate port. -/
theorem portRaw_mid : v.portRaw (j + 1) = v.mid := by
  rw [portRaw_succ, portRaw_le v j le_rfl]
  unfold steps mid
  rw [replaceTwo_at]

/-- At `j + 2` the varied port is the original port: the closure clause. -/
theorem portRaw_end : v.portRaw (j + 2) = w.port (j + 2) := by
  rw [portRaw_succ, portRaw_mid]
  unfold steps mid
  rw [replaceTwo_succ]
  exact v.closed

/-- Beyond the varied pair the ports agree with the original ones. -/
theorem portRaw_ge (n : ℕ) (hn : j + 2 ≤ n) : v.portRaw n = w.port n := by
  obtain ⟨k, rfl⟩ : ∃ k, n = j + 2 + k := ⟨n - (j + 2), by omega⟩
  induction k with
  | zero => exact v.portRaw_end
  | succ k ih =>
    rw [show j + 2 + (k + 1) = (j + 2 + k) + 1 by omega, portRaw_succ, ih (by omega),
      w.port_succ]
    unfold steps
    rw [replaceTwo_other _ _ _ _ _ (by omega) (by omega)]

/-- The varied ports agree with the original ones at every index other
than `j + 1`. -/
theorem portRaw_eq (n : ℕ) (hn : n ≠ j + 1) : v.portRaw n = w.port n := by
  rcases Nat.lt_or_ge n (j + 1) with h | h
  · exact v.portRaw_le n (by omega)
  · exact v.portRaw_ge n (by omega)

/-- The varied worldline: the varied step sequence is admissible at every
step. -/
def worldline : SeamStepWorldline where
  start := w.start
  steps := v.steps
  adm := by
    intro n
    show StepAdmissible (v.portRaw n) (v.steps n)
    by_cases h1 : n = j
    · subst h1
      unfold steps
      rw [replaceTwo_at, v.portRaw_le n le_rfl]
      exact v.adm1
    · by_cases h2 : n = j + 1
      · subst h2
        unfold steps
        rw [replaceTwo_succ, v.portRaw_mid]
        exact v.adm2
      · unfold steps
        rw [replaceTwo_other _ _ _ _ _ h1 h2, v.portRaw_eq n h2]
        exact w.adm n

theorem worldline_steps (n : ℕ) : v.worldline.steps n = replaceTwo w.steps j v.s' v.s'' n := rfl

theorem worldline_port (n : ℕ) : v.worldline.port n = v.portRaw n := rfl

/-- **Varied ports.**  The varied worldline occupies the original port at
every index other than `j + 1`. -/
theorem varied_port_eq (n : ℕ) (hn : n ≠ j + 1) : v.worldline.port n = w.port n :=
  v.portRaw_eq n hn

/-- At `j + 1` the varied worldline occupies the intermediate port. -/
theorem varied_port_mid : v.worldline.port (j + 1) = v.mid := v.portRaw_mid

end ClosedTwoStepVariation

/-! ## The per-step potential pairing of a hop -/

noncomputable section

/-- The potential pairing of one step: `h` times the port potential at the
occupied port minus the seam potential of the next node paired against the
oriented hop indicator.  The minimal-coupling interaction is `q` times its
window sum. -/
def potentialPairing (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ)
    (u v : Fin 12) : ℝ :=
  h * φ n u - realSeamInner (hopIndicator u v) (A (n + 1))

/-- The interaction contribution of one step of a hopping charge. -/
def hopPairing (q h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ)
    (u v : Fin 12) : ℝ :=
  q * potentialPairing h A φ n u v

theorem potentialPairing_rest (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ)
    (u : Fin 12) : potentialPairing h A φ n u u = h * φ n u := by
  unfold potentialPairing
  rw [hopIndicator_self, realSeamInner_zero_left, sub_zero]

theorem potentialPairing_forward (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ)
    (e : Fin 30) :
    potentialPairing h A φ n (seamLeft e) (seamRight e) =
      h * φ n (seamLeft e) - A (n + 1) e := by
  unfold potentialPairing
  rw [hopIndicator_forward, realSeamInner_comm, seamInner_single, one_mul]

theorem potentialPairing_backward (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ)
    (e : Fin 30) :
    potentialPairing h A φ n (seamRight e) (seamLeft e) =
      h * φ n (seamRight e) + A (n + 1) e := by
  unfold potentialPairing
  rw [hopIndicator_backward, realSeamInner_comm, seamInner_single]
  ring

/-- **The interaction as a window sum of hop pairings.**  With `h ≠ 0` the
minimal-coupling interaction of a hopping path is the window sum of the
hop pairings at the occupied port and the hop of each step. -/
theorem interactionA_eq_sum_hopPairing (q h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (γ : HoppingPath) :
    interactionA q h N A φ γ =
      ∑ n ∈ Finset.range (N + 1), hopPairing q h A φ n (γ.port n) (γ.port (n + 1)) := by
  unfold interactionA sourcePairing
  refine Finset.sum_congr rfl fun n _ ↦ ?_
  unfold hoppingCurrent hopPairing potentialPairing
  rw [realSeamInner_comm, seamInner_neg_right, seamInner_smul_right, realPortInner_comm,
    portInner_hoppingLoad, realSeamInner_comm]
  field_simp
  ring

/-- Two functions agreeing off the sites `j` and `j + 1` have window sums
differing by the two-site difference, when both sites lie in the window. -/
theorem sum_two_site_difference (f f' : ℕ → ℝ) (M j : ℕ) (hj : j + 1 < M)
    (hne : ∀ n, n ≠ j → n ≠ j + 1 → f' n = f n) :
    (∑ n ∈ Finset.range M, f' n) - (∑ n ∈ Finset.range M, f n) =
      (f' j + f' (j + 1)) - (f j + f (j + 1)) := by
  have hpt : ∀ n, f' n - f n =
      (if n = j then f' j - f j else 0) + (if n = j + 1 then f' (j + 1) - f (j + 1) else 0) := by
    intro n
    by_cases h1 : n = j
    · subst h1
      simp
    · by_cases h2 : n = j + 1
      · subst h2
        simp
      · rw [hne n h1 h2, if_neg h1, if_neg h2]
        ring
  rw [← Finset.sum_sub_distrib, Finset.sum_congr rfl fun n _ ↦ hpt n, Finset.sum_add_distrib,
    Finset.sum_ite_eq', Finset.sum_ite_eq', if_pos (Finset.mem_range.mpr (by omega)),
    if_pos (Finset.mem_range.mpr hj)]
  ring

/-- The potential pairing along a two-step route `u → p → x` starting at
step `j`. -/
def routePairing (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (j : ℕ)
    (u p x : Fin 12) : ℝ :=
  potentialPairing h A φ j u p + potentialPairing h A φ (j + 1) p x

/-- The clock difference of a closed two-step variation: the original step
norms minus the varied step norms; the unit `τ` cancels. -/
def clockDifference (w : SeamStepWorldline) {j : ℕ} (v : ClosedTwoStepVariation w j) : ℝ :=
  (stepNormSq (w.steps j) + stepNormSq (w.steps (j + 1))) - (stepNormSq v.s' + stepNormSq v.s'')

/-- The interaction difference of a closed two-step variation: `q` times the
route pairing of the alternative route minus that of the original route. -/
def interactionDifference (q h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (w : SeamStepWorldline) {j : ℕ} (v : ClosedTwoStepVariation w j) : ℝ :=
  q * (routePairing h A φ j (w.port j) v.mid (w.port (j + 2)) -
    routePairing h A φ j (w.port j) (w.port (j + 1)) (w.port (j + 2)))

/-- The clock action of the varied generated path minus that of the
original, over a window containing both varied steps. -/
theorem clockAction_closed_variation (τ : ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ClosedTwoStepVariation w j) (M : ℕ) (hj : j + 1 < M) :
    clockAction M (generatedPath τ v.worldline) - clockAction M (generatedPath τ w) =
      clockDifference w v := by
  rw [clockAction_generated, clockAction_generated,
    sum_two_site_difference _ _ M j hj (fun n h1 h2 ↦ by
      rw [v.worldline_steps, replaceTwo_other _ _ _ _ _ h1 h2])]
  rw [v.worldline_steps, v.worldline_steps, replaceTwo_at, replaceTwo_succ]
  unfold clockDifference
  ring

/-- The interaction of the varied hopping path minus that of the original. -/
theorem interactionA_closed_variation (q h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ClosedTwoStepVariation w j) (hj : j + 1 < N + 1) :
    interactionA q h N A φ (hoppingPath v.worldline) -
        interactionA q h N A φ (hoppingPath w) =
      interactionDifference q h A φ w v := by
  rw [interactionA_eq_sum_hopPairing q h hh, interactionA_eq_sum_hopPairing q h hh,
    sum_two_site_difference _ _ (N + 1) j hj (fun n h1 h2 ↦ by
      rw [hoppingPath_port, hoppingPath_port, hoppingPath_port, hoppingPath_port,
        v.varied_port_eq n h2, v.varied_port_eq (n + 1) (by omega)])]
  simp only [hoppingPath_port]
  rw [v.varied_port_eq j (by omega), v.varied_port_mid, v.varied_port_eq (j + 2) (by omega)]
  unfold interactionDifference routePairing hopPairing
  ring

/-- **Exact action difference under a closed two-step variation.**  With the
potentials held fixed and both varied steps inside the window, the
transported action moves by the clock difference plus the interaction
difference; the committed window action cancels identically, and the
declared unit drops out. -/
theorem transportedAction_closed_variation (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ClosedTwoStepVariation w j) (hj : j + 1 < N + 1) :
    transportedAction q h τ N A φ ρ J v.worldline - transportedAction q h τ N A φ ρ J w =
      clockDifference w v + interactionDifference q h A φ w v := by
  rw [← clockAction_closed_variation τ w v (N + 1) hj,
    ← interactionA_closed_variation q h hh N A φ w v hj]
  unfold transportedAction monopoleCoupledAction
  ring

end

/-! ## Stationarity and local minimality under closed two-step variations -/

noncomputable section

/-- Declared single-step stationarity: every closed two-step variation with
both varied steps inside the window leaves the transported action fixed.
The class of variations is declared; the potentials are held fixed. -/
def SingleStepStationary (q h τ : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) : Prop :=
  ∀ j, j + 1 < N + 1 → ∀ v : ClosedTwoStepVariation w j,
    transportedAction q h τ N A φ ρ J v.worldline = transportedAction q h τ N A φ ρ J w

/-- Declared single-step local minimality: no closed two-step variation
inside the window lowers the transported action. -/
def SingleStepLocalMin (q h τ : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) : Prop :=
  ∀ j, j + 1 < N + 1 → ∀ v : ClosedTwoStepVariation w j,
    transportedAction q h τ N A φ ρ J w ≤ transportedAction q h τ N A φ ρ J v.worldline

/-- The scaled electric seam field on `e`, multiplied by `h`: minus the
forward time difference of the seam potential minus `h` times the port
coboundary of the port potential. -/
theorem h_mul_electricFieldScaled (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (n : ℕ) (e : Fin 30) :
    h * electricFieldScaled h A φ n e =
      -(A (n + 1) e - A n e) - h * (φ n (seamRight e) - φ n (seamLeft e)) := by
  unfold electricFieldScaled
  simp only [Pi.sub_apply, Pi.neg_apply, Pi.smul_apply, smul_eq_mul, realCoboundary_apply]
  field_simp

/-! ### The exchange of a crossing and a rest -/

/-- The port at a forward crossing step is the smaller endpoint. -/
theorem port_of_forward (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hf : w.steps j = .forward e) : w.port j = seamLeft e := by
  have h := w.adm j
  rw [hf] at h
  exact h

/-- The port at a backward crossing step is the larger endpoint. -/
theorem port_of_backward (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hb : w.steps j = .backward e) : w.port j = seamRight e := by
  have h := w.adm j
  rw [hb] at h
  exact h

theorem port_succ_of_rest (w : SeamStepWorldline) (j : ℕ) (hr : w.steps j = .rest) :
    w.port (j + 1) = w.port j := by
  rw [w.port_succ, hr]
  rfl

theorem port_succ_of_forward (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hf : w.steps j = .forward e) : w.port (j + 1) = seamRight e := by
  rw [w.port_succ, hf]
  rfl

theorem port_succ_of_backward (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hb : w.steps j = .backward e) : w.port (j + 1) = seamLeft e := by
  rw [w.port_succ, hb]
  rfl

/-- The exchange variation: a forward crossing of `e` at `j` followed by a
rest is replaced by a rest followed by the same crossing. -/
def exchangeVariation (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hf : w.steps j = .forward e) (hr : w.steps (j + 1) = .rest) :
    ClosedTwoStepVariation w j where
  s' := .rest
  s'' := .forward e
  adm1 := trivial
  adm2 := port_of_forward w j e hf
  closed := by
    show seamRight e = w.port (j + 2)
    rw [port_succ_of_rest w (j + 1) hr, port_succ_of_forward w j e hf]

/-- **Exchange difference.**  Delaying a forward crossing of `e` by one
rest moves the transported action by `q h E_{j+1}(e)`, the scaled electric
seam field on `e` at the step of the delayed crossing; the clock difference
is zero. -/
theorem exchange_action_difference (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hf : w.steps j = .forward e) (hr : w.steps (j + 1) = .rest) :
    transportedAction q h τ N A φ ρ J (exchangeVariation w j e hf hr).worldline -
        transportedAction q h τ N A φ ρ J w =
      q * (h * electricFieldScaled h A φ (j + 1) e) ∧
    clockDifference w (exchangeVariation w j e hf hr) = 0 := by
  have hpj := port_of_forward w j e hf
  have hpj1 := port_succ_of_forward w j e hf
  have hpj2 := port_succ_of_rest w (j + 1) hr
  have hmid : (exchangeVariation w j e hf hr).mid = seamLeft e := hpj
  have hclock : clockDifference w (exchangeVariation w j e hf hr) = 0 := by
    unfold clockDifference
    rw [hf, hr]
    show (stepNormSq (.forward e) + stepNormSq .rest) -
      (stepNormSq .rest + stepNormSq (.forward e)) = 0
    ring
  refine ⟨?_, hclock⟩
  rw [transportedAction_closed_variation q h τ hh N A φ ρ J w _ hj, hclock,
    h_mul_electricFieldScaled h hh]
  unfold interactionDifference routePairing
  rw [hmid, hpj2, hpj1, hpj, potentialPairing_rest, potentialPairing_forward,
    potentialPairing_forward, potentialPairing_rest]
  ring

/-- **The exchange balance.**  The interaction difference of the exchange is
`q` times the electromotive pairing `-(A (j+2) e - A (j+1) e) - h (φ (j+1)
(right) - φ (j+1) (left))`, and the clock cost of the two orders is equal:
the balance is pure field. -/
theorem exchange_balance (q h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hf : w.steps j = .forward e) (hr : w.steps (j + 1) = .rest) :
    interactionDifference q h A φ w (exchangeVariation w j e hf hr) =
      q * (-(A (j + 2) e - A (j + 1) e) -
        h * (φ (j + 1) (seamRight e) - φ (j + 1) (seamLeft e))) ∧
    clockDifference w (exchangeVariation w j e hf hr) = 0 := by
  have hpj := port_of_forward w j e hf
  have hpj1 := port_succ_of_forward w j e hf
  have hpj2 := port_succ_of_rest w (j + 1) hr
  have hmid : (exchangeVariation w j e hf hr).mid = seamLeft e := hpj
  constructor
  · unfold interactionDifference routePairing
    rw [hmid, hpj2, hpj1, hpj, potentialPairing_rest, potentialPairing_forward,
      potentialPairing_forward, potentialPairing_rest]
    ring
  · unfold clockDifference
    rw [hf, hr]
    show (stepNormSq (.forward e) + stepNormSq .rest) -
      (stepNormSq .rest + stepNormSq (.forward e)) = 0
    ring

/-- **Exchange stationarity.**  With `q ≠ 0` and `h ≠ 0`, the exchange leaves
the action fixed exactly when the scaled electric seam field on `e`
vanishes at the step of the delayed crossing. -/
theorem exchange_stationary_iff (q h τ : ℝ) (hq : q ≠ 0) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hf : w.steps j = .forward e) (hr : w.steps (j + 1) = .rest) :
    transportedAction q h τ N A φ ρ J (exchangeVariation w j e hf hr).worldline =
        transportedAction q h τ N A φ ρ J w ↔
      electricFieldScaled h A φ (j + 1) e = 0 := by
  rw [← sub_eq_zero, (exchange_action_difference q h τ hh N A φ ρ J w j hj e hf hr).1]
  constructor
  · intro h0
    rcases mul_eq_zero.mp h0 with h0 | h0
    · exact absurd h0 hq
    · rcases mul_eq_zero.mp h0 with h0 | h0
      · exact absurd h0 hh
      · exact h0
  · intro h0
    rw [h0, mul_zero, mul_zero]

/-- The backward exchange variation: a backward crossing of `e` at `j`
followed by a rest is replaced by a rest followed by the same crossing. -/
def exchangeBackwardVariation (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hb : w.steps j = .backward e) (hr : w.steps (j + 1) = .rest) :
    ClosedTwoStepVariation w j where
  s' := .rest
  s'' := .backward e
  adm1 := trivial
  adm2 := port_of_backward w j e hb
  closed := by
    show seamLeft e = w.port (j + 2)
    rw [port_succ_of_rest w (j + 1) hr, port_succ_of_backward w j e hb]

/-- **Backward exchange difference.**  Delaying a backward crossing of `e`
by one rest moves the transported action by `-q h E_{j+1}(e)`. -/
theorem exchange_backward_action_difference (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hb : w.steps j = .backward e) (hr : w.steps (j + 1) = .rest) :
    transportedAction q h τ N A φ ρ J (exchangeBackwardVariation w j e hb hr).worldline -
        transportedAction q h τ N A φ ρ J w =
      -(q * (h * electricFieldScaled h A φ (j + 1) e)) := by
  have hpj := port_of_backward w j e hb
  have hpj1 := port_succ_of_backward w j e hb
  have hpj2 := port_succ_of_rest w (j + 1) hr
  have hmid : (exchangeBackwardVariation w j e hb hr).mid = seamRight e := hpj
  have hclock : clockDifference w (exchangeBackwardVariation w j e hb hr) = 0 := by
    unfold clockDifference
    rw [hb, hr]
    show (stepNormSq (.backward e) + stepNormSq .rest) -
      (stepNormSq .rest + stepNormSq (.backward e)) = 0
    ring
  rw [transportedAction_closed_variation q h τ hh N A φ ρ J w _ hj, hclock,
    h_mul_electricFieldScaled h hh]
  unfold interactionDifference routePairing
  rw [hmid, hpj2, hpj1, hpj, potentialPairing_rest, potentialPairing_backward,
    potentialPairing_backward, potentialPairing_rest]
  ring

end

/-! ### Round trips from rest -/

noncomputable section

/-- The round-trip variation: two rests at the smaller endpoint of `e` are
replaced by a forward crossing of `e` followed by a backward crossing. -/
def roundTripVariation (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hu : w.port j = seamLeft e) (h1 : w.steps j = .rest) (h2 : w.steps (j + 1) = .rest) :
    ClosedTwoStepVariation w j where
  s' := .forward e
  s'' := .backward e
  adm1 := hu
  adm2 := rfl
  closed := by
    show seamLeft e = w.port (j + 2)
    rw [port_succ_of_rest w (j + 1) h2, port_succ_of_rest w j h1, hu]

/-- **Round-trip difference.**  Replacing two rests by the round trip across
`e` moves the transported action by `-8 - q h E_{j+1}(e)`: the clock term
`2 (τ² - 4) - 2 τ² = -8` plus `q` times the electromotive pairing around
the two-step loop, the time difference of `A` on `e` plus `h` times the
port-potential difference at the step between the two crossings. -/
theorem roundTrip_action_difference (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hu : w.port j = seamLeft e) (h1 : w.steps j = .rest) (h2 : w.steps (j + 1) = .rest) :
    transportedAction q h τ N A φ ρ J (roundTripVariation w j e hu h1 h2).worldline -
        transportedAction q h τ N A φ ρ J w =
      -8 - q * (h * electricFieldScaled h A φ (j + 1) e) ∧
    clockDifference w (roundTripVariation w j e hu h1 h2) = -8 ∧
    interactionDifference q h A φ w (roundTripVariation w j e hu h1 h2) =
      q * ((A (j + 2) e - A (j + 1) e) +
        h * (φ (j + 1) (seamRight e) - φ (j + 1) (seamLeft e))) := by
  have hpj1 : w.port (j + 1) = seamLeft e := by rw [port_succ_of_rest w j h1, hu]
  have hpj2 : w.port (j + 2) = seamLeft e := by rw [port_succ_of_rest w (j + 1) h2, hpj1]
  have hmid : (roundTripVariation w j e hu h1 h2).mid = seamRight e := by
    show stepTarget (w.port j) (.forward e) = seamRight e
    rfl
  have hclock : clockDifference w (roundTripVariation w j e hu h1 h2) = -8 := by
    unfold clockDifference
    rw [h1, h2]
    show (stepNormSq .rest + stepNormSq .rest) -
      (stepNormSq (.forward e) + stepNormSq (.backward e)) = -8
    unfold stepNormSq seamNormSq
    ring
  have hint : interactionDifference q h A φ w (roundTripVariation w j e hu h1 h2) =
      q * ((A (j + 2) e - A (j + 1) e) +
        h * (φ (j + 1) (seamRight e) - φ (j + 1) (seamLeft e))) := by
    unfold interactionDifference routePairing
    rw [hmid, hpj2, hpj1, hu, potentialPairing_rest, potentialPairing_forward,
      potentialPairing_backward, potentialPairing_rest]
    ring
  refine ⟨?_, hclock, hint⟩
  rw [transportedAction_closed_variation q h τ hh N A φ ρ J w _ hj, hclock, hint,
    h_mul_electricFieldScaled h hh]
  ring

/-- **Round-trip criterion.**  The round trip lowers the transported action
exactly when `-8 < q h E_{j+1}(e)`. -/
theorem roundTrip_lowers_iff (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hu : w.port j = seamLeft e) (h1 : w.steps j = .rest) (h2 : w.steps (j + 1) = .rest) :
    transportedAction q h τ N A φ ρ J (roundTripVariation w j e hu h1 h2).worldline <
        transportedAction q h τ N A φ ρ J w ↔
      -8 < q * (h * electricFieldScaled h A φ (j + 1) e) := by
  rw [← sub_neg, (roundTrip_action_difference q h τ hh N A φ ρ J w j hj e hu h1 h2).1]
  constructor <;> intro hlt <;> linarith

/-- **Rest is not a local minimum at zero seam field.**  A worldline resting
for two steps at the smaller endpoint of `e` inside the window, with the
scaled electric field on `e` vanishing at the step between, is lowered by
`8` under the round trip; it is neither single-step stationary nor a
single-step local minimum. -/
theorem rest_not_localMin_zero_field (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hu : w.port j = seamLeft e) (h1 : w.steps j = .rest) (h2 : w.steps (j + 1) = .rest)
    (hE : electricFieldScaled h A φ (j + 1) e = 0) :
    transportedAction q h τ N A φ ρ J (roundTripVariation w j e hu h1 h2).worldline =
        transportedAction q h τ N A φ ρ J w - 8 ∧
    ¬ SingleStepLocalMin q h τ N A φ ρ J w ∧
    ¬ SingleStepStationary q h τ N A φ ρ J w := by
  have hd := (roundTrip_action_difference q h τ hh N A φ ρ J w j hj e hu h1 h2).1
  rw [hE, mul_zero, mul_zero, sub_zero] at hd
  refine ⟨by linarith, fun hmin ↦ ?_, fun hst ↦ ?_⟩
  · have := hmin j hj (roundTripVariation w j e hu h1 h2)
    linarith
  · have := hst j hj (roundTripVariation w j e hu h1 h2)
    linarith

/-! ### Static fields and the circulation of the seam potential -/

/-- The oriented hop indicator reverses sign under reversal of the hop. -/
theorem hopIndicator_symm (u v : Fin 12) : hopIndicator v u = -hopIndicator u v := by
  funext e
  have hlt := (seam_table_sound e).1
  simp only [hopIndicator, Pi.neg_apply]
  by_cases h1 : seamLeft e = v ∧ seamRight e = u
  · have h3 : ¬ (seamLeft e = u ∧ seamRight e = v) := fun h3 ↦
      hlt.ne (h1.1.trans h3.2.symm)
    rw [if_pos h1, if_neg h3, if_pos ⟨h1.2, h1.1⟩]
    norm_num
  · by_cases h2 : seamRight e = v ∧ seamLeft e = u
    · have h4 : ¬ (seamRight e = u ∧ seamLeft e = v) := fun h4 ↦
        hlt.ne (h4.2.trans h2.1.symm)
      rw [if_neg h1, if_pos h2, if_pos ⟨h2.2, h2.1⟩]
    · have h3 : ¬ (seamLeft e = u ∧ seamRight e = v) := fun h3 ↦ h2 ⟨h3.2, h3.1⟩
      have h4 : ¬ (seamRight e = u ∧ seamLeft e = v) := fun h4 ↦ h1 ⟨h4.2, h4.1⟩
      rw [if_neg h1, if_neg h2, if_neg h3, if_neg h4]
      norm_num

/-- The circulation of a static seam potential around the port loop
`u → p → x → p' → u`: the sum of the pairings of the four oriented hop
indicators against the potential. -/
def circulation (A₀ : Fin 30 → ℝ) (u p x p' : Fin 12) : ℝ :=
  realSeamInner (hopIndicator u p) A₀ + realSeamInner (hopIndicator p x) A₀ +
    realSeamInner (hopIndicator x p') A₀ + realSeamInner (hopIndicator p' u) A₀

/-- The circulation is the pairing along the route `u → p → x` minus the
pairing along the route `u → p' → x`. -/
theorem circulation_eq_route_difference (A₀ : Fin 30 → ℝ) (u p x p' : Fin 12) :
    circulation A₀ u p x p' =
      (realSeamInner (hopIndicator u p) A₀ + realSeamInner (hopIndicator p x) A₀) -
        (realSeamInner (hopIndicator u p') A₀ + realSeamInner (hopIndicator p' x) A₀) := by
  unfold circulation
  rw [hopIndicator_symm x p', hopIndicator_symm p' u, realSeamInner_comm (-_),
    seamInner_neg_right, realSeamInner_comm (-_), seamInner_neg_right,
    realSeamInner_comm A₀, realSeamInner_comm A₀]
  ring

/-- **Static interaction difference.**  With `A` and `φ` constant in the
step index, the interaction difference of a closed two-step variation is
`q` times `h` times the port-potential difference at the two intermediate
ports plus `q` times the circulation of `A` around the loop through the
original intermediate port and back through the alternative one. -/
theorem static_interaction_difference (q h : ℝ) (A₀ : Fin 30 → ℝ) (φ₀ : Fin 12 → ℝ)
    (w : SeamStepWorldline) {j : ℕ} (v : ClosedTwoStepVariation w j) :
    interactionDifference q h (fun _ ↦ A₀) (fun _ ↦ φ₀) w v =
      q * (h * (φ₀ v.mid - φ₀ (w.port (j + 1))) +
        circulation A₀ (w.port j) (w.port (j + 1)) (w.port (j + 2)) v.mid) := by
  unfold interactionDifference routePairing potentialPairing
  rw [circulation_eq_route_difference, realSeamInner_comm _ A₀, realSeamInner_comm _ A₀,
    realSeamInner_comm _ A₀, realSeamInner_comm _ A₀]
  ring

/-- The two-step round trip across `e` from its smaller endpoint has zero
circulation: forward then backward on one seam. -/
theorem roundTrip_circulation_zero (A₀ : Fin 30 → ℝ) (e : Fin 30) :
    circulation A₀ (seamLeft e) (seamLeft e) (seamLeft e) (seamRight e) = 0 := by
  unfold circulation
  rw [hopIndicator_self, hopIndicator_forward, hopIndicator_backward, realSeamInner_zero_left,
    realSeamInner_comm, seamInner_single, realSeamInner_comm, seamInner_single]
  ring

/-- **Static round trip.**  With static potentials the round trip across `e`
from two rests at its smaller endpoint moves the interaction by
`q h (φ₀ (right) - φ₀ (left))` and the transported action by `-8` plus
that amount. -/
theorem static_roundTrip_difference (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A₀ : Fin 30 → ℝ) (φ₀ : Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hu : w.port j = seamLeft e) (h1 : w.steps j = .rest) (h2 : w.steps (j + 1) = .rest) :
    interactionDifference q h (fun _ ↦ A₀) (fun _ ↦ φ₀) w (roundTripVariation w j e hu h1 h2) =
      q * (h * (φ₀ (seamRight e) - φ₀ (seamLeft e))) ∧
    transportedAction q h τ N (fun _ ↦ A₀) (fun _ ↦ φ₀) ρ J
        (roundTripVariation w j e hu h1 h2).worldline -
      transportedAction q h τ N (fun _ ↦ A₀) (fun _ ↦ φ₀) ρ J w =
      -8 + q * (h * (φ₀ (seamRight e) - φ₀ (seamLeft e))) := by
  have hint := (roundTrip_action_difference q h τ hh N (fun _ ↦ A₀) (fun _ ↦ φ₀) ρ J w j hj e
    hu h1 h2).2.2
  simp only [sub_self, zero_add] at hint
  refine ⟨hint, ?_⟩
  rw [transportedAction_closed_variation q h τ hh N _ _ ρ J w _ hj,
    (roundTrip_action_difference q h τ hh N (fun _ ↦ A₀) (fun _ ↦ φ₀) ρ J w j hj e
      hu h1 h2).2.1, hint]

/-- **Rest is interaction-stationary on an equipotential seam.**  With static
potentials and equal port potential at the two endpoints of `e`, the
round trip leaves the interaction fixed and lowers the transported action
by exactly the clock term `8`. -/
theorem static_rest_interaction_stationary (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A₀ : Fin 30 → ℝ) (φ₀ : Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hu : w.port j = seamLeft e) (h1 : w.steps j = .rest) (h2 : w.steps (j + 1) = .rest)
    (heq : φ₀ (seamRight e) = φ₀ (seamLeft e)) :
    interactionDifference q h (fun _ ↦ A₀) (fun _ ↦ φ₀) w (roundTripVariation w j e hu h1 h2) =
      0 ∧
    transportedAction q h τ N (fun _ ↦ A₀) (fun _ ↦ φ₀) ρ J
        (roundTripVariation w j e hu h1 h2).worldline =
      transportedAction q h τ N (fun _ ↦ A₀) (fun _ ↦ φ₀) ρ J w - 8 := by
  have hd := static_roundTrip_difference q h τ hh N A₀ φ₀ ρ J w j hj e hu h1 h2
  rw [heq, sub_self, mul_zero, mul_zero] at hd
  exact ⟨hd.1, by linarith [hd.2]⟩

/-! ### Non-forcing -/

/-- **Two charges, two balances.**  Distinct charges give distinct exchange
differences whenever the seam field at the delayed crossing is nonzero. -/
theorem two_charges_two_balances (q₁ q₂ h τ : ℝ) (hq : q₁ ≠ q₂) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hf : w.steps j = .forward e) (hr : w.steps (j + 1) = .rest)
    (hE : electricFieldScaled h A φ (j + 1) e ≠ 0) :
    transportedAction q₁ h τ N A φ ρ J (exchangeVariation w j e hf hr).worldline -
        transportedAction q₁ h τ N A φ ρ J w ≠
      transportedAction q₂ h τ N A φ ρ J (exchangeVariation w j e hf hr).worldline -
        transportedAction q₂ h τ N A φ ρ J w := by
  rw [(exchange_action_difference q₁ h τ hh N A φ ρ J w j hj e hf hr).1,
    (exchange_action_difference q₂ h τ hh N A φ ρ J w j hj e hf hr).1]
  intro heq
  exact hq (mul_right_cancel₀ (mul_ne_zero hh hE) heq)

/-- **The unit drops out.**  Every closed two-step action difference is the
same at any two declared units. -/
theorem difference_forgets_unit (q h τ₁ τ₂ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ClosedTwoStepVariation w j) (hj : j + 1 < N + 1) :
    transportedAction q h τ₁ N A φ ρ J v.worldline - transportedAction q h τ₁ N A φ ρ J w =
      transportedAction q h τ₂ N A φ ρ J v.worldline -
        transportedAction q h τ₂ N A φ ρ J w := by
  rw [transportedAction_closed_variation q h τ₁ hh N A φ ρ J w v hj,
    transportedAction_closed_variation q h τ₂ hh N A φ ρ J w v hj]

end

end OPH.TransportedChargeForceLaw

/- Axiom audit: standard axioms only (`propext`, `Classical.choice`,
`Quot.sound`); no `sorry`, no `native_decide`, no project axiom. -/

#print axioms OPH.TransportedChargeForceLaw.ClosedTwoStepVariation.varied_port_eq
#print axioms OPH.TransportedChargeForceLaw.ClosedTwoStepVariation.varied_port_mid
#print axioms OPH.TransportedChargeForceLaw.interactionA_eq_sum_hopPairing
#print axioms OPH.TransportedChargeForceLaw.clockAction_closed_variation
#print axioms OPH.TransportedChargeForceLaw.interactionA_closed_variation
#print axioms OPH.TransportedChargeForceLaw.transportedAction_closed_variation
#print axioms OPH.TransportedChargeForceLaw.exchange_action_difference
#print axioms OPH.TransportedChargeForceLaw.exchange_balance
#print axioms OPH.TransportedChargeForceLaw.exchange_stationary_iff
#print axioms OPH.TransportedChargeForceLaw.exchange_backward_action_difference
#print axioms OPH.TransportedChargeForceLaw.roundTrip_action_difference
#print axioms OPH.TransportedChargeForceLaw.roundTrip_lowers_iff
#print axioms OPH.TransportedChargeForceLaw.rest_not_localMin_zero_field
#print axioms OPH.TransportedChargeForceLaw.hopIndicator_symm
#print axioms OPH.TransportedChargeForceLaw.static_interaction_difference
#print axioms OPH.TransportedChargeForceLaw.roundTrip_circulation_zero
#print axioms OPH.TransportedChargeForceLaw.static_roundTrip_difference
#print axioms OPH.TransportedChargeForceLaw.static_rest_interaction_stationary
#print axioms OPH.TransportedChargeForceLaw.two_charges_two_balances
#print axioms OPH.TransportedChargeForceLaw.difference_forgets_unit
