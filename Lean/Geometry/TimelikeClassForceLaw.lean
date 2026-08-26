import Geometry.TransportedChargeForceLaw
import Geometry.SeamStepSpeedLimit

set_option autoImplicit false

open scoped BigOperators

namespace OPH.TimelikeClassForceLaw

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
open OPH.TransportedChargeForceLaw
open OPH.SeamStepSpeedLimit
open OPH.C1Lorentz (Spatial Herm2 spatialDot spatialNormSq lorentzB lorentzQ)

/-!
# The force law inside the timelike class and the exchange as momentum timing

STATUS.  Candidate module on the coupled-action row.  The closed two-step
variations of `Geometry/TransportedChargeForceLaw.lean` keep every port
other than one.  This module declares two further variation classes of a
seam-step worldline and computes the exact transported-action difference
in each: the one-step variations with shifted ports (one step replaced,
every later port shifted, the later steps kept and assumed admissible at
the shifted ports), and the in-block variations of the rest-diluted class
(the crossing of a block of `k` rests and one crossing moved by one
position).  The in-block variations are closed two-step variations, so the
force law of the closed class restricts to the rest-diluted class exactly.
The exchange is then read on the generated path as a change of timing of
one fixed momentum quantum, and on the field side as a change of timing of
one fixed energy transfer, with the exact factor `h` between the action
difference and the transferred energy.  Every statement is an identity, an
equivalence, or an inequation on the finite complex; the potentials are
held fixed in every action difference.

WHAT IS PROVED.

1. One-step variations with shifted ports.  `replaceOne` replaces step `j`.
   A shifted one-step variation (`ShiftedOneStepVariation`) is an
   admissible alternative step at `w.port j` together with the declared
   hypothesis that every later original step is admissible at the port it
   is shifted to (`admAfter`).  The varied worldline is a seam-step
   worldline (`ShiftedOneStepVariation.worldline`); its ports agree with
   the original ones up to `j` (`shifted_port_le`) and its port at `j + 1`
   is the alternative target (`shifted_port_succ`).  The exact action
   difference over a window containing `j` is the one-step clock
   difference `|s_j|² - |s'|²` plus `q` times the pairing difference at
   step `j` plus `q` times the exact sum over the later window steps of
   the potential pairings at the shifted ports minus those at the original
   ports (`transportedAction_shifted_variation`); nothing is assumed about
   the potentials.  Under the declared shift invariance of the potentials
   (`ShiftInvariant`: every later pairing along the shifted route equals
   the pairing along the original route) the sum vanishes and the
   difference is the one-step term (`shifted_variation_invariant`); for a
   forward crossing replaced by a rest this is `4 + q A (j+1) e`
   (`crossing_to_rest_invariant`), which carries the seam potential itself
   whereas the closed two-step differences carry only field pairings
   (their gauge behavior is stated nowhere in this module).  When
   the shifted route rejoins the original port at `j + 2`
   (`RejoinsAtTwo`), every later port agrees (`shifted_port_rejoin`) and
   the difference is the closed two-step difference of the pair
   `(s', s_{j+1})`, in the form of `transportedAction_closed_variation`
   (`shifted_variation_rejoin`).
2. Rest-diluted class and in-block variations.  `BlockPattern w n k i`
   declares a block of `k + 1` steps from `n` whose only crossing sits at
   position `i`; `RestBlock w n k` is the pattern with `i = k`
   (`restBlock_iff_pattern`); `RestDiluted w k` declares the class of
   worldlines whose every block of length `k + 1` is a rest block.  Delaying
   a forward crossing at position `i < k` by one position is the exchange
   variation of the closed class; the varied worldline carries the pattern
   at position `i + 1` (`delay_pattern`) and every other block is untouched
   (`delay_steps_outside`).  Its action difference is `q h E_{n+i+1}(e)`
   with zero clock difference (`delay_action_difference`); advancing a
   crossing at position `0 < i` gives `-q h E_{n+i}(e)`
   (`advance_action_difference`), and the backward crossing gives the
   opposite signs (`advance_backward_action_difference` here; the delayed
   backward crossing is `TransportedChargeForceLaw.exchange_backward_action_difference`).  `PositionStationary` declares
   stationarity of the crossing position inside its block (no in-block move
   changes the action); with `q ≠ 0`, `h ≠ 0` it holds exactly when the
   scaled seam field on the crossed seam vanishes at the delayed index and
   at the original index (`position_stationary_iff`).  Delay lowers the
   action exactly when `q h E_{n+i+1}(e) < 0`, advance exactly when
   `0 < q h E_{n+i}(e)` (`delay_lowers_iff`, `advance_lowers_iff`).  With a
   static nonzero field on the seam an interior crossing position is never
   stationary and one of the two moves lowers the action, the direction
   fixed by the sign of `q h E₀` (`static_field_drift`): the forward
   crossing moves earlier when `0 < q h E₀`, later when `q h E₀ < 0`.
3. Momentum timing.  The generated path of a closed two-step variation
   agrees with the original at every index other than `j + 1`
   (`generatedPath_closed_variation`), so the discrete momentum over any
   interval avoiding `j + 1` is unchanged (`intervalMomentum_closed_variation`),
   in particular the block momentum `2m ((k+1) τ, s_e)` of every block
   (`blockMomentum_eq`, `blockMomentum_delay`).  The exchange swaps the two
   per-step momenta: the crossing quantum `2m (τ, s_e)` moves from step `j`
   to step `j + 1` and the rest quantum `2m (τ, 0)` the other way
   (`delay_momentum_swap`).  On the field side, along the scaled Ampere
   evolution sourced by the respective hopping current, the field energy
   moves by `(q/2)(E_n e + E_{n+1} e)` at the hop step and by zero at the
   rest step (`hopping_work_energy_forward`, `hopping_work_energy_rest`);
   with the seam field static across the three indices involved the
   transfer is `q E₀` at the hop step and zero at the rest step for both
   timings (`transfer_original`, `transfer_delayed`), and the action
   difference of the delay equals `h` times the transfer at the original
   hop step, equivalently `-h` times the difference of the mid-window
   energy increments between the two timings
   (`exchange_action_eq_h_mul_transfer`).  The factor is exactly `h`: the
   interaction weights are `h ⟨J, A⟩ + h ⟨ρ, φ⟩` while the energy balance
   pairs the current `-(q/h)` against `(h/2)(E_n + E_{n+1})`, so `h`
   cancels on the energy side and survives on the action side.
4. Non-forcing.  Two charges give two in-block balances at a nonzero seam
   field (`two_charges_two_block_balances`); every in-block difference and
   every shifted one-step difference is the same at any two declared units
   (`block_difference_forgets_unit`, `shifted_difference_forgets_unit`).
   The variation classes, the block pattern, the unit `τ`, and the step `h`
   are declared.

PRIOR WORK.  `TransportedChargeForceLaw.transportedAction_closed_variation`,
`exchange_action_difference`, `exchange_backward_action_difference`,
`exchange_stationary_iff`, `two_charges_two_balances`, and
`difference_forgets_unit` are the closed-class results restricted here;
the in-block delay is literally `exchangeVariation`.  `SeamStepSpeedLimit.
RestBlock`, `block_displacement`, `block_lorentzQ`, `block_timelike_iff`
supply the block pattern and its kinematics; the momentum statements use
`block_displacement` through `ray_port_eq_spatial`.  `ChargeFixedInteraction.
discreteMomentum` is the momentum of the weighted clock action;
`coupled_momentum_law` is its balance on the E-paired route and is not
restated.  `PortChargeMinimalCoupling.hopping_work_energy` and its
forward/rest cases supply the field-side balance.  No prior module treats
shifted one-step variations or the in-block position variable.

ROWS TOUCHED.  The coupled-action row (the two variation classes and the
in-block balance are declared here); the source clock and duration row
(the block pattern, the rest count `k`, and the unit `τ` are declared; `τ`
drops out of every difference); the physical spacetime attachment row (no
port is attached to a spacetime point); the light-signal row (no signal is
attached); the laboratory clock and energy calibration import (no unit is
attached to `q`, `h`, `m`, or the seam field); the gravitation-route
energy identification (the clock term is not identified with an energy).
None of these rows is discharged.

NEGATIVES CITED.  The Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`): the clock action
shape, the coupling shape, and both variation classes are declared
enrichments; cited at scope only.

CONVENTIONS.  Signature `(+---)`; a crossing step of the generated path has
Lorentz square `τ² - 4`, a rest `τ²`.  Seam orientation from `seamLeft e`
to `seamRight e`; forward hopping current `-(q / h)` on `e`; scaled field
`E n = -(A (n+1) - A n) / h - d (φ n)`.  Block positions are counted from
the block start; differences are written varied minus original; a delay
moves the crossing to a later index.

FALSIFIER.  The module is wrong if some in-block move carries a clock
term, if the shifted-port sum fails to vanish under shift invariance, if
the block momentum changes under a delay, or if the action-to-transfer
factor differs from `h`.

Axiom audit.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`; no `sorry`, no
`native_decide`, no project axiom.
-/

/-! ## (1) One-step variations with shifted ports -/

/-- Replace step `j` of a step sequence. -/
def replaceOne (steps : ℕ → SeamStep) (j : ℕ) (s' : SeamStep) : ℕ → SeamStep :=
  fun n ↦ if n = j then s' else steps n

theorem replaceOne_at (steps : ℕ → SeamStep) (j : ℕ) (s' : SeamStep) :
    replaceOne steps j s' j = s' := by
  simp [replaceOne]

theorem replaceOne_other (steps : ℕ → SeamStep) (j : ℕ) (s' : SeamStep) (n : ℕ)
    (hn : n ≠ j) : replaceOne steps j s' n = steps n := by
  simp [replaceOne, hn]

/-- The shifted port sequence of a one-step replacement. -/
def shiftedPort (w : SeamStepWorldline) (j : ℕ) (s' : SeamStep) : ℕ → Fin 12 :=
  portSeq w.start (replaceOne w.steps j s')

theorem shiftedPort_succ (w : SeamStepWorldline) (j : ℕ) (s' : SeamStep) (n : ℕ) :
    shiftedPort w j s' (n + 1) = stepTarget (shiftedPort w j s' n) (replaceOne w.steps j s' n) :=
  rfl

/-- Up to `j` the shifted ports are the original ones. -/
theorem shiftedPort_le (w : SeamStepWorldline) (j : ℕ) (s' : SeamStep) (n : ℕ) (hn : n ≤ j) :
    shiftedPort w j s' n = w.port n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    rw [shiftedPort_succ, ih (by omega), w.port_succ, replaceOne_other _ _ _ _ (by omega)]

/-- At `j + 1` the shifted port is the target of the alternative step. -/
theorem shiftedPort_at_succ (w : SeamStepWorldline) (j : ℕ) (s' : SeamStep) :
    shiftedPort w j s' (j + 1) = stepTarget (w.port j) s' := by
  rw [shiftedPort_succ, shiftedPort_le w j s' j le_rfl, replaceOne_at]

/-- Declared shifted one-step variation at `j`: an admissible alternative
step at `w.port j`, with every later original step admissible at the port
it is shifted to.  The admissibility-after-shift clause is a hypothesis of
the class, declared, and is what makes the varied sequence a worldline. -/
structure ShiftedOneStepVariation (w : SeamStepWorldline) (j : ℕ) where
  /-- The alternative step. -/
  s' : SeamStep
  /-- The alternative step is admissible at the port reached at `j`. -/
  adm : StepAdmissible (w.port j) s'
  /-- Every later original step is admissible at its shifted port. -/
  admAfter : ∀ n, j < n → StepAdmissible (shiftedPort w j s' n) (w.steps n)

namespace ShiftedOneStepVariation

variable {w : SeamStepWorldline} {j : ℕ} (v : ShiftedOneStepVariation w j)

/-- The varied worldline. -/
def worldline : SeamStepWorldline where
  start := w.start
  steps := replaceOne w.steps j v.s'
  adm := by
    intro n
    show StepAdmissible (shiftedPort w j v.s' n) (replaceOne w.steps j v.s' n)
    by_cases hn : n = j
    · subst hn
      rw [replaceOne_at, shiftedPort_le w n v.s' n le_rfl]
      exact v.adm
    · rw [replaceOne_other _ _ _ _ hn]
      rcases Nat.lt_or_gt_of_ne hn with h | h
      · rw [shiftedPort_le w j v.s' n h.le]
        exact w.adm n
      · exact v.admAfter n h

theorem worldline_steps (n : ℕ) : v.worldline.steps n = replaceOne w.steps j v.s' n := rfl

theorem worldline_port (n : ℕ) : v.worldline.port n = shiftedPort w j v.s' n := rfl

/-- **Shifted ports below the varied step** agree with the original ones. -/
theorem shifted_port_le (n : ℕ) (hn : n ≤ j) : v.worldline.port n = w.port n :=
  shiftedPort_le w j v.s' n hn

/-- **The shifted port at `j + 1`** is the target of the alternative step. -/
theorem shifted_port_succ : v.worldline.port (j + 1) = stepTarget (w.port j) v.s' :=
  shiftedPort_at_succ w j v.s'

end ShiftedOneStepVariation

noncomputable section

/-- The clock difference of a shifted one-step variation: the original step
norm minus the alternative one; the later steps are kept, so they cancel. -/
def oneStepClockDifference (w : SeamStepWorldline) {j : ℕ} (v : ShiftedOneStepVariation w j) :
    ℝ :=
  stepNormSq (w.steps j) - stepNormSq v.s'

/-- The pairing difference at the varied step itself. -/
def oneStepPairingDifference (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (w : SeamStepWorldline) {j : ℕ} (v : ShiftedOneStepVariation w j) : ℝ :=
  potentialPairing h A φ j (w.port j) (stepTarget (w.port j) v.s') -
    potentialPairing h A φ j (w.port j) (w.port (j + 1))

/-- The exact later-port sum: over the window steps after `j`, the pairing
at the shifted ports minus the pairing at the original ports. -/
def shiftedTail (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (N : ℕ)
    (w : SeamStepWorldline) {j : ℕ} (v : ShiftedOneStepVariation w j) : ℝ :=
  ∑ n ∈ Finset.Ico (j + 1) (N + 1),
    (potentialPairing h A φ n (v.worldline.port n) (v.worldline.port (n + 1)) -
      potentialPairing h A φ n (w.port n) (w.port (n + 1)))

theorem clockAction_shifted_variation (τ : ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ShiftedOneStepVariation w j) (M : ℕ) (hj : j < M) :
    clockAction M (generatedPath τ v.worldline) - clockAction M (generatedPath τ w) =
      oneStepClockDifference w v := by
  rw [clockAction_generated, clockAction_generated, ← Finset.sum_sub_distrib]
  have hpt : ∀ n, (τ ^ 2 - stepNormSq (v.worldline.steps n)) - (τ ^ 2 - stepNormSq (w.steps n)) =
      if n = j then oneStepClockDifference w v else 0 := by
    intro n
    by_cases hn : n = j
    · subst hn
      rw [if_pos rfl, v.worldline_steps, replaceOne_at]
      unfold oneStepClockDifference
      ring
    · rw [if_neg hn, v.worldline_steps, replaceOne_other _ _ _ _ hn]
      ring
  rw [Finset.sum_congr rfl fun n _ ↦ hpt n, Finset.sum_ite_eq',
    if_pos (Finset.mem_range.mpr hj)]

theorem interactionA_shifted_variation (q h : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ShiftedOneStepVariation w j) (hj : j < N + 1) :
    interactionA q h N A φ (hoppingPath v.worldline) - interactionA q h N A φ (hoppingPath w) =
      q * oneStepPairingDifference h A φ w v + q * shiftedTail h A φ N w v := by
  rw [interactionA_eq_sum_hopPairing q h hh, interactionA_eq_sum_hopPairing q h hh,
    ← Finset.sum_sub_distrib, ← Finset.sum_range_add_sum_Ico _ (show j + 1 ≤ N + 1 by omega),
    Finset.sum_range_succ]
  have hlow : ∑ n ∈ Finset.range j,
      (hopPairing q h A φ n ((hoppingPath v.worldline).port n)
          ((hoppingPath v.worldline).port (n + 1)) -
        hopPairing q h A φ n ((hoppingPath w).port n) ((hoppingPath w).port (n + 1))) = 0 := by
    refine Finset.sum_eq_zero fun n hn ↦ ?_
    have hn' := Finset.mem_range.mp hn
    simp only [hoppingPath_port]
    rw [v.shifted_port_le n hn'.le, v.shifted_port_le (n + 1) (by omega), sub_self]
  rw [hlow, zero_add]
  simp only [hoppingPath_port]
  rw [v.shifted_port_le j le_rfl, v.shifted_port_succ]
  unfold shiftedTail oneStepPairingDifference hopPairing
  rw [Finset.mul_sum]
  congr 1
  · ring
  · refine Finset.sum_congr rfl fun n _ ↦ ?_
    ring

/-- **Exact action difference under a shifted one-step variation.**  With
the potentials held fixed and the varied step inside the window, the
transported action moves by the one-step clock difference, plus `q` times
the pairing difference at the varied step, plus `q` times the exact sum of
the later pairing differences at the shifted ports.  Nothing is assumed
about the potentials; the committed window action cancels identically and
the unit drops out. -/
theorem transportedAction_shifted_variation (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ShiftedOneStepVariation w j) (hj : j < N + 1) :
    transportedAction q h τ N A φ ρ J v.worldline - transportedAction q h τ N A φ ρ J w =
      oneStepClockDifference w v + q * oneStepPairingDifference h A φ w v +
        q * shiftedTail h A φ N w v := by
  have hc := clockAction_shifted_variation τ w v (N + 1) hj
  have hi := interactionA_shifted_variation q h hh N A φ w v hj
  unfold transportedAction monopoleCoupledAction
  linarith

end

/-! ### Shift invariance and the rejoining case -/

noncomputable section

/-- Declared shift invariance of the potentials for a shifted one-step
variation: at every later window step the potential pairing along the
shifted route equals the pairing along the original route.  This is a
hypothesis on `A`, `φ`, and the route; it is not derived. -/
def ShiftInvariant (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (N : ℕ)
    (w : SeamStepWorldline) {j : ℕ} (v : ShiftedOneStepVariation w j) : Prop :=
  ∀ n, j + 1 ≤ n → n < N + 1 →
    potentialPairing h A φ n (v.worldline.port n) (v.worldline.port (n + 1)) =
      potentialPairing h A φ n (w.port n) (w.port (n + 1))

theorem shiftedTail_of_invariant (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (N : ℕ)
    (w : SeamStepWorldline) {j : ℕ} (v : ShiftedOneStepVariation w j)
    (hinv : ShiftInvariant h A φ N w v) : shiftedTail h A φ N w v = 0 := by
  unfold shiftedTail
  refine Finset.sum_eq_zero fun n hn ↦ ?_
  have hn' := Finset.mem_Ico.mp hn
  rw [hinv n hn'.1 hn'.2, sub_self]

/-- **Reduction under shift invariance.**  With shift-invariant potentials
the shifted one-step difference is the one-step term: the clock difference
plus `q` times the pairing difference at the varied step. -/
theorem shifted_variation_invariant (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ShiftedOneStepVariation w j) (hj : j < N + 1)
    (hinv : ShiftInvariant h A φ N w v) :
    transportedAction q h τ N A φ ρ J v.worldline - transportedAction q h τ N A φ ρ J w =
      oneStepClockDifference w v + q * oneStepPairingDifference h A φ w v := by
  rw [transportedAction_shifted_variation q h τ hh N A φ ρ J w v hj,
    shiftedTail_of_invariant h A φ N w v hinv, mul_zero, add_zero]

/-- **A crossing replaced by a rest, shift-invariant potentials.**  The
difference is `4 + q A (j+1) e`: the clock gain of one crossing and the
seam potential at the node after the removed crossing.  The seam potential
appears bare, so this one-step difference is not gauge covariant; the
closed class of `TransportedChargeForceLaw` is the covariant one. -/
theorem crossing_to_rest_invariant (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ShiftedOneStepVariation w j) (hj : j < N + 1) (e : Fin 30)
    (hf : w.steps j = .forward e) (hs : v.s' = .rest)
    (hinv : ShiftInvariant h A φ N w v) :
    transportedAction q h τ N A φ ρ J v.worldline - transportedAction q h τ N A φ ρ J w =
      4 + q * A (j + 1) e := by
  rw [shifted_variation_invariant q h τ hh N A φ ρ J w v hj hinv]
  have hpj := port_of_forward w j e hf
  have hpj1 := port_succ_of_forward w j e hf
  unfold oneStepClockDifference oneStepPairingDifference
  rw [hf, hs, hpj1, hpj]
  show (seamNormSq - 0) + q * (potentialPairing h A φ j (seamLeft e) (seamLeft e) -
    potentialPairing h A φ j (seamLeft e) (seamRight e)) = 4 + q * A (j + 1) e
  rw [potentialPairing_rest, potentialPairing_forward]
  unfold seamNormSq
  ring

/-- Declared rejoining: the shifted route occupies the original port at
`j + 2`. -/
def RejoinsAtTwo (w : SeamStepWorldline) {j : ℕ} (v : ShiftedOneStepVariation w j) : Prop :=
  v.worldline.port (j + 2) = w.port (j + 2)

/-- Beyond `j + 2` a rejoining shifted route occupies the original ports. -/
theorem shifted_port_rejoin (w : SeamStepWorldline) {j : ℕ} (v : ShiftedOneStepVariation w j)
    (hr : RejoinsAtTwo w v) (n : ℕ) (hn : j + 2 ≤ n) : v.worldline.port n = w.port n := by
  obtain ⟨k, rfl⟩ : ∃ k, n = j + 2 + k := ⟨n - (j + 2), by omega⟩
  induction k with
  | zero => exact hr
  | succ k ih =>
    rw [show j + 2 + (k + 1) = (j + 2 + k) + 1 by omega, v.worldline.port_succ,
      ih (by omega), w.port_succ, v.worldline_steps, replaceOne_other _ _ _ _ (by omega)]

/-- **Rejoining recovers the closed two-step difference.**  When the shifted
route rejoins at `j + 2` and both steps lie in the window, the difference is
the clock difference of the varied step plus `q` times the route pairing of
`u → stepTarget u s' → w.port (j+2)` minus that of the original route: the
form of `transportedAction_closed_variation` for the pair
`(s', w.steps (j+1))`. -/
theorem shifted_variation_rejoin (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ShiftedOneStepVariation w j) (hj : j + 1 < N + 1) (hr : RejoinsAtTwo w v) :
    transportedAction q h τ N A φ ρ J v.worldline - transportedAction q h τ N A φ ρ J w =
      oneStepClockDifference w v +
        q * (routePairing h A φ j (w.port j) (stepTarget (w.port j) v.s') (w.port (j + 2)) -
          routePairing h A φ j (w.port j) (w.port (j + 1)) (w.port (j + 2))) := by
  rw [transportedAction_shifted_variation q h τ hh N A φ ρ J w v (by omega)]
  have htail : shiftedTail h A φ N w v =
      potentialPairing h A φ (j + 1) (stepTarget (w.port j) v.s') (w.port (j + 2)) -
        potentialPairing h A φ (j + 1) (w.port (j + 1)) (w.port (j + 2)) := by
    unfold shiftedTail
    rw [Finset.sum_eq_sum_Ico_succ_bot hj, v.shifted_port_succ, hr]
    have hzero : ∑ n ∈ Finset.Ico (j + 1 + 1) (N + 1),
        (potentialPairing h A φ n (v.worldline.port n) (v.worldline.port (n + 1)) -
          potentialPairing h A φ n (w.port n) (w.port (n + 1))) = 0 := by
      refine Finset.sum_eq_zero fun n hn ↦ ?_
      have hn' := Finset.mem_Ico.mp hn
      rw [shifted_port_rejoin w v hr n (by omega), shifted_port_rejoin w v hr (n + 1) (by omega),
        sub_self]
    rw [hzero, add_zero]
  rw [htail]
  unfold oneStepPairingDifference routePairing
  ring

end

/-! ## (2) The rest-diluted class and its in-block variations -/

/-- Declared block pattern: `k + 1` steps from `n`, the only crossing at
position `i` (the step `n + i`), every other step of the block a rest. -/
def BlockPattern (w : SeamStepWorldline) (n k i : ℕ) : Prop :=
  (∀ p, p < k + 1 → p ≠ i → w.steps (n + p) = .rest) ∧ w.steps (n + i) ≠ .rest

/-- `RestBlock` is the pattern with the crossing at the last position. -/
theorem restBlock_iff_pattern (w : SeamStepWorldline) (n k : ℕ) :
    RestBlock w n k ↔ BlockPattern w n k k := by
  unfold RestBlock BlockPattern
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨fun p hp hpk ↦ h1 p (by omega), h2⟩
  · rintro ⟨h1, h2⟩
    exact ⟨fun p hp ↦ h1 p (by omega) (by omega), h2⟩

/-- Declared rest-diluted class: every block of `k + 1` steps from a
multiple of `k + 1` is a rest block, i.e. exactly `k` rests between
consecutive crossings, with the crossings at the block ends. -/
def RestDiluted (w : SeamStepWorldline) (k : ℕ) : Prop :=
  ∀ b, RestBlock w (b * (k + 1)) k

/-- The two-rest worldline of `SeamStepSpeedLimit` is rest diluted with
`k = 2`. -/
theorem twoRestWorldline_restDiluted : RestDiluted twoRestWorldline 2 := by
  intro b
  have := twoRestWorldline_block b
  rwa [show 3 * b = b * (2 + 1) by ring] at this

/-- In a block pattern the step after the crossing is a rest when the
crossing is not last. -/
theorem pattern_rest_after (w : SeamStepWorldline) (n k i : ℕ) (hb : BlockPattern w n k i)
    (hi : i < k) : w.steps (n + i + 1) = .rest :=
  hb.1 (i + 1) (by omega) (by omega)

/-- In a block pattern the step before the crossing is a rest when the
crossing is not first. -/
theorem pattern_rest_before (w : SeamStepWorldline) (n k i : ℕ) (hb : BlockPattern w n k (i + 1))
    (hi : i + 1 < k + 1) : w.steps (n + i) = .rest :=
  hb.1 i (by omega) (by omega)

/-! ### The advance variation of the closed class -/

/-- The advance variation: a rest at `j` followed by a forward crossing of
`e` is replaced by the crossing followed by a rest. -/
def advanceVariation (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hr : w.steps j = .rest) (hf : w.steps (j + 1) = .forward e) :
    ClosedTwoStepVariation w j where
  s' := .forward e
  s'' := .rest
  adm1 := by
    have h1 := port_of_forward w (j + 1) e hf
    rw [port_succ_of_rest w j hr] at h1
    exact h1
  adm2 := trivial
  closed := by
    show seamRight e = w.port (j + 2)
    rw [port_succ_of_forward w (j + 1) e hf]

/-- The backward advance variation. -/
def advanceBackwardVariation (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hr : w.steps j = .rest) (hb : w.steps (j + 1) = .backward e) :
    ClosedTwoStepVariation w j where
  s' := .backward e
  s'' := .rest
  adm1 := by
    have h1 := port_of_backward w (j + 1) e hb
    rw [port_succ_of_rest w j hr] at h1
    exact h1
  adm2 := trivial
  closed := by
    show seamLeft e = w.port (j + 2)
    rw [port_succ_of_backward w (j + 1) e hb]

noncomputable section

/-- **Advance difference.**  Advancing a forward crossing of `e` at `j + 1`
by one rest moves the transported action by `-q h E_{j+1}(e)`; the clock
difference is zero. -/
theorem advance_action_difference (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hr : w.steps j = .rest) (hf : w.steps (j + 1) = .forward e) :
    transportedAction q h τ N A φ ρ J (advanceVariation w j e hr hf).worldline -
        transportedAction q h τ N A φ ρ J w =
      -(q * (h * electricFieldScaled h A φ (j + 1) e)) ∧
    clockDifference w (advanceVariation w j e hr hf) = 0 := by
  have hpj1 := port_of_forward w (j + 1) e hf
  have hpj : w.port j = seamLeft e := by rw [← port_succ_of_rest w j hr, hpj1]
  have hpj2 := port_succ_of_forward w (j + 1) e hf
  have hmid : (advanceVariation w j e hr hf).mid = seamRight e := by
    show stepTarget (w.port j) (.forward e) = seamRight e
    rfl
  have hclock : clockDifference w (advanceVariation w j e hr hf) = 0 := by
    unfold clockDifference
    rw [hr, hf]
    show (stepNormSq .rest + stepNormSq (.forward e)) -
      (stepNormSq (.forward e) + stepNormSq .rest) = 0
    ring
  refine ⟨?_, hclock⟩
  rw [transportedAction_closed_variation q h τ hh N A φ ρ J w _ hj, hclock,
    h_mul_electricFieldScaled h hh]
  unfold interactionDifference routePairing
  rw [hmid, hpj2, hpj1, hpj, potentialPairing_forward, potentialPairing_rest,
    potentialPairing_rest, potentialPairing_forward]
  ring

/-- **Backward advance difference.**  Advancing a backward crossing of `e`
at `j + 1` moves the transported action by `q h E_{j+1}(e)`. -/
theorem advance_backward_action_difference (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hr : w.steps j = .rest) (hb : w.steps (j + 1) = .backward e) :
    transportedAction q h τ N A φ ρ J (advanceBackwardVariation w j e hr hb).worldline -
        transportedAction q h τ N A φ ρ J w =
      q * (h * electricFieldScaled h A φ (j + 1) e) := by
  have hpj1 := port_of_backward w (j + 1) e hb
  have hpj : w.port j = seamRight e := by rw [← port_succ_of_rest w j hr, hpj1]
  have hpj2 := port_succ_of_backward w (j + 1) e hb
  have hmid : (advanceBackwardVariation w j e hr hb).mid = seamLeft e := by
    show stepTarget (w.port j) (.backward e) = seamLeft e
    rfl
  have hclock : clockDifference w (advanceBackwardVariation w j e hr hb) = 0 := by
    unfold clockDifference
    rw [hr, hb]
    show (stepNormSq .rest + stepNormSq (.backward e)) -
      (stepNormSq (.backward e) + stepNormSq .rest) = 0
    ring
  rw [transportedAction_closed_variation q h τ hh N A φ ρ J w _ hj, hclock,
    h_mul_electricFieldScaled h hh]
  unfold interactionDifference routePairing
  rw [hmid, hpj2, hpj1, hpj, potentialPairing_backward, potentialPairing_rest,
    potentialPairing_rest, potentialPairing_backward]
  ring

/-! ### In-block moves of the crossing -/

/-- **Delay keeps the pattern.**  Delaying the forward crossing at position
`i < k` of a block gives the pattern at position `i + 1`. -/
theorem delay_pattern (w : SeamStepWorldline) (n k i : ℕ) (hb : BlockPattern w n k i)
    (hi : i < k) (e : Fin 30) (hf : w.steps (n + i) = .forward e) :
    BlockPattern (exchangeVariation w (n + i) e hf (pattern_rest_after w n k i hb hi)).worldline
      n k (i + 1) := by
  refine ⟨fun p hp hpi ↦ ?_, ?_⟩
  · rw [ClosedTwoStepVariation.worldline_steps]
    by_cases hpi' : p = i
    · subst hpi'
      exact replaceTwo_at _ _ _ _
    · rw [replaceTwo_other _ _ _ _ _ (by omega) (by omega)]
      exact hb.1 p hp hpi'
  · rw [ClosedTwoStepVariation.worldline_steps, show n + (i + 1) = n + i + 1 by omega,
      replaceTwo_succ]
    exact SeamStep.noConfusion

/-- **Delay touches one block.**  Outside the block every step is kept. -/
theorem delay_steps_outside (w : SeamStepWorldline) (n k i : ℕ) (hb : BlockPattern w n k i)
    (hi : i < k) (e : Fin 30) (hf : w.steps (n + i) = .forward e) (m : ℕ)
    (hm : m < n ∨ n + k < m) :
    (exchangeVariation w (n + i) e hf (pattern_rest_after w n k i hb hi)).worldline.steps m =
      w.steps m := by
  rw [ClosedTwoStepVariation.worldline_steps, replaceTwo_other _ _ _ _ _ (by omega) (by omega)]

/-- **In-block delay difference.**  Delaying the forward crossing of `e` at
position `i < k` of a block from `n` moves the transported action by
`q h E_{n+i+1}(e)`, with zero clock difference. -/
theorem delay_action_difference (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (n k i : ℕ) (hb : BlockPattern w n k i)
    (hi : i < k) (hN : n + k < N + 1) (e : Fin 30) (hf : w.steps (n + i) = .forward e) :
    transportedAction q h τ N A φ ρ J
        (exchangeVariation w (n + i) e hf (pattern_rest_after w n k i hb hi)).worldline -
        transportedAction q h τ N A φ ρ J w =
      q * (h * electricFieldScaled h A φ (n + i + 1) e) ∧
    clockDifference w (exchangeVariation w (n + i) e hf (pattern_rest_after w n k i hb hi)) = 0 :=
  exchange_action_difference q h τ hh N A φ ρ J w (n + i) (by omega) e hf _

/-- **In-block advance difference.**  Advancing the forward crossing of `e`
at position `i + 1 ≤ k` of a block from `n` moves the transported action by
`-q h E_{n+i+1}(e)`, with zero clock difference. -/
theorem advance_block_action_difference (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (n k i : ℕ)
    (hb : BlockPattern w n k (i + 1)) (hi : i + 1 < k + 1) (hN : n + k < N + 1) (e : Fin 30)
    (hf : w.steps (n + i + 1) = .forward e) :
    transportedAction q h τ N A φ ρ J
        (advanceVariation w (n + i) e (pattern_rest_before w n k i hb hi) hf).worldline -
        transportedAction q h τ N A φ ρ J w =
      -(q * (h * electricFieldScaled h A φ (n + i + 1) e)) ∧
    clockDifference w (advanceVariation w (n + i) e (pattern_rest_before w n k i hb hi) hf) = 0 :=
  advance_action_difference q h τ hh N A φ ρ J w (n + i) (by omega) e _ hf

/-- **Delay lowers the action** exactly when `q h E_{j+1}(e) < 0`. -/
theorem delay_lowers_iff (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hf : w.steps j = .forward e) (hr : w.steps (j + 1) = .rest) :
    transportedAction q h τ N A φ ρ J (exchangeVariation w j e hf hr).worldline <
        transportedAction q h τ N A φ ρ J w ↔
      q * (h * electricFieldScaled h A φ (j + 1) e) < 0 := by
  rw [← sub_neg, (exchange_action_difference q h τ hh N A φ ρ J w j hj e hf hr).1]

/-- **Advance lowers the action** exactly when `0 < q h E_{j+1}(e)`. -/
theorem advance_lowers_iff (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hr : w.steps j = .rest) (hf : w.steps (j + 1) = .forward e) :
    transportedAction q h τ N A φ ρ J (advanceVariation w j e hr hf).worldline <
        transportedAction q h τ N A φ ρ J w ↔
      0 < q * (h * electricFieldScaled h A φ (j + 1) e) := by
  rw [← sub_neg, (advance_action_difference q h τ hh N A φ ρ J w j hj e hr hf).1, neg_neg_iff_pos]

end

/-! ### Stationarity of the crossing position and the drift under a static field -/

noncomputable section

/-- Declared stationarity of the crossing position `i` of the block from `n`
(crossing of `e`, forward): neither the in-block delay (when `i < k`) nor
the in-block advance (when `i = i' + 1`) changes the transported action.
The class of moves is declared; the potentials are held fixed. -/
def PositionStationary (q h τ : ℝ) (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (n k i : ℕ)
    (e : Fin 30) : Prop :=
  (∀ (_ : i < k) (hf : w.steps (n + i) = .forward e) (hr : w.steps (n + i + 1) = .rest),
    transportedAction q h τ N A φ ρ J (exchangeVariation w (n + i) e hf hr).worldline =
      transportedAction q h τ N A φ ρ J w) ∧
  (∀ i', i = i' + 1 →
    ∀ (hr : w.steps (n + i') = .rest) (hf : w.steps (n + i' + 1) = .forward e),
      transportedAction q h τ N A φ ρ J (advanceVariation w (n + i') e hr hf).worldline =
        transportedAction q h τ N A φ ρ J w)

/-- Advance stationarity: the advance leaves the action fixed exactly when
the seam field vanishes at the index of the crossing. -/
theorem advance_stationary_iff (q h τ : ℝ) (hq : q ≠ 0) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hr : w.steps j = .rest) (hf : w.steps (j + 1) = .forward e) :
    transportedAction q h τ N A φ ρ J (advanceVariation w j e hr hf).worldline =
        transportedAction q h τ N A φ ρ J w ↔
      electricFieldScaled h A φ (j + 1) e = 0 := by
  rw [← sub_eq_zero, (advance_action_difference q h τ hh N A φ ρ J w j hj e hr hf).1, neg_eq_zero]
  constructor
  · intro h0
    rcases mul_eq_zero.mp h0 with h0 | h0
    · exact absurd h0 hq
    · rcases mul_eq_zero.mp h0 with h0 | h0
      · exact absurd h0 hh
      · exact h0
  · intro h0
    rw [h0, mul_zero, mul_zero]

/-- **The timelike-class force law.**  Inside a block of the rest-diluted
pattern the forward crossing of `e` at position `i` is stationary exactly
when the scaled seam field on `e` vanishes at the delayed index `n + i + 1`
(when a delay is available) and at the crossing index `n + i` (when an
advance is available). -/
theorem position_stationary_iff (q h τ : ℝ) (hq : q ≠ 0) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (n k i : ℕ) (hb : BlockPattern w n k i)
    (hik : i ≤ k) (hN : n + k < N + 1) (e : Fin 30) (hf : w.steps (n + i) = .forward e) :
    PositionStationary q h τ N A φ ρ J w n k i e ↔
      ((i < k → electricFieldScaled h A φ (n + i + 1) e = 0) ∧
        (∀ i', i = i' + 1 → electricFieldScaled h A φ (n + i' + 1) e = 0)) := by
  constructor
  · rintro ⟨h1, h2⟩
    refine ⟨fun hi ↦ ?_, fun i' hii' ↦ ?_⟩
    · exact (exchange_stationary_iff q h τ hq hh N A φ ρ J w (n + i) (by omega) e hf
        (pattern_rest_after w n k i hb hi)).mp (h1 hi hf _)
    · subst hii'
      exact (advance_stationary_iff q h τ hq hh N A φ ρ J w (n + i') (by omega) e
        (pattern_rest_before w n k i' hb (by omega)) hf).mp (h2 i' rfl _ hf)
  · rintro ⟨h1, h2⟩
    refine ⟨fun hi hf' hr' ↦ ?_, fun i' hii' hr' hf' ↦ ?_⟩
    · exact (exchange_stationary_iff q h τ hq hh N A φ ρ J w (n + i) (by omega) e hf'
        hr').mpr (h1 hi)
    · subst hii'
      exact (advance_stationary_iff q h τ hq hh N A φ ρ J w (n + i') (by omega) e hr'
        hf').mpr (h2 i' rfl)

/-- **Drift under a static seam field.**  With the seam field on `e` equal
to a nonzero `E₀` at the crossing index and at the delayed index, an
interior crossing position is never stationary; the delay lowers the
action when `q h E₀ < 0` and the advance lowers it when `0 < q h E₀`.  With
the committed sign of `electricFieldScaled` a positive `q h E₀` is the
field pointing from `seamLeft e` to `seamRight e` at positive `q h`, so the
forward crossing moves earlier when it runs with the field and later when
it runs against it. -/
theorem static_field_drift (q h τ : ℝ) (hq : q ≠ 0) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (n k i : ℕ)
    (hb : BlockPattern w n k (i + 1)) (hi : i + 1 < k) (hN : n + k < N + 1) (e : Fin 30)
    (hf : w.steps (n + (i + 1)) = .forward e) (E₀ : ℝ) (hE₀ : E₀ ≠ 0)
    (hs1 : electricFieldScaled h A φ (n + i + 1) e = E₀)
    (hs2 : electricFieldScaled h A φ (n + (i + 1) + 1) e = E₀) :
    ¬ PositionStationary q h τ N A φ ρ J w n k (i + 1) e ∧
    (q * (h * E₀) < 0 →
      transportedAction q h τ N A φ ρ J
          (exchangeVariation w (n + (i + 1)) e hf
            (pattern_rest_after w n k (i + 1) hb hi)).worldline <
        transportedAction q h τ N A φ ρ J w) ∧
    (0 < q * (h * E₀) →
      transportedAction q h τ N A φ ρ J
          (advanceVariation w (n + i) e (pattern_rest_before w n k i hb (by omega))
            hf).worldline <
        transportedAction q h τ N A φ ρ J w) := by
  refine ⟨fun hst ↦ ?_, fun hneg ↦ ?_, fun hpos ↦ ?_⟩
  · have := ((position_stationary_iff q h τ hq hh N A φ ρ J w n k (i + 1) hb (by omega) hN e
      hf).mp hst).1 hi
    rw [hs2] at this
    exact hE₀ this
  · rw [delay_lowers_iff q h τ hh N A φ ρ J w (n + (i + 1)) (by omega) e hf _, hs2]
    exact hneg
  · rw [advance_lowers_iff q h τ hh N A φ ρ J w (n + i) (by omega) e _ hf, hs1]
    exact hpos

end

/-! ## (3) Momentum timing and the field-side transfer -/

noncomputable section

/-- One step of the generated path: `(τ, stepVector)`. -/
theorem generatedPath_step (τ : ℝ) (w : SeamStepWorldline) (k : ℕ) :
    generatedPath τ w (k + 1) - generatedPath τ w k = (τ, stepVector (w.steps k)) := by
  refine Prod.ext ?_ ?_
  · rw [Prod.fst_sub, generatedPath_fst, generatedPath_fst]
    push_cast
    ring
  · rw [Prod.snd_sub, generatedPath_increment]

/-- **The generated path of a closed two-step variation** agrees with the
original at every index other than `j + 1`. -/
theorem generatedPath_closed_variation (τ : ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ClosedTwoStepVariation w j) (m : ℕ) (hm : m ≠ j + 1) :
    generatedPath τ v.worldline m = generatedPath τ w m := by
  refine Prod.ext rfl ?_
  rw [← ray_port_eq_spatial τ v.worldline m, ← ray_port_eq_spatial τ w m, v.varied_port_eq m hm]
  rfl

/-- The discrete momentum accumulated over the interval `[a, b]` of the
generated path: `2 m` times the increment. -/
def intervalMomentum (m τ : ℝ) (w : SeamStepWorldline) (a b : ℕ) : Herm2 :=
  (2 * m) • (generatedPath τ w b - generatedPath τ w a)

/-- The interval momentum is the sum of the per-step discrete momenta. -/
theorem intervalMomentum_succ (m τ : ℝ) (w : SeamStepWorldline) (a b : ℕ) :
    intervalMomentum m τ w a (b + 1) =
      intervalMomentum m τ w a b + discreteMomentum m (generatedPath τ w) b := by
  unfold intervalMomentum discreteMomentum
  rw [← smul_add]
  congr 1
  abel

/-- **Interval momentum under a closed variation.**  Over any interval whose
endpoints avoid `j + 1` the momentum is unchanged. -/
theorem intervalMomentum_closed_variation (m τ : ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ClosedTwoStepVariation w j) (a b : ℕ) (ha : a ≠ j + 1) (hb : b ≠ j + 1) :
    intervalMomentum m τ v.worldline a b = intervalMomentum m τ w a b := by
  unfold intervalMomentum
  rw [generatedPath_closed_variation τ w v a ha, generatedPath_closed_variation τ w v b hb]

/-- The block momentum: the interval momentum over the `k + 1` steps of the
block from `n`. -/
def blockMomentum (m τ : ℝ) (w : SeamStepWorldline) (n k : ℕ) : Herm2 :=
  intervalMomentum m τ w n (n + k + 1)

/-- The spatial displacement across a block pattern is the step vector of
its crossing. -/
theorem pattern_displacement (w : SeamStepWorldline) (n k i : ℕ) (hb : BlockPattern w n k i)
    (hik : i ≤ k) :
    spatialAt w (n + k + 1) - spatialAt w n = stepVector (w.steps (n + i)) := by
  have hblock : ∑ p ∈ Finset.range (k + 1), stepVector (w.steps (n + p)) =
      stepVector (w.steps (n + i)) := by
    rw [Finset.sum_eq_single i (fun p hp hpi ↦ ?_) (fun hi ↦ absurd (Finset.mem_range.mpr
      (by omega)) hi)]
    rw [hb.1 p (Finset.mem_range.mp hp) hpi, stepVector_rest]
  unfold spatialAt
  rw [sub_eq_iff_eq_add', ← Finset.sum_range_add_sum_Ico _ (show n ≤ n + k + 1 by omega),
    Finset.sum_Ico_eq_sum_range, show n + k + 1 - n = k + 1 by omega, hblock]

/-- **Block momentum.**  Across a block pattern the momentum is
`2 m ((k + 1) τ, s)` with `s` the step vector of the crossing: one seam
vector, independent of the position of the crossing. -/
theorem blockMomentum_eq (m τ : ℝ) (w : SeamStepWorldline) (n k i : ℕ)
    (hb : BlockPattern w n k i) (hik : i ≤ k) :
    blockMomentum m τ w n k = (2 * m) • (((k + 1 : ℝ) * τ, stepVector (w.steps (n + i))) : Herm2) := by
  unfold blockMomentum intervalMomentum
  congr 1
  refine Prod.ext ?_ ?_
  · rw [Prod.fst_sub, generatedPath_fst, generatedPath_fst]
    push_cast
    ring
  · rw [Prod.snd_sub, generatedPath_snd, generatedPath_snd, pattern_displacement w n k i hb hik]

/-- **Delay keeps the block momentum.** -/
theorem blockMomentum_delay (m τ : ℝ) (w : SeamStepWorldline) (n k i : ℕ)
    (hb : BlockPattern w n k i) (hi : i < k) (e : Fin 30) (hf : w.steps (n + i) = .forward e) :
    blockMomentum m τ (exchangeVariation w (n + i) e hf (pattern_rest_after w n k i hb hi)).worldline
        n k =
      blockMomentum m τ w n k :=
  intervalMomentum_closed_variation m τ w _ n (n + k + 1) (by omega) (by omega)

/-- **The exchange swaps the two per-step momenta.**  The crossing quantum
`2 m (τ, s_e)` moves from step `j` to step `j + 1` and the rest quantum
`2 m (τ, 0)` moves the other way. -/
theorem delay_momentum_swap (m τ : ℝ) (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hf : w.steps j = .forward e) (hr : w.steps (j + 1) = .rest) :
    discreteMomentum m (generatedPath τ (exchangeVariation w j e hf hr).worldline) (j + 1) =
        discreteMomentum m (generatedPath τ w) j ∧
      discreteMomentum m (generatedPath τ (exchangeVariation w j e hf hr).worldline) j =
        discreteMomentum m (generatedPath τ w) (j + 1) ∧
      discreteMomentum m (generatedPath τ w) j = (2 * m) • ((τ, seamVector e) : Herm2) ∧
      discreteMomentum m (generatedPath τ w) (j + 1) = (2 * m) • ((τ, 0) : Herm2) := by
  have hs' : (exchangeVariation w j e hf hr).s' = .rest := rfl
  have hs'' : (exchangeVariation w j e hf hr).s'' = .forward e := rfl
  unfold discreteMomentum
  rw [generatedPath_step, generatedPath_step, generatedPath_step, generatedPath_step,
    ClosedTwoStepVariation.worldline_steps, ClosedTwoStepVariation.worldline_steps,
    replaceTwo_at, replaceTwo_succ, hs', hs'', hf, hr, stepVector_forward, stepVector_rest]
  exact ⟨rfl, rfl, rfl, rfl⟩

/-! ### The field-side transfer at the two timings -/

/-- **Transfer at the original timing.**  Along the scaled Ampere evolution
sourced by the hopping current of `w`, the field energy moves by
`(q/2)(E_j e + E_{j+1} e)` across the crossing step and by zero across the
rest step. -/
theorem transfer_original (q h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hf : w.steps j = .forward e) (hr : w.steps (j + 1) = .rest)
    (hAmp : AmpereEvolutionScaled h A φ (hoppingCurrent q h (hoppingPath w))) :
    fieldEnergyScaled h A φ (j + 1) - fieldEnergyScaled h A φ j =
        (q / 2) * (electricFieldScaled h A φ j e + electricFieldScaled h A φ (j + 1) e) ∧
      fieldEnergyScaled h A φ (j + 2) - fieldEnergyScaled h A φ (j + 1) = 0 := by
  constructor
  · rw [hopping_work_energy_forward q h hh A φ (hoppingPath w) hAmp j e
      (by rw [hoppingPath_port]; exact port_of_forward w j e hf)
      (by rw [hoppingPath_port]; exact port_succ_of_forward w j e hf)]
    ring
  · rw [hopping_work_energy_rest q h hh A φ (hoppingPath w) hAmp (j + 1)
      (by rw [hoppingPath_port, hoppingPath_port]; exact port_succ_of_rest w (j + 1) hr), sub_self]

/-- **Transfer at the delayed timing.**  Along the scaled Ampere evolution
sourced by the hopping current of the delayed worldline, the field energy
moves by zero across the rest step and by `(q/2)(E'_{j+1} e + E'_{j+2} e)`
across the delayed crossing. -/
theorem transfer_delayed (q h : ℝ) (hh : h ≠ 0) (A' : ℕ → Fin 30 → ℝ) (φ' : ℕ → Fin 12 → ℝ)
    (w : SeamStepWorldline) (j : ℕ) (e : Fin 30)
    (hf : w.steps j = .forward e) (hr : w.steps (j + 1) = .rest)
    (hAmp' : AmpereEvolutionScaled h A' φ'
      (hoppingCurrent q h (hoppingPath (exchangeVariation w j e hf hr).worldline))) :
    fieldEnergyScaled h A' φ' (j + 1) - fieldEnergyScaled h A' φ' j = 0 ∧
      fieldEnergyScaled h A' φ' (j + 2) - fieldEnergyScaled h A' φ' (j + 1) =
        (q / 2) * (electricFieldScaled h A' φ' (j + 1) e + electricFieldScaled h A' φ' (j + 2) e) := by
  have hpj : (exchangeVariation w j e hf hr).worldline.port j = seamLeft e := by
    rw [ClosedTwoStepVariation.varied_port_eq _ j (by omega)]
    exact port_of_forward w j e hf
  have hmid : (exchangeVariation w j e hf hr).worldline.port (j + 1) = seamLeft e := by
    rw [ClosedTwoStepVariation.varied_port_mid]
    exact port_of_forward w j e hf
  have hend : (exchangeVariation w j e hf hr).worldline.port (j + 2) = seamRight e := by
    rw [ClosedTwoStepVariation.varied_port_eq _ (j + 2) (by omega),
      port_succ_of_rest w (j + 1) hr, port_succ_of_forward w j e hf]
  constructor
  · rw [hopping_work_energy_rest q h hh A' φ' _ hAmp' j
      (by rw [hoppingPath_port, hoppingPath_port, hmid, hpj]), sub_self]
  · rw [hopping_work_energy_forward q h hh A' φ' _ hAmp' (j + 1) e
      (by rw [hoppingPath_port]; exact hmid) (by rw [hoppingPath_port]; exact hend)]
    ring

/-- **Action difference and field transfer.**  With the seam field on `e`
static across the crossing step of `w` (`E_j e = E_{j+1} e`) and the field
along the Ampere evolution of `w`, the transfer across the crossing step is
`q E_{j+1} e`, and the action difference of the delay is `h` times that
transfer; with a second field along the Ampere evolution of the delayed
worldline, whose transfer across step `j` is zero, the action difference is
`-h` times the difference of the two step-`j` transfers.  The factor is
exactly `h`: the action pairs the load and current with weight `h`, the
energy balance pairs the current `-(q/h)` against `(h/2)(E_j + E_{j+1})`. -/
theorem exchange_action_eq_h_mul_transfer (q h τ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (A' : ℕ → Fin 30 → ℝ) (φ' : ℕ → Fin 12 → ℝ)
    (w : SeamStepWorldline) (j : ℕ) (hj : j + 1 < N + 1) (e : Fin 30)
    (hf : w.steps j = .forward e) (hr : w.steps (j + 1) = .rest)
    (hAmp : AmpereEvolutionScaled h A φ (hoppingCurrent q h (hoppingPath w)))
    (hAmp' : AmpereEvolutionScaled h A' φ'
      (hoppingCurrent q h (hoppingPath (exchangeVariation w j e hf hr).worldline)))
    (hstatic : electricFieldScaled h A φ j e = electricFieldScaled h A φ (j + 1) e) :
    fieldEnergyScaled h A φ (j + 1) - fieldEnergyScaled h A φ j =
        q * electricFieldScaled h A φ (j + 1) e ∧
      transportedAction q h τ N A φ ρ J (exchangeVariation w j e hf hr).worldline -
          transportedAction q h τ N A φ ρ J w =
        h * (fieldEnergyScaled h A φ (j + 1) - fieldEnergyScaled h A φ j) ∧
      transportedAction q h τ N A φ ρ J (exchangeVariation w j e hf hr).worldline -
          transportedAction q h τ N A φ ρ J w =
        -(h * ((fieldEnergyScaled h A' φ' (j + 1) - fieldEnergyScaled h A' φ' j) -
          (fieldEnergyScaled h A φ (j + 1) - fieldEnergyScaled h A φ j))) := by
  have ht := (transfer_original q h hh A φ w j e hf hr hAmp).1
  rw [hstatic] at ht
  have ht' : fieldEnergyScaled h A φ (j + 1) - fieldEnergyScaled h A φ j =
      q * electricFieldScaled h A φ (j + 1) e := by rw [ht]; ring
  have hd := (transfer_delayed q h hh A' φ' w j e hf hr hAmp').1
  have hact := (exchange_action_difference q h τ hh N A φ ρ J w j hj e hf hr).1
  refine ⟨ht', ?_, ?_⟩
  · rw [hact, ht']
    ring
  · rw [hact, hd, ht']
    ring

end

/-! ## (4) Non-forcing -/

noncomputable section

/-- **Two charges, two in-block balances.**  Distinct charges give distinct
in-block delay differences whenever the seam field at the delayed index is
nonzero. -/
theorem two_charges_two_block_balances (q₁ q₂ h τ : ℝ) (hq : q₁ ≠ q₂) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (n k i : ℕ) (hb : BlockPattern w n k i)
    (hi : i < k) (hN : n + k < N + 1) (e : Fin 30) (hf : w.steps (n + i) = .forward e)
    (hE : electricFieldScaled h A φ (n + i + 1) e ≠ 0) :
    transportedAction q₁ h τ N A φ ρ J
        (exchangeVariation w (n + i) e hf (pattern_rest_after w n k i hb hi)).worldline -
        transportedAction q₁ h τ N A φ ρ J w ≠
      transportedAction q₂ h τ N A φ ρ J
        (exchangeVariation w (n + i) e hf (pattern_rest_after w n k i hb hi)).worldline -
        transportedAction q₂ h τ N A φ ρ J w :=
  two_charges_two_balances q₁ q₂ h τ hq hh N A φ ρ J w (n + i) (by omega) e hf _ hE

/-- **The unit drops out of every in-block difference.** -/
theorem block_difference_forgets_unit (q h τ₁ τ₂ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) (n k i : ℕ) (hb : BlockPattern w n k i)
    (hi : i < k) (hN : n + k < N + 1) (e : Fin 30) (hf : w.steps (n + i) = .forward e) :
    transportedAction q h τ₁ N A φ ρ J
        (exchangeVariation w (n + i) e hf (pattern_rest_after w n k i hb hi)).worldline -
        transportedAction q h τ₁ N A φ ρ J w =
      transportedAction q h τ₂ N A φ ρ J
        (exchangeVariation w (n + i) e hf (pattern_rest_after w n k i hb hi)).worldline -
        transportedAction q h τ₂ N A φ ρ J w :=
  difference_forgets_unit q h τ₁ τ₂ hh N A φ ρ J w _ (by omega)

/-- **The unit drops out of every shifted one-step difference.** -/
theorem shifted_difference_forgets_unit (q h τ₁ τ₂ : ℝ) (hh : h ≠ 0) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (w : SeamStepWorldline) {j : ℕ}
    (v : ShiftedOneStepVariation w j) (hj : j < N + 1) :
    transportedAction q h τ₁ N A φ ρ J v.worldline - transportedAction q h τ₁ N A φ ρ J w =
      transportedAction q h τ₂ N A φ ρ J v.worldline -
        transportedAction q h τ₂ N A φ ρ J w := by
  rw [transportedAction_shifted_variation q h τ₁ hh N A φ ρ J w v hj,
    transportedAction_shifted_variation q h τ₂ hh N A φ ρ J w v hj]

end

end OPH.TimelikeClassForceLaw

/- Axiom audit: standard axioms only (`propext`, `Classical.choice`,
`Quot.sound`); no `sorry`, no `native_decide`, no project axiom. -/

#print axioms OPH.TimelikeClassForceLaw.ShiftedOneStepVariation.shifted_port_le
#print axioms OPH.TimelikeClassForceLaw.ShiftedOneStepVariation.shifted_port_succ
#print axioms OPH.TimelikeClassForceLaw.transportedAction_shifted_variation
#print axioms OPH.TimelikeClassForceLaw.shifted_variation_invariant
#print axioms OPH.TimelikeClassForceLaw.crossing_to_rest_invariant
#print axioms OPH.TimelikeClassForceLaw.shifted_port_rejoin
#print axioms OPH.TimelikeClassForceLaw.shifted_variation_rejoin
#print axioms OPH.TimelikeClassForceLaw.restBlock_iff_pattern
#print axioms OPH.TimelikeClassForceLaw.twoRestWorldline_restDiluted
#print axioms OPH.TimelikeClassForceLaw.advance_action_difference
#print axioms OPH.TimelikeClassForceLaw.advance_backward_action_difference
#print axioms OPH.TimelikeClassForceLaw.delay_pattern
#print axioms OPH.TimelikeClassForceLaw.delay_steps_outside
#print axioms OPH.TimelikeClassForceLaw.delay_action_difference
#print axioms OPH.TimelikeClassForceLaw.advance_block_action_difference
#print axioms OPH.TimelikeClassForceLaw.delay_lowers_iff
#print axioms OPH.TimelikeClassForceLaw.advance_lowers_iff
#print axioms OPH.TimelikeClassForceLaw.advance_stationary_iff
#print axioms OPH.TimelikeClassForceLaw.position_stationary_iff
#print axioms OPH.TimelikeClassForceLaw.static_field_drift
#print axioms OPH.TimelikeClassForceLaw.generatedPath_closed_variation
#print axioms OPH.TimelikeClassForceLaw.intervalMomentum_closed_variation
#print axioms OPH.TimelikeClassForceLaw.blockMomentum_eq
#print axioms OPH.TimelikeClassForceLaw.blockMomentum_delay
#print axioms OPH.TimelikeClassForceLaw.delay_momentum_swap
#print axioms OPH.TimelikeClassForceLaw.transfer_original
#print axioms OPH.TimelikeClassForceLaw.transfer_delayed
#print axioms OPH.TimelikeClassForceLaw.exchange_action_eq_h_mul_transfer
#print axioms OPH.TimelikeClassForceLaw.two_charges_two_block_balances
#print axioms OPH.TimelikeClassForceLaw.block_difference_forgets_unit
#print axioms OPH.TimelikeClassForceLaw.shifted_difference_forgets_unit
