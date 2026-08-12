import FiniteConditionalRepair
import FirstLawIdentity
import FluctuationTheorems
import StationaryRealization
import CapFirstLaw

set_option autoImplicit false

namespace OPH.Thermodynamics

/-!
# The composed four-law adequacy surface

One parameterized package for the observation-ledger rows OL-E1 (the
four laws), OL-E2 (the exact fluctuation identities), OL-E3 (the
Landauer bound), and OL-E4 (the arrow of time from repair), V3 issue
#732. Every theorem below is a committed result of the thermodynamics
modules, re-exported through one named premise bundle so that a ledger
row cites a single surface instead of five scattered conditionals.

`RepairLawData` is the premise-register row PR-07, the declared repair
law: one faithful reference law shared by the state and transition
sides, together with the complete repaired-visible datum whose fibres
the conditional-resampling kernel resamples. The B12 and B20 bounded
no-gos (issues #688 and #725; module `CommonReferenceObstruction` with
the compatibility witness `CommonObjectBinding`, and the committed
random-scan preflight `code/b20_random_scan/`) record that the realized
source artifacts do not supply this common reference, and that record
is the justification for carrying the row as declared, disposition
axiomatize. No
theorem in this module asserts a source selection of the reference, of
the visible datum, or of the kernel.

`CalibrationData` is the clock-and-energy calibration register row,
shared with the constants lane: the identification
`-log ref = beta * energy + logZ` of the modular weight with a physical
energy at one positive multiplier. The zeroth-law, third-law, and
Landauer statements consume it; the second-law statements do not.

Composition boundaries, stated exactly:

* the exact first-law split `first_exact_split` is an algebraic matrix
  identity with no repair-law hypothesis, so it does not factor through
  `RepairLawData`; it is re-exported side by side inside the namespace;
* the stationary-kernel second law `second_contraction_of_stationary`
  quantifies over arbitrary stochastic kernels preserving the declared
  reference, covering nonreversible repair chains; detailed balance
  enters only the fluctuation relations;
* the refinement-uniform low-temperature control is the register row
  PR-08; it is consumed by the refinement-family strengthening in
  `LowTemperatureControl` and is a hypothesis of no statement here;
* every remaining theorem takes `RepairLawData`, and where stated
  `CalibrationData`, as its only nonstructural input.

The laws are exact finite theorems under the register. Laboratory
thermodynamics requires the calibration row. Nothing here claims the
axioms force the repair law; that is the register's open disposition.
-/

/-- **Register row PR-07: the declared repair law.** One faithful
reference law on the finite state space, shared by the state and
transition sides, together with the complete repaired-visible datum
`visible` whose fibres the repair kernel resamples. The B12 and B20
bounded no-gos are the record of why this row is declared rather than
derived: the realized source artifacts supply two different finite
objects where one common reference is required. The finite theorems
below consume `ref_pos` and the fibre structure; `ref_law` records that
the declared reference is a probability law and is consumed by
`ref_eq_gibbs`. -/
structure RepairLawData (Ω : Type*) [Fintype Ω] [DecidableEq Ω]
    (B : Type*) [DecidableEq B] where
  /-- The common faithful reference law. -/
  ref : Ω → ℝ
  /-- Faithfulness: every microstate carries positive reference mass. -/
  ref_pos : ∀ x, 0 < ref x
  /-- The reference is a probability law. -/
  ref_law : ∑ x, ref x = 1
  /-- The complete repaired visible datum. -/
  visible : Ω → B

/-- **The clock-and-energy calibration register row**, relative to a
declared repair law: the modular weight of the reference is identified
with a physical energy at one positive multiplier,
`-log ref = beta * energy + logZ`. Shared with the constants lane. The
identification itself is a declared premise; nothing in this module
derives it. -/
structure CalibrationData {Ω : Type*} [Fintype Ω] [DecidableEq Ω]
    {B : Type*} [DecidableEq B] (D : RepairLawData Ω B) where
  /-- The calibrated physical energy observable. -/
  energy : Ω → ℝ
  /-- The inverse-temperature multiplier. -/
  beta : ℝ
  /-- The log-normalization constant of the identification. -/
  logZ : ℝ
  /-- The multiplier is positive. -/
  beta_pos : 0 < beta
  /-- The thermal identification of the modular weight. -/
  thermal : ∀ x, -Real.log (D.ref x) = beta * energy x + logZ

namespace FourLawSurface

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]
variable {B : Type*} [DecidableEq B]
variable (D : RepairLawData Ω B)

/-- The conditional-resampling repair kernel of the declared repair
law: resample the unresolved variables from the reference law inside
the fibre of the complete repaired visible datum. -/
noncomputable def repairKernel : Ω → Ω → ℝ :=
  heatBath D.ref D.visible

theorem repairKernel_def : repairKernel D = heatBath D.ref D.visible :=
  rfl

/-! ## Kernel receipts

The declared repair law determines one kernel, and that kernel is
stochastic, stationary for the reference, reversible, idempotent, and
the information projection onto each visible fibre. These receipts are
the common hypotheses of the four laws below. -/

theorem repairKernel_nonneg (x y : Ω) : 0 ≤ repairKernel D x y :=
  heatBath_nonneg D.ref_pos x y

theorem repairKernel_row_sum (x : Ω) : ∑ y, repairKernel D x y = 1 :=
  heatBath_row_sum D.ref_pos x

theorem repairKernel_stationary (y : Ω) :
    push D.ref (repairKernel D) y = D.ref y :=
  heatBath_stationary D.ref_pos y

theorem repairKernel_detailedBalance (x y : Ω) :
    D.ref x * repairKernel D x y = D.ref y * repairKernel D y x :=
  heatBath_detailedBalance x y

theorem repairKernel_idempotent (x y : Ω) :
    ∑ z, repairKernel D x z * repairKernel D z y = repairKernel D x y :=
  heatBath_idempotent D.ref_pos x y

/-- The kernel row at `x` is the information projection onto the fibre
of the visible datum at `x`: among fibre-supported laws it minimizes
relative entropy to the declared reference. -/
theorem repairKernel_row_optimal (x : Ω) (r : Ω → ℝ)
    (hr0 : ∀ y, 0 ≤ r y) (hr1 : ∑ y, r y = 1)
    (hsupp : ∀ y, D.visible y ≠ D.visible x → r y = 0) :
    kl (repairKernel D x) D.ref ≤ kl r D.ref :=
  heatBath_row_optimal D.ref_pos x r hr0 hr1 hsupp

/-! ## Zeroth law -/

/-- Under the calibration row, the declared reference is the Gibbs
state of the calibrated energy at the calibrated multiplier. This is
the point where `ref_law` is consumed. -/
theorem ref_eq_gibbs [Nonempty Ω] (C : CalibrationData D) (x : Ω) :
    D.ref x = gibbs C.energy C.beta x := by
  have hpt : ∀ y : Ω, D.ref y
      = gibbsWeight C.energy C.beta y * Real.exp (-C.logZ) := by
    intro y
    have hlog : Real.log (D.ref y)
        = -C.beta * C.energy y + -C.logZ := by
      have h := C.thermal y
      linarith
    calc D.ref y
        = Real.exp (Real.log (D.ref y)) :=
          (Real.exp_log (D.ref_pos y)).symm
      _ = Real.exp (-C.beta * C.energy y + -C.logZ) := by rw [hlog]
      _ = gibbsWeight C.energy C.beta y * Real.exp (-C.logZ) := by
          rw [Real.exp_add]
          rfl
  have hZ : partitionZ C.energy C.beta * Real.exp (-C.logZ) = 1 := by
    calc partitionZ C.energy C.beta * Real.exp (-C.logZ)
        = ∑ y, gibbsWeight C.energy C.beta y * Real.exp (-C.logZ) := by
          unfold partitionZ
          rw [Finset.sum_mul]
      _ = ∑ y, D.ref y :=
          Finset.sum_congr rfl fun y _ => (hpt y).symm
      _ = 1 := D.ref_law
  have hZpos := partitionZ_pos C.energy C.beta
  have hexp : Real.exp (-C.logZ) = 1 / partitionZ C.energy C.beta := by
    rw [eq_div_iff (ne_of_gt hZpos)]
    linear_combination hZ
  rw [hpt x, hexp]
  unfold gibbs
  ring

/-- **Zeroth law, thermometer form.** A nondegenerate thermometer
identifies the multiplier: two Gibbs states of the calibrated energy
that agree as distributions have equal `beta`. Contact of `A` with `B`
and of `B` with `C` therefore forces one common temperature through any
thermometer with two distinct calibrated energy levels. -/
theorem zeroth_thermometer [Nonempty Ω] (C : CalibrationData D)
    (beta1 beta2 : ℝ) (i j : Ω) (hij : C.energy i ≠ C.energy j)
    (h : gibbs C.energy beta1 = gibbs C.energy beta2) :
    beta1 = beta2 :=
  gibbs_beta_injective C.energy beta1 beta2 i j hij h

/-- **Zeroth law, multiplier equality.** Two calibration rows over the
same declared reference and the same nondegenerate energy carry the
same multiplier: the declared repair law identifies one inverse
temperature. -/
theorem zeroth_multiplier_unique [Nonempty Ω]
    (C1 C2 : CalibrationData D) (henergy : C1.energy = C2.energy)
    (i j : Ω) (hij : C1.energy i ≠ C1.energy j) :
    C1.beta = C2.beta := by
  have h : gibbs C1.energy C1.beta = gibbs C1.energy C2.beta := by
    funext x
    calc gibbs C1.energy C1.beta x
        = D.ref x := (ref_eq_gibbs D C1 x).symm
      _ = gibbs C2.energy C2.beta x := ref_eq_gibbs D C2 x
      _ = gibbs C1.energy C2.beta x := by rw [henergy]
  exact gibbs_beta_injective C1.energy C1.beta C2.beta i j hij h

/-! ## First law -/

/-- **First law, exact bookkeeping split.** The internal-energy change
across a joint update of state and energy observable splits exactly
into the heat increment, the work increment, and one explicit bilinear
cross term. This is an algebraic matrix identity with no repair-law
hypothesis; it is re-exported side by side because it does not factor
through `RepairLawData`. Physical content additionally requires the
calibration row and a protocol distinguishing controlled changes of the
observable from state changes. -/
theorem first_exact_split {n : ℕ}
    (ρ dρ H dH : Matrix (Fin n) (Fin n) ℂ) :
    internalEnergy (ρ + dρ) (H + dH) - internalEnergy ρ H
      = heatIncrement dρ H + workIncrement ρ dH
        + (Matrix.trace (dρ * dH)).re :=
  firstLaw_split ρ dρ H dH

/-- **First law, modular form with exact remainder.** The entropy
difference to the declared reference equals the modular-energy
difference minus the relative entropy; the first-order first law
`δS = δ⟨K⟩` is this identity with the quadratic-order deficit
dropped. -/
theorem first_cap_exact (p : Ω → ℝ) :
    shannon p - shannon D.ref
      = (∑ x, p x * (-Real.log (D.ref x)))
        - (∑ x, D.ref x * (-Real.log (D.ref x))) - kl p D.ref :=
  cap_firstLaw_exact p D.ref D.ref_pos

/-- **First law, conservation under pure repair.** One repair step
preserves the mean of every observable measurable through the visible
datum. With the calibrated energy fibre-measurable, pure repair
exchanges no energy: the split of `first_exact_split` assigns the whole
change of a repair step to the heat channel at fixed observable, and
that channel vanishes on fibre-measurable energies. -/
theorem first_repair_conserves_fibre_mean (Q : Ω → ℝ)
    (hQ : ∀ x y, D.visible x = D.visible y → Q x = Q y) (p : Ω → ℝ) :
    ∑ y, push p (repairKernel D) y * Q y = ∑ x, p x * Q x :=
  push_heatBath_fixes_mean D.ref_pos Q hQ p

/-! ## Second law -/

/-- **Second law, contraction form.** Relative entropy to the declared
reference is nonincreasing under one repair step. -/
theorem second_contraction (p : Ω → ℝ) (hp : ∀ x, 0 ≤ p x) :
    kl (push p (repairKernel D)) D.ref ≤ kl p D.ref :=
  heatBath_secondLaw D.ref_pos p hp

/-- **Second law without detailed balance.** Any stochastic kernel
preserving the declared reference contracts relative entropy to it.
This covers nonreversible repair chains; the reversible
conditional-resampling kernel is one instance. -/
theorem second_contraction_of_stationary (p : Ω → ℝ) (K : Ω → Ω → ℝ)
    (hp : ∀ x, 0 ≤ p x) (hK0 : ∀ x y, 0 ≤ K x y)
    (hK1 : ∀ x, ∑ y, K x y = 1) (hstat : push D.ref K = D.ref) :
    kl (push p K) D.ref ≤ kl p D.ref :=
  stationary_secondLaw D.ref p K hp D.ref_pos hK0 hK1 hstat

/-- **Second law, Clausius form.** Under one repair step the entropy
change dominates the modular-energy change. On the calibrated thermal
branch this is `ΔS ≥ β Q`. -/
theorem second_clausius (p : Ω → ℝ) (hp : ∀ x, 0 ≤ p x) :
    (∑ x, push p (repairKernel D) x * (-Real.log (D.ref x)))
      - (∑ x, p x * (-Real.log (D.ref x)))
      ≤ shannon (push p (repairKernel D)) - shannon p :=
  clausius p D.ref (repairKernel D) hp D.ref_pos
    (repairKernel_nonneg D) (repairKernel_row_sum D)
    (repairKernel_stationary D)

/-- **Integral fluctuation identity.** The exponential of minus the
fluctuating entropy production of one repair step averages to one
exactly; the second law is the Jensen consequence. -/
theorem second_integral_fluctuation (p : Ω → ℝ)
    (hp : ∀ x, 0 < p x) (hp1 : ∑ x, p x = 1) :
    ∑ x, ∑ y, p x * repairKernel D x y
      * Real.exp
        (-(sigmaEP D.ref p (push p (repairKernel D)) x y)) = 1 :=
  heatBath_integral_fluctuation D.ref_pos p hp hp1

/-- **Detailed fluctuation relation (Crooks form).** The forward step
weight equals the reversed step weight amplified by the exponential of
the entropy production, at every pair of states. -/
theorem second_crooks (p : Ω → ℝ) (hp : ∀ x, 0 < p x) (x y : Ω) :
    p x * repairKernel D x y
      = Real.exp (sigmaEP D.ref p (push p (repairKernel D)) x y)
        * (push p (repairKernel D) y * repairKernel D y x) :=
  heatBath_crooks D.ref_pos p hp x y

/-! ## Arrow of time -/

/-- **Arrow of time from repair.** The mean fluctuating entropy
production of one repair step is nonnegative: it equals the certified
relative-entropy descent of `second_contraction` exactly, so every
repair step points the same way. The emergent rung of this row opens
only if the repair law is promoted into the architecture; this theorem
is the structural rung. -/
theorem arrow_mean_entropy_production_nonneg (p : Ω → ℝ)
    (hp : ∀ x, 0 < p x) :
    0 ≤ ∑ x, ∑ y, p x * repairKernel D x y
      * sigmaEP D.ref p (push p (repairKernel D)) x y := by
  have hmean := sigma_mean_eq_kl_descent D.ref p (repairKernel D) hp
    (repairKernel_row_sum D)
    (heatBath_preserves_pos D.ref_pos p hp)
  have hdesc := second_contraction D p fun x => le_of_lt (hp x)
  rw [hmean]
  linarith

/-! ## Third law -/

/-- **Third law, excited-mass bound.** With ground energy `E0`, gap
`Δ > 0`, and the calibrated energy, the excited Gibbs mass is bounded
by `(d - g₀)/g₀ · exp (-β Δ)`. -/
theorem third_excited_mass_bound [Nonempty Ω] (C : CalibrationData D)
    (beta E0 Δ : ℝ) (hβ : 0 ≤ beta) (hΔ : 0 < Δ)
    (hground : ∀ x, C.energy x = E0 ∨ E0 + Δ ≤ C.energy x)
    (hne : (Finset.univ.filter (fun x => C.energy x = E0)).Nonempty) :
    excitedMass C.energy beta E0
      ≤ ((Finset.univ.filter (fun x => C.energy x ≠ E0)).card : ℝ)
        / ((Finset.univ.filter (fun x => C.energy x = E0)).card : ℝ)
        * Real.exp (-beta * Δ) :=
  excitedMass_le C.energy beta E0 Δ hβ hΔ hground hne

/-- **Third law, quantitative threshold.** Beyond the explicit
threshold `β₀ = log ((d - g₀)/(g₀ ε)) / Δ` the excited mass of the
calibrated Gibbs family is below `ε`. -/
theorem third_excited_mass_threshold [Nonempty Ω]
    (C : CalibrationData D) (beta E0 Δ eps : ℝ)
    (hβ : 0 ≤ beta) (hΔ : 0 < Δ) (heps : 0 < eps)
    (hground : ∀ x, C.energy x = E0 ∨ E0 + Δ ≤ C.energy x)
    (hne : (Finset.univ.filter (fun x => C.energy x = E0)).Nonempty)
    (hthreshold :
      Real.log
        (((Finset.univ.filter (fun x => C.energy x ≠ E0)).card : ℝ)
          / ((Finset.univ.filter
              (fun x => C.energy x = E0)).card : ℝ) / eps)
        < beta * Δ) :
    excitedMass C.energy beta E0 < eps :=
  excitedMass_lt_of_beta_large C.energy beta E0 Δ eps hβ hΔ heps
    hground hne hthreshold

/-- **Third law, finite-step unattainability.** One repair step
extinguishes no atom: a fully supported state stays fully supported,
so finitely many repair steps reach no rank-deficient ground-sector
state. The refinement-uniform strengthening consumes the register row
PR-08 and lives in `LowTemperatureControl`. -/
theorem third_no_step_extinguishes (p : Ω → ℝ) (hp : ∀ x, 0 < p x)
    (y : Ω) : 0 < push p (repairKernel D) y :=
  heatBath_preserves_pos D.ref_pos p hp y

/-! ## Landauer bound -/

/-- **Landauer bound.** Under the calibration row, one repair step that
lowers the entropy by at least `c` expels at least `c / β` of the
calibrated energy: erasing one bit costs at least `log 2 / β`, which is
`k_B T log 2` in physical units. -/
theorem landauer_bound (C : CalibrationData D) (p : Ω → ℝ) (c : ℝ)
    (hp : ∀ x, 0 ≤ p x)
    (herase : shannon (push p (repairKernel D)) - shannon p ≤ -c) :
    (∑ x, p x * C.energy x)
      - (∑ x, push p (repairKernel D) x * C.energy x)
      ≥ c / C.beta :=
  landauer p D.ref (repairKernel D) C.energy C.beta C.logZ c hp
    D.ref_pos (repairKernel_nonneg D) (repairKernel_row_sum D)
    (repairKernel_stationary D) C.beta_pos C.thermal herase

end FourLawSurface

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.FourLawSurface.repairKernel_nonneg
#print axioms OPH.Thermodynamics.FourLawSurface.repairKernel_row_sum
#print axioms OPH.Thermodynamics.FourLawSurface.repairKernel_stationary
#print axioms OPH.Thermodynamics.FourLawSurface.repairKernel_detailedBalance
#print axioms OPH.Thermodynamics.FourLawSurface.repairKernel_idempotent
#print axioms OPH.Thermodynamics.FourLawSurface.repairKernel_row_optimal
#print axioms OPH.Thermodynamics.FourLawSurface.ref_eq_gibbs
#print axioms OPH.Thermodynamics.FourLawSurface.zeroth_thermometer
#print axioms OPH.Thermodynamics.FourLawSurface.zeroth_multiplier_unique
#print axioms OPH.Thermodynamics.FourLawSurface.first_exact_split
#print axioms OPH.Thermodynamics.FourLawSurface.first_cap_exact
#print axioms OPH.Thermodynamics.FourLawSurface.first_repair_conserves_fibre_mean
#print axioms OPH.Thermodynamics.FourLawSurface.second_contraction
#print axioms OPH.Thermodynamics.FourLawSurface.second_contraction_of_stationary
#print axioms OPH.Thermodynamics.FourLawSurface.second_clausius
#print axioms OPH.Thermodynamics.FourLawSurface.second_integral_fluctuation
#print axioms OPH.Thermodynamics.FourLawSurface.second_crooks
#print axioms OPH.Thermodynamics.FourLawSurface.arrow_mean_entropy_production_nonneg
#print axioms OPH.Thermodynamics.FourLawSurface.third_excited_mass_bound
#print axioms OPH.Thermodynamics.FourLawSurface.third_excited_mass_threshold
#print axioms OPH.Thermodynamics.FourLawSurface.third_no_step_extinguishes
#print axioms OPH.Thermodynamics.FourLawSurface.landauer_bound
