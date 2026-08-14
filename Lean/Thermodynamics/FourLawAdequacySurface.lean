import FiniteConditionalRepair
import FirstLawIdentity
import FluctuationTheorems
import StationaryRealization
import CapFirstLaw
import LowTemperatureControl

set_option autoImplicit false

namespace OPH.Thermodynamics

/-!
# The composed four-law adequacy surface

One parameterized package for the observation-ledger rows OL-E1 (the
four laws), OL-E2 (the exact fluctuation identities), OL-E3 (the
Landauer bound), and OL-E4 (the arrow of time from repair), V3 issue
#732. Every theorem in `FourLawSurface` is a committed result of the
thermodynamics modules, stated through the named per-row premise
structures. The composed surface is the one theorem
`fourLaws_composed`: from the single typed antecedent bundle
`FourLawAntecedent`, whose fields are exactly the premise-register
rows PR-07, PR-15, and PR-08, the full conclusion record
`FourLawConclusions` holds for the one kernel `repairKernel A.repair`
and the one calibrated energy `A.calib.energy`. A ledger row cites
that single composed statement.

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
  identity with no repair-law hypothesis; its diagonal (classical)
  instantiation at the repair step and the calibrated energy is the
  threaded clause `first_heat_channel`, and the general matrix
  statement stays a standalone module theorem with no repair-law
  content;
* the stationary-kernel second law `second_contraction_of_stationary`
  quantifies over arbitrary stochastic kernels preserving the declared
  reference, covering nonreversible repair chains; detailed balance
  enters only the fluctuation relations;
* the refinement-uniform low-temperature control is the register row
  PR-08; the bundle field `RefinementAttachment` attaches it to the
  calibrated energy through a carrier identification, and the clauses
  `third_refinement_uniform` and `third_calibrated_member` of the
  composed conclusion consume it; no other clause needs it, and the
  per-clause doc comments name the rows each clause uses;
* the arrow clauses state nonnegativity for every strictly positive
  state, the exact dissipation identity, and strictness: a
  state-changing repair step produces strictly positive mean entropy
  production, and zero production forces the fibre-conditional
  reference form; states with zero atoms carry the nonnegativity
  clause only;
* every theorem of `FourLawSurface` takes `RepairLawData`, and where
  stated `CalibrationData`, as its only nonstructural input; the
  per-row structures stay exported for finer consumers.

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
hypothesis; its diagonal instantiation at the repair step and the
calibrated energy is the threaded theorem `first_heat_channel`, and
the general matrix statement carries no repair-law content. Physical
content additionally requires the calibration row and a protocol
distinguishing controlled changes of the observable from state
changes. -/
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
PR-08; its base theorems live in `LowTemperatureControl`, and the
clauses `third_refinement_uniform` and `third_calibrated_member` of
the composed conclusion thread it to the calibrated energy. -/
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

/-! ## Arrow of time, exact form

The nonnegativity clause above is the Jensen shadow of an exact
identity: the mean fluctuating entropy production of one repair step
equals the relative entropy from the state to its own repaired image.
The identity is classical information theory for conditional-
expectation projections; its role here is that the arrow clause of the
composed surface becomes a forcing statement rather than a bare
inequality. All three theorems consume register row PR-07 only. -/

/-- **Arrow of time, exact dissipation identity.** The mean fluctuating
entropy production of one repair step equals `D(p ‖ push p K)`, the
relative entropy from the state to its repaired image. Stated for
strictly positive states; states with zero atoms carry the
nonnegativity clause only. Consumes register row PR-07. -/
theorem arrow_mean_entropy_production_eq_kl_to_repaired (p : Ω → ℝ)
    (hp : ∀ x, 0 < p x) :
    ∑ x, ∑ y, p x * repairKernel D x y
        * sigmaEP D.ref p (push p (repairKernel D)) x y
      = kl p (push p (repairKernel D)) := by
  have hmean := sigma_mean_eq_kl_descent D.ref p (repairKernel D) hp
    (repairKernel_row_sum D)
    (heatBath_preserves_pos D.ref_pos p hp)
  have hpyth := heatBath_kl_pythagorean (b := D.visible) D.ref_pos p hp
  rw [hmean, repairKernel_def]
  linarith [hpyth]

/-- **Arrow of time, strictness.** A repair step that changes a
strictly positive state produces strictly positive mean entropy
production, and zero mean entropy production forces the state onto the
repaired manifold. Together with the exact identity this excludes a
degenerate always-zero arrow off the fixed states at the structural
rung; the emergent rung is gated on the simulator-export decision
recorded in the register row PR-07 notes. Consumes register row
PR-07. -/
theorem arrow_strict (p : Ω → ℝ) (hp : ∀ x, 0 < p x) :
    push p (repairKernel D) ≠ p ↔
      0 < ∑ x, ∑ y, p x * repairKernel D x y
        * sigmaEP D.ref p (push p (repairKernel D)) x y := by
  rw [arrow_mean_entropy_production_eq_kl_to_repaired D p hp]
  have hq : ∀ y, 0 < push p (repairKernel D) y :=
    heatBath_preserves_pos D.ref_pos p hp
  have hsum : ∑ x, p x = ∑ y, push p (repairKernel D) y :=
    (push_total p (repairKernel D) (repairKernel_row_sum D)).symm
  constructor
  · intro hne
    have hzero : kl p (push p (repairKernel D)) ≠ 0 := fun h =>
      hne ((kl_eq_zero_iff_of_pos p (push p (repairKernel D)) hp hq
        hsum).mp h).symm
    have hnn : 0 ≤ kl p (push p (repairKernel D)) :=
      kl_nonneg p (push p (repairKernel D)) (fun x => (hp x).le)
        (fun y => (hq y).le)
        (fun y hy => absurd hy (ne_of_gt (hq y))) hsum
    exact lt_of_le_of_ne hnn (Ne.symm hzero)
  · intro hpos heq
    have hself : kl p (push p (repairKernel D)) = 0 := by
      rw [heq]
      exact (kl_eq_zero_iff_of_pos p p hp hp rfl).mpr rfl
    rw [hself] at hpos
    exact lt_irrefl 0 hpos

/-- **Zero-dissipation characterization.** A state is fixed by one
repair step exactly when it carries the reference's conditional
structure on every visible fibre: within each fibre the state is the
reference reweighted by the ratio of state fibre mass to reference
fibre mass. Consumes register row PR-07. -/
theorem repair_fixed_iff (p : Ω → ℝ) :
    push p (repairKernel D) = p ↔
      ∀ y, p y = D.ref y * fiberMass p D.visible y
        / fiberMass D.ref D.visible y := by
  constructor
  · intro h y
    have hstep := push_heatBath_eq (π := D.ref) (b := D.visible) p y
    rw [← repairKernel_def, h] at hstep
    exact hstep
  · intro h
    funext y
    rw [repairKernel_def, push_heatBath_eq p y]
    exact (h y).symm

/-! ## First law, diagonal matrix bridge -/

/-- Diagonal embedding of a classical state or observable into the
matrix carrier of the exact first-law identity, along the constructed
equivalence `Fintype.equivFin`. -/
noncomputable def diagR (v : Ω → ℝ) :
    Matrix (Fin (Fintype.card Ω)) (Fin (Fintype.card Ω)) ℂ :=
  Matrix.diagonal (fun i => (v ((Fintype.equivFin Ω).symm i) : ℂ))

/-- **First law, heat channel of one repair step.** Embedding the
state and the calibrated energy as diagonal matrices along
`Fintype.equivFin`, one repair step is a pure-heat stroke of the exact
matrix identity `first_exact_split`: the observable is held fixed
(`dH = 0`), the cross term vanishes, and the heat increment equals the
classical energy flow `∑ x, (push p K x - p x) * energy x`, the
quantity whose negative the Landauer clause bounds from below, at the
same state, kernel, and calibrated energy. The clause instantiates the
algebraic identity on the diagonal (classical) sector only; no claim
is made that the off-diagonal matrix content derives from the repair
law. Consumes register rows PR-07 and PR-15. -/
theorem first_heat_channel (C : CalibrationData D) (p : Ω → ℝ) :
    (internalEnergy (diagR (push p (repairKernel D))) (diagR C.energy)
        - internalEnergy (diagR p) (diagR C.energy)
      = heatIncrement
          (diagR (push p (repairKernel D)) - diagR p)
          (diagR C.energy))
    ∧ heatIncrement
        (diagR (push p (repairKernel D)) - diagR p) (diagR C.energy)
      = ∑ x, (push p (repairKernel D) x - p x) * C.energy x := by
  constructor
  · have hsplit := firstLaw_split (diagR p)
      (diagR (push p (repairKernel D)) - diagR p) (diagR C.energy)
      (0 : Matrix (Fin (Fintype.card Ω)) (Fin (Fintype.card Ω)) ℂ)
    have hadd : diagR p + (diagR (push p (repairKernel D)) - diagR p)
        = diagR (push p (repairKernel D)) := by abel
    rw [hadd, add_zero] at hsplit
    have hwork : workIncrement (diagR p)
        (0 : Matrix (Fin (Fintype.card Ω)) (Fin (Fintype.card Ω)) ℂ)
        = 0 := by
      unfold workIncrement
      rw [Matrix.mul_zero, Matrix.trace_zero]
      exact Complex.zero_re
    have hcross : (Matrix.trace
        ((diagR (push p (repairKernel D)) - diagR p)
          * (0 : Matrix (Fin (Fintype.card Ω))
              (Fin (Fintype.card Ω)) ℂ))).re = 0 := by
      rw [Matrix.mul_zero, Matrix.trace_zero]
      exact Complex.zero_re
    rw [hwork, hcross] at hsplit
    linarith [hsplit]
  · unfold heatIncrement diagR
    simp only [Matrix.diagonal_sub, Matrix.diagonal_mul_diagonal,
      Matrix.trace_diagonal]
    rw [Complex.re_sum]
    have hterm : ∀ i : Fin (Fintype.card Ω),
        (((push p (repairKernel D) ((Fintype.equivFin Ω).symm i) : ℂ)
            - (p ((Fintype.equivFin Ω).symm i) : ℂ))
          * (C.energy ((Fintype.equivFin Ω).symm i) : ℂ)).re
        = (push p (repairKernel D) ((Fintype.equivFin Ω).symm i)
            - p ((Fintype.equivFin Ω).symm i))
          * C.energy ((Fintype.equivFin Ω).symm i) := by
      intro i
      rw [← Complex.ofReal_sub, ← Complex.ofReal_mul,
        Complex.ofReal_re]
    rw [Finset.sum_congr rfl fun i _ => hterm i]
    exact Equiv.sum_comp (Fintype.equivFin Ω).symm
      (fun x => (push p (repairKernel D) x - p x) * C.energy x)

end FourLawSurface

open FourLawSurface

universe u

/-- **Register row PR-08 attachment: the refinement-uniform gap
family, tied to the calibrated energy.** The declared refinement
family of `LowTemperatureControl` lives on an indexed family of
carriers; the distinguished member `r0` is identified with the surface
state space by the equivalence `e`, and `energy_eq` identifies that
member's energy with the calibrated energy of register row PR-15.
These two fields are the bridge between the PR-08 carrier and the
surface carrier; without them the two halves share no object. -/
structure RefinementAttachment {Ω : Type u} [Fintype Ω] [DecidableEq Ω]
    {B : Type*} [DecidableEq B] (D : RepairLawData Ω B)
    (C : CalibrationData D) where
  /-- The refinement index type. -/
  I : Type
  /-- The refinement order on the index type. -/
  [preorderI : Preorder I]
  /-- The member carriers of the declared family. -/
  X : I → Type u
  /-- Every member carrier is finite. -/
  [fintypeX : ∀ r, Fintype (X r)]
  /-- Every member carrier has decidable equality. -/
  [decEqX : ∀ r, DecidableEq (X r)]
  /-- Every member carrier is inhabited. -/
  [nonemptyX : ∀ r, Nonempty (X r)]
  /-- The member energies of the declared family. -/
  E : ∀ r, X r → ℝ
  /-- The declared refinement family with one uniform positive gap
  lower bound and one uniform cardinality bound: register row
  PR-08. -/
  family : UniformGapRefinement X E
  /-- The distinguished member identified with the surface state
  space. -/
  r0 : I
  /-- The carrier identification of the distinguished member with the
  surface state space. -/
  e : X r0 ≃ Ω
  /-- The distinguished member's energy is the calibrated energy of
  register row PR-15 through the carrier identification. -/
  energy_eq : ∀ x, E r0 x = C.energy (e x)

attribute [instance] RefinementAttachment.preorderI
attribute [instance] RefinementAttachment.fintypeX
attribute [instance] RefinementAttachment.decEqX
attribute [instance] RefinementAttachment.nonemptyX

/-- **The composed thermodynamic antecedent.** One typed bundle whose
three fields are exactly the premise-register rows: PR-07, the
declared repair law, disposition axiomatize, carrying the B12 and B20
bounded no-gos as its justification record; PR-15, the clock and
energy calibration import, shared with the constants lane; PR-08, the
refinement-uniform gap family attached to the same calibrated energy
at the distinguished member. Every clause of `FourLawConclusions`
consumes this bundle. Bundling makes clauses that need only one row
formally depend on the whole bundle; each clause doc comment names the
rows it uses, and the per-row structures stay exported for finer
consumers such as the horizon surface. -/
structure FourLawAntecedent (Ω : Type u) [Fintype Ω] [DecidableEq Ω]
    (B : Type*) [DecidableEq B] where
  /-- Register row PR-07: the declared repair law. -/
  repair : RepairLawData Ω B
  /-- Register row PR-15: the clock and energy calibration of the
  declared reference. -/
  calib : CalibrationData repair
  /-- Register row PR-08: the refinement-uniform gap family attached
  to the calibrated energy. -/
  refinement : RefinementAttachment repair calib

/-- The bundle forces an inhabited state space through the carrier
identification of the distinguished refinement member. -/
theorem FourLawAntecedent.nonempty {Ω : Type u} [Fintype Ω]
    [DecidableEq Ω] {B : Type*} [DecidableEq B]
    (A : FourLawAntecedent Ω B) : Nonempty Ω :=
  A.refinement.e.nonempty_congr.mp
    (A.refinement.nonemptyX A.refinement.r0)

/-- **The composed four-law conclusion package.** One `Prop`-valued
record whose clauses are all stated for the one kernel
`repairKernel A.repair`, the one reference `A.repair.ref`, and the one
calibrated energy `A.calib.energy` of the bundle `A`. The
observation-ledger rows OL-E1 through OL-E4 cite `fourLaws_composed`,
which produces this record from the bundle. The `Nonempty` instance is
derivable from the bundle through `FourLawAntecedent.nonempty`. -/
structure FourLawConclusions {Ω : Type u} [Fintype Ω] [DecidableEq Ω]
    [Nonempty Ω] {B : Type*} [DecidableEq B]
    (A : FourLawAntecedent Ω B) : Prop where
  /-- Kernel receipt (row PR-07): entries are nonnegative. -/
  kernel_nonneg : ∀ x y, 0 ≤ repairKernel A.repair x y
  /-- Kernel receipt (row PR-07): rows are normalized. -/
  kernel_row_sum : ∀ x, ∑ y, repairKernel A.repair x y = 1
  /-- Kernel receipt (row PR-07): the declared reference is
  stationary. -/
  kernel_stationary :
    ∀ y, push A.repair.ref (repairKernel A.repair) y = A.repair.ref y
  /-- Kernel receipt (row PR-07): detailed balance with the declared
  reference. -/
  kernel_detailed_balance : ∀ x y,
    A.repair.ref x * repairKernel A.repair x y
      = A.repair.ref y * repairKernel A.repair y x
  /-- Kernel receipt (row PR-07): one full fibre resampling is
  idempotent. -/
  kernel_idempotent : ∀ x y,
    ∑ z, repairKernel A.repair x z * repairKernel A.repair z y
      = repairKernel A.repair x y
  /-- Kernel receipt (row PR-07): each row is the information
  projection onto its visible fibre. -/
  kernel_row_optimal : ∀ (x : Ω) (r : Ω → ℝ), (∀ y, 0 ≤ r y) →
    (∑ y, r y = 1) →
    (∀ y, A.repair.visible y ≠ A.repair.visible x → r y = 0) →
    kl (repairKernel A.repair x) A.repair.ref ≤ kl r A.repair.ref
  /-- Zeroth law (rows PR-07 and PR-15): the declared reference is the
  Gibbs state of the calibrated energy at the calibrated
  multiplier. -/
  zeroth_ref_eq_gibbs :
    ∀ x, A.repair.ref x = gibbs A.calib.energy A.calib.beta x
  /-- Zeroth law, thermometer form (rows PR-07 and PR-15). -/
  zeroth_thermometer : ∀ (beta1 beta2 : ℝ) (i j : Ω),
    A.calib.energy i ≠ A.calib.energy j →
    gibbs A.calib.energy beta1 = gibbs A.calib.energy beta2 →
    beta1 = beta2
  /-- Zeroth law, multiplier uniqueness (rows PR-07 and PR-15): any
  second calibration row over the same declared reference and the same
  nondegenerate energy carries the bundle's multiplier. -/
  zeroth_multiplier_unique : ∀ C2 : CalibrationData A.repair,
    C2.energy = A.calib.energy →
    ∀ i j : Ω, A.calib.energy i ≠ A.calib.energy j →
    A.calib.beta = C2.beta
  /-- First law, modular form with exact remainder (row PR-07). -/
  first_cap_exact : ∀ p : Ω → ℝ,
    shannon p - shannon A.repair.ref
      = (∑ x, p x * (-Real.log (A.repair.ref x)))
        - (∑ x, A.repair.ref x * (-Real.log (A.repair.ref x)))
        - kl p A.repair.ref
  /-- First law, conservation under pure repair (row PR-07). -/
  first_repair_conserves_fibre_mean : ∀ Q : Ω → ℝ,
    (∀ x y, A.repair.visible x = A.repair.visible y → Q x = Q y) →
    ∀ p : Ω → ℝ,
    ∑ y, push p (repairKernel A.repair) y * Q y = ∑ x, p x * Q x
  /-- First law, heat channel of one repair step (rows PR-07 and
  PR-15): the diagonal instantiation of `first_exact_split` at the
  repair step and the calibrated energy; the refinement field is
  unused by this clause. -/
  first_heat_channel : ∀ p : Ω → ℝ,
    (internalEnergy (diagR (push p (repairKernel A.repair)))
          (diagR A.calib.energy)
        - internalEnergy (diagR p) (diagR A.calib.energy)
      = heatIncrement
          (diagR (push p (repairKernel A.repair)) - diagR p)
          (diagR A.calib.energy))
    ∧ heatIncrement
        (diagR (push p (repairKernel A.repair)) - diagR p)
        (diagR A.calib.energy)
      = ∑ x, (push p (repairKernel A.repair) x - p x)
          * A.calib.energy x
  /-- Second law, contraction form (row PR-07). -/
  second_contraction : ∀ p : Ω → ℝ, (∀ x, 0 ≤ p x) →
    kl (push p (repairKernel A.repair)) A.repair.ref
      ≤ kl p A.repair.ref
  /-- Second law without detailed balance (row PR-07): any stochastic
  kernel preserving the declared reference contracts relative entropy
  to it. -/
  second_contraction_of_stationary :
    ∀ (p : Ω → ℝ) (K : Ω → Ω → ℝ), (∀ x, 0 ≤ p x) →
    (∀ x y, 0 ≤ K x y) → (∀ x, ∑ y, K x y = 1) →
    push A.repair.ref K = A.repair.ref →
    kl (push p K) A.repair.ref ≤ kl p A.repair.ref
  /-- Second law, Clausius form (row PR-07). -/
  second_clausius : ∀ p : Ω → ℝ, (∀ x, 0 ≤ p x) →
    (∑ x, push p (repairKernel A.repair) x
        * (-Real.log (A.repair.ref x)))
      - (∑ x, p x * (-Real.log (A.repair.ref x)))
      ≤ shannon (push p (repairKernel A.repair)) - shannon p
  /-- Integral fluctuation identity (row PR-07). -/
  second_integral_fluctuation : ∀ p : Ω → ℝ, (∀ x, 0 < p x) →
    (∑ x, p x = 1) →
    ∑ x, ∑ y, p x * repairKernel A.repair x y
      * Real.exp (-(sigmaEP A.repair.ref p
          (push p (repairKernel A.repair)) x y)) = 1
  /-- Detailed fluctuation relation, Crooks form (row PR-07). -/
  second_crooks : ∀ p : Ω → ℝ, (∀ x, 0 < p x) → ∀ x y,
    p x * repairKernel A.repair x y
      = Real.exp (sigmaEP A.repair.ref p
          (push p (repairKernel A.repair)) x y)
        * (push p (repairKernel A.repair) y
            * repairKernel A.repair y x)
  /-- Arrow of time, nonnegativity (row PR-07). -/
  arrow_mean_entropy_production_nonneg : ∀ p : Ω → ℝ,
    (∀ x, 0 < p x) →
    0 ≤ ∑ x, ∑ y, p x * repairKernel A.repair x y
      * sigmaEP A.repair.ref p (push p (repairKernel A.repair)) x y
  /-- Arrow of time, exact dissipation identity (row PR-07): the mean
  entropy production of one repair step equals the relative entropy
  from the state to its repaired image. Stated for strictly positive
  states. -/
  arrow_mean_ep_eq_kl_to_repaired : ∀ p : Ω → ℝ, (∀ x, 0 < p x) →
    ∑ x, ∑ y, p x * repairKernel A.repair x y
        * sigmaEP A.repair.ref p (push p (repairKernel A.repair)) x y
      = kl p (push p (repairKernel A.repair))
  /-- Arrow of time, strictness (row PR-07): a state-changing repair
  step produces strictly positive mean entropy production, and zero
  production forces the state onto the repaired manifold. Stated for
  strictly positive states. -/
  arrow_strict : ∀ p : Ω → ℝ, (∀ x, 0 < p x) →
    (push p (repairKernel A.repair) ≠ p ↔
      0 < ∑ x, ∑ y, p x * repairKernel A.repair x y
        * sigmaEP A.repair.ref p (push p (repairKernel A.repair)) x y)
  /-- Zero-dissipation characterization (row PR-07): the fixed states
  of one repair step are exactly the states carrying the reference's
  conditional structure on every visible fibre. -/
  repair_fixed_iff : ∀ p : Ω → ℝ,
    (push p (repairKernel A.repair) = p ↔
      ∀ y, p y = A.repair.ref y * fiberMass p A.repair.visible y
        / fiberMass A.repair.ref A.repair.visible y)
  /-- Third law, excited-mass bound (rows PR-07 and PR-15). -/
  third_excited_mass_bound : ∀ beta E0 Δ : ℝ, 0 ≤ beta → 0 < Δ →
    (∀ x, A.calib.energy x = E0 ∨ E0 + Δ ≤ A.calib.energy x) →
    (Finset.univ.filter (fun x => A.calib.energy x = E0)).Nonempty →
    excitedMass A.calib.energy beta E0
      ≤ ((Finset.univ.filter
            (fun x => A.calib.energy x ≠ E0)).card : ℝ)
        / ((Finset.univ.filter
            (fun x => A.calib.energy x = E0)).card : ℝ)
        * Real.exp (-beta * Δ)
  /-- Third law, quantitative threshold (rows PR-07 and PR-15). -/
  third_excited_mass_threshold : ∀ beta E0 Δ eps : ℝ, 0 ≤ beta →
    0 < Δ → 0 < eps →
    (∀ x, A.calib.energy x = E0 ∨ E0 + Δ ≤ A.calib.energy x) →
    (Finset.univ.filter (fun x => A.calib.energy x = E0)).Nonempty →
    Real.log
        (((Finset.univ.filter
            (fun x => A.calib.energy x ≠ E0)).card : ℝ)
          / ((Finset.univ.filter
              (fun x => A.calib.energy x = E0)).card : ℝ) / eps)
      < beta * Δ →
    excitedMass A.calib.energy beta E0 < eps
  /-- Third law, finite-step unattainability (row PR-07). -/
  third_no_step_extinguishes : ∀ p : Ω → ℝ, (∀ x, 0 < p x) →
    ∀ y, 0 < push p (repairKernel A.repair) y
  /-- Third law, refinement-uniform bound (row PR-08): one explicit
  bound controls every member of the declared refinement family at
  every real inverse temperature. -/
  third_refinement_uniform : ∀ (beta : ℝ) (r : A.refinement.I),
    offMinMass (A.refinement.E r) beta
      ≤ (A.refinement.family.cardBound : ℝ)
        * Real.exp (-beta * A.refinement.family.gapBound)
  /-- Third law, calibrated member (rows PR-08 and PR-15): the
  calibrated energy is the distinguished member of the declared
  refinement family through the carrier identification, and inherits
  the refinement-uniform bound at every real inverse temperature. -/
  third_calibrated_member : ∀ beta : ℝ,
    offMinMass A.calib.energy beta
        = offMinMass (A.refinement.E A.refinement.r0) beta
    ∧ offMinMass A.calib.energy beta
        ≤ (A.refinement.family.cardBound : ℝ)
          * Real.exp (-beta * A.refinement.family.gapBound)
  /-- Landauer bound (rows PR-07 and PR-15). -/
  landauer_bound : ∀ (p : Ω → ℝ) (c : ℝ), (∀ x, 0 ≤ p x) →
    shannon (push p (repairKernel A.repair)) - shannon p ≤ -c →
    (∑ x, p x * A.calib.energy x)
      - (∑ x, push p (repairKernel A.repair) x * A.calib.energy x)
      ≥ c / A.calib.beta

/-- **The composed four-law theorem.** From the one typed antecedent
bundle of register rows PR-07, PR-15, and PR-08, the full four-law
conclusion package holds for the one kernel `repairKernel A.repair`
and the one calibrated energy `A.calib.energy`: kernel receipts, both
zeroth-law forms, the fibre-mean-conservation, modular, and
heat-channel first laws, all five second-law and fluctuation
statements, the exact-identity and strict arrow clauses, the finite
and refinement-uniform third laws with the calibrated energy as the
distinguished family member, and the Landauer bound. The laws are
exact finite theorems under the register at the structural rung;
nothing here asserts a source selection of the reference, the visible
datum, or the kernel, and the emergent rung is gated on the
simulator-export decision recorded in the register row PR-07 notes.
The `Nonempty` instance is derivable from the bundle through
`FourLawAntecedent.nonempty`. -/
theorem fourLaws_composed {Ω : Type u} [Fintype Ω] [DecidableEq Ω]
    [Nonempty Ω] {B : Type*} [DecidableEq B]
    (A : FourLawAntecedent Ω B) : FourLawConclusions A := by
  refine
    { kernel_nonneg := repairKernel_nonneg A.repair
      kernel_row_sum := repairKernel_row_sum A.repair
      kernel_stationary := repairKernel_stationary A.repair
      kernel_detailed_balance := repairKernel_detailedBalance A.repair
      kernel_idempotent := repairKernel_idempotent A.repair
      kernel_row_optimal := repairKernel_row_optimal A.repair
      zeroth_ref_eq_gibbs := ref_eq_gibbs A.repair A.calib
      zeroth_thermometer := zeroth_thermometer A.repair A.calib
      zeroth_multiplier_unique := fun C2 hE i j hij =>
        zeroth_multiplier_unique A.repair A.calib C2 hE.symm i j hij
      first_cap_exact := first_cap_exact A.repair
      first_repair_conserves_fibre_mean :=
        first_repair_conserves_fibre_mean A.repair
      first_heat_channel := first_heat_channel A.repair A.calib
      second_contraction := second_contraction A.repair
      second_contraction_of_stationary :=
        second_contraction_of_stationary A.repair
      second_clausius := second_clausius A.repair
      second_integral_fluctuation :=
        second_integral_fluctuation A.repair
      second_crooks := second_crooks A.repair
      arrow_mean_entropy_production_nonneg :=
        arrow_mean_entropy_production_nonneg A.repair
      arrow_mean_ep_eq_kl_to_repaired :=
        arrow_mean_entropy_production_eq_kl_to_repaired A.repair
      arrow_strict := arrow_strict A.repair
      repair_fixed_iff := repair_fixed_iff A.repair
      third_excited_mass_bound :=
        third_excited_mass_bound A.repair A.calib
      third_excited_mass_threshold :=
        third_excited_mass_threshold A.repair A.calib
      third_no_step_extinguishes :=
        third_no_step_extinguishes A.repair
      third_refinement_uniform := fun beta r =>
        A.refinement.family.uniform_bound beta r
      third_calibrated_member := ?_
      landauer_bound := landauer_bound A.repair A.calib }
  intro beta
  have hE : A.refinement.E A.refinement.r0
      = fun x => A.calib.energy (A.refinement.e x) :=
    funext A.refinement.energy_eq
  have hoff : offMinMass (A.refinement.E A.refinement.r0) beta
      = offMinMass A.calib.energy beta := by
    rw [hE]
    exact offMinMass_equiv A.refinement.e A.calib.energy beta
  refine ⟨hoff.symm, ?_⟩
  rw [← hoff]
  exact A.refinement.family.uniform_bound beta A.refinement.r0

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
#print axioms OPH.Thermodynamics.FourLawSurface.arrow_mean_entropy_production_eq_kl_to_repaired
#print axioms OPH.Thermodynamics.FourLawSurface.arrow_strict
#print axioms OPH.Thermodynamics.FourLawSurface.repair_fixed_iff
#print axioms OPH.Thermodynamics.FourLawSurface.first_heat_channel
#print axioms OPH.Thermodynamics.FourLawAntecedent.nonempty
#print axioms OPH.Thermodynamics.fourLaws_composed
