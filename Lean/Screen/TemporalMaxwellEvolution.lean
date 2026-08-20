import LocalFaceMaxwellAction

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.TemporalMaxwellEvolution

open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction

/-!
# Temporal electric-magnetic evolution and charge-current continuity on the committed carrier

WHAT IS PROVED.  Exact finite real linear algebra on the committed
twelve-port, thirty-seam, twenty-oriented-face complex, reusing the exact
committed objects: the port coboundary `realCoboundary` (d), its boundary
transpose `realBoundary` (∂) under the committed equal-weight pairing, the
committed oriented face curvature `faceCurvature` (C) with `C ∘ d = 0`, its
transpose `faceCodifferential` (Cᵀ), and the local normal operator
`localMaxwellOperator` (CᵀC) of `LocalFaceMaxwellAction`.

Fields carry one declared evolution index `n : ℕ` with a stated staggered
convention: the seam potential `A n` and port potential `φ n` sit at integer
steps; the electric seam field `electricField A φ n = -(A (n+1) - A n) - d (φ n)`
sits on the half step between `n` and `n+1`; the magnetic face field
`magneticField A n = C (A n)` sits at integer steps.  One discrete Maxwell
evolution law is declared, the leapfrog Ampere update
`AmpereEvolution`: `E (n+1) - E n = Cᵀ (B (n+1)) - J n`, with `J n` the
declared seam current driving the E-update across integer step `n+1`.
The Faraday law `B (n+1) - B n = -C (E n)` is a theorem of the definitions,
not a second postulate.

The receipts are conditional on the declared Ampere update: (1)
charge-current continuity, both directions: along the
evolution the Gauss constraint `∂ (E n) = ρ n` propagates from step `n` to
`n+1` exactly when `ρ (n+1) - ρ n + ∂ (J n) = 0` (`gauss_step_iff`), and
with the continuity equation it propagates from the initial step to every
step (`gauss_propagation`); (2) an exact per-step identity for the
staggered quadratic form
`(1/2)‖E n‖² + (1/2)⟨B n, B (n+1)⟩` with source work term
`-(1/2)⟨E n + E (n+1), J n⟩` (`energy_balance`), hence exact conservation
of this staggered quadratic form for `J = 0` (`energy_conserved`).  The form
is not proved positive and its conservation is not a stability theorem;
(3) full time-dependent gauge invariance:
`A ↦ A + d (χ n)`, `φ ↦ φ - (χ (n+1) - χ n)` leaves `E`, `B`, the evolution
law, and the quadratic form invariant; (4) the static join: a constant-in-time
zero-current solution has zero magnetic field, a forced-neutral load, and
electric field exactly the canonical Coulomb field `realCoulombField ρ` of
`DiscreteCoulombGreen`/`PositionSpaceMaxwellAction` (`static_join`), and the
canonical Green potential supplies such a static solution for every neutral
load (`coulombStaticSolution`); (5) the exact discrete wave law: zero-current
evolution forces
`A (n+2) - 2 A (n+1) + A n + CᵀC (A (n+1)) + d (φ (n+1) - φ n) = 0`
(`wave_law`).  One composed receipt
`temporalMaxwellEvolution_receipt` carries all clauses from the single typed
antecedent bundle `MaxwellEvolutionBundle`, and `demoBundle` is an explicit
nonstatic inhabitant: zero initial data kicked by one face-boundary cycle,
with `B 0 = 0` and `B 1` equal to three on face zero.

TYPED SOURCE FIREWALL.  The port load `ρ : ℕ → Fin 12 → ℝ` and the seam
current `J : ℕ → Fin 30 → ℝ` stay distinct types throughout.  No definition
identifies them; the only bridge is the proved continuity equivalence
`ρ (n+1) - ρ n + ∂ (J n) = 0` of `gauss_step_iff`.

WHAT IS NOT PROVED.  The step index `n` is a declared evolution parameter,
not physical time; clock and calibration attachment stays open under PR-15.
The sources `ρ` and `J` are declared inputs, not source-produced currents;
the source gauge-field/current/action attachment stays open under PR-54.
No photon, propagation speed, Lorentz covariance, continuum limit, or
laboratory readout is claimed; the physical propagation, frame, and
comparison attachment stays open under PR-53.  The exact unit-step Ampere
update is declared rather than derived from an OPH repair or variational law
and is consumed under PR-66.  The equal-weight pairing coincides with the
PR-20 equal seam-counting selection.  No positive energy, Courant bound, or
stability estimate is proved.  The exact spectrum receipt for `C Cᵀ`
contains eigenvalues above the unit-step leapfrog stability threshold, so a
timestep-scaled update and stability proof remain viable required work under
PR-53 rather than a no-go against temporal propagation.

Axiom audit.  Every proof composes the committed receipts with exact real
linear algebra; the module adds no project axiom and uses no native
decision procedure.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`.
-/

noncomputable section

/-! ## Time-indexed fields on the committed carrier -/

/-- Electric seam field on the half step between `n` and `n+1`:
`E n = -(A (n+1) - A n) - d (φ n)`. -/
def electricField (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    Fin 30 → ℝ :=
  -(A (n + 1) - A n) - realCoboundary (φ n)

/-- Magnetic face field at integer step `n`: `B n = C (A n)`. -/
def magneticField (A : ℕ → Fin 30 → ℝ) (n : ℕ) : Fin 20 → ℝ :=
  faceCurvature (A n)

/-- The declared leapfrog Ampere update: the E-update across integer step
`n+1` is driven by the transpose face incidence applied to `B (n+1)` minus
the declared seam current `J n`. -/
def AmpereEvolution (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) : Prop :=
  ∀ n : ℕ, electricField A φ (n + 1) - electricField A φ n =
    faceCodifferential (magneticField A (n + 1)) - J n

/-- The face curvature of the electric field is the backward magnetic
difference; the scalar-potential term drops out by `C ∘ d = 0`. -/
theorem faceCurvature_electricField (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    faceCurvature (electricField A φ n) =
      magneticField A n - magneticField A (n + 1) := by
  unfold electricField magneticField
  rw [map_sub, faceCurvature_coboundary, sub_zero, map_neg, map_sub, neg_sub]

/-- Faraday law as a theorem of the staggered definitions: no second
postulate is declared. -/
theorem faraday_law (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    magneticField A (n + 1) - magneticField A n =
      -faceCurvature (electricField A φ n) := by
  rw [faceCurvature_electricField, neg_sub]

/-! ## Charge-current continuity, both directions -/

/-- The boundary of the transpose face incidence vanishes: the exact
transpose of the committed `C ∘ d = 0`. -/
theorem realBoundary_faceCodifferential (F : Fin 20 → ℝ) :
    realBoundary (faceCodifferential F) = 0 := by
  have hpair : realPortInner (realBoundary (faceCodifferential F))
      (realBoundary (faceCodifferential F)) = 0 := by
    rw [← realCoboundary_boundary_adjoint,
      ← faceCurvature_codifferential_adjoint, faceCurvature_coboundary]
    simp [faceInner]
  exact (realPortInner_self_eq_zero_iff _).mp hpair

/-- **Continuity, both directions.**  Along the declared evolution, the
Gauss constraint at step `n` propagates to step `n+1` exactly when the
discrete continuity equation `ρ (n+1) - ρ n + ∂ (J n) = 0` holds.  The port
load `ρ` and the seam current `J` are distinct types; this equivalence is
their only bridge. -/
theorem gauss_step_iff (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolution A φ J) (n : ℕ)
    (hn : realBoundary (electricField A φ n) = ρ n) :
    realBoundary (electricField A φ (n + 1)) = ρ (n + 1) ↔
      ρ (n + 1) - ρ n + realBoundary (J n) = 0 := by
  have h : realBoundary (electricField A φ (n + 1) - electricField A φ n) =
      realBoundary (faceCodifferential (magneticField A (n + 1)) - J n) := by
    rw [hAmp n]
  rw [map_sub, map_sub, realBoundary_faceCodifferential, hn, zero_sub] at h
  have hkey : realBoundary (electricField A φ (n + 1)) =
      ρ n - realBoundary (J n) := by
    funext p
    have hp := congrFun h p
    simp only [Pi.sub_apply, Pi.neg_apply] at hp
    simp only [Pi.sub_apply]
    linarith
  rw [hkey]
  constructor
  · intro heq
    funext p
    have hp := congrFun heq p
    simp only [Pi.sub_apply] at hp
    simp only [Pi.add_apply, Pi.sub_apply, Pi.zero_apply]
    linarith
  · intro heq
    funext p
    have hp := congrFun heq p
    simp only [Pi.add_apply, Pi.sub_apply, Pi.zero_apply] at hp
    simp only [Pi.sub_apply]
    linarith

/-- Under the continuity equation, the Gauss constraint propagates from the
initial step to every step. -/
theorem gauss_propagation (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolution A φ J)
    (h0 : realBoundary (electricField A φ 0) = ρ 0)
    (hcont : ∀ n : ℕ, ρ (n + 1) - ρ n + realBoundary (J n) = 0) :
    ∀ n : ℕ, realBoundary (electricField A φ n) = ρ n := by
  intro n
  induction n with
  | zero => exact h0
  | succ m ih => exact (gauss_step_iff A φ J ρ hAmp m ih).mpr (hcont m)

/-! ## The exact staggered quadratic-form identity -/

theorem realSeamInner_sub_right (a b c : Fin 30 → ℝ) :
    realSeamInner a (b - c) = realSeamInner a b - realSeamInner a c := by
  unfold realSeamInner
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun e _ ↦ by
    simp only [Pi.sub_apply]
    ring

theorem faceInner_sub_left (F G H : Fin 20 → ℝ) :
    faceInner (F - G) H = faceInner F H - faceInner G H := by
  unfold faceInner
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun f _ ↦ by
    simp only [Pi.sub_apply]
    ring

theorem realSeamEnergy_sub_eq_inner (x y : Fin 30 → ℝ) :
    realSeamEnergy y - realSeamEnergy x = realSeamInner (y + x) (y - x) := by
  unfold realSeamEnergy realSeamInner
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun e _ ↦ by
    simp only [Pi.add_apply, Pi.sub_apply]
    ring

/-- Staggered discrete quadratic energy form with the committed equal-weight
pairing: `(1/2)‖E n‖² + (1/2)⟨B n, B (n+1)⟩`.  The magnetic cross term at
adjacent integer steps is the exactly conserved form for the leapfrog
update.  No positivity or coercivity property is asserted. -/
def fieldEnergy (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ) : ℝ :=
  (1 / 2) * realSeamEnergy (electricField A φ n) +
    (1 / 2) * faceInner (magneticField A n) (magneticField A (n + 1))

/-- **Exact per-step quadratic-form balance.**  The form's change across integer
step `n+1` is exactly minus the source work of `J n` against the
time-averaged electric field. -/
theorem energy_balance (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (hAmp : AmpereEvolution A φ J) (n : ℕ) :
    fieldEnergy A φ (n + 1) = fieldEnergy A φ n -
      (1 / 2) * realSeamInner
        (electricField A φ n + electricField A φ (n + 1)) (J n) := by
  have hdiff : realSeamEnergy (electricField A φ (n + 1)) -
      realSeamEnergy (electricField A φ n) =
      realSeamInner (electricField A φ (n + 1) + electricField A φ n)
        (electricField A φ (n + 1) - electricField A φ n) :=
    realSeamEnergy_sub_eq_inner _ _
  rw [hAmp n, realSeamInner_sub_right] at hdiff
  have hadj : realSeamInner
      (electricField A φ (n + 1) + electricField A φ n)
      (faceCodifferential (magneticField A (n + 1))) =
      faceInner (faceCurvature
        (electricField A φ (n + 1) + electricField A φ n))
        (magneticField A (n + 1)) :=
    (faceCurvature_codifferential_adjoint _ _).symm
  have hCsum : faceCurvature
      (electricField A φ (n + 1) + electricField A φ n) =
      magneticField A n - magneticField A (n + 1 + 1) := by
    rw [map_add, faceCurvature_electricField, faceCurvature_electricField]
    abel
  rw [hCsum, faceInner_sub_left] at hadj
  rw [hadj] at hdiff
  have hflip : faceInner (magneticField A (n + 1 + 1))
      (magneticField A (n + 1)) =
      faceInner (magneticField A (n + 1)) (magneticField A (n + 1 + 1)) :=
    faceInner_comm _ _
  have hEflip : realSeamInner
      (electricField A φ (n + 1) + electricField A φ n) (J n) =
      realSeamInner
        (electricField A φ n + electricField A φ (n + 1)) (J n) := by
    rw [add_comm (electricField A φ (n + 1)) (electricField A φ n)]
  rw [hflip, hEflip] at hdiff
  unfold fieldEnergy
  linarith [hdiff]

/-- Exact conservation of the staggered quadratic form for zero-current
evolution.  This statement does not imply positivity or stability. -/
theorem energy_conserved (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolution A φ (fun _ ↦ 0)) (n : ℕ) :
    fieldEnergy A φ n = fieldEnergy A φ 0 := by
  induction n with
  | zero => rfl
  | succ m ih =>
    have h : fieldEnergy A φ (m + 1) = fieldEnergy A φ m -
        (1 / 2) * realSeamInner
          (electricField A φ m + electricField A φ (m + 1)) 0 :=
      energy_balance A φ (fun _ ↦ 0) hAmp m
    rw [realSeamInner_zero_right, mul_zero, sub_zero] at h
    rw [h]
    exact ih

/-! ## Full time-dependent gauge invariance -/

/-- Time-dependent vector gauge transformation `A n ↦ A n + d (χ n)`. -/
def gaugeTransformA (A : ℕ → Fin 30 → ℝ) (χ : ℕ → Fin 12 → ℝ) :
    ℕ → Fin 30 → ℝ :=
  fun n ↦ A n + realCoboundary (χ n)

/-- Matching scalar gauge transformation `φ n ↦ φ n - (χ (n+1) - χ n)`. -/
def gaugeTransformPhi (φ : ℕ → Fin 12 → ℝ) (χ : ℕ → Fin 12 → ℝ) :
    ℕ → Fin 12 → ℝ :=
  fun n ↦ φ n - (χ (n + 1) - χ n)

theorem electricField_gauge_invariant (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (χ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    electricField (gaugeTransformA A χ) (gaugeTransformPhi φ χ) n =
      electricField A φ n := by
  unfold electricField gaugeTransformA gaugeTransformPhi
  funext e
  simp only [map_sub, Pi.sub_apply, Pi.neg_apply, Pi.add_apply,
    realCoboundary_apply]
  ring

theorem magneticField_gauge_invariant (A : ℕ → Fin 30 → ℝ)
    (χ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    magneticField (gaugeTransformA A χ) n = magneticField A n := by
  unfold magneticField gaugeTransformA
  rw [map_add, faceCurvature_coboundary, add_zero]

theorem ampereEvolution_gauge_invariant (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (χ : ℕ → Fin 12 → ℝ) :
    AmpereEvolution (gaugeTransformA A χ) (gaugeTransformPhi φ χ) J ↔
      AmpereEvolution A φ J := by
  unfold AmpereEvolution
  constructor
  · intro h n
    have hn := h n
    rwa [electricField_gauge_invariant, electricField_gauge_invariant,
      magneticField_gauge_invariant] at hn
  · intro h n
    rw [electricField_gauge_invariant, electricField_gauge_invariant,
      magneticField_gauge_invariant]
    exact h n

theorem fieldEnergy_gauge_invariant (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (χ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    fieldEnergy (gaugeTransformA A χ) (gaugeTransformPhi φ χ) n =
      fieldEnergy A φ n := by
  unfold fieldEnergy
  rw [electricField_gauge_invariant, magneticField_gauge_invariant,
    magneticField_gauge_invariant]

/-! ## The static join to the committed Coulomb sector -/

/-- A constant-in-time history has constant electric field `-(d (φ 0))`. -/
theorem static_electricField (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (hA : ∀ n, A n = A 0) (hφ : ∀ n, φ n = φ 0) (n : ℕ) :
    electricField A φ n = -realCoboundary (φ 0) := by
  unfold electricField
  rw [hA (n + 1), hA n, hφ n, sub_self, neg_zero, zero_sub]

/-- **Static join, uniqueness direction.**  A constant-in-time zero-current
solution with Gauss load `ρ` has forced-neutral `ρ`, zero
magnetic field, and electric field exactly the canonical Coulomb field of
the committed Green solution. -/
theorem static_join (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (ρ : Fin 12 → ℝ) (hA : ∀ n, A n = A 0) (hφ : ∀ n, φ n = φ 0)
    (hAmp : AmpereEvolution A φ (fun _ ↦ 0))
    (hGauss : realBoundary (electricField A φ 0) = ρ) :
    (∑ p : Fin 12, ρ p) = 0 ∧
      (∀ n, magneticField A n = 0) ∧
      (∀ n, electricField A φ n = realCoulombField ρ) := by
  have hE : ∀ n, electricField A φ n = -realCoboundary (φ 0) :=
    static_electricField A φ hA hφ
  have hneutral : (∑ p : Fin 12, ρ p) = 0 := by
    rw [← hGauss]
    exact realBoundary_total _
  have hCT : faceCodifferential (magneticField A 1) = 0 := by
    have h := hAmp 0
    rw [hE (0 + 1), hE 0, sub_self] at h
    have h' : (0 : Fin 30 → ℝ) =
        faceCodifferential (magneticField A 1) - 0 := h
    rw [sub_zero] at h'
    exact h'.symm
  have hker : A 1 ∈ LinearMap.ker localMaxwellOperator := by
    rw [LinearMap.mem_ker]
    exact hCT
  rw [ker_localMaxwellOperator] at hker
  have hB1 : faceCurvature (A 1) = 0 := LinearMap.mem_ker.mp hker
  have hBn : ∀ n, magneticField A n = 0 := by
    intro n
    unfold magneticField
    rw [hA n, ← hA 1]
    exact hB1
  have hL : realLaplacian (-(φ 0)) = ρ := by
    have h := hGauss
    rw [hE 0, map_neg] at h
    rw [map_neg, realLaplacian_apply]
    exact h
  have hG : realLaplacian (greenPotential ρ) = ρ :=
    greenPotential_stationary ρ hneutral
  have hker2 : realLaplacian (-(φ 0) - greenPotential ρ) = 0 := by
    rw [map_sub, hL, hG, sub_self]
  obtain ⟨c, hc⟩ := (realLaplacian_eq_zero_iff _).mp hker2
  have hsplit : -(φ 0) = greenPotential ρ + fun _ ↦ c := by
    have h1 := sub_eq_iff_eq_add.mp hc
    rw [h1, add_comm]
  have hgrad : realCoboundary (-(φ 0)) = realCoulombField ρ := by
    rw [hsplit, map_add, greenPotential_field,
      (realCoboundary_eq_zero_iff _).mpr ⟨c, rfl⟩, add_zero]
  have hEC : ∀ n, electricField A φ n = realCoulombField ρ := by
    intro n
    rw [hE n, ← hgrad, map_neg]
  exact ⟨hneutral, hBn, hEC⟩

/-- **Static join, existence direction.**  For every neutral load the
canonical Green potential supplies a static solution of the evolution whose
electric field is exactly the committed Coulomb field and whose Gauss load
is exactly `ρ`. -/
theorem coulombStaticSolution (ρ : Fin 12 → ℝ)
    (hρ : (∑ p : Fin 12, ρ p) = 0) :
    AmpereEvolution (fun _ ↦ 0) (fun _ ↦ -greenPotential ρ) (fun _ ↦ 0) ∧
      (∀ n, electricField (fun _ ↦ 0) (fun _ ↦ -greenPotential ρ) n =
        realCoulombField ρ) ∧
      (∀ n, realBoundary
        (electricField (fun _ ↦ 0) (fun _ ↦ -greenPotential ρ) n) = ρ) ∧
      (∀ n, magneticField (fun _ ↦ (0 : Fin 30 → ℝ)) n = 0) := by
  have hE : ∀ n, electricField (fun _ ↦ 0) (fun _ ↦ -greenPotential ρ) n =
      realCoulombField ρ := by
    intro n
    show -((0 : Fin 30 → ℝ) - 0) - realCoboundary (-greenPotential ρ) =
      realCoulombField ρ
    rw [sub_self, neg_zero, zero_sub, map_neg, greenPotential_field, neg_neg]
  have hB : ∀ n, magneticField (fun _ ↦ (0 : Fin 30 → ℝ)) n = 0 := by
    intro n
    show faceCurvature 0 = 0
    exact map_zero _
  refine ⟨?_, hE, ?_, hB⟩
  · intro n
    rw [hE (n + 1), hE n, sub_self, hB (n + 1), map_zero]
    show (0 : Fin 30 → ℝ) = 0 - 0
    rw [sub_zero]
  · intro n
    rw [hE n]
    exact realCoulombField_gauss ρ hρ

/-! ## The exact discrete wave law -/

/-- **Wave law.**  Source-free evolution forces the exact second-order
recursion for `A` with the committed local normal operator `CᵀC`, modulo
the explicit scalar-gradient term. -/
theorem wave_law (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolution A φ (fun _ ↦ 0)) (n : ℕ) :
    A (n + 2) - (2 : ℝ) • A (n + 1) + A n +
      localMaxwellOperator (A (n + 1)) +
      realCoboundary (φ (n + 1) - φ n) = 0 := by
  have h := hAmp n
  have hlm : faceCodifferential (magneticField A (n + 1)) =
      localMaxwellOperator (A (n + 1)) := rfl
  have hE1 : electricField A φ (n + 1) =
      -(A (n + 2) - A (n + 1)) - realCoboundary (φ (n + 1)) := rfl
  have hE0 : electricField A φ n =
      -(A (n + 1) - A n) - realCoboundary (φ n) := rfl
  rw [hlm, hE1, hE0] at h
  funext e
  have he := congrFun h e
  simp only [map_sub, Pi.sub_apply, Pi.neg_apply, Pi.add_apply,
    Pi.smul_apply, smul_eq_mul, realCoboundary_apply, Pi.zero_apply] at he ⊢
  linarith

/-! ## The one citable composed receipt -/

/-- The single typed antecedent bundle: declared histories, declared
sources, the declared Ampere evolution law, the initial Gauss constraint,
and the continuity equation.  The port load `rho` and the seam current `J`
are separate fields with no conversion map. -/
structure MaxwellEvolutionBundle where
  A : ℕ → Fin 30 → ℝ
  phi : ℕ → Fin 12 → ℝ
  rho : ℕ → Fin 12 → ℝ
  J : ℕ → Fin 30 → ℝ
  ampere : AmpereEvolution A phi J
  gauss_init : realBoundary (electricField A phi 0) = rho 0
  continuity : ∀ n : ℕ, rho (n + 1) - rho n + realBoundary (J n) = 0

/-- **The temporal Maxwell evolution receipt (issue #733).**  One typed
conjunction from the single antecedent bundle: Gauss propagation to every
step; the step-local continuity equivalence in both directions for every
candidate load history; the Faraday law of the staggered definitions; the
exact per-step quadratic-form balance with source work term; full
time-dependent gauge invariance of `E`, `B`, the evolution law, and the
quadratic form; exact conservation and the exact discrete wave law in the
zero-current case; and the static join sending every constant zero-current
solution with its declared neutral load to
the committed canonical Coulomb field with forced neutrality and zero
magnetic field.  The step index is a declared evolution parameter and the
sources are declared; PR-66 is consumed and PR-15, PR-53, and PR-54 stay
open.  Positivity and stability are not clauses of the receipt. -/
theorem temporalMaxwellEvolution_receipt (S : MaxwellEvolutionBundle) :
    (∀ n, realBoundary (electricField S.A S.phi n) = S.rho n)
    ∧ (∀ (ρ' : ℕ → Fin 12 → ℝ) (n : ℕ),
        realBoundary (electricField S.A S.phi n) = ρ' n →
        (realBoundary (electricField S.A S.phi (n + 1)) = ρ' (n + 1) ↔
          ρ' (n + 1) - ρ' n + realBoundary (S.J n) = 0))
    ∧ (∀ n, magneticField S.A (n + 1) - magneticField S.A n =
        -faceCurvature (electricField S.A S.phi n))
    ∧ (∀ n, fieldEnergy S.A S.phi (n + 1) = fieldEnergy S.A S.phi n -
        (1 / 2) * realSeamInner
          (electricField S.A S.phi n + electricField S.A S.phi (n + 1))
          (S.J n))
    ∧ (∀ χ : ℕ → Fin 12 → ℝ,
        (∀ n, electricField (gaugeTransformA S.A χ)
            (gaugeTransformPhi S.phi χ) n = electricField S.A S.phi n)
        ∧ (∀ n, magneticField (gaugeTransformA S.A χ) n =
            magneticField S.A n)
        ∧ AmpereEvolution (gaugeTransformA S.A χ)
            (gaugeTransformPhi S.phi χ) S.J
        ∧ (∀ n, fieldEnergy (gaugeTransformA S.A χ)
            (gaugeTransformPhi S.phi χ) n = fieldEnergy S.A S.phi n))
    ∧ (S.J = (fun _ ↦ 0) →
        (∀ n, fieldEnergy S.A S.phi n = fieldEnergy S.A S.phi 0)
        ∧ (∀ n, S.A (n + 2) - (2 : ℝ) • S.A (n + 1) + S.A n +
            localMaxwellOperator (S.A (n + 1)) +
            realCoboundary (S.phi (n + 1) - S.phi n) = 0))
    ∧ ((∀ n, S.A n = S.A 0) → (∀ n, S.phi n = S.phi 0) →
        S.J = (fun _ ↦ 0) →
        ((∑ p : Fin 12, S.rho 0 p) = 0
        ∧ (∀ n, magneticField S.A n = 0)
        ∧ (∀ n, electricField S.A S.phi n = realCoulombField (S.rho 0)))) := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact gauss_propagation S.A S.phi S.J S.rho S.ampere S.gauss_init
      S.continuity
  · intro ρ' n hn
    exact gauss_step_iff S.A S.phi S.J ρ' S.ampere n hn
  · intro n
    exact faraday_law S.A S.phi n
  · exact energy_balance S.A S.phi S.J S.ampere
  · intro χ
    exact ⟨electricField_gauge_invariant S.A S.phi χ,
      magneticField_gauge_invariant S.A χ,
      (ampereEvolution_gauge_invariant S.A S.phi S.J χ).mpr S.ampere,
      fieldEnergy_gauge_invariant S.A S.phi χ⟩
  · intro hJ
    have hAmp0 : AmpereEvolution S.A S.phi (fun _ ↦ 0) := by
      rw [← hJ]
      exact S.ampere
    exact ⟨energy_conserved S.A S.phi hAmp0, wave_law S.A S.phi hAmp0⟩
  · intro hA hφ hJ
    have hAmp0 : AmpereEvolution S.A S.phi (fun _ ↦ 0) := by
      rw [← hJ]
      exact S.ampere
    exact static_join S.A S.phi (S.rho 0) hA hφ hAmp0 S.gauss_init

/-! ## An explicit nonstatic inhabitant -/

/-- Kick field: the signed boundary of committed face zero, a conserved
seam cycle. -/
def demoInitial : Fin 30 → ℝ :=
  faceCodifferential (fun f ↦ if f = 0 then 1 else 0)

theorem demoInitial_apply (e : Fin 30) : demoInitial e = faceIncidenceR 0 e := by
  unfold demoInitial
  rw [faceCodifferential_apply]
  simp [mul_ite, mul_one, mul_zero, Finset.sum_ite_eq']

theorem demoInitial_boundary : realBoundary demoInitial = 0 :=
  realBoundary_faceCodifferential _

/-- Concrete zero-current history: zero initial data kicked by the face-zero
boundary cycle, then evolved by the leapfrog recursion itself. -/
def demoA : ℕ → Fin 30 → ℝ
  | 0 => 0
  | 1 => demoInitial
  | n + 2 => (2 : ℝ) • demoA (n + 1) - demoA n -
      localMaxwellOperator (demoA (n + 1))

def demoPhi : ℕ → Fin 12 → ℝ := fun _ ↦ 0
def demoRho : ℕ → Fin 12 → ℝ := fun _ ↦ 0
def demoJ : ℕ → Fin 30 → ℝ := fun _ ↦ 0

theorem demo_electricField (n : ℕ) :
    electricField demoA demoPhi n = -(demoA (n + 1) - demoA n) := by
  show -(demoA (n + 1) - demoA n) - realCoboundary 0 =
    -(demoA (n + 1) - demoA n)
  rw [map_zero, sub_zero]

theorem demo_ampere : AmpereEvolution demoA demoPhi demoJ := by
  intro n
  rw [demo_electricField (n + 1), demo_electricField n]
  have hlm : faceCodifferential (magneticField demoA (n + 1)) =
      localMaxwellOperator (demoA (n + 1)) := rfl
  rw [hlm, show demoJ n = 0 from rfl, sub_zero]
  have hidx : demoA (n + 1 + 1) = (2 : ℝ) • demoA (n + 1) - demoA n -
      localMaxwellOperator (demoA (n + 1)) := by
    simp only [demoA]
  rw [hidx]
  funext e
  simp only [Pi.sub_apply, Pi.neg_apply, Pi.smul_apply, smul_eq_mul]
  ring

theorem demo_gauss_init :
    realBoundary (electricField demoA demoPhi 0) = demoRho 0 := by
  have hE0 : electricField demoA demoPhi 0 = -demoInitial := by
    rw [demo_electricField 0,
      show demoA (0 + 1) = demoInitial from by simp only [demoA],
      show demoA 0 = (0 : Fin 30 → ℝ) from by simp only [demoA], sub_zero]
  rw [hE0, map_neg, demoInitial_boundary, neg_zero]
  rfl

theorem demo_continuity :
    ∀ n : ℕ, demoRho (n + 1) - demoRho n + realBoundary (demoJ n) = 0 := by
  intro n
  show (0 : Fin 12 → ℝ) - 0 + realBoundary 0 = 0
  rw [map_zero]
  simp

/-- The explicit inhabitant of the antecedent bundle: the evolution is
nonempty. -/
def demoBundle : MaxwellEvolutionBundle where
  A := demoA
  phi := demoPhi
  rho := demoRho
  J := demoJ
  ampere := demo_ampere
  gauss_init := demo_gauss_init
  continuity := demo_continuity

theorem demoA_one_apply_zero : demoA 1 0 = 1 := by
  rw [show demoA 1 = demoInitial from by simp only [demoA],
    demoInitial_apply]
  have hz : faceIncidenceZ 0 0 = 1 := by decide
  show ((faceIncidenceZ 0 0 : ℤ) : ℝ) = 1
  rw [hz]
  norm_num

/-- The inhabitant is nonstatic: the seam potential moves on the first
step. -/
theorem demo_nonstatic : demoA 1 ≠ demoA 0 := by
  intro h
  have h0 := congrFun h 0
  rw [demoA_one_apply_zero] at h0
  have hz : demoA 0 0 = 0 := by
    rw [show demoA 0 = (0 : Fin 30 → ℝ) from by simp only [demoA]]
    rfl
  rw [hz] at h0
  exact one_ne_zero h0

theorem demo_magnetic_start : magneticField demoA 0 = 0 := by
  show faceCurvature (demoA 0) = 0
  rw [show demoA 0 = (0 : Fin 30 → ℝ) from by simp only [demoA], map_zero]

/-- The magnetic field is genuinely excited: it is three on face zero at
step one, against zero at step zero. -/
theorem demo_magnetic_step_one : magneticField demoA 1 0 = 3 := by
  have h1 : demoA 1 = demoInitial := by simp only [demoA]
  show faceCurvature (demoA 1) 0 = 3
  rw [h1, faceCurvature_apply]
  have hterm : ∀ e : Fin 30, faceIncidenceR 0 e * demoInitial e =
      ((faceIncidenceZ 0 e * faceIncidenceZ 0 e : ℤ) : ℝ) := by
    intro e
    rw [demoInitial_apply, Int.cast_mul]
    rfl
  rw [Finset.sum_congr rfl fun e _ ↦ hterm e, ← Int.cast_sum,
    show (∑ e : Fin 30, faceIncidenceZ 0 e * faceIncidenceZ 0 e) = 3 from
      by decide]
  norm_num

theorem demoBundle_nonvacuous :
    demoBundle.A 1 ≠ demoBundle.A 0 ∧
      magneticField demoBundle.A 0 = 0 ∧
      magneticField demoBundle.A 1 0 = 3 :=
  ⟨demo_nonstatic, demo_magnetic_start, demo_magnetic_step_one⟩

end

end OPH.TemporalMaxwellEvolution

/- Axiom audit: the committed receipts and exact real linear algebra only.
Expected axioms per line: at most `propext`, `Classical.choice`,
`Quot.sound`.  No native decision procedure is used. -/

#print axioms OPH.TemporalMaxwellEvolution.faceCurvature_electricField
#print axioms OPH.TemporalMaxwellEvolution.faraday_law
#print axioms OPH.TemporalMaxwellEvolution.realBoundary_faceCodifferential
#print axioms OPH.TemporalMaxwellEvolution.gauss_step_iff
#print axioms OPH.TemporalMaxwellEvolution.gauss_propagation
#print axioms OPH.TemporalMaxwellEvolution.energy_balance
#print axioms OPH.TemporalMaxwellEvolution.energy_conserved
#print axioms OPH.TemporalMaxwellEvolution.electricField_gauge_invariant
#print axioms OPH.TemporalMaxwellEvolution.magneticField_gauge_invariant
#print axioms OPH.TemporalMaxwellEvolution.ampereEvolution_gauge_invariant
#print axioms OPH.TemporalMaxwellEvolution.fieldEnergy_gauge_invariant
#print axioms OPH.TemporalMaxwellEvolution.static_join
#print axioms OPH.TemporalMaxwellEvolution.coulombStaticSolution
#print axioms OPH.TemporalMaxwellEvolution.wave_law
#print axioms OPH.TemporalMaxwellEvolution.temporalMaxwellEvolution_receipt
#print axioms OPH.TemporalMaxwellEvolution.demo_ampere
#print axioms OPH.TemporalMaxwellEvolution.demoBundle_nonvacuous
