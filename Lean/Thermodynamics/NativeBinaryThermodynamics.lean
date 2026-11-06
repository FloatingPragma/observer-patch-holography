import FourLawAdequacySurface

set_option autoImplicit false

namespace OPH.Thermodynamics.NativeBinary

/-!
# Conditional resampling for finite binary swap repairs

The fair mixture of holding a configuration and applying an involution is
exactly conditional resampling from counting measure on its one- or
two-point orbits. On a nonempty finite carrier, normalized counting is a
faithful stationary reference. The existing conditional-repair theorem then
gives relative-entropy decrease and Shannon-entropy increase, as well as
optimality of each repair row within its orbit.

For a fixed-particle-number binary native carrier, endpoint exchange is
such an involution. The separate executable certificate binds its actual
native fair-coin branches to this hold/exchange kernel. This theorem does
not equate the full native quantum instrument with conditional resampling.
A single seam determines uniformity within each orbit, not masses between
orbits; uniqueness of a global reference needs connectivity of the full
repair family. The counting reference has constant modular energy and
supplies no physical energy, calibrated temperature, or physical clock.
-/

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]

noncomputable section

/-- The population transition produced by a fair hold/exchange coin. -/
def fairSwap (τ : Ω → Ω) (x y : Ω) : ℝ :=
  ((if y = x then 1 else 0) + (if y = τ x then 1 else 0)) / 2

omit [Fintype Ω] in
theorem fairSwap_nonneg (τ : Ω → Ω) (x y : Ω) : 0 ≤ fairSwap τ x y := by
  unfold fairSwap
  positivity

/-- The complete invariant of one involution orbit. -/
def orbitLabel (τ : Ω → Ω) (x : Ω) : Finset Ω := {x, τ x}

omit [Fintype Ω] in
theorem orbitLabel_eq_iff {τ : Ω → Ω} (hτ : Function.Involutive τ) (x y : Ω) :
    orbitLabel τ y = orbitLabel τ x ↔ y = x ∨ y = τ x := by
  constructor
  · intro h
    have hy : y ∈ orbitLabel τ x := by rw [← h]; simp [orbitLabel]
    simpa [orbitLabel] using hy
  · rintro (rfl | rfl)
    · rfl
    · simp [orbitLabel, hτ x, Finset.pair_comm]

theorem fiberMass_constant {τ : Ω → Ω} (hτ : Function.Involutive τ)
    (r : ℝ) (x : Ω) :
    fiberMass (fun _ : Ω => r) (orbitLabel τ) x =
      if τ x = x then r else r + r := by
  have hfilter : Finset.univ.filter (fun y => orbitLabel τ y = orbitLabel τ x) =
      orbitLabel τ x := by
    ext y
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    rw [orbitLabel_eq_iff hτ x y]
    simp [orbitLabel]
  rw [fiberMass, hfilter]
  by_cases hx : τ x = x
  · simp [orbitLabel, hx]
  · simp [orbitLabel, hx, Ne.symm hx, two_mul]

/-- The kernel determines conditional uniformity on each orbit for every
positive scaling of counting measure. -/
theorem fairSwap_eq_heatBath_constant {τ : Ω → Ω} (hτ : Function.Involutive τ)
    {r : ℝ} (hr : 0 < r) :
    fairSwap τ = heatBath (fun _ : Ω => r) (orbitLabel τ) := by
  funext x y
  simp only [heatBath, orbitLabel_eq_iff hτ x y, fiberMass_constant hτ]
  by_cases hx : τ x = x
  · by_cases hy : y = x <;> simp [fairSwap, hx, hy, ne_of_gt hr]
  · by_cases hy : y = x
    · subst y
      simp only [fairSwap, hx, Ne.symm hx, if_false, if_true, add_zero,
        true_or]
      field_simp [ne_of_gt (add_pos hr hr)]
      ring
    · by_cases hy' : y = τ x
      · simp only [fairSwap, hx, hy', if_false, if_true, zero_add,
          or_true]
        field_simp [ne_of_gt (add_pos hr hr)]
        ring
      · simp [fairSwap, hy, hy']

/-- Normalized counting on the complete finite carrier. -/
def countingReference (_ : Ω) : ℝ := (Fintype.card Ω : ℝ)⁻¹

omit [DecidableEq Ω] in
theorem countingReference_pos [Nonempty Ω] (x : Ω) :
    0 < countingReference x := by
  exact inv_pos.mpr (Nat.cast_pos.mpr Fintype.card_pos)

omit [DecidableEq Ω] in
theorem countingReference_normalized [Nonempty Ω] :
    ∑ x : Ω, countingReference x = 1 := by
  simp [countingReference, Fintype.card_ne_zero]

/-- The fair exchange kernel is conditional resampling from the
normalized counting reference on its orbit fibers. -/
theorem fairSwap_eq_heatBath [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) :
    fairSwap τ = heatBath (countingReference (Ω := Ω)) (orbitLabel τ) :=
  fairSwap_eq_heatBath_constant hτ
    (inv_pos.mpr (Nat.cast_pos.mpr Fintype.card_pos))

/-- The constructed counting reference and orbit datum instantiate the
repair-law input of the existing four-law surface. No thermal calibration
or refinement field is supplied. -/
def countingRepairLawData [Nonempty Ω] (τ : Ω → Ω) : RepairLawData Ω (Finset Ω) where
  ref := countingReference
  ref_pos := countingReference_pos
  ref_law := countingReference_normalized
  visible := orbitLabel τ

/-- The kernel of the constructed typed repair data is the fair exchange
population kernel. -/
theorem countingRepairLaw_kernel_eq [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) :
    FourLawSurface.repairKernel (countingRepairLawData τ) = fairSwap τ :=
  (fairSwap_eq_heatBath hτ).symm

theorem fairSwap_row_sum [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) (x : Ω) : ∑ y, fairSwap τ x y = 1 := by
  rw [fairSwap_eq_heatBath hτ]
  exact heatBath_row_sum countingReference_pos x

theorem fairSwap_stationary [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) (y : Ω) :
    push countingReference (fairSwap τ) y = countingReference y := by
  rw [fairSwap_eq_heatBath hτ]
  exact heatBath_stationary countingReference_pos y

theorem fairSwap_detailedBalance [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) (x y : Ω) :
    countingReference x * fairSwap τ x y = countingReference y * fairSwap τ y x := by
  rw [fairSwap_eq_heatBath hτ]
  exact heatBath_detailedBalance x y

theorem fairSwap_idempotent [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) (x y : Ω) :
    ∑ z, fairSwap τ x z * fairSwap τ z y = fairSwap τ x y := by
  rw [fairSwap_eq_heatBath hτ]
  exact heatBath_idempotent countingReference_pos x y

/-- The native population rule inherits the transition-side second law
after its hold/exchange identity is supplied. -/
theorem fairSwap_secondLaw [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) (p : Ω → ℝ) (hp : ∀ x, 0 ≤ p x) :
    kl (push p (fairSwap τ)) countingReference ≤ kl p countingReference := by
  simpa only [countingRepairLaw_kernel_eq hτ] using
    FourLawSurface.second_contraction (countingRepairLawData τ) p hp

/-- The exact dissipation theorem of the four-law surface applies to the
constructed repair: its mean entropy production is the relative entropy
from the input to its repaired image. Strict positivity is explicit. -/
theorem fairSwap_exact_dissipation [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) (p : Ω → ℝ) (hp : ∀ x, 0 < p x) :
    ∑ x, ∑ y, p x * fairSwap τ x y *
      sigmaEP countingReference p (push p (fairSwap τ)) x y =
        kl p (push p (fairSwap τ)) := by
  simpa only [countingRepairLaw_kernel_eq hτ] using
    FourLawSurface.arrow_mean_entropy_production_eq_kl_to_repaired
      (countingRepairLawData τ) p hp

/-- The constructed repair satisfies the integral fluctuation identity
on strictly positive normalized inputs. -/
theorem fairSwap_integral_fluctuation [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) (p : Ω → ℝ) (hp : ∀ x, 0 < p x)
    (hp1 : ∑ x, p x = 1) :
    ∑ x, ∑ y, p x * fairSwap τ x y *
      Real.exp (-(sigmaEP countingReference p (push p (fairSwap τ)) x y)) = 1 := by
  simpa only [countingRepairLaw_kernel_eq hτ] using
    FourLawSurface.second_integral_fluctuation (countingRepairLawData τ) p hp hp1

/-- A positive input changes under this repair exactly when its mean
entropy production is strictly positive. -/
theorem fairSwap_strict_entropy_production [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) (p : Ω → ℝ) (hp : ∀ x, 0 < p x) :
    push p (fairSwap τ) ≠ p ↔
      0 < ∑ x, ∑ y, p x * fairSwap τ x y *
        sigmaEP countingReference p (push p (fairSwap τ)) x y := by
  simpa only [countingRepairLaw_kernel_eq hτ] using
    FourLawSurface.arrow_strict (countingRepairLawData τ) p hp

/-- Shannon entropy increases under fair involution repair. The proof
uses the same conditional-resampling law with unnormalized counting. -/
theorem fairSwap_shannon_mono {τ : Ω → Ω} (hτ : Function.Involutive τ)
    (p : Ω → ℝ) (hp : ∀ x, 0 ≤ p x) :
    shannon p ≤ shannon (push p (fairSwap τ)) := by
  have h := heatBath_secondLaw (π := fun _ : Ω => (1 : ℝ))
    (b := orbitLabel τ) (fun _ => by norm_num) p hp
  rw [← fairSwap_eq_heatBath_constant hτ (by norm_num : (0 : ℝ) < 1)] at h
  simpa only [shannon, kl] using neg_le_neg h

/-- Each repair row minimizes relative entropy among laws supported on
its orbit. Orbit masses of a global state are not selected by this result. -/
theorem fairSwap_row_optimal [Nonempty Ω] {τ : Ω → Ω}
    (hτ : Function.Involutive τ) (x : Ω) (p : Ω → ℝ)
    (hp : ∀ y, 0 ≤ p y) (hpsum : ∑ y, p y = 1)
    (hsupp : ∀ y, y ≠ x → y ≠ τ x → p y = 0) :
    kl (fairSwap τ x) countingReference ≤ kl p countingReference := by
  rw [fairSwap_eq_heatBath hτ]
  apply heatBath_row_optimal countingReference_pos x p hp hpsum
  intro y hy
  have hnot : ¬ (y = x ∨ y = τ x) := by
    exact fun h => hy ((orbitLabel_eq_iff hτ x y).mpr h)
  exact hsupp y (fun h => hnot (Or.inl h)) (fun h => hnot (Or.inr h))

omit [DecidableEq Ω] in
/-- The canonical counting reference gives a constant modular energy.
This identity contains no physical temperature calibration. -/
theorem countingReference_modularEnergy (x : Ω) :
    -Real.log (countingReference x) = Real.log (Fintype.card Ω : ℝ) := by
  simp [countingReference, Real.log_inv]

/-- Every calibration with a positive inverse-temperature multiplier compatible with the counting
reference has constant energy. Thus this branch cannot instantiate the
nondegenerate-thermometer premise by choosing a temperature convention. -/
theorem counting_calibration_energy_constant [Nonempty Ω] (τ : Ω → Ω)
    (C : CalibrationData (countingRepairLawData τ)) (x y : Ω) :
    C.energy x = C.energy y := by
  have h : C.beta * C.energy x + C.logZ = C.beta * C.energy y + C.logZ :=
    (C.thermal x).symm.trans (C.thermal y)
  exact mul_left_cancel₀ (ne_of_gt C.beta_pos) (add_right_cancel h)

/-- Relative entropy to normalized counting is the entropy deficit from
the logarithm of the carrier size. -/
theorem kl_counting_eq_log_card_sub_shannon [Nonempty Ω]
    (p : Ω → ℝ) (hpsum : ∑ x, p x = 1) :
    kl p countingReference = Real.log (Fintype.card Ω : ℝ) - shannon p := by
  rw [kl_eq_energy_sub_shannon p countingReference countingReference_pos]
  simp only [countingReference_modularEnergy, ← Finset.sum_mul, hpsum, one_mul]

/-- Any state-independent probability mixture of the exchange repairs
contracts relative entropy to the same normalized counting reference. The
mixture itself need not be idempotent or a one-step conditional expectation. -/
theorem randomScan_secondLaw [Nonempty Ω] {C : Type*} [Fintype C]
    (τ : C → Ω → Ω) (hτ : ∀ c, Function.Involutive (τ c))
    (a : C → ℝ) (ha0 : ∀ c, 0 ≤ a c) (ha1 : ∑ c, a c = 1)
    (p : Ω → ℝ) (hp : ∀ x, 0 ≤ p x) :
    kl (push p (fun x y => ∑ c, a c * fairSwap (τ c) x y)) countingReference ≤
      kl p countingReference := by
  have hstoch := mixture_stochastic a (fun c => fairSwap (τ c)) ha0 ha1
    (fun c => fairSwap_nonneg (τ c)) (fun c => fairSwap_row_sum (hτ c))
  have hstat : push countingReference
      (fun x y => ∑ c, a c * fairSwap (τ c) x y) = countingReference := by
    funext y
    exact mixture_stationary countingReference a (fun c => fairSwap (τ c)) ha1
      (fun c => fairSwap_stationary (hτ c)) y
  have h := kl_push_le p countingReference
    (fun x y => ∑ c, a c * fairSwap (τ c) x y) hp
    (fun x => (countingReference_pos x).le)
    (fun x hx => absurd hx (ne_of_gt (countingReference_pos x))) hstoch.1 hstoch.2
  rwa [hstat] at h

/-- A normalized native population law also has nondecreasing Shannon
entropy under a state-independent random scan of the exchange repairs. -/
theorem randomScan_shannon_mono [Nonempty Ω] {C : Type*} [Fintype C]
    (τ : C → Ω → Ω) (hτ : ∀ c, Function.Involutive (τ c))
    (a : C → ℝ) (ha0 : ∀ c, 0 ≤ a c) (ha1 : ∑ c, a c = 1)
    (p : Ω → ℝ) (hp : ∀ x, 0 ≤ p x) (hpsum : ∑ x, p x = 1) :
    shannon p ≤ shannon (push p (fun x y => ∑ c, a c * fairSwap (τ c) x y)) := by
  have hmass : ∑ y, push p (fun x y => ∑ c, a c * fairSwap (τ c) x y) y = 1 := by
    rw [push_total p _ (mixture_row_sum a (fun c => fairSwap (τ c)) ha1
      (fun c => fairSwap_row_sum (hτ c))), hpsum]
  have h := randomScan_secondLaw τ hτ a ha0 ha1 p hp
  rw [kl_counting_eq_log_card_sub_shannon _ hmass,
    kl_counting_eq_log_card_sub_shannon p hpsum] at h
  linarith

end

#print axioms orbitLabel_eq_iff
#print axioms fairSwap_eq_heatBath_constant
#print axioms countingReference_normalized
#print axioms fairSwap_eq_heatBath
#print axioms countingRepairLaw_kernel_eq
#print axioms fairSwap_stationary
#print axioms fairSwap_detailedBalance
#print axioms fairSwap_idempotent
#print axioms fairSwap_secondLaw
#print axioms fairSwap_exact_dissipation
#print axioms fairSwap_integral_fluctuation
#print axioms fairSwap_strict_entropy_production
#print axioms fairSwap_shannon_mono
#print axioms fairSwap_row_optimal
#print axioms countingReference_modularEnergy
#print axioms counting_calibration_energy_constant
#print axioms kl_counting_eq_log_card_sub_shannon
#print axioms randomScan_secondLaw
#print axioms randomScan_shannon_mono

end OPH.Thermodynamics.NativeBinary
