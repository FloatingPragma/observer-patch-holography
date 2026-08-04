import Mathlib

open scoped BigOperators

namespace OPH.Dynamics

/-!
# Exact dual-channel criterion for a protected finite charge

This module states the reusable linear criterion behind finite charge
conservation. Positivity and trace preservation are separate channel receipts.
No observable is identified with a physical charge here.
-/

variable {State : Type*} [Fintype State]

/-- Pairing of a finite signed state with a real-valued charge observable. -/
noncomputable def chargeExpectation
    (p charge : State → ℝ) : ℝ :=
  ∑ x, p x * charge x

/-- Schrödinger and Heisenberg linear maps are paired as duals when their
finite expectation pairings agree on every signed state and observable. -/
def AreExpectationDual
    (channel dualChannel : (State → ℝ) →ₗ[ℝ] (State → ℝ)) : Prop :=
  ∀ p observable,
    chargeExpectation (channel p) observable =
      chargeExpectation p (dualChannel observable)

/-- A charge expectation is protected on every finite signed state exactly
when the dual channel fixes the charge observable. -/
theorem chargeExpectation_preserved_iff_dual_fixed
    [DecidableEq State]
    (channel dualChannel : (State → ℝ) →ₗ[ℝ] (State → ℝ))
    (hdual : AreExpectationDual channel dualChannel) (charge : State → ℝ) :
    (∀ p, chargeExpectation (channel p) charge = chargeExpectation p charge) ↔
      dualChannel charge = charge := by
  constructor
  · intro hinvariant
    funext x
    let delta : State → ℝ := fun y ↦ if y = x then 1 else 0
    have hpair : chargeExpectation delta (dualChannel charge) =
        chargeExpectation delta charge := by
      calc
        chargeExpectation delta (dualChannel charge) =
            chargeExpectation (channel delta) charge := (hdual delta charge).symm
        _ = chargeExpectation delta charge := hinvariant delta
    simpa [chargeExpectation, delta] using hpair
  · intro hfixed p
    rw [hdual p charge, hfixed]

/-- The forward protected-charge condition in the form consumed by finite
channel constructions. -/
theorem chargeExpectation_preserved_of_dual_fixed
    (channel dualChannel : (State → ℝ) →ₗ[ℝ] (State → ℝ))
    (hdual : AreExpectationDual channel dualChannel) (charge : State → ℝ)
    (hfixed : dualChannel charge = charge) (p : State → ℝ) :
    chargeExpectation (channel p) charge = chargeExpectation p charge := by
  rw [hdual p charge, hfixed]

/-! ## Finite-kernel specialization -/

/-- Schrödinger pushforward through an arbitrary finite real kernel. -/
noncomputable def kernelPush (kernel : State → State → ℝ) :
    (State → ℝ) →ₗ[ℝ] (State → ℝ) where
  toFun p y := ∑ x, p x * kernel x y
  map_add' p q := by
    funext y
    simp only [Pi.add_apply]
    rw [← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro x _
    ring
  map_smul' a p := by
    funext y
    simp only [Pi.smul_apply, RingHom.id_apply, smul_eq_mul]
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro x _
    ring

/-- Heisenberg pullback through the same finite kernel. -/
noncomputable def kernelPull (kernel : State → State → ℝ) :
    (State → ℝ) →ₗ[ℝ] (State → ℝ) where
  toFun observable x := ∑ y, kernel x y * observable y
  map_add' observable other := by
    funext x
    simp only [Pi.add_apply]
    rw [← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro y _
    ring
  map_smul' a observable := by
    funext x
    simp only [Pi.smul_apply, RingHom.id_apply, smul_eq_mul]
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro y _
    ring

/-- Finite kernel pushforward and pullback obey exact push-pull duality. -/
theorem kernelPush_kernelPull_areExpectationDual
    (kernel : State → State → ℝ) :
    AreExpectationDual (kernelPush kernel) (kernelPull kernel) := by
  intro p observable
  change (∑ y, (∑ x, p x * kernel x y) * observable y) =
    ∑ x, p x * ∑ y, kernel x y * observable y
  calc
    (∑ y, (∑ x, p x * kernel x y) * observable y) =
        ∑ y, ∑ x, p x * kernel x y * observable y := by
      apply Finset.sum_congr rfl
      intro y _
      rw [Finset.sum_mul]
    _ = ∑ x, ∑ y, p x * kernel x y * observable y := by
      rw [Finset.sum_comm]
    _ = ∑ x, p x * ∑ y, kernel x y * observable y := by
      apply Finset.sum_congr rfl
      intro x _
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro y _
      ring

/-- For any finite kernel, all-state charge conservation is equivalent to
the pointwise fixed-observable equation for its pullback. -/
theorem kernel_chargeExpectation_preserved_iff_pull_fixed
    [DecidableEq State] (kernel : State → State → ℝ) (charge : State → ℝ) :
    (∀ p, chargeExpectation (kernelPush kernel p) charge =
      chargeExpectation p charge) ↔ kernelPull kernel charge = charge :=
  chargeExpectation_preserved_iff_dual_fixed
    (kernelPush kernel) (kernelPull kernel)
    (kernelPush_kernelPull_areExpectationDual kernel) charge

/-! ## Symmetry is not conservation: an exact negative control -/

/-- A channel is covariant under a declared linear action when the two maps
commute. This is a symmetry condition on the channel, not a protected-charge
condition. -/
def ChannelCovariantUnder
    (channel action : (State → ℝ) →ₗ[ℝ] (State → ℝ)) : Prop :=
  channel.comp action = action.comp channel

/-- Swap the two coordinates of a real function on `Fin 2`. -/
noncomputable def twoStateSwapAction :
    (Fin 2 → ℝ) →ₗ[ℝ] (Fin 2 → ℝ) where
  toFun f i := if i = 0 then f 1 else f 0
  map_add' f g := by
    funext i
    by_cases h : i = 0 <;> simp [h]
  map_smul' a f := by
    funext i
    by_cases h : i = 0 <;> simp [h]

/-- The uniform two-state averaging channel. It erases the odd component of
an arbitrary signed state. -/
noncomputable def twoStateAverageChannel :
    (Fin 2 → ℝ) →ₗ[ℝ] (Fin 2 → ℝ) where
  toFun p _ := (p 0 + p 1) / 2
  map_add' p q := by
    funext i
    simp only [Pi.add_apply]
    ring
  map_smul' a p := by
    funext i
    simp only [Pi.smul_apply, RingHom.id_apply, smul_eq_mul]
    ring

/-- The nonzero odd charge under the two-state swap. -/
def twoStateOddCharge : Fin 2 → ℝ := fun i ↦
  if i = 0 then 1 else -1

theorem twoStateOddCharge_ne_zero : twoStateOddCharge ≠ 0 := by
  intro h
  have h0 := congrFun h 0
  norm_num [twoStateOddCharge] at h0

/-- Positive control: the identity channel protects a nonzero charge. -/
theorem identity_protects_twoStateOddCharge (p : Fin 2 → ℝ) :
    chargeExpectation ((LinearMap.id (R := ℝ) (M := Fin 2 → ℝ)) p)
        twoStateOddCharge =
      chargeExpectation p twoStateOddCharge := by
  rfl

/-- The uniform averaging channel is exactly covariant under the nontrivial
two-state swap. -/
theorem twoStateAverage_covariant_under_swap :
    ChannelCovariantUnder twoStateAverageChannel twoStateSwapAction := by
  apply LinearMap.ext
  intro p
  funext i
  fin_cases i <;>
    simp [twoStateAverageChannel, twoStateSwapAction] <;> ring

/-- Exact negative control: channel covariance under a nontrivial action does
not imply conservation of a charge that transforms under that action. The
covariant averaging channel destroys the nonzero odd charge. -/
theorem channel_covariance_does_not_imply_charge_conservation :
    ChannelCovariantUnder twoStateAverageChannel twoStateSwapAction ∧
      ¬(∀ p : Fin 2 → ℝ,
        chargeExpectation (twoStateAverageChannel p) twoStateOddCharge =
          chargeExpectation p twoStateOddCharge) := by
  refine ⟨twoStateAverage_covariant_under_swap, ?_⟩
  intro h
  let delta : Fin 2 → ℝ := fun i ↦ if i = 0 then 1 else 0
  have hdelta := h delta
  norm_num [chargeExpectation, twoStateAverageChannel, twoStateOddCharge,
    delta, Fin.sum_univ_two] at hdelta

/-! ## Axiom audit -/

#print axioms chargeExpectation_preserved_iff_dual_fixed
#print axioms chargeExpectation_preserved_of_dual_fixed
#print axioms kernelPush_kernelPull_areExpectationDual
#print axioms kernel_chargeExpectation_preserved_iff_pull_fixed
#print axioms twoStateOddCharge_ne_zero
#print axioms identity_protects_twoStateOddCharge
#print axioms twoStateAverage_covariant_under_swap
#print axioms channel_covariance_does_not_imply_charge_conservation

end OPH.Dynamics
