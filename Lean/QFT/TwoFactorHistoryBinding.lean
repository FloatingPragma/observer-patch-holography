import QFT.HistoryActionBinding

/-!
# Two-factor Gibbs/history binding and relative-coupling identifiability

This module extends `HistoryActionBinding` from one certified factor to a
product of two finite additive step groups.  It proves three distinct facts:

1. the Gibbs normalizer and kernel for a sum of two factor costs factorize;
2. the corresponding product history law is bound to that constructed cost;
3. if each factor cost has one nonconstant direction, the two weighted
   coefficients are identifiable modulo exactly one common action-multiplier
   rescaling.  A single multiplier cannot absorb a change in their ratio.

Concrete definitions then transcribe both P blocks and the F-family
`su(3) + so(3)` blocks from `code/e9_kinetic/kinetic_form_v1.json`.  This earns
full-P and one color-bearing F constructed-kernel bindings on finite step
groups of cardinalities `3^6` and `3^11`.

BOUNDARY: every kernel below is constructed from the displayed kinetic cost.
Nothing identifies it with an independently supplied OPH source transition
law, selects the two coefficient values, fixes their ratio, supplies physical
units, or takes a continuum limit. Identifiability says equality inside this
displayed constructed two-factor Gibbs family fixes the relative ratio; it
does not predict or source-select that ratio.
-/

namespace OPH.QFT

open OPH.InformationProjection

/-! ## General two-factor theorem -/

/-- A two-factor cost with independent coefficients. -/
def twoFactorCost {A B : Type*} (qA : A → ℝ) (qB : B → ℝ)
    (a b : ℝ) (z : A × B) : ℝ :=
  a * qA z.1 + b * qB z.2

/-- The product normalizer factorizes exactly. -/
theorem gibbsNorm_twoFactor {A B : Type*} [Fintype A] [Fintype B]
    (qA : A → ℝ) (qB : B → ℝ) (a b : ℝ) :
    gibbsNorm (twoFactorCost qA qB a b) =
      gibbsNorm (fun x ↦ a * qA x) * gibbsNorm (fun y ↦ b * qB y) := by
  unfold gibbsNorm twoFactorCost
  rw [Fintype.sum_prod_type]
  simp_rw [neg_add, Real.exp_add]
  calc
    ∑ x, ∑ y, Real.exp (-(a * qA x)) * Real.exp (-(b * qB y)) =
        ∑ x, Real.exp (-(a * qA x)) * ∑ y, Real.exp (-(b * qB y)) := by
          apply Finset.sum_congr rfl
          intro x _
          rw [Finset.mul_sum]
    _ = (∑ x, Real.exp (-(a * qA x))) * ∑ y, Real.exp (-(b * qB y)) := by
          rw [Finset.sum_mul]

/-- The constructed two-factor Gibbs kernel is the product kernel. -/
theorem gibbsKernel_twoFactor {A B : Type*}
    [AddCommGroup A] [AddCommGroup B] [Fintype A] [Fintype B]
    [Nonempty A] [Nonempty B] (qA : A → ℝ) (qB : B → ℝ)
    (a b : ℝ) (x y : A × B) :
    gibbsKernel (twoFactorCost qA qB a b) x y =
      gibbsKernel (fun u ↦ a * qA u) x.1 y.1 *
        gibbsKernel (fun v ↦ b * qB v) x.2 y.2 := by
  unfold gibbsKernel
  rw [gibbsNorm_twoFactor]
  simp only [twoFactorCost, Prod.fst_sub, Prod.snd_sub]
  change Real.exp (-(a * qA (y.1 - x.1) + b * qB (y.2 - x.2))) /
      (gibbsNorm (fun u ↦ a * qA u) * gibbsNorm (fun v ↦ b * qB v)) = _
  rw [show -(a * qA (y.1 - x.1) + b * qB (y.2 - x.2)) =
      -(a * qA (y.1 - x.1)) + -(b * qB (y.2 - x.2)) by ring,
    Real.exp_add, div_mul_div_comm]

/-- Product-cost history action is the sum of the projected factor actions. -/
theorem incrementAction_twoFactor {A B : Type*}
    [AddCommGroup A] [AddCommGroup B] (qA : A → ℝ) (qB : B → ℝ)
    (a b : ℝ) (n : ℕ) (path : PathSpace (A × B) n) :
    incrementAction (twoFactorCost qA qB a b) n path =
      a * incrementAction qA n (fun i ↦ (path i).1) +
        b * incrementAction qB n (fun i ↦ (path i).2) := by
  unfold incrementAction twoFactorCost
  simp only [Prod.fst_sub, Prod.snd_sub]
  rw [Finset.sum_add_distrib, Finset.mul_sum, Finset.mul_sum]

/-- The general binding theorem instantiated on a two-factor constructed cost. -/
theorem twoFactor_action_reproduces_law {A B : Type*}
    [AddCommGroup A] [AddCommGroup B] [Fintype A] [Fintype B]
    [Nonempty A] [Nonempty B] (pi : A × B → ℝ) (qA : A → ℝ)
    (qB : B → ℝ) (a b : ℝ) (n : ℕ) (hpi : ∀ x, 0 < pi x)
    (hpi1 : ∑ x, pi x = 1) :
    tilt (stepUniformRef pi n) (incrementAction (twoFactorCost qA qB a b) n) 1 =
      markovPathLaw pi (gibbsKernel (twoFactorCost qA qB a b)) n :=
  increment_action_reproduces_law pi (twoFactorCost qA qB a b) n hpi hpi1

/--
Independent variation of the two factor costs identifies the two
multiplier-weighted coefficients, even when the action equality is allowed one
additive gauge constant.
-/
theorem twoFactor_weighted_coefficients_identifiable {A B : Type*}
    (qA : A → ℝ) (qB : B → ℝ) (x₀ x₁ : A) (y₀ y₁ : B)
    (hx : qA x₁ ≠ qA x₀) (hy : qB y₁ ≠ qB y₀)
    (a b a' b' lam lam' gauge : ℝ)
    (h : ∀ z : A × B,
      lam * twoFactorCost qA qB a b z =
        lam' * twoFactorCost qA qB a' b' z + gauge) :
    lam * a = lam' * a' ∧ lam * b = lam' * b' := by
  have hA₁ := h (x₁, y₀)
  have hA₀ := h (x₀, y₀)
  have hB₁ := h (x₀, y₁)
  have hB₀ := h (x₀, y₀)
  have factorA : (lam * a - lam' * a') * (qA x₁ - qA x₀) = 0 := by
    unfold twoFactorCost at hA₁ hA₀
    linear_combination hA₁ - hA₀
  have factorB : (lam * b - lam' * b') * (qB y₁ - qB y₀) = 0 := by
    unfold twoFactorCost at hB₁ hB₀
    linear_combination hB₁ - hB₀
  constructor
  · rcases mul_eq_zero.mp factorA with hcoef | hcost
    · linarith
    · exact False.elim (hx (sub_eq_zero.mp hcost))
  · rcases mul_eq_zero.mp factorB with hcoef | hcost
    · linarith
    · exact False.elim (hy (sub_eq_zero.mp hcost))

/--
Equality of the two *constructed transition kernels* supplies exactly the
additive cost gauge required above, hence identifies both
multiplier-weighted coefficients.  This is the law-level theorem; it still
does not assert that either constructed kernel is the OPH source kernel.
-/
theorem twoFactor_kernel_identifies_weighted_coefficients {A B : Type*}
    [AddCommGroup A] [AddCommGroup B] [Fintype A] [Fintype B]
    [Nonempty A] [Nonempty B]
    (qA : A → ℝ) (qB : B → ℝ) (x₀ x₁ : A) (y₀ y₁ : B)
    (hx : qA x₁ ≠ qA x₀) (hy : qB y₁ ≠ qB y₀)
    (a b a' b' lam lam' : ℝ)
    (hk : ∀ x y : A × B,
      gibbsKernel (fun z ↦ lam * twoFactorCost qA qB a b z) x y =
        gibbsKernel (fun z ↦ lam' * twoFactorCost qA qB a' b' z) x y) :
    lam * a = lam' * a' ∧ lam * b = lam' * b' := by
  let leftCost : A × B → ℝ := fun z ↦ lam * twoFactorCost qA qB a b z
  let rightCost : A × B → ℝ := fun z ↦ lam' * twoFactorCost qA qB a' b' z
  apply twoFactor_weighted_coefficients_identifiable qA qB x₀ x₁ y₀ y₁
    hx hy a b a' b' lam lam'
    (Real.log (gibbsNorm rightCost) - Real.log (gibbsNorm leftCost))
  intro z
  have hlog := congrArg Real.log (hk 0 z)
  rw [gibbsKernel_log leftCost 0 z, gibbsKernel_log rightCost 0 z] at hlog
  simp only [sub_zero] at hlog
  dsimp [leftCost, rightCost] at hlog ⊢
  linarith

/--
With nonzero multipliers, gauge-equivalent two-factor actions differ only by
one common coefficient rescaling.  There is no independent per-factor
multiplier gauge.
-/
theorem twoFactor_only_common_multiplier_scaling {A B : Type*}
    (qA : A → ℝ) (qB : B → ℝ) (x₀ x₁ : A) (y₀ y₁ : B)
    (hx : qA x₁ ≠ qA x₀) (hy : qB y₁ ≠ qB y₀)
    (a b a' b' lam lam' gauge : ℝ) (hlam : lam ≠ 0) (hlam' : lam' ≠ 0)
    (h : ∀ z : A × B,
      lam * twoFactorCost qA qB a b z =
        lam' * twoFactorCost qA qB a' b' z + gauge) :
    ∃ scale : ℝ, scale ≠ 0 ∧ a' = scale * a ∧ b' = scale * b := by
  have hc := twoFactor_weighted_coefficients_identifiable qA qB x₀ x₁ y₀ y₁
    hx hy a b a' b' lam lam' gauge h
  refine ⟨lam / lam', div_ne_zero hlam hlam', ?_, ?_⟩
  · calc
      a' = (lam * a) / lam' := (eq_div_iff hlam').2 (by nlinarith [hc.1])
      _ = (lam / lam') * a := by ring
  · calc
      b' = (lam * b) / lam' := (eq_div_iff hlam').2 (by nlinarith [hc.2])
      _ = (lam / lam') * b := by ring

/-- Kernel equality leaves exactly the same one-dimensional common scaling orbit. -/
theorem twoFactor_kernel_only_common_multiplier_scaling {A B : Type*}
    [AddCommGroup A] [AddCommGroup B] [Fintype A] [Fintype B]
    [Nonempty A] [Nonempty B]
    (qA : A → ℝ) (qB : B → ℝ) (x₀ x₁ : A) (y₀ y₁ : B)
    (hx : qA x₁ ≠ qA x₀) (hy : qB y₁ ≠ qB y₀)
    (a b a' b' lam lam' : ℝ) (hlam : lam ≠ 0) (hlam' : lam' ≠ 0)
    (hk : ∀ x y : A × B,
      gibbsKernel (fun z ↦ lam * twoFactorCost qA qB a b z) x y =
        gibbsKernel (fun z ↦ lam' * twoFactorCost qA qB a' b' z) x y) :
    ∃ scale : ℝ, scale ≠ 0 ∧ a' = scale * a ∧ b' = scale * b := by
  have hc := twoFactor_kernel_identifies_weighted_coefficients qA qB x₀ x₁ y₀ y₁
    hx hy a b a' b' lam lam' hk
  refine ⟨lam / lam', div_ne_zero hlam hlam', ?_, ?_⟩
  · calc
      a' = (lam * a) / lam' := (eq_div_iff hlam').2 (by nlinarith [hc.1])
      _ = (lam / lam') * a := by ring
  · calc
      b' = (lam * b) / lam' := (eq_div_iff hlam').2 (by nlinarith [hc.2])
      _ = (lam / lam') * b := by ring

/-- Consequently the relative ratio is invariant under the only gauge left. -/
theorem twoFactor_relative_ratio_identifiable {A B : Type*}
    (qA : A → ℝ) (qB : B → ℝ) (x₀ x₁ : A) (y₀ y₁ : B)
    (hx : qA x₁ ≠ qA x₀) (hy : qB y₁ ≠ qB y₀)
    (a b a' b' lam lam' gauge : ℝ) (hlam : lam ≠ 0) (hlam' : lam' ≠ 0)
    (hb : b ≠ 0)
    (h : ∀ z : A × B,
      lam * twoFactorCost qA qB a b z =
        lam' * twoFactorCost qA qB a' b' z + gauge) :
    a / b = a' / b' := by
  rcases twoFactor_only_common_multiplier_scaling qA qB x₀ x₁ y₀ y₁ hx hy
      a b a' b' lam lam' gauge hlam hlam' h with ⟨scale, hscale, ha, hbscale⟩
  rw [ha, hbscale]
  field_simp

/-! ## Full P: both certified A1 blocks -/

/-- Galois-conjugate P block on the `three_minus` factor. -/
noncomputable def mirrorSectorKinetic (i j : Fin 3) : ℝ :=
  (sectorKineticRat i j : ℝ) - (sectorKineticIrr i j : ℝ) * Real.sqrt 5

noncomputable def mirrorSectorIncrement (v : SectorStep) : ℝ :=
  ∑ i, ∑ j,
    (sectorLift v i : ℝ) * mirrorSectorKinetic i j * (sectorLift v j : ℝ)

theorem mirrorSectorIncrement_eq_pair (v : SectorStep) :
    mirrorSectorIncrement v =
      (sectorIncrementRat v : ℝ) - (sectorIncrementIrr v : ℝ) * Real.sqrt 5 := by
  unfold mirrorSectorIncrement sectorIncrementRat sectorIncrementIrr mirrorSectorKinetic
  push_cast
  rw [Finset.sum_mul, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.sum_mul, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun j _ => ?_
  ring

theorem mirrorSectorIncrement_zero : mirrorSectorIncrement 0 = 0 := by
  rw [mirrorSectorIncrement_eq_pair,
    show sectorIncrementRat 0 = 0 from by decide,
    show sectorIncrementIrr 0 = 0 from by decide]
  norm_num

theorem mirrorSectorIncrement_unit :
    mirrorSectorIncrement (1, 0, 0) = 1200 + 400 * Real.sqrt 5 := by
  rw [mirrorSectorIncrement_eq_pair,
    show sectorIncrementRat (1, 0, 0) = 1200 from by decide,
    show sectorIncrementIrr (1, 0, 0) = -400 from by decide]
  push_cast
  ring

theorem mirrorSectorIncrement_unit_pos : 0 < mirrorSectorIncrement (1, 0, 0) := by
  rw [mirrorSectorIncrement_unit]
  positivity

abbrev FullPStep := SectorStep × SectorStep

theorem fullPStep_card : Fintype.card FullPStep = 729 := by
  simp only [FullPStep, SectorStep, Fintype.card_prod, ZMod.card]

/-- Both P factors, with independently visible coefficients. -/
noncomputable def fullPCost (a b : ℝ) : FullPStep → ℝ :=
  twoFactorCost sectorIncrement mirrorSectorIncrement a b

/-- Constructed-kernel history binding for the complete six-dimensional P algebra. -/
theorem fullP_action_reproduces_law (pi : FullPStep → ℝ) (a b : ℝ) (n : ℕ)
    (hpi : ∀ x, 0 < pi x) (hpi1 : ∑ x, pi x = 1) :
    tilt (stepUniformRef pi n) (incrementAction (fullPCost a b) n) 1 =
      markovPathLaw pi (gibbsKernel (fullPCost a b)) n :=
  twoFactor_action_reproduces_law pi sectorIncrement mirrorSectorIncrement a b n hpi hpi1

/-- Full-P relative weights are identifiable modulo only common scaling. -/
theorem fullP_relative_coefficients_identifiable
    (a b a' b' lam lam' gauge : ℝ) (hlam : lam ≠ 0) (hlam' : lam' ≠ 0)
    (h : ∀ z : FullPStep,
      lam * fullPCost a b z = lam' * fullPCost a' b' z + gauge) :
    ∃ scale : ℝ, scale ≠ 0 ∧ a' = scale * a ∧ b' = scale * b := by
  apply twoFactor_only_common_multiplier_scaling sectorIncrement mirrorSectorIncrement
    0 (1, 0, 0) 0 (1, 0, 0)
  · rw [sectorIncrement_zero]
    exact ne_of_gt (sectorIncrement_pos _ (by decide))
  · rw [mirrorSectorIncrement_zero]
    exact ne_of_gt mirrorSectorIncrement_unit_pos
  · exact hlam
  · exact hlam'
  · exact h

/-- Equality of full-P constructed kernels leaves only the common scaling orbit. -/
theorem fullP_kernel_relative_coefficients_identifiable
    (a b a' b' lam lam' : ℝ) (hlam : lam ≠ 0) (hlam' : lam' ≠ 0)
    (hk : ∀ x y : FullPStep,
      gibbsKernel (fun z ↦ lam * fullPCost a b z) x y =
        gibbsKernel (fun z ↦ lam' * fullPCost a' b' z) x y) :
    ∃ scale : ℝ, scale ≠ 0 ∧ a' = scale * a ∧ b' = scale * b := by
  apply twoFactor_kernel_only_common_multiplier_scaling
    sectorIncrement mirrorSectorIncrement 0 (1, 0, 0) 0 (1, 0, 0)
  · rw [sectorIncrement_zero]
    exact ne_of_gt (sectorIncrement_pos _ (by decide))
  · rw [mirrorSectorIncrement_zero]
    exact ne_of_gt mirrorSectorIncrement_unit_pos
  · exact hlam
  · exact hlam'
  · intro x y
    simpa [fullPCost] using hk x y

/-! ## F: certified A2 color block plus mirror A1 block -/

abbrev FColorStep := Fin 8 → ZMod 3

def fColorLift (v : FColorStep) (i : Fin 8) : ℤ := centeredLift (v i)

/-- Rational part of the certified eight-dimensional F color block. -/
def fColorKineticRat (i j : Fin 8) : ℤ :=
  if hi : i.val < 3 then
    if hj : j.val < 3 then
      6 * sectorKineticRat ⟨i.val, hi⟩ ⟨j.val, hj⟩
    else 0
  else if _hj : j.val < 3 then 0
  else if i = j then -480 else -240

/-- `sqrt(5)` part of the certified eight-dimensional F color block. -/
def fColorKineticIrr (i j : Fin 8) : ℤ :=
  if hi : i.val < 3 then
    if hj : j.val < 3 then
      6 * sectorKineticIrr ⟨i.val, hi⟩ ⟨j.val, hj⟩
    else 0
  else if _hj : j.val < 3 then 0
  else if i = j then 480 else 240

noncomputable def fColorKinetic (i j : Fin 8) : ℝ :=
  (fColorKineticRat i j : ℝ) + (fColorKineticIrr i j : ℝ) * Real.sqrt 5

noncomputable def fColorIncrement (v : FColorStep) : ℝ :=
  ∑ i, ∑ j, (fColorLift v i : ℝ) * fColorKinetic i j * (fColorLift v j : ℝ)

theorem fColorKinetic_00 : fColorKinetic 0 0 = 7200 - 2400 * Real.sqrt 5 := by
  simp [fColorKinetic, fColorKineticRat, fColorKineticIrr,
    sectorKineticRat, sectorKineticIrr]
  ring

theorem fColorKinetic_33 : fColorKinetic 3 3 = -480 + 480 * Real.sqrt 5 := by
  simp [fColorKinetic, fColorKineticRat, fColorKineticIrr]

theorem fColorIncrement_zero : fColorIncrement 0 = 0 := by
  classical
  norm_num [fColorIncrement, fColorLift, centeredLift]

def fColorUnit : FColorStep := fun i ↦ if i = 0 then 1 else 0

theorem fColorIncrement_unit :
    fColorIncrement fColorUnit = 7200 - 2400 * Real.sqrt 5 := by
  classical
  norm_num [fColorIncrement, fColorUnit, fColorLift, centeredLift]
  rw [fColorKinetic_00]

theorem fColorIncrement_unit_pos : 0 < fColorIncrement fColorUnit := by
  rw [fColorIncrement_unit]
  have hs : Real.sqrt 5 ^ 2 = 5 := by norm_num
  have hspos : 0 < Real.sqrt 5 := by positivity
  nlinarith

abbrev FFamilyStep := FColorStep × SectorStep

theorem fFamilyStep_card : Fintype.card FFamilyStep = 177147 := by
  simp only [FFamilyStep, FColorStep, SectorStep, Fintype.card_prod,
    Fintype.card_fun, Fintype.card_fin, ZMod.card]
  norm_num

/-- The two certified F factors, with independently visible coefficients. -/
noncomputable def fFamilyCost (a b : ℝ) : FFamilyStep → ℝ :=
  twoFactorCost fColorIncrement mirrorSectorIncrement a b

/-- Constructed-kernel binding for the full eleven-dimensional `su(3)+so(3)` F algebra. -/
theorem fFamily_action_reproduces_law (pi : FFamilyStep → ℝ) (a b : ℝ) (n : ℕ)
    (hpi : ∀ x, 0 < pi x) (hpi1 : ∑ x, pi x = 1) :
    tilt (stepUniformRef pi n) (incrementAction (fFamilyCost a b) n) 1 =
      markovPathLaw pi (gibbsKernel (fFamilyCost a b)) n :=
  twoFactor_action_reproduces_law pi fColorIncrement mirrorSectorIncrement a b n hpi hpi1

/-- F-family relative weights are identifiable modulo only common scaling. -/
theorem fFamily_relative_coefficients_identifiable
    (a b a' b' lam lam' gauge : ℝ) (hlam : lam ≠ 0) (hlam' : lam' ≠ 0)
    (h : ∀ z : FFamilyStep,
      lam * fFamilyCost a b z = lam' * fFamilyCost a' b' z + gauge) :
    ∃ scale : ℝ, scale ≠ 0 ∧ a' = scale * a ∧ b' = scale * b := by
  apply twoFactor_only_common_multiplier_scaling fColorIncrement mirrorSectorIncrement
    0 fColorUnit 0 (1, 0, 0)
  · rw [fColorIncrement_zero]
    exact ne_of_gt fColorIncrement_unit_pos
  · rw [mirrorSectorIncrement_zero]
    exact ne_of_gt mirrorSectorIncrement_unit_pos
  · exact hlam
  · exact hlam'
  · exact h

/-- Equality of F-family constructed kernels leaves only the common scaling orbit. -/
theorem fFamily_kernel_relative_coefficients_identifiable
    (a b a' b' lam lam' : ℝ) (hlam : lam ≠ 0) (hlam' : lam' ≠ 0)
    (hk : ∀ x y : FFamilyStep,
      gibbsKernel (fun z ↦ lam * fFamilyCost a b z) x y =
        gibbsKernel (fun z ↦ lam' * fFamilyCost a' b' z) x y) :
    ∃ scale : ℝ, scale ≠ 0 ∧ a' = scale * a ∧ b' = scale * b := by
  apply twoFactor_kernel_only_common_multiplier_scaling
    fColorIncrement mirrorSectorIncrement 0 fColorUnit 0 (1, 0, 0)
  · rw [fColorIncrement_zero]
    exact ne_of_gt fColorIncrement_unit_pos
  · rw [mirrorSectorIncrement_zero]
    exact ne_of_gt mirrorSectorIncrement_unit_pos
  · exact hlam
  · exact hlam'
  · intro x y
    simpa [fFamilyCost] using hk x y

/-! ## Axiom audit -/

#print axioms gibbsNorm_twoFactor
#print axioms gibbsKernel_twoFactor
#print axioms incrementAction_twoFactor
#print axioms twoFactor_weighted_coefficients_identifiable
#print axioms twoFactor_kernel_identifies_weighted_coefficients
#print axioms twoFactor_only_common_multiplier_scaling
#print axioms twoFactor_kernel_only_common_multiplier_scaling
#print axioms twoFactor_relative_ratio_identifiable
#print axioms fullP_action_reproduces_law
#print axioms fullP_relative_coefficients_identifiable
#print axioms fullP_kernel_relative_coefficients_identifiable
#print axioms fFamily_action_reproduces_law
#print axioms fFamily_relative_coefficients_identifiable
#print axioms fFamily_kernel_relative_coefficients_identifiable

end OPH.QFT
