import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Piecewise
import Mathlib.Tactic

/-!
# Finite cell energy conservation commutes with spatial coarse-graining

The actual same-action CAR currents and nearest-face routing are certified in
`model/coarse_energy_ward`. Finite sums prove that an integrated-cell
continuity equation survives arbitrary block maps and their composition.
The quantities are energy increments and integrated outgoing fluxes; no
point-current approximation, symmetric Lorentz tensor or Einstein premise
inhabitant is asserted. The source algebra is not encoded as a Lean value.
-/
namespace OPH.QFT.CoarseEnergyWard

variable {X Y Z M : Type*} [Fintype X] [Fintype Y] [Fintype Z]
  [DecidableEq Y] [DecidableEq Z] [AddCommGroup M]

def push (c : X → Y) (f : X → M) (a : Y) : M :=
  ∑ x, if c x = a then f x else 0

theorem push_total (c : X → Y) (f : X → M) :
    ∑ a, push c f a = ∑ x, f x := by
  simp only [push]
  rw [Finset.sum_comm]
  simp

omit [Fintype Y] in
theorem push_add (c : X → Y) (f g : X → M) (a : Y) :
    push c (fun x => f x + g x) a = push c f a + push c g a := by
  simp only [push, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro x _
  split_ifs <;> simp

theorem ite_sum (p : Prop) [Decidable p] (f : X → M) :
    (if p then ∑ x, f x else 0) = ∑ x, if p then f x else 0 := by
  by_cases h : p <;> simp [h]

omit [Fintype Z] in
theorem push_compose (c : X → Y) (d : Y → Z) (f : X → M) :
    push d (push c f) = push (d ∘ c) f := by
  funext a
  simp only [push, Function.comp_apply]
  simp_rw [ite_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro x _
  rw [Finset.sum_eq_single (c x)]
  · simp
  · intro y _ hy
    by_cases hd : d y = a <;> simp [hd, Ne.symm hy]
  · simp

def coarse (c : X → Y) (J : X → X → M) (a b : Y) : M :=
  push c (fun x => push c (J x) b) a

theorem divergence_natural (c : X → Y) (J : X → X → M) (a : Y) :
    ∑ b, coarse c J a b = push c (fun x => ∑ y, J x y) a := by
  simp only [coarse, push]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro x _
  by_cases h : c x = a
  · simp only [h, ite_true]
    rw [Finset.sum_comm]
    simp
  · simp [h]

/-- Exact for time-integrated flux and endpoint energy increments. -/
theorem coarse_continuity (c : X → Y) (delta : X → M) (J : X → X → M)
    (fine : ∀ x, delta x + ∑ y, J x y = 0) (a : Y) :
    push c delta a + ∑ b, coarse c J a b = 0 := by
  rw [divergence_natural, ← push_add]
  simp [push, fine]

omit [Fintype Z] in
theorem coarse_compose (c : X → Y) (d : Y → Z) (J : X → X → M) :
    coarse d (coarse c J) = coarse (d ∘ c) J := by
  funext a b
  change push d (fun y => push d (coarse c J y) b) a = _
  have hinner : ∀ y, push d (coarse c J y) b =
      push c (fun x => push (d ∘ c) (J x) b) y := by
    intro y
    simp only [push, coarse, Function.comp_apply]
    simp_rw [ite_sum]
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro x _
    by_cases h : c x = y
    · simp only [h, ite_true]
      simpa only [push, Function.comp_apply, ite_sum] using
        congrFun (push_compose c d (J x)) b
    · simp [h]
  simp_rw [hinner]
  exact congrFun (push_compose c d _) a

omit [Fintype Y] in
theorem coarse_antisymmetric (c : X → Y) (J : X → X → M)
    (anti : ∀ x y, J x y = -J y x) (a b : Y) :
    coarse c J a b = -coarse c J b a := by
  simp only [coarse, push, ite_sum, ← Finset.sum_neg_distrib]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro y _
  apply Finset.sum_congr rfl
  intro x _
  by_cases ha : c x = a <;> by_cases hb : c y = b <;> simp only [ha, hb, ite_true, ite_false, neg_zero]
  exact anti x y

omit [Fintype Y] in
/-- For real energy expectations the internal block flux is zero. -/
theorem internal_block_flux_zero (c : X → Y) (J : X → X → ℝ)
    (anti : ∀ x y, J x y = -J y x) (a : Y) :
    coarse c J a a = 0 := by
  have h := coarse_antisymmetric c J anti a a
  linarith

/-- Time-block endpoint increments telescope exactly. Integrals of face
currents add over the same intervals; no endpoint-current rule is used. -/
theorem time_increments_telescope (energy : ℕ → M) (n : ℕ) :
    ∑ k ∈ Finset.range n, (energy (k+1)-energy k) = energy n-energy 0 := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ, ih]
    abel

#print axioms coarse_antisymmetric
#print axioms internal_block_flux_zero
#print axioms time_increments_telescope
#print axioms divergence_natural
#print axioms coarse_continuity
#print axioms push_compose
#print axioms coarse_compose

end OPH.QFT.CoarseEnergyWard
