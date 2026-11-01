import QFT.CoarseEnergyWard

/-!
# Forced momentum cells and their block transport

The actual Weyl/contact densities, kinetic face stresses and contact forces
are attached by the exact Clifford certificate in
`model/lattice_momentum_balance`. The finite sum theorems here retain the
force through coarse-graining, and exclude replacing a nonzero total force
by a closed internal divergence. They do not assert a conserved interacting
Lorentz stress tensor. The source CAR matrices are not encoded as Lean data.
-/
namespace OPH.QFT.LatticeMomentumBalance

open OPH.QFT.CoarseEnergyWard

variable {X Y : Type*} [Fintype X] [Fintype Y] [DecidableEq Y]

theorem coarse_forced_balance (c : X → Y) (delta force : X → ℝ)
    (J : X → X → ℝ) (fine : ∀ x, delta x + ∑ y, J x y = force x) (a : Y) :
    push c delta a + ∑ b, coarse c J a b = push c force a := by
  rw [divergence_natural, ← push_add]
  simp_rw [fine]

omit [Fintype Y] [DecidableEq Y] in
theorem antisymmetric_total_divergence (J : X → X → ℝ)
    (anti : ∀ x y, J x y = -J y x) : ∑ x, ∑ y, J x y = 0 := by
  have h : (∑ x, ∑ y, J x y) = -(∑ x, ∑ y, J x y) := by
    calc
      (∑ x, ∑ y, J x y) = ∑ y, ∑ x, J x y := Finset.sum_comm
      _ = ∑ y, ∑ x, -J y x := by
        apply Finset.sum_congr rfl
        intro y _
        apply Finset.sum_congr rfl
        intro x _
        exact anti x y
      _ = -(∑ x, ∑ y, J x y) := by simp only [Finset.sum_neg_distrib]
  linarith

omit [Fintype Y] [DecidableEq Y] in
theorem total_momentum_change (delta force : X → ℝ) (J : X → X → ℝ)
    (anti : ∀ x y, J x y = -J y x)
    (balance : ∀ x, delta x + ∑ y, J x y = force x) :
    ∑ x, delta x = ∑ x, force x := by
  have h := congrArg (fun f : X → ℝ => ∑ x, f x) (funext balance)
  simpa only [Finset.sum_add_distrib, antisymmetric_total_divergence J anti,
    add_zero] using h

omit [Fintype Y] [DecidableEq Y] in
/-- A nonzero summed force cannot become an internal face flux with closed
external boundary. This excludes relabelling, not a new dressed observable. -/
theorem force_not_closed_divergence (force : X → ℝ)
    (nonzero : ∑ x, force x ≠ 0) :
    ¬ ∃ J : X → X → ℝ, (∀ x y, J x y = -J y x) ∧
        (∀ x, ∑ y, J x y = force x) := by
  rintro ⟨J, anti, h⟩
  have total := antisymmetric_total_divergence J anti
  simp_rw [h] at total
  exact nonzero total

#print axioms coarse_forced_balance
#print axioms antisymmetric_total_divergence
#print axioms total_momentum_change
#print axioms force_not_closed_divergence

end OPH.QFT.LatticeMomentumBalance
