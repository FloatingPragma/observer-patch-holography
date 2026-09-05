import Mathlib.LinearAlgebra.BilinearForm.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic

/-!
# Exact time-prism coefficient and action-defect algebra

For affine magnetic fields on a supplied uniform time grid, this module
proves the polynomial coefficients of their quadratic time averages, the
endpoint-quadrature action defect, and the residual against temporal hat
functions. The pairings may be actual spatial Whitney mass pairings; they
are not assumed to be counting metrics. No positivity is needed for these
identities, while the quadratic identities explicitly require symmetry.

These are algebraic coefficient identities. Analytic integration,
distributional identities, the source-current impulse functional, spatial
Whitney interpolation, physical time and metric attachment, and convergence
are not formalized here. In particular, exact discrete stationarity does
not assert that the exactly integrated interpolant has zero residual.
-/

namespace OPH.WhitneyTimeBridge

open scoped BigOperators
noncomputable section

variable {V : Type*} [AddCommGroup V] [Module ℝ V]

/-- The coefficient functional obtained by assigning `1, 1/2, 1/3` to
the monomials `1, s, s^2`. The name makes no analytic integration claim. -/
def quadraticMeanCoefficient (c₀ c₁ c₂ : ℝ) : ℝ :=
  c₀ + c₁ / 2 + c₂ / 3

/-- Coefficient mean of a symmetric quadratic form on affine endpoint data. -/
def magneticMeanCoefficient (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (x y : V) : ℝ :=
  (β x x + β x y + β y y) / 3

/-- Polarized affine mean; no symmetry assumption is needed. -/
def polarizedMeanCoefficient (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ)
    (x₀ x₁ a₀ a₁ : V) : ℝ :=
  (2 * β x₀ a₀ + β x₀ a₁ + β x₁ a₀ + 2 * β x₁ a₁) / 6

/-- The full scalar polynomial identity before taking its coefficient mean. -/
theorem affine_quadratic_polynomial
    (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ)
    (hβ : ∀ x y, β x y = β y x) (x y : V) (s : ℝ) :
    β ((1 - s) • x + s • y) ((1 - s) • x + s • y) =
      (1 - s) ^ 2 * β x x + 2 * s * (1 - s) * β x y + s ^ 2 * β y y := by
  simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply,
    smul_eq_mul]
  rw [hβ y x]
  ring

/-- Exact symmetric quadratic coefficient mean. -/
theorem affine_quadratic_mean
    (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (x y : V) :
    quadraticMeanCoefficient (β x x)
      (2 * β x y - 2 * β x x)
      (β x x - 2 * β x y + β y y) = magneticMeanCoefficient β x y := by
  unfold quadraticMeanCoefficient magneticMeanCoefficient
  ring

/-- The polarized polynomial has the displayed coefficient mean. -/
theorem affine_polarized_mean
    (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (x₀ x₁ a₀ a₁ : V) :
    quadraticMeanCoefficient (β x₀ a₀)
      (β x₀ a₁ + β x₁ a₀ - 2 * β x₀ a₀)
      (β x₀ a₀ - β x₀ a₁ - β x₁ a₀ + β x₁ a₁) =
      polarizedMeanCoefficient β x₀ x₁ a₀ a₁ := by
  unfold quadraticMeanCoefficient polarizedMeanCoefficient
  ring

/-- Right-endpoint magnetic quadrature minus the affine magnetic mean,
including the factor `1/2` in the field action. -/
theorem magnetic_endpoint_defect
    (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ)
    (hβ : ∀ x y, β x y = β y x) (h : ℝ) (x y : V) :
    h / 2 * β y y - h / 2 * magneticMeanCoefficient β x y =
      h / 4 * (β y y - β x x) + h / 12 * β (y - x) (y - x) := by
  unfold magneticMeanCoefficient
  simp only [map_sub, LinearMap.sub_apply]
  rw [hβ y x]
  ring

/-- Telescoping is proved directly so the action endpoint convention is
independent of an imported interval-sum convention. -/
theorem sum_successive_difference (f : ℕ → ℝ) (N : ℕ) :
    (∑ n ∈ Finset.range N, (f (n + 1) - f n)) = f N - f 0 := by
  induction N with
  | zero => simp
  | succ N ih =>
      rw [Finset.sum_range_succ, ih]
      ring

/-- Exact action-defect sum for `N` slabs, with endpoint data at `0` and `N`.
Electric and source terms cancel only when the two compared actions use
identical electric and source functionals; this theorem is the magnetic
part and makes that cancellation no hidden hypothesis. -/
theorem magnetic_window_action_defect
    (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ)
    (hβ : ∀ x y, β x y = β y x) (h : ℝ) (B : ℕ → V) (N : ℕ) :
    (∑ n ∈ Finset.range N,
      (h / 2 * β (B (n + 1)) (B (n + 1)) -
        h / 2 * magneticMeanCoefficient β (B n) (B (n + 1)))) =
      h / 4 * (β (B N) (B N) - β (B 0) (B 0)) +
        h / 12 * ∑ n ∈ Finset.range N,
          β (B (n + 1) - B n) (B (n + 1) - B n) := by
  simp_rw [magnetic_endpoint_defect β hβ h]
  rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum,
    sum_successive_difference (fun n ↦ β (B n) (B n)) N]

/-- The local magnetic correction is `h^3/12` times the quadratic form of
the divided difference. Only the nonzero-step hypothesis permits division. -/
theorem magnetic_defect_divided_difference
    (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (h : ℝ) (hh : h ≠ 0) (x y : V) :
    h / 12 * β (y - x) (y - x) =
      h ^ 3 / 12 * β (h⁻¹ • (y - x)) (h⁻¹ • (y - x)) := by
  simp only [map_smul, LinearMap.smul_apply, smul_eq_mul]
  field_simp

/-- The trapezoidal polarized pairing differs from the affine mean by
exactly one sixth of the product of endpoint differences. -/
theorem polarized_trapezoid_defect
    (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (h : ℝ) (x₀ x₁ a₀ a₁ : V) :
    h / 2 * (β x₀ a₀ + β x₁ a₁) -
      h * polarizedMeanCoefficient β x₀ x₁ a₀ a₁ =
        h / 6 * β (x₁ - x₀) (a₁ - a₀) := by
  unfold polarizedMeanCoefficient
  simp only [map_sub, LinearMap.sub_apply]
  ring

/-- The local weak correction has the exact `h^2` scaling when written
as `h` times a pairing of divided differences. -/
theorem polarized_defect_divided_difference
    (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (h : ℝ) (hh : h ≠ 0)
    (x₀ x₁ a₀ a₁ : V) :
    h / 6 * β (x₁ - x₀) (a₁ - a₀) =
      h ^ 2 / 6 * (h *
        β (h⁻¹ • (x₁ - x₀)) (h⁻¹ • (a₁ - a₀))) := by
  simp only [map_smul, LinearMap.smul_apply, smul_eq_mul]
  field_simp

/-- Two adjacent affine slabs tested against the hat which is `a` at the
middle node and zero at both neighboring nodes. -/
theorem temporal_hat_mean
    (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (h : ℝ) (x₀ x₁ x₂ a : V) :
    h * polarizedMeanCoefficient β x₀ x₁ 0 a +
      h * polarizedMeanCoefficient β x₁ x₂ a 0 =
        h / 6 * (β x₀ a + 4 * β x₁ a + β x₂ a) := by
  simp only [polarizedMeanCoefficient, map_zero]
  ring

/-- A middle-node impulse minus the two-slab affine hat pairing leaves
the magnetic second-difference residual. This is generally nonzero. -/
theorem temporal_hat_residual
    (β : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (h : ℝ) (x₀ x₁ x₂ a : V) :
    h * β x₁ a - h / 6 * (β x₀ a + 4 * β x₁ a + β x₂ a) =
      h / 6 * β ((2 : ℝ) • x₁ - x₀ - x₂) a := by
  simp only [map_sub, map_smul, LinearMap.sub_apply, LinearMap.smul_apply,
    smul_eq_mul]
  ring

/-- Scalar decomposition of the full weak Ampere coefficient. Each
coordinate of a mass-matrix electric difference, a current covector, and
a magnetic adjoint image may be substituted independently. No discrete
equation and no equality of spatial mass and counting matrices is assumed. -/
theorem ampere_coefficient_decomposition
    (h e₀ e₁ j b₀ b₁ b₂ : ℝ) :
    e₁ - e₀ + h * j - h / 6 * (b₀ + 4 * b₁ + b₂) =
      (e₁ - e₀ - h * (b₁ - j)) + h / 6 * (2 * b₁ - b₀ - b₂) := by
  ring

/-- The same identity for every covector coordinate, including coordinates
outside the carrier's boundary-edge subspace. Applications substitute mass
matrix images and the full magnetic adjoint images into these slots. -/
theorem ampere_covector_decomposition {ι : Type*}
    (h : ℝ) (e₀ e₁ j b₀ b₁ b₂ : ι → ℝ) :
    e₁ - e₀ + h • j - (h / 6) • (b₀ + (4 : ℝ) • b₁ + b₂) =
      (e₁ - e₀ - h • (b₁ - j)) +
        (h / 6) • ((2 : ℝ) • b₁ - b₀ - b₂) := by
  funext i
  simp only [Pi.sub_apply, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  ring

/-- Even exact discrete Ampere leaves the explicit temporal correction. -/
theorem ampere_coefficient_of_discrete
    (h e₀ e₁ j b₀ b₁ b₂ : ℝ)
    (hAmp : e₁ - e₀ = h * (b₁ - j)) :
    e₁ - e₀ + h * j - h / 6 * (b₀ + 4 * b₁ + b₂) =
      h / 6 * (2 * b₁ - b₀ - b₂) := by
  rw [ampere_coefficient_decomposition, hAmp]
  ring

/-- A stable scalar oscillator (`h=1/2`, stiffness `1`) is an explicit
counterexample to promoting discrete Ampere to zero hat residual. -/
theorem scalar_nonzero_residual_control :
    (1 / 2 : ℝ) ^ 2 < 4 ∧
    (1 : ℝ) - 1 = -(1 / 2) * 0 ∧
    (3 / 4 : ℝ) - 1 = -(1 / 2) * (1 / 2) ∧
    (1 / 2 : ℝ) - 0 = (1 / 2) * (1 - 0) ∧
    (1 / 2 : ℝ) - 0 + (1 / 2) * 0 -
      (1 / 2) / 6 * (1 + 4 * 1 + 3 / 4) = 1 / 48 ∧
    (1 / 48 : ℝ) ≠ 0 := by
  norm_num

end

end OPH.WhitneyTimeBridge

#print axioms OPH.WhitneyTimeBridge.affine_quadratic_polynomial
#print axioms OPH.WhitneyTimeBridge.magnetic_window_action_defect
#print axioms OPH.WhitneyTimeBridge.magnetic_defect_divided_difference
#print axioms OPH.WhitneyTimeBridge.polarized_trapezoid_defect
#print axioms OPH.WhitneyTimeBridge.polarized_defect_divided_difference
#print axioms OPH.WhitneyTimeBridge.temporal_hat_mean
#print axioms OPH.WhitneyTimeBridge.temporal_hat_residual
#print axioms OPH.WhitneyTimeBridge.ampere_coefficient_decomposition
#print axioms OPH.WhitneyTimeBridge.ampere_covector_decomposition
#print axioms OPH.WhitneyTimeBridge.scalar_nonzero_residual_control
