import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.LinearCombination

/-! Finite correlation energy and its response on the homogeneous constrained comparison. -/
namespace OPH.EinsteinBranch.OverlapCoupling

theorem gibbs_weight (β E Z : ℝ) (hZ : 0 < Z) :
    -Real.log (Real.exp (-β * E) / Z) = β * E + Real.log Z := by
  rw [Real.log_div (Real.exp_pos _).ne' hZ.ne', Real.log_exp]
  ring

/-- Expectation in the high preparation minus that in the low preparation. -/
noncomputable def cycle (a b c d : ℝ) := (b + c - a - d) / 4

theorem cycle_modular (β a b c d z u v x y : ℝ) :
    cycle (β*a+z+u+x) (β*b+z+u+y) (β*c+z+v+x) (β*d+z+v+y) =
    β * cycle a b c d := by unfold cycle; ring

theorem seam_energy : cycle 0 0 32 0 = 8 ∧ cycle 0 0 32 32 = 0 := by
  unfold cycle
  constructor <;> norm_num

theorem protected_register_energy : cycle 0 0 33 32 = 1/4 := by
  unfold cycle
  norm_num

/-- The mixed contrast is the complete obstruction to splitting a two-by-two read into marginals. -/
theorem zero_cycle_iff_separate (a b c d : ℝ) :
    cycle a b c d = 0 ↔ ∃ r s u v : ℝ, a=r+u ∧ b=r+v ∧ c=s+u ∧ d=s+v := by
  unfold cycle
  constructor
  · intro h
    refine ⟨a,c,0,b-a,?_,?_,?_,?_⟩ <;> linarith
  · rintro ⟨r,s,u,v,h₁,h₂,h₃,h₄⟩
    linarith

theorem volume_constraint (W h κ p E Λ : ℝ) (hW : 0 < W)
    (rate : 6*W*h = -(κ*p))
    (constraint : (κ*p)^2 = 12*(Λ*W^2+κ*W*E)) :
    W*(3*h^2-Λ) = κ*E := by
  have hs : (6*W*h)^2 = 12*(Λ*W^2+κ*W*E) := by rw [rate, neg_sq, constraint]
  have hz : W*(W*(3*h^2-Λ)-κ*E)=0 := by nlinarith only [hs]
  exact sub_eq_zero.mp ((mul_eq_zero.mp hz).resolve_left hW.ne')

theorem held_response (W₁ W₂ G₁ G₂ E₁ E₂ κ Λ β K : ℝ)
    (first : W₁*(G₁-Λ)=κ*E₁) (second : W₂*(G₂-Λ)=κ*E₂)
    (meter : K=β*(E₁-E₂)) :
    β*(W₁*(G₁-Λ)-W₂*(G₂-Λ)) = κ*K := by
  rw [first, second, meter]
  ring

theorem null_response (W G P Λ κ E : ℝ) (constraint : W*(G-Λ)=κ*E)
    (pressure : P = -Λ) : W*(G+P) = κ*E := by
  rw [pressure, ← sub_eq_add_neg]
  exact constraint

/-- Multiply the regional temporal residual by its common area denominator.
The hats use fixed proper volume. The curvature approximation error enters separately. -/
theorem regional_area_residual (π₀ ℓ κ A K T G S : ℝ)
    (area : 4*π₀*ℓ^4*G = -15*A)
    (collar : 8*π₀^2*ℓ^4*S = 15*K) :
    8*π₀^2*ℓ^4*(G-κ*(T+S)) =
      -30*π₀*A-8*π₀^2*ℓ^4*κ*T-15*κ*K := by
  linear_combination 2*π₀*area - κ*collar

/-- A positive regional charge slope and two nonnegative area slopes
give a lower bound for every nondecreasing differentiable area combination. -/
theorem signed_area_lower_bound (charge bound κ scale a b u v : ℝ)
    (hc : bound ≤ charge) (hk : 0 ≤ κ) (hs : 0 ≤ scale)
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hu : 0 ≤ u) (hv : 0 ≤ v) :
    κ*bound ≤ scale*(u*a+v*b)+κ*charge := by
  have h := mul_nonneg hs (add_nonneg (mul_nonneg hu ha) (mul_nonneg hv hb))
  linarith [mul_le_mul_of_nonneg_left hc hk]

/-- A transfer giving opposite area slopes excludes positive rescalings
of both readings into one area derivative. -/
theorem opposite_area_units (a b u v : ℝ)
    (ha : 0 < a) (hb : b < 0) (hu : 0 < u) (hv : 0 < v) : u*a ≠ v*b :=
  ne_of_gt ((mul_neg_of_pos_of_neg hv hb).trans (mul_pos hu ha))

/-- A finite first-exit clock. L is the symmetric interior graph Laplacian
with the exterior absorbing. Its row sums count incoming cut edges. -/
theorem boundary_clock_sum {ι : Type*} [Fintype ι] (L : ι → ι → ℝ)
    (t d b : ι → ℝ) (sym : ∀ i j, L i j = L j i)
    (poisson : ∀ i, ∑ j, L i j * t j = d i)
    (cut : ∀ i, ∑ j, L i j = b i) :
    ∑ i, b i * t i = ∑ i, d i := by
  simp_rw [← cut, Finset.sum_mul, ← poisson]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun i _ =>
    Finset.sum_congr rfl fun j _ => by rw [sym j i]

/-- At equal proper volume the signed area response is encoded in the
difference of the two incoming-flux mean proper chord lengths. -/
theorem equal_volume_chords (A₀ A₁ L₀ L₁ V : ℝ)
    (first : A₀ * L₀ = 4 * V) (second : A₁ * L₁ = 4 * V) :
    (A₁ - A₀) * (L₀ * L₁) = -(4 * V) * (L₁ - L₀) := by
  linear_combination L₀ * second - L₁ * first

/-- Ordered-pair admission is uniformly positive after quartic normalization
if the single-station fraction has a uniform quadratic lower bound. -/
theorem ball_pair_lower (m N α ℓ : ℝ) (hm : 2 ≤ m) (hN : m ≤ N)
    (hα : 0 ≤ α) (density : α * ℓ^2 ≤ m / N) :
    α^2 * ℓ^4 / 2 ≤ m * (m-1) / (N * (N-1)) := by
  have hn : 0 < N := by linarith
  have hsub : 0 < N - 1 := by linarith
  have hm0 : 0 ≤ m := by linarith
  have hratio : m / (2*N) ≤ (m-1)/(N-1) := by
    apply (div_le_div_iff₀ (by positivity) hsub).2
    nlinarith only [mul_nonneg hn.le (sub_nonneg.mpr hm), hm0]
  have hsq := mul_self_le_mul_self (mul_nonneg hα (sq_nonneg ℓ)) density
  calc
    α^2 * ℓ^4 / 2 ≤ (m/N)^2 / 2 := by nlinarith only [hsq]
    _ = (m/N) * (m/(2*N)) := by field_simp
    _ ≤ (m/N) * ((m-1)/(N-1)) :=
      mul_le_mul_of_nonneg_left hratio (div_nonneg hm0 hn.le)
    _ = _ := by field_simp

end OPH.EinsteinBranch.OverlapCoupling
