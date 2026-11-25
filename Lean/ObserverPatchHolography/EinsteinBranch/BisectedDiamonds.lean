import ObserverPatchHolography.EinsteinBranch.RayPairBalance
import ObserverPatchHolography.EinsteinBranch.QuantitativeSmallBall

/-! Coefficient identities for marked, centered four-dimensional causal diamonds.
Clock, curvature expansion and counting accuracy have separate physical hypotheses. -/
namespace OPH.EinsteinBranch.BisectedDiamonds

open RayPairBalance (contraction)

/-- The direction-independent scalar coefficient includes the curvature trace
and the contracted Hessian of the marker intensity. -/
noncomputable def full (h m r b a : ℝ) : ℝ := m + h^2 * (m*r/30 + b + a/40)

/-- Mean normalized count of the two successive half-duration diamonds. -/
noncomputable def halves (h m r b a : ℝ) : ℝ := m + h^2 * ((m*r/30 + b + a/40)/4 + a/32)

noncomputable def corrected (h m r b a : ℝ) : ℝ := full h m r b a - 2/3 * halves h m r b a

theorem cancel_density_hessian (h m r b a : ℝ) :
    corrected h m r b a = m/3 + h^2*(m*r/36 + 5*b/6) := by
  unfold corrected full halves
  ring

/-- Three geodesic directions in one tetrad give the null Ricci contraction.
Each marker Hessian may be arbitrary. The common trace terms also cancel. -/
theorem recover_null_coefficient (h m d j p b a₀ aPos aNeg : ℝ) :
    36 * (-(9/16)*corrected h m (contraction d j p 1 0) b a₀ +
      (81/160)*corrected h m (contraction d j p (5/3) (4/3)) b aPos +
      (9/160)*corrected h m (contraction d j p (5/3) (-4/3)) b aNeg)
      = m*h^2*contraction d j p 1 1 := by
  simp only [cancel_density_hessian, RayPairBalance.contraction]
  ring

/-- Max-norm perturbation of one directional full/half combination. -/
theorem count_error (f s ε : ℝ) (hf : |f| ≤ ε) (hs : |s| ≤ ε) :
    |f - 2/3*s| ≤ 5/3*ε := by
  obtain ⟨a,b⟩ := abs_le.mp hf
  obtain ⟨c,d⟩ := abs_le.mp hs
  apply abs_le.mpr
  constructor <;> linarith

noncomputable def scaledRead (h m κ d j p t f z b u v a : ℝ) : ℝ :=
  36*corrected h m (contraction d j p u v) b a - m*h^2*κ*contraction t f z u v

/-- Compare an independently read energy density on each of the same three
clocks. The geometric and material quadratic forms have unrelated entries. -/
theorem null_residual (h m κ d j p t f z b a₀ aPos aNeg : ℝ) :
    -(9/16)*scaledRead h m κ d j p t f z b 1 0 a₀ +
      (81/160)*scaledRead h m κ d j p t f z b (5/3) (4/3) aPos +
      (9/160)*scaledRead h m κ d j p t f z b (5/3) (-4/3) aNeg =
      m*h^2*(contraction d j p 1 1 - κ*contraction t f z 1 1) := by
  simp only [scaledRead, cancel_density_hessian, RayPairBalance.contraction]
  ring

/-- Arbitrary errors across the three directions may be correlated. -/
theorem three_read_error (x y z ε : ℝ) (hx : |x| ≤ ε) (hy : |y| ≤ ε)
    (hz : |z| ≤ ε) : |-(9/16)*x+(81/160)*y+(9/160)*z| ≤ (9/8)*ε := by
  obtain ⟨a,b⟩ := abs_le.mp hx
  obtain ⟨c,d'⟩ := abs_le.mp hy
  obtain ⟨e,f⟩ := abs_le.mp hz
  apply abs_le.mpr
  constructor <;> linarith

assert_no_sorry null_residual
assert_no_sorry three_read_error

assert_no_sorry cancel_density_hessian
assert_no_sorry recover_null_coefficient
assert_no_sorry count_error

end OPH.EinsteinBranch.BisectedDiamonds
