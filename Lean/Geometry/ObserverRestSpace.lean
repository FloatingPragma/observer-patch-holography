import Geometry.CelestialNullCone
import Geometry.ObserverFrameHyperboloid

/-!
# C1: three-dimensional Euclidean rest spaces

For every point of the internal future unit-timelike hyperboloid, this file
defines its Lorentz-orthogonal rest subspace, proves its real dimension is
three, and proves that minus the ambient Lorentz pairing is positive definite.
It does not formalize a manifold tangent-bundle identification or assert a
physical attachment.
-/

namespace OPH.C1Lorentz

noncomputable section

/-- Lorentz pairing with a fixed left argument, as a real linear functional. -/
def lorentzFunctional (u : Herm2) : Herm2 →ₗ[ℝ] ℝ where
  toFun := lorentzB u
  map_add' := by
    intro v w
    simp only [lorentzB, spatialDot, Prod.fst_add, Prod.snd_add, Pi.add_apply]
    have hsum : (∑ i, u.2 i * (v.2 i + w.2 i)) =
        (∑ i, u.2 i * v.2 i) + ∑ i, u.2 i * w.2 i := by
      rw [← Finset.sum_add_distrib]
      apply Finset.sum_congr rfl
      intro i _
      ring
    rw [hsum]
    ring
  map_smul' := by
    intro a v
    simp only [lorentzB, spatialDot, RingHom.id_apply, Prod.smul_fst,
      Prod.smul_snd, Pi.smul_apply, smul_eq_mul]
    have hsum : (∑ i, u.2 i * (a * v.2 i)) = a * ∑ i, u.2 i * v.2 i := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro i _
      ring
    rw [hsum]
    ring

/-- The Lorentz-orthogonal rest subspace at a frame. -/
def restSubspace (u : FrameHyperboloid) : Submodule ℝ Herm2 :=
  LinearMap.ker (lorentzFunctional u.1)

abbrev RestSpace (u : FrameHyperboloid) := restSubspace u

theorem lorentzFunctional_ne_zero (u : FrameHyperboloid) :
    lorentzFunctional u.1 ≠ 0 := by
  intro hzero
  have h := LinearMap.congr_fun hzero ((1 : ℝ), (0 : Spatial))
  simp [lorentzFunctional, lorentzB, spatialDot] at h
  exact (ne_of_gt u.2.2) h

/-- Every hyperboloid frame has an exactly three-dimensional rest space. -/
theorem finrank_restSpace (u : FrameHyperboloid) :
    Module.finrank ℝ (RestSpace u) = 3 := by
  have h : Module.finrank ℝ (RestSpace u) + 1 = 4 := by
    simpa [RestSpace, restSubspace, finrank_Herm2] using
      (Module.Dual.finrank_ker_add_one_of_ne_zero
        (lorentzFunctional_ne_zero u))
  omega

/-- Cauchy--Schwarz in the explicit three-coordinate realization. -/
theorem spatial_cauchy (x y : Spatial) :
    (spatialDot x y) ^ 2 ≤ spatialNormSq x * spatialNormSq y := by
  simp only [spatialDot, spatialNormSq, Fin.sum_univ_three]
  nlinarith [sq_nonneg (x 0 * y 1 - x 1 * y 0),
    sq_nonneg (x 0 * y 2 - x 2 * y 0),
    sq_nonneg (x 1 * y 2 - x 2 * y 1)]

/-- The positive-sign convention on a rest space: minus the ambient Lorentz
pairing. -/
def restMetric (u : FrameHyperboloid) (v w : RestSpace u) : ℝ :=
  -lorentzB v.1 w.1

theorem restMetric_symm (u : FrameHyperboloid) (v w : RestSpace u) :
    restMetric u v w = restMetric u w v := by
  simp [restMetric, lorentzB_symm]

theorem restMetric_add_right (u : FrameHyperboloid) (v w z : RestSpace u) :
    restMetric u v (w + z) = restMetric u v w + restMetric u v z := by
  simp only [restMetric, Submodule.coe_add, lorentzB_add_right]
  ring

theorem restMetric_smul_right (u : FrameHyperboloid) (a : ℝ)
    (v w : RestSpace u) :
    restMetric u v (a • w) = a * restMetric u v w := by
  simp only [restMetric, SetLike.val_smul, lorentzB_smul_right]
  ring

theorem restMetric_add_left (u : FrameHyperboloid) (v w z : RestSpace u) :
    restMetric u (v + w) z = restMetric u v z + restMetric u w z := by
  simp only [restMetric, Submodule.coe_add, lorentzB_add_left]
  ring

theorem restMetric_smul_left (u : FrameHyperboloid) (a : ℝ)
    (v w : RestSpace u) :
    restMetric u (a • v) w = a * restMetric u v w := by
  simp only [restMetric, SetLike.val_smul, lorentzB_smul_left]
  ring

/-- The sign-reversed Lorentz form is positive definite on every rest space. -/
theorem restMetric_pos (u : FrameHyperboloid) {v : RestSpace u} (hv : v ≠ 0) :
    0 < restMetric u v v := by
  have hut : u.1.1 ≠ 0 := ne_of_gt u.2.2
  have hsp : v.1.2 ≠ 0 := by
    intro hsp0
    apply hv
    apply Subtype.ext
    apply Prod.ext
    · have hort : lorentzFunctional u.1 v.1 = 0 := v.2
      change lorentzB u.1 v.1 = 0 at hort
      change u.1.1 * v.1.1 - spatialDot u.1.2 v.1.2 = 0 at hort
      rw [hsp0] at hort
      simp [spatialDot] at hort
      exact hort.resolve_left hut
    · exact hsp0
  have hrpos : 0 < spatialNormSq v.1.2 := spatialNormSq_pos hsp
  have hpnonneg : 0 ≤ spatialNormSq u.1.2 :=
    Finset.sum_nonneg fun i _ ↦ sq_nonneg (u.1.2 i)
  have hcauchy := spatial_cauchy u.1.2 v.1.2
  have hframe := u.2.1
  have hort : lorentzFunctional u.1 v.1 = 0 := v.2
  change lorentzB u.1 v.1 = 0 at hort
  have hdot : spatialDot u.1.2 v.1.2 = u.1.1 * v.1.1 := by
    simp only [lorentzB] at hort
    linarith
  rw [hdot] at hcauchy
  simp only [lorentzQ] at hframe
  simp only [restMetric, lorentzB_self, lorentzQ]
  nlinarith [sq_nonneg v.1.1]

/-- Zero is the only vector of zero rest norm. -/
theorem restMetric_self_eq_zero_iff (u : FrameHyperboloid) (v : RestSpace u) :
    restMetric u v v = 0 ↔ v = 0 := by
  constructor
  · intro hzero
    by_contra hv
    exact (ne_of_gt (restMetric_pos u hv)) hzero
  · rintro rfl
    simp [restMetric, lorentzB, spatialDot]

#print axioms finrank_restSpace
#print axioms restMetric_add_right
#print axioms restMetric_smul_right
#print axioms restMetric_symm
#print axioms restMetric_pos
#print axioms restMetric_self_eq_zero_iff

end

end OPH.C1Lorentz
