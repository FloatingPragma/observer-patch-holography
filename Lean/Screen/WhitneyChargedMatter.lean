import Mathlib.Analysis.Complex.Trigonometric
import Mathlib.Tactic

set_option autoImplicit false
open scoped BigOperators

namespace OPH.WhitneyChargedMatter

noncomputable section

/-!
# Gauge-equivariant dressed nodal interpolation

For a Whitney one-form on a tetrahedron, the straight-path integral from
vertex `i` to barycentric point `lambda` is `sum_j a i j * lambda j`.
The finite identities below prove the gauge covariance, nodal interpolation,
and common-face algebra of the resulting dressed scalar interpolant.

The domain geometry, straight-path integral identification, analytic H1
conformity, exact integral of the scalar-QED action, variational Noether
identity, and local ODE existence are separate analytic arguments. Neither
the scalar matter content, its charge/mass/coupling, nor its quantum theory
is selected by these algebraic theorems.
-/

variable {ι : Type*}

def unitPhase (x : ℝ) : ℂ := Complex.exp (Complex.I * (x : ℂ))

theorem unitPhase_add (x y : ℝ) : unitPhase (x + y) = unitPhase x * unitPhase y := by
  simp [unitPhase, mul_add, Complex.exp_add]

theorem unitPhase_zero : unitPhase 0 = 1 := by simp [unitPhase]

theorem unitPhase_norm (x : ℝ) : ‖unitPhase x‖ = 1 := by
  simp [unitPhase, Complex.norm_exp]

def gaugeEdges (a : ι → ι → ℝ) (chi : ι → ℝ) (i j : ι) : ℝ :=
  a i j + chi j - chi i

theorem gaugeEdges_diagonal (a : ι → ι → ℝ) (chi : ι → ℝ)
    (ha : ∀ i, a i i = 0) (i : ι) : gaugeEdges a chi i i = 0 := by
  simp [gaugeEdges, ha]

theorem gaugeEdges_antisymmetric (a : ι → ι → ℝ) (chi : ι → ℝ)
    (ha : ∀ i j, a j i = -a i j) (i j : ι) :
    gaugeEdges a chi j i = -gaugeEdges a chi i j := by
  simp only [gaugeEdges, ha j i]
  rw [ha i j]
  ring

def nodalGauge (charge : ℝ) (chi : ι → ℝ) (psi : ι → ℂ) (i : ι) : ℂ :=
  unitPhase (charge * chi i) * psi i

section Finite

variable [Fintype ι]

def nodalValue (lambda chi : ι → ℝ) : ℝ := ∑ j, lambda j * chi j

def pathPhase (a : ι → ι → ℝ) (lambda : ι → ℝ) (i : ι) : ℝ :=
  ∑ j, a i j * lambda j

def interpolate (charge : ℝ) (a : ι → ι → ℝ) (lambda : ι → ℝ)
    (psi : ι → ℂ) : ℂ :=
  ∑ i, (lambda i : ℂ) * unitPhase (charge * pathPhase a lambda i) * psi i

/-- Partition of unity supplies exactly the endpoint gauge correction. -/
theorem pathPhase_gauge (a : ι → ι → ℝ) (lambda chi : ι → ℝ)
    (hlambda : ∑ j, lambda j = 1) (i : ι) :
    pathPhase (gaugeEdges a chi) lambda i =
      pathPhase a lambda i + nodalValue lambda chi - chi i := by
  simp only [pathPhase, gaugeEdges, add_mul, sub_mul,
    Finset.sum_sub_distrib, Finset.sum_add_distrib, ← Finset.mul_sum, hlambda,
    mul_one, nodalValue]
  congr 1
  congr 1
  apply Finset.sum_congr rfl
  intro j _
  ring

/-- The dressed nodal function has the ordinary pointwise U(1) gauge law. -/
theorem interpolate_gauge (charge : ℝ) (a : ι → ι → ℝ)
    (lambda chi : ι → ℝ) (psi : ι → ℂ)
    (hlambda : ∑ j, lambda j = 1) :
    interpolate charge (gaugeEdges a chi) lambda (nodalGauge charge chi psi) =
      unitPhase (charge * nodalValue lambda chi) * interpolate charge a lambda psi := by
  unfold interpolate
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  rw [pathPhase_gauge a lambda chi hlambda i]
  unfold nodalGauge
  have hp : unitPhase (charge * (pathPhase a lambda i + nodalValue lambda chi - chi i)) *
      unitPhase (charge * chi i) =
      unitPhase (charge * nodalValue lambda chi) * unitPhase (charge * pathPhase a lambda i) := by
    rw [← unitPhase_add, ← unitPhase_add]
    congr 1
    ring
  calc
    _ = (lambda i : ℂ) *
        (unitPhase (charge * (pathPhase a lambda i + nodalValue lambda chi - chi i)) *
          unitPhase (charge * chi i)) * psi i := by ring
    _ = _ := by rw [hp]; ring

theorem interpolate_gauge_norm (charge : ℝ) (a : ι → ι → ℝ)
    (lambda chi : ι → ℝ) (psi : ι → ℂ)
    (hlambda : ∑ j, lambda j = 1) :
    ‖interpolate charge (gaugeEdges a chi) lambda (nodalGauge charge chi psi)‖ =
      ‖interpolate charge a lambda psi‖ := by
  rw [interpolate_gauge charge a lambda chi psi hlambda, norm_mul, unitPhase_norm, one_mul]

variable [DecidableEq ι]

/-- At each vertex, the dressed interpolant recovers its nodal value. -/
theorem interpolate_vertex (charge : ℝ) (a : ι → ι → ℝ)
    (ha : ∀ i, a i i = 0) (psi : ι → ℂ) (k : ι) :
    interpolate charge a (fun j => if j = k then 1 else 0) psi = psi k := by
  unfold interpolate
  rw [Finset.sum_eq_single k]
  · simp [pathPhase, ha, unitPhase_zero]
  · intro j _ hj
    simp [hj]
  · simp

/-- Face traces depend only on its vertices and tangential edge coefficients.
The same global vertex carrier allows either adjacent tetrahedron's local
barycentric coordinates to be extended by zero outside that face. -/
theorem interpolate_face_trace (charge : ℝ) (a b : ι → ι → ℝ)
    (lambda : ι → ℝ) (psi phi : ι → ℂ) (face : Finset ι)
    (hlambda : ∀ j, j ∉ face → lambda j = 0)
    (hab : ∀ i ∈ face, ∀ j ∈ face, a i j = b i j)
    (hpsi : ∀ i ∈ face, psi i = phi i) :
    interpolate charge a lambda psi = interpolate charge b lambda phi := by
  unfold interpolate
  apply Finset.sum_congr rfl
  intro i _
  by_cases hi : i ∈ face
  · have hp : pathPhase a lambda i = pathPhase b lambda i := by
      unfold pathPhase
      apply Finset.sum_congr rfl
      intro j _
      by_cases hj : j ∈ face
      · rw [hab i hi j hj]
      · simp [hlambda j hj]
    rw [hp, hpsi i hi]
  · simp [hlambda i hi]

end Finite
end

#print axioms pathPhase_gauge
#print axioms interpolate_gauge
#print axioms interpolate_gauge_norm
#print axioms interpolate_vertex
#print axioms interpolate_face_trace

end OPH.WhitneyChargedMatter
