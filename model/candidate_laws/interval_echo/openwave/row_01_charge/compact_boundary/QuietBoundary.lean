import Mathlib.Analysis.Complex.RealDeriv
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Data.Matrix.Basic
import Mathlib.Tactic

namespace IntervalEcho.QuietBoundary

section Parked

variable {W P : Type*} [Fintype W] [Fintype P] [DecidableEq W] [DecidableEq P]

/-- In the parked aperture, measure.closed_certificate verifies HE = EK in
integer eighth-energy units. E inserts the compact source orbital, and K includes
its scalar energy in addition to the exterior probe Hamiltonian. -/
theorem parked_power (H : Matrix W W ℂ) (E : Matrix W P ℂ) (K : Matrix P P ℂ)
    (parked : H * E = E * K) (n : ℕ) : H ^ n * E = E * K ^ n := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [pow_succ, pow_succ, Matrix.mul_assoc, parked, ← Matrix.mul_assoc, ih,
        Matrix.mul_assoc]

end Parked

/-- A record branch can contain any finite list of operations preserving the
same compact encoding. Intermediate exterior instruments can select branches;
the assumptions have to hold on each enrolled instrument branch. -/
theorem preserved_branch {Tag X Y : Type*} (E : X → Y)
    (small : Tag → X → X) (full : Tag → Y → Y)
    (preserves : ∀ t x, full t (E x) = E (small t x)) (word : List Tag) (x : X) :
    word.foldl (fun y t => full t y) (E x) = E (word.foldl (fun z t => small t z) x) := by
  induction word generalizing x with
  | nil => rfl
  | cons t rest ih =>
      simp only [List.foldl_cons, preserves, ih]

section WeakOpening

variable {W : Type*} [Fintype W]

/-- A sharp load record, or any fixed classical postprocessing or probabilistic
coarsening of it, is diagonal in the complete word readout. -/
def record (weight : W → ℝ) (x : W → ℂ) : ℝ := ∑ w, weight w * (x w).normSq

/-- A source-dependent scalar phase in the encoded exterior amplitude cancels
from each exterior number record. -/
theorem phase_record (weight : W → ℝ) (x : W → ℂ) (phase : ℂ)
    (unitPhase : phase.normSq = 1) :
    record weight (fun w => phase * x w) = record weight x := by
  simp [record, Complex.normSq_mul, unitPhase]

/-- The parked packet and its first variation have disjoint internal-count
supports. On that domain the number readout of the affine tangent is exactly
quadratic. The finite-opening packet can also have a second-order amplitude
on the original support and terms of higher order. -/
theorem diagonal_tangent (weight : W → ℝ) (x v : W → ℂ)
    (separated : ∀ w, x w = 0 ∨ v w = 0) (g : ℝ) :
    record weight (fun w => x w + (g : ℂ) * v w) =
      record weight x + g ^ 2 * record weight v := by
  unfold record
  rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun w _ => ?_
  rcases separated w with hx | hv
  · simp only [hx, zero_add, Complex.normSq_mul, Complex.normSq_ofReal, Complex.normSq_zero,
      mul_zero]
    ring
  · simp [hv]

/-- For differentiable complex word amplitudes the actual number-record
derivative vanishes, independently of the second-order dynamics. -/
theorem zero_first_derivative (weight : W → ℝ) (path : W → ℝ → ℂ) (v : W → ℂ)
    (hv : ∀ w, HasDerivAt (path w) (v w) 0)
    (separated : ∀ w, path w 0 = 0 ∨ v w = 0) :
    HasDerivAt (fun g => record weight (fun w => path w g)) 0 0 := by
  have hw (w : W) :
      HasDerivAt (fun g : ℝ => weight w * (path w g).normSq) 0 0 := by
    have hr : HasDerivAt (fun g => (path w g).re) (v w).re 0 :=
      Complex.reCLM.hasFDerivAt.comp_hasDerivAt 0 (hv w)
    have hi : HasDerivAt (fun g => (path w g).im) (v w).im 0 :=
      Complex.imCLM.hasFDerivAt.comp_hasDerivAt 0 (hv w)
    have h := ((hr.mul hr).add (hi.mul hi)).const_mul (weight w)
    rcases separated w with hx | hvw
    · simpa [Complex.normSq_apply, hx] using h
    · simpa [Complex.normSq_apply, hvw] using h
  simpa only [record, Finset.sum_const_zero] using
    (HasDerivAt.fun_sum (u := Finset.univ) fun w _ => hw w)

end WeakOpening

end IntervalEcho.QuietBoundary

#print axioms IntervalEcho.QuietBoundary.parked_power
#print axioms IntervalEcho.QuietBoundary.preserved_branch
#print axioms IntervalEcho.QuietBoundary.phase_record
#print axioms IntervalEcho.QuietBoundary.diagonal_tangent
#print axioms IntervalEcho.QuietBoundary.zero_first_derivative
