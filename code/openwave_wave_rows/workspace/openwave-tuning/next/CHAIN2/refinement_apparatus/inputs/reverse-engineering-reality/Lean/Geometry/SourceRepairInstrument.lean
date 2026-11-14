import Dynamics.ChoiCPTP
import Geometry.SourceOperationReads

/-!
# Quantum extensions and obstructions for the registered scalar repair

Permutation twirls and diagonal pinching supply completely positive,
trace-preserving extensions of ideal diagonal means. Affinity obstructs the
native threshold and phase-preserving lift. These extensions are explicit
mathematical controls, not an axiom-selected replacement for the simulator.
-/

set_option autoImplicit false

namespace OPH.SourceRepairInstrument
noncomputable section
open Matrix OPH.Dynamics
open scoped ComplexOrder

/-- Conjugation by a permutation, written without a choice of matrix entries. -/
def relabel {n : ℕ} (e : Fin n ≃ Fin n) : CMat n →ₗ[ℂ] CMat n where
  toFun X := X.submatrix e e
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

/-- Every permutation conjugation is completely positive and trace preserving. -/
theorem relabel_cptp {n : ℕ} (e : Fin n ≃ Fin n) : IsCPTP (relabel e) := by
  constructor
  · intro m X hX
    change (X.submatrix (fun p : Fin m × Fin n => (p.1, e p.2))
      (fun p : Fin m × Fin n => (p.1, e p.2))).PosSemidef
    exact hX.submatrix _
  · intro X
    exact e.sum_comp (fun i => X i i)

/-- A retained permutation label supplies the exact inverse on every matrix. -/
theorem relabel_recovery {n : ℕ} (e : Fin n ≃ Fin n) (X : CMat n) :
    relabel e.symm (relabel e X) = X := by
  ext i j
  simp [relabel, Matrix.submatrix]

/-- Recovering each labelled branch before discarding the label preserves input. -/
theorem retained_branch_recovery {n : ℕ} {ι : Type*} [Fintype ι]
    (e : ι → Fin n ≃ Fin n) (w : ι → ℂ) (hw : ∑ i, w i = 1) (X : CMat n) :
    (∑ i, relabel (e i).symm (w i • relabel (e i) X)) = X := by
  simp only [map_smul, relabel_recovery]
  rw [← Finset.sum_smul, hw, one_smul]

/-- Retaining one signed difference reverses a deterministic pair mean. -/
theorem deterministic_pair_recovery (a b : ℝ) :
    ((a+b)/2+(a-b)/2, (a+b)/2-(a-b)/2) = (a,b) := by
  ext <;> ring

/-- Equal mixture of doing nothing and conjugating by the permutation. -/
def twirl {n : ℕ} (e : Fin n ≃ Fin n) : CMat n →ₗ[ℂ] CMat n :=
  (1 / 2 : ℂ) • LinearMap.id + (1 / 2 : ℂ) • relabel e

/-- The ideal repair extension has a genuine CPTP proof. -/
theorem twirl_cptp {n : ℕ} (e : Fin n ≃ Fin n) : IsCPTP (twirl e) := by
  simpa [twirl] using (isCPTP_id (n := n)).convexCombination
    (relabel_cptp e) (a := (1/2 : ℝ)) (b := (1/2 : ℝ))
    (by norm_num) (by norm_num) (by norm_num)

/-- On the diagonal the twirl is exactly the scalar pair mean. -/
theorem twirl_diagonal {n : ℕ} (e : Fin n ≃ Fin n) (X : CMat n) (i : Fin n) :
    twirl e X i i = (X i i + X (e i) (e i))/2 := by
  simp [twirl, relabel, Matrix.submatrix]
  ring

/-- An involutive repair is idempotent on the entire matrix state. -/
theorem twirl_idempotent {n : ℕ} (e : Fin n ≃ Fin n)
    (he : Function.Involutive e) (X : CMat n) :
    twirl e (twirl e X) = twirl e X := by
  ext i j
  simp [twirl, relabel, Matrix.submatrix, he i, he j]
  ring

/-- Erase all off-diagonal entries. -/
def dephase {n : ℕ} : CMat n →ₗ[ℂ] CMat n where
  toFun X := Matrix.diagonal (fun i => X i i)
  map_add' _ _ := by ext i j; by_cases h : i = j <;> simp [Matrix.diagonal, h]
  map_smul' _ _ := by ext i j; by_cases h : i = j <;> simp [Matrix.diagonal, h]

/-- Coordinate projection Kraus operator. -/
def coordinate {n : ℕ} (i : Fin n) : CMat n :=
  Matrix.diagonal (fun j => if j = i then 1 else 0)

/-- A coordinate matrix unit. -/
def matrixUnit {n : ℕ} (i j : Fin n) : CMat n :=
  Matrix.of fun a b => if a = i ∧ b = j then 1 else 0

/-- A dense time-translated port effect exposes every matrix unit. -/
theorem coordinate_sandwich {n : ℕ} (X : CMat n) (i j : Fin n) :
    coordinate i * X * coordinate j = X i j • matrixUnit i j := by
  ext a b
  by_cases ha : a = i <;> by_cases hb : b = j <;>
    simp [coordinate, matrixUnit, Matrix.diagonal, Matrix.mul_apply, ha, hb]

/-- A nonzero off-diagonal entry prevents the port projection being central. -/
theorem coordinate_commutator {n : ℕ} (X : CMat n) (i j : Fin n) (h : i ≠ j) :
    (coordinate i * X - X * coordinate i) i j = X i j := by
  simp [coordinate, Matrix.diagonal, Matrix.mul_apply, Ne.symm h]

/-- The dephasing map has the complete coordinate Kraus family. -/
theorem dephase_kraus {n : ℕ} (X : CMat n) :
    dephase X = ∑ i, coordinate i * X * (coordinate i)ᴴ := by
  classical
  ext a b
  by_cases h : a = b
  · subst b
    simp [dephase, coordinate, Matrix.mul_apply, Matrix.conjTranspose_apply,
      Matrix.diagonal, Matrix.sum_apply]
  · simp [dephase, coordinate, Matrix.mul_apply, Matrix.conjTranspose_apply,
      Matrix.diagonal, Matrix.sum_apply, h, Ne.symm h]

/-- Diagonal erasure is a channel, so it gives a second repair extension. -/
theorem dephase_cptp {n : ℕ} : IsCPTP (dephase (n := n)) := by
  constructor
  · exact isCompletelyPositive_of_kraus _ coordinate dephase_kraus
  · intro X
    simp [dephase, Matrix.trace]

/-- A measured repair is CPTP without assuming a selected instrument. -/
theorem measured_twirl_cptp {n : ℕ} (e : Fin n ≃ Fin n) :
    IsCPTP (dephase ∘ₗ twirl e) := dephase_cptp.comp (twirl_cptp e)

/-- Dephasing leaves every current diagonal read unchanged. -/
theorem dephase_diagonal {n : ℕ} (X : CMat n) (i : Fin n) :
    dephase X i i = X i i := by simp [dephase]

/-- A channel cannot distinguish two equal convex input ensembles. -/
theorem ensemble_identity {V W : Type*} [AddCommGroup V] [Module ℝ V]
    [AddCommGroup W] [Module ℝ W] (f : V →ₗ[ℝ] W)
    (a b c d : V) (h : a+b = c+d) : f a+f b = f c+f d := by
  simpa only [map_add] using congrArg f h

/-- The native no-op band interpreted on real diagonal states. -/
def thresholdRead (epsilon x : ℝ) : ℝ :=
  if |2*x-1| ≤ epsilon then x else 1/2

/-- Two nearby inputs stay unchanged while their doubled displacement is repaired. -/
theorem threshold_nonaffine (epsilon delta : ℝ) (hd : 0 < delta)
    (hlo : 2*delta ≤ epsilon) (hhi : epsilon < 4*delta) :
    thresholdRead epsilon (1/2+delta) ≠
      (thresholdRead epsilon (1/2) + thresholdRead epsilon (1/2+2*delta))/2 := by
  have epos : 0 ≤ epsilon := by linarith
  have h1 : |2*(1/2+delta)-1| ≤ epsilon := by
    rw [show 2*(1/2+delta)-1 = 2*delta by ring, abs_of_pos (by positivity)]
    exact hlo
  have h2 : ¬ |2*(1/2+2*delta)-1| ≤ epsilon := by
    rw [show 2*(1/2+2*delta)-1 = 4*delta by ring, abs_of_pos (by positivity)]
    linarith
  simp only [thresholdRead, if_pos h1, if_neg h2]
  norm_num [epos]
  linarith

/-- A single seam update changes one carrier's trace unless its values agree. -/
theorem local_normalization_change (a b : ℝ) :
    1-a+(a+b)/2 = 1+(b-a)/2 := by ring

/-- The native phase-lift ensemble witness has strictly positive coherence gap. -/
theorem phase_lift_gap {n : ℝ} (hn : 0 < n) :
    1/(8*n) < Real.sqrt 13/(24*n) := by
  have hs : 3 < Real.sqrt 13 := by
    exact (Real.lt_sqrt (by norm_num)).2 (by norm_num)
  apply (div_lt_div_iff₀ (by positivity) (by positivity)).2
  nlinarith

end
end OPH.SourceRepairInstrument
