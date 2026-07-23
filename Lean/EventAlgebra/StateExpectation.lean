import EventAlgebra.Basic

/-!
# The bundled expectation functional of a state

The complex-linear functional `M ↦ Tr(ρ M)` attached to a matrix `ρ`:

* it is additive and homogeneous in the observable;
* it is **positive**: nonnegative on every positive-semidefinite
  observable (not merely on events);
* it is **normalised**: it sends `1` to `1`;
* the Born weight is its restriction to events
  (`bornWeight_eq_expectation`).

## Tagging convention

As in `EventAlgebra.Basic`. Every lemma here is trace-dependent: the
functional *is* the trace pairing.
-/

namespace EventAlgebra

open Matrix
open scoped ComplexOrder MatrixOrder

variable {n : ℕ}

/-- The trace pairing with `ρ`, bundled as a complex-linear functional in
the observable argument. This is the canonical map underlying Born weights. -/
noncomputable def stateExpectationLinearMap
    (ρ : Matrix (Fin n) (Fin n) ℂ) :
    Matrix (Fin n) (Fin n) ℂ →ₗ[ℂ] ℂ where
  toFun := bornWeight ρ
  map_add' M N := bornWeight_add ρ M N
  map_smul' c M := by
    simp only [bornWeight, mul_smul_comm, trace_smul, smul_eq_mul,
      RingHom.id_apply]

/-- Semantic abbreviation for evaluating the bundled state expectation. -/
noncomputable abbrev expectation
    (ρ M : Matrix (Fin n) (Fin n) ℂ) : ℂ :=
  stateExpectationLinearMap ρ M

@[simp]
theorem stateExpectationLinearMap_apply
    (ρ M : Matrix (Fin n) (Fin n) ℂ) :
    stateExpectationLinearMap ρ M = bornWeight ρ M :=
  rfl

/-- **Trace-dependent.** The Born weight of an event is the
expectation functional evaluated at it. -/
theorem bornWeight_eq_expectation (ρ P : Matrix (Fin n) (Fin n) ℂ) :
    bornWeight ρ P = expectation ρ P :=
  rfl

/-- **Trace-dependent.** Normalisation: the expectation of the
identity observable is `1`. -/
theorem expectation_one {ρ : Matrix (Fin n) (Fin n) ℂ} (hρ : IsState ρ) :
    expectation ρ 1 = 1 := by
  change bornWeight ρ 1 = 1
  exact bornWeight_one hρ

/-- **Trace-dependent.** Positivity of the expectation functional
on all positive-semidefinite observables, in the partial order of `ℂ`: by
the spectral theorem `M = V D V⋆` with `D` a nonnegative diagonal, and
cycling the trace, `Tr(ρ M) = Tr((V⋆ ρ V) D)` is a sum of products of
nonnegative diagonal entries. This strengthens `bornWeight_nonneg` from
events to arbitrary positive-semidefinite observables. -/
theorem expectation_nonneg {ρ M : Matrix (Fin n) (Fin n) ℂ}
    (hρ : ρ.PosSemidef) (hM : M.PosSemidef) : 0 ≤ expectation ρ M := by
  classical
  change 0 ≤ bornWeight ρ M
  set V : Matrix (Fin n) (Fin n) ℂ :=
    (hM.1.eigenvectorUnitary : Matrix (Fin n) (Fin n) ℂ) with hV
  set d : Fin n → ℂ := RCLike.ofReal ∘ hM.1.eigenvalues with hd
  -- Spectral decomposition `M = V D V⋆` with nonnegative diagonal `D`.
  have hspec : M = V * diagonal d * star V := by
    rw [hV, hd]
    conv_lhs => rw [hM.1.spectral_theorem, Unitary.conjStarAlgAut_apply]
  -- The sandwiched state `V⋆ ρ V` is positive semidefinite.
  have hsandwich : (star V * ρ * V).PosSemidef := by
    have := hρ.conjTranspose_mul_mul_same V
    rwa [← star_eq_conjTranspose] at this
  -- Cycle the trace: `Tr(ρ M) = Tr((V⋆ ρ V) D)`.
  have hcycle : bornWeight ρ M = ((star V * ρ * V) * diagonal d).trace := by
    rw [bornWeight, hspec,
      show ρ * (V * diagonal d * star V) = (ρ * V * diagonal d) * star V by
        simp only [mul_assoc],
      trace_mul_comm]
    simp only [mul_assoc]
  rw [hcycle]
  -- The trace against a nonnegative diagonal is a sum of nonnegative terms.
  simp only [Matrix.trace, Matrix.diag, Matrix.mul_diagonal]
  refine Finset.sum_nonneg fun i _ => mul_nonneg hsandwich.diag_nonneg ?_
  rw [hd]
  simp only [Function.comp_apply]
  exact RCLike.ofReal_nonneg.mpr (hM.eigenvalues_nonneg i)

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms bornWeight_eq_expectation
#print axioms stateExpectationLinearMap_apply
#print axioms expectation_one
#print axioms expectation_nonneg

end EventAlgebra
