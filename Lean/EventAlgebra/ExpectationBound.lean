import EventAlgebra.StateExpectation
import EventAlgebra.Tsirelson

/-!
# The state-expectation bound and the state-level CHSH corollary

The bridge from the operator-norm layer to the state layer: for a state
`ρ` and any observable `M`,

  `‖Tr(ρ M)‖ ≤ ‖M‖`

in the `L2` operator norm. The proof diagonalizes the state,
`ρ = V D V⋆`, cycles the trace to `Tr(D · V⋆ M V)`, bounds each diagonal
entry of the conjugated observable by its operator norm, and uses that the
nonnegative eigenvalues of a state sum to one.

Composed with the operator-norm Tsirelson bound of
`EventAlgebra.Tsirelson`, this yields the state-level CHSH corollary: for
four projection events with cross-party commutation, the expectation of
the CHSH combination of dichotomic observables under any state has
modulus at most `2√2` (`matrix_state_tsirelson_bound_of_events`).

No attainment claim is made anywhere in this module: both statements are
upper bounds, and nothing here asserts that `2√2` is achieved by any
state.

## Tagging convention

As in `EventAlgebra.Basic`. The entry and unitary norm lemmas are
**algebra-only**; the expectation bound and the CHSH corollary are
**trace-dependent**.
-/

namespace EventAlgebra

open Matrix
open scoped ComplexOrder
open scoped ComplexInnerProductSpace
open scoped Matrix.Norms.L2Operator

variable {n : ℕ}

set_option linter.deprecated false in
/-- **Algebra-only.** Every diagonal entry of a matrix is dominated by its
`L2` operator norm. The proof tests the matrix against a standard basis
vector of the Euclidean space. -/
theorem norm_diag_entry_le_l2_opNorm (B : Matrix (Fin n) (Fin n) ℂ)
    (i : Fin n) : ‖B i i‖ ≤ ‖B‖ := by
  classical
  set x : EuclideanSpace ℂ (Fin n) := EuclideanSpace.single i (1 : ℂ) with hxdef
  set y : EuclideanSpace ℂ (Fin n) :=
    (EuclideanSpace.equiv (Fin n) ℂ).symm (B *ᵥ x) with hydef
  have hxnorm : ‖x‖ = 1 := by
    rw [hxdef, EuclideanSpace.norm_single, norm_one]
  have hyi : y i = B i i := by
    show (B *ᵥ (Pi.single i (1 : ℂ) : Fin n → ℂ)) i = B i i
    simp [Matrix.mulVec, dotProduct_single]
  have hinner : ⟪x, y⟫ = B i i := by
    rw [hxdef, EuclideanSpace.inner_single_left, map_one, one_mul]
    exact hyi
  have hbound : ‖y‖ ≤ ‖B‖ * ‖x‖ := by
    rw [hydef]
    exact B.l2_opNorm_mulVec x
  calc ‖B i i‖ = ‖(⟪x, y⟫ : ℂ)‖ := by rw [hinner]
    _ ≤ ‖x‖ * ‖y‖ := norm_inner_le_norm x y
    _ ≤ ‖x‖ * (‖B‖ * ‖x‖) :=
        mul_le_mul_of_nonneg_left hbound (norm_nonneg x)
    _ = ‖B‖ := by rw [hxnorm, one_mul, mul_one]

/-- **Algebra-only.** A member of the unitary group has `L2` operator norm
one, by the C*-identity. -/
theorem norm_eq_one_of_mem_unitaryGroup [NeZero n]
    {V : Matrix (Fin n) (Fin n) ℂ}
    (hV : V ∈ Matrix.unitaryGroup (Fin n) ℂ) : ‖V‖ = 1 := by
  haveI : Nontrivial (Matrix (Fin n) (Fin n) ℂ) := by
    refine ⟨0, 1, fun h => ?_⟩
    have h00 := Matrix.ext_iff.mpr h ⟨0, Nat.pos_of_ne_zero (NeZero.ne n)⟩
      ⟨0, Nat.pos_of_ne_zero (NeZero.ne n)⟩
    rw [Matrix.zero_apply, Matrix.one_apply_eq] at h00
    exact zero_ne_one h00
  have h1 : star V * V = 1 := Matrix.mem_unitaryGroup_iff'.mp hV
  have hsq : ‖V‖ * ‖V‖ = 1 := by
    rw [← CStarRing.norm_star_mul_self, h1, norm_one]
  rcases mul_self_eq_one_iff.mp hsq with h | hneg
  · exact h
  · exfalso
    have := norm_nonneg V
    rw [hneg] at this
    norm_num at this

/-- **Trace-dependent.** The **state-expectation bound**: the expectation
of any observable under a state is dominated in modulus by the `L2`
operator norm of the observable. -/
theorem norm_expectation_le_l2_opNorm [NeZero n]
    {ρ : Matrix (Fin n) (Fin n) ℂ} (hρ : IsState ρ)
    (M : Matrix (Fin n) (Fin n) ℂ) :
    ‖expectation ρ M‖ ≤ ‖M‖ := by
  classical
  set V : Matrix (Fin n) (Fin n) ℂ :=
    (hρ.1.1.eigenvectorUnitary : Matrix (Fin n) (Fin n) ℂ) with hV
  set d : Fin n → ℂ := RCLike.ofReal ∘ hρ.1.1.eigenvalues with hd
  have hspec : ρ = V * diagonal d * star V := by
    rw [hV, hd]
    conv_lhs => rw [hρ.1.1.spectral_theorem, Unitary.conjStarAlgAut_apply]
  set B : Matrix (Fin n) (Fin n) ℂ := star V * M * V with hB
  have hBnorm : ‖B‖ ≤ ‖M‖ := by
    have hVnorm : ‖V‖ = 1 := by
      refine norm_eq_one_of_mem_unitaryGroup ?_
      rw [hV]
      exact (hρ.1.1.eigenvectorUnitary).2
    have hVstar : ‖star V‖ = 1 := by rw [norm_star, hVnorm]
    calc ‖B‖ = ‖star V * M * V‖ := by rw [hB]
      _ ≤ ‖star V * M‖ * ‖V‖ := l2_opNorm_mul _ _
      _ ≤ ‖star V‖ * ‖M‖ * ‖V‖ :=
          mul_le_mul_of_nonneg_right (l2_opNorm_mul _ _) (norm_nonneg _)
      _ = ‖M‖ := by rw [hVstar, hVnorm, one_mul, mul_one]
  have htr : expectation ρ M = ∑ i, d i * B i i := by
    change bornWeight ρ M = _
    rw [bornWeight, hspec,
      show V * diagonal d * star V * M = V * (diagonal d * (star V * M)) by
        simp only [mul_assoc],
      trace_mul_comm,
      show diagonal d * (star V * M) * V = diagonal d * B by
        rw [hB]; simp only [mul_assoc]]
    simp [Matrix.trace, Matrix.diag, Matrix.diagonal_mul]
  have hsum1 : ∑ i, hρ.1.1.eigenvalues i = 1 := by
    have h := hρ.1.1.trace_eq_sum_eigenvalues
    have h3 : ((∑ i, hρ.1.1.eigenvalues i : ℝ) : ℂ) = ((1 : ℝ) : ℂ) := by
      push_cast
      exact h.symm.trans hρ.2
    exact Complex.ofReal_injective h3
  calc ‖expectation ρ M‖ = ‖∑ i, d i * B i i‖ := by rw [htr]
    _ ≤ ∑ i, ‖d i * B i i‖ := norm_sum_le _ _
    _ ≤ ∑ i, hρ.1.1.eigenvalues i * ‖M‖ := by
        refine Finset.sum_le_sum fun i _ => ?_
        rw [norm_mul]
        have hdi : ‖d i‖ = hρ.1.1.eigenvalues i := by
          rw [hd]
          simp only [Function.comp_apply]
          rw [RCLike.norm_ofReal]
          exact abs_of_nonneg (hρ.1.eigenvalues_nonneg i)
        rw [hdi]
        exact mul_le_mul_of_nonneg_left
          ((norm_diag_entry_le_l2_opNorm B i).trans hBnorm)
          (hρ.1.eigenvalues_nonneg i)
    _ = ‖M‖ := by rw [← Finset.sum_mul, hsum1, one_mul]

/-- **Trace-dependent.** The **state-level CHSH corollary**: for a state
and four projection events with cross-party commutation, the expectation
of the CHSH combination of the dichotomic observables `1 - 2P` has
modulus at most `2√2`. This composes the state-expectation bound with the
operator-norm Tsirelson bound; no attainment claim is made. -/
theorem matrix_state_tsirelson_bound_of_events [NeZero n]
    {ρ : Matrix (Fin n) (Fin n) ℂ} (hρ : IsState ρ)
    (A₀ A₁ B₀ B₁ : Matrix (Fin n) (Fin n) ℂ)
    (hA₀ : IsEvent A₀) (hA₁ : IsEvent A₁)
    (hB₀ : IsEvent B₀) (hB₁ : IsEvent B₁)
    (h₀₀ : A₀ * B₀ = B₀ * A₀) (h₀₁ : A₀ * B₁ = B₁ * A₀)
    (h₁₀ : A₁ * B₀ = B₀ * A₁) (h₁₁ : A₁ * B₁ = B₁ * A₁) :
    ‖expectation ρ
        ((1 - 2 • A₀) * (1 - 2 • B₀) + (1 - 2 • A₀) * (1 - 2 • B₁) +
          (1 - 2 • A₁) * (1 - 2 • B₀) - (1 - 2 • A₁) * (1 - 2 • B₁))‖ ≤
      2 * Real.sqrt 2 :=
  (norm_expectation_le_l2_opNorm hρ _).trans
    (matrix_tsirelson_bound_of_events A₀ A₁ B₀ B₁ hA₀ hA₁ hB₀ hB₁
      h₀₀ h₀₁ h₁₀ h₁₁)

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms norm_diag_entry_le_l2_opNorm
#print axioms norm_eq_one_of_mem_unitaryGroup
#print axioms norm_expectation_le_l2_opNorm
#print axioms matrix_state_tsirelson_bound_of_events

end EventAlgebra
