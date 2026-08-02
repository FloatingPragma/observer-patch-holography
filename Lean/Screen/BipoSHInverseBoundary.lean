import Mathlib
import ObserverPatchHolography.ScalarSeamRepair
import BipoSHTransferInvariant

namespace OPH.BipoSHInverseBoundary

/-!
# The inverse-continuum and operational-response boundary

This file records three exact facts used by the BipoSH continuum gate.

First, a complete uniform seam repair determines its graph stiffness from a
one-tick response:

`L = 2 |E| (I - R)`.

Second, positive finite stiffness operators can converge while their inverses
diverge if no uniform lower spectral bound is supplied.  The two-coordinate
family with eigenvalues `1` and `1/(n+1)` is the minimal counterexample.

Third, scalar transfer cancellation does not extend to an unresolved copy or
radial multiplicity space.  A two-copy readout changes the normalized cross
statistic from one to one half while commuting with the spatial action on the
irreducible factor.

The results use ordinary theorem arguments.  They do not select the refined
equal-seam law, a stochastic source ensemble, a physical intervention, a
screen-to-sky transfer, multiplicity one, or a continuum coercivity margin.
-/

noncomputable section

open ObserverPatchHolography.ScalarSeamRepair

variable {ι ε : Type*} [DecidableEq ι] [Fintype ε] [Nonempty ε]

/-- The exact operational readback identity.  Once the complete finite move
alphabet is uniformly selected, the expected one-tick repair `R` determines
the stiffness without introducing a covariance or an inverse operator. -/
theorem stiffness_recovered_from_uniform_repair
    (left right : ε → ι) :
    graphLaplacian left right =
      (2 * (Fintype.card ε : ℝ)) •
        (LinearMap.id - uniformSeamRepair left right) := by
  rw [uniformSeamRepair_eq_id_sub_graphLaplacian]
  have hcard : (Fintype.card ε : ℝ) ≠ 0 := by
    exact_mod_cast Fintype.card_ne_zero
  have hcoefficient :
      (2 * (Fintype.card ε : ℝ)) *
          (1 / (2 * (Fintype.card ε : ℝ))) = 1 := by
    field_simp
  symm
  calc
    (2 * (Fintype.card ε : ℝ)) •
          (LinearMap.id -
            (LinearMap.id -
              (1 / (2 * (Fintype.card ε : ℝ))) •
                graphLaplacian left right)) =
        (2 * (Fintype.card ε : ℝ)) •
          ((1 / (2 * (Fintype.card ε : ℝ))) •
            graphLaplacian left right) := by module
    _ = ((2 * (Fintype.card ε : ℝ)) *
          (1 / (2 * (Fintype.card ε : ℝ)))) •
            graphLaplacian left right := by rw [smul_smul]
    _ = graphLaplacian left right := by rw [hcoefficient, one_smul]

/-! ## Uniform coercivity is necessary before inversion -/

/-- The collapsing nuisance eigenvalue in the exact two-coordinate control. -/
def collapsingNuisance (n : ℕ) : ℝ := 1 / (n + 1 : ℝ)

/-- Its formal inverse eigenvalue. -/
def inverseNuisance (n : ℕ) : ℝ := n + 1

/-- Every finite member of the counterexample is strictly positive. -/
theorem collapsing_stiffness_positive (n : ℕ) :
    0 < collapsingNuisance n := by
  unfold collapsingNuisance
  positivity

/-- The two declared eigenvalues are exact reciprocals at every finite rung. -/
theorem collapsing_inverse_exact (n : ℕ) :
    collapsingNuisance n * inverseNuisance n = 1 := by
  unfold collapsingNuisance inverseNuisance
  field_simp

/-- The stiffness nuisance eigenvalue can be made smaller than any positive
threshold.  Thus finite positivity supplies no uniform coercivity margin. -/
theorem collapsing_stiffness_arbitrarily_small
    {δ : ℝ} (hδ : 0 < δ) :
    ∃ n : ℕ, collapsingNuisance n < δ := by
  obtain ⟨n, hn⟩ := exists_nat_gt (1 / δ)
  refine ⟨n, ?_⟩
  unfold collapsingNuisance
  have hden : 0 < (n + 1 : ℝ) := by positivity
  apply (div_lt_iff₀ hden).2
  have hn' : 1 / δ < (n + 1 : ℝ) := by
    exact hn.trans (by norm_num)
  have hproduct : 1 < (n + 1 : ℝ) * δ :=
    (div_lt_iff₀ hδ).1 hn'
  simpa [mul_comm] using hproduct

/-- The inverse nuisance eigenvalue is unbounded. -/
theorem inverse_nuisance_unbounded (M : ℝ) :
    ∃ n : ℕ, M < inverseNuisance n := by
  obtain ⟨n, hn⟩ := exists_nat_gt M
  refine ⟨n, ?_⟩
  have hsucc : (n : ℝ) < (n : ℝ) + 1 := by norm_num
  simpa [inverseNuisance] using hn.trans hsucc

/-! ## Copy-space mixing survives rotation equivariance -/

/-- The normalized readout for `A=C=I₂` and `B=diag(1,0)`.  The spatial
irreducible factor is suppressed; `u` and `v` act only on two-dimensional
copy spaces and therefore commute with every spatial rotation. -/
def copyReadoutStatistic (u₁ u₂ v₁ v₂ : ℝ) : ℝ :=
  |u₁ * v₁| /
    Real.sqrt ((u₁ ^ 2 + u₂ ^ 2) * (v₁ ^ 2 + v₂ ^ 2))

@[simp]
theorem copy_readout_first :
    copyReadoutStatistic 1 0 1 0 = 1 := by
  norm_num [copyReadoutStatistic]

@[simp]
theorem copy_readout_mixed :
    copyReadoutStatistic 1 1 1 1 = 1 / 2 := by
  norm_num [copyReadoutStatistic]

/-- Two rotation-compatible copy-space readouts give different normalized
statistics.  Schur scalar cancellation therefore needs multiplicity one (or
a separately selected copy-space map). -/
theorem copy_mixing_changes_statistic :
    copyReadoutStatistic 1 0 1 0 ≠
      copyReadoutStatistic 1 1 1 1 := by
  norm_num

/-! ## Arithmetic admission radius -/

/-- The standard scalar coefficient in the resolvent perturbation bound. -/
def inverseTailRadius (gamma epsilon : ℝ) : ℝ :=
  epsilon / (gamma * (gamma - epsilon))

/-- The resolvent radius is positive only inside a positive coercivity
margin.  The matrix resolvent identity remains a separate operator theorem
input; this lemma verifies the scalar admission arithmetic. -/
theorem inverseTailRadius_pos
    {gamma epsilon : ℝ}
    (hgamma : 0 < gamma) (hepsilon : 0 < epsilon)
    (hgate : epsilon < gamma) :
    0 < inverseTailRadius gamma epsilon := by
  unfold inverseTailRadius
  exact div_pos hepsilon (mul_pos hgamma (sub_pos.mpr hgate))

end

end OPH.BipoSHInverseBoundary

#check OPH.BipoSHInverseBoundary.stiffness_recovered_from_uniform_repair
#check OPH.BipoSHInverseBoundary.collapsing_stiffness_positive
#check OPH.BipoSHInverseBoundary.collapsing_inverse_exact
#check OPH.BipoSHInverseBoundary.collapsing_stiffness_arbitrarily_small
#check OPH.BipoSHInverseBoundary.inverse_nuisance_unbounded
#check OPH.BipoSHInverseBoundary.copy_readout_first
#check OPH.BipoSHInverseBoundary.copy_readout_mixed
#check OPH.BipoSHInverseBoundary.copy_mixing_changes_statistic
#check OPH.BipoSHInverseBoundary.inverseTailRadius_pos

#print axioms OPH.BipoSHInverseBoundary.stiffness_recovered_from_uniform_repair
#print axioms OPH.BipoSHInverseBoundary.collapsing_stiffness_arbitrarily_small
#print axioms OPH.BipoSHInverseBoundary.inverse_nuisance_unbounded
#print axioms OPH.BipoSHInverseBoundary.copy_mixing_changes_statistic
