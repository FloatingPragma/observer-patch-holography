import Mathlib.Analysis.CStarAlgebra.Matrix
import Mathlib.Analysis.SpecialFunctions.Exponential
import Mathlib.Analysis.Calculus.Deriv.Star
import Mathlib.Analysis.Calculus.Deriv.Shift
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.Normed.Module.FiniteDimension
import Mathlib.Data.Matrix.Basis
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Tactic.Module

/-!
# Finite Stone converse on one full matrix block

`PrivateInner.lean` records the forward private direction: a supplied
self-adjoint Hamiltonian generates a unitary real-parameter von Neumann flow,
and every star automorphism of one full finite matrix block is unitarily
inner.  This module proves the converse direction on `Matrix ι ι ℂ` under an
explicit differentiability-at-zero hypothesis.

1. Every complex-linear Leibniz endomorphism of the block is an inner
   derivation, with an explicit matrix-unit witness
   (`matrixDerivation_eq_commutator`).  The proof is the finite matrix-unit
   computation and uses no operator topology.
2. Every star-compatible Leibniz endomorphism is the von Neumann commutator
   action of one self-adjoint generator
   (`starDerivation_exists_selfAdjoint_generator`).
3. That generator is unique up to one additive real scalar multiple of the
   identity (`selfAdjoint_generator_unique_up_to_real_scalar`,
   `stoneFlow_generator_unique_up_to_real_scalar`).
4. A real-parameter star-automorphism group of the block whose orbit maps
   are differentiable at parameter zero is conjugation by the coherent
   unitary family `stonePropagator Ham t = exp (t • ((-I) • Ham))` of one
   time-independent self-adjoint `Ham`
   (`finiteStoneConverse_of_hasDerivAt`).  The family satisfies
   `stonePropagator Ham 0 = 1`, the additive-to-multiplicative group law,
   and unitarity of every member.
5. The differentiability premise is discharged from pointwise continuity
   alone (`continuous_starAutomorphismGroup_exists_hasDerivAt`,
   `finiteStoneConverse_of_continuous`).  The group acts on the matrix-unit
   basis as a norm-continuous matrix family; averaging that family over a
   short parameter interval gives an invertible element, and the
   fundamental theorem of calculus recovers the family from its primitive
   as a differentiable expression.

The results concern one full matrix block.  `CentralBlocks.lean`,
`CenterSpectral.lean`, and `StarWedderburn.lean` extend the construction to
arbitrary finite-dimensional private central-block decompositions.  No
physical clock, source-selected Hamiltonian, or physical time identification
is claimed; the parameter is a real number.
-/

namespace OPH.Dynamics

universe u

variable {ι : Type u} [Fintype ι] [DecidableEq ι]

/-! ## Inner derivations of a finite full matrix block -/

section InnerDerivation

variable (δ : Matrix ι ι ℂ →ₗ[ℂ] Matrix ι ι ℂ)

/-- Matrix-unit witness implementing a Leibniz endomorphism as an inner
derivation.  The reference column `i₀` is arbitrary. -/
noncomputable def innerDerivationWitness (i₀ : ι) : Matrix ι ι ℂ :=
  ∑ i : ι, δ (Matrix.single i i₀ 1) * Matrix.single i₀ i 1

private theorem witness_mul_single (i₀ k l : ι) :
    innerDerivationWitness δ i₀ * Matrix.single k l 1 =
      δ (Matrix.single k i₀ 1) * Matrix.single i₀ l 1 := by
  rw [innerDerivationWitness, Finset.sum_mul,
    Finset.sum_eq_single k
      (fun b _ hbk => by
        rw [mul_assoc, Matrix.single_mul_single_of_ne (h := hbk), mul_zero])
      (fun hk => absurd (Finset.mem_univ k) hk),
    mul_assoc, Matrix.single_mul_single_same, one_mul]

private theorem single_mul_witness
    (hδ : ∀ x y, δ (x * y) = δ x * y + x * δ y) (i₀ k l : ι) :
    Matrix.single k l 1 * innerDerivationWitness δ i₀ =
      δ (Matrix.single k i₀ 1) * Matrix.single i₀ l 1 -
        δ (Matrix.single k l 1) := by
  have hterm : ∀ i : ι,
      Matrix.single k l (1 : ℂ) *
          (δ (Matrix.single i i₀ 1) * Matrix.single i₀ i 1) =
        δ (Matrix.single k l 1 * Matrix.single i i₀ 1) * Matrix.single i₀ i 1
          - δ (Matrix.single k l 1) *
              (Matrix.single i i₀ 1 * Matrix.single i₀ i 1) := by
    intro i
    rw [hδ (Matrix.single k l 1) (Matrix.single i i₀ 1)]
    noncomm_ring
  rw [innerDerivationWitness, Finset.mul_sum,
    Finset.sum_congr rfl fun i _ => hterm i, Finset.sum_sub_distrib]
  congr 1
  · rw [Finset.sum_eq_single l
      (fun b _ hbl => by
        rw [Matrix.single_mul_single_of_ne (h := Ne.symm hbl), map_zero,
          zero_mul])
      (fun hl => absurd (Finset.mem_univ l) hl),
      Matrix.single_mul_single_same, one_mul]
  · rw [← Finset.mul_sum]
    have hone : (∑ i : ι, Matrix.single i i₀ (1 : ℂ) * Matrix.single i₀ i 1)
        = 1 := by
      calc ∑ i : ι, Matrix.single i i₀ (1 : ℂ) * Matrix.single i₀ i 1
          = ∑ i : ι, Matrix.single i i (1 : ℂ) := by
            refine Finset.sum_congr rfl fun i _ => ?_
            rw [Matrix.single_mul_single_same, one_mul]
        _ = 1 := Matrix.sum_single_one
    rw [hone, mul_one]

/-- Every complex-linear Leibniz endomorphism of a finite full matrix block
is the commutator with its matrix-unit witness: derivations of
`Matrix ι ι ℂ` are inner. -/
theorem matrixDerivation_eq_commutator
    (hδ : ∀ x y, δ (x * y) = δ x * y + x * δ y) (i₀ : ι)
    (x : Matrix ι ι ℂ) :
    δ x = innerDerivationWitness δ i₀ * x -
      x * innerDerivationWitness δ i₀ := by
  have hbasis : ∀ k l : ι,
      δ (Matrix.single k l 1) =
        innerDerivationWitness δ i₀ * Matrix.single k l 1 -
          Matrix.single k l 1 * innerDerivationWitness δ i₀ := by
    intro k l
    rw [witness_mul_single, single_mul_witness δ hδ, sub_sub_cancel]
  have hmaps : δ = LinearMap.mulLeft ℂ (innerDerivationWitness δ i₀) -
      LinearMap.mulRight ℂ (innerDerivationWitness δ i₀) := by
    refine Matrix.ext_linearMap ℂ fun k l => LinearMap.ext fun c => ?_
    show δ (Matrix.single k l c) =
      innerDerivationWitness δ i₀ * Matrix.single k l c -
        Matrix.single k l c * innerDerivationWitness δ i₀
    have hc : Matrix.single k l c = c • Matrix.single k l (1 : ℂ) := by
      rw [Matrix.smul_single, smul_eq_mul, mul_one]
    rw [hc, map_smul, hbasis k l, smul_sub, mul_smul_comm, smul_mul_assoc]
  have hx := LinearMap.congr_fun hmaps x
  simpa using hx

end InnerDerivation

/-! ## Self-adjoint generator from a star derivation -/

section StarDerivation

private theorem scalar_eq_smul_one (c : ℂ) :
    Matrix.scalar ι c = c • (1 : Matrix ι ι ℂ) := by
  rw [Matrix.scalar_apply, Matrix.smul_one_eq_diagonal]

variable (δ : Matrix ι ι ℂ →ₗ[ℂ] Matrix ι ι ℂ)

/-- A star-compatible Leibniz endomorphism of a finite full matrix block is
the von Neumann commutator action `x ↦ (-i) • (Ham * x - x * Ham)` of one
self-adjoint generator. -/
theorem starDerivation_exists_selfAdjoint_generator [Nonempty ι]
    (hδ : ∀ x y, δ (x * y) = δ x * y + x * δ y)
    (hstar : ∀ x, δ (star x) = star (δ x)) :
    ∃ Ham : Matrix ι ι ℂ, IsSelfAdjoint Ham ∧
      ∀ x, δ x = (-Complex.I) • (Ham * x - x * Ham) := by
  classical
  obtain ⟨i₀⟩ := ‹Nonempty ι›
  set W := innerDerivationWitness δ i₀ with hWdef
  have hcomm : ∀ x, δ x = W * x - x * W :=
    matrixDerivation_eq_commutator δ hδ i₀
  have hcentral : ∀ z : Matrix ι ι ℂ,
      (W + star W) * z = z * (W + star W) := by
    intro z
    have h1 : δ (star (star z)) = star (δ (star z)) := hstar (star z)
    rw [star_star, hcomm (star z), hcomm z, star_sub, star_mul, star_mul,
      star_star] at h1
    have h4 := sub_eq_sub_iff_add_eq_add.mp h1
    rw [add_mul, mul_add]
    calc W * z + star W * z = z * star W + z * W := h4
      _ = z * W + z * star W := add_comm _ _
  obtain ⟨c, hc⟩ := Matrix.mem_range_scalar_of_commute_single
    (M := W + star W) (fun i j _ => (hcentral (Matrix.single i j 1)).symm)
  have hself : star (W + star W) = W + star W := by
    rw [star_add, star_star, add_comm]
  have hcentry : c = (W + star W) i₀ i₀ := by
    rw [← hc, Matrix.scalar_apply, Matrix.diagonal_apply_eq]
  have hcreal : star c = c := by
    rw [hcentry, ← Matrix.star_apply, hself]
  have hchalf : star (c / 2) = c / 2 := by
    rw [star_div₀, hcreal, star_ofNat]
  set K : Matrix ι ι ℂ := W - (c / 2) • 1 with hKdef
  have hstarW : star W = c • (1 : Matrix ι ι ℂ) - W := by
    have hsum : c • (1 : Matrix ι ι ℂ) = W + star W := by
      rw [← scalar_eq_smul_one, hc]
    rw [hsum]
    module
  have hKskew : star K = -K := by
    have h1 : star K = star W - (c / 2) • (1 : Matrix ι ι ℂ) := by
      rw [hKdef, star_sub, star_smul, star_one, hchalf]
    rw [h1, hstarW, hKdef]
    module
  refine ⟨Complex.I • K, ?_, ?_⟩
  · have h1 : star (Complex.I • K) = (-Complex.I) • (-K) := by
      rw [star_smul, hKskew, Complex.star_def, Complex.conj_I]
    have h2 : ((-Complex.I) • (-K) : Matrix ι ι ℂ) = Complex.I • K := by
      module
    exact h1.trans h2
  · intro x
    have hzx : ((c / 2) • (1 : Matrix ι ι ℂ)) * x =
        x * ((c / 2) • (1 : Matrix ι ι ℂ)) := by
      rw [smul_mul_assoc, one_mul, mul_smul_comm, mul_one]
    have hKcomm : K * x - x * K = W * x - x * W := by
      rw [hKdef, sub_mul, mul_sub, hzx]
      abel
    have h1 : ((Complex.I • K) * x : Matrix ι ι ℂ) = Complex.I • (K * x) :=
      smul_mul_assoc _ _ _
    have h2 : (x * (Complex.I • K) : Matrix ι ι ℂ) = Complex.I • (x * K) :=
      mul_smul_comm _ _ _
    rw [hcomm x, ← hKcomm, h1, h2]
    match_scalars <;> simp

/-- Two self-adjoint generators with one commutator action differ by one
real scalar multiple of the identity. -/
theorem selfAdjoint_generator_unique_up_to_real_scalar [Nonempty ι]
    {H₁ H₂ : Matrix ι ι ℂ} (h₁ : IsSelfAdjoint H₁) (h₂ : IsSelfAdjoint H₂)
    (hsame : ∀ x, (-Complex.I) • (H₁ * x - x * H₁) =
      (-Complex.I) • (H₂ * x - x * H₂)) :
    ∃ r : ℝ, H₁ = H₂ + (r : ℂ) • 1 := by
  classical
  obtain ⟨i₀⟩ := ‹Nonempty ι›
  have hcomm : ∀ x : Matrix ι ι ℂ, (H₁ - H₂) * x = x * (H₁ - H₂) := by
    intro x
    have hcanc : ∀ y : Matrix ι ι ℂ,
        Complex.I • ((-Complex.I) • y) = y := fun y => by
      match_scalars
      simp
    have h'' := congrArg (fun z : Matrix ι ι ℂ => Complex.I • z) (hsame x)
    have h' : H₁ * x - x * H₁ = H₂ * x - x * H₂ := by
      simpa only [hcanc] using h''
    have h4 := sub_eq_sub_iff_add_eq_add.mp h'
    rw [sub_mul, mul_sub, sub_eq_sub_iff_add_eq_add]
    calc H₁ * x + x * H₂ = H₂ * x + x * H₁ := h4
      _ = x * H₁ + H₂ * x := add_comm _ _
  obtain ⟨d, hd⟩ := Matrix.mem_range_scalar_of_commute_single
    (M := H₁ - H₂) (fun i j _ => (hcomm (Matrix.single i j 1)).symm)
  have hself : star (H₁ - H₂) = H₁ - H₂ := by
    rw [star_sub, h₁.star_eq, h₂.star_eq]
  have hdentry : d = (H₁ - H₂) i₀ i₀ := by
    rw [← hd, Matrix.scalar_apply, Matrix.diagonal_apply_eq]
  have hdreal : star d = d := by
    rw [hdentry, ← Matrix.star_apply, hself]
  have hdre : (d.re : ℂ) = d := Complex.conj_eq_iff_re.mp hdreal
  refine ⟨d.re, ?_⟩
  have hsub : H₁ - H₂ = (d.re : ℂ) • 1 := by
    rw [hdre, ← scalar_eq_smul_one, hd]
  rw [← hsub]
  module

end StarDerivation

/-! ## The coherent unitary family and the Stone converse -/

section Stone

open scoped Matrix.Norms.L2Operator

/-- The unitary family generated by a self-adjoint matrix, at the von
Neumann sign convention of `PrivateInner.lean`. -/
noncomputable def stonePropagator (Ham : Matrix ι ι ℂ) (t : ℝ) :
    Matrix ι ι ℂ :=
  NormedSpace.exp (t • ((-Complex.I) • Ham))

theorem stonePropagator_zero (Ham : Matrix ι ι ℂ) :
    stonePropagator Ham 0 = 1 := by
  have hz : ((0 : ℝ) • ((-Complex.I) • Ham) : Matrix ι ι ℂ) = 0 :=
    zero_smul ℝ _
  rw [stonePropagator, hz, NormedSpace.exp_zero]

/-- The propagator family is a one-parameter multiplicative group. -/
theorem stonePropagator_add (Ham : Matrix ι ι ℂ) (s t : ℝ) :
    stonePropagator Ham (s + t) =
      stonePropagator Ham s * stonePropagator Ham t := by
  let +nondep : NormedAlgebra ℚ (Matrix ι ι ℂ) :=
    .restrictScalars ℚ ℂ (Matrix ι ι ℂ)
  have hadd : ((s + t) • ((-Complex.I) • Ham) : Matrix ι ι ℂ) =
      s • ((-Complex.I) • Ham) + t • ((-Complex.I) • Ham) := add_smul s t _
  have hcomm : Commute ((s • ((-Complex.I) • Ham) : Matrix ι ι ℂ))
      (t • ((-Complex.I) • Ham)) :=
    ((Commute.refl _).smul_left s).smul_right t
  rw [stonePropagator, stonePropagator, stonePropagator, hadd,
    NormedSpace.exp_add_of_commute hcomm]

/-- Every member of the propagator family of a self-adjoint generator is
unitary. -/
theorem stonePropagator_mem_unitary
    (Ham : Matrix ι ι ℂ) (hham : IsSelfAdjoint Ham) (t : ℝ) :
    stonePropagator Ham t ∈ unitary (Matrix ι ι ℂ) := by
  let +nondep : NormedAlgebra ℚ (Matrix ι ι ℂ) :=
    .restrictScalars ℚ ℂ (Matrix ι ι ℂ)
  rw [stonePropagator]
  apply NormedSpace.exp_mem_unitary_of_mem_skewAdjoint
  apply skewAdjoint.smul_mem t
  apply hham.smul_mem_skewAdjoint
  rw [skewAdjoint.mem_iff]
  simp

/-- Derivative of the propagator conjugation at parameter zero: the von
Neumann commutator. -/
theorem hasDerivAt_stoneFlow_zero (Ham x : Matrix ι ι ℂ) :
    HasDerivAt
      (fun s : ℝ => stonePropagator Ham s * x * stonePropagator Ham (-s))
      ((-Complex.I) • (Ham * x - x * Ham)) 0 := by
  have hneg : ∀ s : ℝ,
      ((-s : ℝ) • ((-Complex.I) • Ham) : Matrix ι ι ℂ) =
        s • (-((-Complex.I) • Ham)) := fun s =>
    (neg_smul s _).trans (smul_neg s _).symm
  have hfun :
      (fun s : ℝ => stonePropagator Ham s * x * stonePropagator Ham (-s)) =
        fun s : ℝ => NormedSpace.exp (s • ((-Complex.I) • Ham)) * x *
          NormedSpace.exp (s • (-((-Complex.I) • Ham))) := by
    funext s
    rw [stonePropagator, stonePropagator, hneg s]
  rw [hfun]
  have h1 : HasDerivAt
      (fun u : ℝ => NormedSpace.exp (u • ((-Complex.I) • Ham)))
      (NormedSpace.exp ((0 : ℝ) • ((-Complex.I) • Ham)) *
        ((-Complex.I) • Ham)) 0 :=
    hasDerivAt_exp_smul_const _ 0
  have h3 : HasDerivAt
      (fun u : ℝ => NormedSpace.exp (u • (-((-Complex.I) • Ham))))
      ((-((-Complex.I) • Ham)) *
        NormedSpace.exp ((0 : ℝ) • (-((-Complex.I) • Ham)))) 0 :=
    hasDerivAt_exp_smul_const' _ 0
  have h123 : HasDerivAt
      (fun s : ℝ => NormedSpace.exp (s • ((-Complex.I) • Ham)) * x *
        NormedSpace.exp (s • (-((-Complex.I) • Ham))))
      (NormedSpace.exp ((0 : ℝ) • ((-Complex.I) • Ham)) *
          ((-Complex.I) • Ham) * x *
          NormedSpace.exp ((0 : ℝ) • (-((-Complex.I) • Ham))) +
        NormedSpace.exp ((0 : ℝ) • ((-Complex.I) • Ham)) * x *
          ((-((-Complex.I) • Ham)) *
            NormedSpace.exp ((0 : ℝ) • (-((-Complex.I) • Ham))))) 0 :=
    (h1.mul_const x).mul h3
  have hz1 : ((0 : ℝ) • ((-Complex.I) • Ham) : Matrix ι ι ℂ) = 0 :=
    zero_smul ℝ _
  have hz2 : ((0 : ℝ) • (-((-Complex.I) • Ham)) : Matrix ι ι ℂ) = 0 :=
    zero_smul ℝ _
  have hval1 : NormedSpace.exp ((0 : ℝ) • ((-Complex.I) • Ham)) *
          ((-Complex.I) • Ham) * x *
          NormedSpace.exp ((0 : ℝ) • (-((-Complex.I) • Ham))) +
        NormedSpace.exp ((0 : ℝ) • ((-Complex.I) • Ham)) * x *
          ((-((-Complex.I) • Ham)) *
            NormedSpace.exp ((0 : ℝ) • (-((-Complex.I) • Ham)))) =
      ((-Complex.I) • Ham) * x - x * ((-Complex.I) • Ham) := by
    rw [hz1, hz2, NormedSpace.exp_zero]
    noncomm_ring
  have hval2 : (((-Complex.I) • Ham) * x - x * ((-Complex.I) • Ham) :
        Matrix ι ι ℂ) = (-Complex.I) • (Ham * x - x * Ham) := by
    rw [smul_mul_assoc, mul_smul_comm]
    module
  rw [hval1.trans hval2] at h123
  exact h123

/-- Finite Stone converse on one full matrix block.  A real-parameter group
of star automorphisms whose orbit maps are differentiable at parameter zero
is conjugation by the coherent unitary family of one time-independent
self-adjoint Hamiltonian. -/
theorem finiteStoneConverse_of_hasDerivAt [Nonempty ι]
    (α : ℝ → (Matrix ι ι ℂ ≃⋆ₐ[ℂ] Matrix ι ι ℂ))
    (hzero : ∀ x, α 0 x = x)
    (hgroup : ∀ s t : ℝ, ∀ x, α (s + t) x = α s (α t x))
    (hderiv : ∀ x, ∃ dx : Matrix ι ι ℂ,
      HasDerivAt (fun t : ℝ => α t x) dx 0) :
    ∃ Ham : Matrix ι ι ℂ, IsSelfAdjoint Ham ∧
      ∀ (t : ℝ) (x : Matrix ι ι ℂ),
        α t x = stonePropagator Ham t * x * stonePropagator Ham (-t) := by
  classical
  choose d hd using hderiv
  have hadd : ∀ x y, d (x + y) = d x + d y := by
    intro x y
    have h : HasDerivAt (fun t : ℝ => α t x + α t y) (d x + d y) 0 :=
      (hd x).add (hd y)
    have hfun : (fun t : ℝ => α t x + α t y) =
        fun t : ℝ => α t (x + y) := by
      funext s
      rw [map_add]
    rw [hfun] at h
    exact (hd (x + y)).unique h
  have hsmulc : ∀ (a : ℂ) x, d (a • x) = a • d x := by
    intro a x
    have h : HasDerivAt (fun t : ℝ => a • α t x) (a • d x) 0 :=
      (hd x).const_smul a
    have hfun : (fun t : ℝ => a • α t x) = fun t : ℝ => α t (a • x) := by
      funext s
      rw [map_smul]
    rw [hfun] at h
    exact (hd (a • x)).unique h
  let δ : Matrix ι ι ℂ →ₗ[ℂ] Matrix ι ι ℂ :=
    { toFun := d
      map_add' := hadd
      map_smul' := hsmulc }
  have hleib : ∀ x y, d (x * y) = d x * y + x * d y := by
    intro x y
    have h : HasDerivAt (fun t : ℝ => α t x * α t y)
        (d x * α 0 y + α 0 x * d y) 0 := (hd x).mul (hd y)
    rw [hzero x, hzero y] at h
    have hfun : (fun t : ℝ => α t x * α t y) =
        fun t : ℝ => α t (x * y) := by
      funext s
      rw [map_mul]
    rw [hfun] at h
    exact (hd (x * y)).unique h
  have hstar : ∀ x, d (star x) = star (d x) := by
    intro x
    have h : HasDerivAt (fun t : ℝ => star (α t x)) (star (d x)) 0 :=
      (hd x).star
    have hfun : (fun t : ℝ => star (α t x)) =
        fun t : ℝ => α t (star x) := by
      funext s
      rw [map_star]
    rw [hfun] at h
    exact (hd (star x)).unique h
  obtain ⟨Ham, hham, hHam⟩ :=
    starDerivation_exists_selfAdjoint_generator δ hleib hstar
  have hdK : ∀ y : Matrix ι ι ℂ,
      d y = ((-Complex.I) • Ham) * y - y * ((-Complex.I) • Ham) := by
    intro y
    have h : d y = (-Complex.I) • (Ham * y - y * Ham) := hHam y
    have h1 : (((-Complex.I) • Ham) * y : Matrix ι ι ℂ) =
        (-Complex.I) • (Ham * y) := smul_mul_assoc _ _ _
    have h2 : (y * ((-Complex.I) • Ham) : Matrix ι ι ℂ) =
        (-Complex.I) • (y * Ham) := mul_smul_comm _ _ _
    rw [h, h1, h2]
    module
  have hglobal : ∀ (t : ℝ) (x : Matrix ι ι ℂ),
      HasDerivAt (fun s : ℝ => α s x) (d (α t x)) t := by
    intro t x
    have h0 : HasDerivAt (fun h : ℝ => α h (α t x)) (d (α t x))
        (t + -t) := by
      rw [add_neg_cancel]
      exact hd (α t x)
    have h1 := h0.comp_add_const
    have hfun : (fun y : ℝ => α (y + -t) (α t x)) =
        fun y : ℝ => α y x := by
      funext y
      rw [← hgroup (y + -t) t x, add_assoc, neg_add_cancel, add_zero]
    rw [hfun] at h1
    exact h1
  refine ⟨Ham, hham, ?_⟩
  intro t x
  set K : Matrix ι ι ℂ := (-Complex.I) • Ham with hKdef
  set g : ℝ → Matrix ι ι ℂ := fun s =>
    NormedSpace.exp (s • (-K)) * α s x * NormedSpace.exp (s • K) with hgdef
  have hgderiv : ∀ s : ℝ, HasDerivAt g 0 s := by
    intro s
    have h1 : HasDerivAt (fun u : ℝ => NormedSpace.exp (u • (-K)))
        (NormedSpace.exp (s • (-K)) * (-K)) s :=
      hasDerivAt_exp_smul_const (-K) s
    have h2 : HasDerivAt (fun u : ℝ => α u x) (d (α s x)) s := hglobal s x
    have h3 : HasDerivAt (fun u : ℝ => NormedSpace.exp (u • K))
        (K * NormedSpace.exp (s • K)) s := hasDerivAt_exp_smul_const' K s
    have h123 : HasDerivAt g
        ((NormedSpace.exp (s • (-K)) * (-K) * α s x +
            NormedSpace.exp (s • (-K)) * d (α s x)) *
          NormedSpace.exp (s • K) +
            NormedSpace.exp (s • (-K)) * α s x *
              (K * NormedSpace.exp (s • K))) s := (h1.mul h2).mul h3
    have hval : (NormedSpace.exp (s • (-K)) * (-K) * α s x +
          NormedSpace.exp (s • (-K)) * d (α s x)) *
        NormedSpace.exp (s • K) +
          NormedSpace.exp (s • (-K)) * α s x *
            (K * NormedSpace.exp (s • K)) = 0 := by
      rw [hdK (α s x)]
      noncomm_ring
    rw [hval] at h123
    exact h123
  have hconst : g t = g 0 :=
    is_const_of_deriv_eq_zero (fun s => (hgderiv s).differentiableAt)
      (fun s => (hgderiv s).deriv) t 0
  have hz1 : ((0 : ℝ) • (-K) : Matrix ι ι ℂ) = 0 := zero_smul ℝ _
  have hz2 : ((0 : ℝ) • K : Matrix ι ι ℂ) = 0 := zero_smul ℝ _
  have hg0 : g 0 = x := by
    show NormedSpace.exp ((0 : ℝ) • (-K)) * α 0 x *
      NormedSpace.exp ((0 : ℝ) • K) = x
    rw [hz1, hz2, NormedSpace.exp_zero, one_mul, mul_one, hzero]
  have hgt : NormedSpace.exp (t • (-K)) * α t x * NormedSpace.exp (t • K)
      = x := hconst.trans hg0
  have hcomm : Commute ((t • K : Matrix ι ι ℂ)) (t • (-K)) :=
    ((Commute.refl K).neg_right.smul_left t).smul_right t
  have hsm : (t • (-K) : Matrix ι ι ℂ) = -(t • K) := smul_neg t K
  have hsum : (t • K + t • (-K) : Matrix ι ι ℂ) = 0 := by
    rw [hsm]
    exact add_neg_cancel _
  have hcancel : NormedSpace.exp (t • K) * NormedSpace.exp (t • (-K))
      = 1 := by
    let +nondep : NormedAlgebra ℚ (Matrix ι ι ℂ) :=
      .restrictScalars ℚ ℂ (Matrix ι ι ℂ)
    rw [← NormedSpace.exp_add_of_commute hcomm, hsum, NormedSpace.exp_zero]
  have hconj : α t x =
      NormedSpace.exp (t • K) * x * NormedSpace.exp (t • (-K)) := by
    calc α t x
        = NormedSpace.exp (t • K) * NormedSpace.exp (t • (-K)) * α t x *
            (NormedSpace.exp (t • K) * NormedSpace.exp (t • (-K))) := by
          rw [hcancel, one_mul, mul_one]
      _ = NormedSpace.exp (t • K) *
            (NormedSpace.exp (t • (-K)) * α t x *
              NormedSpace.exp (t • K)) *
            NormedSpace.exp (t • (-K)) := by noncomm_ring
      _ = NormedSpace.exp (t • K) * x * NormedSpace.exp (t • (-K)) := by
          rw [hgt]
  have hneg : ((-t : ℝ) • K : Matrix ι ι ℂ) = t • (-K) :=
    (neg_smul t K).trans (smul_neg t K).symm
  rw [stonePropagator, stonePropagator, ← hKdef, hneg]
  exact hconj

/-- Finite expansion of a matrix over the product-indexed matrix units. -/
private theorem matrix_eq_sum_single_prod (x : Matrix ι ι ℂ) :
    x = ∑ q : ι × ι, x q.1 q.2 • Matrix.single q.1 q.2 (1 : ℂ) := by
  have h1 : (∑ q : ι × ι, x q.1 q.2 • Matrix.single q.1 q.2 (1 : ℂ)) =
      ∑ i : ι, ∑ j : ι, x i j • Matrix.single i j (1 : ℂ) := by
    rw [← Finset.univ_product_univ, Finset.sum_product]
  rw [h1]
  calc x = ∑ i : ι, ∑ j : ι, Matrix.single i j (x i j) :=
      Matrix.matrix_eq_sum_single x
    _ = ∑ i : ι, ∑ j : ι, x i j • Matrix.single i j (1 : ℂ) := by
      refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
      rw [Matrix.smul_single, smul_eq_mul, mul_one]

/-- The product-indexed entry matrix of one group element: the matrix of its
linear action in the matrix-unit basis. -/
private noncomputable def groupMatrix
    (α : ℝ → (Matrix ι ι ℂ ≃⋆ₐ[ℂ] Matrix ι ι ℂ)) (t : ℝ) :
    Matrix (ι × ι) (ι × ι) ℂ :=
  Matrix.of fun p q => α t (Matrix.single q.1 q.2 1) p.1 p.2

/-- Differentiability upgrade by interval averaging: every pointwise
continuous real-parameter star-automorphism group of the block has
differentiable orbit maps at parameter zero.  The group is packaged as a
norm-continuous matrix family in the matrix-unit basis, averaged over a
short parameter interval into an invertible element, and recovered from its
primitive through the fundamental theorem of calculus. -/
theorem continuous_starAutomorphismGroup_exists_hasDerivAt
    (α : ℝ → (Matrix ι ι ℂ ≃⋆ₐ[ℂ] Matrix ι ι ℂ))
    (hzero : ∀ x, α 0 x = x)
    (hgroup : ∀ s t : ℝ, ∀ x, α (s + t) x = α s (α t x))
    (hcont : ∀ x, Continuous fun t : ℝ => α t x) :
    ∀ x : Matrix ι ι ℂ, ∃ dx : Matrix ι ι ℂ,
      HasDerivAt (fun t : ℝ => α t x) dx 0 := by
  classical
  haveI : IsScalarTower ℝ ℂ (Matrix ι ι ℂ) :=
    ⟨fun r c N => by
      ext i j
      exact smul_assoc r c (N i j)⟩
  haveI : IsScalarTower ℝ ℂ (Matrix (ι × ι) (ι × ι) ℂ) :=
    ⟨fun r c N => by
      ext p q
      exact smul_assoc r c (N p q)⟩
  set A : ℝ → Matrix (ι × ι) (ι × ι) ℂ := groupMatrix α with hAdef
  have hexp : ∀ (t : ℝ) (x : Matrix ι ι ℂ),
      α t x = ∑ q : ι × ι, x q.1 q.2 • α t (Matrix.single q.1 q.2 1) := by
    intro t x
    conv_lhs => rw [matrix_eq_sum_single_prod x]
    rw [map_sum]
    refine Finset.sum_congr rfl fun q _ => ?_
    rw [map_smul]
  have hentry : ∀ (t : ℝ) (x : Matrix ι ι ℂ) (i j : ι),
      α t x i j = ∑ q : ι × ι, A t (i, j) q * x q.1 q.2 := by
    intro t x i j
    rw [hexp t x, Matrix.sum_apply]
    refine Finset.sum_congr rfl fun q _ => ?_
    rw [Matrix.smul_apply, smul_eq_mul, mul_comm]
    rfl
  have hA : Continuous A := by
    apply continuous_matrix
    intro p q
    exact (hcont (Matrix.single q.1 q.2 1)).matrix_elem p.1 p.2
  have hA0 : A 0 = 1 := by
    ext p q
    show α 0 (Matrix.single q.1 q.2 1) p.1 p.2 =
      (1 : Matrix (ι × ι) (ι × ι) ℂ) p q
    rw [hzero, Matrix.single_apply, Matrix.one_apply]
    by_cases hpq : p = q
    · subst hpq
      simp
    · have hne : ¬(q.1 = p.1 ∧ q.2 = p.2) := fun h =>
        hpq (Prod.ext h.1.symm h.2.symm)
      rw [if_neg hne, if_neg hpq]
  have hAgrp : ∀ s t : ℝ, A (s + t) = A s * A t := by
    intro s t
    ext p q
    show α (s + t) (Matrix.single q.1 q.2 1) p.1 p.2 = _
    rw [hgroup s t, hentry s (α t (Matrix.single q.1 q.2 1)) p.1 p.2,
      Matrix.mul_apply]
    rfl
  have hφ : Continuous fun s : ℝ => ‖A s - 1‖ :=
    (hA.sub continuous_const).norm
  have hφ0 : ‖A 0 - 1‖ = 0 := by rw [hA0, sub_self, norm_zero]
  obtain ⟨δ, hδpos, hδ⟩ := Metric.continuousAt_iff.mp hφ.continuousAt
    (1 / 2) (by norm_num)
  set ε : ℝ := δ / 2 with hεdef
  have hεpos : 0 < ε := by positivity
  have hbound : ∀ s ∈ Set.uIoc (0 : ℝ) ε, ‖A s - 1‖ ≤ 1 / 2 := by
    intro s hs
    rw [Set.uIoc_of_le hεpos.le] at hs
    have hdist : dist s 0 < δ := by
      rw [Real.dist_eq, sub_zero, abs_of_pos hs.1]
      calc s ≤ ε := hs.2
        _ < δ := by rw [hεdef]; linarith
    have hnear := hδ hdist
    rw [Real.dist_eq, hφ0, sub_zero] at hnear
    exact (le_abs_self _).trans hnear.le
  have hAint : ∀ a b : ℝ,
      IntervalIntegrable A MeasureTheory.volume a b := fun a b =>
    hA.intervalIntegrable a b
  set V : Matrix (ι × ι) (ι × ι) ℂ := ∫ s in (0:ℝ)..ε, A s with hVdef
  have hsub : (∫ s in (0:ℝ)..ε, (A s - 1)) =
      V - ε • (1 : Matrix (ι × ι) (ι × ι) ℂ) := by
    rw [intervalIntegral.integral_sub (hAint 0 ε) intervalIntegrable_const,
      intervalIntegral.integral_const, sub_zero, hVdef]
    congr 1
  have hVnorm : ‖V - ε • (1 : Matrix (ι × ι) (ι × ι) ℂ)‖ ≤
      1 / 2 * |ε - 0| := by
    rw [← hsub]
    exact intervalIntegral.norm_integral_le_of_norm_le_const hbound
  have hWnorm : ‖(1 : Matrix (ι × ι) (ι × ι) ℂ) - ε⁻¹ • V‖ < 1 := by
    have h1 : (1 : Matrix (ι × ι) (ι × ι) ℂ) - ε⁻¹ • V =
        -(ε⁻¹ • (V - ε • (1 : Matrix (ι × ι) (ι × ι) ℂ))) := by
      have ha : (ε⁻¹ • (V - ε • (1 : Matrix (ι × ι) (ι × ι) ℂ)) :
            Matrix (ι × ι) (ι × ι) ℂ) =
          ε⁻¹ • V - ε⁻¹ • (ε • (1 : Matrix (ι × ι) (ι × ι) ℂ)) :=
        smul_sub _ _ _
      have hb : (ε⁻¹ • (ε • (1 : Matrix (ι × ι) (ι × ι) ℂ)) :
            Matrix (ι × ι) (ι × ι) ℂ) =
          (ε⁻¹ * ε) • (1 : Matrix (ι × ι) (ι × ι) ℂ) := smul_smul _ _ _
      have hc' : (ε⁻¹ * ε : ℝ) = 1 := inv_mul_cancel₀ hεpos.ne'
      have hd : ((1 : ℝ) • (1 : Matrix (ι × ι) (ι × ι) ℂ) :
          Matrix (ι × ι) (ι × ι) ℂ) = 1 := one_smul ℝ _
      rw [ha, hb, hc', hd]
      exact (neg_sub _ _).symm
    have h3 : ‖(ε⁻¹ : ℝ)‖ = ε⁻¹ := by
      rw [Real.norm_eq_abs, abs_of_pos (by positivity)]
    have h5 : |ε - 0| = ε := by rw [sub_zero, abs_of_pos hεpos]
    rw [h1, norm_neg]
    calc ‖(ε⁻¹ • (V - ε • (1 : Matrix (ι × ι) (ι × ι) ℂ)) :
            Matrix (ι × ι) (ι × ι) ℂ)‖
        ≤ ‖(ε⁻¹ : ℝ)‖ *
            ‖V - ε • (1 : Matrix (ι × ι) (ι × ι) ℂ)‖ := norm_smul_le _ _
      _ = ε⁻¹ * ‖V - ε • (1 : Matrix (ι × ι) (ι × ι) ℂ)‖ := by rw [h3]
      _ ≤ ε⁻¹ * (1 / 2 * |ε - 0|) :=
          mul_le_mul_of_nonneg_left hVnorm (by positivity)
      _ = 1 / 2 := by
          rw [h5, mul_comm (1 / 2 : ℝ) ε, ← mul_assoc,
            inv_mul_cancel₀ hεpos.ne', one_mul]
      _ < 1 := by norm_num
  have hWunit : IsUnit ((ε⁻¹ • V : Matrix (ι × ι) (ι × ι) ℂ)) :=
    sub_sub_self (1 : Matrix (ι × ι) (ι × ι) ℂ) (ε⁻¹ • V) ▸
      (Units.oneSub ((1 : Matrix (ι × ι) (ι × ι) ℂ) - ε⁻¹ • V)
        hWnorm).isUnit
  obtain ⟨u, hu⟩ := hWunit
  set w : Matrix (ι × ι) (ι × ι) ℂ := ↑u⁻¹ with hwdef
  have huw : (↑u : Matrix (ι × ι) (ι × ι) ℂ) * w = 1 := u.mul_inv
  have hshift : ∀ h : ℝ, A h * V = ∫ s in h..(h + ε), A s := by
    intro h
    have h1 : A h * V = ∫ s in (0:ℝ)..ε, A h * A s := by
      have hcc := (ContinuousLinearMap.mul ℂ
          (Matrix (ι × ι) (ι × ι) ℂ) (A h)).intervalIntegral_comp_comm
        (hAint 0 ε)
      exact hcc.symm
    have h2 : (∫ s in (0:ℝ)..ε, A h * A s) = ∫ s in (0:ℝ)..ε, A (h + s) :=
      intervalIntegral.integral_congr fun s _ => (hAgrp h s).symm
    have h3 : (∫ s in (0:ℝ)..ε, A (h + s)) = ∫ s in h..(h + ε), A s := by
      have hca := intervalIntegral.integral_comp_add_left (f := A)
        (a := (0:ℝ)) (b := ε) h
      rw [add_zero] at hca
      exact hca
    rw [h1, h2, h3]
  have hFTC : ∀ b : ℝ,
      HasDerivAt (fun v : ℝ => ∫ s in (0:ℝ)..v, A s) (A b) b := fun b =>
    (hA.integral_hasStrictDerivAt 0 b).hasDerivAt
  have hArep : ∀ h : ℝ, A h =
      (ε⁻¹ • ((∫ s in (0:ℝ)..(h + ε), A s) - ∫ s in (0:ℝ)..h, A s)) *
        w := by
    intro h
    have hsu : (∫ s in (0:ℝ)..(h + ε), A s) - (∫ s in (0:ℝ)..h, A s) =
        ∫ s in h..(h + ε), A s :=
      intervalIntegral.integral_interval_sub_left (hAint 0 (h + ε))
        (hAint 0 h)
    have hmulu : A h * (↑u : Matrix (ι × ι) (ι × ι) ℂ) =
        ε⁻¹ • ((∫ s in (0:ℝ)..(h + ε), A s) - ∫ s in (0:ℝ)..h, A s) := by
      rw [hu]
      have h1 : A h * (ε⁻¹ • V) = ε⁻¹ • (A h * V) := mul_smul_comm _ _ _
      rw [h1, hshift h, hsu]
    calc A h = A h * (↑u : Matrix (ι × ι) (ι × ι) ℂ) * w := by
          rw [mul_assoc, huw, mul_one]
      _ = _ := by rw [hmulu]
  have hd1 : HasDerivAt (fun h : ℝ => ∫ s in (0:ℝ)..(h + ε), A s)
      (A (0 + ε)) 0 := by
    have h0 : HasDerivAt (fun v : ℝ => ∫ s in (0:ℝ)..v, A s) (A (0 + ε))
        (0 + ε) := hFTC (0 + ε)
    exact h0.comp_add_const
  have hd2 : HasDerivAt (fun h : ℝ => ∫ s in (0:ℝ)..h, A s) (A 0) 0 :=
    hFTC 0
  have hd3 : HasDerivAt (fun h : ℝ =>
      (ε⁻¹ • ((∫ s in (0:ℝ)..(h + ε), A s) - ∫ s in (0:ℝ)..h, A s)) * w)
      ((ε⁻¹ • (A (0 + ε) - A 0)) * w) 0 :=
    ((hd1.sub hd2).const_smul ε⁻¹).mul_const _
  have hAd : HasDerivAt A ((ε⁻¹ • (A (0 + ε) - A 0)) * w) 0 := by
    have hfun : (fun h : ℝ =>
        (ε⁻¹ • ((∫ s in (0:ℝ)..(h + ε), A s) - ∫ s in (0:ℝ)..h, A s)) *
          w) = A := funext fun h => (hArep h).symm
    rw [hfun] at hd3
    exact hd3
  intro x
  let Φ : Matrix (ι × ι) (ι × ι) ℂ →ₗ[ℝ] Matrix ι ι ℂ :=
    { toFun := fun M =>
        Matrix.of fun i j => ∑ q : ι × ι, M (i, j) q * x q.1 q.2
      map_add' := fun M N => by
        ext i j
        show (∑ q : ι × ι, (M (i, j) q + N (i, j) q) * x q.1 q.2) =
          (∑ q : ι × ι, M (i, j) q * x q.1 q.2) +
            ∑ q : ι × ι, N (i, j) q * x q.1 q.2
        rw [← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun q _ => ?_
        rw [add_mul]
      map_smul' := fun r M => by
        ext i j
        show (∑ q : ι × ι, (r • M (i, j) q) * x q.1 q.2) =
          r • ∑ q : ι × ι, M (i, j) q * x q.1 q.2
        have hs : (r • ∑ q : ι × ι, M (i, j) q * x q.1 q.2 : ℂ) =
            ∑ q : ι × ι, r • (M (i, j) q * x q.1 q.2) := Finset.smul_sum
        rw [hs]
        refine Finset.sum_congr rfl fun q _ => ?_
        exact smul_mul_assoc r (M (i, j) q) (x q.1 q.2) }
  have hΦc : Continuous Φ := by
    show Continuous fun M : Matrix (ι × ι) (ι × ι) ℂ =>
      Matrix.of fun i j => ∑ q : ι × ι, M (i, j) q * x q.1 q.2
    apply continuous_matrix
    intro i j
    apply continuous_finset_sum
    intro q _
    exact (continuous_id.matrix_elem (i, j) q).mul continuous_const
  have hΦd : HasDerivAt (fun t : ℝ => Φ (A t))
      (Φ ((ε⁻¹ • (A (0 + ε) - A 0)) * w)) 0 :=
    (ContinuousLinearMap.mk Φ hΦc).hasFDerivAt.comp_hasDerivAt 0 hAd
  have hfx : (fun t : ℝ => α t x) = fun t : ℝ => Φ (A t) := by
    funext t
    ext i j
    exact hentry t x i j
  refine ⟨Φ ((ε⁻¹ • (A (0 + ε) - A 0)) * w), ?_⟩
  rw [hfx]
  exact hΦd

/-- Finite Stone converse from pointwise continuity alone: a pointwise
continuous real-parameter star-automorphism group of one full finite matrix
block is conjugation by the coherent unitary family of one time-independent
self-adjoint Hamiltonian. -/
theorem finiteStoneConverse_of_continuous [Nonempty ι]
    (α : ℝ → (Matrix ι ι ℂ ≃⋆ₐ[ℂ] Matrix ι ι ℂ))
    (hzero : ∀ x, α 0 x = x)
    (hgroup : ∀ s t : ℝ, ∀ x, α (s + t) x = α s (α t x))
    (hcont : ∀ x, Continuous fun t : ℝ => α t x) :
    ∃ Ham : Matrix ι ι ℂ, IsSelfAdjoint Ham ∧
      ∀ (t : ℝ) (x : Matrix ι ι ℂ),
        α t x = stonePropagator Ham t * x * stonePropagator Ham (-t) :=
  finiteStoneConverse_of_hasDerivAt α hzero hgroup
    (continuous_starAutomorphismGroup_exists_hasDerivAt α hzero hgroup hcont)

/-- Flow-level uniqueness: two self-adjoint generators with one propagator
conjugation differ by one additive real scalar multiple of the identity. -/
theorem stoneFlow_generator_unique_up_to_real_scalar [Nonempty ι]
    {H₁ H₂ : Matrix ι ι ℂ} (h₁ : IsSelfAdjoint H₁) (h₂ : IsSelfAdjoint H₂)
    (hflow : ∀ (t : ℝ) (x : Matrix ι ι ℂ),
      stonePropagator H₁ t * x * stonePropagator H₁ (-t) =
        stonePropagator H₂ t * x * stonePropagator H₂ (-t)) :
    ∃ r : ℝ, H₁ = H₂ + (r : ℂ) • 1 := by
  refine selfAdjoint_generator_unique_up_to_real_scalar h₁ h₂ fun x => ?_
  have hfx : (fun s : ℝ =>
      stonePropagator H₁ s * x * stonePropagator H₁ (-s)) =
        fun s : ℝ => stonePropagator H₂ s * x * stonePropagator H₂ (-s) :=
    funext fun s => hflow s x
  have h2 : HasDerivAt
      (fun s : ℝ => stonePropagator H₁ s * x * stonePropagator H₁ (-s))
      ((-Complex.I) • (H₂ * x - x * H₂)) 0 := by
    rw [hfx]
    exact hasDerivAt_stoneFlow_zero H₂ x
  exact (hasDerivAt_stoneFlow_zero H₁ x).unique h2

end Stone

#print axioms OPH.Dynamics.matrixDerivation_eq_commutator
#print axioms OPH.Dynamics.starDerivation_exists_selfAdjoint_generator
#print axioms OPH.Dynamics.selfAdjoint_generator_unique_up_to_real_scalar
#print axioms OPH.Dynamics.stonePropagator_add
#print axioms OPH.Dynamics.stonePropagator_mem_unitary
#print axioms OPH.Dynamics.finiteStoneConverse_of_hasDerivAt
#print axioms OPH.Dynamics.continuous_starAutomorphismGroup_exists_hasDerivAt
#print axioms OPH.Dynamics.finiteStoneConverse_of_continuous
#print axioms OPH.Dynamics.stoneFlow_generator_unique_up_to_real_scalar
#print axioms OPH.Dynamics.stonePropagator_zero
#print axioms OPH.Dynamics.hasDerivAt_stoneFlow_zero

end OPH.Dynamics
