import Variational.DiscreteEulerLagrange

namespace OPH.Variational

/-!
# The discrete Noether theorem for invariant local actions

A one-parameter family of transformations that leaves the two-point
Lagrangian invariant yields equality of the two segment momenta adjacent to
each supplied junction that satisfies the discrete Euler-Lagrange equation.
The local chain is finite and exact: differentiating the
invariance in the flow parameter gives the infinitesimal identity
`∂₁L·ξ + ∂₂L·ξ = 0`, and combining it with the Euler-Lagrange
equation at that junction transports the momentum
`p(n) = ∂₂L(γ(n), γ(n+1)) · ξ(γ(n+1))` unchanged across the junction.
When the same minimization hypothesis holds at every interior record, the
chain-wide theorem below applies that local transport at every adjacent pair
and proves that all segment momenta equal one scalar `J` on the whole finite
chain.

The flow, its generator, the Lagrangian derivative data, and the
invariance are named inputs. No physical symmetry, charge, or
continuum Noether current is constructed here; the physical
identification of the resulting scalar chain constant with a laboratory
current is a separate obligation.
-/

variable {N : ℕ}

/-- **Infinitesimal invariance.** Differentiating the exact invariance
`L (T s x) (T s y) = L x y` in the flow parameter at `s = 0` gives
`∂₁L(x,y) · ξ(x) + ∂₂L(x,y) · ξ(y) = 0`, with `ξ` the family
generator. -/
theorem infinitesimal_invariance
    (L : ℝ → ℝ → ℝ) (T : ℝ → ℝ → ℝ) (ξ : ℝ → ℝ)
    (d₁ d₂ : ℝ → ℝ → ℝ)
    (hT0 : ∀ x, T 0 x = x)
    (hTd : ∀ x, HasDerivAt (fun s => T s x) (ξ x) 0)
    (hL : ∀ p : ℝ × ℝ, HasFDerivAt (Function.uncurry L)
      ((d₁ p.1 p.2) • (ContinuousLinearMap.fst ℝ ℝ ℝ)
        + (d₂ p.1 p.2) • (ContinuousLinearMap.snd ℝ ℝ ℝ)) p)
    (hinv : ∀ s x y, L (T s x) (T s y) = L x y)
    (x y : ℝ) :
    d₁ x y * ξ x + d₂ x y * ξ y = 0 := by
  have hg : HasDerivAt (fun s => (T s x, T s y))
      ((ξ x, ξ y) : ℝ × ℝ) 0 := (hTd x).prodMk (hTd y)
  have hLpt : HasFDerivAt (Function.uncurry L)
      ((d₁ x y) • (ContinuousLinearMap.fst ℝ ℝ ℝ)
        + (d₂ x y) • (ContinuousLinearMap.snd ℝ ℝ ℝ))
      ((T 0 x, T 0 y) : ℝ × ℝ) := by
    rw [hT0 x, hT0 y]
    exact hL (x, y)
  have hcomp0 := hLpt.comp_hasDerivAt 0 hg
  have hcomp : HasDerivAt (fun s => L (T s x) (T s y))
      (((d₁ x y) • (ContinuousLinearMap.fst ℝ ℝ ℝ)
        + (d₂ x y) • (ContinuousLinearMap.snd ℝ ℝ ℝ))
          ((ξ x, ξ y) : ℝ × ℝ)) 0 := by
    simpa [Function.comp, Function.uncurry] using hcomp0
  have happ : ((d₁ x y) • (ContinuousLinearMap.fst ℝ ℝ ℝ)
      + (d₂ x y) • (ContinuousLinearMap.snd ℝ ℝ ℝ))
        ((ξ x, ξ y) : ℝ × ℝ)
      = d₁ x y * ξ x + d₂ x y * ξ y := by
    simp [ContinuousLinearMap.add_apply,
      ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst',
      ContinuousLinearMap.coe_snd', smul_eq_mul]
  have hconstfun : (fun s => L (T s x) (T s y))
      = fun _ => L x y := funext fun s => hinv s x y
  have hzero : HasDerivAt (fun s => L (T s x) (T s y)) 0 0 := by
    rw [hconstfun]
    exact hasDerivAt_const 0 (L x y)
  have := hcomp.unique hzero
  rw [happ] at this
  linarith [this]

/-- The segment momentum of a history: the second-slot Lagrangian
derivative on a segment, paired with the flow generator at the
segment's later record. -/
noncomputable def segmentMomentum (d₂ : ℝ → ℝ → ℝ) (ξ : ℝ → ℝ)
    (γ : Fin (N + 1) → ℝ) (n : Fin N) : ℝ :=
  d₂ (γ n.castSucc) (γ n.succ) * ξ (γ n.succ)

/-- **Discrete Noether step.** At a junction where the discrete
Euler-Lagrange equation holds and the flow leaves the Lagrangian
infinitesimally invariant, the segment momentum is transported
unchanged across the junction. -/
theorem noether_step (d₁ d₂ : ℝ → ℝ → ℝ) (ξ : ℝ → ℝ)
    (γ : Fin (N + 1) → ℝ) {k m : Fin N}
    (hkm : k.succ = m.castSucc)
    (hEL : d₂ (γ k.castSucc) (γ k.succ)
      + d₁ (γ m.castSucc) (γ m.succ) = 0)
    (hinf : d₁ (γ m.castSucc) (γ m.succ) * ξ (γ m.castSucc)
      + d₂ (γ m.castSucc) (γ m.succ) * ξ (γ m.succ) = 0) :
    segmentMomentum d₂ ξ γ k = segmentMomentum d₂ ξ γ m := by
  unfold segmentMomentum
  have hpt : γ k.succ = γ m.castSucc := by rw [hkm]
  rw [hpt] at hEL ⊢
  linear_combination ξ (γ m.castSucc) * hEL - hinf

/-- **Discrete Noether theorem, composed form.** A history that
minimizes the local action among single-site variations at a junction,
for a Lagrangian invariant under a differentiable one-parameter family,
transports the segment momentum unchanged across that junction. -/
theorem noether_conserved
    (L : ℝ → ℝ → ℝ) (T : ℝ → ℝ → ℝ) (ξ : ℝ → ℝ)
    (d₁ d₂ : ℝ → ℝ → ℝ)
    (hT0 : ∀ x, T 0 x = x)
    (hTd : ∀ x, HasDerivAt (fun s => T s x) (ξ x) 0)
    (hL : ∀ p : ℝ × ℝ, HasFDerivAt (Function.uncurry L)
      ((d₁ p.1 p.2) • (ContinuousLinearMap.fst ℝ ℝ ℝ)
        + (d₂ p.1 p.2) • (ContinuousLinearMap.snd ℝ ℝ ℝ)) p)
    (hinv : ∀ s x y, L (T s x) (T s y) = L x y)
    (γ : Fin (N + 1) → ℝ) {k m : Fin N}
    (hkm : k.succ = m.castSucc)
    (hmin : ∀ x, localAction L γ
      ≤ localAction L (Function.update γ k.succ x)) :
    segmentMomentum d₂ ξ γ k = segmentMomentum d₂ ξ γ m := by
  have hd₂ : HasDerivAt (fun x => L (γ k.castSucc) x)
      (d₂ (γ k.castSucc) (γ k.succ)) (γ k.succ) := by
    have hL' := hL (γ k.castSucc, γ k.succ)
    have hg : HasDerivAt (fun x => ((γ k.castSucc, x) : ℝ × ℝ))
        ((0, 1) : ℝ × ℝ) (γ k.succ) :=
      (hasDerivAt_const _ _).prodMk (hasDerivAt_id _)
    have := hL'.comp_hasDerivAt (γ k.succ) hg
    simpa [ContinuousLinearMap.add_apply,
      ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst',
      ContinuousLinearMap.coe_snd', smul_eq_mul,
      Function.uncurry] using this
  have hd₁ : HasDerivAt (fun x => L x (γ m.succ))
      (d₁ (γ m.castSucc) (γ m.succ)) (γ m.castSucc) := by
    have hL' := hL (γ m.castSucc, γ m.succ)
    have hg : HasDerivAt (fun x => ((x, γ m.succ) : ℝ × ℝ))
        ((1, 0) : ℝ × ℝ) (γ m.castSucc) :=
      (hasDerivAt_id _).prodMk (hasDerivAt_const _ _)
    have := hL'.comp_hasDerivAt (γ m.castSucc) hg
    simpa [ContinuousLinearMap.add_apply,
      ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst',
      ContinuousLinearMap.coe_snd', smul_eq_mul,
      Function.uncurry] using this
  have hEL := stationary_localAction_discreteEulerLagrange
    L γ hkm _ _ hd₂ hd₁ hmin
  have hinf := infinitesimal_invariance L T ξ d₁ d₂ hT0 hTd hL hinv
    (γ m.castSucc) (γ m.succ)
  exact noether_step d₁ d₂ ξ γ hkm hEL hinf

/-- **Chain-wide finite Noether constancy (scalar and conditional).**
Consider a path with `M + 2` records and hence `M + 1` consecutive segments.
If the invariant differentiable two-point Lagrangian is minimized under every
real single-site replacement at each of the `M` interior records, then there
is one scalar `J` equal to the segment momentum on every segment.

The proof does not assume pairwise momentum equality: it invokes
`noether_conserved` at every adjacent junction and then inducts along the
finite chain.  This remains a conditional scalar theorem; it supplies neither
a vector-valued current nor a physical charge, clock, or continuum limit. -/
theorem noether_current_constant_on_finite_chain
    {M : ℕ}
    (L : ℝ → ℝ → ℝ) (T : ℝ → ℝ → ℝ) (ξ : ℝ → ℝ)
    (d₁ d₂ : ℝ → ℝ → ℝ)
    (hT0 : ∀ x, T 0 x = x)
    (hTd : ∀ x, HasDerivAt (fun s => T s x) (ξ x) 0)
    (hL : ∀ p : ℝ × ℝ, HasFDerivAt (Function.uncurry L)
      ((d₁ p.1 p.2) • (ContinuousLinearMap.fst ℝ ℝ ℝ)
        + (d₂ p.1 p.2) • (ContinuousLinearMap.snd ℝ ℝ ℝ)) p)
    (hinv : ∀ s x y, L (T s x) (T s y) = L x y)
    (γ : Fin (M + 2) → ℝ)
    (hmin : ∀ i : Fin M, ∀ x : ℝ, localAction L γ
      ≤ localAction L (Function.update γ i.castSucc.succ x)) :
    ∃ J : ℝ, ∀ n : Fin (M + 1),
      segmentMomentum d₂ ξ γ n = J := by
  have hadj : ∀ i : Fin M,
      segmentMomentum d₂ ξ γ i.castSucc
        = segmentMomentum d₂ ξ γ i.succ := by
    intro i
    exact noether_conserved L T ξ d₁ d₂ hT0 hTd hL hinv γ
      (k := i.castSucc) (m := i.succ) rfl (hmin i)
  refine ⟨segmentMomentum d₂ ξ γ 0, ?_⟩
  intro n
  induction n using Fin.induction with
  | zero => rfl
  | succ i ih =>
      exact (hadj i).symm.trans ih

/-- Algebraic translation package for the free two-point expression: the
translation fixes every point at parameter zero and leaves `(y - x)^2`
unchanged.  This declaration does not itself contain the derivative packet. -/
theorem free_translation_invariance_package :
    (∀ x : ℝ, (fun s => x + s) 0 = x)
    ∧ (∀ s x y : ℝ, (y + s - (x + s))^2 = (y - x)^2) := by
  constructor
  · intro x
    simp
  · intro s x y
    ring

/-- The exact first-slot derivative used by the free translation witness. -/
def freeDOne (x y : ℝ) : ℝ := 2 * (x - y)

/-- The exact second-slot derivative used by the free translation witness. -/
def freeDTwo (x y : ℝ) : ℝ := 2 * (y - x)

/-- The unit translation generator. -/
def unitTranslationGenerator (_x : ℝ) : ℝ := 1

/-- A three-record affine path `0,1,2`. -/
def affineFreePath3 (i : Fin 3) : ℝ := i.val

/-- **Nonzero local Noether witness.** On the affine free path `0,1,2`, the
two adjacent free-particle segment momenta are equal by `noether_step`, and
their common value is nonzero.  This certifies that the local conservation
identity has nonzero instances.  It supplies no physical clock or
whole-history attachment. -/
theorem free_translation_nonzero_noether_witness :
    segmentMomentum freeDTwo unitTranslationGenerator affineFreePath3
        (0 : Fin 2)
      = segmentMomentum freeDTwo unitTranslationGenerator affineFreePath3
        (1 : Fin 2)
    ∧ segmentMomentum freeDTwo unitTranslationGenerator affineFreePath3
        (0 : Fin 2) ≠ 0 := by
  constructor
  · apply noether_step freeDOne freeDTwo unitTranslationGenerator
      affineFreePath3 (k := (0 : Fin 2)) (m := (1 : Fin 2))
    · decide
    · norm_num [freeDOne, freeDTwo, affineFreePath3]
    · norm_num [freeDOne, freeDTwo, unitTranslationGenerator,
        affineFreePath3]
  · norm_num [segmentMomentum, freeDTwo, unitTranslationGenerator,
      affineFreePath3]

end OPH.Variational

#print axioms OPH.Variational.infinitesimal_invariance
#print axioms OPH.Variational.noether_step
#print axioms OPH.Variational.noether_conserved
#print axioms OPH.Variational.noether_current_constant_on_finite_chain
#print axioms OPH.Variational.free_translation_nonzero_noether_witness
