import FiniteConditionalRepair
import CapFirstLaw

namespace OPH.Thermodynamics

/-!
# Fluctuation identities and reciprocity for conditional repair

The second law of the conditional-repair package is an inequality: the
relative entropy to the reference law is nonincreasing. This module
proves the exact identity underneath it. For one repair step from `p`
to `push p K`, the fluctuating entropy production

`sigmaEP π p q x y = log (p x / π x) - log (q y / π y)`

satisfies the integral fluctuation identity
`∑ p x K x y exp (-sigmaEP) = 1` for every kernel that preserves the
reference law, the pointwise detailed fluctuation relation
`p x K x y = exp (sigmaEP) q y K y x` for every kernel in detailed
balance with the reference law, and its level-set form: forward and
reversed step weights at fixed entropy production `s` differ by the
factor `exp s`. The mean of `sigmaEP` is exactly the certified
relative-entropy descent, so the second law is the Jensen shadow of an
exact identity rather than a primitive inequality. Detailed balance
also yields the finite Onsager relation: equilibrium time-lagged
correlations are symmetric under exchange of the two observables.

For the conditional-resampling kernel the descent itself carries an
exact Pythagorean split, `heatBath_kl_pythagorean`: relative entropy to
the reference decomposes as the relative entropy to the repaired image
plus the relative entropy from the repaired image to the reference. The
proof runs through the fibre-measurable log-ratio of one repair step
and the fibre-mean conservation of the repair kernel; it sums over the
state space only, so the visible alphabet needs no finiteness. The
identity is classical information theory for conditional-expectation
projections, restated for the committed kernel.

All statements are finite and elementary; the conditional-resampling
repair kernel satisfies every hypothesis through the theorems of
`FiniteConditionalRepair`.
-/

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]
variable {B : Type*} [DecidableEq B]

/-- Fluctuating entropy production of one kernel step from `p` (state
`x`) to `q` (state `y`), relative to the reference law `π`: the
difference of the two surprisals `log (p x / π x)` and
`log (q y / π y)`. -/
noncomputable def sigmaEP (pi p q : Ω → ℝ) (x y : Ω) : ℝ :=
  Real.log (p x / pi x) - Real.log (q y / pi y)

/-- **Integral fluctuation identity.** For every kernel with unit row
sums that preserves the positive reference law `π`, and every positive
normalized state `p` with positive image, the exponential of minus the
entropy production averages to one exactly:
`∑ p x K x y exp (-sigmaEP) = 1`. The second law is the Jensen
consequence of this identity. -/
theorem integral_fluctuation (pi p : Ω → ℝ) (K : Ω → Ω → ℝ)
    (hπ : ∀ x, 0 < pi x) (hp : ∀ x, 0 < p x) (hp1 : ∑ x, p x = 1)
    (hK1 : ∀ x, ∑ y, K x y = 1)
    (hstat : ∀ y, push pi K y = pi y)
    (hq : ∀ y, 0 < push p K y) :
    ∑ x, ∑ y, p x * K x y
      * Real.exp (-(sigmaEP pi p (push p K) x y)) = 1 := by
  have hterm : ∀ x y, p x * K x y
      * Real.exp (-(sigmaEP pi p (push p K) x y))
      = pi x * K x y * (push p K y / pi y) := by
    intro x y
    unfold sigmaEP
    rw [neg_sub, Real.exp_sub, Real.exp_log (div_pos (hq y) (hπ y)),
      Real.exp_log (div_pos (hp x) (hπ x))]
    have hx := (hp x).ne'
    have hπx := (hπ x).ne'
    have hπy := (hπ y).ne'
    field_simp
  calc ∑ x, ∑ y, p x * K x y
        * Real.exp (-(sigmaEP pi p (push p K) x y))
      = ∑ x, ∑ y, pi x * K x y * (push p K y / pi y) := by
        exact Finset.sum_congr rfl fun x _ =>
          Finset.sum_congr rfl fun y _ => hterm x y
    _ = ∑ y, ∑ x, pi x * K x y * (push p K y / pi y) :=
        Finset.sum_comm
    _ = ∑ y, push p K y / pi y * ∑ x, pi x * K x y := by
        exact Finset.sum_congr rfl fun y _ => by
          rw [Finset.mul_sum]
          exact Finset.sum_congr rfl fun x _ => by ring
    _ = ∑ y, push p K y / pi y * pi y := by
        exact Finset.sum_congr rfl fun y _ => by
          have hdef : ∑ x, pi x * K x y = push pi K y := rfl
          rw [hdef, hstat y]
    _ = ∑ y, push p K y := by
        exact Finset.sum_congr rfl fun y _ =>
          div_mul_cancel₀ _ (ne_of_gt (hπ y))
    _ = ∑ x, p x := push_total p K hK1
    _ = 1 := hp1

/-- **Pointwise detailed fluctuation relation (Crooks form).** For a
kernel in detailed balance with the positive reference law, the
forward step weight equals the reversed step weight amplified by the
exponential of the entropy production:
`p x K x y = exp (sigmaEP) (push p K) y K y x`, at every pair of
states. -/
theorem crooks_pointwise (pi p : Ω → ℝ) (K : Ω → Ω → ℝ)
    (hπ : ∀ x, 0 < pi x) (hp : ∀ x, 0 < p x)
    (hq : ∀ y, 0 < push p K y)
    (hdb : ∀ x y, pi x * K x y = pi y * K y x) (x y : Ω) :
    p x * K x y
      = Real.exp (sigmaEP pi p (push p K) x y)
        * (push p K y * K y x) := by
  unfold sigmaEP
  rw [Real.exp_sub, Real.exp_log (div_pos (hp x) (hπ x)),
    Real.exp_log (div_pos (hq y) (hπ y))]
  have hx := (hp x).ne'
  have hπx := (hπ x).ne'
  have hπy := (hπ y).ne'
  have hqy := (hq y).ne'
  field_simp
  linear_combination hdb x y

open scoped Classical in
/-- **Level-set fluctuation theorem.** Summed over the pairs with
entropy production exactly `s`, the forward weight equals `exp s`
times the reversed weight. -/
theorem crooks_level_set (pi p : Ω → ℝ) (K : Ω → Ω → ℝ)
    (hπ : ∀ x, 0 < pi x) (hp : ∀ x, 0 < p x)
    (hq : ∀ y, 0 < push p K y)
    (hdb : ∀ x y, pi x * K x y = pi y * K y x) (s : ℝ) :
    ∑ z ∈ Finset.univ.filter
        (fun z : Ω × Ω => sigmaEP pi p (push p K) z.1 z.2 = s),
      p z.1 * K z.1 z.2
      = Real.exp s * ∑ z ∈ Finset.univ.filter
        (fun z : Ω × Ω => sigmaEP pi p (push p K) z.1 z.2 = s),
      push p K z.2 * K z.2 z.1 := by
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro z hz
  have hs := (Finset.mem_filter.mp hz).2
  rw [crooks_pointwise pi p K hπ hp hq hdb z.1 z.2, hs]

/-- The mean entropy production of one step is exactly the certified
relative-entropy descent: `⟨sigmaEP⟩ = D(p ‖ π) - D(push p K ‖ π)`.
The transition-side second law is this identity combined with the
descent theorem. -/
theorem sigma_mean_eq_kl_descent (pi p : Ω → ℝ) (K : Ω → Ω → ℝ)
    (hp : ∀ x, 0 < p x)
    (hK1 : ∀ x, ∑ y, K x y = 1)
    (hq : ∀ y, 0 < push p K y) :
    ∑ x, ∑ y, p x * K x y * sigmaEP pi p (push p K) x y
      = kl p pi - kl (push p K) pi := by
  have hsplit : ∀ x y, p x * K x y * sigmaEP pi p (push p K) x y
      = p x * Real.log (p x / pi x) * K x y
        - p x * K x y * Real.log (push p K y / pi y) := by
    intro x y
    unfold sigmaEP
    ring
  calc ∑ x, ∑ y, p x * K x y * sigmaEP pi p (push p K) x y
      = ∑ x, ∑ y, (p x * Real.log (p x / pi x) * K x y
          - p x * K x y * Real.log (push p K y / pi y)) := by
        exact Finset.sum_congr rfl fun x _ =>
          Finset.sum_congr rfl fun y _ => hsplit x y
    _ = (∑ x, ∑ y, p x * Real.log (p x / pi x) * K x y)
        - ∑ x, ∑ y, p x * K x y
            * Real.log (push p K y / pi y) := by
        rw [← Finset.sum_sub_distrib]
        exact Finset.sum_congr rfl fun x _ => by
          rw [Finset.sum_sub_distrib]
    _ = (∑ x, p x * Real.log (p x / pi x))
        - ∑ y, push p K y * Real.log (push p K y / pi y) := by
        congr 1
        · exact Finset.sum_congr rfl fun x _ => by
            rw [← Finset.mul_sum, hK1 x, mul_one]
        · rw [Finset.sum_comm]
          exact Finset.sum_congr rfl fun y _ => by
            rw [← Finset.sum_mul]
            rfl
    _ = kl p pi - kl (push p K) pi := by
        unfold kl
        congr 1
        · exact (Finset.sum_congr rfl fun x _ =>
            (klTerm_of_pos (ne_of_gt (hp x)) (pi x))).symm
        · exact (Finset.sum_congr rfl fun y _ =>
            (klTerm_of_pos (ne_of_gt (hq y)) (pi y))).symm

omit [DecidableEq Ω] in
/-- **Finite Onsager reciprocity.** In detailed balance, the
equilibrium time-lagged correlation of two observables is symmetric:
`⟨f(X₀) g(X₁)⟩_π = ⟨g(X₀) f(X₁)⟩_π`. Linear-response transport
coefficients built from these correlations inherit the symmetry. -/
theorem correlation_symm (pi : Ω → ℝ) (K : Ω → Ω → ℝ)
    (hdb : ∀ x y, pi x * K x y = pi y * K y x) (f g : Ω → ℝ) :
    ∑ x, ∑ y, pi x * K x y * f x * g y
      = ∑ x, ∑ y, pi x * K x y * g x * f y := by
  calc ∑ x, ∑ y, pi x * K x y * f x * g y
      = ∑ x, ∑ y, pi y * K y x * f x * g y := by
        exact Finset.sum_congr rfl fun x _ =>
          Finset.sum_congr rfl fun y _ => by rw [hdb x y]
    _ = ∑ y, ∑ x, pi y * K y x * f x * g y := Finset.sum_comm
    _ = ∑ x, ∑ y, pi x * K x y * g x * f y := by
        exact Finset.sum_congr rfl fun y _ =>
          Finset.sum_congr rfl fun x _ => by ring

section HeatBathInstantiation

variable {π : Ω → ℝ} {b : Ω → B}

/-- The conditional-resampling repair kernel satisfies the integral
fluctuation identity. -/
theorem heatBath_integral_fluctuation (hπ : ∀ x, 0 < π x)
    (p : Ω → ℝ) (hp : ∀ x, 0 < p x) (hp1 : ∑ x, p x = 1) :
    ∑ x, ∑ y, p x * heatBath π b x y
      * Real.exp
        (-(sigmaEP π p (push p (heatBath π b)) x y)) = 1 :=
  integral_fluctuation π p (heatBath π b) hπ hp hp1
    (heatBath_row_sum hπ) (heatBath_stationary hπ)
    (heatBath_preserves_pos hπ p hp)

/-- The conditional-resampling repair kernel satisfies the pointwise
detailed fluctuation relation. -/
theorem heatBath_crooks (hπ : ∀ x, 0 < π x)
    (p : Ω → ℝ) (hp : ∀ x, 0 < p x) (x y : Ω) :
    p x * heatBath π b x y
      = Real.exp (sigmaEP π p (push p (heatBath π b)) x y)
        * (push p (heatBath π b) y * heatBath π b y x) :=
  crooks_pointwise π p (heatBath π b) hπ hp
    (heatBath_preserves_pos hπ p hp)
    (fun x y => heatBath_detailedBalance x y) x y

/-- The conditional-resampling repair kernel satisfies finite Onsager
reciprocity. -/
theorem heatBath_correlation_symm (f g : Ω → ℝ) :
    ∑ x, ∑ y, π x * heatBath π b x y * f x * g y
      = ∑ x, ∑ y, π x * heatBath π b x y * g x * f y :=
  correlation_symm π (heatBath π b)
    (fun x y => heatBath_detailedBalance x y) f g

/-- **Pythagorean identity for one repair step.** Relative entropy to
the reference law splits exactly across the repaired image:
`D(p ‖ π) = D(p ‖ push p K) + D(push p K ‖ π)` for the
conditional-resampling kernel `K` at every strictly positive state
`p`. The proof uses the closed form of the repair step, whose
log-ratio against the reference is measurable through the visible
datum, together with fibre-mean conservation; every sum ranges over
the state space, so the visible alphabet carries no finiteness
hypothesis. Combined with `sigma_mean_eq_kl_descent` this identifies
the mean entropy production of one repair step with
`D(p ‖ push p K)`. -/
theorem heatBath_kl_pythagorean (hπ : ∀ x, 0 < π x)
    (p : Ω → ℝ) (hp : ∀ x, 0 < p x) :
    kl p π = kl p (push p (heatBath π b))
      + kl (push p (heatBath π b)) π := by
  have hq : ∀ y, 0 < push p (heatBath π b) y :=
    heatBath_preserves_pos hπ p hp
  have hratio : ∀ y, push p (heatBath π b) y / π y
      = fiberMass p b y / fiberMass π b y := by
    intro y
    rw [push_heatBath_eq p y]
    have hπy := (hπ y).ne'
    field_simp
  have hfib : ∀ x y, b x = b y →
      Real.log (push p (heatBath π b) x / π x)
        = Real.log (push p (heatBath π b) y / π y) := by
    intro x y hxy
    rw [hratio x, hratio y, fiberMass_congr (π := p) hxy,
      fiberMass_congr (π := π) hxy]
  have hterm : ∀ x, klTerm (p x) (π x)
      = klTerm (p x) (push p (heatBath π b) x)
        + p x * Real.log (push p (heatBath π b) x / π x) := by
    intro x
    rw [klTerm_of_pos (ne_of_gt (hp x)), klTerm_of_pos (ne_of_gt (hp x)),
      Real.log_div (ne_of_gt (hp x)) (ne_of_gt (hπ x)),
      Real.log_div (ne_of_gt (hp x)) (ne_of_gt (hq x)),
      Real.log_div (ne_of_gt (hq x)) (ne_of_gt (hπ x))]
    ring
  have hmean : ∑ y, push p (heatBath π b) y
        * Real.log (push p (heatBath π b) y / π y)
      = ∑ x, p x * Real.log (push p (heatBath π b) x / π x) :=
    push_heatBath_fixes_mean hπ
      (fun x => Real.log (push p (heatBath π b) x / π x)) hfib p
  have hklq : kl (push p (heatBath π b)) π
      = ∑ y, push p (heatBath π b) y
          * Real.log (push p (heatBath π b) y / π y) :=
    Finset.sum_congr rfl fun y _ =>
      klTerm_of_pos (ne_of_gt (hq y)) (π y)
  calc kl p π
      = ∑ x, (klTerm (p x) (push p (heatBath π b) x)
          + p x * Real.log (push p (heatBath π b) x / π x)) :=
        Finset.sum_congr rfl fun x _ => hterm x
    _ = kl p (push p (heatBath π b))
        + ∑ x, p x * Real.log (push p (heatBath π b) x / π x) :=
        Finset.sum_add_distrib
    _ = kl p (push p (heatBath π b))
        + kl (push p (heatBath π b)) π := by
        rw [← hmean, hklq]

end HeatBathInstantiation

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.integral_fluctuation
#print axioms OPH.Thermodynamics.crooks_pointwise
#print axioms OPH.Thermodynamics.crooks_level_set
#print axioms OPH.Thermodynamics.sigma_mean_eq_kl_descent
#print axioms OPH.Thermodynamics.correlation_symm
#print axioms OPH.Thermodynamics.heatBath_integral_fluctuation
#print axioms OPH.Thermodynamics.heatBath_crooks
#print axioms OPH.Thermodynamics.heatBath_correlation_symm
#print axioms OPH.Thermodynamics.heatBath_kl_pythagorean
