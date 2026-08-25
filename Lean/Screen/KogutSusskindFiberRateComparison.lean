import Mathlib

set_option autoImplicit false

namespace OPH.KogutSusskindFiberRate

/-!
# Kogut--Susskind fiber rates and the variable-rate Dirichlet comparison on an independent product of two-point fibers

The Yang--Mills gap paper's finite `Z2` diagnostic (`code/yang_mills/
z2_finite_transfer_receipt.py`, schema
`oph.yang_mills.z2_finite_transfer_receipt.v2`) evaluates the ground-state
Doob transform `L(o, o') = (H(o, o') - e0 δ) Ω(o') / Ω(o)` of the
Kogut--Susskind Hamiltonian `H = -λ Σ_l X_l - Σ_p U_p` on the gauge-orbit
space, with stationary law `π = Ω²` and one heat-bath collar per spatial
link on the two-point fiber `{o, X_l o}`.  On that fiber the transform is
single-flip with the fiber-dependent rate `c_l(o) = λ (r + 1/r)`,
`r = Ω(o) / Ω(X_l o)`.  This module proves, at the finite exact level:

* the scalar rate identities: `c = λ (r + 1/r)` is the unique rate whose
  heat-bath off-diagonal reproduces the Doob off-diagonal, it is invariant
  under the fiber flip `r ↦ 1/r` (so it is constant on its own two-point
  fiber), and `c ≥ 2λ` with equality exactly at `r = 1` (both directions of
  the arithmetic-geometric mean inequality);
* on an arbitrary finite state space, the variable-rate Dirichlet
  domination: a generator term `c ⬝ (I - E)` with a rate function bounded
  below by `c_*` and constant on the fibers of the `π`-symmetric idempotent
  conditional expectation `E` dominates `c_*` times the unit-rate term
  `I - E` in the Dirichlet pairing, and the same for a finite sum of such
  collar terms;
* on the independent product of `n` two-point fibers with product weights
  and one rate `c_l` per fiber, the full eigenstructure of the product
  generator `Σ_l c_l (I - E_l)`: the Walsh family is an orthogonal
  eigenbasis with eigenvalue `Σ_{l ∈ S} c_l` on the level set `S`, and the
  exact best Dirichlet (Poincaré) constant on mean-zero observables is
  `min_l c_l`; with the Kogut--Susskind rates this is at least `2λ`
  uniformly in every fiber ratio `r_l > 0`;
* concrete two-fiber and three-fiber instances with distinct rational
  ratios, their exact rational rate tables, subset-sum eigenvalue tables,
  and exact gaps, including the `r = 1` equality case.

## Sign and convention re-derivation (recorded, not transcribed)

The committed producer forms `L[i, j] = (H[i, j] - e0 δ_{ij}) ω_j / ω_i` and
reads the rate from the off-diagonal by
`c = -L(o, o') (π(o) + π(o')) / π(o')`.  Re-derivation: the Kogut--Susskind
hopping element is `H(o, X_l o) = -λ`, so the Doob off-diagonal is
`L(o, o') = -λ Ω(o') / Ω(o) = -λ / r` with `r = Ω(o) / Ω(X_l o)`.  The
`π`-preserving two-point heat bath has off-diagonal
`-c π(o') / (π(o) + π(o'))`, and with `π = Ω²` the conditional flip weight
is `π(o') / (π(o) + π(o')) = 1 / (1 + r²)`.  Matching the two off-diagonals
gives `c / (1 + r²) = λ / r`, hence `c = λ (1 + r²) / r = λ (r + 1/r)`,
which agrees with the committed producer's extraction and with the paper's
displayed identity.  Under the flip the ratio inverts, `r ↦ 1/r`, and
`λ (r + 1/r)` is unchanged, so the rate is a fiber constant.  These steps
are proved below as `fiberRate_heatBath_offdiagonal`,
`heatBath_offdiagonal_unique_rate`, and `fiberRate_flip`.

## What is not proved here

The independent-product model is not the gauge-orbit quotient: on the
committed `Z2` orbit space the stationary law `π = Ω²` is not a product
measure and the rates `c_l(o)` vary over orbits, so the product
eigenstructure theorems below do not apply to it directly.  The committed
quotient structure (`Z2GaugeOrbits` in the producer) has no Lean
formalization, and no quotient-space spectral-gap theorem is proved here.
The exact open step is the one the 2026-08-24 deep audit names:
"The next proof is a variable-rate approximate-tensorization comparison on
the gauge-orbit quotient, followed by an anisotropic Wilson-to-Hamiltonian
scan."  The domination lemmas in this file supply the rate-comparison half
of that step under explicit hypotheses (a `π`-symmetric idempotent
conditional expectation and a fiber-constant rate function); the
approximate-tensorization half, which controls the variance of the
non-product orbit law by the sum of fiber variances, is open, and so are
the anisotropic Wilson-to-Hamiltonian limit, refinement and volume
extension, and every continuum statement.  Nothing here is a mass-gap
theorem, a Clay-problem step, an all-coupling result, or evidence for OPH;
no data enters and no comparison is scored.  Self-adjointness of the
product generator with respect to the weighted pairing is not separately
stated; the gap language below is carried entirely by the proved
eigenstructure and the two-sided Poincaré characterization.
-/

noncomputable section

/-! ## The scalar Kogut--Susskind fiber rate -/

/-- The Kogut--Susskind Doob-transform collar rate on one two-point fiber
with ground-state ratio `r = Ω(o) / Ω(X_l o)` and hopping coefficient
`λ`.  The physical range is `r > 0`, `λ > 0`; outside it the formula is a
junk-value extension via Mathlib's `1 / 0 = 0` convention. -/
def fiberRate (lam r : ℝ) : ℝ := lam * (r + 1 / r)

/-- The rate is invariant under the fiber flip `r ↦ 1/r`, so it is constant
on its own two-point fiber.  The identity holds for every real `r` by the
`1 / 0 = 0` convention; it is used on `r > 0`. -/
theorem fiberRate_flip (lam r : ℝ) : fiberRate lam (1 / r) = fiberRate lam r := by
  unfold fiberRate
  rw [one_div_one_div]
  ring

/-- Arithmetic-geometric mean floor: `2λ ≤ λ (r + 1/r)` for `λ > 0`,
`r > 0`. -/
theorem two_mul_le_fiberRate {lam r : ℝ} (hlam : 0 < lam) (hr : 0 < r) :
    2 * lam ≤ fiberRate lam r := by
  unfold fiberRate
  have hne : r ≠ 0 := ne_of_gt hr
  have hident : r + 1 / r - 2 = (r - 1) ^ 2 / r := by
    field_simp
    ring
  have hnn : 0 ≤ (r - 1) ^ 2 / r := div_nonneg (sq_nonneg _) hr.le
  have hge : (2 : ℝ) ≤ r + 1 / r := by linarith
  nlinarith [mul_le_mul_of_nonneg_left hge hlam.le]

/-- Equality case of the floor: `λ (r + 1/r) = 2λ` exactly at `r = 1`
(both directions), for `λ > 0`, `r > 0`. -/
theorem fiberRate_eq_two_mul_iff {lam r : ℝ} (hlam : 0 < lam) (hr : 0 < r) :
    fiberRate lam r = 2 * lam ↔ r = 1 := by
  constructor
  · intro h
    have hne : lam ≠ 0 := ne_of_gt hlam
    have h2 : lam * (r + 1 / r) = lam * 2 := by
      unfold fiberRate at h
      linarith
    have h3 : r + 1 / r = 2 := mul_left_cancel₀ hne h2
    have hident : r + 1 / r - 2 = (r - 1) ^ 2 / r := by
      have hrne : r ≠ 0 := ne_of_gt hr
      field_simp
      ring
    have h4 : (r - 1) ^ 2 / r = 0 := by linarith
    have h5 : (r - 1) ^ 2 = 0 := by
      rcases div_eq_zero_iff.mp h4 with h | h
      · exact h
      · exact absurd h (ne_of_gt hr)
    have h6 : r - 1 = 0 := pow_eq_zero_iff (two_ne_zero).mp h5
    linarith
  · rintro rfl
    unfold fiberRate
    norm_num

/-- Convention receipt, forward direction: the fiber rate times the
heat-bath conditional flip weight `1 / (1 + r²)` reproduces the Doob
off-diagonal magnitude `λ / r`.  With `π = Ω²` the flip weight is
`π(o') / (π(o) + π(o')) = 1 / (1 + r²)`. -/
theorem fiberRate_heatBath_offdiagonal {lam r : ℝ} (hr : 0 < r) :
    fiberRate lam r * (1 / (1 + r ^ 2)) = lam / r := by
  have hrne : r ≠ 0 := ne_of_gt hr
  have hne : (1 : ℝ) + r ^ 2 ≠ 0 := by positivity
  unfold fiberRate
  field_simp
  ring

/-- Convention receipt, uniqueness direction: the fiber rate is the unique
rate whose heat-bath off-diagonal matches the Doob off-diagonal.  This is
the committed producer's extraction
`c = -L(o,o') (π(o) + π(o')) / π(o')` in exact form. -/
theorem heatBath_offdiagonal_unique_rate {lam r c : ℝ} (hr : 0 < r)
    (h : c * (1 / (1 + r ^ 2)) = lam / r) : c = fiberRate lam r := by
  have hrne : r ≠ 0 := ne_of_gt hr
  have hne : (1 : ℝ) + r ^ 2 ≠ 0 := by positivity
  unfold fiberRate
  field_simp at h ⊢
  nlinarith [h]

/-! ## Dirichlet pairing on a finite state space and variable-rate domination

This section is stated on an arbitrary finite state space with an abstract
nonnegative weight, so that it applies verbatim to a gauge-orbit quotient
once that quotient's stationary law and heat-bath conditional expectations
are formalized.  Nothing in it requires a product structure. -/

section GeneralStateSpace

variable {X : Type*} [Fintype X]

/-- Weighted pairing `⟨f, g⟩_π = Σ_x π(x) f(x) g(x)` on a finite state
space.  No normalization or positivity of `π` is built in; hypotheses
carry it where needed. -/
def statePairing (pi : X → ℝ) (f g : X → ℝ) : ℝ := ∑ x, pi x * f x * g x

theorem statePairing_comm (pi : X → ℝ) (f g : X → ℝ) :
    statePairing pi f g = statePairing pi g f := by
  unfold statePairing
  exact Finset.sum_congr rfl fun x _ => by ring

theorem statePairing_zero_left (pi : X → ℝ) (g : X → ℝ) :
    statePairing pi 0 g = 0 := by
  unfold statePairing
  simp

theorem statePairing_sub_left (pi : X → ℝ) (f g h : X → ℝ) :
    statePairing pi (f - g) h = statePairing pi f h - statePairing pi g h := by
  unfold statePairing
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun x _ => by
    simp only [Pi.sub_apply]
    ring

theorem statePairing_sub_right (pi : X → ℝ) (f g h : X → ℝ) :
    statePairing pi f (g - h) = statePairing pi f g - statePairing pi f h := by
  rw [statePairing_comm pi f (g - h), statePairing_sub_left,
    statePairing_comm pi g f, statePairing_comm pi h f]

theorem statePairing_sum_left {ι : Type*} (pi : X → ℝ) (t : Finset ι)
    (f : ι → X → ℝ) (g : X → ℝ) :
    statePairing pi (∑ i ∈ t, f i) g = ∑ i ∈ t, statePairing pi (f i) g := by
  unfold statePairing
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun x _ => by
    rw [Finset.sum_apply, Finset.mul_sum, Finset.sum_mul]

theorem statePairing_sum_right {ι : Type*} (pi : X → ℝ) (t : Finset ι)
    (f : X → ℝ) (g : ι → X → ℝ) :
    statePairing pi f (∑ i ∈ t, g i) = ∑ i ∈ t, statePairing pi f (g i) := by
  rw [statePairing_comm, statePairing_sum_left]
  exact Finset.sum_congr rfl fun i _ => statePairing_comm pi (g i) f

theorem statePairing_smul_left (pi : X → ℝ) (a : ℝ) (f g : X → ℝ) :
    statePairing pi (a • f) g = a * statePairing pi f g := by
  unfold statePairing
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun x _ => by
    simp only [Pi.smul_apply, smul_eq_mul]
    ring

theorem statePairing_smul_right (pi : X → ℝ) (a : ℝ) (f g : X → ℝ) :
    statePairing pi f (a • g) = a * statePairing pi f g := by
  rw [statePairing_comm, statePairing_smul_left, statePairing_comm]

/-- Pointwise multiplication commutes across the pairing:
`⟨u, c ⬝ v⟩ = ⟨c ⬝ u, v⟩`. -/
theorem statePairing_mul_shift (pi : X → ℝ) (c u v : X → ℝ) :
    statePairing pi u (c * v) = statePairing pi (c * u) v := by
  unfold statePairing
  exact Finset.sum_congr rfl fun x _ => by
    simp only [Pi.mul_apply]
    ring

/-- Rate floor inside the pairing: if `π ≥ 0` and `c ≥ c_*` pointwise then
`c_* ⟨g, g⟩ ≤ ⟨g, c ⬝ g⟩`. -/
theorem statePairing_rate_lower (pi : X → ℝ) (hpi : ∀ x, 0 ≤ pi x)
    {c : X → ℝ} {cstar : ℝ} (hc : ∀ x, cstar ≤ c x) (g : X → ℝ) :
    cstar * statePairing pi g g ≤ statePairing pi g (c * g) := by
  unfold statePairing
  rw [Finset.mul_sum, ← sub_nonneg, ← Finset.sum_sub_distrib]
  refine Finset.sum_nonneg fun x _ => ?_
  have hident : pi x * g x * ((c * g) x) - cstar * (pi x * g x * g x)
      = pi x * (c x - cstar) * (g x) ^ 2 := by
    simp only [Pi.mul_apply]
    ring
  rw [hident]
  exact mul_nonneg (mul_nonneg (hpi x) (by linarith [hc x])) (sq_nonneg _)

section CondExp

variable (pi : X → ℝ) (E : (X → ℝ) → (X → ℝ))

/-- For a `π`-symmetric idempotent `E`, the conditional part of `f` is
pairing-orthogonal to the residual `f - E f`. -/
theorem statePairing_condExp_residual
    (hsym : ∀ f g, statePairing pi (E f) g = statePairing pi f (E g))
    (hidem : ∀ f, E (E f) = E f) (f : X → ℝ) :
    statePairing pi (E f) (f - E f) = 0 := by
  rw [statePairing_sub_right, hsym f f, hsym f (E f), hidem f, sub_self]

/-- For a `π`-symmetric idempotent `E` and a rate function constant on the
fibers of `E` (stated operationally: multiplication by `c` commutes with
`E`), the conditional part is also orthogonal to the rate-weighted
residual. -/
theorem statePairing_condExp_rate_residual
    (hsym : ∀ f g, statePairing pi (E f) g = statePairing pi f (E g))
    (hidem : ∀ f, E (E f) = E f)
    {c : X → ℝ} (hcomm : ∀ g, E (c * g) = c * E g) (f : X → ℝ) :
    statePairing pi (E f) (c * (f - E f)) = 0 := by
  rw [statePairing_mul_shift, statePairing_sub_right, ← hsym (c * E f) f,
    hcomm, hidem, sub_self]

/-- The Dirichlet form of a collar term is carried by the residual:
`⟨f, f - E f⟩ = ⟨f - E f, f - E f⟩`. -/
theorem statePairing_residual_split
    (hsym : ∀ f g, statePairing pi (E f) g = statePairing pi f (E g))
    (hidem : ∀ f, E (E f) = E f) (f : X → ℝ) :
    statePairing pi f (f - E f) = statePairing pi (f - E f) (f - E f) := by
  have h0 := statePairing_condExp_residual pi E hsym hidem f
  have hsplit := statePairing_sub_left pi f (E f) (f - E f)
  linarith

/-- Rate-weighted version of `statePairing_residual_split`. -/
theorem statePairing_rate_residual_split
    (hsym : ∀ f g, statePairing pi (E f) g = statePairing pi f (E g))
    (hidem : ∀ f, E (E f) = E f)
    {c : X → ℝ} (hcomm : ∀ g, E (c * g) = c * E g) (f : X → ℝ) :
    statePairing pi f (c * (f - E f))
      = statePairing pi (f - E f) (c * (f - E f)) := by
  have h0 := statePairing_condExp_rate_residual pi E hsym hidem hcomm f
  have hsplit := statePairing_sub_left pi f (E f) (c * (f - E f))
  linarith

/-- Variable-rate Dirichlet domination, one collar: a collar term with a
fiber-constant rate bounded below by `c_*` dominates `c_*` times the
unit-rate collar term in the Dirichlet pairing.  This is the exact form of
the Yang--Mills paper's claim that the finite floor survives
fiber-dependent rates, on any finite state space. -/
theorem dirichlet_domination_single (hpi : ∀ x, 0 ≤ pi x)
    (hsym : ∀ f g, statePairing pi (E f) g = statePairing pi f (E g))
    (hidem : ∀ f, E (E f) = E f)
    {c : X → ℝ} (hcomm : ∀ g, E (c * g) = c * E g)
    {cstar : ℝ} (hc : ∀ x, cstar ≤ c x) (f : X → ℝ) :
    cstar * statePairing pi f (f - E f) ≤ statePairing pi f (c * (f - E f)) := by
  rw [statePairing_residual_split pi E hsym hidem f,
    statePairing_rate_residual_split pi E hsym hidem hcomm f]
  exact statePairing_rate_lower pi hpi hc (f - E f)

end CondExp

/-- Variable-rate Dirichlet domination for a finite collar family: the
variable-rate generator `Σ_l c_l ⬝ (I - E_l)` dominates `c_*` times the
unit-rate heat-bath sum `Σ_l (I - E_l)` in the Dirichlet pairing, under
the stated hypotheses on every collar.  Applying this to a formalized
gauge-orbit quotient is the rate-comparison half of the audit's named next
step; the approximate-tensorization half is open. -/
theorem dirichlet_domination_family {ι : Type*} [Fintype ι]
    (pi : X → ℝ) (hpi : ∀ x, 0 ≤ pi x)
    (E : ι → (X → ℝ) → (X → ℝ))
    (hsym : ∀ l f g, statePairing pi (E l f) g = statePairing pi f (E l g))
    (hidem : ∀ l f, E l (E l f) = E l f)
    (c : ι → X → ℝ) (hcomm : ∀ l g, E l (c l * g) = c l * E l g)
    {cstar : ℝ} (hc : ∀ l x, cstar ≤ c l x) (f : X → ℝ) :
    cstar * statePairing pi f (∑ l, (f - E l f))
      ≤ statePairing pi f (∑ l, c l * (f - E l f)) := by
  rw [statePairing_sum_right, statePairing_sum_right, Finset.mul_sum]
  exact Finset.sum_le_sum fun l _ =>
    dirichlet_domination_single pi (E l) hpi (hsym l) (hidem l) (hcomm l)
      (hc l) f

end GeneralStateSpace

/-! ## Two-point inhabitant of the domination hypotheses

The abstract hypotheses above are inhabited by the genuine single-fiber
object: the two-point conditional expectation onto constants.  With the
Kogut--Susskind rate this gives the exact single-fiber Dirichlet floor
`2λ`. -/

/-- Conditional expectation onto constants on a two-point state space with
weights `p`. -/
def twoPointCondExp (p : Bool → ℝ) : (Bool → ℝ) → (Bool → ℝ) :=
  fun f _ => p false * f false + p true * f true

/-- Idempotence of the two-point conditional expectation under weight
normalization. -/
theorem twoPointCondExp_idem {p : Bool → ℝ} (hsum : p false + p true = 1)
    (f : Bool → ℝ) : twoPointCondExp p (twoPointCondExp p f) = twoPointCondExp p f := by
  funext b
  unfold twoPointCondExp
  linear_combination (p false * f false + p true * f true) * hsum

/-- `p`-symmetry of the two-point conditional expectation with `π = p`. -/
theorem twoPointCondExp_symm (p : Bool → ℝ) (f g : Bool → ℝ) :
    statePairing p (twoPointCondExp p f) g
      = statePairing p f (twoPointCondExp p g) := by
  unfold statePairing twoPointCondExp
  rw [Fintype.sum_bool, Fintype.sum_bool]
  ring

/-- A constant rate commutes with the two-point conditional expectation:
the single two-point fiber carries one rate value, which is the exact
sense in which the Kogut--Susskind rate is constant on its own fiber. -/
theorem twoPointCondExp_const_rate_comm (p : Bool → ℝ) (c0 : ℝ) (g : Bool → ℝ) :
    twoPointCondExp p ((fun _ => c0) * g) = (fun _ => c0) * twoPointCondExp p g := by
  funext b
  simp only [twoPointCondExp, Pi.mul_apply]
  ring

/-- Single-fiber exact Dirichlet floor: on one two-point fiber with
nonnegative normalized weights, the Kogut--Susskind collar term with rate
`λ (r + 1/r)` dominates `2λ` times the unit-rate collar term.  This
inhabits every hypothesis of `dirichlet_domination_single` with the
committed rate. -/
theorem twoPoint_kogutSusskind_domination {p : Bool → ℝ}
    (hp : ∀ b, 0 ≤ p b) (hsum : p false + p true = 1)
    {lam r : ℝ} (hlam : 0 < lam) (hr : 0 < r) (f : Bool → ℝ) :
    (2 * lam) * statePairing p f (f - twoPointCondExp p f)
      ≤ statePairing p f
          ((fun _ => fiberRate lam r) * (f - twoPointCondExp p f)) :=
  dirichlet_domination_single p (twoPointCondExp p) hp
    (twoPointCondExp_symm p) (fun g => twoPointCondExp_idem hsum g)
    (twoPointCondExp_const_rate_comm p (fiberRate lam r))
    (fun _ => two_mul_le_fiberRate hlam hr) f

/-! ## The independent product of two-point fibers

State space `Fin n → Bool`, product weights `π(x) = ∏_l w_l(x_l)`, one
collar per fiber, one rate per fiber.  This is the idealized constant-rate
factorized model; it is not the gauge-orbit quotient. -/

section ProductModel

variable {n : ℕ}

/-- Configurations of `n` two-point fibers. -/
abbrev Config (n : ℕ) := Fin n → Bool

/-- Real observables on the product configuration space. -/
abbrev Obs (n : ℕ) := Config n → ℝ

variable {w : Fin n → Bool → ℝ}

/-- Product weight `π(x) = ∏_l w_l(x_l)`. -/
def prodWeight (w : Fin n → Bool → ℝ) (x : Config n) : ℝ := ∏ l, w l (x l)

/-- The fiber conditional expectation `E_l`: average coordinate `l`
against its fiber weights, leaving the other coordinates fixed. -/
def condExp (w : Fin n → Bool → ℝ) (l : Fin n) : Obs n →ₗ[ℝ] Obs n where
  toFun f := fun x =>
    w l false * f (Function.update x l false) + w l true * f (Function.update x l true)
  map_add' f g := by
    funext x
    simp only [Pi.add_apply]
    ring
  map_smul' a f := by
    funext x
    simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
    ring

@[simp] theorem condExp_apply (l : Fin n) (f : Obs n) (x : Config n) :
    condExp w l f x
      = w l false * f (Function.update x l false)
        + w l true * f (Function.update x l true) := rfl

/-- The variable-rate product generator `Σ_l c_l (I - E_l)` with one
constant rate per fiber. -/
def generator (w : Fin n → Bool → ℝ) (c : Fin n → ℝ) : Obs n →ₗ[ℝ] Obs n :=
  ∑ l, c l • (LinearMap.id - condExp w l)

theorem generator_apply (c : Fin n → ℝ) (f : Obs n) :
    generator w c f = ∑ l, c l • (f - condExp w l f) := by
  unfold generator
  rw [LinearMap.sum_apply]
  exact Finset.sum_congr rfl fun l _ => by
    rw [LinearMap.smul_apply, LinearMap.sub_apply, LinearMap.id_apply]

/-- Signed fiber function with `π_l`-mean zero: `ψ_l(true) = w_l(false)`,
`ψ_l(false) = -w_l(true)`. -/
def fiberSign (w : Fin n → Bool → ℝ) (l : Fin n) (b : Bool) : ℝ :=
  cond b (w l false) (-(w l true))

@[simp] theorem fiberSign_true (l : Fin n) : fiberSign w l true = w l false := rfl

@[simp] theorem fiberSign_false (l : Fin n) : fiberSign w l false = -(w l true) := rfl

/-- Walsh function of the level set `S`: the product of the signed fiber
functions over `S`. -/
def walsh (w : Fin n → Bool → ℝ) (S : Finset (Fin n)) : Obs n :=
  fun x => ∏ l ∈ S, fiberSign w l (x l)

theorem walsh_eq_prod_univ (S : Finset (Fin n)) (x : Config n) :
    walsh w S x = ∏ l, (if l ∈ S then fiberSign w l (x l) else 1) := by
  rw [Finset.prod_ite_mem, Finset.univ_inter]
  rfl

theorem walsh_update_of_notMem {S : Finset (Fin n)} {l : Fin n} (hl : l ∉ S)
    (x : Config n) (b : Bool) :
    walsh w S (Function.update x l b) = walsh w S x := by
  unfold walsh
  refine Finset.prod_congr rfl fun m hm => ?_
  have hml : m ≠ l := fun h => hl (h ▸ hm)
  rw [Function.update_of_ne hml]

theorem walsh_update_of_mem {S : Finset (Fin n)} {l : Fin n} (hl : l ∈ S)
    (x : Config n) (b : Bool) :
    walsh w S (Function.update x l b)
      = fiberSign w l b * ∏ m ∈ S.erase l, fiberSign w m (x m) := by
  unfold walsh
  rw [← Finset.mul_prod_erase S _ hl, Function.update_self]
  congr 1
  refine Finset.prod_congr rfl fun m hm => ?_
  rw [Function.update_of_ne (Finset.ne_of_mem_erase hm)]

/-- A Walsh function over a level set containing `l` is annihilated by
`E_l`: the fiber mean of the signed fiber function vanishes. -/
theorem condExp_walsh_of_mem {S : Finset (Fin n)} {l : Fin n} (hl : l ∈ S) :
    condExp w l (walsh w S) = 0 := by
  funext x
  rw [condExp_apply, walsh_update_of_mem hl, walsh_update_of_mem hl]
  simp only [fiberSign_true, fiberSign_false, Pi.zero_apply]
  ring

/-- A Walsh function over a level set avoiding `l` is fixed by `E_l`
under weight normalization. -/
theorem condExp_walsh_of_notMem {S : Finset (Fin n)} {l : Fin n} (hl : l ∉ S)
    (hsum : w l false + w l true = 1) :
    condExp w l (walsh w S) = walsh w S := by
  funext x
  rw [condExp_apply, walsh_update_of_notMem hl, walsh_update_of_notMem hl,
    ← add_mul, hsum, one_mul]

/-- Explicit eigenstructure: the Walsh function of the level set `S` is an
eigenvector of the product generator with eigenvalue `Σ_{l ∈ S} c_l`. -/
theorem generator_walsh (hsum : ∀ l, w l false + w l true = 1)
    (c : Fin n → ℝ) (S : Finset (Fin n)) :
    generator w c (walsh w S) = (∑ l ∈ S, c l) • walsh w S := by
  rw [generator_apply]
  have h : ∀ l : Fin n,
      c l • (walsh w S - condExp w l (walsh w S))
        = if l ∈ S then c l • walsh w S else 0 := by
    intro l
    by_cases hl : l ∈ S
    · rw [if_pos hl, condExp_walsh_of_mem hl, sub_zero]
    · rw [if_neg hl, condExp_walsh_of_notMem hl (hsum l), sub_self, smul_zero]
  calc
    ∑ l, c l • (walsh w S - condExp w l (walsh w S))
        = ∑ l, (if l ∈ S then c l • walsh w S else 0) :=
          Finset.sum_congr rfl fun l _ => h l
    _ = ∑ l ∈ S, c l • walsh w S := by
          rw [Finset.sum_ite_mem, Finset.univ_inter]
    _ = (∑ l ∈ S, c l) • walsh w S := (Finset.sum_smul).symm

/-- Per-coordinate factor of the pairing of two Walsh functions. -/
def pairFactor (w : Fin n → Bool → ℝ) (S T : Finset (Fin n)) (l : Fin n) : ℝ :=
  ∑ b, w l b * (if l ∈ S then fiberSign w l b else 1)
    * (if l ∈ T then fiberSign w l b else 1)

/-- The pairing of two Walsh functions factorizes over coordinates. -/
theorem statePairing_walsh_walsh (S T : Finset (Fin n)) :
    statePairing (prodWeight w) (walsh w S) (walsh w T)
      = ∏ l, pairFactor w S T l := by
  unfold statePairing
  have hx : ∀ x : Config n,
      prodWeight w x * walsh w S x * walsh w T x
        = ∏ l, (w l (x l) * (if l ∈ S then fiberSign w l (x l) else 1)
            * (if l ∈ T then fiberSign w l (x l) else 1)) := by
    intro x
    rw [walsh_eq_prod_univ S x, walsh_eq_prod_univ T x]
    unfold prodWeight
    rw [← Finset.prod_mul_distrib, ← Finset.prod_mul_distrib]
  rw [Finset.sum_congr rfl fun x _ => hx x]
  unfold pairFactor
  rw [← Finset.sum_prod_piFinset, Fintype.piFinset_univ]

theorem pairFactor_of_notMem_notMem {S T : Finset (Fin n)} {l : Fin n}
    (hS : l ∉ S) (hT : l ∉ T) (hsum : w l false + w l true = 1) :
    pairFactor w S T l = 1 := by
  unfold pairFactor
  rw [Fintype.sum_bool]
  simp only [if_neg hS, if_neg hT]
  linarith

theorem pairFactor_of_mem_notMem {S T : Finset (Fin n)} {l : Fin n}
    (hS : l ∈ S) (hT : l ∉ T) :
    pairFactor w S T l = 0 := by
  unfold pairFactor
  rw [Fintype.sum_bool]
  simp only [if_pos hS, if_neg hT, fiberSign_true, fiberSign_false]
  ring

theorem pairFactor_of_notMem_mem {S T : Finset (Fin n)} {l : Fin n}
    (hS : l ∉ S) (hT : l ∈ T) :
    pairFactor w S T l = 0 := by
  unfold pairFactor
  rw [Fintype.sum_bool]
  simp only [if_neg hS, if_pos hT, fiberSign_true, fiberSign_false]
  ring

theorem pairFactor_of_mem_mem {S T : Finset (Fin n)} {l : Fin n}
    (hS : l ∈ S) (hT : l ∈ T) (hsum : w l false + w l true = 1) :
    pairFactor w S T l = w l false * w l true := by
  unfold pairFactor
  rw [Fintype.sum_bool]
  simp only [if_pos hS, if_pos hT, fiberSign_true, fiberSign_false]
  linear_combination (w l false * w l true) * hsum

/-- Orthogonality of distinct Walsh functions in the product pairing. -/
theorem statePairing_walsh_orthogonal (hsum : ∀ l, w l false + w l true = 1)
    {S T : Finset (Fin n)} (hST : S ≠ T) :
    statePairing (prodWeight w) (walsh w S) (walsh w T) = 0 := by
  obtain ⟨l, hl⟩ : ∃ l, (l ∈ S ∧ l ∉ T) ∨ (l ∈ T ∧ l ∉ S) := by
    by_contra hcon
    push Not at hcon
    refine hST (Finset.ext fun l => ⟨fun h => ?_, fun h => ?_⟩)
    · exact (hcon l).1 h
    · exact (hcon l).2 h
  rw [statePairing_walsh_walsh]
  refine Finset.prod_eq_zero (Finset.mem_univ l) ?_
  rcases hl with ⟨hS, hT⟩ | ⟨hT, hS⟩
  · exact pairFactor_of_mem_notMem hS hT
  · exact pairFactor_of_notMem_mem hS hT

/-- Norm of a Walsh function: `∏_{l ∈ S} w_l(false) w_l(true)`. -/
theorem statePairing_walsh_self (hsum : ∀ l, w l false + w l true = 1)
    (S : Finset (Fin n)) :
    statePairing (prodWeight w) (walsh w S) (walsh w S)
      = ∏ l ∈ S, (w l false * w l true) := by
  rw [statePairing_walsh_walsh]
  calc
    ∏ l, pairFactor w S S l
        = ∏ l, (if l ∈ S then w l false * w l true else 1) := by
          refine Finset.prod_congr rfl fun l _ => ?_
          by_cases hl : l ∈ S
          · rw [if_pos hl, pairFactor_of_mem_mem hl hl (hsum l)]
          · rw [if_neg hl, pairFactor_of_notMem_notMem hl hl (hsum l)]
    _ = ∏ l ∈ S, (w l false * w l true) := by
          rw [Finset.prod_ite_mem, Finset.univ_inter]

theorem statePairing_walsh_self_pos (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1) (S : Finset (Fin n)) :
    0 < statePairing (prodWeight w) (walsh w S) (walsh w S) := by
  rw [statePairing_walsh_self hsum S]
  exact Finset.prod_pos fun l _ => mul_pos (hw l false) (hw l true)

/-- The product weight is a probability weight: total mass one. -/
theorem prodWeight_sum_eq_one (hsum : ∀ l, w l false + w l true = 1) :
    ∑ x : Config n, prodWeight w x = 1 := by
  have h := statePairing_walsh_self hsum (∅ : Finset (Fin n))
  rw [Finset.prod_empty] at h
  calc
    ∑ x : Config n, prodWeight w x
        = statePairing (prodWeight w) (walsh w ∅) (walsh w ∅) := by
          unfold statePairing walsh
          exact Finset.sum_congr rfl fun x _ => by
            rw [Finset.prod_empty, mul_one, mul_one]
    _ = 1 := h

/-- Linear independence of the Walsh family. -/
theorem walsh_linearIndependent (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1) :
    LinearIndependent ℝ (fun S : Finset (Fin n) => walsh w S) := by
  rw [Fintype.linearIndependent_iff]
  intro g hg T
  have h0 : statePairing (prodWeight w) (∑ S, g S • walsh w S) (walsh w T) = 0 := by
    rw [hg, statePairing_zero_left]
  rw [statePairing_sum_left] at h0
  rw [Finset.sum_eq_single T
    (fun S _ hST => by
      rw [statePairing_smul_left, statePairing_walsh_orthogonal hsum hST,
        mul_zero])
    (fun h => absurd (Finset.mem_univ T) h)] at h0
  rw [statePairing_smul_left] at h0
  rcases mul_eq_zero.mp h0 with h | h
  · exact h
  · exact absurd h (ne_of_gt (statePairing_walsh_self_pos hw hsum T))

/-- The Walsh functions as a basis of the observable space: `2^n`
orthogonal nonzero vectors in dimension `2^n`. -/
def walshBasis (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1) :
    Basis (Finset (Fin n)) ℝ (Obs n) :=
  basisOfLinearIndependentOfCardEqFinrank (walsh_linearIndependent hw hsum)
    (by
      rw [Fintype.card_finset, Fintype.card_fin,
        Module.finrank_fintype_fun_eq_card, Fintype.card_fun,
        Fintype.card_bool, Fintype.card_fin])

theorem walshBasis_apply (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1) (S : Finset (Fin n)) :
    walshBasis hw hsum S = walsh w S := by
  unfold walshBasis
  rw [coe_basisOfLinearIndependentOfCardEqFinrank]

/-- Expansion of an observable in the Walsh basis. -/
theorem repr_expansion (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1) (f : Obs n) :
    f = ∑ S, (walshBasis hw hsum).repr f S • walsh w S := by
  conv_lhs => rw [← Basis.sum_repr (walshBasis hw hsum) f]
  exact Finset.sum_congr rfl fun S _ => by rw [walshBasis_apply]

theorem statePairing_walsh_right (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1) (S : Finset (Fin n)) (f : Obs n) :
    statePairing (prodWeight w) (walsh w S) f
      = (walshBasis hw hsum).repr f S
          * statePairing (prodWeight w) (walsh w S) (walsh w S) := by
  conv_lhs => rw [repr_expansion hw hsum f]
  rw [statePairing_sum_right]
  rw [Finset.sum_eq_single S
    (fun T _ hT => by
      rw [statePairing_smul_right,
        statePairing_walsh_orthogonal hsum (Ne.symm hT), mul_zero])
    (fun h => absurd (Finset.mem_univ S) h)]
  rw [statePairing_smul_right]

theorem statePairing_left_expansion (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1) (f g : Obs n) :
    statePairing (prodWeight w) f g
      = ∑ S, (walshBasis hw hsum).repr f S
          * statePairing (prodWeight w) (walsh w S) g := by
  conv_lhs => rw [repr_expansion hw hsum f]
  rw [statePairing_sum_left]
  exact Finset.sum_congr rfl fun S _ => statePairing_smul_left _ _ _ _

/-- Diagonal form of the squared norm in the Walsh basis. -/
theorem statePairing_self_diagonal (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1) (f : Obs n) :
    statePairing (prodWeight w) f f
      = ∑ S, ((walshBasis hw hsum).repr f S) ^ 2
          * statePairing (prodWeight w) (walsh w S) (walsh w S) := by
  rw [statePairing_left_expansion hw hsum f f]
  refine Finset.sum_congr rfl fun S _ => ?_
  rw [statePairing_walsh_right hw hsum S f]
  ring

/-- Diagonal form of the Dirichlet form in the Walsh basis, with the
subset-sum eigenvalues explicit. -/
theorem statePairing_generator_diagonal (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1) (c : Fin n → ℝ) (f : Obs n) :
    statePairing (prodWeight w) f (generator w c f)
      = ∑ S, ((walshBasis hw hsum).repr f S) ^ 2 * (∑ l ∈ S, c l)
          * statePairing (prodWeight w) (walsh w S) (walsh w S) := by
  have hLf : generator w c f
      = ∑ T, ((walshBasis hw hsum).repr f T * (∑ l ∈ T, c l)) • walsh w T := by
    conv_lhs => rw [repr_expansion hw hsum f]
    rw [map_sum]
    exact Finset.sum_congr rfl fun T _ => by
      rw [map_smul, generator_walsh hsum c T, smul_smul]
  rw [hLf, statePairing_sum_right]
  refine Finset.sum_congr rfl fun T _ => ?_
  rw [statePairing_smul_right, statePairing_comm,
    statePairing_walsh_right hw hsum T f]
  ring

/-- The `π`-mean of an observable. -/
def piMean (w : Fin n → Bool → ℝ) (f : Obs n) : ℝ :=
  ∑ x, prodWeight w x * f x

theorem statePairing_walsh_empty_left (f : Obs n) :
    statePairing (prodWeight w) (walsh w ∅) f = piMean w f := by
  unfold statePairing piMean walsh
  exact Finset.sum_congr rfl fun x _ => by rw [Finset.prod_empty, mul_one]

theorem repr_empty_of_meanZero (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1) (f : Obs n)
    (hf : piMean w f = 0) :
    (walshBasis hw hsum).repr f ∅ = 0 := by
  have h := statePairing_walsh_right hw hsum ∅ f
  rw [statePairing_walsh_empty_left, hf] at h
  have hm : statePairing (prodWeight w) (walsh w (∅ : Finset (Fin n)))
      (walsh w (∅ : Finset (Fin n))) = 1 := by
    rw [statePairing_walsh_self hsum, Finset.prod_empty]
  rw [hm, mul_one] at h
  exact h.symm

/-- Poincaré lower bound on the product model: every uniform lower bound
on the fiber rates is a Dirichlet constant on mean-zero observables.
Scope: `n` independent two-point fibers, product weights, one constant
rate per fiber. -/
theorem dirichlet_lower_bound (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1)
    {c : Fin n → ℝ} (hc : ∀ l, 0 ≤ c l)
    {gam : ℝ} (hgam : ∀ l, gam ≤ c l)
    {f : Obs n} (hf : piMean w f = 0) :
    gam * statePairing (prodWeight w) f f
      ≤ statePairing (prodWeight w) f (generator w c f) := by
  rw [statePairing_self_diagonal hw hsum f,
    statePairing_generator_diagonal hw hsum c f, Finset.mul_sum]
  refine Finset.sum_le_sum fun S _ => ?_
  by_cases hS : S = ∅
  · subst hS
    rw [repr_empty_of_meanZero hw hsum f hf]
    simp
  · obtain ⟨l0, hl0⟩ := Finset.nonempty_iff_ne_empty.mpr hS
    have he : gam ≤ ∑ l ∈ S, c l :=
      (hgam l0).trans (Finset.single_le_sum (fun l _ => hc l) hl0)
    have hm : (0 : ℝ) ≤ ((walshBasis hw hsum).repr f S) ^ 2
        * statePairing (prodWeight w) (walsh w S) (walsh w S) :=
      mul_nonneg (sq_nonneg _) (statePairing_walsh_self_pos hw hsum S).le
    calc
      gam * (((walshBasis hw hsum).repr f S) ^ 2
          * statePairing (prodWeight w) (walsh w S) (walsh w S))
          ≤ (∑ l ∈ S, c l) * (((walshBasis hw hsum).repr f S) ^ 2
            * statePairing (prodWeight w) (walsh w S) (walsh w S)) :=
        mul_le_mul_of_nonneg_right he hm
      _ = ((walshBasis hw hsum).repr f S) ^ 2 * (∑ l ∈ S, c l)
          * statePairing (prodWeight w) (walsh w S) (walsh w S) := by ring

/-- Sharpness: any Dirichlet constant on mean-zero observables is at most
every fiber rate; the Walsh function of the singleton `{l}` attains rate
`c_l` exactly. -/
theorem rate_le_of_dirichlet_bound (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1)
    {c : Fin n → ℝ} {gam : ℝ}
    (h : ∀ f : Obs n, piMean w f = 0 →
      gam * statePairing (prodWeight w) f f
        ≤ statePairing (prodWeight w) f (generator w c f))
    (l : Fin n) : gam ≤ c l := by
  have hmz : piMean w (walsh w {l}) = 0 := by
    rw [← statePairing_walsh_empty_left]
    exact statePairing_walsh_orthogonal hsum
      (Ne.symm (Finset.singleton_ne_empty l))
  have hb := h (walsh w {l}) hmz
  have hgen : statePairing (prodWeight w) (walsh w {l})
      (generator w c (walsh w {l}))
        = c l * statePairing (prodWeight w) (walsh w {l}) (walsh w {l}) := by
    rw [generator_walsh hsum c {l}, statePairing_smul_right,
      Finset.sum_singleton]
  rw [hgen] at hb
  exact le_of_mul_le_mul_right hb (statePairing_walsh_self_pos hw hsum {l})

/-- Exact two-sided characterization: on the independent product of `n`
two-point fibers with nonnegative rates, `γ` is a Dirichlet constant on
mean-zero observables exactly when `γ ≤ c_l` for every fiber.  Read
together with `generator_walsh`, this identifies the best constant with
the minimum fiber rate: the spectral gap of the product generator is
`min_l c_l`. -/
theorem dirichlet_constant_iff (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1)
    {c : Fin n → ℝ} (hc : ∀ l, 0 ≤ c l) (gam : ℝ) :
    (∀ f : Obs n, piMean w f = 0 →
        gam * statePairing (prodWeight w) f f
          ≤ statePairing (prodWeight w) f (generator w c f))
      ↔ ∀ l, gam ≤ c l :=
  ⟨fun h l => rate_le_of_dirichlet_bound hw hsum h l,
    fun hgam f hf => dirichlet_lower_bound hw hsum hc hgam hf⟩

/-- The best Dirichlet constant is exactly the minimum fiber rate
(`n ≥ 1`): membership is the lower bound, maximality is sharpness at the
minimizing fiber. -/
theorem bestConstant_isGreatest_min_rate
    (hne : (Finset.univ : Finset (Fin n)).Nonempty)
    (hw : ∀ l b, 0 < w l b) (hsum : ∀ l, w l false + w l true = 1)
    {c : Fin n → ℝ} (hc : ∀ l, 0 ≤ c l) :
    IsGreatest
      {gam : ℝ | ∀ f : Obs n, piMean w f = 0 →
        gam * statePairing (prodWeight w) f f
          ≤ statePairing (prodWeight w) f (generator w c f)}
      (Finset.univ.inf' hne c) := by
  constructor
  · exact (dirichlet_constant_iff hw hsum hc _).mpr fun l =>
      Finset.inf'_le c (Finset.mem_univ l)
  · intro gam hgam
    exact Finset.le_inf' hne c fun l _ =>
      rate_le_of_dirichlet_bound hw hsum hgam l

/-! ## The Kogut--Susskind rates on the product model -/

/-- Uniform Kogut--Susskind floor on the product model: with per-fiber
rates `c_l = λ (r_l + 1/r_l)`, the Dirichlet form dominates `2λ` times
the squared norm of every mean-zero observable, uniformly in every fiber
ratio `r_l > 0`. -/
theorem kogutSusskind_uniform_dirichlet_floor (hw : ∀ l b, 0 < w l b)
    (hsum : ∀ l, w l false + w l true = 1)
    {lam : ℝ} (hlam : 0 < lam) {r : Fin n → ℝ} (hr : ∀ l, 0 < r l)
    {f : Obs n} (hf : piMean w f = 0) :
    (2 * lam) * statePairing (prodWeight w) f f
      ≤ statePairing (prodWeight w) f
          (generator w (fun l => fiberRate lam (r l)) f) :=
  dirichlet_lower_bound hw hsum
    (fun l => le_trans (by linarith) (two_mul_le_fiberRate hlam (hr l)))
    (fun l => two_mul_le_fiberRate hlam (hr l)) hf

/-- Fiber weights induced by ground-state ratios under the `π = Ω²`
convention: on the fiber `{o, X_l o}` with `r = Ω(o)/Ω(X_l o)`, the
conditional law is `(r²/(1+r²), 1/(1+r²))`. -/
def ksWeight (r : Fin n → ℝ) : Fin n → Bool → ℝ := fun l b =>
  cond b ((r l) ^ 2 / (1 + (r l) ^ 2)) (1 / (1 + (r l) ^ 2))

theorem ksWeight_pos {r : Fin n → ℝ} (hr : ∀ l, 0 < r l) :
    ∀ l b, 0 < ksWeight r l b := by
  intro l b
  have h1 : (0 : ℝ) < 1 + (r l) ^ 2 := by positivity
  cases b
  · exact div_pos one_pos h1
  · exact div_pos (pow_pos (hr l) 2) h1

theorem ksWeight_sum (r : Fin n → ℝ) :
    ∀ l, ksWeight r l false + ksWeight r l true = 1 := by
  intro l
  have h1 : (1 : ℝ) + (r l) ^ 2 ≠ 0 := by positivity
  unfold ksWeight
  field_simp

/-- The `ksWeight` fiber odds reproduce the squared ground-state ratio,
matching the committed `π = Ω²` convention. -/
theorem ksWeight_ratio {r : Fin n → ℝ} (hr : ∀ l, 0 < r l) (l : Fin n) :
    ksWeight r l true / ksWeight r l false = (r l) ^ 2 := by
  have h1 : (1 : ℝ) + (r l) ^ 2 ≠ 0 := by positivity
  unfold ksWeight
  field_simp

/-- Kogut--Susskind floor with the induced `π = Ω²` fiber weights. -/
theorem kogutSusskind_ksWeight_floor {lam : ℝ} (hlam : 0 < lam)
    {r : Fin n → ℝ} (hr : ∀ l, 0 < r l)
    {f : Obs n} (hf : piMean (ksWeight r) f = 0) :
    (2 * lam) * statePairing (prodWeight (ksWeight r)) f f
      ≤ statePairing (prodWeight (ksWeight r)) f
          (generator (ksWeight r) (fun l => fiberRate lam (r l)) f) :=
  kogutSusskind_uniform_dirichlet_floor (ksWeight_pos hr) (ksWeight_sum r)
    hlam hr hf

end ProductModel

/-! ## Concrete instances with exact rational data -/

section Instances

/-- Ratios of the two-fiber instance: `r = (2, 3)`, `λ = 1`. -/
def rA : Fin 2 → ℝ := ![2, 3]

theorem rA_pos : ∀ l, 0 < rA l := by
  intro l
  fin_cases l <;> norm_num [rA]

/-- Exact rates of the two-fiber instance: `(5/2, 10/3)`. -/
theorem instanceA_rate_values :
    fiberRate 1 (rA 0) = 5 / 2 ∧ fiberRate 1 (rA 1) = 10 / 3 := by
  constructor <;> norm_num [fiberRate, rA]

/-- Exact subset-sum eigenvalue table of the two-fiber instance:
`0, 5/2, 10/3, 35/6`. -/
theorem instanceA_eigenvalue_table :
    (∑ l ∈ (∅ : Finset (Fin 2)), fiberRate 1 (rA l)) = 0 ∧
    (∑ l ∈ ({0} : Finset (Fin 2)), fiberRate 1 (rA l)) = 5 / 2 ∧
    (∑ l ∈ ({1} : Finset (Fin 2)), fiberRate 1 (rA l)) = 10 / 3 ∧
    (∑ l ∈ (Finset.univ : Finset (Fin 2)), fiberRate 1 (rA l)) = 35 / 6 := by
  refine ⟨by simp, ?_, ?_, ?_⟩
  · rw [Finset.sum_singleton]
    norm_num [fiberRate, rA]
  · rw [Finset.sum_singleton]
    norm_num [fiberRate, rA]
  · rw [Fin.sum_univ_two]
    norm_num [fiberRate, rA]

/-- Exact gap of the two-fiber instance: the best Dirichlet constant on
mean-zero observables is `5/2 = min(5/2, 10/3)`, strictly above the
`2λ = 2` floor because neither ratio is `1`. -/
theorem instanceA_gap :
    IsGreatest
      {gam : ℝ | ∀ f : Obs 2, piMean (ksWeight rA) f = 0 →
        gam * statePairing (prodWeight (ksWeight rA)) f f
          ≤ statePairing (prodWeight (ksWeight rA)) f
              (generator (ksWeight rA) (fun l => fiberRate 1 (rA l)) f)}
      (5 / 2) := by
  have hw := ksWeight_pos rA_pos
  have hs := ksWeight_sum rA
  have hc : ∀ l : Fin 2, 0 ≤ fiberRate 1 (rA l) := fun l =>
    le_trans (by norm_num) (two_mul_le_fiberRate one_pos (rA_pos l))
  constructor
  · refine (dirichlet_constant_iff hw hs hc _).mpr ?_
    intro l
    fin_cases l
    · rw [instanceA_rate_values.1]
    · rw [instanceA_rate_values.2]
      norm_num
  · intro gam hgam
    have h := rate_le_of_dirichlet_bound hw hs hgam 0
    rwa [instanceA_rate_values.1] at h

theorem instanceA_gap_above_floor : (2 : ℝ) * 1 < 5 / 2 := by norm_num

/-- Ratios of the three-fiber instance: `r = (1, 2, 3)`, `λ = 1/2`. -/
def rB : Fin 3 → ℝ := ![1, 2, 3]

theorem rB_pos : ∀ l, 0 < rB l := by
  intro l
  fin_cases l <;> norm_num [rB]

/-- Exact rates of the three-fiber instance: `(1, 5/4, 5/3)`. -/
theorem instanceB_rate_values :
    fiberRate (1 / 2) (rB 0) = 1 ∧ fiberRate (1 / 2) (rB 1) = 5 / 4 ∧
      fiberRate (1 / 2) (rB 2) = 5 / 3 := by
  refine ⟨?_, ?_, ?_⟩ <;> norm_num [fiberRate, rB]

/-- Exact full-level eigenvalue of the three-fiber instance:
`1 + 5/4 + 5/3 = 47/12`. -/
theorem instanceB_top_eigenvalue :
    (∑ l ∈ (Finset.univ : Finset (Fin 3)), fiberRate (1 / 2) (rB l))
      = 47 / 12 := by
  rw [Fin.sum_univ_three]
  norm_num [fiberRate, rB]

/-- Exact gap of the three-fiber instance: the best Dirichlet constant is
`1 = 2λ`; the floor is attained because the first ratio is `1`. -/
theorem instanceB_gap :
    IsGreatest
      {gam : ℝ | ∀ f : Obs 3, piMean (ksWeight rB) f = 0 →
        gam * statePairing (prodWeight (ksWeight rB)) f f
          ≤ statePairing (prodWeight (ksWeight rB)) f
              (generator (ksWeight rB) (fun l => fiberRate (1 / 2) (rB l)) f)}
      1 := by
  have hw := ksWeight_pos rB_pos
  have hs := ksWeight_sum rB
  have hlam : (0 : ℝ) < 1 / 2 := by norm_num
  have hc : ∀ l : Fin 3, 0 ≤ fiberRate (1 / 2) (rB l) := fun l =>
    le_trans (by norm_num) (two_mul_le_fiberRate hlam (rB_pos l))
  constructor
  · refine (dirichlet_constant_iff hw hs hc _).mpr ?_
    intro l
    fin_cases l
    · rw [instanceB_rate_values.1]
    · rw [instanceB_rate_values.2.1]
      norm_num
    · rw [instanceB_rate_values.2.2]
      norm_num
  · intro gam hgam
    have h := rate_le_of_dirichlet_bound hw hs hgam 0
    rwa [instanceB_rate_values.1] at h

/-- The three-fiber instance attains the `2λ` floor through its `r = 1`
fiber, in agreement with the scalar equality case. -/
theorem instanceB_attains_floor :
    fiberRate (1 / 2) (rB 0) = 2 * (1 / 2) ∧ rB 0 = 1 := by
  constructor
  · norm_num [fiberRate, rB]
  · norm_num [rB]

/-- The scalar equality case applied at the attaining fiber: for
`λ = 1/2`, `r = 1` is the only positive ratio whose rate meets `2λ`. -/
theorem instanceB_equality_case (r : ℝ) (hr : 0 < r) :
    fiberRate (1 / 2) r = 2 * (1 / 2) ↔ r = 1 :=
  fiberRate_eq_two_mul_iff (by norm_num) hr

end Instances

end

end OPH.KogutSusskindFiberRate

/- Axiom audit: propext, Classical.choice, and Quot.sound only. -/

#print axioms OPH.KogutSusskindFiberRate.fiberRate_flip
#print axioms OPH.KogutSusskindFiberRate.two_mul_le_fiberRate
#print axioms OPH.KogutSusskindFiberRate.fiberRate_eq_two_mul_iff
#print axioms OPH.KogutSusskindFiberRate.fiberRate_heatBath_offdiagonal
#print axioms OPH.KogutSusskindFiberRate.heatBath_offdiagonal_unique_rate
#print axioms OPH.KogutSusskindFiberRate.statePairing_rate_lower
#print axioms OPH.KogutSusskindFiberRate.statePairing_condExp_residual
#print axioms OPH.KogutSusskindFiberRate.statePairing_condExp_rate_residual
#print axioms OPH.KogutSusskindFiberRate.statePairing_residual_split
#print axioms OPH.KogutSusskindFiberRate.statePairing_rate_residual_split
#print axioms OPH.KogutSusskindFiberRate.dirichlet_domination_single
#print axioms OPH.KogutSusskindFiberRate.dirichlet_domination_family
#print axioms OPH.KogutSusskindFiberRate.twoPointCondExp_idem
#print axioms OPH.KogutSusskindFiberRate.twoPointCondExp_symm
#print axioms OPH.KogutSusskindFiberRate.twoPointCondExp_const_rate_comm
#print axioms OPH.KogutSusskindFiberRate.twoPoint_kogutSusskind_domination
#print axioms OPH.KogutSusskindFiberRate.condExp_walsh_of_mem
#print axioms OPH.KogutSusskindFiberRate.condExp_walsh_of_notMem
#print axioms OPH.KogutSusskindFiberRate.generator_walsh
#print axioms OPH.KogutSusskindFiberRate.statePairing_walsh_walsh
#print axioms OPH.KogutSusskindFiberRate.statePairing_walsh_orthogonal
#print axioms OPH.KogutSusskindFiberRate.statePairing_walsh_self
#print axioms OPH.KogutSusskindFiberRate.statePairing_walsh_self_pos
#print axioms OPH.KogutSusskindFiberRate.prodWeight_sum_eq_one
#print axioms OPH.KogutSusskindFiberRate.walsh_linearIndependent
#print axioms OPH.KogutSusskindFiberRate.walshBasis_apply
#print axioms OPH.KogutSusskindFiberRate.repr_expansion
#print axioms OPH.KogutSusskindFiberRate.statePairing_self_diagonal
#print axioms OPH.KogutSusskindFiberRate.statePairing_generator_diagonal
#print axioms OPH.KogutSusskindFiberRate.repr_empty_of_meanZero
#print axioms OPH.KogutSusskindFiberRate.dirichlet_lower_bound
#print axioms OPH.KogutSusskindFiberRate.rate_le_of_dirichlet_bound
#print axioms OPH.KogutSusskindFiberRate.dirichlet_constant_iff
#print axioms OPH.KogutSusskindFiberRate.bestConstant_isGreatest_min_rate
#print axioms OPH.KogutSusskindFiberRate.kogutSusskind_uniform_dirichlet_floor
#print axioms OPH.KogutSusskindFiberRate.ksWeight_pos
#print axioms OPH.KogutSusskindFiberRate.ksWeight_sum
#print axioms OPH.KogutSusskindFiberRate.ksWeight_ratio
#print axioms OPH.KogutSusskindFiberRate.kogutSusskind_ksWeight_floor
#print axioms OPH.KogutSusskindFiberRate.instanceA_rate_values
#print axioms OPH.KogutSusskindFiberRate.instanceA_eigenvalue_table
#print axioms OPH.KogutSusskindFiberRate.instanceA_gap
#print axioms OPH.KogutSusskindFiberRate.instanceB_rate_values
#print axioms OPH.KogutSusskindFiberRate.instanceB_top_eigenvalue
#print axioms OPH.KogutSusskindFiberRate.instanceB_gap
#print axioms OPH.KogutSusskindFiberRate.instanceB_attains_floor
#print axioms OPH.KogutSusskindFiberRate.instanceB_equality_case
