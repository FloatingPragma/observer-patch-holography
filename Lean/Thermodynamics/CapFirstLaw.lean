import FiniteConditionalRepair

namespace OPH.Thermodynamics

/-!
# The cap first law with exact remainder

The Einstein branch consumes the identity between the entropy
differential and the modular pairing as an explicit premise of its
finite bulk-edge-central first law. This module derives the finite,
exact form of that premise from the thermodynamic package: for any
state `p` and faithful reference `τ` with modular Hamiltonian
`K = -log τ`,

`S(p) - S(τ) = ⟨K⟩_p - ⟨K⟩_τ - D(p ‖ τ)`,

an exact identity whose deficit is the relative entropy itself, so the
first-order first law `δS = δ⟨K⟩` is this identity with the
second-order deficit dropped. On a central split `K = 2π B + Z` the
identity distributes over the two charges, and a `τ`-preserving
channel that conserves the mean central charge obeys the cap Clausius
inequality `2π Δ⟨B⟩ ≤ ΔS`. The conditional-resampling repair kernel
conserves the mean of every fibre-measurable charge, so a central
charge measurable through the repaired visible datum is conserved
automatically and the repair kernel satisfies the cap Clausius
inequality outright. The identification of the abstract split with the
physical cap data stays with the energy-clock receipt.
-/

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]
variable {B : Type*} [DecidableEq B]

/-- Relative entropy of a law to itself vanishes. -/
theorem kl_self (tau : Ω → ℝ) : kl tau tau = 0 := by
  unfold kl
  apply Finset.sum_eq_zero
  intro x _
  by_cases h : tau x = 0
  · simp [h]
  · rw [klTerm_of_pos h, div_self h, Real.log_one, mul_zero]

/-- The Shannon entropy of the faithful reference equals its own
modular energy: `S(τ) = ⟨K⟩_τ` with `K = -log τ`. -/
theorem shannon_ref_eq_modular_energy (tau : Ω → ℝ)
    (hτ : ∀ x, 0 < tau x) :
    shannon tau = ∑ x, tau x * (-Real.log (tau x)) := by
  have h := kl_eq_energy_sub_shannon tau tau hτ
  rw [kl_self] at h
  linarith

/-- **Cap first law with exact remainder.** The entropy difference to
the reference equals the modular-energy difference minus the relative
entropy: `S(p) - S(τ) = ⟨K⟩_p - ⟨K⟩_τ - D(p ‖ τ)`. The first-order
entanglement first law `δS = δ⟨K⟩` is this identity with the
quadratic-order deficit `D(p ‖ τ)` dropped. -/
theorem cap_firstLaw_exact (p tau : Ω → ℝ) (hτ : ∀ x, 0 < tau x) :
    shannon p - shannon tau
      = (∑ x, p x * (-Real.log (tau x)))
        - (∑ x, tau x * (-Real.log (tau x))) - kl p tau := by
  have h1 := kl_eq_energy_sub_shannon p tau hτ
  have h2 := shannon_ref_eq_modular_energy tau hτ
  linarith

/-- Distributing a central split `K = 2π B + Z` of the modular
Hamiltonian over an expectation. -/
theorem modular_energy_split (tau Bc Z : Ω → ℝ) (s : Ω → ℝ)
    (hsplit : ∀ x, -Real.log (tau x) = 2 * Real.pi * Bc x + Z x) :
    ∑ x, s x * (-Real.log (tau x))
      = 2 * Real.pi * (∑ x, s x * Bc x) + ∑ x, s x * Z x := by
  calc ∑ x, s x * (-Real.log (tau x))
      = ∑ x, (2 * Real.pi * (s x * Bc x) + s x * Z x) := by
        apply Finset.sum_congr rfl
        intro x _
        rw [hsplit x]
        ring
    _ = 2 * Real.pi * (∑ x, s x * Bc x) + ∑ x, s x * Z x := by
        rw [Finset.sum_add_distrib, ← Finset.mul_sum]

/-- The cap first law through the central split `K = 2π B + Z`:
`S(p) - S(τ) = 2π Δ⟨B⟩ + Δ⟨Z⟩ - D(p ‖ τ)` exactly. -/
theorem cap_firstLaw_split (p tau Bc Z : Ω → ℝ)
    (hτ : ∀ x, 0 < tau x)
    (hsplit : ∀ x, -Real.log (tau x) = 2 * Real.pi * Bc x + Z x) :
    shannon p - shannon tau
      = 2 * Real.pi * ((∑ x, p x * Bc x) - (∑ x, tau x * Bc x))
        + ((∑ x, p x * Z x) - (∑ x, tau x * Z x))
        - kl p tau := by
  have h := cap_firstLaw_exact p tau hτ
  rw [modular_energy_split tau Bc Z p hsplit,
    modular_energy_split tau Bc Z tau hsplit] at h
  linarith

/-- **Cap Clausius inequality.** A `τ`-preserving stochastic channel
that conserves the mean central charge satisfies
`2π Δ⟨B⟩ ≤ ΔS` for the central split `K = 2π B + Z`. -/
theorem cap_clausius_of_central_conserved (p tau Bc Z : Ω → ℝ)
    (K : Ω → Ω → ℝ)
    (hp : ∀ x, 0 ≤ p x) (hτ : ∀ x, 0 < tau x)
    (hK0 : ∀ x y, 0 ≤ K x y) (hK1 : ∀ x, ∑ y, K x y = 1)
    (hstat : ∀ y, push tau K y = tau y)
    (hsplit : ∀ x, -Real.log (tau x) = 2 * Real.pi * Bc x + Z x)
    (hZ : ∑ x, push p K x * Z x = ∑ x, p x * Z x) :
    2 * Real.pi * ((∑ x, push p K x * Bc x) - (∑ x, p x * Bc x))
      ≤ shannon (push p K) - shannon p := by
  have h := clausius p tau K hp hτ hK0 hK1 hstat
  rw [modular_energy_split tau Bc Z (push p K) hsplit,
    modular_energy_split tau Bc Z p hsplit] at h
  linarith

/-- The repair kernel has nonnegative entries. -/
theorem heatBath_nonneg {π : Ω → ℝ} {b : Ω → B} (hπ : ∀ x, 0 < π x)
    (x y : Ω) : 0 ≤ heatBath π b x y := by
  unfold heatBath
  by_cases h : b y = b x
  · rw [if_pos h]
    exact le_of_lt (div_pos (hπ y) (fiberMass_pos hπ x))
  · rw [if_neg h]

/-- The pushforward through the repair kernel preserves the mean of
every fibre-measurable observable: dynamical conservation of protected
charges at the level of expectations. -/
theorem push_heatBath_fixes_mean {π : Ω → ℝ} {b : Ω → B}
    (hπ : ∀ x, 0 < π x) (Q : Ω → ℝ)
    (hQ : ∀ x y, b x = b y → Q x = Q y) (p : Ω → ℝ) :
    ∑ y, push p (heatBath π b) y * Q y = ∑ x, p x * Q x := by
  unfold push
  calc ∑ y, (∑ x, p x * heatBath π b x y) * Q y
      = ∑ y, ∑ x, p x * (heatBath π b x y * Q y) := by
        apply Finset.sum_congr rfl
        intro y _
        rw [Finset.sum_mul]
        apply Finset.sum_congr rfl
        intro x _
        ring
    _ = ∑ x, ∑ y, p x * (heatBath π b x y * Q y) := Finset.sum_comm
    _ = ∑ x, p x * Q x := by
        apply Finset.sum_congr rfl
        intro x _
        rw [← Finset.mul_sum, heatBath_fixes_fiberObservable hπ Q hQ x]

/-- **Cap Clausius for conditional repair.** When the reference law
carries a central split `K = 2π B + Z` whose central charge is
measurable through the repaired visible datum, one repair step obeys
`2π Δ⟨B⟩ ≤ ΔS` with no further hypothesis: the repair kernel conserves
the central charge automatically and preserves the reference. -/
theorem heatBath_cap_clausius {π : Ω → ℝ} {b : Ω → B}
    (hπ : ∀ x, 0 < π x) (Bc Z : Ω → ℝ) (p : Ω → ℝ)
    (hp : ∀ x, 0 ≤ p x)
    (hsplit : ∀ x, -Real.log (π x) = 2 * Real.pi * Bc x + Z x)
    (hZ : ∀ x y, b x = b y → Z x = Z y) :
    2 * Real.pi * ((∑ x, push p (heatBath π b) x * Bc x)
        - (∑ x, p x * Bc x))
      ≤ shannon (push p (heatBath π b)) - shannon p :=
  cap_clausius_of_central_conserved p π Bc Z (heatBath π b) hp hπ
    (heatBath_nonneg hπ) (heatBath_row_sum hπ)
    (heatBath_stationary hπ) hsplit
    (push_heatBath_fixes_mean hπ Z hZ p)

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.kl_self
#print axioms OPH.Thermodynamics.cap_firstLaw_exact
#print axioms OPH.Thermodynamics.cap_firstLaw_split
#print axioms OPH.Thermodynamics.cap_clausius_of_central_conserved
#print axioms OPH.Thermodynamics.push_heatBath_fixes_mean
#print axioms OPH.Thermodynamics.heatBath_cap_clausius
