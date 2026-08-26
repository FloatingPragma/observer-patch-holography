import GibbsReferenceEnergyIdentification

set_option autoImplicit false

open scoped BigOperators

/-!
# Additivity of the modular-energy ledger and the binding defect
(issues 736, 739, 729)

STATUS.  Exact finite identities for the modular energy ledger
`modularEnergy p tau = ∑ p x * (-log tau x)` of
`Geometry/InternalEnergyInertia.lean` on a product of two finite systems,
against the Gibbs reference `gibbs H β` of
`Thermodynamics/GibbsReferenceEnergyIdentification.lean`.  The ledger of a
composite splits into the ledgers of the parts plus a defect; against an
interacting Gibbs reference the defect is `β` times the interaction
expectation plus a state-independent constant.  Nothing here derives the
identification of the ledger with a laboratory energy, of the reference with
the record reference, or of the interaction term with a physical
interaction; those are declared inputs, named where they enter.

WHAT IS PROVED.

(1) Product additivity (`modularEnergy_productRef`).  For a state `p` on
`Ω₁ × Ω₂` with marginals `marginal₁ p x = ∑ y, p (x, y)` and
`marginal₂ p y = ∑ x, p (x, y)`, and a product reference
`productRef tau₁ tau₂ (x, y) = tau₁ x * tau₂ y` with positive factors,
`modularEnergy p (productRef tau₁ tau₂) = modularEnergy (marginal₁ p) tau₁
+ modularEnergy (marginal₂ p) tau₂` exactly, for every `p`, correlated or
not.  Corollary for product states (`modularEnergy_productState`): the
marginals of `productState p₁ p₂` are the mass-scaled factors
(`marginal₁_productState`), so the ledger is `mass p₂ * E₁ + mass p₁ * E₂`,
which for normalized factors is `E₁ + E₂`.

(2) Binding defect (`defect`).  For a general joint reference `tau₁₂`,
`defect p tau₁₂ tau₁ tau₂ = modularEnergy p tau₁₂ - modularEnergy
(marginal₁ p) tau₁ - modularEnergy (marginal₂ p) tau₂`.  It equals
`∑ p (x, y) * log (tau₁ x * tau₂ y / tau₁₂ (x, y))`
(`defect_eq_sum_log_ratio`) and vanishes for the product reference
(`defect_productRef`).

(3) Gibbs interacting reference (`defect_gibbs`).  With
`tau₁₂ = gibbs (H₁ + H₂ + V) β` and `tau_i = gibbs H_i β`,
`defect = β * expectation p V + mass p * bindingConstant β H₁ H₂ V` where
`bindingConstant = log Z₁₂ - log Z₁ - log Z₂`; for a normalized state the
constant enters with coefficient one (`defect_gibbs_normalized`).  The
constant is the log of the ratio of partition functions
(`bindingConstant_eq_log_ratio`) and vanishes at `V = 0`
(`partitionZ_product`, `bindingConstant_zero`), so at `V = 0` the
defect vanishes (`defect_gibbs_noInteraction`).

(4) Sign (`state_part_nonpos`, `bindingConstant_nonneg`,
`defect_gibbs_le_const`, `defect_gibbs_constInteraction`).  For `V ≤ 0`
pointwise (a declared attractive interaction), `p ≥ 0`, and `β ≥ 0`, the
state-dependent part `β * expectation p V` is `≤ 0`, the binding constant
is `≥ 0`, and the defect of a normalized state is at most the constant.
The sign of the full defect is not fixed by these hypotheses: for a constant
interaction `V = -c` the defect vanishes for every state, and in the
two-point setting the uniform state has positive full defect
(`two_point_defect_pos`, value `-1/4 + log ((3 + e)/4)`, proved through
the bound `(3 + e)/4 > e^(1/4)`).  The two-point example
(`two_point_strict`) on `Bool × Bool` with uniform `p`, `β = 1`, and `V`
equal to `-1` on `(true, true)` has state-dependent part `-1/4 < 0`.

(5) Composition with the inertia precursor (`composite_inertial_ledger`).
Under the declared shape A of `Geometry/InternalEnergyInertia.lean`
(internal energy enters the inertial coefficient with slope one), the
inertial coefficient of the composite whose ledger is the joint modular
energy equals `m₁ + m₂ + E₁ + E₂ + defect`, a definitional identity.  The
shape is declared; the Legendre non-identifiability at its scope is cited
(`legendre_nonidentifiability_recited`): every Lagrangian shape is a
declared enrichment.

(6) Boundary (`defect_is_expectation`).  For every faithful joint
reference the defect is the expectation of the state-independent function
`effectiveInteraction tau₁₂ tau₁ tau₂ = log (tau₁ x * tau₂ y / tau₁₂ (x, y))`
in the state `p`; a normalized finite faithful reference is a Gibbs reference at `β = 1`
for its own modular Hamiltonian, so the defect is always an expectation of a
`p`-independent function.  What is not delivered by any theorem is the
identification of that function with a physical interaction.

ROWS TOUCHED.  Laboratory clock and energy calibration import (the
identification of the ledger with a laboratory energy stays with
`Thermodynamics/PhysicalCalibrationImport.lean`); gravitation-route energy
identification (the identification of `tau₁₂` with the record reference and
of `V` with a physical interaction); coupled-action row (shape A is a
declared composite action).  None discharged.

NEGATIVES CITED.  The Legendre no-go of
`Variational/RealizedHistoryLegendreNoGo.lean`, through
`legendre_nonidentifiability_cited`, at its scope: realized histories do
not select a velocity curvature or Legendre map, so shape A is declared.

CONVENTIONS.  Ledger in nats; `modularEnergy p tau = ∑ p x * (-log tau x)`;
`gibbs H β x = exp (-β * H x) / partitionZ H β`; marginals are
sums; states are arbitrary real functions unless a hypothesis names
nonnegativity or normalization; `defect` is joint minus parts.

FALSIFIER.  The module fails if the product additivity misses a term, if
the defect differs from the displayed log-ratio sum, if the Gibbs defect
differs from `β⟨V⟩ + mass * bindingConstant`, if the binding constant fails
to vanish at `V = 0`, or if the two-point example fails strictness.

Axiom audit.  No project axiom, no native decision procedure; the guard
lines at the end show at most `propext`, `Classical.choice`, `Quot.sound`.
-/

namespace OPH.ModularEnergyAdditivity

open OPH.Thermodynamics OPH.InternalEnergyInertia
open OPH.GibbsReferenceEnergyIdentification

noncomputable section

variable {Ω₁ Ω₂ : Type*} [Fintype Ω₁] [Fintype Ω₂]

/-! ## (1) Marginals, product references, product additivity -/

/-- First marginal `∑ y, p (x, y)`. -/
def marginal₁ (p : Ω₁ × Ω₂ → ℝ) (x : Ω₁) : ℝ := ∑ y, p (x, y)

/-- Second marginal `∑ x, p (x, y)`. -/
def marginal₂ (p : Ω₁ × Ω₂ → ℝ) (y : Ω₂) : ℝ := ∑ x, p (x, y)

/-- Product reference `tau₁ x * tau₂ y`. -/
def productRef (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ) (z : Ω₁ × Ω₂) : ℝ :=
  tau₁ z.1 * tau₂ z.2

/-- Product state `p₁ x * p₂ y`. -/
def productState (p₁ : Ω₁ → ℝ) (p₂ : Ω₂ → ℝ) (z : Ω₁ × Ω₂) : ℝ :=
  p₁ z.1 * p₂ z.2

omit [Fintype Ω₁] [Fintype Ω₂] in
theorem productRef_pos (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ)
    (h₁ : ∀ x, 0 < tau₁ x) (h₂ : ∀ y, 0 < tau₂ y) (z : Ω₁ × Ω₂) :
    0 < productRef tau₁ tau₂ z :=
  mul_pos (h₁ z.1) (h₂ z.2)

theorem mass_marginal₁ (p : Ω₁ × Ω₂ → ℝ) : mass (marginal₁ p) = mass p := by
  unfold mass marginal₁
  rw [Fintype.sum_prod_type]

theorem mass_marginal₂ (p : Ω₁ × Ω₂ → ℝ) : mass (marginal₂ p) = mass p := by
  unfold mass marginal₂
  rw [Fintype.sum_prod_type, Finset.sum_comm]

/-- A joint expectation of a function of the first coordinate is the
marginal expectation. -/
theorem expectation_fst (p : Ω₁ × Ω₂ → ℝ) (f : Ω₁ → ℝ) :
    expectation p (fun z => f z.1) = expectation (marginal₁ p) f := by
  unfold expectation marginal₁
  rw [Fintype.sum_prod_type]
  apply Finset.sum_congr rfl
  intro x _
  rw [Finset.sum_mul]

/-- A joint expectation of a function of the second coordinate is the
marginal expectation. -/
theorem expectation_snd (p : Ω₁ × Ω₂ → ℝ) (g : Ω₂ → ℝ) :
    expectation p (fun z => g z.2) = expectation (marginal₂ p) g := by
  unfold expectation marginal₂
  rw [Fintype.sum_prod_type, Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro y _
  rw [Finset.sum_mul]

theorem expectation_add (p f g : Ω₁ × Ω₂ → ℝ) :
    expectation p (fun z => f z + g z) = expectation p f + expectation p g := by
  unfold expectation
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro z _
  ring

/-- The ledger is the expectation of the modular Hamiltonian. -/
theorem modularEnergy_eq_expectation {Ω : Type*} [Fintype Ω] (p tau : Ω → ℝ) :
    modularEnergy p tau = expectation p (fun x => -Real.log (tau x)) := rfl

/-- **Product additivity.**  Against a product reference the ledger of any
joint state, correlated or not, is the sum of the marginal ledgers. -/
theorem modularEnergy_productRef (p : Ω₁ × Ω₂ → ℝ)
    (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ)
    (h₁ : ∀ x, 0 < tau₁ x) (h₂ : ∀ y, 0 < tau₂ y) :
    modularEnergy p (productRef tau₁ tau₂)
      = modularEnergy (marginal₁ p) tau₁ + modularEnergy (marginal₂ p) tau₂ := by
  rw [modularEnergy_eq_expectation, modularEnergy_eq_expectation,
    modularEnergy_eq_expectation, ← expectation_fst, ← expectation_snd,
    ← expectation_add]
  unfold expectation productRef
  apply Finset.sum_congr rfl
  intro z _
  show p z * (-Real.log (tau₁ z.1 * tau₂ z.2))
    = p z * (-Real.log (tau₁ z.1) + -Real.log (tau₂ z.2))
  rw [Real.log_mul (h₁ z.1).ne' (h₂ z.2).ne']
  ring

omit [Fintype Ω₁] in
theorem marginal₁_productState (p₁ : Ω₁ → ℝ) (p₂ : Ω₂ → ℝ) :
    marginal₁ (productState p₁ p₂) = fun x => mass p₂ * p₁ x := by
  funext x
  unfold marginal₁ productState mass
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro y _
  ring

omit [Fintype Ω₂] in
theorem marginal₂_productState (p₁ : Ω₁ → ℝ) (p₂ : Ω₂ → ℝ) :
    marginal₂ (productState p₁ p₂) = fun y => mass p₁ * p₂ y := by
  funext y
  unfold marginal₂ productState mass
  rw [Finset.sum_mul]

theorem modularEnergy_smul {Ω : Type*} [Fintype Ω] (c : ℝ) (p tau : Ω → ℝ) :
    modularEnergy (fun x => c * p x) tau = c * modularEnergy p tau := by
  unfold modularEnergy
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro x _
  ring

/-- **Product-state corollary.**  For a product state the ledger against a
product reference is `mass p₂ * E₁ + mass p₁ * E₂`. -/
theorem modularEnergy_productState (p₁ : Ω₁ → ℝ) (p₂ : Ω₂ → ℝ)
    (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ)
    (h₁ : ∀ x, 0 < tau₁ x) (h₂ : ∀ y, 0 < tau₂ y) :
    modularEnergy (productState p₁ p₂) (productRef tau₁ tau₂)
      = mass p₂ * modularEnergy p₁ tau₁ + mass p₁ * modularEnergy p₂ tau₂ := by
  rw [modularEnergy_productRef _ _ _ h₁ h₂, marginal₁_productState,
    marginal₂_productState, modularEnergy_smul, modularEnergy_smul]

/-- Normalized factors: the ledger of the product state is `E₁ + E₂`. -/
theorem modularEnergy_productState_normalized (p₁ : Ω₁ → ℝ) (p₂ : Ω₂ → ℝ)
    (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ)
    (h₁ : ∀ x, 0 < tau₁ x) (h₂ : ∀ y, 0 < tau₂ y)
    (n₁ : mass p₁ = 1) (n₂ : mass p₂ = 1) :
    modularEnergy (productState p₁ p₂) (productRef tau₁ tau₂)
      = modularEnergy p₁ tau₁ + modularEnergy p₂ tau₂ := by
  rw [modularEnergy_productState _ _ _ _ h₁ h₂, n₁, n₂, one_mul, one_mul]

/-! ## (2) The binding defect -/

/-- The binding defect: joint ledger minus the two marginal ledgers. -/
def defect (p : Ω₁ × Ω₂ → ℝ) (tau₁₂ : Ω₁ × Ω₂ → ℝ)
    (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ) : ℝ :=
  modularEnergy p tau₁₂ - modularEnergy (marginal₁ p) tau₁
    - modularEnergy (marginal₂ p) tau₂

/-- The state-independent effective interaction of a joint reference
relative to the product of the part references. -/
def effectiveInteraction (tau₁₂ : Ω₁ × Ω₂ → ℝ) (tau₁ : Ω₁ → ℝ)
    (tau₂ : Ω₂ → ℝ) (z : Ω₁ × Ω₂) : ℝ :=
  Real.log (tau₁ z.1 * tau₂ z.2 / tau₁₂ z)

/-- **Defect as a log-ratio sum.**  For faithful references the defect is
`∑ p (x, y) * log (tau₁ x * tau₂ y / tau₁₂ (x, y))`. -/
theorem defect_eq_sum_log_ratio (p : Ω₁ × Ω₂ → ℝ) (tau₁₂ : Ω₁ × Ω₂ → ℝ)
    (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ)
    (h₁₂ : ∀ z, 0 < tau₁₂ z) (h₁ : ∀ x, 0 < tau₁ x) (h₂ : ∀ y, 0 < tau₂ y) :
    defect p tau₁₂ tau₁ tau₂
      = ∑ z, p z * Real.log (tau₁ z.1 * tau₂ z.2 / tau₁₂ z) := by
  unfold defect
  rw [sub_sub, ← modularEnergy_productRef p tau₁ tau₂ h₁ h₂]
  unfold modularEnergy productRef
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro z _
  rw [Real.log_div (mul_pos (h₁ z.1) (h₂ z.2)).ne' (h₁₂ z).ne']
  ring

/-- The defect vanishes for the product reference. -/
theorem defect_productRef (p : Ω₁ × Ω₂ → ℝ) (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ)
    (h₁ : ∀ x, 0 < tau₁ x) (h₂ : ∀ y, 0 < tau₂ y) :
    defect p (productRef tau₁ tau₂) tau₁ tau₂ = 0 := by
  unfold defect
  rw [modularEnergy_productRef p tau₁ tau₂ h₁ h₂]
  ring


/-! ## (3) The Gibbs interacting reference -/

section Gibbs

variable [Nonempty Ω₁] [Nonempty Ω₂] [DecidableEq Ω₁] [DecidableEq Ω₂]

/-- The declared total energy `H₁ + H₂ + V` of the interacting composite. -/
def totalEnergy (H₁ : Ω₁ → ℝ) (H₂ : Ω₂ → ℝ) (V : Ω₁ × Ω₂ → ℝ)
    (z : Ω₁ × Ω₂) : ℝ :=
  H₁ z.1 + H₂ z.2 + V z

/-- The state-independent binding constant `log Z₁₂ - log Z₁ - log Z₂`. -/
def bindingConstant (β : ℝ) (H₁ : Ω₁ → ℝ) (H₂ : Ω₂ → ℝ)
    (V : Ω₁ × Ω₂ → ℝ) : ℝ :=
  Real.log (partitionZ (totalEnergy H₁ H₂ V) β)
    - Real.log (partitionZ H₁ β) - Real.log (partitionZ H₂ β)

omit [Nonempty Ω₁] [Nonempty Ω₂] [DecidableEq Ω₁] [DecidableEq Ω₂] in
theorem expectation_totalEnergy (p : Ω₁ × Ω₂ → ℝ) (H₁ : Ω₁ → ℝ)
    (H₂ : Ω₂ → ℝ) (V : Ω₁ × Ω₂ → ℝ) :
    expectation p (totalEnergy H₁ H₂ V)
      = expectation (marginal₁ p) H₁ + expectation (marginal₂ p) H₂
        + expectation p V := by
  rw [← expectation_fst, ← expectation_snd]
  unfold expectation totalEnergy
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro z _
  ring

/-- **Gibbs defect.**  Against `gibbs (H₁ + H₂ + V) β` with part references
`gibbs H_i β`, the defect is `β * ⟨V⟩_p + (∑ p) * bindingConstant`. -/
theorem defect_gibbs (β : ℝ) (p : Ω₁ × Ω₂ → ℝ) (H₁ : Ω₁ → ℝ) (H₂ : Ω₂ → ℝ)
    (V : Ω₁ × Ω₂ → ℝ) :
    defect p (gibbs (totalEnergy H₁ H₂ V) β) (gibbs H₁ β) (gibbs H₂ β)
      = β * expectation p V + mass p * bindingConstant β H₁ H₂ V := by
  unfold defect bindingConstant
  rw [modularEnergy_gibbs, modularEnergy_gibbs, modularEnergy_gibbs,
    expectation_totalEnergy, mass_marginal₁, mass_marginal₂]
  ring

/-- Normalized state: the defect is `β * ⟨V⟩_p + bindingConstant`. -/
theorem defect_gibbs_normalized (β : ℝ) (p : Ω₁ × Ω₂ → ℝ) (H₁ : Ω₁ → ℝ)
    (H₂ : Ω₂ → ℝ) (V : Ω₁ × Ω₂ → ℝ) (hp : mass p = 1) :
    defect p (gibbs (totalEnergy H₁ H₂ V) β) (gibbs H₁ β) (gibbs H₂ β)
      = β * expectation p V + bindingConstant β H₁ H₂ V := by
  rw [defect_gibbs, hp, one_mul]

/-- The binding constant is the log of the partition-function ratio. -/
theorem bindingConstant_eq_log_ratio (β : ℝ) (H₁ : Ω₁ → ℝ) (H₂ : Ω₂ → ℝ)
    (V : Ω₁ × Ω₂ → ℝ) :
    bindingConstant β H₁ H₂ V
      = Real.log (partitionZ (totalEnergy H₁ H₂ V) β
          / (partitionZ H₁ β * partitionZ H₂ β)) := by
  unfold bindingConstant
  rw [Real.log_div (partitionZ_pos _ _).ne'
      (mul_pos (partitionZ_pos _ _) (partitionZ_pos _ _)).ne',
    Real.log_mul (partitionZ_pos _ _).ne' (partitionZ_pos _ _).ne']
  ring

omit [Nonempty Ω₁] [Nonempty Ω₂] [DecidableEq Ω₁] [DecidableEq Ω₂] in
/-- Without interaction the partition function factorizes. -/
theorem partitionZ_product (β : ℝ) (H₁ : Ω₁ → ℝ) (H₂ : Ω₂ → ℝ) :
    partitionZ (totalEnergy H₁ H₂ (fun _ => 0)) β
      = partitionZ H₁ β * partitionZ H₂ β := by
  unfold partitionZ gibbsWeight totalEnergy
  rw [Finset.sum_mul_sum, Fintype.sum_prod_type]
  apply Finset.sum_congr rfl
  intro x _
  apply Finset.sum_congr rfl
  intro y _
  rw [← Real.exp_add]
  congr 1
  ring

/-- The binding constant vanishes at `V = 0`. -/
theorem bindingConstant_zero (β : ℝ) (H₁ : Ω₁ → ℝ) (H₂ : Ω₂ → ℝ) :
    bindingConstant β H₁ H₂ (fun _ => 0) = 0 := by
  unfold bindingConstant
  rw [partitionZ_product,
    Real.log_mul (partitionZ_pos _ _).ne' (partitionZ_pos _ _).ne']
  ring

/-- At `V = 0` the Gibbs defect vanishes for every state. -/
theorem defect_gibbs_noInteraction (β : ℝ) (p : Ω₁ × Ω₂ → ℝ)
    (H₁ : Ω₁ → ℝ) (H₂ : Ω₂ → ℝ) :
    defect p (gibbs (totalEnergy H₁ H₂ (fun _ => 0)) β) (gibbs H₁ β)
      (gibbs H₂ β) = 0 := by
  rw [defect_gibbs, bindingConstant_zero]
  unfold expectation
  simp

/-! ## (4) Sign for a declared attractive interaction -/

omit [Nonempty Ω₁] [Nonempty Ω₂] in
omit [Nonempty Ω₁] [Nonempty Ω₂] [DecidableEq Ω₁] [DecidableEq Ω₂] in
/-- The state-dependent part `β * ⟨V⟩_p` is nonpositive for `V ≤ 0`,
`p ≥ 0`, `β ≥ 0`. -/
theorem state_part_nonpos (β : ℝ) (p V : Ω₁ × Ω₂ → ℝ)
    (hβ : 0 ≤ β) (hp : ∀ z, 0 ≤ p z) (hV : ∀ z, V z ≤ 0) :
    β * expectation p V ≤ 0 := by
  apply mul_nonpos_of_nonneg_of_nonpos hβ
  unfold expectation
  exact Finset.sum_nonpos fun z _ => mul_nonpos_of_nonneg_of_nonpos (hp z) (hV z)

/-- **Binding constant is nonnegative.**  For `V ≤ 0` and `β ≥ 0`,
`Z₁₂ ≥ Z₁ Z₂`, so `bindingConstant ≥ 0`. -/
theorem bindingConstant_nonneg (β : ℝ) (H₁ : Ω₁ → ℝ) (H₂ : Ω₂ → ℝ)
    (V : Ω₁ × Ω₂ → ℝ) (hβ : 0 ≤ β) (hV : ∀ z, V z ≤ 0) :
    0 ≤ bindingConstant β H₁ H₂ V := by
  rw [bindingConstant_eq_log_ratio, ← partitionZ_product β H₁ H₂]
  apply Real.log_nonneg
  rw [le_div_iff₀ (partitionZ_pos _ _), one_mul]
  unfold partitionZ gibbsWeight
  apply Finset.sum_le_sum
  intro z _
  apply Real.exp_le_exp.mpr
  unfold totalEnergy
  nlinarith [hV z]

/-- **Upper bound by the binding constant.**  For a declared attractive
interaction the defect of a normalized nonnegative state is at most the
state-independent binding constant; the constant is nonnegative
(`bindingConstant_nonneg`), so the sign of the defect is not fixed. -/
theorem defect_gibbs_le_const (β : ℝ) (p : Ω₁ × Ω₂ → ℝ) (H₁ : Ω₁ → ℝ)
    (H₂ : Ω₂ → ℝ) (V : Ω₁ × Ω₂ → ℝ)
    (hβ : 0 ≤ β) (hp : ∀ z, 0 ≤ p z) (hmass : mass p = 1)
    (hV : ∀ z, V z ≤ 0) :
    defect p (gibbs (totalEnergy H₁ H₂ V) β) (gibbs H₁ β) (gibbs H₂ β)
      ≤ bindingConstant β H₁ H₂ V := by
  rw [defect_gibbs_normalized _ _ _ _ _ hmass]
  linarith [state_part_nonpos β p V hβ hp hV]

/-- **Constant interaction.**  For `V = -c` (attractive for `c ≥ 0`) the
joint Gibbs reference is the product reference, so the defect vanishes for
every state: the full defect can be zero with `V < 0`. -/
theorem defect_gibbs_constInteraction (β c : ℝ) (p : Ω₁ × Ω₂ → ℝ)
    (H₁ : Ω₁ → ℝ) (H₂ : Ω₂ → ℝ) :
    defect p (gibbs (totalEnergy H₁ H₂ (fun _ => -c)) β) (gibbs H₁ β)
      (gibbs H₂ β) = 0 := by
  have h : gibbs (totalEnergy H₁ H₂ (fun _ => -c)) β
      = gibbs (totalEnergy H₁ H₂ (fun _ => 0)) β := by
    have := gibbs_add_const (totalEnergy H₁ H₂ (fun _ => 0)) β (-c)
    rw [← this]
    congr 1
    funext z
    simp [totalEnergy]
  rw [h]
  exact defect_gibbs_noInteraction β p H₁ H₂

end Gibbs

/-- Two-point interaction on `Bool × Bool`: `-1` on `(true, true)`. -/
def twoPointV (z : Bool × Bool) : ℝ := if z = (true, true) then -1 else 0

/-- **Two-point strict example.**  Uniform state, `β = 1`: the
state-dependent part is `-1/4 < 0`. -/
theorem two_point_strict :
    (1 : ℝ) * expectation (fun _ : Bool × Bool => (1 / 4 : ℝ)) twoPointV
      = -1 / 4 ∧
    (1 : ℝ) * expectation (fun _ : Bool × Bool => (1 / 4 : ℝ)) twoPointV
      < 0 := by
  have h : (1 : ℝ) * expectation (fun _ : Bool × Bool => (1 / 4 : ℝ)) twoPointV
      = -1 / 4 := by
    unfold expectation twoPointV
    simp
    norm_num
  exact ⟨h, by rw [h]; norm_num⟩

/-- **Two-point full defect.**  In the same example, against the Gibbs
references with `H₁ = H₂ = 0`, `β = 1`, and interaction `twoPointV`, the
full defect of the uniform state is `-1/4 + log ((3 + e) / 4)`, and this is
positive: the binding constant outweighs the state-dependent part. -/
theorem two_point_defect_pos :
    defect (fun _ : Bool × Bool => (1 / 4 : ℝ))
      (gibbs (totalEnergy (fun _ : Bool => (0 : ℝ)) (fun _ => 0) twoPointV) 1)
      (gibbs (fun _ : Bool => (0 : ℝ)) 1) (gibbs (fun _ : Bool => (0 : ℝ)) 1)
      = -1 / 4 + Real.log ((3 + Real.exp 1) / 4) ∧
    0 < defect (fun _ : Bool × Bool => (1 / 4 : ℝ))
      (gibbs (totalEnergy (fun _ : Bool => (0 : ℝ)) (fun _ => 0) twoPointV) 1)
      (gibbs (fun _ : Bool => (0 : ℝ)) 1) (gibbs (fun _ : Bool => (0 : ℝ)) 1) := by
  have hmass : mass (fun _ : Bool × Bool => (1 / 4 : ℝ)) = 1 := by
    unfold mass
    simp
  have hZ12 : partitionZ (totalEnergy (fun _ : Bool => (0 : ℝ)) (fun _ => 0)
      twoPointV) 1 = 3 + Real.exp 1 := by
    unfold partitionZ gibbsWeight totalEnergy twoPointV
    simp [Fintype.sum_prod_type]
    ring
  have hZ1 : partitionZ (fun _ : Bool => (0 : ℝ)) 1 = 2 := by
    unfold partitionZ gibbsWeight
    simp
  have heq : defect (fun _ : Bool × Bool => (1 / 4 : ℝ))
      (gibbs (totalEnergy (fun _ : Bool => (0 : ℝ)) (fun _ => 0) twoPointV) 1)
      (gibbs (fun _ : Bool => (0 : ℝ)) 1) (gibbs (fun _ : Bool => (0 : ℝ)) 1)
      = -1 / 4 + Real.log ((3 + Real.exp 1) / 4) := by
    rw [defect_gibbs_normalized _ _ _ _ _ hmass, two_point_strict.1,
      bindingConstant_eq_log_ratio, hZ12, hZ1]
    norm_num
  refine ⟨heq, ?_⟩
  rw [heq]
  have hpos : (0 : ℝ) < (3 + Real.exp 1) / 4 := by positivity
  have hlt : Real.exp (1 / 4) < (3 + Real.exp 1) / 4 := by
    apply lt_of_pow_lt_pow_left₀ 4 hpos.le
    rw [← Real.exp_nat_mul]
    have he1 := Real.exp_one_gt_d9
    have he2 := Real.exp_one_lt_d9
    have : Real.exp ((4 : ℕ) * (1 / 4)) = Real.exp 1 := by norm_num
    rw [this]
    have h1 : (1.42 : ℝ) < (3 + Real.exp 1) / 4 := by linarith
    calc Real.exp 1 < 2.7182818286 := he2
      _ < (1.42 : ℝ) ^ 4 := by norm_num
      _ < ((3 + Real.exp 1) / 4) ^ 4 :=
          pow_lt_pow_left₀ h1 (by norm_num) (by norm_num)
  have := (Real.lt_log_iff_exp_lt hpos).2 hlt
  linarith

/-! ## (5) Composition with the inertia precursor (shape A, declared) -/

/-- **Composite inertial ledger.**  Under declared shape A the inertial
coefficient of the composite whose ledger is the joint modular energy is
`m₁ + m₂ + E₁ + E₂ + defect`, with `E_i` the marginal ledgers; a
definitional identity. -/
theorem composite_inertial_ledger (m₁ m₂ : ℝ) (p : Ω₁ × Ω₂ → ℝ)
    (tau₁₂ : Ω₁ × Ω₂ → ℝ) (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ) :
    inertialCoefficient (m₁ + m₂) (modularEnergy p tau₁₂)
      = m₁ + m₂ + modularEnergy (marginal₁ p) tau₁
        + modularEnergy (marginal₂ p) tau₂ + defect p tau₁₂ tau₁ tau₂ := by
  unfold inertialCoefficient defect
  ring

/-- The shape-A equation of motion of the composite, with the coefficient
written through the parts and the defect. -/
theorem composite_stationary_iff_ledger (m₁ m₂ : ℝ) (p : Ω₁ × Ω₂ → ℝ)
    (tau₁₂ : Ω₁ × Ω₂ → ℝ) (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ) (N : ℕ)
    (F x : ℕ → OPH.C1Lorentz.Herm2) :
    CompositeStationary (m₁ + m₂) (modularEnergy p tau₁₂) (N + 1) F x ↔
      ∀ j, j < N →
        (m₁ + m₂ + modularEnergy (marginal₁ p) tau₁
          + modularEnergy (marginal₂ p) tau₂ + defect p tau₁₂ tau₁ tau₂)
          • secondDifference x j = impulse F j := by
  rw [← composite_inertial_ledger m₁ m₂ p tau₁₂ tau₁ tau₂]
  exact (composite_inertial_coefficient (m₁ + m₂) (modularEnergy p tau₁₂)
    N F x).1

/-- The Legendre non-identifiability, re-cited at its scope: shape A is a
declared enrichment. -/
theorem legendre_nonidentifiability_recited :
    (¬ ∃ vel : ℝ → ℝ → ℝ, OPH.Variational.SolvesMomentum
      OPH.Variational.chainLogLagrangian vel) ∧
      OPH.Variational.chainCurvedLagrangian 1
        ≠ OPH.Variational.chainCurvedLagrangian 2 :=
  legendre_nonidentifiability_cited

/-! ## (6) Boundary: the defect is always an expectation -/

/-- **Boundary.**  For faithful references the defect is the expectation
of the state-independent effective interaction; a normalized finite
faithful joint reference is a Gibbs reference at `β = 1` for its own modular
Hamiltonian (`faithful_is_gibbs`), so no normalized non-Gibbs joint reference
exists at this scope.  The identification
of `effectiveInteraction` with a physical interaction is declared. -/
theorem defect_is_expectation (tau₁₂ : Ω₁ × Ω₂ → ℝ)
    (tau₁ : Ω₁ → ℝ) (tau₂ : Ω₂ → ℝ)
    (h₁₂ : ∀ z, 0 < tau₁₂ z) (h₁ : ∀ x, 0 < tau₁ x) (h₂ : ∀ y, 0 < tau₂ y) :
    ∀ p, defect p tau₁₂ tau₁ tau₂
      = expectation p (effectiveInteraction tau₁₂ tau₁ tau₂) := fun p =>
  defect_eq_sum_log_ratio p tau₁₂ tau₁ tau₂ h₁₂ h₁ h₂

/-- A normalized faithful reference is its own Gibbs reference at `β = 1`. -/
theorem faithful_is_gibbs {Ω : Type*} [Fintype Ω] [Nonempty Ω] (tau : Ω → ℝ)
    (hτ : ∀ x, 0 < tau x) (hn : ∑ x, tau x = 1) :
    OPH.Thermodynamics.gibbs (fun x => -Real.log (tau x)) 1
      = tau :=
  OPH.GibbsReferenceEnergyIdentification.gibbs_of_neg_log tau hτ hn

end

end OPH.ModularEnergyAdditivity

#print axioms OPH.ModularEnergyAdditivity.modularEnergy_productRef
#print axioms OPH.ModularEnergyAdditivity.modularEnergy_productState_normalized
#print axioms OPH.ModularEnergyAdditivity.defect_eq_sum_log_ratio
#print axioms OPH.ModularEnergyAdditivity.defect_productRef
#print axioms OPH.ModularEnergyAdditivity.defect_gibbs
#print axioms OPH.ModularEnergyAdditivity.defect_gibbs_normalized
#print axioms OPH.ModularEnergyAdditivity.bindingConstant_eq_log_ratio
#print axioms OPH.ModularEnergyAdditivity.bindingConstant_zero
#print axioms OPH.ModularEnergyAdditivity.defect_gibbs_noInteraction
#print axioms OPH.ModularEnergyAdditivity.state_part_nonpos
#print axioms OPH.ModularEnergyAdditivity.bindingConstant_nonneg
#print axioms OPH.ModularEnergyAdditivity.defect_gibbs_le_const
#print axioms OPH.ModularEnergyAdditivity.defect_gibbs_constInteraction
#print axioms OPH.ModularEnergyAdditivity.two_point_strict
#print axioms OPH.ModularEnergyAdditivity.two_point_defect_pos
#print axioms OPH.ModularEnergyAdditivity.composite_inertial_ledger
#print axioms OPH.ModularEnergyAdditivity.composite_stationary_iff_ledger
#print axioms OPH.ModularEnergyAdditivity.legendre_nonidentifiability_recited
#print axioms OPH.ModularEnergyAdditivity.defect_is_expectation
#print axioms OPH.ModularEnergyAdditivity.faithful_is_gibbs
