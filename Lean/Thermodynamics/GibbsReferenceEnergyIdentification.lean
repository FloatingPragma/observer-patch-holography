import Geometry.InternalEnergyInertia
import FourLawAdequacySurface

set_option autoImplicit false

open scoped BigOperators

/-!
# Gibbs reference: the modular energy ledger as an affine function of the
screen field energy (gravitation row issue 729; internal energy ledger of
issues 736 and 739)

STATUS.  Conditional identification, read against a declared register row.
The Einstein composition consumes the modular energy `modularEnergy p tau =
∑ p x * (-log tau x)` of `Geometry/InternalEnergyInertia.lean`; the screen
modules supply a field energy `fieldEnergyScaled`.  The corpus carries the
identification of the modular weight with a physical energy as the declared
clock-and-energy calibration register row
`FourLawAdequacySurface.CalibrationData` (field `thermal : -log ref = beta *
energy + logZ`, with `ref_eq_gibbs` giving the Gibbs form).  This module
supplies an inhabitant of that row built from the screen field energy
(`TemperatureDeclaration.toCalibrationData` with `energy = screenFieldEnergy
h`), reads the ledger against it (affine law, Clausius form, slope), and
bounds its scope (boundary example).  Nothing here derives that the record
reference is Gibbs, and nothing selects `β`; both are declared clauses.

PRIOR ART.  `Thermodynamics/FiniteConditionalRepair.lean` defines
`gibbsWeight`, `partitionZ`, `partitionZ_pos` and `gibbs E beta`, and proves
`gibbs_beta_injective` (two Gibbs states for one energy with two distinct
levels that agree as distributions share `beta`).  `Thermodynamics/
CapFirstLaw.lean` proves `modular_energy_split` (a pointwise split of the
modular Hamiltonian distributes over the ledger).  `Thermodynamics/
FourLawAdequacySurface.lean` defines `CalibrationData` and proves
`ref_eq_gibbs`, `zeroth_thermometer`, `zeroth_multiplier_unique`.  This
module uses those objects and defines no second Gibbs state:
`modularEnergy_gibbs` is the constant-offset instance of
`modular_energy_split`, and `screen_beta_unique` is `gibbs_beta_injective`
applied to the screen instance.

WHAT IS PROVED.

(1) Gibbs reference.  For a finite nonempty `Ω`, an energy `H : Ω → ℝ` and
`β : ℝ`, the corpus state `gibbs H β x = exp (-β * H x) / partitionZ H β` is
positive (`gibbs_pos`), normalized (`gibbs_sum`), and its modular
Hamiltonian is affine in `H` pointwise (`neg_log_gibbs`, the `thermal`
clause of `CalibrationData` with `logZ = log (partitionZ H β)`).  Exact
affine law (`modularEnergy_gibbs`, through `modular_energy_split`):
`modularEnergy p (gibbs H β) = β * expectation p H + mass p * log (partitionZ
H β)`, where `expectation p H = ∑ p x * H x` and `mass p = ∑ p x`; for a
normalized state the mass factor is one (`modularEnergy_gibbs_normalized`).
The ledger is the definition `modularEnergy` of
`Geometry/InternalEnergyInertia.lean`, imported and unfolded, with no
restatement.

(2) Clausius form of the cap first law (`cap_firstLaw_gibbs`): for a
normalized `p`, `shannon p - shannon (gibbs H β) = β * (expectation p H -
expectation (gibbs H β) H) - kl p (gibbs H β)`; the offset `log Z` cancels
exactly.  Between two normalized states the ledger difference is `β` times
the energy expectation difference (`modularEnergy_gibbs_sub`).

(3) Slope and reference uniqueness.  If two Gibbs references for the same
`H` at `β` and `β'` give equal ledgers on two normalized states with
different energy expectations, then `β = β'` (`slope_unique`, a ledger-side
form of `gibbs_beta_injective`); the ledger per unit energy expectation
between two such states is `β` (`slope_eq_beta`).  A faithful reference is
determined by its ledger on all states (`reference_determined_by_ledger`,
via point masses).  The additive freedom of `H` is exact in both
directions: `gibbs (fun x => H x + c) β = gibbs H β` (`gibbs_add_const`),
and two energies with one Gibbs reference at one nonzero `β` differ by a
constant (`gibbs_energy_unique_up_to_const`, with `c = (log Z - log Z') /
β`).  `TemperatureDeclaration`
packages the declared one-parameter identification: one positive number
`beta`; `gibbsRepairLaw H d` is the repair-law datum with reference `gibbs H
d.beta`, and `TemperatureDeclaration.toCalibrationData H d` is the
inhabitant of the `CalibrationData` row it induces, with `energy = H` and
`logZ = log (partitionZ H d.beta)`.

(4) Screen field instance.  `seedPotential e` is the declared potential
history that is zero at step `0` and equals `(e + 1) • Pi.single e 1` at
every later step (the `e`-th seam basis vector scaled by `e + 1`), with
zero port potential; `screenFieldEnergy h e = fieldEnergyScaled h
(seedPotential e) 0 0`.  `screenFieldEnergy_eq` computes it as
`(1 / 2) * h⁻¹ ^ 2 * (e + 1) ^ 2`.  The affine law holds for this instance
(`screen_affine_law`), the Clausius form holds (`screen_cap_firstLaw`),
`screenCalibration h d` is the `CalibrationData` inhabitant with `energy =
screenFieldEnergy h`, and once the family (with its declared `(e + 1)`
scaling) and the nonzero step `h` are declared, `β` is determined by the
reference: two values of `β` giving the same screen Gibbs reference are
equal (`screen_beta_unique`, from `gibbs_beta_injective` with levels `0`
and `1`).

(5) Boundary, a constant-energy example.  On `Fin 2` the reference `![1,
1 / 2]` is not Gibbs for the constant energy `0` at any `β`
(`boundary_not_gibbs`; the exclusion goes through the non-uniformity of
the reference, since a constant energy has the uniform Gibbs state), and
the two point masses have equal energy expectation and different ledgers
(`boundary_equal_expectation`, `boundary_ledger_ne`), so the ledger is not
an affine function of the energy expectation on normalized states
(`boundary_not_affine`).  The identification is specific to the Gibbs
clause for the declared `H`: the clause is the pair `(H, β)`, and `H`
carries the content, since every faithful normalized reference is the
Gibbs state of its own modular Hamiltonian at `β = 1`
(`gibbs_of_neg_log`: `gibbs (fun x => -log (tau x)) 1 = tau` under
positivity and `∑ tau = 1`).  A normalized non-Gibbs reference for a free
`H` therefore does not exist; the boundary example is excluded because a
constant energy has the uniform Gibbs state and the example is not uniform.

(6) Temperature dictionary.  `labInverseTemperature a cal β = β /
hbarOverTau a cal` is the laboratory inverse temperature in joules⁻¹ for
the declared tick and `thermalEnergyJoules a cal β = hbarOverTau a cal / β`
the laboratory thermal energy in joules; `beta_energy_invariant` proves `β
* E = labInverseTemperature a cal β * energyJoules a cal E` exactly,
`labInverseTemperature_mul_thermalEnergy` proves the two readings are
reciprocal for `β ≠ 0`, and `labInverseTemperature_not_forced` proves that
two distinct ticks give two distinct readings of one nonzero `β`, in the
pattern of `composite_joules_not_forced`.  The ledger `β` is dimensionless
per ledger unit; the laboratory reading is a conditional conversion.

ROWS TOUCHED.  The clock-and-energy calibration register row
`FourLawAdequacySurface.CalibrationData.thermal` (declared there; this
module builds an inhabitant from `(β, screenFieldEnergy h)` and reads the
ledger against it; the row stays declared).  The laboratory clock and
energy calibration import (the tick is declared; part (6) converts through
it and forces nothing).  The gravitation-route energy identification (the
identification of the modular energy with the screen field energy is
conditional on the calibration row inhabited by the screen field energy at
one declared `β`; it is not derived).  The source clock and duration row,
the physical spacetime attachment row, and the coupled-action row are
named and untouched: the step `h` is the declared scaled step, no carrier
map is used, and no worldline-field coupling appears.  This module
discharges none of these rows.

NEGATIVES CITED.  The Legendre non-identifiability
(`Variational/RealizedHistoryLegendreNoGo.lean`): realized source histories
select no velocity curvature or Legendre map, so every Lagrangian shape is
a declared enrichment; it is cited because the Gibbs clause is likewise a
declared enrichment of the reference, and nothing here reads `β` or the
reference off a history.  The abstract rate non-identifiability
(`ObserverPatchHolography/RateBridgeObstruction.lean`): its repair-layer
bridge fails, and it forbids nothing about a source-hosted process; no rate
claim is made here.

CONVENTIONS.  Forward differences in the step index; the electric field is
`-(h⁻¹ • (A (n + 1) - A n)) - realCoboundary (φ n)`, so the potential enters
with the negative orientation of `ScaledMaxwellStability`.  The Gibbs
exponent is `-β * H x` (energy lowers weight for `β > 0`), the argument
order is the corpus order `gibbs H β`, the ledger is in nats, and `β` is
per ledger energy unit.

FALSIFIER.  The module fails if the corpus `gibbs` is not positive or not
normalized, if the affine law misses a term, if the offset fails to cancel
in the Clausius form, if the slope is not unique, if two energies with one
Gibbs reference at nonzero `β` fail to differ by a constant, if a faithful
normalized reference is not Gibbs for its own modular Hamiltonian at `β =
1`, if the screen instance
energy differs from `(1 / 2) * h⁻¹ ^ 2 * (e + 1) ^ 2`, if the screen
calibration datum fails the `thermal` clause, if the boundary reference is
of Gibbs form, or if two ticks give one laboratory reading.

Axiom audit.  No project axiom, no native decision procedure; the guard
lines at the end show at most `propext`, `Classical.choice`, `Quot.sound`.
-/

namespace OPH.GibbsReferenceEnergyIdentification

open OPH.Thermodynamics OPH.InternalEnergyInertia
open OPH.PhysicalCalibrationImport OPH.ScaledMaxwellStability
open OPH.DiscreteCoulombGreen OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution

noncomputable section

/-! ## (1) The corpus Gibbs reference and the affine law -/

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]

/-- Energy expectation `∑ p x * H x` of a state `p`. -/
def expectation (p H : Ω → ℝ) : ℝ := ∑ x, p x * H x

/-- Total mass `∑ p x` of a state. -/
def mass (p : Ω → ℝ) : ℝ := ∑ x, p x

/-- The corpus Gibbs reference `gibbs H β` of
`Thermodynamics/FiniteConditionalRepair.lean` is faithful. -/
theorem gibbs_pos [Nonempty Ω] (H : Ω → ℝ) (β : ℝ) (x : Ω) :
    0 < gibbs H β x :=
  div_pos (Real.exp_pos _) (partitionZ_pos H β)

/-- The corpus Gibbs reference is normalized. -/
theorem gibbs_sum [Nonempty Ω] (H : Ω → ℝ) (β : ℝ) :
    ∑ x, gibbs H β x = 1 := by
  unfold gibbs
  rw [← Finset.sum_div]
  exact div_self (partitionZ_pos H β).ne'

/-- The modular Hamiltonian of the Gibbs reference is `β * H + log Z`
pointwise: the `thermal` clause of `CalibrationData` with `logZ = log
(partitionZ H β)`. -/
theorem neg_log_gibbs [Nonempty Ω] (H : Ω → ℝ) (β : ℝ) (x : Ω) :
    -Real.log (gibbs H β x)
      = β * H x + Real.log (partitionZ H β) := by
  unfold gibbs gibbsWeight
  rw [Real.log_div (Real.exp_pos _).ne' (partitionZ_pos H β).ne',
    Real.log_exp]
  ring

/-- **Exact affine law.**  The modular energy ledger of
`Geometry/InternalEnergyInertia.lean` against the Gibbs reference is
`β * ⟨H⟩_p + (∑ p) * log Z`: the constant-offset instance of
`modular_energy_split` with `Bc = β * H / (2π)` and the constant function
`log (partitionZ H β)`. -/
theorem modularEnergy_gibbs [Nonempty Ω] (H : Ω → ℝ) (β : ℝ) (p : Ω → ℝ) :
    modularEnergy p (gibbs H β)
      = β * expectation p H + mass p * Real.log (partitionZ H β) := by
  have hsplit : ∀ x, -Real.log (gibbs H β x)
      = 2 * Real.pi * (β * H x / (2 * Real.pi))
        + (fun _ : Ω => Real.log (partitionZ H β)) x := by
    intro x
    rw [neg_log_gibbs]
    have hπ : (2 * Real.pi) ≠ 0 := by positivity
    field_simp
  have h := modular_energy_split (gibbs H β)
    (fun x => β * H x / (2 * Real.pi)) (fun _ => Real.log (partitionZ H β))
    p hsplit
  unfold modularEnergy expectation mass
  rw [h]
  have hπ : (2 * Real.pi) ≠ 0 := by positivity
  have h1 : ∑ x, p x * (β * H x / (2 * Real.pi))
      = (β / (2 * Real.pi)) * ∑ x, p x * H x := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro x _
    ring
  rw [h1, ← Finset.sum_mul]
  field_simp

/-- The affine law for a normalized state: `β * ⟨H⟩_p + log Z`. -/
theorem modularEnergy_gibbs_normalized [Nonempty Ω] (H : Ω → ℝ) (β : ℝ)
    (p : Ω → ℝ) (hp : ∑ x, p x = 1) :
    modularEnergy p (gibbs H β)
      = β * expectation p H + Real.log (partitionZ H β) := by
  rw [modularEnergy_gibbs]
  unfold mass
  rw [hp, one_mul]

/-- The expectation of `H` in its own Gibbs state, in ledger form. -/
theorem modularEnergy_gibbs_self [Nonempty Ω] (H : Ω → ℝ) (β : ℝ) :
    modularEnergy (gibbs H β) (gibbs H β)
      = β * expectation (gibbs H β) H + Real.log (partitionZ H β) :=
  modularEnergy_gibbs_normalized H β (gibbs H β) (gibbs_sum H β)

/-! ## (2) Clausius form of the cap first law -/

/-- **Clausius form.**  For a normalized state and a Gibbs reference the cap
first law reads `S(p) - S(τ) = β (⟨H⟩_p - ⟨H⟩_τ) - D(p ‖ τ)`; the offset
`log Z` cancels. -/
theorem cap_firstLaw_gibbs [Nonempty Ω] (H : Ω → ℝ) (β : ℝ)
    (p : Ω → ℝ) (hp : ∑ x, p x = 1) :
    shannon p - shannon (gibbs H β)
      = β * (expectation p H - expectation (gibbs H β) H)
        - kl p (gibbs H β) := by
  have h := cap_firstLaw_modular p (gibbs H β) (gibbs_pos H β)
  rw [modularEnergy_gibbs_normalized H β p hp, modularEnergy_gibbs_self] at h
  linarith

/-- Between two normalized states the ledger difference is `β` times the
energy expectation difference. -/
theorem modularEnergy_gibbs_sub [Nonempty Ω] (H : Ω → ℝ) (β : ℝ) (p q : Ω → ℝ)
    (hp : ∑ x, p x = 1) (hq : ∑ x, q x = 1) :
    modularEnergy p (gibbs H β) - modularEnergy q (gibbs H β)
      = β * (expectation p H - expectation q H) := by
  rw [modularEnergy_gibbs_normalized H β p hp,
    modularEnergy_gibbs_normalized H β q hq]
  ring

/-! ## (3) Uniqueness of the slope and of the reference -/

/-- **Slope uniqueness.**  Two Gibbs references for one `H` whose ledgers
agree on two normalized states with different energy expectations share
the inverse temperature.  Ledger-side form of `gibbs_beta_injective`. -/
theorem slope_unique [Nonempty Ω] (H : Ω → ℝ) (β β' : ℝ) (p q : Ω → ℝ)
    (hp : ∑ x, p x = 1) (hq : ∑ x, q x = 1)
    (hne : expectation p H ≠ expectation q H)
    (h1 : modularEnergy p (gibbs H β) = modularEnergy p (gibbs H β'))
    (h2 : modularEnergy q (gibbs H β) = modularEnergy q (gibbs H β')) :
    β = β' := by
  have d1 := modularEnergy_gibbs_sub H β p q hp hq
  have d2 := modularEnergy_gibbs_sub H β' p q hp hq
  have h3 : β * (expectation p H - expectation q H)
      = β' * (expectation p H - expectation q H) := by
    rw [← d1, ← d2, h1, h2]
  exact mul_right_cancel₀ (sub_ne_zero.mpr hne) h3

/-- The ledger per unit energy expectation between two normalized states
with different expectations is `β`. -/
theorem slope_eq_beta [Nonempty Ω] (H : Ω → ℝ) (β : ℝ) (p q : Ω → ℝ)
    (hp : ∑ x, p x = 1) (hq : ∑ x, q x = 1)
    (hne : expectation p H ≠ expectation q H) :
    (modularEnergy p (gibbs H β) - modularEnergy q (gibbs H β))
      / (expectation p H - expectation q H) = β := by
  rw [modularEnergy_gibbs_sub H β p q hp hq]
  exact mul_div_cancel_right₀ β (sub_ne_zero.mpr hne)

/-- The ledger of a point mass is the modular Hamiltonian at that point. -/
theorem modularEnergy_single (tau : Ω → ℝ) (i : Ω) :
    modularEnergy (Pi.single i (1 : ℝ)) tau = -Real.log (tau i) := by
  unfold modularEnergy
  rw [Finset.sum_eq_single i]
  · simp
  · intro b _ hb
    simp [hb]
  · intro h
    exact absurd (Finset.mem_univ i) h

/-- The energy expectation of a point mass is the energy at that point. -/
theorem expectation_single (H : Ω → ℝ) (i : Ω) :
    expectation (Pi.single i (1 : ℝ)) H = H i := by
  unfold expectation
  rw [Finset.sum_eq_single i]
  · simp
  · intro b _ hb
    simp [hb]
  · intro h
    exact absurd (Finset.mem_univ i) h

/-- A point mass is normalized. -/
theorem sum_single (i : Ω) :
    ∑ x, (Pi.single i (1 : ℝ) : Ω → ℝ) x = 1 := by
  rw [Finset.sum_eq_single i]
  · simp
  · intro b _ hb
    simp [hb]
  · intro h
    exact absurd (Finset.mem_univ i) h

/-- **Reference uniqueness.**  A faithful reference is determined by its
ledger on all states: two faithful references with equal ledgers on every
state are equal. -/
theorem reference_determined_by_ledger (tau tau' : Ω → ℝ)
    (hτ : ∀ x, 0 < tau x) (hτ' : ∀ x, 0 < tau' x)
    (h : ∀ p : Ω → ℝ, modularEnergy p tau = modularEnergy p tau') :
    tau = tau' := by
  funext x
  have hx := h (Pi.single x 1)
  rw [modularEnergy_single, modularEnergy_single, neg_inj] at hx
  rw [← Real.exp_log (hτ x), ← Real.exp_log (hτ' x), hx]

omit [DecidableEq Ω] in
/-- **Additive freedom.**  Shifting the energy by a constant leaves the
corpus Gibbs reference unchanged. -/
theorem gibbs_add_const [Nonempty Ω] (H : Ω → ℝ) (β c : ℝ) :
    gibbs (fun x => H x + c) β = gibbs H β := by
  funext x
  simp only [gibbs, gibbsWeight, partitionZ]
  have hZ : ∑ y, Real.exp (-β * (H y + c))
      = Real.exp (-β * c) * ∑ y, Real.exp (-β * H y) := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro y _
    rw [← Real.exp_add]
    congr 1
    ring
  have hx : Real.exp (-β * (H x + c))
      = Real.exp (-β * c) * Real.exp (-β * H x) := by
    rw [← Real.exp_add]
    congr 1
    ring
  rw [hZ, hx, mul_div_mul_left _ _ (Real.exp_pos _).ne']

/-- **Converse of the additive freedom.**  Two energies with one Gibbs
reference at one nonzero `β` differ by a constant: `c = (log Z - log Z') /
β`, from `neg_log_gibbs` at each point. -/
theorem gibbs_energy_unique_up_to_const [Nonempty Ω] (H H' : Ω → ℝ) (β : ℝ)
    (hβ : β ≠ 0) (h : gibbs H β = gibbs H' β) :
    ∃ c : ℝ, ∀ x, H' x = H x + c := by
  refine ⟨(Real.log (partitionZ H β) - Real.log (partitionZ H' β)) / β, ?_⟩
  intro x
  have hx : -Real.log (gibbs H β x) = -Real.log (gibbs H' β x) := by rw [h]
  rw [neg_log_gibbs, neg_log_gibbs] at hx
  field_simp
  linarith

omit [DecidableEq Ω] in
/-- **Every faithful normalized reference is Gibbs for its own modular
Hamiltonian at `β = 1`.**  The clause `(H, β)` carries its content in `H`. -/
theorem gibbs_of_neg_log (tau : Ω → ℝ) (hτ : ∀ x, 0 < tau x)
    (hsum : ∑ x, tau x = 1) :
    gibbs (fun x => -Real.log (tau x)) 1 = tau := by
  have hZ : partitionZ (fun x => -Real.log (tau x)) 1 = 1 := by
    unfold partitionZ gibbsWeight
    simp only [neg_mul, one_mul, neg_neg]
    rw [← hsum]
    exact Finset.sum_congr rfl (fun x _ => Real.exp_log (hτ x))
  funext x
  simp only [gibbs, gibbsWeight, hZ, div_one, neg_mul, one_mul, neg_neg]
  exact Real.exp_log (hτ x)

/-- The declared one-parameter identification: one declared positive
inverse temperature `beta` per ledger energy unit.  The number is declared,
never derived; the energy `H` it is paired with is supplied at use. -/
structure TemperatureDeclaration where
  /-- The declared inverse temperature per ledger energy unit. -/
  beta : ℝ
  /-- The declared multiplier is positive, as `CalibrationData` requires. -/
  beta_pos : 0 < beta

/-- The reference state of a temperature declaration for the energy `H`. -/
def TemperatureDeclaration.reference (H : Ω → ℝ)
    (d : TemperatureDeclaration) : Ω → ℝ :=
  gibbs H d.beta

/-- Under a temperature declaration the ledger per unit energy expectation
is the declared `beta`. -/
theorem TemperatureDeclaration.ledger_slope [Nonempty Ω] (H : Ω → ℝ)
    (d : TemperatureDeclaration) (p q : Ω → ℝ)
    (hp : ∑ x, p x = 1) (hq : ∑ x, q x = 1)
    (hne : expectation p H ≠ expectation q H) :
    (modularEnergy p (d.reference H) - modularEnergy q (d.reference H))
      / (expectation p H - expectation q H) = d.beta :=
  slope_eq_beta H d.beta p q hp hq hne

/-- The declared repair-law datum whose reference is the Gibbs state of `H`
at the declared `beta`, with the trivial visible datum.  This is the
`RepairLawData` object against which the calibration row is stated. -/
def gibbsRepairLaw [Nonempty Ω] (H : Ω → ℝ) (d : TemperatureDeclaration) :
    RepairLawData Ω Unit where
  ref := gibbs H d.beta
  ref_pos := gibbs_pos H d.beta
  ref_law := gibbs_sum H d.beta
  visible := fun _ => ()

/-- **Inhabitant of the clock-and-energy calibration register row.**  A
temperature declaration for the energy `H` gives a `CalibrationData` datum
for `gibbsRepairLaw H d` with `energy = H`, the declared `beta`, and `logZ
= log (partitionZ H beta)`; the `thermal` clause is `neg_log_gibbs`.  The
row stays declared; this supplies one of its inhabitants. -/
def TemperatureDeclaration.toCalibrationData [Nonempty Ω] (H : Ω → ℝ)
    (d : TemperatureDeclaration) : CalibrationData (gibbsRepairLaw H d) where
  energy := H
  beta := d.beta
  logZ := Real.log (partitionZ H d.beta)
  beta_pos := d.beta_pos
  thermal := neg_log_gibbs H d.beta

/-- The reference of the induced calibration datum is the declared
reference, in the form `ref_eq_gibbs` returns it. -/
theorem toCalibrationData_ref [Nonempty Ω] (H : Ω → ℝ)
    (d : TemperatureDeclaration) (x : Ω) :
    (gibbsRepairLaw H d).ref x
      = gibbs (d.toCalibrationData H).energy (d.toCalibrationData H).beta x :=
  FourLawSurface.ref_eq_gibbs (gibbsRepairLaw H d) (d.toCalibrationData H) x

/-! ## (4) The screen field instance -/

/-- Declared potential history for the `e`-th seam: zero at step `0`, and
the `e`-th seam basis vector scaled by `e + 1` at every later step.  The
scaling is a declared choice that makes the family energy non-constant. -/
def seedPotential (e : Fin 30) : ℕ → Fin 30 → ℝ :=
  fun n => if n = 0 then 0 else ((e.val : ℝ) + 1) • (Pi.single e 1 : Fin 30 → ℝ)

theorem seedPotential_zero (e : Fin 30) : seedPotential e 0 = 0 := by
  simp [seedPotential]

theorem seedPotential_succ (e : Fin 30) (n : ℕ) :
    seedPotential e (n + 1) = ((e.val : ℝ) + 1) • (Pi.single e 1 : Fin 30 → ℝ) := by
  simp [seedPotential]

/-- The screen field energy of the `e`-th declared configuration at the
declared step `h`, read at window `0` with zero port potential: the
function `fieldEnergyScaled` of `Screen/ScaledMaxwellStability.lean`. -/
def screenFieldEnergy (h : ℝ) (e : Fin 30) : ℝ :=
  fieldEnergyScaled h (seedPotential e) (fun _ => 0) 0

/-- Closed form of the instance energy: `(1 / 2) * h⁻¹ ^ 2 * (e + 1) ^ 2`
(electric part only; the magnetic part vanishes at window `0`). -/
theorem screenFieldEnergy_eq (h : ℝ) (e : Fin 30) :
    screenFieldEnergy h e = (1 / 2) * h⁻¹ ^ 2 * ((e.val : ℝ) + 1) ^ 2 := by
  have hE : electricFieldScaled h (seedPotential e) (fun _ => 0) 0
      = -(h⁻¹ • (((e.val : ℝ) + 1) • (Pi.single e 1 : Fin 30 → ℝ))) := by
    unfold electricFieldScaled
    rw [seedPotential_succ, seedPotential_zero, sub_zero]
    simp
  have hB : magneticField (seedPotential e) 0 = 0 := by
    unfold magneticField
    rw [seedPotential_zero, map_zero]
  unfold screenFieldEnergy fieldEnergyScaled
  rw [hE, hB]
  have hF : faceInner 0 (magneticField (seedPotential e) 1) = 0 := by
    simp [faceInner]
  rw [hF, mul_zero, add_zero]
  unfold realSeamEnergy
  rw [Finset.sum_eq_single e]
  · simp
    ring
  · intro b _ hb
    simp [hb]
  · intro h
    exact absurd (Finset.mem_univ e) h

/-- The screen Gibbs reference at declared step `h` and declared inverse
temperature `β`: the corpus `gibbs` applied to the screen field energy. -/
def screenGibbs (h β : ℝ) : Fin 30 → ℝ := gibbs (screenFieldEnergy h) β

/-- **Affine law on the screen instance.**  The modular energy ledger against
the screen Gibbs reference is `β` times the expectation of the screen field
energy plus the mass times `log Z`: the same function `fieldEnergyScaled`
appears on the ledger side and on the field side. -/
theorem screen_affine_law (h β : ℝ) (p : Fin 30 → ℝ) :
    modularEnergy p (screenGibbs h β)
      = β * expectation p (screenFieldEnergy h)
        + mass p * Real.log (partitionZ (screenFieldEnergy h) β) :=
  modularEnergy_gibbs (screenFieldEnergy h) β p

/-- Clausius form on the screen instance. -/
theorem screen_cap_firstLaw (h β : ℝ) (p : Fin 30 → ℝ) (hp : ∑ x, p x = 1) :
    shannon p - shannon (screenGibbs h β)
      = β * (expectation p (screenFieldEnergy h)
          - expectation (screenGibbs h β) (screenFieldEnergy h))
        - kl p (screenGibbs h β) :=
  cap_firstLaw_gibbs (screenFieldEnergy h) β p hp

/-- The calibration-row inhabitant of the screen instance: `energy =
screenFieldEnergy h` at the declared `beta`. -/
def screenCalibration (h : ℝ) (d : TemperatureDeclaration) :
    CalibrationData (gibbsRepairLaw (screenFieldEnergy h) d) :=
  d.toCalibrationData (screenFieldEnergy h)

theorem screenCalibration_energy (h : ℝ) (d : TemperatureDeclaration) :
    (screenCalibration h d).energy = screenFieldEnergy h := rfl

/-- The instance energy separates the first two configurations at any
nonzero step. -/
theorem screenFieldEnergy_zero_ne_one (h : ℝ) (hh : h ≠ 0) :
    screenFieldEnergy h 0 ≠ screenFieldEnergy h 1 := by
  rw [screenFieldEnergy_eq, screenFieldEnergy_eq]
  have hne : h⁻¹ ^ 2 ≠ 0 := pow_ne_zero 2 (inv_ne_zero hh)
  have h0 : ((0 : Fin 30).val : ℝ) = 0 := by simp
  have h1 : ((1 : Fin 30).val : ℝ) = 1 := by simp
  rw [h0, h1]
  intro heq
  apply hne
  linarith

/-- **`β` is determined by the reference.**  Once the family and the nonzero
step are declared, two inverse temperatures giving one screen Gibbs
reference are equal: `gibbs_beta_injective` of
`Thermodynamics/FiniteConditionalRepair.lean` at the levels `0` and `1`. -/
theorem screen_beta_unique (h β β' : ℝ) (hh : h ≠ 0)
    (hβ : screenGibbs h β = screenGibbs h β') : β = β' :=
  gibbs_beta_injective (screenFieldEnergy h) β β' 0 1
    (screenFieldEnergy_zero_ne_one h hh) hβ

/-! ## (5) Boundary: a constant-energy, non-Gibbs reference -/

/-- Declared two-point boundary reference `![1, 1 / 2]` (not normalized). -/
def boundaryRef : Fin 2 → ℝ := ![1, 1 / 2]

/-- Declared constant energy on the two-point boundary example. -/
def boundaryEnergy : Fin 2 → ℝ := fun _ => 0

theorem boundaryRef_pos (x : Fin 2) : 0 < boundaryRef x := by
  fin_cases x <;> norm_num [boundaryRef]

/-- The Gibbs reference of a constant energy is uniform. -/
theorem gibbs_boundaryEnergy (β : ℝ) (x : Fin 2) :
    gibbs boundaryEnergy β x = 1 / 2 := by
  simp [gibbs, gibbsWeight, partitionZ, boundaryEnergy]

/-- The boundary reference is not of Gibbs form for the constant energy at
any inverse temperature (the exclusion goes through non-uniformity: a
constant energy has the uniform Gibbs state). -/
theorem boundary_not_gibbs (β : ℝ) : boundaryRef ≠ gibbs boundaryEnergy β := by
  intro h
  have h0 := congrFun h 0
  rw [gibbs_boundaryEnergy] at h0
  norm_num [boundaryRef] at h0

/-- The two point masses have equal energy expectation. -/
theorem boundary_equal_expectation :
    expectation (Pi.single 0 (1 : ℝ)) boundaryEnergy
      = expectation (Pi.single 1 (1 : ℝ)) boundaryEnergy := by
  rw [expectation_single, expectation_single]
  rfl

/-- The two point masses have different ledgers against the boundary
reference. -/
theorem boundary_ledger_ne :
    modularEnergy (Pi.single 0 (1 : ℝ)) boundaryRef
      ≠ modularEnergy (Pi.single 1 (1 : ℝ)) boundaryRef := by
  rw [modularEnergy_single, modularEnergy_single]
  simp only [boundaryRef, Matrix.cons_val_zero, Matrix.cons_val_one]
  rw [Real.log_one, neg_zero]
  intro hlog
  apply Real.log_ne_zero_of_pos_of_ne_one (by norm_num : (0 : ℝ) < 1 / 2)
    (by norm_num)
  linarith

/-- **Boundary.**  Against the boundary reference the ledger is not an affine
function of the energy expectation on normalized states.  For the constant
energy this is the statement that the ledger is non-constant. -/
theorem boundary_not_affine :
    ¬ ∃ a b : ℝ, ∀ p : Fin 2 → ℝ, ∑ x, p x = 1 →
      modularEnergy p boundaryRef = a * expectation p boundaryEnergy + b := by
  rintro ⟨a, b, hab⟩
  apply boundary_ledger_ne
  rw [hab _ (sum_single 0), hab _ (sum_single 1), boundary_equal_expectation]

/-! ## (6) Temperature dictionary through the declared tick -/

/-- Laboratory inverse temperature in inverse joules for the declared tick:
`β / hbarOverTau`.  A conditional conversion of the ledger `β`. -/
def labInverseTemperature (a : SIAnchors) (cal : ClockCalibration) (β : ℝ) : ℝ :=
  β / hbarOverTau a cal

/-- Laboratory thermal energy `k_B T` in joules for the declared tick:
`hbarOverTau / β`.  A conditional conversion of the ledger `β`. -/
def thermalEnergyJoules (a : SIAnchors) (cal : ClockCalibration) (β : ℝ) : ℝ :=
  hbarOverTau a cal / β

/-- **Exact conversion.**  The dimensionless product `β * E` is invariant
under the conversion of both factors to laboratory units. -/
theorem beta_energy_invariant (a : SIAnchors) (cal : ClockCalibration)
    (β E : ℝ) :
    β * E = labInverseTemperature a cal β * energyJoules a cal E := by
  unfold labInverseTemperature energyJoules
  have hk := (hbarOverTau_pos a cal).ne'
  field_simp

/-- The laboratory inverse temperature and thermal energy are reciprocal. -/
theorem labInverseTemperature_mul_thermalEnergy (a : SIAnchors)
    (cal : ClockCalibration) {β : ℝ} (hβ : β ≠ 0) :
    labInverseTemperature a cal β * thermalEnergyJoules a cal β = 1 := by
  unfold labInverseTemperature thermalEnergyJoules
  have hk := (hbarOverTau_pos a cal).ne'
  field_simp

/-- **Non-forcing.**  Two distinct declared ticks give two distinct
laboratory readings of one nonzero ledger `β`. -/
theorem labInverseTemperature_not_forced (a : SIAnchors)
    (cal cal' : ClockCalibration) (hne : cal.tau ≠ cal'.tau) {β : ℝ}
    (hβ : β ≠ 0) :
    labInverseTemperature a cal β ≠ labInverseTemperature a cal' β := by
  intro heq
  unfold labInverseTemperature at heq
  have h1 := (hbarOverTau_pos a cal).ne'
  have h2 := (hbarOverTau_pos a cal').ne'
  rw [div_eq_div_iff h1 h2] at heq
  have key : hbarOverTau a cal = hbarOverTau a cal' :=
    (mul_left_cancel₀ hβ heq).symm
  unfold hbarOverTau at key
  rw [div_eq_div_iff cal.tau_pos.ne' cal'.tau_pos.ne'] at key
  have hc : (a.planckConstant : ℝ) / (2 * Real.pi) ≠ 0 :=
    (div_pos (planckConstant_real_pos a) (mul_pos two_pos Real.pi_pos)).ne'
  exact hne (mul_left_cancel₀ hc key).symm

end

end OPH.GibbsReferenceEnergyIdentification

#print axioms OPH.GibbsReferenceEnergyIdentification.gibbs_pos
#print axioms OPH.GibbsReferenceEnergyIdentification.gibbs_sum
#print axioms OPH.GibbsReferenceEnergyIdentification.neg_log_gibbs
#print axioms OPH.GibbsReferenceEnergyIdentification.modularEnergy_gibbs
#print axioms OPH.GibbsReferenceEnergyIdentification.cap_firstLaw_gibbs
#print axioms OPH.GibbsReferenceEnergyIdentification.modularEnergy_gibbs_sub
#print axioms OPH.GibbsReferenceEnergyIdentification.slope_unique
#print axioms OPH.GibbsReferenceEnergyIdentification.slope_eq_beta
#print axioms OPH.GibbsReferenceEnergyIdentification.reference_determined_by_ledger
#print axioms OPH.GibbsReferenceEnergyIdentification.gibbs_add_const
#print axioms OPH.GibbsReferenceEnergyIdentification.gibbs_energy_unique_up_to_const
#print axioms OPH.GibbsReferenceEnergyIdentification.gibbs_of_neg_log
#print axioms OPH.GibbsReferenceEnergyIdentification.TemperatureDeclaration.toCalibrationData
#print axioms OPH.GibbsReferenceEnergyIdentification.toCalibrationData_ref
#print axioms OPH.GibbsReferenceEnergyIdentification.screenFieldEnergy_eq
#print axioms OPH.GibbsReferenceEnergyIdentification.screen_affine_law
#print axioms OPH.GibbsReferenceEnergyIdentification.screen_cap_firstLaw
#print axioms OPH.GibbsReferenceEnergyIdentification.screenCalibration
#print axioms OPH.GibbsReferenceEnergyIdentification.screen_beta_unique
#print axioms OPH.GibbsReferenceEnergyIdentification.boundary_not_gibbs
#print axioms OPH.GibbsReferenceEnergyIdentification.boundary_not_affine
#print axioms OPH.GibbsReferenceEnergyIdentification.beta_energy_invariant
#print axioms OPH.GibbsReferenceEnergyIdentification.labInverseTemperature_mul_thermalEnergy
#print axioms OPH.GibbsReferenceEnergyIdentification.labInverseTemperature_not_forced
