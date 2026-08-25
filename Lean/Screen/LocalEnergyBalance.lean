import ScaledMaxwellStability

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.LocalEnergyBalance

open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability

/-!
# Local energy balance on the committed carrier: a discrete Poynting law

WHAT IS PROVED.  Exact finite real linear algebra on the committed
twelve-port, thirty-seam, twenty-oriented-face complex, reusing the
committed oriented face curvature `faceCurvature` (C), its transpose
`faceCodifferential` (Cᵀ), the committed scaled staggered packet
`electricFieldScaled`, `magneticField`, `AmpereEvolutionScaled`, and the
committed scaled staggered form `fieldEnergyScaled` of
`ScaledMaxwellStability`.

(A) Face-local densities.  Every committed seam bounds exactly two faces
with opposite induced signs, so `∑ f, C f e ^ 2 = 2` on every seam
(`faceIncidence_sq_sum`).  This is the fact that makes a FACE-local
density exact: the electric seam energy `(1/2) E e ^ 2` splits into equal
quarters over the two faces of the seam, and the staggered magnetic
energy `(1/2) B n f * B (n+1) f` is face-native.  The face energy density
`faceEnergyDensity h A φ n f =
  (1/4) ∑ e, C f e ^ 2 * (E n e) ^ 2 + (1/2) B n f * B (n+1) f`
sums over the twenty faces to the committed `fieldEnergyScaled h A φ n`
exactly (`faceEnergyDensity_sum`).  The equal split is a declared
symmetric convention within a one-parameter-per-seam family of exact
face-local laws: any split `(λ_e, 1 - λ_e)` also sums to the committed
energy and yields a conserved flux with a different local law, and the
module fixes the symmetric member, which is what gives the mean-`B`
Poynting form below.  A port-local density is a live
alternative that this module does not construct: the magnetic face energy
has no port-native home, and the transport term `Cᵀ B` pairs seams with
faces, so the port split would be a further declared choice.

(B) Seam flux.  With the half-step mean `Ē n = (1/2) (E n + E (n+1))`
and the integer-step field `B (n+1)` between them, the flux rate out of
face `f` through its seam `e` is
`seamFluxRate h A φ n f e =
  C f e * Ē n e * B (n+1) f - (1/2) C f e ^ 2 * Ē n e * (Cᵀ B (n+1)) e`.
On the committed carrier this is the Poynting-shaped expression
`C f e * Ē n e * (1/2) (B (n+1) f + B (n+1) f')` with `f'` the opposite
face of the seam (`seamFluxRate_poynting_form`): orientation sign, times
the mean electric field on the seam, times the mean magnetic field across
the seam.  The flux rates of the two faces of a seam cancel
(`seamFluxRate_conserved`), which is the discrete statement that energy
leaving one face through a seam enters the opposite face.

(C) HEADLINE, the local balance (`local_energy_balance`).  Under the
committed scaled Ampere evolution with step `h ≠ 0`, at every face `f`
and every step `n`,
`u (n+1) f - u n f = -h * (div S n f + w n f)`,
with `div S n f = ∑ e, seamFluxRate h A φ n f e` the discrete divergence
of the flux at the face and `w n f = (1/2) ∑ e, C f e ^ 2 * Ē n e * J n e`
the face's share of the seam work `Ē · J`.  Signs: the time difference
is forward across the step `n → n+1`; positive divergence is net outflow;
positive work against a current lowers the field energy.  Summing over
the twenty faces, the divergence sums to zero (`faceFluxDivergence_sum`),
the work sums to `⟨Ē n, J n⟩` (`faceWorkRate_sum`), and the committed
global identity `energy_balance_scaled` is recovered from the local one
(`global_balance_from_local`): the local law refines the committed global
law exactly.

(D) Continuity packaging (`fourDivergence_eq_work`,
`fourDivergence_source_free`).  The pair (density, flux) has discrete
four-divergence `u (n+1) f - u n f + h * div S n f` equal to `-h * w n f`,
hence zero in the source-free case.  This is one conserved current on the
face lattice, the finite precursor of one component of a weak Ward law.

STANDING NEGATIVES, cited at scope.  The step index `n` is a declared
evolution parameter.  The standing result `no_rate_functional` of
`RateNonidentifiability` concerns abstract transition systems; its bridge
to the repair layer fails (`acceptedStep_changes_obs`,
`firstLock_obs_determined` of `RateBridgeObstruction`), so that result
does not forbid a repair-layer clock and is not used here in either
direction.  The kinetic term behind `AmpereEvolutionScaled` is a declared
enrichment: `realizedHistory_legendre_nonidentifiability_receipt` shows
realized source histories do not select a velocity curvature or Legendre
map, so the Lagrangian shape consumed through the scaled packet is
declared, and every statement below is conditional on that declaration.

THE THREE OPEN ROWS.  Source clock and duration: the step `h` and index
`n` are declared; no duration is derived here and that row is untouched.
Physical spacetime attachment: the face lattice is the committed
combinatorial carrier; no metric, no light cone, no continuum limit is
attached, and that row is touched only in that the discrete balance names
what a face-balance limit would have to preserve: the two-face
antisymmetry of the seam flux, the face-sum identity with the committed
global form, the forward-difference and outflow sign conventions, and the
absence of area or Hodge weights; nothing is discharged.
Coupled action: the sources `ρ` and `J` are declared inputs and the
scaled kinetic term is declared; that row is untouched.

BOUNDARY TO THE EINSTEIN COMPOSITION.  The Einstein branch consumes a
same-source symmetric stress tensor, a weak Ward identity for all its
components, locality and covariance, and a discrete-to-continuum face
balance.  This module supplies one component (energy) of one conserved
current on the committed face lattice, for the committed screen-field
energy.  It does not supply momentum or stress components, a symmetric
tensor, covariance, or a continuum limit.  The identification of
`fieldEnergyScaled` with the modular energy of `CapFirstLaw` consumed by
the Einstein branch is a named open join and is not asserted.  No premise
is discharged; no frozen prediction is added.

Axiom audit.  Every proof composes the committed receipts with exact real
linear algebra and two kernel `decide` checks on the committed integer
incidence table (`faceIncidenceZ_sq_sum`, `faceIncidenceZ_opposite`); the module adds no project axiom and uses no native
decision procedure.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`.
-/

noncomputable section

/-! ## (A) The two-faces-per-seam fact in squared form -/

/-- Squared incidence sums to two on every seam: the committed
two-faces-per-seam fact in the form the density split uses. -/
theorem faceIncidenceZ_sq_sum :
    ∀ e : Fin 30, (∑ f : Fin 20, faceIncidenceZ f e ^ 2) = 2 := by
  decide

/-- Real form of `faceIncidenceZ_sq_sum`. -/
theorem faceIncidence_sq_sum (e : Fin 30) :
    (∑ f : Fin 20, faceIncidenceR f e ^ 2) = 2 := by
  have h := faceIncidenceZ_sq_sum e
  unfold faceIncidenceR
  exact_mod_cast h

/-- On a seam, two distinct faces with nonzero incidence carry opposite
signs: the committed opposite-orientation fact in pointwise form. -/
theorem faceIncidenceZ_opposite :
    ∀ (e : Fin 30) (f g : Fin 20), f ≠ g → faceIncidenceZ f e ≠ 0 →
      faceIncidenceZ g e ≠ 0 → faceIncidenceZ g e = -faceIncidenceZ f e := by
  set_option maxRecDepth 16384 in
    decide

/-- A face with nonzero incidence on a seam has an opposite face: the
codifferential at that seam is the two-term sum. -/
theorem faceCodifferential_two_terms (e : Fin 30) (f g : Fin 20)
    (hfg : f ≠ g) (hf : faceIncidenceZ f e ≠ 0) (hg : faceIncidenceZ g e ≠ 0)
    (F : Fin 20 → ℝ) :
    faceCodifferential F e = faceIncidenceR f e * F f + faceIncidenceR g e * F g := by
  rw [faceCodifferential_apply]
  refine Finset.sum_eq_add_of_mem f g (Finset.mem_univ f) (Finset.mem_univ g) hfg ?_
  intro k _ hk
  have hzero : faceIncidenceZ k e = 0 := by
    by_contra hne
    have h1 := faceIncidenceZ_opposite e f k hk.1.symm hf hne
    have h2 := faceIncidenceZ_opposite e g k hk.2.symm hg hne
    have h3 := faceIncidenceZ_opposite e f g hfg hf hg
    omega
  unfold faceIncidenceR
  rw [hzero]
  simp

/-! ## (A) Face-local energy density -/

/-- Half-step mean of the scaled electric field across the step `n → n+1`:
`Ē n = (1/2) (E n + E (n+1))`. -/
def meanElectric (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) : Fin 30 → ℝ :=
  (1 / 2 : ℝ) • (electricFieldScaled h A φ n + electricFieldScaled h A φ (n + 1))

theorem meanElectric_apply (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) (e : Fin 30) :
    meanElectric h A φ n e =
      (1 / 2) * (electricFieldScaled h A φ n e + electricFieldScaled h A φ (n + 1) e) := by
  unfold meanElectric
  simp only [Pi.smul_apply, Pi.add_apply, smul_eq_mul]

/-- Face-local energy density at step `n`: one quarter of the squared
electric field on each of the face's three seams (each seam is shared by
two faces), plus the face-native staggered magnetic term. -/
def faceEnergyDensity (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) (f : Fin 20) : ℝ :=
  (1 / 4) * (∑ e : Fin 30, faceIncidenceR f e ^ 2 * electricFieldScaled h A φ n e ^ 2)
    + (1 / 2) * magneticField A n f * magneticField A (n + 1) f

/-- The face densities sum to the committed scaled staggered form. -/
theorem faceEnergyDensity_sum (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) :
    (∑ f : Fin 20, faceEnergyDensity h A φ n f) = fieldEnergyScaled h A φ n := by
  unfold faceEnergyDensity fieldEnergyScaled realSeamEnergy faceInner
  rw [Finset.sum_add_distrib, ← Finset.mul_sum, Finset.sum_comm]
  have helec : (∑ e : Fin 30, ∑ f : Fin 20,
      faceIncidenceR f e ^ 2 * electricFieldScaled h A φ n e ^ 2) =
      2 * ∑ e : Fin 30, electricFieldScaled h A φ n e ^ 2 := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun e _ ↦ ?_
    rw [← Finset.sum_mul, faceIncidence_sq_sum]
  rw [helec]
  have hmag : (∑ f : Fin 20, (1 / 2) * magneticField A n f * magneticField A (n + 1) f) =
      (1 / 2) * ∑ f : Fin 20, magneticField A n f * magneticField A (n + 1) f := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun f _ ↦ by ring
  rw [hmag]
  ring

/-! ## (B) Seam flux rate and its two-face cancellation -/

/-- Flux rate of energy out of face `f` through its seam `e`, built from
the half-step mean electric field and the integer-step magnetic field
through the committed incidence operators. -/
def seamFluxRate (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) (f : Fin 20) (e : Fin 30) : ℝ :=
  faceIncidenceR f e * meanElectric h A φ n e * magneticField A (n + 1) f
    - (1 / 2) * faceIncidenceR f e ^ 2 * meanElectric h A φ n e
        * faceCodifferential (magneticField A (n + 1)) e

/-- **Two-face cancellation.**  On every seam the flux rates of all faces
sum to zero: energy leaving one face through the seam enters the opposite
face. -/
theorem seamFluxRate_conserved (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) (e : Fin 30) :
    (∑ f : Fin 20, seamFluxRate h A φ n f e) = 0 := by
  unfold seamFluxRate
  rw [Finset.sum_sub_distrib]
  have h1 : (∑ f : Fin 20, faceIncidenceR f e * meanElectric h A φ n e
      * magneticField A (n + 1) f) =
      meanElectric h A φ n e * faceCodifferential (magneticField A (n + 1)) e := by
    rw [faceCodifferential_apply, Finset.mul_sum]
    exact Finset.sum_congr rfl fun f _ ↦ by ring
  have h2 : (∑ f : Fin 20, (1 / 2) * faceIncidenceR f e ^ 2 * meanElectric h A φ n e
      * faceCodifferential (magneticField A (n + 1)) e) =
      (1 / 2) * meanElectric h A φ n e * faceCodifferential (magneticField A (n + 1)) e
        * ∑ f : Fin 20, faceIncidenceR f e ^ 2 := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun f _ ↦ by ring
  rw [h1, h2, faceIncidence_sq_sum]
  ring

/-- **Poynting form on the committed carrier.**  For a face `f` and the
opposite face `g` of its seam `e`, the flux rate is the orientation sign
times the mean electric field on the seam times the mean magnetic field
across the seam. -/
theorem seamFluxRate_poynting_form (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (n : ℕ) (e : Fin 30) (f g : Fin 20)
    (hfg : f ≠ g) (hf : faceIncidenceZ f e ≠ 0) (hg : faceIncidenceZ g e ≠ 0) :
    seamFluxRate h A φ n f e =
      faceIncidenceR f e * meanElectric h A φ n e
        * ((1 / 2) * (magneticField A (n + 1) f + magneticField A (n + 1) g)) := by
  unfold seamFluxRate
  rw [faceCodifferential_two_terms e f g hfg hf hg]
  have hopp : faceIncidenceR g e = -faceIncidenceR f e := by
    unfold faceIncidenceR
    exact_mod_cast faceIncidenceZ_opposite e f g hfg hf hg
  have hsq : faceIncidenceR f e ^ 2 = 1 := by
    unfold faceIncidenceR
    have hcases : faceIncidenceZ f e = 1 ∨ faceIncidenceZ f e = -1 := by
      unfold faceIncidenceZ at hf ⊢
      split_ifs at hf ⊢ <;> simp_all
    rcases hcases with hc | hc <;> rw [hc] <;> norm_num
  rw [hopp]
  linear_combination (-(1 / 2) * faceIncidenceR f e * meanElectric h A φ n e
    * (magneticField A (n + 1) f - magneticField A (n + 1) g)) * hsq
/-! ## (C) Local balance: the headline -/

/-- Discrete divergence of the flux at face `f`: total flux rate out of the
face through its seams (positive means net outflow). -/
def faceFluxDivergence (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) (f : Fin 20) : ℝ :=
  ∑ e : Fin 30, seamFluxRate h A φ n f e

/-- Face share of the seam work rate `Ē · J`: half of each of its seams'
work, since each seam is shared by two faces. -/
def faceWorkRate (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (n : ℕ) (f : Fin 20) : ℝ :=
  (1 / 2) * ∑ e : Fin 30, faceIncidenceR f e ^ 2 * meanElectric h A φ n e * J n e

/-- Pointwise scaled Faraday law at a face, in summed incidence form. -/
theorem faraday_face (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (n : ℕ) (f : Fin 20) :
    magneticField A (n + 1) f - magneticField A n f =
      -(h * ∑ e : Fin 30, faceIncidenceR f e * electricFieldScaled h A φ n e) := by
  have hF := congrFun (faraday_law_scaled h hh A φ n) f
  simp only [Pi.sub_apply, Pi.neg_apply, Pi.smul_apply, smul_eq_mul] at hF
  rw [hF, faceCurvature_apply]

/-- **Local energy balance (discrete Poynting law).**  Under the committed
scaled Ampere evolution with step `h ≠ 0`, the change of the face energy
density across the step `n → n+1` equals minus `h` times the sum of the
flux divergence at the face and the face's share of the work `Ē · J`. -/
theorem local_energy_balance (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ J) (n : ℕ) (f : Fin 20) :
    faceEnergyDensity h A φ (n + 1) f - faceEnergyDensity h A φ n f =
      -(h * (faceFluxDivergence h A φ n f + faceWorkRate h A φ J n f)) := by
  set E := electricFieldScaled h A φ with hEdef
  set B := magneticField A with hBdef
  have hE : ∀ e, E (n + 1) e - E n e =
      h * (faceCodifferential (B (n + 1)) e - J n e) := by
    intro e
    have := congrFun (hAmp n) e
    simpa only [Pi.sub_apply, Pi.smul_apply, smul_eq_mul] using this
  have hsum : (1 / 4) * (∑ e : Fin 30, faceIncidenceR f e ^ 2 * E (n + 1) e ^ 2)
      - (1 / 4) * (∑ e : Fin 30, faceIncidenceR f e ^ 2 * E n e ^ 2)
      + h * (∑ e : Fin 30, seamFluxRate h A φ n f e)
      + h * ((1 / 2) * ∑ e : Fin 30, faceIncidenceR f e ^ 2 * meanElectric h A φ n e * J n e)
      = h * B (n + 1) f * ∑ e : Fin 30, faceIncidenceR f e * meanElectric h A φ n e := by
    simp only [Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun e _ ↦ ?_
    have he := hE e
    unfold seamFluxRate
    rw [meanElectric_apply, ← hEdef, ← hBdef]
    linear_combination ((1 / 4) * faceIncidenceR f e ^ 2 * (E n e + E (n + 1) e)) * he
  have hmean : (∑ e : Fin 30, faceIncidenceR f e * meanElectric h A φ n e) =
      (1 / 2) * ((∑ e : Fin 30, faceIncidenceR f e * E n e)
        + ∑ e : Fin 30, faceIncidenceR f e * E (n + 1) e) := by
    rw [← Finset.sum_add_distrib, Finset.mul_sum]
    refine Finset.sum_congr rfl fun e _ ↦ ?_
    rw [meanElectric_apply, ← hEdef]
    ring
  have hB0 := faraday_face h hh A φ n f
  have hB1 := faraday_face h hh A φ (n + 1) f
  rw [← hEdef, ← hBdef] at hB0 hB1
  unfold faceEnergyDensity faceFluxDivergence faceWorkRate
  rw [← hEdef, ← hBdef]
  linear_combination hsum + (h * B (n + 1) f) * hmean
    + ((1 / 2) * B (n + 1) f) * hB1 + ((1 / 2) * B (n + 1) f) * hB0

/-! ## (C) Summation: the local law refines the committed global law -/

/-- The flux divergence sums to zero over the closed committed carrier. -/
theorem faceFluxDivergence_sum (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) : (∑ f : Fin 20, faceFluxDivergence h A φ n f) = 0 := by
  unfold faceFluxDivergence
  rw [Finset.sum_comm]
  exact Finset.sum_eq_zero fun e _ ↦ seamFluxRate_conserved h A φ n e

/-- The face work rates sum to the seam work `⟨Ē n, J n⟩`. -/
theorem faceWorkRate_sum (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (n : ℕ) :
    (∑ f : Fin 20, faceWorkRate h A φ J n f) =
      realSeamInner (meanElectric h A φ n) (J n) := by
  unfold faceWorkRate realSeamInner
  rw [← Finset.mul_sum, Finset.sum_comm]
  have hinner : (∑ e : Fin 30, ∑ f : Fin 20,
      faceIncidenceR f e ^ 2 * meanElectric h A φ n e * J n e) =
      2 * ∑ e : Fin 30, meanElectric h A φ n e * J n e := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun e _ ↦ ?_
    have hrow : (∑ f : Fin 20, faceIncidenceR f e ^ 2 * meanElectric h A φ n e * J n e) =
        (∑ f : Fin 20, faceIncidenceR f e ^ 2) * (meanElectric h A φ n e * J n e) := by
      rw [Finset.sum_mul]
      exact Finset.sum_congr rfl fun f _ ↦ by ring
    rw [hrow, faceIncidence_sq_sum]
  rw [hinner]
  ring

/-- The seam work in mean form equals the committed time-averaged form. -/
theorem meanElectric_inner (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) (n : ℕ) :
    realSeamInner (meanElectric h A φ n) (J n) =
      (1 / 2) * realSeamInner
        (electricFieldScaled h A φ n + electricFieldScaled h A φ (n + 1)) (J n) := by
  unfold meanElectric realSeamInner
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  simp only [Pi.smul_apply, Pi.add_apply, smul_eq_mul]
  ring

/-- **Global balance from the local law.**  Summing `local_energy_balance`
over the twenty faces reproduces the committed `energy_balance_scaled`
statement exactly; the committed global identity is the face sum of the
local one. -/
theorem global_balance_from_local (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ J) (n : ℕ) :
    fieldEnergyScaled h A φ (n + 1) = fieldEnergyScaled h A φ n -
      (h / 2) * realSeamInner
        (electricFieldScaled h A φ n + electricFieldScaled h A φ (n + 1)) (J n) := by
  have hloc : (∑ f : Fin 20, (faceEnergyDensity h A φ (n + 1) f - faceEnergyDensity h A φ n f))
      = ∑ f : Fin 20, -(h * (faceFluxDivergence h A φ n f + faceWorkRate h A φ J n f)) :=
    Finset.sum_congr rfl fun f _ ↦ local_energy_balance h hh A φ J hAmp n f
  rw [Finset.sum_sub_distrib, faceEnergyDensity_sum, faceEnergyDensity_sum,
    Finset.sum_neg_distrib, ← Finset.mul_sum, Finset.sum_add_distrib,
    faceFluxDivergence_sum, faceWorkRate_sum, meanElectric_inner] at hloc
  linear_combination hloc

/-! ## (D) Continuity packaging: one component of a weak Ward law -/

/-- Discrete four-divergence of the pair (density, flux) at face `f` across
the step `n → n+1`: forward time difference plus `h` times the spatial
flux divergence. -/
def fourDivergence (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (n : ℕ) (f : Fin 20) : ℝ :=
  faceEnergyDensity h A φ (n + 1) f - faceEnergyDensity h A φ n f
    + h * faceFluxDivergence h A φ n f

/-- **Sourced continuity.**  The four-divergence equals minus `h` times the
face work rate: the (energy, flux) current is conserved up to the work of
the declared current. -/
theorem fourDivergence_eq_work (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ J) (n : ℕ) (f : Fin 20) :
    fourDivergence h A φ n f = -(h * faceWorkRate h A φ J n f) := by
  unfold fourDivergence
  have := local_energy_balance h hh A φ J hAmp n f
  linear_combination this

/-- **Source-free continuity.**  For zero current the four-divergence
vanishes at every face and every step: one conserved current on the face
lattice, the finite precursor of one component of a weak Ward law. -/
theorem fourDivergence_source_free (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (hAmp : AmpereEvolutionScaled h A φ (fun _ ↦ 0))
    (n : ℕ) (f : Fin 20) : fourDivergence h A φ n f = 0 := by
  rw [fourDivergence_eq_work h hh A φ (fun _ ↦ 0) hAmp n f]
  unfold faceWorkRate
  simp

/-- The face work rate vanishes on any face all of whose seams carry zero
current; in particular the four-divergence vanishes there even with
sources elsewhere (locality of the source term). -/
theorem fourDivergence_local_source_free (h : ℝ) (hh : h ≠ 0)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ J) (n : ℕ) (f : Fin 20)
    (hJ : ∀ e, faceIncidenceZ f e ≠ 0 → J n e = 0) :
    fourDivergence h A φ n f = 0 := by
  rw [fourDivergence_eq_work h hh A φ J hAmp n f]
  unfold faceWorkRate
  have hzero : (∑ e : Fin 30, faceIncidenceR f e ^ 2 * meanElectric h A φ n e * J n e) = 0 := by
    refine Finset.sum_eq_zero fun e _ ↦ ?_
    by_cases hfe : faceIncidenceZ f e = 0
    · unfold faceIncidenceR
      rw [hfe]
      simp
    · rw [hJ e hfe, mul_zero]
  rw [hzero]
  simp

/-! ## (E) Bundle receipt -/

/-- **Local energy balance receipt.**  From the committed scaled Maxwell
bundle: the face densities sum to the committed form; the local balance
holds at every face and step; the four-divergence equals minus the work;
the face-summed law reproduces the committed global balance.  Everything
is conditional on the bundle's declared step, kinetic term, and sources. -/
theorem localEnergyBalance_receipt (S : ScaledMaxwellBundle) :
    (∀ n, (∑ f : Fin 20, faceEnergyDensity S.h S.A S.phi n f) =
      fieldEnergyScaled S.h S.A S.phi n)
    ∧ (∀ n f, faceEnergyDensity S.h S.A S.phi (n + 1) f
        - faceEnergyDensity S.h S.A S.phi n f =
        -(S.h * (faceFluxDivergence S.h S.A S.phi n f
          + faceWorkRate S.h S.A S.phi S.J n f)))
    ∧ (∀ n f, fourDivergence S.h S.A S.phi n f =
        -(S.h * faceWorkRate S.h S.A S.phi S.J n f))
    ∧ (∀ n, fieldEnergyScaled S.h S.A S.phi (n + 1) = fieldEnergyScaled S.h S.A S.phi n -
        (S.h / 2) * realSeamInner
          (electricFieldScaled S.h S.A S.phi n
            + electricFieldScaled S.h S.A S.phi (n + 1)) (S.J n)) := by
  have hh : S.h ≠ 0 := ne_of_gt S.h_pos
  exact ⟨fun n ↦ faceEnergyDensity_sum S.h S.A S.phi n,
    fun n f ↦ local_energy_balance S.h hh S.A S.phi S.J S.ampere n f,
    fun n f ↦ fourDivergence_eq_work S.h hh S.A S.phi S.J S.ampere n f,
    fun n ↦ global_balance_from_local S.h hh S.A S.phi S.J S.ampere n⟩

/-- The receipt is inhabited by the committed nonstatic demonstration
bundle at `h = 1/2`. -/
theorem demo_local_balance :
    ∀ n f, faceEnergyDensity demoScaledBundle.h demoScaledBundle.A
        demoScaledBundle.phi (n + 1) f
      - faceEnergyDensity demoScaledBundle.h demoScaledBundle.A
        demoScaledBundle.phi n f =
      -(demoScaledBundle.h * (faceFluxDivergence demoScaledBundle.h
        demoScaledBundle.A demoScaledBundle.phi n f
        + faceWorkRate demoScaledBundle.h demoScaledBundle.A
          demoScaledBundle.phi demoScaledBundle.J n f)) :=
  (localEnergyBalance_receipt demoScaledBundle).2.1

end

end OPH.LocalEnergyBalance

#print axioms OPH.LocalEnergyBalance.faceIncidenceZ_sq_sum
#print axioms OPH.LocalEnergyBalance.faceIncidenceZ_opposite
#print axioms OPH.LocalEnergyBalance.faceEnergyDensity_sum
#print axioms OPH.LocalEnergyBalance.seamFluxRate_conserved
#print axioms OPH.LocalEnergyBalance.seamFluxRate_poynting_form
#print axioms OPH.LocalEnergyBalance.local_energy_balance
#print axioms OPH.LocalEnergyBalance.global_balance_from_local
#print axioms OPH.LocalEnergyBalance.fourDivergence_eq_work
#print axioms OPH.LocalEnergyBalance.fourDivergence_source_free
#print axioms OPH.LocalEnergyBalance.fourDivergence_local_source_free
#print axioms OPH.LocalEnergyBalance.localEnergyBalance_receipt
#print axioms OPH.LocalEnergyBalance.demo_local_balance
