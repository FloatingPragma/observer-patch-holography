import GibbsReferenceEnergyIdentification
import ObserverPatchHolography.CollarModularT2
import QFT.StructuralInheritance

set_option autoImplicit false

open scoped BigOperators

/-!
# Collar temperature reading: the collar modular surface pins no inverse
temperature (gravitation row issue 729; calibration import)

STATUS.  Conditional reading; negative on the pinning question.  The collar
modular surface of `ObserverPatchHolography/CollarStates.lean` (with the T2
receipts of `CollarModularT2.lean`, the layer and clause of
`CollarLayer.lean` and `CollarClause.lean`, and the bridge of
`CollarStatesBridge.lean`) carries one object with thermal content, on the realized branch of the
declared `ax:maxent` constraint of `CollarModularT2.lean`: the
realized-branch reference `gibbsM SFam lam = Z⁻¹ • exp (-(∑ lam a • SFam
a))` of the retained family at multipliers `lam`, with matrix logarithm
`-(λ·S) - log Z • 1` (`matLog_gibbsM`) and spectral weights `gibbsProb lam`
on the four joint sectors.  The exponent is the modular Hamiltonian `K = λ·S`
and the multipliers are the only parameters: the surface carries no
separate inverse temperature, no modular period, and no energy scale.  The
KMS row of `QFT/StructuralInheritance.lean` (row 5 of the inheritance
matrix: provable on the finite block, refutable on the net and for
continuous public automorphisms) supplies the finite Gibbs-KMS identity
`gibbsState_kms H β` with `β` a free real parameter.  Reading the collar
reference against the Gibbs ledger of
`Thermodynamics/GibbsReferenceEnergyIdentification.lean` therefore lands in
the modular-flow branch: the collar surface fixes `β` only together with a
declared energy scale.  This module proves that dictionary exactly and the
negative exactly.  The identification of the collar reference with the
record reference of the cap first law is declared, never committed by the
collar modules (their reference is a state on a `4 × 4` model algebra, and
the cap first law reads a classical reference on a finite microstate set);
this module reads the spectral weights as that classical reference and says
so at every use.

WHAT IS PROVED.
(1) The collar reference is the corpus Gibbs state of its own modular
Hamiltonian at `β = 1` in modular units: `gibbsProb lam = gibbs (eigval
lam) 1` (`collarRef_eq_gibbs_one`), with the collar partition function
equal to the corpus one at `β = 1` (`collarZ_eq`).  `collarRepairLaw lam`
is the repair-law datum whose reference is the collar weights.
(2) The modular Hamiltonian alone leaves `β` free.  For every finite
nonempty `Ω`, energy `K` and `β ≠ 0`, `gibbs (fun x => K x / β) β = gibbs K
1` (`gibbs_div_beta`), so `gibbs (fun p => eigval lam p / β) β = gibbsProb
lam` for every `β ≠ 0` (`collar_beta_free`).  Every temperature declaration
`d` induces a `CalibrationData` inhabitant for the collar reference with
`energy = eigval lam / d.beta` (`collarCalibration`), and two declarations
with different `beta` both inhabit the row (`collar_beta_not_pinned`): the
collar surface pins no `β`.  Once an energy `H` with two distinct levels is
declared, `β` is determined (`collar_beta_unique_given_energy`, by
`gibbs_beta_injective`), and once `β ≠ 0` is declared, `H` is determined up
to an additive constant (`collar_energy_unique_given_beta`); the dictionary
`eigval lam = β * H + c` is exact in both directions
(`collar_modular_hamiltonian_split`, `collar_gibbs_of_split`).
(3) Exact collar-side ledger identities at `β = 1`: the affine law
`modularEnergy p (gibbsProb lam) = expectation p (eigval lam) + mass p * log
Z` (`collar_affine_law`), the Clausius form of the cap first law
(`collar_cap_firstLaw`), and the collar reference's ledger against itself
equal to its Shannon entropy (`collar_self_ledger_eq_shannon`, from the
general `modularEnergy_self_eq_shannon` for faithful states).  Under a
declared split `H = K / d.beta` the affine law reads with slope `d.beta`
(`collar_affine_law_declared`).
(4) KMS at parameter one.  The collar density is the normalized finite
Gibbs density of `K = ∑ lam a • SFam a` at `β = 1`
(`collarDensity_eq_gibbsDensity_one`), its functional is `gibbsState K 1`
(`collar_functional_eq_gibbsState`), and it satisfies the finite KMS
identity at parameter one for the inner flow of `K`
(`collar_kms_one`), which is the KMS identity at parameter `β` for the
rescaled generator `β⁻¹ • K` (`collar_kms_beta`, through
`gibbsDensity_inv_smul`): the modular flow parameter and the inverse
temperature trade against the energy scale exactly.
(5) Tick dictionary.  `collarLabInverseTemperature a cal d =
labInverseTemperature a cal d.beta`; two distinct ticks give two distinct
laboratory readings of one declared collar `beta`
(`collar_lab_not_forced`), and two distinct declarations give two distinct
readings for one tick (`collar_lab_declaration_sensitive`).

ROWS TOUCHED.  The gravitation-route energy identification (the collar
reference is Gibbs for its own modular Hamiltonian at `β = 1`; the split of
`K` into `β * H` is declared).  The laboratory clock and energy calibration
import (the tick is declared; part (5) converts through it).  The clock-
and-energy calibration register row `FourLawAdequacySurface.CalibrationData`
(inhabited here for every declared `beta`; the row stays declared).  The
source clock and duration row, the physical spacetime attachment row, the
coupled-action row, and the light-signal row are named and untouched.  This
module discharges none of these rows.

NEGATIVES CITED.  The Legendre non-identifiability
(`Variational/RealizedHistoryLegendreNoGo.lean`): realized histories select
no velocity curvature or Legendre map, so every shape is a declared
enrichment; cited at scope because the split `K = β * H` is likewise a
declared enrichment of the collar reference.  The KMS degeneracy and the
public rigidity of `QFT/StructuralInheritance.lean`
(`supportGradedNet_regionalExpectation_kms_degenerate`,
`publicRecord_no_thermal_flow`): the finite KMS identity used here lives on
the private block and selects no temperature.

CONVENTIONS.  The collar exponent is `-(λ·S)` with `λ` the multipliers
(`gibbsM`); the corpus Gibbs exponent is `-β * H x` in the order `gibbs H
β`; the ledger `modularEnergy p tau = ∑ p x * (-log (tau x))` is in nats;
the KMS point is `z = I * β` for `heisenbergFlow`.  The four collar sectors
`Fin 2 × Fin 2` carry the joint spectral projections `specProjP`.

FALSIFIER.  The module fails if the collar weights differ from the corpus
Gibbs state of `eigval lam` at `β = 1`, if `gibbs (K / β) β` differs from
`gibbs K 1` for some `β ≠ 0`, if a temperature declaration fails to inhabit
the calibration row for the collar reference, if a declared two-level
energy fails to pin `β`, if the self-ledger differs from the Shannon
entropy, if the collar density fails the finite KMS identity at parameter
one, or if two ticks give one laboratory reading.

Axiom audit.  No project axiom, no native decision procedure; the guard
lines at the end show at most `propext`, `Classical.choice`, `Quot.sound`.
-/

namespace OPH.CollarTemperatureReading

open OPH.Thermodynamics OPH.InternalEnergyInertia
open OPH.PhysicalCalibrationImport OPH.GibbsReferenceEnergyIdentification
open OPH.QFT

noncomputable section

/-! ## (1) The collar reference as a corpus Gibbs state at `β = 1` -/

/-- The four joint sectors of the collar witness family. -/
abbrev Sector : Type := Fin 2 × Fin 2

/-- The collar partition function is the corpus partition function of the
modular Hamiltonian `eigval lam` at `β = 1`. -/
theorem collarZ_eq (lam : Fin 2 → ℝ) :
    OPH.partitionZ lam = Thermodynamics.partitionZ (eigval lam) 1 := by
  unfold OPH.partitionZ Thermodynamics.partitionZ Thermodynamics.gibbsWeight
  refine Finset.sum_congr rfl fun p _ => ?_
  ring_nf

/-- **The collar reference is Gibbs for its own modular Hamiltonian at
`β = 1`.**  The spectral weights `gibbsProb lam` of `gibbsM SFam lam` are
the corpus state `gibbs (eigval lam) 1`.  Pointer: this is `gibbs_of_neg_log` of the Gibbs reading composed with `gLog_eq_log_gibbsProb`, proved here by direct unfolding. -/
theorem collarRef_eq_gibbs_one (lam : Fin 2 → ℝ) :
    gibbsProb lam = Thermodynamics.gibbs (eigval lam) 1 := by
  funext p
  unfold gibbsProb Thermodynamics.gibbs Thermodynamics.gibbsWeight
  rw [collarZ_eq]
  ring_nf

/-- The repair-law datum whose reference is the collar weights, with the
trivial visible datum.  Reading the spectral weights as the record
reference of the cap first law is declared here. -/
def collarRepairLaw (lam : Fin 2 → ℝ) : RepairLawData Sector Unit where
  ref := gibbsProb lam
  ref_pos := gibbsProb_pos lam
  ref_law := sum_gibbsProb lam
  visible := fun _ => ()

/-! ## (2) The modular Hamiltonian alone leaves `β` free -/

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]

omit [DecidableEq Ω] in
/-- Rescaling the energy by `β⁻¹` and the inverse temperature by `β`
leaves the corpus partition function unchanged. -/
theorem partitionZ_div_beta (K : Ω → ℝ) {β : ℝ} (hβ : β ≠ 0) :
    Thermodynamics.partitionZ (fun x => K x / β) β
      = Thermodynamics.partitionZ K 1 := by
  unfold Thermodynamics.partitionZ Thermodynamics.gibbsWeight
  refine Finset.sum_congr rfl fun x _ => ?_
  congr 1
  field_simp

omit [DecidableEq Ω] in
/-- **`β` is free against the modular Hamiltonian alone.**  For every
`β ≠ 0`, `gibbs (K / β) β = gibbs K 1`. -/
theorem gibbs_div_beta (K : Ω → ℝ) {β : ℝ} (hβ : β ≠ 0) :
    Thermodynamics.gibbs (fun x => K x / β) β = Thermodynamics.gibbs K 1 := by
  funext x
  unfold Thermodynamics.gibbs
  rw [partitionZ_div_beta K hβ]
  congr 1
  unfold Thermodynamics.gibbsWeight
  congr 1
  field_simp

/-- The collar reference is the Gibbs state of `eigval lam / β` at every
`β ≠ 0`. -/
theorem collar_beta_free (lam : Fin 2 → ℝ) {β : ℝ} (hβ : β ≠ 0) :
    Thermodynamics.gibbs (fun p => eigval lam p / β) β = gibbsProb lam := by
  rw [gibbs_div_beta (eigval lam) hβ, collarRef_eq_gibbs_one]

/-- **Inhabitant of the calibration row for every declared `beta`.**  A
temperature declaration `d` reads the collar reference with `energy = eigval
lam / d.beta` and `logZ = log Z`; the `thermal` clause is `matLog_gibbsM`
on the spectral weights (`gLog_eq_log_gibbsProb`).  The split of `K` into
`d.beta * energy` is declared. -/
def collarCalibration (lam : Fin 2 → ℝ) (d : TemperatureDeclaration) :
    CalibrationData (collarRepairLaw lam) where
  energy := fun p => eigval lam p / d.beta
  beta := d.beta
  logZ := Real.log (OPH.partitionZ lam)
  beta_pos := d.beta_pos
  thermal := by
    intro p
    show -Real.log (gibbsProb lam p) = _
    rw [← gLog_eq_log_gibbsProb, gLog]
    have hb := d.beta_pos.ne'
    field_simp
    ring

/-- **The collar surface pins no `β`.**  Two temperature declarations with
different `beta` both inhabit the calibration row for the collar reference,
with the row's `beta` fields equal to the declared values. -/
theorem collar_beta_not_pinned (lam : Fin 2 → ℝ) (d d' : TemperatureDeclaration)
    (hne : d.beta ≠ d'.beta) :
    (collarCalibration lam d).beta ≠ (collarCalibration lam d').beta ∧
      (∀ p, (collarRepairLaw lam).ref p
        = Thermodynamics.gibbs (collarCalibration lam d).energy
            (collarCalibration lam d).beta p) ∧
      (∀ p, (collarRepairLaw lam).ref p
        = Thermodynamics.gibbs (collarCalibration lam d').energy
            (collarCalibration lam d').beta p) :=
  ⟨hne, FourLawSurface.ref_eq_gibbs _ (collarCalibration lam d),
    FourLawSurface.ref_eq_gibbs _ (collarCalibration lam d')⟩

/-- Existential form of the negative: for every `β ≠ 0` some energy makes
the collar reference Gibbs at `β`. -/
theorem collar_gibbs_at_every_beta (lam : Fin 2 → ℝ) {β : ℝ} (hβ : β ≠ 0) :
    ∃ H : Sector → ℝ, gibbsProb lam = Thermodynamics.gibbs H β :=
  ⟨fun p => eigval lam p / β, (collar_beta_free lam hβ).symm⟩

/-- **A declared two-level energy pins `β`.**  If the collar reference is
Gibbs for one declared `H` with two distinct levels at `β` and at `β'`, then
`β = β'` (`gibbs_beta_injective`). -/
theorem collar_beta_unique_given_energy (lam : Fin 2 → ℝ) (H : Sector → ℝ)
    (i j : Sector) (hij : H i ≠ H j) {β β' : ℝ}
    (h : gibbsProb lam = Thermodynamics.gibbs H β)
    (h' : gibbsProb lam = Thermodynamics.gibbs H β') : β = β' :=
  gibbs_beta_injective H β β' i j hij (h.symm.trans h')

/-- **A declared nonzero `β` pins the energy up to a constant.** -/
theorem collar_energy_unique_given_beta (lam : Fin 2 → ℝ) (H H' : Sector → ℝ)
    {β : ℝ} (hβ : β ≠ 0)
    (h : gibbsProb lam = Thermodynamics.gibbs H β)
    (h' : gibbsProb lam = Thermodynamics.gibbs H' β) :
    ∃ c : ℝ, ∀ p, H' p = H p + c :=
  gibbs_energy_unique_up_to_const H H' β hβ (h.symm.trans h')

/-- **Dictionary, forward.**  If the collar reference is Gibbs for `H` at
`β ≠ 0`, the modular Hamiltonian splits as `eigval lam = β * H + c`. -/
theorem collar_modular_hamiltonian_split (lam : Fin 2 → ℝ) (H : Sector → ℝ)
    {β : ℝ} (hβ : β ≠ 0) (h : gibbsProb lam = Thermodynamics.gibbs H β) :
    ∃ c : ℝ, ∀ p, eigval lam p = β * H p + c := by
  obtain ⟨c, hc⟩ := collar_energy_unique_given_beta lam
    (fun p => eigval lam p / β) H hβ (collar_beta_free lam hβ).symm h
  refine ⟨-(β * c), fun p => ?_⟩
  have := hc p
  simp only at this
  rw [this]
  field_simp
  ring

/-- **Dictionary, backward.**  If `eigval lam = β * H + c` with `β ≠ 0`,
the collar reference is Gibbs for `H` at `β`. -/
theorem collar_gibbs_of_split (lam : Fin 2 → ℝ) (H : Sector → ℝ) {β c : ℝ}
    (hβ : β ≠ 0) (h : ∀ p, eigval lam p = β * H p + c) :
    gibbsProb lam = Thermodynamics.gibbs H β := by
  rw [← collar_beta_free lam hβ]
  have hH : (fun p => eigval lam p / β) = fun p => H p + c / β := by
    funext p
    rw [h p]
    field_simp
  rw [hH, gibbs_add_const]

/-! ## (3) Exact collar-side ledger identities at `β = 1` -/

omit [DecidableEq Ω] in
/-- **The self-ledger of a faithful state is its Shannon entropy.**
`modularEnergy p p = shannon p` when every `p x` is positive.  Pointer: `FiniteConditionalRepair.kl_eq_energy_sub_shannon` at `tau = p` with `kl p p = 0` gives the same identity. -/
theorem modularEnergy_self_eq_shannon (p : Ω → ℝ) (hp : ∀ x, 0 < p x) :
    modularEnergy p p = shannon p := by
  unfold modularEnergy shannon klTerm
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun x _ => ?_
  rw [if_neg (hp x).ne', div_one]
  ring

/-- The collar reference's ledger against itself is its Shannon entropy. -/
theorem collar_self_ledger_eq_shannon (lam : Fin 2 → ℝ) :
    modularEnergy (gibbsProb lam) (gibbsProb lam) = shannon (gibbsProb lam) :=
  modularEnergy_self_eq_shannon (gibbsProb lam) (gibbsProb_pos lam)

/-- **Collar affine law at `β = 1`.**  `modularEnergy p (gibbsProb lam) =
⟨eigval lam⟩_p + (∑ p) * log Z`, the instance of `modularEnergy_gibbs`. -/
theorem collar_affine_law (lam : Fin 2 → ℝ) (p : Sector → ℝ) :
    modularEnergy p (gibbsProb lam)
      = expectation p (eigval lam) + mass p * Real.log (OPH.partitionZ lam) := by
  rw [collarRef_eq_gibbs_one, modularEnergy_gibbs, collarZ_eq, one_mul]

/-- **Collar Clausius form at `β = 1`.**  For a normalized `p`,
`S(p) - S(ω) = ⟨K⟩_p - ⟨K⟩_ω - D(p ‖ ω)` with `K = eigval lam`. -/
theorem collar_cap_firstLaw (lam : Fin 2 → ℝ) (p : Sector → ℝ)
    (hp : ∑ x, p x = 1) :
    shannon p - shannon (gibbsProb lam)
      = (expectation p (eigval lam) - expectation (gibbsProb lam) (eigval lam))
        - kl p (gibbsProb lam) := by
  rw [collarRef_eq_gibbs_one, cap_firstLaw_gibbs (eigval lam) 1 p hp, one_mul]

/-- **Collar affine law under a declared split.**  With `energy = eigval lam
/ d.beta` the slope is the declared `d.beta`; the content is the same
identity as `collar_affine_law`. -/
theorem collar_affine_law_declared (lam : Fin 2 → ℝ) (d : TemperatureDeclaration)
    (p : Sector → ℝ) :
    modularEnergy p (collarRepairLaw lam).ref
      = (collarCalibration lam d).beta
          * expectation p (collarCalibration lam d).energy
        + mass p * (collarCalibration lam d).logZ := by
  show modularEnergy p (gibbsProb lam)
    = d.beta * expectation p (fun q => eigval lam q / d.beta)
      + mass p * Real.log (OPH.partitionZ lam)
  rw [← collar_beta_free lam d.beta_pos.ne', modularEnergy_gibbs,
    partitionZ_div_beta (eigval lam) d.beta_pos.ne', ← collarZ_eq]

/-! ## (4) KMS at parameter one -/

/-- The collar modular Hamiltonian as a matrix: `K = ∑ lam a • SFam a`. -/
def collarK (lam : Fin 2 → ℝ) : CollarC :=
  ∑ a : Fin 2, ((lam a : ℝ) : ℂ) • SFam a

/-- The collar density is the normalized finite Gibbs density of `collarK
lam` at `β = 1`, in the objects of `QFT/StructuralInheritance.lean`. -/
theorem collarDensity_eq_gibbsDensity_one (lam : Fin 2 → ℝ) :
    gibbsM SFam lam
      = (gibbsPartition (collarK lam) 1)⁻¹ • gibbsDensity (collarK lam) 1 := by
  unfold gibbsM gibbsPartition gibbsDensity collarK
  rw [Complex.ofReal_one, one_smul]

/-- The collar functional `a ↦ Tr (ω a)` is `gibbsState (collarK lam) 1`. -/
theorem collar_functional_eq_gibbsState (lam : Fin 2 → ℝ) (a : CollarC) :
    (gibbsM SFam lam * a).trace = gibbsState (collarK lam) 1 a := by
  rw [collarDensity_eq_gibbsDensity_one, gibbsState, Matrix.smul_mul,
    Matrix.trace_smul, smul_eq_mul, inv_mul_eq_div]

/-- **Collar KMS at parameter one.**  The collar functional satisfies the
finite KMS identity `ω (a σ_{i}(b)) = ω (b a)` for the inner flow of
`collarK lam` at `β = 1`. -/
theorem collar_kms_one (lam : Fin 2 → ℝ) (a b : CollarC) :
    (gibbsM SFam lam
        * (a * heisenbergFlow (collarK lam) b (Complex.I * ((1 : ℝ) : ℂ)))).trace
      = (gibbsM SFam lam * (b * a)).trace := by
  rw [collar_functional_eq_gibbsState, collar_functional_eq_gibbsState]
  exact gibbsState_kms (collarK lam) 1 a b

/-- Rescaling the generator by `β⁻¹` and the parameter by `β` leaves the
finite Gibbs density unchanged. -/
theorem gibbsDensity_inv_smul {d : Type*} [Fintype d] [DecidableEq d]
    (H : Matrix d d ℂ) {β : ℝ} (hβ : β ≠ 0) :
    gibbsDensity (((β : ℂ)⁻¹) • H) β = gibbsDensity H 1 := by
  unfold gibbsDensity
  rw [smul_smul, Complex.ofReal_one, one_smul,
    mul_inv_cancel₀ (Complex.ofReal_ne_zero.mpr hβ), one_smul]

/-- **Collar KMS at parameter `β` for the rescaled generator.**  The
modular parameter and the inverse temperature trade against the energy
scale: the collar functional is KMS at `β` for the flow of `β⁻¹ • collarK
lam`, for every `β ≠ 0`. -/
theorem collar_kms_beta (lam : Fin 2 → ℝ) {β : ℝ} (hβ : β ≠ 0) (a b : CollarC) :
    (gibbsM SFam lam
        * (a * heisenbergFlow (((β : ℂ)⁻¹) • collarK lam) b (Complex.I * β))).trace
      = (gibbsM SFam lam * (b * a)).trace := by
  have hstate : ∀ c : CollarC,
      (gibbsM SFam lam * c).trace = gibbsState (((β : ℂ)⁻¹) • collarK lam) β c := by
    intro c
    rw [collar_functional_eq_gibbsState, gibbsState, gibbsState, gibbsPartition,
      gibbsPartition, gibbsDensity_inv_smul (collarK lam) hβ]
  rw [hstate, hstate]
  exact gibbsState_kms _ β a b

/-! ## (5) Tick dictionary -/

/-- Laboratory inverse temperature of a declared collar `beta` for a
declared tick: the composition with `labInverseTemperature`. -/
def collarLabInverseTemperature (a : SIAnchors) (cal : ClockCalibration)
    (d : TemperatureDeclaration) : ℝ :=
  labInverseTemperature a cal d.beta

/-- **Non-forcing.**  Two distinct declared ticks give two distinct
laboratory readings of one declared collar `beta`. -/
theorem collar_lab_not_forced (a : SIAnchors) (cal cal' : ClockCalibration)
    (hne : cal.tau ≠ cal'.tau) (d : TemperatureDeclaration) :
    collarLabInverseTemperature a cal d ≠ collarLabInverseTemperature a cal' d :=
  labInverseTemperature_not_forced a cal cal' hne d.beta_pos.ne'

/-- Two distinct declarations give two distinct laboratory readings for one
tick: the reading carries the declared `beta`, which the collar surface does
not pin. -/
theorem collar_lab_declaration_sensitive (a : SIAnchors) (cal : ClockCalibration)
    (d d' : TemperatureDeclaration) (hne : d.beta ≠ d'.beta) :
    collarLabInverseTemperature a cal d ≠ collarLabInverseTemperature a cal d' := by
  unfold collarLabInverseTemperature labInverseTemperature
  intro h
  exact hne ((div_left_inj' (hbarOverTau_pos a cal).ne').mp h)

end

end OPH.CollarTemperatureReading

#print axioms OPH.CollarTemperatureReading.collarZ_eq
#print axioms OPH.CollarTemperatureReading.collarRef_eq_gibbs_one
#print axioms OPH.CollarTemperatureReading.gibbs_div_beta
#print axioms OPH.CollarTemperatureReading.collar_beta_free
#print axioms OPH.CollarTemperatureReading.collarCalibration
#print axioms OPH.CollarTemperatureReading.collar_beta_not_pinned
#print axioms OPH.CollarTemperatureReading.collar_gibbs_at_every_beta
#print axioms OPH.CollarTemperatureReading.collar_beta_unique_given_energy
#print axioms OPH.CollarTemperatureReading.collar_energy_unique_given_beta
#print axioms OPH.CollarTemperatureReading.collar_modular_hamiltonian_split
#print axioms OPH.CollarTemperatureReading.collar_gibbs_of_split
#print axioms OPH.CollarTemperatureReading.modularEnergy_self_eq_shannon
#print axioms OPH.CollarTemperatureReading.collar_self_ledger_eq_shannon
#print axioms OPH.CollarTemperatureReading.collar_affine_law
#print axioms OPH.CollarTemperatureReading.collar_cap_firstLaw
#print axioms OPH.CollarTemperatureReading.collar_affine_law_declared
#print axioms OPH.CollarTemperatureReading.collarDensity_eq_gibbsDensity_one
#print axioms OPH.CollarTemperatureReading.collar_functional_eq_gibbsState
#print axioms OPH.CollarTemperatureReading.collar_kms_one
#print axioms OPH.CollarTemperatureReading.gibbsDensity_inv_smul
#print axioms OPH.CollarTemperatureReading.collar_kms_beta
#print axioms OPH.CollarTemperatureReading.collar_lab_not_forced
#print axioms OPH.CollarTemperatureReading.collar_lab_declaration_sensitive
