import ObserverPatchHolography.EinsteinBranch.DarkSector

/-!
# Collar-scale saturation for the dark-sector screening interface

The dark-sector module bounds the anomalous collar stress by the collar
recovery defect and proves the screening interface: the stress tends to zero
along any family whose defect tends to zero at fixed radius and bounded
collar constant.  This module supplies the rate side of that interface.  The
collar-recovery claim of the Einstein branch bounds the defect at depth
`delta` by an exponential envelope
`kappa * B * exp (-(delta * rho^(1/3) - r0) / zeta)`, where `rho` is the
local record density, so `rho^(1/3)` is the inverse UV record spacing and
`delta * rho^(1/3)` is the collar depth counted in UV cells.

* `RecoveryEnvelope` types that envelope: prefactor `kappa >= 0`, boundary
  size `B >= 1` in UV cells, offset `r0 >= 0`, decay constant `zeta > 0`,
  record density `rho > 0`, and a defect function `eta` obeying the bound at
  every nonnegative depth.
* Depth saturation: for a fixed envelope the bound and the defect vanish as
  the depth grows (`envelopeBound_tendsto_zero_of_depth`,
  `RecoveryEnvelope.eta_tendsto_zero_of_depth`), and composing with the
  committed screening interface `anomalousStress_tendsto_zero` switches the
  anomalous stress off along any remainder family built on the envelope
  (`eta_and_stress_tendsto_zero_of_depth`).
* Density saturation: at fixed positive depth the bound and the composed
  stress vanish as the record density grows
  (`envelopeBound_tendsto_zero_of_density`,
  `eta_and_stress_tendsto_zero_of_density`).  Where records are dense, a
  source generates no persistent anomalous charge of its own.
* Quantitative rate: any remainder whose defect sits on an envelope at depth
  `delta` has stress bounded by
  `15 C / (8 pi^2 ell^4) * kappa * B * exp (-(delta * rho^(1/3) - r0) / zeta)`
  (`anomalousStress_exponential_rate`, `remainder_stress_rate`).
* Settled-cut threshold: the envelope bound is below `eps` exactly when the
  UV depth `delta * rho^(1/3)` exceeds `r0 + zeta * log (kappa * B / eps)`
  (`envelopeBound_lt_iff`).
* Delimitation in both directions.  The identically zero defect satisfies
  every envelope hypothesis and gives zero stress
  (`recoveryEnvelope_admits_zero_defect`,
  `zeroDefectEnvelope_stress_eq_zero`): the envelope never forces a nonzero
  defect, so galactic-scale sourcing is a premise, not a consequence.  A
  concrete sparsely recorded instance sits at the bound with an order-one
  defect at UV depth equal to the offset
  (`recoveryEnvelope_admits_order_one_defect`,
  `sparseCutRemainder_stress_eq`, `sparseCutRemainder_stress_pos`): the same
  hypotheses admit an order-one defect on sparsely recorded cuts.

What is not proved here: that physical collars satisfy this envelope is the
conditional mixing branch of the collar-recovery claim
(OPH-GR-D2-COLLAR-RECOVERY), taken here as a hypothesis carried by the
structure fields; no physical value of the record density, the prefactor, or
the decay constant is asserted; and nothing here proves that the defect is
nonzero on any galactic cut.  The collar-scale premise of the dark-matter
paper, that the carrying collars are galactic-scale settled cuts, is
supported by the density-saturation clause only in the stated conditional
form.
-/

namespace OPH.EinsteinBranch

open Filter
open scoped Topology

noncomputable section

/-! ## The exponential recovery envelope -/

/-- Exponential recovery envelope for the collar defect.  `kappa` is the
mixing prefactor, `B` the boundary size in UV cells, `r0` the rate offset,
`zeta` the decay constant in UV cells, and `rho` the local record density,
so `rho ^ (1/3)` is the inverse UV record spacing and a collar of depth
`delta` spans `delta * rho ^ (1/3)` UV cells.  `eta delta` is the collar
recovery defect at depth `delta`, bounded by the envelope at every
nonnegative depth.  That physical collars satisfy this bound is the
conditional mixing branch of the collar-recovery claim, carried here as a
premise. -/
structure RecoveryEnvelope where
  kappa : ℝ
  kappa_nonneg : 0 ≤ kappa
  B : ℝ
  B_ge_one : 1 ≤ B
  r0 : ℝ
  r0_nonneg : 0 ≤ r0
  zeta : ℝ
  zeta_pos : 0 < zeta
  rho : ℝ
  rho_pos : 0 < rho
  eta : ℝ → ℝ
  eta_nonneg : ∀ δ, 0 ≤ eta δ
  eta_le : ∀ δ, 0 ≤ δ →
    eta δ ≤ kappa * B * Real.exp (-(δ * rho ^ ((1 : ℝ) / 3) - r0) / zeta)

/-- The envelope bound at depth `δ`. -/
def envelopeBound (E : RecoveryEnvelope) (δ : ℝ) : ℝ :=
  E.kappa * E.B * Real.exp (-(δ * E.rho ^ ((1 : ℝ) / 3) - E.r0) / E.zeta)

/-- The anomaly remainder carried by a collar of depth `δ` whose defect sits
on the envelope `E`, with collar radius `ℓ`, collar constant `c`, and
anomalous energy `energy` obeying the collar bound. -/
def RecoveryEnvelope.remainder (E : RecoveryEnvelope) (ℓ c δ energy : ℝ)
    (hℓ : 0 < ℓ) (hc : 0 ≤ c) (hE : |energy| ≤ c * E.eta δ) :
    AnomalyRemainder where
  ell := ℓ
  ell_pos := hℓ
  C := c
  C_nonneg := hc
  eta := E.eta δ
  eta_nonneg := E.eta_nonneg δ
  anomalousEnergy := energy
  bound := hE

/-! ## Depth saturation -/

/-- The exponential envelope in the UV depth `m` vanishes as `m` grows. -/
theorem expRate_tendsto_zero (kappa B r0 zeta : ℝ) (hzeta : 0 < zeta) :
    Tendsto (fun m : ℝ => kappa * B * Real.exp (-(m - r0) / zeta))
      atTop (𝓝 0) := by
  have h1 : Tendsto (fun m : ℝ => (m - r0) / zeta) atTop atTop := by
    have h0 : Tendsto (fun m : ℝ => m - r0) atTop atTop := by
      simpa [sub_eq_add_neg] using
        tendsto_atTop_add_const_right atTop (-r0) (tendsto_id (α := ℝ))
    exact h0.atTop_div_const hzeta
  have h2 : Tendsto (fun m : ℝ => Real.exp (-(m - r0) / zeta)) atTop (𝓝 0) := by
    have h3 := Real.tendsto_exp_atBot.comp (tendsto_neg_atTop_atBot.comp h1)
    simpa only [Function.comp_def, neg_div] using h3
  have h4 := h2.const_mul (kappa * B)
  rw [mul_zero] at h4
  exact h4

/-- For a fixed envelope the bound vanishes as the depth grows. -/
theorem envelopeBound_tendsto_zero_of_depth (E : RecoveryEnvelope) :
    Tendsto (fun δ : ℝ => envelopeBound E δ) atTop (𝓝 0) := by
  have hc : 0 < E.rho ^ ((1 : ℝ) / 3) := Real.rpow_pos_of_pos E.rho_pos _
  have hm : Tendsto (fun δ : ℝ => δ * E.rho ^ ((1 : ℝ) / 3)) atTop atTop :=
    tendsto_id.atTop_mul_const hc
  have h :=
    (expRate_tendsto_zero E.kappa E.B E.r0 E.zeta E.zeta_pos).comp hm
  simpa only [envelopeBound, Function.comp_def] using h

/-- The defect of an envelope vanishes along any nonnegative depth sequence
tending to infinity. -/
theorem RecoveryEnvelope.eta_tendsto_zero_of_depth (E : RecoveryEnvelope)
    (δseq : ℕ → ℝ) (hnn : ∀ n, 0 ≤ δseq n)
    (hδ : Tendsto δseq atTop atTop) :
    Tendsto (fun n => E.eta (δseq n)) atTop (𝓝 0) := by
  have hb : Tendsto (fun n => envelopeBound E (δseq n)) atTop (𝓝 0) := by
    have h := (envelopeBound_tendsto_zero_of_depth E).comp hδ
    simpa only [Function.comp_def] using h
  have hle : ∀ n, E.eta (δseq n) ≤ envelopeBound E (δseq n) := fun n =>
    E.eta_le _ (hnn n)
  exact squeeze_zero (fun n => E.eta_nonneg _) hle hb

/-- Depth saturation composed with the screening interface: along a family
of anomaly remainders built on one envelope at depths tending to infinity,
with fixed collar radius and bounded collar constant, the recovery defect
and the anomalous stress both tend to zero. -/
theorem eta_and_stress_tendsto_zero_of_depth (E : RecoveryEnvelope)
    (δseq : ℕ → ℝ) (A : ℕ → AnomalyRemainder) (ℓ Cmax : ℝ)
    (hnn : ∀ n, 0 ≤ δseq n) (hδ : Tendsto δseq atTop atTop)
    (hℓ : ∀ n, (A n).ell = ℓ) (hC : ∀ n, (A n).C ≤ Cmax)
    (heta : ∀ n, (A n).eta = E.eta (δseq n)) :
    Tendsto (fun n => (A n).eta) atTop (𝓝 0) ∧
      Tendsto (fun n => anomalousStress (A n)) atTop (𝓝 0) := by
  have h0 : Tendsto (fun n => (A n).eta) atTop (𝓝 0) := by
    have h := E.eta_tendsto_zero_of_depth δseq hnn hδ
    simpa only [heta] using h
  exact ⟨h0, anomalousStress_tendsto_zero A ℓ Cmax hℓ hC h0⟩

/-! ## Density saturation -/

/-- At fixed positive depth the envelope bound vanishes as the record
density grows. -/
theorem envelopeBound_tendsto_zero_of_density (kappa B r0 zeta δ : ℝ)
    (hzeta : 0 < zeta) (hδ : 0 < δ) :
    Tendsto
      (fun rho : ℝ =>
        kappa * B * Real.exp (-(δ * rho ^ ((1 : ℝ) / 3) - r0) / zeta))
      atTop (𝓝 0) := by
  have h13 : (0 : ℝ) < 1 / 3 := by norm_num
  have hm : Tendsto (fun rho : ℝ => δ * rho ^ ((1 : ℝ) / 3)) atTop atTop :=
    (tendsto_rpow_atTop h13).const_mul_atTop hδ
  have h := (expRate_tendsto_zero kappa B r0 zeta hzeta).comp hm
  simpa only [Function.comp_def] using h

/-- Density saturation composed with the screening interface: along a family
of envelopes with common prefactor, boundary size, offset, and decay
constant whose record densities tend to infinity, the defects at any fixed
positive depth and the anomalous stress of any remainder family carried on
them both tend to zero.  This is the compact-source clause: where records
are dense, a source generates no persistent anomalous charge of its own. -/
theorem eta_and_stress_tendsto_zero_of_density (E : ℕ → RecoveryEnvelope)
    (A : ℕ → AnomalyRemainder) (δ ℓ Cmax kappa B r0 zeta : ℝ) (hδ : 0 < δ)
    (hkappa : ∀ n, (E n).kappa = kappa) (hB : ∀ n, (E n).B = B)
    (hr0 : ∀ n, (E n).r0 = r0) (hzeta : ∀ n, (E n).zeta = zeta)
    (hrho : Tendsto (fun n => (E n).rho) atTop atTop)
    (hℓ : ∀ n, (A n).ell = ℓ) (hC : ∀ n, (A n).C ≤ Cmax)
    (heta : ∀ n, (A n).eta = (E n).eta δ) :
    Tendsto (fun n => (A n).eta) atTop (𝓝 0) ∧
      Tendsto (fun n => anomalousStress (A n)) atTop (𝓝 0) := by
  have hzeta0 : 0 < zeta := by
    rw [← hzeta 0]; exact (E 0).zeta_pos
  have hg :
      Tendsto
        (fun n =>
          kappa * B *
            Real.exp (-(δ * (E n).rho ^ ((1 : ℝ) / 3) - r0) / zeta))
        atTop (𝓝 0) := by
    have h :=
      (envelopeBound_tendsto_zero_of_density kappa B r0 zeta δ hzeta0
        hδ).comp hrho
    simpa only [Function.comp_def] using h
  have hle : ∀ n,
      (A n).eta ≤
        kappa * B *
          Real.exp (-(δ * (E n).rho ^ ((1 : ℝ) / 3) - r0) / zeta) := by
    intro n
    rw [heta n]
    have h := (E n).eta_le δ (le_of_lt hδ)
    rwa [hkappa n, hB n, hr0 n, hzeta n] at h
  have h0 : Tendsto (fun n => (A n).eta) atTop (𝓝 0) :=
    squeeze_zero (fun n => (A n).eta_nonneg) hle hg
  exact ⟨h0, anomalousStress_tendsto_zero A ℓ Cmax hℓ hC h0⟩

/-! ## Quantitative rate in the record density -/

/-- Exponential rate on the stress itself: any remainder whose defect sits
on the envelope `E` at nonnegative depth `δ` has anomalous stress bounded by
the small-ball kernel factor times the envelope. -/
theorem anomalousStress_exponential_rate (E : RecoveryEnvelope)
    (A : AnomalyRemainder) (δ : ℝ) (hδ : 0 ≤ δ) (heta : A.eta = E.eta δ) :
    |anomalousStress A| ≤
      15 * A.C / (8 * Real.pi ^ 2 * A.ell ^ 4) *
        (E.kappa * E.B *
          Real.exp (-(δ * E.rho ^ ((1 : ℝ) / 3) - E.r0) / E.zeta)) := by
  have h1 := anomalousStress_abs_le A
  have h2 : A.eta ≤
      E.kappa * E.B *
        Real.exp (-(δ * E.rho ^ ((1 : ℝ) / 3) - E.r0) / E.zeta) := by
    rw [heta]; exact E.eta_le δ hδ
  have hk : 0 ≤ 15 * A.C / (8 * Real.pi ^ 2 * A.ell ^ 4) := by
    have hπ := Real.pi_pos
    have hℓ := A.ell_pos
    have hc := A.C_nonneg
    positivity
  calc |anomalousStress A|
      ≤ 15 * A.C * A.eta / (8 * Real.pi ^ 2 * A.ell ^ 4) := h1
    _ = 15 * A.C / (8 * Real.pi ^ 2 * A.ell ^ 4) * A.eta := by ring
    _ ≤ 15 * A.C / (8 * Real.pi ^ 2 * A.ell ^ 4) *
        (E.kappa * E.B *
          Real.exp (-(δ * E.rho ^ ((1 : ℝ) / 3) - E.r0) / E.zeta)) :=
        mul_le_mul_of_nonneg_left h2 hk

/-- The rate receipt for the composed instance built by
`RecoveryEnvelope.remainder`. -/
theorem remainder_stress_rate (E : RecoveryEnvelope) (ℓ c δ energy : ℝ)
    (hℓ : 0 < ℓ) (hc : 0 ≤ c) (hδ : 0 ≤ δ)
    (hE : |energy| ≤ c * E.eta δ) :
    |anomalousStress (E.remainder ℓ c δ energy hℓ hc hE)| ≤
      15 * c / (8 * Real.pi ^ 2 * ℓ ^ 4) *
        (E.kappa * E.B *
          Real.exp (-(δ * E.rho ^ ((1 : ℝ) / 3) - E.r0) / E.zeta)) :=
  anomalousStress_exponential_rate E _ δ hδ rfl

/-! ## The settled-cut threshold -/

/-- Settled-cut threshold: the envelope bound falls below `ε` exactly when
the UV depth `δ * rho ^ (1/3)` exceeds `r0 + zeta * log (kappa * B / ε)`.
The statement holds for every positive `ε`; for `ε ≤ kappa * B` the
threshold depth is nonnegative and separates saturated collars from settled
cuts. -/
theorem envelopeBound_lt_iff (E : RecoveryEnvelope) (δ ε : ℝ)
    (hkB : 0 < E.kappa * E.B) (hε : 0 < ε) :
    envelopeBound E δ < ε ↔
      E.r0 + E.zeta * Real.log (E.kappa * E.B / ε) <
        δ * E.rho ^ ((1 : ℝ) / 3) := by
  unfold envelopeBound
  set m := δ * E.rho ^ ((1 : ℝ) / 3) with hm
  rw [mul_comm (E.kappa * E.B) (Real.exp (-(m - E.r0) / E.zeta)),
    ← lt_div_iff₀ hkB, ← Real.lt_log_iff_exp_lt (div_pos hε hkB),
    div_lt_iff₀ E.zeta_pos,
    Real.log_div (ne_of_gt hε) (ne_of_gt hkB),
    Real.log_div (ne_of_gt hkB) (ne_of_gt hε)]
  constructor <;> intro h <;> nlinarith [h]

/-! ## Delimitation: the envelope forces nothing -/

/-- Envelope instance with identically zero defect: fully settled records. -/
def zeroDefectEnvelope : RecoveryEnvelope where
  kappa := 3
  kappa_nonneg := by norm_num
  B := 2
  B_ge_one := one_le_two
  r0 := 1
  r0_nonneg := zero_le_one
  zeta := 2
  zeta_pos := two_pos
  rho := 8
  rho_pos := by norm_num
  eta := fun _ => 0
  eta_nonneg := fun _ => le_refl 0
  eta_le := fun δ _ => by
    have h : (0 : ℝ) ≤
        3 * 2 * Real.exp (-(δ * (8 : ℝ) ^ ((1 : ℝ) / 3) - 1) / 2) := by
      positivity
    exact h

/-- The envelope hypotheses admit the identically zero defect: the envelope
never forces a nonzero defect, so galactic-scale sourcing is a premise, not
a consequence. -/
theorem recoveryEnvelope_admits_zero_defect :
    ∃ E : RecoveryEnvelope, ∀ δ : ℝ, E.eta δ = 0 :=
  ⟨zeroDefectEnvelope, fun _ => rfl⟩

/-- Every remainder carried on the zero-defect envelope has zero anomalous
stress. -/
theorem zeroDefectEnvelope_stress_eq_zero (ℓ c δ energy : ℝ)
    (hℓ : 0 < ℓ) (hc : 0 ≤ c)
    (hE : |energy| ≤ c * zeroDefectEnvelope.eta δ) :
    anomalousStress (zeroDefectEnvelope.remainder ℓ c δ energy hℓ hc hE) =
      0 :=
  anomalousStress_eq_zero_of_recovered _ rfl

/-! ## Delimitation: sparse records admit an order-one defect -/

/-- Envelope instance for a sparsely recorded cut: unit record density,
unit decay constant, offset one, and the defect sitting exactly at the
bound value `2 * exp (-(δ - 1))`. -/
def sparseCutEnvelope : RecoveryEnvelope where
  kappa := 1
  kappa_nonneg := zero_le_one
  B := 2
  B_ge_one := one_le_two
  r0 := 1
  r0_nonneg := zero_le_one
  zeta := 1
  zeta_pos := one_pos
  rho := 1
  rho_pos := one_pos
  eta := fun δ => 2 * Real.exp (-(δ - 1))
  eta_nonneg := fun δ => by
    have h : (0 : ℝ) ≤ 2 * Real.exp (-(δ - 1)) := by positivity
    exact h
  eta_le := fun δ _ => by simp [Real.one_rpow]

/-- On the sparse instance the defect at depth one is exactly `2`. -/
theorem sparseCutEnvelope_eta_one : sparseCutEnvelope.eta 1 = 2 := by
  show 2 * Real.exp (-((1 : ℝ) - 1)) = 2
  norm_num

/-- On the sparse instance the UV depth at `δ = 1` equals the offset. -/
theorem sparseCutEnvelope_depth_eq_offset :
    (1 : ℝ) * sparseCutEnvelope.rho ^ ((1 : ℝ) / 3) =
      sparseCutEnvelope.r0 := by
  show (1 : ℝ) * (1 : ℝ) ^ ((1 : ℝ) / 3) = 1
  rw [Real.one_rpow, mul_one]

/-- The envelope hypotheses admit an order-one defect at a depth whose UV
count does not exceed the offset: sparsely recorded cuts are not screened by
the envelope alone. -/
theorem recoveryEnvelope_admits_order_one_defect :
    ∃ (E : RecoveryEnvelope) (δ : ℝ), 0 ≤ δ ∧
      δ * E.rho ^ ((1 : ℝ) / 3) ≤ E.r0 ∧ (1 : ℝ) ≤ E.eta δ := by
  refine ⟨sparseCutEnvelope, 1, zero_le_one,
    le_of_eq sparseCutEnvelope_depth_eq_offset, ?_⟩
  rw [sparseCutEnvelope_eta_one]
  norm_num

/-- Concrete remainder on the sparse instance: unit collar radius, unit
collar constant, depth one, anomalous energy `2` at the collar bound. -/
def sparseCutRemainder : AnomalyRemainder :=
  sparseCutEnvelope.remainder 1 1 1 2 one_pos zero_le_one
    (by rw [sparseCutEnvelope_eta_one, abs_two]; norm_num)

/-- The sparse remainder carries the exact nonzero stress
`15 / (4 pi^2)`. -/
theorem sparseCutRemainder_stress_eq :
    anomalousStress sparseCutRemainder = 15 / (4 * Real.pi ^ 2) := by
  show 15 / (8 * Real.pi ^ 2 * (1 : ℝ) ^ 4) * 2 = 15 / (4 * Real.pi ^ 2)
  have h8 : (8 : ℝ) * Real.pi ^ 2 ≠ 0 := by positivity
  have h4 : (4 : ℝ) * Real.pi ^ 2 ≠ 0 := by positivity
  rw [one_pow, mul_one, div_mul_eq_mul_div, div_eq_div_iff h8 h4]
  ring

/-- The sparse remainder has strictly positive anomalous stress: the same
envelope hypotheses that screen dense regions admit a nonzero dark source on
sparsely recorded cuts. -/
theorem sparseCutRemainder_stress_pos :
    0 < anomalousStress sparseCutRemainder := by
  rw [sparseCutRemainder_stress_eq]
  positivity

/-! ## Per-theorem axiom audit -/

#print axioms expRate_tendsto_zero
#print axioms envelopeBound_tendsto_zero_of_depth
#print axioms RecoveryEnvelope.eta_tendsto_zero_of_depth
#print axioms eta_and_stress_tendsto_zero_of_depth
#print axioms envelopeBound_tendsto_zero_of_density
#print axioms eta_and_stress_tendsto_zero_of_density
#print axioms anomalousStress_exponential_rate
#print axioms remainder_stress_rate
#print axioms envelopeBound_lt_iff
#print axioms recoveryEnvelope_admits_zero_defect
#print axioms zeroDefectEnvelope_stress_eq_zero
#print axioms sparseCutEnvelope_eta_one
#print axioms sparseCutEnvelope_depth_eq_offset
#print axioms recoveryEnvelope_admits_order_one_defect
#print axioms sparseCutRemainder_stress_eq
#print axioms sparseCutRemainder_stress_pos

end

end OPH.EinsteinBranch
