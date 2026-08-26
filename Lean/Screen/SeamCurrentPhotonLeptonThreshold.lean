import SeamCurrentFreePhotonLift

namespace OPH.SeamCurrentPhotonLeptonThreshold

abbrev ThresholdVec3 := OPH.SeamCurrentFreePhotonLift.FreePhotonVec3
abbrev ThresholdEuclideanVec3 :=
  OPH.PrimitivePortFrameQuotient.EuclideanVec3

/-!
# Conditional photon-lepton threshold algebra for FZ-12

This file isolates the leading ultrarelativistic kinematics of
`gamma gamma -> e+ e-` in one preferred frame.  Energy and momentum are
composed additively.  The hard photon and the two charged leptons may carry
independent dimension-six dispersion coefficients.  The soft background
photon is kept at its Lorentz-invariant leading order.

The derivation retains the soft energy in the leading outgoing energy sum and
drops it inside terms suppressed by a mass or dispersion coefficient.
This is the standard leading threshold expansion, not an exact interaction
theorem.  The exact FZ-12 spatial-symbol remainder is proved separately in
`SeamCurrentEdge30Remainder`.

Nothing here selects the electron or positron action, proves additive physical
energy-momentum conservation, supplies a pair-production vertex or cross
section, or identifies the conditional FZ-12 oscillator with a photon.  The
theorems expose the lepton degeneracy that any physical UHE comparison must
resolve.
-/

/-- Leading hard-photon momentum for
`E^2 = p^2 + delta E^4`, expanded to first order in `delta`. -/
noncomputable def photonMomentumLeading (delta energy : ℝ) : ℝ :=
  energy - delta * energy ^ 3 / 2

/-- Leading outgoing charged-lepton momentum.  `share` is its fraction of the
hard energy.  The soft energy is retained only in the leading energy sum,
which is the order needed to recover the head-on `4 epsilon E` term. -/
noncomputable def leptonMomentumLeading
    (delta mass share hardEnergy softEnergy : ℝ) : ℝ :=
  share * (hardEnergy + softEnergy)
    - mass ^ 2 / (2 * share * hardEnergy)
    - delta * (share * hardEnergy) ^ 3 / 2

/-- The single leading coefficient combination visible at a fixed outgoing
energy share. -/
def effectiveThresholdCoefficient
    (deltaGamma deltaPlus deltaMinus share : ℝ) : ℝ :=
  deltaGamma
    - deltaPlus * share ^ 3
    - deltaMinus * (1 - share) ^ 3

/-- Leading head-on threshold residual after additive conservation and the
ultrarelativistic expansion. -/
noncomputable def leadingThresholdResidual
    (deltaGamma deltaPlus deltaMinus share hardEnergy softEnergy mass : ℝ) : ℝ :=
  effectiveThresholdCoefficient deltaGamma deltaPlus deltaMinus share
      * hardEnergy ^ 4
    + 4 * softEnergy * hardEnergy
    - mass ^ 2 / (share * (1 - share))

/-- Additive head-on momentum conservation implies the standard leading
threshold equation while retaining independent charged-lepton coefficients. -/
theorem additive_head_on_implies_leading_threshold
    {deltaGamma deltaPlus deltaMinus share hardEnergy softEnergy mass : ℝ}
    (hEnergy : hardEnergy ≠ 0)
    (hShare : share ≠ 0)
    (hOtherShare : 1 - share ≠ 0)
    (hMomentum :
      photonMomentumLeading deltaGamma hardEnergy - softEnergy =
        leptonMomentumLeading
            deltaPlus mass share hardEnergy softEnergy
          + leptonMomentumLeading
            deltaMinus mass (1 - share) hardEnergy softEnergy) :
    leadingThresholdResidual
      deltaGamma deltaPlus deltaMinus share hardEnergy softEnergy mass = 0 := by
  unfold photonMomentumLeading leptonMomentumLeading at hMomentum
  unfold leadingThresholdResidual effectiveThresholdCoefficient
  field_simp [hEnergy, hShare, hOtherShare] at hMomentum ⊢
  nlinarith

/-- With Lorentz-invariant charged leptons, the coefficient reduces to the
photon coefficient used by the published photon-only threshold model. -/
theorem photon_only_coefficient
    (deltaGamma share : ℝ) :
    effectiveThresholdCoefficient deltaGamma 0 0 share = deltaGamma := by
  simp [effectiveThresholdCoefficient]

/-- At equal energy sharing, the threshold sees the photon coefficient minus
one eighth of the sum of the two charged-lepton coefficients. -/
theorem symmetric_share_coefficient
    (deltaGamma deltaPlus deltaMinus : ℝ) :
    effectiveThresholdCoefficient
      deltaGamma deltaPlus deltaMinus (1 / 2 : ℝ) =
        deltaGamma - (deltaPlus + deltaMinus) / 8 := by
  unfold effectiveThresholdCoefficient
  ring

/-- If charge conjugation supplies one common lepton coefficient, the equal
share threshold depends on `deltaGamma - deltaLepton/4`. -/
theorem common_lepton_coefficient
    (deltaGamma deltaLepton : ℝ) :
    effectiveThresholdCoefficient
      deltaGamma deltaLepton deltaLepton (1 / 2 : ℝ) =
        deltaGamma - deltaLepton / 4 := by
  rw [symmetric_share_coefficient]
  ring

/-- The charge-odd lepton direction is invisible at equal sharing. -/
theorem charge_odd_lepton_direction_invisible
    (deltaGamma deltaLepton chargeOdd : ℝ) :
    effectiveThresholdCoefficient
      deltaGamma (deltaLepton + chargeOdd) (deltaLepton - chargeOdd)
        (1 / 2 : ℝ) =
      deltaGamma - deltaLepton / 4 := by
  rw [symmetric_share_coefficient]
  ring

/-- A second continuous coefficient direction also leaves the equal-share
threshold unchanged.  The photon coefficient and both lepton coefficients
cannot be recovered separately from this one threshold combination. -/
theorem common_shift_degeneracy
    (deltaGamma deltaPlus deltaMinus shift : ℝ) :
    effectiveThresholdCoefficient
      (deltaGamma + shift) (deltaPlus + 4 * shift)
        (deltaMinus + 4 * shift) (1 / 2 : ℝ) =
      effectiveThresholdCoefficient
        deltaGamma deltaPlus deltaMinus (1 / 2 : ℝ) := by
  repeat rw [symmetric_share_coefficient]
  ring

/-- Exchanging electron and positron coefficients while exchanging their
energy shares leaves the leading threshold relation unchanged. -/
theorem charge_exchange_degeneracy
    (deltaGamma deltaPlus deltaMinus share : ℝ) :
    effectiveThresholdCoefficient deltaGamma deltaPlus deltaMinus share =
      effectiveThresholdCoefficient
        deltaGamma deltaMinus deltaPlus (1 - share) := by
  unfold effectiveThresholdCoefficient
  ring

/-- On the equal-share cancellation plane, the leading threshold deformation
vanishes even when all three coefficients are nonzero. -/
theorem symmetric_share_cancellation
    {deltaGamma deltaPlus deltaMinus : ℝ}
    (hCancellation : deltaPlus + deltaMinus = 8 * deltaGamma) :
    effectiveThresholdCoefficient
      deltaGamma deltaPlus deltaMinus (1 / 2 : ℝ) = 0 := by
  rw [symmetric_share_coefficient, hCancellation]
  ring

/-- The frozen FZ-12 photon coefficient with standard charged leptons. -/
theorem fz12_photon_only_coefficient (a : ℝ) :
    effectiveThresholdCoefficient
      (-(a ^ 2) / 20) 0 0 (1 / 2 : ℝ) = -(a ^ 2) / 20 := by
  exact photon_only_coefficient _ _

/-- If the same leading principal symbol is additionally imposed on the
photon, electron, and positron, the visible coefficient is `-3 a^2/80`.
Universality is an explicit premise here, not a theorem of A1 through A3. -/
theorem fz12_universal_principal_symbol_coefficient (a : ℝ) :
    effectiveThresholdCoefficient
      (-(a ^ 2) / 20) (-(a ^ 2) / 20) (-(a ^ 2) / 20)
        (1 / 2 : ℝ) = -(3 * a ^ 2) / 80 := by
  rw [common_lepton_coefficient]
  ring

/-! ## Universal negative-coefficient share envelope

The results in this section still concern only the declared leading,
additive, collinear kinematic residual.  They do not prove a physical
pair-production channel, an interaction rate, an opacity, or a continuum
calibration.  In particular, assigning one common coefficient to the photon
and both charged particles remains an explicit universality premise.
-/

/-- The cubic share combination for a common coefficient depends only on the
product of the two outgoing energy shares. -/
theorem cubic_share_identity (share : ℝ) :
    1 - share ^ 3 - (1 - share) ^ 3 =
      3 * share * (1 - share) := by
  ring

/-- Under the explicit common negative-coefficient premise `delta = -d`, the
whole dispersion coefficient is `-3 d x(1-x)` at every share, not merely at
equal sharing. -/
theorem common_negative_coefficient (d share : ℝ) :
    effectiveThresholdCoefficient (-d) (-d) (-d) share =
      -3 * d * (share * (1 - share)) := by
  unfold effectiveThresholdCoefficient
  ring

/-- The share product used by the universal-coefficient threshold problem. -/
def shareProduct (share : ℝ) : ℝ :=
  share * (1 - share)

/-- A physical open share has a positive share product. -/
theorem shareProduct_pos
    {share : ℝ} (hSharePos : 0 < share) (hShareLt : share < 1) :
    0 < shareProduct share := by
  exact mul_pos hSharePos (sub_pos.mpr hShareLt)

/-- The share product never exceeds one quarter.  The interval hypotheses
needed for physical shares enter separately through `shareProduct_pos`. -/
theorem shareProduct_le_quarter (share : ℝ) :
    shareProduct share ≤ (1 / 4 : ℝ) := by
  unfold shareProduct
  nlinarith [sq_nonneg (2 * share - 1)]

/-- Combined physical-domain bounds for the share product. -/
theorem shareProduct_physical_bounds
    {share : ℝ} (hSharePos : 0 < share) (hShareLt : share < 1) :
    0 < shareProduct share ∧ shareProduct share ≤ (1 / 4 : ℝ) := by
  exact ⟨shareProduct_pos hSharePos hShareLt,
    shareProduct_le_quarter share⟩

/-- One quarter is attained by the share product exactly at equal sharing. -/
theorem shareProduct_eq_quarter_iff (share : ℝ) :
    shareProduct share = (1 / 4 : ℝ) ↔
      share = (1 / 2 : ℝ) := by
  unfold shareProduct
  constructor
  · intro h
    nlinarith [sq_nonneg (2 * share - 1)]
  · intro h
    subst share
    norm_num

/-- For a positive share product `u`, this is the soft energy required by the
common negative-coefficient leading residual.  It has the form `A/u + B*u`.
This is a kinematic definition, not an interaction threshold by itself. -/
noncomputable def universalSoftEnergyAtProduct
    (d mass hardEnergy product : ℝ) : ℝ :=
  (mass ^ 2 / (4 * hardEnergy)) / product
    + (3 * d * hardEnergy ^ 3 / 4) * product

/-- The share-independent lower envelope obtained by minimizing the positive
`A/u + B*u` expression over all positive `u`.  Whether its equality point is
inside the physical range `u ≤ 1/4` is a separate question. -/
noncomputable def universalSoftEnergyEnvelope
    (d mass hardEnergy : ℝ) : ℝ :=
  2 * Real.sqrt
    ((mass ^ 2 / (4 * hardEnergy)) *
      (3 * d * hardEnergy ^ 3 / 4))

/-- Exact common-negative-coefficient residual before solving for the soft
energy. -/
theorem common_negative_residual_eq
    (d share hardEnergy softEnergy mass : ℝ) :
    leadingThresholdResidual (-d) (-d) (-d)
        share hardEnergy softEnergy mass =
      -3 * d * shareProduct share * hardEnergy ^ 4
        + 4 * softEnergy * hardEnergy
        - mass ^ 2 / shareProduct share := by
  unfold leadingThresholdResidual shareProduct
  rw [common_negative_coefficient]

/-- The common-negative residual factors by the difference between the
declared soft energy and the required `A/u + B*u` value. -/
theorem common_negative_residual_factorization
    {d share hardEnergy softEnergy mass : ℝ}
    (hEnergy : hardEnergy ≠ 0)
    (hProduct : shareProduct share ≠ 0) :
    leadingThresholdResidual (-d) (-d) (-d)
        share hardEnergy softEnergy mass =
      4 * hardEnergy *
        (softEnergy - universalSoftEnergyAtProduct
          d mass hardEnergy (shareProduct share)) := by
  rw [common_negative_residual_eq]
  unfold universalSoftEnergyAtProduct
  field_simp [hEnergy, hProduct]
  ring

/-- With nonzero hard energy and share product, vanishing of the leading
residual is equivalent to the explicit required soft energy. -/
theorem common_negative_residual_zero_iff_soft_energy
    {d share hardEnergy softEnergy mass : ℝ}
    (hEnergy : hardEnergy ≠ 0)
    (hProduct : shareProduct share ≠ 0) :
    leadingThresholdResidual (-d) (-d) (-d)
        share hardEnergy softEnergy mass = 0 ↔
      softEnergy = universalSoftEnergyAtProduct
        d mass hardEnergy (shareProduct share) := by
  rw [common_negative_residual_factorization hEnergy hProduct]
  constructor
  · intro h
    have hFactor : 4 * hardEnergy ≠ 0 :=
      mul_ne_zero (by norm_num) hEnergy
    exact sub_eq_zero.mp ((mul_eq_zero.mp h).resolve_left hFactor)
  · intro h
    rw [h]
    ring

/-- Exact square-factor identity behind the reciprocal-linear AM--GM bound.
The nonnegativity hypotheses are precisely those used to split the square
root of the product. -/
theorem reciprocal_linear_gap_eq_square
    {A B product : ℝ}
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hProduct : product ≠ 0) :
    A / product + B * product - 2 * Real.sqrt (A * B) =
      (Real.sqrt A - Real.sqrt B * product) ^ 2 / product := by
  rw [Real.sqrt_mul hA]
  have hSqrtA : Real.sqrt A ^ 2 = A := Real.sq_sqrt hA
  have hSqrtB : Real.sqrt B ^ 2 = B := Real.sq_sqrt hB
  field_simp [hProduct]
  nlinarith

/-- The AM--GM lower bound for the positive reciprocal-linear expression. -/
theorem two_sqrt_le_reciprocal_linear
    {A B product : ℝ}
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hProduct : 0 < product) :
    2 * Real.sqrt (A * B) ≤ A / product + B * product := by
  have hGap := reciprocal_linear_gap_eq_square hA hB hProduct.ne'
  have hSquare :
      0 ≤ (Real.sqrt A - Real.sqrt B * product) ^ 2 / product :=
    div_nonneg (sq_nonneg _) hProduct.le
  linarith

/-- The square-factor bound is sharp: `sqrt A / sqrt B` is an exact equality
witness when both coefficients are positive. -/
theorem reciprocal_linear_balanced_witness
    {A B : ℝ} (hA : 0 < A) (hB : 0 < B) :
    A / (Real.sqrt A / Real.sqrt B)
        + B * (Real.sqrt A / Real.sqrt B) =
      2 * Real.sqrt (A * B) := by
  have hSqrtAPos : 0 < Real.sqrt A := Real.sqrt_pos.2 hA
  have hSqrtBPos : 0 < Real.sqrt B := Real.sqrt_pos.2 hB
  have hProduct : Real.sqrt A / Real.sqrt B ≠ 0 :=
    (div_pos hSqrtAPos hSqrtBPos).ne'
  have hGap := reciprocal_linear_gap_eq_square hA.le hB.le hProduct
  have hBalance :
      Real.sqrt B * (Real.sqrt A / Real.sqrt B) = Real.sqrt A := by
    field_simp [hSqrtBPos.ne']
  rw [hBalance] at hGap
  norm_num at hGap
  linarith

/-- The positive AM--GM equality witness lies in the physical product domain
whenever the reciprocal-linear coefficients obey `16 A ≤ B`. -/
theorem reciprocal_linear_balanced_witness_le_quarter
    {A B : ℝ} (hA : 0 < A) (hB : 0 < B)
    (hInterior : 16 * A ≤ B) :
    Real.sqrt A / Real.sqrt B ≤ (1 / 4 : ℝ) := by
  have hSqrtBPos : 0 < Real.sqrt B := Real.sqrt_pos.2 hB
  have hLinear : 4 * Real.sqrt A ≤ Real.sqrt B := by
    apply (sq_le_sq₀ (by positivity) (Real.sqrt_nonneg B)).mp
    rw [mul_pow, Real.sq_sqrt hA.le, Real.sq_sqrt hB.le]
    norm_num
    exact hInterior
  apply (div_le_iff₀ hSqrtBPos).2
  nlinarith

/-- Every positive product at most one quarter is realized by an open outgoing
share.  This closes the domain step that a bare AM--GM equality witness does
not supply. -/
theorem exists_open_share_with_shareProduct
    {product : ℝ} (hProduct : 0 < product)
    (hProductLe : product ≤ (1 / 4 : ℝ)) :
    ∃ share : ℝ,
      0 < share ∧ share < 1 ∧ shareProduct share = product := by
  let share : ℝ := (1 - Real.sqrt (1 - 4 * product)) / 2
  have hDisc : 0 ≤ 1 - 4 * product := by linarith
  have hDiscLt : 1 - 4 * product < 1 := by linarith
  have hSqrtLt : Real.sqrt (1 - 4 * product) < 1 := by
    simpa using Real.sqrt_lt_sqrt hDisc hDiscLt
  have hSharePos : 0 < share := by
    dsimp [share]
    linarith
  have hShareLt : share < 1 := by
    have hSqrtNonneg := Real.sqrt_nonneg (1 - 4 * product)
    dsimp [share]
    linarith
  have hSqrtSq : Real.sqrt (1 - 4 * product) ^ 2 = 1 - 4 * product :=
    Real.sq_sqrt hDisc
  have hShareProduct : shareProduct share = product := by
    unfold shareProduct
    dsimp [share]
    nlinarith
  exact ⟨share, hSharePos, hShareLt, hShareProduct⟩

/-- The AM--GM equality witness specialized to the common universal
coefficient branch. -/
noncomputable def universalBalancedProduct
    (d mass hardEnergy : ℝ) : ℝ :=
  Real.sqrt (mass ^ 2 / (4 * hardEnergy)) /
    Real.sqrt (3 * d * hardEnergy ^ 3 / 4)

theorem universalBalancedProduct_pos
    {d mass hardEnergy : ℝ}
    (hD : 0 < d) (hMass : 0 < mass) (hEnergy : 0 < hardEnergy) :
    0 < universalBalancedProduct d mass hardEnergy := by
  unfold universalBalancedProduct
  apply div_pos <;> apply Real.sqrt_pos.2 <;> positivity

/-- The reverse endpoint condition is the exact above-transition criterion
for the physical AM--GM witness. -/
theorem universal_interior_transition_iff
    {d mass hardEnergy : ℝ} (hEnergy : 0 < hardEnergy) :
    16 * (mass ^ 2 / (4 * hardEnergy)) ≤
        3 * d * hardEnergy ^ 3 / 4 ↔
      16 * mass ^ 2 ≤ 3 * d * hardEnergy ^ 4 := by
  have hLeft :
      16 * (mass ^ 2 / (4 * hardEnergy)) =
        4 * mass ^ 2 / hardEnergy := by
    field_simp [hEnergy.ne']
    ring
  rw [hLeft, div_le_iff₀ hEnergy]
  constructor <;> intro h <;> nlinarith

theorem universalBalancedProduct_le_quarter_of_supertransition
    {d mass hardEnergy : ℝ}
    (hD : 0 < d) (hMass : 0 < mass) (hEnergy : 0 < hardEnergy)
    (hTransition : 16 * mass ^ 2 ≤ 3 * d * hardEnergy ^ 4) :
    universalBalancedProduct d mass hardEnergy ≤ (1 / 4 : ℝ) := by
  unfold universalBalancedProduct
  apply reciprocal_linear_balanced_witness_le_quarter
  · positivity
  · positivity
  · exact (universal_interior_transition_iff hEnergy).2 hTransition

/-- At the specialized balance product, the universal soft-energy expression
attains the reciprocal-linear lower envelope. -/
theorem universalSoftEnergyAtProduct_balanced
    {d mass hardEnergy : ℝ}
    (hD : 0 < d) (hMass : 0 < mass) (hEnergy : 0 < hardEnergy) :
    universalSoftEnergyAtProduct d mass hardEnergy
        (universalBalancedProduct d mass hardEnergy) =
      universalSoftEnergyEnvelope d mass hardEnergy := by
  unfold universalSoftEnergyAtProduct universalBalancedProduct
    universalSoftEnergyEnvelope
  exact reciprocal_linear_balanced_witness (by positivity) (by positivity)

/-- Every positive-product soft-energy value lies above the universal
share-independent envelope. -/
theorem universalSoftEnergyEnvelope_le_at_product
    {d mass hardEnergy product : ℝ}
    (hD : 0 < d) (hMass : 0 < mass)
    (hEnergy : 0 < hardEnergy) (hProduct : 0 < product) :
    universalSoftEnergyEnvelope d mass hardEnergy ≤
      universalSoftEnergyAtProduct d mass hardEnergy product := by
  unfold universalSoftEnergyEnvelope universalSoftEnergyAtProduct
  apply two_sqrt_le_reciprocal_linear
  · positivity
  · positivity
  · exact hProduct

/-- Above the transition, the reciprocal-linear envelope is not merely an
unconstrained lower bound: it is attained by at least one physical open share
and is therefore the constrained minimum over `0 < x(1-x) ≤ 1/4`. -/
theorem universalSoftEnergyEnvelope_physical_minimum_above_transition
    {d mass hardEnergy : ℝ}
    (hD : 0 < d) (hMass : 0 < mass) (hEnergy : 0 < hardEnergy)
    (hTransition : 16 * mass ^ 2 ≤ 3 * d * hardEnergy ^ 4) :
    (∀ share : ℝ, 0 < share → share < 1 →
      universalSoftEnergyEnvelope d mass hardEnergy ≤
        universalSoftEnergyAtProduct d mass hardEnergy
          (shareProduct share)) ∧
    (∃ share : ℝ, 0 < share ∧ share < 1 ∧
      universalSoftEnergyAtProduct d mass hardEnergy
          (shareProduct share) =
        universalSoftEnergyEnvelope d mass hardEnergy) := by
  constructor
  · intro share hSharePos hShareLt
    exact universalSoftEnergyEnvelope_le_at_product
      hD hMass hEnergy (shareProduct_pos hSharePos hShareLt)
  · have hProductPos := universalBalancedProduct_pos hD hMass hEnergy
    have hProductLe :=
      universalBalancedProduct_le_quarter_of_supertransition
        hD hMass hEnergy hTransition
    obtain ⟨share, hSharePos, hShareLt, hShareProduct⟩ :=
      exists_open_share_with_shareProduct hProductPos hProductLe
    refine ⟨share, hSharePos, hShareLt, ?_⟩
    rw [hShareProduct]
    exact universalSoftEnergyAtProduct_balanced hD hMass hEnergy

/-- Generic endpoint comparison for the physical bound `u ≤ 1/4`.  For
nonnegative `A,B`, the condition `B ≤ 16 A` is exactly the regime in which
the unconstrained AM--GM balance lies at or beyond the equal-share endpoint. -/
theorem reciprocal_linear_quarter_le
    {A B product : ℝ}
    (hB : 0 ≤ B) (hProduct : 0 < product)
    (hProductLe : product ≤ (1 / 4 : ℝ))
    (hBoundary : B ≤ 16 * A) :
    A / (1 / 4 : ℝ) + B * (1 / 4 : ℝ) ≤
      A / product + B * product := by
  have hA : 0 ≤ A := by nlinarith [hB.trans hBoundary]
  have hLeft : 0 ≤ (1 / 4 : ℝ) - product := by linarith
  have hScaled : 16 * A * product ≤ 4 * A := by
    nlinarith [mul_nonneg hA hLeft]
  have hRight : 0 ≤ 4 * A / product - B := by
    have hReciprocal : 16 * A ≤ 4 * A / product :=
      (le_div_iff₀ hProduct).2 hScaled
    linarith
  have hFactor :
      0 ≤ ((1 / 4 : ℝ) - product) * (4 * A / product - B) :=
    mul_nonneg hLeft hRight
  have hIdentity :
      A / product + B * product
          - (A / (1 / 4 : ℝ) + B * (1 / 4 : ℝ)) =
        ((1 / 4 : ℝ) - product) * (4 * A / product - B) := by
    field_simp [hProduct.ne']
    ring
  linarith

/-- In the boundary regime, equal sharing minimizes the required soft energy
over all physical shares.  The transition condition is deliberately stated
in the exact `A,B` coefficients used by the proved optimization. -/
theorem universalSoftEnergy_equal_share_le
    {d mass hardEnergy share : ℝ}
    (hD : 0 < d)
    (hEnergy : 0 < hardEnergy)
    (hSharePos : 0 < share) (hShareLt : share < 1)
    (hBoundary :
      3 * d * hardEnergy ^ 3 / 4 ≤
        16 * (mass ^ 2 / (4 * hardEnergy))) :
    universalSoftEnergyAtProduct d mass hardEnergy (1 / 4 : ℝ) ≤
      universalSoftEnergyAtProduct d mass hardEnergy (shareProduct share) := by
  unfold universalSoftEnergyAtProduct
  apply reciprocal_linear_quarter_le
  · positivity
  · exact shareProduct_pos hSharePos hShareLt
  · exact shareProduct_le_quarter share
  · exact hBoundary

/-- The coefficient-domain boundary condition is exactly the fourth-power
transition `3 d E^4 ≤ 16 m^2` when the hard energy is positive. -/
theorem universal_equal_share_transition_iff
    {d mass hardEnergy : ℝ} (hEnergy : 0 < hardEnergy) :
    3 * d * hardEnergy ^ 3 / 4 ≤
        16 * (mass ^ 2 / (4 * hardEnergy)) ↔
      3 * d * hardEnergy ^ 4 ≤ 16 * mass ^ 2 := by
  have hRight :
      16 * (mass ^ 2 / (4 * hardEnergy)) =
        4 * mass ^ 2 / hardEnergy := by
    field_simp [hEnergy.ne']
    ring
  rw [hRight, le_div_iff₀ hEnergy]
  constructor <;> intro h <;> nlinarith

/-- Equal sharing is therefore the physical-share minimizer below (or at)
the exact fourth-power transition. -/
theorem universalSoftEnergy_equal_share_le_of_subtransition
    {d mass hardEnergy share : ℝ}
    (hD : 0 < d) (hEnergy : 0 < hardEnergy)
    (hSharePos : 0 < share) (hShareLt : share < 1)
    (hTransition :
      3 * d * hardEnergy ^ 4 ≤ 16 * mass ^ 2) :
    universalSoftEnergyAtProduct d mass hardEnergy (1 / 4 : ℝ) ≤
      universalSoftEnergyAtProduct d mass hardEnergy (shareProduct share) := by
  apply universalSoftEnergy_equal_share_le
      hD hEnergy hSharePos hShareLt
  exact (universal_equal_share_transition_iff hEnergy).2 hTransition

/-- The exact equal-share value in the common negative-coefficient branch. -/
theorem universalSoftEnergyAtProduct_quarter
    {d mass hardEnergy : ℝ} (hEnergy : hardEnergy ≠ 0) :
    universalSoftEnergyAtProduct d mass hardEnergy (1 / 4 : ℝ) =
      mass ^ 2 / hardEnergy + 3 * d * hardEnergy ^ 3 / 16 := by
  unfold universalSoftEnergyAtProduct
  field_simp [hEnergy]
  ring

/-- A soft energy strictly below the proved universal envelope cannot solve
the common-negative-coefficient leading residual at any physical open share.
This excludes only this conditional leading equation; it is not a no-go for
other coefficients, non-collinear kinematics, exact finite-symbol dynamics,
or a physical interaction model. -/
theorem no_common_negative_leading_solution_below_envelope
    {d mass hardEnergy softEnergy share : ℝ}
    (hD : 0 < d) (hMass : 0 < mass)
    (hEnergy : 0 < hardEnergy)
    (hSharePos : 0 < share) (hShareLt : share < 1)
    (hBelow : softEnergy <
      universalSoftEnergyEnvelope d mass hardEnergy) :
    leadingThresholdResidual (-d) (-d) (-d)
        share hardEnergy softEnergy mass ≠ 0 := by
  intro hResidual
  have hProduct : 0 < shareProduct share :=
    shareProduct_pos hSharePos hShareLt
  have hSoft :=
    (common_negative_residual_zero_iff_soft_energy
      hEnergy.ne' hProduct.ne').mp hResidual
  have hEnvelope := universalSoftEnergyEnvelope_le_at_product
    hD hMass hEnergy hProduct
  linarith

/-- Equal sharing reduces the complete leading residual to the single visible
coefficient combination and the usual mass threshold. -/
theorem symmetric_share_threshold_residual
    (deltaGamma deltaPlus deltaMinus hardEnergy softEnergy mass : ℝ) :
    leadingThresholdResidual deltaGamma deltaPlus deltaMinus
        (1 / 2 : ℝ) hardEnergy softEnergy mass =
      (deltaGamma - (deltaPlus + deltaMinus) / 8) * hardEnergy ^ 4
        + 4 * softEnergy * hardEnergy - 4 * mass ^ 2 := by
  unfold leadingThresholdResidual
  rw [symmetric_share_coefficient]
  ring

/-! ## Unique equal sharing on the Lorentz-invariant lepton branch -/

/-- Away from equal sharing, the positive-mass term in the leading collinear
threshold is strictly larger than its value at equal sharing.  This theorem is
restricted to `0 < share < 1` and a positive lepton mass. -/
theorem li_lepton_mass_penalty_strict_of_ne_symmetric
    {mass share : ℝ}
    (hMass : 0 < mass)
    (hSharePos : 0 < share)
    (hShareLt : share < 1)
    (hNe : share ≠ (1 / 2 : ℝ)) :
    4 * mass ^ 2 < mass ^ 2 / (share * (1 - share)) := by
  have hDen : 0 < share * (1 - share) :=
    mul_pos hSharePos (sub_pos.mpr hShareLt)
  have hLinear : 2 * share - 1 ≠ 0 := by
    intro h
    apply hNe
    linarith
  have hProduct : 0 < mass ^ 2 * (2 * share - 1) ^ 2 :=
    mul_pos (sq_pos_of_pos hMass) (sq_pos_of_ne_zero hLinear)
  apply (lt_div_iff₀ hDen).2
  nlinarith

/-- Equal sharing globally minimizes the positive-mass penalty in the leading
collinear threshold on the open share interval. -/
theorem li_lepton_mass_penalty_ge_symmetric
    {mass share : ℝ}
    (hMass : 0 < mass)
    (hSharePos : 0 < share)
    (hShareLt : share < 1) :
    4 * mass ^ 2 ≤ mass ^ 2 / (share * (1 - share)) := by
  by_cases hNe : share ≠ (1 / 2 : ℝ)
  · exact (li_lepton_mass_penalty_strict_of_ne_symmetric
      hMass hSharePos hShareLt hNe).le
  · have hEq : share = (1 / 2 : ℝ) := not_ne_iff.mp hNe
    subst share
    norm_num [div_eq_mul_inv]
    nlinarith

/-- The positive-mass penalty reaches its global minimum only at equal
sharing. -/
theorem li_lepton_mass_penalty_eq_symmetric_iff
    {mass share : ℝ}
    (hMass : 0 < mass)
    (hSharePos : 0 < share)
    (hShareLt : share < 1) :
    mass ^ 2 / (share * (1 - share)) = 4 * mass ^ 2 ↔
      share = (1 / 2 : ℝ) := by
  constructor
  · intro hEq
    by_contra hNe
    have hStrict := li_lepton_mass_penalty_strict_of_ne_symmetric
      hMass hSharePos hShareLt hNe
    linarith
  · intro hEq
    subst share
    norm_num
    ring

/-- With Lorentz-invariant charged leptons, equal sharing globally maximizes
the leading head-on collinear threshold residual.  No condition on the common
hard-photon or soft-photon terms is needed because they cancel in the
comparison. -/
theorem photon_only_residual_le_symmetric_share
    {mass share : ℝ}
    (hMass : 0 < mass)
    (hSharePos : 0 < share)
    (hShareLt : share < 1)
    (deltaGamma hardEnergy softEnergy : ℝ) :
    leadingThresholdResidual deltaGamma 0 0 share
        hardEnergy softEnergy mass ≤
      leadingThresholdResidual deltaGamma 0 0 (1 / 2 : ℝ)
        hardEnergy softEnergy mass := by
  have hPenalty := li_lepton_mass_penalty_ge_symmetric
    hMass hSharePos hShareLt
  unfold leadingThresholdResidual effectiveThresholdCoefficient
  norm_num
  linarith

/-- The leading head-on collinear residual on the Lorentz-invariant lepton
branch reaches its global maximum only at equal sharing. -/
theorem photon_only_residual_eq_symmetric_share_iff
    {mass share : ℝ}
    (hMass : 0 < mass)
    (hSharePos : 0 < share)
    (hShareLt : share < 1)
    (deltaGamma hardEnergy softEnergy : ℝ) :
    leadingThresholdResidual deltaGamma 0 0 share
        hardEnergy softEnergy mass =
      leadingThresholdResidual deltaGamma 0 0 (1 / 2 : ℝ)
        hardEnergy softEnergy mass ↔
      share = (1 / 2 : ℝ) := by
  have hPenalty := li_lepton_mass_penalty_eq_symmetric_iff
    hMass hSharePos hShareLt
  unfold leadingThresholdResidual effectiveThresholdCoefficient
  norm_num
  constructor
  · intro hResidual
    apply hPenalty.mp
    linarith
  · intro hShare
    subst share
    norm_num

/-! ## Exact double-root boundary of the leading symmetric threshold -/

/-- Equal-share leading threshold polynomial written in terms of its one
visible effective coefficient. -/
def symmetricThresholdPolynomial
    (deltaEffective hardEnergy softEnergy mass : ℝ) : ℝ :=
  deltaEffective * hardEnergy ^ 4
    + 4 * softEnergy * hardEnergy - 4 * mass ^ 2

/-- Formal derivative of the leading symmetric threshold polynomial with
respect to the hard energy. -/
def symmetricThresholdDerivative
    (deltaEffective hardEnergy softEnergy : ℝ) : ℝ :=
  4 * deltaEffective * hardEnergy ^ 3 + 4 * softEnergy

/-- The hard energy at which the negative-coefficient threshold polynomial
develops its double root. -/
noncomputable def criticalHardEnergy (softEnergy mass : ℝ) : ℝ :=
  4 * mass ^ 2 / (3 * softEnergy)

/-- The corresponding critical effective coefficient. -/
noncomputable def criticalEffectiveCoefficient (softEnergy mass : ℝ) : ℝ :=
  -(27 * softEnergy ^ 4) / (64 * mass ^ 6)

theorem critical_threshold_polynomial_zero
    {softEnergy mass : ℝ}
    (hSoft : softEnergy ≠ 0) :
    symmetricThresholdPolynomial
      (criticalEffectiveCoefficient softEnergy mass)
      (criticalHardEnergy softEnergy mass) softEnergy mass = 0 := by
  unfold symmetricThresholdPolynomial criticalEffectiveCoefficient
    criticalHardEnergy
  field_simp [hSoft]
  ring

theorem critical_threshold_derivative_zero
    {softEnergy mass : ℝ}
    (hMass : mass ≠ 0) :
    symmetricThresholdDerivative
      (criticalEffectiveCoefficient softEnergy mass)
      (criticalHardEnergy softEnergy mass) softEnergy = 0 := by
  unfold symmetricThresholdDerivative criticalEffectiveCoefficient
    criticalHardEnergy
  field_simp [hMass]
  ring

/-! ## Exact consequences of a subluminal photon symbol with LI leptons -/

/-- The repository's explicit Euclidean `L^2` realization of a momentum
coordinate.  The ordinary product norm on `Fin 3 -> R` is not used for
physical momentum magnitudes. -/
noncomputable def euclideanMomentum
    (momentum : ThresholdVec3) : ThresholdEuclideanVec3 :=
  OPH.PrimitivePortFrameQuotient.euclideanEquivVec3.symm momentum

/-- Euclidean momentum magnitude in the metric selected by the exact carrier
Gram completion. -/
noncomputable def euclideanMomentumMagnitude
    (momentum : ThresholdVec3) : ℝ :=
  ‖euclideanMomentum momentum‖

theorem euclideanMomentum_add (p q : ThresholdVec3) :
    euclideanMomentum (p + q) = euclideanMomentum p + euclideanMomentum q := by
  ext d
  rfl

theorem euclideanMomentumMagnitude_nonnegative (p : ThresholdVec3) :
    0 ≤ euclideanMomentumMagnitude p := by
  exact norm_nonneg _

theorem euclideanMomentumMagnitude_add_le
    (p q : ThresholdVec3) :
    euclideanMomentumMagnitude (p + q) ≤
      euclideanMomentumMagnitude p + euclideanMomentumMagnitude q := by
  unfold euclideanMomentumMagnitude
  rw [euclideanMomentum_add]
  exact norm_add_le _ _

theorem radiusSquared_eq_euclideanMomentumMagnitude_sq
    (momentum : ThresholdVec3) :
    OPH.SeamCurrentEdge30Moment.radiusSquared momentum =
      euclideanMomentumMagnitude momentum ^ 2 := by
  unfold OPH.SeamCurrentEdge30Moment.radiusSquared
    OPH.PrimitivePortTranslationBridge.dot euclideanMomentumMagnitude
    euclideanMomentum
  rw [EuclideanSpace.real_norm_sq_eq]
  apply Finset.sum_congr rfl
  intro d _
  rw [OPH.PrimitivePortFrameQuotient.euclideanEquivVec3_symm_apply]
  ring

/-- The normalized thirty-seam support has the exact Euclidean second moment
used by the finite symbol. -/
theorem unit_seam_second_moment_eq (momentum : ThresholdVec3) :
    (∑ e : Fin 30,
      OPH.PrimitivePortTranslationBridge.dot momentum
        (OPH.SeamCurrentHomogeneousAction.unitCarrierSeamDirection e) ^ 2) =
      10 * OPH.SeamCurrentEdge30Moment.radiusSquared momentum := by
  rw [← OPH.SeamCurrentEdge30Moment.seamMoment2_eq momentum]
  unfold OPH.SeamCurrentEdge30Moment.seamMoment2
  apply Finset.sum_congr rfl
  intro e _
  unfold OPH.SeamCurrentHomogeneousAction.unitCarrierSeamDirection
  rw [OPH.PrimitivePortTranslationBridge.dot_smul_right]
  rw [OPH.SeamCurrentFreePhotonLift.dot_comm momentum]
  ring

/-- Global elementary cosine bound used to compare the exact finite symbol
with the Euclidean quadratic form. -/
theorem one_sub_cos_le_half_sq (x : ℝ) :
    1 - Real.cos x ≤ x ^ 2 / 2 := by
  linarith [Real.one_sub_sq_div_two_le_cos (x := x)]

/-- The exact FZ-12 edge-current symbol is bounded above by the Euclidean
quadratic form for every nonzero coordinate scale.  This is a statement about
the complete cosine symbol, not an EFT truncation. -/
theorem exact_fz12_symbol_le_radiusSquared
    {a : ℝ} (ha : a ≠ 0) (momentum : ThresholdVec3) :
    OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol
        a momentum ≤
      OPH.SeamCurrentEdge30Moment.radiusSquared momentum := by
  rw [OPH.SeamCurrentDirichletGenerator.dilatedCompletionFourierSymbol_eq_edgeCurrentCharacterSymbol
    a ha]
  unfold OPH.SeamCurrentHomogeneousAction.edgeCurrentCharacterSymbol
  have hsum :
      (∑ e : Fin 30,
        (1 - Real.cos (a *
          OPH.PrimitivePortTranslationBridge.dot momentum
            (OPH.SeamCurrentHomogeneousAction.unitCarrierSeamDirection e)))) ≤
        ∑ e : Fin 30,
          (a * OPH.PrimitivePortTranslationBridge.dot momentum
            (OPH.SeamCurrentHomogeneousAction.unitCarrierSeamDirection e)) ^ 2 /
              2 := by
    apply Finset.sum_le_sum
    intro e _
    exact one_sub_cos_le_half_sq _
  have hfactor : 0 ≤ 1 / (5 * a ^ 2) := by positivity
  calc
    (1 / (5 * a ^ 2)) *
        ∑ e : Fin 30,
          (1 - Real.cos (a *
            OPH.PrimitivePortTranslationBridge.dot momentum
              (OPH.SeamCurrentHomogeneousAction.unitCarrierSeamDirection e))) ≤
      (1 / (5 * a ^ 2)) *
        ∑ e : Fin 30,
          (a * OPH.PrimitivePortTranslationBridge.dot momentum
            (OPH.SeamCurrentHomogeneousAction.unitCarrierSeamDirection e)) ^ 2 /
              2 := mul_le_mul_of_nonneg_left hsum hfactor
    _ = (1 / (5 * a ^ 2)) * (a ^ 2 / 2) *
        ∑ e : Fin 30,
          OPH.PrimitivePortTranslationBridge.dot momentum
            (OPH.SeamCurrentHomogeneousAction.unitCarrierSeamDirection e) ^ 2 := by
      rw [Finset.mul_sum]
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro e _
      ring
    _ = OPH.SeamCurrentEdge30Moment.radiusSquared momentum := by
      rw [unit_seam_second_moment_eq]
      field_simp [ha]
      ring

/-- The nonnegative frequency of the declared FZ-12 oscillator is bounded by
the Euclidean momentum magnitude.  Physical interpretation of either side
requires the named field, scale, space, and clock attachments. -/
theorem exact_fz12_frequency_le_euclideanMomentumMagnitude
    {a : ℝ} (ha : a ≠ 0) (momentum : ThresholdVec3) :
    OPH.SeamCurrentFreePhotonLift.photonModeFrequency a momentum ≤
      euclideanMomentumMagnitude momentum := by
  apply le_of_sq_le_sq _
    (euclideanMomentumMagnitude_nonnegative momentum)
  rw [OPH.SeamCurrentFreePhotonLift.photonModeFrequency_sq]
  rw [← radiusSquared_eq_euclideanMomentumMagnitude_sq]
  exact exact_fz12_symbol_le_radiusSquared ha momentum

/-- Standard positive-energy relativistic charged-lepton branch. -/
noncomputable def liLeptonEnergy
    (mass : ℝ) (momentum : ThresholdVec3) : ℝ :=
  Real.sqrt (mass ^ 2 + euclideanMomentumMagnitude momentum ^ 2)

/-- A positive lepton mass makes its energy strictly larger than its momentum
norm. -/
theorem euclideanMomentumMagnitude_lt_liLeptonEnergy
    {mass : ℝ} (hMass : 0 < mass) (momentum : ThresholdVec3) :
    euclideanMomentumMagnitude momentum < liLeptonEnergy mass momentum := by
  have hsq : euclideanMomentumMagnitude momentum ^ 2 <
      mass ^ 2 + euclideanMomentumMagnitude momentum ^ 2 := by
    nlinarith [sq_pos_of_pos hMass]
  simpa only [liLeptonEnergy,
    Real.sqrt_sq (euclideanMomentumMagnitude_nonnegative momentum)] using
    Real.sqrt_lt_sqrt (sq_nonneg (euclideanMomentumMagnitude momentum)) hsq

/-- The energy of any LI lepton pair is strictly larger than the norm of its
total momentum. -/
theorem total_euclidean_momentum_lt_li_pair_energy
    {mass : ℝ} (hMass : 0 < mass)
    (pPlus pMinus : ThresholdVec3) :
    euclideanMomentumMagnitude (pPlus + pMinus) <
      liLeptonEnergy mass pPlus + liLeptonEnergy mass pMinus := by
  calc
    euclideanMomentumMagnitude (pPlus + pMinus) ≤
        euclideanMomentumMagnitude pPlus +
          euclideanMomentumMagnitude pMinus :=
      euclideanMomentumMagnitude_add_le pPlus pMinus
    _ < liLeptonEnergy mass pPlus + liLeptonEnergy mass pMinus :=
      add_lt_add
        (euclideanMomentumMagnitude_lt_liLeptonEnergy hMass pPlus)
        (euclideanMomentumMagnitude_lt_liLeptonEnergy hMass pMinus)

/-- A nonnegative frequency whose square is bounded by the squared momentum
norm is subluminal in the only sense used below. -/
theorem frequency_le_momentum_norm_of_sq_le
    {frequency : ℝ} {momentum : ThresholdVec3}
    (hFrequencySq : frequency ^ 2 ≤
      euclideanMomentumMagnitude momentum ^ 2) :
    frequency ≤ euclideanMomentumMagnitude momentum :=
  le_of_sq_le_sq hFrequencySq
    (euclideanMomentumMagnitude_nonnegative momentum)

/-- Under the conditional photon attachment and standard LI leptons, an exact
subluminal photon cannot decay into an electron-positron pair while conserving
ordinary momentum and energy.  The theorem is independent of an EFT
truncation and of the rank-six direction. -/
theorem subluminal_photon_decay_impossible_with_li_leptons
    {mass frequency : ℝ} {photonMomentum pPlus pMinus : ThresholdVec3}
    (hMass : 0 < mass)
    (hMomentum : pPlus + pMinus = photonMomentum)
    (hSubluminal : frequency ≤
      euclideanMomentumMagnitude photonMomentum) :
    frequency <
      liLeptonEnergy mass pPlus + liLeptonEnergy mass pMinus := by
  calc
    frequency ≤ euclideanMomentumMagnitude photonMomentum := hSubluminal
    _ = euclideanMomentumMagnitude (pPlus + pMinus) := by rw [hMomentum]
    _ < liLeptonEnergy mass pPlus + liLeptonEnergy mass pMinus :=
      total_euclidean_momentum_lt_li_pair_energy hMass pPlus pMinus

/-- At fixed incoming momenta, the Lorentz-invariant incoming photon energy is
no smaller than the sum of two subluminal frequencies.  Consequently, every
final-state energy inside the FZ energy budget is also inside the LI energy
budget.  This order statement does not prove that a final state exists or is
dynamically reachable. -/
theorem subluminal_pair_open_implies_li_pair_open
    {omegaOne omegaTwo pairEnergy : ℝ}
    {momentumOne momentumTwo : ThresholdVec3}
    (hOne : omegaOne ≤ euclideanMomentumMagnitude momentumOne)
    (hTwo : omegaTwo ≤ euclideanMomentumMagnitude momentumTwo)
    (hOpen : pairEnergy ≤ omegaOne + omegaTwo) :
    pairEnergy ≤ euclideanMomentumMagnitude momentumOne +
      euclideanMomentumMagnitude momentumTwo := by
  exact hOpen.trans (add_le_add hOne hTwo)

end OPH.SeamCurrentPhotonLeptonThreshold

#print axioms OPH.SeamCurrentPhotonLeptonThreshold.additive_head_on_implies_leading_threshold
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.symmetric_share_coefficient
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.charge_odd_lepton_direction_invisible
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.common_shift_degeneracy
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.charge_exchange_degeneracy
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.fz12_universal_principal_symbol_coefficient
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.common_negative_residual_zero_iff_soft_energy
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.reciprocal_linear_gap_eq_square
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.reciprocal_linear_balanced_witness
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.exists_open_share_with_shareProduct
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.universalSoftEnergyEnvelope_physical_minimum_above_transition
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.universalSoftEnergy_equal_share_le
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.universal_equal_share_transition_iff
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.universalSoftEnergy_equal_share_le_of_subtransition
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.no_common_negative_leading_solution_below_envelope
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.li_lepton_mass_penalty_eq_symmetric_iff
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.photon_only_residual_le_symmetric_share
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.photon_only_residual_eq_symmetric_share_iff
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.critical_threshold_polynomial_zero
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.critical_threshold_derivative_zero
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.exact_fz12_symbol_le_radiusSquared
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.exact_fz12_frequency_le_euclideanMomentumMagnitude
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.subluminal_photon_decay_impossible_with_li_leptons
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.subluminal_pair_open_implies_li_pair_open
