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
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.li_lepton_mass_penalty_eq_symmetric_iff
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.photon_only_residual_le_symmetric_share
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.photon_only_residual_eq_symmetric_share_iff
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.critical_threshold_polynomial_zero
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.critical_threshold_derivative_zero
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.exact_fz12_symbol_le_radiusSquared
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.exact_fz12_frequency_le_euclideanMomentumMagnitude
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.subluminal_photon_decay_impossible_with_li_leptons
#print axioms OPH.SeamCurrentPhotonLeptonThreshold.subluminal_pair_open_implies_li_pair_open
