import NeutralPairCoupledAction

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.NeutralPairJointStationaryWitness

open OPH.SeamCurrentCarrierQuotient
open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.PortChargeMinimalCoupling
open OPH.WorldlineHopTransport
open OPH.TransportedChargeForceLaw
open OPH.SeamChargeContinuity
open OPH.NeutralPairCoupledAction
open OPH.C1Lorentz (lorentzQ)

/-!
# A jointly stationary neutral-pair witness

The Coulomb-started packet in `NeutralPairCoupledAction` is stationary in
the field variables but fails its displayed path exchanges.  That failure
is not a no-go for the same finite action.  This module gives a different
exact field history for the same opposite crossings of seams `0` and `29`.

On the two-step window `N = 1`, the identical neutral-pair action is
stationary in both field variables and under every closed two-step path
variation of either carrier at `j = 0`.  The total load is neutral, the
sources obey continuity, and both crossings are timelike at the separately
declared carrier clock unit `3`.

This is exact finite mathematics on the committed carrier.  The action,
charge, step, clock unit, variation class, and rational field history are
declared.  The theorem neither selects this stationary point nor constructs
a laboratory realization, calibration, continuum limit, or readout.
-/

noncomputable section

/-! ## Exact field history -/

/-- Common denominator of the seam-potential table. -/
def potentialDenominator : ℤ := 9192

/-- Integer numerators of `A 0`, `A 1`, and `A 2` over the common
denominator `9192`. -/
def jointPotentialZ : Fin 3 → Fin 30 → ℤ := ![
  ![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
  ![24006, -6171, -6171, -3534, -3534, 8643, 8643, 3360, 3360, 1265,
    -1323, 2530, 1265, -1323, 2530, 0, 1265, -3534, 0, -1323, 3360,
    1265, -3534, -1323, 3360, -6171, 8643, -6171, 8643, -24006],
  ![24006, -2294, -2294, -7411, -7411, 6591, 6591, 7710, 7710, 2873,
    -2342, 3766, 2873, -2342, 3766, 0, 2873, -7411, 0, -2342, 7710,
    2873, -7411, -2342, 7710, -2294, 6591, -2294, 6591, -24006]
]

/-- Real seam-potential history.  Only slices `0,1,2` enter the window. -/
def jointPotential (n : ℕ) (e : Fin 30) : ℝ :=
  if hn : n < 3 then (jointPotentialZ ⟨n, hn⟩ e : ℝ) / potentialDenominator else 0

/-- Temporal gauge. -/
def jointScalarPotential : ℕ → Fin 12 → ℝ := fun _ _ ↦ 0

/-- Integer pair loads at the two Gauss slices. -/
def jointLoadZ : Fin 2 → Fin 12 → ℤ := ![
  ![1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
  ![0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]
]

/-- `9192` times the pair current at step zero. -/
def jointCurrentZ : Fin 30 → ℤ := fun e ↦
  if e = 0 then -18384 else if e = 29 then 18384 else 0

set_option maxRecDepth 32768 in
/-- Exact Gauss identities for the two electric-field slices. -/
theorem joint_gauss_table :
    ∀ (n : Fin 2) (p : Fin 12),
      -2 * (∑ e : Fin 30, incidenceZ e p *
        (jointPotentialZ ⟨n + 1, by omega⟩ e - jointPotentialZ ⟨n, by omega⟩ e)) =
          potentialDenominator * jointLoadZ n p := by
  decide

set_option maxRecDepth 32768 in
set_option maxHeartbeats 2000000 in
/-- Exact sourced Ampere identity after clearing the denominator. -/
theorem joint_ampere_table :
    ∀ e : Fin 30,
      4 * jointPotentialZ 2 e - 8 * jointPotentialZ 1 e +
        4 * jointPotentialZ 0 e +
        (∑ d : Fin 30, localKineticZ e d * jointPotentialZ 1 d) -
        jointCurrentZ e = 0 := by
  decide

theorem jointPotential_of_lt (n : ℕ) (hn : n < 3) (e : Fin 30) :
    jointPotential n e = (jointPotentialZ ⟨n, hn⟩ e : ℝ) / 9192 := by
  simp [jointPotential, hn, potentialDenominator]

/-! ## Source casts -/

theorem joint_pair_load (n : ℕ) (hn : n < 2) (p : Fin 12) :
    neutralPairLoad 1 (crossingWorldline 0) (crossingWorldline 29) n p =
      (jointLoadZ ⟨n, hn⟩ p : ℝ) := by
  interval_cases n
  · rw [committedNeutralPairLoad_eq, pairLoad_zero_cast]
    fin_cases p <;> rfl
  · rw [committedNeutralPairLoad_eq]
    unfold pairLoad
    rw [familyLoad_apply, Fin.sum_univ_two]
    have hPos : (pairPath 0).port 1 = 1 := by
      change (crossingPath 0).port 1 = 1
      exact crossingPath_one 0
    have hNeg : (pairPath 1).port 1 = 11 := by
      change (crossingPath 29).port 1 = 11
      exact crossingPath_one 29
    rw [hPos, hNeg]
    fin_cases p <;> norm_num [pairCharge, jointLoadZ] <;> decide

theorem joint_pair_current (e : Fin 30) :
    neutralPairCurrent 1 (1 / 2) (crossingWorldline 0) (crossingWorldline 29) 0 e =
      (jointCurrentZ e : ℝ) / 9192 := by
  rw [committedNeutralPairCurrent_eq]
  unfold pairCurrent familyCurrent
  rw [Fin.sum_univ_two]
  simp only [Pi.add_apply]
  rw [show pairPath 0 = crossingPath 0 by rfl,
    show pairPath 1 = crossingPath 29 by rfl,
    show pairCharge 0 = 1 by rfl, show pairCharge 1 = -1 by rfl,
    hoppingCurrent_forward 1 (1 / 2) (crossingPath 0) 0 0
      (crossingPath_zero 0) (crossingPath_one 0),
    hoppingCurrent_forward (-1) (1 / 2) (crossingPath 29) 0 29
      (crossingPath_zero 29) (crossingPath_one 29)]
  by_cases h0 : e = 0
  · subst e
    (simp [jointCurrentZ, show (0 : Fin 30) ≠ 29 by decide]; norm_num)
  · by_cases h29 : e = 29
    · subst e
      (simp [jointCurrentZ, show (29 : Fin 30) ≠ 0 by decide]; norm_num)
    · simp [jointCurrentZ, h0, h29]

/-! ## Real field equations -/

theorem joint_electric_table (n : ℕ) (hn : n < 2) (e : Fin 30) :
    electricFieldScaled (1 / 2) jointPotential jointScalarPotential n e =
      (-2 : ℝ) *
        ((jointPotentialZ ⟨n + 1, by omega⟩ e : ℝ) -
          (jointPotentialZ ⟨n, by omega⟩ e : ℝ)) / 9192 := by
  unfold electricFieldScaled jointScalarPotential
  simp only [Pi.sub_apply, Pi.neg_apply, Pi.smul_apply, smul_eq_mul,
    realCoboundary_apply, sub_self, sub_zero]
  rw [jointPotential_of_lt (n + 1) (by omega), jointPotential_of_lt n (by omega)]
  norm_num
  ring

theorem joint_gauss (n : ℕ) (hn : n < 2) :
    realBoundary (electricFieldScaled (1 / 2) jointPotential jointScalarPotential n) =
      neutralPairLoad 1 (crossingWorldline 0) (crossingWorldline 29) n := by
  funext p
  rw [realBoundary_apply]
  simp_rw [joint_electric_table n hn]
  rw [joint_pair_load n hn]
  have htable := joint_gauss_table ⟨n, hn⟩ p
  have htableR :
      (-2 : ℝ) * (∑ e : Fin 30, (incidenceZ e p : ℝ) *
        ((jointPotentialZ ⟨n + 1, by omega⟩ e : ℝ) -
          (jointPotentialZ ⟨n, by omega⟩ e : ℝ))) =
        9192 * (jointLoadZ ⟨n, hn⟩ p : ℝ) := by
    exact_mod_cast htable
  calc
    (∑ e : Fin 30,
        ((if p = seamRight e then
            -2 * ((jointPotentialZ ⟨n + 1, by omega⟩ e : ℝ) -
              (jointPotentialZ ⟨n, by omega⟩ e : ℝ)) / 9192 else 0) -
          (if p = seamLeft e then
            -2 * ((jointPotentialZ ⟨n + 1, by omega⟩ e : ℝ) -
              (jointPotentialZ ⟨n, by omega⟩ e : ℝ)) / 9192 else 0))) =
      ((-2 : ℝ) / 9192) *
        (∑ e : Fin 30, (incidenceZ e p : ℝ) *
          ((jointPotentialZ ⟨n + 1, by omega⟩ e : ℝ) -
            (jointPotentialZ ⟨n, by omega⟩ e : ℝ))) := by
              rw [Finset.mul_sum]
              refine Finset.sum_congr rfl fun e _ ↦ ?_
              unfold incidenceZ
              split_ifs <;> ring
    _ = (jointLoadZ ⟨n, hn⟩ p : ℝ) := by
      linear_combination (1 / (9192 : ℝ)) * htableR

theorem joint_ampere :
    electricFieldScaled (1 / 2) jointPotential jointScalarPotential 1 -
        electricFieldScaled (1 / 2) jointPotential jointScalarPotential 0 =
      (1 / 2 : ℝ) •
        (localMaxwellOperator (jointPotential 1) -
          neutralPairCurrent 1 (1 / 2) (crossingWorldline 0)
            (crossingWorldline 29) 0) := by
  funext e
  simp only [Pi.sub_apply, Pi.smul_apply, smul_eq_mul]
  rw [joint_electric_table 1 (by norm_num), joint_electric_table 0 (by norm_num),
    localMaxwellOperator_apply, joint_pair_current]
  simp_rw [jointPotential_of_lt 1 (by norm_num)]
  have hfin0 : (⟨0, by omega⟩ : Fin 3) = 0 := Fin.ext (by norm_num)
  have hfin1 : (⟨0 + 1, by omega⟩ : Fin 3) = 1 := Fin.ext (by norm_num)
  have hfin2 : (⟨1 + 1, by omega⟩ : Fin 3) = 2 := Fin.ext (by norm_num)
  simp only [hfin0, hfin1, hfin2]
  have htable := joint_ampere_table e
  have htableR :
      4 * (jointPotentialZ 2 e : ℝ) - 8 * (jointPotentialZ 1 e : ℝ) +
        4 * (jointPotentialZ 0 e : ℝ) +
        (∑ d : Fin 30, (localKineticZ e d : ℝ) * (jointPotentialZ 1 d : ℝ)) -
        (jointCurrentZ e : ℝ) = 0 := by
    exact_mod_cast htable
  rw [show (∑ x : Fin 30, (localKineticZ e x : ℝ) *
      ((jointPotentialZ 1 x : ℝ) / 9192)) =
      (∑ x : Fin 30, (localKineticZ e x : ℝ) *
        (jointPotentialZ 1 x : ℝ)) / 9192 by
        rw [Finset.sum_div]
        exact Finset.sum_congr rfl fun x _ ↦ by ring]
  norm_num
  linear_combination (-1 / (18384 : ℝ)) * htableR

/-- Field Euler--Lagrange stationarity of the exact history in the same
neutral-pair action, for arbitrary declared clock units. -/
theorem joint_field_stationary (tauPos tauNeg : ℝ) :
    NeutralPairFieldStationary 1 (1 / 2) tauPos tauNeg 1 jointPotential
      jointScalarPotential (crossingWorldline 0) (crossingWorldline 29) := by
  unfold NeutralPairFieldStationary
  apply (neutralPair_field_equations 1 (1 / 2) tauPos tauNeg (by norm_num) 1
    jointPotential jointScalarPotential (crossingWorldline 0) (crossingWorldline 29)).mpr
  exact ⟨fun m hm ↦ by interval_cases m; exact joint_ampere,
    fun n hn ↦ joint_gauss n hn⟩

/-! ## Exhaustion of the closed two-step variations -/

/-- Exact neighbor list of port `0`, including a rest. -/
theorem hop_zero_cases (p : Fin 12) (h : Hop 0 p) :
    p = 0 ∨ p = 1 ∨ p = 2 ∨ p = 3 ∨ p = 4 ∨ p = 6 := by
  rcases h with h | ⟨e, hl, hr⟩ | ⟨e, hr, hl⟩
  · exact Or.inl h.symm
  · fin_cases e <;> simp_all [seamLeft, seamRight]
  · fin_cases e <;> simp_all [seamRight]

/-- Exact neighbor list of port `1`, including a rest, read in the incoming
argument order used below. -/
theorem hop_to_one_cases (p : Fin 12) (h : Hop p 1) :
    p = 0 ∨ p = 1 ∨ p = 2 ∨ p = 3 ∨ p = 5 ∨ p = 7 := by
  rcases h with h | ⟨e, hl, hr⟩ | ⟨e, hr, hl⟩
  · exact Or.inr (Or.inl h)
  · fin_cases e <;> simp_all [seamLeft, seamRight]
  · fin_cases e <;> simp_all [seamLeft, seamRight]

/-- Exact neighbor list of port `10`, including a rest. -/
theorem hop_ten_cases (p : Fin 12) (h : Hop 10 p) :
    p = 4 ∨ p = 6 ∨ p = 8 ∨ p = 9 ∨ p = 10 ∨ p = 11 := by
  rcases h with h | ⟨e, hl, hr⟩ | ⟨e, hr, hl⟩
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inl h.symm))))
  · fin_cases e <;> simp_all [seamLeft, seamRight]
  · fin_cases e <;> simp_all [seamLeft, seamRight]

/-- Exact neighbor list of port `11`, including a rest, in incoming order. -/
theorem hop_to_eleven_cases (p : Fin 12) (h : Hop p 11) :
    p = 5 ∨ p = 7 ∨ p = 8 ∨ p = 9 ∨ p = 10 ∨ p = 11 := by
  rcases h with h | ⟨e, hl, hr⟩ | ⟨e, hr, hl⟩
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr h))))
  · fin_cases e <;> simp_all [seamLeft, seamRight]
  · fin_cases e <;> simp_all [seamLeft, seamRight]

/-- The only ports connected to both endpoints of seam `0` are its two
endpoints and the two adjacent triangle vertices `2,3`. -/
theorem seam_zero_two_hop_mids (p : Fin 12) (h0 : Hop 0 p) (h1 : Hop p 1) :
    p = 0 ∨ p = 1 ∨ p = 2 ∨ p = 3 := by
  rcases hop_zero_cases p h0 with h | h | h | h | h | h
  · exact Or.inl h
  · exact Or.inr (Or.inl h)
  · exact Or.inr (Or.inr (Or.inl h))
  · exact Or.inr (Or.inr (Or.inr h))
  · subst p
    rcases hop_to_one_cases 4 h1 with h | h | h | h | h | h <;> omega
  · subst p
    rcases hop_to_one_cases 6 h1 with h | h | h | h | h | h <;> omega

/-- The corresponding four ports for seam `29` are `10,11,8,9`. -/
theorem seam_29_two_hop_mids (p : Fin 12) (h10 : Hop 10 p) (h11 : Hop p 11) :
    p = 10 ∨ p = 11 ∨ p = 8 ∨ p = 9 := by
  rcases hop_ten_cases p h10 with h | h | h | h | h | h
  · subst p
    rcases hop_to_eleven_cases 4 h11 with h | h | h | h | h | h <;> omega
  · subst p
    rcases hop_to_eleven_cases 6 h11 with h | h | h | h | h | h <;> omega
  · exact Or.inr (Or.inr (Or.inl h))
  · exact Or.inr (Or.inr (Or.inr h))
  · exact Or.inl h
  · exact Or.inr (Or.inl h)

/-- Every closed two-step variation of the positive carrier has one of the
four enumerated intermediate ports. -/
theorem positive_variation_mid_cases
    (v : ClosedTwoStepVariation (crossingWorldline 0) 0) :
    v.mid = 0 ∨ v.mid = 1 ∨ v.mid = 2 ∨ v.mid = 3 := by
  apply seam_zero_two_hop_mids
  · simpa [ClosedTwoStepVariation.mid] using
      hop_of_admissible ((crossingWorldline 0).port 0) v.s' v.adm1
  · have h := hop_of_admissible v.mid v.s'' v.adm2
    have hend : stepTarget v.mid v.s'' = 1 := by
      change stepTarget (stepTarget ((crossingWorldline 0).port 0) v.s') v.s'' = 1
      rw [v.closed]
      rfl
    rw [hend] at h
    exact h

/-- Every closed two-step variation of the negative carrier has one of the
four enumerated intermediate ports. -/
theorem negative_variation_mid_cases
    (v : ClosedTwoStepVariation (crossingWorldline 29) 0) :
    v.mid = 10 ∨ v.mid = 11 ∨ v.mid = 8 ∨ v.mid = 9 := by
  apply seam_29_two_hop_mids
  · simpa [ClosedTwoStepVariation.mid, crossingWorldline, seamLeft] using
      hop_of_admissible ((crossingWorldline 29).port 0) v.s' v.adm1
  · have h := hop_of_admissible v.mid v.s'' v.adm2
    have hend : stepTarget v.mid v.s'' = 11 := by
      change stepTarget (stepTarget ((crossingWorldline 29).port 0) v.s') v.s'' = 11
      rw [v.closed]
      rfl
    rw [hend] at h
    exact h

/-- On an admissible carrier step, the squared spatial step norm is zero
exactly for a rest and four otherwise, expressed through its target. -/
theorem stepNormSq_of_admissible (u : Fin 12) (s : SeamStep)
    (hs : StepAdmissible u s) :
    stepNormSq s = if stepTarget u s = u then 0 else 4 := by
  cases s with
  | rest => simp [stepNormSq, stepTarget]
  | forward e =>
      have hu : u = seamLeft e := hs
      simp [stepNormSq, stepTarget, seamNormSq, hu, (seam_table_sound e).1.ne']
  | backward e =>
      have hu : u = seamRight e := hs
      simp [stepNormSq, stepTarget, seamNormSq, hu, (seam_table_sound e).1.ne]

/-- Clock difference of any closed replacement of a crossing followed by a
rest, written only in terms of its intermediate port. -/
theorem crossing_clockDifference (e : Fin 30)
    (v : ClosedTwoStepVariation (crossingWorldline e) 0) :
    clockDifference (crossingWorldline e) v =
      4 - ((if v.mid = seamLeft e then 0 else 4) +
        (if seamRight e = v.mid then 0 else 4)) := by
  have hs1 := stepNormSq_of_admissible ((crossingWorldline e).port 0) v.s' v.adm1
  have hs2 := stepNormSq_of_admissible v.mid v.s'' v.adm2
  unfold clockDifference
  rw [show (crossingWorldline e).steps 0 = .forward e by rfl,
    show (crossingWorldline e).steps (0 + 1) = .rest by rfl,
    hs1, hs2]
  simp only [stepNormSq]
  rw [show stepTarget ((crossingWorldline e).port 0) v.s' = v.mid by rfl]
  have hend : stepTarget v.mid v.s'' = seamRight e := by
    change stepTarget (stepTarget ((crossingWorldline e).port 0) v.s') v.s'' = seamRight e
    rw [v.closed]
    rfl
  rw [hend]
  simp [seamNormSq, crossingWorldline_port, crossingPath]

/-! ## Exact route balances selected by the field table -/

/-- Positive-charge interaction difference through an intermediate port. -/
def positiveRouteDelta (p : Fin 12) : ℝ :=
  routePairing (1 / 2) jointPotential jointScalarPotential 0 0 p 1 -
    routePairing (1 / 2) jointPotential jointScalarPotential 0 0 1 1

/-- Negative-charge route-pairing difference before multiplication by
charge `-1`. -/
def negativeRouteDelta (p : Fin 12) : ℝ :=
  routePairing (1 / 2) jointPotential jointScalarPotential 0 10 p 11 -
    routePairing (1 / 2) jointPotential jointScalarPotential 0 10 11 11

/-- Kernel-checked integer identities behind all nontrivial route balances. -/
theorem route_balance_table :
    jointPotentialZ 1 0 - jointPotentialZ 2 0 = 0 ∧
    (-jointPotentialZ 1 1 + jointPotentialZ 2 5 + jointPotentialZ 1 0 =
      4 * potentialDenominator) ∧
    (-jointPotentialZ 1 2 + jointPotentialZ 2 6 + jointPotentialZ 1 0 =
      4 * potentialDenominator) ∧
    (jointPotentialZ 1 29 - jointPotentialZ 2 29 = 0) ∧
    (jointPotentialZ 1 25 - jointPotentialZ 2 26 + jointPotentialZ 1 29 =
      -4 * potentialDenominator) ∧
    (jointPotentialZ 1 27 - jointPotentialZ 2 28 + jointPotentialZ 1 29 =
      -4 * potentialDenominator) := by
  decide

theorem pairing_0_1 (n : ℕ) :
    potentialPairing (1 / 2) jointPotential jointScalarPotential n 0 1 =
      (1 / 2) * jointScalarPotential n 0 - jointPotential (n + 1) 0 := by
  simpa [seamLeft, seamRight] using
    potentialPairing_forward (1 / 2) jointPotential jointScalarPotential n 0

theorem pairing_0_2 (n : ℕ) :
    potentialPairing (1 / 2) jointPotential jointScalarPotential n 0 2 =
      (1 / 2) * jointScalarPotential n 0 - jointPotential (n + 1) 1 := by
  simpa [seamLeft, seamRight] using
    potentialPairing_forward (1 / 2) jointPotential jointScalarPotential n 1

theorem pairing_2_1 (n : ℕ) :
    potentialPairing (1 / 2) jointPotential jointScalarPotential n 2 1 =
      (1 / 2) * jointScalarPotential n 2 + jointPotential (n + 1) 5 := by
  simpa [seamLeft, seamRight] using
    potentialPairing_backward (1 / 2) jointPotential jointScalarPotential n 5

theorem pairing_0_3 (n : ℕ) :
    potentialPairing (1 / 2) jointPotential jointScalarPotential n 0 3 =
      (1 / 2) * jointScalarPotential n 0 - jointPotential (n + 1) 2 := by
  simpa [seamLeft, seamRight] using
    potentialPairing_forward (1 / 2) jointPotential jointScalarPotential n 2

theorem pairing_3_1 (n : ℕ) :
    potentialPairing (1 / 2) jointPotential jointScalarPotential n 3 1 =
      (1 / 2) * jointScalarPotential n 3 + jointPotential (n + 1) 6 := by
  simpa [seamLeft, seamRight] using
    potentialPairing_backward (1 / 2) jointPotential jointScalarPotential n 6

theorem pairing_10_11 (n : ℕ) :
    potentialPairing (1 / 2) jointPotential jointScalarPotential n 10 11 =
      (1 / 2) * jointScalarPotential n 10 - jointPotential (n + 1) 29 := by
  simpa [seamLeft, seamRight] using
    potentialPairing_forward (1 / 2) jointPotential jointScalarPotential n 29

theorem pairing_10_8 (n : ℕ) :
    potentialPairing (1 / 2) jointPotential jointScalarPotential n 10 8 =
      (1 / 2) * jointScalarPotential n 10 + jointPotential (n + 1) 25 := by
  simpa [seamLeft, seamRight] using
    potentialPairing_backward (1 / 2) jointPotential jointScalarPotential n 25

theorem pairing_8_11 (n : ℕ) :
    potentialPairing (1 / 2) jointPotential jointScalarPotential n 8 11 =
      (1 / 2) * jointScalarPotential n 8 - jointPotential (n + 1) 26 := by
  simpa [seamLeft, seamRight] using
    potentialPairing_forward (1 / 2) jointPotential jointScalarPotential n 26

theorem pairing_10_9 (n : ℕ) :
    potentialPairing (1 / 2) jointPotential jointScalarPotential n 10 9 =
      (1 / 2) * jointScalarPotential n 10 + jointPotential (n + 1) 27 := by
  simpa [seamLeft, seamRight] using
    potentialPairing_backward (1 / 2) jointPotential jointScalarPotential n 27

theorem pairing_9_11 (n : ℕ) :
    potentialPairing (1 / 2) jointPotential jointScalarPotential n 9 11 =
      (1 / 2) * jointScalarPotential n 9 - jointPotential (n + 1) 28 := by
  simpa [seamLeft, seamRight] using
    potentialPairing_forward (1 / 2) jointPotential jointScalarPotential n 28

theorem positiveRouteDelta_zero : positiveRouteDelta 0 = 0 := by
  unfold positiveRouteDelta routePairing
  rw [potentialPairing_rest, pairing_0_1, pairing_0_1, potentialPairing_rest]
  norm_num [jointPotential, jointPotentialZ, potentialDenominator, jointScalarPotential]

theorem positiveRouteDelta_one : positiveRouteDelta 1 = 0 := by
  unfold positiveRouteDelta routePairing
  rw [pairing_0_1, potentialPairing_rest]
  ring

theorem positiveRouteDelta_two : positiveRouteDelta 2 = 4 := by
  unfold positiveRouteDelta routePairing
  rw [pairing_0_2, pairing_2_1, pairing_0_1, potentialPairing_rest]
  simp only [jointScalarPotential, mul_zero, zero_add]
  simp_rw [jointPotential_of_lt 1 (by norm_num), jointPotential_of_lt 2 (by norm_num)]
  have hfin1 : (⟨1, by omega⟩ : Fin 3) = 1 := Fin.ext rfl
  have hfin2 : (⟨2, by omega⟩ : Fin 3) = 2 := Fin.ext rfl
  simp only [hfin1, hfin2]
  have hR :
      -(jointPotentialZ 1 1 : ℝ) + (jointPotentialZ 2 5 : ℝ) +
        (jointPotentialZ 1 0 : ℝ) = 4 * 9192 := by
    exact_mod_cast route_balance_table.2.1
  linear_combination (1 / (9192 : ℝ)) * hR

theorem positiveRouteDelta_three : positiveRouteDelta 3 = 4 := by
  unfold positiveRouteDelta routePairing
  rw [pairing_0_3, pairing_3_1, pairing_0_1, potentialPairing_rest]
  simp only [jointScalarPotential, mul_zero, zero_add]
  simp_rw [jointPotential_of_lt 1 (by norm_num), jointPotential_of_lt 2 (by norm_num)]
  have hfin1 : (⟨1, by omega⟩ : Fin 3) = 1 := Fin.ext rfl
  have hfin2 : (⟨2, by omega⟩ : Fin 3) = 2 := Fin.ext rfl
  simp only [hfin1, hfin2]
  have hR :
      -(jointPotentialZ 1 2 : ℝ) + (jointPotentialZ 2 6 : ℝ) +
        (jointPotentialZ 1 0 : ℝ) = 4 * 9192 := by
    exact_mod_cast route_balance_table.2.2.1
  linear_combination (1 / (9192 : ℝ)) * hR

theorem negativeRouteDelta_ten : negativeRouteDelta 10 = 0 := by
  unfold negativeRouteDelta routePairing
  rw [potentialPairing_rest, pairing_10_11, pairing_10_11, potentialPairing_rest]
  simp only [jointScalarPotential, mul_zero, zero_add]
  simp_rw [jointPotential_of_lt 1 (by norm_num), jointPotential_of_lt 2 (by norm_num)]
  have hfin1 : (⟨1, by omega⟩ : Fin 3) = 1 := Fin.ext rfl
  have hfin2 : (⟨2, by omega⟩ : Fin 3) = 2 := Fin.ext rfl
  simp only [hfin1, hfin2]
  have hR : (jointPotentialZ 1 29 : ℝ) - (jointPotentialZ 2 29 : ℝ) = 0 := by
    exact_mod_cast route_balance_table.2.2.2.1
  linear_combination hR / (9192 : ℝ)

theorem negativeRouteDelta_eleven : negativeRouteDelta 11 = 0 := by
  unfold negativeRouteDelta routePairing
  rw [pairing_10_11, potentialPairing_rest]
  ring

theorem negativeRouteDelta_eight : negativeRouteDelta 8 = -4 := by
  unfold negativeRouteDelta routePairing
  rw [pairing_10_8, pairing_8_11, pairing_10_11, potentialPairing_rest]
  simp only [jointScalarPotential, mul_zero, zero_add]
  simp_rw [jointPotential_of_lt 1 (by norm_num), jointPotential_of_lt 2 (by norm_num)]
  have hfin1 : (⟨1, by omega⟩ : Fin 3) = 1 := Fin.ext rfl
  have hfin2 : (⟨2, by omega⟩ : Fin 3) = 2 := Fin.ext rfl
  simp only [hfin1, hfin2]
  have hR :
      (jointPotentialZ 1 25 : ℝ) - (jointPotentialZ 2 26 : ℝ) +
        (jointPotentialZ 1 29 : ℝ) = -4 * 9192 := by
    exact_mod_cast route_balance_table.2.2.2.2.1
  linear_combination (1 / (9192 : ℝ)) * hR

theorem negativeRouteDelta_nine : negativeRouteDelta 9 = -4 := by
  unfold negativeRouteDelta routePairing
  rw [pairing_10_9, pairing_9_11, pairing_10_11, potentialPairing_rest]
  simp only [jointScalarPotential, mul_zero, zero_add]
  simp_rw [jointPotential_of_lt 1 (by norm_num), jointPotential_of_lt 2 (by norm_num)]
  have hfin1 : (⟨1, by omega⟩ : Fin 3) = 1 := Fin.ext rfl
  have hfin2 : (⟨2, by omega⟩ : Fin 3) = 2 := Fin.ext rfl
  simp only [hfin1, hfin2]
  have hR :
      (jointPotentialZ 1 27 : ℝ) - (jointPotentialZ 2 28 : ℝ) +
        (jointPotentialZ 1 29 : ℝ) = -4 * 9192 := by
    exact_mod_cast route_balance_table.2.2.2.2.2
  linear_combination (1 / (9192 : ℝ)) * hR

/-! ## All closed two-step path variations -/

/-- For the positive crossing, the interaction part of an arbitrary closed
two-step variation is exactly the route table indexed by its intermediate
port. -/
theorem positive_interaction_eq_routeDelta
    (v : ClosedTwoStepVariation (crossingWorldline 0) 0) :
    interactionDifference 1 (1 / 2) jointPotential jointScalarPotential
        (crossingWorldline 0) v = positiveRouteDelta v.mid := by
  unfold interactionDifference positiveRouteDelta
  simp only [one_mul]
  rw [crossingWorldline_port, crossingWorldline_port, crossingWorldline_port]
  simp [crossingPath, seamLeft, seamRight]

/-- For the negative crossing, charge `-1` reverses the corresponding route
table. -/
theorem negative_interaction_eq_routeDelta
    (v : ClosedTwoStepVariation (crossingWorldline 29) 0) :
    interactionDifference (-1) (1 / 2) jointPotential jointScalarPotential
        (crossingWorldline 29) v = -negativeRouteDelta v.mid := by
  unfold interactionDifference negativeRouteDelta
  rw [crossingWorldline_port, crossingWorldline_port, crossingWorldline_port]
  simp [crossingPath, seamLeft, seamRight]

/-- Every allowed closed two-step variation of the positive carrier has
zero clock-plus-interaction first difference in the declared field. -/
theorem positive_all_closed_balance
    (v : ClosedTwoStepVariation (crossingWorldline 0) 0) :
    clockDifference (crossingWorldline 0) v +
        interactionDifference 1 (1 / 2) jointPotential jointScalarPotential
          (crossingWorldline 0) v = 0 := by
  rw [crossing_clockDifference, positive_interaction_eq_routeDelta,
    seam_zero_endpoints.1, seam_zero_endpoints.2]
  rcases positive_variation_mid_cases v with hm | hm | hm | hm
  · rw [hm, positiveRouteDelta_zero]
    norm_num
  · rw [hm, positiveRouteDelta_one]
    norm_num
  · rw [hm, positiveRouteDelta_two]
    simp [show (2 : Fin 12) ≠ 0 by decide, show (1 : Fin 12) ≠ 2 by decide]
  · rw [hm, positiveRouteDelta_three]
    simp [show (3 : Fin 12) ≠ 0 by decide, show (1 : Fin 12) ≠ 3 by decide]

/-- Every allowed closed two-step variation of the negative carrier has
zero clock-plus-interaction first difference in the declared field. -/
theorem negative_all_closed_balance
    (v : ClosedTwoStepVariation (crossingWorldline 29) 0) :
    clockDifference (crossingWorldline 29) v +
        interactionDifference (-1) (1 / 2) jointPotential jointScalarPotential
          (crossingWorldline 29) v = 0 := by
  rw [crossing_clockDifference, negative_interaction_eq_routeDelta,
    seam_29_endpoints.1, seam_29_endpoints.2]
  rcases negative_variation_mid_cases v with hm | hm | hm | hm
  · rw [hm, negativeRouteDelta_ten]
    simp [show (11 : Fin 12) ≠ 10 by decide]
  · rw [hm, negativeRouteDelta_eleven]
    simp [show (11 : Fin 12) ≠ 10 by decide]
  · rw [hm, negativeRouteDelta_eight]
    simp [show (8 : Fin 12) ≠ 10 by decide, show (11 : Fin 12) ≠ 8 by decide]
  · rw [hm, negativeRouteDelta_nine]
    simp [show (9 : Fin 12) ≠ 10 by decide, show (11 : Fin 12) ≠ 9 by decide]

/-- Exact stationarity of the same finite action under an arbitrary closed
two-step variation of the positive path in the `N = 1` window. -/
theorem positive_all_closed_stationary (tauPos tauNeg : ℝ)
    (v : ClosedTwoStepVariation (crossingWorldline 0) 0) :
    neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg 1 jointPotential
          jointScalarPotential v.worldline (crossingWorldline 29) =
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg 1 jointPotential
          jointScalarPotential (crossingWorldline 0) (crossingWorldline 29) := by
  have hVariation := positivePath_closed_variation 1 (1 / 2) tauPos tauNeg
    (by norm_num) 1 jointPotential jointScalarPotential (crossingWorldline 0)
    (crossingWorldline 29) v (by norm_num)
  rw [positive_all_closed_balance v] at hVariation
  linarith

/-- Exact stationarity of the same finite action under an arbitrary closed
two-step variation of the negative path in the `N = 1` window. -/
theorem negative_all_closed_stationary (tauPos tauNeg : ℝ)
    (v : ClosedTwoStepVariation (crossingWorldline 29) 0) :
    neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg 1 jointPotential
          jointScalarPotential (crossingWorldline 0) v.worldline =
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg 1 jointPotential
          jointScalarPotential (crossingWorldline 0) (crossingWorldline 29) := by
  have hVariation := negativePath_closed_variation 1 (1 / 2) tauPos tauNeg
    (by norm_num) 1 jointPotential jointScalarPotential (crossingWorldline 0)
    (crossingWorldline 29) v (by norm_num)
  rw [negative_all_closed_balance v] at hVariation
  linarith

/-- Stationarity in every positive-path two-step slot contained in the
one-step interior window.  The index bound forces `j = 0`, but leaves the
closed replacement itself arbitrary. -/
def PositiveWindowPathStationary (tauPos tauNeg : ℝ) : Prop :=
  ∀ (j : ℕ) (v : ClosedTwoStepVariation (crossingWorldline 0) j),
    j + 1 < 1 + 1 →
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg 1 jointPotential
            jointScalarPotential v.worldline (crossingWorldline 29) =
        neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg 1 jointPotential
            jointScalarPotential (crossingWorldline 0) (crossingWorldline 29)

/-- Stationarity in every negative-path two-step slot contained in the
one-step interior window. -/
def NegativeWindowPathStationary (tauPos tauNeg : ℝ) : Prop :=
  ∀ (j : ℕ) (v : ClosedTwoStepVariation (crossingWorldline 29) j),
    j + 1 < 1 + 1 →
      neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg 1 jointPotential
            jointScalarPotential (crossingWorldline 0) v.worldline =
        neutralPairCoupledAction 1 (1 / 2) tauPos tauNeg 1 jointPotential
            jointScalarPotential (crossingWorldline 0) (crossingWorldline 29)

theorem positive_window_path_stationary (tauPos tauNeg : ℝ) :
    PositiveWindowPathStationary tauPos tauNeg := by
  intro j v hj
  have hj0 : j = 0 := by omega
  subst j
  exact positive_all_closed_stationary tauPos tauNeg v

theorem negative_window_path_stationary (tauPos tauNeg : ℝ) :
    NegativeWindowPathStationary tauPos tauNeg := by
  intro j v hj
  have hj0 : j = 0 := by omega
  subst j
  exact negative_all_closed_stationary tauPos tauNeg v

/-! ## Same-instance receipt and scope contacts -/

/-- A typed receipt that the one declared field history and crossing pair
are stationary in both field slots and in every in-window closed two-step
path slot, while their exact sources remain neutral and continuous. -/
structure JointStationaryReceipt (tauPos tauNeg : ℝ) : Prop where
  field : NeutralPairFieldStationary 1 (1 / 2) tauPos tauNeg 1 jointPotential
    jointScalarPotential (crossingWorldline 0) (crossingWorldline 29)
  positivePath : PositiveWindowPathStationary tauPos tauNeg
  negativePath : NegativeWindowPathStationary tauPos tauNeg
  neutral : ∀ n, (∑ p : Fin 12,
    neutralPairLoad 1 (crossingWorldline 0) (crossingWorldline 29) n p) = 0
  continuity : ∀ n,
    neutralPairLoad 1 (crossingWorldline 0) (crossingWorldline 29) (n + 1) -
      neutralPairLoad 1 (crossingWorldline 0) (crossingWorldline 29) n +
        (1 / 2 : ℝ) • realBoundary
          (neutralPairCurrent 1 (1 / 2) (crossingWorldline 0)
            (crossingWorldline 29) n) = 0

/-- The rational packet is jointly stationary for arbitrary declared clock
units.  This does not assert that the packet is unique or dynamically
selected. -/
theorem joint_stationary_receipt (tauPos tauNeg : ℝ) :
    JointStationaryReceipt tauPos tauNeg where
  field := joint_field_stationary tauPos tauNeg
  positivePath := positive_window_path_stationary tauPos tauNeg
  negativePath := negative_window_path_stationary tauPos tauNeg
  neutral := neutralPairLoad_total 1 (crossingWorldline 0) (crossingWorldline 29)
  continuity := neutralPair_continuity 1 (1 / 2) (by norm_num)
    (crossingWorldline 0) (crossingWorldline 29)

/-- The field equations are genuinely sourced on both crossed seams, with
opposite current and load signs. -/
theorem joint_sources_nonzero :
    neutralPairCurrent 1 (1 / 2) (crossingWorldline 0)
        (crossingWorldline 29) 0 0 = -2 ∧
      neutralPairCurrent 1 (1 / 2) (crossingWorldline 0)
        (crossingWorldline 29) 0 29 = 2 ∧
      neutralPairLoad 1 (crossingWorldline 0) (crossingWorldline 29) 0 0 = 1 ∧
      neutralPairLoad 1 (crossingWorldline 0) (crossingWorldline 29) 0 10 = -1 := by
  constructor
  · rw [joint_pair_current]
    norm_num [jointCurrentZ]
  constructor
  · rw [joint_pair_current]
    norm_num [jointCurrentZ, show (29 : Fin 30) ≠ 0 by decide]
  constructor
  · rw [joint_pair_load 0 (by norm_num)]
    norm_num [jointLoadZ]
  · rw [joint_pair_load 0 (by norm_num)]
    change ((-1 : ℤ) : ℝ) = -1
    norm_num

/-- At carrier clock unit `3`, the jointly stationary source crossings are
timelike in the existing finite Lorentz module. -/
theorem joint_crossings_timelike :
    0 < lorentzQ (generatedPath 3 (crossingWorldline 0) (0 + 1) -
      generatedPath 3 (crossingWorldline 0) 0) ∧
    0 < lorentzQ (generatedPath 3 (crossingWorldline 29) (0 + 1) -
      generatedPath 3 (crossingWorldline 29) 0) :=
  committedNeutralPair_initial_crossings_timelike

/-- **Exact same-instance joint-stationarity witness.**  One finite action,
one nonzero neutral source pair, and one rational field history satisfy the
sourced field Euler--Lagrange equations and every in-window closed two-step
path stationarity equation.  The source crossings are timelike at the
separately declared clock unit `3`. -/
theorem exact_joint_stationary_timelike_sourced :
    JointStationaryReceipt 3 3 ∧
      (0 < lorentzQ (generatedPath 3 (crossingWorldline 0) (0 + 1) -
          generatedPath 3 (crossingWorldline 0) 0) ∧
        0 < lorentzQ (generatedPath 3 (crossingWorldline 29) (0 + 1) -
          generatedPath 3 (crossingWorldline 29) 0)) ∧
      (neutralPairCurrent 1 (1 / 2) (crossingWorldline 0)
          (crossingWorldline 29) 0 0 = -2 ∧
        neutralPairCurrent 1 (1 / 2) (crossingWorldline 0)
          (crossingWorldline 29) 0 29 = 2 ∧
        neutralPairLoad 1 (crossingWorldline 0)
          (crossingWorldline 29) 0 0 = 1 ∧
        neutralPairLoad 1 (crossingWorldline 0)
          (crossingWorldline 29) 0 10 = -1) := by
  exact ⟨joint_stationary_receipt 3 3, joint_crossings_timelike,
    joint_sources_nonzero⟩

end

end OPH.NeutralPairJointStationaryWitness

#print axioms OPH.NeutralPairJointStationaryWitness.joint_gauss_table
#print axioms OPH.NeutralPairJointStationaryWitness.joint_ampere_table
#print axioms OPH.NeutralPairJointStationaryWitness.joint_field_stationary
#print axioms OPH.NeutralPairJointStationaryWitness.positive_variation_mid_cases
#print axioms OPH.NeutralPairJointStationaryWitness.negative_variation_mid_cases
#print axioms OPH.NeutralPairJointStationaryWitness.positive_all_closed_stationary
#print axioms OPH.NeutralPairJointStationaryWitness.negative_all_closed_stationary
#print axioms OPH.NeutralPairJointStationaryWitness.exact_joint_stationary_timelike_sourced
