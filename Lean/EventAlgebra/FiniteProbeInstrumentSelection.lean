import EventAlgebra.InstrumentSelectionByCertainStates

set_option autoImplicit false

namespace EventAlgebra

/-!
# Finite pure-state probes select the identity channel

A finite Kraus family fixing every coordinate pure-state matrix and the
all-ones rank-one matrix consists of scalar identity operators, and its
channel is the identity map. These are
`n + 1` algebraic fixation conditions in dimension `n`; normalization of the
all-ones input gives the same condition by linearity.

The probability form assumes a normalized Kraus channel and exact unit
Born return probability for each normalized probe. It selects the Lüders
map of the sure event in the finite complex matrix semantics of
`InstrumentSelectionByCertainStates`. No finite-data inference, complete
family of labeled blocks, or quantitative stability bound is proved here.
-/

open Matrix
open OPH.QFT
open OPH.Dynamics (CMat)
open scoped ComplexOrder

noncomputable section

variable {n : ℕ}

/-- Conjugating a pure rank-one matrix gives the rank-one matrix of the
transformed vector. -/
theorem kraus_mul_pure_matrix (K : CMat n) (ψ : Fin n → ℂ) :
    K * vecMulVec ψ (star ψ) * Kᴴ =
      vecMulVec (K *ᵥ ψ) (star (K *ᵥ ψ)) := by
  rw [mul_vecMulVec, vecMulVec_mul, star_mulVec]

/-- A Kraus sum fixing a pure rank-one matrix sends its defining vector
into its own ray in every Kraus branch. -/
theorem kraus_mulVec_collinear_of_fixes_pure {ι : Type*} [Fintype ι]
    {K : ι → CMat n} {ψ : Fin n → ℂ}
    (hfix : ∑ c, K c * vecMulVec ψ (star ψ) * (K c)ᴴ =
      vecMulVec ψ (star ψ)) (c : ι) :
    ∃ lam : ℂ, K c *ᵥ ψ = lam • ψ := by
  simp only [kraus_mul_pure_matrix] at hfix
  exact rank_one_decomposition_collinear hfix c

/-- Coordinate pure-state fixation and one coherent pure-state fixation
force every Kraus operator to be a scalar identity. -/
theorem kraus_eq_smul_one_of_finite_pure_probes {ι : Type*} [Fintype ι]
    {K : ι → CMat n}
    (hbasis : ∀ j : Fin n,
      ∑ c, K c * vecMulVec (Pi.single j 1) (star (Pi.single j 1)) * (K c)ᴴ =
        vecMulVec (Pi.single j 1) (star (Pi.single j 1)))
    (hcoherent : ∑ c, K c * vecMulVec (fun _ => 1) (star (fun _ => 1)) * (K c)ᴴ =
      vecMulVec (fun _ => 1) (star (fun _ => 1))) (c : ι) :
    ∃ lam : ℂ, K c = lam • (1 : CMat n) := by
  obtain ⟨lam, hlam⟩ := kraus_mulVec_collinear_of_fixes_pure hcoherent c
  have hoff : ∀ i j : Fin n, i ≠ j → K c i j = 0 := by
    intro i j hij
    obtain ⟨a, ha⟩ := kraus_mulVec_collinear_of_fixes_pure (hbasis j) c
    have hi := congrFun ha i
    simpa [mulVec_single_one, Pi.single_eq_of_ne hij] using hi
  have hdiag : ∀ i : Fin n, K c i i = lam := by
    intro i
    have hi := congrFun hlam i
    simp only [mulVec, dotProduct, mul_one, Pi.smul_apply, smul_eq_mul] at hi
    rw [Finset.sum_eq_single i] at hi
    · exact hi
    · intro j _ hji
      exact hoff i j (Ne.symm hji)
    · simp
  refine ⟨lam, ?_⟩
  ext i j
  by_cases hij : i = j
  · subst j
    simp [hdiag]
  · simp [hoff i j hij, hij]

/-- The finite Kraus channel selected by the coordinate and
coherent probes acts as the identity on every matrix. -/
theorem kraus_channel_eq_id_of_finite_pure_probes {ι : Type*} [Fintype ι]
    {K : ι → CMat n}
    (hbasis : ∀ j : Fin n,
      ∑ c, K c * vecMulVec (Pi.single j 1) (star (Pi.single j 1)) * (K c)ᴴ =
        vecMulVec (Pi.single j 1) (star (Pi.single j 1)))
    (hcoherent : ∑ c, K c * vecMulVec (fun _ => 1) (star (fun _ => 1)) * (K c)ᴴ =
      vecMulVec (fun _ => 1) (star (fun _ => 1))) (X : CMat n) :
    ∑ c, K c * X * (K c)ᴴ = X := by
  choose lam hlam using kraus_eq_smul_one_of_finite_pure_probes hbasis hcoherent
  have hexp (Y : CMat n) : ∑ c, K c * Y * (K c)ᴴ =
      (∑ c, lam c * star (lam c)) • Y := by
    rw [Finset.sum_smul]
    apply Finset.sum_congr rfl
    intro c _
    rw [hlam c]
    simp only [conjTranspose_smul, conjTranspose_one, smul_mul_assoc,
      mul_smul_comm, one_mul, mul_one, smul_smul]
    rw [mul_comm]
  ext i j
  have hi := congrFun (congrFun (hbasis i) i) i
  rw [hexp] at hi
  simp only [Matrix.smul_apply, vecMulVec_apply, Pi.single_eq_same,
    Pi.star_apply, star_one, mul_one, smul_eq_mul] at hi
  rw [hexp, hi, one_smul]

/-- In the committed finite matrix semantics, a Kraus-form map fixing
the finite probes equals the Lüders map of the sure event. -/
theorem finite_pure_probes_select_lueders_one
    {Φ : CMat n →ₗ[ℂ] CMat n} (hK : HasKrausFormWithEffect Φ 1)
    (hbasis : ∀ j : Fin n,
      Φ (vecMulVec (Pi.single j 1) (star (Pi.single j 1))) =
        vecMulVec (Pi.single j 1) (star (Pi.single j 1)))
    (hcoherent : Φ (vecMulVec (fun _ => 1) (star (fun _ => 1))) =
      vecMulVec (fun _ => 1) (star (fun _ => 1))) :
    Φ = luedersOutcomeMap (1 : CMat n) := by
  obtain ⟨m, K, hΦ, _⟩ := hK
  apply LinearMap.ext
  intro X
  rw [hΦ, kraus_channel_eq_id_of_finite_pure_probes
    (fun j => by rw [← hΦ]; exact hbasis j)
    (by rw [← hΦ]; exact hcoherent)]
  simp [luedersOutcomeMap_apply]

/-- The coherent probe can be supplied as the trace-one uniform pure
state. Together with the coordinate pure states it selects the identity
channel among maps with the committed Kraus form and sure effect. -/
theorem normalized_finite_pure_probes_select_lueders_one
    (hn : n ≠ 0) {Φ : CMat n →ₗ[ℂ] CMat n}
    (hK : HasKrausFormWithEffect Φ 1)
    (hbasis : ∀ j : Fin n,
      Φ (vecMulVec (Pi.single j 1) (star (Pi.single j 1))) =
        vecMulVec (Pi.single j 1) (star (Pi.single j 1)))
    (hcoherent : Φ ((n : ℂ)⁻¹ • vecMulVec (fun _ => 1) (star (fun _ => 1))) =
      (n : ℂ)⁻¹ • vecMulVec (fun _ => 1) (star (fun _ => 1))) :
    Φ = luedersOutcomeMap (1 : CMat n) := by
  apply finite_pure_probes_select_lueders_one hK hbasis
  have hnC : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn
  apply smul_right_injective (CMat n) (inv_ne_zero hnC)
  simpa only [map_smul] using hcoherent

/-- A scalar multiple of a pure matrix compresses any matrix to its
Born weight times that same pure matrix. -/
theorem scaled_pure_sandwich (a : ℂ) (ψ : Fin n → ℂ) (X : CMat n) :
    (a • vecMulVec ψ (star ψ)) * X * (a • vecMulVec ψ (star ψ)) =
      bornWeight X (a • vecMulVec ψ (star ψ)) • (a • vecMulVec ψ (star ψ)) := by
  have hpure : vecMulVec ψ (star ψ) * X * vecMulVec ψ (star ψ) =
      bornWeight X (vecMulVec ψ (star ψ)) • vecMulVec ψ (star ψ) := by
    rw [vecMulVec_mul, vecMulVec_mul_vecMulVec, vecMulVec_smul,
      bornWeight, mul_vecMulVec, trace_vecMulVec,
      dotProduct_comm (X *ᵥ ψ) (star ψ), dotProduct_mulVec]
  have hweight : bornWeight X (a • vecMulVec ψ (star ψ)) =
      a * bornWeight X (vecMulVec ψ (star ψ)) := by
    simp only [bornWeight, mul_smul_comm, trace_smul, smul_eq_mul]
  simp only [smul_mul_assoc, mul_smul_comm, hpure, hweight, smul_smul]
  congr 1
  ring

/-- A state certain of a rank-one event is that pure state itself. -/
theorem state_eq_scaled_pure_of_bornWeight_one {σ P : CMat n}
    (hσ : IsState σ) (hP : IsEvent P) (a : ℂ) (ψ : Fin n → ℂ)
    (hform : P = a • vecMulVec ψ (star ψ)) (hreturn : bornWeight σ P = 1) :
    σ = P := by
  obtain ⟨hright, hleft⟩ := mul_eq_self_of_bornWeight_one hσ hP hreturn
  calc
    σ = P * σ * P := by rw [hleft, hright]
    _ = bornWeight σ P • P := by rw [hform, scaled_pure_sandwich]
    _ = P := by rw [hreturn, one_smul]

/-- A normalized Kraus channel with return probability one on a pure
state fixes that state. -/
theorem fixes_scaled_pure_of_return_probability_one
    {Φ : CMat n →ₗ[ℂ] CMat n} (hK : HasKrausFormWithEffect Φ 1)
    {P : CMat n} (hstate : IsState P) (hevent : IsEvent P)
    (a : ℂ) (ψ : Fin n → ℂ) (hform : P = a • vecMulVec ψ (star ψ))
    (hreturn : bornWeight (Φ P) P = 1) : Φ P = P := by
  have houtput : IsState (Φ P) :=
    ⟨hK.isCompletelyPositive.isPositiveMap _ hstate.1,
      by rw [hK.trace_apply, bornWeight_one hstate]⟩
  exact state_eq_scaled_pure_of_bornWeight_one houtput hevent a ψ hform hreturn

/-- A rank-one state is an event. -/
theorem isEvent_of_state_scaled_pure {P : CMat n} (hstate : IsState P)
    (a : ℂ) (ψ : Fin n → ℂ) (hform : P = a • vecMulVec ψ (star ψ)) :
    IsEvent P := by
  refine ⟨hstate.1.isHermitian, ?_⟩
  calc
    P * P = P * 1 * P := by rw [mul_one]
    _ = bornWeight 1 P • P := by rw [hform, scaled_pure_sandwich]
    _ = P := by rw [bornWeight, one_mul, hstate.2, one_smul]

/-- For a normalized finite Kraus channel, unit return probability on
the coordinate pure states and the uniform pure state selects the identity
channel. Preparations and readback use the committed Born pairing. -/
theorem finite_pure_return_probabilities_select_lueders_one
    (hn : n ≠ 0) {Φ : CMat n →ₗ[ℂ] CMat n}
    (hK : HasKrausFormWithEffect Φ 1)
    (hbasis : ∀ j : Fin n,
      bornWeight (Φ (vecMulVec (Pi.single j 1) (star (Pi.single j 1))))
        (vecMulVec (Pi.single j 1) (star (Pi.single j 1))) = 1)
    (hcoherent : bornWeight
      (Φ ((n : ℂ)⁻¹ • vecMulVec (fun _ => 1) (star (fun _ => 1))))
      ((n : ℂ)⁻¹ • vecMulVec (fun _ => 1) (star (fun _ => 1))) = 1) :
    Φ = luedersOutcomeMap (1 : CMat n) := by
  apply normalized_finite_pure_probes_select_lueders_one hn hK
  · intro j
    have hstate : IsState (vecMulVec (Pi.single j 1) (star (Pi.single j 1))) := by
      refine ⟨posSemidef_vecMulVec_self_star _, ?_⟩
      rw [trace_vecMulVec, dotProduct, Finset.sum_eq_single j]
      · simp
      · intro i _ hij
        simp [Pi.single_eq_of_ne hij]
      · simp
    have hform : vecMulVec (Pi.single j (1 : ℂ)) (star (Pi.single j 1)) =
        (1 : ℂ) • vecMulVec (Pi.single j 1) (star (Pi.single j 1)) := by
      rw [one_smul]
    exact fixes_scaled_pure_of_return_probability_one hK hstate
      (isEvent_of_state_scaled_pure hstate 1 _ hform) 1 _ hform (hbasis j)
  · have hne : (fun _ : Fin n => (1 : ℂ)) ≠ 0 := by
      intro h
      have hi := congrFun h ⟨0, Nat.pos_of_ne_zero hn⟩
      simp at hi
    have hstate : IsState
        ((n : ℂ)⁻¹ • vecMulVec (fun _ : Fin n => 1) (star (fun _ => 1))) := by
      have h := smul_vecMulVec_mem_certainStates isEvent_one
        (one_mulVec (fun _ : Fin n => (1 : ℂ))) hne
      simpa [dotProduct] using h.1
    exact fixes_scaled_pure_of_return_probability_one hK hstate
      (isEvent_of_state_scaled_pure hstate (n : ℂ)⁻¹ _ rfl) (n : ℂ)⁻¹ _ rfl hcoherent

end

#print axioms kraus_mul_pure_matrix
#print axioms kraus_mulVec_collinear_of_fixes_pure
#print axioms kraus_eq_smul_one_of_finite_pure_probes
#print axioms kraus_channel_eq_id_of_finite_pure_probes
#print axioms finite_pure_probes_select_lueders_one
#print axioms normalized_finite_pure_probes_select_lueders_one
#print axioms scaled_pure_sandwich
#print axioms state_eq_scaled_pure_of_bornWeight_one
#print axioms fixes_scaled_pure_of_return_probability_one
#print axioms isEvent_of_state_scaled_pure
#print axioms finite_pure_return_probabilities_select_lueders_one

end EventAlgebra
