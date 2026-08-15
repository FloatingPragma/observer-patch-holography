import DiscreteCoulombGreen

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.PositionSpaceMaxwellAction

open OPH.SeamCurrentCarrierQuotient
open OPH.SeamU1HolonomyClassification
open OPH.DiscreteCoulombGreen

/-!
# The local real position-space action on the committed seams (issue #733)

WHAT IS PROVED.  Finite exact real linear algebra on the committed
thirty-seam table: real position-space fields `A : Fin 30 → ℝ` on the
committed seams, real port gauge `χ : Fin 12 → ℝ` acting by
`A + realCoboundary χ` through the cast of the committed rational
incidence.  The field-strength projector
`fieldProjector = 1 - d ∘ G ∘ ∂` (with `G` the real cast of the exact
Green matrix of `DiscreteCoulombGreen`) is idempotent
(`fieldProjector_idempotent`) and self-adjoint for the equal-weight
seam pairing (`fieldProjector_selfAdjoint`); its range is the kernel of
the boundary, the nineteen-dimensional cycle space
(`range_fieldProjector`, `cycleSpace_finrank`), and its kernel is the
image of the coboundary (`ker_fieldProjector`).  Hence
`fieldStrength A = fieldProjector A` is gauge-invariant
(`fieldStrength_gauge_invariant`) and determines `A` exactly up to
gauge (`fieldStrength_determines_up_to_gauge`); the real quotient of
position fields by gauge is linearly equivalent to the cycle space
(`positionQuotientEquivCycles`), the real de Rham companion of the
committed finite `U(1)^19` classification
`seamConnectionQuotientMulEquiv`, with the dimension agreeing with the
committed chord count and the committed rational cycle rank, nineteen
(`cycleSpace_finrank_eq_chord_count`, `cycleSpace_finrank_eq_rational`,
`chord_count_eq_gauss_cycle_rank`).  The sourced quadratic action
`sourcedAction J A = (1/2) Σ (P A)² - Σ J·A` is gauge-invariant exactly
when the source is conserved, `realBoundary J = 0`, both directions
(`sourcedAction_gauge_invariant_iff`); the exact quadratic expansion
holds with linear term `⟨P A - J, h⟩` and quadratic term
`(1/2) Σ (P h)²` (`sourcedAction_expansion`); stationarity holds
exactly at `P A = J` (`stationary_iff`); solvability holds exactly for
conserved sources, with `A = J` a solution (`stationary_solvable_iff`,
`conserved_source_fixed`); any two solutions differ exactly by gauge
(`solutions_differ_by_gauge`); every solution is a global minimum, with
equality exactly on the gauge orbit (`sourced_global_minimum`).  For a
neutral rational charge the real cast of the Part-A Coulomb field is
the unique minimal-energy solution of the real Gauss sector
(`electrostatic_join`, composing the Part-A Thomson principle through
`realCoulombField_cast`).  One composed receipt
`positionSpaceMaxwellCoulomb_receipt` carries these clauses as a single
typed conjunction.

PREMISES.  Every theorem in this file is premise-free finite
mathematics on the committed seam carrier, EXCEPT for one selection
remark: the equal-weight inner products `realSeamInner` and
`realPortInner` used throughout coincide with the PR-20 equal
seam-counting selection — the same forced equal weight per seam that
the register row PR-20 selects for the light-signal surface.  No field
of any premise bundle and no premise-register row is consumed by any
proof here; the coincidence is stated so that an instrument lane
consuming PR-20 reads the same pairing.

WHAT IS NOT PROVED.  Calling `sourcedAction` a Maxwell action,
identifying `J` with a physical current, any amplitude with a photon,
the seam table with physical space, a continuum limit, or a laboratory
readout is NOT proved and remains open under PR-53 and PR-54.  The
word Coulomb names the canonical discrete Gauss solution only.  No
dynamics, no quantization, no covariance, and no measurement statement
appears in this module.

Axiom audit.  Every proof composes the `DiscreteCoulombGreen` receipts
with exact real linear algebra; the module adds no project axiom and
uses no native decision procedure.  The audit lines at the end of the
file show at most `propext`, `Classical.choice`, and `Quot.sound`.
-/

noncomputable section

/-! ## Pairing bookkeeping for the equal-weight inner product -/

theorem realSeamInner_comm (c d : Fin 30 → ℝ) :
    realSeamInner c d = realSeamInner d c :=
  Finset.sum_congr rfl fun _ _ ↦ mul_comm _ _

theorem realPortInner_comm (x y : Fin 12 → ℝ) :
    realPortInner x y = realPortInner y x :=
  Finset.sum_congr rfl fun _ _ ↦ mul_comm _ _

theorem realSeamInner_add_right (a b c : Fin 30 → ℝ) :
    realSeamInner a (b + c) = realSeamInner a b + realSeamInner a c := by
  unfold realSeamInner
  rw [← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun e _ ↦ mul_add _ _ _

theorem realSeamInner_sub_left (a b c : Fin 30 → ℝ) :
    realSeamInner (a - b) c = realSeamInner a c - realSeamInner b c := by
  unfold realSeamInner
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun e _ ↦ sub_mul _ _ _

theorem realSeamInner_zero_left (c : Fin 30 → ℝ) :
    realSeamInner 0 c = 0 := by
  unfold realSeamInner
  simp

theorem realSeamInner_zero_right (c : Fin 30 → ℝ) :
    realSeamInner c 0 = 0 := by
  unfold realSeamInner
  simp

theorem realSeamEnergy_zero : realSeamEnergy (0 : Fin 30 → ℝ) = 0 := by
  unfold realSeamEnergy
  simp

/-! ## The field-strength projector -/

/-- The field-strength projector `1 - d ∘ G ∘ ∂` on real position-space
seam fields, built from the committed incidence and the exact real
Green matrix. -/
def fieldProjector : (Fin 30 → ℝ) →ₗ[ℝ] (Fin 30 → ℝ) :=
  LinearMap.id - realCoboundary ∘ₗ greenMatrixR.mulVecLin ∘ₗ realBoundary

theorem fieldProjector_apply (A : Fin 30 → ℝ) :
    fieldProjector A =
      A - realCoboundary (greenMatrixR.mulVec (realBoundary A)) := by
  simp [fieldProjector, LinearMap.sub_apply, LinearMap.comp_apply,
    LinearMap.id_apply]

/-- The projected field is source-free: its boundary vanishes. -/
theorem realBoundary_fieldProjector (A : Fin 30 → ℝ) :
    realBoundary (fieldProjector A) = 0 := by
  rw [fieldProjector_apply, map_sub]
  have hc : realBoundary (realCoboundary (greenMatrixR.mulVec
      (realBoundary A))) = realLaplacian (greenMatrixR.mulVec
      (realBoundary A)) := rfl
  rw [hc, realLaplacian_eq_mulVec, Matrix.mulVec_mulVec,
    laplacian_mul_green_real, Matrix.sub_mulVec, Matrix.one_mulVec,
    Matrix.smul_mulVec, allOnesR_mulVec]
  funext p
  simp [realBoundary_total A]

/-- Cycle fields are fixed by the projector. -/
theorem fieldProjector_fixes (A : Fin 30 → ℝ)
    (hA : realBoundary A = 0) : fieldProjector A = A := by
  rw [fieldProjector_apply, hA, Matrix.mulVec_zero, map_zero, sub_zero]

/-- **(M1)** The projector is idempotent. -/
theorem fieldProjector_idempotent (A : Fin 30 → ℝ) :
    fieldProjector (fieldProjector A) = fieldProjector A :=
  fieldProjector_fixes _ (realBoundary_fieldProjector A)

/-- Pure-gauge fields are annihilated by the projector; the proof uses
the reversed exact Green identity `G L = 1 - (1/12) J`. -/
theorem fieldProjector_coboundary (χ : Fin 12 → ℝ) :
    fieldProjector (realCoboundary χ) = 0 := by
  rw [fieldProjector_apply]
  have hb : realBoundary (realCoboundary χ) = realLaplacian χ := rfl
  rw [hb, realLaplacian_eq_mulVec, Matrix.mulVec_mulVec,
    green_mul_laplacian_real, Matrix.sub_mulVec, Matrix.one_mulVec,
    Matrix.smul_mulVec, allOnesR_mulVec, map_sub, map_smul]
  have hconst : realCoboundary (fun _ ↦ ∑ q : Fin 12, χ q) = 0 :=
    (realCoboundary_eq_zero_iff _).mpr ⟨_, rfl⟩
  rw [hconst, smul_zero, sub_zero, sub_self]

/-- **(M1)** The kernel of the projector is exactly the gauge
directions, the image of the coboundary. -/
theorem ker_fieldProjector :
    LinearMap.ker fieldProjector = LinearMap.range realCoboundary := by
  ext A
  rw [LinearMap.mem_ker, LinearMap.mem_range]
  constructor
  · intro h
    refine ⟨greenMatrixR.mulVec (realBoundary A), ?_⟩
    have h2 : A - realCoboundary (greenMatrixR.mulVec (realBoundary A))
        = 0 := by
      rw [← fieldProjector_apply]
      exact h
    exact (sub_eq_zero.mp h2).symm
  · rintro ⟨χ, rfl⟩
    exact fieldProjector_coboundary χ

/-- **(M1)** The range of the projector is exactly the cycle space, the
kernel of the boundary. -/
theorem range_fieldProjector :
    LinearMap.range fieldProjector = LinearMap.ker realBoundary := by
  ext A
  rw [LinearMap.mem_range, LinearMap.mem_ker]
  constructor
  · rintro ⟨B, rfl⟩
    exact realBoundary_fieldProjector B
  · intro hA
    exact ⟨A, fieldProjector_fixes A hA⟩

/-- **(M1)** The projector is self-adjoint for the equal-weight seam
pairing, by the exact transposition and the symmetry of the Green
matrix. -/
theorem fieldProjector_selfAdjoint (A B : Fin 30 → ℝ) :
    realSeamInner (fieldProjector A) B =
      realSeamInner A (fieldProjector B) := by
  rw [fieldProjector_apply A, fieldProjector_apply B,
    realSeamInner_sub_left, realSeamInner_comm A
      (B - realCoboundary (greenMatrixR.mulVec (realBoundary B))),
    realSeamInner_sub_left, realCoboundary_boundary_adjoint,
    realCoboundary_boundary_adjoint, realPortInner_green_mulVec,
    realSeamInner_comm B A]
  rw [realPortInner_comm (realBoundary A)
    (greenMatrixR.mulVec (realBoundary B))]

/-! ## Gauge invariance of the field strength and the real quotient -/

/-- The field strength of a real position-space field. -/
def fieldStrength (A : Fin 30 → ℝ) : Fin 30 → ℝ := fieldProjector A

/-- **(M1)** The field strength is gauge-invariant. -/
theorem fieldStrength_gauge_invariant (A : Fin 30 → ℝ)
    (χ : Fin 12 → ℝ) :
    fieldStrength (A + realCoboundary χ) = fieldStrength A := by
  unfold fieldStrength
  rw [map_add, fieldProjector_coboundary, add_zero]

/-- **(M1)** Two fields have the same projector image exactly when they
differ by gauge: the field is determined by its field strength exactly
up to gauge. -/
theorem fieldProjector_eq_iff (A B : Fin 30 → ℝ) :
    fieldProjector A = fieldProjector B ↔
      A - B ∈ LinearMap.range realCoboundary := by
  rw [← ker_fieldProjector, LinearMap.mem_ker, map_sub, sub_eq_zero]

theorem fieldStrength_determines_up_to_gauge (A B : Fin 30 → ℝ) :
    fieldStrength A = fieldStrength B ↔
      A - B ∈ LinearMap.range realCoboundary :=
  fieldProjector_eq_iff A B

/-- **(M1)** The real de Rham companion of the committed `U(1)^19`
classification: position-space fields modulo gauge are linearly
equivalent to the nineteen-dimensional cycle space. -/
def positionQuotientEquivCycles :
    ((Fin 30 → ℝ) ⧸ LinearMap.range realCoboundary) ≃ₗ[ℝ]
      LinearMap.ker realBoundary :=
  (Submodule.quotEquivOfEq _ _ ker_fieldProjector.symm).trans
    ((LinearMap.quotKerEquivRange fieldProjector).trans
      (LinearEquiv.ofEq _ _ range_fieldProjector))

/-- The kernel of the real coboundary is the span of the constant
load. -/
theorem realCoboundary_ker_eq_span :
    LinearMap.ker realCoboundary =
      Submodule.span ℝ {(fun _ ↦ 1 : Fin 12 → ℝ)} := by
  ext φ
  rw [LinearMap.mem_ker, Submodule.mem_span_singleton,
    realCoboundary_eq_zero_iff]
  constructor
  · rintro ⟨c, rfl⟩
    exact ⟨c, funext fun p ↦ by simp⟩
  · rintro ⟨a, rfl⟩
    exact ⟨a, funext fun p ↦ by simp⟩

/-- The gauge directions form an eleven-dimensional subspace. -/
theorem gauge_finrank :
    Module.finrank ℝ (LinearMap.range realCoboundary) = 11 := by
  have h := LinearMap.finrank_range_add_finrank_ker realCoboundary
  rw [realCoboundary_ker_eq_span, finrank_span_singleton] at h
  · simp at h
    omega
  · intro hzero
    have h0 := congrFun hzero 0
    norm_num at h0

/-- **(M1)** The cycle space has exact real dimension nineteen. -/
theorem cycleSpace_finrank :
    Module.finrank ℝ (LinearMap.ker realBoundary) = 19 := by
  have h := LinearMap.finrank_range_add_finrank_ker fieldProjector
  rw [range_fieldProjector, ker_fieldProjector, gauge_finrank] at h
  simp at h
  omega

/-- **(M1)** Dimension agreement with the committed chord count of the
`U(1)^19` classification. -/
theorem cycleSpace_finrank_eq_chord_count :
    Module.finrank ℝ (LinearMap.ker realBoundary) =
      ((treeSeamsᶜ : Finset (Fin 30)).card : ℕ) := by
  rw [cycleSpace_finrank, treeSeams_compl_card]

/-- **(M1)** Dimension agreement with the committed rational Gauss
cycle rank. -/
theorem cycleSpace_finrank_eq_rational :
    Module.finrank ℝ (LinearMap.ker realBoundary) =
      Module.finrank ℚ (LinearMap.ker rationalSeamBoundary) := by
  rw [cycleSpace_finrank, OPH.DiscreteGauss.gauss_cycle_space_finrank]

/-! ## (M2) The sourced action and source conservation -/

/-- The sourced quadratic action on real position-space fields.  The
name records shape only; identifying it with a laboratory Maxwell
action is exactly the open attachment PR-53/PR-54. -/
def sourcedAction (J A : Fin 30 → ℝ) : ℝ :=
  (1 / 2) * realSeamEnergy (fieldProjector A) - realSeamInner J A

/-- **(M2)** The sourced action is gauge-invariant exactly when the
source is conserved, in both directions. -/
theorem sourcedAction_gauge_invariant_iff (J : Fin 30 → ℝ) :
    (∀ (A : Fin 30 → ℝ) (χ : Fin 12 → ℝ),
        sourcedAction J (A + realCoboundary χ) = sourcedAction J A) ↔
      realBoundary J = 0 := by
  constructor
  · intro h
    have hpair : ∀ χ : Fin 12 → ℝ,
        realSeamInner J (realCoboundary χ) = 0 := by
      intro χ
      have h0 := h 0 χ
      rw [zero_add] at h0
      unfold sourcedAction at h0
      rw [fieldProjector_coboundary, map_zero, realSeamInner_zero_right]
        at h0
      linarith
    have h1 := hpair (realBoundary J)
    rw [realSeamInner_comm J, realCoboundary_boundary_adjoint] at h1
    exact (realPortInner_self_eq_zero_iff _).mp h1
  · intro hJ A χ
    unfold sourcedAction
    have hP : fieldProjector (A + realCoboundary χ) = fieldProjector A := by
      rw [map_add, fieldProjector_coboundary, add_zero]
    rw [hP, realSeamInner_add_right]
    have hpair : realSeamInner J (realCoboundary χ) = 0 := by
      rw [realSeamInner_comm, realCoboundary_boundary_adjoint, hJ]
      unfold realPortInner
      simp
    rw [hpair, add_zero]

/-! ## (M3) The exact quadratic expansion and the variational package -/

/-- **(M3)** Exact quadratic expansion of the sourced action: linear
term `⟨P A - J, h⟩`, quadratic term `(1/2) Σ (P h)²`, no remainder. -/
theorem sourcedAction_expansion (J A h : Fin 30 → ℝ) :
    sourcedAction J (A + h) = sourcedAction J A +
      realSeamInner (fieldProjector A - J) h +
      (1 / 2) * realSeamEnergy (fieldProjector h) := by
  unfold sourcedAction
  rw [map_add]
  have hsq : realSeamEnergy (fieldProjector A + fieldProjector h) =
      realSeamEnergy (fieldProjector A) +
        2 * realSeamInner (fieldProjector A) (fieldProjector h) +
        realSeamEnergy (fieldProjector h) := by
    unfold realSeamEnergy realSeamInner
    rw [Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun e _ ↦ ?_
    simp only [Pi.add_apply]
    ring
  have hcross : realSeamInner (fieldProjector A) (fieldProjector h) =
      realSeamInner (fieldProjector A) h := by
    rw [fieldProjector_selfAdjoint A (fieldProjector h),
      fieldProjector_idempotent, ← fieldProjector_selfAdjoint A h]
  have hJadd : realSeamInner J (A + h) =
      realSeamInner J A + realSeamInner J h := realSeamInner_add_right J A h
  have hsub : realSeamInner (fieldProjector A - J) h =
      realSeamInner (fieldProjector A) h - realSeamInner J h :=
    realSeamInner_sub_left _ J h
  rw [hsq, hcross, hJadd, hsub]
  ring

/-- **(M3)** A conserved source is fixed by the projector, so the
linear pairing of the expansion is gauge-covariant. -/
theorem conserved_source_fixed (J : Fin 30 → ℝ)
    (hJ : realBoundary J = 0) : fieldProjector J = J :=
  fieldProjector_fixes J hJ

/-- **(M3)** Stationarity of the sourced action holds exactly at
`P A = J`. -/
theorem stationary_iff (J A : Fin 30 → ℝ) :
    (∀ h : Fin 30 → ℝ, sourcedAction J (A + h) = sourcedAction J A +
        (1 / 2) * realSeamEnergy (fieldProjector h)) ↔
      fieldProjector A = J := by
  constructor
  · intro hs
    have hlin : ∀ h : Fin 30 → ℝ,
        realSeamInner (fieldProjector A - J) h = 0 := by
      intro h
      have h1 := sourcedAction_expansion J A h
      have h2 := hs h
      linarith
    have hself := hlin (fieldProjector A - J)
    rw [realSeamInner_self_eq_energy] at hself
    exact sub_eq_zero.mp ((realSeamEnergy_eq_zero_iff _).mp hself)
  · intro hPA h
    rw [sourcedAction_expansion, hPA, sub_self, realSeamInner_zero_left]
    ring

/-- **(M3)** The stationarity equation is solvable exactly for
conserved sources, and then `A = J` is a solution. -/
theorem stationary_solvable_iff (J : Fin 30 → ℝ) :
    (∃ A : Fin 30 → ℝ, fieldProjector A = J) ↔ realBoundary J = 0 := by
  constructor
  · rintro ⟨A, rfl⟩
    exact realBoundary_fieldProjector A
  · intro hJ
    exact ⟨J, conserved_source_fixed J hJ⟩

/-- **(M3)** Any two solutions of the stationarity equation differ
exactly by gauge. -/
theorem solutions_differ_by_gauge (J A B : Fin 30 → ℝ)
    (hA : fieldProjector A = J) :
    fieldProjector B = J ↔ B - A ∈ LinearMap.range realCoboundary := by
  rw [← hA]
  exact fieldProjector_eq_iff B A

/-- **(M3)** Every solution is a global minimum of the sourced action,
with equality exactly on its gauge orbit. -/
theorem sourced_global_minimum (J A : Fin 30 → ℝ)
    (hA : fieldProjector A = J) (B : Fin 30 → ℝ) :
    sourcedAction J A ≤ sourcedAction J B ∧
      (sourcedAction J B = sourcedAction J A ↔
        B - A ∈ LinearMap.range realCoboundary) := by
  have hab : A + (B - A) = B := by abel
  have hB : sourcedAction J B = sourcedAction J A +
      (1 / 2) * realSeamEnergy (fieldProjector (B - A)) := by
    have h1 := sourcedAction_expansion J A (B - A)
    rw [hab] at h1
    rw [h1, hA, sub_self, realSeamInner_zero_left]
    ring
  refine ⟨?_, ?_, ?_⟩
  · rw [hB]
    have hnn := realSeamEnergy_nonneg (fieldProjector (B - A))
    linarith
  · intro heq
    have hz : realSeamEnergy (fieldProjector (B - A)) = 0 := by
      rw [hB] at heq
      linarith
    have hP0 : fieldProjector (B - A) = 0 :=
      (realSeamEnergy_eq_zero_iff _).mp hz
    rw [← ker_fieldProjector]
    exact LinearMap.mem_ker.mpr hP0
  · intro hmem
    have hmem' : B - A ∈ LinearMap.ker fieldProjector := by
      rw [ker_fieldProjector]
      exact hmem
    rw [hB, LinearMap.mem_ker.mp hmem', realSeamEnergy_zero]
    ring

/-! ## (M4) The electrostatic join with the Part-A Coulomb field -/

/-- The real Coulomb field of a real port load: the same canonical
construction as Part A after scalar extension. -/
def realCoulombField (ρ : Fin 12 → ℝ) : Fin 30 → ℝ :=
  realCoboundary (greenMatrixR.mulVec ρ)

theorem realCoulombField_gauss (ρ : Fin 12 → ℝ)
    (hρ : (∑ p : Fin 12, ρ p) = 0) :
    realBoundary (realCoulombField ρ) = ρ := by
  have h1 : realBoundary (realCoulombField ρ) =
      realLaplacian (greenMatrixR.mulVec ρ) := rfl
  rw [h1, realLaplacian_eq_mulVec, Matrix.mulVec_mulVec,
    laplacian_mul_green_real, Matrix.sub_mulVec, Matrix.one_mulVec,
    Matrix.smul_mulVec, allOnesR_mulVec]
  funext p
  simp [hρ]

theorem realCoulombField_orthogonal (ρ : Fin 12 → ℝ) (c : Fin 30 → ℝ)
    (hc : realBoundary c = 0) :
    realSeamInner (realCoulombField ρ) c = 0 := by
  unfold realCoulombField
  rw [realCoboundary_boundary_adjoint, hc]
  unfold realPortInner
  simp

/-- Real Thomson decomposition: every real Gauss solution splits as the
real Coulomb field plus a cycle with exactly additive energy. -/
theorem realThomson_decomposition (ρ : Fin 12 → ℝ)
    (hρ : (∑ p : Fin 12, ρ p) = 0) (E : Fin 30 → ℝ)
    (hE : realBoundary E = ρ) :
    realSeamEnergy E = realSeamEnergy (realCoulombField ρ) +
      realSeamEnergy (E - realCoulombField ρ) := by
  have hFg : realBoundary (realCoulombField ρ) = ρ :=
    realCoulombField_gauss ρ hρ
  have hcyc : realBoundary (E - realCoulombField ρ) = 0 := by
    rw [map_sub, hE, hFg, sub_self]
  have horth : realSeamInner (realCoulombField ρ)
      (E - realCoulombField ρ) = 0 :=
    realCoulombField_orthogonal ρ _ hcyc
  unfold realSeamEnergy
  have hsplit : ∀ e : Fin 30, E e ^ 2 =
      realCoulombField ρ e ^ 2 + (E - realCoulombField ρ) e ^ 2 +
        2 * (realCoulombField ρ e * (E - realCoulombField ρ) e) := by
    intro e
    simp only [Pi.sub_apply]
    ring
  rw [Finset.sum_congr rfl fun e _ ↦ hsplit e, Finset.sum_add_distrib,
    Finset.sum_add_distrib, ← Finset.mul_sum]
  have hzero : (∑ e : Fin 30,
      realCoulombField ρ e * (E - realCoulombField ρ) e) = 0 := horth
  rw [hzero]
  ring

/-- Real Thomson minimality with the equality case. -/
theorem realThomson_minimum (ρ : Fin 12 → ℝ)
    (hρ : (∑ p : Fin 12, ρ p) = 0) (E : Fin 30 → ℝ)
    (hE : realBoundary E = ρ) :
    realSeamEnergy (realCoulombField ρ) ≤ realSeamEnergy E ∧
      (realSeamEnergy E = realSeamEnergy (realCoulombField ρ) ↔
        E = realCoulombField ρ) := by
  have hdec := realThomson_decomposition ρ hρ E hE
  constructor
  · rw [hdec]
    exact le_add_of_nonneg_right (realSeamEnergy_nonneg _)
  · constructor
    · intro heq
      have hz : realSeamEnergy (E - realCoulombField ρ) = 0 := by linarith
      exact sub_eq_zero.mp ((realSeamEnergy_eq_zero_iff _).mp hz)
    · rintro rfl
      rfl

/-- The real Coulomb field of a cast rational load is the cast of the
Part-A rational Coulomb field: the same canonical object. -/
theorem realCoulombField_cast (ρ : RationalPortLoad) :
    realCoulombField (fun p ↦ ((ρ p : ℚ) : ℝ)) =
      fun s ↦ ((coulombField ρ s : ℚ) : ℝ) := by
  funext s
  have hd : realCoulombField (fun p ↦ ((ρ p : ℚ) : ℝ)) s =
      greenMatrixR.mulVec (fun p ↦ ((ρ p : ℚ) : ℝ)) (seamRight s) -
        greenMatrixR.mulVec (fun p ↦ ((ρ p : ℚ) : ℝ)) (seamLeft s) := rfl
  have hcast : ∀ p : Fin 12,
      greenMatrixR.mulVec (fun q ↦ ((ρ q : ℚ) : ℝ)) p =
        ((greenMatrix.mulVec ρ p : ℚ) : ℝ) := by
    intro p
    have h1 : greenMatrixR.mulVec (fun q ↦ ((ρ q : ℚ) : ℝ)) p =
        ∑ q : Fin 12, greenMatrixR p q * ((ρ q : ℚ) : ℝ) := rfl
    have h2 : greenMatrix.mulVec ρ p =
        ∑ q : Fin 12, greenMatrix p q * ρ q := rfl
    rw [h1, h2]
    push_cast
    exact Finset.sum_congr rfl fun q _ ↦ by rw [greenMatrixR_eq_cast]
  have hq : coulombField ρ s =
      greenMatrix.mulVec ρ (seamRight s) - greenMatrix.mulVec ρ (seamLeft s) :=
    rfl
  rw [hd, hcast, hcast, hq]
  push_cast
  ring

/-- **(M4)** Electrostatic join: for a neutral rational charge, the
cast of the Part-A Coulomb field solves the real Gauss problem, is
orthogonal to every real cycle, and is the unique minimal-energy real
Gauss solution — the Part-A Thomson principle carried into the
position-space sector. -/
theorem electrostatic_join (ρ : RationalPortLoad)
    (hρ : rationalPortTotal ρ = 0) :
    realBoundary (fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) =
        (fun p ↦ ((ρ p : ℚ) : ℝ)) ∧
      (∀ c : Fin 30 → ℝ, realBoundary c = 0 →
        realSeamInner (fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) c = 0) ∧
      ∀ E : Fin 30 → ℝ,
        realBoundary E = (fun p ↦ ((ρ p : ℚ) : ℝ)) →
          realSeamEnergy (fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) ≤
              realSeamEnergy E ∧
            (realSeamEnergy E =
                realSeamEnergy (fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) ↔
              E = fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) := by
  have hneutral : (∑ p : Fin 12, ((ρ p : ℚ) : ℝ)) = 0 := by
    rw [← Rat.cast_sum]
    have h : (∑ p : Fin 12, ρ p) = 0 := hρ
    rw [h, Rat.cast_zero]
  have hcast : (fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) =
      realCoulombField (fun p ↦ ((ρ p : ℚ) : ℝ)) :=
    (realCoulombField_cast ρ).symm
  refine ⟨?_, ?_, ?_⟩
  · rw [hcast]
    exact realCoulombField_gauss _ hneutral
  · intro c hc
    rw [hcast]
    exact realCoulombField_orthogonal _ c hc
  · intro E hE
    rw [hcast]
    exact realThomson_minimum _ hneutral E hE

/-! ## The one citable composed receipt -/

/-- **The position-space Maxwell/Coulomb receipt (issue #733).**  One
typed conjunction of the M1-M4 clauses: the projector is idempotent
and self-adjoint with kernel the gauge directions and range the
nineteen-dimensional cycle space; the field strength is gauge-invariant
and determines the field exactly up to gauge, with the real quotient
equivalent to the cycle space at the committed chord count; the sourced
action is gauge-invariant exactly for conserved sources; the exact
quadratic expansion, stationarity exactly at `P A = J`, solvability
exactly for conserved sources, gauge-uniqueness of solutions, and the
global-minimum property; and the electrostatic join with the Part-A
Coulomb field.  Premise-free finite mathematics; PR-53 and PR-54 stay
open, and the equal-weight pairing coincides with the PR-20 selection
without consuming it. -/
theorem positionSpaceMaxwellCoulomb_receipt :
    (∀ A : Fin 30 → ℝ,
      fieldProjector (fieldProjector A) = fieldProjector A)
    ∧ (∀ A B : Fin 30 → ℝ,
      realSeamInner (fieldProjector A) B =
        realSeamInner A (fieldProjector B))
    ∧ LinearMap.ker fieldProjector = LinearMap.range realCoboundary
    ∧ LinearMap.range fieldProjector = LinearMap.ker realBoundary
    ∧ (∀ (A : Fin 30 → ℝ) (χ : Fin 12 → ℝ),
      fieldStrength (A + realCoboundary χ) = fieldStrength A)
    ∧ (∀ A B : Fin 30 → ℝ,
      fieldStrength A = fieldStrength B ↔
        A - B ∈ LinearMap.range realCoboundary)
    ∧ Nonempty (((Fin 30 → ℝ) ⧸ LinearMap.range realCoboundary) ≃ₗ[ℝ]
        LinearMap.ker realBoundary)
    ∧ Module.finrank ℝ (LinearMap.ker realBoundary) = 19
    ∧ Module.finrank ℝ (LinearMap.ker realBoundary) =
        ((treeSeamsᶜ : Finset (Fin 30)).card : ℕ)
    ∧ Module.finrank ℝ (LinearMap.ker realBoundary) =
        Module.finrank ℚ (LinearMap.ker rationalSeamBoundary)
    ∧ (∀ J : Fin 30 → ℝ,
      (∀ (A : Fin 30 → ℝ) (χ : Fin 12 → ℝ),
          sourcedAction J (A + realCoboundary χ) = sourcedAction J A) ↔
        realBoundary J = 0)
    ∧ (∀ J A h : Fin 30 → ℝ,
      sourcedAction J (A + h) = sourcedAction J A +
        realSeamInner (fieldProjector A - J) h +
        (1 / 2) * realSeamEnergy (fieldProjector h))
    ∧ (∀ J A : Fin 30 → ℝ,
      (∀ h : Fin 30 → ℝ, sourcedAction J (A + h) = sourcedAction J A +
          (1 / 2) * realSeamEnergy (fieldProjector h)) ↔
        fieldProjector A = J)
    ∧ (∀ J : Fin 30 → ℝ,
      (∃ A : Fin 30 → ℝ, fieldProjector A = J) ↔ realBoundary J = 0)
    ∧ (∀ J A B : Fin 30 → ℝ, fieldProjector A = J →
      (fieldProjector B = J ↔ B - A ∈ LinearMap.range realCoboundary))
    ∧ (∀ J A : Fin 30 → ℝ, fieldProjector A = J → ∀ B : Fin 30 → ℝ,
      sourcedAction J A ≤ sourcedAction J B ∧
        (sourcedAction J B = sourcedAction J A ↔
          B - A ∈ LinearMap.range realCoboundary))
    ∧ ∀ ρ : RationalPortLoad, rationalPortTotal ρ = 0 →
      realBoundary (fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) =
          (fun p ↦ ((ρ p : ℚ) : ℝ)) ∧
        (∀ c : Fin 30 → ℝ, realBoundary c = 0 →
          realSeamInner (fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) c = 0) ∧
        ∀ E : Fin 30 → ℝ,
          realBoundary E = (fun p ↦ ((ρ p : ℚ) : ℝ)) →
            realSeamEnergy (fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) ≤
                realSeamEnergy E ∧
              (realSeamEnergy E =
                  realSeamEnergy (fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) ↔
                E = fun s ↦ ((coulombField ρ s : ℚ) : ℝ)) :=
  ⟨fieldProjector_idempotent, fieldProjector_selfAdjoint,
    ker_fieldProjector, range_fieldProjector,
    fieldStrength_gauge_invariant, fieldStrength_determines_up_to_gauge,
    ⟨positionQuotientEquivCycles⟩, cycleSpace_finrank,
    cycleSpace_finrank_eq_chord_count, cycleSpace_finrank_eq_rational,
    sourcedAction_gauge_invariant_iff, sourcedAction_expansion,
    stationary_iff, stationary_solvable_iff, solutions_differ_by_gauge,
    fun J A hA B ↦ sourced_global_minimum J A hA B, electrostatic_join⟩

end

end OPH.PositionSpaceMaxwellAction

/- Axiom audit: the `DiscreteCoulombGreen` receipts and exact real
linear algebra only.  Expected axioms per line: at most `propext`,
`Classical.choice`, `Quot.sound`.  No native decision procedure is
used. -/

#print axioms OPH.PositionSpaceMaxwellAction.fieldProjector_apply
#print axioms OPH.PositionSpaceMaxwellAction.realBoundary_fieldProjector
#print axioms OPH.PositionSpaceMaxwellAction.fieldProjector_fixes
#print axioms OPH.PositionSpaceMaxwellAction.fieldProjector_idempotent
#print axioms OPH.PositionSpaceMaxwellAction.fieldProjector_coboundary
#print axioms OPH.PositionSpaceMaxwellAction.ker_fieldProjector
#print axioms OPH.PositionSpaceMaxwellAction.range_fieldProjector
#print axioms OPH.PositionSpaceMaxwellAction.fieldProjector_selfAdjoint
#print axioms OPH.PositionSpaceMaxwellAction.fieldStrength_gauge_invariant
#print axioms OPH.PositionSpaceMaxwellAction.fieldProjector_eq_iff
#print axioms OPH.PositionSpaceMaxwellAction.fieldStrength_determines_up_to_gauge
#print axioms OPH.PositionSpaceMaxwellAction.positionQuotientEquivCycles
#print axioms OPH.PositionSpaceMaxwellAction.realCoboundary_ker_eq_span
#print axioms OPH.PositionSpaceMaxwellAction.gauge_finrank
#print axioms OPH.PositionSpaceMaxwellAction.cycleSpace_finrank
#print axioms OPH.PositionSpaceMaxwellAction.cycleSpace_finrank_eq_chord_count
#print axioms OPH.PositionSpaceMaxwellAction.cycleSpace_finrank_eq_rational
#print axioms OPH.PositionSpaceMaxwellAction.sourcedAction_gauge_invariant_iff
#print axioms OPH.PositionSpaceMaxwellAction.sourcedAction_expansion
#print axioms OPH.PositionSpaceMaxwellAction.conserved_source_fixed
#print axioms OPH.PositionSpaceMaxwellAction.stationary_iff
#print axioms OPH.PositionSpaceMaxwellAction.stationary_solvable_iff
#print axioms OPH.PositionSpaceMaxwellAction.solutions_differ_by_gauge
#print axioms OPH.PositionSpaceMaxwellAction.sourced_global_minimum
#print axioms OPH.PositionSpaceMaxwellAction.realCoulombField_gauss
#print axioms OPH.PositionSpaceMaxwellAction.realCoulombField_orthogonal
#print axioms OPH.PositionSpaceMaxwellAction.realThomson_decomposition
#print axioms OPH.PositionSpaceMaxwellAction.realThomson_minimum
#print axioms OPH.PositionSpaceMaxwellAction.realCoulombField_cast
#print axioms OPH.PositionSpaceMaxwellAction.electrostatic_join
#print axioms OPH.PositionSpaceMaxwellAction.positionSpaceMaxwellCoulomb_receipt
