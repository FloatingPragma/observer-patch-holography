import ObserverPatchHolography.EinsteinBranch.SmallBall

/-!
# Generic dark sector and the modular-anomaly source

The rest-frame Einstein relation of `SmallBall.lean` reads the stress
contraction as one real number.  On the Einstein branch that number is the
modular-charge stress of the screen: every source carrying modular charge
enters it, and the universal-coupling clause says that nothing else
distinguishes one source from another in the geometry response.  This module
records the algebraic consequences.

* A finite source family splits into luminous and non-luminous members.  Under
  universal coupling the geometry responds to the total modular charge, so the
  luminous-only Einstein relation holds exactly when the non-luminous charge
  vanishes (`luminousOnlyRelation_iff_darkStress_eq_zero`).  A geometric
  excess over the luminous prediction forces a non-luminous source with
  nonzero charge (`exists_dark_source_of_geometric_excess`).  This is the
  sense in which a dark sector is generic rather than optional.
* The candidate dark source is the anomalous modular energy carried on overlap
  collars, the Einstein branch's controlled finite-stage remainder.  Its
  rest-frame stress is bounded by the collar recovery defect
  (`anomalousStress_abs_le`), vanishes exactly at full recovery
  (`anomalousStress_eq_zero_of_recovered`), and tends to zero along any
  family whose recovery defect tends to zero at fixed radius and bounded
  collar constant (`anomalousStress_tendsto_zero`).  The last statement is the
  screening interface: where recovery saturates, the dark source switches off.
* A conserved comoving charge dilutes as the inverse cube of the scale factor
  (`comovingDensity_mul_scale_cubed`), the cold-background scaling.

What is not proved here: that the physical stress contraction is the total
modular charge (universal coupling is a premise field), that the anomalous
modular energy is localized on collars with a quotient-invariant receipt, the
size of the collar constant, the radial profile of the anomalous charge, and
that the recovery defect vanishes in any physical regime.  Those are the
open targets named in the dark-matter paper.
-/

namespace OPH.EinsteinBranch

open Filter
open scoped Topology BigOperators

noncomputable section

/-! ## Luminous and non-luminous modular charge -/

/-- A finite family of screen sources.  `charge i` is the rest-frame
modular-charge contribution of source `i` to the stress contraction;
`luminous i` records whether the source carries the electromagnetic readout. -/
structure SourceFamily (ι : Type) [Fintype ι] where
  charge : ι → ℝ
  luminous : ι → Bool

variable {ι : Type} [Fintype ι]

/-- Total modular charge of the family. -/
def totalStress (S : SourceFamily ι) : ℝ := ∑ i, S.charge i

/-- Modular charge carried by luminous sources. -/
def luminousStress (S : SourceFamily ι) : ℝ :=
  ∑ i ∈ Finset.univ.filter (fun i => S.luminous i = true), S.charge i

/-- Modular charge carried by non-luminous sources. -/
def darkStress (S : SourceFamily ι) : ℝ :=
  ∑ i ∈ Finset.univ.filter (fun i => S.luminous i = false), S.charge i

theorem totalStress_eq_luminous_add_dark (S : SourceFamily ι) :
    totalStress S = luminousStress S + darkStress S := by
  unfold totalStress luminousStress darkStress
  rw [← Finset.sum_filter_add_sum_filter_not Finset.univ
    (fun i => S.luminous i = true)]
  congr 1
  apply Finset.sum_congr _ (fun _ _ => rfl)
  ext i
  simp [Bool.not_eq_true]

/-- Universal coupling as a premise: the stress contraction entering the
rest-frame Einstein relation is the total modular charge of the family,
with no readout-dependent weight. -/
def UniversalCoupling (P : SmallBallPremises) (S : SourceFamily ι) : Prop :=
  P.stressContraction = totalStress S

/-- Under universal coupling the geometry responds to luminous and
non-luminous charge alike. -/
theorem geometry_responds_to_all_charge (P : SmallBallPremises)
    (S : SourceFamily ι) (h : UniversalCoupling P S) :
    P.geometryContraction =
      8 * Real.pi * P.G * (luminousStress S + darkStress S) := by
  rw [restFrameEinsteinRelation P, h, totalStress_eq_luminous_add_dark]

/-- The luminous-only Einstein relation holds exactly when the non-luminous
modular charge vanishes. -/
theorem luminousOnlyRelation_iff_darkStress_eq_zero (P : SmallBallPremises)
    (S : SourceFamily ι) (h : UniversalCoupling P S) :
    P.geometryContraction = 8 * Real.pi * P.G * luminousStress S ↔
      darkStress S = 0 := by
  rw [geometry_responds_to_all_charge P S h]
  have hpos : (0 : ℝ) < 8 * Real.pi * P.G := by
    have := Real.pi_pos
    have := P.G_pos
    positivity
  constructor
  · intro heq
    have : 8 * Real.pi * P.G * darkStress S = 0 := by linarith
    rcases mul_eq_zero.mp this with h0 | h0
    · exact absurd h0 (ne_of_gt hpos)
    · exact h0
  · intro h0
    rw [h0, add_zero]

/-- A geometric excess over the luminous prediction forces a non-luminous
source with nonzero modular charge. -/
theorem exists_dark_source_of_geometric_excess (P : SmallBallPremises)
    (S : SourceFamily ι) (h : UniversalCoupling P S)
    (hex : P.geometryContraction ≠ 8 * Real.pi * P.G * luminousStress S) :
    ∃ i, S.luminous i = false ∧ S.charge i ≠ 0 := by
  by_contra hnone
  push Not at hnone
  apply hex
  rw [luminousOnlyRelation_iff_darkStress_eq_zero P S h]
  unfold darkStress
  apply Finset.sum_eq_zero
  intro i hi
  exact hnone i (Finset.mem_filter.mp hi).2

/-- Nonnegative non-luminous charges give a nonnegative, hence attractive,
dark stress in the weak-field equation. -/
theorem darkStress_nonneg_of_charges_nonneg (S : SourceFamily ι)
    (h : ∀ i, S.luminous i = false → 0 ≤ S.charge i) :
    0 ≤ darkStress S := by
  unfold darkStress
  apply Finset.sum_nonneg
  intro i hi
  exact h i (Finset.mem_filter.mp hi).2

/-! ## The anomalous modular energy as a bounded source -/

/-- The finite-stage anomalous modular energy carried on an overlap collar,
with the controlled bound of the Einstein branch: `|E| ≤ C * eta`, where
`eta` is the collar recovery defect and `C` the collar constant. -/
structure AnomalyRemainder where
  ell : ℝ
  ell_pos : 0 < ell
  C : ℝ
  C_nonneg : 0 ≤ C
  eta : ℝ
  eta_nonneg : 0 ≤ eta
  anomalousEnergy : ℝ
  bound : |anomalousEnergy| ≤ C * eta

/-- Rest-frame anomalous stress: the small-ball kernel inverse applied to the
anomalous modular energy. -/
def anomalousStress (A : AnomalyRemainder) : ℝ :=
  15 / (8 * Real.pi ^ 2 * A.ell ^ 4) * A.anomalousEnergy

theorem smallBallKernelInverse_pos (A : AnomalyRemainder) :
    0 < 15 / (8 * Real.pi ^ 2 * A.ell ^ 4) := by
  have := Real.pi_pos
  have := A.ell_pos
  positivity

/-- The anomalous stress is bounded by the recovery defect. -/
theorem anomalousStress_abs_le (A : AnomalyRemainder) :
    |anomalousStress A| ≤ 15 * A.C * A.eta / (8 * Real.pi ^ 2 * A.ell ^ 4) := by
  unfold anomalousStress
  rw [abs_mul, abs_of_pos (smallBallKernelInverse_pos A)]
  have hk := smallBallKernelInverse_pos A
  calc 15 / (8 * Real.pi ^ 2 * A.ell ^ 4) * |A.anomalousEnergy|
      ≤ 15 / (8 * Real.pi ^ 2 * A.ell ^ 4) * (A.C * A.eta) :=
        mul_le_mul_of_nonneg_left A.bound (le_of_lt hk)
    _ = 15 * A.C * A.eta / (8 * Real.pi ^ 2 * A.ell ^ 4) := by ring

/-- Exact recovery switches the dark source off. -/
theorem anomalousStress_eq_zero_of_recovered (A : AnomalyRemainder)
    (h : A.eta = 0) : anomalousStress A = 0 := by
  have hb := A.bound
  rw [h, mul_zero] at hb
  have hE : A.anomalousEnergy = 0 := abs_nonpos_iff.mp hb
  unfold anomalousStress
  rw [hE, mul_zero]

/-- Screening interface: along a family with fixed radius and bounded collar
constant, the anomalous stress tends to zero whenever the recovery defect
does.  That the recovery defect vanishes in any physical regime is not
proved here. -/
theorem anomalousStress_tendsto_zero (A : ℕ → AnomalyRemainder) (ℓ Cmax : ℝ)
    (hℓ : ∀ n, (A n).ell = ℓ) (hC : ∀ n, (A n).C ≤ Cmax)
    (heta : Tendsto (fun n => (A n).eta) atTop (𝓝 0)) :
    Tendsto (fun n => anomalousStress (A n)) atTop (𝓝 0) := by
  have hℓpos : 0 < ℓ := by
    have := (A 0).ell_pos
    rwa [hℓ 0] at this
  have hCmax : 0 ≤ Cmax := le_trans (A 0).C_nonneg (hC 0)
  set k : ℝ := 15 / (8 * Real.pi ^ 2 * ℓ ^ 4) with hk
  have hkpos : 0 < k := by
    have := Real.pi_pos
    positivity
  -- squeeze: |stress n| ≤ k * Cmax * eta n → 0
  have hbound : ∀ n, |anomalousStress (A n)| ≤ k * Cmax * (A n).eta := by
    intro n
    have h1 := anomalousStress_abs_le (A n)
    rw [hℓ n] at h1
    calc |anomalousStress (A n)|
        ≤ 15 * (A n).C * (A n).eta / (8 * Real.pi ^ 2 * ℓ ^ 4) := h1
      _ = k * (A n).C * (A n).eta := by rw [hk]; ring
      _ ≤ k * Cmax * (A n).eta := by
          apply mul_le_mul_of_nonneg_right _ (A n).eta_nonneg
          exact mul_le_mul_of_nonneg_left (hC n) (le_of_lt hkpos)
  have hlim : Tendsto (fun n => k * Cmax * (A n).eta) atTop (𝓝 0) := by
    have := heta.const_mul (k * Cmax)
    simpa using this
  rw [tendsto_zero_iff_abs_tendsto_zero]
  refine squeeze_zero (fun n => abs_nonneg _) hbound hlim

/-! ## Cold scaling of a conserved comoving charge -/

/-- A conserved modular charge `Q` in a comoving cell whose physical volume
scales as `a ^ 3`. -/
structure ComovingCharge where
  Q : ℝ
  a : ℕ → ℝ
  a_pos : ∀ n, 0 < a n

/-- Physical charge density at step `n`. -/
def comovingDensity (X : ComovingCharge) (n : ℕ) : ℝ := X.Q / X.a n ^ 3

/-- Conserved charge dilutes as the inverse cube of the scale factor. -/
theorem comovingDensity_mul_scale_cubed (X : ComovingCharge) (n : ℕ) :
    comovingDensity X n * X.a n ^ 3 = X.Q := by
  unfold comovingDensity
  have h : X.a n ^ 3 ≠ 0 := pow_ne_zero 3 (ne_of_gt (X.a_pos n))
  exact div_mul_cancel₀ X.Q h

/-! ## Per-theorem axiom audit -/

#print axioms totalStress_eq_luminous_add_dark
#print axioms geometry_responds_to_all_charge
#print axioms luminousOnlyRelation_iff_darkStress_eq_zero
#print axioms exists_dark_source_of_geometric_excess
#print axioms darkStress_nonneg_of_charges_nonneg
#print axioms anomalousStress_abs_le
#print axioms anomalousStress_eq_zero_of_recovered
#print axioms anomalousStress_tendsto_zero
#print axioms comovingDensity_mul_scale_cubed

end

end OPH.EinsteinBranch
