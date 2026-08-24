import ObserverPatchHolography.EinsteinBranch.DeepProfileClosure

/-!
# Per-cut collar composition and the acceleration-constant dictionary

`DeepProfileClosure.lean` characterizes the deep-regime enclosed-mass
profile from two atomic premises, taking quadrature closure as given.
This module derives quadrature closure, and the enclosed-mass law itself,
from a more primitive per-cut collar model, and extracts the magnitude
dictionary for the acceleration constant.

* `PerCutCollarModel` packages Newton's constant `G`, a cut linear
  density `n` (settled codimension-one cuts per unit radius around a
  point source, so `n * r` cuts are nested within radius `r`), and a
  per-cut defect variance function `V` with four premises: variances of
  independent per-source collar defects add (`V_additive`), the variance
  does not decrease with the source (`V_mono`), it is nonnegative on
  positive sources (`V_nonneg`), and it is somewhere nonzero
  (`V_nondeg`).  The per-cut anomalous amplitude is
  `cutAmplitude Mb = sqrt (V Mb)` and the induced enclosed mass is
  `inducedM Mb r = n * r * cutAmplitude Mb`: each of the nested cuts
  carries the amplitude once.
* `varianceSlope`: by the Cauchy lemma
  `linear_on_pos_of_additive_monotone` the variance function is linear on
  positive sources, `V m = c * m` for a unique slope `c`, strictly
  positive by `V_nondeg`.
* `perCutLaw` and `perCut_isDeepRegimeLaw`: the induced enclosed-mass
  law inhabits `DeepRegimeLaw` at the same `G`.  Quadrature closure, an
  atomic premise there, is derived here from `V_additive` through the
  identity `inducedM Mb r ^ 2 = (n * r) ^ 2 * V Mb` (`inducedM_sq`):
  squared induced masses are variances up to the fixed factor
  `(n * r) ^ 2`, and variances add.  Scale covariance holds by
  construction, and the remaining premises follow from `V_mono`,
  `V_nonneg`, `V_nondeg`.
* `perCut_a0_dictionary` and `perCut_a0_unique`: for a model with
  variance slope `c`, the induced law realizes the committed profile with
  constant `a0 = G * n ^ 2 * c` exactly, and by
  `deepProfile_characterization` this is the unique constant of the
  induced law.
* `perCut_nonidentifiability`: two explicit models with different cut
  densities (`thinCutModel` with `n = 1`, slope `4` and `denseCutModel`
  with `n = 2`, slope `1`, both at `G = 1`) have equal product
  `G * n ^ 2 * c` and induce the same enclosed-mass law
  (`thinCutModel_varianceSlope`, `denseCutModel_varianceSlope`).  The
  dictionary therefore determines only the product `n ^ 2 * c`; neither
  factor is separately fixed by the deep law.
* `unitCollarModel` is a full explicit inhabitant with `G = 1`, `n = 1`,
  `V m = m`; its induced law is the enclosed-mass profile of `unitLaw`
  from `DeepProfileClosure.lean` (`unitCollarModel_matches_unitLaw`) and
  its dictionary constant is `G * n ^ 2 * 1`
  (`unitCollarModel_dictionary`).

What is not proved here: variance additivity of independent per-source
collar defects and the linear cut count are premises of the model, not
derived from the axioms of the event-algebra corpus.  No numeric value of
`n`, `c`, or `a0` is asserted or derived; the dictionary reduces the
magnitude question to the collar density and the per-cut variance
constant, both open.  No comparison with data is made here.
-/

namespace OPH.EinsteinBranch

noncomputable section

/-! ## The per-cut collar model -/

/-- A per-cut collar model around a point source.  `n` is the cut linear
density: the number of settled codimension-one collar cuts per unit
radius, so `n * r` cuts are nested within radius `r`.  `V m` is the
variance of the collar defect contributed by an independent source of
mass `m`.  The premises: variances of independent per-source defects add
(`V_additive`), the variance does not decrease with the source
(`V_mono`), it is nonnegative on positive sources (`V_nonneg`), and it
is somewhere nonzero (`V_nondeg`). -/
structure PerCutCollarModel where
  G : ℝ
  G_pos : 0 < G
  n : ℝ
  n_pos : 0 < n
  V : ℝ → ℝ
  V_additive : ∀ m1 m2 : ℝ, 0 < m1 → 0 < m2 → V (m1 + m2) = V m1 + V m2
  V_mono : ∀ m1 m2 : ℝ, 0 < m1 → m1 ≤ m2 → V m1 ≤ V m2
  V_nonneg : ∀ m : ℝ, 0 < m → 0 ≤ V m
  V_nondeg : ∃ m : ℝ, 0 < m ∧ 0 < V m

/-- Per-cut anomalous amplitude: the standard deviation of the collar
defect sourced by a baryonic mass `Mb`. -/
def cutAmplitude (P : PerCutCollarModel) (Mb : ℝ) : ℝ :=
  Real.sqrt (P.V Mb)

/-- Induced enclosed anomalous mass: the `n * r` cuts nested within
radius `r` each carry the per-cut amplitude once. -/
def inducedM (P : PerCutCollarModel) (Mb r : ℝ) : ℝ :=
  P.n * r * cutAmplitude P Mb

/-- The squared induced mass is the per-cut variance scaled by the
squared cut count.  This identity turns variance additivity into
quadrature closure. -/
theorem inducedM_sq (P : PerCutCollarModel) (Mb r : ℝ) (hMb : 0 < Mb) :
    inducedM P Mb r ^ 2 = (P.n * r) ^ 2 * P.V Mb := by
  unfold inducedM cutAmplitude
  rw [mul_pow, Real.sq_sqrt (P.V_nonneg Mb hMb)]

/-! ## The variance slope -/

/-- The variance function of a per-cut collar model is linear on positive
sources with a unique, strictly positive slope: `V m = c * m`.
Existence and linearity are the Cauchy lemma
`linear_on_pos_of_additive_monotone` applied to `V_additive` and
`V_mono`; strict positivity of the slope is `V_nondeg`. -/
theorem varianceSlope (P : PerCutCollarModel) :
    ∃! c : ℝ, 0 < c ∧ ∀ m : ℝ, 0 < m → P.V m = c * m := by
  obtain ⟨c, _, hlin⟩ :=
    linear_on_pos_of_additive_monotone P.V P.V_additive P.V_mono
  have hcpos : 0 < c := by
    obtain ⟨m0, hm0, hVm0⟩ := P.V_nondeg
    rw [hlin m0 hm0] at hVm0
    by_contra hle
    push Not at hle
    nlinarith [mul_le_mul_of_nonneg_right hle hm0.le]
  refine ⟨c, ⟨hcpos, hlin⟩, ?_⟩
  rintro c' ⟨-, hlin'⟩
  have h1 := hlin 1 one_pos
  have h1' := hlin' 1 one_pos
  rw [mul_one] at h1 h1'
  rw [← h1']
  exact h1

/-! ## The induced law is a deep-regime law -/

/-- The deep-regime law induced by a per-cut collar model: same `G`,
enclosed mass `inducedM`.  Quadrature closure is derived from variance
additivity via `inducedM_sq`; scale covariance is linearity of the cut
count in the radius. -/
def perCutLaw (P : PerCutCollarModel) : DeepRegimeLaw where
  G := P.G
  G_pos := P.G_pos
  M := inducedM P
  nonneg := fun Mb r _ hr => by
    have hn := P.n_pos
    unfold inducedM cutAmplitude
    positivity
  scaleCovariant := fun Mb lam r _ _ _ => by
    unfold inducedM
    ring
  quadratureClosed := fun m1 m2 r h1 h2 _ => by
    rw [inducedM_sq P (m1 + m2) r (by linarith), inducedM_sq P m1 r h1,
      inducedM_sq P m2 r h2, P.V_additive m1 m2 h1 h2]
    ring
  monoMass := fun m1 m2 r hm1 hm12 hr => by
    unfold inducedM cutAmplitude
    exact mul_le_mul_of_nonneg_left
      (Real.sqrt_le_sqrt (P.V_mono m1 m2 hm1 hm12))
      (mul_nonneg P.n_pos.le hr.le)
  nondeg := by
    obtain ⟨m0, hm0, hVm0⟩ := P.V_nondeg
    exact ⟨m0, 1, hm0, one_pos, by
      unfold inducedM cutAmplitude
      rw [mul_one]
      exact mul_pos P.n_pos (Real.sqrt_pos.mpr hVm0)⟩

/-- The induced deep-regime law carries exactly the induced enclosed mass
and the model's `G`. -/
theorem perCut_isDeepRegimeLaw (P : PerCutCollarModel) :
    (perCutLaw P).M = inducedM P ∧ (perCutLaw P).G = P.G :=
  ⟨rfl, rfl⟩

/-! ## The magnitude dictionary -/

/-- Magnitude dictionary: for a model with variance slope `c`, the
induced law realizes the committed deep profile with constant
`a0 = G * n ^ 2 * c`, which is positive. -/
theorem perCut_a0_dictionary (P : PerCutCollarModel) (c : ℝ) (hc : 0 < c)
    (hlin : ∀ m : ℝ, 0 < m → P.V m = c * m) :
    0 < P.G * P.n ^ 2 * c ∧
      ∀ Mb r : ℝ, 0 < Mb → 0 < r →
        inducedM P Mb r = r * Real.sqrt (Mb * (P.G * P.n ^ 2 * c) / P.G) := by
  constructor
  · have hG := P.G_pos
    have hn := P.n_pos
    positivity
  · intro Mb r hMb _
    have hGne : P.G ≠ 0 := ne_of_gt P.G_pos
    have key : Mb * (P.G * P.n ^ 2 * c) / P.G = P.n ^ 2 * (c * Mb) := by
      rw [div_eq_iff hGne]
      ring
    rw [key, Real.sqrt_mul (sq_nonneg P.n) (c * Mb),
      Real.sqrt_sq P.n_pos.le]
    unfold inducedM cutAmplitude
    rw [hlin Mb hMb]
    ring

/-- Uniqueness corollary: by `deepProfile_characterization`, any constant
realizing the committed profile for the induced law equals
`G * n ^ 2 * c`.  The characterized acceleration constant of the induced
law is therefore exactly the dictionary value. -/
theorem perCut_a0_unique (P : PerCutCollarModel) (c : ℝ) (hc : 0 < c)
    (hlin : ∀ m : ℝ, 0 < m → P.V m = c * m) (a0 : ℝ) (ha0 : 0 < a0)
    (hchar : ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      (perCutLaw P).M Mb r = r * Real.sqrt (Mb * a0 / (perCutLaw P).G)) :
    a0 = P.G * P.n ^ 2 * c := by
  obtain ⟨b, -, huniq⟩ := deepProfile_characterization (perCutLaw P)
  have h1 : a0 = b := huniq a0 ⟨ha0, hchar⟩
  have h2 : P.G * P.n ^ 2 * c = b := by
    obtain ⟨hpos, hdict⟩ := perCut_a0_dictionary P c hc hlin
    refine huniq _ ⟨hpos, ?_⟩
    intro Mb r hMb hr
    show inducedM P Mb r = r * Real.sqrt (Mb * (P.G * P.n ^ 2 * c) / P.G)
    exact hdict Mb r hMb hr
  exact h1.trans h2.symm

/-! ## Non-identifiability of the factors -/

/-- Model with cut density `n = 1` and variance `V m = 4 * m`. -/
def thinCutModel : PerCutCollarModel where
  G := 1
  G_pos := one_pos
  n := 1
  n_pos := one_pos
  V := fun m => 4 * m
  V_additive := by
    intro m1 m2 _ _
    show 4 * (m1 + m2) = 4 * m1 + 4 * m2
    ring
  V_mono := by
    intro m1 m2 _ h12
    show 4 * m1 ≤ 4 * m2
    linarith
  V_nonneg := by
    intro m hm
    show 0 ≤ 4 * m
    linarith
  V_nondeg := ⟨1, one_pos, by norm_num⟩

/-- Model with cut density `n = 2` and variance `V m = m`. -/
def denseCutModel : PerCutCollarModel where
  G := 1
  G_pos := one_pos
  n := 2
  n_pos := two_pos
  V := fun m => m
  V_additive := fun _ _ _ _ => rfl
  V_mono := fun _ _ _ h => h
  V_nonneg := fun _ hm => hm.le
  V_nondeg := ⟨1, one_pos, one_pos⟩

/-- The thin model has variance slope `4`. -/
theorem thinCutModel_varianceSlope :
    ∀ m : ℝ, 0 < m → thinCutModel.V m = 4 * m :=
  fun _ _ => rfl

/-- The dense model has variance slope `1`. -/
theorem denseCutModel_varianceSlope :
    ∀ m : ℝ, 0 < m → denseCutModel.V m = 1 * m :=
  fun m _ => (one_mul m).symm

/-- Non-identifiability receipt: the two models have different cut
densities and equal dictionary product `G * n ^ 2 * c`, and they induce
the same enclosed-mass law.  The deep law fixes only the product
`n ^ 2 * c`; neither factor separately. -/
theorem perCut_nonidentifiability :
    thinCutModel.n ≠ denseCutModel.n ∧
      thinCutModel.G * thinCutModel.n ^ 2 * 4 =
        denseCutModel.G * denseCutModel.n ^ 2 * 1 ∧
      ∀ Mb r : ℝ, inducedM thinCutModel Mb r = inducedM denseCutModel Mb r := by
  refine ⟨?_, ?_, ?_⟩
  · show (1 : ℝ) ≠ 2
    norm_num
  · show (1 : ℝ) * 1 ^ 2 * 4 = 1 * 2 ^ 2 * 1
    norm_num
  · intro Mb r
    show 1 * r * Real.sqrt (4 * Mb) = 2 * r * Real.sqrt Mb
    have h42 : Real.sqrt 4 = 2 := by
      rw [show (4 : ℝ) = 2 ^ 2 by norm_num,
        Real.sqrt_sq (by norm_num : (0 : ℝ) ≤ 2)]
    rw [Real.sqrt_mul (by norm_num : (0 : ℝ) ≤ 4) Mb, h42]
    ring

/-! ## Explicit unit inhabitant -/

/-- Full explicit inhabitant: `G = 1`, `n = 1`, `V m = m`. -/
def unitCollarModel : PerCutCollarModel where
  G := 1
  G_pos := one_pos
  n := 1
  n_pos := one_pos
  V := fun m => m
  V_additive := fun _ _ _ _ => rfl
  V_mono := fun _ _ _ h => h
  V_nonneg := fun _ hm => hm.le
  V_nondeg := ⟨1, one_pos, one_pos⟩

/-- The unit collar model induces the enclosed-mass profile of `unitLaw`
from `DeepProfileClosure.lean`. -/
theorem unitCollarModel_matches_unitLaw (Mb r : ℝ) :
    inducedM unitCollarModel Mb r = unitLaw.M Mb r := by
  show 1 * r * Real.sqrt Mb = r * Real.sqrt Mb
  ring

/-- Dictionary instance for the unit collar model: its induced law
realizes the committed profile with constant `G * n ^ 2 * 1`. -/
theorem unitCollarModel_dictionary :
    ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      inducedM unitCollarModel Mb r =
        r * Real.sqrt (Mb *
          (unitCollarModel.G * unitCollarModel.n ^ 2 * 1) /
          unitCollarModel.G) :=
  (perCut_a0_dictionary unitCollarModel 1 one_pos
    (fun m _ => (one_mul m).symm)).2

/-! ## Per-theorem axiom audit -/

#print axioms inducedM_sq
#print axioms varianceSlope
#print axioms perCut_isDeepRegimeLaw
#print axioms perCut_a0_dictionary
#print axioms perCut_a0_unique
#print axioms thinCutModel_varianceSlope
#print axioms denseCutModel_varianceSlope
#print axioms perCut_nonidentifiability
#print axioms unitCollarModel_matches_unitLaw
#print axioms unitCollarModel_dictionary

end

end OPH.EinsteinBranch
