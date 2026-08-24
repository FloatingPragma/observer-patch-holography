import ObserverPatchHolography.EinsteinBranch.DarkSector

/-!
# Deep-profile characterization from atomic premises

`DarkSector.lean` commits to the linear enclosed-mass profile
`M_A(r) = r * sqrt(Mb * a0 / G)` as a single functional-form premise
(premise `pre:linear` of the dark-matter paper).  This module replaces the
functional form by two atomic premises and derives the profile from them.

* `DeepRegimeLaw` packages a deep-regime enclosed-mass law
  `M : Mb -> r -> M_A(r)` with two structural premises: degree-one scale
  covariance in the radius (`scaleCovariant`: a sparsely recorded deep cut
  retains no length scale, so rescaling the radius rescales the enclosed
  anomalous mass by the same factor) and quadrature composition in the
  baryonic source (`quadratureClosed`: per-source collar defects are
  independent, variances add, so anomalous amplitudes add in quadrature),
  together with the regularity premises `nonneg`, `monoMass`, `nondeg`.
* `linear_on_pos_of_additive_monotone` is a self-contained Cauchy lemma: a
  function additive and monotone on the positive reals is linear there,
  with a nonnegative slope.
* `deepProfile_characterization`: every `DeepRegimeLaw` carries a unique
  constant `a0 > 0` with `M Mb r = r * sqrt(Mb * a0 / G)` for all positive
  `Mb` and `r`.  The committed profile is forced by the two atomic premises
  together with the stated regularity premises; `monoMass` is load-bearing
  in the Cauchy step and `nonneg`/`nondeg` in the square-root extraction.
* Bridge theorems transfer `deep_radial_acceleration`,
  `circularSpeedSq_const` and `baryonic_tully_fisher` from `DeepProfile` to
  any `DeepRegimeLaw` (`law_deep_radial_acceleration`,
  `law_circularSpeedSq_const`, `law_baryonic_tully_fisher`,
  `law_baryonic_tully_fisher_exists`).
* `matchingRadius` `r_M = sqrt(G * Mb / a0)`: the baryonic and anomalous
  accelerations both equal `a0` there
  (`baryonicAcceleration_matchingRadius`,
  `anomalousAcceleration_matchingRadius`), and the anomalous acceleration
  dominates the baryonic one exactly outside `r_M`
  (`anomalous_dominates_iff`).
* Load-bearing receipts: `additiveLaw` satisfies every premise except
  `quadratureClosed`, `sqrtLaw` satisfies every premise except
  `scaleCovariant`, and `unitLaw` is a full explicit inhabitant with
  `G = 1` and constant `a0 = 1` (`unitLaw_char`).  Each of the two atomic
  premises is therefore necessary, and the premise bundle is jointly
  satisfiable.

What is not proved here: `scaleCovariant` and `quadratureClosed` are
premises, motivated by variance additivity of independent collar defects
and by the absence of a retained scale on sparsely recorded deep cuts;
neither is derived from the axioms of the event-algebra corpus.  The value
of `a0` is not determined: the characterization fixes the functional form
and the uniqueness of the constant, not its magnitude.  The transition
region near and inside `r_M` is not modelled.  Comparable deep-limit
scale-invariance and square-root laws in the literature are discussed in
the owning paper; no novelty claim is made here.
-/

namespace OPH.EinsteinBranch

noncomputable section

/-! ## The deep-regime law: two atomic premises -/

/-- A deep-regime enclosed-mass law.  `M Mb r` is the anomalous mass enclosed
within radius `r` around a baryonic source of mass `Mb`.  The two atomic
premises are `scaleCovariant` (no retained length scale in the deep regime)
and `quadratureClosed` (independent per-source collar defects compose in
variance, so amplitudes add in quadrature).  `nonneg`, `monoMass` and
`nondeg` are regularity premises: an enclosed mass is nonnegative, does not
decrease with the source, and is somewhere nonzero. -/
structure DeepRegimeLaw where
  G : ℝ
  G_pos : 0 < G
  M : ℝ → ℝ → ℝ
  nonneg : ∀ Mb r : ℝ, 0 < Mb → 0 < r → 0 ≤ M Mb r
  scaleCovariant : ∀ Mb lam r : ℝ, 0 < Mb → 0 < lam → 0 < r →
    M Mb (lam * r) = lam * M Mb r
  quadratureClosed : ∀ m1 m2 r : ℝ, 0 < m1 → 0 < m2 → 0 < r →
    M (m1 + m2) r ^ 2 = M m1 r ^ 2 + M m2 r ^ 2
  monoMass : ∀ m1 m2 r : ℝ, 0 < m1 → m1 ≤ m2 → 0 < r → M m1 r ≤ M m2 r
  nondeg : ∃ m r : ℝ, 0 < m ∧ 0 < r ∧ 0 < M m r

/-! ## A Cauchy lemma on the positive reals -/

/-- A function additive and monotone on the positive reals is linear there,
with a nonnegative slope.  Self-contained Cauchy argument: additivity gives
linearity at positive rational points, monotonicity squeezes the value at
every positive real between nearby rational values. -/
theorem linear_on_pos_of_additive_monotone (f : ℝ → ℝ)
    (hadd : ∀ x y : ℝ, 0 < x → 0 < y → f (x + y) = f x + f y)
    (hmono : ∀ x y : ℝ, 0 < x → x ≤ y → f x ≤ f y) :
    ∃ c : ℝ, 0 ≤ c ∧ ∀ x : ℝ, 0 < x → f x = c * x := by
  -- nonnegativity on positives from additivity and monotonicity
  have hnonneg : ∀ x : ℝ, 0 < x → 0 ≤ f x := by
    intro x hx
    have h1 : f x ≤ f (x + x) := hmono x (x + x) hx (by linarith)
    have h2 : f (x + x) = f x + f x := hadd x x hx hx
    linarith
  -- natural-number scaling
  have hnat : ∀ n : ℕ, 1 ≤ n → ∀ x : ℝ, 0 < x →
      f ((n : ℝ) * x) = (n : ℝ) * f x := by
    intro n hn
    induction n, hn using Nat.le_induction with
    | base => intro x hx; norm_num
    | succ k hk ih =>
      intro x hx
      have hk1 : (1 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk
      have hkx : 0 < (k : ℝ) * x :=
        mul_pos (lt_of_lt_of_le zero_lt_one hk1) hx
      have hcast : ((k + 1 : ℕ) : ℝ) * x = (k : ℝ) * x + x := by
        push_cast; ring
      rw [hcast, hadd _ x hkx hx, ih x hx]
      push_cast; ring
  have hnatpt : ∀ n : ℕ, 1 ≤ n → f ((n : ℝ)) = (n : ℝ) * f 1 := by
    intro n hn
    have h := hnat n hn 1 one_pos
    rwa [mul_one] at h
  -- rational scaling
  have hrat : ∀ q : ℚ, 0 < q → f ((q : ℝ)) = f 1 * (q : ℝ) := by
    intro q hq
    have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
    have hnum : 0 < q.num := Rat.num_pos.mpr hq
    have hdpos : 0 < q.den := q.den_pos
    have hdR : (0 : ℝ) < (q.den : ℝ) := by exact_mod_cast hdpos
    have hkey : (q.den : ℝ) * (q : ℝ) = (q.num : ℝ) := by
      rw [Rat.cast_def]
      field_simp
    have hfd : f ((q.den : ℝ) * (q : ℝ)) = (q.den : ℝ) * f ((q : ℝ)) :=
      hnat q.den (by omega) _ hqR
    have hfnum : f ((q.num : ℝ)) = (q.num : ℝ) * f 1 := by
      have h1 : (1 : ℕ) ≤ q.num.toNat := by omega
      have hcast : ((q.num.toNat : ℕ) : ℝ) = (q.num : ℝ) := by
        exact_mod_cast Int.toNat_of_nonneg hnum.le
      have h2 := hnatpt q.num.toNat h1
      rwa [hcast] at h2
    rw [hkey, hfnum] at hfd
    have hfq : f ((q : ℝ)) = (q.num : ℝ) * f 1 / (q.den : ℝ) := by
      rw [eq_div_iff (ne_of_gt hdR), mul_comm]
      linarith [hfd]
    rw [hfq, Rat.cast_def]
    ring
  refine ⟨f 1, hnonneg 1 one_pos, ?_⟩
  intro x hx
  by_contra hne
  rcases lt_or_gt_of_ne hne with hlt | hgt
  · -- f x < f 1 * x is excluded by a rational point just below x
    have hfx0 : 0 ≤ f x := hnonneg x hx
    have hcpos : 0 < f 1 := by
      rcases (hnonneg 1 one_pos).eq_or_lt with heq | h
      · exfalso
        rw [← heq, zero_mul] at hlt
        linarith
      · exact h
    have hdivlt : f x / f 1 < x := (div_lt_iff₀ hcpos).mpr (by linarith)
    obtain ⟨q, hq1, hq2⟩ := exists_rat_btwn hdivlt
    have hq0R : (0 : ℝ) < (q : ℝ) :=
      lt_of_le_of_lt (div_nonneg hfx0 hcpos.le) hq1
    have hq0 : (0 : ℚ) < q := by exact_mod_cast hq0R
    have hmono1 : f ((q : ℝ)) ≤ f x := hmono _ x hq0R hq2.le
    rw [hrat q hq0] at hmono1
    have h3 := (div_lt_iff₀ hcpos).mp hq1
    linarith
  · -- f x > f 1 * x is excluded by a rational point just above x
    rcases (hnonneg 1 one_pos).eq_or_lt with heq | hcpos
    · obtain ⟨q, hq⟩ := exists_rat_gt x
      have hq0R : (0 : ℝ) < (q : ℝ) := lt_trans hx hq
      have hq0 : (0 : ℚ) < q := by exact_mod_cast hq0R
      have hmono1 : f x ≤ f ((q : ℝ)) := hmono x _ hx hq.le
      rw [hrat q hq0, ← heq, zero_mul] at hmono1
      rw [← heq, zero_mul] at hgt
      linarith
    · have hdivgt : x < f x / f 1 := (lt_div_iff₀ hcpos).mpr (by linarith)
      obtain ⟨q, hq1, hq2⟩ := exists_rat_btwn hdivgt
      have hq0R : (0 : ℝ) < (q : ℝ) := lt_trans hx hq1
      have hq0 : (0 : ℚ) < q := by exact_mod_cast hq0R
      have hmono1 : f x ≤ f ((q : ℝ)) := hmono x _ hx hq1.le
      rw [hrat q hq0] at hmono1
      have h3 := (lt_div_iff₀ hcpos).mp hq2
      linarith

/-! ## The characterization theorem -/

/-- Every deep-regime law carries a unique acceleration constant `a0 > 0`
realizing the committed linear profile `M Mb r = r * sqrt(Mb * a0 / G)`.
Scale covariance factorizes the law as `r` times a radial coefficient,
quadrature closure makes the squared coefficient additive in the source,
and the Cauchy lemma forces the square root. -/
theorem deepProfile_characterization (L : DeepRegimeLaw) :
    ∃! a0 : ℝ, 0 < a0 ∧ ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      L.M Mb r = r * Real.sqrt (Mb * a0 / L.G) := by
  have hGne : L.G ≠ 0 := ne_of_gt L.G_pos
  -- radial factorization from scale covariance
  have hfact : ∀ Mb r : ℝ, 0 < Mb → 0 < r → L.M Mb r = r * L.M Mb 1 := by
    intro Mb r hMb hr
    have h := L.scaleCovariant Mb r 1 hMb hr one_pos
    rwa [mul_one] at h
  -- the squared radial coefficient is additive and monotone on positives
  have hgadd : ∀ x y : ℝ, 0 < x → 0 < y →
      L.M (x + y) 1 ^ 2 = L.M x 1 ^ 2 + L.M y 1 ^ 2 := fun x y hx hy =>
    L.quadratureClosed x y 1 hx hy one_pos
  have hgmono : ∀ x y : ℝ, 0 < x → x ≤ y → L.M x 1 ^ 2 ≤ L.M y 1 ^ 2 := by
    intro x y hx hxy
    have h1 : 0 ≤ L.M x 1 := L.nonneg x 1 hx one_pos
    have h2 : L.M x 1 ≤ L.M y 1 := L.monoMass x y 1 hx hxy one_pos
    rw [pow_two, pow_two]
    exact mul_self_le_mul_self h1 h2
  obtain ⟨lam, hlam0, hglin0⟩ :=
    linear_on_pos_of_additive_monotone (fun m => L.M m 1 ^ 2) hgadd hgmono
  have hglin : ∀ x : ℝ, 0 < x → L.M x 1 ^ 2 = lam * x := fun x hx =>
    hglin0 x hx
  -- nondegeneracy forces a strictly positive slope
  obtain ⟨m0, r0, hm0, hr0, hM0⟩ := L.nondeg
  have hf0 : 0 < L.M m0 1 := by
    have heq := hfact m0 r0 hm0 hr0
    rw [heq] at hM0
    by_contra hle
    push Not at hle
    nlinarith [mul_le_mul_of_nonneg_left hle hr0.le]
  have hlampos : 0 < lam := by
    have hgm0 : 0 < L.M m0 1 ^ 2 := pow_pos hf0 2
    rw [hglin m0 hm0] at hgm0
    by_contra hle
    push Not at hle
    nlinarith [mul_le_mul_of_nonneg_right hle hm0.le]
  refine ⟨lam * L.G, ⟨mul_pos hlampos L.G_pos, ?_⟩, ?_⟩
  · intro Mb r hMb hr
    have hsq : Mb * (lam * L.G) / L.G = lam * Mb := by
      rw [div_eq_iff hGne]; ring
    rw [hfact Mb r hMb hr, hsq, ← hglin Mb hMb,
      Real.sqrt_sq (L.nonneg Mb 1 hMb one_pos)]
  · rintro a0' ⟨ha0', hchar'⟩
    have hM11 : L.M 1 1 ^ 2 = lam := by
      have h := hglin 1 one_pos
      rwa [mul_one] at h
    have h1 := hchar' 1 1 one_pos one_pos
    simp only [one_mul] at h1
    have h2 : L.M 1 1 ^ 2 = a0' / L.G := by
      rw [h1, Real.sq_sqrt (div_nonneg ha0'.le L.G_pos.le)]
    have h3 : a0' / L.G = lam := by rw [← h2, hM11]
    rw [div_eq_iff hGne] at h3
    exact h3

/-! ## Bridge to the committed profile -/

/-- The committed `DeepProfile` carrying the characterized constant of a
deep-regime law at a given baryonic mass. -/
def toDeepProfile (L : DeepRegimeLaw) (a0 Mb : ℝ) (ha0 : 0 < a0)
    (hMb : 0 < Mb) : DeepProfile :=
  ⟨L.G, L.G_pos, a0, ha0, Mb, hMb⟩

/-- Under the characterization the law's enclosed mass is the committed
`enclosedMass` of the bridged profile. -/
theorem toDeepProfile_enclosedMass (L : DeepRegimeLaw) (a0 : ℝ)
    (ha0 : 0 < a0)
    (hchar : ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      L.M Mb r = r * Real.sqrt (Mb * a0 / L.G))
    (Mb r : ℝ) (hMb : 0 < Mb) (hr : 0 < r) :
    L.M Mb r = enclosedMass (toDeepProfile L a0 Mb ha0 hMb) r := by
  rw [hchar Mb r hMb hr]
  rfl

/-- Deep radial-acceleration relation transferred to a deep-regime law:
`G M(Mb, r) / r^2 = sqrt((G Mb / r^2) a0)`. -/
theorem law_deep_radial_acceleration (L : DeepRegimeLaw) (a0 : ℝ)
    (ha0 : 0 < a0)
    (hchar : ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      L.M Mb r = r * Real.sqrt (Mb * a0 / L.G))
    (Mb r : ℝ) (hMb : 0 < Mb) (hr : 0 < r) :
    L.G * L.M Mb r / r ^ 2 = Real.sqrt (L.G * Mb / r ^ 2 * a0) := by
  have hD := deep_radial_acceleration (toDeepProfile L a0 Mb ha0 hMb) r hr
  unfold anomalousAcceleration baryonicAcceleration at hD
  rw [← toDeepProfile_enclosedMass L a0 ha0 hchar Mb r hMb hr] at hD
  exact hD

/-- Flat rotation curve transferred to a deep-regime law: the anomalous
circular speed squared is radius independent. -/
theorem law_circularSpeedSq_const (L : DeepRegimeLaw) (a0 : ℝ)
    (ha0 : 0 < a0)
    (hchar : ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      L.M Mb r = r * Real.sqrt (Mb * a0 / L.G))
    (Mb r : ℝ) (hMb : 0 < Mb) (hr : 0 < r) :
    L.G * L.M Mb r / r ^ 2 * r = Real.sqrt (L.G * Mb * a0) := by
  have hD := circularSpeedSq_const (toDeepProfile L a0 Mb ha0 hMb) r hr
  unfold circularSpeedSq anomalousAcceleration at hD
  rw [← toDeepProfile_enclosedMass L a0 ha0 hchar Mb r hMb hr] at hD
  exact hD

/-- Baryonic Tully--Fisher relation transferred to a deep-regime law:
`v^4 = G Mb a0`. -/
theorem law_baryonic_tully_fisher (L : DeepRegimeLaw) (a0 : ℝ)
    (ha0 : 0 < a0)
    (hchar : ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      L.M Mb r = r * Real.sqrt (Mb * a0 / L.G))
    (Mb r : ℝ) (hMb : 0 < Mb) (hr : 0 < r) :
    (L.G * L.M Mb r / r ^ 2 * r) ^ 2 = L.G * Mb * a0 := by
  have hD := baryonic_tully_fisher (toDeepProfile L a0 Mb ha0 hMb) r hr
  unfold circularSpeedSq anomalousAcceleration at hD
  rw [← toDeepProfile_enclosedMass L a0 ha0 hchar Mb r hMb hr] at hD
  exact hD

/-- Composed form: every deep-regime law satisfies the baryonic
Tully--Fisher relation for some positive constant. -/
theorem law_baryonic_tully_fisher_exists (L : DeepRegimeLaw) :
    ∃ a0 : ℝ, 0 < a0 ∧ ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      (L.G * L.M Mb r / r ^ 2 * r) ^ 2 = L.G * Mb * a0 := by
  obtain ⟨a0, ⟨ha0, hchar⟩, -⟩ := deepProfile_characterization L
  exact ⟨a0, ha0, fun Mb r hMb hr =>
    law_baryonic_tully_fisher L a0 ha0 hchar Mb r hMb hr⟩

/-! ## The matching radius -/

/-- Matching radius `r_M = sqrt(G Mb / a0)` of a committed profile: the
scale at which the baryonic and anomalous accelerations exchange
dominance. -/
def matchingRadius (D : DeepProfile) : ℝ := Real.sqrt (D.G * D.Mb / D.a0)

/-- The matching radius of a committed profile is positive. -/
theorem matchingRadius_pos (D : DeepProfile) : 0 < matchingRadius D := by
  have hG := D.G_pos; have hMb := D.Mb_pos; have ha0 := D.a0_pos
  unfold matchingRadius
  exact Real.sqrt_pos.mpr (by positivity)

/-- The squared matching radius is `G Mb / a0` exactly. -/
theorem matchingRadius_sq (D : DeepProfile) :
    matchingRadius D ^ 2 = D.G * D.Mb / D.a0 := by
  have hG := D.G_pos; have hMb := D.Mb_pos; have ha0 := D.a0_pos
  unfold matchingRadius
  exact Real.sq_sqrt (by positivity)

/-- The flat-curve speed scale factorizes through the matching radius:
`sqrt(G Mb a0) = a0 * r_M`. -/
theorem sqrt_eq_a0_mul_matchingRadius (D : DeepProfile) :
    Real.sqrt (D.G * D.Mb * D.a0) = D.a0 * matchingRadius D := by
  have hG := D.G_pos; have hMb := D.Mb_pos; have ha0 := D.a0_pos
  have ha0ne : D.a0 ≠ 0 := ne_of_gt ha0
  unfold matchingRadius
  rw [show D.G * D.Mb * D.a0 = D.a0 ^ 2 * (D.G * D.Mb / D.a0) by
        field_simp,
      Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq ha0.le]

/-- At the matching radius the baryonic acceleration equals `a0`. -/
theorem baryonicAcceleration_matchingRadius (D : DeepProfile) :
    baryonicAcceleration D (matchingRadius D) = D.a0 := by
  have hG := D.G_pos; have hMb := D.Mb_pos; have ha0 := D.a0_pos
  have hGMbne : D.G * D.Mb ≠ 0 := by positivity
  have ha0ne : D.a0 ≠ 0 := ne_of_gt ha0
  unfold baryonicAcceleration
  rw [matchingRadius_sq]
  field_simp

/-- At the matching radius the anomalous acceleration equals `a0`. -/
theorem anomalousAcceleration_matchingRadius (D : DeepProfile) :
    anomalousAcceleration D (matchingRadius D) = D.a0 := by
  rw [anomalousAcceleration_eq D _ (matchingRadius_pos D),
    sqrt_eq_a0_mul_matchingRadius, mul_div_assoc,
    div_self (ne_of_gt (matchingRadius_pos D)), mul_one]

/-- Deep and Newtonian regimes exchange dominance exactly at the matching
radius: the anomalous acceleration exceeds the baryonic one at `r` if and
only if `r_M < r`. -/
theorem anomalous_dominates_iff (D : DeepProfile) (r : ℝ) (hr : 0 < r) :
    anomalousAcceleration D r > baryonicAcceleration D r ↔
      matchingRadius D < r := by
  have hG := D.G_pos; have hMb := D.Mb_pos; have ha0 := D.a0_pos
  have hrM := matchingRadius_pos D
  have ha0ne : D.a0 ≠ 0 := ne_of_gt ha0
  have hGMb : D.G * D.Mb = D.a0 * matchingRadius D ^ 2 := by
    rw [matchingRadius_sq]
    field_simp
  rw [anomalousAcceleration_eq D r hr, sqrt_eq_a0_mul_matchingRadius]
  unfold baryonicAcceleration
  rw [hGMb, gt_iff_lt,
    div_lt_div_iff₀ (by positivity : (0 : ℝ) < r ^ 2) hr]
  constructor
  · intro h
    by_contra hle
    push Not at hle
    nlinarith [mul_le_mul_of_nonneg_left
      (mul_le_mul_of_nonneg_right hle hr.le) (mul_pos ha0 hrM).le]
  · intro h
    nlinarith [mul_lt_mul_of_pos_left h (mul_pos (mul_pos ha0 hrM) hr)]

/-! ## Load-bearing receipts

`additiveLaw` shows `quadratureClosed` is necessary, `sqrtLaw` shows
`scaleCovariant` is necessary, and `unitLaw` shows the premise bundle is
jointly satisfiable. -/

/-- Enclosed mass proportional to the source: satisfies every premise of
`DeepRegimeLaw` except quadrature closure. -/
def additiveLaw (c : ℝ) : ℝ → ℝ → ℝ := fun Mb r => r * (Mb * c)

theorem additiveLaw_scaleCovariant (c : ℝ) :
    ∀ Mb lam r : ℝ, 0 < Mb → 0 < lam → 0 < r →
      additiveLaw c Mb (lam * r) = lam * additiveLaw c Mb r := by
  intro Mb lam r _ _ _
  unfold additiveLaw
  ring

theorem additiveLaw_nonneg (c : ℝ) (hc : 0 < c) :
    ∀ Mb r : ℝ, 0 < Mb → 0 < r → 0 ≤ additiveLaw c Mb r := by
  intro Mb r hMb hr
  unfold additiveLaw
  positivity

theorem additiveLaw_monoMass (c : ℝ) (hc : 0 < c) :
    ∀ m1 m2 r : ℝ, 0 < m1 → m1 ≤ m2 → 0 < r →
      additiveLaw c m1 r ≤ additiveLaw c m2 r := by
  intro m1 m2 r hm1 hm12 hr
  unfold additiveLaw
  exact mul_le_mul_of_nonneg_left
    (mul_le_mul_of_nonneg_right hm12 hc.le) hr.le

theorem additiveLaw_nondeg (c : ℝ) (hc : 0 < c) :
    ∃ m r : ℝ, 0 < m ∧ 0 < r ∧ 0 < additiveLaw c m r :=
  ⟨1, 1, one_pos, one_pos, by simpa [additiveLaw] using hc⟩

/-- The additive law violates quadrature closure: amplitudes that add
linearly do not add in quadrature. -/
theorem additiveLaw_not_quadratureClosed (c : ℝ) (hc : 0 < c) :
    ¬ ∀ m1 m2 r : ℝ, 0 < m1 → 0 < m2 → 0 < r →
      additiveLaw c (m1 + m2) r ^ 2 =
        additiveLaw c m1 r ^ 2 + additiveLaw c m2 r ^ 2 := by
  intro h
  have h1 := h 1 1 1 one_pos one_pos one_pos
  unfold additiveLaw at h1
  nlinarith [h1, pow_pos hc 2]

/-- Enclosed mass growing as `r^2`: satisfies every premise of
`DeepRegimeLaw` except scale covariance. -/
def sqrtLaw (lam : ℝ) : ℝ → ℝ → ℝ := fun Mb r => r ^ 2 * Real.sqrt (lam * Mb)

theorem sqrtLaw_quadratureClosed (lam : ℝ) (hlam : 0 < lam) :
    ∀ m1 m2 r : ℝ, 0 < m1 → 0 < m2 → 0 < r →
      sqrtLaw lam (m1 + m2) r ^ 2 =
        sqrtLaw lam m1 r ^ 2 + sqrtLaw lam m2 r ^ 2 := by
  intro m1 m2 r h1 h2 hr
  unfold sqrtLaw
  have e0 : Real.sqrt (lam * (m1 + m2)) ^ 2 = lam * (m1 + m2) :=
    Real.sq_sqrt (by positivity)
  have e1 : Real.sqrt (lam * m1) ^ 2 = lam * m1 :=
    Real.sq_sqrt (by positivity)
  have e2 : Real.sqrt (lam * m2) ^ 2 = lam * m2 :=
    Real.sq_sqrt (by positivity)
  rw [mul_pow, mul_pow, mul_pow, e0, e1, e2]
  ring

theorem sqrtLaw_nonneg (lam : ℝ) (hlam : 0 < lam) :
    ∀ Mb r : ℝ, 0 < Mb → 0 < r → 0 ≤ sqrtLaw lam Mb r := by
  intro Mb r hMb hr
  unfold sqrtLaw
  positivity

theorem sqrtLaw_monoMass (lam : ℝ) (hlam : 0 < lam) :
    ∀ m1 m2 r : ℝ, 0 < m1 → m1 ≤ m2 → 0 < r →
      sqrtLaw lam m1 r ≤ sqrtLaw lam m2 r := by
  intro m1 m2 r hm1 hm12 hr
  unfold sqrtLaw
  exact mul_le_mul_of_nonneg_left
    (Real.sqrt_le_sqrt (mul_le_mul_of_nonneg_left hm12 hlam.le))
    (sq_nonneg r)

theorem sqrtLaw_nondeg (lam : ℝ) (hlam : 0 < lam) :
    ∃ m r : ℝ, 0 < m ∧ 0 < r ∧ 0 < sqrtLaw lam m r :=
  ⟨1, 1, one_pos, one_pos, by
    unfold sqrtLaw
    rw [one_pow, one_mul, mul_one]
    exact Real.sqrt_pos.mpr hlam⟩

/-- The `r^2` law violates scale covariance: doubling the radius
quadruples the enclosed mass instead of doubling it. -/
theorem sqrtLaw_not_scaleCovariant (lam : ℝ) (hlam : 0 < lam) :
    ¬ ∀ Mb l r : ℝ, 0 < Mb → 0 < l → 0 < r →
      sqrtLaw lam Mb (l * r) = l * sqrtLaw lam Mb r := by
  intro h
  have h1 := h 1 2 1 one_pos two_pos one_pos
  unfold sqrtLaw at h1
  have hs : 0 < Real.sqrt (lam * 1) := Real.sqrt_pos.mpr (by rwa [mul_one])
  nlinarith [h1, hs]

/-- Full explicit inhabitant of `DeepRegimeLaw`: `G = 1` and
`M Mb r = r * sqrt Mb`.  The premise bundle is jointly satisfiable. -/
def unitLaw : DeepRegimeLaw where
  G := 1
  G_pos := one_pos
  M := fun Mb r => r * Real.sqrt Mb
  nonneg := fun Mb r hMb hr => by positivity
  scaleCovariant := fun Mb lam r _ _ _ => by ring
  quadratureClosed := fun m1 m2 r h1 h2 hr => by
    have e0 : Real.sqrt (m1 + m2) ^ 2 = m1 + m2 := Real.sq_sqrt (by positivity)
    have e1 : Real.sqrt m1 ^ 2 = m1 := Real.sq_sqrt h1.le
    have e2 : Real.sqrt m2 ^ 2 = m2 := Real.sq_sqrt h2.le
    rw [mul_pow, mul_pow, mul_pow, e0, e1, e2]
    ring
  monoMass := fun m1 m2 r hm1 hm12 hr =>
    mul_le_mul_of_nonneg_left (Real.sqrt_le_sqrt hm12) hr.le
  nondeg := ⟨1, 1, one_pos, one_pos, by norm_num [Real.sqrt_one]⟩

/-- The inhabitant realizes the characterization with `a0 = 1`. -/
theorem unitLaw_char :
    ∀ Mb r : ℝ, 0 < Mb → 0 < r →
      unitLaw.M Mb r = r * Real.sqrt (Mb * 1 / unitLaw.G) := by
  intro Mb r _ _
  show r * Real.sqrt Mb = r * Real.sqrt (Mb * 1 / 1)
  rw [mul_one, div_one]

/-! ## Per-theorem axiom audit -/

#print axioms linear_on_pos_of_additive_monotone
#print axioms deepProfile_characterization
#print axioms toDeepProfile_enclosedMass
#print axioms law_deep_radial_acceleration
#print axioms law_circularSpeedSq_const
#print axioms law_baryonic_tully_fisher
#print axioms law_baryonic_tully_fisher_exists
#print axioms matchingRadius_pos
#print axioms matchingRadius_sq
#print axioms sqrt_eq_a0_mul_matchingRadius
#print axioms baryonicAcceleration_matchingRadius
#print axioms anomalousAcceleration_matchingRadius
#print axioms anomalous_dominates_iff
#print axioms additiveLaw_scaleCovariant
#print axioms additiveLaw_nonneg
#print axioms additiveLaw_monoMass
#print axioms additiveLaw_nondeg
#print axioms additiveLaw_not_quadratureClosed
#print axioms sqrtLaw_quadratureClosed
#print axioms sqrtLaw_nonneg
#print axioms sqrtLaw_monoMass
#print axioms sqrtLaw_nondeg
#print axioms sqrtLaw_not_scaleCovariant
#print axioms unitLaw_char

end

end OPH.EinsteinBranch
