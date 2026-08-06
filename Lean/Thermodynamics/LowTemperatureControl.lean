import FiniteConditionalRepair

namespace OPH.Thermodynamics

/-!
# Refinement-uniform low-temperature control

This module proves the finite low-temperature control receipt of the
conditional thermodynamics package at the representation level:

* the mass of the finite Gibbs family off the exact minimum set of the
  energy is bounded by `card X * exp (-beta * gap)` at every real
  inverse temperature, where `gap` is the exact least positive energy
  excess over the minimum;
* the off-minimum mass tends to zero as `beta` tends to infinity;
* a declared refinement family with minimum-preserving maps, one
  uniform positive gap lower bound, and one cardinality bound is
  controlled by the single bound `cardBound * exp (-beta * gapBound)`
  across every member simultaneously, giving uniform concentration;
* the conditional-resampling reference law is the constant-energy
  member of the reference-relative Gibbs family at every `beta`, and
  the repair kernel fixes it; this is the degenerate instance receipt
  on the four-law package;
* a two-level witness carries a nonconstant energy with exact gap and
  an exact off-minimum mass formula;
* a negative control: a refinement family whose gaps shrink to zero
  admits no uniform gap lower bound and violates the uniform
  concentration conclusion, exhibited at `beta = 1` on member `r = 2`.

Imported premises, not established here:

* the identification of the energy with the realized global objective
  and of the reference with the shared optimizer reference (B12
  receipts 1 and 2);
* the identification of `beta` with a physical inverse temperature and
  the energy-clock calibration (E5, #703);
* every continuum and physical-realization statement.
-/

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]
variable {B : Type*} [DecidableEq B]

section MinimumSet

variable (E : Ω → ℝ)

/-- The exact minimum energy value over the finite state space. -/
noncomputable def minEnergy [Nonempty Ω] : ℝ :=
  Finset.univ.inf' Finset.univ_nonempty E

omit [DecidableEq Ω] in
theorem minEnergy_le [Nonempty Ω] (x : Ω) : minEnergy E ≤ E x :=
  Finset.inf'_le E (Finset.mem_univ x)

omit [DecidableEq Ω] in
theorem exists_minEnergy [Nonempty Ω] : ∃ x, E x = minEnergy E := by
  obtain ⟨x, _, hx⟩ := Finset.exists_mem_eq_inf' Finset.univ_nonempty E
  exact ⟨x, hx.symm⟩

/-- The set of states off the exact minimum set. -/
noncomputable def offMin [Nonempty Ω] : Finset Ω :=
  Finset.univ.filter (fun x => E x ≠ minEnergy E)

omit [DecidableEq Ω] in
@[simp] theorem mem_offMin [Nonempty Ω] {x : Ω} :
    x ∈ offMin E ↔ E x ≠ minEnergy E := by
  unfold offMin
  simp

/-- Gibbs mass off the exact minimum set at inverse temperature
`beta`. -/
noncomputable def offMinMass [Nonempty Ω] (beta : ℝ) : ℝ :=
  excitedMass E beta (minEnergy E)

/-- The exact least positive energy excess over the minimum. -/
noncomputable def energyGap [Nonempty Ω] (h : (offMin E).Nonempty) : ℝ :=
  (offMin E).inf' h (fun x => E x - minEnergy E)

omit [DecidableEq Ω] in
theorem energyGap_pos [Nonempty Ω] (h : (offMin E).Nonempty) :
    0 < energyGap E h := by
  unfold energyGap
  rw [Finset.lt_inf'_iff]
  intro x hx
  have hne : E x ≠ minEnergy E := (mem_offMin E).mp hx
  have hlt : minEnergy E < E x :=
    lt_of_le_of_ne (minEnergy_le E x) (Ne.symm hne)
  linarith

omit [DecidableEq Ω] in
theorem energyGap_le [Nonempty Ω] (h : (offMin E).Nonempty) {x : Ω}
    (hx : x ∈ offMin E) : energyGap E h ≤ E x - minEnergy E :=
  Finset.inf'_le _ hx

omit [DecidableEq Ω] in
/-- The gap dichotomy: every state sits at the minimum or at least one
gap above it. -/
theorem minEnergy_add_gap_le [Nonempty Ω] (h : (offMin E).Nonempty)
    (x : Ω) :
    E x = minEnergy E ∨ minEnergy E + energyGap E h ≤ E x := by
  by_cases hx : E x = minEnergy E
  · exact Or.inl hx
  · right
    have hmem : x ∈ offMin E := (mem_offMin E).mpr hx
    have := energyGap_le E h hmem
    linarith

omit [Fintype Ω] [DecidableEq Ω] in
theorem gibbsWeight_pos (beta : ℝ) (x : Ω) :
    0 < gibbsWeight E beta x :=
  Real.exp_pos _

theorem offMinMass_nonneg [Nonempty Ω] (beta : ℝ) :
    0 ≤ offMinMass E beta := by
  unfold offMinMass excitedMass
  apply div_nonneg
  · exact Finset.sum_nonneg fun x _ => (gibbsWeight_pos E beta x).le
  · exact (partitionZ_pos E beta).le

theorem offMinMass_le_one [Nonempty Ω] (beta : ℝ) :
    offMinMass E beta ≤ 1 := by
  unfold offMinMass excitedMass
  rw [div_le_one (partitionZ_pos E beta)]
  unfold partitionZ
  apply Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
  exact fun x _ _ => (gibbsWeight_pos E beta x).le

omit [DecidableEq Ω] in
/-- Constant energy makes the off-minimum set empty and its mass zero:
the degenerate concentration case. -/
theorem offMinMass_eq_zero_of_constant [Nonempty Ω] (c beta : ℝ)
    (hc : ∀ x, E x = c) : offMinMass E beta = 0 := by
  have hmin : minEnergy E = c := by
    obtain ⟨x0, hx0⟩ := exists_minEnergy E
    rw [← hx0, hc x0]
  unfold offMinMass excitedMass
  have hempty : Finset.univ.filter (fun x => E x ≠ minEnergy E) = ∅ := by
    rw [Finset.eq_empty_iff_forall_notMem]
    intro x hx
    exact (Finset.mem_filter.mp hx).2 (by rw [hc x, hmin])
  rw [hempty, Finset.sum_empty, zero_div]

end MinimumSet

section ExplicitBound

variable (E : Ω → ℝ)

/-- **Cardinality-uniform low-temperature bound.** Any positive gap
below the exact gap and any cardinality bound above the exact
cardinality control the off-minimum mass at every real inverse
temperature: `offMinMass beta ≤ n * exp (-beta * gap)`. -/
theorem offMinMass_le_card_mul [Nonempty Ω] (beta gap : ℝ)
    (hgap : 0 < gap)
    (hdich : ∀ x, E x = minEnergy E ∨ minEnergy E + gap ≤ E x)
    (n : ℕ) (hn : Fintype.card Ω ≤ n) :
    offMinMass E beta ≤ (n : ℝ) * Real.exp (-beta * gap) := by
  have hcard1 : 1 ≤ n := le_trans Fintype.card_pos hn
  rcases le_or_gt 0 beta with hβ | hβ
  · obtain ⟨x0, hx0⟩ := exists_minEnergy E
    have hnonempty :
        (Finset.univ.filter (fun x => E x = minEnergy E)).Nonempty :=
      ⟨x0, Finset.mem_filter.mpr ⟨Finset.mem_univ x0, hx0⟩⟩
    have hb := excitedMass_le E beta (minEnergy E) gap hβ hgap hdich
      hnonempty
    refine le_trans hb ?_
    apply mul_le_mul_of_nonneg_right _ (Real.exp_pos (-beta * gap)).le
    have hGpos : (0 : ℝ) <
        ((Finset.univ.filter (fun x => E x = minEnergy E)).card : ℝ) := by
      exact_mod_cast Finset.card_pos.mpr hnonempty
    rw [div_le_iff₀ hGpos]
    have hXle :
        ((Finset.univ.filter (fun x => E x ≠ minEnergy E)).card : ℝ)
          ≤ (n : ℝ) := by
      have hle : (Finset.univ.filter
          (fun x => E x ≠ minEnergy E)).card ≤ Fintype.card Ω := by
        rw [← Finset.card_univ]
        exact Finset.card_le_card (Finset.filter_subset _ _)
      exact_mod_cast le_trans hle hn
    have hG1 : (1 : ℝ) ≤
        ((Finset.univ.filter (fun x => E x = minEnergy E)).card : ℝ) := by
      exact_mod_cast Finset.card_pos.mpr hnonempty
    nlinarith
  · have h1 := offMinMass_le_one E beta
    have hexp1 : (1 : ℝ) ≤ Real.exp (-beta * gap) := by
      have hx : (0 : ℝ) ≤ -beta * gap := by nlinarith
      have := Real.add_one_le_exp (-beta * gap)
      linarith
    have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hcard1
    nlinarith

/-- **Explicit finite low-temperature bound at every real inverse
temperature**, with the exact cardinality and the exact gap:
`offMinMass beta ≤ card Ω * exp (-beta * energyGap)`. -/
theorem offMinMass_le [Nonempty Ω] (h : (offMin E).Nonempty) (beta : ℝ) :
    offMinMass E beta
      ≤ (Fintype.card Ω : ℝ) * Real.exp (-beta * energyGap E h) :=
  offMinMass_le_card_mul E beta (energyGap E h) (energyGap_pos E h)
    (minEnergy_add_gap_le E h) (Fintype.card Ω) le_rfl

theorem tendsto_const_mul_exp_neg_gap (c gap : ℝ) (hgap : 0 < gap) :
    Filter.Tendsto (fun beta : ℝ => c * Real.exp (-beta * gap))
      Filter.atTop (nhds 0) := by
  have h1 : Filter.Tendsto (fun beta : ℝ => -beta * gap)
      Filter.atTop Filter.atBot := by
    have := (Filter.tendsto_neg_atTop_atBot (G := ℝ)).atBot_mul_const hgap
    simpa using this
  have h2 : Filter.Tendsto (fun beta : ℝ => Real.exp (-beta * gap))
      Filter.atTop (nhds 0) := Real.tendsto_exp_atBot.comp h1
  simpa using h2.const_mul c

/-- **Low-temperature concentration limit.** As `beta` tends to
infinity the off-minimum Gibbs mass tends to zero. -/
theorem offMinMass_tendsto_zero [Nonempty Ω] (h : (offMin E).Nonempty) :
    Filter.Tendsto (fun beta => offMinMass E beta)
      Filter.atTop (nhds 0) :=
  squeeze_zero (offMinMass_nonneg E) (offMinMass_le E h)
    (tendsto_const_mul_exp_neg_gap (Fintype.card Ω) (energyGap E h)
      (energyGap_pos E h))

end ExplicitBound

section RefinementUniform

/-- A declared refinement family: indexed finite energies with
minimum-preserving refinement maps, one uniform positive gap lower
bound, and one uniform cardinality bound. -/
structure UniformGapRefinement {I : Type*} [Preorder I]
    (X : I → Type*) [∀ r, Fintype (X r)] [∀ r, DecidableEq (X r)]
    [∀ r, Nonempty (X r)] (E : ∀ r, X r → ℝ) where
  /-- Coarse-graining map from a finer member to a coarser one. -/
  refineMap : ∀ {r s : I}, r ≤ s → X s → X r
  /-- Refinement maps carry minimizers to minimizers. -/
  refine_min : ∀ {r s : I} (hrs : r ≤ s) (x : X s),
    E s x = minEnergy (E s) → E r (refineMap hrs x) = minEnergy (E r)
  /-- The uniform gap lower bound. -/
  gapBound : ℝ
  gapBound_pos : 0 < gapBound
  /-- Every member sits at its minimum or at least `gapBound` above
  it. -/
  gap_dichotomy : ∀ (r : I) (x : X r),
    E r x = minEnergy (E r) ∨ minEnergy (E r) + gapBound ≤ E r x
  /-- The uniform cardinality bound. -/
  cardBound : ℕ
  card_le : ∀ r : I, Fintype.card (X r) ≤ cardBound

variable {I : Type*} [Preorder I] {X : I → Type*}
  [∀ r, Fintype (X r)] [∀ r, DecidableEq (X r)] [∀ r, Nonempty (X r)]
  {E : ∀ r, X r → ℝ}

/-- **Refinement-uniform low-temperature bound.** One explicit bound
controls every member of the declared family simultaneously at every
real inverse temperature. -/
theorem UniformGapRefinement.uniform_bound
    (F : UniformGapRefinement X E) (beta : ℝ) (r : I) :
    offMinMass (E r) beta
      ≤ (F.cardBound : ℝ) * Real.exp (-beta * F.gapBound) :=
  offMinMass_le_card_mul (E r) beta F.gapBound F.gapBound_pos
    (F.gap_dichotomy r) F.cardBound (F.card_le r)

/-- **Refinement-uniform concentration.** For every positive tolerance
there is one inverse-temperature threshold past which every member of
the family holds its off-minimum mass below the tolerance. -/
theorem UniformGapRefinement.uniform_concentration
    (F : UniformGapRefinement X E) :
    ∀ eps : ℝ, 0 < eps → ∃ beta0 : ℝ, ∀ beta : ℝ, beta0 ≤ beta →
      ∀ r : I, offMinMass (E r) beta < eps := by
  intro eps heps
  have htend := tendsto_const_mul_exp_neg_gap (F.cardBound : ℝ)
    F.gapBound F.gapBound_pos
  have hev := htend.eventually_lt_const heps
  obtain ⟨beta0, hb⟩ := Filter.eventually_atTop.mp hev
  exact ⟨beta0, fun beta hbeta r =>
    lt_of_le_of_lt (F.uniform_bound beta r) (hb beta hbeta)⟩

end RefinementUniform

section FourLawInstance

variable (π : Ω → ℝ) (b : Ω → B)

/-- Reference-relative Gibbs weight over the reference law `π`. -/
noncomputable def refGibbsWeight (E : Ω → ℝ) (beta : ℝ) (x : Ω) : ℝ :=
  π x * Real.exp (-beta * E x)

/-- Reference-relative partition function. -/
noncomputable def refPartitionZ (E : Ω → ℝ) (beta : ℝ) : ℝ :=
  ∑ x, refGibbsWeight π E beta x

/-- Normalized reference-relative Gibbs state. -/
noncomputable def refGibbs (E : Ω → ℝ) (beta : ℝ) (x : Ω) : ℝ :=
  refGibbsWeight π E beta x / refPartitionZ π E beta

omit [DecidableEq Ω] in
/-- **Degenerate instance receipt, state side.** At constant energy the
reference-relative Gibbs state equals the reference law itself at every
inverse temperature. -/
theorem refGibbs_const_energy (hπ1 : ∑ x, π x = 1) (c beta : ℝ) :
    refGibbs π (fun _ => c) beta = π := by
  funext x
  simp only [refGibbs, refGibbsWeight, refPartitionZ]
  rw [← Finset.sum_mul, hπ1, one_mul, mul_div_assoc,
    div_self (Real.exp_ne_zero _), mul_one]

/-- **Degenerate instance receipt, transition side.** The
conditional-resampling repair kernel of the four-law package fixes the
constant-energy Gibbs member at every inverse temperature, because that
member is the reference law. -/
theorem refGibbs_const_stationary (hπ : ∀ x, 0 < π x)
    (hπ1 : ∑ x, π x = 1) (c beta : ℝ) :
    push (refGibbs π (fun _ => c) beta) (heatBath π b)
      = refGibbs π (fun _ => c) beta := by
  rw [refGibbs_const_energy π hπ1 c beta]
  funext y
  exact heatBath_stationary hπ y

omit [DecidableEq Ω] in
/-- At constant energy the off-minimum mass vanishes identically: the
degenerate member concentrates at every inverse temperature. -/
theorem constEnergy_offMinMass [Nonempty Ω] (c beta : ℝ) :
    offMinMass (fun _ : Ω => c) beta = 0 :=
  offMinMass_eq_zero_of_constant _ c beta (fun _ => rfl)

end FourLawInstance

section TwoLevelWitness

/-- The nonconstant two-level witness energy: level `0` on `false` and
level `δ` on `true`. -/
noncomputable def twoLevelEnergy (δ : ℝ) : Bool → ℝ :=
  fun s => if s then δ else 0

theorem twoLevel_minEnergy (δ : ℝ) (hδ : 0 ≤ δ) :
    minEnergy (twoLevelEnergy δ) = 0 := by
  apply le_antisymm
  · simpa [twoLevelEnergy] using minEnergy_le (twoLevelEnergy δ) false
  · apply Finset.le_inf'
    intro x _
    cases x <;> simp [twoLevelEnergy, hδ]

theorem twoLevel_offMin (δ : ℝ) (hδ : 0 < δ) :
    offMin (twoLevelEnergy δ) = {true} := by
  ext x
  rw [mem_offMin, twoLevel_minEnergy δ hδ.le]
  cases x <;> simp [twoLevelEnergy, hδ.ne']

theorem twoLevel_offMin_nonempty (δ : ℝ) (hδ : 0 < δ) :
    (offMin (twoLevelEnergy δ)).Nonempty := by
  rw [twoLevel_offMin δ hδ]
  exact ⟨true, Finset.mem_singleton_self true⟩

/-- The two-level witness carries the exact gap `δ`. -/
theorem twoLevel_energyGap (δ : ℝ) (hδ : 0 < δ)
    (h : (offMin (twoLevelEnergy δ)).Nonempty) :
    energyGap (twoLevelEnergy δ) h = δ := by
  apply le_antisymm
  · have htrue : true ∈ offMin (twoLevelEnergy δ) := by
      rw [twoLevel_offMin δ hδ]
      exact Finset.mem_singleton_self true
    have := energyGap_le (twoLevelEnergy δ) h htrue
    simpa [twoLevelEnergy, twoLevel_minEnergy δ hδ.le] using this
  · apply Finset.le_inf'
    intro x hx
    have hxne := (mem_offMin (twoLevelEnergy δ)).mp hx
    cases x
    · exact absurd (by simp [twoLevelEnergy, twoLevel_minEnergy δ hδ.le])
        hxne
    · simp [twoLevelEnergy, twoLevel_minEnergy δ hδ.le]

/-- Exact off-minimum mass of the two-level witness. -/
theorem twoLevel_offMinMass (δ : ℝ) (hδ : 0 < δ) (beta : ℝ) :
    offMinMass (twoLevelEnergy δ) beta
      = Real.exp (-beta * δ) / (1 + Real.exp (-beta * δ)) := by
  unfold offMinMass excitedMass partitionZ
  rw [twoLevel_minEnergy δ hδ.le]
  have hfilter :
      Finset.univ.filter (fun x : Bool => twoLevelEnergy δ x ≠ 0)
        = {true} := by
    ext x
    cases x <;> simp [twoLevelEnergy, hδ.ne']
  rw [hfilter, Finset.sum_singleton]
  have hZ : ∑ x : Bool, gibbsWeight (twoLevelEnergy δ) beta x
      = Real.exp (-beta * δ) + 1 := by
    rw [Fintype.sum_bool]
    unfold gibbsWeight twoLevelEnergy
    simp
  rw [hZ]
  unfold gibbsWeight twoLevelEnergy
  rw [add_comm]
  simp

/-- Explicit-bound instance on the two-level witness: off-minimum mass
at most `2 * exp (-beta * δ)` at every real inverse temperature. -/
theorem twoLevel_bound (δ : ℝ) (hδ : 0 < δ) (beta : ℝ) :
    offMinMass (twoLevelEnergy δ) beta ≤ 2 * Real.exp (-beta * δ) := by
  have h := offMinMass_le (twoLevelEnergy δ)
    (twoLevel_offMin_nonempty δ hδ) beta
  rw [twoLevel_energyGap δ hδ (twoLevel_offMin_nonempty δ hδ)] at h
  simpa [Fintype.card_bool] using h

/-- The constant family over the two-level witness inhabits the
refinement structure with exact bounds `cardBound = 2` and
`gapBound = δ`. -/
noncomputable def twoLevelTower (δ : ℝ) (hδ : 0 < δ) :
    UniformGapRefinement (fun _ : ℕ => Bool)
      (fun _ => twoLevelEnergy δ) where
  refineMap _ x := x
  refine_min _ _ hx := hx
  gapBound := δ
  gapBound_pos := hδ
  gap_dichotomy _ x := by
    cases x
    · left
      simp [twoLevelEnergy, twoLevel_minEnergy δ hδ.le]
    · right
      simp [twoLevelEnergy, twoLevel_minEnergy δ hδ.le]
  cardBound := 2
  card_le _ := by simp [Fintype.card_bool]

end TwoLevelWitness

section NegativeControl

/-- Off-minimum mass lower bound on the two-level witness in the
small-`beta * δ` regime: below the threshold `beta * δ ≤ 1/2` at least
one quarter of the mass sits off the minimum. -/
theorem twoLevel_offMinMass_ge (δ : ℝ) (hδ : 0 < δ) (beta : ℝ)
    (hsmall : beta * δ ≤ 1 / 2) :
    (1 : ℝ) / 4 ≤ offMinMass (twoLevelEnergy δ) beta := by
  rw [twoLevel_offMinMass δ hδ beta]
  have hepos : 0 < Real.exp (-beta * δ) := Real.exp_pos _
  have hege : (1 : ℝ) / 2 ≤ Real.exp (-beta * δ) := by
    have h1 := Real.add_one_le_exp (-beta * δ)
    nlinarith
  rw [le_div_iff₀ (by positivity)]
  nlinarith

/-- The shrinking-gap family: two-level members with gaps
`1 / (r + 1)` tending to zero. -/
noncomputable def shrinkingGapEnergy (r : ℕ) : Bool → ℝ :=
  twoLevelEnergy (1 / ((r : ℝ) + 1))

theorem shrinkingGap_delta_pos (r : ℕ) :
    (0 : ℝ) < 1 / ((r : ℝ) + 1) := by positivity

/-- Member `r` of the shrinking-gap family carries the exact gap
`1 / (r + 1)`. -/
theorem shrinkingGap_energyGap (r : ℕ)
    (h : (offMin (shrinkingGapEnergy r)).Nonempty) :
    energyGap (shrinkingGapEnergy r) h = 1 / ((r : ℝ) + 1) :=
  twoLevel_energyGap _ (shrinkingGap_delta_pos r) h

/-- The shrinking-gap family admits no uniform positive gap lower
bound. -/
theorem shrinkingGap_no_uniform_gap :
    ∀ g : ℝ, 0 < g → ∃ (r : ℕ)
      (h : (offMin (shrinkingGapEnergy r)).Nonempty),
      energyGap (shrinkingGapEnergy r) h < g := by
  intro g hg
  obtain ⟨n, hn⟩ := exists_nat_gt (1 / g)
  refine ⟨n, twoLevel_offMin_nonempty _ (shrinkingGap_delta_pos n), ?_⟩
  rw [shrinkingGap_energyGap]
  rw [div_lt_iff₀ (by positivity : (0 : ℝ) < (n : ℝ) + 1)]
  rw [div_lt_iff₀ hg] at hn
  nlinarith

/-- At every inverse temperature some member of the shrinking-gap
family keeps at least one quarter of its mass off the minimum. -/
theorem shrinkingGap_mass_lower (beta : ℝ) :
    ∃ r : ℕ, (1 : ℝ) / 4 ≤ offMinMass (shrinkingGapEnergy r) beta := by
  obtain ⟨n, hn⟩ := exists_nat_ge (2 * beta)
  refine ⟨n, ?_⟩
  apply twoLevel_offMinMass_ge _ (shrinkingGap_delta_pos n)
  rw [mul_one_div, div_le_iff₀ (by positivity : (0 : ℝ) < (n : ℝ) + 1)]
  linarith

/-- The declared failure point: at `beta = 1`, member `r = 2` of the
shrinking-gap family keeps at least one quarter of its mass off the
minimum. -/
theorem shrinkingGap_failure_at_beta_one :
    (1 : ℝ) / 4 ≤ offMinMass (shrinkingGapEnergy 2) 1 := by
  apply twoLevel_offMinMass_ge _ (shrinkingGap_delta_pos 2)
  norm_num

/-- **Negative control.** The shrinking-gap family violates the uniform
concentration conclusion of `UniformGapRefinement`: no
inverse-temperature threshold controls every member below the
tolerance `1/4`. -/
theorem shrinkingGap_uniform_concentration_fails :
    ¬ ∀ eps : ℝ, 0 < eps → ∃ beta0 : ℝ, ∀ beta : ℝ, beta0 ≤ beta →
      ∀ r : ℕ, offMinMass (shrinkingGapEnergy r) beta < eps := by
  intro hcon
  obtain ⟨beta0, hb⟩ := hcon (1 / 4) (by norm_num)
  obtain ⟨r, hr⟩ := shrinkingGap_mass_lower (max beta0 0)
  exact absurd (hb (max beta0 0) (le_max_left _ _) r) (not_lt.mpr hr)

end NegativeControl

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.offMinMass_le_card_mul
#print axioms OPH.Thermodynamics.offMinMass_le
#print axioms OPH.Thermodynamics.offMinMass_tendsto_zero
#print axioms OPH.Thermodynamics.UniformGapRefinement.uniform_bound
#print axioms OPH.Thermodynamics.UniformGapRefinement.uniform_concentration
#print axioms OPH.Thermodynamics.refGibbs_const_energy
#print axioms OPH.Thermodynamics.refGibbs_const_stationary
#print axioms OPH.Thermodynamics.twoLevel_energyGap
#print axioms OPH.Thermodynamics.twoLevel_offMinMass
#print axioms OPH.Thermodynamics.twoLevel_bound
#print axioms OPH.Thermodynamics.shrinkingGap_no_uniform_gap
#print axioms OPH.Thermodynamics.shrinkingGap_uniform_concentration_fails
