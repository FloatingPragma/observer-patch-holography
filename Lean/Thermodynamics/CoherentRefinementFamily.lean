import LowTemperatureControl
import FourLawAdequacySurface

set_option autoImplicit false

namespace OPH.Thermodynamics

/-!
# Coherent directed refinement families with a uniform spectral-tail envelope

Targeted PR-08 work (V3 issue #739). This module is finite mathematics
throughout; every theorem is stated on explicit finite state spaces and
uses the exact Gibbs-mass definitions of `LowTemperatureControl`
(`offMinMass`, `minEnergy`, `gibbsWeight`, `partitionZ`).  Its companion
result for the four-law surface is described precisely below.

What IS proved here:

* `CoherentRefinementFamily`: an indexed family of finite energies over
  a preorder with a directedness field, coarse-graining maps satisfying
  the identity and composition coherence laws, an energy-compatibility
  law used by `minEnergy_mono`,
  a designated ground datum carried to the stage minimum and preserved
  by the maps, and one uniform spectral-tail envelope: a threshold
  `beta0`, a function `Env` bounding the off-minimum Gibbs mass of
  every stage simultaneously beyond the threshold, with `Env` tending
  to zero at infinite `beta`;
* `uniform_concentration`: for every positive tolerance one
  inverse-temperature threshold takes the off-minimum mass of every
  stage below the tolerance; the companion theorem exports this beside
  the old four-law conclusion record;
* exact field use: `minEnergy_mono` uses energy compatibility,
  `refine_ground_min` uses the ground-minimum and ground-preservation
  fields, and `refine_energy_le_chain` uses composition and energy
  compatibility.  The current concentration and four-law companion do
  not use directedness, identity coherence, or those auxiliary lemmas;
* retention: `ofUniformGap` constructs a coherent family from
  `UniformGapRefinement` data only when directedness, coherence, energy
  compatibility, and coherent ground data are additionally supplied.
  `ofUniformGapConst` uses the old gap and cardinality bounds but installs
  identity maps rather than preserving arbitrary old refine maps.
  `twoLevelCoherent` retains the committed two-level gap witness
  `twoLevelTower` as a degenerate constant-carrier instance;
* cardinality separation: the ladder family `ladderFamily` has stage `n`
  carrying energy levels `0, δ, 2δ, …, nδ` with level `k` of
  degeneracy `2^k`, nontrivial stage-truncation refine maps satisfying
  the coherence laws, and the geometric spectral-tail envelope
  `r / (1 - r)` with `r = 2 · exp (-beta · δ)` beyond the explicit
  threshold `beta0 = log 4 / δ` (`geom_tail_sum_le`,
  `geom_tail_tsum`, `ladder_offMinMass_le`); its stage cardinality is
  unbounded (`ladder_card_unbounded`), so the uniform-cardinality
  hypothesis of `UniformGapRefinement` fails on it:
  `ladder_uniformGap_isEmpty` proves `UniformGapRefinement` is
  uninhabited on the ladder data while `ladderFamily` inhabits the
  coherent structure.  This separates the spectral-envelope condition
  from uniform total cardinality on the same carriers and energies; it
  is not a proof that either complete structure contains the other;
* `constUniformGap`: every finite energy carries the
  gap-plus-cardinality data on a one-member constant family (exact gap
  when the off-minimum set is inhabited, vacuous gap otherwise);
* the companion hook on the four-law surface:
  `CoherentRefinementAttachment` ties the coherent family to the
  calibrated energy exactly as `RefinementAttachment` ties the
  uniform-gap form, `FourLawCoherentAntecedent` bundles rows PR-07 and
  PR-15 unchanged with row PR-08 in the narrowed form, and
  `fourLaws_composed_coherent` conjoins two results.  First it produces
  `FourLawConclusions` through `toAntecedent`, which carries rows PR-07
  and PR-15 over definitionally but fills the old PR-08 field with a
  newly constructed degenerate constant family.  Separately it proves
  the supplied directed family's envelope, concentration, and
  calibrated-member clauses.  It does not substitute that supplied
  family into the old `FourLawConclusions` type.
  The hook lives in this module rather than in
  `FourLawAdequacySurface`, keeping every committed declaration of the
  surface untouched; this module is registered as an
  `OPHThermodynamics` lakefile root, so the hook is elaborated by
  `lake build OPHThermodynamics` and by CI.

What is NOT proved here:

* cofinality: `directed` says only that two indices have a common upper
  bound.  There is no ambient regulator preorder, stage map into such a
  preorder, cofinality predicate, strict-extension/no-max condition, or
  continuum limit.  Constant and singleton index systems are allowed;
* inclusion between the complete `UniformGapRefinement` and
  `CoherentRefinementFamily` structures.  The former lacks several
  fields required by the latter, while the latter need not have a
  uniform cardinality bound.  The ladder proves only the stated
  same-data cardinality separation;
* replacement of the old four-law refinement member by the supplied
  directed family.  The theorem `fourLaws_composed_coherent` is the
  companion conjunction described above;
* the ground-degeneracy entropy limit (Gibbs entropy tending to
  `log g0`): no committed normalization or Gibbs-entropy limit lemma
  exists in the thermodynamics modules, and this file adds none;
* any identification of `beta` with a physical inverse temperature or
  of the energies with a realized objective: register rows PR-15 and
  PR-07 stay open exactly as recorded in the four-law surface;
* any claim that a realized source artifact supplies a coherent directed
  refinement family: register row PR-08 remains a registered premise.
-/

/-! ## Geometric tail bounds

The finite partial sums of the shifted geometric series are bounded by
`r / (1 - r)`, and the full series sums to that value. These are the
envelope receipts of the ladder family. -/

/-- Finite geometric tail bound: for `0 ≤ r < 1` every partial sum
`∑_{j < n} r^(j+1)` is at most `r / (1 - r)`. -/
theorem geom_tail_sum_le (radius : ℝ) (h0 : 0 ≤ radius)
    (h1 : radius < 1) (n : ℕ) :
    ∑ j ∈ Finset.range n, radius ^ (j + 1) ≤ radius / (1 - radius) := by
  have h1r : (0 : ℝ) < 1 - radius := by linarith
  have hS : ∑ j ∈ Finset.range n, radius ^ (j + 1)
      = radius * ∑ j ∈ Finset.range n, radius ^ j := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun j _ => by ring
  have hgeom : (∑ j ∈ Finset.range n, radius ^ j) * (1 - radius)
      = 1 - radius ^ n := by
    have h := geom_sum_mul radius n
    have h2 : (∑ j ∈ Finset.range n, radius ^ j) * (1 - radius)
        = -((∑ j ∈ Finset.range n, radius ^ j) * (radius - 1)) := by
      ring
    rw [h2, h]
    ring
  rw [hS, le_div_iff₀ h1r]
  calc radius * (∑ j ∈ Finset.range n, radius ^ j) * (1 - radius)
      = radius * ((∑ j ∈ Finset.range n, radius ^ j) * (1 - radius)) := by
        ring
    _ = radius * (1 - radius ^ n) := by rw [hgeom]
    _ ≤ radius * 1 := by
        have hpow : (0 : ℝ) ≤ radius ^ n := pow_nonneg h0 n
        exact mul_le_mul_of_nonneg_left (by linarith) h0
    _ = radius := mul_one radius

/-- Exact geometric tail sum: for `0 ≤ r < 1`,
`∑_{j ≥ 0} r^(j+1) = r / (1 - r)`. The finite stage sums of the ladder
family are partial sums of this series. -/
theorem geom_tail_tsum (radius : ℝ) (h0 : 0 ≤ radius)
    (h1 : radius < 1) :
    ∑' j : ℕ, radius ^ (j + 1) = radius / (1 - radius) := by
  calc ∑' j : ℕ, radius ^ (j + 1)
      = ∑' j : ℕ, radius * radius ^ j := tsum_congr fun j => by ring
    _ = radius * ∑' j : ℕ, radius ^ j := tsum_mul_left
    _ = radius * (1 - radius)⁻¹ := by
        rw [tsum_geometric_of_lt_one h0 h1]
    _ = radius / (1 - radius) := (div_eq_mul_inv _ _).symm

/-! ## The coherent directed refinement family -/

/-- **A coherent directed refinement family with a uniform spectral-tail
envelope.** An indexed family of finite energies over a preorder with:
a directedness field; coarse-graining maps satisfying the identity and
composition coherence laws; energy compatibility, which
`minEnergy_mono` uses; a
designated ground datum realizing each stage minimum and preserved by
the maps; and one envelope `Env` with threshold `beta0` bounding the
off-minimum Gibbs mass of every stage simultaneously, with `Env`
tending to zero at infinite `beta`.  `ofUniformGap` constructs this
structure from gap-plus-cardinality data only after the additional
directedness, coherence, energy-compatibility, and ground hypotheses are
supplied.  The ladder shows that this structure's envelope need not come
with a uniform total-cardinality bound; it does not establish inclusion
between the two complete structures. -/
structure CoherentRefinementFamily {I : Type*} [Preorder I]
    (X : I → Type*) [∀ r, Fintype (X r)] [∀ r, DecidableEq (X r)]
    [∀ r, Nonempty (X r)] (E : ∀ r, X r → ℝ) where
  /-- Directedness: every two stages admit a common refinement. -/
  directed : ∀ r s : I, ∃ t : I, r ≤ t ∧ s ≤ t
  /-- Coarse-graining map from a finer stage to a coarser stage. -/
  refineMap : ∀ {r s : I}, r ≤ s → X s → X r
  /-- Identity coherence law. -/
  refine_id : ∀ (r : I) (x : X r), refineMap (le_refl r) x = x
  /-- Composition coherence law. -/
  refine_comp : ∀ {r s t : I} (hrs : r ≤ s) (hst : s ≤ t) (x : X t),
    refineMap hrs (refineMap hst x) = refineMap (le_trans hrs hst) x
  /-- Energy compatibility: coarse-graining does not raise the energy.
  The auxiliary theorem `minEnergy_mono` consumes this field. -/
  refine_energy_le : ∀ {r s : I} (hrs : r ≤ s) (x : X s),
    E r (refineMap hrs x) ≤ E s x
  /-- The designated ground datum at each stage. -/
  ground : ∀ r : I, X r
  /-- The ground datum realizes the exact stage minimum. -/
  ground_min : ∀ r : I, E r (ground r) = minEnergy (E r)
  /-- Coarse-graining carries the ground datum to the ground datum. -/
  refine_ground : ∀ {r s : I} (hrs : r ≤ s),
    refineMap hrs (ground s) = ground r
  /-- The envelope threshold. -/
  beta0 : ℝ
  /-- The uniform spectral-tail envelope. -/
  Env : ℝ → ℝ
  /-- Beyond the threshold, one envelope value bounds the off-minimum
  Gibbs mass of every stage simultaneously. -/
  env_bound : ∀ beta : ℝ, beta0 ≤ beta → ∀ r : I,
    offMinMass (E r) beta ≤ Env beta
  /-- The envelope tends to zero at infinite inverse temperature. -/
  env_tendsto : Filter.Tendsto Env Filter.atTop (nhds 0)

namespace CoherentRefinementFamily

variable {I : Type*} [Preorder I] {X : I → Type*}
  [∀ r, Fintype (X r)] [∀ r, DecidableEq (X r)] [∀ r, Nonempty (X r)]
  {E : ∀ r, X r → ℝ}

/-- Stage minima are monotone along refinement. This theorem consumes
the energy-compatibility law through the refine maps. -/
theorem minEnergy_mono (F : CoherentRefinementFamily X E) {r s : I}
    (hrs : r ≤ s) : minEnergy (E r) ≤ minEnergy (E s) := by
  obtain ⟨x, hx⟩ := exists_minEnergy (E s)
  calc minEnergy (E r) ≤ E r (F.refineMap hrs x) :=
        minEnergy_le (E r) (F.refineMap hrs x)
    _ ≤ E s x := F.refine_energy_le hrs x
    _ = minEnergy (E s) := hx

/-- The refine maps carry the ground datum of a finer stage onto the
exact minimum of the coarser stage. Consumes `refine_ground` and
`ground_min`. -/
theorem refine_ground_min (F : CoherentRefinementFamily X E) {r s : I}
    (hrs : r ≤ s) :
    E r (F.refineMap hrs (F.ground s)) = minEnergy (E r) := by
  rw [F.refine_ground hrs, F.ground_min r]

/-- Two-step coarse-graining lowers energy stagewise. Consumes the
composition coherence law. -/
theorem refine_energy_le_chain (F : CoherentRefinementFamily X E)
    {r s t : I} (hrs : r ≤ s) (hst : s ≤ t) (x : X t) :
    E r (F.refineMap (le_trans hrs hst) x)
      ≤ E s (F.refineMap hst x) := by
  rw [← F.refine_comp hrs hst x]
  exact F.refine_energy_le hrs (F.refineMap hst x)

/-- Envelope restatement with the threshold hypothesis explicit. -/
theorem envelope_bound (F : CoherentRefinementFamily X E) {beta : ℝ}
    (hbeta : F.beta0 ≤ beta) (r : I) :
    offMinMass (E r) beta ≤ F.Env beta :=
  F.env_bound beta hbeta r

/-- **Family-uniform concentration.** For every positive tolerance
there is one inverse-temperature threshold past which every stage of
the family holds its off-minimum Gibbs mass below the tolerance.  The
four-law companion theorem exports this conclusion beside, rather than
inside, the old gap-plus-cardinality conclusion record. -/
theorem uniform_concentration (F : CoherentRefinementFamily X E) :
    ∀ eps : ℝ, 0 < eps → ∃ beta1 : ℝ, ∀ beta : ℝ, beta1 ≤ beta →
      ∀ r : I, offMinMass (E r) beta < eps := by
  intro eps heps
  have hev := F.env_tendsto.eventually_lt_const heps
  obtain ⟨b, hb⟩ := Filter.eventually_atTop.mp hev
  refine ⟨max F.beta0 b, fun beta hbeta r => ?_⟩
  exact lt_of_le_of_lt
    (F.env_bound beta (le_trans (le_max_left _ _) hbeta) r)
    (hb beta (le_trans (le_max_right _ _) hbeta))

/-! ## Retention of the committed gap-plus-cardinality data -/

/-- **Retention, general form.** The committed `UniformGapRefinement`
data yields a coherent family whenever its maps additionally satisfy
directedness, the two coherence laws, energy compatibility, and ground
preservation — exactly the content the committed structure does not
carry. The envelope is the committed bound
`cardBound * exp (-beta * gapBound)` with threshold `0`. -/
noncomputable def ofUniformGap (F : UniformGapRefinement X E)
    (hdir : ∀ r s : I, ∃ t : I, r ≤ t ∧ s ≤ t)
    (hid : ∀ (r : I) (x : X r), F.refineMap (le_refl r) x = x)
    (hcomp : ∀ {r s t : I} (hrs : r ≤ s) (hst : s ≤ t) (x : X t),
      F.refineMap hrs (F.refineMap hst x)
        = F.refineMap (le_trans hrs hst) x)
    (hle : ∀ {r s : I} (hrs : r ≤ s) (x : X s),
      E r (F.refineMap hrs x) ≤ E s x)
    (g : ∀ r : I, X r)
    (hgmin : ∀ r : I, E r (g r) = minEnergy (E r))
    (hg : ∀ {r s : I} (hrs : r ≤ s), F.refineMap hrs (g s) = g r) :
    CoherentRefinementFamily X E where
  directed := hdir
  refineMap := F.refineMap
  refine_id := hid
  refine_comp := hcomp
  refine_energy_le := hle
  ground := g
  ground_min := hgmin
  refine_ground := hg
  beta0 := 0
  Env := fun beta => (F.cardBound : ℝ) * Real.exp (-beta * F.gapBound)
  env_bound := fun beta _ r => F.uniform_bound beta r
  env_tendsto :=
    tendsto_const_mul_exp_neg_gap (F.cardBound : ℝ) F.gapBound
      F.gapBound_pos

/-- **Retention, constant-carrier form.** The gap and cardinality bounds
of a `UniformGapRefinement` on a constant family over a directed preorder
yield a coherent family with newly installed identity refine maps and a
ground datum chosen from `exists_minEnergy`.  This is a degenerate
constant-carrier witness; it does not preserve arbitrary refine maps from
the supplied `UniformGapRefinement`. -/
noncomputable def ofUniformGapConst {Ω : Type*} [Fintype Ω]
    [DecidableEq Ω] [Nonempty Ω] {E0 : Ω → ℝ}
    (hdir : ∀ r s : I, ∃ t : I, r ≤ t ∧ s ≤ t)
    (F : UniformGapRefinement (fun _ : I => Ω) (fun _ => E0)) :
    CoherentRefinementFamily (fun _ : I => Ω) (fun _ => E0) where
  directed := hdir
  refineMap := fun _ x => x
  refine_id := fun _ _ => rfl
  refine_comp := fun _ _ _ => rfl
  refine_energy_le := fun _ _ => le_rfl
  ground := fun _ => (exists_minEnergy E0).choose
  ground_min := fun _ => (exists_minEnergy E0).choose_spec
  refine_ground := fun _ => rfl
  beta0 := 0
  Env := fun beta => (F.cardBound : ℝ) * Real.exp (-beta * F.gapBound)
  env_bound := fun beta _ r => F.uniform_bound beta r
  env_tendsto :=
    tendsto_const_mul_exp_neg_gap (F.cardBound : ℝ) F.gapBound
      F.gapBound_pos

end CoherentRefinementFamily

/-- **Retention of the committed gap witness.** The committed constant
tower `twoLevelTower` over the two-level energy inhabits the coherent
structure: identity refine maps, envelope `2 * exp (-beta * δ)`,
threshold `0`. -/
noncomputable def twoLevelCoherent (δ : ℝ) (hδ : 0 < δ) :
    CoherentRefinementFamily (fun _ : ℕ => Bool)
      (fun _ => twoLevelEnergy δ) :=
  CoherentRefinementFamily.ofUniformGapConst
    (fun r s => ⟨max r s, le_max_left r s, le_max_right r s⟩)
    (twoLevelTower δ hδ)

/-! ## The one-member constant family on any finite energy

Route used by the four-law surface to carry a coherent antecedent
through the committed composed statement: any finite energy carries
the gap-plus-cardinality data on the one-member constant family, with
the exact gap when the off-minimum set is inhabited and a vacuous gap
otherwise. -/

/-- Every finite energy admits the committed gap-plus-cardinality data
on a constant family: exact gap `energyGap` when the off-minimum set
is inhabited, gap `1` vacuously otherwise, cardinality bound the exact
cardinality. -/
noncomputable def constUniformGap (I : Type*) [Preorder I]
    (Ω : Type*) [Fintype Ω] [DecidableEq Ω] [Nonempty Ω] (E : Ω → ℝ) :
    UniformGapRefinement (fun _ : I => Ω) (fun _ => E) where
  refineMap := fun _ x => x
  refine_min := fun _ _ hx => hx
  gapBound := if h : (offMin E).Nonempty then energyGap E h else 1
  gapBound_pos := by
    split
    · exact energyGap_pos E ‹_›
    · exact one_pos
  gap_dichotomy := by
    intro r x
    split
    · exact minEnergy_add_gap_le E ‹_› x
    · rename_i h
      refine Or.inl ?_
      by_contra hne
      exact h ⟨x, (mem_offMin E).mpr hne⟩
  cardBound := Fintype.card Ω
  card_le := fun _ => le_rfl

/-! ## Cardinality separation: the exponential ladder family

Stage `n` carries the energy levels `0, δ, 2δ, …, nδ` with level `k`
of degeneracy `2^k`. The stage cardinality `2^(n+1) - 1` is unbounded,
so the committed uniform-cardinality hypothesis fails on this family
(`ladder_uniformGap_isEmpty`), while the geometric spectral tail
`r / (1 - r)`, `r = 2 * exp (-beta * δ)`, is a uniform envelope beyond
the explicit threshold `beta0 = log 4 / δ`. -/

/-- Stage `n` of the ladder: level `k ≤ n` with degeneracy `2^k`. -/
abbrev LadderState (n : ℕ) : Type :=
  Σ k : Fin (n + 1), Fin (2 ^ (k : ℕ))

/-- The designated ground datum of stage `n`: level `0`, copy `0`. -/
def ladderGround (n : ℕ) : LadderState n :=
  ⟨⟨0, Nat.succ_pos n⟩, ⟨0, pow_pos (by norm_num) _⟩⟩

instance ladderState_nonempty (n : ℕ) : Nonempty (LadderState n) :=
  ⟨ladderGround n⟩

/-- The ladder energy: level `k` sits at energy `k * δ`. -/
noncomputable def ladderEnergy (δ : ℝ) (n : ℕ) : LadderState n → ℝ :=
  fun x => ((x.1 : ℕ) : ℝ) * δ

theorem ladder_ground_energy (δ : ℝ) (n : ℕ) :
    ladderEnergy δ n (ladderGround n) = 0 := by
  show ((0 : ℕ) : ℝ) * δ = 0
  simp

theorem ladder_minEnergy (δ : ℝ) (hδ : 0 ≤ δ) (n : ℕ) :
    minEnergy (ladderEnergy δ n) = 0 := by
  apply le_antisymm
  · have h := minEnergy_le (ladderEnergy δ n) (ladderGround n)
    rwa [ladder_ground_energy] at h
  · apply Finset.le_inf'
    intro x _
    exact mul_nonneg (Nat.cast_nonneg _) hδ

/-- Stage truncation: keep a state whose level fits in the coarser
stage, send everything above to the ground datum. -/
def ladderTrunc (r : ℕ) {s : ℕ} (x : LadderState s) : LadderState r :=
  if h : (x.1 : ℕ) ≤ r then ⟨⟨(x.1 : ℕ), Nat.lt_succ_of_le h⟩, x.2⟩
  else ladderGround r

theorem ladderTrunc_id (r : ℕ) (x : LadderState r) :
    ladderTrunc r x = x := by
  unfold ladderTrunc
  exact dif_pos (Nat.lt_succ_iff.mp x.1.isLt)

theorem ladderTrunc_ground (r s : ℕ) :
    ladderTrunc r (ladderGround s) = ladderGround r := by
  unfold ladderTrunc
  exact dif_pos (Nat.zero_le r)

theorem ladderTrunc_comp {r s t : ℕ} (hrs : r ≤ s)
    (x : LadderState t) :
    ladderTrunc r (ladderTrunc s x) = ladderTrunc r x := by
  by_cases hs : (x.1 : ℕ) ≤ s
  · have h1 : ladderTrunc s x
        = (⟨⟨(x.1 : ℕ), Nat.lt_succ_of_le hs⟩, x.2⟩ : LadderState s) :=
      dif_pos hs
    rw [h1]
    by_cases hr : (x.1 : ℕ) ≤ r
    · have h2 : ladderTrunc r
          (⟨⟨(x.1 : ℕ), Nat.lt_succ_of_le hs⟩, x.2⟩ : LadderState s)
          = (⟨⟨(x.1 : ℕ), Nat.lt_succ_of_le hr⟩, x.2⟩ : LadderState r) :=
        dif_pos hr
      have h3 : ladderTrunc r x
          = (⟨⟨(x.1 : ℕ), Nat.lt_succ_of_le hr⟩, x.2⟩ : LadderState r) :=
        dif_pos hr
      rw [h2, h3]
    · have h2 : ladderTrunc r
          (⟨⟨(x.1 : ℕ), Nat.lt_succ_of_le hs⟩, x.2⟩ : LadderState s)
          = ladderGround r := dif_neg hr
      have h3 : ladderTrunc r x = ladderGround r := dif_neg hr
      rw [h2, h3]
  · have hr : ¬ (x.1 : ℕ) ≤ r := fun h => hs (le_trans h hrs)
    have h1 : ladderTrunc s x = ladderGround s := dif_neg hs
    have h3 : ladderTrunc r x = ladderGround r := dif_neg hr
    rw [h1, ladderTrunc_ground, h3]

theorem ladderTrunc_energy_le (δ : ℝ) (hδ : 0 ≤ δ) (r : ℕ) {s : ℕ}
    (x : LadderState s) :
    ladderEnergy δ r (ladderTrunc r x) ≤ ladderEnergy δ s x := by
  by_cases h : (x.1 : ℕ) ≤ r
  · have h1 : ladderTrunc r x
        = (⟨⟨(x.1 : ℕ), Nat.lt_succ_of_le h⟩, x.2⟩ : LadderState r) :=
      dif_pos h
    rw [h1]
    exact le_of_eq rfl
  · have h1 : ladderTrunc r x = ladderGround r := dif_neg h
    rw [h1, ladder_ground_energy]
    exact mul_nonneg (Nat.cast_nonneg _) hδ

theorem ladderGround_min (δ : ℝ) (hδ : 0 ≤ δ) (n : ℕ) :
    ladderEnergy δ n (ladderGround n) = minEnergy (ladderEnergy δ n) := by
  rw [ladder_ground_energy, ladder_minEnergy δ hδ n]

theorem ladder_partitionZ_ge_one (δ : ℝ) (n : ℕ) (beta : ℝ) :
    1 ≤ partitionZ (ladderEnergy δ n) beta := by
  have hg : gibbsWeight (ladderEnergy δ n) beta (ladderGround n) = 1 := by
    unfold gibbsWeight
    rw [ladder_ground_energy, mul_zero, Real.exp_zero]
  calc (1 : ℝ)
      = gibbsWeight (ladderEnergy δ n) beta (ladderGround n) := hg.symm
    _ ≤ ∑ x, gibbsWeight (ladderEnergy δ n) beta x :=
        Finset.single_le_sum
          (fun x _ => (gibbsWeight_pos (ladderEnergy δ n) beta x).le)
          (Finset.mem_univ (ladderGround n))
    _ = partitionZ (ladderEnergy δ n) beta := rfl

/-- The off-minimum weight sum of ladder stage `n` is exactly the
finite geometric tail `∑_{j < n} r^(j+1)` with
`r = 2 * exp (-beta * δ)`: level `k ≥ 1` contributes
`2^k * exp (-beta * k * δ) = r^k` exactly. -/
theorem ladder_excitedSum_eq (δ : ℝ) (hδ : 0 < δ) (n : ℕ) (beta : ℝ) :
    (∑ x ∈ Finset.univ.filter
        (fun x : LadderState n => ladderEnergy δ n x ≠ 0),
      gibbsWeight (ladderEnergy δ n) beta x)
      = ∑ j ∈ Finset.range n, (2 * Real.exp (-beta * δ)) ^ (j + 1) := by
  have hterm : ∀ k : Fin (n + 1),
      (∑ i : Fin (2 ^ (k : ℕ)),
        if ladderEnergy δ n ⟨k, i⟩ ≠ 0
        then gibbsWeight (ladderEnergy δ n) beta ⟨k, i⟩ else 0)
      = (fun m : ℕ =>
          if m ≠ 0 then (2 * Real.exp (-beta * δ)) ^ m else 0)
          (k : ℕ) := by
    intro k
    show (∑ i : Fin (2 ^ (k : ℕ)),
        if ladderEnergy δ n ⟨k, i⟩ ≠ 0
        then gibbsWeight (ladderEnergy δ n) beta ⟨k, i⟩ else 0)
      = if (k : ℕ) ≠ 0 then (2 * Real.exp (-beta * δ)) ^ (k : ℕ) else 0
    by_cases hk : (k : ℕ) = 0
    · rw [if_neg (fun h => h hk)]
      apply Finset.sum_eq_zero
      intro i _
      apply if_neg
      intro hcon
      apply hcon
      show (((k : ℕ) : ℝ)) * δ = 0
      rw [hk]
      simp
    · rw [if_pos hk]
      have hall : ∀ i ∈ (Finset.univ : Finset (Fin (2 ^ (k : ℕ)))),
          (if ladderEnergy δ n ⟨k, i⟩ ≠ 0
           then gibbsWeight (ladderEnergy δ n) beta ⟨k, i⟩ else 0)
          = Real.exp (-beta * (((k : ℕ) : ℝ) * δ)) := by
        intro i _
        rw [if_pos (show ladderEnergy δ n ⟨k, i⟩ ≠ 0 from
          mul_ne_zero (Nat.cast_ne_zero.mpr hk) hδ.ne')]
        rfl
      rw [Finset.sum_congr rfl hall, Finset.sum_const,
        Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_pow,
        ← Real.exp_nat_mul,
        show -beta * (((k : ℕ) : ℝ) * δ)
          = ((k : ℕ) : ℝ) * (-beta * δ) by ring]
      push_cast
      ring
  calc (∑ x ∈ Finset.univ.filter
        (fun x : LadderState n => ladderEnergy δ n x ≠ 0),
      gibbsWeight (ladderEnergy δ n) beta x)
      = ∑ x : LadderState n,
          if ladderEnergy δ n x ≠ 0
          then gibbsWeight (ladderEnergy δ n) beta x else 0 :=
        Finset.sum_filter _ _
    _ = ∑ k : Fin (n + 1), ∑ i : Fin (2 ^ (k : ℕ)),
          if ladderEnergy δ n ⟨k, i⟩ ≠ 0
          then gibbsWeight (ladderEnergy δ n) beta ⟨k, i⟩ else 0 :=
        Fintype.sum_sigma _
    _ = ∑ k : Fin (n + 1),
          (fun m : ℕ =>
            if m ≠ 0 then (2 * Real.exp (-beta * δ)) ^ m else 0)
            (k : ℕ) :=
        Finset.sum_congr rfl fun k _ => hterm k
    _ = ∑ m ∈ Finset.range (n + 1),
          (if m ≠ 0 then (2 * Real.exp (-beta * δ)) ^ m else 0) :=
        Fin.sum_univ_eq_sum_range
          (fun m : ℕ =>
            if m ≠ 0 then (2 * Real.exp (-beta * δ)) ^ m else 0)
          (n + 1)
    _ = ∑ j ∈ Finset.range n, (2 * Real.exp (-beta * δ)) ^ (j + 1) := by
        rw [Finset.sum_range_succ']
        simp

/-- The geometric envelope receipt of the ladder: whenever
`r = 2 * exp (-beta * δ) < 1`, every stage `n` holds its off-minimum
Gibbs mass below `r / (1 - r)`, uniformly in `n`. -/
theorem ladder_offMinMass_le (δ : ℝ) (hδ : 0 < δ) (n : ℕ) (beta : ℝ)
    (hlt : 2 * Real.exp (-beta * δ) < 1) :
    offMinMass (ladderEnergy δ n) beta
      ≤ 2 * Real.exp (-beta * δ) / (1 - 2 * Real.exp (-beta * δ)) := by
  have hR0 : (0 : ℝ) ≤ 2 * Real.exp (-beta * δ) := by positivity
  have hnumnn : (0 : ℝ) ≤ ∑ j ∈ Finset.range n,
      (2 * Real.exp (-beta * δ)) ^ (j + 1) :=
    Finset.sum_nonneg fun j _ => pow_nonneg hR0 (j + 1)
  have hZ := ladder_partitionZ_ge_one δ n beta
  unfold offMinMass excitedMass
  rw [ladder_minEnergy δ hδ.le n, ladder_excitedSum_eq δ hδ n beta]
  calc (∑ j ∈ Finset.range n, (2 * Real.exp (-beta * δ)) ^ (j + 1))
        / partitionZ (ladderEnergy δ n) beta
      ≤ ∑ j ∈ Finset.range n, (2 * Real.exp (-beta * δ)) ^ (j + 1) :=
        div_le_self hnumnn hZ
    _ ≤ 2 * Real.exp (-beta * δ)
          / (1 - 2 * Real.exp (-beta * δ)) :=
        geom_tail_sum_le _ hR0 hlt n

/-- The explicit threshold: at and beyond `beta0 = log 4 / δ` the ratio
`r = 2 * exp (-beta * δ)` is at most `1/2`, hence strictly below `1`. -/
theorem ladder_r_lt_one (δ : ℝ) (hδ : 0 < δ) (beta : ℝ)
    (hbeta : Real.log 4 / δ ≤ beta) :
    2 * Real.exp (-beta * δ) < 1 := by
  have hβδ : Real.log 4 ≤ beta * δ := by
    rw [div_le_iff₀ hδ] at hbeta
    linarith
  have h1 : Real.exp (-(beta * δ)) ≤ Real.exp (-Real.log 4) :=
    Real.exp_le_exp.mpr (by linarith)
  have h2 : Real.exp (-Real.log 4) = 1 / 4 := by
    rw [Real.exp_neg, Real.exp_log (by norm_num : (0 : ℝ) < 4)]
    norm_num
  have h3 : -beta * δ = -(beta * δ) := by ring
  rw [h3]
  linarith

/-- **The growing-cardinality envelope instance.** The ladder family inhabits the
coherent structure: stage truncations are nontrivial coherent refine
maps, and the geometric spectral tail `r / (1 - r)` with
`r = 2 * exp (-beta * δ)` is a uniform envelope beyond the explicit
threshold `log 4 / δ`. Its stage cardinality is unbounded
(`ladder_card_unbounded`), so the committed gap-plus-cardinality
structure is uninhabited on the same data
(`ladder_uniformGap_isEmpty`).  Thus an envelope does not imply a
uniform total-cardinality bound on the same carriers and energies.  This
is a separation of those two conditions, not a structure-inclusion
theorem. -/
noncomputable def ladderFamily (δ : ℝ) (hδ : 0 < δ) :
    CoherentRefinementFamily LadderState (ladderEnergy δ) where
  directed := fun r s => ⟨max r s, le_max_left r s, le_max_right r s⟩
  refineMap := fun {r s} _ x => ladderTrunc r x
  refine_id := fun r x => ladderTrunc_id r x
  refine_comp := fun {r s t} hrs _ x => ladderTrunc_comp hrs x
  refine_energy_le := fun {r s} _ x => ladderTrunc_energy_le δ hδ.le r x
  ground := ladderGround
  ground_min := fun n => ladderGround_min δ hδ.le n
  refine_ground := fun {r s} _ => ladderTrunc_ground r s
  beta0 := Real.log 4 / δ
  Env := fun beta =>
    2 * Real.exp (-beta * δ) / (1 - 2 * Real.exp (-beta * δ))
  env_bound := fun beta hbeta r =>
    ladder_offMinMass_le δ hδ r beta (ladder_r_lt_one δ hδ beta hbeta)
  env_tendsto := by
    have hnum : Filter.Tendsto
        (fun beta : ℝ => 2 * Real.exp (-beta * δ))
        Filter.atTop (nhds 0) :=
      tendsto_const_mul_exp_neg_gap 2 δ hδ
    have hden : Filter.Tendsto
        (fun beta : ℝ => 1 - 2 * Real.exp (-beta * δ))
        Filter.atTop (nhds 1) := by
      have h := hnum.const_sub 1
      simpa using h
    have h := hnum.div hden one_ne_zero
    simpa using h

/-- Stage `n` of the ladder has at least `2^n` states: level `n` alone
carries `2^n` copies. -/
theorem ladder_card_lower (n : ℕ) :
    2 ^ n ≤ Fintype.card (LadderState n) := by
  rw [Fintype.card_sigma]
  calc 2 ^ n
      = Fintype.card
          (Fin (2 ^ ((⟨n, Nat.lt_succ_self n⟩ : Fin (n + 1)) : ℕ))) := by
        rw [Fintype.card_fin]
    _ ≤ ∑ k : Fin (n + 1), Fintype.card (Fin (2 ^ (k : ℕ))) :=
        Finset.single_le_sum
          (f := fun k : Fin (n + 1) => Fintype.card (Fin (2 ^ (k : ℕ))))
          (fun k _ => Nat.zero_le _)
          (Finset.mem_univ (⟨n, Nat.lt_succ_self n⟩ : Fin (n + 1)))

/-- The ladder stage cardinality is unbounded. -/
theorem ladder_card_unbounded (m : ℕ) :
    m < Fintype.card (LadderState m) :=
  lt_of_lt_of_le Nat.lt_two_pow_self (ladder_card_lower m)

/-- No uniform cardinality bound exists across the ladder stages: the
committed uniform-cardinality hypothesis fails on this family. -/
theorem ladder_no_uniform_card :
    ¬ ∃ cardBound : ℕ, ∀ n : ℕ,
      Fintype.card (LadderState n) ≤ cardBound := by
  rintro ⟨cB, hcB⟩
  exact absurd (hcB cB) (not_le.mpr (ladder_card_unbounded cB))

/-- **Cardinality-separation receipt.** The committed gap-plus-cardinality
structure is uninhabited on the ladder data: any instance would carry
a uniform cardinality bound, and `ladder_card_unbounded` refutes it.
Together with `ladderFamily`, this proves that the spectral-envelope
condition does not imply a uniform total-cardinality bound on the same
carriers and energies. -/
theorem ladder_uniformGap_isEmpty (δ : ℝ) :
    IsEmpty (UniformGapRefinement LadderState (ladderEnergy δ)) :=
  ⟨fun F => absurd (F.card_le F.cardBound)
    (not_le.mpr (ladder_card_unbounded F.cardBound))⟩

/-! ## Coherent-envelope companion to the composed four-law surface

Register row PR-08, narrowed form.  The declarations below leave every
committed declaration of `FourLawAdequacySurface` untouched.  They pair
the old conclusion record, reached through a separately constructed
constant-family witness, with consequences of the supplied directed
family. -/

universe u

/-- **Register row PR-08, coherent-envelope attachment.** The coherent
directed refinement family, tied to the calibrated energy of register
row PR-15 at the distinguished member through a carrier
identification, exactly as `RefinementAttachment` ties the uniform-gap
form. -/
structure CoherentRefinementAttachment {Ω : Type u} [Fintype Ω]
    [DecidableEq Ω] {B : Type*} [DecidableEq B] (D : RepairLawData Ω B)
    (C : CalibrationData D) where
  /-- The refinement index type. -/
  I : Type
  /-- The refinement order on the index type. -/
  [preorderI : Preorder I]
  /-- The member carriers of the declared family. -/
  X : I → Type u
  /-- Every member carrier is finite. -/
  [fintypeX : ∀ r, Fintype (X r)]
  /-- Every member carrier has decidable equality. -/
  [decEqX : ∀ r, DecidableEq (X r)]
  /-- Every member carrier is inhabited. -/
  [nonemptyX : ∀ r, Nonempty (X r)]
  /-- The member energies of the declared family. -/
  E : ∀ r, X r → ℝ
  /-- The declared coherent directed family with a uniform spectral-tail
  envelope: register row PR-08 in the narrowed form. -/
  family : CoherentRefinementFamily X E
  /-- The distinguished member identified with the surface state
  space. -/
  r0 : I
  /-- The carrier identification of the distinguished member with the
  surface state space. -/
  e : X r0 ≃ Ω
  /-- The distinguished member's energy is the calibrated energy of
  register row PR-15 through the carrier identification. -/
  energy_eq : ∀ x, E r0 x = C.energy (e x)

attribute [instance] CoherentRefinementAttachment.preorderI
attribute [instance] CoherentRefinementAttachment.fintypeX
attribute [instance] CoherentRefinementAttachment.decEqX
attribute [instance] CoherentRefinementAttachment.nonemptyX

/-- **The coherent composed antecedent.** One typed bundle whose
fields are the premise-register rows PR-07 and PR-15 exactly as in
`FourLawAntecedent`, with row PR-08 in the narrowed coherent-envelope
form. -/
structure FourLawCoherentAntecedent (Ω : Type u) [Fintype Ω]
    [DecidableEq Ω] (B : Type*) [DecidableEq B] where
  /-- Register row PR-07: the declared repair law. -/
  repair : RepairLawData Ω B
  /-- Register row PR-15: the clock and energy calibration of the
  declared reference. -/
  calib : CalibrationData repair
  /-- Register row PR-08, narrowed form: the coherent directed family
  with a uniform spectral-tail envelope, attached to the calibrated
  energy. -/
  refinement : CoherentRefinementAttachment repair calib

/-- The coherent bundle forces an inhabited state space through the
carrier identification of the distinguished member. -/
theorem FourLawCoherentAntecedent.nonempty {Ω : Type u} [Fintype Ω]
    [DecidableEq Ω] {B : Type*} [DecidableEq B]
    (A : FourLawCoherentAntecedent Ω B) : Nonempty Ω :=
  A.refinement.e.nonempty_congr.mp
    (A.refinement.nonemptyX A.refinement.r0)

/-- The route from the coherent antecedent to the committed
antecedent: the PR-07 and PR-15 rows are carried over unchanged
(`toAntecedent_repair`, `toAntecedent_calib`, both definitional), and
the committed PR-08 field is filled by the degenerate one-member
constant family on the calibrated energy (`constUniformGap`), whose
gap is the exact gap of the calibrated energy when the off-minimum set
is inhabited and vacuous otherwise. -/
noncomputable def FourLawCoherentAntecedent.toAntecedent {Ω : Type u}
    [Fintype Ω] [DecidableEq Ω] {B : Type*} [DecidableEq B]
    (A : FourLawCoherentAntecedent Ω B) : FourLawAntecedent Ω B :=
  have : Nonempty Ω := A.nonempty
  { repair := A.repair
    calib := A.calib
    refinement :=
      { I := ℕ
        X := fun _ => Ω
        fintypeX := fun _ => inferInstance
        decEqX := fun _ => inferInstance
        nonemptyX := fun _ => inferInstance
        E := fun _ => A.calib.energy
        family := constUniformGap ℕ Ω A.calib.energy
        r0 := 0
        e := Equiv.refl Ω
        energy_eq := fun _ => rfl } }

/-- The coherent route carries register row PR-07 unchanged. -/
theorem FourLawCoherentAntecedent.toAntecedent_repair {Ω : Type u}
    [Fintype Ω] [DecidableEq Ω] {B : Type*} [DecidableEq B]
    (A : FourLawCoherentAntecedent Ω B) :
    A.toAntecedent.repair = A.repair := rfl

/-- The coherent route carries register row PR-15 unchanged. -/
theorem FourLawCoherentAntecedent.toAntecedent_calib {Ω : Type u}
    [Fintype Ω] [DecidableEq Ω] {B : Type*} [DecidableEq B]
    (A : FourLawCoherentAntecedent Ω B) :
    A.toAntecedent.calib = A.calib := rfl

/-- **Four-law companion conjunction for the coherent-envelope row.**
From the bundled data, the full old conclusion package holds along the
separately constructed degenerate one-member route (`toAntecedent`;
rows PR-07 and PR-15 are carried over definitionally).  The supplied
directed family is not installed in that old conclusion record.
Separately, its uniform spectral-tail envelope
controls every member of the declared coherent family beyond the
declared threshold; for every positive tolerance one
inverse-temperature threshold takes every member's off-minimum mass
below the tolerance; the calibrated energy is the distinguished member
of the coherent family through the carrier identification; and the
calibrated energy carries the envelope bound beyond the declared
threshold. The `Nonempty` instance is derivable from the bundle
through `FourLawCoherentAntecedent.nonempty`. -/
theorem fourLaws_composed_coherent {Ω : Type u} [Fintype Ω]
    [DecidableEq Ω] [Nonempty Ω] {B : Type*} [DecidableEq B]
    (A : FourLawCoherentAntecedent Ω B) :
    FourLawConclusions A.toAntecedent
    ∧ (∀ beta : ℝ, A.refinement.family.beta0 ≤ beta →
        ∀ r : A.refinement.I,
          offMinMass (A.refinement.E r) beta
            ≤ A.refinement.family.Env beta)
    ∧ (∀ eps : ℝ, 0 < eps → ∃ beta1 : ℝ, ∀ beta : ℝ, beta1 ≤ beta →
        ∀ r : A.refinement.I,
          offMinMass (A.refinement.E r) beta < eps)
    ∧ (∀ beta : ℝ,
        offMinMass A.calib.energy beta
          = offMinMass (A.refinement.E A.refinement.r0) beta)
    ∧ (∀ beta : ℝ, A.refinement.family.beta0 ≤ beta →
        offMinMass A.calib.energy beta
          ≤ A.refinement.family.Env beta) := by
  have hmember : ∀ beta : ℝ,
      offMinMass A.calib.energy beta
        = offMinMass (A.refinement.E A.refinement.r0) beta := by
    intro beta
    have hE : A.refinement.E A.refinement.r0
        = fun x => A.calib.energy (A.refinement.e x) :=
      funext A.refinement.energy_eq
    rw [hE]
    exact (offMinMass_equiv A.refinement.e A.calib.energy beta).symm
  refine ⟨fourLaws_composed A.toAntecedent,
    A.refinement.family.env_bound,
    A.refinement.family.uniform_concentration,
    hmember, fun beta hbeta => ?_⟩
  rw [hmember beta]
  exact A.refinement.family.env_bound beta hbeta A.refinement.r0

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.geom_tail_sum_le
#print axioms OPH.Thermodynamics.geom_tail_tsum
#print axioms OPH.Thermodynamics.CoherentRefinementFamily.minEnergy_mono
#print axioms OPH.Thermodynamics.CoherentRefinementFamily.refine_ground_min
#print axioms OPH.Thermodynamics.CoherentRefinementFamily.refine_energy_le_chain
#print axioms OPH.Thermodynamics.CoherentRefinementFamily.envelope_bound
#print axioms OPH.Thermodynamics.CoherentRefinementFamily.uniform_concentration
#print axioms OPH.Thermodynamics.CoherentRefinementFamily.ofUniformGap
#print axioms OPH.Thermodynamics.CoherentRefinementFamily.ofUniformGapConst
#print axioms OPH.Thermodynamics.twoLevelCoherent
#print axioms OPH.Thermodynamics.constUniformGap
#print axioms OPH.Thermodynamics.ladder_ground_energy
#print axioms OPH.Thermodynamics.ladder_minEnergy
#print axioms OPH.Thermodynamics.ladderTrunc_id
#print axioms OPH.Thermodynamics.ladderTrunc_ground
#print axioms OPH.Thermodynamics.ladderTrunc_comp
#print axioms OPH.Thermodynamics.ladderTrunc_energy_le
#print axioms OPH.Thermodynamics.ladderGround_min
#print axioms OPH.Thermodynamics.ladder_partitionZ_ge_one
#print axioms OPH.Thermodynamics.ladder_excitedSum_eq
#print axioms OPH.Thermodynamics.ladder_offMinMass_le
#print axioms OPH.Thermodynamics.ladder_r_lt_one
#print axioms OPH.Thermodynamics.ladderFamily
#print axioms OPH.Thermodynamics.ladder_card_lower
#print axioms OPH.Thermodynamics.ladder_card_unbounded
#print axioms OPH.Thermodynamics.ladder_no_uniform_card
#print axioms OPH.Thermodynamics.ladder_uniformGap_isEmpty
#print axioms OPH.Thermodynamics.FourLawCoherentAntecedent.nonempty
#print axioms OPH.Thermodynamics.FourLawCoherentAntecedent.toAntecedent
#print axioms OPH.Thermodynamics.FourLawCoherentAntecedent.toAntecedent_repair
#print axioms OPH.Thermodynamics.FourLawCoherentAntecedent.toAntecedent_calib
#print axioms OPH.Thermodynamics.fourLaws_composed_coherent
-- Expected axioms for every declaration above:
-- propext, Classical.choice, Quot.sound. No native_decide is used.
