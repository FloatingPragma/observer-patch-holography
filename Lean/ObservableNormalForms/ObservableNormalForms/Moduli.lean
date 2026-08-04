import ObservableNormalForms.Stability

/-!
# Intrinsic finite stability moduli

This module replaces the *assumed* certificates of `Stability.lean`
(`ErrorBoundWitness`, `InverseObservationBound`) with *defined* finite
maxima on a finite metric state space, and proves the manuscript's
`prop:minimal-moduli`, `cor:rate-transfer`, and `thm:sharpness`
(`extra/observable_normal_forms.tex`, Section "Observational stability
moduli").

The three intrinsic objects are

* `residualModulus` (ηΦ): the largest distance to the consistent set among
  states of residual at most `t`;
* `inverseModulus` (ωB): the largest state distance among consistent pairs
  whose observations differ by at most `r`;
* `settledModulus` (Θ): the largest state distance among pairs that are both
  `δ`-settled and `r`-close in observation.

Every extremum is a `Finset.sup'`/`Finset.inf'` over a set whose
nonemptiness is *proved* at each use site, never assumed: the sublevel set
at `t ≥ 0` contains every consistent state, the inverse and settled sets at
nonnegative radii contain the diagonal pair on a consistent state, and the
nearest-point infimum carries the nonemptiness of `C` as an explicit
argument.  The `dite` fallback value `0` is unreachable in every theorem
below; each proof routes through an explicit `*_of_nonempty` rewrite.

The sharpness witness (three points of ℝ) and the negative control (two
points of ℝ) compute their moduli exactly rather than assuming them.
-/

namespace ObservableNormalForms

universe u v

section Intrinsic

variable {Q : Type u} {𝓑 : Type v}
variable [MetricSpace Q] [Fintype Q] [PseudoMetricSpace 𝓑]

/-! ### Distance to a finite consistent set -/

/-- Distance from `x` to the nonempty finite set `C`, as an attained finite
minimum.  Nonemptiness is an explicit argument, not a default value. -/
noncomputable def distToFinset (C : Finset Q) (hC : C.Nonempty) (x : Q) : ℝ :=
  C.inf' hC fun c => dist x c

omit [Fintype Q] in
lemma distToFinset_nonneg (C : Finset Q) (hC : C.Nonempty) (x : Q) :
    0 ≤ distToFinset C hC x :=
  Finset.le_inf' hC _ fun _ _ => dist_nonneg

omit [Fintype Q] in
lemma distToFinset_le (C : Finset Q) (hC : C.Nonempty) {c : Q}
    (hc : c ∈ C) (x : Q) :
    distToFinset C hC x ≤ dist x c :=
  Finset.inf'_le _ hc

omit [Fintype Q] in
/-- The nearest-point minimum is attained. -/
lemma exists_distToFinset_eq (C : Finset Q) (hC : C.Nonempty) (x : Q) :
    ∃ c ∈ C, distToFinset C hC x = dist x c :=
  Finset.exists_mem_eq_inf' hC _

omit [Fintype Q] in
lemma distToFinset_eq_zero_of_mem (C : Finset Q) (hC : C.Nonempty)
    {x : Q} (hx : x ∈ C) :
    distToFinset C hC x = 0 :=
  le_antisymm (by simpa using distToFinset_le C hC hx x)
    (distToFinset_nonneg C hC x)

/-! ### The residual error-bound modulus ηΦ -/

open Classical in
/-- The residual sublevel set `{x : Φ x ≤ t}`. -/
noncomputable def residualSet (Φ : Q → ℝ) (t : ℝ) : Finset Q :=
  Finset.univ.filter fun x => Φ x ≤ t

omit [MetricSpace Q] in
lemma mem_residualSet {Φ : Q → ℝ} {t : ℝ} {x : Q} :
    x ∈ residualSet Φ t ↔ Φ x ≤ t := by
  simp [residualSet]

omit [MetricSpace Q] in
lemma residualSet_mono {Φ : Q → ℝ} {s t : ℝ} (hst : s ≤ t) :
    residualSet Φ s ⊆ residualSet Φ t := fun _ hx =>
  mem_residualSet.mpr ((mem_residualSet.mp hx).trans hst)

omit [MetricSpace Q] in
/-- Nonemptiness of the sublevel set at any nonnegative threshold, witnessed
by a consistent state.  This is the discharged witness that keeps every
downstream ηΦ statement out of the vacuous `dite` branch. -/
lemma residualSet_nonempty {C : Finset Q} (hC : C.Nonempty) {Φ : Q → ℝ}
    (hCzero : ∀ c ∈ C, Φ c = 0) {t : ℝ} (ht : 0 ≤ t) :
    (residualSet Φ t).Nonempty := by
  obtain ⟨c, hc⟩ := hC
  exact ⟨c, mem_residualSet.mpr (by rw [hCzero c hc]; exact ht)⟩

/-- The intrinsic residual error-bound modulus
`ηΦ(t) = max {dist(x, C) : Φ x ≤ t}`. -/
noncomputable def residualModulus (C : Finset Q) (hC : C.Nonempty)
    (Φ : Q → ℝ) (t : ℝ) : ℝ :=
  if h : (residualSet Φ t).Nonempty then
    (residualSet Φ t).sup' h (distToFinset C hC)
  else 0

lemma residualModulus_of_nonempty {C : Finset Q} (hC : C.Nonempty)
    {Φ : Q → ℝ} {t : ℝ} (h : (residualSet Φ t).Nonempty) :
    residualModulus C hC Φ t = (residualSet Φ t).sup' h (distToFinset C hC) :=
  dif_pos h

lemma residualModulus_of_empty {C : Finset Q} (hC : C.Nonempty)
    {Φ : Q → ℝ} {t : ℝ} (h : ¬(residualSet Φ t).Nonempty) :
    residualModulus C hC Φ t = 0 :=
  dif_neg h

lemma residualModulus_nonneg (C : Finset Q) (hC : C.Nonempty)
    (Φ : Q → ℝ) (t : ℝ) :
    0 ≤ residualModulus C hC Φ t := by
  by_cases h : (residualSet Φ t).Nonempty
  · obtain ⟨x, hx⟩ := h
    rw [residualModulus_of_nonempty hC ⟨x, hx⟩]
    exact (distToFinset_nonneg C hC x).trans (Finset.le_sup' _ hx)
  · rw [residualModulus_of_empty hC h]

/-- Any state of residual at most `t` lies within `ηΦ(t)` of the consistent
set.  The sublevel set is nonempty because the state itself belongs to it. -/
lemma distToFinset_le_residualModulus {C : Finset Q} (hC : C.Nonempty)
    {Φ : Q → ℝ} {t : ℝ} {x : Q} (hx : Φ x ≤ t) :
    distToFinset C hC x ≤ residualModulus C hC Φ t := by
  have hmem : x ∈ residualSet Φ t := mem_residualSet.mpr hx
  rw [residualModulus_of_nonempty hC ⟨x, hmem⟩]
  exact Finset.le_sup' _ hmem

/-- `prop:minimal-moduli` (monotonicity of ηΦ). -/
theorem residualModulus_mono (C : Finset Q) (hC : C.Nonempty) (Φ : Q → ℝ) :
    Monotone (residualModulus C hC Φ) := by
  intro s t hst
  by_cases h : (residualSet Φ s).Nonempty
  · rw [residualModulus_of_nonempty hC h,
      residualModulus_of_nonempty hC (h.mono (residualSet_mono hst))]
    exact Finset.sup'_le _ _ fun x hx =>
      Finset.le_sup' _ (residualSet_mono hst hx)
  · rw [residualModulus_of_empty hC h]
    exact residualModulus_nonneg C hC Φ t

/-- `prop:minimal-moduli` (i): `ηΦ(0) = 0`. -/
theorem residualModulus_zero {C : Finset Q} (hC : C.Nonempty) {Φ : Q → ℝ}
    (hΦnn : ∀ x, 0 ≤ Φ x) (hCiff : ∀ x, x ∈ C ↔ Φ x = 0) :
    residualModulus C hC Φ 0 = 0 := by
  have hCzero : ∀ c ∈ C, Φ c = 0 := fun c hc => (hCiff c).mp hc
  have hne := residualSet_nonempty hC hCzero le_rfl
  rw [residualModulus_of_nonempty hC hne]
  refine le_antisymm (Finset.sup'_le _ _ fun x hx => ?_) ?_
  · have hx0 : Φ x = 0 := le_antisymm (mem_residualSet.mp hx) (hΦnn x)
    exact le_of_eq (distToFinset_eq_zero_of_mem C hC ((hCiff x).mpr hx0))
  · obtain ⟨c, hc⟩ := id hC
    have hmem : c ∈ residualSet Φ 0 :=
      mem_residualSet.mpr (le_of_eq (hCzero c hc))
    exact (distToFinset_nonneg C hC c).trans (Finset.le_sup' _ hmem)

/-- The intrinsic residual modulus satisfies the reusable certificate of
`Stability.lean`.  This turns the previously *assumed* `ErrorBoundWitness`
into a theorem about a defined object. -/
theorem errorBoundWitness_residualModulus (C : Finset Q) (hC : C.Nonempty)
    (Φ : Q → ℝ) :
    ErrorBoundWitness (↑C : Set Q) Φ (residualModulus C hC Φ) := by
  intro x δ hx
  obtain ⟨c, hcC, hcd⟩ := exists_distToFinset_eq C hC x
  refine ⟨c, Finset.mem_coe.mpr hcC, ?_⟩
  calc dist x c = distToFinset C hC x := hcd.symm
    _ ≤ residualModulus C hC Φ δ := distToFinset_le_residualModulus hC hx

/-- `prop:minimal-moduli` (ii): ηΦ is the least error-bound certificate.  Any
`η` satisfying `ErrorBoundWitness` dominates ηΦ at every nonnegative
threshold; no monotonicity of `η` is required. -/
theorem residualModulus_le_of_errorBoundWitness {C : Finset Q}
    (hC : C.Nonempty) {Φ : Q → ℝ} {η : ℝ → ℝ}
    (hη : ErrorBoundWitness (↑C : Set Q) Φ η)
    (hCzero : ∀ c ∈ C, Φ c = 0) {t : ℝ} (ht : 0 ≤ t) :
    residualModulus C hC Φ t ≤ η t := by
  have hne := residualSet_nonempty hC hCzero ht
  rw [residualModulus_of_nonempty hC hne]
  refine Finset.sup'_le _ _ fun z hz => ?_
  obtain ⟨c, hcC, hcd⟩ := hη (mem_residualSet.mp hz)
  exact (distToFinset_le C hC (Finset.mem_coe.mp hcC) z).trans hcd

/-! ### The inverse-observation modulus ωB -/

open Classical in
/-- Consistent pairs whose observations are `r`-close. -/
noncomputable def inverseSet (C : Finset Q) (B : Q → 𝓑) (r : ℝ) :
    Finset (Q × Q) :=
  (C ×ˢ C).filter fun p => dist (B p.1) (B p.2) ≤ r

omit [MetricSpace Q] [Fintype Q] in
lemma mem_inverseSet {C : Finset Q} {B : Q → 𝓑} {r : ℝ} {p : Q × Q} :
    p ∈ inverseSet C B r ↔
      p.1 ∈ C ∧ p.2 ∈ C ∧ dist (B p.1) (B p.2) ≤ r := by
  simp [inverseSet, Finset.mem_product, and_assoc]

omit [MetricSpace Q] [Fintype Q] in
lemma inverseSet_mono {C : Finset Q} {B : Q → 𝓑} {r s : ℝ} (hrs : r ≤ s) :
    inverseSet C B r ⊆ inverseSet C B s := by
  intro p hp
  rw [mem_inverseSet] at hp ⊢
  exact ⟨hp.1, hp.2.1, hp.2.2.trans hrs⟩

omit [MetricSpace Q] [Fintype Q] in
/-- Nonemptiness of the inverse-observation set at any nonnegative radius,
witnessed by the diagonal pair on a consistent state. -/
lemma inverseSet_nonempty {C : Finset Q} (hC : C.Nonempty) {B : Q → 𝓑}
    {r : ℝ} (hr : 0 ≤ r) :
    (inverseSet C B r).Nonempty := by
  obtain ⟨c, hc⟩ := hC
  exact ⟨(c, c), mem_inverseSet.mpr ⟨hc, hc, by simpa using hr⟩⟩

/-- The intrinsic inverse-observation modulus
`ωB(r) = max {dist(c, c') : c, c' ∈ C, dist(Bc, Bc') ≤ r}`. -/
noncomputable def inverseModulus (C : Finset Q) (B : Q → 𝓑) (r : ℝ) : ℝ :=
  if h : (inverseSet C B r).Nonempty then
    (inverseSet C B r).sup' h fun p => dist p.1 p.2
  else 0

omit [Fintype Q] in
lemma inverseModulus_of_nonempty {C : Finset Q} {B : Q → 𝓑} {r : ℝ}
    (h : (inverseSet C B r).Nonempty) :
    inverseModulus C B r =
      (inverseSet C B r).sup' h fun p => dist p.1 p.2 :=
  dif_pos h

omit [Fintype Q] in
lemma inverseModulus_of_empty {C : Finset Q} {B : Q → 𝓑} {r : ℝ}
    (h : ¬(inverseSet C B r).Nonempty) :
    inverseModulus C B r = 0 :=
  dif_neg h

omit [Fintype Q] in
lemma inverseModulus_nonneg (C : Finset Q) (B : Q → 𝓑) (r : ℝ) :
    0 ≤ inverseModulus C B r := by
  by_cases h : (inverseSet C B r).Nonempty
  · obtain ⟨p, hp⟩ := h
    rw [inverseModulus_of_nonempty ⟨p, hp⟩]
    exact dist_nonneg.trans
      (Finset.le_sup' (fun q : Q × Q => dist q.1 q.2) hp)
  · rw [inverseModulus_of_empty h]

omit [Fintype Q] in
/-- Any consistent pair is bounded by ωB at its own observation distance.
The inverse-observation set is nonempty because the pair itself belongs. -/
lemma dist_le_inverseModulus {C : Finset Q} {B : Q → 𝓑} {c d : Q}
    (hc : c ∈ C) (hd : d ∈ C) :
    dist c d ≤ inverseModulus C B (dist (B c) (B d)) := by
  have hmem : (c, d) ∈ inverseSet C B (dist (B c) (B d)) :=
    mem_inverseSet.mpr ⟨hc, hd, le_rfl⟩
  rw [inverseModulus_of_nonempty ⟨(c, d), hmem⟩]
  exact Finset.le_sup' (fun q : Q × Q => dist q.1 q.2) hmem

omit [Fintype Q] in
/-- `prop:minimal-moduli` (monotonicity of ωB). -/
theorem inverseModulus_mono (C : Finset Q) (B : Q → 𝓑) :
    Monotone (inverseModulus C B) := by
  intro r s hrs
  by_cases h : (inverseSet C B r).Nonempty
  · rw [inverseModulus_of_nonempty h,
      inverseModulus_of_nonempty (h.mono (inverseSet_mono hrs))]
    exact Finset.sup'_le _ (fun q : Q × Q => dist q.1 q.2) fun p hp =>
      Finset.le_sup' (fun q : Q × Q => dist q.1 q.2)
        (inverseSet_mono hrs hp)
  · rw [inverseModulus_of_empty h]
    exact inverseModulus_nonneg C B s

omit [Fintype Q] in
/-- The intrinsic inverse modulus satisfies the reusable certificate of
`Stability.lean`.  This turns the previously *assumed*
`InverseObservationBound` into a theorem about a defined object. -/
theorem inverseObservationBound_inverseModulus (C : Finset Q) (B : Q → 𝓑) :
    InverseObservationBound (↑C : Set Q) B (inverseModulus C B) :=
  fun _ _ hc hd =>
    dist_le_inverseModulus (Finset.mem_coe.mp hc) (Finset.mem_coe.mp hd)

omit [Fintype Q] in
/-- `prop:minimal-moduli` (iv): ωB is the least monotone inverse-observation
certificate at every nonnegative radius. -/
theorem inverseModulus_le_of_inverseObservationBound {C : Finset Q}
    (hC : C.Nonempty) {B : Q → 𝓑} {ω : ℝ → ℝ}
    (hω : InverseObservationBound (↑C : Set Q) B ω) (hmono : Monotone ω)
    {r : ℝ} (hr : 0 ≤ r) :
    inverseModulus C B r ≤ ω r := by
  have hne := inverseSet_nonempty hC (B := B) hr
  rw [inverseModulus_of_nonempty hne]
  refine Finset.sup'_le _ _ fun p hp => ?_
  rw [mem_inverseSet] at hp
  exact (hω (Finset.mem_coe.mpr hp.1) (Finset.mem_coe.mpr hp.2.1)).trans
    (hmono hp.2.2)

/-! ### The settled-output modulus Θ -/

open Classical in
/-- Pairs that are `δ`-settled in residual and `r`-close in observation. -/
noncomputable def settledSet (Φ : Q → ℝ) (B : Q → 𝓑) (δ r : ℝ) :
    Finset (Q × Q) :=
  Finset.univ.filter fun p : Q × Q =>
    Φ p.1 ≤ δ ∧ Φ p.2 ≤ δ ∧ dist (B p.1) (B p.2) ≤ r

omit [MetricSpace Q] in
lemma mem_settledSet {Φ : Q → ℝ} {B : Q → 𝓑} {δ r : ℝ} {p : Q × Q} :
    p ∈ settledSet Φ B δ r ↔
      Φ p.1 ≤ δ ∧ Φ p.2 ≤ δ ∧ dist (B p.1) (B p.2) ≤ r := by
  simp [settledSet]

omit [MetricSpace Q] in
/-- Nonemptiness of the settled set at nonnegative parameters, witnessed by
the diagonal pair on a consistent state. -/
lemma settledSet_nonempty {C : Finset Q} (hC : C.Nonempty) {Φ : Q → ℝ}
    {B : Q → 𝓑} (hCzero : ∀ c ∈ C, Φ c = 0) {δ r : ℝ}
    (hδ : 0 ≤ δ) (hr : 0 ≤ r) :
    (settledSet Φ B δ r).Nonempty := by
  obtain ⟨c, hc⟩ := hC
  refine ⟨(c, c), mem_settledSet.mpr ⟨?_, ?_, by simpa using hr⟩⟩ <;>
    simpa [hCzero c hc] using hδ

/-- The intrinsic settled-output modulus
`Θ(δ, r) = max {dist(x, y) : Φ x, Φ y ≤ δ, dist(Bx, By) ≤ r}`. -/
noncomputable def settledModulus (Φ : Q → ℝ) (B : Q → 𝓑) (δ r : ℝ) : ℝ :=
  if h : (settledSet Φ B δ r).Nonempty then
    (settledSet Φ B δ r).sup' h fun p => dist p.1 p.2
  else 0

lemma settledModulus_of_nonempty {Φ : Q → ℝ} {B : Q → 𝓑} {δ r : ℝ}
    (h : (settledSet Φ B δ r).Nonempty) :
    settledModulus Φ B δ r =
      (settledSet Φ B δ r).sup' h fun p => dist p.1 p.2 :=
  dif_pos h

lemma settledModulus_of_empty {Φ : Q → ℝ} {B : Q → 𝓑} {δ r : ℝ}
    (h : ¬(settledSet Φ B δ r).Nonempty) :
    settledModulus Φ B δ r = 0 :=
  dif_neg h

lemma settledModulus_nonneg (Φ : Q → ℝ) (B : Q → 𝓑) (δ r : ℝ) :
    0 ≤ settledModulus Φ B δ r := by
  by_cases h : (settledSet Φ B δ r).Nonempty
  · obtain ⟨p, hp⟩ := h
    rw [settledModulus_of_nonempty ⟨p, hp⟩]
    exact dist_nonneg.trans
      (Finset.le_sup' (fun q : Q × Q => dist q.1 q.2) hp)
  · rw [settledModulus_of_empty h]

/-! ### Intrinsic corollaries of the master estimate -/

/-- `thm:master-bound` with the intrinsic moduli.  The previously assumed
certificates are here *instantiated by theorems*, so the estimate holds for
the defined finite maxima with no residual hypotheses about them. -/
theorem intrinsic_two_output_estimate {C : Finset Q} (hC : C.Nonempty)
    {Φ : Q → ℝ} {B : Q → 𝓑} {L δx δy ε : ℝ} {x y : Q}
    (hL : 0 ≤ L) (hLip : LipschitzBound L B)
    (hx : Φ x ≤ δx) (hy : Φ y ≤ δy)
    (hBxy : dist (B x) (B y) ≤ ε) :
    dist x y ≤
      residualModulus C hC Φ δx + residualModulus C hC Φ δy +
        inverseModulus C B
          (ε + L * (residualModulus C hC Φ δx + residualModulus C hC Φ δy)) :=
  heterogeneous_two_output_estimate hL
    (errorBoundWitness_residualModulus C hC Φ)
    (inverseObservationBound_inverseModulus C B)
    (inverseModulus_mono C B) hLip hx hy hBxy

/-- `eq:symmetric-bound` with the intrinsic moduli. -/
theorem intrinsic_symmetric_estimate {C : Finset Q} (hC : C.Nonempty)
    {Φ : Q → ℝ} {B : Q → 𝓑} {L δ ε : ℝ} {x y : Q}
    (hL : 0 ≤ L) (hLip : LipschitzBound L B)
    (hx : Φ x ≤ δ) (hy : Φ y ≤ δ)
    (hBxy : dist (B x) (B y) ≤ ε) :
    dist x y ≤
      2 * residualModulus C hC Φ δ +
        inverseModulus C B (ε + 2 * L * residualModulus C hC Φ δ) :=
  symmetric_two_output_estimate hL
    (errorBoundWitness_residualModulus C hC Φ)
    (inverseObservationBound_inverseModulus C B)
    (inverseModulus_mono C B) hLip hx hy hBxy

/-- The settled-output modulus obeys the symmetric certificate bound.  This
is the inequality whose coefficient `2` the sharpness witness shows to be
optimal. -/
theorem settledModulus_le_symmetric_bound {C : Finset Q} (hC : C.Nonempty)
    {Φ : Q → ℝ} {B : Q → 𝓑} {L δ r : ℝ}
    (hL : 0 ≤ L) (hLip : LipschitzBound L B) :
    settledModulus Φ B δ r ≤
      2 * residualModulus C hC Φ δ +
        inverseModulus C B (r + 2 * L * residualModulus C hC Φ δ) := by
  by_cases h : (settledSet Φ B δ r).Nonempty
  · rw [settledModulus_of_nonempty h]
    refine Finset.sup'_le _ _ fun p hp => ?_
    rw [mem_settledSet] at hp
    exact intrinsic_symmetric_estimate hC hL hLip hp.1 hp.2.1 hp.2.2
  · rw [settledModulus_of_empty h]
    have h1 := residualModulus_nonneg C hC Φ δ
    have h2 := inverseModulus_nonneg C B
      (r + 2 * L * residualModulus C hC Φ δ)
    linarith

/-- `cor:rate-transfer`: power-rate bounds on the intrinsic moduli transfer
to a compound output certificate.  The exponents are real (`Real.rpow`). -/
theorem rate_transfer {C : Finset Q} (hC : C.Nonempty)
    {Φ : Q → ℝ} {B : Q → 𝓑} {L a b p q δx δy ε : ℝ} {x y : Q}
    (hL : 0 ≤ L) (hLip : LipschitzBound L B) (hΦnn : ∀ z, 0 ≤ Φ z)
    (hηrate : ∀ t, 0 ≤ t → residualModulus C hC Φ t ≤ a * t ^ p)
    (hωrate : ∀ r, 0 ≤ r → inverseModulus C B r ≤ b * r ^ q)
    (hx : Φ x ≤ δx) (hy : Φ y ≤ δy)
    (hBxy : dist (B x) (B y) ≤ ε) :
    dist x y ≤
      a * (δx ^ p + δy ^ p) +
        b * (ε + L * a * (δx ^ p + δy ^ p)) ^ q := by
  have hδx : 0 ≤ δx := (hΦnn x).trans hx
  have hδy : 0 ≤ δy := (hΦnn y).trans hy
  have hε : 0 ≤ ε := dist_nonneg.trans hBxy
  have hmaster := intrinsic_two_output_estimate hC hL hLip hx hy hBxy
  have hηx := hηrate δx hδx
  have hηy := hηrate δy hδy
  have hηxnn := residualModulus_nonneg C hC Φ δx
  have hηynn := residualModulus_nonneg C hC Φ δy
  have hsum : residualModulus C hC Φ δx + residualModulus C hC Φ δy ≤
      a * (δx ^ p + δy ^ p) := by
    have := add_le_add hηx hηy
    linarith
  have hL' := mul_le_mul_of_nonneg_left hsum hL
  have harg : ε + L * (residualModulus C hC Φ δx + residualModulus C hC Φ δy)
      ≤ ε + L * a * (δx ^ p + δy ^ p) := by linarith
  have hargnn : 0 ≤ ε + L * a * (δx ^ p + δy ^ p) := by
    have h1 : 0 ≤ a * (δx ^ p + δy ^ p) := le_trans (by linarith) hsum
    have h2 := mul_nonneg hL h1
    linarith
  have hωbound :
      inverseModulus C B
          (ε + L * (residualModulus C hC Φ δx + residualModulus C hC Φ δy))
        ≤ b * (ε + L * a * (δx ^ p + δy ^ p)) ^ q :=
    le_trans (inverseModulus_mono C B harg) (hωrate _ hargnn)
  linarith

/-- `thm:sharpness` (i): at zero residual the settled-output modulus *equals*
the inverse-observation modulus, so the observation term of the symmetric
bound is exact.  The two defining maxima range over the same finite set of
pairs, which is nonempty at every `r ≥ 0`. -/
theorem settledModulus_zero_residual {C : Finset Q} (hC : C.Nonempty)
    {Φ : Q → ℝ} {B : Q → 𝓑}
    (hΦnn : ∀ x, 0 ≤ Φ x) (hCiff : ∀ x, x ∈ C ↔ Φ x = 0)
    {r : ℝ} (hr : 0 ≤ r) :
    settledModulus Φ B 0 r = inverseModulus C B r := by
  have hsets : settledSet Φ B 0 r = inverseSet C B r := by
    ext p
    rw [mem_settledSet, mem_inverseSet]
    constructor
    · rintro ⟨h1, h2, h3⟩
      exact ⟨(hCiff p.1).mpr (le_antisymm h1 (hΦnn p.1)),
        (hCiff p.2).mpr (le_antisymm h2 (hΦnn p.2)), h3⟩
    · rintro ⟨h1, h2, h3⟩
      exact ⟨le_of_eq ((hCiff p.1).mp h1), le_of_eq ((hCiff p.2).mp h2), h3⟩
  have hne : (inverseSet C B r).Nonempty := inverseSet_nonempty hC hr
  have hne' : (settledSet Φ B 0 r).Nonempty := by rw [hsets]; exact hne
  rw [settledModulus_of_nonempty hne', inverseModulus_of_nonempty hne]
  exact Finset.sup'_congr hne' hsets fun _ _ => rfl

end Intrinsic

/-! ### Injectivity boundary and finite separation radius -/

section InjectivityBoundary

variable {Q : Type u} {𝓑 : Type v}
variable [MetricSpace Q] [Fintype Q] [MetricSpace 𝓑]

omit [Fintype Q] in
/-- `prop:minimal-moduli` (iii): `ωB(0) = 0` iff `B` is injective on the
consistent set. -/
theorem inverseModulus_zero_iff_injOn {C : Finset Q} (hC : C.Nonempty)
    {B : Q → 𝓑} :
    inverseModulus C B 0 = 0 ↔ Set.InjOn B (↑C : Set Q) := by
  constructor
  · intro h0 c hc d hd hBcd
    have h1 := dist_le_inverseModulus (B := B)
      (Finset.mem_coe.mp hc) (Finset.mem_coe.mp hd)
    rw [hBcd, dist_self, h0] at h1
    exact eq_of_dist_eq_zero (le_antisymm h1 dist_nonneg)
  · intro hinj
    have hne := inverseSet_nonempty hC (B := B) le_rfl
    rw [inverseModulus_of_nonempty hne]
    refine le_antisymm (Finset.sup'_le _ _ fun p hp => ?_) ?_
    · rw [mem_inverseSet] at hp
      have hB : B p.1 = B p.2 :=
        eq_of_dist_eq_zero (le_antisymm hp.2.2 dist_nonneg)
      have hp12 : p.1 = p.2 :=
        hinj (Finset.mem_coe.mpr hp.1) (Finset.mem_coe.mpr hp.2.1) hB
      simp [hp12]
    · obtain ⟨c, hc⟩ := id hC
      have hmem : (c, c) ∈ inverseSet C B (0 : ℝ) :=
        mem_inverseSet.mpr ⟨hc, hc, by simp⟩
      exact dist_nonneg.trans
        (Finset.le_sup' (fun q : Q × Q => dist q.1 q.2) hmem)

omit [Fintype Q] in
/-- `prop:minimal-moduli` (finite separation radius): if `B` is injective on
the finite consistent set, then ωB vanishes on a whole nonempty interval
`[0, r₀)`.  When distinct consistent points exist, `r₀` is the attained
minimum of their observation separations; otherwise any positive radius
works — an explicit case split, not a default value. -/
theorem exists_separation_radius {C : Finset Q} (hC : C.Nonempty)
    {B : Q → 𝓑} (hinj : Set.InjOn B (↑C : Set Q)) :
    ∃ r₀ : ℝ, 0 < r₀ ∧
      ∀ r : ℝ, 0 ≤ r → r < r₀ → inverseModulus C B r = 0 := by
  classical
  set D := (C ×ˢ C).filter (fun p : Q × Q => p.1 ≠ p.2) with hD
  have memD : ∀ {p : Q × Q}, p ∈ D ↔ p.1 ∈ C ∧ p.2 ∈ C ∧ p.1 ≠ p.2 := by
    intro p
    rw [hD, Finset.mem_filter, Finset.mem_product]
    tauto
  have key : ∀ r : ℝ, 0 ≤ r →
      (∀ p ∈ inverseSet C B r, p.1 = p.2) → inverseModulus C B r = 0 := by
    intro r hr hdiag
    have hne := inverseSet_nonempty hC (B := B) hr
    rw [inverseModulus_of_nonempty hne]
    refine le_antisymm (Finset.sup'_le _ _ fun p hp => ?_) ?_
    · rw [hdiag p hp]; simp
    · obtain ⟨c, hc⟩ := id hC
      have hmem : (c, c) ∈ inverseSet C B r :=
        mem_inverseSet.mpr ⟨hc, hc, by simpa using hr⟩
      exact dist_nonneg.trans
        (Finset.le_sup' (fun q : Q × Q => dist q.1 q.2) hmem)
  by_cases hDne : D.Nonempty
  · obtain ⟨p₀, hp₀, hp₀eq⟩ :=
      Finset.exists_mem_eq_inf' hDne fun p : Q × Q => dist (B p.1) (B p.2)
    set r₀ := D.inf' hDne fun p : Q × Q => dist (B p.1) (B p.2) with hr₀
    have hp₀' := memD.mp hp₀
    have hr₀pos : 0 < r₀ := by
      rw [hp₀eq]
      refine dist_pos.mpr fun hEq => hp₀'.2.2 ?_
      exact hinj (Finset.mem_coe.mpr hp₀'.1) (Finset.mem_coe.mpr hp₀'.2.1) hEq
    refine ⟨r₀, hr₀pos, fun r hr hlt => key r hr fun p hp => ?_⟩
    by_contra hne'
    rw [mem_inverseSet] at hp
    have hpD : p ∈ D := memD.mpr ⟨hp.1, hp.2.1, hne'⟩
    have hle : r₀ ≤ dist (B p.1) (B p.2) :=
      Finset.inf'_le (fun p : Q × Q => dist (B p.1) (B p.2)) hpD
    exact absurd (hle.trans hp.2.2) (not_le.mpr hlt)
  · refine ⟨1, one_pos, fun r hr _ => key r hr fun p hp => ?_⟩
    by_contra hne'
    rw [mem_inverseSet] at hp
    exact hDne ⟨p, memD.mpr ⟨hp.1, hp.2.1, hne'⟩⟩

end InjectivityBoundary

/-! ### The three-point sharpness witness

Three points of ℝ: `{-1, 0, 1}` with the induced metric, consistent set
`{0}`, residual `Φ = |·|`, constant observation `B ≡ 0`.  All three moduli
are *computed* from their defining maxima:
`ηΦ(1) = 1`, `ωB(0) = 0`, `Θ(1, 0) = 2`, so `Θ(1, 0) = 2·ηΦ(1) + ωB(0)`
and no coefficient `a < 2` can replace the `2`. -/

namespace SharpnessWitness

noncomputable def carrier : Finset ℝ := {-1, 0, 1}

noncomputable def pNeg : ↥carrier := ⟨-1, by simp [carrier]⟩
noncomputable def pZero : ↥carrier := ⟨0, by simp [carrier]⟩
noncomputable def pPos : ↥carrier := ⟨1, by simp [carrier]⟩

noncomputable def Cw : Finset ↥carrier := {pZero}
noncomputable def Φw : ↥carrier → ℝ := fun z => |z.val|
def Bw : ↥carrier → ℝ := fun _ => 0

lemma Cw_nonempty : Cw.Nonempty := ⟨pZero, Finset.mem_singleton_self _⟩

lemma Φw_nonneg : ∀ z, 0 ≤ Φw z := fun z => abs_nonneg z.val

lemma mem_Cw_iff : ∀ z : ↥carrier, z ∈ Cw ↔ Φw z = 0 := by
  intro z
  constructor
  · intro hz
    rw [Cw, Finset.mem_singleton] at hz
    subst hz
    simp [Φw, pZero]
  · intro hz
    have hv : z.val = 0 := abs_eq_zero.mp hz
    rw [Cw, Finset.mem_singleton]
    exact Subtype.ext hv

lemma point_cases (z : ↥carrier) : z = pNeg ∨ z = pZero ∨ z = pPos := by
  obtain ⟨v, hv⟩ := z
  have hv' : v = -1 ∨ v = 0 ∨ v = 1 := by simpa [carrier] using hv
  rcases hv' with h | h | h
  · exact Or.inl (Subtype.ext h)
  · exact Or.inr (Or.inl (Subtype.ext h))
  · exact Or.inr (Or.inr (Subtype.ext h))

lemma dist_val (z w : ↥carrier) : dist z w = |z.val - w.val| := by
  rw [Subtype.dist_eq, Real.dist_eq]

lemma distToFinset_Cw (z : ↥carrier) :
    distToFinset Cw Cw_nonempty z = dist z pZero := by
  simp [distToFinset, Cw]

/-- Computed: `ηΦ(1) = 1`, attained at the state `1`. -/
theorem residual_at_one : residualModulus Cw Cw_nonempty Φw 1 = 1 := by
  have hmem : pPos ∈ residualSet Φw 1 :=
    mem_residualSet.mpr (by norm_num [Φw, pPos])
  have hne : (residualSet Φw 1).Nonempty := ⟨pPos, hmem⟩
  rw [residualModulus_of_nonempty Cw_nonempty hne]
  apply le_antisymm
  · refine Finset.sup'_le _ _ fun z _ => ?_
    rw [distToFinset_Cw, dist_val]
    rcases point_cases z with h | h | h <;> subst h <;>
      norm_num [pNeg, pZero, pPos]
  · have h1 : distToFinset Cw Cw_nonempty pPos = 1 := by
      rw [distToFinset_Cw, dist_val]
      norm_num [pPos, pZero]
    calc (1 : ℝ) = distToFinset Cw Cw_nonempty pPos := h1.symm
      _ ≤ _ := Finset.le_sup' _ hmem

/-- Computed: `ωB(0) = 0` (the consistent set is a singleton, so `B` is
injective on it). -/
theorem inverse_at_zero : inverseModulus Cw Bw 0 = 0 := by
  refine (inverseModulus_zero_iff_injOn Cw_nonempty).mpr ?_
  intro u hu v hv _
  have hu' : u = pZero := by
    have := Finset.mem_coe.mp hu
    simpa [Cw] using this
  have hv' : v = pZero := by
    have := Finset.mem_coe.mp hv
    simpa [Cw] using this
  rw [hu', hv']

/-- Computed: `Θ(1, 0) = 2`, attained by the pair `(1, -1)`. -/
theorem settled_at_one_zero : settledModulus Φw Bw 1 0 = 2 := by
  have hmem : (pPos, pNeg) ∈ settledSet Φw Bw 1 0 := by
    rw [mem_settledSet]
    norm_num [Φw, Bw, pPos, pNeg]
  have hne : (settledSet Φw Bw 1 0).Nonempty := ⟨(pPos, pNeg), hmem⟩
  rw [settledModulus_of_nonempty hne]
  apply le_antisymm
  · refine Finset.sup'_le _ _ fun p _ => ?_
    rcases point_cases p.1 with h1 | h1 | h1 <;>
      rcases point_cases p.2 with h2 | h2 | h2 <;>
        rw [h1, h2, dist_val] <;> norm_num [pNeg, pZero, pPos]
  · have h2 : dist pPos pNeg = 2 := by
      rw [dist_val]; norm_num [pPos, pNeg]
    calc (2 : ℝ) = dist pPos pNeg := h2.symm
      _ ≤ _ := Finset.le_sup' (fun p : ↥carrier × ↥carrier => dist p.1 p.2) hmem

/-- `thm:sharpness` (ii): on the three-point witness the symmetric bound is
an equality, so its coefficient `2` cannot be lowered anywhere. -/
theorem symmetric_bound_attained :
    settledModulus Φw Bw 1 0 =
      2 * residualModulus Cw Cw_nonempty Φw 1 + inverseModulus Cw Bw 0 := by
  rw [residual_at_one, inverse_at_zero, settled_at_one_zero]
  norm_num

/-- `thm:sharpness` (ii), quantitative form: every coefficient `a < 2` fails
on the computed witness values. -/
theorem coefficient_two_sharp {a : ℝ} (ha : a < 2) :
    a * residualModulus Cw Cw_nonempty Φw 1 + inverseModulus Cw Bw 0 <
      settledModulus Φw Bw 1 0 := by
  rw [residual_at_one, inverse_at_zero, settled_at_one_zero]
  linarith

end SharpnessWitness

/-! ### Negative control

Two points of ℝ: `{0, 1}`, consistent set `{0}`, `Φ = |·|`, `B ≡ 0`.
Computed: `ηΦ(1) = 1`, `ωB(0) = 0`, but `Θ(1, 0) = 1 < 2 = 2·ηΦ(1) + ωB(0)`.
The symmetric certificate bound is *strictly slack* here, so the sharpness
statement above has content: attainment is a property of the three-point
witness, not a tautology of the definitions. -/

namespace NegativeControl

noncomputable def carrier : Finset ℝ := {0, 1}

noncomputable def pZero : ↥carrier := ⟨0, by simp [carrier]⟩
noncomputable def pOne : ↥carrier := ⟨1, by simp [carrier]⟩

noncomputable def Cn : Finset ↥carrier := {pZero}
noncomputable def Φn : ↥carrier → ℝ := fun z => |z.val|
def Bn : ↥carrier → ℝ := fun _ => 0

lemma Cn_nonempty : Cn.Nonempty := ⟨pZero, Finset.mem_singleton_self _⟩

lemma point_cases (z : ↥carrier) : z = pZero ∨ z = pOne := by
  obtain ⟨v, hv⟩ := z
  have hv' : v = 0 ∨ v = 1 := by simpa [carrier] using hv
  rcases hv' with h | h
  · exact Or.inl (Subtype.ext h)
  · exact Or.inr (Subtype.ext h)

lemma dist_val (z w : ↥carrier) : dist z w = |z.val - w.val| := by
  rw [Subtype.dist_eq, Real.dist_eq]

lemma distToFinset_Cn (z : ↥carrier) :
    distToFinset Cn Cn_nonempty z = dist z pZero := by
  simp [distToFinset, Cn]

/-- Computed: `ηΦ(1) = 1`, attained at the state `1`. -/
theorem residual_at_one : residualModulus Cn Cn_nonempty Φn 1 = 1 := by
  have hmem : pOne ∈ residualSet Φn 1 :=
    mem_residualSet.mpr (by norm_num [Φn, pOne])
  have hne : (residualSet Φn 1).Nonempty := ⟨pOne, hmem⟩
  rw [residualModulus_of_nonempty Cn_nonempty hne]
  apply le_antisymm
  · refine Finset.sup'_le _ _ fun z _ => ?_
    rw [distToFinset_Cn, dist_val]
    rcases point_cases z with h | h <;> subst h <;> norm_num [pZero, pOne]
  · have h1 : distToFinset Cn Cn_nonempty pOne = 1 := by
      rw [distToFinset_Cn, dist_val]
      norm_num [pOne, pZero]
    calc (1 : ℝ) = distToFinset Cn Cn_nonempty pOne := h1.symm
      _ ≤ _ := Finset.le_sup' _ hmem

/-- Computed: `ωB(0) = 0`. -/
theorem inverse_at_zero : inverseModulus Cn Bn 0 = 0 := by
  refine (inverseModulus_zero_iff_injOn Cn_nonempty).mpr ?_
  intro u hu v hv _
  have hu' : u = pZero := by
    have := Finset.mem_coe.mp hu
    simpa [Cn] using this
  have hv' : v = pZero := by
    have := Finset.mem_coe.mp hv
    simpa [Cn] using this
  rw [hu', hv']

/-- Computed: `Θ(1, 0) = 1` — strictly below the certificate bound `2`. -/
theorem settled_at_one_zero : settledModulus Φn Bn 1 0 = 1 := by
  have hmem : (pOne, pZero) ∈ settledSet Φn Bn 1 0 := by
    rw [mem_settledSet]
    norm_num [Φn, Bn, pOne, pZero]
  have hne : (settledSet Φn Bn 1 0).Nonempty := ⟨(pOne, pZero), hmem⟩
  rw [settledModulus_of_nonempty hne]
  apply le_antisymm
  · refine Finset.sup'_le _ _ fun p _ => ?_
    rcases point_cases p.1 with h1 | h1 <;>
      rcases point_cases p.2 with h2 | h2 <;>
        rw [h1, h2, dist_val] <;> norm_num [pZero, pOne]
  · have h1 : dist pOne pZero = 1 := by
      rw [dist_val]; norm_num [pOne, pZero]
    calc (1 : ℝ) = dist pOne pZero := h1.symm
      _ ≤ _ := Finset.le_sup' (fun p : ↥carrier × ↥carrier => dist p.1 p.2) hmem

/-- Negative control: the symmetric bound is NOT attained on the two-point
system.  Sharpness of the coefficient `2` is therefore a property of the
three-point witness, not an artifact of the definitions. -/
theorem bound_not_attained :
    settledModulus Φn Bn 1 0 <
      2 * residualModulus Cn Cn_nonempty Φn 1 + inverseModulus Cn Bn 0 := by
  rw [residual_at_one, inverse_at_zero, settled_at_one_zero]
  norm_num

end NegativeControl

end ObservableNormalForms
