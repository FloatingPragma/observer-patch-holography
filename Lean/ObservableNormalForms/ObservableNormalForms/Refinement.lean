import ObservableNormalForms.Stability

/-!
# Naturality and refinement certificates

This module formalizes the one-step approximate-naturality theorem and the
metric telescope/projective-comparison core of the manuscript's refinement
results.  It remains generic and imports only the preceding mathematical
certificate layer.
-/

namespace ObservableNormalForms

universe uₛ uᵣ vₛ vᵣ

section OneStep

variable {Qₛ : Type uₛ} {Qᵣ : Type uᵣ}
variable {𝓑ₛ : Type vₛ} {𝓑ᵣ : Type vᵣ}

/-- One-step approximate naturality.  Semantic observation and consistency defects
generate a bound on the failure of naturality of exact normalizers. -/
theorem one_step_approximate_naturality
    [PseudoMetricSpace Qᵣ] [PseudoMetricSpace 𝓑ᵣ]
    {Cₛ : Set Qₛ} {Cᵣ : Set Qᵣ}
    {Bₛ : Qₛ → 𝓑ₛ} {Bᵣ : Qᵣ → 𝓑ᵣ}
    {Φᵣ : Qᵣ → ℝ} {ηᵣ ωᵣ : ℝ → ℝ}
    {Nₛ : Qₛ → Qₛ} {Nᵣ : Qᵣ → Qᵣ}
    {f : Qₛ → Qᵣ} {h : 𝓑ₛ → 𝓑ᵣ}
    {Lᵣ δC εB : ℝ}
    (hL : 0 ≤ Lᵣ)
    (hη : ErrorBoundWitness Cᵣ Φᵣ ηᵣ)
    (hω : InverseObservationBound Cᵣ Bᵣ ωᵣ)
    (hωmono : Monotone ωᵣ)
    (hLip : LipschitzBound Lᵣ Bᵣ)
    (hNₛC : ∀ q : Qₛ, Nₛ q ∈ Cₛ)
    (hNₛobs : ∀ q : Qₛ, Bₛ (Nₛ q) = Bₛ q)
    (hNᵣC : ∀ q : Qᵣ, Nᵣ q ∈ Cᵣ)
    (hNᵣobs : ∀ q : Qᵣ, Bᵣ (Nᵣ q) = Bᵣ q)
    (hConsistencyDefect : ∀ c ∈ Cₛ, Φᵣ (f c) ≤ δC)
    (hObservationDefect : ∀ q : Qₛ,
      dist (Bᵣ (f q)) (h (Bₛ q)) ≤ εB)
    (q : Qₛ) :
    dist (f (Nₛ q)) (Nᵣ (f q)) ≤
      ηᵣ δC + ωᵣ (2 * εB + Lᵣ * ηᵣ δC) := by
  rcases hη (hConsistencyDefect (Nₛ q) (hNₛC q)) with
    ⟨c, hcC, hclose⟩
  have hYObs : dist (Bᵣ (f (Nₛ q))) (h (Bₛ q)) ≤ εB := by
    simpa [hNₛobs q] using hObservationDefect (Nₛ q)
  have hQObs : dist (h (Bₛ q)) (Bᵣ (f q)) ≤ εB := by
    simpa [dist_comm] using hObservationDefect q
  have hObs :
      dist (Bᵣ c) (Bᵣ (Nᵣ (f q))) ≤
        2 * εB + Lᵣ * ηᵣ δC := by
    rw [hNᵣobs]
    calc
      dist (Bᵣ c) (Bᵣ (f q))
          ≤ dist (Bᵣ c) (Bᵣ (f (Nₛ q))) +
              dist (Bᵣ (f (Nₛ q))) (Bᵣ (f q)) :=
            dist_triangle _ _ _
      _ ≤ dist (Bᵣ c) (Bᵣ (f (Nₛ q))) +
            (dist (Bᵣ (f (Nₛ q))) (h (Bₛ q)) +
              dist (h (Bₛ q)) (Bᵣ (f q))) := by
            gcongr
            exact dist_triangle _ _ _
      _ ≤ Lᵣ * dist c (f (Nₛ q)) + (εB + εB) := by
            exact add_le_add (hLip c (f (Nₛ q)))
              (add_le_add hYObs hQObs)
      _ ≤ Lᵣ * ηᵣ δC + (εB + εB) := by
            exact add_le_add
              (mul_le_mul_of_nonneg_left
                (by simpa [dist_comm] using hclose) hL)
              le_rfl
      _ = 2 * εB + Lᵣ * ηᵣ δC := by ring
  have hConsistentDistance :
      dist c (Nᵣ (f q)) ≤ ωᵣ (2 * εB + Lᵣ * ηᵣ δC) :=
    (hω hcC (hNᵣC (f q))).trans (hωmono hObs)
  calc
    dist (f (Nₛ q)) (Nᵣ (f q))
        ≤ dist (f (Nₛ q)) c + dist c (Nᵣ (f q)) :=
          dist_triangle _ _ _
    _ ≤ ηᵣ δC + ωᵣ (2 * εB + Lᵣ * ηᵣ δC) := by
          exact add_le_add hclose hConsistentDistance

/-- Exact-naturality corollary.  Exact preservation of consistency and
observations makes normalizers natural whenever coarse consistent fibers are
observable-determined. -/
theorem exact_naturality_from_uniqueness
    {Cₛ : Set Qₛ} {Cᵣ : Set Qᵣ}
    {Bₛ : Qₛ → 𝓑ₛ} {Bᵣ : Qᵣ → 𝓑ᵣ}
    {Nₛ : Qₛ → Qₛ} {Nᵣ : Qᵣ → Qᵣ}
    {f : Qₛ → Qᵣ} {h : 𝓑ₛ → 𝓑ᵣ}
    (hDetermined : Set.InjOn Bᵣ Cᵣ)
    (hNₛC : ∀ q : Qₛ, Nₛ q ∈ Cₛ)
    (hNₛobs : ∀ q : Qₛ, Bₛ (Nₛ q) = Bₛ q)
    (hNᵣC : ∀ q : Qᵣ, Nᵣ q ∈ Cᵣ)
    (hNᵣobs : ∀ q : Qᵣ, Bᵣ (Nᵣ q) = Bᵣ q)
    (hfC : ∀ c ∈ Cₛ, f c ∈ Cᵣ)
    (hcommute : ∀ q : Qₛ, Bᵣ (f q) = h (Bₛ q)) :
    f ∘ Nₛ = Nᵣ ∘ f := by
  funext q
  apply hDetermined (hfC (Nₛ q) (hNₛC q)) (hNᵣC (f q))
  rw [hcommute, hNₛobs, ← hcommute, hNᵣobs]

end OneStep

section Telescope

variable {Q : Type uᵣ} [PseudoMetricSpace Q]

/-- The metric telescope used in the refinement-tower theorem.  The paper obtains
the step bound by applying a restriction-map Lipschitz constant to each
one-step naturality defect.  This theorem performs the nontrivial arbitrary-
depth accumulation and introduces no positivity assumptions on the supplied
bounds beyond what each step proof itself entails. -/
theorem telescoping_refinement_bound
    (p : ℕ → Q) (K a : ℕ → ℝ)
    (hstep : ∀ j : ℕ, dist (p (j + 1)) (p j) ≤ K j * a j)
    (n k : ℕ) :
    dist (p (n + k)) (p n) ≤
      ∑ j ∈ Finset.range k, K (n + j) * a (n + j) := by
  induction k with
  | zero => simp
  | succ k ih =>
      calc
        dist (p (n + (k + 1))) (p n)
            ≤ dist (p (n + (k + 1))) (p (n + k)) +
                dist (p (n + k)) (p n) := dist_triangle _ _ _
        _ ≤ K (n + k) * a (n + k) +
              ∑ j ∈ Finset.range k, K (n + j) * a (n + j) := by
              gcongr
              simpa [Nat.add_assoc] using hstep (n + k)
        _ = ∑ j ∈ Finset.range (k + 1),
              K (n + j) * a (n + j) := by
              rw [Finset.sum_range_succ]
              ring

/-- Metric core of the fine-to-coarse solver certificate, after the exact
tower has provided its telescoping receipt.  The theorem proves how a solver
receipt and the exact tower defect combine under a Lipschitz restriction map. -/
theorem projective_implementation_bound_from_tower_receipt
    {Qfine Qcoarse : Type*}
    [PseudoMetricSpace Qfine] [PseudoMetricSpace Qcoarse]
    {ρ : Qfine → Qcoarse}
    {approx exactFine : Qfine} {exactCoarse : Qcoarse}
    {K e tail : ℝ}
    (hK : 0 ≤ K)
    (hρ : LipschitzBound K ρ)
    (hsolver : dist approx exactFine ≤ e)
    (htower : dist (ρ exactFine) exactCoarse ≤ tail) :
    dist (ρ approx) exactCoarse ≤ K * e + tail := by
  calc
    dist (ρ approx) exactCoarse
        ≤ dist (ρ approx) (ρ exactFine) +
            dist (ρ exactFine) exactCoarse := dist_triangle _ _ _
    _ ≤ K * dist approx exactFine + tail :=
      add_le_add (hρ approx exactFine) htower
    _ ≤ K * e + tail :=
      add_le_add (mul_le_mul_of_nonneg_left hsolver hK) le_rfl

/-- Metric core of the manuscript's same-level implementation-agreement
proposition.  Both implementations are compared with one exact state before
the common restriction is applied, so no refinement-tail term appears. -/
theorem same_level_implementation_agreement
    {Qfine Qcoarse : Type*}
    [PseudoMetricSpace Qfine] [PseudoMetricSpace Qcoarse]
    {ρ : Qfine → Qcoarse}
    {approx₁ approx₂ exact : Qfine}
    {K e₁ e₂ : ℝ}
    (hK : 0 ≤ K)
    (hρ : LipschitzBound K ρ)
    (h₁ : dist approx₁ exact ≤ e₁)
    (h₂ : dist approx₂ exact ≤ e₂) :
    dist (ρ approx₁) (ρ approx₂) ≤ K * (e₁ + e₂) := by
  calc
    dist (ρ approx₁) (ρ approx₂)
        ≤ dist (ρ approx₁) (ρ exact) +
            dist (ρ exact) (ρ approx₂) := dist_triangle _ _ _
    _ ≤ K * dist approx₁ exact + K * dist exact approx₂ :=
      add_le_add (hρ approx₁ exact) (hρ exact approx₂)
    _ ≤ K * e₁ + K * e₂ :=
      add_le_add (mul_le_mul_of_nonneg_left h₁ hK)
        (mul_le_mul_of_nonneg_left (by simpa [dist_comm] using h₂) hK)
    _ = K * (e₁ + e₂) := by ring

/-- Metric core of the anchored cross-level comparison theorem.  The five
points are the restricted approximate output, restricted exact output,
first anchored normal form, second anchored normal form, second restricted
exact output, and second restricted approximate output.  All Lipschitz and
modulus amplification has already been discharged into the five scalar
receipts. -/
theorem anchored_cross_level_metric_core
    {X : Type*} [PseudoMetricSpace X]
    {approxM exactM anchorM anchorL exactL approxL : X}
    {solverM pathM anchorMismatch pathL solverL : ℝ}
    (hsolverM : dist approxM exactM ≤ solverM)
    (hpathM : dist exactM anchorM ≤ pathM)
    (hanchor : dist anchorM anchorL ≤ anchorMismatch)
    (hpathL : dist anchorL exactL ≤ pathL)
    (hsolverL : dist exactL approxL ≤ solverL) :
    dist approxM approxL ≤
      solverM + pathM + anchorMismatch + pathL + solverL := by
  have h₀ := dist_triangle approxM exactM approxL
  have h₁ := dist_triangle exactM anchorM approxL
  have h₂ := dist_triangle anchorM anchorL approxL
  have h₃ := dist_triangle anchorL exactL approxL
  linarith

/-- Metric core of the nested-compatible-level corollary.  Compatibility
identifies the anchor observations, eliminating both the first path segment
and the inverse-observation mismatch term. -/
theorem nested_compatible_levels_metric_core
    {X : Type*} [PseudoMetricSpace X]
    {approxCoarse exactCoarse exactFine approxFine : X}
    {coarseReceipt pathSegment fineReceipt : ℝ}
    (hcoarse : dist approxCoarse exactCoarse ≤ coarseReceipt)
    (hpath : dist exactCoarse exactFine ≤ pathSegment)
    (hfine : dist exactFine approxFine ≤ fineReceipt) :
    dist approxCoarse approxFine ≤
      coarseReceipt + pathSegment + fineReceipt := by
  have h₀ := dist_triangle approxCoarse exactCoarse approxFine
  have h₁ := dist_triangle exactCoarse exactFine approxFine
  linarith

/-- Supporting precursor retained for comparison with the earlier draft:
two implementations are compared with a common exact coarse point after each
has paid the same tower receipt.  The current manuscript uses the sharper
same-level and anchored/nested statements above. -/
theorem two_implementations_bound_from_tower_receipts
    {Qfine Qcoarse : Type*}
    [PseudoMetricSpace Qfine] [PseudoMetricSpace Qcoarse]
    {ρ : Qfine → Qcoarse}
    {approx₁ approx₂ exactFine : Qfine} {exactCoarse : Qcoarse}
    {K e₁ e₂ tail : ℝ}
    (hK : 0 ≤ K)
    (hρ : LipschitzBound K ρ)
    (hsolver₁ : dist approx₁ exactFine ≤ e₁)
    (hsolver₂ : dist approx₂ exactFine ≤ e₂)
    (htower : dist (ρ exactFine) exactCoarse ≤ tail) :
    dist (ρ approx₁) (ρ approx₂) ≤
      K * (e₁ + e₂) + 2 * tail := by
  have h₁ := projective_implementation_bound_from_tower_receipt
    hK hρ hsolver₁ htower
  have h₂ := projective_implementation_bound_from_tower_receipt
    hK hρ hsolver₂ htower
  calc
    dist (ρ approx₁) (ρ approx₂)
        ≤ dist (ρ approx₁) exactCoarse +
            dist exactCoarse (ρ approx₂) := dist_triangle _ _ _
    _ ≤ (K * e₁ + tail) + (K * e₂ + tail) :=
      add_le_add h₁ (by simpa [dist_comm] using h₂)
    _ = K * (e₁ + e₂) + 2 * tail := by ring

end Telescope

end ObservableNormalForms
