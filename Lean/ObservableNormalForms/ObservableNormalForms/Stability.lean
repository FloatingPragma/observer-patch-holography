import Mathlib

/-!
# Quantitative observable stability

This module formalizes the proof-bearing core of the manuscript's principal
stability estimate and its schedule corollary.  The hypotheses
are stated as reusable certificates rather than as finite maxima.  On a
finite metric space the manuscript's intrinsic `η` and `ω` supply exactly
these certificates.
-/

namespace ObservableNormalForms

universe u v w

section Certificates

variable {Q : Type u} {𝓑 : Type v}
variable [PseudoMetricSpace Q] [PseudoMetricSpace 𝓑]

/-- A residual-to-consistency certificate with attained witnesses.  For the
finite systems in the manuscript, choose a nearest consistent state and use
the intrinsic residual modulus. -/
def ErrorBoundWitness (C : Set Q) (Φ : Q → ℝ) (η : ℝ → ℝ) : Prop :=
  ∀ ⦃x : Q⦄ ⦃δ : ℝ⦄, Φ x ≤ δ →
    ∃ c : Q, c ∈ C ∧ dist x c ≤ η δ

/-- An inverse-observation certificate on the consistent subset. -/
def InverseObservationBound (C : Set Q) (B : Q → 𝓑) (ω : ℝ → ℝ) : Prop :=
  ∀ ⦃c d : Q⦄, c ∈ C → d ∈ C →
    dist c d ≤ ω (dist (B c) (B d))

/-- An explicit real-valued Lipschitz bound.  This avoids importing any
nonnegative-coefficient coercions into the theorem statement. -/
def LipschitzBound (L : ℝ) (B : Q → 𝓑) : Prop :=
  ∀ x y : Q, dist (B x) (B y) ≤ L * dist x y

/-- Heterogeneous two-output estimate.  This is a
complete proof of the load-bearing triangle argument, not a restatement of
its conclusion as an assumption. -/
theorem heterogeneous_two_output_estimate
    {C : Set Q} {Φ : Q → ℝ} {B : Q → 𝓑}
    {η ω : ℝ → ℝ} {L δx δy ε : ℝ} {x y : Q}
    (hL : 0 ≤ L)
    (hη : ErrorBoundWitness C Φ η)
    (hω : InverseObservationBound C B ω)
    (hωmono : Monotone ω)
    (hLip : LipschitzBound L B)
    (hx : Φ x ≤ δx)
    (hy : Φ y ≤ δy)
    (hBxy : dist (B x) (B y) ≤ ε) :
    dist x y ≤
      η δx + η δy + ω (ε + L * (η δx + η δy)) := by
  rcases hη hx with ⟨cx, hcxC, hcx⟩
  rcases hη hy with ⟨cy, hcyC, hcy⟩
  have hObs : dist (B cx) (B cy) ≤ ε + L * (η δx + η δy) := by
    calc
      dist (B cx) (B cy)
          ≤ dist (B cx) (B x) + dist (B x) (B cy) :=
            dist_triangle _ _ _
      _ ≤ dist (B cx) (B x) +
            (dist (B x) (B y) + dist (B y) (B cy)) := by
            gcongr
            exact dist_triangle _ _ _
      _ ≤ L * dist cx x + (ε + L * dist y cy) := by
            exact add_le_add (hLip cx x)
              (add_le_add hBxy (hLip y cy))
      _ ≤ L * η δx + (ε + L * η δy) := by
            exact add_le_add
              (mul_le_mul_of_nonneg_left (by simpa [dist_comm] using hcx) hL)
              (add_le_add_right (mul_le_mul_of_nonneg_left hcy hL) ε)
      _ = ε + L * (η δx + η δy) := by ring
  have hMiddle : dist cx cy ≤ ω (ε + L * (η δx + η δy)) :=
    (hω hcxC hcyC).trans (hωmono hObs)
  calc
    dist x y ≤ dist x cx + dist cx y := dist_triangle _ _ _
    _ ≤ dist x cx + (dist cx cy + dist cy y) := by
      gcongr
      exact dist_triangle _ _ _
    _ ≤ η δx +
        (ω (ε + L * (η δx + η δy)) + η δy) := by
      exact add_le_add hcx
        (add_le_add hMiddle (by simpa [dist_comm] using hcy))
    _ = η δx + η δy + ω (ε + L * (η δx + η δy)) := by ring

/-- The equal-residual specialization of the heterogeneous estimate. -/
theorem symmetric_two_output_estimate
    {C : Set Q} {Φ : Q → ℝ} {B : Q → 𝓑}
    {η ω : ℝ → ℝ} {L δ ε : ℝ} {x y : Q}
    (hL : 0 ≤ L)
    (hη : ErrorBoundWitness C Φ η)
    (hω : InverseObservationBound C B ω)
    (hωmono : Monotone ω)
    (hLip : LipschitzBound L B)
    (hx : Φ x ≤ δ)
    (hy : Φ y ≤ δ)
    (hBxy : dist (B x) (B y) ≤ ε) :
    dist x y ≤ 2 * η δ + ω (ε + 2 * L * η δ) := by
  have h := heterogeneous_two_output_estimate hL hη hω hωmono hLip hx hy hBxy
  simpa [two_mul, mul_add, add_mul, add_assoc] using h

/-- Approximate-schedule-independence corollary: two approximately settled schedules with a
common initial observation have close endpoints.  No relation between their
rewrite paths is assumed. -/
theorem approximate_schedule_independence
    {C : Set Q} {Φ : Q → ℝ} {B : Q → 𝓑}
    {η ω : ℝ → ℝ} {L δ₁ δ₂ β₁ β₂ : ℝ}
    {initial z₁ z₂ : Q}
    (hL : 0 ≤ L)
    (hη : ErrorBoundWitness C Φ η)
    (hω : InverseObservationBound C B ω)
    (hωmono : Monotone ω)
    (hLip : LipschitzBound L B)
    (hz₁ : Φ z₁ ≤ δ₁)
    (hz₂ : Φ z₂ ≤ δ₂)
    (hobs₁ : dist (B z₁) (B initial) ≤ β₁)
    (hobs₂ : dist (B z₂) (B initial) ≤ β₂) :
    dist z₁ z₂ ≤
      η δ₁ + η δ₂ +
        ω (β₁ + β₂ + L * (η δ₁ + η δ₂)) := by
  have hEndpoints : dist (B z₁) (B z₂) ≤ β₁ + β₂ := by
    calc
      dist (B z₁) (B z₂)
          ≤ dist (B z₁) (B initial) + dist (B initial) (B z₂) :=
            dist_triangle _ _ _
      _ ≤ β₁ + β₂ := by
        gcongr
        simpa [dist_comm] using hobs₂
  exact heterogeneous_two_output_estimate
    hL hη hω hωmono hLip hz₁ hz₂ hEndpoints

/-- Observation enrichment cannot worsen any monotone inverse-observation
certificate when the enriched metric dominates the old observation metric.
This is the certificate-level content of the sensor-enrichment corollary. -/
theorem inverse_bound_of_sensor_enrichment
    {𝓓 : Type w} [PseudoMetricSpace 𝓓]
    {C : Set Q} {B : Q → 𝓑} {Bplus : Q → 𝓓} {ω : ℝ → ℝ}
    (hω : InverseObservationBound C B ω)
    (hωmono : Monotone ω)
    (hdominates : ∀ c d : Q,
      dist (B c) (B d) ≤ dist (Bplus c) (Bplus d)) :
    InverseObservationBound C Bplus ω := by
  intro c d hc hd
  exact (hω hc hd).trans (hωmono (hdominates c d))

end Certificates

end ObservableNormalForms
