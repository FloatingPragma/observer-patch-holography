import Mathlib

/-!
# Coupling-Algebra Part A: Bridge Equivalence

The algebraic equivalence at the core of the OPH coupling bridge: for real
`P > 0`, `alphaU > 0`, `N > π`,

```text
N = π · exp(6π / (P · alphaU))   ↔   (24π) / (alphaU · log(N/π)) = 4 · P.
```

The left-hand side is the count form of the bridge; the right-hand side is
the tick-projection form. The literal `6` in the exponent is not primitive:
it is the ratio `m_rep / β_EW` with `m_rep = 24` and `β_EW = 4`, and the
generalised statement `bridge_equivalence_beta` keeps that dependence
explicit by proving the equivalence for every `β > 0` with the exponent
written as `(m_rep / β) · π / (P · alphaU)`.

**Scope warning.** This module formalises the ALGEBRAIC layer of the
coupling theorem only: it is a log/exp manipulation over `ℝ` and carries no
physical-derivation content. The physical identities I1/I2 (which give the
two sides their physical readings) are outside the formalised set. Numeric
enclosures stay in the Python certificates; nothing here is floating-point.
-/

namespace OPH.BridgeEquivalence

/-- Replica multiplicity `m_rep = 24`. Declared as a named constant so the
    provenance of the literal `24` in `tickProjection` is explicit. -/
def mRep : ℝ := 24

/-- Electroweak beta slot `β_EW = 4`. Declared as a named constant so the
    provenance of the literal `4` on the tick side is explicit. -/
def betaEW : ℝ := 4

/-- The literal `6` in the bridge exponent is the ratio `m_rep / β_EW`. -/
theorem six_eq_mRep_div_betaEW : (6 : ℝ) = mRep / betaEW := by
  norm_num [mRep, betaEW]

/-- The literal `24` decomposes as `6 · β_EW`. -/
theorem mRep_eq_six_mul_betaEW : mRep = 6 * betaEW := by
  norm_num [mRep, betaEW]

/-- The tick projection `(m_rep · π) / (alphaU · log(N/π))` with
    `m_rep = 24`. -/
noncomputable def tickProjection (alphaU N : ℝ) : ℝ :=
  (mRep * Real.pi) / (alphaU * Real.log (N / Real.pi))

/-- **Bridge equivalence, general beta form.** For every `β > 0` the count
    form with exponent `(m_rep / β) · π / (P · alphaU)` is equivalent to the
    tick form `tickProjection alphaU N = β · P`. The dependence of the
    exponent's numerator on `m_rep / β` (specialising to `6 = 24 / 4` at
    `β = β_EW`) is explicit in the statement. -/
theorem bridge_equivalence_beta {beta P alphaU N : ℝ}
    (hbeta : 0 < beta) (hP : 0 < P) (ha : 0 < alphaU) (hN : Real.pi < N) :
    N = Real.pi * Real.exp ((mRep / beta) * Real.pi / (P * alphaU)) ↔
      tickProjection alphaU N = beta * P := by
  have hpi : (0 : ℝ) < Real.pi := Real.pi_pos
  have hNpos : (0 : ℝ) < N := hpi.trans hN
  have hL : 0 < Real.log (N / Real.pi) := Real.log_pos ((one_lt_div hpi).mpr hN)
  have hmRep : (0 : ℝ) < mRep := by norm_num [mRep]
  unfold tickProjection
  constructor
  · intro hEq
    have hlog : Real.log (N / Real.pi) = (mRep / beta) * Real.pi / (P * alphaU) := by
      rw [hEq, mul_div_cancel_left₀ _ hpi.ne', Real.log_exp]
    rw [hlog]
    field_simp
  · intro hTick
    have hlog : Real.log (N / Real.pi) = (mRep / beta) * Real.pi / (P * alphaU) := by
      have hden : alphaU * Real.log (N / Real.pi) ≠ 0 := (mul_pos ha hL).ne'
      field_simp at hTick ⊢
      linear_combination -hTick
    calc N = Real.pi * (N / Real.pi) := by field_simp
      _ = Real.pi * Real.exp (Real.log (N / Real.pi)) := by
          rw [Real.exp_log (div_pos hNpos hpi)]
      _ = Real.pi * Real.exp ((mRep / beta) * Real.pi / (P * alphaU)) := by
          rw [hlog]

/-- **Bridge equivalence, literal form.** For `P > 0`, `alphaU > 0`,
    `N > π`:

    `N = π · exp(6π / (P · alphaU)) ↔ (24π) / (alphaU · log(N/π)) = 4 · P`.

    Specialisation of `bridge_equivalence_beta` at `β = β_EW = 4`, where
    `m_rep / β_EW = 24 / 4 = 6`. -/
theorem bridge_equivalence {P alphaU N : ℝ}
    (hP : 0 < P) (ha : 0 < alphaU) (hN : Real.pi < N) :
    N = Real.pi * Real.exp (6 * Real.pi / (P * alphaU)) ↔
      (24 * Real.pi) / (alphaU * Real.log (N / Real.pi)) = 4 * P := by
  have h := bridge_equivalence_beta (beta := betaEW)
    (by norm_num [betaEW]) hP ha hN
  rw [show mRep / betaEW = (6 : ℝ) by norm_num [mRep, betaEW]] at h
  rw [tickProjection, mRep, betaEW] at h
  exact h

/-- **Bridge equivalence, tick-projection corollary.** The count form is
    equivalent to `tickProjection alphaU N = β · P` when `β = 4` and
    `24 = 6 · β`; the hypothesis `h24` records in the statement itself that
    the exponent numerator `6` is `m_rep / β_EW` with `m_rep = 24`. -/
theorem bridge_equivalence_tick {beta P alphaU N : ℝ}
    (hb : beta = 4) (h24 : mRep = 6 * beta)
    (hP : 0 < P) (ha : 0 < alphaU) (hN : Real.pi < N) :
    N = Real.pi * Real.exp (6 * Real.pi / (P * alphaU)) ↔
      tickProjection alphaU N = beta * P := by
  subst hb
  have h6 : mRep / (4 : ℝ) = 6 := by rw [h24]; norm_num
  have h := bridge_equivalence_beta (beta := (4 : ℝ)) (by norm_num) hP ha hN
  rwa [h6] at h

-- Axiom audit: these must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms six_eq_mRep_div_betaEW
#print axioms mRep_eq_six_mul_betaEW
#print axioms bridge_equivalence_beta
#print axioms bridge_equivalence
#print axioms bridge_equivalence_tick

end OPH.BridgeEquivalence
