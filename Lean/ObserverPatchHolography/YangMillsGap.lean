import Mathlib
import ObserverPatchHolography.YangMillsLemma72
import ObserverPatchHolography.YangMillsProp81

/-!
# Yang–Mills finite repair-gap — legacy commuting-projection branch

This module preserves a special conditional result from an earlier version of
B. Müller, *Explaining the Yang–Mills Mass Gap with Observer-Patch Repair
Dynamics*: a finite-gap assembly for a family of pairwise commuting
orthogonal projections.  It is **not** a formalisation of the current paper's
finite-stage theorem, whose collar conditional expectations need not commute
and whose gap argument instead uses a Dobrushin/approximate-tensorisation
hypothesis.

## Honest scope (read this first)

This module proves an **implication**, nothing more:

> IF the repair generator on a finite collar family has the special structural
> form assumed here — each collar acts by an orthogonal projection `E_C`, the
> collars **mutually commute**, their joint fixed space is represented by
> `P₀`, and there are **finitely many** active collar types each with a
> **strictly positive** relaxation rate `c_C` (from Lemma 7.2) — THEN the
> finite generator `L_r^rep = ∑_C c_C · (I − E_C)` dominates `c_* · (I − P₀)`
> with `c_* = min_C c_C > 0` (Proposition 8.1): a strictly positive
> finite-stage spectral gap.

Everything to the left of `THEN` is a **hypothesis** here.  In particular,
pairwise commutation is not derived from the present Yang–Mills construction
and is discharged **nowhere** in Lean.  This file proves only the legacy/special
implication `commuting structure ⇒ finite gap`.

The result must not be substituted for, or described as, the current paper's
noncommuting Dobrushin finite-stage theorem.  Nor can it be composed into the
current continuum receipt merely because both conclusions are called a finite
gap: doing so would require a separate source-level proof that the current
collar expectations commute and that this commuting branch satisfies the
continuum certificate's compatibility hypotheses.  No such bridge is present
here.

One thing is deliberately **NOT** proved and **NOT** claimed:

* `Δ_YM = Δ_rep`. This is Müller's **Assumption 9.2** — the continuum
  certificate (Schwinger-function convergence, reflection positivity,
  Osterwalder–Schrader reconstruction, non-triviality). It is the genuine open
  problem and is **untouched**. The deliverable here is `Δ_rep`, the
  finite-stage *representation* gap — **not** the physical Yang–Mills mass gap.

The two keystone lemmas this legacy branch *assembles* are proved in sibling modules and
**imported** here (no `sorry`s remain):

* `lemma_7_2` — scalar relaxation on a uniform hidden fiber (each `c_C > 0`);
  proved in `ObserverPatchHolography.YangMillsLemma72` (extracted verbatim from
  the retired single-file artifact `RepairGap.lean`, Part I — that file is now
  consolidated into these modules; see git history for the original);
* `prop_8_1`  — commuting-color finite-stage gap (`commuting ⇒ gap`);
  proved in `ObserverPatchHolography.YangMillsProp81`.

This file is the historical **setup + Theorem 7.3 / Lemma 7.4 assembly**; with the
keystone imports wired, `thm_7_3_finite_gap` carries zero `sorry`s and no
project-level axioms. The conditional continuum chain lives in the sibling
`ObserverPatchHolography.RepairGapChain`.

SCOPE: machine-checked only on the special pairwise-commuting branch: the
finite representation-gap implication `Δ_rep ≥ c_* > 0` (Lemma 7.2 / Lemma
7.4 / Prop 8.1 / Thm 7.3 assembly).  This is neither the current paper's
noncommuting finite-stage theorem nor a proof of any continuum-certificate
premise.
-/

namespace ObserverPatchHolography.YangMillsGap

/-! ## §7.2 keystone — imported statement (uniform hidden fiber) -/

/-- `E_F`, the expectation onto the constants of a uniform hidden fiber `F`,
    realised as the matrix `|F|⁻¹ · J` where `J` is the all-ones matrix. -/
noncomputable def EF (F : Type*) [Fintype F] : Matrix F F ℝ :=
  (Fintype.card F : ℝ)⁻¹ • Matrix.of (fun _ _ => (1 : ℝ))

/-- **Lemma 7.2 (imported keystone).** A positive-semidefinite relaxation `D` on
    a uniform hidden fiber `F` with `|F| ≥ 2`, commuting with the full symmetric
    group `S_F` (every permutation matrix) and with kernel exactly the constants,
    is a *strictly positive scalar* multiple of `I − E_F`.

    Proved in the sibling keystone module `YangMillsLemma72` via the commutant
    of the permutation representation (`commutant_perm_two_valued` + Schur on
    the standard rep); discharged here by direct import — this file is the
    assembly, not the keystone. (The two `EF` definitions are token-identical,
    hence definitionally equal.) -/
theorem lemma_7_2 {F : Type*} [Fintype F] [DecidableEq F]
    (hF : 2 ≤ Fintype.card F) (D : Matrix F F ℝ)
    (hPSD : D.PosSemidef)
    (hComm : ∀ σ : Equiv.Perm F, Commute (σ.permMatrix ℝ) D)
    (hKer : ∀ v : F → ℝ, D.mulVec v = 0 ↔ ∃ c : ℝ, v = fun _ => c) :
    ∃ cF : ℝ, 0 < cF ∧ D = cF • ((1 : Matrix F F ℝ) - EF F) :=
  ObserverPatchHolography.YangMillsLemma72.lemma_7_2 hF D hPSD hComm hKer

/-- Convenience projection of `lemma_7_2` down to its positivity component.

    **Read the statement, not the proof.** The conclusion `∃ cF, 0 < cF` mentions none of the
    hypotheses and is provable by `⟨1, by norm_num⟩` alone, so as a *statement* this carries no
    information: it is a weaker sibling of the re-export directly above, which keeps the real content
    `D = cF • (1 - EF F)`. The proof does invoke `lemma_7_2`, but a reader auditing statements would
    be misled by that. Nothing in the assembly consumes this lemma (an earlier docstring claimed it
    did, which was false); it is retained only as a named entry point for `#print axioms`. Prefer
    `lemma_7_2` — or the full-strength re-export above — for any real use. -/
theorem collar_rate_pos {F : Type*} [Fintype F] [DecidableEq F]
    (hF : 2 ≤ Fintype.card F) (D : Matrix F F ℝ)
    (hPSD : D.PosSemidef)
    (hComm : ∀ σ : Equiv.Perm F, Commute (σ.permMatrix ℝ) D)
    (hKer : ∀ v : F → ℝ, D.mulVec v = 0 ↔ ∃ c : ℝ, v = fun _ => c) :
    ∃ cF : ℝ, 0 < cF :=
  let ⟨cF, hpos, _⟩ := lemma_7_2 hF D hPSD hComm hKer
  ⟨cF, hpos⟩

/-! ## Operator setup: a real Hilbert space

(The paper's application is finite-dimensional, whence complete; only
completeness is consumed by the Loewner/positivity machinery, so we
hypothesize exactly that — the finite-dimensional case is an instance.) -/

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
  [CompleteSpace E]

/-- **The repair generator** `L_r^rep = ∑_{C ∈ s} c_C · (I − E_C)` over a finite
    active-collar index set `s`, with per-collar orthogonal projections `Ec` and
    relaxation rates `rate`. -/
noncomputable def repairGenerator {ι : Type*} (s : Finset ι)
    (Ec : ι → (E →L[ℝ] E)) (rate : ι → ℝ) : E →L[ℝ] E :=
  ∑ a ∈ s, rate a • ((1 : E →L[ℝ] E) - Ec a)

/-! ## §8.1 legacy gap engine — imported commuting-projection statement -/

/-- **Proposition 8.1 (imported).** For a finite family of **mutually commuting**
    orthogonal (star) projections `Ec` whose non-commutative product equals the
    projection `P₀` onto the joint fixed (constants) space, the constant-rate
    generator `∑ c_* · (I − E_C)` dominates `c_* · (I − P₀)`.

    Proved in the sibling gap module `YangMillsProp81` via `1 − ∏ Eₐ ≤ ∑ (1 − Eₐ)`
    for a commuting family of star projections; discharged here by direct
    import. -/
theorem prop_8_1 {ι : Type*} (s : Finset ι) (Ec : ι → (E →L[ℝ] E)) (P0 : E →L[ℝ] E)
    (hE : ∀ a ∈ s, IsStarProjection (Ec a))
    (hc : (↑s : Set ι).Pairwise (Function.onFun Commute Ec))
    (hprod : s.noncommProd Ec hc = P0)
    {cstar : ℝ} (hcpos : 0 < cstar) :
    cstar • ((1 : E →L[ℝ] E) - P0) ≤ ∑ a ∈ s, cstar • ((1 : E →L[ℝ] E) - Ec a) :=
  ObserverPatchHolography.YangMillsProp81.prop_8_1 s Ec P0 hE hc hprod hcpos

/-! ## §7.4 uniform floor -/

/-- **Lemma 7.4 (uniform floor).** Finitely many strictly positive relaxation
    rates over a nonempty active-collar set have a strictly positive minimum
    `c_*` that bounds them all from below. -/
theorem uniform_floor {ι : Type*} (s : Finset ι) (hne : s.Nonempty) (rate : ι → ℝ)
    (hpos : ∀ a ∈ s, 0 < rate a) :
    ∃ cstar : ℝ, 0 < cstar ∧ ∀ a ∈ s, cstar ≤ rate a := by
  obtain ⟨a₀, ha₀, hmin⟩ := s.exists_min_image rate hne
  exact ⟨rate a₀, hpos a₀ ha₀, hmin⟩

/-! ## §7.3 legacy assembly — the special commuting finite-gap deliverable -/

/-- **Legacy Theorem 7.3 / 7.4 (commuting finite representation gap).** Combine Lemma 7.2 (each
    collar rate `rate a > 0`), Proposition 8.1 (commuting colors ⇒ constant-rate
    gap), and the uniform floor (Lemma 7.4): the repair generator dominates
    `c_* · (I − P₀)` with `c_* > 0`.

    This is `Δ_rep ≥ c_* > 0` only on the special pairwise-commuting branch.
    It is not the current paper's noncommuting Dobrushin theorem and supplies
    no premise of the current continuum receipt without an additional
    compatibility proof.  It says nothing about `Δ_YM`; the continuum bridge
    is not stated, assumed, or proved here.

    Proof shape: `c_* · (I − P₀) ≤ ∑ c_* · (I − E_C) ≤ ∑ c_C · (I − E_C)`, where
    the first `≤` is Prop 8.1 and the second is rate-monotonicity — each summand
    grows because `c_C ≥ c_*` (uniform floor) and `I − E_C ⪰ 0` (star
    projection), so the difference is a sum of positive operators. -/
theorem thm_7_3_finite_gap {ι : Type*} (s : Finset ι) (hne : s.Nonempty)
    (Ec : ι → (E →L[ℝ] E)) (P0 : E →L[ℝ] E)
    (hE : ∀ a ∈ s, IsStarProjection (Ec a))
    (hc : (↑s : Set ι).Pairwise (Function.onFun Commute Ec))
    (hprod : s.noncommProd Ec hc = P0)
    (rate : ι → ℝ) (hrate : ∀ a ∈ s, 0 < rate a) :
    ∃ cstar : ℝ, 0 < cstar ∧
      cstar • ((1 : E →L[ℝ] E) - P0) ≤ repairGenerator s Ec rate := by
  -- Lemma 7.4: a strictly positive uniform floor c_* = min_C c_C.
  obtain ⟨cstar, hcstar_pos, hfloor⟩ := uniform_floor s hne rate hrate
  refine ⟨cstar, hcstar_pos, ?_⟩
  -- Step 1 (Prop 8.1): the constant-rate gap  c_* (I − P₀) ≤ ∑ c_* (I − E_C).
  have step1 :
      cstar • ((1 : E →L[ℝ] E) - P0) ≤ ∑ a ∈ s, cstar • ((1 : E →L[ℝ] E) - Ec a) :=
    prop_8_1 s Ec P0 hE hc hprod hcstar_pos
  -- Step 2: actual rates dominate the floor rate collar-by-collar.
  have step2 :
      (∑ a ∈ s, cstar • ((1 : E →L[ℝ] E) - Ec a))
        ≤ ∑ a ∈ s, rate a • ((1 : E →L[ℝ] E) - Ec a) := by
    have hdiff :
        (∑ a ∈ s, (rate a - cstar) • ((1 : E →L[ℝ] E) - Ec a))
          = (∑ a ∈ s, rate a • ((1 : E →L[ℝ] E) - Ec a))
              - ∑ a ∈ s, cstar • ((1 : E →L[ℝ] E) - Ec a) := by
      simp only [sub_smul, Finset.sum_sub_distrib]
    refine (ContinuousLinearMap.le_def _ _).mpr ?_
    rw [← hdiff]
    refine ContinuousLinearMap.isPositive_sum s (fun a ha => ?_)
    exact ContinuousLinearMap.IsPositive.smul_of_nonneg
      (ContinuousLinearMap.IsPositive.of_isStarProjection (hE a ha).one_sub)
      (sub_nonneg.mpr (hfloor a ha))
  -- Assemble: c_* (I − P₀) ≤ ∑ c_* (I − E_C) ≤ ∑ c_C (I − E_C) = L_r^rep.
  show cstar • ((1 : E →L[ℝ] E) - P0) ≤ ∑ a ∈ s, rate a • ((1 : E →L[ℝ] E) - Ec a)
  exact le_trans step1 step2

/-! ## Axiom self-audit (build-log visible)

Expected report for every theorem below: exactly
`[propext, Classical.choice, Quot.sound]` — the three standard Mathlib axioms,
no `sorryAx`, no project-level axiom. -/

#print axioms lemma_7_2
#print axioms collar_rate_pos
#print axioms prop_8_1
#print axioms uniform_floor
#print axioms thm_7_3_finite_gap

end ObserverPatchHolography.YangMillsGap
