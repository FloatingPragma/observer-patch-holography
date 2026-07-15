# Mathlib notes — EventAlgebra development (2026-07-15)

Running log of every Mathlib gap, friction point, and affordance hit while
building the `EventAlgebra` library (toolchain `leanprover/lean4:v4.29.1`,
Mathlib pinned per `lake-manifest.json`). Kept as raw material for the
paper's formalisation-experience section.

## Scoped instances (friction, easily missed)

1. **The partial order on `ℂ` is scoped.** `Matrix.PosSemidef` over `ℂ`
   and all `0 ≤ z` statements need `open scoped ComplexOrder`
   (`Mathlib/Analysis/Complex/Order.lean`, `scoped[ComplexOrder] attribute
   [instance] Complex.partialOrder`). Without it, goals mentioning the
   order fail to elaborate with unhelpful instance errors.
2. **The Loewner order on matrices is scoped** (`MatrixOrder`,
   `Mathlib/Analysis/Matrix/Order.lean`), as is the ordered-ring structure
   of `ℂ`/`RCLike` (`RCLike.toIsStrictOrderedRing`, scoped in
   `ComplexOrder`) that provides `mul_nonneg` in the ℂ order.
3. **The C*-structure on matrices is scoped**:
   `open scoped Matrix.Norms.L2Operator` activates the L2 operator norm,
   `NormedRing`, and `CStarRing` instances
   (`Mathlib/Analysis/CStarAlgebra/Matrix.lean`). There is no global norm
   on `Matrix`, so the abstract C*-statement of the Tsirelson bound and
   its matrix instantiation must live behind this scope. `Nontrivial
   (Matrix (Fin n) (Fin n) ℂ)` for `n ≠ 0` also had to be constructed by
   hand (entrywise, via `Matrix.one_apply_eq`); no instance fired.

## Genuine gaps and workarounds

4. **No `Tr(AB) ≥ 0` for a product of two PSD matrices.**
   `Matrix.PosSemidef.trace_nonneg` exists for a single PSD matrix, and
   sandwiches `B * A * Bᴴ` are covered
   (`PosSemidef.mul_mul_conjTranspose_same`), but the bilinear statement
   is absent. For Born weights this is harmless (`Tr(ρP) = Tr(PρP)` via
   `trace_mul_cycle` for idempotent `P`); for the general expectation
   functional (`expectation_nonneg`) it forced a real detour (next item).
5. **Gram factorisation `0 ≤ M ↔ ∃ B, M = Bᴴ B` is CFC-gated and did not
   synthesize.** The natural route,
   `CStarAlgebra.nonneg_iff_eq_star_mul_self` + `Matrix.PosSemidef.nonneg`
   (Loewner order), fails on
   `Matrix (Fin n) (Fin n) ℂ` with two unsynthesized classes:
   `NonUnitalContinuousFunctionalCalculus ℝ (Matrix (Fin n) (Fin n) ℂ)
   IsSelfAdjoint` and `NonnegSpectrumClass ℝ (Matrix (Fin n) (Fin n) ℂ)` —
   even with `ComplexOrder` and `MatrixOrder` open. (Mathlib's own
   `Analysis/Matrix/Order.lean` uses the lemma successfully inside its own
   section context; downstream, with the Pi topology on matrices, the
   instance chain did not resolve for us.) **Workaround:** the elementary
   spectral route — `Matrix.IsHermitian.spectral_theorem`
   (`M = U D U⋆` via `Unitary.conjStarAlgAut_apply`),
   `PosSemidef.eigenvalues_nonneg`, cycle the trace, and sum
   termwise-nonnegative diagonal products. Robust and about a dozen lines.
6. **No noncommutative `linear_combination`.** The CHSH squaring identity
   `S² = 4·1 − [a₀,a₁][b₀,b₁]` is pure noncommutative-ring algebra with
   side relations (involutivity, `aᵢbⱼ = bⱼaᵢ`); `noncomm_ring` cannot use
   hypotheses and `linear_combination` is commutative-only.
   **Workaround:** a hand-rolled confluent rewrite system as a `simp only`
   set — ∀-quantified association-compatible forms
   `b·(a·x) = a·(b·x)` and `a·(a·x) = x` derived from the hypotheses by
   reassociation — followed by `abel`. Terminates because each rewrite
   strictly decreases the number of letters out of normal order.
7. **`rw` on a projection-heavy goal picks the wrong occurrence.**
   Rewriting with the definition equation of `bornWeight` unfolds the
   *first* syntactic occurrence (often the normalising scalar
   `(bornWeight ρ P)⁻¹` rather than the intended trace argument).
   **Workaround:** isolate each trace identity as a standalone `have`
   with explicit arguments to `trace_mul_comm`/`trace_mul_cycle`.
8. **Unicode identifier surprise:** `𝒫` (script capital P, U+1D4AB) is
   not a legal identifier character in Lean 4.29; the partition variable
   had to be plain ASCII.

## Affordances that made things short

9. **`Matrix.PosSemidef` API is definition-agnostic.** The pinned Mathlib
   defines `PosSemidef` via finitely-supported vectors (`n →₀ R`); we
   never had to unfold it — `mul_mul_conjTranspose_same`,
   `conjTranspose_mul_mul_same`, `posSemidef_sum`, `PosSemidef.smul`,
   `diag_nonneg`, `trace_nonneg`, `trace_eq_zero_iff`,
   `posSemidef_conjTranspose_mul_self` covered every need.
   `PosSemidef.smul` accepts a complex scalar with `0 ≤ a` in the ℂ order
   (`PosSMulMono ℂ ℂ` resolves under `ComplexOrder`), so the Lüders
   normalisation is one lemma application.
10. **Trace faithfulness is in the box:**
    `Matrix.trace_mul_conjTranspose_self_eq_zero_iff` powers the
    uniqueness of the pinching, and
    `PosSemidef.trace_eq_zero_iff` + `PosSemidef.dotProduct_mulVec_zero_iff`
    power the certainty-support lemma (`σP = σ` for states certain of `P`).
11. **C*-norm machinery:** `CStarRing.norm_star_mul_self`,
    `IsSelfAdjoint.norm_mul_self` (`‖x·x‖ = ‖x‖²`), the `NormOneClass`
    instance for nontrivial C*-rings, and `norm_nsmul_le` reduce the
    analytic half of the Tsirelson bound to ~40 lines over bare
    `[NormedRing A] [StarRing A] [CStarRing A] [Nontrivial A]` — no
    completeness, no `CStarAlgebra` bundle needed.
12. **Mathlib prior art:** `Mathlib/Algebra/Star/CHSH.lean` already
    contains `IsCHSHTuple` and the *order-form* Tsirelson inequality
    `S ≤ √2³ • 1` in a star-ordered ℝ-algebra (Kim Morrison). Our
    `tsirelson_bound` is the complementary *norm-form* statement in a
    C*-ring; `tsirelson_bound_of_isCHSHTuple` consumes Mathlib's
    hypothesis bundle directly, so the two results are interoperable.

## Comparison notes vs. the Isabelle/HOL line (Echenim–Mhalla)

For the paper's related-work section; their JAR 2024 development and AFP
entries prove CHSH/Tsirelson matrix-first over a bespoke density-matrix and
L2-norm stack.

- **Shorter in Lean/Mathlib:** (a) the norm-form Tsirelson bound is proved
  *abstractly* for any unital C*-ring — the C*-identity, `norm_mul_le`,
  and `NormOneClass` are typeclass-generic, so the finite matrix case is a
  two-line instantiation rather than the theorem; (b) Born-weight
  positivity/upper bounds live directly in the partial order of `ℂ`
  (`StarOrderedRing`), which is strictly stronger than the `re`-projected
  statements and needed no bespoke real-part bookkeeping; (c) PSD
  closure lemmas (sandwich, sum, scalar) came from the general
  `Matrix.PosSemidef` API.
- **Longer in Lean/Mathlib:** (a) the PSD Gram factorisation that a
  concrete matrix stack gets from its spectral toolbox is CFC-gated in
  Mathlib and the instance chain failed for the Pi-topology matrix
  algebra, forcing the manual spectral-theorem detour of item 5;
  (b) scoped-instance discipline (`ComplexOrder`, `MatrixOrder`,
  `Matrix.Norms.L2Operator`) is a real usability tax that a bespoke stack
  does not pay.
- **Organisational difference (the paper's claim):** the present bundle
  separates algebra-only content from trace-dependent content lemma-by-lemma
  (machine-visible doc tags), and routes classicality through an explicit
  conditional-expectation projector onto a distinguished block-diagonal
  commutant (`ProjectivePartition` / `pinchingExpectation`), with
  existence/uniqueness/contractivity and the Lüders-compatibility law
  `𝔼 ∘ L_P = L_P ∘ 𝔼` for commutant `P`. That organisation, not the CHSH
  arithmetic, is the delta over the matrix-first prior art.
