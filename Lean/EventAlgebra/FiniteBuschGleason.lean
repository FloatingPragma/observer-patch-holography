import EventAlgebra.Basic
import EventAlgebra.StateExpectation
import EventAlgebra.FiniteEffectClosureBoundary

/-!
# The finite Busch--Gleason representation theorem

This module proves the effect-algebra representation theorem in finite
dimension `d`: every real-valued assignment `v` on the effects of
`Matrix (Fin d) (Fin d) ℂ` that

* takes values in `[0, 1]`,
* sends the sure effect `1` to `1`, and
* is **additive across contexts** (`v (E + F) = v E + v F` whenever `E`,
  `F`, and `E + F` are all effects)

is the Born functional `E ↦ Re (Tr (ρ E))` of a unique density matrix `ρ`.
Together with the converse (`isEffectValuation_born`), Born weights are
exactly the noncontextually additive effect assignments
(`isEffectValuation_iff_born`).

## Proof route

1. Additivity forces `v 0 = 0` and monotonicity along effect splittings.
2. **Rational homogeneity**: splitting an effect `E` into `n` equal effect
   summands `n⁻¹ • E` gives `v ((p / n) • E) = (p / n) * v E`.
3. **Real homogeneity**: for `c ∈ [0, 1]` the floor fractions `⌊c n⌋ / n`
   sandwich `v (c • E)` between `(c - 1/n) * v E` and `(c + 1/n) * v E`;
   a finite Archimedean squeeze (no continuity axiom) yields
   `v (c • E) = c * v E`.
4. **Extension to Hermitian observables**: every positive semidefinite `A`
   satisfies `A ≤ (Re (Tr A) + 1) • 1`, so `posValue v A := t * v (t⁻¹ • A)`
   is well defined for any admissible bound `t`, additive, and positively
   homogeneous; every Hermitian `H` admits a shift with `H + s • 1`
   positive semidefinite, exhibiting `H = t • E - s • 1` with `E` an
   effect, and `hermValue v H := posValue v (H + s • 1) - s` is a
   real-linear extension of `v`.
5. **Finite Riesz representation**: expanding Hermitian matrices in the
   basis of symmetrised and antisymmetrised matrix units determines a
   matrix `gleasonState v` with `Tr (gleasonState v * H) = hermValue v H`
   for every Hermitian `H`.
6. **State properties**: positivity of `v` on rank-one projectors makes
   `gleasonState v` positive semidefinite, and `v 1 = 1` gives unit trace;
   testing on the Hermitian basis gives uniqueness.

## Negative control

The cube response of `EventAlgebra.FiniteEffectClosureBoundary` extends to
qubit effects as `cubeResponse`: it takes values in `[0, 1]` on effects,
sends `1` to `1`, and satisfies per-context normalisation
`cubeResponse E + cubeResponse (1 - E) = 1`, yet it is not of Born form.
The failure point is exhibited exactly: on the unsharp three-effect
decomposition `1 = halfUp + halfUp + downProj` in dimension two the cube
response sums to `5/8` instead of `1`, so cross-context additivity fails.
This closes the loop with the dense-test/affinity boundary module: range,
normalisation, and continuity do not force Born form; additivity does.

## Tagging convention

As in `EventAlgebra.Basic`: **algebra-only** lemmas consume no valuation
and no trace functional; **valuation-dependent** lemmas consume the
hypotheses `IsEffectValuation`; **trace-dependent** lemmas pass through
the trace pairing.
-/

namespace EventAlgebra

open Matrix
open scoped ComplexOrder

variable {d : ℕ}

/-! ## Effects -/

/-- An **effect** of the finite algebra `Matrix (Fin d) (Fin d) ℂ`: a
matrix `E` with `0 ≤ E ≤ 1` in the positive-semidefinite order. -/
def IsEffect (E : Matrix (Fin d) (Fin d) ℂ) : Prop :=
  E.PosSemidef ∧ (1 - E).PosSemidef

/-- **Algebra-only.** Every projection event is an effect. -/
theorem IsEvent.isEffect {P : Matrix (Fin d) (Fin d) ℂ} (hP : IsEvent P) :
    IsEffect P :=
  ⟨hP.posSemidef, hP.compl.posSemidef⟩

/-- **Algebra-only.** The zero matrix is an effect. -/
theorem isEffect_zero : IsEffect (0 : Matrix (Fin d) (Fin d) ℂ) :=
  isEvent_zero.isEffect

/-- **Algebra-only.** The identity matrix is an effect. -/
theorem isEffect_one : IsEffect (1 : Matrix (Fin d) (Fin d) ℂ) :=
  isEvent_one.isEffect

/-- **Algebra-only.** Effects are Hermitian. -/
theorem IsEffect.isHermitian {E : Matrix (Fin d) (Fin d) ℂ}
    (hE : IsEffect E) : E.IsHermitian :=
  hE.1.isHermitian

/-- **Algebra-only.** The complement `1 - E` of an effect is an effect. -/
theorem IsEffect.compl {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E) :
    IsEffect (1 - E) :=
  ⟨hE.2, by simpa using hE.1⟩

/-- **Algebra-only.** A real scalar acts on a complex matrix as its
coercion: the bridge between the real-linear statements of this module and
the complex-scalar Mathlib API. -/
theorem real_smul_eq_coe_smul (c : ℝ) (E : Matrix (Fin d) (Fin d) ℂ) :
    c • E = ((c : ℂ)) • E := by
  ext i j
  simp [Matrix.smul_apply, Complex.real_smul]

/-- **Algebra-only.** Nonnegative real scalings preserve positive
semidefiniteness. -/
theorem posSemidef_real_smul {E : Matrix (Fin d) (Fin d) ℂ}
    (hE : E.PosSemidef) {c : ℝ} (hc : 0 ≤ c) : (c • E).PosSemidef := by
  rw [real_smul_eq_coe_smul]
  exact hE.smul (Complex.zero_le_real.mpr hc)

/-- **Algebra-only.** Scaling an effect by `c ∈ [0, 1]` yields an effect:
`1 - c • E = c • (1 - E) + (1 - c) • 1` is a sum of positive matrices. -/
theorem IsEffect.smul {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E)
    {c : ℝ} (h0 : 0 ≤ c) (h1 : c ≤ 1) : IsEffect (c • E) := by
  refine ⟨posSemidef_real_smul hE.1 h0, ?_⟩
  have key : (1 : Matrix (Fin d) (Fin d) ℂ) - c • E =
      c • (1 - E) + (1 - c) • 1 := by module
  rw [key]
  exact (posSemidef_real_smul hE.2 h0).add
    (posSemidef_real_smul Matrix.PosSemidef.one (by linarith))

/-- **Algebra-only.** The convex splitting of an effect into complementary
scalings: `E = c • E + (1 - c) • E`. -/
theorem effect_convex_split (E : Matrix (Fin d) (Fin d) ℂ) (c : ℝ) :
    E = c • E + (1 - c) • E := by module

/-! ## Additive effect valuations -/

/-- The hypotheses of the finite Busch--Gleason theorem: a real assignment
on matrices that maps effects into `[0, 1]`, is normalised on the sure
effect, and is additive whenever the sum of two effects is an effect.
The assignment is unconstrained away from the effects. -/
structure IsEffectValuation (v : Matrix (Fin d) (Fin d) ℂ → ℝ) : Prop where
  /-- Effects receive nonnegative values. -/
  nonneg : ∀ {E : Matrix (Fin d) (Fin d) ℂ}, IsEffect E → 0 ≤ v E
  /-- Effects receive values at most one. -/
  le_one : ∀ {E : Matrix (Fin d) (Fin d) ℂ}, IsEffect E → v E ≤ 1
  /-- The sure effect receives value one. -/
  map_one : v 1 = 1
  /-- Additivity across contexts: whenever `E`, `F`, and `E + F` are all
  effects, the value of the sum is the sum of the values. -/
  additive : ∀ {E F : Matrix (Fin d) (Fin d) ℂ}, IsEffect E → IsEffect F →
    IsEffect (E + F) → v (E + F) = v E + v F

variable {v : Matrix (Fin d) (Fin d) ℂ → ℝ}

/-- **Valuation-dependent.** Additivity at `0 + 0` forces `v 0 = 0`. -/
theorem IsEffectValuation.map_zero (hv : IsEffectValuation v) : v 0 = 0 := by
  have h := hv.additive isEffect_zero isEffect_zero (by simpa using isEffect_zero)
  simp only [add_zero] at h
  linarith

/-- **Valuation-dependent.** The complement rule `v (1 - E) = 1 - v E`. -/
theorem IsEffectValuation.map_compl (hv : IsEffectValuation v)
    {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E) :
    v (1 - E) = 1 - v E := by
  have h := hv.additive hE hE.compl (by simpa using isEffect_one)
  rw [add_sub_cancel, hv.map_one] at h
  linarith

/-- **Valuation-dependent.** Rational homogeneity, numerator induction: the
`j`-th partial sum of `n` equal effect summands has value `j` times the
value of one summand. -/
theorem IsEffectValuation.natDiv_smul (hv : IsEffectValuation v)
    {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E) {n : ℕ} (hn : 0 < n) :
    ∀ j : ℕ, j ≤ n →
      v (((j : ℝ) / n) • E) = j * v (((1 : ℝ) / n) • E) := by
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  intro j
  induction j with
  | zero =>
    intro _
    have h0 : ((((0 : ℕ) : ℝ)) / n) • E = 0 := by push_cast; module
    rw [h0, hv.map_zero]
    simp
  | succ k ih =>
    intro hk
    have hk' : k ≤ n := (Nat.le_succ k).trans hk
    have hkR : (k : ℝ) + 1 ≤ n := by exact_mod_cast hk
    have hsplit : (((k + 1 : ℕ) : ℝ) / n) • E =
        ((k : ℝ) / n) • E + ((1 : ℝ) / n) • E := by
      push_cast
      module
    have hEk : IsEffect (((k : ℝ) / n) • E) :=
      hE.smul (by positivity) (by rw [div_le_one hnR]; linarith)
    have hE1 : IsEffect (((1 : ℝ) / n) • E) :=
      hE.smul (by positivity) (by rw [div_le_one hnR]; exact_mod_cast hn)
    have hEsum : IsEffect ((((k + 1 : ℕ) : ℝ) / n) • E) :=
      hE.smul (by positivity) (by rw [div_le_one hnR]; exact_mod_cast hk)
    have h := hv.additive hEk hE1 (hsplit ▸ hEsum)
    rw [hsplit, h, ih hk']
    push_cast
    ring

/-- **Valuation-dependent.** Rational homogeneity: for `p ≤ n`,
`v ((p / n) • E) = (p / n) * v E`. -/
theorem IsEffectValuation.ratDiv_smul (hv : IsEffectValuation v)
    {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E) {n p : ℕ}
    (hn : 0 < n) (hp : p ≤ n) :
    v (((p : ℝ) / n) • E) = ((p : ℝ) / n) * v E := by
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  have hbase : v E = n * v (((1 : ℝ) / n) • E) := by
    have h := hv.natDiv_smul hE hn n le_rfl
    have hone : (((n : ℝ)) / n) • E = E := by
      rw [div_self hnR.ne']
      module
    rwa [hone] at h
  rw [hv.natDiv_smul hE hn p hp, hbase]
  field_simp

/-- **Valuation-dependent.** The floor-fraction lower bound behind the
Archimedean squeeze: for `c ∈ [0, 1]` and `0 < n`,
`c * v E - 1/n ≤ v (c • E)`. -/
theorem IsEffectValuation.smul_ge (hv : IsEffectValuation v)
    {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E) {c : ℝ}
    (h0 : 0 ≤ c) (h1 : c ≤ 1) {n : ℕ} (hn : 0 < n) :
    c * v E - ((n : ℝ))⁻¹ ≤ v (c • E) := by
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn
  set p : ℕ := ⌊c * n⌋₊ with hp
  have hpn : p ≤ n := by
    have : c * n ≤ n := by nlinarith
    calc p ≤ ⌊(n : ℝ)⌋₊ := Nat.floor_le_floor this
    _ = n := Nat.floor_natCast n
  have hple : (p : ℝ) / n ≤ c := by
    rw [div_le_iff₀ hnR]
    exact Nat.floor_le (by positivity)
  have hplt : c - ((n : ℝ))⁻¹ < (p : ℝ) / n := by
    rw [lt_div_iff₀ hnR]
    have hcancel : ((n : ℝ))⁻¹ * n = 1 := inv_mul_cancel₀ hnR.ne'
    have hfl : c * n < (p : ℝ) + 1 := by
      exact_mod_cast Nat.lt_floor_add_one (c * n)
    nlinarith
  have hsplit : c • E = ((p : ℝ) / n) • E + (c - (p : ℝ) / n) • E := by module
  have hEp : IsEffect (((p : ℝ) / n) • E) :=
    hE.smul (by positivity) ((div_le_one hnR).mpr (by exact_mod_cast hpn))
  have hpn0 : (0 : ℝ) ≤ (p : ℝ) / n := by positivity
  have hErest : IsEffect ((c - (p : ℝ) / n) • E) :=
    hE.smul (by linarith) (by linarith)
  have hEc : IsEffect (c • E) := hE.smul h0 h1
  have hadd := hv.additive hEp hErest (hsplit ▸ hEc)
  have hcv : v (c • E) =
      v (((p : ℝ) / n) • E) + v ((c - (p : ℝ) / n) • E) := by
    rw [hsplit]
    exact hadd
  have hval : v (((p : ℝ) / n) • E) = ((p : ℝ) / n) * v E :=
    hv.ratDiv_smul hE hn hpn
  have hnonneg : 0 ≤ v ((c - (p : ℝ) / n) • E) := hv.nonneg hErest
  have hvE0 : 0 ≤ v E := hv.nonneg hE
  have hvE1 : v E ≤ 1 := hv.le_one hE
  have hinv : (0 : ℝ) < (n : ℝ)⁻¹ := inv_pos.mpr hnR
  nlinarith [mul_le_mul_of_nonneg_right hplt.le hvE0,
    mul_le_mul_of_nonneg_left hvE1 hinv.le]

/-- **Valuation-dependent.** Real homogeneity on effects, by the
Archimedean sandwich between floor fractions: for `c ∈ [0, 1]`,
`v (c • E) = c * v E`. No continuity hypothesis enters; only additivity,
positivity, and the archimedean property of `ℝ`. -/
theorem IsEffectValuation.smul_eq (hv : IsEffectValuation v)
    {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E) {c : ℝ}
    (h0 : 0 ≤ c) (h1 : c ≤ 1) :
    v (c • E) = c * v E := by
  have key : ∀ n : ℕ, 0 < n → |v (c • E) - c * v E| ≤ ((n : ℝ))⁻¹ := by
    intro n hn
    have hlow := hv.smul_ge hE h0 h1 hn
    -- upper bound via the complement splitting `E = c • E + (1 - c) • E`
    have hEc : IsEffect (c • E) := hE.smul h0 h1
    have hErest : IsEffect ((1 - c) • E) := hE.smul (by linarith) (by linarith)
    have hadd := hv.additive hEc hErest ((effect_convex_split E c) ▸ hE)
    rw [← effect_convex_split E c] at hadd
    have hlow' := hv.smul_ge hE (by linarith : (0:ℝ) ≤ 1 - c)
      (by linarith : (1:ℝ) - c ≤ 1) hn
    have hhigh : v (c • E) ≤ c * v E + (n : ℝ)⁻¹ := by
      have : v (c • E) = v E - v ((1 - c) • E) := by linarith
      rw [this]
      have hvE1 : v E ≤ 1 := hv.le_one hE
      nlinarith
    rw [abs_le]
    constructor <;> linarith
  by_contra hne
  have habs : 0 < |v (c • E) - c * v E| := abs_pos.mpr (sub_ne_zero.mpr hne)
  obtain ⟨n, hn⟩ := exists_nat_one_div_lt habs
  have hkey := key (n + 1) n.succ_pos
  rw [one_div] at hn
  push_cast at hkey hn
  linarith

/-! ## Spectral shift bounds

Every Hermitian matrix sits between real multiples of the identity; a
positive semidefinite matrix is dominated by one more than its real
trace times the identity. These bounds drive the decomposition
`H = t • E - s • 1` of step (4). -/

/-- **Algebra-only.** Real smul composes by multiplication of the
scalars. -/
theorem real_smul_smul (a b : ℝ) (M : Matrix (Fin d) (Fin d) ℂ) :
    a • b • M = (a * b) • M := by
  ext i j
  simp [Matrix.smul_apply, Complex.real_smul]
  ring

/-- **Algebra-only.** Matrix multiplication commutes with real
scalings. -/
theorem mul_real_smul (M : Matrix (Fin d) (Fin d) ℂ) (c : ℝ)
    (N : Matrix (Fin d) (Fin d) ℂ) : M * (c • N) = c • (M * N) := by
  rw [real_smul_eq_coe_smul, real_smul_eq_coe_smul, Matrix.mul_smul]

/-- **Algebra-only.** Real scalings commute with matrix
multiplication. -/
theorem real_smul_mul (c : ℝ) (M N : Matrix (Fin d) (Fin d) ℂ) :
    (c • M) * N = c • (M * N) := by
  rw [real_smul_eq_coe_smul, real_smul_eq_coe_smul, Matrix.smul_mul]

/-- **Trace-dependent.** The trace of a real scaling. -/
theorem trace_real_smul (c : ℝ) (A : Matrix (Fin d) (Fin d) ℂ) :
    (c • A).trace = (c : ℂ) * A.trace := by
  rw [real_smul_eq_coe_smul, trace_smul, smul_eq_mul]

/-- **Algebra-only.** Conjugated form of the spectral theorem:
`H = V D V⋆` with `V` the eigenvector unitary and `D` the real diagonal
of eigenvalues. -/
theorem spectral_conj {H : Matrix (Fin d) (Fin d) ℂ}
    (hH : H.IsHermitian) :
    H = (hH.eigenvectorUnitary : Matrix (Fin d) (Fin d) ℂ) *
      diagonal (RCLike.ofReal ∘ hH.eigenvalues) *
      star (hH.eigenvectorUnitary : Matrix (Fin d) (Fin d) ℂ) := by
  conv_lhs => rw [hH.spectral_theorem, Unitary.conjStarAlgAut_apply]

/-- **Algebra-only.** Upper spectral shift: if every eigenvalue of the
Hermitian `H` is at most `t`, then `t • 1 - H` is positive
semidefinite. -/
theorem smul_one_sub_posSemidef {H : Matrix (Fin d) (Fin d) ℂ}
    (hH : H.IsHermitian) {t : ℝ} (ht : ∀ i, hH.eigenvalues i ≤ t) :
    (t • (1 : Matrix (Fin d) (Fin d) ℂ) - H).PosSemidef := by
  classical
  set V : Matrix (Fin d) (Fin d) ℂ :=
    (hH.eigenvectorUnitary : Matrix (Fin d) (Fin d) ℂ) with hV
  have hVstar : V * star V = 1 := Unitary.coe_mul_star_self _
  have hspec : H = V * diagonal (RCLike.ofReal ∘ hH.eigenvalues) * star V := by
    rw [hV]
    exact spectral_conj hH
  have hone : V * (((t : ℝ) : ℂ) • 1) * star V =
      t • (1 : Matrix (Fin d) (Fin d) ℂ) := by
    rw [Matrix.mul_smul, Matrix.smul_mul, mul_one, hVstar,
      ← real_smul_eq_coe_smul]
  have hdiag : ((t : ℂ)) • (1 : Matrix (Fin d) (Fin d) ℂ)
      - diagonal (RCLike.ofReal ∘ hH.eigenvalues)
      = diagonal (fun i => ((t - hH.eigenvalues i : ℝ) : ℂ)) := by
    ext i j
    by_cases h : i = j
    · subst h
      simp only [Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply_eq,
        diagonal_apply_eq, Function.comp_apply, smul_eq_mul, mul_one]
      push_cast
      rfl
    · simp only [Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply_ne h,
        diagonal_apply_ne _ h, smul_zero, sub_zero]
  have hkey : t • (1 : Matrix (Fin d) (Fin d) ℂ) - H =
      V * (diagonal (fun i => ((t - hH.eigenvalues i : ℝ) : ℂ))) * star V := by
    rw [← hdiag, mul_sub, sub_mul, hone, ← hspec]
  rw [hkey, Matrix.star_eq_conjTranspose]
  refine PosSemidef.mul_mul_conjTranspose_same ?_ V
  exact posSemidef_diagonal_iff.mpr fun i =>
    Complex.zero_le_real.mpr (sub_nonneg.mpr (ht i))

/-- **Algebra-only.** Lower spectral shift: if every eigenvalue of the
Hermitian `H` is at least `s`, then `H - s • 1` is positive
semidefinite. -/
theorem sub_smul_one_posSemidef {H : Matrix (Fin d) (Fin d) ℂ}
    (hH : H.IsHermitian) {s : ℝ} (hs : ∀ i, s ≤ hH.eigenvalues i) :
    (H - s • (1 : Matrix (Fin d) (Fin d) ℂ)).PosSemidef := by
  classical
  set V : Matrix (Fin d) (Fin d) ℂ :=
    (hH.eigenvectorUnitary : Matrix (Fin d) (Fin d) ℂ) with hV
  have hVstar : V * star V = 1 := Unitary.coe_mul_star_self _
  have hspec : H = V * diagonal (RCLike.ofReal ∘ hH.eigenvalues) * star V := by
    rw [hV]
    exact spectral_conj hH
  have hone : V * (((s : ℝ) : ℂ) • 1) * star V =
      s • (1 : Matrix (Fin d) (Fin d) ℂ) := by
    rw [Matrix.mul_smul, Matrix.smul_mul, mul_one, hVstar,
      ← real_smul_eq_coe_smul]
  have hdiag : diagonal (RCLike.ofReal ∘ hH.eigenvalues)
      - ((s : ℂ)) • (1 : Matrix (Fin d) (Fin d) ℂ)
      = diagonal (fun i => ((hH.eigenvalues i - s : ℝ) : ℂ)) := by
    ext i j
    by_cases h : i = j
    · subst h
      simp only [Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply_eq,
        diagonal_apply_eq, Function.comp_apply, smul_eq_mul, mul_one]
      push_cast
      rfl
    · simp only [Matrix.sub_apply, Matrix.smul_apply, Matrix.one_apply_ne h,
        diagonal_apply_ne _ h, smul_zero, sub_zero]
  have hkey : H - s • (1 : Matrix (Fin d) (Fin d) ℂ) =
      V * (diagonal (fun i => ((hH.eigenvalues i - s : ℝ) : ℂ))) * star V := by
    rw [← hdiag, mul_sub, sub_mul, hone, ← hspec]
  rw [hkey, Matrix.star_eq_conjTranspose]
  refine PosSemidef.mul_mul_conjTranspose_same ?_ V
  exact posSemidef_diagonal_iff.mpr fun i =>
    Complex.zero_le_real.mpr (sub_nonneg.mpr (hs i))

/-- **Algebra-only.** Every Hermitian matrix admits a nonnegative shift
into the positive cone: `H + s • 1` is positive semidefinite for
`s = ∑ |eigenvalues| + 1`. This is the existence half of the effect
decomposition `H = t • E - s • 1`. -/
theorem exists_add_smul_one_posSemidef
    {H : Matrix (Fin d) (Fin d) ℂ} (hH : H.IsHermitian) :
    ∃ s : ℝ, (H + s • 1).PosSemidef := by
  refine ⟨(∑ i, |hH.eigenvalues i|) + 1, ?_⟩
  have hbound : ∀ i, -((∑ i, |hH.eigenvalues i|) + 1) ≤ hH.eigenvalues i := by
    intro i
    have h1 : |hH.eigenvalues i| ≤ ∑ j, |hH.eigenvalues j| :=
      Finset.single_le_sum (f := fun j => |hH.eigenvalues j|)
        (fun j _ => abs_nonneg _) (Finset.mem_univ i)
    have h2 : -|hH.eigenvalues i| ≤ hH.eigenvalues i := neg_abs_le _
    linarith
  have h := sub_smul_one_posSemidef hH hbound
  have key : H - (-((∑ i, |hH.eigenvalues i|) + 1)) • (1 : Matrix (Fin d) (Fin d) ℂ) =
      H + ((∑ i, |hH.eigenvalues i|) + 1) • 1 := by module
  rwa [key] at h

/-! ## The extension to the positive cone -/

/-- The canonical dominating scalar of a positive matrix: one more than
its real trace. Every eigenvalue is at most the trace, so
`A ≤ posBound A • 1` strictly. -/
noncomputable def posBound (A : Matrix (Fin d) (Fin d) ℂ) : ℝ :=
  A.trace.re + 1

/-- **Algebra-only.** The canonical bound is positive on the positive
cone. -/
theorem posBound_pos {A : Matrix (Fin d) (Fin d) ℂ} (hA : A.PosSemidef) :
    0 < posBound A := by
  have h := hA.trace_nonneg
  have hre : 0 ≤ A.trace.re := by
    simpa using (Complex.le_def.mp h).1
  unfold posBound
  linarith

/-- **Algebra-only.** The canonical bound dominates:
`posBound A • 1 - A` is positive semidefinite. -/
theorem posBound_smul_one_sub {A : Matrix (Fin d) (Fin d) ℂ}
    (hA : A.PosSemidef) :
    (posBound A • (1 : Matrix (Fin d) (Fin d) ℂ) - A).PosSemidef := by
  apply smul_one_sub_posSemidef hA.1
  intro i
  have htr : A.trace.re = ∑ j, hA.1.eigenvalues j := by
    rw [hA.1.trace_eq_sum_eigenvalues, Complex.re_sum]
    exact Finset.sum_congr rfl fun j _ => by
      simp [Complex.coe_algebraMap]
  have hle : hA.1.eigenvalues i ≤ ∑ j, hA.1.eigenvalues j :=
    Finset.single_le_sum (fun j _ => hA.eigenvalues_nonneg j)
      (Finset.mem_univ i)
  unfold posBound
  linarith

/-- **Algebra-only.** Normalising a positive matrix by any dominating
positive scalar produces an effect. -/
theorem isEffect_inv_smul_of_bound {A : Matrix (Fin d) (Fin d) ℂ}
    (hA : A.PosSemidef) {t : ℝ} (ht : 0 < t)
    (hb : (t • (1 : Matrix (Fin d) (Fin d) ℂ) - A).PosSemidef) :
    IsEffect (t⁻¹ • A) := by
  refine ⟨posSemidef_real_smul hA (inv_nonneg.mpr ht.le), ?_⟩
  have key : (1 : Matrix (Fin d) (Fin d) ℂ) - t⁻¹ • A =
      t⁻¹ • (t • (1 : Matrix (Fin d) (Fin d) ℂ) - A) := by
    rw [show t⁻¹ • (t • (1 : Matrix (Fin d) (Fin d) ℂ) - A)
        = (t⁻¹ * t) • (1 : Matrix (Fin d) (Fin d) ℂ) - t⁻¹ • A from by module,
      inv_mul_cancel₀ ht.ne']
    module
  rw [key]
  exact posSemidef_real_smul hb (inv_nonneg.mpr ht.le)

/-- **Algebra-only.** The canonical normalisation of a positive matrix is
an effect. -/
theorem isEffect_posBound_inv_smul {A : Matrix (Fin d) (Fin d) ℂ}
    (hA : A.PosSemidef) : IsEffect ((posBound A)⁻¹ • A) :=
  isEffect_inv_smul_of_bound hA (posBound_pos hA) (posBound_smul_one_sub hA)

/-- The extension of an effect valuation to the positive cone: rescale by
the canonical bound, evaluate on the resulting effect, and scale back. -/
noncomputable def posValue (v : Matrix (Fin d) (Fin d) ℂ → ℝ)
    (A : Matrix (Fin d) (Fin d) ℂ) : ℝ :=
  posBound A * v ((posBound A)⁻¹ • A)

/-- **Valuation-dependent.** Well-definedness of the positive-cone
extension: any dominating positive scalar computes the same value. The
proof is the real homogeneity of step (3) applied to the ratio of two
bounds. -/
theorem posValue_eq_of_bound (hv : IsEffectValuation v)
    {A : Matrix (Fin d) (Fin d) ℂ} (hA : A.PosSemidef) {t : ℝ} (ht : 0 < t)
    (hb : (t • (1 : Matrix (Fin d) (Fin d) ℂ) - A).PosSemidef) :
    posValue v A = t * v (t⁻¹ • A) := by
  have main : ∀ s u : ℝ, 0 < s → 0 < u → s ≤ u →
      ((s • (1 : Matrix (Fin d) (Fin d) ℂ) - A)).PosSemidef →
      u * v (u⁻¹ • A) = s * v (s⁻¹ • A) := by
    intro s u hs hu hsu hbs
    have hEs : IsEffect (s⁻¹ • A) := isEffect_inv_smul_of_bound hA hs hbs
    have hcalc : u⁻¹ • A = (s / u) • (s⁻¹ • A) := by
      rw [real_smul_smul]
      congr 1
      field_simp
    have h01 : 0 ≤ s / u := by positivity
    have h11 : s / u ≤ 1 := (div_le_one hu).mpr hsu
    rw [hcalc, hv.smul_eq hEs h01 h11]
    field_simp
  rcases le_total (posBound A) t with h | h
  · exact (main (posBound A) t (posBound_pos hA) ht h
      (posBound_smul_one_sub hA)).symm
  · exact main t (posBound A) ht (posBound_pos hA) h hb

/-- **Valuation-dependent.** The positive-cone extension restricts to the
valuation on effects. -/
theorem posValue_of_isEffect (hv : IsEffectValuation v)
    {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E) :
    posValue v E = v E := by
  have hb : ((1 : ℝ) • (1 : Matrix (Fin d) (Fin d) ℂ) - E).PosSemidef := by
    have key : (1 : ℝ) • (1 : Matrix (Fin d) (Fin d) ℂ) - E = 1 - E := by
      module
    rw [key]
    exact hE.2
  rw [posValue_eq_of_bound hv hE.1 one_pos hb, inv_one]
  have h1 : (1 : ℝ) • E = E := by module
  rw [h1, one_mul]

/-- **Valuation-dependent.** Additivity of the positive-cone extension:
both summands and the sum are computed with the common bound
`posBound A + posBound B`, where additivity of the valuation applies. -/
theorem posValue_add (hv : IsEffectValuation v)
    {A B : Matrix (Fin d) (Fin d) ℂ} (hA : A.PosSemidef)
    (hB : B.PosSemidef) :
    posValue v (A + B) = posValue v A + posValue v B := by
  have hs := posBound_pos hA
  have ht := posBound_pos hB
  have hst : 0 < posBound A + posBound B := by linarith
  have hboundA : ((posBound A + posBound B) • (1 : Matrix (Fin d) (Fin d) ℂ)
      - A).PosSemidef := by
    have key : (posBound A + posBound B) • (1 : Matrix (Fin d) (Fin d) ℂ) - A
        = (posBound A • 1 - A) + posBound B • 1 := by module
    rw [key]
    exact (posBound_smul_one_sub hA).add
      (posSemidef_real_smul Matrix.PosSemidef.one ht.le)
  have hboundB : ((posBound A + posBound B) • (1 : Matrix (Fin d) (Fin d) ℂ)
      - B).PosSemidef := by
    have key : (posBound A + posBound B) • (1 : Matrix (Fin d) (Fin d) ℂ) - B
        = posBound A • 1 + (posBound B • 1 - B) := by module
    rw [key]
    exact (posSemidef_real_smul Matrix.PosSemidef.one hs.le).add
      (posBound_smul_one_sub hB)
  have hboundAB : ((posBound A + posBound B) • (1 : Matrix (Fin d) (Fin d) ℂ)
      - (A + B)).PosSemidef := by
    have key : (posBound A + posBound B) • (1 : Matrix (Fin d) (Fin d) ℂ)
        - (A + B) = (posBound A • 1 - A) + (posBound B • 1 - B) := by module
    rw [key]
    exact (posBound_smul_one_sub hA).add (posBound_smul_one_sub hB)
  have hEA : IsEffect ((posBound A + posBound B)⁻¹ • A) :=
    isEffect_inv_smul_of_bound hA hst hboundA
  have hEB : IsEffect ((posBound A + posBound B)⁻¹ • B) :=
    isEffect_inv_smul_of_bound hB hst hboundB
  have hEAB : IsEffect ((posBound A + posBound B)⁻¹ • (A + B)) :=
    isEffect_inv_smul_of_bound (hA.add hB) hst hboundAB
  have hsplit : (posBound A + posBound B)⁻¹ • (A + B) =
      (posBound A + posBound B)⁻¹ • A + (posBound A + posBound B)⁻¹ • B := by
    module
  rw [posValue_eq_of_bound hv (hA.add hB) hst hboundAB,
    posValue_eq_of_bound hv hA hst hboundA,
    posValue_eq_of_bound hv hB hst hboundB, hsplit,
    hv.additive hEA hEB (hsplit ▸ hEAB)]
  ring

/-- **Valuation-dependent.** Positive homogeneity of the positive-cone
extension. -/
theorem posValue_smul (hv : IsEffectValuation v)
    {A : Matrix (Fin d) (Fin d) ℂ} (hA : A.PosSemidef) {c : ℝ}
    (hc : 0 ≤ c) : posValue v (c • A) = c * posValue v A := by
  rcases eq_or_lt_of_le hc with rfl | hc'
  · have h0 : (0 : ℝ) • A = (0 : Matrix (Fin d) (Fin d) ℂ) := by module
    have h0' : ∀ t : ℝ, t • (0 : Matrix (Fin d) (Fin d) ℂ) = 0 := by
      intro t
      module
    rw [h0]
    simp [posValue, h0', hv.map_zero]
  · have hs := posBound_pos hA
    have hbound : ((c * posBound A) • (1 : Matrix (Fin d) (Fin d) ℂ)
        - c • A).PosSemidef := by
      have key : (c * posBound A) • (1 : Matrix (Fin d) (Fin d) ℂ) - c • A
          = c • (posBound A • 1 - A) := by module
      rw [key]
      exact posSemidef_real_smul (posBound_smul_one_sub hA) hc
    rw [posValue_eq_of_bound hv (posSemidef_real_smul hA hc)
      (mul_pos hc' hs) hbound]
    have harg : (c * posBound A)⁻¹ • (c • A) = (posBound A)⁻¹ • A := by
      rw [real_smul_smul]
      congr 1
      field_simp
    rw [harg]
    unfold posValue
    ring

/-- **Valuation-dependent.** The positive-cone extension is
nonnegative. -/
theorem posValue_nonneg (hv : IsEffectValuation v)
    {A : Matrix (Fin d) (Fin d) ℂ} (hA : A.PosSemidef) :
    0 ≤ posValue v A :=
  mul_nonneg (posBound_pos hA).le (hv.nonneg (isEffect_posBound_inv_smul hA))

/-! ## The extension to Hermitian observables -/

/-- **Algebra-only.** Real scalings preserve hermiticity. -/
theorem isHermitian_real_smul {H : Matrix (Fin d) (Fin d) ℂ}
    (hH : H.IsHermitian) (c : ℝ) : (c • H).IsHermitian :=
  hH.smul (IsSelfAdjoint.all c)

open Classical in
/-- The chosen shift of a matrix into the positive cone; for Hermitian
matrices `H + hermShift H • 1` is positive semidefinite
(`hermShift_posSemidef`). -/
noncomputable def hermShift (H : Matrix (Fin d) (Fin d) ℂ) : ℝ :=
  if h : ∃ s : ℝ, (H + s • (1 : Matrix (Fin d) (Fin d) ℂ)).PosSemidef then
    h.choose
  else 0

/-- **Algebra-only.** The chosen shift works for Hermitian matrices. -/
theorem hermShift_posSemidef {H : Matrix (Fin d) (Fin d) ℂ}
    (hH : H.IsHermitian) :
    (H + hermShift H • (1 : Matrix (Fin d) (Fin d) ℂ)).PosSemidef := by
  have h := exists_add_smul_one_posSemidef hH
  rw [hermShift, dif_pos h]
  exact h.choose_spec

/-- The extension of an effect valuation to Hermitian observables: shift
into the positive cone, extend, and subtract the shift. In the language of
step (4), `H = t • E - s • 1` with `s = hermShift H`,
`t = posBound (H + s • 1)`, and `E` the canonical effect of the shifted
matrix; `hermValue v H = t * v E - s`. -/
noncomputable def hermValue (v : Matrix (Fin d) (Fin d) ℂ → ℝ)
    (H : Matrix (Fin d) (Fin d) ℂ) : ℝ :=
  posValue v (H + hermShift H • 1) - hermShift H

/-- **Valuation-dependent.** Shifting a positive matrix by a nonnegative
multiple of the identity adds that multiple to the positive-cone
extension. -/
theorem posValue_add_smul_one (hv : IsEffectValuation v)
    {A : Matrix (Fin d) (Fin d) ℂ} (hA : A.PosSemidef) {c : ℝ}
    (hc : 0 ≤ c) :
    posValue v (A + c • 1) = posValue v A + c := by
  have hone : posValue v (1 : Matrix (Fin d) (Fin d) ℂ) = 1 := by
    rw [posValue_of_isEffect hv isEffect_one, hv.map_one]
  have hcone : posValue v (c • (1 : Matrix (Fin d) (Fin d) ℂ)) = c := by
    rw [posValue_smul hv Matrix.PosSemidef.one hc, hone, mul_one]
  rw [posValue_add hv hA (posSemidef_real_smul Matrix.PosSemidef.one hc),
    hcone]

/-- **Valuation-dependent.** Well-definedness of the Hermitian extension:
any shift into the positive cone computes the same value. Two shifts are
compared through their common refinement, using only additivity and
homogeneity of the positive-cone extension. -/
theorem hermValue_eq_of_shift (hv : IsEffectValuation v)
    {H : Matrix (Fin d) (Fin d) ℂ} (hH : H.IsHermitian) {s : ℝ}
    (hs : (H + s • (1 : Matrix (Fin d) (Fin d) ℂ)).PosSemidef) :
    hermValue v H = posValue v (H + s • 1) - s := by
  have hσ := hermShift_posSemidef hH
  have hsu : 0 ≤ s + (|s| + |hermShift H|) := by
    have := neg_abs_le s
    have := abs_nonneg (hermShift H)
    linarith
  have hσu : 0 ≤ hermShift H + (|s| + |hermShift H|) := by
    have := neg_abs_le (hermShift H)
    have := abs_nonneg s
    linarith
  have hkey : (H + s • 1) +
      (hermShift H + (|s| + |hermShift H|)) • (1 : Matrix (Fin d) (Fin d) ℂ)
      = (H + hermShift H • 1) + (s + (|s| + |hermShift H|)) • 1 := by
    module
  have h1 := posValue_add_smul_one hv hs hσu
  have h2 := posValue_add_smul_one hv hσ hsu
  rw [hkey, h2] at h1
  unfold hermValue
  linarith

/-- **Valuation-dependent.** The Hermitian extension vanishes at zero. -/
theorem hermValue_zero (hv : IsEffectValuation v) :
    hermValue v (0 : Matrix (Fin d) (Fin d) ℂ) = 0 := by
  have key : (0 : Matrix (Fin d) (Fin d) ℂ) + (0 : ℝ) • 1 = 0 := by module
  have h0 : ((0 : Matrix (Fin d) (Fin d) ℂ) +
      (0 : ℝ) • (1 : Matrix (Fin d) (Fin d) ℂ)).PosSemidef := by
    rw [key]
    exact Matrix.PosSemidef.zero
  rw [hermValue_eq_of_shift hv isHermitian_zero h0, key,
    posValue_of_isEffect hv isEffect_zero, hv.map_zero]
  ring

/-- **Valuation-dependent.** The Hermitian extension restricts to the
valuation on effects. -/
theorem hermValue_of_isEffect (hv : IsEffectValuation v)
    {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E) :
    hermValue v E = v E := by
  have key : E + (0 : ℝ) • (1 : Matrix (Fin d) (Fin d) ℂ) = E := by module
  have h0 : (E + (0 : ℝ) • (1 : Matrix (Fin d) (Fin d) ℂ)).PosSemidef := by
    rw [key]
    exact hE.1
  rw [hermValue_eq_of_shift hv hE.isHermitian h0, key,
    posValue_of_isEffect hv hE]
  ring

/-- **Valuation-dependent.** Additivity of the Hermitian extension. -/
theorem hermValue_add (hv : IsEffectValuation v)
    {H K : Matrix (Fin d) (Fin d) ℂ} (hH : H.IsHermitian)
    (hK : K.IsHermitian) :
    hermValue v (H + K) = hermValue v H + hermValue v K := by
  obtain ⟨s, hs⟩ := exists_add_smul_one_posSemidef hH
  obtain ⟨t, ht⟩ := exists_add_smul_one_posSemidef hK
  have key : (H + K) + (s + t) • (1 : Matrix (Fin d) (Fin d) ℂ)
      = (H + s • 1) + (K + t • 1) := by module
  have hsum : ((H + K) +
      (s + t) • (1 : Matrix (Fin d) (Fin d) ℂ)).PosSemidef := by
    rw [key]
    exact hs.add ht
  rw [hermValue_eq_of_shift hv (hH.add hK) hsum,
    hermValue_eq_of_shift hv hH hs, hermValue_eq_of_shift hv hK ht, key,
    posValue_add hv hs ht]
  ring

/-- **Valuation-dependent.** The Hermitian extension is odd. -/
theorem hermValue_neg (hv : IsEffectValuation v)
    {H : Matrix (Fin d) (Fin d) ℂ} (hH : H.IsHermitian) :
    hermValue v (-H) = - hermValue v H := by
  have h := hermValue_add hv hH hH.neg
  rw [add_neg_cancel, hermValue_zero hv] at h
  linarith

/-- **Valuation-dependent.** Nonnegative real homogeneity of the
Hermitian extension. -/
theorem hermValue_smul_of_nonneg (hv : IsEffectValuation v)
    {H : Matrix (Fin d) (Fin d) ℂ} (hH : H.IsHermitian) {c : ℝ}
    (hc : 0 ≤ c) :
    hermValue v (c • H) = c * hermValue v H := by
  obtain ⟨s, hs⟩ := exists_add_smul_one_posSemidef hH
  have key : (c • H) + (c * s) • (1 : Matrix (Fin d) (Fin d) ℂ)
      = c • (H + s • 1) := by module
  have hcs : ((c • H) +
      (c * s) • (1 : Matrix (Fin d) (Fin d) ℂ)).PosSemidef := by
    rw [key]
    exact posSemidef_real_smul hs hc
  rw [hermValue_eq_of_shift hv (isHermitian_real_smul hH c) hcs, key,
    posValue_smul hv hs hc, hermValue_eq_of_shift hv hH hs]
  ring

/-- **Valuation-dependent.** Real homogeneity of the Hermitian
extension. -/
theorem hermValue_smul (hv : IsEffectValuation v)
    {H : Matrix (Fin d) (Fin d) ℂ} (hH : H.IsHermitian) (c : ℝ) :
    hermValue v (c • H) = c * hermValue v H := by
  rcases le_total 0 c with hc | hc
  · exact hermValue_smul_of_nonneg hv hH hc
  · have key : c • H = (-c) • (-H) := by module
    rw [key, hermValue_smul_of_nonneg hv hH.neg (by linarith),
      hermValue_neg hv hH]
    ring

/-- **Valuation-dependent.** The Hermitian extension commutes with finite
sums of Hermitian matrices. -/
theorem hermValue_sum (hv : IsEffectValuation v) {ι : Type*}
    (s : Finset ι) (f : ι → Matrix (Fin d) (Fin d) ℂ)
    (hf : ∀ i ∈ s, (f i).IsHermitian) :
    hermValue v (∑ i ∈ s, f i) = ∑ i ∈ s, hermValue v (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using hermValue_zero hv
  | insert a s ha ih =>
    have hha : (f a).IsHermitian := hf a (Finset.mem_insert_self a s)
    have hhs : ∀ i ∈ s, (f i).IsHermitian := fun i hi =>
      hf i (Finset.mem_insert_of_mem hi)
    have hsum : (∑ i ∈ s, f i).IsHermitian :=
      Finset.sum_induction f _ (fun x y hx hy => hx.add hy)
        isHermitian_zero hhs
    rw [Finset.sum_insert ha, hermValue_add hv hha hsum, ih hhs,
      Finset.sum_insert ha]

/-! ## The Hermitian matrix units and the finite Riesz representation -/

/-- The symmetrised matrix unit `E_{jk} + E_{kj}`: the real part of the
`(j,k)` matrix entry, as a Hermitian observable. -/
noncomputable def symUnit (j k : Fin d) : Matrix (Fin d) (Fin d) ℂ :=
  single j k 1 + single k j 1

/-- The antisymmetrised matrix unit `i E_{jk} - i E_{kj}`: the imaginary
part of the `(j,k)` matrix entry, as a Hermitian observable. -/
noncomputable def skewUnit (j k : Fin d) : Matrix (Fin d) (Fin d) ℂ :=
  Complex.I • single j k 1 - Complex.I • single k j 1

/-- **Algebra-only.** The symmetrised units are Hermitian. -/
theorem symUnit_isHermitian (j k : Fin d) : (symUnit j k).IsHermitian := by
  unfold symUnit Matrix.IsHermitian
  rw [conjTranspose_add, conjTranspose_single, conjTranspose_single]
  simp [add_comm]

/-- **Algebra-only.** The antisymmetrised units are Hermitian. -/
theorem skewUnit_isHermitian (j k : Fin d) : (skewUnit j k).IsHermitian := by
  unfold skewUnit Matrix.IsHermitian
  rw [conjTranspose_sub, conjTranspose_smul, conjTranspose_smul,
    conjTranspose_single, conjTranspose_single]
  simp [Complex.conj_I, neg_smul]
  module

/-- **Algebra-only.** Index symmetry of the symmetrised units. -/
theorem symUnit_swap (j k : Fin d) : symUnit k j = symUnit j k := by
  unfold symUnit
  exact add_comm _ _

/-- **Algebra-only.** Index antisymmetry of the antisymmetrised units. -/
theorem skewUnit_swap (j k : Fin d) : skewUnit k j = - skewUnit j k := by
  unfold skewUnit
  rw [neg_sub]

/-- **Trace-dependent.** The trace pairing against a symmetrised unit
reads off the symmetrised matrix entry. -/
theorem trace_mul_symUnit (M : Matrix (Fin d) (Fin d) ℂ) (j k : Fin d) :
    (M * symUnit j k).trace = M k j + M j k := by
  unfold symUnit
  rw [mul_add, trace_add, trace_mul_single, trace_mul_single]
  simp

/-- **Trace-dependent.** The trace pairing against an antisymmetrised
unit reads off the antisymmetrised matrix entry. -/
theorem trace_mul_skewUnit (M : Matrix (Fin d) (Fin d) ℂ) (j k : Fin d) :
    (M * skewUnit j k).trace = Complex.I * M k j - Complex.I * M j k := by
  unfold skewUnit
  rw [mul_sub, trace_sub, Matrix.mul_smul, Matrix.mul_smul, trace_smul,
    trace_smul, trace_mul_single, trace_mul_single]
  simp

/-- **Algebra-only.** Every Hermitian matrix is the real span of the
Hermitian matrix units: the finite Riesz expansion behind step (5). Each
unordered index pair is counted twice, once from each order, whence the
factor `1/2`. -/
theorem isHermitian_eq_sum {H : Matrix (Fin d) (Fin d) ℂ}
    (hH : H.IsHermitian) :
    H = (1/2 : ℝ) • ∑ j, ∑ k,
      ((H j k).re • symUnit j k + (H j k).im • skewUnit j k) := by
  ext p q
  have hqp : H q p = star (H p q) := by
    rw [← hH.apply p q, star_star]
  have hre : (H q p).re = (H p q).re := by
    rw [hqp]
    simp
  have him : (H q p).im = -(H p q).im := by
    rw [hqp]
    simp
  rw [Matrix.smul_apply, Matrix.sum_apply]
  simp only [Matrix.sum_apply, Matrix.add_apply, Matrix.smul_apply, symUnit,
    skewUnit, Matrix.sub_apply, Matrix.single_apply, Complex.real_smul,
    smul_eq_mul]
  simp only [mul_add, mul_sub, mul_ite, mul_one, mul_zero, ite_and,
    Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_ite_irrel,
    Finset.sum_const_zero, Finset.sum_ite_eq, Finset.sum_ite_eq',
    Finset.mem_univ, if_true]
  rw [hre, him]
  conv_lhs => rw [← Complex.re_add_im (H p q)]
  push_cast
  ring

/-- The candidate density matrix of the valuation, read off from the
Hermitian extension on the matrix units: the diagonal from the
symmetrised units, the off-diagonal phases from both unit families. -/
noncomputable def gleasonState (v : Matrix (Fin d) (Fin d) ℂ → ℝ) :
    Matrix (Fin d) (Fin d) ℂ :=
  Matrix.of fun j k =>
    (((hermValue v (symUnit j k) : ℝ) : ℂ) +
      Complex.I * ((hermValue v (skewUnit j k) : ℝ) : ℂ)) / 2

theorem gleasonState_apply (v : Matrix (Fin d) (Fin d) ℂ → ℝ)
    (j k : Fin d) :
    gleasonState v j k =
      (((hermValue v (symUnit j k) : ℝ) : ℂ) +
        Complex.I * ((hermValue v (skewUnit j k) : ℝ) : ℂ)) / 2 :=
  rfl

/-- **Trace-dependent.** The candidate state reproduces the extension on
the symmetrised units. -/
theorem trace_gleasonState_mul_symUnit (hv : IsEffectValuation v)
    (j k : Fin d) :
    (gleasonState v * symUnit j k).trace =
      ((hermValue v (symUnit j k) : ℝ) : ℂ) := by
  rw [trace_mul_symUnit, gleasonState_apply, gleasonState_apply,
    symUnit_swap j k, skewUnit_swap j k,
    hermValue_neg hv (skewUnit_isHermitian j k)]
  push_cast
  ring

/-- **Trace-dependent.** The candidate state reproduces the extension on
the antisymmetrised units. -/
theorem trace_gleasonState_mul_skewUnit (hv : IsEffectValuation v)
    (j k : Fin d) :
    (gleasonState v * skewUnit j k).trace =
      ((hermValue v (skewUnit j k) : ℝ) : ℂ) := by
  rw [trace_mul_skewUnit, gleasonState_apply, gleasonState_apply,
    symUnit_swap j k, skewUnit_swap j k,
    hermValue_neg hv (skewUnit_isHermitian j k)]
  push_cast
  linear_combination
    (-((hermValue v (skewUnit j k) : ℝ) : ℂ)) * Complex.I_mul_I

/-- **Trace-dependent.** The finite Riesz representation, step (5): the
trace pairing against the candidate state computes the Hermitian
extension of the valuation on every Hermitian observable. -/
theorem trace_gleasonState_mul_isHermitian (hv : IsEffectValuation v)
    {H : Matrix (Fin d) (Fin d) ℂ} (hH : H.IsHermitian) :
    (gleasonState v * H).trace = ((hermValue v H : ℝ) : ℂ) := by
  set W : Matrix (Fin d) (Fin d) ℂ :=
    ∑ j, ∑ k, ((H j k).re • symUnit j k + (H j k).im • skewUnit j k)
    with hW
  have hterm : ∀ j k : Fin d,
      ((H j k).re • symUnit j k + (H j k).im • skewUnit j k).IsHermitian :=
    fun j k => (isHermitian_real_smul (symUnit_isHermitian j k) _).add
      (isHermitian_real_smul (skewUnit_isHermitian j k) _)
  have hrow : ∀ j : Fin d,
      (∑ k, ((H j k).re • symUnit j k + (H j k).im • skewUnit j k)).IsHermitian :=
    fun j => Finset.sum_induction _ _ (fun x y hx hy => hx.add hy)
      isHermitian_zero fun k _ => hterm j k
  have hWherm : W.IsHermitian := by
    rw [hW]
    exact Finset.sum_induction _ _ (fun x y hx hy => hx.add hy)
      isHermitian_zero fun j _ => hrow j
  have htermtrace : ∀ j k : Fin d,
      (gleasonState v *
        ((H j k).re • symUnit j k + (H j k).im • skewUnit j k)).trace =
      ((hermValue v
        ((H j k).re • symUnit j k + (H j k).im • skewUnit j k) : ℝ) : ℂ) := by
    intro j k
    rw [mul_add, mul_real_smul, mul_real_smul, trace_add, trace_real_smul,
      trace_real_smul, trace_gleasonState_mul_symUnit hv,
      trace_gleasonState_mul_skewUnit hv,
      hermValue_add hv (isHermitian_real_smul (symUnit_isHermitian j k) _)
        (isHermitian_real_smul (skewUnit_isHermitian j k) _),
      hermValue_smul hv (symUnit_isHermitian j k),
      hermValue_smul hv (skewUnit_isHermitian j k)]
    push_cast
    ring
  have hsum : (gleasonState v * W).trace = ((hermValue v W : ℝ) : ℂ) := by
    rw [hW, Finset.mul_sum, trace_sum,
      hermValue_sum hv _ _ fun j _ => hrow j]
    push_cast
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Finset.mul_sum, trace_sum,
      hermValue_sum hv _ _ fun k _ => hterm j k]
    push_cast
    exact Finset.sum_congr rfl fun k _ => htermtrace j k
  have hdecomp : H = (1/2 : ℝ) • W := by
    rw [hW]
    exact isHermitian_eq_sum hH
  rw [hdecomp, mul_real_smul, trace_real_smul, hermValue_smul hv hWherm,
    hsum]
  push_cast
  ring

/-! ## State properties of the represented matrix -/

/-- **Trace-dependent.** The candidate state is Hermitian: index symmetry
of the symmetrised units and antisymmetry of the antisymmetrised units
conjugate the entry formula into itself. -/
theorem gleasonState_isHermitian (hv : IsEffectValuation v) :
    (gleasonState v).IsHermitian := by
  show (gleasonState v)ᴴ = gleasonState v
  ext i j
  rw [Matrix.conjTranspose_apply, gleasonState_apply, gleasonState_apply,
    symUnit_swap i j, skewUnit_swap i j,
    hermValue_neg hv (skewUnit_isHermitian i j)]
  push_cast
  simp only [star_div₀, star_add, star_mul', star_ofNat, RCLike.star_def,
    Complex.conj_ofReal, Complex.conj_I, map_neg]
  ring

/-- **Trace-dependent.** The quadratic form of the trace pairing: for any
matrix `M` and vector `x`, `x⋆ M x = Tr (M · x x⋆)`. -/
theorem dotProduct_mulVec_eq_trace (M : Matrix (Fin d) (Fin d) ℂ)
    (x : Fin d → ℂ) :
    star x ⬝ᵥ (M *ᵥ x) = (M * vecMulVec x (star x)).trace := by
  simp only [dotProduct, mulVec, Matrix.trace, Matrix.diag, Matrix.mul_apply,
    vecMulVec_apply, Pi.star_apply, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun l _ => ?_
  ring

/-- **Algebra-only.** The outer square identity
`(x x⋆) (x x⋆) = ‖x‖² (x x⋆)`. -/
theorem vecMulVec_mul_self (x : Fin d → ℂ) :
    vecMulVec x (star x) * vecMulVec x (star x) =
      (∑ i, Complex.normSq (x i)) • vecMulVec x (star x) := by
  ext i j
  simp only [Matrix.mul_apply, vecMulVec_apply, Pi.star_apply,
    Matrix.smul_apply, Complex.real_smul]
  push_cast
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [← Complex.mul_conj]
  ring

/-- **Algebra-only.** The outer product is Hermitian. -/
theorem vecMulVec_star_isHermitian (x : Fin d → ℂ) :
    (vecMulVec x (star x)).IsHermitian := by
  show (vecMulVec x (star x))ᴴ = vecMulVec x (star x)
  ext i j
  simp only [Matrix.conjTranspose_apply, vecMulVec_apply, Pi.star_apply,
    star_mul', star_star]
  ring

/-- **Algebra-only.** The normalised outer product of a nonzero vector is
a projection event. -/
theorem isEvent_normalized_vecMulVec {x : Fin d → ℂ} (hx : x ≠ 0) :
    IsEvent ((∑ i, Complex.normSq (x i))⁻¹ • vecMulVec x (star x)) := by
  have hr : 0 < ∑ i, Complex.normSq (x i) := by
    obtain ⟨i, hi⟩ := Function.ne_iff.mp hx
    refine Finset.sum_pos' (fun j _ => Complex.normSq_nonneg _) ?_
    exact ⟨i, Finset.mem_univ i, Complex.normSq_pos.mpr hi⟩
  constructor
  · exact isHermitian_real_smul (vecMulVec_star_isHermitian x) _
  · rw [real_smul_mul, mul_real_smul, vecMulVec_mul_self, real_smul_smul,
      real_smul_smul]
    congr 1
    field_simp

/-- **Trace-dependent.** Positivity, step (6): the value on normalised
rank-one projectors is nonnegative, so the quadratic form of the
candidate state is nonnegative. -/
theorem gleasonState_posSemidef (hv : IsEffectValuation v) :
    (gleasonState v).PosSemidef := by
  refine PosSemidef.of_dotProduct_mulVec_nonneg (gleasonState_isHermitian hv)
    fun x => ?_
  rcases eq_or_ne x 0 with rfl | hx
  · simp
  · set r : ℝ := ∑ i, Complex.normSq (x i) with hr
    have hrpos : 0 < r := by
      obtain ⟨i, hi⟩ := Function.ne_iff.mp hx
      refine Finset.sum_pos' (fun j _ => Complex.normSq_nonneg _) ?_
      exact ⟨i, Finset.mem_univ i, Complex.normSq_pos.mpr hi⟩
    have hQ : IsEvent (r⁻¹ • vecMulVec x (star x)) :=
      isEvent_normalized_vecMulVec hx
    have hval := trace_gleasonState_mul_isHermitian hv hQ.1
    rw [hermValue_of_isEffect hv hQ.isEffect] at hval
    have hsplit : vecMulVec x (star x) =
        r • (r⁻¹ • vecMulVec x (star x)) := by
      rw [real_smul_smul, mul_inv_cancel₀ hrpos.ne']
      module
    rw [dotProduct_mulVec_eq_trace, hsplit, mul_real_smul, trace_real_smul,
      hval, ← Complex.ofReal_mul]
    exact Complex.zero_le_real.mpr
      (mul_nonneg hrpos.le (hv.nonneg hQ.isEffect))

/-- **Trace-dependent.** Normalisation, step (6): the candidate state has
unit trace because `v 1 = 1`. -/
theorem gleasonState_trace (hv : IsEffectValuation v) :
    (gleasonState v).trace = 1 := by
  have h := trace_gleasonState_mul_isHermitian hv isHermitian_one
  rw [mul_one, hermValue_of_isEffect hv isEffect_one, hv.map_one] at h
  simpa using h

/-- **Trace-dependent.** The candidate is a state. -/
theorem gleasonState_isState (hv : IsEffectValuation v) :
    IsState (gleasonState v) :=
  ⟨gleasonState_posSemidef hv, gleasonState_trace hv⟩

/-- **Trace-dependent.** The candidate state represents the valuation on
every effect through the Born weight. -/
theorem gleasonState_represents (hv : IsEffectValuation v)
    {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E) :
    v E = (bornWeight (gleasonState v) E).re := by
  have h := trace_gleasonState_mul_isHermitian hv hE.isHermitian
  rw [hermValue_of_isEffect hv hE] at h
  rw [bornWeight, h, Complex.ofReal_re]

/-! ## Uniqueness -/

/-- **Algebra-only.** The effect decomposition of step (4), packaged: any
Hermitian matrix is a nonnegative scaling of an effect minus a real
multiple of the identity. -/
theorem exists_effect_decomposition {H : Matrix (Fin d) (Fin d) ℂ}
    (hH : H.IsHermitian) :
    ∃ (t s : ℝ) (E : Matrix (Fin d) (Fin d) ℂ),
      IsEffect E ∧ 0 ≤ t ∧ H = t • E - s • 1 := by
  obtain ⟨s, hs⟩ := exists_add_smul_one_posSemidef hH
  refine ⟨posBound (H + s • 1), s, (posBound (H + s • 1))⁻¹ • (H + s • 1),
    isEffect_posBound_inv_smul hs, (posBound_pos hs).le, ?_⟩
  rw [real_smul_smul, mul_inv_cancel₀ (posBound_pos hs).ne']
  module

/-- **Trace-dependent.** Two states with the same Born weights on all
effects agree on all Hermitian observables, hence on the Hermitian matrix
units, hence entrywise: the uniqueness half of the theorem. -/
theorem state_ext_of_forall_effect_re {ρ ρ' : Matrix (Fin d) (Fin d) ℂ}
    (hρ : IsState ρ) (hρ' : IsState ρ')
    (h : ∀ E, IsEffect E → (bornWeight ρ E).re = (bornWeight ρ' E).re) :
    ρ = ρ' := by
  have hC : ∀ E, IsEffect E → bornWeight ρ E = bornWeight ρ' E := by
    intro E hE
    rw [← bornWeight_eq_re hρ.1.1 hE.isHermitian,
      ← bornWeight_eq_re hρ'.1.1 hE.isHermitian, h E hE]
  have hsm : ∀ (σ : Matrix (Fin d) (Fin d) ℂ) (c : ℝ)
      (M : Matrix (Fin d) (Fin d) ℂ),
      bornWeight σ (c • M) = (c : ℂ) * bornWeight σ M := by
    intro σ c M
    rw [bornWeight, bornWeight, mul_real_smul, trace_real_smul]
  have hHerm : ∀ H : Matrix (Fin d) (Fin d) ℂ, H.IsHermitian →
      bornWeight ρ H = bornWeight ρ' H := by
    intro H hHh
    obtain ⟨t, s, E, hE, ht, hdec⟩ := exists_effect_decomposition hHh
    rw [hdec, bornWeight_sub, bornWeight_sub, hsm, hsm, hsm, hsm,
      hC E hE, bornWeight_one hρ, bornWeight_one hρ']
  ext j k
  have hS := hHerm (symUnit j k) (symUnit_isHermitian j k)
  have hA := hHerm (skewUnit j k) (skewUnit_isHermitian j k)
  rw [bornWeight, bornWeight, trace_mul_symUnit, trace_mul_symUnit] at hS
  rw [bornWeight, bornWeight, trace_mul_skewUnit, trace_mul_skewUnit,
    ← mul_sub, ← mul_sub] at hA
  have hA' := mul_left_cancel₀ Complex.I_ne_zero hA
  have h2 : (2 : ℂ) * ρ j k = 2 * ρ' j k := by
    linear_combination hS - hA'
  exact mul_left_cancel₀ two_ne_zero h2

/-! ## The representation theorem -/

/-- **Trace-dependent.** The finite Busch--Gleason representation
theorem: every additive effect valuation on `Matrix (Fin d) (Fin d) ℂ`
is the Born functional `E ↦ Re (Tr (ρ E))` of a unique density
matrix `ρ`. -/
theorem finite_busch_gleason (hv : IsEffectValuation v) :
    ∃! ρ : Matrix (Fin d) (Fin d) ℂ,
      IsState ρ ∧ ∀ E, IsEffect E → v E = (bornWeight ρ E).re := by
  refine ⟨gleasonState v,
    ⟨gleasonState_isState hv, fun E hE => gleasonState_represents hv hE⟩,
    ?_⟩
  rintro ρ ⟨hρ, hrep⟩
  exact state_ext_of_forall_effect_re hρ (gleasonState_isState hv)
    fun E hE => by rw [← hrep E hE, gleasonState_represents hv hE]

/-- **Trace-dependent.** The easy converse: the Born functional of any
state is an additive effect valuation. Nonnegativity is positivity of the
expectation functional; the upper bound passes through the complement
effect. -/
theorem isEffectValuation_born {ρ : Matrix (Fin d) (Fin d) ℂ}
    (hρ : IsState ρ) :
    IsEffectValuation (fun E => (bornWeight ρ E).re) := by
  constructor
  · intro E hE
    have h := expectation_nonneg hρ.1 hE.1
    simpa using (Complex.le_def.mp h).1
  · intro E hE
    have h := expectation_nonneg hρ.1 hE.2
    have hsub : expectation ρ (1 - E) = 1 - bornWeight ρ E := by
      rw [← bornWeight_eq_expectation, bornWeight_sub, bornWeight_one hρ]
    rw [hsub] at h
    have := (Complex.le_def.mp h).1
    simpa using this
  · rw [bornWeight_one hρ]
    simp
  · intro E F hE hF hEF
    rw [bornWeight_add]
    simp

/-- **Trace-dependent.** The exact characterisation: the additive effect
valuations are precisely the Born functionals of density matrices. This
is the mathematical half of the Born derivation; Born weights are exactly
the noncontextually additive effect assignments. -/
theorem isEffectValuation_iff_born :
    IsEffectValuation v ↔
      ∃ ρ, IsState ρ ∧ ∀ E, IsEffect E → v E = (bornWeight ρ E).re := by
  constructor
  · intro hv
    obtain ⟨ρ, hρ, -⟩ := finite_busch_gleason hv
    exact ⟨ρ, hρ⟩
  · rintro ⟨ρ, hρ, hrep⟩
    have hb := isEffectValuation_born hρ
    constructor
    · intro E hE
      rw [hrep E hE]
      exact hb.nonneg hE
    · intro E hE
      rw [hrep E hE]
      exact hb.le_one hE
    · rw [hrep 1 isEffect_one]
      exact hb.map_one
    · intro E F hE hF hEF
      rw [hrep _ hEF, hrep _ hE, hrep _ hF]
      exact hb.additive hE hF hEF

/-- **Trace-dependent.** Corollary for projective decompositions: any
finite family of projection events summing to the identity receives total
value one, with each weight of Born form under the represented state. -/
theorem sum_valuation_projective_decomposition (hv : IsEffectValuation v)
    {ι : Type*} (s : Finset ι) (P : ι → Matrix (Fin d) (Fin d) ℂ)
    (hP : ∀ i ∈ s, IsEvent (P i)) (hsum : ∑ i ∈ s, P i = 1) :
    ∑ i ∈ s, v (P i) = 1 := by
  obtain ⟨ρ, ⟨hρ, hrep⟩, -⟩ := finite_busch_gleason hv
  calc ∑ i ∈ s, v (P i)
      = ∑ i ∈ s, (bornWeight ρ (P i)).re :=
        Finset.sum_congr rfl fun i hi => hrep _ (hP i hi).isEffect
    _ = (∑ i ∈ s, bornWeight ρ (P i)).re := by rw [Complex.re_sum]
    _ = (bornWeight ρ (∑ i ∈ s, P i)).re := by rw [bornWeight_sum]
    _ = 1 := by
        rw [hsum, bornWeight_one hρ]
        simp

end EventAlgebra
