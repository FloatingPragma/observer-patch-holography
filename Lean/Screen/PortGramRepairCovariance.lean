import Mathlib
import PortGramRepairBand

open scoped BigOperators Matrix
open Filter Topology

namespace OPH.PortGramRepairCovariance

open OPH.A5FamilyBand

/-!
# Exact diffusion kernel of the twelve-port repair mean

This file separates a finite deterministic calculation from its possible
physical interpretation.  The registered adjacency `A` supplies the declared
one-step conservative repair mean

`T = I - (5 I - A) / 60`.

On the nonconstant sector `Q = I - P0`, the exact kernel
`Q T^(2 n) Q` is the sum of three repair bands.  Their response coefficients
are strictly ordered, and normalization by the slow coefficient converges
entrywise to `4 P3`, which is the pinned unit-diagonal port Gram matrix.

Everything through the finite-step formula is exact arithmetic in
`Q(sqrt(5))`.  The final limit is an ordinary real, finite-dimensional
analytic statement.  Neither statement selects this repair mean from the
axioms, identifies its completed readback topology with physical space,
chooses a physical field sector or clock, nor supplies a laboratory
observable.  At every finite step all three nonconstant response coefficients
remain positive; the rank-three form appears only after the normalized
infinite-step limit.  Any quotient/completion by `G` therefore applies to the
limiting carrier-position readback, while record order and repair cost remain
separate.  Such a completion is one abstract local carrier, canonical only up
to port-label-preserving isometry.  A faithful isometric action on it and
overlap/refinement gluing are not formalized here.

TRUST BOUNDARY.  Seven large exact rational table identities in this file are
closed with `native_decide`; their printed axiom receipts contain generated
native decision axioms, which place those identities and downstream theorems
that consume them outside the kernel-only trust profile.  Replacing them
requires explicit normalized-rational table certificates or entrywise proofs;
plain kernel `decide` does not reduce these quotient-valued matrix equalities.
-/

/-- Exact quadratic coefficient field represented by `x + y sqrt(5)`. -/
abbrev Q5 := QuadraticAlgebra ℚ 5 0

abbrev PortMatrix := Matrix (Fin 12) (Fin 12) Q5

/-- Exact coefficient constructor. -/
def q5 (x y : ℚ) : Q5 := ⟨x, y⟩

/-- The constant projector `P0`. -/
def pZero : PortMatrix := fun i j ↦
  q5 ((entry M1t i j : ℚ) / 120) 0

/-- The positive rank-three projector `P3`. -/
def pLow : PortMatrix := fun i j ↦
  q5 ((entry Xt i j : ℚ) / 40) ((entry Yt i j : ℚ) / 40)

/-- The even quintet projector `P5`. -/
def pFive : PortMatrix := fun i j ↦
  q5 ((entry M5t i j : ℚ) / 24) 0

/-- The Galois rank-three projector `P3'`. -/
def pHigh : PortMatrix := fun i j ↦
  q5 ((entry Xt i j : ℚ) / 40) (-(entry Yt i j : ℚ) / 40)

/-- Projector onto the complete nonconstant response sector. -/
def qNonconstant : PortMatrix := 1 - pZero

/-- Registered twelve-port adjacency in the exact coefficient field. -/
def adjacency : PortMatrix := fun i j ↦ q5 (entry At i j) 0

/-- Mean of one uniformly sampled conservative seam repair. -/
def repairMean : PortMatrix :=
  1 - q5 (1 / 60) 0 • (q5 5 0 • (1 : PortMatrix) - adjacency)

/-- Response on the positive rank-three band. -/
def lowResponse : Q5 := q5 (55 / 60) (1 / 60)

/-- Response on the quintet. -/
def fiveResponse : Q5 := q5 (9 / 10) 0

/-- Response on the Galois rank-three band. -/
def highResponse : Q5 := q5 (55 / 60) (-1 / 60)

set_option maxHeartbeats 4000000 in
/-- The four pinned tables are exact projectors. -/
@[simp] theorem projector_squares :
    pZero * pZero = pZero ∧ pLow * pLow = pLow ∧
      pFive * pFive = pFive ∧ pHigh * pHigh = pHigh := by
  native_decide

set_option maxHeartbeats 4000000 in
/-- Every two distinct pinned bands are mutually orthogonal, in both orders. -/
@[simp] theorem projector_cross_products :
    pZero * pLow = 0 ∧ pLow * pZero = 0 ∧
      pZero * pFive = 0 ∧ pFive * pZero = 0 ∧
      pZero * pHigh = 0 ∧ pHigh * pZero = 0 ∧
      pLow * pFive = 0 ∧ pFive * pLow = 0 ∧
      pLow * pHigh = 0 ∧ pHigh * pLow = 0 ∧
      pFive * pHigh = 0 ∧ pHigh * pFive = 0 := by
  native_decide

set_option maxHeartbeats 4000000 in
/-- The four exact bands resolve the identity. -/
@[simp] theorem projector_resolution :
    pZero + pLow + pFive + pHigh = 1 := by
  native_decide

/-- A four-band spectral combination. -/
def bandSum (a0 a3 a5 a3p : Q5) : PortMatrix :=
  a0 • pZero + a3 • pLow + a5 • pFive + a3p • pHigh

set_option maxHeartbeats 4000000 in
/-- The source-backed repair mean has the exact four-band decomposition. -/
theorem repairMean_band_decomposition :
    repairMean = bandSum 1 lowResponse fiveResponse highResponse := by
  native_decide

/-- Multiplication of spectral combinations is coefficientwise. -/
theorem bandSum_mul
    (a0 a3 a5 a3p b0 b3 b5 b3p : Q5) :
    bandSum a0 a3 a5 a3p * bandSum b0 b3 b5 b3p =
      bandSum (a0 * b0) (a3 * b3) (a5 * b5) (a3p * b3p) := by
  rcases projector_squares with ⟨h00, h33, h55, hpp⟩
  rcases projector_cross_products with
    ⟨h03, h30, h05, h50, h0p, hp0, h35, h53, h3p, hp3, h5p, hp5⟩
  simp only [bandSum, Matrix.add_mul, Matrix.mul_add, Matrix.smul_mul,
    Matrix.mul_smul, h00, h33, h55, hpp, h03, h30, h05, h50, h0p,
    hp0, h35, h53, h3p, hp3, h5p, hp5, smul_zero]
  module

/-- Exact finite-step spectral formula for every natural power of `T`. -/
theorem repairMean_pow (n : ℕ) :
    repairMean ^ n =
      bandSum 1 (lowResponse ^ n) (fiveResponse ^ n) (highResponse ^ n) := by
  induction n with
  | zero =>
      simp only [pow_zero, one_smul, bandSum]
      exact projector_resolution.symm
  | succ n ih =>
      rw [pow_succ, ih, repairMean_band_decomposition, bandSum_mul]
      simp only [one_mul, pow_succ]

/-- The nonconstant projector is exactly the sum of the three positive-cost
bands. -/
theorem qNonconstant_band_decomposition :
    qNonconstant = bandSum 0 1 1 1 := by
  unfold qNonconstant bandSum
  rw [← projector_resolution]
  module

/-- Deterministic two-leg diffusion kernel at finite step `n`. -/
def diffusionKernel (n : ℕ) : PortMatrix :=
  qNonconstant * repairMean ^ (2 * n) * qNonconstant

/-- Exact finite-step kernel formula.  No stochastic ensemble or asymptotic
argument enters this theorem. -/
theorem diffusionKernel_exact (n : ℕ) :
    diffusionKernel n =
      bandSum 0 (lowResponse ^ (2 * n)) (fiveResponse ^ (2 * n))
        (highResponse ^ (2 * n)) := by
  rw [diffusionKernel, repairMean_pow, qNonconstant_band_decomposition,
    bandSum_mul, bandSum_mul]
  simp

/-! ## Exact Gram binding and the real asymptotic statement -/

/-- The normalized twelve-port Gram matrix, evaluated in the same exact
quadratic coefficient field. -/
def portGram : PortMatrix := fun i j ↦
  q5 ((OPH.PortFrameGram.g5 i j).1 / 5) ((OPH.PortFrameGram.g5 i j).2 / 5)

set_option maxHeartbeats 4000000 in
/-- Exact source-table identity `G = 4 P3`. -/
theorem portGram_eq_four_pLow :
    portGram = q5 4 0 • pLow := by
  native_decide

/-- The full twelve-port Gram table respects the signed antipodal quotient:
replacing either port by its antipode negates the corresponding row or
column.  Thus the relation `e_(antipode i) = -e_i` is compatible with `G`,
instead of being imposed only after choosing six representatives. -/
theorem portGram_antipodal_quotient :
    (∀ i j : Fin 12,
      OPH.PortFrameGram.g5 (OPH.PortFrameGram.antipode i) j =
        (-(OPH.PortFrameGram.g5 i j).1,
          -(OPH.PortFrameGram.g5 i j).2)) ∧
    (∀ i j : Fin 12,
      OPH.PortFrameGram.g5 i (OPH.PortFrameGram.antipode j) =
        (-(OPH.PortFrameGram.g5 i j).1,
          -(OPH.PortFrameGram.g5 i j).2)) := by
  constructor <;> decide

private theorem sqrt5_sq : Real.sqrt 5 ^ 2 = 5 := by norm_num

/-- The real embedding of the exact quadratic coefficient field. -/
noncomputable def evalReal : Q5 →+* ℝ where
  toFun z := z.re + z.im * Real.sqrt 5
  map_zero' := by simp
  map_one' := by
    simp only [QuadraticAlgebra.re_one, QuadraticAlgebra.im_one]
    norm_num
  map_add' x y := by
    simp only [QuadraticAlgebra.re_add, QuadraticAlgebra.im_add]
    push_cast
    ring
  map_mul' x y := by
    simp only [QuadraticAlgebra.re_mul, QuadraticAlgebra.im_mul]
    push_cast
    ring_nf
    rw [sqrt5_sq]

/-- Entrywise real evaluation of an exact port matrix. -/
noncomputable def realMatrix (M : PortMatrix) : Matrix (Fin 12) (Fin 12) ℝ :=
  M.map evalReal

noncomputable abbrev pZeroR := realMatrix pZero
noncomputable abbrev pLowR := realMatrix pLow
noncomputable abbrev pFiveR := realMatrix pFive
noncomputable abbrev pHighR := realMatrix pHigh

/-! ## Intrinsic local carrier in the twelve-port counting space

The selected carrier can be presented without choosing Cartesian axes.  It is
the range of `P3` acting on the real twelve-port counting space.  The
coordinate frame in `PrimitivePortFrameQuotient` is therefore an isometric
chart and density witness, not part of this definition.
-/

abbrev PortVector := Fin 12 → ℝ

/-- Counting-space scalar product on the twelve declared port labels. -/
noncomputable def portDot (x y : PortVector) : ℝ :=
  ∑ i : Fin 12, x i * y i

/-- The selected projector as a linear endomorphism of port counting space. -/
noncomputable def pLowLinear : PortVector →ₗ[ℝ] PortVector :=
  Matrix.toLin' pLowR

/-- Coordinate-free local carrier: the image of the selected low band inside
the twelve-port counting space. -/
noncomputable def IntrinsicCarrier : Submodule ℝ PortVector :=
  LinearMap.range pLowLinear

/-- Unit record on one declared port label. -/
def portBasis (p : Fin 12) : PortVector :=
  fun i ↦ if i = p then 1 else 0

/-- Intrinsic port generator `v_p = 2 P3 e_p`. -/
noncomputable def intrinsicPortVector (p : Fin 12) : PortVector :=
  pLowLinear (2 • portBasis p)

/-- The intrinsic generator, packaged as a member of the projector range. -/
noncomputable def intrinsicPortGenerator (p : Fin 12) : IntrinsicCarrier :=
  ⟨intrinsicPortVector p, ⟨2 • portBasis p, rfl⟩⟩

set_option maxHeartbeats 4000000 in
/-- The exact low-band table is symmetric. -/
theorem pLow_symmetric : pLow.transpose = pLow := by
  native_decide

/-- Real evaluation preserves the projector identity. -/
theorem pLowR_idempotent : pLowR * pLowR = pLowR := by
  change pLow.map evalReal * pLow.map evalReal = pLow.map evalReal
  rw [← Matrix.map_mul, projector_squares.2.1]

/-- Real evaluation preserves symmetry. -/
theorem pLowR_symmetric : pLowR.transpose = pLowR := by
  have h := congrArg realMatrix pLow_symmetric
  simpa [realMatrix] using h

/-- Each intrinsic generator is twice the corresponding projector column. -/
theorem intrinsicPortVector_apply (p i : Fin 12) :
    intrinsicPortVector p i = 2 * pLowR i p := by
  simp [intrinsicPortVector, pLowLinear, portBasis, Matrix.toLin'_apply,
    Matrix.mulVec, dotProduct]
  ring

set_option maxHeartbeats 4000000 in
/-- Antipodal columns of the selected projector differ by sign. -/
theorem pLow_antipode_columns :
    ∀ i p : Fin 12,
      pLow i (OPH.PortFrameGram.antipode p) = -pLow i p := by
  native_decide

/-- The intrinsic carrier identifies antipodal ports as opposite generators. -/
theorem intrinsicPortVector_antipode (p : Fin 12) :
    intrinsicPortVector (OPH.PortFrameGram.antipode p) =
      -intrinsicPortVector p := by
  funext i
  simp only [Pi.neg_apply]
  rw [intrinsicPortVector_apply, intrinsicPortVector_apply]
  have h := congrArg evalReal (pLow_antipode_columns i p)
  simp only [map_neg] at h
  change pLowR i (OPH.PortFrameGram.antipode p) = -pLowR i p at h
  rw [h]
  ring

theorem intrinsicPortGenerator_antipode (p : Fin 12) :
    intrinsicPortGenerator (OPH.PortFrameGram.antipode p) =
      -intrinsicPortGenerator p := by
  apply Subtype.ext
  exact intrinsicPortVector_antipode p

/-- The intrinsic generators reproduce the selected Gram table without a
Cartesian-coordinate premise. -/
theorem intrinsicPortGenerator_gram (p q : Fin 12) :
    portDot (intrinsicPortGenerator p) (intrinsicPortGenerator q) =
      realMatrix portGram p q := by
  have hsym : ∀ i j, pLowR i j = pLowR j i := by
    intro i j
    have h := congrFun (congrFun pLowR_symmetric i) j
    simpa using h.symm
  have hidem := congrFun (congrFun pLowR_idempotent p) q
  rw [Matrix.mul_apply] at hidem
  simp only [portDot, intrinsicPortGenerator, intrinsicPortVector_apply]
  have hsum :
      (∑ i : Fin 12, pLowR i p * pLowR i q) = pLowR p q := by
    calc
      (∑ i : Fin 12, pLowR i p * pLowR i q) =
          ∑ i : Fin 12, pLowR p i * pLowR i q := by
            apply Finset.sum_congr rfl
            intro i _
            rw [hsym]
      _ = pLowR p q := hidem
  have hgram := congrFun (congrFun (congrArg realMatrix
    portGram_eq_four_pLow) p) q
  have hfour : evalReal (q5 4 0) = (4 : ℝ) := by
    norm_num [evalReal, q5]
  simp only [realMatrix, Matrix.map_apply, Matrix.smul_apply, smul_eq_mul,
    map_mul, hfour] at hgram
  calc
    (∑ x : Fin 12, 2 * pLowR x p * (2 * pLowR x q)) =
        4 * ∑ x : Fin 12, pLowR x p * pLowR x q := by
          rw [Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro x _
          ring
    _ = 4 * pLowR p q := by rw [hsum]
    _ = realMatrix portGram p q := by
      exact hgram.symm

/-- Real four-band spectral combination. -/
noncomputable def realBandSum (a0 a3 a5 a3p : ℝ) :
    Matrix (Fin 12) (Fin 12) ℝ :=
  a0 • pZeroR + a3 • pLowR + a5 • pFiveR + a3p • pHighR

/-- Exact responses in the chosen real embedding. -/
noncomputable def lowResponseR : ℝ := evalReal lowResponse
noncomputable def fiveResponseR : ℝ := evalReal fiveResponse
noncomputable def highResponseR : ℝ := evalReal highResponse

theorem response_real_values :
    lowResponseR = (55 + Real.sqrt 5) / 60 ∧
      fiveResponseR = 9 / 10 ∧
      highResponseR = (55 - Real.sqrt 5) / 60 := by
  constructor
  · norm_num [lowResponseR, evalReal, lowResponse, q5]
    ring
  constructor
  · norm_num [fiveResponseR, evalReal, fiveResponse, q5]
  · norm_num [highResponseR, evalReal, highResponse, q5]
    ring

/-- The slow rank-three response is strictly above the quintet response,
which is strictly above the Galois response; all lie strictly between zero
and one. -/
theorem response_strict_order :
    0 < highResponseR ∧ highResponseR < fiveResponseR ∧
      fiveResponseR < lowResponseR ∧ lowResponseR < 1 := by
  rcases response_real_values with ⟨hlow, hfive, hhigh⟩
  rw [hlow, hfive, hhigh]
  have hsnonneg : 0 ≤ Real.sqrt 5 := Real.sqrt_nonneg 5
  have hs_gt_one : 1 < Real.sqrt 5 := by
    nlinarith [sqrt5_sq]
  have hs_lt_five : Real.sqrt 5 < 5 := by
    nlinarith [sqrt5_sq]
  constructor
  · nlinarith
  constructor
  · nlinarith
  constructor <;> nlinarith

/-- No nonconstant repair band disappears at a finite step.  This is a
coefficient statement; the rank and physical-sector interpretations are
separate. -/
theorem finite_response_coefficients_positive (n : ℕ) :
    0 < lowResponseR ^ (2 * n) ∧
      0 < fiveResponseR ^ (2 * n) ∧
      0 < highResponseR ^ (2 * n) := by
  rcases response_strict_order with ⟨hhigh0, hhighfive, hfivelow, hlow1⟩
  exact ⟨pow_pos (hhigh0.trans (hhighfive.trans hfivelow)) _,
    pow_pos (hhigh0.trans hhighfive) _, pow_pos hhigh0 _⟩

/-- Real evaluation commutes with a spectral combination. -/
theorem realMatrix_bandSum (a0 a3 a5 a3p : Q5) :
    realMatrix (bandSum a0 a3 a5 a3p) =
      realBandSum (evalReal a0) (evalReal a3) (evalReal a5) (evalReal a3p) := by
  ext i j
  simp [realMatrix, bandSum, realBandSum, pZeroR, pLowR, pFiveR, pHighR]

/-- The exact finite kernel transported into the real embedding. -/
noncomputable def realDiffusionKernel (n : ℕ) :
    Matrix (Fin 12) (Fin 12) ℝ :=
  realMatrix (diffusionKernel n)

theorem realDiffusionKernel_exact (n : ℕ) :
    realDiffusionKernel n =
      realBandSum 0 (lowResponseR ^ (2 * n))
        (fiveResponseR ^ (2 * n)) (highResponseR ^ (2 * n)) := by
  rw [realDiffusionKernel, diffusionKernel_exact, realMatrix_bandSum]
  simp [lowResponseR, fiveResponseR, highResponseR]

/-- The quintet response relative to the unique slow band. -/
noncomputable def fiveRatio : ℝ := fiveResponseR / lowResponseR

/-- The Galois response relative to the unique slow band. -/
noncomputable def highRatio : ℝ := highResponseR / lowResponseR

/-- Both non-leading response ratios are strictly between zero and one, with
the Galois ratio strictly smaller than the quintet ratio. -/
theorem response_ratio_strict_order :
    0 < highRatio ∧ highRatio < fiveRatio ∧ fiveRatio < 1 := by
  rcases response_strict_order with ⟨hhigh0, hhighfive, hfivelow, hlow1⟩
  have hlow0 : 0 < lowResponseR :=
    hhigh0.trans (hhighfive.trans hfivelow)
  simp only [highRatio, fiveRatio]
  constructor
  · exact div_pos hhigh0 hlow0
  constructor
  · exact (div_lt_div_iff_of_pos_right hlow0).2 hhighfive
  · exact (div_lt_one hlow0).2 hfivelow

/-- Low-response normalization defined directly from the source kernel.  The
factor four is the unit-diagonal normalization of its eventual leading band. -/
noncomputable def normalizedKernel (n : ℕ) :
    Matrix (Fin 12) (Fin 12) ℝ :=
  (4 / lowResponseR ^ (2 * n)) • realDiffusionKernel n

/-- Exact expansion of the source-defined normalized kernel.  The target
projector occurs here as a proved spectral consequence, not in the definition
of `normalizedKernel`. -/
theorem normalizedKernel_exact (n : ℕ) :
    normalizedKernel n =
      4 • (pLowR + fiveRatio ^ (2 * n) • pFiveR +
        highRatio ^ (2 * n) • pHighR) := by
  rcases response_strict_order with ⟨hhigh0, hhighfive, hfivelow, hlow1⟩
  have hlow0 : 0 < lowResponseR :=
    hhigh0.trans (hhighfive.trans hfivelow)
  have hlowne : lowResponseR ≠ 0 := ne_of_gt hlow0
  rw [normalizedKernel, realDiffusionKernel_exact]
  ext i j
  simp only [realBandSum, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul]
  simp [fiveRatio, highRatio, div_pow]
  field_simp [pow_ne_zero _ hlowne]

/-- The two exact relative coefficients decay to zero. -/
theorem response_ratios_tendsto_zero :
    Tendsto (fun n : ℕ ↦ fiveRatio ^ (2 * n)) atTop (nhds 0) ∧
      Tendsto (fun n : ℕ ↦ highRatio ^ (2 * n)) atTop (nhds 0) := by
  rcases response_ratio_strict_order with ⟨hhigh0, hhighfive, hfive1⟩
  have hhigh1 : highRatio < 1 := hhighfive.trans hfive1
  have hfiveSq : 0 ≤ fiveRatio ^ 2 ∧ fiveRatio ^ 2 < 1 := by
    constructor
    · positivity
    · nlinarith
  have hhighSq : 0 ≤ highRatio ^ 2 ∧ highRatio ^ 2 < 1 := by
    constructor
    · positivity
    · nlinarith
  constructor
  · simpa [pow_mul] using
      (tendsto_pow_atTop_nhds_zero_of_lt_one hfiveSq.1 hfiveSq.2)
  · simpa [pow_mul] using
      (tendsto_pow_atTop_nhds_zero_of_lt_one hhighSq.1 hhighSq.2)

/-- Entrywise asymptotic selection of the unit-diagonal Gram matrix.  This is
an analytic consequence of the exact finite kernel, not a physical promotion
of the repair readback. -/
theorem normalizedKernel_tendsto_portGram (i j : Fin 12) :
    Tendsto (fun n : ℕ ↦ normalizedKernel n i j) atTop
      (nhds (realMatrix portGram i j)) := by
  rcases response_ratios_tendsto_zero with ⟨hfive, hhigh⟩
  have hlimit :
      Tendsto (fun n : ℕ ↦ normalizedKernel n i j) atTop
        (nhds ((4 : ℝ) * pLowR i j)) := by
    have hconst : Tendsto (fun _ : ℕ ↦ pLowR i j) atTop
        (nhds (pLowR i j)) := tendsto_const_nhds
    have h := (hconst.add
      ((hfive.mul_const (pFiveR i j)).add
        (hhigh.mul_const (pHighR i j)))).const_mul (4 : ℝ)
    convert h using 1
    all_goals
      simp [normalizedKernel_exact, Matrix.add_apply, Matrix.smul_apply,
        smul_eq_mul]
    all_goals ring_nf
  have hgram := congrArg realMatrix portGram_eq_four_pLow
  have hentry := congrFun (congrFun hgram i) j
  have hfour : evalReal (q5 4 0) = (4 : ℝ) := by
    norm_num [evalReal, q5]
  have htarget : realMatrix portGram i j = (4 : ℝ) * pLowR i j := by
    simpa [realMatrix, hfour] using hentry
  rw [htarget]
  exact hlimit

#print axioms projector_squares
#print axioms projector_cross_products
#print axioms projector_resolution
#print axioms repairMean_band_decomposition
#print axioms repairMean_pow
#print axioms diffusionKernel_exact
#print axioms portGram_eq_four_pLow
#print axioms portGram_antipodal_quotient
#print axioms pLowR_idempotent
#print axioms pLowR_symmetric
#print axioms intrinsicPortGenerator_antipode
#print axioms intrinsicPortGenerator_gram
#print axioms response_strict_order
#print axioms finite_response_coefficients_positive
#print axioms response_ratio_strict_order
#print axioms normalizedKernel_exact
#print axioms normalizedKernel_tendsto_portGram

end OPH.PortGramRepairCovariance
