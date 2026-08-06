import EventAlgebra.StateExpectation
import EventAlgebra.FiniteEffectClosureBoundary
import Tower.ConsensusTower

/-!
# The convex-affinity bridge for celestial binary responses

This issue-scoped module records the positive counterpart to the
`FiniteEffectClosureBoundary` obstruction.  That module proves that
continuity, range, and antipodal normalization alone do not force the
affine/Born form of a celestial binary response.  The bridge here derives
the missing affinity for one declared class of producers: preparation
families closed under convex mixing.

`ConvexPreparationFamily` packages a parameter type, a mixing operation
with parameter in `[0,1]`, a density-matrix assignment certified by
`EventAlgebra.IsState`, and the mixing law stating that the assigned
state of a mixture is the convex combination of the assigned states.
Convex mixtures of certified states are certified states
(`isState_convexMix`), and for such a family the outcome functional
`p ↦ Tr(ρ_p E)` on every fixed effect matrix is affine in the mixing
parameter (`ConvexPreparationFamily.outcome_affine`,
`ConvexPreparationFamily.outcome_re_affine`): the trace pairing is
linear in the state, so affinity is inherited from the declared convex
structure rather than assumed.  `segmentFamily` realizes the interface
on the segment between any two certified states, and
`towerSegmentFamily` instantiates it on the selected states of two
observers at one regulator of a `ConsensusTower`.

`CelestialEffectScheme` attaches one operator per spatial coordinate and
sends each spatial vector `v` to the operator `2⁻¹ • 1 + ∑ i, v i • C i`.
The induced response of any certified state is exactly the centered
affine ansatz `affineBinaryWeight` at the scheme's coefficient vector
(`CelestialEffectScheme.response_eq_affineBinaryWeight`) and obeys exact
antipodal normalization (`CelestialEffectScheme.response_antipodal`).
The composition theorem
`CelestialEffectScheme.dense_probability_tests_force_closed_unit_ball`
instantiates the affinity hypothesis of the attained dense-test theorem
`EventAlgebra.FiniteEffectClosureBoundary.dense_affine_probability_tests_force_closed_unit_ball`
with the derived affine form: probability bounds on any dense direction
family force the scheme's coefficient vector into the closed unit ball.
`affine_response_dense_tests_force_closed_unit_ball` is the same adaptor
for an abstract response with a supplied affine representation.

On the negative side the `(1 + n_z ^ 3) / 2` countermodel is restated
against this producer class.  No dimension, scheme, and certified state
induces the cube response
(`CelestialEffectScheme.cube_response_not_scheme_induced`), and the
disagreement is pinned on the four explicit directions of the attained
nonrepresentation proof: the three coordinate axes and `(3/5, 0, 4/5)`
(`cube_response_four_direction_disagreement`).

The witness layer instantiates everything on `M₂`: a two-parameter
convex family over the simplex `{(a, b) | 0 ≤ a, 0 ≤ b, a + b ≤ 1}`
mixing three projection vertices, two projection effects with exact
rational outcome forms, the halved-Pauli scheme whose operators at the
`z` and `x` axes are the two witness effects, the exact coefficient
vector `(1 - a - b, 0, a - b)`, and its closed-unit-ball receipt.

## Claim boundary

The bridge derives affinity for source families with a declared convex
mixing law.  Whether the OPH source produces such a family with
interlocking public contexts is the open B13 core; shared effects across
contexts, Gleason/Busch hypotheses, and a Born-rule closure are out of
scope here.  The schemes, families, and witness states are supplied by
hand; no physical instrument or measurement attachment is claimed.
-/

namespace OPH.QFT

open EventAlgebra
open EventAlgebra.FiniteEffectClosureBoundary
open OPH.C1Lorentz
open Matrix
open scoped ComplexOrder

noncomputable section

universe u

/-! ## Convex mixtures of certified states -/

/-- Convex mixtures of certified density matrices are certified density
matrices: positivity and unit trace pass through the mixture. -/
theorem isState_convexMix {d : ℕ} {ρ σ : Matrix (Fin d) (Fin d) ℂ}
    (hρ : IsState ρ) (hσ : IsState σ) {l : ℝ} (h0 : 0 ≤ l) (h1 : l ≤ 1) :
    IsState ((l : ℂ) • ρ + ((1 - l : ℝ) : ℂ) • σ) := by
  constructor
  · exact (hρ.1.smul (by exact_mod_cast h0)).add
      (hσ.1.smul (by
        have : (0 : ℝ) ≤ 1 - l := by linarith
        exact_mod_cast this))
  · rw [trace_add, trace_smul, trace_smul, hρ.2, hσ.2, smul_eq_mul,
      smul_eq_mul]
    push_cast
    ring

/-! ## Convex preparation families -/

/-- A preparation family with a declared convex mixing operation.  The
parameter type `P` labels preparations, `mix` combines two preparations
with a weight in `[0,1]`, `state` assigns a certified density matrix, and
`state_mix` is the mixing law: the state of a mixture is the convex
combination of the states.  The mixing law is declared structure data;
no source realization is claimed by the interface itself. -/
structure ConvexPreparationFamily (d : ℕ) (P : Type u) where
  /-- The convex mixing operation on preparation labels. -/
  mix : (l : ℝ) → l ∈ Set.Icc (0 : ℝ) 1 → P → P → P
  /-- The density-matrix assignment. -/
  state : P → Matrix (Fin d) (Fin d) ℂ
  /-- Every assigned matrix is a certified state. -/
  state_isState : ∀ p, IsState (state p)
  /-- The declared mixing law. -/
  state_mix : ∀ (l : ℝ) (hl : l ∈ Set.Icc (0 : ℝ) 1) (p q : P),
    state (mix l hl p q) = (l : ℂ) • state p + ((1 - l : ℝ) : ℂ) • state q

namespace ConvexPreparationFamily

variable {d : ℕ} {P : Type u}

/-- Affinity of the outcome functional in the mixing parameter, for every
fixed effect matrix `E`: the trace pairing is linear in the state, so the
mixing law transports convex combinations of preparations to convex
combinations of outcomes.  This is the positive affinity principle the
`FiniteEffectClosureBoundary` obstruction leaves open; here it holds by
construction for the declared producer class. -/
theorem outcome_affine (F : ConvexPreparationFamily d P)
    (E : Matrix (Fin d) (Fin d) ℂ) (l : ℝ) (hl : l ∈ Set.Icc (0 : ℝ) 1)
    (p q : P) :
    bornWeight (F.state (F.mix l hl p q)) E
      = (l : ℂ) * bornWeight (F.state p) E
        + ((1 - l : ℝ) : ℂ) * bornWeight (F.state q) E := by
  rw [F.state_mix l hl p q]
  simp only [bornWeight, Matrix.add_mul, Matrix.smul_mul, trace_add,
    trace_smul, smul_eq_mul]

/-- Real-part form of outcome affinity. -/
theorem outcome_re_affine (F : ConvexPreparationFamily d P)
    (E : Matrix (Fin d) (Fin d) ℂ) (l : ℝ) (hl : l ∈ Set.Icc (0 : ℝ) 1)
    (p q : P) :
    (bornWeight (F.state (F.mix l hl p q)) E).re
      = l * (bornWeight (F.state p) E).re
        + (1 - l) * (bornWeight (F.state q) E).re := by
  rw [F.outcome_affine E l hl p q]
  simp [Complex.add_re, Complex.mul_re]

/-- Outcomes of a convex preparation family on projection events are
genuine probabilities. -/
theorem outcome_mem_unitInterval (F : ConvexPreparationFamily d P)
    {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEvent E) (p : P) :
    (bornWeight (F.state p) E).re ∈ Set.Icc (0 : ℝ) 1 :=
  ⟨bornWeight_re_nonneg (F.state_isState p).1 hE,
    bornWeight_re_le_one (F.state_isState p) hE⟩

end ConvexPreparationFamily

/-- The segment family between two certified states: parameters are the
points of `[0,1]`, mixing is affine interpolation of parameters, and the
mixing law is the module identity collecting coefficients. -/
def segmentFamily {d : ℕ} (ρ σ : Matrix (Fin d) (Fin d) ℂ)
    (hρ : IsState ρ) (hσ : IsState σ) :
    ConvexPreparationFamily d {t : ℝ // t ∈ Set.Icc (0 : ℝ) 1} where
  mix l hl t s :=
    ⟨l * t.1 + (1 - l) * s.1, by
      have h1l : (0 : ℝ) ≤ 1 - l := by linarith [hl.2]
      constructor
      · exact add_nonneg (mul_nonneg hl.1 t.2.1) (mul_nonneg h1l s.2.1)
      · nlinarith [mul_le_mul_of_nonneg_left t.2.2 hl.1,
          mul_le_mul_of_nonneg_left s.2.2 h1l]⟩
  state t := (t.1 : ℂ) • ρ + ((1 - t.1 : ℝ) : ℂ) • σ
  state_isState t := isState_convexMix hρ hσ t.2.1 t.2.2
  state_mix l hl t s := by
    push_cast
    module

/-- The segment family on the selected states of two observers at one
regulator of a consensus tower.  The tower's `state_isState` receipts
certify both endpoints, so every mixture along the segment is a certified
state and every outcome functional on it is affine in the mixing
parameter. -/
def towerSegmentFamily {ι : Type u} [Preorder ι]
    (T : OPH.Tower.ConsensusTower ι) (r : ι) (o o' : T.Observer r) :
    ConvexPreparationFamily (T.dim r) {t : ℝ // t ∈ Set.Icc (0 : ℝ) 1} :=
  segmentFamily (T.state r o) (T.state r o')
    (T.state_isState r o) (T.state_isState r o')

/-! ## Celestial effect schemes -/

/-- A celestial effect scheme in matrix dimension `d`: one component
operator per spatial coordinate.  The scheme sends a spatial vector to a
centered operator pencil; no hermiticity, positivity, or effect-algebra
structure is required for the affinity computation. -/
structure CelestialEffectScheme (d : ℕ) where
  /-- The three component operators. -/
  component : Fin 3 → Matrix (Fin d) (Fin d) ℂ

namespace CelestialEffectScheme

variable {d : ℕ}

/-- The operator assigned to the spatial vector `v`. -/
def effect (S : CelestialEffectScheme d) (v : Spatial) :
    Matrix (Fin d) (Fin d) ℂ :=
  ((2⁻¹ : ℝ) : ℂ) • 1 + ∑ i, (v i : ℂ) • S.component i

/-- Exact antipodal normalization at the operator level. -/
theorem effect_antipodal (S : CelestialEffectScheme d) (v : Spatial) :
    S.effect v + S.effect (-v) = 1 := by
  simp only [effect, Pi.neg_apply, Complex.ofReal_neg, neg_smul,
    Finset.sum_neg_distrib]
  rw [add_add_add_comm, add_neg_cancel, add_zero, ← add_smul]
  push_cast
  norm_num

/-- The trace pairing against a scheme operator splits into the constant
and the coordinate-linear part. -/
theorem bornWeight_effect (S : CelestialEffectScheme d)
    (ρ : Matrix (Fin d) (Fin d) ℂ) (v : Spatial) :
    bornWeight ρ (S.effect v)
      = ((2⁻¹ : ℝ) : ℂ) * ρ.trace
        + ∑ i, (v i : ℂ) * bornWeight ρ (S.component i) := by
  simp only [effect, bornWeight, Matrix.mul_add, Matrix.trace_add,
    Matrix.mul_smul, Matrix.trace_smul, Matrix.mul_one, Finset.mul_sum,
    Matrix.trace_sum, smul_eq_mul]

/-- The real-valued response of a state through a scheme. -/
def response (S : CelestialEffectScheme d) (ρ : Matrix (Fin d) (Fin d) ℂ)
    (v : Spatial) : ℝ :=
  (bornWeight ρ (S.effect v)).re

/-- The coefficient vector of a state through a scheme: twice the real
part of the trace pairing against each component operator. -/
def blochCoefficient (S : CelestialEffectScheme d)
    (ρ : Matrix (Fin d) (Fin d) ℂ) : Spatial :=
  fun i => 2 * (bornWeight ρ (S.component i)).re

/-- The response of a certified state through a scheme is exactly the
centered affine ansatz of the `FiniteEffectClosureBoundary` module at the
scheme's coefficient vector.  Affinity of the celestial response is
derived from the operator pencil and the linearity of the trace pairing;
for this producer class it is a theorem rather than a hypothesis. -/
theorem response_eq_affineBinaryWeight (S : CelestialEffectScheme d)
    {ρ : Matrix (Fin d) (Fin d) ℂ} (hρ : IsState ρ) (v : Spatial) :
    S.response ρ v = affineBinaryWeight (S.blochCoefficient ρ) v := by
  simp only [response, affineBinaryWeight, spatialDot, blochCoefficient,
    bornWeight_effect, hρ.2, mul_one]
  rw [Fin.sum_univ_three, Fin.sum_univ_three]
  simp [Complex.add_re, Complex.mul_re]
  ring

/-- Exact antipodal normalization of the induced response. -/
theorem response_antipodal (S : CelestialEffectScheme d)
    {ρ : Matrix (Fin d) (Fin d) ℂ} (hρ : IsState ρ) (v : Spatial) :
    S.response ρ (-v) = 1 - S.response ρ v := by
  have hsum : bornWeight ρ (S.effect v) + bornWeight ρ (S.effect (-v)) = 1 := by
    rw [← bornWeight_add, S.effect_antipodal v, bornWeight_one hρ]
  have hre := congrArg Complex.re hsum
  simp only [Complex.add_re, Complex.one_re] at hre
  simp only [response]
  linarith

/-- Affinity of the response in the preparation, for a convex preparation
family: the response to every fixed direction is a convex combination of
the endpoint responses. -/
theorem response_mix_affine {P : Type u} (S : CelestialEffectScheme d)
    (F : ConvexPreparationFamily d P) (l : ℝ)
    (hl : l ∈ Set.Icc (0 : ℝ) 1) (p q : P) (v : Spatial) :
    S.response (F.state (F.mix l hl p q)) v
      = l * S.response (F.state p) v + (1 - l) * S.response (F.state q) v :=
  F.outcome_re_affine (S.effect v) l hl p q

/-- Composition with the attained dense-test theorem: the derived affine
form discharges the affinity hypothesis, so probability bounds on any
dense family of directions force the scheme's coefficient vector into the
closed unit ball.  The dense-test theorem is consumed as proved; no part
of it is restated. -/
theorem dense_probability_tests_force_closed_unit_ball
    (S : CelestialEffectScheme d) {ρ : Matrix (Fin d) (Fin d) ℂ}
    (hρ : IsState ρ) (D : Set CelestialSphere) (hDense : Dense D)
    (hProb : ∀ m ∈ D, S.response ρ m.1 ∈ Set.Icc (0 : ℝ) 1) :
    spatialNormSq (S.blochCoefficient ρ) ≤ 1 := by
  refine dense_affine_probability_tests_force_closed_unit_ball D hDense _ ?_
  intro m hm
  rw [← S.response_eq_affineBinaryWeight hρ m.1]
  exact hProb m hm

end CelestialEffectScheme

/-- The abstract adaptor: any celestial response with a supplied affine
representation and dense probability bounds has its coefficient vector in
the closed unit ball.  This is the exact interface shape of the attained
dense-test theorem, restated for responses rather than raw coefficient
vectors. -/
theorem affine_response_dense_tests_force_closed_unit_ball
    (F : CelestialSphere → ℝ) (q : Spatial)
    (hAffine : ∀ m : CelestialSphere, F m = affineBinaryWeight q m.1)
    (D : Set CelestialSphere) (hDense : Dense D)
    (hProb : ∀ m ∈ D, F m ∈ Set.Icc (0 : ℝ) 1) :
    spatialNormSq q ≤ 1 := by
  refine dense_affine_probability_tests_force_closed_unit_ball D hDense q ?_
  intro m hm
  rw [← hAffine m]
  exact hProb m hm

/-- End-to-end composition on tower data: dense probability tests on the
response of any segment mixture of two tower observers' selected states
force the scheme coefficient vector into the closed unit ball. -/
theorem towerSegment_dense_tests_force_closed_unit_ball {ι : Type u}
    [Preorder ι] (T : OPH.Tower.ConsensusTower ι) (r : ι)
    (o o' : T.Observer r) (S : CelestialEffectScheme (T.dim r))
    (t : {t : ℝ // t ∈ Set.Icc (0 : ℝ) 1}) (D : Set CelestialSphere)
    (hDense : Dense D)
    (hProb : ∀ m ∈ D,
      S.response ((towerSegmentFamily T r o o').state t) m.1
        ∈ Set.Icc (0 : ℝ) 1) :
    spatialNormSq
      (S.blochCoefficient ((towerSegmentFamily T r o o').state t)) ≤ 1 :=
  S.dense_probability_tests_force_closed_unit_ball
    ((towerSegmentFamily T r o o').state_isState t) D hDense hProb

/-! ## The cube countermodel against the producer class -/

/-- The positive `x` axis. -/
def axisX : Spatial := ![1, 0, 0]

/-- The positive `y` axis. -/
def axisY : Spatial := ![0, 1, 0]

/-- The positive `z` axis. -/
def axisZ : Spatial := ![0, 0, 1]

/-- The probe direction `(3/5, 0, 4/5)` separating `z ^ 3` from `z`. -/
def probeDirection : Spatial := ![3 / 5, 0, 4 / 5]

theorem axisX_unit : spatialNormSq axisX = 1 := by
  norm_num [axisX, spatialNormSq, Fin.sum_univ_succ]

theorem axisY_unit : spatialNormSq axisY = 1 := by
  norm_num [axisY, spatialNormSq, Fin.sum_univ_succ]

theorem axisZ_unit : spatialNormSq axisZ = 1 := by
  norm_num [axisZ, spatialNormSq, Fin.sum_univ_succ]

theorem probeDirection_unit : spatialNormSq probeDirection = 1 := by
  norm_num [probeDirection, spatialNormSq, Fin.sum_univ_succ]

/-- Finite-support sharpening of the nonrepresentation theorem: no
coefficient vector matches the cube response on the four explicit
directions alone.  The three axes pin the coefficient vector to
`(0, 0, 1)` and the probe direction separates `z ^ 3` from `z`. -/
theorem affineBinaryWeight_four_direction_disagreement (q : Spatial) :
    ¬ (nonlinearBinaryWeight axisX = affineBinaryWeight q axisX ∧
        nonlinearBinaryWeight axisY = affineBinaryWeight q axisY ∧
        nonlinearBinaryWeight axisZ = affineBinaryWeight q axisZ ∧
        nonlinearBinaryWeight probeDirection
          = affineBinaryWeight q probeDirection) := by
  rintro ⟨hx, hy, hz, hw⟩
  have hq0 : q 0 = 0 := by
    simp [axisX, nonlinearBinaryWeight, affineBinaryWeight, spatialDot,
      Fin.sum_univ_succ, Matrix.cons_val_two] at hx
    norm_num at hx ⊢
    exact hx
  have hq1 : q 1 = 0 := by
    simp [axisY, nonlinearBinaryWeight, affineBinaryWeight, spatialDot,
      Fin.sum_univ_succ, Matrix.cons_val_two] at hy
    norm_num at hy ⊢
    exact hy
  have hq2 : q 2 = 1 := by
    simp [axisZ, nonlinearBinaryWeight, affineBinaryWeight, spatialDot,
      Fin.sum_univ_succ, Matrix.cons_val_two] at hz
    linarith
  simp [probeDirection, nonlinearBinaryWeight, affineBinaryWeight,
    spatialDot, Fin.sum_univ_succ, Matrix.cons_val_two] at hw
  simp [hq0, hq2] at hw
  norm_num at hw

/-- The cube response disagrees with the response of every certified
state through every scheme in every matrix dimension on the four explicit
directions. -/
theorem cube_response_four_direction_disagreement
    {d : ℕ} (S : CelestialEffectScheme d)
    {ρ : Matrix (Fin d) (Fin d) ℂ} (hρ : IsState ρ) :
    ¬ (nonlinearBinaryWeight axisX = S.response ρ axisX ∧
        nonlinearBinaryWeight axisY = S.response ρ axisY ∧
        nonlinearBinaryWeight axisZ = S.response ρ axisZ ∧
        nonlinearBinaryWeight probeDirection = S.response ρ probeDirection) := by
  rintro ⟨hx, hy, hz, hw⟩
  exact affineBinaryWeight_four_direction_disagreement (S.blochCoefficient ρ)
    ⟨hx.trans (S.response_eq_affineBinaryWeight hρ axisX),
      hy.trans (S.response_eq_affineBinaryWeight hρ axisY),
      hz.trans (S.response_eq_affineBinaryWeight hρ axisZ),
      hw.trans (S.response_eq_affineBinaryWeight hρ probeDirection)⟩

/-- The countermodel boundary restated for the producer class: the cube
response `(1 + n_z ^ 3) / 2` is induced by no certified state through any
scheme in any finite matrix dimension. -/
theorem cube_response_not_scheme_induced :
    ¬ ∃ (d : ℕ) (S : CelestialEffectScheme d)
        (ρ : Matrix (Fin d) (Fin d) ℂ),
        IsState ρ ∧ ∀ v : Spatial, spatialNormSq v = 1 →
          nonlinearBinaryWeight v = S.response ρ v := by
  rintro ⟨d, S, ρ, hρ, hEq⟩
  exact cube_response_four_direction_disagreement S hρ
    ⟨hEq axisX axisX_unit, hEq axisY axisY_unit, hEq axisZ axisZ_unit,
      hEq probeDirection probeDirection_unit⟩

/-! ## The `M₂` witness family -/

/-- The projection onto the first basis vector. -/
def vertexZero : Matrix (Fin 2) (Fin 2) ℂ := !![1, 0; 0, 0]

/-- The projection onto the second basis vector. -/
def vertexOne : Matrix (Fin 2) (Fin 2) ℂ := !![0, 0; 0, 1]

/-- The projection onto the equal-superposition vector. -/
def vertexPlus : Matrix (Fin 2) (Fin 2) ℂ :=
  !![((2⁻¹ : ℝ) : ℂ), ((2⁻¹ : ℝ) : ℂ); ((2⁻¹ : ℝ) : ℂ), ((2⁻¹ : ℝ) : ℂ)]

theorem isEvent_vertexZero : IsEvent vertexZero := by
  constructor
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [vertexZero, Matrix.conjTranspose_apply]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [vertexZero, Matrix.mul_apply, Fin.sum_univ_two]

theorem isEvent_vertexOne : IsEvent vertexOne := by
  constructor
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [vertexOne, Matrix.conjTranspose_apply]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [vertexOne, Matrix.mul_apply, Fin.sum_univ_two]

theorem isEvent_vertexPlus : IsEvent vertexPlus := by
  constructor
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [vertexPlus, Matrix.conjTranspose_apply]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      (simp [vertexPlus, Matrix.mul_apply, Fin.sum_univ_two]; norm_num)

theorem trace_vertexZero : vertexZero.trace = 1 := by
  simp [vertexZero, Matrix.trace_fin_two]

theorem trace_vertexOne : vertexOne.trace = 1 := by
  simp [vertexOne, Matrix.trace_fin_two]

theorem trace_vertexPlus : vertexPlus.trace = 1 := by
  simp only [vertexPlus, Matrix.trace_fin_two]
  push_cast
  norm_num

/-- Parameters of the witness family: the closed planar simplex. -/
def WitnessParam : Type :=
  {ab : ℝ × ℝ // 0 ≤ ab.1 ∧ 0 ≤ ab.2 ∧ ab.1 + ab.2 ≤ 1}

/-- The witness state at simplex parameter `(a, b)`: the convex
combination `a • |0⟩⟨0| + b • |1⟩⟨1| + (1 - a - b) • |+⟩⟨+|`. -/
def affinityWitnessState (p : WitnessParam) : Matrix (Fin 2) (Fin 2) ℂ :=
  (p.1.1 : ℂ) • vertexZero + (p.1.2 : ℂ) • vertexOne
    + ((1 - p.1.1 - p.1.2 : ℝ) : ℂ) • vertexPlus

/-- The witness state written entrywise. -/
theorem affinityWitnessState_eq (p : WitnessParam) :
    affinityWitnessState p =
      !![(((1 + p.1.1 - p.1.2) / 2 : ℝ) : ℂ),
          (((1 - p.1.1 - p.1.2) / 2 : ℝ) : ℂ);
        (((1 - p.1.1 - p.1.2) / 2 : ℝ) : ℂ),
          (((1 - p.1.1 + p.1.2) / 2 : ℝ) : ℂ)] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    (simp [affinityWitnessState, vertexZero, vertexOne, vertexPlus]; ring)

/-- Every witness parameter yields a certified state. -/
theorem affinityWitnessState_isState (p : WitnessParam) : IsState (affinityWitnessState p) := by
  obtain ⟨⟨a, b⟩, ha, hb, hab⟩ := p
  constructor
  · refine Matrix.PosSemidef.add (Matrix.PosSemidef.add ?_ ?_) ?_
    · exact isEvent_vertexZero.posSemidef.smul (by exact_mod_cast ha)
    · exact isEvent_vertexOne.posSemidef.smul (by exact_mod_cast hb)
    · refine isEvent_vertexPlus.posSemidef.smul ?_
      have : (0 : ℝ) ≤ 1 - a - b := by linarith
      exact_mod_cast this
  · simp only [affinityWitnessState, trace_add, trace_smul, trace_vertexZero,
      trace_vertexOne, trace_vertexPlus, smul_eq_mul]
    push_cast
    ring

/-- The witness family: the two-parameter simplex family on `M₂` with
componentwise affine mixing. -/
def witnessFamily : ConvexPreparationFamily 2 WitnessParam where
  mix l hl p q :=
    ⟨(l * p.1.1 + (1 - l) * q.1.1, l * p.1.2 + (1 - l) * q.1.2), by
      have h1l : (0 : ℝ) ≤ 1 - l := by linarith [hl.2]
      refine ⟨add_nonneg (mul_nonneg hl.1 p.2.1) (mul_nonneg h1l q.2.1),
        add_nonneg (mul_nonneg hl.1 p.2.2.1) (mul_nonneg h1l q.2.2.1), ?_⟩
      nlinarith [mul_le_mul_of_nonneg_left p.2.2.2 hl.1,
        mul_le_mul_of_nonneg_left q.2.2.2 h1l]⟩
  state := affinityWitnessState
  state_isState := affinityWitnessState_isState
  state_mix l hl p q := by
    simp only [affinityWitnessState]
    push_cast
    module

/-- Exact outcome of the witness family on the first effect. -/
theorem witness_outcome_vertexZero (p : WitnessParam) :
    bornWeight (affinityWitnessState p) vertexZero
      = (((1 + p.1.1 - p.1.2) / 2 : ℝ) : ℂ) := by
  rw [affinityWitnessState_eq p]
  simp [bornWeight, vertexZero, Matrix.trace_fin_two]

/-- Exact outcome of the witness family on the second effect. -/
theorem witness_outcome_vertexPlus (p : WitnessParam) :
    bornWeight (affinityWitnessState p) vertexPlus
      = (((2 - p.1.1 - p.1.2) / 2 : ℝ) : ℂ) := by
  rw [affinityWitnessState_eq p]
  simp [bornWeight, vertexPlus, Matrix.trace_fin_two]
  ring

/-- Affinity of the first witness outcome in the mixing parameter,
written with the exact closed forms on both sides. -/
theorem witness_outcome_vertexZero_affine (l : ℝ)
    (hl : l ∈ Set.Icc (0 : ℝ) 1) (p q : WitnessParam) :
    bornWeight (witnessFamily.state (witnessFamily.mix l hl p q)) vertexZero
      = (l : ℂ) * (((1 + p.1.1 - p.1.2) / 2 : ℝ) : ℂ)
        + ((1 - l : ℝ) : ℂ) * (((1 + q.1.1 - q.1.2) / 2 : ℝ) : ℂ) := by
  rw [witnessFamily.outcome_affine vertexZero l hl p q]
  simp only [witnessFamily, witness_outcome_vertexZero]

/-- Affinity of the second witness outcome in the mixing parameter,
written with the exact closed forms on both sides. -/
theorem witness_outcome_vertexPlus_affine (l : ℝ)
    (hl : l ∈ Set.Icc (0 : ℝ) 1) (p q : WitnessParam) :
    bornWeight (witnessFamily.state (witnessFamily.mix l hl p q)) vertexPlus
      = (l : ℂ) * (((2 - p.1.1 - p.1.2) / 2 : ℝ) : ℂ)
        + ((1 - l : ℝ) : ℂ) * (((2 - q.1.1 - q.1.2) / 2 : ℝ) : ℂ) := by
  rw [witnessFamily.outcome_affine vertexPlus l hl p q]
  simp only [witnessFamily, witness_outcome_vertexPlus]

/-- The halved-Pauli scheme on `M₂`. -/
def pauliScheme : CelestialEffectScheme 2 where
  component :=
    ![!![0, ((2⁻¹ : ℝ) : ℂ); ((2⁻¹ : ℝ) : ℂ), 0],
      !![0, -(((2⁻¹ : ℝ) : ℂ) * Complex.I); ((2⁻¹ : ℝ) : ℂ) * Complex.I, 0],
      !![((2⁻¹ : ℝ) : ℂ), 0; 0, -((2⁻¹ : ℝ) : ℂ)]]

/-- The Pauli-scheme operator at the `z` axis is the first witness
effect. -/
theorem pauliScheme_effect_axisZ : pauliScheme.effect axisZ = vertexZero := by
  ext i j
  fin_cases i <;> fin_cases j
  · simp [CelestialEffectScheme.effect, pauliScheme, axisZ, vertexZero,
      Fin.sum_univ_three]
    norm_num
  · simp [CelestialEffectScheme.effect, pauliScheme, axisZ, vertexZero,
      Fin.sum_univ_three]
  · simp [CelestialEffectScheme.effect, pauliScheme, axisZ, vertexZero,
      Fin.sum_univ_three]
  · simp [CelestialEffectScheme.effect, pauliScheme, axisZ, vertexZero,
      Fin.sum_univ_three]

/-- The Pauli-scheme operator at the `x` axis is the second witness
effect. -/
theorem pauliScheme_effect_axisX : pauliScheme.effect axisX = vertexPlus := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [CelestialEffectScheme.effect, pauliScheme, axisX, vertexPlus,
      Fin.sum_univ_three]

/-- Exact coefficient vector of the witness state through the Pauli
scheme. -/
theorem witness_blochCoefficient (p : WitnessParam) :
    pauliScheme.blochCoefficient (affinityWitnessState p)
      = ![1 - p.1.1 - p.1.2, 0, p.1.1 - p.1.2] := by
  rw [show pauliScheme.blochCoefficient (affinityWitnessState p)
      = fun i => 2 * (bornWeight (affinityWitnessState p)
          (pauliScheme.component i)).re from rfl]
  funext i
  rw [affinityWitnessState_eq p]
  fin_cases i
  · simp [pauliScheme, bornWeight, Matrix.trace_fin_two]
    ring
  · simp [pauliScheme, bornWeight, Matrix.trace_fin_two]
  · simp [pauliScheme, bornWeight, Matrix.trace_fin_two]
    ring

/-- The witness coefficient vectors fill out receipts of the composition
theorem: they lie in the closed unit ball for every simplex parameter. -/
theorem witness_blochCoefficient_mem_closed_unit_ball (p : WitnessParam) :
    spatialNormSq (pauliScheme.blochCoefficient (affinityWitnessState p)) ≤ 1 := by
  rw [witness_blochCoefficient p]
  obtain ⟨⟨a, b⟩, ha, hb, hab⟩ := p
  simp only [spatialNormSq, Fin.sum_univ_three, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two,
    Matrix.tail_cons]
  nlinarith [mul_nonneg ha (by linarith : (0 : ℝ) ≤ 1 - a),
    mul_nonneg hb (by linarith : (0 : ℝ) ≤ 1 - b)]

end

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms isState_convexMix
#print axioms ConvexPreparationFamily.outcome_affine
#print axioms ConvexPreparationFamily.outcome_re_affine
#print axioms ConvexPreparationFamily.outcome_mem_unitInterval
#print axioms segmentFamily
#print axioms towerSegmentFamily
#print axioms CelestialEffectScheme.effect_antipodal
#print axioms CelestialEffectScheme.bornWeight_effect
#print axioms CelestialEffectScheme.response_eq_affineBinaryWeight
#print axioms CelestialEffectScheme.response_antipodal
#print axioms CelestialEffectScheme.response_mix_affine
#print axioms CelestialEffectScheme.dense_probability_tests_force_closed_unit_ball
#print axioms affine_response_dense_tests_force_closed_unit_ball
#print axioms towerSegment_dense_tests_force_closed_unit_ball
#print axioms affineBinaryWeight_four_direction_disagreement
#print axioms cube_response_four_direction_disagreement
#print axioms cube_response_not_scheme_induced
#print axioms isEvent_vertexZero
#print axioms isEvent_vertexOne
#print axioms isEvent_vertexPlus
#print axioms affinityWitnessState_isState
#print axioms witnessFamily
#print axioms witness_outcome_vertexZero
#print axioms witness_outcome_vertexPlus
#print axioms witness_outcome_vertexZero_affine
#print axioms witness_outcome_vertexPlus_affine
#print axioms pauliScheme_effect_axisZ
#print axioms pauliScheme_effect_axisX
#print axioms witness_blochCoefficient
#print axioms witness_blochCoefficient_mem_closed_unit_ball

end OPH.QFT
