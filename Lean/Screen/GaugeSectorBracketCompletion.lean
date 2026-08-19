import AssembledActionComposition
import OrientedFaceInvariantMetric

/-!
# Gauge-sector completion of the assembled action at the committed bracket level (V3, issue #735, OL-G9)

The committed assembled structural action of
`AssembledActionComposition` carries its gauge kinetic sector as a
declared diagonal pairing: one declared strength magnitude per
invariant carrier sector (`3+`, `3-`, `5`), squared against its
committed weight.  This module upgrades that sector to the committed
invariant-form level as an extension.  It imports the committed
modules and edits none of them.

The carrier.  The pinned bracket carrier splits multiplicity-free
into the sectors `1 + 3 + 3' + 5` (the commutant-dimension receipt of
the invariant-metric certificate
`code/b14_jacobi/invariant_metric_phase.certificate.json`; the
sector-balanced point of that phase diagram selects the committed
compact family `G`, re-exported below from
`OrientedFaceInvariantMetric`).  The committed weight space
`ThreeSectorWeights` of `GaugeKineticInvariantForms` carries one
weight per nontrivial sector; the trivial sector carries no committed
weight and the committed distance forms never involve its scale.  The
typed gauge-field slot of this module is therefore one coordinate
vector per weighted sector: `Fin 3 → ℝ`, `Fin 3 → ℝ`, `Fin 5 → ℝ`.

The kinetic functional.  For a committed weight vector `w`, the
quadratic functional is the sector-weighted sum of coordinate squared
norms, with its symmetric bilinear polarization.  Proved here:

1. **Exact identification with the committed diagonal pairing.**  On
   the axis slice (one basis coordinate per sector, the committed
   basis of the diagonal pairing) the functional is exactly the
   committed `gaugeKineticEnergy`; conversely the functional factors
   exactly through the committed pairing at the sector magnitudes.
   The committed assembled action embeds in the completed one
   (`completed_embeds_committed`), and the completed action factors
   exactly through the committed action
   (`completed_factors_through_committed`).  Agreement on the axis
   slice pins the weights (`axis_restriction_pins`).
2. **Invariance at the committed level, exact grade: infinitesimal.**
   The invariance carried here is the polarized first-order identity
   `B_w(L x, y) + B_w(x, L y) = 0` for typed linear generators `L` on
   the carrier.  Sector-internal skew-adjoint generators leave every
   weight vector invariant (the multiplicity-free commutant grade;
   this clause consumes no committed constraint row).  The committed
   constraint rows are realized exactly as cross-sector invariance
   rows: for the exhibited `3-`/`5` cross generator at coefficient
   ratio `(1, -sqrt 5)`, infinitesimal invariance holds exactly when
   the committed `G` row `GAdInvariant` holds
   (`gMirrorGen_invariance_iff`); the mirror `3+`/`5` generator
   realizes the committed `F` row (`fMirrorGen_invariance_iff`); and
   simultaneous infinitesimal invariance under both exhibited
   generators leaves exactly the committed one-ray mirror family
   (`mirror_pair_invariance_one_ray`).  No group action, no bracket
   representation, and no derivation of the generators from the
   pinned face bracket is claimed in Lean: the sector pattern and the
   coefficient ratio of the cross generators are declared here so
   that their invariance rows are exactly the committed constraint
   rows, and the derivation of those rows from the pinned bracket
   tensors is the independently replayed certificate computation of
   `code/e9_kinetic/gauge_kinetic_invariant_forms.certificate.json`.
3. **The completed assembled action.**  One real function of one
   completed field tuple: the invariant-form kinetic functional at
   the committed two-parameter `G` normal form, plus the committed
   quartic potential, gauge mass form, and realified Yukawa line, on
   the unchanged committed premise bundle
   `AssembledActionPremiseData`.  Exact restriction receipts recover
   the committed action under the axis identification and recover
   each committed sector functional on its slice; the global abelian
   phase invariance and the committed-charge invariance persist; the
   parameter census is unchanged: the invariant-form coefficients are
   exactly the two declared normal-form parameters of the committed
   family (`completed_gauge_parameter_count` via the committed
   `g_exact_two_parameter`), and the completed action introduces no
   coefficient beyond the committed census
   (`completed_parameter_census`).  The committed inhabitant extends,
   with one exact off-axis witness evaluation at a carrier point
   outside the axis slice.
4. **Load-bearing receipts.**  The control weight vector `(1, 1, 1)`
   violates the committed `G` row and breaks the cross-generator
   invariance at an exhibited carrier pair with exact value
   `1 - sqrt 5`; the control generator at ratio `(1, 1)` breaks the
   invariance of the committed normal form; the axis restriction pins
   the extension against silently changing the committed pairing.

Consumed register rows.  The bundle field `ew` carries the committed
rows of the assembled action: PR-09, PR-10, PR-11, PR-59 through its
base, declared parameter rows PR-48, PR-49, PR-50, PR-54.  The two
gauge weight parameters stay declared under the classification row
PR-11, matching inventory rows SML-01 and SML-02.  The carrier sector
split and the balanced-point selection anchor are committed content
of rows PR-09, PR-10, PR-11 through `OrientedFaceInvariantMetric`.

Boundary and nonclaims.  The invariant-form coefficients are the
declared committed normal-form parameters and select nothing.  There
is no source action: register row PR-54 stays open, and no
derivative, spacetime carrier, or field-strength two-form appears;
the carrier coordinates are algebraic field slots exactly as the
committed strength magnitudes are.  There is no Spin or chirality
attachment (row PR-47 open), no global form (rows PR-35, PR-46), no
three-family structure (row PR-36), no measured coupling, and no
quantization.  The invariance grade is infinitesimal and typed; no
gauge group action on the carrier is committed or claimed.  OL-G9
stays partial, with the gauge kinetic sector entry upgraded from the
declared diagonal pairing to the committed invariant-form level.

Falsifier.  The completion fails if the axis evaluation differs from
the committed diagonal pairing, if the sector-magnitude factorization
fails, if the cross-generator invariance row differs from the
committed constraint row on either mirror side, if the exhibited
control weight vector fails to break invariance at the exhibited
carrier pair, if any restriction receipt differs from the committed
functional, or if the census equality fails on agreeing parameter
lists.
-/

namespace OPH.GaugeSectorBracketCompletion

open OPH.GaugeKineticInvariantForms (ThreeSectorWeights FAdInvariant
  GAdInvariant gNormalForm mirrorCommonRay
  mirror_common_extra_premise_one_ray g_exact_two_parameter)
open OPH.AssembledActionComposition (SectorStrength FieldTuple
  AssembledActionPremiseData assembledAction baseTuple scaleTuple
  gaugeKineticEnergy gaugeKineticEnergy_zero assembledAction_base
  assembled_abelian_invariance committed_charge_abelian_invariance
  committedAssembledActionPremiseData committed_inhabitant_value
  assembled_gauge_coefficients_G_invariant
  assembled_gauge_energy_nonneg gaugeMassForm_zero_direction
  yukawaLine_zero_doublet)
open OPH.ElectroweakBreakingComposition (Doublet normSqSum vacuum
  potential potential_vacuum gaugeMassForm gaugeMassForm_exact
  EWParam yukawaLine epsPair wDirection scalarChargeSix
  committedEWBreakingPremiseData)

/-! ## The typed carrier: one coordinate vector per weighted sector -/

/-- The typed gauge-field slot on the committed bracket carrier: one
real coordinate vector per weighted invariant sector of the committed
multiplicity-free split `1 + 3 + 3' + 5`.  The trivial sector carries
no committed weight and no slot. -/
@[ext]
structure CarrierField where
  /-- Coordinates on the `3+` sector. -/
  threePlus : Fin 3 → ℝ
  /-- Coordinates on the `3-` sector. -/
  threeMinus : Fin 3 → ℝ
  /-- Coordinates on the `5` sector. -/
  five : Fin 5 → ℝ

/-- The coordinate pairing of one sector. -/
def sectorPair {n : ℕ} (x y : Fin n → ℝ) : ℝ := ∑ i, x i * y i

/-- The coordinate squared norm of one sector. -/
def sectorNormSq {n : ℕ} (x : Fin n → ℝ) : ℝ := sectorPair x x

theorem sectorNormSq_nonneg {n : ℕ} (x : Fin n → ℝ) :
    0 ≤ sectorNormSq x := by
  unfold sectorNormSq sectorPair
  exact Finset.sum_nonneg fun i _ => mul_self_nonneg (x i)

/-- The symmetric bilinear polarization of the invariant-form family:
one committed weight per sector against the sector coordinate
pairing. -/
def invariantKineticBilinear (w : ThreeSectorWeights)
    (x y : CarrierField) : ℝ :=
  w.threePlus * sectorPair x.threePlus y.threePlus
    + w.threeMinus * sectorPair x.threeMinus y.threeMinus
    + w.five * sectorPair x.five y.five

/-- **The quadratic kinetic functional at the committed invariant-form
level**: the sector-weighted sum of coordinate squared norms, the
diagonal of the bilinear polarization. -/
def invariantKineticForm (w : ThreeSectorWeights) (x : CarrierField) :
    ℝ :=
  invariantKineticBilinear w x x

/-! ## Exact identification with the committed diagonal pairing -/

/-- The axis slice: the committed basis of the diagonal pairing, one
strength magnitude on the first coordinate of each sector. -/
def axisCarrier (s : SectorStrength) : CarrierField where
  threePlus := ![s.threePlus, 0, 0]
  threeMinus := ![s.threeMinus, 0, 0]
  five := ![s.five, 0, 0, 0, 0]

/-- **The committed diagonal pairing is the axis restriction of the
invariant-form functional.**  An exact identification under the
committed basis, for every weight vector. -/
theorem invariantKineticForm_axis (w : ThreeSectorWeights)
    (s : SectorStrength) :
    invariantKineticForm w (axisCarrier s) = gaugeKineticEnergy w s := by
  simp only [invariantKineticForm, invariantKineticBilinear,
    sectorPair, axisCarrier, gaugeKineticEnergy, Fin.sum_univ_three,
    Fin.sum_univ_five, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.cons_val_three, Matrix.cons_val_four,
    Matrix.head_cons, Matrix.tail_cons]
  ring

/-- The sector magnitudes of a carrier field: the committed strength
tuple through which the invariant-form functional factors. -/
noncomputable def sectorMagnitudes (x : CarrierField) :
    SectorStrength :=
  ⟨Real.sqrt (sectorNormSq x.threePlus),
    Real.sqrt (sectorNormSq x.threeMinus),
    Real.sqrt (sectorNormSq x.five)⟩

/-- **The invariant-form functional factors exactly through the
committed diagonal pairing** at the sector magnitudes. -/
theorem gaugeKineticEnergy_sectorMagnitudes (w : ThreeSectorWeights)
    (x : CarrierField) :
    gaugeKineticEnergy w (sectorMagnitudes x)
      = invariantKineticForm w x := by
  simp only [gaugeKineticEnergy, sectorMagnitudes,
    invariantKineticForm, invariantKineticBilinear]
  rw [Real.sq_sqrt (sectorNormSq_nonneg x.threePlus),
    Real.sq_sqrt (sectorNormSq_nonneg x.threeMinus),
    Real.sq_sqrt (sectorNormSq_nonneg x.five)]
  rfl

/-- **The axis restriction pins the weights.**  Two weight vectors
whose functionals agree on the axis slice are equal: the extension
cannot silently change the committed diagonal pairing. -/
theorem axis_restriction_pins (w w' : ThreeSectorWeights)
    (h : ∀ s : SectorStrength,
      invariantKineticForm w (axisCarrier s)
        = invariantKineticForm w' (axisCarrier s)) :
    w = w' := by
  have h1 := h ⟨1, 0, 0⟩
  have h2 := h ⟨0, 1, 0⟩
  have h3 := h ⟨0, 0, 1⟩
  rw [invariantKineticForm_axis, invariantKineticForm_axis] at h1 h2 h3
  simp only [gaugeKineticEnergy] at h1 h2 h3
  norm_num at h1 h2 h3
  ext
  · exact h1
  · exact h2
  · exact h3

/-! ## Invariance at the committed level: exact grade, infinitesimal -/

/-- Infinitesimal invariance of the weighted form under a typed
carrier generator: the polarized first-order identity.  This is the
exact grade the committed structure supports; no group action is
claimed. -/
def InfinitesimallyInvariant (w : ThreeSectorWeights)
    (L : CarrierField → CarrierField) : Prop :=
  ∀ x y : CarrierField,
    invariantKineticBilinear w (L x) y
      + invariantKineticBilinear w x (L y) = 0

/-- Sector-internal skew-adjointness of a carrier generator: the
polarized pairing identity within each weighted sector. -/
def SectorInternalSkew (L : CarrierField → CarrierField) : Prop :=
  (∀ x y : CarrierField,
      sectorPair (L x).threePlus y.threePlus
        + sectorPair x.threePlus (L y).threePlus = 0)
    ∧ (∀ x y : CarrierField,
      sectorPair (L x).threeMinus y.threeMinus
        + sectorPair x.threeMinus (L y).threeMinus = 0)
    ∧ (∀ x y : CarrierField,
      sectorPair (L x).five y.five
        + sectorPair x.five (L y).five = 0)

/-- Every sector-internal skew-adjoint generator leaves every weight
vector infinitesimally invariant: the multiplicity-free commutant
grade.  This clause consumes no committed constraint row; the
committed rows are exactly the cross-sector conditions below. -/
theorem sector_internal_skew_invariant (w : ThreeSectorWeights)
    (L : CarrierField → CarrierField) (hL : SectorInternalSkew L) :
    InfinitesimallyInvariant w L := by
  intro x y
  have h1 := hL.1 x y
  have h2 := hL.2.1 x y
  have h3 := hL.2.2 x y
  simp only [invariantKineticBilinear]
  linear_combination w.threePlus * h1 + w.threeMinus * h2
    + w.five * h3

/-- One concrete nonzero sector-internal generator: an infinitesimal
rotation in the first two coordinates of the `3+` sector. -/
def rotPlusGen (x : CarrierField) : CarrierField where
  threePlus := ![-(x.threePlus 1), x.threePlus 0, 0]
  threeMinus := fun _ => 0
  five := fun _ => 0

theorem rotPlusGen_sectorInternalSkew : SectorInternalSkew rotPlusGen := by
  refine ⟨fun x y => ?_, fun x y => ?_, fun x y => ?_⟩ <;>
    simp only [sectorPair, rotPlusGen, Fin.sum_univ_three,
      Fin.sum_univ_five, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons] <;>
    ring

/-- The `3-`/`5` cross generator at declared coefficients: the sector
pattern of the committed `G` constraint row.  The `3-` slot carries
the leading `5` coordinates at coefficient `μ`; the `5` slot carries
the `3-` coordinates at coefficient `ν` on its leading block. -/
def crossMinusFiveGen (μ ν : ℝ) (x : CarrierField) : CarrierField where
  threePlus := fun _ => 0
  threeMinus := ![μ * x.five 0, μ * x.five 1, μ * x.five 2]
  five := ![ν * x.threeMinus 0, ν * x.threeMinus 1,
    ν * x.threeMinus 2, 0, 0]

/-- The `3+`/`5` mirror cross generator: the sector pattern of the
committed `F` constraint row. -/
def crossPlusFiveGen (μ ν : ℝ) (x : CarrierField) : CarrierField where
  threePlus := ![μ * x.five 0, μ * x.five 1, μ * x.five 2]
  threeMinus := fun _ => 0
  five := ![ν * x.threePlus 0, ν * x.threePlus 1,
    ν * x.threePlus 2, 0, 0]

/-- The cross pairing of the `3-`/`5` sectors on the shared leading
block. -/
def crossPairingMinus (x y : CarrierField) : ℝ :=
  (x.five 0 * y.threeMinus 0 + x.five 1 * y.threeMinus 1
      + x.five 2 * y.threeMinus 2)
    + (x.threeMinus 0 * y.five 0 + x.threeMinus 1 * y.five 1
      + x.threeMinus 2 * y.five 2)

/-- The cross pairing of the `3+`/`5` sectors on the shared leading
block. -/
def crossPairingPlus (x y : CarrierField) : ℝ :=
  (x.five 0 * y.threePlus 0 + x.five 1 * y.threePlus 1
      + x.five 2 * y.threePlus 2)
    + (x.threePlus 0 * y.five 0 + x.threePlus 1 * y.five 1
      + x.threePlus 2 * y.five 2)

/-- Exact polarized value of the `3-`/`5` cross generator: one linear
row in the weights times the cross pairing. -/
theorem crossMinusFive_bilinear (μ ν : ℝ) (w : ThreeSectorWeights)
    (x y : CarrierField) :
    invariantKineticBilinear w (crossMinusFiveGen μ ν x) y
        + invariantKineticBilinear w x (crossMinusFiveGen μ ν y)
      = (w.threeMinus * μ + w.five * ν) * crossPairingMinus x y := by
  simp only [invariantKineticBilinear, sectorPair, crossMinusFiveGen,
    crossPairingMinus, Fin.sum_univ_three, Fin.sum_univ_five,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.cons_val_four, Matrix.head_cons,
    Matrix.tail_cons]
  ring

/-- Exact polarized value of the `3+`/`5` cross generator. -/
theorem crossPlusFive_bilinear (μ ν : ℝ) (w : ThreeSectorWeights)
    (x y : CarrierField) :
    invariantKineticBilinear w (crossPlusFiveGen μ ν x) y
        + invariantKineticBilinear w x (crossPlusFiveGen μ ν y)
      = (w.threePlus * μ + w.five * ν) * crossPairingPlus x y := by
  simp only [invariantKineticBilinear, sectorPair, crossPlusFiveGen,
    crossPairingPlus, Fin.sum_univ_three, Fin.sum_univ_five,
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.cons_val_four, Matrix.head_cons,
    Matrix.tail_cons]
  ring

/-- The carrier pair exhibiting the cross pairing: one unit `5`
coordinate against one unit `3-` coordinate. -/
def witnessFive : CarrierField :=
  ⟨fun _ => 0, fun _ => 0, ![1, 0, 0, 0, 0]⟩

/-- The mirror-side exhibiting pair member on the `3-` sector. -/
def witnessMinus : CarrierField :=
  ⟨fun _ => 0, ![1, 0, 0], fun _ => 0⟩

/-- The exhibiting pair member on the `3+` sector. -/
def witnessPlus : CarrierField :=
  ⟨![1, 0, 0], fun _ => 0, fun _ => 0⟩

/-- **Infinitesimal invariance under the `3-`/`5` cross generator is
exactly one linear row in the weights.**  The forward direction
evaluates the polarized identity at the exhibited carrier pair. -/
theorem crossMinusFive_invariance_iff (μ ν : ℝ)
    (w : ThreeSectorWeights) :
    InfinitesimallyInvariant w (crossMinusFiveGen μ ν)
      ↔ w.threeMinus * μ + w.five * ν = 0 := by
  constructor
  · intro hinv
    have h := hinv witnessFive witnessMinus
    rw [crossMinusFive_bilinear] at h
    simpa [crossPairingMinus, witnessFive, witnessMinus,
      Matrix.cons_val_two, Matrix.tail_cons] using h
  · intro h0 x y
    rw [crossMinusFive_bilinear, h0, zero_mul]

/-- **Infinitesimal invariance under the `3+`/`5` cross generator is
exactly one linear row in the weights.** -/
theorem crossPlusFive_invariance_iff (μ ν : ℝ)
    (w : ThreeSectorWeights) :
    InfinitesimallyInvariant w (crossPlusFiveGen μ ν)
      ↔ w.threePlus * μ + w.five * ν = 0 := by
  constructor
  · intro hinv
    have h := hinv witnessFive witnessPlus
    rw [crossPlusFive_bilinear] at h
    simpa [crossPairingPlus, witnessFive, witnessPlus,
      Matrix.cons_val_two, Matrix.tail_cons] using h
  · intro h0 x y
    rw [crossPlusFive_bilinear, h0, zero_mul]

/-- The declared realization of the committed `G` row: the `3-`/`5`
cross generator at coefficient ratio `(1, -sqrt 5)`.  The ratio is
declared so that the invariance row below is exactly the committed
constraint row; its derivation from the pinned bracket tensors is
certificate content. -/
noncomputable def gMirrorGen : CarrierField → CarrierField :=
  crossMinusFiveGen 1 (-Real.sqrt 5)

/-- The declared realization of the committed `F` row: the `3+`/`5`
cross generator at coefficient ratio `(1, -sqrt 5)`. -/
noncomputable def fMirrorGen : CarrierField → CarrierField :=
  crossPlusFiveGen 1 (-Real.sqrt 5)

/-- **The committed `G` constraint row is exactly the infinitesimal
invariance row of the exhibited cross generator.**  The committed
ad-invariance receipt of `GaugeKineticInvariantForms` is upgraded
from a coefficient row to a typed invariance statement of the
kinetic functional, at the exact infinitesimal grade. -/
theorem gMirrorGen_invariance_iff (w : ThreeSectorWeights) :
    InfinitesimallyInvariant w gMirrorGen ↔ GAdInvariant w := by
  unfold gMirrorGen
  rw [crossMinusFive_invariance_iff]
  constructor
  · intro h
    show w.threeMinus = Real.sqrt 5 * w.five
    linear_combination h
  · intro h
    have h' : w.threeMinus = Real.sqrt 5 * w.five := h
    linear_combination h'

/-- **The committed `F` constraint row is exactly the infinitesimal
invariance row of the mirror cross generator.** -/
theorem fMirrorGen_invariance_iff (w : ThreeSectorWeights) :
    InfinitesimallyInvariant w fMirrorGen ↔ FAdInvariant w := by
  unfold fMirrorGen
  rw [crossPlusFive_invariance_iff]
  constructor
  · intro h
    show w.threePlus = Real.sqrt 5 * w.five
    linear_combination h
  · intro h
    have h' : w.threePlus = Real.sqrt 5 * w.five := h
    linear_combination h'

/-- **Simultaneous infinitesimal invariance under both exhibited
mirror generators leaves exactly the committed one-ray family**: the
committed extra-premise one-ray theorem restated at the invariance
grade. -/
theorem mirror_pair_invariance_one_ray (w : ThreeSectorWeights) :
    (InfinitesimallyInvariant w fMirrorGen
        ∧ InfinitesimallyInvariant w gMirrorGen)
      ↔ ∃! sc : ℝ, w = mirrorCommonRay sc := by
  rw [fMirrorGen_invariance_iff, gMirrorGen_invariance_iff]
  exact mirror_common_extra_premise_one_ray w

/-- The committed two-parameter normal form is infinitesimally
invariant under the exhibited `G` cross generator, for every
declared parameter pair. -/
theorem completed_gauge_coefficients_invariant (p f : ℝ) :
    InfinitesimallyInvariant (gNormalForm p f) gMirrorGen :=
  (gMirrorGen_invariance_iff (gNormalForm p f)).mpr
    (assembled_gauge_coefficients_G_invariant p f)

/-! ## Mutation receipts: the invariance is load-bearing -/

theorem sqrt_five_ne_one : Real.sqrt 5 ≠ 1 := by
  intro h
  have h5 : Real.sqrt 5 ^ 2 = 5 := Real.sq_sqrt (by norm_num)
  rw [h] at h5
  norm_num at h5

/-- The control weight vector: it violates the committed `G` row. -/
def controlWeights : ThreeSectorWeights := ⟨1, 1, 1⟩

theorem control_form_not_G_invariant : ¬ GAdInvariant controlWeights := by
  intro h
  have h1 : (1 : ℝ) = Real.sqrt 5 * 1 := h
  rw [mul_one] at h1
  exact sqrt_five_ne_one h1.symm

/-- **Exhibited break.**  At the control weight vector the polarized
identity of the exhibited `G` cross generator evaluates to the
nonzero value `1 - sqrt 5` on the exhibited carrier pair. -/
theorem control_form_breaks_invariance_witness :
    invariantKineticBilinear controlWeights (gMirrorGen witnessFive)
          witnessMinus
        + invariantKineticBilinear controlWeights witnessFive
          (gMirrorGen witnessMinus)
      = 1 - Real.sqrt 5
    ∧ (1 - Real.sqrt 5 : ℝ) ≠ 0 := by
  constructor
  · show invariantKineticBilinear controlWeights
        (crossMinusFiveGen 1 (-Real.sqrt 5) witnessFive) witnessMinus
      + invariantKineticBilinear controlWeights witnessFive
        (crossMinusFiveGen 1 (-Real.sqrt 5) witnessMinus)
      = 1 - Real.sqrt 5
    rw [crossMinusFive_bilinear]
    simp only [crossPairingMinus, witnessFive, witnessMinus,
      controlWeights, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
    ring
  · intro h
    have h1 : Real.sqrt 5 = 1 := by linarith
    exact sqrt_five_ne_one h1

/-- The control weight vector is not infinitesimally invariant under
the exhibited `G` cross generator. -/
theorem control_form_not_invariant :
    ¬ InfinitesimallyInvariant controlWeights gMirrorGen := by
  intro hinv
  have h := hinv witnessFive witnessMinus
  rw [control_form_breaks_invariance_witness.1] at h
  exact control_form_breaks_invariance_witness.2 h

/-- The control generator at ratio `(1, 1)` breaks the invariance of
the committed normal form at the all-ones parameters: the exhibited
ratio `(1, -sqrt 5)` is load-bearing. -/
theorem control_generator_breaks_committed_form :
    ¬ InfinitesimallyInvariant (gNormalForm 1 1)
        (crossMinusFiveGen 1 1) := by
  rw [crossMinusFive_invariance_iff]
  intro h
  have h5 : (0 : ℝ) ≤ Real.sqrt 5 := Real.sqrt_nonneg 5
  simp only [gNormalForm] at h
  norm_num at h
  linarith

/-! ## The completed assembled action -/

/-- The completed field tuple: the committed field tuple with its
gauge strength slot upgraded to the typed carrier slot. -/
structure CompletedFieldTuple where
  /-- Gauge-field slot on the committed bracket carrier. -/
  gauge : CarrierField
  /-- Scalar doublet slot (committed carrier). -/
  scalar : Doublet
  /-- Constant electroweak direction slot. -/
  direction : EWParam
  /-- Yukawa matter doublet slot. -/
  matterDoublet : Doublet
  /-- Yukawa matter singlet slot. -/
  matterSinglet : ℂ

/-- **The completed assembled action.**  The invariant-form kinetic
functional at the committed two-parameter `G` normal form, plus the
committed quartic potential, gauge mass form at the minimum, and
realified Yukawa line, on the unchanged committed premise bundle.
No parameter is added: the invariant-form coefficients are the two
declared bundle weights. -/
noncomputable def completedAction (D : AssembledActionPremiseData)
    (Ψ : CompletedFieldTuple) : ℝ :=
  invariantKineticForm (gNormalForm D.gaugeFreePlus D.gaugeFive)
      Ψ.gauge
    + potential D.ew.lam D.ew.vev Ψ.scalar
    + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev Ψ.direction
    + (yukawaLine D.ew.yuk Ψ.scalar Ψ.matterDoublet
        Ψ.matterSinglet).re

/-- The axis embedding of committed field tuples into completed
ones. -/
def completedOfFieldTuple (Φ : FieldTuple) : CompletedFieldTuple :=
  ⟨axisCarrier Φ.strength, Φ.scalar, Φ.direction, Φ.matterDoublet,
    Φ.matterSinglet⟩

/-- The committed tuple of a completed tuple: sector magnitudes on
the gauge slot, identity elsewhere. -/
noncomputable def fieldTupleOfCompleted (Ψ : CompletedFieldTuple) :
    FieldTuple :=
  ⟨sectorMagnitudes Ψ.gauge, Ψ.scalar, Ψ.direction, Ψ.matterDoublet,
    Ψ.matterSinglet⟩

/-- **The committed assembled action embeds exactly in the completed
one** under the axis identification. -/
theorem completed_embeds_committed (D : AssembledActionPremiseData)
    (Φ : FieldTuple) :
    completedAction D (completedOfFieldTuple Φ)
      = assembledAction D Φ := by
  show invariantKineticForm (gNormalForm D.gaugeFreePlus D.gaugeFive)
        (axisCarrier Φ.strength)
      + potential D.ew.lam D.ew.vev Φ.scalar
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev Φ.direction
      + (yukawaLine D.ew.yuk Φ.scalar Φ.matterDoublet
          Φ.matterSinglet).re = _
  rw [invariantKineticForm_axis]
  rfl

/-- **The completed action factors exactly through the committed
action** at the sector magnitudes. -/
theorem completed_factors_through_committed
    (D : AssembledActionPremiseData) (Ψ : CompletedFieldTuple) :
    completedAction D Ψ
      = assembledAction D (fieldTupleOfCompleted Ψ) := by
  show completedAction D Ψ
      = gaugeKineticEnergy (gNormalForm D.gaugeFreePlus D.gaugeFive)
          (sectorMagnitudes Ψ.gauge)
        + potential D.ew.lam D.ew.vev Ψ.scalar
        + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev Ψ.direction
        + (yukawaLine D.ew.yuk Ψ.scalar Ψ.matterDoublet
            Ψ.matterSinglet).re
  rw [gaugeKineticEnergy_sectorMagnitudes]
  rfl

/-- The completed base tuple: the axis embedding of the committed
base tuple. -/
noncomputable def completedBaseTuple (D : AssembledActionPremiseData) :
    CompletedFieldTuple :=
  completedOfFieldTuple (baseTuple D)

/-- The completed action vanishes at the completed base tuple. -/
theorem completedAction_base (D : AssembledActionPremiseData) :
    completedAction D (completedBaseTuple D) = 0 := by
  unfold completedBaseTuple
  rw [completed_embeds_committed]
  exact assembledAction_base D

/-- **Gauge restriction.**  On the gauge slice through the completed
base tuple the completed action is exactly the invariant-form
kinetic functional. -/
theorem completed_restrict_gauge (D : AssembledActionPremiseData)
    (A : CarrierField) :
    completedAction D { completedBaseTuple D with gauge := A }
      = invariantKineticForm
          (gNormalForm D.gaugeFreePlus D.gaugeFive) A := by
  show invariantKineticForm (gNormalForm D.gaugeFreePlus D.gaugeFive)
        A
      + potential D.ew.lam D.ew.vev (vacuum D.ew.vev)
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev ⟨0, 0, 0, 0⟩
      + (yukawaLine D.ew.yuk (vacuum D.ew.vev) (0 : Doublet)
          (0 : ℂ)).re = _
  rw [potential_vacuum, gaugeMassForm_zero_direction,
    yukawaLine_zero_doublet]
  simp

/-- **Scalar restriction persists.**  On the scalar slice the
completed action is exactly the committed scalar sector. -/
theorem completed_restrict_scalar (D : AssembledActionPremiseData)
    (h : Doublet) (A : EWParam) :
    completedAction D
        { completedBaseTuple D with scalar := h, direction := A }
      = potential D.ew.lam D.ew.vev h
        + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev A := by
  show invariantKineticForm (gNormalForm D.gaugeFreePlus D.gaugeFive)
        (axisCarrier ⟨0, 0, 0⟩)
      + potential D.ew.lam D.ew.vev h
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev A
      + (yukawaLine D.ew.yuk h (0 : Doublet) (0 : ℂ)).re = _
  rw [invariantKineticForm_axis, gaugeKineticEnergy_zero,
    yukawaLine_zero_doublet]
  simp

/-- **Yukawa restriction persists.**  On the matter slice the
completed action is exactly the realified committed Yukawa line at
the chosen minimum. -/
theorem completed_restrict_yukawa (D : AssembledActionPremiseData)
    (ψ : Doublet) (χ : ℂ) :
    completedAction D
        { completedBaseTuple D with
            matterDoublet := ψ, matterSinglet := χ }
      = (yukawaLine D.ew.yuk (vacuum D.ew.vev) ψ χ).re := by
  show invariantKineticForm (gNormalForm D.gaugeFreePlus D.gaugeFive)
        (axisCarrier ⟨0, 0, 0⟩)
      + potential D.ew.lam D.ew.vev (vacuum D.ew.vev)
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev ⟨0, 0, 0, 0⟩
      + (yukawaLine D.ew.yuk (vacuum D.ew.vev) ψ χ).re = _
  rw [invariantKineticForm_axis, gaugeKineticEnergy_zero,
    potential_vacuum, gaugeMassForm_zero_direction]
  simp

/-! ## Explicit display, census, and positivity -/

/-- **The explicit completed sum.**  Every coefficient is a declared
bundle field; the gauge weights are the committed `G` normal form of
the two declared parameters. -/
theorem completedAction_explicit (D : AssembledActionPremiseData)
    (Ψ : CompletedFieldTuple) :
    completedAction D Ψ
      = D.gaugeFreePlus * sectorNormSq Ψ.gauge.threePlus
        + Real.sqrt 5 * D.gaugeFive * sectorNormSq Ψ.gauge.threeMinus
        + D.gaugeFive * sectorNormSq Ψ.gauge.five
        + D.ew.lam * (normSqSum Ψ.scalar - D.ew.vev ^ 2) ^ 2
        + D.ew.vev ^ 2 / 4
            * (D.ew.gW ^ 2
                * (Ψ.direction.a1 ^ 2 + Ψ.direction.a2 ^ 2)
              + (D.ew.gY * Ψ.direction.b
                  - D.ew.gW * Ψ.direction.a3) ^ 2)
        + ((D.ew.yuk : ℂ) * epsPair Ψ.matterDoublet Ψ.scalar
            * Ψ.matterSinglet).re := by
  show invariantKineticForm (gNormalForm D.gaugeFreePlus D.gaugeFive)
        Ψ.gauge
      + potential D.ew.lam D.ew.vev Ψ.scalar
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev Ψ.direction
      + (yukawaLine D.ew.yuk Ψ.scalar Ψ.matterDoublet
          Ψ.matterSinglet).re = _
  rw [gaugeMassForm_exact]
  simp only [invariantKineticForm, invariantKineticBilinear,
    gNormalForm, potential, yukawaLine, sectorNormSq]

/-- **The parameter census is unchanged.**  Two bundles agreeing on
the committed declared coefficient list (the two gauge weights, the
quartic, the scale, the two couplings, the Yukawa coefficient) have
equal completed actions: the completion introduces no coefficient
beyond the committed census. -/
theorem completed_parameter_census
    (D D' : AssembledActionPremiseData)
    (hplus : D.gaugeFreePlus = D'.gaugeFreePlus)
    (hfive : D.gaugeFive = D'.gaugeFive)
    (hlam : D.ew.lam = D'.ew.lam)
    (hvev : D.ew.vev = D'.ew.vev)
    (hgW : D.ew.gW = D'.ew.gW)
    (hgY : D.ew.gY = D'.ew.gY)
    (hyuk : D.ew.yuk = D'.ew.yuk) :
    ∀ Ψ : CompletedFieldTuple,
      completedAction D Ψ = completedAction D' Ψ := by
  intro Ψ
  rw [completedAction_explicit, completedAction_explicit, hplus,
    hfive, hlam, hvev, hgW, hgY, hyuk]

/-- **Exact parameter count of the completed gauge sector.**  The
invariant-form coefficient vector of the completed action lies in
the committed two-parameter family with exactly one parameter pair:
the two declared bundle weights, per the committed
`g_exact_two_parameter`. -/
theorem completed_gauge_parameter_count
    (D : AssembledActionPremiseData) :
    ∃! p : ℝ × ℝ,
      gNormalForm D.gaugeFreePlus D.gaugeFive
        = gNormalForm p.1 p.2 :=
  (g_exact_two_parameter (gNormalForm D.gaugeFreePlus D.gaugeFive)).mp
    (assembled_gauge_coefficients_G_invariant D.gaugeFreePlus
      D.gaugeFive)

/-- The completed gauge sector is nonnegative under the declared
positivity clauses, through the exact factorization. -/
theorem completed_gauge_energy_nonneg (D : AssembledActionPremiseData)
    (A : CarrierField) :
    0 ≤ invariantKineticForm
        (gNormalForm D.gaugeFreePlus D.gaugeFive) A := by
  rw [← gaugeKineticEnergy_sectorMagnitudes]
  exact assembled_gauge_energy_nonneg D (sectorMagnitudes A)

/-- Kinetic-level invariance transports to the completed action: two
gauge slots of equal kinetic value give equal completed actions. -/
theorem completedAction_gauge_transport (D : AssembledActionPremiseData)
    (Ψ : CompletedFieldTuple) (A : CarrierField)
    (h : invariantKineticForm
          (gNormalForm D.gaugeFreePlus D.gaugeFive) A
        = invariantKineticForm
          (gNormalForm D.gaugeFreePlus D.gaugeFive) Ψ.gauge) :
    completedAction D { Ψ with gauge := A } = completedAction D Ψ := by
  show invariantKineticForm (gNormalForm D.gaugeFreePlus D.gaugeFive)
        A
      + potential D.ew.lam D.ew.vev Ψ.scalar
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev Ψ.direction
      + (yukawaLine D.ew.yuk Ψ.scalar Ψ.matterDoublet
          Ψ.matterSinglet).re
    = invariantKineticForm (gNormalForm D.gaugeFreePlus D.gaugeFive)
        Ψ.gauge
      + potential D.ew.lam D.ew.vev Ψ.scalar
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev Ψ.direction
      + (yukawaLine D.ew.yuk Ψ.scalar Ψ.matterDoublet
          Ψ.matterSinglet).re
  rw [h]

/-! ## Abelian phase invariance persists -/

/-- The abelian phase rescaling of the charged slots of a completed
tuple: the gauge carrier and direction slots are uncharged and stay
fixed. -/
def completedScaleTuple (u v w : ℂ) (Ψ : CompletedFieldTuple) :
    CompletedFieldTuple where
  gauge := Ψ.gauge
  scalar := fun i => u * Ψ.scalar i
  direction := Ψ.direction
  matterDoublet := fun i => v * Ψ.matterDoublet i
  matterSinglet := w * Ψ.matterSinglet

theorem fieldTupleOfCompleted_scale (u v w : ℂ)
    (Ψ : CompletedFieldTuple) :
    fieldTupleOfCompleted (completedScaleTuple u v w Ψ)
      = scaleTuple u v w (fieldTupleOfCompleted Ψ) := rfl

/-- **Exact abelian invariance of the completed sum**, inherited
through the exact factorization from the committed invariance. -/
theorem completed_abelian_invariance (D : AssembledActionPremiseData)
    (u v w : ℂ) (hu : Complex.normSq u = 1) (hprod : u * v * w = 1)
    (Ψ : CompletedFieldTuple) :
    completedAction D (completedScaleTuple u v w Ψ)
      = completedAction D Ψ := by
  rw [completed_factors_through_committed,
    completed_factors_through_committed, fieldTupleOfCompleted_scale]
  exact assembled_abelian_invariance D u v w hu hprod
    (fieldTupleOfCompleted Ψ)

/-- **Committed-charge abelian invariance of the completed sum**: the
committed integer charge balance of the bundled selection mask leaves
the completed action exactly fixed at every unit phase. -/
theorem completed_committed_charge_abelian_invariance
    (D : AssembledActionPremiseData) :
    ∃ (cH : ℤ) (dRow sRow : Fin 10),
      (cH = scalarChargeSix ∨ cH = -scalarChargeSix)
        ∧ OPH.ExteriorSelection.mem D.ew.base.selectionMask.val dRow
            = true
        ∧ OPH.ExteriorSelection.mem D.ew.base.selectionMask.val sRow
            = true
        ∧ cH + OPH.ExteriorSelection.charge dRow
            + OPH.ExteriorSelection.charge sRow = 0
        ∧ ∀ p : ℂ, Complex.normSq p = 1 →
            ∀ Ψ : CompletedFieldTuple,
              completedAction D
                  (completedScaleTuple (p ^ cH)
                    (p ^ OPH.ExteriorSelection.charge dRow)
                    (p ^ OPH.ExteriorSelection.charge sRow) Ψ)
                = completedAction D Ψ := by
  obtain ⟨cH, dRow, sRow, hcH, hd, hs, hbal, hinv⟩ :=
    committed_charge_abelian_invariance D
  refine ⟨cH, dRow, sRow, hcH, hd, hs, hbal, fun p hp Ψ => ?_⟩
  rw [completed_factors_through_committed,
    completed_factors_through_committed, fieldTupleOfCompleted_scale]
  exact hinv p hp (fieldTupleOfCompleted Ψ)

/-! ## The balanced-point anchor, re-exported -/

/-- The committed selection anchor of the carrier: at the
sector-balanced reference metric the committed compact family `G` is
the unique nearest classified family (committed content of
`OrientedFaceInvariantMetric`, register rows PR-09, PR-10, PR-11). -/
theorem carrier_family_selected_at_balanced_metric
    (family : OPH.OrientedFaceBracketSelector.CompactFamily)
    (h : family ≠ OPH.OrientedFaceBracketSelector.CompactFamily.G) :
    OPH.OrientedFaceInvariantMetric.dG2 1 1 1
      < OPH.OrientedFaceInvariantMetric.squaredDistanceAt 1 1 1
          family :=
  OPH.OrientedFaceInvariantMetric.balanced_unique_nearest_G
    one_pos one_pos family h

/-! ## Nonvacuity: the completed committed inhabitant -/

/-- The completed premise bundle is the unchanged committed bundle;
it is inhabited by the committed inhabitant. -/
theorem completedPremiseData_nonvacuous :
    Nonempty AssembledActionPremiseData :=
  ⟨committedAssembledActionPremiseData⟩

/-- The committed witness evaluation extends exactly: the axis
embedding of the committed witness tuple evaluates to the committed
value `9/4`. -/
theorem completed_committed_inhabitant_embedded :
    completedAction committedAssembledActionPremiseData
        (completedOfFieldTuple
          ⟨⟨1, 0, 0⟩, vacuum 1, wDirection, ![1, 0], 1⟩) = 9 / 4 := by
  rw [completed_embeds_committed]
  exact committed_inhabitant_value

/-- One off-axis witness: the all-ones `3+` carrier vector, outside
the axis slice, with the committed scalar and matter witness
slots. -/
noncomputable def offAxisWitness : CompletedFieldTuple :=
  ⟨⟨![1, 1, 1], fun _ => 0, fun _ => 0⟩, vacuum 1, wDirection,
    ![1, 0], 1⟩

/-- Exact witness evaluation at the completed committed inhabitant:
the off-axis `3+` carrier vector contributes `3`, the potential `0`,
the mass form `1/4`, the Yukawa line `1`.  A witness value of the
declared parameters, not a prediction. -/
theorem completed_off_axis_witness_value :
    completedAction committedAssembledActionPremiseData offAxisWitness
      = 17 / 4 := by
  rw [completedAction_explicit]
  norm_num [committedAssembledActionPremiseData,
    committedEWBreakingPremiseData, offAxisWitness, sectorNormSq,
    sectorPair, Fin.sum_univ_three, Fin.sum_univ_five, vacuum,
    normSqSum, epsPair, wDirection, Matrix.cons_val_two,
    Matrix.tail_cons, OPH.ElectroweakBreakingComposition.wDirection,
    OPH.AssembledActionComposition.committedAssembledActionPremiseData]

/-- The off-axis witness carrier lies outside the axis slice: the
completion is a proper extension of the committed gauge slot. -/
theorem off_axis_witness_not_axis :
    ∀ s : SectorStrength, offAxisWitness.gauge ≠ axisCarrier s := by
  intro s h
  have h1 := congrArg (fun A : CarrierField => A.threePlus 1) h
  simp only [offAxisWitness, axisCarrier, Matrix.cons_val_one] at h1
  norm_num at h1

/-! ## The composed receipt -/

/-- **The composed gauge-sector completion receipt (issue #735,
OL-G9 partial).**  For every committed premise bundle `D`: the
completed action is one explicit sum with every coefficient a
declared bundle field; the committed assembled action embeds exactly
under the axis identification and the completed action factors
exactly through the committed one at the sector magnitudes; the
gauge slice recovers the invariant-form functional and the scalar
and Yukawa slices recover their committed functionals; the completed
gauge coefficients satisfy the committed `G` row, are infinitesimally
invariant under the exhibited cross generator realizing that row,
collapse to the committed one-ray mirror family under simultaneous
invariance for the mirror generator, and are invariant under every
sector-internal skew-adjoint generator; the abelian phase invariance
and the committed-charge invariance persist; the axis restriction
pins the weights; the exhibited control weight vector breaks the
invariance; the completed base tuple evaluates to zero; and the
completed gauge sector is nonnegative.  Consumed premise rows through
the bundle: PR-09, PR-10, PR-11, PR-59 (committed base), PR-48,
PR-49, PR-50, PR-54 (declared parameter fields; the rows stay open).
Rows PR-35, PR-36, PR-46, PR-47 are untouched. -/
theorem gaugeSectorBracketCompletion_receipt
    (D : AssembledActionPremiseData) :
    (∀ Φ : FieldTuple,
        completedAction D (completedOfFieldTuple Φ)
          = assembledAction D Φ)
    ∧ (∀ Ψ : CompletedFieldTuple,
        completedAction D Ψ
          = assembledAction D (fieldTupleOfCompleted Ψ))
    ∧ ((∀ A : CarrierField,
          completedAction D { completedBaseTuple D with gauge := A }
            = invariantKineticForm
                (gNormalForm D.gaugeFreePlus D.gaugeFive) A)
      ∧ (∀ (h : Doublet) (A : EWParam),
          completedAction D
              { completedBaseTuple D with
                  scalar := h, direction := A }
            = potential D.ew.lam D.ew.vev h
              + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev A)
      ∧ (∀ (ψ : Doublet) (χ : ℂ),
          completedAction D
              { completedBaseTuple D with
                  matterDoublet := ψ, matterSinglet := χ }
            = (yukawaLine D.ew.yuk (vacuum D.ew.vev) ψ χ).re)
      ∧ completedAction D (completedBaseTuple D) = 0)
    ∧ (GAdInvariant (gNormalForm D.gaugeFreePlus D.gaugeFive)
      ∧ InfinitesimallyInvariant
          (gNormalForm D.gaugeFreePlus D.gaugeFive) gMirrorGen
      ∧ (InfinitesimallyInvariant
            (gNormalForm D.gaugeFreePlus D.gaugeFive) fMirrorGen →
          ∃! sc : ℝ,
            gNormalForm D.gaugeFreePlus D.gaugeFive
              = mirrorCommonRay sc)
      ∧ (∀ L : CarrierField → CarrierField, SectorInternalSkew L →
          InfinitesimallyInvariant
            (gNormalForm D.gaugeFreePlus D.gaugeFive) L))
    ∧ ((∀ (u v w : ℂ), Complex.normSq u = 1 → u * v * w = 1 →
          ∀ Ψ : CompletedFieldTuple,
            completedAction D (completedScaleTuple u v w Ψ)
              = completedAction D Ψ)
      ∧ (∃ (cH : ℤ) (dRow sRow : Fin 10),
          (cH = scalarChargeSix ∨ cH = -scalarChargeSix)
            ∧ OPH.ExteriorSelection.mem D.ew.base.selectionMask.val
                dRow = true
            ∧ OPH.ExteriorSelection.mem D.ew.base.selectionMask.val
                sRow = true
            ∧ cH + OPH.ExteriorSelection.charge dRow
                + OPH.ExteriorSelection.charge sRow = 0
            ∧ ∀ p : ℂ, Complex.normSq p = 1 →
                ∀ Ψ : CompletedFieldTuple,
                  completedAction D
                      (completedScaleTuple (p ^ cH)
                        (p ^ OPH.ExteriorSelection.charge dRow)
                        (p ^ OPH.ExteriorSelection.charge sRow) Ψ)
                    = completedAction D Ψ))
    ∧ ((∀ w w' : ThreeSectorWeights,
          (∀ s : SectorStrength,
            invariantKineticForm w (axisCarrier s)
              = invariantKineticForm w' (axisCarrier s)) → w = w')
      ∧ ¬ InfinitesimallyInvariant controlWeights gMirrorGen
      ∧ (∀ A : CarrierField,
          0 ≤ invariantKineticForm
              (gNormalForm D.gaugeFreePlus D.gaugeFive) A)) :=
  ⟨fun Φ => completed_embeds_committed D Φ,
    fun Ψ => completed_factors_through_committed D Ψ,
    ⟨fun A => completed_restrict_gauge D A,
      fun h A => completed_restrict_scalar D h A,
      fun ψ χ => completed_restrict_yukawa D ψ χ,
      completedAction_base D⟩,
    ⟨assembled_gauge_coefficients_G_invariant D.gaugeFreePlus
        D.gaugeFive,
      completed_gauge_coefficients_invariant D.gaugeFreePlus
        D.gaugeFive,
      fun hF => (mirror_pair_invariance_one_ray
          (gNormalForm D.gaugeFreePlus D.gaugeFive)).mp
        ⟨hF, completed_gauge_coefficients_invariant D.gaugeFreePlus
          D.gaugeFive⟩,
      fun L hL => sector_internal_skew_invariant
        (gNormalForm D.gaugeFreePlus D.gaugeFive) L hL⟩,
    ⟨fun u v w hu hprod Ψ =>
        completed_abelian_invariance D u v w hu hprod Ψ,
      completed_committed_charge_abelian_invariance D⟩,
    ⟨fun w w' h => axis_restriction_pins w w' h,
      control_form_not_invariant,
      fun A => completed_gauge_energy_nonneg D A⟩⟩

end OPH.GaugeSectorBracketCompletion

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.GaugeSectorBracketCompletion.invariantKineticForm_axis
#print axioms OPH.GaugeSectorBracketCompletion.gaugeKineticEnergy_sectorMagnitudes
#print axioms OPH.GaugeSectorBracketCompletion.axis_restriction_pins
#print axioms OPH.GaugeSectorBracketCompletion.sector_internal_skew_invariant
#print axioms OPH.GaugeSectorBracketCompletion.rotPlusGen_sectorInternalSkew
#print axioms OPH.GaugeSectorBracketCompletion.crossMinusFive_bilinear
#print axioms OPH.GaugeSectorBracketCompletion.crossPlusFive_bilinear
#print axioms OPH.GaugeSectorBracketCompletion.crossMinusFive_invariance_iff
#print axioms OPH.GaugeSectorBracketCompletion.crossPlusFive_invariance_iff
#print axioms OPH.GaugeSectorBracketCompletion.gMirrorGen_invariance_iff
#print axioms OPH.GaugeSectorBracketCompletion.fMirrorGen_invariance_iff
#print axioms OPH.GaugeSectorBracketCompletion.mirror_pair_invariance_one_ray
#print axioms OPH.GaugeSectorBracketCompletion.completed_gauge_coefficients_invariant
#print axioms OPH.GaugeSectorBracketCompletion.control_form_not_G_invariant
#print axioms OPH.GaugeSectorBracketCompletion.control_form_breaks_invariance_witness
#print axioms OPH.GaugeSectorBracketCompletion.control_form_not_invariant
#print axioms OPH.GaugeSectorBracketCompletion.control_generator_breaks_committed_form
#print axioms OPH.GaugeSectorBracketCompletion.completed_embeds_committed
#print axioms OPH.GaugeSectorBracketCompletion.completed_factors_through_committed
#print axioms OPH.GaugeSectorBracketCompletion.completedAction_base
#print axioms OPH.GaugeSectorBracketCompletion.completed_restrict_gauge
#print axioms OPH.GaugeSectorBracketCompletion.completed_restrict_scalar
#print axioms OPH.GaugeSectorBracketCompletion.completed_restrict_yukawa
#print axioms OPH.GaugeSectorBracketCompletion.completedAction_explicit
#print axioms OPH.GaugeSectorBracketCompletion.completed_parameter_census
#print axioms OPH.GaugeSectorBracketCompletion.completed_gauge_parameter_count
#print axioms OPH.GaugeSectorBracketCompletion.completed_gauge_energy_nonneg
#print axioms OPH.GaugeSectorBracketCompletion.completedAction_gauge_transport
#print axioms OPH.GaugeSectorBracketCompletion.completed_abelian_invariance
#print axioms OPH.GaugeSectorBracketCompletion.completed_committed_charge_abelian_invariance
#print axioms OPH.GaugeSectorBracketCompletion.carrier_family_selected_at_balanced_metric
#print axioms OPH.GaugeSectorBracketCompletion.completedPremiseData_nonvacuous
#print axioms OPH.GaugeSectorBracketCompletion.completed_committed_inhabitant_embedded
#print axioms OPH.GaugeSectorBracketCompletion.completed_off_axis_witness_value
#print axioms OPH.GaugeSectorBracketCompletion.off_axis_witness_not_axis
#print axioms OPH.GaugeSectorBracketCompletion.gaugeSectorBracketCompletion_receipt
