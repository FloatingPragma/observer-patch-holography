import GaugeKineticInvariantForms
import ElectroweakBreakingComposition

/-!
# The first assembled structural action (V3, issue #735, OL-G9)

One typed action object on one committed carrier bundle.  The
committed term inventory of `tracking/sm_lagrangian_correspondence.json`
records sectors one by one; this module sums three of them as one
real-valued function of one field tuple and proves the composition
receipts on one typed antecedent bundle.

The three summed sectors:

* **Gauge kinetic sector.**  The coefficient object is the committed
  `ThreeSectorWeights` of `GaugeKineticInvariantForms`, taken in the
  committed two-parameter ad-invariant normal form `gNormalForm` of
  the measure-selected compact family `G` (the family that
  `SMStructureComposition` selects under register rows PR-09, PR-10,
  PR-11).  The two normal-form parameters are declared bundle fields
  with declared positivity.  The pairing with a field slot is the
  diagonal sector pairing: one declared strength magnitude per
  invariant carrier sector (`3+`, `3-`, `5`), each squared against
  its committed weight.  The pairing convention is an assembly
  declaration of this file; the committed module supplies the weight
  space, the constraint rows, and the normal forms, and this file
  reuses them without restating the constraint content.
* **Scalar sector.**  The committed quartic potential and the
  committed gauge mass form of `ElectroweakBreakingComposition`: the
  potential is evaluated on the doublet slot, and the mass form (the
  exact constant-direction coefficient of a covariant kinetic term at
  the chosen minimum) is evaluated on an electroweak direction slot.
* **Yukawa sector.**  The committed one-line `yukawaLine`, realified:
  the assembled action carries its real part.  The complex line is
  recovered exactly from two assembled evaluations
  (`yukawa_complex_recovery`), so the realification loses nothing.

Proved on the one bundle `AssembledActionPremiseData`:

1. the assembled action is one well-typed real function of one field
   tuple on the committed carriers, with every coefficient a declared
   bundle field (`assembledAction_explicit` displays the sum with all
   coefficients);
2. three exact restriction theorems: on each sector's field slice
   (the other slots held at the base tuple) the assembled action
   equals exactly that sector's committed functional;
3. invariance of the sum, stated exactly at the strength the
   committed structure supports: the assembled gauge coefficients
   satisfy the committed `G` ad-invariance constraint row, the extra
   declared `F` row collapses them to the committed one-ray mirror
   form, and the whole sum is exactly invariant under the global
   abelian phase action whose integer charges are the committed
   exterior-table charges at the committed charge balance
   (`committed_charge_abelian_invariance`).  No local gauge
   invariance and no nonabelian invariance of the sum is claimed:
   the committed invariant-form receipts are constraint rows on
   coefficients, and the committed exterior grammar supplies abelian
   integer charges only;
4. an explicit inhabitant with every declared parameter equal to one,
   a nonvacuity receipt, and one exact witness evaluation;
5. one composed receipt (`assembledActionComposition_receipt`) whose
   antecedent bundle names every consumed premise row.

Consumed register rows.  The bundle field `ew` is the committed
`EWBreakingPremiseData`, which carries the composed
Standard-Model-structure premises of register rows PR-09, PR-10,
PR-11, PR-59 through its base, the declared scalar potential
parameters of row PR-49, the declared coupling constants whose
selection stays on rows PR-48 and PR-54, and the declared Yukawa
coefficient of row PR-50.  The two gauge weight parameters are
declared here under the classification row PR-11, matching the
committed boundary of inventory rows SML-01 and SML-02: the
Yang-Mills shape is attained on the classified branch and the
coefficient values are free.

Boundary and nonclaims.  This is an assembled finite structural
action on the committed coarse grammar.  It is not the physical
Standard Model action.  There is no spacetime, no field strength
two-form, no covariant derivative, no source-produced connection or
current (register row PR-54 stays open), no quantization, no measured
coefficient, and no three-family structure (row PR-36 untouched).
The gauge kinetic sector enters through one declared magnitude per
invariant sector; the full Yang-Mills contraction over a spacetime
carrier is exactly the structure that row PR-54 names and is absent
here.  The scalar covariant kinetic term enters through its committed
constant-direction coefficient at the chosen minimum; a kinetic term
on varying fields is likewise on row PR-54.  Global-form, physical
chirality, and Spin attachments stay on rows PR-35, PR-46, PR-47.
The declared parameters select nothing.  OL-G9 moves from owed to
partial at most, with the assembly named as structural and
conditional.

Falsifier.  The composition fails if any sector restriction differs
from the committed sector functional, if the assembled gauge
coefficients violate the committed `G` constraint row, if the abelian
phase action at the committed charge balance changes the value of the
assembled action, or if the base tuple fails to evaluate to zero.
-/

namespace OPH.AssembledActionComposition

open OPH.GaugeKineticInvariantForms (ThreeSectorWeights FAdInvariant
  GAdInvariant gNormalForm mirrorCommonRay
  mirror_common_extra_premise_one_ray)
open OPH.ElectroweakBreakingComposition (Doublet normSqSum vacuum
  potential gaugeMassForm EWParam yukawaLine epsPair
  EWBreakingPremiseData committedEWBreakingPremiseData
  potential_vacuum vacuum_normSqSum gaugeMassForm_exact
  yukawa_mass_at_vacuum wDirection wMass_value scalarChargeSix
  selected_sector_supports_yukawa_line)

/-! ## The committed carrier bundle: one field tuple -/

/-- One declared field-strength magnitude per committed invariant
carrier sector (`3+`, `3-`, `5`).  The sector split is the committed
split of `ThreeSectorWeights`; the magnitudes are the finite stand-in
for a field strength, since no spacetime field strength exists on the
committed structure (register row PR-54 open). -/
@[ext]
structure SectorStrength where
  /-- Magnitude on the `3+` carrier sector. -/
  threePlus : ℝ
  /-- Magnitude on the `3-` carrier sector. -/
  threeMinus : ℝ
  /-- Magnitude on the `5` carrier sector. -/
  five : ℝ

/-- The diagonal sector pairing of a committed weight vector with a
sector strength: each squared magnitude against its committed weight.
This pairing convention is the assembly declaration of this file. -/
def gaugeKineticEnergy (w : ThreeSectorWeights) (s : SectorStrength) :
    ℝ :=
  w.threePlus * s.threePlus ^ 2 + w.threeMinus * s.threeMinus ^ 2
    + w.five * s.five ^ 2

/-- The one committed field tuple of the assembled action: gauge
sector strengths, the scalar doublet carrier, one constant
electroweak direction, and the Yukawa matter pair (doublet carrier
and singlet) on the committed coarse grammar. -/
structure FieldTuple where
  /-- Gauge kinetic sector slot. -/
  strength : SectorStrength
  /-- Scalar doublet slot (committed carrier `Fin 2 → ℂ`). -/
  scalar : Doublet
  /-- Constant electroweak direction slot for the covariant-kinetic
  mass form at the minimum. -/
  direction : EWParam
  /-- Yukawa matter doublet slot. -/
  matterDoublet : Doublet
  /-- Yukawa matter singlet slot. -/
  matterSinglet : ℂ

/-! ## The typed antecedent bundle -/

/-- The assembled-action premise bundle.  The field `ew` is the
committed electroweak-breaking bundle (register rows PR-09, PR-10,
PR-11, PR-59 through its base; declared rows PR-48, PR-49, PR-50,
PR-54 parameters).  The two gauge fields are the declared parameters
of the committed two-parameter ad-invariant normal form of the
selected family `G` (classification row PR-11; coefficient values
free per inventory rows SML-01 and SML-02), with declared
positivity as the kinetic-positivity clause.  Every field is a
declared hypothesis; none is derived here. -/
structure AssembledActionPremiseData where
  /-- The committed electroweak-breaking premise bundle. -/
  ew : EWBreakingPremiseData
  /-- Declared free `3+` weight of the `G` ad-invariant normal
  form. -/
  gaugeFreePlus : ℝ
  /-- Declared positivity of the `3+` weight. -/
  gaugeFreePlus_pos : 0 < gaugeFreePlus
  /-- Declared free `5` weight of the `G` ad-invariant normal
  form. -/
  gaugeFive : ℝ
  /-- Declared positivity of the `5` weight. -/
  gaugeFive_pos : 0 < gaugeFive

/-! ## The assembled action: one function of one field tuple -/

/-- **The assembled structural action.**  One real-valued function of
one field tuple: the committed `G`-invariant gauge kinetic pairing,
plus the committed quartic potential, plus the committed gauge mass
form at the chosen minimum, plus the realified committed Yukawa
line.  Every coefficient is a declared bundle field. -/
noncomputable def assembledAction (D : AssembledActionPremiseData)
    (Φ : FieldTuple) : ℝ :=
  gaugeKineticEnergy (gNormalForm D.gaugeFreePlus D.gaugeFive)
      Φ.strength
    + potential D.ew.lam D.ew.vev Φ.scalar
    + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev Φ.direction
    + (yukawaLine D.ew.yuk Φ.scalar Φ.matterDoublet
        Φ.matterSinglet).re

/-- The base tuple: zero strengths, the chosen scalar minimum, the
zero direction, and zero matter slots.  Each sector's committed
functional vanishes here, so it is the common base point of the three
sector slices. -/
noncomputable def baseTuple (D : AssembledActionPremiseData) :
    FieldTuple where
  strength := ⟨0, 0, 0⟩
  scalar := vacuum D.ew.vev
  direction := ⟨0, 0, 0, 0⟩
  matterDoublet := 0
  matterSinglet := 0

/-! ## Sector vanishing lemmas at the base point -/

theorem gaugeKineticEnergy_zero (w : ThreeSectorWeights) :
    gaugeKineticEnergy w ⟨0, 0, 0⟩ = 0 := by
  simp [gaugeKineticEnergy]

theorem gaugeMassForm_zero_direction (g g' vev : ℝ) :
    gaugeMassForm g g' vev ⟨0, 0, 0, 0⟩ = 0 := by
  rw [gaugeMassForm_exact]
  norm_num

theorem yukawaLine_zero_doublet (y : ℝ) (h : Doublet) (χ : ℂ) :
    yukawaLine y h (0 : Doublet) χ = 0 := by
  simp [yukawaLine, epsPair]

/-! ## Sector decomposition: three exact restriction theorems -/

/-- **Gauge restriction.**  On the gauge field slice the assembled
action is exactly the committed `G`-invariant kinetic pairing. -/
theorem restrict_gauge (D : AssembledActionPremiseData)
    (s : SectorStrength) :
    assembledAction D { baseTuple D with strength := s }
      = gaugeKineticEnergy (gNormalForm D.gaugeFreePlus D.gaugeFive)
          s := by
  show gaugeKineticEnergy (gNormalForm D.gaugeFreePlus D.gaugeFive) s
      + potential D.ew.lam D.ew.vev (vacuum D.ew.vev)
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev ⟨0, 0, 0, 0⟩
      + (yukawaLine D.ew.yuk (vacuum D.ew.vev) (0 : Doublet)
          (0 : ℂ)).re = _
  rw [potential_vacuum, gaugeMassForm_zero_direction,
    yukawaLine_zero_doublet]
  simp

/-- **Scalar restriction.**  On the scalar field slice the assembled
action is exactly the committed wave-1 scalar sector: the quartic
potential plus the gauge mass form at the minimum. -/
theorem restrict_scalar (D : AssembledActionPremiseData)
    (h : Doublet) (A : EWParam) :
    assembledAction D { baseTuple D with scalar := h, direction := A }
      = potential D.ew.lam D.ew.vev h
        + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev A := by
  show gaugeKineticEnergy (gNormalForm D.gaugeFreePlus D.gaugeFive)
        ⟨0, 0, 0⟩
      + potential D.ew.lam D.ew.vev h
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev A
      + (yukawaLine D.ew.yuk h (0 : Doublet) (0 : ℂ)).re = _
  rw [gaugeKineticEnergy_zero, yukawaLine_zero_doublet]
  simp

/-- **Yukawa restriction.**  On the matter field slice the assembled
action is exactly the real part of the committed Yukawa line at the
chosen minimum. -/
theorem restrict_yukawa (D : AssembledActionPremiseData)
    (ψ : Doublet) (χ : ℂ) :
    assembledAction D
        { baseTuple D with matterDoublet := ψ, matterSinglet := χ }
      = (yukawaLine D.ew.yuk (vacuum D.ew.vev) ψ χ).re := by
  show gaugeKineticEnergy (gNormalForm D.gaugeFreePlus D.gaugeFive)
        ⟨0, 0, 0⟩
      + potential D.ew.lam D.ew.vev (vacuum D.ew.vev)
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev ⟨0, 0, 0, 0⟩
      + (yukawaLine D.ew.yuk (vacuum D.ew.vev) ψ χ).re = _
  rw [gaugeKineticEnergy_zero, potential_vacuum,
    gaugeMassForm_zero_direction]
  simp

/-- The Yukawa restriction carries the exact mass bilinear
coefficient `yuk * vev` of the committed line. -/
theorem restrict_yukawa_mass_coefficient
    (D : AssembledActionPremiseData) (ψ : Doublet) (χ : ℂ) :
    assembledAction D
        { baseTuple D with matterDoublet := ψ, matterSinglet := χ }
      = D.ew.yuk * D.ew.vev * (ψ 0 * χ).re := by
  rw [restrict_yukawa, yukawa_mass_at_vacuum]
  rw [Complex.re_ofReal_mul]

/-- The assembled action vanishes at the base tuple: the common
vacuum normalization of the three sectors. -/
theorem assembledAction_base (D : AssembledActionPremiseData) :
    assembledAction D (baseTuple D) = 0 := by
  show gaugeKineticEnergy (gNormalForm D.gaugeFreePlus D.gaugeFive)
        ⟨0, 0, 0⟩
      + potential D.ew.lam D.ew.vev (vacuum D.ew.vev)
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev ⟨0, 0, 0, 0⟩
      + (yukawaLine D.ew.yuk (vacuum D.ew.vev) (0 : Doublet)
          (0 : ℂ)).re = 0
  rw [gaugeKineticEnergy_zero, potential_vacuum,
    gaugeMassForm_zero_direction, yukawaLine_zero_doublet]
  simp

/-! ## Well-typedness display: the one sum with declared
coefficients -/

/-- **The explicit assembled sum.**  The assembled action written out
as one polynomial-plus-realified-line expression, with every
coefficient a declared bundle field: the gauge weights
`gaugeFreePlus`, `sqrt 5 * gaugeFive`, `gaugeFive` (the committed `G`
normal form), the quartic `lam`, the couplings `gW`, `gY` with the
scale `vev`, and the Yukawa coefficient `yuk`. -/
theorem assembledAction_explicit (D : AssembledActionPremiseData)
    (Φ : FieldTuple) :
    assembledAction D Φ
      = D.gaugeFreePlus * Φ.strength.threePlus ^ 2
        + Real.sqrt 5 * D.gaugeFive * Φ.strength.threeMinus ^ 2
        + D.gaugeFive * Φ.strength.five ^ 2
        + D.ew.lam * (normSqSum Φ.scalar - D.ew.vev ^ 2) ^ 2
        + D.ew.vev ^ 2 / 4
            * (D.ew.gW ^ 2
                * (Φ.direction.a1 ^ 2 + Φ.direction.a2 ^ 2)
              + (D.ew.gY * Φ.direction.b
                  - D.ew.gW * Φ.direction.a3) ^ 2)
        + ((D.ew.yuk : ℂ) * epsPair Φ.matterDoublet Φ.scalar
            * Φ.matterSinglet).re := by
  show gaugeKineticEnergy (gNormalForm D.gaugeFreePlus D.gaugeFive)
        Φ.strength
      + potential D.ew.lam D.ew.vev Φ.scalar
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev Φ.direction
      + (yukawaLine D.ew.yuk Φ.scalar Φ.matterDoublet
          Φ.matterSinglet).re = _
  rw [gaugeMassForm_exact]
  simp only [gaugeKineticEnergy, gNormalForm, potential, yukawaLine]

/-- The assembled gauge kinetic sector is nonnegative under the
declared positivity clauses. -/
theorem assembled_gauge_energy_nonneg (D : AssembledActionPremiseData)
    (s : SectorStrength) :
    0 ≤ gaugeKineticEnergy (gNormalForm D.gaugeFreePlus D.gaugeFive)
        s := by
  have h5 : (0 : ℝ) ≤ Real.sqrt 5 := Real.sqrt_nonneg 5
  simp only [gaugeKineticEnergy, gNormalForm]
  have t1 : 0 ≤ D.gaugeFreePlus * s.threePlus ^ 2 :=
    mul_nonneg D.gaugeFreePlus_pos.le (sq_nonneg _)
  have t2 : 0 ≤ Real.sqrt 5 * D.gaugeFive * s.threeMinus ^ 2 :=
    mul_nonneg (mul_nonneg h5 D.gaugeFive_pos.le) (sq_nonneg _)
  have t3 : 0 ≤ D.gaugeFive * s.five ^ 2 :=
    mul_nonneg D.gaugeFive_pos.le (sq_nonneg _)
  linarith

/-! ## Invariance of the sum, at the committed strength -/

/-- The assembled gauge coefficients satisfy the committed `G`
ad-invariance constraint row: the invariance the gauge sector carries
individually persists for the coefficients of the sum. -/
theorem assembled_gauge_coefficients_G_invariant (p f : ℝ) :
    GAdInvariant (gNormalForm p f) := by
  simp [GAdInvariant, gNormalForm]

/-- Under the extra declared `F` constraint row, the assembled gauge
coefficients collapse to the committed one-ray mirror form.  The
simultaneous premise is declared, exactly as in the committed
module. -/
theorem assembled_gauge_mirror_under_F (p f : ℝ)
    (hF : FAdInvariant (gNormalForm p f)) :
    ∃! sc : ℝ, gNormalForm p f = mirrorCommonRay sc :=
  (mirror_common_extra_premise_one_ray (gNormalForm p f)).mp
    ⟨hF, assembled_gauge_coefficients_G_invariant p f⟩

/-- On the mirror ray the gauge kinetic sector has the committed
one-scale shape. -/
theorem gaugeKineticEnergy_mirror_ray (sc : ℝ) (s : SectorStrength) :
    gaugeKineticEnergy (mirrorCommonRay sc) s
      = sc * (Real.sqrt 5 * (s.threePlus ^ 2 + s.threeMinus ^ 2)
          + s.five ^ 2) := by
  simp only [gaugeKineticEnergy, mirrorCommonRay]
  ring

/-- The abelian phase rescaling of the charged slots of a field
tuple: one unit scalar per charged carrier (scalar doublet, matter
doublet, matter singlet).  The gauge strength and direction slots are
uncharged under the abelian phase and stay fixed. -/
def scaleTuple (u v w : ℂ) (Φ : FieldTuple) : FieldTuple where
  strength := Φ.strength
  scalar := fun i => u * Φ.scalar i
  direction := Φ.direction
  matterDoublet := fun i => v * Φ.matterDoublet i
  matterSinglet := w * Φ.matterSinglet

theorem normSqSum_scale (u : ℂ) (h : Doublet) :
    normSqSum (fun i => u * h i)
      = Complex.normSq u * normSqSum h := by
  simp only [normSqSum, Complex.normSq_mul]
  ring

theorem potential_scale (lam vev : ℝ) (u : ℂ)
    (hu : Complex.normSq u = 1) (h : Doublet) :
    potential lam vev (fun i => u * h i) = potential lam vev h := by
  simp only [potential, normSqSum_scale, hu, one_mul]

theorem epsPair_scale (u v : ℂ) (ψ h : Doublet) :
    epsPair (fun i => v * ψ i) (fun i => u * h i)
      = u * v * epsPair ψ h := by
  simp only [epsPair]
  ring

theorem yukawaLine_scale (y : ℝ) (u v w : ℂ) (h ψ : Doublet)
    (χ : ℂ) :
    yukawaLine y (fun i => u * h i) (fun i => v * ψ i)
        (w * χ)
      = u * v * w * yukawaLine y h ψ χ := by
  simp only [yukawaLine, epsPair_scale]
  ring

/-- **Exact abelian invariance of the sum.**  For any unit scalar on
the scalar slot and any scalar pair on the matter slots whose triple
product is one, the assembled action is exactly invariant.  This is
the invariance the committed structure supports for the sum: a global
abelian phase statement.  No local gauge invariance and no nonabelian
invariance of the sum is claimed. -/
theorem assembled_abelian_invariance (D : AssembledActionPremiseData)
    (u v w : ℂ) (hu : Complex.normSq u = 1) (hprod : u * v * w = 1)
    (Φ : FieldTuple) :
    assembledAction D (scaleTuple u v w Φ) = assembledAction D Φ := by
  show gaugeKineticEnergy (gNormalForm D.gaugeFreePlus D.gaugeFive)
        Φ.strength
      + potential D.ew.lam D.ew.vev (fun i => u * Φ.scalar i)
      + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev Φ.direction
      + (yukawaLine D.ew.yuk (fun i => u * Φ.scalar i)
          (fun i => v * Φ.matterDoublet i)
          (w * Φ.matterSinglet)).re = _
  rw [potential_scale D.ew.lam D.ew.vev u hu, yukawaLine_scale,
    hprod, one_mul]
  rfl

/-- Integer-charge specialization: powers of one unit phase, with
charges summing to zero, leave the assembled action fixed. -/
theorem zpow_charge_invariance (D : AssembledActionPremiseData)
    (a b c : ℤ) (hsum : a + b + c = 0) (p : ℂ)
    (hp : Complex.normSq p = 1) (Φ : FieldTuple) :
    assembledAction D (scaleTuple (p ^ a) (p ^ b) (p ^ c) Φ)
      = assembledAction D Φ := by
  have hp0 : p ≠ 0 := by
    intro h
    rw [h] at hp
    simp at hp
  have hu : Complex.normSq (p ^ a) = 1 := by
    rw [map_zpow₀ Complex.normSq p a, hp, one_zpow]
  have hprod : p ^ a * p ^ b * p ^ c = 1 := by
    rw [← zpow_add₀ hp0, ← zpow_add₀ hp0, hsum, zpow_zero]
  exact assembled_abelian_invariance D _ _ _ hu hprod Φ

/-- **Committed-charge abelian invariance.**  The bundled selection
mask supplies a matter doublet row and singlet row of the committed
exterior table whose integer charges balance the scalar row charge or
its negative; at that committed charge assignment, every unit phase
acting through the integer charges leaves the assembled action
exactly fixed.  The charges are the committed `q = 6Y` table entries;
the scalar slot carries plus or minus the committed scalar row
charge. -/
theorem committed_charge_abelian_invariance
    (D : AssembledActionPremiseData) :
    ∃ (cH : ℤ) (dRow sRow : Fin 10),
      (cH = scalarChargeSix ∨ cH = -scalarChargeSix)
        ∧ OPH.ExteriorSelection.mem D.ew.base.selectionMask.val dRow
            = true
        ∧ OPH.ExteriorSelection.mem D.ew.base.selectionMask.val sRow
            = true
        ∧ cH + OPH.ExteriorSelection.charge dRow
            + OPH.ExteriorSelection.charge sRow = 0
        ∧ ∀ p : ℂ, Complex.normSq p = 1 → ∀ Φ : FieldTuple,
            assembledAction D
                (scaleTuple (p ^ cH)
                  (p ^ OPH.ExteriorSelection.charge dRow)
                  (p ^ OPH.ExteriorSelection.charge sRow) Φ)
              = assembledAction D Φ := by
  obtain ⟨dRow, sRow, hmd, hms, -, -, -, -, -, -, hbal⟩ :=
    selected_sector_supports_yukawa_line D.ew
  rcases hbal with hbal | hbal
  · exact ⟨scalarChargeSix, dRow, sRow, Or.inl rfl, hmd, hms,
      by omega,
      fun p hp Φ => zpow_charge_invariance D _ _ _ (by omega) p hp Φ⟩
  · exact ⟨-scalarChargeSix, dRow, sRow, Or.inr rfl, hmd, hms,
      by omega,
      fun p hp Φ => zpow_charge_invariance D _ _ _ (by omega) p hp Φ⟩

/-! ## Exact complex recovery of the realified Yukawa line -/

/-- The realification loses nothing: the imaginary part of the
committed complex Yukawa line is itself an assembled evaluation, at
the singlet slot rotated by a quarter phase.  Together with
`restrict_yukawa` this recovers the full committed complex line from
two assembled values. -/
theorem yukawa_complex_recovery (D : AssembledActionPremiseData)
    (ψ : Doublet) (χ : ℂ) :
    assembledAction D
        { baseTuple D with
            matterDoublet := ψ, matterSinglet := -Complex.I * χ }
      = (yukawaLine D.ew.yuk (vacuum D.ew.vev) ψ χ).im := by
  rw [restrict_yukawa]
  have h1 : yukawaLine D.ew.yuk (vacuum D.ew.vev) ψ
        (-Complex.I * χ)
      = -Complex.I
          * yukawaLine D.ew.yuk (vacuum D.ew.vev) ψ χ := by
    simp only [yukawaLine]
    ring
  rw [h1]
  simp [Complex.mul_re]

/-! ## Nonvacuity: the committed inhabitant -/

/-- The committed inhabitant: the committed electroweak-breaking
instance with every declared gauge weight set to one.  The instance
is a nonvacuity control; the parameter values select nothing. -/
noncomputable def committedAssembledActionPremiseData :
    AssembledActionPremiseData where
  ew := committedEWBreakingPremiseData
  gaugeFreePlus := 1
  gaugeFreePlus_pos := one_pos
  gaugeFive := 1
  gaugeFive_pos := one_pos

/-- The premise bundle is satisfiable. -/
theorem assembledActionPremiseData_nonvacuous :
    Nonempty AssembledActionPremiseData :=
  ⟨committedAssembledActionPremiseData⟩

/-- Exact witness evaluation at the committed inhabitant: unit `3+`
strength, the chosen minimum, the charged direction, and a unit
matter pair give the value `9/4` (gauge `1`, potential `0`, mass form
`1/4`, Yukawa `1`).  A witness value of the declared parameters, not
a prediction. -/
theorem committed_inhabitant_value :
    assembledAction committedAssembledActionPremiseData
        ⟨⟨1, 0, 0⟩, vacuum 1, wDirection, ![1, 0], 1⟩ = 9 / 4 := by
  rw [assembledAction_explicit]
  norm_num [committedAssembledActionPremiseData,
    committedEWBreakingPremiseData, wDirection, vacuum, normSqSum,
    epsPair,
    OPH.ElectroweakBreakingComposition.wDirection]

/-! ## The composed receipt -/

/-- **The composed assembled-action receipt (issue #735, OL-G9
partial).**  For every premise bundle `D`: the assembled action is
one explicit sum with every coefficient a declared bundle field; its
restriction to each sector's field slice is exactly that sector's
committed functional (gauge kinetic pairing, potential plus gauge
mass form, realified Yukawa line at the minimum); the assembled gauge
coefficients satisfy the committed `G` ad-invariance constraint row
and collapse to the committed one-ray mirror form under the extra
declared `F` row; the sum is exactly invariant under the global
abelian phase action at unit scalars with product one, and in
particular at the committed integer charge balance of the bundled
selection mask; the imaginary part of the committed complex Yukawa
line is recovered as an assembled evaluation; the base tuple
evaluates to zero; and the gauge sector is nonnegative under the
declared positivity clauses.  Every clause consumes the bundle.
Consumed premise rows through the bundle: PR-09, PR-10, PR-11, PR-59
(committed base), PR-48, PR-49, PR-50, PR-54 (declared parameter
fields; the rows stay open).  Rows PR-35, PR-36, PR-46, PR-47 are
untouched. -/
theorem assembledActionComposition_receipt
    (D : AssembledActionPremiseData) :
    (∀ Φ : FieldTuple,
        assembledAction D Φ
          = D.gaugeFreePlus * Φ.strength.threePlus ^ 2
            + Real.sqrt 5 * D.gaugeFive * Φ.strength.threeMinus ^ 2
            + D.gaugeFive * Φ.strength.five ^ 2
            + D.ew.lam * (normSqSum Φ.scalar - D.ew.vev ^ 2) ^ 2
            + D.ew.vev ^ 2 / 4
                * (D.ew.gW ^ 2
                    * (Φ.direction.a1 ^ 2 + Φ.direction.a2 ^ 2)
                  + (D.ew.gY * Φ.direction.b
                      - D.ew.gW * Φ.direction.a3) ^ 2)
            + ((D.ew.yuk : ℂ) * epsPair Φ.matterDoublet Φ.scalar
                * Φ.matterSinglet).re)
    ∧ ((∀ s : SectorStrength,
          assembledAction D { baseTuple D with strength := s }
            = gaugeKineticEnergy
                (gNormalForm D.gaugeFreePlus D.gaugeFive) s)
      ∧ (∀ (h : Doublet) (A : EWParam),
          assembledAction D
              { baseTuple D with scalar := h, direction := A }
            = potential D.ew.lam D.ew.vev h
              + gaugeMassForm D.ew.gW D.ew.gY D.ew.vev A)
      ∧ (∀ (ψ : Doublet) (χ : ℂ),
          assembledAction D
              { baseTuple D with
                  matterDoublet := ψ, matterSinglet := χ }
            = (yukawaLine D.ew.yuk (vacuum D.ew.vev) ψ χ).re))
    ∧ (GAdInvariant (gNormalForm D.gaugeFreePlus D.gaugeFive)
      ∧ (FAdInvariant (gNormalForm D.gaugeFreePlus D.gaugeFive) →
          ∃! sc : ℝ,
            gNormalForm D.gaugeFreePlus D.gaugeFive
              = mirrorCommonRay sc)
      ∧ (∀ (u v w : ℂ), Complex.normSq u = 1 → u * v * w = 1 →
          ∀ Φ : FieldTuple,
            assembledAction D (scaleTuple u v w Φ)
              = assembledAction D Φ)
      ∧ (∃ (cH : ℤ) (dRow sRow : Fin 10),
          (cH = scalarChargeSix ∨ cH = -scalarChargeSix)
            ∧ OPH.ExteriorSelection.mem D.ew.base.selectionMask.val
                dRow = true
            ∧ OPH.ExteriorSelection.mem D.ew.base.selectionMask.val
                sRow = true
            ∧ cH + OPH.ExteriorSelection.charge dRow
                + OPH.ExteriorSelection.charge sRow = 0
            ∧ ∀ p : ℂ, Complex.normSq p = 1 → ∀ Φ : FieldTuple,
                assembledAction D
                    (scaleTuple (p ^ cH)
                      (p ^ OPH.ExteriorSelection.charge dRow)
                      (p ^ OPH.ExteriorSelection.charge sRow) Φ)
                  = assembledAction D Φ))
    ∧ ((∀ (ψ : Doublet) (χ : ℂ),
          assembledAction D
              { baseTuple D with
                  matterDoublet := ψ, matterSinglet := -Complex.I * χ }
            = (yukawaLine D.ew.yuk (vacuum D.ew.vev) ψ χ).im)
      ∧ assembledAction D (baseTuple D) = 0
      ∧ (∀ s : SectorStrength,
          0 ≤ gaugeKineticEnergy
              (gNormalForm D.gaugeFreePlus D.gaugeFive) s)) :=
  ⟨fun Φ => assembledAction_explicit D Φ,
    ⟨fun s => restrict_gauge D s,
      fun h A => restrict_scalar D h A,
      fun ψ χ => restrict_yukawa D ψ χ⟩,
    ⟨assembled_gauge_coefficients_G_invariant D.gaugeFreePlus
        D.gaugeFive,
      fun hF => assembled_gauge_mirror_under_F D.gaugeFreePlus
        D.gaugeFive hF,
      fun u v w hu hprod Φ =>
        assembled_abelian_invariance D u v w hu hprod Φ,
      committed_charge_abelian_invariance D⟩,
    ⟨fun ψ χ => yukawa_complex_recovery D ψ χ,
      assembledAction_base D,
      fun s => assembled_gauge_energy_nonneg D s⟩⟩

end OPH.AssembledActionComposition

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.AssembledActionComposition.restrict_gauge
#print axioms OPH.AssembledActionComposition.restrict_scalar
#print axioms OPH.AssembledActionComposition.restrict_yukawa
#print axioms OPH.AssembledActionComposition.restrict_yukawa_mass_coefficient
#print axioms OPH.AssembledActionComposition.assembledAction_base
#print axioms OPH.AssembledActionComposition.assembledAction_explicit
#print axioms OPH.AssembledActionComposition.assembled_gauge_energy_nonneg
#print axioms OPH.AssembledActionComposition.assembled_gauge_coefficients_G_invariant
#print axioms OPH.AssembledActionComposition.assembled_gauge_mirror_under_F
#print axioms OPH.AssembledActionComposition.gaugeKineticEnergy_mirror_ray
#print axioms OPH.AssembledActionComposition.assembled_abelian_invariance
#print axioms OPH.AssembledActionComposition.zpow_charge_invariance
#print axioms OPH.AssembledActionComposition.committed_charge_abelian_invariance
#print axioms OPH.AssembledActionComposition.yukawa_complex_recovery
#print axioms OPH.AssembledActionComposition.assembledActionPremiseData_nonvacuous
#print axioms OPH.AssembledActionComposition.committed_inhabitant_value
#print axioms OPH.AssembledActionComposition.assembledActionComposition_receipt
