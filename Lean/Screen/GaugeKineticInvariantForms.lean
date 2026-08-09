import Mathlib

/-!
# Exact carrier/ad-invariant gauge-form controls

This file is the small proof-kernel companion to
`code/e9_kinetic/gauge_kinetic_invariant_forms.certificate.json`.  The exact
finite replay supplies the constraint rows.  Here we prove what those rows
mean: the generic `F` and `G` carrier form spaces fall from three parameters
to two, the plane control remains two-to-two, and imposing **both** mirror
constraints on one form leaves one ray.

The last premise is deliberately explicit.  Nothing here says that the OPH
source chooses simultaneous mirror invariance, a bracket, an overall scale,
or a relative product-factor coefficient.
-/

namespace OPH.GaugeKineticInvariantForms

/-- Projector weights on the `3+`, `3-`, and `5` carrier sectors. -/
@[ext]
structure ThreeSectorWeights where
  threePlus : ℝ
  threeMinus : ℝ
  five : ℝ

/-- The single exact constraint replayed for the supplied generic `F` bracket. -/
def FAdInvariant (w : ThreeSectorWeights) : Prop :=
  w.threePlus = Real.sqrt 5 * w.five

/-- The single exact constraint replayed for the supplied generic `G` bracket. -/
def GAdInvariant (w : ThreeSectorWeights) : Prop :=
  w.threeMinus = Real.sqrt 5 * w.five

/-- A two-parameter normal form for the `F` ad-invariant weights. -/
noncomputable def fNormalForm (freeMinus five : ℝ) : ThreeSectorWeights :=
  ⟨Real.sqrt 5 * five, freeMinus, five⟩

/-- A two-parameter normal form for the `G` ad-invariant weights. -/
noncomputable def gNormalForm (freePlus five : ℝ) : ThreeSectorWeights :=
  ⟨freePlus, Real.sqrt 5 * five, five⟩

/-- The replayed `F` row is exactly `w₃₊ - √5 w₅ = 0`. -/
theorem f_constraint_row_iff (w : ThreeSectorWeights) :
    FAdInvariant w ↔ w.threePlus - Real.sqrt 5 * w.five = 0 := by
  constructor <;> intro h
  · exact sub_eq_zero.mpr h
  · exact sub_eq_zero.mp h

/-- The replayed `G` row is exactly `w₃₋ - √5 w₅ = 0`. -/
theorem g_constraint_row_iff (w : ThreeSectorWeights) :
    GAdInvariant w ↔ w.threeMinus - Real.sqrt 5 * w.five = 0 := by
  constructor <;> intro h
  · exact sub_eq_zero.mpr h
  · exact sub_eq_zero.mp h

/-- Exact `F` carrier-three to ad-invariant-two result, with no hidden quotient. -/
theorem f_exact_two_parameter (w : ThreeSectorWeights) :
    FAdInvariant w ↔
      ∃! p : ℝ × ℝ, w = fNormalForm p.1 p.2 := by
  constructor
  · intro h
    refine ⟨(w.threeMinus, w.five), ?_, ?_⟩
    · ext <;> simp [fNormalForm, FAdInvariant] at h ⊢
      exact h
    · intro p hp
      apply Prod.ext
      · exact (congrArg ThreeSectorWeights.threeMinus hp).symm
      · exact (congrArg ThreeSectorWeights.five hp).symm
  · rintro ⟨p, hp, -⟩
    rw [hp]
    simp [FAdInvariant, fNormalForm]

/-- Exact `G` carrier-three to ad-invariant-two result, with no hidden quotient. -/
theorem g_exact_two_parameter (w : ThreeSectorWeights) :
    GAdInvariant w ↔
      ∃! p : ℝ × ℝ, w = gNormalForm p.1 p.2 := by
  constructor
  · intro h
    refine ⟨(w.threePlus, w.five), ?_, ?_⟩
    · ext <;> simp [gNormalForm, GAdInvariant] at h ⊢
      exact h
    · intro p hp
      apply Prod.ext
      · exact (congrArg ThreeSectorWeights.threePlus hp).symm
      · exact (congrArg ThreeSectorWeights.five hp).symm
  · rintro ⟨p, hp, -⟩
    rw [hp]
    simp [GAdInvariant, gNormalForm]

/-- The plane family has only the two triplet weights. -/
@[ext]
structure PlaneWeights where
  threePlus : ℝ
  threeMinus : ℝ

/-- There is no additional ad-invariance row in the exact plane control. -/
def PAdInvariant (_w : PlaneWeights) : Prop := True

def pNormalForm (p : ℝ × ℝ) : PlaneWeights := ⟨p.1, p.2⟩

/-- Exact plane carrier-two to ad-invariant-two control. -/
theorem p_exact_two_parameter (w : PlaneWeights) :
    PAdInvariant w ↔ ∃! p : ℝ × ℝ, w = pNormalForm p := by
  constructor
  · intro _
    refine ⟨(w.threePlus, w.threeMinus), ?_, ?_⟩
    · rfl
    · intro p hp
      apply Prod.ext
      · exact (congrArg PlaneWeights.threePlus hp).symm
      · exact (congrArg PlaneWeights.threeMinus hp).symm
  · intro _
    trivial

/-- The ray left by the *extra* simultaneous `F`/`G` premise. -/
noncomputable def mirrorCommonRay (scale : ℝ) : ThreeSectorWeights :=
  ⟨Real.sqrt 5 * scale, Real.sqrt 5 * scale, scale⟩

/--
If one explicitly requires the same form to be ad-invariant for both supplied
mirror brackets, it is uniquely parametrized by one scale.  This theorem does
not derive that simultaneous premise from an OPH source.
-/
theorem mirror_common_extra_premise_one_ray (w : ThreeSectorWeights) :
    (FAdInvariant w ∧ GAdInvariant w) ↔
      ∃! scale : ℝ, w = mirrorCommonRay scale := by
  constructor
  · rintro ⟨hF, hG⟩
    refine ⟨w.five, ?_, ?_⟩
    · ext <;> simp [mirrorCommonRay, FAdInvariant, GAdInvariant] at hF hG ⊢
      · exact hF
      · exact hG
    · intro scale hs
      exact (congrArg ThreeSectorWeights.five hs).symm
  · rintro ⟨scale, hs, -⟩
    rw [hs]
    exact ⟨rfl, rfl⟩

#print axioms f_exact_two_parameter
#print axioms g_exact_two_parameter
#print axioms p_exact_two_parameter
#print axioms mirror_common_extra_premise_one_ray

end OPH.GaugeKineticInvariantForms
