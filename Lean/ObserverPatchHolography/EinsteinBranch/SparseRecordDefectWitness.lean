import ObserverPatchHolography.EinsteinBranch.DarkSector

/-!
# Sparse-record defect witness for the dark-sector screening interface

The collar-scale saturation module proves that the recovery envelope never
forces a nonzero defect.  This module supplies the complementary positive
half: an explicit finite record model in which sparse recording provably
forces an order-one recovery defect, computed exactly, vanishing as the
record density saturates.

The model.  A source bit `A` uniform on `Fin 2` and a correlate `D`
perfectly correlated with it: `sourceCorrelateWeight` puts mass `1/2` on
each diagonal pair and none off the diagonal, so `D = A` almost surely and
predicting `D` from a record is predicting `A`.  A collar record of density
`p` carries the bit with probability `p` and an erasure with probability
`1 - p`: `recordWeight p` is the joint weight on `Fin 2 × Option (Fin 2)`
with mass `p/2` on `(a, some a)` and `(1 - p)/2` on `(a, none)`.  A
deterministic predictor `g : Option (Fin 2) → Fin 2` reads the record and
guesses the correlate; `predErr p g` is the probability of a wrong guess.

* The record law is a probability weight: nonnegative for `p ∈ [0, 1]`,
  total mass one, source marginal `1/2` matching the correlate law,
  erasure mass `1 - p`, carried mass `p` (`recordWeight_nonneg`,
  `recordWeight_sum`, `recordWeight_marginal`, `recordWeight_consistent`,
  `recordWeight_erasure`, `recordWeight_carried`).
* Optimality.  Any predictor that reads the carried bit errs exactly
  `(1 - p)/2` (`predErr_of_informed`), every deterministic predictor errs
  at least `(1 - p)/2` (`recordDefect_le_predErr`), and the model defect
  `recordDefect p = (1 - p)/2` is the least achievable error
  (`recordDefect_isLeast`).  A randomized predictor, a probability mixture
  of the eight deterministic ones, has error a convex combination of
  deterministic errors and obeys the same bound with the same least value
  (`recordDefect_le_randErr`, `recordDefect_isLeast_randomized`).
* Mechanism.  `recordDefect 0 = 1/2` and every predictor of the
  unrecorded cut errs at least `1/2` (`unrecordedCut_error_ge_half`);
  `recordDefect 1 = 0` and the zero error is attained
  (`saturatedRecord_error_eq_zero`); the defect is strictly decreasing and
  continuous in the density (`recordDefect_strictAnti`,
  `recordDefect_continuous`) and falls below `ε` exactly when
  `p > 1 - 2ε` (`recordDefect_lt_iff`), the mechanism-side counterpart of
  the envelope's settled-cut threshold.
* Composition with the dark-sector interface.  `recordRemainder` carries
  the model defect as an `AnomalyRemainder` with anomalous energy at the
  collar bound.  The sparse instance at `p = 0` with unit radius and unit
  collar constant has anomalous stress exactly `15 / (16 * π ^ 2) > 0`
  (`sparseRecordRemainder_stress_eq`, `sparseRecordRemainder_stress_pos`),
  and along any density family `p n → 1` the defect and the stress tend
  to zero (`recordRemainder_eta_and_stress_tendsto_zero`).  Within this
  record model, sparsely recorded cuts carry an order-one defect by
  computation and densely recorded cuts carry none.
* The identification of the model defect with the collar recovery defect
  is the `identification` field of `RecordedCollar`, a named modeling
  premise.  It forces the density bound `p ≤ 1`
  (`RecordedCollar.p_le_one`) and yields the exact stress bound and the
  saturation switch-off (`RecordedCollar.stress_abs_le`,
  `RecordedCollar.stress_eq_zero_of_saturated`,
  `RecordedCollar.eta_of_unrecorded`).

What is not proved here: the model is one record process, the erasure of a
single perfectly correlated bit, not the physical collar record process;
that the collar recovery defect of a physical cut equals the model defect
at some density is the modeling premise carried by the `identification`
field of `RecordedCollar`, not a consequence; and no physical value of the
record density is asserted.  Which density a galactic-scale settled cut
realizes is open.
-/

namespace OPH.EinsteinBranch

open Filter
open scoped Topology BigOperators

noncomputable section

/-! ## The finite record model -/

/-- Case dichotomy for `Fin 2`. -/
theorem finTwoCases : ∀ i : Fin 2, i = 0 ∨ i = 1 := by decide

/-- Joint law of the source bit `A` and its correlate `D`: uniform marginal,
perfect correlation.  All mass sits on the diagonal, so `D = A` almost
surely. -/
def sourceCorrelateWeight (a d : Fin 2) : ℝ := if a = d then 1 / 2 else 0

/-- Perfect correlation: off-diagonal pairs carry no mass. -/
theorem sourceCorrelateWeight_eq_zero_of_ne (a d : Fin 2) (h : a ≠ d) :
    sourceCorrelateWeight a d = 0 := if_neg h

/-- The source marginal of the correlate law is uniform. -/
theorem sourceCorrelateWeight_marginal (a : Fin 2) :
    ∑ d : Fin 2, sourceCorrelateWeight a d = 1 / 2 := by
  rcases finTwoCases a with rfl | rfl <;>
    simp [sourceCorrelateWeight]

/-- The correlate law is a probability weight. -/
theorem sourceCorrelateWeight_sum :
    ∑ a : Fin 2, ∑ d : Fin 2, sourceCorrelateWeight a d = 1 := by
  rw [Fin.sum_univ_two, sourceCorrelateWeight_marginal,
    sourceCorrelateWeight_marginal]
  norm_num

/-- Joint law of the source bit and the collar record at record density `p`:
with probability `p` the record carries the bit (`some a`, weight `p / 2`
per value), with probability `1 - p` it is an erasure (`none`, weight
`(1 - p) / 2` per value). -/
def recordWeight (p : ℝ) (a : Fin 2) : Option (Fin 2) → ℝ
  | some b => if a = b then p / 2 else 0
  | none => (1 - p) / 2

/-- For `p ∈ [0, 1]` the record law is nonnegative. -/
theorem recordWeight_nonneg (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1)
    (a : Fin 2) (r : Option (Fin 2)) : 0 ≤ recordWeight p a r := by
  cases r with
  | none =>
      show 0 ≤ (1 - p) / 2
      linarith
  | some b =>
      show 0 ≤ if a = b then p / 2 else 0
      by_cases h : a = b
      · rw [if_pos h]; linarith
      · rw [if_neg h]

/-- The source marginal of the record law is uniform at every density. -/
theorem recordWeight_marginal (p : ℝ) (a : Fin 2) :
    ∑ r : Option (Fin 2), recordWeight p a r = 1 / 2 := by
  rcases finTwoCases a with rfl | rfl <;>
    simp [recordWeight, Fintype.sum_option] <;> ring

/-- The algebraic weights sum to one for every real `p`. Together with
`recordWeight_nonneg`, they form a probability law only for `p ∈ [0,1]`. -/
theorem recordWeight_sum (p : ℝ) :
    ∑ a : Fin 2, ∑ r : Option (Fin 2), recordWeight p a r = 1 := by
  rw [Fin.sum_univ_two, recordWeight_marginal, recordWeight_marginal]
  norm_num

/-- The record law extends the correlate law: both have the same source
marginal. -/
theorem recordWeight_consistent (p : ℝ) (a : Fin 2) :
    ∑ r : Option (Fin 2), recordWeight p a r =
      ∑ d : Fin 2, sourceCorrelateWeight a d := by
  rw [recordWeight_marginal, sourceCorrelateWeight_marginal]

/-- The erasure mass is `1 - p`. -/
theorem recordWeight_erasure (p : ℝ) :
    ∑ a : Fin 2, recordWeight p a none = 1 - p := by
  rw [Fin.sum_univ_two]
  show (1 - p) / 2 + (1 - p) / 2 = 1 - p
  ring

/-- The carried mass is `p`. -/
theorem recordWeight_carried (p : ℝ) :
    ∑ a : Fin 2, recordWeight p a (some a) = p := by
  rw [Fin.sum_univ_two]
  show (if (0 : Fin 2) = 0 then p / 2 else 0) +
    (if (1 : Fin 2) = 1 then p / 2 else 0) = p
  norm_num

/-! ## Deterministic predictors and the model defect -/

/-- Error of the deterministic predictor `g` at record density `p`: the
probability that the guess `g` makes from the record differs from the
correlate, which equals the source bit by perfect correlation. -/
def predErr (p : ℝ) (g : Option (Fin 2) → Fin 2) : ℝ :=
  ∑ a : Fin 2, ∑ r : Option (Fin 2),
    recordWeight p a r * (if g r = a then 0 else 1)

/-- The model defect at record density `p`. -/
def recordDefect (p : ℝ) : ℝ := (1 - p) / 2

/-- For `p ≤ 1` the model defect is nonnegative. -/
theorem recordDefect_nonneg (p : ℝ) (hp : p ≤ 1) : 0 ≤ recordDefect p := by
  simp only [recordDefect]
  linarith

/-- The predictor that reads the carried bit and guesses `0` on erasure. -/
def informedPredictor : Option (Fin 2) → Fin 2
  | some a => a
  | none => 0

/-- Any predictor that reads the carried bit errs exactly `(1 - p) / 2`,
whatever it guesses on erasure. -/
theorem predErr_of_informed (p : ℝ) (g : Option (Fin 2) → Fin 2)
    (hg : ∀ a, g (some a) = a) : predErr p g = (1 - p) / 2 := by
  rcases finTwoCases (g none) with hn | hn <;>
    simp [predErr, recordWeight, Fintype.sum_option, Fin.sum_univ_two, hg,
      hn]

/-- The informed predictor errs exactly the model defect. -/
theorem predErr_informedPredictor (p : ℝ) :
    predErr p informedPredictor = recordDefect p :=
  predErr_of_informed p informedPredictor fun _ => rfl

/-- Algebraic optimality lower bound for `p ≥ 0`. Its interpretation as a
prediction-error bound additionally uses the probability range `p ≤ 1`. -/
theorem recordDefect_le_predErr (p : ℝ) (hp0 : 0 ≤ p)
    (g : Option (Fin 2) → Fin 2) : recordDefect p ≤ predErr p g := by
  rcases finTwoCases (g none) with hn | hn <;>
    rcases finTwoCases (g (some 0)) with h0 | h0 <;>
      rcases finTwoCases (g (some 1)) with h1 | h1 <;>
        simp [predErr, recordDefect, recordWeight, Fintype.sum_option,
          Fin.sum_univ_two, hn, h0, h1] <;>
        linarith

/-- Algebraically, the model defect is the least deterministic error
functional for `p ≥ 0`. It is a probability-error minimum when also `p ≤ 1`. -/
theorem recordDefect_isLeast (p : ℝ) (hp0 : 0 ≤ p) :
    IsLeast (Set.range fun g : Option (Fin 2) → Fin 2 => predErr p g)
      (recordDefect p) :=
  ⟨⟨informedPredictor, predErr_informedPredictor p⟩, by
    rintro e ⟨g, rfl⟩
    exact recordDefect_le_predErr p hp0 g⟩

/-! ## Randomized predictors -/

/-- A randomized predictor: a probability weight over the eight
deterministic predictors.  Its error is by definition a convex combination
of deterministic errors. -/
structure RandomizedPredictor where
  w : (Option (Fin 2) → Fin 2) → ℝ
  w_nonneg : ∀ g, 0 ≤ w g
  w_sum : ∑ g, w g = 1

/-- Error of a randomized predictor at record density `p`. -/
def randErr (p : ℝ) (R : RandomizedPredictor) : ℝ :=
  ∑ g, R.w g * predErr p g

/-- The point mass on the informed predictor. -/
def informedPoint : RandomizedPredictor where
  w := fun g => if g = informedPredictor then 1 else 0
  w_nonneg := fun g => by
    by_cases h : g = informedPredictor <;> simp [h]
  w_sum := by simp

/-- The point mass on the informed predictor errs exactly the informed
error. -/
theorem randErr_informedPoint (p : ℝ) :
    randErr p informedPoint = predErr p informedPredictor := by
  simp [randErr, informedPoint, ite_mul]

/-- Convexity reduction: every randomized predictor errs at least the model
defect, because its error is a convex combination of deterministic errors
each at least the defect. -/
theorem recordDefect_le_randErr (p : ℝ) (hp0 : 0 ≤ p)
    (R : RandomizedPredictor) : recordDefect p ≤ randErr p R := by
  calc recordDefect p
      = ∑ g : Option (Fin 2) → Fin 2, R.w g * recordDefect p := by
        rw [← Finset.sum_mul, R.w_sum, one_mul]
    _ ≤ ∑ g : Option (Fin 2) → Fin 2, R.w g * predErr p g :=
        Finset.sum_le_sum fun g _ =>
          mul_le_mul_of_nonneg_left (recordDefect_le_predErr p hp0 g)
            (R.w_nonneg g)

/-- Algebraically, the model defect is also the least randomized error
functional for `p ≥ 0`; the probability interpretation additionally needs
`p ≤ 1`. -/
theorem recordDefect_isLeast_randomized (p : ℝ) (hp0 : 0 ≤ p) :
    IsLeast (Set.range fun R : RandomizedPredictor => randErr p R)
      (recordDefect p) :=
  ⟨⟨informedPoint, by
      show randErr p informedPoint = recordDefect p
      rw [randErr_informedPoint, predErr_informedPredictor]⟩, by
    rintro e ⟨R, rfl⟩
    exact recordDefect_le_randErr p hp0 R⟩

/-! ## The mechanism -/

/-- Unrecorded cut: the model defect is one half. -/
theorem recordDefect_zero : recordDefect 0 = 1 / 2 := by
  simp only [recordDefect]
  norm_num

/-- Saturated record: the model defect is zero. -/
theorem recordDefect_one : recordDefect 1 = 0 := by
  simp only [recordDefect]
  norm_num

/-- On the unrecorded cut every deterministic predictor errs at least one
half: the order-one defect is forced, not chosen. -/
theorem unrecordedCut_error_ge_half (g : Option (Fin 2) → Fin 2) :
    1 / 2 ≤ predErr 0 g := by
  have h := recordDefect_le_predErr 0 le_rfl g
  rwa [recordDefect_zero] at h

/-- On the saturated record the informed predictor errs zero. -/
theorem saturatedRecord_error_eq_zero :
    predErr 1 informedPredictor = 0 := by
  rw [predErr_informedPredictor, recordDefect_one]

/-- The model defect is strictly decreasing in the record density. -/
theorem recordDefect_strictAnti : StrictAnti recordDefect := by
  intro x y hxy
  simp only [recordDefect]
  linarith

/-- The model defect is continuous in the record density. -/
theorem recordDefect_continuous : Continuous recordDefect := by
  unfold recordDefect
  fun_prop

/-- Exact density threshold: the model defect is below `ε` exactly when the
record density exceeds `1 - 2ε`.  This is the mechanism-side counterpart of
the envelope's settled-cut threshold. -/
theorem recordDefect_lt_iff (p ε : ℝ) :
    recordDefect p < ε ↔ 1 - 2 * ε < p := by
  simp only [recordDefect]
  constructor <;> intro h <;> linarith

/-! ## Composition with the dark-sector interface -/

/-- The anomaly remainder carried by a collar whose recovery defect is the
model defect at record density `p`, with collar radius `ℓ`, collar constant
`c`, and anomalous energy sitting exactly at the collar bound. -/
def recordRemainder (p ℓ c : ℝ) (hp : p ≤ 1) (hℓ : 0 < ℓ) (hc : 0 ≤ c) :
    AnomalyRemainder where
  ell := ℓ
  ell_pos := hℓ
  C := c
  C_nonneg := hc
  eta := recordDefect p
  eta_nonneg := recordDefect_nonneg p hp
  anomalousEnergy := c * recordDefect p
  bound :=
    le_of_eq (abs_of_nonneg (mul_nonneg hc (recordDefect_nonneg p hp)))

/-- The sparse instance: unrecorded cut, unit collar radius, unit collar
constant. -/
def sparseRecordRemainder : AnomalyRemainder :=
  recordRemainder 0 1 1 zero_le_one one_pos zero_le_one

/-- The sparse instance carries the order-one defect `1 / 2`. -/
theorem sparseRecordRemainder_eta : sparseRecordRemainder.eta = 1 / 2 :=
  recordDefect_zero

/-- The sparse instance carries the exact anomalous stress
`15 / (16 π ^ 2)`. -/
theorem sparseRecordRemainder_stress_eq :
    anomalousStress sparseRecordRemainder = 15 / (16 * Real.pi ^ 2) := by
  show 15 / (8 * Real.pi ^ 2 * (1 : ℝ) ^ 4) * (1 * recordDefect 0) =
    15 / (16 * Real.pi ^ 2)
  rw [recordDefect_zero]
  have h8 : (8 : ℝ) * Real.pi ^ 2 ≠ 0 := by positivity
  have h16 : (16 : ℝ) * Real.pi ^ 2 ≠ 0 := by positivity
  rw [one_pow, mul_one, one_mul, div_mul_eq_mul_div,
    div_eq_div_iff h8 h16]
  ring

/-- The sparse instance has strictly positive anomalous stress: within this
record model an unrecorded cut is a nonzero dark source. -/
theorem sparseRecordRemainder_stress_pos :
    0 < anomalousStress sparseRecordRemainder := by
  rw [sparseRecordRemainder_stress_eq]
  positivity

/-- Density saturation composed with the screening interface: along any
family of record densities tending to one, the model defect and the
anomalous stress of the carried remainders both tend to zero. -/
theorem recordRemainder_eta_and_stress_tendsto_zero (pseq : ℕ → ℝ)
    (ℓ c : ℝ) (hp : ∀ n, pseq n ≤ 1) (hℓ : 0 < ℓ) (hc : 0 ≤ c)
    (hlim : Tendsto pseq atTop (𝓝 1)) :
    Tendsto (fun n => (recordRemainder (pseq n) ℓ c (hp n) hℓ hc).eta)
        atTop (𝓝 0) ∧
      Tendsto
        (fun n => anomalousStress (recordRemainder (pseq n) ℓ c (hp n) hℓ hc))
        atTop (𝓝 0) := by
  have heta : Tendsto (fun n => recordDefect (pseq n)) atTop (𝓝 0) := by
    have h1 : Tendsto (fun n => (1 - pseq n) / 2) atTop
        (𝓝 ((1 - 1) / 2)) := (hlim.const_sub 1).div_const 2
    simpa [recordDefect] using h1
  refine ⟨heta, ?_⟩
  exact anomalousStress_tendsto_zero _ ℓ c (fun _ => rfl) (fun _ => le_rfl)
    heta

/-- A collar under the record-model premise: the `identification` field
asserts that the recovery defect of the carried remainder is the model
defect at record density `p`.  This identification is the modeling premise
of the sparse-record mechanism, not a consequence of the model. -/
structure RecordedCollar where
  p : ℝ
  p_nonneg : 0 ≤ p
  remainder : AnomalyRemainder
  identification : remainder.eta = recordDefect p

/-- The identification forces the density bound `p ≤ 1`: a recovery defect
is nonnegative, so the identified density cannot exceed one. -/
theorem RecordedCollar.p_le_one (X : RecordedCollar) : X.p ≤ 1 := by
  have h := X.remainder.eta_nonneg
  rw [X.identification] at h
  simp only [recordDefect] at h
  linarith

/-- Exact stress bound under the record-model premise: the anomalous stress
is bounded by the small-ball kernel factor times `(1 - p) / 2`. -/
theorem RecordedCollar.stress_abs_le (X : RecordedCollar) :
    |anomalousStress X.remainder| ≤
      15 * X.remainder.C * ((1 - X.p) / 2) /
        (8 * Real.pi ^ 2 * X.remainder.ell ^ 4) := by
  have h := anomalousStress_abs_le X.remainder
  rwa [X.identification, recordDefect] at h

/-- Saturation switch-off under the record-model premise: at density one
the anomalous stress vanishes. -/
theorem RecordedCollar.stress_eq_zero_of_saturated (X : RecordedCollar)
    (h : X.p = 1) : anomalousStress X.remainder = 0 := by
  apply anomalousStress_eq_zero_of_recovered
  rw [X.identification, h, recordDefect_one]

/-- Unrecorded cut under the record-model premise: the recovery defect is
exactly one half. -/
theorem RecordedCollar.eta_of_unrecorded (X : RecordedCollar)
    (h : X.p = 0) : X.remainder.eta = 1 / 2 := by
  rw [X.identification, h, recordDefect_zero]

/-- The sparse collar: the unrecorded cut carried as a recorded collar.
The hypotheses of `RecordedCollar` are jointly satisfiable with an
order-one defect and strictly positive stress. -/
def sparseRecordedCollar : RecordedCollar where
  p := 0
  p_nonneg := le_rfl
  remainder := sparseRecordRemainder
  identification := rfl

/-- The sparse collar has strictly positive anomalous stress. -/
theorem sparseRecordedCollar_stress_pos :
    0 < anomalousStress sparseRecordedCollar.remainder :=
  sparseRecordRemainder_stress_pos

/-! ## Per-theorem axiom audit -/

#print axioms finTwoCases
#print axioms sourceCorrelateWeight_eq_zero_of_ne
#print axioms sourceCorrelateWeight_marginal
#print axioms sourceCorrelateWeight_sum
#print axioms recordWeight_nonneg
#print axioms recordWeight_marginal
#print axioms recordWeight_sum
#print axioms recordWeight_consistent
#print axioms recordWeight_erasure
#print axioms recordWeight_carried
#print axioms predErr_of_informed
#print axioms predErr_informedPredictor
#print axioms recordDefect_nonneg
#print axioms recordDefect_le_predErr
#print axioms recordDefect_isLeast
#print axioms randErr_informedPoint
#print axioms recordDefect_le_randErr
#print axioms recordDefect_isLeast_randomized
#print axioms recordDefect_zero
#print axioms recordDefect_one
#print axioms unrecordedCut_error_ge_half
#print axioms saturatedRecord_error_eq_zero
#print axioms recordDefect_strictAnti
#print axioms recordDefect_continuous
#print axioms recordDefect_lt_iff
#print axioms sparseRecordRemainder_eta
#print axioms sparseRecordRemainder_stress_eq
#print axioms sparseRecordRemainder_stress_pos
#print axioms recordRemainder_eta_and_stress_tendsto_zero
#print axioms RecordedCollar.p_le_one
#print axioms RecordedCollar.stress_abs_le
#print axioms RecordedCollar.stress_eq_zero_of_saturated
#print axioms RecordedCollar.eta_of_unrecorded
#print axioms sparseRecordedCollar_stress_pos

end

end OPH.EinsteinBranch
