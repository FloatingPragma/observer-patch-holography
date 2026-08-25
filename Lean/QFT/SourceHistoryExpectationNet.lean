import QFT.SourceHistoryThreeSlotLocalGNS

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# Finite-probability record conditioning on the source-history carrier

The eight retained length-three histories carry three binary slots.  On the
rational record layer `Fin 8 → ℚ`, this module constructs, for every slot
region `S ⊆ {0,1,2}`, a finite-probability conditioning family `E_S`
with respect to the empirical window law.  At a history `g`, `E_S f` averages
`f` over histories agreeing with `g` on every slot of `S`, weighted by the
committed window counts `![1149, 96, 2, 87, 3, 2, 3, 412]`.  Each `E_S` is linear over
the rationals, unital, pointwise positive, localized, idempotent, tower-compatible,
and preserves the empirical mean.  The two domain walls are fixed points of
their own regional conditioning maps, and every `E_S` preserves the exact
bond means `94/1754`, `103/1754`, and action mean `197/1754`.  The diagonal
embedding relates these rational record observables to represented diagonal
elements; it does not extend `E_S` to the full matrix algebra.

The matrix no-go has a fixed and narrow scope: there is no complex-linear
map from the full `M₈(ℂ)` carrier into the specific slot-`(0,1)` algebra
`algebra01` that simultaneously has range in `algebra01`, obeys the
`algebra01`-bimodule law, and preserves this committed state.  The slot-two
conditional weight is `383/415` in the `(s₀,s₁) = (0,0)` cell (`1149 : 96`)
but `2/89` in the `(0,1)` cell (`2 : 87`), whereas such a fixed-target
bimodule projection would force equality.

This theorem does not rule out ordinary AQFT inclusions with state
restriction, a different or enlarged target algebra, boundary memory or
ancillas, generalized expectations, or non-bimodule positive/completely
positive channels.  Those remain viable enriched routes.  Slots are
retained-window positions with no calibrated physical clock; no Lorentzian
localization, spectrum condition, continuum limit, fields, particles,
scattering, or detector readback is constructed here.
-/

namespace OPH.QFT.SourceHistoryExpectationNet

open Matrix
open OPH.QFT
open OPH.InformationProjection
open OPH.QFT.FiniteTwoSiteIsingField
open OPH.QFT.SourceHistoryGNSDynamics
open OPH.QFT.SourceHistoryThreeSlotLocalGNS

open scoped ComplexOrder InnerProductSpace

/-! ## The record layer: slot regions and agreement cells -/

/-- The state read out at slot `i` of the retained history `g`. -/
def slotState (i : Fin 3) : Fin 8 → Fin 2 :=
  if i = 0 then sourceState0 else if i = 1 then sourceState1 else sourceState2

theorem slotState_zero : slotState 0 = sourceState0 := rfl
theorem slotState_one : slotState 1 = sourceState1 := rfl
theorem slotState_two : slotState 2 = sourceState2 := rfl

/-- Two histories agree on every slot of the region `S`. -/
def agreesOn (S : Finset (Fin 3)) (g h : Fin 8) : Prop :=
  ∀ i ∈ S, slotState i g = slotState i h

instance (S : Finset (Fin 3)) (g h : Fin 8) : Decidable (agreesOn S g h) :=
  Finset.decidableDforallFinset

theorem agreesOn_refl (S : Finset (Fin 3)) (g : Fin 8) : agreesOn S g g :=
  fun _ _ => rfl

theorem agreesOn_symm {S : Finset (Fin 3)} {g h : Fin 8}
    (hgh : agreesOn S g h) : agreesOn S h g :=
  fun i hi => (hgh i hi).symm

theorem agreesOn_trans {S : Finset (Fin 3)} {g h k : Fin 8}
    (h₁ : agreesOn S g h) (h₂ : agreesOn S h k) : agreesOn S g k :=
  fun i hi => (h₁ i hi).trans (h₂ i hi)

theorem agreesOn_mono {S T : Finset (Fin 3)} (hST : S ⊆ T) {g h : Fin 8}
    (hgh : agreesOn T g h) : agreesOn S g h :=
  fun i hi => hgh i (hST hi)

theorem agreesOn_empty (g h : Fin 8) : agreesOn ∅ g h :=
  fun i hi => absurd hi (Finset.notMem_empty i)

/-- The three slots separate histories: agreement on the full region is
equality, by the committed bit encoding `g = 4 s₀ + 2 s₁ + s₂`. -/
theorem agreesOn_univ_iff {g h : Fin 8} :
    agreesOn Finset.univ g h ↔ g = h := by
  constructor
  · intro hgh
    have h0 : sourceState0 g = sourceState0 h := hgh 0 (Finset.mem_univ 0)
    have h1 : sourceState1 g = sourceState1 h := hgh 1 (Finset.mem_univ 1)
    have h2 : sourceState2 g = sourceState2 h := hgh 2 (Finset.mem_univ 2)
    apply Fin.ext
    rw [← sourceState_encoding g, ← sourceState_encoding h, h0, h1, h2]
  · rintro rfl
    exact agreesOn_refl _ _

/-! ## Cell weights from the committed window counts -/

/-- The committed window count of a history, as a rational weight. -/
def windowWeightQ (h : Fin 8) : ℚ := (sourceWindowCount h : ℚ)

theorem windowWeightQ_pos (h : Fin 8) : 0 < windowWeightQ h := by
  unfold windowWeightQ
  exact_mod_cast sourceWindowCount_pos h

/-- The indicator of the agreement cell of `g` for the region `S`. -/
def cellIndicator (S : Finset (Fin 3)) (g h : Fin 8) : ℚ :=
  if agreesOn S g h then 1 else 0

theorem cellIndicator_of_agrees {S : Finset (Fin 3)} {g h : Fin 8}
    (hgh : agreesOn S g h) : cellIndicator S g h = 1 := if_pos hgh

theorem cellIndicator_of_not_agrees {S : Finset (Fin 3)} {g h : Fin 8}
    (hgh : ¬ agreesOn S g h) : cellIndicator S g h = 0 := if_neg hgh

theorem cellIndicator_nonneg (S : Finset (Fin 3)) (g h : Fin 8) :
    0 ≤ cellIndicator S g h := by
  unfold cellIndicator
  split <;> norm_num

theorem cellIndicator_self (S : Finset (Fin 3)) (g : Fin 8) :
    cellIndicator S g g = 1 :=
  cellIndicator_of_agrees (agreesOn_refl S g)

theorem cellIndicator_congr_left {S : Finset (Fin 3)} {g g' : Fin 8}
    (hgg' : agreesOn S g g') (h : Fin 8) :
    cellIndicator S g h = cellIndicator S g' h := by
  by_cases hgh : agreesOn S g h
  · rw [cellIndicator_of_agrees hgh,
      cellIndicator_of_agrees (agreesOn_trans (agreesOn_symm hgg') hgh)]
  · rw [cellIndicator_of_not_agrees hgh,
      cellIndicator_of_not_agrees fun hc => hgh (agreesOn_trans hgg' hc)]

theorem cellIndicator_congr_right {S : Finset (Fin 3)} {h h' : Fin 8}
    (hhh' : agreesOn S h h') (g : Fin 8) :
    cellIndicator S g h = cellIndicator S g h' := by
  by_cases hgh : agreesOn S g h
  · rw [cellIndicator_of_agrees hgh,
      cellIndicator_of_agrees (agreesOn_trans hgh hhh')]
  · rw [cellIndicator_of_not_agrees hgh,
      cellIndicator_of_not_agrees
        fun hc => hgh (agreesOn_trans hc (agreesOn_symm hhh'))]

/-- The window weight of the agreement cell of `g` for the region `S`. -/
def cellWeight (S : Finset (Fin 3)) (g : Fin 8) : ℚ :=
  ∑ h, cellIndicator S g h * windowWeightQ h

theorem cellWeight_pos (S : Finset (Fin 3)) (g : Fin 8) :
    0 < cellWeight S g := by
  refine Finset.sum_pos'
    (fun h _ => mul_nonneg (cellIndicator_nonneg S g h)
      (le_of_lt (windowWeightQ_pos h)))
    ⟨g, Finset.mem_univ g, ?_⟩
  rw [cellIndicator_self, one_mul]
  exact windowWeightQ_pos g

theorem cellWeight_ne_zero (S : Finset (Fin 3)) (g : Fin 8) :
    cellWeight S g ≠ 0 :=
  ne_of_gt (cellWeight_pos S g)

theorem cellWeight_congr {S : Finset (Fin 3)} {g g' : Fin 8}
    (hgg' : agreesOn S g g') : cellWeight S g = cellWeight S g' :=
  Finset.sum_congr rfl fun h _ => by rw [cellIndicator_congr_left hgg' h]

/-! ## Finite-probability regional record conditioning -/

/-- Numerator of the regional record-conditioning map: the cell-restricted
window-weighted sum of `f`. -/
def regionalNum (S : Finset (Fin 3)) (f : Fin 8 → ℚ) (g : Fin 8) : ℚ :=
  ∑ h, cellIndicator S g h * windowWeightQ h * f h

theorem regionalNum_congr {S : Finset (Fin 3)} {g g' : Fin 8}
    (hgg' : agreesOn S g g') (f : Fin 8 → ℚ) :
    regionalNum S f g = regionalNum S f g' :=
  Finset.sum_congr rfl fun h _ => by rw [cellIndicator_congr_left hgg' h]

/-- Finite-probability conditioning of the rational record observable `f`
under the empirical window law, given the slots in `S`. -/
def regionalExpectation (S : Finset (Fin 3)) (f : Fin 8 → ℚ) (g : Fin 8) : ℚ :=
  regionalNum S f g / cellWeight S g

/-- `E_S f` is measurable for its own region: it is constant on each
agreement cell. -/
theorem regionalExpectation_localized {S : Finset (Fin 3)} {g g' : Fin 8}
    (f : Fin 8 → ℚ) (hgg' : agreesOn S g g') :
    regionalExpectation S f g = regionalExpectation S f g' := by
  unfold regionalExpectation
  rw [regionalNum_congr hgg' f, cellWeight_congr hgg']

/-- Additivity of the regional record-conditioning map. -/
theorem regionalExpectation_add (S : Finset (Fin 3)) (f f' : Fin 8 → ℚ) :
    regionalExpectation S (f + f') =
      regionalExpectation S f + regionalExpectation S f' := by
  funext g
  simp only [regionalExpectation, regionalNum, Pi.add_apply, mul_add]
  rw [Finset.sum_add_distrib, add_div]

/-- Rational homogeneity of the regional record-conditioning map. -/
theorem regionalExpectation_smul (S : Finset (Fin 3)) (c : ℚ) (f : Fin 8 → ℚ) :
    regionalExpectation S (c • f) = c • regionalExpectation S f := by
  funext g
  simp only [regionalExpectation, regionalNum, Pi.smul_apply, smul_eq_mul]
  rw [Finset.sum_congr rfl fun h _ =>
    show cellIndicator S g h * windowWeightQ h * (c * f h)
        = c * (cellIndicator S g h * windowWeightQ h * f h) by ring]
  rw [← Finset.mul_sum, mul_div_assoc]

/-- The regional record-conditioning map is unital. -/
theorem regionalExpectation_one (S : Finset (Fin 3)) :
    regionalExpectation S 1 = 1 := by
  funext g
  simp only [regionalExpectation, regionalNum, Pi.one_apply, mul_one]
  exact div_self (cellWeight_ne_zero S g)

/-- The regional record-conditioning map is positive on nonnegative record
observables. -/
theorem regionalExpectation_nonneg (S : Finset (Fin 3)) (f : Fin 8 → ℚ)
    (hf : ∀ g, 0 ≤ f g) (g : Fin 8) : 0 ≤ regionalExpectation S f g :=
  div_nonneg
    (Finset.sum_nonneg fun h _ =>
      mul_nonneg
        (mul_nonneg (cellIndicator_nonneg S g h)
          (le_of_lt (windowWeightQ_pos h)))
        (hf h))
    (le_of_lt (cellWeight_pos S g))

/-- A record observable that is constant on the agreement cells of `S` is
a fixed point of `E_S`. -/
theorem regionalExpectation_of_localized {S : Finset (Fin 3)} {f : Fin 8 → ℚ}
    (hf : ∀ g h, agreesOn S g h → f g = f h) :
    regionalExpectation S f = f := by
  funext g
  have hnum : regionalNum S f g = f g * cellWeight S g := by
    unfold regionalNum cellWeight
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun h _ => ?_
    by_cases hgh : agreesOn S g h
    · rw [hf g h hgh]
      ring
    · rw [cellIndicator_of_not_agrees hgh]
      ring
  unfold regionalExpectation
  rw [hnum, mul_div_cancel_right₀ _ (cellWeight_ne_zero S g)]

/-- On the full region the record-conditioning map is the identity: the
three slots separate the eight histories. -/
theorem regionalExpectation_univ (f : Fin 8 → ℚ) :
    regionalExpectation Finset.univ f = f :=
  regionalExpectation_of_localized fun g h hgh => by
    rw [agreesOn_univ_iff.mp hgh]

/-! ## The empirical mean and the tower law -/

/-- The empirical mean of a record observable under the committed window
law. -/
def recordMean (f : Fin 8 → ℚ) : ℚ := ∑ g, sourceTauEmpQ g * f g

theorem recordMean_eq_weighted (f : Fin 8 → ℚ) :
    recordMean f = (∑ h, windowWeightQ h * f h) / 1754 := by
  unfold recordMean
  rw [Finset.sum_div]
  refine Finset.sum_congr rfl fun h _ => ?_
  rw [sourceTauEmpQ_eq_counts h]
  unfold windowWeightQ
  ring

theorem cellWeight_empty (g : Fin 8) : cellWeight ∅ g = 1754 := by
  unfold cellWeight
  rw [Finset.sum_congr rfl fun h _ => by
    rw [cellIndicator_of_agrees (agreesOn_empty g h), one_mul]]
  rw [Fin.sum_univ_eight]
  simp only [windowWeightQ, sourceWindowCount_apply_0,
    sourceWindowCount_apply_1, sourceWindowCount_apply_2,
    sourceWindowCount_apply_3, sourceWindowCount_apply_4,
    sourceWindowCount_apply_5, sourceWindowCount_apply_6,
    sourceWindowCount_apply_7]
  norm_num

/-- The empty-region record-conditioning map is the constant empirical
mean. -/
theorem regionalExpectation_empty (f : Fin 8 → ℚ) :
    regionalExpectation ∅ f = fun _ => recordMean f := by
  funext g
  unfold regionalExpectation regionalNum
  rw [cellWeight_empty, recordMean_eq_weighted]
  congr 1
  exact Finset.sum_congr rfl fun h _ => by
    rw [cellIndicator_of_agrees (agreesOn_empty g h), one_mul]

theorem regionalNum_tower {S T : Finset (Fin 3)} (hST : S ⊆ T)
    (f : Fin 8 → ℚ) (g : Fin 8) :
    regionalNum S (regionalExpectation T f) g = regionalNum S f g := by
  unfold regionalNum
  simp only [regionalExpectation]
  calc
    ∑ h, cellIndicator S g h * windowWeightQ h *
        (regionalNum T f h / cellWeight T h)
      = ∑ h, ∑ k, cellIndicator S g h * windowWeightQ h *
          (cellIndicator T h k * windowWeightQ k * f k / cellWeight T h) := by
        refine Finset.sum_congr rfl fun h _ => ?_
        unfold regionalNum
        rw [Finset.sum_div, Finset.mul_sum]
    _ = ∑ k, ∑ h, cellIndicator S g h * windowWeightQ h *
          (cellIndicator T h k * windowWeightQ k * f k / cellWeight T h) :=
        Finset.sum_comm
    _ = ∑ k, cellIndicator S g k * windowWeightQ k * f k := by
        refine Finset.sum_congr rfl fun k _ => ?_
        have hkey : ∀ h : Fin 8,
            cellIndicator S g h * windowWeightQ h *
              (cellIndicator T h k * windowWeightQ k * f k / cellWeight T h)
            = cellIndicator S g k *
                (windowWeightQ k * f k / cellWeight T k) *
                (cellIndicator T k h * windowWeightQ h) := by
          intro h
          by_cases hTk : agreesOn T h k
          · rw [cellIndicator_of_agrees hTk,
              cellIndicator_of_agrees (agreesOn_symm hTk),
              cellWeight_congr hTk,
              cellIndicator_congr_right (agreesOn_mono hST hTk) g]
            ring
          · rw [cellIndicator_of_not_agrees hTk,
              cellIndicator_of_not_agrees
                fun hc => hTk (agreesOn_symm hc)]
            ring
        rw [Finset.sum_congr rfl fun h _ => hkey h, ← Finset.mul_sum]
        have hW : (∑ h, cellIndicator T k h * windowWeightQ h)
            = cellWeight T k := rfl
        have hW0 : cellWeight T k ≠ 0 := cellWeight_ne_zero T k
        rw [hW]
        field_simp

/-- **Tower law.**  Conditioning on a region and then on a subregion is
conditioning on the subregion. -/
theorem regionalExpectation_tower {S T : Finset (Fin 3)} (hST : S ⊆ T)
    (f : Fin 8 → ℚ) :
    regionalExpectation S (regionalExpectation T f) =
      regionalExpectation S f := by
  funext g
  show regionalNum S (regionalExpectation T f) g / cellWeight S g
      = regionalNum S f g / cellWeight S g
  rw [regionalNum_tower hST f g]

/-- Idempotence of every regional record-conditioning map. -/
theorem regionalExpectation_idem (S : Finset (Fin 3)) (f : Fin 8 → ℚ) :
    regionalExpectation S (regionalExpectation S f) =
      regionalExpectation S f :=
  regionalExpectation_tower (subset_refl S) f

/-- **Mean preservation.**  Every regional record-conditioning map preserves
the empirical mean: the tower law at the empty region. -/
theorem recordMean_regionalExpectation (S : Finset (Fin 3)) (f : Fin 8 → ℚ) :
    recordMean (regionalExpectation S f) = recordMean f := by
  have h := congrFun
    (regionalExpectation_tower (Finset.empty_subset S) f) 0
  rwa [congrFun (regionalExpectation_empty (regionalExpectation S f)) 0,
    congrFun (regionalExpectation_empty f) 0] at h

/-! ## The bundled net and its inhabitant -/

/-- A legacy-named bundle for the finite-probability record-conditioning
net on `Fin 8 → ℚ`: the maps are rational-linear, unital, pointwise
positive, localized, tower-compatible, and preserve the empirical mean. -/
structure StatePreservingRegionalNet where
  expect : Finset (Fin 3) → (Fin 8 → ℚ) → Fin 8 → ℚ
  map_add : ∀ S f f', expect S (f + f') = expect S f + expect S f'
  map_smul : ∀ (S) (c : ℚ) (f), expect S (c • f) = c • expect S f
  unital : ∀ S, expect S 1 = 1
  nonneg : ∀ S f, (∀ g, 0 ≤ f g) → ∀ g, 0 ≤ expect S f g
  localized : ∀ S f g g', agreesOn S g g' → expect S f g = expect S f g'
  tower : ∀ S T, S ⊆ T → ∀ f, expect S (expect T f) = expect S f
  statePreserving : ∀ S f, recordMean (expect S f) = recordMean f

/-- The committed counts inhabit this finite-probability conditioning
bundle: the rational record maps satisfy every listed clause jointly. -/
def sourceRegionalNet : StatePreservingRegionalNet where
  expect := regionalExpectation
  map_add := regionalExpectation_add
  map_smul := regionalExpectation_smul
  unital := regionalExpectation_one
  nonneg := regionalExpectation_nonneg
  localized := fun _ f _ _ hgg' => regionalExpectation_localized f hgg'
  tower := fun _ _ hST => regionalExpectation_tower hST
  statePreserving := recordMean_regionalExpectation

theorem sourceRegionalNet_expect :
    sourceRegionalNet.expect = regionalExpectation := rfl

/-! ## Intervals, walls, and exact conditional witnesses -/

/-- The adjacent slot interval `(0,1)`. -/
def interval01 : Finset (Fin 3) := {0, 1}

/-- The adjacent slot interval `(1,2)`. -/
def interval12 : Finset (Fin 3) := {1, 2}

/-- The early domain wall as a record observable. -/
def bond01RecordQ : Fin 8 → ℚ := fun g => (bond01Energy g : ℚ)

/-- The late domain wall as a record observable. -/
def bond12RecordQ : Fin 8 → ℚ := fun g => (bond12Energy g : ℚ)

/-- The repair action as a record observable. -/
def actionRecordQ : Fin 8 → ℚ := fun g => (sourceAction g : ℚ)

/-- The slot-two indicator `s₂ = 0` as a record observable. -/
def lateRecordQ : Fin 8 → ℚ := fun g => if sourceState2 g = 0 then 1 else 0

theorem lateRecordQ_apply_0 : lateRecordQ 0 = 1 := rfl
theorem lateRecordQ_apply_1 : lateRecordQ 1 = 0 := rfl
theorem lateRecordQ_apply_2 : lateRecordQ 2 = 1 := rfl
theorem lateRecordQ_apply_3 : lateRecordQ 3 = 0 := rfl

theorem bond01Energy_localized :
    ∀ g h : Fin 8, agreesOn interval01 g h →
      bond01Energy g = bond01Energy h := by decide

theorem bond12Energy_localized :
    ∀ g h : Fin 8, agreesOn interval12 g h →
      bond12Energy g = bond12Energy h := by decide

/-- The early wall is measurable for its own interval: `E_(0,1)` fixes
it. -/
theorem regionalExpectation_fixes_bond01 :
    regionalExpectation interval01 bond01RecordQ = bond01RecordQ :=
  regionalExpectation_of_localized fun g h hgh => by
    simp only [bond01RecordQ]
    exact_mod_cast bond01Energy_localized g h hgh

/-- The late wall is measurable for its own interval: `E_(1,2)` fixes
it. -/
theorem regionalExpectation_fixes_bond12 :
    regionalExpectation interval12 bond12RecordQ = bond12RecordQ :=
  regionalExpectation_of_localized fun g h hgh => by
    simp only [bond12RecordQ]
    exact_mod_cast bond12Energy_localized g h hgh

theorem recordMean_bond01 : recordMean bond01RecordQ = 94 / 1754 :=
  sourceTauEmpQ_bond01_mean

theorem recordMean_bond12 : recordMean bond12RecordQ = 103 / 1754 :=
  sourceTauEmpQ_bond12_mean

theorem recordMean_action : recordMean actionRecordQ = 197 / 1754 :=
  sourceTauEmpQ_meanAction

/-- The repair action is the sum of the two wall observables. -/
theorem sourceAction_eq_bond_sum :
    ∀ g : Fin 8, sourceAction g = bond01Energy g + bond12Energy g := by
  decide

theorem actionRecordQ_eq_bond_sum :
    actionRecordQ = bond01RecordQ + bond12RecordQ := by
  funext g
  simp only [actionRecordQ, bond01RecordQ, bond12RecordQ, Pi.add_apply]
  exact_mod_cast sourceAction_eq_bond_sum g

/-- Every regional record-conditioning map preserves the first-edge wall
mean `94/1754`. -/
theorem recordMean_regional_bond01 (S : Finset (Fin 3)) :
    recordMean (regionalExpectation S bond01RecordQ) = 94 / 1754 := by
  rw [recordMean_regionalExpectation, recordMean_bond01]

/-- Every regional record-conditioning map preserves the second-edge wall
mean `103/1754`. -/
theorem recordMean_regional_bond12 (S : Finset (Fin 3)) :
    recordMean (regionalExpectation S bond12RecordQ) = 103 / 1754 := by
  rw [recordMean_regionalExpectation, recordMean_bond12]

/-- Every regional record-conditioning map preserves the two-edge action-
incidence mean `197/1754`: the position pairing `94 + 103 = 197`. -/
theorem recordMean_regional_action (S : Finset (Fin 3)) :
    recordMean (regionalExpectation S actionRecordQ) = 197 / 1754 := by
  rw [recordMean_regionalExpectation, recordMean_action]

/-- The slot-two conditional weight in the `(s₀,s₁) = (0,0)` cell: the
committed counts `1149 : 96` give exactly `383/415`. -/
theorem regionalExpectation_late_cell00 :
    regionalExpectation interval01 lateRecordQ 0 = 383 / 415 := by
  show regionalNum interval01 lateRecordQ 0 / cellWeight interval01 0
      = 383 / 415
  unfold regionalNum cellWeight
  rw [Fin.sum_univ_eight, Fin.sum_univ_eight,
    cellIndicator_of_agrees (show agreesOn interval01 0 0 by decide),
    cellIndicator_of_agrees (show agreesOn interval01 0 1 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 0 2 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 0 3 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 0 4 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 0 5 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 0 6 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 0 7 by decide),
    lateRecordQ_apply_0, lateRecordQ_apply_1]
  norm_num [windowWeightQ, sourceWindowCount_apply_0,
    sourceWindowCount_apply_1]

/-- The slot-two conditional weight in the `(s₀,s₁) = (0,1)` cell: the
committed counts `2 : 87` give exactly `2/89`. -/
theorem regionalExpectation_late_cell01 :
    regionalExpectation interval01 lateRecordQ 2 = 2 / 89 := by
  show regionalNum interval01 lateRecordQ 2 / cellWeight interval01 2
      = 2 / 89
  unfold regionalNum cellWeight
  rw [Fin.sum_univ_eight, Fin.sum_univ_eight,
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 2 0 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 2 1 by decide),
    cellIndicator_of_agrees (show agreesOn interval01 2 2 by decide),
    cellIndicator_of_agrees (show agreesOn interval01 2 3 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 2 4 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 2 5 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 2 6 by decide),
    cellIndicator_of_not_agrees (show ¬ agreesOn interval01 2 7 by decide),
    lateRecordQ_apply_2, lateRecordQ_apply_3]
  norm_num [windowWeightQ, sourceWindowCount_apply_2,
    sourceWindowCount_apply_3]

/-- The two conditional weights differ: the slot-two distribution depends
on the interval cell. -/
theorem lateConditional_ne : (383 / 415 : ℚ) ≠ 2 / 89 := by norm_num

/-- The net genuinely averages: the slot-two indicator is not fixed by
`E_(0,1)`. -/
theorem regionalExpectation_late_not_fixed :
    regionalExpectation interval01 lateRecordQ ≠ lateRecordQ := by
  intro hcontra
  have h0 := congrFun hcontra 0
  rw [regionalExpectation_late_cell00, lateRecordQ_apply_0] at h0
  norm_num at h0

noncomputable section

/-! ## The diagonal embedding into the matrix carrier -/

/-- A record observable embedded as a diagonal matrix on the
source-counted carrier. -/
def recordDiagonal (f : Fin 8 → ℚ) : HistoryMatrix :=
  Matrix.diagonal fun g => (f g : ℂ)

/-- The diagonal embedding is multiplicative. -/
theorem recordDiagonal_mul (f f' : Fin 8 → ℚ) :
    recordDiagonal f * recordDiagonal f' = recordDiagonal (f * f') := by
  unfold recordDiagonal
  rw [Matrix.diagonal_mul_diagonal]
  congr 1
  funext g
  rw [Pi.mul_apply]
  push_cast
  ring

/-- The committed density trace of an embedded record observable is its
empirical mean. -/
theorem recordDiagonal_expectation (f : Fin 8 → ℚ) :
    EventAlgebra.expectation sourceHistoryDensity (recordDiagonal f)
      = ((recordMean f : ℚ) : ℂ) := by
  change (sourceHistoryDensity * recordDiagonal f).trace = _
  rw [sourceHistoryDensity, recordDiagonal, Matrix.diagonal_mul_diagonal,
    Matrix.trace_diagonal]
  rw [Finset.sum_congr rfl fun g _ =>
    show ((sourceTauEmpR g : ℝ) : ℂ) * ((f g : ℚ) : ℂ)
        = ((sourceTauEmpQ g * f g : ℚ) : ℂ) by
      simp only [sourceTauEmpR]
      push_cast
      ring]
  rw [← Rat.cast_sum]
  rfl

/-- State preservation of the net, expressed through the committed density
trace on the matrix carrier. -/
theorem recordDiagonal_regional_expectation (S : Finset (Fin 3))
    (f : Fin 8 → ℚ) :
    EventAlgebra.expectation sourceHistoryDensity
        (recordDiagonal (regionalExpectation S f))
      = EventAlgebra.expectation sourceHistoryDensity (recordDiagonal f) := by
  rw [recordDiagonal_expectation, recordDiagonal_expectation,
    recordMean_regionalExpectation]

theorem recordDiagonal_bond01 : recordDiagonal bond01RecordQ = bond01 := by
  rw [bond01_eq_diagonal]
  unfold recordDiagonal
  congr 1

theorem recordDiagonal_bond12 : recordDiagonal bond12RecordQ = bond12 := by
  rw [bond12_eq_diagonal]
  unfold recordDiagonal
  congr 1

theorem recordDiagonal_action :
    recordDiagonal actionRecordQ = sourceHistoryHamiltonian := by
  unfold recordDiagonal sourceHistoryHamiltonian
  congr 1

/-- Every regional conditioning of the two-edge action keeps the committed
mean `197/1754` under the density trace. -/
theorem regional_action_matrix_expectation (S : Finset (Fin 3)) :
    EventAlgebra.expectation sourceHistoryDensity
        (recordDiagonal (regionalExpectation S actionRecordQ))
      = (197 / 1754 : ℂ) := by
  rw [recordDiagonal_regional_expectation, recordDiagonal_action]
  exact sourceHistoryDensity_meanEnergy

/-- The interval expectation of the early wall is the wall itself on the
carrier, with its committed expectation `94/1754`. -/
theorem regional_bond01_matrix_expectation (S : Finset (Fin 3)) :
    EventAlgebra.expectation sourceHistoryDensity
        (recordDiagonal (regionalExpectation S bond01RecordQ))
      = (94 / 1754 : ℂ) := by
  rw [recordDiagonal_regional_expectation, recordDiagonal_bond01]
  exact sourceHistoryDensity_bond01_expectation

theorem regional_bond12_matrix_expectation (S : Finset (Fin 3)) :
    EventAlgebra.expectation sourceHistoryDensity
        (recordDiagonal (regionalExpectation S bond12RecordQ))
      = (103 / 1754 : ℂ) := by
  rw [recordDiagonal_regional_expectation, recordDiagonal_bond12]
  exact sourceHistoryDensity_bond12_expectation

theorem recordDiagonal_regional_bond01_fixed :
    recordDiagonal (regionalExpectation interval01 bond01RecordQ) = bond01 := by
  rw [regionalExpectation_fixes_bond01, recordDiagonal_bond01]

theorem recordDiagonal_regional_bond12_fixed :
    recordDiagonal (regionalExpectation interval12 bond12RecordQ) = bond12 := by
  rw [regionalExpectation_fixes_bond12, recordDiagonal_bond12]

/-! ## The matrix-level obstruction: setup -/

/-- The committed density trace of a matrix unit. -/
theorem expectation_single (a b : Fin 8) (z : ℂ) :
    EventAlgebra.expectation sourceHistoryDensity (Matrix.single a b z)
      = if b = a then ((sourceWindowCount a : ℂ) / 1754) * z else 0 := by
  change (sourceHistoryDensity * Matrix.single a b z).trace = _
  rw [Matrix.trace_mul_single, op_smul_eq_mul]
  by_cases hba : b = a
  · subst hba
    rw [if_pos rfl, sourceHistoryDensity_eq_windowCount]
  · rw [if_neg hba, sourceHistoryDensity,
      Matrix.diagonal_apply_ne _ hba, zero_mul]

/-- The `(0,1)` pair embedding on the exact source carrier. -/
def pair01 : Matrix (Spin × Spin) (Spin × Spin) ℂ →⋆ₐ[ℂ] HistoryMatrix :=
  (historyStageEquiv : ThreeSlotMatrix →⋆ₐ[ℂ] HistoryMatrix).comp
    (slotLeft Spin (α := Spin × Spin))

theorem slot0_eq_pair01 (A : SpinMatrix) :
    slot0 A = pair01 (slotLeft Spin (α := Spin) A) := rfl

theorem slot1_eq_pair01 (B : SpinMatrix) :
    slot1 B = pair01 (slotRight Spin (β := Spin) B) := rfl

theorem algebra01_le_pair01_range : algebra01 ≤ pair01.range := by
  show algebra0 ⊔ algebra1 ≤ pair01.range
  refine sup_le ?_ ?_
  · intro X hX
    obtain ⟨A, rfl⟩ := hX
    exact ⟨slotLeft Spin (α := Spin) A, rfl⟩
  · intro X hX
    obtain ⟨B, rfl⟩ := hX
    exact ⟨slotRight Spin (β := Spin) B, rfl⟩

theorem pair01_apply_two_two (B : Matrix (Spin × Spin) (Spin × Spin) ℂ) :
    pair01 B 2 2 = B ((0 : Spin), (1 : Spin)) ((0 : Spin), (1 : Spin)) := by
  norm_num [pair01, historyStageEquiv, historyEquiv, slotLeft,
    Matrix.reindex_apply, Matrix.kroneckerMap_apply, Matrix.one_apply,
    (by decide : sourceState0 2 = 0), (by decide : sourceState1 2 = 1),
    (by decide : sourceState2 2 = 0)]

theorem pair01_apply_three_three (B : Matrix (Spin × Spin) (Spin × Spin) ℂ) :
    pair01 B 3 3 = B ((0 : Spin), (1 : Spin)) ((0 : Spin), (1 : Spin)) := by
  norm_num [pair01, historyStageEquiv, historyEquiv, slotLeft,
    Matrix.reindex_apply, Matrix.kroneckerMap_apply, Matrix.one_apply,
    (by decide : sourceState0 3 = 0), (by decide : sourceState1 3 = 1),
    (by decide : sourceState2 3 = 1)]

/-- Histories `2` and `3` differ only at slot two, so every element of the
interval algebra has equal diagonal entries there. -/
theorem algebra01_wall_diagonal_pair {M : HistoryMatrix}
    (hM : M ∈ algebra01) : M 2 2 = M 3 3 := by
  obtain ⟨B, hB⟩ := algebra01_le_pair01_range hM
  have hB' : pair01 B = M := hB
  rw [← hB', pair01_apply_two_two, pair01_apply_three_three]

/-! ## The matrix-level obstruction: witnesses -/

/-- The slot-two projection onto `s₂ = 0`. -/
def latePoint : HistoryMatrix := slot2 (Matrix.single 0 0 1)

/-- The interval partial isometry moving the `(s₀,s₁)` cell `(0,1)` onto
`(0,0)`. -/
def cellShift : HistoryMatrix :=
  slot0 (Matrix.single 0 0 1) * slot1 (Matrix.single 0 1 1)

/-- The adjoint interval partial isometry. -/
def cellShiftStar : HistoryMatrix :=
  slot0 (Matrix.single 0 0 1) * slot1 (Matrix.single 1 0 1)

/-- The interval projection onto the `(s₀,s₁)` cell `(0,1)`. -/
def lateCellProj : HistoryMatrix :=
  slot0 (Matrix.single 0 0 1) * slot1 (Matrix.single 1 1 1)

theorem cellShift_mem_algebra01 : cellShift ∈ algebra01 :=
  mul_mem (algebra0_le_algebra01 ⟨Matrix.single 0 0 1, rfl⟩)
    (algebra1_le_algebra01 ⟨Matrix.single 0 1 1, rfl⟩)

theorem cellShiftStar_mem_algebra01 : cellShiftStar ∈ algebra01 :=
  mul_mem (algebra0_le_algebra01 ⟨Matrix.single 0 0 1, rfl⟩)
    (algebra1_le_algebra01 ⟨Matrix.single 1 0 1, rfl⟩)

theorem lateCellProj_mem_algebra01 : lateCellProj ∈ algebra01 :=
  mul_mem (algebra0_le_algebra01 ⟨Matrix.single 0 0 1, rfl⟩)
    (algebra1_le_algebra01 ⟨Matrix.single 1 1 1, rfl⟩)

theorem slot0_single00_explicit :
    slot0 (Matrix.single 0 0 1) =
      Matrix.single (0 : Fin 8) 0 (1 : ℂ) + Matrix.single 1 1 1 +
        Matrix.single 2 2 1 + Matrix.single 3 3 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [slot0, rawSlot0, historyStageEquiv, historyEquiv,
      sourceState0, sourceState1, sourceState2, slotLeft,
      Matrix.reindex_apply, Matrix.kroneckerMap_apply, Matrix.one_apply,
      Matrix.single_apply, Fin.ext_iff]

theorem slot1_single01_explicit :
    slot1 (Matrix.single 0 1 1) =
      Matrix.single (0 : Fin 8) 2 (1 : ℂ) + Matrix.single 1 3 1 +
        Matrix.single 4 6 1 + Matrix.single 5 7 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [slot1, rawSlot1, historyStageEquiv, historyEquiv,
      sourceState0, sourceState1, sourceState2, slotLeft, slotRight,
      Matrix.reindex_apply, Matrix.kroneckerMap_apply, Matrix.one_apply,
      Matrix.single_apply, Fin.ext_iff]

theorem slot1_single10_explicit :
    slot1 (Matrix.single 1 0 1) =
      Matrix.single (2 : Fin 8) 0 (1 : ℂ) + Matrix.single 3 1 1 +
        Matrix.single 6 4 1 + Matrix.single 7 5 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [slot1, rawSlot1, historyStageEquiv, historyEquiv,
      sourceState0, sourceState1, sourceState2, slotLeft, slotRight,
      Matrix.reindex_apply, Matrix.kroneckerMap_apply, Matrix.one_apply,
      Matrix.single_apply, Fin.ext_iff]

theorem slot1_single11_explicit :
    slot1 (Matrix.single 1 1 1) =
      Matrix.single (2 : Fin 8) 2 (1 : ℂ) + Matrix.single 3 3 1 +
        Matrix.single 6 6 1 + Matrix.single 7 7 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [slot1, rawSlot1, historyStageEquiv, historyEquiv,
      sourceState0, sourceState1, sourceState2, slotLeft, slotRight,
      Matrix.reindex_apply, Matrix.kroneckerMap_apply, Matrix.one_apply,
      Matrix.single_apply, Fin.ext_iff]

theorem slot2_single00_explicit :
    slot2 (Matrix.single 0 0 1) =
      Matrix.single (0 : Fin 8) 0 (1 : ℂ) + Matrix.single 2 2 1 +
        Matrix.single 4 4 1 + Matrix.single 6 6 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [slot2, rawSlot2, historyStageEquiv, historyEquiv,
      sourceState0, sourceState1, sourceState2, slotRight,
      Matrix.reindex_apply, Matrix.kroneckerMap_apply, Matrix.one_apply,
      Matrix.single_apply, Fin.ext_iff]

theorem latePoint_explicit :
    latePoint =
      Matrix.single (0 : Fin 8) 0 (1 : ℂ) + Matrix.single 2 2 1 +
        Matrix.single 4 4 1 + Matrix.single 6 6 1 :=
  slot2_single00_explicit

theorem cellShift_explicit :
    cellShift = Matrix.single (0 : Fin 8) 2 (1 : ℂ) + Matrix.single 1 3 1 := by
  rw [cellShift, slot0_single00_explicit, slot1_single01_explicit]
  simp [add_mul, mul_add]

theorem cellShiftStar_explicit :
    cellShiftStar =
      Matrix.single (2 : Fin 8) 0 (1 : ℂ) + Matrix.single 3 1 1 := by
  rw [cellShiftStar, slot0_single00_explicit, slot1_single10_explicit]
  simp [add_mul, mul_add]

theorem lateCellProj_explicit :
    lateCellProj =
      Matrix.single (2 : Fin 8) 2 (1 : ℂ) + Matrix.single 3 3 1 := by
  rw [lateCellProj, slot0_single00_explicit, slot1_single11_explicit]
  simp [add_mul, mul_add]

/-- Conjugating the slot-two projection by the interval shift produces the
single-history unit at `g = 0`. -/
theorem cellShift_latePoint_conj :
    cellShift * latePoint * cellShiftStar =
      Matrix.single (0 : Fin 8) 0 (1 : ℂ) := by
  rw [cellShift_explicit, latePoint_explicit, cellShiftStar_explicit]
  simp [add_mul, mul_add]

/-- Compressing the slot-two projection into the `(0,1)` cell produces the
single-history unit at `g = 2`. -/
theorem lateCellProj_latePoint_conj :
    lateCellProj * latePoint * lateCellProj =
      Matrix.single (2 : Fin 8) 2 (1 : ℂ) := by
  rw [lateCellProj_explicit, latePoint_explicit]
  simp [add_mul, mul_add]

theorem expectation_density_add (X Y : HistoryMatrix) :
    EventAlgebra.expectation sourceHistoryDensity (X + Y)
      = EventAlgebra.expectation sourceHistoryDensity X
        + EventAlgebra.expectation sourceHistoryDensity Y :=
  map_add (EventAlgebra.stateExpectationLinearMap sourceHistoryDensity) X Y

/-- The density trace of a shift-conjugated observable reads the diagonal
of the `(0,1)` cell with the `(0,0)`-cell weights `1149` and `96`. -/
theorem expectation_cellShift_conj (Y : HistoryMatrix) :
    EventAlgebra.expectation sourceHistoryDensity
        (cellShift * Y * cellShiftStar)
      = ((1149 : ℂ) / 1754) * Y 2 2 + ((96 : ℂ) / 1754) * Y 3 3 := by
  rw [cellShift_explicit, cellShiftStar_explicit]
  simp only [add_mul, mul_add]
  simp only [Matrix.single_mul_mul_single, one_mul, mul_one]
  simp only [expectation_density_add, expectation_single]
  norm_num [sourceWindowCount_apply_0, sourceWindowCount_apply_1,
    Fin.ext_iff]

/-- The density trace of a cell-compressed observable reads the same
diagonal with the `(0,1)`-cell weights `2` and `87`. -/
theorem expectation_lateCellProj_conj (Y : HistoryMatrix) :
    EventAlgebra.expectation sourceHistoryDensity
        (lateCellProj * Y * lateCellProj)
      = ((2 : ℂ) / 1754) * Y 2 2 + ((87 : ℂ) / 1754) * Y 3 3 := by
  rw [lateCellProj_explicit]
  simp only [add_mul, mul_add]
  simp only [Matrix.single_mul_mul_single, one_mul, mul_one]
  simp only [expectation_density_add, expectation_single]
  norm_num [sourceWindowCount_apply_2, sourceWindowCount_apply_3,
    Fin.ext_iff]

/-! ## The matrix-level obstruction -/

/-- The identity map satisfies the bimodule and state-preservation clauses
alone: the obstruction below is produced by compression into the interval
algebra, not by either clause separately. -/
theorem identity_bimodule_statePreserving :
    (∀ B ∈ algebra01, ∀ B' ∈ algebra01, ∀ X : HistoryMatrix,
        (LinearMap.id : HistoryMatrix →ₗ[ℂ] HistoryMatrix) (B * X * B')
          = B * (LinearMap.id : HistoryMatrix →ₗ[ℂ] HistoryMatrix) X * B') ∧
      (∀ X : HistoryMatrix,
        EventAlgebra.expectation sourceHistoryDensity
            ((LinearMap.id : HistoryMatrix →ₗ[ℂ] HistoryMatrix) X)
          = EventAlgebra.expectation sourceHistoryDensity X) :=
  ⟨fun _ _ _ _ _ => rfl, fun _ => rfl⟩

/-- **Fixed-target obstruction.**  No complex-linear map from the full
source-counted matrix carrier into the specific slot-`(0,1)` algebra
`algebra01` satisfies the `algebra01`-bimodule law and preserves the
committed source state.  Such a map would give the slot-two projection a
single conditional weight over `algebra01`, but the `(0,0)` cell demands
`1149/1245 = 383/415` while the `(0,1)` cell demands `2/89`.  This theorem
does not address other targets, enlarged carriers, ancillary memory,
generalized expectations, or non-bimodule positive/CP channels. -/
theorem no_statePreserving_bimodule_projection_onto_interval01 :
    ¬ ∃ E : HistoryMatrix →ₗ[ℂ] HistoryMatrix,
        (∀ X : HistoryMatrix, E X ∈ algebra01) ∧
        (∀ B ∈ algebra01, ∀ B' ∈ algebra01, ∀ X : HistoryMatrix,
          E (B * X * B') = B * E X * B') ∧
        (∀ X : HistoryMatrix,
          EventAlgebra.expectation sourceHistoryDensity (E X)
            = EventAlgebra.expectation sourceHistoryDensity X) := by
  rintro ⟨E, hrange, hbimod, hstate⟩
  have hdiag : E latePoint 2 2 = E latePoint 3 3 :=
    algebra01_wall_diagonal_pair (hrange latePoint)
  have h1 : EventAlgebra.expectation sourceHistoryDensity
      (cellShift * E latePoint * cellShiftStar) = (1149 / 1754 : ℂ) := by
    rw [← hbimod cellShift cellShift_mem_algebra01 cellShiftStar
      cellShiftStar_mem_algebra01 latePoint, hstate,
      cellShift_latePoint_conj, expectation_single]
    norm_num [sourceWindowCount_apply_0]
  have h2 : EventAlgebra.expectation sourceHistoryDensity
      (lateCellProj * E latePoint * lateCellProj) = (2 / 1754 : ℂ) := by
    rw [← hbimod lateCellProj lateCellProj_mem_algebra01 lateCellProj
      lateCellProj_mem_algebra01 latePoint, hstate,
      lateCellProj_latePoint_conj, expectation_single]
    norm_num [sourceWindowCount_apply_2]
  have e1 := (expectation_cellShift_conj (E latePoint)).symm.trans h1
  have e2 := (expectation_lateCellProj_conj (E latePoint)).symm.trans h2
  rw [hdiag] at e1 e2
  have hcontra : (99771 : ℂ) = 0 := by
    linear_combination (1245 * 1754 : ℂ) * e2 - (89 * 1754 : ℂ) * e1
  norm_num at hcontra

/-- **No state-preserving conditional expectation onto fixed
`algebra01`.**  The clause list (range in `algebra01`, identity on it,
`algebra01`-bimodule law, state preservation) is jointly unsatisfiable on
the committed density; the identity clause is not needed for the
contradiction.  This does not forbid a net of inclusions with ordinary
state restriction or an enriched expectation construction. -/
theorem no_statePreserving_conditional_expectation_onto_interval01 :
    ¬ ∃ E : HistoryMatrix →ₗ[ℂ] HistoryMatrix,
        (∀ X : HistoryMatrix, E X ∈ algebra01) ∧
        (∀ X ∈ algebra01, E X = X) ∧
        (∀ B ∈ algebra01, ∀ B' ∈ algebra01, ∀ X : HistoryMatrix,
          E (B * X * B') = B * E X * B') ∧
        (∀ X : HistoryMatrix,
          EventAlgebra.expectation sourceHistoryDensity (E X)
            = EventAlgebra.expectation sourceHistoryDensity X) := by
  rintro ⟨E, hrange, _hid, hbimod, hstate⟩
  exact no_statePreserving_bimodule_projection_onto_interval01
    ⟨E, hrange, hbimod, hstate⟩

/-! ## Selected-GNS composition -/

/-- The matrix coefficient of the selected cyclic vector is the committed
density trace, for every carrier observable. -/
theorem stage_matrixCoefficient (X : HistoryMatrix) :
    ⟪cyclicUnit, stageRepresentation X cyclicUnit⟫_ℂ =
      EventAlgebra.expectation sourceHistoryDensity X := by
  rw [cyclicUnit, stageRepresentation]
  change ⟪colimitGNSUnitClass selectedFunctional,
    colimitGNSRepresentation selectedFunctional
      (stageToCompletion sourceHistoryTower () X)
      (colimitGNSUnitClass selectedFunctional)⟫_ℂ = _
  rw [colimitGNSUnitClass_expectation, selectedColimitFunctional_stage]
  rfl

/-- Mean preservation of the rational record-conditioning map, expressed
through the selected GNS functional on represented diagonal elements. -/
theorem represented_regional_statePreservation (S : Finset (Fin 3))
    (f : Fin 8 → ℚ) :
    ⟪cyclicUnit,
        stageRepresentation
          (recordDiagonal (regionalExpectation S f)) cyclicUnit⟫_ℂ =
      ⟪cyclicUnit, stageRepresentation (recordDiagonal f) cyclicUnit⟫_ℂ := by
  rw [stage_matrixCoefficient, stage_matrixCoefficient,
    recordDiagonal_regional_expectation]

/-- The interval expectation of the early wall keeps its committed matrix
coefficient `94/1754` on the selected GNS space. -/
theorem represented_regional_bond01_coefficient :
    ⟪cyclicUnit,
        stageRepresentation
          (recordDiagonal (regionalExpectation interval01 bond01RecordQ))
          cyclicUnit⟫_ℂ = (94 / 1754 : ℂ) := by
  rw [recordDiagonal_regional_bond01_fixed, stage_matrixCoefficient]
  exact sourceHistoryDensity_bond01_expectation

/-- The interval expectation of the late wall keeps its committed matrix
coefficient `103/1754` on the selected GNS space. -/
theorem represented_regional_bond12_coefficient :
    ⟪cyclicUnit,
        stageRepresentation
          (recordDiagonal (regionalExpectation interval12 bond12RecordQ))
          cyclicUnit⟫_ℂ = (103 / 1754 : ℂ) := by
  rw [recordDiagonal_regional_bond12_fixed, stage_matrixCoefficient]
  exact sourceHistoryDensity_bond12_expectation

/-! ## Bundled receipt -/

/-- One bundled receipt for finite-probability record conditioning, its
committed wall means and exact conditional witnesses, together with the
narrow matrix no-go: the fixed target `algebra01` admits no
state-preserving `algebra01`-bimodule projection from this carrier. -/
theorem sourceHistoryExpectationNetAttachment :
    (∀ S T : Finset (Fin 3), S ⊆ T → ∀ f,
        regionalExpectation S (regionalExpectation T f)
          = regionalExpectation S f) ∧
      (∀ S f, recordMean (regionalExpectation S f) = recordMean f) ∧
      (∀ f, regionalExpectation Finset.univ f = f) ∧
      regionalExpectation interval01 bond01RecordQ = bond01RecordQ ∧
      regionalExpectation interval12 bond12RecordQ = bond12RecordQ ∧
      recordMean bond01RecordQ = 94 / 1754 ∧
      recordMean bond12RecordQ = 103 / 1754 ∧
      (∀ S, recordMean (regionalExpectation S actionRecordQ) = 197 / 1754) ∧
      regionalExpectation interval01 lateRecordQ 0 = 383 / 415 ∧
      regionalExpectation interval01 lateRecordQ 2 = 2 / 89 ∧
      ((383 / 415 : ℚ) ≠ 2 / 89) ∧
      ¬ ∃ E : HistoryMatrix →ₗ[ℂ] HistoryMatrix,
        (∀ X : HistoryMatrix, E X ∈ algebra01) ∧
        (∀ B ∈ algebra01, ∀ B' ∈ algebra01, ∀ X : HistoryMatrix,
          E (B * X * B') = B * E X * B') ∧
        (∀ X : HistoryMatrix,
          EventAlgebra.expectation sourceHistoryDensity (E X)
            = EventAlgebra.expectation sourceHistoryDensity X) :=
  ⟨fun _ _ hST => regionalExpectation_tower hST,
    recordMean_regionalExpectation,
    regionalExpectation_univ,
    regionalExpectation_fixes_bond01,
    regionalExpectation_fixes_bond12,
    recordMean_bond01,
    recordMean_bond12,
    recordMean_regional_action,
    regionalExpectation_late_cell00,
    regionalExpectation_late_cell01,
    lateConditional_ne,
    no_statePreserving_bimodule_projection_onto_interval01⟩

end

end OPH.QFT.SourceHistoryExpectationNet

-- Axiom audit: all declarations must remain on the standard Mathlib basis.
#print axioms OPH.QFT.SourceHistoryExpectationNet.agreesOn_univ_iff
#print axioms OPH.QFT.SourceHistoryExpectationNet.cellWeight_pos
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_localized
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_add
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_smul
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_one
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_nonneg
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_of_localized
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_univ
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_empty
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalNum_tower
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_tower
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_idem
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordMean_regionalExpectation
#print axioms OPH.QFT.SourceHistoryExpectationNet.sourceRegionalNet
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_fixes_bond01
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_fixes_bond12
#print axioms OPH.QFT.SourceHistoryExpectationNet.actionRecordQ_eq_bond_sum
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordMean_regional_action
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_late_cell00
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_late_cell01
#print axioms OPH.QFT.SourceHistoryExpectationNet.regionalExpectation_late_not_fixed
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordDiagonal_expectation
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordDiagonal_regional_expectation
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordDiagonal_bond01
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordDiagonal_bond12
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordDiagonal_action
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordDiagonal_regional_bond01_fixed
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordDiagonal_regional_bond12_fixed
#print axioms OPH.QFT.SourceHistoryExpectationNet.cellWeight_empty
#print axioms OPH.QFT.SourceHistoryExpectationNet.lateConditional_ne
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordMean_bond01
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordMean_bond12
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordMean_action
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordMean_regional_bond01
#print axioms OPH.QFT.SourceHistoryExpectationNet.recordMean_regional_bond12
#print axioms OPH.QFT.SourceHistoryExpectationNet.sourceAction_eq_bond_sum
#print axioms OPH.QFT.SourceHistoryExpectationNet.regional_action_matrix_expectation
#print axioms OPH.QFT.SourceHistoryExpectationNet.regional_bond01_matrix_expectation
#print axioms OPH.QFT.SourceHistoryExpectationNet.regional_bond12_matrix_expectation
#print axioms OPH.QFT.SourceHistoryExpectationNet.expectation_single
#print axioms OPH.QFT.SourceHistoryExpectationNet.algebra01_le_pair01_range
#print axioms OPH.QFT.SourceHistoryExpectationNet.algebra01_wall_diagonal_pair
#print axioms OPH.QFT.SourceHistoryExpectationNet.cellShift_latePoint_conj
#print axioms OPH.QFT.SourceHistoryExpectationNet.lateCellProj_latePoint_conj
#print axioms OPH.QFT.SourceHistoryExpectationNet.expectation_cellShift_conj
#print axioms OPH.QFT.SourceHistoryExpectationNet.expectation_lateCellProj_conj
#print axioms OPH.QFT.SourceHistoryExpectationNet.identity_bimodule_statePreserving
#print axioms OPH.QFT.SourceHistoryExpectationNet.no_statePreserving_bimodule_projection_onto_interval01
#print axioms OPH.QFT.SourceHistoryExpectationNet.no_statePreserving_conditional_expectation_onto_interval01
#print axioms OPH.QFT.SourceHistoryExpectationNet.stage_matrixCoefficient
#print axioms OPH.QFT.SourceHistoryExpectationNet.represented_regional_statePreservation
#print axioms OPH.QFT.SourceHistoryExpectationNet.represented_regional_bond01_coefficient
#print axioms OPH.QFT.SourceHistoryExpectationNet.represented_regional_bond12_coefficient
#print axioms OPH.QFT.SourceHistoryExpectationNet.sourceHistoryExpectationNetAttachment
