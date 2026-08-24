import ObserverPatchHolography.EinsteinBranch.DarkSector

/-!
# Finite cap-generator split and localization receipt

The dark-matter paper's source section uses the local split of the cap
modular generator, `K_C = 2 pi B_C + K_C^anom + const`, with the anomalous
part carried on the overlap collar.  This module proves the finite
counterpart of that split on a finite carrier of microstates.

* A `CapState` is a finite carrier labelled bulk or collar, with strictly
  positive unnormalized weights and a declared boost datum vanishing on the
  collar.  The modular generator is `K x = -log (w x)`.
* Existence: with the bulk-mean normalization clause (the bulk-weighted mean
  of the anomalous part vanishes), the split
  `K = 2 pi B + K_anom + c 1` holds exactly (`modularGenerator_split`,
  `anomalousPart_bulk_normalized`), and any split satisfying the clause
  coincides with the canonical one (`split_unique`, `split_eq_canonical`).
  Without the clause the split is non-unique
  (`split_nonunique_without_normalization`), so the clause is load-bearing.
* Quotient invariance: under the state rescaling `w -> exp (-t) * w`, which
  shifts `K` by `t`, the constant moves to `c + t` and the anomalous part is
  unchanged (`normalizingConstant_rescale`, `anomalousPart_rescale`); the
  collar expectation of the anomalous part is invariant as well
  (`collarAnomalousEnergy_rescale`).
* Localization: if the bulk weights are exactly Gibbs for `2 pi B`, the
  anomalous part vanishes on the bulk, so it is carried entirely on the
  collar (`anomalousPart_eq_zero_on_bulk`,
  `anomalousPart_supported_on_collar`), and the state expectation of the
  anomalous part equals its collar expectation
  (`anomalousExpectation_eq_collar`).
* Interface: a defect hypothesis `|<K_anom>_collar| <= C * eta` builds the
  `AnomalyRemainder` consumed by the dark-sector module, with the remainder
  energy equal to the collar expectation of the split by definitional
  unfolding (`collarRemainder`, `collarRemainder_energy_eq`,
  `collarRemainder_stress_abs_le`).
* A fully explicit inhabitant on `Fin 3` (two bulk states of weight `1`, one
  collar state of weight `exp (-1)`, zero boost) carries the exact numbers:
  normalizing constant `0`, anomalous part `(0, 0, 1)`, collar expectation
  `exp (-1) / (2 + exp (-1))`, and a remainder instance with `C = 1/2`,
  `eta = 1`.

What is not proved here: no continuum cap or type-III modular theory enters;
the boost datum is declared, not derived from a physical modular flow; the
bulk-mean normalization clause, the collar support of the boost, and the
Gibbs-on-bulk condition are model premises; no claim is made that the
physical cap generator satisfies these clauses; the collar constant and the
recovery defect in the interface step are inputs.  The receipt is the finite
counterpart of the paper's split, stated at exactly that strength.
-/

namespace OPH.EinsteinBranch

open scoped BigOperators

universe u

noncomputable section

/-- Finite cap state: microstates labelled bulk or collar, strictly positive
unnormalized weights, and a declared boost datum vanishing on the collar.
The bulk sector is required to be inhabited so the bulk-mean normalization
below is well posed. -/
structure CapState (α : Type u) [Fintype α] where
  isBulk : α → Bool
  bulk_nonempty : ∃ x, isBulk x = true
  weight : α → ℝ
  weight_pos : ∀ x, 0 < weight x
  boost : α → ℝ
  boost_collar_zero : ∀ x, isBulk x = false → boost x = 0

variable {α : Type u} [Fintype α]

/-- Bulk microstates. -/
def bulkSet (D : CapState α) : Finset α :=
  Finset.univ.filter (fun x => D.isBulk x = true)

/-- Collar microstates. -/
def collarSet (D : CapState α) : Finset α :=
  Finset.univ.filter (fun x => D.isBulk x = false)

/-- Modular generator of the state, `K x = -log (w x)`. -/
def modularGenerator (D : CapState α) : α → ℝ :=
  fun x => -Real.log (D.weight x)

/-- Total weight of the bulk sector. -/
def bulkMass (D : CapState α) : ℝ := ∑ x ∈ bulkSet D, D.weight x

/-- Total weight of the carrier. -/
def totalMass (D : CapState α) : ℝ := ∑ x, D.weight x

/-- Bulk-weighted sum of `K - 2 pi B`, the numerator of the normalizing
constant. -/
def bulkBoostGap (D : CapState α) : ℝ :=
  ∑ x ∈ bulkSet D, D.weight x *
    (modularGenerator D x - 2 * Real.pi * D.boost x)

/-- The normalizing constant `c`: the bulk-weighted mean of `K - 2 pi B`. -/
def normalizingConstant (D : CapState α) : ℝ :=
  bulkBoostGap D / bulkMass D

/-- The anomalous part `K_anom = K - 2 pi B - c 1`. -/
def anomalousPart (D : CapState α) : α → ℝ :=
  fun x => modularGenerator D x - 2 * Real.pi * D.boost x -
    normalizingConstant D

/-- The normalization clause: the bulk-weighted mean of `A` vanishes. -/
def BulkNormalized (D : CapState α) (A : α → ℝ) : Prop :=
  ∑ x ∈ bulkSet D, D.weight x * A x = 0

theorem bulkMass_pos (D : CapState α) : 0 < bulkMass D := by
  obtain ⟨x0, hx0⟩ := D.bulk_nonempty
  refine Finset.sum_pos (fun x _ => D.weight_pos x) ⟨x0, ?_⟩
  simp [bulkSet, hx0]

theorem totalMass_pos (D : CapState α) : 0 < totalMass D := by
  obtain ⟨x0, _⟩ := D.bulk_nonempty
  exact Finset.sum_pos (fun x _ => D.weight_pos x) ⟨x0, Finset.mem_univ x0⟩

/-! ## Existence of the split -/

/-- The split exists: `K = 2 pi B + K_anom + c 1` with `c` the explicit
normalizing constant. -/
theorem modularGenerator_split (D : CapState α) (x : α) :
    modularGenerator D x =
      2 * Real.pi * D.boost x + anomalousPart D x + normalizingConstant D := by
  unfold anomalousPart
  ring

/-- The canonical anomalous part satisfies the normalization clause. -/
theorem anomalousPart_bulk_normalized (D : CapState α) :
    BulkNormalized D (anomalousPart D) := by
  have hb : bulkMass D ≠ 0 := ne_of_gt (bulkMass_pos D)
  unfold BulkNormalized
  have h : ∀ x ∈ bulkSet D,
      D.weight x * anomalousPart D x
        = D.weight x * (modularGenerator D x - 2 * Real.pi * D.boost x)
          - D.weight x * normalizingConstant D := by
    intro x _
    unfold anomalousPart
    ring
  rw [Finset.sum_congr rfl h, Finset.sum_sub_distrib, ← Finset.sum_mul]
  show bulkBoostGap D - bulkMass D * normalizingConstant D = 0
  unfold normalizingConstant
  field_simp
  ring

/-! ## Quotient-invariant uniqueness -/

/-- Uniqueness: two splits of the same generator over the same boost datum
that both satisfy the normalization clause coincide, anomalous part and
constant alike. -/
theorem split_unique (D : CapState α) (A A' : α → ℝ) (c c' : ℝ)
    (hA : ∀ x, modularGenerator D x = 2 * Real.pi * D.boost x + A x + c)
    (hA' : ∀ x, modularGenerator D x = 2 * Real.pi * D.boost x + A' x + c')
    (hN : BulkNormalized D A) (hN' : BulkNormalized D A') :
    A = A' ∧ c = c' := by
  have hdiff : ∀ x, A x = A' x + (c' - c) := by
    intro x
    have h1 := hA x
    have h2 := hA' x
    linarith
  have hsum : ∑ x ∈ bulkSet D, D.weight x * A x
      = ∑ x ∈ bulkSet D, D.weight x * A' x + bulkMass D * (c' - c) := by
    have h : ∀ x ∈ bulkSet D, D.weight x * A x
        = D.weight x * A' x + D.weight x * (c' - c) := by
      intro x _
      rw [hdiff x]
      ring
    rw [Finset.sum_congr rfl h, Finset.sum_add_distrib, ← Finset.sum_mul]
    rfl
  rw [hN, hN'] at hsum
  have hb := bulkMass_pos D
  have hzero : bulkMass D * (c' - c) = 0 := by linarith
  have hcc : c' - c = 0 := by
    rcases mul_eq_zero.mp hzero with h | h
    · exact absurd h (ne_of_gt hb)
    · exact h
  refine ⟨funext fun x => ?_, by linarith⟩
  rw [hdiff x, hcc, add_zero]

/-- Any split satisfying the normalization clause is the canonical one. -/
theorem split_eq_canonical (D : CapState α) (A : α → ℝ) (c : ℝ)
    (hA : ∀ x, modularGenerator D x = 2 * Real.pi * D.boost x + A x + c)
    (hN : BulkNormalized D A) :
    A = anomalousPart D ∧ c = normalizingConstant D :=
  split_unique D A (anomalousPart D) c (normalizingConstant D) hA
    (modularGenerator_split D) hN (anomalousPart_bulk_normalized D)

/-- The normalization clause is load-bearing: without it, any nonzero
constant can be shuffled between the anomalous part and the constant term,
so the split is non-unique. -/
theorem split_nonunique_without_normalization (D : CapState α) (t : ℝ)
    (ht : t ≠ 0) :
    ∃ A' c', (∀ x, modularGenerator D x
        = 2 * Real.pi * D.boost x + A' x + c') ∧ A' ≠ anomalousPart D := by
  refine ⟨fun x => anomalousPart D x + t, normalizingConstant D - t,
    fun x => ?_, fun hcontra => ?_⟩
  · have h := modularGenerator_split D x
    show modularGenerator D x
      = 2 * Real.pi * D.boost x + (anomalousPart D x + t)
        + (normalizingConstant D - t)
    linarith
  · obtain ⟨x0, _⟩ := D.bulk_nonempty
    have h2 : anomalousPart D x0 + t = anomalousPart D x0 :=
      congrFun hcontra x0
    exact ht (by linarith)

/-! ## Quotient invariance under state rescaling -/

/-- State rescaling `w -> exp (-t) * w`.  It shifts the modular generator by
`t` and realizes the additive-constant quotient of the split. -/
def rescale (D : CapState α) (t : ℝ) : CapState α where
  isBulk := D.isBulk
  bulk_nonempty := D.bulk_nonempty
  weight := fun x => Real.exp (-t) * D.weight x
  weight_pos := fun x => mul_pos (Real.exp_pos (-t)) (D.weight_pos x)
  boost := D.boost
  boost_collar_zero := D.boost_collar_zero

theorem modularGenerator_rescale (D : CapState α) (t : ℝ) (x : α) :
    modularGenerator (rescale D t) x = modularGenerator D x + t := by
  show -Real.log (Real.exp (-t) * D.weight x) = -Real.log (D.weight x) + t
  rw [Real.log_mul (Real.exp_ne_zero (-t)) (ne_of_gt (D.weight_pos x)),
    Real.log_exp]
  ring

theorem bulkMass_rescale (D : CapState α) (t : ℝ) :
    bulkMass (rescale D t) = Real.exp (-t) * bulkMass D := by
  calc bulkMass (rescale D t)
      = ∑ x ∈ bulkSet D, Real.exp (-t) * D.weight x := rfl
    _ = Real.exp (-t) * ∑ x ∈ bulkSet D, D.weight x := by
        rw [Finset.mul_sum]
    _ = Real.exp (-t) * bulkMass D := rfl

theorem bulkBoostGap_rescale (D : CapState α) (t : ℝ) :
    bulkBoostGap (rescale D t)
      = Real.exp (-t) * (bulkBoostGap D + t * bulkMass D) := by
  have h : ∀ x ∈ bulkSet D,
      Real.exp (-t) * D.weight x *
          (modularGenerator (rescale D t) x - 2 * Real.pi * D.boost x)
        = Real.exp (-t) *
            (D.weight x * (modularGenerator D x - 2 * Real.pi * D.boost x)
              + t * D.weight x) := by
    intro x _
    rw [modularGenerator_rescale]
    ring
  calc bulkBoostGap (rescale D t)
      = ∑ x ∈ bulkSet D, Real.exp (-t) * D.weight x *
          (modularGenerator (rescale D t) x - 2 * Real.pi * D.boost x) := rfl
    _ = ∑ x ∈ bulkSet D, Real.exp (-t) *
          (D.weight x * (modularGenerator D x - 2 * Real.pi * D.boost x)
            + t * D.weight x) := Finset.sum_congr rfl h
    _ = Real.exp (-t) * ∑ x ∈ bulkSet D,
          (D.weight x * (modularGenerator D x - 2 * Real.pi * D.boost x)
            + t * D.weight x) := by rw [Finset.mul_sum]
    _ = Real.exp (-t) * (bulkBoostGap D + t * bulkMass D) := by
        rw [Finset.sum_add_distrib, ← Finset.mul_sum]
        rfl

/-- Under the rescaling the normalizing constant absorbs the shift. -/
theorem normalizingConstant_rescale (D : CapState α) (t : ℝ) :
    normalizingConstant (rescale D t) = normalizingConstant D + t := by
  have hb : bulkMass D ≠ 0 := ne_of_gt (bulkMass_pos D)
  have he : Real.exp (-t) ≠ 0 := Real.exp_ne_zero (-t)
  unfold normalizingConstant
  rw [bulkBoostGap_rescale, bulkMass_rescale]
  field_simp

/-- The anomalous part is quotient-invariant: the rescaling changes the
constant and nothing else. -/
theorem anomalousPart_rescale (D : CapState α) (t : ℝ) :
    anomalousPart (rescale D t) = anomalousPart D := by
  funext x
  show modularGenerator (rescale D t) x - 2 * Real.pi * D.boost x
      - normalizingConstant (rescale D t)
    = anomalousPart D x
  rw [modularGenerator_rescale, normalizingConstant_rescale]
  unfold anomalousPart
  ring

/-! ## Collar expectation and localization -/

/-- Normalized state expectation of the anomalous part over the collar
sector.  This is the finite stand-in for the anomalous cap energy. -/
def collarAnomalousEnergy (D : CapState α) : ℝ :=
  (∑ x ∈ collarSet D, D.weight x * anomalousPart D x) / totalMass D

/-- The collar expectation of the anomalous part is quotient-invariant. -/
theorem collarAnomalousEnergy_rescale (D : CapState α) (t : ℝ) :
    collarAnomalousEnergy (rescale D t) = collarAnomalousEnergy D := by
  have he : Real.exp (-t) ≠ 0 := Real.exp_ne_zero (-t)
  have hnum : ∑ x ∈ collarSet (rescale D t),
      (rescale D t).weight x * anomalousPart (rescale D t) x
      = Real.exp (-t) * ∑ x ∈ collarSet D, D.weight x * anomalousPart D x := by
    rw [anomalousPart_rescale]
    calc ∑ x ∈ collarSet D,
        Real.exp (-t) * D.weight x * anomalousPart D x
        = ∑ x ∈ collarSet D,
            Real.exp (-t) * (D.weight x * anomalousPart D x) := by
          apply Finset.sum_congr rfl
          intro x _
          ring
      _ = Real.exp (-t) * ∑ x ∈ collarSet D,
            D.weight x * anomalousPart D x := by rw [Finset.mul_sum]
  have htot : totalMass (rescale D t) = Real.exp (-t) * totalMass D := by
    calc totalMass (rescale D t)
        = ∑ x, Real.exp (-t) * D.weight x := rfl
      _ = Real.exp (-t) * ∑ x, D.weight x := by rw [Finset.mul_sum]
      _ = Real.exp (-t) * totalMass D := rfl
  unfold collarAnomalousEnergy
  rw [hnum, htot, mul_div_mul_left _ _ he]

/-- Gibbs-on-bulk condition: the bulk weights are exactly the Gibbs weights
of `2 pi B` at partition constant `logZ`. -/
def GibbsOnBulk (D : CapState α) (logZ : ℝ) : Prop :=
  ∀ x, D.isBulk x = true →
    D.weight x = Real.exp (-(2 * Real.pi * D.boost x) - logZ)

theorem boostGap_eq_logZ_on_bulk (D : CapState α) (logZ : ℝ)
    (h : GibbsOnBulk D logZ) (x : α) (hx : D.isBulk x = true) :
    modularGenerator D x - 2 * Real.pi * D.boost x = logZ := by
  unfold modularGenerator
  rw [h x hx, Real.log_exp]
  ring

/-- Under the Gibbs-on-bulk condition the normalizing constant is the
partition constant. -/
theorem normalizingConstant_eq_of_gibbs (D : CapState α) (logZ : ℝ)
    (h : GibbsOnBulk D logZ) : normalizingConstant D = logZ := by
  have hb : bulkMass D ≠ 0 := ne_of_gt (bulkMass_pos D)
  have hnum : bulkBoostGap D = bulkMass D * logZ := by
    unfold bulkBoostGap bulkMass
    rw [Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro x hx
    have hxB : D.isBulk x = true := by
      simpa [bulkSet] using hx
    rw [boostGap_eq_logZ_on_bulk D logZ h x hxB]
  unfold normalizingConstant
  rw [hnum, mul_comm, mul_div_assoc, div_self hb, mul_one]

/-- Localization on the bulk: Gibbs bulk weights make the anomalous part
vanish there. -/
theorem anomalousPart_eq_zero_on_bulk (D : CapState α) (logZ : ℝ)
    (h : GibbsOnBulk D logZ) (x : α) (hx : D.isBulk x = true) :
    anomalousPart D x = 0 := by
  unfold anomalousPart
  rw [normalizingConstant_eq_of_gibbs D logZ h]
  have hg := boostGap_eq_logZ_on_bulk D logZ h x hx
  linarith

/-- Localization on the collar: under the Gibbs-on-bulk condition every
microstate carrying nonzero anomalous part is a collar state. -/
theorem anomalousPart_supported_on_collar (D : CapState α) (logZ : ℝ)
    (h : GibbsOnBulk D logZ) (x : α) (hx : anomalousPart D x ≠ 0) :
    D.isBulk x = false := by
  cases hB : D.isBulk x with
  | false => rfl
  | true => exact absurd (anomalousPart_eq_zero_on_bulk D logZ h x hB) hx

/-- Under the Gibbs-on-bulk condition the full state expectation of the
anomalous part is its collar expectation: the anomalous energy is carried
entirely on the collar. -/
theorem anomalousExpectation_eq_collar (D : CapState α) (logZ : ℝ)
    (h : GibbsOnBulk D logZ) :
    (∑ x, D.weight x * anomalousPart D x) / totalMass D
      = collarAnomalousEnergy D := by
  unfold collarAnomalousEnergy
  congr 1
  rw [← Finset.sum_filter_add_sum_filter_not Finset.univ
    (fun x => D.isBulk x = true) (fun x => D.weight x * anomalousPart D x)]
  have hzero : ∑ x ∈ Finset.univ.filter (fun x => D.isBulk x = true),
      D.weight x * anomalousPart D x = 0 := by
    apply Finset.sum_eq_zero
    intro x hx
    have hxB : D.isBulk x = true := (Finset.mem_filter.mp hx).2
    rw [anomalousPart_eq_zero_on_bulk D logZ h x hxB, mul_zero]
  rw [hzero, zero_add]
  apply Finset.sum_congr _ (fun _ _ => rfl)
  ext x
  simp [collarSet, Bool.not_eq_true]

/-! ## The dark-sector interface -/

/-- The finite localization receipt: a defect hypothesis on the collar
expectation of the anomalous part builds exactly the `AnomalyRemainder`
consumed by the dark-sector interface. -/
def collarRemainder (D : CapState α) (ℓ C η : ℝ) (hℓ : 0 < ℓ)
    (hC : 0 ≤ C) (hη : 0 ≤ η)
    (hbound : |collarAnomalousEnergy D| ≤ C * η) : AnomalyRemainder where
  ell := ℓ
  ell_pos := hℓ
  C := C
  C_nonneg := hC
  eta := η
  eta_nonneg := hη
  anomalousEnergy := collarAnomalousEnergy D
  bound := hbound

/-- The remainder energy is the collar expectation of the split, by
definitional unfolding. -/
theorem collarRemainder_energy_eq (D : CapState α) (ℓ C η : ℝ)
    (hℓ : 0 < ℓ) (hC : 0 ≤ C) (hη : 0 ≤ η)
    (hbound : |collarAnomalousEnergy D| ≤ C * η) :
    (collarRemainder D ℓ C η hℓ hC hη hbound).anomalousEnergy
      = collarAnomalousEnergy D := rfl

/-- The dark-sector stress bound applies verbatim to the collar remainder of
a finite split. -/
theorem collarRemainder_stress_abs_le (D : CapState α) (ℓ C η : ℝ)
    (hℓ : 0 < ℓ) (hC : 0 ≤ C) (hη : 0 ≤ η)
    (hbound : |collarAnomalousEnergy D| ≤ C * η) :
    |anomalousStress (collarRemainder D ℓ C η hℓ hC hη hbound)|
      ≤ 15 * C * η / (8 * Real.pi ^ 2 * ℓ ^ 4) :=
  anomalousStress_abs_le _

/-! ## Explicit inhabitant -/

/-- Three-state inhabitant: two bulk states of weight `1`, one collar state
of weight `exp (-1)`, zero boost. -/
def exampleState : CapState (Fin 3) where
  isBulk := ![true, true, false]
  bulk_nonempty := ⟨0, rfl⟩
  weight := ![1, 1, Real.exp (-1)]
  weight_pos := by
    intro x
    fin_cases x <;> norm_num [Real.exp_pos]
  boost := fun _ => 0
  boost_collar_zero := fun _ _ => rfl

/-- The example satisfies the Gibbs-on-bulk condition at `logZ = 0`. -/
theorem exampleState_gibbs : GibbsOnBulk exampleState 0 := by
  intro x hx
  fin_cases x
  · show (1 : ℝ) = Real.exp (-(2 * Real.pi * 0) - 0)
    norm_num
  · show (1 : ℝ) = Real.exp (-(2 * Real.pi * 0) - 0)
    norm_num
  · simp [exampleState] at hx

/-- Exact normalizing constant of the example. -/
theorem exampleState_normalizingConstant :
    normalizingConstant exampleState = 0 := by
  have hgap : bulkBoostGap exampleState = 0 := by
    unfold bulkBoostGap bulkSet
    rw [Finset.sum_filter, Fin.sum_univ_three]
    norm_num [exampleState, modularGenerator, Matrix.cons_val_two,
      Matrix.tail_cons, Matrix.head_cons]
  unfold normalizingConstant
  rw [hgap, zero_div]

/-- Exact anomalous part of the example: zero on the bulk, `1` on the
collar. -/
theorem exampleState_anomalousPart :
    anomalousPart exampleState = ![0, 0, 1] := by
  funext x
  unfold anomalousPart
  rw [exampleState_normalizingConstant]
  fin_cases x <;> norm_num [modularGenerator, exampleState,
    Matrix.cons_val_two, Matrix.tail_cons, Matrix.head_cons]

/-- Exact collar expectation of the example. -/
theorem exampleState_collarEnergy :
    collarAnomalousEnergy exampleState
      = Real.exp (-1) / (2 + Real.exp (-1)) := by
  have hnum : ∑ x ∈ collarSet exampleState,
      exampleState.weight x * anomalousPart exampleState x
      = Real.exp (-1) := by
    unfold collarSet
    rw [Finset.sum_filter, Fin.sum_univ_three, exampleState_anomalousPart]
    norm_num [exampleState, Matrix.cons_val_two, Matrix.tail_cons,
      Matrix.head_cons]
  have htot : totalMass exampleState = 2 + Real.exp (-1) := by
    unfold totalMass
    rw [Fin.sum_univ_three]
    norm_num [exampleState, Matrix.cons_val_two, Matrix.tail_cons,
      Matrix.head_cons]
  unfold collarAnomalousEnergy
  rw [hnum, htot]

/-- The example's collar expectation is bounded by `(1/2) * 1`. -/
theorem exampleState_energy_le :
    |collarAnomalousEnergy exampleState| ≤ (1 / 2 : ℝ) * 1 := by
  rw [exampleState_collarEnergy]
  have h1 : (0 : ℝ) < Real.exp (-1) := Real.exp_pos (-1)
  have h2 : Real.exp (-1) < 1 := by
    rw [← Real.exp_zero]
    exact Real.exp_lt_exp.mpr (by norm_num)
  have hpos : (0 : ℝ) < 2 + Real.exp (-1) := by linarith
  rw [abs_of_pos (div_pos h1 hpos), div_le_iff₀ hpos]
  nlinarith

/-- Explicit remainder instance built from the example split, with
`ell = 1`, `C = 1/2`, `eta = 1`. -/
def exampleRemainder : AnomalyRemainder :=
  collarRemainder exampleState 1 (1 / 2) 1 one_pos (by norm_num)
    zero_le_one exampleState_energy_le

/-- The example remainder carries exactly the collar expectation of the
example split. -/
theorem exampleRemainder_energy :
    exampleRemainder.anomalousEnergy
      = Real.exp (-1) / (2 + Real.exp (-1)) :=
  exampleState_collarEnergy

/-- On collar states the boost datum vanishes, so the anomalous part reads
the modular generator directly: `K_anom = K - c` on the collar.  This is
the collar-support clause of the boost datum made explicit. -/
theorem anomalousPart_collar_boost_free (D : CapState α) (x : α)
    (hx : D.isBulk x = false) :
    anomalousPart D x = modularGenerator D x - normalizingConstant D := by
  unfold anomalousPart
  rw [D.boost_collar_zero x hx]
  ring

/-- Three-state inhabitant with a nonzero boost datum: bulk weights Gibbs
for the boost `(1, 0)` at `logZ = 0`, collar weight `exp (-1)`.  The
premise bundle is jointly satisfiable with a nontrivial boost. -/
def boostedState : CapState (Fin 3) where
  isBulk := ![true, true, false]
  bulk_nonempty := ⟨0, rfl⟩
  weight := ![Real.exp (-(2 * Real.pi)), 1, Real.exp (-1)]
  weight_pos := by
    intro x
    fin_cases x <;> norm_num [Real.exp_pos]
  boost := ![1, 0, 0]
  boost_collar_zero := by
    intro x hx
    fin_cases x
    · simp at hx
    · simp at hx
    · rfl

/-- The boosted inhabitant satisfies the Gibbs-on-bulk condition at
`logZ = 0` with its nonzero boost. -/
theorem boostedState_gibbs : GibbsOnBulk boostedState 0 := by
  intro x hx
  fin_cases x
  · show Real.exp (-(2 * Real.pi)) = Real.exp (-(2 * Real.pi * 1) - 0)
    norm_num
  · show (1 : ℝ) = Real.exp (-(2 * Real.pi * 0) - 0)
    norm_num
  · simp [boostedState] at hx

/-- The boosted inhabitant carries the same exact anomalous part
`(0, 0, 1)`: the split localizes on the collar with a nontrivial boost in
play. -/
theorem boostedState_anomalousPart :
    anomalousPart boostedState = ![0, 0, 1] := by
  funext x
  unfold anomalousPart
  rw [normalizingConstant_eq_of_gibbs boostedState 0 boostedState_gibbs]
  fin_cases x <;>
    norm_num [modularGenerator, boostedState, Matrix.cons_val_two,
      Matrix.tail_cons, Matrix.head_cons]

/-! ## Per-theorem axiom audit -/

#print axioms bulkMass_pos
#print axioms totalMass_pos
#print axioms modularGenerator_split
#print axioms anomalousPart_bulk_normalized
#print axioms split_unique
#print axioms split_eq_canonical
#print axioms split_nonunique_without_normalization
#print axioms modularGenerator_rescale
#print axioms bulkMass_rescale
#print axioms bulkBoostGap_rescale
#print axioms normalizingConstant_rescale
#print axioms anomalousPart_rescale
#print axioms collarAnomalousEnergy_rescale
#print axioms boostGap_eq_logZ_on_bulk
#print axioms normalizingConstant_eq_of_gibbs
#print axioms anomalousPart_eq_zero_on_bulk
#print axioms anomalousPart_supported_on_collar
#print axioms anomalousExpectation_eq_collar
#print axioms collarRemainder_energy_eq
#print axioms collarRemainder_stress_abs_le
#print axioms exampleState_gibbs
#print axioms exampleState_normalizingConstant
#print axioms exampleState_anomalousPart
#print axioms exampleState_collarEnergy
#print axioms exampleState_energy_le
#print axioms exampleRemainder_energy
#print axioms anomalousPart_collar_boost_free
#print axioms boostedState_gibbs
#print axioms boostedState_anomalousPart

end

end OPH.EinsteinBranch
