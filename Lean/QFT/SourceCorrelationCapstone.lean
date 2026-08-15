import QFT.TwoSlotCPNetWitness

/-!
# The source correlation capstone on a fully disjoint pair

The independent closure audit of issue #692 records two facts this
module answers.  The designated pair (86, 88) has overlapping full
supports, so its slot split does not reflect source disjointness; and
the two regional conditional expectations are not jointly injective on
the committed pair's carrier (`slotExpectations_not_jointly_injective`),
so this pair of marginal maps does not reconstruct an arbitrary joint state
and the CP net needs a marginal-coherent correlation datum; the same erasure is proved on
this module's carrier below.

The committed source-operator payload (schema v2, sha256
`765979faf7166126e8eb0c76078986fd5d2d35f2e40ded8bd3fce4b6e07d0a10`)
exports the verbatim 96-node full supports and the step-aligned joint
label paths.  Observers 86 and 247 have **fully disjoint** supports;
this module transcribes both support lists and kernel-checks their
disjointness, so the regional split of this pair is grounded in a
source datum rather than by declaration alone.

On the joint carrier of that pair, the source-counted joint empirical
state is the correlation receipt:

* the 32 aligned joint labels give a diagonal state with exact rational
  masses, certified positive semidefinite and trace one;
* its two partial traces equal the marginal empirical states exactly,
  through kernel-checked count marginalisation;
* it differs from the product of its marginals at an explicit entry, so
  the counted joint empirical state carries genuine correlation content;
* the two regional conditional expectations applied to it reproduce
  exactly the normalized marginal states.

Together: the expectations recover the marginals of the source state
and, by the on-carrier erasure theorem proved below, map the counted
state and its marginal product to the same expectation pair while the
two states differ; the joint state itself is therefore the
marginal-coherent correlation datum of the committed diamond on this
carrier, and it is source-counted, not declared.  The diamond itself is
instantiated on this carrier below, with the factor identifications of
both observers' source-generated algebras.  The separate constant-tower
transport is proved for the different 86/88 carrier in
`QFT.TowerAnchoredDiamond`; `QFT.TowerAnchoredCorrelation` anchors this
86/247 diamond to its own tower stage by the same committed pattern and
relates the two carriers through the shared observer 86, without
identifying the carriers.

**Boundary.**  The label conventions and the reading of two observers
as the two factors stay declared postprocessors; the support
disjointness, joint labels, counts, and marginals are committed source
data.  The extraction is post-hoc and ineligible as validation.  No
region lattice beyond the committed diamond, tower anchoring,
nonconstant tower,
channel semantics, instrument, clock, or physical claim is attached.
Of the reopened issue's three residual items, the justified regional
construction and the marginal-coherent correlation receipt are
addressed here at the stated strength; issue #728 owns the
nonconstant source realization or its scoped no-go.
-/

set_option maxRecDepth 65536

namespace OPH.QFT

open Matrix
open Kronecker
open scoped ComplexOrder

noncomputable section

/-! ## The kernel-checked source disjointness of the pair -/

/-- The verbatim 96-node full support of observer 86. -/
def support86Full : Fin 96 → ℕ :=
  ![86,44,52,65,73,78,94,99,107,120,128,141,154,10,18,23,31,36,57,26,39,
    60,81,115,70,91,112,133,146,102,136,149,167,162,175,188,209,170,183,
    196,217,230,222,243,0,2,5,7,13,15,47,28,49,41,62,21,34,68,89,123,157,
    204,83,104,125,180,201,159,214,235,144,178,191,225,238,256,251,264,
    277,298,259,272,285,306,319,290,311,332,1,3,4,6,8,11,16,12]

/-- The verbatim 96-node full support of observer 247. -/
def support247Full : Fin 96 → ℕ :=
  ![247,158,179,192,213,226,234,260,268,281,302,315,336,90,103,116,124,
    137,145,171,111,166,200,205,239,294,221,255,289,323,273,328,349,357,
    378,370,391,412,383,404,425,459,446,480,48,56,69,77,82,98,132,61,95,
    74,129,150,184,153,187,218,252,307,362,208,242,276,310,344,399,433,
    467,341,396,417,438,472,493,501,488,522,514,535,556,527,548,569,603,
    590,624,14,19,22,27,35,40,30]

/-- **Source disjointness.**  The two full supports share no node, so
the support disjointness backing the regional split of this pair is a
source datum. -/
theorem supports_disjoint :
    ∀ i j : Fin 96, support86Full i ≠ support247Full j := by decide

/-! ## The source-counted joint law -/

/-- The step-aligned joint label path of the pair `(86, 247)`,
transcribed from the committed payload. -/
def jointPath247 : Fin 32 → Fin 13 × Fin 13 :=
  ![(1,5),(2,4),(2,6),(7,6),(8,8),(10,10),(8,10),(9,9),(7,7),(2,7),
    (2,7),(8,3),(8,5),(6,1),(6,4),(7,4),(5,2),(8,3),(4,4),(4,4),(3,5),
    (0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(12,12),(12,12),(12,11),
    (11,11)]

/-- Transcription cross-check: the joint path's marginals are the
committed per-observer paths. -/
theorem jointPath247_marginals :
    ∀ t : Fin 32, (jointPath247 t).1 = obs86.path t ∧
      (jointPath247 t).2 = obs247.path t := by decide

/-- The joint occupation count of one joint label. -/
def jointCount (p : Fin 13 × Fin 13) : ℕ :=
  (Finset.univ.filter fun t => jointPath247 t = p).card

/-- The marginal occupation counts. -/
def count86 (i : Fin 13) : ℕ :=
  (Finset.univ.filter fun t => obs86.path t = i).card

def count247 (j : Fin 13) : ℕ :=
  (Finset.univ.filter fun t => obs247.path t = j).card

theorem jointCount_total : ∑ p : Fin 13 × Fin 13, jointCount p = 32 := by
  decide

theorem jointCount_marginal_left :
    ∀ i : Fin 13, ∑ j : Fin 13, jointCount (i, j) = count86 i := by decide

theorem jointCount_marginal_right :
    ∀ j : Fin 13, ∑ i : Fin 13, jointCount (i, j) = count247 j := by decide

/-- The counted table is not symmetric under exchanging the two observer
labels: the ordered cell `(1,5)` occurs once. -/
theorem jointCount_one_five : jointCount (1, 5) = 1 := by decide

/-- The transposed ordered cell `(5,1)` does not occur. -/
theorem jointCount_five_one : jointCount (5, 1) = 0 := by decide

/-- The source-counted joint empirical state. -/
def correlationState : Matrix (Fin 13 × Fin 13) (Fin 13 × Fin 13) ℂ :=
  Matrix.diagonal fun p => (jointCount p : ℂ) / 32

/-- The marginal empirical states. -/
def marginal86State : Matrix (Fin 13) (Fin 13) ℂ :=
  Matrix.diagonal fun i => (count86 i : ℂ) / 32

def marginal247State : Matrix (Fin 13) (Fin 13) ℂ :=
  Matrix.diagonal fun j => (count247 j : ℂ) / 32

/-- The joint state is positive semidefinite. -/
theorem correlationState_posSemidef : correlationState.PosSemidef := by
  refine Matrix.posSemidef_diagonal_iff.mpr fun p => ?_
  have h : (0 : ℝ) ≤ (jointCount p : ℝ) / 32 := by positivity
  rw [show ((jointCount p : ℂ) / 32) =
      (((jointCount p : ℝ) / 32 : ℝ) : ℂ) by push_cast; ring]
  exact_mod_cast h

/-- The joint state has trace one. -/
theorem correlationState_trace : correlationState.trace = 1 := by
  unfold correlationState
  rw [Matrix.trace_diagonal, ← Finset.sum_div]
  rw [show (∑ p : Fin 13 × Fin 13, (jointCount p : ℂ)) =
      ((∑ p : Fin 13 × Fin 13, jointCount p : ℕ) : ℂ) by push_cast; rfl]
  rw [jointCount_total]
  norm_num

/-! ## Partial traces of diagonal joint states -/

variable {α β : Type*} [Fintype α] [Fintype β] [DecidableEq α] [DecidableEq β]

omit [Fintype α] in
/-- The second partial trace of a diagonal joint matrix is the diagonal
of the fibre sums. -/
theorem ptraceSnd_diagonal (f : α × β → ℂ) :
    ptraceSnd (Matrix.diagonal f) =
      Matrix.diagonal fun a => ∑ b : β, f (a, b) := by
  ext a a'
  unfold ptraceSnd
  by_cases h : a = a'
  · subst h
    simp [Matrix.diagonal]
  · have : ∀ b : β, Matrix.diagonal f (a, b) (a', b) = 0 := by
      intro b
      apply Matrix.diagonal_apply_ne
      intro hc
      exact h (congrArg Prod.fst hc)
    simp [Matrix.diagonal_apply_ne _ h, this]

omit [Fintype β] in
/-- The first partial trace of a diagonal joint matrix is the diagonal
of the fibre sums. -/
theorem ptraceFst_diagonal (f : α × β → ℂ) :
    OPH.Locality.ptraceFst (Matrix.diagonal f) =
      Matrix.diagonal fun b => ∑ a : α, f (a, b) := by
  ext b b'
  unfold OPH.Locality.ptraceFst
  by_cases h : b = b'
  · subst h
    simp [Matrix.diagonal]
  · have : ∀ a : α, Matrix.diagonal f (a, b) (a, b') = 0 := by
      intro a
      apply Matrix.diagonal_apply_ne
      intro hc
      exact h (congrArg Prod.snd hc)
    simp [Matrix.diagonal_apply_ne _ h, this]

/-! ## Marginal consistency and the correlation receipts -/

/-- The left marginal of the source state is the left empirical
state. -/
theorem correlationState_marginal_left :
    ptraceSnd correlationState = marginal86State := by
  unfold correlationState marginal86State
  rw [ptraceSnd_diagonal]
  congr 1
  funext i
  rw [← Finset.sum_div]
  congr 1
  rw [show (∑ b : Fin 13, (jointCount (i, b) : ℂ)) =
      ((∑ b : Fin 13, jointCount (i, b) : ℕ) : ℂ) by push_cast; rfl]
  rw [jointCount_marginal_left i]

/-- The right marginal of the source state is the right empirical
state. -/
theorem correlationState_marginal_right :
    OPH.Locality.ptraceFst correlationState = marginal247State := by
  unfold correlationState marginal247State
  rw [ptraceFst_diagonal]
  congr 1
  funext j
  rw [← Finset.sum_div]
  congr 1
  rw [show (∑ a : Fin 13, (jointCount (a, j) : ℂ)) =
      ((∑ a : Fin 13, jointCount (a, j) : ℕ) : ℂ) by push_cast; rfl]
  rw [jointCount_marginal_right j]

/-- The count-level correlation witness at the freezeout label. -/
theorem correlation_witness :
    jointCount (0, 0) * 32 ≠ count86 0 * count247 0 := by decide

/-- **The counted joint empirical state is not the product of its
marginals**: it differs from the marginal product at the freezeout
entry, so the counted state carries genuine correlation content. -/
theorem correlationState_ne_product :
    correlationState ≠ marginal86State ⊗ₖ marginal247State := by
  intro h
  have hentry := congrFun (congrFun h ((0, 0) : Fin 13 × Fin 13))
    ((0, 0) : Fin 13 × Fin 13)
  unfold correlationState marginal86State marginal247State at hentry
  rw [Matrix.diagonal_apply_eq] at hentry
  simp only [Matrix.kroneckerMap_apply, Matrix.diagonal_apply_eq] at hentry
  have hj : jointCount (0, 0) = 7 := by decide
  have h86 : count86 0 = 7 := by decide
  have h247 : count247 0 = 7 := by decide
  rw [hj, h86, h247] at hentry
  norm_num at hentry

/-- **CP consumption.**  The right regional conditional expectation of
the source state is the normalized-left tensor of the right marginal:
the expectation reproduces the source marginal exactly. -/
theorem rightSlotExpectation_correlationState :
    rightSlotExpectation correlationState =
      (13 : ℂ)⁻¹ •
        ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ marginal247State) := by
  have h : rightSlotExpectation (α := Fin 13) (β := Fin 13)
      correlationState =
      ((Fintype.card (Fin 13) : ℂ))⁻¹ •
        ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
          OPH.Locality.ptraceFst correlationState) := rfl
  rw [h, correlationState_marginal_right]
  norm_num

/-- The left regional conditional expectation reproduces the left
marginal. -/
theorem leftSlotExpectation_correlationState :
    leftSlotExpectation correlationState =
      (13 : ℂ)⁻¹ •
        (marginal86State ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ)) := by
  have h : leftSlotExpectation (α := Fin 13) (β := Fin 13)
      correlationState =
      ((Fintype.card (Fin 13) : ℂ))⁻¹ •
        (ptraceSnd correlationState ⊗ₖ
          (1 : Matrix (Fin 13) (Fin 13) ℂ)) := rfl
  rw [h, correlationState_marginal_left]
  norm_num

/-! ## The committed diamond on the justified pair -/

/-- The (86, 247) carrier. -/
abbrev PairIndex247 := Fin 13 × Fin 13

/-- Observer 86's generated algebra, lifted into the left slot of the
justified carrier, is exactly the left factor. -/
theorem obs86_lifted_eq_leftSlot247 :
    StarSubalgebra.map (slotLeft (Fin 13)) obs86.sourceAlgebra =
      (slotLeft (Fin 13) (α := Fin 13)).range := by
  rw [obs86_sourceAlgebra_eq_top]
  exact (StarAlgHom.range_eq_map_top _).symm

/-- Observer 247's generated algebra, lifted into the right slot, is
exactly the right factor. -/
theorem obs247_lifted_eq_rightSlot247 :
    StarSubalgebra.map (slotRight (Fin 13)) obs247.sourceAlgebra =
      (slotRight (Fin 13) (β := Fin 13)).range := by
  rw [obs247_sourceAlgebra_eq_top]
  exact (StarAlgHom.range_eq_map_top _).symm

/-- The regional algebra assignment of the justified diamond. -/
def twoSlotAlgebra247 :
    TwoSlotRegion → StarSubalgebra ℂ (Matrix PairIndex247 PairIndex247 ℂ)
  | TwoSlotRegion.bot => ⊥
  | TwoSlotRegion.left => (slotLeft (Fin 13) (α := Fin 13)).range
  | TwoSlotRegion.right => (slotRight (Fin 13) (β := Fin 13)).range
  | TwoSlotRegion.top => ⊤

/-- The conditional-expectation assignment of the justified diamond. -/
def twoSlotExpect247 : TwoSlotRegion →
    Matrix PairIndex247 PairIndex247 ℂ →ₗ[ℂ] Matrix PairIndex247 PairIndex247 ℂ
  | TwoSlotRegion.bot => scalarExpectation PairIndex247
  | TwoSlotRegion.left => leftSlotExpectation
  | TwoSlotRegion.right => rightSlotExpectation
  | TwoSlotRegion.top => LinearMap.id

theorem leftSlotExpectation_mem247 (M : Matrix PairIndex247 PairIndex247 ℂ) :
    leftSlotExpectation M ∈ (slotLeft (Fin 13) (α := Fin 13)).range := by
  have h : leftSlotExpectation (α := Fin 13) (β := Fin 13) M =
      (Fintype.card (Fin 13) : ℂ)⁻¹ •
        (ptraceSnd M ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ)) := rfl
  rw [h]
  exact SMulMemClass.smul_mem _ ⟨ptraceSnd M, rfl⟩

theorem rightSlotExpectation_mem247 (M : Matrix PairIndex247 PairIndex247 ℂ) :
    rightSlotExpectation M ∈ (slotRight (Fin 13) (β := Fin 13)).range := by
  have h : rightSlotExpectation (α := Fin 13) (β := Fin 13) M =
      (Fintype.card (Fin 13) : ℂ)⁻¹ •
        ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ OPH.Locality.ptraceFst M) := rfl
  rw [h]
  exact SMulMemClass.smul_mem _ ⟨OPH.Locality.ptraceFst M, rfl⟩

/-- **The justified-pair diamond.**  The (86, 247) carrier inhabits the
conditional-expectation net interface with the identified
source-generated slot algebras as regions. -/
def twoSlotNet247 : CPRegionalNet PairIndex247 where
  Region := TwoSlotRegion
  regionFintype := inferInstance
  regionNonempty := ⟨TwoSlotRegion.bot⟩
  regionLE U V := TwoSlotRegion.le U V = true
  regionLE_refl := TwoSlotRegion.le_refl'
  regionLE_trans {U V W} := TwoSlotRegion.le_trans' U V W
  regionLE_antisymm {U V} := TwoSlotRegion.le_antisymm' U V
  overlap := TwoSlotRegion.meet
  overlap_le_left := TwoSlotRegion.meet_le_left'
  overlap_le_right := TwoSlotRegion.meet_le_right'
  le_overlap {W U V} := TwoSlotRegion.le_meet' W U V
  disjoint U V := TwoSlotRegion.disj U V = true
  disjoint_symm {U V} := TwoSlotRegion.disj_symm' U V
  disjoint_irrefl := TwoSlotRegion.disj_irrefl'
  localAlgebra := twoSlotAlgebra247
  isotony {U V} hUV := by
    cases U <;> cases V <;> (try (simp [TwoSlotRegion.le] at hUV))
    · exact le_refl _
    · exact bot_le
    · exact bot_le
    · exact bot_le
    · exact le_refl _
    · exact le_top
    · exact le_refl _
    · exact le_top
    · exact le_refl _
  locality {U V} hUV := by
    cases U <;> cases V <;> (try (simp [TwoSlotRegion.disj] at hUV))
    · exact fun X hX Y hY => slot_commute hX hY
    · exact fun X hX Y hY => (slot_commute hY hX).symm
  expect := twoSlotExpect247
  expect_mem U X := by
    cases U
    · exact scalarExpectation_mem_bot X
    · exact leftSlotExpectation_mem247 X
    · exact rightSlotExpectation_mem247 X
    · trivial
  expect_fixes U X hX := by
    cases U
    · exact scalarExpectation_fixes X hX
    · obtain ⟨A, rfl⟩ := hX
      exact leftSlotExpectation_fixes_left A
    · obtain ⟨B, rfl⟩ := hX
      exact rightSlotExpectation_fixes_right B
    · rfl
  expect_posSemidef U {X} hX := by
    cases U
    · exact scalarExpectation_posSemidef hX
    · exact leftSlotExpectation_posSemidef hX
    · exact rightSlotExpectation_posSemidef hX
    · exact hX
  expect_trace U X := by
    cases U
    · exact scalarExpectation_trace X
    · exact leftSlotExpectation_trace X
    · exact rightSlotExpectation_trace X
    · rfl
  expect_tower {U V} hUV X := by
    cases U <;> cases V <;> (try (simp [TwoSlotRegion.le] at hUV))
    · exact scalarExpectation_of_trace_eq (scalarExpectation_trace X)
    · exact scalarExpectation_of_trace_eq (leftSlotExpectation_trace X)
    · exact scalarExpectation_of_trace_eq (rightSlotExpectation_trace X)
    · rfl
    · exact leftSlotExpectation_idem X
    · rfl
    · exact rightSlotExpectation_idem X
    · rfl
    · rfl
  generating := {TwoSlotRegion.left, TwoSlotRegion.right}
  coverage := by
    have hset : (⋃ U ∈ ({TwoSlotRegion.left, TwoSlotRegion.right} :
        Finset TwoSlotRegion),
          (twoSlotAlgebra247 U : Set (Matrix PairIndex247 PairIndex247 ℂ))) =
        (((slotLeft (Fin 13) (α := Fin 13)).range :
            Set (Matrix PairIndex247 PairIndex247 ℂ)) ∪
          ((slotRight (Fin 13) (β := Fin 13)).range :
            Set (Matrix PairIndex247 PairIndex247 ℂ))) := by
      ext x
      simp only [Finset.mem_insert, Finset.mem_singleton, Set.mem_iUnion,
        Set.mem_union, exists_prop]
      constructor
      · rintro ⟨U, (rfl | rfl), hx⟩
        · exact Or.inl hx
        · exact Or.inr hx
      · rintro (hx | hx)
        · exact ⟨TwoSlotRegion.left, Or.inl rfl, hx⟩
        · exact ⟨TwoSlotRegion.right, Or.inr rfl, hx⟩
    rw [hset]
    exact slot_ranges_generate_top

/-! ## The on-carrier erasure theorem -/

/-- Partial second trace of a Kronecker product. -/
theorem ptraceSnd_kronecker {α β : Type*} [Fintype α] [Fintype β]
    (A : Matrix α α ℂ) (B : Matrix β β ℂ) :
    ptraceSnd (A ⊗ₖ B) = B.trace • A := by
  ext a a'
  unfold ptraceSnd
  simp only [Matrix.kroneckerMap_apply, Matrix.smul_apply, smul_eq_mul]
  have htr : B.trace = ∑ b : β, B b b := rfl
  rw [htr, Finset.sum_mul]
  exact Finset.sum_congr rfl fun b _ => mul_comm _ _

/-- Partial first trace of a Kronecker product. -/
theorem ptraceFst_kronecker {α β : Type*} [Fintype α] [Fintype β]
    (A : Matrix α α ℂ) (B : Matrix β β ℂ) :
    OPH.Locality.ptraceFst (A ⊗ₖ B) = A.trace • B := by
  ext b b'
  unfold OPH.Locality.ptraceFst
  simp only [Matrix.kroneckerMap_apply, Matrix.smul_apply, smul_eq_mul]
  have htr : A.trace = ∑ a : α, A a a := rfl
  rw [htr, Finset.sum_mul]

/-- The marginal empirical states have trace one. -/
theorem marginal86State_trace : marginal86State.trace = 1 := by
  unfold marginal86State
  rw [Matrix.trace_diagonal, ← Finset.sum_div]
  rw [show (∑ i : Fin 13, (count86 i : ℂ)) =
      ((∑ i : Fin 13, count86 i : ℕ) : ℂ) by push_cast; rfl]
  have h : ∑ i : Fin 13, count86 i = 32 := by decide
  rw [h]
  norm_num

theorem marginal247State_trace : marginal247State.trace = 1 := by
  unfold marginal247State
  rw [Matrix.trace_diagonal, ← Finset.sum_div]
  rw [show (∑ j : Fin 13, (count247 j : ℂ)) =
      ((∑ j : Fin 13, count247 j : ℕ) : ℂ) by push_cast; rfl]
  have h : ∑ j : Fin 13, count247 j = 32 := by decide
  rw [h]
  norm_num

/-- **On-carrier erasure.**  The two regional expectations map the
counted joint state and the product of its marginals to the same
expectation pair, while the two states differ: the expectation pair is
not injective at exactly the source correlation, so the joint state is
the datum the net must carry. -/
theorem slotExpectations_erase_source_correlation :
    (leftSlotExpectation correlationState,
        rightSlotExpectation correlationState) =
      (leftSlotExpectation (marginal86State ⊗ₖ marginal247State),
        rightSlotExpectation (marginal86State ⊗ₖ marginal247State)) ∧
      correlationState ≠ marginal86State ⊗ₖ marginal247State := by
  refine ⟨?_, correlationState_ne_product⟩
  have hL : leftSlotExpectation (marginal86State ⊗ₖ marginal247State) =
      leftSlotExpectation correlationState := by
    have h1 : leftSlotExpectation (α := Fin 13) (β := Fin 13)
        (marginal86State ⊗ₖ marginal247State) =
        (Fintype.card (Fin 13) : ℂ)⁻¹ •
          (ptraceSnd (marginal86State ⊗ₖ marginal247State) ⊗ₖ
            (1 : Matrix (Fin 13) (Fin 13) ℂ)) := rfl
    have h2 : leftSlotExpectation (α := Fin 13) (β := Fin 13)
        correlationState =
        (Fintype.card (Fin 13) : ℂ)⁻¹ •
          (ptraceSnd correlationState ⊗ₖ
            (1 : Matrix (Fin 13) (Fin 13) ℂ)) := rfl
    rw [h1, h2, ptraceSnd_kronecker, marginal247State_trace,
      correlationState_marginal_left, one_smul]
  have hR : rightSlotExpectation (marginal86State ⊗ₖ marginal247State) =
      rightSlotExpectation correlationState := by
    have h1 : rightSlotExpectation (α := Fin 13) (β := Fin 13)
        (marginal86State ⊗ₖ marginal247State) =
        (Fintype.card (Fin 13) : ℂ)⁻¹ •
          ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
            OPH.Locality.ptraceFst (marginal86State ⊗ₖ marginal247State)) := rfl
    have h2 : rightSlotExpectation (α := Fin 13) (β := Fin 13)
        correlationState =
        (Fintype.card (Fin 13) : ℂ)⁻¹ •
          ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
            OPH.Locality.ptraceFst correlationState) := rfl
    rw [h1, h2, ptraceFst_kronecker, marginal86State_trace,
      correlationState_marginal_right, one_smul]
  rw [hL, hR]

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.supports_disjoint
#print axioms OPH.QFT.jointPath247_marginals
#print axioms OPH.QFT.jointCount_total
#print axioms OPH.QFT.jointCount_marginal_left
#print axioms OPH.QFT.jointCount_marginal_right
#print axioms OPH.QFT.ptraceFst_diagonal
#print axioms OPH.QFT.correlationState_posSemidef
#print axioms OPH.QFT.correlationState_trace
#print axioms OPH.QFT.ptraceSnd_diagonal
#print axioms OPH.QFT.correlationState_marginal_left
#print axioms OPH.QFT.correlationState_marginal_right
#print axioms OPH.QFT.correlation_witness
#print axioms OPH.QFT.correlationState_ne_product
#print axioms OPH.QFT.rightSlotExpectation_correlationState
#print axioms OPH.QFT.leftSlotExpectation_correlationState
#print axioms OPH.QFT.obs86_lifted_eq_leftSlot247
#print axioms OPH.QFT.obs247_lifted_eq_rightSlot247
#print axioms OPH.QFT.twoSlotNet247
#print axioms OPH.QFT.ptraceSnd_kronecker
#print axioms OPH.QFT.ptraceFst_kronecker
#print axioms OPH.QFT.slotExpectations_erase_source_correlation
