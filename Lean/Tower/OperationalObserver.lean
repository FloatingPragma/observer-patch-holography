import Tower.ConsensusTower
import Time.ObserverHistory

/-!
# Operational observerhood receipts over consensus towers

This module implements completion-plan issue `#711`.  An
`OperationalObserver` over a `ConsensusTower` packages the seven operational
observerhood clauses as typed fields: a bounded accessible interface, a
self-readback map, durable records, record-conditioned control, predictive
boundary use, checkpoint and refinement continuation, and overlap-readable
evidence.  Every clause is a finite statement about declared tower data;
nothing is inferred from a limit, a dynamics, or a physical interpretation.

`witnessObserver` is an explicit inhabitant over the explicit one-regulator
tower `witnessTower`: the private algebra is `M₂(ℂ)`, the public algebra is
the span of the two diagonal coordinate projectors, three records carry the
strict chain order, the readback is the partition average, and the
prediction map advances records.  The tower generator is the `i`-scaled
commutator with the diagonal Hamiltonian `diag(1, 2)`, which annihilates
every diagonal record element and moves the raising matrix unit, so record
durability is a computed fixed-point fact rather than a vanishing-generator
triviality.  Every record-conditioned intervention is conjugation by the
diagonal involution `diag(1, -1)`, which fixes every diagonal record
element and negates the raising matrix unit, so evidence invariance under
control is likewise computed against a map distinct from the identity.
`witnessObserver_nontrivial` collects the facts separating this witness
from the constant adaptor: two-dimensional private algebra, strictly
ordered records, a nonconstant non-identity prediction map, a record
element distinct from the zero and identity matrices, a nonzero generator,
and a control distinct from the identity map.  Three negative controls show
that an erasing readback, an identity generator, and a resetting prediction
each violate exactly one clause on the same concrete data.

## Claim boundary

The receipt is operational and finite.  It selects no unique observer over
a given tower, attaches no physical instrument or laboratory system, makes
no claim about consciousness, and supplies no source realization of the
tower data.  The E2 source-realized observer family quantifies over this
receipt, and the E6 chart-owner typing may consume it; neither consumer is
constructed here.
-/

namespace OPH.Tower

open EventAlgebra
open Matrix
open OPH.D1Time

universe u

variable {ι : Type u} [Preorder ι]

/-- Operational observerhood receipt over a consensus tower.

An inhabitant certifies, for one refinement-compatible observer label
family, the seven operational clauses of completion-plan issue `#711`.  The
receipt consumes only declared tower data.  It asserts no uniqueness, no
physical attachment, and no source realization. -/
structure OperationalObserver (T : ConsensusTower ι) where
  /-- Clause 6 (checkpoint and refinement continuation), identity part: the
  observer is a regulator-indexed family of tower observer labels. -/
  label : ∀ r : ι, T.Observer r
  /-- Clause 6 (checkpoint and refinement continuation): the label family
  commutes with the declared observer refinement maps, so the observer
  persists across tower refinement. -/
  label_refine : ∀ {r s : ι} (hrs : r ≤ s),
    T.observerRefine hrs (label r) = label s
  /-- Clause 1 (bounded accessible interface): a declared accessible star
  subalgebra of the finite private matrix algebra at each regulator. -/
  accessible : ∀ r : ι, StarSubalgebra ℂ (T.PrivateAlgebra r)
  /-- Clause 1 (bounded accessible interface): the observer's own public
  record algebra lies inside the declared interface, so own records are
  accessible to their owner. -/
  public_le_accessible : ∀ r : ι,
    T.publicAlgebra r (label r) ≤ accessible r
  /-- Clause 2 (self-readback): a declared readback map on the private
  algebra.  The map is total; its receipt content is the two laws below,
  which constrain it on the accessible interface and on public elements. -/
  readback : ∀ r : ι, T.PrivateAlgebra r → T.PrivateAlgebra r
  /-- Clause 2 (self-readback): readback of an accessible element lands in
  the observer's own public record algebra. -/
  readback_mem_public : ∀ (r : ι) (X : T.PrivateAlgebra r),
    X ∈ accessible r → readback r X ∈ T.publicAlgebra r (label r)
  /-- Clause 2 (self-readback): reading an own public record back returns
  it. -/
  readback_fixes_public : ∀ (r : ι) (X : T.PrivateAlgebra r),
    X ∈ T.publicAlgebra r (label r) → readback r X = X
  /-- Clause 3 (durable records): every record element is annihilated by
  the observer's declared generator.  Records are fixed points of the
  declared persistence law. -/
  record_durable : ∀ (r : ι) (x : T.Record r),
    T.generator r (label r) (T.recordElement r (label r) x).1 = 0
  /-- Clause 4 (record-conditioned control): a record-indexed family of
  declared linear interventions on the private algebra. -/
  control : ∀ r : ι, T.Record r →
    T.PrivateAlgebra r →ₗ[ℂ] T.PrivateAlgebra r
  /-- Clause 4 (record-conditioned control): every declared intervention
  preserves the accessible interface. -/
  control_preserves_accessible : ∀ (r : ι) (x : T.Record r)
      (X : T.PrivateAlgebra r),
    X ∈ accessible r → control r x X ∈ accessible r
  /-- Clause 5 (predictive boundary use): a declared record-to-record
  prediction map. -/
  predict : ∀ r : ι, T.Record r → T.Record r
  /-- Clause 5 (predictive boundary use): the prediction output never
  strictly precedes its input in the observer's own record order. -/
  predict_no_regress : ∀ (r : ι) (x : T.Record r),
    ¬ (T.recordOrder r (label r)).precedes (predict r x) x
  /-- Clause 6 (checkpoint and refinement continuation): prediction
  commutes with the declared record refinement maps. -/
  predict_refine : ∀ {r s : ι} (hrs : r ≤ s) (x : T.Record r),
    T.recordRefine hrs (predict r x) = predict s (T.recordRefine hrs x)
  /-- Clause 7 (overlap-readable evidence): a declared readable value for
  each record.  The record element itself is typed into the observer's own
  public algebra by the tower field `recordElement`. -/
  readout : ∀ r : ι, T.Record r → ℂ
  /-- Clause 7 (overlap-readable evidence): the declared readable value is
  the trace pairing of the record element against the observer's own
  state. -/
  readout_eq_pairing : ∀ (r : ι) (x : T.Record r),
    readout r x =
      (T.state r (label r) * (T.recordElement r (label r) x).1).trace
  /-- Clause 7 (overlap-readable evidence): every record element is fixed
  by every declared intervention, so committed evidence survives the
  observer's own record-conditioned control. -/
  record_fixed_by_control : ∀ (r : ι) (x y : T.Record r),
    control r y (T.recordElement r (label r) x).1 =
      (T.recordElement r (label r) x).1

namespace OperationalObserver

variable {T : ConsensusTower ι} (O : OperationalObserver T)

/-- Reading an own record element back returns it: the record-element
special case of the self-readback law. -/
theorem readback_recordElement (r : ι) (x : T.Record r) :
    O.readback r (T.recordElement r (O.label r) x).1 =
      (T.recordElement r (O.label r) x).1 :=
  O.readback_fixes_public r _ (T.recordElement r (O.label r) x).2

/-- The readback is idempotent on the accessible interface. -/
theorem readback_idem (r : ι) (X : T.PrivateAlgebra r)
    (hX : X ∈ O.accessible r) :
    O.readback r (O.readback r X) = O.readback r X :=
  O.readback_fixes_public r _ (O.readback_mem_public r X hX)

/-- The declared readable value of a record is stable under refinement:
the refined record read at the finer regulator carries the coarse value.
This combines the continuation clause with the evidence clause. -/
theorem readout_refine {r s : ι} (hrs : r ≤ s) (x : T.Record r) :
    O.readout s (T.recordRefine hrs x) = O.readout r x := by
  rw [O.readout_eq_pairing, O.readout_eq_pairing, ← O.label_refine hrs,
    ← T.record_element_natural hrs (O.label r) x,
    T.state_natural hrs (O.label r)]

end OperationalObserver

/-! ## An explicit nontrivial witness

The witness data are exact and small: private algebra `M₂(ℂ)`, the diagonal
two-projector partition as public algebra, three chain-ordered records, the
first coordinate projector as state, the scaled Hamiltonian commutator as
generator, the partition average as readback, conjugation by a diagonal
involution as record-conditioned control, and a record-advancing prediction
map. -/

/-- The diagonal two-projector partition of `M₂(ℂ)`: the coordinate
projectors `diag(1,0)` and `diag(0,1)`. -/
noncomputable def witnessPartition : ProjectivePartition 2 2 where
  proj i := Matrix.diagonal (Pi.single i 1)
  isEvent i := by
    refine ⟨?_, ?_⟩
    · show (Matrix.diagonal (Pi.single i (1 : ℂ)))ᴴ =
        Matrix.diagonal (Pi.single i 1)
      rw [Matrix.diagonal_conjTranspose]
      congr 1
      funext j
      rcases eq_or_ne j i with h | h <;>
        simp [h, Pi.single_apply]
    · rw [Matrix.diagonal_mul_diagonal]
      congr 1
      funext j
      rcases eq_or_ne j i with h | h <;>
        simp [h, Pi.single_apply]
  orthogonal i j hij := by
    rw [Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_zero]
    congr 1
    funext a
    rcases eq_or_ne a i with h | h
    · subst h
      simp [hij]
    · simp [h]
  complete := by
    ext a b
    rw [Matrix.sum_apply]
    rcases eq_or_ne a b with h | h
    · subst h
      simp [Matrix.diagonal_apply_eq, Pi.single_apply, Matrix.one_apply_eq]
    · simp [h]

/-- The first witness projector has unit trace. -/
theorem witnessPartition_proj_zero_trace :
    (witnessPartition.proj 0).trace = 1 := by
  show (Matrix.diagonal (Pi.single (0 : Fin 2) (1 : ℂ))).trace = 1
  rw [Matrix.trace_diagonal]
  simp [Pi.single_apply]

/-- The first witness projector is a density matrix. -/
theorem witnessState_isState : IsState (witnessPartition.proj 0) :=
  ⟨(witnessPartition.isEvent 0).posSemidef, witnessPartition_proj_zero_trace⟩

/-- The first witness projector is nonzero. -/
theorem witnessPartition_proj_zero_ne_zero : witnessPartition.proj 0 ≠ 0 := by
  intro h
  have h00 : (1 : ℂ) = 0 := congrFun (congrFun h 0) 0
  exact one_ne_zero h00

/-- The first witness projector differs from the identity. -/
theorem witnessPartition_proj_zero_ne_one : witnessPartition.proj 0 ≠ 1 := by
  intro h
  have h11 : (0 : ℂ) = 1 := congrFun (congrFun h 1) 1
  exact zero_ne_one h11

/-- Record labels represented on the partition: the first record on the
first projector, the two later records on the second projector.
Injectivity is neither claimed nor needed. -/
def witnessRecordIndex : Fin 3 → Fin 2 := fun x => if x = 0 then 0 else 1

/-- The witness prediction map: advance each record by one step and hold
the final record fixed. -/
def witnessPredict : Fin 3 → Fin 3 := fun x => if x = 2 then 2 else x + 1

/-- The declared witness readouts: unit weight on the first record, zero
weight on the later records. -/
def witnessReadout : Fin 3 → ℂ := fun x => if x = 0 then 1 else 0

/-- The witness Hamiltonian `diag(1, 2)`.  Its distinct diagonal entries
make the commutator generator vanish exactly on the diagonal matrices. -/
noncomputable def witnessHamiltonian : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.diagonal ![1, 2]

/-- The witness generator: the `i`-scaled commutator with the witness
Hamiltonian.  It annihilates every diagonal matrix and moves the raising
matrix unit. -/
noncomputable def witnessGenerator :
    Matrix (Fin 2) (Fin 2) ℂ →ₗ[ℂ] Matrix (Fin 2) (Fin 2) ℂ where
  toFun X := Complex.I • (witnessHamiltonian * X - X * witnessHamiltonian)
  map_add' X Y := by
    rw [mul_add, add_mul, add_sub_add_comm, smul_add]
  map_smul' c X := by
    rw [RingHom.id_apply, mul_smul_comm, smul_mul_assoc, ← smul_sub,
      smul_comm]

@[simp] theorem witnessGenerator_apply (X : Matrix (Fin 2) (Fin 2) ℂ) :
    witnessGenerator X =
      Complex.I • (witnessHamiltonian * X - X * witnessHamiltonian) :=
  rfl

/-- The witness generator annihilates every diagonal matrix: diagonal
matrices commute with the diagonal witness Hamiltonian. -/
theorem witnessGenerator_diagonal (d : Fin 2 → ℂ) :
    witnessGenerator (Matrix.diagonal d) = 0 := by
  rw [witnessGenerator_apply, witnessHamiltonian,
    Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  have hcomm : (fun i => (![1, 2] : Fin 2 → ℂ) i * d i) =
      fun i => d i * (![1, 2] : Fin 2 → ℂ) i := by
    funext i
    exact mul_comm _ _
  rw [hcomm, sub_self, smul_zero]

/-- The witness intervention conjugator `diag(1, -1)`: a diagonal
involution. -/
noncomputable def witnessConjugator : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.diagonal ![1, -1]

/-- The witness conjugator squares to the identity, so conjugation by it is
conjugation by a self-inverse element. -/
theorem witnessConjugator_involution :
    witnessConjugator * witnessConjugator = 1 := by
  rw [witnessConjugator, Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]
  congr 1
  funext i
  fin_cases i <;> norm_num

/-- The witness record-conditioned control: conjugation by the diagonal
involution.  With the conjugator self-inverse, `X ↦ U * X * U` is exactly
`X ↦ U * X * U⁻¹`.  It fixes every diagonal matrix and negates the raising
matrix unit. -/
noncomputable def witnessControl :
    Matrix (Fin 2) (Fin 2) ℂ →ₗ[ℂ] Matrix (Fin 2) (Fin 2) ℂ where
  toFun X := witnessConjugator * X * witnessConjugator
  map_add' X Y := by
    rw [mul_add, add_mul]
  map_smul' c X := by
    rw [RingHom.id_apply, mul_smul_comm, smul_mul_assoc]

@[simp] theorem witnessControl_apply (X : Matrix (Fin 2) (Fin 2) ℂ) :
    witnessControl X = witnessConjugator * X * witnessConjugator :=
  rfl

/-- The witness control fixes every diagonal matrix: diagonal matrices
commute with the conjugator and the conjugator squares to the identity. -/
theorem witnessControl_diagonal (d : Fin 2 → ℂ) :
    witnessControl (Matrix.diagonal d) = Matrix.diagonal d := by
  rw [witnessControl_apply, witnessConjugator,
    Matrix.diagonal_mul_diagonal, Matrix.diagonal_mul_diagonal]
  congr 1
  funext i
  fin_cases i <;> simp

/-- The raising matrix unit `|0><1|`, restated over the witness algebra as
the off-diagonal probe for the generator and control facts. -/
def witnessRaise : Matrix (Fin 2) (Fin 2) ℂ :=
  fun i j => if i = 0 ∧ j = 1 then 1 else 0

/-- The witness generator negates and `i`-scales the raising matrix unit. -/
theorem witnessGenerator_raise :
    witnessGenerator witnessRaise = -(Complex.I • witnessRaise) := by
  rw [witnessGenerator_apply]
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [witnessHamiltonian, witnessRaise, Matrix.mul_apply,
      Matrix.diagonal_apply]

/-- The witness control negates the raising matrix unit. -/
theorem witnessControl_raise :
    witnessControl witnessRaise = -witnessRaise := by
  rw [witnessControl_apply]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [witnessConjugator, witnessRaise, Matrix.mul_apply,
      Fin.sum_univ_two, Matrix.diagonal_apply]

/-- An explicit one-regulator consensus tower over `M₂(ℂ)` with three
chain-ordered records represented on the diagonal partition, the first
coordinate projector as state, and the scaled Hamiltonian commutator as
generator. -/
noncomputable def witnessTower : ConsensusTower Unit where
  indexNonempty := ⟨()⟩
  commonRefinement := fun _ _ => ⟨(), by trivial, by trivial⟩
  dim := fun _ => 2
  Observer := fun _ => Unit
  observerFinite := fun _ => inferInstance
  Record := fun _ => Fin 3
  recordFinite := fun _ => inferInstance
  recordOrder := fun _ _ => finRecordOrder 3
  publicAlgebra := fun _ _ => witnessPartition.publicSubalgebra
  public_mul_comm := fun _ _ X Y =>
    witnessPartition.publicSubalgebra_mul_comm X Y
  recordElement := fun _ _ x =>
    ⟨witnessPartition.proj (witnessRecordIndex x),
      witnessPartition.proj_mem_span (witnessRecordIndex x)⟩
  state := fun _ _ => witnessPartition.proj 0
  state_isState := fun _ _ => witnessState_isState
  generator := fun _ _ => witnessGenerator
  observerRefine := fun _ o => o
  recordRefine := fun _ x => x
  algebraRefine := fun _ => StarAlgHom.id ℂ _
  observer_refine_refl := by intros; rfl
  observer_refine_trans := by intros; rfl
  record_refine_refl := by intros; rfl
  record_refine_trans := by intros; rfl
  algebra_refine_refl := by intros; rfl
  algebra_refine_trans := by intros; rfl
  record_order_natural := by
    intro r s hrs o x y h
    exact h
  public_natural := by
    intro r s hrs o X hX
    exact hX
  record_element_natural := by intros; rfl
  state_natural := by intros; rfl
  generator_natural := by intros; rfl

/-- The operational witness observer over `witnessTower`: full accessible
interface, partition-average readback, conjugation interventions by the
diagonal involution, advancing prediction, and explicit declared
readouts. -/
noncomputable def witnessObserver : OperationalObserver witnessTower where
  label := fun _ => ()
  label_refine := by intros; rfl
  accessible := fun _ => ⊤
  public_le_accessible := fun _ => le_top
  readback := fun _ X => partitionAverage witnessPartition X
  readback_mem_public := fun _ X _ =>
    partitionAverage_mem_span witnessPartition X
  readback_fixes_public := fun _ X hX =>
    partitionAverage_fixes witnessPartition hX
  record_durable := by
    intro r x
    exact witnessGenerator_diagonal (Pi.single (witnessRecordIndex x) 1)
  control := fun _ _ => witnessControl
  control_preserves_accessible := fun _ _ _ _ => StarSubalgebra.mem_top
  predict := fun _ => witnessPredict
  predict_no_regress := by
    intro r
    show ∀ x : Fin 3, ¬ witnessPredict x < x
    decide
  predict_refine := by intros; rfl
  readout := fun _ => witnessReadout
  readout_eq_pairing := by
    intro r
    show ∀ x : Fin 3,
      witnessReadout x =
        (witnessPartition.proj 0 *
          witnessPartition.proj (witnessRecordIndex x)).trace
    intro x
    fin_cases x
    · show (1 : ℂ) =
        (witnessPartition.proj 0 * witnessPartition.proj 0).trace
      rw [(witnessPartition.isEvent 0).2, witnessPartition_proj_zero_trace]
    · show (0 : ℂ) =
        (witnessPartition.proj 0 * witnessPartition.proj 1).trace
      rw [witnessPartition.orthogonal 0 1 (by decide), Matrix.trace_zero]
    · show (0 : ℂ) =
        (witnessPartition.proj 0 * witnessPartition.proj 1).trace
      rw [witnessPartition.orthogonal 0 1 (by decide), Matrix.trace_zero]
  record_fixed_by_control := by
    intro r x y
    exact witnessControl_diagonal (Pi.single (witnessRecordIndex x) 1)

/-- Two witness records are strictly ordered. -/
theorem witnessTower_strict_precedence :
    (witnessTower.recordOrder () ()).precedes (0 : Fin 3) (1 : Fin 3) := by
  show (0 : Fin 3) < 1
  decide

/-- The witness prediction advances the first record. -/
theorem witnessObserver_predict_zero :
    witnessObserver.predict () (0 : Fin 3) = (1 : Fin 3) := rfl

/-- The witness prediction map is nonconstant. -/
theorem witnessObserver_predict_nonconstant :
    witnessObserver.predict () (0 : Fin 3) ≠
      witnessObserver.predict () (1 : Fin 3) := by
  intro h
  have h' : (1 : Fin 3) = 2 := h
  exact absurd h' (by decide)

/-- The witness prediction map differs from the identity map. -/
theorem witnessObserver_predict_ne_id :
    witnessObserver.predict () ≠ (id : Fin 3 → Fin 3) := by
  intro h
  have h0 : (1 : Fin 3) = 0 := congrFun h (0 : Fin 3)
  exact absurd h0 (by decide)

/-- The witness generator is nonzero on the raising matrix unit, so record
durability is a fixed-point property of the diagonal records rather than a
property of a vanishing generator. -/
theorem witnessGenerator_raise_ne_zero :
    witnessGenerator witnessRaise ≠ 0 := by
  rw [witnessGenerator_raise]
  intro h
  have h01 := congrFun (congrFun h 0) 1
  simp [witnessRaise] at h01

/-- The witness generator is nonzero as a linear map. -/
theorem witnessGenerator_ne_zero :
    witnessGenerator ≠
      (0 : Matrix (Fin 2) (Fin 2) ℂ →ₗ[ℂ] Matrix (Fin 2) (Fin 2) ℂ) := by
  intro h
  exact witnessGenerator_raise_ne_zero
    ((congrArg (fun F => F witnessRaise) h).trans
      (LinearMap.zero_apply witnessRaise))

/-- The witness control differs from the identity map: it negates the
raising matrix unit while fixing every diagonal record element. -/
theorem witnessControl_ne_id :
    witnessControl ≠
      (LinearMap.id :
        Matrix (Fin 2) (Fin 2) ℂ →ₗ[ℂ] Matrix (Fin 2) (Fin 2) ℂ) := by
  intro h
  have hr : witnessControl witnessRaise = witnessRaise :=
    (LinearMap.congr_fun h witnessRaise).trans
      (LinearMap.id_apply witnessRaise)
  rw [witnessControl_raise] at hr
  have h01 := congrFun (congrFun hr 0) 1
  norm_num [witnessRaise] at h01

/-- Nontriviality receipt for the witness observer: two-dimensional
private algebra, strictly ordered records, a nonconstant non-identity
prediction map, a record element distinct from the zero and identity
matrices, a nonzero record-durability generator, and a record-conditioned
control distinct from the identity map. -/
theorem witnessObserver_nontrivial :
    witnessTower.dim () = 2 ∧
      (witnessTower.recordOrder () ()).precedes (0 : Fin 3) (1 : Fin 3) ∧
      witnessObserver.predict () (0 : Fin 3) ≠
        witnessObserver.predict () (1 : Fin 3) ∧
      witnessObserver.predict () ≠ (id : Fin 3 → Fin 3) ∧
      (witnessTower.recordElement () () (0 : Fin 3)).1 ≠ 0 ∧
      (witnessTower.recordElement () () (0 : Fin 3)).1 ≠ 1 ∧
      witnessTower.generator () () ≠ 0 ∧
      witnessObserver.control () (0 : Fin 3) ≠ LinearMap.id :=
  ⟨rfl, witnessTower_strict_precedence, witnessObserver_predict_nonconstant,
    witnessObserver_predict_ne_id, witnessPartition_proj_zero_ne_zero,
    witnessPartition_proj_zero_ne_one, witnessGenerator_ne_zero,
    witnessControl_ne_id⟩

/-! ## Negative controls

Each control presents explicit candidate data over the same witness tower
and derives a contradiction from the corresponding clause hypothesis, in
the exact shape of the structure field it violates. -/

/-- Negative-control candidate: a readback that erases every input. -/
def erasingReadback :
    witnessTower.PrivateAlgebra () → witnessTower.PrivateAlgebra () :=
  fun _ => 0

/-- Negative control (a): the erasing readback moves the public identity
element, so it cannot satisfy the self-readback law
`readback_fixes_public` over the witness tower. -/
theorem erasingReadback_not_self_readback :
    ¬ ∀ X : witnessTower.PrivateAlgebra (),
        X ∈ witnessTower.publicAlgebra () () → erasingReadback X = X := by
  intro h
  have h1 : erasingReadback 1 = 1 := h 1 (one_mem _)
  have h00 : (0 : ℂ) = 1 := congrFun (congrFun h1 (0 : Fin 2)) (0 : Fin 2)
  exact zero_ne_one h00

/-- Negative control (b): the identity generator fails to annihilate the
first record element, so it cannot satisfy the durability law
`record_durable` over the witness tower. -/
theorem idGenerator_not_durable :
    ¬ ∀ x : Fin 3,
        (LinearMap.id : witnessTower.PrivateAlgebra () →ₗ[ℂ]
            witnessTower.PrivateAlgebra ())
          (witnessTower.recordElement () () x).1 = 0 := by
  intro h
  exact witnessPartition_proj_zero_ne_zero (h 0)

/-- Negative-control candidate: a prediction that resets every record to
the first record. -/
def resettingPredict : Fin 3 → Fin 3 := fun _ => 0

/-- Negative control (c): the resetting prediction sends the second record
strictly backwards, so it cannot satisfy the no-regress law
`predict_no_regress` over the witness tower. -/
theorem resettingPredict_not_predictive :
    ¬ ∀ x : Fin 3,
        ¬ (witnessTower.recordOrder () ()).precedes (resettingPredict x) x := by
  intro h
  exact h 1 (show (0 : Fin 3) < 1 by decide)

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.Tower.OperationalObserver.readback_recordElement
#print axioms OPH.Tower.OperationalObserver.readback_idem
#print axioms OPH.Tower.OperationalObserver.readout_refine
#print axioms OPH.Tower.witnessGenerator_diagonal
#print axioms OPH.Tower.witnessConjugator_involution
#print axioms OPH.Tower.witnessControl_diagonal
#print axioms OPH.Tower.witnessGenerator_raise
#print axioms OPH.Tower.witnessControl_raise
#print axioms OPH.Tower.witnessGenerator_raise_ne_zero
#print axioms OPH.Tower.witnessGenerator_ne_zero
#print axioms OPH.Tower.witnessControl_ne_id
#print axioms OPH.Tower.witnessObserver_nontrivial
#print axioms OPH.Tower.erasingReadback_not_self_readback
#print axioms OPH.Tower.idGenerator_not_durable
#print axioms OPH.Tower.resettingPredict_not_predictive

end OPH.Tower
