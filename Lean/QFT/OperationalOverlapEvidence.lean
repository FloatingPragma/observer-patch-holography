import QFT.ObserverEventWorld
import Tower.OperationalObserver

/-!
# Cross-observer readable evidence at one finite access cut

This module supplies the cross-observer half of the operational-observerhood
receipt at one declared regulator.  Two `OperationalObserver`s over the same
consensus tower are compared through an E6 `ObserverAccessCut`.  Each tower
record is first regarded as a local section on its owner's region and then
restricted, through the net's typed restriction maps, to the meet of the two
owner regions.  `OperationalOverlapEvidence` requires distinct observer
labels, exact compatibility of both operational accessible algebras with the
cut, equality of every committed record after those two restrictions, equality
of the corresponding own readouts, and at least one nonzero common restricted
record.

The receipt is deliberately fixed-regulator.  It does not construct a family
of cuts, refinement-natural cross-observer evidence, triple-overlap coherence,
source realization, or a physical instrument.  Those remain downstream E2
obligations.

The synthetic witness below makes the restriction map load-bearing.  Its two
observer regions are distinct proper regions meeting in a scalar bottom
region.  Observer-specific diagonal record representatives differ before
restriction but agree after the two character restrictions, and one common
restriction is nonzero.  A corner-mutated tower keeps the same carriers and
all operational-observer fields but makes one overlap restriction and its own
readout disagree; it is an exact negative control.

Finally, `operationalObservers_share_visible_event_record` is an E6 consumer
independent of the synthetic record witness.  When one event is visible in
charts owned by two distinct operational observers, its one packet-derived
record restricts equally from both owner regions and is fixed by both
self-readback maps.
-/

namespace OPH.QFT

open Matrix
open EventAlgebra
open OPH.Tower
open OPH.D1Time

universe u

variable {ι : Type u} [Preorder ι]
variable {T : ConsensusTower ι} {N : FiniteCausalObserverNet T} {r : ι}

namespace ObserverAccessCut

variable (A : ObserverAccessCut T N r)

/-- A tower record, represented in one observer's public algebra, regarded as
a local section on that observer's declared access region.  Public containment
and localization of the access cut supply the membership proof. -/
def recordAtObserverRegion (o : T.Observer r) (x : T.Record r) :
    N.localAlgebra r (A.observerRegion o) :=
  ⟨(T.recordElement r o x).1,
    A.accessible_le_region o
      (A.public_le_accessible o (T.recordElement r o x).2)⟩

/-- Restrict the left observer's representative of a record to the declared
meet of the left and right observer regions. -/
noncomputable def recordRestrictionLeft (o p : T.Observer r) (x : T.Record r) :
    N.localAlgebra r
      (N.overlap r (A.observerRegion o) (A.observerRegion p)) :=
  N.restrict r
    (N.overlap_le_left r (A.observerRegion o) (A.observerRegion p))
    (A.recordAtObserverRegion o x)

/-- Restrict the right observer's representative of a record to the same
oriented meet used by `recordRestrictionLeft`. -/
noncomputable def recordRestrictionRight (o p : T.Observer r) (x : T.Record r) :
    N.localAlgebra r
      (N.overlap r (A.observerRegion o) (A.observerRegion p)) :=
  N.restrict r
    (N.overlap_le_right r (A.observerRegion o) (A.observerRegion p))
    (A.recordAtObserverRegion p x)

/-- The left restriction remains accessible to its source observer. -/
theorem recordRestrictionLeft_mem_accessible (o p : T.Observer r)
    (x : T.Record r) :
    (A.recordRestrictionLeft o p x : ConsensusTower.PrivateAlgebra T r) ∈
      A.accessibleAlgebra o :=
  A.restrict_preserves o
    (N.regionLE_refl r (A.observerRegion o))
    (N.overlap_le_left r (A.observerRegion o) (A.observerRegion p))
    (A.recordAtObserverRegion o x)
    (A.public_le_accessible o (T.recordElement r o x).2)

/-- The right restriction remains accessible to its source observer. -/
theorem recordRestrictionRight_mem_accessible (o p : T.Observer r)
    (x : T.Record r) :
    (A.recordRestrictionRight o p x : ConsensusTower.PrivateAlgebra T r) ∈
      A.accessibleAlgebra p :=
  A.restrict_preserves p
    (N.regionLE_refl r (A.observerRegion p))
    (N.overlap_le_right r (A.observerRegion o) (A.observerRegion p))
    (A.recordAtObserverRegion p x)
    (A.public_le_accessible p (T.recordElement r p x).2)

end ObserverAccessCut

/-- Cross-observer readable evidence for two operational observers at one E6
access cut.  Every record label is committed: the two observer-indexed public
representatives must agree after typed restriction to the owner-region meet,
and the two own readouts must agree.  The nonzero field prevents an empty or
all-zero discharge. -/
structure OperationalOverlapEvidence (A : ObserverAccessCut T N r)
    (O P : OperationalObserver T) : Prop where
  /-- The two operational observers have different tower labels at the cut. -/
  labels_ne : O.label r ≠ P.label r
  /-- The first operational interface is exactly the first cut interface. -/
  left_accessible_eq :
    O.accessible r = A.accessibleAlgebra (O.label r)
  /-- The second operational interface is exactly the second cut interface. -/
  right_accessible_eq :
    P.accessible r = A.accessibleAlgebra (P.label r)
  /-- All committed record representatives agree on the shared overlap. -/
  record_restriction_eq : ∀ x : T.Record r,
    A.recordRestrictionLeft (O.label r) (P.label r) x =
      A.recordRestrictionRight (O.label r) (P.label r) x
  /-- Both observers assign the same own readable value to every committed
  record label.  This is independent of equality of the restricted matrix
  representatives because the tower states are observer-indexed. -/
  readout_eq : ∀ x : T.Record r, O.readout r x = P.readout r x
  /-- At least one common restricted record is nonzero. -/
  exists_nonzero_restriction : ∃ x : T.Record r,
    A.recordRestrictionLeft (O.label r) (P.label r) x ≠ 0

namespace OperationalOverlapEvidence

variable {A : ObserverAccessCut T N r}
variable {O P : OperationalObserver T}

/-- The common overlap section is accessible to both operational observers.
The right-hand membership uses exact overlap equality, not an ambient cast. -/
theorem commonRecord_mem_both (E : OperationalOverlapEvidence A O P)
    (x : T.Record r) :
    (A.recordRestrictionLeft (O.label r) (P.label r) x :
        ConsensusTower.PrivateAlgebra T r) ∈ O.accessible r ∧
      (A.recordRestrictionLeft (O.label r) (P.label r) x :
        ConsensusTower.PrivateAlgebra T r) ∈ P.accessible r := by
  constructor
  · rw [E.left_accessible_eq]
    exact A.recordRestrictionLeft_mem_accessible _ _ x
  · rw [E.right_accessible_eq]
    have h := A.recordRestrictionRight_mem_accessible
      (O.label r) (P.label r) x
    rw [← E.record_restriction_eq x] at h
    exact h

end OperationalOverlapEvidence

/-! ## E6 event-record consumer -/

namespace FiniteEventGeometryWorld

variable {A : ObserverAccessCut T N r} {sig : EventWorldSignature}
variable {S : OPH.C2Soldering.EventFrameSoldering sig.Event sig.Chart}

/-- The record of a visible event, promoted from its packet support to the
seeing chart owner's full observer region. -/
def recordAtOwnerRegion (W : FiniteEventGeometryWorld A sig S)
    {i : sig.Chart} {e : sig.Event} (hi : W.activeChart i)
    (he : W.activeEvent e) (hv : S.visible i e) :
    N.localAlgebra r (A.observerRegion (W.chartOwner i)) :=
  ⟨W.recordRead (W.packetOf e),
    N.isotony r
      (N.regionLE_trans r (W.visible_support hi he hv)
        (W.chartRegion_le_owner i hi))
      (W.record_local (W.packetOf e))⟩

/-- E6 cross-observer consumer.  If one active event is visible in two active
charts owned by distinct operational observers, then its one packet-derived
record restricts to the same section from the two owner regions.  Moreover,
the record is public to both owners and hence fixed by both operational
self-readback maps.

The theorem consumes two `OperationalObserver` values but does not claim that
the world or either observer was source-realized. -/
theorem operationalObservers_share_visible_event_record
    (W : FiniteEventGeometryWorld A sig S) (O P : OperationalObserver T)
    {i j : sig.Chart} {e : sig.Event}
    (hi : W.activeChart i) (hj : W.activeChart j)
    (he : W.activeEvent e) (hvi : S.visible i e) (hvj : S.visible j e)
    (hoi : W.chartOwner i = O.label r)
    (hoj : W.chartOwner j = P.label r)
    (_hne : O.label r ≠ P.label r) :
    N.restrict r
        (N.overlap_le_left r
          (A.observerRegion (W.chartOwner i))
          (A.observerRegion (W.chartOwner j)))
        (W.recordAtOwnerRegion hi he hvi) =
      N.restrict r
        (N.overlap_le_right r
          (A.observerRegion (W.chartOwner i))
          (A.observerRegion (W.chartOwner j)))
        (W.recordAtOwnerRegion hj he hvj) ∧
      O.readback r (W.recordRead (W.packetOf e)) =
        W.recordRead (W.packetOf e) ∧
      P.readback r (W.recordRead (W.packetOf e)) =
        W.recordRead (W.packetOf e) := by
  have hsi : N.regionLE r (W.supportRead (W.packetOf e))
      (A.observerRegion (W.chartOwner i)) :=
    N.regionLE_trans r (W.visible_support hi he hvi)
      (W.chartRegion_le_owner i hi)
  have hsj : N.regionLE r (W.supportRead (W.packetOf e))
      (A.observerRegion (W.chartOwner j)) :=
    N.regionLE_trans r (W.visible_support hj he hvj)
      (W.chartRegion_le_owner j hj)
  have hso : N.regionLE r (W.supportRead (W.packetOf e))
      (N.overlap r (A.observerRegion (W.chartOwner i))
        (A.observerRegion (W.chartOwner j))) :=
    N.le_overlap r hsi hsj
  let Xo : N.localAlgebra r
      (N.overlap r (A.observerRegion (W.chartOwner i))
        (A.observerRegion (W.chartOwner j))) :=
    ⟨W.recordRead (W.packetOf e),
      N.isotony r hso (W.record_local (W.packetOf e))⟩
  have hleft : W.recordAtOwnerRegion hi he hvi =
      N.localInclusion
        (N.overlap_le_left r (A.observerRegion (W.chartOwner i))
          (A.observerRegion (W.chartOwner j))) Xo := by
    apply Subtype.ext
    rfl
  have hright : W.recordAtOwnerRegion hj he hvj =
      N.localInclusion
        (N.overlap_le_right r (A.observerRegion (W.chartOwner i))
          (A.observerRegion (W.chartOwner j))) Xo := by
    apply Subtype.ext
    rfl
  refine ⟨?_, ?_, ?_⟩
  · rw [hleft, hright]
    have hl := N.restrict_inclusion r
      (N.overlap_le_left r (A.observerRegion (W.chartOwner i))
        (A.observerRegion (W.chartOwner j))) Xo
    have hr := N.restrict_inclusion r
      (N.overlap_le_right r (A.observerRegion (W.chartOwner i))
        (A.observerRegion (W.chartOwner j))) Xo
    have hl' :
        N.restrict r
            (N.overlap_le_left r (A.observerRegion (W.chartOwner i))
              (A.observerRegion (W.chartOwner j)))
            (N.localInclusion
              (N.overlap_le_left r (A.observerRegion (W.chartOwner i))
                (A.observerRegion (W.chartOwner j))) Xo) = Xo := by
      simpa [FiniteCausalObserverNet.localInclusion] using hl
    have hr' :
        N.restrict r
            (N.overlap_le_right r (A.observerRegion (W.chartOwner i))
              (A.observerRegion (W.chartOwner j)))
            (N.localInclusion
              (N.overlap_le_right r (A.observerRegion (W.chartOwner i))
                (A.observerRegion (W.chartOwner j))) Xo) = Xo := by
      simpa [FiniteCausalObserverNet.localInclusion] using hr
    exact hl'.trans hr'.symm
  · have hpub := W.record_public hi he hvi
    rw [hoi] at hpub
    exact O.readback_fixes_public r _ hpub
  · have hpub := W.record_public hj he hvj
    rw [hoj] at hpub
    exact P.readback_fixes_public r _ hpub

end FiniteEventGeometryWorld

/-! ## A proper-meet two-observer witness and a corner-mutation control -/

section OverlapWitness

/-- Three region labels: two distinct observer regions and their proper
bottom meet. -/
inductive OperationalOverlapRegion : Type
  | bot
  | left
  | right
  deriving DecidableEq, Fintype

namespace OperationalOverlapRegion

/-- Region inclusion: reflexivity plus the bottom region below both observer
regions. -/
def le (U V : OperationalOverlapRegion) : Prop := U = V ∨ U = .bot

instance : DecidableRel le := fun U V =>
  decidable_of_iff (U = V ∨ U = .bot) Iff.rfl

/-- The declared meet: equal regions meet in themselves and the two distinct
observer regions meet in the bottom region. -/
def meet (U V : OperationalOverlapRegion) : OperationalOverlapRegion :=
  if U = V then U else .bot

theorem le_refl : ∀ U, le U U := by decide
theorem le_trans : ∀ U V W, le U V → le V W → le U W := by decide
theorem le_antisymm : ∀ U V, le U V → le V U → U = V := by decide
theorem meet_le_left : ∀ U V, le (meet U V) U := by decide
theorem meet_le_right : ∀ U V, le (meet U V) V := by decide
theorem le_meet : ∀ W U V, le W U → le W V → le W (meet U V) := by decide

/-- The region selected by each of the two observer labels. -/
def ofObserver (o : Fin 2) : OperationalOverlapRegion :=
  if o = 0 then .left else .right

theorem ofObserver_zero : ofObserver 0 = .left := rfl
theorem ofObserver_one : ofObserver 1 = .right := rfl
theorem left_ne_right : OperationalOverlapRegion.left ≠ .right := by decide
theorem bot_lt_left : le .bot .left ∧ ¬ le .left .bot := by decide
theorem bot_lt_right : le .bot .right ∧ ¬ le .right .bot := by decide

end OperationalOverlapRegion

/-- The common corner of an observer-indexed record.  In the good tower it
depends only on the record label.  In the mutated tower the second observer's
second record changes the shared corner from zero to one. -/
def operationalOverlapCorner (mutated : Bool) (o x : Fin 2) : ℂ :=
  if x = 0 then 1 else if mutated = true ∧ o = 1 then 1 else 0

/-- An observer-specific private diagonal entry.  This ensures the two good
record representatives differ before restriction even though their shared
corner agrees. -/
def operationalOverlapTail (o x : Fin 2) : ℂ :=
  ((2 * o.val + x.val + 1 : ℕ) : ℂ)

/-- The observer-indexed diagonal representative of a record label. -/
def operationalOverlapRecord (mutated : Bool) (o x : Fin 2) :
    Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.diagonal ![operationalOverlapCorner mutated o x,
    operationalOverlapTail o x]

theorem operationalOverlapRecord_mem_diagonal (mutated : Bool)
    (o x : Fin 2) :
    operationalOverlapRecord mutated o x ∈
      (diagonalPartition 2).publicSubalgebra :=
  diagonalPartition_mem_iff.mpr
    ⟨![operationalOverlapCorner mutated o x,
      operationalOverlapTail o x], rfl⟩

/-- Character restriction from the diagonal algebra to the scalar bottom
algebra, by evaluation at the common corner. -/
noncomputable def operationalOverlapCharacter :
    (diagonalPartition 2).publicSubalgebra →⋆ₐ[ℂ]
      (oneBlockPartition 2).publicSubalgebra where
  toFun X := ⟨(X : Matrix (Fin 2) (Fin 2) ℂ) 0 0 • 1,
    ((oneBlockPartition 2).publicSubalgebra).smul_mem (one_mem _) _⟩
  map_zero' := by
    apply Subtype.ext
    simp
  map_one' := by
    apply Subtype.ext
    simp
  map_add' X Y := by
    apply Subtype.ext
    simp [Matrix.add_apply, add_smul]
  map_mul' X Y := by
    apply Subtype.ext
    obtain ⟨c, hc⟩ := diagonalPartition_mem_iff.mp X.2
    obtain ⟨d, hd⟩ := diagonalPartition_mem_iff.mp Y.2
    simp only [MulMemClass.coe_mul]
    rw [← hc, ← hd]
    simp [smul_smul, Algebra.mul_smul_comm, mul_comm]
  commutes' c := by
    apply Subtype.ext
    simp [Algebra.algebraMap_eq_smul_one, Matrix.smul_apply]
  map_star' X := by
    apply Subtype.ext
    simp [star_smul]

/-- Regional algebras for the proper-meet witness: scalars at the bottom and
the full diagonal algebra on each observer region. -/
noncomputable def operationalOverlapRegionAlgebra :
    OperationalOverlapRegion →
      StarSubalgebra ℂ (Matrix (Fin 2) (Fin 2) ℂ)
  | .bot => (oneBlockPartition 2).publicSubalgebra
  | .left => (diagonalPartition 2).publicSubalgebra
  | .right => (diagonalPartition 2).publicSubalgebra

/-- The proper-meet restriction system: identities on equal regions and the
corner character on either drop to the bottom. -/
noncomputable def operationalOverlapRestrict :
    ∀ U V : OperationalOverlapRegion, OperationalOverlapRegion.le V U →
      operationalOverlapRegionAlgebra U →⋆ₐ[ℂ]
        operationalOverlapRegionAlgebra V
  | .bot, .bot, _ => StarAlgHom.id ℂ _
  | .left, .left, _ => StarAlgHom.id ℂ _
  | .right, .right, _ => StarAlgHom.id ℂ _
  | .left, .bot, _ => operationalOverlapCharacter
  | .right, .bot, _ => operationalOverlapCharacter
  | .bot, .left, h => absurd h (by decide)
  | .bot, .right, h => absurd h (by decide)
  | .left, .right, h => absurd h (by decide)
  | .right, .left, h => absurd h (by decide)

/-- The good and bad towers share all carriers and operations.  Only the
second observer's second public record is corner-mutated when `mutated=true`.
The nonzero commutator generator is inherited from the B17 witness. -/
noncomputable def operationalOverlapTower (mutated : Bool) :
    ConsensusTower Unit where
  indexNonempty := ⟨()⟩
  commonRefinement := fun _ _ => ⟨(), le_rfl, le_rfl⟩
  dim := fun _ => 2
  Observer := fun _ => Fin 2
  observerFinite := fun _ => inferInstance
  Record := fun _ => Fin 2
  recordFinite := fun _ => inferInstance
  recordOrder := fun _ _ => finRecordOrder 2
  publicAlgebra := fun _ _ => (diagonalPartition 2).publicSubalgebra
  public_mul_comm := fun _ _ X Y =>
    (diagonalPartition 2).publicSubalgebra_mul_comm X Y
  recordElement := fun _ o x =>
    ⟨operationalOverlapRecord mutated o x,
      operationalOverlapRecord_mem_diagonal mutated o x⟩
  state := fun _ _ => witnessState.1
  state_isState := fun _ _ => witnessState.2
  generator := fun _ _ => OPH.Tower.witnessGenerator
  observerRefine := fun _ o => o
  recordRefine := fun _ x => x
  algebraRefine := fun _ => StarAlgHom.id ℂ _
  observer_refine_refl := by intros; rfl
  observer_refine_trans := by intros; rfl
  record_refine_refl := by intros; rfl
  record_refine_trans := by intros; rfl
  algebra_refine_refl := by intros; rfl
  algebra_refine_trans := by intros; rfl
  record_order_natural := by intros; assumption
  public_natural := by intros; assumption
  record_element_natural := by intros; rfl
  state_natural := by intros; rfl
  generator_natural := by intros; rfl

/-- The proper-meet E6 net over either overlap tower. -/
noncomputable def operationalOverlapNet (mutated : Bool) :
    FiniteCausalObserverNet (operationalOverlapTower mutated) where
  Region := fun _ => OperationalOverlapRegion
  regionFintype := fun _ => inferInstance
  regionNonempty := fun _ => ⟨.bot⟩
  regionLE := fun _ => OperationalOverlapRegion.le
  regionLE_refl := fun _ U => OperationalOverlapRegion.le_refl U
  regionLE_trans := fun _ {U V W} hUV hVW =>
    OperationalOverlapRegion.le_trans U V W hUV hVW
  regionLE_antisymm := fun _ {U V} hUV hVU =>
    OperationalOverlapRegion.le_antisymm U V hUV hVU
  overlap := fun _ => OperationalOverlapRegion.meet
  overlap_le_left := fun _ U V => OperationalOverlapRegion.meet_le_left U V
  overlap_le_right := fun _ U V => OperationalOverlapRegion.meet_le_right U V
  le_overlap := fun _ {W U V} hWU hWV =>
    OperationalOverlapRegion.le_meet W U V hWU hWV
  disjoint := fun _ _ _ => False
  disjoint_symm := by intros; assumption
  disjoint_irrefl := by intro r U h; exact h
  localAlgebra := fun _ => operationalOverlapRegionAlgebra
  isotony := by
    intro r U V hUV
    rcases hUV with rfl | rfl
    · exact le_rfl
    · exact oneBlockPartition_publicSubalgebra_le _
  locality := by intros; contradiction
  restrict := fun _ {U V} h => operationalOverlapRestrict U V h
  restrict_refl := by intro r U X; cases U <;> rfl
  restrict_trans := by
    intro r U V W hVU hWV X
    cases U <;> cases V <;> cases W <;>
      first
        | rfl
        | exact absurd hVU (by decide)
        | exact absurd hWV (by decide)
  restrict_inclusion := by
    intro r U V hUV X
    cases U with
    | bot =>
        cases V with
        | bot => rfl
        | left =>
            apply Subtype.ext
            exact oneBlockPartition_corner_eval X (0 : Fin 2)
        | right =>
            apply Subtype.ext
            exact oneBlockPartition_corner_eval X (0 : Fin 2)
    | left =>
        cases V with
        | bot => exact absurd hUV (by decide)
        | left => rfl
        | right => exact absurd hUV (by decide)
    | right =>
        cases V with
        | bot => exact absurd hUV (by decide)
        | left => exact absurd hUV (by decide)
        | right => rfl
  regionRefine := fun _ U => U
  region_refine_refl := by intros; rfl
  region_refine_trans := by intros; rfl
  region_refine_mono := by intros; assumption
  overlap_natural := by intros; rfl
  disjoint_natural := by intros; assumption
  localAlgebra_natural := by intros; assumption
  repair := fun _ _ => LinearMap.id
  repair_idempotent := by intros; rfl
  repair_fixes_region := by intros; rfl
  repair_fixes_disjoint := by intros; rfl
  repair_natural := by intros; rfl

/-- The proper-meet access cut.  Each observer accesses the full diagonal
algebra on its own proper region. -/
noncomputable def operationalOverlapCut (mutated : Bool) :
    ObserverAccessCut (operationalOverlapTower mutated)
      (operationalOverlapNet mutated) () where
  observerRegion := OperationalOverlapRegion.ofObserver
  accessibleAlgebra := fun _ => (diagonalPartition 2).publicSubalgebra
  accessible_le_region := by intro o; fin_cases o <;> exact le_rfl
  public_le_accessible := fun _ => le_rfl
  restrict_preserves := by
    intro o U V hU hVU X hX
    fin_cases o <;>
      exact (operationalOverlapNet mutated).isotony ()
        ((operationalOverlapNet mutated).regionLE_trans () hVU hU)
        ((operationalOverlapNet mutated).restrict () hVU X).2

/-- The record-advancing prediction on the two-record chain. -/
def operationalOverlapPredict : Fin 2 → Fin 2 := fun _ => 1

/-- One of the two operational observers over either tower.  Its readback is
the diagonal partition expectation, its generator is the nonzero diagonal
Hamiltonian commutator, its control is nonidentity diagonal conjugation, and
its prediction advances to the final record. -/
noncomputable def operationalOverlapObserver (mutated : Bool) (o : Fin 2) :
    OperationalObserver (operationalOverlapTower mutated) where
  label := fun _ => o
  label_refine := by intros; rfl
  accessible := fun _ => (diagonalPartition 2).publicSubalgebra
  public_le_accessible := fun _ => le_rfl
  readback := fun _ X => partitionAverage (diagonalPartition 2) X
  readback_mem_public := fun _ X _ =>
    partitionAverage_mem_span (diagonalPartition 2) X
  readback_fixes_public := fun _ X hX =>
    partitionAverage_fixes (diagonalPartition 2) hX
  record_durable := by
    intro r x
    exact OPH.Tower.witnessGenerator_diagonal
      ![operationalOverlapCorner mutated o x,
        operationalOverlapTail o x]
  control := fun _ _ => OPH.Tower.witnessControl
  control_preserves_accessible := by
    intro r x X hX
    obtain ⟨d, hd⟩ := diagonalPartition_mem_iff.mp hX
    rw [← hd]
    have hfix : OPH.Tower.witnessControl (Matrix.diagonal d) =
        Matrix.diagonal d :=
      OPH.Tower.witnessControl_diagonal d
    have hmem : Matrix.diagonal d ∈
        (diagonalPartition 2).publicSubalgebra :=
      diagonalPartition_mem_iff.mpr ⟨d, rfl⟩
    exact hfix.symm ▸ hmem
  predict := fun _ => operationalOverlapPredict
  predict_no_regress := by
    intro r
    show ∀ x : Fin 2, ¬ operationalOverlapPredict x < x
    decide
  predict_refine := by intros; rfl
  readout := fun r x =>
    ((operationalOverlapTower mutated).state r o *
      ((operationalOverlapTower mutated).recordElement r o x).1).trace
  predict_readout_congruent := by intros; rfl
  readout_eq_pairing := by intros; rfl
  record_fixed_by_control := by
    intro r x y
    exact OPH.Tower.witnessControl_diagonal
      ![operationalOverlapCorner mutated o x,
        operationalOverlapTail o x]

/-- The two good observer-specific representatives of the first record are
different before restriction.  Their agreement below is therefore caused by
the overlap characters rather than by reusing one ambient matrix. -/
theorem operationalOverlap_good_records_differ_before_restriction :
    ((operationalOverlapTower false).recordElement () (0 : Fin 2)
        (0 : Fin 2)).1 ≠
      ((operationalOverlapTower false).recordElement () (1 : Fin 2)
        (0 : Fin 2)).1 := by
  change operationalOverlapRecord false 0 0 ≠
    operationalOverlapRecord false 1 0
  intro h
  have h11 := congrFun (congrFun h (1 : Fin 2)) (1 : Fin 2)
  norm_num [operationalOverlapRecord,
    operationalOverlapCorner, operationalOverlapTail] at h11

/-- The two observer regions in the good witness are distinct and their
declared meet is the bottom region, proper below each. -/
theorem operationalOverlap_good_regions_proper :
    OperationalOverlapRegion.ofObserver (0 : Fin 2) ≠
        OperationalOverlapRegion.ofObserver (1 : Fin 2) ∧
      OperationalOverlapRegion.meet
          (OperationalOverlapRegion.ofObserver (0 : Fin 2))
          (OperationalOverlapRegion.ofObserver (1 : Fin 2)) =
        OperationalOverlapRegion.bot ∧
      (¬ OperationalOverlapRegion.le OperationalOverlapRegion.left
        OperationalOverlapRegion.bot) ∧
      (¬ OperationalOverlapRegion.le OperationalOverlapRegion.right
        OperationalOverlapRegion.bot) := by
  decide

/-- In the good tower, every common record label has equal typed restrictions
from the two distinct owner regions. -/
theorem operationalOverlap_good_restrictions_agree (x : Fin 2) :
    (operationalOverlapCut false).recordRestrictionLeft
        (0 : Fin 2) (1 : Fin 2) x =
      (operationalOverlapCut false).recordRestrictionRight
        (0 : Fin 2) (1 : Fin 2) x := by
  apply Subtype.ext
  change
    operationalOverlapCorner false 0 x •
        (1 : Matrix (Fin 2) (Fin 2) ℂ) =
      operationalOverlapCorner false 1 x •
        (1 : Matrix (Fin 2) (Fin 2) ℂ)
  fin_cases x <;> rfl

/-- The first good restricted record is the nonzero scalar sure event on the
bottom overlap. -/
theorem operationalOverlap_good_restriction_nonzero :
    (operationalOverlapCut false).recordRestrictionLeft
        (0 : Fin 2) (1 : Fin 2) (0 : Fin 2) ≠ 0 := by
  intro h
  have hv := congrArg Subtype.val h
  change
    operationalOverlapCorner false 0 0 •
        (1 : Matrix (Fin 2) (Fin 2) ℂ) = 0 at hv
  have h00 := congrFun (congrFun hv (0 : Fin 2)) (0 : Fin 2)
  norm_num [operationalOverlapCorner, Matrix.smul_apply,
    Matrix.one_apply_eq] at h00

/-- The good observers assign the same own readout to both common record
labels. -/
theorem operationalOverlap_good_readouts_agree (x : Fin 2) :
    (operationalOverlapObserver false (0 : Fin 2)).readout () x =
      (operationalOverlapObserver false (1 : Fin 2)).readout () x := by
  change (witnessState.1 * operationalOverlapRecord false 0 x).trace =
    (witnessState.1 * operationalOverlapRecord false 1 x).trace
  fin_cases x <;>
    norm_num [operationalOverlapRecord, operationalOverlapCorner,
      operationalOverlapTail, witnessState, Matrix.trace_fin_two,
      Matrix.mul_apply, Matrix.diagonal_apply, Fin.sum_univ_two]

/-- The full nonvacuous cross-observer receipt on the good proper-meet
witness. -/
theorem operationalOverlap_good_evidence :
    OperationalOverlapEvidence (operationalOverlapCut false)
      (operationalOverlapObserver false (0 : Fin 2))
      (operationalOverlapObserver false (1 : Fin 2)) := by
  constructor
  · show (0 : Fin 2) ≠ 1
    decide
  · rfl
  · rfl
  · exact operationalOverlap_good_restrictions_agree
  · exact operationalOverlap_good_readouts_agree
  · change ∃ x : Fin 2,
      (operationalOverlapCut false).recordRestrictionLeft
        (0 : Fin 2) (1 : Fin 2) x ≠ 0
    exact ⟨0, operationalOverlap_good_restriction_nonzero⟩

/-- Evaluated common-access theorem: every good common overlap record belongs
to both operational accessible algebras. -/
theorem operationalOverlap_good_common_mem_both (x : Fin 2) :
    ((operationalOverlapCut false).recordRestrictionLeft
        (0 : Fin 2) (1 : Fin 2) x).1 ∈
        (operationalOverlapObserver false (0 : Fin 2)).accessible () ∧
      ((operationalOverlapCut false).recordRestrictionLeft
        (0 : Fin 2) (1 : Fin 2) x).1 ∈
        (operationalOverlapObserver false (1 : Fin 2)).accessible () :=
  operationalOverlap_good_evidence.commonRecord_mem_both x

/-- Corner-mutation negative control: on the otherwise identical bad tower,
the second committed record restricts to zero from the left and to one from
the right. -/
theorem operationalOverlap_bad_restrictions_disagree :
    (operationalOverlapCut true).recordRestrictionLeft
        (0 : Fin 2) (1 : Fin 2) (1 : Fin 2) ≠
      (operationalOverlapCut true).recordRestrictionRight
        (0 : Fin 2) (1 : Fin 2) (1 : Fin 2) := by
  intro h
  have hv := congrArg Subtype.val h
  change
    operationalOverlapCorner true 0 1 •
        (1 : Matrix (Fin 2) (Fin 2) ℂ) =
      operationalOverlapCorner true 1 1 •
        (1 : Matrix (Fin 2) (Fin 2) ℂ) at hv
  have h00 := congrFun (congrFun hv (0 : Fin 2)) (0 : Fin 2)
  norm_num [operationalOverlapCorner, Matrix.smul_apply,
    Matrix.one_apply_eq] at h00

/-- The same corner mutation also makes the two own readouts of the second
record disagree. -/
theorem operationalOverlap_bad_readouts_disagree :
    (operationalOverlapObserver true (0 : Fin 2)).readout () (1 : Fin 2) ≠
      (operationalOverlapObserver true (1 : Fin 2)).readout () (1 : Fin 2) := by
  change (witnessState.1 * operationalOverlapRecord true 0 1).trace ≠
    (witnessState.1 * operationalOverlapRecord true 1 1).trace
  norm_num [operationalOverlapRecord, operationalOverlapCorner,
    operationalOverlapTail, witnessState, Matrix.trace_fin_two,
    Matrix.mul_apply, Matrix.diagonal_apply, Fin.sum_univ_two]

/-- Exact negative receipt: the corner-mutated pair cannot inhabit the
cross-observer evidence structure. -/
theorem operationalOverlap_bad_not_evidence :
    ¬ OperationalOverlapEvidence (operationalOverlapCut true)
      (operationalOverlapObserver true (0 : Fin 2))
      (operationalOverlapObserver true (1 : Fin 2)) := by
  intro E
  exact operationalOverlap_bad_restrictions_disagree
    (E.record_restriction_eq (1 : Fin 2))

end OverlapWitness

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.OperationalOverlapEvidence.commonRecord_mem_both
#print axioms OPH.QFT.FiniteEventGeometryWorld.operationalObservers_share_visible_event_record
#print axioms OPH.QFT.operationalOverlap_good_records_differ_before_restriction
#print axioms OPH.QFT.operationalOverlap_good_regions_proper
#print axioms OPH.QFT.operationalOverlap_good_restrictions_agree
#print axioms OPH.QFT.operationalOverlap_good_restriction_nonzero
#print axioms OPH.QFT.operationalOverlap_good_readouts_agree
#print axioms OPH.QFT.operationalOverlap_good_evidence
#print axioms OPH.QFT.operationalOverlap_good_common_mem_both
#print axioms OPH.QFT.operationalOverlap_bad_restrictions_disagree
#print axioms OPH.QFT.operationalOverlap_bad_readouts_disagree
#print axioms OPH.QFT.operationalOverlap_bad_not_evidence

end OPH.QFT
