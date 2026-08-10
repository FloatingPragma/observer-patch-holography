import QFT.RichFibreWitness
import QFT.CoverageReceipt

/-!
# The rich-fibre regional net: four noncommutative source regions

This module builds the finite causal observer net over the rich-fibre
witness of `QFT/RichFibreWitness.lean` (the second preregistered run
`e1_prereg2_64k_20260810`; every preregistered gate passed).  Each of the
four selected observers realizes at least three split fibres on its
twenty-node window, so the uniform enrichment rule of
`QFT/SourceRegionalNet.lean` yields a genuinely noncommutative regional
algebra at **every** observer region.  The windows are pairwise disjoint
in the source (`richSupport_pairwise_disjoint`, the preregistered gate G3
on the literals), so the ambient block model mirrors four separated
physical supports.

Contents: the block-relation machinery over the `Fin 81` ambient basis;
the enrichment rule from the payload literals (two ambient indices are
related exactly when they lie in one realized split fibre of one
observer); the region lattice; the block regional algebras with isotony,
locality on declared-disjoint regions (a genuine block-versus-block
commutation theorem, since two enriched regions are now disjoint), and
the character/block restriction system; the tower and the
`FiniteCausalObserverNet` instance; and the receipts (noncommutativity at
all four regions, the realized two-by-two matrix factor of each
designated split fibre, region separation, and the computed coverage law
for the four-window family under the top region).

**Boundary.**  The net mirrors the committed payload of one preregistered
bounded run; the enrichment rule, region lattice, restriction system, and
state are declared readings of those literals through the committed E1
interfaces.  No physical causality, locality, clock, instrument, or
prediction is claimed; source-produced CP/CPTP instrument provenance
remains with E2 and continuum structure with E3.
-/

namespace OPH.QFT

open OPH.Tower
open Matrix

set_option maxRecDepth 65536

/-! ## Ambient owner facts and collapse maps -/

/-- The base point is the one ambient index owned by no observer. -/
theorem richNodeOwner_none_iff : ∀ i : Fin 81,
    richNodeOwner i = none ↔ i = 0 := by decide

/-- The collapse of the ambient basis onto one observer window: fixed on
the window, the base point elsewhere. -/
def richCollapse (o : RichObserver) : Fin 81 → Fin 81 := fun i =>
  if richNodeOwner i = some o then i else 0

theorem richCollapse_of_owner {o : RichObserver} {i : Fin 81}
    (h : richNodeOwner i = some o) : richCollapse o i = i := by
  simp [richCollapse, h]

theorem richCollapse_of_not_owner {o : RichObserver} {i : Fin 81}
    (h : ¬ richNodeOwner i = some o) : richCollapse o i = 0 := by
  simp [richCollapse, h]

/-! ## Class maps for the tower's commutative layers -/

/-- The bottom class map: one class. -/
def richBot : Fin 81 → Unit := fun _ => ()

/-- The public record class map of one observer: the record class on the
window, one ambient class elsewhere. -/
def richPubMap (o : RichObserver) : Fin 81 → Option (Fin 32) := fun i =>
  if richNodeOwner i = some o then some (richNodeRecordClass i) else none

/-! ## Diagonal spans over the ambient basis -/

theorem richDiagonal_const_eq_algebraMap (z : ℂ) :
    Matrix.diagonal (fun _ : Fin 81 => z) =
      algebraMap ℂ (Matrix (Fin 81) (Fin 81) ℂ) z := by
  rw [Matrix.algebraMap_eq_diagonal]
  rfl

/-- The diagonal span induced by a class function `f` on the ambient
basis. -/
def richDiagSpan {α : Type} (f : Fin 81 → α) :
    StarSubalgebra ℂ (Matrix (Fin 81) (Fin 81) ℂ) where
  carrier := {X | ∃ c : α → ℂ, Matrix.diagonal (fun i => c (f i)) = X}
  zero_mem' := ⟨fun _ => 0, by simp⟩
  one_mem' := ⟨fun _ => 1, by simp⟩
  add_mem' := by
    rintro X Y ⟨c, rfl⟩ ⟨d, rfl⟩
    exact ⟨fun a => c a + d a, by simp [Matrix.diagonal_add]⟩
  mul_mem' := by
    rintro X Y ⟨c, rfl⟩ ⟨d, rfl⟩
    exact ⟨fun a => c a * d a, by simp [Matrix.diagonal_mul_diagonal]⟩
  algebraMap_mem' := fun z => ⟨fun _ => z, richDiagonal_const_eq_algebraMap z⟩
  star_mem' := by
    rintro X ⟨c, rfl⟩
    exact ⟨fun a => star (c a), by
      simp [Matrix.star_eq_conjTranspose, Matrix.diagonal_conjTranspose]⟩

theorem mem_richDiagSpan_iff {α : Type} {f : Fin 81 → α}
    {X : Matrix (Fin 81) (Fin 81) ℂ} :
    X ∈ richDiagSpan f ↔
      ∃ c : α → ℂ, Matrix.diagonal (fun i => c (f i)) = X :=
  Iff.rfl

/-- Every diagonal span is commutative. -/
theorem richDiagSpan_mul_comm {α : Type} (f : Fin 81 → α)
    (X Y : richDiagSpan f) : X * Y = Y * X := by
  apply Subtype.ext
  obtain ⟨c, hc⟩ := X.2
  obtain ⟨d, hd⟩ := Y.2
  show (X : Matrix (Fin 81) (Fin 81) ℂ) * Y =
    (Y : Matrix (Fin 81) (Fin 81) ℂ) * X
  rw [← hc, ← hd, Matrix.diagonal_mul_diagonal,
    Matrix.diagonal_mul_diagonal]
  congr 1
  funext i
  exact mul_comm _ _

/-! ## The enrichment rule from the payload literals -/

/-- The split marker: `true` exactly on the payload's listed split
record classes of the observer. -/
def richSplitRecord (o : RichObserver) (r : Fin 32) : Bool :=
  (RichObserver.splitFibreClasses o).contains (r : ℕ)

/-- The first ambient index of one observer's window carrying a given
record class; the base point when the class is absent. -/
def richFirstOfClass (o : RichObserver) (c : Fin 32) : Fin 81 :=
  match (List.finRange 81).find? (fun k =>
      decide (richNodeOwner k = some o) &&
        decide (richNodeRecordClass k = c)) with
  | some k => k
  | none => 0

/-- The representative of an ambient index under one observer's
enrichment: indices of one realized split fibre share the fibre's first
window index; everything else represents itself. -/
def richRep (o : RichObserver) : Fin 81 → Fin 81 := fun i =>
  if richNodeOwner i = some o ∧
      richSplitRecord o (richNodeRecordClass i) = true then
    richFirstOfClass o (richNodeRecordClass i)
  else i

/-- The enrichment relation of one observer: the kernel of its
representative map. -/
def richRelObs (o : RichObserver) : Fin 81 → Fin 81 → Bool :=
  fun i j => decide (richRep o i = richRep o j)

/-- The top representative: the owner's representative on every owned
index, identity at the base point. -/
def richRepTop : Fin 81 → Fin 81 := fun i =>
  match richNodeOwner i with
  | none => i
  | some o => richRep o i

/-- The top enrichment relation: the kernel of the top representative. -/
def richRelTop : Fin 81 → Fin 81 → Bool :=
  fun i j => decide (richRepTop i = richRepTop j)

/-- The bottom relation: equality. -/
def richRelBot : Fin 81 → Fin 81 → Bool := fun i j => decide (i = j)

/-! ## Block relations and their subalgebras -/

private theorem bool_eq_false' {b : Bool} (h : ¬ b = true) : b = false := by
  cases b
  · rfl
  · exact absurd rfl h

private theorem ne_true_of_eq_false' {b : Bool} (h : b = false) :
    ¬ b = true := fun h2 => by
  rw [h] at h2
  exact Bool.noConfusion h2

/-- A Boolean block relation on the ambient basis, with a marker of
ambient indices whose diagonal must be constant. -/
structure RichBlockRelation where
  rel : Fin 81 → Fin 81 → Bool
  amb : Fin 81 → Bool
  rel_refl : ∀ i, rel i i = true
  rel_symm : ∀ i j, rel i j = rel j i
  rel_trans : ∀ i j k, rel i j = true → rel j k = true → rel i k = true
  amb_singleton : ∀ i j, amb i = true → rel i j = true → i = j

namespace RichBlockRelation

/-- Build a block relation from a representative-map characterization. -/
def ofRep (rel : Fin 81 → Fin 81 → Bool) (amb : Fin 81 → Bool)
    (f : Fin 81 → Fin 81) (hchar : ∀ i j, rel i j = decide (f i = f j))
    (hamb : ∀ i j, amb i = true → rel i j = true → i = j) :
    RichBlockRelation where
  rel := rel
  amb := amb
  rel_refl i := by rw [hchar]; exact decide_eq_true rfl
  rel_symm i j := by
    rw [hchar, hchar]
    exact decide_eq_decide.mpr ⟨Eq.symm, Eq.symm⟩
  rel_trans i j k h1 h2 := by
    rw [hchar] at h1 h2 ⊢
    exact decide_eq_true ((of_decide_eq_true h1).trans (of_decide_eq_true h2))
  amb_singleton := hamb

/-- The block subalgebra of a block relation: matrices supported on
related pairs, constant on the diagonal across ambient indices. -/
def subalgebra (B : RichBlockRelation) :
    StarSubalgebra ℂ (Matrix (Fin 81) (Fin 81) ℂ) where
  carrier := {X | (∀ i j, B.rel i j = false → X i j = 0) ∧
    (∀ i j, B.amb i = true → B.amb j = true → X i i = X j j)}
  zero_mem' := ⟨fun _ _ _ => rfl, fun _ _ _ _ => rfl⟩
  one_mem' := by
    refine ⟨fun i j h => ?_, fun i j _ _ => ?_⟩
    · have hij : i ≠ j := fun he => by
        subst he
        exact Bool.noConfusion ((B.rel_refl i).symm.trans h)
      exact Matrix.one_apply_ne hij
    · rw [Matrix.one_apply_eq, Matrix.one_apply_eq]
  add_mem' := by
    rintro X Y ⟨hX0, hXs⟩ ⟨hY0, hYs⟩
    refine ⟨fun i j h => ?_, fun i j hi hj => ?_⟩
    · rw [Matrix.add_apply, hX0 i j h, hY0 i j h, add_zero]
    · rw [Matrix.add_apply, Matrix.add_apply, hXs i j hi hj, hYs i j hi hj]
  mul_mem' := by
    rintro X Y ⟨hX0, hXs⟩ ⟨hY0, hYs⟩
    refine ⟨fun i j h => ?_, fun i j hi hj => ?_⟩
    · rw [Matrix.mul_apply]
      refine Finset.sum_eq_zero fun k _ => ?_
      by_cases hik : B.rel i k = true
      · have hkj : B.rel k j = false := by
          cases hkj' : B.rel k j
          · rfl
          · exact absurd (B.rel_trans i k j hik hkj')
              (fun htrue => Bool.noConfusion (h.symm.trans htrue))
        rw [hY0 k j hkj, mul_zero]
      · rw [hX0 i k (bool_eq_false' hik), zero_mul]
    · rw [Matrix.mul_apply, Matrix.mul_apply]
      rw [Finset.sum_eq_single_of_mem i (Finset.mem_univ i)
          (fun k _ hk => by
            have h0 : B.rel i k = false := by
              cases hik : B.rel i k
              · rfl
              · exact absurd (B.amb_singleton i k hi hik).symm hk
            rw [hX0 i k h0, zero_mul]),
        Finset.sum_eq_single_of_mem j (Finset.mem_univ j)
          (fun k _ hk => by
            have h0 : B.rel j k = false := by
              cases hjk : B.rel j k
              · rfl
              · exact absurd (B.amb_singleton j k hj hjk).symm hk
            rw [hX0 j k h0, zero_mul]),
        hXs i j hi hj, hYs i j hi hj]
  algebraMap_mem' := by
    intro z
    refine ⟨fun i j h => ?_, fun i j _ _ => ?_⟩
    · have hij : i ≠ j := fun he => by
        subst he
        exact Bool.noConfusion ((B.rel_refl i).symm.trans h)
      rw [← richDiagonal_const_eq_algebraMap, Matrix.diagonal_apply_ne _ hij]
    · rw [← richDiagonal_const_eq_algebraMap, Matrix.diagonal_apply_eq,
        Matrix.diagonal_apply_eq]
  star_mem' := by
    rintro X ⟨hX0, hXs⟩
    refine ⟨fun i j h => ?_, fun i j hi hj => ?_⟩
    · rw [Matrix.star_apply, hX0 j i ((B.rel_symm j i).trans h), star_zero]
    · rw [Matrix.star_apply, Matrix.star_apply, hXs i j hi hj]

theorem mem_subalgebra_iff (B : RichBlockRelation)
    {X : Matrix (Fin 81) (Fin 81) ℂ} :
    X ∈ B.subalgebra ↔
      (∀ i j, B.rel i j = false → X i j = 0) ∧
        (∀ i j, B.amb i = true → B.amb j = true → X i i = X j j) :=
  Iff.rfl

end RichBlockRelation

/-! ## The concrete block relations and kernel facts -/

/-- The ambient-singleton law of one observer's block relation. -/
theorem richObsBlock_amb : ∀ (o : RichObserver) (i j : Fin 81),
    decide (¬ richNodeOwner i = some o) = true →
      richRelObs o i j = true → i = j := by decide

/-- The block relation of one observer region: the enrichment relation
on the window, constant diagonal off the window. -/
def richObsBlock (o : RichObserver) : RichBlockRelation :=
  .ofRep (richRelObs o) (fun i => decide (¬ richNodeOwner i = some o))
    (richRep o) (fun _ _ => rfl) (richObsBlock_amb o)

/-- The block relation of the top region: no ambient constraint. -/
def richTopBlock : RichBlockRelation :=
  .ofRep richRelTop (fun _ => false) richRepTop (fun _ _ => rfl)
    (fun _ _ h _ => Bool.noConfusion h)

/-- The block relation of the bottom region: equality with a global
ambient constraint, so its subalgebra is the scalars. -/
def richBotBlock : RichBlockRelation :=
  .ofRep richRelBot (fun _ => true) id (fun _ _ => rfl)
    (fun _ _ _ h => of_decide_eq_true h)

/-- Off-diagonal relatedness under one observer forces double
ownership. -/
theorem richRelObs_offdiag_owner : ∀ (o : RichObserver) (i j : Fin 81),
    richRelObs o i j = true → ¬ i = j →
      richNodeOwner i = some o ∧ richNodeOwner j = some o := by decide

/-- Relatedness under one observer is closed over its ownership. -/
theorem richRelObs_owner_closed : ∀ (o : RichObserver) (i k : Fin 81),
    richRelObs o i k = true → richNodeOwner i = some o →
      richNodeOwner k = some o := by decide

/-- Off the window, one observer's relation is the identity. -/
theorem richRelObs_not_owner_singleton : ∀ (o : RichObserver) (i j : Fin 81),
    ¬ richNodeOwner i = some o → richRelObs o i j = true → i = j := by
  decide

/-- The base point is a singleton of every observer relation. -/
theorem richRelObs_zero_singleton : ∀ (o : RichObserver) (k : Fin 81),
    richRelObs o 0 k = true → k = 0 := by decide

/-- The base point is a singleton of the top relation. -/
theorem richRelTop_zero_singleton : ∀ k : Fin 81,
    richRelTop 0 k = true → k = 0 := by decide

/-- On an owned index, the top relation restricts to the owner's
relation. -/
theorem richRelTop_of_owner : ∀ (o : RichObserver) (i j : Fin 81),
    richNodeOwner i = some o → richRelTop i j = richRelObs o i j := by
  decide

/-- Each observer's relation refines the top relation. -/
theorem richRelObs_le_top : ∀ (o : RichObserver) (i j : Fin 81),
    richRelObs o i j = true → richRelTop i j = true := by decide

/-- Off-diagonal top relatedness happens away from the base point. -/
theorem richRelTop_offdiag_owner : ∀ i j : Fin 81,
    richRelTop i j = true → ¬ i = j → richNodeOwner i ≠ none := by decide

/-! ## The region lattice -/

/-- Region labels: bottom, the four windows, top. -/
inductive RichRegion : Type
  | bot
  | r86
  | r88
  | r247
  | r384
  | top
  deriving DecidableEq, Fintype

namespace RichRegion

/-- Region inclusion: bottom below everything, everything below top,
windows pairwise incomparable. -/
def le : RichRegion → RichRegion → Prop
  | _, .top => True
  | .bot, _ => True
  | .r86, .r86 => True
  | .r88, .r88 => True
  | .r247, .r247 => True
  | .r384, .r384 => True
  | _, _ => False

instance : (U V : RichRegion) → Decidable (le U V) := fun U V => by
  cases U <;> cases V <;>
    first
      | exact (inferInstance : Decidable True)
      | exact (inferInstance : Decidable False)

theorem le_refl : ∀ U, le U U := by decide

theorem le_trans : ∀ U V W, le U V → le V W → le U W := by decide

theorem le_antisymm : ∀ U V, le U V → le V U → U = V := by decide

/-- The declared overlap: the meet in the region lattice. -/
def meet (U V : RichRegion) : RichRegion :=
  if le U V then U else if le V U then V else .bot

theorem meet_le_left : ∀ U V, le (meet U V) U := by decide

theorem meet_le_right : ∀ U V, le (meet U V) V := by decide

theorem le_meet : ∀ W U V, le W U → le W V → le W (meet U V) := by decide

/-- The support region of each observer. -/
def regionOf : RichObserver → RichRegion
  | .o86 => .r86
  | .o88 => .r88
  | .o247 => .r247
  | .o384 => .r384

/-- Declared disjointness: two distinct support regions.  The literal
receipt `richSupport_pairwise_disjoint` (the preregistered gate G3)
backs this declaration in the source. -/
def disj (U V : RichRegion) : Prop :=
  U ≠ V ∧ U ≠ .bot ∧ U ≠ .top ∧ V ≠ .bot ∧ V ≠ .top

end RichRegion

/-- The block relation of a region. -/
def richRegionBlock : RichRegion → RichBlockRelation
  | .bot => richBotBlock
  | .r86 => richObsBlock .o86
  | .r88 => richObsBlock .o88
  | .r247 => richObsBlock .o247
  | .r384 => richObsBlock .o384
  | .top => richTopBlock

/-- The regional algebra assignment: the block subalgebra of the
region's relation.  Scalars at the bottom, one enriched block algebra
per window, the block-diagonal join at the top. -/
def richRegionAlgebra (U : RichRegion) :
    StarSubalgebra ℂ (Matrix (Fin 81) (Fin 81) ℂ) :=
  (richRegionBlock U).subalgebra

theorem richRegionAlgebra_regionOf (o : RichObserver) :
    richRegionAlgebra (RichRegion.regionOf o) =
      (richObsBlock o).subalgebra := by
  cases o <;> rfl

/-! ## Isotony -/

/-- Bottom members are scalar diagonals. -/
theorem richBot_mem_iff {X : Matrix (Fin 81) (Fin 81) ℂ} :
    X ∈ richRegionAlgebra .bot ↔
      X = Matrix.diagonal (fun _ => X 0 0) := by
  constructor
  · rintro ⟨h0, hc⟩
    ext i j
    by_cases hij : i = j
    · subst hij
      rw [Matrix.diagonal_apply_eq]
      exact hc i 0 rfl rfl
    · rw [Matrix.diagonal_apply_ne _ hij]
      exact h0 i j (decide_eq_false hij)
  · intro h
    refine ⟨fun i j hrel => ?_, fun i j _ _ => ?_⟩
    · have hij : i ≠ j := fun he => by
        subst he
        exact Bool.noConfusion
          ((richBotBlock.rel_refl i).symm.trans hrel)
      conv_lhs => rw [h]
      rw [Matrix.diagonal_apply_ne _ hij]
    · conv_lhs => rw [h]
      conv_rhs => rw [h]
      rw [Matrix.diagonal_apply_eq, Matrix.diagonal_apply_eq]

theorem richRegionAlgebra_bot_le (U : RichRegion) :
    richRegionAlgebra .bot ≤ richRegionAlgebra U := by
  intro X hX
  have hX' := richBot_mem_iff.mp hX
  refine ⟨fun i j hrel => ?_, fun i j _ _ => ?_⟩
  · have hij : i ≠ j := fun he => by
      subst he
      exact Bool.noConfusion
        (((richRegionBlock U).rel_refl i).symm.trans hrel)
    conv_lhs => rw [hX']
    rw [Matrix.diagonal_apply_ne _ hij]
  · conv_lhs => rw [hX']
    conv_rhs => rw [hX']
    rw [Matrix.diagonal_apply_eq, Matrix.diagonal_apply_eq]

/-- An observer relation is false wherever the top relation is false. -/
theorem richRelObs_false_of_top_false (o : RichObserver) (i j : Fin 81)
    (hrel : richRelTop i j = false) : richRelObs o i j = false := by
  cases hobs : richRelObs o i j
  · rfl
  · exact absurd (richRelObs_le_top o i j hobs)
      (ne_true_of_eq_false' hrel)

theorem richRegionAlgebra_le_top (U : RichRegion) :
    richRegionAlgebra U ≤ richRegionAlgebra .top := by
  cases U
  case bot => exact richRegionAlgebra_bot_le .top
  case top => exact le_rfl
  all_goals
    rintro X ⟨h0, _⟩
    exact ⟨fun i j hrel =>
        h0 i j (richRelObs_false_of_top_false _ i j hrel),
      fun i j hi _ => Bool.noConfusion hi⟩

/-- Isotony of the regional assignment along region inclusion. -/
theorem richRegionAlgebra_isotony {U V : RichRegion}
    (hUV : RichRegion.le U V) :
    richRegionAlgebra U ≤ richRegionAlgebra V := by
  cases U <;> cases V <;>
    first
      | exact le_rfl
      | exact absurd hUV (by decide)
      | exact richRegionAlgebra_bot_le _
      | exact richRegionAlgebra_le_top _

/-! ## Locality: two disjoint enriched regions commute -/

/-- Structure of a block member: an off-diagonal nonzero entry lies
inside the observer's window. -/
theorem richBlock_offdiag_window {o : RichObserver}
    {X : Matrix (Fin 81) (Fin 81) ℂ}
    (hX : X ∈ (richObsBlock o).subalgebra) {i j : Fin 81}
    (hne : ¬ i = j) (hval : X i j ≠ 0) :
    richNodeOwner i = some o ∧ richNodeOwner j = some o := by
  by_cases hrel : richRelObs o i j = true
  · exact richRelObs_offdiag_owner o i j hrel hne
  · exact absurd (hX.1 i j (bool_eq_false' hrel)) hval

/-- Structure of a block member: the diagonal is the base-point scalar
off the observer's window. -/
theorem richBlock_offwindow_diag {o : RichObserver}
    {X : Matrix (Fin 81) (Fin 81) ℂ}
    (hX : X ∈ (richObsBlock o).subalgebra) {i : Fin 81}
    (hi : ¬ richNodeOwner i = some o) : X i i = X 0 0 := by
  have h0 : ¬ richNodeOwner (0 : Fin 81) = some o := by
    rw [(richNodeOwner_none_iff (0 : Fin 81)).mpr rfl]
    simp
  exact hX.2 i 0 (decide_eq_true hi) (decide_eq_true h0)

/-- Members of two block algebras over distinct observers commute: the
windows are separated, each block sees the other's diagonal as a
constant, and cross-window products vanish. -/
theorem richBlock_commute_disjoint {o o' : RichObserver} (hne : o ≠ o')
    {X Y : Matrix (Fin 81) (Fin 81) ℂ}
    (hX : X ∈ (richObsBlock o).subalgebra)
    (hY : Y ∈ (richObsBlock o').subalgebra) :
    Commute X Y := by
  have howner : ∀ i : Fin 81,
      ¬ (richNodeOwner i = some o ∧ richNodeOwner i = some o') := by
    rintro i ⟨h1, h2⟩
    rw [h1] at h2
    exact hne (Option.some.inj h2)
  show X * Y = Y * X
  ext i j
  rw [Matrix.mul_apply, Matrix.mul_apply]
  by_cases hij : i = j
  · -- diagonal entry: only the k = i term survives on both sides
    subst hij
    rw [Finset.sum_eq_single_of_mem i (Finset.mem_univ i)
        (fun k _ hk => by
          by_cases hXik : X i k = 0
          · rw [hXik, zero_mul]
          · obtain ⟨hio, hko⟩ :=
              richBlock_offdiag_window hX (fun h => hk h.symm) hXik
            by_cases hYki : Y k i = 0
            · rw [hYki, mul_zero]
            · obtain ⟨hko', _⟩ :=
                richBlock_offdiag_window hY hk hYki
              exact absurd ⟨hko, hko'⟩ (howner k)),
      Finset.sum_eq_single_of_mem i (Finset.mem_univ i)
        (fun k _ hk => by
          by_cases hYik : Y i k = 0
          · rw [hYik, zero_mul]
          · obtain ⟨hio', hko'⟩ :=
              richBlock_offdiag_window hY (fun h => hk h.symm) hYik
            by_cases hXki : X k i = 0
            · rw [hXki, mul_zero]
            · obtain ⟨hko, _⟩ :=
                richBlock_offdiag_window hX hk hXki
              exact absurd ⟨hko, hko'⟩ (howner k))]
    exact mul_comm _ _
  · by_cases hwin : richNodeOwner i = some o ∧ richNodeOwner j = some o
    · -- both endpoints in o's window: Y contributes its scalar
      have hYii : Y i i = Y 0 0 :=
        richBlock_offwindow_diag hY (fun h => howner i ⟨hwin.1, h⟩)
      have hYjj : Y j j = Y 0 0 :=
        richBlock_offwindow_diag hY (fun h => howner j ⟨hwin.2, h⟩)
      rw [Finset.sum_eq_single_of_mem j (Finset.mem_univ j)
          (fun k _ hk => by
            by_cases hYkj : Y k j = 0
            · rw [hYkj, mul_zero]
            · obtain ⟨_, hjo'⟩ := richBlock_offdiag_window hY hk hYkj
              exact absurd ⟨hwin.2, hjo'⟩ (howner j)),
        Finset.sum_eq_single_of_mem i (Finset.mem_univ i)
          (fun k _ hk => by
            by_cases hYik : Y i k = 0
            · rw [hYik, zero_mul]
            · obtain ⟨hio', _⟩ :=
                richBlock_offdiag_window hY (fun h => hk h.symm) hYik
              exact absurd ⟨hwin.1, hio'⟩ (howner i)),
        hYjj, hYii]
      exact mul_comm _ _
    · by_cases hwin' : richNodeOwner i = some o' ∧ richNodeOwner j = some o'
      · -- both endpoints in o''s window: X contributes its scalar
        have hXii : X i i = X 0 0 :=
          richBlock_offwindow_diag hX (fun h => howner i ⟨h, hwin'.1⟩)
        have hXjj : X j j = X 0 0 :=
          richBlock_offwindow_diag hX (fun h => howner j ⟨h, hwin'.2⟩)
        rw [Finset.sum_eq_single_of_mem i (Finset.mem_univ i)
            (fun k _ hk => by
              by_cases hXik : X i k = 0
              · rw [hXik, zero_mul]
              · obtain ⟨hio, _⟩ :=
                  richBlock_offdiag_window hX (fun h => hk h.symm) hXik
                exact absurd ⟨hio, hwin'.1⟩ (howner i)),
          Finset.sum_eq_single_of_mem j (Finset.mem_univ j)
            (fun k _ hk => by
              by_cases hXkj : X k j = 0
              · rw [hXkj, mul_zero]
              · obtain ⟨_, hjo⟩ := richBlock_offdiag_window hX hk hXkj
                exact absurd ⟨hjo, hwin'.2⟩ (howner j))]
        rw [hXii, hXjj]
        exact mul_comm _ _
      · -- mixed or outside: every product term vanishes on both sides
        have hXij : X i j = 0 := by
          by_contra hval
          exact hwin (richBlock_offdiag_window hX hij hval)
        have hYij : Y i j = 0 := by
          by_contra hval
          exact hwin' (richBlock_offdiag_window hY hij hval)
        have hL : ∑ k, X i k * Y k j = 0 := by
          refine Finset.sum_eq_zero fun k _ => ?_
          by_cases hik : i = k
          · subst hik
            rw [hYij, mul_zero]
          · by_cases hXik : X i k = 0
            · rw [hXik, zero_mul]
            · obtain ⟨hio, hko⟩ :=
                richBlock_offdiag_window hX hik hXik
              by_cases hkj : k = j
              · subst hkj
                exact absurd hXij hXik
              · by_cases hYkj : Y k j = 0
                · rw [hYkj, mul_zero]
                · obtain ⟨hko', _⟩ :=
                    richBlock_offdiag_window hY hkj hYkj
                  exact absurd ⟨hko, hko'⟩ (howner k)
        have hR : ∑ k, Y i k * X k j = 0 := by
          refine Finset.sum_eq_zero fun k _ => ?_
          by_cases hik : i = k
          · subst hik
            rw [hXij, mul_zero]
          · by_cases hYik : Y i k = 0
            · rw [hYik, zero_mul]
            · obtain ⟨hio', hko'⟩ :=
                richBlock_offdiag_window hY hik hYik
              by_cases hkj : k = j
              · subst hkj
                exact absurd hYij hYik
              · by_cases hXkj : X k j = 0
                · rw [hXkj, mul_zero]
                · obtain ⟨hko, _⟩ :=
                    richBlock_offdiag_window hX hkj hXkj
                  exact absurd ⟨hko, hko'⟩ (howner k)
        rw [hL, hR]

/-! ## The restriction system -/

/-- Pull the diagonal back along a basis self-map. -/
def richDiagPull (σ : Fin 81 → Fin 81) (X : Matrix (Fin 81) (Fin 81) ℂ) :
    Matrix (Fin 81) (Fin 81) ℂ :=
  Matrix.diagonal fun i => X (σ i) (σ i)

theorem richDiagPull_one (σ : Fin 81 → Fin 81) :
    richDiagPull σ (1 : Matrix (Fin 81) (Fin 81) ℂ) = 1 := by
  simp [richDiagPull]

theorem richDiagPull_add (σ : Fin 81 → Fin 81)
    (X Y : Matrix (Fin 81) (Fin 81) ℂ) :
    richDiagPull σ (X + Y) = richDiagPull σ X + richDiagPull σ Y := by
  simp [richDiagPull, Matrix.diagonal_add]

theorem richDiagPull_zero (σ : Fin 81 → Fin 81) :
    richDiagPull σ (0 : Matrix (Fin 81) (Fin 81) ℂ) = 0 := by
  simp [richDiagPull]

theorem richDiagPull_star (σ : Fin 81 → Fin 81)
    (X : Matrix (Fin 81) (Fin 81) ℂ) :
    richDiagPull σ (star X) = star (richDiagPull σ X) := by
  simp [richDiagPull, Matrix.star_eq_conjTranspose,
    Matrix.diagonal_conjTranspose]

theorem richDiagPull_algebraMap (σ : Fin 81 → Fin 81) (z : ℂ) :
    richDiagPull σ (algebraMap ℂ (Matrix (Fin 81) (Fin 81) ℂ) z) =
      algebraMap ℂ (Matrix (Fin 81) (Fin 81) ℂ) z := by
  rw [← richDiagonal_const_eq_algebraMap]
  simp [richDiagPull]

/-- The constant collapse of any matrix lands in the bottom scalars. -/
theorem richConstPull_mem_bot (X : Matrix (Fin 81) (Fin 81) ℂ) :
    richDiagPull (fun _ => 0) X ∈ richRegionAlgebra .bot := by
  refine ⟨fun i j hrel => ?_, fun i j _ _ => ?_⟩
  · have hij : i ≠ j := fun he => by
      subst he
      exact Bool.noConfusion ((richBotBlock.rel_refl i).symm.trans hrel)
    exact Matrix.diagonal_apply_ne _ hij
  · rw [richDiagPull, Matrix.diagonal_apply_eq, Matrix.diagonal_apply_eq]

/-- The character onto the bottom scalars from a block subalgebra whose
base point is a relation singleton. -/
noncomputable def richCharToBot (B : RichBlockRelation)
    (h0 : ∀ k, B.rel 0 k = true → k = 0) :
    B.subalgebra →⋆ₐ[ℂ] richRegionAlgebra .bot where
  toFun X := ⟨richDiagPull (fun _ => 0) X.1, richConstPull_mem_bot X.1⟩
  map_one' := Subtype.ext (richDiagPull_one _)
  map_mul' X Y := Subtype.ext (by
    show richDiagPull (fun _ => 0) (X.1 * Y.1) =
      richDiagPull (fun _ => 0) X.1 * richDiagPull (fun _ => 0) Y.1
    unfold richDiagPull
    rw [Matrix.diagonal_mul_diagonal]
    congr 1
    funext i
    rw [Matrix.mul_apply]
    refine Finset.sum_eq_single_of_mem 0 (Finset.mem_univ _)
      fun k _ hk => ?_
    have hrel : B.rel 0 k = false := by
      cases hbk : B.rel 0 k
      · rfl
      · exact absurd (h0 k hbk) hk
    rw [X.2.1 0 k hrel, zero_mul])
  map_zero' := Subtype.ext (richDiagPull_zero _)
  map_add' X Y := Subtype.ext (richDiagPull_add _ X.1 Y.1)
  map_star' X := Subtype.ext (richDiagPull_star _ X.1)
  commutes' z := Subtype.ext (richDiagPull_algebraMap _ z)

/-- The matrix value of the block restriction onto one observer window:
block entries kept on the window, the base-point scalar on the
off-window diagonal, zero elsewhere. -/
def richBlockPull (o : RichObserver) (X : Matrix (Fin 81) (Fin 81) ℂ) :
    Matrix (Fin 81) (Fin 81) ℂ :=
  fun i j => if richRelObs o i j = true then
    X (richCollapse o i) (richCollapse o j) else 0

theorem richBlockPull_mem (o : RichObserver)
    (X : Matrix (Fin 81) (Fin 81) ℂ) :
    richBlockPull o X ∈ richRegionAlgebra (RichRegion.regionOf o) := by
  rw [richRegionAlgebra_regionOf]
  refine ⟨fun i j h0 => ?_, fun i j hi hj => ?_⟩
  · have h0' : richRelObs o i j = false := h0
    simp only [richBlockPull]
    rw [if_neg (ne_true_of_eq_false' h0')]
  · have hoi : ¬ richNodeOwner i = some o := of_decide_eq_true hi
    have hoj : ¬ richNodeOwner j = some o := of_decide_eq_true hj
    have hri : richRelObs o i i = true := (richObsBlock o).rel_refl i
    have hrj : richRelObs o j j = true := (richObsBlock o).rel_refl j
    simp only [richBlockPull]
    rw [if_pos hri, if_pos hrj,
      richCollapse_of_not_owner hoi, richCollapse_of_not_owner hoj]

/-- The block restriction fixes every member of the window algebra. -/
theorem richBlockPull_fixes {o : RichObserver}
    {X : Matrix (Fin 81) (Fin 81) ℂ}
    (hX : X ∈ richRegionAlgebra (RichRegion.regionOf o)) :
    richBlockPull o X = X := by
  rw [richRegionAlgebra_regionOf] at hX
  ext i j
  simp only [richBlockPull]
  by_cases hrel : richRelObs o i j = true
  · rw [if_pos hrel]
    by_cases hij : i = j
    · subst hij
      by_cases ho : richNodeOwner i = some o
      · rw [richCollapse_of_owner ho]
      · rw [richCollapse_of_not_owner ho]
        exact (richBlock_offwindow_diag hX ho).symm
    · obtain ⟨hoi, hoj⟩ := richRelObs_offdiag_owner o i j hrel hij
      rw [richCollapse_of_owner hoi, richCollapse_of_owner hoj]
  · rw [if_neg hrel]
    exact (hX.1 i j (bool_eq_false' hrel)).symm

/-- Entry support of a top member: entries vanish off the top
relation. -/
theorem richTop_entry_zero {X : Matrix (Fin 81) (Fin 81) ℂ}
    (hX : X ∈ richRegionAlgebra .top) {i j : Fin 81}
    (h : richRelTop i j = false) : X i j = 0 :=
  hX.1 i j h

/-- On an owned index, a top member's entries vanish off the owner's
relation. -/
theorem richTop_entry_zero_of_owner {X : Matrix (Fin 81) (Fin 81) ℂ}
    (hX : X ∈ richRegionAlgebra .top) {o : RichObserver} {i j : Fin 81}
    (hio : richNodeOwner i = some o) (h : richRelObs o i j = false) :
    X i j = 0 :=
  hX.1 i j ((richRelTop_of_owner o i j hio).trans h)

/-- The restriction from the top block algebra onto one observer
window: a star-algebra homomorphism keeping the realized blocks. -/
noncomputable def richPullback (o : RichObserver) :
    richRegionAlgebra .top →⋆ₐ[ℂ]
      richRegionAlgebra (RichRegion.regionOf o) where
  toFun X := ⟨richBlockPull o X.1, richBlockPull_mem o X.1⟩
  map_one' := Subtype.ext (by
    show richBlockPull o 1 = 1
    ext i j
    simp only [richBlockPull]
    by_cases hrel : richRelObs o i j = true
    · rw [if_pos hrel]
      by_cases hij : i = j
      · subst hij
        rw [Matrix.one_apply_eq, Matrix.one_apply_eq]
      · obtain ⟨hoi, hoj⟩ := richRelObs_offdiag_owner o i j hrel hij
        rw [richCollapse_of_owner hoi, richCollapse_of_owner hoj]
    · have hij : i ≠ j := fun he => by
        subst he
        exact hrel ((richObsBlock o).rel_refl i)
      rw [if_neg hrel, Matrix.one_apply_ne hij])
  map_mul' X Y := Subtype.ext (by
    show richBlockPull o (X.1 * Y.1) =
      richBlockPull o X.1 * richBlockPull o Y.1
    ext i j
    by_cases hrel : richRelObs o i j = true
    · by_cases hoi : richNodeOwner i = some o
      · have hoj : richNodeOwner j = some o := by
          by_cases hij : i = j
          · subst hij; exact hoi
          · exact (richRelObs_offdiag_owner o i j hrel hij).2
        have hL : richBlockPull o (X.1 * Y.1) i j =
            ∑ k, X.1 i k * Y.1 k j := by
          simp only [richBlockPull]
          rw [if_pos hrel, richCollapse_of_owner hoi,
            richCollapse_of_owner hoj, Matrix.mul_apply]
        have hR : ∀ k, richBlockPull o X.1 i k * richBlockPull o Y.1 k j =
            X.1 i k * Y.1 k j := by
          intro k
          by_cases hik : richRelObs o i k = true
          · have hok : richNodeOwner k = some o :=
              richRelObs_owner_closed o i k hik hoi
            by_cases hkj : richRelObs o k j = true
            · simp only [richBlockPull]
              rw [if_pos hik, if_pos hkj, richCollapse_of_owner hoi,
                richCollapse_of_owner hok, richCollapse_of_owner hoj]
            · simp only [richBlockPull]
              rw [if_neg hkj, mul_zero,
                richTop_entry_zero_of_owner Y.2
                  (richRelObs_owner_closed o i k hik hoi)
                  (bool_eq_false' hkj), mul_zero]
          · simp only [richBlockPull]
            rw [if_neg hik, zero_mul,
              richTop_entry_zero_of_owner X.2 hoi (bool_eq_false' hik),
              zero_mul]
        rw [hL, Matrix.mul_apply]
        exact Finset.sum_congr rfl fun k _ => (hR k).symm
      · have hij : i = j := richRelObs_not_owner_singleton o i j hoi hrel
        subst hij
        have hσ : richCollapse o i = 0 := richCollapse_of_not_owner hoi
        have h0i : ∀ m, m ≠ 0 → X.1 0 m = 0 := fun m hm =>
          richTop_entry_zero X.2 (by
            cases hrel0 : richRelTop 0 m
            · rfl
            · exact absurd (richRelTop_zero_singleton m hrel0) hm)
        have hL : richBlockPull o (X.1 * Y.1) i i = X.1 0 0 * Y.1 0 0 := by
          simp only [richBlockPull]
          rw [if_pos hrel, hσ, Matrix.mul_apply]
          refine Finset.sum_eq_single_of_mem 0 (Finset.mem_univ _)
            fun m _ hm => ?_
          rw [h0i m hm, zero_mul]
        have hR : (richBlockPull o X.1 * richBlockPull o Y.1) i i =
            X.1 0 0 * Y.1 0 0 := by
          rw [Matrix.mul_apply]
          rw [Finset.sum_eq_single_of_mem i (Finset.mem_univ _)
            (fun k _ hk => by
              have hik : richRelObs o i k = false := by
                cases htrue : richRelObs o i k
                · rfl
                · exact absurd
                    (richRelObs_not_owner_singleton o i k hoi htrue).symm hk
              simp only [richBlockPull]
              rw [if_neg (ne_true_of_eq_false' hik), zero_mul])]
          simp only [richBlockPull]
          rw [if_pos hrel, if_pos hrel, hσ]
        exact hL.trans hR.symm
    · have hL : richBlockPull o (X.1 * Y.1) i j = 0 := by
        simp only [richBlockPull]
        rw [if_neg hrel]
      have hR : (richBlockPull o X.1 * richBlockPull o Y.1) i j = 0 := by
        rw [Matrix.mul_apply]
        refine Finset.sum_eq_zero fun k _ => ?_
        by_cases hik : richRelObs o i k = true
        · by_cases hkj : richRelObs o k j = true
          · exact absurd ((richObsBlock o).rel_trans i k j hik hkj) hrel
          · simp only [richBlockPull]
            rw [if_neg hkj, mul_zero]
        · simp only [richBlockPull]
          rw [if_neg hik, zero_mul]
      exact hL.trans hR.symm)
  map_zero' := Subtype.ext (by
    show richBlockPull o 0 = 0
    ext i j
    simp [richBlockPull])
  map_add' X Y := Subtype.ext (by
    show richBlockPull o (X.1 + Y.1) =
      richBlockPull o X.1 + richBlockPull o Y.1
    ext i j
    simp only [richBlockPull, Matrix.add_apply]
    by_cases hrel : richRelObs o i j = true
    · rw [if_pos hrel, if_pos hrel, if_pos hrel]
    · rw [if_neg hrel, if_neg hrel, if_neg hrel, add_zero])
  map_star' X := Subtype.ext (by
    show richBlockPull o (star X.1) = star (richBlockPull o X.1)
    ext i j
    have hsym : richRelObs o i j = richRelObs o j i :=
      (richObsBlock o).rel_symm i j
    simp only [richBlockPull, Matrix.star_apply]
    by_cases hrel : richRelObs o i j = true
    · rw [if_pos hrel, if_pos (hsym.symm.trans hrel)]
    · rw [if_neg hrel, if_neg fun h => hrel (hsym.trans h), star_zero])
  commutes' z := Subtype.ext (by
    show richBlockPull o (algebraMap ℂ (Matrix (Fin 81) (Fin 81) ℂ) z) =
      algebraMap ℂ (Matrix (Fin 81) (Fin 81) ℂ) z
    rw [← richDiagonal_const_eq_algebraMap]
    ext i j
    simp only [richBlockPull]
    by_cases hrel : richRelObs o i j = true
    · rw [if_pos hrel]
      by_cases hij : i = j
      · subst hij
        rw [Matrix.diagonal_apply_eq, Matrix.diagonal_apply_eq]
      · obtain ⟨hoi, hoj⟩ := richRelObs_offdiag_owner o i j hrel hij
        rw [richCollapse_of_owner hoi, richCollapse_of_owner hoj]
    · have hij : i ≠ j := fun he => by
        subst he
        exact hrel ((richObsBlock o).rel_refl i)
      rw [if_neg hrel, Matrix.diagonal_apply_ne _ hij])

/-- The supplied restriction system: identities on equal regions,
characters onto the bottom, and block restrictions from the top. -/
noncomputable def richRestrict : ∀ U V : RichRegion,
    RichRegion.le V U →
      (richRegionAlgebra U →⋆ₐ[ℂ] richRegionAlgebra V)
  | .bot, .bot, _ => StarAlgHom.id ℂ _
  | .r86, .r86, _ => StarAlgHom.id ℂ _
  | .r88, .r88, _ => StarAlgHom.id ℂ _
  | .r247, .r247, _ => StarAlgHom.id ℂ _
  | .r384, .r384, _ => StarAlgHom.id ℂ _
  | .top, .top, _ => StarAlgHom.id ℂ _
  | .r86, .bot, _ => richCharToBot _ (richRelObs_zero_singleton .o86)
  | .r88, .bot, _ => richCharToBot _ (richRelObs_zero_singleton .o88)
  | .r247, .bot, _ => richCharToBot _ (richRelObs_zero_singleton .o247)
  | .r384, .bot, _ => richCharToBot _ (richRelObs_zero_singleton .o384)
  | .top, .bot, _ => richCharToBot _ richRelTop_zero_singleton
  | .top, .r86, _ => richPullback .o86
  | .top, .r88, _ => richPullback .o88
  | .top, .r247, _ => richPullback .o247
  | .top, .r384, _ => richPullback .o384
  | .bot, .r86, h => nomatch h
  | .bot, .r88, h => nomatch h
  | .bot, .r247, h => nomatch h
  | .bot, .r384, h => nomatch h
  | .bot, .top, h => nomatch h
  | .r86, .r88, h => nomatch h
  | .r86, .r247, h => nomatch h
  | .r86, .r384, h => nomatch h
  | .r86, .top, h => nomatch h
  | .r88, .r86, h => nomatch h
  | .r88, .r247, h => nomatch h
  | .r88, .r384, h => nomatch h
  | .r88, .top, h => nomatch h
  | .r247, .r86, h => nomatch h
  | .r247, .r88, h => nomatch h
  | .r247, .r384, h => nomatch h
  | .r247, .top, h => nomatch h
  | .r384, .r86, h => nomatch h
  | .r384, .r88, h => nomatch h
  | .r384, .r247, h => nomatch h
  | .r384, .top, h => nomatch h


/-! ## The tower -/

open EventAlgebra in
/-- The selected tower state: the declared pure state at the base
point.  The payload pins no state; this field is declared data. -/
noncomputable def richFibreState : EventAlgebra.StateMatrix 81 :=
  ⟨Matrix.diagonal (Pi.single 0 1),
    ((diagonalPartition 81).isEvent 0).posSemidef, by
      rw [Matrix.trace_diagonal]
      simp [Pi.single_apply]⟩

/-- The rich-fibre consensus tower: `Record () = Fin 32` is the
payload's record-class type, `Observer () = RichObserver` carries the
four greedy-disjoint observers, and the public algebra of each observer
is the diagonal span of its record classes on its window. -/
noncomputable def richFibreTower : ConsensusTower Unit where
  indexNonempty := ⟨()⟩
  commonRefinement := fun _ _ => ⟨(), le_rfl, le_rfl⟩
  dim := fun _ => 81
  Observer := fun _ => RichObserver
  observerFinite := fun _ => inferInstance
  Record := fun _ => Fin 32
  recordFinite := fun _ => inferInstance
  recordOrder := fun _ _ =>
    OPH.TimeOrderLedger.ObserverRecordOrder.discrete (Fin 32)
  publicAlgebra := fun _ o => richDiagSpan (richPubMap o)
  public_mul_comm := fun _ o X Y => richDiagSpan_mul_comm (richPubMap o) X Y
  recordElement := fun _ o c =>
    ⟨Matrix.diagonal fun i => if richPubMap o i = some c then 1 else 0,
      fun x => if x = some c then 1 else 0, rfl⟩
  state := fun _ _ => richFibreState.1
  state_isState := fun _ _ => richFibreState.2
  generator := fun _ _ => 0
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
    simp [OPH.TimeOrderLedger.ObserverRecordOrder.discrete]
  public_natural := by
    intro r s hrs o X hX
    exact hX
  record_element_natural := by intros; rfl
  state_natural := by intros; rfl
  generator_natural := by intros; simp

/-! ## Locality on the region lattice -/

/-- Window members commute across distinct observers. -/
theorem richWindow_commute {o o' : RichObserver} (hne : o ≠ o')
    (X : richRegionAlgebra (RichRegion.regionOf o))
    (Y : richRegionAlgebra (RichRegion.regionOf o'))  :
    Commute (X : Matrix (Fin 81) (Fin 81) ℂ)
      (Y : Matrix (Fin 81) (Fin 81) ℂ) := by
  have hX : (X : Matrix (Fin 81) (Fin 81) ℂ) ∈
      (richObsBlock o).subalgebra := by
    rw [← richRegionAlgebra_regionOf]
    exact X.2
  have hY : (Y : Matrix (Fin 81) (Fin 81) ℂ) ∈
      (richObsBlock o').subalgebra := by
    rw [← richRegionAlgebra_regionOf]
    exact Y.2
  exact richBlock_commute_disjoint hne hX hY

/-- Members of two declared-disjoint regions commute elementwise. -/
theorem richRegional_locality {U V : RichRegion}
    (h : RichRegion.disj U V)
    (X : richRegionAlgebra U) (Y : richRegionAlgebra V) :
    Commute (X : Matrix (Fin 81) (Fin 81) ℂ)
      (Y : Matrix (Fin 81) (Fin 81) ℂ) := by
  obtain ⟨hne, hUb, hUt, hVb, hVt⟩ := h
  cases U
  case bot => exact absurd rfl hUb
  case top => exact absurd rfl hUt
  case r86 =>
    cases V
    case bot => exact absurd rfl hVb
    case top => exact absurd rfl hVt
    case r86 => exact absurd rfl hne
    case r88 => exact richWindow_commute (o := .o86) (o' := .o88) (by decide) X Y
    case r247 => exact richWindow_commute (o := .o86) (o' := .o247) (by decide) X Y
    case r384 => exact richWindow_commute (o := .o86) (o' := .o384) (by decide) X Y
  case r88 =>
    cases V
    case bot => exact absurd rfl hVb
    case top => exact absurd rfl hVt
    case r86 => exact richWindow_commute (o := .o88) (o' := .o86) (by decide) X Y
    case r88 => exact absurd rfl hne
    case r247 => exact richWindow_commute (o := .o88) (o' := .o247) (by decide) X Y
    case r384 => exact richWindow_commute (o := .o88) (o' := .o384) (by decide) X Y
  case r247 =>
    cases V
    case bot => exact absurd rfl hVb
    case top => exact absurd rfl hVt
    case r86 => exact richWindow_commute (o := .o247) (o' := .o86) (by decide) X Y
    case r88 => exact richWindow_commute (o := .o247) (o' := .o88) (by decide) X Y
    case r247 => exact absurd rfl hne
    case r384 => exact richWindow_commute (o := .o247) (o' := .o384) (by decide) X Y
  case r384 =>
    cases V
    case bot => exact absurd rfl hVb
    case top => exact absurd rfl hVt
    case r86 => exact richWindow_commute (o := .o384) (o' := .o86) (by decide) X Y
    case r88 => exact richWindow_commute (o := .o384) (o' := .o88) (by decide) X Y
    case r247 => exact richWindow_commute (o := .o384) (o' := .o247) (by decide) X Y
    case r384 => exact absurd rfl hne

/-! ## The rich-fibre finite causal observer net -/

set_option maxHeartbeats 2000000 in
/-- The rich-fibre net: the region lattice `RichRegion`, one genuinely
noncommutative block algebra at every observer window, the
character-and-block restriction system, distinct-window disjointness,
and identity repair. -/
noncomputable def richFibreNet : FiniteCausalObserverNet richFibreTower where
  Region := fun _ => RichRegion
  regionFintype := fun _ => inferInstance
  regionNonempty := fun _ => ⟨.bot⟩
  regionLE := fun _ => RichRegion.le
  regionLE_refl := fun _ U => RichRegion.le_refl U
  regionLE_trans := fun _ {U V W} h1 h2 => RichRegion.le_trans U V W h1 h2
  regionLE_antisymm := fun _ {U V} h1 h2 => RichRegion.le_antisymm U V h1 h2
  overlap := fun _ => RichRegion.meet
  overlap_le_left := fun _ U V => RichRegion.meet_le_left U V
  overlap_le_right := fun _ U V => RichRegion.meet_le_right U V
  le_overlap := fun _ {W U V} h1 h2 => RichRegion.le_meet W U V h1 h2
  disjoint := fun _ => RichRegion.disj
  disjoint_symm := fun _ {U V} h =>
    ⟨h.1.symm, h.2.2.2.1, h.2.2.2.2, h.2.1, h.2.2.1⟩
  disjoint_irrefl := fun _ U h => h.1 rfl
  localAlgebra := fun _ => richRegionAlgebra
  isotony := fun _ {U V} h => richRegionAlgebra_isotony h
  locality := fun _ {U V} h X Y => richRegional_locality h X Y
  restrict := fun _ {U V} h => richRestrict U V h
  restrict_refl := by
    intro r U X
    cases U <;> rfl
  restrict_trans := by
    intro r U V W hVU hWV X
    cases U <;> cases V <;> cases W <;>
      first
        | rfl
        | exact absurd hVU (by decide)
        | exact absurd hWV (by decide)
  restrict_inclusion := by
    intro r U V hUV X
    cases U <;> cases V <;>
      first
        | rfl
        | exact absurd hUV (by decide)
        | exact Subtype.ext (richBlockPull_fixes X.2)
        | exact Subtype.ext (richBot_mem_iff.mp X.2).symm
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

/-! ## Region separation receipts -/

/-- Every pair of distinct observer windows is declared disjoint; the
kernel receipt `richSupport_pairwise_disjoint` (preregistered gate G3)
backs the declaration in the source. -/
theorem richFibreNet_disjoint_supports :
    ∀ o o' : RichObserver, o ≠ o' →
      richFibreNet.disjoint () (RichRegion.regionOf o)
        (RichRegion.regionOf o') := by
  intro o o' h
  cases o <;> cases o' <;>
    first
      | exact absurd rfl h
      | exact ⟨by decide, by decide, by decide, by decide, by decide⟩

/-- Separation: regional observables of distinct windows commute. -/
theorem richFibreNet_separation (o o' : RichObserver) (h : o ≠ o')
    (X : richFibreNet.localAlgebra () (RichRegion.regionOf o))
    (Y : richFibreNet.localAlgebra () (RichRegion.regionOf o')) :
    Commute X.1 Y.1 :=
  richFibreNet.locality () (richFibreNet_disjoint_supports o o' h) X Y

/-! ## The realized matrix units and per-region noncommutativity -/

/-- The matrix unit at one ambient entry. -/
def richFibreUnit (a b : Fin 81) : Matrix (Fin 81) (Fin 81) ℂ :=
  fun k l => if k = a ∧ l = b then 1 else 0

/-- A matrix unit on a related window pair lies in the window's block
algebra. -/
theorem richFibreUnit_mem {o : RichObserver} {a b : Fin 81}
    (ha : richNodeOwner a = some o) (hb : richNodeOwner b = some o)
    (hrel : richRelObs o a b = true) :
    richFibreUnit a b ∈ richRegionAlgebra (RichRegion.regionOf o) := by
  rw [richRegionAlgebra_regionOf]
  refine ⟨fun i j h0 => ?_, fun i j hi hj => ?_⟩
  · by_cases h : i = a ∧ j = b
    · obtain ⟨rfl, rfl⟩ := h
      exact absurd hrel (ne_true_of_eq_false' h0)
    · simp [richFibreUnit, h]
  · have hia : ¬ i = a := fun he => by
      subst he
      exact of_decide_eq_true hi ha
    have hja : ¬ j = a := fun he => by
      subst he
      exact of_decide_eq_true hj ha
    simp [richFibreUnit, hia, hja]

/-- The product of opposite fibre units is the initial diagonal unit. -/
theorem richFibreUnit_mul (a b : Fin 81) (_hab : a ≠ b) :
    richFibreUnit a b * richFibreUnit b a = richFibreUnit a a := by
  ext k l
  rw [Matrix.mul_apply]
  rw [Finset.sum_eq_single_of_mem b (Finset.mem_univ b)
    (fun m _ hm => by simp [richFibreUnit, hm])]
  by_cases hk : k = a
  · by_cases hl : l = a
    · simp [richFibreUnit, hk, hl]
    · simp [richFibreUnit, hk, hl]
  · simp [richFibreUnit, hk]

/-- Distinct diagonal units differ. -/
theorem richFibreUnit_diag_ne (a b : Fin 81) (hab : a ≠ b) :
    richFibreUnit a a ≠ richFibreUnit b b := by
  intro h
  have := congrFun (congrFun h a) a
  simp [richFibreUnit, hab] at this

/-- **Every window region of the rich-fibre net is noncommutative**: the
realized designated split fibre of each observer supplies a
noncommuting pair of regional observables.  This is the receipt the
truncation-8 net could earn at one region only. -/
theorem richRegional_noncommutative_all (o : RichObserver) :
    ∃ X Y : richRegionAlgebra (RichRegion.regionOf o),
      X.1 * Y.1 ≠ Y.1 * X.1 := by
  cases o
  case o86 =>
    refine ⟨⟨richFibreUnit 3 7,
        richFibreUnit_mem (by decide) (by decide) (by decide)⟩,
      ⟨richFibreUnit 7 3,
        richFibreUnit_mem (by decide) (by decide) (by decide)⟩, ?_⟩
    rw [richFibreUnit_mul 3 7 (by decide),
      richFibreUnit_mul 7 3 (by decide)]
    exact richFibreUnit_diag_ne 3 7 (by decide)
  case o88 =>
    refine ⟨⟨richFibreUnit 24 35,
        richFibreUnit_mem (by decide) (by decide) (by decide)⟩,
      ⟨richFibreUnit 35 24,
        richFibreUnit_mem (by decide) (by decide) (by decide)⟩, ?_⟩
    rw [richFibreUnit_mul 24 35 (by decide),
      richFibreUnit_mul 35 24 (by decide)]
    exact richFibreUnit_diag_ne 24 35 (by decide)
  case o247 =>
    refine ⟨⟨richFibreUnit 42 57,
        richFibreUnit_mem (by decide) (by decide) (by decide)⟩,
      ⟨richFibreUnit 57 42,
        richFibreUnit_mem (by decide) (by decide) (by decide)⟩, ?_⟩
    rw [richFibreUnit_mul 42 57 (by decide),
      richFibreUnit_mul 57 42 (by decide)]
    exact richFibreUnit_diag_ne 42 57 (by decide)
  case o384 =>
    refine ⟨⟨richFibreUnit 65 66,
        richFibreUnit_mem (by decide) (by decide) (by decide)⟩,
      ⟨richFibreUnit 66 65,
        richFibreUnit_mem (by decide) (by decide) (by decide)⟩, ?_⟩
    rw [richFibreUnit_mul 65 66 (by decide),
      richFibreUnit_mul 66 65 (by decide)]
    exact richFibreUnit_diag_ne 65 66 (by decide)

/-- The regional-factor localization receipt: at every observer, the
four matrix units of its designated split fibre lie in its regional
algebra and multiply as a full two-by-two matrix factor. -/
theorem richDesignatedFactor (o : RichObserver) :
    ∃ a b : Fin 81, a ≠ b ∧
      richNodeOwner a = some o ∧ richNodeOwner b = some o ∧
      richFibreUnit a b ∈ richRegionAlgebra (RichRegion.regionOf o) ∧
      richFibreUnit b a ∈ richRegionAlgebra (RichRegion.regionOf o) ∧
      richFibreUnit a b * richFibreUnit b a = richFibreUnit a a ∧
      richFibreUnit b a * richFibreUnit a b = richFibreUnit b b := by
  cases o
  case o86 =>
    exact ⟨3, 7, by decide, by decide, by decide,
      richFibreUnit_mem (by decide) (by decide) (by decide),
      richFibreUnit_mem (by decide) (by decide) (by decide),
      richFibreUnit_mul 3 7 (by decide), richFibreUnit_mul 7 3 (by decide)⟩
  case o88 =>
    exact ⟨24, 35, by decide, by decide, by decide,
      richFibreUnit_mem (by decide) (by decide) (by decide),
      richFibreUnit_mem (by decide) (by decide) (by decide),
      richFibreUnit_mul 24 35 (by decide),
      richFibreUnit_mul 35 24 (by decide)⟩
  case o247 =>
    exact ⟨42, 57, by decide, by decide, by decide,
      richFibreUnit_mem (by decide) (by decide) (by decide),
      richFibreUnit_mem (by decide) (by decide) (by decide),
      richFibreUnit_mul 42 57 (by decide),
      richFibreUnit_mul 57 42 (by decide)⟩
  case o384 =>
    exact ⟨65, 66, by decide, by decide, by decide,
      richFibreUnit_mem (by decide) (by decide) (by decide),
      richFibreUnit_mem (by decide) (by decide) (by decide),
      richFibreUnit_mul 65 66 (by decide),
      richFibreUnit_mul 66 65 (by decide)⟩


/-! ## The coverage computation for the window family under the top -/

/-- The diagonal projector onto one observer window. -/
def richSupportProjector (o : RichObserver) :
    Matrix (Fin 81) (Fin 81) ℂ :=
  Matrix.diagonal fun i => if richNodeOwner i = some o then 1 else 0

/-- The window projector lies in every window's block algebra whose
observer owns it, and in particular in its own. -/
theorem richSupportProjector_mem (o : RichObserver) :
    richSupportProjector o ∈
      richRegionAlgebra (RichRegion.regionOf o) := by
  rw [richRegionAlgebra_regionOf]
  refine ⟨fun i j h0 => ?_, fun i j hi hj => ?_⟩
  · have hij : i ≠ j := fun he => by
      subst he
      exact Bool.noConfusion (((richObsBlock o).rel_refl i).symm.trans h0)
    exact Matrix.diagonal_apply_ne _ hij
  · simp only [richSupportProjector, Matrix.diagonal_apply_eq]
    rw [if_neg (of_decide_eq_true hi), if_neg (of_decide_eq_true hj)]

/-- The diagonal unit at the base point. -/
def richAmbientUnit : Matrix (Fin 81) (Fin 81) ℂ :=
  Matrix.diagonal fun i => if i = 0 then 1 else 0

/-- The base-point unit is the identity minus the four window
projectors: the windows partition the ambient basis away from the base
point. -/
theorem richAmbientUnit_eq :
    richAmbientUnit = 1 - (richSupportProjector .o86 +
      richSupportProjector .o88 + richSupportProjector .o247 +
        richSupportProjector .o384) := by
  unfold richAmbientUnit richSupportProjector
  rw [← Matrix.diagonal_one, Matrix.diagonal_add, Matrix.diagonal_add,
    Matrix.diagonal_add, Matrix.diagonal_sub]
  congr 1
  funext i
  rcases h : richNodeOwner i with _ | o
  · have h0 : i = 0 := (richNodeOwner_none_iff i).mp h
    subst h0
    simp
  · have h0 : i ≠ 0 := by
      intro he
      rw [(richNodeOwner_none_iff i).mpr he] at h
      simp at h
    cases o <;> simp [h0]

/-- The truncation of a matrix onto one observer window: entries on
top-related pairs led by the window, zero elsewhere. -/
def richSupportTruncation (o : RichObserver)
    (X : Matrix (Fin 81) (Fin 81) ℂ) : Matrix (Fin 81) (Fin 81) ℂ :=
  fun i j => if richNodeOwner i = some o ∧ richRelTop i j = true then
    X i j else 0

/-- Every window truncation lies in the window's block algebra. -/
theorem richSupportTruncation_mem (o : RichObserver)
    (X : Matrix (Fin 81) (Fin 81) ℂ) :
    richSupportTruncation o X ∈
      richRegionAlgebra (RichRegion.regionOf o) := by
  rw [richRegionAlgebra_regionOf]
  refine ⟨fun i j h0 => ?_, fun i j hi hj => ?_⟩
  · simp only [richSupportTruncation]
    rw [if_neg]
    rintro ⟨h1, h2⟩
    rw [richRelTop_of_owner o i j h1] at h2
    exact Bool.noConfusion (h0.symm.trans h2)
  · simp only [richSupportTruncation]
    rw [if_neg fun hc => (of_decide_eq_true hi) hc.1,
      if_neg fun hc => (of_decide_eq_true hj) hc.1]

/-- Exact decomposition of a top observable into the four window
truncations and the base-point scalar term. -/
theorem richTop_decomposition {X : Matrix (Fin 81) (Fin 81) ℂ}
    (hX : X ∈ richRegionAlgebra .top) :
    X = richSupportTruncation .o86 X + richSupportTruncation .o88 X +
      richSupportTruncation .o247 X + richSupportTruncation .o384 X +
        X 0 0 • richAmbientUnit := by
  ext i j
  simp only [Matrix.add_apply, Matrix.smul_apply, richSupportTruncation,
    richAmbientUnit]
  by_cases hij : i = j
  · subst hij
    rw [Matrix.diagonal_apply_eq]
    rcases h : richNodeOwner i with _ | o
    · have h0 : i = 0 := (richNodeOwner_none_iff i).mp h
      subst h0
      simp
    · have h0 : i ≠ 0 := by
        intro he
        rw [(richNodeOwner_none_iff i).mpr he] at h
        simp at h
      have hrefl : richRelTop i i = true := richTopBlock.rel_refl i
      cases o <;> simp [h0, hrefl]
  · rw [Matrix.diagonal_apply_ne _ hij]
    by_cases hrel : richRelTop i j = true
    · have hne : richNodeOwner i ≠ none :=
        richRelTop_offdiag_owner i j hrel (fun he => hij he)
      rcases h : richNodeOwner i with _ | o
      · exact absurd h hne
      · cases o <;> simp [hrel]
    · have h0 : X i j = 0 := hX.1 i j (bool_eq_false' hrel)
      simp [hrel, h0]

/-- The four-window family under the top region. -/
noncomputable def richWindowCover :
    richFibreNet.FiniteCover () RichRegion.top where
  regions := ({.r86, .r88, .r247, .r384} : Finset RichRegion)
  nonempty := by
    refine ⟨RichRegion.r86, ?_⟩
    change RichRegion.r86 ∈
      ({.r86, .r88, .r247, .r384} : Finset RichRegion)
    simp
  subregion := fun U _ => by cases U <;> trivial

/-- Every top observable is generated by the four window algebras. -/
theorem richTop_mem_familyJoin {X : Matrix (Fin 81) (Fin 81) ℂ}
    (hX : X ∈ richRegionAlgebra .top) :
    X ∈ richFibreNet.familyJoin richWindowCover := by
  have h86 : RichRegion.r86 ∈ richWindowCover.regions := by
    change RichRegion.r86 ∈
      ({.r86, .r88, .r247, .r384} : Finset RichRegion)
    simp
  have h88 : RichRegion.r88 ∈ richWindowCover.regions := by
    change RichRegion.r88 ∈
      ({.r86, .r88, .r247, .r384} : Finset RichRegion)
    simp
  have h247 : RichRegion.r247 ∈ richWindowCover.regions := by
    change RichRegion.r247 ∈
      ({.r86, .r88, .r247, .r384} : Finset RichRegion)
    simp
  have h384 : RichRegion.r384 ∈ richWindowCover.regions := by
    change RichRegion.r384 ∈
      ({.r86, .r88, .r247, .r384} : Finset RichRegion)
    simp
  have hT86 := richFibreNet.le_familyJoin richWindowCover h86
    (richSupportTruncation_mem .o86 X)
  have hT88 := richFibreNet.le_familyJoin richWindowCover h88
    (richSupportTruncation_mem .o88 X)
  have hT247 := richFibreNet.le_familyJoin richWindowCover h247
    (richSupportTruncation_mem .o247 X)
  have hT384 := richFibreNet.le_familyJoin richWindowCover h384
    (richSupportTruncation_mem .o384 X)
  have hP86 := richFibreNet.le_familyJoin richWindowCover h86
    (richSupportProjector_mem .o86)
  have hP88 := richFibreNet.le_familyJoin richWindowCover h88
    (richSupportProjector_mem .o88)
  have hP247 := richFibreNet.le_familyJoin richWindowCover h247
    (richSupportProjector_mem .o247)
  have hP384 := richFibreNet.le_familyJoin richWindowCover h384
    (richSupportProjector_mem .o384)
  have hAmb : richAmbientUnit ∈
      richFibreNet.familyJoin richWindowCover := by
    rw [richAmbientUnit_eq]
    exact sub_mem (one_mem _)
      (add_mem (add_mem (add_mem hP86 hP88) hP247) hP384)
  rw [richTop_decomposition hX]
  refine add_mem (add_mem (add_mem (add_mem hT86 hT88) hT247) hT384) ?_
  have hsmul : X 0 0 • richAmbientUnit =
      algebraMap ℂ (Matrix (Fin 81) (Fin 81) ℂ) (X 0 0) *
        richAmbientUnit :=
    Algebra.smul_def _ _
  rw [hsmul]
  exact mul_mem
    ((richFibreNet.familyJoin richWindowCover).algebraMap_mem (X 0 0)) hAmb

/-- **The coverage law, computed**: the four-window family generates the
top algebra. -/
theorem richWindowCover_coverageLaw :
    richFibreNet.CoverageLaw richWindowCover :=
  ⟨fun _X hX => richTop_mem_familyJoin hX⟩

/-- The family join of the four windows equals the top algebra. -/
theorem richWindowCover_familyJoin :
    richFibreNet.familyJoin richWindowCover =
      richFibreNet.localAlgebra () RichRegion.top :=
  richWindowCover_coverageLaw.join_eq

/-! ## Every window is load-bearing for coverage -/

/-- The top relation away from one observer: the representative map
fixes that observer's window pointwise. -/
def richRepDrop (o : RichObserver) : Fin 81 → Fin 81 := fun i =>
  if richNodeOwner i = some o then i else richRepTop i

/-- The dropped-window relation. -/
def richRelDrop (o : RichObserver) : Fin 81 → Fin 81 → Bool :=
  fun i j => decide (richRepDrop o i = richRepDrop o j)

/-- The dropped-window block relation: no ambient constraint. -/
def richDropBlock (o : RichObserver) : RichBlockRelation :=
  .ofRep (richRelDrop o) (fun _ => false) (richRepDrop o) (fun _ _ => rfl)
    (fun _ _ h _ => Bool.noConfusion h)

theorem richRepDrop_of_not_owner {o : RichObserver} {i : Fin 81}
    (h : ¬ richNodeOwner i = some o) : richRepDrop o i = richRepTop i := by
  simp [richRepDrop, h]

theorem richRepTop_of_owner {o : RichObserver} {i : Fin 81}
    (h : richNodeOwner i = some o) : richRepTop i = richRep o i := by
  simp [richRepTop, h]

/-- Away from its own window, every other observer's relation refines
the dropped-window relation. -/
theorem richRelObs_le_drop {o o' : RichObserver} (hne : o' ≠ o)
    (i j : Fin 81) (h : richRelObs o' i j = true) :
    richRelDrop o i j = true := by
  by_cases hij : i = j
  · subst hij
    exact (richDropBlock o).rel_refl i
  · obtain ⟨hio, hjo⟩ := richRelObs_offdiag_owner o' i j h hij
    have hio' : ¬ richNodeOwner i = some o := by
      rw [hio]
      exact fun hc => hne (Option.some.inj hc)
    have hjo' : ¬ richNodeOwner j = some o := by
      rw [hjo]
      exact fun hc => hne (Option.some.inj hc)
    have hrep : richRepDrop o i = richRepDrop o j := by
      rw [richRepDrop_of_not_owner hio', richRepDrop_of_not_owner hjo',
        richRepTop_of_owner hio, richRepTop_of_owner hjo]
      exact of_decide_eq_true h
    exact decide_eq_true hrep

/-- Every remaining region algebra embeds into the dropped-window block
algebra. -/
theorem richRegionAlgebra_le_dropBlock (o : RichObserver)
    (U : RichRegion) (hU : U ≠ RichRegion.regionOf o)
    (hUt : U ≠ RichRegion.top) :
    richRegionAlgebra U ≤ (richDropBlock o).subalgebra := by
  have key : ∀ (o' : RichObserver), o' ≠ o →
      richRegionAlgebra (RichRegion.regionOf o') ≤
        (richDropBlock o).subalgebra := by
    intro o' hne X hX
    rw [richRegionAlgebra_regionOf] at hX
    obtain ⟨h0, _⟩ := hX
    refine ⟨fun i j hd => ?_, fun i j hi _ => Bool.noConfusion hi⟩
    refine h0 i j ?_
    show richRelObs o' i j = false
    cases hV : richRelObs o' i j
    · rfl
    · exact absurd (richRelObs_le_drop hne i j hV)
        (ne_true_of_eq_false' hd)
  cases U
  case top => exact absurd rfl hUt
  case bot =>
    intro X hX
    have hX' := richBot_mem_iff.mp hX
    refine ⟨fun i j hrel => ?_, fun i j hi _ => Bool.noConfusion hi⟩
    have hij : i ≠ j := fun he => by
      subst he
      exact Bool.noConfusion (((richDropBlock o).rel_refl i).symm.trans hrel)
    conv_lhs => rw [hX']
    rw [Matrix.diagonal_apply_ne _ hij]
  case r86 =>
    exact key .o86 (fun he => hU (he ▸ rfl))
  case r88 =>
    exact key .o88 (fun he => hU (he ▸ rfl))
  case r247 =>
    exact key .o247 (fun he => hU (he ▸ rfl))
  case r384 =>
    exact key .o384 (fun he => hU (he ▸ rfl))

/-- The designated fibre unit of each observer is outside the
dropped-window block algebra: its window pair is unrelated once the
window is dropped. -/
theorem richDesignated_not_in_drop (o : RichObserver) :
    ∃ a b : Fin 81,
      richFibreUnit a b ∈ richRegionAlgebra (RichRegion.regionOf o) ∧
      richFibreUnit a b ∈ richRegionAlgebra .top ∧
      richFibreUnit a b ∉ (richDropBlock o).subalgebra := by
  have hmem : ∀ (a b : Fin 81), richRelDrop o a b = false →
      richFibreUnit a b ∉ (richDropBlock o).subalgebra := by
    intro a b hd hin
    have h := hin.1 a b hd
    simp [richFibreUnit] at h
  have htop : ∀ (a b : Fin 81), richNodeOwner a = some o →
      richNodeOwner b = some o → richRelObs o a b = true →
      richFibreUnit a b ∈ richRegionAlgebra .top := by
    intro a b ha hb hrel
    refine ⟨fun i j h0 => ?_, fun i j hi _ => Bool.noConfusion hi⟩
    by_cases h : i = a ∧ j = b
    · obtain ⟨rfl, rfl⟩ := h
      have h0' : richRelTop i j = false := h0
      rw [richRelTop_of_owner o i j ha] at h0'
      exact absurd hrel (ne_true_of_eq_false' h0')
    · simp [richFibreUnit, h]
  cases o
  case o86 =>
    exact ⟨3, 7,
      richFibreUnit_mem (by decide) (by decide) (by decide),
      htop 3 7 (by decide) (by decide) (by decide),
      hmem 3 7 (by decide)⟩
  case o88 =>
    exact ⟨24, 35,
      richFibreUnit_mem (by decide) (by decide) (by decide),
      htop 24 35 (by decide) (by decide) (by decide),
      hmem 24 35 (by decide)⟩
  case o247 =>
    exact ⟨42, 57,
      richFibreUnit_mem (by decide) (by decide) (by decide),
      htop 42 57 (by decide) (by decide) (by decide),
      hmem 42 57 (by decide)⟩
  case o384 =>
    exact ⟨65, 66,
      richFibreUnit_mem (by decide) (by decide) (by decide),
      htop 65 66 (by decide) (by decide) (by decide),
      hmem 65 66 (by decide)⟩

/-- **Negative coverage control at every window**: any family avoiding
one window and the top region itself fails the coverage law, because the
dropped window carries its own realized noncommutative block.  The
top-exclusion hypothesis is necessary: the singleton family at the top
trivially covers. -/
theorem richDropCover_not_coverageLaw (o : RichObserver)
    (C : richFibreNet.FiniteCover () RichRegion.top)
    (hC : ∀ U ∈ C.regions,
      U ≠ RichRegion.regionOf o ∧ U ≠ RichRegion.top) :
    ¬ richFibreNet.CoverageLaw C := by
  intro hcov
  obtain ⟨a, b, _hreg, htop, hdrop⟩ := richDesignated_not_in_drop o
  have hjoin : richFibreNet.familyJoin C ≤ (richDropBlock o).subalgebra :=
    richFibreNet.familyJoin_le fun U hU =>
      richRegionAlgebra_le_dropBlock o U (hC U hU).1 (hC U hU).2
  exact hdrop (hjoin (hcov.covered_le_join htop))

/-! ## Unique restriction gluing on the four-window cover -/

private theorem regionOf_mem_cover (o : RichObserver) :
    RichRegion.regionOf o ∈ richWindowCover.regions := by
  cases o <;>
    · change _ ∈ ({.r86, .r88, .r247, .r384} : Finset RichRegion)
      decide

/-- The window section of a compatible family, as an ambient matrix. -/
private noncomputable def secVal
    (F : richWindowCover.CompatibleFamily) (o : RichObserver) :
    Matrix (Fin 81) (Fin 81) ℂ :=
  (F.localSection (RichRegion.regionOf o) (regionOf_mem_cover o)).val

private theorem secVal_mem (F : richWindowCover.CompatibleFamily)
    (o : RichObserver) : secVal F o ∈ (richObsBlock o).subalgebra := by
  cases o
  case o86 =>
    exact (F.localSection (RichRegion.regionOf .o86)
      (regionOf_mem_cover .o86)).property
  case o88 =>
    exact (F.localSection (RichRegion.regionOf .o88)
      (regionOf_mem_cover .o88)).property
  case o247 =>
    exact (F.localSection (RichRegion.regionOf .o247)
      (regionOf_mem_cover .o247)).property
  case o384 =>
    exact (F.localSection (RichRegion.regionOf .o384)
      (regionOf_mem_cover .o384)).property

/-- Sections of one compatible family share their base-point scalar:
every window meet is the bottom region and both bottom characters
evaluate there. -/
private theorem family_scalar_eq
    (F : richWindowCover.CompatibleFamily) (o o' : RichObserver) :
    secVal F o 0 0 = secVal F o' 0 0 := by
  by_cases hoo : o = o'
  · subst hoo
    rfl
  · have hcomp := F.overlap_compatible
      (RichRegion.regionOf o) (regionOf_mem_cover o)
      (RichRegion.regionOf o') (regionOf_mem_cover o')
    have hval := congrArg Subtype.val hcomp
    revert hval
    cases o <;> cases o' <;>
      first
        | exact absurd rfl hoo
        | (intro hval
           exact congrFun (congrFun hval (0 : Fin 81)) (0 : Fin 81))

/-- The glued matrix of a compatible family on the four-window cover. -/
private noncomputable def familyGlue
    (F : richWindowCover.CompatibleFamily) :
    Matrix (Fin 81) (Fin 81) ℂ := fun i j =>
  match richNodeOwner i with
  | none => if i = j then secVal F .o86 0 0 else 0
  | some o =>
      if richNodeOwner j = some o ∧ richRelObs o i j = true then
        secVal F o i j
      else 0

private theorem familyGlue_entry_owner
    (F : richWindowCover.CompatibleFamily) {o : RichObserver}
    {i j : Fin 81} (hio : richNodeOwner i = some o) :
    familyGlue F i j =
      if richNodeOwner j = some o ∧ richRelObs o i j = true then
        secVal F o i j
      else 0 := by
  unfold familyGlue
  rw [hio]

private theorem familyGlue_entry_base
    (F : richWindowCover.CompatibleFamily) {i j : Fin 81}
    (hi : richNodeOwner i = none) :
    familyGlue F i j = if i = j then secVal F .o86 0 0 else 0 := by
  unfold familyGlue
  rw [hi]

/-- **Unique restriction gluing on the four-window cover**: every
compatible family has exactly one top observable restricting to it.
Together with the coverage law this is the reconstruction receipt the
committed descent interface consumes. -/
theorem richWindowCover_hasUniqueDescent :
    richWindowCover.HasUniqueDescent := by
  intro F
  have howner0 : richNodeOwner (0 : Fin 81) = none :=
    (richNodeOwner_none_iff 0).mpr rfl
  -- membership of the glue in the top algebra
  have hGmem : familyGlue F ∈ richRegionAlgebra .top := by
    refine ⟨fun i j h0 => ?_, fun i j hi _ => Bool.noConfusion hi⟩
    rcases howner : richNodeOwner i with _ | o
    · have hi0 : i = 0 := (richNodeOwner_none_iff i).mp howner
      subst hi0
      rw [familyGlue_entry_base F howner]
      by_cases hij : (0 : Fin 81) = j
      · subst hij
        exact absurd (richTopBlock.rel_refl 0) (ne_true_of_eq_false' h0)
      · rw [if_neg hij]
    · rw [familyGlue_entry_owner F howner, if_neg]
      rintro ⟨_hjo, hrel⟩
      have h0' : richRelObs o i j = false := by
        have := richRelTop_of_owner o i j howner
        exact this ▸ h0
      exact Bool.noConfusion (h0'.symm.trans hrel)
  -- the glue restricts to every section
  have key : ∀ o : RichObserver,
      richBlockPull o (familyGlue F) = secVal F o := by
    intro o
    ext i j
    simp only [richBlockPull]
    by_cases hrel : richRelObs o i j = true
    · rw [if_pos hrel]
      by_cases hio : richNodeOwner i = some o
      · have hjo := richRelObs_owner_closed o i j hrel hio
        rw [richCollapse_of_owner hio, richCollapse_of_owner hjo,
          familyGlue_entry_owner F hio, if_pos ⟨hjo, hrel⟩]
      · have hij := richRelObs_not_owner_singleton o i j hio hrel
        subst hij
        rw [richCollapse_of_not_owner hio]
        have hb : familyGlue F (0 : Fin 81) (0 : Fin 81) =
            secVal F .o86 0 0 := by
          rw [familyGlue_entry_base F howner0]
          simp
        rw [hb, family_scalar_eq F .o86 o]
        exact (richBlock_offwindow_diag (secVal_mem F o) hio).symm
    · rw [if_neg hrel]
      exact ((secVal_mem F o).1 i j (bool_eq_false' hrel)).symm
  have hGlob : richWindowCover.Globalizes F ⟨familyGlue F, hGmem⟩ := by
    intro U hU
    have hU' : U = .r86 ∨ U = .r88 ∨ U = .r247 ∨ U = .r384 := by
      have hU2 : U ∈ ({.r86, .r88, .r247, .r384} : Finset RichRegion) := hU
      cases U
      case bot => exact absurd hU2 (by decide)
      case top => exact absurd hU2 (by decide)
      case r86 => exact Or.inl rfl
      case r88 => exact Or.inr (Or.inl rfl)
      case r247 => exact Or.inr (Or.inr (Or.inl rfl))
      case r384 => exact Or.inr (Or.inr (Or.inr rfl))
    obtain h | h | h | h := hU' <;> subst h
    · exact Subtype.ext (key .o86)
    · exact Subtype.ext (key .o88)
    · exact Subtype.ext (key .o247)
    · exact Subtype.ext (key .o384)
  refine ⟨⟨familyGlue F, hGmem⟩, hGlob, ?_⟩
  -- uniqueness: a globalizing top observable is the glue entrywise
  rintro ⟨Y, hYmem⟩ hY
  apply Subtype.ext
  show (Y : Matrix (Fin 81) (Fin 81) ℂ) = familyGlue F
  have hYrestrict : ∀ o : RichObserver,
      richBlockPull o Y = secVal F o := fun o => by
    have h := hY (RichRegion.regionOf o) (regionOf_mem_cover o)
    have hval := congrArg Subtype.val h
    cases o <;> exact hval
  refine Matrix.ext fun (i : Fin 81) (j : Fin 81) => ?_
  rcases howner : richNodeOwner i with _ | o
  · have hi0 : i = 0 := (richNodeOwner_none_iff i).mp howner
    subst hi0
    by_cases hij : (0 : Fin 81) = j
    · subst hij
      have hb : familyGlue F (0 : Fin 81) (0 : Fin 81) =
          secVal F .o86 0 0 := by
        rw [familyGlue_entry_base F howner0]
        simp
      have hcalc : richBlockPull .o86 Y (0 : Fin 81) (0 : Fin 81) =
          Y (0 : Fin 81) (0 : Fin 81) := by
        have hri : richRelObs .o86 0 0 = true :=
          (richObsBlock .o86).rel_refl 0
        simp only [richBlockPull]
        rw [if_pos hri, richCollapse_of_not_owner (by
          rw [howner0]
          simp)]
      have h := congrFun (congrFun (hYrestrict .o86) (0 : Fin 81))
        (0 : Fin 81)
      rw [hb, ← h, hcalc]
    · have hb : familyGlue F (0 : Fin 81) j = 0 := by
        rw [familyGlue_entry_base F howner0, if_neg hij]
      rw [hb]
      refine hYmem.1 0 j ?_
      show richRelTop 0 j = false
      cases hrel : richRelTop 0 j
      · rfl
      · exact absurd (richRelTop_zero_singleton j hrel).symm hij
  · by_cases hcond : richNodeOwner j = some o ∧ richRelObs o i j = true
    · obtain ⟨hjo, hrel⟩ := hcond
      have hb : familyGlue F i j = secVal F o i j := by
        rw [familyGlue_entry_owner F howner, if_pos ⟨hjo, hrel⟩]
      have hcalc : richBlockPull o Y i j = Y i j := by
        simp only [richBlockPull]
        rw [if_pos hrel, richCollapse_of_owner howner,
          richCollapse_of_owner hjo]
      have h := congrFun (congrFun (hYrestrict o) i) j
      rw [hb, ← h, hcalc]
    · have hb : familyGlue F i j = 0 := by
        rw [familyGlue_entry_owner F howner, if_neg hcond]
      rw [hb]
      refine hYmem.1 i j ?_
      show richRelTop i j = false
      cases hrel : richRelTop i j
      · rfl
      · exfalso
        have hrel' : richRelObs o i j = true := by
          have := richRelTop_of_owner o i j howner
          exact this ▸ hrel
        exact hcond ⟨richRelObs_owner_closed o i j hrel' howner, hrel'⟩

/-- The descent packet of the rich-fibre net: every declared cover with
a proved unique-descent receipt is admitted, and the singleton family
supplies the universal existence branch. -/
noncomputable def richFibreDescent : richFibreNet.ObserverNetDescent where
  declaredCover := fun C => C.HasUniqueDescent
  cover_exists := fun r W =>
    ⟨richFibreNet.singletonCover r W,
      richFibreNet.singletonCover_hasUniqueDescent r W⟩
  unique_descent := fun _C hC => hC

/-- **Reconstruction on the four-window cover**: every compatible family
has exactly one globalization, and it lies in the family-generated
join. -/
theorem richWindowCover_reconstruction
    (F : richWindowCover.CompatibleFamily) :
    ∃! X : richFibreNet.localAlgebra () RichRegion.top,
      richWindowCover.Globalizes F X ∧
        X.val ∈ richFibreNet.familyJoin richWindowCover :=
  richWindowCover_coverageLaw.reconstruction richFibreDescent
    richWindowCover_hasUniqueDescent F

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.richBlock_commute_disjoint
#print axioms OPH.QFT.richRegional_locality
#print axioms OPH.QFT.richFibreNet_disjoint_supports
#print axioms OPH.QFT.richFibreNet_separation
#print axioms OPH.QFT.richRegional_noncommutative_all
#print axioms OPH.QFT.richDesignatedFactor
#print axioms OPH.QFT.richWindowCover_coverageLaw
#print axioms OPH.QFT.richWindowCover_familyJoin
#print axioms OPH.QFT.richDropCover_not_coverageLaw
#print axioms OPH.QFT.richWindowCover_hasUniqueDescent
#print axioms OPH.QFT.richWindowCover_reconstruction
