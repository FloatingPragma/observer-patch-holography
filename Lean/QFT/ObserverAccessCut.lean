import QFT.ObserverNetDescent

/-!
# Observer access cuts on finite causal observer nets

This module implements the access-cut layer of completion-plan issue `#712`
over the committed E1 interface.  An `ObserverAccessCut` assigns to every A3
observer one declared support region and one accessible star subalgebra of
the tower's private matrix algebra.  The accessible algebra contains the
observer's commutative public record algebra, localizes inside the regional
algebra of the declared observer region, and is stable under the declared
contravariant restrictions below that region.  `RepairStableAccessCut` is
the optional strengthening receipt for stability under the local repair
idempotents.

The observer-local algebra at a region is the meet of the regional local
algebra with the accessible algebra.  The inheritance theorems transport
the E1 receipts to this meet: isotony, elementwise commutation on declared
disjoint regions, restriction stability, selected-state expectation
compatibility with the ambient regional expectation, regional-repair fixed
points, and remote nondisturbance.

`AccessibleGlueClosure` is the deliverable-8 receipt: on every declared
descent family under the observer region, any ambient globalization of a
compatible family of accessible sections is accessible.  The theorem
`observerLocal_hasUniqueDescent` consumes exactly this receipt together
with an ambient `ObserverNetDescent` packet and returns unique descent
into the observer-local algebra.  The receipt is a separate premise: no
theorem in this module derives accessibility from descent, and whether
descent could force accessibility stays open until a countermodel is on
record.  The
universal singleton-family packet satisfies the receipt for every cut,
which keeps the receipt type inhabited without a nontrivial gluing claim.
A countermodel in which a declared family glues to an inaccessible element
requires a net whose restrictions fail to capture the glue algebraically
and is out of scope here.

The witness section instantiates the cut on a dimension-two model.  The
one-block partition drives the constant tower, so the public algebra is
the scalar span of the identity; the full diagonal partition supplies
every regional algebra; the accessible algebra is the diagonal algebra.
Both interposition facts are strict: the public scalars sit strictly
inside the accessible diagonal algebra, which sits strictly inside the
private matrix algebra.  In that model the accessible algebra equals
every regional algebra, so every observer-local meet collapses onto its
regional component.  The proper-meet section supplies a dimension-four
cut in the block-algebra shape of `QFT/NoncommutativeWitness.lean`: the
accessible algebra is the full diagonal algebra while the block region
carries the noncommutative algebra of matrices `diag(a, a, A)`, so the
observer-local meet at the block region is the computed intersection of
its two components, with separating elements on both sides.

Claim boundary: every construction here is a finite conditional structure
over the committed E1 interface.  The ambient private algebra is a
relational provenance object; no theorem grants any observer access to it,
and no physical region, coverage, instrument, spacetime, causal cone, or
continuum object is claimed.  E1 retains the noncommutative source
witness; source realization belongs to E2.
-/

namespace OPH.QFT

open OPH.Tower

universe u

variable {ι : Type u} [Preorder ι]

/-- A per-observer access cut on one finite causal observer net.

`observerRegion` is declared finite-region data and `accessibleAlgebra` is
declared algebraic data.  The fields squeeze the accessible algebra between
the observer's public record algebra and the local algebra of the declared
region, and make the declared restrictions below the observer region
preserve accessibility.  No field asserts that the cut is maximal,
canonical, or induced by an instrument. -/
structure ObserverAccessCut (T : ConsensusTower ι)
    (N : FiniteCausalObserverNet T) (r : ι) where
  /-- The declared support region of each observer. -/
  observerRegion : T.Observer r → N.Region r
  /-- The declared accessible star subalgebra of each observer. -/
  accessibleAlgebra : T.Observer r →
    StarSubalgebra ℂ (ConsensusTower.PrivateAlgebra T r)
  /-- Accessibility localizes in the declared observer region. -/
  accessible_le_region : ∀ o,
    accessibleAlgebra o ≤ N.localAlgebra r (observerRegion o)
  /-- The observer's public record algebra is accessible. -/
  public_le_accessible : ∀ o,
    T.publicAlgebra r o ≤ accessibleAlgebra o
  /-- Declared restrictions below the observer region preserve
  accessibility, stated at the element level. -/
  restrict_preserves : ∀ (o : T.Observer r) {U V : N.Region r},
    N.regionLE r U (observerRegion o) →
    ∀ (hVU : N.regionLE r V U) (X : N.localAlgebra r U),
      (X : ConsensusTower.PrivateAlgebra T r) ∈ accessibleAlgebra o →
        (N.restrict r hVU X : ConsensusTower.PrivateAlgebra T r) ∈
          accessibleAlgebra o

/-- Optional strengthening receipt: the local repair idempotents below the
observer region preserve accessibility.  This is a separate mixin so that
access cuts without repair stability remain expressible. -/
structure RepairStableAccessCut {T : ConsensusTower ι}
    {N : FiniteCausalObserverNet T} {r : ι}
    (A : ObserverAccessCut T N r) : Prop where
  repair_preserves : ∀ (o : T.Observer r) {U : N.Region r},
    N.regionLE r U (A.observerRegion o) →
    ∀ X : ConsensusTower.PrivateAlgebra T r,
      X ∈ A.accessibleAlgebra o →
        N.repair r U X ∈ A.accessibleAlgebra o

namespace ObserverAccessCut

variable {T : ConsensusTower ι} {N : FiniteCausalObserverNet T} {r : ι}
variable (A : ObserverAccessCut T N r)

/-! ## The observer-local algebra net -/

/-- The observer-local algebra at a region: the meet of the regional local
algebra with the observer's accessible algebra.  The meet is defined for
every region; below the declared observer region it is the observer-local
net of issue `#712`. -/
noncomputable def observerLocalAlgebra (o : T.Observer r) (U : N.Region r) :
    StarSubalgebra ℂ (ConsensusTower.PrivateAlgebra T r) :=
  N.localAlgebra r U ⊓ A.accessibleAlgebra o

@[simp]
theorem mem_observerLocalAlgebra {o : T.Observer r} {U : N.Region r}
    {X : ConsensusTower.PrivateAlgebra T r} :
    X ∈ A.observerLocalAlgebra o U ↔
      X ∈ N.localAlgebra r U ∧ X ∈ A.accessibleAlgebra o :=
  Iff.rfl

theorem observerLocal_le_local (o : T.Observer r) (U : N.Region r) :
    A.observerLocalAlgebra o U ≤ N.localAlgebra r U :=
  inf_le_left

theorem observerLocal_le_accessible (o : T.Observer r) (U : N.Region r) :
    A.observerLocalAlgebra o U ≤ A.accessibleAlgebra o :=
  inf_le_right

/-- At the declared observer region the meet is the whole accessible
algebra: localization makes the second component redundant there. -/
theorem observerLocalAlgebra_observerRegion (o : T.Observer r) :
    A.observerLocalAlgebra o (A.observerRegion o) = A.accessibleAlgebra o :=
  inf_eq_right.mpr (A.accessible_le_region o)

/-! ## Inheritance of the E1 receipts -/

/-- Inheritance, isotony: region inclusion gives inclusion of
observer-local algebras.  The statement holds for every region pair and in
particular below the declared observer region. -/
theorem observerLocal_isotony (o : T.Observer r) {U V : N.Region r}
    (hUV : N.regionLE r U V) :
    A.observerLocalAlgebra o U ≤ A.observerLocalAlgebra o V :=
  inf_le_inf (N.isotony r hUV) le_rfl

/-- Inheritance, locality: observer-local elements of declared-disjoint
regions commute elementwise, directly from the ambient E1 locality. -/
theorem observerLocal_commute (o : T.Observer r) {U V : N.Region r}
    (hUV : N.disjoint r U V) (X : A.observerLocalAlgebra o U)
    (Y : A.observerLocalAlgebra o V) :
    Commute (X : ConsensusTower.PrivateAlgebra T r) Y :=
  N.locality r hUV ⟨X.1, (A.mem_observerLocalAlgebra.mp X.2).1⟩
    ⟨Y.1, (A.mem_observerLocalAlgebra.mp Y.2).1⟩

/-- Inheritance, restriction stability: the declared restriction of an
observer-local element below the observer region is observer-local. -/
theorem restrict_mem_observerLocal (o : T.Observer r) {U V : N.Region r}
    (hU : N.regionLE r U (A.observerRegion o)) (hVU : N.regionLE r V U)
    (X : N.localAlgebra r U)
    (hX : (X : ConsensusTower.PrivateAlgebra T r) ∈ A.accessibleAlgebra o) :
    (N.restrict r hVU X : ConsensusTower.PrivateAlgebra T r) ∈
      A.observerLocalAlgebra o V :=
  A.mem_observerLocalAlgebra.mpr
    ⟨(N.restrict r hVU X).2, A.restrict_preserves o hU hVU X hX⟩

/-- The observer expectation of an observer-local element: the selected A3
density state paired with the underlying private observable. -/
noncomputable def observerExpectation (o : T.Observer r) (U : N.Region r)
    (X : A.observerLocalAlgebra o U) : ℂ :=
  (T.state r o * (X : ConsensusTower.PrivateAlgebra T r)).trace

/-- Inheritance, expectation compatibility: under the meet inclusion into
the regional algebra the observer expectation is the ambient regional
expectation.  Both sides are the same finite trace pairing, so the
identity holds definitionally. -/
@[simp]
theorem regionalExpectation_observerLocal (o : T.Observer r)
    (U : N.Region r) (X : A.observerLocalAlgebra o U) :
    N.regionalExpectation r o U
        (StarSubalgebra.inclusion (A.observerLocal_le_local o U) X) =
      A.observerExpectation o U X :=
  rfl

/-- Observer expectations are invariant under the observer-local isotony
inclusions, again definitionally. -/
@[simp]
theorem observerExpectation_inclusion (o : T.Observer r) {U V : N.Region r}
    (hUV : N.regionLE r U V) (X : A.observerLocalAlgebra o U) :
    A.observerExpectation o V
        (StarSubalgebra.inclusion (A.observerLocal_isotony o hUV) X) =
      A.observerExpectation o U X :=
  rfl

/-- Inheritance, repair fixed points: the regional repair fixes every
observer-local element of its own region. -/
theorem repair_fixes_observerLocal (o : T.Observer r) (U : N.Region r)
    {X : ConsensusTower.PrivateAlgebra T r}
    (hX : X ∈ A.observerLocalAlgebra o U) :
    N.repair r U X = X :=
  N.repair_fixes_region r U ⟨X, (A.mem_observerLocalAlgebra.mp hX).1⟩

/-- Inheritance, remote nondisturbance: repair in a declared-disjoint
region fixes every observer-local element. -/
theorem repair_fixes_observerLocal_of_disjoint (o : T.Observer r)
    {U V : N.Region r} (hUV : N.disjoint r U V)
    {X : ConsensusTower.PrivateAlgebra T r}
    (hX : X ∈ A.observerLocalAlgebra o V) :
    N.repair r U X = X :=
  N.repair_fixes_disjoint r hUV ⟨X, (A.mem_observerLocalAlgebra.mp hX).1⟩

/-! ## The accessible-glue receipt and observer-local descent -/

/-- Deliverable-8 receipt: on every declared descent family inside a region
below the observer region, any ambient globalization of a compatible
family of accessible sections is accessible.  Under a declared descent
packet the globalization is unique, so this says exactly that the selected
ambient glue lies in the accessible algebra.  The receipt is a premise
independent of ambient descent: `QFT.GlueCountermodel` exhibits a net
and cut where unique descent holds while the unique glue of an
accessible family escapes the accessible algebra, and the conditional
theorem below only consumes the receipt. -/
structure AccessibleGlueClosure (o : T.Observer r)
    (D : N.ObserverNetDescent) : Prop where
  glue_accessible : ∀ {W : N.Region r},
    N.regionLE r W (A.observerRegion o) →
    ∀ C : N.FiniteCover r W, D.declaredCover C →
    ∀ F : C.CompatibleFamily,
      (∀ U (hU : U ∈ C.regions),
        (F.localSection U hU : ConsensusTower.PrivateAlgebra T r) ∈
          A.accessibleAlgebra o) →
      ∀ X : N.localAlgebra r W, C.Globalizes F X →
        (X : ConsensusTower.PrivateAlgebra T r) ∈ A.accessibleAlgebra o

/-- Observer-local unique descent, conditional on the accessible-glue
receipt: a compatible family of accessible sections on a declared family
below the observer region has exactly one ambient globalization, and that
globalization is observer-local at the ambient region.  Existence in the
observer-local algebra comes from the receipt; uniqueness is inherited
from ambient uniqueness. -/
theorem observerLocal_hasUniqueDescent {o : T.Observer r}
    {D : N.ObserverNetDescent} (G : A.AccessibleGlueClosure o D)
    {W : N.Region r} (hW : N.regionLE r W (A.observerRegion o))
    (C : N.FiniteCover r W) (hC : D.declaredCover C)
    (F : C.CompatibleFamily)
    (hF : ∀ U (hU : U ∈ C.regions),
      (F.localSection U hU : ConsensusTower.PrivateAlgebra T r) ∈
        A.accessibleAlgebra o) :
    ∃! X : N.localAlgebra r W, C.Globalizes F X ∧
      (X : ConsensusTower.PrivateAlgebra T r) ∈
        A.observerLocalAlgebra o W := by
  obtain ⟨X, hX, huniq⟩ := D.unique_descent C hC F
  refine ⟨X, ⟨hX, A.mem_observerLocalAlgebra.mpr
    ⟨X.2, G.glue_accessible hW C hC F hF X hX⟩⟩, ?_⟩
  rintro Y ⟨hY, -⟩
  exact huniq Y hY

/-- The universal singleton descent packet satisfies the accessible-glue
receipt for every access cut: the singleton glue is the declared section
itself.  This inhabits the receipt type; it carries no nontrivial gluing
content. -/
theorem accessibleGlueClosure_singleton (o : T.Observer r) :
    A.AccessibleGlueClosure o N.singletonDescent := by
  constructor
  intro W hW C hC F hF X hX
  have hCs : C = N.singletonCover r W := hC
  subst hCs
  have hWmem : W ∈ (N.singletonCover r W).regions := by
    simp [FiniteCausalObserverNet.singletonCover]
  have hXW : X = F.localSection W hWmem :=
    (N.restrict_refl r W X).symm.trans (hX W hWmem)
  rw [hXW]
  exact hF W hWmem

end ObserverAccessCut

namespace RepairStableAccessCut

variable {T : ConsensusTower ι} {N : FiniteCausalObserverNet T} {r : ι}
variable {A : ObserverAccessCut T N r}

/-- Given the repair-stability receipt, every B2 relaxation of a repair
below the observer region preserves accessibility: the relaxed element is
an algebraic combination of the element and its repair. -/
theorem relaxedRepair_mem (hA : RepairStableAccessCut A) (o : T.Observer r)
    {U : N.Region r} (hU : N.regionLE r U (A.observerRegion o)) (a : ℂ)
    {X : ConsensusTower.PrivateAlgebra T r}
    (hX : X ∈ A.accessibleAlgebra o) :
    N.relaxedRepair a r U X ∈ A.accessibleAlgebra o := by
  have hrep : N.repair r U X ∈ A.accessibleAlgebra o :=
    hA.repair_preserves o hU X hX
  have hgoal : N.repair r U X + a • (X - N.repair r U X) ∈
      A.accessibleAlgebra o :=
    add_mem hrep ((A.accessibleAlgebra o).smul_mem (sub_mem hX hrep) a)
  simpa [FiniteCausalObserverNet.relaxedRepair, EventAlgebra.publicRelax]
    using hgoal

end RepairStableAccessCut

/-! ## A strictly interposed dimension-two witness

The witness data: the one-block partition supplies the tower and its
scalar public algebra, the full diagonal partition supplies the constant
regional algebra of the net, and the access cut declares the diagonal
algebra accessible.  The two strictness lemmas exhibit the proper chain
public < accessible < private. -/

section AccessWitness

open Matrix
open EventAlgebra

/-- The one-block projective partition: the single sure event. -/
noncomputable def oneBlockPartition (n : ℕ) : ProjectivePartition n 1 where
  proj := fun _ => 1
  isEvent := fun _ => isEvent_one
  orthogonal := fun i j hij => absurd (Subsingleton.elim i j) hij
  complete := by simp

/-- The full diagonal projective partition: one rank-one diagonal
projector per basis index. -/
noncomputable def diagonalPartition (n : ℕ) : ProjectivePartition n n where
  proj := fun i => Matrix.diagonal (Pi.single i 1)
  isEvent := fun i => by
    refine ⟨?_, ?_⟩
    · show (Matrix.diagonal (Pi.single i 1))ᴴ = Matrix.diagonal (Pi.single i 1)
      rw [Matrix.diagonal_conjTranspose]
      congr 1
      funext j
      by_cases h : j = i <;> simp [Pi.single_apply, h]
    · rw [Matrix.diagonal_mul_diagonal]
      congr 1
      funext j
      by_cases h : j = i <;> simp [Pi.single_apply, h]
  orthogonal := fun i j hij => by
    rw [Matrix.diagonal_mul_diagonal]
    have hsingle : (fun x => Pi.single i (1 : ℂ) x * Pi.single j 1 x) =
        fun _ : Fin n => (0 : ℂ) := by
      funext x
      by_cases hxi : x = i
      · subst hxi
        simp [hij]
      · simp [Pi.single_apply, hxi]
    rw [hsingle]
    exact Matrix.diagonal_zero
  complete := by
    ext i j
    by_cases h : i = j
    · subst h
      simp [Matrix.sum_apply, Matrix.diagonal_apply_eq, Pi.single_apply]
    · simp [Matrix.sum_apply, Matrix.diagonal_apply_ne _ h,
        Matrix.one_apply_ne h]

/-- Membership in the one-block public algebra is exactly being a scalar
multiple of the identity. -/
theorem oneBlockPartition_mem_iff {n : ℕ} (X : Matrix (Fin n) (Fin n) ℂ) :
    X ∈ (oneBlockPartition n).publicSubalgebra ↔
      ∃ c : ℂ, c • (1 : Matrix (Fin n) (Fin n) ℂ) = X := by
  constructor
  · intro hX
    have hle : (oneBlockPartition n).span ≤
        Submodule.span ℂ {(1 : Matrix (Fin n) (Fin n) ℂ)} :=
      Submodule.span_mono (by rintro Y ⟨i, rfl⟩; exact rfl)
    exact Submodule.mem_span_singleton.mp (hle hX)
  · rintro ⟨c, rfl⟩
    exact ((oneBlockPartition n).publicSubalgebra).smul_mem
      (one_mem _) c

/-- Every star subalgebra contains the one-block public scalars. -/
theorem oneBlockPartition_publicSubalgebra_le {n : ℕ}
    (S : StarSubalgebra ℂ (Matrix (Fin n) (Fin n) ℂ)) :
    (oneBlockPartition n).publicSubalgebra ≤ S := by
  intro X hX
  obtain ⟨c, rfl⟩ := (oneBlockPartition_mem_iff X).mp hX
  exact S.smul_mem (one_mem S) c

/-- The selected witness state: the pure state on the first basis
direction, certified through the diagonal partition's first projector. -/
noncomputable def witnessState : StateMatrix 2 :=
  ⟨Matrix.diagonal (Pi.single 0 1),
    ((diagonalPartition 2).isEvent 0).posSemidef, by
      rw [Matrix.trace_diagonal]
      simp [Pi.single_apply]⟩

/-- A finite region-indexed net on the constant tower of one partition
with every regional algebra supplied by a second partition.  Restrictions
are the identity and repair is the second partition's average.  The two
partitions are independent parameters, so the tower's public algebra and
the regional algebras can differ; the access-cut witness uses exactly this
freedom. -/
noncomputable def twoPartitionCausalNet {n k m : ℕ}
    (partPub : ProjectivePartition n k) (partLoc : ProjectivePartition n m)
    (rho : StateMatrix n) :
    FiniteCausalObserverNet
      (ConsensusTower.constantConsensusTower partPub rho) where
  Region := fun _ => Finset (Fin 2)
  regionFintype := fun _ => by classical exact inferInstance
  regionNonempty := fun _ => ⟨∅⟩
  regionLE := fun _ U V => U ⊆ V
  regionLE_refl := fun _ _ => Finset.Subset.rfl
  regionLE_trans := by
    intro r U V W hUV hVW
    exact fun x hx => hVW (hUV hx)
  regionLE_antisymm := by
    intro r U V hUV hVU
    exact hUV.antisymm hVU
  overlap := fun _ U V => U ∩ V
  overlap_le_left := fun _ _ _ => Finset.inter_subset_left
  overlap_le_right := fun _ _ _ => Finset.inter_subset_right
  le_overlap := fun _ _ _ _ hWU hWV => Finset.subset_inter hWU hWV
  disjoint := fun _ U V => U.Nonempty ∧ V.Nonempty ∧ Disjoint U V
  disjoint_symm := by
    intro r U V h
    exact ⟨h.2.1, h.1, h.2.2.symm⟩
  disjoint_irrefl := by
    intro r U h
    obtain ⟨x, hx⟩ := h.1
    exact (Finset.disjoint_left.mp h.2.2 hx) hx
  localAlgebra := fun _ _ => partLoc.publicSubalgebra
  isotony := fun _ _ _ _ X hX => hX
  locality := by
    intro r U V hUV X Y
    exact congrArg Subtype.val (partLoc.publicSubalgebra_mul_comm X Y)
  restrict := fun _ {_ _} _ => StarAlgHom.id ℂ partLoc.publicSubalgebra
  restrict_refl := by intros; rfl
  restrict_trans := by intros; rfl
  restrict_inclusion := by intros; rfl
  regionRefine := fun _ U => U
  region_refine_refl := by intros; rfl
  region_refine_trans := by intros; rfl
  region_refine_mono := by intros; assumption
  overlap_natural := by intros; rfl
  disjoint_natural := by intros; assumption
  localAlgebra_natural := by intros; assumption
  repair := fun _ _ => partitionAverageLinearMap partLoc
  repair_idempotent := by
    intro r U X
    exact partitionAverage_idem partLoc X
  repair_fixes_region := by
    intro r U X
    exact partitionAverage_fixes partLoc X.2
  repair_fixes_disjoint := by
    intro r U V hUV X
    exact partitionAverage_fixes partLoc X.2
  repair_natural := by intros; rfl

/-- The witness tower: constant tower of the one-block partition at
dimension two, so the public algebra is the scalar span of the
identity. -/
noncomputable def witnessTower : ConsensusTower Unit :=
  ConsensusTower.constantConsensusTower (oneBlockPartition 2) witnessState

/-- The witness net over the witness tower: every regional algebra is the
full diagonal algebra of the diagonal partition. -/
noncomputable def witnessNet : FiniteCausalObserverNet witnessTower :=
  twoPartitionCausalNet (oneBlockPartition 2) (diagonalPartition 2)
    witnessState

/-- The witness net keeps a declared pair of distinct nonempty disjoint
region labels. -/
theorem witnessNet_has_disjoint_pair :
    witnessNet.disjoint () ({0} : Finset (Fin 2)) ({1} : Finset (Fin 2)) := by
  simp [witnessNet, twoPartitionCausalNet]

/-- The witness access cut: the observer region is the full two-label
region and the accessible algebra is the diagonal algebra. -/
noncomputable def witnessAccessCut : ObserverAccessCut witnessTower witnessNet () where
  observerRegion := fun _ => ({0, 1} : Finset (Fin 2))
  accessibleAlgebra := fun _ => (diagonalPartition 2).publicSubalgebra
  accessible_le_region := fun _ => le_rfl
  public_le_accessible := fun _ =>
    oneBlockPartition_publicSubalgebra_le _
  restrict_preserves := fun _ _ _ _ _ _ hX => hX

/-- The first diagonal projector is accessible to the witness observer. -/
theorem witness_proj_mem_accessible :
    (diagonalPartition 2).proj 0 ∈ (diagonalPartition 2).publicSubalgebra :=
  (diagonalPartition 2).proj_mem_span 0

/-- The first diagonal projector is outside the scalar public algebra. -/
theorem witness_proj_not_mem_public :
    (diagonalPartition 2).proj 0 ∉ (oneBlockPartition 2).publicSubalgebra := by
  intro hmem
  obtain ⟨c, hc⟩ := (oneBlockPartition_mem_iff _).mp hmem
  have h00 := congrFun (congrFun hc 0) 0
  have h11 := congrFun (congrFun hc 1) 1
  simp [diagonalPartition, Matrix.smul_apply] at h00 h11
  rw [h00] at h11
  exact one_ne_zero h11

/-- The raising matrix unit is outside the diagonal accessible algebra. -/
theorem witness_m2Raise_not_mem_accessible :
    FiniteCausalObserverNet.m2Raise ∉
      (diagonalPartition 2).publicSubalgebra := by
  intro hmem
  have hmem' : FiniteCausalObserverNet.m2Raise ∈
      Submodule.span ℂ (Set.range (diagonalPartition 2).proj) := hmem
  obtain ⟨c, hc⟩ :=
    (Submodule.mem_span_range_iff_exists_fun ℂ).mp hmem'
  have h01 := congrFun (congrFun hc 0) 1
  simp [diagonalPartition, Matrix.sum_apply, Matrix.smul_apply,
    FiniteCausalObserverNet.m2Raise] at h01

/-- Strictness, lower interposition: the tower public algebra is strictly
inside the witness accessible algebra. -/
theorem witnessAccessCut_public_lt_accessible
    (o : witnessTower.Observer ()) :
    witnessTower.publicAlgebra () o < witnessAccessCut.accessibleAlgebra o := by
  rw [SetLike.lt_iff_le_and_exists]
  exact ⟨witnessAccessCut.public_le_accessible o,
    (diagonalPartition 2).proj 0, witness_proj_mem_accessible,
    witness_proj_not_mem_public⟩

/-- Strictness, upper interposition: the witness accessible algebra is a
proper subalgebra of the private matrix algebra. -/
theorem witnessAccessCut_accessible_lt_top
    (o : witnessTower.Observer ()) :
    witnessAccessCut.accessibleAlgebra o <
      (⊤ : StarSubalgebra ℂ (ConsensusTower.PrivateAlgebra witnessTower ())) := by
  rw [SetLike.lt_iff_le_and_exists]
  exact ⟨le_top, FiniteCausalObserverNet.m2Raise,
    StarSubalgebra.mem_top, witness_m2Raise_not_mem_accessible⟩

/-- The witness cut satisfies the optional repair-stability receipt: the
diagonal partition average maps every private element into the diagonal
algebra. -/
theorem witnessAccessCut_repairStable :
    RepairStableAccessCut witnessAccessCut := by
  constructor
  intro o U _hU X _hX
  exact partitionAverage_mem_span (diagonalPartition 2) X

/-- Every nonempty subregion family of the witness net has unique descent:
all regional algebras and restriction maps in the witness model
coincide. -/
noncomputable def witnessAllCoversDescent : witnessNet.ObserverNetDescent where
  declaredCover := fun _ => True
  cover_exists := by
    intro r W
    exact ⟨witnessNet.singletonCover r W, trivial⟩
  unique_descent := by
    intro r W C _hC F
    obtain ⟨U, hU⟩ := C.nonempty
    refine ⟨F.localSection U hU, ?_, ?_⟩
    · intro V hV
      have hcompat := F.overlap_compatible U hU V hV
      simpa [witnessNet, twoPartitionCausalNet] using hcompat
    · intro X hX
      have hx := hX U hU
      simpa [witnessNet, twoPartitionCausalNet] using hx

/-- The witness cut satisfies the accessible-glue receipt over the
all-families descent packet of the witness net: every globalization is
one of the declared accessible sections. -/
theorem witnessAccessCut_accessibleGlueClosure
    (o : witnessTower.Observer ()) :
    witnessAccessCut.AccessibleGlueClosure o witnessAllCoversDescent := by
  constructor
  intro W hW C hC F hF X hX
  obtain ⟨U, hU⟩ := C.nonempty
  have hXU : X = F.localSection U hU := hX U hU
  rw [hXU]
  exact hF U hU

end AccessWitness

/-! ## A dimension-four witness with a proper observer-local meet

The dimension-two cut above declares the accessible algebra equal to
every regional algebra, so its observer-local meets collapse onto the
regional algebras.  The witness in this section separates the two meet
components on dimension four, in the block-algebra shape of the character
net of `QFT/NoncommutativeWitness.lean`; that module imports this one, so
the block data is declared here under its own names.  The tower is the
constant tower of the full diagonal partition, so the public algebra is
the diagonal span and nonscalar public records exist.  The net declares
three regions: a scalar bottom region, the diagonal observer region, and
a block region carrying the noncommutative algebra of matrices
`diag(a, a, A)`.  The cut declares the diagonal algebra accessible.  The
observer-local meet at the block region is computed exactly: it is the
algebra of diagonal matrices with equal first two entries, strictly below
the block regional algebra and strictly below the accessible diagonal
algebra. -/

section ProperMeetWitness

open Matrix
open EventAlgebra

/-- The sum of coefficient-weighted diagonal projectors is the diagonal
matrix of the coefficient vector. -/
theorem sum_smul_diagonalPartition_proj {n : ℕ} (c : Fin n → ℂ) :
    ∑ i, c i • (diagonalPartition n).proj i = Matrix.diagonal c := by
  ext a b
  by_cases hab : a = b
  · subst hab
    simp [diagonalPartition, Matrix.sum_apply, Matrix.smul_apply,
      Pi.single_apply]
  · simp [diagonalPartition, Matrix.sum_apply, Matrix.smul_apply,
      Matrix.diagonal_apply_ne _ hab]

/-- Membership in the diagonal-partition public algebra is exactly being
a diagonal matrix. -/
theorem diagonalPartition_mem_iff {n : ℕ} {X : Matrix (Fin n) (Fin n) ℂ} :
    X ∈ (diagonalPartition n).publicSubalgebra ↔
      ∃ c : Fin n → ℂ, Matrix.diagonal c = X := by
  constructor
  · intro hX
    have hX' : X ∈ Submodule.span ℂ (Set.range (diagonalPartition n).proj) :=
      hX
    obtain ⟨c, hc⟩ := (Submodule.mem_span_range_iff_exists_fun ℂ).mp hX'
    exact ⟨c, (sum_smul_diagonalPartition_proj c).symm.trans hc⟩
  · rintro ⟨c, rfl⟩
    rw [← sum_smul_diagonalPartition_proj c]
    exact Submodule.sum_mem _ fun i _ =>
      Submodule.smul_mem _ _ ((diagonalPartition n).proj_mem_span i)

/-- Diagonal evaluation of a scalar public element reproduces the
element. -/
theorem oneBlockPartition_corner_eval {n : ℕ}
    (X : (oneBlockPartition n).publicSubalgebra) (i : Fin n) :
    (X : Matrix (Fin n) (Fin n) ℂ) i i • (1 : Matrix (Fin n) (Fin n) ℂ) =
      (X : Matrix (Fin n) (Fin n) ℂ) := by
  obtain ⟨c, hc⟩ := (oneBlockPartition_mem_iff _).mp X.2
  rw [← hc]
  simp [Matrix.smul_apply]

/-- The block algebra of matrices `diag(a, a, A)` inside `M₄(ℂ)`: scalar
upper block, free lower two-by-two block, vanishing mixed blocks.  The
carrier repeats the left-block shape of `QFT/NoncommutativeWitness.lean`,
which imports this module and therefore cannot be imported here. -/
def properMeetBlockAlgebra : StarSubalgebra ℂ (Matrix (Fin 4) (Fin 4) ℂ) where
  carrier := {X | X 0 1 = 0 ∧ X 1 0 = 0 ∧ X 1 1 = X 0 0 ∧
    X 0 2 = 0 ∧ X 0 3 = 0 ∧ X 1 2 = 0 ∧ X 1 3 = 0 ∧
    X 2 0 = 0 ∧ X 2 1 = 0 ∧ X 3 0 = 0 ∧ X 3 1 = 0}
  zero_mem' := by simp [Set.mem_setOf_eq]
  one_mem' := by simp [Set.mem_setOf_eq]
  add_mem' := by
    intro X Y hX hY
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ := hX
    obtain ⟨g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11⟩ := hY
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
      simp [Matrix.add_apply, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11,
        g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11]
  mul_mem' := by
    intro X Y hX hY
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ := hX
    obtain ⟨g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11⟩ := hY
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
      simp [Matrix.mul_apply, Fin.sum_univ_four, h1, h2, h3, h4, h5, h6, h7,
        h8, h9, h10, h11, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11]
  algebraMap_mem' := by
    intro c
    simp [Set.mem_setOf_eq, Algebra.algebraMap_eq_smul_one,
      Matrix.smul_apply]
  star_mem' := by
    intro X hX
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ := hX
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
      simp [h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11]

theorem mem_properMeetBlockAlgebra_iff {X : Matrix (Fin 4) (Fin 4) ℂ} :
    X ∈ properMeetBlockAlgebra ↔
      (X 0 1 = 0 ∧ X 1 0 = 0 ∧ X 1 1 = X 0 0 ∧
        X 0 2 = 0 ∧ X 0 3 = 0 ∧ X 1 2 = 0 ∧ X 1 3 = 0 ∧
        X 2 0 = 0 ∧ X 2 1 = 0 ∧ X 3 0 = 0 ∧ X 3 1 = 0) :=
  Iff.rfl

/-- The lower-block raising unit `|2><3|`. -/
def properMeetRaise : Matrix (Fin 4) (Fin 4) ℂ :=
  fun k l => if k = 2 ∧ l = 3 then 1 else 0

/-- The lower-block lowering unit `|3><2|`. -/
def properMeetLower : Matrix (Fin 4) (Fin 4) ℂ :=
  fun k l => if k = 3 ∧ l = 2 then 1 else 0

theorem properMeetRaise_mem_block :
    properMeetRaise ∈ properMeetBlockAlgebra := by
  simp [mem_properMeetBlockAlgebra_iff, properMeetRaise]

theorem properMeetLower_mem_block :
    properMeetLower ∈ properMeetBlockAlgebra := by
  simp [mem_properMeetBlockAlgebra_iff, properMeetLower]

/-- The block regional algebra is noncommutative: its two lower-block
matrix units have distinct products. -/
theorem properMeetBlockAlgebra_noncommutative :
    ∃ X Y : properMeetBlockAlgebra,
      ¬ Commute (X : Matrix (Fin 4) (Fin 4) ℂ)
        (Y : Matrix (Fin 4) (Fin 4) ℂ) := by
  refine ⟨⟨properMeetRaise, properMeetRaise_mem_block⟩,
    ⟨properMeetLower, properMeetLower_mem_block⟩, fun h => ?_⟩
  have h22 := congrFun (congrFun h.eq 2) 2
  simp [properMeetRaise, properMeetLower, Matrix.mul_apply] at h22

/-- The lower-block raising unit is outside the diagonal algebra. -/
theorem properMeetRaise_not_mem_diagonal :
    properMeetRaise ∉ (diagonalPartition 4).publicSubalgebra := by
  intro hmem
  obtain ⟨c, hc⟩ := diagonalPartition_mem_iff.mp hmem
  have h23 := congrFun (congrFun hc 2) 3
  simp [properMeetRaise] at h23

/-- The second diagonal projector is outside the block algebra. -/
theorem diagonalProj_not_mem_properMeetBlock :
    (diagonalPartition 4).proj 1 ∉ properMeetBlockAlgebra := by
  intro hmem
  have h := (mem_properMeetBlockAlgebra_iff.mp hmem).2.2.1
  simp [diagonalPartition] at h

/-- The selected dimension-four witness state: the pure state on the
first basis direction. -/
noncomputable def properMeetState : StateMatrix 4 :=
  ⟨Matrix.diagonal (Pi.single 0 1),
    ((diagonalPartition 4).isEvent 0).posSemidef, by
      rw [Matrix.trace_diagonal]
      simp [Pi.single_apply]⟩

/-- The dimension-four witness tower: the constant tower of the full
diagonal partition, so the public algebra is the whole diagonal span and
nonscalar public records exist. -/
noncomputable def properMeetTower : ConsensusTower Unit :=
  ConsensusTower.constantConsensusTower (diagonalPartition 4) properMeetState

/-- Region labels for the proper-meet net: a scalar bottom region, the
diagonal observer region, and the noncommutative block region. -/
inductive ProperMeetRegion : Type
  | bot
  | diag
  | block
  deriving DecidableEq, Fintype

namespace ProperMeetRegion

/-- Region inclusion: reflexive, and the bottom region is below every
region. -/
def le (U V : ProperMeetRegion) : Prop := U = V ∨ U = bot

instance : DecidableRel le := fun U V =>
  decidable_of_iff (U = V ∨ U = bot) Iff.rfl

/-- Declared overlap: the common region on the diagonal of the relation
and the bottom region for distinct pairs. -/
def meet (U V : ProperMeetRegion) : ProperMeetRegion :=
  if U = V then U else bot

theorem le_refl : ∀ U, le U U := by decide

theorem le_trans : ∀ U V W, le U V → le V W → le U W := by decide

theorem le_antisymm : ∀ U V, le U V → le V U → U = V := by decide

theorem meet_le_left : ∀ U V, le (meet U V) U := by decide

theorem meet_le_right : ∀ U V, le (meet U V) V := by decide

theorem le_meet : ∀ W U V, le W U → le W V → le W (meet U V) := by decide

/-- The block region is outside the cone below the diagonal region, so
the containment laws stated at the diagonal region are falsifiable. -/
theorem block_not_le_diag : ¬ le block diag := by decide

end ProperMeetRegion

/-- The regional algebra assignment of the proper-meet net: scalars at
the bottom region, the diagonal algebra at the diagonal region, and the
block algebra at the block region. -/
noncomputable def properMeetRegionAlgebra :
    ProperMeetRegion → StarSubalgebra ℂ (Matrix (Fin 4) (Fin 4) ℂ)
  | .bot => (oneBlockPartition 4).publicSubalgebra
  | .diag => (diagonalPartition 4).publicSubalgebra
  | .block => properMeetBlockAlgebra

/-- Scalar multiples of the identity are one-block public elements. -/
theorem properMeet_smul_one_mem (c : ℂ) :
    c • (1 : Matrix (Fin 4) (Fin 4) ℂ) ∈
      (oneBlockPartition 4).publicSubalgebra :=
  (oneBlockPartition_mem_iff _).mpr ⟨c, rfl⟩

/-- The character of the diagonal algebra: evaluation of the first
diagonal entry, returned as a scalar public element. -/
noncomputable def properMeetDiagCharacter :
    (diagonalPartition 4).publicSubalgebra →⋆ₐ[ℂ]
      (oneBlockPartition 4).publicSubalgebra where
  toFun X := ⟨(X : Matrix (Fin 4) (Fin 4) ℂ) 0 0 • 1,
    properMeet_smul_one_mem _⟩
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

/-- The character of the block algebra: evaluation of the scalar block,
returned as a scalar public element. -/
noncomputable def properMeetBlockCharacter :
    properMeetBlockAlgebra →⋆ₐ[ℂ] (oneBlockPartition 4).publicSubalgebra where
  toFun X := ⟨(X : Matrix (Fin 4) (Fin 4) ℂ) 0 0 • 1,
    properMeet_smul_one_mem _⟩
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
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11⟩ :=
      mem_properMeetBlockAlgebra_iff.mp X.2
    simp [Matrix.mul_apply, Fin.sum_univ_four, h1, h4, h5, smul_smul,
      mul_comm]
  commutes' c := by
    apply Subtype.ext
    simp [Algebra.algebraMap_eq_smul_one, Matrix.smul_apply]
  map_star' X := by
    apply Subtype.ext
    simp [star_smul]

/-- The supplied restriction system of the proper-meet net: identities on
equal regions and the two characters on the drops to the bottom region.
The remaining constructor pairs are excluded by the region order. -/
noncomputable def properMeetRestrict :
    ∀ U V : ProperMeetRegion, ProperMeetRegion.le V U →
      (properMeetRegionAlgebra U →⋆ₐ[ℂ] properMeetRegionAlgebra V)
  | .bot, .bot, _ => StarAlgHom.id ℂ _
  | .diag, .diag, _ => StarAlgHom.id ℂ _
  | .block, .block, _ => StarAlgHom.id ℂ _
  | .diag, .bot, _ => properMeetDiagCharacter
  | .block, .bot, _ => properMeetBlockCharacter
  | .bot, .diag, h => absurd h (by decide)
  | .bot, .block, h => absurd h (by decide)
  | .diag, .block, h => absurd h (by decide)
  | .block, .diag, h => absurd h (by decide)

/-- The proper-meet causal net over the dimension-four tower: scalar
bottom region, diagonal region, and noncommutative block region, with
character restrictions onto the bottom region.  The declared disjointness
relation is empty and the repair maps are the identity; the interface
laws hold for them, and nontrivial repair semantics belong to E2. -/
noncomputable def properMeetNet : FiniteCausalObserverNet properMeetTower where
  Region := fun _ => ProperMeetRegion
  regionFintype := fun _ => inferInstance
  regionNonempty := fun _ => ⟨ProperMeetRegion.bot⟩
  regionLE := fun _ => ProperMeetRegion.le
  regionLE_refl := fun _ U => ProperMeetRegion.le_refl U
  regionLE_trans := fun _ {U V W} hUV hVW =>
    ProperMeetRegion.le_trans U V W hUV hVW
  regionLE_antisymm := fun _ {U V} hUV hVU =>
    ProperMeetRegion.le_antisymm U V hUV hVU
  overlap := fun _ => ProperMeetRegion.meet
  overlap_le_left := fun _ U V => ProperMeetRegion.meet_le_left U V
  overlap_le_right := fun _ U V => ProperMeetRegion.meet_le_right U V
  le_overlap := fun _ {W U V} hWU hWV =>
    ProperMeetRegion.le_meet W U V hWU hWV
  disjoint := fun _ _ _ => False
  disjoint_symm := by
    intro r U V h
    exact h.elim
  disjoint_irrefl := fun _ _ h => h
  localAlgebra := fun _ => properMeetRegionAlgebra
  isotony := by
    intro r U V hUV
    rcases hUV with rfl | rfl
    · exact le_rfl
    · exact oneBlockPartition_publicSubalgebra_le _
  locality := by
    intro r U V h
    exact h.elim
  restrict := fun _ {U V} h => properMeetRestrict U V h
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
    cases U with
    | bot =>
        cases V with
        | bot => rfl
        | diag =>
            apply Subtype.ext
            exact oneBlockPartition_corner_eval X (0 : Fin 4)
        | block =>
            apply Subtype.ext
            exact oneBlockPartition_corner_eval X (0 : Fin 4)
    | diag =>
        cases V with
        | bot => exact absurd hUV (by decide)
        | diag => rfl
        | block => exact absurd hUV (by decide)
    | block =>
        cases V with
        | bot => exact absurd hUV (by decide)
        | diag => exact absurd hUV (by decide)
        | block => rfl
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

/-- The proper-meet access cut: the observer region is the diagonal
region and the accessible algebra is the diagonal algebra.  At the block
region both meet components differ from the meet; the two strictness
receipts below certify this. -/
noncomputable def properMeetAccessCut :
    ObserverAccessCut properMeetTower properMeetNet () where
  observerRegion := fun _ => ProperMeetRegion.diag
  accessibleAlgebra := fun _ => (diagonalPartition 4).publicSubalgebra
  accessible_le_region := fun _ => le_rfl
  public_le_accessible := fun _ => le_rfl
  restrict_preserves := by
    intro o U V hU hVU X hX
    cases U with
    | bot =>
        cases V with
        | bot => exact hX
        | diag =>
            exact absurd
              (show ProperMeetRegion.le ProperMeetRegion.diag
                ProperMeetRegion.bot from hVU) (by decide)
        | block =>
            exact absurd
              (show ProperMeetRegion.le ProperMeetRegion.block
                ProperMeetRegion.bot from hVU) (by decide)
    | diag =>
        cases V with
        | bot =>
            exact ((diagonalPartition 4).publicSubalgebra).smul_mem
              (one_mem _) _
        | diag => exact hX
        | block =>
            exact absurd
              (show ProperMeetRegion.le ProperMeetRegion.block
                ProperMeetRegion.diag from hVU)
              ProperMeetRegion.block_not_le_diag
    | block =>
        exact absurd
          (show ProperMeetRegion.le ProperMeetRegion.block
            ProperMeetRegion.diag from hU)
          ProperMeetRegion.block_not_le_diag

/-- The observer-local meet at the block region, computed exactly: a
matrix lies in the meet precisely when it is diagonal with equal first
two entries, the intersection `diag(a, a, d, e)` of the diagonal algebra
with the block algebra. -/
theorem properMeetAccessCut_mem_block_meet_iff
    (o : properMeetTower.Observer ()) (X : Matrix (Fin 4) (Fin 4) ℂ) :
    X ∈ properMeetAccessCut.observerLocalAlgebra o ProperMeetRegion.block ↔
      ∃ a d e : ℂ, Matrix.diagonal ![a, a, d, e] = X := by
  constructor
  · intro hX
    obtain ⟨hblock, hdiag⟩ :=
      properMeetAccessCut.mem_observerLocalAlgebra.mp hX
    obtain ⟨c, hc⟩ := diagonalPartition_mem_iff.mp hdiag
    have h11 : X 1 1 = X 0 0 :=
      (mem_properMeetBlockAlgebra_iff.mp hblock).2.2.1
    have hc1 : c 1 = c 0 := by
      rw [← hc] at h11
      simpa using h11
    refine ⟨c 0, c 2, c 3, ?_⟩
    rw [← hc]
    congr 1
    funext i
    fin_cases i <;> simp [hc1]
  · rintro ⟨a, d, e, rfl⟩
    refine properMeetAccessCut.mem_observerLocalAlgebra.mpr ⟨?_, ?_⟩
    · refine mem_properMeetBlockAlgebra_iff.mpr ?_
      refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp
    · exact diagonalPartition_mem_iff.mpr ⟨![a, a, d, e], rfl⟩

/-- Strictness, regional side: the observer-local meet at the block
region is strictly below the block regional algebra; the lower-block
raising unit separates them. -/
theorem properMeetAccessCut_meet_lt_regional
    (o : properMeetTower.Observer ()) :
    properMeetAccessCut.observerLocalAlgebra o ProperMeetRegion.block <
      properMeetNet.localAlgebra () ProperMeetRegion.block := by
  rw [SetLike.lt_iff_le_and_exists]
  refine ⟨properMeetAccessCut.observerLocal_le_local o ProperMeetRegion.block,
    properMeetRaise, properMeetRaise_mem_block, fun hmem => ?_⟩
  exact properMeetRaise_not_mem_diagonal
    (properMeetAccessCut.mem_observerLocalAlgebra.mp hmem).2

/-- Strictness, accessible side: the observer-local meet at the block
region is strictly below the accessible diagonal algebra; the second
diagonal projector separates them. -/
theorem properMeetAccessCut_meet_lt_accessible
    (o : properMeetTower.Observer ()) :
    properMeetAccessCut.observerLocalAlgebra o ProperMeetRegion.block <
      properMeetAccessCut.accessibleAlgebra o := by
  rw [SetLike.lt_iff_le_and_exists]
  refine
    ⟨properMeetAccessCut.observerLocal_le_accessible o ProperMeetRegion.block,
      (diagonalPartition 4).proj 1, (diagonalPartition 4).proj_mem_span 1,
      fun hmem => ?_⟩
  exact diagonalProj_not_mem_properMeetBlock
    (properMeetAccessCut.mem_observerLocalAlgebra.mp hmem).1

end ProperMeetWitness

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.ObserverAccessCut.observerLocalAlgebra_observerRegion
#print axioms OPH.QFT.ObserverAccessCut.observerLocal_isotony
#print axioms OPH.QFT.ObserverAccessCut.observerLocal_commute
#print axioms OPH.QFT.ObserverAccessCut.restrict_mem_observerLocal
#print axioms OPH.QFT.ObserverAccessCut.regionalExpectation_observerLocal
#print axioms OPH.QFT.ObserverAccessCut.observerExpectation_inclusion
#print axioms OPH.QFT.ObserverAccessCut.repair_fixes_observerLocal
#print axioms OPH.QFT.ObserverAccessCut.repair_fixes_observerLocal_of_disjoint
#print axioms OPH.QFT.ObserverAccessCut.observerLocal_hasUniqueDescent
#print axioms OPH.QFT.ObserverAccessCut.accessibleGlueClosure_singleton
#print axioms OPH.QFT.RepairStableAccessCut.relaxedRepair_mem
#print axioms OPH.QFT.witnessNet_has_disjoint_pair
#print axioms OPH.QFT.witnessAccessCut_public_lt_accessible
#print axioms OPH.QFT.witnessAccessCut_accessible_lt_top
#print axioms OPH.QFT.witnessAccessCut_repairStable
#print axioms OPH.QFT.witnessAccessCut_accessibleGlueClosure
#print axioms OPH.QFT.diagonalPartition_mem_iff
#print axioms OPH.QFT.properMeetBlockAlgebra_noncommutative
#print axioms OPH.QFT.properMeetAccessCut_mem_block_meet_iff
#print axioms OPH.QFT.properMeetAccessCut_meet_lt_regional
#print axioms OPH.QFT.properMeetAccessCut_meet_lt_accessible
