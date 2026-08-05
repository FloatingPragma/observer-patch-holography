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
into the observer-local algebra.  Ambient unique descent alone does not
place the glue in the accessible algebra; the receipt is a separate
premise, and no theorem in this module derives it from descent.  The
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
private matrix algebra.

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
ambient glue lies in the accessible algebra.  The receipt is a premise:
ambient unique descent alone does not place the glue in the accessible
algebra, and the conditional theorem below only consumes the receipt. -/
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

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.ObserverAccessCut.observerLocalAlgebra_observerRegion
#print axioms OPH.QFT.ObserverAccessCut.observerLocal_isotony
#print axioms OPH.QFT.ObserverAccessCut.observerLocal_commute
#print axioms OPH.QFT.ObserverAccessCut.restrict_mem_observerLocal
#print axioms OPH.QFT.ObserverAccessCut.regionalExpectation_observerLocal
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
