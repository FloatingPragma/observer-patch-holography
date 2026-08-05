import Tower.ConsensusTower
import EventAlgebra.TwoScalePublicRepair
import ObserverPatchHolography.Locality.NoSignalling

/-!
# Finite causal observer nets

This file implements a proof-carrying finite algebraic packet for
completion-plan issue `#692`.  A `FiniteCausalObserverNet` decorates one A3 `ConsensusTower` with a
finite region system at every regulator, region-indexed star subalgebras of
the tower's private matrix algebra, isotony, algebraic locality on *declared*
disjoint regions, contravariant restriction maps, and idempotent local repair
maps.  The restriction maps are an extra presheaf-like receipt: ordinary
AQFT isotony does not provide a homomorphic retraction from a larger local
algebra to a smaller one.

The repair maps are deliberately only complex-linear idempotents.  Their
locality is an explicit field: a repair in `U` fixes every observable in a
declared-disjoint region `V`.  The B2 `publicRelax` construction then gives an
exact multiplicative-amplitude relaxation law, and remote observables remain
fixed throughout that relaxation.  No CP or CPTP predicate is inferred from
this algebraic statement.

The A3 selected density matrices restrict to compatible expectation
functionals on the local algebras.  The compatibility theorem follows from
the tower's exact trace-pairing naturality.  A generic B4 basis-reindexed
marginal identity is exposed through a supplied tensor-split receipt, but no
theorem identifies its factors with the regional algebras or uses regional
disjointness.

Nothing here constructs spacetime, a Lorentzian causal relation, an adaptive
scheduler, a continuum net, a physical clock, or a source-selected regional
factorization.
-/

namespace OPH.QFT

open Matrix
open OPH.Tower

universe u

variable {ι : Type u} [Preorder ι]

/-- A proof-carrying enriched finite algebraic observer net decorating one
timeless consensus tower.

`regionLE`, `overlap`, and `disjoint` are declared finite-region data.  They
are not interpreted as spacetime relations.  The local algebras all live in
the same finite private matrix algebra at a fixed regulator. -/
structure FiniteCausalObserverNet (T : ConsensusTower ι) where
  /-- Finite observer-accessible regions at each regulator. -/
  Region : ι → Type u
  regionFintype : ∀ r, Fintype (Region r)
  regionNonempty : ∀ r, Nonempty (Region r)
  /-- Region inclusion, bundled with its finite-poset laws. -/
  regionLE : ∀ r, Region r → Region r → Prop
  regionLE_refl : ∀ r U, regionLE r U U
  regionLE_trans : ∀ r {U V W},
    regionLE r U V → regionLE r V W → regionLE r U W
  regionLE_antisymm : ∀ r {U V},
    regionLE r U V → regionLE r V U → U = V
  /-- A declared overlap region with the meet universal property. -/
  overlap : ∀ r, Region r → Region r → Region r
  overlap_le_left : ∀ r U V, regionLE r (overlap r U V) U
  overlap_le_right : ∀ r U V, regionLE r (overlap r U V) V
  le_overlap : ∀ r {W U V},
    regionLE r W U → regionLE r W V → regionLE r W (overlap r U V)
  /-- Declared disjointness; no geometric or relativistic reading is built in. -/
  disjoint : ∀ r, Region r → Region r → Prop
  disjoint_symm : ∀ r {U V}, disjoint r U V → disjoint r V U
  disjoint_irrefl : ∀ r U, ¬ disjoint r U U
  /-- Region-indexed private star subalgebras of the A3 matrix algebra. -/
  localAlgebra : ∀ r, Region r →
    StarSubalgebra ℂ (ConsensusTower.PrivateAlgebra T r)
  /-- Isotony: inclusion of regions gives inclusion of local algebras. -/
  isotony : ∀ r {U V}, regionLE r U V →
    localAlgebra r U ≤ localAlgebra r V
  /-- Algebraic locality on the declared disjointness relation. -/
  locality : ∀ r {U V}, disjoint r U V →
    ∀ (X : localAlgebra r U) (Y : localAlgebra r V),
      Commute (X : ConsensusTower.PrivateAlgebra T r) Y
  /-- Contravariant restriction along a region inclusion.  This is a
  nonstandard extra receipt, not a consequence of AQFT isotony. -/
  restrict : ∀ r {U V}, regionLE r V U →
    localAlgebra r U →⋆ₐ[ℂ] localAlgebra r V
  restrict_refl : ∀ r U (X : localAlgebra r U),
    restrict r (regionLE_refl r U) X = X
  restrict_trans : ∀ r {U V W}
      (hVU : regionLE r V U) (hWV : regionLE r W V)
      (X : localAlgebra r U),
    restrict r hWV (restrict r hVU X) =
      restrict r (regionLE_trans r hWV hVU) X
  /-- Restriction retracts the isotony inclusion. -/
  restrict_inclusion : ∀ r {U V} (hUV : regionLE r U V)
      (X : localAlgebra r U),
    restrict r hUV
      ((StarSubalgebra.inclusion (isotony r hUV)) X) = X
  /-- Covariant region refinement along the A3 regulator preorder. -/
  regionRefine : ∀ {r s}, r ≤ s → Region r → Region s
  region_refine_refl : ∀ r U,
    regionRefine (show r ≤ r from le_rfl) U = U
  region_refine_trans : ∀ {r s t} (hrs : r ≤ s) (hst : s ≤ t) U,
    regionRefine hst (regionRefine hrs U) =
      regionRefine (hrs.trans hst) U
  region_refine_mono : ∀ {r s} (hrs : r ≤ s) {U V},
    regionLE r U V →
      regionLE s (regionRefine hrs U) (regionRefine hrs V)
  overlap_natural : ∀ {r s} (hrs : r ≤ s) U V,
    regionRefine hrs (overlap r U V) =
      overlap s (regionRefine hrs U) (regionRefine hrs V)
  disjoint_natural : ∀ {r s} (hrs : r ≤ s) {U V},
    disjoint r U V → disjoint s (regionRefine hrs U) (regionRefine hrs V)
  /-- A3 algebra refinement maps every local algebra into its refined region. -/
  localAlgebra_natural : ∀ {r s} (hrs : r ≤ s) U X,
    X ∈ localAlgebra r U →
      T.algebraRefine hrs X ∈ localAlgebra s (regionRefine hrs U)
  /-- An algebraic, idempotent repair map localized by explicit fixed-point
  conditions.  It need not be positive, multiplicative, CP, or CPTP. -/
  repair : ∀ r, Region r →
    ConsensusTower.PrivateAlgebra T r →ₗ[ℂ]
      ConsensusTower.PrivateAlgebra T r
  repair_idempotent : ∀ r U X, repair r U (repair r U X) = repair r U X
  repair_fixes_region : ∀ r U
      (X : localAlgebra r U), repair r U X = X
  repair_fixes_disjoint : ∀ r {U V}, disjoint r U V →
    ∀ (X : localAlgebra r V), repair r U X = X
  /-- Repair naturality is an exact intertwining law, not a continuum limit. -/
  repair_natural : ∀ {r s} (hrs : r ≤ s) U X,
    T.algebraRefine hrs (repair r U X) =
      repair s (regionRefine hrs U) (T.algebraRefine hrs X)

attribute [instance] FiniteCausalObserverNet.regionFintype
attribute [instance] FiniteCausalObserverNet.regionNonempty

namespace FiniteCausalObserverNet

variable {T : ConsensusTower ι} (N : FiniteCausalObserverNet T)

/-- Isotony bundled as a star-algebra homomorphism. -/
noncomputable def localInclusion {r : ι} {U V : N.Region r}
    (hUV : N.regionLE r U V) :
    N.localAlgebra r U →⋆ₐ[ℂ] N.localAlgebra r V :=
  StarSubalgebra.inclusion (N.isotony r hUV)

/-- A3 algebra refinement restricted to a regional algebra. -/
noncomputable def localRefine {r s : ι} (hrs : r ≤ s)
    (U : N.Region r) :
    N.localAlgebra r U →⋆ₐ[ℂ]
      N.localAlgebra s (N.regionRefine hrs U) where
  toFun X := ⟨T.algebraRefine hrs X.1,
    N.localAlgebra_natural hrs U X.1 X.2⟩
  map_zero' := by
    apply Subtype.ext
    exact map_zero (T.algebraRefine hrs)
  map_one' := by
    apply Subtype.ext
    exact map_one (T.algebraRefine hrs)
  map_add' X Y := by
    apply Subtype.ext
    exact map_add (T.algebraRefine hrs) X.1 Y.1
  map_mul' X Y := by
    apply Subtype.ext
    exact map_mul (T.algebraRefine hrs) X.1 Y.1
  map_star' X := by
    apply Subtype.ext
    exact map_star (T.algebraRefine hrs) X.1
  commutes' c := by
    apply Subtype.ext
    exact (T.algebraRefine hrs).commutes c

/-- Algebraic locality restated at the ambient-matrix level. -/
theorem commute_of_disjoint {r : ι} {U V : N.Region r}
    (hUV : N.disjoint r U V)
    (X : N.localAlgebra r U) (Y : N.localAlgebra r V) :
    (X : ConsensusTower.PrivateAlgebra T r) * Y =
      Y * (X : ConsensusTower.PrivateAlgebra T r) :=
  (N.locality r hUV X Y).eq

/-- The selected A3 density matrix evaluated on a regional observable. -/
noncomputable def regionalExpectation (r : ι) (o : T.Observer r)
    (U : N.Region r) (X : N.localAlgebra r U) : ℂ :=
  (T.state r o * (X : ConsensusTower.PrivateAlgebra T r)).trace

/-- Isotony does not change the selected-state expectation of an observable. -/
@[simp]
theorem regionalExpectation_inclusion {r : ι} (o : T.Observer r)
    {U V : N.Region r} (hUV : N.regionLE r U V)
    (X : N.localAlgebra r U) :
    N.regionalExpectation r o V (N.localInclusion hUV X) =
      N.regionalExpectation r o U X :=
  rfl

/-- A3 state naturality makes regional expectations compatible with both
regulator and region refinement. -/
theorem regionalExpectation_refine {r s : ι} (hrs : r ≤ s)
    (o : T.Observer r) (U : N.Region r)
    (X : N.localAlgebra r U) :
    N.regionalExpectation s (T.observerRefine hrs o)
        (N.regionRefine hrs U) (N.localRefine hrs U X) =
      N.regionalExpectation r o U X :=
  T.state_natural hrs o X.1

/-- B2 relaxation applied to one declared local repair idempotent. -/
noncomputable def relaxedRepair (a : ℂ) (r : ι) (U : N.Region r)
    (X : ConsensusTower.PrivateAlgebra T r) :
    ConsensusTower.PrivateAlgebra T r :=
  EventAlgebra.publicRelax a (N.repair r U) X

/-- Local repair amplitudes compose multiplicatively, exactly as in B2. -/
theorem relaxedRepair_compose (a b : ℂ) (r : ι) (U : N.Region r)
    (X : ConsensusTower.PrivateAlgebra T r) :
    N.relaxedRepair a r U (N.relaxedRepair b r U X) =
      N.relaxedRepair (a * b) r U X := by
  apply EventAlgebra.publicRelax_compose
  apply LinearMap.ext
  intro Y
  exact N.repair_idempotent r U Y

/-- Every observable in the repaired region is fixed throughout the B2
relaxation. -/
@[simp]
theorem relaxedRepair_fixes_region (a : ℂ) (r : ι) (U : N.Region r)
    (X : N.localAlgebra r U) :
    N.relaxedRepair a r U X = X := by
  simp [relaxedRepair, EventAlgebra.publicRelax, N.repair_fixes_region r U X]

/-- Algebraic no-signalling: repair in `U` fixes every observable in a
declared-disjoint `V`, for every relaxation amplitude.  This is a
Heisenberg-observable identity, not a bundled channel theorem. -/
theorem relaxedRepair_fixes_disjoint {r : ι} {U V : N.Region r}
    (hUV : N.disjoint r U V) (a : ℂ)
    (X : N.localAlgebra r V) :
    N.relaxedRepair a r U X = X := by
  simp [relaxedRepair, EventAlgebra.publicRelax,
    N.repair_fixes_disjoint r hUV X]

/-- Consequently every matrix state, including but not limited to the A3
selected states, assigns the same expectation to a remote observable before
and after the declared local relaxation. -/
theorem relaxedRepair_remote_expectation {r : ι} {U V : N.Region r}
    (hUV : N.disjoint r U V) (a : ℂ)
    (rho : ConsensusTower.PrivateAlgebra T r)
    (X : N.localAlgebra r V) :
    (rho * N.relaxedRepair a r U X).trace = (rho * X.1).trace := by
  rw [N.relaxedRepair_fixes_disjoint hUV a X]

/-! ## B4 basis-split helper carried alongside a region pair -/

/-- A declared basis factorization carried alongside two disjoint regions.
The current fields do not identify either regional local algebra with one of
the two matrix factors, and the marginal theorem below does not use regional
disjointness.  It is therefore a typed home for the generic B4 helper, not yet
a regional-factor adapter. -/
structure TensorSplitReceipt (r : ι) (U V : N.Region r) where
  leftDim : ℕ
  rightDim : ℕ
  basis : Fin (T.dim r) ≃ Fin leftDim × Fin rightDim
  regions_disjoint : N.disjoint r U V

namespace TensorSplitReceipt

variable {N}

/-- Reindex an A3 private matrix in the declared tensor-product basis. -/
noncomputable def splitMatrix {r : ι} {U V : N.Region r}
    (S : N.TensorSplitReceipt r U V)
    (X : ConsensusTower.PrivateAlgebra T r) :
    Matrix (Fin S.leftDim × Fin S.rightDim)
      (Fin S.leftDim × Fin S.rightDim) ℂ :=
  Matrix.reindex S.basis S.basis X

/-- The generic B4 partial-trace identity after the supplied basis reindexing.
Kraus completeness is used; the region labels and `regions_disjoint` field are
not.  A physical or regional no-signalling use needs an additional theorem
identifying the left matrix factor with operations localized in `U` and the
retained factor with `V`. -/
theorem kraus_remote_marginal_invariant {r : ι} {U V : N.Region r}
    (S : N.TensorSplitReceipt r U V)
    {κ : Type*} (s : Finset κ)
    (Kr : κ → Matrix (Fin S.leftDim) (Fin S.leftDim) ℂ)
    (hK : ∑ i ∈ s, (Kr i)ᴴ * Kr i = 1)
    (X : ConsensusTower.PrivateAlgebra T r) :
    OPH.Locality.ptraceFst
        (∑ i ∈ s,
          OPH.Locality.liftLeft (β := Fin S.rightDim) (Kr i) *
            S.splitMatrix X *
              (OPH.Locality.liftLeft (β := Fin S.rightDim) (Kr i))ᴴ) =
      OPH.Locality.ptraceFst (S.splitMatrix X) :=
  OPH.Locality.ptraceFst_local_kraus s Kr hK (S.splitMatrix X)

end TensorSplitReceipt

/-! ## Independence controls -/

/-- Algebraic locality for an unbundled finite region assignment. -/
def IsAlgebraicallyLocal {R : Type*}
    (Disjoint : R → R → Prop)
    (A : R → StarSubalgebra ℂ (Matrix (Fin 2) (Fin 2) ℂ)) : Prop :=
  ∀ {U V}, Disjoint U V → ∀ (X : A U) (Y : A V),
    Commute (X : Matrix (Fin 2) (Fin 2) ℂ) Y

/-- The matrix unit `|0><1|`. -/
def m2Raise : Matrix (Fin 2) (Fin 2) ℂ :=
  fun i j => if i = 0 ∧ j = 1 then 1 else 0

/-- The matrix unit `|1><0|`. -/
def m2Lower : Matrix (Fin 2) (Fin 2) ℂ :=
  fun i j => if i = 1 ∧ j = 0 then 1 else 0

theorem m2Raise_m2Lower_not_commute :
    ¬ Commute m2Raise m2Lower := by
  intro h
  have h00 := congr_fun (congr_fun h.eq (0 : Fin 2)) (0 : Fin 2)
  simp [m2Raise, m2Lower, Matrix.mul_apply] at h00

/-- **Negative control.** Two distinct labels, both assigned the full `M₂`
algebra and declared disjoint, violate algebraic locality.  Thus finiteness,
region labels, and isotony cannot manufacture the E1 locality receipt. -/
theorem fullM2_distinct_regions_not_local :
    ¬ IsAlgebraicallyLocal (fun U V : Bool => U ≠ V)
      (fun _ => (⊤ : StarSubalgebra ℂ (Matrix (Fin 2) (Fin 2) ℂ))) := by
  intro hlocal
  have hc := hlocal Bool.false_ne_true
      (⟨m2Raise, Set.mem_univ _⟩ :
        (⊤ : StarSubalgebra ℂ (Matrix (Fin 2) (Fin 2) ℂ)))
      (⟨m2Lower, Set.mem_univ _⟩ :
        (⊤ : StarSubalgebra ℂ (Matrix (Fin 2) (Fin 2) ℂ)))
  exact m2Raise_m2Lower_not_commute hc

/-- **Negative control.** Idempotence alone does not imply local
non-disturbance: the zero map is idempotent but does not fix the remote unit
observable. -/
theorem idempotence_does_not_force_remote_fix :
    ∃ E : Matrix (Fin 1) (Fin 1) ℂ →ₗ[ℂ]
        Matrix (Fin 1) (Fin 1) ℂ,
      E.comp E = E ∧ E 1 ≠ 1 := by
  refine ⟨0, ?_, ?_⟩
  · ext X
    simp
  · simp

/-! ## A finite nonvacuity model from the A3/B2 stack -/

section PartitionPublicModel

variable {n k : ℕ}

/-- A finite region-indexed model on the constant A3 tower of a projective
partition.  Regions are subsets of a two-site label set.  Every region is
assigned the same commutative public record algebra, and repair is the B2
partition average.  Given the supplied partition and state, this intentionally
degenerate assignment gives a parameterized consistency construction for the
E1 interface; it is neither a closed witness nor a source-derived spatial
factorization. -/
noncomputable def partitionPublicCausalNet
    (part : EventAlgebra.ProjectivePartition n k)
    (rho : EventAlgebra.StateMatrix n) :
    FiniteCausalObserverNet
      (ConsensusTower.constantConsensusTower part rho) where
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
  localAlgebra := fun _ _ => part.publicSubalgebra
  isotony := fun _ _ _ _ X hX => hX
  locality := by
    intro r U V hUV X Y
    exact congrArg Subtype.val (part.publicSubalgebra_mul_comm X Y)
  restrict := fun _ {_ _} _ => StarAlgHom.id ℂ part.publicSubalgebra
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
  repair := fun _ _ => EventAlgebra.partitionAverageLinearMap part
  repair_idempotent := by
    intro r U X
    exact EventAlgebra.partitionAverage_idem part X
  repair_fixes_region := by
    intro r U X
    exact EventAlgebra.partitionAverage_fixes part X.2
  repair_fixes_disjoint := by
    intro r U V hUV X
    exact EventAlgebra.partitionAverage_fixes part X.2
  repair_natural := by intros; rfl

/-- The nonvacuity model contains a genuinely declared pair of distinct,
nonempty disjoint region labels.  Their assigned algebras are nevertheless
the same commutative algebra, so this is not a tensor-factor construction. -/
theorem partitionPublicCausalNet_has_disjoint_pair
    (part : EventAlgebra.ProjectivePartition n k)
    (rho : EventAlgebra.StateMatrix n) :
    (partitionPublicCausalNet part rho).disjoint ()
      ({0} : Finset (Fin 2)) ({1} : Finset (Fin 2)) := by
  simp [partitionPublicCausalNet]

end PartitionPublicModel

end FiniteCausalObserverNet

-- Axiom audit: no project axiom or admission is used.
#print axioms OPH.QFT.FiniteCausalObserverNet.commute_of_disjoint
#print axioms OPH.QFT.FiniteCausalObserverNet.regionalExpectation_refine
#print axioms OPH.QFT.FiniteCausalObserverNet.relaxedRepair_compose
#print axioms OPH.QFT.FiniteCausalObserverNet.relaxedRepair_remote_expectation
#print axioms OPH.QFT.FiniteCausalObserverNet.TensorSplitReceipt.kraus_remote_marginal_invariant
#print axioms OPH.QFT.FiniteCausalObserverNet.fullM2_distinct_regions_not_local
#print axioms OPH.QFT.FiniteCausalObserverNet.idempotence_does_not_force_remote_fix
#print axioms OPH.QFT.FiniteCausalObserverNet.partitionPublicCausalNet_has_disjoint_pair

end OPH.QFT
