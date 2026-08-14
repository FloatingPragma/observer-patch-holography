import QFT.TowerAnchoredDiamond
import QFT.SourceCorrelationCapstone

/-!
# The structural-net adequacy surface

One display surface for the committed structural-field-theory receipts
of the quantum lane (observation ledger row OL-C6, issue #730).  A
`StructuralNetDeclaredPremises` record bundles the three declared
premises of the committed regional construction: the Cartesian joint
carrier (the joint index of an observer pair is the Cartesian product
of the two label alphabets), the region map (the four-region diamond
assignment of one star subalgebra of the joint matrix algebra to each
region), and the slot split (the two observers enter the joint algebra
as the two tensor factors through the left and right slot embeddings).
Each field is a declared premise of the V3 program: a choice the
source theorems do not produce.  No theorem below derives any of the
three, and each stays declared here.

Two committed instances carry the surface side by side, on different
joint carriers, and their receipts do not compose into one net on one
carrier:

* the designated pair (86, 88) on `Fin 13 × Fin 14` carries net
  existence with the genuine joint coverage law, the full
  conditional-expectation diamond, and the constant-tower transport:
  the diamond transports along the ambient basis equivalence onto one
  constant consensus-tower stage at dimension 182, with the anchoring
  receipts carried;
* the support-disjoint pair (86, 247) on `Fin 13 × Fin 13` carries the
  source-grounded receipts: kernel-checked support disjointness behind
  the split, the factor identification of both observers'
  source-generated algebras with the declared regions, net existence
  with coverage, the source-counted correlation state with exact
  marginal consistency, and the on-carrier erasure of the correlation
  content by the regional expectation pair.

Every theorem in this module is a re-export or a direct conjunction of
committed theorems from `QFT.CPRestrictionNet`,
`QFT.TwoSlotCPNetWitness`, `QFT.TowerAnchoredDiamond`,
`QFT.JointSlotFactorisation`, and `QFT.SourceCorrelationCapstone`,
presented as composition, with one exception: the slot-swap section
proves the finite Z/2 swap receipts of the support-disjoint carrier
(`slotSwap`, `supportDisjoint_swap_covariance`).  The theorem
`supportDisjoint_composed_adequacy` states the support-disjoint
receipts as one conjunction consuming `supportDisjointPairPremises`,
with every clause about the one net `supportDisjointPairNet` and the
one state `correlationState`; the conjunction is presentation-level
composition of the committed receipts, and the clause list is the
derivation record.  Where the committed receipts live on different
carriers they are re-exported side by side rather than composed, and
that separation is part of the committed level.

**What OL-C6 lacks at the committed level.**  The time-slice property
is not attained: the committed region systems carry no time parameter,
no dynamics on the net, and no statement that the algebra of a slice
region determines the top algebra; the coverage law is an algebraic
generation statement over a declared generating family and is not a
time-slice statement.  The exact block: no registered premise supplies
net-compatible dynamics, PR-43's flow lives on one full matrix block
with no region compatibility, and PR-52 owns the physical time and
causal reading.  Genuine local covariance beyond the committed receipts
is not attained: the receipts in this direction are the transport of
one diamond along one basis equivalence onto one constant tower stage,
and the Z/2 slot-swap action of `supportDisjoint_swap_covariance`, an
internal algebra symmetry of the declared two-slot split with
equivariant regional expectations.  The swap is not a spacetime
symmetry; no naturality across tower stages and no nonconstant stage
is attached.  Every statement is finite-dimensional matrix algebra; no
continuum limit, no vacuum state, and no local quantum field is
constructed, and no continuum quantum field theory is claimed.

**Boundary.**  The slot assembly and the label conventions stay
declared postprocessors over post-hoc payload transcriptions,
ineligible as validation; the support disjointness, joint labels,
counts, and marginals of the support-disjoint pair are committed
source data.  No physical channel, instrument, clock, causal relation,
or prediction is attached.
-/

namespace OPH.QFT

open Matrix
open Kronecker
open OPH.Tower
open scoped ComplexOrder

noncomputable section

/-! ## The declared premises, bundled -/

/-- The three declared premises of the committed regional construction,
bundled as one record.  Each field is a declared premise of the V3
program; the source theorems produce none of them. -/
structure StructuralNetDeclaredPremises where
  /-- Declared premise, Cartesian joint carrier (left alphabet): the
  label alphabet of the left observer.  The joint carrier is the
  Cartesian product of this type with `RightLabels`; the product form
  itself is the declared choice. -/
  LeftLabels : Type
  /-- Declared premise, Cartesian joint carrier (right alphabet): the
  label alphabet of the right observer. -/
  RightLabels : Type
  [leftLabelsFintype : Fintype LeftLabels]
  [leftLabelsDecEq : DecidableEq LeftLabels]
  [rightLabelsFintype : Fintype RightLabels]
  [rightLabelsDecEq : DecidableEq RightLabels]
  /-- Declared premise, region map: the assignment of one star
  subalgebra of the joint matrix algebra to each region of the
  two-slot diamond. -/
  regionMap : TwoSlotRegion →
    StarSubalgebra ℂ (Matrix (LeftLabels × RightLabels)
      (LeftLabels × RightLabels) ℂ)
  /-- Declared premise, slot split (left factor): the left observer
  enters the joint algebra through this star embedding. -/
  splitLeft : Matrix LeftLabels LeftLabels ℂ →⋆ₐ[ℂ]
    Matrix (LeftLabels × RightLabels) (LeftLabels × RightLabels) ℂ
  /-- Declared premise, slot split (right factor): the right observer
  enters the joint algebra through this star embedding. -/
  splitRight : Matrix RightLabels RightLabels ℂ →⋆ₐ[ℂ]
    Matrix (LeftLabels × RightLabels) (LeftLabels × RightLabels) ℂ

attribute [instance] StructuralNetDeclaredPremises.leftLabelsFintype
attribute [instance] StructuralNetDeclaredPremises.leftLabelsDecEq
attribute [instance] StructuralNetDeclaredPremises.rightLabelsFintype
attribute [instance] StructuralNetDeclaredPremises.rightLabelsDecEq

/-- The declared Cartesian joint carrier of a premise record. -/
abbrev StructuralNetDeclaredPremises.Carrier
    (P : StructuralNetDeclaredPremises) : Type :=
  P.LeftLabels × P.RightLabels

/-- The committed premise instance of the designated pair (86, 88):
carrier `Fin 13 × Fin 14`, the diamond region map of the two-slot
witness, and the two slot embeddings.  The pair's supports overlap, so
this instance carries no support-disjointness receipt. -/
def designatedPairPremises : StructuralNetDeclaredPremises where
  LeftLabels := Fin 13
  RightLabels := Fin 14
  leftLabelsFintype := inferInstance
  leftLabelsDecEq := inferInstance
  rightLabelsFintype := inferInstance
  rightLabelsDecEq := inferInstance
  regionMap := twoSlotAlgebra
  splitLeft := slotLeft (Fin 14)
  splitRight := slotRight (Fin 13)

/-- The committed premise instance of the support-disjoint pair
(86, 247): carrier `Fin 13 × Fin 13`, the diamond region map of the
justified-pair witness, and the two slot embeddings.  The split of
this pair is backed by the kernel-checked support disjointness
re-exported in `supportDisjoint_split_grounding`; the split itself
stays a declared premise. -/
def supportDisjointPairPremises : StructuralNetDeclaredPremises where
  LeftLabels := Fin 13
  RightLabels := Fin 13
  leftLabelsFintype := inferInstance
  leftLabelsDecEq := inferInstance
  rightLabelsFintype := inferInstance
  rightLabelsDecEq := inferInstance
  regionMap := twoSlotAlgebra247
  splitLeft := slotLeft (Fin 13)
  splitRight := slotRight (Fin 13)

/-! ## Net existence with coverage -/

/-- **Net existence, designated pair.**  The committed two-slot diamond
of the designated pair, typed at the declared carrier.  This
definition is the existence witness; it is the committed `twoSlotNet`
re-exported through the premise record. -/
def designatedPairNet : CPRegionalNet designatedPairPremises.Carrier :=
  twoSlotNet

/-- **Net existence, support-disjoint pair.**  The committed
justified-pair diamond, typed at the declared carrier.  This
definition is the existence witness; it is the committed
`twoSlotNet247` re-exported through the premise record. -/
def supportDisjointPairNet :
    CPRegionalNet supportDisjointPairPremises.Carrier :=
  twoSlotNet247

/-- **Net existence with coverage, designated pair.**  The region type
of the committed net is the diamond, its regional algebras are exactly
the declared region map, and the declared generating pair of slot
regions generates the top algebra (the genuine joint coverage law).
All three components are re-exports of interface fields of
`twoSlotNet`. -/
theorem designatedPair_net_exists_with_coverage :
    (designatedPairNet.Region = TwoSlotRegion) ∧
    (∀ U : TwoSlotRegion,
      designatedPairNet.localAlgebra U =
        designatedPairPremises.regionMap U) ∧
    (StarAlgebra.adjoin ℂ
        (⋃ U ∈ ({TwoSlotRegion.left, TwoSlotRegion.right} :
            Finset TwoSlotRegion),
          (designatedPairPremises.regionMap U :
            Set (Matrix designatedPairPremises.Carrier
              designatedPairPremises.Carrier ℂ))) = ⊤) :=
  ⟨rfl, fun _ => rfl, designatedPairNet.coverage⟩

/-- **Net existence with coverage, support-disjoint pair.**  The same
three receipts for the justified pair, re-exported from
`twoSlotNet247`. -/
theorem supportDisjoint_net_exists_with_coverage :
    (supportDisjointPairNet.Region = TwoSlotRegion) ∧
    (∀ U : TwoSlotRegion,
      supportDisjointPairNet.localAlgebra U =
        supportDisjointPairPremises.regionMap U) ∧
    (StarAlgebra.adjoin ℂ
        (⋃ U ∈ ({TwoSlotRegion.left, TwoSlotRegion.right} :
            Finset TwoSlotRegion),
          (supportDisjointPairPremises.regionMap U :
            Set (Matrix supportDisjointPairPremises.Carrier
              supportDisjointPairPremises.Carrier ℂ))) = ⊤) :=
  ⟨rfl, fun _ => rfl, supportDisjointPairNet.coverage⟩

/-! ## The conditional-expectation diamond -/

/-- **The conditional-expectation diamond, designated pair.**  The
committed restriction data of the net: one conditional expectation per
region, with membership in the declared region map, pointwise fixing
of the region, positivity, trace preservation, the tower law along
region inclusion, unitality, and idempotence.  Every component is an
interface field or interface theorem of `twoSlotNet`, re-exported
through the premise record. -/
theorem designatedPair_conditional_expectation_diamond :
    (∀ (U : TwoSlotRegion)
        (X : Matrix designatedPairPremises.Carrier
          designatedPairPremises.Carrier ℂ),
      designatedPairNet.expect U X ∈ designatedPairPremises.regionMap U) ∧
    (∀ U : TwoSlotRegion,
      ∀ X ∈ designatedPairPremises.regionMap U,
        designatedPairNet.expect U X = X) ∧
    (∀ (U : TwoSlotRegion)
        (X : Matrix designatedPairPremises.Carrier
          designatedPairPremises.Carrier ℂ),
      X.PosSemidef → (designatedPairNet.expect U X).PosSemidef) ∧
    (∀ (U : TwoSlotRegion)
        (X : Matrix designatedPairPremises.Carrier
          designatedPairPremises.Carrier ℂ),
      (designatedPairNet.expect U X).trace = X.trace) ∧
    (∀ U V : TwoSlotRegion, designatedPairNet.regionLE U V →
      ∀ X, designatedPairNet.expect U (designatedPairNet.expect V X) =
        designatedPairNet.expect U X) ∧
    (∀ U : TwoSlotRegion, designatedPairNet.expect U 1 = 1) ∧
    (∀ (U : TwoSlotRegion)
        (X : Matrix designatedPairPremises.Carrier
          designatedPairPremises.Carrier ℂ),
      designatedPairNet.expect U (designatedPairNet.expect U X) =
        designatedPairNet.expect U X) :=
  ⟨designatedPairNet.expect_mem, designatedPairNet.expect_fixes,
    fun U _ hX => designatedPairNet.expect_posSemidef U hX,
    designatedPairNet.expect_trace,
    fun _ _ hUV X => designatedPairNet.expect_tower hUV X,
    fun U => CPRegionalNet.expect_one (N := designatedPairNet) U,
    fun U X => CPRegionalNet.expect_idem (N := designatedPairNet) U X⟩

/-- **The conditional-expectation diamond, support-disjoint pair.**
Five of the designated-pair receipts for the justified pair
(membership, fixing, positivity, trace, tower), re-exported from
`twoSlotNet247`; unitality and idempotence hold generically for every
`CPRegionalNet` and are not restated here. -/
theorem supportDisjoint_conditional_expectation_diamond :
    (∀ (U : TwoSlotRegion)
        (X : Matrix supportDisjointPairPremises.Carrier
          supportDisjointPairPremises.Carrier ℂ),
      supportDisjointPairNet.expect U X ∈
        supportDisjointPairPremises.regionMap U) ∧
    (∀ U : TwoSlotRegion,
      ∀ X ∈ supportDisjointPairPremises.regionMap U,
        supportDisjointPairNet.expect U X = X) ∧
    (∀ (U : TwoSlotRegion)
        (X : Matrix supportDisjointPairPremises.Carrier
          supportDisjointPairPremises.Carrier ℂ),
      X.PosSemidef → (supportDisjointPairNet.expect U X).PosSemidef) ∧
    (∀ (U : TwoSlotRegion)
        (X : Matrix supportDisjointPairPremises.Carrier
          supportDisjointPairPremises.Carrier ℂ),
      (supportDisjointPairNet.expect U X).trace = X.trace) ∧
    (∀ U V : TwoSlotRegion, supportDisjointPairNet.regionLE U V →
      ∀ X, supportDisjointPairNet.expect U
          (supportDisjointPairNet.expect V X) =
        supportDisjointPairNet.expect U X) :=
  ⟨supportDisjointPairNet.expect_mem, supportDisjointPairNet.expect_fixes,
    fun U _ hX => supportDisjointPairNet.expect_posSemidef U hX,
    supportDisjointPairNet.expect_trace,
    fun _ _ hUV X => supportDisjointPairNet.expect_tower hUV X⟩

/-! ## Constant-tower transport -/

/-- **The constant-tower transport receipt.**  The committed diamond of the designated pair transports
along the ambient basis star-equivalence onto the private algebra of
one constant consensus-tower stage at dimension 182, and the anchoring
receipts hold: the anchored net is definitionally the transported
diamond, the tower stage's private algebra is the anchor ambient,
every member of the tower's public checkpoint partition lies in the
anchored left region, and the anchored left region is a proper
subalgebra.  This transport along one equivalence onto one constant
stage is the whole committed content in this direction; no
symmetry-group action on regions and no nonconstant stage is
attached. -/
theorem designatedPair_constant_tower_transport :
    (anchoredNet = twoSlotNet.transport ambientEquiv
        (fun hX => posSemidef_reindex_equiv pairFin hX)
        (fun hY => posSemidef_reindex_equiv pairFin.symm hY)
        (fun X => trace_reindex_equiv pairFin X)) ∧
    (ConsensusTower.PrivateAlgebra anchoredTower () =
        Matrix (Fin 182) (Fin 182) ℂ) ∧
    (∀ v : Fin 4, anchoredCheckpointPartition.proj v ∈
        anchoredNet.localAlgebra
          (show anchoredNet.Region from TwoSlotRegion.left)) ∧
    (anchoredNet.localAlgebra
        (show anchoredNet.Region from TwoSlotRegion.left) ≠ ⊤) :=
  ⟨rfl, anchoredTower_privateAlgebra, partition_members_in_left_region,
    anchoredNet_left_ne_top⟩

/-! ## The source-grounded split of the support-disjoint pair -/

/-- **Split grounding.**  The kernel-checked source disjointness of the
two full supports, and the identification of both observers'
source-generated algebras with the declared regions through the
declared split embeddings.  The disjointness is a committed source
datum; the split through which it is read stays declared. -/
theorem supportDisjoint_split_grounding :
    (∀ i j : Fin 96, support86Full i ≠ support247Full j) ∧
    (StarSubalgebra.map supportDisjointPairPremises.splitLeft
        obs86.sourceAlgebra =
      supportDisjointPairPremises.regionMap TwoSlotRegion.left) ∧
    (StarSubalgebra.map supportDisjointPairPremises.splitRight
        obs247.sourceAlgebra =
      supportDisjointPairPremises.regionMap TwoSlotRegion.right) :=
  ⟨supports_disjoint, obs86_lifted_eq_leftSlot247,
    obs247_lifted_eq_rightSlot247⟩

/-- Factor identification of the designated pair through its declared
split embeddings.  This pair carries no support-disjointness receipt;
its supports overlap. -/
theorem designatedPair_split_identification :
    (StarSubalgebra.map designatedPairPremises.splitLeft
        obs86.sourceAlgebra =
      designatedPairPremises.regionMap TwoSlotRegion.left) ∧
    (StarSubalgebra.map designatedPairPremises.splitRight
        obs88.sourceAlgebra =
      designatedPairPremises.regionMap TwoSlotRegion.right) :=
  ⟨obs86_lifted_eq_leftSlot, obs88_lifted_eq_rightSlot⟩

/-! ## The correlation-carrying counted state -/

/-- **The correlation-carrying counted state with marginal
consistency.**  The source-counted joint empirical state of the
support-disjoint pair: total count 32, positive semidefinite, trace
one, both partial traces equal to the marginal empirical states
exactly, and different from the product of its marginals, so the
counted state carries genuine correlation content.  All components are
re-exports from the source correlation capstone. -/
theorem supportDisjoint_correlation_counted_state :
    ((∑ p : Fin 13 × Fin 13, jointCount p) = 32) ∧
    correlationState.PosSemidef ∧
    correlationState.trace = 1 ∧
    ptraceSnd correlationState = marginal86State ∧
    OPH.Locality.ptraceFst correlationState = marginal247State ∧
    correlationState ≠ marginal86State ⊗ₖ marginal247State :=
  ⟨jointCount_total, correlationState_posSemidef, correlationState_trace,
    correlationState_marginal_left, correlationState_marginal_right,
    correlationState_ne_product⟩

/-- The regional conditional expectations of the committed net
reproduce the two marginal empirical states of the counted joint state
exactly, phrased through the net's expectation fields. -/
theorem supportDisjointNet_expect_recovers_marginals :
    (supportDisjointPairNet.expect TwoSlotRegion.left correlationState =
        (13 : ℂ)⁻¹ •
          (marginal86State ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ))) ∧
    (supportDisjointPairNet.expect TwoSlotRegion.right correlationState =
        (13 : ℂ)⁻¹ •
          ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ marginal247State)) :=
  ⟨leftSlotExpectation_correlationState,
    rightSlotExpectation_correlationState⟩

/-! ## The on-carrier erasure -/

/-- **On-carrier erasure.**  The two regional expectations of the
committed net map the counted joint state and the product of its
marginals to the same expectation pair while the two states differ:
the expectation pair is not injective at exactly the source
correlation, so the joint state is the datum the net must carry.
Re-export of the capstone erasure theorem, phrased through the net's
expectation fields. -/
theorem supportDisjoint_on_carrier_erasure :
    ((supportDisjointPairNet.expect TwoSlotRegion.left correlationState,
        supportDisjointPairNet.expect TwoSlotRegion.right
          correlationState) =
      (supportDisjointPairNet.expect TwoSlotRegion.left
          (marginal86State ⊗ₖ marginal247State),
        supportDisjointPairNet.expect TwoSlotRegion.right
          (marginal86State ⊗ₖ marginal247State))) ∧
    correlationState ≠ marginal86State ⊗ₖ marginal247State :=
  slotExpectations_erase_source_correlation

/-! ## The support-disjoint receipts as one composed conjunction -/

/-- **The composed support-disjoint adequacy statement.**  One
conjunction consuming the premise bundle `supportDisjointPairPremises`
(register rows PR-32, PR-33, PR-34), with every clause about the one
net `supportDisjointPairNet` and the one state `correlationState`: net
existence with the joint coverage law, the conditional-expectation
diamond, the source-grounded split, the correlation-carrying counted
state with exact marginal consistency, the marginal recovery through
the net's expectations, and the on-carrier erasure.  Every clause is a
committed receipt of this module; the conjunction is the
presentation-level composition that lets the ledger row cite one
statement on the shared objects, and it adds no derivation. -/
theorem supportDisjoint_composed_adequacy :
    ((supportDisjointPairNet.Region = TwoSlotRegion) ∧
      (∀ U : TwoSlotRegion,
        supportDisjointPairNet.localAlgebra U =
          supportDisjointPairPremises.regionMap U) ∧
      (StarAlgebra.adjoin ℂ
          (⋃ U ∈ ({TwoSlotRegion.left, TwoSlotRegion.right} :
              Finset TwoSlotRegion),
            (supportDisjointPairPremises.regionMap U :
              Set (Matrix supportDisjointPairPremises.Carrier
                supportDisjointPairPremises.Carrier ℂ))) = ⊤)) ∧
    ((∀ (U : TwoSlotRegion)
        (X : Matrix supportDisjointPairPremises.Carrier
          supportDisjointPairPremises.Carrier ℂ),
      supportDisjointPairNet.expect U X ∈
        supportDisjointPairPremises.regionMap U) ∧
      (∀ U : TwoSlotRegion,
        ∀ X ∈ supportDisjointPairPremises.regionMap U,
          supportDisjointPairNet.expect U X = X) ∧
      (∀ (U : TwoSlotRegion)
          (X : Matrix supportDisjointPairPremises.Carrier
            supportDisjointPairPremises.Carrier ℂ),
        X.PosSemidef → (supportDisjointPairNet.expect U X).PosSemidef) ∧
      (∀ (U : TwoSlotRegion)
          (X : Matrix supportDisjointPairPremises.Carrier
            supportDisjointPairPremises.Carrier ℂ),
        (supportDisjointPairNet.expect U X).trace = X.trace) ∧
      (∀ U V : TwoSlotRegion, supportDisjointPairNet.regionLE U V →
        ∀ X, supportDisjointPairNet.expect U
            (supportDisjointPairNet.expect V X) =
          supportDisjointPairNet.expect U X)) ∧
    ((∀ i j : Fin 96, support86Full i ≠ support247Full j) ∧
      (StarSubalgebra.map supportDisjointPairPremises.splitLeft
          obs86.sourceAlgebra =
        supportDisjointPairPremises.regionMap TwoSlotRegion.left) ∧
      (StarSubalgebra.map supportDisjointPairPremises.splitRight
          obs247.sourceAlgebra =
        supportDisjointPairPremises.regionMap TwoSlotRegion.right)) ∧
    (((∑ p : Fin 13 × Fin 13, jointCount p) = 32) ∧
      correlationState.PosSemidef ∧
      correlationState.trace = 1 ∧
      ptraceSnd correlationState = marginal86State ∧
      OPH.Locality.ptraceFst correlationState = marginal247State ∧
      correlationState ≠ marginal86State ⊗ₖ marginal247State) ∧
    ((supportDisjointPairNet.expect TwoSlotRegion.left correlationState =
        (13 : ℂ)⁻¹ •
          (marginal86State ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ))) ∧
      (supportDisjointPairNet.expect TwoSlotRegion.right correlationState =
        (13 : ℂ)⁻¹ •
          ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ marginal247State))) ∧
    (((supportDisjointPairNet.expect TwoSlotRegion.left correlationState,
        supportDisjointPairNet.expect TwoSlotRegion.right
          correlationState) =
      (supportDisjointPairNet.expect TwoSlotRegion.left
          (marginal86State ⊗ₖ marginal247State),
        supportDisjointPairNet.expect TwoSlotRegion.right
          (marginal86State ⊗ₖ marginal247State))) ∧
      correlationState ≠ marginal86State ⊗ₖ marginal247State) :=
  ⟨supportDisjoint_net_exists_with_coverage,
    supportDisjoint_conditional_expectation_diamond,
    supportDisjoint_split_grounding,
    supportDisjoint_correlation_counted_state,
    supportDisjointNet_expect_recovers_marginals,
    supportDisjoint_on_carrier_erasure⟩

/-! ## The Z/2 slot-swap action on the support-disjoint diamond -/

/-- The label-swap map on the two-slot diamond: exchanges the left and
right regions and fixes bottom and top. -/
def TwoSlotRegion.swap : TwoSlotRegion → TwoSlotRegion
  | TwoSlotRegion.bot => TwoSlotRegion.bot
  | TwoSlotRegion.left => TwoSlotRegion.right
  | TwoSlotRegion.right => TwoSlotRegion.left
  | TwoSlotRegion.top => TwoSlotRegion.top

/-- The region swap is an involution. -/
theorem TwoSlotRegion.swap_swap : ∀ U : TwoSlotRegion, U.swap.swap = U := by
  decide

/-- The region swap preserves the diamond order. -/
theorem TwoSlotRegion.swap_monotone : ∀ U V : TwoSlotRegion,
    TwoSlotRegion.le U V = true →
      TwoSlotRegion.le U.swap V.swap = true := by
  decide

/-- The slot-swap star-algebra equivalence of the support-disjoint
carrier: matrix reindexing along the exchange of the two tensor
factors.  This is a constructed object, not a declared premise; it is
an internal algebra symmetry of the declared two-slot split, not a
spacetime symmetry. -/
def slotSwap :
    Matrix PairIndex247 PairIndex247 ℂ ≃⋆ₐ[ℂ]
      Matrix PairIndex247 PairIndex247 ℂ :=
  StarAlgEquiv.ofAlgEquiv
    (Matrix.reindexAlgEquiv ℂ ℂ (Equiv.prodComm (Fin 13) (Fin 13)))
    (fun M => by
      simp only [Matrix.reindexAlgEquiv_apply, Matrix.star_eq_conjTranspose]
      exact reindex_conjTranspose _ M)

theorem slotSwap_apply (M : Matrix PairIndex247 PairIndex247 ℂ)
    (z w : PairIndex247) :
    slotSwap M z w = M (z.2, z.1) (w.2, w.1) :=
  rfl

/-- The slot swap exchanges the two Kronecker factors. -/
theorem slotSwap_kronecker (A B : Matrix (Fin 13) (Fin 13) ℂ) :
    slotSwap (A ⊗ₖ B) = B ⊗ₖ A := by
  ext z w
  rw [slotSwap_apply]
  simp only [Matrix.kroneckerMap_apply]
  exact mul_comm _ _

/-- The slot swap is an involution. -/
theorem slotSwap_involutive (M : Matrix PairIndex247 PairIndex247 ℂ) :
    slotSwap (slotSwap M) = M :=
  rfl

/-- The swap exchanges the two partial traces (first form). -/
theorem ptraceFst_slotSwap (M : Matrix PairIndex247 PairIndex247 ℂ) :
    OPH.Locality.ptraceFst (slotSwap M) = ptraceSnd M :=
  rfl

/-- The swap exchanges the two partial traces (second form). -/
theorem ptraceSnd_slotSwap (M : Matrix PairIndex247 PairIndex247 ℂ) :
    ptraceSnd (slotSwap M) = OPH.Locality.ptraceFst M :=
  rfl

/-- **Region covariance of the swap.**  The slot swap maps each
regional algebra of the support-disjoint diamond onto the regional
algebra of the swapped region. -/
theorem slotSwap_map_localAlgebra (U : TwoSlotRegion) :
    StarSubalgebra.map
        (slotSwap : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
          Matrix PairIndex247 PairIndex247 ℂ)
        (twoSlotAlgebra247 U) =
      twoSlotAlgebra247 U.swap := by
  cases U
  · -- bottom: the scalars map onto the scalars.
    apply le_antisymm
    · rintro x ⟨y, hy, rfl⟩
      obtain ⟨c, rfl⟩ := StarSubalgebra.mem_bot.mp hy
      exact ⟨c, (AlgHomClass.commutes
        (slotSwap : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
          Matrix PairIndex247 PairIndex247 ℂ) c).symm⟩
    · intro x hx
      obtain ⟨c, rfl⟩ := StarSubalgebra.mem_bot.mp hx
      exact ⟨algebraMap ℂ _ c, ⟨c, rfl⟩,
        AlgHomClass.commutes
          (slotSwap : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
            Matrix PairIndex247 PairIndex247 ℂ) c⟩
  · -- left: the left factor maps onto the right factor.
    apply le_antisymm
    · rintro x ⟨y, ⟨A, rfl⟩, rfl⟩
      refine ⟨A, ?_⟩
      show (1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ A =
        slotSwap (A ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ))
      exact (slotSwap_kronecker A 1).symm
    · rintro x ⟨B, rfl⟩
      refine ⟨slotLeft (Fin 13) B, ⟨B, rfl⟩, ?_⟩
      show slotSwap (B ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ)) =
        (1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ B
      exact slotSwap_kronecker B 1
  · -- right: the right factor maps onto the left factor.
    apply le_antisymm
    · rintro x ⟨y, ⟨B, rfl⟩, rfl⟩
      refine ⟨B, ?_⟩
      show B ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ) =
        slotSwap ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ B)
      exact (slotSwap_kronecker 1 B).symm
    · rintro x ⟨A, rfl⟩
      refine ⟨slotRight (Fin 13) A, ⟨A, rfl⟩, ?_⟩
      show slotSwap ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ A) =
        A ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ)
      exact slotSwap_kronecker 1 A
  · -- top: a star-algebra equivalence is surjective.
    show StarSubalgebra.map
        (slotSwap : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
          Matrix PairIndex247 PairIndex247 ℂ) ⊤ = ⊤
    rw [← StarAlgHom.range_eq_map_top, StarSubalgebra.eq_top_iff]
    intro x
    exact ⟨slotSwap.symm x, slotSwap.apply_symm_apply x⟩

/-- **Expectation equivariance of the swap.**  The regional conditional
expectations of the support-disjoint net intertwine the slot swap: the
expectation of the swapped region evaluated on the swapped matrix is
the swap of the original expectation. -/
theorem slotSwap_expect (U : TwoSlotRegion)
    (M : Matrix PairIndex247 PairIndex247 ℂ) :
    supportDisjointPairNet.expect U.swap (slotSwap M) =
      slotSwap (supportDisjointPairNet.expect U M) := by
  have hswap_smul_kron : ∀ (c : ℂ) (A B : Matrix (Fin 13) (Fin 13) ℂ),
      slotSwap (c • (A ⊗ₖ B)) = c • (B ⊗ₖ A) := by
    intro c A B
    ext z w
    rw [slotSwap_apply]
    simp only [Matrix.smul_apply, Matrix.kroneckerMap_apply, smul_eq_mul]
    ring
  cases U
  · -- bottom onto bottom: the scalar expectation.
    show scalarExpectation PairIndex247 (slotSwap M) =
      slotSwap (scalarExpectation PairIndex247 M)
    have htr : (slotSwap M).trace = M.trace :=
      trace_reindex_equiv (Equiv.prodComm (Fin 13) (Fin 13)) M
    have h1 : scalarExpectation PairIndex247 (slotSwap M) =
        (Fintype.card PairIndex247 : ℂ)⁻¹ • (slotSwap M).trace •
          (1 : Matrix PairIndex247 PairIndex247 ℂ) := rfl
    have h2 : scalarExpectation PairIndex247 M =
        (Fintype.card PairIndex247 : ℂ)⁻¹ • M.trace •
          (1 : Matrix PairIndex247 PairIndex247 ℂ) := rfl
    rw [h1, h2, htr]
    ext z w
    rw [slotSwap_apply]
    simp only [Matrix.smul_apply, smul_eq_mul, Matrix.one_apply]
    by_cases h : z = w
    · subst h
      simp
    · have h' : ¬ ((z.2, z.1) = (w.2, w.1)) := by
        intro hc
        exact h (Prod.ext (congrArg Prod.snd hc) (congrArg Prod.fst hc))
      simp [h, h']
  · -- left onto right.
    show rightSlotExpectation (slotSwap M) =
      slotSwap (leftSlotExpectation M)
    have h1 : rightSlotExpectation (α := Fin 13) (β := Fin 13)
        (slotSwap M) =
        (Fintype.card (Fin 13) : ℂ)⁻¹ •
          ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
            OPH.Locality.ptraceFst (slotSwap M)) := rfl
    have h2 : leftSlotExpectation (α := Fin 13) (β := Fin 13) M =
        (Fintype.card (Fin 13) : ℂ)⁻¹ •
          (ptraceSnd M ⊗ₖ (1 : Matrix (Fin 13) (Fin 13) ℂ)) := rfl
    rw [h1, h2, ptraceFst_slotSwap, hswap_smul_kron]
  · -- right onto left.
    show leftSlotExpectation (slotSwap M) =
      slotSwap (rightSlotExpectation M)
    have h1 : leftSlotExpectation (α := Fin 13) (β := Fin 13)
        (slotSwap M) =
        (Fintype.card (Fin 13) : ℂ)⁻¹ •
          (ptraceSnd (slotSwap M) ⊗ₖ
            (1 : Matrix (Fin 13) (Fin 13) ℂ)) := rfl
    have h2 : rightSlotExpectation (α := Fin 13) (β := Fin 13) M =
        (Fintype.card (Fin 13) : ℂ)⁻¹ •
          ((1 : Matrix (Fin 13) (Fin 13) ℂ) ⊗ₖ
            OPH.Locality.ptraceFst M) := rfl
    rw [h1, h2, ptraceSnd_slotSwap, hswap_smul_kron]
  · -- top onto top: the identity expectation.
    rfl

/-- **The swap-covariance receipt.**  The Z/2 slot-swap action on the
support-disjoint diamond, in one statement: the swap maps each regional
algebra onto the swapped region's algebra, the regional expectations
intertwine the swap, the swap preserves the region order, and the swap
is an involution.  This is a symmetry-group action on the region system
of the committed net, strictly beyond the constant-tower transport; it
is an internal algebra symmetry of the declared two-slot split and
carries no spacetime reading. -/
theorem supportDisjoint_swap_covariance :
    (∀ U : TwoSlotRegion,
      StarSubalgebra.map
          (slotSwap : Matrix PairIndex247 PairIndex247 ℂ →⋆ₐ[ℂ]
            Matrix PairIndex247 PairIndex247 ℂ)
          (twoSlotAlgebra247 U) =
        twoSlotAlgebra247 U.swap) ∧
    (∀ (U : TwoSlotRegion) (M : Matrix PairIndex247 PairIndex247 ℂ),
      supportDisjointPairNet.expect U.swap (slotSwap M) =
        slotSwap (supportDisjointPairNet.expect U M)) ∧
    (∀ U V : TwoSlotRegion, supportDisjointPairNet.regionLE U V →
      supportDisjointPairNet.regionLE U.swap V.swap) ∧
    (∀ M : Matrix PairIndex247 PairIndex247 ℂ,
      slotSwap (slotSwap M) = M) :=
  ⟨slotSwap_map_localAlgebra, slotSwap_expect,
    fun U V hUV => TwoSlotRegion.swap_monotone U V hUV,
    slotSwap_involutive⟩

end

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.designatedPairPremises
#print axioms OPH.QFT.supportDisjointPairPremises
#print axioms OPH.QFT.designatedPairNet
#print axioms OPH.QFT.supportDisjointPairNet
#print axioms OPH.QFT.designatedPair_net_exists_with_coverage
#print axioms OPH.QFT.supportDisjoint_net_exists_with_coverage
#print axioms OPH.QFT.designatedPair_conditional_expectation_diamond
#print axioms OPH.QFT.supportDisjoint_conditional_expectation_diamond
#print axioms OPH.QFT.designatedPair_constant_tower_transport
#print axioms OPH.QFT.supportDisjoint_split_grounding
#print axioms OPH.QFT.designatedPair_split_identification
#print axioms OPH.QFT.supportDisjoint_correlation_counted_state
#print axioms OPH.QFT.supportDisjointNet_expect_recovers_marginals
#print axioms OPH.QFT.supportDisjoint_on_carrier_erasure
#print axioms OPH.QFT.supportDisjoint_composed_adequacy
#print axioms OPH.QFT.TwoSlotRegion.swap_swap
#print axioms OPH.QFT.TwoSlotRegion.swap_monotone
#print axioms OPH.QFT.slotSwap
#print axioms OPH.QFT.slotSwap_apply
#print axioms OPH.QFT.slotSwap_kronecker
#print axioms OPH.QFT.slotSwap_involutive
#print axioms OPH.QFT.ptraceFst_slotSwap
#print axioms OPH.QFT.ptraceSnd_slotSwap
#print axioms OPH.QFT.slotSwap_map_localAlgebra
#print axioms OPH.QFT.slotSwap_expect
#print axioms OPH.QFT.supportDisjoint_swap_covariance
