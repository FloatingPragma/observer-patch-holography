import TrichotomyCases
import PortDualMetricSelection
import ExteriorSelection
import Z6Descent
import QuantumMatterIntegration
import A5FamilyBand

open scoped BigOperators

/-!
# The composed Standard Model structure surface (V3, issue #734)

One citable package for the observation-ledger rows OL-G1 through
OL-G5: the forced gauge Lie type, the measure-conditional selection of
the compact family `G`, the exterior one-generation matter selection
with anomaly cancellation and its parity grading, the `Z6` kernel with
its conditional global-form status, and the family-band receipts
behind the three-family row.  This module adds no mathematical
content; every theorem re-exports a committed result through the
premise structure below.  Where the committed hypotheses do not
literally compose, the re-exports stand side by side inside the
premise-parameterized namespace and this header names the gap.

`SMStructurePremiseData` bundles the supplied objects the committed
theorems consume, in three blocks.  The port-dual weighting block
carries a vertex share pinned to one third and a constant normalized
face weight on the twenty reconstructed faces: the equal-face-weight
and barycentric one-third rules of register row PR-09.  The
classification block carries the centre dimension and the simple-ideal
dimension multiset with the hypotheses of the committed trichotomy:
the rational-centre dimension list, the compact-simple dimension list
below twelve, the twelve-dimensional total, the centre bound from the
S5 centre receipt, and the exclusion of centre dimension zero (the
declared fixed-space argument against `su(2)^4`).  The Lie-theoretic
bridge behind those lists is the compact-Lie classification import of
register row PR-11.  The matter block carries a selection mask
constrained exactly by the exterior selection boundary recorded under register row
PR-12: nonempty, chiral, and anomaly free over the ten exterior
components.  The measure-to-metric and nearest-point repair rules of
register row PR-10 enter through the committed definitions rather
than through a field: the diagonal reading substitutes the port-dual
weight into all three sector-scale slots of the distance forms, and
the discriminator verdict is the nearest classified compact family.

Two further physical premises are registered as PR-35 and PR-36, but
they are not fields of `SMStructurePremiseData` and are not discharged
by the composite theorem below.  The
same-source loop-to-kernel identity (the physical carrier-loop class
of the source is the six-axis gluing class, and its image under the
committed intertwiner generates the tensor-action kernel) gates the
global-form row OL-G2; the theorems here supply the abstract
intertwiner and the four-quotient menu only, so the row stays
conditional.  The family band-selection rules (the family
multiplicity object is a single complete faithful band of the
twelve-port adjacency with dimension in the pinned window, and
admissible bands are compared by the operational seam-cost order)
gate the three-family row OL-G4; the theorems here supply the
selection arithmetic on the committed carrier, and exactly those
rules plus the realization clause of the executable certificate
`family_band_attachment_certificate.py` separate the row from a
derived family triplication.

Composition gaps.  The measure-conditional selection picks the
compact family `G` in the classified locus, and the trichotomy forces
the dimension data `(1, {3, 8})`; the identification of the selected
family with the forced type passes through the pinned compact-locus
classification of row PR-11, so the two theorems stand side by side
here and no composition across that bridge is claimed.  The committed
gap theorem of `A5OPH` (a `G`-module isomorphism determines no Lie
bracket) scopes the type claim.  The matter selection and the band
selection concern different committed objects (the exterior component
table and the adjacency spectrum) and are re-exported beside one
another without a joint statement.

Boundary.  No physical current, laboratory charge, coupling,
Lagrangian assembly, symmetry breaking, or spacetime Spin attachment
is claimed.  The parity grading is a scan output on the declared
component table; a physical chirality claim consumes a Spin
attachment, and the committed receipt
`Z6Descent.no_universal_fermion_minus_one` is the recorded
obstruction to reading one off the finite central arithmetic.  The
global form is conditional: the realized local weights descend
through all four admissible quotients, so the local table selects no
physical form.
-/

namespace OPH.SMStructureAdequacySurface

open OPH.PrimitivePortDualMeasure (FaceTriple baseFaces portDualWeight)
open OPH.OrientedFaceInvariantMetric (dG2 squaredDistanceAt)
open OPH.OrientedFaceBracketSelector (CompactFamily)

/-- The supplied objects of the composed Standard Model structure
surface.  The first block carries the port-dual weighting rules of
register row PR-09: a vertex share `vertexShare` pinned to the
barycentric one third and a face weight `faceWeight` that is constant
and normalized on the twenty reconstructed faces.  The second block
carries the classification inputs of the committed trichotomy, with
the Lie-theoretic bridge behind them the compact-Lie classification
import of register row PR-11: the centre dimension `centreDim` in the
rational-centre list, the simple-ideal dimension multiset
`simpleIdealDims` drawn from the compact-simple dimension list below
twelve, the twelve-dimensional total, the S5 centre bound, and the
exclusion of centre dimension zero.  The third block carries the
exterior selection boundary recorded under register row PR-12: a
selection mask over the ten exterior components that is nonempty,
chiral, and anomaly free.  The measure-to-metric and nearest-point repair rules of
register row PR-10 are consumed through the committed distance forms,
so they carry no field here. -/
structure SMStructurePremiseData where
  /-- The declared per-vertex share of a face weight. -/
  vertexShare : ℚ
  /-- The declared weight on reconstructed faces. -/
  faceWeight : FaceTriple → ℚ
  /-- The share is barycentric: one third per incident port
  (register row PR-09). -/
  vertexShare_barycentric : vertexShare = 1 / 3
  /-- The face weight is constant across the twenty reconstructed
  faces (register row PR-09). -/
  faceWeight_constant : ∀ f ∈ baseFaces, ∀ g ∈ baseFaces,
    faceWeight f = faceWeight g
  /-- The face weight is normalized (register row PR-09). -/
  faceWeight_normalized : ∑ f ∈ baseFaces, faceWeight f = 1
  /-- The centre dimension of the compact coefficient algebra. -/
  centreDim : ℕ
  /-- The multiset of simple-ideal dimensions. -/
  simpleIdealDims : Multiset ℕ
  /-- The centre dimension lies in the rational-centre subset-sum
  list (the `A5CharacterField` receipt; the torus/cocharacter step is
  part of register row PR-11). -/
  centreDim_mem : centreDim ∈ ({0, 1, 5, 6, 7, 11, 12} : Finset ℕ)
  /-- Every simple-ideal dimension lies in the compact-simple list
  below twelve (classification input of register row PR-11). -/
  simpleIdealDims_mem : ∀ d ∈ simpleIdealDims, d ∈ OPHSelection.SimpleDims
  /-- Centre and semisimple dimensions total the port-module
  dimension twelve. -/
  dims_total : centreDim + simpleIdealDims.sum = 12
  /-- The S5 centre bound of the `A5PortModule` receipt. -/
  centreDim_le_one : centreDim ≤ 1
  /-- The declared fixed-space exclusion of centre dimension zero
  (part of register row PR-11). -/
  centreDim_ne_zero : centreDim ≠ 0
  /-- The matter selection mask over the ten exterior components
  (register row PR-12). -/
  selectionMask : Fin 1024
  /-- The selection is nonempty (register row PR-12). -/
  selection_nonempty : selectionMask.val ≠ 0
  /-- The selection is chiral: no component appears with its
  conjugate (register row PR-12). -/
  selection_chiral :
    OPH.ExteriorSelection.chiral selectionMask.val = true
  /-- All four anomaly forms vanish on the selection (register row
  PR-12). -/
  selection_anomalyFree :
    OPH.ExteriorSelection.anomalyFree selectionMask.val = true

namespace SMStructurePremiseData

variable (D : SMStructurePremiseData)

/-! ## The forced gauge Lie type (OL-G1) -/

/-- The four-case enumeration behind the forced type: under the
bundled classification inputs the centre-versus-semisimple dimension
data is abelian, the dimension-six branch, the Standard-Model type,
or `su(2)^4`.  Re-export of
`OPH.TrichotomyCases.centre_case_enumeration`. -/
theorem centre_semisimple_enumeration :
    (D.centreDim = 12 ∧ D.simpleIdealDims = 0)
      ∨ (D.centreDim = 6 ∧ D.simpleIdealDims = {3, 3})
      ∨ (D.centreDim = 1 ∧ D.simpleIdealDims = {3, 8})
      ∨ (D.centreDim = 0 ∧ D.simpleIdealDims = {3, 3, 3, 3}) :=
  OPH.TrichotomyCases.centre_case_enumeration D.centreDim
    D.simpleIdealDims D.centreDim_mem D.simpleIdealDims_mem D.dims_total

/-- **The forced gauge Lie type (OL-G1).**  Under the bundled
classification inputs the centre dimension is one and the simple-ideal
dimension multiset is `{3, 8}`.  With the compact-simple dimension
list of register row PR-11 this dimension data reads as the type
`u(1) + su(2) + su(3)`; the reading is part of the row PR-11 import,
and this theorem is dimension arithmetic only.  Re-export of
`OPH.TrichotomyCases.inner_action_selection` at the bundled
premises. -/
theorem forced_gauge_type :
    D.centreDim = 1 ∧ D.simpleIdealDims = {3, 8} :=
  OPH.TrichotomyCases.inner_action_selection D.centreDim
    D.simpleIdealDims D.simpleIdealDims_mem D.dims_total
    D.centreDim_le_one D.centreDim_ne_zero

/-! ## The measure-conditional selection of the family G (OL-G1) -/

/-- The bundled weighting rules pin the port-dual weight to `1/12` at
every port.  Re-export of
`OPH.PrimitivePortDualMeasure.equal_face_weights_imply_port_dual_one_twelfth`. -/
theorem port_dual_weight_pinned (p : Fin 12) :
    portDualWeight D.vertexShare D.faceWeight p = 1 / 12 :=
  OPH.PrimitivePortDualMeasure.equal_face_weights_imply_port_dual_one_twelfth
    D.vertexShare D.faceWeight D.vertexShare_barycentric
    D.faceWeight_constant D.faceWeight_normalized p

/-- **The measure-conditional family selection (OL-G1).**  Under the
bundled weighting rules and the declared measure-to-metric and
nearest-point rules of register row PR-10, the compact family `G` is
strictly nearer to the pinned face bracket than every other classified
family, at the metric induced by the port-dual weight of any port.
The compact locus is conditional on the pinned classification
certificate (register row PR-11).  Re-export of
`OPH.PortDualMetricSelection.port_dual_weight_selects_G` at the
bundled premises. -/
theorem measure_selects_G (p : Fin 12) (family : CompactFamily)
    (h : family ≠ .G) :
    dG2 ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)
        ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)
        ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ) <
      squaredDistanceAt
        ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)
        ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)
        ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)
        family :=
  OPH.PortDualMetricSelection.port_dual_weight_selects_G
    D.vertexShare D.faceWeight D.vertexShare_barycentric
    D.faceWeight_constant D.faceWeight_normalized p family h

/-! ## The exterior one-generation selection (OL-G3, OL-G5) -/

/-- **The one-generation selection (OL-G3).**  The bundled selection
mask, constrained exactly by the matter candidate grammar of register
row PR-12, is one of the two conjugate fermionic-parity sectors of the
ten-component exterior table.  The theorem does not choose between
them.  Re-export of `OPH.ExteriorSelection.selection_unique` at the
bundled premises. -/
theorem matter_selection_is_parity_sector :
    D.selectionMask.val = OPH.ExteriorSelection.evenMask
      ∨ D.selectionMask.val = OPH.ExteriorSelection.oddMask :=
  OPH.ExteriorSelection.selection_unique D.selectionMask
    D.selection_nonempty D.selection_chiral D.selection_anomalyFree

/-- **The parity grading of the selection (OL-G5, attained part).**
Componentwise, the bundled selection agrees with the even-parity table
or with the odd-parity table of the declared component grammar.  The
grading is table provenance: the parity column is part of register row
PR-12, and the scan outputs which sectors survive.  A physical
chirality claim consumes a spacetime Spin attachment that no theorem
here supplies. -/
theorem selection_parity_graded :
    (∀ i : Fin 10, OPH.ExteriorSelection.mem D.selectionMask.val i
        = !OPH.ExteriorSelection.odd i)
      ∨ (∀ i : Fin 10, OPH.ExteriorSelection.mem D.selectionMask.val i
        = OPH.ExteriorSelection.odd i) := by
  rcases D.matter_selection_is_parity_sector with h | h
  · refine Or.inl fun i => ?_
    rw [h]
    exact OPH.ExteriorSelection.evenMask_is_even_sector i
  · refine Or.inr fun i => ?_
    rw [h]
    exact OPH.ExteriorSelection.oddMask_is_odd_sector i

end SMStructurePremiseData

/-! ## Committed generation receipts beside the premise namespace -/

/-- Both parity sectors survive the scan, and charge conjugation
exchanges them.  Re-export of
`OPH.ExteriorSelection.parity_sectors_survive` and
`OPH.ExteriorSelection.conj_exchanges_survivors`. -/
theorem smStructure_generation_witness :
    (OPH.ExteriorSelection.chiral OPH.ExteriorSelection.evenMask = true
        ∧ OPH.ExteriorSelection.anomalyFree
            OPH.ExteriorSelection.evenMask = true
        ∧ OPH.ExteriorSelection.chiral OPH.ExteriorSelection.oddMask = true
        ∧ OPH.ExteriorSelection.anomalyFree
            OPH.ExteriorSelection.oddMask = true)
      ∧ ∀ i : Fin 10,
        OPH.ExteriorSelection.mem OPH.ExteriorSelection.evenMask i
          = OPH.ExteriorSelection.mem OPH.ExteriorSelection.oddMask
              (OPH.ExteriorSelection.conj i) :=
  ⟨OPH.ExteriorSelection.parity_sectors_survive,
    OPH.ExteriorSelection.conj_exchanges_survivors⟩

/-- In Standard-Model field order, the even exterior rows reproduce
the five matter weights of the descent theorem.  Re-export of
`OPH.QuantumMatterIntegration.even_component_weights_eq_matterWeights`. -/
theorem smStructure_generation_weights :
    [OPH.QuantumMatterIntegration.componentWeight 3,
      OPH.QuantumMatterIntegration.componentWeight 2,
      OPH.QuantumMatterIntegration.componentWeight 4,
      OPH.QuantumMatterIntegration.componentWeight 8,
      OPH.QuantumMatterIntegration.componentWeight 9] =
      OPH.Z6Descent.matterWeights :=
  OPH.QuantumMatterIntegration.even_component_weights_eq_matterWeights

/-! ## The Z6 kernel and the conditional global form (OL-G2) -/

/-- **The Z6 kernel (OL-G2, attained part).**  The common central
stabilizer of the realized matter weights is exactly the tensor-action
kernel; with the committed cyclicity receipt the kernel is the
six-element group generated by `(omega_3 I_3, -I_2, e^{i pi/3})`.
Re-export of `OPH.Z6Descent.kernel_on_realized_weights`. -/
theorem smStructure_z6_kernel :
    ∀ c : OPH.TraceBalancedKernel.C,
      (∀ w ∈ OPH.Z6Descent.realizedWeights,
        OPH.Z6Descent.phase c w = 0) ↔
        OPH.TraceBalancedKernel.tensorCharFun c = 0 :=
  OPH.Z6Descent.kernel_on_realized_weights

/-- **The conditional global-form status (OL-G2, boundary part).**
The selected diagonal kernel has exactly four subgroups, and every
realized local weight descends through all four, so the local tensor
table determines a maximal effective quotient without selecting a
physical global form.  The leading undischarged clause toward the
`S(U(3) x U(2))/Z6` quotient is the same-source loop-to-kernel
identity, a declared premise of the V3 program; the committed record
additionally names a source-derived complete character category,
instanton normalization, and laboratory flux attachment as premises of
the physical selection, and the global-form row is conditional on all
of them.  Re-export of
`OPH.Z6Descent.four_admissible_global_forms`. -/
theorem smStructure_global_form_menu :
    (Finset.univ.filter fun s : Finset (ZMod 6) =>
      (0 : ZMod 6) ∈ s ∧ ∀ a ∈ s, ∀ b ∈ s, a + b ∈ s).card = 4 :=
  OPH.Z6Descent.four_admissible_global_forms

/-- The committed abstract half of the loop-to-kernel premise: the
intertwiner from the six-axis gluing class onto the tensor-action
kernel carries the axis generator to the kernel generator and
antipodal reversal to conjugation, is injective, and has image exactly
the kernel.  The identification of the six-axis class with a physical
carrier-loop class of the same source is the declared premise; this
theorem supplies the abstract identification only.  Re-export of the
`OPH.Z6Descent.sixAxisToKernel` receipts. -/
theorem smStructure_loop_kernel_intertwiner :
    OPH.Z6Descent.sixAxisToKernel
        (QuotientAddGroup.mk (fun i => if i = 0 then 1 else 0))
          = OPH.TraceBalancedKernel.gen
      ∧ (∀ x : OPH.Z6Exact.Lattice ⧸ OPH.Z6Exact.gauge,
          OPH.Z6Descent.sixAxisToKernel (-x)
            = -OPH.Z6Descent.sixAxisToKernel x)
      ∧ Function.Injective OPH.Z6Descent.sixAxisToKernel
      ∧ ∀ c : OPH.TraceBalancedKernel.C,
          OPH.TraceBalancedKernel.tensorCharFun c = 0 ↔
            ∃ x, OPH.Z6Descent.sixAxisToKernel x = c :=
  ⟨OPH.Z6Descent.sixAxis_generator_maps_to_kernel_generator,
    OPH.Z6Descent.sixAxisToKernel_intertwines_involutions,
    OPH.Z6Descent.sixAxisToKernel_injective,
    OPH.Z6Descent.sixAxisToKernel_range⟩

/-! ## The family-band receipts (OL-G4) -/

/-- **The family-band selection (OL-G4, conditional).**  Under the
family band-selection rules, declared premises of the V3 program (the
family multiplicity object is a single complete faithful band of the
twelve-port adjacency with dimension in the pinned window, and
admissible bands are compared by the operational seam-cost order over
`Z[sqrt5]`), the admissible candidates number three and the rank-three
band is the unique strict cost minimizer.  The three-family row is
conditional on those registered but undischarged rules, and exactly those
rules plus the realization clause of the executable certificate
separate it from a derived family triplication.  Re-export of
`OPH.A5FamilyBand.family_band_selected`. -/
theorem smStructure_family_band_selected :
    (OPH.A5FamilyBand.candidates.filter
        OPH.A5FamilyBand.admissible).length = 3
      ∧ ((OPH.A5FamilyBand.candidates.filter
          OPH.A5FamilyBand.admissible).all fun c =>
            c.1 == "3" || OPH.A5FamilyBand.z5lt (5, -1) c.2.2.1) = true :=
  OPH.A5FamilyBand.family_band_selected

/-- The exact seam-cost order is strict, and without the faithfulness
clause the cost minimizer over all four bands is the trivial band at
cost zero: the clause is load-bearing, and cost minimization alone
selects nothing physical.  Re-export of
`OPH.A5FamilyBand.cost_order_strict` and
`OPH.A5FamilyBand.trivial_band_without_faithfulness`. -/
theorem smStructure_family_band_cost_receipts :
    (OPH.A5FamilyBand.z5lt (5, -1) (6, 0) = true
        ∧ OPH.A5FamilyBand.z5lt (6, 0) (5, 1) = true)
      ∧ (OPH.A5FamilyBand.candidates.all fun c =>
          c.1 == "1" || OPH.A5FamilyBand.z5lt (0, 0) c.2.2.1) = true :=
  ⟨OPH.A5FamilyBand.cost_order_strict,
    OPH.A5FamilyBand.trivial_band_without_faithfulness⟩

/-- The three-copy generation arithmetic on the selected band: the
integer-scaled anomaly forms of three copies vanish, and the diagonal
`Z6` generator fixes every field of the fifteen-state generation and
hence all three copies.  These are table receipts under register row
PR-12 and the declared band rules; no derived triplication is claimed.
Re-export of `OPH.A5FamilyBand.three_copy_anomaly_forms` and
`OPH.A5FamilyBand.diagonal_z6_fixes_three_copies`. -/
theorem smStructure_three_family_receipts :
    (3 * (6 * (1 : Int)^3 + 3 * (-4)^3 + 6^3 + 3 * 2^3
          + 2 * (-3)^3) = 0
        ∧ 3 * (6 * (1 : Int) + 3 * (-4) + 6 + 3 * 2 + 2 * (-3)) = 0
        ∧ 3 * ((-4 : Int) + 2 * 1 + 2) = 0
        ∧ 3 * (3 * (1 : Int) + (-3)) = 0)
      ∧ ((2 * (1 : Int) + 3 * 1 + 1) % 6 = 0
        ∧ (2 * (2 : Int) + 3 * 0 - 4) % 6 = 0
        ∧ (2 * (0 : Int) + 3 * 0 + 6) % 6 = 0
        ∧ (2 * (2 : Int) + 3 * 0 + 2) % 6 = 0
        ∧ (2 * (0 : Int) + 3 * 1 - 3) % 6 = 0
        ∧ 3 * 15 = 45) :=
  ⟨OPH.A5FamilyBand.three_copy_anomaly_forms,
    OPH.A5FamilyBand.diagonal_z6_fixes_three_copies⟩

/-! ## Boundary receipts: the surface carries its own limits -/

/-- **The chirality and global-form boundary (OL-G5, OL-G2
scoping).**  No central parameter acts by `-1` on all five matter
multiplets, so the finite central arithmetic supplies no spacetime
Spin or fermion-parity attachment; and charge conjugation exchanges
the two surviving parity sectors, so the scan does not choose between
them.  Re-export of `OPH.Z6Descent.no_universal_fermion_minus_one`
and `OPH.ExteriorSelection.conj_exchanges_survivors`. -/
theorem smStructure_boundary :
    (∀ c : OPH.TraceBalancedKernel.C,
        ¬ (∀ w ∈ OPH.Z6Descent.matterWeights,
          OPH.Z6Descent.phase c w = 3))
      ∧ ∀ i : Fin 10,
        OPH.ExteriorSelection.mem OPH.ExteriorSelection.evenMask i
          = OPH.ExteriorSelection.mem OPH.ExteriorSelection.oddMask
              (OPH.ExteriorSelection.conj i) :=
  ⟨OPH.Z6Descent.no_universal_fermion_minus_one,
    OPH.ExteriorSelection.conj_exchanges_survivors⟩

/-! ## The committed instance of the premise data

The barycentric one-third share, the equal face weight `1/20`, the
Standard-Model dimension data, and the even-parity selection mask
discharge every field with exact literals.  The instance is a
nonvacuity control for the premise structure; it selects nothing by
itself. -/

/-- The committed premise instance. -/
def committedSMStructurePremiseData : SMStructurePremiseData where
  vertexShare := 1 / 3
  faceWeight := fun _ => 1 / 20
  vertexShare_barycentric := rfl
  faceWeight_constant := fun _ _ _ _ => rfl
  faceWeight_normalized := by
    rw [Finset.sum_const, OPH.PrimitivePortDualMeasure.base_face_count,
      nsmul_eq_mul]
    norm_num
  centreDim := 1
  simpleIdealDims := {3, 8}
  centreDim_mem := by decide
  simpleIdealDims_mem := by
    intro d hd
    simp only [Multiset.insert_eq_cons, Multiset.mem_cons,
      Multiset.mem_singleton] at hd
    rcases hd with rfl | rfl <;> decide
  dims_total := by decide
  centreDim_le_one := le_refl 1
  centreDim_ne_zero := one_ne_zero
  selectionMask := ⟨OPH.ExteriorSelection.evenMask, by decide⟩
  selection_nonempty := by decide
  selection_chiral := OPH.ExteriorSelection.parity_sectors_survive.1
  selection_anomalyFree :=
    OPH.ExteriorSelection.parity_sectors_survive.2.1

/-- The pinned weight theorem instantiated at the committed premise
instance. -/
theorem committedSMStructurePremiseData_weight :
    ∀ p : Fin 12,
      portDualWeight committedSMStructurePremiseData.vertexShare
        committedSMStructurePremiseData.faceWeight p = 1 / 12 :=
  committedSMStructurePremiseData.port_dual_weight_pinned

/-! ## Composite receipt -/

/-- **The Standard Model structure surface receipt (issue #734).**
For every premise bundle: the classification inputs force centre
dimension one with simple-ideal dimensions `{3, 8}`; the weighting
rules pin the port-dual weight to `1/12`, and at that measure point
the family `G` is the unique nearest classified compact family; the
selection mask is one of the two conjugate parity sectors; the
realized-weight stabilizer is exactly the `Z6` tensor kernel; the
selected kernel admits exactly four global forms; and the rank-three
band is the unique strict cost minimizer among admissible bands.  The
first four conjuncts consume the bundled premises; the last three are
committed arithmetic receipts re-exported beside them.  The theorem type
does not consume or discharge PR-35 or PR-36.  Physical global-form and
family claims remain conditional on those rows and on the additional
attachment rows recorded in the premise register. -/
theorem smStructureAdequacySurface_receipt (D : SMStructurePremiseData) :
    (D.centreDim = 1 ∧ D.simpleIdealDims = {3, 8})
      ∧ (∀ p : Fin 12,
          portDualWeight D.vertexShare D.faceWeight p = 1 / 12)
      ∧ (∀ (p : Fin 12) (family : CompactFamily), family ≠ .G →
          dG2 ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)
              ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)
              ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ) <
            squaredDistanceAt
              ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)
              ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)
              ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)
              family)
      ∧ (D.selectionMask.val = OPH.ExteriorSelection.evenMask
          ∨ D.selectionMask.val = OPH.ExteriorSelection.oddMask)
      ∧ (∀ c : OPH.TraceBalancedKernel.C,
          (∀ w ∈ OPH.Z6Descent.realizedWeights,
            OPH.Z6Descent.phase c w = 0) ↔
            OPH.TraceBalancedKernel.tensorCharFun c = 0)
      ∧ ((Finset.univ.filter fun s : Finset (ZMod 6) =>
          (0 : ZMod 6) ∈ s ∧ ∀ a ∈ s, ∀ b ∈ s, a + b ∈ s).card = 4)
      ∧ ((OPH.A5FamilyBand.candidates.filter
            OPH.A5FamilyBand.admissible).length = 3
          ∧ ((OPH.A5FamilyBand.candidates.filter
              OPH.A5FamilyBand.admissible).all fun c =>
                c.1 == "3"
                  || OPH.A5FamilyBand.z5lt (5, -1) c.2.2.1) = true) :=
  ⟨D.forced_gauge_type,
    D.port_dual_weight_pinned,
    fun p family h => D.measure_selects_G p family h,
    D.matter_selection_is_parity_sector,
    OPH.Z6Descent.kernel_on_realized_weights,
    OPH.Z6Descent.four_admissible_global_forms,
    OPH.A5FamilyBand.family_band_selected⟩

end OPH.SMStructureAdequacySurface

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.SMStructureAdequacySurface.SMStructurePremiseData.centre_semisimple_enumeration
#print axioms OPH.SMStructureAdequacySurface.SMStructurePremiseData.forced_gauge_type
#print axioms OPH.SMStructureAdequacySurface.SMStructurePremiseData.port_dual_weight_pinned
#print axioms OPH.SMStructureAdequacySurface.SMStructurePremiseData.measure_selects_G
#print axioms OPH.SMStructureAdequacySurface.SMStructurePremiseData.matter_selection_is_parity_sector
#print axioms OPH.SMStructureAdequacySurface.SMStructurePremiseData.selection_parity_graded
#print axioms OPH.SMStructureAdequacySurface.smStructure_generation_witness
#print axioms OPH.SMStructureAdequacySurface.smStructure_generation_weights
#print axioms OPH.SMStructureAdequacySurface.smStructure_z6_kernel
#print axioms OPH.SMStructureAdequacySurface.smStructure_global_form_menu
#print axioms OPH.SMStructureAdequacySurface.smStructure_loop_kernel_intertwiner
#print axioms OPH.SMStructureAdequacySurface.smStructure_family_band_selected
#print axioms OPH.SMStructureAdequacySurface.smStructure_family_band_cost_receipts
#print axioms OPH.SMStructureAdequacySurface.smStructure_three_family_receipts
#print axioms OPH.SMStructureAdequacySurface.smStructure_boundary
#print axioms OPH.SMStructureAdequacySurface.committedSMStructurePremiseData
#print axioms OPH.SMStructureAdequacySurface.committedSMStructurePremiseData_weight
#print axioms OPH.SMStructureAdequacySurface.smStructureAdequacySurface_receipt
