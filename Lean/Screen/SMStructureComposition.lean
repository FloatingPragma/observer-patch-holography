import SMStructureAdequacySurface

open scoped BigOperators

/-!
# The composed Standard Model structure theorem (V3, issue #734)

One theorem on one premise bundle.  `Screen/SMStructureAdequacySurface.lean`
supplies the typed bundle `SMStructurePremiseData` (register rows
PR-09, PR-11, PR-12, with the measure-to-metric and nearest-point
rules of row PR-10 entering through the committed distance forms) and
proves its clauses separately.  This module closes the composition
contract of the issue: the measure-conditional selection of the
compact family `G` and the forced dimension type `(1, {3, 8})` meet
in a single uniqueness statement, the matter clauses are stated on
the bundled selection mask, and a two-model theorem identifies every
supplied family/mask pair that satisfies the composed predicate.

`genericTypeData` is the committed generic-stratum dimension column
of the pinned compact-locus classification certificate
`code/b14_jacobi/b14_compact_locus.certificate.json`
(`classification.compact_families`): the generic algebra of the
plane `P` is `so(3) + so(3) + R^6`, dimension data `(6, {3, 3})`,
and the generic algebra of each of `F` and `G` is `su(3) + so(3) + R`,
dimension data `(1, {3, 8})`.  The table is consumed under register
row PR-11 exactly as the component tables of `ExteriorSelection` are
consumed under register row PR-12: a committed literal reading of a
pinned certificate, with its receipts beside it.  The mirror tie
`genericTypeData .F = genericTypeData .G` is one of those receipts:
the type reading alone separates nothing inside the compact
`su(3) + so(3)` stratum, the measure clause of rows PR-09/PR-10 is
load-bearing, and `F_wins_at_witness` exhibits an invariant carrier
metric at which `F` is nearer.

Boundary.  The composed theorem forces dimension data and identifies
the nearest classified family at the bundled measure; it does not
construct a source bracket, current, coupling, action, Lagrangian
assembly, symmetry breaking, or laboratory identification.  The
committed gap theorem `OPHGap.gModule_iso_does_not_determine_bracket`
scopes the type claim: a `G`-module isomorphism determines no Lie
bracket, and the source current and action attachment stays on the
register as row PR-54.  The global-form clauses of OL-G2 stay
outside every theorem here; they are conditional on register rows
PR-35 and PR-46 and are stated on the surface only.  The two-model
theorem quantifies over the classified compact locus and the
1024-mask grammar under the bundled premises; no wider model class
is claimed, and the compact locus is conditional on the pinned
classification certificate of row PR-11.
-/

namespace OPH.SMStructureComposition

open OPH.SMStructureAdequacySurface
open OPH.PrimitivePortDualMeasure (portDualWeight)
open OPH.OrientedFaceInvariantMetric (dG2 squaredDistanceAt)
open OPH.OrientedFaceBracketSelector (CompactFamily)

/-- The generic-stratum dimension data `(centre dimension,
simple-ideal dimension multiset)` of each classified compact family:
the `classification.compact_families` column of the pinned
certificate `code/b14_jacobi/b14_compact_locus.certificate.json`,
committed as a literal table under register row PR-11. -/
def genericTypeData : CompactFamily → ℕ × Multiset ℕ
  | .P => (6, {3, 3})
  | .F => (1, {3, 8})
  | .G => (1, {3, 8})

/-- No-cheat receipt for the committed table: `F` and `G` share the
generic dimension data, the two families are distinct, and there is
an invariant carrier metric at which `F` is strictly nearer than
`G`.  The type reading alone selects nothing; the bundled measure
clause of rows PR-09/PR-10 is load-bearing. -/
theorem genericTypeData_mirror_tie :
    genericTypeData .F = genericTypeData .G
      ∧ (CompactFamily.F ≠ CompactFamily.G)
      ∧ OPH.OrientedFaceInvariantMetric.dF2 8 1 1 <
          OPH.OrientedFaceInvariantMetric.dG2 8 1 1 :=
  ⟨rfl, by decide, OPH.OrientedFaceInvariantMetric.F_wins_at_witness⟩

/-- Classification-locus/enumeration compatibility: the committed
dimension data of every classified family satisfies the hypotheses
of the committed centre-versus-semisimple enumeration (register row
PR-11). -/
theorem genericTypeData_satisfies_enumeration_premises
    (family : CompactFamily) :
    (genericTypeData family).1 ∈ ({0, 1, 5, 6, 7, 11, 12} : Finset ℕ)
      ∧ (∀ d ∈ (genericTypeData family).2, d ∈ OPHSelection.SimpleDims)
      ∧ (genericTypeData family).1 + (genericTypeData family).2.sum
          = 12 := by
  cases family <;>
    exact ⟨by decide, fun d hd => by
      simp only [genericTypeData, Multiset.insert_eq_cons,
        Multiset.mem_cons, Multiset.mem_singleton] at hd
      rcases hd with rfl | rfl <;> decide, by decide⟩

/-- The committed table lands inside the four survivors of the
centre-versus-semisimple enumeration: the classification locus and
the enumeration arithmetic agree row by row. -/
theorem genericTypeData_in_enumeration (family : CompactFamily) :
    ((genericTypeData family).1 = 12 ∧ (genericTypeData family).2 = 0)
      ∨ ((genericTypeData family).1 = 6
          ∧ (genericTypeData family).2 = {3, 3})
      ∨ ((genericTypeData family).1 = 1
          ∧ (genericTypeData family).2 = {3, 8})
      ∨ ((genericTypeData family).1 = 0
          ∧ (genericTypeData family).2 = {3, 3, 3, 3}) := by
  have h := genericTypeData_satisfies_enumeration_premises family
  exact OPH.TrichotomyCases.centre_case_enumeration _ _ h.1 h.2.1 h.2.2

/-- The real port-dual weight of the bundle at a port: the bundled
PR-09 weighting read through the diagonal measure-to-metric rule of
register row PR-10. -/
noncomputable def bundleWeight (D : SMStructurePremiseData)
    (p : Fin 12) : ℝ :=
  ((portDualWeight D.vertexShare D.faceWeight p : ℚ) : ℝ)

/-- Strict nearest-point predicate at the bundled measure: at every
port, the family is strictly nearer to the pinned face bracket than
every other classified family under the metric induced by the
bundled port-dual weight (register rows PR-09 and PR-10; the compact
locus is conditional on the pinned classification certificate of row
PR-11). -/
def StrictlyNearestAt (D : SMStructurePremiseData)
    (family : CompactFamily) : Prop :=
  ∀ p : Fin 12, ∀ other : CompactFamily, other ≠ family →
    squaredDistanceAt (bundleWeight D p) (bundleWeight D p)
        (bundleWeight D p) family <
      squaredDistanceAt (bundleWeight D p) (bundleWeight D p)
        (bundleWeight D p) other

/-- The family `G` is strictly nearest at the bundled measure. -/
theorem strictlyNearestAt_G (D : SMStructurePremiseData) :
    StrictlyNearestAt D .G := by
  intro p other hother
  exact D.measure_selects_G p other hother

/-- Only `G` is strictly nearest at the bundled measure: two distinct
strict minimizers at the same port contradict each other. -/
theorem strictlyNearestAt_eq_G (D : SMStructurePremiseData)
    (family : CompactFamily) (h : StrictlyNearestAt D family) :
    family = .G := by
  by_contra hne
  have hG : dG2 (bundleWeight D 0) (bundleWeight D 0)
      (bundleWeight D 0) <
      squaredDistanceAt (bundleWeight D 0) (bundleWeight D 0)
        (bundleWeight D 0) family :=
    D.measure_selects_G 0 family hne
  have hf : squaredDistanceAt (bundleWeight D 0) (bundleWeight D 0)
      (bundleWeight D 0) family <
      dG2 (bundleWeight D 0) (bundleWeight D 0) (bundleWeight D 0) :=
    h 0 .G (Ne.symm hne)
  exact lt_asymm hG hf

/-- **The composed family selection (OL-G1).**  Under the bundle
there is exactly one classified compact family that is both strictly
nearest to the pinned face bracket at the bundled port-dual measure
(rows PR-09, PR-10) and whose committed generic dimension data
equals the forced classification data of the bundle (row PR-11): the
family `G`.  The metric clause and the type clause share the single
witness. -/
theorem composed_family_unique (D : SMStructurePremiseData) :
    ∃! family : CompactFamily,
      StrictlyNearestAt D family ∧
        genericTypeData family = (D.centreDim, D.simpleIdealDims) := by
  obtain ⟨h1, h2⟩ := D.forced_gauge_type
  refine ⟨.G, ⟨strictlyNearestAt_G D, ?_⟩, ?_⟩
  · rw [h1, h2]; rfl
  · rintro family ⟨hnear, -⟩
    exact strictlyNearestAt_eq_G D family hnear

/-! ## The two-model capstone -/

/-- The composed Standard-Model structure predicate at a supplied
family/mask pair, relative to the bundle: the family is strictly
nearest at the bundled measure and carries the forced generic
dimension data, and the mask satisfies the register row PR-12
grammar (nonempty, chiral, anomaly free).  The mask constraints
mirror the bundle's own PR-12 fields, restated on the supplied
mask. -/
def IsSMStructureModel (D : SMStructurePremiseData)
    (family : CompactFamily) (m : Fin 1024) : Prop :=
  StrictlyNearestAt D family
    ∧ genericTypeData family = (D.centreDim, D.simpleIdealDims)
    ∧ m.val ≠ 0
    ∧ OPH.ExteriorSelection.chiral m.val = true
    ∧ OPH.ExteriorSelection.anomalyFree m.val = true

/-- **The two-model theorem.**  Within the classified compact locus
and the 1024-mask grammar, under the bundled premises PR-09 through
PR-12, the composed predicate has exactly two models: the family `G`
with the even-parity mask and the family `G` with the odd-parity
mask.  No wider model class is claimed, and the compact locus is
conditional on the pinned classification certificate of row
PR-11. -/
theorem smStructure_two_models (D : SMStructurePremiseData) :
    ∀ (family : CompactFamily) (m : Fin 1024),
      IsSMStructureModel D family m ↔
        ((family = .G ∧ m.val = OPH.ExteriorSelection.evenMask)
          ∨ (family = .G ∧ m.val = OPH.ExteriorSelection.oddMask)) := by
  intro family m
  constructor
  · rintro ⟨hnear, -, hne, hchi, hanom⟩
    have hfam := strictlyNearestAt_eq_G D family hnear
    rcases OPH.ExteriorSelection.selection_unique m hne hchi hanom
      with h | h
    · exact Or.inl ⟨hfam, h⟩
    · exact Or.inr ⟨hfam, h⟩
  · rintro (⟨rfl, hm⟩ | ⟨rfl, hm⟩) <;>
      obtain ⟨h1, h2⟩ := D.forced_gauge_type
    · refine ⟨strictlyNearestAt_G D, by rw [h1, h2]; rfl, ?_, ?_, ?_⟩
      · rw [hm]; decide
      · rw [hm]; exact OPH.ExteriorSelection.parity_sectors_survive.1
      · rw [hm]; exact OPH.ExteriorSelection.parity_sectors_survive.2.1
    · refine ⟨strictlyNearestAt_G D, by rw [h1, h2]; rfl, ?_, ?_, ?_⟩
      · rw [hm]; decide
      · rw [hm]; exact OPH.ExteriorSelection.parity_sectors_survive.2.2.1
      · rw [hm]; exact OPH.ExteriorSelection.parity_sectors_survive.2.2.2

/-- Both models are realized and charge conjugation exchanges their
masks componentwise: the composed surface is unique up to charge
conjugation within the stated locus and grammar.  The theorem does
not choose between the two models. -/
theorem smStructure_conj_exchanges_models (D : SMStructurePremiseData) :
    IsSMStructureModel D .G
        ⟨OPH.ExteriorSelection.evenMask, by decide⟩
      ∧ IsSMStructureModel D .G
          ⟨OPH.ExteriorSelection.oddMask, by decide⟩
      ∧ ∀ i : Fin 10,
          OPH.ExteriorSelection.mem OPH.ExteriorSelection.evenMask i
            = OPH.ExteriorSelection.mem OPH.ExteriorSelection.oddMask
                (OPH.ExteriorSelection.conj i) :=
  ⟨(smStructure_two_models D .G _).mpr (Or.inl ⟨rfl, rfl⟩),
    (smStructure_two_models D .G _).mpr (Or.inr ⟨rfl, rfl⟩),
    OPH.ExteriorSelection.conj_exchanges_survivors⟩

/-! ## The composed receipt -/

/-- **The composed Standard Model structure receipt (issue #734).**
For every premise bundle `D` (register rows PR-09 through PR-12):
exactly one classified compact family is strictly nearest at the
bundled measure with generic dimension data equal to the bundle's
forced classification data; that data is `(1, {3, 8})`; the bundled
selection mask is one of the two conjugate parity sectors; its
selected components number fifteen in total dimension; its four
anomaly forms vanish; it carries the parity grading; and the common
central stabilizer of its selected component weights is exactly the
`Z6` tensor-action kernel.  Every clause consumes the bundle.  The
theorem does not consume or discharge PR-35, PR-36, PR-46, PR-47, or
PR-54; physical global-form, family, chirality, and action claims
stay conditional on those rows. -/
theorem smStructureComposition_receipt (D : SMStructurePremiseData) :
    (∃! family : CompactFamily,
        StrictlyNearestAt D family ∧
          genericTypeData family = (D.centreDim, D.simpleIdealDims))
      ∧ (D.centreDim, D.simpleIdealDims)
          = ((1 : ℕ), ({3, 8} : Multiset ℕ))
      ∧ (D.selectionMask.val = OPH.ExteriorSelection.evenMask
          ∨ D.selectionMask.val = OPH.ExteriorSelection.oddMask)
      ∧ OPH.ExteriorSelection.maskSum
          (fun i => OPH.ExteriorSelection.colorDim i *
            OPH.ExteriorSelection.weakDim i)
          D.selectionMask.val = 15
      ∧ (OPH.ExteriorSelection.grav D.selectionMask.val = 0
          ∧ OPH.ExteriorSelection.su3 D.selectionMask.val = 0
          ∧ OPH.ExteriorSelection.su2 D.selectionMask.val = 0
          ∧ OPH.ExteriorSelection.cube D.selectionMask.val = 0)
      ∧ ((∀ i : Fin 10,
            OPH.ExteriorSelection.mem D.selectionMask.val i
              = !OPH.ExteriorSelection.odd i)
          ∨ (∀ i : Fin 10,
            OPH.ExteriorSelection.mem D.selectionMask.val i
              = OPH.ExteriorSelection.odd i))
      ∧ (∀ c : OPH.TraceBalancedKernel.C,
          (∀ i : Fin 10,
            OPH.ExteriorSelection.mem D.selectionMask.val i = true →
              OPH.Z6Descent.phase c
                (OPH.QuantumMatterIntegration.componentWeight i) = 0) ↔
            OPH.TraceBalancedKernel.tensorCharFun c = 0) :=
  ⟨composed_family_unique D,
    by rw [D.forced_gauge_type.1, D.forced_gauge_type.2],
    D.matter_selection_is_parity_sector,
    D.masked_generation_count,
    D.masked_anomaly_forms,
    D.selection_parity_graded,
    D.masked_z6_kernel⟩

end OPH.SMStructureComposition

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.SMStructureComposition.genericTypeData_mirror_tie
#print axioms OPH.SMStructureComposition.genericTypeData_satisfies_enumeration_premises
#print axioms OPH.SMStructureComposition.genericTypeData_in_enumeration
#print axioms OPH.SMStructureComposition.strictlyNearestAt_G
#print axioms OPH.SMStructureComposition.strictlyNearestAt_eq_G
#print axioms OPH.SMStructureComposition.composed_family_unique
#print axioms OPH.SMStructureComposition.smStructure_two_models
#print axioms OPH.SMStructureComposition.smStructure_conj_exchanges_models
#print axioms OPH.SMStructureComposition.smStructureComposition_receipt
