import QFT.ColimitNormedCompletion
import Mathlib.Analysis.CStarAlgebra.GelfandNaimarkSegal

set_option autoImplicit false
set_option relaxedAutoImplicit false

/-!
# A conditional GNS representation of the completed observer colimit

`QFT.ColimitNormedCompletion` constructs a C*-algebra by completing the
filtered colimit of the finite private matrix algebras.  This module attaches
that completed algebra, and each of its completed local algebras, to the GNS
Hilbert space of an **explicitly supplied** positive linear functional for
Mathlib's canonical C*-spectral order on the completion.

The input functional is not derived here.  In particular, this module does
not select a physical state, prove that it is a vacuum, prove faithfulness of the GNS
representation, or attach the Hilbert space to physical spacetime, time,
energy, momentum, fields, particles, or scattering.  It does not discharge
PR-52 or PR-58 and it does not promote OL-C6.  The Cauchy statement transported
below is only the committed order-theoretic equality of completed local
algebras; it is not a physical time-slice theorem.

Mathlib's GNS construction supplies the Hilbert completion and the unital star
representation.  The contribution here is the typed adapter to the committed
OPH colimit and regional-net objects, with the existing isotony, locality, and
order-theoretic Cauchy equality carried through that adapter.
-/

namespace OPH.QFT

open OPH.Tower
open UniformSpace

open scoped ComplexOrder InnerProductSpace OPH.QFT.CompletionStar

universe u

variable {ι : Type u} [Preorder ι]
variable {T : ConsensusTower ι}

/-- The completed colimit algebra, named locally for the GNS adapter. -/
abbrev ColimitCompletion (T : ConsensusTower ι) :=
  Completion (TowerColimit T)

/- Mathlib deliberately does not register a global order on every C*-algebra.
For this adapter we use its canonical spectral order, locally, so the choice
cannot leak into unrelated modules. -/
noncomputable local instance colimitCompletionPartialOrder :
    PartialOrder (ColimitCompletion T) :=
  CStarAlgebra.spectralOrder (ColimitCompletion T)

local instance colimitCompletionStarOrderedRing :
    StarOrderedRing (ColimitCompletion T) :=
  CStarAlgebra.spectralOrderedRing (ColimitCompletion T)

/-- The GNS Hilbert space of an explicitly supplied positive functional on
the completed colimit. -/
abbrev ColimitGNS (ω : ColimitCompletion T →ₚ[ℂ] ℂ) := ω.GNS

/-- The bounded operators on the conditional GNS Hilbert space. -/
abbrev ColimitGNSOperators (ω : ColimitCompletion T →ₚ[ℂ] ℂ) :=
  ColimitGNS ω →L[ℂ] ColimitGNS ω

/-- The conditional unital star representation of the completed colimit on
the GNS Hilbert space of the supplied positive functional. -/
noncomputable def colimitGNSRepresentation
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ) :
    ColimitCompletion T →⋆ₐ[ℂ] ColimitGNSOperators ω :=
  ω.gnsStarAlgHom

@[simp]
theorem colimitGNSRepresentation_apply
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ) (a : ColimitCompletion T) :
    colimitGNSRepresentation ω a = ω.gnsStarAlgHom a :=
  rfl

/-- Normalization is an additional hypothesis on the supplied functional;
positivity alone does not provide it. -/
def IsNormalizedColimitFunctional
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ) : Prop :=
  ω 1 = 1

/-- The completion image of the unit class in the GNS pre-Hilbert space.
It is deliberately not called a vacuum: no dynamics or energy operator has
been supplied. -/
noncomputable def colimitGNSUnitClass
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ) : ColimitGNS ω :=
  (ω.toPreGNS (1 : ColimitCompletion T) : ColimitGNS ω)

/-- Acting on the unit class recovers the canonical dense copy of the source
algebra in its GNS completion. -/
@[simp]
theorem colimitGNSRepresentation_apply_unitClass
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ) (a : ColimitCompletion T) :
    colimitGNSRepresentation ω a (colimitGNSUnitClass ω) =
      (ω.toPreGNS a : ColimitGNS ω) := by
  change ω.gnsNonUnitalStarAlgHom a
      (ω.toPreGNS (1 : ColimitCompletion T) : ColimitGNS ω) =
    (ω.toPreGNS a : ColimitGNS ω)
  rw [PositiveLinearMap.gnsNonUnitalStarAlgHom_apply_coe]
  simp [PositiveLinearMap.leftMulMapPreGNS]

/-- The supplied functional is recovered as a matrix coefficient of the unit
class.  This is an algebraic GNS fact, not a vacuum-expectation claim. -/
theorem colimitGNSUnitClass_expectation
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ) (a : ColimitCompletion T) :
    ⟪colimitGNSUnitClass ω,
      colimitGNSRepresentation ω a (colimitGNSUnitClass ω)⟫_ℂ = ω a := by
  rw [colimitGNSRepresentation_apply_unitClass]
  simp [colimitGNSUnitClass, PositiveLinearMap.preGNS_inner_def]

/-- Normalization of the supplied functional makes the unit class a unit
vector.  Positivity alone is intentionally insufficient for this theorem. -/
theorem norm_colimitGNSUnitClass
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ)
    (hω : IsNormalizedColimitFunctional ω) :
    ‖colimitGNSUnitClass ω‖ = 1 := by
  change ω 1 = 1 at hω
  have hinner := colimitGNSUnitClass_expectation ω
    (1 : ColimitCompletion T)
  have hsq : ‖colimitGNSUnitClass ω‖ ^ 2 = 1 := by
    rw [norm_sq_eq_re_inner (𝕜 := ℂ)]
    simpa [hω] using congrArg Complex.re hinner
  nlinarith [norm_nonneg (colimitGNSUnitClass ω)]

/-- The orbit of the unit class under the whole completed colimit is dense in
the GNS completion.  This mathematical cyclicity does not select a physical
vacuum or supply dynamics. -/
theorem denseRange_colimitGNSRepresentation_apply_unitClass
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ) :
    DenseRange (fun a : ColimitCompletion T =>
      colimitGNSRepresentation ω a (colimitGNSUnitClass ω)) := by
  have horbit :
      (fun a : ColimitCompletion T =>
          colimitGNSRepresentation ω a (colimitGNSUnitClass ω)) =
        ((↑) : ω.PreGNS → ColimitGNS ω) ∘ ω.toPreGNS := by
    funext a
    exact colimitGNSRepresentation_apply_unitClass ω a
  rw [horbit]
  exact Completion.denseRange_coe.comp ω.toPreGNS.surjective.denseRange
    (Completion.continuous_coe ω.PreGNS)

/-! ## The represented regional net -/

/-- The image of a completed regional algebra in the conditional GNS
representation.  Taking the image (rather than silently identifying it with
the source) does not assume that the representation is faithful. -/
noncomputable def representedCompletionLocalAlgebra
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ)
    (N : FiniteCausalObserverNet T) (U : EventRegion N) :
    StarSubalgebra ℂ (ColimitGNSOperators ω) :=
  (completionLocalAlgebra N U).map (colimitGNSRepresentation ω)

/-- A represented observable belongs to the represented regional algebra
whenever its source observable belongs to the completed regional algebra. -/
theorem colimitGNSRepresentation_mem_representedCompletionLocalAlgebra
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ)
    {N : FiniteCausalObserverNet T} {U : EventRegion N}
    {a : ColimitCompletion T} (ha : a ∈ completionLocalAlgebra N U) :
    colimitGNSRepresentation ω a ∈ representedCompletionLocalAlgebra ω N U :=
  StarSubalgebra.mem_map.mpr ⟨a, ha, rfl⟩

/-- Isotony of completed regional algebras survives passage to their images
under the conditional GNS representation. -/
theorem representedCompletionLocalAlgebra_mono
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ)
    {N : FiniteCausalObserverNet T} {U V : EventRegion N} (f : U ⟶ V) :
    representedCompletionLocalAlgebra ω N U ≤
      representedCompletionLocalAlgebra ω N V :=
  StarSubalgebra.map_mono (completionLocalAlgebra_mono f)

/-- Refinement independence survives passage to represented regional
algebras. -/
theorem representedCompletionLocalAlgebra_refine
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ)
    (N : FiniteCausalObserverNet T) {r s : ι} (hrs : r ≤ s) (U : N.Region r) :
    representedCompletionLocalAlgebra ω N ⟨s, N.regionRefine hrs U⟩ =
      representedCompletionLocalAlgebra ω N ⟨r, U⟩ := by
  unfold representedCompletionLocalAlgebra
  rw [completionLocalAlgebra_refine hrs U]

/-- Elementwise locality survives passage to the represented regional
algebras.  This conclusion needs no faithfulness hypothesis. -/
theorem representedCompletionLocalAlgebra_commute
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ)
    {N : FiniteCausalObserverNet T} {r : ι} {U V : N.Region r}
    (hUV : N.disjoint r U V) {x y : ColimitGNSOperators ω}
    (hx : x ∈ representedCompletionLocalAlgebra ω N ⟨r, U⟩)
    (hy : y ∈ representedCompletionLocalAlgebra ω N ⟨r, V⟩) :
    Commute x y := by
  obtain ⟨a, ha, rfl⟩ := StarSubalgebra.mem_map.mp hx
  obtain ⟨b, hb, rfl⟩ := StarSubalgebra.mem_map.mp hy
  exact (completionLocalAlgebra_commute hUV ha hb).map
    (colimitGNSRepresentation ω)

/-- The committed order-theoretic Cauchy equality survives passage to the
represented regional algebras.  This is not a physical time-slice theorem. -/
theorem representedCompletionLocalAlgebra_eq_of_isCauchyEmbedding
    (ω : ColimitCompletion T →ₚ[ℂ] ℂ)
    {N : FiniteCausalObserverNet T} {U V : EventRegion N}
    (f : U ⟶ V) (hf : IsCauchyEmbedding f) :
    representedCompletionLocalAlgebra ω N U =
      representedCompletionLocalAlgebra ω N V := by
  unfold representedCompletionLocalAlgebra
  rw [completionLocalAlgebra_eq_of_isCauchyEmbedding f hf]

end OPH.QFT

-- Axiom audit: this adapter must stay on the standard Mathlib basis.
#print axioms OPH.QFT.colimitGNSRepresentation
#print axioms OPH.QFT.colimitGNSRepresentation_apply
#print axioms OPH.QFT.colimitGNSRepresentation_apply_unitClass
#print axioms OPH.QFT.colimitGNSUnitClass_expectation
#print axioms OPH.QFT.norm_colimitGNSUnitClass
#print axioms OPH.QFT.denseRange_colimitGNSRepresentation_apply_unitClass
#print axioms OPH.QFT.representedCompletionLocalAlgebra_mono
#print axioms OPH.QFT.representedCompletionLocalAlgebra_refine
#print axioms OPH.QFT.representedCompletionLocalAlgebra_commute
#print axioms OPH.QFT.representedCompletionLocalAlgebra_eq_of_isCauchyEmbedding
