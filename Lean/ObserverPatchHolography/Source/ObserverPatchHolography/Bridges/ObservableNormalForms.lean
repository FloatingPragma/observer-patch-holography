import ObserverPatchHolography.Primitives
import ObservableNormalForms.ObserverConfluence

/-!
# Bridge to observation-determined normal forms

This module records the exact relationship between Jonathan Hill's concrete
local-repair development and the substrate-neutral endpoint theorem used by
the Observation-Determined Normal Forms paper.  It does not assert that any
particular physical boundary map identifies the concrete consistent quotient;
that remains a separate application theorem.
-/

namespace OPH

universe u

open ObservableNormalForms

section LocalRepairBridge

variable {C : OPHCarrier}
variable (lr : C.Patch → Records C → Records C)
variable
  (H1 : ∀ (i : C.Patch) (x : Records C) (j : C.Patch),
    j ≠ i → (lr i x) j = x j)
  (H2 : ∀ (i : C.Patch) (x : Records C),
    lr i x ≠ x ↔
      ∃ e : C.Edge,
        (C.src e = i ∨ C.tgt e = i) ∧ ¬ edgeConsistentAt e x)
  (H3 : ∀ (i : C.Patch) (x : Records C),
    lr i x ≠ x →
      ∀ e : C.Edge,
        (C.src e = i ∨ C.tgt e = i) → edgeConsistentAt e (lr i x))

/-- A boundary invariant for every local repair is an observation-preservation
certificate for the induced rewrite relation. -/
theorem acceptedStepLR_observationPreserving
    {β : Type u} (B : Records C → β)
    (HB : ∀ (i : C.Patch) (x : Records C), B (lr i x) = B x) :
    ObservationPreserving (acceptedStepLR lr) B := by
  intro x y hxy
  obtain ⟨i, rfl, _⟩ := hxy
  exact (HB i x).symm

/-- Jonathan's H1--H3 completeness theorem is exactly the `CompleteFor`
premise of the neutral normal-form theorem. -/
theorem acceptedStepLR_completeFor :
    CompleteFor (acceptedStepLR lr) {x | Consistent C x} := by
  intro x
  change NormalFormLR lr x ↔ Consistent C x
  exact completeness lr H1 H2 H3 x

/-- Boundary identification on the concrete consistent records is equivalent
to cross-source endpoint uniqueness for the concrete local-repair relation.

This recovers the logical core of `boundary_fiber_observer_unique` as the
forward implication and additionally exposes the converse.  It deliberately
leaves the boundary map and its injectivity proof as explicit inputs. -/
theorem boundaryIdentifiesModulo_iff_observerEndpointUniqueModuloLR
    {β : Type u} (B : Records C → β)
    (HB : ∀ (i : C.Patch) (x : Records C), B (lr i x) = B x) :
    BoundaryIdentifiesModulo {x | Consistent C x} B (gaugeEquiv C) ↔
      ObserverEndpointUniqueModulo (acceptedStepLR lr) B (gaugeEquiv C) :=
  boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo
    (acceptedStepLR_observationPreserving lr B HB)
    (acceptedStepLR_completeFor lr H1 H2 H3)

end LocalRepairBridge

#print axioms boundaryIdentifiesModulo_iff_observerEndpointUniqueModuloLR

end OPH
