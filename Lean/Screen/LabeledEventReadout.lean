import Mathlib
import ObserverPatchHolography.ScalarSeamRepair

namespace OPH.LabeledEventReadout

open ObserverPatchHolography.ScalarSeamRepair

noncomputable section

/-!
# Labeled event readout of one refinement layer

A refinement layer contains twelve inherited slots and thirty midpoint slots.
The base map reads the inherited slots.  For each midpoint, a declared labeled
event pair-averages that midpoint with one chosen inherited parent and exposes
the parent's updated value through the base map.

This file proves the exact decoder for any choice of parent.  It gives an
equivalence between the 42 scalar field coordinates and the twelve baseline
responses plus thirty labeled event responses.  The geometry-specific carrier
enumeration, A5 transport checks, exact rank calculation, and independent
replay are supplied by the OPH-FPE labeled-event certificate.

The event responses must all refer to the same pre-event field.  This theorem
does not construct repeatable preparation, reset or checkpoint access, a
nondestructive instrument, or a physical observable.  It also makes no claim
about passive averaged-semigroup readout protocols.
-/

/-- The first refinement layer, typed as inherited and midpoint slots. -/
abbrev Carrier := Sum (Fin 12) (Fin 30)

/-- A scalar field on the typed first-refinement carrier. -/
abbrev Field := Carrier → ℝ

/-- Twelve baseline values and thirty selected labeled-event responses. -/
abbrev Readout := (Fin 12 → ℝ) × (Fin 30 → ℝ)

/-- Restriction to the twelve inherited carrier slots. -/
def baseline (x : Field) : Fin 12 → ℝ :=
  fun u ↦ x (.inl u)

/-- The labeled intervention for midpoint `m`: pair-average it with the
chosen inherited parent and fix every other carrier coordinate. -/
def labeledEvent (parent : Fin 30 → Fin 12) (m : Fin 30) :
    Field →ₗ[ℝ] Field :=
  pairAverage (.inl (parent m)) (.inr m)

/-- The response coordinate exposed at the chosen inherited parent after the
labeled midpoint event. -/
def eventResponse (parent : Fin 30 → Fin 12) (x : Field) (m : Fin 30) : ℝ :=
  baseline (labeledEvent parent m x) (parent m)

@[simp]
theorem eventResponse_eq_average
    (parent : Fin 30 → Fin 12) (x : Field) (m : Fin 30) :
    eventResponse parent x m =
      (x (.inl (parent m)) + x (.inr m)) / 2 := by
  simp [eventResponse, baseline, labeledEvent]

/-- The complete minimal labeled readout: twelve baseline coordinates and one
chosen event response for each of the thirty midpoint slots. -/
def measure (parent : Fin 30 → Fin 12) (x : Field) : Readout :=
  (baseline x, eventResponse parent x)

/-- Exact decoder.  Inherited coordinates come directly from the baseline;
each midpoint is twice its event response minus its chosen parent's baseline. -/
def reconstruct (parent : Fin 30 → Fin 12) (y : Readout) : Field
  | .inl u => y.1 u
  | .inr m => 2 * y.2 m - y.1 (parent m)

@[simp]
theorem reconstruct_measure
    (parent : Fin 30 → Fin 12) (x : Field) :
    reconstruct parent (measure parent x) = x := by
  funext slot
  cases slot with
  | inl u => simp [reconstruct, measure, baseline]
  | inr m =>
      simp [reconstruct, measure, baseline, eventResponse_eq_average]
      ring

@[simp]
theorem measure_reconstruct
    (parent : Fin 30 → Fin 12) (y : Readout) :
    measure parent (reconstruct parent y) = y := by
  rcases y with ⟨base, events⟩
  apply Prod.ext
  · funext u
    simp [measure, baseline, reconstruct]
  · funext m
    simp [measure, eventResponse_eq_average, reconstruct]

/-- The minimal labeled readout is exactly invertible for every parent choice.
The parent selector need not be canonical or symmetry invariant. -/
def readoutEquiv (parent : Fin 30 → Fin 12) : Field ≃ Readout where
  toFun := measure parent
  invFun := reconstruct parent
  left_inv := reconstruct_measure parent
  right_inv := measure_reconstruct parent

theorem measure_injective (parent : Fin 30 → Fin 12) :
    Function.Injective (measure parent) :=
  (readoutEquiv parent).injective

theorem measure_surjective (parent : Fin 30 → Fin 12) :
    Function.Surjective (measure parent) :=
  (readoutEquiv parent).surjective

/-- The typed carrier contains exactly the twelve inherited and thirty
midpoint slots used by the executable rank certificate. -/
theorem carrier_card : Fintype.card Carrier = 42 := by
  simp [Carrier]

#print axioms OPH.LabeledEventReadout.eventResponse_eq_average
#print axioms OPH.LabeledEventReadout.reconstruct_measure
#print axioms OPH.LabeledEventReadout.measure_reconstruct
#print axioms OPH.LabeledEventReadout.measure_injective
#print axioms OPH.LabeledEventReadout.carrier_card

end

end OPH.LabeledEventReadout
