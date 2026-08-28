import Tower.AdaptiveFixedPointEndpoint
import ObserverPatchHolography.Execution.CumulativeAttemptCapacity

/-!
# Capacity-bounded adaptive public endpoints

This module is the downstream consumer of cumulative attempt capacity. A
`BoundedWaste` scheduler and a sufficient attempt budget reach a stable normal
form before exhaustion. Completeness places that record in the existing
`FixedPointObject`; confluence identifies it with the canonical raw and public
endpoint.

No theorem selects the scheduler, supplies the budget, or attaches the attempt
unit to physical time, rate, energy, bandwidth, fees, or hardware.
-/

namespace OPH.Tower.CumulativeCapacityEndpoint

open OPH
open OPH.Execution
open OPH.Locality
open OPH.Tower.AdaptiveEndpoint
open OPH.Tower.PublicWorldPresentation
open OPH.Tower.PublicWorldPresentation.FiniteRepairSystem
open OPH.Tower.PublicWorldPresentation.OPHPrimitiveEndpoint

noncomputable section

variable (C : OPHCarrier)
variable [FinitePatchStates C]

theorem boundedWaste_endpoint_within_attempt_budget
    (seed start : Records C) (sigma : AdaptiveScheduler C)
    (waste budget : Nat)
    (hcomplete : Completeness C)
    (hbounded : BoundedWaste C waste sigma)
    (hbudget : (waste + 1) * mismatchCount start ≤ budget) :
    ∃ N, cumulativeAttemptCost N ≤ budget ∧
      Consistent C (adaptiveRun N sigma start) ∧
      (∀ n, N ≤ n →
        adaptiveRun n sigma start = adaptiveRun N sigma start) ∧
      ∃ w : FixedPointObject (presentation C seed) (localSystem C seed),
        w.1 = (presentation C seed).toPublicWorld
          (adaptiveRun N sigma start) := by
  obtain ⟨N, hN, hnormal, hconstant⟩ :=
    boundedWaste_eventually_normal C waste start sigma hbounded
  have hwithin : cumulativeAttemptCost N ≤ budget := by
    simpa [cumulativeAttemptCost] using le_trans hN hbudget
  have hconsistent : Consistent C (adaptiveRun N sigma start) :=
    (hcomplete (adaptiveRun N sigma start)).mp hnormal
  let w : FixedPointObject (presentation C seed) (localSystem C seed) :=
    ⟨(presentation C seed).toPublicWorld (adaptiveRun N sigma start),
      ⟨⟨adaptiveRun N sigma start,
        normalForm_isFixed C seed hnormal⟩, rfl⟩⟩
  exact ⟨N, hwithin, hconsistent, hconstant, w, rfl⟩

theorem boundedWaste_public_endpoint_exists_unique
    (seed start : Records C) (waste budget : Nat)
    (hcomplete : Completeness C)
    (hconf : Confluence C)
    (hbudget : (waste + 1) * mismatchCount start ≤ budget) :
    ∃ w : FixedPointObject (presentation C seed) (localSystem C seed),
      ∀ sigma : AdaptiveScheduler C,
        BoundedWaste C waste sigma →
        ∃ N, cumulativeAttemptCost N ≤ budget ∧
          Consistent C (adaptiveRun N sigma start) ∧
          (∀ n, N ≤ n →
            adaptiveRun n sigma start = adaptiveRun N sigma start) ∧
          (presentation C seed).toPublicWorld
              (adaptiveRun N sigma start) = w.1 := by
  obtain ⟨w, hw⟩ := canonicalEndpoint_mem_fixedPointObject C seed start
  refine ⟨w, ?_⟩
  intro sigma hbounded
  obtain ⟨N, hN, hnormal, hconstant⟩ :=
    boundedWaste_eventually_normal C waste start sigma hbounded
  have hwithin : cumulativeAttemptCost N ≤ budget := by
    simpa [cumulativeAttemptCost] using le_trans hN hbudget
  have hconsistent : Consistent C (adaptiveRun N sigma start) :=
    (hcomplete (adaptiveRun N sigma start)).mp hnormal
  have hadaptive : Relation.ReflTransGen (acceptedStep C) start
      (adaptiveRun N sigma start) :=
    adaptiveRun_reachable C N sigma start
  have hcanonical : Relation.ReflTransGen (acceptedStep C) start
      (Repair C start) := Repair_reachable C start
  have heq : adaptiveRun N sigma start = Repair C start :=
    AbstractRewriting.unique_normal_form (acceptedStep C) hconf
      ⟨hadaptive, hnormal⟩
      ⟨hcanonical, Repair_normalForm C start⟩
  refine ⟨N, hwithin, hconsistent, hconstant, ?_⟩
  calc
    (presentation C seed).toPublicWorld (adaptiveRun N sigma start) =
        (presentation C seed).toPublicWorld (Repair C start) :=
      congrArg (presentation C seed).toPublicWorld heq
    _ = w.1 := hw.symm

#print axioms boundedWaste_endpoint_within_attempt_budget
#print axioms boundedWaste_public_endpoint_exists_unique

end

end OPH.Tower.CumulativeCapacityEndpoint
