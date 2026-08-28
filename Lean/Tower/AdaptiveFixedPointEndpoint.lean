import Tower.FixedPointEndpoint
import ObserverPatchHolography.Execution.AdaptiveRunStratification

/-!
# Adaptive canonical-repair endpoints

This module is a new consumer of the source `adaptiveRun` stratification.
It does not rewire the existing finite `CompletedSchedule` results. A
pathwise weak-fair adaptive scheduler reaches a stable normal form. Under
the separately named completeness and confluence premises, all such runs
from one start reach the canonical `Repair` endpoint and hence the same
public fixed-point object.

No finite stochastic kernel, active-source law, almost-sure first-hit
claim, lumpability, physical scheduler, or private/shared provenance is
produced here.
-/

namespace OPH.Tower.AdaptiveEndpoint

open OPH
open OPH.Execution
open OPH.Locality
open OPH.Tower.PublicWorldPresentation
open OPH.Tower.PublicWorldPresentation.FiniteRepairSystem
open OPH.Tower.PublicWorldPresentation.OPHPrimitiveEndpoint

noncomputable section

variable (C : OPHCarrier)
variable [FinitePatchStates C]

theorem normalForm_isFixed (seed : Records C) {x : Records C}
    (hnormal : NormalForm C x) :
    IsFixed (presentation C seed) (localSystem C seed) x := by
  intro i
  by_contra hfire
  exact hnormal (localRepair C i x) ⟨i, rfl, hfire⟩

theorem adaptive_endpoint_mem_fixedPointObject
    (seed start : Records C) (sigma : AdaptiveScheduler C)
    (hcomplete : Completeness C)
    (hfair : PathwiseWeakFair C sigma start) :
    ∃ N,
      Consistent C (adaptiveRun N sigma start) ∧
      (∀ n, N ≤ n →
        adaptiveRun n sigma start = adaptiveRun N sigma start) ∧
      ∃ w : FixedPointObject (presentation C seed) (localSystem C seed),
        w.1 = (presentation C seed).toPublicWorld
          (adaptiveRun N sigma start) := by
  obtain ⟨N, hnormal, hconstant⟩ :=
    pathwiseWeakFair_eventual_normalForm C sigma start hfair
  have hconsistent : Consistent C (adaptiveRun N sigma start) :=
    (hcomplete (adaptiveRun N sigma start)).mp hnormal
  let w : FixedPointObject (presentation C seed) (localSystem C seed) :=
    ⟨(presentation C seed).toPublicWorld (adaptiveRun N sigma start),
      ⟨⟨adaptiveRun N sigma start,
        normalForm_isFixed C seed hnormal⟩, rfl⟩⟩
  exact ⟨N, hconsistent, hconstant, w, rfl⟩

theorem adaptive_public_endpoint_exists_unique
    (seed start : Records C)
    (hcomplete : Completeness C)
    (hconf : Confluence C) :
    ∃ w : FixedPointObject (presentation C seed) (localSystem C seed),
      ∀ sigma : AdaptiveScheduler C,
        PathwiseWeakFair C sigma start →
        ∃ N,
          Consistent C (adaptiveRun N sigma start) ∧
          (∀ n, N ≤ n →
            adaptiveRun n sigma start = adaptiveRun N sigma start) ∧
          (presentation C seed).toPublicWorld
              (adaptiveRun N sigma start) = w.1 := by
  obtain ⟨w, hw⟩ := canonicalEndpoint_mem_fixedPointObject C seed start
  refine ⟨w, ?_⟩
  intro sigma hfair
  obtain ⟨N, hnormal, hconstant⟩ :=
    pathwiseWeakFair_eventual_normalForm C sigma start hfair
  have hconsistent : Consistent C (adaptiveRun N sigma start) :=
    (hcomplete (adaptiveRun N sigma start)).mp hnormal
  have hadaptive : Relation.ReflTransGen (acceptedStep C) start
      (adaptiveRun N sigma start) :=
    adaptiveRun_reachable C N sigma start
  have hcanonical : Relation.ReflTransGen (acceptedStep C) start
      (Repair C start) :=
    Repair_reachable C start
  have heq : adaptiveRun N sigma start = Repair C start :=
    AbstractRewriting.unique_normal_form (acceptedStep C) hconf
      ⟨hadaptive, hnormal⟩
      ⟨hcanonical, Repair_normalForm C start⟩
  refine ⟨N, hconsistent, hconstant, ?_⟩
  calc
    (presentation C seed).toPublicWorld (adaptiveRun N sigma start) =
        (presentation C seed).toPublicWorld (Repair C start) :=
      congrArg (presentation C seed).toPublicWorld heq
    _ = w.1 := hw.symm

#print axioms adaptive_endpoint_mem_fixedPointObject
#print axioms adaptive_public_endpoint_exists_unique

end

end OPH.Tower.AdaptiveEndpoint
