import SeamCurrentCarrierQuotient

open scoped BigOperators

namespace OPH.RegionalContinuity

open OPH.PrimitivePortFrameQuotient
open OPH.RepairWordCarrierReadout
open OPH.SeamCurrentCarrierQuotient

/-!
# Regional continuity on the exact twelve-port seam graph

`seamBoundary` uses the incidence convention `-n` at the smaller endpoint
and `+n` at the larger endpoint.  Thus it records net inflow.  The outward
flux below uses the opposite regional sign, so the finite balance law is

`qNext - q = source + seamBoundary current`

and regional summation gives

`change in stored load = regional source - outward flux`.

The statements are exact finite incidence identities.  They do not identify
the integer load with a physical charge, the seam current with a stress
current, or the finite update index with time.
-/

/-- Stored integer load in a finite set of ports. -/
def regionLoad (region : Finset (Fin 12)) (q : PortLoad) : ℤ :=
  ∑ p ∈ region, q p

/-- Outward flux carried by one canonically oriented seam.  Positive current
points from `seamLeft e` to `seamRight e`. -/
def seamOutwardFlux (region : Finset (Fin 12)) (e : Fin 30) (n : ℤ) : ℤ :=
  (if seamLeft e ∈ region then n else 0) -
    (if seamRight e ∈ region then n else 0)

/-- Total signed outward flux through the boundary of a port region.  Seams
with both endpoints inside or both endpoints outside contribute zero. -/
def regionFlux (region : Finset (Fin 12)) (current : SeamCurrent) : ℤ :=
  ∑ e : Fin 30, seamOutwardFlux region e (current e)

/-- Regional summation of one directed incidence boundary is the negative of
its outward flux. -/
theorem directedBoundary_region_eq_neg_outwardFlux
    (region : Finset (Fin 12)) (u v : Fin 12) (n : ℤ) :
    (∑ p ∈ region, directedBoundary u v n p) =
      -((if u ∈ region then n else 0) -
        (if v ∈ region then n else 0)) := by
  classical
  simp [directedBoundary]

/-- The incidence sum over a region is minus its signed outward boundary
flux.  Internal seam contributions cancel exactly. -/
theorem seamBoundary_region_eq_neg_flux
    (region : Finset (Fin 12)) (current : SeamCurrent) :
    regionLoad region (seamBoundary current) = -regionFlux region current := by
  classical
  unfold regionLoad seamBoundary regionFlux seamOutwardFlux
  calc
    (∑ p ∈ region, ∑ e : Fin 30,
        directedBoundary (seamLeft e) (seamRight e) (current e) p) =
        ∑ e : Fin 30, ∑ p ∈ region,
          directedBoundary (seamLeft e) (seamRight e) (current e) p := by
      rw [Finset.sum_comm]
    _ = ∑ e : Fin 30,
        -((if seamLeft e ∈ region then current e else 0) -
          (if seamRight e ∈ region then current e else 0)) := by
      apply Finset.sum_congr rfl
      intro e _
      exact directedBoundary_region_eq_neg_outwardFlux
        region (seamLeft e) (seamRight e) (current e)
    _ = -∑ e : Fin 30,
        ((if seamLeft e ∈ region then current e else 0) -
          (if seamRight e ∈ region then current e else 0)) := by
      rw [Finset.sum_neg_distrib]

/-- No seam crosses the boundary of the complete twelve-port region. -/
theorem regionFlux_univ (current : SeamCurrent) :
    regionFlux Finset.univ current = 0 := by
  simp [regionFlux, seamOutwardFlux]

/-- One finite continuity update in the incidence convention of
`SeamCurrentCarrierQuotient`: boundary load is net inflow. -/
def SatisfiesContinuityUpdate
    (q qNext source : PortLoad) (current : SeamCurrent) : Prop :=
  qNext - q = source + seamBoundary current

/-- Typed finite object consumed by a later Ward-limit manifest. It contains
the actual B5 pointwise continuity premise rather than a detached source
label; nonzero activity is imposed separately by the promotion manifest. -/
structure FiniteContinuityWitness where
  q : PortLoad
  qNext : PortLoad
  source : PortLoad
  current : SeamCurrent
  satisfiesUpdate : SatisfiesContinuityUpdate q qNext source current

/-- Exact regional continuity: stored-load change equals the regional source
minus signed outward seam flux. -/
theorem regional_continuity
    (region : Finset (Fin 12)) (q qNext source : PortLoad)
    (current : SeamCurrent)
    (hupdate : SatisfiesContinuityUpdate q qNext source current) :
    regionLoad region qNext - regionLoad region q =
      regionLoad region source - regionFlux region current := by
  classical
  unfold SatisfiesContinuityUpdate at hupdate
  have hpoint : ∀ p : Fin 12,
      qNext p - q p = source p + seamBoundary current p := by
    intro p
    have hp := congrFun hupdate p
    simpa using hp
  calc
    regionLoad region qNext - regionLoad region q =
        ∑ p ∈ region, (qNext p - q p) := by
      simp [regionLoad, Finset.sum_sub_distrib]
    _ = ∑ p ∈ region, (source p + seamBoundary current p) := by
      apply Finset.sum_congr rfl
      intro p _
      exact hpoint p
    _ = regionLoad region source + regionLoad region (seamBoundary current) := by
      simp [regionLoad, Finset.sum_add_distrib]
    _ = regionLoad region source - regionFlux region current := by
      rw [seamBoundary_region_eq_neg_flux]
      ring

/-- Global continuity on the closed seam graph: total-load change is exactly
the total source. -/
theorem global_continuity
    (q qNext source : PortLoad) (current : SeamCurrent)
    (hupdate : SatisfiesContinuityUpdate q qNext source current) :
    regionLoad Finset.univ qNext - regionLoad Finset.univ q =
      regionLoad Finset.univ source := by
  rw [regional_continuity Finset.univ q qNext source current hupdate,
    regionFlux_univ, sub_zero]

/-- A source with zero total load preserves the closed graph's total stored
load. -/
theorem global_conservation_of_zero_total_source
    (q qNext source : PortLoad) (current : SeamCurrent)
    (hupdate : SatisfiesContinuityUpdate q qNext source current)
    (hsource : regionLoad Finset.univ source = 0) :
    regionLoad Finset.univ qNext = regionLoad Finset.univ q := by
  have h := global_continuity q qNext source current hupdate
  rw [hsource] at h
  exact sub_eq_zero.mp h

/-- Every finite continuity witness carries the exact regional balance law. -/
theorem FiniteContinuityWitness.regionalBalance
    (witness : FiniteContinuityWitness) (region : Finset (Fin 12)) :
    regionLoad region witness.qNext - regionLoad region witness.q =
      regionLoad region witness.source - regionFlux region witness.current :=
  regional_continuity region witness.q witness.qNext witness.source
    witness.current witness.satisfiesUpdate

/-! ## Axiom audit -/

#print axioms directedBoundary_region_eq_neg_outwardFlux
#print axioms seamBoundary_region_eq_neg_flux
#print axioms regionFlux_univ
#print axioms regional_continuity
#print axioms global_continuity
#print axioms global_conservation_of_zero_total_source
#print axioms FiniteContinuityWitness.regionalBalance

end OPH.RegionalContinuity
