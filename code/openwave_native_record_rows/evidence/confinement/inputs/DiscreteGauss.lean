import RegionalContinuity

namespace OPH.DiscreteGauss

open OPH.SeamCurrentCarrierQuotient

/-!
# Exact finite Gauss solvability on the twelve-port seam graph

The existing incidence module supplies an explicit spanning-tree section and
proves that the rational seam-boundary image is the zero-total hyperplane.
This module packages those results as a Gauss interface and identifies the
complete solution ambiguity with the nineteen-dimensional cycle kernel.

No theorem in this file identifies the rational seam current with an
electromagnetic field or the port load with a physical charge density.
-/

/-- The exact incidence range is the subspace of neutral rational port
loads. -/
theorem incidence_range_eq_zeroTotal :
    LinearMap.range rationalSeamBoundary = LinearMap.ker rationalPortTotal :=
  rationalSeamBoundary_range

/-- A rational port load admits a seam-flux solution exactly when its total
load vanishes. -/
theorem gauss_solution_exists_iff_total_zero (charge : RationalPortLoad) :
    (∃ field : RationalSeamCurrent, rationalSeamBoundary field = charge) ↔
      rationalPortTotal charge = 0 := by
  constructor
  · rintro ⟨field, rfl⟩
    have hrange : rationalSeamBoundary field ∈
        LinearMap.range rationalSeamBoundary := ⟨field, rfl⟩
    rw [incidence_range_eq_zeroTotal, LinearMap.mem_ker] at hrange
    exact hrange
  · intro hneutral
    have hker : charge ∈ LinearMap.ker rationalPortTotal := by
      rw [LinearMap.mem_ker]
      exact hneutral
    rw [← incidence_range_eq_zeroTotal] at hker
    exact hker

/-- The existing port-zero spanning-tree section is one explicit declared
Gauss solution for every neutral load. No naturality or canonical-choice
claim is made. -/
theorem rationalBoundarySection_is_gauss_solution
    (charge : RationalPortLoad) (hneutral : rationalPortTotal charge = 0) :
    rationalSeamBoundary (rationalBoundarySection charge) = charge :=
  rationalSeamBoundary_section charge hneutral

/-- Any two flux fields solving the same Gauss equation differ by a
source-free cycle current. -/
theorem gauss_solutions_differ_by_cycle
    {field field' : RationalSeamCurrent} {charge : RationalPortLoad}
    (hfield : rationalSeamBoundary field = charge)
    (hfield' : rationalSeamBoundary field' = charge) :
    field' - field ∈ LinearMap.ker rationalSeamBoundary := by
  rw [LinearMap.mem_ker, LinearMap.map_sub, hfield', hfield, sub_self]

/-- Conversely, translating one solution by a cycle current preserves its
Gauss boundary. -/
theorem add_cycle_preserves_gauss_solution
    {field cycle : RationalSeamCurrent} {charge : RationalPortLoad}
    (hfield : rationalSeamBoundary field = charge)
    (hcycle : cycle ∈ LinearMap.ker rationalSeamBoundary) :
    rationalSeamBoundary (field + cycle) = charge := by
  rw [LinearMap.map_add, hfield, LinearMap.mem_ker.mp hcycle, add_zero]

/-- The full solution fibre over a fixed charge is one affine translate of
the cycle kernel. -/
theorem gauss_solution_iff_difference_is_cycle
    {field field' : RationalSeamCurrent} {charge : RationalPortLoad}
    (hfield : rationalSeamBoundary field = charge) :
    rationalSeamBoundary field' = charge ↔
      field' - field ∈ LinearMap.ker rationalSeamBoundary := by
  constructor
  · intro hfield'
    exact gauss_solutions_differ_by_cycle hfield hfield'
  · intro hcycle
    rw [LinearMap.mem_ker, LinearMap.map_sub, hfield] at hcycle
    exact sub_eq_zero.mp hcycle

/-- The source-free ambiguity has exact rational dimension nineteen. -/
theorem gauss_cycle_space_finrank :
    Module.finrank ℚ (LinearMap.ker rationalSeamBoundary) = 19 :=
  rationalSeamBoundary_kernel_finrank

/-! ## Axiom audit -/

#print axioms incidence_range_eq_zeroTotal
#print axioms gauss_solution_exists_iff_total_zero
#print axioms rationalBoundarySection_is_gauss_solution
#print axioms gauss_solutions_differ_by_cycle
#print axioms add_cycle_preserves_gauss_solution
#print axioms gauss_solution_iff_difference_is_cycle
#print axioms gauss_cycle_space_finrank

end OPH.DiscreteGauss
