import Mathlib

/-!
# Register-Clebsch light-family ratio boundary

This module formalizes the algebraic part of the rejected down-type
register-Clebsch lane.  It does not formalize the empirical FLAG comparison
or promote a quark-mass prediction.

The first theorem states the exact boundary identity implied by

* `y_s = y_mu / 3`, and
* `y_d = 3 * y_e`.

The second theorem states the limited transport result used by the numerical
lane: one common nonzero multiplicative transport factor cancels from the
light-family ratio.  Flavor-dependent threshold or matching factors are
outside this theorem and would define a modified lane.

The final theorem evaluates all six assignments of the conditionally selected,
runtime-target-free unordered weight set `{1/3, 1, 3}`.  The list contains
every permutation once.  Its light-family factors are therefore restricted to
`1/9`, `1/3`, `3`, and `9`.
-/

namespace OPH.ClebschRatio

/-- The exact light-family identity at the register boundary. -/
theorem boundary_light_ratio
    {K : Type*} [Field K]
    {yMu yE yS yD : K}
    (hyE : yE ≠ 0)
    (hyS : yS = yMu / 3)
    (hyD : yD = 3 * yE) :
    yS / yD = (yMu / yE) / 9 := by
  rw [hyS, hyD]
  field_simp
  ring

/-- A shared nonzero multiplicative transport preserves the boundary ratio. -/
theorem common_transport_light_ratio
    {K : Type*} [Field K]
    {yMu yE transport : K}
    (hyE : yE ≠ 0)
    (htransport : transport ≠ 0) :
    (transport * (yMu / 3)) / (transport * (3 * yE)) =
      (yMu / yE) / 9 := by
  field_simp
  ring

/-- One assignment of the unordered register weights to `(b, s, d)`. -/
structure Assignment where
  b : ℚ
  s : ℚ
  d : ℚ
  deriving DecidableEq, Repr

/-- Every permutation of `{1/3, 1, 3}`, listed once. -/
def assignments : List Assignment :=
  [
    ⟨1, 1 / 3, 3⟩,
    ⟨3, 1 / 3, 1⟩,
    ⟨3, 1, 1 / 3⟩,
    ⟨1 / 3, 1, 3⟩,
    ⟨1, 3, 1 / 3⟩,
    ⟨1 / 3, 3, 1⟩
  ]

/-- The Clebsch factor multiplying the lepton ratio in `m_s/m_d`. -/
def lightRatioFactor (assignment : Assignment) : ℚ :=
  assignment.s / assignment.d

/-- The exact light-ratio factors for all six generation assignments. -/
theorem assignment_light_ratio_factors :
    assignments.map lightRatioFactor =
      [1 / 9, 1 / 3, 3, 1 / 3, 9, 3] := by
  norm_num [assignments, lightRatioFactor]

/-- Every declared row is a permutation of the three distinct weights. -/
theorem every_assignment_is_weight_permutation :
    ∀ assignment ∈ assignments,
      [assignment.b, assignment.s, assignment.d].Perm [1 / 3, 1, 3] := by
  intro assignment h
  simp only [assignments, List.mem_cons, List.not_mem_nil, or_false] at h
  rcases h with rfl | rfl | rfl | rfl | rfl | rfl
  · exact List.Perm.swap _ _ _
  · exact (List.Perm.swap _ _ _).trans
      (List.Perm.cons _ (List.Perm.swap _ _ _))
  · exact (List.Perm.cons _ (List.Perm.swap _ _ _)).trans
      ((List.Perm.swap _ _ _).trans
        (List.Perm.cons _ (List.Perm.swap _ _ _)))
  · exact List.Perm.refl _
  · exact (List.Perm.cons _ (List.Perm.swap _ _ _)).trans
      (List.Perm.swap _ _ _)
  · exact List.Perm.cons _ (List.Perm.swap _ _ _)

/-- No generation assignment occurs twice in the declared table. -/
theorem assignments_nodup : assignments.Nodup := by
  norm_num [assignments]

/-- The assignment table has the cardinality of all permutations of three
distinct weights. Together with the preceding membership and no-duplicate
theorems, this certifies the exhaustive six-row enumeration. -/
theorem assignment_count : assignments.length = 6 := by
  decide

end OPH.ClebschRatio

#print axioms OPH.ClebschRatio.boundary_light_ratio
#print axioms OPH.ClebschRatio.common_transport_light_ratio
#print axioms OPH.ClebschRatio.assignment_light_ratio_factors
#print axioms OPH.ClebschRatio.every_assignment_is_weight_permutation
#print axioms OPH.ClebschRatio.assignments_nodup
#print axioms OPH.ClebschRatio.assignment_count
