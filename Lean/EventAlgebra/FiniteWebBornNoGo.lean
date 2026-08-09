import EventAlgebra.InterlockingContexts

/-!
# A planar cubic counterexample on the current finite context web

The finite unsharp battery in `InterlockingContexts` uses only directions in
the `x-z` plane.  Consequently it excludes the committed `z`-cubic response
but leaves the transverse `y`-cubic response invisible on every added
unsharp slot.  This module packages that response as a normalized,
probability-valued, noncontextual assignment on the entire current
`extendedWeb` and proves that its sharp binary restriction is not affine.
Therefore the present `FiniteBuschGleasonInterface` is false.

Boundary.  This is a finite-web incompleteness theorem, not a no-go for
Busch's theorem on all effects.  Adding an unsharp battery outside the common
plane, or assuming an additive valuation on the full effect algebra, can
remove this counterexample.
-/

namespace EventAlgebra.InterlockingContexts

open OPH.C1Lorentz
open EventAlgebra.FiniteEffectClosureBoundary

noncomputable section

/-! ## The transverse response -/

/-- The cubic binary response in the coordinate transverse to every axis in
the current finite unsharp battery. -/
def planarCubicWeight (n : Spatial) : ℝ :=
  (1 + (n 1) ^ 3) / 2

/-- A single effect response.  Sharp projectors have scalar coordinate
`1/2`, in which case the spatial coordinate recovers their direction; all
unsharp effects in the finite battery are valued by their scalar coordinate.
-/
def planarCubicEffectResponse (E : Herm2) : ℝ :=
  if E.1 = 1 / 2 then (1 + (2 * E.2 1) ^ 3) / 2 else E.1

theorem planarCubicWeight_neg (n : Spatial) :
    planarCubicWeight (-n) = 1 - planarCubicWeight n := by
  simp only [planarCubicWeight, Pi.neg_apply]
  ring

theorem planarCubicWeight_antipodal_sum (n : Spatial) :
    planarCubicWeight n + planarCubicWeight (-n) = 1 := by
  rw [planarCubicWeight_neg]
  ring

/-- On a sharp projector the global effect response is exactly the
transverse cubic binary response. -/
@[simp] theorem planarCubicEffectResponse_proj (n : Spatial) :
    planarCubicEffectResponse (proj n) = planarCubicWeight n := by
  simp [planarCubicEffectResponse, planarCubicWeight, proj, scaledProj,
    Pi.smul_apply, smul_eq_mul]

@[simp] theorem planarCubicEffectResponse_grain (n : Spatial) :
    planarCubicEffectResponse (grain n) = 1 / 16 := by
  norm_num [planarCubicEffectResponse, grain, scaledProj]

@[simp] theorem planarCubicEffectResponse_scaled_three_quarters (n : Spatial) :
    planarCubicEffectResponse (scaledProj (3 / 4) n) = 3 / 8 := by
  norm_num [planarCubicEffectResponse, scaledProj]

@[simp] theorem planarCubicEffectResponse_scaled_five_eighths (n : Spatial) :
    planarCubicEffectResponse (scaledProj (5 / 8) n) = 5 / 16 := by
  norm_num [planarCubicEffectResponse, scaledProj]

@[simp] theorem planarCubicWeight_u1 : planarCubicWeight u1 = 1 / 2 := by
  norm_num [planarCubicWeight, u1]

@[simp] theorem planarCubicWeight_u2 : planarCubicWeight u2 = 1 / 2 := by
  norm_num [planarCubicWeight, u2]

@[simp] theorem planarCubicWeight_u3 : planarCubicWeight u3 = 1 / 2 := by
  norm_num [planarCubicWeight, u3]

@[simp] theorem planarCubicWeight_neg_u1 : planarCubicWeight (-u1) = 1 / 2 := by
  rw [planarCubicWeight_neg, planarCubicWeight_u1]
  norm_num

@[simp] theorem planarCubicWeight_neg_u2 : planarCubicWeight (-u2) = 1 / 2 := by
  rw [planarCubicWeight_neg, planarCubicWeight_u2]
  norm_num

@[simp] theorem planarCubicWeight_neg_u3 : planarCubicWeight (-u3) = 1 / 2 := by
  rw [planarCubicWeight_neg, planarCubicWeight_u3]
  norm_num

theorem planarCubicEffectResponse_bin_sum (n : Spatial)
    (hn : spatialNormSq n = 1) :
    ∑ j, planarCubicEffectResponse ((binContext n hn).member j) = 1 := by
  unfold binContext
  rw [Fin.sum_univ_two]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one,
    planarCubicEffectResponse_proj]
  exact planarCubicWeight_antipodal_sum n

theorem planarCubicEffectResponse_trine_sum :
    ∑ j, planarCubicEffectResponse (trineContext.member j) = 1 := by
  unfold trineContext
  norm_num [trineE1, trineE2, trineE3, Fin.sum_univ_succ]

theorem planarCubicEffectResponse_rep_u1_sum :
    ∑ j, planarCubicEffectResponse ((repContext u1 u1_unit).member j) = 1 := by
  unfold repContext
  norm_num [Fin.sum_univ_succ]

theorem planarCubicEffectResponse_rep_u2_sum :
    ∑ j, planarCubicEffectResponse ((repContext u2 u2_unit).member j) = 1 := by
  unfold repContext
  norm_num [Fin.sum_univ_succ]

theorem planarCubicEffectResponse_rep_u3_sum :
    ∑ j, planarCubicEffectResponse ((repContext u3 u3_unit).member j) = 1 := by
  unfold repContext
  norm_num [Fin.sum_univ_succ]

theorem planarCubicEffectResponse_splitA_u1_sum :
    ∑ j, planarCubicEffectResponse ((splitContextA u1 u1_unit).member j) = 1 := by
  unfold splitContextA
  norm_num [Fin.sum_univ_succ]

theorem planarCubicEffectResponse_splitB_u2_sum :
    ∑ j, planarCubicEffectResponse ((splitContextB u2 u2_unit).member j) = 1 := by
  unfold splitContextB
  norm_num [Fin.sum_univ_succ]

theorem planarCubicEffectResponse_splitB_u3_sum :
    ∑ j, planarCubicEffectResponse ((splitContextB u3 u3_unit).member j) = 1 := by
  unfold splitContextB
  norm_num [Fin.sum_univ_succ]

/-! ## A normalized noncontextual assignment on the full finite web -/

/-- The transverse cubic response, evaluated slotwise on the current
extended web. -/
def planarCubicAssignment : Assignment extendedWeb where
  value := fun C j => planarCubicEffectResponse (C.member j)
  normalized := by
    intro C hC
    rcases hC with hBinary | hFinite
    · rcases hBinary with ⟨n, hn, rfl⟩
      exact planarCubicEffectResponse_bin_sum n hn
    · simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hFinite
      rcases hFinite with hC | hC | hC | hC | hC | hC | hC
      · subst C
        exact planarCubicEffectResponse_trine_sum
      · subst C
        exact planarCubicEffectResponse_rep_u1_sum
      · subst C
        exact planarCubicEffectResponse_rep_u2_sum
      · subst C
        exact planarCubicEffectResponse_rep_u3_sum
      · subst C
        exact planarCubicEffectResponse_splitA_u1_sum
      · subst C
        exact planarCubicEffectResponse_splitB_u2_sum
      · subst C
        exact planarCubicEffectResponse_splitB_u3_sum

/-- The assignment factors through one effect response and is therefore
noncontextual on all interlocking slots. -/
theorem planarCubicAssignment_noncontextual :
    Noncontextual extendedWeb planarCubicAssignment :=
  (noncontextual_iff_factorsThroughEffects extendedWeb
    planarCubicAssignment).mpr
    ⟨planarCubicEffectResponse, fun _ _ _ => rfl⟩

/-- The transverse cubic is probability-valued on every unit direction. -/
theorem planarCubicWeight_mem_Icc {n : Spatial}
    (hn : spatialNormSq n = 1) : planarCubicWeight n ∈ Set.Icc (0 : ℝ) 1 := by
  have hy := coordinate_mem_unitInterval hn (1 : Fin 3)
  have hquadLower : 0 ≤ (n 1) ^ 2 - n 1 + 1 := by
    nlinarith [sq_nonneg (n 1 - (1 / 2 : ℝ))]
  have hquadUpper : 0 ≤ (n 1) ^ 2 + n 1 + 1 := by
    nlinarith [sq_nonneg (n 1 + (1 / 2 : ℝ))]
  have hcubedLower : 0 ≤ (n 1 + 1) * ((n 1) ^ 2 - n 1 + 1) :=
    mul_nonneg (by linarith [hy.1]) hquadLower
  have hcubedUpper : 0 ≤ (1 - n 1) * ((n 1) ^ 2 + n 1 + 1) :=
    mul_nonneg (by linarith [hy.2]) hquadUpper
  constructor <;>
    simp only [planarCubicWeight] <;>
    nlinarith [hcubedLower, hcubedUpper]

theorem planarCubicEffectResponse_bin_mem_Icc (n : Spatial)
    (hn : spatialNormSq n = 1) (j : Fin 2) :
    planarCubicEffectResponse ((binContext n hn).member j) ∈
      Set.Icc (0 : ℝ) 1 := by
  unfold binContext
  fin_cases j
  · simpa using planarCubicWeight_mem_Icc hn
  · simpa using planarCubicWeight_mem_Icc
      (by rw [spatialNormSq_neg, hn])

theorem planarCubicEffectResponse_trine_mem_Icc (j : Fin 3) :
    planarCubicEffectResponse (trineContext.member j) ∈ Set.Icc (0 : ℝ) 1 := by
  unfold trineContext
  fin_cases j <;>
    norm_num [trineE1, trineE2, trineE3]

theorem planarCubicEffectResponse_rep_u1_mem_Icc (j : Fin 9) :
    planarCubicEffectResponse ((repContext u1 u1_unit).member j) ∈
      Set.Icc (0 : ℝ) 1 := by
  unfold repContext
  fin_cases j <;> norm_num

theorem planarCubicEffectResponse_rep_u2_mem_Icc (j : Fin 9) :
    planarCubicEffectResponse ((repContext u2 u2_unit).member j) ∈
      Set.Icc (0 : ℝ) 1 := by
  unfold repContext
  fin_cases j <;> norm_num

theorem planarCubicEffectResponse_rep_u3_mem_Icc (j : Fin 9) :
    planarCubicEffectResponse ((repContext u3 u3_unit).member j) ∈
      Set.Icc (0 : ℝ) 1 := by
  unfold repContext
  fin_cases j <;> norm_num

theorem planarCubicEffectResponse_splitA_u1_mem_Icc (j : Fin 4) :
    planarCubicEffectResponse ((splitContextA u1 u1_unit).member j) ∈
      Set.Icc (0 : ℝ) 1 := by
  unfold splitContextA
  fin_cases j <;> norm_num

theorem planarCubicEffectResponse_splitB_u2_mem_Icc (j : Fin 5) :
    planarCubicEffectResponse ((splitContextB u2 u2_unit).member j) ∈
      Set.Icc (0 : ℝ) 1 := by
  unfold splitContextB
  fin_cases j <;> norm_num

theorem planarCubicEffectResponse_splitB_u3_mem_Icc (j : Fin 5) :
    planarCubicEffectResponse ((splitContextB u3 u3_unit).member j) ∈
      Set.Icc (0 : ℝ) 1 := by
  unfold splitContextB
  fin_cases j <;> norm_num

/-- Every slot of the counterexample on the current extended web lies in
the probability interval. -/
theorem planarCubicAssignment_value_mem_Icc
    (C : Context) (hC : C ∈ extendedWeb) (j : Fin C.size) :
    planarCubicAssignment.value C j ∈ Set.Icc (0 : ℝ) 1 := by
  rcases hC with hBinary | hFinite
  · rcases hBinary with ⟨n, hn, rfl⟩
    exact planarCubicEffectResponse_bin_mem_Icc n hn j
  · simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hFinite
    rcases hFinite with hC | hC | hC | hC | hC | hC | hC
    · subst C
      exact planarCubicEffectResponse_trine_mem_Icc j
    · subst C
      exact planarCubicEffectResponse_rep_u1_mem_Icc j
    · subst C
      exact planarCubicEffectResponse_rep_u2_mem_Icc j
    · subst C
      exact planarCubicEffectResponse_rep_u3_mem_Icc j
    · subst C
      exact planarCubicEffectResponse_splitA_u1_mem_Icc j
    · subst C
      exact planarCubicEffectResponse_splitB_u2_mem_Icc j
    · subst C
      exact planarCubicEffectResponse_splitB_u3_mem_Icc j

/-! ## Exact non-affinity and interface failure -/

/-- The transverse cubic response is not an affine/Born binary functional.
The three coordinate axes fix `q = (0,1,0)` and the exact unit direction
`(0,4/5,3/5)` then distinguishes `y^3` from `y`. -/
theorem planarCubicWeight_not_affine :
    ¬ ∃ q : Spatial, ∀ n : Spatial, spatialNormSq n = 1 →
      planarCubicWeight n = affineBinaryWeight q n := by
  rintro ⟨q, hq⟩
  have hx := hq (![(1 : ℝ), 0, 0] : Spatial) (by
    norm_num [spatialNormSq, Fin.sum_univ_succ])
  have hy := hq (![(0 : ℝ), 1, 0] : Spatial) (by
    norm_num [spatialNormSq, Fin.sum_univ_succ])
  have hz := hq (![(0 : ℝ), 0, 1] : Spatial) (by
    norm_num [spatialNormSq, Fin.sum_univ_succ])
  have hq0 : q 0 = 0 := by
    simp [planarCubicWeight, affineBinaryWeight, spatialDot,
      Fin.sum_univ_succ] at hx
    linarith
  have hq1 : q 1 = 1 := by
    simp [planarCubicWeight, affineBinaryWeight, spatialDot,
      Fin.sum_univ_succ] at hy
    linarith
  have hq2 : q 2 = 0 := by
    simp [planarCubicWeight, affineBinaryWeight, spatialDot,
      Fin.sum_univ_succ] at hz
    linarith
  have hw := hq (![(0 : ℝ), 4 / 5, 3 / 5] : Spatial) (by
    norm_num [spatialNormSq, Fin.sum_univ_succ])
  simp [planarCubicWeight, affineBinaryWeight, spatialDot,
    Fin.sum_univ_succ, hq0, hq1, hq2] at hw
  norm_num at hw

/-- The consumption interface currently declared in
`InterlockingContexts` is false. -/
theorem finiteBuschGleasonInterface_false :
    ¬ FiniteBuschGleasonInterface := by
  intro hInterface
  obtain ⟨q, _, hq⟩ :=
    hInterface planarCubicAssignment planarCubicAssignment_noncontextual
  apply planarCubicWeight_not_affine
  refine ⟨q, fun n hn => ?_⟩
  have h := hq n hn
  simpa [planarCubicAssignment, binContext] using h

/-- Summary: a probability-valued, normalized, noncontextual assignment on
the current finite web has no affine binary representation. -/
theorem current_finite_web_born_no_go :
    ∃ v : Assignment extendedWeb,
      Noncontextual extendedWeb v ∧
      (∀ C, C ∈ extendedWeb → ∀ j, v.value C j ∈ Set.Icc (0 : ℝ) 1) ∧
      ¬ ∃ q : Spatial, ∀ (n : Spatial) (hn : spatialNormSq n = 1),
        v.value (binContext n hn) (0 : Fin 2) = affineBinaryWeight q n := by
  refine ⟨planarCubicAssignment, planarCubicAssignment_noncontextual,
    planarCubicAssignment_value_mem_Icc, ?_⟩
  rintro ⟨q, hq⟩
  apply planarCubicWeight_not_affine
  exact ⟨q, fun n hn => by
    simpa [planarCubicAssignment, binContext] using hq n hn⟩

#print axioms planarCubicAssignment_noncontextual
#print axioms planarCubicAssignment_value_mem_Icc
#print axioms planarCubicWeight_not_affine
#print axioms finiteBuschGleasonInterface_false
#print axioms current_finite_web_born_no_go

end


end EventAlgebra.InterlockingContexts
