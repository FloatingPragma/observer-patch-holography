import Mathlib.NumberTheory.Real.Irrational
import Mathlib.Algebra.BigOperators.Ring.Finset

/-!
# Native finite-tape references and an exact sampling boundary

For any finite deterministic record map, retain every complete-tape
multiplicity, including unused draws. This applies to bounded adaptive
native sharp-register programs after their fair-read reduction.

Finite conditioning preserves rational masses. The native three-record
reference has a mean-one entropy projection with mass at outcome zero
(4-sqrt 2)/7. Bounded fair-read programs therefore cannot sample this
optimizer exactly. The experiment's integer threshold controller gives a
bounded approximation and an exact sampler with unbounded stopping time.
-/

namespace OPH.InformationProjection.NativePath

open scoped BigOperators

variable {T X Y : Type*} [Fintype T] [Fintype X] [Fintype Y]
  [DecidableEq X] [DecidableEq Y]

/-- Pushforward with all tape multiplicities retained. -/
def pushQ (f : T → X) (p : T → ℚ) (x : X) : ℚ :=
  ∑ t, if f t = x then p t else 0

/-- The record law of a uniform complete tape. -/
def tapeReference (f : T → X) (x : X) : ℚ :=
  pushQ f (fun _ => (Fintype.card T : ℚ)⁻¹) x

omit [Fintype Y] in
/-- All coarse references are pushforwards of the fine reference. -/
theorem pushQ_comp (f : T → X) (g : X → Y) (p : T → ℚ) :
    pushQ g (pushQ f p) = pushQ (g ∘ f) p := by
  funext y
  unfold pushQ
  calc
    ∑ x, (if g x = y then ∑ t, if f t = x then p t else 0 else 0)
      = ∑ x, ∑ t, if f t = x then (if g x = y then p t else 0) else 0 := by
        apply Finset.sum_congr rfl
        intro x _
        by_cases h : g x = y <;> simp [h]
    _ = ∑ t, ∑ x, if f t = x then (if g x = y then p t else 0) else 0 :=
      Finset.sum_comm
    _ = ∑ t, if (g ∘ f) t = y then p t else 0 := by
      apply Finset.sum_congr rfl
      intro t _
      simp [eq_comm]

/-- Normalization is inherited from the complete tape law. -/
theorem pushQ_total (f : T → X) (p : T → ℚ) :
    ∑ x, pushQ f p x = ∑ t, p t := by
  unfold pushQ
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro t _
  simp [eq_comm]

theorem tapeReference_total [Nonempty T] (f : T → X) :
    ∑ x, tapeReference f x = 1 := by
  unfold tapeReference
  rw [pushQ_total]
  simp [Fintype.card_ne_zero]

omit [Fintype X] in
/-- Each record receives one equal contribution per complete tape in its fibre. -/
theorem tapeReference_fibre (f : T → X) (x : X) :
    tapeReference f x = ((Finset.univ.filter (fun t => f t = x)).card : ℚ) / Fintype.card T := by
  simp [tapeReference, pushQ, ← Finset.sum_filter, div_eq_mul_inv]

omit [Fintype X] in
/-- Padding each tape with the same finite set of unused suffixes preserves its record law. -/
theorem tapeReference_padding {B : Type*} [Fintype B] [Nonempty B] (f : T → X) :
    tapeReference (fun tb : T × B => f tb.1) = tapeReference f := by
  funext x
  unfold tapeReference pushQ
  rw [Fintype.sum_prod_type]
  apply Finset.sum_congr rfl
  intro t _
  by_cases h : f t = x
  · simp [h, Fintype.card_prod, mul_inv_rev, Fintype.card_ne_zero]
  · simp [h]

omit [Fintype X] in
/-- Selecting an event visible to a reader commutes with pushing to that reader. -/
theorem visible_condition_push (f : T → X) (p : T → ℚ) (s : Finset X) (mass : ℚ) :
    pushQ f (fun t => if f t ∈ s then p t / mass else 0) =
      fun x => if x ∈ s then pushQ f p x / mass else 0 := by
  funext x
  by_cases hx : x ∈ s
  · simp only [hx, if_true, pushQ, div_eq_mul_inv, Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro t _
    by_cases h : f t = x <;> simp [h, hx]
  · simp only [hx, if_false, pushQ]
    apply Finset.sum_eq_zero
    intro t _
    by_cases h : f t = x <;> simp [h, hx]

/-- Any event after finite postselection has a rational mass. -/
def conditionedEventMass [DecidableEq T] (accepted event : Finset T) : ℚ :=
  ((accepted ∩ event).card : ℚ) / accepted.card

theorem optimum_zero_irrational : Irrational ((4 - Real.sqrt 2) / 7) := by
  simpa using (irrational_sqrt_two.ratCast_sub (4 : ℚ)).div_ratCast
    (by norm_num : (7 : ℚ) ≠ 0)

omit [Fintype T] in
/-- A universal bounded fair-read exact A3 sampler is impossible. -/
theorem no_finite_conditioned_sampler [DecidableEq T]
    (accepted event : Finset T) :
    (conditionedEventMass accepted event : ℝ) ≠ (4 - Real.sqrt 2) / 7 := by
  exact (optimum_zero_irrational.ne_rat _).symm

/-- A cover containing the full record separates rational classical joint laws. -/
theorem full_record_separates (p q : X → ℚ) (h : pushQ id p = pushQ id q) : p = q := by
  funext x
  simpa [pushQ, eq_comm] using congrFun h x

end OPH.InformationProjection.NativePath

#print axioms OPH.InformationProjection.NativePath.pushQ_comp
#print axioms OPH.InformationProjection.NativePath.tapeReference_total
#print axioms OPH.InformationProjection.NativePath.tapeReference_fibre
#print axioms OPH.InformationProjection.NativePath.tapeReference_padding
#print axioms OPH.InformationProjection.NativePath.visible_condition_push
#print axioms OPH.InformationProjection.NativePath.no_finite_conditioned_sampler
#print axioms OPH.InformationProjection.NativePath.full_record_separates
