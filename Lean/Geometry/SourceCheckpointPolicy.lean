import Mathlib.Analysis.SpecialFunctions.Log.Basic

/-!
# A finite checkpoint law and its causal implementation

Backward masses on the complete native move tree construct the transition
law. Positive trajectories are exactly terminally accepted trajectories;
the law never samples a failed trajectory and retries it. The finite move
grammar, reference weights and terminal checkpoint are explicit data.
-/

set_option autoImplicit false

namespace OPH.SourceCheckpointPolicy
noncomputable section

variable {S A : Type*} [Fintype A]

def mass (step : S → A → S) (weight : S → A → ℝ) (terminal : S → ℝ) : ℕ → S → ℝ
  | 0, s => terminal s
  | n+1, s => ∑ a, weight s a * mass step weight terminal n (step s a)

def kernel (step : S → A → S) (weight : S → A → ℝ) (terminal : S → ℝ)
    (n : ℕ) (s : S) (a : A) : ℝ :=
  weight s a * mass step weight terminal n (step s a) /
    mass step weight terminal (n+1) s

def pathWeight (step : S → A → S) (weight : S → A → ℝ) : S → List A → ℝ
  | _, [] => 1
  | s, a::as => weight s a * pathWeight step weight (step s a) as

def finish (step : S → A → S) : S → List A → S
  | s, [] => s
  | s, a::as => finish step (step s a) as

def pathLaw (step : S → A → S) (weight : S → A → ℝ) (terminal : S → ℝ) :
    S → List A → ℝ
  | _, [] => 1
  | s, a::as => kernel step weight terminal as.length s a *
      pathLaw step weight terminal (step s a) as

theorem mass_nonneg (step : S → A → S) (weight : S → A → ℝ) (terminal : S → ℝ)
    (hw : ∀ s a, 0 ≤ weight s a) (ht : ∀ s, 0 ≤ terminal s) (n : ℕ) (s : S) :
    0 ≤ mass step weight terminal n s := by
  induction n generalizing s with
  | zero => exact ht s
  | succ n ih => exact Finset.sum_nonneg (fun a _ => mul_nonneg (hw s a) (ih _))

theorem kernel_nonneg (step : S → A → S) (weight : S → A → ℝ) (terminal : S → ℝ)
    (hw : ∀ s a, 0 ≤ weight s a) (ht : ∀ s, 0 ≤ terminal s)
    (n : ℕ) (s : S) (a : A) : 0 ≤ kernel step weight terminal n s a :=
  div_nonneg (mul_nonneg (hw s a) (mass_nonneg step weight terminal hw ht n _))
    (mass_nonneg step weight terminal hw ht (n+1) s)

theorem kernel_normalized (step : S → A → S) (weight : S → A → ℝ)
    (terminal : S → ℝ) (n : ℕ) (s : S)
    (h : mass step weight terminal (n+1) s ≠ 0) :
    ∑ a, kernel step weight terminal n s a = 1 := by
  simp only [kernel, div_eq_mul_inv, ← Finset.sum_mul]
  exact mul_inv_cancel₀ h

/-- The old proposal probability is retained at a move exactly when that
move leaves the continuation mass unchanged. A fixed deadline generally
couples successive moves, even for an IID reference. -/
theorem kernel_eq_reference_iff (step : S → A → S) (weight : S → A → ℝ)
    (terminal : S → ℝ) (n : ℕ) (s : S) (a : A)
    (hw : weight s a ≠ 0) (hm : mass step weight terminal (n+1) s ≠ 0) :
    kernel step weight terminal n s a = weight s a ↔
      mass step weight terminal n (step s a) = mass step weight terminal (n+1) s := by
  unfold kernel
  rw [div_eq_iff hm]
  constructor
  · exact mul_left_cancel₀ hw
  · intro h
    rw [h]

theorem mass_child_le (step : S → A → S) (weight : S → A → ℝ) (terminal : S → ℝ)
    (hw : ∀ s a, 0 ≤ weight s a) (ht : ∀ s, 0 ≤ terminal s)
    (n : ℕ) (s : S) (a : A) :
    weight s a * mass step weight terminal n (step s a) ≤
      mass step weight terminal (n+1) s := by
  exact Finset.single_le_sum (fun b _ => mul_nonneg (hw s b)
    (mass_nonneg step weight terminal hw ht n _)) (Finset.mem_univ a)

theorem positive_path_mass (step : S → A → S) (weight : S → A → ℝ)
    (terminal : S → ℝ) (hw : ∀ s a, 0 < weight s a) (ht : ∀ s, 0 ≤ terminal s)
    (s : S) (as : List A) (h : 0 < terminal (finish step s as)) :
    0 < mass step weight terminal as.length s := by
  induction as generalizing s with
  | nil => exact h
  | cons a as ih =>
    exact lt_of_lt_of_le (mul_pos (hw s a) (ih _ h))
      (mass_child_le step weight terminal (fun s a => (hw s a).le) ht as.length s a)

/-- Positive mass is exactly existence of a native word in the complete
finite tree. Decidability of the terminal predicate is a separate question. -/
theorem mass_positive_iff (step : S → A → S) (weight : S → A → ℝ)
    (terminal : S → ℝ) (hw : ∀ s a, 0 < weight s a) (ht : ∀ s, 0 ≤ terminal s)
    (n : ℕ) (s : S) :
    0 < mass step weight terminal n s ↔
      ∃ as : List A, as.length = n ∧ 0 < terminal (finish step s as) := by
  induction n generalizing s with
  | zero =>
    simp [mass,finish]
  | succ n ih =>
    constructor
    · intro h
      obtain ⟨a, _, ha⟩ := (Finset.sum_pos_iff_of_nonneg
        (fun a _ => mul_nonneg (hw s a).le
          (mass_nonneg step weight terminal (fun s a => (hw s a).le) ht n _))).mp h
      have hc : 0 < mass step weight terminal n (step s a) :=
        (mul_pos_iff_of_pos_left (hw s a)).mp ha
      obtain ⟨as, hl, he⟩ := (ih (step s a)).mp hc
      exact ⟨a::as, by simp [hl], he⟩
    · rintro ⟨as, hl, he⟩
      simpa [hl] using positive_path_mass step weight terminal hw ht s as he

private def consEquiv (n : ℕ) : A × (Fin n → A) ≃ (Fin (n+1) → A) where
  toFun p := Fin.cons p.1 p.2
  invFun w := (w 0, Fin.tail w)
  left_inv _ := by simp
  right_inv w := Fin.cons_self_tail w

/-- The backward recursion sums every native word exactly once. -/
theorem mass_eq_word_sum (step : S → A → S) (weight : S → A → ℝ)
    (terminal : S → ℝ) (n : ℕ) (s : S) :
    mass step weight terminal n s = ∑ w : Fin n → A,
      pathWeight step weight s (List.ofFn w) * terminal (finish step s (List.ofFn w)) := by
  induction n generalizing s with
  | zero => simp [mass,pathWeight,finish]
  | succ n ih =>
    rw [← (consEquiv (A := A) n).sum_comp
      (fun w => pathWeight step weight s (List.ofFn w) * terminal (finish step s (List.ofFn w))),
      Fintype.sum_prod_type]
    simp only [consEquiv,Equiv.coe_fn_mk,List.ofFn_succ,Fin.cons_zero,Fin.cons_succ,
      pathWeight,finish,mass]
    simp_rw [mul_assoc, ← Finset.mul_sum, ← ih]

omit [Fintype A] in
theorem pathWeight_positive (step : S → A → S) (weight : S → A → ℝ)
    (hw : ∀ s a, 0 < weight s a) (s : S) (as : List A) :
    0 < pathWeight step weight s as := by
  induction as generalizing s with
  | nil => exact zero_lt_one
  | cons a as ih => exact mul_pos (hw s a) (ih _)

/-- Products of causal transition probabilities telescope to the globally
conditioned reference. Unreachable children carry zero probability. -/
theorem pathLaw_formula (step : S → A → S) (weight : S → A → ℝ)
    (terminal : S → ℝ) (hw : ∀ s a, 0 < weight s a) (ht : ∀ s, 0 ≤ terminal s)
    (s : S) (as : List A) (hm : 0 < mass step weight terminal as.length s) :
    pathLaw step weight terminal s as =
      pathWeight step weight s as * terminal (finish step s as) /
        mass step weight terminal as.length s := by
  induction as generalizing s with
  | nil => simp only [pathLaw,pathWeight,finish,List.length_nil,mass,one_mul] at *
           exact (div_self (ne_of_gt hm)).symm
  | cons a as ih =>
    have hn := mass_nonneg step weight terminal (fun s a => (hw s a).le) ht as.length (step s a)
    by_cases hz : mass step weight terminal as.length (step s a) = 0
    · have hend : terminal (finish step (step s a) as) = 0 := by
        by_contra he
        have hp := positive_path_mass step weight terminal hw ht (step s a) as
          (lt_of_le_of_ne (ht _) (Ne.symm he))
        rw [hz] at hp
        exact (lt_irrefl 0 hp)
      simp [pathLaw,kernel,pathWeight,finish,hz,hend]
    · have hp : 0 < mass step weight terminal as.length (step s a) := lt_of_le_of_ne hn (Ne.symm hz)
      simp only [pathLaw,kernel,pathWeight,finish,List.length_cons]
      rw [ih _ hp]
      field_simp

/-- Terminal certainty follows from the derived law, without a fairness or
independent-proposal hypothesis. This does not derive the terminal demand. -/
theorem pathLaw_positive_iff (step : S → A → S) (weight : S → A → ℝ)
    (terminal : S → ℝ) (hw : ∀ s a, 0 < weight s a) (ht : ∀ s, 0 ≤ terminal s)
    (s : S) (as : List A) (hm : 0 < mass step weight terminal as.length s) :
    0 < pathLaw step weight terminal s as ↔ 0 < terminal (finish step s as) := by
  rw [pathLaw_formula step weight terminal hw ht s as hm]
  exact (div_pos_iff_of_pos_right hm).trans
    (mul_pos_iff_of_pos_left (pathWeight_positive step weight hw s as))

/-- Full-word normalization includes the zero probability of every failed
history. It is not normalization of a retrospectively thinned sample. -/
theorem pathLaw_normalized (step : S → A → S) (weight : S → A → ℝ)
    (terminal : S → ℝ) (hw : ∀ s a, 0 < weight s a) (ht : ∀ s, 0 ≤ terminal s)
    (s : S) (n : ℕ) (hm : 0 < mass step weight terminal n s) :
    ∑ w : Fin n → A, pathLaw step weight terminal s (List.ofFn w) = 1 := by
  have he : ∀ w : Fin n → A,
      pathLaw step weight terminal s (List.ofFn w) =
        pathWeight step weight s (List.ofFn w) * terminal (finish step s (List.ofFn w)) /
          mass step weight terminal n s := by
    intro w
    simpa using pathLaw_formula step weight terminal hw ht s (List.ofFn w) (by simpa using hm)
  simp_rw [he,div_eq_mul_inv]
  rw [← Finset.sum_mul,← mass_eq_word_sum]
  exact mul_inv_cancel₀ (ne_of_gt hm)

end
end OPH.SourceCheckpointPolicy
