import Mathlib

/-!
# Finite-cut survival does not select a global-capacity action

One normalized connected-cut datum supplies a local blocked-event expectation
`q` and local survival `1 - q`.  This file first exhibits two positive actions
on a global capacity and then separates three blocked-event semantics: no
capacity action, one-class blocking, and six-class-total blocking.  Every
action composes under disconnected-cut grouping, while the alternatives
disagree after one cut on their stated positive domains.

The theorem is a typed non-identifiability statement.  It does not select a
physical blocked event, a global capacity, or a cosmological value.
-/

namespace OPH.GlobalCapacityAttachment

/-- A local one-class survival factor. -/
def survival (q : ℝ) : ℝ := 1 - q

/-- A completion in which the local cut record has no action on capacity. -/
def neutralAction (_q : ℝ) (_cuts : ℕ) (capacity : ℝ) : ℝ := capacity

/-- A completion in which every disconnected cut contributes one factor. -/
def multiplicativeAction (q : ℝ) (cuts : ℕ) (capacity : ℝ) : ℝ :=
  survival q ^ cuts * capacity

/-- The neutral completion has the identity at zero cuts. -/
theorem neutral_zero (q capacity : ℝ) :
    neutralAction q 0 capacity = capacity := rfl

/-- The multiplicative completion has the identity at zero cuts. -/
theorem multiplicative_zero (q capacity : ℝ) :
    multiplicativeAction q 0 capacity = capacity := by
  simp [multiplicativeAction]

/-- The neutral completion composes under every grouping of disconnected cuts. -/
theorem neutral_compose (q : ℝ) (m n : ℕ) (capacity : ℝ) :
    neutralAction q (m + n) capacity =
      neutralAction q m (neutralAction q n capacity) := rfl

/-- The multiplicative completion composes under every grouping of cuts. -/
theorem multiplicative_compose (q : ℝ) (m n : ℕ) (capacity : ℝ) :
    multiplicativeAction q (m + n) capacity =
      multiplicativeAction q m (multiplicativeAction q n capacity) := by
  simp only [multiplicativeAction, pow_add]
  ring

/-- Regrouping a finite partition of a cut count does not change the neutral action. -/
theorem neutral_refinement_partition (q : ℝ) (pieces : List ℕ) (capacity : ℝ) :
    neutralAction q pieces.sum capacity =
      pieces.foldr (fun cuts current => neutralAction q cuts current) capacity := by
  simp [neutralAction]

/-- Regrouping a finite partition of a cut count does not change the multiplicative action. -/
theorem multiplicative_refinement_partition
    (q : ℝ) (pieces : List ℕ) (capacity : ℝ) :
    multiplicativeAction q pieces.sum capacity =
      pieces.foldr (fun cuts current => multiplicativeAction q cuts current) capacity := by
  induction pieces with
  | nil => simp [multiplicativeAction]
  | cons head tail ih =>
      simp only [List.sum_cons, List.foldr_cons]
      rw [multiplicative_compose, ih]

/-- Both completion actions preserve positivity in the physical domain. -/
theorem multiplicative_positive {q capacity : ℝ}
    (_hq0 : 0 < q) (hq1 : q < 1) (hcapacity : 0 < capacity) (cuts : ℕ) :
    0 < multiplicativeAction q cuts capacity := by
  have hsurvival : 0 < survival q := by
    dsimp [survival]
    linarith
  exact mul_pos (pow_pos hsurvival cuts) hcapacity

/-- The two lawful completions disagree after one cut. -/
theorem completions_disagree_at_one {q capacity : ℝ}
    (hq : 0 < q) (hcapacity : 0 < capacity) :
    neutralAction q 1 capacity ≠ multiplicativeAction q 1 capacity := by
  simp only [neutralAction, multiplicativeAction, survival, pow_one]
  intro heq
  have : q * capacity = 0 := by linarith
  exact (mul_ne_zero (ne_of_gt hq) (ne_of_gt hcapacity)) this

/-- A monoid-action requirement leaves two distinct global completions. -/
theorem local_survival_does_not_select_global_action {q capacity : ℝ}
    (hq0 : 0 < q) (_hq1 : q < 1) (hcapacity : 0 < capacity) :
    (∀ m n, neutralAction q (m + n) capacity =
        neutralAction q m (neutralAction q n capacity)) ∧
    (∀ m n, multiplicativeAction q (m + n) capacity =
        multiplicativeAction q m (multiplicativeAction q n capacity)) ∧
    neutralAction q 1 capacity ≠ multiplicativeAction q 1 capacity := by
  exact ⟨fun m n => neutral_compose q m n capacity,
    fun m n => multiplicative_compose q m n capacity,
    completions_disagree_at_one hq0 hcapacity⟩

/-- One-class blocked-event factor. -/
def oneClassFactor (q : ℝ) : ℝ := 1 - q

/-- Six mutually exclusive equal classes, when that total projector exists. -/
def sixClassFactor (q : ℝ) : ℝ := 1 - 6 * q

/-- No action of the local blocked event on the global object. -/
def noCapacityActionFactor : ℝ := 1

/-- Every scalar factor defines a cut-count action on global capacity. -/
def factorAction (factor : ℝ) (cuts : ℕ) (capacity : ℝ) : ℝ :=
  factor ^ cuts * capacity

/-- Factor actions compose under disconnected-cut grouping. -/
theorem factorAction_compose (factor : ℝ) (m n : ℕ) (capacity : ℝ) :
    factorAction factor (m + n) capacity =
      factorAction factor m (factorAction factor n capacity) := by
  simp only [factorAction, pow_add]
  ring

/-- Factor actions are unchanged by finite regrouping of a cut count. -/
theorem factorAction_refinement_partition
    (factor : ℝ) (pieces : List ℕ) (capacity : ℝ) :
    factorAction factor pieces.sum capacity =
      pieces.foldr (fun cuts current => factorAction factor cuts current) capacity := by
  induction pieces with
  | nil => simp [factorAction]
  | cons head tail ih =>
      simp only [List.sum_cons, List.foldr_cons]
      rw [factorAction_compose, ih]

/-- A positive factor acts positively on a positive capacity. -/
theorem factorAction_positive {factor capacity : ℝ}
    (hfactor : 0 < factor) (hcapacity : 0 < capacity) (cuts : ℕ) :
    0 < factorAction factor cuts capacity := by
  exact mul_pos (pow_pos hfactor cuts) hcapacity

/-- The three candidate semantics are exactly distinct on `0 < q < 1/6`. -/
theorem blocked_event_factors_pairwise_distinct {q : ℝ}
    (hq0 : 0 < q) (_hq6 : q < 1 / 6) :
    oneClassFactor q ≠ sixClassFactor q ∧
      oneClassFactor q ≠ noCapacityActionFactor ∧
      sixClassFactor q ≠ noCapacityActionFactor := by
  constructor
  · dsimp [oneClassFactor, sixClassFactor]
    intro h
    linarith
  constructor
  · dsimp [oneClassFactor, noCapacityActionFactor]
    intro h
    linarith
  · dsimp [sixClassFactor, noCapacityActionFactor]
    intro h
    linarith

/-- The three blocked-event actions are positive on the declared domain. -/
theorem blocked_event_actions_positive {q capacity : ℝ}
    (hq0 : 0 < q) (hq6 : q < 1 / 6) (hcapacity : 0 < capacity) (cuts : ℕ) :
    0 < factorAction noCapacityActionFactor cuts capacity ∧
      0 < factorAction (oneClassFactor q) cuts capacity ∧
      0 < factorAction (sixClassFactor q) cuts capacity := by
  have hOne : 0 < oneClassFactor q := by
    dsimp [oneClassFactor]
    linarith
  have hSix : 0 < sixClassFactor q := by
    dsimp [sixClassFactor]
    linarith
  exact ⟨factorAction_positive (by norm_num [noCapacityActionFactor]) hcapacity cuts,
    factorAction_positive hOne hcapacity cuts,
    factorAction_positive hSix hcapacity cuts⟩

/-- The three lawful blocked-event completions disagree after one cut. -/
theorem blocked_event_actions_pairwise_distinct_at_one {q capacity : ℝ}
    (hq0 : 0 < q) (hq6 : q < 1 / 6) (hcapacity : 0 < capacity) :
    factorAction noCapacityActionFactor 1 capacity ≠
        factorAction (oneClassFactor q) 1 capacity ∧
      factorAction noCapacityActionFactor 1 capacity ≠
        factorAction (sixClassFactor q) 1 capacity ∧
      factorAction (oneClassFactor q) 1 capacity ≠
        factorAction (sixClassFactor q) 1 capacity := by
  have hfactors := blocked_event_factors_pairwise_distinct hq0 hq6
  simp only [factorAction, pow_one]
  constructor
  · intro h
    apply hfactors.2.1
    exact mul_right_cancel₀ (ne_of_gt hcapacity) h.symm
  constructor
  · intro h
    apply hfactors.2.2
    exact mul_right_cancel₀ (ne_of_gt hcapacity) h.symm
  · intro h
    apply hfactors.1
    exact mul_right_cancel₀ (ne_of_gt hcapacity) h

-- Axiom audit: all results use only standard Mathlib axioms.
#print axioms neutral_compose
#print axioms multiplicative_compose
#print axioms neutral_refinement_partition
#print axioms multiplicative_refinement_partition
#print axioms multiplicative_positive
#print axioms completions_disagree_at_one
#print axioms local_survival_does_not_select_global_action
#print axioms blocked_event_factors_pairwise_distinct
#print axioms factorAction_compose
#print axioms factorAction_refinement_partition
#print axioms factorAction_positive
#print axioms blocked_event_actions_positive
#print axioms blocked_event_actions_pairwise_distinct_at_one

end OPH.GlobalCapacityAttachment
