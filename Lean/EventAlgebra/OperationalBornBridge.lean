import EventAlgebra.FiniteBuschGleason

set_option autoImplicit false

namespace EventAlgebra

/-!
# Born representation from complete record probabilities

Each experiment has a finite complete record space and a normalized
nonnegative record law. A separate positive, complete family of matrix
effects represents its records. An event's probability and effect are
their respective finite sums; their disjoint-union laws are proved from
these definitions, without a trace definition of probability. The collection
of experiments is not assumed finite.

Two operational hypotheses connect different experiments: equal effects
have equal event probabilities, and every coexistent pair of complex matrix
effects has disjoint event implementations in one experiment. Pairing any
effect with zero also gives universal effect realizability. The chosen-implementation valuation is
then well defined on effects and additive. The finite Busch--Gleason theorem
produces a unique density matrix representing every implemented event.

These universal hypotheses concern the full complex effect algebra. A
classical native record construction verifies them only on its represented
commutative domain. Finite complex tomography supplies neither universal
probability equivalence nor joint refinement for all coexistent effects.
No trace-defined runtime is used as independent evidence for these clauses.
-/

open Matrix
open scoped ComplexOrder

universe u v

variable {C : Type u} {Ω : C → Type v} [∀ c, Fintype (Ω c)]

/-- A complete probability law for every finite record experiment. -/
structure CompleteRecordWeights (C : Type u) (Ω : C → Type v)
    [∀ c, Fintype (Ω c)] where
  weight : ∀ c, Ω c → ℝ
  nonneg : ∀ c x, 0 ≤ weight c x
  total : ∀ c, ∑ x, weight c x = 1

/-- A separate positive matrix representation of complete record outcomes.
No probability or trace pairing occurs in these fields. -/
structure CompleteRecordEffects (d : ℕ) (C : Type u) (Ω : C → Type v)
    [∀ c, Fintype (Ω c)] where
  atom : ∀ c, Ω c → Matrix (Fin d) (Fin d) ℂ
  positive : ∀ c x, (atom c x).PosSemidef
  complete : ∀ c, ∑ x, atom c x = 1

variable {d : ℕ}

noncomputable section

/-- The probability of a public record event is the sum of its atomic
record probabilities. -/
def recordEventWeight (W : CompleteRecordWeights C Ω) (c : C)
    (S : Finset (Ω c)) : ℝ := ∑ x ∈ S, W.weight c x

/-- Coarse-graining records sums their represented effects. -/
def recordEventEffect (R : CompleteRecordEffects d C Ω) (c : C)
    (S : Finset (Ω c)) : Matrix (Fin d) (Fin d) ℂ := ∑ x ∈ S, R.atom c x

theorem recordEventWeight_nonneg (W : CompleteRecordWeights C Ω)
    (c : C) (S : Finset (Ω c)) : 0 ≤ recordEventWeight W c S :=
  Finset.sum_nonneg fun x _ => W.nonneg c x

theorem recordEventWeight_univ (W : CompleteRecordWeights C Ω) (c : C) :
    recordEventWeight W c Finset.univ = 1 := W.total c

theorem recordEventWeight_le_one (W : CompleteRecordWeights C Ω)
    (c : C) (S : Finset (Ω c)) : recordEventWeight W c S ≤ 1 := by
  rw [← recordEventWeight_univ W c]
  exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ S)
    (fun x _ _ => W.nonneg c x)

theorem recordEventEffect_univ (R : CompleteRecordEffects d C Ω) (c : C) :
    recordEventEffect R c Finset.univ = 1 := R.complete c

variable [∀ c, DecidableEq (Ω c)]

/-- Disjoint record unions are additive because their finite sums have no
overlapping terms. This precedes any effect quotient or Born rule. -/
theorem recordEventWeight_union (W : CompleteRecordWeights C Ω)
    (c : C) {S T : Finset (Ω c)} (h : Disjoint S T) :
    recordEventWeight W c (S ∪ T) = recordEventWeight W c S + recordEventWeight W c T :=
  Finset.sum_union h

theorem recordEventEffect_union (R : CompleteRecordEffects d C Ω)
    (c : C) {S T : Finset (Ω c)} (h : Disjoint S T) :
    recordEventEffect R c (S ∪ T) = recordEventEffect R c S + recordEventEffect R c T :=
  Finset.sum_union h

/-- Every represented record event is an effect: its matrix and its
complement are sums of positive outcome matrices. -/
theorem recordEventEffect_isEffect (R : CompleteRecordEffects d C Ω)
    (c : C) (S : Finset (Ω c)) : IsEffect (recordEventEffect R c S) := by
  have hpos (T : Finset (Ω c)) : (recordEventEffect R c T).PosSemidef :=
    posSemidef_sum T (fun x _ => R.positive c x)
  refine ⟨hpos S, ?_⟩
  have hpartition : recordEventEffect R c S +
      recordEventEffect R c (Finset.univ \ S) = 1 := by
    have hdisj : Disjoint S (Finset.univ \ S) :=
      Finset.disjoint_left.mpr (fun _ hS hdiff => (Finset.mem_sdiff.mp hdiff).2 hS)
    rw [← recordEventEffect_union R c hdisj,
      Finset.union_sdiff_of_subset (Finset.subset_univ S), recordEventEffect_univ]
  have hcompl : 1 - recordEventEffect R c S =
      recordEventEffect R c (Finset.univ \ S) := by
    rw [← hpartition]
    abel
  rw [hcompl]
  exact hpos _

/-- An implementation is an experiment and one event of its complete
finite public record space. -/
abbrev ImplementedRecordEvent (C : Type u) (Ω : C → Type v) :=
  (c : C) × Finset (Ω c)

/-- Operational conditions required to pass from implemented record events
to the entire complex matrix effect algebra. The refinement clause supplies
one common experiment for each coexistent pair; it asserts no probability
equation. Probability additivity is derived from disjoint record sums. -/
structure OperationalEffectQuotient (W : CompleteRecordWeights C Ω)
    (R : CompleteRecordEffects d C Ω) : Prop where
  implementation_independent : ∀ (c c' : C) (S : Finset (Ω c)) (T : Finset (Ω c')),
    recordEventEffect R c S = recordEventEffect R c' T →
      recordEventWeight W c S = recordEventWeight W c' T
  coexistent_refinement : ∀ {E F}, IsEffect E → IsEffect F → IsEffect (E + F) →
    ∃ (c : C) (S T : Finset (Ω c)), Disjoint S T ∧
      recordEventEffect R c S = E ∧ recordEventEffect R c T = F

omit [∀ c, DecidableEq (Ω c)] in
/-- Joint refinement applied to an effect and zero supplies an
implementation of every effect; no separate realizability premise is needed. -/
theorem OperationalEffectQuotient.realizes
    {W : CompleteRecordWeights C Ω} {R : CompleteRecordEffects d C Ω}
    (Q : OperationalEffectQuotient W R) {E : Matrix (Fin d) (Fin d) ℂ}
    (hE : IsEffect E) :
    ∃ e : ImplementedRecordEvent C Ω, recordEventEffect R e.1 e.2 = E := by
  obtain ⟨c, S, _, _, hS, _⟩ := Q.coexistent_refinement
    (E := E) (F := 0) hE isEffect_zero (by simpa using hE)
  exact ⟨⟨c, S⟩, hS⟩

variable (W : CompleteRecordWeights C Ω) (R : CompleteRecordEffects d C Ω)
    (Q : OperationalEffectQuotient W R)

/-- The valuation chooses any implementing record event. Values outside
the effect interval are set to zero and play no role in the theorem. -/
def operationalEffectValuation (E : Matrix (Fin d) (Fin d) ℂ) : ℝ := by
  classical
  exact if hE : IsEffect E then
    let e := Classical.choose (Q.realizes hE)
    recordEventWeight W e.1 e.2
  else 0

/-- The chosen valuation agrees with every implementation, so it is
independent of the representative selected by choice. -/
theorem operationalEffectValuation_eq_eventWeight (c : C) (S : Finset (Ω c)) :
    operationalEffectValuation W R Q (recordEventEffect R c S) = recordEventWeight W c S := by
  rw [operationalEffectValuation, dif_pos (recordEventEffect_isEffect R c S)]
  exact Q.implementation_independent _ _ _ _
    (Classical.choose_spec (Q.realizes (recordEventEffect_isEffect R c S)))

theorem operationalEffectValuation_eq_of_realizes
    {E : Matrix (Fin d) (Fin d) ℂ} (c : C) (S : Finset (Ω c))
    (h : recordEventEffect R c S = E) :
    operationalEffectValuation W R Q E = recordEventWeight W c S := by
  rw [← h]
  exact operationalEffectValuation_eq_eventWeight W R Q c S

/-- The operational quotient is an additive effect valuation. Range and
normalization descend from complete record laws. The nontrivial additivity
step uses a common disjoint refinement and the proved finite-sum union law. -/
theorem operationalEffectValuation_isEffectValuation :
    IsEffectValuation (operationalEffectValuation W R Q) := by
  constructor
  · intro E hE
    obtain ⟨⟨c, S⟩, h⟩ := Q.realizes hE
    rw [operationalEffectValuation_eq_of_realizes W R Q c S h]
    exact recordEventWeight_nonneg W c S
  · intro E hE
    obtain ⟨⟨c, S⟩, h⟩ := Q.realizes hE
    rw [operationalEffectValuation_eq_of_realizes W R Q c S h]
    exact recordEventWeight_le_one W c S
  · obtain ⟨⟨c, _⟩, _⟩ := Q.realizes (isEffect_one (d := d))
    rw [operationalEffectValuation_eq_of_realizes W R Q c Finset.univ
      (recordEventEffect_univ R c), recordEventWeight_univ]
  · intro E F hE hF hEF
    obtain ⟨c, S, T, hdisj, hS, hT⟩ := Q.coexistent_refinement hE hF hEF
    have hsum : recordEventEffect R c (S ∪ T) = E + F := by
      rw [recordEventEffect_union R c hdisj, hS, hT]
    rw [operationalEffectValuation_eq_of_realizes W R Q c (S ∪ T) hsum,
      operationalEffectValuation_eq_of_realizes W R Q c S hS,
      operationalEffectValuation_eq_of_realizes W R Q c T hT,
      recordEventWeight_union W c hdisj]

/-- The existing finite Busch--Gleason theorem supplies the density matrix
after the valuation has been derived from record probabilities. -/
theorem operational_busch_gleason :
    ∃! ρ : Matrix (Fin d) (Fin d) ℂ,
      IsState ρ ∧ ∀ E, IsEffect E →
        operationalEffectValuation W R Q E = (bornWeight ρ E).re :=
  finite_busch_gleason (operationalEffectValuation_isEffectValuation W R Q)

include Q in
/-- One unique density matrix represents the independently supplied
probabilities of every implemented record event. The Born formula occurs
in this conclusion, rather than in the record-law definitions or quotient
hypotheses. -/
theorem operational_record_probabilities_unique_born :
    ∃! ρ : Matrix (Fin d) (Fin d) ℂ,
      IsState ρ ∧ ∀ (c : C) (S : Finset (Ω c)),
        recordEventWeight W c S = (bornWeight ρ (recordEventEffect R c S)).re := by
  obtain ⟨ρ, ⟨hρ, hrep⟩, huniq⟩ := operational_busch_gleason W R Q
  refine ⟨ρ, ⟨hρ, ?_⟩, ?_⟩
  · intro c S
    rw [← operationalEffectValuation_eq_eventWeight W R Q c S]
    exact hrep _ (recordEventEffect_isEffect R c S)
  · rintro σ ⟨hσ, hrecords⟩
    apply huniq σ
    refine ⟨hσ, ?_⟩
    intro E hE
    obtain ⟨⟨c, S⟩, h⟩ := Q.realizes hE
    rw [← h, operationalEffectValuation_eq_eventWeight W R Q c S]
    exact hrecords c S

end

/-! ## The diagonal native grammar cannot supply the full quotient

Complete records of weighted partial-injection branches have diagonal
effects. Summing records preserves that property. This blocks the universal
quotient whenever its target algebra contains an off-diagonal effect. The
separate native-interface certificate proves the partial-injection premise
for the actual C/E branches and its closure under finite classical feedback;
the coherent proper-code mixer lies outside that grammar.
-/

/-- Diagonal outcome effects stay diagonal under record coarse graining. -/
theorem recordEventEffect_offDiagonal_zero
    (R : CompleteRecordEffects d C Ω)
    (hdiag : ∀ c x (i j : Fin d), i ≠ j → R.atom c x i j = 0)
    (c : C) (S : Finset (Ω c)) (i j : Fin d) (hij : i ≠ j) :
    recordEventEffect R c S i j = 0 := by
  rw [recordEventEffect, Matrix.sum_apply]
  exact Finset.sum_eq_zero fun x _ => hdiag c x i j hij

/-- A grammar with diagonal record effects cannot implement the entire
complex effect algebra when it contains an off-diagonal effect. This is a
failure of effect reachability, independent of the probability assignment. -/
theorem no_operationalEffectQuotient_of_offDiagonal_effect
    (W : CompleteRecordWeights C Ω) (R : CompleteRecordEffects d C Ω)
    (hdiag : ∀ c x (i j : Fin d), i ≠ j → R.atom c x i j = 0)
    {E : Matrix (Fin d) (Fin d) ℂ} (hE : IsEffect E)
    (i j : Fin d) (hij : i ≠ j) (hEij : E i j ≠ 0) :
    ¬ OperationalEffectQuotient W R := by
  intro Q
  obtain ⟨⟨c, S⟩, h⟩ := Q.realizes hE
  apply hEij
  rw [← h]
  exact recordEventEffect_offDiagonal_zero R hdiag c S i j hij

/-- A non-diagonal projection in the smallest complex quantum carrier. -/
noncomputable def diagonalObstructionEffect : Matrix (Fin 2) (Fin 2) ℂ :=
  !![1 / 2, 1 / 2; 1 / 2, 1 / 2]

theorem diagonalObstructionEffect_isEvent : IsEvent diagonalObstructionEffect := by
  constructor
  · ext i j
    fin_cases i <;> fin_cases j <;> norm_num [diagonalObstructionEffect]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [diagonalObstructionEffect, Matrix.mul_apply, Fin.sum_univ_succ]

/-- Even a qubit's full operational quotient is incompatible with a
diagonal-only record-effect family. -/
theorem diagonal_qubit_records_no_operationalEffectQuotient
    (W : CompleteRecordWeights C Ω) (R : CompleteRecordEffects 2 C Ω)
    (hdiag : ∀ c x (i j : Fin 2), i ≠ j → R.atom c x i j = 0) :
    ¬ OperationalEffectQuotient W R := by
  apply no_operationalEffectQuotient_of_offDiagonal_effect W R hdiag
    diagonalObstructionEffect_isEvent.isEffect 0 1
  · decide
  · norm_num [diagonalObstructionEffect]

#print axioms recordEventWeight_union
#print axioms recordEventEffect_union
#print axioms recordEventEffect_isEffect
#print axioms operationalEffectValuation_eq_eventWeight
#print axioms operationalEffectValuation_isEffectValuation
#print axioms operational_busch_gleason
#print axioms operational_record_probabilities_unique_born
#print axioms recordEventEffect_offDiagonal_zero
#print axioms no_operationalEffectQuotient_of_offDiagonal_effect
#print axioms diagonalObstructionEffect_isEvent
#print axioms diagonal_qubit_records_no_operationalEffectQuotient

end EventAlgebra
