import EventAlgebra.OperationalBornBridge

set_option autoImplicit false
namespace EventAlgebra.NativePublicEventGrammar
open Matrix
open scoped ComplexOrder
noncomputable section
abbrev Mat (d : ℕ) := Matrix (Fin d) (Fin d) ℂ


/-- The source's declared complete classical preparation law supplies a
positive normalized density, including every sharp initialized word. -/
def initialDensity {d : ℕ} (w : Fin d → ℝ) : Mat d :=
  Matrix.diagonal (fun i => (w i : ℂ))

theorem initialDensity_isState {d : ℕ} (w : Fin d → ℝ)
    (hpos : ∀ i, 0 ≤ w i) (htotal : ∑ i, w i = 1) : IsState (initialDensity w) := by
  constructor
  · exact Matrix.posSemidef_diagonal_iff.mpr (fun i => Complex.zero_le_real.mpr (hpos i))
  · simp only [initialDensity, Matrix.trace_diagonal]
    exact_mod_cast htotal

/-- Finite decision trees with complete primitive quantum branches. Classical
coins are scalar-identity branches. Early stops and refusals are retained leaves. -/
inductive Program (d : ℕ) where
  | halt : Program d
  | step (m : ℕ) (K : Fin m → Mat d)
      (complete : ∑ i, (K i)ᴴ * K i = 1)
      (next : Fin m → Program d) : Program d

def Record {d : ℕ} : Program d → Type
  | .halt => Unit
  | .step m _ _ next => (i : Fin m) × Record (next i)

instance recordFintype {d : ℕ} (p : Program d) : Fintype (Record p) := by
  induction p with
  | halt => exact inferInstanceAs (Fintype Unit)
  | step m K h next ih =>
    letI : ∀ i, Fintype (Record (next i)) := ih
    exact inferInstanceAs (Fintype ((i : Fin m) × Record (next i)))

instance recordDecidableEq {d : ℕ} (p : Program d) : DecidableEq (Record p) :=
  Classical.decEq _

/-- Actual chronological post-measurement maps, with every record retained. -/
def branch {d : ℕ} : (p : Program d) → Record p → Mat d
  | .halt, _ => 1
  | .step _ K _ next, ⟨i, r⟩ => branch (next i) r * K i

/-- Composite completeness for every finite adaptive program. -/
theorem branch_complete {d : ℕ} (p : Program d) :
    ∑ r, (branch p r)ᴴ * branch p r = 1 := by
  induction p with
  | halt => simp [Record, branch]
  | step m K h next ih =>
    change (∑ r : (i : Fin m) × Record (next i),
      (branch (next r.1) r.2 * K r.1)ᴴ * (branch (next r.1) r.2 * K r.1)) = 1
    rw [Fintype.sum_sigma]
    simp only [conjTranspose_mul]
    calc
      (∑ i, ∑ r, (K i)ᴴ * (branch (next i) r)ᴴ * (branch (next i) r * K i))
          = ∑ i, (K i)ᴴ * (∑ r, (branch (next i) r)ᴴ * branch (next i) r) * K i := by
            apply Finset.sum_congr rfl
            intro i _
            simp only [Finset.mul_sum, Finset.sum_mul, Matrix.mul_assoc]
      _ = 1 := by simp only [ih, mul_one]; exact h

/-- Positive complete event atoms derived from the actual history maps. -/
def effects {d : ℕ} (p : Program d) : CompleteRecordEffects d Unit (fun _ => Record p) where
  atom _ r := (branch p r)ᴴ * branch p r
  positive _ _ := posSemidef_conjTranspose_mul_self _
  complete _ := branch_complete p

/-- Primitive runtime quantum semantics, before defining public event effects. -/
def weight {d : ℕ} (p : Program d) (ρ : Mat d) (r : Record p) : ℝ :=
  (Matrix.trace (branch p r * ρ * (branch p r)ᴴ)).re

theorem weight_eq_pairing {d : ℕ} (p : Program d) (ρ : Mat d) (r : Record p) :
    weight p ρ r = (bornWeight ρ ((branch p r)ᴴ * branch p r)).re := by
  unfold weight bornWeight
  rw [trace_mul_cycle, trace_mul_comm]

/-- A normalized positive initial density gives a complete public record law. -/
def weights {d : ℕ} (p : Program d) (ρ : Mat d) (hρ : IsState ρ) :
    CompleteRecordWeights Unit (fun _ => Record p) where
  weight _ r := weight p ρ r
  nonneg _ r := by
    have h := (hρ.1.mul_mul_conjTranspose_same (branch p r)).trace_nonneg
    exact (Complex.le_def.mp h).1
  total _ := by
    simp only [weight_eq_pairing]
    rw [← Complex.re_sum, ← bornWeight_sum, branch_complete, bornWeight_one hρ]
    rfl

def eventEffect {d : ℕ} (p : Program d) (S : Finset (Record p)) : Mat d :=
  recordEventEffect (effects p) () S

def eventWeight {d : ℕ} (p : Program d) (ρ : Mat d) (S : Finset (Record p)) : ℝ :=
  ∑ r ∈ S, weight p ρ r

/-- One positive normalized state represents every generated public effect.
The quantum branch semantics and initial state are supplied, not inferred. -/
theorem eventWeight_eq_pairing {d : ℕ} (p : Program d) (ρ : Mat d)
    (S : Finset (Record p)) :
    eventWeight p ρ S = (bornWeight ρ (eventEffect p S)).re := by
  simp only [eventWeight, eventEffect, recordEventEffect, effects, weight_eq_pairing]
  rw [← Complex.re_sum, bornWeight_sum]

/-- Independent operator equality certifies reading equivalence, without
comparing probabilities on one target state. -/
theorem same_effect_same_probability {d : ℕ} (p q : Program d) (ρ : Mat d)
    (S : Finset (Record p)) (T : Finset (Record q))
    (h : eventEffect p S = eventEffect q T) :
    eventWeight p ρ S = eventWeight q ρ T := by
  rw [eventWeight_eq_pairing, eventWeight_eq_pairing, h]

theorem disjoint_union {d : ℕ} (p : Program d) (ρ : Mat d)
    (S T : Finset (Record p)) (h : Disjoint S T) :
    eventWeight p ρ (S ∪ T) = eventWeight p ρ S + eventWeight p ρ T := by
  exact Finset.sum_union h

theorem public_effect {d : ℕ} (p : Program d) (S : Finset (Record p)) :
    IsEffect (eventEffect p S) := recordEventEffect_isEffect (effects p) () S

/-- Every label of a record-map refinement retains its mass. -/
theorem labelled_refinement {d : ℕ} {β : Type*} [Fintype β] [DecidableEq β]
    (p : Program d) (ρ : Mat d) (label : Record p → β) :
    ∑ b, eventWeight p ρ (Finset.univ.filter (fun r => label r = b)) =
      eventWeight p ρ Finset.univ := by
  simp only [eventWeight, Finset.sum_filter]
  rw [Finset.sum_comm]
  simp


/-- A complete operator signature independent of the target state. The
symmetric/antisymmetric matrix units are real linear combinations of
coordinate and pair-superposition probe projectors. -/
theorem operator_signature_ext {d : ℕ} (E F : Mat d)
    (hs : ∀ j k, (E * symUnit j k).trace = (F * symUnit j k).trace)
    (ha : ∀ j k, (E * skewUnit j k).trace = (F * skewUnit j k).trace) : E = F := by
  ext j k
  have hS := hs j k
  have hA := ha j k
  rw [trace_mul_symUnit, trace_mul_symUnit] at hS
  rw [trace_mul_skewUnit, trace_mul_skewUnit, ← mul_sub, ← mul_sub] at hA
  have hA' := mul_left_cancel₀ Complex.I_ne_zero hA
  have h2 : (2 : ℂ) * E j k = 2 * F j k := by
    linear_combination hS - hA'
  exact mul_left_cancel₀ two_ne_zero h2

/-- Context normalization does not force refinement consistency. Squaring
primitive probabilities and renormalizing each complete record produces
fine weights(1/6,1/6,2/3), coarse weights(1/2,1/2). Both contexts sum to one,
but the two events representing the same scalar effect I/2 disagree. -/
theorem normalized_context_sampler_refinement_counterexample :
    (0 ≤ (1/6 : ℚ) ∧ 0 ≤ (2/3 : ℚ) ∧
     (1/6 : ℚ) + 1/6 + 2/3 = 1 ∧ (1/2 : ℚ) + 1/2 = 1) ∧
    (1/4 : ℚ) + 1/4 = 1/2 ∧
    (1/6 : ℚ) + 1/6 ≠ 1/2 := by norm_num

#print axioms initialDensity_isState
#print axioms branch_complete
#print axioms eventWeight_eq_pairing
#print axioms same_effect_same_probability
#print axioms public_effect
#print axioms labelled_refinement
#print axioms operator_signature_ext
#print axioms normalized_context_sampler_refinement_counterexample
end
end EventAlgebra.NativePublicEventGrammar
