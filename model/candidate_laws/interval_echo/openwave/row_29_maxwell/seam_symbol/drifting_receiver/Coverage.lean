import Mathlib.Tactic
import Mathlib.Analysis.Convex.SpecificFunctions.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Series

/-! A finite receiver tree may choose fresh dark and bright response rates
at every history. The eight branches keep every arm and click/miss outcome.
The fourth arm may be any other calibration or signal. -/

namespace IntervalEcho.FiniteReceiver

open Real

structure Response where
  dark : ℝ
  bright : ℝ
  dark_mem : dark ∈ Set.Icc 0 1
  bright_mem : bright ∈ Set.Icc 0 1
  fourth : ℝ
  fourth_mem : fourth ∈ Set.Icc 0 1

noncomputable def prob (r : Response) (p : ℝ) : Fin 8 → ℝ :=
  let q := (1-p)*r.dark + p*r.bright
  ![r.dark/4, (1-r.dark)/4, r.bright/4, (1-r.bright)/4,
    q/4, (1-q)/4, r.fourth/4, (1-r.fourth)/4]

def score (p : ℝ) : Fin 8 → ℝ := ![-(1-p), 0, -p, 0, 1, 0, 0, 0]

lemma prob_nonneg (r : Response) {p : ℝ} (hp : p ∈ Set.Icc 0 1) (i : Fin 8) :
    0 ≤ prob r p i := by
  obtain ⟨hp₀, hp₁⟩ := hp
  obtain ⟨ha, ha'⟩ := r.dark_mem
  obtain ⟨hb, hb'⟩ := r.bright_mem
  obtain ⟨hc, hc'⟩ := r.fourth_mem
  have hx := add_nonneg (mul_nonneg (sub_nonneg.mpr hp₁) ha) (mul_nonneg hp₀ hb)
  have hy : (1-p) * r.dark + p * r.bright ≤ 1 := by
    nlinarith [mul_nonneg (sub_nonneg.mpr hp₁) (sub_nonneg.mpr ha'),
               mul_nonneg hp₀ (sub_nonneg.mpr hb')]
  fin_cases i <;> simp [prob] <;> linarith

lemma exp_chord {x : ℝ} (hx : x ∈ Set.Icc (-1) 1) (t : ℝ) :
    exp (t*x) ≤ (1-x)/2 * exp (-t) + (1+x)/2 * exp t := by
  obtain ⟨lo, hi⟩ := hx
  have h := convexOn_exp.2 (Set.mem_univ (-t)) (Set.mem_univ t)
    (show 0 ≤ (1-x)/2 by linarith) (show 0 ≤ (1+x)/2 by linarith) (by ring)
  simp only [smul_eq_mul] at h
  convert h using 2
  ring

lemma oneStep (r : Response) {p : ℝ} (hp : p ∈ Set.Icc 0 1) (t : ℝ) :
    ∑ i, prob r p i * exp (t * score p i) ≤ exp (t^2/2) := by
  calc
    _ ≤ ∑ i, prob r p i *
        ((1-score p i)/2 * exp (-t) + (1+score p i)/2 * exp t) := by
      apply Finset.sum_le_sum
      intro i _
      apply mul_le_mul_of_nonneg_left _ (prob_nonneg r hp i)
      apply exp_chord
      obtain ⟨h₀, h₁⟩ := hp
      fin_cases i <;> simp [score, Set.mem_Icc] <;> constructor <;> linarith
    _ = cosh t := by rw [Real.cosh_eq]; simp [Fin.sum_univ_succ, prob, score]; ring
    _ ≤ _ := Real.cosh_le_exp_half_sq t

/-- Expected payoff at the leaves. Each node can select both response
endpoints after observing the complete path to that node. -/
noncomputable def expect (r : List (Fin 8) → Response) (p : ℝ) (pay : ℝ → ℝ) :
    ℕ → List (Fin 8) → ℝ → ℝ
  | 0, _, s => pay s
  | n+1, h, s => ∑ i, prob (r h) p i * expect r p pay n (i::h) (s + score p i)

lemma exp_tree (r : List (Fin 8) → Response) {p : ℝ} (hp : p ∈ Set.Icc 0 1) (t : ℝ)
    (n : ℕ) (h : List (Fin 8)) (s : ℝ) :
    expect r p (fun x => exp (t*x)) n h s ≤ exp (t*s + n*(t^2/2)) := by
  induction n generalizing h s with
  | zero => simp [expect]
  | succ n ih =>
    calc
      _ ≤ ∑ i, prob (r h) p i * exp (t*(s + score p i) + n*(t^2/2)) :=
        Finset.sum_le_sum fun i _ => mul_le_mul_of_nonneg_left (ih _ _) (prob_nonneg _ hp _)
      _ = exp (t*s + n*(t^2/2)) * ∑ i, prob (r h) p i * exp (t*score p i) := by
        rw [Finset.mul_sum]; apply Finset.sum_congr rfl; intro i _
        rw [show t*(s+score p i)+n*(t^2/2) = (t*s+n*(t^2/2))+t*score p i by ring, exp_add]
        ring
      _ ≤ exp (t*s + n*(t^2/2)) * exp (t^2/2) :=
        mul_le_mul_of_nonneg_left (oneStep _ hp _) (exp_pos _).le
      _ = _ := by rw [← exp_add, Nat.cast_succ]; congr 1; ring

lemma expect_dominate (r : List (Fin 8) → Response) {p : ℝ} (hp : p ∈ Set.Icc 0 1)
    {f g : ℝ → ℝ} {a b : ℝ} (bound : ∀ x, f x ≤ a * g x + b * g (-x))
    (n : ℕ) (h : List (Fin 8)) (s : ℝ) :
    expect r p f n h s ≤ a * expect r p g n h s +
      b * expect r p (fun x => g (-x)) n h s := by
  induction n generalizing h s with
  | zero => exact bound s
  | succ n ih =>
    simp only [expect, Finset.mul_sum, ← Finset.sum_add_distrib]
    apply Finset.sum_le_sum
    intro i _
    exact (mul_le_mul_of_nonneg_left (ih (i::h) (s+score p i)) (prob_nonneg _ hp i)).trans_eq (by ring)

/-- Two-sided finite-sample deviation of the actual predictable receiver,
uniformly over the choice of response at every node. -/
theorem twoSided (r : List (Fin 8) → Response) {p ε : ℝ} (hp : p ∈ Set.Icc 0 1)
    (he : 0 ≤ ε) (n : ℕ) :
    expect r p (fun x => if (n:ℝ)*ε ≤ |x| then 1 else 0) n [] 0 ≤
      2 * exp (-(n:ℝ)*ε^2/2) := by
  let a := exp (-(n:ℝ)*ε^2)
  have hbound (x : ℝ) : (if (n:ℝ)*ε ≤ |x| then (1:ℝ) else 0) ≤
      a * exp (ε*x) + a * exp (ε*(-x)) := by
    dsimp [a]
    split_ifs with hx
    · rw [← exp_add, ← exp_add]
      rcases le_abs.mp hx with h | h
      · have h' := exp_le_exp.mpr (show 0 ≤ -(n:ℝ)*ε^2+ε*x by nlinarith [mul_le_mul_of_nonneg_left h he])
        rw [exp_zero] at h'; linarith [exp_pos (-(n:ℝ)*ε^2+ε* -x)]
      · have h' := exp_le_exp.mpr (show 0 ≤ -(n:ℝ)*ε^2+ε* -x by nlinarith [mul_le_mul_of_nonneg_left h he])
        rw [exp_zero] at h'; linarith [exp_pos (-(n:ℝ)*ε^2+ε*x)]
    · positivity
  have hpos := exp_tree r hp ε n [] 0
  have hneg := exp_tree r hp (-ε) n [] 0
  have ha : 0 ≤ a := (exp_pos _).le
  simp only [mul_zero, zero_add, neg_sq] at hpos hneg
  have step := expect_dominate r hp hbound n [] 0
  simp only [mul_neg, ← neg_mul] at step
  refine step.trans ((add_le_add (mul_le_mul_of_nonneg_left hpos ha)
    (mul_le_mul_of_nonneg_left hneg ha)).trans_eq ?_)
  dsimp [a]
  rw [← exp_add, ← two_mul]
  congr 2
  ring

/-- Sum of the marginal failure bounds for the enrolled score family.
The finite union bound requires no cross-bank independence. -/
theorem jointBudget {n m α : ℝ} (hn : 0 < n) (hm : 0 < m) (ha : 0 < α)
    (hr : 1 ≤ 2*m/α) :
    m * (2 * exp (-n * (sqrt (2 * log (2*m/α) / n))^2 / 2)) = α := by
  have hlog := Real.log_nonneg hr
  rw [Real.sq_sqrt (by positivity)]
  rw [show -n * (2 * log (2*m/α) / n) / 2 = -log (2*m/α) by field_simp,
    Real.exp_neg, Real.exp_log (by positivity)]
  field_simp

end IntervalEcho.FiniteReceiver
