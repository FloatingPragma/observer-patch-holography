import Geometry.SourceReusableBus

/-!
# Arbitrary finite paired-rail chains

This is a conditional construction for the candidate scalar pair-mean law.
The rail embedding, initial preparation and controller are supplied.
Cells 0,...,n form the bus; cells greater than n are untouched archives.
Only cell 0 needs a physical reset rung. Each adjacent copy needs both rail
edges. A weighted bound gives polynomial sufficient cleanup work; it does
not remove the attenuation of a single forward read.
-/

set_option autoImplicit false

namespace OPH.SourceBusScaling
noncomputable section
open OPH.SourceEncodedMemory OPH.SourceReusableBus

def forward : ℕ → List (Instruction ℕ)
  | 0 => []
  | n+1 => forward n ++ [.copy n (n+1)]

def sweep (n : ℕ) : List (Instruction ℕ) := .clear 0 :: forward n

def carry (a : ℕ → ℝ) : ℕ → ℝ
  | 0 => a 0
  | j+1 => (carry a j + a (j+1))/2

theorem forward_later (n j : ℕ) (a : ℕ → ℝ) (h : n < j) :
    execute (forward n) a j = a j := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [forward, execute_append, execute, step]
    simp only [copied, show j ≠ n by omega, show j ≠ n+1 by omega, or_self,
      if_false]
    exact ih (by omega)

theorem forward_last (n : ℕ) (a : ℕ → ℝ) :
    execute (forward n) a n = carry a n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [forward, execute_append, execute, step, copied, or_true, if_true]
    rw [ih, forward_later n (n+1) a (by omega)]
    rfl

theorem carry_error (n : ℕ) (a : ℕ → ℝ) (reference E : ℝ)
    (h0 : |a 0-reference| ≤ E) (ha : ∀ j, 0 < j → j ≤ n → |a j| ≤ E) :
    |carry a n-reference/(2:ℝ)^n| ≤ E := by
  induction n with
  | zero => simpa [carry] using h0
  | succ n ih =>
    have hc := ih (fun j hj hk => ha j hj (by omega))
    have h := export_error (carry a n) (a (n+1)) (reference/(2:ℝ)^n) E E
      hc (ha (n+1) (by omega) (by omega))
    simpa [carry, pow_succ, div_div] using h

theorem forward_read_error (n : ℕ) (a : ℕ → ℝ) (reference E : ℝ)
    (h0 : |a 0-reference| ≤ E) (ha : ∀ j, 0 < j → j ≤ n → |a j| ≤ E) :
    |execute (forward n) a n-reference/(2:ℝ)^n| ≤ E := by
  rw [forward_last]
  exact carry_error n a reference E h0 ha

theorem forward_before (n j : ℕ) (a : ℕ → ℝ) (h : j < n) :
    execute (forward n) a j = carry a (j+1) := by
  induction n with
  | zero => omega
  | succ n ih =>
    simp only [forward, execute_append, execute, step]
    by_cases hj : j = n
    · subst j
      simp only [copied, true_or, if_true]
      rw [forward_last, forward_later n (n+1) a (by omega)]
      rfl
    · simp only [copied, hj, show j ≠ n+1 by omega, or_self, if_false]
      exact ih (by omega)

theorem forward_length (n : ℕ) : (compile (forward n)).length = 2*n := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp [forward, compile_append, compile, word, copyWord, ih]
    omega

theorem sweep_length (n : ℕ) : (compile (sweep n)).length = 1+2*n := by
  simp [sweep, compile, word, clearWord, forward_length]
  omega

theorem sweep_archives (n j : ℕ) (a : ℕ → ℝ) (h : n < j) :
    execute (sweep n) a j = a j := by
  change execute (forward n) (cleared 0 a) j = a j
  rw [forward_later n j _ h]
  simp [cleared, show j ≠ 0 by omega]

/-- Weight and reciprocal spectral margin. Subtraction is in ℝ. -/
def mass (n : ℕ) : ℝ := ((n:ℝ)+1)^2
def weight (n i : ℕ) : ℝ := mass n - ((n:ℝ)-(i:ℝ))^2
def rate (n : ℕ) : ℝ := 1-1/mass n

theorem mass_pos (n : ℕ) : 0 < mass n := by unfold mass; positivity

theorem mass_ge_four (n : ℕ) (hn : 1 ≤ n) : 4 ≤ mass n := by
  have h : (1:ℝ) ≤ n := by exact_mod_cast hn
  unfold mass
  nlinarith

theorem weight_zero (n : ℕ) : weight n 0 = 2*(n:ℝ)+1 := by
  simp [weight, mass]
  ring

theorem weight_bounds (n i : ℕ) (hi : i ≤ n) :
    2*(n:ℝ)+1 ≤ weight n i ∧ weight n i ≤ mass n := by
  have hi' : (i:ℝ) ≤ n := by exact_mod_cast hi
  have hn : (0:ℝ) ≤ n := Nat.cast_nonneg n
  have hj : (0:ℝ) ≤ i := Nat.cast_nonneg i
  dsimp [weight, mass]
  constructor <;> nlinarith [sq_nonneg ((n:ℝ)-(i:ℝ)),
    mul_nonneg hj (show 0 ≤ 2*(n:ℝ)-(i:ℝ) by linarith)]

theorem weight_mono (n i j : ℕ) (hi : i ≤ j) (hj : j ≤ n) :
    weight n i ≤ weight n j := by
  have hi' : (i:ℝ) ≤ j := by exact_mod_cast hi
  have hj' : (j:ℝ) ≤ n := by exact_mod_cast hj
  unfold weight
  nlinarith [mul_nonneg (show 0 ≤ (j:ℝ)-(i:ℝ) by linarith)
    (show 0 ≤ 2*(n:ℝ)-(i:ℝ)-(j:ℝ) by linarith)]

theorem rate_bounds (n : ℕ) (hn : 1 ≤ n) : 0 ≤ rate n ∧ rate n < 1 := by
  have hm := mass_ge_four n hn
  have hp := mass_pos n
  have hi : 0 < 1/mass n := by positivity
  have hu : 1/mass n ≤ 1 := (div_le_one hp).mpr (by linarith)
  unfold rate
  constructor <;> linarith

theorem first_weight (n : ℕ) :
    weight n 1 / 2 ≤ rate n * weight n 0 := by
  have hp := mass_pos n
  apply (mul_le_mul_iff_left₀ hp).mp
  have he : rate n * weight n 0 = (mass n-1)*weight n 0/mass n := by
    unfold rate
    field_simp
  rw [he, div_mul_cancel₀ _ (ne_of_gt hp)]
  dsimp [weight, mass]
  norm_num
  nlinarith [sq_nonneg (n:ℝ)]

theorem step_weight (n j : ℕ) (hn : 1 ≤ n) (hj : j+2 ≤ n) :
    (rate n * weight n j + weight n (j+2))/2 ≤ rate n * weight n (j+1) := by
  have hp := mass_pos n
  have hm := mass_ge_four n hn
  have hw := (weight_bounds n (j+2) hj).2
  have hc : weight n j + weight n (j+2) = 2*weight n (j+1)-2 := by
    simp [weight, Nat.cast_add, Nat.cast_one]
    ring
  have hr : rate n * mass n = mass n-1 := by
    unfold rate
    field_simp
  have hid :
      (2*rate n*weight n (j+1)-rate n*weight n j-weight n (j+2))*mass n =
        2*mass n-2-weight n (j+2) := by
    linear_combination (2*weight n (j+1)-weight n j)*hr - (mass n-1)*hc
  have : 0 ≤ (2*rate n*weight n (j+1)-rate n*weight n j-weight n (j+2))*mass n := by
    rw [hid]; linarith
  have := (nonneg_of_mul_nonneg_left this hp)
  linarith

def WeightedBound (n : ℕ) (E : ℝ) (a : ℕ → ℝ) : Prop :=
  ∀ i ≤ n, |a i| ≤ E * weight n i

theorem carry_bound (n k : ℕ) (a : ℕ → ℝ) (E : ℝ)
    (hn : 1 ≤ n) (hE : 0 ≤ E) (ha : WeightedBound n E a) (hk : k+1 ≤ n) :
    |carry (cleared 0 a) (k+1)| ≤ E * rate n * weight n k := by
  induction k with
  | zero =>
    have h := ha 1 hk
    have h' := abs_le.mp h
    have hw := mul_le_mul_of_nonneg_left (first_weight n) hE
    norm_num [carry, cleared]
    apply abs_le.mpr
    constructor <;> nlinarith
  | succ k ih =>
    have hc := abs_le.mp (ih (by omega))
    have ha' := abs_le.mp (ha (k+2) hk)
    have hw := mul_le_mul_of_nonneg_left (step_weight n k hn hk) hE
    change |(carry (cleared 0 a) (k+1) + cleared 0 a (k+2))/2| ≤ _
    simp only [cleared, show k+2 ≠ 0 by omega, if_false]
    apply abs_le.mpr
    constructor <;> nlinarith

theorem sweep_contraction (n : ℕ) (a : ℕ → ℝ) (E : ℝ)
    (hn : 1 ≤ n) (hE : 0 ≤ E) (ha : WeightedBound n E a) :
    WeightedBound n (E*rate n) (execute (sweep n) a) := by
  intro i hi
  change |execute (forward n) (cleared 0 a) i| ≤ _
  by_cases h : i < n
  · rw [forward_before n i _ h]
    exact carry_bound n i a E hn hE ha (by omega)
  · have hin : i = n := by omega
    subst i
    obtain ⟨k, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : n ≠ 0)
    rw [forward_last]
    exact (carry_bound (k+1) k a E hn hE ha (by omega)).trans
      (mul_le_mul_of_nonneg_left (weight_mono (k+1) k (k+1) (by omega) (by omega))
        (mul_nonneg hE (rate_bounds (k+1) hn).1))

def cleaning (n : ℕ) : ℕ → List (Instruction ℕ)
  | 0 => []
  | r+1 => sweep n ++ cleaning n r

theorem cleaning_length (n r : ℕ) :
    (compile (cleaning n r)).length = (1+2*n)*r := by
  induction r with
  | zero => rfl
  | succ r ih =>
    simp only [cleaning, compile_append, List.length_append, sweep_length, ih]
    ring

def readCycle (source n r : ℕ) : List (Instruction ℕ) :=
  [.copy source 0] ++ forward n ++ cleaning n r

theorem readCycle_length (source n r : ℕ) :
    (compile (readCycle source n r)).length = 2+2*n+(1+2*n)*r := by
  simp [readCycle, compile_append, compile, word, copyWord,
    forward_length, cleaning_length]
  omega

theorem cleaning_contraction (n r : ℕ) (a : ℕ → ℝ) (E : ℝ)
    (hn : 1 ≤ n) (hE : 0 ≤ E) (ha : WeightedBound n E a) :
    WeightedBound n (E*rate n^r) (execute (cleaning n r) a) := by
  induction r generalizing a E with
  | zero => simpa [cleaning, execute] using ha
  | succ r ih =>
    have hb := ih (execute (sweep n) a) (E*rate n)
      (mul_nonneg hE (rate_bounds n hn).1) (sweep_contraction n a E hn hE ha)
    simpa [cleaning, execute_append, pow_succ, mul_assoc, mul_comm, mul_left_comm] using hb

theorem cleanup_bound (n r i : ℕ) (a : ℕ → ℝ) (A : ℝ)
    (hn : 1 ≤ n) (hi : i ≤ n) (hA : 0 ≤ A)
    (ha : ∀ j ≤ n, |a j| ≤ A) :
    |execute (cleaning n r) a i| ≤ A*mass n/(2*(n:ℝ)+1)*rate n^r := by
  have hd : 0 < 2*(n:ℝ)+1 := by positivity
  have he : 0 ≤ A/(2*(n:ℝ)+1) := by positivity
  have hb : WeightedBound n (A/(2*(n:ℝ)+1)) a := by
    intro j hj
    apply (ha j hj).trans
    calc
      A = A/(2*(n:ℝ)+1)*(2*(n:ℝ)+1) := by field_simp
      _ ≤ _ := mul_le_mul_of_nonneg_left (weight_bounds n j hj).1 he
  have hc := cleaning_contraction n r a _ hn he hb i hi
  have hw := mul_le_mul_of_nonneg_left (weight_bounds n i hi).2
    (mul_nonneg he (pow_nonneg (rate_bounds n hn).1 r))
  calc
    _ ≤ (A/(2*(n:ℝ)+1)*rate n^r)*mass n := hc.trans hw
    _ = _ := by ring

theorem sweep_zero (a : ℕ → ℝ) : execute (sweep 0) a = cleared 0 a := rfl

theorem cleaning_archives (n r j : ℕ) (a : ℕ → ℝ) (h : n < j) :
    execute (cleaning n r) a j = a j := by
  induction r generalizing a with
  | zero => rfl
  | succ r ih =>
    rw [cleaning, execute_append, ih, sweep_archives n j a h]

theorem rate_power_bound (n r : ℕ) (hn : 1 ≤ n) :
    rate n^r ≤ mass n/(mass n+(r:ℝ)) := by
  have hp := mass_pos n
  have hr := (rate_bounds n hn).1
  induction r with
  | zero => simp [ne_of_gt hp]
  | succ r ih =>
    have hd : 0 < mass n+(r:ℝ) := by positivity
    have hd' : 0 < mass n+(r:ℝ)+1 := by positivity
    calc
      rate n^(r+1) = rate n^r*rate n := pow_succ _ _
      _ ≤ (mass n/(mass n+(r:ℝ)))*rate n := mul_le_mul_of_nonneg_right ih hr
      _ ≤ mass n/(mass n+(r+1:ℕ)) := by
        unfold rate
        push_cast
        rw [show mass n + ((r:ℝ)+1) = mass n+(r:ℝ)+1 by ring]
        field_simp
        nlinarith [show (0:ℝ) ≤ r from Nat.cast_nonneg r]

/-- One block of (n+1)^2 sweeps at least halves the weighted norm. -/
theorem block_rate (n : ℕ) (hn : 1 ≤ n) :
    rate n^((n+1)^2) ≤ (1/2:ℝ) := by
  have h := rate_power_bound n ((n+1)^2) hn
  have he : (((n+1)^2:ℕ):ℝ) = mass n := by simp [mass]
  rw [he] at h
  have hh : mass n/(mass n+mass n) = (1/2:ℝ) := by
    have hp := ne_of_gt (mass_pos n)
    field_simp
    ring
  simpa only [hh] using h

/-- Polynomial operation count in chain length and the requested binary
cleanup accuracy: (1+2*n)*(n+1)^2*k scalar means. -/
theorem cleanup_blocks (n k i : ℕ) (a : ℕ → ℝ) (A : ℝ)
    (hn : 1 ≤ n) (hi : i ≤ n) (hA : 0 ≤ A)
    (ha : ∀ j ≤ n, |a j| ≤ A) :
    |execute (cleaning n ((n+1)^2*k)) a i| ≤
      A*mass n/(2*(n:ℝ)+1)*(1/2:ℝ)^k := by
  have hc := cleanup_bound n ((n+1)^2*k) i a A hn hi hA ha
  rw [pow_mul] at hc
  exact hc.trans (mul_le_mul_of_nonneg_left
    (pow_le_pow_left₀ (pow_nonneg (rate_bounds n hn).1 _) (block_rate n hn) k)
    (div_nonneg (mul_nonneg hA (le_of_lt (mass_pos n))) (by positivity)))

end
end OPH.SourceBusScaling
