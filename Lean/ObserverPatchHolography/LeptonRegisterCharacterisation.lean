import ObserverPatchHolography.KoideCirculant

/-!
# Conditional characterisations of charged-lepton register choices

The assumed branch weights for activity `p` are idle `1 - p`, agree `p / 2`,
and each oriented repair `p / 4`. For activities `pc` and `pe`, the squared
coefficients are `a^2 = (1 - pc) * (pe / 2) / 3` and
`b^2 = (pc / 4) * (1 - pe) / 3`. Equal activities give `(b/a)^2 = 1/2`
throughout the open unit interval. The signed Koide formula is reused from
`OPH.KoideCirculant`.

Each result is conditional on its stated flag rule or weight formula. These
results do not establish those weights for the 12-state seam object or the
432-state tick. They do not select a physical phase, identify physical masses,
or assert that OPH forces a construction choice.
-/

namespace OPH.LeptonRegisterCharacterisation

/-- Under the blank-pair and single-flag flip conditions, only exclusive OR remains.
This is a Boolean characterisation, not a proof of the weight formula for the
12-state object. It does not assert that OPH forces the combination rule. -/
theorem xor_unique (f : Bool → Bool → Bool)
    (h0 : f false false = false)
    (hx : ∀ x y, f (!x) y = !(f x y))
    (hy : ∀ x y, f x (!y) = !(f x y)) :
    f = xor := by
  have h10 : f true false = true := by simpa [h0] using hx false false
  have h01 : f false true = true := by simpa [h0] using hy false false
  have h11 : f true true = false := by simpa [h01] using hx false true
  funext x y
  cases x <;> cases y <;> simp [h0, h10, h01, h11]

/-- For the stated squared weights, equal activity in `(0,1)` gives `b^2/a^2 = 1/2`.
This is an identity for the weight formula, not a proof of the weights for the
12-state object. It does not assert that OPH forces equal activity. -/
theorem balance_of_equal_activity {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    ((p / 4) * (1 - p) / 3) / ((1 - p) * (p / 2) / 3) = 1 / 2 := by
  have h1 : (1 - p) ≠ 0 := by linarith
  have h2 : p ≠ 0 := ne_of_gt hp0
  field_simp
  ring

/-- For the stated squared weights and activities in `(0,1)`, balance holds exactly
when the activities are equal. This characterises the weight formula, not the
12-state object. It does not assert that OPH forces equal preparation. -/
theorem balance_iff_equal_activity {pc pe : ℝ}
    (_hc0 : 0 < pc) (hc1 : pc < 1) (he0 : 0 < pe) (he1 : pe < 1) :
    ((pc / 4) * (1 - pe) / 3) / ((1 - pc) * (pe / 2) / 3) = 1 / 2 ↔ pc = pe := by
  have hc : 0 < 1 - pc := by linarith
  rw [div_eq_iff (by positivity)]
  constructor
  · intro h
    nlinarith
  · intro h
    subst h
    ring

/-- If `a^2` and `b^2` have the equal-activity weights, their modulus ratio is balanced.
This uses the weight formula as a hypothesis; it does not establish the 12-state
object or assert that OPH forces these coefficients. -/
theorem coefficient_balance_of_equal_activity {a b p : ℝ}
    (hp0 : 0 < p) (hp1 : p < 1)
    (ha : a ^ 2 = (1 - p) * (p / 2) / 3)
    (hb : b ^ 2 = (p / 4) * (1 - p) / 3) :
    (b / a) ^ 2 = 1 / 2 := by
  rw [div_pow, ha, hb]
  exact balance_of_equal_activity hp0 hp1

/-- With the assumed equal-activity weights, the signed Koide quotient is `2/3`.
This uses `KoideCirculant.koide_formula`. It does not establish the weight formula
for the 12-state object, nonnegative roots, a physical mass reading, or a choice
forced by OPH. -/
theorem signed_koide_of_equal_activity {a b p : ℝ}
    (hp0 : 0 < p) (hp1 : p < 1)
    (ha : a ^ 2 = (1 - p) * (p / 2) / 3)
    (hb : b ^ 2 = (p / 4) * (1 - p) / 3) :
    (3 * a ^ 2 + 6 * b ^ 2) / (3 * a) ^ 2 = 2 / 3 := by
  have ha0 : a ≠ 0 := by
    have hp : 0 < 1 - p := by linarith
    have hpos : 0 < a ^ 2 := by rw [ha]; positivity
    exact fun h => by simp [h] at hpos
  rw [KoideCirculant.koide_formula ha0,
    coefficient_balance_of_equal_activity hp0 hp1 ha hb]
  norm_num

/-- For independent flags with probabilities `p/4` and `p/2`, the exclusive count
has the stated polynomial form. This is an identity for assumed probabilities,
not a proof of the weight formula for the 12-state object. It does not assert
that OPH forces independent flags or exclusive OR. -/
theorem exclusive_count (p : ℝ) :
    p / 4 + p / 2 - 2 * ((p / 4) * (p / 2)) = 3 * p / 4 - p ^ 2 / 4 := by
  ring

/-- Substitution of the assumed activity `4/n` gives the stated count for nonzero `n`.
This does not prove a port count, the weight formula for the 12-state object,
or that OPH forces `p = 4/n` or a carrier with twelve ports. -/
theorem exclusive_count_ports {n : ℝ} (hn : n ≠ 0) :
    3 * (4 / n) / 4 - (4 / n) ^ 2 / 4 = 3 / n - 4 / n ^ 2 := by
  field_simp

/-- Under the polynomial count formula on `[0,1]`, count `2/9` leaves only activity
`1/3`. This does not establish the weight formula for the 12-state object,
identify the count with a physical phase, or assert that OPH forces the activity. -/
theorem exclusive_count_eq_two_ninths_iff {p : ℝ} (_hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    3 * p / 4 - p ^ 2 / 4 = 2 / 9 ↔ p = 1 / 3 := by
  constructor
  · intro h
    have hfac : (3 * p - 1) * (3 * p - 8) = 0 := by nlinarith
    rcases mul_eq_zero.mp hfac with h1 | h2
    · linarith
    · linarith
  · intro h
    subst h
    norm_num

/-- The four disjoint port bins of the marked-edge reader furnish one instance
of the equal-activity weights. The orientation is part of the reader. -/
theorem marked_edge_branch_weights :
    (8 / 12 : ℝ) = 1 - 1 / 3 ∧ (2 / 12 : ℝ) = (1 / 3) / 2 ∧
    (1 / 12 : ℝ) = (1 / 3) / 4 ∧
    (1 / 12 : ℝ) + 2 / 12 - 2 * (1 / 12) * (2 / 12) = 2 / 9 := by
  norm_num

/-- The amplitude of the symmetric endpoint pair is sqrt(2) times either
tip amplitude. Three-character Fourier readout therefore has balanced
intensities in the positive-root chamber. -/
theorem marked_edge_intensity_quotient {x y z : ℝ}
    (hs : x + y + z = 0) (hq : x ^ 2 + y ^ 2 + z ^ 2 = 3 / 2) :
    let a := Real.sqrt 2 / 6
    let b := (1 / 6 : ℝ)
    ((a + 2*b*x)^2 + (a + 2*b*y)^2 + (a + 2*b*z)^2) /
      ((a + 2*b*x) + (a + 2*b*y) + (a + 2*b*z))^2 = 2 / 3 := by
  dsimp
  rw [KoideCirculant.root_square_sum hs hq, KoideCirculant.root_sum hs]
  have h2 := Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)
  have hn : Real.sqrt 2 ≠ 0 := by positivity
  field_simp
  nlinarith

/-- On the marked seam (0,1), the free-generator rows at ports 4, 6, 8 are
proportional to y, z, x+z on a word with endpoint amplitudes x and tip
amplitudes y,z. Vanishing leakage forces the zero word. -/
theorem marked_edge_no_free_mode {x y z : ℂ} (h4 : y = 0) (h6 : z = 0)
    (h8 : x + z = 0) : x = 0 ∧ y = 0 ∧ z = 0 := by
  exact ⟨by simpa [h6] using h8, h4, h6⟩

#print axioms marked_edge_branch_weights
#print axioms marked_edge_intensity_quotient
#print axioms marked_edge_no_free_mode
#print axioms xor_unique
#print axioms balance_of_equal_activity
#print axioms balance_iff_equal_activity
#print axioms coefficient_balance_of_equal_activity
#print axioms signed_koide_of_equal_activity
#print axioms exclusive_count
#print axioms exclusive_count_ports
#print axioms exclusive_count_eq_two_ninths_iff

end OPH.LeptonRegisterCharacterisation
