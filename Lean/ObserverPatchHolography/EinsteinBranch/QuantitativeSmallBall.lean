import ObserverPatchHolography.EinsteinBranch.SmallBall
import ObserverPatchHolography.EinsteinBranch.Tensor
import Mathlib.Util.AssertNoSorry

/-!
# Quantitative small-ball tensor reconstruction

Ten explicit future unit timelike vectors determine a symmetric four-tensor.
The exact inverse has sharp coordinate infinity norm `17/8`. The small-ball
residual identity uses the same bulk entropy, area, stress and geometry inputs
in each frame. These theorems supply finite algebra and a conditional limit;
they construct no regional entropy or responsive metric source.
-/

namespace OPH.EinsteinBranch.Quantitative

noncomputable section

abbrev Coords := Fin 10 → ℝ

def matrix (x : Coords) : Mat 3 :=
  ![![x 0, x 1, x 2, x 3],
    ![x 1, x 4, x 5, x 6],
    ![x 2, x 5, x 7, x 8],
    ![x 3, x 6, x 8, x 9]]

def coordinates (B : Mat 3) : Coords :=
  ![B 0 0, B 0 1, B 0 2, B 0 3, B 1 1,
    B 1 2, B 1 3, B 2 2, B 2 3, B 3 3]

def probes : Fin 10 → V 3 :=
  ![![1, 0, 0, 0],
    ![5/3, 4/3, 0, 0], ![5/3, -4/3, 0, 0],
    ![5/3, 0, 4/3, 0], ![5/3, 0, -4/3, 0],
    ![5/3, 0, 0, 4/3], ![5/3, 0, 0, -4/3],
    ![3, 2, 2, 0], ![3, 2, 0, 2], ![3, 0, 2, 2]]

theorem probes_unit_future (r : Fin 10) :
    quadOf (eta 3) (probes r) = -1 ∧ 0 < probes r 0 := by
  fin_cases r <;>
    norm_num [probes, quadOf, bilinOf, eta, Fin.sum_univ_succ]

def design (x : Coords) : Coords :=
  ![x 0,
    25/9*x 0 + 40/9*x 1 + 16/9*x 4,
    25/9*x 0 - 40/9*x 1 + 16/9*x 4,
    25/9*x 0 + 40/9*x 2 + 16/9*x 7,
    25/9*x 0 - 40/9*x 2 + 16/9*x 7,
    25/9*x 0 + 40/9*x 3 + 16/9*x 9,
    25/9*x 0 - 40/9*x 3 + 16/9*x 9,
    9*x 0 + 12*x 1 + 12*x 2 + 4*x 4 + 8*x 5 + 4*x 7,
    9*x 0 + 12*x 1 + 12*x 3 + 4*x 4 + 8*x 6 + 4*x 9,
    9*x 0 + 12*x 2 + 12*x 3 + 4*x 7 + 8*x 8 + 4*x 9]

def decoder (q : Coords) : Coords :=
  ![q 0,
    9/80*(q 1-q 2), 9/80*(q 3-q 4), 9/80*(q 5-q 6),
    9/32*(q 1+q 2)-25/16*q 0,
    q 7/8+7/16*q 0-99/320*(q 1+q 3)+9/320*(q 2+q 4),
    q 8/8+7/16*q 0-99/320*(q 1+q 5)+9/320*(q 2+q 6),
    9/32*(q 3+q 4)-25/16*q 0,
    q 9/8+7/16*q 0-99/320*(q 3+q 5)+9/320*(q 4+q 6),
    9/32*(q 5+q 6)-25/16*q 0]

theorem matrix_coordinates (B : Mat 3) (hB : ∀ i j, B i j = B j i) :
    matrix (coordinates B) = B := by
  funext i j
  fin_cases i <;> fin_cases j <;> simp [matrix, coordinates]
  all_goals apply hB

theorem design_eq_charge (x : Coords) (r : Fin 10) :
    design x r = quadOf (matrix x) (probes r) := by
  fin_cases r <;>
    simp [design, probes, matrix, quadOf, bilinOf, Fin.sum_univ_succ] <;> ring

theorem decoder_design (x : Coords) : decoder (design x) = x := by
  funext i
  fin_cases i <;> simp [decoder, design] <;> ring

theorem design_decoder (q : Coords) : design (decoder q) = q := by
  funext i
  fin_cases i <;> simp [decoder, design] <;> ring

theorem decoder_bound (q : Coords) (ε : ℝ) (hε : 0 ≤ ε)
    (hq : ∀ r, |q r| ≤ ε) (i : Fin 10) :
    |decoder q i| ≤ 17/8*ε := by
  have h0 := abs_le.mp (hq 0)
  have h1 := abs_le.mp (hq 1)
  have h2 := abs_le.mp (hq 2)
  have h3 := abs_le.mp (hq 3)
  have h4 := abs_le.mp (hq 4)
  have h5 := abs_le.mp (hq 5)
  have h6 := abs_le.mp (hq 6)
  have h7 := abs_le.mp (hq 7)
  have h8 := abs_le.mp (hq 8)
  have h9 := abs_le.mp (hq 9)
  apply abs_le.mpr
  fin_cases i <;> simp [decoder] <;> constructor <;> linarith

def sharpCharges : Coords := ![-1, 1, 1, 0, 0, 0, 0, 0, 0, 0]

theorem decoder_bound_sharp :
    (∀ r, |sharpCharges r| ≤ 1) ∧ decoder sharpCharges 4 = 17/8 := by
  constructor
  · intro r; fin_cases r <;> norm_num [sharpCharges]
  · change (9/32:ℝ)*(1+1)-25/16*(-1) = 17/8
    norm_num

theorem tensor_bound (B : Mat 3) (hB : ∀ i j, B i j = B j i)
    (ε : ℝ) (hε : 0 ≤ ε)
    (hq : ∀ r, |quadOf B (probes r)| ≤ ε) (i j : Fin 4) :
    |B i j| ≤ 17/8*ε := by
  let x := coordinates B
  have hx : matrix x = B := matrix_coordinates B hB
  have hd : ∀ r, |design x r| ≤ ε := by
    intro r
    rw [design_eq_charge, hx]
    exact hq r
  have hb : ∀ k, |x k| ≤ 17/8*ε := by
    intro k
    have h := decoder_bound (design x) ε hε hd k
    rw [decoder_design] at h
    exact h
  rw [← hx]
  fin_cases i <;> fin_cases j <;> simp [matrix] <;> apply hb

def coefficient (ℓ : ℝ) : ℝ := 4*Real.pi*ℓ^4/15

theorem coefficient_pos (ℓ : ℝ) (hℓ : 0 < ℓ) : 0 < coefficient ℓ := by
  unfold coefficient
  positivity

/-- The scalar residual uses one entropy/area pair. `Λ` is the common
isotropic vacuum-reference term, not a direction-dependent fitted value. -/
theorem residual_identity (s a G T Λ η ℓ : ℝ) (hη : η ≠ 0) :
    η*coefficient ℓ*(G-Λ-(2*Real.pi/η)*T) =
      (s-2*Real.pi*coefficient ℓ*T) +
        η*(a+coefficient ℓ*(G-Λ))-(s+η*a) := by
  field_simp
  ring

theorem residual_bound (s a G T Λ η ℓ eb ea es : ℝ)
    (hη : 0 < η) (hℓ : 0 < ℓ)
    (hb : |s-2*Real.pi*coefficient ℓ*T| ≤ eb)
    (ha : |a+coefficient ℓ*(G-Λ)| ≤ ea)
    (hs : |s+η*a| ≤ es) :
    |G-Λ-(2*Real.pi/η)*T| ≤ (eb+η*ea+es)/(η*coefficient ℓ) := by
  have hc : 0 < η*coefficient ℓ := mul_pos hη (coefficient_pos ℓ hℓ)
  apply (le_div_iff₀ hc).mpr
  rw [mul_comm, ← abs_of_pos hc, ← abs_mul, residual_identity s a G T Λ η ℓ hη.ne']
  calc
    |(s-2*Real.pi*coefficient ℓ*T)+η*(a+coefficient ℓ*(G-Λ))-(s+η*a)|
      ≤ |s-2*Real.pi*coefficient ℓ*T| +
          |η*(a+coefficient ℓ*(G-Λ))| + |s+η*a| :=
        (abs_sub _ _).trans (add_le_add (abs_add_le _ _) (le_refl _))
    _ ≤ eb+η*ea+es := by
      rw [abs_mul, abs_of_pos hη]
      exact add_le_add (add_le_add hb (mul_le_mul_of_nonneg_left ha hη.le)) hs

def residual (G T : Mat 3) (Λ η : ℝ) : Mat 3 :=
  fun i j => G i j+Λ*EinsteinBranch.eta 3 i j-(2*Real.pi/η)*T i j

theorem residual_charge (G T : Mat 3) (Λ η : ℝ) (r : Fin 10) :
    quadOf (residual G T Λ η) (probes r) =
      quadOf G (probes r)-Λ-(2*Real.pi/η)*quadOf T (probes r) := by
  have he : (fun i j => G i j+Λ*EinsteinBranch.eta 3 i j) =
      fun i j => G i j-(-Λ)*EinsteinBranch.eta 3 i j := by
    funext i j; ring
  unfold residual
  rw [quadOf_sub_smul, he, quadOf_sub_smul, (probes_unit_future r).1]
  ring

/-- One common entropy/area family supplies every tensor charge. There is
no separate null-balance hypothesis in this conclusion. -/
theorem small_ball_tensor_bound (G T : Mat 3)
    (hG : ∀ i j, G i j = G j i) (hT : ∀ i j, T i j = T j i)
    (s a : Fin 10 → ℝ) (Λ η ℓ eb ea es : ℝ)
    (hη : 0 < η) (hℓ : 0 < ℓ) (heb : 0 ≤ eb) (hea : 0 ≤ ea) (hes : 0 ≤ es)
    (hb : ∀ r, |s r-2*Real.pi*coefficient ℓ*quadOf T (probes r)| ≤ eb)
    (ha : ∀ r, |a r+coefficient ℓ*(quadOf G (probes r)-Λ)| ≤ ea)
    (hs : ∀ r, |s r+η*a r| ≤ es) (i j : Fin 4) :
    |residual G T Λ η i j| ≤
      17/8*((eb+η*ea+es)/(η*coefficient ℓ)) := by
  apply tensor_bound
  · intro k l
    simp only [residual]
    rw [hG k l, hT k l, eta_symm k l]
  · exact div_nonneg (by positivity) (mul_pos hη (coefficient_pos ℓ hℓ)).le
  · intro r
    rw [residual_charge]
    exact residual_bound (s r) (a r) (quadOf G (probes r))
      (quadOf T (probes r)) Λ η ℓ eb ea es hη hℓ (hb r) (ha r) (hs r)

theorem normalized_error_identity (eb ea es η ℓ : ℝ)
    (hη : η ≠ 0) (hℓ : ℓ ≠ 0) :
    (eb+η*ea+es)/(η*coefficient ℓ) =
      15/(4*Real.pi)*((eb/ℓ^4)/η+ea/ℓ^4+(es/ℓ^4)/η) := by
  unfold coefficient
  field_simp

/-- Each error is measured on the same shrinking family. The geometric
scale and positive area coefficient cannot be canceled at unrelated indices. -/
theorem normalized_error_tendsto (eb ea es η ℓ : ℕ → ℝ) (ηstar : ℝ)
    (hη : ∀ n, η n ≠ 0) (hℓ : ∀ n, ℓ n ≠ 0) (hstar : ηstar ≠ 0)
    (heta : Filter.Tendsto η Filter.atTop (nhds ηstar))
    (hb : Filter.Tendsto (fun n => eb n/(ℓ n)^4) Filter.atTop (nhds 0))
    (ha : Filter.Tendsto (fun n => ea n/(ℓ n)^4) Filter.atTop (nhds 0))
    (hs : Filter.Tendsto (fun n => es n/(ℓ n)^4) Filter.atTop (nhds 0)) :
    Filter.Tendsto (fun n => (eb n+η n*ea n+es n)/(η n*coefficient (ℓ n)))
      Filter.atTop (nhds 0) := by
  have h := (((hb.div heta hstar).add ha).add (hs.div heta hstar)).const_mul
    (15/(4*Real.pi))
  simp only [zero_div, zero_add, mul_zero] at h
  convert h using 1
  funext n
  exact normalized_error_identity _ _ _ _ _ (hη n) (hℓ n)

/-- A bound tending to zero forces every tensor component to tend to zero.
The small-ball theorem supplies `ε`; its required scaling remains explicit. -/
theorem tensor_tendsto_zero (B : ℕ → Mat 3)
    (hB : ∀ n i j, B n i j = B n j i)
    (ε : ℕ → ℝ) (hε : ∀ n, 0 ≤ ε n)
    (hq : ∀ n r, |quadOf (B n) (probes r)| ≤ ε n)
    (hlim : Filter.Tendsto ε Filter.atTop (nhds 0)) (i j : Fin 4) :
    Filter.Tendsto (fun n => B n i j) Filter.atTop (nhds 0) := by
  have hb : ∀ n, |B n i j| ≤ 17/8*ε n :=
    fun n => tensor_bound (B n) (hB n) (ε n) (hε n) (hq n) i j
  have he : Filter.Tendsto (fun n => (17/8:ℝ)*ε n)
      Filter.atTop (nhds 0) := by
    simpa using hlim.const_mul (17/8:ℝ)
  have habs : Filter.Tendsto (fun n => |B n i j|) Filter.atTop (nhds 0) :=
    squeeze_zero (fun n => abs_nonneg _) hb he
  exact tendsto_zero_iff_norm_tendsto_zero.mpr habs

/-- A common cosmological coefficient is constant under the stated
conservation identities and one spacetime-independent coupling. This uses
the exact finite-difference Ward/Bianchi interface of `Tensor.lean`. -/
theorem lambda_constant {P : Type*} (step : Fin 4 → P → P) (base : P)
    (connected : ∀ q, SymmReachable step base q)
    (G T : P → Mat 3) (Λ : P → ℝ) (κ : ℝ)
    (heq : ∀ p i j, G p i j+Λ p*eta 3 i j = κ*T p i j)
    (bianchi : ∀ p j, ddiv step G p j = 0)
    (ward : ∀ p j, ddiv step T p j = 0) :
    ∀ p, Λ p = Λ base := by
  have heq' : ∀ p i j, G p i j = κ*T p i j+(-Λ p)*eta 3 i j := by
    intro p i j
    have h := heq p i j
    linarith
  have hstep := step_invariant_of_divergence_free step G T (fun p => -Λ p) κ
    heq' bianchi ward
  intro p
  have h := lam_eq_of_symmReachable step (fun p => -Λ p) hstep (connected p)
  linarith

assert_no_sorry probes_unit_future
assert_no_sorry decoder_design
assert_no_sorry design_decoder
assert_no_sorry decoder_bound
assert_no_sorry decoder_bound_sharp
assert_no_sorry tensor_bound
assert_no_sorry residual_identity
assert_no_sorry residual_bound
assert_no_sorry small_ball_tensor_bound
assert_no_sorry normalized_error_tendsto
assert_no_sorry tensor_tendsto_zero
assert_no_sorry lambda_constant

#print axioms tensor_bound
#print axioms residual_bound
#print axioms small_ball_tensor_bound
#print axioms normalized_error_tendsto
#print axioms tensor_tendsto_zero
#print axioms lambda_constant

end
end OPH.EinsteinBranch.Quantitative
