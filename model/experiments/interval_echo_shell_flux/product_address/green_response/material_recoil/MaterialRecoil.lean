import Mathlib.Tactic
import Mathlib.RingTheory.RootsOfUnity.Complex
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv

namespace MaterialRecoil

/-- Averaging three independent fair contact switches. Signs describe the
    coherence of one test load in the left/right material basis. -/
theorem switch_average (x y z : ℝ) :
    (1 + x + y + z + x*y + x*z + y*z + x*y*z)/8 =
    ((1+x)/2)*((1+y)/2)*((1+z)/2) := by ring

theorem active_phase_filter (a b c : Bool) :
    let sign : Bool → ℝ := fun hit => if hit then 1 else -1
    ((1+sign a)/2)*((1+sign b)/2)*((1+sign c)/2) =
      if a && b && c then 1 else 0 := by
  cases a <;> cases b <;> cases c <;> norm_num

/-- The free flight time is a half-period. Equal mirror weights divide
    the native displacement response by two. -/
theorem current_to_load (near far : ℝ) :
    let d : ℝ := 65/2
    d/(2*(1+d^2))*(far-near) = -(520/4229)*((near-far)/8) := by
  dsimp
  ring

/-- The two bare diagonal energies on a board of arbitrary native degree
    differ by one common constant. The off-diagonal swap is shared. -/
theorem bare_shift (k deg : ℝ) :
    ((6-deg)*k*2^2 + deg*(k*2^2+1/2)) -
    ((6-deg)*k + deg*(k+1/2)) = 18*k := by ring

namespace CommonReceiver

/-- Twice the contrast of an occupied or vacant sensing pad against a blank.
    The probe is load one when `two` is false and load two when it is true. -/
def twiceContrast (k : ℕ) (two : Bool) (s : Fin 3) : ℕ :=
  let q := if two then 2 else 1
  if s.val = 0 then 0 else 4*k*q*s.val + 1 + (if q = s.val then 1 else 0)

private def contact (k q p : ℕ) : ℤ :=
  2*(k:ℤ)*((q:ℤ)-(p:ℤ))^2 + (if q = p then 0 else 1)

/-- The twice-contrast equals the difference of signed splittings computed
    from four native contact energies, including when the sensed pad is vacant. -/
theorem contrast_from_contacts (k : ℕ) (two : Bool) (s : Fin 3) :
    let q := if two then 2 else 1
    (twiceContrast k two s : ℤ) =
      contact k q 0 - contact k 0 0 - contact k q s.val + contact k 0 s.val := by
  cases two <;> fin_cases s <;> norm_num [twiceContrast, contact] <;> ring

/-- The equally likely zero-, one- and two-period thirds. A reference
    coefficient of two represents minus one in the three-cycle. -/
noncomputable def window (z : ℂ) (k : ℕ) (two : Bool) (s p : Fin 3) : ℂ :=
  let n := twiceContrast k two s + 2*twiceContrast k two p
  (1 + z^(n%3) + z^((n*2)%3))/3

private theorem reduce_power {z : ℂ} (hz : z^3=1) (n : ℕ) : z^n = z^(n%3) := by
  conv_lhs => rw [← Nat.mod_add_div n 3, pow_add, pow_mul, hz, one_pow, mul_one]

/-- Directly averaging the sensed phase divided by the reference phase
    gives the reduced window, before any species-specific arithmetic. -/
theorem phase_average {z : ℂ} (hz : z^3 = 1) (k : ℕ) (two : Bool) (s p : Fin 3) :
    let f := z^(twiceContrast k two s) / z^(twiceContrast k two p)
    (1 + f + f^2)/3 = window z k two s p := by
  let a := twiceContrast k two s
  let b := twiceContrast k two p
  have h₀ : z ≠ 0 := by rintro rfl; norm_num at hz
  have signed : z^(a+2*b) = z^a/z^b := by
    apply (eq_div_iff (pow_ne_zero b h₀)).mpr
    rw [← pow_add, show a+2*b+b = a+3*b by omega, pow_add, pow_mul, hz, one_pow, mul_one]
  change (1 + z^a/z^b + (z^a/z^b)^2)/3 =
    (1 + z^((a+2*b)%3) + z^(((a+2*b)*2)%3))/3
  rw [← signed, ← pow_mul, reduce_power hz ((a+2*b)*2), reduce_power hz (a+2*b)]

/-- Both readings of the phase orientation are covered by an arbitrary
    primitive cubic root. The physical positive orientation is exhibited below. -/
theorem ternary_filter {z : ℂ} (hz : IsPrimitiveRoot z 3)
    (two : Bool) (s p : Fin 3) :
    window z 8 two s p = if s = p then 1 else 0 := by
  have h : (z-1)*(1+z+z^2) = 0 := by linear_combination hz.pow_eq_one
  have h₁ := (mul_eq_zero.mp h).resolve_left (sub_ne_zero.mpr (hz.ne_one (by decide)))
  have h₂ : 1+z^2+z=0 := by linear_combination h₁
  cases two <;> fin_cases s <;> fin_cases p <;> norm_num [window, twiceContrast, h₁, h₂]

theorem physical_root :
    IsPrimitiveRoot (Complex.exp (2 * Real.pi * Complex.I / 3)) 3 := by
  simpa using Complex.isPrimitiveRoot_exp 3 (by decide)

/-- At stiffness four the very same shutters pass every material label. -/
theorem stiffness_four (z : ℂ) (two : Bool) (s p : Fin 3) :
    window z 4 two s p = 1 := by
  cases two <;> fin_cases s <;> fin_cases p <;> norm_num [window, twiceContrast]

noncomputable def blind (z : ℂ) (two : Bool) (s : Fin 3) : ℂ :=
  (window z 8 two s 1 + window z 8 two s 2)/2

theorem blind_coordinate {z : ℂ} (hz : IsPrimitiveRoot z 3) (two : Bool) (s : Fin 3) :
    blind z two s = if s = 0 then 0 else 1/2 := by
  unfold blind
  rw [ternary_filter hz, ternary_filter hz]
  fin_cases s <;> norm_num <;> decide

theorem blind_three {z : ℂ} (hz : IsPrimitiveRoot z 3) (two : Bool) (a b c : Fin 3) :
    blind z two a * blind z two b * blind z two c =
    if a ≠ 0 ∧ b ≠ 0 ∧ c ≠ 0 then 1/8 else 0 := by
  simp only [blind_coordinate hz]
  by_cases ha : a = 0 <;> by_cases hb : b = 0 <;> by_cases hc : c = 0 <;>
    norm_num [ha, hb, hc]

/-- Degree completion leaves the same signed blank-link splitting for every
    mobile load. The off-diagonal exchange element is held at minus one half. -/
theorem blank_degree_detuning (k q : ℝ) :
    let E := fun d : ℝ => (6-d)*k*q^2 + d*(k*q^2+1/2)
    E 6 - E 5 = 1/2 := by dsimp; ring

/-- Mean displacement from the two equal face arms, with negative real
    coherence at the near arm and positive coherence at the far arm. -/
noncomputable def displacement (w near far t : ℝ) : ℝ :=
  (w*(far-near)/10) * (1-Real.cos (t*(Real.sqrt 5/2)))

theorem signed_face (w near far t : ℝ) :
    let D : ℝ := 1/2
    let a := D / (2*(1+D^2)) * (1-Real.cos (t*(Real.sqrt 5/2)))
    ((-w)*near*a + w*far*a)/2 = displacement w near far t := by
  dsimp [displacement]
  ring

theorem full_shutter (w near far : ℝ) :
    displacement w near far (Real.pi / (Real.sqrt 5/2)) =
    -w*(8/5)*((near-far)/8) := by
  have h : Real.sqrt 5 / 2 ≠ 0 := by positivity
  rw [displacement, div_mul_cancel₀ _ h, Real.cos_pi]
  ring

/-- Initial acceleration is an actual second derivative of the signed
    displacement. Weight one gives matched current transfer; weight one eighth
    gives blind transfer on occupied triples. -/
theorem initial_acceleration (w near far : ℝ) :
    deriv (deriv (displacement w near far)) 0 = -w*((near-far)/8) := by
  let A := w*(far-near)/10
  let v := Real.sqrt 5/2
  have speed (t : ℝ) :
      HasDerivAt (displacement w near far) (A*v*Real.sin (t*v)) t := by
    convert ((hasDerivAt_const t (1:ℝ)).sub (((hasDerivAt_id t).mul_const v).cos)).const_mul A
      using 1
    simp only [id_eq]
    ring
  rw [show deriv (displacement w near far) = (fun t => A*v*Real.sin (t*v))
    from funext (fun t => (speed t).deriv)]
  have accel := ((((hasDerivAt_id (0:ℝ)).mul_const v).sin).const_mul (A*v)).deriv
  convert accel using 1
  dsimp [A, v]
  rw [zero_mul, Real.cos_zero]
  nlinarith [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 5),
    congrArg (fun x : ℝ => w*(far-near)*x)
      (Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 5))]

/-- End-to-end blind transfer for either moving species and any occupied
    source triple. The gain is the measured-coordinate coherence product. -/
theorem occupied_blind_transfer {z : ℂ} (hz : IsPrimitiveRoot z 3) (two : Bool)
    (a b c : Fin 3) (occ : a ≠ 0 ∧ b ≠ 0 ∧ c ≠ 0) (near far : ℝ) :
    let w := (blind z two a * blind z two b * blind z two c).re
    let J := (near-far)/8
    displacement w near far (Real.pi / (Real.sqrt 5/2)) = -J/5 ∧
      deriv (deriv (displacement w near far)) 0 = -J/8 := by
  constructor <;> simp only [blind_three hz, if_pos occ, full_shutter, initial_acceleration] <;>
    norm_num [Complex.div_re] <;> ring

/-- A standard made of the enrolled guide materials passes all three
    matches for either probe. The two independent shutters read the same J. -/
theorem matched_transfer {z : ℂ} (hz : IsPrimitiveRoot z 3) (two : Bool)
    (a b c : Fin 3) (near far : ℝ) :
    let w := (window z 8 two a a * window z 8 two b b * window z 8 two c c).re
    let J := (near-far)/8
    displacement w near far (Real.pi / (Real.sqrt 5/2)) = -(8/5)*J ∧
      deriv (deriv (displacement w near far)) 0 = -J := by
  simp [ternary_filter hz, full_shutter, initial_acceleration]

end CommonReceiver
end MaterialRecoil
