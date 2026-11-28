import Mathlib.Tactic

/-! Fixed-marginal source contrasts, calibrated native contact pilots, and
the finite reservoir needed for a compensated positive aperture. -/

namespace OPH.EinsteinBranch.ClosedOverlap

open scoped BigOperators

/-- Any linear read on two classical inputs is determined by its four
conditional outputs. Their connected contribution survives matching both
one-head marginals, also after a common quantum evolution. -/
theorem covariance_read (l r j a b c d : ℝ) :
    ((1-l)*(1-r)*a + (1-l)*r*b + l*(1-r)*c + l*r*d) -
    ((1-l-r+j)*a + (r-j)*b + (l-j)*c + j*d) =
    (j-l*r)*(b+c-a-d) := by ring

/-- In a closed quantum flight, a unitary commuting with its full energy
leaves the energy of every prepared input unchanged. -/
theorem closed_flight {A : Type*} [Monoid A] (U Uadj H : A)
    (unitary : Uadj * U = 1) (conserved : H * U = U * H) :
    Uadj * (H * U) = H := by rw [conserved, ← mul_assoc, unitary, one_mul]

/-- Measuring balanced empty/double against a symmetric one-load state costs
eight native units; antisymmetric against symmetric costs one. A masked
sharp record reference gives identical modular readings on both parities.
The two stationary pilot banks below agree on mean energy and each endpoint. -/
theorem matched_pilot (kEqual kUnlike : ℝ) :
    ((8:ℝ)/9*9 + (1/9)*0 = 8) ∧
    kUnlike - ((8/9)*kUnlike + (1/9)*kEqual) = (kUnlike-kEqual)/9 := by
  constructor <;> ring

theorem energy_match_all_pilots (kEqual kEven kOdd beta : ℝ)
    (ordinary : kEven-kEqual = beta*8) :
    kEven - ((8/9)*kOdd+(1/9)*kEqual) = 0 ↔ kOdd-kEven=beta := by
  constructor <;> intro h <;> linarith

/-- Applying the calibration from 0 to energy 8 to the masked log produces
a spurious 8/9 energy units (one in 576 tickets) in the equal-energy banks. -/
theorem masked_excess (kEqual kUnlike beta : ℝ) (hb : beta ≠ 0)
    (cal : kUnlike-kEqual=8*beta) :
    (kUnlike-((8/9)*kUnlike+(1/9)*kEqual))/beta = 8/9 := by
  apply (div_eq_iff hb).mpr
  linarith

/-- Independent ordinary shutters keep a preparation with probability s^e.
Normalizing the survivors gives a reference log fixed on every native level. -/
theorem shutter_level (s Z : ℝ) (hs : 0 < s) (hz : 0 < Z) (e : ℕ) :
    -Real.log (s^e/Z) = Real.log Z + (-Real.log s)*(e : ℝ) := by
  rw [Real.log_div (pow_ne_zero _ (ne_of_gt hs)) (ne_of_gt hz), Real.log_pow]
  ring

/-- This transfer holds simultaneously for all source changes on the measured
spectrum, including changes of parity and same-energy stationary banks. -/
theorem shutter_charge {A : Type*} (S : Finset A) (e : A → ℕ) (delta : A → ℝ)
    (total : ∑ i ∈ S, delta i = 0) (s Z : ℝ) (hs : 0 < s) (hz : 0 < Z) :
    (∑ i ∈ S, delta i * (-Real.log (s^(e i)/Z))) =
      (-Real.log s) * ∑ i ∈ S, delta i * (e i : ℝ) := by
  simp_rw [shutter_level s Z hs hz, mul_add]
  rw [Finset.sum_add_distrib, ← Finset.sum_mul, total, zero_mul, zero_add, Finset.mul_sum]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- The annealed sharp mismatch contrast in any reference diagonal on the
four energy eigenwords. Its ordinary calibrated residual equals one half
of the parity-temperature error times the change in mismatch probability. -/
theorem sharp_reference_residual (equal even odd d beta : ℝ)
    (cal : even-equal=8*beta) :
    d*((even+odd)/2-equal) - beta*(17/2*d) = d*(odd-even-beta)/2 := by
  rw [show equal = even-8*beta by linarith]
  ring

/-- For any positive aperture response, the relative multiplier in an outer
shell must explain all excess over the original outer-to-inner bound. -/
theorem aperture_gain (inner outer gainIn gainOut bound obs : ℝ)
    (hi : 0 < inner) (ho : 0 ≤ gainOut)
    (profile : outer ≤ (bound-1)*inner)
    (needed : (obs-1)*(gainIn*inner) ≤ gainOut*outer) :
    (obs-1)*gainIn ≤ gainOut*(bound-1) := by
  apply le_of_mul_le_mul_right (a := inner) _ hi
  calc
    (obs-1)*gainIn*inner ≤ gainOut*outer := by simpa [mul_assoc] using needed
    _ ≤ (gainOut*(bound-1))*inner := by
      simpa [mul_assoc] using mul_le_mul_of_nonneg_left profile ho

/-- Population changes in a closed binary contact bank. The primitive move
replaces 18 symmetric one-load states by 16 antisymmetric states and two
zero-energy states, one empty and one doubly occupied. -/
theorem population_resonance (a b c d : ℤ)
    (number : a+b+c+d=0) (load : b+c+2*d=0) (energy : 8*b+9*c=0) :
    ∃ t : ℤ, a=t ∧ b=(-18)*t ∧ c=16*t ∧ d=t := by
  exact ⟨d, by omega, by omega, by omega, rfl⟩

theorem no_small_population_change (a b c d : ℤ)
    (number : a+b+c+d=0) (load : b+c+2*d=0) (energy : 8*b+9*c=0)
    (boundLo : -18 < b) (boundHi : b < 18) :
    a=0 ∧ b=0 ∧ c=0 ∧ d=0 := by
  obtain ⟨t, ha, hb, hc, hd⟩ := population_resonance a b c d number load energy
  omega

/-- A positive enclosed excess requires an exterior debit of at least as much
available initial stock. The hypothesis is a common positive energy-to-mass
conversion before and after redistribution. -/
theorem outside_budget {V : Type*} [Fintype V] [DecidableEq V]
    (S : Finset V) (initial final : V → ℝ)
    (same : ∑ x, (final x-initial x)=0)
    (positive : ∀ x, 0 ≤ final x) :
    ∑ x ∈ S, (final x-initial x) ≤ ∑ x ∈ Sᶜ, initial x := by
  have split := Finset.sum_add_sum_compl S (fun x => final x-initial x)
  rw [same, Finset.sum_sub_distrib] at split
  have hn : 0 ≤ ∑ x ∈ Sᶜ, final x := Finset.sum_nonneg fun x _ => positive x
  simp only [Finset.sum_sub_distrib] at split ⊢
  linarith

/-- If only an exterior shell can donate, its *initial* stock bounds the
positive aperture excess for arbitrary internal redistribution. -/
theorem reachable_budget {V : Type*} [Fintype V] [DecidableEq V]
    (S donors : Finset V) (initial final : V → ℝ) (hd : donors ⊆ Sᶜ)
    (same : ∑ x, (final x-initial x)=0)
    (positive : ∀ x, 0 ≤ final x)
    (unchanged : ∀ x ∈ Sᶜ \ donors, final x=initial x) :
    ∑ x ∈ S, (final x-initial x) ≤ ∑ x ∈ donors, initial x := by
  have split := Finset.sum_add_sum_compl S (fun x => final x-initial x)
  rw [same] at split
  have outside := Finset.sum_sdiff hd (f := fun x => final x - initial x)
  have hzero : ∑ x ∈ Sᶜ \ donors, (final x-initial x)=0 :=
    Finset.sum_eq_zero fun x hx => sub_eq_zero.mpr (unchanged x hx)
  rw [hzero, zero_add] at outside
  rw [← outside, Finset.sum_sub_distrib] at split
  have hn : 0 ≤ ∑ x ∈ donors, final x := Finset.sum_nonneg fun x _ => positive x
  simp only [Finset.sum_sub_distrib] at split ⊢
  linarith

end OPH.EinsteinBranch.ClosedOverlap
