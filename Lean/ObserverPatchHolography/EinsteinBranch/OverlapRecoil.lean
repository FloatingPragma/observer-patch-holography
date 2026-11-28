import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

/-! A finite comparison seam and the flux of its ordinary-energy tickets.
A common fair flip masks the sign of both binary records. The apparatus can
install the second load in either orientation. Currents obey node balance in
the interior. The exterior may absorb arbitrary current. -/

namespace OPH.EinsteinBranch.OverlapRecoil

open scoped BigOperators

/-- Independent sharp signs with means b and c differ with this probability. -/
noncomputable def unlike (b c : ℝ) : ℝ := ((1+b)/2)*((1-c)/2) + ((1-b)/2)*((1+c)/2)

theorem unlike_eq (b c : ℝ) : unlike b c = (1-b*c)/2 := by unfold unlike; ring

/-- Either single head remains fair after the common flip, even for a biased word. -/
theorem fair_marginal (b : ℝ) : (1/2 : ℝ)*((1+b)/2) + (1/2)*((1-b)/2) = 1/2 := by ring

/-- Each of the ordered unlike outputs occurs half as often as an unlike pair.
The two unlike entries of the density matrix therefore commute with the
native hopping term. This identity holds for any hopping matrix element. -/
theorem exchange_commutes (w hop : ℝ) : (w/2)*hop - hop*(w/2) = 0 := by ring

/-- The unlike load has equal symmetric and antisymmetric weights. -/
noncomputable def seamEnergy (b c : ℝ) : ℝ := 8*(unlike b c/2) + 9*(unlike b c/2)

noncomputable def ticket (b c : ℝ) : ℝ := seamEnergy b c / 512

theorem ordinary_calibration : ((8 : ℝ)/2 + 9/2)/512 = 17/1024 := by norm_num

/-- For native 0/1 words the energy is 17/2 times the mismatch
probability. At fixed single-head margins, erasing two original heads
replaces the joint occupied probability by the product of the margins.
Both routes must deliver to retain their original covariance. -/
theorem original_head_energy (left right both hL hR : ℝ) :
    (17/2 : ℝ)*(left+right-2*(left*right)) -
    (17/2 : ℝ)*(left+right-2*(left*right+hL*hR*(both-left*right))) =
    17*hL*hR*(both-left*right) := by ring

/-- Independent local binary readers give an exact signed covariance gain.
Here a0,a1 (and b0,b1) are probabilities of one at each local output. -/
theorem local_reader_covariance (l r j a0 a1 b0 b1 : ℝ) :
    (a0*b0*(1-l-r+j) + a1*b0*(l-j) + a0*b1*(r-j) + a1*b1*j)
    - ((1-l)*a0+l*a1)*((1-r)*b0+r*b1) =
    (a1-a0)*(b1-b0)*(j-l*r) := by ring

/-- Reversing both calibrations or adding independent identical read noise
preserves strict attraction for a positively correlated pair. -/
theorem matched_contact_sign (a0 a1 cov hL hR : ℝ)
    (sensitive : a1 ≠ a0) (hc : 0 < cov) (hl : 0 < hL) (hr : 0 < hR) :
    0 < 17*hL*hR*(a1-a0)^2*cov := by
  have hs : 0 < (a1-a0)^2 := sq_pos_of_ne_zero (sub_ne_zero.mpr sensitive)
  positivity

/-- A positive-temperature canonical contact has equal populations on
the two zero-energy states and weights 0<u,v<1 on the unlike levels.
Sharp reading preserves this sign after destroying all phase information. -/
theorem passive_sharp_energy (u v : ℝ) (hu : 0 < u) (hv : 0 < v)
    (hop : u < 1) (hodd : v < 1) :
    0 < (17/4 : ℝ) - (17/2)*(u+v)/(2+u+v) := by
  have h : 0 < 2+u+v := by positivity
  apply sub_pos.mpr
  rw [div_lt_iff₀ h]
  linarith

theorem erased_minus_reference (b s : ℝ) :
    seamEnergy 0 0 - seamEnergy b (s*b) = 17/4*s*b^2 := by
  simp only [seamEnergy, unlike_eq]
  ring

theorem ticket_difference (b s : ℝ) :
    ticket 0 0 - ticket b (s*b) = (17/2048)*s*b^2 := by
  unfold ticket
  rw [← sub_div, erased_minus_reference]
  ring

theorem energy_residual (b c d e : ℝ) :
    ticket b c - ticket d e - (seamEnergy b c - seamEnergy d e)/512 = 0 := by
  unfold ticket
  ring

/-- Probabilities of the four sharp states after applying the fair common flip. -/
noncomputable def population (v : ℝ) (x : Bool × Bool) : ℝ :=
  if x.1 = x.2 then (1+v)/4 else (1-v)/4

/-- Every off-diagonal entry of the native hopping term connects a sharp
state to its swap. Thus every entry of its commutator with the complete
masked diagonal state vanishes. The contrast term is itself diagonal. -/
theorem native_stationarity (v : ℝ) (x y : Bool × Bool) (hop : ℝ)
    (support : x = y ∨ x = (y.2, y.1) ∨ hop = 0) :
    hop * (population v y - population v x) = 0 := by
  rcases support with rfl | rfl | rfl
  · simp
  · simp [population, eq_comm]
  · simp

/-- Rates for every fixed downstream measurement are affine in the prepared
state. The erased state is the midpoint of the two port orientations. -/
theorem opposite_effects (v : ℝ) (read : Bool × Bool → ℝ) :
    (∑ x, read x * (population 0 x - population v x)) =
    -(∑ x, read x * (population 0 x - population (-v) x)) := by
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro x _
  by_cases h : x.1 = x.2 <;> simp [population, h] <;> ring

/-- If a closed energy budget and a single additive ticket rate cover both
the reader and its supply, an energy transfer changes neither the sum of
the two rates nor the flux of their joint enclosure. -/
theorem closed_energy_budget (reader supply reader' supply' : ℝ)
    (conserved : reader + supply = reader' + supply') :
    (reader'/512 - reader/512) + (supply'/512 - supply/512) = 0 := by linarith

theorem independent_erasure (b s e : ℝ) :
    seamEnergy ((1-e)*b) (s*((1-e)*b)) - seamEnergy b (s*b) =
    (17/4)*s*b^2*(2*e-e^2) := by
  simp only [seamEnergy, unlike_eq]
  ring

theorem shared_erasure (b s e : ℝ) :
    (1-e)*seamEnergy b (s*b) + e*seamEnergy 0 0 - seamEnergy b (s*b) =
    (17/4)*s*b^2*e := by
  simp only [seamEnergy, unlike_eq]
  ring

/-- Count visits to a target in disjoint first-hit cohorts. Each hit starts
the same expected subsequent visit count, including the first appointment. -/
theorem first_hit_factor {T : Type*} (times : Finset T) (p : T → ℝ)
    (diagonal allVisits : ℝ) (hd : 0 < diagonal)
    (restart : allVisits = ∑ t ∈ times, p t * diagonal) :
    allVisits / diagonal = ∑ t ∈ times, p t := by
  rw [restart, ← Finset.sum_mul]
  exact mul_div_cancel_right₀ _ (ne_of_gt hd)

/-- Antisymmetry cancels every internal bond, including self currents. -/
theorem internal_flux {V : Type*} (S : Finset V) (J : V → V → ℝ)
    (antisym : ∀ x y, J x y = -J y x) :
    ∑ x ∈ S, ∑ y ∈ S, J x y = 0 := by
  have h : (∑ x ∈ S, ∑ y ∈ S, J x y) = -(∑ x ∈ S, ∑ y ∈ S, J x y) := by
    calc
      _ = ∑ x ∈ S, ∑ y ∈ S, -J y x := by
        exact Finset.sum_congr rfl fun x _ =>
          Finset.sum_congr rfl fun y _ => antisym x y
      _ = -(∑ y ∈ S, ∑ x ∈ S, J y x) := by
        simp only [Finset.sum_neg_distrib]
        rw [Finset.sum_comm]
  linarith

/-- Includes edges to the drain. Only interior vertices require balance. -/
theorem finite_flux {V : Type*} [Fintype V] [DecidableEq V] (S : Finset V)
    (J : V → V → ℝ) (rate : V → ℝ)
    (antisym : ∀ x y, J x y = -J y x)
    (balance : ∀ x ∈ S, ∑ y, J x y = rate x) :
    ∑ x ∈ S, ∑ y ∈ Finset.univ \ S, J x y = ∑ x ∈ S, rate x := by
  have split (x : V) :
      (∑ y ∈ S, J x y) + (∑ y ∈ Finset.univ \ S, J x y) = ∑ y, J x y := by
    exact Finset.sum_add_sum_compl S (J x)
  have total : (∑ x ∈ S, ∑ y ∈ S, J x y) +
      (∑ x ∈ S, ∑ y ∈ Finset.univ \ S, J x y) =
      ∑ x ∈ S, rate x := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun x hx => (split x).trans (balance x hx)
  rw [internal_flux S J antisym, zero_add] at total
  exact total

theorem shell_flux {V : Type*} [Fintype V] [DecidableEq V] (S : Finset V)
    (J : V → V → ℝ) (b : V → ℝ) (s : ℝ)
    (antisym : ∀ x y, J x y = -J y x)
    (balance : ∀ x ∈ S, ∑ y, J x y = ticket 0 0 - ticket (b x) (s*b x)) :
    ∑ x ∈ S, ∑ y ∈ Finset.univ \ S, J x y = (17/2048)*s*∑ x ∈ S, (b x)^2 := by
  rw [finite_flux S J _ antisym balance, Finset.mul_sum]
  exact Finset.sum_congr rfl fun x _ => ticket_difference (b x) s

theorem original_head_shell {V : Type*} [Fintype V] [DecidableEq V] (S : Finset V)
    (J : V → V → ℝ) (hit : V → ℝ) (cov faces : ℝ)
    (hc : 0 < cov) (hf : 0 < faces) (seen : ∃ x ∈ S, hit x ≠ 0)
    (antisym : ∀ x y, J x y = -J y x)
    (balance : ∀ x ∈ S, ∑ y, J x y = 17/512*cov*(hit x)^2) :
    0 < (8/5 : ℝ)*(∑ x ∈ S, ∑ y ∈ Finset.univ \ S, J x y)/faces := by
  rw [finite_flux S J _ antisym balance, ← Finset.mul_sum]
  obtain ⟨x, hx, hhit⟩ := seen
  have hp : 0 < ∑ y ∈ S, (hit y)^2 :=
    Finset.sum_pos' (fun y _ => sq_nonneg _) ⟨x, hx, sq_pos_of_ne_zero hhit⟩
  positivity

/-- The standard receiver reads 8/5 times outward current as inward recoil.
Reverse encoding changes the sign for exactly the same admission profile. -/
theorem mean_inward_sign (stock faces : ℝ) (h : 0 < stock) (hf : 0 < faces) :
    0 < (8/5 : ℝ)*(17/2048)*stock/faces ∧
    (8/5 : ℝ)*(17/2048)*(-1)*stock/faces < 0 := by
  constructor <;> first | positivity | apply div_neg_of_neg_of_pos <;> nlinarith

/-- On a shell containing an admitted address, the material recoil is strictly
inward with matching port orientation. Neither a radial power law nor a flux
sign enters as an assumption. Node balances can include drains outside S. -/
theorem inward_from_admission {V : Type*} [Fintype V] [DecidableEq V]
    (S : Finset V) (J : V → V → ℝ) (hit : V → ℝ) (z faces : ℝ)
    (hz : z ≠ 0) (hf : 0 < faces) (visible : ∃ x ∈ S, hit x ≠ 0)
    (antisym : ∀ x y, J x y = -J y x)
    (balance : ∀ x ∈ S, ∑ y, J x y =
      ticket 0 0 - ticket (z*hit x) (1*(z*hit x))) :
    0 < (8/5 : ℝ) * (∑ x ∈ S, ∑ y ∈ Finset.univ \ S, J x y) / faces := by
  rw [shell_flux S J (fun x => z * hit x) 1 antisym balance]
  obtain ⟨x, hx, hhit⟩ := visible
  have hp : 0 < ∑ y ∈ S, (z*hit y)^2 :=
    Finset.sum_pos' (fun y _ => sq_nonneg (z*hit y))
      ⟨x, hx, sq_pos_of_ne_zero (mul_ne_zero hz hhit)⟩
  positivity

end OPH.EinsteinBranch.OverlapRecoil
