import ObserverPatchHolography.DirectedSeamRepair
import PortFrameGram

namespace ObserverPatchHolography.DirectedSeamRepair

/-!
# Finite-path progress for directed seam repair

The local strict-descent theorem in `DirectedSeamRepair` does not by itself
exclude a state in which every adjacent load difference is at most one while
distant vertices differ by two or more.  On a connected carrier, a directed
repair can first transport a maximum load across edges whose endpoint loads
are `M` and `M - 1`, without changing the quadratic load, and then balance it
against the first endpoint at least two below `M`.

This file proves the arithmetic part of that argument for paths of length at
most three.  Three is the diameter of the twelve-vertex icosahedral carrier.
The graph-theoretic premise that the listed vertices form such a shortest
path remains explicit; these theorems do not select a stochastic schedule.
-/

/-- Every pair of ports in the concrete twelve-port carrier is connected by
a path with at most three seams. -/
theorem carrier_path_length_at_most_three :
    ∀ source target : Fin 12,
      source = target ∨
      OPH.PortFrameGram.adj source target = true ∨
      (∃ middle : Fin 12,
        OPH.PortFrameGram.adj source middle = true ∧
        OPH.PortFrameGram.adj middle target = true) ∨
      (∃ middle₁ middle₂ : Fin 12,
        OPH.PortFrameGram.adj source middle₁ = true ∧
        OPH.PortFrameGram.adj middle₁ middle₂ = true ∧
        OPH.PortFrameGram.adj middle₂ target = true) := by
  decide

/-! ## The global minimum in every signed total sector -/

/-- Total signed load on the concrete twelve-port carrier. -/
def totalLoad12 (x : Fin 12 → ℤ) : ℤ :=
  ∑ port, x port

/-- Global quadratic load on the concrete twelve-port carrier. -/
def quadraticLoad12 (x : Fin 12 → ℤ) : ℤ :=
  ∑ port, x port ^ 2

/-- The nonnegative integer defect from the two-level shell based at `q`. -/
def shellDefect12 (x : Fin 12 → ℤ) (q : ℤ) : ℤ :=
  ∑ port, (x port - q) * (x port - q - 1)

theorem integer_shell_term_nonnegative (z : ℤ) :
    0 ≤ z * (z - 1) := by
  by_cases hz : z ≤ 0
  · nlinarith
  · have : 1 ≤ z := by omega
    nlinarith

/-- The shell defect is nonnegative for arbitrary signed loads. -/
theorem shellDefect12_nonnegative (x : Fin 12 → ℤ) (q : ℤ) :
    0 ≤ shellDefect12 x q := by
  exact Finset.sum_nonneg fun port _ =>
    integer_shell_term_nonnegative (x port - q)

/-- Vanishing shell defect is equivalent to every port carrying one of the
two nearest integer loads. -/
theorem shellDefect12_eq_zero_iff (x : Fin 12 → ℤ) (q : ℤ) :
    shellDefect12 x q = 0 ↔
      ∀ port, x port = q ∨ x port = q + 1 := by
  constructor
  · intro hzero port
    have hall :
        ∀ port ∈ (Finset.univ : Finset (Fin 12)),
          (x port - q) * (x port - q - 1) = 0 :=
      (Finset.sum_eq_zero_iff_of_nonneg
        (fun port _ =>
          integer_shell_term_nonnegative (x port - q))).mp hzero
    have hport := hall port (Finset.mem_univ port)
    rcases mul_eq_zero.mp hport with hleft | hright
    · left
      omega
    · right
      omega
  · intro hshell
    apply (Finset.sum_eq_zero_iff_of_nonneg
      (fun port _ =>
        integer_shell_term_nonnegative (x port - q))).mpr
    intro port _
    rcases hshell port with hlow | hupp
    · simp [hlow]
    · simp [hupp]

/-- Exact completion-of-squares identity in a sector whose total is
`12*q + r`.  No nonnegativity premise is imposed on the port loads. -/
theorem quadraticLoad12_eq_shellEnergy_add_defect
    (x : Fin 12 → ℤ) (q r : ℤ)
    (htotal : totalLoad12 x = 12 * q + r) :
    quadraticLoad12 x =
      (12 - r) * q ^ 2 + r * (q + 1) ^ 2 + shellDefect12 x q := by
  simp only [quadraticLoad12, shellDefect12, totalLoad12] at htotal ⊢
  have hdefect :
      (∑ port : Fin 12, (x port - q) * (x port - q - 1)) =
        (∑ port : Fin 12, x port ^ 2) -
          (2 * q + 1) * (∑ port : Fin 12, x port) +
            12 * q * (q + 1) := by
    calc
      (∑ port : Fin 12, (x port - q) * (x port - q - 1)) =
          ∑ port : Fin 12,
            (x port ^ 2 - (2 * q + 1) * x port + q * (q + 1)) := by
              apply Finset.sum_congr rfl
              intro port _
              ring
      _ = (∑ port : Fin 12, x port ^ 2) -
            (2 * q + 1) * (∑ port : Fin 12, x port) +
              12 * q * (q + 1) := by
            rw [Finset.sum_add_distrib, Finset.sum_sub_distrib,
              ← Finset.mul_sum]
            simp
            ring
  rw [hdefect, htotal]
  ring

/-- The balanced two-level energy is a lower bound throughout every signed
sector represented as `12*q + r`. -/
theorem shellEnergy_is_global_lower_bound
    (x : Fin 12 → ℤ) (q r : ℤ)
    (htotal : totalLoad12 x = 12 * q + r) :
    (12 - r) * q ^ 2 + r * (q + 1) ^ 2 ≤ quadraticLoad12 x := by
  rw [quadraticLoad12_eq_shellEnergy_add_defect x q r htotal]
  exact le_add_of_nonneg_right (shellDefect12_nonnegative x q)

/-- Equality in the global lower bound occurs exactly on the balanced shell.
For Euclidean division `s = 12*q + r`, the total equation then fixes exactly
`r` ports at `q + 1`. -/
theorem shellEnergy_eq_iff_balanced
    (x : Fin 12 → ℤ) (q r : ℤ)
    (htotal : totalLoad12 x = 12 * q + r) :
    quadraticLoad12 x =
        (12 - r) * q ^ 2 + r * (q + 1) ^ 2 ↔
      ∀ port, x port = q ∨ x port = q + 1 := by
  rw [quadraticLoad12_eq_shellEnergy_add_defect x q r htotal]
  have hnonnegative := shellDefect12_nonnegative x q
  constructor
  · intro heq
    apply (shellDefect12_eq_zero_iff x q).mp
    nlinarith
  · intro hshell
    rw [(shellDefect12_eq_zero_iff x q).mpr hshell]
    ring

/-- Nearest balancing never increases the local quadratic load. -/
theorem quadraticLoad_forward_le (a b : ℤ) :
    quadraticLoad (balanceForward a b) ≤ quadraticLoad (a, b) := by
  by_cases hgap : 1 < (a - b) ^ 2
  · exact (quadraticLoad_forward_strict a b hgap).le
  · have hsmall : (a - b) ^ 2 ≤ 1 := le_of_not_gt hgap
    have hdrop := two_mul_quadratic_drop a b
    rcases output_gap_zero_or_one a b with hzero | hone
    · rw [hzero] at hdrop
      norm_num at hdrop
      nlinarith [sq_nonneg (a - b)]
    · have hne : a - b ≠ 0 := by
        intro hab
        have hsum := lower_add_upper a b
        omega
      have hpositive : 0 < (a - b) ^ 2 := sq_pos_of_ne_zero hne
      rw [hone] at hdrop
      norm_num at hdrop
      nlinarith

/-- The same nonincrease statement holds for the reverse direction. -/
theorem quadraticLoad_reverse_le (a b : ℤ) :
    quadraticLoad (balanceReverse a b) ≤ quadraticLoad (a, b) := by
  rw [quadraticLoad_reverse_eq_forward]
  exact quadraticLoad_forward_le a b

/-- A forward repair transports a maximum `M` through an endpoint carrying
either `M` or `M - 1`, without changing the two-load multiset. -/
theorem balanceForward_transports_high
    (M z : ℤ) (hz : z = M ∨ z = M - 1) :
    balanceForward M z = (z, M) := by
  rcases hz with rfl | rfl
  · apply Prod.ext <;> simp [balanceForward, lowerLoad, upperLoad] <;> omega
  · apply Prod.ext <;> simp [balanceForward, lowerLoad, upperLoad] <;> omega

/-- A one-edge path from `M` to a load at most `M - 2` strictly lowers the
quadratic load. -/
theorem path_one_forward_strict
    (M endpoint : ℤ) (hend : endpoint ≤ M - 2) :
    quadraticLoad (balanceForward M endpoint) <
      quadraticLoad (M, endpoint) := by
  apply quadraticLoad_forward_strict
  nlinarith [sq_nonneg (M - endpoint - 2)]

/-- Along a two-edge path, one neutral transport step followed by the first
two-level mismatch strictly lowers the total quadratic load on the path. -/
theorem path_two_forward_strict
    (M middle endpoint : ℤ)
    (hmiddle : middle = M ∨ middle = M - 1)
    (hend : endpoint ≤ M - 2) :
    let first := balanceForward M middle
    let second := balanceForward first.2 endpoint
    first.1 ^ 2 + quadraticLoad second <
      M ^ 2 + middle ^ 2 + endpoint ^ 2 := by
  rw [balanceForward_transports_high M middle hmiddle]
  simp only
  have hstrict := path_one_forward_strict M endpoint hend
  simp only [quadraticLoad] at hstrict ⊢
  nlinarith

/-- Along a three-edge path, two neutral transport steps followed by the
first two-level mismatch strictly lower the total quadratic load on the
path. -/
theorem path_three_forward_strict
    (M middle₁ middle₂ endpoint : ℤ)
    (hmiddle₁ : middle₁ = M ∨ middle₁ = M - 1)
    (hmiddle₂ : middle₂ = M ∨ middle₂ = M - 1)
    (hend : endpoint ≤ M - 2) :
    let first := balanceForward M middle₁
    let second := balanceForward first.2 middle₂
    let third := balanceForward second.2 endpoint
    first.1 ^ 2 + second.1 ^ 2 + quadraticLoad third <
      M ^ 2 + middle₁ ^ 2 + middle₂ ^ 2 + endpoint ^ 2 := by
  rw [balanceForward_transports_high M middle₁ hmiddle₁]
  simp only
  rw [balanceForward_transports_high M middle₂ hmiddle₂]
  simp only
  have hstrict := path_one_forward_strict M endpoint hend
  simp only [quadraticLoad] at hstrict ⊢
  nlinarith

end ObserverPatchHolography.DirectedSeamRepair

/- Axiom audit: exact integer arithmetic with standard axioms. -/

#print axioms ObserverPatchHolography.DirectedSeamRepair.quadraticLoad_forward_le
#print axioms ObserverPatchHolography.DirectedSeamRepair.carrier_path_length_at_most_three
#print axioms ObserverPatchHolography.DirectedSeamRepair.shellEnergy_eq_iff_balanced
#print axioms ObserverPatchHolography.DirectedSeamRepair.path_three_forward_strict
