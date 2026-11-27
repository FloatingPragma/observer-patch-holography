import Mathlib
import A5PortAction

namespace OPH.StaticFamily

open OPH.A5PortAction (app natAdj natNeighbors perms)

/-- Both orientations of a marked edge with its two face tips. -/
def markedEdgeTiles : List (List Nat) :=
  (List.range 12).flatMap fun a => (natNeighbors a).map fun b =>
    a :: b :: (natNeighbors a).filter (natAdj b)

/-- Free order-three rotations from the exact port action. -/
def orderThree : List (List Nat) :=
  perms.filter fun p => (List.range 12).all fun i =>
    app p i != i && app p (app p (app p i)) == i

/-- Input rails of disjoint rotated copies of a marked-edge tile. -/
def disjointInputRails : List (List Nat) :=
  markedEdgeTiles.flatMap fun tile => orderThree.filterMap fun p =>
    let first := tile.map (app p)
    let second := first.map (app p)
    if (tile ++ first ++ second).eraseDups.length = 12 then
      some [tile.headD 0, first.headD 0, second.headD 0]
    else none

set_option maxRecDepth 10000 in
set_option maxHeartbeats 10000000 in
/-- No disjoint analyzer tiling puts two input rails on a common native seam. -/
theorem analyzer_input_rails_independent :
    disjointInputRails.length = 240 ∧
      disjointInputRails.all (fun rail => rail.all fun i => rail.all fun j => !natAdj i j) = true := by
  decide

/-- A projection and its complement have opposite reflection signs. -/
theorem paired_couplings_preserve_ports {A : Type*} [Ring A] (p q : A)
    (hp : p * p = p) :
    p * (q + (p + p - 1) * q * (p + p - 1)) * (1 - p) = 0 := by
  have left : p * (p + p - 1) = p := by
    rw [mul_sub, mul_add, hp, mul_one]
    abel
  have right : (p + p - 1) * (1 - p) = -(1 - p) := by
    have hz : p * (1 - p) = 0 := by rw [mul_sub, mul_one, hp, sub_self]
    rw [sub_mul, add_mul, hz, one_mul]
    abel
  rw [mul_add, add_mul]
  simp only [← mul_assoc, left]
  rw [mul_assoc (p * q), right, mul_neg, add_neg_cancel]

/-- The single-projector smaller roots increase every ratio to the light root. -/
theorem single_projector_ratio {l m s q : ℝ} (hl : 0 < l) (hm : l < m) (hs : m < s)
    (roots : m * (s - m) = q * l * (s - l)) : q < m / l := by
  apply (lt_div_iff₀ hl).2
  have hsl : 0 < s - l := by linarith
  have hml : 0 < m * (m - l) := mul_pos (by linarith) (sub_pos.mpr hm)
  exact lt_of_mul_lt_mul_right (by nlinarith only [hml, roots]) (le_of_lt hsl)

/-- Conservative margins from the two ratio intervals exclude a common virtual gap. -/
theorem single_virtual_gap_excluded {eta u v qu qv : ℝ} (heta : 0 ≤ eta)
    (hu0 : 206 ≤ u) (hu : u ≤ 206.768292) (hqu : 206.770315 ≤ qu)
    (hv : 3477.012 ≤ v) (hqv : qv ≤ 3477.473) :
    ¬(u * (1 + eta * u) = qu * (1 + eta) ∧
      v * (1 + eta * v) = qv * (1 + eta)) := by
  rintro ⟨mu, tau⟩
  have u2 : u ^ 2 ≤ 42849 := by
    have h := mul_nonneg (show 0 ≤ 207 - u by linarith) (show 0 ≤ 207 + u by linarith)
    nlinarith only [h]
  have v2 : 12089529 ≤ v ^ 2 := by
    have h := mul_nonneg (show 0 ≤ v - 3477 by linarith) (show 0 ≤ v + 3477 by linarith)
    nlinarith only [h]
  have eta_lo : 0.000000045 < eta := by
    by_contra! upper
    have hprod := mul_nonneg heta (show 0 ≤ 42849 - (u ^ 2 - qu) by linarith)
    nlinarith only [hprod, upper, mu, hu, hqu]
  have hprod := mul_nonneg heta (show 0 ≤ v ^ 2 - qv - 12085000 by linarith)
  nlinarith only [hprod, eta_lo, tau, hv, hqv]

/-- The squared Wilson regulator leaves a strict one-valley derivative. -/
theorem squared_wilson_derivative {u s am : ℝ} (hu : 0 ≤ u) (hus : u ≤ s) (ham : 0 ≤ am) :
    0 < 1 - u + 2 * s * (am + s ^ 2) := by
  have hs : 0 ≤ s := le_trans hu hus
  have hc : 0 ≤ 2 * s * (am + s ^ 2) := by positivity
  by_cases h : u < 1
  · linarith
  · have h1 : 1 ≤ u := le_of_not_gt h
    have h2 : 1 ≤ s ^ 2 := by nlinarith [h1, hus, sq_nonneg (s - 1)]
    have h3 := mul_nonneg hs (show 0 ≤ am + s ^ 2 - 1 by linarith)
    nlinarith only [h3, hus, hs]

/-- A circulant two-current kernel preserves the individual family tags exactly
when both nonconstant Fourier coefficients vanish. -/
theorem coherent_current_iff {K : Type*} [Field K] {w a b c : K}
    (root : w ^ 2 + w + 1 = 0) (distinct : w ≠ w ^ 2) :
    (a + b * w + c * w ^ 2 = 0 ∧ a + b * w ^ 2 + c * w = 0) ↔
      (a = b ∧ b = c) := by
  constructor
  · rintro ⟨first, second⟩
    have product : (b - c) * (w - w ^ 2) = 0 := by
      linear_combination first - second
    have bc := sub_eq_zero.mp ((mul_eq_zero.mp product).resolve_right (sub_ne_zero.mpr distinct))
    refine ⟨?_, bc⟩
    rw [← bc] at first
    linear_combination first - b * root
  · rintro ⟨ab, bc⟩
    subst b
    subst c
    constructor <;> linear_combination a * root

end OPH.StaticFamily
