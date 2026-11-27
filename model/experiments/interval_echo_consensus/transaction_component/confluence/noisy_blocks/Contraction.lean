import Mathlib.Tactic

namespace OPH.NoisyTriangleRepair

abbrev Address := Fin 6
abbrev Offer := Fin 8
abbrev Flags := Address → Bool

/-- Order 01, 02, 12, 012 on each triangle, interleaved between triangles. -/
def hit (o : Offer) (i : Address) : Bool :=
  decide (i.val / 3 = o.val % 2) &&
    match o.val / 2 with
    | 0 => decide (i.val % 3 ≠ 2)
    | 1 => decide (i.val % 3 ≠ 1)
    | 2 => decide (i.val % 3 ≠ 0)
    | _ => true

def force (o : Offer) (b : Bool) (s : Flags) : Flags := fun i => if hit o i then b else s i
def bit (b : Bool) : ℚ := if b then 1 else 0
def distance (s : Flags) : ℚ := ∑ i, bit (s i)

/-- The three branch weights apply once to the complete offered mask. -/
def expect : List Offer → (Flags → ℚ) → Flags → ℚ
  | [], f, s => f s
  | o :: rest, f, s =>
    12/16 * expect rest f (force o false s) + 3/16 * expect rest f s +
      1/16 * expect rest f (force o true s)

def visits (w : List Offer) (i : Address) := (w.filter (fun o => hit o i)).length
def marginal (k : ℕ) (p : ℚ) : ℚ := (p-1/13) * (3/16)^k + 1/13

theorem coordinate (w : List Offer) (s : Flags) (i : Address) :
    expect w (fun q => bit (q i)) s = marginal (visits w i) (bit (s i)) := by
  induction w generalizing s with
  | nil => simp [expect, visits, marginal]
  | cons o w ih =>
    simp only [expect, ih]
    cases hh : hit o i <;> cases hq : s i <;>
      simp [visits, hh, force, hq, marginal, bit, pow_succ] <;> ring

theorem sum_expect (w : List Offer) (s : Flags) (f : Address → Flags → ℚ) :
    expect w (fun q => ∑ i, f i q) s = ∑ i, expect w (f i) s := by
  induction w generalizing s with
  | nil => rfl
  | cons o w ih =>
    simp only [expect, ih, Finset.sum_add_distrib, Finset.mul_sum]

theorem expected_distance (w : List Offer) (s : Flags) :
    expect w distance s = ∑ i, marginal (visits w i) (bit (s i)) := by
  unfold distance
  rw [sum_expect w s (fun i q => bit (q i))]
  simp only [coordinate]

def Fair (w : List Offer) := w.Perm (List.finRange 8)

theorem three_exposures {w : List Offer} (hw : Fair w) (i : Address) : visits w i = 3 := by
  have ref : ∀ a : Address, visits (List.finRange 8) a = 3 := by decide
  exact ((hw.filter (fun o => hit o i)).length_eq).trans (ref i)

theorem block_mean {w : List Offer} (hw : Fair w) (s : Flags) :
    expect w distance s = (27/4096)*distance s + 1878/4096 := by
  rw [expected_distance, distance, Finset.mul_sum]
  simp only [three_exposures hw]
  have term (i : Address) : marginal 3 (bit (s i)) = (27/4096)*bit (s i)+313/4096 := by
    unfold marginal; ring
  simp only [term, Finset.sum_add_distrib]
  norm_num

theorem split_visits (u v : List Offer) (i : Address) :
    visits (u++v) i = visits u i + visits v i := by simp [visits]

def offset (last : Offer) (i : Address) : ℚ := if hit last i then 304/4096 else 313/4096

set_option maxHeartbeats 1000000 in
theorem offset_budget (p : Offer) : (∑ i : Address, offset p i) ≤ 1860/4096 := by
  fin_cases p <;> norm_num [Fin.sum_univ_succ, offset, hit]

set_option maxHeartbeats 1000000 in
theorem prefix_bound (w tail : List Offer) (o : Offer)
    (h : Fair (w ++ o :: tail)) (s : Flags) :
    expect w distance s ≤ distance s + 1860/4096 := by
  have bound (i : Address) :
      marginal (visits w i) (bit (s i)) ≤ bit (s i) + offset o i := by
    have total := three_exposures h i
    rw [split_visits] at total
    generalize hn : visits w i = k at total ⊢
    cases hh : hit o i with
    | false =>
      have hk : k ≤ 3 := by omega
      interval_cases k <;> cases s i <;> norm_num [marginal, bit, offset, hh]
    | true =>
      have pos : 0 < visits (o::tail) i := by simp [visits, hh]
      have hk : k ≤ 2 := by omega
      interval_cases k <;> cases s i <;> norm_num [marginal, bit, offset, hh]
  rw [expected_distance]
  calc
    _ ≤ ∑ i, (bit (s i) + offset o i) := Finset.sum_le_sum (fun i _ => bound i)
    _ = distance s + ∑ i, offset o i := by rw [Finset.sum_add_distrib, distance]
    _ ≤ _ := by linarith [offset_budget o]

/-- Fresh block-start choices all have the same conditional first moment.
The sequence recursion here is the resulting law of total expectation. -/
theorem tube (x : ℕ → ℚ)
    (step : ∀ m, x (m+1) = (27/4096)*x m+1878/4096) (m : ℕ) :
    x m = 6/13 + (27/4096)^m * (x 0 - 6/13) := by
  induction m with
  | zero => ring
  | succ m ih => rw [step, ih, pow_succ]; ring

end OPH.NoisyTriangleRepair
