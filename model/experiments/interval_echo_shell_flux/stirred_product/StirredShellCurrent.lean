import Mathlib.Algebra.BigOperators.Field
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

open scoped BigOperators

namespace OPH.StirredShellCurrent

variable {V : Type*} [Fintype V]

/-- A shell is the fiber of the independently ranged integer layer. -/
def shell (layer : V → ℕ) (n : ℕ) : Finset V :=
  Finset.univ.filter fun u => layer u = n

/-- Flight target selected by a fresh uniform accepted word on the same shell.
The input mass counts active arrivals per full renewal cycle, so collecting
addresses have zero active input mass. -/
def targetMass (layer : V → ℕ) (mu : V → ℚ) (x : V) : ℚ :=
  ∑ u, mu u * (if layer u = layer x then 1 / ((shell layer (layer x)).card : ℚ) else 0)

/-- Equal target mass within a shell follows by summing over every possible
old word, even if the arrivals are angularly nonuniform. -/
theorem targetMass_is_average (layer : V → ℕ) (mu : V → ℚ) (x : V) :
    targetMass layer mu x = (∑ u ∈ shell layer (layer x), mu u) /
      ((shell layer (layer x)).card : ℚ) := by
  classical
  unfold targetMass
  simp only [mul_ite, mul_zero, div_eq_mul_inv, one_mul]
  rw [Finset.sum_mul]
  simp only [shell, Finset.sum_filter]

/-- Signed current of one spatially offered ordered bond from independently
chosen target and direction. Neither term uses the total aperture count. -/
def edgeCurrent (layer : V → ℕ) (mu : V → ℚ) (offer : V → V → ℚ) (i j : V) : ℚ :=
  targetMass layer mu i * offer i j - targetMass layer mu j * offer j i

def perWord (layer : V → ℕ) (mu : V → ℚ) (n : ℕ) : ℚ :=
  (∑ u ∈ shell layer n, mu u) / ((shell layer n).card : ℚ)

/-- Every actual outward bond between the same two layers has this same signed
mean. The assumption on the tag concerns the contact probability conditional
on an active word. A collecting word has zero active mass and sends no return. -/
theorem edgeCurrent_from_tags (layer : V → ℕ) (mu : V → ℚ) (offer : V → V → ℚ)
    (n : ℕ) (i j : V) (hi : layer i = n) (hj : layer j = n+1)
    (forward : offer i j = 1/8)
    (reverse : offer j i = 1/8 ∨ targetMass layer mu j = 0) :
    edgeCurrent layer mu offer i j = (perWord layer mu n - perWord layer mu (n+1))/8 := by
  have ret : targetMass layer mu j * offer j i = targetMass layer mu j / 8 := by
    rcases reverse with fair | inactive
    · rw [fair]; ring
    · rw [inactive]; ring
  rw [edgeCurrent, forward, ret, targetMass_is_average, targetMass_is_average, hi, hj]
  unfold perWord
  ring

/-- A signed crossing is evaluated only from the two addresses actually read. -/
def crossing (layer : V → ℕ) (n : ℕ) (i j : V) : ℚ :=
  (if layer i ≤ n then 1 else 0) - (if layer j ≤ n then 1 else 0)

omit [Fintype V] in
theorem tangential_crossing_zero (layer : V → ℕ) (n : ℕ) (i j : V)
    (h : layer i = layer j) : crossing layer n i j = 0 := by
  simp [crossing, h]

/-- Swaps on the address route contribute only their actual boundary crossings;
on every recorded path the signed sum telescopes, including repeated visits. -/
def pathCurrent (layer : V → ℕ) (n : ℕ) (start : V) : List V → ℚ
  | [] => 0
  | next :: tail => crossing layer n start next + pathCurrent layer n next tail

def lastAddress (start : V) : List V → V
  | [] => start
  | next :: tail => lastAddress next tail

omit [Fintype V] in
theorem pathCurrent_endpoints (layer : V → ℕ) (n : ℕ) (start : V) (tail : List V) :
    pathCurrent layer n start tail = crossing layer n start (lastAddress start tail) := by
  induction tail generalizing start with
  | nil => exact (sub_self _).symm
  | cons v tail ih =>
    rw [pathCurrent, lastAddress, ih]
    unfold crossing
    ring

omit [Fintype V] in
/-- A cycle's flight starts inside and ends outside the measured shell. No
fictitious collector-to-source return edge enters the signed count. -/
theorem cycle_net_one (layer : V → ℕ) (n : ℕ) (start : V) (tail : List V)
    (hs : layer start ≤ n) (hf : ¬layer (lastAddress start tail) ≤ n) :
    pathCurrent layer n start tail = 1 := by
  rw [pathCurrent_endpoints]
  simp [crossing, hs, hf]

/-- The independent input is the total signed crossing per cycle (one for
paths ending at the collector). Conditional fair tags and uniform stirring
give equal shares on every actual ordered edge. The measured cut cardinality
and common clock then give the current in units per appointment. -/
theorem cubical_current_per_appointment
    (layer : V → ℕ) (mu : V → ℚ) (offer : V → V → ℚ)
    (edges : Finset (V × V)) (n : ℕ) (hn : 0 < n) (T : ℚ) (hT : 0 < T)
    (geometry : edges.card = 24*n^2)
    (sides : ∀ e ∈ edges, layer e.1 = n ∧ layer e.2 = n+1)
    (tags : ∀ e ∈ edges, offer e.1 e.2 = 1/8 ∧
      (offer e.2 e.1 = 1/8 ∨ targetMass layer mu e.2 = 0))
    (total : ∑ e ∈ edges, edgeCurrent layer mu offer e.1 e.2 = 1)
    (e : V × V) (he : e ∈ edges) :
    edgeCurrent layer mu offer e.1 e.2 / T = (1/T) / (24*(n:ℚ)^2) := by
  let j := (perWord layer mu n - perWord layer mu (n+1))/8
  have eqj : ∀ b ∈ edges, edgeCurrent layer mu offer b.1 b.2 = j := by
    intro b hb
    exact edgeCurrent_from_tags layer mu offer n b.1 b.2
      (sides b hb).1 (sides b hb).2 (tags b hb).1 (tags b hb).2
  have hsum : (edges.card : ℚ) * j = 1 := by
    rw [← total, Finset.sum_congr rfl eqj, Finset.sum_const, nsmul_eq_mul]
  rw [geometry] at hsum
  push_cast at hsum
  rw [eqj e he]
  have hnq : (n:ℚ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
  have htq : T ≠ 0 := ne_of_gt hT
  field_simp [htq, hnq]
  nlinarith only [hsum]

/-- Balance of visits in a complete renewal. Index k denotes measured layer k+1;
the final outer layer may collect and need not have a balance equation.
The two counts of visits come from arrivals and departures respectively. -/
structure LocalVisitBalance (last : ℕ) (visits stays fwd back : ℕ → ℚ) : Prop where
  firstArrival : visits 0 = 1 + stays 0 + back 0
  firstDeparture : visits 0 = stays 0 + fwd 0
  arrival : ∀ k, k < last → visits (k+1) = stays (k+1) + fwd k + back (k+1)
  departure : ∀ k, k < last → visits (k+1) = stays (k+1) + back k + fwd (k+1)

/-- Conservation follows by cancelling the visits and within-shell stays at
successive layers. The outgoing and returning crossing rates may differ at
every boundary. -/
theorem local_balance_transmits_one {last : ℕ} {visits stays fwd back : ℕ → ℚ}
    (h : LocalVisitBalance last visits stays fwd back) :
    ∀ k, k ≤ last → fwd k - back k = 1 := by
  intro k
  induction k with
  | zero =>
    intro _
    linarith only [h.firstArrival, h.firstDeparture]
  | succ k ih =>
    intro hk
    have hlt : k < last := hk
    have prior := ih (Nat.le_of_lt hlt)
    linarith only [prior, h.arrival k hlt, h.departure k hlt]

/-- The combined instrument implication does not assume a common net current
or constant shell total. The two directional aperture sums are the means of
actual target/tag events. Local renewal visit balance supplies their difference;
the independent geometry census fixes the squared denominator. -/
theorem cubical_current_from_local_visits
    (layer : V → ℕ) (mu : V → ℚ) (offer : V → V → ℚ)
    (last k : ℕ) (hk : k ≤ last) (visits stays fwd back : ℕ → ℚ)
    (balance : LocalVisitBalance last visits stays fwd back)
    (edges : Finset (V × V)) (T : ℚ) (hT : 0 < T)
    (geometry : edges.card = 24*(k+1)^2)
    (sides : ∀ e ∈ edges, layer e.1 = k+1 ∧ layer e.2 = k+1+1)
    (tags : ∀ e ∈ edges, offer e.1 e.2 = 1/8 ∧
      (offer e.2 e.1 = 1/8 ∨ targetMass layer mu e.2 = 0))
    (outRead : fwd k = ∑ e ∈ edges, targetMass layer mu e.1 * offer e.1 e.2)
    (inRead : back k = ∑ e ∈ edges, targetMass layer mu e.2 * offer e.2 e.1)
    (e : V × V) (he : e ∈ edges) :
    edgeCurrent layer mu offer e.1 e.2 / T = (1/T) / (24*((k+1:ℕ):ℚ)^2) := by
  have total : ∑ b ∈ edges, edgeCurrent layer mu offer b.1 b.2 = 1 := by
    simp only [edgeCurrent, Finset.sum_sub_distrib, ← outRead, ← inRead]
    exact local_balance_transmits_one balance k hk
  exact cubical_current_per_appointment layer mu offer edges (k+1)
    (Nat.succ_pos k) T hT geometry sides tags total e he

/-- Positive cycle length and a nonempty positive-radius aperture give a
strictly positive prediction. -/
theorem cubical_coefficient_positive (T : ℚ) (hT : 0 < T) (n : ℕ) (hn : 0 < n) :
    0 < (1/T) / (24*(n:ℚ)^2) := by
  have radius_pos : (0 : ℚ) < n := Nat.cast_pos.mpr hn
  exact div_pos (div_pos zero_lt_one hT) (mul_pos (by decide) (pow_pos radius_pos 2))

end OPH.StirredShellCurrent
