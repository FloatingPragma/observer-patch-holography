import Std.Tactic

/-! Integer height certificate for timed material echoes across opposite base
faces of a fourfold triangle refinement. The native base incidence gives
constant heights zero and two on the two chosen faces. Each nonbacktracking
crossing uses a fresh port with at least one internal transfer between
crossings. These two input checks connect the inequalities to the apparatus. -/

namespace OPH.RefinedTriangleEcho

def Near (a b : Nat) : Prop := a ≤ b + 1 ∧ b ≤ a + 1

def Small (a b c : Nat) : Prop := Near a b ∧ Near b c ∧ Near c a

/-- Old vertices double and each edge midpoint receives the endpoint sum. -/
theorem children_small {a b c : Nat} (h : Small a b c) :
    Small (2*a) (a+b) (c+a) ∧ Small (2*b) (b+c) (a+b) ∧
    Small (2*c) (c+a) (b+c) ∧ Small (a+b) (b+c) (c+a) := by
  simp only [Small, Near] at *
  omega

/-- Two triangles sharing an edge have height sums at most two apart. -/
theorem adjacent_sum {a b c d : Nat} (h₁ : Small a b c) (h₂ : Small a b d) :
    a + b + d ≤ a + b + c + 2 := by
  simp only [Small, Near] at *
  omega

theorem sum_on_path (n : Nat) (f : Nat → Nat)
    (steps : ∀ i, i < n → f (i+1) ≤ f i + 2) :
    f n ≤ f 0 + 2*n := by
  induction n with
  | zero => simp
  | succ k ih =>
    have h := ih (fun i hi => steps i (by omega))
    have hk := steps k (by omega)
    omega

/-- Orient each leg's list of crossed faces from the low to the high face,
reversing the return list. Flight durations are in material-transfer ticks. -/
theorem opposite_echo_bound (L n m outbound inbound : Nat)
    (forward backward : Nat → Nat)
    (fo : forward 0 = 0) (bo : backward 0 = 0)
    (fe : forward n = 6 * 2^L) (be : backward m = 6 * 2^L)
    (fs : ∀ i, i < n → forward (i+1) ≤ forward i + 2)
    (bs : ∀ i, i < m → backward (i+1) ≤ backward i + 2)
    (ft : 2*n ≤ outbound + 1) (bt : 2*m ≤ inbound + 1) :
    12 * 2^L - 2 ≤ outbound + inbound := by
  have h₁ := sum_on_path n forward fs
  have h₂ := sum_on_path m backward bs
  omega

end OPH.RefinedTriangleEcho
