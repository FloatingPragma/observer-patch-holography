import OphGap.SignedGap

/-!
# A certificate checker for the hypotheses of `gram_posDef_of_negative_cycle`, with soundness

Data format. A graph on vertices `0, …, n-1` is a list `es : List Edge` of triples
`(src, tgt, neg)`; `neg = true` means sign `-1`, `neg = false` means sign `+1`.
A certificate consists of the *ordering* of `es` and a cycle description:

* the first `n-1` edges form a spanning tree rooted at vertex `0`, in the sense that
  edge number `i` (0-based, `i < n-1`) has `tgt = i+1` and `src < i+1`;
* `cyc : List (Nat × Bool)` lists edge indices with a direction flag (`true` = traverse
  `src → tgt`), describing a closed walk starting and ending at vertex `0` whose
  sign product is `-1`.

`check n es cyc = true` is decidable by the kernel. The theorem `posDef_of_check`
shows it implies the Gram matrix `Dᵀ D` of the signed derivative built from `es`
is positive definite. The checks use kernel reduction only.
-/

open Matrix

namespace GapScout

abbrev Edge := Nat × Nat × Bool

/-- Every edge has both endpoints `< n` and is not a loop. -/
def edgesOk (n : Nat) : List Edge → Bool
  | [] => true
  | (a, b, _) :: es => (decide (a < n) && decide (b < n) && decide (a ≠ b)) && edgesOk n es

/-- `treeOk k i es`: the next `k` edges, starting at index `i`, are tree edges:
edge at index `j` has `tgt = j+1` and `src < j+1`. -/
def treeOk : Nat → Nat → List Edge → Bool
  | 0, _, _ => true
  | _ + 1, _, [] => false
  | k + 1, i, (a, b, _) :: es => (decide (b = i + 1) && decide (a < i + 1)) && treeOk k (i + 1) es

/-- Walk the cycle: returns the final vertex and the sign parity (`true` = negative). -/
def walk (es : List Edge) : Nat → Bool → List (Nat × Bool) → Option (Nat × Bool)
  | cur, par, [] => some (cur, par)
  | cur, par, (i, fwd) :: rest =>
    match es[i]? with
    | none => none
    | some (a, b, neg) =>
      if fwd then (if a = cur then walk es b (par ^^ neg) rest else none)
      else (if b = cur then walk es a (par ^^ neg) rest else none)

def check (n : Nat) (es : List Edge) (cyc : List (Nat × Bool)) : Bool :=
  decide (0 < n) && edgesOk n es && treeOk (n - 1) 0 es &&
    (match walk es 0 false cyc with
      | some (c, p) => decide (c = 0) && p
      | none => false)

/-! ## From the data to `src`, `tgt`, `σ` -/

theorem edgesOk_mem {n : Nat} {es : List Edge} (h : edgesOk n es = true) :
    ∀ x ∈ es, x.1 < n ∧ x.2.1 < n ∧ x.1 ≠ x.2.1 := by
  induction es with
  | nil => simp
  | cons x es ih =>
    obtain ⟨a, b, s⟩ := x
    simp only [edgesOk, Bool.and_eq_true, decide_eq_true_eq] at h
    intro y hy
    simp only [List.mem_cons] at hy
    rcases hy with rfl | hy
    · exact ⟨h.1.1.1, h.1.1.2, h.1.2⟩
    · exact ih h.2 y hy

def srcOf (n : Nat) (es : List Edge) (hok : edgesOk n es = true) (e : Fin es.length) : Fin n :=
  ⟨(es.get e).1, (edgesOk_mem hok _ (List.get_mem es e)).1⟩

def tgtOf (n : Nat) (es : List Edge) (hok : edgesOk n es = true) (e : Fin es.length) : Fin n :=
  ⟨(es.get e).2.1, (edgesOk_mem hok _ (List.get_mem es e)).2.1⟩

def sgn (b : Bool) : ℝ := if b then -1 else 1

def σOf (es : List Edge) (e : Fin es.length) : ℝ := sgn (es.get e).2.2

theorem sgn_xor (p q : Bool) : sgn p * sgn q = sgn (p ^^ q) := by
  cases p <;> cases q <;> simp [sgn]

theorem loop_free (n : Nat) (es : List Edge) (hok : edgesOk n es = true) :
    ∀ e, srcOf n es hok e ≠ tgtOf n es hok e := by
  intro e h
  exact (edgesOk_mem hok _ (List.get_mem es e)).2.2 (congrArg Fin.val h)

theorem signs_pm (es : List Edge) : ∀ e, σOf es e = 1 ∨ σOf es e = -1 := by
  intro e
  unfold σOf sgn
  split <;> simp

/-! ## Spanning tree ⇒ every vertex reachable from `0` -/

theorem treeOk_get {k i : Nat} {es : List Edge} (h : treeOk k i es = true) :
    ∀ j, j < k → ∃ hj : j < es.length,
      (es[j]'hj).2.1 = i + j + 1 ∧ (es[j]'hj).1 < i + j + 1 := by
  induction k generalizing i es with
  | zero => intro j hj; omega
  | succ k ih =>
    intro j hj
    cases es with
    | nil => simp [treeOk] at h
    | cons x es =>
      obtain ⟨a, b, s⟩ := x
      simp only [treeOk, Bool.and_eq_true, decide_eq_true_eq] at h
      cases j with
      | zero =>
        refine ⟨by simp, ?_, ?_⟩ <;> simp [h.1.1, h.1.2]
      | succ j =>
        obtain ⟨hj', h1, h2⟩ := ih h.2 j (by omega)
        refine ⟨by simp; omega, ?_, ?_⟩
        · simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using h1
        · simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using h2

theorem reach_all {n : Nat} {es : List Edge} (hok : edgesOk n es = true)
    (htree : treeOk (n - 1) 0 es = true) (h0 : 0 < n) :
    ∀ v : Fin n, ∃ s, SignedReach (srcOf n es hok) (tgtOf n es hok) (σOf es) ⟨0, h0⟩ v s := by
  suffices H : ∀ k (hk : k < n), ∃ s,
      SignedReach (srcOf n es hok) (tgtOf n es hok) (σOf es) ⟨0, h0⟩ ⟨k, hk⟩ s by
    intro v; exact H v.1 v.2
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
    intro hk
    cases k with
    | zero => exact ⟨1, SignedReach.refl _⟩
    | succ k =>
      obtain ⟨hj, htgt, hsrc⟩ := treeOk_get htree k (by omega)
      obtain ⟨s, hs⟩ := ih (es[k]'hj).1 (by omega) (by omega)
      refine ⟨s * σOf es ⟨k, hj⟩, ?_⟩
      have step := SignedReach.fwd (src := srcOf n es hok) (tgt := tgtOf n es hok)
        (σ := σOf es) ⟨k, hj⟩ hs
      have e : tgtOf n es hok ⟨k, hj⟩ = ⟨k + 1, hk⟩ := Fin.ext (by simp [tgtOf]; omega)
      rw [e] at step
      exact step

/-! ## The walk ⇒ a negative closed walk at `0` -/

theorem walk_sound {n : Nat} {es : List Edge} (hok : edgesOk n es = true) (h0 : 0 < n) :
    ∀ (rest : List (Nat × Bool)) (cur : Nat) (par : Bool) (c : Nat) (p : Bool),
      walk es cur par rest = some (c, p) →
      ∀ (hcur : cur < n),
        SignedReach (srcOf n es hok) (tgtOf n es hok) (σOf es) ⟨0, h0⟩ ⟨cur, hcur⟩ (sgn par) →
        ∃ hc : c < n,
          SignedReach (srcOf n es hok) (tgtOf n es hok) (σOf es) ⟨0, h0⟩ ⟨c, hc⟩ (sgn p) := by
  intro rest
  induction rest with
  | nil =>
    intro cur par c p hw hcur hr
    simp only [walk, Option.some.injEq, Prod.mk.injEq] at hw
    obtain ⟨rfl, rfl⟩ := hw
    exact ⟨hcur, hr⟩
  | cons x rest ih =>
    obtain ⟨i, fwd⟩ := x
    intro cur par c p hw hcur hr
    simp only [walk] at hw
    split at hw
    · exact absurd hw (by simp)
    · rename_i a b neg hget
      obtain ⟨hi, hget'⟩ := List.getElem?_eq_some_iff.mp hget
      have hmem := edgesOk_mem hok _ (List.getElem_mem hi)
      rw [hget'] at hmem
      split at hw
      · -- forward: a = cur, move to b
        split at hw
        · rename_i ha
          subst ha
          have hr' : SignedReach (srcOf n es hok) (tgtOf n es hok) (σOf es) ⟨0, h0⟩
              (srcOf n es hok ⟨i, hi⟩) (sgn par) := by
            have e : srcOf n es hok ⟨i, hi⟩ = ⟨a, hcur⟩ := Fin.ext (by simp [srcOf, hget'])
            rw [e]; exact hr
          have step := SignedReach.fwd (src := srcOf n es hok) (tgt := tgtOf n es hok)
            (σ := σOf es) ⟨i, hi⟩ hr'
          have e2 : tgtOf n es hok ⟨i, hi⟩ = ⟨b, hmem.2.1⟩ := Fin.ext (by simp [tgtOf, hget'])
          have e3 : σOf es ⟨i, hi⟩ = sgn neg := by simp [σOf, hget']
          rw [e2, e3, sgn_xor] at step
          exact ih b (par ^^ neg) c p hw hmem.2.1 step
        · exact absurd hw (by simp)
      · -- backward: b = cur, move to a
        split at hw
        · rename_i hb
          subst hb
          have hr' : SignedReach (srcOf n es hok) (tgtOf n es hok) (σOf es) ⟨0, h0⟩
              (tgtOf n es hok ⟨i, hi⟩) (sgn par) := by
            have e : tgtOf n es hok ⟨i, hi⟩ = ⟨b, hcur⟩ := Fin.ext (by simp [tgtOf, hget'])
            rw [e]; exact hr
          have step := SignedReach.bwd (src := srcOf n es hok) (tgt := tgtOf n es hok)
            (σ := σOf es) ⟨i, hi⟩ hr'
          have e2 : srcOf n es hok ⟨i, hi⟩ = ⟨a, hmem.1⟩ := Fin.ext (by simp [srcOf, hget'])
          have e3 : σOf es ⟨i, hi⟩ = sgn neg := by simp [σOf, hget']
          rw [e2, e3, sgn_xor] at step
          exact ih a (par ^^ neg) c p hw hmem.1 step
        · exact absurd hw (by simp)

/-! ## Unpacking `check` -/

theorem check_pos {n es cyc} (h : check n es cyc = true) : 0 < n := by
  unfold check at h
  simp only [Bool.and_eq_true, decide_eq_true_eq] at h
  exact h.1.1.1

theorem check_edgesOk {n es cyc} (h : check n es cyc = true) : edgesOk n es = true := by
  unfold check at h
  simp only [Bool.and_eq_true] at h
  exact h.1.1.2

theorem check_treeOk {n es cyc} (h : check n es cyc = true) : treeOk (n - 1) 0 es = true := by
  unfold check at h
  simp only [Bool.and_eq_true] at h
  exact h.1.2

theorem check_walk {n es cyc} (h : check n es cyc = true) :
    walk es 0 false cyc = some (0, true) := by
  unfold check at h
  simp only [Bool.and_eq_true] at h
  have h2 := h.2
  split at h2
  · rename_i c p hw
    simp only [Bool.and_eq_true, decide_eq_true_eq] at h2
    obtain ⟨rfl, rfl⟩ := h2
    exact hw
  · simp at h2

theorem negative_cycle_of_check {n es cyc} (h : check n es cyc = true) :
    SignedReach (srcOf n es (check_edgesOk h)) (tgtOf n es (check_edgesOk h)) (σOf es)
      ⟨0, check_pos h⟩ ⟨0, check_pos h⟩ (-1) := by
  have hrefl : SignedReach (srcOf n es (check_edgesOk h)) (tgtOf n es (check_edgesOk h)) (σOf es)
      ⟨0, check_pos h⟩ ⟨0, check_pos h⟩ (sgn false) := by
    simpa [sgn] using SignedReach.refl (src := srcOf n es (check_edgesOk h))
      (tgt := tgtOf n es (check_edgesOk h)) (σ := σOf es) (⟨0, check_pos h⟩ : Fin n)
  obtain ⟨hc, hr⟩ := walk_sound (check_edgesOk h) (check_pos h) cyc 0 false 0 true
    (check_walk h) (check_pos h) hrefl
  simpa [sgn] using hr

theorem connected_of_check {n es cyc} (h : check n es cyc = true) :
    ∀ v : Fin n, ∃ s, SignedReach (srcOf n es (check_edgesOk h)) (tgtOf n es (check_edgesOk h))
      (σOf es) ⟨0, check_pos h⟩ v s :=
  reach_all (check_edgesOk h) (check_treeOk h) (check_pos h)

/-- SOUNDNESS: a passing certificate yields a positive definite Gram matrix. -/
theorem posDef_of_check {n : Nat} {es : List Edge} {cyc : List (Nat × Bool)}
    (h : check n es cyc = true) :
    ((signedDerivative (srcOf n es (check_edgesOk h)) (tgtOf n es (check_edgesOk h)) (σOf es))ᵀ *
      signedDerivative (srcOf n es (check_edgesOk h)) (tgtOf n es (check_edgesOk h)) (σOf es)).PosDef :=
  gram_posDef_of_negative_cycle _ _ _ (loop_free n es _) (signs_pm es) _
    (negative_cycle_of_check h) (connected_of_check h)

/-- Kernel triviality, stated directly. -/
theorem kernel_trivial_of_check {n : Nat} {es : List Edge} {cyc : List (Nat × Bool)}
    (h : check n es cyc = true) (x : Fin n → ℝ)
    (hx : signedDerivative (srcOf n es (check_edgesOk h)) (tgtOf n es (check_edgesOk h)) (σOf es) *ᵥ x = 0) :
    x = 0 :=
  kernel_trivial_of_negative_cycle _ _ _ (loop_free n es _) (signs_pm es) _
    (negative_cycle_of_check h) (connected_of_check h) x hx

/-! ## Opaque statement wrappers

Elaborating the statements above against a concrete 11,816-edge list makes the elaborator
try to unfold the list. These `Prop`-valued wrappers keep the concrete statements
syntactic; unfold them (`unfold GramPosDef`) to see exactly what is asserted. -/

/-- `Dᵀ D` for the graph `es` is positive definite. -/
def GramPosDef (n : Nat) (es : List Edge) (cyc : List (Nat × Bool)) (h : check n es cyc = true) : Prop :=
  ((signedDerivative (srcOf n es (check_edgesOk h)) (tgtOf n es (check_edgesOk h)) (σOf es))ᵀ *
    signedDerivative (srcOf n es (check_edgesOk h)) (tgtOf n es (check_edgesOk h)) (σOf es)).PosDef

theorem gramPosDef_of_check {n es cyc} (h : check n es cyc = true) : GramPosDef n es cyc h :=
  posDef_of_check h

/-- Every eigenvalue of that `Dᵀ D` is strictly positive (so the least one, the gap, is). -/
def EigenvaluesPos (n : Nat) (es : List Edge) (cyc : List (Nat × Bool)) (h : check n es cyc = true) : Prop :=
  ∀ i : Fin n, 0 < (posDef_of_check h).1.eigenvalues i

theorem eigenvaluesPos_of_check {n es cyc} (h : check n es cyc = true) : EigenvaluesPos n es cyc h :=
  fun i => (posDef_of_check h).eigenvalues_pos i

/-- The signed derivative has trivial kernel. -/
def KernelTrivial (n : Nat) (es : List Edge) (cyc : List (Nat × Bool)) (h : check n es cyc = true) : Prop :=
  ∀ x : Fin n → ℝ,
    signedDerivative (srcOf n es (check_edgesOk h)) (tgtOf n es (check_edgesOk h)) (σOf es) *ᵥ x = 0 →
      x = 0

theorem kernelTrivial_of_check' {n es cyc} (h : check n es cyc = true) : KernelTrivial n es cyc h :=
  fun x hx => kernel_trivial_of_check h x hx

/-- Negative-control helper: flip the sign of the edge at index `i`. -/
def flipAt : Nat → List Edge → List Edge
  | _, [] => []
  | 0, (a, b, s) :: es => (a, b, !s) :: es
  | i + 1, e :: es => e :: flipAt i es

end GapScout
