import Mathlib.Data.Complex.Basic
import Mathlib.Data.Real.Sqrt
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Tactic

/-! The physical basis of the reachable output consists of port loads and the
sign/coin label of each addressed collar. The indexing by the input used in
the finite sums below is justified by a left inverse on actual output records. -/

namespace StrictQuantum

open scoped BigOperators
noncomputable section

inductive Load | empty | single | double deriving DecidableEq, Fintype
open Load
abbrev Sign := Bool

/-- true records the oriented contrast +2; false records −2. -/
def sign : Load → Load → Option Sign
  | double, empty => some true
  | empty, double => some false
  | _, _ => none

def leftBefore (b : Sign) : Load := if b then double else empty
def rightBefore (b : Sign) : Load := if b then empty else double

theorem sign_determines {a b : Load} {s : Sign} (h : sign a b = some s) :
    a = leftBefore s ∧ b = rightBefore s := by
  cases a <;> cases b <;> cases s <;> simp_all [sign, leftBefore, rightBefore]

variable {V E : Type} [DecidableEq V]

structure Screen (E V : Type) where
  left : E → V
  right : E → V
  distinct : ∀ e, left e ≠ right e

abbrev Loads (V : Type) := V → Load

def putPair (q : Screen E V) (e : E) (x : Loads V) (a b : Load) : Loads V :=
  fun v => if v = q.left e then a else if v = q.right e then b else x v

def tag (q : Screen E V) (e : E) (x : Loads V) :=
  sign (x (q.left e)) (x (q.right e))

def fire (q : Screen E V) (e : E) (x : Loads V) : Loads V :=
  if (tag q e x).isSome then putPair q e x single single else x

def undo (q : Screen E V) (e : E) (r : Option Sign) (y : Loads V) : Loads V :=
  match r with
  | none => y
  | some s => putPair q e y (leftBefore s) (rightBefore s)

theorem undo_fire (q : Screen E V) (e : E) (x : Loads V) :
    undo q e (tag q e x) (fire q e x) = x := by
  cases h : tag q e x with
  | none => simp [undo, fire, h]
  | some s =>
    have hh := sign_determines (by exact h : sign (x (q.left e)) (x (q.right e)) = some s)
    funext v
    by_cases hl : v = q.left e
    · subst v; simp [undo, putPair, hh.1]
    · by_cases hr : v = q.right e
      · subst v; simp [undo, putPair, hl, hh.2]
      · simp [undo, fire, putPair, h, hl, hr]

/-- For a fixed visit list the record list is indexed by position in that list.
On a duplicate-free list each position is a distinct physical collar. -/
abbrev SharpOutput (V : Type) := Loads V × List (Option Sign)

def scan (q : Screen E V) : List E → Loads V → SharpOutput V
  | [], x => (x, [])
  | e :: es, x =>
    let tail := scan q es (fire q e x)
    (tail.1, tag q e x :: tail.2)

def recover (q : Screen E V) : List E → List (Option Sign) → Loads V → Loads V
  | e :: es, r :: rs, y => undo q e r (recover q es rs y)
  | _, _, y => y

theorem recover_scan (q : Screen E V) (es : List E) (x : Loads V) :
    recover q es (scan q es x).2 (scan q es x).1 = x := by
  induction es generalizing x with
  | nil => rfl
  | cons e es ih => simp only [scan, recover, ih, undo_fire]

/-- Generate a specified classical output branch from the actual strict action.
Each quiet or hit decision is checked *at its service*, on the current loads. -/
def selected (q : Screen E V) (M : E → Bool) : List E → Loads V → Option (SharpOutput V)
  | [], x => some (x, [])
  | e :: es, x =>
    if (tag q e x).isSome = M e then
      (selected q M es (fire q e x)).map fun y => (y.1, tag q e x :: y.2)
    else none

def flags (r : List (Option Sign)) := r.map Option.isSome

theorem selected_scan (q : Screen E V) (M : E → Bool) (es : List E) (x : Loads V) :
    selected q M es x =
      if flags (scan q es x).2 = es.map M then some (scan q es x) else none := by
  induction es generalizing x with
  | nil => rfl
  | cons e es ih =>
    by_cases h : (tag q e x).isSome = M e
    · simp only [selected, h, ↓reduceIte, ih, scan, flags, List.map_cons,
        List.cons.injEq, true_and]
      split <;> simp_all
    · simp [selected, scan, flags, h]

/-- An independent coin qubit is enrolled in each occupied collar only. -/
def Coins : List (Option Sign) → Type
  | [] => Unit
  | none :: rs => Coins rs
  | some _ :: rs => Bool × Coins rs

instance coinFintype : (rs : List (Option Sign)) → Fintype (Coins rs)
  | [] => inferInstanceAs (Fintype Unit)
  | none :: rs => coinFintype rs
  | some _ :: rs => @instFintypeProd _ _ inferInstance (coinFintype rs)

/-- The exact two equal positive amplitudes of the native write. -/
def alpha : ℝ := (Real.sqrt 2)⁻¹

theorem alpha_pos : 0 < alpha := by
  unfold alpha
  positivity

theorem two_alpha_sq : 2 * (alpha * alpha) = 1 := by
  have h : (Real.sqrt 2) ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hpos : 0 < Real.sqrt 2 := by positivity
  unfold alpha
  field_simp
  nlinarith

def weight : List (Option Sign) → ℝ
  | [] => 1
  | none :: rs => weight rs
  | some _ :: rs => alpha * weight rs

theorem weight_pos (rs : List (Option Sign)) : 0 < weight rs := by
  induction rs with
  | nil => simp [weight]
  | cons r rs ih => cases r <;> simp [weight, ih, mul_pos alpha_pos ih]

def labels : (rs : List (Option Sign)) → Coins rs → List (Option (Sign × Bool))
  | [], _ => []
  | none :: rs, c => none :: labels rs c
  | some s :: rs, c => some (s, c.1) :: labels rs c.2

theorem labels_signs (rs : List (Option Sign)) (c : Coins rs) :
    (labels rs c).map (Option.map Prod.fst) = rs := by
  induction rs with
  | nil => rfl
  | cons r rs ih => cases r <;> simp [labels, ih]

theorem labels_injective (rs : List (Option Sign)) : Function.Injective (labels rs) := by
  induction rs with
  | nil => intro a b _; exact Subsingleton.elim (α := Unit) a b
  | cons r rs ih =>
    cases r with
    | none => intro a b h; exact ih (List.cons.inj h).2
    | some s =>
      intro a b h
      exact Prod.ext (by simpa [labels] using (List.cons.inj h).1) (ih (List.cons.inj h).2)

theorem normalized (rs : List (Option Sign)) :
    (∑ _ : Coins rs, weight rs * weight rs) = 1 := by
  induction rs with
  | nil => simp [Coins, weight]
  | cons r rs ih =>
    cases r with
    | none => exact ih
    | some s =>
      change (∑ _ : Bool × Coins rs, (alpha * weight rs) * (alpha * weight rs)) = 1
      rw [Fintype.sum_prod_type]
      simp only [Fintype.sum_bool]
      have he : (∑ _ : Coins rs, (alpha * weight rs) * (alpha * weight rs)) = alpha * alpha := by
        calc
          _ = ∑ _ : Coins rs, (alpha * alpha) * (weight rs * weight rs) :=
            Finset.sum_congr rfl (fun _ _ => by ring)
          _ = _ := by rw [← Finset.mul_sum, ih, mul_one]
      rw [he]
      linarith [two_alpha_sq]

/-- Distinct source and coin labels are orthogonal physical output basis states:
the decoder reads only final loads and actual sign/coin records. -/
def physical (q : Screen E V) (es : List E)
    (z : (x : Loads V) × Coins (scan q es x).2) :=
  ((scan q es z.1).1, labels (scan q es z.1).2 z.2)

theorem physical_injective (q : Screen E V) (es : List E) :
    Function.Injective (physical q es) := by
  rintro ⟨x,c⟩ ⟨y,d⟩ hh
  have hl := congrArg Prod.fst hh
  have hr := congrArg Prod.snd hh
  change (scan q es x).1 = (scan q es y).1 at hl
  change labels (scan q es x).2 c = labels (scan q es y).2 d at hr
  have hs := congrArg (List.map (Option.map Prod.fst)) hr
  simp only [labels_signs] at hs
  have he : x = y := by
    calc
      x = recover q es (scan q es x).2 (scan q es x).1 := (recover_scan ..).symm
      _ = recover q es (scan q es y).2 (scan q es y).1 := by rw [hs, hl]
      _ = y := recover_scan ..
  subst y
  have hc := labels_injective _ hr
  subst d
  rfl

abbrev StateVector (V : Type) := Loads V → ℂ

/-- V in compressed coordinates, using precisely the reachable physical image.
Coordinates outside that image have zero amplitude. -/
def Vout (q : Screen E V) (es : List E) (ψ : StateVector V)
    (x : Loads V) (_ : Coins (scan q es x).2) : ℂ :=
  ψ x * (weight (scan q es x).2 : ℝ)

/-- Sparse multiplication of the addressed Kraus columns. A positive column
contributes α for the chosen coin basis value and a quiet column contributes
one. A rejected branch contributes zero at the first disagreement. -/
def branchCoefficient (q : Screen E V) (M : E → Bool) : List E → Loads V → ℝ
  | [], _ => 1
  | e :: es, x =>
    let r := tag q e x
    if r.isSome = M e then
      (if r.isSome then alpha else 1) * branchCoefficient q M es (fire q e x)
    else 0

theorem coefficient_scan (q : Screen E V) (M : E → Bool) (es : List E) (x : Loads V) :
    branchCoefficient q M es x =
      if flags (scan q es x).2 = es.map M then weight (scan q es x).2 else 0 := by
  induction es generalizing x with
  | nil => rfl
  | cons e es ih =>
    cases hr : tag q e x <;> cases hm : M e <;>
      simp [branchCoefficient, scan, flags, hr, hm, ih, weight]

def project (q : Screen E V) (es : List E) (M : E → Bool) (ψ : StateVector V)
    (x : Loads V) : ℂ :=
  if flags (scan q es x).2 = es.map M then ψ x else 0

def kraus (q : Screen E V) (es : List E) (M : E → Bool) (ψ : StateVector V)
    (x : Loads V) (_ : Coins (scan q es x).2) : ℂ :=
  ψ x * (branchCoefficient q M es x : ℝ)

/-- Derived by running the sequentially selected hit/quiet branch, using the
blank native column normalization at every accepted address. -/
theorem kraus_eq_V_project (q : Screen E V) (es : List E) (M : E → Bool) (ψ : StateVector V) :
    kraus q es M ψ = Vout q es (project q es M ψ) := by
  funext x c
  rw [kraus, coefficient_scan]
  by_cases h : flags (scan q es x).2 = es.map M <;> simp [h, project, Vout]

theorem Vout_norm (q : Screen E V) (es : List E) (ψ : StateVector V) (x : Loads V) :
    (∑ c : Coins (scan q es x).2, Complex.normSq (Vout q es ψ x c)) = Complex.normSq (ψ x) := by
  simp only [Vout, Complex.normSq_mul, Complex.normSq_ofReal]
  rw [← Finset.mul_sum, normalized, mul_one]

theorem branch_probability [Fintype V] (q : Screen E V) (es : List E)
    (M : E → Bool) (ψ : StateVector V) :
    (∑ x, ∑ c : Coins (scan q es x).2, Complex.normSq (kraus q es M ψ x c)) =
      ∑ x, if flags (scan q es x).2 = es.map M then Complex.normSq (ψ x) else 0 := by
  rw [kraus_eq_V_project]
  apply Finset.sum_congr rfl
  intro x _
  rw [Vout_norm, project]
  split <;> simp_all

end
end StrictQuantum
