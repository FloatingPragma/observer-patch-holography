import ObserverPatchHolography.ScalarSeamRepair

/-! # Scalar readback of a resolved half exchange
The terminal records at the two touched addresses have the same distribution.
Its expectation can be evaluated for every assignment of real scores to the
three port bins. A prefix calendar transfers the mass of every nonzero branch
exactly; a zero-mass prefix assigns zero to every allowed child.
-/
namespace OPH.SeamWordMeter
open ObserverPatchHolography.ScalarSeamRepair
open scoped BigOperators
noncomputable section

variable {α : Type*} [Fintype α]
/-- Joint sharp terminal bins before the addressed offer. -/
def firstRead (p : α → α → ℝ) (f : α → ℝ) : ℝ :=
  ∑ a, ∑ b, p a b * f a

def secondRead (p : α → α → ℝ) (f : α → ℝ) : ℝ :=
  ∑ a, ∑ b, p a b * f b

def resolvedRead (p : α → α → ℝ) (f : α → ℝ) : ℝ :=
  ∑ a, ∑ b, p a b * ((f a + f b) / 2)

theorem sharp_score (p : α → α → ℝ) (f : α → ℝ) :
    resolvedRead p f = (firstRead p f + secondRead p f) / 2 := by
  simp only [resolvedRead, firstRead, secondRead, mul_div_assoc', mul_add,
    ← Finset.sum_div, Finset.sum_add_distrib]

variable {ι : Type*} [DecidableEq ι]

theorem calibrate (i j : ι) (q : ι → ℝ) (offset gain : ℝ) :
    pairAverage i j (fun v => offset + gain * q v) =
      fun v => offset + gain * pairAverage i j q v := by
  funext v
  by_cases h : v = i ∨ v = j <;> simp [pairAverage, h]
  ring

def run (w : List (ι × ι)) (q : ι → ℝ) : ι → ℝ :=
  w.foldl (fun x e => pairAverage e.1 e.2 x) q

theorem whole_word (w : List (ι × ι)) (q : ι → ℝ) (offset gain : ℝ) :
    run w (fun v => offset + gain * q v) = fun v => offset + gain * run w q v := by
  induction w generalizing q with
  | nil => rfl
  | cons e rest ih =>
    change run rest (pairAverage e.1 e.2 _) = fun v => offset + gain * run rest (pairAverage e.1 e.2 q) v
    rw [calibrate, ih]

/-- Either output of the calibrated binary contact has its prescribed mass. -/
theorem split_mass (left right : ℝ) (hl : 0 ≤ left) (hr : 0 ≤ right) :
    let p := if left + right = 0 then 0 else right / (left + right)
    (left + right) * (1 - p) = left ∧ (left + right) * p = right := by
  dsimp only
  split_ifs with h
  · constructor <;> linarith
  · constructor <;> field_simp [h]
    ring

def conditional (p : List Bool → ℝ) (w : List Bool) (b : Bool) : ℝ :=
  let l := p (w ++ [false])
  let r := p (w ++ [true])
  let q := if l + r = 0 then 0 else r / (l + r)
  if b then q else 1 - q

theorem transferred (p : List Bool → ℝ) (hn : ∀ w, 0 ≤ p w)
    (mass : ∀ w, p (w ++ [false]) + p (w ++ [true]) = p w)
    (w : List Bool) (b : Bool) : p w * conditional p w b = p (w ++ [b]) := by
  have hs := split_mass _ _ (hn (w ++ [false])) (hn (w ++ [true]))
  rw [← mass w]
  cases b
  · exact hs.1
  · exact hs.2

/-- The source may assign dependent masses to the successive edges of a word. -/
theorem leaf_mass (p : List Bool → ℝ) (hn : ∀ w, 0 ≤ p w)
    (mass : ∀ w, p (w ++ [false]) + p (w ++ [true]) = p w)
    (start : List Bool) (bits : List Bool) :
    (bits.foldl (fun s b => (s.1 ++ [b], s.2 * conditional p s.1 b))
      (start, p start)).2 = p (start ++ bits) := by
  induction bits generalizing start with
  | nil => simp
  | cons b bs ih =>
    simp only [List.foldl_cons]
    rw [transferred p hn mass, ih]
    simp

end
end OPH.SeamWordMeter
#print axioms OPH.SeamWordMeter.sharp_score
#print axioms OPH.SeamWordMeter.whole_word
#print axioms OPH.SeamWordMeter.split_mass
#print axioms OPH.SeamWordMeter.leaf_mass
