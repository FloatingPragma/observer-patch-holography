import Mathlib

namespace OPH.ScreenCount

abbrev F := ZMod 2

variable {ι : Type} [Fintype ι] [DecidableEq ι]

def digit (i : ι) : ι → F := Pi.single i 1

def moves (G : SimpleGraph ι) : Submodule F (ι → F) :=
  Submodule.span F {v | ∃ i j, G.Adj i j ∧ v = digit i - digit j}

omit [Fintype ι] in
private theorem walk_move (G : SimpleGraph ι) {i j : ι} (p : G.Walk i j) :
    digit i - digit j ∈ moves G := by
  induction p with
  | nil => simp
  | @cons i j k edge _ ih =>
    have h : digit i - digit j ∈ moves G :=
      Submodule.subset_span ⟨i,j,edge,rfl⟩
    convert (moves G).add_mem h ih using 1; abel

/-- One fair orientation bit per bridge supplies its incidence column. A connected
region realizes all and only zero-total changes of its owner parity word. -/
theorem moves_iff (G : SimpleGraph ι) (conn : G.Connected) (root : ι) (v : ι → F) :
    v ∈ moves G ↔ ∑ i, v i = 0 := by
  constructor
  · intro hv
    induction hv using Submodule.span_induction with
    | mem x hx =>
      obtain ⟨a,b,_,rfl⟩ := hx
      simp [digit, Pi.single_apply, Finset.sum_sub_distrib]
    | zero => simp
    | add a b _ _ ha hb => simpa [Finset.sum_add_distrib] using congrArg₂ (·+·) ha hb
    | smul t a _ ha => simp [Pi.smul_apply, smul_eq_mul, ← Finset.mul_sum, ha]
  · intro hz
    have cert : ∑ i, v i • (digit i - digit root) = v := by
      ext j
      simp [digit, Pi.single_apply, Finset.sum_apply, Pi.smul_apply, mul_sub,
        Finset.sum_sub_distrib, hz, mul_ite]
    rw [← cert]
    exact Submodule.sum_mem _ fun i _ => (moves G).smul_mem _ (walk_move G (conn i root).some)

/-- Fix the observed total parity and delete any one owner. Every strict
subscreen word has exactly one completion to the chosen parity. -/
noncomputable def deleteOwner {n : ℕ} (a : Fin (n+1)) (p : F) :
    {w : Fin (n+1) → F // ∑ i, w i = p} ≃ (Fin n → F) where
  toFun w := a.removeNth w.1
  invFun w := ⟨a.insertNth (p - ∑ i, w i) w, by simp [Fin.sum_insertNth]⟩
  right_inv w := by ext j; simp [Fin.removeNth]
  left_inv w := by
    apply Subtype.ext
    have h := Fin.sum_univ_succAbove (f := w.1) a
    have entry : p - ∑ i, a.removeNth w.1 i = w.1 a := by
      rw [w.2] at h
      simpa [Fin.removeNth] using (sub_eq_iff_eq_add.mpr h)
    simp only [entry, Fin.insertNth_self_removeNth]

theorem supportCount (n : ℕ) (p : F) :
    Fintype.card {w : Fin (n+1) → F // ∑ i, w i = p} = 2^n := by
  simpa using Fintype.card_congr (deleteOwner (0 : Fin (n+1)) p)

end OPH.ScreenCount
