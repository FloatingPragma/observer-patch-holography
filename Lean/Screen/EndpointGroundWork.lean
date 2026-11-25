import Mathlib.Tactic.Linarith
import Mathlib.Algebra.Order.Group.Abs
import Mathlib.Data.Real.Basic

/-!
# Static work with local endpoints on a common state space

The trial states may have arbitrarily many excitations or be entangled with an
environment. A minimum is used only where a state actually attains it, as in
the finite nonempty closed state families of the experiment.
-/
namespace OPH.EndpointGroundWork

structure Minimum {S : Type*} (energy : S → ℝ) (value : ℝ) : Prop where
  lower : ∀ s, value ≤ energy s
  attained : ∃ s, energy s = value

theorem insertion {S : Type*} {e v : S → ℝ} {x y l u : ℝ}
    (a : Minimum e x) (b : Minimum (fun s => e s + v s) y)
    (width : ∀ s, l ≤ v s ∧ v s ≤ u) : l ≤ y-x ∧ y-x ≤ u := by
  obtain ⟨s, hs⟩ := a.attained
  obtain ⟨t, ht⟩ := b.attained
  have ht' : e t + v t = y := ht
  have hs' : y ≤ e s + v s := b.lower s
  constructor
  · linarith [(width t).1, a.lower t]
  · linarith [(width s).2]

/-- Adding one pin in the presence and absence of a second pin on the same
family. The connected ground work cannot exceed the smaller field width. -/
theorem connected_work {S : Type*} {e p q : S → ℝ} {e0 ep eq epq lp up lq uq : ℝ}
    (h0 : Minimum e e0) (hp : Minimum (fun s => e s + p s) ep)
    (hq : Minimum (fun s => e s + q s) eq)
    (hpq : Minimum (fun s => e s + p s + q s) epq)
    (bp : ∀ s, lp ≤ p s ∧ p s ≤ up) (bq : ∀ s, lq ≤ q s ∧ q s ≤ uq) :
    |epq - ep - eq + e0| ≤ min (up-lp) (uq-lq) := by
  have hqp : Minimum (fun s => e s + q s + p s) epq := by
    simpa only [add_right_comm] using hpq
  obtain ⟨a,b⟩ := insertion h0 hp bp
  obtain ⟨c,d⟩ := insertion hq hqp bp
  obtain ⟨f,g⟩ := insertion h0 hq bq
  obtain ⟨i,j⟩ := insertion hp hpq bq
  apply le_min <;> apply abs_le.mpr <;> constructor <;> linarith

/-- A trial-state bijection between two hard pin positions. -/
theorem relocation {S T : Type*} (swap : S ≃ T) {e : S → ℝ} {f : T → ℝ}
    {x y bound : ℝ} (a : Minimum e x) (b : Minimum f y)
    (displacement : ∀ s, |f (swap s) - e s| ≤ bound) : |y-x| ≤ bound := by
  obtain ⟨s, hs⟩ := a.attained
  obtain ⟨t, ht⟩ := b.attained
  have hs' := (abs_le.mp (displacement s)).2
  have ht' := (abs_le.mp (displacement (swap.symm t))).1
  rw [swap.apply_symm_apply] at ht'
  apply abs_le.mpr
  constructor
  · linarith [a.lower (swap.symm t)]
  · linarith [b.lower (swap s)]

end OPH.EndpointGroundWork
