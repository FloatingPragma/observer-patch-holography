import Mathlib.Data.Complex.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Logic.Relation

/-! Endpoint operators on the connected, contrast-zero blank ground sector.
The full bond proof supplies the edge-mismatch support condition. Locality
here means a kernel vanishes between records differing on an untouched site. -/

namespace OPH.BlankChainEndpoint

variable {V : Type*}
abbrev Load := Fin 3

theorem along {R : V → V → Prop} {w : V → Load}
    (agree : ∀ ⦃a b⦄, R a b → w a = w b) {u v : V}
    (path : Relation.ReflTransGen R u v) : w u = w v := by
  induction path with
  | refl => rfl
  | tail _ e ih => exact ih.trans (agree e)

theorem nonzero_record_is_constant {R : V → V → Prop}
    (connected : ∀ u v, Relation.ReflTransGen R u v) (base : V)
    (amp : (V → Load) → ℂ) (w : V → Load)
    (edgeZero : ∀ ⦃a b⦄, R a b → w a ≠ w b → amp w = 0)
    (occupied : amp w ≠ 0) : ∃ q, ∀ v, w v = q := by
  refine ⟨w base, fun v => along ?_ (connected v base)⟩
  intro a b e
  by_contra unequal
  exact occupied (edgeZero e unequal)

def Uniform (q : Load) : V → Load := fun _ => q

theorem untouched_diagonal (site : V) (op : (V → Load) → (V → Load) → ℂ)
    (locality : ∀ u v, u site ≠ v site → op u v = 0)
    {a b : Load} (distinct : a ≠ b) : op (Uniform a) (Uniform b) = 0 :=
  locality _ _ distinct

noncomputable def compressed (op : (V → Load) → (V → Load) → ℂ)
    (psi : Load → ℂ) (a : Load) : ℂ :=
  ∑ b, op (Uniform a) (Uniform b) * psi b

theorem compressed_apply (site : V) (op : (V → Load) → (V → Load) → ℂ)
    (locality : ∀ u v, u site ≠ v site → op u v = 0)
    (psi : Load → ℂ) (a : Load) :
    compressed op psi a = op (Uniform a) (Uniform a) * psi a := by
  unfold compressed
  apply Finset.sum_eq_single a
  · intro b _ hb
    rw [untouched_diagonal site op locality (Ne.symm hb), zero_mul]
  · simp

theorem local_orders_agree (left right : V) (a b : (V → Load) → (V → Load) → ℂ)
    (ha : ∀ u v, u left ≠ v left → a u v = 0)
    (hb : ∀ u v, u right ≠ v right → b u v = 0)
    (psi : Load → ℂ) :
    compressed a (compressed b psi) = compressed b (compressed a psi) := by
  funext q
  simp only [compressed_apply left a ha, compressed_apply right b hb]
  exact mul_left_comm _ _ _

end OPH.BlankChainEndpoint
