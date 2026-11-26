import Mathlib

namespace OPH.PublicationTree

variable {I A : Type*} [DecidableEq I]

/-- The addressee takes the retained word of its transmitting neighbour. -/
def send (s : I → A) (p c : I) : I → A := Function.update s c (s p)

def transmit (s : I → A) : List (I × I) → I → A
  | [] => s
  | (p,c) :: rest => transmit (send s p c) rest

def reached (seen : Finset I) : List (I × I) → Finset I
  | [] => seen
  | (_,c) :: rest => reached (insert c seen) rest

def scheduled (seen : Finset I) : List (I × I) → Prop
  | [] => True
  | (p,c) :: rest => p ∈ seen ∧ scheduled (insert c seen) rest

theorem all_readers (steps : List (I × I)) (seen : Finset I) (s : I → A) (a : A)
    (known : ∀ i ∈ seen, s i = a) (valid : scheduled seen steps) :
    ∀ i ∈ reached seen steps, transmit s steps i = a := by
  induction steps generalizing seen s with
  | nil => exact known
  | cons edge rest ih =>
    obtain ⟨p,c⟩ := edge
    obtain ⟨hp, hv⟩ := valid
    apply ih (insert c seen) (send s p c) _ hv
    intro i hi
    by_cases eq : i = c
    · subst i
      simpa [send] using known p hp
    · simpa [send, Function.update_of_ne eq] using known i ((Finset.mem_insert.mp hi).resolve_left eq)

variable {R : Type*} [Fintype A] [DecidableEq A]

/-- Matrix elements after mapping the sharp symbol a to two copies (a,a),
while retaining an arbitrary reference index. -/
def copied (rho : (R × A) → (R × A) → ℂ) (u v : R × A) (b c : A) : ℂ :=
  if u.2 = b then (if v.2 = c then rho u v else 0) else 0

theorem discard_reader (rho : (R × A) → (R × A) → ℂ) (u v : R × A) :
    ∑ b, copied rho u v b b = (if v.2 = u.2 then rho u v else 0) := by
  classical
  simp [copied]

theorem retains_reference_iff (rho : (R × A) → (R × A) → ℂ) :
    (∀ u v, ∑ b, copied rho u v b b = rho u v) ↔
    (∀ u v, v.2 ≠ u.2 → rho u v = 0) := by
  constructor
  · intro h u v huv
    have this := h u v
    simpa [discard_reader, huv] using this.symm
  · intro h u v
    by_cases huv : v.2 = u.2
    · simp [discard_reader, huv]
    · simp [discard_reader, huv, h u v huv]

end OPH.PublicationTree
