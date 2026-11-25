import Mathlib

/-! The charge-raising port instrument has one target rail and any finite
number of address rails. An arriving vacancy selects a phase on the target.
These lemmas relate the three ternary echo outcomes, the coherent sum over
one-vacancy addresses, and the signed regional charge change. -/

namespace SourceLaw.ArrivalChargeLadder

/-- After an echo on (x,2) both output words have the same number of unit
loads. The inverse exposure on the unencoded control removes the x=1 term. -/
theorem ternary_select (x : Fin 3) :
    (2 - (x.val : ℤ)) - (if x.val = 1 then 1 else 0) =
      (if x.val = 0 then 2 else 0) := by
  fin_cases x <;> decide

/-- Every vacancy address receives the same total phase when each port is
visited once with a common, restored control. An unread quantum reference
can be included in the coefficients. -/
theorem one_vacancy_aperture {A : Type*} [Fintype A] [DecidableEq A]
    (port : A) (φ : ℝ) :
    ∑ addr : A, (if addr = port then φ else 0) = φ := by
  simp

inductive Code where
  | charged | neutral | singlet
  deriving DecidableEq

/-- Charge in the address and target region relative to a unit-filled host. -/
def materialCharge : Code → ℤ
  | .charged => -1
  | .neutral => 0
  | .singlet => 1

/-- The oriented, reservoir-resolved conversion acts identically at every address. -/
def raiseRead {A : Type*} (v : Code → A → ℂ) : Code → A → ℂ
  | .neutral, a => v .charged a
  | _, _ => 0

/-- Multiplication by the regional charge on each output record. -/
def chargeRead {A : Type*} (v : Code → A → ℂ) (s : Code) (a : A) : ℂ :=
  (materialCharge s : ℂ) * v s a

/-- The forward branch raises regional charge by one even when the input is
correlated with an arbitrary untouched address or reference. -/
theorem charge_step {A : Type*} (v : Code → A → ℂ) (s : Code) (a : A) :
    chargeRead (raiseRead v) s a - raiseRead (chargeRead v) s a = raiseRead v s a := by
  cases s <;> simp [chargeRead, raiseRead, materialCharge]

/-- The charged singlet input has zero forward amplitude at every output. -/
theorem singlet_zero {A : Type*} (v : A → ℂ) (s : Code) (a : A) :
    raiseRead (fun slot addr => if slot = .singlet then v addr else 0) s a = 0 := by
  cases s <;> simp [raiseRead]

#print axioms ternary_select
#print axioms one_vacancy_aperture
#print axioms charge_step
#print axioms singlet_zero

end SourceLaw.ArrivalChargeLadder
