import Mathlib.Tactic

/-! Signed delivery as an observer-record function in the native ternary head.
    The two arguments are the actually written coin and old-load contrast. -/
namespace IntervalEcho.PortPhases

def delivery (coin : Fin 2) (contrast : ℤ) : ℤ :=
  (contrast + (1 - 2 * (coin.val : ℤ)) * (contrast % 2)) / 2

def leftAfter (a b : Fin 3) (coin : Fin 2) : ℤ :=
  ((a.val + b.val) / 2 + coin.val * ((a.val + b.val) % 2) : ℕ)

def rightAfter (a b : Fin 3) (coin : Fin 2) : ℤ :=
  ((a.val + b.val) / 2 + (1 - coin.val) * ((a.val + b.val) % 2) : ℕ)

theorem measuredDelivery (a b : Fin 3) (c : Fin 2) :
    let j := delivery c ((a.val : ℤ) - b.val)
    leftAfter a b c + j = a.val ∧ rightAfter a b c - j = b.val := by
  decide +revert

theorem reversal (c : Fin 2) (k : Fin 5) :
    let d := (k.val : ℤ) - 2
    delivery c d - delivery c (-d) = d := by
  decide +revert

end IntervalEcho.PortPhases
