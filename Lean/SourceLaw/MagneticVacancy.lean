import Mathlib.Data.Fin.VecNotation
import Mathlib.Tactic

/-!
# A charged orbital in one prescribed unit of face flux

The integer tables describe one perpendicular half-turn pair in the root-zero
gauge on twelve material ports. Each table stores the exponent of a twentieth
root of unity. The executable derives all sixty lifts from the link phases.
The trace identity reduces the finite rotation average to three.
-/

namespace SourceLaw.MagneticVacancy

def turnR : Fin 12 → Fin 12 := ![1,0,3,2,6,7,4,5,11,10,9,8]
def turnS : Fin 12 → Fin 12 := ![2,3,0,1,7,6,5,4,9,8,11,10]
def phaseR : Fin 12 → Fin 20 := ![0,10,5,5,6,6,4,4,6,7,3,4]
def phaseS : Fin 12 → Fin 20 := ![0,15,6,11,9,9,17,17,15,11,0,6]

/-- On every sharp vacancy address, the two orders have the same destination
and differ by ten twentieths of a full phase cycle. -/
theorem order_phase :
    ∀ i, turnR (turnS i) = turnS (turnR i) ∧
      phaseS i + phaseR (turnS i) = phaseR i + phaseS (turnR i) + 10 := by
  decide

/-- Averaging the squared character with the enumerated coefficients gives
three commuting directions. Here w is the fourth power of the primitive
twentieth root used by the link instrument. -/
theorem character_average {w : ℂ} (hw : w ^ 5 = 1) (hn : w ≠ 1) :
    (192 + 12 * (w + w^2 + w^3 + w^4)) / 60 = 3 := by
  have hsum : 1 + w + w^2 + w^3 + w^4 = 0 := by
    have hz : (w - 1) * (1 + w + w^2 + w^3 + w^4) = 0 := by
      calc
        _ = w ^ 5 - 1 := by ring
        _ = 0 := by rw [hw]; ring
    exact (mul_eq_zero.mp hz).resolve_left (sub_ne_zero.mpr hn)
  linear_combination hsum / 5

#print axioms order_phase
#print axioms character_average

end SourceLaw.MagneticVacancy
