import Mathlib.Algebra.Module.LinearMap.Defs
import Mathlib.Data.Complex.Basic

/-! An observer's receiver may read a coherent sum of two energy components.
At a relative half period, opposite half-ring translation characters turn
the free evolution into translation on every accepted input. -/

namespace OPH.InteractingChain

variable {V : Type*} [AddCommGroup V] [Module ℂ V]

/-- The common phase multiplies the translated sum for arbitrary components,
including components of an input correlated with a reference memory. -/
theorem halfBeat_translation (flight shift : V →ₗ[ℂ] V) (phase : ℂ) (even odd : V)
    (flightEven : flight even = phase • even) (flightOdd : flight odd = -(phase • odd))
    (shiftEven : shift even = even) (shiftOdd : shift odd = -odd) :
    flight (even + odd) = phase • shift (even + odd) := by
  simp [map_add, flightEven, flightOdd, shiftEven, shiftOdd, smul_add]

end OPH.InteractingChain
