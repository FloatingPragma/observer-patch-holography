import Mathlib

namespace OPH.Thermodynamics

/-!
# The exact first-law bookkeeping identity

For internal energy `U(ρ, H) = Re tr(ρ H)`, the change across a joint
update of state and Hamiltonian splits exactly into the heat increment
`Re tr(dρ H)`, the work increment `Re tr(ρ dH)`, and one explicit
bilinear cross term `Re tr(dρ dH)`. The identity is algebraic. A
differentiable protocol sends the cross term to zero at first order; a
finite ordered stroke keeps it explicit, so no silent path choice enters
the heat/work split. Physical content additionally requires the
energy-clock identification receipt: an identified physical `H` and a
protocol distinguishing controlled changes of `H` from state changes.
Dynamical conservation of a fibre-constant energy under the repair
kernel is `heatBath_fixes_fiberObservable` in the conditional-repair
module.
-/

variable {n : ℕ}

/-- Internal energy of a state against an energy observable. -/
noncomputable def internalEnergy (ρ H : Matrix (Fin n) (Fin n) ℂ) : ℝ :=
  (Matrix.trace (ρ * H)).re

/-- Heat increment: state change at fixed energy observable. -/
noncomputable def heatIncrement (dρ H : Matrix (Fin n) (Fin n) ℂ) : ℝ :=
  (Matrix.trace (dρ * H)).re

/-- Work increment: controlled change of the energy observable. -/
noncomputable def workIncrement (ρ dH : Matrix (Fin n) (Fin n) ℂ) : ℝ :=
  (Matrix.trace (ρ * dH)).re

/-- **Exact first law.** The energy change across a joint update splits
into heat, work, and the explicit bilinear cross term. -/
theorem firstLaw_split (ρ dρ H dH : Matrix (Fin n) (Fin n) ℂ) :
    internalEnergy (ρ + dρ) (H + dH) - internalEnergy ρ H
      = heatIncrement dρ H + workIncrement ρ dH
        + (Matrix.trace (dρ * dH)).re := by
  unfold internalEnergy heatIncrement workIncrement
  rw [add_mul, mul_add, mul_add, Matrix.trace_add, Matrix.trace_add,
    Matrix.trace_add]
  simp [Complex.add_re]
  ring

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.firstLaw_split
