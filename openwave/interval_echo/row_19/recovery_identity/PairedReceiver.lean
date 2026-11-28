import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.LinearAlgebra.Dimension.Constructions
import Mathlib.Data.Real.Basic
import Mathlib.Tactic

namespace PairedReceiver

/-- The source coordinates are the nine joint Pauli directions in each of two
parity blocks. A single registered two-qubit output has at most 15 traceless
real coordinates. The channel and its evolution time may be arbitrary. -/
theorem invisible_directions (m : ℕ) (F : (Fin 18 → ℝ) →ₗ[ℝ] (Fin m → ℝ)) :
    18 ≤ m + Module.finrank ℝ (LinearMap.ker F) := by
  have h := LinearMap.finrank_range_add_finrank_ker F
  have hr := Submodule.finrank_le (LinearMap.range F)
  simp at hr h
  omega

theorem paired_head_kernel (F : (Fin 18 → ℝ) →ₗ[ℝ] (Fin 15 → ℝ)) :
    3 ≤ Module.finrank ℝ (LinearMap.ker F) := by
  have := invisible_directions 15 F
  omega

/-- The same output statistics give the same readout even for a nonlinear
processing rule. This statement uses the available state at one station. -/
theorem hidden_in_every_read {m : ℕ} {α : Type*}
    (F : (Fin 18 → ℝ) →ₗ[ℝ] (Fin m → ℝ)) (g : (Fin m → ℝ) → α)
    (x d : Fin 18 → ℝ) (hd : d ∈ LinearMap.ker F) :
    g (F (x + d)) = g (F x) := by
  have h : F d = 0 := hd
  simp [h]

end PairedReceiver
