import Mathlib
import A5PortModule

namespace OPH.A5PrimitivePortPrediction

/-!
# Primitive twelve-port propagation arithmetic

The proper icosahedral carrier action is transitive on the twelve ports.
Consequently, a scalar coefficient field invariant under that action is
constant.  Once a physical branch declares the complete primitive port orbit
as its sole directional support and normalizes the quadratic continuum term,
the equal-weight cosine symbol has fixed long-wavelength coefficients.

This file proves the transitivity-to-equal-weight implication and the exact
coefficient ratios used by the frozen FZ-11 primitive-port prediction.

Boundary: no theorem here identifies this scalar symbol with a photon or any
other physical sector.  The physical sector, polarization, coherent-frame,
and exclusivity bridges remain premises of the prediction contract and are
owned by issue #655.
-/

/-- Every scalar port coefficient invariant under all sixty proper carrier
rotations is constant on the primitive twelve-port orbit. -/
theorem invariant_port_weights_are_equal (w : Fin 12 → ℚ)
    (h : ∀ g ∈ OPH.A5PortModule.P60, OPH.A5PortModule.pAct g w = w) :
    ∃ c : ℚ, w = c • OPH.A5PortModule.ones := by
  have hw : w ∈ OPH.A5PortModule.Fixed := h
  rw [OPH.A5PortModule.fixed_eq_span, Submodule.mem_span_singleton] at hw
  rcases hw with ⟨c, hc⟩
  exact ⟨c, hc.symm⟩

/-- Coefficient of `a^2 k^4` in the normalized primitive-port symbol. -/
def C4 : ℚ := -1 / 20

/-- Isotropic coefficient of `a^4 k^6`. -/
def B0 : ℚ := 1 / 840

/-- Coefficient of `a^4 k^6 I6`. -/
def B6 : ℚ := 2 / 7875

/-- The anisotropic coefficient is fixed by the square of the lower-order
isotropic coefficient. -/
theorem b6_over_c4_squared : B6 / C4 ^ 2 = 32 / 315 := by
  norm_num [B6, C4]

/-- The isotropic sixth-order coefficient is fixed by the same lower-order
coefficient. -/
theorem b0_over_c4_squared : B0 / C4 ^ 2 = 10 / 21 := by
  norm_num [B0, C4]

/-- The anisotropic and isotropic sixth-order coefficients have one fixed
ratio. -/
theorem b6_over_b0 : B6 / B0 = 16 / 75 := by
  norm_num [B6, B0]

/-- Binary refinement suppresses the primitive-port spin-six coefficient by
sixteen at fixed physical momentum. -/
theorem binary_refinement (a : ℚ) :
    B6 * (a / 2) ^ 4 = (B6 * a ^ 4) / 16 := by
  ring

#print axioms invariant_port_weights_are_equal
#print axioms b6_over_c4_squared
#print axioms b0_over_c4_squared
#print axioms b6_over_b0
#print axioms binary_refinement

end OPH.A5PrimitivePortPrediction
