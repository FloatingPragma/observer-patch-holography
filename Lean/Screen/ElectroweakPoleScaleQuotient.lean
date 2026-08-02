import Mathlib

namespace OPH.ElectroweakPoleScaleQuotient

/-!
# Exact scale quotient of the strict W/Z pole consumer

Write `S = g*v_F/2`, `t = g'/g`, and normalize the declared strict one-loop
self-energy values by the corresponding tree squared masses:

`pW = 1 + Delta_WW(w)/w`,
`pZ = 1 + Delta_ZZ(z)/z`.

The convention-separated strict poles then have the algebraic form

`sW = S^2*pW`,
`sZ = S^2*(1+t^2)*pZ`.

This file proves that the explicit common scale cancels from their ratio at
fixed normalized corrections and constructs an unconstrained algebraic
counterfamily attaining every complex ratio in the factored coordinates.
A passive common-unit rescaling is covered. The theorem is not a total source
derivative under an active vacuum-scale change, since threshold ratios can move
`pW` and `pZ`. The counterfamily is not a physical EFT completion. It proves
that a numerical ratio requires source-selected `t`, `pW`, and `pZ`. No
measured W/Z quantity occurs in this module.
-/

/-- Strict charged pole in the factored consumer coordinates. -/
def chargedPole (S pW : ℂ) : ℂ :=
  S ^ 2 * pW

/-- Strict neutral pole in the factored consumer coordinates. -/
def neutralPole (S t pZ : ℂ) : ℂ :=
  S ^ 2 * (1 + t ^ 2) * pZ

/-- The explicit common pole scale cancels globally from the complex-pole
    ratio at fixed normalized corrections. The remaining inputs are the
    coupling ratio and the two normalized self-energy factors. -/
theorem common_scale_cancels
    (S t pW pZ : ℂ)
    (hS : S ≠ 0) (ht : 1 + t ^ 2 ≠ 0) (hpZ : pZ ≠ 0) :
    chargedPole S pW / neutralPole S t pZ =
      pW / ((1 + t ^ 2) * pZ) := by
  simp only [chargedPole, neutralPole]
  field_simp

/-- Rescaling the common pole scale leaves the complex-pole ratio unchanged. -/
theorem common_rescaling_invariant
    (S t pW pZ scaleFactor : ℂ)
    (hS : S ≠ 0) (hScaleFactor : scaleFactor ≠ 0)
    (ht : 1 + t ^ 2 ≠ 0) (hpZ : pZ ≠ 0) :
    chargedPole (scaleFactor * S) pW /
        neutralPole (scaleFactor * S) t pZ =
      chargedPole S pW / neutralPole S t pZ := by
  rw [common_scale_cancels (scaleFactor * S) t pW pZ
    (mul_ne_zero hScaleFactor hS) ht hpZ]
  rw [common_scale_cancels S t pW pZ hS ht hpZ]

/-- With vanishing normalized self-energy corrections, the quotient reduces
    to the tree relation. This is a conditional consumer identity, not a
    physical prediction of `t`. -/
theorem born_ratio
    (S t : ℂ) (hS : S ≠ 0) (ht : 1 + t ^ 2 ≠ 0) :
    chargedPole S 1 / neutralPole S t 1 = 1 / (1 + t ^ 2) := by
  simpa using common_scale_cancels S t 1 1 hS ht one_ne_zero

/-- The open normalized W coordinate can realize any desired algebraic ratio
    while `t` and the normalized Z coordinate are held fixed. This theorem is
    on the unconstrained factored coordinates. It does not assert that the
    chosen correction passes the decaying-sheet or perturbative-domain guards,
    or is emitted by an OPH source or a physical quantum field theory. -/
theorem normalized_w_correction_counterfamily
    (S t pZ target : ℂ)
    (hS : S ≠ 0) (ht : 1 + t ^ 2 ≠ 0) (hpZ : pZ ≠ 0) :
    chargedPole S (target * ((1 + t ^ 2) * pZ)) /
        neutralPole S t pZ = target := by
  rw [common_scale_cancels S t
    (target * ((1 + t ^ 2) * pZ)) pZ hS ht hpZ]
  field_simp

/-! ## Convention-separated mass and width output vector -/

/-- The dimensionless mass-width vocabulary of the strict pole consumer.
Residues and asymmetries are not members of this output type. -/
noncomputable def massWidthVector
    (massW widthW massZ widthZ : ℝ) : ℝ × ℝ × ℝ :=
  (massW / massZ, widthW / massW, widthZ / massZ)

/-- One passive common change of positive dimensional units cancels globally
from all three strict mass-width ratios.  Nonzero scale is enough for the
algebra; positivity is the physical unit convention. -/
theorem common_scale_cancels_massWidthVector
    {scale massW widthW massZ widthZ : ℝ}
    (hScale : scale ≠ 0) (hMassW : massW ≠ 0) (hMassZ : massZ ≠ 0) :
    massWidthVector
      (scale * massW) (scale * widthW)
      (scale * massZ) (scale * widthZ) =
    massWidthVector massW widthW massZ widthZ := by
  unfold massWidthVector
  apply Prod.ext
  · field_simp [hScale, hMassZ]
  · apply Prod.ext
    · field_simp [hScale, hMassW]
    · field_simp [hScale, hMassZ]

#print axioms common_scale_cancels
#print axioms common_rescaling_invariant
#print axioms born_ratio
#print axioms normalized_w_correction_counterfamily
#print axioms common_scale_cancels_massWidthVector

end OPH.ElectroweakPoleScaleQuotient
