import Mathlib

namespace OPH.BipoSHFrameInvariant

/-!
# Frame-quotient invariant for the rank-six BipoSH statistic

The frozen issue-659 statistic uses the norm of the complete `L = 6`
coefficient vector.  A coherent change of angular frame acts by a linear
isometry on that vector, so its norm and the normalized squared statistic do
not depend on the unknown carrier orientation.

This result removes the numerical value of a common global orientation from
the comparison contract.  It does not construct a screen-to-observable map,
prove that the map is equivariant or isometric, exclude radial copy-space
mixing, or identify a physical covariance.  Those remain physical premises.
-/

/-- Squared normalized cross-band statistic.  `cross` is the complete
rank-six coefficient vector and `d2`, `d4` are the two diagonal scalar
normalizations. -/
noncomputable def normalizedCrossPower {V : Type*} [SeminormedAddCommGroup V]
    (cross : V) (d2 d4 : ℝ) : ℝ :=
  ‖cross‖ ^ 2 / (d2 * d4)

/-- The frozen unsquared statistic used by the certificate surface.  The
absolute value accommodates the certificate's fail-closed real convention;
physical covariance diagonal blocks are positive on an attached branch. -/
noncomputable def normalizedCrossAmplitude {V : Type*}
    [SeminormedAddCommGroup V] (cross : V) (d2 d4 : ℝ) : ℝ :=
  ‖cross‖ / Real.sqrt |d2 * d4|

/-- A coherent frame isometry preserves the frozen unsquared statistic
exactly, without selecting an absolute carrier orientation. -/
theorem normalizedCrossAmplitude_frame_invariant
    {𝕜 V : Type*} [RCLike 𝕜] [SeminormedAddCommGroup V]
    [NormedSpace 𝕜 V]
    (frame : V ≃ₗᵢ[𝕜] V) (cross : V) (d2 d4 : ℝ) :
    normalizedCrossAmplitude (frame cross) d2 d4 =
      normalizedCrossAmplitude cross d2 d4 := by
  simp [normalizedCrossAmplitude]

/-- Any coherent orthogonal or unitary frame change preserves the statistic
exactly. -/
theorem normalizedCrossPower_frame_invariant
    {𝕜 V : Type*} [RCLike 𝕜] [SeminormedAddCommGroup V]
    [NormedSpace 𝕜 V]
    (frame : V ≃ₗᵢ[𝕜] V) (cross : V) (d2 d4 : ℝ) :
    normalizedCrossPower (frame cross) d2 d4 =
      normalizedCrossPower cross d2 d4 := by
  simp [normalizedCrossPower]

/-- The orientation quotient and the multiplicity-one scalar transfer can be
taken together.  Positive band amplitudes rescale the cross block and its two
diagonal normalizations in the matching way, while the frame isometry leaves
the coefficient norm fixed. -/
theorem normalizedCrossPower_frame_and_band_scale_invariant
    {V : Type*} [SeminormedAddCommGroup V] [NormedSpace ℝ V]
    (frame : V ≃ₗᵢ[ℝ] V) (cross : V)
    (d2 d4 u2 u4 : ℝ)
    (hd2 : d2 ≠ 0) (hd4 : d4 ≠ 0)
    (hu2 : 0 < u2) (hu4 : 0 < u4) :
    normalizedCrossPower ((u2 * u4) • frame cross)
        (u2 ^ 2 * d2) (u4 ^ 2 * d4) =
      normalizedCrossPower cross d2 d4 := by
  unfold normalizedCrossPower
  rw [norm_smul, frame.norm_map, Real.norm_eq_abs,
    abs_of_pos (mul_pos hu2 hu4)]
  field_simp [hd2, hd4, ne_of_gt hu2, ne_of_gt hu4]

/-- Complex transfer phases are harmless as well.  Observable covariance
normalizations acquire squared norms of the two band amplitudes, so every
nonzero complex phase and magnitude cancels from the squared statistic. -/
theorem normalizedCrossPower_complex_band_scale_invariant
    {𝕜 V : Type*} [RCLike 𝕜] [SeminormedAddCommGroup V]
    [NormedSpace 𝕜 V]
    (frame : V ≃ₗᵢ[𝕜] V) (cross : V)
    (d2 d4 : ℝ) (u2 u4 : 𝕜)
    (hd2 : d2 ≠ 0) (hd4 : d4 ≠ 0)
    (hu2 : u2 ≠ 0) (hu4 : u4 ≠ 0) :
    normalizedCrossPower ((u2 * star u4) • frame cross)
        (‖u2‖ ^ 2 * d2) (‖u4‖ ^ 2 * d4) =
      normalizedCrossPower cross d2 d4 := by
  unfold normalizedCrossPower
  rw [norm_smul, frame.norm_map, norm_mul, norm_star]
  field_simp [hd2, hd4, norm_ne_zero_iff.mpr hu2, norm_ne_zero_iff.mpr hu4]

/-- A frame choice can rotate the coefficient vector while leaving its norm
fixed.  Therefore a nonzero source norm cannot be tuned to zero by changing
the coherent global orientation. -/
theorem nonzero_survives_frame_change
    {𝕜 V : Type*} [RCLike 𝕜] [NormedAddCommGroup V]
    [NormedSpace 𝕜 V]
    (frame : V ≃ₗᵢ[𝕜] V) (cross : V) (hcross : cross ≠ 0) :
    frame cross ≠ 0 := by
  simpa using frame.injective.ne hcross

/-- The invariant also survives composition of any finite chain of coherent
frame changes. -/
theorem normalizedCrossPower_two_frames
    {𝕜 V : Type*} [RCLike 𝕜] [SeminormedAddCommGroup V]
    [NormedSpace 𝕜 V]
    (frame₁ frame₂ : V ≃ₗᵢ[𝕜] V) (cross : V) (d2 d4 : ℝ) :
    normalizedCrossPower (frame₂ (frame₁ cross)) d2 d4 =
      normalizedCrossPower cross d2 d4 := by
  rw [normalizedCrossPower_frame_invariant frame₂,
    normalizedCrossPower_frame_invariant frame₁]

#print axioms normalizedCrossPower_frame_invariant
#print axioms normalizedCrossAmplitude_frame_invariant
#print axioms normalizedCrossPower_frame_and_band_scale_invariant
#print axioms normalizedCrossPower_complex_band_scale_invariant
#print axioms nonzero_survives_frame_change
#print axioms normalizedCrossPower_two_frames

end OPH.BipoSHFrameInvariant
