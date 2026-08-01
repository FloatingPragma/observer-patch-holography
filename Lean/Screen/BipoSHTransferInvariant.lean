import Mathlib

namespace OPH.BipoSHTransferInvariant

/-!
# Scale-free BipoSH transfer invariant

Let a scalar source field with one copy of each harmonic representation be
carried to an observable sphere with the same multiplicity-one bands by a
linear rotation-equivariant transfer.  Schur's lemma makes the transfer
scalar on each irreducible harmonic band.  Consequently an `ell = 2` by `ell = 4`
covariance block is multiplied by `u2 * u4`, while its two scalar diagonal
normalizations are multiplied by `u2^2` and `u4^2`.

This file proves that the squared amplitude-free statistic is unchanged by
those band scalings.  The theorem removes the transfer amplitudes from the
candidate statistic.  It does not prove that an OPH screen field is a
physical covariance, that a quotient-visible common frame exists, or that a
screen-to-observable map is rotation-equivariant and multiplicity one, or
that a radial transfer has no additional copy-space mixing.  Those remain
physical premises of issue 659.
-/

/-- Squared form of the amplitude-free cross-band statistic.  Here `n` is
the squared norm of the cross-band BipoSH vector, while `d2` and `d4` are the
two nonzero scalar diagonal normalizations. -/
noncomputable def scaleFreeSquared (n d2 d4 : ℝ) : ℝ := n / (d2 * d4)

/-- A nonzero rotation-equivariant scalar transfer on the two harmonic bands
leaves the squared amplitude-free statistic exactly unchanged. -/
theorem scale_free_squared_invariant
    (n d2 d4 u2 u4 : ℝ)
    (hd2 : d2 ≠ 0) (hd4 : d4 ≠ 0)
    (hu2 : u2 ≠ 0) (hu4 : u4 ≠ 0) :
    scaleFreeSquared ((u2 * u4) ^ 2 * n)
        (u2 ^ 2 * d2) (u4 ^ 2 * d4) =
      scaleFreeSquared n d2 d4 := by
  unfold scaleFreeSquared
  field_simp [hd2, hd4, hu2, hu4]

/-- The same cancellation applies when the cross-band norm itself, rather
than its square, is retained and all scale factors are positive. -/
theorem positive_norm_ratio_invariant
    (n d2 d4 u2 u4 : ℝ)
    (hd2 : 0 < d2) (hd4 : 0 < d4)
    (hu2 : 0 < u2) (hu4 : 0 < u4) :
    (u2 * u4 * n) /
        Real.sqrt ((u2 ^ 2 * d2) * (u4 ^ 2 * d4)) =
      n / Real.sqrt (d2 * d4) := by
  rw [show (u2 ^ 2 * d2) * (u4 ^ 2 * d4) =
      (u2 * u4) ^ 2 * (d2 * d4) by ring]
  rw [Real.sqrt_mul (sq_nonneg (u2 * u4))]
  rw [Real.sqrt_sq_eq_abs, abs_of_pos (mul_pos hu2 hu4)]
  field_simp [ne_of_gt hu2, ne_of_gt hu4,
    Real.sqrt_ne_zero'.mpr (mul_pos hd2 hd4)]

#print axioms scale_free_squared_invariant
#print axioms positive_norm_ratio_invariant

end OPH.BipoSHTransferInvariant
