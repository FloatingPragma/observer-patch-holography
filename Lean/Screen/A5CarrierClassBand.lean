import Mathlib

namespace OPH.A5CarrierClassBand

/-!
# Rational skeleton of the carrier-class dispersion band

The carrier-class dispersion certificate
(`code/a5_fingerprint/carrier_class_dispersion_certificate.py`, receipt
`runtime/carrier_class_dispersion_receipt.json`) proves that every member of
the declared positive-weight scalar cosine class shares one exact test
surface: a negative quartic coefficient, the isotropic floor `B0/C4^2 >= 10/21`
saturated exactly on single-radius members, empty anisotropic ranks one
through five, and the rank-six band
`B6/B0 = (16/75) <I6(seed)> in [-16/135, 16/75]`.

This file checks the rational arithmetic of that surface: the
kernel-forced orbit multiples, the three-branch coefficient table with the
frozen FZ-11 and FZ-12 rays as the vertex and edge points, the band
endpoint and interior values, the tuned vertex-face zero, the two-shell
Lagrange identity with its sign and equality case, the general-member
control ratio, and the three-value weighted-mean band confinement.

Imported premises, not established here:

* the kernel factorization
  `sum_g ((g u).n)^6 = 60/7 + (64/35) I6(u) I6(n)` and its completeness
  argument, machine-checked in exact arithmetic by the Python certificate;
* the degree-eight factorization
  `sum_g ((g u).n)^8 = 20/3 + (256/75) I6(u) I6(n)` and the vanishing
  harmonic angular-rank-eight multiplicity, with even-rank tuple
  `(m0,m2,m4,m6,m8) = (1,0,0,1,0)`, likewise checked by the exact certificate;
* the seed values `I6 = 1, -5/16, -5/9` on the vertex, edge, and face
  orbits and the range `[-5/9, 1]` from the 62-direction stationary
  census of the fixed-point packet;
* every physical sector, frame, scale, and exclusivity premise of the
  frozen branch predictions. The target-clean producer and physical bridge
  are scoped to the electromagnetism lane (#733); validation follows the
  owning lane's preregistered comparison contract.

The arbitrary finite-mixture formulas for `D0` and `D6` are proved in the
written invariant argument and checked by the exact Python certificate. This
file remains the rational coefficient skeleton rather than a general finite-sum
formalization.

No measured value appears in this file.
-/

/-- Kernel-forced orbit-sum rank-six multiple: `(|O|/60)(64/35) I6(seed)`. -/
def orbitMultiple (size i6seed : ℚ) : ℚ :=
  size / 60 * (64 / 35) * i6seed

/-- Vertex orbit multiple `64/175`. -/
theorem orbitMultiple_vertex : orbitMultiple 12 1 = 64 / 175 := by
  norm_num [orbitMultiple]

/-- Edge orbit multiple `-2/7`, the moment stated in the FZ-12 ray. -/
theorem orbitMultiple_edge : orbitMultiple 30 (-5 / 16) = -2 / 7 := by
  norm_num [orbitMultiple]

/-- Face orbit multiple `-64/189`. -/
theorem orbitMultiple_face : orbitMultiple 20 (-5 / 9) = -64 / 189 := by
  norm_num [orbitMultiple]

/-- Single-orbit branch `B6` in units of `a^4`: `beta / (120 |O|)`. -/
def b6OverA4 (size beta : ℚ) : ℚ :=
  beta / (120 * size)

/-- Vertex branch `B6 = 2 a^4/7875`, the frozen FZ-11 value. -/
theorem b6_vertex : b6OverA4 12 (64 / 175) = 2 / 7875 := by
  norm_num [b6OverA4]

/-- Edge branch `B6 = -a^4/12600`, the frozen FZ-12 value. -/
theorem b6_edge : b6OverA4 30 (-2 / 7) = -1 / 12600 := by
  norm_num [b6OverA4]

/-- Face branch `B6 = -2 a^4/14175`, the completed third point. -/
theorem b6_face : b6OverA4 20 (-64 / 189) = -2 / 14175 := by
  norm_num [b6OverA4]

/-- Single-orbit rank-six-to-isotropic ratio: `7 beta / |O|`. -/
def b6OverB0 (size beta : ℚ) : ℚ :=
  7 * beta / size

/-- The ratio map is `(16/75) I6(seed)` on every orbit. -/
theorem b6OverB0_eq_scaled_seed (size i6seed : ℚ) (hsize : size ≠ 0) :
    b6OverB0 size (orbitMultiple size i6seed) = 16 / 75 * i6seed := by
  unfold b6OverB0 orbitMultiple
  field_simp
  ring

/-- Vertex point `16/75`, the FZ-11 band endpoint. -/
theorem band_vertex : b6OverB0 12 (64 / 175) = 16 / 75 := by
  norm_num [b6OverB0]

/-- Edge point `-1/15`, the FZ-12 interior value. -/
theorem band_edge : b6OverB0 30 (-2 / 7) = -1 / 15 := by
  norm_num [b6OverB0]

/-- Face point `-16/135`, the other band endpoint. -/
theorem band_face : b6OverB0 20 (-64 / 189) = -16 / 135 := by
  norm_num [b6OverB0]

/-- Band endpoints from the census range `[-5/9, 1]`. -/
theorem band_endpoints :
    (16 : ℚ) / 75 * (-5 / 9) = -16 / 135 ∧ (16 : ℚ) / 75 * 1 = 16 / 75 := by
  norm_num

/-- The 25:27 vertex-face mixture cancels the rank-six content. -/
theorem tuned_zero :
    25 * (64 / 175 : ℚ) + 27 * (-64 / 189) = 0 := by
  norm_num

/-- Two-shell Lagrange identity for the radial moments
    `mu2 mu6 - mu4^2 = w1 w2 s t (s - t)^2` with `s, t` the squared
    radii. -/
theorem lagrange_two_shell (w1 w2 s t : ℚ) :
    (w1 * s + w2 * t) * (w1 * s ^ 3 + w2 * t ^ 3)
      - (w1 * s ^ 2 + w2 * t ^ 2) ^ 2
      = w1 * w2 * s * t * (s - t) ^ 2 := by
  ring

/-- Positive weights and radii give the nonnegative moment gap, so the
    isotropic ratio sits on or above the floor. -/
theorem moment_gap_nonneg (w1 w2 s t : ℚ)
    (hw1 : 0 < w1) (hw2 : 0 < w2) (hs : 0 < s) (ht : 0 < t) :
    0 ≤ (w1 * s + w2 * t) * (w1 * s ^ 3 + w2 * t ^ 3)
      - (w1 * s ^ 2 + w2 * t ^ 2) ^ 2 := by
  rw [lagrange_two_shell]
  positivity

/-- Floor saturation forces one radius: a zero moment gap with positive
    weights and radii gives equal squared radii. -/
theorem gap_zero_iff_single_radius (w1 w2 s t : ℚ)
    (hw1 : 0 < w1) (hw2 : 0 < w2) (hs : 0 < s) (ht : 0 < t) :
    (w1 * s + w2 * t) * (w1 * s ^ 3 + w2 * t ^ 3)
        - (w1 * s ^ 2 + w2 * t ^ 2) ^ 2 = 0
      ↔ s = t := by
  rw [lagrange_two_shell]
  constructor
  · intro h
    have hprod : w1 * w2 * s * t ≠ 0 := by positivity
    have hsq : (s - t) ^ 2 = 0 := by
      rcases mul_eq_zero.mp h with h' | h'
      · exact absurd h' hprod
      · exact h'
    have := pow_eq_zero_iff (n := 2) (by norm_num) |>.mp hsq
    linarith [sub_eq_zero.mp this]
  · intro h
    subst h
    ring

/-- General-member control: the three-shell mixed-orbit member of the
    certificate evaluates to `-383632/5682975` by the direct route. -/
theorem general_member_direct :
    7 * ((1 : ℚ) * (64 / 175) * 1 ^ 3 + 1 * (-2 / 7) * 4 ^ 3
        + 3 / 7 * (-64 / 189) * (9 / 4) ^ 3)
      / (1 * 12 * 1 ^ 3 + 1 * 30 * 4 ^ 3 + 3 / 7 * 20 * (9 / 4) ^ 3)
      = -383632 / 5682975 := by
  norm_num

/-- General-member control: the weighted-mean route gives the same value. -/
theorem general_member_mean :
    (16 : ℚ) / 75
      * (1 * 12 * 1 ^ 3 * 1 + 1 * 30 * 4 ^ 3 * (-5 / 16)
        + 3 / 7 * 20 * (9 / 4) ^ 3 * (-5 / 9))
      / (1 * 12 * 1 ^ 3 + 1 * 30 * 4 ^ 3 + 3 / 7 * 20 * (9 / 4) ^ 3)
      = -383632 / 5682975 := by
  norm_num

/-- The control value sits inside the band. -/
theorem general_member_in_band :
    (-16 : ℚ) / 135 ≤ -383632 / 5682975
      ∧ (-383632 : ℚ) / 5682975 ≤ 16 / 75 := by
  norm_num

/-- Three-value weighted-mean confinement: nonnegative weights with
    positive total and values inside `[lo, hi]` keep the mean inside
    `[lo, hi]`.  This is the band argument for three-orbit members. -/
theorem weighted_mean_three_in_band (a b c x y z lo hi : ℚ)
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : 0 ≤ c) (htot : 0 < a + b + c)
    (hx : lo ≤ x ∧ x ≤ hi) (hy : lo ≤ y ∧ y ≤ hi)
    (hz : lo ≤ z ∧ z ≤ hi) :
    lo ≤ (a * x + b * y + c * z) / (a + b + c)
      ∧ (a * x + b * y + c * z) / (a + b + c) ≤ hi := by
  constructor
  · rw [le_div_iff₀ htot]
    nlinarith [mul_le_mul_of_nonneg_left hx.1 ha,
      mul_le_mul_of_nonneg_left hy.1 hb,
      mul_le_mul_of_nonneg_left hz.1 hc]
  · rw [div_le_iff₀ htot]
    nlinarith [mul_le_mul_of_nonneg_left hx.2 ha,
      mul_le_mul_of_nonneg_left hy.2 hb,
      mul_le_mul_of_nonneg_left hz.2 hc]

/-- Eighth-order ratio map: `9 m8 / |O|` with
    `m8 = (|O|/60)(256/75) I6(seed)` the degree-eight orbit multiple. -/
def d6OverD0 (size m8 : ℚ) : ℚ :=
  9 * m8 / size

/-- Degree-eight orbit multiple. -/
def orbitMultiple8 (size i6seed : ℚ) : ℚ :=
  size / 60 * (256 / 75) * i6seed

/-- The eighth-order ratio map is `(64/125) I6(seed)` on every orbit. -/
theorem d6OverD0_eq_scaled_seed (size i6seed : ℚ) (hsize : size ≠ 0) :
    d6OverD0 size (orbitMultiple8 size i6seed) = 64 / 125 * i6seed := by
  unfold d6OverD0 orbitMultiple8
  field_simp
  ring

/-- Vertex eighth-order point `64/125`. -/
theorem eighth_vertex : d6OverD0 12 (orbitMultiple8 12 1) = 64 / 125 := by
  norm_num [d6OverD0, orbitMultiple8]

/-- Edge eighth-order point `-4/25`. -/
theorem eighth_edge :
    d6OverD0 30 (orbitMultiple8 30 (-5 / 16)) = -4 / 25 := by
  norm_num [d6OverD0, orbitMultiple8]

/-- Face eighth-order point `-64/225`, the other endpoint. -/
theorem eighth_face :
    d6OverD0 20 (orbitMultiple8 20 (-5 / 9)) = -64 / 225 := by
  norm_num [d6OverD0, orbitMultiple8]

/-- The cross-order lock: the eighth-order and sixth-order ratio maps
    differ by the universal factor `12/5` at every seed. -/
theorem cross_order_lock (size i6seed : ℚ) (hsize : size ≠ 0) :
    d6OverD0 size (orbitMultiple8 size i6seed)
      = 12 / 5 * b6OverB0 size (orbitMultiple size i6seed) := by
  rw [d6OverD0_eq_scaled_seed size i6seed hsize,
    b6OverB0_eq_scaled_seed size i6seed hsize]
  ring

/-- Genuinely division-free form of the single-radius cross-order lock.
    If the rank-six and rank-eight anisotropies are the certified multiples
    of their isotropic coefficients, the polynomial relation remains
    meaningful even when both anisotropies vanish. -/
theorem cross_order_polynomial (b0 d0 i6seed : ℚ) :
    5 * ((64 / 125 * i6seed) * d0) * b0
      = 12 * ((16 / 75 * i6seed) * b0) * d0 := by
  ring

/-- Exact two-radius control. Equal per-direction vertex and face weights
    at radii one and two violate the single-radius ratio lock by a nonzero
    rational amount. This prevents widening the theorem to all radial
    mixtures. -/
theorem multi_radius_negative_control :
    (64 / 125 : ℚ)
          * (12 * 1 ^ 4 * 1 + 20 * 4 ^ 4 * (-5 / 9))
          / (12 * 1 ^ 4 + 20 * 4 ^ 4)
        - 12 / 5
          * ((16 / 75 : ℚ)
            * (12 * 1 ^ 3 * 1 + 20 * 4 ^ 3 * (-5 / 9))
            / (12 * 1 ^ 3 + 20 * 4 ^ 3))
      = -57344 / 10360225
      ∧ (-57344 : ℚ) / 10360225 ≠ 0 := by
  norm_num

/-- Vertex eighth-order coefficient `-a^6/118125` in units of `a^6`,
    matching the certified through-eighth fixed-point template whose
    `x^8` coefficient is `-2/118125` at symbol normalization two. -/
theorem eighth_vertex_coefficient :
    -(orbitMultiple8 12 1) / (6720 * 12) = -1 / 118125
      ∧ (-2 : ℚ) / 118125 = 2 * (-1 / 118125) := by
  constructor
  · norm_num [orbitMultiple8]
  · norm_num

end OPH.A5CarrierClassBand

#print axioms OPH.A5CarrierClassBand.b6OverB0_eq_scaled_seed
#print axioms OPH.A5CarrierClassBand.moment_gap_nonneg
#print axioms OPH.A5CarrierClassBand.gap_zero_iff_single_radius
#print axioms OPH.A5CarrierClassBand.weighted_mean_three_in_band
#print axioms OPH.A5CarrierClassBand.cross_order_lock
#print axioms OPH.A5CarrierClassBand.cross_order_polynomial
#print axioms OPH.A5CarrierClassBand.multi_radius_negative_control
