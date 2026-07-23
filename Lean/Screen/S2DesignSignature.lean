import Mathlib

namespace OPH.S2DesignSignature

/-! This file checks the degree-six Legendre arithmetic for an assumed
icosahedral inner-product multiset.  It does not select that geometry. -/

/-- The sixth Legendre polynomial, expressed as a polynomial in `x²`.
This lets the icosahedral inner product be certified using only `x² = 1/5`. -/
def legendreSixFromSquare (x2 : ℚ) : ℚ :=
  (231 * x2 ^ 3 - 315 * x2 ^ 2 + 105 * x2 - 5) / 16

/-- At either icosahedral off-diagonal inner product `x = ±1/√5`,
the exact sixth Legendre value is `41/125`. -/
theorem legendreSix_at_icosahedral_square :
    legendreSixFromSquare (1 / 5) = 41 / 125 := by
  norm_num [legendreSixFromSquare]

/-- For any fixed icosahedron vertex, the twelve dot products consist of
one `1`, one `-1`, five `1/√5`, and five `-1/√5`. Since `P₆` is even,
its normalized row (and hence all-pairs) average is `(2 + 10 P₆(1/√5))/12`. -/
def normalizedPairAverageSix : ℚ :=
  (2 + 10 * legendreSixFromSquare (1 / 5)) / 12

theorem normalizedPairAverageSix_eq :
    normalizedPairAverageSix = 11 / 25 := by
  norm_num [normalizedPairAverageSix, legendreSixFromSquare]

/-- The same arithmetic exposed as a one-line consequence of the exact
Legendre value, keeping the combinatorial multiplicities visible. -/
theorem pair_average_from_legendre_value
    (h : legendreSixFromSquare (1 / 5) = 41 / 125) :
    (2 + 10 * legendreSixFromSquare (1 / 5)) / 12 = 11 / 25 := by
  rw [h]
  norm_num

#print axioms legendreSix_at_icosahedral_square
#print axioms normalizedPairAverageSix_eq
#print axioms pair_average_from_legendre_value

end OPH.S2DesignSignature
