import Mathlib
import PortFrameGram
import A5FamilyBand

namespace OPH.PortGramRepairBand

open OPH.A5FamilyBand

/-!
# The port Gram form is the low-cost repair band

`PortFrameGram.g5` is the exact pair-valued table for `5 G`, where `G`
is the normalized twelve-port Gram matrix.  `A5FamilyBand.Xt` and `Yt`
are the exact tables for `40 P3 = X + sqrt(5) Y`.  Entrywise,

`2 (5 G) = X + sqrt(5) Y = 40 P3`,

so `G = 4 P3`.  The Galois-conjugated Gram table similarly equals
`4 P3'`.  This file checks those identities directly against the pinned
tables and records their exact action under the declared repair operator
`L = 5 I - A`:

* `P3` has cost `5 - sqrt(5)`;
* the quintet has cost `6`;
* `P3'` has cost `5 + sqrt(5)`.

Thus the actual positive port Gram form occupies the unique lowest positive
eigenband of this repair operator; its Galois conjugate is the maximal-cost
faithful band.  The conclusion is conditional on using the pinned port Gram
form and `5 I - A` as the relevant response data.  It does not derive a
physical kinetic action, family attachment, scale, or laboratory observable.
-/

/-- Twice the exact `5 G` table. -/
def twiceG5 (i j : Fin 12) : ℤ × ℤ :=
  (2 * (OPH.PortFrameGram.g5 i j).1,
    2 * (OPH.PortFrameGram.g5 i j).2)

/-- Twice the Galois conjugate of the exact `5 G` table. -/
def twiceG5Conjugate (i j : Fin 12) : ℤ × ℤ :=
  (2 * (OPH.PortFrameGram.g5 i j).1,
    -2 * (OPH.PortFrameGram.g5 i j).2)

/-- Exact table for `40 P3`. -/
def scaledP3 (i j : Fin 12) : ℤ × ℤ :=
  (entry Xt i j, entry Yt i j)

/-- Exact table for `40 P3'`. -/
def scaledP3Prime (i j : Fin 12) : ℤ × ℤ :=
  (entry Xt i j, -entry Yt i j)

set_option maxHeartbeats 4000000 in
/-- `2 (5 G) = 40 P3`, equivalently `G = 4 P3`, entry for entry. -/
theorem portGram_eq_four_p3 :
    ∀ i j : Fin 12, twiceG5 i j = scaledP3 i j := by
  decide

set_option maxHeartbeats 4000000 in
/-- The Galois-conjugated identity `2 sigma(5 G) = 40 P3'`. -/
theorem portGram_galois_eq_four_p3prime :
    ∀ i j : Fin 12, twiceG5Conjugate i j = scaledP3Prime i j := by
  decide

/-! ## Exact `5 I - A` eigenband packet

The two components below encode multiplication in `Z[sqrt(5)]`.  For
example, multiplying `(X,Y)` by the cost `(5,-1)` gives
`(5X-5Y, 5Y-X)`.
-/

set_option maxHeartbeats 4000000 in
/-- The port Gram band is an exact `5 - sqrt(5)` eigenband of `5 I - A`. -/
theorem p3_repair_eigen :
    matEq
        (fun i j => 5 * entry Xt i j - mulEntry At Xt i j)
        (fun i j => 5 * entry Xt i j - 5 * entry Yt i j) = true ∧
      matEq
        (fun i j => 5 * entry Yt i j - mulEntry At Yt i j)
        (fun i j => 5 * entry Yt i j - entry Xt i j) = true := by
  constructor <;> decide

set_option maxHeartbeats 4000000 in
/-- The Galois control is an exact `5 + sqrt(5)` eigenband of `5 I - A`. -/
theorem p3prime_repair_eigen :
    matEq
        (fun i j => 5 * entry Xt i j - mulEntry At Xt i j)
        (fun i j => 5 * entry Xt i j - 5 * entry Yt i j) = true ∧
      matEq
        (fun i j => -5 * entry Yt i j + mulEntry At Yt i j)
        (fun i j => entry Xt i j - 5 * entry Yt i j) = true := by
  constructor <;> decide

set_option maxHeartbeats 4000000 in
/-- The trivial and quintet controls have repair costs zero and six. -/
theorem outer_repair_eigen :
    matEq
        (fun i j => 5 * entry M1t i j - mulEntry At M1t i j)
        (fun _ _ => 0) = true ∧
      matEq
        (fun i j => 5 * entry M5t i j - mulEntry At M5t i j)
        (fun i j => 6 * entry M5t i j) = true := by
  constructor <;> decide

/-- The Gram band is strictly positive and uniquely below the other faithful
repair bands.  The Galois band lies strictly above the quintet and is therefore
the maximal-cost faithful control. -/
theorem portGram_unique_lowest_positive_galois_maximal :
    z5lt (0, 0) (5, -1) = true ∧
      z5lt (5, -1) (6, 0) = true ∧
      z5lt (6, 0) (5, 1) = true := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

/-- One theorem packet for the exact, non-arbitrary branch identification. -/
theorem portGram_repair_band_packet :
    (∀ i j : Fin 12, twiceG5 i j = scaledP3 i j) ∧
      (∀ i j : Fin 12, twiceG5Conjugate i j = scaledP3Prime i j) ∧
      z5lt (0, 0) (5, -1) = true ∧
      z5lt (5, -1) (6, 0) = true ∧
      z5lt (6, 0) (5, 1) = true := by
  exact ⟨portGram_eq_four_p3, portGram_galois_eq_four_p3prime,
    portGram_unique_lowest_positive_galois_maximal⟩

/-- The band selected by the complete faithful in-window cost theorem is the
same band carried by the pinned port Gram form. -/
theorem selected_family_band_is_port_gram :
    (∀ i j : Fin 12, twiceG5 i j = scaledP3 i j) ∧
      (candidates.filter admissible).length = 3 ∧
      ((candidates.filter admissible).all fun c =>
        c.1 == "3" || z5lt (5, -1) c.2.2.1) = true := by
  exact ⟨portGram_eq_four_p3, family_band_selected⟩

#print axioms portGram_eq_four_p3
#print axioms portGram_galois_eq_four_p3prime
#print axioms p3_repair_eigen
#print axioms p3prime_repair_eigen
#print axioms outer_repair_eigen
#print axioms portGram_unique_lowest_positive_galois_maximal
#print axioms portGram_repair_band_packet
#print axioms selected_family_band_is_port_gram

end OPH.PortGramRepairBand
