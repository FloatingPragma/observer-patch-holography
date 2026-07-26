import Mathlib
import PortFrameGram

namespace OPH.A5IncidenceResponse

open OPH.PortFrameGram (adj antipode distance_two_witness)

/-!
# What the twelve-port incidence determines about inverse-port response

This file separates the theorem forced by the certified icosahedral incidence
from the physical interpretation of that theorem.

The incidence has a unique port at graph distance three from each source
(`distance_three_partner_unique`).  Its permutation matrix is exactly the
degree-three adjacency polynomial

`J = (A^3 - 4 A^2 - 5 A + 10 I) / 10`

(`antipode_polynomial`).  Thus a perturb/readback protocol that defines
inverse-port readback as the unique maximal-distance echo computes `J` without
using gauge labels or selecting four sector signs.

BOUNDARY.  Incidence alone does not say that a physical carrier must use this
inverse-port protocol.  The physical response statement additionally requires
the separately named operational response contract.  Nor does this file
identify laboratory gauge currents.
-/

/-- A pair is at graph distance at least three when it is distinct,
non-adjacent, and has no common neighbor. -/
def DistanceAtLeastThree (i j : Fin 12) : Prop :=
  i ≠ j ∧ adj i j = false ∧
    ∀ k : Fin 12, ¬(adj i k = true ∧ adj k j = true)

/-- The incidence determines the inverse port: any port at graph distance at
least three from `i` is the declared antipode of `i`. -/
theorem distance_three_partner_unique (i j : Fin 12)
    (h : DistanceAtLeastThree i j) :
    j = antipode i := by
  by_contra hne
  obtain ⟨k, hik, hkj⟩ := distance_two_witness i j h.1 h.2.1 hne
  exact h.2.2 k ⟨hik, hkj⟩

/-- The integral identity matrix on the twelve ports. -/
def identityMatrix (i j : Fin 12) : ℤ := if i = j then 1 else 0

/-- The integral adjacency matrix of the certified incidence. -/
def adjacencyMatrix (i j : Fin 12) : ℤ := if adj i j then 1 else 0

/-- The permutation matrix of the incidence-determined antipode. -/
def antipodeMatrix (i j : Fin 12) : ℤ := if j = antipode i then 1 else 0

/-- Matrix multiplication for integral kernels on the twelve ports. -/
def matMul (M N : Fin 12 → Fin 12 → ℤ) (i j : Fin 12) : ℤ :=
  ((List.finRange 12).map fun k => M i k * N k j).sum

/-- The square of the adjacency matrix. -/
def adjacencySq : Fin 12 → Fin 12 → ℤ :=
  matMul adjacencyMatrix adjacencyMatrix

/-- The cube of the adjacency matrix. -/
def adjacencyCube : Fin 12 → Fin 12 → ℤ :=
  matMul adjacencySq adjacencyMatrix

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Exact incidence identity for the inverse-port permutation matrix:
`10J = A^3 - 4A^2 - 5A + 10I`. -/
theorem antipode_polynomial :
    ∀ i j : Fin 12,
      10 * antipodeMatrix i j =
        adjacencyCube i j - 4 * adjacencySq i j -
          5 * adjacencyMatrix i j + 10 * identityMatrix i j := by
  decide

end OPH.A5IncidenceResponse

/- Axiom audit: kernel computation with standard axioms only. -/

#print axioms OPH.A5IncidenceResponse.distance_three_partner_unique
#print axioms OPH.A5IncidenceResponse.antipode_polynomial
