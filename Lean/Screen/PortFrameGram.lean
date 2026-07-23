import Mathlib

open scoped BigOperators

namespace OPH.PortFrameGram

/-! # The icosahedral port-frame Gram identity, exactly

This file formalizes the port-frame layer of the twelve-port screen packet
(issue #568): on the explicit icosahedral incidence structure, the declared
Gram rule with entries `1, 1/√5, -1/√5, -1` at graph distances `0,1,2,3`
satisfies `G² = 4·G` and is symmetric with trace `12`.  Arithmetic is exact:
the file works with the scaled matrix `5·G`, whose entries live in `ℤ(√5)`
represented as integer pairs `(a, b) ↦ a + b·√5`, and every identity is
checked by kernel `decide` with no additional axioms.

It also verifies the combinatorial spine those entries rely on: the encoded
graph is 5-regular, the antipode map is a fixed-point-free involution, every
non-adjacent non-antipodal pair has a common neighbor (distance two), and no
port shares a neighbor with its antipode (distance three).

Boundary: this file does not derive the icosahedral carrier from OPH axioms,
does not construct the `A5` action or classify automorphisms, and proves
nothing about physical currents.  Rank three of `G` follows on paper from
`(G/4)² = G/4` and `tr G = 12`; the idempotency and trace facts proved here
are the two inputs to that argument. -/

/-- Neighbor lists of the icosahedral incidence used by the screen packet.
Vertex `i` is antipodal to `11 - i`. -/
def neighbors : Fin 12 → List (Fin 12)
  | 0 => [1, 2, 3, 4, 6]
  | 1 => [0, 2, 3, 5, 7]
  | 2 => [0, 1, 4, 5, 8]
  | 3 => [0, 1, 6, 7, 9]
  | 4 => [0, 2, 6, 8, 10]
  | 5 => [1, 2, 7, 8, 11]
  | 6 => [0, 3, 4, 9, 10]
  | 7 => [1, 3, 5, 9, 11]
  | 8 => [2, 4, 5, 10, 11]
  | 9 => [3, 6, 7, 10, 11]
  | 10 => [4, 6, 8, 9, 11]
  | 11 => [5, 7, 8, 9, 10]

/-- Adjacency relation read off the neighbor lists. -/
def adj (i j : Fin 12) : Bool := j ∈ neighbors i

/-- The antipode of port `i`. -/
def antipode (i : Fin 12) : Fin 12 := 11 - i

/-! ## Combinatorial spine -/

/-- The incidence is symmetric. -/
theorem adj_symm : ∀ i j : Fin 12, adj i j = adj j i := by decide

/-- No port is its own neighbor. -/
theorem adj_irrefl : ∀ i : Fin 12, adj i i = false := by decide

/-- Every port has exactly five neighbors. -/
theorem degree_five :
    ∀ i : Fin 12, (Finset.univ.filter (fun j => adj i j)).card = 5 := by decide

/-- The antipode map is an involution. -/
theorem antipode_involutive : ∀ i : Fin 12, antipode (antipode i) = i := by decide

/-- The antipode map has no fixed point. -/
theorem antipode_fixed_point_free : ∀ i : Fin 12, antipode i ≠ i := by decide

/-- A port is never adjacent to its antipode. -/
theorem antipode_not_adjacent : ∀ i : Fin 12, adj i (antipode i) = false := by decide

/-- A port and its antipode share no neighbor: the antipode sits at graph
distance three, not two. -/
theorem antipode_no_common_neighbor :
    ∀ i : Fin 12, ∀ k : Fin 12, ¬(adj i k = true ∧ adj k (antipode i) = true) := by
  decide

/-- Every non-adjacent, non-antipodal, distinct pair has a common neighbor:
those pairs sit at graph distance exactly two. -/
theorem distance_two_witness :
    ∀ i j : Fin 12, i ≠ j → adj i j = false → j ≠ antipode i →
      ∃ k : Fin 12, adj i k = true ∧ adj k j = true := by decide

/-! ## The Gram identities

The declared Gram rule has entries `1, 1/√5, -1/√5, -1` at graph distances
`0,1,2,3`.  This file works with the scaled matrix `g5 = 5·G`, whose entries
lie in `ℤ(√5)` (pairs `(a,b)` denoting `a + b·√5`), so every arithmetic step
in the load-bearing identity is exact integer arithmetic that the kernel
evaluates directly.  In that normalization `G² = 4·G` reads `g5² = 20·g5`,
and `tr G = 12` reads `tr g5 = 60`. -/

/-- The scaled Gram matrix `5·G`, with integer components in `ℤ(√5)`:
`(5,0)` on the diagonal, `(0,1)` at distance one, `(0,-1)` at distance two,
`(-5,0)` at the antipode. -/
def g5 (i j : Fin 12) : ℤ × ℤ :=
  if i = j then (5, 0)
  else if adj i j then (0, 1)
  else if j = antipode i then (-5, 0)
  else (0, -1)

/-- Multiplication in `ℤ(√5)`. -/
def mulZ5 (x y : ℤ × ℤ) : ℤ × ℤ :=
  (x.1 * y.1 + 5 * x.2 * y.2, x.1 * y.2 + x.2 * y.1)

/-- One entry of the matrix product `(5G)·(5G)`. -/
def g5sq (i j : Fin 12) : ℤ × ℤ :=
  (((List.finRange 12).map fun k => mulZ5 (g5 i k) (g5 k j))).foldr
    (fun x s => (x.1 + s.1, x.2 + s.2)) (0, 0)

/-- The Gram matrix is symmetric. -/
theorem g5_symm : ∀ i j : Fin 12, g5 i j = g5 j i := by decide

/-- The trace of `5·G` is `60`, i.e. `tr G = 12`. -/
theorem g5_trace :
    (((List.finRange 12).map fun i => g5 i i)).foldr
      (fun x s => (x.1 + s.1, x.2 + s.2)) (0, 0) = ((60 : ℤ), (0 : ℤ)) := by
  decide

/-- The port-frame Gram identity, exactly: `(5G)² = 20·(5G)`, equivalently
`G² = 4·G`.  With `g5_trace` this gives the rank-three statement: `G/4` is
idempotent with trace `3`, and the rank of an idempotent equals its trace
(that final step is on paper). -/
theorem gram_sq :
    ∀ i j : Fin 12, g5sq i j = (20 * (g5 i j).1, 20 * (g5 i j).2) := by decide

end OPH.PortFrameGram
