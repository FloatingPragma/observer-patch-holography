import Mathlib

namespace OPH.A5PortAction

/-! # The icosahedral rotation action on the twelve ports, explicitly

Companion to `PortFrameGram` for issue #568: the sixty proper rotations of
the icosahedron, listed explicitly as permutations of the twelve ports.  The
theorems verify, by kernel `decide` alone, that the listed set is a group of
order sixty (closed, containing the identity, with two-sided inverses inside
the list), that every element preserves the port incidence of
`PortFrameGram.neighbors` and commutes with the antipode `i ↦ 11 - i`, and
that the action is transitive on ports.

Boundary: this file exhibits the order-sixty incidence-automorphism group and
its properties; the identification of that abstract group with `A5`
(simplicity, presentation, or classification) and everything downstream
(current lifts, physical descent) stay on paper or in later files. -/

/-- The sixty proper icosahedral rotations as port permutations
(row `p` sends port `k` to `p.getD k 0`). -/
def perms : List (List Nat) := [
  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
  [0, 2, 4, 1, 6, 8, 3, 5, 10, 7, 9, 11],
  [0, 3, 1, 6, 2, 7, 4, 9, 5, 10, 8, 11],
  [0, 4, 6, 2, 3, 10, 1, 8, 9, 5, 7, 11],
  [0, 6, 3, 4, 1, 9, 2, 10, 7, 8, 5, 11],
  [1, 0, 3, 2, 7, 6, 5, 4, 9, 8, 11, 10],
  [1, 2, 0, 5, 3, 4, 7, 8, 6, 11, 9, 10],
  [1, 3, 7, 0, 5, 9, 2, 6, 11, 4, 8, 10],
  [1, 5, 2, 7, 0, 8, 3, 11, 4, 9, 6, 10],
  [1, 7, 5, 3, 2, 11, 0, 9, 8, 6, 4, 10],
  [2, 0, 1, 4, 5, 3, 8, 6, 7, 10, 11, 9],
  [2, 1, 5, 0, 8, 7, 4, 3, 11, 6, 10, 9],
  [2, 4, 0, 8, 1, 6, 5, 10, 3, 11, 7, 9],
  [2, 5, 8, 1, 4, 11, 0, 7, 10, 3, 6, 9],
  [2, 8, 4, 5, 0, 10, 1, 11, 6, 7, 3, 9],
  [3, 0, 6, 1, 9, 4, 7, 2, 10, 5, 11, 8],
  [3, 1, 0, 7, 6, 2, 9, 5, 4, 11, 10, 8],
  [3, 6, 9, 0, 7, 10, 1, 4, 11, 2, 5, 8],
  [3, 7, 1, 9, 0, 5, 6, 11, 2, 10, 4, 8],
  [3, 9, 7, 6, 1, 11, 0, 10, 5, 4, 2, 8],
  [4, 0, 2, 6, 8, 1, 10, 3, 5, 9, 11, 7],
  [4, 2, 8, 0, 10, 5, 6, 1, 11, 3, 9, 7],
  [4, 6, 0, 10, 2, 3, 8, 9, 1, 11, 5, 7],
  [4, 8, 10, 2, 6, 11, 0, 5, 9, 1, 3, 7],
  [4, 10, 6, 8, 0, 9, 2, 11, 3, 5, 1, 7],
  [5, 1, 7, 2, 11, 3, 8, 0, 9, 4, 10, 6],
  [5, 2, 1, 8, 7, 0, 11, 4, 3, 10, 9, 6],
  [5, 7, 11, 1, 8, 9, 2, 3, 10, 0, 4, 6],
  [5, 8, 2, 11, 1, 4, 7, 10, 0, 9, 3, 6],
  [5, 11, 8, 7, 2, 10, 1, 9, 4, 3, 0, 6],
  [6, 0, 4, 3, 10, 2, 9, 1, 8, 7, 11, 5],
  [6, 3, 0, 9, 4, 1, 10, 7, 2, 11, 8, 5],
  [6, 4, 10, 0, 9, 8, 3, 2, 11, 1, 7, 5],
  [6, 9, 3, 10, 0, 7, 4, 11, 1, 8, 2, 5],
  [6, 10, 9, 4, 3, 11, 0, 8, 7, 2, 1, 5],
  [7, 1, 3, 5, 9, 0, 11, 2, 6, 8, 10, 4],
  [7, 3, 9, 1, 11, 6, 5, 0, 10, 2, 8, 4],
  [7, 5, 1, 11, 3, 2, 9, 8, 0, 10, 6, 4],
  [7, 9, 11, 3, 5, 10, 1, 6, 8, 0, 2, 4],
  [7, 11, 5, 9, 1, 8, 3, 10, 2, 6, 0, 4],
  [8, 2, 5, 4, 11, 1, 10, 0, 7, 6, 9, 3],
  [8, 4, 2, 10, 5, 0, 11, 6, 1, 9, 7, 3],
  [8, 5, 11, 2, 10, 7, 4, 1, 9, 0, 6, 3],
  [8, 10, 4, 11, 2, 6, 5, 9, 0, 7, 1, 3],
  [8, 11, 10, 5, 4, 9, 2, 7, 6, 1, 0, 3],
  [9, 3, 6, 7, 10, 0, 11, 1, 4, 5, 8, 2],
  [9, 6, 10, 3, 11, 4, 7, 0, 8, 1, 5, 2],
  [9, 7, 3, 11, 6, 1, 10, 5, 0, 8, 4, 2],
  [9, 10, 11, 6, 7, 8, 3, 4, 5, 0, 1, 2],
  [9, 11, 7, 10, 3, 5, 6, 8, 1, 4, 0, 2],
  [10, 4, 8, 6, 11, 2, 9, 0, 5, 3, 7, 1],
  [10, 6, 4, 9, 8, 0, 11, 3, 2, 7, 5, 1],
  [10, 8, 11, 4, 9, 5, 6, 2, 7, 0, 3, 1],
  [10, 9, 6, 11, 4, 3, 8, 7, 0, 5, 2, 1],
  [10, 11, 9, 8, 6, 7, 4, 5, 3, 2, 0, 1],
  [11, 5, 7, 8, 9, 1, 10, 2, 3, 4, 6, 0],
  [11, 7, 9, 5, 10, 3, 8, 1, 6, 2, 4, 0],
  [11, 8, 5, 10, 7, 2, 9, 4, 1, 6, 3, 0],
  [11, 9, 10, 7, 8, 6, 5, 3, 4, 1, 2, 0],
  [11, 10, 8, 9, 5, 4, 7, 6, 2, 3, 1, 0]]

/-- Apply a listed permutation. -/
def app (p : List Nat) (k : Nat) : Nat := p.getD k 0

/-- Compose two listed permutations (`p` after `q`). -/
def comp (p q : List Nat) : List Nat := (List.range 12).map fun k => app p (app q k)

/-- Port adjacency at the `Nat` level, matching `PortFrameGram.neighbors`. -/
def natNeighbors : Nat → List Nat
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
  | _ => [5, 7, 8, 9, 10]

/-- Adjacency test. -/
def natAdj (i j : Nat) : Bool := j ∈ natNeighbors i

/-- The list has exactly sixty entries. -/
theorem perms_length : perms.length = 60 := by decide

/-- No permutation is listed twice. -/
theorem perms_nodup : perms.Nodup := by decide

/-- Every row is a genuine permutation of `0..11`. -/
theorem perms_are_permutations :
    perms.all (fun p => p.length = 12 && (List.range 12).all (fun k => p.contains k)) = true := by
  decide

/-- The identity is listed. -/
theorem id_mem : List.range 12 ∈ perms := by decide

set_option maxHeartbeats 4000000 in
/-- The listed set is closed under composition. -/
theorem comp_closed :
    (perms.all fun p => perms.all fun q => perms.contains (comp p q)) = true := by
  decide

set_option maxHeartbeats 4000000 in
/-- Every listed permutation has a listed two-sided inverse. -/
theorem inverses_listed :
    (perms.all fun p => perms.any fun q =>
      comp p q == List.range 12 && comp q p == List.range 12) = true := by
  decide

set_option maxHeartbeats 4000000 in
/-- Every listed permutation preserves port adjacency. -/
theorem adjacency_preserved :
    (perms.all fun p => (List.range 12).all fun i => (List.range 12).all fun j =>
      natAdj i j == natAdj (app p i) (app p j)) = true := by
  decide

/-- Every listed permutation commutes with the antipode `k ↦ 11 - k`. -/
theorem antipode_equivariant :
    (perms.all fun p => (List.range 12).all fun k =>
      app p (11 - k) == 11 - app p k) = true := by
  decide

/-- The action is transitive: every port is the image of port `0`. -/
theorem transitive_on_ports :
    ((List.range 12).all fun t => perms.any fun p => app p 0 == t) = true := by
  decide

end OPH.A5PortAction
